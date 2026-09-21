#!/usr/bin/env python3
"""Slice LINCS 2020 Level 5 signatures for CLDN4 and TACSTD2.

The shRNA gctx is ~12 GB and the CRISPR gctx is ~6.5 GB. This script downloads
metadata, then HTTP-range-reads only the signature rows that match. It does not
write the full archives.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import io
import time
from collections import Counter
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import requests

ROOT = Path(__file__).resolve().parent
CACHE = ROOT / "data" / "cache"
DATA = ROOT / "data"
TABLES = ROOT / "tables"

SIGINFO_URL = "https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/siginfo_beta.txt"
GENEINFO_URL = "https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/geneinfo_beta.txt"
CELLINFO_URL = "https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/cellinfo_beta.txt"
XPR_GCTX_URL = (
    "https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/level5/"
    "level5_beta_trt_xpr_n142901x12328.gctx"
)
PHASE1_SIG_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/"
    "GSE92742_Broad_LINCS_sig_info.txt.gz"
)
PHASE1_PERT_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/"
    "GSE92742_Broad_LINCS_pert_info.txt.gz"
)
PHASE2_SIG_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/"
    "GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz"
)

GENES = ("CLDN4", "TACSTD2")
# Loss-of-function classes. trt_sh is shRNA. trt_xpr is CRISPR. trt_si is siRNA.
LOF_TYPES = ("trt_sh", "trt_sh.cgs", "trt_sh.css", "trt_xpr", "trt_si")


class HTTPRangeFile:
    """Read-only file object over HTTP Range requests, with a block cache."""

    def __init__(self, url: str, block: int = 4 * 1024 * 1024, max_blocks: int = 48):
        self.session = requests.Session()
        head = self.session.head(url, timeout=60, allow_redirects=True)
        head.raise_for_status()
        self.url = head.url
        self.size = int(head.headers["Content-Length"])
        self.block = block
        self.max_blocks = max_blocks
        self.pos = 0
        self.cache: dict[int, bytes] = {}
        self.order: list[int] = []
        self.n_req = 0
        self.n_bytes = 0

    def seek(self, offset: int, whence: int = 0) -> int:
        if whence == 0:
            self.pos = offset
        elif whence == 1:
            self.pos += offset
        elif whence == 2:
            self.pos = self.size + offset
        else:
            raise ValueError(whence)
        if self.pos < 0:
            raise ValueError("negative seek")
        return self.pos

    def tell(self) -> int:
        return self.pos

    def seekable(self) -> bool:
        return True

    def readable(self) -> bool:
        return True

    def writable(self) -> bool:
        return False

    def close(self) -> None:
        return None

    def flush(self) -> None:
        return None

    def read(self, n: int = -1) -> bytes:
        if n is None or n < 0:
            n = self.size - self.pos
        if n <= 0 or self.pos >= self.size:
            return b""
        end = min(self.pos + n, self.size)
        out = bytearray()
        while self.pos < end:
            data = self._block(self.pos // self.block)
            inner = self.pos % self.block
            take = min(len(data) - inner, end - self.pos)
            if take <= 0:
                break
            out += data[inner : inner + take]
            self.pos += take
        return bytes(out)

    def _block(self, bidx: int) -> bytes:
        hit = self.cache.get(bidx)
        if hit is not None:
            return hit
        start = bidx * self.block
        end = min(start + self.block, self.size) - 1
        last_err: Exception | None = None
        for attempt in range(4):
            try:
                resp = self.session.get(
                    self.url,
                    headers={"Range": f"bytes={start}-{end}"},
                    timeout=180,
                )
                if resp.status_code not in (200, 206):
                    raise RuntimeError(f"HTTP {resp.status_code} for bytes {start}-{end}")
                data = resp.content
                break
            except Exception as exc:  # noqa: BLE001 — retry transient range failures
                last_err = exc
                time.sleep(2 ** attempt)
        else:
            raise RuntimeError(f"range read failed at {start}-{end}") from last_err
        self.n_req += 1
        self.n_bytes += len(data)
        self.cache[bidx] = data
        self.order.append(bidx)
        while len(self.order) > self.max_blocks:
            self.cache.pop(self.order.pop(0), None)
        return data


def _download(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        print(f"cache hit {dest.name} ({dest.stat().st_size} bytes)")
        return dest
    print(f"GET {url}")
    with requests.get(url, stream=True, timeout=300) as resp:
        resp.raise_for_status()
        tmp = dest.with_suffix(dest.suffix + ".partial")
        with tmp.open("wb") as fh:
            for chunk in resp.iter_content(chunk_size=8 * 1024 * 1024):
                if chunk:
                    fh.write(chunk)
        tmp.replace(dest)
    print(f"wrote {dest.name} ({dest.stat().st_size} bytes)")
    return dest


def _open_text(path: Path):
    if path.suffix == ".gz":
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", newline="")
    return path.open(newline="")


def _count_names(path: Path, name_field: str) -> dict:
    """Count signature rows and CLDN4/TACSTD2 hits. Streamed, no full frame."""
    n = 0
    by_type: Counter[str] = Counter()
    genes_by_type: dict[str, set[str]] = {}
    hits = []
    with _open_text(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = reader.fieldnames or []
        # Phase I uses pert_iname; LINCS2020 uses cmap_name.
        if name_field not in fields:
            raise SystemExit(f"{path.name} missing {name_field}; columns={fields[:12]}")
        type_field = "pert_type" if "pert_type" in fields else None
        for row in reader:
            n += 1
            pt = row[type_field] if type_field else ""
            by_type[pt] += 1
            name = row[name_field]
            genes_by_type.setdefault(pt, set()).add(name)
            if name in GENES or any(g in name for g in GENES):
                hits.append(row)
    return {
        "n": n,
        "by_type": by_type,
        "n_genes_by_type": {k: len(v) for k, v in genes_by_type.items()},
        "hits": hits,
        "fields": fields,
    }


def _decode_fixed(arr) -> list[str]:
    out = []
    for item in arr:
        if isinstance(item, bytes):
            out.append(item.decode("utf-8").rstrip("\x00"))
        else:
            out.append(str(item).rstrip("\x00"))
    return out


def slice_tacstd2(
    meta: pd.DataFrame, geneinfo: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if meta.empty:
        raise SystemExit("No TACSTD2 trt_xpr rows in siginfo; nothing to slice.")
    remote = HTTPRangeFile(XPR_GCTX_URL)
    print(f"xpr gctx bytes {remote.size}")
    t0 = time.time()
    handle = h5py.File(remote, "r", rdcc_nbytes=128 * 1024 * 1024, rdcc_nslots=20011)
    matrix = handle["0/DATA/0/matrix"]
    col_ids = _decode_fixed(handle["0/META/COL/id"][:])
    row_ids = _decode_fixed(handle["0/META/ROW/id"][:])
    print(
        f"opened gctx shape={matrix.shape} chunks={matrix.chunks} "
        f"ids={len(col_ids)} genes={len(row_ids)} "
        f"http_MB={remote.n_bytes/1e6:.1f} s={time.time()-t0:.1f}"
    )
    if matrix.shape != (len(col_ids), len(row_ids)):
        raise SystemExit(f"unexpected gctx shape {matrix.shape}")
    id_to_row = {sid: i for i, sid in enumerate(col_ids)}
    missing = [s for s in meta["sig_id"] if s not in id_to_row]
    if missing:
        raise SystemExit(f"{len(missing)} sig_ids missing from gctx, e.g. {missing[:3]}")
    geneinfo = geneinfo.copy()
    geneinfo["gene_id"] = geneinfo["gene_id"].astype(str)
    if set(row_ids) != set(geneinfo["gene_id"]):
        raise SystemExit("gctx row ids do not match geneinfo_beta gene_id")
    if row_ids != geneinfo["gene_id"].tolist():
        geneinfo = geneinfo.set_index("gene_id").loc[row_ids].reset_index()
    order = np.array([id_to_row[s] for s in meta["sig_id"]], dtype=int)
    data = np.empty((len(order), matrix.shape[1]), dtype=np.float32)
    for out_i, src_i in enumerate(order):
        data[out_i, :] = matrix[int(src_i), :]
        print(
            f"  sliced {out_i+1}/{len(order)} row {int(src_i)} "
            f"http_MB={remote.n_bytes/1e6:.1f}"
        )
    if not np.isfinite(data).all():
        n_bad = int((~np.isfinite(data)).sum())
        raise SystemExit(f"non-finite z-scores in slice: {n_bad}")
    handle.close()
    wide = pd.DataFrame(data.T, index=row_ids)
    wide.columns = meta["sig_id"].tolist()
    wide.insert(0, "gene_symbol", geneinfo["gene_symbol"].tolist())
    wide.insert(1, "feature_space", geneinfo["feature_space"].tolist())
    wide.index.name = "gene_id"
    meta = meta.copy()
    meta["gctx_row"] = order
    meta["gctx_bytes_fetched"] = remote.n_bytes
    meta["gctx_http_requests"] = remote.n_req
    print(
        f"slice done http_MB={remote.n_bytes/1e6:.1f} "
        f"requests={remote.n_req} of gctx_MB={remote.size/1e6:.1f} "
        f"s={time.time()-t0:.1f}"
    )
    return wide, meta


def main() -> None:
    CACHE.mkdir(parents=True, exist_ok=True)
    TABLES.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)

    sig_path = _download(SIGINFO_URL, CACHE / "siginfo_beta.txt")
    gene_path = _download(GENEINFO_URL, CACHE / "geneinfo_beta.txt")
    cell_path = _download(CELLINFO_URL, CACHE / "cellinfo_beta.txt")
    p1_sig = _download(PHASE1_SIG_URL, CACHE / "GSE92742_sig_info.txt.gz")
    p1_pert = _download(PHASE1_PERT_URL, CACHE / "GSE92742_pert_info.txt.gz")
    p2_sig = _download(PHASE2_SIG_URL, CACHE / "GSE70138_sig_info.txt.gz")

    print("scanning LINCS2020 siginfo")
    lincs = _count_names(sig_path, "cmap_name")
    print("scanning GSE92742 sig_info")
    phase1 = _count_names(p1_sig, "pert_iname")
    print("scanning GSE92742 pert_info")
    phase1_pert = _count_names(p1_pert, "pert_iname")
    print("scanning GSE70138 sig_info")
    phase2 = _count_names(p2_sig, "pert_iname")

    geneinfo = pd.read_csv(gene_path, sep="\t", dtype=str)
    cells = pd.read_csv(cell_path, sep="\t", dtype=str)
    cell_keep = cells[
        ["cell_iname", "cell_lineage", "primary_disease", "subtype", "cell_type"]
    ].drop_duplicates("cell_iname")

    hits = pd.DataFrame(lincs["hits"])
    if hits.empty:
        hits = pd.DataFrame(columns=["cmap_name", "pert_type", "sig_id"])
    inventory_rows = []
    inventory_rows.append(
        {
            "source": "LINCS2020 siginfo_beta",
            "n_rows": lincs["n"],
            "cldn4_rows": int((hits["cmap_name"] == "CLDN4").sum()) if len(hits) else 0,
            "tacstd2_rows": int((hits["cmap_name"] == "TACSTD2").sum()) if len(hits) else 0,
            "note": "cmap_name exact",
        }
    )
    for label, scanned, name_col in (
        ("GSE92742 sig_info (phase I)", phase1, "pert_iname"),
        ("GSE92742 pert_info (phase I)", phase1_pert, "pert_iname"),
        ("GSE70138 sig_info (phase II, 2017-03-06)", phase2, "pert_iname"),
    ):
        h = pd.DataFrame(scanned["hits"])
        inventory_rows.append(
            {
                "source": label,
                "n_rows": scanned["n"],
                "cldn4_rows": int((h[name_col] == "CLDN4").sum()) if len(h) else 0,
                "tacstd2_rows": int((h[name_col] == "TACSTD2").sum()) if len(h) else 0,
                "note": f"{name_col} exact",
            }
        )
    inv = pd.DataFrame(inventory_rows)
    type_rows = []
    for pt, n in sorted(lincs["by_type"].items(), key=lambda kv: -kv[1]):
        type_rows.append(
            {
                "pert_type": pt,
                "n_sig": n,
                "n_cmap_name": lincs["n_genes_by_type"].get(pt, 0),
                "cldn4": int(((hits["pert_type"] == pt) & (hits["cmap_name"] == "CLDN4")).sum())
                if len(hits)
                else 0,
                "tacstd2": int(((hits["pert_type"] == pt) & (hits["cmap_name"] == "TACSTD2")).sum())
                if len(hits)
                else 0,
            }
        )
    types = pd.DataFrame(type_rows)
    inv.to_csv(TABLES / "inventory_sources.tsv", sep="\t", index=False)
    types.to_csv(TABLES / "inventory_pert_type.tsv", sep="\t", index=False)
    print(inv.to_string(index=False))
    print(types.to_string(index=False))

    if (hits["cmap_name"] == "CLDN4").any():
        raise SystemExit("CLDN4 rows appeared; this script's empty-CLDN4 branch is stale")

    tac = hits[(hits["cmap_name"] == "TACSTD2") & (hits["pert_type"] == "trt_xpr")].copy()
    other = hits[hits["cmap_name"] == "TACSTD2"]["pert_type"].value_counts().to_dict()
    if other != {"trt_xpr": len(tac)}:
        # Still slice CRISPR, but do not hide another class if one appears later.
        print("TACSTD2 pert_type counts", other)
    if tac.empty:
        raise SystemExit("TACSTD2 trt_xpr not found")

    num_cols = [
        "nsample",
        "cc_q75",
        "ss_ngene",
        "tas",
        "pct_self_rank_q25",
        "pert_time",
        "is_hiq",
        "qc_pass",
        "is_exemplar_sig",
        "is_ncs_sig",
        "is_null_sig",
    ]
    for col in num_cols:
        tac[col] = pd.to_numeric(tac[col], errors="coerce")
    tac = tac.merge(cell_keep, on="cell_iname", how="left")
    if tac["cell_lineage"].isna().any():
        missing_cells = tac.loc[tac["cell_lineage"].isna(), "cell_iname"].unique()
        raise SystemExit(f"cellinfo missing {missing_cells}")
    tac = tac.sort_values(["cell_iname", "pert_id", "sig_id"]).reset_index(drop=True)

    keep_cols = [
        "sig_id",
        "pert_id",
        "pert_type",
        "cmap_name",
        "cell_iname",
        "cell_lineage",
        "primary_disease",
        "subtype",
        "cell_type",
        "pert_itime",
        "pert_time",
        "nsample",
        "cc_q75",
        "ss_ngene",
        "tas",
        "pct_self_rank_q25",
        "is_hiq",
        "qc_pass",
        "is_exemplar_sig",
        "is_ncs_sig",
        "is_null_sig",
        "bead_batch",
        "det_plates",
        "distil_ids",
        "project_code",
        "build_name",
    ]
    tac = tac[keep_cols]
    wide, tac = slice_tacstd2(tac, geneinfo)
    # slice_tacstd2 returns meta with gctx fields appended; keep column order stable.
    meta_cols = keep_cols + ["gctx_row", "gctx_bytes_fetched", "gctx_http_requests"]
    tac = tac[meta_cols]

    meta_path = DATA / "tacstd2_xpr_sig_meta.tsv"
    tac.to_csv(meta_path, sep="\t", index=False)
    z_path = DATA / "tacstd2_xpr_level5_z.tsv.gz"
    wide.to_csv(z_path, sep="\t", compression="gzip")
    digest = hashlib.sha256(z_path.read_bytes()).hexdigest()
    manifest = pd.DataFrame(
        [
            {
                "file": z_path.name,
                "n_genes": wide.shape[0],
                "n_sig": wide.shape[1] - 2,
                "sha256": digest,
                "bytes": z_path.stat().st_size,
                "gctx_bytes_fetched": int(tac["gctx_bytes_fetched"].iloc[0]),
                "gctx_http_requests": int(tac["gctx_http_requests"].iloc[0]),
                "gctx_url": XPR_GCTX_URL,
            }
        ]
    )
    manifest.to_csv(DATA / "slice_manifest.tsv", sep="\t", index=False)
    print(manifest.to_string(index=False))
    print(f"wrote {meta_path}")
    print(f"wrote {z_path}")


if __name__ == "__main__":
    main()
