#!/usr/bin/env python3
"""Download and catalog DepMap 24Q2 RNA + Nusinow 2020 CCLE protein for B2.

Does not invent accessions. Figshare article 25880521 and MassIVE MSV000085836
were queried over HTTPS before transfer. The 461 MB RNA matrix is hashed and
column-subsetted in one pass; the full matrix is not retained.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
import ssl
import sys
import time
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

USER_AGENT = "B2-CCLE-all/1.0 (research; +https://github.com/jinxuanhong1-blip/sdaxcge)"
FIGSHARE_ARTICLE = "https://api.figshare.com/v2/articles/25880521"
MASSIVE_PROXI = (
    "https://massive.ucsd.edu/ProteoSAFe/proxi/v0.1/datasets"
    "?accession=MSV000085836&resultType=full"
)
GYGI_PROTEIN = "https://gygi.hms.harvard.edu/data/ccle/protein_quant_current_normalized.csv.gz"
GYGI_SAMPLE = "https://gygi.hms.harvard.edu/data/ccle/Table_S1_Sample_Information.xlsx"
CLDN_RE = re.compile(r"^CLDN\d+$")
TACSTD2_RE = re.compile(r"^TACSTD2$")


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def stamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def ctx() -> ssl.SSLContext:
    return ssl.create_default_context()


def open_url(url: str, timeout: int = 600):
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    return urllib.request.urlopen(req, timeout=timeout, context=ctx())


class HashingRaw(io.RawIOBase):
    def __init__(self, raw):
        self.raw = raw
        self.h = hashlib.sha256()
        self.n = 0

    def read(self, size: int = -1) -> bytes:
        b = self.raw.read(size)
        if b:
            self.h.update(b)
            self.n += len(b)
        return b

    def readinto(self, b) -> int:
        data = self.read(len(b))
        n = len(data)
        b[:n] = data
        return n

    def readable(self) -> bool:
        return True

    def writable(self) -> bool:
        return False

    def seekable(self) -> bool:
        return False


def parse_symbol(col: str) -> str:
    col = col.strip()
    if not col:
        return ""
    # DepMap: "TACSTD2 (4070)"
    m = re.match(r"^([A-Za-z0-9\-._]+)(?:\s+\(\d+\))?$", col)
    return m.group(1) if m else col


def save_json(url: str, dest: Path) -> dict:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with open_url(url, timeout=60) as resp:
        raw = resp.read()
    dest.write_bytes(raw)
    return json.loads(raw.decode("utf-8"))


def download_to(url: str, dest: Path) -> tuple[int, str, int]:
    dest.parent.mkdir(parents=True, exist_ok=True)
    hasher = hashlib.sha256()
    n = 0
    http_len = None
    with open_url(url) as resp:
        cl = resp.headers.get("Content-Length")
        if cl and cl.isdigit():
            http_len = int(cl)
        with dest.open("wb") as out:
            while True:
                chunk = resp.read(1024 * 1024)
                if not chunk:
                    break
                hasher.update(chunk)
                n += len(chunk)
                out.write(chunk)
    return n, hasher.hexdigest(), http_len if http_len is not None else n


def stream_extract_rna(url: str, out_tsv: Path, means_json: Path) -> dict:
    """Hash the full CSV while writing TACSTD2 + CLDN\\d+ columns and gene means."""
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    with open_url(url) as resp:
        cl = resp.headers.get("Content-Length")
        http_len = int(cl) if cl and cl.isdigit() else None
        hashed = HashingRaw(resp)
        text = io.TextIOWrapper(hashed, encoding="utf-8", newline="")
        reader = csv.reader(text)
        header = next(reader)
        symbols = [parse_symbol(c) for c in header]
        keep_idx = []
        keep_sym = []
        for i, s in enumerate(symbols):
            if i == 0:
                continue
            if TACSTD2_RE.match(s) or CLDN_RE.match(s):
                keep_idx.append(i)
                keep_sym.append(s)
        if "TACSTD2" not in keep_sym:
            raise SystemExit("TACSTD2 column not found in DepMap RNA matrix")
        if "CLDN4" not in keep_sym:
            raise SystemExit("CLDN4 column not found in DepMap RNA matrix")

        n_genes = len(header) - 1
        sums = [0.0] * n_genes
        counts = [0] * n_genes
        n_lines = 0
        with out_tsv.open("w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["ModelID"] + keep_sym)
            for row in reader:
                if not row:
                    continue
                mid = row[0]
                out = [mid]
                for i in keep_idx:
                    val = row[i] if i < len(row) else ""
                    out.append(val)
                w.writerow(out)
                for gi in range(n_genes):
                    j = gi + 1
                    if j >= len(row) or row[j] == "":
                        continue
                    try:
                        x = float(row[j])
                    except ValueError:
                        continue
                    sums[gi] += x
                    counts[gi] += 1
                n_lines += 1
                if n_lines % 200 == 0:
                    print(f"  RNA rows {n_lines} bytes={hashed.n}", file=sys.stderr)

    means = {}
    for gi in range(n_genes):
        s = symbols[gi + 1]
        if not s or counts[gi] == 0:
            continue
        means[s] = {
            "mean": sums[gi] / counts[gi],
            "n": counts[gi],
        }
    means_json.write_text(json.dumps({"n_lines": n_lines, "n_genes": len(means), "means": means}))
    print(
        f"  RNA done lines={n_lines} keep={keep_sym} bytes={hashed.n} s={time.time()-t0:.1f}",
        file=sys.stderr,
    )
    return {
        "bytes": hashed.n,
        "sha256": hashed.h.hexdigest(),
        "http_len": http_len if http_len is not None else hashed.n,
        "n_lines": n_lines,
        "keep": keep_sym,
        "n_genes_with_mean": len(means),
    }


def catalog_row(**kwargs) -> dict:
    keys = [
        "sample_id",
        "accession",
        "database",
        "source_url",
        "verified_at_utc",
        "verification_evidence",
        "file_name",
        "bytes",
        "decision",
        "refusal_reason",
        "sha256",
        "organism",
        "assay",
        "layout",
        "read_pair",
        "registry_checksum",
        "local_path",
        "downloaded_at_utc",
        "notes",
    ]
    row = {k: kwargs.get(k, "NA") for k in keys}
    return row


def write_catalog(path: Path, rows: list[dict]) -> None:
    keys = list(rows[0].keys())
    rows = sorted(rows, key=lambda r: (r["sample_id"], r["accession"], r["file_name"]))
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", default="/tmp/ccle_data")
    ap.add_argument("--outdir", default="results/w200/B2_CCLE_all")
    args = ap.parse_args()

    work = Path(args.workdir)
    out = Path(args.outdir)
    extracted = work / "extracted"
    repro = out / "repro" / "accessions"
    work.mkdir(parents=True, exist_ok=True)
    extracted.mkdir(parents=True, exist_ok=True)
    repro.mkdir(parents=True, exist_ok=True)

    verified_at = utc_now()
    fig_path = repro / f"figshare.25880521.{stamp()}.json"
    print("Verifying Figshare 25880521", file=sys.stderr)
    article = save_json(FIGSHARE_ARTICLE, fig_path)
    if str(article.get("id")) != "25880521":
        raise SystemExit(f"Figshare id mismatch: {article.get('id')}")
    if "DepMap 24Q2" not in str(article.get("title", "")):
        raise SystemExit(f"Unexpected Figshare title: {article.get('title')}")

    files = {f["name"]: f for f in article.get("files", [])}
    need = [
        "Model.csv",
        "OmicsExpressionProteinCodingGenesTPMLogp1.csv",
        "README.txt",
    ]
    for n in need:
        if n not in files:
            raise SystemExit(f"Missing Figshare file {n}")

    msv_path = repro / f"MSV000085836.{stamp()}.json"
    print("Verifying MassIVE MSV000085836", file=sys.stderr)
    msv = save_json(MASSIVE_PROXI, msv_path)
    rec = msv[0] if isinstance(msv, list) else msv
    acc_vals = []
    for a in rec.get("accession", []):
        if isinstance(a, dict):
            acc_vals.append(str(a.get("value", "")))
    if "MSV000085836" not in acc_vals:
        raise SystemExit(f"MassIVE accession mismatch: {acc_vals}")
    if "Cancer Cell Line Encyclopedia" not in rec.get("title", ""):
        raise SystemExit(f"Unexpected MassIVE title: {rec.get('title')}")

    rows = []
    # Model.csv
    meta = files["Model.csv"]
    dest = work / "Model.csv"
    print("Downloading Model.csv", file=sys.stderr)
    n, sha, http_len = download_to(meta["download_url"], dest)
    if n != int(meta["size"]):
        print(f"WARNING Model.csv size {n} != registry {meta['size']}", file=sys.stderr)
    rows.append(
        catalog_row(
            sample_id="DepMap-24Q2-models",
            accession="10.25452/figshare.plus.25880521",
            database="Figshare+",
            source_url=meta["download_url"],
            verified_at_utc=verified_at,
            verification_evidence=str(fig_path.as_posix().replace("/workspace/", "")),
            file_name="Model.csv",
            bytes=str(int(meta["size"])),
            decision="accepted",
            refusal_reason="NA",
            sha256=sha,
            organism="Homo sapiens",
            assay="cell_line_metadata",
            layout="NA",
            read_pair="NA",
            registry_checksum=meta.get("computed_md5") or "NA",
            local_path=str(dest),
            downloaded_at_utc=utc_now(),
            notes=f"figshare file id {meta['id']}; downloaded_bytes={n}; http_len={http_len}",
        )
    )

    # README
    meta = files["README.txt"]
    dest = work / "README.txt"
    print("Downloading README.txt", file=sys.stderr)
    n, sha, http_len = download_to(meta["download_url"], dest)
    rows.append(
        catalog_row(
            sample_id="DepMap-24Q2-readme",
            accession="10.25452/figshare.plus.25880521",
            database="Figshare+",
            source_url=meta["download_url"],
            verified_at_utc=verified_at,
            verification_evidence=str(fig_path.as_posix().replace("/workspace/", "")),
            file_name="README.txt",
            bytes=str(int(meta["size"])),
            decision="accepted",
            refusal_reason="NA",
            sha256=sha,
            organism="Homo sapiens",
            assay="release_readme",
            layout="NA",
            read_pair="NA",
            registry_checksum=meta.get("computed_md5") or "NA",
            local_path=str(dest),
            downloaded_at_utc=utc_now(),
            notes=f"figshare file id {meta['id']}; downloaded_bytes={n}",
        )
    )

    # RNA matrix: stream extract
    meta = files["OmicsExpressionProteinCodingGenesTPMLogp1.csv"]
    print("Streaming RNA matrix (hash + extract TACSTD2/CLDN*)", file=sys.stderr)
    rna_info = stream_extract_rna(
        meta["download_url"],
        extracted / "rna_tacstd2_cldn.tsv",
        extracted / "rna_gene_means.json",
    )
    if rna_info["bytes"] != int(meta["size"]):
        print(
            f"WARNING RNA size {rna_info['bytes']} != registry {meta['size']}",
            file=sys.stderr,
        )
    rows.append(
        catalog_row(
            sample_id="DepMap-24Q2-RNA",
            accession="10.25452/figshare.plus.25880521",
            database="Figshare+",
            source_url=meta["download_url"],
            verified_at_utc=verified_at,
            verification_evidence=str(fig_path.as_posix().replace("/workspace/", "")),
            file_name="OmicsExpressionProteinCodingGenesTPMLogp1.csv",
            bytes=str(int(meta["size"])),
            decision="accepted",
            refusal_reason="NA",
            sha256=rna_info["sha256"],
            organism="Homo sapiens",
            assay="RNAseq_log2_TPM_plus1_protein_coding",
            layout="NA",
            read_pair="NA",
            registry_checksum=meta.get("computed_md5") or "NA",
            local_path=str(extracted / "rna_tacstd2_cldn.tsv"),
            downloaded_at_utc=utc_now(),
            notes=(
                f"figshare file id {meta['id']}; full file hashed, not retained; "
                f"extracted n_lines={rna_info['n_lines']} genes={','.join(rna_info['keep'])}"
            ),
        )
    )

    # Gygi protein (processed table cited by MassIVE record as the CCLE portal table)
    print("Downloading Nusinow protein table from Gygi lab", file=sys.stderr)
    dest = work / "protein_quant_current_normalized.csv.gz"
    n, sha, http_len = download_to(GYGI_PROTEIN, dest)
    rows.append(
        catalog_row(
            sample_id="CCLE-Nusinow2020-protein",
            accession="MSV000085836",
            database="MassIVE",
            source_url=GYGI_PROTEIN,
            verified_at_utc=verified_at,
            verification_evidence=str(msv_path.as_posix().replace("/workspace/", "")),
            file_name="protein_quant_current_normalized.csv.gz",
            bytes=str(n),
            decision="accepted",
            refusal_reason="NA",
            sha256=sha,
            organism="Homo sapiens",
            assay="TMT_MS3_proteomics_normalized",
            layout="NA",
            read_pair="NA",
            registry_checksum="NA",
            local_path=str(dest),
            downloaded_at_utc=utc_now(),
            notes=(
                "Processed gene-level table hosted by Gygi lab "
                "(https://gygi.hms.harvard.edu/publications/ccle.html). "
                f"Raw spectra: MassIVE MSV000085836. http_len={http_len}."
            ),
        )
    )

    print("Downloading Nusinow Table S1 sample information", file=sys.stderr)
    dest = work / "Table_S1_Sample_Information.xlsx"
    n, sha, http_len = download_to(GYGI_SAMPLE, dest)
    rows.append(
        catalog_row(
            sample_id="CCLE-Nusinow2020-samples",
            accession="MSV000085836",
            database="MassIVE",
            source_url=GYGI_SAMPLE,
            verified_at_utc=verified_at,
            verification_evidence=str(msv_path.as_posix().replace("/workspace/", "")),
            file_name="Table_S1_Sample_Information.xlsx",
            bytes=str(n),
            decision="accepted",
            refusal_reason="NA",
            sha256=sha,
            organism="Homo sapiens",
            assay="sample_metadata",
            layout="NA",
            read_pair="NA",
            registry_checksum="NA",
            local_path=str(dest),
            downloaded_at_utc=utc_now(),
            notes=f"Gygi Table S1; http_len={http_len}",
        )
    )

    cat = out / "catalog.tsv"
    write_catalog(cat, rows)
    notes_cat = Path("notes/w200/B2_CCLE_all/catalog.tsv")
    notes_cat.parent.mkdir(parents=True, exist_ok=True)
    notes_cat.write_bytes(cat.read_bytes())
    print(f"Wrote {cat} rows={len(rows)}", file=sys.stderr)


if __name__ == "__main__":
    main()
