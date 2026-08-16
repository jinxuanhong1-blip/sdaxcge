#!/usr/bin/env python3
"""Extract the gene panel + per-cell library size from public GEO objects.

Outputs one cells-x-genes TSV.gz per cohort under --datadir.
Does not load full dense matrices when a stream is enough.
"""
from __future__ import annotations

import argparse
import gzip
import shutil
import tarfile
import tempfile
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent


def load_panel(path: Path) -> list[str]:
    return list(pd.read_csv(path, sep="\t")["gene"])


def write_panel_df(barcodes, total_umi, kept: dict[str, np.ndarray], out: Path) -> None:
    df = pd.DataFrame({"barcode": barcodes, "total_umi": np.asarray(total_umi, dtype=np.int64)})
    for g in sorted(kept):
        df[g] = np.asarray(kept[g], dtype=np.float32)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, sep="\t", index=False, compression="gzip")
    print(f"wrote {out} rows={len(df)} genes={len(kept)}", flush=True)


def extract_gse207422(datadir: Path, panel: set[str]) -> None:
    matrix = datadir / "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"
    out = datadir / "gse207422_panel_cells.tsv.gz"
    if out.exists():
        print("skip existing", out, flush=True)
        return
    print("=== extract GSE207422 ===", flush=True)
    with gzip.open(matrix, "rb") as fh:
        header = fh.readline().decode("ascii").rstrip("\n")
        cells = header.split("\t")[1:]
        n_cells = len(cells)
        print(f"cells {n_cells}", flush=True)
        total = np.zeros(n_cells, dtype=np.float64)
        kept: dict[str, np.ndarray] = {}
        n_genes = 0
        for raw in fh:
            n_genes += 1
            tab = raw.find(b"\t")
            gene = raw[:tab].decode("ascii")
            arr = np.fromstring(raw[tab + 1 :].rstrip(b"\r\n"), sep="\t", dtype=np.float32)
            if arr.size != n_cells:
                raise ValueError(f"{gene}: got {arr.size}, expected {n_cells}")
            total += arr
            if gene in panel:
                kept[gene] = arr.copy()
            if n_genes % 2000 == 0:
                print(f"  genes={n_genes} kept={len(kept)}", flush=True)
    print(f"done genes={n_genes} kept={len(kept)} missing={sorted(panel - set(kept))}", flush=True)
    write_panel_df(cells, total, kept, out)


def _ungzip_until_rds(src: Path, dest: Path) -> None:
    current = src
    tmp_paths = []
    for _ in range(3):
        with open(current, "rb") as fh:
            magic = fh.read(2)
        if magic == b"\x1f\x8b":
            nxt = dest.with_suffix(dest.suffix + f".pass{len(tmp_paths)}")
            with gzip.open(current, "rb") as zin, open(nxt, "wb") as zout:
                shutil.copyfileobj(zin, zout, 16 * 1024 * 1024)
            tmp_paths.append(nxt)
            current = nxt
            continue
        break
    if current != dest:
        shutil.copyfile(current, dest)
    for p in tmp_paths:
        if p.exists() and p != dest:
            p.unlink(missing_ok=True)


def extract_gse205335(datadir: Path, panel: list[str]) -> None:
    import rdata
    from scipy import sparse

    rds_gz = datadir / "GSE205335_Lung_IO_UMI_matrix.rds.gz"
    out = datadir / "gse205335_panel_cells.tsv.gz"
    if out.exists():
        print("skip existing", out, flush=True)
        return
    print("=== extract GSE205335 ===", flush=True)
    with tempfile.TemporaryDirectory(prefix="gse205335-") as tmp:
        rds_path = Path(tmp) / "matrix.rds"
        _ungzip_until_rds(rds_gz, rds_path)
        print(f"rds bytes={rds_path.stat().st_size}", flush=True)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message='Missing constructor for R class "dgCMatrix"')
            obj = rdata.read_rds(rds_path)
    expected = {"i", "p", "Dim", "Dimnames", "x"}
    if not expected.issubset(vars(obj)):
        raise TypeError(f"RDS is not dgCMatrix; keys={list(vars(obj))}")
    genes = np.asarray(obj.Dimnames[0], dtype=str)
    barcodes = np.asarray(obj.Dimnames[1], dtype=str)
    matrix = sparse.csc_matrix((obj.x, obj.i, obj.p), shape=tuple(obj.Dim), copy=False)
    print(f"matrix {matrix.shape[0]} genes x {matrix.shape[1]} cells", flush=True)
    gene_index = {g: i for i, g in enumerate(genes)}
    present = [g for g in panel if g in gene_index]
    missing = [g for g in panel if g not in gene_index]
    print(f"panel present={len(present)} missing={missing}", flush=True)
    total_umi = np.asarray(matrix.sum(axis=0)).ravel()
    kept = {g: np.asarray(matrix.getrow(int(gene_index[g])).toarray()).ravel() for g in present}
    write_panel_df(barcodes, total_umi, kept, out)


def _read_tsv_names(path: Path, col: int = 0) -> list[str]:
    names = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            names.append(parts[col] if len(parts) > col else parts[0])
    return names


def stream_mtx_panel(mtx_path: Path, features: list[str], n_cells: int, panel: set[str]) -> dict[str, np.ndarray]:
    """Stream a genes x cells MatrixMarket file; keep panel gene rows only."""
    wanted_rows = {}
    for i, g in enumerate(features, start=1):
        if g in panel and g not in wanted_rows:
            wanted_rows[i] = g
    kept = {g: np.zeros(n_cells, dtype=np.float32) for g in wanted_rows.values()}
    print(f"  streaming {mtx_path.name} want_rows={len(wanted_rows)} n_cells={n_cells}", flush=True)
    n_kept_entries = 0
    with gzip.open(mtx_path, "rt") as fh:
        header = None
        for line in fh:
            if line.startswith("%"):
                continue
            header = line.strip().split()
            break
        if header is None:
            raise ValueError(f"empty MTX {mtx_path}")
        n_rows, n_cols, n_entries = map(int, header)
        if n_cols != n_cells:
            raise ValueError(f"{mtx_path}: MTX cols {n_cols} != barcodes {n_cells}")
        print(f"  MTX {n_rows} x {n_cols} nnz={n_entries}", flush=True)
        for i, line in enumerate(fh, start=1):
            r_s, c_s, v_s = line.split()
            r = int(r_s)
            if r in wanted_rows:
                kept[wanted_rows[r]][int(c_s) - 1] = float(v_s)
                n_kept_entries += 1
            if i % 50_000_000 == 0:
                print(f"    entries {i}/{n_entries} kept={n_kept_entries}", flush=True)
    print(f"  kept_entries={n_kept_entries} genes={list(kept)}", flush=True)
    return kept


def extract_gse241934_split(datadir: Path, panel: set[str], split: str) -> None:
    if split == "IIT":
        mtx = datadir / "GSE241934_IIT_Matrix.mtx.gz"
        feat = datadir / "GSE241934_IIT_features.tsv.gz"
        bc = datadir / "GSE241934_IIT_barcodes.tsv.gz"
        meta_p = datadir / "GSE241934_IIT_Meta.txt.gz"
        out = datadir / "gse241934_iit_panel_cells.tsv.gz"
    else:
        mtx = datadir / "GSE241934_Real_Matrix.mtx.gz"
        feat = datadir / "GSE241934_RWC_features.tsv.gz"
        bc = datadir / "GSE241934_RWC_barcodes.tsv.gz"
        meta_p = datadir / "GSE241934_Real_Meta.txt.gz"
        out = datadir / "gse241934_rwc_panel_cells.tsv.gz"
    if out.exists():
        print("skip existing", out, flush=True)
        return
    print(f"=== extract GSE241934 {split} ===", flush=True)
    features = _read_tsv_names(feat, 0)
    barcodes = _read_tsv_names(bc, 0)
    kept = stream_mtx_panel(mtx, features, len(barcodes), panel)
    meta = pd.read_csv(meta_p, sep="\t", low_memory=False)
    meta = meta.set_index("cellID")
    # author nCount_RNA is the library size after their QC
    total = meta.reindex(barcodes)["nCount_RNA"].to_numpy(dtype=np.float64)
    if np.isnan(total).any():
        n_miss = int(np.isnan(total).sum())
        print(f"  WARNING {n_miss} barcodes missing from meta; filling total_umi from panel sum", flush=True)
        panel_sum = np.zeros(len(barcodes), dtype=np.float64)
        for arr in kept.values():
            panel_sum += arr
        total = np.where(np.isnan(total), panel_sum, total)
    write_panel_df(barcodes, total, kept, out)
    # also write a slim meta join key
    keep_cols = [
        c for c in [
            "orig.ident", "sampleID", "major.cell.type", "major_cell_type",
            "Pathological Response", "EGFR", "PD1", "Histology",
            "nCount_RNA", "nFeature_RNA", "percent.mt",
        ] if c in meta.columns
    ]
    slim = meta.reindex(barcodes)[keep_cols].reset_index()
    if "cellID" in slim.columns:
        slim = slim.rename(columns={"cellID": "barcode"})
    elif "index" in slim.columns:
        slim = slim.rename(columns={"index": "barcode"})
    slim_out = datadir / f"gse241934_{split.lower()}_meta_join.tsv.gz"
    slim.to_csv(slim_out, sep="\t", index=False, compression="gzip")
    print("wrote", slim_out, flush=True)


def _extract_10x_tar_samples(
    tar_path: Path,
    dest_dir: Path,
    panel: set[str],
    samples: list[tuple[str, dict]],
    out_name: str,
) -> None:
    out = dest_dir / out_name
    if out.exists():
        print("skip existing", out, flush=True)
        return
    dest_dir.mkdir(parents=True, exist_ok=True)
    extract_dir = dest_dir / (tar_path.stem + "_untar")
    extract_dir.mkdir(exist_ok=True)
    needed = []
    for prefix, _meta in samples:
        for suf in ("_matrix.mtx.gz", "_barcodes.tsv.gz", "_features.tsv.gz"):
            needed.append(prefix + suf)
    have = {p.name for p in extract_dir.iterdir()} if extract_dir.exists() else set()
    if not set(needed).issubset(have):
        print(f"untar {tar_path.name}", flush=True)
        with tarfile.open(tar_path) as tf:
            members = [m for m in tf.getmembers() if Path(m.name).name in set(needed)]
            tf.extractall(extract_dir, members=members)
        # flatten if nested
        for p in extract_dir.rglob("*"):
            if p.is_file() and p.parent != extract_dir:
                target = extract_dir / p.name
                if not target.exists():
                    shutil.move(str(p), str(target))

    import scipy.io

    frames = []
    for prefix, extra in samples:
        print(f"  10x {prefix}", flush=True)
        mtx = scipy.io.mmread(str(extract_dir / f"{prefix}_matrix.mtx.gz")).tocsr().astype(np.float32)
        barcodes = _read_tsv_names(extract_dir / f"{prefix}_barcodes.tsv.gz", 0)
        genes = []
        with gzip.open(extract_dir / f"{prefix}_features.tsv.gz", "rt") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                genes.append(parts[1] if len(parts) > 1 else parts[0])
        if mtx.shape[0] != len(genes) and mtx.shape[1] == len(genes):
            mtx = mtx.T.tocsr()
        if mtx.shape[0] != len(genes) or mtx.shape[1] != len(barcodes):
            raise ValueError(f"{prefix}: mtx {mtx.shape} genes {len(genes)} bc {len(barcodes)}")
        gidx = {}
        for i, g in enumerate(genes):
            if g not in gidx:
                gidx[g] = i
        total = np.asarray(mtx.sum(axis=0)).ravel()
        rec = {"barcode": [f"{prefix}__{b}" for b in barcodes], "total_umi": total.astype(np.int64)}
        rec.update(extra)
        for g in sorted(panel):
            if g in gidx:
                rec[g] = np.asarray(mtx[gidx[g], :].todense()).ravel().astype(np.float32)
        frames.append(pd.DataFrame(rec))
        print(f"    cells={len(barcodes)} genes_kept={sum(g in gidx for g in panel)}", flush=True)
    df = pd.concat(frames, ignore_index=True)
    df.to_csv(out, sep="\t", index=False, compression="gzip")
    print(f"wrote {out} rows={len(df)}", flush=True)


def extract_gse291670(datadir: Path, panel: set[str]) -> None:
    samples = [
        ("GSM8839599_MPR-1", {"sample": "MPR-1", "patient": "MPR-1", "response": "MPR"}),
        ("GSM8839600_MPR-2", {"sample": "MPR-2", "patient": "MPR-2", "response": "MPR"}),
        ("GSM8839601_MPR-3", {"sample": "MPR-3", "patient": "MPR-3", "response": "MPR"}),
        ("GSM8839602_Non-MPR-1", {"sample": "Non-MPR-1", "patient": "Non-MPR-1", "response": "NMPR"}),
        ("GSM8839603_Non-MPR-2", {"sample": "Non-MPR-2", "patient": "Non-MPR-2", "response": "NMPR"}),
        ("GSM8839604_Non-MPR-3", {"sample": "Non-MPR-3", "patient": "Non-MPR-3", "response": "NMPR"}),
    ]
    print("=== extract GSE291670 ===", flush=True)
    _extract_10x_tar_samples(datadir / "GSE291670_RAW.tar", datadir, panel, samples, "gse291670_panel_cells.tsv.gz")


def extract_gse233203(datadir: Path, panel: set[str]) -> None:
    samples = [
        ("GSM7412612_NCCLu_162", {"sample": "NCCLu_162", "patient": "NCCLu_162", "response": "R"}),
        ("GSM7412613_NCCLu_185", {"sample": "NCCLu_185", "patient": "NCCLu_185", "response": "NR"}),
        ("GSM7412614_NCCLu_327", {"sample": "NCCLu_327", "patient": "NCCLu_327", "response": "NR"}),
        ("GSM7412615_NCCLu_334", {"sample": "NCCLu_334", "patient": "NCCLu_334", "response": "NR"}),
        ("GSM7412616_NCCLu_376", {"sample": "NCCLu_376", "patient": "NCCLu_376", "response": "R"}),
        ("GSM7412617_NCCLu_383", {"sample": "NCCLu_383", "patient": "NCCLu_383", "response": "R"}),
        ("GSM7412618_NCCLu_397", {"sample": "NCCLu_397", "patient": "NCCLu_397", "response": "NR"}),
    ]
    print("=== extract GSE233203 leftover ===", flush=True)
    _extract_10x_tar_samples(datadir / "GSE233203_RAW.tar", datadir, panel, samples, "gse233203_panel_cells.tsv.gz")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--datadir", type=Path, default=Path("/tmp/scrna_exh_meta"))
    ap.add_argument("--panel", type=Path, default=HERE / "gene_panel.tsv")
    args = ap.parse_args()
    panel_list = load_panel(args.panel)
    panel = set(panel_list)
    extract_gse207422(args.datadir, panel)
    extract_gse205335(args.datadir, panel_list)
    extract_gse241934_split(args.datadir, panel, "IIT")
    extract_gse241934_split(args.datadir, panel, "RWC")
    extract_gse291670(args.datadir, panel)
    extract_gse233203(args.datadir, panel)
    print("extract complete")


if __name__ == "__main__":
    main()
