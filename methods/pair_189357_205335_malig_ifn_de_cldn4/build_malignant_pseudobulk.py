#!/usr/bin/env python3
"""UMI-sum malignant cells to patient pseudobulk.

GSE189357: marker-malignant (EPCAM|KRT8|KRT18|KRT19)>0 and PTPRC==0.
GSE205335: author lineage.sub == 'Malignant cells'.

PR #459 locked patients only (n=9 + n=22). No dual-high. No T/NK.
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
from scipy import sparse
from scipy.io import mmread

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
EPI = ["EPCAM", "KRT8", "KRT18", "KRT19"]
DEFAULT_GEO = Path("/tmp/geo_pair_189357_205335")


def _upper_map(names: list[str]) -> dict[str, int]:
    return {n.upper(): i for i, n in enumerate(names)}


def _mal_mask(counts: np.ndarray, idx: dict[str, int]) -> np.ndarray:
    def col(g: str) -> np.ndarray:
        if g not in idx:
            return np.zeros(counts.shape[0], dtype=float)
        return counts[:, idx[g]].astype(float)

    epi = np.zeros(counts.shape[0], dtype=bool)
    for g in EPI:
        epi |= col(g) > 0
    return epi & (col("PTPRC") == 0)


def parse_geo_soft(path: Path) -> pd.DataFrame:
    opener = gzip.open if path.suffix == ".gz" else open
    records: list[dict[str, str]] = []
    current: dict[str, str] | None = None
    descriptions: list[str] = []
    titles: list[str] = []
    with opener(path, "rt", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current is not None:
                    current["description"] = descriptions[0] if descriptions else ""
                    current["title"] = titles[0] if titles else ""
                    records.append(current)
                current = {"gsm": line.split(" = ", 1)[1]}
                descriptions, titles = [], []
            elif current is not None and line.startswith("!Sample_title = "):
                titles.append(line.split(" = ", 1)[1])
            elif current is not None and line.startswith("!Sample_description = "):
                descriptions.append(line.split(" = ", 1)[1])
            elif current is not None and line.startswith("!Sample_characteristics_ch1 = "):
                value = line.split(" = ", 1)[1]
                if ": " in value:
                    key, item = value.split(": ", 1)
                    current[key] = item
        if current is not None:
            current["description"] = descriptions[0] if descriptions else ""
            current["title"] = titles[0] if titles else ""
            records.append(current)
    metadata = pd.DataFrame(records)
    read_end = metadata["platform"].str.extract(r"Single Cell ([35])'")[0]
    metadata["orig.ident"] = (
        metadata["description"].str.replace("_", "-", regex=False) + "-" + read_end + "P"
    )
    return metadata.rename(
        columns={"tumor stage": "tumor_stage", "cancer subtype": "cancer_subtype"}
    )


def _read_10x_features(tf: tarfile.TarFile, members: dict, sample: str) -> list[str]:
    name = next(n for n in members if f"_{sample}_features" in n or f"_{sample}_genes" in n)
    genes = []
    with gzip.GzipFile(fileobj=tf.extractfile(members[name])) as f:
        for line in f:
            p = line.decode().strip().split("\t")
            genes.append((p[1] if len(p) > 1 else p[0]).upper())
    return genes


def _extract_mtx(tf: tarfile.TarFile, members: dict, sample: str, dest: Path) -> Path:
    name = next(n for n in members if f"_{sample}_matrix.mtx" in n)
    dest.mkdir(parents=True, exist_ok=True)
    out = dest / f"{sample}_matrix.mtx"
    if out.exists() and out.stat().st_size > 1_000_000:
        return out
    with gzip.GzipFile(fileobj=tf.extractfile(members[name])) as src, out.open("wb") as fh:
        while True:
            chunk = src.read(8 * 1024 * 1024)
            if not chunk:
                break
            fh.write(chunk)
    return out


def build_gse189357(tar_path: Path, scratch: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    units = pd.read_csv(DATA / "GSE189357_marker_units.tsv", sep="\t")
    keep = set(units.loc[units["eligible"].astype(str).str.lower() == "true", "patient"].astype(str))
    pb: dict[str, pd.Series] = {}
    meta_rows = []
    with tarfile.open(tar_path) as tf:
        members = {m.name: m for m in tf.getmembers()}
        for sample in [f"TD{i}" for i in range(1, 10)]:
            if sample not in keep:
                print(f"  skip {sample} (not eligible)", flush=True)
                continue
            print(f"  reading {sample}", flush=True)
            genes = _read_10x_features(tf, members, sample)
            mtx_path = _extract_mtx(tf, members, sample, scratch)
            mat = mmread(mtx_path).tocsc()
            if mat.shape[0] != len(genes):
                raise RuntimeError(f"{sample}: features {len(genes)} != mtx rows {mat.shape[0]}")
            idx = _upper_map(genes)
            panel = ["PTPRC"] + EPI
            keep_idx = [idx[g] for g in panel if g in idx]
            keep_names = [g for g in panel if g in idx]
            cell_panel = np.zeros((mat.shape[1], len(keep_names)), dtype=float)
            for j, gi in enumerate(keep_idx):
                cell_panel[:, j] = np.asarray(mat[gi, :].todense()).ravel()
            mal = _mal_mask(cell_panel, {g: i for i, g in enumerate(keep_names)})
            n_mal = int(mal.sum())
            if n_mal == 0:
                print(f"  WARN {sample}: 0 marker-malignant cells", flush=True)
                continue
            sums_raw = np.asarray(mat[:, mal].sum(axis=1)).ravel()
            ser = pd.Series(sums_raw, index=genes, dtype=float).groupby(level=0).sum()
            pb[sample] = ser
            meta_rows.append(
                {
                    "patient": sample,
                    "cohort": "GSE189357",
                    "unit": "patient",
                    "file": f"{sample}_matrix.mtx",
                    "n_cells": int(mat.shape[1]),
                    "n_malignant_summed": n_mal,
                    "n_genes": int(len(ser)),
                    "libsize": float(ser.sum()),
                }
            )
            print(f"  {sample}: cells={mat.shape[1]} mal={n_mal} lib={ser.sum():.0f}", flush=True)
            del mat
    if not pb:
        raise SystemExit("GSE189357: no patient pseudobulks")
    genes = sorted(set().union(*[set(s.index) for s in pb.values()]))
    counts = pd.DataFrame({k: v.reindex(genes).fillna(0.0) for k, v in pb.items()}, index=genes)
    return counts, pd.DataFrame(meta_rows)


def build_gse205335(ident_path: Path, matrix_path: Path, soft_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    locked = pd.read_csv(DATA / "GSE205335_patients.tsv", sep="\t")
    keep = set(locked["patient"].astype(str))
    identities = pd.read_csv(ident_path, sep="\t")
    if identities["barcode"].duplicated().any():
        raise ValueError("GSE205335 identity barcodes are not unique")
    metadata = parse_geo_soft(soft_path)

    print("read GSE205335 RDS (author-malignant UMI-sum)", flush=True)
    import rdata

    with tempfile.TemporaryDirectory(prefix="gse205335-") as tmp:
        rds_path = matrix_path
        if matrix_path.suffix == ".gz":
            rds_path = Path(tmp) / matrix_path.stem
            print(f"decompress {matrix_path.name}", flush=True)
            with gzip.open(matrix_path, "rb") as source, rds_path.open("wb") as dest:
                shutil.copyfileobj(source, dest, 16 * 1024 * 1024)
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", message='Missing constructor for R class "dgCMatrix"'
            )
            obj = rdata.read_rds(rds_path)

    expected = {"i", "p", "Dim", "Dimnames", "x"}
    if not expected.issubset(vars(obj)):
        raise TypeError("RDS is not the expected Matrix::dgCMatrix object")
    genes = np.asarray(obj.Dimnames[0], dtype=str)
    barcodes = np.asarray(obj.Dimnames[1], dtype=str)
    matrix = sparse.csc_matrix((obj.x, obj.i, obj.p), shape=tuple(obj.Dim), copy=False)
    print(f"  matrix {matrix.shape[0]} genes x {matrix.shape[1]} cells", flush=True)

    indexed = identities.set_index("barcode")
    missing = pd.Index(barcodes).difference(indexed.index)
    extra = indexed.index.difference(pd.Index(barcodes))
    if len(missing) or len(extra):
        raise ValueError(
            f"GSE205335 matrix/identity mismatch: {len(missing)} missing, {len(extra)} extra"
        )
    cells = indexed.loc[barcodes].reset_index()
    cells = cells.merge(
        metadata[["orig.ident", "gsm", "patient", "tissue", "recist", "cancer_subtype"]],
        on="orig.ident",
        how="left",
        validate="many_to_one",
    )
    if cells["patient"].isna().any():
        raise ValueError("GSE205335 identity samples did not match GEO metadata")
    cells["patient"] = cells["patient"].astype(str)
    mal = cells["lineage.sub"].eq("Malignant cells").to_numpy()
    print(f"  author malignant cells: {int(mal.sum())}", flush=True)

    gene_names = pd.Index(pd.Index(genes).astype(str).str.upper())
    pb: dict[str, pd.Series] = {}
    meta_rows = []
    for patient in sorted(keep):
        idx = np.flatnonzero(mal & cells["patient"].eq(patient).to_numpy())
        n_mal = int(idx.size)
        n_cells = int(cells["patient"].eq(patient).sum())
        if n_mal == 0:
            print(f"  skip {patient}: 0 author-malignant (should already be out of n=22)", flush=True)
            continue
        sums = np.asarray(matrix[:, idx].sum(axis=1)).ravel()
        ser = pd.Series(sums, index=gene_names, dtype=float).groupby(level=0).sum()
        pb[patient] = ser
        meta_rows.append(
            {
                "patient": patient,
                "cohort": "GSE205335",
                "unit": "patient",
                "file": "GSE205335_Lung_IO_UMI_matrix.rds.gz",
                "n_cells": n_cells,
                "n_malignant_summed": n_mal,
                "n_genes": int(len(ser)),
                "libsize": float(ser.sum()),
                "malig_def": "author_malig",
            }
        )
        print(f"  {patient}: cells={n_cells} mal={n_mal} lib={ser.sum():.0f}", flush=True)
    if not pb:
        raise SystemExit("GSE205335: no patient pseudobulks")
    all_genes = sorted(set().union(*[set(s.index) for s in pb.values()]))
    counts = pd.DataFrame({k: v.reindex(all_genes).fillna(0.0) for k, v in pb.items()}, index=all_genes)
    return counts, pd.DataFrame(meta_rows)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--geo", type=Path, default=DEFAULT_GEO)
    p.add_argument("--scratch", type=Path, default=DEFAULT_GEO / "extract")
    p.add_argument("--rebuild-189357", action="store_true")
    args = p.parse_args()
    DATA.mkdir(parents=True, exist_ok=True)

    out189 = DATA / "GSE189357_malignant_counts.tsv.gz"
    if args.rebuild_189357 or not out189.exists():
        tar189 = args.geo / "gse189357" / "GSE189357_RAW.tar"
        if not tar189.exists():
            raise SystemExit(f"missing {tar189}; run download.py (or keep committed counts)")
        print("Building GSE189357 patient malignant UMI-sum", flush=True)
        c189, m189 = build_gse189357(tar189, args.scratch)
        c189.to_csv(out189, sep="\t", compression="gzip")
        m189.to_csv(DATA / "GSE189357_malignant_meta.tsv", sep="\t", index=False)
        print(f"wrote GSE189357 {c189.shape[0]} genes x {c189.shape[1]} patients", flush=True)
    else:
        print(f"keep existing {out189}", flush=True)

    ident = args.geo / "gse205335" / "GSE205335_Lung_IO_CellIdentity.txt.gz"
    matrix = args.geo / "gse205335" / "GSE205335_Lung_IO_UMI_matrix.rds.gz"
    soft = args.geo / "gse205335" / "GSE205335_family.soft.gz"
    for path in (ident, matrix, soft):
        if not path.exists():
            raise SystemExit(f"missing {path}; run download.py")
    print("Building GSE205335 patient malignant UMI-sum", flush=True)
    c205, m205 = build_gse205335(ident, matrix, soft)
    c205.to_csv(DATA / "GSE205335_malignant_counts.tsv.gz", sep="\t", compression="gzip")
    m205.to_csv(DATA / "GSE205335_malignant_meta.tsv", sep="\t", index=False)
    print(f"wrote GSE205335 {c205.shape[0]} genes x {c205.shape[1]} patients", flush=True)


if __name__ == "__main__":
    main()
