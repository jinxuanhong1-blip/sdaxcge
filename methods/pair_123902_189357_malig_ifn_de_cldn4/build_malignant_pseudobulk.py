#!/usr/bin/env python3
"""UMI-sum marker-malignant cells to patient/donor pseudobulk.

Same gate as PR #459: (EPCAM|KRT8|KRT18|KRT19)>0 and PTPRC==0.
GSE123902: tumor/met only; donor is the unit (normals dropped).
GSE189357: one 10x sample per patient.
"""
from __future__ import annotations

import argparse
import gzip
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
EPI = ["EPCAM", "KRT8", "KRT18", "KRT19"]
DEFAULT_TAR = Path("/tmp/geo_pair_123902_189357")


def _upper_map(names: list[str]) -> dict[str, int]:
    return {n.upper(): i for i, n in enumerate(names)}


def _mal_mask(counts: np.ndarray, idx: dict[str, int]) -> np.ndarray:
    """counts: cells x genes (or cells x panel)."""
    def col(g: str) -> np.ndarray:
        if g not in idx:
            return np.zeros(counts.shape[0], dtype=float)
        return counts[:, idx[g]].astype(float)

    epi = np.zeros(counts.shape[0], dtype=bool)
    for g in EPI:
        epi |= col(g) > 0
    return epi & (col("PTPRC") == 0)


def build_gse123902(tar_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    units = pd.read_csv(DATA / "GSE123902_marker_units.tsv", sep="\t")
    tumor = units[units["tissue"].isin(["PRIMARY", "METASTASIS"])].copy()
    tumor = tumor.sort_values(["patient", "tissue"]).drop_duplicates("patient", keep="first")
    keep_files = set(tumor["file"].astype(str))
    keep_donors = set(tumor.loc[tumor["eligible"] == True, "patient"].astype(str))

    pb: dict[str, pd.Series] = {}
    meta_rows = []
    with tarfile.open(tar_path) as tf:
        for m in tf.getmembers():
            if not m.name.endswith(".csv.gz"):
                continue
            name = Path(m.name).name
            if name not in keep_files:
                print(f"  skip {name} (not locked tumor/met unit)", flush=True)
                continue
            parts = name.split("_")
            donor = parts[2] if len(parts) > 2 else name
            if donor not in keep_donors:
                print(f"  skip {name} (donor {donor} not eligible)", flush=True)
                continue
            print(f"  reading {name}", flush=True)
            raw = gzip.GzipFile(fileobj=tf.extractfile(m))
            df = pd.read_csv(raw, index_col=0)
            df.columns = [str(c).upper() for c in df.columns]
            if df.columns.duplicated().any():
                df = df.T.groupby(level=0).sum().T
            idx = {g: i for i, g in enumerate(df.columns)}
            mat = df.to_numpy(dtype=float)
            mal = _mal_mask(mat, idx)
            n_mal = int(mal.sum())
            if n_mal == 0:
                print(f"  WARN {donor}: 0 marker-malignant cells", flush=True)
                continue
            sums = pd.Series(mat[mal].sum(axis=0), index=df.columns, dtype=float)
            pb[donor] = sums
            meta_rows.append(
                {
                    "patient": donor,
                    "cohort": "GSE123902",
                    "unit": "donor",
                    "file": name,
                    "n_cells": int(mat.shape[0]),
                    "n_malignant_summed": n_mal,
                    "n_genes": int(len(sums)),
                    "libsize": float(sums.sum()),
                }
            )
            print(f"  {donor}: cells={mat.shape[0]} mal={n_mal} lib={sums.sum():.0f}", flush=True)
    if not pb:
        raise SystemExit("GSE123902: no donor pseudobulks")
    genes = sorted(set().union(*[set(s.index) for s in pb.values()]))
    counts = pd.DataFrame({k: v.reindex(genes).fillna(0.0) for k, v in pb.items()}, index=genes)
    return counts, pd.DataFrame(meta_rows)


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
    dest.parent.mkdir(parents=True, exist_ok=True)
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
    keep = set(units.loc[units["eligible"] == True, "patient"].astype(str))
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
            mat = mmread(mtx_path).tocsc()  # genes x cells
            if mat.shape[0] != len(genes):
                raise RuntimeError(f"{sample}: features {len(genes)} != mtx rows {mat.shape[0]}")
            idx = _upper_map(genes)
            # cells x panel for the gate
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


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--tars", type=Path, default=DEFAULT_TAR)
    p.add_argument("--scratch", type=Path, default=Path("/tmp/geo_pair_123902_189357/extract"))
    args = p.parse_args()
    DATA.mkdir(parents=True, exist_ok=True)

    p123 = args.tars / "GSE123902_RAW.tar"
    p189 = args.tars / "GSE189357_RAW.tar"
    if not p123.exists():
        raise SystemExit(f"missing {p123}; run download.py")
    if not p189.exists():
        raise SystemExit(f"missing {p189}; run download.py")

    print("Building GSE123902 donor malignant UMI-sum", flush=True)
    c123, m123 = build_gse123902(p123)
    c123.to_csv(DATA / "GSE123902_malignant_counts.tsv.gz", sep="\t", compression="gzip")
    m123.to_csv(DATA / "GSE123902_malignant_meta.tsv", sep="\t", index=False)
    print(f"wrote GSE123902 {c123.shape[0]} genes x {c123.shape[1]} donors", flush=True)

    print("Building GSE189357 patient malignant UMI-sum", flush=True)
    c189, m189 = build_gse189357(p189, args.scratch)
    c189.to_csv(DATA / "GSE189357_malignant_counts.tsv.gz", sep="\t", compression="gzip")
    m189.to_csv(DATA / "GSE189357_malignant_meta.tsv", sep="\t", index=False)
    print(f"wrote GSE189357 {c189.shape[0]} genes x {c189.shape[1]} patients", flush=True)


if __name__ == "__main__":
    main()
