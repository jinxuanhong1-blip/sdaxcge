#!/usr/bin/env python3
"""Stream GSE241934 MTX for CLDN4 in author Epi cells. IIT and Real kept separate."""
from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd

DATA = Path("/tmp/gse241934")
OUT = Path(__file__).resolve().parents[1] / "inputs"
OUT.mkdir(parents=True, exist_ok=True)


def gene_index(features_path: Path, gene: str) -> int:
    with gzip.open(features_path, "rt") as f:
        for i, line in enumerate(f, start=1):
            parts = line.rstrip("\n").split("\t")
            if gene in parts[:2]:
                return i
    raise KeyError(gene)


def extract_gene_from_mtx(mtx_path: Path, gene_row: int, n_cells: int) -> np.ndarray:
    vals = np.zeros(n_cells, dtype=np.float32)
    with gzip.open(mtx_path, "rt") as f:
        for line in f:
            if line.startswith("%"):
                continue
            break
        for line in f:
            r, c, v = line.split()
            if int(r) == gene_row:
                vals[int(c) - 1] = float(v)
    return vals


def group_of(x) -> str:
    x = str(x)
    if x in {"non-MPR", "NMPR", "nonMPR"}:
        return "NMPR"
    if x in {"MPR", "pCR", "MPR (pCR)"}:
        return "MPR"
    return ""


def analyze(name: str, mtx: Path, features: Path, barcodes: Path, meta: Path) -> pd.DataFrame:
    print(f"=== {name} CLDN4 ===", flush=True)
    md = pd.read_csv(meta, sep="\t", dtype=str, low_memory=False)
    bc = pd.read_csv(barcodes, sep="\t", header=None)
    n_cells = len(bc)
    gi = gene_index(features, "CLDN4")
    print("CLDN4 mtx row", gi, "n_cells", n_cells, flush=True)
    cldn = extract_gene_from_mtx(mtx, gi, n_cells)
    print("CLDN4 nonzero", int((cldn > 0).sum()), "max", float(cldn.max()), flush=True)
    if len(md) == n_cells:
        md = md.reset_index(drop=True)
        md["CLDN4"] = cldn
    else:
        raise RuntimeError(f"{name}: meta {len(md)} vs barcodes {n_cells}")
    major = md["major.cell.type"] if "major.cell.type" in md.columns else md["major_cell_type"]
    md["is_epi"] = major.eq("Epi")
    md["group"] = md["Pathological Response"].map(group_of)
    rows = []
    for sid, sub in md.groupby(md["sampleID"].astype(str)):
        n_epi = int(sub["is_epi"].sum())
        if n_epi > 0:
            mean_log1p = float(np.log1p(sub.loc[sub["is_epi"], "CLDN4"].astype(float)).mean())
            frac_pos = float((sub.loc[sub["is_epi"], "CLDN4"] > 0).mean())
        else:
            mean_log1p = np.nan
            frac_pos = np.nan
        rows.append(
            {
                "cohort": name,
                "sampleID": sid,
                "n_epi": n_epi,
                "malignant_CLDN4_mean_log1p": mean_log1p,
                "malignant_CLDN4_frac_pos": frac_pos,
                "group": sub["group"].iloc[0],
                "pathologic_response": sub["Pathological Response"].iloc[0],
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    needed = [
        DATA / "GSE241934_IIT_Matrix.mtx.gz",
        DATA / "GSE241934_IIT_features.tsv.gz",
        DATA / "GSE241934_IIT_barcodes.tsv.gz",
        DATA / "GSE241934_IIT_Meta.txt.gz",
        DATA / "GSE241934_Real_Matrix.mtx.gz",
        DATA / "GSE241934_RWC_features.tsv.gz",
        DATA / "GSE241934_RWC_barcodes.tsv.gz",
        DATA / "GSE241934_Real_Meta.txt.gz",
    ]
    missing = [p for p in needed if not p.exists()]
    if missing:
        print("MISSING", [p.name for p in missing], flush=True)
        raise SystemExit(2)
    iit = analyze(
        "IIT_EGFRmut",
        DATA / "GSE241934_IIT_Matrix.mtx.gz",
        DATA / "GSE241934_IIT_features.tsv.gz",
        DATA / "GSE241934_IIT_barcodes.tsv.gz",
        DATA / "GSE241934_IIT_Meta.txt.gz",
    )
    real = analyze(
        "REAL_WT",
        DATA / "GSE241934_Real_Matrix.mtx.gz",
        DATA / "GSE241934_RWC_features.tsv.gz",
        DATA / "GSE241934_RWC_barcodes.tsv.gz",
        DATA / "GSE241934_Real_Meta.txt.gz",
    )
    per = pd.concat([iit, real], ignore_index=True)
    per.to_csv(OUT / "GSE241934_cldn4_per_sample.csv", index=False)
    print(per.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
