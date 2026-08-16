#!/usr/bin/env python3
"""Cell-level supplement for GSE131907 malignant and GSE207422 epithelial cells.

Sample-pseudobulk GSEA is underpowered (n=14–21). This streams only the
pre-specified genes and scores TACSTD2 vs focal genes / signatures among cells.
Cells in one tumor are not independent; treat as a powered but clustered check.
"""
from __future__ import annotations

import gzip
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
OUT = ROOT / "results" / "hunt_tj_gsea"
FOCAL = ["CLDN1", "CLDN4", "CLDN7", "F11R", "PARD3"]


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def bh(p: pd.Series) -> pd.Series:
    valid = p.dropna()
    n = len(valid)
    if n == 0:
        return p
    order = valid.sort_values().index
    adj = valid.loc[order] * n / np.arange(1, n + 1)
    adj = adj[::-1].cummin()[::-1].clip(upper=1)
    q = p.copy()
    q.loc[order] = adj
    return q


def stream_genes(path: Path, keep_cells: list[str], genes: set[str]) -> pd.DataFrame:
    keep = set(keep_cells)
    found: dict[str, np.ndarray] = {}
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        cells = header[1:]
        idx = [i for i, c in enumerate(cells) if c in keep]
        names = [cells[i] for i in idx]
        log(f"  {path.name}: keeping {len(names)} / {len(cells)} cells, hunting {len(genes)} genes")
        for line in fh:
            gene = line.split("\t", 1)[0]
            if gene not in genes:
                continue
            parts = line.rstrip("\n").split("\t")
            vals = np.array([float(parts[i + 1]) for i in idx], dtype=np.float32)
            found[gene] = vals
            if len(found) == len(genes):
                break
    log(f"  found {len(found)} genes")
    return pd.DataFrame(found, index=names).T


def zmean(df: pd.DataFrame, genes: list[str]) -> pd.Series:
    present = [g for g in genes if g in df.index]
    if len(present) < 3:
        return pd.Series(np.nan, index=df.columns)
    sub = np.log1p(df.loc[present])
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, np.nan), axis=0)
    return z.mean(axis=0)


def correlate(x: pd.Series, mat: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for g in mat.index:
        if g == "TACSTD2":
            continue
        rho, p = stats.spearmanr(x, mat.loc[g], nan_policy="omit")
        rows.append({"gene": g, "rho": float(rho), "p": float(p)})
    out = pd.DataFrame(rows)
    out["q"] = bh(out["p"])
    return out


def main() -> int:
    payload = json.loads((ROOT / "data" / "genesets" / "primary_sets.json").read_text())
    sets = payload["sets"]
    want = set(["TACSTD2"] + FOCAL)
    for gs in sets.values():
        want.update(gs)

    # GSE131907 malignant, tumor/met sites
    ann = pd.read_csv(RAW / "geo" / "GSE131907_Lung_Cancer_cell_annotation.txt.gz", sep="\t")
    mal = ann[(ann["Cell_subtype"] == "Malignant cells") & (~ann["Sample_Origin"].isin(["nLung", "nLN"]))]
    log(f"GSE131907 malignant tumor/met cells {len(mal)}")
    m131 = stream_genes(
        RAW / "geo" / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz",
        mal["Index"].tolist(),
        want,
    )
    x = np.log1p(m131.loc["TACSTD2"])
    foc = correlate(x, np.log1p(m131.reindex(FOCAL).dropna(how="all")))
    foc["cohort"] = "GSE131907_malignant_cells"
    sig_rows = []
    for term, genes in sets.items():
        sc = zmean(m131, genes)
        rho, p = stats.spearmanr(x, sc, nan_policy="omit")
        sig_rows.append({"term": term, "rho": float(rho), "p": float(p), "cohort": "GSE131907_malignant_cells"})
    sig = pd.DataFrame(sig_rows)
    sig["q"] = bh(sig["p"])

    # GSE207422 epithelial gate (same rule as main hunt)
    path = RAW / "geo" / "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"
    markers = {"EPCAM", "KRT8", "KRT18", "KRT19", "PTPRC", "PECAM1", "COL1A1", "COL3A1"}
    cells, mark = _stream_markers(path, markers)
    epi_score = mark.reindex(["EPCAM", "KRT8", "KRT18", "KRT19"]).fillna(0).sum(axis=0)
    ptprc = mark.reindex(["PTPRC"]).fillna(0).iloc[0]
    pecam = mark.reindex(["PECAM1"]).fillna(0).iloc[0]
    col = mark.reindex(["COL1A1", "COL3A1"]).fillna(0).sum(axis=0)
    epi_cells = [c for c, keep in zip(cells, ((epi_score >= 1) & (ptprc == 0) & (pecam == 0) & (col == 0)).to_numpy()) if keep]
    log(f"GSE207422 epithelial cells {len(epi_cells)}")
    m207 = stream_genes(path, epi_cells, want)
    x2 = np.log1p(m207.loc["TACSTD2"])
    foc2 = correlate(x2, np.log1p(m207.reindex(FOCAL).dropna(how="all")))
    foc2["cohort"] = "GSE207422_epithelial_cells"
    sig_rows2 = []
    for term, genes in sets.items():
        sc = zmean(m207, genes)
        rho, p = stats.spearmanr(x2, sc, nan_policy="omit")
        sig_rows2.append({"term": term, "rho": float(rho), "p": float(p), "cohort": "GSE207422_epithelial_cells"})
    sig2 = pd.DataFrame(sig_rows2)
    sig2["q"] = bh(sig2["p"])

    foc_all = pd.concat([foc, foc2], ignore_index=True)
    sig_all = pd.concat([sig, sig2], ignore_index=True)
    foc_all.to_csv(OUT / "celllevel_focal.tsv", sep="\t", index=False)
    sig_all.to_csv(OUT / "celllevel_signatures.tsv", sep="\t", index=False)
    log("wrote cell-level tables")
    print(foc_all.to_string(index=False))
    print(sig_all.to_string(index=False))
    return 0


def _stream_markers(path: Path, markers: set[str]) -> tuple[list[str], pd.DataFrame]:
    found: dict[str, np.ndarray] = {}
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        cells = header[1:]
        for line in fh:
            gene = line.split("\t", 1)[0]
            if gene not in markers:
                continue
            found[gene] = np.fromstring(line.split("\t", 1)[1], sep="\t", dtype=np.float32)
            if len(found) == len(markers):
                break
    return cells, pd.DataFrame(found, index=cells).T


if __name__ == "__main__":
    raise SystemExit(main())
