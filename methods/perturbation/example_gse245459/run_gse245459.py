#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_gse245459.py
================
Public TACSTD2 shRNA RNA-seq (verified GEO accession GSE245459).

Title  : Exploring the Mechanism of TACSTD2 Causing Cisplatin Resistance
Model  : SKOV3 human ovarian cancer (n=12 = 4 groups x 3 biological reps)
Design : shNC / sh (TACSTD2 KD)  x  vehicle / cisplatin (DDP)
Suppl  : GSE245459_fpkm.anno.txt.gz  (~14.0 MB FPKM; public FTP)

PRIMARY CONTRAST: sh vs shNC  (TACSTD2 knockdown, no cisplatin)
SECONDARY      : shDDP vs shNCDDP  (same KD under cisplatin)

CAVEAT: the GEO supplement is FPKM, not raw integer counts, so this is NOT a
DESeq2/edgeR run. With n=3 we report Welch t-tests on log2(FPKM+1) as an
exploratory ranking aid and state that FDR here is not a count-model FDR.
See playbook.md §4 and §8.

Opposite-gene question (reverse of the CLDN4-KO example):
    TACSTD2 KD → does CLDN4 go UP?  (test, do not assume)
"""
from __future__ import annotations
import gzip
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
RES = HERE / "results"
GMT = HERE.parent / "gene_sets" / "playbook_sets.gmt"
RES.mkdir(parents=True, exist_ok=True)
DATA.mkdir(parents=True, exist_ok=True)

SUPPL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE245nnn/GSE245459/suppl/GSE245459_fpkm.anno.txt.gz"
MATRIX = DATA / "GSE245459_fpkm.anno.txt.gz"

GROUPS = {
    "shNC": ["shNC1", "shNC2", "shNC3"],
    "sh": ["sh1", "sh2", "sh3"],
    "shNCDDP": ["shNCDDP1", "shNCDDP2", "shNCDDP3"],
    "shDDP": ["shDDP1", "shDDP2", "shDDP3"],
}
CONTRASTS = [
    ("sh_vs_shNC", "sh", "shNC"),           # primary: TACSTD2 KD
    ("shDDP_vs_shNCDDP", "shDDP", "shNCDDP"),
]
FOCUS = [
    "TACSTD2", "CLDN4", "EPCAM",
    "CLDN1", "CLDN3", "CLDN7",
    "TJP1", "TJP2", "OCLN", "F11R",
    "CDH1", "VIM", "SNAI1", "SNAI2", "ZEB1",
    "STAT1", "STAT2", "IRF1", "ISG15", "MX1",
]


def log(msg: str) -> None:
    print(msg, flush=True)


def ensure_matrix() -> None:
    if MATRIX.exists():
        log(f"[data] cached {MATRIX.name} ({MATRIX.stat().st_size:,} B)")
        return
    log(f"[data] downloading {SUPPL}")
    urllib.request.urlretrieve(SUPPL, MATRIX)


def load() -> pd.DataFrame:
    with gzip.open(MATRIX, "rt") as fh:
        df = pd.read_csv(fh, sep="\t", low_memory=False)
    sample_cols = [c for g in GROUPS.values() for c in g]
    df = df.rename(columns={"GeneName": "gene"})
    df[sample_cols] = df[sample_cols].apply(pd.to_numeric, errors="coerce")
    df = df.dropna(subset=["gene"])
    df = df.groupby("gene", as_index=False)[sample_cols].max()
    return df


def bh_fallback(p: np.ndarray) -> np.ndarray:
    """BH FDR without statsmodels (always available)."""
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    pv = p[ok]
    n = len(pv)
    order = np.argsort(pv)
    ranked = pv[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    tmp = np.empty(n)
    tmp[order] = q
    out[ok] = tmp
    return out


def run_gsea(fc: pd.DataFrame, tag: str):
    import gseapy as gp
    rnk = fc[["gene", "log2FC"]].dropna()
    rnk = rnk[np.isfinite(rnk["log2FC"])].drop_duplicates("gene")
    rnk = rnk.sort_values("log2FC", ascending=False)
    rnk_path = RES / f"{tag}.rnk"
    rnk.to_csv(rnk_path, sep="\t", header=False, index=False)
    if not GMT.exists():
        log("[gsea] missing playbook_sets.gmt — skip")
        return None
    pre = gp.prerank(
        rnk=rnk_path.as_posix(), gene_sets=GMT.as_posix(),
        min_size=5, max_size=1000, permutation_num=1000, seed=42,
        outdir=None, no_plot=True, verbose=False,
    )
    res = pre.res2d.copy()
    res.to_csv(RES / f"{tag}_gsea.tsv", sep="\t", index=False)
    return res


def main() -> None:
    ensure_matrix()
    df = load()
    log(f"[data] genes x samples = {df.shape[0]} x {df.shape[1]-1}")

    def _contrast(frame, treat, ctrl):
        t = np.log2(frame[treat].to_numpy(float) + 1.0)
        c = np.log2(frame[ctrl].to_numpy(float) + 1.0)
        lfc = t.mean(axis=1) - c.mean(axis=1)
        tstat, p = stats.ttest_ind(t, c, axis=1, equal_var=False, nan_policy="omit")
        out = pd.DataFrame({
            "gene": frame["gene"], "log2FC": lfc, "t": tstat, "pvalue": p,
            "mean_treat": frame[treat].mean(axis=1).to_numpy(),
            "mean_ctrl": frame[ctrl].mean(axis=1).to_numpy(),
        })
        out["padj"] = bh_fallback(out["pvalue"].to_numpy())
        return out.sort_values("log2FC", ascending=False).reset_index(drop=True)

    lines = ["=" * 72, "GSE245459  SKOV3 shTACSTD2  (public FPKM, n=3)", "=" * 72,
             "NOT DESeq2: FPKM + Welch t on log2(FPKM+1). Exploratory only.", ""]
    focus_rows = []
    for tag, treat_name, ctrl_name in CONTRASTS:
        treat, ctrl = GROUPS[treat_name], GROUPS[ctrl_name]
        fc = _contrast(df, treat, ctrl)
        fc.to_csv(RES / f"{tag}_de.tsv.gz", sep="\t", index=False)
        sub = fc[fc["gene"].isin(FOCUS)].copy()
        sub.insert(0, "contrast", tag)
        focus_rows.append(sub)
        lines.append(f"--- {tag}  ({treat_name} vs {ctrl_name}) ---")
        lines.append(f"{'gene':<9}{'log2FC':>9}{'padj':>12}{'mean_KD':>10}{'mean_ctrl':>11}")
        for _, r in sub.set_index("gene").reindex(FOCUS).dropna(how="all").iterrows():
            lines.append(f"{r.name:<9}{r['log2FC']:>9.2f}{r['padj']:>12.3g}"
                         f"{r['mean_treat']:>10.2f}{r['mean_ctrl']:>11.2f}")
        t2 = fc[fc["gene"] == "TACSTD2"]
        c4 = fc[fc["gene"] == "CLDN4"]
        if len(t2):
            r = t2.iloc[0]
            lines.append(f"[on-target] TACSTD2 log2FC={r['log2FC']:+.2f} padj={r['padj']:.3g}")
        if len(c4):
            r = c4.iloc[0]
            got = "UP" if r["log2FC"] > 0 else "DOWN"
            lines.append(f"[opposite ] CLDN4   log2FC={r['log2FC']:+.2f} padj={r['padj']:.3g} -> {got} "
                         f"(hypothesis was UP on TACSTD2 loss)")
        try:
            gsea = run_gsea(fc, tag)
            if gsea is not None:
                lines.append("GSEA prerank (positive NES = up in KD):")
                for _, r in gsea.iterrows():
                    lines.append(f"  {r['Term']}: NES={r['NES']:.2f}  FDR={r['FDR q-val']:.3g}")
        except Exception as e:
            lines.append(f"[gsea] skipped ({e!r})")
        lines.append("")

    pd.concat(focus_rows, ignore_index=True).to_csv(RES / "focus_genes.tsv", sep="\t", index=False)
    text = "\n".join(lines)
    (RES / "run_summary.txt").write_text(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
