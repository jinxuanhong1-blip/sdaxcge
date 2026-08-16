#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_gse207704.py
================
Worked example for the KD/KO RNA-seq playbook.

Dataset : GSE207704  (REAL, verified GEO accession)
Title   : "Claudin-4-adhesion signaling drives breast cancer metabolism and
           progression via liver X receptor beta"
Design  : CLDN4 CRISPR knockout (KO) vs wild-type (WT) in two human breast
          cancer cell lines (MCF7, T47D). Human (Homo sapiens).
Suppl   : GSE207704_CLDN4_RNAseq.txt.gz  (~1.0 MB, Cufflinks FPKM matrix)

IMPORTANT CAVEAT (read the playbook section "n=2-3 pitfalls"):
    The *processed* supplementary matrix distributed on GEO for this series is
    an FPKM table with exactly ONE column per condition per cell line
    (MCF7_WT, MCF7_CLDN4KO, T47D_WT, T47D_CLDN4KO). There are NO within-group
    replicates in the processed file, so we CANNOT compute a valid DESeq2/edgeR
    p-value from it. A dispersion-based negative-binomial test needs raw integer
    counts AND >= 2 replicates per group.

    Therefore this script demonstrates the parts of the playbook that a
    replicate-free FPKM matrix legitimately supports:
        1. per-cell-line log2 fold-change (KO vs WT)
        2. an "on-target" sanity check on CLDN4 itself
        3. the "opposite gene" check: does CLDN4 loss coincide with TACSTD2
           (TROP2) up-regulation, and what happens to other tight-junction /
           EMT markers?
        4. GSEA (pre-ranked) against Hallmark IFN-alpha, IFN-gamma, EMT and a
           KEGG Tight-junction set, ranking genes by the mean KO-vs-WT log2FC.

    To reproduce a full DESeq2/edgeR run you must start from raw counts with
    replicates (see ../templates/). This series does not ship those, so the
    proper route is re-quantifying the FASTQs from SRA (out of the <2 GB scope).

Outputs (written to results/):
    de_log2fc_table.tsv        full per-gene log2FC table
    focus_genes.tsv            CLDN4/TACSTD2/claudins/TJ/EMT focus panel
    gsea_prerank_report.tsv    GSEA pre-ranked enrichment results
    rank.rnk                   the ranking vector used for GSEA
    run_summary.txt            human-readable summary printed at the end
"""
from __future__ import annotations
import os
import sys
import gzip
import json
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
RES = HERE / "results"
GMT_DIR = HERE.parent / "gene_sets"
RES.mkdir(parents=True, exist_ok=True)
GMT_DIR.mkdir(parents=True, exist_ok=True)

GSE = "GSE207704"
SUPPL_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/"
    "GSE207704/suppl/GSE207704_CLDN4_RNAseq.txt.gz"
)
MATRIX = DATA / "GSE207704_CLDN4_RNAseq.txt.gz"

# Enrichr library endpoints (real gene sets, fetched once and cached as .gmt)
ENRICHR = "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName={lib}"
WANTED = {
    "MSigDB_Hallmark_2020": [
        "Interferon Alpha Response",
        "Interferon Gamma Response",
        "Epithelial Mesenchymal Transition",
    ],
    "KEGG_2021_Human": [
        "Tight junction",
    ],
}


def log(msg: str) -> None:
    print(msg, flush=True)


def ensure_matrix() -> None:
    if MATRIX.exists():
        log(f"[data] using cached {MATRIX.name} ({MATRIX.stat().st_size:,} bytes)")
        return
    DATA.mkdir(parents=True, exist_ok=True)
    log(f"[data] downloading {SUPPL_URL}")
    urllib.request.urlretrieve(SUPPL_URL, MATRIX)
    log(f"[data] saved {MATRIX} ({MATRIX.stat().st_size:,} bytes)")


def load_fpkm() -> pd.DataFrame:
    with gzip.open(MATRIX, "rt") as fh:
        df = pd.read_csv(fh, sep="\t")
    # Columns of interest carry a "(fpkm)" suffix in the header.
    fpkm_cols = {c: c.split(" ")[0] for c in df.columns if "(fpkm)" in c.lower()}
    df = df.rename(columns=fpkm_cols)
    keep = ["gene_short_name"] + list(fpkm_cols.values())
    df = df[keep].copy()
    df = df.rename(columns={"gene_short_name": "gene"})
    # Collapse duplicate gene symbols by max FPKM (isoform/locus level -> gene).
    num_cols = [c for c in df.columns if c != "gene"]
    df[num_cols] = df[num_cols].apply(pd.to_numeric, errors="coerce")
    df = df.dropna(subset=["gene"])
    df = df.groupby("gene", as_index=False)[num_cols].max()
    return df


def compute_log2fc(df: pd.DataFrame) -> pd.DataFrame:
    pc = 1.0  # pseudocount on FPKM before ratio
    out = pd.DataFrame({"gene": df["gene"]})
    out["MCF7_log2FC"] = np.log2((df["MCF7_CLDN4KO_FPKM"] + pc) / (df["MCF7_WT_FPKM"] + pc))
    out["T47D_log2FC"] = np.log2((df["T47D_CLDN4KO_FPKM"] + pc) / (df["T47D_WT_FPKM"] + pc))
    out["mean_log2FC"] = out[["MCF7_log2FC", "T47D_log2FC"]].mean(axis=1)
    # concordance: same sign in both cell lines
    out["concordant"] = np.sign(out["MCF7_log2FC"]) == np.sign(out["T47D_log2FC"])
    for c in ["MCF7", "T47D"]:
        out[f"{c}_WT_FPKM"] = df[f"{c}_WT_FPKM"].values
        out[f"{c}_CLDN4KO_FPKM"] = df[f"{c}_CLDN4KO_FPKM"].values
    return out.sort_values("mean_log2FC", ascending=False).reset_index(drop=True)


def fetch_gene_sets() -> Path:
    """Fetch requested gene sets from Enrichr, cache a combined .gmt, return path."""
    gmt = GMT_DIR / "playbook_sets.gmt"
    if gmt.exists():
        log(f"[gsea] using cached gene sets {gmt.name}")
        return gmt
    lines = []
    for lib, terms in WANTED.items():
        url = ENRICHR.format(lib=lib)
        log(f"[gsea] fetching {lib}")
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        text = urllib.request.urlopen(req, timeout=60).read().decode("utf-8", "replace")
        wanted_lower = {t.lower() for t in terms}
        for row in text.strip().split("\n"):
            parts = row.split("\t")
            name = parts[0].strip()
            base = name.split(" (")[0].strip().lower()
            if base in wanted_lower or name.strip().lower() in wanted_lower:
                genes = [g.split(",")[0].strip() for g in parts[2:] if g.strip()]
                tag = name.replace(" ", "_")
                lines.append("\t".join([f"{lib}::{tag}", ""] + genes))
    gmt.write_text("\n".join(lines) + "\n")
    log(f"[gsea] wrote {gmt} ({len(lines)} sets)")
    return gmt


def focus_report(fc: pd.DataFrame) -> pd.DataFrame:
    panel = [
        "CLDN4",            # the knocked-out gene (on-target check)
        "TACSTD2",          # TROP2 - the "opposite gene" hypothesis
        "EPCAM",            # TROP2 paralog / epithelial marker
        "CLDN1", "CLDN3", "CLDN7",   # other claudins (compensation?)
        "TJP1", "TJP2", "OCLN", "F11R",  # tight-junction scaffold
        "CDH1", "VIM", "SNAI1", "SNAI2", "ZEB1",  # EMT axis
        "STAT1", "STAT2", "IRF1", "ISG15", "MX1",  # IFN response
    ]
    sub = fc[fc["gene"].isin(panel)].copy()
    sub["gene"] = pd.Categorical(sub["gene"], categories=panel, ordered=True)
    return sub.sort_values("gene")


def run_gsea(fc: pd.DataFrame, gmt: Path):
    import gseapy as gp
    rnk = fc[["gene", "mean_log2FC"]].dropna()
    rnk = rnk[np.isfinite(rnk["mean_log2FC"])]
    rnk = rnk.drop_duplicates("gene").sort_values("mean_log2FC", ascending=False)
    rnk_path = RES / "rank.rnk"
    rnk.to_csv(rnk_path, sep="\t", header=False, index=False)
    log(f"[gsea] ranked {len(rnk)} genes -> {rnk_path.name}")
    pre = gp.prerank(
        rnk=rnk_path.as_posix(),
        gene_sets=gmt.as_posix(),
        min_size=5,
        max_size=1000,
        permutation_num=1000,
        seed=42,
        outdir=None,
        no_plot=True,
        verbose=False,
    )
    res = pre.res2d.copy()
    res.to_csv(RES / "gsea_prerank_report.tsv", sep="\t", index=False)
    return res


def main() -> None:
    ensure_matrix()
    df = load_fpkm()
    log(f"[data] genes x samples = {df.shape[0]} x {df.shape[1]-1}")
    fc = compute_log2fc(df)
    fc.to_csv(RES / "de_log2fc_table.tsv", sep="\t", index=False)

    focus = focus_report(fc)
    focus.to_csv(RES / "focus_genes.tsv", sep="\t", index=False)

    gmt = fetch_gene_sets()
    try:
        gsea = run_gsea(fc, gmt)
        gsea_ok = True
    except Exception as e:  # network / gseapy hiccup should not lose the DE table
        gsea = None
        gsea_ok = False
        log(f"[gsea] SKIPPED ({e!r})")

    # ---- human-readable summary -------------------------------------------
    lines = []
    def emit(s=""):
        lines.append(s)
        log(s)

    emit("=" * 72)
    emit(f"{GSE}  CLDN4 KO vs WT  (MCF7, T47D)  -  worked example")
    emit("=" * 72)
    emit("Processed matrix is FPKM with 1 sample/condition -> no p-values.")
    emit("Metric = mean(log2((KO+1)/(WT+1))) across MCF7 and T47D.")
    emit("")
    emit("On-target + opposite-gene focus panel:")
    emit(f"{'gene':<9}{'MCF7_l2fc':>11}{'T47D_l2fc':>11}{'mean':>9}  concordant")
    for _, r in focus.iterrows():
        emit(f"{r['gene']:<9}{r['MCF7_log2FC']:>11.2f}{r['T47D_log2FC']:>11.2f}"
             f"{r['mean_log2FC']:>9.2f}  {bool(r['concordant'])}")
    emit("")
    cldn4 = focus[focus["gene"] == "CLDN4"]
    trop2 = focus[focus["gene"] == "TACSTD2"]
    if len(cldn4):
        r = cldn4.iloc[0]
        emit(f"[on-target] CLDN4 mean log2FC = {r['mean_log2FC']:+.2f} "
             f"(expect strongly negative in a KO).")
    if len(trop2):
        r = trop2.iloc[0]
        direction = "UP" if r["mean_log2FC"] > 0 else "DOWN"
        emit(f"[opposite ] TACSTD2/TROP2 mean log2FC = {r['mean_log2FC']:+.2f} "
             f"-> {direction} on CLDN4 loss "
             f"(concordant={bool(r['concordant'])}).")
    emit("")
    if gsea_ok and gsea is not None:
        emit("GSEA pre-ranked (KO high -> positive NES):")
        cols = [c for c in ["Term", "NES", "NOM p-val", "FDR q-val"] if c in gsea.columns]
        emit("  " + "  ".join(cols))
        for _, r in gsea.iterrows():
            emit("  " + "  ".join(str(r[c]) for c in cols))
    else:
        emit("GSEA: skipped (see log).")
    emit("")
    emit("Files written to results/: de_log2fc_table.tsv, focus_genes.tsv,")
    emit("gsea_prerank_report.tsv, rank.rnk, run_summary.txt")

    (RES / "run_summary.txt").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
