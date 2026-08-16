#!/usr/bin/env python3
"""scRNA-seq pseudobulk + single-cell analysis of Tacstd2 / Cldn4 vs anti-PD-1.

Cohorts (both mouse lung, 10x):
  GSE129297  Ctrl / PD1 / YKL(CDK7i) / combo   (SCLC lung, unfiltered droplet matrices)
  GSE133604  Ctrl / Ctrl+PD1 / Asf1a-KO / KO+PD1 (KP lung, cell-filtered matrices)

For each condition we compute, from real counts:
  * n cells, pseudobulk log2 CP10K of Tacstd2 & Cldn4
  * % cells expressing each gene
  * per-cell Spearman of the gene vs a cytotoxic/immune-module score (single-cell immune axis)
Output only under results/mouse/.
"""
import gzip
import json
import os

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
DATA = os.path.join(ROOT, "notes", "mouse", "data")
RES = os.path.join(ROOT, "results", "mouse")

GMAP = json.load(open(os.path.join(ROOT, "notes", "mouse", "gene_map.json")))
SYM2ENS = GMAP["symbol_to_ensembl"]
TARGETS = ["Tacstd2", "Cldn4"]
IMMUNE = ["Cd8a", "Cd8b1", "Gzmb", "Gzmk", "Prf1", "Ifng", "Nkg7",
          "Cxcl9", "Cxcl10", "Cd3e", "Pdcd1", "Ptprc"]

# Each sample entry: (matrix_file, features_file, condition)
COHORTS = {
    "GSE129297": dict(dir="GSE129297", min_umi=1000,
                      model="SCLC lung (RPM); anti-PD-1",
                      samples=[
                          ("GSM3704194_Ctrl.matrix.mtx.gz", "../GSE129297_features.tsv.gz", "Ctrl"),
                          ("GSM3704195_PD1.matrix.mtx.gz", "../GSE129297_features.tsv.gz", "antiPD1"),
                          ("GSM3704196_YKL.matrix.mtx.gz", "../GSE129297_features.tsv.gz", "CDK7i_YKL"),
                          ("GSM3704197_combo.matrix.mtx.gz", "../GSE129297_features.tsv.gz", "antiPD1_CDK7i"),
                      ]),
    "GSE133604": dict(dir="GSE133604", min_umi=500,
                      model="KP lung (KrasG12D;p53-/-); anti-PD-1",
                      samples=[
                          ("GSM3912860_1_Ctrl_matrix.mtx.gz", "../GSE133604_genes.tsv.gz", "Ctrl"),
                          ("GSM3912861_2_CtrlplusPD1_matrix.mtx.gz", "../GSE133604_genes.tsv.gz", "antiPD1"),
                          ("GSM3912862_3_ko_matrix.mtx.gz", "../GSE133604_genes.tsv.gz", "Asf1a_KO"),
                          ("GSM3912863_4_KOplusPD1_matrix.mtx.gz", "../GSE133604_genes.tsv.gz", "Asf1aKO_antiPD1"),
                      ]),
    "GSE297632": dict(dir="GSE297632", min_umi=500,
                      model="LLC subcutaneous; anti-PD-1 tolerant (scRNA companion of GSE297630)",
                      samples=[
                          ("GSM8995911_Control_GEX_matrix.mtx.gz", "GSM8995911_Control_GEX_features.tsv.gz", "Control"),
                          ("GSM8995912_PD_1_treatment_GEX_matrix.mtx.gz", "GSM8995912_PD_1_treatment_GEX_features.tsv.gz", "antiPD1"),
                      ]),
    "GSE222158": dict(dir="GSE222158", min_umi=500,
                      model="Oncogene-driven NSCLC; sorted CD45+/CD3+ immune cells; DC-CCL21 + anti-PD-1",
                      samples=[
                          ("GSM6915764_A_matrix.mtx.gz", "GSM6915764_A_features.tsv.gz", "CD45_Ctrl"),
                          ("GSM6915765_B_matrix.mtx.gz", "GSM6915765_B_features.tsv.gz", "CD45_antiPD1"),
                          ("GSM6915766_C_matrix.mtx.gz", "GSM6915766_C_features.tsv.gz", "CD45_DCvax"),
                          ("GSM6915767_D_matrix.mtx.gz", "GSM6915767_D_features.tsv.gz", "CD45_DCvax_antiPD1"),
                          ("GSM6915768_AT_matrix.mtx.gz", "GSM6915768_AT_features.tsv.gz", "CD3_Ctrl"),
                          ("GSM6915769_DT_matrix.mtx.gz", "GSM6915769_DT_features.tsv.gz", "CD3_DCvax_antiPD1"),
                      ]),
}


def load_feature_rows(path):
    """Return {ensembl_id: 1-based row index}."""
    rows = {}
    with gzip.open(path, "rt") as f:
        for i, line in enumerate(f, start=1):
            ens = line.split("\t")[0]
            rows[ens] = i
    return rows


def analyze_sample(mtx_path, gene_rows, min_umi):
    """gene_rows: dict symbol -> 1-based feature row. Returns per-cell dict."""
    df = pd.read_csv(mtx_path, sep=r"\s+", skiprows=3, header=None,
                     names=["gene", "cell", "count"], dtype=np.int32)
    ncell = int(df["cell"].max())
    col_total = np.bincount(df["cell"].values, weights=df["count"].values,
                            minlength=ncell + 1).astype(np.float64)
    cells = np.where(col_total >= min_umi)[0]
    cell_set = set(cells.tolist())
    # remap cell id -> position in filtered array
    pos = {c: i for i, c in enumerate(cells)}
    n = len(cells)
    tot = col_total[cells]

    want = {row: sym for sym, row in gene_rows.items()}
    sub = df[df["gene"].isin(want.keys())]
    # build gene x cell dense (few genes)
    mat = {sym: np.zeros(n, dtype=np.float64) for sym in gene_rows}
    for gene_row, cell_id, cnt in sub.itertuples(index=False):
        if cell_id in cell_set:
            mat[want[gene_row]][pos[cell_id]] = cnt
    return dict(n=n, tot=tot, mat=mat)


def main():
    pb_rows, sc_rows = [], []
    for cohort, cfg in COHORTS.items():
        cdir = os.path.join(DATA, cfg["dir"])
        for fname, feat, cond in cfg["samples"]:
            frows = load_feature_rows(os.path.join(cdir, feat))
            gene_rows = {}
            for sym in TARGETS + IMMUNE:
                ens = SYM2ENS.get(sym)
                if ens in frows:
                    gene_rows[sym] = frows[ens]
            res = analyze_sample(os.path.join(cdir, fname), gene_rows, cfg["min_umi"])
            n, tot, mat = res["n"], res["tot"], res["mat"]
            # log-normalized per-cell (CP10K)
            lognorm = {sym: np.log1p(mat[sym] / tot * 1e4) for sym in mat}
            # immune score
            imm_syms = [s for s in IMMUNE if s in mat]
            Z = []
            for s in imm_syms:
                v = lognorm[s]
                sd = v.std()
                Z.append((v - v.mean()) / sd if sd > 0 else np.zeros_like(v))
            score = np.mean(Z, axis=0) if Z else None
            for gene in TARGETS:
                if gene not in mat:
                    continue
                counts = mat[gene]
                pb_cpm = counts.sum() / tot.sum() * 1e6
                pct = float((counts > 0).mean() * 100)
                pb_rows.append(dict(dataset=cohort, model=cfg["model"], condition=cond,
                                    gene=gene, n_cells=n,
                                    pseudobulk_log2CPM=round(float(np.log2(pb_cpm + 1)), 4),
                                    pct_cells_expressing=round(pct, 2)))
                if score is not None:
                    rho, p = stats.spearmanr(lognorm[gene], score)
                    sc_rows.append(dict(dataset=cohort, condition=cond, gene=gene,
                                        n_cells=n, spearman_rho=round(float(rho), 4),
                                        spearman_p=float(p),
                                        immune_genes_used=";".join(imm_syms)))
            print(f"[{cohort}/{cond}] cells={n}")

    pb = pd.DataFrame(pb_rows)
    sc = pd.DataFrame(sc_rows)
    pb.to_csv(os.path.join(RES, "scrna_pseudobulk.csv"), index=False)
    sc.to_csv(os.path.join(RES, "scrna_singlecell_immune_corr.csv"), index=False)
    print("\n== pseudobulk ==")
    print(pb.to_string(index=False))
    print("\n== single-cell immune correlation ==")
    print(sc.to_string(index=False))


if __name__ == "__main__":
    main()
