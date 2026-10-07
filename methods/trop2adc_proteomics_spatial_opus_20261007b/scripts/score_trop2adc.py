#!/usr/bin/env python3
"""Score tumor-cell immune signatures in TROP2-ADC treated-vs-control RNA-seq.

Qualifying datasets (HARD CONSTRAINT 1: TROP2-ADC only):
  GSE304294  KYSE30 ESCC cells        IMMU132 (sacituzumab govitecan) vs Control
  GSE312098  CX-1 CRC cells           IMMU132 vs Control
  GSE311016  CRC patient-derived xeno IMMU132 vs Control (paired by PDX model)

IACS010759 (OXPHOS inhibitor) and GSK2606414 (PERK inhibitor) are side-by-side
internal controls in the SAME experiments. They are NOT TROP2-ADC and are scored
ONLY as specificity context, never as TROP2-ADC evidence. Combination arms contain
IMMU132 but are confounded by the partner drug and are reported as secondary.

Method: log2(x+1) transform of FPKM; log2FC = mean_log2(treated) - mean_log2(control);
per-gene Wilcoxon rank-sum (Mann-Whitney U); per-sample signature score = mean
log2(x+1) over detected genes in the set, compared treated-vs-control with Welch t
(and paired t where the design is paired). No fabricated numbers; null/opposite
results reported as-is.
"""
import gzip, io, os, json, sys
import numpy as np
import pandas as pd
from scipy import stats

DATA = os.path.join(os.path.dirname(__file__), "..", "data")
OUT = os.path.join(os.path.dirname(__file__), "..")

# ---- Signature gene sets (with aliases to match Ensembl/HGNC gene_name col) ----
SIGNATURES = {
    "APM": ["B2M", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F",
            "TAP1", "TAP2", "PSMB8", "PSMB9", "NLRC5", "TAPBP"],
    "IFN": ["STAT1", "IRF1", "CXCL9", "CXCL10", "CXCL11", "GBP1", "GBP2", "GBP4",
            "IRF7", "ISG15", "MX1", "OAS1", "IFIT1", "IFIT3", "IDO1"],
    "STING": ["CGAS", "STING1", "TBK1", "IKBKE", "IRF3", "IFNB1"],
    "SANITY_DOWN": ["CLDN4"],
    "TARGET": ["TACSTD2"],
}
# expected direction of change on TROP2-ADC treatment (per user hypothesis)
EXPECTED = {"APM": "up", "IFN": "up", "STING": "up",
            "SANITY_DOWN": "down", "TARGET": "sanity"}

ALIASES = {
    "CGAS": ["CGAS", "MB21D1", "C6orf150"],
    "STING1": ["STING1", "TMEM173"],
    "NLRC5": ["NLRC5"],
    "TAPBP": ["TAPBP"],
    "ISG15": ["ISG15", "G1P2"],
    "MX1": ["MX1"],
    "IFIT1": ["IFIT1"],
    "IFIT3": ["IFIT3"],
}


def read_fpkm(path):
    """Read a GEO gene_fpkm.txt.gz with unknown encoding (utf-8/utf-16)."""
    with gzip.open(path, "rb") as fh:
        raw = fh.read()
    for enc in ("utf-8", "utf-16", "utf-16-le", "latin-1"):
        try:
            txt = raw.decode(enc)
            if "gene_id" in txt[:200] and "gene_name" in txt[:4000]:
                return pd.read_csv(io.StringIO(txt), sep="\t")
        except Exception:
            continue
    # last resort
    return pd.read_csv(io.StringIO(raw.decode("latin-1")), sep="\t")


def expr_matrix(df, sample_cols):
    """Collapse to gene_name x samples (max FPKM across duplicate symbols)."""
    sub = df[["gene_name"] + sample_cols].copy()
    for c in sample_cols:
        sub[c] = pd.to_numeric(sub[c], errors="coerce")
    sub = sub.dropna(subset=["gene_name"])
    g = sub.groupby("gene_name")[sample_cols].max()
    return g


def gene_candidates(gene):
    return ALIASES.get(gene, [gene])


def score_dataset(name, df, ctrl_cols, trt_cols, paired=False):
    mat = expr_matrix(df, ctrl_cols + trt_cols)
    log2 = np.log2(mat + 1.0)

    per_gene_rows = []
    sig_present = {s: [] for s in SIGNATURES}
    for sig, genes in SIGNATURES.items():
        for gene in genes:
            found = None
            for cand in gene_candidates(gene):
                if cand in log2.index:
                    found = cand
                    break
            if found is None:
                per_gene_rows.append(dict(dataset=name, signature=sig, gene=gene,
                                          detected=False, mean_log2_ctrl=np.nan,
                                          mean_log2_trt=np.nan, log2FC=np.nan,
                                          wilcoxon_p=np.nan))
                continue
            c = log2.loc[found, ctrl_cols].astype(float).values
            t = log2.loc[found, trt_cols].astype(float).values
            log2fc = float(np.mean(t) - np.mean(c))
            try:
                _, wp = stats.mannwhitneyu(t, c, alternative="two-sided")
            except ValueError:
                wp = np.nan
            per_gene_rows.append(dict(dataset=name, signature=sig, gene=gene,
                                      detected=True, matched_symbol=found,
                                      mean_log2_ctrl=float(np.mean(c)),
                                      mean_log2_trt=float(np.mean(t)),
                                      log2FC=log2fc, wilcoxon_p=float(wp) if wp==wp else np.nan))
            sig_present[sig].append(found)

    # per-sample signature scores + t-test
    sig_rows = []
    for sig, genes in SIGNATURES.items():
        present = sig_present[sig]
        if not present:
            continue
        ctrl_score = log2.loc[present, ctrl_cols].mean(axis=0).values
        trt_score = log2.loc[present, trt_cols].mean(axis=0).values
        delta = float(np.mean(trt_score) - np.mean(ctrl_score))
        tt = stats.ttest_ind(trt_score, ctrl_score, equal_var=False)
        row = dict(dataset=name, signature=sig, n_genes=len(present),
                   expected=EXPECTED[sig],
                   score_ctrl_mean=float(np.mean(ctrl_score)),
                   score_trt_mean=float(np.mean(trt_score)),
                   delta_log2=delta, welch_t=float(tt.statistic),
                   welch_p=float(tt.pvalue))
        if paired and len(ctrl_score) == len(trt_score):
            pt = stats.ttest_rel(trt_score, ctrl_score)
            row["paired_t"] = float(pt.statistic)
            row["paired_p"] = float(pt.pvalue)
        sig_rows.append(row)

    return pd.DataFrame(per_gene_rows), pd.DataFrame(sig_rows)


def main():
    datasets = []

    # ---------- GSE304294 (ESCC KYSE30) ----------
    df = read_fpkm(os.path.join(DATA, "GSE304294_gene_fpkm.txt.gz"))
    # Library names: OX1=Control(1-3), OX2=IMMU132(1-2), OX3=IACS010759(1-3), OX4=Combination(1-3)
    ds = dict(name="GSE304294",
              ctrl=["OX1_1", "OX1_2", "OX1_3"],
              trt=["OX2_1", "OX2_2"],
              paired=False,
              context={"IACS010759(OXPHOSi,specificity-ctrl)": ["OX3_1", "OX3_2", "OX3_3"],
                       "Combination(IMMU132+IACS)": ["OX4_1", "OX4_2", "OX4_3"]})
    datasets.append((df, ds))

    # ---------- GSE312098 (CRC CX-1) ----------
    df = read_fpkm(os.path.join(DATA, "GSE312098_gene_fpkm.txt.gz"))
    df.columns = [c.strip() for c in df.columns]
    ds = dict(name="GSE312098",
              ctrl=["X_1", "X_2", "X_3"],
              trt=["X_4", "X_5", "X_6"],
              paired=False,
              context={"GSK2606414(PERKi,specificity-ctrl)": ["X_7", "X_8", "X_9"],
                       "Combination(IMMU132+GSK)": ["X_10", "X_11", "X_12"]})
    datasets.append((df, ds))

    # ---------- GSE311016 (CRC PDX, paired) ----------
    df = read_fpkm(os.path.join(DATA, "GSE311016_gene_fpkm.txt.gz"))
    df.columns = [c.strip() for c in df.columns]
    ds = dict(name="GSE311016",
              ctrl=["C_114", "C_36", "C_82", "C_83", "C_196"],
              trt=["T_114", "T_36", "T_82", "T_83", "T_196"],
              paired=True, context={})
    datasets.append((df, ds))

    all_gene, all_sig, all_ctx = [], [], []
    for df, ds in datasets:
        pg, ps = score_dataset(ds["name"], df, ds["ctrl"], ds["trt"], ds["paired"])
        all_gene.append(pg); all_sig.append(ps)
        print(f"\n===== {ds['name']}  IMMU132(n={len(ds['trt'])}) vs Control(n={len(ds['ctrl'])})"
              f"{' [paired]' if ds['paired'] else ''} =====")
        print(ps[["signature", "n_genes", "expected", "delta_log2", "welch_p"] +
                 (["paired_p"] if "paired_p" in ps.columns else [])].to_string(index=False))

        # specificity-context arms (NOT scored as TROP2 evidence)
        for label, cols in ds.get("context", {}).items():
            _, cs = score_dataset(ds["name"], df, ds["ctrl"], cols, False)
            cs.insert(1, "context_arm", label)
            all_ctx.append(cs)

    gene_df = pd.concat(all_gene, ignore_index=True)
    sig_df = pd.concat(all_sig, ignore_index=True)
    gene_df.to_csv(os.path.join(OUT, "per_gene_log2fc.csv"), index=False)
    sig_df.to_csv(os.path.join(OUT, "signature_scores.csv"), index=False)
    if all_ctx:
        pd.concat(all_ctx, ignore_index=True).to_csv(
            os.path.join(OUT, "specificity_context_arms.csv"), index=False)
    print("\nWrote per_gene_log2fc.csv, signature_scores.csv, specificity_context_arms.csv")


if __name__ == "__main__":
    main()
