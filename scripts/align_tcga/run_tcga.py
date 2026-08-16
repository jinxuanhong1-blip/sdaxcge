"""
USER-ALIGN: replicate the user's bulk finding in TCGA LUAD+LUSC (NSCLC).

Method (matching the user's description):
  1. Spearman correlation of TACSTD2 (TROP-2) vs immune / cytotoxic /
     exhaustion signatures.
  2. Partial Spearman correlation after controlling for tumour purity
     (ESTIMATE and ABSOLUTE; CPE as the Aran consensus, reported alongside).

Data (public, open access):
  * Expression: UCSC Xena TCGA LUAD & LUSC HiSeqV2 (RSEM, log2(norm+1)).
  * Purity: Aran et al. Nat Commun 2015, Supplementary Data 1
    (per-sample ESTIMATE, ABSOLUTE, CPE).

Outputs go to results/align_tcga/.
"""
from __future__ import annotations

import os
import pandas as pd

import common as C

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "..", "results", "align_tcga", "data")
OUT = os.path.join(HERE, "..", "..", "results", "align_tcga")

PURITY_METHODS = ["ESTIMATE", "ABSOLUTE", "CPE"]
SIG_FEATURES = list(C.SIGNATURES.keys())
IMMUNE_PREFIXES = tuple(SIG_FEATURES) + ("gene_",)


def load_purity() -> pd.DataFrame:
    xl = pd.ExcelFile(os.path.join(DATA, "Aran_CPE_purity.xlsx"))
    raw = xl.parse(xl.sheet_names[0], header=None)
    hdr = raw.index[raw.iloc[:, 0].astype(str).str.strip() == "Sample ID"][0]
    names = [str(v).strip() for v in raw.iloc[hdr].tolist()]
    df = raw.iloc[hdr + 1:].copy()
    df.columns = names
    df = df.rename(columns={"Sample ID": "SampleID",
                            "Cancer type": "CancerType"})
    df = df[df["SampleID"].astype(str).str.startswith("TCGA")].copy()
    df["barcode15"] = df["SampleID"].astype(str).str[:15]
    for col in PURITY_METHODS:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.drop_duplicates("barcode15")
    return df.set_index("barcode15")[["CancerType"] + PURITY_METHODS]


def primary_tumors(expr: pd.DataFrame) -> pd.DataFrame:
    cols = [c for c in expr.columns if str(c).endswith("-01")]
    return expr[cols]


def build_expression() -> tuple[pd.DataFrame, pd.Series]:
    luad = C.load_xena_matrix(os.path.join(DATA, "LUAD.HiSeqV2.gz"))
    lusc = C.load_xena_matrix(os.path.join(DATA, "LUSC.HiSeqV2.gz"))
    luad_t = primary_tumors(luad)
    lusc_t = primary_tumors(lusc)
    genes = luad_t.index.intersection(lusc_t.index)
    luad_t = luad_t.loc[genes]
    lusc_t = lusc_t.loc[genes]
    pooled = pd.concat([luad_t, lusc_t], axis=1)
    cohort = pd.Series(["LUAD"] * luad_t.shape[1] + ["LUSC"] * lusc_t.shape[1],
                       index=list(luad_t.columns) + list(lusc_t.columns))
    return pooled, cohort


def compute_features(expr: pd.DataFrame) -> pd.DataFrame:
    feats = {C.TARGET_GENE: expr.loc[C.TARGET_GENE].astype(float)}
    for name, genes in C.SIGNATURES.items():
        score, _ = C.signature_score(expr, genes)
        feats[name] = score
    avail = set(expr.index)
    for g in C.SINGLE_IMMUNE_GENES:
        r = C.resolve_symbol(g, avail)
        if r is not None:
            feats["gene_" + g] = expr.loc[r].astype(float)
    for g in C.TJ_GENES:
        r = C.resolve_symbol(g, avail)
        if r is not None:
            feats["TJ_" + g] = expr.loc[r].astype(float)
    return pd.DataFrame(feats)


def is_immune(feature: str) -> bool:
    return feature.startswith("gene_") or feature in SIG_FEATURES


def run_block(feat: pd.DataFrame, purity: pd.DataFrame, tag: str,
              purity_cols: list[str] | None = None) -> pd.DataFrame:
    """Spearman + partial-Spearman of TACSTD2 vs every other feature."""
    if purity_cols is None:
        purity_cols = PURITY_METHODS
    df = feat.join(purity[purity_cols], how="left")
    x = df[C.TARGET_GENE].to_numpy(float)
    rows = []
    targets = [c for c in feat.columns if c != C.TARGET_GENE]
    for col in targets:
        y = df[col].to_numpy(float)
        r, p, n = C.spearman(x, y)
        rec = {"cohort": tag, "feature": col, "class":
               "TJ" if col.startswith("TJ_") else "immune",
               "spearman_rho": r, "spearman_p": p, "n": n}
        for pm in purity_cols:
            z = df[pm].to_numpy(float)
            pr, pp, pn = C.partial_spearman(x, y, z)
            rec[f"partial_rho_{pm}"] = pr
            rec[f"partial_p_{pm}"] = pp
            rec[f"partial_n_{pm}"] = pn
        rows.append(rec)
    out = pd.DataFrame(rows)
    # BH-FDR within this block, immune features only (user's claim is about
    # immune/cytotoxic/exhaustion; TJ genes are a separate positive control).
    imm = out["class"] == "immune"
    for pcol, qcol in [("spearman_p", "spearman_q"),
                       ("partial_p_ESTIMATE", "partial_q_ESTIMATE"),
                       ("partial_p_ABSOLUTE", "partial_q_ABSOLUTE"),
                       ("partial_p_CPE", "partial_q_CPE")]:
        if pcol not in out.columns:
            continue
        q = [float("nan")] * len(out)
        idx = list(out.index[imm])
        qs = C.bh_fdr(out.loc[idx, pcol].tolist())
        for i, qi in zip(idx, qs):
            q[i] = qi
        out[qcol] = q
    return out


def tacstd2_vs_purity(feat: pd.DataFrame, purity: pd.DataFrame,
                      cohort: pd.Series) -> pd.DataFrame:
    rows = []
    frames = {"NSCLC_pooled": feat}
    for ck in ["LUAD", "LUSC"]:
        frames[ck] = feat.loc[cohort[cohort == ck].index]
    for tag, sub in frames.items():
        df = sub.join(purity[PURITY_METHODS], how="left")
        x = df[C.TARGET_GENE].to_numpy(float)
        for pm in PURITY_METHODS:
            r, p, n = C.spearman(x, df[pm].to_numpy(float))
            rows.append({"cohort": tag, "purity_method": pm,
                         "spearman_rho": r, "spearman_p": p, "n": n})
    return pd.DataFrame(rows)


def fmt_csv(df: pd.DataFrame) -> pd.DataFrame:
    """Round ρ to 4 d.p.; leave p/q in full float (scientific in writeup)."""
    out = df.copy()
    for c in out.columns:
        if c.endswith("_rho"):
            out[c] = out[c].round(4)
    return out


def write_direction_tally(res: pd.DataFrame, path: str) -> None:
    """Count how many immune features stay negative / significant."""
    rows = []
    for cohort, g in res.groupby("cohort", sort=False):
        imm = g[g["class"] == "immune"]
        rec = {"cohort": cohort, "n_immune_features": len(imm),
               "spearman_neg": int((imm.spearman_rho < 0).sum()),
               "spearman_neg_p05": int(((imm.spearman_rho < 0)
                                        & (imm.spearman_p < 0.05)).sum()),
               "spearman_neg_q05": int(((imm.spearman_rho < 0)
                                        & (imm.spearman_q < 0.05)).sum())}
        for pm in PURITY_METHODS:
            rho, p, q = (f"partial_rho_{pm}", f"partial_p_{pm}",
                         f"partial_q_{pm}")
            if rho not in imm.columns:
                continue
            rec[f"{pm}_neg"] = int((imm[rho] < 0).sum())
            rec[f"{pm}_neg_p05"] = int(((imm[rho] < 0)
                                        & (imm[p] < 0.05)).sum())
            rec[f"{pm}_neg_q05"] = int(((imm[rho] < 0)
                                        & (imm[q] < 0.05)).sum())
        rows.append(rec)
    pd.DataFrame(rows).to_csv(path, index=False)


def main() -> None:
    pooled, cohort = build_expression()
    purity = load_purity()
    feat = compute_features(pooled)
    feat["cohort"] = cohort

    common_ids = feat.index.intersection(purity.index)
    n_cpe = int(purity.loc[common_ids, "CPE"].notna().sum())
    n_est = int(purity.loc[common_ids, "ESTIMATE"].notna().sum())
    n_abs = int(purity.loc[common_ids, "ABSOLUTE"].notna().sum())
    both_mask = (purity.loc[common_ids, "ESTIMATE"].notna()
                 & purity.loc[common_ids, "ABSOLUTE"].notna())
    n_both = int(both_mask.sum())
    print(f"[TCGA] pooled primary tumours: {pooled.shape[1]} "
          f"(LUAD={int((cohort=='LUAD').sum())}, "
          f"LUSC={int((cohort=='LUSC').sum())})")
    print(f"[TCGA] purity overlap: CPE={n_cpe} ESTIMATE={n_est} "
          f"ABSOLUTE={n_abs} ESTIMATE∩ABSOLUTE={n_both}")

    results = []
    results.append(run_block(feat.drop(columns="cohort"), purity, "NSCLC_pooled"))
    for ck in ["LUAD", "LUSC"]:
        sub = feat[feat["cohort"] == ck].drop(columns="cohort")
        results.append(run_block(sub, purity, ck))
    res = pd.concat(results, ignore_index=True)

    # Same-n complete case: only samples with BOTH ESTIMATE and ABSOLUTE.
    # This is the fair head-to-head the user asked for (ESTIMATE vs ABSOLUTE).
    both_ids = common_ids[both_mask]
    feat_both = feat.loc[both_ids].drop(columns="cohort")
    cc_nsclc = run_block(feat_both, purity, "NSCLC_ESTIMATE_ABSOLUTE_complete",
                         purity_cols=["ESTIMATE", "ABSOLUTE"])
    cc_parts = [cc_nsclc]
    for ck in ["LUAD", "LUSC"]:
        ids = feat.index[(feat["cohort"] == ck) & feat.index.isin(both_ids)]
        cc_parts.append(run_block(feat.loc[ids].drop(columns="cohort"),
                                  purity, f"{ck}_ESTIMATE_ABSOLUTE_complete",
                                  purity_cols=["ESTIMATE", "ABSOLUTE"]))
    cc = pd.concat(cc_parts, ignore_index=True)

    vs_purity = tacstd2_vs_purity(feat.drop(columns="cohort"), purity, cohort)

    fmt_csv(res).to_csv(os.path.join(OUT, "tcga_tacstd2_correlations.csv"),
                        index=False)
    fmt_csv(cc).to_csv(os.path.join(OUT, "tcga_estimate_vs_absolute_completen.csv"),
                       index=False)
    fmt_csv(vs_purity).to_csv(os.path.join(OUT, "tcga_tacstd2_vs_purity.csv"),
                              index=False)
    write_direction_tally(res, os.path.join(OUT, "tcga_direction_tally.csv"))
    write_direction_tally(cc, os.path.join(OUT, "tcga_direction_tally_completen.csv"))

    # Compact signature-only table for the writeup (pooled + histologies).
    sig = res[res.feature.isin(SIG_FEATURES)].copy()
    keep = ["cohort", "feature", "n",
            "spearman_rho", "spearman_p", "spearman_q",
            "partial_rho_ESTIMATE", "partial_p_ESTIMATE", "partial_q_ESTIMATE",
            "partial_n_ESTIMATE",
            "partial_rho_ABSOLUTE", "partial_p_ABSOLUTE", "partial_q_ABSOLUTE",
            "partial_n_ABSOLUTE",
            "partial_rho_CPE", "partial_p_CPE", "partial_q_CPE",
            "partial_n_CPE"]
    fmt_csv(sig[keep]).to_csv(os.path.join(OUT, "tcga_signatures_estimate_absolute.csv"),
                              index=False)

    with open(os.path.join(OUT, "tcga_signature_genes.txt"), "w") as fh:
        for name, genes in C.SIGNATURES.items():
            fh.write(f"{name}: {', '.join(genes)}\n")
        fh.write(f"TARGET: {C.TARGET_GENE}\n")
        fh.write(f"TJ_GENES: {', '.join(C.TJ_GENES)}\n")

    with open(os.path.join(OUT, "tcga_sample_counts.txt"), "w") as fh:
        fh.write(f"pooled_primary_tumors\t{pooled.shape[1]}\n")
        fh.write(f"LUAD_primary\t{int((cohort=='LUAD').sum())}\n")
        fh.write(f"LUSC_primary\t{int((cohort=='LUSC').sum())}\n")
        fh.write(f"with_purity_row\t{len(common_ids)}\n")
        fh.write(f"nonmissing_CPE\t{n_cpe}\n")
        fh.write(f"nonmissing_ESTIMATE\t{n_est}\n")
        fh.write(f"nonmissing_ABSOLUTE\t{n_abs}\n")
        fh.write(f"ESTIMATE_and_ABSOLUTE\t{n_both}\n")
        fh.write(f"LUAD_ESTIMATE_and_ABSOLUTE\t"
                 f"{int(((cohort=='LUAD') & feat.index.isin(both_ids)).sum())}\n")
        fh.write(f"LUSC_ESTIMATE_and_ABSOLUTE\t"
                 f"{int(((cohort=='LUSC') & feat.index.isin(both_ids)).sum())}\n")

    print("\n=== TACSTD2 vs purity ===")
    print(fmt_csv(vs_purity).to_string(index=False))
    print("\n=== POOLED signatures: Spearman / ESTIMATE / ABSOLUTE / CPE ===")
    show = ["feature", "spearman_rho", "partial_rho_ESTIMATE",
            "partial_p_ESTIMATE", "partial_rho_ABSOLUTE",
            "partial_p_ABSOLUTE", "partial_rho_CPE"]
    print(fmt_csv(sig[sig.cohort == "NSCLC_pooled"][show]).to_string(index=False))
    print("\n=== Complete-case (ESTIMATE∩ABSOLUTE) pooled signatures ===")
    show2 = ["feature", "n", "spearman_rho",
             "partial_rho_ESTIMATE", "partial_p_ESTIMATE",
             "partial_rho_ABSOLUTE", "partial_p_ABSOLUTE"]
    print(fmt_csv(cc_nsclc[cc_nsclc.feature.isin(SIG_FEATURES)][show2]
                  ).to_string(index=False))
    print("\n=== Direction tally (immune features) ===")
    print(pd.read_csv(os.path.join(OUT, "tcga_direction_tally.csv")).to_string(index=False))


if __name__ == "__main__":
    main()
