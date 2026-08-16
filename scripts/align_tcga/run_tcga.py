"""
USER-ALIGN: replicate the user's bulk finding in TCGA LUAD+LUSC (NSCLC).

Method (matching the user's description):
  1. Spearman correlation of TACSTD2 (TROP-2) vs immune / cytotoxic /
     exhaustion signatures.
  2. Partial Spearman correlation after controlling for tumour purity
     (ESTIMATE, ABSOLUTE and CPE from Aran et al. 2015).

Data (public, open access):
  * Expression: UCSC Xena TCGA LUAD & LUSC HiSeqV2 (RSEM, log2(norm+1)).
  * Purity: Aran et al. Nat Commun 2015, Supplementary Data 1 (CPE table).

Outputs go to results/align_tcga/.
"""
from __future__ import annotations

import os
import numpy as np
import pandas as pd

import common as C

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "..", "results", "align_tcga", "data")
OUT = os.path.join(HERE, "..", "..", "results", "align_tcga")


def load_purity() -> pd.DataFrame:
    xl = pd.ExcelFile(os.path.join(DATA, "Aran_CPE_purity.xlsx"))
    raw = xl.parse(xl.sheet_names[0], header=None)
    # Locate the header row containing "Sample ID".
    hdr = raw.index[raw.iloc[:, 0].astype(str).str.strip() == "Sample ID"][0]
    names = [str(v).strip() for v in raw.iloc[hdr].tolist()]
    df = raw.iloc[hdr + 1:].copy()
    df.columns = names
    df = df.rename(columns={"Sample ID": "SampleID",
                            "Cancer type": "CancerType"})
    df = df[df["SampleID"].astype(str).str.startswith("TCGA")].copy()
    # Barcode in purity table has a vial letter (TCGA-XX-XXXX-01A); trim to 15.
    df["barcode15"] = df["SampleID"].astype(str).str[:15]
    for col in ["ESTIMATE", "ABSOLUTE", "CPE"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.drop_duplicates("barcode15")
    return df.set_index("barcode15")[["CancerType", "ESTIMATE",
                                      "ABSOLUTE", "CPE"]]


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


def compute_features(expr: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Return a samples x features frame: TACSTD2, signatures, single genes."""
    feats = {}
    used_map = {}
    # target gene
    feats[C.TARGET_GENE] = expr.loc[C.TARGET_GENE].astype(float)
    # signatures
    for name, genes in C.SIGNATURES.items():
        score, used = C.signature_score(expr, genes)
        feats[name] = score
        used_map[name] = used
    # single immune genes
    for g in C.SINGLE_IMMUNE_GENES:
        r = C.resolve_symbol(g, set(expr.index))
        if r is not None:
            feats["gene_" + g] = expr.loc[r].astype(float)
    # tight-junction genes
    for g in C.TJ_GENES:
        r = C.resolve_symbol(g, set(expr.index))
        if r is not None:
            feats["TJ_" + g] = expr.loc[r].astype(float)
    return pd.DataFrame(feats), used_map


def run_block(feat: pd.DataFrame, purity: pd.DataFrame, tag: str) -> pd.DataFrame:
    """Spearman + partial-Spearman of TACSTD2 vs every other feature."""
    df = feat.join(purity[["ESTIMATE", "ABSOLUTE", "CPE"]], how="left")
    x = df[C.TARGET_GENE].to_numpy(float)
    rows = []
    targets = [c for c in feat.columns if c != C.TARGET_GENE]
    for col in targets:
        y = df[col].to_numpy(float)
        r, p, n = C.spearman(x, y)
        rec = {"cohort": tag, "feature": col,
               "spearman_rho": r, "spearman_p": p, "n": n}
        for pm in ["ESTIMATE", "ABSOLUTE", "CPE"]:
            z = df[pm].to_numpy(float)
            pr, pp, pn = C.partial_spearman(x, y, z)
            rec[f"partial_rho_{pm}"] = pr
            rec[f"partial_p_{pm}"] = pp
            rec[f"partial_n_{pm}"] = pn
        rows.append(rec)
    return pd.DataFrame(rows)


def main() -> None:
    pooled, cohort = build_expression()
    purity = load_purity()

    print(f"[TCGA] pooled primary tumours: {pooled.shape[1]} "
          f"(LUAD={int((cohort=='LUAD').sum())}, "
          f"LUSC={int((cohort=='LUSC').sum())}); genes={pooled.shape[0]}")

    feat, used_map = compute_features(pooled)
    feat["cohort"] = cohort

    # overlap with purity
    common_ids = feat.index.intersection(purity.index)
    print(f"[TCGA] samples with purity (CPE table): {len(common_ids)}")
    n_cpe = purity.loc[common_ids, "CPE"].notna().sum()
    n_est = purity.loc[common_ids, "ESTIMATE"].notna().sum()
    n_abs = purity.loc[common_ids, "ABSOLUTE"].notna().sum()
    print(f"[TCGA] non-missing purity -> CPE={n_cpe}, ESTIMATE={n_est}, "
          f"ABSOLUTE={n_abs}")

    results = []
    # pooled NSCLC
    results.append(run_block(feat.drop(columns="cohort"), purity, "NSCLC_pooled"))
    # per cohort robustness
    for ck in ["LUAD", "LUSC"]:
        sub = feat[feat["cohort"] == ck].drop(columns="cohort")
        results.append(run_block(sub, purity, ck))

    res = pd.concat(results, ignore_index=True)
    res = res.round(4)
    out_csv = os.path.join(OUT, "tcga_tacstd2_correlations.csv")
    res.to_csv(out_csv, index=False)
    print(f"[TCGA] wrote {out_csv} ({len(res)} rows)")

    # signature membership record (transparency)
    with open(os.path.join(OUT, "tcga_signature_genes.txt"), "w") as fh:
        for name, genes in C.SIGNATURES.items():
            fh.write(f"{name}: {', '.join(genes)}\n")
        fh.write(f"TARGET: {C.TARGET_GENE}\n")
        fh.write(f"TJ_GENES: {', '.join(C.TJ_GENES)}\n")

    # cohort n record
    with open(os.path.join(OUT, "tcga_sample_counts.txt"), "w") as fh:
        fh.write(f"pooled_primary_tumors\t{pooled.shape[1]}\n")
        fh.write(f"LUAD_primary\t{int((cohort=='LUAD').sum())}\n")
        fh.write(f"LUSC_primary\t{int((cohort=='LUSC').sum())}\n")
        fh.write(f"with_purity_row\t{len(common_ids)}\n")
        fh.write(f"nonmissing_CPE\t{int(n_cpe)}\n")
        fh.write(f"nonmissing_ESTIMATE\t{int(n_est)}\n")
        fh.write(f"nonmissing_ABSOLUTE\t{int(n_abs)}\n")

    # concise console summary of the pooled result
    pooled_res = res[res["cohort"] == "NSCLC_pooled"]
    print("\n=== POOLED NSCLC: TACSTD2 vs signatures ===")
    show = ["feature", "spearman_rho", "spearman_p",
            "partial_rho_CPE", "partial_p_CPE", "partial_n_CPE"]
    print(pooled_res[show].to_string(index=False))


if __name__ == "__main__":
    main()
