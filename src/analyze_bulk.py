"""A11 bulk analysis: TACSTD2 vs the CD47 / Galectin / Nectin / TGFB axes in
TCGA-LUAD and TCGA-LUSC, with and without adjustment for tumour purity.

Implements analyses A1-A4 of the pre-registered plan.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C
from stats_utils import (spearman, partial_spearman, cliffs_delta, bh_fdr,
                         zscore_module)

DERIVED = C.DATA / "derived"
TAB = C.OUT / "tables"
TAB.mkdir(parents=True, exist_ok=True)


def load(cohort):
    mat = pd.read_parquet(DERIVED / f"bulk_{cohort}_panel.parquet")
    samp = pd.read_parquet(DERIVED / f"bulk_{cohort}_samples.parquet")
    return mat, samp


def add_scores(mat):
    """Compartment scores, axis ligand module scores, and the TGFB activity score."""
    sc = pd.DataFrame(index=mat.index)
    used = {}
    for name, genes in C.COMPARTMENT.items():
        sc[f"score_{name}"], used[name] = zscore_module(mat, genes)
    for axis, genes in C.AXIS_LIGANDS.items():
        sc[f"module_{axis}"], used[f"module_{axis}"] = zscore_module(mat, genes)
    sc["module_TGFB_activity_FTBRS"], used["FTBRS"] = zscore_module(mat, C.F_TBRS)
    return sc, used


# ---------------------------------------------------------------------------
# A1 + A2: gene-level naive and purity-adjusted association
# ---------------------------------------------------------------------------
def gene_level(mat, samp, cohort):
    tac = mat[C.TARGET]
    purity = samp["purity"]
    axis_of = C.gene_to_axis()

    panel = C.primary_panel()
    controls = [g for g in C.NEG_CONTROL if g in mat.columns]

    # Pre-specified tertile split on TACSTD2.
    lo_cut, hi_cut = tac.quantile([1 / 3, 2 / 3])
    grp = pd.Series(np.where(tac >= hi_cut, "high",
                    np.where(tac <= lo_cut, "low", "mid")), index=tac.index)

    rows = []
    for gene in panel + controls:
        if gene not in mat.columns or gene == C.TARGET:
            continue
        y = mat[gene]
        naive = spearman(tac, y)
        adj = partial_spearman(tac, y, [purity])
        d = cliffs_delta(y[grp == "high"], y[grp == "low"])

        # Sensitivity: adjust for an RNA-based epithelial score instead of
        # DNA-based purity. Over-adjustment is possible here because TACSTD2 is
        # itself an epithelial gene, so this is reported as a bound, not a
        # preferred estimate.
        adj_epi = partial_spearman(tac, y, [samp["score_epithelial"]])

        rows.append(dict(
            cohort=cohort, gene=gene,
            axis=axis_of.get(gene, "negative_control"),
            is_control=gene in controls,
            mean_log2tpm=float(y.mean()),
            rho_naive=naive["rho"], rho_naive_lo=naive["lo"],
            rho_naive_hi=naive["hi"], p_naive=naive["p"], n=naive["n"],
            rho_adj_purity=adj["rho"], rho_adj_lo=adj["lo"],
            rho_adj_hi=adj["hi"], p_adj=adj["p"], n_adj=adj["n"],
            rho_adj_epithelial=adj_epi["rho"], p_adj_epithelial=adj_epi["p"],
            cliffs_delta_hi_vs_lo=d["delta"], p_tertile=d["p"],
            n_high=d["n1"], n_low=d["n2"],
        ))

    df = pd.DataFrame(rows)
    # FDR is computed separately within the primary panel and within the
    # negative controls, so controls neither dilute nor inflate the correction.
    for mask, tag in [(~df.is_control, "panel"), (df.is_control, "control")]:
        for pcol, qcol in [("p_naive", "q_naive"), ("p_adj", "q_adj"),
                           ("p_tertile", "q_tertile")]:
            df.loc[mask, qcol] = bh_fdr(df.loc[mask, pcol].values)

    df["attenuation"] = 1 - (df["rho_adj_purity"].abs() / df["rho_naive"].abs())
    df["sign_flip"] = np.sign(df["rho_adj_purity"]) != np.sign(df["rho_naive"])
    return df.sort_values(["axis", "rho_naive"], ascending=[True, False])


# ---------------------------------------------------------------------------
# A2 (module level) + confounder diagnostics
# ---------------------------------------------------------------------------
def module_level(mat, samp, cohort):
    tac = mat[C.TARGET]
    purity = samp["purity"]
    cols = [c for c in samp.columns if c.startswith(("module_", "score_"))]

    rows = []
    for col in cols + ["purity"]:
        y = samp[col]
        naive = spearman(tac, y)
        adj = (dict(rho=np.nan, p=np.nan, n=np.nan, lo=np.nan, hi=np.nan)
               if col == "purity" else partial_spearman(tac, y, [purity]))
        rows.append(dict(cohort=cohort, feature=col,
                         rho_naive=naive["rho"], lo=naive["lo"], hi=naive["hi"],
                         p_naive=naive["p"], n=naive["n"],
                         rho_adj_purity=adj["rho"], p_adj=adj["p"]))
    df = pd.DataFrame(rows)
    df["q_naive"] = bh_fdr(df["p_naive"].values)
    df["q_adj"] = bh_fdr(df["p_adj"].values)
    return df


# ---------------------------------------------------------------------------
# A3: association within strata of tumour purity
# ---------------------------------------------------------------------------
def purity_strata(mat, samp, cohort):
    tac = mat[C.TARGET]
    pur = samp["purity"]
    ok = pur.notna()
    q1, q2 = pur[ok].quantile([1 / 3, 2 / 3])
    strat = pd.Series(np.where(pur >= q2, "high_purity",
                      np.where(pur <= q1, "low_purity", "mid_purity")),
                      index=pur.index).where(ok)

    feats = [c for c in samp.columns if c.startswith("module_")]
    rows = []
    for s in ["low_purity", "mid_purity", "high_purity"]:
        idx = strat[strat == s].index
        for f in feats:
            r = spearman(tac.loc[idx], samp.loc[idx, f])
            rows.append(dict(cohort=cohort, stratum=s, feature=f,
                             purity_range=f"{pur.loc[idx].min():.2f}-{pur.loc[idx].max():.2f}",
                             rho=r["rho"], lo=r["lo"], hi=r["hi"],
                             p=r["p"], n=r["n"]))
    df = pd.DataFrame(rows)
    df["q"] = bh_fdr(df["p"].values)
    return df


def tumor_vs_normal(mat_all, samp_all, cohort):
    """Context: is TACSTD2 even tumour-enriched in lung?"""
    rows = []
    for gene in [C.TARGET, "EPCAM", "CD47", "LGALS1", "LGALS3", "LGALS9",
                 "PVR", "NECTIN2", "NECTIN4", "TGFB1", "TGFB2", "TGFB3"]:
        if gene not in mat_all.columns:
            continue
        t = mat_all.loc[samp_all.sample_type == "PrimaryTumor", gene]
        n = mat_all.loc[samp_all.sample_type == "NormalAdjacent", gene]
        d = cliffs_delta(t, n)
        rows.append(dict(cohort=cohort, gene=gene,
                         mean_tumor=float(t.mean()), mean_normal=float(n.mean()),
                         log2fc_tumor_vs_normal=float(t.mean() - n.mean()),
                         cliffs_delta=d["delta"], p=d["p"],
                         n_tumor=d["n1"], n_normal=d["n2"]))
    df = pd.DataFrame(rows)
    df["q"] = bh_fdr(df["p"].values)
    return df


def main():
    gene_all, mod_all, strat_all, tn_all = [], [], [], []
    used_report = {}

    for cohort in ["LUAD", "LUSC"]:
        mat_all, samp_all = load(cohort)
        sc_all, used = add_scores(mat_all)
        samp_all = samp_all.join(sc_all)
        used_report[cohort] = used

        tn_all.append(tumor_vs_normal(mat_all, samp_all, cohort))

        keep = samp_all.sample_type == "PrimaryTumor"
        mat, samp = mat_all[keep], samp_all[keep]
        print(f"[{cohort}] primary tumours n={len(mat)}, "
              f"with purity n={samp['purity'].notna().sum()}", flush=True)

        gene_all.append(gene_level(mat, samp, cohort))
        mod_all.append(module_level(mat, samp, cohort))
        strat_all.append(purity_strata(mat, samp, cohort))

    pd.concat(gene_all).to_csv(TAB / "10_bulk_gene_level.csv", index=False)
    pd.concat(mod_all).to_csv(TAB / "11_bulk_module_level.csv", index=False)
    pd.concat(strat_all).to_csv(TAB / "12_bulk_purity_strata.csv", index=False)
    pd.concat(tn_all).to_csv(TAB / "13_bulk_tumor_vs_normal.csv", index=False)
    print("wrote bulk tables to", TAB, flush=True)


if __name__ == "__main__":
    main()
