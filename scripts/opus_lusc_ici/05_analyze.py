#!/usr/bin/env python3
"""Step 5 - TACSTD2 / CLDN4 versus ICI response and immune context, split by histology.

Primary analysis set
    pathology-confirmed LUSC, pre-treatment, with a binary benefit label.

Sensitivity set
    primary plus high-confidence inferred LUSC (GSE135222 / GSE126044).

Contrast set
    pathology-confirmed non-LUSC (mostly LUAD), same endpoints, reported so the
    histology split is visible rather than mixed.

Nothing is fabricated: every row written below is computed from the public
matrices in ``results/opus_lusc_ici/data``.  Cohorts that lack a gene or an
endpoint are omitted from that comparison and counted as missing.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from common import (benjamini_hochberg, hedges_g, load_signatures, mann_whitney,
                    partial_spearman, permutation_pvalue, random_effects_meta,
                    signature_coverage, signature_score, spearman, write_table)
from config import (COMPANION_GENES, DATA_DIR, LOGS_DIR, RANDOM_SEED,
                    TABLES_DIR, TARGET_GENES)

N_PERM = 10000
HOUSEKEEPING = ["ACTB", "GAPDH", "PPIA", "RPLP0", "HPRT1", "TBP"]
NEGCTRL_IMMUNE = ["ACTB", "GAPDH", "PPIA", "RPLP0"]  # should not track CYT / IFNG
IMMUNE_FOCUS = ["CYT", "IFNG_6", "TIS_18", "CD8_T_CELL", "CHECKPOINT",
                "ANTIGEN_PRESENTATION", "TGFB_FTBRS", "EPITHELIAL", "PROLIFERATION"]


def load_expr(cohort: str) -> pd.DataFrame:
    return pd.read_csv(DATA_DIR / f"{cohort}_expr.tsv.gz", sep="\t", index_col=0)


def build_sample_table(clinical: pd.DataFrame, sigs: dict) -> pd.DataFrame:
    rows = []
    coverage_rows = []
    genes_wanted = TARGET_GENES + COMPANION_GENES + HOUSEKEEPING
    for cohort, sub in clinical.groupby("cohort", sort=False):
        expr = load_expr(cohort)
        for gene in genes_wanted:
            if gene in expr.index:
                sub = sub.copy()
                sub[gene] = expr.loc[gene, sub["sample_id"]].to_numpy()
            else:
                sub[gene] = np.nan
        for name, spec in sigs.items():
            present, total = signature_coverage(expr, spec["genes"])
            coverage_rows.append(dict(cohort=cohort, signature=name,
                                      genes_present=present, genes_total=total,
                                      fraction=present / total if total else np.nan))
            if present >= 2:
                scores = signature_score(expr, spec["genes"])
                sub[f"sig_{name}"] = scores.reindex(sub["sample_id"]).to_numpy()
            else:
                sub[f"sig_{name}"] = np.nan
        rows.append(sub)
    write_table(pd.DataFrame(coverage_rows), TABLES_DIR / "signature_coverage.csv")
    return pd.concat(rows, ignore_index=True)


def fisher_z_meta(rhos, ns):
    """Random-effects meta-analysis of Spearman rho via Fisher z."""
    z, var = [], []
    for rho, n in zip(rhos, ns):
        if pd.isna(rho) or n is None or n < 6 or abs(rho) >= 1:
            continue
        zi = np.arctanh(np.clip(rho, -0.999999, 0.999999))
        z.append(zi)
        var.append(1 / (n - 3))
    pooled = random_effects_meta(z, var)
    if pd.isna(pooled["estimate"]):
        return pooled | {"rho": np.nan, "rho_ci_low": np.nan, "rho_ci_high": np.nan}
    rho = float(np.tanh(pooled["estimate"]))
    return pooled | {
        "rho": rho,
        "rho_ci_low": float(np.tanh(pooled["ci_low"])),
        "rho_ci_high": float(np.tanh(pooled["ci_high"])),
    }


def response_rows(frame: pd.DataFrame, gene: str, stratum: str, analysis_set: str) -> list[dict]:
    out = []
    for cohort, sub in frame.groupby("cohort", sort=False):
        x = sub.loc[sub["benefit"] == 1, gene]
        y = sub.loc[sub["benefit"] == 0, gene]
        if x.notna().sum() < 2 or y.notna().sum() < 2:
            continue
        mw = mann_whitney(x, y)
        g, se = hedges_g(x, y)
        seed = (RANDOM_SEED + hash((cohort, gene, stratum, analysis_set)) % 100000) % (2**31)
        p_perm = permutation_pvalue(x, y, lambda a, b: hedges_g(a, b)[0],
                                    n_perm=N_PERM, seed=seed)
        out.append(dict(
            analysis_set=analysis_set, stratum=stratum, cohort=cohort, gene=gene,
            n_benefit=int(x.notna().sum()), n_nobenefit=int(y.notna().sum()),
            median_benefit=float(x.median()), median_nobenefit=float(y.median()),
            mean_benefit=float(x.mean()), mean_nobenefit=float(y.mean()),
            hedges_g=g, hedges_g_se=se,
            mw_u=mw["u"], mw_p=mw["p"], auc=mw["auc"], rank_biserial=mw["rbc"],
            permutation_p=p_perm,
        ))
    return out


def meta_from_response(rows: list[dict], analysis_set: str, stratum: str, gene: str) -> dict:
    use = [r for r in rows if r["analysis_set"] == analysis_set
           and r["stratum"] == stratum and r["gene"] == gene]
    if not use:
        return dict(analysis_set=analysis_set, stratum=stratum, gene=gene, k=0)
    pooled = random_effects_meta([r["hedges_g"] for r in use],
                                 [r["hedges_g_se"] ** 2 for r in use])
    return dict(analysis_set=analysis_set, stratum=stratum, gene=gene,
                n_total=int(sum(r["n_benefit"] + r["n_nobenefit"] for r in use)),
                n_benefit=int(sum(r["n_benefit"] for r in use)),
                n_nobenefit=int(sum(r["n_nobenefit"] for r in use)),
                **pooled)


def cox_pfs(sub: pd.DataFrame, gene: str) -> dict:
    try:
        from lifelines import CoxPHFitter
    except ImportError:
        return dict(n=len(sub), hr=np.nan, hr_lo=np.nan, hr_hi=np.nan, p=np.nan, note="lifelines missing")
    frame = sub[["pfs_time_months", "pfs_event", gene]].dropna()
    frame = frame[frame["pfs_time_months"] > 0]
    if len(frame) < 8 or frame["pfs_event"].sum() < 3:
        return dict(n=len(frame), events=int(frame["pfs_event"].sum()) if len(frame) else 0,
                    hr=np.nan, hr_lo=np.nan, hr_hi=np.nan, p=np.nan, note="too few events")
    # Standardise the gene so the HR is per +1 SD within the analysed subset.
    frame = frame.copy()
    sd = frame[gene].std(ddof=1)
    if not sd or sd == 0:
        return dict(n=len(frame), events=int(frame["pfs_event"].sum()),
                    hr=np.nan, hr_lo=np.nan, hr_hi=np.nan, p=np.nan, note="zero variance")
    frame[gene] = (frame[gene] - frame[gene].mean()) / sd
    cph = CoxPHFitter()
    try:
        cph.fit(frame, duration_col="pfs_time_months", event_col="pfs_event")
    except Exception as exc:  # noqa: BLE001 - report and move on
        return dict(n=len(frame), events=int(frame["pfs_event"].sum()),
                    hr=np.nan, hr_lo=np.nan, hr_hi=np.nan, p=np.nan, note=str(exc)[:120])
    s = cph.summary.loc[gene]
    return dict(n=int(len(frame)), events=int(frame["pfs_event"].sum()),
                hr=float(s["exp(coef)"]), hr_lo=float(s["exp(coef) lower 95%"]),
                hr_hi=float(s["exp(coef) upper 95%"]), p=float(s["p"]), note="ok")


def immune_rows(frame: pd.DataFrame, gene: str, signatures: list[str],
                stratum: str, analysis_set: str) -> list[dict]:
    out = []
    for cohort, sub in frame.groupby("cohort", sort=False):
        for sig in signatures:
            col = f"sig_{sig}"
            if col not in sub.columns:
                continue
            sp = spearman(sub[gene], sub[col])
            if sp["n"] < 6:
                continue
            cov = {}
            if "sig_EPITHELIAL" in sub.columns and sig != "EPITHELIAL":
                cov = partial_spearman(sub[gene], sub[col],
                                       sub[["sig_EPITHELIAL"]].rename(columns={"sig_EPITHELIAL": "epithelial"}))
            out.append(dict(
                analysis_set=analysis_set, stratum=stratum, cohort=cohort,
                gene=gene, signature=sig, n=sp["n"], rho=sp["rho"], p=sp["p"],
                partial_rho_epithelial=cov.get("rho", np.nan),
                partial_p_epithelial=cov.get("p", np.nan),
                partial_n=cov.get("n", np.nan),
            ))
    return out


def meta_immune(rows: list[dict], analysis_set: str, stratum: str,
                gene: str, signature: str) -> dict:
    use = [r for r in rows if r["analysis_set"] == analysis_set and r["stratum"] == stratum
           and r["gene"] == gene and r["signature"] == signature]
    if not use:
        return dict(analysis_set=analysis_set, stratum=stratum, gene=gene,
                    signature=signature, k=0)
    pooled = fisher_z_meta([r["rho"] for r in use], [r["n"] for r in use])
    return dict(analysis_set=analysis_set, stratum=stratum, gene=gene, signature=signature,
                n_total=int(sum(r["n"] for r in use)), **{k: pooled[k] for k in
                ("k", "estimate", "se", "ci_low", "ci_high", "p", "tau2", "q", "i2",
                 "rho", "rho_ci_low", "rho_ci_high")})


def leave_one_out(rows: list[dict], analysis_set: str, stratum: str, gene: str) -> list[dict]:
    use = [r for r in rows if r["analysis_set"] == analysis_set
           and r["stratum"] == stratum and r["gene"] == gene]
    out = []
    for drop in {r["cohort"] for r in use}:
        keep = [r for r in use if r["cohort"] != drop]
        pooled = random_effects_meta([r["hedges_g"] for r in keep],
                                     [r["hedges_g_se"] ** 2 for r in keep])
        out.append(dict(analysis_set=analysis_set, stratum=stratum, gene=gene,
                        dropped=drop, **pooled))
    return out


def tcga_immune(sigs: dict) -> pd.DataFrame:
    expr = pd.read_csv(DATA_DIR / "TCGA_LUSC_expr.tsv.gz", sep="\t", index_col=0)
    rows = []
    for gene in TARGET_GENES:
        if gene not in expr.index:
            continue
        for name in IMMUNE_FOCUS:
            present, total = signature_coverage(expr, sigs[name]["genes"])
            if present < 2:
                continue
            score = signature_score(expr, sigs[name]["genes"])
            sp = spearman(expr.loc[gene], score)
            cov = {}
            if name != "EPITHELIAL":
                epi = signature_score(expr, sigs["EPITHELIAL"]["genes"])
                cov = partial_spearman(expr.loc[gene], score,
                                       pd.DataFrame({"epithelial": epi}))
            rows.append(dict(
                dataset="TCGA-LUSC", gene=gene, signature=name,
                genes_present=present, genes_total=total,
                n=sp["n"], rho=sp["rho"], p=sp["p"],
                partial_rho_epithelial=cov.get("rho", np.nan),
                partial_p_epithelial=cov.get("p", np.nan),
            ))
    frame = pd.DataFrame(rows)
    if len(frame):
        frame["q"] = benjamini_hochberg(frame["p"])
    return frame


def main() -> int:
    clinical = pd.read_csv(TABLES_DIR / "cohort_clinical_with_histology.csv")
    clinical["histology"] = clinical["histology"].fillna("NA")
    clinical["histology_final"] = clinical["histology_final"].fillna("NA")
    sigs = load_signatures()

    samples = build_sample_table(clinical, sigs)
    write_table(samples, TABLES_DIR / "sample_level_scores.csv")

    # Analysis sets.  Pre-treatment only: post-treatment samples would leak
    # on-treatment immune changes into a "predictor" analysis.
    pre = samples[samples["timepoint"] == "pre-treatment"].copy()
    sets = {
        "primary_pathology_LUSC": pre[(pre["histology_tier"] == "pathology")
                                      & (pre["histology_final"] == "LUSC")],
        "sensitivity_LUSC_plus_inferred": pre[pre["histology_final"] == "LUSC"],
        "contrast_pathology_nonLUSC": pre[(pre["histology_tier"] == "pathology")
                                          & (pre["histology_final"] == "non-LUSC")],
    }
    SET_STRATUM = {
        "primary_pathology_LUSC": "LUSC",
        "sensitivity_LUSC_plus_inferred": "LUSC",
        "contrast_pathology_nonLUSC": "non-LUSC",
    }

    inventory = []
    for name, frame in sets.items():
        for cohort, sub in frame.groupby("cohort"):
            inventory.append(dict(
                analysis_set=name, cohort=cohort, n=len(sub),
                n_with_benefit=int(sub["benefit"].notna().sum()),
                n_benefit=int((sub["benefit"] == 1).sum()),
                n_nobenefit=int((sub["benefit"] == 0).sum()),
                n_with_pfs=int(sub["pfs_time_months"].notna().sum()),
                n_TACSTD2=int(sub["TACSTD2"].notna().sum()),
                n_CLDN4=int(sub["CLDN4"].notna().sum()),
            ))
    write_table(pd.DataFrame(inventory), TABLES_DIR / "analysis_set_inventory.csv")
    print("\nAnalysis-set inventory")
    print(pd.DataFrame(inventory).to_string(index=False))

    # ----- response -----
    resp_rows = []
    for set_name, frame in sets.items():
        usable = frame[frame["benefit"].notna()]
        stratum = SET_STRATUM[set_name]
        for gene in TARGET_GENES + HOUSEKEEPING[:3]:
            resp_rows.extend(response_rows(usable, gene, stratum, set_name))
    resp = pd.DataFrame(resp_rows)
    write_table(resp, TABLES_DIR / "response_by_cohort.csv")

    meta_rows = []
    for set_name in sets:
        stratum = SET_STRATUM[set_name]
        for gene in TARGET_GENES + HOUSEKEEPING[:3]:
            meta_rows.append(meta_from_response(resp_rows, set_name, stratum, gene))
    meta = pd.DataFrame(meta_rows)
    write_table(meta, TABLES_DIR / "response_meta.csv")

    loo_rows = []
    for set_name in ("primary_pathology_LUSC", "sensitivity_LUSC_plus_inferred"):
        for gene in TARGET_GENES:
            loo_rows.extend(leave_one_out(resp_rows, set_name, "LUSC", gene))
    write_table(pd.DataFrame(loo_rows), TABLES_DIR / "response_leave_one_cohort_out.csv")

    # ----- PFS Cox -----
    cox_rows = []
    for set_name, frame in sets.items():
        stratum = SET_STRATUM[set_name]
        for gene in TARGET_GENES:
            for cohort, sub in frame.groupby("cohort"):
                if sub["pfs_time_months"].notna().sum() < 8:
                    continue
                fit = cox_pfs(sub, gene)
                cox_rows.append(dict(analysis_set=set_name, stratum=stratum,
                                     cohort=cohort, gene=gene, **fit))
            # pooled-within-set Cox on within-cohort z-scores (avoids mixing platforms)
            pieces = []
            for cohort, sub in frame.groupby("cohort"):
                z = sub[["pfs_time_months", "pfs_event", gene]].dropna()
                if len(z) < 4 or z[gene].std(ddof=1) == 0:
                    continue
                z = z.copy()
                z[gene] = (z[gene] - z[gene].mean()) / z[gene].std(ddof=1)
                z["cohort"] = cohort
                pieces.append(z)
            if pieces:
                pooled = pd.concat(pieces)
                fit = cox_pfs(pooled, gene)
                cox_rows.append(dict(analysis_set=set_name, stratum=stratum,
                                     cohort="POOLED_within_cohort_z", gene=gene, **fit))
    write_table(pd.DataFrame(cox_rows), TABLES_DIR / "pfs_cox.csv")

    # ----- immune -----
    imm_rows = []
    for set_name, frame in sets.items():
        stratum = SET_STRATUM[set_name]
        for gene in TARGET_GENES + NEGCTRL_IMMUNE:
            imm_rows.extend(immune_rows(frame, gene, IMMUNE_FOCUS, stratum, set_name))
    imm = pd.DataFrame(imm_rows)
    if len(imm):
        imm["q"] = np.nan
        for (aset, gene), idx in imm.groupby(["analysis_set", "gene"]).groups.items():
            imm.loc[idx, "q"] = benjamini_hochberg(imm.loc[idx, "p"])
    write_table(imm, TABLES_DIR / "immune_by_cohort.csv")

    imm_meta = []
    for set_name in sets:
        stratum = SET_STRATUM[set_name]
        for gene in TARGET_GENES + NEGCTRL_IMMUNE:
            for sig in IMMUNE_FOCUS:
                imm_meta.append(meta_immune(imm_rows, set_name, stratum, gene, sig))
    imm_meta_df = pd.DataFrame(imm_meta)
    if len(imm_meta_df):
        for aset, idx in imm_meta_df.groupby("analysis_set").groups.items():
            imm_meta_df.loc[idx, "q"] = benjamini_hochberg(imm_meta_df.loc[idx, "p"])
    write_table(imm_meta_df, TABLES_DIR / "immune_meta.csv")

    tcga = tcga_immune(sigs)
    write_table(tcga, TABLES_DIR / "tcga_lusc_immune.csv")

    # ----- gene-gene -----
    gg_rows = []
    for set_name, frame in sets.items():
        stratum = SET_STRATUM[set_name]
        for cohort, sub in frame.groupby("cohort"):
            sp = spearman(sub["TACSTD2"], sub["CLDN4"])
            if sp["n"] >= 6:
                gg_rows.append(dict(analysis_set=set_name, stratum=stratum,
                                    cohort=cohort, **sp))
    write_table(pd.DataFrame(gg_rows), TABLES_DIR / "tacstd2_cldn4_correlation.csv")

    # ----- verification summary (computed, not asserted) -----
    primary_meta = meta[(meta["analysis_set"] == "primary_pathology_LUSC")
                        & (meta["gene"].isin(TARGET_GENES))]
    hk_meta = meta[(meta["analysis_set"] == "primary_pathology_LUSC")
                   & (meta["gene"].isin(HOUSEKEEPING[:3]))]
    verification = {
        "n_perm": N_PERM,
        "primary_response_meta": primary_meta.to_dict(orient="records"),
        "housekeeping_response_meta": hk_meta.to_dict(orient="records"),
        "notes": [
            "Permutation p-values randomise benefit labels within each cohort.",
            "Housekeeping genes (ACTB, GAPDH, PPIA) are negative controls for the response contrast.",
            "Leave-one-cohort-out re-estimates the random-effects meta after dropping each cohort.",
            "Partial Spearman removes a linear epithelial-content score before correlating with immune signatures.",
            "TCGA-LUSC is treatment-naive surgical resection and is used only to replicate the immune axis, not ICI response.",
        ],
    }
    (TABLES_DIR / "verification_bundle.json").write_text(
        json.dumps(verification, indent=2, default=float) + "\n")

    print("\nPrimary pathology-LUSC response meta")
    cols = [c for c in ("gene", "k", "n_total", "n_benefit", "n_nobenefit",
                        "estimate", "ci_low", "ci_high", "p", "i2") if c in primary_meta.columns]
    if len(primary_meta):
        print(primary_meta[cols].to_string(index=False))
    else:
        print("  (no rows — check inventory)")

    (LOGS_DIR / "05_analyze.done").write_text("ok\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
