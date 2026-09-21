#!/usr/bin/env python3
"""Human CLDN4 protein vs DNA-PK subunits and STING-pathway proteins.

Public source: Mirhadi et al., Nat Commun 2022 (PMID 35383172),
Supplementary Data 1, sheet "Normalized log2 Prot Quant".
TMT proteome of NSCLC PDX tumors (PRIDE PXD016579).

Human rows only. LUAD and LUSC are tested separately. Pairwise-complete
Spearman. CLDN4 NA is not imputed. Technical replicate pairs are averaged
to one model before the test.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

PROT_SHEET = "Normalized log2 Prot Quant"
CLIN_SHEET = "Clinical & omic-subtype info"
XLSX_NAME = "41467_2022_29444_MOESM4_ESM.xlsx"

# Pre-specified human tumor proteins.
DNA_PK = ["PRKDC", "XRCC5", "XRCC6"]
STING = ["STING1", "CGAS", "TBK1", "IRF3"]
PRIMARY = DNA_PK + STING
# Secondary sensors / IFN nodes. Not part of the DNA-PK or STING-core FDR family.
EXTENDED = [
    "IFI16",
    "STAT1",
    "STAT2",
    "IRF9",
    "MAVS",
    "DDX58",
    "IFIH1",
    "IRF7",
    "NFKB1",
    "RELA",
    "JAK1",
    "TYK2",
]
ALL_GENES = ["CLDN4"] + PRIMARY + EXTENDED
# Mouse rows inspected so they are not silently used as the human protein.
MOUSE_CONTRAST = ["PRKDC", "XRCC5", "XRCC6", "STING1", "TBK1", "STAT1"]

REPLICATE_PAIRS = {
    "PHLC113": ["PHLC113-X1", "PHLC113-X2"],
    "PHLC116": ["PHLC116-X1", "PHLC116-X2"],
    "PHLC277": ["PHLC277-X1", "PHLC277-X2"],
}

FAMILY = {g: "DNA-PK" for g in DNA_PK}
FAMILY.update({g: "STING" for g in STING})
FAMILY.update({g: "extended" for g in EXTENDED})


def bh(pvals: list[float]) -> list[float]:
    arr = np.asarray(pvals, dtype=float)
    out = np.full(arr.shape, np.nan)
    ok = np.isfinite(arr)
    if ok.any():
        out[ok] = stats.false_discovery_control(arr[ok], method="bh")
    return [float(x) if np.isfinite(x) else np.nan for x in out]


def spearman(x: pd.Series, y: pd.Series) -> dict:
    d = pd.concat([x, y], axis=1).dropna()
    n = int(len(d))
    rec = {"n_pairwise": n, "rho": np.nan, "p": np.nan, "tested": False, "note": "not tested"}
    if n == 0:
        rec["note"] = "no pairwise-complete models"
        return rec
    if n < 4:
        rec["note"] = f"n={n} < 4; Spearman not computed"
        return rec
    if d.iloc[:, 0].nunique() < 2 or d.iloc[:, 1].nunique() < 2:
        rec["note"] = f"n={n} but no variance"
        return rec
    rho, p = stats.spearmanr(d.iloc[:, 0], d.iloc[:, 1])
    rec.update(
        {
            "rho": float(rho),
            "p": float(p),
            "tested": True,
            "note": "two-sided Spearman, pairwise-complete log2 protein",
        }
    )
    return rec


def partial_spearman(y: pd.Series, x: pd.Series, cov: pd.Series) -> dict:
    """Spearman partial via Pearson correlation of rank residuals."""
    d = pd.concat([y, x, cov], axis=1).dropna()
    d.columns = ["y", "x", "c"]
    n = int(len(d))
    rec = {"n_pairwise": n, "rho": np.nan, "p": np.nan, "tested": False, "note": "not tested"}
    if n < 8 or d["y"].nunique() < 2 or d["x"].nunique() < 2 or d["c"].nunique() < 2:
        rec["note"] = f"n={n}; partial Spearman not computed"
        return rec

    def rank_resid(a: np.ndarray, b: np.ndarray) -> np.ndarray:
        ra = stats.rankdata(a).astype(float)
        rb = stats.rankdata(b).astype(float)
        slope, intercept = np.polyfit(rb, ra, 1)
        return ra - (slope * rb + intercept)

    ry = rank_resid(d["y"].to_numpy(), d["c"].to_numpy())
    rx = rank_resid(d["x"].to_numpy(), d["c"].to_numpy())
    rho, p = stats.pearsonr(ry, rx)
    rec.update(
        {
            "rho": float(rho),
            "p": float(p),
            "tested": True,
            "note": "partial Spearman, covariate = median log2 human proteome",
        }
    )
    return rec


def glass_rbc(detected: pd.Series, missing: pd.Series) -> dict:
    a = detected.dropna()
    b = missing.dropna()
    rec = {
        "n_cldn4_quantified": int(len(a)),
        "n_cldn4_missing": int(len(b)),
        "median_quantified": float(a.median()) if len(a) else np.nan,
        "median_missing": float(b.median()) if len(b) else np.nan,
        "rank_biserial": np.nan,
        "p": np.nan,
        "tested": False,
        "note": "not tested",
    }
    if len(a) < 3 or len(b) < 3:
        rec["note"] = "too few models on one side of CLDN4 detection"
        return rec
    if a.nunique() < 1 or b.nunique() < 1:
        rec["note"] = "no variance"
        return rec
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    rbc = (2.0 * float(u)) / (len(a) * len(b)) - 1.0
    rec.update(
        {
            "rank_biserial": float(rbc),
            "p": float(p),
            "tested": True,
            "note": "Mann-Whitney; Glass rank-biserial > 0 means higher protein when CLDN4 was quantified. NA is not imputed.",
        }
    )
    return rec


def load_protein_rows(prot: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    meta_cols = ["Gene Name", "ENTREZ", "HGNC", "H/M", "Entry", "Protein names"]
    samples = [c for c in prot.columns if c not in meta_cols]
    prot = prot.copy()
    prot["Gene Name"] = prot["Gene Name"].astype(str).str.strip().str.upper()
    prot["H/M"] = prot["H/M"].astype(str).str.strip()

    def grab(genes: list[str], species: str) -> pd.DataFrame:
        rows = []
        presence = []
        for gene in genes:
            hit = prot[(prot["Gene Name"] == gene) & (prot["H/M"] == species)]
            presence.append(
                {
                    "gene": gene,
                    "species_filter": species,
                    "n_rows": int(len(hit)),
                    "entry": ",".join(hit["Entry"].astype(str)) if len(hit) else "",
                    "protein_names": str(hit["Protein names"].iloc[0])[:180] if len(hit) else "",
                    "n_columns": len(samples),
                    "n_observed_columns": int(pd.to_numeric(hit.iloc[0, prot.columns.get_loc(samples[0]) :], errors="coerce").notna().sum())
                    if len(hit) == 1
                    else 0,
                }
            )
            if len(hit) != 1:
                continue
            s = pd.to_numeric(hit.iloc[0][samples], errors="coerce")
            s.name = gene
            rows.append(s)
        mat = pd.DataFrame(rows).T if rows else pd.DataFrame(index=samples)
        mat.index.name = "sample_id"
        return mat, pd.DataFrame(presence)

    human, presence_h = grab(ALL_GENES, "H")
    mouse, presence_m = grab(MOUSE_CONTRAST, "M")
    presence = pd.concat([presence_h, presence_m], ignore_index=True)
    return human, mouse, presence, samples


def collapse_replicates(human: pd.DataFrame, mouse: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Average log2 across the three paper QA replicate pairs. Skip NA in the mean."""
    rep_cols = {c for pair in REPLICATE_PAIRS.values() for c in pair}
    missing = sorted(rep_cols - set(human.index))
    if missing:
        raise SystemExit(f"expected technical replicates not in matrix: {missing}")

    def collapse(mat: pd.DataFrame) -> pd.DataFrame:
        keep = mat.drop(index=list(rep_cols), errors="ignore")
        extra = []
        audit = []
        for base, cols in REPLICATE_PAIRS.items():
            block = mat.loc[cols]
            averaged = block.mean(axis=0, skipna=True)
            averaged.name = base
            extra.append(averaged)
            for gene in mat.columns:
                vals = pd.to_numeric(block[gene], errors="coerce")
                audit.append(
                    {
                        "model": base,
                        "gene": gene,
                        "n_replicates_observed": int(vals.notna().sum()),
                        "replicate_values": ";".join(
                            f"{idx}={vals.loc[idx]:.4f}" if pd.notna(vals.loc[idx]) else f"{idx}=NA"
                            for idx in cols
                        ),
                    }
                )
        out = pd.concat([keep, pd.DataFrame(extra)])
        out.index.name = "model"
        return out, pd.DataFrame(audit)

    human_m, audit = collapse(human)
    mouse_m, _ = collapse(mouse)
    return human_m, mouse_m, audit


def attach_clinical(models: pd.DataFrame, clin: pd.DataFrame) -> pd.DataFrame:
    clin = clin.copy()
    clin["phlcid"] = clin["phlcid"].astype(str).str.strip()
    keep = clin.set_index("phlcid")[["histology", "proteome subtypes", "model", "sex", "pathStage"]]
    keep = keep.rename(columns={"proteome subtypes": "proteotype", "model": "pdx_stability"})
    out = models.join(keep, how="left")
    if out["histology"].isna().any():
        missing = out.index[out["histology"].isna()].tolist()
        raise SystemExit(f"models without clinical histology: {missing}")
    return out


def median_human_proteome(prot: pd.DataFrame, samples: list[str]) -> pd.Series:
    human = prot[prot["H/M"].astype(str).str.strip() == "H"]
    mat = human[samples].apply(pd.to_numeric, errors="coerce")
    med = mat.median(axis=0, skipna=True)
    med.index = samples
    # collapse replicates the same way
    rep_cols = [c for pair in REPLICATE_PAIRS.values() for c in pair]
    base = med.drop(index=rep_cols)
    for model, cols in REPLICATE_PAIRS.items():
        base.loc[model] = float(np.nanmean(med.loc[cols].to_numpy(dtype=float)))
    base.name = "median_human_log2"
    return base


def zscore(s: pd.Series) -> pd.Series:
    v = pd.to_numeric(s, errors="coerce")
    mu = v.mean(skipna=True)
    sd = v.std(skipna=True, ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return pd.Series(np.nan, index=s.index)
    return (v - mu) / sd


def add_scores(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    z_dn = pd.DataFrame({g: zscore(out[g]) for g in DNA_PK})
    z_st = pd.DataFrame({g: zscore(out[g]) for g in STING})
    out["DNAPK_score"] = z_dn.mean(axis=1, skipna=False)
    n_st = z_st.notna().sum(axis=1)
    sting = z_st.mean(axis=1, skipna=True)
    sting = sting.where(n_st >= 2)
    out["STING_score"] = sting
    out["n_sting_members"] = n_st
    return out


def annotate_q(df: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    parts = []
    for _, sub in df.groupby(group_cols, dropna=False):
        sub = sub.copy()
        sub["q_bh"] = bh(sub["p"].tolist())
        parts.append(sub)
    return pd.concat(parts, ignore_index=True) if parts else df


def run_spearman_block(
    df: pd.DataFrame,
    cohort: str,
    analysis: str,
    genes: list[str],
    partial: bool,
    include_scores: bool,
) -> list[dict]:
    rows = []
    for gene in genes:
        if partial:
            rec = partial_spearman(df["CLDN4"], df[gene], df["median_human_log2"])
        else:
            rec = spearman(df["CLDN4"], df[gene])
        rec.update(
            {
                "cohort": cohort,
                "analysis": analysis,
                "family": FAMILY[gene],
                "partner": gene,
                "species": "human",
            }
        )
        rows.append(rec)
    if not include_scores or partial:
        return rows
    for score, family in (("DNAPK_score", "DNA-PK"), ("STING_score", "STING")):
        rec = spearman(df["CLDN4"], df[score])
        rec.update(
            {
                "cohort": cohort,
                "analysis": analysis,
                "family": family,
                "partner": score,
                "species": "human",
                "note": rec["note"] + "; mean of within-cohort z-scores (composite, not in protein FDR)",
            }
        )
        rows.append(rec)
    return rows


def run_detection_block(df: pd.DataFrame, cohort: str, genes: list[str]) -> list[dict]:
    rows = []
    det = df["CLDN4"].notna()
    for gene in genes:
        rec = glass_rbc(df.loc[det, gene], df.loc[~det, gene])
        rec.update(
            {
                "cohort": cohort,
                "analysis": "detection_contrast",
                "family": FAMILY[gene],
                "partner": gene,
                "species": "human",
            }
        )
        rows.append(rec)
    return rows


def pearson(x: pd.Series, y: pd.Series) -> dict:
    d = pd.concat([x, y], axis=1).dropna()
    n = int(len(d))
    rec = {"n_pairwise": n, "rho": np.nan, "p": np.nan, "tested": False, "note": "not tested"}
    if n < 4 or d.iloc[:, 0].nunique() < 2 or d.iloc[:, 1].nunique() < 2:
        rec["note"] = f"n={n}; Pearson not computed"
        return rec
    r, p = stats.pearsonr(d.iloc[:, 0], d.iloc[:, 1])
    rec.update(
        {
            "rho": float(r),
            "p": float(p),
            "tested": True,
            "note": "two-sided Pearson, same complete cases as the Spearman",
        }
    )
    return rec


def tukey_mask(s: pd.Series, k: float = 1.5) -> pd.Series:
    q1 = s.quantile(0.25)
    q3 = s.quantile(0.75)
    iqr = q3 - q1
    return (s < q1 - k * iqr) | (s > q3 + k * iqr)


def cooks_flag(y: pd.Series, x: pd.Series) -> pd.Series:
    """Cook's distance from OLS y ~ x. Flag if D > 4/n. Index follows y/x after dropna alignment."""
    d = pd.concat([y, x], axis=1).dropna()
    d.columns = ["y", "x"]
    n = len(d)
    flag = pd.Series(False, index=d.index)
    if n < 6:
        return flag
    X = np.column_stack([np.ones(n), d["x"].to_numpy(dtype=float)])
    yy = d["y"].to_numpy(dtype=float)
    beta, *_ = np.linalg.lstsq(X, yy, rcond=None)
    resid = yy - X @ beta
    pdim = 2
    dof = n - pdim
    if dof <= 0:
        return flag
    mse = float(np.sum(resid**2) / dof)
    if mse <= 0:
        return flag
    hat = np.diag(X @ np.linalg.pinv(X.T @ X) @ X.T)
    cook = resid**2 / (pdim * mse) * (hat / np.maximum(1 - hat, 1e-12) ** 2)
    flag.loc[:] = cook > (4.0 / n)
    return flag


def luad_dnapk_sweep(luad: pd.DataFrame, clin: pd.DataFrame) -> pd.DataFrame:
    """One sweep of the LUAD CLDN4–DNA-PK trend. Does not impute.

    Pre-specified tests that decide the call (counts_toward_final_call):
    complete-case Spearman, complete-case Pearson, Tukey 1.5-IQR deletion,
    Cook's D > 4/n deletion, and within-label Spearman only when that label
    has n>=10. Smaller subtype cells and leave-one-out minima are written
    out and do not decide the call.
    """
    clin = clin.copy()
    clin["phlcid"] = clin["phlcid"].astype(str).str.strip()
    clin = clin.set_index("phlcid")
    partners = DNA_PK + ["DNAPK_score"]
    cc = luad.dropna(subset=["CLDN4", *DNA_PK]).copy()
    cc["transcriptome_subtype"] = clin.reindex(cc.index.astype(str))["Transcriptome subtype"]
    rows: list[dict] = []

    def add(rec: dict, **extra) -> None:
        rec = dict(rec)
        rec.update(extra)
        rows.append(rec)

    # WHO histologic pattern is not in this workbook. Record the absence.
    add(
        {
            "n_pairwise": 0,
            "rho": np.nan,
            "p": np.nan,
            "tested": False,
            "note": "WHO growth pattern (lepidic/acinar/papillary/micropapillary/solid) is not a column in Supplementary Data 1; not tested",
        },
        sweep="histology_subtype",
        partner="",
        stratum="WHO_pattern",
        n_stratum=0,
        n_dropped=0,
        dropped_models="",
        counts_toward_final_call=False,
    )

    for partner in partners:
        sp = spearman(cc["CLDN4"], cc[partner])
        add(
            sp,
            sweep="complete_case_spearman",
            partner=partner,
            stratum="LUAD",
            n_stratum=int(len(cc)),
            n_dropped=0,
            dropped_models="",
            counts_toward_final_call=True,
        )
        pe = pearson(cc["CLDN4"], cc[partner])
        add(
            pe,
            sweep="complete_case_pearson",
            partner=partner,
            stratum="LUAD",
            n_stratum=int(len(cc)),
            n_dropped=0,
            dropped_models="",
            counts_toward_final_call=True,
        )

        flag = tukey_mask(cc["CLDN4"]) | tukey_mask(cc[partner])
        dropped = cc.index[flag].astype(str).tolist()
        kept = cc.loc[~flag]
        sp = spearman(kept["CLDN4"], kept[partner])
        sp["note"] = sp["note"] + "; dropped Tukey 1.5 IQR outliers on CLDN4 or the partner"
        add(
            sp,
            sweep="outlier_tukey",
            partner=partner,
            stratum="LUAD",
            n_stratum=int(len(cc)),
            n_dropped=int(flag.sum()),
            dropped_models=",".join(dropped),
            counts_toward_final_call=True,
        )

        cflag = cooks_flag(cc["CLDN4"], cc[partner])
        dropped = cflag.index[cflag].astype(str).tolist()
        kept = cc.drop(index=cflag.index[cflag])
        sp = spearman(kept["CLDN4"], kept[partner])
        sp["note"] = sp["note"] + "; Cook's D > 4/n from OLS CLDN4 ~ partner, none imputed"
        add(
            sp,
            sweep="outlier_cook",
            partner=partner,
            stratum="LUAD",
            n_stratum=int(len(cc)),
            n_dropped=int(cflag.sum()),
            dropped_models=",".join(dropped),
            counts_toward_final_call=True,
        )

        loo = []
        for model in cc.index:
            sub = cc.drop(index=model)
            rec = spearman(sub["CLDN4"], sub[partner])
            loo.append((rec["p"], rec["rho"], str(model)))
        loo_ok = [t for t in loo if np.isfinite(t[0])]
        loo_ok.sort()
        n_under = int(sum(t[0] <= 0.05 for t in loo_ok))
        min_p, min_rho, min_model = loo_ok[0]
        add(
            {
                "n_pairwise": int(len(cc) - 1),
                "rho": float(min_rho),
                "p": float(min_p),
                "tested": True,
                "note": (
                    f"minimum leave-one-out Spearman p across {len(loo_ok)} deletions; "
                    f"{n_under} deletions have p<=0.05. This minimum is not a test and does not decide the call."
                ),
            },
            sweep="leave_one_out_min_p",
            partner=partner,
            stratum="LUAD",
            n_stratum=int(len(cc)),
            n_dropped=1,
            dropped_models=min_model,
            counts_toward_final_call=False,
            n_loo_p_le_005=n_under,
        )

    # Published within-LUAD labels. n>=10 counts; smaller cells are disclosed only.
    for col, sweep_name in (
        ("proteotype", "proteotype"),
        ("transcriptome_subtype", "transcriptome_subtype"),
    ):
        for level, sub in cc.groupby(cc[col].fillna("NA").astype(str)):
            counts = int(len(sub)) >= 10
            for partner in DNA_PK:
                sp = spearman(sub["CLDN4"], sub[partner])
                sp["note"] = sp["note"] + (
                    "; within-LUAD label, n>=10, counts toward the call"
                    if counts
                    else "; within-LUAD label, n<10, disclosed and not used to decide the call"
                )
                add(
                    sp,
                    sweep=sweep_name,
                    partner=partner,
                    stratum=str(level),
                    n_stratum=int(len(sub)),
                    n_dropped=0,
                    dropped_models="",
                    counts_toward_final_call=counts,
                )

    out = pd.DataFrame(rows)
    # BH within each disclosed subtype family, proteins only.
    for sweep_name in ("proteotype", "transcriptome_subtype"):
        mask = out["sweep"] == sweep_name
        if mask.any():
            out.loc[mask, "q_bh"] = bh(out.loc[mask, "p"].tolist())
    if "q_bh" not in out.columns:
        out["q_bh"] = np.nan
    out["q_bh"] = out["q_bh"].astype(float)
    return out


def final_null_from_sweep(sweep: pd.DataFrame) -> dict:
    deciding = sweep[sweep["counts_toward_final_call"] & sweep["tested"]].copy()
    if deciding.empty:
        return {"final_call": "not called", "n_deciding_tests": 0, "min_p": None}
    min_row = deciding.loc[deciding["p"].idxmin()]
    all_above = bool((deciding["p"] > 0.05).all())
    return {
        "final_call": "FINAL null" if all_above else "not null",
        "n_deciding_tests": int(len(deciding)),
        "min_p": float(min_row["p"]),
        "min_p_sweep": str(min_row["sweep"]),
        "min_p_partner": str(min_row["partner"]),
        "min_p_stratum": str(min_row["stratum"]),
        "min_p_n": int(min_row["n_pairwise"]),
        "min_p_rho": float(min_row["rho"]),
        "rule": (
            "FINAL null only if every pre-specified deciding test has p>0.05: "
            "LUAD complete-case Spearman, complete-case Pearson, Tukey outlier deletion, "
            "Cook outlier deletion, and within-LUAD label Spearman with n>=10. "
            "Leave-one-out minima and n<10 subtype cells are disclosed and do not decide."
        ),
    }


def plot(models: pd.DataFrame, pairs: pd.DataFrame, path: Path) -> None:
    panels = [
        ("DNAPK_score", "DNA-PK score"),
        ("PRKDC", "PRKDC"),
        ("STING_score", "STING score"),
        ("STING1", "STING1"),
    ]
    fig, axes = plt.subplots(2, 4, figsize=(11.2, 6.2), sharex=False, sharey=False)
    for i, cohort in enumerate(("LUAD", "LUSC")):
        d = models[models["histology"] == cohort]
        for j, (col, label) in enumerate(panels):
            ax = axes[i, j]
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            sub = d[["CLDN4", col]].dropna()
            row = pairs[
                (pairs.cohort == cohort)
                & (pairs.analysis == "primary_histology")
                & (pairs.partner == col)
            ]
            rho = p = n = np.nan
            if len(row):
                rho = row.iloc[0]["rho"]
                p = row.iloc[0]["p"]
                n = int(row.iloc[0]["n_pairwise"])
            if len(sub) == 0:
                ax.text(0.5, 0.5, "no pairwise models", ha="center", va="center", transform=ax.transAxes)
            else:
                ax.scatter(sub["CLDN4"], sub[col], s=18, c="#3C5488", alpha=0.85, linewidths=0)
            title_p = f"{p:.3g}" if np.isfinite(p) else "NA"
            title_r = f"{rho:+.2f}" if np.isfinite(rho) else "NA"
            ax.set_title(f"{cohort}  CLDN4 vs {label}\nρ={title_r}  p={title_p}  n={n}", fontsize=8)
            if i == 1:
                ax.set_xlabel("CLDN4 log2 protein", fontsize=8)
            if j == 0:
                ax.set_ylabel(label if col.endswith("score") else f"{label} log2", fontsize=8)
            ax.tick_params(labelsize=7)
    fig.suptitle(
        "Mirhadi NSCLC PDX TMT — human tumor protein (pairwise complete)",
        fontsize=11,
        y=1.02,
    )
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def n_table(models: pd.DataFrame, presence: pd.DataFrame) -> pd.DataFrame:
    rows = []
    cohorts = {
        "LUAD": models["histology"] == "LUAD",
        "LUSC": models["histology"] == "LUSC",
        "LUAD_plus_LUSC": models["histology"].isin(["LUAD", "LUSC"]),
        "other_histology": ~models["histology"].isin(["LUAD", "LUSC"]),
        "all_models": pd.Series(True, index=models.index),
    }
    for name, mask in cohorts.items():
        sub = models.loc[mask]
        rec = {
            "cohort": name,
            "n_models": int(len(sub)),
            "n_CLDN4_quantified": int(sub["CLDN4"].notna().sum()),
            "n_CLDN4_missing": int(sub["CLDN4"].isna().sum()),
        }
        for gene in PRIMARY + EXTENDED:
            rec[f"n_{gene}"] = int(sub[gene].notna().sum())
            both = sub[["CLDN4", gene]].dropna()
            rec[f"n_CLDN4_and_{gene}"] = int(len(both))
        rows.append(rec)
    # column-level presence before collapse is already in presence.tsv
    _ = presence
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, default=Path("data/pdx_cldn4_dnapk_sting"))
    ap.add_argument("--outdir", type=Path, default=Path("methods/pdx_cldn4_dnapk_sting"))
    args = ap.parse_args()
    xlsx = args.data / XLSX_NAME
    if not xlsx.exists():
        raise SystemExit(f"missing {xlsx}; run download.py first")

    prot = pd.read_excel(xlsx, sheet_name=PROT_SHEET)
    clin = pd.read_excel(xlsx, sheet_name=CLIN_SHEET)
    human_cols, mouse_cols, presence, samples = load_protein_rows(prot)
    if "CLDN4" not in human_cols.columns:
        raise SystemExit("human CLDN4 row missing — stop, do not substitute RNA or mouse CLDN4")
    for gene in PRIMARY:
        if gene not in human_cols.columns:
            raise SystemExit(f"human {gene} row missing")

    human, mouse, audit = collapse_replicates(human_cols, mouse_cols)
    med = median_human_proteome(prot, samples)
    models = attach_clinical(human, clin)
    models = models.join(med)
    models = add_scores(models)
    # mouse STING / DNA-PK kept beside the model for the stroma contrast only
    mouse = mouse.add_prefix("mouse_")
    models = models.join(mouse)

    # Scores were fit on all histologies together above. Refit within the
    # analysis cohort so LUAD z-scores do not use LUSC means.
    def with_cohort_scores(df: pd.DataFrame) -> pd.DataFrame:
        return add_scores(df.drop(columns=["DNAPK_score", "STING_score", "n_sting_members"], errors="ignore"))

    spearman_rows: list[dict] = []
    detect_rows: list[dict] = []
    for cohort, hist in (("LUAD", "LUAD"), ("LUSC", "LUSC")):
        sub = with_cohort_scores(models[models["histology"] == hist].copy())
        spearman_rows += run_spearman_block(
            sub, cohort, "primary_histology", PRIMARY, partial=False, include_scores=True
        )
        spearman_rows += run_spearman_block(
            sub, cohort, "partial_median", PRIMARY, partial=True, include_scores=False
        )
        spearman_rows += run_spearman_block(
            sub, cohort, "extended", EXTENDED, partial=False, include_scores=False
        )
        detect_rows += run_detection_block(sub, cohort, PRIMARY + ["IFI16", "STAT1"])

    pooled = with_cohort_scores(models[models["histology"].isin(["LUAD", "LUSC"])].copy())
    spearman_rows += run_spearman_block(
        pooled, "LUAD_plus_LUSC", "pooled_sensitivity", PRIMARY, partial=False, include_scores=True
    )

    # Stroma contrast: mouse protein vs human CLDN4. Not a primary test.
    stroma_rows = []
    for cohort, df in (
        ("LUAD", models[models.histology == "LUAD"]),
        ("LUSC", models[models.histology == "LUSC"]),
        ("LUAD_plus_LUSC", models[models.histology.isin(["LUAD", "LUSC"])]),
    ):
        for gene in MOUSE_CONTRAST:
            col = f"mouse_{gene}"
            if col not in df.columns:
                rec = {
                    "cohort": cohort,
                    "partner": gene,
                    "species": "mouse",
                    "n_pairwise": 0,
                    "rho": np.nan,
                    "p": np.nan,
                    "tested": False,
                    "note": "mouse row absent; not used as a human stand-in",
                }
            else:
                rec = spearman(df["CLDN4"], df[col])
                rec.update(
                    {
                        "cohort": cohort,
                        "partner": gene,
                        "species": "mouse",
                        "note": rec["note"] + "; murine stroma protein, excluded from the human test",
                    }
                )
            stroma_rows.append(rec)

    spearman_df = pd.DataFrame(spearman_rows)
    # FDR within cohort × analysis, proteins only (scores excluded).
    protein_mask = ~spearman_df["partner"].str.endswith("_score")
    scored = spearman_df[protein_mask].copy()
    scored = annotate_q(scored, ["cohort", "analysis"])
    scores = spearman_df[~protein_mask].copy()
    scores["q_bh"] = np.nan
    spearman_df = pd.concat([scored, scores], ignore_index=True)

    detect_df = pd.DataFrame(detect_rows)
    detect_df = annotate_q(detect_df, ["cohort"])

    # sample table: human proteins + scores fit within histology for LUAD/LUSC,
    # and the all-model z-score only as a stored column for other histologies.
    scored_parts = []
    for hist, sub in models.groupby("histology"):
        if hist in ("LUAD", "LUSC"):
            scored_parts.append(with_cohort_scores(sub.copy()))
        else:
            scored_parts.append(sub.copy())
    sample = pd.concat(scored_parts).sort_index()
    sample_out = sample.reset_index().rename(columns={"index": "model"})
    if "model" not in sample_out.columns:
        sample_out = sample.reset_index()
    # join may name the index 'model' already
    if sample_out.columns[0] != "model":
        sample_out = sample_out.rename(columns={sample_out.columns[0]: "model"})

    presence_model = []
    for gene in ALL_GENES:
        s = models[gene]
        presence_model.append(
            {
                "gene": gene,
                "species": "human",
                "family": "predictor" if gene == "CLDN4" else FAMILY[gene],
                "n_models": int(len(models)),
                "n_observed_models": int(s.notna().sum()),
                "n_missing_models": int(s.isna().sum()),
                "entry": presence.loc[
                    (presence.gene == gene) & (presence.species_filter == "H"), "entry"
                ].iloc[0],
            }
        )
    presence_out = pd.DataFrame(presence_model)

    outdir = args.outdir / "results"
    outdir.mkdir(parents=True, exist_ok=True)
    spearman_df.to_csv(outdir / "spearman.tsv", sep="\t", index=False)
    detect_df.to_csv(outdir / "detection.tsv", sep="\t", index=False)
    pd.DataFrame(stroma_rows).to_csv(outdir / "stroma_contrast.tsv", sep="\t", index=False)
    presence_out.to_csv(outdir / "presence.tsv", sep="\t", index=False)
    n_table(models, presence).to_csv(outdir / "n_table.tsv", sep="\t", index=False)
    audit.to_csv(outdir / "replicate_audit.tsv", sep="\t", index=False)
    keep_cols = [
        "model",
        "histology",
        "proteotype",
        "pdx_stability",
        "median_human_log2",
        "CLDN4",
        *PRIMARY,
        *EXTENDED,
        "DNAPK_score",
        "STING_score",
        "n_sting_members",
    ]
    sample_out[keep_cols].to_csv(outdir / "sample_scores.tsv", sep="\t", index=False)

    luad_scored = with_cohort_scores(models[models["histology"] == "LUAD"].copy())
    sweep = luad_dnapk_sweep(luad_scored, clin)
    verdict = final_null_from_sweep(sweep)
    sweep.to_csv(outdir / "luad_dnapk_sweep.tsv", sep="\t", index=False)

    primary = spearman_df[spearman_df.analysis == "primary_histology"]
    plot(pd.concat([
        with_cohort_scores(models[models.histology == "LUAD"].copy()),
        with_cohort_scores(models[models.histology == "LUSC"].copy()),
    ]), primary, outdir / "fig_cldn4_vs_dnapk_sting.png")

    def pack(df: pd.DataFrame) -> list[dict]:
        rows = []
        for rec in df.to_dict(orient="records"):
            clean = {}
            for k, v in rec.items():
                if isinstance(v, (np.floating, float)):
                    clean[k] = None if not np.isfinite(v) else float(v)
                elif isinstance(v, (np.integer,)):
                    clean[k] = int(v)
                elif isinstance(v, (np.bool_,)):
                    clean[k] = bool(v)
                else:
                    clean[k] = v
            rows.append(clean)
        return rows

    summary = {
        "dataset": "Mirhadi et al. Nat Commun 2022 Supplementary Data 1",
        "accession": "PXD016579",
        "doi": "10.1038/s41467-022-29444-9",
        "matrix": "Normalized log2 Prot Quant, human rows, author tumor/stroma normalization",
        "n_tmt_columns": len(samples),
        "n_models_after_replicate_average": int(len(models)),
        "replicate_pairs_averaged": list(REPLICATE_PAIRS),
        "primary_partners": PRIMARY,
        "extended_partners": EXTENDED,
        "n_CLDN4_LUAD": int(models.loc[models.histology == "LUAD", "CLDN4"].notna().sum()),
        "n_CLDN4_LUSC": int(models.loc[models.histology == "LUSC", "CLDN4"].notna().sum()),
        "n_models_LUAD": int((models.histology == "LUAD").sum()),
        "n_models_LUSC": int((models.histology == "LUSC").sum()),
        "primary_spearman": pack(primary),
        "luad_dnapk_final": verdict,
        "rule": "No RNA stand-in. No mouse-stroma stand-in. No imputation of CLDN4 NA. LUAD and LUSC not pooled in the primary test.",
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: summary[k] for k in summary if k != "primary_spearman"}, indent=2))
    deciding = sweep[sweep["counts_toward_final_call"] & sweep["tested"]]
    print(deciding[["sweep", "partner", "stratum", "n_pairwise", "rho", "p"]].to_string(index=False))
    print("VERDICT", verdict["final_call"], "min_p", verdict.get("min_p"))
    show = primary[primary.partner.isin(PRIMARY + ["DNAPK_score", "STING_score"])][
        ["cohort", "partner", "n_pairwise", "rho", "p", "q_bh"]
    ]
    print(show.to_string(index=False))


if __name__ == "__main__":
    main()
