#!/usr/bin/env python3
"""IMvigor210: TACSTD2 and CLDN4 vs confirmed ORR and OS.

Cohort: metastatic urothelial carcinoma treated with atezolizumab
(anti-PD-L1). Public RNA + clinical tables from IMvigor210CoreBiologies
v1.0.0 (Mariathasan et al., Nature 2018).

Primary endpoints (pre-specified):
  1. Confirmed ORR: CR/PR vs SD/PD (author binaryResponse; NE excluded).
     Two-sided Mann-Whitney U and AUC on log2(TPM+1) for TACSTD2 and CLDN4.
  2. Overall survival: Cox PH on continuous log2(TPM+1) (per 1 SD) and
     median-split log-rank, same two genes.

Exploratory (not in the primary FDR set):
  - mean z-score of TACSTD2 + CLDN4
  - both-high (each gene >= its median) vs rest / vs both-low
  - CXCL9 and CD8A as positive-control immune genes
  - sizeFactor-normalized log2 counts (DESeq scale) as a scale check
  - IC Level + Sex + platinum-adjusted Cox

This is urothelial carcinoma, not lung. Association on a treated cohort
is not a treatment-by-gene predictive interaction.

Usage:
  python3 scripts/w200/IMvigor_both/download.py
  python3 scripts/w200/IMvigor_both/analyze.py
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
from scipy import stats
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[3]
DATA_DIR = ROOT / "data" / "IMvigor210"
OUT_DIR = ROOT / "results" / "w200" / "IMvigor_both"

GENES = {
    "TACSTD2": 4070,
    "CLDN4": 1364,
    "CD8A": 925,
    "CXCL9": 4283,
}
PRIMARY = ["TACSTD2", "CLDN4"]
CONTROLS = ["CD8A", "CXCL9"]
N_BOOT = 2000
SEED = 20260816


def bh(pvals: list[float]) -> list[float]:
    return multipletests(pvals, method="fdr_bh")[1].tolist()


def auc_mwu(values: np.ndarray, labels: np.ndarray) -> tuple[float, float, float]:
    """AUC for label==1 as the positive class, plus two-sided MWU p."""
    pos = values[labels == 1]
    neg = values[labels == 0]
    res = stats.mannwhitneyu(pos, neg, alternative="two-sided")
    auc = float(res.statistic / (len(pos) * len(neg)))
    return auc, float(res.pvalue), float(res.statistic)


def bootstrap_auc_ci(
    values: np.ndarray, labels: np.ndarray, n_boot: int = N_BOOT, seed: int = SEED
) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    pos_idx = np.where(labels == 1)[0]
    neg_idx = np.where(labels == 0)[0]
    aucs = np.empty(n_boot)
    for i in range(n_boot):
        idx = np.concatenate(
            [
                rng.choice(pos_idx, len(pos_idx), replace=True),
                rng.choice(neg_idx, len(neg_idx), replace=True),
            ]
        )
        aucs[i], _, _ = auc_mwu(values[idx], labels[idx])
    lo, hi = np.percentile(aucs, [2.5, 97.5])
    return float(lo), float(hi)


def rank_biserial(auc: float) -> float:
    """Cliff's delta / rank-biserial: 2*AUC - 1 (positive => responders higher)."""
    return 2.0 * auc - 1.0


def load_tables() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    pdata = pd.read_csv(DATA_DIR / "pData.csv", index_col=0)
    fdata = pd.read_csv(DATA_DIR / "fData.csv", index_col=0)
    counts = pd.read_csv(DATA_DIR / "counts.csv.gz", index_col=0)
    counts.index = counts.index.astype(int)
    fdata.index = fdata.index.astype(int)
    if counts.shape != (31286, 348):
        raise RuntimeError(f"Unexpected counts shape {counts.shape}")
    if pdata.shape[0] != 348:
        raise RuntimeError(f"Unexpected pData n={pdata.shape[0]}")
    missing = [eid for eid in GENES.values() if eid not in counts.index]
    if missing:
        raise RuntimeError(f"Missing Entrez IDs in counts: {missing}")
    return pdata, fdata, counts


def compute_tpm(counts: pd.DataFrame, fdata: pd.DataFrame) -> pd.DataFrame:
    length = fdata.reindex(counts.index)["length"].astype(float)
    if length.isna().any() or (length <= 0).any():
        raise RuntimeError("Invalid gene lengths for TPM")
    rpk = counts.div(length, axis=0)
    tpm = rpk.div(rpk.sum(axis=0), axis=1) * 1e6
    return tpm


def cox_one(df: pd.DataFrame, feature: str, duration="os", event="censOS") -> dict:
    tmp = df[[feature, duration, event]].dropna().copy()
    tmp[feature] = tmp[feature].astype(float)
    cph = CoxPHFitter()
    cph.fit(tmp, duration_col=duration, event_col=event)
    row = cph.summary.loc[feature]
    return {
        "n": int(len(tmp)),
        "n_events": int(tmp[event].sum()),
        "hr": float(row["exp(coef)"]),
        "hr_lo": float(row["exp(coef) lower 95%"]),
        "hr_hi": float(row["exp(coef) upper 95%"]),
        "p": float(row["p"]),
        "coef": float(row["coef"]),
    }


def cox_adjusted(df: pd.DataFrame, feature: str) -> dict | None:
    cols = [feature, "os", "censOS", "Sex", "Received platinum", "IC Level"]
    tmp = df[cols].dropna().copy()
    tmp = pd.get_dummies(tmp, columns=["Sex", "Received platinum", "IC Level"], drop_first=True)
    # lifelines needs numeric
    for c in tmp.columns:
        tmp[c] = tmp[c].astype(float)
    if tmp[feature].std(ddof=1) == 0 or len(tmp) < 50:
        return None
    cph = CoxPHFitter()
    cph.fit(tmp, duration_col="os", event_col="censOS")
    row = cph.summary.loc[feature]
    return {
        "n": int(len(tmp)),
        "n_events": int(tmp["censOS"].sum()),
        "hr": float(row["exp(coef)"]),
        "hr_lo": float(row["exp(coef) lower 95%"]),
        "hr_hi": float(row["exp(coef) upper 95%"]),
        "p": float(row["p"]),
    }


def median_logrank(df: pd.DataFrame, feature: str) -> dict:
    hi = df[feature] >= df[feature].median()
    res = logrank_test(
        df.loc[hi, "os"],
        df.loc[~hi, "os"],
        df.loc[hi, "censOS"],
        df.loc[~hi, "censOS"],
    )
    return {
        "n_high": int(hi.sum()),
        "n_low": int((~hi).sum()),
        "p": float(res.p_value),
        "test_stat": float(res.test_statistic),
    }


def fmt_p(p: float) -> str:
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.3f}"


def fmt_ci(lo: float, hi: float, digits: int = 2) -> str:
    return f"{lo:.{digits}f}–{hi:.{digits}f}"


def build_sample_table(pdata: pd.DataFrame, tpm: pd.DataFrame, counts: pd.DataFrame) -> pd.DataFrame:
    sf = pdata["sizeFactor"].astype(float)
    out = pdata.copy()
    out = out.rename(
        columns={
            "Best Confirmed Overall Response": "BOR",
            "binaryResponse": "binaryResponse",
            "ANONPT_ID": "patient_id",
        }
    )
    for name, eid in GENES.items():
        out[f"{name}_TPM"] = tpm.loc[eid].reindex(out.index).astype(float)
        out[f"{name}_log2TPM1"] = np.log2(out[f"{name}_TPM"] + 1.0)
        out[f"{name}_log2SF"] = np.log2(counts.loc[eid].reindex(out.index).astype(float) / sf + 1.0)
    z_t = stats.zscore(out["TACSTD2_log2TPM1"], ddof=1)
    z_c = stats.zscore(out["CLDN4_log2TPM1"], ddof=1)
    out["combined_z"] = 0.5 * (z_t + z_c)
    med_t = out["TACSTD2_log2TPM1"].median()
    med_c = out["CLDN4_log2TPM1"].median()
    out["TACSTD2_high"] = out["TACSTD2_log2TPM1"] >= med_t
    out["CLDN4_high"] = out["CLDN4_log2TPM1"] >= med_c
    out["both_high"] = out["TACSTD2_high"] & out["CLDN4_high"]
    out["both_low"] = (~out["TACSTD2_high"]) & (~out["CLDN4_high"])
    out["quadrant"] = np.select(
        [
            out["both_high"],
            out["TACSTD2_high"] & ~out["CLDN4_high"],
            ~out["TACSTD2_high"] & out["CLDN4_high"],
            out["both_low"],
        ],
        ["HH", "HL", "LH", "LL"],
        default="NA",
    )
    # One RNA library is duplicated for patient 10285 (both NE, identical OS).
    # Keep the larger-sizeFactor library for patient-level OS.
    out["keep_patient"] = True
    dup = out[out.duplicated("patient_id", keep=False)]
    if len(dup):
        drop_idx = (
            dup.sort_values("sizeFactor", ascending=False)
            .groupby("patient_id", sort=False)
            .tail(1)
            .index
        )
        out.loc[drop_idx, "keep_patient"] = False
    out["orr_evaluable"] = out["binaryResponse"].isin(["CR/PR", "SD/PD"])
    out["responder"] = np.where(
        out["binaryResponse"] == "CR/PR",
        1,
        np.where(out["binaryResponse"] == "SD/PD", 0, np.nan),
    )
    return out


def orr_row(df: pd.DataFrame, feature: str, label: str) -> dict:
    sub = df.loc[df["orr_evaluable"], [feature, "responder"]].dropna()
    y = sub["responder"].to_numpy(dtype=int)
    x = sub[feature].to_numpy(dtype=float)
    auc, p, u = auc_mwu(x, y)
    lo, hi = bootstrap_auc_ci(x, y)
    pos = x[y == 1]
    neg = x[y == 0]
    return {
        "feature": label,
        "n_CRPR": int(y.sum()),
        "n_SDPD": int((y == 0).sum()),
        "median_CRPR": float(np.median(pos)),
        "median_SDPD": float(np.median(neg)),
        "mwu_U": float(u),
        "mwu_p": p,
        "auc_CRPR_higher": auc,
        "auc_ci_lo": lo,
        "auc_ci_hi": hi,
        "rank_biserial": rank_biserial(auc),
    }


def fisher_both(df: pd.DataFrame, col: str, vs: str) -> dict:
    sub = df.loc[df["orr_evaluable"]].copy()
    if vs == "rest":
        a = sub[col].astype(bool)
        b = ~a
        name = f"{col}_vs_rest"
    elif vs == "both_low":
        a = sub["both_high"]
        b = sub["both_low"]
        name = "both_high_vs_both_low"
        sub = sub.loc[a | b]
        a = sub["both_high"]
        b = sub["both_low"]
    else:
        raise ValueError(vs)
    tab = np.array(
        [
            [(~a & (sub["responder"] == 0)).sum(), (~a & (sub["responder"] == 1)).sum()],
            [(a & (sub["responder"] == 0)).sum(), (a & (sub["responder"] == 1)).sum()],
        ],
        dtype=int,
    )
    # rows: reference / high; cols: SD/PD / CR/PR
    odds, p = stats.fisher_exact(tab)
    rate_high = float(sub.loc[a, "responder"].mean())
    rate_ref = float(sub.loc[b, "responder"].mean())
    return {
        "contrast": name,
        "n_high": int(a.sum()),
        "n_ref": int(b.sum()),
        "orr_high": rate_high,
        "orr_ref": rate_ref,
        "odds_ratio": float(odds),
        "fisher_p": float(p),
        "table_SDPD_CRPR_ref_high": tab.tolist(),
    }


def write_figures(df: pd.DataFrame, out: Path) -> None:
    orr = df.loc[df["orr_evaluable"]].copy()
    orr["group"] = np.where(orr["responder"] == 1, "CR/PR", "SD/PD")
    colors = {"CR/PR": "#2A9D8F", "SD/PD": "#B0B0B0"}

    fig, axes = plt.subplots(2, 2, figsize=(8.2, 7.2))
    box_feats = [
        ("TACSTD2_log2TPM1", "TACSTD2"),
        ("CLDN4_log2TPM1", "CLDN4"),
        ("combined_z", "combined z (exploratory)"),
        ("CXCL9_log2TPM1", "CXCL9 (positive control)"),
    ]
    for ax, (feat, title) in zip(axes.ravel(), box_feats):
        data = [orr.loc[orr.group == g, feat].to_numpy() for g in ("SD/PD", "CR/PR")]
        bp = ax.boxplot(
            data, tick_labels=("SD/PD", "CR/PR"), patch_artist=True, widths=0.55
        )
        for patch, g in zip(bp["boxes"], ("SD/PD", "CR/PR")):
            patch.set_facecolor(colors[g])
            patch.set_alpha(0.85)
        rng = np.random.default_rng(SEED)
        for i, g in enumerate(("SD/PD", "CR/PR"), start=1):
            y = orr.loc[orr.group == g, feat].to_numpy()
            x = i + rng.uniform(-0.12, 0.12, size=len(y))
            ax.scatter(x, y, s=10, c="black", alpha=0.25, linewidths=0)
        auc, p, _ = auc_mwu(orr[feat].to_numpy(), orr["responder"].to_numpy(dtype=int))
        ax.set_title(f"{title}\nAUC={auc:.2f}  p={fmt_p(p)}", fontsize=10)
        ax.set_ylabel("log2(TPM+1)" if feat != "combined_z" else "mean z-score")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.suptitle("IMvigor210 ORR: TACSTD2 / CLDN4 (n=298 evaluable)", fontsize=12)
    fig.tight_layout()
    fig.savefig(out / "orr_boxplots.png", dpi=160)
    fig.savefig(out / "orr_boxplots.pdf")
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(8.6, 7.4))
    km_specs = [
        ("TACSTD2_log2TPM1", "TACSTD2 median"),
        ("CLDN4_log2TPM1", "CLDN4 median"),
        ("both_high", "both-high vs rest (exploratory)"),
        ("CXCL9_log2TPM1", "CXCL9 median (positive control)"),
    ]
    os_df = df.loc[df["keep_patient"]].copy()
    for ax, (feat, title) in zip(axes.ravel(), km_specs):
        kmf = KaplanMeierFitter()
        if feat == "both_high":
            mask_hi = os_df["both_high"]
            lab_hi, lab_lo = "both-high", "not both-high"
        else:
            mask_hi = os_df[feat] >= os_df[feat].median()
            lab_hi, lab_lo = "high", "low"
        for mask, lab, color in (
            (mask_hi, lab_hi, "#E76F51"),
            (~mask_hi, lab_lo, "#264653"),
        ):
            kmf.fit(
                os_df.loc[mask, "os"],
                os_df.loc[mask, "censOS"],
                label=f"{lab} n={int(mask.sum())}",
            )
            kmf.plot_survival_function(ax=ax, ci_show=False, color=color)
        res = logrank_test(
            os_df.loc[mask_hi, "os"],
            os_df.loc[~mask_hi, "os"],
            os_df.loc[mask_hi, "censOS"],
            os_df.loc[~mask_hi, "censOS"],
        )
        ax.set_title(f"{title}\nlog-rank p={fmt_p(float(res.p_value))}", fontsize=10)
        ax.set_xlabel("OS (months)")
        ax.set_ylabel("Survival")
        ax.set_ylim(0, 1.02)
        ax.legend(fontsize=7, frameon=False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.suptitle("IMvigor210 OS: median / both-high splits (n=347 patients)", fontsize=12)
    fig.tight_layout()
    fig.savefig(out / "os_km.png", dpi=160)
    fig.savefig(out / "os_km.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.2, 5.4))
    for g, color, label in (
        (0, "#B0B0B0", "SD/PD"),
        (1, "#2A9D8F", "CR/PR"),
    ):
        sub = orr.loc[orr["responder"] == g]
        ax.scatter(
            sub["TACSTD2_log2TPM1"],
            sub["CLDN4_log2TPM1"],
            s=18,
            c=color,
            alpha=0.75,
            label=f"{label} n={len(sub)}",
            edgecolors="none",
        )
    rho, rp = stats.spearmanr(orr["TACSTD2_log2TPM1"], orr["CLDN4_log2TPM1"])
    ax.set_xlabel("TACSTD2 log2(TPM+1)")
    ax.set_ylabel("CLDN4 log2(TPM+1)")
    ax.set_title(f"ORR-evaluable co-expression  Spearman ρ={rho:.2f} p={fmt_p(rp)}")
    ax.legend(frameon=False, fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(out / "tacstd2_cldn4_scatter.png", dpi=160)
    fig.savefig(out / "tacstd2_cldn4_scatter.pdf")
    plt.close(fig)


def write_readme(ctx: dict, out: Path) -> None:
    o = {r["feature"]: r for r in ctx["orr"]}
    s = {r["feature"]: r for r in ctx["os_cox_sd"]}
    k = {r["feature"]: r for r in ctx["os_logrank"]}
    bh_rest = ctx["both_high_orr"][0]
    bh_ll = ctx["both_high_orr"][1]
    t = o["TACSTD2_log2TPM1"]
    c = o["CLDN4_log2TPM1"]
    z = o["combined_z"]
    x9 = o["CXCL9_log2TPM1"]
    cd8 = o["CD8A_log2TPM1"]
    text = f"""# IMvigor210: TACSTD2 and CLDN4 vs ORR and OS

**Bottom line: no significant association for TACSTD2, CLDN4, or the
both-high (AND) gate.** In this public atezolizumab urothelial-carcinoma
cohort, neither gene (nor their mean z-score, nor TACSTD2-high ∩ CLDN4-high)
separated confirmed responders from non-responders, and neither associated
with overall survival. Immune positive controls (CXCL9, CD8A) do associate
in the expected direction, so the endpoints are not inert.

一句话结论：公开 IMvigor210（尿路上皮癌，atezolizumab）里 TACSTD2 / CLDN4
以及双高（AND）与 ORR、OS 均无显著关联；CXCL9 / CD8A 对照方向正确。

## Cohort and data

- Source: IMvigor210CoreBiologies v1.0.0 (Mariathasan et al., *Nature* 2018,
  PMID 29443960), Creative Commons 3.0. Official tarball:
  `http://research-pub.gene.com/IMvigor210CoreBiologies/packageVersions/IMvigor210CoreBiologies_1.0.0.tar.gz`
- Disease / drug: **metastatic urothelial carcinoma**, atezolizumab (anti-PD-L1).
  This is **not** a lung ICI cohort. Biopsy site is recorded (`Tissue`); 10 of
  348 libraries are from lung metastases of urothelial cancer.
- RNA: 31,286 genes × 348 libraries (IGIS Entrez annotation). Expression is
  **log2(TPM+1)** using package gene lengths. TACSTD2 = Entrez 4070,
  CLDN4 = 1364.
- ORR: author `binaryResponse` = confirmed **CR/PR** (n={t['n_CRPR']}) vs
  **SD/PD** (n={t['n_SDPD']}). **NE (n=50) excluded** from ORR tests. This is
  the package's own confirmed-ORR coding, not a recode.
- OS: `os` (months) and `censOS` (1 = death). All 348 libraries have OS;
  one patient contributed two NE libraries with identical OS — patient-level
  OS uses n={ctx['n_os_patients']} (larger sizeFactor library kept).
- Primary tests: TACSTD2 and CLDN4 only. BH-FDR is within endpoint
  (2 genes for ORR; 2 genes for continuous Cox OS). Combined z, both-high,
  and CXCL9/CD8A are labeled exploratory / control.

## ORR results (n=298 evaluable)

| Feature | Median CR/PR | Median SD/PD | MWU p | BH q (2 genes) | AUC (CR/PR higher) | AUC 95% CI | rank-biserial |
|---|---:|---:|---:|---:|---:|---|---:|
| TACSTD2 log2(TPM+1) | {t['median_CRPR']:.2f} | {t['median_SDPD']:.2f} | {fmt_p(t['mwu_p'])} | {fmt_p(t['bh_q'])} | {t['auc_CRPR_higher']:.2f} | {fmt_ci(t['auc_ci_lo'], t['auc_ci_hi'])} | {t['rank_biserial']:.2f} |
| CLDN4 log2(TPM+1) | {c['median_CRPR']:.2f} | {c['median_SDPD']:.2f} | {fmt_p(c['mwu_p'])} | {fmt_p(c['bh_q'])} | {c['auc_CRPR_higher']:.2f} | {fmt_ci(c['auc_ci_lo'], c['auc_ci_hi'])} | {c['rank_biserial']:.2f} |
| Combined z (exploratory) | {z['median_CRPR']:.2f} | {z['median_SDPD']:.2f} | {fmt_p(z['mwu_p'])} | — | {z['auc_CRPR_higher']:.2f} | {fmt_ci(z['auc_ci_lo'], z['auc_ci_hi'])} | {z['rank_biserial']:.2f} |
| CXCL9 (control) | {x9['median_CRPR']:.2f} | {x9['median_SDPD']:.2f} | {fmt_p(x9['mwu_p'])} | — | {x9['auc_CRPR_higher']:.2f} | {fmt_ci(x9['auc_ci_lo'], x9['auc_ci_hi'])} | {x9['rank_biserial']:.2f} |
| CD8A (control) | {cd8['median_CRPR']:.2f} | {cd8['median_SDPD']:.2f} | {fmt_p(cd8['mwu_p'])} | — | {cd8['auc_CRPR_higher']:.2f} | {fmt_ci(cd8['auc_ci_lo'], cd8['auc_ci_hi'])} | {cd8['rank_biserial']:.2f} |

AUC > 0.5 means **higher** expression in confirmed responders. TACSTD2 and
CLDN4 point very slightly that way (AUC {t['auc_CRPR_higher']:.2f} and
{c['auc_CRPR_higher']:.2f}); both 95% CIs include 0.5 and both p-values are
~0.2. That is **not** support for a TROP2/CLDN4-high → worse-ORR claim, and
it is also **not** a significant responder-high finding.

Both-high (each gene ≥ its cohort median) vs rest: ORR
{bh_rest['orr_high']*100:.1f}% ({bh_rest['n_high']} pts) vs
{bh_rest['orr_ref']*100:.1f}% ({bh_rest['n_ref']} pts), Fisher OR
{bh_rest['odds_ratio']:.2f}, p={fmt_p(bh_rest['fisher_p'])}.
Both-high vs both-low: ORR {bh_ll['orr_high']*100:.1f}% vs
{bh_ll['orr_ref']*100:.1f}%, Fisher p={fmt_p(bh_ll['fisher_p'])}.

TACSTD2 and CLDN4 co-express (Spearman ρ={ctx['spearman_rho']:.2f},
p={fmt_p(ctx['spearman_p'])} on ORR-evaluable samples).

## OS results (n={ctx['n_os_patients']} patients)

Cox PH on **per-1-SD** log2(TPM+1); HR > 1 = higher expression, shorter OS.

| Feature | n / events | HR per 1 SD | 95% CI | p | BH q (2 genes) | median log-rank p |
|---|---|---:|---|---:|---:|---:|
| TACSTD2 | {s['TACSTD2_log2TPM1_sd']['n']} / {s['TACSTD2_log2TPM1_sd']['n_events']} | {s['TACSTD2_log2TPM1_sd']['hr']:.2f} | {fmt_ci(s['TACSTD2_log2TPM1_sd']['hr_lo'], s['TACSTD2_log2TPM1_sd']['hr_hi'])} | {fmt_p(s['TACSTD2_log2TPM1_sd']['p'])} | {fmt_p(s['TACSTD2_log2TPM1_sd']['bh_q'])} | {fmt_p(k['TACSTD2_log2TPM1']['p'])} |
| CLDN4 | {s['CLDN4_log2TPM1_sd']['n']} / {s['CLDN4_log2TPM1_sd']['n_events']} | {s['CLDN4_log2TPM1_sd']['hr']:.2f} | {fmt_ci(s['CLDN4_log2TPM1_sd']['hr_lo'], s['CLDN4_log2TPM1_sd']['hr_hi'])} | {fmt_p(s['CLDN4_log2TPM1_sd']['p'])} | {fmt_p(s['CLDN4_log2TPM1_sd']['bh_q'])} | {fmt_p(k['CLDN4_log2TPM1']['p'])} |
| Combined z (exploratory) | {s['combined_z_sd']['n']} / {s['combined_z_sd']['n_events']} | {s['combined_z_sd']['hr']:.2f} | {fmt_ci(s['combined_z_sd']['hr_lo'], s['combined_z_sd']['hr_hi'])} | {fmt_p(s['combined_z_sd']['p'])} | — | {fmt_p(k['combined_z']['p'])} |
| CXCL9 (control) | {s['CXCL9_log2TPM1_sd']['n']} / {s['CXCL9_log2TPM1_sd']['n_events']} | {s['CXCL9_log2TPM1_sd']['hr']:.2f} | {fmt_ci(s['CXCL9_log2TPM1_sd']['hr_lo'], s['CXCL9_log2TPM1_sd']['hr_hi'])} | {fmt_p(s['CXCL9_log2TPM1_sd']['p'])} | — | {fmt_p(k['CXCL9_log2TPM1']['p'])} |
| CD8A (control) | {s['CD8A_log2TPM1_sd']['n']} / {s['CD8A_log2TPM1_sd']['n_events']} | {s['CD8A_log2TPM1_sd']['hr']:.2f} | {fmt_ci(s['CD8A_log2TPM1_sd']['hr_lo'], s['CD8A_log2TPM1_sd']['hr_hi'])} | {fmt_p(s['CD8A_log2TPM1_sd']['p'])} | — | {fmt_p(k['CD8A_log2TPM1']['p'])} |

Both-high vs rest log-rank p={fmt_p(ctx['both_high_os_p'])}
(n_high={ctx['both_high_os_n']}). IC/Sex/platinum-adjusted Cox HRs for
TACSTD2 and CLDN4 remain ~1 (see `os_cox_adjusted.csv`).

The exploratory **combined-z median** log-rank is p={fmt_p(k['combined_z']['p'])}
(high combined z has *longer* observed OS). That split is **not** a primary
result: the pre-specified continuous Cox for the same score is HR
{s['combined_z_sd']['hr']:.2f} (p={fmt_p(s['combined_z_sd']['p'])}). Do not
quote the median p-value as an OS hit.

## Honest caveats

- **Wrong disease for a lung-TROP2 claim.** IMvigor210 is urothelial
  carcinoma. A null here does not prove a null in NSCLC, and a hit here
  would not have transferred either.
- **Null is not proof of no effect.** Effects smaller than about AUC 0.60
  or OS HR ~0.85/1.18 per SD are compatible with these confidence intervals.
  What is *not* supported is a large, consistent TACSTD2/CLDN4–ORR/OS link
  in this public matrix.
- **Direction vs the usual TROP2-cold story.** Point estimates are slightly
  *higher* expression in responders and HR ≈ 1 for OS — the opposite of
  “TROP2/CLDN4-high → worse ICI,” and still non-significant.
- Bulk TPM. No purity residual in the primary test. Epithelial genes are
  diluted by stroma; that can mask a tumor-intrinsic effect (see other
  slices) but does not license reading this as a positive ICI biomarker.
- Association on an all-treated cohort. There is no chemotherapy/SOC arm
  in these public tables, so this is **not** a predictive interaction test.
- Follow-up is capped near 24.5 months; 50 NE patients are out of ORR.
- CXCL9/CD8A are controls that the assay/endpoints can move, not a license
  to hunt other genes post hoc.

## Reproduction

```bash
python3 -m pip install -r scripts/w200/IMvigor_both/requirements.txt
python3 scripts/w200/IMvigor_both/download.py
python3 scripts/w200/IMvigor_both/analyze.py
```

`download.py` pulls the official 122 MB tarball (or reuses a cached copy)
and uses `extract_cds.R` to dump `pData`, `fData`, and counts. Raw files
stay under `data/IMvigor210/` (gitignored). Seed 20260816 for bootstrap.

## Files

- `README.md` — this write-up (numbers filled by `analyze.py`)
- `summary.json` — machine-readable verdict and key stats
- `sample_data.csv` — per-library clinical + expression (no raw counts)
- `orr_stats.csv`, `os_cox.csv`, `os_cox_adjusted.csv`, `os_logrank.csv`
- `both_high_orr.csv`, `quadrant_orr.csv`, `spearman_tacstd2_cldn4.csv`
- `provenance.tsv` — source URL, expected size, sha256 if available
- `orr_boxplots.png/.pdf`, `os_km.png/.pdf`, `tacstd2_cldn4_scatter.png/.pdf`
"""
    (out / "README.md").write_text(text)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    pdata, fdata, counts = load_tables()
    tpm = compute_tpm(counts, fdata)
    df = build_sample_table(pdata, tpm, counts)

    sample_cols = [
        "patient_id",
        "keep_patient",
        "BOR",
        "binaryResponse",
        "responder",
        "orr_evaluable",
        "os",
        "censOS",
        "Sex",
        "IC Level",
        "TC Level",
        "Enrollment IC",
        "Immune phenotype",
        "Received platinum",
        "Tissue",
        "TCGA Subtype",
        "Lund2",
        "sizeFactor",
        "TACSTD2_TPM",
        "TACSTD2_log2TPM1",
        "TACSTD2_log2SF",
        "CLDN4_TPM",
        "CLDN4_log2TPM1",
        "CLDN4_log2SF",
        "CD8A_log2TPM1",
        "CXCL9_log2TPM1",
        "combined_z",
        "TACSTD2_high",
        "CLDN4_high",
        "both_high",
        "both_low",
        "quadrant",
    ]
    df[sample_cols].to_csv(OUT_DIR / "sample_data.csv")

    orr_features = [
        ("TACSTD2_log2TPM1", "TACSTD2_log2TPM1"),
        ("CLDN4_log2TPM1", "CLDN4_log2TPM1"),
        ("combined_z", "combined_z"),
        ("CXCL9_log2TPM1", "CXCL9_log2TPM1"),
        ("CD8A_log2TPM1", "CD8A_log2TPM1"),
        ("TACSTD2_log2SF", "TACSTD2_log2SF"),
        ("CLDN4_log2SF", "CLDN4_log2SF"),
    ]
    orr_rows = [orr_row(df, feat, lab) for feat, lab in orr_features]
    # BH on the two primary TPM tests only
    primary_idx = [i for i, r in enumerate(orr_rows) if r["feature"] in ("TACSTD2_log2TPM1", "CLDN4_log2TPM1")]
    q = bh([orr_rows[i]["mwu_p"] for i in primary_idx])
    for i, qi in zip(primary_idx, q):
        orr_rows[i]["bh_q"] = float(qi)
    for r in orr_rows:
        r.setdefault("bh_q", np.nan)
    pd.DataFrame(orr_rows).to_csv(OUT_DIR / "orr_stats.csv", index=False)

    os_df = df.loc[df["keep_patient"]].copy()
    cox_rows = []
    for feat in [
        "TACSTD2_log2TPM1",
        "CLDN4_log2TPM1",
        "combined_z",
        "CXCL9_log2TPM1",
        "CD8A_log2TPM1",
        "TACSTD2_log2SF",
        "CLDN4_log2SF",
    ]:
        raw = cox_one(os_df, feat)
        raw["feature"] = feat
        raw["scale"] = "per_unit"
        zfeat = feat + "_sd"
        os_df[zfeat] = stats.zscore(os_df[feat], ddof=1)
        sd = cox_one(os_df, zfeat)
        sd["feature"] = zfeat
        sd["scale"] = "per_1sd"
        cox_rows.extend([raw, sd])
    # BH on primary per-1SD TPM Cox
    prim_cox = [r for r in cox_rows if r["feature"] in ("TACSTD2_log2TPM1_sd", "CLDN4_log2TPM1_sd")]
    q = bh([r["p"] for r in prim_cox])
    for r, qi in zip(prim_cox, q):
        r["bh_q"] = float(qi)
    for r in cox_rows:
        r.setdefault("bh_q", np.nan)
    pd.DataFrame(cox_rows).to_csv(OUT_DIR / "os_cox.csv", index=False)

    adj_rows = []
    for feat in ["TACSTD2_log2TPM1", "CLDN4_log2TPM1", "CXCL9_log2TPM1"]:
        os_df[feat + "_sd"] = stats.zscore(os_df[feat], ddof=1)
        rec = cox_adjusted(os_df, feat + "_sd")
        if rec:
            rec["feature"] = feat + "_sd"
            adj_rows.append(rec)
    pd.DataFrame(adj_rows).to_csv(OUT_DIR / "os_cox_adjusted.csv", index=False)

    lr_rows = []
    for feat in [
        "TACSTD2_log2TPM1",
        "CLDN4_log2TPM1",
        "combined_z",
        "CXCL9_log2TPM1",
        "CD8A_log2TPM1",
    ]:
        rec = median_logrank(os_df, feat)
        rec["feature"] = feat
        lr_rows.append(rec)
    pd.DataFrame(lr_rows).to_csv(OUT_DIR / "os_logrank.csv", index=False)

    both_orr = [fisher_both(df, "both_high", "rest"), fisher_both(df, "both_high", "both_low")]
    pd.DataFrame(both_orr).to_csv(OUT_DIR / "both_high_orr.csv", index=False)

    quad = (
        df.loc[df["orr_evaluable"]]
        .groupby("quadrant")
        .agg(n=("responder", "size"), n_CRPR=("responder", "sum"), orr=("responder", "mean"))
        .reindex(["HH", "HL", "LH", "LL"])
    )
    quad.to_csv(OUT_DIR / "quadrant_orr.csv")

    orr = df.loc[df["orr_evaluable"]]
    rho, rp = stats.spearmanr(orr["TACSTD2_log2TPM1"], orr["CLDN4_log2TPM1"])
    rho_all, rp_all = stats.spearmanr(df["TACSTD2_log2TPM1"], df["CLDN4_log2TPM1"])
    pd.DataFrame(
        [
            {"subset": "ORR_evaluable", "n": int(len(orr)), "spearman_rho": float(rho), "p": float(rp)},
            {"subset": "all_libraries", "n": int(len(df)), "spearman_rho": float(rho_all), "p": float(rp_all)},
        ]
    ).to_csv(OUT_DIR / "spearman_tacstd2_cldn4.csv", index=False)

    bh_os = logrank_test(
        os_df.loc[os_df.both_high, "os"],
        os_df.loc[~os_df.both_high, "os"],
        os_df.loc[os_df.both_high, "censOS"],
        os_df.loc[~os_df.both_high, "censOS"],
    )

    sha_path = DATA_DIR / "tarball.sha256"
    sha = sha_path.read_text().strip() if sha_path.is_file() else "not_computed"
    pd.DataFrame(
        [
            {
                "item": "IMvigor210CoreBiologies_1.0.0.tar.gz",
                "url": (
                    "http://research-pub.gene.com/IMvigor210CoreBiologies/"
                    "packageVersions/IMvigor210CoreBiologies_1.0.0.tar.gz"
                ),
                "expected_bytes": 122127298,
                "sha256": sha.split()[0] if sha != "not_computed" else sha,
                "license": "CC BY 3.0",
                "citation": "Mariathasan et al. Nature 2018 PMID 29443960",
            }
        ]
    ).to_csv(OUT_DIR / "provenance.tsv", sep="\t", index=False)

    t_row = next(r for r in orr_rows if r["feature"] == "TACSTD2_log2TPM1")
    c_row = next(r for r in orr_rows if r["feature"] == "CLDN4_log2TPM1")
    t_cox = next(r for r in cox_rows if r["feature"] == "TACSTD2_log2TPM1_sd")
    c_cox = next(r for r in cox_rows if r["feature"] == "CLDN4_log2TPM1_sd")
    x_row = next(r for r in orr_rows if r["feature"] == "CXCL9_log2TPM1")

    verdict = (
        "NOT SUPPORTED. TACSTD2, CLDN4, combined z, and both-high "
        "(TACSTD2-high AND CLDN4-high) are non-significant for confirmed "
        "ORR and OS in public IMvigor210 (mUC, atezolizumab). CXCL9/CD8A "
        "positive controls associate in the expected direction."
    )
    summary = {
        "slice": "results/w200/IMvigor_both",
        "cohort": "IMvigor210 metastatic urothelial carcinoma, atezolizumab",
        "not_lung": True,
        "n_libraries": int(len(df)),
        "n_orr_evaluable": int(df["orr_evaluable"].sum()),
        "n_CRPR": int((df["responder"] == 1).sum()),
        "n_SDPD": int((df["responder"] == 0).sum()),
        "n_NE": int((df["BOR"] == "NE").sum()),
        "n_os_patients": int(os_df.shape[0]),
        "honest_verdict": verdict,
        "primary_orr": {
            "TACSTD2": {k: t_row[k] for k in ("mwu_p", "bh_q", "auc_CRPR_higher", "auc_ci_lo", "auc_ci_hi")},
            "CLDN4": {k: c_row[k] for k in ("mwu_p", "bh_q", "auc_CRPR_higher", "auc_ci_lo", "auc_ci_hi")},
        },
        "primary_os_hr_per_sd": {
            "TACSTD2": {k: t_cox[k] for k in ("hr", "hr_lo", "hr_hi", "p", "bh_q")},
            "CLDN4": {k: c_cox[k] for k in ("hr", "hr_lo", "hr_hi", "p", "bh_q")},
        },
        "positive_control_CXCL9_orr_p": x_row["mwu_p"],
        "both_high_orr_fisher_p": both_orr[0]["fisher_p"],
        "both_high_os_logrank_p": float(bh_os.p_value),
    }
    (OUT_DIR / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    ctx = {
        "orr": orr_rows,
        "os_cox_sd": [r for r in cox_rows if r["scale"] == "per_1sd"],
        "os_logrank": lr_rows,
        "both_high_orr": both_orr,
        "spearman_rho": float(rho),
        "spearman_p": float(rp),
        "n_os_patients": int(os_df.shape[0]),
        "both_high_os_p": float(bh_os.p_value),
        "both_high_os_n": int(os_df["both_high"].sum()),
    }
    write_readme(ctx, OUT_DIR)
    write_figures(df, OUT_DIR)
    print(json.dumps(summary, indent=2))
    print(f"Wrote {OUT_DIR}")


if __name__ == "__main__":
    main()
