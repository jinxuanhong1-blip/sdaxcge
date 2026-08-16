#!/usr/bin/env python3
"""Gide 2019 melanoma ICI: public RNA CLDN4 versus response (OR / n / p)."""

from __future__ import annotations

import hashlib
import json
import math
import urllib.request
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rdata
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test
from scipy import stats
from statsmodels.formula.api import logit
from statsmodels.stats.contingency_tables import Table2x2


ROOT = Path(__file__).resolve().parent
RAW = ROOT / "raw"
DATA = ROOT / "data"
TABLES = ROOT / "tables"
FIGURES = ROOT / "figures"
SEED = 2019
CLDN4_SYMBOL = "CLDN4"
CLDN4_ENSG = "ENSG00000189143.10"

GIDE_RDA = RAW / "Gide_PD1highCD8Tscore.rda"
GIDE_RDA_URL = (
    "https://raw.githubusercontent.com/Liulab/PD1highCD8Tscore/master/data/Gide.rda"
)
GIDE_RDA_SHA256 = (
    "170a6079505310ddc7cd88a1283bb415331d89f5a144675fcc96d1602e83f9a5"
)
ORCESTRA_CSV = DATA / "orcestra_icb_gide_cldn4.csv"
ORCESTRA_RDS_URL = "https://zenodo.org/record/7332096/files/ICB_Gide.rds?download=1"
ORCESTRA_RDS_SHA256 = (
    "4fba93efefa2f18bcc952d398eb34e5a085e71032d34cd0482f4329dfbd7b536"
)
ORCESTRA_DOI = "10.5281/zenodo.7332096"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_gide_rda() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    if not GIDE_RDA.exists():
        print(f"Downloading {GIDE_RDA_URL}")
        urllib.request.urlretrieve(GIDE_RDA_URL, GIDE_RDA)
    digest = sha256(GIDE_RDA)
    if digest != GIDE_RDA_SHA256:
        raise ValueError(f"Unexpected Gide.rda sha256: {digest}")


def load_liulab() -> pd.DataFrame:
    converted = rdata.conversion.convert(rdata.parser.parse_file(GIDE_RDA))
    bundle = converted["Gide"]
    rna = bundle["RNA"]
    clinic = bundle["clinic"].copy()
    genes = [str(g) for g in rna.coords["dim_0"].values]
    samples = [str(s) for s in rna.coords["dim_1"].values]
    hits = [i for i, g in enumerate(genes) if g == CLDN4_SYMBOL]
    if len(hits) != 1:
        raise ValueError(f"Expected one CLDN4 row; found {hits}")
    if list(clinic.index.astype(str)) != samples:
        raise ValueError("Liulab RNA sample order does not match clinic rownames")
    tpm = np.asarray(rna.isel(dim_0=hits[0]).values, dtype=float)
    if tpm.shape != (len(samples),):
        raise ValueError("Unexpected CLDN4 vector shape")
    out = clinic.copy()
    out.index = samples
    out.index.name = "sample_id"
    out["cldn4_tpm"] = tpm
    out["cldn4_log2_tpm_plus_1"] = np.log2(tpm + 1.0)
    out["patient_id"] = out["patient"].astype(str).str.replace(
        r"_(PRE|EDT)$", "", regex=True
    )
    out["visit"] = out["Alias"].astype(str)
    out["drug"] = out["drug"].astype(str)
    out["recist"] = out["Best.RECIST.response"].astype(str)
    out["paper_response"] = out["response"].astype(str)
    out["treatment_class"] = np.where(
        out["drug"].isin(["Pembrolizumab", "Nivolumab"]),
        "anti-PD-1",
        "anti-PD-1 + anti-CTLA-4",
    )
    return out


def load_orcestra() -> pd.DataFrame:
    if not ORCESTRA_CSV.exists():
        raise FileNotFoundError(ORCESTRA_CSV)
    orc = pd.read_csv(ORCESTRA_CSV)
    if orc["cldn4_gene_id"].nunique() != 1 or orc["cldn4_gene_id"].iloc[0] != CLDN4_ENSG:
        raise ValueError("ORCESTRA extract is not a single CLDN4 ENSG row")
    if len(orc) != 41:
        raise ValueError(f"Expected 41 ORCESTRA patients; found {len(orc)}")
    orc = orc.set_index("sample_id")
    keep = [
        "recist",
        "response",
        "cldn4_orcestra_log2tpm",
        "survival_time_pfs",
        "event_occurred_pfs",
        "survival_time_os",
        "event_occurred_os",
        "sex",
        "age",
        "treatment",
        "Site.of.PRE.Biopsy",
    ]
    return orc[keep].rename(
        columns={
            "recist": "orcestra_recist",
            "response": "orcestra_response",
            "sex": "orcestra_sex",
            "age": "orcestra_age",
            "treatment": "orcestra_treatment",
        }
    )


def fisher_or(table: np.ndarray) -> dict[str, float | str]:
    """table rows = [low, high], cols = [nonresponder, responder]."""
    table = np.asarray(table, dtype=int)
    if table.shape != (2, 2):
        raise ValueError("Need a 2x2 table")
    fisher = stats.fisher_exact(table, alternative="two-sided")
    t2 = Table2x2(table)
    lo, hi = t2.oddsratio_confint()
    n_high_r = int(table[1, 1])
    n_high = int(table[1].sum())
    n_low_r = int(table[0, 1])
    n_low = int(table[0].sum())
    return {
        "fisher_or": float(fisher.statistic),
        "fisher_or_ci_low": float(lo),
        "fisher_or_ci_high": float(hi),
        "fisher_p": float(fisher.pvalue),
        "n": int(table.sum()),
        "n_high_resp": n_high_r,
        "n_high": n_high,
        "n_low_resp": n_low_r,
        "n_low": n_low,
        "orr_high": n_high_r / n_high if n_high else math.nan,
        "orr_low": n_low_r / n_low if n_low else math.nan,
        "high_over_n": f"{n_high_r}/{n_high}",
        "low_over_n": f"{n_low_r}/{n_low}",
        "table_low_nr_r_high_nr_r": json.dumps(table.tolist()),
    }


def median_split_or(df: pd.DataFrame, expr: str, y: str) -> dict[str, float | str]:
    work = df.dropna(subset=[expr, y]).copy()
    work[y] = work[y].astype(bool)
    median = float(work[expr].median())
    high = work[expr] >= median
    table = np.array(
        [
            [
                int(((~high) & (~work[y])).sum()),
                int(((~high) & work[y]).sum()),
            ],
            [
                int((high & (~work[y])).sum()),
                int((high & work[y]).sum()),
            ],
        ],
        dtype=int,
    )
    out = fisher_or(table)
    out["median"] = median
    out["n_high_group"] = int(high.sum())
    out["n_low_group"] = int((~high).sum())
    return out


def mannwhitney(df: pd.DataFrame, expr: str, y: str) -> dict[str, float | int]:
    work = df.dropna(subset=[expr, y])
    responders = work.loc[work[y].astype(bool), expr].to_numpy(dtype=float)
    nonresponders = work.loc[~work[y].astype(bool), expr].to_numpy(dtype=float)
    test = stats.mannwhitneyu(
        responders, nonresponders, alternative="two-sided", method="asymptotic"
    )
    auc = float(test.statistic / (len(responders) * len(nonresponders)))
    rng = np.random.default_rng(SEED)
    estimates = np.empty(20_000)
    for i in range(20_000):
        a = rng.choice(responders, len(responders), replace=True)
        b = rng.choice(nonresponders, len(nonresponders), replace=True)
        estimates[i] = stats.mannwhitneyu(a, b).statistic / (len(a) * len(b))
    lo, hi = np.quantile(estimates, [0.025, 0.975])
    return {
        "n_responder": int(len(responders)),
        "n_nonresponder": int(len(nonresponders)),
        "median_responder": float(np.median(responders)),
        "median_nonresponder": float(np.median(nonresponders)),
        "mann_whitney_u": float(test.statistic),
        "mann_whitney_p": float(test.pvalue),
        "auc_responder_higher": auc,
        "auc_bootstrap_95ci_low": float(lo),
        "auc_bootstrap_95ci_high": float(hi),
    }


def logistic_per_sd(df: pd.DataFrame, expr: str, y: str) -> dict[str, float | int]:
    work = df.dropna(subset=[expr, y]).copy()
    work[y] = work[y].astype(int)
    work["z"] = (work[expr] - work[expr].mean()) / work[expr].std(ddof=1)
    model = logit(f"{y} ~ z", work).fit(disp=0)
    coef = float(model.params["z"])
    ci = model.conf_int().loc["z"].to_numpy(dtype=float)
    return {
        "logit_or_per_sd": float(math.exp(coef)),
        "logit_or_ci_low": float(math.exp(ci[0])),
        "logit_or_ci_high": float(math.exp(ci[1])),
        "logit_p": float(model.pvalues["z"]),
        "n": int(len(work)),
    }


def cox_hr(df: pd.DataFrame, time: str, event: str, cov: str) -> dict[str, float | int]:
    work = df[[time, event, cov]].dropna().copy()
    work[cov] = work[cov].astype(float)
    cph = CoxPHFitter()
    cph.fit(work, duration_col=time, event_col=event)
    row = cph.summary.loc[cov]
    return {
        "n": int(len(work)),
        "n_event": int(work[event].sum()),
        "hr": float(row["exp(coef)"]),
        "hr_ci_low": float(row["exp(coef) lower 95%"]),
        "hr_ci_high": float(row["exp(coef) upper 95%"]),
        "cox_p": float(row["p"]),
    }


def detection_or(df: pd.DataFrame, tpm: str, y: str) -> dict[str, float | str]:
    work = df.dropna(subset=[tpm, y]).copy()
    detected = work[tpm] > 0
    yb = work[y].astype(bool)
    table = np.array(
        [
            [
                int(((~detected) & (~yb)).sum()),
                int(((~detected) & yb).sum()),
            ],
            [
                int((detected & (~yb)).sum()),
                int((detected & yb).sum()),
            ],
        ],
        dtype=int,
    )
    out = fisher_or(table)
    out["n_detected"] = int(detected.sum())
    out["n_undetected"] = int((~detected).sum())
    return out


def analyze_binary(df: pd.DataFrame, expr: str, y: str, tpm: str | None = None) -> dict:
    out = {}
    out.update({f"split_{k}": v for k, v in median_split_or(df, expr, y).items()})
    out.update({f"mw_{k}": v for k, v in mannwhitney(df, expr, y).items()})
    out.update({f"logit_{k}": v for k, v in logistic_per_sd(df, expr, y).items()})
    if tpm is not None:
        out.update({f"detect_{k}": v for k, v in detection_or(df, tpm, y).items()})
    return out


def write_tsv(path: Path, rows: list[dict], fieldnames: list[str] | None = None) -> None:
    frame = pd.DataFrame(rows)
    if fieldnames:
        frame = frame[fieldnames]
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False)


def make_figures(primary: pd.DataFrame, all_pre: pd.DataFrame, summaries: dict) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    plt.style.use("seaborn-v0_8-whitegrid")
    rng = np.random.default_rng(SEED)

    fig, axes = plt.subplots(1, 3, figsize=(12.2, 4.2))
    panels = [
        (
            "Anti-PD-1 PRE (primary)",
            primary,
            "cldn4_log2_tpm_plus_1",
            "orr_crpr_vs_sdpd",
            "log2(TPM + 1)",
            "primary_pre_pd1_orr",
        ),
        (
            "Anti-PD-1 PRE, CR/PR vs PD",
            primary[primary["recist"].isin(["CR", "PR", "PD"])],
            "cldn4_log2_tpm_plus_1",
            "orr_crpr_vs_pd",
            "log2(TPM + 1)",
            "pre_pd1_crpr_vs_pd",
        ),
        (
            "All PRE (PD-1 + combo)",
            all_pre,
            "cldn4_log2_tpm_plus_1",
            "orr_crpr_vs_sdpd",
            "log2(TPM + 1)",
            "all_pre_orr",
        ),
    ]
    colors = ["#6b7280", "#2563eb"]
    for ax, (title, frame, expr, y, ylabel, key) in zip(axes, panels, strict=True):
        yb = frame[y].astype(bool)
        groups = [
            frame.loc[~yb, expr].to_numpy(dtype=float),
            frame.loc[yb, expr].to_numpy(dtype=float),
        ]
        bp = ax.boxplot(
            groups,
            positions=[1, 2],
            widths=0.55,
            patch_artist=True,
            showfliers=False,
            medianprops={"color": "black", "linewidth": 1.5},
        )
        for box, color in zip(bp["boxes"], colors, strict=True):
            box.set_facecolor(color)
            box.set_alpha(0.35)
        for pos, values, color in zip([1, 2], groups, colors, strict=True):
            jitter = rng.uniform(-0.12, 0.12, len(values))
            ax.scatter(
                pos + jitter,
                values,
                s=23,
                color=color,
                alpha=0.8,
                edgecolor="white",
                linewidth=0.35,
                zorder=3,
            )
        ax.set_xticks([1, 2], [f"NR\nn={len(groups[0])}", f"R\nn={len(groups[1])}"])
        ax.set_title(title, fontsize=10, weight="bold")
        ax.set_ylabel(ylabel)
        p = summaries[key]["mw_mann_whitney_p"]
        ax.text(
            0.5,
            0.98,
            f"MW p={p:.3g}",
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=9,
        )
    fig.suptitle(
        "Gide 2019 melanoma: pretreatment CLDN4 versus ICI response", y=1.02
    )
    fig.tight_layout()
    fig.savefig(FIGURES / "cldn4_vs_response.png", dpi=220, bbox_inches="tight")
    fig.savefig(FIGURES / "cldn4_vs_response.svg", bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.2, 4.1))
    prim = summaries["primary_pre_pd1_orr"]
    rates = [prim["split_orr_low"], prim["split_orr_high"]]
    counts = [prim["split_low_over_n"], prim["split_high_over_n"]]
    bars = ax.bar(
        ["CLDN4-low\n(< median)", "CLDN4-high\n(≥ median)"],
        rates,
        color=["#6b7280", "#2563eb"],
        width=0.62,
        alpha=0.85,
    )
    for bar, rate, count in zip(bars, rates, counts, strict=True):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            rate + 0.02,
            f"{rate:.0%}\n({count})",
            ha="center",
            va="bottom",
            fontsize=9,
        )
    ax.set_ylim(0, 0.85)
    ax.set_ylabel("ORR (CR/PR)")
    ax.set_title(
        f"Primary median split  OR={prim['split_fisher_or']:.2f}  "
        f"n={prim['split_n']}  p={prim['split_fisher_p']:.2f}"
    )
    fig.tight_layout()
    fig.savefig(FIGURES / "orr_bar_median.png", dpi=220, bbox_inches="tight")
    fig.savefig(FIGURES / "orr_bar_median.svg", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)
    download_gide_rda()
    liulab = load_liulab()
    orcestra = load_orcestra()
    merged = liulab.join(orcestra, how="left")

    recist_mismatch = merged.dropna(subset=["orcestra_recist"])
    if not (recist_mismatch["recist"] == recist_mismatch["orcestra_recist"]).all():
        raise ValueError("RECIST disagrees between Liulab and ORCESTRA")

    merged["orr_crpr_vs_sdpd"] = merged["recist"].isin(["CR", "PR"])
    merged["orr_crpr_vs_pd"] = np.where(
        merged["recist"].isin(["CR", "PR"]),
        True,
        np.where(merged["recist"].eq("PD"), False, np.nan),
    )
    merged["paper_yes"] = merged["paper_response"].eq("Yes")
    merged["orcestra_r"] = np.where(
        merged["orcestra_response"].eq("R"),
        True,
        np.where(merged["orcestra_response"].eq("NR"), False, np.nan),
    )
    merged["cldn4_high_primary"] = pd.Series(pd.NA, index=merged.index, dtype="boolean")

    pre = merged[merged["visit"].eq("PRE")].copy()
    if pre["patient_id"].duplicated().any():
        raise ValueError("Duplicate pretreatment patients")
    primary = pre[pre["treatment_class"].eq("anti-PD-1")].copy()
    if len(primary) != 41:
        raise ValueError(f"Expected 41 pretreatment anti-PD-1 patients; got {len(primary)}")
    if primary["orcestra_recist"].notna().sum() != 41:
        raise ValueError("Primary set is not the complete ORCESTRA ICB_Gide object")
    primary_median = float(primary["cldn4_log2_tpm_plus_1"].median())
    primary["cldn4_high_primary"] = primary["cldn4_log2_tpm_plus_1"] >= primary_median
    merged.loc[primary.index, "cldn4_high_primary"] = primary["cldn4_high_primary"]

    combo = pre[pre["treatment_class"].eq("anti-PD-1 + anti-CTLA-4")].copy()
    edt = merged[merged["visit"].eq("EDT")].copy()
    orcestra_bin = primary.dropna(subset=["orcestra_r"]).copy()
    orcestra_bin["orcestra_r"] = orcestra_bin["orcestra_r"].astype(bool)
    strict = primary[primary["recist"].isin(["CR", "PR", "PD"])].copy()
    strict["orr_crpr_vs_pd"] = strict["orr_crpr_vs_pd"].astype(bool)

    analyses = {
        "primary_pre_pd1_orr": (
            primary,
            "cldn4_log2_tpm_plus_1",
            "orr_crpr_vs_sdpd",
            "cldn4_tpm",
        ),
        "pre_pd1_crpr_vs_pd": (
            strict,
            "cldn4_log2_tpm_plus_1",
            "orr_crpr_vs_pd",
            "cldn4_tpm",
        ),
        "pre_pd1_paper_yesno": (
            primary,
            "cldn4_log2_tpm_plus_1",
            "paper_yes",
            "cldn4_tpm",
        ),
        "pre_pd1_orcestra_rnr_liulab": (
            orcestra_bin,
            "cldn4_log2_tpm_plus_1",
            "orcestra_r",
            "cldn4_tpm",
        ),
        "pre_pd1_orcestra_rnr_orcestra": (
            orcestra_bin,
            "cldn4_orcestra_log2tpm",
            "orcestra_r",
            None,
        ),
        "all_pre_orr": (
            pre,
            "cldn4_log2_tpm_plus_1",
            "orr_crpr_vs_sdpd",
            "cldn4_tpm",
        ),
        "pre_combo_orr": (
            combo,
            "cldn4_log2_tpm_plus_1",
            "orr_crpr_vs_sdpd",
            "cldn4_tpm",
        ),
        "edt_orr": (
            edt,
            "cldn4_log2_tpm_plus_1",
            "orr_crpr_vs_sdpd",
            "cldn4_tpm",
        ),
    }
    summaries = {
        name: analyze_binary(frame, expr, y, tpm)
        for name, (frame, expr, y, tpm) in analyses.items()
    }

    primary = primary.copy()
    primary["pfs_months"] = primary["survival_time_pfs"].astype(float)
    primary["pfs_event"] = primary["event_occurred_pfs"].astype(int)
    primary["os_months"] = primary["survival_time_os"].astype(float)
    primary["os_event"] = primary["event_occurred_os"].astype(int)
    primary["high"] = primary["cldn4_high_primary"].astype(int)
    primary["z"] = (
        primary["cldn4_log2_tpm_plus_1"] - primary["cldn4_log2_tpm_plus_1"].mean()
    ) / primary["cldn4_log2_tpm_plus_1"].std(ddof=1)
    survival = {
        "pfs_median_split": cox_hr(primary, "pfs_months", "pfs_event", "high"),
        "os_median_split": cox_hr(primary, "os_months", "os_event", "high"),
        "pfs_per_sd": cox_hr(primary, "pfs_months", "pfs_event", "z"),
        "os_per_sd": cox_hr(primary, "os_months", "os_event", "z"),
    }
    hi = primary[primary["high"].eq(1)]
    lo = primary[primary["high"].eq(0)]
    survival["pfs_median_split"]["logrank_p"] = float(
        logrank_test(hi.pfs_months, lo.pfs_months, hi.pfs_event, lo.pfs_event).p_value
    )
    survival["os_median_split"]["logrank_p"] = float(
        logrank_test(hi.os_months, lo.os_months, hi.os_event, lo.os_event).p_value
    )

    spearman = stats.spearmanr(
        primary["cldn4_log2_tpm_plus_1"], primary["cldn4_orcestra_log2tpm"]
    )

    prim = summaries["primary_pre_pd1_orr"]
    headline = [
        {"item": "primary_cohort", "value": "Gide 2019 pretreatment anti-PD-1 RNA"},
        {"item": "primary_endpoint", "value": "ORR CR/PR vs SD/PD"},
        {"item": "primary_contrast", "value": "CLDN4-high vs low (median log2(TPM+1))"},
        {"item": "primary_ORR_fisher_OR_high_vs_low", "value": f"{prim['split_fisher_or']:.3f}"},
        {
            "item": "primary_ORR_fisher_OR_95CI",
            "value": (
                f"{prim['split_fisher_or_ci_low']:.3f}-"
                f"{prim['split_fisher_or_ci_high']:.3f}"
            ),
        },
        {"item": "primary_ORR_fisher_p", "value": f"{prim['split_fisher_p']:.4g}"},
        {"item": "primary_ORR_n_eval", "value": str(prim["split_n"])},
        {
            "item": "primary_ORR_n_high_resp_over_n_high",
            "value": prim["split_high_over_n"],
        },
        {
            "item": "primary_ORR_n_low_resp_over_n_low",
            "value": prim["split_low_over_n"],
        },
        {"item": "primary_ORR_rate_high", "value": f"{prim['split_orr_high']:.3f}"},
        {"item": "primary_ORR_rate_low", "value": f"{prim['split_orr_low']:.3f}"},
        {
            "item": "primary_MW_p",
            "value": f"{prim['mw_mann_whitney_p']:.4g}",
        },
        {
            "item": "primary_logit_OR_per_SD",
            "value": f"{prim['logit_logit_or_per_sd']:.3f}",
        },
        {
            "item": "primary_logit_p",
            "value": f"{prim['logit_logit_p']:.4g}",
        },
        {
            "item": "primary_PFS_HR_high_vs_low",
            "value": f"{survival['pfs_median_split']['hr']:.3f}",
        },
        {
            "item": "primary_PFS_HR_95CI",
            "value": (
                f"{survival['pfs_median_split']['hr_ci_low']:.3f}-"
                f"{survival['pfs_median_split']['hr_ci_high']:.3f}"
            ),
        },
        {
            "item": "primary_PFS_cox_p",
            "value": f"{survival['pfs_median_split']['cox_p']:.4g}",
        },
        {
            "item": "primary_OS_HR_high_vs_low",
            "value": f"{survival['os_median_split']['hr']:.3f}",
        },
        {
            "item": "primary_OS_cox_p",
            "value": f"{survival['os_median_split']['cox_p']:.4g}",
        },
        {"item": "CLDN4_median_log2_tpm_plus_1", "value": f"{primary_median:.6f}"},
        {"item": "n_primary_patients", "value": "41"},
        {"item": "n_CLDN4_zero_primary", "value": str(int((primary['cldn4_tpm'] == 0).sum()))},
        {
            "item": "spearman_liulab_vs_orcestra",
            "value": f"{float(spearman.statistic):.3f}",
        },
    ]
    write_tsv(TABLES / "headline_OR_n_p.tsv", headline, ["item", "value"])

    write_tsv(
        TABLES / "orr_2x2_primary.tsv",
        [
            {
                "CLDN4": "low",
                "SD_PD": int((~primary["cldn4_high_primary"] & ~primary["orr_crpr_vs_sdpd"]).sum()),
                "CR_PR": int((~primary["cldn4_high_primary"] & primary["orr_crpr_vs_sdpd"]).sum()),
            },
            {
                "CLDN4": "high",
                "SD_PD": int((primary["cldn4_high_primary"] & ~primary["orr_crpr_vs_sdpd"]).sum()),
                "CR_PR": int((primary["cldn4_high_primary"] & primary["orr_crpr_vs_sdpd"]).sum()),
            },
        ],
    )

    summary_rows = []
    for name, values in summaries.items():
        row = {"analysis": name}
        row.update(values)
        summary_rows.append(row)
    pd.DataFrame(summary_rows).to_csv(TABLES / "summary.csv", index=False)

    survival_rows = [{"analysis": name, **values} for name, values in survival.items()]
    pd.DataFrame(survival_rows).to_csv(TABLES / "survival.csv", index=False)

    sample_cols = [
        "sample_id",
        "patient_id",
        "visit",
        "drug",
        "treatment_class",
        "recist",
        "paper_response",
        "orr_crpr_vs_sdpd",
        "cldn4_tpm",
        "cldn4_log2_tpm_plus_1",
        "cldn4_high_primary",
        "orcestra_response",
        "cldn4_orcestra_log2tpm",
        "survival_time_pfs",
        "event_occurred_pfs",
        "survival_time_os",
        "event_occurred_os",
        "age",
        "sex",
        "Site.of.PRE.Biopsy",
        "Progression.Free.Survival..Days.",
        "Overall.Survival..Days.",
        "Progressed",
        "status",
    ]
    sample = merged.reset_index()
    sample.to_csv(TABLES / "sample_level.csv", index=False, columns=sample_cols)

    with (DATA / "provenance.tsv").open("w") as handle:
        handle.write("file\tsha256\turl\n")
        handle.write(f"{GIDE_RDA.name}\t{sha256(GIDE_RDA)}\t{GIDE_RDA_URL}\n")
        handle.write(
            f"{ORCESTRA_CSV.name}\t{sha256(ORCESTRA_CSV)}\t"
            f"{ORCESTRA_RDS_URL}#CLDN4_extract\n"
        )
        handle.write(
            f"ICB_Gide.rds\t{ORCESTRA_RDS_SHA256}\t{ORCESTRA_RDS_URL}\n"
        )

    make_figures(primary, pre, summaries)

    print("PRIMARY OR/n/p")
    print(
        f"OR={prim['split_fisher_or']:.3f} "
        f"[{prim['split_fisher_or_ci_low']:.3f}-{prim['split_fisher_or_ci_high']:.3f}] "
        f"n={prim['split_n']} p={prim['split_fisher_p']:.4g} "
        f"high={prim['split_high_over_n']} low={prim['split_low_over_n']}"
    )
    print(
        f"MW p={prim['mw_mann_whitney_p']:.4g} "
        f"logit OR/SD={prim['logit_logit_or_per_sd']:.3f} "
        f"p={prim['logit_logit_p']:.4g}"
    )
    print(
        "PFS HR",
        f"{survival['pfs_median_split']['hr']:.3f} p={survival['pfs_median_split']['cox_p']:.4g}",
    )
    print("Liulab vs ORCESTRA Spearman", float(spearman.statistic), "p", float(spearman.pvalue))


if __name__ == "__main__":
    main()
