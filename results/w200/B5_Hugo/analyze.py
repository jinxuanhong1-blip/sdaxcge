#!/usr/bin/env python3
"""Hugo 2016 melanoma anti-PD-1: TACSTD2 and CLDN4 versus response.

Primary analysis (prespecified):
  - Pretreatment biopsies only (exclude Pt16, the sole on-treatment sample).
  - One row per patient. Pt27A and Pt27B (same patient, both CR, both
    pretreatment) are collapsed by mean FPKM.
  - Responder = CR or PR; nonresponder = PD. This cohort has no SD.
  - Expression scale: log2(FPKM + 1).
  - Two-sided Mann-Whitney U and responder-higher AUC with a bootstrap CI.

Sensitivity analyses (also prespecified, not data-driven fishing):
  - Drop Pt28, the scalp subcutaneous biopsy with TACSTD2/CLDN4 ~100-1000x
    the rest of the cohort.
  - Keep both Pt27 biopsies (sample-level).
  - Add Pt16 back (on-treatment PD).
  - Three-level Kruskal-Wallis on CR / PR / PD.
  - Cohort-median split odds ratio of response (the B5 claim scale).
"""

from __future__ import annotations

import hashlib
import math
import warnings
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
FPKM_XLSX = DATA / "GSE78220_PatientFPKM.xlsx"
MATRIX = DATA / "GSE78220_series_matrix.txt.gz"
URLS = {
    FPKM_XLSX: (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE78nnn/GSE78220/suppl/"
        "GSE78220_PatientFPKM.xlsx"
    ),
    MATRIX: (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE78nnn/GSE78220/matrix/"
        "GSE78220_series_matrix.txt.gz"
    ),
}
GENES = ("TACSTD2", "CLDN4")
SEED = 78220
DETECT_FPKM = 1.0
# GEO last-update lines change; pin only the author-supplied FPKM workbook.
FPKM_SHA256 = "548c3d627df111c107d83b3f7d33acf0"  # md5 recorded 2026-08-16
# The line above is MD5. SHA-256 is computed at runtime and written to
# provenance.tsv. We verify MD5 because that is what was recorded on download.
FPKM_MD5 = "548c3d627df111c107d83b3f7d33acf0"


def download_inputs() -> None:
    import urllib.request

    DATA.mkdir(parents=True, exist_ok=True)
    for path, url in URLS.items():
        if not path.exists():
            print(f"Downloading {url}")
            urllib.request.urlretrieve(url, path)
    digest = hashlib.md5(FPKM_XLSX.read_bytes()).hexdigest()
    if digest != FPKM_MD5:
        raise ValueError(f"Unexpected MD5 for {FPKM_XLSX.name}: {digest}")


def parse_metadata() -> pd.DataFrame:
    import csv
    import gzip

    sample_rows: list[list[str]] = []
    with gzip.open(MATRIX, "rt") as handle:
        for line in handle:
            if line.startswith("!Sample_"):
                sample_rows.append(next(csv.reader([line], delimiter="\t")))

    def one_row(name: str) -> list[str]:
        matches = [row[1:] for row in sample_rows if row[0] == name]
        if len(matches) != 1:
            raise ValueError(f"Expected one {name} row; found {len(matches)}")
        return matches[0]

    titles = one_row("!Sample_title")
    accessions = one_row("!Sample_geo_accession")
    descriptions = one_row("!Sample_description")
    characteristics = [
        row[1:] for row in sample_rows if row[0] == "!Sample_characteristics_ch1"
    ]
    if len(titles) != 28:
        raise ValueError(f"Expected 28 GEO samples; found {len(titles)}")

    rows = []
    for i, title in enumerate(titles):
        kv: dict[str, str] = {}
        for values in characteristics:
            raw = values[i]
            if raw and ":" in raw:
                key, val = raw.split(":", 1)
                kv[key.strip()] = val.strip()
        response = kv["anti-pd-1 response"]
        if response == "Complete Response":
            recist = "CR"
        elif response == "Partial Response":
            recist = "PR"
        elif response == "Progressive Disease":
            recist = "PD"
        else:
            raise ValueError(f"Unexpected response {response!r}")
        rows.append(
            {
                "sample_title": title,
                "gsm": accessions[i],
                "description": descriptions[i],
                "patient": kv["patient id"],
                "response_full": response,
                "recist": recist,
                "binary_group": "Responder" if recist in {"CR", "PR"} else "Nonresponder",
                "biopsy_time": kv.get("biopsy time", ""),
                "tissue": kv.get("tissue", ""),
                "study_site": kv.get("study site", ""),
                "gender": kv.get("gender", ""),
                "age_yrs": kv.get("age (yrs)", ""),
                "disease_status": kv.get("disease status", ""),
                "os_days": kv.get("overall survival (days)", ""),
                "vital_status": kv.get("vital status", ""),
                "previous_mapki": kv.get("previous mapki", ""),
                "anatomical_location": kv.get("anatomical location", ""),
            }
        )
    meta = pd.DataFrame(rows)
    if (meta["tissue"] != "Melanoma biopsies").any():
        raise ValueError("Non-melanoma tissue present; this analysis is melanoma-only")
    return meta


def parse_expression() -> pd.DataFrame:
    warnings.filterwarnings(
        "ignore",
        message="Unknown extension is not supported and will be removed",
    )
    raw = pd.read_excel(FPKM_XLSX, sheet_name="FPKM")
    if raw.shape != (25268, 29):
        raise ValueError(f"Unexpected FPKM shape {raw.shape}")
    if raw["Gene"].duplicated().any():
        raise ValueError("Duplicate gene symbols in FPKM matrix")
    missing = [g for g in GENES if g not in set(raw["Gene"])]
    if missing:
        raise ValueError(f"Missing genes: {missing}")
    expr = raw.set_index("Gene").loc[list(GENES)].T
    expr.index.name = "fpkm_column"
    expr = expr.astype(float)
    return expr


def column_to_title(column: str) -> str:
    if column.endswith(".baseline"):
        return column[: -len(".baseline")]
    if column.endswith(".OnTx"):
        return column[: -len(".OnTx")]
    raise ValueError(f"Unrecognized FPKM column {column!r}")


def log2p1(x: float) -> float:
    return math.log2(x + 1.0)


def bootstrap_auc_ci(
    responders: np.ndarray, nonresponders: np.ndarray, draws: int = 20_000
) -> tuple[float, float]:
    rng = np.random.default_rng(SEED)
    estimates = np.empty(draws)
    n_r, n_nr = len(responders), len(nonresponders)
    for i in range(draws):
        a = rng.choice(responders, n_r, replace=True)
        b = rng.choice(nonresponders, n_nr, replace=True)
        estimates[i] = stats.mannwhitneyu(a, b, alternative="two-sided").statistic / (
            n_r * n_nr
        )
    return tuple(np.quantile(estimates, [0.025, 0.975]))


def odds_ratio_ci(a: int, b: int, c: int, d: int) -> tuple[float, float, float]:
    """Response OR for high vs low, Haldane-Anscombe 0.5 correction, Wald CI."""
    aa, bb, cc, dd = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    odds = (aa * dd) / (bb * cc)
    se = math.sqrt(1 / aa + 1 / bb + 1 / cc + 1 / dd)
    lo = math.exp(math.log(odds) - 1.96 * se)
    hi = math.exp(math.log(odds) + 1.96 * se)
    return odds, lo, hi


def summarize_continuous(values_r: np.ndarray, values_nr: np.ndarray) -> dict:
    test = stats.mannwhitneyu(
        values_r, values_nr, alternative="two-sided", method="asymptotic"
    )
    n_r, n_nr = len(values_r), len(values_nr)
    auc = float(test.statistic / (n_r * n_nr))
    auc_lo, auc_hi = bootstrap_auc_ci(values_r, values_nr)
    return {
        "n_responder": n_r,
        "n_nonresponder": n_nr,
        "median_responder": float(np.median(values_r)),
        "median_nonresponder": float(np.median(values_nr)),
        "mean_responder": float(np.mean(values_r)),
        "mean_nonresponder": float(np.mean(values_nr)),
        "mann_whitney_u": float(test.statistic),
        "mann_whitney_p": float(test.pvalue),
        "auc_responder_higher": auc,
        "auc_bootstrap_95ci_low": float(auc_lo),
        "auc_bootstrap_95ci_high": float(auc_hi),
        "cliffs_delta": float(2 * auc - 1),
    }


def summarize_detection(fpkm_r: np.ndarray, fpkm_nr: np.ndarray) -> dict:
    det_r = int(np.sum(fpkm_r > DETECT_FPKM))
    det_nr = int(np.sum(fpkm_nr > DETECT_FPKM))
    table = [
        [det_r, len(fpkm_r) - det_r],
        [det_nr, len(fpkm_nr) - det_nr],
    ]
    fisher = stats.fisher_exact(table, alternative="two-sided")
    return {
        "detected_responder": det_r,
        "detected_nonresponder": det_nr,
        "detection_threshold_fpkm": DETECT_FPKM,
        "detection_fisher_odds_ratio": float(fisher.statistic),
        "detection_fisher_p": float(fisher.pvalue),
    }


def summarize_median_split(fpkm_r: np.ndarray, fpkm_nr: np.ndarray) -> dict:
    all_fpkm = np.concatenate([fpkm_r, fpkm_nr])
    cut = float(np.median(all_fpkm))
    high_r = int(np.sum(fpkm_r > cut))
    low_r = int(np.sum(fpkm_r <= cut))
    high_nr = int(np.sum(fpkm_nr > cut))
    low_nr = int(np.sum(fpkm_nr <= cut))
    fisher = stats.fisher_exact(
        [[high_r, high_nr], [low_r, low_nr]], alternative="two-sided"
    )
    or_hat, or_lo, or_hi = odds_ratio_ci(high_r, high_nr, low_r, low_nr)
    return {
        "median_split_cut_fpkm": cut,
        "high_responder": high_r,
        "high_nonresponder": high_nr,
        "low_responder": low_r,
        "low_nonresponder": low_nr,
        "response_or_high_vs_low": or_hat,
        "response_or_95ci_low": or_lo,
        "response_or_95ci_high": or_hi,
        "median_split_fisher_p": float(fisher.pvalue),
    }


def group_arrays(frame: pd.DataFrame, gene: str, value_col: str) -> tuple[np.ndarray, np.ndarray]:
    r = frame.loc[frame["binary_group"] == "Responder", value_col].to_numpy(float)
    nr = frame.loc[frame["binary_group"] == "Nonresponder", value_col].to_numpy(float)
    return r, nr


def collapse_patients(sample_frame: pd.DataFrame) -> pd.DataFrame:
    """Mean FPKM across biopsies from the same patient."""
    num = [f"{g}_fpkm" for g in GENES]
    first = sample_frame.drop(columns=num).groupby("patient", as_index=False).first()
    means = sample_frame.groupby("patient", as_index=False)[num].mean()
    out = first.merge(means, on="patient")
    for gene in GENES:
        out[f"{gene}_log2_fpkm_plus_1"] = out[f"{gene}_fpkm"].map(log2p1)
    out["n_biopsies"] = (
        sample_frame.groupby("patient").size().reindex(out["patient"]).to_numpy()
    )
    return out


def analysis_bundle(frame: pd.DataFrame, analysis: str, gene: str) -> dict:
    log_col = f"{gene}_log2_fpkm_plus_1"
    fpkm_col = f"{gene}_fpkm"
    r_log, nr_log = group_arrays(frame, gene, log_col)
    r_fpkm, nr_fpkm = group_arrays(frame, gene, fpkm_col)
    recist = frame["recist"].to_numpy()
    log_vals = frame[log_col].to_numpy(float)
    groups = [log_vals[recist == lab] for lab in ("CR", "PR", "PD")]
    nonempty = [g for g in groups if len(g)]
    if len(nonempty) >= 2 and all(len(g) >= 1 for g in nonempty):
        kw = stats.kruskal(*nonempty)
        kw_stat, kw_p = float(kw.statistic), float(kw.pvalue)
    else:
        kw_stat, kw_p = float("nan"), float("nan")
    row = {
        "analysis": analysis,
        "gene": gene,
        "n_patients_or_samples": len(frame),
        **summarize_continuous(r_log, nr_log),
        **summarize_detection(r_fpkm, nr_fpkm),
        **summarize_median_split(r_fpkm, nr_fpkm),
        "kruskal_cr_pr_pd_h": kw_stat,
        "kruskal_cr_pr_pd_p": kw_p,
    }
    return row


def write_csv(path: Path, rows: list[dict]) -> None:
    pd.DataFrame(rows).to_csv(path, index=False)


def jitter(rng: np.random.Generator, n: int, width: float = 0.12) -> np.ndarray:
    return rng.uniform(-width, width, n)


def box_panel(ax, frame: pd.DataFrame, gene: str, title: str, pvalue: float) -> None:
    rng = np.random.default_rng(SEED)
    log_col = f"{gene}_log2_fpkm_plus_1"
    ordered = [
        frame.loc[frame["binary_group"] == "Nonresponder", log_col].to_numpy(float),
        frame.loc[frame["binary_group"] == "Responder", log_col].to_numpy(float),
    ]
    colors = ["#6b7280", "#2563eb"]
    bp = ax.boxplot(
        ordered,
        positions=[1, 2],
        widths=0.55,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "black", "linewidth": 1.5},
    )
    for box, color in zip(bp["boxes"], colors, strict=True):
        box.set_facecolor(color)
        box.set_alpha(0.35)
    for pos, values, color in zip([1, 2], ordered, colors, strict=True):
        ax.scatter(
            pos + jitter(rng, len(values)),
            values,
            s=28,
            color=color,
            alpha=0.85,
            edgecolor="white",
            linewidth=0.4,
            zorder=3,
        )
    # Label the epithelial-like scalp biopsy if present.
    hit = frame.loc[frame["patient"] == "Pt28"]
    if len(hit):
        y = float(hit.iloc[0][log_col])
        x = 2 if hit.iloc[0]["binary_group"] == "Responder" else 1
        ax.annotate(
            "Pt28",
            (x, y),
            textcoords="offset points",
            xytext=(8, 4),
            fontsize=8,
            color="#b45309",
        )
    ax.set_xticks([1, 2], [f"NR\nn={len(ordered[0])}", f"R\nn={len(ordered[1])}"])
    ax.set_title(title, fontsize=10, weight="bold")
    ax.set_ylabel(f"{gene}  log2(FPKM + 1)")
    ax.text(
        0.5,
        0.98,
        f"MW p={pvalue:.3g}",
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=9,
    )


def make_figures(
    primary: pd.DataFrame,
    no_pt28: pd.DataFrame,
    summaries: pd.DataFrame,
) -> None:
    plt.style.use("seaborn-v0_8-whitegrid")
    pmap = {
        (row.analysis, row.gene): row.mann_whitney_p
        for row in summaries.itertuples(index=False)
    }

    fig, axes = plt.subplots(2, 2, figsize=(8.6, 7.2))
    box_panel(
        axes[0, 0],
        primary,
        "TACSTD2",
        "TACSTD2 pretreatment (primary)",
        pmap[("primary_patient_pretreatment", "TACSTD2")],
    )
    box_panel(
        axes[0, 1],
        no_pt28,
        "TACSTD2",
        "TACSTD2, Pt28 dropped",
        pmap[("sensitivity_drop_Pt28", "TACSTD2")],
    )
    box_panel(
        axes[1, 0],
        primary,
        "CLDN4",
        "CLDN4 pretreatment (primary)",
        pmap[("primary_patient_pretreatment", "CLDN4")],
    )
    box_panel(
        axes[1, 1],
        no_pt28,
        "CLDN4",
        "CLDN4, Pt28 dropped",
        pmap[("sensitivity_drop_Pt28", "CLDN4")],
    )
    fig.suptitle(
        "Hugo 2016 melanoma: TACSTD2 and CLDN4 versus anti-PD-1 response",
        y=1.01,
    )
    fig.tight_layout()
    fig.savefig(ROOT / "tacstd2_cldn4_vs_response.png", dpi=220, bbox_inches="tight")
    fig.savefig(ROOT / "tacstd2_cldn4_vs_response.svg", bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.6, 4.6))
    for group, color, label in (
        ("Nonresponder", "#6b7280", "PD"),
        ("Responder", "#2563eb", "CR/PR"),
    ):
        sub = primary.loc[primary["binary_group"] == group]
        ax.scatter(
            sub["TACSTD2_fpkm"],
            sub["CLDN4_fpkm"],
            s=36,
            color=color,
            alpha=0.85,
            edgecolor="white",
            linewidth=0.4,
            label=f"{label} n={len(sub)}",
        )
    for _, row in primary.iterrows():
        if row["TACSTD2_fpkm"] >= 5 or row["CLDN4_fpkm"] >= 5:
            ax.annotate(
                row["patient"],
                (row["TACSTD2_fpkm"], row["CLDN4_fpkm"]),
                textcoords="offset points",
                xytext=(6, 4),
                fontsize=8,
            )
    ax.set_xscale("symlog", linthresh=1)
    ax.set_yscale("symlog", linthresh=1)
    ax.set_xlabel("TACSTD2 FPKM")
    ax.set_ylabel("CLDN4 FPKM")
    ax.set_title("Almost all melanomas sit near zero; Pt28 does not")
    ax.legend(frameon=True, fontsize=8)
    fig.tight_layout()
    fig.savefig(ROOT / "fpkm_scatter.png", dpi=220, bbox_inches="tight")
    fig.savefig(ROOT / "fpkm_scatter.svg", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    download_inputs()
    meta = parse_metadata()
    expr = parse_expression()
    expr = expr.reset_index()
    expr["sample_title"] = expr["fpkm_column"].map(column_to_title)
    if set(expr["sample_title"]) != set(meta["sample_title"]):
        raise ValueError("FPKM columns and GEO sample titles do not match")

    samples = meta.merge(expr, on="sample_title", how="inner")
    if len(samples) != 28:
        raise ValueError(f"Expected 28 joined samples; got {len(samples)}")
    for gene in GENES:
        samples[f"{gene}_fpkm"] = samples[gene].astype(float)
        samples[f"{gene}_log2_fpkm_plus_1"] = samples[f"{gene}_fpkm"].map(log2p1)
    samples = samples.drop(columns=list(GENES))
    samples["is_pretreatment"] = samples["biopsy_time"] == "pre-treatment"
    samples["is_pt28"] = samples["patient"] == "Pt28"

    pretreat = samples.loc[samples["is_pretreatment"]].copy()
    primary = collapse_patients(pretreat)
    if len(primary) != 26:
        raise ValueError(f"Expected 26 pretreatment patients; got {len(primary)}")
    no_pt28 = primary.loc[primary["patient"] != "Pt28"].copy()
    sample_level_pre = pretreat.copy()
    all_patients = collapse_patients(samples)

    summary_rows = []
    for gene in GENES:
        summary_rows.append(analysis_bundle(primary, "primary_patient_pretreatment", gene))
        summary_rows.append(analysis_bundle(no_pt28, "sensitivity_drop_Pt28", gene))
        summary_rows.append(
            analysis_bundle(sample_level_pre, "sensitivity_sample_level_pretreatment", gene)
        )
        summary_rows.append(
            analysis_bundle(all_patients, "sensitivity_include_on_treatment_Pt16", gene)
        )

    summaries = pd.DataFrame(summary_rows)

    rho, rho_p = stats.spearmanr(primary["TACSTD2_fpkm"], primary["CLDN4_fpkm"])
    rho_no, rho_no_p = stats.spearmanr(no_pt28["TACSTD2_fpkm"], no_pt28["CLDN4_fpkm"])
    corr_rows = [
        {
            "analysis": "primary_patient_pretreatment",
            "spearman_tacstd2_cldn4": float(rho),
            "spearman_p": float(rho_p),
            "n": len(primary),
        },
        {
            "analysis": "sensitivity_drop_Pt28",
            "spearman_tacstd2_cldn4": float(rho_no),
            "spearman_p": float(rho_no_p),
            "n": len(no_pt28),
        },
    ]

    write_csv(ROOT / "summary.csv", summary_rows)
    write_csv(ROOT / "gene_correlation.csv", corr_rows)
    sample_cols = [
        "patient",
        "sample_title",
        "fpkm_column",
        "gsm",
        "recist",
        "binary_group",
        "biopsy_time",
        "anatomical_location",
        "previous_mapki",
        "study_site",
        "gender",
        "age_yrs",
        "disease_status",
        "os_days",
        "vital_status",
        "TACSTD2_fpkm",
        "CLDN4_fpkm",
        "TACSTD2_log2_fpkm_plus_1",
        "CLDN4_log2_fpkm_plus_1",
        "description",
    ]
    samples.loc[:, sample_cols].sort_values(["patient", "sample_title"]).to_csv(
        ROOT / "sample_level.csv", index=False
    )
    patient_cols = [
        "patient",
        "n_biopsies",
        "recist",
        "binary_group",
        "anatomical_location",
        "previous_mapki",
        "study_site",
        "TACSTD2_fpkm",
        "CLDN4_fpkm",
        "TACSTD2_log2_fpkm_plus_1",
        "CLDN4_log2_fpkm_plus_1",
    ]
    primary.loc[:, patient_cols].sort_values("patient").to_csv(
        ROOT / "patient_level_primary.csv", index=False
    )
    make_figures(primary, no_pt28, summaries)

    with (ROOT / "provenance.tsv").open("w") as handle:
        handle.write("file\tmd5\tsha256\turl\n")
        for path, url in URLS.items():
            raw = path.read_bytes()
            handle.write(
                f"{path.name}\t{hashlib.md5(raw).hexdigest()}\t"
                f"{hashlib.sha256(raw).hexdigest()}\t{url}\n"
            )

    print("GEO RECIST:", dict(Counter(samples["recist"])))
    print("Pretreatment patients:", len(primary), dict(Counter(primary["recist"])))
    print(
        "Pt28 FPKM TACSTD2/CLDN4:",
        float(primary.loc[primary["patient"] == "Pt28", "TACSTD2_fpkm"].iloc[0]),
        float(primary.loc[primary["patient"] == "Pt28", "CLDN4_fpkm"].iloc[0]),
    )
    for row in summary_rows:
        if row["analysis"] in {
            "primary_patient_pretreatment",
            "sensitivity_drop_Pt28",
        }:
            print(
                f"{row['analysis']} {row['gene']}: "
                f"R={row['n_responder']} NR={row['n_nonresponder']} "
                f"medR={row['median_responder']:.3f} medNR={row['median_nonresponder']:.3f} "
                f"AUC={row['auc_responder_higher']:.3f} "
                f"MW p={row['mann_whitney_p']:.4g} "
                f"OR={row['response_or_high_vs_low']:.3f} "
                f"Fisher p={row['median_split_fisher_p']:.4g}"
            )
    print(f"Spearman TACSTD2-CLDN4 primary rho={rho:.3f} p={rho_p:.4g}")
    print(f"Spearman TACSTD2-CLDN4 no Pt28 rho={rho_no:.3f} p={rho_no_p:.4g}")


if __name__ == "__main__":
    main()
