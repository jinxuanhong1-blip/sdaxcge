#!/usr/bin/env python3
"""B5 analog: public GBM/glioma ICI RNA, CLDN4 versus clinical outcome.

Primary analysis uses GSE121810 (Cloughesy 2019) because it is the only
located public bulk tumor RNA-seq matrix that is both strand-consistent
and joinable to a public per-patient survival table.

GSE264695 and GSE318645 are retained as QC-failed sensitivity analyses.
No public RECIST/ORR labels were found for these RNA cohorts.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import re
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
from scipy.stats import mannwhitneyu


H_STRAND_MT = [
    "MT-CO1",
    "MT-CO2",
    "MT-CO3",
    "MT-CYB",
    "MT-ND1",
    "MT-ND2",
    "MT-ND3",
    "MT-ND4",
    "MT-ND5",
    "MT-ATP6",
]

SOURCES = {
    "GSE121810_counts": {
        "url": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE121nnn/GSE121810/"
            "suppl/GSE121810_Prins.PD1NeoAdjv.Jul2018.HUGO.PtID.xlsx"
        ),
        "sha256": "078aeda14e23168400f6d50c3fa92dd611072c50da91613f34cf3e75b0a0c288",
        "filename": "GSE121810_Prins.PD1NeoAdjv.Jul2018.HUGO.PtID.xlsx",
    },
    "GSE264695_counts": {
        "url": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE264nnn/GSE264695/suppl/"
            "GSE264695_Prins.PD1NeoAdjv.Combo.Merck.Nov2022.htseqCount.SampleRenamedForGEO.xlsx"
        ),
        "sha256": "0adde00580f394c23a0ae64f3da42168c5e0cc132bd22027e92507545a8fbaf7",
        "filename": "GSE264695_Prins.PD1NeoAdjv.Combo.Merck.Nov2022.htseqCount.SampleRenamedForGEO.xlsx",
    },
    "GSE318645_counts": {
        "url": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE318nnn/GSE318645/"
            "suppl/GSE318645_IpiNivo.Tumor.htseqCount.RenamedForGEO.txt.gz"
        ),
        "sha256": "bfa2ac335f93579176e394fd2e90f7f8fa0ad37087cd64306696f40653587508",
        "filename": "GSE318645_IpiNivo.Tumor.htseqCount.RenamedForGEO.txt.gz",
    },
    "GSE318645_matrix": {
        "url": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE318nnn/GSE318645/"
            "matrix/GSE318645_series_matrix.txt.gz"
        ),
        "sha256": "ddc80257127e2533e15394c3e10af78b8211933d5bafb8ca3791628dad8fb33e",
        "filename": "GSE318645_series_matrix.txt.gz",
    },
    "NC2024_clinical": {
        "url": (
            "https://static-content.springer.com/esm/"
            "art%3A10.1038%2Fs41467-024-54326-7/MediaObjects/"
            "41467_2024_54326_MOESM4_ESM.xlsx"
        ),
        "sha256": "2b216cbba0c280ab8c40dc36153703262ddbd59c41e5f9e143866d2da60d4311",
        "filename": "41467_2024_54326_MOESM4_ESM.xlsx",
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def obtain(source: dict, cache_dir: Path) -> Path:
    path = cache_dir / source["filename"]
    if path.exists() and sha256(path) == source["sha256"]:
        return path
    cache_dir.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".download")
    print(f"Downloading {source['url']}")
    urllib.request.urlretrieve(source["url"], temporary)
    observed = sha256(temporary)
    if observed != source["sha256"]:
        temporary.unlink(missing_ok=True)
        raise RuntimeError(
            f"Checksum changed for {source['filename']}: {observed}"
        )
    temporary.replace(path)
    return path


def load_count_matrix(path: Path) -> pd.DataFrame:
    if path.suffix == ".gz":
        frame = pd.read_csv(path, sep="\t", index_col=0)
    else:
        frame = pd.read_excel(path, index_col=0)
    frame.index = frame.index.astype(str)
    frame = frame[~frame.index.str.startswith("__")]
    frame = frame.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    return frame.groupby(level=0).sum()


def library_metrics(counts: pd.DataFrame) -> dict:
    assigned = counts.sum(axis=0)
    antisense = counts.index.str.endswith(("-AS1", "-AS2"))
    mt_h = counts.reindex([gene for gene in H_STRAND_MT if gene in counts.index]).sum(
        axis=0
    )
    nd6 = counts.loc["MT-ND6"] if "MT-ND6" in counts.index else pd.Series(np.nan, index=counts.columns)
    nd6_ratio = float((nd6 / mt_h.replace(0, np.nan)).median())
    as_frac = float((counts.loc[antisense].sum(axis=0) / assigned.replace(0, np.nan)).median())
    # Wrong-strand HTSeq matrices put most biological signal into antisense
    # and L-strand MT-ND6. GSE121810 is the only matrix that looks stranded.
    strand_ok = bool(nd6_ratio < 0.05 and as_frac < 0.02)
    return {
        "n_genes": int(counts.shape[0]),
        "n_samples": int(counts.shape[1]),
        "assigned_reads_median": float(assigned.median()),
        "assigned_reads_min": float(assigned.min()),
        "assigned_reads_max": float(assigned.max()),
        "genes_detected_median": float((counts > 0).sum(axis=0).median()),
        "mt_nd6_over_h_strand_mt_median": nd6_ratio,
        "antisense_fraction_median": as_frac,
        "strand_ok": strand_ok,
    }


def cpm_table(counts: pd.DataFrame) -> pd.DataFrame:
    lib = counts.sum(axis=0)
    return counts.div(lib, axis=1) * 1e6


def gene_detectability(counts: pd.DataFrame, gene: str) -> dict:
    if gene not in counts.index:
        return {"present_in_matrix": False}
    raw = counts.loc[gene]
    cpm = cpm_table(counts).loc[gene]
    return {
        "present_in_matrix": True,
        "nonzero_n": int((raw > 0).sum()),
        "n": int(raw.size),
        "raw_median": float(raw.median()),
        "raw_max": float(raw.max()),
        "cpm_median": float(cpm.median()),
        "cpm_q1": float(cpm.quantile(0.25)),
        "cpm_q3": float(cpm.quantile(0.75)),
        "cpm_max": float(cpm.max()),
        "n_cpm_ge_1": int((cpm >= 1).sum()),
        "n_cpm_ge_5": int((cpm >= 5).sum()),
    }


def quantiles(values: np.ndarray) -> dict:
    values = np.asarray(values, dtype=float)
    q1, median, q3 = np.quantile(values, [0.25, 0.5, 0.75])
    return {
        "n": int(values.size),
        "median": float(median),
        "q1": float(q1),
        "q3": float(q3),
    }


def bootstrap_median_difference(
    group_a: np.ndarray, group_b: np.ndarray, seed: int = 20260816
) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    differences = np.empty(20_000)
    for index in range(differences.size):
        differences[index] = np.median(
            rng.choice(group_a, group_a.size, replace=True)
        ) - np.median(rng.choice(group_b, group_b.size, replace=True))
    lower, upper = np.quantile(differences, [0.025, 0.975])
    return float(lower), float(upper)


def cox_continuous(time: np.ndarray, event: np.ndarray, x: np.ndarray) -> dict:
    frame = pd.DataFrame({"os_days": time, "os_event": event, "cldn4_log2cpm": x})
    model = CoxPHFitter()
    model.fit(frame, duration_col="os_days", event_col="os_event")
    summary = model.summary.loc["cldn4_log2cpm"]
    return {
        "n": int(len(frame)),
        "events": int(frame["os_event"].sum()),
        "hr_per_log2cpm": float(summary["exp(coef)"]),
        "hr_95_ci": [float(summary["exp(coef) lower 95%"]), float(summary["exp(coef) upper 95%"])],
        "p_value": float(summary["p"]),
        "concordance": float(model.concordance_index_),
    }


def cox_binary(time: np.ndarray, event: np.ndarray, high: np.ndarray) -> dict:
    frame = pd.DataFrame({"os_days": time, "os_event": event, "cldn4_high": high.astype(int)})
    model = CoxPHFitter()
    model.fit(frame, duration_col="os_days", event_col="os_event")
    summary = model.summary.loc["cldn4_high"]
    return {
        "n": int(len(frame)),
        "n_high": int(high.sum()),
        "n_low": int((~high).sum()),
        "events": int(frame["os_event"].sum()),
        "hr_high_vs_low": float(summary["exp(coef)"]),
        "hr_95_ci": [float(summary["exp(coef) lower 95%"]), float(summary["exp(coef) upper 95%"])],
        "p_value": float(summary["p"]),
        "concordance": float(model.concordance_index_),
    }


def binary_expression_test(high_outcome: np.ndarray, low_outcome: np.ndarray) -> dict:
    test = mannwhitneyu(high_outcome, low_outcome, alternative="two-sided", method="asymptotic")
    median_difference = float(np.median(high_outcome) - np.median(low_outcome))
    ci_lower, ci_upper = bootstrap_median_difference(high_outcome, low_outcome)
    return {
        "group_a": quantiles(high_outcome),
        "group_b": quantiles(low_outcome),
        "median_difference_a_minus_b": median_difference,
        "median_difference_bootstrap_95_ci": [ci_lower, ci_upper],
        "mann_whitney_u": float(test.statistic),
        "p_value_two_sided": float(test.pvalue),
        "common_language_effect_probability": float(
            test.statistic / (high_outcome.size * low_outcome.size)
        ),
    }


def spearman_os(records: list[dict], time_key: str) -> dict:
    from scipy.stats import spearmanr

    usable = [
        record
        for record in records
        if record.get(time_key) is not None
        and not (isinstance(record[time_key], float) and math.isnan(record[time_key]))
    ]
    rho, p_value = spearmanr(
        [record["cldn4_log2cpm"] for record in usable],
        [record[time_key] for record in usable],
    )
    return {"n": len(usable), "spearman_rho": float(rho), "p_value": float(p_value)}


def analyze_survival(
    records: list[dict], time_key: str, event_key: str, landmark_days: int | None = None
) -> dict:
    usable = [
        record
        for record in records
        if record.get(time_key) is not None
        and record.get(event_key) is not None
        and not (isinstance(record[time_key], float) and math.isnan(record[time_key]))
    ]
    time = np.array([record[time_key] for record in usable], dtype=float)
    event = np.array([record[event_key] for record in usable], dtype=int)
    x = np.array([record["cldn4_log2cpm"] for record in usable], dtype=float)
    high = x >= np.median(x)
    logrank = logrank_test(time[high], time[~high], event[high], event[~high])
    result = {
        "n_with_endpoint": len(usable),
        "continuous_cox": cox_continuous(time, event, x),
        "median_split_cox": cox_binary(time, event, high),
        "median_split_logrank_p": float(logrank.p_value),
        "spearman_time": spearman_os(usable, time_key),
    }
    if landmark_days is not None:
        longer = x[time >= landmark_days]
        shorter = x[time < landmark_days]
        result["os365_expression_test"] = {
            "definition": (
                "Exploratory binary outcome OS>=365 days versus OS<365 days. "
                "This is not RECIST/ORR."
            ),
            "longer_survivor": quantiles(longer) if longer.size else None,
            "shorter_survivor": quantiles(shorter) if shorter.size else None,
            **(
                {
                    "test": binary_expression_test(longer, shorter),
                    "group_a_label": "OS>=365d",
                    "group_b_label": "OS<365d",
                }
                if longer.size and shorter.size
                else {"test": None}
            ),
        }
    return result


def cox_arm_adjusted(records: list[dict]) -> dict:
    frame = pd.DataFrame(
        {
            "os_days": [record["os_days"] for record in records],
            "os_event": [record["os_event"] for record in records],
            "cldn4_log2cpm": [record["cldn4_log2cpm"] for record in records],
            "neoadjuvant": [int(record["arm"] == "neoadjuvant") for record in records],
        }
    )
    model = CoxPHFitter()
    model.fit(frame, duration_col="os_days", event_col="os_event")
    cldn4 = model.summary.loc["cldn4_log2cpm"]
    arm = model.summary.loc["neoadjuvant"]
    return {
        "n": int(len(frame)),
        "events": int(frame["os_event"].sum()),
        "cldn4_hr_per_log2cpm": float(cldn4["exp(coef)"]),
        "cldn4_hr_95_ci": [
            float(cldn4["exp(coef) lower 95%"]),
            float(cldn4["exp(coef) upper 95%"]),
        ],
        "cldn4_p_value": float(cldn4["p"]),
        "neoadjuvant_hr": float(arm["exp(coef)"]),
        "neoadjuvant_p_value": float(arm["p"]),
        "concordance": float(model.concordance_index_),
    }


def load_nc2024_clinical(path: Path) -> pd.DataFrame:
    clinical = pd.read_excel(path, sheet_name="SFig4C", header=2)
    clinical = clinical.dropna(subset=["SampleID"]).copy()
    clinical["SampleID"] = clinical["SampleID"].astype(str)
    rename = {
        "Treatment Group (B: Adjuvant; A: Neoadjuvant)": "treatment_group",
    }
    clinical = clinical.rename(columns=rename)
    return clinical


def load_gse318645_clinical(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        text = handle.read().split("!series_matrix_table_begin")[0]
    rows: dict[str, list[str]] = {}
    for line in text.splitlines():
        if not line.startswith("!Sample_"):
            continue
        tag = line.split("\t", 1)[0]
        values = [value.strip('"') for value in line.split("\t")[1:]]
        if tag == "!Sample_title":
            rows["title"] = values
        elif tag == "!Sample_geo_accession":
            rows["gsm"] = values
        elif tag == "!Sample_source_name_ch1":
            rows["source"] = values
        elif tag.startswith("!Sample_characteristics"):
            keys = {value.split(":", 1)[0] for value in values if ":" in value}
            if len(keys) == 1:
                key = next(iter(keys))
                rows[key] = [
                    value.split(":", 1)[1].strip() if ":" in value else value
                    for value in values
                ]
    frame = pd.DataFrame(rows)
    tumor = frame[frame["source"] == "Tumor"].copy()
    tumor["sample_id"] = tumor["title"].str.extract(r"Patient (\d+)", expand=False)
    tumor["sample_id"] = tumor["sample_id"].map(lambda value: f"P{int(value):03d}")
    tumor["os_days"] = tumor["os"].astype(int)
    tumor["os_event"] = tumor["os.event"].astype(int)
    tumor["age"] = tumor["age"].astype(float)
    tumor["arm"] = tumor["arm"].astype(str)
    return tumor


def records_from_prins(
    counts: pd.DataFrame,
    clinical: pd.DataFrame,
    geo: str,
    qc_pass: bool,
) -> list[dict]:
    cpm = cpm_table(counts)
    clinical = clinical.set_index("SampleID")
    records = []
    for sample_id in counts.columns:
        if sample_id not in clinical.index:
            raise RuntimeError(f"{geo} sample {sample_id} missing from NC2024 clinical table")
        row = clinical.loc[sample_id]
        raw = float(counts.loc["CLDN4", sample_id])
        value = float(cpm.loc["CLDN4", sample_id])
        records.append(
            {
                "cohort": geo,
                "sample_id": sample_id,
                "subject_id": re.sub(r"_.*$", "", sample_id.replace("Pt", "Pt")),
                "qc_pass": qc_pass,
                "arm": "neoadjuvant" if str(row["treatment_group"]).strip() == "A" else "adjuvant",
                "rna_timing": (
                    "on_treatment_after_one_pembro_dose"
                    if str(row["treatment_group"]).strip() == "A"
                    else "pretreatment"
                ),
                "os_days": None if pd.isna(row["OS"]) else float(row["OS"]),
                "os_event": None if pd.isna(row["OS.Event"]) else int(row["OS.Event"]),
                "pfs_days": None if pd.isna(row["PFS"]) else float(row["PFS"]),
                "pfs_event": None if pd.isna(row["PFS.Event"]) else int(row["PFS.Event"]),
                "age": None if pd.isna(row["Age"]) else float(row["Age"]),
                "sex": None if pd.isna(row["Gender"]) else str(row["Gender"]),
                "mgmt": None if pd.isna(row["MGMT"]) else str(row["MGMT"]),
                "idh": None if pd.isna(row["IDH"]) else str(row["IDH"]),
                "cldn4_raw": raw,
                "cldn4_cpm": value,
                "cldn4_log2cpm": float(np.log2(value + 1.0)),
                "orr_available": False,
            }
        )
    records.sort(key=lambda record: record["sample_id"])
    return records


def records_from_gse318645(counts: pd.DataFrame, clinical: pd.DataFrame, qc_pass: bool) -> list[dict]:
    cpm = cpm_table(counts)
    clinical = clinical.set_index("sample_id")
    missing = set(counts.columns) - set(clinical.index)
    if missing:
        raise RuntimeError(f"GSE318645 count columns missing clinical rows: {sorted(missing)}")
    arm_map = {
        "1": "nivolumab_ipilimumab_neoadjuvant",
        "2": "nivolumab_neoadjuvant",
        "3": "placebo_neoadjuvant_then_dual_adjuvant",
    }
    records = []
    for sample_id in counts.columns:
        row = clinical.loc[sample_id]
        raw = float(counts.loc["CLDN4", sample_id])
        value = float(cpm.loc["CLDN4", sample_id])
        records.append(
            {
                "cohort": "GSE318645",
                "sample_id": sample_id,
                "subject_id": sample_id,
                "qc_pass": qc_pass,
                "arm": arm_map.get(str(row["arm"]), str(row["arm"])),
                "rna_timing": (
                    "pretreatment"
                    if str(row["arm"]) == "3"
                    else "on_treatment_after_neoadjuvant_dose"
                ),
                "os_days": float(row["os_days"]),
                "os_event": int(row["os_event"]),
                "pfs_days": None,
                "pfs_event": None,
                "age": float(row["age"]),
                "sex": str(row["Sex"]),
                "mgmt": None,
                "idh": None,
                "cldn4_raw": raw,
                "cldn4_cpm": value,
                "cldn4_log2cpm": float(np.log2(value + 1.0)),
                "orr_available": False,
            }
        )
    records.sort(key=lambda record: record["sample_id"])
    return records


def write_samples(records: list[dict], path: Path) -> None:
    fieldnames = [
        "cohort",
        "sample_id",
        "subject_id",
        "qc_pass",
        "arm",
        "rna_timing",
        "os_days",
        "os_event",
        "pfs_days",
        "pfs_event",
        "age",
        "sex",
        "mgmt",
        "idh",
        "cldn4_raw",
        "cldn4_cpm",
        "cldn4_log2cpm",
        "orr_available",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)


def stable_jitter(key: str) -> float:
    digest = hashlib.sha256(key.encode()).digest()
    return (int.from_bytes(digest[:4], "big") / (2**32 - 1) - 0.5) * 0.30


def make_detectability_figure(counts: pd.DataFrame, path: Path) -> None:
    genes = ["CLDN4", "TACSTD2", "CLDN7", "EPCAM", "CD274", "PTPRC", "ACTB", "GFAP"]
    cpm = cpm_table(counts)
    fig, axis = plt.subplots(figsize=(7.2, 4.6))
    positions = []
    values = []
    labels = []
    for index, gene in enumerate(genes):
        if gene not in cpm.index:
            continue
        series = np.log2(cpm.loc[gene].to_numpy() + 1.0)
        values.append(series)
        positions.append(index)
        labels.append(gene)
        axis.scatter(
            np.full(series.size, index) + np.array([stable_jitter(f"{gene}-{i}") for i in range(series.size)]),
            series,
            s=18,
            alpha=0.55,
            color="#C65A46" if gene == "CLDN4" else "#68788A",
            edgecolors="white",
            linewidths=0.3,
            zorder=3,
        )
    box = axis.boxplot(
        values,
        positions=positions,
        widths=0.55,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "#111111", "linewidth": 1.4},
        whiskerprops={"color": "#333333"},
        capprops={"color": "#333333"},
    )
    for patch, gene in zip(box["boxes"], labels):
        color = "#C65A46" if gene == "CLDN4" else "#68788A"
        patch.set_facecolor(color)
        patch.set_alpha(0.16)
        patch.set_edgecolor(color)
    axis.axhline(np.log2(2.0), color="#888888", linestyle="--", linewidth=0.8)
    axis.set_xticks(positions, labels, rotation=30, ha="right")
    axis.set_ylabel("log2(CPM + 1)")
    axis.set_title(
        "GSE121810 recurrent GBM (pembrolizumab): CLDN4 is barely expressed",
        loc="left",
        fontsize=11,
    )
    axis.spines[["top", "right"]].set_visible(False)
    axis.grid(axis="y", color="#E5E5E5", linewidth=0.8)
    axis.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def make_km_figure(records: list[dict], results: dict, path: Path) -> None:
    time = np.array([record["os_days"] for record in records], dtype=float)
    event = np.array([record["os_event"] for record in records], dtype=int)
    x = np.array([record["cldn4_log2cpm"] for record in records], dtype=float)
    high = x >= np.median(x)
    fig, axis = plt.subplots(figsize=(5.6, 4.6))
    km = KaplanMeierFitter()
    for mask, label, color in [
        (~high, f"CLDN4-low (n={(~high).sum()})", "#68788A"),
        (high, f"CLDN4-high (n={high.sum()})", "#C65A46"),
    ]:
        km.fit(time[mask], event[mask], label=label)
        km.plot(ax=axis, ci_show=False, color=color, linewidth=2.0)
    p_value = results["os"]["median_split_logrank_p"]
    hr = results["os"]["median_split_cox"]["hr_high_vs_low"]
    axis.set_title(
        (
            "GSE121810: OS by tumor CLDN4 median split\n"
            f"log-rank p={p_value:.3f}; Cox HR(high vs low)={hr:.2f}"
        ),
        loc="left",
        fontsize=11,
        linespacing=1.4,
    )
    axis.set_xlabel("Overall survival (days)")
    axis.set_ylabel("Survival probability")
    axis.spines[["top", "right"]].set_visible(False)
    axis.grid(axis="y", color="#E5E5E5", linewidth=0.8)
    axis.set_axisbelow(True)
    axis.legend(frameon=False, loc="upper right")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def make_os365_figure(records: list[dict], results: dict, path: Path) -> None:
    groups = [
        ("OS < 365d", 0, "#68788A"),
        ("OS ≥ 365d", 1, "#C65A46"),
    ]
    fig, axis = plt.subplots(figsize=(5.2, 4.6))
    values = []
    for label, binary, color in groups:
        selected = [
            record
            for record in records
            if (record["os_days"] >= 365) == bool(binary)
        ]
        group_values = [record["cldn4_log2cpm"] for record in selected]
        values.append(group_values)
        axis.scatter(
            [binary + stable_jitter(record["sample_id"]) for record in selected],
            group_values,
            s=22,
            alpha=0.65,
            color=color,
            edgecolors="white",
            linewidths=0.35,
            zorder=3,
        )
    box = axis.boxplot(
        values,
        positions=[0, 1],
        widths=0.48,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "#111111", "linewidth": 1.6},
        whiskerprops={"color": "#333333"},
        capprops={"color": "#333333"},
    )
    for patch, (_, _, color) in zip(box["boxes"], groups):
        patch.set_facecolor(color)
        patch.set_alpha(0.18)
        patch.set_edgecolor(color)
    test = results["os"]["os365_expression_test"]["test"]
    axis.set_title(
        (
            "GSE121810: CLDN4 vs 1-year OS (not RECIST)\n"
            f"Two-sided Mann–Whitney p={test['p_value_two_sided']:.3f}; n={len(records)}"
        ),
        loc="left",
        fontsize=11,
        linespacing=1.4,
    )
    axis.set_ylabel("CLDN4 log2(CPM + 1)")
    axis.set_xticks([0, 1], [f"OS < 365d\nn={len(values[0])}", f"OS ≥ 365d\nn={len(values[1])}"])
    axis.spines[["top", "right"]].set_visible(False)
    axis.grid(axis="y", color="#E5E5E5", linewidth=0.8)
    axis.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def inventory() -> list[dict]:
    return [
        {
            "accession": "GSE121810",
            "citation": "Cloughesy TF et al. Nat Med 2019. doi:10.1038/s41591-018-0337-7",
            "n_tumor_rna": 29,
            "tissue": "recurrent GBM tumor bulk RNA-seq",
            "treatment": "pembrolizumab neoadjuvant vs adjuvant",
            "expression_public": True,
            "orr_public": False,
            "survival_public": True,
            "cldn4_in_matrix": True,
            "strand_ok": True,
            "used_as": "primary",
            "note": "Only analysis-grade public bulk tumor matrix with public OS/PFS join.",
        },
        {
            "accession": "GSE264695",
            "citation": "de Groot J et al. Nat Commun 2024. doi:10.1038/s41467-024-54326-7",
            "n_tumor_rna": 20,
            "tissue": "recurrent GBM tumor bulk RNA-seq",
            "treatment": "neoadjuvant pembrolizumab expansion",
            "expression_public": True,
            "orr_public": False,
            "survival_public": True,
            "cldn4_in_matrix": True,
            "strand_ok": False,
            "used_as": "qc_failed_sensitivity",
            "note": "HTSeq counts look reverse-stranded (high MT-ND6 / antisense fraction).",
        },
        {
            "accession": "GSE318645",
            "citation": "GEO GSE318645; NCT04606316 dual ICB window trial",
            "n_tumor_rna": 47,
            "tissue": "recurrent GBM tumor bulk RNA-seq",
            "treatment": "nivo+ipi vs nivo vs placebo neoadjuvant, then ICB",
            "expression_public": True,
            "orr_public": False,
            "survival_public": True,
            "cldn4_in_matrix": True,
            "strand_ok": False,
            "used_as": "qc_failed_sensitivity",
            "note": "OS is in GEO metadata, but gene-assigned fraction ~2% and strand QC fails.",
        },
        {
            "accession": "GSE226976",
            "citation": "CAPTIVE / DNX-2401 + pembrolizumab NanoString",
            "n_tumor_rna": 48,
            "tissue": "recurrent glioma NanoString PanCancer Immune",
            "treatment": "oncolytic adenovirus + pembrolizumab",
            "expression_public": True,
            "orr_public": False,
            "survival_public": False,
            "cldn4_in_matrix": False,
            "strand_ok": None,
            "used_as": "excluded",
            "note": "770-gene immune panel; CLDN4 is absent.",
        },
        {
            "accession": "Zhao2019_PRJNA482620",
            "citation": "Zhao J et al. Nat Med 2019. doi:10.1038/s41591-019-0349-y",
            "n_tumor_rna": 66,
            "tissue": "GBM tumor RNA-seq (raw SRA only)",
            "treatment": "nivolumab or pembrolizumab",
            "expression_public": False,
            "orr_public": False,
            "survival_public": False,
            "cldn4_in_matrix": None,
            "strand_ok": None,
            "used_as": "excluded",
            "note": "Raw SRA only; paper says processed data on request. Sample-response map is not a public joinable table.",
        },
        {
            "accession": "GSE154795",
            "citation": "Lee AH et al. Nat Med 2021 neoadjuvant PD-1 GBM scRNA",
            "n_tumor_rna": 40,
            "tissue": "GBM scRNA-seq",
            "treatment": "neoadjuvant PD-1",
            "expression_public": True,
            "orr_public": False,
            "survival_public": False,
            "cldn4_in_matrix": None,
            "strand_ok": None,
            "used_as": "excluded",
            "note": "Single-cell, not bulk ICI-response matrix; processed object is >1 GB.",
        },
        {
            "accession": "GSE287371",
            "citation": "Newly diagnosed GBM: RT + nivolumab + IDO1 inhibitor",
            "n_tumor_rna": 0,
            "tissue": "peripheral blood leukocytes",
            "treatment": "radiation + nivolumab + BMS-986205",
            "expression_public": True,
            "orr_public": False,
            "survival_public": False,
            "cldn4_in_matrix": None,
            "strand_ok": None,
            "used_as": "excluded",
            "note": "Blood, not tumor. CLDN4 is an epithelial/tumor gene.",
        },
        {
            "accession": "TCGA-GBM",
            "citation": "TCGA glioblastoma",
            "n_tumor_rna": None,
            "tissue": "untreated/standard-care GBM",
            "treatment": "not ICI",
            "expression_public": True,
            "orr_public": False,
            "survival_public": True,
            "cldn4_in_matrix": True,
            "strand_ok": True,
            "used_as": "excluded",
            "note": "Not an ICI-response cohort. Not used to test B5.",
        },
    ]


def write_inventory(rows: list[dict], path: Path) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache-dir", type=Path, default=Path(".cache/B5_GBM"))
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    paths = {name: obtain(source, args.cache_dir) for name, source in SOURCES.items()}
    counts_121810 = load_count_matrix(paths["GSE121810_counts"])
    counts_264695 = load_count_matrix(paths["GSE264695_counts"])
    counts_318645 = load_count_matrix(paths["GSE318645_counts"])
    clinical_nc = load_nc2024_clinical(paths["NC2024_clinical"])
    clinical_318645 = load_gse318645_clinical(paths["GSE318645_matrix"])

    qc = {
        "GSE121810": library_metrics(counts_121810),
        "GSE264695": library_metrics(counts_264695),
        "GSE318645": library_metrics(counts_318645),
    }
    if not qc["GSE121810"]["strand_ok"]:
        raise RuntimeError("GSE121810 unexpectedly failed strand QC")
    if qc["GSE264695"]["strand_ok"] or qc["GSE318645"]["strand_ok"]:
        raise RuntimeError("Expected GSE264695 and GSE318645 to fail strand QC")

    detect = {
        geo: {gene: gene_detectability(counts, gene) for gene in ["CLDN4", "TACSTD2", "CLDN7", "GFAP", "ACTB"]}
        for geo, counts in [
            ("GSE121810", counts_121810),
            ("GSE264695", counts_264695),
            ("GSE318645", counts_318645),
        ]
    }

    records_121810 = records_from_prins(
        counts_121810, clinical_nc, "GSE121810", qc_pass=True
    )
    records_264695 = records_from_prins(
        counts_264695, clinical_nc, "GSE264695", qc_pass=False
    )
    records_318645 = records_from_gse318645(
        counts_318645, clinical_318645, qc_pass=False
    )
    if len({record["sample_id"] for record in records_121810}) != 29:
        raise RuntimeError("Expected 29 GSE121810 samples")
    if any(record["orr_available"] for record in records_121810 + records_264695 + records_318645):
        raise RuntimeError("ORR labels were not supposed to be present")

    adjuvant = [record for record in records_121810 if record["arm"] == "adjuvant"]
    neoadjuvant = [record for record in records_121810 if record["arm"] == "neoadjuvant"]

    cpm_ge1 = [record for record in records_121810 if record["cldn4_cpm"] >= 1]
    drop_top2 = sorted(records_121810, key=lambda record: record["cldn4_cpm"])[:-2]
    drop_cpm_ge1 = [record for record in records_121810 if record["cldn4_cpm"] < 1]
    primary = {
        "all_samples": {
            "os": analyze_survival(
                records_121810, "os_days", "os_event", landmark_days=365
            ),
            "pfs": analyze_survival(records_121810, "pfs_days", "pfs_event"),
            "arm_adjusted_os_cox": cox_arm_adjusted(records_121810),
        },
        "pretreatment_adjuvant_arm_only": {
            "note": (
                "Closest pretreatment predictive test. RNA is from surgery before "
                "adjuvant pembrolizumab. n is small."
            ),
            "os": analyze_survival(adjuvant, "os_days", "os_event", landmark_days=365),
        },
        "on_treatment_neoadjuvant_arm_only": {
            "note": "RNA is after one pembrolizumab dose. Not a pretreatment biomarker test.",
            "os": analyze_survival(
                neoadjuvant, "os_days", "os_event", landmark_days=365
            ),
        },
        "detectability_gated_sensitivity": {
            "note": (
                "CLDN4 CPM>=1 in only 3/29 tumors. The unadjusted continuous Cox "
                "is therefore a test of a few samples above the detection floor, "
                "not of a graded epithelial program."
            ),
            "samples_cpm_ge_1": [
                {
                    "sample_id": record["sample_id"],
                    "arm": record["arm"],
                    "cldn4_cpm": record["cldn4_cpm"],
                    "os_days": record["os_days"],
                    "os_event": record["os_event"],
                }
                for record in sorted(cpm_ge1, key=lambda record: -record["cldn4_cpm"])
            ],
            "exclude_cpm_ge_1": analyze_survival(
                drop_cpm_ge1, "os_days", "os_event", landmark_days=365
            ),
            "exclude_two_highest_cldn4": analyze_survival(
                drop_top2, "os_days", "os_event", landmark_days=365
            ),
        },
    }
    sensitivity = {
        "GSE264695_qc_failed": {
            "os": analyze_survival(records_264695, "os_days", "os_event", landmark_days=365)
        },
        "GSE121810_plus_GSE264695_ignore_qc": {
            "note": "Combined Nat Commun 2024 RNA set. Do not interpret; extension matrix fails strand QC.",
            "os": analyze_survival(
                records_121810 + records_264695, "os_days", "os_event", landmark_days=365
            ),
        },
        "GSE318645_qc_failed": {
            "os": analyze_survival(records_318645, "os_days", "os_event", landmark_days=365)
        },
    }

    inventory_rows = inventory()
    results = {
        "claim_tested": (
            "B5 analog in GBM/glioma: higher tumor CLDN4 associates with worse "
            "ICI outcome. User-stated pan-cancer OR=0.42 is not tested here."
        ),
        "endpoint_honesty": {
            "radiographic_orr_found": False,
            "endpoint_used": "overall survival (primary), PFS (secondary), OS>=365d exploratory binary",
            "reason": (
                "GEO sample records and the Nat Commun 2024 clinical workbook "
                "provide OS/PFS, not RECIST CR/PR/SD/PD."
            ),
        },
        "source": {
            "GSE121810": SOURCES["GSE121810_counts"],
            "GSE264695": SOURCES["GSE264695_counts"],
            "GSE318645_counts": SOURCES["GSE318645_counts"],
            "GSE318645_matrix": SOURCES["GSE318645_matrix"],
            "clinical": {
                **SOURCES["NC2024_clinical"],
                "sheet": "SFig4C",
                "citation": "de Groot J et al. Nat Commun 2024;15:10792. doi:10.1038/s41467-024-54326-7",
            },
        },
        "qc": qc,
        "detectability": detect,
        "analysis_set": {
            "primary_n": len(records_121810),
            "primary_definition": (
                "All 29 GSE121810 recurrent GBM tumors with public OS. "
                "Responder/nonresponder RECIST labels are not available."
            ),
            "adjuvant_pretreatment_n": len(adjuvant),
            "neoadjuvant_on_treatment_n": len(neoadjuvant),
        },
        "primary_GSE121810": primary,
        "sensitivity_qc_failed": sensitivity,
        "inventory": inventory_rows,
        "conclusion": {
            "supports_B5_in_GBM": False,
            "reason": (
                "No public RECIST/ORR labels exist for GBM ICI RNA. In GSE121810, "
                "CLDN4 is below 1 CPM in 26/29 tumors. An unadjusted continuous Cox "
                "OS association (HR 2.52, p=0.019) is leverage from two samples "
                "above that floor and is gone after they are removed (p=0.55). "
                "Median-split OS, Spearman, OS>=365d, pretreatment adjuvant arm, "
                "and two larger wrong-strand matrices do not confirm a usable "
                "CLDN4 response biomarker in glioma."
            ),
        },
    }

    write_samples(
        records_121810 + records_264695 + records_318645,
        args.output_dir / "analysis_samples.csv",
    )
    write_inventory(inventory_rows, args.output_dir / "cohort_inventory.csv")
    with (args.output_dir / "results.json").open("w") as handle:
        json.dump(results, handle, indent=2)
        handle.write("\n")
    make_detectability_figure(counts_121810, args.output_dir / "CLDN4_detectability.png")
    make_km_figure(records_121810, primary["all_samples"], args.output_dir / "CLDN4_vs_OS_GSE121810.png")
    make_os365_figure(records_121810, primary["all_samples"], args.output_dir / "CLDN4_vs_OS365.png")
    print(json.dumps(
        {
            "qc_strand_ok": {k: v["strand_ok"] for k, v in qc.items()},
            "cldn4_cpm_median_GSE121810": detect["GSE121810"]["CLDN4"]["cpm_median"],
            "primary_os_cox": primary["all_samples"]["os"]["continuous_cox"],
            "primary_os_logrank_p": primary["all_samples"]["os"]["median_split_logrank_p"],
            "primary_os365_p": primary["all_samples"]["os"]["os365_expression_test"]["test"]["p_value_two_sided"],
        },
        indent=2,
    ))


if __name__ == "__main__":
    main()
