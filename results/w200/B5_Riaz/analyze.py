#!/usr/bin/env python3
"""Reproduce the Riaz 2017 CLDN4-response analysis from public GEO files."""

from __future__ import annotations

import csv
import gzip
import hashlib
import math
import re
import urllib.request
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
FPKM = DATA / "GSE91061_fpkm.csv.gz"
MATRIX = DATA / "GSE91061_series_matrix.txt.gz"
URLS = {
    FPKM: (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE91nnn/GSE91061/suppl/"
        "GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv.gz"
    ),
    MATRIX: (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE91nnn/GSE91061/matrix/"
        "GSE91061_series_matrix.txt.gz"
    ),
}
CLDN4_ENTREZ = "1364"
KNOWN_RESPONSES = {"PRCR", "SD", "PD"}
SEED = 91061


def download_inputs() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    for path, url in URLS.items():
        if not path.exists():
            print(f"Downloading {url}")
            urllib.request.urlretrieve(url, path)


def parse_metadata() -> dict[str, dict[str, str]]:
    sample_rows: list[list[str]] = []
    with gzip.open(MATRIX, "rt") as handle:
        for line in handle:
            if line.startswith("!Sample_"):
                sample_rows.append(next(csv.reader([line], delimiter="\t")))

    def one_row(name: str) -> list[str]:
        matches = [row[1:] for row in sample_rows if row[0] == name]
        if len(matches) != 1:
            raise ValueError(f"Expected exactly one {name} row; found {len(matches)}")
        return matches[0]

    titles = one_row("!Sample_title")
    accessions = one_row("!Sample_geo_accession")
    characteristics = [
        row[1:] for row in sample_rows if row[0] == "!Sample_characteristics_ch1"
    ]
    characteristic_map: dict[str, list[str]] = {}
    for values in characteristics:
        key = values[0].split(":", 1)[0]
        characteristic_map[key] = [value.split(": ", 1)[1] for value in values]

    if not (len(titles) == len(accessions) == 109):
        raise ValueError("Unexpected GEO sample count")

    metadata = {}
    for title, gsm, visit, response, tissue in zip(
        titles,
        accessions,
        characteristic_map["visit (pre or on treatment)"],
        characteristic_map["response"],
        characteristic_map["tissue"],
        strict=True,
    ):
        match = re.match(r"Pt([^_]+)_", title)
        if not match:
            raise ValueError(f"Cannot derive patient ID from {title}")
        metadata[title] = {
            "sample": title,
            "gsm": gsm,
            "patient": f"Pt{match.group(1)}",
            "visit": visit,
            "response": response,
            "tissue": tissue,
        }
    return metadata


def parse_cldn4() -> dict[str, float]:
    with gzip.open(FPKM, "rt", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        matches = [row for row in reader if row and row[0] == CLDN4_ENTREZ]
    if len(matches) != 1:
        raise ValueError(f"Expected one Entrez {CLDN4_ENTREZ} row; found {len(matches)}")
    values = [float(value) for value in matches[0][1:]]
    if len(header[1:]) != 109 or len(values) != 109:
        raise ValueError("Unexpected expression matrix dimensions")
    return dict(zip(header[1:], values, strict=True))


def bootstrap_auc_ci(
    responders: np.ndarray, nonresponders: np.ndarray, draws: int = 20_000
) -> tuple[float, float]:
    rng = np.random.default_rng(SEED)
    estimates = np.empty(draws)
    for i in range(draws):
        a = rng.choice(responders, len(responders), replace=True)
        b = rng.choice(nonresponders, len(nonresponders), replace=True)
        estimates[i] = stats.mannwhitneyu(a, b).statistic / (len(a) * len(b))
    return tuple(np.quantile(estimates, [0.025, 0.975]))


def summarize(groups: dict[str, np.ndarray]) -> dict[str, float | int | str]:
    responders = groups["Responder"]
    nonresponders = groups["Nonresponder"]
    test = stats.mannwhitneyu(
        responders, nonresponders, alternative="two-sided", method="asymptotic"
    )
    auc = float(test.statistic / (len(responders) * len(nonresponders)))
    auc_lo, auc_hi = bootstrap_auc_ci(responders, nonresponders)
    detected_r = int(np.sum(responders > 0))
    detected_nr = int(np.sum(nonresponders > 0))
    fisher = stats.fisher_exact(
        [
            [detected_r, len(responders) - detected_r],
            [detected_nr, len(nonresponders) - detected_nr],
        ],
        alternative="two-sided",
    )
    return {
        "n_responder": len(responders),
        "n_nonresponder": len(nonresponders),
        "median_responder": float(np.median(responders)),
        "median_nonresponder": float(np.median(nonresponders)),
        "mann_whitney_u": float(test.statistic),
        "mann_whitney_p": float(test.pvalue),
        "auc_responder_higher": auc,
        "auc_bootstrap_95ci_low": auc_lo,
        "auc_bootstrap_95ci_high": auc_hi,
        "detected_responder": detected_r,
        "detected_nonresponder": detected_nr,
        "detection_fisher_odds_ratio": float(fisher.statistic),
        "detection_fisher_p": float(fisher.pvalue),
    }


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def make_figure(
    baseline: dict[str, np.ndarray],
    on_treatment: dict[str, np.ndarray],
    changes: dict[str, np.ndarray],
    summaries: dict[str, dict],
) -> None:
    plt.style.use("seaborn-v0_8-whitegrid")
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.1))
    rng = np.random.default_rng(SEED)
    panels = [
        ("Pretreatment (primary)", baseline, "log2(FPKM + 1)", "pretreatment"),
        ("On-treatment (exploratory)", on_treatment, "log2(FPKM + 1)", "on_treatment"),
        ("Paired on − pre (exploratory)", changes, "Δ log2(FPKM + 1)", "paired_change"),
    ]
    colors = ["#6b7280", "#2563eb"]
    for ax, (title, groups, ylabel, key) in zip(axes, panels, strict=True):
        ordered = [groups["Nonresponder"], groups["Responder"]]
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
        ax.set_xticks([1, 2], [f"NR\nn={len(ordered[0])}", f"R\nn={len(ordered[1])}"])
        ax.set_title(title, fontsize=10, weight="bold")
        ax.set_ylabel(ylabel)
        ax.text(
            0.5,
            0.98,
            f"MW p={summaries[key]['mann_whitney_p']:.3g}",
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=9,
        )
    fig.suptitle("Riaz 2017 melanoma: CLDN4 expression versus nivolumab response", y=1.02)
    fig.tight_layout()
    fig.savefig(ROOT / "cldn4_vs_response.png", dpi=220, bbox_inches="tight")
    fig.savefig(ROOT / "cldn4_vs_response.svg", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    download_inputs()
    metadata = parse_metadata()
    expression = parse_cldn4()
    if set(expression) != set(metadata):
        raise ValueError("Expression and metadata sample names do not match")

    sample_rows = []
    for sample, fpkm in expression.items():
        row = dict(metadata[sample])
        row["cldn4_fpkm"] = fpkm
        row["cldn4_log2_fpkm_plus_1"] = math.log2(fpkm + 1)
        row["binary_group"] = (
            "Responder"
            if row["response"] == "PRCR"
            else "Nonresponder"
            if row["response"] in {"SD", "PD"}
            else ""
        )
        sample_rows.append(row)
    sample_rows.sort(key=lambda row: (row["patient"], row["visit"]))

    eligible = [
        row for row in sample_rows if row["response"] in KNOWN_RESPONSES
    ]

    def visit_groups(visit: str) -> dict[str, np.ndarray]:
        return {
            group: np.array(
                [
                    row["cldn4_log2_fpkm_plus_1"]
                    for row in eligible
                    if row["visit"] == visit and row["binary_group"] == group
                ]
            )
            for group in ("Responder", "Nonresponder")
        }

    baseline = visit_groups("Pre")
    on_treatment = visit_groups("On")

    by_patient: dict[str, dict[str, dict]] = {}
    for row in eligible:
        by_patient.setdefault(row["patient"], {})[row["visit"]] = row
    change_rows = []
    for patient, visits in by_patient.items():
        if "Pre" not in visits or "On" not in visits:
            continue
        pre, on = visits["Pre"], visits["On"]
        if pre["response"] != on["response"]:
            raise ValueError(f"Inconsistent response for paired patient {patient}")
        change_rows.append(
            {
                "patient": patient,
                "response": pre["response"],
                "binary_group": pre["binary_group"],
                "pre_log2_fpkm_plus_1": pre["cldn4_log2_fpkm_plus_1"],
                "on_log2_fpkm_plus_1": on["cldn4_log2_fpkm_plus_1"],
                "change_on_minus_pre": (
                    on["cldn4_log2_fpkm_plus_1"]
                    - pre["cldn4_log2_fpkm_plus_1"]
                ),
            }
        )
    changes = {
        group: np.array(
            [
                row["change_on_minus_pre"]
                for row in change_rows
                if row["binary_group"] == group
            ]
        )
        for group in ("Responder", "Nonresponder")
    }

    summaries = {
        "pretreatment": summarize(baseline),
        "on_treatment": summarize(on_treatment),
        "paired_change": summarize(changes),
    }
    summary_rows = [
        {"analysis": analysis, **values} for analysis, values in summaries.items()
    ]
    write_csv(ROOT / "summary.csv", summary_rows, list(summary_rows[0]))
    write_csv(ROOT / "sample_level.csv", sample_rows, list(sample_rows[0]))
    write_csv(ROOT / "paired_changes.csv", change_rows, list(change_rows[0]))
    make_figure(baseline, on_treatment, changes, summaries)

    with (ROOT / "provenance.tsv").open("w") as handle:
        handle.write("file\tsha256\turl\n")
        for path, url in URLS.items():
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            handle.write(f"{path.name}\t{digest}\t{url}\n")

    print("GEO response counts:", dict(Counter(row["response"] for row in sample_rows)))
    for analysis, values in summaries.items():
        print(
            f"{analysis}: R={values['n_responder']}, "
            f"NR={values['n_nonresponder']}, "
            f"AUC={values['auc_responder_higher']:.3f}, "
            f"MW p={values['mann_whitney_p']:.4g}"
        )


if __name__ == "__main__":
    main()
