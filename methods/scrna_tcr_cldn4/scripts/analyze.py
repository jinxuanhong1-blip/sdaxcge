#!/usr/bin/env python3
"""Patient-level TCR clone expansion vs author-epithelial CLDN4 on GSE241934.

Public GEO only. Unit of inference = patient. GSE207422 has no TCR/BCR.
"""
from __future__ import annotations

import csv
import gzip
import json
import math
import os
import re
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
INPUTS = ROOT / "inputs"
RESULTS = ROOT / "results"
TCR_DIR = Path(os.environ.get("TCR_DIR", "/tmp/gse241934/tcr"))
MIN_TCR = 50
MIN_EPI_PRIMARY = 20


def shannon_clonality(freqs: list[int]) -> float:
    n = sum(freqs)
    k = len(freqs)
    if n <= 0 or k <= 1:
        return 0.0
    p = np.asarray(freqs, dtype=float) / n
    h = float(-(p * np.log(p)).sum())
    return 1.0 - h / math.log(k)


def cliffs_delta(a: np.ndarray, b: np.ndarray) -> float:
    """δ = P(a>b) - P(a<b). Positive => a larger than b."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.size == 0 or b.size == 0:
        return float("nan")
    gt = 0
    lt = 0
    for x in a:
        gt += int(np.sum(x > b))
        lt += int(np.sum(x < b))
    return (gt - lt) / (a.size * b.size)


def load_tcr_patient(clonotypes_path: Path, contig_path: Path) -> dict:
    with gzip.open(clonotypes_path, "rt") as f:
        rows = list(csv.DictReader(f))
    freqs = [int(r["frequency"]) for r in rows]
    n_cells = int(sum(freqs))
    n_clones = len(freqs)
    exp2 = [x for x in freqs if x >= 2]
    exp3 = [x for x in freqs if x >= 3]
    top = max(freqs) if freqs else 0

    chains: Counter[str] = Counter()
    barcodes: set[str] = set()
    paired: set[str] = set()
    tra: set[str] = set()
    trb: set[str] = set()
    if contig_path.exists():
        with gzip.open(contig_path, "rt") as f:
            for r in csv.DictReader(f):
                chains[r["chain"]] += 1
                bc = r["barcode"]
                barcodes.add(bc)
                if r["chain"] == "TRA":
                    tra.add(bc)
                elif r["chain"] == "TRB":
                    trb.add(bc)
        paired = tra & trb

    return {
        "n_tcr_cells": n_cells,
        "n_clonotypes": n_clones,
        "n_expanded_clones_ge2": len(exp2),
        "n_expanded_clones_ge3": len(exp3),
        "expanded_cell_frac_ge2": (sum(exp2) / n_cells) if n_cells else float("nan"),
        "expanded_cell_frac_ge3": (sum(exp3) / n_cells) if n_cells else float("nan"),
        "expanded_clones_per_1k_ge2": (1000.0 * len(exp2) / n_cells) if n_cells else float("nan"),
        "top_clone_frac": (top / n_cells) if n_cells else float("nan"),
        "clonality_shannon": shannon_clonality(freqs),
        "n_contig_barcodes": len(barcodes),
        "n_paired_tra_trb": len(paired),
        "n_tra_contigs": int(chains.get("TRA", 0)),
        "n_trb_contigs": int(chains.get("TRB", 0)),
        "n_igh_igk_igl_contigs": int(
            chains.get("IGH", 0) + chains.get("IGK", 0) + chains.get("IGL", 0)
        ),
        "contig_chains": dict(chains),
    }


def load_all_tcr() -> pd.DataFrame:
    recs = []
    for p in sorted(TCR_DIR.glob("*_clonotypes.csv.gz")):
        m = re.search(r"(GSM\d+)_(P\d+)_clonotypes", p.name)
        if not m:
            raise ValueError(f"unparsed TCR filename: {p.name}")
        gsm, patient = m.group(1), m.group(2)
        contig = p.with_name(f"{gsm}_{patient}_filtered_contig_annotations.csv.gz")
        rec = load_tcr_patient(p, contig)
        rec["sampleID"] = patient
        rec["gsm"] = gsm
        recs.append(rec)
    return pd.DataFrame(recs)


def spearman_row(x, y, test: str, n_note: str) -> dict:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    n = int(x.size)
    if n < 3:
        return {
            "test": test,
            "n": n,
            "rho": float("nan"),
            "p": float("nan"),
            "note": n_note + "; n<3",
        }
    rho, p = stats.spearmanr(x, y)
    return {
        "test": test,
        "n": n,
        "rho": float(rho),
        "p": float(p),
        "note": n_note,
    }


def mwu_row(a, b, test: str, note: str) -> dict:
    a = np.asarray(a, dtype=float)
    a = a[np.isfinite(a)]
    b = np.asarray(b, dtype=float)
    b = b[np.isfinite(b)]
    if a.size == 0 or b.size == 0:
        return {
            "test": test,
            "n_a": int(a.size),
            "n_b": int(b.size),
            "median_a": float("nan"),
            "median_b": float("nan"),
            "U": float("nan"),
            "p": float("nan"),
            "cliffs_delta": float("nan"),
            "note": note,
        }
    U, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    return {
        "test": test,
        "n_a": int(a.size),
        "n_b": int(b.size),
        "median_a": float(np.median(a)),
        "median_b": float(np.median(b)),
        "U": float(U),
        "p": float(p),
        "cliffs_delta": float(cliffs_delta(a, b)),
        "note": note,
    }


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    tcr = load_all_tcr()
    cldn = pd.read_csv(INPUTS / "GSE241934_author_epi_cldn4.tsv", sep="\t")
    df = tcr.merge(cldn, on="sampleID", how="outer", indicator=True)
    df["has_tcr"] = df["n_tcr_cells"].notna()
    df["has_cldn4"] = df["malignant_CLDN4_mean_log1p"].notna()
    df["join"] = df["_merge"].map(
        {
            "both": "tcr+cldn4_row",
            "left_only": "tcr_only",
            "right_only": "cldn4_only",
        }
    )
    df = df.drop(columns=["_merge"])
    df["pass_tcr"] = df["n_tcr_cells"].fillna(0) >= MIN_TCR
    df["pass_epi20"] = df["n_epi"].fillna(0) >= MIN_EPI_PRIMARY
    df["pass_primary"] = df["pass_tcr"] & df["pass_epi20"] & df["has_cldn4"]

    df.to_csv(RESULTS / "patient_tcr_cldn4.tsv", sep="\t", index=False)

    primary = df[df["pass_primary"]].copy()
    tests = []
    tests.append(
        spearman_row(
            primary["malignant_CLDN4_mean_log1p"],
            primary["expanded_cell_frac_ge2"],
            "CLDN4 vs expanded_cell_frac (clone≥2)",
            f"author Epi ≥{MIN_EPI_PRIMARY}; TCR cells ≥{MIN_TCR}; GSE241934 IIT+Real",
        )
    )
    tests.append(
        spearman_row(
            primary["malignant_CLDN4_mean_log1p"],
            primary["expanded_clones_per_1k_ge2"],
            "CLDN4 vs expanded_clones_per_1k (clone≥2)",
            "same primary floor",
        )
    )
    tests.append(
        spearman_row(
            primary["malignant_CLDN4_mean_log1p"],
            primary["clonality_shannon"],
            "CLDN4 vs Shannon clonality",
            "same primary floor",
        )
    )
    tests.append(
        spearman_row(
            primary["malignant_CLDN4_mean_log1p"],
            primary["top_clone_frac"],
            "CLDN4 vs top_clone_frac",
            "same primary floor",
        )
    )
    tests.append(
        spearman_row(
            primary["malignant_CLDN4_mean_log1p"],
            primary["expanded_cell_frac_ge3"],
            "CLDN4 vs expanded_cell_frac (clone≥3)",
            "same primary floor; expansion cutoff ≥3",
        )
    )
    tests.append(
        spearman_row(
            primary["malignant_TACSTD2_mean_log1p"],
            primary["expanded_cell_frac_ge2"],
            "TACSTD2 vs expanded_cell_frac (clone≥2)",
            "secondary gene; same primary floor",
        )
    )

    for floor in (10, 5):
        sub = df[
            df["pass_tcr"]
            & df["has_cldn4"]
            & (df["n_epi"].fillna(0) >= floor)
        ]
        tests.append(
            spearman_row(
                sub["malignant_CLDN4_mean_log1p"],
                sub["expanded_cell_frac_ge2"],
                f"CLDN4 vs expanded_cell_frac (clone≥2), epi≥{floor}",
                "sensitivity floor",
            )
        )

    for cohort, label in (
        ("IIT_EGFRmut", "IIT EGFR-mut"),
        ("REAL_WT", "Real-world EGFR-WT"),
    ):
        sub = primary[primary["cohort"] == cohort]
        tests.append(
            spearman_row(
                sub["malignant_CLDN4_mean_log1p"],
                sub["expanded_cell_frac_ge2"],
                f"CLDN4 vs expanded_cell_frac (clone≥2), {label}",
                "primary floor, one cohort",
            )
        )

    mpr = primary[primary["group"] == "MPR"]
    nmpr = primary[primary["group"] == "NMPR"]
    tests.append(
        mwu_row(
            nmpr["expanded_cell_frac_ge2"],
            mpr["expanded_cell_frac_ge2"],
            "expanded_cell_frac NMPR vs MPR",
            "primary floor; n_a=NMPR, n_b=MPR; secondary (not the CLDN4 test)",
        )
    )
    tests.append(
        mwu_row(
            nmpr["malignant_CLDN4_mean_log1p"],
            mpr["malignant_CLDN4_mean_log1p"],
            "author-Epi CLDN4 NMPR vs MPR",
            "primary floor; n_a=NMPR, n_b=MPR; secondary",
        )
    )

    pd.DataFrame(tests).to_csv(RESULTS / "stats.tsv", sep="\t", index=False)

    hunt = [
        {
            "accession": "GSE207422",
            "public_tcr_bcr": "no",
            "malignant_or_epi_CLDN4": "yes (UMI; 15 patients)",
            "same_patient_tcr_and_cldn4": "no",
            "why": "BD Rhapsody WTA + bulk RNA only. GEO suppl: UMI matrix + sample xlsx + bulk TPM. No contig/clonotype/VDJ. Paper methods: Rhapsody WTA, no TCR/BCR library. Raw not deposited.",
        },
        {
            "accession": "GSE241934",
            "public_tcr_bcr": "yes TCR (43 patients; TRA/TRB). No BCR (0 IGH/IGK/IGL contigs).",
            "malignant_or_epi_CLDN4": "yes (author Epi on IIT+Real MTX)",
            "same_patient_tcr_and_cldn4": "yes",
            "why": "Leftover neoadjuvant IO+chemo scRNA. RAW.tar = Cell Ranger clonotypes + contig CSVs. Same sampleIDs as author-Epi CLDN4 table.",
        },
        {
            "accession": "GSE243013",
            "public_tcr_bcr": "yes TCR",
            "malignant_or_epi_CLDN4": "no (CD45+ immune MTX only)",
            "same_patient_tcr_and_cldn4": "no",
            "why": "Immune-only. Residual TACSTD2/CLDN4 in CD45+ is not malignant epithelium.",
        },
        {
            "accession": "GSE179994",
            "public_tcr_bcr": "yes TCR",
            "malignant_or_epi_CLDN4": "no (T cells only)",
            "same_patient_tcr_and_cldn4": "no",
            "why": "Public tables are T-cell scTCR. Cannot score malignant CLDN4.",
        },
        {
            "accession": "GSE176021 / GSE176022",
            "public_tcr_bcr": "yes (VDJ tars / bulk culture TCR)",
            "malignant_or_epi_CLDN4": "no",
            "same_patient_tcr_and_cldn4": "no",
            "why": "GEO processed = lymphocytes. GSE176022 is bulk MANAFEST TCR, not scRNA.",
        },
        {
            "accession": "GSE185204 / GSE185206 / GSE186446",
            "public_tcr_bcr": "yes TCR (n=3)",
            "malignant_or_epi_CLDN4": "no (sorted T cells)",
            "same_patient_tcr_and_cldn4": "no",
            "why": "CD3+ TIL atlas. n=3. No malignant compartment.",
        },
        {
            "accession": "GSE280232",
            "public_tcr_bcr": "yes (paired GEX+TCR GSMs)",
            "malignant_or_epi_CLDN4": "no",
            "same_patient_tcr_and_cldn4": "no",
            "why": "GEO characteristic cell type = Sorted T cells. TIL-only.",
        },
        {
            "accession": "GSE229353",
            "public_tcr_bcr": "no public contig/clonotype",
            "malignant_or_epi_CLDN4": "no (CD45+ MTX)",
            "same_patient_tcr_and_cldn4": "no",
            "why": "Post-neoadjuvant CD45+ GEX only.",
        },
        {
            "accession": "GSE291670",
            "public_tcr_bcr": "no",
            "malignant_or_epi_CLDN4": "yes (snRNA MTX)",
            "same_patient_tcr_and_cldn4": "no",
            "why": "RAW.tar is barcodes/features/matrix only. No VDJ.",
        },
        {
            "accession": "GSE205335",
            "public_tcr_bcr": "no",
            "malignant_or_epi_CLDN4": "yes (UMI RDS + CellIdentity)",
            "same_patient_tcr_and_cldn4": "no",
            "why": "GEO suppl: UMI RDS + identity table. No TCR.",
        },
        {
            "accession": "GSE267108 / GSE274595 / GSE337519 / E-MTAB-13526",
            "public_tcr_bcr": "no",
            "malignant_or_epi_CLDN4": "epithelial scores exist on some",
            "same_patient_tcr_and_cldn4": "no",
            "why": "2024–26 leftover epithelial matrices without public VDJ.",
        },
        {
            "accession": "GSE308745",
            "public_tcr_bcr": "yes TCR+ADT",
            "malignant_or_epi_CLDN4": "no (PBMC)",
            "same_patient_tcr_and_cldn4": "no",
            "why": "Blood only.",
        },
    ]
    pd.DataFrame(hunt).to_csv(RESULTS / "hunt_tcr_cldn4.tsv", sep="\t", index=False)

    summary = {
        "gse207422_public_tcr_bcr": False,
        "additive_series": "GSE241934",
        "n_tcr_patients": int(tcr.shape[0]),
        "n_cldn4_table_patients": int(cldn.shape[0]),
        "n_join_both": int((df["join"] == "tcr+cldn4_row").sum()),
        "n_tcr_with_cldn4_value": int((df["has_tcr"] & df["has_cldn4"]).sum()),
        "n_primary": int(primary.shape[0]),
        "n_primary_mpr": int((primary["group"] == "MPR").sum()),
        "n_primary_nmpr": int((primary["group"] == "NMPR").sum()),
        "n_primary_iit": int((primary["cohort"] == "IIT_EGFRmut").sum()),
        "n_primary_real": int((primary["cohort"] == "REAL_WT").sum()),
        "bcr_contigs_total": int(tcr["n_igh_igk_igl_contigs"].sum()),
        "min_tcr_cells": MIN_TCR,
        "min_epi_primary": MIN_EPI_PRIMARY,
        "cldn4_source": "author epithelial cells on public GSE241934 IIT+Real MTX; mean log1p from methods/scrna_meta_mpr/results/per_patient_GSE241934.tsv",
        "primary_test": tests[0],
    }
    (RESULTS / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    fig, ax = plt.subplots(figsize=(6.2, 4.8))
    colors = {"MPR": "#2a6f97", "NMPR": "#c1121f"}
    markers = {"IIT_EGFRmut": "o", "REAL_WT": "s"}
    for (grp, cohort), sub in primary.groupby(["group", "cohort"]):
        ax.scatter(
            sub["malignant_CLDN4_mean_log1p"],
            sub["expanded_cell_frac_ge2"],
            c=colors.get(grp, "0.4"),
            marker=markers.get(cohort, "o"),
            s=42,
            alpha=0.85,
            edgecolors="white",
            linewidths=0.4,
            label=f"{grp} / {cohort.replace('_', ' ')} (n={len(sub)})",
        )
    rho = tests[0]["rho"]
    p = tests[0]["p"]
    ax.set_xlabel("Author-epithelial CLDN4 (mean log1p)")
    ax.set_ylabel("Expanded TCR-cell fraction (clone size ≥2)")
    ax.set_title(
        f"GSE241934  n={tests[0]['n']}  Spearman ρ={rho:.3f}  p={p:.3g}"
    )
    ax.legend(fontsize=7, loc="best", frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(RESULTS / "fig_cldn4_vs_tcr_expansion.png", dpi=160)
    fig.savefig(RESULTS / "fig_cldn4_vs_tcr_expansion.pdf")
    plt.close(fig)

    print(json.dumps(summary, indent=2))
    print(pd.DataFrame(tests).to_string(index=False))


if __name__ == "__main__":
    main()
