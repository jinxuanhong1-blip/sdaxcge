#!/usr/bin/env python3
"""Build the patient-level Spearman registry for the scRNA TACSTD2/CLDN4 vs T/NK meta.

GSE207422 A3 is ingested as given (no re-audit). Other must-try series use
public processed-GEO patient tables (this repo) or a fresh GSE241934 download.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from lib_stats import spearman

HERE = Path(__file__).resolve().parent
DATA = HERE / "data" / "existing"
COHORTS = HERE / "results" / "cohorts"
OUT = HERE / "results"


def _row(**kwargs) -> dict:
    return kwargs


def from_xy(cohort, gene, metric, immune, x, y, **extra) -> dict:
    rho, p, n = spearman(x, y)
    d = _row(
        cohort=cohort,
        gene=gene,
        metric=metric,
        immune=immune,
        n=n,
        rho=rho,
        p=p,
        **extra,
    )
    return d


def gse207422_given() -> list[dict]:
    given = json.loads((HERE / "given_a3.json").read_text())
    df = pd.read_csv(DATA / "GSE207422_drmref_patients.tsv", sep="\t")
    # Confirm TACSTD2 matches the given slide; do not replace it.
    rho, p, n = spearman(df["malig_TACSTD2_mean"], df["frac_tnk"])
    rows = [
        _row(
            cohort="GSE207422",
            gene="TACSTD2",
            metric="malig_mean_log1p_cp10k",
            immune="frac_tnk",
            n=given["n_patients"],
            rho=given["spearman_rho_TACSTD2_mean"],
            p=given["spearman_p_TACSTD2_mean"],
            primary=True,
            source="A3_GIVEN",
            malignant_def="DRMref Malignant cells (taken as given)",
            note="User A3 slide taken as given. Table check ρ={:.6f} p={:.6g} n={}.".format(rho, p, n),
        ),
        _row(
            cohort="GSE207422",
            gene="TACSTD2",
            metric="malig_pct_pos",
            immune="frac_tnk",
            n=given["n_patients"],
            rho=given["spearman_rho_TACSTD2_pct"],
            p=given["spearman_p_TACSTD2_pct"],
            primary=False,
            source="A3_GIVEN_secondary",
            malignant_def="DRMref Malignant cells (taken as given)",
            note="Secondary %pos on the same given 12-patient table.",
        ),
    ]
    # CLDN4 is not on the A3 slide; compute from the same given table.
    rows.append(
        from_xy(
            "GSE207422",
            "CLDN4",
            "malig_mean_log1p_cp10k",
            "frac_tnk",
            df["malig_CLDN4_mean"],
            df["frac_tnk"],
            primary=True,
            source="A3_table_CLDN4_not_on_slide",
            malignant_def="DRMref Malignant cells (same given table)",
            note="CLDN4 was not on the A3 TACSTD2 slide; computed from the given 12-patient table.",
        )
    )
    return rows


def gse205335() -> list[dict]:
    df = pd.read_csv(DATA / "GSE205335_patients.tsv", sep="\t")
    # All deposited patients with malignant + T/NK (n=22). Author malignant labels.
    rows = []
    for gene, col, primary in (
        ("TACSTD2", "mal_TACSTD2_mean", True),
        ("TACSTD2", "mal_TACSTD2_pct_pos", False),
        ("CLDN4", "mal_CLDN4_mean", True),
        ("CLDN4", "mal_CLDN4_pct_pos", False),
    ):
        rows.append(
            from_xy(
                "GSE205335",
                gene,
                col,
                "frac_tnk",
                df[col],
                df["frac_tnk"],
                primary=primary,
                source="processed_GEO_author_malignant",
                malignant_def="author lineage.sub == Malignant cells",
                note="n=22 patients with captured malignant cells; mixed histology/sites; RECIST not required.",
            )
        )
    nsclc = df[df["cancer_subtype"].isin(["ADC", "SQ"])]
    for gene, col in (("TACSTD2", "mal_TACSTD2_mean"), ("CLDN4", "mal_CLDN4_mean")):
        rows.append(
            from_xy(
                "GSE205335_NSCLC",
                gene,
                col,
                "frac_tnk",
                nsclc[col],
                nsclc["frac_tnk"],
                primary=False,
                source="processed_GEO_author_malignant",
                malignant_def="author malignant; ADC+SQ only",
                note="Sensitivity: drop SCLC/NUT.",
            )
        )
    return rows


def gse131907() -> list[dict]:
    df = pd.read_csv(DATA / "GSE131907_samples.tsv", sep="\t")
    tumor = {"tLung", "tL/B", "mLN", "PE", "mBrain"}
    tdf = df[df["origin"].isin(tumor)].copy()
    epi = tdf[tdf["n_epithelial"] >= 20]
    mal = tdf[tdf["n_malignant"] >= 20]
    rows = []
    # Primary: tumor-site epithelium (author malignant calls are mets-heavy; tLung often n_mal=0).
    for gene, col, primary in (
        ("TACSTD2", "epi_TACSTD2_mean", True),
        ("TACSTD2", "epi_TACSTD2_pct", False),
        ("CLDN4", "epi_CLDN4_mean", True),
        ("CLDN4", "epi_CLDN4_pct", False),
    ):
        rows.append(
            from_xy(
                "GSE131907_tumor_sites",
                gene,
                col,
                "frac_tnk",
                epi[col],
                epi["frac_tnk"],
                primary=primary,
                source="processed_GEO_author_annotation",
                malignant_def="author Epithelial cells at tumor sites (tLung/tL-B/mLN/PE/mBrain)",
                note="Tumor sites only. Author 'Malignant cells' subtype is mostly mets/PE; tLung n_mal=0 so epithelium is the usable score.",
            )
        )
    for gene, col in (
        ("TACSTD2", "mal_TACSTD2_mean"),
        ("CLDN4", "mal_CLDN4_mean"),
    ):
        rows.append(
            from_xy(
                "GSE131907_malignant_subtype",
                gene,
                col,
                "frac_tnk",
                mal[col],
                mal["frac_tnk"],
                primary=False,
                source="processed_GEO_author_malignant_subtype",
                malignant_def="author Cell_subtype == Malignant cells",
                note="Sensitivity; mostly metastases / PE, not primary tLung.",
            )
        )
    return rows


def gse253013() -> list[dict]:
    df = pd.read_csv(DATA / "GSE253013_patients.tsv", sep="\t")
    t = df[(df["tissue"] == "Tumor") & (df["eligible_malig"])].copy()
    rows = []
    for gene, col, primary in (
        ("TACSTD2", "TACSTD2_mean_log1p", True),
        ("TACSTD2", "TACSTD2_pct_pos", False),
        ("CLDN4", "CLDN4_mean_log1p", True),
        ("CLDN4", "CLDN4_pct_pos", False),
    ):
        rows.append(
            from_xy(
                "GSE253013",
                gene,
                col,
                "frac_tnk",
                t[col],
                t["tnk_fraction"],
                primary=primary,
                source="processed_GEO_malig_like",
                malignant_def="epithelial AND not normal-lung marker program (treatment-naive LUAD tumors)",
                note="9 tumor patients; no ICI labels.",
            )
        )
    return rows


def gse291670() -> list[dict]:
    df = pd.read_csv(DATA / "GSE291670_patients.tsv", sep="\t")
    rows = []
    for gene, col, primary in (
        ("TACSTD2", "mal_TACSTD2_mean_log1p_cp10k", True),
        ("TACSTD2", "mal_TACSTD2_pct_pos", False),
        ("CLDN4", "mal_CLDN4_mean_log1p_cp10k", True),
        ("CLDN4", "mal_CLDN4_pct_pos", False),
    ):
        rows.append(
            from_xy(
                "GSE291670",
                gene,
                col,
                "frac_tnk",
                df[col],
                df["frac_lineage_tnk"],
                primary=primary,
                source="processed_GEO_marker_malignant",
                malignant_def="marker-score malignant (n=6 post anlotinib+PD-1)",
                note="n=6 (3 MPR / 3 NMPR). Lineage T/NK fraction.",
            )
        )
    return rows


def gse241934() -> list[dict]:
    rows = []
    for cohort, label, note in (
        (
            "GSE241934_IIT",
            "GSE241934_IIT",
            "EGFR-mutant NEOTIDE IIT residual tumors; author Epi vs T+NK; all 11 patients eligible.",
        ),
        (
            "GSE241934_Real",
            "GSE241934_Real",
            "EGFR-WT real-world residual tumors; author Epi vs T+NK; ≥20 Epi and ≥20 T/NK.",
        ),
    ):
        pdf = pd.read_csv(COHORTS / f"{cohort}_patients.tsv", sep="\t")
        elig = pdf[pdf["eligible"] == 1]
        for gene, col, primary in (
            ( "TACSTD2", "TACSTD2_log1p_cp10k", True),
            ( "TACSTD2", "TACSTD2_pct_pos", False),
            ( "CLDN4", "CLDN4_log1p_cp10k", True),
            ( "CLDN4", "CLDN4_pct_pos", False),
        ):
            rows.append(
                from_xy(
                    label,
                    gene,
                    col,
                    "frac_tnk",
                    elig[col],
                    elig["frac_tnk"],
                    primary=primary,
                    source="processed_GEO_author_Epi",
                    malignant_def="author major.cell.type == Epi (residual epithelium; no CNV call)",
                    note=note,
                )
            )
    return rows


def gse325414() -> list[dict]:
    df = pd.read_csv(DATA / "GSE325414_donors.csv")
    rows = []
    for gene, col in (("TACSTD2", "malignant_TACSTD2_mean"), ("CLDN4", "malignant_CLDN4_mean")):
        rows.append(
            from_xy(
                "GSE325414",
                gene,
                col,
                "frac_tnk",
                df[col],
                df["frac_TNK"],
                primary=True,
                source="leftover_2026_author_labels",
                malignant_def="author malignant labels; PEF ablation treat-and-resect (not ICI)",
                note="2026 leftover lung tumor scRNA. Donor-level (n=25).",
            )
        )
    return rows


def skip_rows() -> list[dict]:
    return [
        _row(
            cohort="GSE146100",
            gene="TACSTD2",
            metric="NA",
            immune="NA",
            n=1,
            rho=float("nan"),
            p=float("nan"),
            primary=False,
            source="skip",
            malignant_def="cluster epithelial (1 patient, 3 nodules)",
            note="n=1 patient. Not usable for a patient-level Spearman (n>1 required).",
        ),
        _row(
            cohort="GSE243013",
            gene="both",
            metric="NA",
            immune="NA",
            n=0,
            rho=float("nan"),
            p=float("nan"),
            primary=False,
            source="skip",
            malignant_def="immune-cell-only atlas",
            note="Deposited matrix is immune cells, not malignant/epithelial. No T/NK-vs-malignant score.",
        ),
        _row(
            cohort="EGA/dbGaP raw (GSE205335 raw EGAD00001008703; GSE207422 HRA001033)",
            gene="both",
            metric="NA",
            immune="NA",
            n=0,
            rho=float("nan"),
            p=float("nan"),
            primary=False,
            source="skip",
            malignant_def="NA",
            note="Controlled raw skipped by policy. Processed GEO used instead.",
        ),
    ]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    rows += gse207422_given()
    rows += gse205335()
    rows += gse241934()
    rows += gse291670()
    rows += gse253013()
    rows += gse131907()
    rows += gse325414()
    rows += skip_rows()
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "cohort_effects.tsv", sep="\t", index=False)
    print(df[df["primary"] == True][["cohort", "gene", "n", "rho", "p"]].to_string(index=False))
    print("wrote", OUT / "cohort_effects.tsv")


if __name__ == "__main__":
    main()
