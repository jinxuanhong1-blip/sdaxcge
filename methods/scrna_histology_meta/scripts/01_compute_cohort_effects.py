#!/usr/bin/env python3
"""Compute histology-stratified malignant/epithelial TACSTD2/CLDN4 vs T/NK Spearman.

Reads harvested sibling sample tables plus optional GSE241934 extract.
Does not invent histology. ASC / SCLC / NUT / UNKNOWN stay out of LUAD and LUSC pools.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from lib import map_histology, rho_ci, spearman_rho_p


ROOT = Path(__file__).resolve().parents[1]
HARVESTED = ROOT / "harvested"


def effect_row(
    *,
    series: str,
    cohort: str,
    histology: str,
    gene: str,
    compartment: str,
    n: int,
    rho: float,
    p: float,
    min_epi: int,
    min_tnk: int,
    note: str,
    primary: bool,
) -> dict:
    lo, hi = rho_ci(rho, n)
    return {
        "series": series,
        "cohort": cohort,
        "histology": histology,
        "gene": gene,
        "compartment": compartment,
        "n": n,
        "rho": rho,
        "p": p,
        "ci_low": lo,
        "ci_high": hi,
        "min_epi_or_mal": min_epi,
        "min_tnk": min_tnk,
        "primary_meta": int(bool(primary and n >= 6 and histology in {"LUAD", "LUSC"})),
        "note": note,
    }


def add_gene_effects(rows, **common):
    x = common.pop("x")
    y = common.pop("y")
    rho, p, n = spearman_rho_p(x, y)
    rows.append(effect_row(n=n, rho=rho, p=p, **common))


def gse207422(rows: list[dict]) -> None:
    df = pd.read_csv(HARVESTED / "gse207422_sample_table.tsv", sep="\t")
    post = df[(df["timing"] == "post") & (df["n_malignant_like"] >= 10)].copy()
    post["histology"] = post["Pathology"].map(map_histology)
    for hist, sub in post.groupby("histology"):
        for gene, col in (("TACSTD2", "mal_TACSTD2_mean"), ("CLDN4", "mal_CLDN4_mean")):
            add_gene_effects(
                rows,
                series="GSE207422",
                cohort="GSE207422_post_nmal10",
                histology=hist,
                gene=gene,
                compartment="malignant_like",
                x=sub[col],
                y=sub["frac_tnk"],
                min_epi=10,
                min_tnk=0,
                note="Hu 2023 neoadjuvant PD-1+chemo; post-surgery; marker malignant-like; author Pathology",
                primary=True,
            )


def gse205335(rows: list[dict]) -> None:
    df = pd.read_csv(HARVESTED / "gse205335_patient_table.tsv", sep="\t")
    df = df[df["n_malignant"] >= 20].copy()
    df["histology"] = df["cancer_subtype"].map(map_histology)
    for hist, sub in df.groupby("histology"):
        for gene, col in (("TACSTD2", "mal_TACSTD2_mean"), ("CLDN4", "mal_CLDN4_mean")):
            add_gene_effects(
                rows,
                series="GSE205335",
                cohort="GSE205335_nmal20",
                histology=hist,
                gene=gene,
                compartment="malignant",
                x=sub[col],
                y=sub["frac_tnk"],
                min_epi=20,
                min_tnk=0,
                note="Park/Ahn lung ICI atlas; author lineage.sub Malignant cells; patient unit; includes NE",
                primary=True,
            )


def gse131907(rows: list[dict]) -> None:
    df = pd.read_csv(HARVESTED / "gse131907_sample_table.tsv", sep="\t")
    # Kim 2020 is a LUAD atlas. origin tLung = primary tumor.
    specs = [
        ("tLung", "tLung", "epithelial", "epi_TACSTD2_mean", "epi_CLDN4_mean", 20, True, "primary tumor tLung; author Epithelial cells"),
        ("tumor_sites", None, "epithelial", "epi_TACSTD2_mean", "epi_CLDN4_mean", 20, False, "all tumor-bearing sites; author Epithelial cells"),
        ("tumor_sites_mal", None, "malignant", "mal_TACSTD2_mean", "mal_CLDN4_mean", 20, False, "author Malignant cells subtype; mostly mets/PE"),
    ]
    for cohort, origin, compartment, tcol, ccol, min_epi, primary, note in specs:
        if origin == "tLung":
            sub = df[(df["origin"] == "tLung") & (df["n_epithelial"] >= min_epi)].copy()
        elif cohort == "tumor_sites":
            sub = df[(df["origin"].isin(["tLung", "tL/B", "mLN", "mBrain", "PE"])) & (df["n_epithelial"] >= min_epi)].copy()
        else:
            sub = df[(df["n_malignant"] >= min_epi)].copy()
        for gene, col in (("TACSTD2", tcol), ("CLDN4", ccol)):
            add_gene_effects(
                rows,
                series="GSE131907",
                cohort=f"GSE131907_{cohort}",
                histology="LUAD",
                gene=gene,
                compartment=compartment,
                x=sub[col],
                y=sub["frac_tnk"],
                min_epi=min_epi,
                min_tnk=0,
                note="Kim 2020 LUAD atlas; no ICI labels; " + note,
                primary=primary,
            )


def gse253013(rows: list[dict]) -> None:
    df = pd.read_csv(HARVESTED / "gse253013_per_patient.tsv", sep="\t")
    sub = df[(df["tissue"] == "Tumor") & (df["eligible_malig"])].copy()
    for gene, col in (("TACSTD2", "TACSTD2_mean_log1p"), ("CLDN4", "CLDN4_mean_log1p")):
        add_gene_effects(
            rows,
            series="GSE253013",
            cohort="GSE253013_tumor_malig",
            histology="LUAD",
            gene=gene,
            compartment="malignant_like",
            x=sub[col],
            y=sub["tnk_fraction"],
            min_epi=10,
            min_tnk=20,
            note="Sze/Xiang 2024 treatment-naive LUAD; patient unit; no ICI labels",
            primary=True,
        )


def gse241934(rows: list[dict], path: Path) -> None:
    if not path.exists():
        return
    df = pd.read_csv(path, sep="\t")
    df["histology"] = df["histology_raw"].map(map_histology)
    keep = df[(df["n_epi"] >= 10) & (df["n_tnk"] >= 20)].copy()
    for cohort, sub0 in keep.groupby("cohort"):
        for hist, sub in sub0.groupby("histology"):
            for gene in ("TACSTD2", "CLDN4"):
                add_gene_effects(
                    rows,
                    series="GSE241934",
                    cohort=str(cohort),
                    histology=hist,
                    gene=gene,
                    compartment="epithelial",
                    x=sub[f"epi_{gene}_mean_log1p_cp10k"],
                    y=sub["frac_tnk"],
                    min_epi=10,
                    min_tnk=20,
                    note="Zhang/Zhong NEOTIDE+RWC; author Epi vs T+NK; IIT=EGFR-mut, RWC=WT LUAD/ASC",
                    primary=True,
                )
    # Combined LUAD across IIT+RWC as sensitivity (independent patients, same series)
    luad = keep[keep["histology"] == "LUAD"]
    for gene in ("TACSTD2", "CLDN4"):
        add_gene_effects(
            rows,
            series="GSE241934",
            cohort="GSE241934_LUAD_combined",
            histology="LUAD",
            gene=gene,
            compartment="epithelial",
            x=luad[f"epi_{gene}_mean_log1p_cp10k"],
            y=luad["frac_tnk"],
            min_epi=10,
            min_tnk=20,
            note="IIT+RWC LUAD patients combined; sensitivity only; not used in primary meta",
            primary=False,
        )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gse241934", default=str(HARVESTED / "gse241934_per_sample.tsv"))
    ap.add_argument("--out-dir", default=str(ROOT / "results"))
    args = ap.parse_args()
    rows: list[dict] = []
    gse207422(rows)
    gse205335(rows)
    gse131907(rows)
    gse253013(rows)
    gse241934(rows, Path(args.gse241934))
    out = pd.DataFrame(rows)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_dir / "cohort_effects.tsv", sep="\t", index=False)
    summary = {
        "n_effect_rows": int(len(out)),
        "primary_rows": int(out["primary_meta"].sum()),
        "histology_counts": out.groupby("histology").size().to_dict(),
        "series": sorted(out["series"].unique().tolist()),
    }
    (out_dir / "cohort_effects_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(out.to_string(index=False))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
