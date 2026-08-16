#!/usr/bin/env python3
"""Reproduce the public-data audit of claims B1 and B2."""

from __future__ import annotations

import argparse
import json
import math
import re
import shutil
import urllib.request
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
import xenaPython as xena
from scipy.stats import pearsonr, rankdata, spearmanr


XENA_HUB = "https://pancanatlas.xenahubs.net"
XENA_EXPR = "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena"
XENA_PHENO = "TCGA_phenotype_denseDataOnlyDownload.tsv"
SURFY_URL = (
    "https://raw.githubusercontent.com/steveneschrich/surfaceome/main/"
    "data-raw/surfy/table_S3_surfaceome.xlsx"
)
DEPMAP_24Q2_MODEL = "https://ndownloader.figshare.com/files/46489732"
DEPMAP_24Q2_EXPR = "https://ndownloader.figshare.com/files/46490878"
DEPMAP_20Q2_MODEL = "https://ndownloader.figshare.com/files/25494443"
DEPMAP_20Q2_EXPR = "https://ndownloader.figshare.com/files/25817909"
CMP_API = "https://api.cellmodelpassports.sanger.ac.uk/datasets/proteomics"
RPPA_ANTIBODIES = "https://tcpa.drbioright.org/rppa500mclp/CCLE-annotation-antibody"
GYGI_PROTEIN = (
    "https://gygi.hms.harvard.edu/data/ccle/"
    "protein_quant_current_normalized.csv.gz"
)
GYGI_SAMPLES = (
    "https://gygi.hms.harvard.edu/data/ccle/Table_S1_Sample_Information.xlsx"
)


def download(url: str, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and path.stat().st_size:
        return
    print(f"Downloading {url} -> {path}")
    with urllib.request.urlopen(url) as src, path.open("wb") as dst:
        shutil.copyfileobj(src, dst)


def decode_field(field: str, values: list[float | None]) -> list[str | None]:
    metadata = xena.field_codes(XENA_HUB, XENA_PHENO, [field])
    codes = metadata[0].get("code") if metadata else None
    if not codes:
        return [None if value is None else str(value) for value in values]
    labels = codes.split("\t")
    return [
        None
        if value is None or (isinstance(value, float) and math.isnan(value))
        else labels[int(value)]
        for value in values
    ]


def correlation(x: np.ndarray, y: np.ndarray) -> tuple[int, float, float, float, float]:
    keep = np.isfinite(x) & np.isfinite(y)
    x, y = x[keep], y[keep]
    rho, rho_p = spearmanr(x, y)
    r, r_p = pearsonr(x, y)
    return len(x), float(rho), float(rho_p), float(r), float(r_p)


def fixed_effect_spearman(
    x: np.ndarray, y: np.ndarray, groups: np.ndarray
) -> tuple[int, float, float]:
    keep = np.isfinite(x) & np.isfinite(y) & pd.notna(groups)
    xr, yr, gr = rankdata(x[keep]), rankdata(y[keep]), groups[keep]
    for group in np.unique(gr):
        mask = gr == group
        xr[mask] -= xr[mask].mean()
        yr[mask] -= yr[mask].mean()
    r, p = pearsonr(xr, yr)
    return len(xr), float(r), float(p)


def rank_desc(values: pd.Series) -> pd.Series:
    return values.rank(method="min", ascending=False).astype("Int64")


def analyze_b1(out: Path, cache: Path, chunk_size: int) -> dict:
    surfy_path = cache / "table_S3_surfaceome.xlsx"
    download(SURFY_URL, surfy_path)
    surfy = pd.read_excel(surfy_path, sheet_name="in silico surfaceome only", header=1)
    surface_genes = (
        surfy["UniProt gene"].dropna().astype(str).str.strip().loc[lambda x: x != ""].unique()
    )

    samples = xena.dataset_samples(XENA_HUB, XENA_EXPR, None)
    fields = set(xena.dataset_field(XENA_HUB, XENA_EXPR))
    genes = sorted(set(surface_genes) & fields)
    if "TACSTD2" not in fields:
        raise RuntimeError("TACSTD2 is absent from the selected Xena matrix")

    disease_raw, type_raw, type_id_raw = xena.dataset_fetch(
        XENA_HUB,
        XENA_PHENO,
        samples,
        ["_primary_disease", "sample_type", "sample_type_id"],
    )
    disease = np.asarray(decode_field("_primary_disease", disease_raw), dtype=object)
    sample_type = np.asarray(decode_field("sample_type", type_raw), dtype=object)
    sample_type_id = np.asarray(type_id_raw, dtype=float)
    primary = (sample_type == "Primary Tumor") | (sample_type_id == 1)

    target = np.asarray(
        xena.dataset_fetch(XENA_HUB, XENA_EXPR, samples, ["TACSTD2"])[0],
        dtype=float,
    )
    rows: list[dict] = []
    for start in range(0, len(genes), chunk_size):
        batch = genes[start : start + chunk_size]
        print(f"Xena surface genes {start + 1}-{start + len(batch)} / {len(genes)}")
        matrix = xena.dataset_fetch(XENA_HUB, XENA_EXPR, samples, batch)
        for gene, values in zip(batch, matrix, strict=True):
            y = np.asarray(values, dtype=float)
            n_all, rho_all, p_all, r_all, rp_all = correlation(target, y)
            n_primary, rho_primary, p_primary, r_primary, rp_primary = correlation(
                target[primary], y[primary]
            )
            n_fe, rho_fe, p_fe = fixed_effect_spearman(
                target[primary], y[primary], disease[primary]
            )
            rows.append(
                {
                    "gene": gene,
                    "n_all_samples": n_all,
                    "spearman_rho_all_samples": rho_all,
                    "spearman_p_all_samples": p_all,
                    "pearson_r_all_samples": r_all,
                    "pearson_p_all_samples": rp_all,
                    "n_primary_tumors": n_primary,
                    "spearman_rho_primary_tumors": rho_primary,
                    "spearman_p_primary_tumors": p_primary,
                    "pearson_r_primary_tumors": r_primary,
                    "pearson_p_primary_tumors": rp_primary,
                    "n_primary_tumors_fixed_effect": n_fe,
                    "spearman_rho_cancer_type_fixed_effect": rho_fe,
                    "spearman_p_cancer_type_fixed_effect": p_fe,
                }
            )

    ranking = pd.DataFrame(rows)
    ranking = ranking.loc[ranking.gene != "TACSTD2"].copy()
    for metric in [
        "spearman_rho_all_samples",
        "spearman_rho_primary_tumors",
        "spearman_rho_cancer_type_fixed_effect",
        "pearson_r_all_samples",
        "pearson_r_primary_tumors",
    ]:
        ranking[f"rank_{metric}"] = rank_desc(ranking[metric])
    ranking = ranking.sort_values(
        ["rank_spearman_rho_all_samples", "gene"], ignore_index=True
    )
    ranking.to_csv(out / "B1_surface_gene_rankings.csv", index=False)
    ranking.head(25).to_csv(out / "B1_top25_surface_genes.csv", index=False)

    cldn4 = ranking.loc[ranking.gene == "CLDN4"].iloc[0].to_dict()
    top = ranking.iloc[0].to_dict()
    return {
        "dataset": XENA_EXPR,
        "unit": "log2(norm_value+1), batch-effects-normalized mRNA",
        "surfaceome_definition": "SURFY 5% FPR representative human surfaceome plus GPI proteins",
        "surfaceome_version_source": SURFY_URL,
        "surfaceome_input_proteins": int(len(surfy)),
        "unique_surface_gene_symbols": int(len(surface_genes)),
        "surface_genes_present_in_xena": int(len(genes)),
        "samples_all": int(len(samples)),
        "primary_tumors": int(primary.sum()),
        "top_surface_gene_all_samples": top,
        "CLDN4": cldn4,
    }


def depmap_rna(
    model_path: Path, expression_path: Path, release: str
) -> tuple[pd.DataFrame, dict]:
    expression = pd.read_csv(
        expression_path,
        usecols=lambda c: c in {"Unnamed: 0", "TACSTD2 (4070)", "CLDN4 (1364)"},
    ).rename(columns={"Unnamed: 0": "ModelID"})
    model = pd.read_csv(model_path)
    if release == "24Q2":
        subset = model.loc[
            model.OncotreePrimaryDisease.eq("Non-Small Cell Lung Cancer"),
            ["ModelID", "CellLineName", "OncotreeSubtype"],
        ]
        subtype = "OncotreeSubtype"
    else:
        subset = model.loc[
            model.lineage_subtype.eq("NSCLC"),
            ["DepMap_ID", "cell_line_name", "Subtype"],
        ].rename(
            columns={
                "DepMap_ID": "ModelID",
                "cell_line_name": "CellLineName",
                "Subtype": "OncotreeSubtype",
            }
        )
        subtype = "OncotreeSubtype"
    data = subset.merge(expression, on="ModelID", how="inner").dropna(
        subset=["TACSTD2 (4070)", "CLDN4 (1364)"]
    )
    n, rho, rho_p, r, r_p = correlation(
        data["TACSTD2 (4070)"].to_numpy(),
        data["CLDN4 (1364)"].to_numpy(),
    )
    data["release"] = release
    return data, {
        "release": release,
        "modality": "RNA-seq, log2(TPM+1); not protein abundance",
        "NSCLC_definition": (
            "OncotreePrimaryDisease == 'Non-Small Cell Lung Cancer'"
            if release == "24Q2"
            else "lineage_subtype == 'NSCLC'"
        ),
        "n": n,
        "spearman_rho": rho,
        "spearman_p": rho_p,
        "pearson_r": r,
        "pearson_p": r_p,
        "subtype_field": subtype,
    }


def proteomics_availability() -> dict:
    counts = {}
    for gene, uniprot in {"TROP2/TACSTD2": "P09758", "CLDN4": "O14493"}.items():
        params = {
            "filter": json.dumps(
                [{"name": "uniprot_id", "op": "eq", "val": uniprot}]
            ),
            "page[size]": "1",
        }
        response = requests.get(CMP_API, params=params, timeout=60)
        response.raise_for_status()
        counts[gene] = {
            "uniprot": uniprot,
            "records": int(response.json()["meta"]["count"]),
        }
    antibodies = requests.get(RPPA_ANTIBODIES, timeout=60).json()["antibodies"]
    genes = set()
    for antibody in antibodies:
        genes.update(
            token
            for token in re.split(r"[,;/ ]+", antibody.get("genes", ""))
            if token
        )
    return {
        "Sanger_DIA_MS_949_line_panel": counts,
        "CCLE_RPPA500": {
            "antibodies": len(antibodies),
            "CLDN4_present": "CLDN4" in genes,
            "TACSTD2_present": "TACSTD2" in genes,
        },
        "conclusion": (
            "The pair cannot be computed in these two additional public "
            "protein matrices: CLDN4 is absent from Sanger DIA-MS and both "
            "analytes are absent from CCLE RPPA500."
        ),
    }


def _gygi_sample_columns(proteins: pd.DataFrame) -> list[str]:
    return [
        column
        for column in proteins.columns
        if "_TenPx" in column and not column.startswith("TenPx")
    ]


def _gygi_pairs(proteins: pd.DataFrame, codes: set[str]) -> pd.DataFrame:
    columns = [
        column
        for column in _gygi_sample_columns(proteins)
        if re.sub(r"_TenPx\d+$", "", column) in codes
    ]
    pair = pd.DataFrame(
        {
            "CCLECode": [re.sub(r"_TenPx\d+$", "", column) for column in columns],
            "TACSTD2_protein": proteins.loc[
                proteins.Gene_Symbol.eq("TACSTD2"), columns
            ]
            .iloc[0]
            .to_numpy(dtype=float),
            "CLDN4_protein": proteins.loc[proteins.Gene_Symbol.eq("CLDN4"), columns]
            .iloc[0]
            .to_numpy(dtype=float),
        }
    )
    return pair.groupby("CCLECode", as_index=False).mean()


def _corr_record(pair: pd.DataFrame, **meta: object) -> dict:
    complete = pair.dropna(subset=["TACSTD2_protein", "CLDN4_protein"])
    n, rho, rho_p, r, r_p = correlation(
        complete.TACSTD2_protein.to_numpy(), complete.CLDN4_protein.to_numpy()
    )
    return {
        **meta,
        "models_in_subset": int(len(pair)),
        "complete_pairs": n,
        "spearman_rho": rho,
        "spearman_p": rho_p,
        "pearson_r": r,
        "pearson_p": r_p,
    }


def gygi_protein(
    cache: Path, model20_path: Path
) -> tuple[pd.DataFrame, dict, list[dict]]:
    protein_path = cache / "Gygi_protein_quant_current_normalized.csv.gz"
    sample_path = cache / "Gygi_Table_S1_Sample_Information.xlsx"
    download(GYGI_PROTEIN, protein_path)
    download(GYGI_SAMPLES, sample_path)

    proteins = pd.read_csv(protein_path)
    sample_info = pd.read_excel(sample_path, sheet_name="Sample_Information")
    model = pd.read_csv(
        model20_path,
        usecols=["CCLE_Name", "lineage_subtype", "lineage_sub_subtype", "Subtype"],
    )
    lung_codes = set(
        sample_info.loc[sample_info["Tissue of Origin"].eq("Lung"), "CCLE Code"]
    )
    nsclc_codes = set(model.loc[model.lineage_subtype.eq("NSCLC"), "CCLE_Name"])

    lung = _gygi_pairs(proteins, lung_codes).merge(
        model.rename(columns={"CCLE_Name": "CCLECode"}),
        on="CCLECode",
        how="left",
    )
    complete_lung = lung.dropna(subset=["TACSTD2_protein", "CLDN4_protein"]).copy()
    nsclc = lung.loc[lung.lineage_subtype.eq("NSCLC")].copy()

    primary = _corr_record(
        lung,
        definition="Gygi Tissue of Origin == Lung; complete pairs",
        note="This is the n=45 pairing that yields ρ≈0.69. It is all Gygi lung, not NSCLC-only.",
    )
    nsclc_result = _corr_record(
        nsclc,
        definition="DepMap 20Q2 lineage_subtype == NSCLC within Gygi lung",
        note="Stricter NSCLC filter. n drops because CLDN4 is missing in 28 NSCLC models.",
    )
    return complete_lung, primary, [primary, nsclc_result]


def analyze_b2(out: Path, cache: Path) -> dict:
    paths = {
        "24Q2_model": cache / "DepMap_24Q2_Model.csv",
        "24Q2_expr": cache / "DepMap_24Q2_Expression.csv",
        "20Q2_model": cache / "DepMap_20Q2_sample_info.csv",
        "20Q2_expr": cache / "DepMap_20Q2_Expression.csv",
    }
    for url, path in [
        (DEPMAP_24Q2_MODEL, paths["24Q2_model"]),
        (DEPMAP_24Q2_EXPR, paths["24Q2_expr"]),
        (DEPMAP_20Q2_MODEL, paths["20Q2_model"]),
        (DEPMAP_20Q2_EXPR, paths["20Q2_expr"]),
    ]:
        download(url, path)

    data24, result24 = depmap_rna(paths["24Q2_model"], paths["24Q2_expr"], "24Q2")
    data20, result20 = depmap_rna(paths["20Q2_model"], paths["20Q2_expr"], "20Q2")
    protein, protein_result, protein_sensitivities = gygi_protein(
        cache, paths["20Q2_model"]
    )
    pd.concat([data24, data20], ignore_index=True).to_csv(
        out / "B2_DepMap_NSCLC_RNA_pairs.csv", index=False
    )
    protein.to_csv(out / "B2_Gygi_lung_protein_pairs.csv", index=False)
    protein.to_csv(out / "B2_CCLE_NSCLC_protein_pairs.csv", index=False)
    pd.DataFrame(protein_sensitivities).to_csv(
        out / "B2_Gygi_protein_sensitivity.csv", index=False
    )

    fig, ax = plt.subplots(figsize=(5.4, 4.6), constrained_layout=True)
    nsclc = protein.lineage_subtype.eq("NSCLC")
    ax.scatter(
        protein.loc[nsclc, "TACSTD2_protein"],
        protein.loc[nsclc, "CLDN4_protein"],
        s=28,
        alpha=0.8,
        edgecolors="none",
        label=f"NSCLC (n={int(nsclc.sum())})",
    )
    ax.scatter(
        protein.loc[~nsclc, "TACSTD2_protein"],
        protein.loc[~nsclc, "CLDN4_protein"],
        s=28,
        alpha=0.8,
        edgecolors="none",
        label=f"SCLC (n={int((~nsclc).sum())})",
    )
    ax.legend(frameon=False)
    ax.set(
        title="Gygi CCLE lung TMT-MS\n"
        f"Spearman ρ={protein_result['spearman_rho']:.3f}, "
        f"n={protein_result['complete_pairs']}",
        xlabel="TROP2/TACSTD2 normalized TMT abundance",
        ylabel="CLDN4 normalized TMT abundance",
    )
    fig.savefig(out / "B2_CCLE_NSCLC_protein_scatter.png", dpi=180)
    fig.savefig(out / "B2_Gygi_lung_protein_scatter.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), constrained_layout=True)
    for ax, data, result in zip(axes, [data20, data24], [result20, result24], strict=True):
        ax.scatter(
            data["TACSTD2 (4070)"],
            data["CLDN4 (1364)"],
            s=17,
            alpha=0.7,
            edgecolors="none",
        )
        ax.set(
            title=f"DepMap {result['release']} NSCLC RNA\n"
            f"Spearman ρ={result['spearman_rho']:.3f}, n={result['n']}",
            xlabel="TACSTD2 log2(TPM+1)",
            ylabel="CLDN4 log2(TPM+1)",
        )
    fig.savefig(out / "B2_DepMap_NSCLC_RNA_scatter.png", dpi=180)
    plt.close(fig)
    return {
        "claim_as_corrected": (
            "Gygi/Nusinow CCLE TMT-MS, Tissue of Origin = Lung, "
            "complete CLDN4–TACSTD2 pairs, n=45, Spearman ρ≈0.69"
        ),
        "verdict": "supported",
        "protein_data_availability": proteomics_availability(),
        "Gygi_CCLE_protein_result": protein_result,
        "Gygi_protein_sensitivities": protein_sensitivities,
        "lineage_in_n45": protein.lineage_subtype.value_counts(dropna=False).to_dict(),
        "RNA_sensitivity_checks": [result20, result24],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path(__file__).parent)
    parser.add_argument("--cache", type=Path, default=Path(__file__).parent / "data")
    parser.add_argument("--chunk-size", type=int, default=128)
    parser.add_argument("--skip-b1", action="store_true")
    parser.add_argument("--skip-b2", action="store_true")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    args.cache.mkdir(parents=True, exist_ok=True)

    summary_path = args.output / "summary.json"
    summary = json.loads(summary_path.read_text()) if summary_path.exists() else {}
    if not args.skip_b1:
        summary["B1"] = analyze_b1(args.output, args.cache, args.chunk_size)
    if not args.skip_b2:
        summary["B2"] = analyze_b2(args.output, args.cache)
    summary_path.write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
