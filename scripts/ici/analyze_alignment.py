#!/usr/bin/env python3
"""Targeted follow-up: malignant-cell scoring and purity-adjusted bulk tests."""

from __future__ import annotations

import gzip
from pathlib import Path
import tempfile
import urllib.request

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from scipy.stats import mannwhitneyu, pearsonr, rankdata, spearmanr, t


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "ici"
SCRNA_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/"
    "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"
)
SCMETA_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/"
    "GSE207422_NSCLC_scRNAseq_metadata.xlsx"
)
PURITY_URL = (
    "https://api.gdc.cancer.gov/data/"
    "4f277128-f793-4354-a13d-30cc7fe9f6b5"
)
CBIO = "https://www.cbioportal.org/api"
GENES = {
    4070: "TACSTD2",
    1364: "CLDN4",
    915: "CD3D",
    916: "CD3E",
    925: "CD8A",
    4818: "NKG7",
    10578: "GNLY",
    3824: "KLRD1",
}
SC_GENES = {
    "TACSTD2",
    "CLDN4",
    "EPCAM",
    "KRT8",
    "KRT18",
    "KRT19",
    "KRT7",
    "PTPRC",
    "CD3D",
    "CD3E",
    "CD8A",
    "NKG7",
    "GNLY",
    "KLRD1",
    "HLA-A",
    "HLA-B",
    "HLA-C",
    "B2M",
    "STAT1",
    "IRF1",
    "IFITM1",
    "CXCL10",
}


def download(url: str, path: Path) -> Path:
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_suffix(path.suffix + ".part")
        urllib.request.urlretrieve(url, temporary)
        temporary.replace(path)
    return path


def extract_scrna(path: Path) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    selected: dict[str, np.ndarray] = {}
    with gzip.open(path, "rt") as handle:
        cells = np.asarray(handle.readline().rstrip("\n").split("\t")[1:])
        for line in handle:
            gene = line.split("\t", 1)[0]
            if gene in SC_GENES:
                selected[gene] = np.fromstring(
                    line[len(gene) + 1 :], sep="\t", dtype=np.int32
                )
    missing = SC_GENES.difference(selected)
    if missing:
        raise RuntimeError(f"Missing marker rows: {sorted(missing)}")
    return cells, selected


def scrna_patient_scores(
    matrix_path: Path, metadata_path: Path
) -> tuple[pd.DataFrame, list[dict]]:
    cells, counts = extract_scrna(matrix_path)
    samples = np.asarray(["_".join(cell.split("_")[:2]) for cell in cells])
    epithelial = sum(
        counts[gene] for gene in ("EPCAM", "KRT8", "KRT18", "KRT19", "KRT7")
    )
    # Operational tumor-epithelial gate. GEO supplies no cell-level malignant calls.
    malignant_gate = (epithelial > 0) & (counts["PTPRC"] == 0)
    tnk_gate = (counts["PTPRC"] > 0) & (
        sum(counts[gene] for gene in ("CD3D", "CD3E", "NKG7", "GNLY", "KLRD1"))
        > 0
    )
    mhc_ifn = ("HLA-A", "HLA-B", "HLA-C", "B2M", "STAT1", "IRF1", "IFITM1", "CXCL10")

    rows = []
    for sample in sorted(set(samples)):
        sample_cells = samples == sample
        malignant_cells = sample_cells & malignant_gate
        record = {
            "sample": sample,
            "n_cells": int(sample_cells.sum()),
            "n_epithelial_ptprc_negative": int(malignant_cells.sum()),
            "tnk_fraction_all_cells": float(tnk_gate[sample_cells].mean()),
        }
        for gene in ("TACSTD2", "CLDN4"):
            values = counts[gene][malignant_cells]
            record[f"{gene}_malignant_mean_log1p"] = float(np.log1p(values).mean())
            record[f"{gene}_malignant_detection_fraction"] = float((values > 0).mean())
            record[f"{gene}_malignant_umi_per_cell"] = float(values.mean())
        record["malignant_ifn_mhci_mean_log1p"] = float(
            np.mean([np.log1p(counts[gene][malignant_cells]).mean() for gene in mhc_ifn])
        )
        rows.append(record)

    scores = pd.DataFrame(rows)
    metadata = pd.read_excel(metadata_path)
    metadata = metadata[metadata["Sample"].astype(str).str.startswith("BD_")].rename(
        columns={
            "Sample": "sample",
            "Patient": "patient",
            "Resource": "resource",
            "Pathologic Response": "pathologic_response",
            "RECIST": "recist",
        }
    )
    scores = scores.merge(
        metadata[["sample", "patient", "resource", "pathologic_response", "recist"]],
        on="sample",
    )
    scores["mpr_or_pcr"] = scores["pathologic_response"].astype(str).str.match(
        r"^(MPR|pCR)"
    )
    post = scores[scores["resource"].astype(str).str.contains("Post")].copy()

    statistics: list[dict] = []
    for gene in ("TACSTD2", "CLDN4"):
        column = f"{gene}_malignant_mean_log1p"
        nmpr = post.loc[~post["mpr_or_pcr"], column]
        benefit = post.loc[post["mpr_or_pcr"], column]
        _, p_value = mannwhitneyu(nmpr, benefit, alternative="two-sided")
        statistics.append(
            {
                "dataset": "GSE207422 scRNA",
                "analysis": "post-treatment epithelial/PTPRC-negative gate; NMPR vs MPR/pCR",
                "feature": gene,
                "n": len(post),
                "n_group1": len(nmpr),
                "n_group2": len(benefit),
                "effect": float(nmpr.median() - benefit.median()),
                "effect_definition": "median(NMPR)-median(MPR/pCR), mean log1p UMI/cell",
                "rho": np.nan,
                "p": p_value,
            }
        )

    for feature, outcome in (
        ("TACSTD2_malignant_mean_log1p", "tnk_fraction_all_cells"),
        ("CLDN4_malignant_mean_log1p", "tnk_fraction_all_cells"),
        ("TACSTD2_malignant_mean_log1p", "malignant_ifn_mhci_mean_log1p"),
        ("CLDN4_malignant_mean_log1p", "malignant_ifn_mhci_mean_log1p"),
    ):
        rho, p_value = spearmanr(post[feature], post[outcome])
        statistics.append(
            {
                "dataset": "GSE207422 scRNA",
                "analysis": f"post-treatment patient-level {feature} vs {outcome}",
                "feature": feature.split("_")[0],
                "n": len(post),
                "n_group1": np.nan,
                "n_group2": np.nan,
                "effect": np.nan,
                "effect_definition": "Spearman rho",
                "rho": rho,
                "p": p_value,
            }
        )
    return scores, statistics


def partial_spearman(
    x: pd.Series, y: pd.Series, covariate: pd.Series
) -> tuple[float, float, np.ndarray, np.ndarray]:
    design = np.column_stack([np.ones(len(x)), rankdata(covariate)])
    x_rank = rankdata(x)
    y_rank = rankdata(y)
    x_residual = x_rank - design @ np.linalg.lstsq(design, x_rank, rcond=None)[0]
    y_residual = y_rank - design @ np.linalg.lstsq(design, y_rank, rcond=None)[0]
    rho = pearsonr(x_residual, y_residual).statistic
    degrees_freedom = len(x) - 3
    statistic = rho * np.sqrt(degrees_freedom / (1 - rho**2))
    p_value = 2 * t.sf(abs(statistic), degrees_freedom)
    return float(rho), float(p_value), x_residual, y_residual


def tcga_expression(study: str) -> pd.DataFrame:
    headers = {"Accept": "application/json", "Content-Type": "application/json"}
    samples = requests.get(
        f"{CBIO}/sample-lists/{study}_all", headers=headers, timeout=60
    ).json()["sampleIds"]
    response = requests.post(
        f"{CBIO}/molecular-profiles/{study}_rna_seq_v2_mrna/"
        "molecular-data/fetch?projection=DETAILED",
        headers=headers,
        json={"entrezGeneIds": list(GENES), "sampleIds": samples},
        timeout=120,
    )
    response.raise_for_status()
    records = [
        {
            "sample": item["sampleId"],
            "gene": item["gene"]["hugoGeneSymbol"],
            "value": item["value"],
        }
        for item in response.json()
    ]
    return pd.DataFrame(records).pivot(index="sample", columns="gene", values="value")


def tcga_analysis(purity_path: Path) -> tuple[pd.DataFrame, list[dict], dict]:
    purity = pd.read_csv(purity_path, sep="\t").set_index("array")["purity"]
    outputs = []
    statistics = []
    plot_data = {}
    immune_genes = ["CD3D", "CD3E", "CD8A", "NKG7", "GNLY", "KLRD1"]
    for histology in ("luad", "lusc"):
        study = f"{histology}_tcga_pan_can_atlas_2018"
        expression = tcga_expression(study).dropna()
        data = np.log2(expression + 1).join(purity).dropna()
        immune = data[immune_genes].apply(
            lambda column: (column - column.mean()) / column.std(ddof=1)
        ).mean(axis=1)
        data["tnk_score"] = immune
        data["histology"] = histology.upper()
        outputs.append(data.reset_index())
        for gene in ("TACSTD2", "CLDN4"):
            rho, p_value, x_residual, y_residual = partial_spearman(
                data[gene], immune, data["purity"]
            )
            statistics.append(
                {
                    "dataset": f"TCGA-{histology.upper()}",
                    "analysis": "partial Spearman; adjusted for ABSOLUTE purity",
                    "feature": gene,
                    "n": len(data),
                    "n_group1": np.nan,
                    "n_group2": np.nan,
                    "effect": np.nan,
                    "effect_definition": "partial Spearman rho",
                    "rho": rho,
                    "p": p_value,
                }
            )
            plot_data[(histology, gene)] = (x_residual, y_residual, rho, p_value)
        rho, p_value, _, _ = partial_spearman(
            data["TACSTD2"], data["CLDN4"], data["purity"]
        )
        statistics.append(
            {
                "dataset": f"TCGA-{histology.upper()}",
                "analysis": "partial Spearman; adjusted for ABSOLUTE purity",
                "feature": "TACSTD2~CLDN4",
                "n": len(data),
                "n_group1": np.nan,
                "n_group2": np.nan,
                "effect": np.nan,
                "effect_definition": "partial Spearman rho",
                "rho": rho,
                "p": p_value,
            }
        )
    return pd.concat(outputs, ignore_index=True), statistics, plot_data


def oncosg_analysis() -> tuple[pd.DataFrame, list[dict]]:
    study = "luad_oncosg_2020"
    profile = f"{study}_rna_seq_v2_mrna_median_all_sample_Zscores"
    headers = {"Accept": "application/json", "Content-Type": "application/json"}
    samples = requests.get(
        f"{CBIO}/sample-lists/{study}_rna_seq_v2_mrna",
        headers=headers,
        timeout=60,
    ).json()["sampleIds"]
    response = requests.post(
        f"{CBIO}/molecular-profiles/{profile}/molecular-data/fetch?projection=DETAILED",
        headers=headers,
        json={"entrezGeneIds": list(GENES), "sampleIds": samples},
        timeout=120,
    )
    response.raise_for_status()
    expression = pd.DataFrame(
        [
            {
                "sample": item["sampleId"],
                "gene": item["gene"]["hugoGeneSymbol"],
                "value": item["value"],
            }
            for item in response.json()
        ]
    ).pivot(index="sample", columns="gene", values="value").dropna()
    immune = expression[["CD3D", "CD3E", "CD8A", "NKG7", "GNLY", "KLRD1"]].mean(
        axis=1
    )
    statistics = []
    for feature, outcome, label in (
        ("TACSTD2", immune, "T/NK z-score"),
        ("CLDN4", immune, "T/NK z-score"),
        ("TACSTD2", expression["CLDN4"], "CLDN4"),
    ):
        rho, p_value = spearmanr(expression[feature], outcome)
        statistics.append(
            {
                "dataset": "OncoSG-LUAD",
                "analysis": f"unadjusted Spearman; {feature} vs {label}",
                "feature": feature if label != "CLDN4" else "TACSTD2~CLDN4",
                "n": len(expression),
                "n_group1": np.nan,
                "n_group2": np.nan,
                "effect": np.nan,
                "effect_definition": "Spearman rho",
                "rho": rho,
                "p": p_value,
            }
        )
    return expression.reset_index(), statistics


def bh_adjust(values: pd.Series) -> pd.Series:
    p = values.astype(float).to_numpy()
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    adjusted = np.empty(len(p))
    adjusted[order] = np.minimum.accumulate(ranked[::-1])[::-1]
    return pd.Series(np.minimum(adjusted, 1), index=values.index)


def plot_alignment(
    scores: pd.DataFrame, tcga_plot: dict, statistics: pd.DataFrame
) -> None:
    post = scores[scores["resource"].astype(str).str.contains("Post")]
    fig, axes = plt.subplots(2, 2, figsize=(9, 8))
    rng = np.random.default_rng(20260816)
    for ax, gene in zip(axes[0], ("TACSTD2", "CLDN4")):
        column = f"{gene}_malignant_mean_log1p"
        groups = [
            post.loc[~post["mpr_or_pcr"], column],
            post.loc[post["mpr_or_pcr"], column],
        ]
        ax.boxplot(groups, tick_labels=[f"NMPR (n={len(groups[0])})", f"MPR/pCR (n={len(groups[1])})"])
        for x_position, values in enumerate(groups, 1):
            ax.scatter(
                x_position + rng.uniform(-0.06, 0.06, len(values)),
                values,
                color="#276FBF",
            )
        row = statistics[
            (statistics["dataset"] == "GSE207422 scRNA")
            & statistics["analysis"].str.contains("NMPR vs")
            & (statistics["feature"] == gene)
        ].iloc[0]
        ax.set(
            title=f"GSE207422 epithelial/PTPRC−: {gene}\nNMPR−MPR={row.effect:.3f}, p={row.p:.3g}",
            ylabel="Mean log1p UMI per gated cell",
        )
    for ax, histology in zip(axes[1], ("luad", "lusc")):
        x, y, rho, p_value = tcga_plot[(histology, "TACSTD2")]
        ax.scatter(x, y, s=12, alpha=0.5, color="#9C2C77")
        ax.set(
            title=f"TCGA-{histology.upper()} TACSTD2 vs T/NK\npurity-adjusted ρ={rho:.3f}, p={p_value:.3g}",
            xlabel="TACSTD2 rank residual",
            ylabel="T/NK-score rank residual",
        )
    fig.tight_layout()
    fig.savefig(OUT / "user_alignment.png", dpi=180)
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    cache = Path(tempfile.gettempdir()) / "ici_alignment"
    matrix = download(SCRNA_URL, cache / SCRNA_URL.rsplit("/", 1)[-1])
    metadata = download(SCMETA_URL, cache / SCMETA_URL.rsplit("/", 1)[-1])
    purity = download(PURITY_URL, cache / "TCGA_mastercalls.abs_tables_JSedit.fixed.txt")

    scores, scrna_statistics = scrna_patient_scores(matrix, metadata)
    scores.to_csv(OUT / "gse207422_scrna_patient_scores.tsv", sep="\t", index=False)
    tcga, tcga_statistics, plot_data = tcga_analysis(purity)
    tcga.to_csv(OUT / "tcga_lung_selected_expression.tsv", sep="\t", index=False)
    oncosg, oncosg_statistics = oncosg_analysis()
    oncosg.to_csv(OUT / "oncosg_luad_selected_expression.tsv", sep="\t", index=False)
    statistics = pd.DataFrame(
        scrna_statistics + tcga_statistics + oncosg_statistics
    )
    statistics["q_bh"] = bh_adjust(statistics["p"])
    statistics.to_csv(OUT / "alignment_statistics.tsv", sep="\t", index=False)
    plot_alignment(scores, plot_data, statistics)


if __name__ == "__main__":
    main()
