#!/usr/bin/env python3
"""Reproduce the public-data audit of the CLDN4-loss IFN/MHC-I claim."""

from __future__ import annotations

import gzip
import hashlib
import io
import math
import urllib.parse
import urllib.request
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parent
RAW = ROOT / "data" / "raw"
PROCESSED = ROOT / "data" / "processed"
FIGURES = ROOT / "figures"

HUMAN_SIGNATURE = ["IFI27", "OAS2", "IFIT1", "MX1", "ISG15", "HLA-A"]
MOUSE_MAP = {
    "IFI27": "Ifi27l2a",
    "OAS2": "Oas2",
    "IFIT1": "Ifit1",
    "MX1": "Mx1",
    "ISG15": "Isg15",
    "HLA-A": "H2-K1",
}

URLS = {
    "GSE22493_family.soft.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE22nnn/GSE22493/soft/GSE22493_family.soft.gz",
    "GSE22493_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE22nnn/GSE22493/matrix/GSE22493_series_matrix.txt.gz",
    "GSE50927_Cldn4lungWTvsKOgenes.csv.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50927/suppl/GSE50927_Cldn4lungWTvsKOgenes.csv.gz",
    "GSE50927_VILIwtkoloGenes.csv.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50927/suppl/GSE50927_VILIwtkoloGenes.csv.gz",
    "GSE50927_VILIwtkohiGenes.csv.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50927/suppl/GSE50927_VILIwtkohiGenes.csv.gz",
    "GSE207704_CLDN4_RNAseq.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207704/suppl/GSE207704_CLDN4_RNAseq.txt.gz",
    "LINCS_CRISPR_KO_consensus.gmt": "https://cfde-drc.s3.amazonaws.com/LINCS/XMT/2022-12-13/LINCS_XMT_2022-12-13_LINCS_L1000_CRISPR_KO_Consensus_Sigs.gmt",
}


def download() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    PROCESSED.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    records = []
    for name, url in URLS.items():
        path = RAW / name
        if not path.exists():
            urllib.request.urlretrieve(url, path)
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        records.append(
            {"file": name, "url": url, "bytes": path.stat().st_size, "sha256": digest}
        )
    pd.DataFrame(records).to_csv(PROCESSED / "source_manifest.tsv", sep="\t", index=False)


def soft_table(path: Path, begin: str, end: str) -> pd.DataFrame:
    with gzip.open(path, "rt", errors="replace") as handle:
        lines = handle.readlines()
    start = lines.index(begin + "\n") + 1
    stop = lines.index(end + "\n")
    return pd.read_csv(io.StringIO("".join(lines[start:stop])), sep="\t", dtype=str)


def array_effects() -> list[dict]:
    """Read submitted normalized log2(channel 2 KD / channel 1 control) ratios."""
    platform = soft_table(
        RAW / "GSE22493_family.soft.gz", "!platform_table_begin", "!platform_table_end"
    )
    matrix = soft_table(
        RAW / "GSE22493_series_matrix.txt.gz",
        "!series_matrix_table_begin",
        "!series_matrix_table_end",
    )
    matrix = matrix.rename(columns={"ID_REF": "ID"})
    sample_cols = ["GSM558700", "GSM558701", "GSM558702"]
    for col in sample_cols:
        matrix[col] = pd.to_numeric(matrix[col], errors="coerce")

    # ORF is incomplete on this old array. Recover unambiguous SYMBOL-- descriptions.
    platform["symbol"] = platform["ORF"].fillna("")
    desc_symbol = platform["DESCRIPTION"].fillna("").str.extract(
        r"^([A-Za-z0-9.-]+)--", expand=False
    )
    platform.loc[platform["symbol"].eq(""), "symbol"] = desc_symbol
    platform["symbol"] = platform["symbol"].replace({"G1P2": "ISG15"})
    joined = matrix.merge(platform[["ID", "symbol"]], on="ID", how="left")

    rows = []
    for gene in ["CLDN4"] + HUMAN_SIGNATURE:
        probes = joined.loc[joined["symbol"].eq(gene), sample_cols]
        if probes.empty:
            rows.append(effect_row("GSE22493", "SKOV-3 KD", gene, None, note="not on array"))
            continue
        # Multiple printed probes are collapsed within each paired array.
        values = probes.median(axis=0, skipna=True).dropna().to_numpy(dtype=float)
        if not len(values):
            rows.append(effect_row("GSE22493", "SKOV-3 KD", gene, None, note="all missing"))
            continue
        pvalue = (
            stats.ttest_1samp(values, 0, alternative="greater").pvalue
            if len(values) >= 2
            else np.nan
        )
        rows.append(
            effect_row(
                "GSE22493",
                "SKOV-3 KD",
                gene,
                float(np.mean(values)),
                pvalue=float(pvalue),
                n=len(values),
                note=f"{len(probes)} probe(s); paired two-colour arrays",
            )
        )
    return rows


def mouse_effects() -> list[dict]:
    contrasts = {
        "baseline KO vs WT": "GSE50927_Cldn4lungWTvsKOgenes.csv.gz",
        "VILI-low KO vs WT VILI": "GSE50927_VILIwtkoloGenes.csv.gz",
        "VILI-high KO vs WT VILI": "GSE50927_VILIwtkohiGenes.csv.gz",
    }
    rows = []
    for contrast, filename in contrasts.items():
        table = pd.read_csv(RAW / filename)
        for human_gene, mouse_gene in {"CLDN4": "Cldn4", **MOUSE_MAP}.items():
            hit = table.loc[table["Marker.Symbol"].eq(mouse_gene)]
            if hit.empty:
                rows.append(
                    effect_row(
                        "GSE50927", contrast, human_gene, None, ortholog=mouse_gene
                    )
                )
                continue
            record = hit.iloc[0]
            rows.append(
                effect_row(
                    "GSE50927",
                    contrast,
                    human_gene,
                    float(record.logFC),
                    pvalue=float(record.PValue),
                    padj=float(record.FDR),
                    ortholog=mouse_gene,
                    n=2,
                    note="submitter edgeR result; whole mouse lung",
                )
            )
    return rows


def breast_effects() -> list[dict]:
    """Use the only submitted processed values (condition-mean FPKM, no replicate values)."""
    table = pd.read_csv(RAW / "GSE207704_CLDN4_RNAseq.txt.gz", sep="\t")
    rows = []
    columns = {
        "MCF-7 KO vs WT": ("MCF7_CLDN4KO_FPKM (fpkm)", "MCF7_WT_FPKM (fpkm)"),
        "T47D KO vs WT": ("T47D_CLDN4KO_FPKM (fpkm)", "T47D_WT_FPKM (fpkm)"),
    }
    for contrast, (ko_col, wt_col) in columns.items():
        for gene in ["CLDN4"] + HUMAN_SIGNATURE:
            hits = table.loc[table["gene_short_name"].eq(gene), [ko_col, wt_col]]
            if hits.empty:
                rows.append(
                    effect_row(
                        "GSE207704",
                        contrast,
                        gene,
                        None,
                        note="omitted from submitted processed table",
                    )
                )
                continue
            # Cufflinks contains two CLDN4 loci; sum FPKM before computing the ratio.
            ko, wt = hits.sum(axis=0)
            logfc = math.log2((float(ko) + 0.5) / (float(wt) + 0.5))
            rows.append(
                effect_row(
                    "GSE207704",
                    contrast,
                    gene,
                    logfc,
                    n=2,
                    note="log2 ratio of submitted mean FPKM; replicate values unavailable",
                )
            )
    return rows


def effect_row(
    dataset: str,
    contrast: str,
    gene: str,
    log2fc: float | None,
    *,
    pvalue: float | None = None,
    padj: float | None = None,
    ortholog: str = "",
    n: int | None = None,
    note: str = "",
) -> dict:
    return {
        "dataset": dataset,
        "contrast": contrast,
        "gene": gene,
        "assayed_symbol": ortholog or gene,
        "log2FC_KD_or_KO_vs_control": log2fc,
        "pvalue": pvalue,
        "FDR": padj,
        "nominal_n_per_condition": n,
        "note": note,
    }


def adjust_array_pvalues(effects: pd.DataFrame) -> pd.DataFrame:
    mask = effects["dataset"].eq("GSE22493") & effects["pvalue"].notna()
    if mask.any():
        effects.loc[mask, "FDR"] = multipletests(
            effects.loc[mask, "pvalue"], method="fdr_bh"
        )[1]
    return effects


def summarize(effects: pd.DataFrame) -> pd.DataFrame:
    rows = []
    signature = effects.loc[effects["gene"].isin(HUMAN_SIGNATURE)].copy()
    for (dataset, contrast), group in signature.groupby(["dataset", "contrast"], sort=False):
        observed = group.dropna(subset=["log2FC_KD_or_KO_vs_control"])
        values = observed["log2FC_KD_or_KO_vs_control"].to_numpy()
        positives = int((values > 0).sum())
        count = len(values)
        sign_p = (
            float(stats.binomtest(positives, count, 0.5, alternative="greater").pvalue)
            if count
            else np.nan
        )
        cldn = effects.loc[
            effects["dataset"].eq(dataset)
            & effects["contrast"].eq(contrast)
            & effects["gene"].eq("CLDN4"),
            "log2FC_KD_or_KO_vs_control",
        ]
        rows.append(
            {
                "dataset": dataset,
                "contrast": contrast,
                "CLDN4_log2FC": cldn.iloc[0] if len(cldn) else np.nan,
                "signature_genes_observed_of_6": count,
                "signature_genes_positive": positives,
                "median_signature_log2FC": np.median(values) if count else np.nan,
                "one_sided_exact_sign_p": sign_p,
                "strict_all_observed_positive": bool(count == 6 and positives == 6),
            }
        )
    return pd.DataFrame(rows)


def check_lincs() -> pd.DataFrame:
    matches = []
    with open(RAW / "LINCS_CRISPR_KO_consensus.gmt", encoding="utf-8") as handle:
        for line in handle:
            name = line.split("\t", 1)[0]
            if name.upper() in {"CLDN4 UP", "CLDN4 DOWN"}:
                matches.append(name)
    return pd.DataFrame(
        [
            {
                "resource": "LINCS L1000 CRISPR KO consensus",
                "release": "CFDE 2022-12-13 / Harmonizome 2023-09-05",
                "CLDN4_signature_count": len(matches) // 2,
                "result": "No CLDN4 perturbation signature" if not matches else "; ".join(matches),
            }
        ]
    )


def plot_heatmap(effects: pd.DataFrame) -> None:
    signature = effects.loc[effects["gene"].isin(HUMAN_SIGNATURE)].copy()
    signature["label"] = signature["dataset"] + "\n" + signature["contrast"]
    labels = signature["label"].drop_duplicates().tolist()
    matrix = (
        signature.pivot(index="label", columns="gene", values="log2FC_KD_or_KO_vs_control")
        .reindex(index=labels, columns=HUMAN_SIGNATURE)
        .to_numpy()
    )
    fig, ax = plt.subplots(figsize=(8.4, 4.7))
    image = ax.imshow(matrix, cmap="RdBu_r", vmin=-2, vmax=2, aspect="auto")
    ax.set_xticks(range(len(HUMAN_SIGNATURE)), HUMAN_SIGNATURE)
    ax.set_yticks(range(len(labels)), labels)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            text = "NA" if np.isnan(matrix[i, j]) else f"{matrix[i, j]:.2f}"
            ax.text(j, i, text, ha="center", va="center", fontsize=8)
    ax.set_title("CLDN4 loss: prespecified IFN/MHC-I genes (log2 fold change)")
    fig.colorbar(image, ax=ax, label="KD/KO vs matched control")
    fig.tight_layout()
    fig.savefig(FIGURES / "claim_C4_heatmap.png", dpi=180)
    fig.savefig(FIGURES / "claim_C4_heatmap.svg")
    plt.close(fig)


def write_discovery() -> None:
    records = [
        ["GSE22493", "GEO", "included", "Direct CLDN4 lentiviral KD; paired two-colour expression arrays in SKOV-3"],
        ["GSE50927", "GEO", "included", "Direct germline Cldn4 KO; whole-lung RNA-seq at baseline and after VILI"],
        ["GSE207704", "GEO", "included", "Direct CRISPR CLDN4 KO; RNA-seq in MCF-7 and T47D"],
        ["LINCS CRISPR KO consensus", "LINCS/CFDE", "excluded", "No CLDN4 perturbation signature in the downloaded consensus GMT"],
        ["GSE99415/16/17", "GEO", "excluded", "CLDN4-related ncRNA study, but deposited assays are not direct CLDN4-loss transcriptomes"],
        ["GSE60885", "GEO", "excluded", "CLDN4 is a study result; deposited assay is DNA methylation, not CLDN4-loss transcriptomics"],
        ["GSE84742", "GEO", "excluded", "CLDN4 is mentioned in a broader differentiation study, not directly perturbed"],
        ["Kashiwagi et al. 2025, PMID 41016339", "publication/web search", "excluded", "H1688 CLDN4-KO RNA-seq reported, but no public expression accession or matrix was located"],
        ["Other CLDN4 KD/KO papers", "PubMed/web search", "excluded", "Targeted or phenotypic assays only; no public genome-wide expression matrix located"],
    ]
    pd.DataFrame(records, columns=["record", "source", "status", "reason"]).to_csv(
        PROCESSED / "discovery_audit.tsv", sep="\t", index=False
    )


def main() -> None:
    download()
    effects = pd.DataFrame(array_effects() + mouse_effects() + breast_effects())
    effects = adjust_array_pvalues(effects)
    effects.to_csv(PROCESSED / "gene_effects.tsv", sep="\t", index=False, na_rep="NA")
    summarize(effects).to_csv(
        PROCESSED / "contrast_summary.tsv", sep="\t", index=False, na_rep="NA"
    )
    check_lincs().to_csv(PROCESSED / "lincs_coverage.tsv", sep="\t", index=False)
    write_discovery()
    plot_heatmap(effects)


if __name__ == "__main__":
    main()
