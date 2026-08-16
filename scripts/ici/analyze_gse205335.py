#!/usr/bin/env python3
"""Patient-level malignant TACSTD2/CLDN4 tests in GSE205335 core ICI samples."""

from __future__ import annotations

from pathlib import Path
import subprocess
import tempfile
import urllib.request

from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "ici"
IDENTITY_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/"
    "GSE205335_Lung_IO_CellIdentity.txt.gz"
)
MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/"
    "GSE205335_Lung_IO_UMI_matrix.rds.gz"
)
CLINICAL_URL = (
    "https://cdn.elifesciences.org/articles/98366/elife-98366-supp1-v1.xlsx"
)
MIN_MALIGNANT = 10
IFN_MHC = ("HLA-A", "HLA-B", "HLA-C", "B2M", "STAT1", "IRF1", "IFITM1", "CXCL10")


def download(url: str, path: Path) -> Path:
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_suffix(path.suffix + ".part")
        urllib.request.urlretrieve(url, temporary)
        temporary.replace(path)
    return path


def extract_genes(matrix: Path, destination: Path) -> Path:
    if destination.exists() and destination.stat().st_size > 0:
        return destination
    script = ROOT / "scripts" / "ici" / "extract_gse205335_genes.R"
    subprocess.run(
        ["Rscript", str(script), str(matrix), str(destination)],
        check=True,
    )
    return destination


def sample_key(orig_ident: str) -> str:
    return orig_ident.rsplit("-", 1)[0].replace("-", "_")


def bh_adjust(values: pd.Series) -> pd.Series:
    p = values.astype(float).to_numpy()
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    adjusted = np.empty(len(p))
    adjusted[order] = np.minimum.accumulate(ranked[::-1])[::-1]
    return pd.Series(np.minimum(adjusted, 1), index=values.index)


def load_clinical(path: Path) -> pd.DataFrame:
    raw = pd.read_excel(path, sheet_name=0, engine="calamine", header=None)
    table = raw.iloc[3:].copy()
    table.columns = raw.iloc[2].tolist()
    table = table.dropna(subset=["Sample"])
    table = table[table["Sample"].astype(str).str.match(r"^[A-Z]")]
    table["sample"] = table["Sample"].astype(str)
    return table


def patient_scores(identities: pd.DataFrame, genes: pd.DataFrame) -> pd.DataFrame:
    identities = identities.copy()
    identities["sample"] = identities["orig.ident"].map(sample_key)
    merged = identities.merge(genes, on="barcode", how="inner")
    if len(merged) != len(identities):
        raise RuntimeError(
            f"Barcode mismatch: identities={len(identities)} merged={len(merged)}"
        )
    malignant = merged["lineage.sub"].eq("Malignant cells")
    tnk = merged["lineage.total"].eq("T/NK cells")
    rows = []
    for sample, cells in merged.groupby("sample", sort=True):
        malignant_cells = cells[cells["lineage.sub"].eq("Malignant cells")]
        non_malignant = cells[cells["lineage.sub"].ne("Malignant cells")]
        record = {
            "sample": sample,
            "core_add": cells["core.patient"].iloc[0],
            "n_cells": len(cells),
            "n_malignant": int(len(malignant_cells)),
            "tnk_fraction_all_cells": float(cells["lineage.total"].eq("T/NK cells").mean()),
            "tnk_fraction_non_malignant": (
                float(non_malignant["lineage.total"].eq("T/NK cells").mean())
                if len(non_malignant)
                else np.nan
            ),
        }
        for gene in ("TACSTD2", "CLDN4"):
            values = malignant_cells[gene].to_numpy(float) if len(malignant_cells) else np.array([])
            record[f"{gene}_malignant_mean_log1p"] = (
                float(np.log1p(values).mean()) if len(values) else np.nan
            )
            record[f"{gene}_malignant_detection_fraction"] = (
                float((values > 0).mean()) if len(values) else np.nan
            )
        if len(malignant_cells):
            record["malignant_ifn_mhci_mean_log1p"] = float(
                np.mean([np.log1p(malignant_cells[gene]).mean() for gene in IFN_MHC])
            )
        else:
            record["malignant_ifn_mhci_mean_log1p"] = np.nan
        rows.append(record)
    _ = malignant, tnk
    return pd.DataFrame(rows)


def collapse_patients(scores: pd.DataFrame, clinical: pd.DataFrame) -> pd.DataFrame:
    merged = scores.merge(clinical, on="sample", how="left")
    core = merged[merged["Class"].eq("Core")].copy()
    if core["Patient"].nunique() != 11:
        raise RuntimeError(f"Expected 11 core patients, found {core['Patient'].nunique()}")
    rows = []
    for patient, samples in core.groupby("Patient", sort=True):
        record = {
            "patient": patient,
            "n_samples": len(samples),
            "samples": ",".join(samples["sample"]),
            "recist": samples["RECIST"].iloc[0],
            "response": samples["Response"].iloc[0],
            "histology": samples["Cancer subtype"].iloc[0],
            "immunotherapy": samples["Immunotherapy"].iloc[0],
            "tissue_origins": ",".join(samples["Tissue origin"].astype(str)),
            "pfs_months": float(samples["PFS"].iloc[0]),
            "pfs_event": int(samples["PFS_event"].iloc[0]),
            "os_months": float(samples["OS"].iloc[0]),
            "os_event": int(samples["OS_event"].iloc[0]),
            "n_cells": int(samples["n_cells"].sum()),
            "n_malignant": int(samples["n_malignant"].sum()),
        }
        malignant_weights = samples["n_malignant"].clip(lower=1)
        cell_weights = samples["n_cells"].clip(lower=1)
        non_malignant_weights = (samples["n_cells"] - samples["n_malignant"]).clip(lower=1)
        for column in (
            "TACSTD2_malignant_mean_log1p",
            "CLDN4_malignant_mean_log1p",
            "TACSTD2_malignant_detection_fraction",
            "CLDN4_malignant_detection_fraction",
            "malignant_ifn_mhci_mean_log1p",
        ):
            record[column] = float(np.average(samples[column], weights=malignant_weights))
        record["tnk_fraction_all_cells"] = float(
            np.average(samples["tnk_fraction_all_cells"], weights=cell_weights)
        )
        record["tnk_fraction_non_malignant"] = float(
            np.average(
                samples["tnk_fraction_non_malignant"], weights=non_malignant_weights
            )
        )
        rows.append(record)
    patients = pd.DataFrame(rows)
    patients["responder"] = patients["response"].eq("Responder")
    patients["included_min10_malignant"] = patients["n_malignant"] >= MIN_MALIGNANT
    return core, patients


def statistics_table(patients: pd.DataFrame) -> pd.DataFrame:
    analysis = patients[patients["included_min10_malignant"]].copy()
    rows = []
    for gene in ("TACSTD2", "CLDN4"):
        column = f"{gene}_malignant_mean_log1p"
        non = analysis.loc[~analysis["responder"], column]
        yes = analysis.loc[analysis["responder"], column]
        _, p_value = mannwhitneyu(non, yes, alternative="two-sided")
        rows.append(
            {
                "analysis": "core malignant score; non-responder vs PR",
                "feature": gene,
                "n": len(analysis),
                "n_group1": len(non),
                "n_group2": len(yes),
                "effect": float(non.median() - yes.median()),
                "effect_definition": "median(SD+PD)-median(PR), mean log1p UMI/cell",
                "rho": np.nan,
                "hr": np.nan,
                "p": p_value,
            }
        )
        pd_only = analysis[analysis["recist"].isin(["PR", "PD"])]
        high = pd_only.loc[pd_only["recist"].eq("PD"), column]
        low = pd_only.loc[pd_only["recist"].eq("PR"), column]
        _, p_pd = mannwhitneyu(high, low, alternative="two-sided")
        rows.append(
            {
                "analysis": "core malignant score; PD vs PR",
                "feature": gene,
                "n": len(pd_only),
                "n_group1": len(high),
                "n_group2": len(low),
                "effect": float(high.median() - low.median()),
                "effect_definition": "median(PD)-median(PR), mean log1p UMI/cell",
                "rho": np.nan,
                "hr": np.nan,
                "p": p_pd,
            }
        )
        for outcome in (
            "tnk_fraction_all_cells",
            "tnk_fraction_non_malignant",
            "malignant_ifn_mhci_mean_log1p",
        ):
            rho, p_rho = spearmanr(analysis[column], analysis[outcome])
            rows.append(
                {
                    "analysis": f"core malignant {gene} vs {outcome}",
                    "feature": gene,
                    "n": len(analysis),
                    "n_group1": np.nan,
                    "n_group2": np.nan,
                    "effect": np.nan,
                    "effect_definition": "Spearman rho",
                    "rho": rho,
                    "hr": np.nan,
                    "p": p_rho,
                }
            )
        cox = analysis[[column, "pfs_months", "pfs_event"]].dropna().copy()
        cox["z"] = (cox[column] - cox[column].mean()) / cox[column].std(ddof=1)
        fitter = CoxPHFitter()
        fitter.fit(cox[["z", "pfs_months", "pfs_event"]], "pfs_months", "pfs_event")
        rows.append(
            {
                "analysis": "core malignant score; PFS Cox per SD",
                "feature": gene,
                "n": len(cox),
                "n_group1": int(cox["pfs_event"].sum()),
                "n_group2": np.nan,
                "effect": np.nan,
                "effect_definition": "Cox HR per 1 SD",
                "rho": np.nan,
                "hr": float(np.exp(fitter.params_["z"])),
                "p": float(fitter.summary.loc["z", "p"]),
            }
        )
        high = analysis[column] >= analysis[column].median()
        logrank = logrank_test(
            analysis.loc[high, "pfs_months"],
            analysis.loc[~high, "pfs_months"],
            analysis.loc[high, "pfs_event"],
            analysis.loc[~high, "pfs_event"],
        )
        rows.append(
            {
                "analysis": "core malignant score; PFS median-split log-rank",
                "feature": gene,
                "n": len(analysis),
                "n_group1": int(high.sum()),
                "n_group2": int((~high).sum()),
                "effect": np.nan,
                "effect_definition": "log-rank high vs low",
                "rho": np.nan,
                "hr": np.nan,
                "p": float(logrank.p_value),
            }
        )
    rho, p_value = spearmanr(
        analysis["TACSTD2_malignant_mean_log1p"],
        analysis["CLDN4_malignant_mean_log1p"],
    )
    rows.append(
        {
            "analysis": "core malignant TACSTD2 vs CLDN4",
            "feature": "TACSTD2~CLDN4",
            "n": len(analysis),
            "n_group1": np.nan,
            "n_group2": np.nan,
            "effect": np.nan,
            "effect_definition": "Spearman rho",
            "rho": rho,
            "hr": np.nan,
            "p": p_value,
        }
    )
    table = pd.DataFrame(rows)
    table["q_bh"] = bh_adjust(table["p"])
    return table


def plot_response(patients: pd.DataFrame, statistics: pd.DataFrame) -> None:
    analysis = patients[patients["included_min10_malignant"]]
    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    rng = np.random.default_rng(20260816)
    for ax, gene in zip(axes, ("TACSTD2", "CLDN4")):
        column = f"{gene}_malignant_mean_log1p"
        groups = [
            analysis.loc[~analysis["responder"], column],
            analysis.loc[analysis["responder"], column],
        ]
        ax.boxplot(
            groups,
            tick_labels=[
                f"SD+PD (n={len(groups[0])})",
                f"PR (n={len(groups[1])})",
            ],
        )
        for position, group in enumerate(groups, 1):
            ax.scatter(
                position + rng.uniform(-0.06, 0.06, len(group)),
                group,
                color="#276FBF",
            )
        row = statistics[
            statistics["analysis"].eq("core malignant score; non-responder vs PR")
            & statistics["feature"].eq(gene)
        ].iloc[0]
        ax.set(
            title=f"GSE205335 core: {gene}\n(SD+PD)−PR={row.effect:.3f}, p={row.p:.3g}",
            ylabel="Malignant-cell mean log1p UMI",
        )
    fig.tight_layout()
    fig.savefig(OUT / "gse205335_response.png", dpi=180)
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    cache = Path(tempfile.gettempdir()) / "ici_next"
    cache.mkdir(parents=True, exist_ok=True)
    identity_path = download(IDENTITY_URL, cache / IDENTITY_URL.rsplit("/", 1)[-1])
    matrix_path = download(MATRIX_URL, cache / MATRIX_URL.rsplit("/", 1)[-1])
    clinical_path = download(CLINICAL_URL, cache / "elife-98366-supp1-v1.xlsx")
    gene_path = extract_genes(matrix_path, cache / "GSE205335_selected_genes.tsv")

    identities = pd.read_csv(identity_path, sep="\t")
    genes = pd.read_csv(gene_path, sep="\t")
    clinical = load_clinical(clinical_path)
    scores = patient_scores(identities, genes)
    sample_table, patients = collapse_patients(scores, clinical)
    sample_table.to_csv(OUT / "gse205335_sample_scores.tsv", sep="\t", index=False)
    patients.to_csv(OUT / "gse205335_patient_scores.tsv", sep="\t", index=False)
    statistics = statistics_table(patients)
    statistics.to_csv(OUT / "gse205335_statistics.tsv", sep="\t", index=False)
    plot_response(patients, statistics)


if __name__ == "__main__":
    main()
