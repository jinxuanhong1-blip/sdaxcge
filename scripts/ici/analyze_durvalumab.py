#!/usr/bin/env python3
"""Pretreatment TACSTD2/CLDN4 tests in NCT02904954 durvalumab ± SBRT RNA-seq."""

from __future__ import annotations

import gzip
from pathlib import Path
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
PRE_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE253nnn/GSE253564/suppl/"
    "GSE253564_Pre-treatment_Samples_Pubs_FPKMs.txt.gz"
)
SOFT_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE253nnn/GSE253564/soft/"
    "GSE253564_family.soft.gz"
)
SOURCE_URL = (
    "https://static-content.springer.com/esm/"
    "art%3A10.1038%2Fs41467-023-44195-x/MediaObjects/"
    "41467_2023_44195_MOESM6_ESM.xlsx"
)
GENES = ("TACSTD2", "CLDN4")
IMMUNE = ("CD3D", "CD3E", "CD8A")


def download(url: str, path: Path) -> Path:
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_suffix(path.suffix + ".part")
        urllib.request.urlretrieve(url, temporary)
        temporary.replace(path)
    return path


def parse_soft(path: Path) -> pd.DataFrame:
    records: list[dict[str, str]] = []
    current: dict[str, str] | None = None
    with gzip.open(path, "rt", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current:
                    records.append(current)
                current = {"gsm": line.split(" = ", 1)[1]}
            elif current is not None and line.startswith("!Sample_title = "):
                current["title"] = line.split(" = ", 1)[1]
            elif current is not None and line.startswith(
                "!Sample_characteristics_ch1 = "
            ):
                value = line.split(" = ", 1)[1]
                if ": " in value:
                    key, item = value.split(": ", 1)
                    current[key.strip().lower()] = item.strip()
    if current:
        records.append(current)
    return pd.DataFrame(records)


def bh_adjust(values: pd.Series) -> pd.Series:
    p = values.astype(float).to_numpy()
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    adjusted = np.empty(len(p))
    adjusted[order] = np.minimum.accumulate(ranked[::-1])[::-1]
    return pd.Series(np.minimum(adjusted, 1), index=values.index)


def load_source_outcomes(path: Path) -> pd.DataFrame:
    mpr = pd.read_excel(path, sheet_name="Figure 2e")
    mpr = mpr.dropna(subset=["Reference number"]).copy()
    mpr["patient"] = mpr["Reference number"].astype(str).str.replace(
        "Durva", "durva", regex=False
    )
    mpr["mpr"] = mpr["Path Response_2Group"].eq("Major")
    pfs = pd.read_excel(path, sheet_name="Figure 2b")
    pfs = pfs.dropna(subset=["Reference number"]).copy()
    pfs["patient"] = pfs["Reference number"].astype(str).str.replace(
        "Durva", "durva", regex=False
    )
    status = pfs["Progression Status"].astype(str)
    pfs["pfs_months"] = pd.to_numeric(
        pfs["Progression-free survival in Months"], errors="coerce"
    )
    pfs["pfs_event"] = status.str.contains(
        r"recurrence|Dead|Death", case=False, regex=True
    ).astype(int)
    pfs["source_arm"] = pfs["ARM"].astype(str)
    return mpr[["patient", "mpr"]].merge(
        pfs[["patient", "pfs_months", "pfs_event", "source_arm", "Progression Status"]],
        on="patient",
        how="outer",
    )


def infer_arm_label(geo_arm: str) -> str:
    # GEO only stores Arm1/Arm2. Published MPR imbalance (2 vs 16) maps Arm2 to dual therapy.
    return {
        "Arm1": "durvalumab_monotherapy",
        "Arm2": "durvalumab_SBRT",
    }.get(geo_arm, geo_arm)


def statistics_table(data: pd.DataFrame) -> pd.DataFrame:
    rows = []
    subsets = {
        "all pretreatment": data,
        "durvalumab monotherapy": data[data["arm_label"].eq("durvalumab_monotherapy")],
        "durvalumab + SBRT": data[data["arm_label"].eq("durvalumab_SBRT")],
    }
    for subset_name, subset in subsets.items():
        labeled = subset.dropna(subset=["mpr"])
        if labeled["mpr"].nunique() < 2:
            continue
        for gene in GENES:
            nmpr = labeled.loc[~labeled["mpr"], gene]
            mpr = labeled.loc[labeled["mpr"], gene]
            _, p_value = mannwhitneyu(nmpr, mpr, alternative="two-sided")
            rows.append(
                {
                    "analysis": f"{subset_name}; NMPR vs MPR",
                    "feature": gene,
                    "n": len(labeled),
                    "n_group1": len(nmpr),
                    "n_group2": len(mpr),
                    "effect": float(nmpr.median() - mpr.median()),
                    "effect_definition": "median(NMPR)-median(MPR), log2(FPKM+1)",
                    "rho": np.nan,
                    "hr": np.nan,
                    "p": p_value,
                }
            )
            rho, p_rho = spearmanr(labeled[gene], labeled["tnk_score"])
            rows.append(
                {
                    "analysis": f"{subset_name}; {gene} vs T-cell proxy",
                    "feature": gene,
                    "n": len(labeled),
                    "n_group1": np.nan,
                    "n_group2": np.nan,
                    "effect": np.nan,
                    "effect_definition": "Spearman rho",
                    "rho": rho,
                    "hr": np.nan,
                    "p": p_rho,
                }
            )
        rho, p_value = spearmanr(labeled["TACSTD2"], labeled["CLDN4"])
        rows.append(
            {
                "analysis": f"{subset_name}; TACSTD2 vs CLDN4",
                "feature": "TACSTD2~CLDN4",
                "n": len(labeled),
                "n_group1": np.nan,
                "n_group2": np.nan,
                "effect": np.nan,
                "effect_definition": "Spearman rho",
                "rho": rho,
                "hr": np.nan,
                "p": p_value,
            }
        )
        survival = labeled.dropna(subset=["pfs_months", "pfs_event"])
        if survival["pfs_event"].sum() >= 3:
            for gene in GENES:
                cox = survival[[gene, "pfs_months", "pfs_event"]].copy()
                cox["z"] = (cox[gene] - cox[gene].mean()) / cox[gene].std(ddof=1)
                fitter = CoxPHFitter()
                fitter.fit(cox[["z", "pfs_months", "pfs_event"]], "pfs_months", "pfs_event")
                rows.append(
                    {
                        "analysis": f"{subset_name}; PFS Cox per SD",
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
                high = survival[gene] >= survival[gene].median()
                logrank = logrank_test(
                    survival.loc[high, "pfs_months"],
                    survival.loc[~high, "pfs_months"],
                    survival.loc[high, "pfs_event"],
                    survival.loc[~high, "pfs_event"],
                )
                rows.append(
                    {
                        "analysis": f"{subset_name}; PFS median-split log-rank",
                        "feature": gene,
                        "n": len(survival),
                        "n_group1": int(high.sum()),
                        "n_group2": int((~high).sum()),
                        "effect": np.nan,
                        "effect_definition": "log-rank high vs low",
                        "rho": np.nan,
                        "hr": np.nan,
                        "p": float(logrank.p_value),
                    }
                )
    table = pd.DataFrame(rows)
    table["q_bh"] = bh_adjust(table["p"])
    return table


def plot_mpr(data: pd.DataFrame, statistics: pd.DataFrame) -> None:
    labeled = data.dropna(subset=["mpr"])
    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    rng = np.random.default_rng(20260816)
    for ax, gene in zip(axes, GENES):
        groups = [
            labeled.loc[~labeled["mpr"], gene],
            labeled.loc[labeled["mpr"], gene],
        ]
        ax.boxplot(
            groups,
            tick_labels=[f"NMPR (n={len(groups[0])})", f"MPR (n={len(groups[1])})"],
        )
        for position, group in enumerate(groups, 1):
            ax.scatter(
                position + rng.uniform(-0.06, 0.06, len(group)),
                group,
                color="#B23A48",
            )
        row = statistics[
            statistics["analysis"].eq("all pretreatment; NMPR vs MPR")
            & statistics["feature"].eq(gene)
        ].iloc[0]
        ax.set(
            title=f"GSE253564 pretreatment: {gene}\nNMPR−MPR={row.effect:.3f}, p={row.p:.3g}",
            ylabel="log2(FPKM+1)",
        )
    fig.tight_layout()
    fig.savefig(OUT / "durvalumab_mpr.png", dpi=180)
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    cache = Path(tempfile.gettempdir()) / "ici_next"
    cache.mkdir(parents=True, exist_ok=True)
    pre_path = download(PRE_URL, cache / PRE_URL.rsplit("/", 1)[-1])
    soft_path = download(SOFT_URL, cache / SOFT_URL.rsplit("/", 1)[-1])
    source_path = download(SOURCE_URL, cache / SOURCE_URL.rsplit("/", 1)[-1])

    expression = pd.read_csv(pre_path, sep="\t")
    wanted = list(GENES + IMMUNE)
    missing = set(wanted) - set(expression["gene"])
    if missing:
        raise RuntimeError(f"Missing genes in GSE253564: {sorted(missing)}")
    matrix = (
        expression.set_index("gene")
        .loc[wanted]
        .drop(columns=["Entrez.ID"])
        .apply(pd.to_numeric)
    )
    log_expression = np.log2(matrix.astype(float) + 1).T
    log_expression.index.name = "patient"
    log_expression["tnk_score"] = log_expression[list(IMMUNE)].apply(
        lambda column: (column - column.mean()) / column.std(ddof=1)
    ).mean(axis=1)

    soft = parse_soft(soft_path)
    soft["patient"] = soft["title"].astype(str)
    soft["arm_label"] = soft["treatment"].map(infer_arm_label)
    outcomes = load_source_outcomes(source_path)
    data = (
        log_expression.reset_index()
        .merge(soft[["patient", "gsm", "treatment", "cell type"]], on="patient")
        .merge(outcomes, on="patient", how="left")
    )
    arm_from_source = data["source_arm"].map(
        {
            "Durvalumab": "durvalumab_monotherapy",
            "SBRT + durva": "durvalumab_SBRT",
        }
    )
    data["arm_label"] = arm_from_source.fillna(data["treatment"].map(infer_arm_label))
    data.to_csv(OUT / "gse253564_patient_scores.tsv", sep="\t", index=False)
    statistics = statistics_table(data)
    statistics.to_csv(OUT / "gse253564_statistics.tsv", sep="\t", index=False)
    plot_mpr(data, statistics)


if __name__ == "__main__":
    main()
