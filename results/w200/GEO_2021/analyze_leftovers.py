#!/usr/bin/env python3
"""Leftover 2021 GEO lung ICI analysis: TACSTD2/CLDN4 plus deposited controls."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import urllib.request
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
from scipy.stats import mannwhitneyu, spearmanr

BASE190 = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn"
FILES = {
    "GSE190265_TPM_France3.csv.gz": f"{BASE190}/GSE190265/suppl/GSE190265_TPM_France3.csv.gz",
    "GSE190265_samples_info_France3.csv.gz": f"{BASE190}/GSE190265/suppl/GSE190265_samples_info_France3.csv.gz",
    "GSE190265_series_matrix.txt.gz": f"{BASE190}/GSE190265/matrix/GSE190265_series_matrix.txt.gz",
    "GSE190266_TPM_France4.csv.gz": f"{BASE190}/GSE190266/suppl/GSE190266_TPM_France4.csv.gz",
    "GSE190266_series_matrix.txt.gz": f"{BASE190}/GSE190266/matrix/GSE190266_series_matrix.txt.gz",
}

PRIMARY = ("TACSTD2", "CLDN4")
CONTROLS = ("CXCL10", "CD274", "OPTN", "TLR9")


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=Path("/tmp/GEO_2021"))
    parser.add_argument("--out-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--no-download", action="store_true")
    return parser.parse_args()


def obtain_files(data_dir: Path, allow_download: bool) -> None:
    data_dir.mkdir(parents=True, exist_ok=True)
    for filename, url in FILES.items():
        destination = data_dir / filename
        if destination.is_file() and destination.stat().st_size:
            continue
        if not allow_download:
            raise FileNotFoundError(destination)
        print(f"Downloading {url}")
        urllib.request.urlretrieve(url, destination)


def parse_series_matrix(path: Path) -> pd.DataFrame:
    titles: list[str] | None = None
    characteristics: list[list[str]] = []
    with gzip.open(path, "rt") as source:
        for line in source:
            fields = [value.strip('"') for value in line.rstrip("\n").split("\t")]
            if fields[0] == "!Sample_title":
                titles = fields[1:]
            elif fields[0] == "!Sample_characteristics_ch1":
                characteristics.append(fields[1:])
    if titles is None:
        raise ValueError(f"No sample titles in {path}")
    frame = pd.DataFrame(index=titles)
    for row in characteristics:
        key = row[0].split(": ", 1)[0]
        frame[key] = [value.split(": ", 1)[1] if ": " in value else value for value in row]
    return frame


def normalize_histology(value: object) -> str:
    text = str(value).strip().lower()
    if text in {"non-squamous", "nonsquamous"}:
        return "non-squamous"
    if text == "squamous":
        return "squamous"
    return "unknown"


def dcb_label(time: float, event: int) -> str:
    if time >= 6:
        return "DCB"
    if event == 1:
        return "NDB"
    return "excluded"


def load_france3(data_dir: Path) -> pd.DataFrame:
    expression = pd.read_csv(data_dir / "GSE190265_TPM_France3.csv.gz", sep=";", index_col=0)
    clinical = pd.read_csv(
        data_dir / "GSE190265_samples_info_France3.csv.gz", sep=";"
    ).set_index("sample")
    matrix = parse_series_matrix(data_dir / "GSE190265_series_matrix.txt.gz")
    frame = clinical.rename(columns={"time_PFS": "time", "evtPFS": "event"}).copy()
    frame["histology"] = matrix.reindex(frame.index)["disease state"].map(normalize_histology)
    frame["histology"] = frame["histology"].fillna("unknown")
    frame["pfs_capped"] = False
    frame["dataset"] = "GSE190265"
    genes = [gene for gene in PRIMARY + CONTROLS if gene in expression.columns]
    return frame.join(expression[genes], how="inner")


def load_france4(data_dir: Path) -> pd.DataFrame:
    expression = pd.read_csv(
        data_dir / "GSE190266_TPM_France4.csv.gz",
        sep=";",
        decimal=",",
        index_col=0,
    )
    matrix = parse_series_matrix(data_dir / "GSE190266_series_matrix.txt.gz")
    frame = pd.DataFrame(
        {
            "time": pd.to_numeric(matrix["pfs_time (6 months)"], errors="coerce"),
            "event": pd.to_numeric(matrix["pfs_evt (6 months)"], errors="coerce"),
            "histology": matrix["disease state"].map(normalize_histology),
        },
        index=matrix.index,
    )
    frame["pfs_capped"] = True
    frame["dataset"] = "GSE190266"
    genes = [gene for gene in PRIMARY + CONTROLS if gene in expression.columns]
    return frame.join(expression[genes], how="inner")


def gene_inventory(data_dir: Path) -> pd.DataFrame:
    france3 = pd.read_csv(data_dir / "GSE190265_TPM_France3.csv.gz", sep=";", index_col=0)
    france4 = pd.read_csv(
        data_dir / "GSE190266_TPM_France4.csv.gz", sep=";", decimal=",", index_col=0
    )
    rows = []
    for dataset, matrix in (("GSE190265", france3), ("GSE190266", france4)):
        rows.append(
            {
                "dataset": dataset,
                "n_samples": matrix.shape[0],
                "n_genes": matrix.shape[1],
                "first_gene": matrix.columns[0],
                "last_gene": matrix.columns[-1],
                "TACSTD2": "TACSTD2" in matrix.columns,
                "CLDN4": "CLDN4" in matrix.columns,
                "CXCL10": "CXCL10" in matrix.columns,
                "CD274": "CD274" in matrix.columns,
                "OPTN": "OPTN" in matrix.columns,
                "TLR9": "TLR9" in matrix.columns,
                "note": (
                    "Excel-style 16384-column cap; TACSTD2 sorts after T and is absent"
                    if dataset == "GSE190266"
                    else "full deposited transcriptome"
                ),
            }
        )
    return pd.DataFrame(rows)


def analyze_gene(frame: pd.DataFrame, gene: str) -> dict[str, object]:
    usable = frame[["time", "event", gene, "histology", "pfs_capped"]].dropna(
        subset=["time", "event", gene]
    )
    usable["log2"] = np.log2(usable[gene] + 1)
    usable["dcb"] = [dcb_label(time, int(event)) for time, event in zip(usable.time, usable.event)]
    cox = CoxPHFitter().fit(usable[["time", "event", "log2"]], duration_col="time", event_col="event")
    coefficient = cox.summary.loc["log2"]
    median = usable["log2"].median()
    high = usable["log2"] > median
    logrank = logrank_test(
        usable.loc[high, "time"],
        usable.loc[~high, "time"],
        usable.loc[high, "event"],
        usable.loc[~high, "event"],
    )
    rho, rho_p = spearmanr(usable["log2"], usable["time"])
    dcb = usable.loc[usable.dcb == "DCB", "log2"]
    ndb = usable.loc[usable.dcb == "NDB", "log2"]
    if len(dcb) and len(ndb):
        wilcoxon_p = mannwhitneyu(dcb, ndb, alternative="two-sided").pvalue
        auc = (np.array([(a > b) for a in dcb for b in ndb]).mean()
               + 0.5 * np.array([(a == b) for a in dcb for b in ndb]).mean())
    else:
        wilcoxon_p = np.nan
        auc = np.nan
    histology_p = np.nan
    labeled = usable[usable.histology.isin(["non-squamous", "squamous"])].copy()
    if labeled["histology"].nunique() == 2 and labeled["histology"].value_counts().min() >= 5:
        labeled["nonsquamous"] = (labeled.histology == "non-squamous").astype(int)
        adjusted = CoxPHFitter().fit(
            labeled[["time", "event", "log2", "nonsquamous"]],
            duration_col="time",
            event_col="event",
        )
        histology_p = float(adjusted.summary.loc["log2", "p"])
    return {
        "dataset": usable.attrs.get("dataset", frame["dataset"].iloc[0]),
        "gene": gene,
        "n": len(usable),
        "events": int(usable.event.sum()),
        "n_DCB": int((usable.dcb == "DCB").sum()),
        "n_NDB": int((usable.dcb == "NDB").sum()),
        "n_excluded_from_DCB": int((usable.dcb == "excluded").sum()),
        "median_TPM_DCB": float(usable.loc[usable.dcb == "DCB", gene].median()) if len(dcb) else np.nan,
        "median_TPM_NDB": float(usable.loc[usable.dcb == "NDB", gene].median()) if len(ndb) else np.nan,
        "auc_DCB_gt_NDB": auc,
        "wilcoxon_p_DCB": wilcoxon_p,
        "spearman_rho_vs_PFS": rho,
        "spearman_p": rho_p,
        "cox_hr_per_log2_tpm_plus_1": float(coefficient["exp(coef)"]),
        "cox_ci95_low": float(coefficient["exp(coef) lower 95%"]),
        "cox_ci95_high": float(coefficient["exp(coef) upper 95%"]),
        "cox_p": float(coefficient["p"]),
        "cox_p_histology_adjusted": histology_p,
        "median_split_logrank_p": float(logrank.p_value),
        "pfs_capped": bool(usable.pfs_capped.iloc[0]),
        "n_with_histology": int(usable.histology.isin(["non-squamous", "squamous"]).sum()),
    }


def add_bh(rows: list[dict[str, object]], p_field: str, q_field: str) -> None:
    valid = [
        (index, float(row[p_field]))
        for index, row in enumerate(rows)
        if p_field in row and pd.notna(row[p_field])
    ]
    if not valid:
        for row in rows:
            row[q_field] = np.nan
        return
    order = np.argsort([p_value for _, p_value in valid])
    p_values = np.array([valid[index][1] for index in order])
    ranked = p_values * len(p_values) / np.arange(1, len(p_values) + 1)
    adjusted = np.minimum.accumulate(ranked[::-1])[::-1].clip(max=1)
    mapping = {valid[order[i]][0]: adjusted[i] for i in range(len(order))}
    for index, row in enumerate(rows):
        row[q_field] = mapping.get(index, np.nan)


def write_per_sample(frame: pd.DataFrame, destination: Path) -> None:
    keep = ["dataset", "time", "event", "histology", "pfs_capped"]
    genes = [gene for gene in PRIMARY + CONTROLS if gene in frame.columns]
    out = frame[keep + genes].copy()
    out["dcb"] = [dcb_label(time, int(event)) for time, event in zip(out.time, out.event)]
    for gene in genes:
        out[f"log2_{gene}_plus_1"] = np.log2(out[gene] + 1)
    out.to_csv(destination)


def plot_km(frame: pd.DataFrame, gene: str, destination: Path) -> None:
    usable = frame[["time", "event", gene]].dropna()
    usable["log2"] = np.log2(usable[gene] + 1)
    high = usable["log2"] > usable["log2"].median()
    figure, axis = plt.subplots(figsize=(5.2, 3.8))
    fitter = KaplanMeierFitter()
    for mask, label in ((high, "above median"), (~high, "at/below median")):
        fitter.fit(usable.loc[mask, "time"], usable.loc[mask, "event"], label=f"{gene} {label}")
        fitter.plot(ax=axis, ci_show=False)
    axis.set_xlabel("PFS (months)")
    axis.set_ylabel("Progression-free probability")
    axis.set_title(f"{frame['dataset'].iloc[0]} {gene}")
    figure.tight_layout()
    figure.savefig(destination, dpi=140)
    plt.close(figure)


def plot_dcb_boxes(frames: list[pd.DataFrame], destination: Path) -> None:
    figure, axes = plt.subplots(1, 3, figsize=(9.6, 3.6))
    panels = [
        (frames[0], "TACSTD2"),
        (frames[0], "CLDN4"),
        (frames[1], "CLDN4"),
    ]
    for axis, (frame, gene) in zip(axes, panels):
        usable = frame[["time", "event", gene, "dataset"]].dropna()
        usable["dcb"] = [dcb_label(time, int(event)) for time, event in zip(usable.time, usable.event)]
        usable = usable[usable.dcb.isin(["DCB", "NDB"])]
        usable["log2"] = np.log2(usable[gene] + 1)
        data = [usable.loc[usable.dcb == label, "log2"] for label in ("NDB", "DCB")]
        axis.boxplot(data, tick_labels=["NDB", "DCB"], widths=0.55)
        axis.set_title(f"{usable.dataset.iloc[0]} {gene}")
        axis.set_ylabel("log2(TPM+1)")
    figure.tight_layout()
    figure.savefig(destination, dpi=140)
    plt.close(figure)


def write_manifest(data_dir: Path, destination: Path) -> None:
    with destination.open("w", newline="") as output:
        writer = csv.writer(output)
        writer.writerow(["file", "bytes", "sha256", "url"])
        for filename, url in FILES.items():
            path = data_dir / filename
            writer.writerow(
                [
                    filename,
                    path.stat().st_size,
                    hashlib.sha256(path.read_bytes()).hexdigest(),
                    url,
                ]
            )


def main() -> None:
    args = arguments()
    obtain_files(args.data_dir, allow_download=not args.no_download)
    out = args.out_dir
    figures = out / "figures"
    figures.mkdir(parents=True, exist_ok=True)

    france3 = load_france3(args.data_dir)
    france4 = load_france4(args.data_dir)
    inventory = gene_inventory(args.data_dir)
    inventory.to_csv(out / "gene_inventory.csv", index=False)

    write_per_sample(france3, out / "per_sample_GSE190265.csv")
    write_per_sample(france4, out / "per_sample_GSE190266.csv")

    primary_rows: list[dict[str, object]] = []
    for gene in PRIMARY:
        if gene in france3.columns:
            primary_rows.append(analyze_gene(france3, gene))
        else:
            primary_rows.append(
                {"dataset": "GSE190265", "gene": gene, "note": "absent from deposited matrix"}
            )
        if gene in france4.columns:
            primary_rows.append(analyze_gene(france4, gene))
        else:
            primary_rows.append(
                {
                    "dataset": "GSE190266",
                    "gene": gene,
                    "n": 70,
                    "note": "absent from deposited TPM matrix (last gene MTMR14; 16383 columns)",
                }
            )
    add_bh(primary_rows, "cox_p", "cox_bh_q_primary")
    add_bh(primary_rows, "wilcoxon_p_DCB", "wilcoxon_bh_q_primary")
    pd.DataFrame(primary_rows).to_csv(out / "survival_results.csv", index=False)

    control_rows: list[dict[str, object]] = []
    for frame in (france3, france4):
        for gene in CONTROLS:
            if gene in frame.columns:
                control_rows.append(analyze_gene(frame, gene))
    add_bh(control_rows, "wilcoxon_p_DCB", "wilcoxon_bh_q_controls")
    pd.DataFrame(control_rows).to_csv(out / "control_results.csv", index=False)

    plot_km(france3, "TACSTD2", figures / "GSE190265_TACSTD2_KM.png")
    plot_km(france3, "CLDN4", figures / "GSE190265_CLDN4_KM.png")
    plot_km(france4, "CLDN4", figures / "GSE190266_CLDN4_KM.png")
    plot_dcb_boxes([france3, france4], figures / "DCB_boxplots.png")
    write_manifest(args.data_dir, out / "input_manifest.csv")
    print(f"Wrote leftover analyses to {out}")


if __name__ == "__main__":
    main()
