#!/usr/bin/env python3
"""Score TACSTD2 against published Leader/Salcher signatures on public processed data."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


TOIL_DEFAULT = Path("/tmp/xena/tcga_RSEM_gene_tpm.gz")
ENSEMBL_DEFAULT = Path("/tmp/leader-upstream/input_tables/ensemble_ids.tsv")
THORSSON_DEFAULT = Path(
    "/tmp/leader-upstream/input_tables/Thorsson_et_al_TCGA_tumor_immunology_metrics.csv"
)
ESTIMATE_DEFAULT = Path("/tmp/leader-upstream/input_tables/TCGA_ESTIMATE.txt")
IMMUNE_EP_DEFAULT = Path("/tmp/leader-upstream/input_tables/immune_vs_ep_de.csv")
LCAM_DE_DEFAULT = Path(
    "/tmp/leader-upstream/input_tables/DE_LCAMhi_vs_LCAMlo_pseudobulk.rd"
)
SIGNATURES_DEFAULT = Path(__file__).resolve().with_name("signatures.tsv")

LCAM_HI = ["IgG_plasma", "SPP1_mac", "T_activated"]
LCAM_LO = ["B", "AM", "cDC2", "AZU1_mac", "cDC1"]
LCAM_SUBTYPE_ROWS = {
    "Leader_LCAM_hi_IgG_plasma": "IgG_plasma",
    "Leader_LCAM_hi_SPP1_mac": "SPP1_mac",
    "Leader_LCAM_hi_T_activated": "T_activated",
    "Leader_LCAM_lo_B": "B",
    "Leader_LCAM_lo_AM": "AM",
    "Leader_LCAM_lo_cDC2": "cDC2",
    "Leader_LCAM_lo_AZU1_mac": "AZU1_mac",
    "Leader_LCAM_lo_cDC1": "cDC1",
}
COMPARATOR_GENES = [
    "TACSTD2",
    "EPCAM",
    "KRT8",
    "CD3D",
    "CXCL13",
    "SPP1",
    "CXCR2",
    "FCGR3B",
]


def load_signatures(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def load_ensembl_map(path: Path) -> dict[str, str]:
    table = pd.read_csv(path, sep="\t", header=None, names=["symbol", "ensembl"])
    return dict(zip(table["ensembl"], table["symbol"]))


def needed_symbols(signatures: pd.DataFrame) -> set[str]:
    genes = {"TACSTD2"}
    for row in signatures.itertuples(index=False):
        genes.update(row.genes.split(","))
    return genes


def extract_toil(
    toil_path: Path, ensembl_to_symbol: dict[str, str], symbols: set[str]
) -> pd.DataFrame:
    wanted_ensembl = {
        ensembl: symbol
        for ensembl, symbol in ensembl_to_symbol.items()
        if symbol in symbols
    }
    rows: dict[str, np.ndarray] = {}
    with gzip.open(toil_path, "rt") as handle:
        samples = handle.readline().rstrip("\n").split("\t")[1:]
        for line in handle:
            ensembl_version, *values = line.rstrip("\n").split("\t")
            ensembl = ensembl_version.split(".")[0]
            symbol = wanted_ensembl.get(ensembl)
            if symbol is None or symbol in rows:
                continue
            rows[symbol] = np.asarray(values, dtype=np.float64)
    if "TACSTD2" not in rows:
        raise RuntimeError("TACSTD2 was not found in the TOIL matrix")
    return pd.DataFrame(rows, index=samples).T


def unlog_toil(log2_tpm: pd.DataFrame) -> pd.DataFrame:
    linear = np.power(2.0, log2_tpm) - 0.001
    return pd.DataFrame(np.clip(linear, 0.0, None), index=log2_tpm.index, columns=log2_tpm.columns)


def zscore_rows(matrix: pd.DataFrame) -> pd.DataFrame:
    values = matrix.to_numpy(dtype=np.float64)
    means = np.nanmean(values, axis=1, keepdims=True)
    stds = np.nanstd(values, axis=1, ddof=0, keepdims=True)
    stds[stds == 0] = np.nan
    return pd.DataFrame((values - means) / stds, index=matrix.index, columns=matrix.columns)


def mean_z_score(linear: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, list[str]]:
    present = [gene for gene in genes if gene in linear.index]
    if not present:
        return pd.Series(np.nan, index=linear.columns), present
    logged = np.log1p(linear.loc[present])
    return zscore_rows(logged).mean(axis=0), present


def lcam_scores(linear: pd.DataFrame, signatures: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    subtype_genes: dict[str, list[str]] = {}
    present: dict[str, list[str]] = {}
    for row_name, subtype in LCAM_SUBTYPE_ROWS.items():
        genes = signatures.loc[signatures["signature"] == row_name, "genes"].iloc[0].split(",")
        subtype_genes[subtype] = genes
        present[subtype] = [gene for gene in genes if gene in linear.index]

    fractions = linear / linear.sum(axis=0).replace(0, np.nan)
    logged = np.log10(1e-6 + fractions)
    gene_z = zscore_rows(logged)

    subtype_mat = {}
    for subtype, genes in present.items():
        if not genes:
            subtype_mat[subtype] = pd.Series(np.nan, index=linear.columns)
            continue
        subtype_mat[subtype] = gene_z.loc[genes].mean(axis=0)
    subtype_df = pd.DataFrame(subtype_mat).T
    subtype_z = zscore_rows(subtype_df)

    hi = [name for name in LCAM_HI if present[name]]
    lo = [name for name in LCAM_LO if present[name]]
    scores = pd.DataFrame(
        {
            "LCAMhi": subtype_z.loc[hi].mean(axis=0) if hi else np.nan,
            "LCAMlo": subtype_z.loc[lo].mean(axis=0) if lo else np.nan,
        }
    )
    scores["LCAM"] = scores["LCAMhi"] - scores["LCAMlo"]
    return scores, present


def primary_nsclc(samples: list[str], thorsson: pd.DataFrame) -> pd.DataFrame:
    study = thorsson.set_index("TCGA Participant Barcode")["TCGA Study"]
    rows = []
    for sample in samples:
        if not sample.startswith("TCGA-") or len(sample) < 15:
            continue
        if sample[13:15] != "01":
            continue
        participant = sample[:12]
        disease = study.get(participant)
        if disease in {"LUAD", "LUSC"}:
            rows.append(
                {
                    "sample": sample,
                    "participant": participant,
                    "disease": disease,
                    "leukocyte_fraction": thorsson.set_index("TCGA Participant Barcode")
                    .loc[participant, "Leukocyte Fraction"]
                    if participant in thorsson["TCGA Participant Barcode"].values
                    else np.nan,
                }
            )
    return pd.DataFrame(rows)


def spearman(x: pd.Series, y: pd.Series) -> dict[str, float | int]:
    aligned = pd.concat([x, y], axis=1).dropna()
    if len(aligned) < 3:
        return {"n": int(len(aligned)), "rho": np.nan, "p": np.nan}
    rho, p_value = stats.spearmanr(aligned.iloc[:, 0], aligned.iloc[:, 1])
    return {"n": int(len(aligned)), "rho": float(rho), "p": float(p_value)}


def extract_leader_tables(output_dir: Path) -> None:
    immune_ep = pd.read_csv(IMMUNE_EP_DEFAULT, index_col=0)
    immune_ep.index = immune_ep.index.astype(str)
    rows = []
    for gene in COMPARATOR_GENES:
        if gene not in immune_ep.index:
            continue
        row = immune_ep.loc[gene]
        rows.append(
            {
                "source": "Leader_immune_vs_epithelial",
                "gene": gene,
                "l2fc_immune_over_epithelial": float(row["l2fc"]),
                "immune_expression": float(row["fg_exprs"]),
                "epithelial_expression": float(row["bg_exprs"]),
                "note": "columns are expression, not p-values",
            }
        )
    immune_out = pd.DataFrame(rows)
    immune_out.to_csv(
        output_dir / "leader_immune_vs_epithelial.tsv",
        sep="\t",
        index=False,
        float_format="%.6g",
    )

    import rdata

    de = rdata.read_rda(LCAM_DE_DEFAULT)["DE_total"]
    de.index = de.index.astype(str)
    if "TACSTD2" not in de.index:
        raise RuntimeError("TACSTD2 missing from Leader LCAM-hi vs LCAM-lo DE")
    row = de.loc["TACSTD2"]
    pd.DataFrame(
        [
            {
                "source": "Leader_LCAMhi_vs_LCAMlo_pseudobulk",
                "gene": "TACSTD2",
                "log2_FC_hi_over_lo": float(row["log2_FC"]),
                "p_value": float(row["p.value"]),
                "adj_p_value": float(row["adj.p.value"]),
                "freq_hi": float(row["freq_fg"]),
                "freq_lo": float(row["freq_bg"]),
            }
        ]
    ).to_csv(
        output_dir / "leader_lcam_de_tacstd2.tsv",
        sep="\t",
        index=False,
        float_format="%.6g",
    )


def score_tcga(output_dir: Path, toil_path: Path) -> None:
    signatures = load_signatures(SIGNATURES_DEFAULT)
    symbols = needed_symbols(signatures)
    ensembl_to_symbol = load_ensembl_map(ENSEMBL_DEFAULT)
    log2_tpm = extract_toil(toil_path, ensembl_to_symbol, symbols)
    linear = unlog_toil(log2_tpm)

    thorsson = pd.read_csv(THORSSON_DEFAULT)
    meta = primary_nsclc(list(linear.columns), thorsson)
    linear = linear.loc[:, meta["sample"]]
    log2_tpm = log2_tpm.loc[:, meta["sample"]]
    meta = meta.set_index("sample")

    estimate = pd.read_csv(ESTIMATE_DEFAULT, sep="\t")
    estimate = estimate.set_index("ID")
    meta = meta.join(estimate[["Immune_score"]], how="left")

    lcam, lcam_present = lcam_scores(linear, signatures)
    trn_genes = signatures.loc[signatures["signature"] == "Salcher_TRN", "genes"].iloc[0].split(",")
    tan_genes = signatures.loc[signatures["signature"] == "Salcher_TAN", "genes"].iloc[0].split(",")
    nan_genes = signatures.loc[signatures["signature"] == "Salcher_NAN", "genes"].iloc[0].split(",")
    neut_genes = signatures.loc[
        signatures["signature"] == "Salcher_major_neutrophils", "genes"
    ].iloc[0].split(",")
    trn, trn_present = mean_z_score(linear, trn_genes)
    tan, tan_present = mean_z_score(linear, tan_genes)
    nan, nan_present = mean_z_score(linear, nan_genes)
    neut, neut_present = mean_z_score(linear, neut_genes)

    scores = meta.copy()
    scores["TACSTD2_log2tpm"] = log2_tpm.loc["TACSTD2"]
    scores["LCAM"] = lcam["LCAM"]
    scores["LCAMhi"] = lcam["LCAMhi"]
    scores["LCAMlo"] = lcam["LCAMlo"]
    scores["TRN"] = trn
    scores["TAN"] = tan
    scores["NAN"] = nan
    scores["major_neutrophils"] = neut
    scores.to_csv(output_dir / "tcga_primary_nsclc_scores.tsv", sep="\t")

    coverage_rows = [
        {"signature": "TACSTD2", "n_requested": 1, "n_present": 1, "present_genes": "TACSTD2", "missing_genes": ""},
        {
            "signature": "Salcher_TRN",
            "n_requested": len(trn_genes),
            "n_present": len(trn_present),
            "present_genes": ",".join(trn_present),
            "missing_genes": ",".join(g for g in trn_genes if g not in trn_present),
        },
        {
            "signature": "Salcher_TAN",
            "n_requested": len(tan_genes),
            "n_present": len(tan_present),
            "present_genes": ",".join(tan_present),
            "missing_genes": ",".join(g for g in tan_genes if g not in tan_present),
        },
        {
            "signature": "Salcher_NAN",
            "n_requested": len(nan_genes),
            "n_present": len(nan_present),
            "present_genes": ",".join(nan_present),
            "missing_genes": ",".join(g for g in nan_genes if g not in nan_present),
        },
        {
            "signature": "Salcher_major_neutrophils",
            "n_requested": len(neut_genes),
            "n_present": len(neut_present),
            "present_genes": ",".join(neut_present),
            "missing_genes": ",".join(g for g in neut_genes if g not in neut_present),
        },
    ]
    for subtype, genes in (
        (name, signatures.loc[signatures["signature"] == row, "genes"].iloc[0].split(","))
        for row, name in LCAM_SUBTYPE_ROWS.items()
    ):
        have = lcam_present[subtype]
        coverage_rows.append(
            {
                "signature": f"Leader_{subtype}",
                "n_requested": len(genes),
                "n_present": len(have),
                "present_genes": ",".join(have),
                "missing_genes": ",".join(g for g in genes if g not in have),
            }
        )
    pd.DataFrame(coverage_rows).to_csv(
        output_dir / "tcga_signature_coverage.tsv", sep="\t", index=False
    )

    pairs = [
        ("TACSTD2_log2tpm", "LCAM", "all_primary"),
        ("TACSTD2_log2tpm", "TRN", "all_primary"),
        ("TACSTD2_log2tpm", "TAN", "all_primary"),
        ("TACSTD2_log2tpm", "NAN", "all_primary"),
        ("TACSTD2_log2tpm", "major_neutrophils", "all_primary"),
        ("TACSTD2_log2tpm", "leukocyte_fraction", "all_primary"),
        ("LCAM", "TRN", "all_primary"),
        ("LCAM", "leukocyte_fraction", "all_primary"),
        ("TRN", "leukocyte_fraction", "all_primary"),
    ]
    corr_rows = []
    for left, right, subset in pairs:
        result = spearman(scores[left], scores[right])
        corr_rows.append({"subset": subset, "x": left, "y": right, **result})

    for disease in ["LUAD", "LUSC"]:
        part = scores[scores["disease"] == disease]
        for left, right in [
            ("TACSTD2_log2tpm", "LCAM"),
            ("TACSTD2_log2tpm", "TRN"),
            ("LCAM", "TRN"),
        ]:
            result = spearman(part[left], part[right])
            corr_rows.append({"subset": disease, "x": left, "y": right, **result})

    leuk = scores["leukocyte_fraction"].dropna()
    low_cut = leuk.quantile(0.25)
    high_cut = leuk.quantile(0.75)
    trop_cut = scores["TACSTD2_log2tpm"].quantile(0.75)
    immune_low = scores[scores["leukocyte_fraction"] <= low_cut]
    immune_high = scores[scores["leukocyte_fraction"] >= high_cut]
    trop_high = scores[scores["TACSTD2_log2tpm"] >= trop_cut]
    immune_low_trop_high = scores[
        (scores["leukocyte_fraction"] <= low_cut)
        & (scores["TACSTD2_log2tpm"] >= trop_cut)
    ]
    exclude_lowest10 = scores[
        scores["leukocyte_fraction"] >= scores["leukocyte_fraction"].quantile(0.10)
    ]

    for label, frame in [
        ("immune_low_Q1_leukocyte", immune_low),
        ("immune_high_Q4_leukocyte", immune_high),
        ("trop2_high_Q4", trop_high),
        ("immune_low_and_trop2_high", immune_low_trop_high),
        ("exclude_lowest10pct_leukocyte", exclude_lowest10),
    ]:
        for left, right in [
            ("TACSTD2_log2tpm", "LCAM"),
            ("TACSTD2_log2tpm", "TRN"),
            ("LCAM", "TRN"),
        ]:
            result = spearman(frame[left], frame[right])
            corr_rows.append({"subset": label, "x": left, "y": right, **result})

    pd.DataFrame(corr_rows).to_csv(
        output_dir / "tcga_spearman.tsv", sep="\t", index=False, float_format="%.6g"
    )

    scores = scores.copy()
    scores["leukocyte_quartile"] = pd.qcut(
        scores["leukocyte_fraction"], 4, labels=["Q1_low", "Q2", "Q3", "Q4_high"]
    )
    scores["tacstd2_quartile"] = pd.qcut(
        scores["TACSTD2_log2tpm"], 4, labels=["Q1_low", "Q2", "Q3", "Q4_high"]
    )
    strata = (
        scores.dropna(subset=["leukocyte_quartile"])
        .groupby(["leukocyte_quartile", "tacstd2_quartile"], observed=True)
        .agg(
            n=("TACSTD2_log2tpm", "size"),
            TACSTD2_log2tpm_mean=("TACSTD2_log2tpm", "mean"),
            LCAM_mean=("LCAM", "mean"),
            TRN_mean=("TRN", "mean"),
            leukocyte_fraction_mean=("leukocyte_fraction", "mean"),
        )
        .reset_index()
    )
    strata.to_csv(output_dir / "tcga_quartile_strata.tsv", sep="\t", index=False, float_format="%.6f")

    summary = pd.DataFrame(
        [
            {
                "stratum": "all_primary_LUAD_LUSC",
                "n": int(len(scores)),
                "n_LUAD": int((scores["disease"] == "LUAD").sum()),
                "n_LUSC": int((scores["disease"] == "LUSC").sum()),
                "leukocyte_Q1_cutoff": float(low_cut),
                "TACSTD2_Q4_cutoff": float(trop_cut),
                "n_immune_low_Q1": int(len(immune_low)),
                "n_trop2_high_Q4": int(len(trop_high)),
                "n_immune_low_and_trop2_high": int(len(immune_low_trop_high)),
            }
        ]
    )
    summary.to_csv(output_dir / "tcga_stratum_counts.tsv", sep="\t", index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--toil", type=Path, default=TOIL_DEFAULT)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    extract_leader_tables(args.output_dir)
    score_tcga(args.output_dir, args.toil)


if __name__ == "__main__":
    main()
