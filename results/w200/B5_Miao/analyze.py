#!/usr/bin/env python3
"""Miao 2018 RCC ICI: public RNA CLDN4 and TACSTD2 versus response.

Downloads the BHK Lab / ORCESTRA ICB_Miao1 TSV bundle (CC-BY-4.0), which is a
re-annotated public extract of Miao et al., Science 2018. Tests the two
prespecified genes against RECIST objective response, with batch, regimen,
and clinical-benefit sensitivities. Does not use the curated PFS event field
(see data-integrity checks).
"""

from __future__ import annotations

import hashlib
import json
import math
import urllib.request
import zipfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from scipy import stats


ROOT = Path(__file__).resolve().parent
CACHE = ROOT / ".cache"
ZIP_PATH = CACHE / "ICB_Miao1.zip"
EXTRACT = CACHE
SOURCE_URL = "https://zenodo.org/api/records/7058399/files/ICB_Miao1.zip/content"
SOURCE_SHA256 = "121c7a843aa3ceb95c8ec58354b92d131d998d56d029762ea8a3ab46bca07ffe"
SOURCE_MD5 = "9fdb0d89c6725410db09bad4d06a27f1"
SOURCE_DOI = "10.5281/zenodo.7058399"
GENES = {
    "CLDN4": "ENSG00000189143.8",
    "TACSTD2": "ENSG00000184292.5",
}
FLOOR = -9.96578428466209
SEED = 2018
N_BOOT = 20_000


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def obtain_bundle() -> None:
    CACHE.mkdir(parents=True, exist_ok=True)
    if ZIP_PATH.exists() and sha256(ZIP_PATH) == SOURCE_SHA256:
        pass
    else:
        temporary = ZIP_PATH.with_suffix(".download")
        print(f"Downloading {SOURCE_URL}")
        urllib.request.urlretrieve(SOURCE_URL, temporary)
        observed = sha256(temporary)
        if observed != SOURCE_SHA256:
            temporary.unlink(missing_ok=True)
            raise RuntimeError(
                f"Source checksum changed: expected {SOURCE_SHA256}, observed {observed}"
            )
        if md5(temporary) != SOURCE_MD5:
            temporary.unlink(missing_ok=True)
            raise RuntimeError("MD5 does not match the Zenodo record listing")
        temporary.replace(ZIP_PATH)
    EXTRACT.mkdir(parents=True, exist_ok=True)
    needed = [
        EXTRACT / "ICB_Miao1_metadata.tsv",
        EXTRACT / "ICB_Miao1_expr.tsv",
        EXTRACT / "ICB_Miao1_expr_genes.tsv",
    ]
    if not all(path.exists() for path in needed):
        with zipfile.ZipFile(ZIP_PATH) as archive:
            archive.extractall(CACHE)


def load_tables() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    metadata = pd.read_csv(EXTRACT / "ICB_Miao1_metadata.tsv", sep="\t")
    expression = pd.read_csv(EXTRACT / "ICB_Miao1_expr.tsv", sep="\t", index_col=0)
    genes = pd.read_csv(EXTRACT / "ICB_Miao1_expr_genes.tsv", sep="\t")
    if metadata.shape[0] != 52:
        raise RuntimeError(f"Expected 52 metadata rows; found {metadata.shape[0]}")
    if expression.shape != (57820, 33):
        raise RuntimeError(f"Unexpected expression shape {expression.shape}")
    # The TSV has a leading sample-ID field with no header; pandas uses it as index.
    if metadata.index.astype(str).duplicated().any():
        raise RuntimeError("Duplicate metadata sample IDs")
    return metadata, expression, genes


def analysis_frame(
    metadata: pd.DataFrame, expression: pd.DataFrame, genes: pd.DataFrame
) -> pd.DataFrame:
    rna = metadata.loc[metadata["rna"].notna()].copy()
    if len(rna) != 33:
        raise RuntimeError(f"Expected 33 RNA samples; found {len(rna)}")
    sample_ids = rna.index.astype(str)
    if set(sample_ids) != set(expression.columns.astype(str)):
        raise RuntimeError("RNA metadata IDs do not match expression columns")

    for symbol, gene_id in GENES.items():
        hits = genes.loc[genes["gene_name"] == symbol, "gene_id"]
        if list(hits) != [gene_id]:
            raise RuntimeError(f"{symbol} mapping changed: {list(hits)}")
        if gene_id not in expression.index:
            raise RuntimeError(f"{gene_id} missing from expression matrix")

    frame = pd.DataFrame(
        {
            "sample_id": sample_ids,
            "patientid": rna["patientid"].to_numpy(),
            "sex": rna["sex"].to_numpy(),
            "age": rna["age"].to_numpy(),
            "recist": rna["recist"].to_numpy(),
            "response_curated": rna["response"].to_numpy(),
            "response_category": rna["response_category"].to_numpy(),
            "drug": rna["drug"].to_numpy(),
            "first_line": rna["first_line"].to_numpy(),
            "survival_time_pfs": rna["survival_time_pfs"].to_numpy(),
            "event_occurred_pfs": rna["event_occurred_pfs"].to_numpy(),
            "survival_time_os": rna["survival_time_os"].to_numpy(),
            "event_occurred_os": rna["event_occurred_os"].to_numpy(),
        }
    )
    frame["id_prefix"] = np.where(frame["sample_id"].str.startswith("RCC_"), "RCC_", "PD_")
    frame["nivo_monotherapy"] = frame["drug"].eq("nivolumab")
    frame["orr_responder"] = np.where(
        frame["recist"].isin(["CR", "PR"]),
        "Responder",
        np.where(frame["recist"].isin(["SD", "PD"]), "Nonresponder", ""),
    )
    frame["cb_vs_ncb"] = np.where(
        frame["response_category"].eq("clinical benefit"),
        "CB",
        np.where(frame["response_category"].eq("no clinical benefit"), "NCB", ""),
    )
    for symbol, gene_id in GENES.items():
        values = expression.loc[gene_id, frame["sample_id"]].to_numpy(dtype=float)
        frame[symbol] = values
        frame[f"{symbol}_detected"] = values > (FLOOR + 1e-8)
        cohort_mean = frame.groupby("id_prefix")[symbol].transform("mean")
        frame[f"{symbol}_resid"] = values - cohort_mean
    if frame["sample_id"].duplicated().any():
        raise RuntimeError("Duplicate sample IDs")
    return frame.reset_index(drop=True)


def bootstrap_auc_ci(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    rng = np.random.default_rng(SEED)
    estimates = np.empty(N_BOOT)
    for i in range(N_BOOT):
        aa = rng.choice(a, len(a), replace=True)
        bb = rng.choice(b, len(b), replace=True)
        estimates[i] = stats.mannwhitneyu(aa, bb, alternative="two-sided").statistic / (
            len(aa) * len(bb)
        )
    return tuple(np.quantile(estimates, [0.025, 0.975]))


def bootstrap_median_diff_ci(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    rng = np.random.default_rng(SEED)
    estimates = np.empty(N_BOOT)
    for i in range(N_BOOT):
        aa = rng.choice(a, len(a), replace=True)
        bb = rng.choice(b, len(b), replace=True)
        estimates[i] = float(np.median(aa) - np.median(bb))
    return tuple(np.quantile(estimates, [0.025, 0.975]))


def compare_groups(high: np.ndarray, low: np.ndarray) -> dict:
    if len(high) == 0 or len(low) == 0:
        return {
            "n_high": int(len(high)),
            "n_low": int(len(low)),
            "median_high": math.nan,
            "median_low": math.nan,
            "median_diff": math.nan,
            "median_diff_ci_low": math.nan,
            "median_diff_ci_high": math.nan,
            "mann_whitney_u": math.nan,
            "mann_whitney_p": math.nan,
            "auc_high_higher": math.nan,
            "auc_bootstrap_95ci_low": math.nan,
            "auc_bootstrap_95ci_high": math.nan,
            "cliffs_delta": math.nan,
        }
    test = stats.mannwhitneyu(high, low, alternative="two-sided", method="asymptotic")
    auc = float(test.statistic / (len(high) * len(low)))
    auc_lo, auc_hi = bootstrap_auc_ci(high, low)
    med_lo, med_hi = bootstrap_median_diff_ci(high, low)
    return {
        "n_high": int(len(high)),
        "n_low": int(len(low)),
        "median_high": float(np.median(high)),
        "median_low": float(np.median(low)),
        "median_diff": float(np.median(high) - np.median(low)),
        "median_diff_ci_low": float(med_lo),
        "median_diff_ci_high": float(med_hi),
        "mann_whitney_u": float(test.statistic),
        "mann_whitney_p": float(test.pvalue),
        "auc_high_higher": auc,
        "auc_bootstrap_95ci_low": float(auc_lo),
        "auc_bootstrap_95ci_high": float(auc_hi),
        "cliffs_delta": float(2 * auc - 1),
    }


def run_contrast(frame: pd.DataFrame, gene: str, label_col: str, high: str, low: str, mask: pd.Series) -> dict:
    subset = frame.loc[mask & frame[label_col].isin([high, low])]
    values_high = subset.loc[subset[label_col].eq(high), gene].to_numpy()
    values_low = subset.loc[subset[label_col].eq(low), gene].to_numpy()
    result = compare_groups(values_high, values_low)
    if gene.endswith("_resid"):
        result["detected_high"] = math.nan
        result["detected_low"] = math.nan
        result["detection_fisher_p"] = math.nan
    else:
        detected_high = int(subset.loc[subset[label_col].eq(high), f"{gene}_detected"].sum())
        detected_low = int(subset.loc[subset[label_col].eq(low), f"{gene}_detected"].sum())
        table = [
            [detected_high, len(values_high) - detected_high],
            [detected_low, len(values_low) - detected_low],
        ]
        fisher = stats.fisher_exact(table, alternative="two-sided")
        result["detected_high"] = detected_high
        result["detected_low"] = detected_low
        result["detection_fisher_p"] = float(fisher.pvalue)
    result["n_total"] = int(len(subset))
    return result


def pca_batch(expression: pd.DataFrame, sample_ids: list[str]) -> dict:
    matrix = expression.loc[:, sample_ids].to_numpy(dtype=float)
    keep = (matrix > FLOOR).mean(axis=1) > 0.5
    kept = matrix[keep]
    top = np.argsort(kept.var(axis=1))[-2000:]
    centered = kept[top] - kept[top].mean(axis=1, keepdims=True)
    _u, singular, vt = np.linalg.svd(centered, full_matrices=False)
    explained = (singular**2) / (singular**2).sum()
    pc1 = vt[0]
    prefix = np.where(pd.Index(sample_ids).str.startswith("RCC_"), "RCC_", "PD_")
    a = pc1[prefix == "PD_"]
    b = pc1[prefix == "RCC_"]
    return {
        "genes_used": int(keep.sum()),
        "pc1_variance_fraction": float(explained[0]),
        "pc2_variance_fraction": float(explained[1]),
        "pc1_pd_mean": float(a.mean()),
        "pc1_rcc_mean": float(b.mean()),
        "pc1_prefix_mw_p": float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue),
    }


def cox_os(frame: pd.DataFrame, gene: str, adjust_prefix: bool) -> dict:
    cols = [gene, "survival_time_os", "event_occurred_os"]
    data = frame[cols].copy()
    if adjust_prefix:
        data["prefix_rcc"] = frame["id_prefix"].eq("RCC_").astype(int)
    model = CoxPHFitter()
    model.fit(data, duration_col="survival_time_os", event_col="event_occurred_os")
    summary = model.summary.loc[gene]
    return {
        "n": int(len(data)),
        "n_events": int(frame["event_occurred_os"].sum()),
        "hr_per_log2_unit": float(summary["exp(coef)"]),
        "hr_ci_low": float(summary["exp(coef) lower 95%"]),
        "hr_ci_high": float(summary["exp(coef) upper 95%"]),
        "p": float(summary["p"]),
        "adjusted_for_id_prefix": bool(adjust_prefix),
    }


def make_figure(frame: pd.DataFrame, summaries: pd.DataFrame) -> None:
    plt.style.use("seaborn-v0_8-whitegrid")
    fig, axes = plt.subplots(2, 3, figsize=(11.6, 7.2), sharey=False)
    rng = np.random.default_rng(SEED)
    colors = ["#6b7280", "#2563eb"]
    panels = [
        ("pooled", "All ICI (primary)", pd.Series(True, index=frame.index)),
        ("PD_", "ID prefix PD_", frame["id_prefix"].eq("PD_")),
        ("RCC_", "ID prefix RCC_", frame["id_prefix"].eq("RCC_")),
    ]
    pmap = {
        (row.gene, row.stratum): row.mann_whitney_p
        for row in summaries.itertuples()
        if row.endpoint == "RECIST_ORR" and row.expression == "raw"
    }
    for row_i, gene in enumerate(("CLDN4", "TACSTD2")):
        for col_i, (stratum, title, mask) in enumerate(panels):
            ax = axes[row_i, col_i]
            subset = frame.loc[mask & frame["orr_responder"].ne("")]
            groups = [
                subset.loc[subset["orr_responder"].eq("Nonresponder"), gene].to_numpy(),
                subset.loc[subset["orr_responder"].eq("Responder"), gene].to_numpy(),
            ]
            bp = ax.boxplot(
                groups,
                positions=[1, 2],
                widths=0.55,
                patch_artist=True,
                showfliers=False,
                medianprops={"color": "black", "linewidth": 1.5},
            )
            for box, color in zip(bp["boxes"], colors, strict=True):
                box.set_facecolor(color)
                box.set_alpha(0.35)
            for pos, values, color in zip([1, 2], groups, colors, strict=True):
                jitter = rng.uniform(-0.12, 0.12, len(values))
                ax.scatter(
                    pos + jitter,
                    values,
                    s=22,
                    color=color,
                    alpha=0.85,
                    edgecolor="white",
                    linewidth=0.35,
                    zorder=3,
                )
            ax.set_xticks([1, 2], [f"NR\nn={len(groups[0])}", f"R\nn={len(groups[1])}"])
            ax.set_title(f"{gene}: {title}", fontsize=10, weight="bold")
            if col_i == 0:
                ax.set_ylabel("log2(TPM), floor ≈ log2(0.001)")
            pvalue = pmap.get((gene, stratum), math.nan)
            ax.text(
                0.5,
                0.98,
                f"MW p={pvalue:.3g}" if pd.notna(pvalue) else "MW p=NA",
                transform=ax.transAxes,
                ha="center",
                va="top",
                fontsize=9,
            )
    fig.suptitle(
        "Miao 2018 RCC ICI: CLDN4 and TACSTD2 versus RECIST CR/PR vs SD/PD",
        y=1.02,
    )
    fig.tight_layout()
    fig.savefig(ROOT / "cldn4_tacstd2_vs_response.png", dpi=220, bbox_inches="tight")
    fig.savefig(ROOT / "cldn4_tacstd2_vs_response.svg", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    obtain_bundle()
    metadata, expression, genes = load_tables()
    frame = analysis_frame(metadata, expression, genes)
    batch = pca_batch(expression, list(frame["sample_id"]))

    pfs_pd = frame.loc[frame["recist"].eq("PD"), "event_occurred_pfs"]
    integrity = {
        "n_rna": int(len(frame)),
        "n_recist_pd": int(frame["recist"].eq("PD").sum()),
        "n_recist_pd_pfs_event": int(pfs_pd.sum()),
        "pfs_event_rate": float(frame["event_occurred_pfs"].mean()),
        "os_event_rate": float(frame["event_occurred_os"].mean()),
        "pfs_field_unusable": bool(pfs_pd.sum() == 0 and (frame["recist"].eq("PD")).any()),
        "orr_responders_by_drug": frame.loc[frame["orr_responder"].eq("Responder"), "drug"]
        .value_counts()
        .to_dict(),
        "pca": batch,
        "cldn4_tacstd2_spearman": {
            "rho": float(stats.spearmanr(frame["CLDN4"], frame["TACSTD2"]).statistic),
            "p": float(stats.spearmanr(frame["CLDN4"], frame["TACSTD2"]).pvalue),
        },
    }

    contrasts = []
    contrast_specs = [
        ("RECIST_ORR", "orr_responder", "Responder", "Nonresponder", "raw"),
        ("curated_R_vs_NR", "response_curated", "R", "NR", "raw"),
        ("CB_vs_NCB", "cb_vs_ncb", "CB", "NCB", "raw"),
    ]
    strata = [
        ("pooled", pd.Series(True, index=frame.index)),
        ("PD_", frame["id_prefix"].eq("PD_")),
        ("RCC_", frame["id_prefix"].eq("RCC_")),
        ("nivo_monotherapy", frame["nivo_monotherapy"]),
    ]
    for endpoint, column, high, low, _scale in contrast_specs:
        for gene in GENES:
            for stratum, mask in strata:
                row = run_contrast(frame, gene, column, high, low, mask)
                row.update(
                    {
                        "gene": gene,
                        "endpoint": endpoint,
                        "stratum": stratum,
                        "expression": "raw",
                        "high_label": high,
                        "low_label": low,
                    }
                )
                contrasts.append(row)
            residual = run_contrast(
                frame, f"{gene}_resid", column, high, low, pd.Series(True, index=frame.index)
            )
            residual.update(
                {
                    "gene": gene,
                    "endpoint": endpoint,
                    "stratum": "pooled_prefix_residualized",
                    "expression": "prefix_residual",
                    "high_label": high,
                    "low_label": low,
                }
            )
            contrasts.append(residual)

    summary = pd.DataFrame(contrasts)
    cox_rows = []
    for gene in GENES:
        for adjust in (False, True):
            row = cox_os(frame, gene, adjust)
            row["gene"] = gene
            cox_rows.append(row)
    cox = pd.DataFrame(cox_rows)

    frame.to_csv(ROOT / "sample_level.csv", index=False)
    summary.to_csv(ROOT / "summary.csv", index=False)
    cox.to_csv(ROOT / "os_cox.csv", index=False)
    make_figure(frame, summary)

    provenance = pd.DataFrame(
        [
            {
                "file": ZIP_PATH.name,
                "sha256": sha256(ZIP_PATH),
                "md5": md5(ZIP_PATH),
                "url": SOURCE_URL,
                "doi": SOURCE_DOI,
                "license": "CC-BY-4.0",
            }
        ]
    )
    provenance.to_csv(ROOT / "provenance.tsv", sep="\t", index=False)

    payload = {
        "source": {
            "paper": "Miao et al. Science 2018;359:801-806",
            "pubmed": "29301960",
            "curation": "BHK Lab ICB_Miao1 TSV extract, Zenodo 7058399",
            "license": "CC-BY-4.0",
        },
        "integrity": integrity,
        "primary": summary.loc[
            (summary["endpoint"] == "RECIST_ORR")
            & (summary["stratum"] == "pooled")
            & (summary["expression"] == "raw")
        ].to_dict(orient="records"),
        "os_cox": cox.to_dict(orient="records"),
    }
    (ROOT / "results.json").write_text(json.dumps(payload, indent=2) + "\n")

    print("RNA n=", len(frame))
    print("RECIST", frame["recist"].value_counts().to_dict())
    print("ORR R/NR", frame["orr_responder"].value_counts().to_dict())
    print("PFS events among RECIST PD:", int(pfs_pd.sum()), "/", int(len(pfs_pd)))
    print("PC1 variance", round(batch["pc1_variance_fraction"], 3), "prefix MW p", f"{batch['pc1_prefix_mw_p']:.2e}")
    for row in payload["primary"]:
        print(
            f"PRIMARY {row['gene']}: R={row['n_high']} NR={row['n_low']} "
            f"AUC={row['auc_high_higher']:.3f} MW p={row['mann_whitney_p']:.4g}"
        )


if __name__ == "__main__":
    main()
