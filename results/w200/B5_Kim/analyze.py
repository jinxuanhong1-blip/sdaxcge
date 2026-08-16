#!/usr/bin/env python3
"""Kim 2018 gastric pembrolizumab: public CLDN4 versus ICI response.

B5 analog of results/w200/B5_Riaz. Uses the open ORCESTRA/PredictIO ICB_Kim
TSV release (Zenodo 10.5281/zenodo.7058399), not FASTQ.
"""

from __future__ import annotations

import csv
import hashlib
import math
import zipfile
from pathlib import Path
from urllib.request import urlretrieve

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
ZIP_PATH = DATA / "ICB_Kim.zip"
META_PATH = DATA / "ICB_Kim_metadata.tsv"
EXPR_PATH = DATA / "ICB_Kim_expr_gene_tpm.tsv"
GENE_PATH = DATA / "ICB_Kim_expr_gene_tpm_genes.tsv"
ZIP_URL = "https://zenodo.org/api/records/7058399/files/ICB_Kim.zip/content"
ZIP_MD5 = "8ec1669e80c0f365172982ff39f8f6ca"
NEEDED = (
    "ICB_Kim_metadata.tsv",
    "ICB_Kim_expr_gene_tpm.tsv",
    "ICB_Kim_expr_gene_tpm_genes.tsv",
)
GENES = {
    "CLDN4": "ENSG00000189143.10",
    "TACSTD2": "ENSG00000184292.7",
    "CD274": "ENSG00000120217.14",
    "CD8A": "ENSG00000153563.16",
    "CXCL9": "ENSG00000138755.6",
    "GZMB": "ENSG00000100453.14",
    "EPCAM": "ENSG00000119888.11",
}
SEED = 25780
BOOTSTRAP_DRAWS = 20_000
EXPECTED_RECIST = {"CR": 13, "SD": 14, "PD": 18}


def md5sum(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256sum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_inputs() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    if not ZIP_PATH.exists():
        print(f"Downloading {ZIP_URL}")
        urlretrieve(ZIP_URL, ZIP_PATH)
    got = md5sum(ZIP_PATH)
    if got != ZIP_MD5:
        raise ValueError(f"ICB_Kim.zip MD5 mismatch: {got} != {ZIP_MD5}")
    with zipfile.ZipFile(ZIP_PATH) as archive:
        for name in NEEDED:
            dest = DATA / name
            if not dest.exists():
                archive.extract(name, DATA)


def read_r_tsv(path: Path) -> list[dict[str, str]]:
    """Read an R write.table TSV that has a rownames column and no header for it."""
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = [item.strip('"') for item in next(reader)]
        rows = []
        for raw in reader:
            raw = [item.strip('"') for item in raw]
            if len(raw) == len(header) + 1:
                raw = raw[1:]
            if len(raw) != len(header):
                raise ValueError(
                    f"{path.name}: expected {len(header)} fields, found {len(raw)}"
                )
            rows.append(dict(zip(header, raw, strict=True)))
    return rows


def read_metadata() -> list[dict[str, str]]:
    rows = read_r_tsv(META_PATH)
    if len(rows) != 45:
        raise ValueError(f"Expected 45 Kim RNA patients; found {len(rows)}")
    counts: dict[str, int] = {}
    for row in rows:
        recist = row["recist"]
        counts[recist] = counts.get(recist, 0) + 1
    if counts != EXPECTED_RECIST:
        raise ValueError(f"Unexpected RECIST counts: {counts}")
    return rows


def read_expression() -> tuple[list[str], dict[str, list[float]]]:
    with EXPR_PATH.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        samples = [item.strip('"') for item in next(reader)]
        matrix: dict[str, list[float]] = {}
        for row in reader:
            gene_id = row[0].strip('"')
            matrix[gene_id] = [float(value) for value in row[1:]]
    if len(samples) != 45:
        raise ValueError(f"Expected 45 expression columns; found {len(samples)}")
    return samples, matrix


def read_gene_names() -> dict[str, str]:
    mapping = {row["gene_name"]: row["gene_id"] for row in read_r_tsv(GENE_PATH)}
    for symbol, gene_id in GENES.items():
        if mapping.get(symbol) != gene_id:
            raise ValueError(f"{symbol} gene_id mismatch: {mapping.get(symbol)}")
    return mapping


def bootstrap_auc_ci(
    responders: np.ndarray, nonresponders: np.ndarray, draws: int = BOOTSTRAP_DRAWS
) -> tuple[float, float]:
    rng = np.random.default_rng(SEED)
    estimates = np.empty(draws)
    for i in range(draws):
        a = rng.choice(responders, len(responders), replace=True)
        b = rng.choice(nonresponders, len(nonresponders), replace=True)
        estimates[i] = stats.mannwhitneyu(a, b).statistic / (len(a) * len(b))
    return tuple(np.quantile(estimates, [0.025, 0.975]))


def woolf_or_ci(table: np.ndarray) -> tuple[float, float, float]:
    a, b, c, d = (float(x) for x in table.ravel())
    if min(a, b, c, d) == 0:
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
        corrected = True
    else:
        corrected = False
    log_or = math.log((a * d) / (b * c))
    se = math.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    lo = math.exp(log_or - 1.96 * se)
    hi = math.exp(log_or + 1.96 * se)
    point = math.exp(log_or)
    if corrected:
        return point, lo, hi
    return point, lo, hi


def fisher_or(table: np.ndarray) -> tuple[float, float]:
    result = stats.fisher_exact(table, alternative="two-sided")
    return float(result.statistic), float(result.pvalue)


def median_split_table(values: np.ndarray, responder: np.ndarray) -> np.ndarray:
    high = values >= np.median(values)
    return np.array(
        [
            [int(np.sum(high & responder)), int(np.sum(high & ~responder))],
            [int(np.sum(~high & responder)), int(np.sum(~high & ~responder))],
        ],
        dtype=int,
    )


def summarize_continuous(responders: np.ndarray, nonresponders: np.ndarray) -> dict:
    test = stats.mannwhitneyu(
        responders, nonresponders, alternative="two-sided", method="asymptotic"
    )
    auc = float(test.statistic / (len(responders) * len(nonresponders)))
    auc_lo, auc_hi = bootstrap_auc_ci(responders, nonresponders)
    return {
        "n_responder": int(len(responders)),
        "n_nonresponder": int(len(nonresponders)),
        "n_total": int(len(responders) + len(nonresponders)),
        "median_responder": float(np.median(responders)),
        "median_nonresponder": float(np.median(nonresponders)),
        "mann_whitney_u": float(test.statistic),
        "mann_whitney_p": float(test.pvalue),
        "auc_responder_higher": auc,
        "auc_bootstrap_95ci_low": float(auc_lo),
        "auc_bootstrap_95ci_high": float(auc_hi),
    }


def summarize_or(values: np.ndarray, responder: np.ndarray) -> dict:
    table = median_split_table(values, responder)
    or_hat, p_value = fisher_or(table)
    or_woolf, lo, hi = woolf_or_ci(table.astype(float))
    return {
        "n_total": int(len(values)),
        "n_responder": int(np.sum(responder)),
        "n_nonresponder": int(np.sum(~responder)),
        "n_high": int(table[0].sum()),
        "n_low": int(table[1].sum()),
        "high_responder": int(table[0, 0]),
        "high_nonresponder": int(table[0, 1]),
        "low_responder": int(table[1, 0]),
        "low_nonresponder": int(table[1, 1]),
        "split_median": float(np.median(values)),
        "odds_ratio_fisher": or_hat,
        "odds_ratio_woolf": or_woolf,
        "or_95ci_low": lo,
        "or_95ci_high": hi,
        "fisher_p": p_value,
    }


def write_csv(path: Path, rows: list[dict]) -> None:
    fieldnames: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for key in row:
            if key not in seen:
                seen.add(key)
                fieldnames.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def make_figure(
    sample_rows: list[dict],
    summaries: dict[str, dict],
) -> None:
    plt.style.use("seaborn-v0_8-whitegrid")
    fig, axes = plt.subplots(1, 3, figsize=(12.2, 4.2))
    rng = np.random.default_rng(SEED)

    def values_for(gene: str, recists: set[str]) -> np.ndarray:
        return np.array(
            [
                float(row[gene])
                for row in sample_rows
                if row["recist"] in recists
            ]
        )

    panels = [
        (
            "CLDN4, CR vs SD+PD (primary)",
            [values_for("CLDN4", {"SD", "PD"}), values_for("CLDN4", {"CR"})],
            "log2 TPM",
            summaries["cldn4_orr_mw"]["mann_whitney_p"],
            (
                f"OR {summaries['cldn4_orr_or']['odds_ratio_fisher']:.2f}, "
                f"p={summaries['cldn4_orr_or']['fisher_p']:.2g}, "
                f"n={summaries['cldn4_orr_or']['n_total']}"
            ),
        ),
        (
            "CLDN4, CR vs PD (exclude SD)",
            [values_for("CLDN4", {"PD"}), values_for("CLDN4", {"CR"})],
            "log2 TPM",
            summaries["cldn4_cr_pd_mw"]["mann_whitney_p"],
            (
                f"OR {summaries['cldn4_cr_pd_or']['odds_ratio_fisher']:.2f}, "
                f"p={summaries['cldn4_cr_pd_or']['fisher_p']:.2g}, "
                f"n={summaries['cldn4_cr_pd_or']['n_total']}"
            ),
        ),
        (
            "CD274 control, CR vs SD+PD",
            [values_for("CD274", {"SD", "PD"}), values_for("CD274", {"CR"})],
            "log2 TPM",
            summaries["cd274_orr_mw"]["mann_whitney_p"],
            (
                f"OR {summaries['cd274_orr_or']['odds_ratio_fisher']:.2f}, "
                f"p={summaries['cd274_orr_or']['fisher_p']:.2g}, "
                f"n={summaries['cd274_orr_or']['n_total']}"
            ),
        ),
    ]
    colors = ["#6b7280", "#2563eb"]
    for ax, (title, groups, ylabel, p_value, extra) in zip(axes, panels, strict=True):
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
                s=23,
                color=color,
                alpha=0.85,
                edgecolor="white",
                linewidth=0.35,
                zorder=3,
            )
        ax.set_xticks([1, 2], [f"NR\nn={len(groups[0])}", f"R\nn={len(groups[1])}"])
        ax.set_title(title, fontsize=9.5, weight="bold")
        ax.set_ylabel(ylabel)
        ax.text(
            0.5,
            0.98,
            f"MW p={p_value:.3g}\n{extra}",
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=8.5,
        )
    fig.suptitle(
        "Kim 2018 gastric pembrolizumab: CLDN4 versus ICI response",
        y=1.03,
    )
    fig.tight_layout()
    fig.savefig(ROOT / "cldn4_vs_response.png", dpi=220, bbox_inches="tight")
    fig.savefig(ROOT / "cldn4_vs_response.svg", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    download_inputs()
    metadata = read_metadata()
    samples, matrix = read_expression()
    read_gene_names()
    meta_ids = [row["patientid"] for row in metadata]
    if samples != meta_ids:
        raise ValueError("Expression columns and metadata patient IDs differ")

    sample_rows = []
    for idx, row in enumerate(metadata):
        out = {
            "patientid": row["patientid"],
            "recist": row["recist"],
            "response": row["response"],
            "orr_responder": "yes" if row["recist"] == "CR" else "no",
        }
        for symbol, gene_id in GENES.items():
            out[symbol] = matrix[gene_id][idx]
        sample_rows.append(out)

    def gene_split(symbol: str, responder_labels: set[str], allowed: set[str]):
        values = np.array(
            [
                float(row[symbol])
                for row in sample_rows
                if row["recist"] in allowed
            ]
        )
        responder = np.array(
            [
                row["recist"] in responder_labels
                for row in sample_rows
                if row["recist"] in allowed
            ]
        )
        return values[responder], values[~responder], values, responder

    cldn4_r, cldn4_nr, cldn4_all, cldn4_is_r = gene_split(
        "CLDN4", {"CR"}, {"CR", "SD", "PD"}
    )
    cldn4_cr, cldn4_pd, cldn4_crpd, cldn4_crpd_r = gene_split(
        "CLDN4", {"CR"}, {"CR", "PD"}
    )
    cldn4_dcb, cldn4_pd_only, cldn4_dcb_all, cldn4_dcb_r = gene_split(
        "CLDN4", {"CR", "SD"}, {"CR", "SD", "PD"}
    )
    cd274_r, cd274_nr, cd274_all, cd274_is_r = gene_split(
        "CD274", {"CR"}, {"CR", "SD", "PD"}
    )
    tac_r, tac_nr, tac_all, tac_is_r = gene_split(
        "TACSTD2", {"CR"}, {"CR", "SD", "PD"}
    )
    epcam_r, epcam_nr, epcam_all, epcam_is_r = gene_split(
        "EPCAM", {"CR"}, {"CR", "SD", "PD"}
    )

    summaries = {
        "cldn4_orr_or": summarize_or(cldn4_all, cldn4_is_r),
        "cldn4_orr_mw": summarize_continuous(cldn4_r, cldn4_nr),
        "cldn4_cr_pd_or": summarize_or(cldn4_crpd, cldn4_crpd_r),
        "cldn4_cr_pd_mw": summarize_continuous(cldn4_cr, cldn4_pd),
        "cldn4_dcb_or": summarize_or(cldn4_dcb_all, cldn4_dcb_r),
        "cldn4_dcb_mw": summarize_continuous(cldn4_dcb, cldn4_pd_only),
        "cd274_orr_or": summarize_or(cd274_all, cd274_is_r),
        "cd274_orr_mw": summarize_continuous(cd274_r, cd274_nr),
        "tacstd2_orr_or": summarize_or(tac_all, tac_is_r),
        "tacstd2_orr_mw": summarize_continuous(tac_r, tac_nr),
        "epcam_orr_or": summarize_or(epcam_all, epcam_is_r),
        "epcam_orr_mw": summarize_continuous(epcam_r, epcam_nr),
    }

    analysis_meta = {
        "cldn4_orr_or": ("CLDN4", "CR vs SD+PD"),
        "cldn4_orr_mw": ("CLDN4", "CR vs SD+PD"),
        "cldn4_cr_pd_or": ("CLDN4", "CR vs PD"),
        "cldn4_cr_pd_mw": ("CLDN4", "CR vs PD"),
        "cldn4_dcb_or": ("CLDN4", "CR+SD vs PD"),
        "cldn4_dcb_mw": ("CLDN4", "CR+SD vs PD"),
        "cd274_orr_or": ("CD274", "CR vs SD+PD"),
        "cd274_orr_mw": ("CD274", "CR vs SD+PD"),
        "tacstd2_orr_or": ("TACSTD2", "CR vs SD+PD"),
        "tacstd2_orr_mw": ("TACSTD2", "CR vs SD+PD"),
        "epcam_orr_or": ("EPCAM", "CR vs SD+PD"),
        "epcam_orr_mw": ("EPCAM", "CR vs SD+PD"),
    }
    summary_rows = []
    for analysis, values in summaries.items():
        gene, contrast_name = analysis_meta[analysis]
        row = {"analysis": analysis, "gene": gene, "contrast": contrast_name, **values}
        summary_rows.append(row)

    write_csv(ROOT / "summary.csv", summary_rows)
    write_csv(ROOT / "sample_level.csv", sample_rows)
    make_figure(sample_rows, summaries)

    or_rows = [
        {
            "gene": row["gene"],
            "contrast": row["contrast"],
            "odds_ratio": row["odds_ratio_fisher"],
            "or_95ci_low": row["or_95ci_low"],
            "or_95ci_high": row["or_95ci_high"],
            "n": row["n_total"],
            "n_responder": row["n_responder"],
            "n_nonresponder": row["n_nonresponder"],
            "high_R": row["high_responder"],
            "high_NR": row["high_nonresponder"],
            "low_R": row["low_responder"],
            "low_NR": row["low_nonresponder"],
            "fisher_p": row["fisher_p"],
        }
        for row in summary_rows
        if row["analysis"].endswith("_or")
    ]
    write_csv(ROOT / "or_n_p.csv", or_rows)

    with (ROOT / "provenance.tsv").open("w") as handle:
        handle.write("file\tmd5\tsha256\turl_or_note\n")
        handle.write(
            f"{ZIP_PATH.name}\t{md5sum(ZIP_PATH)}\t{sha256sum(ZIP_PATH)}\t{ZIP_URL}\n"
        )
        for path in (META_PATH, EXPR_PATH, GENE_PATH):
            handle.write(
                f"{path.name}\t{md5sum(path)}\t{sha256sum(path)}\t"
                "extracted from ICB_Kim.zip (Zenodo 10.5281/zenodo.7058399)\n"
            )

    primary = summaries["cldn4_orr_or"]
    print(
        "PRIMARY CLDN4 median-split CR vs SD+PD: "
        f"OR={primary['odds_ratio_fisher']:.3f} "
        f"n={primary['n_total']} "
        f"p={primary['fisher_p']:.4g} "
        f"table={primary['high_responder']}/{primary['high_nonresponder']} vs "
        f"{primary['low_responder']}/{primary['low_nonresponder']}"
    )
    mw = summaries["cldn4_orr_mw"]
    print(
        "SUPPORT CLDN4 MW: "
        f"AUC={mw['auc_responder_higher']:.3f} "
        f"p={mw['mann_whitney_p']:.4g} "
        f"nR={mw['n_responder']} nNR={mw['n_nonresponder']}"
    )
    control = summaries["cd274_orr_or"]
    print(
        "CONTROL CD274 median-split: "
        f"OR={control['odds_ratio_fisher']:.3f} "
        f"n={control['n_total']} "
        f"p={control['fisher_p']:.4g}"
    )


if __name__ == "__main__":
    main()
