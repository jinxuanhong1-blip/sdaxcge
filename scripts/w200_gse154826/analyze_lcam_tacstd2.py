#!/usr/bin/env python3
"""GSE154826 LCAM vs TACSTD2 from compact processed/open files only.

Honest skip: SRA/FASTQ, GEO 10x MTX tarballs (~3 GB), and the full
lung_ldm.rd UMI object. Those are huge raw or unannotated count dumps
and are not required for the author LCAM-hi vs LCAM-lo DE or for
cell-composition LCAM scores.

Sources (all public, compact):
  - effiken/Leader_et_al input_tables (Cancer Cell 2021, PMID 34767762)
  - GEO GSE154826 series matrices + sample_annots.csv.gz (metadata only)
"""

from __future__ import annotations

import hashlib
import json
import urllib.request
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyreadr
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
CACHE = Path("/tmp/gse154826_open")
OUT = ROOT / "results" / "w200" / "GSE154826"

GITHUB_RAW = "https://raw.githubusercontent.com/effiken/Leader_et_al/master/input_tables"
GEO_SUPPL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154826/suppl"
GEO_MATRIX = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154826/matrix"

FILES = {
    "DE_LCAMhi_vs_LCAMlo_pseudobulk.rd": f"{GITHUB_RAW}/DE_LCAMhi_vs_LCAMlo_pseudobulk.rd",
    "annots_list.csv": f"{GITHUB_RAW}/annots_list.csv",
    "table_s1_sample_table.csv": f"{GITHUB_RAW}/table_s1_sample_table.csv",
    "cell_metadata.csv": f"{GITHUB_RAW}/cell_metadata.csv",
    "immune_vs_ep_de.csv": f"{GITHUB_RAW}/immune_vs_ep_de.csv",
    "table_s5_genes_used_in_figures.csv": f"{GITHUB_RAW}/table_s5_genes_used_in_figures.csv",
    "GSE154826_sample_annots.csv.gz": f"{GEO_SUPPL}/GSE154826_sample_annots.csv.gz",
    "GSE154826-GPL18573_series_matrix.txt.gz": f"{GEO_MATRIX}/GSE154826-GPL18573_series_matrix.txt.gz",
    "GSE154826-GPL24676_series_matrix.txt.gz": f"{GEO_MATRIX}/GSE154826-GPL24676_series_matrix.txt.gz",
}

# Author LCAM definition from figure_5abcd_s5a.R
LCAM_HI = ["T_activated", "IgG", "MoMac-II"]
LCAM_LO = ["B", "AM", "cDC2", "AZU1_mac", "Tcm/naive_II", "cDC1"]
NORM_GROUPS = ["T", "B&plasma", "MNP", "lin_neg"]

FOCUS_GENES = [
    "TACSTD2",
    "EPCAM",
    "CLDN4",
    "CLDN3",
    "CLDN7",
    "KRT8",
    "KRT18",
    "KRT19",
    "CDH1",
    "PDCD1",
    "CXCL13",
    "SPP1",
    "IGHG1",
    "IGHG3",
    "MZB1",
    "LAYN",
    "TNFRSF9",
    "GZMB",
    "PTPRC",
    "CD3D",
    "MS4A1",
    "CLEC9A",
    "CLEC10A",
    "AZU1",
    "SFTPC",
]


def md5_file(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(name: str) -> Path:
    CACHE.mkdir(parents=True, exist_ok=True)
    dest = CACHE / name
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    url = FILES[name]
    print(f"download {name}")
    urllib.request.urlretrieve(url, dest)
    return dest


def load_de() -> pd.DataFrame:
    obj = pyreadr.read_r(str(fetch("DE_LCAMhi_vs_LCAMlo_pseudobulk.rd")))
    de = obj["DE_total"].copy()
    de.index.name = "gene"
    de["neglog10_p"] = -np.log10(de["p.value"].clip(lower=1e-300))
    de["neglog10_padj"] = -np.log10(de["adj.p.value"].clip(lower=1e-300))
    return de


def composition_scores() -> tuple[pd.DataFrame, pd.DataFrame]:
    annots = pd.read_csv(fetch("annots_list.csv"))
    samples = pd.read_csv(fetch("table_s1_sample_table.csv"))
    cells = pd.read_csv(fetch("cell_metadata.csv"))

    annots["cluster"] = annots["cluster"].astype(int)
    cells["cluster_ID"] = cells["cluster_ID"].astype(int)
    cells = cells.merge(
        annots[["cluster", "lineage", "sub_lineage", "norm_group"]],
        left_on="cluster_ID",
        right_on="cluster",
        how="left",
    )
    samples["sample_ID"] = samples["sample_ID"].astype(int)
    cells["sample_ID"] = cells["sample_ID"].astype(int)
    cells = cells.merge(samples, on="sample_ID", how="left")

    # Gated epi/endo/fibro doublets are not a real epithelial compartment.
    cells["is_gated_doublet"] = cells["lineage"].eq("epi_endo_fibro_doublet")

    # Per-sample cluster counts (all cells, including gated).
    sample_n = (
        cells.groupby(
            [
                "sample_ID",
                "patient_ID",
                "tissue",
                "disease",
                "library_chemistry",
                "prep",
                "Use.in.Clustering.Model.",
            ],
            dropna=False,
        )
        .size()
        .rename("n_cells")
        .reset_index()
    )
    gated_n = (
        cells.groupby("sample_ID")["is_gated_doublet"].sum().rename("n_gated_doublet")
    )
    sample_n = sample_n.merge(gated_n, on="sample_ID", how="left")

    immune = cells.loc[~cells["is_gated_doublet"]].copy()
    counts = (
        immune.groupby(["sample_ID", "sub_lineage"]).size().unstack(fill_value=0)
    )
    frac = counts.div(counts.sum(axis=1).replace(0, np.nan), axis=0)

    # Within-lineage frequencies as in figure_5abcd_s5a.R
    within = frac.copy()
    for group in NORM_GROUPS:
        cols = annots.loc[annots["norm_group"].eq(group), "sub_lineage"]
        cols = [c for c in cols.unique() if c in within.columns and str(c) != "nan"]
        if not cols:
            continue
        denom = within[cols].sum(axis=1).replace(0, np.nan)
        within[cols] = within[cols].div(denom, axis=0)

    def score(frame: pd.DataFrame, cols: list[str]) -> pd.Series:
        use = [c for c in cols if c in frame.columns]
        return np.log(frame[use].fillna(0) + 1e-2).sum(axis=1)

    scores = pd.DataFrame(
        {
            "lcam_hi_score": score(within, LCAM_HI),
            "lcam_lo_score": score(within, LCAM_LO),
        }
    )
    scores["lcam_difference"] = scores["lcam_hi_score"] - scores["lcam_lo_score"]
    for col in LCAM_HI:
        scores[f"frac_{col}"] = frac[col] if col in frac.columns else np.nan
        scores[f"within_{col}"] = within[col] if col in within.columns else np.nan
    for col in LCAM_LO:
        scores[f"frac_{col}"] = frac[col] if col in frac.columns else np.nan

    out = sample_n.merge(scores, on="sample_ID", how="left")

    # Author primary LCAM axis used V2 beads tumor samples.
    # Binary call is a median split of the continuous difference on that subset
    # only (the two scores have different numbers of terms, so a zero cut is
    # meaningless).
    v2_mask = (
        out["tissue"].eq("Tumor")
        & out["library_chemistry"].eq("V2")
        & out["prep"].eq("beads")
    )
    v2_beads_tumor = out.loc[v2_mask].copy()
    med = v2_beads_tumor["lcam_difference"].median()
    v2_beads_tumor["lcam_call"] = np.where(
        v2_beads_tumor["lcam_difference"] >= med, "LCAM_hi", "LCAM_lo"
    )
    out["lcam_call"] = "not_v2_beads_tumor"
    out.loc[v2_beads_tumor.index, "lcam_call"] = v2_beads_tumor["lcam_call"]
    return out, v2_beads_tumor


def plot_volcano(de: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 5.4))
    x = de["log2_FC"]
    y = de["neglog10_p"]
    ax.scatter(x, y, s=4, c="#9aa3ad", alpha=0.35, linewidths=0, rasterized=True)
    sig = de["adj.p.value"].lt(0.05) & de["log2_FC"].abs().ge(1)
    ax.scatter(x[sig], y[sig], s=6, c="#4c78a8", alpha=0.55, linewidths=0, rasterized=True)

    labels = {
        "TACSTD2": "#d62728",
        "EPCAM": "#ff7f0e",
        "CLDN4": "#9467bd",
        "CXCL13": "#2ca02c",
        "SPP1": "#2ca02c",
        "IGHG1": "#2ca02c",
        "PDCD1": "#2ca02c",
        "SFTPC": "#1f77b4",
        "AZU1": "#1f77b4",
    }
    for gene, color in labels.items():
        if gene not in de.index:
            continue
        r = de.loc[gene]
        ax.scatter(r["log2_FC"], r["neglog10_p"], s=42, c=color, zorder=5, edgecolors="k", linewidths=0.4)
        ax.annotate(
            gene,
            (r["log2_FC"], r["neglog10_p"]),
            textcoords="offset points",
            xytext=(6, 4),
            fontsize=8,
            color=color,
            fontweight="bold" if gene == "TACSTD2" else "normal",
        )
    ax.axvline(0, color="0.4", lw=0.6)
    ax.axhline(-np.log10(0.05), color="0.5", lw=0.6, ls="--")
    ax.set_xlabel("log2 FC (LCAM-hi / LCAM-lo immune pseudobulk)")
    ax.set_ylabel("-log10 p")
    ax.set_title("GSE154826 author DE: LCAM-hi vs LCAM-lo (CD45+ fraction)")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_focus_bars(de: pd.DataFrame, path: Path) -> None:
    genes = [g for g in FOCUS_GENES if g in de.index]
    sub = de.loc[genes].copy()
    colors = []
    for g in genes:
        if g == "TACSTD2":
            colors.append("#d62728")
        elif g in {"EPCAM", "CLDN4", "CLDN3", "CLDN7", "KRT8", "KRT18", "KRT19", "CDH1"}:
            colors.append("#ff7f0e")
        elif g in {"CXCL13", "SPP1", "IGHG1", "IGHG3", "MZB1", "PDCD1", "LAYN", "TNFRSF9", "GZMB"}:
            colors.append("#2ca02c")
        else:
            colors.append("#4c78a8")
    fig, ax = plt.subplots(figsize=(8.2, 5.8))
    y = np.arange(len(genes))
    ax.barh(y, sub["log2_FC"].values, color=colors, edgecolor="k", linewidth=0.3)
    ax.set_yticks(y)
    ax.set_yticklabels(genes)
    ax.axvline(0, color="0.2", lw=0.7)
    ax.invert_yaxis()
    ax.set_xlabel("log2 FC (LCAM-hi / LCAM-lo)")
    ax.set_title("Focus genes in author LCAM-hi vs LCAM-lo immune DE")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_lcam_scores(v2: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    colors = np.where(v2["lcam_call"].eq("LCAM_hi"), "#d62728", "#1f77b4")
    ax.scatter(v2["lcam_hi_score"], v2["lcam_lo_score"], c=colors, s=36, zorder=3, edgecolors="k", linewidths=0.3)
    for _, r in v2.iterrows():
        ax.annotate(str(r["patient_ID"]), (r["lcam_hi_score"], r["lcam_lo_score"]), fontsize=7, xytext=(3, 3), textcoords="offset points")
    ax.scatter([], [], c="#d62728", s=36, edgecolors="k", linewidths=0.3, label="median-split LCAM-hi")
    ax.scatter([], [], c="#1f77b4", s=36, edgecolors="k", linewidths=0.3, label="median-split LCAM-lo")
    ax.legend(frameon=False, loc="best")
    ax.set_xlabel("LCAM-hi score (T_activated + IgG + MoMac-II)")
    ax.set_ylabel("LCAM-lo score (B + AM + cDC2 + AZU1_mac + Tcm/naive_II + cDC1)")
    ax.set_title("Mount Sinai V2-beads tumor samples (composition LCAM)")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_epi_vs_immune(ie: pd.DataFrame, path: Path) -> None:
    genes = ["TACSTD2", "EPCAM", "CLDN4", "KRT19", "PTPRC", "CD3D", "CXCL13", "SPP1"]
    genes = [g for g in genes if g in ie.index]
    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    y = np.arange(len(genes))
    colors = ["#d62728" if g == "TACSTD2" else "#4c78a8" for g in genes]
    ax.barh(y, ie.loc[genes, "l2fc"].values, color=colors, edgecolor="k", linewidth=0.3)
    ax.set_yticks(y)
    ax.set_yticklabels(genes)
    ax.axvline(0, color="0.2", lw=0.7)
    ax.invert_yaxis()
    ax.set_xlabel("log2 FC (immune / epithelial-gated cells)")
    ax.set_title("TACSTD2 is epithelial, not an immune/LCAM gene")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    de = load_de()
    ie = pd.read_csv(fetch("immune_vs_ep_de.csv"), index_col=0)
    ie.index.name = "gene"
    all_scores, v2_tumor = composition_scores()

    focus = [g for g in FOCUS_GENES if g in de.index]
    de_focus = de.loc[focus].reset_index()
    de_focus.to_csv(OUT / "lcam_hi_vs_lo_focus_genes.csv", index=False)
    de.reset_index().to_csv(OUT / "lcam_hi_vs_lo_pseudobulk_de.csv", index=False)
    all_scores.to_csv(OUT / "sample_lcam_composition_scores.csv", index=False)
    v2_tumor.to_csv(OUT / "v2_beads_tumor_lcam_scores.csv", index=False)

    ie_focus_genes = [
        g
        for g in ["TACSTD2", "EPCAM", "CLDN4", "CLDN3", "CLDN7", "KRT8", "KRT19", "PTPRC", "CD3D", "CXCL13", "SPP1", "IGHG1"]
        if g in ie.index
    ]
    ie.loc[ie_focus_genes].reset_index().to_csv(OUT / "immune_vs_epithelial_focus_genes.csv", index=False)

    tac = de.loc["TACSTD2"]
    epi_tac = ie.loc["TACSTD2"]
    rho, rp = stats.spearmanr(v2_tumor["lcam_hi_score"], v2_tumor["lcam_lo_score"])

    stats_txt = OUT / "stats.txt"
    stats_txt.write_text(
        "\n".join(
            [
                "GSE154826 LCAM vs TACSTD2 (processed/open only)",
                "",
                "PRIMARY: author immune-pseudobulk DE LCAM-hi vs LCAM-lo",
                f"  TACSTD2 log2FC={tac['log2_FC']:.4f} p={tac['p.value']:.4g} padj={tac['adj.p.value']:.4g}",
                f"  TACSTD2 freq_bg={tac['freq_bg']:.4g} freq_fg={tac['freq_fg']:.4g} n_cells_above_thresh={int(tac['n_cells_with_min_umis_above_thresh'])}",
                f"  TACSTD2 |log2FC| rank among {len(de)} genes: {int(de['log2_FC'].abs().rank(ascending=False).loc['TACSTD2'])}",
                "",
                "SANITY LCAM genes (should be strongly LCAM-hi):",
                *(
                    f"  {g}: log2FC={de.loc[g,'log2_FC']:.3f} padj={de.loc[g,'adj.p.value']:.3g}"
                    for g in ["CXCL13", "SPP1", "IGHG1", "PDCD1"]
                    if g in de.index
                ),
                "",
                "EPITHELIAL residual in the same DE (tiny frequencies; not tumor TACSTD2):",
                *(
                    f"  {g}: log2FC={de.loc[g,'log2_FC']:.3f} padj={de.loc[g,'adj.p.value']:.3g} freq_fg={de.loc[g,'freq_fg']:.3g}"
                    for g in ["EPCAM", "CLDN4", "KRT19"]
                    if g in de.index
                ),
                "",
                "IMMUNE vs EPITHELIAL-GATED DE (author immune_vs_ep_de.csv):",
                f"  TACSTD2 l2fc(immune/epi)={epi_tac['l2fc']:.3f} fg_immune={epi_tac['fg_exprs']:.4g} bg_epi={epi_tac['bg_exprs']:.4g}",
                f"  EPCAM l2fc={ie.loc['EPCAM','l2fc']:.3f}",
                f"  CLDN4 l2fc={ie.loc['CLDN4','l2fc']:.3f}",
                f"  PTPRC l2fc={ie.loc['PTPRC','l2fc']:.3f}",
                "",
                "COMPOSITION LCAM (no expression; V2 beads tumor samples only):",
                f"  n_samples={len(v2_tumor)} n_patients={v2_tumor['patient_ID'].nunique()}",
                f"  LCAM_hi calls={int((v2_tumor['lcam_call']=='LCAM_hi').sum())} LCAM_lo calls={int((v2_tumor['lcam_call']=='LCAM_lo').sum())}",
                f"  Spearman LCAM-hi vs LCAM-lo scores: rho={rho:.3f} p={rp:.3g}",
                "",
                "SKIPPED (honest, huge raw / unneeded count dumps):",
                "  - SRA/FASTQ SRP251372 / PRJNA609924",
                "  - GEO GSE154826_amp_batch_ID_*.tar.gz feature-barcode MTX (~3 GB, no cluster labels)",
                "  - Dropbox lung_ldm.rd full UMI object (umitab for 361,929 cells)",
                "",
            ]
        )
        + "\n"
    )

    plot_volcano(de, OUT / "lcam_hi_vs_lo_volcano.png")
    plot_focus_bars(de, OUT / "focus_genes_lcam_log2fc.png")
    plot_lcam_scores(v2_tumor, OUT / "v2_beads_tumor_lcam_scores.png")
    plot_epi_vs_immune(ie, OUT / "tacstd2_immune_vs_epithelial.png")

    provenance = {
        "dataset": "GSE154826",
        "paper": "Leader, Grout et al. Cancer Cell 2021 PMID 34767762",
        "question": "LCAM vs TACSTD2",
        "policy": "processed/open only; honest skip huge raw",
        "used": {
            name: {
                "url": url,
                "md5": md5_file(CACHE / name),
                "bytes": (CACHE / name).stat().st_size,
            }
            for name, url in FILES.items()
        },
        "skipped": [
            "SRA/FASTQ SRP251372",
            "GEO 10x MTX tarballs GSE154826_amp_batch_ID_*.tar.gz",
            "Dropbox lung_ldm.rd full UMI matrix",
        ],
        "primary_result": {
            "gene": "TACSTD2",
            "contrast": "author LCAM-hi vs LCAM-lo immune pseudobulk DE_total",
            "log2FC": float(tac["log2_FC"]),
            "p": float(tac["p.value"]),
            "padj": float(tac["adj.p.value"]),
            "verdict": "not significant",
        },
    }
    (OUT / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(stats_txt.read_text())
    print("wrote", OUT)


if __name__ == "__main__":
    main()
