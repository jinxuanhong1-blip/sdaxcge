#!/usr/bin/env python3
"""Honest overlap of the A8/A9 panel and keratin/TJ GSEA across datasets.

No genome-wide fishing. The question is: which of
  CLDN1, CLDN4, CLDN7, F11R, PARD3
are up in TACSTD2-high samples in each dataset, and which GSEA sets
replicate.

Primary bulk pair (the claim as stated): TCGA-LUAD + GSE207422.
TCGA-LUSC is a histology-matched squamous check.
GSE131907 is supporting scRNA-seq (tLung epithelial pseudobulk is the
comparable view; cell-level is exploratory and labelled as such).
"""
import json
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import lib_analysis as L

OUT = "results/claim_A8A9"
FDR = L.FDR_THRESHOLD

# dataset_id -> (label, panel_stats path, gsea path, role)
DATASETS = [
    ("tcga_luad", "TCGA-LUAD bulk",
     f"{OUT}/tcga_luad/panel_stats.csv", f"{OUT}/tcga_luad/gsea_prerank.csv",
     "primary"),
    ("gse207422", "GSE207422 NSCLC bulk",
     f"{OUT}/gse207422/panel_stats.csv", f"{OUT}/gse207422/gsea_prerank.csv",
     "primary"),
    ("tcga_lusc", "TCGA-LUSC bulk",
     f"{OUT}/tcga_lusc/panel_stats.csv", f"{OUT}/tcga_lusc/gsea_prerank.csv",
     "histology_check"),
    ("gse131907_tlung_pb", "GSE131907 tLung epi pseudobulk",
     f"{OUT}/gse131907/tlung_pseudobulk/panel_stats.csv",
     f"{OUT}/gse131907/tlung_pseudobulk/gsea_prerank.csv",
     "supporting_pseudobulk"),
    ("gse131907_alltumor_pb", "GSE131907 all-tumour-site epi pseudobulk",
     f"{OUT}/gse131907/alltumor_pseudobulk/panel_stats.csv",
     f"{OUT}/gse131907/alltumor_pseudobulk/gsea_prerank.csv",
     "sensitivity"),
    ("gse131907_tlung_cell", "GSE131907 tLung epi cell-level (exploratory)",
     f"{OUT}/gse131907/tlung_cell/panel_stats.csv",
     f"{OUT}/gse131907/tlung_cell/gsea_prerank.csv",
     "exploratory_cell"),
]


def _md_table(df):
    """Markdown table without requiring the tabulate extra."""
    df = df.copy()
    cols = list(df.columns)
    def fmt(v):
        if isinstance(v, (bool, np.bool_)):
            return "True" if v else "False"
        if isinstance(v, (float, np.floating)):
            if np.isnan(v):
                return ""
            return f"{v:.4g}"
        return str(v)
    lines = ["| " + " | ".join(cols) + " |",
             "| " + " | ".join("---" for _ in cols) + " |"]
    for _, row in df.iterrows():
        lines.append("| " + " | ".join(fmt(row[c]) for c in cols) + " |")
    return "\n".join(lines)


def load_existing():
    rows = []
    gsea_rows = []
    for did, label, ppath, gpath, role in DATASETS:
        if not os.path.exists(ppath):
            print(f"SKIP missing {ppath}")
            continue
        p = pd.read_csv(ppath, index_col=0)
        for gene in L.PANEL:
            if gene not in p.index:
                rows.append({
                    "dataset_id": did, "dataset": label, "role": role, "gene": gene,
                    "present": False, "up_in_trop2_high": False,
                })
                continue
            r = p.loc[gene]
            rows.append({
                "dataset_id": did, "dataset": label, "role": role, "gene": gene,
                "present": True,
                "n_total": int(r["n_total"]), "n_high": int(r["n_high"]), "n_low": int(r["n_low"]),
                "log2FC_high_vs_low": float(r["log2FC_high_vs_low"]),
                "spearman_rho": float(r["spearman_rho"]),
                "spearman_p": float(r["spearman_p"]),
                "spearman_fdr": float(r["spearman_fdr"]),
                "welch_p": float(r["welch_p"]),
                "welch_fdr": float(r["welch_fdr"]),
                "mwu_p": float(r["mwu_p"]),
                "up_in_trop2_high": bool(r["up_in_trop2_high"]),
                "same_direction": bool((r["spearman_rho"] > 0) and (r["log2FC_high_vs_low"] > 0)),
            })
        if os.path.exists(gpath):
            g = pd.read_csv(gpath)
            for _, gr in g.iterrows():
                gsea_rows.append({
                    "dataset_id": did, "dataset": label, "role": role,
                    "Term": gr["Term"],
                    "NES": float(gr["NES"]),
                    "NOM_p": float(gr["NOM p-val"]),
                    "FDR_q": float(gr["FDR q-val"]),
                    "ES": float(gr["ES"]),
                    "significant_fdr05": bool(gr["FDR q-val"] < FDR),
                    "positive_NES": bool(gr["NES"] > 0),
                })
    return pd.DataFrame(rows), pd.DataFrame(gsea_rows)


def intersection_table(panel_long):
    """Wide call matrix + explicit intersection sets."""
    call = panel_long.pivot_table(index="gene", columns="dataset_id",
                                  values="up_in_trop2_high", aggfunc="first")
    direction = panel_long.pivot_table(index="gene", columns="dataset_id",
                                       values="same_direction", aggfunc="first")
    recs = []
    for gene in L.PANEL:
        luad = bool(call.loc[gene, "tcga_luad"]) if "tcga_luad" in call.columns and gene in call.index else False
        gse = bool(call.loc[gene, "gse207422"]) if "gse207422" in call.columns and gene in call.index else False
        lusc = bool(call.loc[gene, "tcga_lusc"]) if "tcga_lusc" in call.columns and gene in call.index else False
        pb = bool(call.loc[gene, "gse131907_tlung_pb"]) if "gse131907_tlung_pb" in call.columns and gene in call.index else False
        cell = bool(call.loc[gene, "gse131907_tlung_cell"]) if "gse131907_tlung_cell" in call.columns and gene in call.index else False
        # direction-only for underpowered GSE131907 pseudobulk
        pb_dir = False
        if "gse131907_tlung_pb" in direction.columns and gene in direction.index:
            pb_dir = bool(direction.loc[gene, "gse131907_tlung_pb"])
        recs.append({
            "gene": gene,
            "pass_TCGA_LUAD": luad,
            "pass_GSE207422": gse,
            "pass_TCGA_LUSC": lusc,
            "pass_GSE131907_tLung_pseudobulk": pb,
            "same_direction_GSE131907_tLung_pseudobulk": pb_dir,
            "pass_GSE131907_tLung_cell_exploratory": cell,
            "intersection_LUAD_and_GSE207422": luad and gse,
            "intersection_LUAD_GSE207422_and_GSE131907_pb": luad and gse and pb,
            "intersection_LUAD_GSE207422_and_GSE131907_pb_direction": luad and gse and pb_dir,
            "intersection_LUAD_LUSC_GSE207422": luad and lusc and gse,
        })
    return pd.DataFrame(recs).set_index("gene"), call


def plot_panel(panel_long, path):
    ds_order = [d[0] for d in DATASETS if d[0] in set(panel_long["dataset_id"])]
    labels = {d[0]: d[1] for d in DATASETS}
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, col, title, cmap, v in zip(
        axes,
        ["spearman_rho", "log2FC_high_vs_low"],
        ["Spearman rho vs TACSTD2", "log2FC (TACSTD2 high vs low tertile)"],
        ["RdBu_r", "RdBu_r"],
        [0.6, 1.5],
    ):
        mat = panel_long.pivot_table(index="gene", columns="dataset_id", values=col, aggfunc="first")
        mat = mat.reindex(index=L.PANEL, columns=ds_order)
        im = ax.imshow(mat.to_numpy(dtype=float), cmap=cmap, vmin=-v, vmax=v, aspect="auto")
        ax.set_xticks(range(len(ds_order)))
        ax.set_xticklabels([labels[d] for d in ds_order], rotation=35, ha="right", fontsize=8)
        ax.set_yticks(range(len(L.PANEL)))
        ax.set_yticklabels(L.PANEL)
        ax.set_title(title, fontsize=10)
        # mark failures of the strict call
        call = panel_long.pivot_table(index="gene", columns="dataset_id",
                                      values="up_in_trop2_high", aggfunc="first")
        call = call.reindex(index=L.PANEL, columns=ds_order)
        for i, g in enumerate(L.PANEL):
            for j, d in enumerate(ds_order):
                val = mat.iloc[i, j]
                passed = bool(call.iloc[i, j]) if pd.notna(call.iloc[i, j]) else False
                txt = "·" if pd.isna(val) else f"{val:.2f}"
                ax.text(j, i, txt, ha="center", va="center",
                        fontsize=7, color="white" if passed else "black",
                        fontweight="bold" if passed else "normal")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.suptitle("A8/A9 panel vs TACSTD2  |  bold white = up_in_trop2_high (rho>0, log2FC>0, Welch FDR<0.05)",
                 fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_gsea(gsea_long, path):
    if gsea_long.empty:
        return
    ds_order = [d[0] for d in DATASETS if d[0] in set(gsea_long["dataset_id"])]
    labels = {d[0]: d[1] for d in DATASETS}
    terms = list(gsea_long["Term"].drop_duplicates())
    mat = gsea_long.pivot_table(index="Term", columns="dataset_id", values="NES", aggfunc="first")
    mat = mat.reindex(index=terms, columns=ds_order)
    fdr = gsea_long.pivot_table(index="Term", columns="dataset_id", values="FDR_q", aggfunc="first")
    fdr = fdr.reindex(index=terms, columns=ds_order)
    fig, ax = plt.subplots(figsize=(11, 5.2))
    im = ax.imshow(mat.to_numpy(dtype=float), cmap="RdBu_r", vmin=-2.5, vmax=2.5, aspect="auto")
    ax.set_xticks(range(len(ds_order)))
    ax.set_xticklabels([labels[d] for d in ds_order], rotation=35, ha="right", fontsize=8)
    ax.set_yticks(range(len(terms)))
    ax.set_yticklabels(terms, fontsize=8)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            nes = mat.iloc[i, j]
            q = fdr.iloc[i, j]
            if pd.isna(nes):
                continue
            star = "*" if pd.notna(q) and q < FDR else ""
            ax.text(j, i, f"{nes:.2f}{star}", ha="center", va="center", fontsize=7)
    ax.set_title("Pre-ranked GSEA NES (ranking = Spearman vs TACSTD2). * FDR q < 0.05")
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="NES")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def write_markdown(panel_long, gsea_long, inter, path):
    lines = []
    lines.append("# Claim A8/A9 — TROP2-high keratin / tight-junction GSEA")
    lines.append("")
    lines.append("Computed numbers only. No fabricated statistics.")
    lines.append("")
    lines.append("## Question")
    lines.append("")
    lines.append("Are keratin and tight-junction (TJ) programs enriched among genes")
    lines.append("positively associated with **TACSTD2 (TROP2)** in lung cancer, and do")
    lines.append("the TJ genes **CLDN1, CLDN4, CLDN7, F11R, PARD3** themselves track TROP2-high?")
    lines.append("")
    lines.append("## Method (identical in every dataset)")
    lines.append("")
    lines.append("- Ranking metric for pre-ranked GSEA: Spearman correlation of each gene with TACSTD2.")
    lines.append("- TROP2-high vs TROP2-low: top vs bottom tertile of TACSTD2.")
    lines.append("- A panel gene is called `up_in_trop2_high` only if **all** of: rho > 0,")
    lines.append("  log2FC(high vs low) > 0, Welch t-test BH-FDR < 0.05.")
    lines.append("- Gene sets: Enrichr GO BP/CC 2021, KEGG 2021 Human, Reactome 2022")
    lines.append("  (local GMT in `data/genesets/`; see provenance JSON).")
    lines.append("- GSEA: gseapy prerank, 1000 permutations, min_size=3.")
    lines.append("")
    lines.append("## Datasets")
    lines.append("")
    lines.append("| id | description | role |")
    lines.append("|---|---|---|")
    lines.append("| tcga_luad | TCGA-LUAD STAR log2(TPM+1), primary tumours, 1/patient | primary bulk |")
    lines.append("| gse207422 | GSE207422 NSCLC pre-treatment bulk log2TPM (n=24) | primary bulk |")
    lines.append("| tcga_lusc | TCGA-LUSC, same processing | squamous histology check |")
    lines.append("| gse131907_tlung_pb | GSE131907 tLung epithelial sample pseudobulk | supporting, small n |")
    lines.append("| gse131907_alltumor_pb | + tL/B, mLN, mBrain, PE epithelial pseudobulk | sensitivity |")
    lines.append("| gse131907_tlung_cell | tLung epithelial cells, cell-level | exploratory (cells not independent) |")
    lines.append("")
    lines.append("## Panel gene calls")
    lines.append("")
    show = panel_long.copy()
    cols = ["dataset", "role", "gene", "n_total", "log2FC_high_vs_low",
            "spearman_rho", "welch_fdr", "up_in_trop2_high", "same_direction"]
    lines.append(_md_table(show[cols]))
    lines.append("")
    lines.append("## Honest intersection of the five genes")
    lines.append("")
    lines.append("Intersection is **not** a genome-wide overlap hunt. It is the subset of")
    lines.append("{CLDN1, CLDN4, CLDN7, F11R, PARD3} that pass the same call in the named datasets.")
    lines.append("")
    lines.append(_md_table(inter.reset_index()))
    lines.append("")
    primary_pass = inter.index[inter["intersection_LUAD_and_GSE207422"]].tolist()
    plus_pb = inter.index[inter["intersection_LUAD_GSE207422_and_GSE131907_pb"]].tolist()
    plus_dir = inter.index[inter["intersection_LUAD_GSE207422_and_GSE131907_pb_direction"]].tolist()
    luad_lusc_gse = inter.index[inter["intersection_LUAD_LUSC_GSE207422"]].tolist()
    lines.append(f"- **LUAD ∩ GSE207422 (primary):** {primary_pass or 'none'}")
    lines.append(f"- **LUAD ∩ LUSC ∩ GSE207422:** {luad_lusc_gse or 'none'}")
    lines.append(f"- **LUAD ∩ GSE207422 ∩ GSE131907 tLung pseudobulk (strict):** {plus_pb or 'none'}")
    lines.append(f"- **LUAD ∩ GSE207422 ∩ GSE131907 tLung pseudobulk (direction only):** {plus_dir or 'none'}")
    lines.append("")
    lines.append("PARD3 is listed because it is in the a priori panel, not because it passed.")
    lines.append("A False in the table is a real negative under this criterion, not a missing value.")
    lines.append("")
    lines.append("## GSEA")
    lines.append("")
    if gsea_long.empty:
        lines.append("(no GSEA tables found)")
    else:
        gshow = gsea_long[["dataset", "role", "Term", "NES", "NOM_p", "FDR_q", "significant_fdr05"]]
        lines.append(_md_table(gshow))
    lines.append("")
    lines.append("## What this does **not** show")
    lines.append("")
    lines.append("- Causality or that TROP2 drives TJ/keratin transcription.")
    lines.append("- Protein-level TROP2 or claudin status.")
    lines.append("- A genome-wide intersecting signature beyond the five a priori genes.")
    lines.append("- Independent samples in the GSE131907 cell-level view.")
    lines.append("- Adequately powered sample-level tests in GSE131907 tLung (n is small).")
    lines.append("")
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def main():
    panel_long, gsea_long = load_existing()
    if panel_long.empty:
        raise SystemExit("no panel_stats found — run the dataset scripts first")
    panel_long.to_csv(f"{OUT}/panel_long.csv", index=False)
    if not gsea_long.empty:
        gsea_long.to_csv(f"{OUT}/gsea_long.csv", index=False)

    inter, call = intersection_table(panel_long)
    inter.to_csv(f"{OUT}/panel_intersection.csv")
    call.to_csv(f"{OUT}/panel_call_matrix.csv")

    plot_panel(panel_long, f"{OUT}/fig_panel_heatmap.png")
    plot_gsea(gsea_long, f"{OUT}/fig_gsea_nes.png")
    write_markdown(panel_long, gsea_long, inter, f"{OUT}/REPORT.md")

    summary = {
        "primary_intersection_LUAD_GSE207422": inter.index[inter["intersection_LUAD_and_GSE207422"]].tolist(),
        "intersection_LUAD_LUSC_GSE207422": inter.index[inter["intersection_LUAD_LUSC_GSE207422"]].tolist(),
        "intersection_plus_GSE131907_pb_strict": inter.index[inter["intersection_LUAD_GSE207422_and_GSE131907_pb"]].tolist(),
        "intersection_plus_GSE131907_pb_direction": inter.index[inter["intersection_LUAD_GSE207422_and_GSE131907_pb_direction"]].tolist(),
        "criterion": "up_in_trop2_high = spearman_rho>0 AND log2FC_high_vs_low>0 AND welch_fdr<0.05",
        "panel": L.PANEL,
        "note": "Empty lists are real negatives, not missing analyses.",
    }
    with open(f"{OUT}/intersection_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(summary, indent=2))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
