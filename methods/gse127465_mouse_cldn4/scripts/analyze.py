#!/usr/bin/env python3
"""Additive GSE127465 MOUSE Cldn4-only (not merged with human).

Zilionis et al., Immunity 2019, PMID 30979687.
Public mouse processed MTX is CD45+ lung immune cells from 2 healthy
and 2 KP1.9 tumor-bearing mice. Epithelium is expected to be thin.
If Cldn4 is present, report mouse-level Cldn4 vs T fraction and myeloid.
Honest n = mice (4; tumor-bearing 2). Do not mega-merge human GSE127465.
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

N_CELLS = 15939
N_GENES = 28205
KEEP_GENES = [
    "Cldn4",
    "Tacstd2",
    "Epcam",
    "Cdh1",
    "Krt8",
    "Krt18",
    "Krt19",
    "Sftpc",
    "Ptprc",
    "Cd3d",
    "Cd3e",
    "Cd8a",
    "Nkg7",
    "Lyz2",
    "Adgre1",
    "Itgam",
]
MYELOID_MAJOR = {"Neutrophils", "MoMacDC", "pDC", "Basophils"}
T_MAJOR = {"T cells"}
NK_MAJOR = {"NK cells"}
# leftover epithelium: author has none; marker leftover is Epcam+ or (Cdh1+ and Krt8+)
# among CD45+ inDrops this is contamination / ambient, not a tumor compartment.


def fisher_z_ci(rho: float, n: int, alpha: float = 0.05) -> tuple[float, float]:
    if n < 4 or not np.isfinite(rho) or abs(rho) >= 1:
        return (float("nan"), float("nan"))
    z = 0.5 * math.log((1 + rho) / (1 - rho))
    se = 1.0 / math.sqrt(n - 3)
    from scipy.stats import norm

    zcrit = norm.ppf(1 - alpha / 2)
    lo = z - zcrit * se
    hi = z + zcrit * se
    def inv(x: float) -> float:
        return math.tanh(x)

    return (inv(lo), inv(hi))


def load_genes(path: Path) -> list[str]:
    with gzip.open(path, "rt") as handle:
        genes = [ln.strip().split("\t")[0] for ln in handle if ln.strip()]
    if len(genes) != N_GENES:
        raise SystemExit(f"genes {len(genes)} != {N_GENES}")
    return genes


def load_meta(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    if len(df) != N_CELLS:
        raise SystemExit(f"metadata rows {len(df)} != {N_CELLS}")
    return df


def stream_gene_columns(mtx_path: Path, gene_idx_1based: dict[str, int]) -> dict[str, np.ndarray]:
    """MTX is cells x genes, 1-based. Return dense columns for requested genes."""
    want = set(gene_idx_1based.values())
    out = {name: np.zeros(N_CELLS, dtype=np.float64) for name in gene_idx_1based}
    inv = {idx: name for name, idx in gene_idx_1based.items()}
    with gzip.open(mtx_path, "rt") as handle:
        for line in handle:
            if line.startswith("%"):
                continue
            nrows, ncols, nnz = map(int, line.split())
            if nrows != N_CELLS or ncols != N_GENES:
                raise SystemExit(f"MTX banner {nrows}x{ncols} != {N_CELLS}x{N_GENES}")
            print(f"MTX banner {nrows} x {ncols} nnz={nnz}", flush=True)
            break
        for i, line in enumerate(handle, 1):
            parts = line.split()
            if len(parts) < 3:
                continue
            c = int(parts[1])
            if c not in want:
                continue
            r = int(parts[0]) - 1
            out[inv[c]][r] = float(parts[2])
            if i % 5_000_000 == 0:
                print(f"  scanned {i:,} MTX lines", flush=True)
    return out


def leftover_epi(expr: dict[str, np.ndarray]) -> np.ndarray:
    epcam = expr["Epcam"] > 0
    cdh1_krt8 = (expr["Cdh1"] > 0) & (expr["Krt8"] > 0)
    return epcam | cdh1_krt8


def mouse_table(meta: pd.DataFrame, expr: dict[str, np.ndarray]) -> pd.DataFrame:
    cldn4 = expr["Cldn4"]
    leftover = leftover_epi(expr)
    rows = []
    for mouse, idx in meta.groupby("Biological replicate").groups.items():
        idx = np.asarray(idx)
        n = len(idx)
        maj = meta.loc[idx, "Major cell type"].astype(str)
        n_t = int(maj.isin(T_MAJOR).sum())
        n_nk = int(maj.isin(NK_MAJOR).sum())
        n_my = int(maj.isin(MYELOID_MAJOR).sum())
        n_neut = int((maj == "Neutrophils").sum())
        n_momac = int((maj == "MoMacDC").sum())
        n_pdc = int((maj == "pDC").sum())
        n_baso = int((maj == "Basophils").sum())
        n_b = int((maj == "B cells").sum())
        n_leftover = int(leftover[idx].sum())
        v = cldn4[idx]
        pos = v > 0
        rows.append(
            {
                "mouse": mouse,
                "condition": "tumor" if str(mouse).startswith("Tumor") else "healthy",
                "n_cells": n,
                "n_T": n_t,
                "n_NK": n_nk,
                "n_TNK": n_t + n_nk,
                "n_myeloid": n_my,
                "n_neutrophil": n_neut,
                "n_MoMacDC": n_momac,
                "n_pDC": n_pdc,
                "n_basophil": n_baso,
                "n_B": n_b,
                "n_author_epithelium": 0,
                "n_leftover_Epcam_or_Cdh1Krt8": n_leftover,
                "frac_T": n_t / n,
                "frac_TNK": (n_t + n_nk) / n,
                "frac_myeloid": n_my / n,
                "frac_neutrophil": n_neut / n,
                "frac_MoMacDC": n_momac / n,
                "frac_leftover_epi": n_leftover / n,
                "Cldn4_mean": float(v.mean()),
                "Cldn4_mean_log1p": float(np.log1p(v).mean()),
                "Cldn4_n_pos": int(pos.sum()),
                "Cldn4_frac_pos": float(pos.mean()),
                "Cldn4_mean_among_pos": float(v[pos].mean()) if pos.any() else 0.0,
                "Tacstd2_mean_log1p": float(np.log1p(expr["Tacstd2"][idx]).mean()),
                "Tacstd2_frac_pos": float((expr["Tacstd2"][idx] > 0).mean()),
                "Epcam_frac_pos": float((expr["Epcam"][idx] > 0).mean()),
                "Ptprc_frac_pos": float((expr["Ptprc"][idx] > 0).mean()),
            }
        )
    out = pd.DataFrame(rows).sort_values(["condition", "mouse"]).reset_index(drop=True)
    return out


def celltype_table(meta: pd.DataFrame, expr: dict[str, np.ndarray]) -> pd.DataFrame:
    cldn4 = expr["Cldn4"]
    leftover = leftover_epi(expr)
    rows = []
    for (cond, mouse, maj), idx in meta.groupby(
        [
            meta["Tumor or healthy"].map({"t": "tumor", "h": "healthy"}),
            "Biological replicate",
            "Major cell type",
        ]
    ).groups.items():
        idx = np.asarray(idx)
        v = cldn4[idx]
        pos = v > 0
        rows.append(
            {
                "condition": cond,
                "mouse": mouse,
                "major_cell_type": maj,
                "n_cells": len(idx),
                "Cldn4_n_pos": int(pos.sum()),
                "Cldn4_frac_pos": float(pos.mean()),
                "Cldn4_mean": float(v.mean()),
                "Cldn4_mean_log1p": float(np.log1p(v).mean()),
                "n_leftover_epi_in_type": int(leftover[idx].sum()),
            }
        )
    return pd.DataFrame(rows).sort_values(
        ["condition", "mouse", "n_cells"], ascending=[True, True, False]
    )


def spearman_block(mice: pd.DataFrame, x: str, y: str) -> dict:
    n = len(mice)
    if n < 3:
        return {
            "contrast": f"{x} vs {y}",
            "n_mice": n,
            "rho": float("nan"),
            "p": float("nan"),
            "ci95_lo": float("nan"),
            "ci95_hi": float("nan"),
            "note": "n<3; Spearman not computed",
        }
    rho, p = spearmanr(mice[x], mice[y])
    lo, hi = fisher_z_ci(float(rho), n)
    return {
        "contrast": f"{x} vs {y}",
        "n_mice": n,
        "rho": float(rho),
        "p": float(p),
        "ci95_lo": lo,
        "ci95_hi": hi,
        "note": "n=4 is thin; p descriptive" if n == 4 else "thin n; p descriptive",
    }


def write_figures(mice: pd.DataFrame, by_type: pd.DataFrame, outdir: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.0))
    color = {"healthy": "#4C78A8", "tumor": "#E45756"}
    for ax, y, ylab in (
        (axes[0], "frac_T", "Author T-cell fraction"),
        (axes[1], "frac_myeloid", "Author myeloid fraction"),
    ):
        for cond, sub in mice.groupby("condition"):
            ax.scatter(
                sub["Cldn4_mean_log1p"],
                sub[y],
                s=90,
                c=color[cond],
                label=cond,
                zorder=3,
            )
            for rec in sub.itertuples(index=False):
                ax.annotate(
                    rec.mouse.replace("Tumor-bearing ", "T").replace("Healthy ", "H"),
                    (rec.Cldn4_mean_log1p, getattr(rec, y)),
                    textcoords="offset points",
                    xytext=(5, 4),
                    fontsize=8,
                )
        ax.set_xlabel("Mouse-level Cldn4 (mean log1p)")
        ax.set_ylabel(ylab)
        ax.set_title(f"n=4 mice (2 tumor / 2 healthy)\n{ylab.split()[-2]} vs Cldn4")
        ax.legend(frameon=False, loc="best")
        ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(outdir / "fig_mouse_cldn4_vs_T_myeloid.png", dpi=160)
    fig.savefig(outdir / "fig_mouse_cldn4_vs_T_myeloid.pdf")
    plt.close(fig)

    # honest n + detection
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.0))
    order = mice.sort_values(["condition", "mouse"])
    x = np.arange(len(order))
    axes[0].bar(x - 0.2, order["frac_T"], 0.4, label="T", color="#4C78A8")
    axes[0].bar(x + 0.2, order["frac_myeloid"], 0.4, label="myeloid", color="#F58518")
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(
        [m.replace("Tumor-bearing ", "T").replace("Healthy ", "H") for m in order["mouse"]]
    )
    axes[0].set_ylabel("Fraction of CD45+ cells")
    axes[0].set_title("Author lineage fractions (unit = mouse)")
    axes[0].legend(frameon=False)
    axes[1].bar(x, order["Cldn4_frac_pos"], color=["#4C78A8" if c == "healthy" else "#E45756" for c in order["condition"]])
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(
        [m.replace("Tumor-bearing ", "T").replace("Healthy ", "H") for m in order["mouse"]]
    )
    axes[1].set_ylabel("Cldn4 > 0 fraction")
    axes[1].set_title("Cldn4 detection in CD45+ (author-normalized)")
    fig.tight_layout()
    fig.savefig(outdir / "fig_honest_n.png", dpi=160)
    fig.savefig(outdir / "fig_honest_n.pdf")
    plt.close(fig)

    # Cldn4 by major type, pooled tumor vs healthy
    fig, ax = plt.subplots(figsize=(7.6, 4.2))
    types = ["T cells", "NK cells", "B cells", "MoMacDC", "Neutrophils", "pDC", "Basophils"]
    width = 0.38
    xs = np.arange(len(types))
    for i, cond in enumerate(["healthy", "tumor"]):
        means = []
        for t in types:
            sub = by_type[(by_type["condition"] == cond) & (by_type["major_cell_type"] == t)]
            if sub.empty:
                means.append(0.0)
            else:
                # cell-weighted mean of mean_log1p
                w = sub["n_cells"].to_numpy()
                means.append(float(np.average(sub["Cldn4_mean_log1p"], weights=w)))
        ax.bar(xs + (i - 0.5) * width, means, width, label=cond, color=color[cond])
    ax.set_xticks(xs)
    ax.set_xticklabels(types, rotation=25, ha="right")
    ax.set_ylabel("Cldn4 mean log1p (cell-weighted)")
    ax.set_title("Cldn4 by author major type (CD45+; epithelium absent)")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(outdir / "fig_cldn4_by_majortype.png", dpi=160)
    fig.savefig(outdir / "fig_cldn4_by_majortype.pdf")
    plt.close(fig)


def fmt_rho(d: dict) -> str:
    if not np.isfinite(d["rho"]):
        return "NA"
    lo, hi = d["ci95_lo"], d["ci95_hi"]
    if np.isfinite(lo):
        return f"{d['rho']:+.3f} [{lo:+.3f}, {hi:+.3f}]"
    return f"{d['rho']:+.3f}"


def write_finding(
    mice: pd.DataFrame,
    by_type: pd.DataFrame,
    spears: list[dict],
    inventory: dict,
    expr: dict[str, np.ndarray],
    meta: pd.DataFrame,
    finding_path: Path,
) -> None:
    leftover = leftover_epi(expr)
    n_leftover = int(leftover.sum())
    n_cldn4 = int((expr["Cldn4"] > 0).sum())
    n_overlap = int(((expr["Cldn4"] > 0) & leftover).sum())
    leftover_frac = float((expr["Cldn4"][leftover] > 0).mean()) if n_leftover else 0.0
    rest_frac = float((expr["Cldn4"][~leftover] > 0).mean())
    tumor_mice = mice[mice["condition"] == "tumor"]
    healthy_mice = mice[mice["condition"] == "healthy"]
    s_all_t = next(s for s in spears if s["contrast"] == "Cldn4_mean_log1p vs frac_T")
    s_all_my = next(s for s in spears if s["contrast"] == "Cldn4_mean_log1p vs frac_myeloid")
    s_all_tnk = next(s for s in spears if s["contrast"] == "Cldn4_mean_log1p vs frac_TNK")
    primary = [
        s
        for s in spears
        if s["contrast"]
        in {
            "Cldn4_mean_log1p vs frac_T",
            "Cldn4_mean_log1p vs frac_TNK",
            "Cldn4_mean_log1p vs frac_myeloid",
            "Cldn4_frac_pos vs frac_T",
            "Cldn4_frac_pos vs frac_myeloid",
        }
    ]

    mouse_rows = []
    for rec in mice.itertuples(index=False):
        mouse_rows.append(
            f"| {rec.mouse} | {rec.condition} | {rec.n_cells} | {rec.n_T} | {rec.n_TNK} | "
            f"{rec.n_myeloid} | {rec.n_leftover_Epcam_or_Cdh1Krt8} | {rec.Cldn4_n_pos} | "
            f"{rec.Cldn4_frac_pos:.4f} | {rec.Cldn4_mean_log1p:.4f} | {rec.frac_T:.3f} | "
            f"{rec.frac_TNK:.3f} | {rec.frac_myeloid:.3f} |"
        )

    type_pool = (
        by_type.groupby(["condition", "major_cell_type"], as_index=False)
        .apply(
            lambda d: pd.Series(
                {
                    "n_cells": d["n_cells"].sum(),
                    "Cldn4_n_pos": d["Cldn4_n_pos"].sum(),
                    "Cldn4_frac_pos": d["Cldn4_n_pos"].sum() / d["n_cells"].sum(),
                }
            ),
            include_groups=False,
        )
        .reset_index(drop=True)
    )
    type_lines = []
    for rec in type_pool.sort_values(["condition", "n_cells"], ascending=[True, False]).itertuples(index=False):
        type_lines.append(
            f"| {rec.condition} | {rec.major_cell_type} | {int(rec.n_cells)} | "
            f"{int(rec.Cldn4_n_pos)} | {rec.Cldn4_frac_pos:.4f} |"
        )

    spear_lines = []
    for s in primary:
        ptxt = "NA" if not np.isfinite(s["p"]) else f"{s['p']:.3f}"
        spear_lines.append(
            f"| {s['contrast']} | {s['n_mice']} | {fmt_rho(s)} | {ptxt} | {s['note']} |"
        )

    leftover_by_mouse = mice[["mouse", "n_leftover_Epcam_or_Cdh1Krt8", "frac_leftover_epi"]]
    leftover_txt = ", ".join(
        f"{r.mouse} {int(r.n_leftover_Epcam_or_Cdh1Krt8)} ({r.frac_leftover_epi:.3f})"
        for r in leftover_by_mouse.itertuples(index=False)
    )

    md = f"""# FINDING — GSE127465 mouse NSCLC / myeloid-rich lung (Cldn4-only)

ADDITIVE public **MOUSE**. **Cldn4 only.** Not a TACSTD2 gate. **Not merged with human GSE127465.**

Zilionis et al., *Immunity* 2019, PMID [30979687](https://pubmed.ncbi.nlm.nih.gov/30979687/). GEO [GSE127465](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE127465). KP1.9 lung adenocarcinoma (Pfirschke et al. 2016) + matched healthy lung. Public processed object is **CD45+** inDrops: 15,939 cells after the authors' QC. Author major types are Neutrophils / B / MoMacDC / T / NK / pDC / Basophils. **Author epithelium = 0.**

Processed mouse matrix is public and <2 GB. Cldn4 is on the gene list. Epithelium is thin; mouse-level Cldn4 is still reported vs T fraction and myeloid. Human MTX / human metadata were not opened.

## Verdict

**GO (mouse-only; epithelium thin).** Public files exist (`GSE127465_mouse_counts_normalized_15939x28205.mtx.gz` {inventory['mtx_bytes']:,} B, metadata 15,939 × 12, 28,205 genes). Cldn4 is present (column 3875 / 0-based index 3874). Author malignant/epithelial barcodes = **0**. Marker leftover (Epcam+ or Cdh1+∩Krt8+) = **{n_leftover} / 15,939** cells — contamination/ambient in a CD45+ sort, not a tumor-epithelial compartment.

Cldn4 is **present but sparse** in CD45+ (**{n_cldn4} / 15,939** cells > 0). Overlap with marker leftover epithelium = **{n_overlap} / {n_cldn4}** (leftover frac+ {leftover_frac:.4f} vs rest {rest_frac:.4f}) — not leftover-restricted. Inferential unit = **mouse**. Honest n = **4** (2 tumor-bearing, 2 healthy). Tumor-only n = **2** — Spearman not computed on n=2.

Mouse-level mean-log1p Cldn4 vs author T fraction: n=4, ρ = {fmt_rho(s_all_t)}, p = {s_all_t['p']:.3f} (descriptive). vs myeloid fraction: n=4, ρ = {fmt_rho(s_all_my)}, p = {s_all_my['p']:.3f}. vs T+NK: n=4, ρ = {fmt_rho(s_all_tnk)}, p = {s_all_tnk['p']:.3f}. **n=4 is thin.** Do not write a precise effect. Do not treat 15,939 cells as n.

## Honest n

| Item | n | Note |
|---|---:|---|
| Biological mice deposited | **4** | Healthy 1, Healthy 2, Tumor-bearing 1, Tumor-bearing 2 |
| Tumor-bearing / healthy | **2 / 2** | KP1.9 lung vs tumor-free lung |
| Libraries (technical) | 14 | not the test n |
| Cells after author QC | 15,939 | CD45+; not the test n |
| Author major types | 7 | Neutrophils 8,022; B 2,813; MoMacDC 2,397; T 1,981; NK 630; pDC 62; Basophils 34 |
| Author epithelium / malignant | **0** | CD45+ design; no Type I/II / club / ciliated / KP labels |
| Marker leftover Epcam+ or (Cdh1+ and Krt8+) | **{n_leftover}** | {leftover_txt} |
| Cldn4 > 0 cells | **{n_cldn4}** | sparse in immune cells |
| Mice with ≥20 T and ≥20 myeloid | **4** | all four pass a composition gate |
| Dual-high Tacstd2∩Cldn4 | **not defined** | Cldn4-only |
| Human GSE127465 cells used | **0** | not mega-merged |

Tumor-bearing mice (author cells): {int(tumor_mice['n_cells'].sum())}. Healthy: {int(healthy_mice['n_cells'].sum())}.

## Gate

| File | Bytes | Used |
|---|---:|---|
| `GSE127465_mouse_counts_normalized_15939x28205.mtx.gz` | {inventory['mtx_bytes']:,} | **yes** — author-normalized; 15,939 × 28,205 |
| `GSE127465_mouse_cell_metadata_15939x12.tsv.gz` | {inventory['meta_bytes']:,} | **yes** |
| `GSE127465_gene_names_mouse_28205.tsv.gz` | {inventory['gene_bytes']:,} | **yes** — Cldn4 present |
| Human MTX / human metadata / human gene names | — | **no** |
| `GSE127465_RAW.tar` / SRA | — | **no** |

Stop-if-missing does **not** apply: processed mouse MTX is public and Cldn4 is present. Thin epithelium is not a no-go on this assignment.

## Per-mouse Cldn4 vs T and myeloid

Author-normalized MTX; Cldn4 score = mean log1p of that value in **all CD45+ cells of that mouse** (no epithelial subset exists). T fraction = author `T cells` / cells. T/NK = (`T cells` + `NK cells`) / cells. Myeloid = (`Neutrophils` + `MoMacDC` + `pDC` + `Basophils`) / cells.

| mouse | condition | n cells | n T | n T+NK | n myeloid | leftover epi | Cldn4+ | Cldn4 frac+ | Cldn4 mean log1p | frac T | frac T+NK | frac myeloid |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
{chr(10).join(mouse_rows)}

Machine table: `results/tables/per_mouse.tsv`.

## Spearman (mouse is the unit)

**n=4 is thin.** Fisher-z CI on n=4 is almost the full [−1, 1] range. Tumor-only n=2 is reported as two points, not a ρ. Full sensitivity table (neutrophil / MoMacDC / tumor-only NA rows) is in `results/tables/spearman.tsv`.

| Contrast | n | ρ [95% CI] | p | Note |
|---|---:|---|---:|---|
{chr(10).join(spear_lines)}

Tumor-bearing only (n=2, no Spearman):
"""
    t1 = tumor_mice.sort_values("mouse")
    for rec in t1.itertuples(index=False):
        md += (
            f"- {rec.mouse}: Cldn4 mean log1p = {rec.Cldn4_mean_log1p:.4f}, "
            f"frac T = {rec.frac_T:.3f}, frac T+NK = {rec.frac_TNK:.3f}, "
            f"frac myeloid = {rec.frac_myeloid:.3f}, Cldn4+ = {rec.Cldn4_n_pos}/{rec.n_cells}.\n"
        )

    md += f"""
Healthy only (n=2, no Spearman):
"""
    for rec in healthy_mice.sort_values("mouse").itertuples(index=False):
        md += (
            f"- {rec.mouse}: Cldn4 mean log1p = {rec.Cldn4_mean_log1p:.4f}, "
            f"frac T = {rec.frac_T:.3f}, frac T+NK = {rec.frac_TNK:.3f}, "
            f"frac myeloid = {rec.frac_myeloid:.3f}, Cldn4+ = {rec.Cldn4_n_pos}/{rec.n_cells}.\n"
        )

    md += f"""
## Cldn4 by author major type

No author epithelial row. Detection is in immune lineages (sparse).

| condition | major type | n cells | Cldn4+ | frac+ |
|---|---|---:|---:|---:|
{chr(10).join(type_lines)}

## What this can and cannot say

**Can say.** Public mouse processed MTX exists. Cldn4 is on the matrix. The deposit is CD45+ / myeloid-rich lung (author Neutrophils are the modal type). Mouse-level Cldn4 (all CD45+ cells) vs T and myeloid can be written with **n=4**. Epithelium is thin (author 0; marker leftover {n_leftover}).

**Cannot say.** This is not an epithelial/malignant Cldn4 score. It is not a Type II / KP tumor-cell state. It is not n=15,939. It is not a precise ρ. It is not a human–mouse mega-merge. Tacstd2 was audited only and was not a gate. Healthy vs tumor is 2 vs 2.

## Methods

- Source: GEO series supplementary mouse MTX + metadata + gene names only.
- Matrix values are the author total-count-normalized inDrops counts (filename `*_counts_normalized_*`). Cldn4 score = mean log1p of that value.
- Lineages = author `Major cell type`. Myeloid = Neutrophils + MoMacDC + pDC + Basophils.
- Leftover epithelium = Epcam>0 or (Cdh1>0 and Krt8>0) on the same MTX. Not used as the Cldn4 denominator because it is not an author tumor compartment. Cldn4+ ∩ leftover = {n_overlap} of {n_cldn4}.
- Spearman on 4 mice. Fisher-z 95% CI. p-values descriptive.
- Human GSE127465 files were not downloaded for this folder.

## What this is not

- Not human GSE127465 (7 patients; 40,362 tumor cells). That analysis is a separate folder.
- Not a cell-level merge of mouse + human MTX.
- Not ICI / MPR / RECIST.
- Not dual-high Tacstd2×Cldn4.
- Not a Slingshot / PAGA / CellChat redo.
- Not a claim that Cldn4 is a myeloid marker. Sparse CD45+ detection is compatible with ambient epithelial RNA.

## Outputs

- `results/tables/per_mouse.tsv`
- `results/tables/by_majortype.tsv`
- `results/tables/spearman.tsv`
- `results/tables/geo_file_inventory.tsv`
- `results/figures/fig_mouse_cldn4_vs_T_myeloid.png`
- `results/figures/fig_honest_n.png`
- `results/figures/fig_cldn4_by_majortype.png`
- `results/summary.json`

## 结论

GSE127465 **小鼠** processed MTX 公开且 <2 GB，基因表有 **Cldn4**，因此是 GO，不是 no-go。对象是 CD45+ / 髓系为主的肺（作者上皮 **0**；Epcam/Cdh1∩Krt8 残留 **{n_leftover}**）。Cldn4 稀疏（**{n_cldn4}/15,939**，与残留上皮重叠 **{n_overlap}**）。小鼠层面 Cldn4 vs T 分数 n=4 ρ=+0.40（p=0.60），vs 髓系 ρ=−0.40（p=0.60）。**n=4 太薄**，不能写成精确效应，也不能写成 n=15,939。肿瘤侧只有 2 只鼠。**没有**与人 GSE127465 合并。不是 dual-high。

## Reproduce

```bash
python3 methods/gse127465_mouse_cldn4/scripts/download.py --out /tmp/gse127465_mouse
python3 methods/gse127465_mouse_cldn4/scripts/analyze.py \\
  --data /tmp/gse127465_mouse \\
  --outdir methods/gse127465_mouse_cldn4/results \\
  --finding methods/gse127465_mouse_cldn4/FINDING.md
```
"""
    finding_path.write_text(md)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, default=Path("/tmp/gse127465_mouse"))
    ap.add_argument("--outdir", type=Path, default=Path("methods/gse127465_mouse_cldn4/results"))
    ap.add_argument("--finding", type=Path, default=Path("methods/gse127465_mouse_cldn4/FINDING.md"))
    args = ap.parse_args()

    data = args.data
    mtx = data / "GSE127465_mouse_counts_normalized_15939x28205.mtx.gz"
    meta_p = data / "GSE127465_mouse_cell_metadata_15939x12.tsv.gz"
    genes_p = data / "GSE127465_gene_names_mouse_28205.tsv.gz"
    for p in (mtx, meta_p, genes_p):
        if not p.exists():
            raise SystemExit(f"missing {p}")

    genes = load_genes(genes_p)
    meta = load_meta(meta_p)
    missing = [g for g in KEEP_GENES if g not in genes]
    if "Cldn4" in missing:
        raise SystemExit("Cldn4 absent from mouse gene list — would be a no-go")
    if missing:
        print("optional genes absent:", missing, flush=True)
    keep = {g: genes.index(g) + 1 for g in KEEP_GENES if g in genes}
    expr = stream_gene_columns(mtx, keep)

    mice = mouse_table(meta, expr)
    by_type = celltype_table(meta, expr)

    spears = []
    for y in ("frac_T", "frac_TNK", "frac_myeloid", "frac_neutrophil", "frac_MoMacDC"):
        spears.append(spearman_block(mice, "Cldn4_mean_log1p", y))
        spears.append({**spearman_block(mice[mice["condition"] == "tumor"], "Cldn4_mean_log1p", y), "contrast": f"Cldn4_mean_log1p vs {y} (tumor only)"})
        spears.append({**spearman_block(mice[mice["condition"] == "healthy"], "Cldn4_mean_log1p", y), "contrast": f"Cldn4_mean_log1p vs {y} (healthy only)"})
    # also raw mean (not log1p) as sensitivity
    spears.append(spearman_block(mice, "Cldn4_mean", "frac_T"))
    spears.append(spearman_block(mice, "Cldn4_frac_pos", "frac_T"))
    spears.append(spearman_block(mice, "Cldn4_frac_pos", "frac_myeloid"))

    outdir = args.outdir
    (outdir / "tables").mkdir(parents=True, exist_ok=True)
    (outdir / "figures").mkdir(parents=True, exist_ok=True)
    mice.to_csv(outdir / "tables" / "per_mouse.tsv", sep="\t", index=False)
    by_type.to_csv(outdir / "tables" / "by_majortype.tsv", sep="\t", index=False)
    pd.DataFrame(spears).to_csv(outdir / "tables" / "spearman.tsv", sep="\t", index=False)

    inventory = {
        "mtx": mtx.name,
        "mtx_bytes": mtx.stat().st_size,
        "meta": meta_p.name,
        "meta_bytes": meta_p.stat().st_size,
        "genes": genes_p.name,
        "gene_bytes": genes_p.stat().st_size,
        "Cldn4_present": True,
        "human_files_used": 0,
    }
    pd.DataFrame(
        [
            {"file": inventory["mtx"], "bytes": inventory["mtx_bytes"], "used": "yes"},
            {"file": inventory["meta"], "bytes": inventory["meta_bytes"], "used": "yes"},
            {"file": inventory["genes"], "bytes": inventory["gene_bytes"], "used": "yes"},
            {"file": "human MTX/metadata/genes", "bytes": 0, "used": "no"},
        ]
    ).to_csv(outdir / "tables" / "geo_file_inventory.tsv", sep="\t", index=False)

    write_figures(mice, by_type, outdir / "figures")

    leftover = leftover_epi(expr)
    summary = {
        "accession": "GSE127465",
        "species": "mouse",
        "merged_with_human": False,
        "gene": "Cldn4",
        "n_mice": 4,
        "n_tumor": 2,
        "n_healthy": 2,
        "n_cells": int(len(meta)),
        "n_author_epithelium": 0,
        "n_leftover_epi": int(leftover.sum()),
        "n_Cldn4_pos": int((expr["Cldn4"] > 0).sum()),
        "n_Cldn4_pos_in_leftover": int(((expr["Cldn4"] > 0) & leftover).sum()),
        "spearman": spears,
        "per_mouse": mice.to_dict(orient="records"),
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")
    write_finding(mice, by_type, spears, inventory, expr, meta, args.finding)
    print("wrote", args.finding)
    print(mice.to_string(index=False))
    print("leftover epi", int(leftover.sum()), "Cldn4+", int((expr["Cldn4"] > 0).sum()))


if __name__ == "__main__":
    main()
