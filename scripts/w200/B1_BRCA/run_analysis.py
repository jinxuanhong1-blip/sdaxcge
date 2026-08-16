#!/usr/bin/env python3
"""B1 analog in TCGA-BRCA: is CLDN4 the top surface-gene partner of TACSTD2?

Honest ranking + robustness, matching the depth of the PAAD analog
(``results/w200/B1_PAAD``) while keeping the surfaceome-universe test that
defines claim B1.

Primary question
----------------
Among cell-surface genes (Bausch-Fluck 2018 in-silico surfaceome), does
CLDN4 rank #1 by Spearman co-expression with TACSTD2 (TROP2) in TCGA-BRCA?

Analyses
--------
1. Primary tumours (barcode sample-type ``01``), one sample per patient.
2. Surfaceome-wide Spearman + Pearson ranking of every surface gene vs
   TACSTD2 (the B1 test). Also rank CLDN4 transcriptome-wide.
3. Bootstrap 95% CI on the TACSTD2–CLDN4 Spearman rho.
4. Sensitivity: all primary vials (no patient collapse) — first-pass
   definition, so the published #4 rank stays comparable.
5. Sensitivity: PAM50Call_RNAseq strata (LumA / LumB / Her2 / Basal /
   Normal-like) — BRCA's analog of the PAAD neuroendocrine exclusion.
6. Sensitivity: rank-based partial correlation given ABSOLUTE purity.
7. Adjacent-normal (sample-type ``11``) TACSTD2–CLDN4 rho, as a check
   that the tumour number is not just a normal-epithelium program.
8. Claudin-family ranking (is CLDN4 even the top *claudin* partner?).

Outputs land in ``results/w200/B1_BRCA/``.
"""
from __future__ import annotations

import json
import os
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[3]
DATA_DIR = Path(os.environ.get("B1_BRCA_DATA", ROOT / "data"))
OUT_DIR = ROOT / "results" / "w200" / "B1_BRCA"
RNG = np.random.default_rng(20260816)

ANCHOR = "TACSTD2"
FOCUS = "CLDN4"
WINDOW = 200
CLAUDINS = [
    "CLDN1", "CLDN2", "CLDN3", "CLDN4", "CLDN5", "CLDN6", "CLDN7", "CLDN8",
    "CLDN9", "CLDN10", "CLDN11", "CLDN12", "CLDN14", "CLDN15", "CLDN16",
    "CLDN17", "CLDN18", "CLDN19", "CLDN20", "CLDN22", "CLDN23", "CLDN24",
    "CLDN25",
]
PAM50_ORDER = ["LumA", "LumB", "Her2", "Basal", "Normal"]
PAM50_COLORS = {
    "LumA": "#4c78a8",
    "LumB": "#f58518",
    "Her2": "#e45756",
    "Basal": "#72b7b2",
    "Normal": "#54a24b",
    "unknown": "#9e9e9e",
}


def load_expression(path: Path) -> pd.DataFrame:
    expr = pd.read_csv(path, sep="\t", index_col=0)
    expr.index = expr.index.astype(str)
    if expr.index.duplicated().any():
        expr = expr.groupby(level=0).mean()
    return expr


def load_surfaceome(path: Path) -> list[str]:
    surf = pd.read_excel(
        path, sheet_name="in silico surfaceome only", engine="openpyxl", header=1
    )
    genes = surf["UniProt gene"].dropna().astype(str)
    return sorted({g for g in genes if g not in ("nan", "None", "")})


def load_clinical(path: Path) -> pd.DataFrame:
    cl = pd.read_csv(path, sep="\t")
    cl = cl.rename(columns={"sampleID": "sample"})
    cl["patient"] = cl["sample"].astype(str).str[:12]
    return cl.set_index("sample")


def load_purity(path: Path) -> pd.Series:
    pur = pd.read_csv(path, sep="\t")
    pur["sample15"] = pur["array"].astype(str).str[:15]
    s = pur.set_index("sample15")["purity"]
    return s[~s.index.duplicated()]


def sample_type(barcode: str) -> str:
    parts = barcode.split("-")
    return parts[3][:2] if len(parts) >= 4 else ""


def one_per_patient(barcodes: list[str]) -> list[str]:
    """Keep the lexicographically first vial per patient (TCGA-XX-XXXX)."""
    first: dict[str, str] = {}
    for b in sorted(barcodes):
        first.setdefault(b[:12], b)
    return sorted(first.values())


def spearman_ci(x: np.ndarray, y: np.ndarray, n_boot: int = 2000) -> tuple[float, float]:
    n = len(x)
    boots = np.empty(n_boot)
    for i in range(n_boot):
        idx = RNG.integers(0, n, n)
        boots[i] = stats.spearmanr(x[idx], y[idx]).statistic
    return float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))


def partial_spearman(x, y, z) -> tuple[float, float]:
    """Spearman partial correlation of x,y given z (residualised ranks)."""
    rx, ry, rz = (stats.rankdata(v) for v in (x, y, z))

    def resid(a, b):
        beta = np.polyfit(b, a, 1)
        return a - np.polyval(beta, b)

    r, p = stats.pearsonr(resid(rx, rz), resid(ry, rz))
    return float(r), float(p)


def _corr_rows_vs_vector(mat: np.ndarray, vec: np.ndarray) -> np.ndarray:
    vc = vec - vec.mean()
    mc = mat - mat.mean(axis=1, keepdims=True)
    num = mc @ vc
    den = np.sqrt((mc**2).sum(axis=1) * (vc**2).sum())
    with np.errstate(invalid="ignore", divide="ignore"):
        return num / den


def _pvalues(r: np.ndarray, n: int) -> np.ndarray:
    r = np.clip(r, -0.999999999, 0.999999999)
    with np.errstate(invalid="ignore", divide="ignore"):
        t = r * np.sqrt((n - 2) / (1.0 - r**2))
    return 2 * stats.t.sf(np.abs(t), df=n - 2)


def _bh_fdr(p: np.ndarray) -> np.ndarray:
    q = np.full_like(p, np.nan, dtype=float)
    ok = ~np.isnan(p)
    pv = p[ok]
    m = pv.size
    order = np.argsort(pv)
    ranked = pv[order]
    adj = ranked * m / (np.arange(m) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.clip(adj, 0, 1)
    q[ok] = out
    return q


def rank_partners(expr: pd.DataFrame, universe: list[str], samples: list[str]) -> pd.DataFrame:
    sub = expr.loc[universe, samples]
    anchor_vec = expr.loc[ANCHOR, samples].to_numpy(dtype=float)
    mat = sub.to_numpy(dtype=float)
    n = len(samples)
    pear_r = _corr_rows_vs_vector(mat, anchor_vec)
    mat_rank = np.apply_along_axis(stats.rankdata, 1, mat)
    spear_r = _corr_rows_vs_vector(mat_rank, stats.rankdata(anchor_vec))
    df = pd.DataFrame(
        {
            "gene": universe,
            "spearman_r": spear_r,
            "spearman_p": _pvalues(spear_r, n),
            "pearson_r": pear_r,
            "pearson_p": _pvalues(pear_r, n),
            "mean_log2_expr": mat.mean(axis=1),
            "pct_expressed": (mat > 0).mean(axis=1) * 100.0,
        }
    )
    df["spearman_q"] = _bh_fdr(df["spearman_p"].to_numpy())
    df["pearson_q"] = _bh_fdr(df["pearson_p"].to_numpy())
    df = df.sort_values("spearman_r", ascending=False, kind="mergesort", na_position="last").reset_index(drop=True)
    df.insert(1, "spearman_rank", np.arange(1, len(df) + 1))
    df["pearson_rank"] = (
        df["pearson_r"].rank(ascending=False, method="min", na_option="bottom").astype("Int64")
    )
    return df


def corr_block(x: np.ndarray, y: np.ndarray, label: str, n_boot: int = 2000) -> dict:
    pr, pp = stats.pearsonr(x, y)
    sr, sp = stats.spearmanr(x, y)
    lo, hi = spearman_ci(x, y, n_boot)
    return {
        "subset": label,
        "n": int(len(x)),
        "pearson_r": round(float(pr), 4),
        "pearson_p": float(pp),
        "spearman_rho": round(float(sr), 4),
        "spearman_p": float(sp),
        "spearman_rho_ci95_lo": round(lo, 4),
        "spearman_rho_ci95_hi": round(hi, 4),
    }


def focus_row(df: pd.DataFrame, n_universe: int) -> dict:
    row = df[df["gene"] == FOCUS]
    if row.empty:
        return {"gene": FOCUS, "in_universe": False}
    r = row.iloc[0]
    s_rank = int(r["spearman_rank"])
    p_rank = int(r["pearson_rank"])
    return {
        "gene": FOCUS,
        "in_universe": True,
        "is_top_spearman": s_rank == 1,
        "is_top_pearson": p_rank == 1,
        "spearman_rank": s_rank,
        "spearman_r": float(r["spearman_r"]),
        "spearman_q": float(r["spearman_q"]),
        "spearman_percentile": round(100 * (1 - (s_rank - 1) / n_universe), 3),
        "pearson_rank": p_rank,
        "pearson_r": float(r["pearson_r"]),
        "pearson_q": float(r["pearson_q"]),
        "pearson_percentile": round(100 * (1 - (p_rank - 1) / n_universe), 3),
        "top1_spearman": str(df.iloc[0]["gene"]),
        "top1_spearman_r": float(df.iloc[0]["spearman_r"]),
    }


def make_plots(
    expr: pd.DataFrame,
    samples: list[str],
    pam50: pd.Series,
    surf_rank: pd.DataFrame,
    subtype_rho: pd.DataFrame,
    subtype_rank: pd.DataFrame,
    main_block: dict,
) -> None:
    x = expr.loc[ANCHOR, samples].to_numpy(dtype=float)
    y = expr.loc[FOCUS, samples].to_numpy(dtype=float)
    labels = pam50.reindex(samples).fillna("unknown")

    fig, ax = plt.subplots(figsize=(5.6, 5.2), dpi=150)
    for lab in list(PAM50_ORDER) + ["unknown"]:
        m = labels.values == lab
        if not m.any():
            continue
        ax.scatter(
            x[m], y[m], s=14, alpha=0.65, edgecolors="none",
            color=PAM50_COLORS.get(lab, "#9e9e9e"),
            label=f"{lab} (n={int(m.sum())})",
        )
    ax.set_xlabel("TACSTD2  log2(norm_count+1)")
    ax.set_ylabel("CLDN4  log2(norm_count+1)")
    ax.set_title(
        f"TCGA-BRCA primary tumors (n={main_block['n']}, 1/patient)\n"
        f"Spearman ρ = {main_block['spearman_rho']:.2f} "
        f"[{main_block['spearman_rho_ci95_lo']:.2f}, "
        f"{main_block['spearman_rho_ci95_hi']:.2f}]",
        fontsize=10,
    )
    ax.legend(fontsize=7, loc="lower right", frameon=False)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig1_scatter_tacstd2_cldn4.png")
    plt.close(fig)

    topn = surf_rank.head(20).iloc[::-1]
    colors = ["#d62728" if g == FOCUS else "#4c78a8" for g in topn["gene"]]
    fig, ax = plt.subplots(figsize=(6.2, 7.0), dpi=150)
    ax.barh(topn["gene"], topn["spearman_r"], color=colors)
    ax.set_xlabel("Spearman ρ with TACSTD2")
    ax.set_title("Top 20 surface-gene partners of TACSTD2\nTCGA-BRCA (CLDN4 in red)")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig2_top_surface_partners.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.4, 4.2), dpi=150)
    order = subtype_rho["subset"].tolist()
    y_pos = np.arange(len(order))
    ax.errorbar(
        subtype_rho["spearman_rho"],
        y_pos,
        xerr=[
            subtype_rho["spearman_rho"] - subtype_rho["spearman_rho_ci95_lo"],
            subtype_rho["spearman_rho_ci95_hi"] - subtype_rho["spearman_rho"],
        ],
        fmt="o",
        color="#2b6a99",
        capsize=3,
    )
    ax.axvline(0, color="#bbbbbb", lw=0.8)
    ax.set_yticks(y_pos)
    ax.set_yticklabels([f"{s} (n={n})" for s, n in zip(subtype_rho["subset"], subtype_rho["n"])])
    ax.set_xlabel("Spearman ρ  TACSTD2 × CLDN4  (bootstrap 95% CI)")
    ax.set_title("TACSTD2–CLDN4 co-expression by PAM50 / tissue")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig3_subtype_rho.png")
    plt.close(fig)

    if not subtype_rank.empty:
        fig, ax = plt.subplots(figsize=(6.0, 4.0), dpi=150)
        sr = subtype_rank.sort_values("cldn4_surface_rank")
        ax.barh(sr["subset"], sr["cldn4_surface_rank"], color="#4c78a8")
        ax.axvline(1, color="#d62728", ls="--", lw=1, label="rank #1")
        ax.set_xlabel("CLDN4 Spearman rank among surface genes (lower = closer to top)")
        ax.set_title("Does CLDN4 become #1 inside any PAM50 stratum?")
        ax.invert_yaxis()
        ax.legend(fontsize=8, frameon=False)
        fig.tight_layout()
        fig.savefig(OUT_DIR / "fig4_cldn4_rank_by_subtype.png")
        plt.close(fig)


def write_readme(payload: dict) -> None:
    fr = payload["focus_surfaceome_one_per_patient"]
    main = next(r for r in payload["tacstd2_cldn4_rho"] if r["subset"] == "primary_one_per_patient")
    allv = next(r for r in payload["tacstd2_cldn4_rho"] if r["subset"] == "primary_all_vials")
    adj = next((r for r in payload["tacstd2_cldn4_rho"] if r["subset"] == "adjacent_normal"), None)
    part = payload["purity"]
    tw = payload["transcriptome_wide"]
    lines = [
        "# B1 analog — TCGA-BRCA: is CLDN4 the top surface partner of TACSTD2?",
        "",
        "Honest ranking of every surfaceome gene against `TACSTD2` in TCGA-BRCA,",
        "plus the robustness checks used in the PAAD analog (bootstrap CI, one",
        "sample per patient, purity-partial correlation, subtype strata).",
        "",
        "## TL;DR",
        "",
        "**中文**：在 TCGA-BRCA 原发肿瘤（每患者一个 -01 样本，n="
        f"{main['n']}）中，CLDN4 **不是** TACSTD2 的第一表面共表达伙伴。",
        f"按 Spearman，CLDN4 在 {payload['n_surface_universe']} 个表面组基因中排名",
        f"**第 {fr['spearman_rank']}**（ρ = {fr['spearman_r']:.3f}，FDR q ≈ "
        f"{fr['spearman_q']:.1e}，{fr['spearman_percentile']} 百分位）。",
        f"排在它前面的是 {', '.join(payload['genes_ahead_of_cldn4'])}。",
        f"TACSTD2–CLDN4 本身的 Spearman ρ = **{main['spearman_rho']:.3f}**",
        f"（bootstrap 95% CI {main['spearman_rho_ci95_lo']:.2f}–"
        f"{main['spearman_rho_ci95_hi']:.2f}）。全转录组排名第 "
        f"**{tw['cldn4_rank']}** / {tw['n_expressed_genes']}。",
        "PAM50 分层后，CLDN4 在任何一个亚型里都没有变成 #1（最接近的是 Basal，第 2）。",
        "ABSOLUTE 纯度偏相关几乎不改变 ρ。",
        "癌旁正常组织的 ρ（0.80）**高于**肿瘤，提示这是正常上皮程序在肿瘤中被稀释，而非肿瘤特异耦合。",
        "",
        f"**English**: In TCGA-BRCA primary tumours (n={main['n']}, one -01",
        f"sample per patient), CLDN4 is **not** the top surface-gene partner of",
        f"TACSTD2. It ranks **#{fr['spearman_rank']} of {payload['n_surface_universe']}**",
        f"by Spearman (ρ = {fr['spearman_r']:.3f}, FDR q ≈ {fr['spearman_q']:.1e},",
        f"{fr['spearman_percentile']}th percentile). Genes ahead of it:",
        f"{', '.join(payload['genes_ahead_of_cldn4'])}.",
        f"The TACSTD2–CLDN4 pair itself is ρ = **{main['spearman_rho']:.3f}**",
        f"(bootstrap 95% CI {main['spearman_rho_ci95_lo']:.2f}–"
        f"{main['spearman_rho_ci95_hi']:.2f}). Transcriptome-wide, CLDN4 ranks",
        f"**#{tw['cldn4_rank']} of {tw['n_expressed_genes']}** expressed genes.",
        "CLDN4 does not become #1 inside any PAM50 stratum (closest: Basal, #2).",
        "Rank-based partial correlation given ABSOLUTE purity leaves ρ unchanged.",
        "Adjacent-normal ρ is *higher* (0.80) than the tumour pair — the coupling",
        "looks like a normal-epithelium program diluted in tumours, not a",
        "tumour-specific link.",
        "",
        f"Headline vs first pass: all-vial ranking (n={allv['n']}) was also",
        f"#{payload['focus_surfaceome_all_vials']['spearman_rank']}. This HiSeqV2",
        "extract already has one primary vial per patient, so collapsing is a no-op.",
        "",
        "## Data (all open access)",
        "",
        "| File | Source |",
        "| --- | --- |",
        "| `TCGA-BRCA.HiSeqV2.gz` | [Xena TCGA hub HiSeqV2](https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.BRCA.sampleMap%2FHiSeqV2.gz) — log2(norm_count+1), HGNC symbols |",
        "| `BRCA_clinicalMatrix` | [Xena BRCA clinicalMatrix](https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.BRCA.sampleMap%2FBRCA_clinicalMatrix) — `PAM50Call_RNAseq` |",
        "| `tcga_absolute_purity.txt` | [GDC open ABSOLUTE](https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5) |",
        "| `table_S3_surfaceome.xlsx` | Bausch-Fluck 2018, GitHub mirror `steveneschrich/surfaceome` |",
        "",
        "## Methods",
        "",
        "- Primary tumours only (barcode sample-type `01`); main analysis keeps",
        "  the lexicographically first vial per patient.",
        "- Surface-gene universe: in-silico human surfaceome (Bausch-Fluck et al.,",
        f"  *PNAS* 2018), intersected with the matrix, anchor removed → "
        f"{payload['n_surface_universe']} genes.",
        "- Primary metric: Spearman; Pearson reported alongside. 2,000-fold",
        "  case-resampling bootstrap (seed 20260816) for the pair CI.",
        "- Transcriptome-wide rank: genes with log2(norm+1) > 1 in ≥ 20% of",
        "  tumours, duplicate symbols collapsed to the highest-mean row.",
        "- PAM50 strata use `PAM50Call_RNAseq` (better coverage than the 2012",
        "  Nature PAM50 column).",
        "- Purity: rank-based partial Spearman of TACSTD2 vs CLDN4 given",
        "  ABSOLUTE purity.",
        "",
        "## Results",
        "",
        "### TACSTD2 × CLDN4 pair",
        "",
        "| subset | n | Spearman ρ [95% CI] | Pearson r |",
        "| --- | --- | --- | --- |",
    ]
    for r in payload["tacstd2_cldn4_rho"]:
        lines.append(
            f"| {r['subset']} | {r['n']} | **{r['spearman_rho']:.3f}** "
            f"[{r['spearman_rho_ci95_lo']:.2f}, {r['spearman_rho_ci95_hi']:.2f}] "
            f"| {r['pearson_r']:.3f} |"
        )
    lines += [
        "",
        "### Surfaceome ranking of CLDN4 (the B1 test)",
        "",
        "| analysis | CLDN4 Spearman rank | ρ | #1 gene |",
        "| --- | --- | --- | --- |",
        f"| one-per-patient (main) | **#{fr['spearman_rank']} / {payload['n_surface_universe']}** "
        f"| {fr['spearman_r']:.3f} | {fr['top1_spearman']} |",
        f"| all primary vials | #{payload['focus_surfaceome_all_vials']['spearman_rank']} / "
        f"{payload['n_surface_universe']} | {payload['focus_surfaceome_all_vials']['spearman_r']:.3f} "
        f"| {payload['focus_surfaceome_all_vials']['top1_spearman']} |",
    ]
    for r in payload["cldn4_rank_by_pam50"]:
        lines.append(
            f"| {r['subset']} | #{r['cldn4_surface_rank']} / {r['n_universe']} "
            f"| {r['cldn4_spearman_r']:.3f} | {r['top1_gene']} |"
        )
    lines += [
        "",
        "### Transcriptome-wide rank",
        "",
        f"CLDN4 is the **{tw['cldn4_rank']}th** strongest TACSTD2 partner of",
        f"{tw['n_expressed_genes']} expressed genes ({tw['cldn4_percentile']}th",
        "percentile). The genes immediately ahead of it are themselves",
        "epithelial / junctional / keratin-program genes — a shared",
        "malignant-epithelial module, not a CLDN4-unique link. See",
        "`top25_tacstd2_partners_allgenes.csv`.",
        "",
        "### Claudin family",
        "",
        f"Among measured claudins, CLDN4 ranks **#{payload['claudin_family']['cldn4_rank_among_claudins']}**",
        f"as a TACSTD2 partner (top claudin: {payload['claudin_family']['top_claudin']}).",
        "See `claudin_family_rank.csv`.",
        "",
        "### Purity",
        "",
        "| quantity | n | ρ |",
        "| --- | --- | --- |",
    ]
    for r in part:
        lines.append(f"| {r['quantity']} | {r['n']} | {r['rho']:.3f} |")
    if adj is not None:
        lines += [
            "",
            "### Adjacent normal",
            "",
            f"In {adj['n']} adjacent-normal samples, TACSTD2–CLDN4 Spearman ρ =",
            f"**{adj['spearman_rho']:.3f}** [{adj['spearman_rho_ci95_lo']:.2f},",
            f"{adj['spearman_rho_ci95_hi']:.2f}] — *higher* than the tumour pair",
            f"({main['spearman_rho']:.3f}). The tumour number is therefore a",
            "weaker echo of a normal-epithelium program, not a tumour-acquired",
            "coupling. Residual normal epithelium cannot be invoked to *inflate*",
            "the tumour ρ; if anything it would dilute it.",
        ]
    lines += [
        "",
        "### Comparison to other B1 analogs (Spearman ρ, primary tumours)",
        "",
        "| cohort | TACSTD2–CLDN4 ρ | CLDN4 surfaceome rank | source |",
        "| --- | --- | --- | --- |",
        "| TCGA-LUAD | 0.53 | (pair only) | `results/fable_tcga` |",
        "| TCGA-LUSC | 0.39 | (pair only) | `results/fable_tcga` |",
        "| TCGA-PAAD | 0.71 | 13th of 20,282 *all* genes | `results/w200/B1_PAAD` |",
        f"| **TCGA-BRCA** | **{main['spearman_rho']:.2f}** | "
        f"**#{fr['spearman_rank']} of {payload['n_surface_universe']} surface** / "
        f"#{tw['cldn4_rank']} of {tw['n_expressed_genes']} all | this analysis |",
        "",
        "BRCA's pair-wise ρ is *weaker* than PAAD and LUAD and close to LUSC;",
        "the surfaceome rank (#4) is still near-top, just not #1.",
        "",
        "## Honest caveats",
        "",
        "1. **CLDN4 is near-top, not top.** EFNA1 is a clear #1 (ρ ≈ 0.43);",
        "   PVRL4 (= NECTIN4) and EFNA4 sit just above CLDN4. Quoting CLDN4",
        "   as *the* top TACSTD2 surface partner in BRCA is not literally true.",
        "2. **Shared epithelial program.** Genes ahead of and around CLDN4",
        "   (ephrins, nectin-4, TM4SF1, ITGB4, CD151) are themselves",
        "   epithelial / junctional. This is a module, not a CLDN4-specific",
        "   coupling.",
        "3. **Bulk mRNA ≠ protein, ≠ single cell.** ABSOLUTE (DNA) purity",
        "   adjustment barely moves ρ and neither gene tracks DNA purity",
        "   strongly, but DNA purity is an imperfect proxy for the epithelial",
        "   mRNA fraction. Compartment confounding is reduced, not eliminated.",
        "4. **PAM50 does not rescue the #1 claim.** CLDN4's surfaceome rank",
        "   stays off #1 in LumA (#23) / LumB (#54) / Her2 (#28) / Basal (#2) /",
        "   Normal-like (#28). Basal is the closest (behind SLC39A2) and has",
        "   the strongest tumour-pair ρ (0.47). PAM50 Normal-like n=23 has a",
        "   wide CI and should not be over-read.",
        "5. **Adjacent-normal ρ > tumour ρ.** The pair is tighter in 114",
        "   adjacent-normal samples (ρ = 0.80) than in tumours (ρ = 0.35).",
        "   That is the opposite of a tumour-specific coupling.",
        "6. **HiSeqV2 is gene-level RNA-seq** (log2 norm_count+1), not the",
        "   GDC STAR TPM matrix used in the PAAD analog. Pair-wise ρ is",
        "   therefore not strictly interchangeable across those two matrices;",
        "   the *rank* answer (not #1) is the robust claim.",
        "7. TCGA-BRCA is a resected, mostly untreated cohort. Nothing here",
        "   speaks to TROP2-ADC or ICI response.",
        "",
        "## Outputs",
        "",
        "| file | contents |",
        "| --- | --- |",
        "| `coexpression_TACSTD2_surfaceome.csv` | full surfaceome ranking (one-per-patient) |",
        "| `top200.csv` | top-`w200` window |",
        "| `top25_tacstd2_partners_allgenes.csv` | transcriptome-wide top 25 |",
        "| `tacstd2_cldn4_rho.csv` | pair-wise ρ across subsets |",
        "| `purity_adjusted.csv` | ABSOLUTE partial correlation |",
        "| `cldn4_rank_by_pam50.csv` | surfaceome rank of CLDN4 inside each PAM50 |",
        "| `claudin_family_rank.csv` | TACSTD2 vs every measured claudin |",
        "| `summary.json` / `summary.md` | machine + short human verdict |",
        "| `fig1_scatter_tacstd2_cldn4.png` | PAM50-coloured scatter |",
        "| `fig2_top_surface_partners.png` | top-20 surface bar chart |",
        "| `fig3_subtype_rho.png` | forest plot of pair ρ |",
        "| `fig4_cldn4_rank_by_subtype.png` | CLDN4 rank by PAM50 |",
        "",
        "## Reproduce",
        "",
        "```bash",
        "pip install -r requirements.txt",
        "python3 scripts/w200/B1_BRCA/download_data.py   # → data/  (or $B1_BRCA_DATA)",
        "python3 scripts/w200/B1_BRCA/run_analysis.py    # → results/w200/B1_BRCA/",
        "```",
        "",
        f"Generated {payload['generated_utc']}.",
        "",
    ]
    (OUT_DIR / "README.md").write_text("\n".join(lines))


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    expr = load_expression(DATA_DIR / "TCGA-BRCA.HiSeqV2.gz")
    surf_genes = load_surfaceome(DATA_DIR / "table_S3_surfaceome.xlsx")
    clin = load_clinical(DATA_DIR / "BRCA_clinicalMatrix")
    purity = load_purity(DATA_DIR / "tcga_absolute_purity.txt")

    if ANCHOR not in expr.index or FOCUS not in expr.index:
        raise SystemExit(f"{ANCHOR} or {FOCUS} missing from expression matrix")

    primary_all = [c for c in expr.columns if sample_type(c) == "01"]
    primary = one_per_patient(primary_all)
    adjacent = one_per_patient([c for c in expr.columns if sample_type(c) == "11"])

    pam50 = clin["PAM50Call_RNAseq"]
    pam50 = pam50.reindex(expr.columns)

    universe = [g for g in surf_genes if g in expr.index and g != ANCHOR]
    n_universe = len(universe)

    print(f"[data] genes x samples = {expr.shape[0]} x {expr.shape[1]}")
    print(f"[data] primary all-vials={len(primary_all)}  one/patient={len(primary)}  adjacent={len(adjacent)}")
    print(f"[data] surfaceome universe={n_universe}")

    # --- pair-wise TACSTD2 x CLDN4 -----------------------------------------
    pair_rows = []
    x_main = expr.loc[ANCHOR, primary].to_numpy(dtype=float)
    y_main = expr.loc[FOCUS, primary].to_numpy(dtype=float)
    pair_rows.append(corr_block(x_main, y_main, "primary_one_per_patient"))
    pair_rows.append(
        corr_block(
            expr.loc[ANCHOR, primary_all].to_numpy(dtype=float),
            expr.loc[FOCUS, primary_all].to_numpy(dtype=float),
            "primary_all_vials",
        )
    )
    if adjacent:
        pair_rows.append(
            corr_block(
                expr.loc[ANCHOR, adjacent].to_numpy(dtype=float),
                expr.loc[FOCUS, adjacent].to_numpy(dtype=float),
                "adjacent_normal",
            )
        )

    subtype_rank_rows = []
    for lab in PAM50_ORDER:
        cols = [c for c in primary if str(pam50.get(c, "")) == lab]
        if len(cols) < 20:
            continue
        pair_rows.append(
            corr_block(
                expr.loc[ANCHOR, cols].to_numpy(dtype=float),
                expr.loc[FOCUS, cols].to_numpy(dtype=float),
                f"PAM50_{lab}",
            )
        )
        rdf = rank_partners(expr, universe, cols)
        fr = focus_row(rdf, len(universe))
        subtype_rank_rows.append(
            {
                "subset": f"PAM50_{lab}",
                "n": len(cols),
                "n_universe": len(universe),
                "cldn4_surface_rank": fr.get("spearman_rank"),
                "cldn4_spearman_r": round(fr.get("spearman_r", float("nan")), 4),
                "top1_gene": fr.get("top1_spearman"),
                "top1_r": round(fr.get("top1_spearman_r", float("nan")), 4),
            }
        )

    pair_df = pd.DataFrame(pair_rows)
    pair_df.to_csv(OUT_DIR / "tacstd2_cldn4_rho.csv", index=False)
    subtype_rank_df = pd.DataFrame(subtype_rank_rows)
    subtype_rank_df.to_csv(OUT_DIR / "cldn4_rank_by_pam50.csv", index=False)

    # --- surfaceome rankings -----------------------------------------------
    surf_main = rank_partners(expr, universe, primary)
    surf_all = rank_partners(expr, universe, primary_all)
    surf_main.to_csv(OUT_DIR / "coexpression_TACSTD2_surfaceome.csv", index=False)
    surf_main.head(WINDOW).to_csv(OUT_DIR / "top200.csv", index=False)
    focus_main = focus_row(surf_main, n_universe)
    focus_all = focus_row(surf_all, n_universe)
    genes_ahead = surf_main.loc[surf_main["spearman_rank"] < focus_main["spearman_rank"], "gene"].tolist()

    # --- transcriptome-wide rank -------------------------------------------
    mat = expr.loc[:, primary]
    expressed = mat[(mat > 1.0).mean(axis=1) >= 0.20]
    expressed = expressed[expressed.index != ANCHOR]
    expressed = expressed.copy()
    expressed["_mean"] = expressed.mean(axis=1)
    expressed = expressed.sort_values("_mean", ascending=False).drop(columns="_mean")
    expressed = expressed[~expressed.index.duplicated()]
    xr = stats.rankdata(x_main)
    xr_c = xr - xr.mean()
    m_ranks = np.apply_along_axis(stats.rankdata, 1, expressed.to_numpy(dtype=float))
    m_c = m_ranks - m_ranks.mean(axis=1, keepdims=True)
    num = m_c @ xr_c
    den = np.sqrt((m_c**2).sum(axis=1) * (xr_c**2).sum())
    rhos = pd.Series(num / den, index=expressed.index).sort_values(ascending=False)
    cldn4_tw_rank = int(rhos.index.get_loc(FOCUS)) + 1
    top_all = rhos.head(25).round(4).reset_index()
    top_all.columns = ["gene", "spearman_rho_vs_TACSTD2"]
    top_all.insert(0, "rank", range(1, len(top_all) + 1))
    top_all.to_csv(OUT_DIR / "top25_tacstd2_partners_allgenes.csv", index=False)
    tw = {
        "n_expressed_genes": int(len(rhos)),
        "cldn4_rank": cldn4_tw_rank,
        "cldn4_percentile": round(100 * (1 - (cldn4_tw_rank - 1) / len(rhos)), 3),
        "cldn4_rho": float(rhos.loc[FOCUS]),
        "top1_gene": str(rhos.index[0]),
        "top1_rho": float(rhos.iloc[0]),
    }

    # --- claudin family ----------------------------------------------------
    cldn_present = [g for g in CLAUDINS if g in expr.index]
    cldn_rank = rank_partners(expr, cldn_present, primary)
    cldn_rank.to_csv(OUT_DIR / "claudin_family_rank.csv", index=False)
    cldn_focus = cldn_rank[cldn_rank["gene"] == FOCUS].iloc[0]
    claudin_family = {
        "n_claudins_measured": len(cldn_present),
        "cldn4_rank_among_claudins": int(cldn_focus["spearman_rank"]),
        "cldn4_r": float(cldn_focus["spearman_r"]),
        "top_claudin": str(cldn_rank.iloc[0]["gene"]),
        "top_claudin_r": float(cldn_rank.iloc[0]["spearman_r"]),
    }

    # --- purity ------------------------------------------------------------
    pvec = purity.reindex([c[:15] for c in primary]).to_numpy(dtype=float)
    mask = ~np.isnan(pvec)
    r_part, p_part = partial_spearman(x_main[mask], y_main[mask], pvec[mask])
    sr_sub, sp_sub = stats.spearmanr(x_main[mask], y_main[mask])
    purity_rows = [
        {
            "quantity": "spearman_unadjusted_purity_subset",
            "n": int(mask.sum()),
            "rho": round(float(sr_sub), 4),
            "p": float(sp_sub),
        },
        {
            "quantity": "spearman_partial_given_ABSOLUTE_purity",
            "n": int(mask.sum()),
            "rho": round(r_part, 4),
            "p": p_part,
        },
        {
            "quantity": "spearman_TACSTD2_vs_purity",
            "n": int(mask.sum()),
            "rho": round(float(stats.spearmanr(x_main[mask], pvec[mask]).statistic), 4),
            "p": float(stats.spearmanr(x_main[mask], pvec[mask]).pvalue),
        },
        {
            "quantity": "spearman_CLDN4_vs_purity",
            "n": int(mask.sum()),
            "rho": round(float(stats.spearmanr(y_main[mask], pvec[mask]).statistic), 4),
            "p": float(stats.spearmanr(y_main[mask], pvec[mask]).pvalue),
        },
    ]
    pd.DataFrame(purity_rows).to_csv(OUT_DIR / "purity_adjusted.csv", index=False)

    # --- plots -------------------------------------------------------------
    subtype_rho_plot = pair_df[
        pair_df["subset"].isin(
            ["primary_one_per_patient", "adjacent_normal"] + [f"PAM50_{s}" for s in PAM50_ORDER]
        )
    ].copy()
    make_plots(
        expr, primary, pam50, surf_main, subtype_rho_plot, subtype_rank_df, pair_rows[0]
    )

    # --- summaries ---------------------------------------------------------
    verdict_yes = bool(focus_main.get("is_top_spearman"))
    summary = {
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "cohort": "TCGA-BRCA",
        "expression_dataset": "UCSC Xena TCGA.BRCA.sampleMap/HiSeqV2 (log2 norm_count+1)",
        "surfaceome_source": "Bausch-Fluck et al. 2018 in-silico surfaceome (table S3)",
        "anchor_gene": ANCHOR,
        "focus_gene": FOCUS,
        "n_primary_all_vials": len(primary_all),
        "n_primary_one_per_patient": len(primary),
        "n_adjacent_normal": len(adjacent),
        "n_with_ABSOLUTE_purity": int(mask.sum()),
        "n_surface_universe": n_universe,
        "correlation_primary": "spearman",
        "window": WINDOW,
        "question": "Is CLDN4 the top surface-gene co-expression partner of TACSTD2 in TCGA-BRCA?",
        "answer": "YES" if verdict_yes else "NO",
        "focus_surfaceome_one_per_patient": focus_main,
        "focus_surfaceome_all_vials": focus_all,
        "genes_ahead_of_cldn4": genes_ahead,
        "tacstd2_cldn4_rho": pair_rows,
        "purity": purity_rows,
        "transcriptome_wide": tw,
        "claudin_family": claudin_family,
        "cldn4_rank_by_pam50": subtype_rank_rows,
        "top10_surface_spearman": surf_main.head(10)[
            ["gene", "spearman_rank", "spearman_r", "spearman_q"]
        ].to_dict("records"),
    }
    (OUT_DIR / "summary.json").write_text(json.dumps(summary, indent=2))

    md = [
        f"# {FOCUS} vs {ANCHOR} co-expression in TCGA-BRCA (honest ranking)",
        "",
        f"**Question:** Is `{FOCUS}` the top co-expression partner of `{ANCHOR}` among surface genes?",
        "",
        f"**Answer:** {'YES' if verdict_yes else 'NO'} — {FOCUS} ranks "
        f"**#{focus_main['spearman_rank']} of {n_universe}** surface genes "
        f"(Spearman ρ = {focus_main['spearman_r']:.3f}). "
        f"#{1} is {focus_main['top1_spearman']}. Genes ahead of {FOCUS}: "
        f"{', '.join(genes_ahead) if genes_ahead else '(none)'}.",
        "",
        f"- Main cohort: {len(primary)} primary tumours, one sample per patient "
        f"(all-vial n = {len(primary_all)} gives the same rank "
        f"#{focus_all['spearman_rank']}).",
        f"- Pair-wise TACSTD2–CLDN4 Spearman ρ = {pair_rows[0]['spearman_rho']:.3f} "
        f"[{pair_rows[0]['spearman_rho_ci95_lo']:.2f}, {pair_rows[0]['spearman_rho_ci95_hi']:.2f}].",
        f"- Transcriptome-wide rank: #{tw['cldn4_rank']} of {tw['n_expressed_genes']}.",
        f"- Claudin-family rank: #{claudin_family['cldn4_rank_among_claudins']} "
        f"(top claudin = {claudin_family['top_claudin']}).",
        f"- Purity-partial ρ = {r_part:.3f} (unadjusted on the same n = {int(mask.sum())}: "
        f"{float(sr_sub):.3f}).",
        "",
        "See `README.md` in this folder for the full bilingual write-up.",
        "",
        "## Top 15 surface-gene partners (one-per-patient, Spearman)",
        "",
        "| rank | gene | spearman_rho | pearson_rho | FDR q |",
        "| --- | --- | --- | --- | --- |",
    ]
    for _, r in surf_main.head(15).iterrows():
        star = "  <-- focus" if r["gene"] == FOCUS else ""
        md.append(
            f"| {int(r['spearman_rank'])} | {r['gene']}{star} | "
            f"{r['spearman_r']:.3f} | {r['pearson_r']:.3f} | {r['spearman_q']:.2e} |"
        )
    (OUT_DIR / "summary.md").write_text("\n".join(md) + "\n")

    write_readme(summary)

    # drop first-pass figures that the new names supersede
    for stale in (
        "scatter_TACSTD2_vs_CLDN4.png",
        "top_partners_TACSTD2.png",
    ):
        p = OUT_DIR / stale
        if p.exists():
            p.unlink()

    print(json.dumps({k: summary[k] for k in (
        "answer", "n_primary_one_per_patient", "n_surface_universe",
        "focus_surfaceome_one_per_patient", "genes_ahead_of_cldn4",
        "transcriptome_wide", "claudin_family",
    )}, indent=2))
    print(pair_df.to_string(index=False))
    print(pd.DataFrame(purity_rows).to_string(index=False))
    print(subtype_rank_df.to_string(index=False))
    print("\nTop 10 surface partners:")
    print(surf_main.head(10)[["spearman_rank", "gene", "spearman_r"]].to_string(index=False))
    print(f"\n[done] {OUT_DIR}")


if __name__ == "__main__":
    main()
