#!/usr/bin/env python3
"""Fill REPORT.md from the grid TSVs after analyze.py."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent.parent
RHO_CUT = -0.35


def fmt_rho(r) -> str:
    if r is None or (isinstance(r, float) and pd.isna(r)):
        return "NA"
    return f"{float(r):+.3f}"


def fmt_p(p) -> str:
    if p is None or (isinstance(p, float) and pd.isna(p)):
        return "NA"
    p = float(p)
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.3f}"


def fmt_num(x, nd=3) -> str:
    if x is None or (isinstance(x, float) and pd.isna(x)):
        return "NA"
    return f"{float(x):.{nd}f}"


def main() -> None:
    spear = pd.read_csv(HERE / "grid_spearman.tsv", sep="\t")
    nmpr = pd.read_csv(HERE / "grid_nmpr_mpr.tsv", sep="\t")
    per = pd.read_csv(HERE / "per_sample.tsv", sep="\t")
    summary = json.loads((HERE / "summary.json").read_text())
    post = per[per["is_post"]].copy()

    hi = spear[spear["highlight_rho_le_neg035"] == True]  # noqa: E712
    hi_n = nmpr[nmpr["nmpr_gt_mpr"] == True]  # noqa: E712

    lines = []
    a = lines.append
    a("# GSE207422 public-UMI definition grid")
    a("")
    a("ADDITIVE extra. Public GEO UMI only. A3 slide numbers taken as given (malignant TACSTD2 higher in NMPR than MPR; per-patient Spearman vs T/NK ρ ≈ −0.40 to −0.50). This slice does not re-litigate the slide. Author CopyKAT / epithelium RDS barcodes are **not public** and are not in the numeric grid.")
    a("")
    a(f"Matrix: **{summary['n_cells_public_matrix']:,}** cells × **{summary['n_genes_matrix']:,}** genes (user PPT ~90,652; deposited matrix is 92,330). Unit = **12 post-treatment patients** (MPR n=4 including pCR P06; NMPR n=8). Extra figure: `fig_def_grid.png`. Gold boxes = ρ ≤ −0.35 or NMPR>MPR.")
    a("")
    a("## Cell counts by definition")
    a("")
    a("| Axis | Definition | n cells |")
    a("|---|---|---|")
    a(f"| lineage | Hu-mean epithelial | {summary['n_epithelial_hu_mean']:,} |")
    a(f"| lineage | Hu-%pos epithelial | {summary['n_epithelial_hu_pctpos']:,} |")
    a(f"| lineage | stromal (fibroblast+endothelial) | {summary['n_stromal']:,} |")
    a(f"| immune | T | {summary['n_T']:,} |")
    a(f"| immune | NK | {summary['n_NK']:,} |")
    a(f"| immune | CD8 only | {summary['n_CD8']:,} |")
    a(f"| immune | CXCL13+ T/NK | {summary['n_CXCL13_TNK']:,} |")
    a(f"| immune | cytotoxicity-high T | {summary['n_cyto_high_T']:,} |")
    for k, v in summary["malignant_n"].items():
        a(f"| malignant | {k} | {int(v):,} |")
    a("| malignant | author CopyKAT | not public |")
    a("")
    a("UCell module = {TACSTD2, CLDN4, EPCAM}, maxRank=1500. CopyKAT-like / inferCNV-like use stromal reference; they are window-smoothed expression scores, not the R packages.")
    a("")
    a("## Highlight: ρ ≤ −0.35 vs T+NK")
    a("")
    hi_tnk = hi[hi["immune_def"] == "T_NK"] if not hi.empty else hi
    if hi_tnk.empty:
        a("No public combination has Spearman ρ ≤ −0.35 vs T+NK.")
    else:
        a("Requested flag: TACSTD2 vs T/NK fraction, ρ ≤ −0.35. n=12 keeps all four MPR samples; n<12 drops MPR samples with zero malignant cells.")
        a("")
        a("| Malignant | Score | ρ | p | n |")
        a("|---|---|---|---|---|")
        for _, r in hi_tnk.sort_values(["n", "rho"], ascending=[False, True]).iterrows():
            a(
                f"| {r['malignant_label']} | {r['score_label']} | "
                f"**{fmt_rho(r['rho'])}** | {fmt_p(r['p'])} | {int(r['n'])} |"
            )
    a("")
    a("## Highlight: ρ ≤ −0.35 (all five immune defs)")
    a("")
    if hi.empty:
        a("No public combination on this grid has Spearman ρ ≤ −0.35 vs any of the five immune fractions.")
    else:
        a(f"{len(hi)} combination(s). Strongest n=12 rows are EPCAM+KRT (pos) or DRMref public UCell vs CXCL13+ / cytotoxicity-high T.")
        a("")
        a("| Malignant | Immune | Score | ρ | p | n |")
        a("|---|---|---|---|---|---|")
        for _, r in hi.sort_values("rho").iterrows():
            a(
                f"| {r['malignant_label']} | {r['immune_label']} | {r['score_label']} | "
                f"**{fmt_rho(r['rho'])}** | {fmt_p(r['p'])} | {int(r['n'])} |"
            )
    a("")
    a("## Highlight: NMPR>MPR (Δ > 0)")
    a("")
    if hi_n.empty:
        a("No public combination has NMPR mean > MPR mean.")
    else:
        a(f"{len(hi_n)} combination(s). Exact Wilcoxon. None of these is a retuning of the slide.")
        a("")
        a("| Malignant | Score | mean NMPR | mean MPR | Δ | U | p | n (NMPR vs MPR) |")
        a("|---|---|---|---|---|---|---|---|")
        for _, r in hi_n.sort_values("p").iterrows():
            a(
                f"| {r['malignant_label']} | {r['score_label']} | {fmt_num(r['mean_NMPR'])} | "
                f"{fmt_num(r['mean_MPR'])} | {fmt_num(r['delta_NMPR_minus_MPR'])} | "
                f"{fmt_num(r['U'], 1)} | {fmt_p(r['p'])} | {int(r['n_NMPR'])} vs {int(r['n_MPR'])} |"
            )
    a("")
    a("## Full Spearman grid (every ρ / p / n)")
    a("")
    a("Post-treatment patients only. Samples with zero malignant cells under that definition drop out (n < 12).")
    a("")
    for score, sg in spear.groupby("score", sort=False):
        label = sg.iloc[0]["score_label"]
        a(f"### Score: {label}")
        a("")
        a("| Malignant | T+NK ρ (p, n) | CD8 only ρ (p, n) | CD8+NK ρ (p, n) | CXCL13+ ρ (p, n) | cyto-high T ρ (p, n) |")
        a("|---|---|---|---|---|---|")
        for mal, mg in sg.groupby("malignant_def", sort=False):
            cells = []
            lab = mg.iloc[0]["malignant_label"]
            for imm in ["T_NK", "CD8_only", "CD8_NK", "CXCL13_pos", "cyto_high_T"]:
                hit = mg[mg["immune_def"] == imm]
                if hit.empty:
                    cells.append("NA")
                    continue
                r = hit.iloc[0]
                mark = "**" if bool(r["highlight_rho_le_neg035"]) else ""
                cells.append(f"{mark}{fmt_rho(r['rho'])}{mark} ({fmt_p(r['p'])}, n={int(r['n'])})")
            a(f"| {lab} | " + " | ".join(cells) + " |")
        a("")
    a("## Full NMPR vs MPR grid (every Δ / p / n)")
    a("")
    a("| Malignant | Score | mean NMPR | mean MPR | Δ | U | p | p one-sided NMPR>MPR | n |")
    a("|---|---|---|---|---|---|---|---|---|")
    for _, r in nmpr.iterrows():
        star = "**" if bool(r["nmpr_gt_mpr"]) else ""
        a(
            f"| {r['malignant_label']} | {r['score_label']} | {fmt_num(r['mean_NMPR'])} | "
            f"{fmt_num(r['mean_MPR'])} | {star}{fmt_num(r['delta_NMPR_minus_MPR'])}{star} | "
            f"{fmt_num(r['U'], 1)} | {fmt_p(r['p'])} | {fmt_p(r['p_onesided_NMPR_gt_MPR'])} | "
            f"{int(r['n_NMPR'])} vs {int(r['n_MPR'])} |"
        )
    a("")
    a("## Per-sample malignant n and TACSTD2 (post-treatment)")
    a("")
    mal_cols = [c for c in post.columns if c.endswith("_n") and not c.startswith("n_")]
    a("Malignant cell counts:")
    a("")
    keep = ["patient", "response"] + mal_cols
    a("| " + " | ".join(keep) + " |")
    a("|" + "|".join(["---"] * len(keep)) + "|")
    for _, r in post.sort_values("patient").iterrows():
        a("| " + " | ".join(str(r[c]) for c in keep) + " |")
    a("")
    a("Immune fractions:")
    a("")
    imm = ["patient", "response", "frac_T_NK", "frac_CD8_only", "frac_CD8_NK", "frac_CXCL13_pos", "frac_cyto_high_T"]
    a("| " + " | ".join(imm) + " |")
    a("|" + "|".join(["---"] * len(imm)) + "|")
    for _, r in post.sort_values("patient").iterrows():
        vals = [str(r["patient"]), str(r["response"])] + [fmt_num(r[c], 3) for c in imm[2:]]
        a("| " + " | ".join(vals) + " |")
    a("")
    a("## Notes")
    a("")
    a("- Author CopyKAT IDs remain unpublished. CopyKAT-like and inferCNV-like rows are window-smoothed expression CNV vs stromal cells on the public UMI.")
    a("- DRMref public labels are third-party marker annotations (30,877 cells), not Hu CopyKAT. DRMref-like is a 12-type marker rebuild.")
    a("- n=12 is small. A true ρ = −0.45 has two-sided Spearman p ≈ 0.14 at n=12. Every p is reported as computed.")
    a("- Highlight threshold ρ ≤ −0.35 is the requested flag for this extra figure, not a new claim cutoff.")
    a("")
    (HERE / "REPORT.md").write_text("\n".join(lines) + "\n")
    print("wrote", HERE / "REPORT.md", "lines", len(lines))


if __name__ == "__main__":
    main()
