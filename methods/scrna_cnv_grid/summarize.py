#!/usr/bin/env python3
"""Combine GSE207422 + GSE241934 grids and write the direction-recovery FINDING."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path("results/scrna_cnv_grid")


def load_primary() -> pd.DataFrame:
    frames = []
    for p in [
        ROOT / "GSE207422" / "direction_table.tsv",
        ROOT / "GSE241934" / "direction_table.tsv",
    ]:
        if p.exists():
            frames.append(pd.read_csv(p, sep="\t"))
    if not frames:
        raise SystemExit("no direction tables yet")
    return pd.concat(frames, ignore_index=True)


def fmt_p(x) -> str:
    if pd.isna(x):
        return "NA"
    x = float(x)
    if x < 0.001:
        return f"{x:.2e}"
    return f"{x:.3f}"


def fmt_rho(x) -> str:
    if pd.isna(x):
        return "NA"
    return f"{float(x):+.3f}"


def fmt_mean(x) -> str:
    if pd.isna(x):
        return "NA"
    return f"{float(x):.3f}"


def main() -> None:
    df = load_primary()
    # Prefer one T/NK def per dataset for the lead table.
    prefer = {
        ("GSE207422", "drmref_TNK"): 0,
        ("GSE207422", "lineage_TNK"): 1,
        ("GSE241934", "author_TNK"): 0,
        ("GSE241934", "lineage_TNK"): 1,
    }
    df["pref"] = [
        prefer.get((r.dataset, r.tnk_def), 2) for r in df.itertuples(index=False)
    ]
    lead = (
        df.sort_values(["dataset", "cohort", "malignant_def", "pref"])
        .groupby(["dataset", "cohort", "malignant_def"], as_index=False)
        .first()
    )
    lead.to_csv(ROOT / "direction_lead.tsv", sep="\t", index=False)
    df.to_csv(ROOT / "grid_stats_primary.tsv", sep="\t", index=False)

    both = lead[lead["recovers_user_direction"] == "both"].copy()
    nmpr_only = lead[lead["recovers_user_direction"] == "NMPR>MPR_only"]
    rho_only = lead[lead["recovers_user_direction"] == "rho_neg_only"]

    lines = []
    lines.append("# Malignant-definition grid — which public calls recover the user direction")
    lines.append("")
    lines.append("User A3 direction (taken as given): **malignant TACSTD2 higher in NMPR than MPR**, and **negative Spearman vs T/NK fraction**. This slice grids public malignant definitions on GSE207422 and GSE241934 and reports honest n / ρ / p. It does not re-state the original claim.")
    lines.append("")
    lines.append("Primary metric: mean `log1p(CP10k)` TACSTD2 in the called-malignant compartment. pCR counted as MPR. Sample is the unit.")
    lines.append("")
    lines.append("## Definitions that recover both directions")
    lines.append("")
    if both.empty:
        lines.append("None of the lead (preferred T/NK) rows recovered both directions. See the per-T/NK table below and `grid_stats_primary.tsv` for every pairing.")
    else:
        lines.append("| Dataset | Cohort | Malignant def | T/NK def | n (NMPR/MPR) | mean NMPR | mean MPR | p (NMPR vs MPR) | ρ | p(ρ) |")
        lines.append("|---|---|---|---|---|---|---|---|---|---|")
        for r in both.itertuples(index=False):
            lines.append(
                f"| {r.dataset} | {r.cohort} | `{r.malignant_def}` | {r.tnk_def} | "
                f"{int(r.n)} ({int(r.n_NMPR)}/{int(r.n_MPR)}) | {fmt_mean(r.mean_NMPR)} | {fmt_mean(r.mean_MPR)} | "
                f"{fmt_p(r.p_nmpr)} | {fmt_rho(r.rho)} | {fmt_p(r.p_rho)} |"
            )
    lines.append("")
    lines.append("## One-half recoveries (lead T/NK)")
    lines.append("")
    lines.append(f"- NMPR>MPR only: {len(nmpr_only)} definition×cohort rows")
    lines.append(f"- ρ<0 only: {len(rho_only)} definition×cohort rows")
    lines.append("")
    if not nmpr_only.empty:
        lines.append("NMPR>MPR only: " + ", ".join(f"{r.dataset}/{r.cohort}/{r.malignant_def} (p={fmt_p(r.p_nmpr)})" for r in nmpr_only.itertuples(index=False)))
        lines.append("")
    if not rho_only.empty:
        lines.append("ρ<0 only: " + ", ".join(f"{r.dataset}/{r.cohort}/{r.malignant_def} (ρ={fmt_rho(r.rho)}, p={fmt_p(r.p_rho)})" for r in rho_only.itertuples(index=False)))
        lines.append("")

    lines.append("## Full lead table (preferred T/NK per dataset)")
    lines.append("")
    lines.append("| Dataset | Cohort | Malignant def | T/NK | n | NMPR | MPR | Δ | p | ρ | p(ρ) | direction |")
    lines.append("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in lead.sort_values(["dataset", "cohort", "malignant_def"]).itertuples(index=False):
        lines.append(
            f"| {r.dataset} | {r.cohort} | `{r.malignant_def}` | {r.tnk_def} | {int(r.n)} | "
            f"{fmt_mean(r.mean_NMPR)} | {fmt_mean(r.mean_MPR)} | {fmt_mean(r.delta_NMPR_minus_MPR)} | "
            f"{fmt_p(r.p_nmpr)} | {fmt_rho(r.rho)} | {fmt_p(r.p_rho)} | {r.recovers_user_direction} |"
        )
    lines.append("")
    lines.append("## Notes")
    lines.append("")
    lines.append("- GSE207422 has **no public author CopyKAT IDs**. `thirdparty_drmref` is Liu *NAR* 2024 DRMref, 30,877 / 92,330 cells.")
    lines.append("- GSE241934 `author_epi` is the deposited `major.cell.type==Epi` compartment (IIT 1,699 cells; paper CopyKAT malignant N=1,669).")
    lines.append("- CNV is a window-smoothed expression score (win=25, stromal reference), not the unpublished R `copykat` object.")
    lines.append("- n=11–12 (IIT / GSE207422) is small; a true ρ=−0.45 has two-sided Spearman p≈0.14 at n=12. Direction and p are reported separately.")
    lines.append("- Every pairing (all T/NK defs, %positive, pseudobulk, min≥5 cells) is in the per-dataset `grid_stats.tsv` files.")
    lines.append("")
    (ROOT / "FINDING.md").write_text("\n".join(lines) + "\n")
    print("wrote", ROOT / "FINDING.md", "both=", len(both), "lead=", len(lead), flush=True)


if __name__ == "__main__":
    main()
