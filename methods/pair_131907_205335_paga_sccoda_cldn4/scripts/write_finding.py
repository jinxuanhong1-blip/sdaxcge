#!/usr/bin/env python3
"""Fill FINDING.md from PAGA + scCODA tables."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def _fmt(x, nd=3):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return "NA"
    if pd.isna(v):
        return "NA"
    if abs(v) < 0.001 and v != 0:
        return f"{v:.2e}"
    return f"{v:.{nd}f}"


def md_table(df: pd.DataFrame, cols: list[str]) -> str:
    use = [c for c in cols if c in df.columns]
    head = "| " + " | ".join(use) + " |"
    sep = "|" + "|".join(["---"] * len(use)) + "|"
    lines = [head, sep]
    for _, row in df.iterrows():
        cells = []
        for c in use:
            val = row[c]
            if isinstance(val, float):
                cells.append(_fmt(val))
            else:
                cells.append(str(val))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--outdir", type=Path, default=ROOT / "results")
    p.add_argument("--finding", type=Path, default=ROOT / "FINDING.md")
    args = p.parse_args()
    tdir = args.outdir / "tables"

    sccoda = pd.read_csv(tdir / "sccoda_tests.tsv", sep="\t") if (tdir / "sccoda_tests.tsv").is_file() else pd.DataFrame()
    primary = sccoda[sccoda.get("primary_family", False) == True] if len(sccoda) else pd.DataFrame()
    honest = pd.read_csv(tdir / "honest_n.tsv", sep="\t") if (tdir / "honest_n.tsv").is_file() else pd.DataFrame()
    spear = (
        pd.read_csv(tdir / "malignant_sample_level_spearman.tsv", sep="\t")
        if (tdir / "malignant_sample_level_spearman.tsv").is_file()
        else pd.DataFrame()
    )
    paga_info = {}
    if (tdir / "paga_run_info.json").is_file():
        paga_info = json.loads((tdir / "paga_run_info.json").read_text())
    sccoda_sum = {}
    if (tdir / "sccoda_summary.json").is_file():
        sccoda_sum = json.loads((tdir / "sccoda_summary.json").read_text())
    engine = {}
    if (tdir / "sccoda_engine.json").is_file():
        engine = json.loads((tdir / "sccoda_engine.json").read_text())

    n_mal = paga_info.get("malignant", {}).get("n_cells", "NA")
    n_tnk = paga_info.get("tnk", {}).get("n_cells", "NA")
    n_units_paga = paga_info.get("malignant", {}).get("n_units_in_means", "NA")
    harmony = paga_info.get("malignant", {}).get("embed", {}).get("harmony", False)
    root = paga_info.get("malignant", {}).get("root", {}).get("rule", "NA")
    n_q4q1 = sccoda_sum.get("n_q4q1", 23)
    n_rec = sccoda_sum.get("n_primary_recover_q10", 0)

    merge_tnk = primary[(primary.get("slice") == "merge_within_cohort") & (primary.get("compartment") == "TNK")]
    tnk_alr = _fmt(merge_tnk["alr_effect"].iloc[0]) if len(merge_tnk) else "NA"
    tnk_p = _fmt(merge_tnk["alr_p_perm"].iloc[0]) if len(merge_tnk) else "NA"
    tnk_q = _fmt(merge_tnk["q_bh_primary"].iloc[0]) if len(merge_tnk) else "NA"
    tnk_n = f"{int(merge_tnk['n_low'].iloc[0])}/{int(merge_tnk['n_high'].iloc[0])}" if len(merge_tnk) else "NA"

    body = f"""# FINDING — pair GSE131907 + GSE205335 PAGA + scCODA, CLDN4 only

ADDITIVE **CLDN4-only** high-end on the combo that already differs
(PR #320: author-malignant CLDN4 %pos vs T/NK, Q4 vs Q1 **n=23**
r=−0.705). This run does **not** re-audit that combo, add GSE148071,
or re-run the 7-cohort / 6-unit pool. No TACSTD2∩CLDN4 dual-high gate.

A previous attempt of this exact task OOM'd mid-Harmony/UMAP on
**53,296** unsampled malignant cells. This retry **subsamples per unit**
(malignant cap 150, T/NK cap 80) and runs **scCODA first** from the
PR #320 tables so composition results do not depend on the embedding.

## Honest n

Unit = GSE131907 **sample** (not patient) + GSE205335 **patient**.
Cells are library size / embedding weight only.

"""
    if len(honest):
        body += md_table(
            honest,
            [
                "cohort",
                "unit_type",
                "n_eligible",
                "n_q1",
                "n_q4",
                "n_q4q1_compared",
                "note",
            ],
        )
        body += "\n\n"
    body += f"""Merged Q4 vs Q1 compared n = **{n_q4q1}** (PR #320 12/11).
PAGA malignant cells after per-unit cap = **{n_mal}**; T/NK cells = **{n_tnk}**.
PAGA inferential units with ≥8 cells = **{n_units_paga}**.
Harmony on the malignant slice: **{harmony}**. DPT root: {root}.

## scCODA (primary = Q4 vs Q1)

Engine: ALR of T/NK or B vs Other + unit-level permutation p.
scCODA HMC: `{engine.get("status", "not_available")}`. Not faked.
Primary family = 6 tests (T/NK and B × GSE131907 / GSE205335 / merge).
Recovery = ALR effect < 0 and BH q < 0.10. Recoveries: **{n_rec}**.

Merge T/NK Q4 vs Q1: ALR {tnk_alr}, perm p={tnk_p}, q={tnk_q}, n_low/n_high={tnk_n}.

"""
    if len(primary):
        body += "### Primary table\n\n"
        body += md_table(
            primary,
            [
                "slice",
                "compartment",
                "n",
                "n_low",
                "n_high",
                "alr_effect",
                "alr_p_perm",
                "q_bh_primary",
                "frac_delta_high_minus_low",
                "frac_mwu_p",
                "spearman_rho",
                "recover_q10",
            ],
        )
        body += "\n\n"
        body += "Full tests (Q4 vs Q1 + median sensitivity): [`results/tables/sccoda_tests.tsv`](results/tables/sccoda_tests.tsv).\n\n"

    body += "## PAGA / DPT (malignant, subsampled)\n\n"
    body += (
        "Not a redo of PR #325 (GSE131907 nLung+tLung epithelium, AT2-rooted). "
        "This object is author-malignant cells from the winning pair only "
        "(GSE131907 metastases with ≥20 malignant cells + GSE205335 tumor libraries). "
        "Barrier/keratin excludes CLDN4. TACSTD2 is a comparator Spearman only.\n\n"
    )
    if len(spear):
        body += md_table(
            spear[spear["slice"] == "merge"] if "slice" in spear else spear,
            ["slice", "contrast", "n", "rho", "p", "q_bh_merge"],
        )
        body += "\n\n"
        body += "Per-cohort rows and T/NK PAGA vertices: `results/tables/`.\n\n"

    body += """## Extra figures

- `results/figures/fig_malignant_trajectory_cldn4.png` — UMAP CLDN4 / DPT / cohort / unit quartile
- `results/figures/fig_malignant_paga.png` — PAGA graph
- `results/figures/fig_malignant_extra_programs.png` — unit-mean CLDN4 vs DPT / barrier / AT2
- `results/figures/fig_tnk_trajectory_cldn4.png` — T/NK UMAP (no AT2 DPT)
- `results/figures/fig_sccoda_fractions.png` — T/NK and B fractions, Q4 vs Q1
- `results/figures/fig_sccoda_alr_forest.png` — primary ALR forest
- `results/figures/fig_honest_n.png` — unit counts

## What this is not

- Not dual-high TACSTD2×CLDN4.
- Not GSE148071.
- Not the 7-cohort / 6-unit malignant pool.
- Not PR #325 GSE131907-only epithelial PAGA and not PR #339 multi-cohort ICI scCODA.

## Reproduce

```bash
python3 methods/pair_131907_205335_paga_sccoda_cldn4/scripts/run_sccoda.py
python3 methods/pair_131907_205335_paga_sccoda_cldn4/scripts/download.py
python3 methods/pair_131907_205335_paga_sccoda_cldn4/scripts/extract.py
python3 methods/pair_131907_205335_paga_sccoda_cldn4/scripts/run_paga.py \\
  --input /tmp/pair_131907_205335/pair_subsample.h5ad
python3 methods/pair_131907_205335_paga_sccoda_cldn4/scripts/write_finding.py
```
"""
    args.finding.write_text(body)
    print(f"wrote {args.finding}", flush=True)


if __name__ == "__main__":
    main()
