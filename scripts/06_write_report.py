#!/usr/bin/env python3
"""Fill results/hunt_tf/REPORT.md from the scored tables. No new statistics."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

HEADER = """# Hunt: ELF3 / GRHL1 / KLF4 / TFAP2A with TACSTD2 and CLDN4, NKX2-1 opposite

Honest public-data search. Co-expression is not evidence of direct transcriptional regulation.

## Claim being tested

In public **human and mouse lung RNA**, the user TF set **ELF3, GRHL1, KLF4, TFAP2A** co-correlates with **TACSTD2** and **CLDN4**, and **NKX2-1** anti-correlates.

## What would count as support (fixed before opening expression matrices)

Pair (raw Spearman on log2(CPM+1), n ≥ 12 after filters):

- positive pair: ρ ≥ 0.30 and two-sided p < 0.05
- NKX2-1 pair: ρ ≤ −0.20 and two-sided p < 0.05

Cohort call:

- **SUPPORTED** — ≥ 6 of 8 positive pairs **and** both NKX2-1 pairs anti
- **PARTIAL** — ≥ 4 of 8 positive pairs (NKX2-1 not required)
- **NOT_SUPPORTED** — otherwise
- **UNINFORMATIVE** — targets/TFs mostly undetected, or n < 12

The same rules are applied after rank-partialling **EPCAM** (epithelial content) and, separately, **SFTPC** (AT2 / alveolar content). A raw SUPPORTED call that becomes NOT_SUPPORTED after EPCAM residualization is **composition-confounded**, not a validated network.

Score (continuous, used for ranking only):
`mean ρ(4 TFs × 2 targets) − mean ρ(NKX2-1 × 2 targets)`.

A 300-shuffle permutation of TF sample labels gives an empirical p for that score. It is a within-cohort sanity check, not a study-wide FDR.

## Search space

Source: [recount3](https://rna.recount.bio/) open data (human Gencode v26 / G026, mouse Gencode vM23 / M023). Same pipeline for both species.

1. All SRA studies with ≥ 6 runs in the recount3 project index (6,831 human, 8,195 mouse).
2. Per-study SRA metadata text is searched for lung / airway / pulmonary disease terms.
3. Runs are labelled tissue / culture-organoid / cell line / sorted-or-fluid. Single-cell studies are flagged and not downloaded.
4. Download list (pre-specified, not chosen by looking at correlations):
   - **Anchors (always):** GTEx lung, TCGA LUAD, TCGA LUSC
   - **SRA tissue:** ≥ 12 tissue lung RNA-seq runs, single-cell fraction < 0.25, median spots ≥ 5×10^5, and either ≥ 40% of runs are lung-tagged or the lung tag is run-level
   - **SRA culture / cell line:** culture/organoid with the same n/depth filters; cell-line panels (n ≥ 24 or a multi-line title). Single-line nanoparticle time courses are mostly excluded.
   - Sorted cells, blood, BAL, platelets, and isolated immune studies are excluded

This is a screen of public lung RNA that recount3 already processed. It is not every GEO series, and it is not single-cell ATAC or spatial data.

## Why bulk lung can fake this network

TACSTD2 and CLDN4 are high in airway / regenerating epithelium. NKX2-1 and SFTPC are high in AT2 / alveolar lineage. ELF3, GRHL1, KLF4 and TFAP2A are epithelial transcription factors with airway-leaning published roles. In mixed lung tissue the same sample-to-sample swing in **airway vs alveolar fraction** will raise ELF3/GRHL1/KLF4/TFAP2A/TACSTD2/CLDN4 together and lower NKX2-1, without any of those TFs writing TACSTD2 or CLDN4.

That is why EPCAM- and SFTPC-partial correlations, and the culture/cell-line slice, are required before calling the network supported.

"""

FOOTER = """
## How to rerun

```bash
# metadata cache (once): scripts/00_fetch_metadata.sh
python3 scripts/01_select_lung_studies.py
python3 scripts/02_choose_cohorts.py
python3 scripts/03_fetch_matrices.py
python3 scripts/04_analyze.py
python3 scripts/05_figures.py
python3 scripts/06_write_report.py
```

Gene IDs are in `scripts/genes.json`.
"""


def fmt_rho(x) -> str:
    if pd.isna(x):
        return "NA"
    return f"{x:.2f}"


def call_counts(df: pd.DataFrame, col: str) -> str:
    vc = df[col].value_counts()
    return ", ".join(f"{k} {int(v)}" for k, v in vc.items())


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scores", default="results/hunt_tf/tables/cohort_scores.tsv")
    ap.add_argument("--pairs", default="results/hunt_tf/tables/pair_correlations.tsv")
    ap.add_argument("--screen-summary", default="results/hunt_tf/tables/studies_screened_summary.json")
    ap.add_argument("--out", default="results/hunt_tf/REPORT.md")
    args = ap.parse_args()

    scores = pd.read_csv(args.scores, sep="\t")
    pairs = pd.read_csv(args.pairs, sep="\t")
    informative = scores[scores["call_raw"] != "UNINFORMATIVE"].copy()

    lines = [HEADER, "## Results\n"]
    lines.append(
        f"Scored **{len(scores)}** cohort slices "
        f"({int((scores['organism']=='human').sum())} human, "
        f"{int((scores['organism']=='mouse').sum())} mouse). "
        f"**{len(informative)}** were informative under the pre-specified filters.\n"
    )
    lines.append(f"- Raw calls: {call_counts(informative, 'call_raw') if len(informative) else 'none'}")
    lines.append(f"- EPCAM-residual calls: {call_counts(informative, 'call_epcam') if len(informative) else 'none'}")
    lines.append(f"- SFTPC-residual calls: {call_counts(informative, 'call_sftpc') if len(informative) else 'none'}\n")

    raw_sup = informative[informative["call_raw"] == "SUPPORTED"]
    robust = informative[(informative["call_raw"] == "SUPPORTED") & (informative["call_epcam"] == "SUPPORTED")]
    confounded = informative[(informative["call_raw"] == "SUPPORTED") & (informative["call_epcam"] != "SUPPORTED")]

    lines.append("### Verdict\n")
    lines.append(
        "**The pre-specified 4-TF + NKX2-1-anti package is not supported** in public human or mouse lung RNA.\n"
    )
    lines.append(
        f"{len(raw_sup)} of {len(informative)} informative slices meet the raw SUPPORTED rule; "
        f"{len(robust)} still meets it after EPCAM residualization. "
        "That is not a validation. The two largest clean human matrices fail the full rule, "
        "mouse has zero SUPPORTED calls, and the few raw hits are small epithelial-cell / EMT / single-lab tumor series, "
        "not independent bulk parenchyma.\n"
    )
    lines.append("What actually holds, and what does not:\n")
    lines.append(
        "- **ELF3** is the most consistent positive correlate of TACSTD2 and especially CLDN4 "
        "(GTEx, TCGA LUAD/LUSC, CCLE, LUAD line catalogue)."
    )
    lines.append(
        "- **GRHL1** tracks both targets in TCGA tumors and cell-line panels, but is weakly *negative* in GTEx lung. "
        "It is not a universal lung co-correlate."
    )
    lines.append("- **KLF4** is mixed: present in GTEx and some line panels, absent or reversed in TCGA LUAD vs CLDN4.")
    lines.append("- **TFAP2A** is not part of this module. It is near-zero in GTEx and often fails in tumors and CCLE.")
    lines.append(
        "- **NKX2-1 anti-correlation is rejected.** In GTEx lung NKX2-1 is *positively* correlated with TACSTD2 (ρ=0.56) "
        "and CLDN4 (ρ=0.71). In a 160-line lung cancer panel it tracks the same genes (ρ≈0.90). "
        "TCGA tumors are near zero / slightly positive. Apparent anti-correlation shows up in infection, EMT, "
        "and injury time courses where AT2 / NKX2-1 programs drop — that is biology of damage, not a TF network."
    )
    lines.append(
        "- A smaller claim — ELF3 ± GRHL1 ± KLF4 co-vary with TACSTD2/CLDN4 as an epithelial/barrier state — "
        "is visible in tumors and cell-line panels. That is still co-expression, not regulation.\n"
    )

    # Anchors
    lines.append("### Anchor cohorts (always downloaded)\n")
    lines.append("| cohort | n | raw call | EPCAM call | mean ρ TFs | mean ρ NKX2-1 | score | p_perm |")
    lines.append("|---|---:|---|---|---:|---:|---:|---:|")
    for cid in ["GTEx_LUNG", "TCGA_LUAD_tumor", "TCGA_LUAD_normal", "TCGA_LUSC_tumor", "TCGA_LUSC_normal"]:
        hit = scores[scores["cohort_id"] == cid]
        if hit.empty:
            continue
        r = hit.iloc[0]
        lines.append(
            f"| {cid} | {int(r.n_used)} | {r.call_raw} | {r.call_epcam} | "
            f"{fmt_rho(r.mean_rho_pos)} | {fmt_rho(r.mean_rho_nkx)} | {fmt_rho(r.score)} | {fmt_rho(r.p_perm)} |"
        )
    lines.append("")

    def pair_table(cohort_id: str) -> None:
        sub = pairs[(pairs["cohort_id"] == cohort_id) & (pairs["scale"] == "raw")]
        if sub.empty:
            return
        lines.append(f"#### {cohort_id} pairwise raw Spearman\n")
        lines.append("| TF | TACSTD2 ρ (p) | CLDN4 ρ (p) |")
        lines.append("|---|---|---|")
        for tf in ["ELF3", "GRHL1", "KLF4", "TFAP2A", "NKX2-1"]:
            a = sub[(sub["tf"] == tf) & (sub["target"] == "TACSTD2")]
            b = sub[(sub["tf"] == tf) & (sub["target"] == "CLDN4")]
            def cell(df):
                if df.empty or pd.isna(df.iloc[0]["rho"]):
                    return "NA"
                return f"{df.iloc[0]['rho']:.2f} ({df.iloc[0]['p']:.1e})"
            lines.append(f"| {tf} | {cell(a)} | {cell(b)} |")
        lines.append("")

    for cid in ["GTEx_LUNG", "TCGA_LUAD_tumor", "TCGA_LUSC_tumor", "SRP186687", "DRP001919"]:
        pair_table(cid)

    # Per-TF pass rates on informative tissue/tumor/cell-line
    lines.append("### Per-TF pair pass rates (informative tissue, tumor, cell-line)\n")
    focus = pairs[
        (pairs["scale"] == "raw")
        & (pairs["tf"].isin(["ELF3", "GRHL1", "KLF4", "TFAP2A", "NKX2-1"]))
        & (pairs["target"].isin(["TACSTD2", "CLDN4"]))
        & (pairs["cohort_id"].isin(informative.loc[informative["material"].isin(["tissue", "tumor", "cell_line"]), "cohort_id"]))
    ]
    lines.append("| TF | target | n pairs | median ρ | fraction ρ≥0.30, p<0.05 | fraction ρ≤−0.20, p<0.05 |")
    lines.append("|---|---|---:|---:|---:|---:|")
    for tf in ["ELF3", "GRHL1", "KLF4", "TFAP2A", "NKX2-1"]:
        for tgt in ["TACSTD2", "CLDN4"]:
            sub = focus[(focus["tf"] == tf) & (focus["target"] == tgt)].dropna(subset=["rho"])
            if sub.empty:
                continue
            pos = ((sub["rho"] >= 0.30) & (sub["p"] < 0.05)).mean()
            neg = ((sub["rho"] <= -0.20) & (sub["p"] < 0.05)).mean()
            lines.append(f"| {tf} | {tgt} | {len(sub)} | {sub['rho'].median():.2f} | {pos:.2f} | {neg:.2f} |")
    lines.append("")

    # Cross-species
    lines.append("### Human vs mouse tissue\n")
    for org in ["human", "mouse"]:
        sub = informative[(informative["organism"] == org) & (informative["material"].isin(["tissue", "tumor", "adjacent_normal"]))]
        if sub.empty:
            lines.append(f"- **{org} tissue/tumor:** no informative cohorts")
            continue
        lines.append(
            f"- **{org} tissue/tumor:** n={len(sub)}, raw {call_counts(sub, 'call_raw')}; "
            f"EPCAM {call_counts(sub, 'call_epcam')}; "
            f"median score {sub['score'].median():.2f}"
        )
    lines.append("")

    cul = informative[informative["material"].isin(["culture_or_organoid", "cell_line"])]
    lines.append("### Culture / cell-line slice (less composition)\n")
    if cul.empty:
        lines.append("No informative culture or cell-line cohorts.\n")
    else:
        lines.append(
            f"{len(cul)} informative culture/cell-line slices. "
            f"Raw: {call_counts(cul, 'call_raw')}. EPCAM: {call_counts(cul, 'call_epcam')}.\n"
        )
        show = cul.sort_values("score", ascending=False).head(8)
        lines.append("| cohort | material | n | raw | EPCAM | score | title |")
        lines.append("|---|---|---:|---|---|---:|---|")
        for r in show.itertuples():
            title = str(r.study_title)[:80].replace("|", "/")
            lines.append(
                f"| {r.cohort_id} | {r.material} | {int(r.n_used)} | {r.call_raw} | {r.call_epcam} | {fmt_rho(r.score)} | {title} |"
            )
        lines.append("")

    # Composition diagnostic
    lines.append("### Composition diagnostic\n")
    cols = [
        ("TACSTD2 vs SCGB1A1 (club)", "comp_TACSTD2__SCGB1A1"),
        ("TACSTD2 vs SFTPC (AT2)", "comp_TACSTD2__SFTPC"),
        ("CLDN4 vs SCGB1A1", "comp_CLDN4__SCGB1A1"),
        ("CLDN4 vs SFTPC", "comp_CLDN4__SFTPC"),
        ("NKX2-1 vs SFTPC", "comp_NKX2-1__SFTPC"),
        ("NKX2-1 vs SCGB1A1", "comp_NKX2-1__SCGB1A1"),
    ]
    lines.append("Mean raw Spearman across informative cohorts:\n")
    lines.append("| pair | mean ρ |")
    lines.append("|---|---:|")
    for label, col in cols:
        if col in informative.columns:
            lines.append(f"| {label} | {fmt_rho(informative[col].mean())} |")
    lines.append("")
    lines.append(
        "If TACSTD2/CLDN4 track SCGB1A1 and anti-track SFTPC, while NKX2-1 tracks SFTPC, "
        "the hunt recovered the airway-vs-alveolar axis, not a new TF circuit.\n"
    )

    if len(raw_sup):
        lines.append("### The three raw SUPPORTED hits (do not over-read them)\n")
        lines.append("| cohort | n | EPCAM | note |")
        lines.append("|---|---:|---|---|")
        notes = {
            "SRP066794": "EMT time course in culture (mislabelled tissue). TFAP2A goes the *wrong* way (ρ≈−0.8). ELF3/GRHL1/KLF4 collapse together during EMT.",
            "SRP076732": "Title is human lung epithelial cells, n=18. ELF3 vs TACSTD2 is 0.01. Only double-SUPPORTED hit; too small and not parenchyma.",
            "SRP223534": "29 lung tumor RNAs from one lab, multiple aliquots per patient. Batch/patient structure can inflate ρ.",
        }
        for r in raw_sup.itertuples():
            lines.append(f"| {r.cohort_id} | {int(r.n_used)} | {r.call_epcam} | {notes.get(r.cohort_id, str(r.study_title)[:100])} |")
        lines.append("")

    ep_any = informative[informative["call_epcam"] == "SUPPORTED"]
    if len(ep_any):
        lines.append("### EPCAM-residual SUPPORTED (includes cohorts that were only PARTIAL on raw)\n")
        lines.append("| cohort | raw | n | score | title |")
        lines.append("|---|---|---:|---:|---|")
        for r in ep_any.sort_values("score", ascending=False).itertuples():
            title = str(r.study_title)[:80].replace("|", "/")
            lines.append(f"| {r.cohort_id} | {r.call_raw} | {int(r.n_used)} | {fmt_rho(r.score)} | {title} |")
        lines.append("")
        lines.append(
            "None of these is a large independent bulk-lung parenchyma series. "
            "They do not rescue the 4-TF + NKX2-1 package.\n"
        )

    lines.append("## What this is not\n")
    lines.append("- Not ChIP, CUT&RUN, motif, or reporter evidence.")
    lines.append("- Not single-cell (those studies were excluded; recount3 gene sums of scRNA-seq are the wrong object).")
    lines.append("- Not a claim that GRHL1 or ELF3 is *the* TACSTD2/CLDN4 factor. Sister PR A10 tests GRHL1 in TCGA with a different matrix.")
    lines.append("- recount3 ends in 2020-era SRA. Later GEO lung series are not in this hunt.\n")
    lines.append(FOOTER)

    Path(args.out).write_text("\n".join(lines))
    print("wrote", args.out)


if __name__ == "__main__":
    main()
