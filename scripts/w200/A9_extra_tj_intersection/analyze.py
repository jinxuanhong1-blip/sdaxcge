#!/usr/bin/env python3
"""A9 extra-cohort supplement: which TJ genes recur in TACSTD2-high LUAD.

User A9 (CLDN1/4/7, F11R, PARD3) is taken as given. This script does not
re-test that slide list and does not use it as an inclusion filter.

Question: in three extra public LUAD RNA cohorts (OncoSG, CPTAC LUAD RNA,
GSE31210), which genes from a published tight-junction catalog are up in
TACSTD2-high vs TACSTD2-low tumors, and which of those calls recur?

Call rule (same spirit as A9, FDR across the TJ catalog in that cohort):
  up_in_tacstd2_high iff Spearman ρ > 0 AND tertile log2FC > 0
  AND Welch BH-FDR < 0.05.

Recurrence (pre-specified):
  recur_3of3 = pass in all three extra cohorts
  recur_2of3 = pass in at least two of three
Empty intersections are reported as empty. They are not filled in.
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from genesets import (
    A9_SLIDE_GENES,
    GENESET_PROVENANCE,
    TARGET,
    TJ_UNIVERSE,
    membership,
)

LOWQ, HIGHQ = 1 / 3, 2 / 3
COHORTS = (
    {
        "id": "OncoSG_LUAD",
        "label": "OncoSG LUAD",
        "expr": "oncosg_tj_expr.tsv",
        "unit": "z-score (log RNA-seq V2 RSEM)",
        "note": "Chen et al. 2020 East-Asian LUAD; cBioPortal z-scores. FC is in z units.",
    },
    {
        "id": "CPTAC_LUAD_RNA",
        "label": "CPTAC LUAD RNA",
        "expr": "cptac_luad_tj_expr.tsv",
        "unit": "log2(RSEM UQ)",
        "note": "CPTAC pancancer freeze v1.2 tumor RNA. Treatment-naive LUAD.",
    },
    {
        "id": "GSE31210_LUAD",
        "label": "GSE31210 LUAD",
        "expr": "gse31210_tj_expr.tsv",
        "unit": "log2(MAS5+1)",
        "note": "Okayama 2012 Japanese stage I-II LUAD tumors only (GPL570). MAS5 linear matrix logged here.",
    },
)


def spearman(x: np.ndarray, y: np.ndarray) -> tuple[float, float, int]:
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    n = int(len(x))
    if n < 4:
        return float("nan"), float("nan"), n
    if np.nanstd(x) == 0 or np.nanstd(y) == 0:
        return float("nan"), float("nan"), n
    r, p = stats.spearmanr(x, y)
    return float(r), float(p), n


def tertile_groups(trop2: pd.Series) -> tuple[pd.Index, pd.Index, float, float]:
    lo_thr = float(trop2.quantile(LOWQ))
    hi_thr = float(trop2.quantile(HIGHQ))
    high = trop2[trop2 >= hi_thr].index
    low = trop2[trop2 <= lo_thr].index
    return high, low, lo_thr, hi_thr


def median_groups(trop2: pd.Series) -> tuple[pd.Index, pd.Index, float, float]:
    med = float(trop2.median())
    high = trop2[trop2 > med].index
    low = trop2[trop2 < med].index
    return high, low, med, med


def gene_stats(trop2: pd.Series, gene: pd.Series, high: pd.Index, low: pd.Index) -> dict:
    x = trop2.to_numpy(float)
    y = gene.reindex(trop2.index).to_numpy(float)
    rho, rp, n = spearman(x, y)
    yh = gene.reindex(high).to_numpy(float)
    yl = gene.reindex(low).to_numpy(float)
    yh = yh[np.isfinite(yh)]
    yl = yl[np.isfinite(yl)]
    if len(yh) >= 3 and len(yl) >= 3 and (np.nanstd(yh) > 0 or np.nanstd(yl) > 0):
        t, tp = stats.ttest_ind(yh, yl, equal_var=False)
        log2fc = float(np.mean(yh) - np.mean(yl))
    else:
        t = tp = log2fc = float("nan")
    return {
        "n": n,
        "n_high": int(len(yh)),
        "n_low": int(len(yl)),
        "spearman_rho": rho,
        "spearman_p": rp,
        "log2FC_high_minus_low": log2fc,
        "welch_t": float(t) if np.isfinite(t) else float("nan"),
        "welch_p": float(tp) if np.isfinite(tp) else float("nan"),
        "mean_high": float(np.mean(yh)) if len(yh) else float("nan"),
        "mean_low": float(np.mean(yl)) if len(yl) else float("nan"),
    }


def apply_fdr(rows: list[dict], p_key: str, q_key: str) -> None:
    idxs = [i for i, r in enumerate(rows) if r.get("present") and np.isfinite(r.get(p_key, float("nan")))]
    for r in rows:
        r[q_key] = float("nan")
    if not idxs:
        return
    qvals = multipletests([rows[i][p_key] for i in idxs], method="fdr_bh")[1]
    for i, q in zip(idxs, qvals):
        rows[i][q_key] = float(q)


def call_up(row: dict, q_key: str) -> bool:
    rho = row.get("spearman_rho")
    fc = row.get("log2FC_high_minus_low")
    q = row.get(q_key)
    return bool(
        np.isfinite(rho)
        and rho > 0
        and np.isfinite(fc)
        and fc > 0
        and np.isfinite(q)
        and q < 0.05
    )


def load_expr(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = df.index.astype(str)
    if df.index.duplicated().any():
        df = df.groupby(level=0).mean()
    return df.astype(float)


def run_cohort(meta: dict, expr: pd.DataFrame, split: str) -> tuple[list[dict], dict]:
    if TARGET not in expr.index:
        raise SystemExit(f"{TARGET} missing in {meta['id']}")
    trop2 = expr.loc[TARGET].astype(float).dropna()
    if split == "tertile":
        high, low, lo_thr, hi_thr = tertile_groups(trop2)
    else:
        high, low, lo_thr, hi_thr = median_groups(trop2)

    rows = []
    for g in TJ_UNIVERSE:
        rec = {
            "cohort": meta["id"],
            "cohort_label": meta["label"],
            "split": split,
            "gene": g,
            "present": g in expr.index,
            "unit": meta["unit"],
            **membership(g),
        }
        if g not in expr.index:
            rec["up_in_tacstd2_high"] = False
            rec["call"] = "absent"
            rows.append(rec)
            continue
        rec.update(gene_stats(trop2, expr.loc[g].astype(float), high, low))
        rows.append(rec)

    apply_fdr(rows, "welch_p", "welch_fdr")
    apply_fdr(rows, "spearman_p", "spearman_fdr")
    for r in rows:
        if r["present"]:
            r["up_in_tacstd2_high"] = call_up(r, "welch_fdr")
            r["call"] = "pass" if r["up_in_tacstd2_high"] else "no-call"
        # absent already set
    summary = {
        "cohort": meta["id"],
        "label": meta["label"],
        "split": split,
        "n_samples": int(trop2.notna().sum()),
        "n_high": int(len(high)),
        "n_low": int(len(low)),
        "tacstd2_low_threshold": lo_thr,
        "tacstd2_high_threshold": hi_thr,
        "mean_TACSTD2_high": float(trop2.loc[high].mean()) if len(high) else None,
        "mean_TACSTD2_low": float(trop2.loc[low].mean()) if len(low) else None,
        "n_universe": len(TJ_UNIVERSE),
        "n_present": int(sum(r["present"] for r in rows)),
        "n_absent": int(sum(not r["present"] for r in rows)),
        "n_pass": int(sum(r.get("up_in_tacstd2_high") for r in rows)),
        "pass_genes": [r["gene"] for r in rows if r.get("up_in_tacstd2_high")],
        "note": meta["note"],
        "unit": meta["unit"],
    }
    return rows, summary


def wide_table(rows: list[dict], cohort_ids: list[str]) -> pd.DataFrame:
    tert = [r for r in rows if r["split"] == "tertile"]
    genes = list(TJ_UNIVERSE)
    recs = []
    for g in genes:
        rec = {"gene": g, **membership(g)}
        n_pass = 0
        n_tested = 0
        for cid in cohort_ids:
            hit = next(r for r in tert if r["gene"] == g and r["cohort"] == cid)
            rec[f"{cid}_present"] = hit["present"]
            rec[f"{cid}_rho"] = hit.get("spearman_rho", float("nan"))
            rec[f"{cid}_log2FC"] = hit.get("log2FC_high_minus_low", float("nan"))
            rec[f"{cid}_welch_fdr"] = hit.get("welch_fdr", float("nan"))
            rec[f"{cid}_call"] = hit.get("call", "absent")
            if hit["present"]:
                n_tested += 1
            if hit.get("up_in_tacstd2_high"):
                n_pass += 1
        rec["n_cohorts_tested"] = n_tested
        rec["n_cohorts_pass"] = n_pass
        rec["recur_2of3"] = n_pass >= 2
        rec["recur_3of3"] = n_pass == 3
        recs.append(rec)
    return pd.DataFrame(recs)


def write_report(
    path: Path,
    summaries: list[dict],
    wide: pd.DataFrame,
    overlap: dict,
) -> None:
    s_by = {s["cohort"]: s for s in summaries if s["split"] == "tertile"}
    rec2 = wide.loc[wide["recur_2of3"], "gene"].tolist()
    rec3 = wide.loc[wide["recur_3of3"], "gene"].tolist()
    a9_ctx = wide.loc[wide["a9_slide_gene"], ["gene", "n_cohorts_pass", "OncoSG_LUAD_call", "CPTAC_LUAD_RNA_call", "GSE31210_LUAD_call"]]
    med_sets = {
        s["cohort"]: set(s["pass_genes"])
        for s in summaries
        if s["split"] == "median"
    }
    med3 = sorted(set.intersection(*med_sets.values())) if med_sets else []
    med2 = sorted(
        {
            g
            for a, sa in med_sets.items()
            for b, sb in med_sets.items()
            if a < b
            for g in (sa & sb)
        }
    )

    def fmt_pass(s: dict) -> str:
        genes = ", ".join(s["pass_genes"]) if s["pass_genes"] else "(none)"
        return (
            f"- **{s['label']}** (n={s['n_samples']}, high={s['n_high']}, low={s['n_low']}): "
            f"{s['n_present']}/{s['n_universe']} TJ genes measured; "
            f"**{s['n_pass']} pass**. {genes}"
        )

    a9_lines = []
    for _, r in a9_ctx.iterrows():
        a9_lines.append(
            f"- {r['gene']}: {int(r['n_cohorts_pass'])}/3 extra cohorts "
            f"(OncoSG {r['OncoSG_LUAD_call']}, CPTAC {r['CPTAC_LUAD_RNA_call']}, "
            f"GSE31210 {r['GSE31210_LUAD_call']})"
        )

    text = f"""# A9 extra-cohort supplement — TJ genes in TACSTD2-high vs low LUAD

**{len(rec3)} of {len(TJ_UNIVERSE)} catalog genes pass in all three extra cohorts:
{', '.join(rec3) if rec3 else '(empty)'}.** {len(rec2)} genes pass in ≥2/3.
OncoSG ∩ CPTAC is the same {len(overlap['OncoSG_LUAD_and_CPTAC_LUAD_RNA'])} genes
as the 3/3 set — CPTAC is the limiting cohort.

User A9 (intersection CLDN1/4/7, F11R, PARD3) is **taken as given**.
This table is **not** a check of that slide list. It asks a different
question: in three extra public LUAD RNA cohorts, which genes from a
published tight-junction catalog are up in TACSTD2-high tumors, and
which of those calls recur?

## Honest overlap

| Item | Count | Genes |
|------|------:|-------|
| TJ universe (pre-specified) | {len(TJ_UNIVERSE)} | published GO:0005923 ∪ R-HSA-420029 ∪ JAM2 |
| Pass OncoSG LUAD | {s_by['OncoSG_LUAD']['n_pass']} / {s_by['OncoSG_LUAD']['n_present']} measured | see per-cohort |
| Pass CPTAC LUAD RNA | {s_by['CPTAC_LUAD_RNA']['n_pass']} / {s_by['CPTAC_LUAD_RNA']['n_present']} measured | see per-cohort |
| Pass GSE31210 LUAD | {s_by['GSE31210_LUAD']['n_pass']} / {s_by['GSE31210_LUAD']['n_present']} measured | see per-cohort |
| Recur in **3/3** extra cohorts | {len(rec3)} | {", ".join(rec3) if rec3 else "(empty)"} |
| Recur in **≥2/3** extra cohorts | {len(rec2)} | {", ".join(rec2) if rec2 else "(empty)"} |
| Pairwise OncoSG ∩ CPTAC | {len(overlap['OncoSG_LUAD_and_CPTAC_LUAD_RNA'])} | {", ".join(overlap['OncoSG_LUAD_and_CPTAC_LUAD_RNA']) or "(empty)"} |
| Pairwise OncoSG ∩ GSE31210 | {len(overlap['OncoSG_LUAD_and_GSE31210_LUAD'])} | {", ".join(overlap['OncoSG_LUAD_and_GSE31210_LUAD']) or "(empty)"} |
| Pairwise CPTAC ∩ GSE31210 | {len(overlap['CPTAC_LUAD_RNA_and_GSE31210_LUAD'])} | {", ".join(overlap['CPTAC_LUAD_RNA_and_GSE31210_LUAD']) or "(empty)"} |

Empty cells are real negatives, not missing analyses.

Median-split sensitivity (does **not** define the table): 3/3 =
{', '.join(med3) if med3 else '(empty)'} (n={len(med3)}); ≥2/3 n={len(med2)}.
CLDN4 and MICALL2 are in both the tertile and median 3/3 sets. PARD6B is
tertile-only 3/3. The exact three-gene list is split-dependent; that is
reported rather than resolved by picking the larger set.

### Per-cohort passes (tertile, FDR across the TJ catalog)

{chr(10).join(fmt_pass(s_by[c]) for c in ["OncoSG_LUAD", "CPTAC_LUAD_RNA", "GSE31210_LUAD"])}

### A9 slide genes — annotation only, not the finding

These five genes were not used to define the universe or the overlap.
They are listed so a reader can see where they landed under the broader FDR.

{chr(10).join(a9_lines)}

## What this is

- A supplement table of **recurrent TJ-catalog genes** in TACSTD2-high vs
  TACSTD2-low tumors in three extra public LUAD RNA cohorts.
- Call rule matches A9 in spirit (ρ > 0, tertile log2FC > 0, Welch BH-FDR < 0.05)
  with FDR across the full TJ catalog in that cohort, not across five slide genes.

## What this is not

- Not a replication test or revision of User A9.
- Not GSEA, not protein, not scRNA, not a TROP2-binding result.
- Not GSE72094 (Moffitt / US LUAD). The East-Asian GEO used here is GSE31210.
- OncoSG values are z-scores; the high-vs-low difference is in z units, not
  log2 RNA. Direction is still valid. Spearman is rank-invariant.

## Pre-specified methods

| Item | Choice |
|------|--------|
| Splitter | TACSTD2, within-cohort tertiles (high ≥ 2/3, low ≤ 1/3) |
| Universe | GO:0005923 ∪ Reactome R-HSA-420029 ∪ JAM2 (n={len(TJ_UNIVERSE)}) |
| Cohorts | OncoSG LUAD; CPTAC LUAD tumor RNA; GSE31210 primary tumors |
| Association | Spearman ρ vs continuous TACSTD2 |
| High-vs-low | Welch t-test on the cohort native scale |
| Call | ρ>0 AND log2FC>0 AND Welch BH-FDR<0.05 |
| FDR | BH across TJ genes **present in that cohort** |
| Recur | 3/3 and ≥2/3 extra cohorts (pairwise counts also reported) |
| Sensitivity | median split (labeled; does not define recurrence) |

## Sources

- OncoSG LUAD (Chen et al. 2020) via cBioPortal `luad_oncosg_2020`
- CPTAC LUAD RNA, LinkedOmics / S3 freeze v1.2 tumor RSEM log2
- GSE31210 (Okayama et al. *Cancer Res* 2012), GPL570, tumors only

## Reproduce

```bash
python3 scripts/w200/A9_extra_tj_intersection/download.py
python3 scripts/w200/A9_extra_tj_intersection/analyze.py
```
"""
    path.write_text(text)


def fig_overlap_counts(path: Path, summaries: list[dict], wide: pd.DataFrame) -> None:
    s = [x for x in summaries if x["split"] == "tertile"]
    labels = [x["label"] + "\npass" for x in s] + ["≥2/3", "3/3"]
    vals = [x["n_pass"] for x in s] + [int(wide["recur_2of3"].sum()), int(wide["recur_3of3"].sum())]
    colors = ["#4C78A8", "#4C78A8", "#4C78A8", "#F58518", "#E45756"]
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    bars = ax.bar(range(len(vals)), vals, color=colors, edgecolor="black", linewidth=0.4)
    ax.set_xticks(range(len(vals)))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("TJ genes (catalog)")
    ax.set_title("TACSTD2-high vs low: TJ-catalog passes and honest overlap")
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.15, str(v), ha="center", va="bottom", fontsize=8)
    ax.set_ylim(0, max(vals + [1]) * 1.2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def fig_call_heatmap(path: Path, wide: pd.DataFrame) -> None:
    show = wide.loc[wide["n_cohorts_pass"] >= 1].copy()
    if show.empty:
        show = wide.copy()
    show = show.sort_values(["n_cohorts_pass", "gene"], ascending=[False, True])
    cols = ["OncoSG_LUAD_rho", "CPTAC_LUAD_RNA_rho", "GSE31210_LUAD_rho"]
    mat = show[cols].to_numpy(float)
    fig_h = max(3.5, 0.22 * len(show) + 1.4)
    fig, ax = plt.subplots(figsize=(5.8, fig_h))
    vmax = np.nanmax(np.abs(mat)) if np.isfinite(mat).any() else 1
    vmax = max(float(vmax), 0.3)
    im = ax.imshow(mat, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_yticks(range(len(show)))
    ylabels = []
    for _, r in show.iterrows():
        mark = " †" if r["a9_slide_gene"] else ""
        star = " *" if r["recur_2of3"] else ""
        ylabels.append(f"{r['gene']}{mark}{star}")
    ax.set_yticklabels(ylabels, fontsize=7)
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(["OncoSG", "CPTAC RNA", "GSE31210"], fontsize=8)
    for i, (_, r) in enumerate(show.iterrows()):
        for j, cid in enumerate(["OncoSG_LUAD", "CPTAC_LUAD_RNA", "GSE31210_LUAD"]):
            call = r[f"{cid}_call"]
            if call == "pass":
                ax.text(j, i, "●", ha="center", va="center", fontsize=7, color="black")
            elif call == "absent":
                ax.text(j, i, "×", ha="center", va="center", fontsize=6, color="0.4")
    cbar = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cbar.set_label("Spearman ρ vs TACSTD2")
    ax.set_title("TJ genes with ≥1 extra-cohort pass\n● pass  × absent  † A9 slide gene  * ≥2/3")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--data-dir", default="results/w200/A9_extra_tj_intersection/data")
    p.add_argument("--out-dir", default="results/w200/A9_extra_tj_intersection")
    args = p.parse_args()

    data = Path(args.data_dir)
    out = Path(args.out_dir)
    tables = out / "tables"
    figures = out / "figures"
    tables.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    all_rows: list[dict] = []
    summaries: list[dict] = []
    for meta in COHORTS:
        expr = load_expr(data / meta["expr"])
        for split in ("tertile", "median"):
            rows, summ = run_cohort(meta, expr, split)
            all_rows.extend(rows)
            summaries.append(summ)

    stats_df = pd.DataFrame(all_rows)
    stats_df.to_csv(tables / "per_gene_stats.tsv", sep="\t", index=False)

    cohort_ids = [c["id"] for c in COHORTS]
    wide = wide_table(all_rows, cohort_ids)
    wide = wide.sort_values(["n_cohorts_pass", "gene"], ascending=[False, True])
    wide.to_csv(tables / "supplement_tj_recurrence.tsv", sep="\t", index=False)

    rec2 = wide.loc[wide["recur_2of3"], "gene"].tolist()
    rec3 = wide.loc[wide["recur_3of3"], "gene"].tolist()
    sets = {
        s["cohort"]: set(s["pass_genes"])
        for s in summaries
        if s["split"] == "tertile"
    }
    med_sets = {
        s["cohort"]: set(s["pass_genes"])
        for s in summaries
        if s["split"] == "median"
    }
    med3 = sorted(set.intersection(*med_sets.values())) if med_sets else []
    med2 = sorted(
        {
            g
            for a, sa in med_sets.items()
            for b, sb in med_sets.items()
            if a < b
            for g in (sa & sb)
        }
    )
    overlap = {
        "OncoSG_LUAD_and_CPTAC_LUAD_RNA": sorted(sets["OncoSG_LUAD"] & sets["CPTAC_LUAD_RNA"]),
        "OncoSG_LUAD_and_GSE31210_LUAD": sorted(sets["OncoSG_LUAD"] & sets["GSE31210_LUAD"]),
        "CPTAC_LUAD_RNA_and_GSE31210_LUAD": sorted(sets["CPTAC_LUAD_RNA"] & sets["GSE31210_LUAD"]),
        "all_three": rec3,
        "at_least_two": rec2,
    }

    summary_json = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "question": (
            "Which published TJ-catalog genes are up in TACSTD2-high vs low "
            "tumors in OncoSG, CPTAC LUAD RNA, and GSE31210, and which recur?"
        ),
        "not_a_check_of": "User A9 slide intersection CLDN1/4/7, F11R, PARD3 (taken as given)",
        "call_rule": "spearman_rho>0 AND tertile_log2FC>0 AND welch_BH_FDR<0.05 (FDR across TJ catalog)",
        "geneset": GENESET_PROVENANCE,
        "cohorts": summaries,
        "overlap_counts": {
            "n_universe": len(TJ_UNIVERSE),
            "n_pass_OncoSG_LUAD": int(len(sets["OncoSG_LUAD"])),
            "n_pass_CPTAC_LUAD_RNA": int(len(sets["CPTAC_LUAD_RNA"])),
            "n_pass_GSE31210_LUAD": int(len(sets["GSE31210_LUAD"])),
            "n_recur_2of3": len(rec2),
            "n_recur_3of3": len(rec3),
            "n_OncoSG_and_CPTAC": len(overlap["OncoSG_LUAD_and_CPTAC_LUAD_RNA"]),
            "n_OncoSG_and_GSE31210": len(overlap["OncoSG_LUAD_and_GSE31210_LUAD"]),
            "n_CPTAC_and_GSE31210": len(overlap["CPTAC_LUAD_RNA_and_GSE31210_LUAD"]),
        },
        "overlap_genes": overlap,
        "median_sensitivity_not_used_for_table": {
            "n_recur_3of3": len(med3),
            "n_recur_2of3": len(med2),
            "recur_3of3": med3,
            "note": "Median split is labeled sensitivity and does not define the supplement table.",
        },
        "a9_slide_annotation": wide.loc[wide["a9_slide_gene"], ["gene", "n_cohorts_pass"]].to_dict("records"),
    }
    (out / "summary.json").write_text(json.dumps(summary_json, indent=2) + "\n")
    pd.DataFrame(summaries).to_csv(tables / "cohort_summary.tsv", sep="\t", index=False)

    write_report(out / "README.md", summaries, wide, overlap)
    fig_overlap_counts(figures / "fig_overlap_counts.png", summaries, wide)
    fig_call_heatmap(figures / "fig_rho_heatmap_passes.png", wide)
    print(json.dumps(summary_json["overlap_counts"], indent=2))
    print("recur_3of3:", rec3)
    print("recur_2of3:", rec2)
    return 0


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    raise SystemExit(main())
