#!/usr/bin/env python3
"""Correlate TACSTD2 / CLDN4 with immune scores; optional ICI endpoint tests.

Family-wise BH-FDR is applied inside (target x {exclusion, inflamed, other}).
A batch / purity residualisation is run when the corresponding phenotype
columns exist.

Example:

    python scripts/04_correlate.py --indir results/GSE126044 \\
        --group-col response --positive responder \\
        --batch-col tissue_preservation --purity-col estimate:ESTIMATEScore
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from bulkimmune.batch import batch_association, residualise  # noqa: E402
from bulkimmune.stats import (  # noqa: E402
    adjust_pvalues,
    cox_ph,
    partial_spearman,
    spearman_table,
    two_group,
)


def _strip_prefix(col: str) -> str:
    return col.split(":", 1)[-1]


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--indir", type=Path, required=True)
    p.add_argument("--outdir", type=Path, default=None)
    p.add_argument("--group-col", default=None, help="binary phenotype column (e.g. response)")
    p.add_argument("--positive", default=None)
    p.add_argument("--time-col", default=None)
    p.add_argument("--event-col", default=None)
    p.add_argument("--batch-col", default=None)
    p.add_argument("--purity-col", default=None, help="score column used as a purity covariate")
    p.add_argument("--strata-col", default=None, help="run correlations inside each stratum (e.g. cohort)")
    args = p.parse_args()

    outdir = args.outdir or args.indir
    outdir.mkdir(parents=True, exist_ok=True)

    scores = pd.read_csv(args.indir / "scores.tsv", sep="\t", index_col=0)
    targets = pd.read_csv(args.indir / "targets.tsv", sep="\t", index_col=0)
    pheno_path = args.indir / "pheno.tsv"
    pheno = pd.read_csv(pheno_path, sep="\t", index_col=0) if pheno_path.exists() else pd.DataFrame(index=scores.index)

    # Use unprefixed names for family assignment, keep prefixed names in the table
    scores_for_family = scores.copy()
    scores_for_family.columns = [_strip_prefix(c) for c in scores.columns]
    # but we need unique columns -- keep prefixed in the reported table
    corr = spearman_table(targets, scores)
    corr["score_short"] = corr["score"].map(_strip_prefix)
    corr["family"] = None  # filled by adjust_pvalues via score_short
    # map family from the short name
    from bulkimmune.stats import assign_family

    corr["family"] = corr["score_short"].map(assign_family)
    corr = adjust_pvalues(corr, p_col="p", group_cols=("target",), family_col="family")
    corr = corr.sort_values(["target", "family", "p"])
    corr.to_csv(outdir / "correlation_spearman.tsv", sep="\t", index=False)

    # partial correlations
    cov_parts = []
    if args.batch_col and args.batch_col in pheno.columns:
        cov_parts.append(pheno[[args.batch_col]])
    if args.purity_col and args.purity_col in scores.columns:
        cov_parts.append(scores[[args.purity_col]].rename(columns={args.purity_col: "purity"}))
    if args.strata_col and args.strata_col in pheno.columns:
        # include cohort as a covariate when not stratifying the whole analysis
        pass
    partial_rows = []
    if cov_parts:
        cov = pd.concat(cov_parts, axis=1)
        for tcol in targets.columns:
            for scol in scores.columns:
                res = partial_spearman(targets[tcol], scores[scol], cov)
                partial_rows.append({"target": tcol, "score": scol, **res})
        part = pd.DataFrame(partial_rows)
        part["score_short"] = part["score"].map(_strip_prefix)
        part["family"] = part["score_short"].map(assign_family)
        part = adjust_pvalues(part, p_col="p", group_cols=("target",), family_col="family")
        part.to_csv(outdir / "correlation_partial.tsv", sep="\t", index=False)

    if args.strata_col and args.strata_col in pheno.columns:
        pieces = []
        for level, idx in pheno.groupby(args.strata_col).groups.items():
            sub = spearman_table(targets.loc[idx], scores.loc[idx])
            sub["stratum"] = level
            sub["score_short"] = sub["score"].map(_strip_prefix)
            sub["family"] = sub["score_short"].map(assign_family)
            sub = adjust_pvalues(sub, p_col="p", group_cols=("target",), family_col="family")
            pieces.append(sub)
        pd.concat(pieces).to_csv(outdir / "correlation_by_stratum.tsv", sep="\t", index=False)

    if args.batch_col and args.batch_col in pheno.columns:
        ba = batch_association(scores, pheno[args.batch_col])
        ba = adjust_pvalues(ba, p_col="p", group_cols=(), family_col=None)
        ba.to_csv(outdir / "batch_association.tsv", sep="\t", index=False)

    if args.group_col and args.group_col in pheno.columns:
        tg = two_group(pd.concat([targets, scores], axis=1), pheno[args.group_col], positive=args.positive)
        tg["score_short"] = tg["score"].map(lambda s: _strip_prefix(s) if ":" in str(s) else s)
        tg["family"] = tg["score_short"].map(assign_family)
        tg = adjust_pvalues(tg, p_col="p", group_cols=(), family_col="family")
        tg.to_csv(outdir / "endpoint_twogroup.tsv", sep="\t", index=False)

    if args.time_col and args.event_col and args.time_col in pheno.columns:
        surv = cox_ph(
            pd.concat([targets, scores], axis=1),
            pheno[args.time_col],
            pheno[args.event_col],
        )
        surv["score_short"] = surv["score"].map(lambda s: _strip_prefix(s) if ":" in str(s) else s)
        surv["family"] = surv["score_short"].map(assign_family)
        surv = adjust_pvalues(surv, p_col="p", group_cols=(), family_col="family")
        surv.to_csv(outdir / "endpoint_survival.tsv", sep="\t", index=False)

    # compact highlight table: exclusion + inflamed, both targets
    highlight = corr[corr["family"].isin(["exclusion", "inflamed"])].copy()
    highlight.to_csv(outdir / "correlation_highlight.tsv", sep="\t", index=False)

    summary = {
        "n_samples": int(scores.shape[0]),
        "n_scores": int(scores.shape[1]),
        "n_tests": int(corr.shape[0]),
        "n_fdr_lt_0.05": int((corr["p_adj"] < 0.05).sum()) if "p_adj" in corr else 0,
        "top_by_target": {},
    }
    for tcol, sub in corr.groupby("target"):
        top = sub.dropna(subset=["spearman_r"]).assign(abs_r=lambda d: d["spearman_r"].abs())
        top = top.sort_values("abs_r", ascending=False).head(8)
        summary["top_by_target"][tcol] = top[["score", "family", "n", "spearman_r", "p", "p_adj"]].to_dict(
            orient="records"
        )
    (outdir / "correlate_summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps(summary, indent=2, default=str))


if __name__ == "__main__":
    main()
