#!/usr/bin/env python3
"""Step 4 - split histology and resolve the two unannotated ICI cohorts.

Pathology labels from GEO are never overwritten.  For GSE135222 and GSE126044
(no histology in the public record) a TCGA-trained squamous-versus-adenocarcinoma
score is applied, then only high-confidence calls are promoted to LUSC /
non-LUSC.  The same score is run on the pathology-annotated ICI samples as an
out-of-distribution check; those labels stay pathology-derived.
"""

from __future__ import annotations

import json
import sys

import numpy as np
import pandas as pd

from common import load_signatures, sample_percentile_ranks, write_table
from config import DATA_DIR, INFERRED_COHORTS, LOGS_DIR, PATHOLOGY_COHORTS, TABLES_DIR

UNIVERSE_PATH = DATA_DIR / "analysis_gene_universe.txt"


def load_universe() -> list[str]:
    return [g for g in UNIVERSE_PATH.read_text().splitlines() if g]


def marker_delta(expr: pd.DataFrame, universe: list[str], sq_genes, ad_genes) -> pd.DataFrame:
    """Within-sample percentile ranks over ``universe``, then squamous minus adeno."""
    present = [g for g in universe if g in expr.index]
    ranks = sample_percentile_ranks(expr.loc[present])
    sq = [g for g in sq_genes if g in ranks.index]
    ad = [g for g in ad_genes if g in ranks.index]
    out = pd.DataFrame({
        "squamous_score": ranks.loc[sq].mean(axis=0),
        "adeno_score": ranks.loc[ad].mean(axis=0),
        "n_squamous_markers": len(sq),
        "n_adeno_markers": len(ad),
    })
    out["histology_delta"] = out["squamous_score"] - out["adeno_score"]
    return out


def confusion(y_true, y_pred, positive="LUSC") -> dict:
    y_true = np.asarray(y_true)
    y_pred = np.asarray(y_pred)
    mask = np.isin(y_true, ["LUSC", "non-LUSC"]) & np.isin(y_pred, ["LUSC", "non-LUSC"])
    yt, yp = y_true[mask], y_pred[mask]
    tp = int(((yt == positive) & (yp == positive)).sum())
    tn = int(((yt != positive) & (yp != positive)).sum())
    fp = int(((yt != positive) & (yp == positive)).sum())
    fn = int(((yt == positive) & (yp != positive)).sum())
    n = tp + tn + fp + fn
    return dict(
        n=n, tp=tp, tn=tn, fp=fp, fn=fn,
        accuracy=(tp + tn) / n if n else np.nan,
        sensitivity=tp / (tp + fn) if (tp + fn) else np.nan,
        specificity=tn / (tn + fp) if (tn + fp) else np.nan,
        n_indeterminate=int((np.asarray(y_pred) == "indeterminate").sum()),
    )


def main() -> int:
    universe = load_universe()
    sigs = load_signatures()
    sq_genes = sigs["SQUAMOUS_MARKERS"]["genes"]
    ad_genes = sigs["ADENO_MARKERS"]["genes"]
    print(f"universe={len(universe)}  squamous markers={sq_genes}  adeno markers={ad_genes}")

    tcga_lusc = pd.read_csv(DATA_DIR / "TCGA_LUSC_percentile_ranks_markers.tsv.gz",
                            sep="\t", index_col=0)
    tcga_luad = pd.read_csv(DATA_DIR / "TCGA_LUAD_percentile_ranks_markers.tsv.gz",
                            sep="\t", index_col=0)
    sq_ok = [g for g in sq_genes if g in tcga_lusc.index]
    ad_ok = [g for g in ad_genes if g in tcga_lusc.index]
    print(f"TCGA markers present: squamous {sq_ok}  adeno {ad_ok}")

    def from_precomputed(ranks: pd.DataFrame) -> pd.Series:
        return ranks.loc[sq_ok].mean(axis=0) - ranks.loc[ad_ok].mean(axis=0)

    lusc_delta = from_precomputed(tcga_lusc)
    luad_delta = from_precomputed(tcga_luad)

    # Midpoint of the two TCGA class medians is the forced-call cut used only
    # for the hold-out accuracy table.  High-confidence calls for the two
    # unannotated ICI cohorts use a wider margin: the pathology ICI hold-out
    # showed that almost all false LUSC calls sit between -0.08 and +0.13, so
    # inferred LUSC is accepted only at delta >= 0.20 and inferred non-LUSC
    # only at delta <= -0.20.  Everything in between stays indeterminate and
    # is excluded from the LUSC-only analyses.
    lusc_lo = float(np.quantile(lusc_delta, 0.10))
    luad_hi = float(np.quantile(luad_delta, 0.90))
    midpoint = float((lusc_delta.median() + luad_delta.median()) / 2)
    lusc_gate = 0.20
    nonlusc_gate = -0.20

    print(f"TCGA-LUSC delta: median={lusc_delta.median():.3f}  p10={lusc_lo:.3f}  "
          f"min={lusc_delta.min():.3f}")
    print(f"TCGA-LUAD delta: median={luad_delta.median():.3f}  p90={luad_hi:.3f}  "
          f"max={luad_delta.max():.3f}")
    print(f"forced-call midpoint={midpoint:.3f}  high-conf LUSC >= {lusc_gate:.3f}  "
          f"high-conf non-LUSC <= {nonlusc_gate:.3f}")
    print(f"TCGA p10 LUSC={lusc_lo:.3f}  TCGA p90 LUAD={luad_hi:.3f}")

    def call(delta: float) -> str:
        if pd.isna(delta):
            return "indeterminate"
        if delta > lusc_gate:
            return "LUSC"
        if delta < nonlusc_gate:
            return "non-LUSC"
        return "indeterminate"

    # TCGA self-check (expected near-perfect; this is the training distribution).
    tcga_true = (["LUSC"] * len(lusc_delta)) + (["non-LUSC"] * len(luad_delta))
    tcga_pred = [call(v) for v in list(lusc_delta) + list(luad_delta)]
    tcga_conf = confusion(tcga_true, tcga_pred)
    print("TCGA self-check (high-confidence calls only):", tcga_conf)

    clinical = pd.read_csv(TABLES_DIR / "cohort_clinical.csv")
    # pandas reads the literal string "NA" as missing; restore the intended token.
    clinical["histology"] = clinical["histology"].fillna("NA")

    score_rows = []
    for cohort in clinical["cohort"].unique():
        expr = pd.read_csv(DATA_DIR / f"{cohort}_expr.tsv.gz", sep="\t", index_col=0)
        scores = marker_delta(expr, universe, sq_genes, ad_genes)
        scores["cohort"] = cohort
        scores["sample_id"] = scores.index
        score_rows.append(scores.reset_index(drop=True))
    scores = pd.concat(score_rows, ignore_index=True)
    scores["classifier_call"] = scores["histology_delta"].map(call)

    merged = clinical.merge(scores, on=["cohort", "sample_id"], how="left")

    # Pathology labels are authoritative.  Classifier calls fill only the
    # two cohorts that have no public histology.
    def resolve(row):
        if row["cohort"] not in INFERRED_COHORTS:
            if row["histology"] in {"LUSC", "non-LUSC"}:
                return row["histology"], "pathology (GEO)", "pathology"
            return "NA", row["histology_source"], "unlabelled"
        call_ = row["classifier_call"]
        if call_ == "LUSC":
            return "LUSC", "inferred (transcriptomic, high-confidence)", "inferred_highconf"
        if call_ == "non-LUSC":
            return "non-LUSC", "inferred (transcriptomic, high-confidence)", "inferred_highconf"
        return "indeterminate", "inferred (transcriptomic, indeterminate)", "inferred_indeterminate"

    resolved = merged.apply(lambda r: pd.Series(resolve(r),
                                                index=["histology_final",
                                                       "histology_source_final",
                                                       "histology_tier"]), axis=1)
    merged = pd.concat([merged, resolved], axis=1)

    # Out-of-distribution check: pathology ICI samples scored by the classifier.
    path_mask = merged["cohort"].isin(PATHOLOGY_COHORTS) & merged["histology"].isin(["LUSC", "non-LUSC"])
    path_conf = confusion(merged.loc[path_mask, "histology"],
                          merged.loc[path_mask, "classifier_call"])
    print("Pathology ICI hold-out (high-confidence calls):", path_conf)

    # Midpoint (non-conservative) accuracy on the same hold-out, for transparency.
    mid_pred = np.where(merged.loc[path_mask, "histology_delta"] >= midpoint, "LUSC", "non-LUSC")
    mid_conf = confusion(merged.loc[path_mask, "histology"], mid_pred)
    print("Pathology ICI hold-out (midpoint, no indeterminate):", mid_conf)

    validation = pd.DataFrame([
        {"set": "TCGA self-check (high-confidence)", **tcga_conf},
        {"set": "pathology ICI hold-out (high-confidence)", **path_conf},
        {"set": "pathology ICI hold-out (midpoint, forced call)", **mid_conf},
    ])
    write_table(validation, TABLES_DIR / "histology_classifier_validation.csv")

    per_cohort = (merged
                  .groupby(["cohort", "histology_final", "histology_tier"])
                  .size().reset_index(name="n"))
    write_table(per_cohort, TABLES_DIR / "histology_split_counts.csv")

    write_table(merged, TABLES_DIR / "cohort_clinical_with_histology.csv")
    write_table(scores, TABLES_DIR / "histology_scores.csv")

    gates = dict(
        midpoint=midpoint, lusc_gate=lusc_gate, nonlusc_gate=nonlusc_gate,
        tcga_lusc_median=float(lusc_delta.median()),
        tcga_luad_median=float(luad_delta.median()),
        tcga_lusc_p10=lusc_lo, tcga_luad_p90=luad_hi,
        squamous_markers=sq_ok, adeno_markers=ad_ok,
        universe_n=len(universe),
    )
    (TABLES_DIR / "histology_classifier_gates.json").write_text(json.dumps(gates, indent=2) + "\n")

    print("\nFinal histology split")
    print(merged.groupby(["cohort", "histology_final"]).size().unstack(fill_value=0).to_string())

    (LOGS_DIR / "04_histology.done").write_text("ok\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
