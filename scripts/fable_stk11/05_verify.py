#!/usr/bin/env python3
"""Independent verification of the fable_stk11 results.

Checks performed
  1. Data-integrity: expected sample counts, no duplicate samples, WT/mut universe.
  2. Biological sanity: STK11/KEAP1/KRAS mutation frequencies within known ranges.
  3. Independent recomputation: recompute the headline TCGA STK11->TACSTD2/CLDN4
     Mann-Whitney p-values straight from the processed matrices and confirm they
     match the analysis-table values within tolerance.
  4. Determinism: re-run analysis logic and confirm identical key numbers.
  5. Response-direction sanity: in GSE135222, immune markers (CD8A/CD274) must
     trend protective (Cox HR < 1) - guards against event-coding sign errors.
  6. Footprint: processed data < 2 GB.
Writes results/fable_stk11/tables/verification_report.json and
notes/fable_stk11/verification.md.
"""
from __future__ import annotations

import os

import numpy as np
import pandas as pd
from scipy import stats

import fable_common as fc

REPORT = {"checks": [], "passed": 0, "failed": 0}


def check(name, ok, detail=""):
    REPORT["checks"].append({"check": name, "status": "PASS" if ok else "FAIL",
                             "detail": detail})
    REPORT["passed" if ok else "failed"] += 1
    print(f"  [{'PASS' if ok else 'FAIL'}] {name} :: {detail}")


def verify_counts():
    print("\n-- 1. data integrity --")
    exp = {"TCGA-LUAD": (566, 510), "TCGA-LUSC": (484, 481),
           "Rizvi2015": (35, None), "Hellmann2018": (75, None),
           "Rizvi2018": (240, None)}
    for lab, (nseq, nexpr) in exp.items():
        geno = pd.read_csv(os.path.join(fc.PROC, f"{lab}_genotype.csv"))
        check(f"{lab} genotype rows≈{nseq}", abs(geno.shape[0] - nseq) <= 3,
              f"got {geno.shape[0]}")
        check(f"{lab} no duplicate samples", geno.sampleId.is_unique,
              f"{geno.sampleId.duplicated().sum()} dups")
        if nexpr is not None:
            expr = pd.read_csv(os.path.join(fc.PROC, f"{lab}_expression.csv"))
            check(f"{lab} expression rows≈{nexpr}", abs(expr.shape[0] - nexpr) <= 3,
                  f"got {expr.shape[0]}")
            check(f"{lab} TACSTD2 & CLDN4 present",
                  {"TACSTD2", "CLDN4"}.issubset(expr.columns), "")


def verify_frequencies():
    print("\n-- 2. mutation frequency sanity (LUAD) --")
    geno = pd.read_csv(os.path.join(fc.PROC, "TCGA-LUAD_genotype.csv"))
    n = geno.shape[0]
    freqs = {g: geno[f"{g}_mut"].mean() for g in ["STK11", "KEAP1", "KRAS", "TP53"]}
    # literature ranges for TCGA LUAD
    ranges = {"STK11": (0.08, 0.22), "KEAP1": (0.10, 0.24),
              "KRAS": (0.25, 0.38), "TP53": (0.40, 0.60)}
    for g, (lo, hi) in ranges.items():
        check(f"LUAD {g} freq in [{lo},{hi}]", lo <= freqs[g] <= hi,
              f"{freqs[g]:.3f} (n={n})")


def verify_recompute():
    print("\n-- 3. independent recomputation (TCGA-LUAD STK11 effect) --")
    expr = pd.read_csv(os.path.join(fc.PROC, "TCGA-LUAD_expression.csv"))
    geno = pd.read_csv(os.path.join(fc.PROC, "TCGA-LUAD_genotype.csv"))
    df = expr.merge(geno, on="sampleId")
    tab = pd.read_csv(os.path.join(fc.TAB, "tcga_genotype_vs_expression.csv"))
    for tgt in ["TACSTD2", "CLDN4"]:
        v = np.log2(df[tgt].clip(lower=0) + 1)
        mut = v[df.STK11_mut == 1]
        wt = v[df.STK11_mut == 0]
        _, p = stats.mannwhitneyu(mut, wt, alternative="two-sided")
        ref = tab[(tab.cohort == "TCGA-LUAD") & (tab.genotype == "STK11") &
                  (tab.target == tgt)]["mwu_p"].iloc[0]
        check(f"recompute STK11->{tgt} p matches table",
              np.isclose(p, ref, rtol=1e-6, atol=1e-12),
              f"recomputed={p:.3e} table={ref:.3e}")
        check(f"STK11->{tgt} direction is lower-in-mutant",
              mut.median() < wt.median(),
              f"med_mut={mut.median():.3f} med_wt={wt.median():.3f}")


def verify_determinism():
    print("\n-- 4. determinism (recompute co-mutation p) --")
    expr = pd.read_csv(os.path.join(fc.PROC, "TCGA-LUAD_expression.csv"))
    geno = pd.read_csv(os.path.join(fc.PROC, "TCGA-LUAD_genotype.csv"))
    df = expr.merge(geno, on="sampleId")
    kras = df[df.KRAS_mut == 1]
    tab = pd.read_csv(os.path.join(fc.TAB, "tcga_kras_comutation_vs_expression.csv"))
    for tgt in ["TACSTD2", "CLDN4"]:
        v = np.log2(kras[tgt].clip(lower=0) + 1)
        _, p = stats.mannwhitneyu(v[kras.STK11_mut == 1], v[kras.STK11_mut == 0],
                                  alternative="two-sided")
        ref = tab[(tab.co_mutation == "STK11") & (tab.target == tgt) &
                  (tab.cohort == "TCGA-LUAD")]["mwu_p"].iloc[0]
        check(f"KL co-mutation ->{tgt} reproducible",
              np.isclose(p, ref, rtol=1e-6, atol=1e-12),
              f"{p:.3e} vs {ref:.3e}")


def verify_response_direction():
    print("\n-- 5. response-direction sanity (GSE135222 immune markers) --")
    cox = pd.read_csv(os.path.join(fc.TAB, "geo_expression_vs_pfs_cox.csv"))
    for g in ["CD8A", "CD274"]:
        hr = cox[cox.target == g]["cox_hr_per_log2"].iloc[0]
        check(f"{g} Cox HR<1 (protective direction)", hr < 1.0, f"HR={hr:.3f}")
    # ICI: STK11 within KRAS should reduce DCB odds (pooled OR<1)
    kp = pd.read_csv(os.path.join(fc.TAB, "ici_kras_restricted_dcb.csv"))
    # recompute pooled MH from tables column
    import ast
    tabs = [ast.literal_eval(t) for t in kp[kp.co_mutation == "STK11"]["table"]]
    import fable_stats as fs
    orr = fs.mantel_haenszel_or(tabs)
    check("KRAS-restricted STK11 pooled OR<1 (worse ICI benefit)", orr < 1.0,
          f"MH-OR={orr:.3f}")


def verify_footprint():
    print("\n-- 6. footprint --")
    mb = fc.dir_size_mb(fc.DATA)
    check("processed+raw data < 2 GB", mb < 2048, f"{mb:.1f} MB")


def main():
    verify_counts()
    verify_frequencies()
    verify_recompute()
    verify_determinism()
    verify_response_direction()
    verify_footprint()

    REPORT["summary"] = f"{REPORT['passed']} passed / {REPORT['failed']} failed"
    fc.save_json(REPORT, os.path.join(fc.TAB, "verification_report.json"))

    lines = ["# Verification report (fable_stk11)", "",
             f"**{REPORT['summary']}**", "",
             "| Check | Status | Detail |", "|---|---|---|"]
    for c in REPORT["checks"]:
        lines.append(f"| {c['check']} | {c['status']} | {c['detail']} |")
    with open(os.path.join(fc.NOTES, "verification.md"), "w") as f:
        f.write("\n".join(lines) + "\n")

    print(f"\n{REPORT['summary']}")
    if REPORT["failed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
