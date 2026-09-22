#!/usr/bin/env python3
"""Patient-level TACSTD2, junction score, and T/NK for the concordant-4 lock.

Three relations are estimated separately. A junction association is not
treated as evidence for either of the other two.

Inclusion (locked before looking at these correlations):
  Human NSCLC scRNA cohorts in the concordant-4 unit lock
  (GSE123902, GSE131907, GSE205335, GSE189357). A unit enters only if it
  is in that lock and has a malignant pseudobulk column. Units are not
  dropped for the sign or the P value of any correlation.

  GSE123902: eligible primary or metastasis donors; marker-defined malignant.
  GSE131907: tumour-origin samples with at least 20 author-defined malignant cells.
  GSE205335: patients in the locked table; author-defined malignant.
  GSE189357: eligible patients; marker-defined malignant.

Junction score (same definition in every cohort):
  MSigDB Hallmark apical junction (HALLMARK_APICAL_JUNCTION), frozen in
  input/a8_sets.json. TACSTD2 is removed if present (it is not a member).
  Malignant UMI sums, genes with count >= 10 in >= 3 units, TMM,
  log2(CPM + 1). Score = unweighted mean of genes from that set that
  remain after the filter. Sensitivity: the same mean after also removing
  CLDN4. A KEGG tight-junction mean uses the same formula and is stored,
  not substituted for the Hallmark score.

Unit of analysis: patient, donor, or sample. Not cells.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
DATA = HERE / "input"
OUT = HERE / "source_data"
OUT.mkdir(parents=True, exist_ok=True)

COHORTS = ["GSE123902", "GSE131907", "GSE205335", "GSE189357"]
TUMOR_ORIGINS = {"tLung", "tL/B", "mLN", "PE", "mBrain"}
ZCRIT = 1.96


def load_units() -> pd.DataFrame:
    rows = []
    d = pd.read_csv(DATA / "GSE123902_marker_units.tsv", sep="\t")
    tumor = d[d["tissue"].isin(["PRIMARY", "METASTASIS"])].copy()
    tumor = tumor.sort_values(["patient", "tissue"]).drop_duplicates("patient", keep="first")
    el = tumor[tumor["eligible"].astype(str).str.lower() == "true"].copy()
    for _, r in el.iterrows():
        rows.append(
            {
                "patient": str(r["patient"]),
                "cohort": "GSE123902",
                "unit": "donor",
                "malig_def": "marker_malig",
                "tissue": str(r["tissue"]),
                "n_cells": int(r["n_cells"]),
                "n_malignant": int(r["n_malignant"]),
                "n_tnk": int(r["n_tnk"]),
                "frac_tnk": float(r["frac_tnk"]),
            }
        )

    d = pd.read_csv(DATA / "GSE131907_samples.tsv", sep="\t")
    mal = d[d["origin"].isin(TUMOR_ORIGINS) & (d["n_malignant"] >= 20)].copy()
    for _, r in mal.iterrows():
        rows.append(
            {
                "patient": str(r["sample"]),
                "cohort": "GSE131907",
                "unit": "sample",
                "malig_def": "author_malig",
                "tissue": str(r["origin"]),
                "n_cells": int(r["n_cells"]),
                "n_malignant": int(r["n_malignant"]),
                "n_tnk": int(r["n_tnk"]),
                "frac_tnk": float(r["frac_tnk"]),
            }
        )

    d = pd.read_csv(DATA / "GSE205335_patients.tsv", sep="\t")
    for _, r in d.iterrows():
        rows.append(
            {
                "patient": str(r["patient"]),
                "cohort": "GSE205335",
                "unit": "patient",
                "malig_def": "author_malig",
                "tissue": str(r["tissue"]),
                "n_cells": int(r["n_cells"]),
                "n_malignant": int(r["n_malignant"]),
                "n_tnk": int(r["n_tnk"]),
                "frac_tnk": float(r["frac_tnk"]),
            }
        )

    d = pd.read_csv(DATA / "GSE189357_marker_units.tsv", sep="\t")
    el = d[d["eligible"].astype(str).str.lower() == "true"].copy()
    for _, r in el.iterrows():
        rows.append(
            {
                "patient": str(r["patient"]),
                "cohort": "GSE189357",
                "unit": "patient",
                "malig_def": "marker_malig",
                "tissue": str(r["tissue"]),
                "n_cells": int(r["n_cells"]),
                "n_malignant": int(r["n_malignant"]),
                "n_tnk": int(r["n_tnk"]),
                "frac_tnk": float(r["frac_tnk"]),
            }
        )
    return pd.DataFrame(rows)


def read_counts(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = df.index.astype(str).str.upper()
    df.columns = df.columns.astype(str)
    return df.groupby(df.index).sum()


def tmm_norm_factors(counts: pd.DataFrame) -> pd.Series:
    lib = counts.sum(axis=0).astype(float).replace(0, np.nan)
    rel = counts.div(lib, axis=1)
    f75 = rel.quantile(0.75, axis=0)
    ref = (f75 - f75.mean()).abs().idxmin()
    ref_c = counts[ref].astype(float)
    ref_lib = float(lib[ref])
    factors = {}
    for col in counts.columns:
        obs = counts[col].astype(float)
        obs_lib = float(lib[col])
        keep = (obs > 0) & (ref_c > 0)
        if int(keep.sum()) < 50:
            factors[col] = 1.0
            continue
        m = np.log2((obs[keep] / obs_lib) / (ref_c[keep] / ref_lib))
        a = 0.5 * np.log2((obs[keep] / obs_lib) * (ref_c[keep] / ref_lib))
        w = (obs_lib - obs[keep]) / (obs_lib * obs[keep]) + (ref_lib - ref_c[keep]) / (
            ref_lib * ref_c[keep]
        )
        ok = np.isfinite(m) & np.isfinite(a) & np.isfinite(w) & (w > 0)
        m, a, w = m[ok], a[ok], w[ok]
        if len(m) < 50:
            factors[col] = 1.0
            continue
        lo_m, hi_m = np.quantile(m, [0.30, 0.70])
        lo_a, hi_a = np.quantile(a, [0.05, 0.95])
        trim = (m >= lo_m) & (m <= hi_m) & (a >= lo_a) & (a <= hi_a)
        if int(trim.sum()) < 20:
            factors[col] = 1.0
            continue
        tmm = float(np.average(m[trim], weights=1.0 / w[trim]))
        factors[col] = 2 ** tmm
    fac = pd.Series(factors)
    return fac / fac.mean()


def log_cpm(counts: pd.DataFrame, factors: pd.Series) -> pd.DataFrame:
    lib = counts.sum(axis=0).astype(float) * factors.reindex(counts.columns).astype(float)
    return np.log2(counts.div(lib, axis=1) * 1e6 + 1.0)


def filter_genes(counts: pd.DataFrame, min_count: int = 10, min_samples: int = 3) -> pd.DataFrame:
    keep = (counts >= min_count).sum(axis=1) >= min_samples
    return counts.loc[keep]


def spearman(x: pd.Series, y: pd.Series) -> tuple[float, float, int]:
    mask = x.notna() & y.notna()
    n = int(mask.sum())
    if n < 3:
        return float("nan"), float("nan"), n
    rho, p = stats.spearmanr(x[mask].to_numpy(float), y[mask].to_numpy(float))
    return float(rho), float(p), n


def fisher_z_dl(rhos: list[float], ns: list[int]) -> dict:
    """DerSimonian–Laird random-effects pool of Fisher-z Spearman coefficients.

    Variance of z is 1/(n-3). Interval uses 1.96 SE. Two-sided normal P.
    This is the pooler that reproduces the Part 2 CLDN4 figure's published
    pooled ρ when given that figure's cohort coefficients.
    """
    rho = np.asarray(rhos, float)
    n = np.asarray(ns, int)
    if np.any(n <= 3):
        raise ValueError("Fisher z requires n > 3")
    if np.any(np.abs(rho) >= 1):
        raise ValueError("Fisher z is undefined at |ρ| = 1")
    z = np.arctanh(rho)
    v = 1.0 / (n - 3)
    w = 1.0 / v
    z_fe = float(np.sum(w * z) / np.sum(w))
    q = float(np.sum(w * (z - z_fe) ** 2))
    df = len(rho) - 1
    c = float(np.sum(w) - np.sum(w**2) / np.sum(w))
    tau2 = max(0.0, (q - df) / c) if c > 0 else 0.0
    w_re = 1.0 / (v + tau2)
    z_re = float(np.sum(w_re * z) / np.sum(w_re))
    se = float(1.0 / np.sqrt(np.sum(w_re)))
    # Hartung–Knapp scale, reported alongside the pre-specified normal P.
    q_hk = float(np.sum(w_re * (z - z_re) ** 2) / df) if df > 0 else float("nan")
    se_hk = se * math_sqrt(q_hk) if q_hk > 0 else float("nan")
    i2 = max(0.0, (q - df) / q) * 100.0 if q > 0 else 0.0
    return {
        "k": int(len(rho)),
        "n": int(np.sum(n)),
        "rho": float(np.tanh(z_re)),
        "ci95_lo": float(np.tanh(z_re - ZCRIT * se)),
        "ci95_hi": float(np.tanh(z_re + ZCRIT * se)),
        "p": float(2.0 * stats.norm.sf(abs(z_re / se))),
        "p_hk_t": float(2.0 * stats.t.sf(abs(z_re / se_hk), df)) if se_hk == se_hk else float("nan"),
        "I2": float(i2),
        "tau2": float(tau2),
        "Q": q,
        "se_z": se,
        "z": z_re,
    }


def math_sqrt(x: float) -> float:
    return float(np.sqrt(x))


def cohort_ci(rho: float, n: int) -> tuple[float, float]:
    if n <= 3 or not np.isfinite(rho) or abs(rho) >= 1:
        return float("nan"), float("nan")
    z = float(np.arctanh(rho))
    se = 1.0 / np.sqrt(n - 3)
    return float(np.tanh(z - ZCRIT * se)), float(np.tanh(z + ZCRIT * se))


def effects_for(units: pd.DataFrame, xcol: str, ycol: str, relation: str, score_name: str) -> list[dict]:
    rows = []
    rhos, ns, members = [], [], []
    for cohort in COHORTS:
        m = units.loc[units["cohort"] == cohort]
        rho, p, n = spearman(m[xcol], m[ycol])
        lo, hi = cohort_ci(rho, n)
        rows.append(
            {
                "kind": "cohort",
                "relation": relation,
                "score": score_name,
                "x": xcol,
                "y": ycol,
                "cohort": cohort,
                "n": n,
                "rho": rho,
                "p": p,
                "ci95_lo": lo,
                "ci95_hi": hi,
                "I2": "",
                "k": 1,
                "tau2": "",
                "p_hk_t": "",
            }
        )
        if n > 3 and np.isfinite(rho) and abs(rho) < 1:
            rhos.append(rho)
            ns.append(n)
            members.append(cohort)
    pool = fisher_z_dl(rhos, ns)
    rows.append(
        {
            "kind": "pooled_dl_fisher_z",
            "relation": relation,
            "score": score_name,
            "x": xcol,
            "y": ycol,
            "cohort": "+".join(members),
            "n": pool["n"],
            "rho": pool["rho"],
            "p": pool["p"],
            "ci95_lo": pool["ci95_lo"],
            "ci95_hi": pool["ci95_hi"],
            "I2": pool["I2"],
            "k": pool["k"],
            "tau2": pool["tau2"],
            "p_hk_t": pool["p_hk_t"],
            "Q": pool["Q"],
            "member_rhos": ",".join(f"{c}:{r:.6f}" for c, r in zip(members, rhos)),
            "member_ns": ",".join(f"{c}:{n}" for c, n in zip(members, ns)),
        }
    )
    return rows


def mean_score(lc: pd.DataFrame, genes: list[str]) -> pd.Series:
    present = [g for g in genes if g in lc.index]
    if not present:
        raise SystemExit("junction score has no genes in the filtered matrix")
    return lc.loc[present].astype(float).mean(axis=0)


def main() -> None:
    a8 = json.loads((DATA / "a8_sets.json").read_text())
    hallmark = [str(g).upper() for g in a8["sets"]["HALLMARK_APICAL_JUNCTION"]]
    kegg = [str(g).upper() for g in a8["sets"]["KEGG_TIGHT_JUNCTION"]]
    if len(hallmark) != len(set(hallmark)):
        raise SystemExit("duplicate Hallmark gene symbols")
    hallmark_primary = [g for g in hallmark if g != "TACSTD2"]
    hallmark_nocldn4 = [g for g in hallmark_primary if g != "CLDN4"]
    kegg_primary = [g for g in kegg if g != "TACSTD2"]
    kegg_nocldn4 = [g for g in kegg_primary if g != "CLDN4"]

    units = load_units()
    if len(units) != 65:
        raise SystemExit(f"locked concordant-4 table changed: n={len(units)}, expected 65")

    mats = [read_counts(DATA / f"{c}_malignant_counts.tsv.gz") for c in COHORTS]
    all_genes = sorted(set().union(*[set(m.index) for m in mats]))
    pieces = []
    for c, m in zip(COHORTS, mats):
        want = [p for p in units.loc[units["cohort"] == c, "patient"] if p in m.columns]
        if not want:
            continue
        pieces.append(m.reindex(all_genes).fillna(0.0)[want])
    counts_raw = pd.concat(pieces, axis=1)
    counts_raw = counts_raw.loc[:, ~counts_raw.columns.duplicated()]
    units["in_count_matrix"] = units["patient"].isin(counts_raw.columns)
    included = units.loc[units["in_count_matrix"]].copy()
    counts = filter_genes(counts_raw)
    factors = tmm_norm_factors(counts)
    lc = log_cpm(counts, factors)
    if "TACSTD2" not in lc.index:
        raise SystemExit("TACSTD2 absent after the expression filter")

    included = included.set_index("patient")
    included["tacstd2"] = lc.loc["TACSTD2", included.index].astype(float)
    included["junction_hallmark"] = mean_score(lc, hallmark_primary).reindex(included.index)
    included["junction_hallmark_no_cldn4"] = mean_score(lc, hallmark_nocldn4).reindex(included.index)
    included["junction_kegg"] = mean_score(lc, kegg_primary).reindex(included.index)
    included["junction_kegg_no_cldn4"] = mean_score(lc, kegg_nocldn4).reindex(included.index)
    included = included.reset_index()

    # Calibration against the prior malignant TACSTD2 log2(TMM-CPM+1) column.
    prior_path = DATA / "prior_tacstd2_log2tmm.tsv"
    if prior_path.exists():
        prior = pd.read_csv(prior_path, sep="\t")
        merged = included.merge(prior[["patient", "tacstd2"]], on="patient", suffixes=("", "_prior"))
        delta = (merged["tacstd2"] - merged["tacstd2_prior"]).abs()
        if len(merged) != len(included) or float(delta.max()) > 1e-8:
            raise SystemExit(
                f"TACSTD2 calibration failed: overlap={len(merged)} max|Δ|={float(delta.max())}"
            )
        cal = pd.DataFrame(
            [
                {
                    "check": "tacstd2_log2_tmm_cpm_vs_prior_pr741",
                    "n": int(len(merged)),
                    "max_abs_diff": float(delta.max()),
                    "status": "match",
                }
            ]
        )
        cal.to_csv(OUT / "calibration_tacstd2.tsv", sep="\t", index=False)

    hm_present = [g for g in hallmark_primary if g in lc.index]
    hm_nocl = [g for g in hallmark_nocldn4 if g in lc.index]
    kg_present = [g for g in kegg_primary if g in lc.index]
    kg_nocl = [g for g in kegg_nocldn4 if g in lc.index]

    gene_rows = []
    for g in sorted(set(hallmark) | set(kegg) | {"TACSTD2", "CLDN4"}):
        gene_rows.append(
            {
                "gene": g,
                "hallmark_apical_junction": g in set(hallmark),
                "kegg_tight_junction": g in set(kegg),
                "in_filtered_matrix": g in lc.index,
                "used_in_hallmark_score": g in set(hm_present),
                "used_in_hallmark_score_no_cldn4": g in set(hm_nocl),
                "used_in_kegg_score": g in set(kg_present),
                "used_in_kegg_score_no_cldn4": g in set(kg_nocl),
            }
        )
    pd.DataFrame(gene_rows).to_csv(OUT / "genes_used.tsv", sep="\t", index=False)

    units_out = units.merge(
        included[
            [
                "patient",
                "tacstd2",
                "junction_hallmark",
                "junction_hallmark_no_cldn4",
                "junction_kegg",
                "junction_kegg_no_cldn4",
            ]
        ],
        on="patient",
        how="left",
    )
    units_out["included"] = units_out["in_count_matrix"]
    units_out["exclude_reason"] = np.where(
        units_out["included"],
        "",
        "locked unit has no malignant pseudobulk column",
    )
    units_out.to_csv(OUT / "units.tsv", sep="\t", index=False)

    effect_rows: list[dict] = []
    specs = [
        ("tacstd2_vs_junction", "junction_hallmark", "tacstd2", "junction_hallmark"),
        ("tacstd2_vs_tnk", "junction_hallmark", "tacstd2", "frac_tnk"),
        ("junction_vs_tnk", "junction_hallmark", "junction_hallmark", "frac_tnk"),
        (
            "tacstd2_vs_junction",
            "junction_hallmark_no_cldn4",
            "tacstd2",
            "junction_hallmark_no_cldn4",
        ),
        (
            "junction_vs_tnk",
            "junction_hallmark_no_cldn4",
            "junction_hallmark_no_cldn4",
            "frac_tnk",
        ),
        ("tacstd2_vs_junction", "junction_kegg", "tacstd2", "junction_kegg"),
        ("junction_vs_tnk", "junction_kegg", "junction_kegg", "frac_tnk"),
        ("tacstd2_vs_junction", "junction_kegg_no_cldn4", "tacstd2", "junction_kegg_no_cldn4"),
        ("junction_vs_tnk", "junction_kegg_no_cldn4", "junction_kegg_no_cldn4", "frac_tnk"),
    ]
    for relation, score_name, xcol, ycol in specs:
        effect_rows.extend(effects_for(included, xcol, ycol, relation, score_name))
    effects = pd.DataFrame(effect_rows)
    effects.to_csv(OUT / "effects.tsv", sep="\t", index=False)

    # Part 2 pooler check (not plotted): published CLDN4 %pos cohort ρ values.
    part2 = fisher_z_dl(
        [-0.6593406593406593, -0.5220779220779221, -0.43534726143421804, -0.6],
        [13, 21, 22, 9],
    )
    if abs(part2["rho"] - (-0.5311678045689989)) > 1e-12:
        raise SystemExit("Fisher-z pooler does not reproduce the Part 2 CLDN4 ρ")

    key = {
        "question": "TROP2 marks a junction-enriched tumour-cell state associated with lymphocyte paucity.",
        "estimands": [
            "Tumour TACSTD2 vs junction program score",
            "Tumour TACSTD2 vs T/NK fraction",
            "Tumour junction program score vs T/NK fraction",
        ],
        "separation": "The three Spearman relations are estimated separately. TACSTD2–junction and TACSTD2–T/NK do not imply junction–T/NK.",
        "part2_not_this_figure": {
            "statement": "Prior CLDN4 percent-positive versus T/NK (Fig. 5 / PR #780, n=65 units, pooled Spearman ρ=-0.531) is Part 2. It is not retitled as TROP2 evidence and is not drawn here.",
            "pooler_reproduction_rho": part2["rho"],
            "pooler_reproduction_p": part2["p"],
            "pooler_reproduction_I2": part2["I2"],
        },
        "inclusion": {
            "cohorts": COHORTS,
            "rule": "concordant-4 locked human NSCLC scRNA units with a malignant pseudobulk column",
            "locked_units": int(len(units)),
            "included_units": int(included.shape[0]),
            "excluded_units": units_out.loc[~units_out["included"], ["patient", "cohort", "exclude_reason"]].to_dict(
                "records"
            ),
            "not_in_this_lock": [
                "GSE148071",
                "GSE127465",
                "GSE154826",
                "GSE200563",
                "E-MTAB-13526",
            ],
            "unit": "patient/donor/sample, never cell n",
            "cohort_n": included.groupby("cohort").size().reindex(COHORTS).astype(int).to_dict(),
        },
        "junction_score": {
            "primary": "MSigDB HALLMARK_APICAL_JUNCTION",
            "source": "input/a8_sets.json sets.HALLMARK_APICAL_JUNCTION (repo freeze of MSigDB Hallmark)",
            "n_genes_in_set": len(hallmark),
            "tacstd2_in_set": "TACSTD2" in hallmark,
            "cldn4_in_set": "CLDN4" in hallmark,
            "n_genes_in_primary_score": len(hm_present),
            "n_genes_in_no_cldn4_score": len(hm_nocl),
            "formula": "mean log2(TMM-CPM+1) of member genes present after count>=10 in >=3 units; TACSTD2 removed",
            "same_definition_across_cohorts": True,
            "sensitivity_no_cldn4": "same mean after removing CLDN4",
            "kegg_sensitivity": {
                "set": "KEGG_TIGHT_JUNCTION",
                "n_genes_in_set": len(kegg),
                "n_genes_in_score": len(kg_present),
                "n_genes_no_cldn4": len(kg_nocl),
                "role": "stored sensitivity; not the primary score",
            },
        },
        "expression": {
            "compartment": "malignant pseudobulk UMI sum",
            "transform": "log2(TMM-CPM+1)",
            "gene_filter": "count >= 10 in >= 3 units, applied once to the merged matrix",
            "tacstd2_calibration": "exact match to prior concordant-4 TACSTD2 log2(TMM-CPM+1)",
        },
        "tnk": "locked frac_tnk on the same unit (T/NK cell fraction). Not recomputed.",
        "statistics": {
            "cohort_effect": "two-sided Spearman ρ (scipy.stats.spearmanr)",
            "pool": "DerSimonian-Laird random effects on Fisher z, variance 1/(n-3), back-transformed to ρ",
            "ci": "1.96 * SE on the z scale, then tanh",
            "p_pool": "two-sided normal",
            "p_hk_t": "Hartung-Knapp t with k-1 degrees of freedom, stored, not drawn",
            "I2": "max(0, (Q-(k-1))/Q)*100 on the Fisher-z scale",
            "multiple_testing": "three pre-specified relations; P values are not adjusted",
            "effect_axis": "Spearman ρ only",
        },
        "effects_primary": effects.loc[effects["score"] == "junction_hallmark"]
        .drop(columns=["member_rhos", "member_ns"], errors="ignore")
        .to_dict("records"),
    }
    # Keep member strings in a compact block.
    key["effects_all_pooled"] = (
        effects.loc[effects["kind"] == "pooled_dl_fisher_z"]
        .assign(rho=lambda d: d["rho"].astype(float))
        .to_dict("records")
    )
    def json_ready(obj):
        if isinstance(obj, dict):
            return {k: json_ready(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [json_ready(v) for v in obj]
        if isinstance(obj, float) and not math.isfinite(obj):
            return None
        return obj

    (HERE / "KEY_STATS.json").write_text(json.dumps(json_ready(key), indent=2, allow_nan=False) + "\n")

    show = effects.loc[
        effects["score"].isin(["junction_hallmark", "junction_hallmark_no_cldn4", "junction_kegg"])
    ]
    cols = ["kind", "relation", "score", "cohort", "n", "rho", "p", "ci95_lo", "ci95_hi", "I2"]
    print(show[cols].to_string(index=False))
    print("included", len(included), "locked", len(units))
    print("excluded", units_out.loc[~units_out["included"], ["patient", "cohort"]].to_dict("records"))
    print("hallmark genes used", len(hm_present), "no cldn4", len(hm_nocl), "kegg", len(kg_present))


if __name__ == "__main__":
    main()
