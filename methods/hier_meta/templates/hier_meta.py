#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
hier_meta.py -- Bayesian / hierarchical pooling template for TACSTD2 (TROP2) and
CLDN4 effects across many small immune-checkpoint-inhibitor (ICI) GEO cohorts.

OPEN COHORTS ONLY. Rows whose access is EGA / dbGaP / DAC-controlled are never
pooled. They belong in the missing-cohort registry and enter the analysis only as
MNAR sensitivity, never as downloaded data.

METHODS TEMPLATE. This script contains no findings. It reads per-cohort effect
estimates from `results/*.csv` **if they are present** and, when they are not,
prints the required input contract and exits without doing anything.

What it does when inputs are present
------------------------------------
  1. Harmonises per-cohort estimates (effect scale, unit, sign convention, overlap).
  2. Fits a normal-normal hierarchical (random-effects) model by exact numerical
     integration over tau -- no MCMC, fully deterministic and reproducible.
  3. Reports partial pooling / shrinkage per cohort, tau, and a prediction interval.
  4. Runs frequentist companions (FE, DL, PM, REML, Hartung-Knapp) as cross-checks.
  5. Runs small-study / selection-bias diagnostics (Egger, Begg, trim-and-fill,
     PET-PEESE, step-function selection model).
  6. Runs the missing-cohort (EGA / dbGaP data-access-controlled trials) sensitivity
     analysis: MNAR delta sweep with a tipping point, Manski-style bounds, and the
     number of unobserved null cohorts required to overturn the conclusion.
  7. Writes every intermediate quantity to CSV plus a run manifest for reproducibility.

Design notes
------------
  * Conditional on tau, the model is Gaussian-conjugate, so the posterior of mu, of each
    cohort effect theta_i, and of a future cohort effect are exact finite mixtures over a
    tau grid. That removes sampler variance from a problem where K is small and the
    answer is prior-sensitive: the reported prior sensitivity is then a real property of
    the model rather than Monte Carlo noise.
  * The level-1 covariance V is allowed to be non-diagonal so that cohorts sharing
    patients (e.g. pre-treatment and on-treatment series from one trial) can be modelled
    as correlated rather than silently double-counted.
  * Everything that a reviewer would ask to see -- the tau posterior, per-cohort
    shrinkage weights, leave-one-out refits, the prior sensitivity table -- is written
    out, not just printed.

Dependencies: numpy, scipy (stdlib csv/json/argparse otherwise). matplotlib optional.

Usage
-----
    python3 hier_meta.py --results-dir results --outdir methods/hier_meta/out
    python3 hier_meta.py --write-config my_config.json     # dump the default config
    python3 hier_meta.py --config my_config.json --results-dir results
    python3 hier_meta.py --demo                            # run on the bundled synthetic example
"""

from __future__ import annotations

import argparse
import copy
import csv
import glob
import hashlib
import json
import math
import os
import platform as _platform
import sys
from dataclasses import dataclass, field, asdict
from datetime import datetime, timezone
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
from scipy import optimize, stats

__version__ = "0.1.0"

# numpy >= 2 renamed trapz; keep the template runnable on either.
_trapz = getattr(np, "trapezoid", None) or np.trapz


def _quad_weights(x: np.ndarray) -> np.ndarray:
    """Trapezoid quadrature weights, normalised to sum to 1."""
    w = np.zeros_like(x)
    d = np.diff(x)
    w[:-1] += d / 2.0
    w[1:] += d / 2.0
    return w

# --------------------------------------------------------------------------------------
# 0. Configuration
# --------------------------------------------------------------------------------------

DEFAULT_CONFIG: Dict[str, Any] = {
    # One run pools one gene. Pass several via --genes TACSTD2,CLDN4 (writes one
    # subdirectory per gene). Defaults cover both project genes.
    "gene": "TACSTD2",
    "genes": ["TACSTD2", "CLDN4"],
    # Endpoint actually pooled. Effects from different endpoint families are NOT pooled
    # together; run the script once per endpoint.
    "primary_endpoint": "DCB",
    # "smd": Hedges' g / Cohen's d of expression (responder minus non-responder).
    # "log_ratio": log HR or log OR per unit (survival / logistic association).
    # These two families are not on the same scale and must never be pooled together.
    "effect_family": "smd",
    # Effect unit that the pooled estimate is expressed in (log_ratio family).
    "effect_unit": "per_sd",
    # "harm": positive = higher gene expression associated with WORSE ICI outcome.
    # For SMD rows stored as responder-minus-nonresponder this flips the sign.
    # Response-type log-OR rows are also flipped so every endpoint shares one direction.
    "orientation": "harm",
    # "open_only" (default): drop EGA/dbGaP/controlled rows from the likelihood.
    # "catalog": keep them (not recommended; reserved for a DAC-approved re-run).
    "access_policy": "open_only",
    # Convert per-log2 / dichotomised effects onto the per-SD scale (see UNIT_TO_SD).
    # Setting this to false is the conservative choice: non-conforming rows are dropped.
    "allow_unit_conversion": True,
    # Which adjustment set to pool. Rows whose `model_adjust` differs are dropped.
    "model_adjust": "univariable",
    # Cohorts that share patients: "drop" keeps one row per `sample_group` by
    # `overlap_precedence`; "correlate" keeps all and gives them within-group
    # correlation `overlap_rho` in the level-1 covariance.
    "overlap_policy": "drop",
    "overlap_rho": 0.5,
    "overlap_precedence": ["pre_treatment", "baseline", "on_treatment", "post_treatment"],
    "prior": {
        "mu_mean": 0.0,
        "mu_sd": 1.0,
        "tau_family": "half_normal",
        "tau_scale": 0.5,
    },
    # Every prior in this list is refitted and reported side by side. With K small the
    # tau prior is doing real work and this table is a required output, not an extra.
    "prior_sensitivity": [
        {"label": "primary_hn0.5", "mu_mean": 0.0, "mu_sd": 1.0,
         "tau_family": "half_normal", "tau_scale": 0.5},
        {"label": "tight_hn0.25", "mu_mean": 0.0, "mu_sd": 1.0,
         "tau_family": "half_normal", "tau_scale": 0.25},
        {"label": "wide_hn1.0", "mu_mean": 0.0, "mu_sd": 1.0,
         "tau_family": "half_normal", "tau_scale": 1.0},
        {"label": "half_cauchy0.5", "mu_mean": 0.0, "mu_sd": 1.0,
         "tau_family": "half_cauchy", "tau_scale": 0.5},
        {"label": "uniform_tau", "mu_mean": 0.0, "mu_sd": 1.0,
         "tau_family": "uniform", "tau_scale": 2.0},
        {"label": "skeptical_mu", "mu_mean": 0.0, "mu_sd": 0.35,
         "tau_family": "half_normal", "tau_scale": 0.5},
        {"label": "flat_mu", "mu_mean": 0.0, "mu_sd": 100.0,
         "tau_family": "half_normal", "tau_scale": 0.5},
        # Empirical priors from published predictive distributions for tau^2
        # (Turner et al. for log-OR, Rhodes et al. for continuous/HR-type outcomes)
        # must be filled in from the relevant table cell before use; the entry below is
        # a placeholder that is skipped unless `enabled` is set to true.
        {"label": "empirical_tau2_lognormal", "enabled": False, "mu_mean": 0.0,
         "mu_sd": 1.0, "tau_family": "lognormal_tau2",
         "log_tau2_mean": None, "log_tau2_sd": None},
    ],
    # Region of practical equivalence for the pooled effect. For effect_family=smd
    # this is a small-effect Hedges' g (0.2). For log_ratio it is log(1.1) ≈ 0.095.
    "rope": 0.2,
    "tau_grid": {"max": None, "n": 801},
    "publication_bias": {
        # "sei" = classic Egger; "inv_n" = sample-size based variant (Peters/Macaskill
        # style), which is preferable for log-OR because log-OR and its SE are
        # structurally correlated.
        "egger_predictor": "sei",
        "selection_alpha": 0.025,
        "omega_grid": [1.0, 0.5, 0.25, 0.1],
        "trimfill_estimator": "L0",
        # Optional CSV of pooled results for a null / negative-control gene panel used to
        # calibrate the gene-selection winner's curse (see playbook section 7.3).
        "null_gene_results": "",
    },
    "missingness": {
        # Cohort registry listing every known cohort, accessible or not.
        "registry": "",
        "delta_grid": [-1.0, -0.8, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0,
                       0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0],
        "decision_threshold": 0.95,
        # "auto": SE for an unobserved cohort from events (time-to-event) or from n and
        # the response rate (binary); otherwise the median observed SE is used.
        "se_rule": "auto",
        "max_failsafe_k": 200,
    },
    # Exploratory only. Rule of thumb: at most one moderator per ~10 cohorts.
    "meta_regression": {"moderators": []},
    "quantiles": [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975],
}

# Multiplier converting a dichotomised / categorised contrast to a per-1-SD slope,
# assuming the underlying expression is approximately normal and the effect is linear in
# expression. For a standard normal X, the mean difference between the groups defined by
# a median split is 2*phi(0)/0.5 = 1.596 SD, etc. Dividing the contrast by this factor
# puts it on the per-SD scale; the SE is divided by the same factor, so the z-statistic
# (and hence the efficiency loss from dichotomising) is preserved.
UNIT_TO_SD: Dict[str, float] = {
    "per_sd": 1.0,
    "per_z": 1.0,
    "hedges_g": 1.0,
    "smd": 1.0,
    "cohens_d": 1.0,
    "high_vs_low": 1.0 / 1.5957691216057308,
    "median_split": 1.0 / 1.5957691216057308,
    "tertile_top_vs_bottom": 1.0 / 2.181916,
    "quartile_top_vs_bottom": 1.0 / 2.541560,
}
# Units needing a cohort-specific SD of the expression variable (column `x_sd`).
UNIT_NEEDS_XSD = {"per_log2", "per_unit", "per_tpm", "per_log2_tpm"}

RESPONSE_ENDPOINTS = {"ORR", "RESPONSE", "DCB", "CB", "RECIST", "PR_CR", "OBJECTIVE_RESPONSE",
                      "MPR", "pCR", "PCR", "MAJOR_PATHOLOGIC_RESPONSE"}
SURVIVAL_ENDPOINTS = {"OS", "PFS", "EFS", "DFS", "TTP", "TTF"}

# Access tokens treated as open vs restricted. Matching is punctuation-insensitive.
OPEN_ACCESS = {
    "open", "public", "geo", "arrayexpress", "ae", "zenodo", "pride",
    "processedopen", "openprocessed", "ftp", "supplement",
}
RESTRICTED_ACCESS = {
    "ega", "egad", "egas", "dbgap", "phs", "controlled", "restricted",
    "dac", "controlledaccess", "accesscontrolled", "egacontrolled",
}
OPEN_ID_PREFIXES = ("GSE", "GDS", "GPL", "GSM", "E-MTAB-", "E-GEOD-", "E-MEXP-",
                    "E-TABM-", "SYNTH_OPEN")

# --------------------------------------------------------------------------------------
# 1. Input contract
# --------------------------------------------------------------------------------------

# Canonical column -> accepted aliases (case-insensitive, punctuation-insensitive).
COLUMN_ALIASES: Dict[str, Sequence[str]] = {
    "cohort_id": ("cohort_id", "cohort", "dataset", "series", "gse", "accession",
                  "study_id", "id"),
    "study": ("study", "first_author", "author", "reference", "citation"),
    "gene": ("gene", "gene_symbol", "symbol", "feature", "hgnc_symbol"),
    "endpoint": ("endpoint", "outcome", "outcome_type", "response_variable", "event"),
    "unit": ("unit", "effect_unit", "scale", "exposure_unit", "predictor_unit"),
    "effect_measure": ("effect_measure", "measure", "effect_type", "statistic_type"),
    "yi": ("yi", "estimate", "beta", "coef", "coefficient", "log_hr", "loghr",
           "log_or", "logor", "logodds", "effect", "effect_size", "b",
           "hedges_g", "hedgesg", "smd", "cohens_d", "cohensd"),
    "sei": ("sei", "se", "std_error", "stderr", "standard_error", "se_log_hr",
            "se_beta", "se_log_or", "sd_error"),
    "ci_low": ("ci_low", "ci_lower", "lower", "lcl", "l95", "ci95_low", "conf_low",
               "hr_lower", "or_lower", "ci_low_95"),
    "ci_high": ("ci_high", "ci_upper", "upper", "ucl", "u95", "ci95_high", "conf_high",
                "hr_upper", "or_upper", "ci_high_95"),
    "ci_level": ("ci_level", "conf_level", "ci_width", "level"),
    "hr": ("hr", "hazard_ratio", "hazardratio"),
    "or": ("or", "odds_ratio", "oddsratio"),
    "pvalue": ("pvalue", "p_value", "p", "pval", "p_val"),
    "ni": ("ni", "n", "n_total", "sample_size", "n_patients", "nsamples", "n_samples"),
    "events": ("events", "n_events", "nevent", "event_count", "deaths", "n_death"),
    "responders": ("responders", "n_responders", "n_response", "n_cr_pr", "responder_n"),
    "x_sd": ("x_sd", "sd_x", "expr_sd", "sd_expression", "sd_log2", "predictor_sd"),
    "model_adjust": ("model_adjust", "model", "adjustment", "covariates", "adjusted"),
    "cancer_type": ("cancer_type", "tumour_type", "tumor_type", "indication", "cancer"),
    "ici_class": ("ici_class", "ici", "drug_class", "agent_class", "treatment_class"),
    "ici_target": ("ici_target", "target", "checkpoint", "drug_target"),
    "treatment_line": ("treatment_line", "line", "line_of_therapy", "therapy_line"),
    "platform": ("platform", "assay", "technology", "seq_platform", "array"),
    "tissue_timing": ("tissue_timing", "timing", "biopsy_timing", "sample_timing"),
    "sample_group": ("sample_group", "patient_group", "overlap_group", "cluster",
                     "trial", "parent_trial"),
    "access": ("access", "access_type", "data_access", "availability", "repo",
               "repository", "source_access"),
    "notes": ("notes", "note", "comment", "comments"),
}

REQUIRED_MESSAGE = """
No usable per-cohort estimates were found.

This template pools estimates that another step of the pipeline has already produced.
It expects one row per (cohort x gene x endpoint) in one or more CSV files under the
results directory. Minimum contract:

  cohort_id   cohort / GEO series identifier                       (required)
  gene        gene symbol, e.g. TACSTD2 or CLDN4                   (required)
  endpoint    DCB | ORR | OS | PFS | MPR ...                       (required)
  unit        hedges_g | smd | per_sd | per_log2 | high_vs_low     (required)
  yi          effect (Hedges' g, or log HR / log OR)               (required*)
  sei         standard error of yi                                 (required*)
  access      open | geo | ega | dbgap | controlled                (recommended)

  * instead of (yi, sei) you may supply any of:
      hr / or  + ci_low + ci_high        (ratio scale, back-transformed for you)
      yi       + ci_low + ci_high        (log scale)
      yi       + pvalue                  (two-sided; least preferred, recorded as such)

Recommended additional columns, used by the heterogeneity, meta-regression and
missing-cohort modules:

  ni, events, responders, x_sd, model_adjust, cancer_type, ici_class, ici_target,
  treatment_line, platform, tissue_timing, sample_group, study, notes

`x_sd` (the within-cohort SD of the expression variable actually used in the model) is
what makes a per-log2 slope convertible to the per-SD scale. Without it, per-log2 rows
cannot be pooled with per-SD rows and are dropped.

Nothing was written. Re-run once the per-cohort step has produced results/*.csv, or run
with --demo to exercise the pipeline on the bundled synthetic example.
"""


def _norm_key(s: str) -> str:
    return "".join(ch for ch in s.strip().lower() if ch.isalnum())


_ALIAS_LOOKUP: Dict[str, str] = {}
for _canon, _aliases in COLUMN_ALIASES.items():
    for _a in _aliases:
        _ALIAS_LOOKUP[_norm_key(_a)] = _canon


def _to_float(x: Any) -> Optional[float]:
    if x is None:
        return None
    s = str(x).strip()
    if s == "" or s.lower() in {"na", "nan", "null", "none", "n/a", ".", "-"}:
        return None
    s = s.replace(",", "")
    try:
        v = float(s)
    except ValueError:
        return None
    return v if math.isfinite(v) else None


def _clean_str(x: Any) -> str:
    return "" if x is None else str(x).strip()


def classify_access(access_raw: str, cohort_id: str) -> str:
    """Return 'open', 'restricted', or 'unknown'.

    Restricted tokens always win, even if the cohort_id looks like a GEO series
    (e.g. GSE205335 raw listed as EGAD). An empty access field is treated as open
    only when the identifier itself is a public GEO / ArrayExpress / synthetic-open
    accession; otherwise it is unknown and `open_only` drops it.
    """
    tok = _norm_key(access_raw)
    if tok in RESTRICTED_ACCESS or any(tok.startswith(p) for p in ("ega", "egad", "egas", "phs", "dbgap")):
        return "restricted"
    if tok in OPEN_ACCESS:
        return "open"
    cid = (cohort_id or "").strip().upper()
    if any(cid.startswith(p) for p in OPEN_ID_PREFIXES):
        return "open"
    if cid.startswith("EGA") or cid.startswith("PHS"):
        return "restricted"
    return "unknown"


@dataclass
class Cohort:
    """One harmonised per-cohort estimate."""
    cohort_id: str
    gene: str = ""
    endpoint: str = ""
    unit: str = ""
    yi: float = float("nan")
    sei: float = float("nan")
    yi_raw: float = float("nan")
    sei_raw: float = float("nan")
    unit_factor: float = 1.0
    se_source: str = ""
    sign_flipped: bool = False
    ni: Optional[float] = None
    events: Optional[float] = None
    responders: Optional[float] = None
    x_sd: Optional[float] = None
    study: str = ""
    model_adjust: str = ""
    cancer_type: str = ""
    ici_class: str = ""
    ici_target: str = ""
    treatment_line: str = ""
    platform: str = ""
    tissue_timing: str = ""
    sample_group: str = ""
    access: str = ""
    notes: str = ""
    source_file: str = ""
    source_row: int = -1


def read_result_files(paths: Sequence[str]) -> List[Dict[str, Any]]:
    """Read CSVs and map their headers onto canonical column names."""
    rows: List[Dict[str, Any]] = []
    for path in paths:
        with open(path, newline="", encoding="utf-8-sig") as fh:
            reader = csv.DictReader(fh)
            if reader.fieldnames is None:
                continue
            mapping = {}
            for raw in reader.fieldnames:
                canon = _ALIAS_LOOKUP.get(_norm_key(raw))
                if canon is not None and canon not in mapping.values():
                    mapping[raw] = canon
            for idx, raw_row in enumerate(reader, start=2):  # header is line 1
                rec: Dict[str, Any] = {"__file__": path, "__line__": idx}
                for raw, canon in mapping.items():
                    rec[canon] = raw_row.get(raw)
                rows.append(rec)
    return rows


def _derive_effect(rec: Dict[str, Any], audit: Dict[str, str],
                   effect_family: str = "log_ratio") -> Tuple[Optional[float], Optional[float], str]:
    """Return (yi, sei, se_source) from whatever the row provides.

    For effect_family=log_ratio the value is forced onto the log-HR / log-OR scale.
    For effect_family=smd the value is left as a standardized mean difference; ratio
    columns are ignored so a stray HR cannot silently enter an SMD pool.
    """
    yi = _to_float(rec.get("yi"))
    sei = _to_float(rec.get("sei"))
    lo, hi = _to_float(rec.get("ci_low")), _to_float(rec.get("ci_high"))
    level = _to_float(rec.get("ci_level")) or 0.95
    if level > 1:  # given as 95 rather than 0.95
        level /= 100.0
    zc = stats.norm.ppf(0.5 + level / 2.0)

    if effect_family == "smd":
        ratio = None
        ratio_scale = False
        measure = _clean_str(rec.get("effect_measure")).lower()
        if measure in {"hr", "or", "hazard_ratio", "odds_ratio", "log_hr", "log_or"}:
            audit["yi"] = "skipped: ratio measure in smd family"
            return None, None, ""
    else:
        ratio = None
        for key in ("hr", "or"):
            v = _to_float(rec.get(key))
            if v is not None and v > 0:
                ratio = v
                break
        measure = _clean_str(rec.get("effect_measure")).lower()
        ratio_scale = ratio is not None or measure in {"hr", "or", "hazard_ratio", "odds_ratio"}

    if yi is None and ratio is not None:
        yi = math.log(ratio)
        if lo is not None and hi is not None and lo > 0 and hi > 0:
            lo, hi = math.log(lo), math.log(hi)
        audit["yi"] = "log(ratio)"
    elif yi is not None and ratio_scale and yi > 0 and measure in {"hr", "or"}:
        # A column named `estimate` alongside effect_measure=HR is on the ratio scale.
        yi = math.log(yi)
        if lo is not None and hi is not None and lo > 0 and hi > 0:
            lo, hi = math.log(lo), math.log(hi)
        audit["yi"] = "log(ratio)"

    if yi is None:
        return None, None, ""

    if sei is not None and sei > 0:
        return yi, sei, "reported_se"
    if lo is not None and hi is not None and hi > lo:
        return yi, (hi - lo) / (2.0 * zc), f"ci{int(round(level * 100))}"
    p = _to_float(rec.get("pvalue"))
    if p is not None and 0 < p < 1 and yi != 0:
        z = stats.norm.ppf(1.0 - p / 2.0)
        if z > 0:
            return yi, abs(yi) / z, "pvalue"
    return yi, None, ""


def harmonise(rows: Sequence[Dict[str, Any]], cfg: Dict[str, Any]) -> Tuple[List[Cohort], List[Dict[str, str]]]:
    """Filter and convert raw rows into a list of comparable cohort estimates.

    Every input row is recorded in the audit trail with the reason it was kept or
    dropped, so the count that reaches the model can be reconciled with the count that
    left the per-cohort step.
    """
    gene_want = _norm_key(cfg["gene"])
    endpoint_want = _norm_key(cfg["primary_endpoint"])
    adjust_want = _norm_key(cfg.get("model_adjust") or "")
    kept: List[Cohort] = []
    audit: List[Dict[str, str]] = []

    for rec in rows:
        a = {"file": os.path.basename(str(rec.get("__file__"))),
             "line": str(rec.get("__line__")),
             "cohort_id": _clean_str(rec.get("cohort_id")),
             "gene": _clean_str(rec.get("gene")),
             "endpoint": _clean_str(rec.get("endpoint")),
             "unit": _clean_str(rec.get("unit")),
             "access": _clean_str(rec.get("access")),
             "kept": "no", "reason": ""}

        if not a["cohort_id"]:
            a["reason"] = "missing cohort_id"
            audit.append(a); continue
        access_class = classify_access(_clean_str(rec.get("access")), a["cohort_id"])
        if cfg.get("access_policy", "open_only") == "open_only" and access_class != "open":
            a["reason"] = (f"access={access_class} rejected by open_only policy "
                           "(EGA/dbGaP/controlled rows stay in the missingness registry)")
            audit.append(a); continue
        if gene_want and _norm_key(a["gene"]) != gene_want:
            a["reason"] = f"gene != {cfg['gene']}"
            audit.append(a); continue
        if endpoint_want and _norm_key(a["endpoint"]) != endpoint_want:
            a["reason"] = f"endpoint != {cfg['primary_endpoint']}"
            audit.append(a); continue
        if adjust_want:
            got = _norm_key(_clean_str(rec.get("model_adjust")))
            if got and got != adjust_want:
                a["reason"] = f"model_adjust != {cfg['model_adjust']}"
                audit.append(a); continue

        conv: Dict[str, str] = {}
        yi, sei, se_source = _derive_effect(rec, conv, cfg.get("effect_family", "log_ratio"))
        if yi is None:
            a["reason"] = "no effect estimate"
            audit.append(a); continue
        if sei is None or not (sei > 0):
            a["reason"] = "no usable standard error"
            audit.append(a); continue

        unit = _norm_key(_clean_str(rec.get("unit")) or cfg["effect_unit"])
        x_sd = _to_float(rec.get("x_sd"))
        factor = None
        if unit in {_norm_key(u) for u in UNIT_TO_SD}:
            key = next(u for u in UNIT_TO_SD if _norm_key(u) == unit)
            factor = UNIT_TO_SD[key]
            if factor != 1.0 and not cfg.get("allow_unit_conversion", True):
                a["reason"] = f"unit {unit} needs conversion (disabled)"
                audit.append(a); continue
        elif unit in {_norm_key(u) for u in UNIT_NEEDS_XSD}:
            if not cfg.get("allow_unit_conversion", True):
                a["reason"] = f"unit {unit} needs conversion (disabled)"
                audit.append(a); continue
            if x_sd is None or x_sd <= 0:
                a["reason"] = f"unit {unit} requires x_sd to reach per-SD scale"
                audit.append(a); continue
            factor = x_sd
        if factor is None:
            a["reason"] = f"unrecognised unit '{unit}'"
            audit.append(a); continue

        yi_c, sei_c = yi * factor, sei * factor

        flipped = False
        ep = _clean_str(rec.get("endpoint")).upper()
        if cfg.get("orientation") == "harm" and ep in RESPONSE_ENDPOINTS:
            # Positive must mean "worse ICI outcome". A positive SMD (responder minus
            # non-responder) or a positive log-OR for response means *better* outcome.
            yi_c = -yi_c
            flipped = True

        kept.append(Cohort(
            cohort_id=a["cohort_id"],
            gene=a["gene"], endpoint=a["endpoint"], unit=unit,
            yi=yi_c, sei=sei_c, yi_raw=yi, sei_raw=sei, unit_factor=factor,
            se_source=se_source, sign_flipped=flipped,
            ni=_to_float(rec.get("ni")), events=_to_float(rec.get("events")),
            responders=_to_float(rec.get("responders")), x_sd=x_sd,
            study=_clean_str(rec.get("study")),
            model_adjust=_clean_str(rec.get("model_adjust")),
            cancer_type=_clean_str(rec.get("cancer_type")),
            ici_class=_clean_str(rec.get("ici_class")),
            ici_target=_clean_str(rec.get("ici_target")),
            treatment_line=_clean_str(rec.get("treatment_line")),
            platform=_clean_str(rec.get("platform")),
            tissue_timing=_clean_str(rec.get("tissue_timing")),
            sample_group=_clean_str(rec.get("sample_group")),
            access=access_class,
            notes=_clean_str(rec.get("notes")),
            source_file=str(rec.get("__file__")), source_row=int(rec.get("__line__", -1)),
        ))
        a["kept"] = "yes"
        a["reason"] = f"unit_factor={factor:.4g}; se_source={se_source}" + ("; sign flipped" if flipped else "")
        audit.append(a)

    kept, audit = _resolve_duplicates(kept, audit, cfg)
    return kept, audit


def _resolve_duplicates(cohorts: List[Cohort], audit: List[Dict[str, str]],
                        cfg: Dict[str, Any]) -> Tuple[List[Cohort], List[Dict[str, str]]]:
    """Keep one row per cohort_id, and apply the overlapping-samples policy."""
    prec = [_norm_key(p) for p in cfg.get("overlap_precedence", [])]

    def rank(c: Cohort) -> Tuple[int, float]:
        t = _norm_key(c.tissue_timing)
        r = prec.index(t) if t in prec else len(prec)
        return (r, c.sei)  # then most precise

    by_id: Dict[str, List[Cohort]] = {}
    for c in cohorts:
        by_id.setdefault(c.cohort_id, []).append(c)
    deduped: List[Cohort] = []
    for cid, group in by_id.items():
        group.sort(key=rank)
        deduped.append(group[0])
        for c in group[1:]:
            audit.append({"file": os.path.basename(c.source_file), "line": str(c.source_row),
                          "cohort_id": cid, "gene": c.gene, "endpoint": c.endpoint,
                          "unit": c.unit, "kept": "no",
                          "reason": f"duplicate cohort_id; kept {group[0].tissue_timing or 'most precise'}"})

    if cfg.get("overlap_policy") == "drop":
        by_group: Dict[str, List[Cohort]] = {}
        singles: List[Cohort] = []
        for c in deduped:
            (by_group.setdefault(c.sample_group, []) if c.sample_group else singles).append(c)
        out = list(singles)
        for gname, group in by_group.items():
            group.sort(key=rank)
            out.append(group[0])
            for c in group[1:]:
                audit.append({"file": os.path.basename(c.source_file), "line": str(c.source_row),
                              "cohort_id": c.cohort_id, "gene": c.gene, "endpoint": c.endpoint,
                              "unit": c.unit, "kept": "no",
                              "reason": f"overlapping samples with {group[0].cohort_id} "
                                        f"(sample_group={gname}); policy=drop"})
        deduped = out

    deduped.sort(key=lambda c: c.cohort_id)
    return deduped, audit


def build_covariance(cohorts: Sequence[Cohort], cfg: Dict[str, Any]) -> np.ndarray:
    """Level-1 covariance V. Diagonal unless overlapping cohorts are kept and correlated."""
    s = np.array([c.sei for c in cohorts], float)
    V = np.diag(s ** 2)
    if cfg.get("overlap_policy") == "correlate":
        rho = float(cfg.get("overlap_rho", 0.0))
        groups = [c.sample_group for c in cohorts]
        for i in range(len(cohorts)):
            for j in range(i + 1, len(cohorts)):
                if groups[i] and groups[i] == groups[j]:
                    V[i, j] = V[j, i] = rho * s[i] * s[j]
    return V


# --------------------------------------------------------------------------------------
# 2. Bayesian normal-normal hierarchical model, exact by quadrature over tau
# --------------------------------------------------------------------------------------

@dataclass
class MixtureNormal:
    """A finite mixture of normals: sum_k w_k N(mean_k, var_k). Weights sum to 1."""
    w: np.ndarray
    mean: np.ndarray
    var: np.ndarray

    def moments(self) -> Tuple[float, float]:
        m = float(np.sum(self.w * self.mean))
        v = float(np.sum(self.w * (self.var + self.mean ** 2)) - m ** 2)
        return m, max(v, 0.0)

    def cdf(self, x: float) -> float:
        return float(np.sum(self.w * stats.norm.cdf(x, loc=self.mean, scale=np.sqrt(self.var))))

    def pdf(self, x: float) -> float:
        return float(np.sum(self.w * stats.norm.pdf(x, loc=self.mean, scale=np.sqrt(self.var))))

    def quantile(self, q: float) -> float:
        m, v = self.moments()
        sd = math.sqrt(v) if v > 0 else 1e-6
        lo, hi = m - 10 * sd, m + 10 * sd
        for _ in range(12):
            if self.cdf(lo) < q < self.cdf(hi):
                break
            lo, hi = m - 2 * (m - lo), m + 2 * (hi - m)
        return float(optimize.brentq(lambda x: self.cdf(x) - q, lo, hi, xtol=1e-12, rtol=1e-14))

    def prob_gt(self, x: float) -> float:
        return 1.0 - self.cdf(x)

    def summary(self, quantiles: Sequence[float]) -> Dict[str, float]:
        m, v = self.moments()
        out = {"mean": m, "sd": math.sqrt(v)}
        for q in quantiles:
            out[f"q{q:g}"] = self.quantile(q)
        return out


def log_tau_prior(tau: np.ndarray, prior: Dict[str, Any]) -> np.ndarray:
    fam = prior.get("tau_family", "half_normal")
    scale = float(prior.get("tau_scale", 0.5) or 0.5)
    t = np.maximum(tau, 0.0)
    with np.errstate(divide="ignore"):
        if fam == "half_normal":
            lp = stats.halfnorm.logpdf(t, scale=scale)
        elif fam == "half_cauchy":
            lp = stats.halfcauchy.logpdf(t, scale=scale)
        elif fam == "exponential":
            lp = stats.expon.logpdf(t, scale=scale)
        elif fam == "uniform":
            lp = np.where(t <= scale, -math.log(scale), -np.inf)
        elif fam == "lognormal_tau2":
            m, s = prior.get("log_tau2_mean"), prior.get("log_tau2_sd")
            if m is None or s is None:
                raise ValueError("lognormal_tau2 prior requires log_tau2_mean and log_tau2_sd; "
                                 "fill these in from the published predictive distribution "
                                 "for the relevant outcome/intervention type.")
            # log(tau^2) ~ N(m, s^2); change of variables to tau adds log|d log tau^2/d tau|
            tt = np.maximum(t, 1e-12)
            lp = stats.norm.logpdf(2 * np.log(tt), loc=float(m), scale=float(s)) + np.log(2.0 / tt)
        else:
            raise ValueError(f"unknown tau_family '{fam}'")
    return lp


@dataclass
class BayesFit:
    tau: np.ndarray
    tau_post: np.ndarray           # normalised density over the tau grid
    mu: MixtureNormal
    theta: List[MixtureNormal]
    pred: MixtureNormal            # effect in a new, exchangeable cohort
    shrinkage: np.ndarray          # posterior mean of s_i^2/(s_i^2+tau^2)
    log_marg: float                # log p(y) under this prior
    prior: Dict[str, Any]
    edge_mass: float               # posterior mass in the top 1% of the tau grid

    def tau_weights(self) -> np.ndarray:
        """Normalised quadrature weights for the tau posterior."""
        w = self.tau_post * _quad_weights(self.tau)
        return w / w.sum()


def fit_bayes(y: np.ndarray, V: np.ndarray, prior: Dict[str, Any],
              grid_cfg: Dict[str, Any]) -> BayesFit:
    """Exact posterior for the normal-normal hierarchical model by 1-D quadrature.

    Model:  y | theta ~ N(theta, V) with V known
            theta_i   ~ N(mu, tau^2)          (exchangeable cohorts)
            mu        ~ N(m0, v0),  tau ~ prior
    Conditional on tau everything is Gaussian, so p(mu | y, tau) and p(theta | y, tau) are
    available in closed form and the marginal posterior of tau is obtained on a grid.
    """
    K = len(y)
    I = np.eye(K)
    one = np.ones(K)
    m0 = float(prior.get("mu_mean", 0.0))
    v0 = float(prior.get("mu_sd", 1.0)) ** 2

    s = np.sqrt(np.diag(V))
    tmax = grid_cfg.get("max")
    if tmax is None:
        spread = float(np.std(y, ddof=1)) if K > 1 else float(np.abs(y).max())
        tmax = max(1.0, 4.0 * spread, 2.0 * float(np.max(s)))
    n_grid = int(grid_cfg.get("n", 801))
    tau = np.linspace(0.0, float(tmax), n_grid)

    log_lik = np.empty(n_grid)
    mu_mean = np.empty(n_grid)
    mu_var = np.empty(n_grid)
    th_mean = np.empty((n_grid, K))
    th_var = np.empty((n_grid, K))

    for k, t in enumerate(tau):
        A = V + (t ** 2) * I
        Ainv = np.linalg.inv(A)
        # marginal of y with mu integrated out: y ~ N(m0 * 1, A + v0 * 11')
        Sig = A + v0 * np.outer(one, one)
        sign, logdet = np.linalg.slogdet(Sig)
        r = y - m0 * one
        sol = np.linalg.solve(Sig, r)
        log_lik[k] = -0.5 * (K * math.log(2 * math.pi) + logdet + float(r @ sol))
        # p(mu | y, tau)
        P = 1.0 / v0 + float(one @ Ainv @ one)
        b = m0 / v0 + float(one @ Ainv @ y)
        mu_mean[k], mu_var[k] = b / P, 1.0 / P
        # p(theta | y, tau, mu) then integrate mu
        if t == 0.0:
            th_mean[k] = mu_mean[k]
            th_var[k] = mu_var[k]
        else:
            C = np.linalg.inv(np.linalg.inv(V) + I / t ** 2)
            g = C @ (np.ones(K) / t ** 2)              # d(theta_hat)/d(mu)
            th_mean[k] = C @ (np.linalg.solve(V, y) + one * mu_mean[k] / t ** 2)
            th_var[k] = np.diag(C) + (g ** 2) * mu_var[k]

    lp = log_lik + log_tau_prior(tau, prior)
    lp_max = float(np.max(lp[np.isfinite(lp)]))
    dens = np.exp(np.where(np.isfinite(lp), lp - lp_max, -np.inf))
    norm_const = float(_trapz(dens, tau))
    if not (norm_const > 0):
        raise RuntimeError("degenerate tau posterior; check inputs and prior")
    tau_post = dens / norm_const

    w = tau_post * _quad_weights(tau)
    w = w / w.sum()
    edge_mass = float(w[int(0.99 * n_grid):].sum())

    mu_mix = MixtureNormal(w, mu_mean, mu_var)
    theta_mix = [MixtureNormal(w, th_mean[:, i], th_var[:, i]) for i in range(K)]
    pred_mix = MixtureNormal(w, mu_mean, mu_var + tau ** 2)
    # Shrinkage weight in the diagonal case: B_i = s_i^2 / (s_i^2 + tau^2), averaged over
    # the tau posterior. B_i = 0 means no pooling, B_i = 1 means the cohort is replaced by
    # the grand mean.
    shrink = np.array([float(np.sum(w * (s[i] ** 2 / (s[i] ** 2 + tau ** 2)))) for i in range(K)])
    log_marg = float(math.log(norm_const) + lp_max)

    return BayesFit(tau=tau, tau_post=tau_post, mu=mu_mix, theta=theta_mix, pred=pred_mix,
                    shrinkage=shrink, log_marg=log_marg, prior=prior, edge_mass=edge_mass)


def savage_dickey_bf01(fit: BayesFit) -> float:
    """BF for mu = 0 vs mu != 0 (values > 1 favour the null). Requires a proper mu prior."""
    prior_sd = float(fit.prior.get("mu_sd", 1.0))
    if prior_sd > 50:
        return float("nan")  # not defined against an effectively flat prior
    prior_dens = stats.norm.pdf(0.0, loc=float(fit.prior.get("mu_mean", 0.0)), scale=prior_sd)
    post_dens = fit.mu.pdf(0.0)
    return post_dens / prior_dens if prior_dens > 0 else float("nan")


# --------------------------------------------------------------------------------------
# 3. Frequentist companions
# --------------------------------------------------------------------------------------

def _fe(y: np.ndarray, v: np.ndarray) -> Tuple[float, float]:
    w = 1.0 / v
    m = float(np.sum(w * y) / np.sum(w))
    return m, float(math.sqrt(1.0 / np.sum(w)))


def cochran_q(y: np.ndarray, v: np.ndarray) -> Tuple[float, int, float]:
    m, _ = _fe(y, v)
    Q = float(np.sum((y - m) ** 2 / v))
    df = len(y) - 1
    p = float(stats.chi2.sf(Q, df)) if df > 0 else float("nan")
    return Q, df, p


def tau2_dl(y: np.ndarray, v: np.ndarray) -> float:
    w = 1.0 / v
    Q, df, _ = cochran_q(y, v)
    c = float(np.sum(w) - np.sum(w ** 2) / np.sum(w))
    return max(0.0, (Q - df) / c) if c > 0 else 0.0


def tau2_pm(y: np.ndarray, v: np.ndarray) -> float:
    """Paule-Mandel: solve sum w_i(tau^2)(y_i - mu_w)^2 = K - 1."""
    df = len(y) - 1
    if df <= 0:
        return 0.0

    def gen_q(t2: float) -> float:
        w = 1.0 / (v + t2)
        m = np.sum(w * y) / np.sum(w)
        return float(np.sum(w * (y - m) ** 2) - df)

    if gen_q(0.0) <= 0:
        return 0.0
    hi = 1.0
    for _ in range(80):
        if gen_q(hi) < 0:
            break
        hi *= 2.0
    else:
        return hi
    return float(optimize.brentq(gen_q, 0.0, hi, xtol=1e-12))


def tau2_reml(y: np.ndarray, v: np.ndarray) -> float:
    """REML for the intercept-only random-effects model."""
    def neg_ll(t2: float) -> float:
        t2 = max(t2, 0.0)
        w = 1.0 / (v + t2)
        sw = np.sum(w)
        m = np.sum(w * y) / sw
        return float(0.5 * (np.sum(np.log(v + t2)) + np.sum(w * (y - m) ** 2) + math.log(sw)))

    res = optimize.minimize_scalar(neg_ll, bounds=(0.0, max(1.0, 10 * float(np.var(y, ddof=1) if len(y) > 1 else 1.0))),
                                   method="bounded", options={"xatol": 1e-10})
    return float(max(0.0, res.x))


def random_effects(y: np.ndarray, v: np.ndarray, tau2: float,
                   knapp_hartung: bool = False) -> Dict[str, float]:
    K = len(y)
    w = 1.0 / (v + tau2)
    m = float(np.sum(w * y) / np.sum(w))
    se = float(math.sqrt(1.0 / np.sum(w)))
    if knapp_hartung and K > 1:
        # Hartung-Knapp-Sidik-Jonkman: t reference distribution and a variance estimate
        # driven by the observed scatter. Recommended default for small K; the ad hoc
        # truncation keeps it from being anti-conservative.
        se_hk = float(math.sqrt(np.sum(w * (y - m) ** 2) / ((K - 1) * np.sum(w))))
        se_used = max(se_hk, se)
        crit = float(stats.t.ppf(0.975, K - 1))
    else:
        se_used, crit = se, float(stats.norm.ppf(0.975))
    lo, hi = m - crit * se_used, m + crit * se_used
    z = m / se_used
    pval = 2 * (1 - stats.t.cdf(abs(z), K - 1)) if knapp_hartung else 2 * (1 - stats.norm.cdf(abs(z)))
    out = {"estimate": m, "se": se_used, "ci_low": lo, "ci_high": hi, "pvalue": float(pval),
           "tau2": tau2, "tau": math.sqrt(tau2)}
    if K > 2 and tau2 > 0:
        # Higgins-Thompson-Spiegelhalter prediction interval
        crit_p = float(stats.t.ppf(0.975, K - 2))
        half = crit_p * math.sqrt(tau2 + se ** 2)
        out["pi_low"], out["pi_high"] = m - half, m + half
    else:
        out["pi_low"] = out["pi_high"] = float("nan")
    return out


def heterogeneity(y: np.ndarray, v: np.ndarray) -> Dict[str, float]:
    Q, df, p = cochran_q(y, v)
    I2 = max(0.0, (Q - df) / Q) * 100 if Q > 0 else 0.0
    H2 = Q / df if df > 0 else float("nan")
    return {"Q": Q, "df": float(df), "Q_pvalue": p, "I2_percent": I2, "H2": H2,
            "tau2_DL": tau2_dl(y, v), "tau2_PM": tau2_pm(y, v), "tau2_REML": tau2_reml(y, v)}


# --------------------------------------------------------------------------------------
# 4. Small-study effects and selection
# --------------------------------------------------------------------------------------

def egger_test(y: np.ndarray, s: np.ndarray, ni: Optional[np.ndarray] = None,
               predictor: str = "sei") -> Dict[str, float]:
    """Regression test for funnel asymmetry.

    predictor="sei"   classic Egger (radial regression, intercept = asymmetry).
    predictor="inv_n" sample-size based variant, preferable for log-OR where the
                      estimate and its SE are structurally correlated.
    """
    K = len(y)
    if K < 3:
        return {"intercept": float("nan"), "se": float("nan"), "t": float("nan"),
                "df": float(max(K - 2, 0)), "pvalue": float("nan"), "predictor": predictor}
    if predictor == "inv_n":
        if ni is None or np.any(~np.isfinite(ni)) or np.any(ni <= 0):
            return {"intercept": float("nan"), "se": float("nan"), "t": float("nan"),
                    "df": float(K - 2), "pvalue": float("nan"), "predictor": predictor}
        x = 1.0 / ni
        w = 1.0 / s ** 2
        X = np.column_stack([np.ones(K), x])
        W = np.diag(w)
        beta = np.linalg.solve(X.T @ W @ X, X.T @ W @ y)
        resid = y - X @ beta
        dof = K - 2
        sigma2 = float(resid @ W @ resid) / dof
        cov = sigma2 * np.linalg.inv(X.T @ W @ X)
        b0, se0 = float(beta[0]), float(math.sqrt(cov[0, 0]))
        # Here the intercept is the bias-free effect and the slope carries the asymmetry.
        b1, se1 = float(beta[1]), float(math.sqrt(cov[1, 1]))
        tstat = b1 / se1
        return {"intercept": b0, "se": se1, "t": tstat, "df": float(dof),
                "pvalue": float(2 * (1 - stats.t.cdf(abs(tstat), dof))),
                "predictor": predictor, "slope": b1}
    z, xx = y / s, 1.0 / s
    X = np.column_stack([np.ones(K), xx])
    beta, *_ = np.linalg.lstsq(X, z, rcond=None)
    resid = z - X @ beta
    dof = K - 2
    sigma2 = float(resid @ resid) / dof
    cov = sigma2 * np.linalg.inv(X.T @ X)
    b0, se0 = float(beta[0]), float(math.sqrt(cov[0, 0]))
    tstat = b0 / se0
    return {"intercept": b0, "se": se0, "t": tstat, "df": float(dof),
            "pvalue": float(2 * (1 - stats.t.cdf(abs(tstat), dof))),
            "predictor": predictor, "slope": float(beta[1])}


def begg_test(y: np.ndarray, v: np.ndarray) -> Dict[str, float]:
    """Begg & Mazumdar rank correlation between standardised effect and variance."""
    K = len(y)
    if K < 3:
        return {"kendall_tau": float("nan"), "pvalue": float("nan")}
    m, _ = _fe(y, v)
    v_star = v - 1.0 / np.sum(1.0 / v)
    v_star = np.where(v_star > 0, v_star, np.nan)
    y_star = (y - m) / np.sqrt(v_star)
    ok = np.isfinite(y_star)
    if ok.sum() < 3:
        return {"kendall_tau": float("nan"), "pvalue": float("nan")}
    tau_k, p = stats.kendalltau(y_star[ok], v[ok])
    return {"kendall_tau": float(tau_k), "pvalue": float(p)}


def trim_and_fill(y: np.ndarray, v: np.ndarray, estimator: str = "L0",
                  side: Optional[str] = None, max_iter: int = 50) -> Dict[str, Any]:
    """Duval & Tweedie trim-and-fill. Reported as a *what-if*, never as a corrected truth."""
    K = len(y)
    if K < 3:
        return {"k0": 0, "side": side or "none", "estimate": float("nan"),
                "se": float("nan"), "ci_low": float("nan"), "ci_high": float("nan")}
    if side is None:
        eg = egger_test(y, np.sqrt(v))
        side = "left" if eg["intercept"] > 0 else "right"
    flip = 1.0 if side == "left" else -1.0
    yy, vv = flip * y, v

    def k0_estimate(yc: np.ndarray, center: float) -> int:
        d = yc - center
        n = len(d)
        order = np.argsort(np.abs(d), kind="mergesort")
        ranks = np.empty(n)
        ranks[order] = np.arange(1, n + 1)
        signs = np.sign(d)
        if estimator == "R0":
            # k0 = gamma* - 1, with gamma* the length of the rightmost run (largest
            # absolute deviations) of same-signed values on the over-represented side.
            srt = signs[order]
            run = 0
            for sgn in srt[::-1]:
                if sgn > 0:
                    run += 1
                else:
                    break
            return max(0, run - 1)
        Tn = float(np.sum(ranks[signs > 0]))
        return int(max(0, round((4 * Tn - n * (n + 1)) / (2 * n - 1))))

    # Iterate: trim the k0 most extreme studies on the over-represented side, re-centre,
    # re-estimate k0 from the full set relative to the trimmed centre, until k0 is stable.
    k0 = 0
    keep = np.ones(len(yy), bool)
    for _ in range(max_iter):
        keep = np.ones(len(yy), bool)
        if k0 > 0:
            m_tmp, _ = _fe(yy, vv)
            keep[np.argsort(yy - m_tmp)[-k0:]] = False   # trim the largest positive deviations
        m_c, _ = _fe(yy[keep], vv[keep])
        k_new = k0_estimate(yy, m_c)
        if k_new == k0:
            break
        k0 = k_new
    m_c, _ = _fe(yy[keep], vv[keep])
    if k0 > 0:
        order = np.argsort(yy - m_c)
        idx = order[-k0:]                            # mirror the k0 most extreme positives
        y_fill = 2 * m_c - yy[idx]
        y_aug = np.concatenate([yy, y_fill])
        v_aug = np.concatenate([vv, vv[idx]])
    else:
        y_aug, v_aug = yy, vv
    t2 = tau2_dl(y_aug, v_aug)
    res = random_effects(y_aug, v_aug, t2)
    return {"k0": int(k0), "side": side, "estimate": flip * res["estimate"],
            "se": res["se"], "ci_low": flip * res["ci_high"] if flip < 0 else res["ci_low"],
            "ci_high": flip * res["ci_low"] if flip < 0 else res["ci_high"],
            "tau2": t2}


def pet_peese(y: np.ndarray, s: np.ndarray) -> Dict[str, float]:
    """PET and PEESE WLS regressions; the conditional rule is reported, not enforced."""
    K = len(y)
    out: Dict[str, float] = {}
    if K < 3:
        return {"pet_estimate": float("nan"), "pet_pvalue": float("nan"),
                "peese_estimate": float("nan"), "selected": "insufficient_K"}
    w = 1.0 / s ** 2
    for name, x in (("pet", s), ("peese", s ** 2)):
        X = np.column_stack([np.ones(K), x])
        W = np.diag(w)
        XtWX = X.T @ W @ X
        beta = np.linalg.solve(XtWX, X.T @ W @ y)
        resid = y - X @ beta
        dof = K - 2
        sigma2 = float(resid @ W @ resid) / dof
        cov = sigma2 * np.linalg.inv(XtWX)
        out[f"{name}_estimate"] = float(beta[0])
        out[f"{name}_se"] = float(math.sqrt(cov[0, 0]))
        tstat = float(beta[0] / math.sqrt(cov[0, 0]))
        out[f"{name}_t"] = tstat
        out[f"{name}_pvalue"] = float(2 * (1 - stats.t.cdf(abs(tstat), dof)))
    # Conditional rule: if PET rejects a null bias-adjusted effect (one-sided), use PEESE.
    one_sided = out["pet_pvalue"] / 2 if out["pet_estimate"] > 0 else 1 - out["pet_pvalue"] / 2
    out["pet_one_sided_p"] = float(one_sided)
    out["selected"] = "peese" if one_sided < 0.05 else "pet"
    out["selected_estimate"] = out["peese_estimate"] if out["selected"] == "peese" else out["pet_estimate"]
    return out


def selection_model_step(y: np.ndarray, s: np.ndarray, alpha: float = 0.025,
                         omega_fixed: Optional[float] = None) -> Dict[str, float]:
    """One-step weight-function selection model (Iyengar-Greenhouse / Vevea-Hedges).

    Studies with a one-sided p below `alpha` in the favoured direction have weight 1;
    the rest have weight omega. omega < 1 means non-significant cohorts are less likely
    to be observed. With K small this is weakly identified: report it as a sensitivity
    curve over fixed omega, and treat the free-omega MLE as descriptive.
    """
    K = len(y)
    crit = float(stats.norm.ppf(1 - alpha))
    favoured = 1.0 if np.sum(y / s) >= 0 else -1.0
    yy = favoured * y
    c = crit * s                                     # threshold for "significant"

    def neg_ll(par: np.ndarray) -> float:
        mu, log_tau = par[0], par[1]
        omega = omega_fixed if omega_fixed is not None else math.exp(par[2])
        tau2 = math.exp(2 * log_tau)
        sd = np.sqrt(s ** 2 + tau2)
        p_sig = 1 - stats.norm.cdf((c - mu) / sd)
        A = p_sig + omega * (1 - p_sig)
        wts = np.where(yy > c, 1.0, omega)
        if np.any(A <= 0) or omega <= 0:
            return 1e12
        ll = np.sum(np.log(wts) + stats.norm.logpdf(yy, loc=mu, scale=sd) - np.log(A))
        return -ll if math.isfinite(ll) else 1e12

    x0 = [float(np.mean(yy)), math.log(max(1e-3, float(np.std(yy, ddof=1)) if K > 1 else 0.1))]
    if omega_fixed is None:
        x0 = x0 + [0.0]
    res = optimize.minimize(neg_ll, x0=np.array(x0), method="Nelder-Mead",
                            options={"maxiter": 5000, "xatol": 1e-8, "fatol": 1e-10})
    mu_hat = favoured * float(res.x[0])
    tau_hat = float(math.exp(res.x[1]))
    omega_hat = omega_fixed if omega_fixed is not None else float(math.exp(res.x[2]))
    # LRT against omega = 1 (no selection)
    ll_sel = -float(res.fun)
    res1 = optimize.minimize(lambda p: selection_model_step_ll(p, yy, s, c, 1.0),
                             x0=np.array(x0[:2]), method="Nelder-Mead",
                             options={"maxiter": 5000})
    ll_null = -float(res1.fun)
    lrt = 2 * (ll_sel - ll_null)
    return {"estimate": mu_hat, "tau": tau_hat, "omega": omega_hat,
            "loglik": ll_sel, "lrt_vs_no_selection": float(lrt),
            "lrt_pvalue": float(stats.chi2.sf(max(lrt, 0.0), 1)) if omega_fixed is None else float("nan"),
            "converged": bool(res.success)}


def selection_model_step_ll(par: np.ndarray, yy: np.ndarray, s: np.ndarray,
                            c: np.ndarray, omega: float) -> float:
    mu, log_tau = par[0], par[1]
    tau2 = math.exp(2 * log_tau)
    sd = np.sqrt(s ** 2 + tau2)
    p_sig = 1 - stats.norm.cdf((c - mu) / sd)
    A = p_sig + omega * (1 - p_sig)
    wts = np.where(yy > c, 1.0, omega)
    if np.any(A <= 0):
        return 1e12
    ll = np.sum(np.log(wts) + stats.norm.logpdf(yy, loc=mu, scale=sd) - np.log(A))
    return -ll if math.isfinite(ll) else 1e12


# --------------------------------------------------------------------------------------
# 5. Missing cohorts (EGA / dbGaP data-access-controlled trials)
# --------------------------------------------------------------------------------------

@dataclass
class MissingCohort:
    cohort_id: str
    status: str = "missing"
    ni: Optional[float] = None
    events: Optional[float] = None
    response_rate: Optional[float] = None
    sei: Optional[float] = None
    notes: str = ""


def read_registry(path: str) -> Tuple[List[MissingCohort], List[Dict[str, str]]]:
    """Read the cohort registry; rows whose status is not `included` count as unobserved."""
    missing: List[MissingCohort] = []
    all_rows: List[Dict[str, str]] = []
    with open(path, newline="", encoding="utf-8-sig") as fh:
        for row in csv.DictReader(fh):
            clean = { _norm_key(k): v for k, v in row.items() if k }
            all_rows.append({k: _clean_str(v) for k, v in row.items() if k})
            status = _clean_str(clean.get("status")).lower()
            if status in {"included", "include", "analysed", "analyzed"}:
                continue
            missing.append(MissingCohort(
                cohort_id=_clean_str(clean.get("cohortid") or clean.get("accession") or "?"),
                status=status or "unknown",
                ni=_to_float(clean.get("n") or clean.get("ntotal") or clean.get("samplesize")),
                events=_to_float(clean.get("events") or clean.get("nevents")),
                response_rate=_to_float(clean.get("responserate") or clean.get("orr")),
                sei=_to_float(clean.get("sei") or clean.get("se")),
                notes=_clean_str(clean.get("notes")),
            ))
    return missing, all_rows


def impute_se(mc: MissingCohort, endpoint: str, median_se: float, rule: str = "auto") -> Tuple[float, str]:
    """Plausible SE for an unobserved cohort.

    For a Cox model with a single standardised covariate, Var(log HR) ~= 1/E where E is
    the number of events. For logistic regression with a standardised covariate,
    Var(log OR) ~= 1/(n * p * (1-p)). Both ignore covariate adjustment and assume a small
    true effect, which is adequate for a sensitivity analysis but not for inference.
    """
    if mc.sei and mc.sei > 0:
        return mc.sei, "registry_se"
    if rule == "auto":
        ep = endpoint.upper()
        if ep in SURVIVAL_ENDPOINTS and mc.events and mc.events > 0:
            return 1.0 / math.sqrt(mc.events), "1/sqrt(events)"
        if ep in RESPONSE_ENDPOINTS and mc.ni and mc.ni > 0:
            p = mc.response_rate if (mc.response_rate and 0 < mc.response_rate < 1) else 0.25
            return 1.0 / math.sqrt(mc.ni * p * (1 - p)), "1/sqrt(n*p*(1-p))"
        if mc.ni and mc.ni > 0:
            # fall back on a 40% event fraction, a common order of magnitude in ICI trials
            return 1.0 / math.sqrt(0.4 * mc.ni), "1/sqrt(0.4n)"
    return median_se, "median_observed_se"


def delta_sweep(y: np.ndarray, V: np.ndarray, se_missing: Sequence[float],
                prior: Dict[str, Any], grid_cfg: Dict[str, Any],
                deltas: Sequence[float], quantiles: Sequence[float]) -> List[Dict[str, float]]:
    """Pattern-mixture MNAR sensitivity.

    Unobserved cohorts are assumed to have effects centred at mu + delta. delta = 0 is
    MAR: unobserved cohorts then carry no information about mu beyond precision, so the
    pooled estimate is unchanged and only the interval narrows. Everything interesting is
    at delta != 0.
    """
    out: List[Dict[str, float]] = []
    K = len(y)
    m_obs = float(np.sum(y / np.diag(V)) / np.sum(1.0 / np.diag(V)))
    for d in deltas:
        y_aug = np.concatenate([y, np.full(len(se_missing), m_obs + d)])
        V_aug = np.zeros((K + len(se_missing), K + len(se_missing)))
        V_aug[:K, :K] = V
        for j, sm in enumerate(se_missing):
            V_aug[K + j, K + j] = sm ** 2
        fit = fit_bayes(y_aug, V_aug, prior, grid_cfg)
        s = fit.mu.summary(quantiles)
        out.append({"delta": float(d), "n_missing": len(se_missing),
                    "mu_mean": s["mean"], "mu_sd": s["sd"],
                    "ci_low": s[f"q{quantiles[0]:g}"], "ci_high": s[f"q{quantiles[-1]:g}"],
                    "prob_gt_0": fit.mu.prob_gt(0.0),
                    "tau_mean": float(np.sum(fit.tau * fit.tau_weights()))})
    return out


def manski_bounds(y: np.ndarray, V: np.ndarray, se_missing: Sequence[float],
                  prior: Dict[str, Any], grid_cfg: Dict[str, Any],
                  bound_low: float, bound_high: float,
                  quantiles: Sequence[float]) -> Dict[str, float]:
    """Worst-case envelope: all unobserved cohorts at the low, then the high, bound."""
    K = len(y)
    res: Dict[str, float] = {"bound_low_value": bound_low, "bound_high_value": bound_high}
    for tag, val in (("low", bound_low), ("high", bound_high)):
        y_aug = np.concatenate([y, np.full(len(se_missing), val)])
        V_aug = np.zeros((K + len(se_missing), K + len(se_missing)))
        V_aug[:K, :K] = V
        for j, sm in enumerate(se_missing):
            V_aug[K + j, K + j] = sm ** 2
        fit = fit_bayes(y_aug, V_aug, prior, grid_cfg)
        s = fit.mu.summary(quantiles)
        res[f"mu_{tag}"] = s["mean"]
        res[f"ci_low_{tag}"] = s[f"q{quantiles[0]:g}"]
        res[f"ci_high_{tag}"] = s[f"q{quantiles[-1]:g}"]
        res[f"prob_gt_0_{tag}"] = fit.mu.prob_gt(0.0)
    res["envelope_low"] = min(res["ci_low_low"], res["ci_low_high"])
    res["envelope_high"] = max(res["ci_high_low"], res["ci_high_high"])
    return res


def failsafe_k(y: np.ndarray, V: np.ndarray, prior: Dict[str, Any], grid_cfg: Dict[str, Any],
               threshold: float = 0.95, se_null: Optional[float] = None,
               max_k: int = 200) -> Dict[str, float]:
    """Smallest number of unobserved *null* cohorts of typical precision that would push
    the posterior probability of a non-zero effect below `threshold`.

    This is the Bayesian replacement for Rosenthal's fail-safe N, which is discredited
    because it counts studies needed to reach p > 0.05 under an implausible exactly-null
    file drawer and ignores precision.
    """
    K = len(y)
    s_obs = np.sqrt(np.diag(V))
    se_null = float(np.median(s_obs)) if se_null is None else se_null
    base = fit_bayes(y, V, prior, grid_cfg)
    p0 = base.mu.prob_gt(0.0)
    direction = 1.0 if p0 >= 0.5 else -1.0
    prob = p0 if direction > 0 else 1 - p0
    if prob < threshold:
        return {"failsafe_k": 0.0, "baseline_prob": p0, "se_null": se_null,
                "note": "baseline already below threshold"}
    for k in range(1, max_k + 1):
        y_aug = np.concatenate([y, np.zeros(k)])
        V_aug = np.zeros((K + k, K + k))
        V_aug[:K, :K] = V
        V_aug[K:, K:] = np.eye(k) * se_null ** 2
        fit = fit_bayes(y_aug, V_aug, prior, grid_cfg)
        pk = fit.mu.prob_gt(0.0)
        pk = pk if direction > 0 else 1 - pk
        if pk < threshold:
            return {"failsafe_k": float(k), "baseline_prob": p0, "se_null": se_null,
                    "prob_at_k": float(pk), "note": ""}
    return {"failsafe_k": float("inf"), "baseline_prob": p0, "se_null": se_null,
            "note": f"not reached within max_k={max_k}"}


# --------------------------------------------------------------------------------------
# 6. Meta-regression (exploratory) and leave-one-out
# --------------------------------------------------------------------------------------

def meta_regression(y: np.ndarray, v: np.ndarray, X: np.ndarray,
                    names: Sequence[str]) -> Dict[str, Any]:
    """Random-effects meta-regression by WLS with a moment-based residual tau^2."""
    K, p = X.shape
    if K <= p:
        return {"error": f"K={K} cohorts cannot support {p} coefficients"}

    def resid_tau2(t2: float) -> float:
        W = np.diag(1.0 / (v + t2))
        H = X @ np.linalg.solve(X.T @ W @ X, X.T @ W)
        r = y - H @ y
        return float(r @ W @ r) - (K - p)

    t2 = 0.0
    if resid_tau2(0.0) > 0:
        hi = 1.0
        for _ in range(60):
            if resid_tau2(hi) < 0:
                break
            hi *= 2
        t2 = float(optimize.brentq(resid_tau2, 0.0, hi, xtol=1e-12))
    W = np.diag(1.0 / (v + t2))
    XtWX = X.T @ W @ X
    beta = np.linalg.solve(XtWX, X.T @ W @ y)
    cov = np.linalg.inv(XtWX)
    se = np.sqrt(np.diag(cov))
    z = beta / se
    return {"tau2_resid": t2,
            "coefficients": [
                {"term": names[j], "estimate": float(beta[j]), "se": float(se[j]),
                 "z": float(z[j]), "pvalue": float(2 * (1 - stats.norm.cdf(abs(z[j])))),
                 "ci_low": float(beta[j] - 1.96 * se[j]),
                 "ci_high": float(beta[j] + 1.96 * se[j])}
                for j in range(p)],
            "beta": beta, "cov": cov}


def leave_one_out(cohorts: Sequence[Cohort], V: np.ndarray, prior: Dict[str, Any],
                  grid_cfg: Dict[str, Any], quantiles: Sequence[float]) -> List[Dict[str, Any]]:
    y = np.array([c.yi for c in cohorts])
    out = []
    for i in range(len(cohorts)):
        idx = [j for j in range(len(cohorts)) if j != i]
        fit = fit_bayes(y[idx], V[np.ix_(idx, idx)], prior, grid_cfg)
        s = fit.mu.summary(quantiles)
        out.append({"omitted": cohorts[i].cohort_id, "mu_mean": s["mean"], "mu_sd": s["sd"],
                    "ci_low": s[f"q{quantiles[0]:g}"], "ci_high": s[f"q{quantiles[-1]:g}"],
                    "prob_gt_0": fit.mu.prob_gt(0.0)})
    return out


# --------------------------------------------------------------------------------------
# 7. Output helpers
# --------------------------------------------------------------------------------------

def write_csv(path: str, rows: Sequence[Dict[str, Any]], fieldnames: Optional[Sequence[str]] = None) -> None:
    if not rows:
        rows = []
    if fieldnames is None:
        seen: List[str] = []
        for r in rows:
            for k in r:
                if k not in seen:
                    seen.append(k)
        fieldnames = seen or ["empty"]
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(fieldnames), extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: _fmt(r.get(k)) for k in fieldnames})


def _fmt(v: Any) -> Any:
    if isinstance(v, float):
        if not math.isfinite(v):
            return "NA" if math.isnan(v) else ("Inf" if v > 0 else "-Inf")
        return f"{v:.10g}"
    return v


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def maybe_plot(outdir: str, cohorts: Sequence[Cohort], fit: BayesFit,
               deltas: Sequence[Dict[str, float]]) -> List[str]:
    """Figures are a convenience; every underlying number is already in the CSVs."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return []
    made = []
    y = np.array([c.yi for c in cohorts])
    s = np.array([c.sei for c in cohorts])
    labels = [c.cohort_id for c in cohorts]
    order = np.argsort(y)

    fig, ax = plt.subplots(figsize=(7, 0.38 * len(cohorts) + 2.6))
    pos = np.arange(len(cohorts))
    for rank, i in enumerate(order):
        ax.errorbar(y[i], rank + 0.16, xerr=1.96 * s[i], fmt="o", color="0.35",
                    ms=4, lw=1, capsize=2)
        th = fit.theta[i].summary([0.025, 0.975])
        ax.errorbar(th["mean"], rank - 0.16,
                    xerr=[[th["mean"] - th["q0.025"]], [th["q0.975"] - th["mean"]]],
                    fmt="s", color="#1f77b4", ms=4, lw=1.4, capsize=2)
    mu = fit.mu.summary([0.025, 0.975])
    ax.axvline(0, color="0.7", lw=0.8)
    ax.axvspan(mu["q0.025"], mu["q0.975"], color="#1f77b4", alpha=0.12)
    ax.axvline(mu["mean"], color="#1f77b4", lw=1.2)
    ax.set_yticks(pos)
    ax.set_yticklabels([labels[i] for i in order], fontsize=8)
    ax.set_xlabel("effect (log scale, per SD)  --  grey: observed, blue: partially pooled")
    ax.set_title("Forest / shrinkage")
    fig.tight_layout()
    p = os.path.join(outdir, "fig_forest_shrinkage.png")
    fig.savefig(p, dpi=150); plt.close(fig); made.append(p)

    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    ax.scatter(y, s, s=22, color="0.3")
    ax.invert_yaxis()
    ax.axvline(mu["mean"], color="#1f77b4", lw=1)
    smax = float(np.max(s)) * 1.05
    for z, alpha in ((1.96, 0.10), (2.58, 0.06)):
        ax.fill_betweenx([0, smax], [mu["mean"], mu["mean"] - z * smax],
                         [mu["mean"], mu["mean"] + z * smax], color="#1f77b4", alpha=alpha)
    ax.set_xlabel("effect"); ax.set_ylabel("standard error"); ax.set_title("Funnel")
    fig.tight_layout()
    p = os.path.join(outdir, "fig_funnel.png")
    fig.savefig(p, dpi=150); plt.close(fig); made.append(p)

    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    ax.plot(fit.tau, fit.tau_post, color="#1f77b4")
    ax.set_xlabel("tau (between-cohort SD)"); ax.set_ylabel("posterior density")
    ax.set_title("Heterogeneity")
    fig.tight_layout()
    p = os.path.join(outdir, "fig_tau_posterior.png")
    fig.savefig(p, dpi=150); plt.close(fig); made.append(p)

    if deltas:
        d = np.array([r["delta"] for r in deltas])
        m = np.array([r["mu_mean"] for r in deltas])
        lo = np.array([r["ci_low"] for r in deltas])
        hi = np.array([r["ci_high"] for r in deltas])
        fig, ax = plt.subplots(figsize=(5.6, 3.8))
        ax.fill_between(d, lo, hi, color="#d62728", alpha=0.15)
        ax.plot(d, m, color="#d62728")
        ax.axhline(0, color="0.6", lw=0.8)
        ax.set_xlabel("delta: shift of unobserved (EGA/dbGaP) cohorts, log scale")
        ax.set_ylabel("pooled effect")
        ax.set_title("MNAR sensitivity")
        fig.tight_layout()
        p = os.path.join(outdir, "fig_delta_sensitivity.png")
        fig.savefig(p, dpi=150); plt.close(fig); made.append(p)
    return made


# --------------------------------------------------------------------------------------
# 8. Driver
# --------------------------------------------------------------------------------------

def deep_update(base: Dict[str, Any], other: Dict[str, Any]) -> Dict[str, Any]:
    out = copy.deepcopy(base)
    for k, v in other.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = deep_update(out[k], v)
        else:
            out[k] = v
    return out


def find_result_files(results_dir: str, recursive: bool) -> List[str]:
    pattern = os.path.join(results_dir, "**", "*.csv") if recursive else os.path.join(results_dir, "*.csv")
    return sorted(p for p in glob.glob(pattern, recursive=recursive) if os.path.isfile(p))


def run(cfg: Dict[str, Any], files: Sequence[str], outdir: str) -> int:
    os.makedirs(outdir, exist_ok=True)
    quantiles = cfg["quantiles"]
    q_lo, q_hi = quantiles[0], quantiles[-1]

    rows = read_result_files(files)
    cohorts, audit = harmonise(rows, cfg)
    write_csv(os.path.join(outdir, "input_audit.csv"), audit,
              ["file", "line", "cohort_id", "gene", "endpoint", "unit", "access", "kept", "reason"])

    K = len(cohorts)
    print(f"[hier_meta] {len(rows)} input rows from {len(files)} file(s); "
          f"{K} cohort(s) retained for gene={cfg['gene']} endpoint={cfg['primary_endpoint']}")
    if K == 0:
        print(REQUIRED_MESSAGE)
        return 0
    if K < 3:
        print("[hier_meta] WARNING: fewer than 3 cohorts. Between-cohort variance is not "
              "estimable in any useful sense; the tau prior will determine the answer. "
              "Report the cohorts individually instead of pooling.")

    y = np.array([c.yi for c in cohorts], float)
    s = np.array([c.sei for c in cohorts], float)
    v = s ** 2
    V = build_covariance(cohorts, cfg)
    ni = np.array([c.ni if c.ni else np.nan for c in cohorts], float)

    # --- Bayesian primary ------------------------------------------------------------
    fit = fit_bayes(y, V, cfg["prior"], cfg["tau_grid"])
    if fit.edge_mass > 0.01:
        print(f"[hier_meta] WARNING: {fit.edge_mass:.1%} of the tau posterior sits at the top "
              f"of the grid; increase tau_grid.max.")
    mu_s = fit.mu.summary(quantiles)
    pred_s = fit.pred.summary(quantiles)
    tau_w = fit.tau_weights()
    tau_mean = float(np.sum(fit.tau * tau_w))
    tau_med = float(fit.tau[np.searchsorted(np.cumsum(tau_w), 0.5)])
    tau_q975 = float(fit.tau[np.searchsorted(np.cumsum(tau_w), 0.975)])
    rope = float(cfg["rope"])
    is_smd = cfg.get("effect_family", "log_ratio") == "smd"

    def _ratio(x: float) -> float:
        return float("nan") if is_smd or not math.isfinite(x) else math.exp(x)

    pooled: List[Dict[str, Any]] = [{
        "method": "bayes_hierarchical",
        "prior": cfg["prior"].get("tau_family") + f"({cfg['prior'].get('tau_scale')})",
        "k": K, "estimate": mu_s["mean"], "se": mu_s["sd"],
        "ci_low": mu_s[f"q{q_lo:g}"], "ci_high": mu_s[f"q{q_hi:g}"],
        "ratio": _ratio(mu_s["mean"]),
        "ratio_low": _ratio(mu_s[f"q{q_lo:g}"]), "ratio_high": _ratio(mu_s[f"q{q_hi:g}"]),
        "effect_family": cfg.get("effect_family", "log_ratio"),
        "access_policy": cfg.get("access_policy", "open_only"),
        "pi_low": pred_s[f"q{q_lo:g}"], "pi_high": pred_s[f"q{q_hi:g}"],
        "tau_mean": tau_mean, "tau_median": tau_med, "tau_q97.5": tau_q975,
        "prob_gt_0": fit.mu.prob_gt(0.0),
        "prob_outside_rope": fit.mu.prob_gt(rope) + fit.mu.cdf(-rope),
        "prob_in_rope": 1 - (fit.mu.prob_gt(rope) + fit.mu.cdf(-rope)),
        "bf01_savage_dickey": savage_dickey_bf01(fit),
        "note": "posterior mean; CrI = credible interval; PI = new-cohort prediction interval",
    }]

    # --- Frequentist companions ------------------------------------------------------
    fe_m, fe_se = _fe(y, v)
    zc = stats.norm.ppf(0.975)
    pooled.append({"method": "fixed_effect", "k": K, "estimate": fe_m, "se": fe_se,
                   "ci_low": fe_m - zc * fe_se, "ci_high": fe_m + zc * fe_se,
                   "ratio": _ratio(fe_m), "ratio_low": _ratio(fe_m - zc * fe_se),
                   "ratio_high": _ratio(fe_m + zc * fe_se),
                   "note": "common-effect; assumes no between-cohort heterogeneity"})
    het = heterogeneity(y, v)
    for label, t2, hk in (("random_DL", het["tau2_DL"], False),
                          ("random_PM", het["tau2_PM"], False),
                          ("random_REML", het["tau2_REML"], False),
                          ("random_REML_HKSJ", het["tau2_REML"], True)):
        r = random_effects(y, v, t2, knapp_hartung=hk)
        pooled.append({"method": label, "k": K, "estimate": r["estimate"], "se": r["se"],
                       "ci_low": r["ci_low"], "ci_high": r["ci_high"],
                       "ratio": _ratio(r["estimate"]),
                       "ratio_low": _ratio(r["ci_low"]), "ratio_high": _ratio(r["ci_high"]),
                       "pi_low": r["pi_low"], "pi_high": r["pi_high"],
                       "tau_mean": math.sqrt(t2), "pvalue": r["pvalue"],
                       "note": "Hartung-Knapp: t reference, recommended for small K" if hk else ""})

    # --- Shrinkage per cohort --------------------------------------------------------
    shrink_rows = []
    w_re = 1.0 / (v + het["tau2_REML"])
    for i, c in enumerate(cohorts):
        th = fit.theta[i].summary(quantiles)
        shrink_rows.append({
            "cohort_id": c.cohort_id, "study": c.study, "cancer_type": c.cancer_type,
            "ici_target": c.ici_target, "platform": c.platform, "n": c.ni, "events": c.events,
            "unit_in": c.unit, "unit_factor": c.unit_factor, "se_source": c.se_source,
            "sign_flipped": c.sign_flipped, "access": c.access,
            "yi_observed": c.yi, "sei_observed": c.sei,
            "ratio_observed": _ratio(c.yi),
            "ci_low_observed": c.yi - zc * c.sei, "ci_high_observed": c.yi + zc * c.sei,
            "theta_post_mean": th["mean"], "theta_post_sd": th["sd"],
            "theta_ci_low": th[f"q{q_lo:g}"], "theta_ci_high": th[f"q{q_hi:g}"],
            "shrinkage_B": fit.shrinkage[i],
            "shrinkage_pct": 100 * fit.shrinkage[i],
            "weight_re_pct": 100 * float(w_re[i] / w_re.sum()),
            "z_observed": c.yi / c.sei,
            "std_residual_vs_pooled": (c.yi - mu_s["mean"]) / math.sqrt(c.sei ** 2 + tau_mean ** 2),
            "prob_theta_gt_0": fit.theta[i].prob_gt(0.0),
        })
    write_csv(os.path.join(outdir, "cohort_shrinkage.csv"), shrink_rows)

    write_csv(os.path.join(outdir, "tau_posterior.csv"),
              [{"tau": float(t), "density": float(d)} for t, d in zip(fit.tau, fit.tau_post)])
    mu_grid = np.linspace(mu_s["mean"] - 6 * mu_s["sd"], mu_s["mean"] + 6 * mu_s["sd"], 601)
    write_csv(os.path.join(outdir, "mu_posterior.csv"),
              [{"mu": float(x), "density": fit.mu.pdf(float(x)),
                "predictive_density": fit.pred.pdf(float(x))} for x in mu_grid])

    het_out = dict(het)
    het_out.update({"tau_post_mean": tau_mean, "tau_post_median": tau_med,
                    "tau_post_q97.5": tau_q975,
                    "prediction_interval_low": pred_s[f"q{q_lo:g}"],
                    "prediction_interval_high": pred_s[f"q{q_hi:g}"],
                    "prob_new_cohort_gt_0": fit.pred.prob_gt(0.0)})
    write_csv(os.path.join(outdir, "heterogeneity.csv"), [het_out])

    # --- Prior sensitivity -----------------------------------------------------------
    prior_rows = []
    for p in cfg.get("prior_sensitivity", []):
        if p.get("enabled") is False:
            prior_rows.append({"label": p.get("label", "?"), "status": "skipped (not configured)"})
            continue
        try:
            f2 = fit_bayes(y, V, p, cfg["tau_grid"])
        except Exception as exc:  # a badly specified sensitivity prior must not kill the run
            prior_rows.append({"label": p.get("label", "?"), "status": f"error: {exc}"})
            continue
        s2 = f2.mu.summary(quantiles)
        tw = f2.tau_weights()
        prior_rows.append({
            "label": p.get("label", "?"), "status": "ok",
            "mu_prior": f"N({p.get('mu_mean', 0)},{p.get('mu_sd', 1)}^2)",
            "tau_prior": f"{p.get('tau_family')}({p.get('tau_scale')})",
            "estimate": s2["mean"], "sd": s2["sd"],
            "ci_low": s2[f"q{q_lo:g}"], "ci_high": s2[f"q{q_hi:g}"],
                    "ratio": _ratio(s2["mean"]),
            "tau_post_mean": float(np.sum(f2.tau * tw)),
            "prob_gt_0": f2.mu.prob_gt(0.0),
            "log_marginal_likelihood": f2.log_marg,
            "edge_mass": f2.edge_mass,
        })
    write_csv(os.path.join(outdir, "prior_sensitivity.csv"), prior_rows)

    # --- Leave-one-out ---------------------------------------------------------------
    write_csv(os.path.join(outdir, "leave_one_out.csv"),
              leave_one_out(cohorts, V, cfg["prior"], cfg["tau_grid"], quantiles))

    # --- Small-study effects ---------------------------------------------------------
    pb_cfg = cfg["publication_bias"]
    bias_rows: List[Dict[str, Any]] = []
    if K >= 3:
        eg = egger_test(y, s, ni if np.all(np.isfinite(ni)) else None, pb_cfg["egger_predictor"])
        bias_rows.append({"test": f"egger[{eg['predictor']}]", "statistic": eg["t"],
                          "pvalue": eg["pvalue"], "estimate": eg["intercept"],
                          "note": "intercept != 0 indicates funnel asymmetry; "
                                  "underpowered when K < 10"})
        bg = begg_test(y, v)
        bias_rows.append({"test": "begg_rank", "statistic": bg["kendall_tau"],
                          "pvalue": bg["pvalue"], "estimate": float("nan"),
                          "note": "rank correlation; very low power at small K"})
        tf = trim_and_fill(y, v, pb_cfg.get("trimfill_estimator", "L0"))
        bias_rows.append({"test": f"trim_and_fill[{tf['side']}]", "statistic": float(tf["k0"]),
                          "pvalue": float("nan"), "estimate": tf["estimate"],
                          "ci_low": tf["ci_low"], "ci_high": tf["ci_high"],
                          "note": "k0 imputed cohorts; a what-if, not a corrected estimate"})
        pp = pet_peese(y, s)
        bias_rows.append({"test": "PET", "statistic": pp["pet_t"], "pvalue": pp["pet_pvalue"],
                          "estimate": pp["pet_estimate"], "note": "WLS intercept at SE=0"})
        bias_rows.append({"test": "PEESE", "statistic": float("nan"), "pvalue": pp["peese_pvalue"],
                          "estimate": pp["peese_estimate"],
                          "note": f"conditional rule selects: {pp['selected']}"})
        try:
            sm = selection_model_step(y, s, pb_cfg["selection_alpha"])
            bias_rows.append({"test": "selection_model_mle", "statistic": sm["omega"],
                              "pvalue": sm["lrt_pvalue"], "estimate": sm["estimate"],
                              "note": "step weight function; statistic = omega (relative "
                                      "probability that a non-significant cohort is observed)"})
        except Exception as exc:
            bias_rows.append({"test": "selection_model_mle", "statistic": float("nan"),
                              "pvalue": float("nan"), "estimate": float("nan"),
                              "note": f"did not converge: {exc}"})
        for om in pb_cfg.get("omega_grid", []):
            try:
                smf = selection_model_step(y, s, pb_cfg["selection_alpha"], omega_fixed=float(om))
                bias_rows.append({"test": f"selection_model_fixed_omega={om}",
                                  "statistic": float(om), "pvalue": float("nan"),
                                  "estimate": smf["estimate"],
                                  "note": "adjusted effect if non-significant cohorts were "
                                          f"{1/float(om):.3g}x less likely to be available"})
            except Exception:
                continue
    else:
        bias_rows.append({"test": "all", "statistic": float("nan"), "pvalue": float("nan"),
                          "estimate": float("nan"),
                          "note": "K < 3: small-study diagnostics not computed"})
    write_csv(os.path.join(outdir, "smallstudy_bias.csv"), bias_rows)

    # Gene-level null calibration (winner's curse on the gene axis)
    null_path = pb_cfg.get("null_gene_results") or ""
    if null_path and os.path.exists(null_path):
        null_rows = read_result_files([null_path])
        by_gene: Dict[str, List[Tuple[float, float]]] = {}
        for r in null_rows:
            yv, sv, _ = _derive_effect(r, {}, cfg.get("effect_family", "log_ratio"))
            g = _clean_str(r.get("gene"))
            if g and yv is not None and sv and sv > 0:
                by_gene.setdefault(g, []).append((yv, sv))
        null_stats = []
        for g, pairs in by_gene.items():
            if len(pairs) < 2:
                continue
            gy = np.array([a for a, _ in pairs]); gs = np.array([b for _, b in pairs])
            t2 = tau2_pm(gy, gs ** 2)
            r = random_effects(gy, gs ** 2, t2)
            null_stats.append({"gene": g, "k": len(pairs), "estimate": r["estimate"],
                               "se": r["se"], "z": r["estimate"] / r["se"]})
        if null_stats:
            obs_z = pooled[1]["estimate"] / pooled[1]["se"]
            zs = np.array([abs(r["z"]) for r in null_stats])
            emp_p = float((np.sum(zs >= abs(obs_z)) + 1) / (len(zs) + 1))
            write_csv(os.path.join(outdir, "gene_null_calibration.csv"),
                      null_stats + [{"gene": "__TARGET__", "k": K,
                                     "estimate": pooled[1]["estimate"], "se": pooled[1]["se"],
                                     "z": obs_z, "empirical_p_vs_null_panel": emp_p}])
            bias_rows.append({"test": "gene_panel_empirical_p", "statistic": obs_z,
                              "pvalue": emp_p, "estimate": pooled[1]["estimate"],
                              "note": f"target gene vs {len(zs)} negative-control genes "
                                      "pooled identically"})
            write_csv(os.path.join(outdir, "smallstudy_bias.csv"), bias_rows)

    # --- Missing cohorts -------------------------------------------------------------
    miss_cfg = cfg["missingness"]
    registry_path = miss_cfg.get("registry") or ""
    missing: List[MissingCohort] = []
    if registry_path and os.path.exists(registry_path):
        missing, reg_rows = read_registry(registry_path)
        write_csv(os.path.join(outdir, "registry_snapshot.csv"), reg_rows)
    if missing:
        med_se = float(np.median(s))
        se_missing, se_notes = [], []
        for mc in missing:
            se_i, how = impute_se(mc, cfg["primary_endpoint"], med_se, miss_cfg.get("se_rule", "auto"))
            se_missing.append(se_i)
            se_notes.append({"cohort_id": mc.cohort_id, "status": mc.status, "n": mc.ni,
                             "events": mc.events, "imputed_se": se_i, "se_rule": how,
                             "notes": mc.notes})
        write_csv(os.path.join(outdir, "missing_cohorts.csv"), se_notes)

        sweep = delta_sweep(y, V, se_missing, cfg["prior"], cfg["tau_grid"],
                            miss_cfg["delta_grid"], quantiles)
        # tipping point: smallest |delta| at which the credible interval covers 0
        tip = None
        for r in sorted(sweep, key=lambda r: abs(r["delta"])):
            if r["ci_low"] <= 0 <= r["ci_high"]:
                tip = r["delta"]
                break
        for r in sweep:
            r["is_tipping_point"] = (tip is not None and r["delta"] == tip)
        write_csv(os.path.join(outdir, "missingness_delta_sweep.csv"), sweep)

        lo_b = float(np.min(y) - 2 * tau_mean)
        hi_b = float(np.max(y) + 2 * tau_mean)
        bounds = manski_bounds(y, V, se_missing, cfg["prior"], cfg["tau_grid"], lo_b, hi_b, quantiles)
        fs = failsafe_k(y, V, cfg["prior"], cfg["tau_grid"],
                        miss_cfg.get("decision_threshold", 0.95),
                        max_k=int(miss_cfg.get("max_failsafe_k", 200)))
        summary_rows = [
            {"quantity": "n_cohorts_observed", "value": K, "note": ""},
            {"quantity": "n_cohorts_unobserved", "value": len(missing),
             "note": "; ".join(sorted({m.status for m in missing}))},
            {"quantity": "fraction_unobserved", "value": len(missing) / (K + len(missing)),
             "note": "cohort count, not patient count"},
            {"quantity": "patients_observed",
             "value": float(np.nansum(ni)) if np.any(np.isfinite(ni)) else float("nan"), "note": ""},
            {"quantity": "patients_unobserved",
             "value": float(sum(m.ni for m in missing if m.ni)) or float("nan"), "note": ""},
            {"quantity": "tipping_point_delta", "value": tip if tip is not None else float("nan"),
             "note": "smallest shift of unobserved cohorts (same scale as yi) at which the "
                     "credible interval covers 0; NA = not reached on the delta grid"},
            {"quantity": "tipping_point_ratio",
             "value": _ratio(tip) if tip is not None else float("nan"),
             "note": "exp(delta) when effect_family=log_ratio; NA for SMD"},
            {"quantity": "manski_envelope_low", "value": bounds["envelope_low"], "note": ""},
            {"quantity": "manski_envelope_high", "value": bounds["envelope_high"], "note": ""},
            {"quantity": "failsafe_k_null_cohorts", "value": fs["failsafe_k"],
             "note": f"cohorts with true effect 0 and SE={fs['se_null']:.3g} needed to push "
                     f"P(direction) below {miss_cfg.get('decision_threshold', 0.95)}"},
        ]
        write_csv(os.path.join(outdir, "missingness_summary.csv"), summary_rows)
    else:
        sweep = []
        write_csv(os.path.join(outdir, "missingness_summary.csv"),
                  [{"quantity": "registry", "value": "not provided",
                    "note": "set missingness.registry to a cohort registry CSV to run the "
                            "MNAR sensitivity analysis for EGA/dbGaP-restricted trials"}])

    # --- Meta-regression -------------------------------------------------------------
    mods = cfg["meta_regression"].get("moderators", [])
    if mods:
        design_cols = [np.ones(K)]
        names = ["intercept"]
        for m in mods:
            vals = [getattr(c, m, "") for c in cohorts]
            numeric = [_to_float(x) for x in vals]
            if all(x is not None for x in numeric):
                arr = np.array(numeric, float)
                design_cols.append((arr - arr.mean()) / (arr.std() or 1.0))
                names.append(f"{m}[z]")
            else:
                levels = sorted({str(x) for x in vals if str(x)})
                for lev in levels[1:]:
                    design_cols.append(np.array([1.0 if str(x) == lev else 0.0 for x in vals]))
                    names.append(f"{m}={lev}")
        X = np.column_stack(design_cols)
        mr = meta_regression(y, v, X, names)
        if "error" in mr:
            write_csv(os.path.join(outdir, "meta_regression.csv"),
                      [{"term": "ERROR", "estimate": float("nan"), "note": mr["error"]}])
        else:
            rows_mr = mr["coefficients"]
            for r in rows_mr:
                r["tau2_resid"] = mr["tau2_resid"]
                r["note"] = ("EXPLORATORY: cohort-level moderator; aggregation bias applies "
                             "and this is not a patient-level interaction")
            write_csv(os.path.join(outdir, "meta_regression.csv"), rows_mr)

    write_csv(os.path.join(outdir, "pooled_estimates.csv"), pooled)

    figs = maybe_plot(outdir, cohorts, fit, sweep)

    manifest = {
        "template": "methods/hier_meta/templates/hier_meta.py",
        "template_version": __version__,
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "python": sys.version.split()[0],
        "platform": _platform.platform(),
        "numpy": np.__version__,
        "scipy": __import__("scipy").__version__,
        "config": cfg,
        "inputs": [{"path": f, "sha256": sha256_file(f), "bytes": os.path.getsize(f)} for f in files],
        "n_input_rows": len(rows),
        "n_cohorts_pooled": K,
        "cohorts_pooled": [c.cohort_id for c in cohorts],
        "figures": figs,
        "warnings": ([f"tau grid edge mass {fit.edge_mass:.3f}"] if fit.edge_mass > 0.01 else []),
    }
    with open(os.path.join(outdir, "run_manifest.json"), "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2, default=str)

    scale_note = "SMD (Hedges' g, harm orientation)" if is_smd else (
        f"log-scale; ratio {_ratio(mu_s['mean']):.3f}")
    print(f"[hier_meta] {cfg['gene']} {cfg['primary_endpoint']} "
          f"pooled ({cfg['prior']['tau_family']} prior): "
          f"{mu_s['mean']:.3f} [{mu_s[f'q{q_lo:g}']:.3f}, {mu_s[f'q{q_hi:g}']:.3f}] "
          f"{scale_note}; tau {tau_mean:.3f}; P(>0)={fit.mu.prob_gt(0.0):.3f}")
    print(f"[hier_meta] wrote {len(os.listdir(outdir))} file(s) to {outdir}")
    print("[hier_meta] REMINDER: this is a methods run. Interpretation belongs in the "
          "prespecified report, and the missing-cohort sensitivity bounds the claim.")
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="Bayesian hierarchical pooling of a single-gene ICI effect across small cohorts.")
    ap.add_argument("--results-dir", default="results",
                    help="directory searched for *.csv per-cohort estimates (default: results)")
    ap.add_argument("--recursive", action="store_true", help="search results-dir recursively")
    ap.add_argument("--config", default=None, help="JSON config; see --write-config")
    ap.add_argument("--outdir", default=os.path.join("methods", "hier_meta", "out"),
                    help="output directory")
    ap.add_argument("--gene", default=None, help="override config gene (single)")
    ap.add_argument("--genes", default=None,
                    help="comma-separated genes (default: TACSTD2,CLDN4); writes one subdir each")
    ap.add_argument("--endpoint", default=None, help="override config primary_endpoint")
    ap.add_argument("--effect-family", default=None, choices=("smd", "log_ratio"),
                    help="override config effect_family")
    ap.add_argument("--access-policy", default=None, choices=("open_only", "catalog"),
                    help="open_only drops EGA/dbGaP rows (default)")
    ap.add_argument("--registry", default=None, help="override missingness.registry path")
    ap.add_argument("--write-config", default=None, metavar="PATH",
                    help="write the default config to PATH and exit")
    ap.add_argument("--demo", action="store_true",
                    help="run on the bundled synthetic example instead of results/")
    ap.add_argument("--strict", action="store_true",
                    help="exit non-zero when no usable inputs are found")
    args = ap.parse_args(argv)

    if args.write_config:
        with open(args.write_config, "w", encoding="utf-8") as fh:
            json.dump(DEFAULT_CONFIG, fh, indent=2)
        print(f"wrote default config to {args.write_config}")
        return 0

    cfg = copy.deepcopy(DEFAULT_CONFIG)
    if args.config:
        with open(args.config, encoding="utf-8") as fh:
            cfg = deep_update(cfg, json.load(fh))
    if args.gene:
        cfg["gene"] = args.gene
        cfg["genes"] = [args.gene]
    if args.genes:
        cfg["genes"] = [g.strip() for g in args.genes.split(",") if g.strip()]
    if args.endpoint:
        cfg["primary_endpoint"] = args.endpoint
    if args.effect_family:
        cfg["effect_family"] = args.effect_family
    if args.access_policy:
        cfg["access_policy"] = args.access_policy
    if args.registry:
        cfg["missingness"]["registry"] = args.registry

    results_dir = args.results_dir
    if args.demo:
        here = os.path.dirname(os.path.abspath(__file__))
        results_dir = os.path.join(os.path.dirname(here), "example", "results")
        if not cfg["missingness"].get("registry"):
            reg = os.path.join(os.path.dirname(here), "example", "cohort_registry.csv")
            if os.path.exists(reg):
                cfg["missingness"]["registry"] = reg
        print(f"[hier_meta] DEMO MODE: synthetic inputs from {results_dir}. "
              "Numbers produced here are simulated and mean nothing.")

    if not os.path.isdir(results_dir):
        print(f"[hier_meta] results directory '{results_dir}' does not exist.")
        print(REQUIRED_MESSAGE)
        return 1 if args.strict else 0
    files = find_result_files(results_dir, args.recursive)
    if not files:
        print(f"[hier_meta] no CSV files under '{results_dir}'.")
        print(REQUIRED_MESSAGE)
        return 1 if args.strict else 0

    genes = cfg.get("genes") or [cfg["gene"]]
    rc = 0
    for gene in genes:
        one = copy.deepcopy(cfg)
        one["gene"] = gene
        out = args.outdir if len(genes) == 1 else os.path.join(args.outdir, gene)
        print(f"[hier_meta] === {gene} / {one['primary_endpoint']} / "
              f"{one.get('effect_family')} / access={one.get('access_policy')} ===")
        rc = max(rc, run(one, files, out))
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
