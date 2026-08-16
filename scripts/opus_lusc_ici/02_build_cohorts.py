#!/usr/bin/env python3
"""Step 2 - harmonise every cohort into a common expression + clinical schema.

Each cohort is reduced to:
  * ``<KEY>_expr.tsv.gz``  genes (HGNC symbols) x samples, log2 scale
  * one row per sample in ``cohort_clinical.csv`` with a shared column schema

Histology is taken from the GEO pathology annotation where it exists; cohorts
without it are left as ``NA`` here and resolved in step 3.
"""

from __future__ import annotations

import gzip
import sys

import numpy as np
import pandas as pd

from common import (collapse_to_symbols, counts_to_logcpm, load_probemap,
                    read_series_matrix, to_log2p1, write_table)
from config import COHORTS_BY_KEY, DATA_DIR, LOGS_DIR, RAW_DIR, TABLES_DIR

# Shared clinical schema. ``benefit`` is the harmonised binary endpoint used by
# the response meta-analysis; ``benefit_definition`` records how it was derived.
CLINICAL_COLUMNS = [
    "cohort", "sample_id", "gsm", "histology", "histology_source", "setting",
    "regimen", "timepoint", "benefit", "benefit_definition", "response_raw",
    "pfs_time_months", "pfs_event", "sex", "age", "stage",
]


def empty_clinical(index) -> pd.DataFrame:
    return pd.DataFrame({c: pd.Series([np.nan] * len(index), index=index) for c in CLINICAL_COLUMNS})


def norm_histology(value: str) -> str:
    if not isinstance(value, str):
        return "NA"
    v = value.strip().lower()
    if v in {"squamous", "sqcc", "squamous cell carcinoma", "lung squamous cell carcinoma"}:
        return "LUSC"
    if v in {"adeno", "adenocarcinoma", "ac", "non-squamous"}:
        return "non-LUSC"
    if v in {"", "na", "nan", "nsclc", "other", "nos", "sarcomatoid",
             "carcinosarcomatoid", "large cell neuroendocarine carcinoma"}:
        return "NA"
    return "NA"


# ---------------------------------------------------------------------------
# Cohort-specific loaders. Each returns (expr [genes x samples, log2], clinical)
# ---------------------------------------------------------------------------


def build_france(key: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """GSE190265 / GSE190266 - Dijon anti-PD-1 monotherapy, TPM, samples in rows."""
    cohort = COHORTS_BY_KEY[key]
    filename = cohort.supplementary[0]
    with gzip.open(RAW_DIR / filename, "rt", encoding="utf8", errors="replace") as fh:
        text = fh.read()
    # France-4 uses a comma decimal separator; both use ';' as field separator.
    decimal = "," if key == "FRANCE4" else "."
    frame = pd.read_csv(pd.io.common.StringIO(text), sep=";", decimal=decimal,
                        index_col=0, low_memory=False)
    frame.index = frame.index.astype(str)
    if key == "FRANCE3":
        # France-3 has no index header: the first column holds the sample id.
        frame.index.name = "sample"
    expr = to_log2p1(frame.T.apply(pd.to_numeric, errors="coerce"))
    expr = collapse_to_symbols(expr)

    meta = read_series_matrix(cohort.gse)
    meta["sample_id"] = meta["title"].astype(str)
    meta = meta.set_index("sample_id")

    clinical = empty_clinical(expr.columns)
    clinical["cohort"] = key
    clinical["sample_id"] = expr.columns
    clinical["gsm"] = [meta["gsm"].get(s, np.nan) for s in expr.columns]
    clinical["histology"] = [norm_histology(meta["disease state"].get(s)) for s in expr.columns]
    clinical["histology_source"] = "pathology (GEO)"
    clinical["setting"] = cohort.setting
    clinical["regimen"] = cohort.regimen
    clinical["timepoint"] = "pre-treatment"

    if key == "FRANCE3":
        info = pd.read_csv(RAW_DIR / "GSE190265_samples_info_France3.csv.gz", sep=";")
        info["sample"] = info["sample"].astype(str)
        info = info.set_index("sample")
        clinical["pfs_time_months"] = [pd.to_numeric(info["time_PFS"].get(s), errors="coerce")
                                       for s in expr.columns]
        clinical["pfs_event"] = [pd.to_numeric(info["evtPFS"].get(s), errors="coerce")
                                 for s in expr.columns]
    else:
        clinical["pfs_time_months"] = [pd.to_numeric(meta["pfs_time (6 months)"].get(s), errors="coerce")
                                       for s in expr.columns]
        clinical["pfs_event"] = [pd.to_numeric(meta["pfs_evt (6 months)"].get(s), errors="coerce")
                                 for s in expr.columns]

    # Durable clinical benefit = progression-free at the 6-month landmark.
    # Patients censored before 6 months cannot be classified and stay NA.
    benefit = []
    for t, e in zip(clinical["pfs_time_months"], clinical["pfs_event"]):
        if pd.isna(t) or pd.isna(e):
            benefit.append(np.nan)
        elif t >= 6:
            benefit.append(1)
        elif e == 1:
            benefit.append(0)
        else:
            benefit.append(np.nan)
    clinical["benefit"] = benefit
    clinical["benefit_definition"] = "progression-free at 6 months (landmark)"
    clinical["response_raw"] = [f"PFS={t}m event={e}" for t, e in
                                zip(clinical["pfs_time_months"], clinical["pfs_event"])]
    return expr, clinical


def build_neochemo() -> tuple[pd.DataFrame, pd.DataFrame]:
    """GSE207422 - neoadjuvant anti-PD-1 + chemotherapy, bulk log2 TPM."""
    cohort = COHORTS_BY_KEY["NEOCHEMO"]
    expr = pd.read_csv(RAW_DIR / cohort.supplementary[0], sep="\t", index_col=0)
    expr = collapse_to_symbols(expr.apply(pd.to_numeric, errors="coerce"))

    meta = read_series_matrix(cohort.gse).set_index("title")
    extra = pd.read_excel(RAW_DIR / cohort.supplementary[1]).set_index("Sample")

    clinical = empty_clinical(expr.columns)
    clinical["cohort"] = "NEOCHEMO"
    clinical["sample_id"] = expr.columns
    clinical["gsm"] = [meta["gsm"].get(s, np.nan) for s in expr.columns]
    clinical["histology"] = [norm_histology(meta["pathologic_sub"].get(s)) for s in expr.columns]
    clinical["histology_source"] = "pathology (GEO)"
    clinical["setting"] = cohort.setting
    clinical["regimen"] = cohort.regimen
    clinical["timepoint"] = ["pre-treatment" if str(meta["sampling_time"].get(s, "")).startswith("Pre")
                             else "post-treatment" for s in expr.columns]
    clinical["sex"] = [meta["Sex"].get(s, np.nan) for s in expr.columns]
    clinical["age"] = [pd.to_numeric(meta["age"].get(s), errors="coerce") for s in expr.columns]
    clinical["stage"] = [meta["clinical_stage"].get(s, np.nan) for s in expr.columns]

    raw = [str(meta["pathologic_response"].get(s, "")) for s in expr.columns]
    clinical["response_raw"] = [f"{r} / RECIST {extra['RECIST'].get(s, 'NA')}"
                                for r, s in zip(raw, expr.columns)]
    clinical["benefit"] = [1 if r.startswith("MPR") else (0 if r == "NMPR" else np.nan) for r in raw]
    clinical["benefit_definition"] = "major pathological response (MPR, including pCR)"
    return expr, clinical


def build_pla() -> tuple[pd.DataFrame, pd.DataFrame]:
    """GSE283829 - advanced NSCLC on ICI, raw counts with Ensembl ids."""
    cohort = COHORTS_BY_KEY["PLA"]
    counts = pd.read_csv(RAW_DIR / cohort.supplementary[0], sep="\t", index_col=0)
    counts.columns = [c.lstrip("S") for c in counts.columns]
    expr = collapse_to_symbols(counts_to_logcpm(counts.apply(pd.to_numeric, errors="coerce")),
                               load_probemap())

    meta = read_series_matrix(cohort.gse).set_index("title")
    clinical = empty_clinical(expr.columns)
    clinical["cohort"] = "PLA"
    clinical["sample_id"] = expr.columns
    clinical["gsm"] = [meta["gsm"].get(s, np.nan) for s in expr.columns]
    clinical["histology"] = [norm_histology(meta["tumor type"].get(s)) for s in expr.columns]
    clinical["histology_source"] = "pathology (GEO)"
    clinical["setting"] = cohort.setting
    clinical["regimen"] = cohort.regimen
    clinical["timepoint"] = "pre-treatment"
    clinical["sex"] = [meta["Sex"].get(s, np.nan) for s in expr.columns]

    raw = [str(meta["disease stage"].get(s, "")) for s in expr.columns]
    clinical["response_raw"] = raw
    # 'disease stage' in this series records RECIST best response (CR/SD/PD).
    clinical["benefit"] = [1 if r == "CR" else (0 if r in {"SD", "PD"} else np.nan) for r in raw]
    clinical["benefit_definition"] = "complete response (CR) vs SD/PD"
    return expr, clinical


def build_durvart() -> tuple[pd.DataFrame, pd.DataFrame]:
    """GSE253564 - neoadjuvant durvalumab +/- SBRT, pre-treatment FPKM."""
    cohort = COHORTS_BY_KEY["DURVART"]
    frame = pd.read_csv(RAW_DIR / cohort.supplementary[0], sep="\t")
    frame = frame.drop(columns=[c for c in ("Entrez.ID",) if c in frame.columns])
    frame = frame.set_index("gene")
    expr = collapse_to_symbols(to_log2p1(frame.apply(pd.to_numeric, errors="coerce")))

    meta = read_series_matrix(cohort.gse).set_index("title")
    clinical = empty_clinical(expr.columns)
    clinical["cohort"] = "DURVART"
    clinical["sample_id"] = expr.columns
    clinical["gsm"] = [meta["gsm"].get(s, np.nan) for s in expr.columns]
    clinical["histology"] = [norm_histology(meta["cell type"].get(s)) for s in expr.columns]
    clinical["histology_source"] = "pathology (GEO)"
    clinical["setting"] = cohort.setting
    clinical["regimen"] = cohort.regimen
    clinical["timepoint"] = "pre-treatment"
    clinical["response_raw"] = [meta["treatment"].get(s, np.nan) for s in expr.columns]
    clinical["benefit_definition"] = "not available in the public record"
    return expr, clinical


def build_smc() -> tuple[pd.DataFrame, pd.DataFrame]:
    """GSE135222 - SMC Korea anti-PD-1/PD-L1, TPM with Ensembl ids, no histology in GEO."""
    cohort = COHORTS_BY_KEY["SMC"]
    frame = pd.read_csv(RAW_DIR / cohort.supplementary[0], sep="\t", index_col=0)
    expr = collapse_to_symbols(to_log2p1(frame.apply(pd.to_numeric, errors="coerce")),
                               load_probemap())

    meta = read_series_matrix(cohort.gse)
    meta["key"] = meta["title"].str.replace(" ", "", regex=False)
    meta = meta.set_index("key")

    clinical = empty_clinical(expr.columns)
    clinical["cohort"] = "SMC"
    clinical["sample_id"] = expr.columns
    clinical["gsm"] = [meta["gsm"].get(s, np.nan) for s in expr.columns]
    clinical["histology"] = "NA"          # resolved by the step-3 classifier
    clinical["histology_source"] = "inferred (transcriptomic)"
    clinical["setting"] = cohort.setting
    clinical["regimen"] = cohort.regimen
    clinical["timepoint"] = "pre-treatment"
    clinical["sex"] = [meta["gender"].get(s, np.nan) for s in expr.columns]
    clinical["age"] = [pd.to_numeric(meta["age"].get(s), errors="coerce") for s in expr.columns]

    days = pd.Series([pd.to_numeric(meta["pfs.time"].get(s), errors="coerce") for s in expr.columns],
                     index=expr.columns)
    event = pd.Series([pd.to_numeric(meta["progression-free survival (pfs)"].get(s), errors="coerce")
                       for s in expr.columns], index=expr.columns)
    clinical["pfs_time_months"] = days / 30.44
    clinical["pfs_event"] = event.values
    benefit = []
    for t, e in zip(clinical["pfs_time_months"], clinical["pfs_event"]):
        if pd.isna(t) or pd.isna(e):
            benefit.append(np.nan)
        elif t >= 6:
            benefit.append(1)
        elif e == 1:
            benefit.append(0)
        else:
            benefit.append(np.nan)
    clinical["benefit"] = benefit
    clinical["benefit_definition"] = "progression-free at 6 months (landmark)"
    clinical["response_raw"] = [f"PFS={d}d event={e}" for d, e in zip(days, event)]
    return expr, clinical


def build_yuhs() -> tuple[pd.DataFrame, pd.DataFrame]:
    """GSE126044 - Yonsei anti-PD-1, raw counts by symbol, no histology in GEO."""
    cohort = COHORTS_BY_KEY["YUHS"]
    counts = pd.read_csv(RAW_DIR / cohort.supplementary[0], sep="\t", index_col=0)
    expr = collapse_to_symbols(counts_to_logcpm(counts.apply(pd.to_numeric, errors="coerce")))

    meta = read_series_matrix(cohort.gse)
    meta["key"] = meta["title"].str.replace("RNA-seq_", "", regex=False)
    meta = meta.set_index("key")

    clinical = empty_clinical(expr.columns)
    clinical["cohort"] = "YUHS"
    clinical["sample_id"] = expr.columns
    clinical["gsm"] = [meta["gsm"].get(s, np.nan) for s in expr.columns]
    clinical["histology"] = "NA"
    clinical["histology_source"] = "inferred (transcriptomic)"
    clinical["setting"] = cohort.setting
    clinical["regimen"] = cohort.regimen
    clinical["timepoint"] = "pre-treatment"
    raw = [str(meta["patient response"].get(s, "")) for s in expr.columns]
    clinical["response_raw"] = raw
    clinical["benefit"] = [1 if r == "responder" else (0 if r == "non-responder" else np.nan)
                           for r in raw]
    clinical["benefit_definition"] = "RECIST responder as reported by the submitters"
    return expr, clinical


BUILDERS = {
    "FRANCE3": lambda: build_france("FRANCE3"),
    "FRANCE4": lambda: build_france("FRANCE4"),
    "NEOCHEMO": build_neochemo,
    "PLA": build_pla,
    "DURVART": build_durvart,
    "SMC": build_smc,
    "YUHS": build_yuhs,
}

# Deposited matrices differ enormously in how many genes they carry (16k for the
# reduced France-4 panel, 59k for the Ensembl-complete ones), so library
# complexity is judged relative to each cohort's own median rather than against
# one absolute cut-off.
QC_ABSOLUTE_FLOOR = 2000
QC_RELATIVE_FRACTION = 0.4


def main() -> int:
    all_clinical, qc_rows = [], []

    for key, builder in BUILDERS.items():
        print(f"[{key}]")
        expr, clinical = builder()

        detected = (expr > 0).sum(axis=0)
        threshold = max(QC_ABSOLUTE_FLOOR, QC_RELATIVE_FRACTION * float(np.median(detected)))
        keep = detected >= threshold
        dropped = list(expr.columns[~keep])
        if dropped:
            print(f"  QC: dropping {len(dropped)} sample(s) below {threshold:.0f} "
                  f"detected genes: {dropped}")
        expr = expr.loc[:, keep]
        clinical = clinical.loc[keep]

        expr = expr.round(4)
        out = DATA_DIR / f"{key}_expr.tsv.gz"
        expr.to_csv(out, sep="\t", compression={"method": "gzip", "mtime": 0})
        print(f"  {expr.shape[0]} genes x {expr.shape[1]} samples -> {out.name} "
              f"({out.stat().st_size/1e6:.1f} MB)")

        qc_rows.append(dict(
            cohort=key, genes=expr.shape[0], samples=expr.shape[1],
            dropped_samples=len(dropped),
            qc_threshold=int(threshold),
            median_detected_genes=int(np.median(detected[keep])) if keep.any() else 0,
            lusc_pathology=int((clinical["histology"] == "LUSC").sum()),
            has_TACSTD2=int("TACSTD2" in expr.index),
            has_CLDN4=int("CLDN4" in expr.index),
        ))
        all_clinical.append(clinical)

    clinical = pd.concat(all_clinical, ignore_index=True)
    write_table(clinical, TABLES_DIR / "cohort_clinical.csv")
    write_table(pd.DataFrame(qc_rows), TABLES_DIR / "cohort_qc.csv")

    print("\nPer-cohort summary")
    print(pd.DataFrame(qc_rows).to_string(index=False))
    print("\nHistology (pathology-annotated cohorts)")
    print(clinical.groupby(["cohort", "histology"]).size().to_string())

    (LOGS_DIR / "02_build_cohorts.done").write_text("ok\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
