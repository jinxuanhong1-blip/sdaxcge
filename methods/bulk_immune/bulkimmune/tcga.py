"""TCGA LUAD / LUSC loaders (UCSC Xena / GDC-hub public mirrors).

TCGA is not an ICI cohort. It is used here as a large, treatment-naive NSCLC
expression atlas to estimate the *direction and magnitude* of the
TACSTD2/CLDN4 ~ immune-score correlations with enough n to survive
multiple-testing, and to stratify LUAD vs LUSC. ICI response is then tested
separately in GSE126044 / GSE135222. Do not quote a TCGA correlation as
evidence that TACSTD2 predicts ICI outcome.

Files (public, no login):

* Expression (log2 RSEM, HiSeqV2, gene symbols):
  ``https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.{COHORT}.sampleMap/HiSeqV2.gz``
* Clinical matrix:
  ``https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.{COHORT}.sampleMap/{COHORT}_clinicalMatrix.gz``

Primary-tumour samples are kept (barcode sample type ``01``).
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .geo import download

__all__ = ["XENA_URLS", "fetch_tcga_cohort", "load_tcga_nsclc", "barcode_sample_type"]

XENA_URLS = {
    "LUAD": {
        "expr": "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LUAD.sampleMap%2FHiSeqV2.gz",
        # clinical matrix is stored uncompressed on this mirror (.gz returns 403)
        "clin": "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LUAD.sampleMap%2FLUAD_clinicalMatrix",
    },
    "LUSC": {
        "expr": "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LUSC.sampleMap%2FHiSeqV2.gz",
        "clin": "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LUSC.sampleMap%2FLUSC_clinicalMatrix",
    },
}


def barcode_sample_type(barcode: str) -> str:
    """TCGA sample-type code (01 = primary solid tumour, 11 = solid tissue normal)."""
    parts = str(barcode).split("-")
    if len(parts) < 4:
        return ""
    return parts[3][:2]


def fetch_tcga_cohort(cohort: str, dest_dir: str | Path) -> tuple[Path, Path]:
    if cohort not in XENA_URLS:
        raise ValueError(f"unknown cohort {cohort}")
    dest_dir = Path(dest_dir)
    expr = download(XENA_URLS[cohort]["expr"], dest_dir / f"TCGA_{cohort}_HiSeqV2.gz")
    clin = download(XENA_URLS[cohort]["clin"], dest_dir / f"TCGA_{cohort}_clinicalMatrix")
    return expr, clin


def load_tcga_nsclc(
    dest_dir: str | Path,
    cohorts: tuple[str, ...] = ("LUAD", "LUSC"),
    sample_types: tuple[str, ...] = ("01",),
    fetch: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (log2-RSEM genes x samples, phenotype) for the requested cohorts.

    Phenotype columns: cohort, sample_type, sex, age, pathologic_stage,
    os_time, os_event (when present in the Xena clinical matrix).
    """
    dest_dir = Path(dest_dir)
    expr_frames = []
    pheno_frames = []
    for cohort in cohorts:
        if fetch:
            epath, cpath = fetch_tcga_cohort(cohort, dest_dir)
        else:
            epath = dest_dir / f"TCGA_{cohort}_HiSeqV2.gz"
            cpath = dest_dir / f"TCGA_{cohort}_clinicalMatrix"
        expr = pd.read_csv(epath, sep="\t", index_col=0, compression="infer")
        expr.index = expr.index.astype(str)
        clin = pd.read_csv(cpath, sep="\t", index_col=0, compression="infer", dtype=str)
        keep = [c for c in expr.columns if barcode_sample_type(c) in sample_types]
        expr = expr[keep]
        clin = clin.reindex(keep)
        clin = clin.copy()
        clin["cohort"] = cohort
        clin["sample_type"] = [barcode_sample_type(i) for i in clin.index]
        # Xena column names vary slightly across cohorts
        def _first(*cands):
            lower = {c.lower(): c for c in clin.columns}
            for cand in cands:
                if cand.lower() in lower:
                    return lower[cand.lower()]
            for c in clin.columns:
                if any(cand.lower() in c.lower() for cand in cands):
                    return c
            return None

        sex_c = _first("gender", "sex")
        age_c = _first("age_at_initial_pathologic_diagnosis", "age")
        stage_c = _first("pathologic_stage", "ajcc_pathologic_tumor_stage")
        os_t = _first("OS.time", "os_time", "vital_status")
        # Xena often uses _OS_IND / _OS
        os_time_c = _first("OS.time", "_OS", "os_time")
        os_event_c = _first("OS", "_OS_IND", "vital_status")
        pheno = pd.DataFrame(
            {
                "cohort": cohort,
                "sample_type": clin["sample_type"],
                "sex": clin[sex_c].str.lower() if sex_c else pd.NA,
                "age": pd.to_numeric(clin[age_c], errors="coerce") if age_c else pd.NA,
                "pathologic_stage": clin[stage_c] if stage_c else pd.NA,
                "os_time": pd.to_numeric(clin[os_time_c], errors="coerce") if os_time_c else pd.NA,
                "os_event": pd.to_numeric(clin[os_event_c], errors="coerce") if os_event_c else pd.NA,
            },
            index=expr.columns,
        )
        expr_frames.append(expr)
        pheno_frames.append(pheno)

    # Inner-join genes so LUAD and LUSC share a symbol universe
    genes = set(expr_frames[0].index)
    for fr in expr_frames[1:]:
        genes &= set(fr.index)
    genes = [g for g in expr_frames[0].index if g in genes]
    expr_all = pd.concat([fr.loc[genes] for fr in expr_frames], axis=1)
    pheno_all = pd.concat(pheno_frames, axis=0)
    pheno_all = pheno_all.loc[expr_all.columns]
    return expr_all, pheno_all
