#!/usr/bin/env python3
"""Step 3 - build the TCGA reference used for calibration and replication.

TCGA-LUSC provides a large (n~500), pathology-confirmed, treatment-naive lung
squamous reference in which the TACSTD2/CLDN4 immune associations can be
estimated with real precision; TCGA-LUAD is used only to calibrate the
histology classifier of step 4, so only its summary scores are kept.
"""

from __future__ import annotations

import sys

import numpy as np
import pandas as pd

from common import collapse_to_symbols, load_probemap, sample_percentile_ranks, write_table
from config import DATA_DIR, LOGS_DIR, RAW_DIR, TABLES_DIR

# Xena GDC "star_tpm" matrices are already log2(TPM + 1).
MIN_DETECTION_FRACTION = 0.10
MIN_DETECTION_LEVEL = 0.10

UNIVERSE_COHORTS = ["FRANCE3", "NEOCHEMO", "PLA", "DURVART", "SMC", "YUHS"]


def load_tcga(filename: str) -> pd.DataFrame:
    frame = pd.read_csv(RAW_DIR / filename, sep="\t", index_col=0)
    primary = [c for c in frame.columns if c.split("-")[3].startswith("01")]
    frame = frame[primary]
    # One RNA-seq aliquot per patient.
    patients = pd.Series([c[:12] for c in frame.columns], index=frame.columns)
    frame = frame.loc[:, ~patients.duplicated()]
    frame.columns = [c[:12] for c in frame.columns]
    return collapse_to_symbols(frame, load_probemap())


def main() -> int:
    print("[TCGA-LUSC] loading")
    lusc = load_tcga("TCGA-LUSC.star_tpm.tsv.gz")
    print(f"  {lusc.shape[0]} genes x {lusc.shape[1]} primary tumours")

    print("[TCGA-LUAD] loading")
    luad = load_tcga("TCGA-LUAD.star_tpm.tsv.gz")
    print(f"  {luad.shape[0]} genes x {luad.shape[1]} primary tumours")

    shared = lusc.index.intersection(luad.index)
    print(f"  {len(shared)} shared symbols")

    # A within-sample percentile rank only means the same thing across datasets
    # if it is taken over the same set of genes, so the classifier of step 4
    # works on one explicit universe shared by TCGA and the ICI cohorts.
    # France-4's deposited matrix carries a reduced 16k-gene panel that omits
    # several histology markers; it is pathology-annotated and therefore does
    # not need the classifier, so it is left out of the universe.
    universe = set(shared)
    for key in UNIVERSE_COHORTS:
        genes = pd.read_csv(DATA_DIR / f"{key}_expr.tsv.gz", sep="\t", usecols=[0]).iloc[:, 0]
        universe &= set(genes)
        print(f"  universe after {key}: {len(universe)}")
    universe = sorted(universe)
    (DATA_DIR / "analysis_gene_universe.txt").write_text("\n".join(universe) + "\n")
    print(f"  wrote analysis_gene_universe.txt ({len(universe)} genes)")

    from common import load_signatures

    sigs = load_signatures()
    markers = [g for g in sigs["SQUAMOUS_MARKERS"]["genes"] + sigs["ADENO_MARKERS"]["genes"]
               if g in universe]
    missing = [g for g in sigs["SQUAMOUS_MARKERS"]["genes"] + sigs["ADENO_MARKERS"]["genes"]
               if g not in universe]
    if missing:
        print(f"  histology markers absent from the universe: {missing}")

    for name, frame in (("LUSC", lusc), ("LUAD", luad)):
        ranks = sample_percentile_ranks(frame.loc[universe])
        out = DATA_DIR / f"TCGA_{name}_percentile_ranks_markers.tsv.gz"
        ranks.loc[markers].round(5).to_csv(out, sep="\t", compression={"method": "gzip", "mtime": 0})
        print(f"  wrote {out.name} ({len(markers)} marker genes x {frame.shape[1]} samples)")

    # Keep the LUSC expression matrix itself: it powers the replication of the
    # immune-context analysis. Drop genes that are essentially never detected.
    detected = (lusc > MIN_DETECTION_LEVEL).mean(axis=1)
    keep = detected >= MIN_DETECTION_FRACTION
    lusc_out = lusc.loc[keep].round(3)
    print(f"  keeping {keep.sum()} / {len(keep)} genes detected in >= "
          f"{MIN_DETECTION_FRACTION:.0%} of tumours")

    path = DATA_DIR / "TCGA_LUSC_expr.tsv.gz"
    lusc_out.to_csv(path, sep="\t", compression={"method": "gzip", "mtime": 0})
    print(f"  wrote {path.name} ({path.stat().st_size/1e6:.1f} MB)")

    clinical = pd.read_csv(RAW_DIR / "TCGA-LUSC.clinical.tsv.gz", sep="\t", low_memory=False)
    clinical["patient"] = clinical["sample"].str[:12]
    columns = {
        "patient": "patient",
        "gender.demographic": "sex",
        "age_at_index.demographic": "age",
        "ajcc_pathologic_stage.diagnoses": "stage",
        "primary_diagnosis.diagnoses": "primary_diagnosis",
        "vital_status.demographic": "vital_status",
        "days_to_death.demographic": "days_to_death",
        "days_to_last_follow_up.diagnoses": "days_to_last_follow_up",
    }
    available = {k: v for k, v in columns.items() if k in clinical.columns}
    slim = clinical[list(available)].rename(columns=available).drop_duplicates("patient")
    slim = slim[slim["patient"].isin(lusc_out.columns)]
    slim["os_event"] = (slim["vital_status"].astype(str).str.lower() == "dead").astype(int)
    slim["os_days"] = pd.to_numeric(slim.get("days_to_death"), errors="coerce").fillna(
        pd.to_numeric(slim.get("days_to_last_follow_up"), errors="coerce"))
    write_table(slim, TABLES_DIR / "tcga_lusc_clinical.csv")

    summary = pd.DataFrame([
        dict(dataset="TCGA-LUSC", genes=lusc_out.shape[0], samples=lusc_out.shape[1],
             file_mb=round(path.stat().st_size / 1e6, 1)),
        dict(dataset="TCGA-LUAD (markers only)", genes=len(shared), samples=luad.shape[1],
             file_mb=np.nan),
    ])
    write_table(summary, TABLES_DIR / "tcga_reference_summary.csv")

    (LOGS_DIR / "03_tcga_reference.done").write_text("ok\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
