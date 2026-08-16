#!/usr/bin/env python3
"""Replace the broken ORCESTRA IMvigor210 CLDN4 vector with official counts.

ORCESTRA ICB_Mariathasan_expr.tsv stores CLDN4 (ENSG00000189143) as a constant
floor (-9.96578) in all 348 samples. Housekeeping and other epithelial genes in
the same object vary, so this is a gene-specific quantification failure, not a
row-alignment error.

Official counts come from IMvigor210CoreBiologies (Mariathasan et al. Nature
2018), extracted copy:
  https://github.com/mimifp/tfm_mUC/tree/main/2_data/1_IMvigor210
which is a CSV dump of counts(cds) / fData / pData from
  http://research-pub.gene.com/IMvigor210CoreBiologies/

Patient IDs match PredictIO as P{ANONPT_ID}.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "imvigor210"
DER = ROOT / "data" / "derived"
OUT = ROOT / "results" / "claim_B5"


def main() -> None:
    f = pd.read_csv(RAW / "fData_IMvigor210.csv")
    iloc = int(np.flatnonzero(f["symbol"].eq("CLDN4"))[0])
    row = pd.read_csv(
        RAW / "exmat_censored_IMvigor210.csv",
        skiprows=lambda i: i != 0 and i != (iloc + 1),
    )
    assert str(row.iloc[0, 0]) == f"gene_{iloc + 1}", (row.iloc[0, 0], iloc)
    vals = row.iloc[0, 1:].astype(float)
    p = pd.read_csv(RAW / "pData_IMvigor210.csv")
    p["CLDN4_counts"] = p["X"].map(vals)
    p["CLDN4_log2"] = np.log2(p["CLDN4_counts"] + 1.0)
    # One ANONPT_ID (10285) has two NE samples; PredictIO labels them P10285 / P10285_1.
    p = p.sort_values(["ANONPT_ID", "X"]).copy()
    p["dup"] = p.groupby("ANONPT_ID").cumcount()
    p["patientid"] = "P" + p["ANONPT_ID"].astype(str)
    p.loc[p["dup"] > 0, "patientid"] = p.loc[p["dup"] > 0, "patientid"] + "_" + p.loc[p["dup"] > 0, "dup"].astype(str)

    clin = pd.read_csv(DER / "cldn4_clinical.tsv", sep="\t")
    mask = clin["study"].eq("Mariathasan")
    before = clin.loc[mask, "CLDN4"].copy()
    mapped = clin.loc[mask, "patientid"].map(p.set_index("patientid")["CLDN4_log2"])
    n_hit = int(mapped.notna().sum())
    if n_hit != int(mask.sum()):
        raise SystemExit(f"ID match failed: {n_hit}/{int(mask.sum())}")
    clin.loc[mask, "CLDN4"] = mapped.to_numpy()
    clin["cldn4_source"] = "ORCESTRA"
    clin.loc[mask, "cldn4_source"] = "IMvigor210CoreBiologies_counts_log2p1"

    clin.to_csv(DER / "cldn4_clinical.tsv", sep="\t", index=False)
    clin.to_csv(OUT / "cldn4_clinical.tsv", sep="\t", index=False)

    qc = {
        "n": int(mask.sum()),
        "orcestra_unique_values": int(before.nunique()),
        "orcestra_value": float(before.iloc[0]) if len(before) else None,
        "official_n_matched": n_hit,
        "official_log2_min": float(mapped.min()),
        "official_log2_median": float(mapped.median()),
        "official_log2_max": float(mapped.max()),
        "official_log2_iqr": float(mapped.quantile(0.75) - mapped.quantile(0.25)),
        "official_n_unique": int(mapped.nunique()),
    }
    pd.Series(qc).to_csv(OUT / "imvigor210_cldn4_rescue.tsv", sep="\t", header=["value"])
    print(qc)


if __name__ == "__main__":
    main()
