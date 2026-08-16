#!/usr/bin/env python3
"""Recompute GSE243013 patient-level TACSTD2/CLDN4 leak vs MPR / RECIST.

Expression: two-gene nonzero entries extracted from the public immune MTX
(GSE243013_NSCLC_immune_scRNA_counts.mtx.gz; 1,254,749 cells × 31,831 genes).
Labels and library sizes: independently downloaded GEO metadata CSV.

This matrix is the deposited CD45+ immune atlas after neoadjuvant chemo-IO.
It is not malignant-cell RNA. Residual tumor / ambient epithelial transcripts
are a live alternative explanation.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd

from config import GSE243013_EXTRACT, GSE243013_GENES, GSE243013_META, RESULTS
from stats_lib import two_group


def _load_meta() -> pd.DataFrame:
    meta = pd.read_csv(GSE243013_META)
    if len(meta) != 1_254_749:
        raise SystemExit(f"unexpected metadata rows: {len(meta)}")
    meta = meta.reset_index(drop=True)
    meta["cell_index_1based"] = np.arange(1, len(meta) + 1)
    return meta


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    genes = pd.read_csv(GSE243013_GENES)
    gene_col = genes.columns[0]
    symbols = genes[gene_col].astype(str).tolist()
    for g in ("TACSTD2", "CLDN4"):
        if symbols.count(g) != 1:
            raise SystemExit(f"{g} count={symbols.count(g)}")

    meta = _load_meta()
    ext = pd.read_csv(GSE243013_EXTRACT, sep="\t")
    if ext["cell_index_1based"].max() > len(meta):
        raise SystemExit("extract index exceeds metadata")

    wide = (
        ext.pivot_table(index="cell_index_1based", columns="gene", values="count", aggfunc="sum")
        .reindex(columns=["TACSTD2", "CLDN4"])
        .fillna(0)
        .astype(int)
    )
    meta = meta.merge(wide, how="left", left_on="cell_index_1based", right_index=True)
    meta[["TACSTD2", "CLDN4"]] = meta[["TACSTD2", "CLDN4"]].fillna(0).astype(int)

    def _agg(g: pd.DataFrame) -> pd.Series:
        n = len(g)
        lib = g["total_counts"].to_numpy(dtype=float)
        out = {
            "n_cells": n,
            "sum_total_counts": float(lib.sum()),
            "pathological_response": g["pathological_response"].iloc[0],
            "pathological_response_rate": g["pathological_response_rate"].iloc[0],
            "radiological_response": g["radiological_response"].iloc[0],
            "cancer_type": g["cancer_type"].iloc[0],
            "gender": g["gender"].iloc[0],
            "age": g["age"].iloc[0],
            "anti_PD1_therapy": g["anti-PD1_therapy"].iloc[0],
            "chemotherapy": g["chemotherapy"].iloc[0],
        }
        for gene in ("TACSTD2", "CLDN4"):
            vec = g[gene].to_numpy(dtype=float)
            out[f"{gene}_n_pos"] = int(np.sum(vec > 0))
            out[f"{gene}_frac_pos"] = float(np.mean(vec > 0))
            out[f"{gene}_sum"] = int(vec.sum())
            out[f"{gene}_mean_count"] = float(vec.mean())
            out[f"{gene}_mean_cpm"] = float(np.mean(np.where(lib > 0, 1e6 * vec / lib, 0.0)))
        return pd.Series(out)

    patient = (
        meta.groupby("sampleID", sort=False)
        .apply(_agg)
        .reset_index()
    )
    if "level_1" in patient.columns:
        patient = patient.drop(columns=["level_1"])
    patient["dataset"] = "GSE243013"
    patient["compartment"] = "CD45+_deposited_immune_atlas"
    patient["not_malignant_rna"] = True
    patient["mpr_any"] = np.where(
        patient["pathological_response"].isin(["pCR", "MPR"]), "yes",
        np.where(patient["pathological_response"] == "non-MPR", "no", "unknown"),
    )
    # RECIST-like deposited labels
    recist_map = {
        "CR": "CR", "PR": "PR", "SD": "SD", "PD": "PD",
        "pCR": "CR",  # do not invent; only if deposited this way
    }
    patient["recist"] = patient["radiological_response"].astype(str)
    patient_path = RESULTS / "gse243013_patient_level.tsv"
    patient.to_csv(patient_path, sep="\t", index=False)

    analyzed = patient[patient["mpr_any"].isin(["yes", "no"])].copy()
    tests = []
    for gene in ("TACSTD2", "CLDN4"):
        metric = f"{gene}_frac_pos"
        a = analyzed.loc[analyzed["mpr_any"] == "yes", metric].to_numpy()
        b = analyzed.loc[analyzed["mpr_any"] == "no", metric].to_numpy()
        t = two_group(a, b, "MPR-any", "non-MPR")
        t.update({
            "dataset": "GSE243013",
            "gene": gene,
            "metric": metric,
            "contrast": "mpr_any_vs_nonMPR",
            "note": "GEO pCR+MPR vs non-MPR; exclude unknowm",
        })
        tests.append(t)

        # RECIST: CR+PR vs SD+PD if those strings exist
        rec = analyzed.copy()
        rec["recist_bin"] = rec["radiological_response"].astype(str).str.upper()
        responder = rec["recist_bin"].isin(["CR", "PR", "OR", "RESPONSE"])
        nonresp = rec["recist_bin"].isin(["SD", "PD", "NR", "NON-RESPONSE", "NONRESPONSE"])
        if responder.sum() >= 2 and nonresp.sum() >= 2:
            t2 = two_group(rec.loc[responder, metric], rec.loc[nonresp, metric], "CR/PR", "SD/PD")
            t2.update({
                "dataset": "GSE243013",
                "gene": gene,
                "metric": metric,
                "contrast": "RECIST_CRPR_vs_SDPD",
                "note": f"deposited radiological_response; n_resp={int(responder.sum())} n_nr={int(nonresp.sum())}",
            })
            tests.append(t2)

    tests_df = pd.DataFrame(tests)
    tests_df.to_csv(RESULTS / "gse243013_stats.tsv", sep="\t", index=False)

    recist_counts = analyzed["radiological_response"].value_counts(dropna=False).to_dict()
    mpr_counts = analyzed["pathological_response"].value_counts(dropna=False).to_dict()
    summary = {
        "dataset": "GSE243013",
        "compartment": "CD45+ deposited immune scRNA (not malignant RNA)",
        "n_cells_metadata": int(len(meta)),
        "n_samples_geo": int(patient["sampleID"].nunique()),
        "n_analyzed_mpr_known": int(len(analyzed)),
        "n_mpr_any": int((analyzed["mpr_any"] == "yes").sum()),
        "n_nonmpr": int((analyzed["mpr_any"] == "no").sum()),
        "n_unknown_mpr": int((patient["mpr_any"] == "unknown").sum()),
        "pathological_response_counts": {str(k): int(v) for k, v in mpr_counts.items()},
        "radiological_response_counts": {str(k): int(v) for k, v in recist_counts.items()},
        "tacstd2_cells_pos": int((meta["TACSTD2"] > 0).sum()),
        "cldn4_cells_pos": int((meta["CLDN4"] > 0).sum()),
        "extract_nonzero_rows": int(len(ext)),
        "mtx_url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE243nnn/GSE243013/suppl/GSE243013_NSCLC_immune_scRNA_counts.mtx.gz",
        "tests": tests_df.to_dict(orient="records"),
    }
    (RESULTS / "gse243013_summary.json").write_text(json.dumps(summary, indent=2))
    print("n_analyzed", len(analyzed), "mpr_any", (analyzed.mpr_any == "yes").sum(),
          "nonMPR", (analyzed.mpr_any == "no").sum())
    print("RECIST labels:", recist_counts)
    print(tests_df[["gene", "contrast", "n_a", "n_b", "g_median_a", "g_median_b", "mw_p", "g_hedges_g"]].to_string(index=False))


if __name__ == "__main__":
    main()
