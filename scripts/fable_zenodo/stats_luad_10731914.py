#!/usr/bin/env python3
"""Honest n / Spearman ρ / p for TACSTD2 vs CLDN4 in Zenodo 10731914.

Uses the already-downloaded human 8-LUAD count matrix. No fabricated stats.
Writes tables under results/noskip/zenodo/.
"""
import csv
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
MAT = ROOT / "results" / "fable_zenodo" / "downloads" / "10731914" / "cell_matrix_sce_8LUADs.csv.gz"
ANN = ROOT / "results" / "fable_zenodo" / "downloads" / "10731914" / "cell_annotation_sce_8LUADs.csv"
OUT = ROOT / "results" / "noskip" / "zenodo"
OUT.mkdir(parents=True, exist_ok=True)

GENES = ["TACSTD2", "CLDN4", "EPCAM"]


def stream_rows(path, wanted):
    found = {}
    with gzip.open(path, "rt") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        barcodes = header[1:]
        for row in reader:
            if row[0] in wanted:
                found[row[0]] = np.array(row[1:], dtype=np.float64)
                if len(found) == len(wanted):
                    break
    return found, barcodes


def spearman_block(x, y, label, extra=None):
    n = int(len(x))
    if n < 3:
        return {"subset": label, "n": n, "rho": None, "p": None, "note": "n<3"}
    rho, p = stats.spearmanr(x, y, nan_policy="omit")
    out = {
        "subset": label,
        "n": n,
        "rho": None if np.isnan(rho) else float(rho),
        "p": None if np.isnan(p) else float(p),
        "method": "scipy.stats.spearmanr two-sided",
    }
    if extra:
        out.update(extra)
    return out


def mw_block(a, b, label):
    n_a, n_b = int(len(a)), int(len(b))
    if n_a < 1 or n_b < 1:
        return {"contrast": label, "n_a": n_a, "n_b": n_b, "U": None, "p": None}
    U, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    return {
        "contrast": label,
        "n_a": n_a,
        "n_b": n_b,
        "median_a": float(np.median(a)),
        "median_b": float(np.median(b)),
        "U": float(U),
        "p": float(p),
        "method": "scipy.stats.mannwhitneyu two-sided",
    }


def main():
    genes, barcodes = stream_rows(MAT, set(GENES))
    missing = [g for g in GENES if g not in genes]
    if missing:
        raise SystemExit(f"missing genes: {missing}")

    ann = pd.read_csv(ANN)
    ann = ann.rename(columns={ann.columns[0]: "barcode"})
    ann["barcode"] = ann["barcode"].astype(str)
    meta = ann.set_index("barcode").loc[list(barcodes)]
    ct = meta["Cell_type.unimodel"].astype(str)
    lib = meta["nCount_RNA"].astype(float).replace(0, np.nan)

    raw = {g: genes[g] for g in GENES}
    cp10k = {g: raw[g] / lib.to_numpy() * 1e4 for g in GENES}

    # Presence
    presence = []
    for g in GENES:
        presence.append({
            "dataset": "zenodo.10731914 human 8-LUAD scRNA-seq",
            "gene": g,
            "present": True,
            "n_cells": int(len(raw[g])),
            "n_expressing": int((raw[g] > 0).sum()),
            "pct_expressing": float((raw[g] > 0).mean() * 100),
            "mean_raw": float(raw[g].mean()),
            "mean_cp10k": float(np.nanmean(cp10k[g])),
        })
    pd.DataFrame(presence).to_csv(OUT / "10731914_gene_presence.csv", index=False)

    # Spearman TACSTD2 vs CLDN4
    x, y = cp10k["TACSTD2"], cp10k["CLDN4"]
    # drop cells with missing libsize
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    ct_ok = ct.to_numpy()[ok]
    epi = ct_ok == "Epithelial cells"
    both_pos = (raw["TACSTD2"][ok] > 0) & (raw["CLDN4"][ok] > 0)

    corrs = [
        spearman_block(x, y, "all_cells_cp10k"),
        spearman_block(x[epi], y[epi], "epithelial_cp10k"),
        spearman_block(x[~epi], y[~epi], "non_epithelial_cp10k"),
        spearman_block(x[both_pos], y[both_pos], "both_genes_detected_cp10k"),
        spearman_block(raw["TACSTD2"][ok], raw["CLDN4"][ok], "all_cells_raw_UMI"),
    ]
    pd.DataFrame(corrs).to_csv(OUT / "10731914_spearman_tacstd2_cldn4.csv", index=False)

    # Mann-Whitney epithelial vs rest
    mw = [
        mw_block(x[epi], x[~epi], "TACSTD2_cp10k epithelial vs rest"),
        mw_block(y[epi], y[~epi], "CLDN4_cp10k epithelial vs rest"),
    ]
    pd.DataFrame(mw).to_csv(OUT / "10731914_mw_epithelial_vs_rest.csv", index=False)

    summary = {
        "doi": "10.5281/zenodo.10731914",
        "n_cells": int(ok.sum()),
        "n_epithelial": int(epi.sum()),
        "n_non_epithelial": int((~epi).sum()),
        "TACSTD2_present": True,
        "CLDN4_present": True,
        "spearman_all_cells": corrs[0],
        "spearman_epithelial": corrs[1],
        "mw_TACSTD2_epi_vs_rest_p": mw[0]["p"],
        "mw_CLDN4_epi_vs_rest_p": mw[1]["p"],
        "note": "ρ is Spearman on CP10K (or raw UMI where labelled). p is two-sided. No multiple-testing correction applied to the five Spearman tests.",
    }
    (OUT / "10731914_stats.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    print("\nSpearman:")
    print(pd.DataFrame(corrs).to_string(index=False))
    print("\nMann-Whitney:")
    print(pd.DataFrame(mw).to_string(index=False))


if __name__ == "__main__":
    main()
