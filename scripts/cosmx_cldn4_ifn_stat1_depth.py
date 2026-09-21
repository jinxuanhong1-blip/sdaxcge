#!/usr/bin/env python3
"""Depth and within-FOV checks for the CosMx CLDN4–IFN/STAT1 correlation.

Complements cosmx_cldn4_ifn_stat1_spatial.py. Same tumor gate, same core
module, same 50-count QC. Asks whether the positive same-cell ρ survives
(1) residualizing on log library size alone and (2) stratification into
within-section count quintiles, and whether it is present inside FOVs.
"""

from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse, stats

import importlib.util

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "cosmx_ifn", ROOT / "scripts" / "cosmx_cldn4_ifn_stat1_spatial.py"
)
base = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(base)

OUT = ROOT / "results" / "cosmx_cldn4_ifn_stat1" / "tables"


def spearman(x, y):
    return base.spearman(np.asarray(x, float), np.asarray(y, float))


def main() -> None:
    meta = base.load_matrix()
    cell_type = meta["cell_type"]
    tumor = np.array([str(x).startswith("tumor") for x in cell_type])
    qc = meta["n_counts"] >= base.MIN_COUNTS
    genes = meta["genes"]
    gix = {g: i for i, g in enumerate(genes)}
    core = meta["inventory"]["ifn_stat1_core_present"]
    core_idx = np.array([gix[g] for g in core])
    rows = []
    gene_rows = []
    for sample in list(dict.fromkeys(meta["sample"].tolist())):
        m = (meta["sample"] == sample) & tumor & qc
        fov = meta["fov"][m]
        counts_fov = pd.Series(fov).value_counts()
        keep_fov = set(counts_fov[counts_fov >= base.MIN_FOV_TUMOR].index.tolist())
        m = m & np.isin(meta["fov"], list(keep_fov))
        ncount = meta["n_counts"][m]
        raw = meta["counts"][m]
        logn = np.log1p(ncount)
        lognorm = np.log1p(raw * (1e4 / np.maximum(ncount, 1.0))[:, None])
        cldn4 = lognorm[:, gix["CLDN4"]]
        detect = (raw[:, core_idx] > 0).mean(axis=0)
        score, _ = base.module_from_lognorm(lognorm[:, core_idx], detect)
        if score is None:
            raise SystemExit(sample)
        patient = str(meta["patient"][m][0])
        fov = meta["fov"][m]
        # Quintiles of library size inside the section.
        q = pd.qcut(ncount, 5, labels=False, duplicates="drop")
        q_rhos = []
        for qi in sorted(set(q)):
            sel = q == qi
            q_rhos.append(spearman(cldn4[sel], score[sel]))
        # Within-FOV Spearman, FOVs with >=40 tumor cells (already true).
        fov_rhos = []
        fov_rhos_lib = []
        for f in np.unique(fov):
            sel = fov == f
            if sel.sum() < 40:
                continue
            fov_rhos.append(spearman(cldn4[sel], score[sel]))
            fov_rhos_lib.append(
                base.residual_spearman(
                    cldn4[sel], score[sel], logn[sel].reshape(-1, 1)
                )
            )
        rows.append(
            {
                "sample": sample,
                "patient": patient,
                "n": int(m.sum()),
                "rho_cldn4_logn": spearman(cldn4, logn),
                "rho_ifn_logn": spearman(score, logn),
                "rho_same": spearman(cldn4, score),
                "rho_same_lib": base.residual_spearman(cldn4, score, logn.reshape(-1, 1)),
                "rho_quintile_mean": float(np.nanmean(q_rhos)),
                "rho_quintiles": ",".join(f"{r:.3f}" for r in q_rhos),
                "rho_within_fov_mean": float(np.nanmean(fov_rhos)),
                "rho_within_fov_lib_mean": float(np.nanmean(fov_rhos_lib)),
                "n_fov": int(len(fov_rhos)),
                "n_fov_neg": int(np.sum(np.array(fov_rhos) < 0)),
            }
        )
        for j, gene in enumerate(core):
            vec = lognorm[:, core_idx[j]]
            gene_rows.append(
                {
                    "sample": sample,
                    "patient": patient,
                    "gene": gene,
                    "rho_lib": base.residual_spearman(cldn4, vec, logn.reshape(-1, 1)),
                }
            )
        print(sample, rows[-1]["rho_same"], rows[-1]["rho_same_lib"], rows[-1]["rho_quintile_mean"], flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "depth_within_fov.tsv", sep="\t", index=False)
    gdf = pd.DataFrame(gene_rows)
    gdf.to_csv(OUT / "gene_spearman_lib.tsv", sep="\t", index=False)
    # Patient means.
    num = [
        "rho_cldn4_logn",
        "rho_ifn_logn",
        "rho_same",
        "rho_same_lib",
        "rho_quintile_mean",
        "rho_within_fov_mean",
        "rho_within_fov_lib_mean",
    ]
    pat = df.groupby("patient")[num].mean()
    summary = {
        "patient_means": {c: {k: float(v) for k, v in pat[c].items()} for c in num},
        "patient_mean_of_means": {c: float(pat[c].mean()) for c in num},
        "patient_n_neg": {c: int((pat[c] < 0).sum()) for c in num},
    }
    gene_pat = gdf.groupby(["patient", "gene"])["rho_lib"].mean().reset_index()
    gene_sum = gene_pat.groupby("gene")["rho_lib"].agg(["mean", lambda s: int((s < 0).sum())])
    gene_sum.columns = ["mean_rho_lib", "n_patients_neg"]
    summary["gene_lib"] = gene_sum.reset_index().to_dict(orient="records")
    (OUT / "depth_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary["patient_mean_of_means"], indent=2))
    print(json.dumps(summary["patient_n_neg"], indent=2))


if __name__ == "__main__":
    main()
