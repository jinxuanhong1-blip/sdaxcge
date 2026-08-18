#!/usr/bin/env python3
"""Mouse-level Cldn4 vs T/NK and IFN/MHC in E-MTAB-5311 (LLC tumors; TISMO-class)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from gene_sets import CLDN4, IFN_MHC, T_NK  # noqa: E402
from score_utils import association_table, mean_score, present_genes, zscore_rows  # noqa: E402

RAW = Path("/tmp/geo_dl")
OUT = Path("/workspace/analysis/E-MTAB-5311")


def load_sdrf(path: Path) -> pd.DataFrame:
    sdrf = pd.read_csv(path, sep="\t")
    # Comment[ENA_RUN] may be present; otherwise map via condensed file.
    cols = {c.lower(): c for c in sdrf.columns}
    run_col = None
    for key in ("comment[ena_run]", "comment[ena run]", "comment[fastqtouri]"):
        if key in cols:
            run_col = cols[key]
            break
    # Find ERR columns by scanning values if needed
    meta = pd.DataFrame()
    meta["source"] = sdrf["Source Name"]
    meta["env"] = sdrf["Characteristics[environmental history]"]
    meta["disease"] = sdrf["Characteristics[disease]"]
    # ENA run is often in a Comment column
    err_series = None
    for c in sdrf.columns:
        vals = sdrf[c].astype(str)
        if vals.str.match(r"ERR\d+").any():
            err_series = vals
            break
    if err_series is None:
        raise RuntimeError("Could not find ERR run accessions in SDRF")
    # If Fastq URI, extract ERR
    extracted = err_series.str.extract(r"(ERR\d+)", expand=False)
    meta["run"] = extracted
    return meta.dropna(subset=["run"]).drop_duplicates("run")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    counts = pd.read_csv(RAW / "E-MTAB-5311-raw-counts.tsv", sep="\t")
    counts = counts.drop_duplicates("Gene Name")
    counts = counts.set_index("Gene Name")
    run_cols = [c for c in counts.columns if c.startswith("ERR")]
    counts = counts[run_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)

    # log2 CPM
    lib = counts.sum(axis=0)
    cpm = counts.div(lib, axis=1) * 1e6
    logcpm = np.log2(cpm + 1.0)

    sdrf = load_sdrf(RAW / "E-MTAB-5311.sdrf.txt")
    sdrf = sdrf.set_index("run").reindex(run_cols)

    present = {
        "Cldn4": present_genes(logcpm.index, CLDN4),
        "T_NK": present_genes(logcpm.index, T_NK),
        "IFN_MHC": present_genes(logcpm.index, IFN_MHC),
    }
    z = zscore_rows(logcpm)
    mouse = pd.DataFrame(index=run_cols)
    mouse["run"] = run_cols
    mouse["source_name"] = sdrf["source"].values
    mouse["exercise"] = sdrf["env"].values
    mouse["model"] = "LLC"
    mouse["cldn4_log2cpm"] = logcpm.loc["Cldn4"] if "Cldn4" in logcpm.index else np.nan
    mouse["cldn4_cpm"] = cpm.loc["Cldn4"] if "Cldn4" in cpm.index else np.nan
    mouse["tnk_score"] = mean_score(z, present["T_NK"])
    mouse["ifn_mhc_score"] = mean_score(z, present["IFN_MHC"])
    mouse["library_size"] = lib.values
    mouse["tnk_genes_used"] = ",".join(present["T_NK"])
    mouse["ifn_mhc_genes_used"] = ",".join(present["IFN_MHC"])
    # One tumor biopsy per mouse in this study (n=16).
    mouse["mouse_id"] = mouse["source_name"].fillna(mouse["run"])

    assoc_all = association_table(mouse, "cldn4_log2cpm", ["tnk_score", "ifn_mhc_score"])
    assoc_all.insert(0, "subset", "all_LLC_mice")
    parts = [assoc_all]
    for env, sub in mouse.groupby("exercise"):
        a = association_table(sub, "cldn4_log2cpm", ["tnk_score", "ifn_mhc_score"])
        a.insert(0, "subset", f"exercise={env}")
        parts.append(a)
    assoc = pd.concat(parts, ignore_index=True)

    genes_used = pd.DataFrame(
        {
            "set": ["Cldn4", "T_NK", "IFN_MHC"],
            "requested": [",".join(CLDN4), ",".join(T_NK), ",".join(IFN_MHC)],
            "present": [
                ",".join(present["Cldn4"]),
                ",".join(present["T_NK"]),
                ",".join(present["IFN_MHC"]),
            ],
            "n_present": [len(present["Cldn4"]), len(present["T_NK"]), len(present["IFN_MHC"])],
        }
    )

    mouse.to_csv(OUT / "mouse_level_scores.tsv", sep="\t", index=False)
    assoc.to_csv(OUT / "high_vs_low_summary.tsv", sep="\t", index=False)
    genes_used.to_csv(OUT / "genes_used.tsv", sep="\t", index=False)

    print(
        f"E-MTAB-5311: n_mice={len(mouse)} Cldn4={bool(present['Cldn4'])} "
        f"T/NK={len(present['T_NK'])} IFN/MHC={len(present['IFN_MHC'])}"
    )
    print(mouse["exercise"].value_counts().to_string())
    print(assoc.to_string(index=False))


if __name__ == "__main__":
    main()
