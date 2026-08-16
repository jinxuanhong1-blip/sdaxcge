#!/usr/bin/env python3
"""Leftover open lung ICI scTCR: GSE179994. Catalog GSE185204 (n=3, skip 900 MB).

GSE179994 GEO tables have scTCR + Tex cluster but no response / MPR column.
Combinatorial expanded vs non-expanded Tex is estimable. Vs response is not,
unless a public response table is added later.
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from stats_util import mw_row, paired_mw, write_json, write_tsv


def parse_series(path: Path) -> pd.DataFrame:
    titles, patients, conditions, gsms = [], [], [], []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("!Sample_geo_accession"):
                gsms = [x.strip().strip('"') for x in line.split("\t")[1:]]
            elif line.startswith("!Sample_title"):
                titles = [x.strip().strip('"') for x in line.split("\t")[1:]]
            elif line.startswith("!Sample_characteristics") and "patient id:" in line:
                patients = [
                    x.strip().strip('"').split("patient id:")[-1].strip()
                    for x in line.split("\t")[1:]
                ]
            elif line.startswith("!Sample_characteristics") and "condition:" in line:
                conditions = [
                    x.strip().strip('"').split("condition:")[-1].strip()
                    for x in line.split("\t")[1:]
                ]
    n = max(len(gsms), len(titles))
    return pd.DataFrame(
        {
            "gsm": gsms[:n] if gsms else [""] * n,
            "title": titles[:n] if titles else [""] * n,
            "patient_geo": patients[:n] if patients else [""] * n,
            "condition": conditions[:n] if conditions else [""] * n,
        }
    )


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--tcr", type=Path, required=True)
    p.add_argument("--meta", type=Path, required=True)
    p.add_argument("--series", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--min-tcr", type=int, default=50)
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    tcr = pd.read_csv(args.tcr, sep="\t", low_memory=False)
    meta = pd.read_csv(args.meta, sep="\t")
    series = parse_series(args.series)
    write_tsv(series, args.out / "gse179994_geo_sample_map.tsv")

    # join on cell id (TCR CellName vs meta cellid are not identical schemes)
    tcr = tcr.copy()
    tcr["clone_id"] = tcr["clone.id"].fillna(
        tcr["CDR3(Alpha1)"].fillna("") + "|" + tcr["CDR3(Beta1)"].fillna("")
    )
    tcr = tcr[tcr["clone_id"].astype(str).str.len() > 1]
    tcr["clone_size"] = tcr.groupby(["sample", "clone_id"], observed=True)["CellName"].transform("size")
    tcr["expanded"] = tcr["clone_size"] >= 2

    # metadata cellid examples: P1.ut.AAACCTGAGGCACATG-1 ; TCR: P1.tr.1.AAACCTGCACATGACT-1
    meta = meta.copy()
    meta["barcode"] = meta["cellid"].astype(str).str.split(".").str[-1]
    tcr["barcode"] = tcr["CellName"].astype(str).str.split(".").str[-1]
    # GEO meta sample is already P1.pre / P1.post.1 — do not strip the capture suffix
    j = tcr.merge(
        meta[["cellid", "patient", "sample", "celltype", "cluster", "barcode"]].rename(
            columns={"sample": "meta_sample"}
        ),
        left_on=["barcode", "sample"],
        right_on=["barcode", "meta_sample"],
        how="left",
    )
    j["is_tex"] = j["cluster"].astype(str).eq("Tex")
    j["is_cd8"] = j["celltype"].astype(str).eq("CD8")
    join_rate = float(j["cluster"].notna().mean())

    rows = []
    for sid, g in j.groupby("sample", observed=True):
        n = int(len(g))
        n_exp = int(g["expanded"].sum())
        n_non = n - n_exp
        n_tex = int(g["is_tex"].sum()) if g["cluster"].notna().any() else 0
        n_exp_tex = int((g["is_tex"] & g["expanded"]).sum()) if g["cluster"].notna().any() else 0
        rows.append(
            {
                "sample": sid,
                "patient": g["patient"].dropna().iloc[0] if g["patient"].notna().any() else sid.split(".")[0],
                "meta_sample": g["meta_sample"].dropna().iloc[0] if g["meta_sample"].notna().any() else "",
                "n_tcr": n,
                "n_joined_meta": int(g["cluster"].notna().sum()),
                "n_clones": int(g["clone_id"].nunique()),
                "frac_exp_cells": n_exp / n if n else np.nan,
                "n_tex": n_tex,
                "frac_tex": n_tex / n if n else np.nan,
                "frac_tex_among_expanded": (n_exp_tex / n_exp) if n_exp else np.nan,
                "frac_tex_among_nonexpanded": (
                    (n_tex - n_exp_tex) / n_non if n_non else np.nan
                ),
                "timepoint": "post" if ".post" in str(sid) else ("pre" if ".pre" in str(sid) else "unknown"),
                "pass_depth": n >= args.min_tcr,
            }
        )
    end = pd.DataFrame(rows)
    write_tsv(end, args.out / "patient_sample_endpoints.tsv")

    work = end[end["pass_depth"]].copy()
    joined_ok = work[work["n_joined_meta"] >= 20].copy()
    stats = []
    stats.append(
        paired_mw(
            joined_ok["frac_tex_among_expanded"],
            joined_ok["frac_tex_among_nonexpanded"],
            "Tex_expanded_vs_nonexpanded",
        )
    )
    stats.append(
        mw_row(
            work.loc[work["timepoint"] == "post", "frac_exp_cells"],
            work.loc[work["timepoint"] == "pre", "frac_exp_cells"],
            "frac_exp_cells",
            "post_vs_pre",
        )
    )
    stats.append(
        mw_row(
            joined_ok.loc[joined_ok["timepoint"] == "post", "frac_tex"],
            joined_ok.loc[joined_ok["timepoint"] == "pre", "frac_tex"],
            "frac_tex",
            "post_vs_pre_joined_only",
        )
    )
    write_tsv(pd.DataFrame(stats), args.out / "stats_combinatorial_and_time.tsv")

    catalog = pd.DataFrame(
        [
            {
                "accession": "GSE179994",
                "paper": "Liu et al. Nat Cancer 2022 PMID 35121991",
                "public_tcr": "GSE179994_all.scTCR.tsv.gz",
                "n_patients_meta": int(meta["patient"].nunique()),
                "n_tcr_cells": int(len(tcr)),
                "response_in_geo": False,
                "mpr_in_geo": False,
                "tacstd2": "not in T-cell tables",
                "use": "leftover open TCR; combinatorial Tex vs expansion; not vs MPR",
            },
            {
                "accession": "GSE185204",
                "paper": "Pai et al. Cancer Cell 2023 PMID 36787796 (lung ICB regions)",
                "public_tcr": "per-region MTX + TCR in GSE185204_RAW.tar (~900 MB)",
                "n_patients_meta": 3,
                "n_tcr_cells": None,
                "response_in_geo": False,
                "mpr_in_geo": False,
                "tacstd2": "possible in GEX MTX; n=3 descriptive only",
                "use": "SKIPPED download (n=3; 900 MB). Listed only.",
            },
            {
                "accession": "GSE236581",
                "paper": "CRC ICB spatiotemporal",
                "public_tcr": "GSE236581_VDJ_merge.txt.gz",
                "n_patients_meta": 22,
                "n_tcr_cells": None,
                "response_in_geo": None,
                "mpr_in_geo": False,
                "tacstd2": "CRC, not lung",
                "use": "not lung; excluded",
            },
            {
                "accession": "GSE162025",
                "paper": "NPC tumour-blood pairs",
                "public_tcr": "filtered_contig_annotations_TCR.csv.gz",
                "n_patients_meta": 10,
                "n_tcr_cells": None,
                "response_in_geo": False,
                "mpr_in_geo": False,
                "tacstd2": "NPC, not lung ICI",
                "use": "not lung ICI; excluded",
            },
        ]
    )
    write_tsv(catalog, args.out / "leftover_tcr_catalog.tsv")
    (args.out / "response_status.txt").write_text(
        "GSE179994 GEO series matrix has patient id and pre/on-treatment only. "
        "No RECIST/MPR column. Paper Supplementary Table 1 has response but is "
        "not redistributed here. Vs-response is not estimable from the public GEO tables.\n"
    )
    spec = {
        "dataset": "GSE179994 leftover + catalog",
        "n_tcr_cells": int(len(tcr)),
        "n_samples": int(end.shape[0]),
        "n_pass_depth": int(end["pass_depth"].sum()),
        "meta_join_rate": join_rate,
        "response_in_geo": False,
        "e2": "not estimable",
        "gse185204": "listed only; n=3; download skipped",
    }
    write_json(spec, args.out / "run_info.json")
    print(f"leftover wrote {args.out}  samples={end.shape[0]} join_rate={join_rate:.3f}")


if __name__ == "__main__":
    main()
