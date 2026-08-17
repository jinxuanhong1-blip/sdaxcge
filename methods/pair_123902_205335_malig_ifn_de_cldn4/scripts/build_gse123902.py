#!/usr/bin/env python3
"""Donor-level marker-malignant UMI-sum for GSE123902 (Laughney 2020).

Locked gates from PR #459:
  malignant = (EPCAM|KRT8|KRT18|KRT19)>0 and PTPRC==0
  tumor/met only; one donor (PRIMARY preferred over METASTASIS); normals dropped.
Patient/donor is the unit. No T/NK matrix. No GSE148071.
"""
from __future__ import annotations

import argparse
import gzip
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
DATA = HERE / "data"
EPI = ["EPCAM", "KRT8", "KRT18", "KRT19"]


def locked_donors() -> pd.DataFrame:
    d = pd.read_csv(DATA / "GSE123902_marker_units.tsv", sep="\t")
    tumor = d[d["tissue"].isin(["PRIMARY", "METASTASIS"])].copy()
    tumor = tumor.sort_values(["patient", "tissue"]).drop_duplicates("patient", keep="first")
    el = tumor[tumor["eligible"].astype(str).str.lower() == "true"].copy()
    return el


def tissue_of(name: str) -> str:
    if "NORMAL" in name:
        return "NORMAL"
    if "METASTASIS" in name:
        return "METASTASIS"
    if "PRIMARY" in name:
        return "PRIMARY"
    return "OTHER"


def marker_mal(df: pd.DataFrame) -> np.ndarray:
    cols = {c.upper(): c for c in df.columns}
    epi = np.zeros(len(df), dtype=bool)
    for g in EPI:
        if g in cols:
            epi |= df[cols[g]].to_numpy(dtype=float) > 0
    ptprc = df[cols["PTPRC"]].to_numpy(dtype=float) if "PTPRC" in cols else np.zeros(len(df))
    return epi & (ptprc == 0)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--tar",
        type=Path,
        default=Path("/tmp/pair_123902_205335_malig/GSE123902/GSE123902_RAW.tar"),
    )
    args = ap.parse_args()
    if not args.tar.exists():
        raise SystemExit(f"missing {args.tar}; run scripts/download_gse123902.py")

    locked = locked_donors()
    keep = set(locked["patient"].astype(str))
    print(f"locked donors n={len(keep)}: {sorted(keep)}", flush=True)

    members = []
    with tarfile.open(args.tar) as tf:
        for m in tf.getmembers():
            if not m.name.endswith(".csv.gz"):
                continue
            name = Path(m.name).name
            parts = name.split("_")
            donor = parts[2] if len(parts) > 2 else name
            tissue = tissue_of(name)
            if tissue not in {"PRIMARY", "METASTASIS"}:
                continue
            if donor not in keep:
                continue
            members.append((donor, tissue, m, name))
        members.sort(key=lambda x: (x[0], 0 if x[1] == "PRIMARY" else 1))

        sums: dict[str, pd.Series] = {}
        meta_rows = []
        chosen: dict[str, str] = {}
        for donor, tissue, m, name in members:
            if donor in chosen:
                continue
            chosen[donor] = tissue
            raw = gzip.GzipFile(fileobj=tf.extractfile(m))
            df = pd.read_csv(raw, index_col=0)
            df.columns = [str(c).upper() for c in df.columns]
            mal = marker_mal(df)
            n_mal = int(mal.sum())
            if n_mal < 1:
                print(f"  SKIP {name}: 0 marker-malignant", flush=True)
                continue
            umi = df.loc[mal].sum(axis=0).astype(float)
            umi.index = umi.index.astype(str)
            sums[donor] = umi
            rec = locked.loc[locked["patient"].astype(str) == donor].iloc[0]
            meta_rows.append(
                {
                    "patient": donor,
                    "cohort": "GSE123902",
                    "unit": "donor",
                    "tissue": tissue,
                    "file": name,
                    "n_cells": int(len(df)),
                    "n_malignant": n_mal,
                    "n_malignant_locked": int(rec["n_malignant"]),
                    "mal_CLDN4_pct": float(rec["mal_CLDN4_pct"]),
                    "mal_CLDN4_mean": float(rec["mal_CLDN4_mean"]),
                    "malignant_rule": "marker_malig EPCAM|KRT8|KRT18|KRT19>0 and PTPRC==0; tumor/met only",
                }
            )
            print(
                f"  {name}: cells={len(df)} mal={n_mal} locked_mal={int(rec['n_malignant'])} genes={umi.size}",
                flush=True,
            )

    if not sums:
        raise SystemExit("no GSE123902 donor sums")
    genes = sorted(set().union(*[set(s.index) for s in sums.values()]))
    mat = pd.DataFrame({p: s.reindex(genes).fillna(0.0) for p, s in sums.items()}, index=genes)
    mat.index.name = ""
    DATA.mkdir(parents=True, exist_ok=True)
    out_c = DATA / "GSE123902_malignant_counts.tsv.gz"
    out_m = DATA / "GSE123902_malignant_meta.tsv"
    mat.to_csv(out_c, sep="\t")
    pd.DataFrame(meta_rows).to_csv(out_m, sep="\t", index=False)
    print(f"wrote {out_c} genes={mat.shape[0]} donors={mat.shape[1]}", flush=True)
    print(f"wrote {out_m}", flush=True)


if __name__ == "__main__":
    main()
