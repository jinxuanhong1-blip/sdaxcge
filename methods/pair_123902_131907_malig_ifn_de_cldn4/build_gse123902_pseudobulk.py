#!/usr/bin/env python3
"""Marker-malignant UMI-sum per GSE123902 donor (Laughney 2020).

Gate (locked, PR #459): (EPCAM|KRT8|KRT18|KRT19)>0 and PTPRC==0.
Tumor/met only. Patient/donor is the unit. Normals dropped.
"""
from __future__ import annotations

import gzip
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
TAR = Path("/tmp/geo_dl/GSE123902_RAW.tar")
EPI = ["EPCAM", "KRT8", "KRT18", "KRT19"]


def _is_tumor(name: str) -> bool:
    return ("NORMAL" not in name) and ("PRIMARY" in name or "METASTASIS" in name)


def build(tar_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    units = pd.read_csv(DATA / "GSE123902_marker_units.tsv", sep="\t")
    tumor = units[units["tissue"].isin(["PRIMARY", "METASTASIS"])].copy()
    tumor = tumor.sort_values(["patient", "tissue"]).drop_duplicates("patient", keep="first")
    elig = tumor[tumor["eligible"].astype(str).str.lower() == "true"]
    keep_donors = set(elig["patient"].astype(str))

    sums: dict[str, pd.Series] = {}
    audit = []

    with tarfile.open(tar_path) as tf:
        for m in tf.getmembers():
            if not m.name.endswith(".csv.gz"):
                continue
            name = Path(m.name).name
            if not _is_tumor(name):
                continue
            parts = name.split("_")
            donor = parts[2] if len(parts) > 2 else name
            if donor not in keep_donors:
                continue
            raw = gzip.GzipFile(fileobj=tf.extractfile(m))
            header = raw.readline().decode().strip().split(",")
            genes = [g.strip().upper() for g in header[1:]]
            idx = {g: i for i, g in enumerate(genes)}
            epi_i = [idx[g] for g in EPI if g in idx]
            ptprc_i = idx.get("PTPRC")
            acc = np.zeros(len(genes), dtype=np.float64)
            n_mal = 0
            n_cells = 0
            n_genes = len(genes)
            for line in raw:
                bits = line.decode().strip().split(",")
                vec = np.fromiter((float(x) for x in bits[1:]), dtype=np.float64, count=n_genes)
                n_cells += 1
                epi = False
                for i in epi_i:
                    if vec[i] > 0:
                        epi = True
                        break
                if not epi:
                    continue
                if ptprc_i is not None and vec[ptprc_i] != 0:
                    continue
                acc += vec
                n_mal += 1
            series = pd.Series(acc, index=genes).groupby(level=0).sum()
            sums[donor] = series if donor not in sums else sums[donor].add(series, fill_value=0.0)
            audit.append(
                {
                    "file": name,
                    "patient": donor,
                    "n_cells": n_cells,
                    "n_malignant_built": n_mal,
                    "n_genes": int(len(series)),
                    "libsize": float(acc.sum()),
                }
            )
            print(f"  {name}: cells={n_cells} mal={n_mal} genes={len(series)} lib={acc.sum():.0f}", flush=True)

    if not sums:
        raise SystemExit("no GSE123902 tumor matrices scored")
    mat = pd.DataFrame(sums).fillna(0.0)
    mat.index = mat.index.astype(str)
    mat = mat.groupby(mat.index).sum()
    mat = mat.reindex(sorted(mat.columns), axis=1)
    return mat, pd.DataFrame(audit)


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    if not TAR.exists():
        raise SystemExit(f"missing {TAR}; run download.sh")
    print("Building GSE123902 marker-malignant donor UMI-sum")
    mat, audit = build(TAR)
    out = DATA / "GSE123902_malignant_counts.tsv.gz"
    mat.to_csv(out, sep="\t", compression="gzip")
    audit.to_csv(DATA / "GSE123902_build_audit.tsv", sep="\t", index=False)
    print(f"wrote {out}  genes={mat.shape[0]} donors={mat.shape[1]}")
    print(audit.to_string(index=False))


if __name__ == "__main__":
    main()
