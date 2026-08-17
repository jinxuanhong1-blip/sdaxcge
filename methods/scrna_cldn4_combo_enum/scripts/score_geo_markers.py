#!/usr/bin/env python3
"""Marker-malignant CLDN4 vs T/NK on GEO processed matrices <2GB.

GSE123902: 17 dense CSV (Laughney 2020). Tumor samples only.
GSE189357: 9× 10x MTX (Zhu/Wang AIS–IAC atlas). One sample per patient.

Locked marker gates (same for both):
  malignant-like = (EPCAM>0 or KRT8>0 or KRT18>0 or KRT19>0) and PTPRC==0
  T/NK           = (CD3D>0 or CD3E>0 or CD8A>0 or NKG7>0 or GNLY>0 or KLRD1>0)
                   and not malignant-like
CLDN4 mean = mean log1p(count) in malignant-like cells; %pos = fraction count>0.
Eligible unit: n_mal>=20 and n_tnk>=20. Patient/donor is the unit.
"""
from __future__ import annotations

import gzip
import io
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
DATA = HERE / "data"
CACHE = Path("/tmp/geo_dl")

EPI = ["EPCAM", "KRT8", "KRT18", "KRT19"]
TNK = ["CD3D", "CD3E", "CD8A", "NKG7", "GNLY", "KLRD1"]
PANEL = ["CLDN4", "PTPRC"] + EPI + TNK


def _upper_map(names: list[str]) -> dict[str, int]:
    return {n.upper(): i for i, n in enumerate(names)}


def _score_counts(counts: np.ndarray, idx: dict[str, int]) -> dict:
    """counts: cells x genes (raw UMI)."""
    def col(g):
        if g not in idx:
            return np.zeros(counts.shape[0], dtype=float)
        return counts[:, idx[g]].astype(float)

    epi = np.zeros(counts.shape[0], dtype=bool)
    for g in EPI:
        epi |= col(g) > 0
    ptprc = col("PTPRC")
    mal = epi & (ptprc == 0)
    tnk = np.zeros(counts.shape[0], dtype=bool)
    for g in TNK:
        tnk |= col(g) > 0
    tnk = tnk & (~mal)
    cldn4 = col("CLDN4")
    n_mal = int(mal.sum())
    n_tnk = int(tnk.sum())
    if n_mal == 0:
        mean = pct = float("nan")
    else:
        x = np.log1p(cldn4[mal])
        mean = float(x.mean())
        pct = float((cldn4[mal] > 0).mean())
    return {
        "n_cells": int(counts.shape[0]),
        "n_malignant": n_mal,
        "n_tnk": n_tnk,
        "frac_tnk": n_tnk / counts.shape[0] if counts.shape[0] else float("nan"),
        "mal_CLDN4_mean": mean,
        "mal_CLDN4_pct": pct,
        "eligible": n_mal >= 20 and n_tnk >= 20,
    }


def score_gse123902(tar_path: Path) -> pd.DataFrame:
    rows = []
    with tarfile.open(tar_path) as tf:
        for m in tf.getmembers():
            if not m.name.endswith(".csv.gz"):
                continue
            name = Path(m.name).name
            # GSM3516662_MSK_LX653_PRIMARY_TUMOUR_dense.csv.gz
            parts = name.split("_")
            gsm = parts[0]
            donor = parts[2] if len(parts) > 2 else name
            tissue = "NORMAL" if "NORMAL" in name else (
                "METASTASIS" if "METASTASIS" in name else (
                    "PRIMARY" if "PRIMARY" in name else "OTHER"
                )
            )
            raw = gzip.GzipFile(fileobj=tf.extractfile(m))
            # cells x genes; header gene names; col0 cell id
            header = raw.readline().decode().strip().split(",")
            genes = [g.upper() for g in header[1:]]
            imap = _upper_map(genes)
            keep = [imap[g] for g in PANEL if g in imap]
            keep_names = [g for g in PANEL if g in imap]
            cols = []
            n = 0
            for line in raw:
                bits = line.decode().strip().split(",")
                vec = np.array([float(bits[j + 1]) for j in keep], dtype=float)
                cols.append(vec)
                n += 1
            counts = np.vstack(cols) if cols else np.zeros((0, len(keep)))
            idx = {g: i for i, g in enumerate(keep_names)}
            rec = _score_counts(counts, idx)
            rec.update(
                {
                    "cohort": "GSE123902",
                    "unit_id": donor,
                    "patient": donor,
                    "gsm": gsm,
                    "tissue": tissue,
                    "file": name,
                    "malig_def": "marker_malig",
                    "unit": "donor",
                }
            )
            rows.append(rec)
            print(f"  {name}: n={rec['n_cells']} mal={rec['n_malignant']} tnk={rec['n_tnk']} elig={rec['eligible']}", flush=True)
    return pd.DataFrame(rows)


def _read_10x_features(tf: tarfile.TarFile, members: dict, sample: str) -> list[str]:
    name = next(n for n in members if f"_{sample}_features" in n or f"_{sample}_genes" in n)
    genes = []
    with gzip.GzipFile(fileobj=tf.extractfile(members[name])) as f:
        for line in f:
            p = line.decode().strip().split("\t")
            genes.append((p[1] if len(p) > 1 else p[0]).upper())
    return genes


def _read_10x_nbarcodes(tf: tarfile.TarFile, members: dict, sample: str) -> int:
    name = next(n for n in members if f"_{sample}_barcodes" in n)
    n = 0
    with gzip.GzipFile(fileobj=tf.extractfile(members[name])) as f:
        for _ in f:
            n += 1
    return n


def score_gse189357(tar_path: Path) -> pd.DataFrame:
    rows = []
    with tarfile.open(tar_path) as tf:
        members = {m.name: m for m in tf.getmembers()}
        samples = [f"TD{i}" for i in range(1, 10)]
        for sample in samples:
            genes = _read_10x_features(tf, members, sample)
            imap = _upper_map(genes)
            n_cells = _read_10x_nbarcodes(tf, members, sample)
            keep_idx = {imap[g] for g in PANEL if g in imap}  # 0-based gene
            # accumulate cells x panel
            keep_names = [g for g in PANEL if g in imap]
            gene_to_col = {g: i for i, g in enumerate(keep_names)}
            counts = np.zeros((n_cells, len(keep_names)), dtype=float)
            mtx_name = next(n for n in members if f"_{sample}_matrix.mtx" in n)
            with gzip.GzipFile(fileobj=tf.extractfile(members[mtx_name])) as f:
                # skip comments + size line
                for line in f:
                    s = line.decode()
                    if s.startswith("%"):
                        continue
                    # dims
                    break
                for line in f:
                    a, b, v = line.decode().split()
                    gi = int(a) - 1
                    ci = int(b) - 1
                    if gi in keep_idx:
                        g = genes[gi]
                        counts[ci, gene_to_col[g]] = float(v)
            rec = _score_counts(counts, gene_to_col)
            rec.update(
                {
                    "cohort": "GSE189357",
                    "unit_id": sample,
                    "patient": sample,
                    "tissue": "TUMOR",
                    "malig_def": "marker_malig",
                    "unit": "patient",
                }
            )
            rows.append(rec)
            print(f"  {sample}: n={rec['n_cells']} mal={rec['n_malignant']} tnk={rec['n_tnk']} elig={rec['eligible']}", flush=True)
    return pd.DataFrame(rows)


def main():
    DATA.mkdir(parents=True, exist_ok=True)
    p123 = CACHE / "GSE123902_RAW.tar"
    p189 = CACHE / "GSE189357_RAW.tar"
    if p123.exists():
        print("Scoring GSE123902")
        d = score_gse123902(p123)
        d.to_csv(DATA / "GSE123902_marker_units.tsv", sep="\t", index=False)
        print(d[["unit_id", "tissue", "n_cells", "n_malignant", "n_tnk", "eligible", "mal_CLDN4_mean"]].to_string())
    else:
        print("SKIP GSE123902: tar missing")
    if p189.exists():
        print("Scoring GSE189357")
        d = score_gse189357(p189)
        d.to_csv(DATA / "GSE189357_marker_units.tsv", sep="\t", index=False)
        print(d[["unit_id", "n_cells", "n_malignant", "n_tnk", "eligible", "mal_CLDN4_mean"]].to_string())
    else:
        print("SKIP GSE189357: tar missing")


if __name__ == "__main__":
    main()
