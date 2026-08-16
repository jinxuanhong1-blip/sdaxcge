#!/usr/bin/env python3
"""Score E-MTAB-13526 tumor h5ad remotely: author cell types + TACSTD2/CLDN4 columns.

Streams CSR X over HTTPS and keeps only TACSTD2 and CLDN4. Public ArrayExpress only.
"""
from __future__ import annotations

import json
from pathlib import Path

import fsspec
import h5py
import numpy as np
import pandas as pd
from scipy import stats

OUT = Path(__file__).resolve().parents[1] / "tables"
META = Path(__file__).resolve().parents[1] / "metadata"
OUT.mkdir(parents=True, exist_ok=True)
URL = "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/526/E-MTAB-13526/Files/10X_Lung_Tumour_Annotated_v2.h5ad"

EPI_KEYS = ("epithelial", "tumour", "tumor", "at2", "at1", "basal", "club", "goblet", "ciliated", "ionocyte", "tuft", "malignant")
TNK_KEYS = ("t cell", "t-cell", "cd4", "cd8", "nk", "ilc", "nkt")


def decode(ds):
    x = ds[...]
    if getattr(x, "dtype", None) is not None and x.dtype.kind in {"S", "O"}:
        return np.array([i.decode() if isinstance(i, bytes) else str(i) for i in x])
    return x


def cat_col(obs, name):
    codes = np.asarray(obs[name][...])
    cats = decode(obs["__categories"][name])
    return np.array([str(cats[int(c)]) if 0 <= int(c) < len(cats) else "NA" for c in codes])


def is_epi(label: str) -> bool:
    s = label.lower()
    return any(k in s for k in EPI_KEYS)


def is_tnk(label: str) -> bool:
    s = label.lower()
    if "macrophage" in s or "monocyte" in s or "dendritic" in s:
        return False
    return any(k in s for k in TNK_KEYS)


def main():
    fs = fsspec.filesystem("https")
    f = fs.open(URL, "rb", block_size=16 * 2**20)
    h = h5py.File(f, "r")
    obs = h["obs"]
    patient = cat_col(obs, "patient")
    ctype = cat_col(obs, "Cell types")
    sorting = cat_col(obs, "sorting") if "sorting" in obs else np.array(["NA"] * len(patient))
    ttype = cat_col(obs, "tumour type") if "tumour type" in obs else np.array(["NA"] * len(patient))
    sample = cat_col(obs, "patient_sample") if "patient_sample" in obs else cat_col(obs, "sample")
    n = len(patient)
    print("n_cells", n, "n_patients", len(set(patient)))
    print("cell types", sorted(set(ctype))[:40], "...", len(set(ctype)))
    print("sorting", sorted(set(sorting)))
    print("tumour type", sorted(set(ttype)))

    var_idx = decode(h["var"]["_index"])
    want = {"TACSTD2": None, "CLDN4": None}
    for i, g in enumerate(var_idx):
        if g in want:
            want[g] = i
    print("gene cols", want)
    if any(v is None for v in want.values()):
        raise SystemExit("missing genes")

    X = h["X"]
    indptr = np.asarray(X["indptr"][...])  # n+1
    tac = np.zeros(n, dtype=np.float32)
    cld = np.zeros(n, dtype=np.float32)
    nnz = int(X["data"].shape[0])
    chunk = 8_000_000
    # map nnz offset -> cell via searchsorted on indptr
    for start in range(0, nnz, chunk):
        stop = min(start + chunk, nnz)
        idx = np.asarray(X["indices"][start:stop])
        dat = np.asarray(X["data"][start:stop])
        # cell id for each nnz
        cells = np.searchsorted(indptr, np.arange(start, stop, dtype=np.int64), side="right") - 1
        m_t = idx == want["TACSTD2"]
        m_c = idx == want["CLDN4"]
        if m_t.any():
            tac[cells[m_t]] = dat[m_t]
        if m_c.any():
            cld[cells[m_c]] = dat[m_c]
        print(f"stream {stop}/{nnz} ({100*stop/nnz:.1f}%)", flush=True)

    epi = np.array([is_epi(x) for x in ctype])
    tnk = np.array([is_tnk(x) for x in ctype])
    print("n_epi", int(epi.sum()), "n_tnk", int(tnk.sum()))
    cells = pd.DataFrame(
        {
            "patient": patient,
            "cell_type": ctype,
            "sorting": np.char.strip(sorting.astype(str)),
            "tumour_type": ttype,
            "sample": sample,
            "epi": epi,
            "tnk": tnk,
            "TACSTD2": tac,
            "CLDN4": cld,
        }
    )
    cells.to_csv("/tmp/scrna_leftover_2426/EMTAB13526_per_cell_genes.tsv.gz", sep="\t", index=False)

    def agg(d, subset):
        rows = []
        for p, g in d.groupby("patient"):
            e = g[g["epi"]]
            t = g[g["tnk"]]
            rows.append(
                {
                    "dataset": "E-MTAB-13526",
                    "subset": subset,
                    "patient": p,
                    "n_cells": int(len(g)),
                    "n_epi": int(len(e)),
                    "n_tnk": int(len(t)),
                    "frac_tnk": float(g["tnk"].mean()) if len(g) else np.nan,
                    "frac_epi": float(g["epi"].mean()) if len(g) else np.nan,
                    "epi_TACSTD2": float(e["TACSTD2"].mean()) if len(e) else np.nan,
                    "epi_CLDN4": float(e["CLDN4"].mean()) if len(e) else np.nan,
                    "tumour_type": ",".join(sorted(set(g["tumour_type"]))),
                    "sorting": ",".join(sorted(set(g["sorting"]))),
                    "samples": ",".join(sorted(set(g["sample"]))),
                }
            )
        return pd.DataFrame(rows)

    df_all = agg(cells, "all_sorts")
    df_cd = agg(cells[cells["sorting"] == "CD235a-"], "CD235a_minus")
    df = pd.concat([df_all, df_cd], ignore_index=True)
    df.to_csv(OUT / "EMTAB13526_per_patient.tsv", sep="\t", index=False)

    # cell-type inventory
    ct = (
        pd.DataFrame({"cell_type": ctype})
        .value_counts()
        .rename("n")
        .reset_index()
    )
    ct.to_csv(OUT / "EMTAB13526_celltype_counts.tsv", sep="\t", index=False)

    stats_rows = []
    for subset, primary in [("CD235a_minus", True), ("all_sorts", False)]:
        elig = df[(df["subset"] == subset) & (df["n_epi"] >= 20) & (df["n_tnk"] >= 20)]
        for gene in ["epi_TACSTD2", "epi_CLDN4"]:
            d = elig.dropna(subset=[gene, "frac_tnk"])
            if len(d) >= 4:
                rho, p = stats.spearmanr(d[gene], d["frac_tnk"])
            else:
                rho, p = np.nan, np.nan
            stats_rows.append(
                {
                    "dataset": "E-MTAB-13526",
                    "analysis": "spearman_patient_" + subset,
                    "comparison": f"{gene} vs frac_tnk",
                    "metric": gene,
                    "n": int(len(d)),
                    "stat": "spearman_rho",
                    "value": float(rho) if rho == rho else np.nan,
                    "p_value": float(p) if p == p else np.nan,
                    "note": (
                        "Cvejic/De Zuani NSCLC atlas, tumor h5ad, author Cell types; not ICI; extra n. "
                        + ("PRIMARY: CD235a- lanes only (has epithelium + immune). " if primary else "SENSITIVITY: all sorts including myeloid/MDSC. ")
                        + "X streamed remotely (author-normalized CSR)."
                    ),
                }
            )
    pd.DataFrame(stats_rows).to_csv(OUT / "EMTAB13526_stats.tsv", sep="\t", index=False)
    print(df.to_string(index=False))
    print(pd.DataFrame(stats_rows).to_string(index=False))
    # save a small audit
    (META / "E-MTAB-13526_extract_audit.json").write_text(
        json.dumps(
            {
                "n_cells": n,
                "n_patients": int(df.shape[0]),
                "n_epi": int(epi.sum()),
                "n_tnk": int(tnk.sum()),
                "cell_types": sorted(set(ctype)),
                "url": URL,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
