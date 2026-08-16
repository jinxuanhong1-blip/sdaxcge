#!/usr/bin/env python3
"""Probe E-MTAB-13526 tumor h5ad over HTTPS (obs/var/X encoding only)."""
from __future__ import annotations

import json
from pathlib import Path

import fsspec
import h5py
import numpy as np

OUT = Path(__file__).resolve().parents[1] / "metadata"
URL = "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/526/E-MTAB-13526/Files/10X_Lung_Tumour_Annotated_v2.h5ad"


def decode(ds):
    if hasattr(ds, "asstr"):
        try:
            return ds.asstr()[...]
        except Exception:
            pass
    x = ds[...]
    if getattr(x, "dtype", None) is not None and x.dtype.kind in {"S", "O"}:
        return np.array([i.decode() if isinstance(i, bytes) else str(i) for i in x])
    return x


def main():
    fs = fsspec.filesystem("https")
    f = fs.open(URL, "rb", block_size=8 * 2**20)
    h = h5py.File(f, "r")
    info = {"keys": list(h.keys()), "obs": list(h["obs"].keys()), "var": list(h["var"].keys())}
    x = h["X"]
    info["X_keys"] = list(x.keys())
    for k in x.keys():
        d = x[k]
        info[f"X_{k}_shape"] = list(d.shape)
        info[f"X_{k}_dtype"] = str(d.dtype)
    # encoding
    info["X_attrs"] = {k: str(v) for k, v in x.attrs.items()}
    # var index
    var_idx = decode(h["var"]["_index"])
    info["n_var"] = int(len(var_idx))
    info["var_example"] = [str(v) for v in var_idx[:8]]
    wanted = ["TACSTD2", "CLDN4", "EPCAM", "CD3D", "NKG7"]
    info["gene_hits"] = {g: int(np.where(var_idx == g)[0][0]) for g in wanted if g in set(var_idx)}
    # obs useful columns
    for col in [
        "Cell types",
        "Cell types v25",
        "batch",
        "environment",
        "exp",
        "cancer stage",
        "n_counts",
        "patient",
        "sample",
        "donor",
    ]:
        if col in h["obs"]:
            ds = h["obs"][col]
            info[f"obs_{col}_shape"] = list(ds.shape)
            info[f"obs_{col}_dtype"] = str(ds.dtype)
            try:
                vals = decode(ds)
                if getattr(vals, "dtype", None) is not None and vals.dtype.kind in {"U", "O"}:
                    uniq, cnt = np.unique(vals, return_counts=True)
                    info[f"obs_{col}_nunique"] = int(len(uniq))
                    info[f"obs_{col}_top"] = {str(a): int(b) for a, b in zip(uniq[np.argsort(-cnt)][:20], cnt[np.argsort(-cnt)][:20])}
                else:
                    info[f"obs_{col}_min"] = float(np.min(vals))
                    info[f"obs_{col}_max"] = float(np.max(vals))
            except Exception as e:
                info[f"obs_{col}_err"] = str(e)
    # categories
    if "__categories" in h["obs"]:
        info["obs_categories"] = list(h["obs"]["__categories"].keys())
    (OUT / "E-MTAB-13526_h5ad_probe.json").write_text(json.dumps(info, indent=2, default=str))
    print(json.dumps(info, indent=2, default=str)[:4000])


if __name__ == "__main__":
    main()
