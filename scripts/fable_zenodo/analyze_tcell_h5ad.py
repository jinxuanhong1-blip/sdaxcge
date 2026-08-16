#!/usr/bin/env python3
"""Extract TACSTD2/CLDN4 from Fig6 T-cell h5ad objects and compute Spearman if both vary."""
import json
import zipfile
from pathlib import Path

import h5py
import numpy as np
from scipy import stats
from scipy.sparse import csr_matrix, csc_matrix

ROOT = Path(__file__).resolve().parents[2]
ZPATH = ROOT / "results" / "noskip" / "zenodo" / "downloads" / "13947395" / "Fig6_processed_scRNAseq.zip"
OUT = ROOT / "results" / "noskip" / "zenodo"


def fmt_p(p):
    p = float(p)
    return "<1e-300" if p == 0.0 else p


def decode(arr):
    return [x.decode() if isinstance(x, bytes) else str(x) for x in arr]


def read_var_names(f):
    for key in ["var/_index", "var/gene_ids", "var/gene_symbols", "var/features", "raw/var/_index"]:
        if key in f:
            return decode(f[key][:]), key
    # fallback: any 1d string dataset under var
    if "var" in f:
        for k, v in f["var"].items():
            if isinstance(v, h5py.Dataset) and v.dtype.kind in ("S", "O", "U") and v.ndim == 1:
                return decode(v[:]), f"var/{k}"
    return [], None


def read_X_gene_columns(f, gene_idx):
    """Return dense arrays (n_obs,) for each gene. AnnData X is cells x genes."""
    X = f["X"]
    if isinstance(X, h5py.Dataset):
        return np.asarray(X[:, gene_idx], dtype=np.float64).T
    if "data" in X and "indices" in X and "indptr" in X:
        data = X["data"][:]
        indices = X["indices"][:]
        indptr = X["indptr"][:]
        attrs = dict(X.attrs)
        encoding = str(attrs.get("encoding-type") or attrs.get("h5sparse_format") or "")
        shape_attr = attrs.get("shape", attrs.get("h5sparse_shape", None))
        shape = tuple(int(x) for x in np.asarray(shape_attr).ravel()) if shape_attr is not None else None
        if shape is None:
            if "csr" in encoding:
                shape = (len(indptr) - 1, int(indices.max()) + 1 if len(indices) else 0)
            else:
                shape = (int(indices.max()) + 1 if len(indices) else 0, len(indptr) - 1)
        if "csc" in encoding:
            mat = csc_matrix((data, indices, indptr), shape=shape)
        else:
            mat = csr_matrix((data, indices, indptr), shape=shape)
        return np.asarray(mat[:, gene_idx].todense(), dtype=np.float64).T
    raise RuntimeError("unrecognized X layout")


def analyze_one(path):
    with h5py.File(path, "r") as f:
        print(path.name, "top", list(f.keys()))
        if "var" in f:
            print("  var keys", list(f["var"].keys())[:20])
        genes, src = read_var_names(f)
        n_obs = None
        if "obs" in f and "_index" in f["obs"]:
            n_obs = int(f["obs/_index"].shape[0])
        rec = {
            "file": path.name,
            "var_source": src,
            "n_var": len(genes),
            "n_obs": n_obs,
            "TACSTD2_present": "TACSTD2" in genes,
            "CLDN4_present": "CLDN4" in genes,
        }
        if rec["TACSTD2_present"] and rec["CLDN4_present"] and "X" in f:
            i = [genes.index("TACSTD2"), genes.index("CLDN4")]
            rows = read_X_gene_columns(f, i)
            x, y = rows[0], rows[1]
            rec["n"] = int(len(x))
            rec["TACSTD2_pct"] = float((x > 0).mean() * 100)
            rec["CLDN4_pct"] = float((y > 0).mean() * 100)
            rec["TACSTD2_mean"] = float(x.mean())
            rec["CLDN4_mean"] = float(y.mean())
            if np.unique(x).size > 1 and np.unique(y).size > 1:
                rho, p = stats.spearmanr(x, y, nan_policy="omit")
                rec["rho"] = float(rho)
                rec["p"] = fmt_p(p)
                rec["p_raw_float"] = float(p)
            else:
                rec["rho"] = None
                rec["p"] = None
                rec["note"] = "at least one gene is constant; Spearman undefined"
        return rec


def main():
    recs = []
    tmp = OUT / "_tmp_h5ad"
    tmp.mkdir(exist_ok=True)
    with zipfile.ZipFile(ZPATH) as z:
        for name in z.namelist():
            if not name.endswith(".h5ad") or "site-packages" in name:
                continue
            dest = tmp / Path(name).name
            with z.open(name) as src, open(dest, "wb") as fh:
                while True:
                    c = src.read(1 << 20)
                    if not c:
                        break
                    fh.write(c)
            rec = analyze_one(dest)
            print(rec)
            recs.append(rec)
            dest.unlink()
    (OUT / "13947395_tcell_expression.json").write_text(json.dumps(recs, indent=2))


if __name__ == "__main__":
    main()
