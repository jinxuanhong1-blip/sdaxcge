#!/usr/bin/env python3
"""Compute honest TACSTD2 / CLDN4 presence and Spearman n/ρ/p across downloaded matrices.

No fabricated statistics. If a gene is absent from a panel, ρ is not computed.
Writes per-dataset CSVs and a master summary under results/noskip/zenodo/.
"""
import csv
import gzip
import json
import tarfile
import zipfile
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import stats
from scipy.sparse import csc_matrix

ROOT = Path(__file__).resolve().parents[2]
DL = ROOT / "results" / "noskip" / "zenodo" / "downloads"
OUT = ROOT / "results" / "noskip" / "zenodo"
OUT.mkdir(parents=True, exist_ok=True)

HUMAN = {"TACSTD2", "CLDN4"}
MOUSE = {"Tacstd2", "Cldn4"}


def fmt_p(p):
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return None
    p = float(p)
    if p == 0.0:
        return "<1e-300"  # scipy underflows
    return p


def spearman(x, y, subset):
    n = int(len(x))
    if n < 3:
        return {"subset": subset, "n": n, "rho": None, "p": None, "note": "n<3"}
    rho, p = stats.spearmanr(x, y, nan_policy="omit")
    return {
        "subset": subset,
        "n": n,
        "rho": None if np.isnan(rho) else float(rho),
        "p": fmt_p(p),
        "p_raw_float": float(p) if not np.isnan(p) else None,
        "method": "scipy.stats.spearmanr two-sided",
    }


def extract_10x_h5(path, genes=("TACSTD2", "CLDN4")):
    """Return dict gene->1d counts and n_cells from a 10x filtered_feature_bc_matrix.h5."""
    with h5py.File(path, "r") as f:
        names = [x.decode() if isinstance(x, bytes) else str(x) for x in f["matrix/features/name"][:]]
        data = f["matrix/data"][:]
        indices = f["matrix/indices"][:]
        indptr = f["matrix/indptr"][:]
        shape = tuple(int(x) for x in f["matrix/shape"][:])
    # 10x: CSC, shape = (n_genes, n_cells)
    mat = csc_matrix((data, indices, indptr), shape=shape)
    out = {}
    for g in genes:
        if g in names:
            i = names.index(g)
            out[g] = np.asarray(mat[i, :].todense()).ravel().astype(np.float64)
    return out, names, shape


def analyze_visium_8417887():
    ddir = DL / "8417887"
    rows = []
    xs, ys, sample_ids = [], [], []
    for h5 in sorted(ddir.glob("P*_filtered_feature_bc_matrix.h5")):
        genes, names, shape = extract_10x_h5(h5)
        present = {g: g in genes for g in HUMAN}
        n_spots = int(shape[1])
        rec = {
            "sample": h5.stem.replace("_filtered_feature_bc_matrix", ""),
            "n_spots": n_spots,
            "n_features": int(shape[0]),
            "TACSTD2_present": present["TACSTD2"],
            "CLDN4_present": present["CLDN4"],
        }
        if present["TACSTD2"] and present["CLDN4"]:
            x, y = genes["TACSTD2"], genes["CLDN4"]
            rec["TACSTD2_pct"] = float((x > 0).mean() * 100)
            rec["CLDN4_pct"] = float((y > 0).mean() * 100)
            rec["TACSTD2_mean"] = float(x.mean())
            rec["CLDN4_mean"] = float(y.mean())
            sp = spearman(x, y, rec["sample"])
            rec.update({"rho": sp["rho"], "p": sp["p"], "n": sp["n"]})
            xs.append(x)
            ys.append(y)
            sample_ids.append(np.full(n_spots, rec["sample"]))
        rows.append(rec)
        print("visium", rec)
    tab = pd.DataFrame(rows)
    tab.to_csv(OUT / "8417887_visium_per_sample.csv", index=False)
    pooled = None
    if xs:
        X = np.concatenate(xs)
        Y = np.concatenate(ys)
        pooled = spearman(X, Y, "all_visium_spots_raw_counts")
        pooled["TACSTD2_pct"] = float((X > 0).mean() * 100)
        pooled["CLDN4_pct"] = float((Y > 0).mean() * 100)
        (OUT / "8417887_visium_pooled_spearman.json").write_text(json.dumps(pooled, indent=2))
        print("visium pooled", pooled)
    return {"per_sample": rows, "pooled": pooled}


def stream_gene_rows_txt(path, wanted):
    found = {}
    with open(path, "rt") as fh:
        header = fh.readline()
        n_cells = len(header.rstrip("\n").split("\t"))
        for line in fh:
            gene, _, rest = line.partition("\t")
            if gene in wanted:
                found[gene] = np.fromstring(rest, sep="\t", dtype=np.float64)
                if len(found) == len(wanted):
                    break
    return found, n_cells


def analyze_txt_8417887():
    out = {}
    for label, fname, meta_name, meta_col in [
        ("malignant", "Malignant_cell_data.txt", "Malignant_cell_meta_data.txt", "tissue_type"),
        ("immune", "Immune_cell_data.txt", "Immune_cell_meta_data.txt", "celltype"),
    ]:
        genes, n_cells = stream_gene_rows_txt(DL / "8417887" / fname, HUMAN)
        rec = {
            "table": label,
            "n_cells_header": n_cells,
            "TACSTD2_present": "TACSTD2" in genes,
            "CLDN4_present": "CLDN4" in genes,
            "genes_found": sorted(genes),
        }
        if HUMAN <= set(genes):
            x, y = genes["TACSTD2"], genes["CLDN4"]
            rec["n"] = int(len(x))
            rec["TACSTD2_pct"] = float((x > 0).mean() * 100)
            rec["CLDN4_pct"] = float((y > 0).mean() * 100)
            rec["TACSTD2_mean"] = float(x.mean())
            rec["CLDN4_mean"] = float(y.mean())
            rec["spearman"] = spearman(x, y, f"{label}_raw")
            # optional tissue split
            meta = pd.read_csv(DL / "8417887" / meta_name, sep="\t")
            rec["meta_n"] = int(len(meta))
            rec["meta_cols"] = list(meta.columns)
        print(label, rec)
        out[label] = rec
    # JSON-safe
    def conv(o):
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, (np.integer,)):
            return int(o)
        return o
    (OUT / "8417887_scrna_tables.json").write_text(json.dumps(out, indent=2, default=conv))
    return out


def analyze_11205626():
    tar_path = DL / "11205626" / "24samples.h5.tar.gz"
    meta = pd.read_csv(DL / "11205626" / "24samples_metadata.csv")
    # Use Cell (`_all_`) only to avoid double-counting SN / depletion of the same sample
    wanted = meta[meta["method"] == "Cell"].copy()
    tmp = OUT / "_tmp_h5"
    tmp.mkdir(exist_ok=True)
    rows = []
    xs_n, ys_n, xs_t, ys_t = [], [], [], []
    with tarfile.open(tar_path, "r:gz") as tar:
        for _, row in wanted.iterrows():
            member = f"{row['filenames']}/filtered_feature_bc_matrix.h5"
            try:
                src = tar.extractfile(member)
            except KeyError:
                print("missing", member)
                continue
            dest = tmp / (row["filenames"] + ".h5")
            dest.write_bytes(src.read())
            genes, names, shape = extract_10x_h5(dest)
            dest.unlink()
            rec = {
                "sample": row["filenames"],
                "tissue": row["tissue"],
                "patient": str(row["patient_original_label"]),
                "method": row["method"],
                "n_cells": int(shape[1]),
                "TACSTD2_present": "TACSTD2" in genes,
                "CLDN4_present": "CLDN4" in genes,
            }
            if HUMAN <= set(genes):
                x, y = genes["TACSTD2"], genes["CLDN4"]
                rec["TACSTD2_pct"] = float((x > 0).mean() * 100)
                rec["CLDN4_pct"] = float((y > 0).mean() * 100)
                rec["TACSTD2_mean"] = float(x.mean())
                rec["CLDN4_mean"] = float(y.mean())
                sp = spearman(x, y, rec["sample"])
                rec["rho"] = sp["rho"]
                rec["p"] = sp["p"]
                rec["n"] = sp["n"]
                if row["tissue"] == "Normal":
                    xs_n.append(x)
                    ys_n.append(y)
                else:
                    xs_t.append(x)
                    ys_t.append(y)
            rows.append(rec)
            print("11205626", rec)
    pd.DataFrame(rows).to_csv(OUT / "11205626_cell_protocol_per_sample.csv", index=False)
    pooled = {}
    if xs_n:
        X, Y = np.concatenate(xs_n), np.concatenate(ys_n)
        pooled["normal"] = spearman(X, Y, "normal_cell_protocol")
        pooled["normal"]["TACSTD2_pct"] = float((X > 0).mean() * 100)
        pooled["normal"]["CLDN4_pct"] = float((Y > 0).mean() * 100)
    if xs_t:
        X, Y = np.concatenate(xs_t), np.concatenate(ys_t)
        pooled["tumor"] = spearman(X, Y, "tumor_cell_protocol")
        pooled["tumor"]["TACSTD2_pct"] = float((X > 0).mean() * 100)
        pooled["tumor"]["CLDN4_pct"] = float((Y > 0).mean() * 100)
    if xs_n and xs_t:
        X = np.concatenate(xs_n + xs_t)
        Y = np.concatenate(ys_n + ys_t)
        pooled["all"] = spearman(X, Y, "all_cell_protocol")
        pooled["all"]["TACSTD2_pct"] = float((X > 0).mean() * 100)
        pooled["all"]["CLDN4_pct"] = float((Y > 0).mean() * 100)
    (OUT / "11205626_pooled_spearman.json").write_text(json.dumps(pooled, indent=2))
    print("11205626 pooled", pooled)
    return {"per_sample": rows, "pooled": pooled}


def analyze_mouse():
    path = DL / "10731914" / "mouse_scRNAseq_adata_raw_counts.csv.gz"
    # cells x genes; first column is barcode
    usecols = None
    with gzip.open(path, "rt") as fh:
        header = next(csv.reader(fh))
    genes = [c for c in header if c in MOUSE or c in {"Epcam"}]
    print("mouse columns found:", genes)
    if not ({"Tacstd2", "Cldn4"} <= set(genes)):
        rec = {"TACSTD2_present": "Tacstd2" in genes, "CLDN4_present": "Cldn4" in genes, "columns": genes}
        (OUT / "10731914_mouse_stats.json").write_text(json.dumps(rec, indent=2))
        return rec
    df = pd.read_csv(path, usecols=["Tacstd2", "Cldn4"])
    meta = pd.read_csv(DL / "10731914" / "mouse_scRNAseq_adata_metadata.csv")
    x = df["Tacstd2"].to_numpy(dtype=float)
    y = df["Cldn4"].to_numpy(dtype=float)
    rec = {
        "n_cells": int(len(x)),
        "TACSTD2_present": True,
        "CLDN4_present": True,
        "symbol_note": "mouse symbols Tacstd2 / Cldn4",
        "Tacstd2_pct": float((x > 0).mean() * 100),
        "Cldn4_pct": float((y > 0).mean() * 100),
        "spearman_all": spearman(x, y, "mouse_all_raw"),
        "meta_n": int(len(meta)),
        "meta_cols": list(meta.columns),
    }
    if "treatmentGroup" in meta.columns and len(meta) == len(x):
        by = []
        for g, idx in meta.groupby("treatmentGroup").groups.items():
            xi, yi = x[list(idx)], y[list(idx)]
            by.append({
                "treatmentGroup": str(g),
                "n": int(len(idx)),
                "Tacstd2_pct": float((xi > 0).mean() * 100),
                "Cldn4_pct": float((yi > 0).mean() * 100),
                **spearman(xi, yi, str(g)),
            })
        rec["by_treatment"] = by
        pd.DataFrame(by).to_csv(OUT / "10731914_mouse_by_treatment.csv", index=False)
    (OUT / "10731914_mouse_stats.json").write_text(json.dumps(rec, indent=2))
    print("mouse", rec)
    return rec


def analyze_nanostring():
    learn = pd.read_csv(DL / "2635194" / "learn-data.csv")
    val = pd.read_csv(DL / "2635194" / "validate-data.csv")
    genes = [c for c in learn.columns if c != "Unnamed: 0"]
    rec = {
        "doi": "10.5281/zenodo.2635194",
        "assay": "Nanostring nCounter, advanced NSCLC anti-PD-1",
        "n_learn": int(learn.shape[0]),
        "n_validate": int(val.shape[0]),
        "n_genes_in_panel": int(len(genes)),
        "TACSTD2_present": "TACSTD2" in genes,
        "CLDN4_present": "CLDN4" in genes,
        "TROP2_present": "TROP2" in genes,
        "rho": None,
        "p": None,
        "note": "TACSTD2 and CLDN4 are not on this immune-focused nCounter panel. Spearman was not computed.",
    }
    (OUT / "2635194_nanostring_absent.json").write_text(json.dumps(rec, indent=2))
    print("nanostring", rec)
    return rec


def analyze_i3lung():
    rec = {
        "doi": "10.64898/2026.01.16.25342913",
        "files": ["cb.csv", "outcomes.csv", "genomics.csv", "fmrad.csv", "digital_pathology.csv"],
        "TACSTD2_present": False,
        "CLDN4_present": False,
        "rho": None,
        "p": None,
        "note": "I3LUNG_DATA.zip is clinical / radiomics / driver-mutation tables (KRAS/P53/STK11). No expression matrix. results.zip (1.9 GB) was not needed after this negative.",
    }
    (OUT / "17535424_i3lung_absent.json").write_text(json.dumps(rec, indent=2))
    print("i3lung", rec)
    return rec


def analyze_fig6_tcells():
    """T-cell-only h5ad objects: check whether TACSTD2/CLDN4 exist in var."""
    zpath = DL / "13947395" / "Fig6_processed_scRNAseq.zip"
    recs = []
    with zipfile.ZipFile(zpath) as z:
        for name in z.namelist():
            if not name.endswith(".h5ad") or "site-packages" in name:
                continue
            tmp = OUT / "_tmp_h5ad"
            tmp.mkdir(exist_ok=True)
            dest = tmp / Path(name).name
            with z.open(name) as src, open(dest, "wb") as fh:
                while True:
                    chunk = src.read(1 << 20)
                    if not chunk:
                        break
                    fh.write(chunk)
            genes = []
            n_obs = None
            with h5py.File(dest, "r") as f:
                # anndata var names usually /var/_index or /var/gene_ids
                for key in ["var/_index", "var/gene_ids", "var/features", "var/gene_symbols"]:
                    if key in f:
                        genes = [x.decode() if isinstance(x, bytes) else str(x) for x in f[key][:]]
                        break
                if "obs/_index" in f:
                    n_obs = int(f["obs/_index"].shape[0])
            dest.unlink()
            rec = {
                "file": Path(name).name,
                "n_obs": n_obs,
                "n_var_checked": len(genes),
                "TACSTD2_present": "TACSTD2" in genes,
                "CLDN4_present": "CLDN4" in genes,
            }
            recs.append(rec)
            print("fig6", rec)
    (OUT / "13947395_tcell_h5ad_gene_check.json").write_text(json.dumps(recs, indent=2))
    return recs


def main():
    summary = {
        "10731914_human": json.loads((OUT / "10731914_stats.json").read_text()) if (OUT / "10731914_stats.json").exists() else None,
        "2635194_nanostring": analyze_nanostring(),
        "17535424_i3lung": analyze_i3lung(),
        "8417887_visium": analyze_visium_8417887(),
        "8417887_scrna": analyze_txt_8417887(),
        "11205626_paired_luad": analyze_11205626(),
        "10731914_mouse": analyze_mouse(),
        "13947395_tcells": analyze_fig6_tcells(),
        "8041882_imc": {
            "doi": "10.5281/zenodo.8041882",
            "TACSTD2_present": False,
            "CLDN4_present": False,
            "n_cells": 99659,
            "n_markers": 41,
            "rho": None,
            "p": None,
            "note": "IMC protein panel; TACSTD2/CLDN4 not among 41 antibodies. No ρ computed.",
        },
    }
    (OUT / "master_summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print("\nWrote", OUT / "master_summary.json")


if __name__ == "__main__":
    main()
