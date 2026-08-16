#!/usr/bin/env python3
"""Score epithelial TACSTD2/CLDN4 vs T/NK fraction in TISCH2 NSCLC objects.

Unit of analysis is Sample when CellMetainfo has a real sample column,
otherwise Patient. Tumor-like tissue only. Honest skips are written, not
silenced. Cell-level p-values are not used as primary evidence.
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

from config import (
    DATASETS,
    EPITHELIAL_NONMALIGNANT,
    EXCLUDE_TISSUE_TOKENS,
    GENES,
    MALIGNANT,
    MIN_EPI_CELLS,
    MIN_N_FOR_SPEARMAN,
    MIN_TNK_CELLS,
    TARGET_GENES,
    TNK_LINEAGES,
    TUMOR_TOKENS,
)


def _as_str_array(arr) -> np.ndarray:
    out = []
    for x in arr:
        if isinstance(x, bytes):
            out.append(x.decode())
        else:
            out.append(str(x))
    return np.array(out, dtype=object)


def load_tisch_genes(path: Path, wanted: list[str]) -> tuple[pd.DataFrame, list[str], dict]:
    """Return barcode-indexed gene table + present genes + h5 audit."""
    audit = {"path": str(path), "keys": [], "shape": None, "n_genes": None, "n_cells": None}
    with h5py.File(path, "r") as f:
        audit["keys"] = list(f.keys())
        grp = None
        if all(k in f for k in ("data", "indices", "indptr")):
            grp = f
        else:
            for key in f.keys():
                g = f[key]
                if isinstance(g, h5py.Group) and all(k in g for k in ("data", "indices", "indptr")):
                    grp = g
                    break
        if grp is None:
            raise ValueError(f"no sparse matrix in {path}: {audit['keys']}")

        def _names(candidates):
            for c in candidates:
                node = grp.get(c) if hasattr(grp, "get") else None
                if node is None:
                    node = f.get(c)
                if node is None:
                    continue
                if isinstance(node, h5py.Dataset):
                    return _as_str_array(node[:])
                if isinstance(node, h5py.Group):
                    for sub in ("name", "id", "gene_names"):
                        if sub in node:
                            return _as_str_array(node[sub][:])
            return None

        genes = _names(["features", "gene_names", "genes", "rownames"])
        barcodes = _names(["barcodes", "cell_names", "colnames"])
        shape = tuple(int(x) for x in grp["shape"][:]) if "shape" in grp else None
        data = grp["data"][:]
        indices = grp["indices"][:]
        indptr = grp["indptr"][:]

    if shape is None:
        raise ValueError(f"{path}: missing shape")
    n0, n1 = int(shape[0]), int(shape[1])
    # 10x/MAESTRO: CSC genes x cells. Same buffers can be read as CSR cells x genes.
    if genes is not None and barcodes is not None and len(genes) == n1 and len(barcodes) == n0:
        n_genes, n_cells = n1, n0
        csc = sp.csc_matrix((data, indices, indptr), shape=(n_cells, n_genes))
        # columns = genes
        def col(i: int) -> np.ndarray:
            return np.asarray(csc[:, i].todense()).ravel()
    else:
        n_genes, n_cells = n0, n1
        csc = sp.csc_matrix((data, indices, indptr), shape=(n_genes, n_cells))

        def col(i: int) -> np.ndarray:
            # row i of CSC (n_genes, n_cells) == gene i across cells
            return np.asarray(csc.getrow(i).todense()).ravel()

    if genes is None or barcodes is None:
        raise ValueError(f"{path}: missing gene or barcode names")
    if len(genes) != n_genes or len(barcodes) != n_cells:
        raise ValueError(
            f"{path}: dim mismatch genes={len(genes)} barcodes={len(barcodes)} "
            f"shape=({n_genes},{n_cells})"
        )

    audit["shape"] = [n_genes, n_cells]
    audit["n_genes"] = n_genes
    audit["n_cells"] = n_cells
    present = [g for g in wanted if g in set(genes)]
    audit["genes_present"] = present
    audit["genes_absent"] = [g for g in wanted if g not in set(genes)]
    gene_index = {g: int(np.where(genes == g)[0][0]) for g in present}
    table = {"barcode": barcodes}
    for g, i in gene_index.items():
        table[g] = col(i)
    return pd.DataFrame(table).set_index("barcode"), present, audit


def _norm_label(x) -> str:
    return str(x).strip()


def is_tumor_tissue(val: str) -> bool:
    s = str(val).strip().lower()
    if not s or s in {"nan", "none", "na"}:
        return False
    if any(tok in s for tok in EXCLUDE_TISSUE_TOKENS):
        return False
    return any(tok in s for tok in TUMOR_TOKENS)


def pick_unit_columns(meta: pd.DataFrame) -> tuple[str | None, str | None, str]:
    """Return (patient_col, sample_col, unit_kind)."""
    cols = {c.lower(): c for c in meta.columns}
    patient = cols.get("patient")
    sample = cols.get("sample")
    # GSE151537 stores a per-cell id in Sample (n_sample == n_cells).
    if sample is not None and meta[sample].nunique(dropna=True) >= 0.9 * len(meta):
        sample = None
    if sample is not None:
        return patient, sample, "sample"
    if patient is not None:
        return patient, None, "patient"
    return None, None, "none"


def pick_tissue_col(meta: pd.DataFrame) -> str | None:
    for name in ("Source", "Tissue", "source", "tissue"):
        if name in meta.columns:
            return name
    return None


def ici_columns(meta: pd.DataFrame) -> list[str]:
    keys = ("response", "treatment", "therapy", "ici", "icb", "recist", "mpr", "pdl1", "pd-1", "pd1")
    out = []
    for c in meta.columns:
        cl = c.lower()
        if any(k in cl for k in keys):
            out.append(c)
    return out


def spearman_row(x: np.ndarray, y: np.ndarray, contrast: str, dataset: str, unit: str) -> dict:
    n = int(len(x))
    rec = {
        "dataset": dataset,
        "unit": unit,
        "contrast": contrast,
        "n": n,
        "rho": np.nan,
        "p": np.nan,
        "ok": False,
        "note": "",
    }
    if n < MIN_N_FOR_SPEARMAN:
        rec["note"] = f"n={n}<{MIN_N_FOR_SPEARMAN}; Spearman not computed"
        return rec
    if np.nanstd(x) == 0 or np.nanstd(y) == 0:
        rec["note"] = "zero variance; Spearman undefined"
        return rec
    rho, p = stats.spearmanr(x, y, nan_policy="omit")
    rec["rho"] = float(rho)
    rec["p"] = float(p)
    rec["ok"] = True
    return rec


def analyze_dataset(ds: str, data_dir: Path) -> dict:
    spec = DATASETS[ds]
    meta_path = data_dir / f"{ds}_CellMetainfo_table.tsv"
    h5_path = data_dir / f"{ds}_expression.h5"
    rec: dict = {
        "dataset": ds,
        "pmid": spec.get("pmid"),
        "paper": spec.get("paper"),
        "ici_gallery": spec.get("ici_gallery"),
        "priority": spec.get("priority"),
        "included_in_pool": False,
        "skip_reason": None,
    }
    if not meta_path.exists():
        rec["skip_reason"] = "CellMetainfo missing"
        return rec, pd.DataFrame(), pd.DataFrame()

    meta = pd.read_csv(meta_path, sep="\t", low_memory=False)
    cell_col = meta.columns[0]
    meta = meta.set_index(cell_col)
    meta.index = meta.index.astype(str)
    lineage_col = next((c for c in meta.columns if "major-lineage" in c.lower()), None)
    if lineage_col is None:
        rec["skip_reason"] = "no major-lineage column"
        return rec, pd.DataFrame(), pd.DataFrame()
    meta["_lineage"] = meta[lineage_col].map(_norm_label)
    rec["n_cells_meta"] = int(len(meta))
    rec["lineages"] = meta["_lineage"].value_counts().to_dict()
    rec["n_malignant"] = int((meta["_lineage"].isin(MALIGNANT)).sum())
    rec["n_epithelial_nonmalig"] = int((meta["_lineage"].isin(EPITHELIAL_NONMALIGNANT)).sum())
    rec["n_tnk"] = int((meta["_lineage"].isin(TNK_LINEAGES)).sum())
    rec["ici_metainfo_columns"] = ici_columns(meta)
    rec["ici_metainfo_values"] = {
        c: meta[c].astype(str).value_counts().head(12).to_dict() for c in rec["ici_metainfo_columns"]
    }
    patient_col, sample_col, unit_kind = pick_unit_columns(meta)
    rec["patient_col"] = patient_col
    rec["sample_col"] = sample_col
    rec["unit"] = unit_kind
    tissue_col = pick_tissue_col(meta)
    rec["tissue_col"] = tissue_col
    if tissue_col:
        rec["tissue_values"] = meta[tissue_col].astype(str).value_counts().to_dict()

    if rec["n_malignant"] + rec["n_epithelial_nonmalig"] == 0:
        rec["skip_reason"] = spec.get("skip_reason_expected") or "no epithelial/malignant cells"
        return rec, pd.DataFrame(), pd.DataFrame()
    if rec["n_tnk"] == 0:
        rec["skip_reason"] = "no T/NK cells"
        return rec, pd.DataFrame(), pd.DataFrame()
    if unit_kind == "none":
        rec["skip_reason"] = spec.get("skip_reason_expected") or "no Patient/Sample column"
        return rec, pd.DataFrame(), pd.DataFrame()
    if not spec.get("download_h5", False) or not h5_path.exists():
        rec["skip_reason"] = "expression.h5 not downloaded"
        return rec, pd.DataFrame(), pd.DataFrame()

    expr, present, h5_audit = load_tisch_genes(h5_path, GENES)
    rec["h5"] = h5_audit
    rec["genes_present"] = present
    if "TACSTD2" not in present:
        rec["skip_reason"] = "TACSTD2 absent from TISCH h5"
        return rec, pd.DataFrame(), pd.DataFrame()

    # Align
    common = meta.index.intersection(expr.index)
    rec["n_cells_aligned"] = int(len(common))
    if len(common) < 50:
        rec["skip_reason"] = f"barcode alignment too small (n={len(common)})"
        return rec, pd.DataFrame(), pd.DataFrame()
    meta = meta.loc[common]
    expr = expr.loc[common]

    if tissue_col is not None:
        tumor_mask = meta[tissue_col].map(is_tumor_tissue)
        rec["n_tumor_like_cells"] = int(tumor_mask.sum())
        if tumor_mask.sum() >= 50:
            meta = meta.loc[tumor_mask]
            expr = expr.loc[meta.index]
        else:
            rec["tissue_filter"] = "kept all cells; tumor-like filter left <50 cells"
    else:
        rec["tissue_filter"] = "no Source/Tissue column; all cells kept"

    unit_col = sample_col if unit_kind == "sample" else patient_col
    meta["_unit"] = meta[unit_col].astype(str)
    meta["_is_malig"] = meta["_lineage"].isin(MALIGNANT)
    meta["_is_epi"] = meta["_lineage"].isin(MALIGNANT | EPITHELIAL_NONMALIGNANT)
    meta["_is_tnk"] = meta["_lineage"].isin(TNK_LINEAGES)

    rows = []
    for unit, sub in meta.groupby("_unit", sort=True):
        if unit in {"nan", "None", ""}:
            continue
        n_cells = int(len(sub))
        n_mal = int(sub["_is_malig"].sum())
        n_epi = int(sub["_is_epi"].sum())
        n_tnk = int(sub["_is_tnk"].sum())
        # Prefer malignant when enough cells; else all epithelial-like.
        if n_mal >= MIN_EPI_CELLS:
            score_mask = sub["_is_malig"].to_numpy()
            epi_def = "Malignant"
        else:
            score_mask = sub["_is_epi"].to_numpy()
            epi_def = "epithelial_like"
        n_score = int(score_mask.sum())
        eligible = n_score >= MIN_EPI_CELLS and n_tnk >= MIN_TNK_CELLS
        row = {
            "dataset": ds,
            "unit_kind": unit_kind,
            "unit_id": unit,
            "patient": str(sub[patient_col].iloc[0]) if patient_col else "",
            "n_cells": n_cells,
            "n_malignant": n_mal,
            "n_epithelial_like": n_epi,
            "n_tnk": n_tnk,
            "frac_tnk": n_tnk / n_cells if n_cells else np.nan,
            "epi_definition": epi_def,
            "n_scored_epi": n_score,
            "eligible": eligible,
            "ici_gallery": spec.get("ici_gallery"),
        }
        idx = sub.index
        for g in present:
            v = expr.loc[idx, g].to_numpy()
            scored = v[score_mask]
            row[f"{g}_epi_mean"] = float(np.mean(scored)) if n_score else np.nan
            row[f"{g}_epi_pctpos"] = float(np.mean(scored > 0) * 100) if n_score else np.nan
        rows.append(row)

    units = pd.DataFrame(rows)
    rec["n_units_total"] = int(len(units))
    rec["n_units_eligible"] = int(units["eligible"].sum()) if len(units) else 0
    elig = units.loc[units["eligible"]].copy() if len(units) else units

    stats_rows = []
    for g in [x for x in TARGET_GENES if x in present]:
        for metric, col in (("mean", f"{g}_epi_mean"), ("pctpos", f"{g}_epi_pctpos")):
            stats_rows.append(
                spearman_row(
                    elig[col].to_numpy() if len(elig) else np.array([]),
                    elig["frac_tnk"].to_numpy() if len(elig) else np.array([]),
                    f"{g}_{metric}_vs_fracTNK",
                    ds,
                    unit_kind,
                )
            )
    # Controls
    for g in ("EPCAM", "PTPRC"):
        if g in present and len(elig):
            stats_rows.append(
                spearman_row(
                    elig[f"{g}_epi_mean"].to_numpy(),
                    elig["frac_tnk"].to_numpy(),
                    f"{g}_mean_vs_fracTNK",
                    ds,
                    unit_kind,
                )
            )

    stats_df = pd.DataFrame(stats_rows)
    rec["statistics"] = stats_df.to_dict(orient="records")
    if rec["n_units_eligible"] >= MIN_N_FOR_SPEARMAN and "TACSTD2" in present:
        rec["included_in_pool"] = True
    else:
        if rec["skip_reason"] is None:
            rec["skip_reason"] = (
                f"eligible {unit_kind}s n={rec['n_units_eligible']} "
                f"(need >={MIN_N_FOR_SPEARMAN} for Spearman / pool)"
            )
    rec["has_ici_labels"] = bool(rec["ici_metainfo_columns"]) or spec.get("ici_gallery") == "Immunotherapy"
    rec["ici_label_honest"] = (
        "TISCH gallery=Immunotherapy; CellMetainfo has no response/RECIST/MPR column"
        if spec.get("ici_gallery") == "Immunotherapy" and not rec["ici_metainfo_columns"]
        else (
            "CellMetainfo treatment/response columns present"
            if rec["ici_metainfo_columns"]
            else "no ICI labels in TISCH gallery or CellMetainfo"
        )
    )
    return rec, units, stats_df


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="data/tisch_nsclc_pool")
    ap.add_argument("--out-dir", default="methods/tisch_nsclc_pool")
    args = ap.parse_args()
    data_dir = Path(args.data_dir)
    out_dir = Path(args.out_dir)
    tables = out_dir / "tables"
    tables.mkdir(parents=True, exist_ok=True)

    audits = []
    all_units = []
    all_stats = []
    for ds in DATASETS:
        print(f"== {ds}", flush=True)
        rec, units, stats_df = analyze_dataset(ds, data_dir)
        audits.append(rec)
        if len(units):
            all_units.append(units)
        if len(stats_df):
            all_stats.append(stats_df)
        print(
            f"   skip={rec.get('skip_reason')} eligible={rec.get('n_units_eligible')} "
            f"pool={rec.get('included_in_pool')}",
            flush=True,
        )

    units_df = pd.concat(all_units, ignore_index=True) if all_units else pd.DataFrame()
    stats_df = pd.concat(all_stats, ignore_index=True) if all_stats else pd.DataFrame()
    units_path = tables / "per_unit_metrics.tsv"
    stats_path = tables / "per_dataset_spearman.tsv"
    audit_path = out_dir / "dataset_audit.json"
    units_df.to_csv(units_path, sep="\t", index=False)
    stats_df.to_csv(stats_path, sep="\t", index=False)
    payload = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "n_datasets_catalog": len(DATASETS),
        "n_included_in_pool": sum(1 for r in audits if r.get("included_in_pool")),
        "datasets": audits,
    }
    audit_path.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"wrote {units_path} {stats_path} {audit_path}")


if __name__ == "__main__":
    main()
