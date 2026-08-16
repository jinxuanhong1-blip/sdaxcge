#!/usr/bin/env python3
"""Build per-patient lineage counts + malignant TACSTD2 from public GEO files.

Cohorts
-------
GSE207422  Hu et al. 2023 neoadjuvant PD-1 + chemo. Public dense UMI TSV.
           Marker hierarchy (author CopyKAT barcodes are not on GEO).
GSE241934  Zhang et al. 2024 NEOTIDE/CTONG2104 + real-world. Author
           major.cell.type in the GEO metadata; TACSTD2 streamed from MTX.
GSE291670  Xia et al. 2025 anlotinib + camrelizumab. Public 10x MTX.
           Same marker hierarchy as GSE207422.

Malignant TACSTD2 = mean log1p(UMI) in epithelial cells of that patient.
pCR is grouped with MPR. Patient/sample is the unit.
"""
from __future__ import annotations

import gzip
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUT = ROOT / "results"

MARKER_SETS = {
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT7"],
    "T": ["CD3D", "CD3E", "CD3G", "CD2", "CD8A", "CD4", "IL7R"],
    "NK": ["NKG7", "GNLY", "KLRD1", "NCAM1"],
    "B": ["MS4A1", "CD79A"],
    "Myeloid": ["LYZ", "CD68", "CD14", "FCGR3A"],
    "Endo": ["PECAM1", "VWF"],
    "Fibro": ["COL1A1", "DCN", "PDGFRA"],
    "Mast": ["TPSAB1", "KIT"],
}
# Assignment priority (first match wins). T/NK before other immune.
PRIORITY = ["T", "NK", "B", "Myeloid", "Mast", "Epithelial", "Endo", "Fibro"]


def _group_of(x) -> str:
    s = str(x)
    if s in {"non-MPR", "NMPR", "nonMPR"} or s.startswith("NMPR"):
        return "NMPR"
    if s in {"MPR", "pCR", "MPR (pCR)"} or s.startswith("MPR") or s == "pCR" or "pCR" in s:
        return "MPR"
    return ""


def _collapse_lineage(name: str) -> str:
    if name in {"T", "NK", "T_NK"}:
        return "T_NK"
    if name in {"Epi", "Epithelial"}:
        return "Epithelial"
    if name in {"B"}:
        return "B"
    if name in {"Myeloid"}:
        return "Myeloid"
    if name in {"Fibro", "Endo", "Mast", "Stromal"}:
        return "Stromal"
    return "Other"


def _empty_counts() -> dict[str, int]:
    return {k: 0 for k in ["Epithelial", "T_NK", "B", "Myeloid", "Stromal", "Other"]}


def assign_markers(df: pd.DataFrame) -> pd.Series:
    """Exclusive marker hierarchy on raw UMI presence."""
    ptprc = df["PTPRC"] > 0 if "PTPRC" in df.columns else pd.Series(False, index=df.index)
    hits = {}
    for name, genes in MARKER_SETS.items():
        present = [g for g in genes if g in df.columns]
        hits[name] = df[present].gt(0).any(axis=1) if present else pd.Series(False, index=df.index)
    # Immune lineages require PTPRC+; epithelial/stromal require PTPRC-.
    immune = {"T", "NK", "B", "Myeloid", "Mast"}
    assigned = pd.Series("Other", index=df.index)
    remaining = pd.Series(True, index=df.index)
    for name in PRIORITY:
        ok = hits[name] & remaining
        if name in immune:
            ok = ok & ptprc
        else:
            ok = ok & ~ptprc
        assigned[ok] = name
        remaining = remaining & ~ok
    return assigned


def gene_row_index(features_path: Path, genes: list[str]) -> dict[str, int]:
    out: dict[str, int] = {}
    with gzip.open(features_path, "rt") as f:
        for i, line in enumerate(f, start=1):
            parts = line.rstrip("\n").split("\t")
            for g in genes:
                if g in parts[:2] and g not in out:
                    out[g] = i
    return out


def stream_mtx_genes(mtx_path: Path, rows: dict[str, int], n_cells: int) -> dict[str, np.ndarray]:
    want = {v: k for k, v in rows.items()}
    vals = {g: np.zeros(n_cells, dtype=np.float32) for g in rows}
    with gzip.open(mtx_path, "rt") as f:
        for line in f:
            if line.startswith("%"):
                continue
            break
        for line in f:
            r_s, c_s, v_s = line.split()
            r = int(r_s)
            if r in want:
                vals[want[r]][int(c_s) - 1] = float(v_s)
    return vals


# ---------------------------------------------------------------------------
# GSE207422
# ---------------------------------------------------------------------------
def build_gse207422() -> pd.DataFrame:
    marker_path = DATA / "GSE207422" / "GSE207422_markers.tsv.gz"
    raw = pd.read_csv(marker_path, sep="\t", index_col=0)
    raw = raw[~raw.index.duplicated(keep="first")]
    expr = raw.T.astype(np.float32)
    expr.index.name = "cell"
    expr["sample"] = [c.rsplit("_", 1)[0] for c in expr.index]
    expr["lineage_fine"] = assign_markers(expr)
    expr["lineage"] = expr["lineage_fine"].map(_collapse_lineage)
    expr["is_epi"] = expr["lineage"].eq("Epithelial")
    expr["tac_log1p"] = np.log1p(expr["TACSTD2"].astype(float))

    meta = pd.read_excel(DATA / "GSE207422" / "GSE207422_NSCLC_scRNAseq_metadata.xlsx")
    meta = meta.dropna(subset=["Sample"]).copy()
    meta = meta[meta["Sample"].astype(str).str.startswith("BD_")]
    meta["group"] = meta["Pathologic Response"].map(_group_of)
    meta = meta.set_index("Sample")

    rows = []
    for sample, sub in expr.groupby("sample"):
        counts = _empty_counts()
        for lin, n in sub["lineage"].value_counts().items():
            counts[lin] = int(n)
        n_epi = int(sub["is_epi"].sum())
        if n_epi:
            tac = float(sub.loc[sub["is_epi"], "tac_log1p"].mean())
            tac_pos = float((sub.loc[sub["is_epi"], "TACSTD2"] > 0).mean())
        else:
            tac = tac_pos = np.nan
        rec = {
            "dataset": "GSE207422",
            "cohort": "GSE207422",
            "sample_id": sample,
            "patient_id": meta.loc[sample, "Patient"] if sample in meta.index else sample,
            "n_cells": int(len(sub)),
            "n_epi": n_epi,
            "malignant_TACSTD2_mean_log1p": tac,
            "malignant_TACSTD2_frac_pos": tac_pos,
            "annotation": "marker_hierarchy",
            **counts,
        }
        if sample in meta.index:
            rec.update(
                {
                    "group": meta.loc[sample, "group"],
                    "pathologic_response": str(meta.loc[sample, "Pathologic Response"]),
                    "timepoint": str(meta.loc[sample, "Resource"]),
                    "histology": str(meta.loc[sample, "Pathology"]),
                    "pd1": str(meta.loc[sample, "PD1 Antibody"]),
                }
            )
        else:
            rec.update({"group": "", "pathologic_response": "", "timepoint": "", "histology": "", "pd1": ""})
        rows.append(rec)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# GSE241934
# ---------------------------------------------------------------------------
def _gse241934_cohort(name: str, meta_path: Path, mtx_path: Path, feat_path: Path, bc_path: Path) -> pd.DataFrame:
    print(f"GSE241934 {name}: streaming TACSTD2 from {mtx_path.name}", flush=True)
    meta = pd.read_csv(meta_path, sep="\t", dtype=str, low_memory=False)
    major = meta["major.cell.type"] if "major.cell.type" in meta.columns else meta["major_cell_type"]
    meta = meta.copy()
    meta["lineage"] = major.map(_collapse_lineage)
    meta["is_epi"] = major.eq("Epi")
    with gzip.open(bc_path, "rt") as f:
        barcodes = [line.rstrip("\n") for line in f]
    n_cells = len(barcodes)
    rows = gene_row_index(feat_path, ["TACSTD2"])
    if "TACSTD2" not in rows:
        raise RuntimeError(f"TACSTD2 missing in {feat_path}")
    tac = stream_mtx_genes(mtx_path, rows, n_cells)["TACSTD2"]
    # barcodes == cellID and same order as meta (verified)
    if "cellID" in meta.columns:
        if list(meta["cellID"].astype(str)) != barcodes:
            meta = meta.set_index("cellID").reindex(barcodes).reset_index()
    meta["TACSTD2"] = tac
    meta["tac_log1p"] = np.log1p(meta["TACSTD2"].astype(float))
    meta["group"] = meta["Pathological Response"].map(_group_of)
    out = []
    for sid, sub in meta.groupby("sampleID"):
        counts = _empty_counts()
        for lin, n in sub["lineage"].value_counts().items():
            counts[lin] = int(n)
        n_epi = int(sub["is_epi"].sum())
        if n_epi:
            tmean = float(sub.loc[sub["is_epi"], "tac_log1p"].mean())
            tpos = float((sub.loc[sub["is_epi"], "TACSTD2"] > 0).mean())
        else:
            tmean = tpos = np.nan
        out.append(
            {
                "dataset": "GSE241934",
                "cohort": name,
                "sample_id": str(sid),
                "patient_id": str(sid),
                "n_cells": int(len(sub)),
                "n_epi": n_epi,
                "malignant_TACSTD2_mean_log1p": tmean,
                "malignant_TACSTD2_frac_pos": tpos,
                "group": sub["group"].iloc[0],
                "pathologic_response": str(sub["Pathological Response"].iloc[0]),
                "timepoint": "Post-treatment surgery",
                "histology": str(sub["Histology"].iloc[0]) if "Histology" in sub.columns else "",
                "pd1": str(sub["PD1"].iloc[0]) if "PD1" in sub.columns else "",
                "annotation": "author_major.cell.type",
                **counts,
            }
        )
    return pd.DataFrame(out)


def build_gse241934() -> pd.DataFrame:
    d = DATA / "GSE241934"
    iit = _gse241934_cohort(
        "GSE241934_IIT",
        d / "GSE241934_IIT_Meta.txt.gz",
        d / "GSE241934_IIT_Matrix.mtx.gz",
        d / "GSE241934_IIT_features.tsv.gz",
        d / "GSE241934_IIT_barcodes.tsv.gz",
    )
    real = _gse241934_cohort(
        "GSE241934_REAL",
        d / "GSE241934_Real_Meta.txt.gz",
        d / "GSE241934_Real_Matrix.mtx.gz",
        d / "GSE241934_RWC_features.tsv.gz",
        d / "GSE241934_RWC_barcodes.tsv.gz",
    )
    return pd.concat([iit, real], ignore_index=True)


# ---------------------------------------------------------------------------
# GSE291670
# ---------------------------------------------------------------------------
def build_gse291670() -> pd.DataFrame:
    d = DATA / "GSE291670"
    tar = d / "GSE291670_RAW.tar"
    if tar.exists() and not any(d.glob("*_matrix.mtx.gz")):
        with tarfile.open(tar) as tf:
            tf.extractall(d)
    genes = ["TACSTD2", "PTPRC"] + [g for gs in MARKER_SETS.values() for g in gs]
    rows = []
    for mtx in sorted(d.glob("*_matrix.mtx.gz")):
        stem = mtx.name.replace("_matrix.mtx.gz", "")
        label = stem.split("_", 1)[1]  # MPR-1 / Non-MPR-1
        group = "MPR" if label.startswith("MPR") else "NMPR"
        feat = d / f"{stem}_features.tsv.gz"
        bc = d / f"{stem}_barcodes.tsv.gz"
        n_cells = sum(1 for _ in gzip.open(bc, "rt"))
        idx = gene_row_index(feat, genes)
        print(f"GSE291670 {label}: {n_cells} cells, {len(idx)} markers", flush=True)
        vals = stream_mtx_genes(mtx, idx, n_cells)
        expr = pd.DataFrame(vals)
        lin_fine = assign_markers(expr)
        lin = lin_fine.map(_collapse_lineage)
        is_epi = lin.eq("Epithelial")
        counts = _empty_counts()
        for name, n in lin.value_counts().items():
            counts[name] = int(n)
        n_epi = int(is_epi.sum())
        if n_epi:
            tac = float(np.log1p(expr.loc[is_epi, "TACSTD2"]).mean())
            tpos = float((expr.loc[is_epi, "TACSTD2"] > 0).mean())
        else:
            tac = tpos = np.nan
        rows.append(
            {
                "dataset": "GSE291670",
                "cohort": "GSE291670",
                "sample_id": label,
                "patient_id": label,
                "n_cells": int(n_cells),
                "n_epi": n_epi,
                "malignant_TACSTD2_mean_log1p": tac,
                "malignant_TACSTD2_frac_pos": tpos,
                "group": group,
                "pathologic_response": group,
                "timepoint": "Post-treatment surgery",
                "histology": "NSCLC",
                "pd1": "Camrelizumab",
                "annotation": "marker_hierarchy",
                **counts,
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    parts = []
    print("=== GSE207422 ===", flush=True)
    a = build_gse207422()
    a.to_csv(OUT / "composition_GSE207422.tsv", sep="\t", index=False)
    parts.append(a)
    print(a[["sample_id", "group", "n_cells", "Epithelial", "T_NK", "malignant_TACSTD2_mean_log1p"]].to_string(index=False))

    print("=== GSE241934 ===", flush=True)
    b = build_gse241934()
    b.to_csv(OUT / "composition_GSE241934.tsv", sep="\t", index=False)
    parts.append(b)
    print(b[["cohort", "sample_id", "group", "n_cells", "Epithelial", "T_NK", "malignant_TACSTD2_mean_log1p"]].to_string(index=False))

    print("=== GSE291670 ===", flush=True)
    c = build_gse291670()
    c.to_csv(OUT / "composition_GSE291670.tsv", sep="\t", index=False)
    parts.append(c)
    print(c[["sample_id", "group", "n_cells", "Epithelial", "T_NK", "malignant_TACSTD2_mean_log1p"]].to_string(index=False))

    allp = pd.concat(parts, ignore_index=True)
    allp.to_csv(OUT / "composition_all.tsv", sep="\t", index=False)
    print("wrote", OUT / "composition_all.tsv", "n=", len(allp), flush=True)


if __name__ == "__main__":
    main()
