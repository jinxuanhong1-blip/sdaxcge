#!/usr/bin/env python3
"""Build patient-level malignant (or tumor-Epi) UMI-sum matrices.

Public processed GEO only. Writes genes x patient count tables plus a
per-patient inventory. GSE253013 is inventory-only (9.3 GB RDS not loaded).
"""
from __future__ import annotations

import gzip
import json
import os
import shutil
import sys
import tempfile
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

DATA = Path(os.environ.get("SCRNA_PB_DATA", "/tmp/scrna_pb_gsea_data"))
ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "scrna_pseudobulk_gsea_meta"
PB = OUT / "pseudobulk"
PB.mkdir(parents=True, exist_ok=True)

MIN_UMI = 200
MIN_MALIG = 30

LINEAGE = {
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT5", "KRT7", "KRT17", "ELF3", "CDH1", "MUC1"],
    "T/NK": ["CD3D", "CD3E", "CD3G", "TRAC", "CD2", "NKG7", "GNLY", "KLRD1"],
    "B/Plasma": ["CD79A", "CD79B", "MS4A1", "JCHAIN", "MZB1", "IGHM"],
    "Myeloid": ["LYZ", "CD68", "CD14", "C1QA", "C1QB", "FCN1", "ITGAX"],
    "Mast": ["TPSAB1", "TPSB2", "CPA3", "MS4A2"],
    "Endothelial": ["PECAM1", "VWF", "CLDN5", "CDH5", "RAMP2"],
    "Fibroblast": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "TAGLN"],
}
NORMAL_LUNG = [
    "SFTPA1", "SFTPA2", "SFTPB", "SFTPC", "SFTPD", "AGER", "NAPSA",
    "SCGB1A1", "SCGB3A2", "TPPP3", "FOXJ1", "CAPS",
]
PANEL = sorted({g for vs in LINEAGE.values() for g in vs} | set(NORMAL_LUNG) | {"TACSTD2", "PTPRC"})

GSE291670_SAMPLES = [
    ("GSM8839599_MPR-1", "MPR", "MPR-1"),
    ("GSM8839600_MPR-2", "MPR", "MPR-2"),
    ("GSM8839601_MPR-3", "MPR", "MPR-3"),
    ("GSM8839602_Non-MPR-1", "NMPR", "Non-MPR-1"),
    ("GSM8839603_Non-MPR-2", "NMPR", "Non-MPR-2"),
    ("GSM8839604_Non-MPR-3", "NMPR", "Non-MPR-3"),
]


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def present(genes, columns) -> list[str]:
    return [g for g in genes if g in columns]


def log1p_cp10k(counts: pd.DataFrame, total: pd.Series) -> pd.DataFrame:
    tot = np.asarray(total, float)
    tot[tot <= 0] = np.nan
    return pd.DataFrame(
        np.log1p(np.asarray(counts, float) / tot[:, None] * 1e4),
        index=counts.index,
        columns=counts.columns,
    )


def marker_malignant(panel_df: pd.DataFrame, total: pd.Series) -> pd.Series:
    """Epithelial lineage AND normal-lung score ≤ 75th percentile of epithelial."""
    cols = [c for c in panel_df.columns if c in PANEL]
    ln = log1p_cp10k(panel_df[cols], total)
    scores = pd.DataFrame(
        {lin: ln[present(gs, ln.columns)].mean(axis=1) for lin, gs in LINEAGE.items()},
        index=panel_df.index,
    )
    lineage = scores.idxmax(axis=1)
    nl_cols = present(NORMAL_LUNG, ln.columns)
    normal = ln[nl_cols].mean(axis=1) if nl_cols else pd.Series(0.0, index=panel_df.index)
    epi = lineage == "Epithelial"
    cut = float(normal.loc[epi].quantile(0.75)) if epi.any() else 0.0
    return epi & (normal <= cut)


def save_counts(cohort: str, counts: pd.DataFrame, meta: pd.DataFrame) -> None:
    counts = counts.groupby(counts.index.astype(str)).sum()
    counts = counts.loc[counts.sum(axis=1) > 0]
    path = PB / f"{cohort}_counts.tsv.gz"
    counts.to_csv(path, sep="\t", compression="gzip")
    meta.to_csv(PB / f"{cohort}_patient_meta.tsv", sep="\t", index=False)
    log(f"  wrote {path.name} genes={counts.shape[0]} patients={counts.shape[1]}")


def stream_txt_header(path: Path) -> tuple[list[str], object]:
    fh = gzip.open(path, "rt")
    header = fh.readline().rstrip("\n")
    cells = header.split("\t")[1:]
    return cells, fh


def stream_gene_row(fh, n_cells: int):
    for raw in fh:
        tab = raw.find("\t")
        gene = raw[:tab]
        arr = np.fromstring(raw[tab + 1 :].rstrip("\r\n"), sep="\t", dtype=np.float32)
        if arr.size != n_cells:
            raise ValueError(f"{gene}: got {arr.size}, expected {n_cells}")
        yield gene, arr


def accumulate_from_txt(path: Path, cell_to_patient: dict[str, str], keep_cells: np.ndarray | None = None):
    """Stream genes x cells txt; sum UMI into genes x patients for keep_cells."""
    cells, fh = stream_txt_header(path)
    n = len(cells)
    if keep_cells is None:
        keep_cells = np.array([c in cell_to_patient for c in cells])
    patients = sorted({cell_to_patient[c] for c, k in zip(cells, keep_cells) if k})
    pidx = {p: i for i, p in enumerate(patients)}
    col_patient = np.full(n, -1, dtype=np.int32)
    for i, c in enumerate(cells):
        if keep_cells[i] and c in cell_to_patient:
            col_patient[i] = pidx[cell_to_patient[c]]
    use = col_patient >= 0
    use_i = np.where(use)[0]
    assign = col_patient[use]
    log(f"  stream {path.name}: cells={n} kept={int(use.sum())} patients={len(patients)}")
    genes: list[str] = []
    rows: list[np.ndarray] = []
    n_genes = 0
    try:
        for gene, arr in stream_gene_row(fh, n):
            n_genes += 1
            acc = np.zeros(len(patients), dtype=np.float64)
            np.add.at(acc, assign, arr[use_i])
            if acc.sum() > 0:
                genes.append(gene)
                rows.append(acc)
            if n_genes % 4000 == 0:
                log(f"    genes={n_genes} kept_nonzero={len(genes)}")
    finally:
        fh.close()
    mat = pd.DataFrame(np.vstack(rows) if rows else np.zeros((0, len(patients))), index=genes, columns=patients)
    return mat


# ---------------------------------------------------------------------------
# GSE207422
# ---------------------------------------------------------------------------
def run_gse207422() -> dict:
    log("=== GSE207422 ===")
    matrix = DATA / "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"
    meta_x = pd.read_excel(DATA / "GSE207422_NSCLC_scRNAseq_metadata.xlsx")
    meta_x = meta_x.loc[meta_x["Sample"].astype(str).str.startswith("BD_immune")].copy()
    sample_to_pt = dict(zip(meta_x["Sample"].astype(str), meta_x["Patient"].astype(str)))

    cells, fh = stream_txt_header(matrix)
    n = len(cells)
    samples = [c.rsplit("_", 1)[0] for c in cells]
    total = np.zeros(n, dtype=np.float64)
    kept: dict[str, np.ndarray] = {}
    n_genes = 0
    log(f"  pass1 panel+totals cells={n}")
    try:
        for gene, arr in stream_gene_row(fh, n):
            n_genes += 1
            total += arr
            if gene in PANEL:
                kept[gene] = arr.copy()
            if n_genes % 4000 == 0:
                log(f"    pass1 genes={n_genes} panel={len(kept)}")
    finally:
        fh.close()
    panel_df = pd.DataFrame(kept)
    mal = marker_malignant(panel_df, pd.Series(total))
    ok_umi = total >= MIN_UMI
    mal = mal.to_numpy() & ok_umi
    log(f"  malignant-like {int(mal.sum())} / {int(ok_umi.sum())} (UMI>={MIN_UMI}); genes_scanned={n_genes}")

    cell_to_pt = {}
    keep = np.zeros(n, dtype=bool)
    for i, (cell, samp) in enumerate(zip(cells, samples)):
        if mal[i] and samp in sample_to_pt:
            cell_to_pt[cell] = sample_to_pt[samp]
            keep[i] = True
    counts = accumulate_from_txt(matrix, cell_to_pt, keep)
    n_mal = pd.Series(0, index=sorted(set(cell_to_pt.values())), dtype=int)
    for cell, pt in cell_to_pt.items():
        n_mal[pt] += 1
    keep_pt = n_mal[n_mal >= MIN_MALIG].index.tolist()
    counts = counts.loc[:, keep_pt]
    pmeta = meta_x.loc[meta_x["Patient"].astype(str).isin(keep_pt)].copy()
    pmeta = pmeta.rename(columns={"Patient": "patient"})
    pmeta["n_malignant"] = pmeta["patient"].map(n_mal)
    pmeta["malignant_rule"] = "marker_epithelial_low_normal_lung"
    pmeta["cohort"] = "GSE207422"
    save_counts("GSE207422", counts, pmeta)
    return {
        "cohort": "GSE207422",
        "usable": len(keep_pt) >= 6,
        "n_patients_gated": len(keep_pt),
        "n_malignant_cells": int(n_mal.loc[keep_pt].sum()) if keep_pt else 0,
        "n_genes": int(counts.shape[0]),
        "rule": "marker_epithelial_low_normal_lung",
        "reason": "ok" if len(keep_pt) >= 6 else "n_patients<6",
    }


# ---------------------------------------------------------------------------
# GSE131907
# ---------------------------------------------------------------------------
def run_gse131907() -> dict:
    log("=== GSE131907 ===")
    annot = pd.read_csv(DATA / "GSE131907" / "cell_annotation.txt.gz", sep="\t")
    tumor = ~annot["Sample_Origin"].isin(["nLung", "nLN"])
    mal = annot["Cell_subtype"].astype(str).eq("Malignant cells") & tumor
    annot = annot.loc[mal, ["Index", "Sample", "Sample_Origin"]].copy()
    n_mal = annot.groupby("Sample").size()
    keep_pt = n_mal[n_mal >= MIN_MALIG].index
    annot = annot.loc[annot["Sample"].isin(keep_pt)]
    cell_to_pt = dict(zip(annot["Index"].astype(str), annot["Sample"].astype(str)))
    counts = accumulate_from_txt(DATA / "GSE131907" / "raw_UMI_matrix.txt.gz", cell_to_pt)
    counts = counts.loc[:, [c for c in keep_pt if c in counts.columns]]
    pmeta = (
        annot.groupby("Sample")
        .agg(n_malignant=("Index", "size"), Sample_Origin=("Sample_Origin", "first"))
        .reset_index()
        .rename(columns={"Sample": "patient"})
    )
    pmeta["malignant_rule"] = "author_Cell_subtype==Malignant cells; drop nLung/nLN"
    pmeta["cohort"] = "GSE131907"
    save_counts("GSE131907", counts, pmeta)
    return {
        "cohort": "GSE131907",
        "usable": counts.shape[1] >= 6,
        "n_patients_gated": int(counts.shape[1]),
        "n_malignant_cells": int(pmeta["n_malignant"].sum()),
        "n_genes": int(counts.shape[0]),
        "rule": "author_malignant_tumor_met",
        "reason": "ok",
    }


# ---------------------------------------------------------------------------
# GSE241934
# ---------------------------------------------------------------------------
def _stream_mtx(mtx_gz: Path, features: list[str], col_patient: np.ndarray, patients: list[str]) -> pd.DataFrame:
    pidx_max = len(patients)
    acc: dict[int, np.ndarray] = {}
    n_lines = 0
    with gzip.open(mtx_gz, "rt") as fh:
        header = fh.readline()
        if not header.startswith("%%MatrixMarket"):
            raise ValueError(f"not MTX: {mtx_gz}")
        dims = fh.readline()
        while dims.startswith("%"):
            dims = fh.readline()
        n_genes, n_cells, nnz = map(int, dims.split())
        log(f"  MTX {mtx_gz.name}: {n_genes} x {n_cells} nnz={nnz}")
        if n_cells != len(col_patient):
            raise ValueError(f"barcode/meta length {len(col_patient)} != MTX cols {n_cells}")
        for line in fh:
            n_lines += 1
            a, b, c = line.split()
            col = int(b) - 1
            pid = col_patient[col]
            if pid < 0:
                continue
            row = int(a) - 1
            if row not in acc:
                acc[row] = np.zeros(pidx_max, dtype=np.float64)
            acc[row][pid] += float(c)
            if n_lines % 20_000_000 == 0:
                log(f"    MTX lines={n_lines} genes_hit={len(acc)}")
    genes = []
    rows = []
    for i in sorted(acc):
        if acc[i].sum() > 0 and i < len(features):
            genes.append(features[i])
            rows.append(acc[i])
    return pd.DataFrame(np.vstack(rows) if rows else np.zeros((0, pidx_max)), index=genes, columns=patients)


def _gse241934_one(prefix: str, meta_name: str, feat_name: str, bc_name: str, mtx_name: str, tag: str):
    d = DATA / "GSE241934"
    meta = pd.read_csv(d / meta_name, sep="\t", low_memory=False)
    feat = pd.read_csv(d / feat_name, sep="\t", header=None)
    genes = feat.iloc[:, 1].astype(str).tolist() if feat.shape[1] > 1 else feat.iloc[:, 0].astype(str).tolist()
    bc = pd.read_csv(d / bc_name, header=None)[0].astype(str).tolist()
    meta = meta.set_index(meta["cellID"].astype(str))
    epi = meta["major.cell.type"].astype(str).eq("Epi")
    n_epi = meta.loc[epi].groupby("sampleID").size()
    keep_pt = n_epi[n_epi >= MIN_MALIG].index.astype(str)
    patients = sorted(keep_pt)
    pidx = {p: i for i, p in enumerate(patients)}
    col_patient = np.full(len(bc), -1, dtype=np.int32)
    n_keep = 0
    for i, b in enumerate(bc):
        if b in meta.index and bool(epi.get(b, False)):
            pt = str(meta.at[b, "sampleID"])
            if pt in pidx:
                col_patient[i] = pidx[pt]
                n_keep += 1
    log(f"  {tag}: Epi cells used={n_keep} patients={len(patients)}")
    counts = _stream_mtx(d / mtx_name, genes, col_patient, patients)
    pmeta = (
        meta.loc[epi & meta["sampleID"].astype(str).isin(patients)]
        .groupby("sampleID")
        .agg(
            n_malignant=("cellID", "size"),
            Histology=("Histology", "first"),
            Pathological_Response=("Pathological Response", "first"),
            EGFR=("EGFR", "first"),
            PD1=("PD1", "first"),
        )
        .reset_index()
        .rename(columns={"sampleID": "patient"})
    )
    pmeta["slice"] = tag
    pmeta["malignant_rule"] = "author_major.cell.type==Epi"
    pmeta["cohort"] = "GSE241934"
    return counts, pmeta


def run_gse241934() -> dict:
    log("=== GSE241934 ===")
    c1, m1 = _gse241934_one("IIT", "IIT_Meta.txt.gz", "IIT_features.tsv.gz", "IIT_barcodes.tsv.gz", "IIT_Matrix.mtx.gz", "IIT")
    c2, m2 = _gse241934_one("RWC", "Real_Meta.txt.gz", "RWC_features.tsv.gz", "RWC_barcodes.tsv.gz", "Real_Matrix.mtx.gz", "RWC")
    # patient IDs do not overlap between IIT and RWC
    genes = sorted(set(c1.index) | set(c2.index))
    counts = pd.concat(
        [c1.reindex(genes, fill_value=0), c2.reindex(genes, fill_value=0)],
        axis=1,
    )
    pmeta = pd.concat([m1, m2], ignore_index=True)
    save_counts("GSE241934", counts, pmeta)
    return {
        "cohort": "GSE241934",
        "usable": counts.shape[1] >= 6,
        "n_patients_gated": int(counts.shape[1]),
        "n_malignant_cells": int(pmeta["n_malignant"].sum()),
        "n_genes": int(counts.shape[0]),
        "rule": "author_Epi_IIT+RWC",
        "reason": "ok",
        "n_IIT": int(m1.shape[0]),
        "n_RWC": int(m2.shape[0]),
    }


# ---------------------------------------------------------------------------
# GSE205335
# ---------------------------------------------------------------------------
def _ungzip_until_rds(src: Path, dest: Path) -> None:
    current = src
    tmp_paths = []
    for _ in range(3):
        with open(current, "rb") as fh:
            magic = fh.read(2)
        if magic == b"\x1f\x8b":
            nxt = dest.with_suffix(dest.suffix + f".pass{len(tmp_paths)}")
            with gzip.open(current, "rb") as zin, open(nxt, "wb") as zout:
                shutil.copyfileobj(zin, zout, 16 * 1024 * 1024)
            tmp_paths.append(nxt)
            current = nxt
            continue
        break
    if current != dest:
        shutil.copyfile(current, dest)
    for p in tmp_paths:
        if p.exists() and p != dest:
            p.unlink(missing_ok=True)


def _soft_gse205335() -> pd.DataFrame:
    recs = []
    cur = None
    with gzip.open(DATA / "GSE205335_family.soft.gz", "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("^SAMPLE = "):
                if cur:
                    recs.append(cur)
                cur = {"gsm": line.split(" = ", 1)[1].strip()}
            elif cur is not None:
                if line.startswith("!Sample_title = "):
                    cur["title"] = line.split(" = ", 1)[1].strip()
                elif line.startswith("!Sample_description = ") and "description" not in cur:
                    cur["description"] = line.split(" = ", 1)[1].strip()
                elif line.startswith("!Sample_characteristics_ch1 = "):
                    v = line.split(" = ", 1)[1].strip()
                    if ": " in v:
                        k, x = v.split(": ", 1)
                        cur[k.strip()] = x.strip()
        if cur:
            recs.append(cur)
    return pd.DataFrame(recs)


def run_gse205335() -> dict:
    log("=== GSE205335 ===")
    import rdata
    from scipy import sparse

    ident = pd.read_csv(DATA / "GSE205335_Lung_IO_CellIdentity.txt.gz", sep="\t")
    soft = _soft_gse205335()
    soft["orig_key"] = soft["description"].str.replace("_", "-", regex=False)
    ident["orig_key"] = ident["orig.ident"].str.replace(r"-[35]P$", "", regex=True)
    smap = ident[["orig.ident", "orig_key"]].drop_duplicates().merge(soft, on="orig_key", how="left")
    ident = ident.merge(smap[["orig.ident", "patient", "tissue", "recist", "cancer subtype"]], on="orig.ident", how="left")
    normal = {"Normal Brain", "Normal LN", "Normal Lung"}
    mal = (
        ident["lineage.sub"].eq("Malignant cells")
        & ~ident["tissue"].isin(normal)
        & ident["patient"].notna()
    )
    ident_m = ident.loc[mal].copy()
    n_mal = ident_m.groupby("patient").size()
    keep_pt = n_mal[n_mal >= MIN_MALIG].index
    ident_m = ident_m.loc[ident_m["patient"].isin(keep_pt)]
    barcode_to_pt = dict(zip(ident_m["barcode"].astype(str), ident_m["patient"].astype(str)))
    log(f"  author malignant tumor cells={len(ident_m)} patients={len(keep_pt)}")

    with tempfile.TemporaryDirectory(prefix="gse205335-") as tmp:
        rds_path = Path(tmp) / "matrix.rds"
        log("  decompressing RDS...")
        _ungzip_until_rds(DATA / "GSE205335_Lung_IO_UMI_matrix.rds.gz", rds_path)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message='Missing constructor for R class "dgCMatrix"')
            obj = rdata.read_rds(rds_path)
    genes = np.asarray(obj.Dimnames[0], dtype=str)
    barcodes = np.asarray(obj.Dimnames[1], dtype=str)
    matrix = sparse.csc_matrix((obj.x, obj.i, obj.p), shape=tuple(obj.Dim), copy=False)
    log(f"  RDS {matrix.shape[0]} genes x {matrix.shape[1]} cells")
    patients = sorted(keep_pt.astype(str))
    pidx = {p: i for i, p in enumerate(patients)}
    cols_by_pt: dict[str, list[int]] = {p: [] for p in patients}
    for j, bc in enumerate(barcodes):
        pt = barcode_to_pt.get(str(bc))
        if pt is not None:
            cols_by_pt[pt].append(j)
    acc = np.zeros((matrix.shape[0], len(patients)), dtype=np.float64)
    n_used = 0
    for pt, cols in cols_by_pt.items():
        if not cols:
            continue
        acc[:, pidx[pt]] = np.asarray(matrix[:, cols].sum(axis=1)).ravel()
        n_used += len(cols)
    log(f"  summed columns={n_used}")
    counts = pd.DataFrame(acc, index=genes, columns=patients)
    pmeta = (
        ident_m.groupby("patient")
        .agg(
            n_malignant=("barcode", "size"),
            recist=("recist", "first"),
            tissue=("tissue", lambda s: ",".join(sorted(set(map(str, s))))),
            cancer_subtype=("cancer subtype", "first"),
        )
        .reset_index()
    )
    pmeta["malignant_rule"] = "author_lineage.sub==Malignant cells; drop normal tissues"
    pmeta["cohort"] = "GSE205335"
    save_counts("GSE205335", counts, pmeta)
    return {
        "cohort": "GSE205335",
        "usable": counts.shape[1] >= 6,
        "n_patients_gated": int(counts.shape[1]),
        "n_malignant_cells": int(pmeta["n_malignant"].sum()),
        "n_genes": int(counts.shape[0]),
        "rule": "author_malignant_tumor",
        "reason": "ok",
    }


# ---------------------------------------------------------------------------
# GSE291670
# ---------------------------------------------------------------------------
def run_gse291670() -> dict:
    log("=== GSE291670 ===")
    import scipy.io

    d = DATA / "GSE291670"
    parts = []
    metas = []
    for prefix, resp, name in GSE291670_SAMPLES:
        mtx = scipy.io.mmread(str(d / f"{prefix}_matrix.mtx.gz")).tocsr().astype(np.float32)
        with gzip.open(d / f"{prefix}_barcodes.tsv.gz", "rt") as fh:
            barcodes = [line.strip().split("\t")[0] for line in fh]
        genes = []
        with gzip.open(d / f"{prefix}_features.tsv.gz", "rt") as fh:
            for line in fh:
                parts_l = line.rstrip("\n").split("\t")
                genes.append(parts_l[1] if len(parts_l) > 1 else parts_l[0])
        if mtx.shape[0] != len(genes) and mtx.shape[1] == len(genes):
            mtx = mtx.T.tocsr()
        gidx = {}
        for i, g in enumerate(genes):
            if g not in gidx:
                gidx[g] = i
        total = np.asarray(mtx.sum(axis=0)).ravel()
        panel_cols = {}
        for g in PANEL:
            if g in gidx:
                panel_cols[g] = np.asarray(mtx[gidx[g], :].todense()).ravel()
        panel_df = pd.DataFrame(panel_cols)
        mal = marker_malignant(panel_df, pd.Series(total)).to_numpy() & (total >= MIN_UMI)
        n_mal = int(mal.sum())
        log(f"  {name}: cells={mtx.shape[1]} malignant-like={n_mal}")
        metas.append({"patient": name, "response": resp, "n_malignant": n_mal, "n_cells": int(mtx.shape[1])})
        if n_mal < MIN_MALIG:
            continue
        vec = np.asarray(mtx[:, mal].sum(axis=1)).ravel()
        parts.append(pd.Series(vec, index=genes, name=name))
    if not parts:
        counts = pd.DataFrame()
    else:
        counts = pd.concat(parts, axis=1).fillna(0)
    pmeta = pd.DataFrame(metas)
    pmeta["malignant_rule"] = "marker_epithelial_low_normal_lung"
    pmeta["cohort"] = "GSE291670"
    gated = pmeta.loc[pmeta["n_malignant"] >= MIN_MALIG, "patient"].tolist()
    if len(gated):
        counts = counts.loc[:, [c for c in gated if c in counts.columns]]
    save_counts("GSE291670", counts, pmeta)
    return {
        "cohort": "GSE291670",
        "usable": counts.shape[1] >= 6,
        "n_patients_gated": int(counts.shape[1]),
        "n_malignant_cells": int(pmeta.loc[pmeta["n_malignant"] >= MIN_MALIG, "n_malignant"].sum()),
        "n_genes": int(counts.shape[0]),
        "rule": "marker_epithelial_low_normal_lung",
        "reason": "ok" if counts.shape[1] >= 6 else f"n_patients={counts.shape[1]}<6",
    }


def run_gse253013_inventory() -> dict:
    log("=== GSE253013 inventory-only ===")
    titles = []
    path = DATA / "GSE253013_series_matrix.txt.gz"
    if path.exists():
        with gzip.open(path, "rt", errors="replace") as fh:
            for line in fh:
                if line.startswith("!Sample_title"):
                    titles = [t.strip().strip('"') for t in line.split("\t")[1:]]
                    break
    patients = sorted({t.split("_")[0] for t in titles if t.startswith("MRC")})
    note = (
        "9.3 GB Garnett RDS not loaded (16 GB RAM). Public series: 9 treatment-naive "
        "LUAD (tumor+adjacent), no ICI/MPR labels. Prior public marker extract: "
        "MRC004/007 <30 malignant-like cells → ≤7 gated patients, below meta n gate."
    )
    inv = pd.DataFrame({"patient": patients, "cohort": "GSE253013", "gsea_run": False, "note": note})
    inv.to_csv(PB / "GSE253013_patient_meta.tsv", sep="\t", index=False)
    return {
        "cohort": "GSE253013",
        "usable": False,
        "n_patients_gated": 0,
        "n_malignant_cells": None,
        "n_genes": 0,
        "rule": "not_run",
        "reason": "rds_9.3GB_not_loaded_and_n_patients_below_meta_gate",
        "n_patients_series": len(patients),
        "note": note,
    }


def main() -> None:
    inventory = []
    inventory.append(run_gse207422())
    inventory.append(run_gse131907())
    inventory.append(run_gse241934())
    inventory.append(run_gse205335())
    inventory.append(run_gse291670())
    inventory.append(run_gse253013_inventory())
    (OUT / "tables").mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(inventory)
    df.to_csv(OUT / "tables" / "cohort_inventory.tsv", sep="\t", index=False)
    (OUT / "pseudobulk_inventory.json").write_text(json.dumps(inventory, indent=2) + "\n")
    log("done inventory")
    print(df.to_string(index=False))


if __name__ == "__main__":
    (OUT / "tables").mkdir(parents=True, exist_ok=True)
    main()
