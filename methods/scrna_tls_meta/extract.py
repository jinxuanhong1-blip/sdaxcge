#!/usr/bin/env python3
"""Stream selected genes from public lung scRNA matrices to per-cell tables.

Does not load a full dense matrix. Author labels preferred; marker lineage
only when barcodes are unlabeled. Filters are not tuned to a claimed ρ.
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import tarfile
from pathlib import Path

import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from genes import LINEAGE, NORMAL_LUNG, PANEL

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "scrna_tls_meta"
EXTRACTED = DATA / "extracted"


def _norm_symbol(g: str) -> str:
    return str(g).split(".")[0].strip()


def marker_lineage(expr: dict[str, np.ndarray], n: int) -> np.ndarray:
    names = list(LINEAGE)
    scores = []
    for genes in LINEAGE.values():
        mats = [np.log1p(expr[g].astype(np.float32)) for g in genes if g in expr]
        if mats:
            scores.append(np.mean(np.vstack(mats), axis=0))
        else:
            scores.append(np.zeros(n, dtype=np.float32))
    scores = np.vstack(scores)
    best = scores.argmax(axis=0)
    best_val = scores.max(axis=0)
    second = np.partition(scores, -2, axis=0)[-2]
    assigned = np.array(names, dtype=object)[best]
    assigned[(best_val < 0.15) | ((best_val - second) < 0.05)] = "other"
    return assigned


def malignant_like(expr: dict[str, np.ndarray], lineage: np.ndarray) -> np.ndarray:
    n = len(lineage)
    mats = [np.log1p(expr[g].astype(np.float32)) for g in NORMAL_LUNG if g in expr]
    if not mats:
        return lineage == "epithelial"
    nl = np.max(np.vstack(mats), axis=0)
    return (lineage == "epithelial") & (nl <= 0.05)


def stream_gene_cell_tsv(
    path: Path,
    genes: set[str],
    skip_index_name: bool = True,
) -> tuple[list[str], dict[str, np.ndarray], np.ndarray]:
    """Genes × cells TSV.gz. Accumulate library size while keeping selected rows."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if skip_index_name:
            cells = header[1:]
        else:
            cells = header[1:] if header[0] in {"Index", "gene", "Gene", ""} else header
        n = len(cells)
        lib = np.zeros(n, dtype=np.float64)
        keep: dict[str, np.ndarray] = {}
        for line in fh:
            if not line:
                continue
            gene, rest = line.split("\t", 1)
            gene = _norm_symbol(gene)
            vals = np.fromstring(rest, sep="\t", dtype=np.float32)
            if vals.size != n:
                continue
            lib += vals
            if gene in genes:
                keep[gene] = vals
    return cells, keep, lib.astype(np.float32)


def stream_mtx(
    mtx: Path,
    features: Path,
    barcodes: Path,
    genes: set[str],
    lib_from_meta: np.ndarray | None = None,
) -> tuple[list[str], dict[str, np.ndarray], np.ndarray]:
    feats = pd.read_csv(features, sep="\t", header=None, dtype=str)
    if feats.shape[1] >= 2:
        symbols = feats[1].map(_norm_symbol)
        alt = feats[0].map(_norm_symbol)
        symbols = np.where(symbols.isin(genes), symbols, alt)
    else:
        symbols = feats[0].map(_norm_symbol).to_numpy()
    symbols = np.asarray(symbols)
    want_rows = {i for i, s in enumerate(symbols) if s in genes}
    bcs = pd.read_csv(barcodes, sep="\t", header=None, dtype=str)[0].tolist()
    n = len(bcs)
    keep = {s: np.zeros(n, dtype=np.float32) for s in set(symbols[list(want_rows)])} if want_rows else {}
    lib = np.zeros(n, dtype=np.float64) if lib_from_meta is None else None
    opener = gzip.open if str(mtx).endswith(".gz") else open
    with opener(mtx, "rt") as fh:
        for line in fh:
            if line.startswith("%"):
                continue
            # dims
            break
        for line in fh:
            parts = line.split()
            if len(parts) < 3:
                continue
            r = int(parts[0]) - 1
            c = int(parts[1]) - 1
            v = float(parts[2])
            if lib is not None:
                lib[c] += v
            if r in want_rows:
                keep[symbols[r]][c] += v
    if lib is None:
        lib = lib_from_meta.astype(np.float32)
    else:
        lib = lib.astype(np.float32)
    return bcs, keep, lib


def write_cells(
    dest: Path,
    cohort: str,
    cells: list[str],
    patient: list[str] | np.ndarray,
    sample: list[str] | np.ndarray,
    tissue: list[str] | np.ndarray,
    lineage: np.ndarray,
    is_mal: np.ndarray,
    lib: np.ndarray,
    expr: dict[str, np.ndarray],
    extra: dict[str, np.ndarray] | None = None,
) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(
        {
            "cohort": cohort,
            "cell": cells,
            "patient": patient,
            "sample": sample,
            "tissue": tissue,
            "lineage": lineage,
            "is_malignant": is_mal.astype(bool),
            "n_umi": lib,
        }
    )
    for g, arr in expr.items():
        df[g] = arr
    if extra:
        for k, v in extra.items():
            df[k] = v
    df.to_parquet(dest, index=False)
    print(f"wrote {dest} n={len(df)} genes={sorted(expr)}", flush=True)


def extract_gse131907(data: Path, out: Path) -> None:
    annot = pd.read_csv(
        data / "GSE131907" / "GSE131907_Lung_Cancer_cell_annotation.txt.gz",
        sep="\t",
    )
    series_path = data / "GSE131907" / "GSE131907_series_matrix.txt.gz"
    sample_to_patient: dict[str, str] = {}
    cur_title = None
    with gzip.open(series_path, "rt") as fh:
        for line in fh:
            if line.startswith("!Sample_title"):
                titles = [x.strip().strip('"') for x in line.split("\t")[1:]]
                titles_iter = titles
            elif line.startswith("!Sample_characteristics_ch1"):
                fields = [x.strip().strip('"') for x in line.split("\t")[1:]]
                if fields and fields[0].lower().startswith("patient id"):
                    for t, f in zip(titles_iter, fields):
                        sample_to_patient[t] = f.split(":", 1)[-1].strip()
    cells, expr, lib = stream_gene_cell_tsv(
        data / "GSE131907" / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz",
        set(PANEL),
    )
    annot = annot.set_index("Index").reindex(cells)
    sample = annot["Sample"].fillna("NA").astype(str).to_numpy()
    origin = annot["Sample_Origin"].fillna("").astype(str).to_numpy()
    ctype = annot["Cell_type"].fillna("").astype(str).to_numpy()
    subtype = annot["Cell_subtype"].fillna("").astype(str).to_numpy()
    patient = [sample_to_patient.get(s, s) for s in sample]
    lin_map = {
        "Epithelial cells": "epithelial",
        "T lymphocytes": "T",
        "NK cells": "NK",
        "B lymphocytes": "B",
        "Myeloid cells": "myeloid",
        "MAST cells": "myeloid",
        "Fibroblasts": "fibroblast",
        "Endothelial cells": "endothelial",
    }
    lineage = np.array([lin_map.get(x, "other") for x in ctype], dtype=object)
    tumor_origins = {"tLung", "tL/B", "mLN", "mBrain", "PE"}
    tissue = np.array(
        ["tumor" if o in tumor_origins else ("normal" if o in {"nLung", "nLN"} else "other") for o in origin],
        dtype=object,
    )
    mal_sub = {"Malignant cells", "tS1", "tS2", "tS3"}
    is_mal = (lineage == "epithelial") & np.array([s in mal_sub for s in subtype])
    # tumor epithelial without a named malignant subtype still counts as epithelial
    write_cells(
        out / "GSE131907.parquet",
        "GSE131907",
        cells,
        patient,
        sample,
        tissue,
        lineage,
        is_mal,
        lib,
        expr,
        extra={"cell_type": ctype, "cell_subtype": subtype, "sample_origin": origin},
    )


def extract_gse207422(data: Path, out: Path) -> None:
    meta = pd.read_excel(data / "GSE207422" / "GSE207422_NSCLC_scRNAseq_metadata.xlsx")
    cells, expr, lib = stream_gene_cell_tsv(
        data / "GSE207422" / "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz",
        set(PANEL),
    )
    sample = [c.rsplit("_", 1)[0] for c in cells]
    meta = meta.copy()
    meta["Sample"] = meta["Sample"].astype(str)
    smap = meta.set_index("Sample")
    patient = [str(smap.loc[s, "Patient"]) if s in smap.index else s for s in sample]
    resource = [str(smap.loc[s, "Resource"]) if s in smap.index else "" for s in sample]
    tissue = np.array(["tumor"] * len(cells), dtype=object)
    n = len(cells)
    lineage = marker_lineage(expr, n)
    is_mal = malignant_like(expr, lineage)
    write_cells(
        out / "GSE207422.parquet",
        "GSE207422",
        cells,
        patient,
        sample,
        tissue,
        lineage,
        is_mal,
        lib,
        expr,
        extra={"resource": np.array(resource, dtype=object)},
    )


def _map_241934_lineage(major: str, fine: str = "") -> str:
    fine_l = str(fine).lower()
    if "plasma" in fine_l:
        return "plasma"
    t = str(major).strip().lower()
    if t in {"epi", "epithelial", "cancer", "tumor", "malig", "malignant"}:
        return "epithelial"
    if t in {"b", "bcell", "b_cell", "b-cell"} or t.startswith("b"):
        return "B"
    if t in {"t", "tcell", "t_cell"} or t.startswith("t"):
        return "T"
    if t in {"nk", "nkcell"}:
        return "NK"
    if any(x in t for x in ("mye", "mono", "mac", "dc", "neut", "mast")):
        return "myeloid"
    if "fibro" in t:
        return "fibroblast"
    if "endo" in t:
        return "endothelial"
    return "other"


def extract_gse241934_split(data: Path, out: Path, split: str) -> None:
    d = data / "GSE241934"
    if split == "IIT":
        mtx, feat, bc, meta_p = (
            d / "GSE241934_IIT_Matrix.mtx.gz",
            d / "GSE241934_IIT_features.tsv.gz",
            d / "GSE241934_IIT_barcodes.tsv.gz",
            d / "GSE241934_IIT_Meta.txt.gz",
        )
        cohort = "GSE241934_IIT"
    else:
        mtx, feat, bc, meta_p = (
            d / "GSE241934_Real_Matrix.mtx.gz",
            d / "GSE241934_RWC_features.tsv.gz",
            d / "GSE241934_RWC_barcodes.tsv.gz",
            d / "GSE241934_Real_Meta.txt.gz",
        )
        cohort = "GSE241934_RWC"
    meta = pd.read_csv(meta_p, sep="\t")
    # barcodes may be raw 10x or already prefixed
    bcs_file = pd.read_csv(bc, sep="\t", header=None, dtype=str)[0].tolist()
    if "cellID" in meta.columns:
        key = "cellID"
    elif "barcode" in meta.columns:
        key = "barcode"
    else:
        key = meta.columns[0]
    # align matrix columns to metadata order if possible
    meta_ids = meta[key].astype(str).tolist()
    if set(bcs_file) == set(meta_ids) or all(b in set(meta_ids) for b in bcs_file[:5]):
        order = {b: i for i, b in enumerate(bcs_file)}
        # extract in matrix order then reindex? easier: extract then map
        pass
    lib_meta = None
    if "nCount_RNA" in meta.columns:
        # will reindex after extract
        lib_lookup = dict(zip(meta[key].astype(str), meta["nCount_RNA"].astype(float)))
    else:
        lib_lookup = {}
    cells, expr, lib = stream_mtx(mtx, feat, bc, set(PANEL), lib_from_meta=None)
    # if metadata uses prefixed IDs, try to join
    def _find_row(cell: str):
        if cell in lib_lookup:
            return cell
        # try suffix
        return None

    major_col = None
    for c in ("major_cell_type", "major.cell.type", "cell.type", "cell_type"):
        if c in meta.columns:
            major_col = c
            break
    fine_col = "cell.type" if "cell.type" in meta.columns else None
    meta_idx = meta.set_index(meta[key].astype(str))
    # also index by barcode suffix
    if "cellID" in meta.columns:
        meta["_bc"] = meta["cellID"].astype(str).str.split("_").str[-1]
    patient = []
    sample = []
    lineage = []
    is_mal = []
    lib2 = []
    major_vals = []
    sid_col = "sampleID" if "sampleID" in meta.columns else ("orig.ident" if "orig.ident" in meta.columns else key)
    for cell in cells:
        row = None
        if cell in meta_idx.index:
            row = meta_idx.loc[cell]
        else:
            # try match last token
            tok = cell.split("_")[-1]
            hits = meta[meta[key].astype(str).str.endswith(tok)] if "cellID" in meta.columns else None
            if hits is not None and len(hits) == 1:
                row = hits.iloc[0]
        if row is None:
            patient.append("unknown")
            sample.append(cell.split("_")[0] if "_" in cell else cell)
            lineage.append("other")
            is_mal.append(False)
            lib2.append(lib[len(patient) - 1] if False else 0)
            major_vals.append("")
            continue
        if isinstance(row, pd.DataFrame):
            row = row.iloc[0]
        patient.append(str(row[sid_col]))
        sample.append(str(row[sid_col]))
        maj = str(row[major_col]) if major_col else ""
        fine = str(row[fine_col]) if fine_col else ""
        major_vals.append(maj)
        lin = _map_241934_lineage(maj, fine)
        lineage.append(lin)
        is_mal.append(lin == "epithelial")
        lib2.append(float(row["nCount_RNA"]) if "nCount_RNA" in row.index else float("nan"))
    lib_use = np.array(lib2, dtype=np.float32)
    lib_use = np.where(np.isfinite(lib_use) & (lib_use > 0), lib_use, lib)
    n = len(cells)
    # refine malignant with normal-lung markers among epithelial
    lin_arr = np.array(lineage, dtype=object)
    is_mal_arr = malignant_like(expr, lin_arr)
    # if author already called epithelial, keep author epithelial as malignant-like
    # unless normal-lung high
    write_cells(
        out / f"{cohort}.parquet",
        cohort,
        cells,
        patient,
        sample,
        np.array(["tumor"] * n, dtype=object),
        lin_arr,
        is_mal_arr,
        lib_use,
        expr,
        extra={"author_major": np.array(major_vals, dtype=object)},
    )


def extract_gse148071(data: Path, out: Path) -> None:
    tar_p = data / "GSE148071" / "GSE148071_RAW.tar"
    frames = []
    with tarfile.open(tar_p, "r") as tf:
        members = [m for m in tf.getmembers() if m.name.endswith(".txt.gz") or m.name.endswith(".txt")]
        for m in members:
            pat = re.search(r"(P\d+)", m.name)
            patient = pat.group(1) if pat else Path(m.name).stem
            import io

            raw = tf.extractfile(m).read()
            if m.name.endswith(".gz"):
                text_fh = io.TextIOWrapper(gzip.GzipFile(fileobj=io.BytesIO(raw)), encoding="utf-8")
            else:
                text_fh = io.TextIOWrapper(io.BytesIO(raw), encoding="utf-8")
            first = text_fh.readline().rstrip("\n")
            cols = first.split("\t")
            # Wu GSE148071: header is barcodes only (no gene column).
            if cols and cols[0] not in {"", "Index", "gene", "Gene", "GENE", "symbol"}:
                cells = cols
            else:
                cells = cols[1:]
            n = len(cells)
            lib = np.zeros(n, dtype=np.float64)
            expr: dict[str, np.ndarray] = {}
            for line in text_fh:
                gene, rest = line.split("\t", 1)
                gene = _norm_symbol(gene)
                vals = np.fromstring(rest, sep="\t", dtype=np.float32)
                if vals.size != n:
                    continue
                lib += vals
                if gene in PANEL:
                    expr[gene] = vals
            lin = marker_lineage(expr, n)
            is_mal = malignant_like(expr, lin)
            df = pd.DataFrame(
                {
                    "cohort": "GSE148071",
                    "cell": [f"{patient}_{c}" for c in cells],
                    "patient": patient,
                    "sample": patient,
                    "tissue": "tumor",
                    "lineage": lin,
                    "is_malignant": is_mal,
                    "n_umi": lib.astype(np.float32),
                }
            )
            for g, arr in expr.items():
                df[g] = arr
            frames.append(df)
            print(f"  GSE148071 {patient} n={n}", flush=True)
    dest = out / "GSE148071.parquet"
    pd.concat(frames, ignore_index=True).to_parquet(dest, index=False)
    print(f"wrote {dest}", flush=True)


def extract_gse154826(data: Path, out: Path) -> None:
    lead = data / "leader"
    cm = pd.read_csv(lead / "cell_metadata.csv")
    s1 = pd.read_csv(lead / "table_s1_sample_table.csv")
    annots = pd.read_csv(lead / "annots_list.csv")
    cm["barcode"] = cm["cell_ID"].astype(str).str.split("_", n=1).str[1]
    cm["sample_ID"] = cm["sample_ID"].astype(str)
    s1["sample_ID"] = s1["sample_ID"].astype(str)
    cm = cm.merge(s1, on="sample_ID", how="left")
    annots["cluster"] = annots["cluster"].astype(int)
    cm["cluster_ID"] = cm["cluster_ID"].astype(int)
    cm = cm.merge(annots, left_on="cluster_ID", right_on="cluster", how="left")
    lin_map = {
        "T": "T",
        "NK": "NK",
        "B&plasma": "B",
        "MNP": "myeloid",
        "mast": "myeloid",
        "pDC": "myeloid",
        "epi_endo_fibro_doublet": "gate",
    }
    # split true B vs plasma via sub_lineage
    lineage = []
    for lin, sub in zip(cm["lineage"].astype(str), cm["sub_lineage"].fillna("").astype(str)):
        if lin == "B&plasma":
            lineage.append("B" if sub == "B" else "plasma")
        else:
            lineage.append(lin_map.get(lin, "other"))
    cm["_lineage"] = lineage
    cm["_tissue"] = np.where(cm["tissue"].astype(str).str.lower().eq("tumor"), "tumor", "normal")
    # stream each batch tar
    batch_dir = data / "GSE154826"
    tars = sorted(batch_dir.glob("GSE154826_amp_batch_ID_*.tar.gz"))
    frames = []
    for tar_p in tars:
        m = re.search(r"amp_batch_ID_(\d+)", tar_p.name)
        batch = int(m.group(1)) if m else -1
        sub = cm[cm["amp_batch_ID"] == batch]
        if sub.empty:
            continue
        want_bc = set(sub["barcode"].astype(str))
        with tarfile.open(tar_p, "r:gz") as tf:
            names = tf.getnames()
            fname = [n for n in names if n.endswith("features.tsv") or n.endswith("genes.tsv")][0]
            bname = [n for n in names if n.endswith("barcodes.tsv")][0]
            mname = [n for n in names if n.endswith("matrix.mtx")][0]
            feats = pd.read_csv(tf.extractfile(fname), sep="\t", header=None, dtype=str)
            if feats.shape[1] >= 2:
                symbols = feats[1].map(_norm_symbol).to_numpy()
            else:
                symbols = feats[0].map(_norm_symbol).to_numpy()
            barcodes = pd.read_csv(tf.extractfile(bname), header=None, dtype=str)[0].str.replace("-1$", "", regex=True)
            # keep original too
            barcodes_raw = pd.read_csv(tf.extractfile(bname), header=None, dtype=str)[0]
            # reopen mtx
            fh = tf.extractfile(mname)
            while True:
                line = fh.readline().decode()
                if line.startswith("%"):
                    continue
                dims = [int(x) for x in line.split()]
                break
            n_cell = int(dims[1])
            keep_idx = {}
            for i, bc in enumerate(barcodes_raw.tolist()):
                bc2 = bc.replace("-1", "") if False else bc
                if bc in want_bc:
                    keep_idx[i] = bc
                elif bc2 in want_bc:
                    keep_idx[i] = bc2
                else:
                    # barcode may include -1
                    if str(bc).split("-")[0] in want_bc:
                        keep_idx[i] = str(bc)
            # rebuild want with exact barcode strings from matrix
            # map annotated barcodes
            ann_bc = set(sub["barcode"].astype(str))
            col_keep = {}
            for i, bc in enumerate(barcodes_raw.astype(str).tolist()):
                if bc in ann_bc:
                    col_keep[i] = bc
                elif bc.replace("-1", "") in ann_bc:
                    col_keep[i] = bc.replace("-1", "")
                elif (bc + "-1") in ann_bc:
                    col_keep[i] = bc + "-1"
            want_rows = {i: symbols[i] for i in range(len(symbols)) if symbols[i] in PANEL}
            gene_acc = {g: {} for g in set(want_rows.values())}
            lib_acc: dict[int, float] = {i: 0.0 for i in col_keep}
            for line in fh:
                parts = line.decode().split()
                if len(parts) < 3:
                    continue
                r = int(parts[0]) - 1
                c = int(parts[1]) - 1
                v = float(parts[2])
                if c in lib_acc:
                    lib_acc[c] += v
                    if r in want_rows:
                        gene_acc[want_rows[r]][c] = gene_acc[want_rows[r]].get(c, 0.0) + v
        if not col_keep:
            print(f"  GSE154826 batch {batch}: no barcode overlap", flush=True)
            continue
        cols = sorted(col_keep)
        bcs = [col_keep[c] for c in cols]
        sub_i = sub.set_index("barcode")
        rows = []
        for c, bc in zip(cols, bcs):
            if bc not in sub_i.index:
                continue
            r = sub_i.loc[bc]
            if isinstance(r, pd.DataFrame):
                r = r.iloc[0]
            rec = {
                "cohort": "GSE154826",
                "cell": str(r["cell_ID"]),
                "patient": str(r["patient_ID"]),
                "sample": str(r["sample_ID"]),
                "tissue": str(r["_tissue"]),
                "lineage": str(r["_lineage"]),
                "is_malignant": bool(r["_lineage"] == "gate"),
                "n_umi": float(lib_acc.get(c, 0.0)),
                "prep": str(r.get("prep", "")),
                "disease": str(r.get("disease", "")),
                "author_lineage": str(r.get("lineage", "")),
                "sub_lineage": str(r.get("sub_lineage", "")),
            }
            for g in gene_acc:
                rec[g] = float(gene_acc[g].get(c, 0.0))
            rows.append(rec)
        if rows:
            frames.append(pd.DataFrame(rows))
            print(f"  GSE154826 batch {batch}: {len(rows)} annotated cells", flush=True)
    dest = out / "GSE154826.parquet"
    if frames:
        pd.concat(frames, ignore_index=True).to_parquet(dest, index=False)
        print(f"wrote {dest}", flush=True)
    else:
        print("GSE154826: no cells extracted", flush=True)


def extract_gse253013_from_assembled(data: Path, out: Path) -> None:
    """Build a per-cell parquet from the streamed RDS gene_panel.npz."""
    ex = data / "GSE253013" / "extracted"
    expr_p = ex / "gene_panel.npz"
    meta_p = ex / "cell_metadata.tsv"
    if not expr_p.exists() or not meta_p.exists():
        print("GSE253013 assembled files missing; skip", flush=True)
        return
    packed = np.load(expr_p)
    expr = {k: packed[k] for k in packed.files}
    meta = pd.read_csv(meta_p, sep="\t")
    n = len(meta)
    lineage = marker_lineage(expr, n)
    is_mal = malignant_like(expr, lineage)
    # author labels if present
    if "cell_type" in meta.columns:
        author = meta["cell_type"].astype(str).to_numpy()
    else:
        author = np.array([""] * n, dtype=object)
    patient = meta["patient"].astype(str) if "patient" in meta.columns else meta.iloc[:, 0].astype(str)
    tissue_raw = meta["tissue_raw"].astype(str) if "tissue_raw" in meta.columns else pd.Series(["T"] * n)
    tissue = np.where(tissue_raw.str.upper().isin({"T", "TUMOR"}), "tumor", "normal")
    lib = meta["nCount_RNA"].to_numpy(dtype=np.float32) if "nCount_RNA" in meta.columns else np.zeros(n, np.float32)
    if "total" in meta.columns:
        lib = np.where(lib > 0, lib, meta["total"].to_numpy(dtype=np.float32))
    cells = meta["barcode"].astype(str) if "barcode" in meta.columns else [f"c{i}" for i in range(n)]
    sample = meta["library_id"].astype(str) if "library_id" in meta.columns else patient
    write_cells(
        out / "GSE253013.parquet",
        "GSE253013",
        list(cells),
        list(patient),
        list(sample),
        tissue,
        lineage,
        is_mal,
        lib,
        {k: v for k, v in expr.items() if k in PANEL or k in expr},
        extra={"author_cell_type": author},
    )


COHORT_FUNCS = {
    "GSE131907": extract_gse131907,
    "GSE207422": extract_gse207422,
    "GSE241934_IIT": lambda d, o: extract_gse241934_split(d, o, "IIT"),
    "GSE241934_RWC": lambda d, o: extract_gse241934_split(d, o, "RWC"),
    "GSE148071": extract_gse148071,
    "GSE154826": extract_gse154826,
    "GSE253013": extract_gse253013_from_assembled,
}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, default=DATA)
    ap.add_argument("--outdir", type=Path, default=EXTRACTED)
    ap.add_argument("--cohorts", default=",".join(COHORT_FUNCS))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    for name in [c.strip() for c in args.cohorts.split(",") if c.strip()]:
        print("EXTRACT", name, flush=True)
        COHORT_FUNCS[name](args.data, args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
