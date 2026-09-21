#!/usr/bin/env python3
"""Subsample Q1/Q4 cells and store a shared gene panel for Augur.

Raw counts and library size are saved. CLDN4 stays in the matrix so the
classifier can drop it. GSE205335 is added when the RDS extract exists.
"""
from __future__ import annotations

import gzip
import io
import tarfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import io as spio
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
GEO = Path("/tmp/geo_c4")
OUT = Path("/tmp/geo_c4/augur_extract")
OUT.mkdir(parents=True, exist_ok=True)

PER_UNIT = 25
SEED = 1
MIN_CELLS = 8

MARKER_MAL = ("EPCAM", "KRT8", "KRT18", "KRT19")
MARKER_T = ("CD3D", "CD3E", "CD3G", "CD4", "CD8A")
MARKER_NK = ("NKG7", "GNLY", "KLRD1", "NCAM1", "NCR1")
MARKER_MYE = ("LYZ", "CD68", "CD14", "FCGR3A", "CD163", "C1QA", "MARCO", "AIF1")
MARKER_B = ("MS4A1", "CD79A", "CD19")


def panel_genes() -> list[str]:
    genes = []
    seen = set()
    for line in (DATA / "gene_panel.txt").read_text().splitlines():
        line = line.split("#", 1)[0].strip().upper()
        if line and line not in seen:
            seen.add(line)
            genes.append(line)
    return genes


def classify_arrays(present: dict[str, np.ndarray], n: int) -> np.ndarray:
    def any_of(genes: tuple[str, ...]) -> np.ndarray:
        acc = np.zeros(n, dtype=bool)
        for g in genes:
            if g in present:
                acc |= present[g]
        return acc

    ptprc = present.get("PTPRC", np.zeros(n, dtype=bool))
    mal = any_of(MARKER_MAL) & ~ptprc
    t = ~mal & any_of(MARKER_T)
    nk = ~mal & ~t & any_of(MARKER_NK)
    mye = ~mal & ~t & ~nk & any_of(MARKER_MYE)
    b = ~mal & ~t & ~nk & ~mye & any_of(MARKER_B)
    labels = np.array(["other"] * n, dtype=object)
    labels[b] = "B"
    labels[mye] = "myeloid"
    labels[nk] = "NK"
    labels[t] = "T"
    labels[mal] = "malignant"
    return labels


def subsample(labels: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    keep = []
    for lab in ("T", "NK", "myeloid", "malignant", "B"):
        idx = np.flatnonzero(labels == lab)
        if len(idx) < MIN_CELLS:
            continue
        take = idx if len(idx) <= PER_UNIT else rng.choice(idx, PER_UNIT, replace=False)
        keep.append(np.sort(take))
    if not keep:
        return np.array([], dtype=int)
    return np.concatenate(keep)


def q_units(tab: pd.DataFrame, dataset: str) -> pd.DataFrame:
    sub = tab[(tab.dataset == dataset) & (tab.cldn4_quartile.isin(["Q1", "Q4"]))]
    return sub


def store_unit(bucket: list, dataset: str, unit: str, quartile: str, labels, keep, counts, lib, genes):
    # counts: (n_genes, n_keep)
    for j, cell_i in enumerate(keep):
        bucket.append(
            {
                "dataset": dataset,
                "unit_id": unit,
                "quartile": quartile,
                "cell_type": labels[cell_i],
                "lib_size": float(lib[j]),
                "expr": counts[:, j].astype(np.float32),
            }
        )


def extract_123(tab: pd.DataFrame, genes: list[str], bucket: list) -> None:
    sub = q_units(tab, "GSE123902")
    marker = pd.read_csv(GEO / "GSE123902_marker_units.tsv", sep="\t")
    marker = marker[marker.tissue.isin(["PRIMARY", "METASTASIS"])]
    marker = marker.sort_values(["patient", "tissue"]).drop_duplicates("patient")
    file_of = dict(zip(marker.patient, marker.file))
    want = {file_of[u]: u for u in sub.unit_id}
    q_of = dict(zip(sub.unit_id, sub.cldn4_quartile))
    gene_set = set(genes)
    with tarfile.open(GEO / "GSE123902_RAW.tar") as tar:
        for member in tar.getmembers():
            if member.name not in want:
                continue
            unit = want[member.name]
            raw = gzip.GzipFile(fileobj=tar.extractfile(member)).read()
            df = pd.read_csv(io.BytesIO(raw), index_col=0)
            # cells x genes
            df.columns = df.columns.astype(str).str.upper()
            # sum duplicate gene columns
            if df.columns.duplicated().any():
                df = df.T.groupby(level=0).sum().T
            n = df.shape[0]
            present = {}
            for g in set(MARKER_MAL + MARKER_T + MARKER_NK + MARKER_MYE + MARKER_B + ("PTPRC",)):
                if g in df.columns:
                    present[g] = df[g].to_numpy() > 0
            labels = classify_arrays(present, n)
            rng = np.random.default_rng(SEED)
            keep = subsample(labels, rng)
            if len(keep) == 0:
                continue
            use = [g for g in genes if g in df.columns]
            mat = df.iloc[keep][use].to_numpy(dtype=np.float64)
            lib = df.iloc[keep].to_numpy(dtype=np.float64).sum(axis=1)
            full = np.zeros((len(genes), len(keep)), dtype=np.float32)
            for j, g in enumerate(use):
                full[genes.index(g)] = mat[:, j]
            store_unit(bucket, "GSE123902", unit, q_of[unit], labels, keep, full, lib, genes)
            print(f"  123902 {unit} kept={len(keep)}", flush=True)


def extract_189(tab: pd.DataFrame, genes: list[str], bucket: list) -> None:
    sub = q_units(tab, "GSE189357")
    q_of = dict(zip(sub.unit_id, sub.cldn4_quartile))
    need = set(genes) | set(MARKER_MAL + MARKER_T + MARKER_NK + MARKER_MYE + MARKER_B + ("PTPRC",))
    with tarfile.open(GEO / "GSE189357_RAW.tar") as tar:
        names = tar.getnames()
        for unit, quartile in q_of.items():
            hits = [n for n in names if n.endswith(f"_{unit}_matrix.mtx.gz")]
            prefix = hits[0][: -len("_matrix.mtx.gz")]

            def member(suffix: str):
                return gzip.GzipFile(fileobj=tar.extractfile(f"{prefix}{suffix}"))

            mtx = spio.mmread(member("_matrix.mtx.gz")).tocsc()
            features = pd.read_csv(member("_features.tsv.gz"), sep="\t", header=None)
            symbols = features.iloc[:, 1 if features.shape[1] > 1 else 0].astype(str).str.upper()
            groups: dict[str, list[int]] = defaultdict(list)
            for i, g in enumerate(symbols.tolist()):
                if g in need:
                    groups[g].append(i)
            n = mtx.shape[1]
            lib = np.asarray(mtx.sum(axis=0)).ravel()
            present = {}
            expr_rows = {}
            for g, idxs in groups.items():
                if len(idxs) == 1:
                    row = np.asarray(mtx[idxs[0]].todense()).ravel()
                else:
                    row = np.asarray(mtx[idxs].sum(axis=0)).ravel()
                expr_rows[g] = row
                if g in set(MARKER_MAL + MARKER_T + MARKER_NK + MARKER_MYE + MARKER_B + ("PTPRC",)):
                    present[g] = row > 0
            labels = classify_arrays(present, n)
            rng = np.random.default_rng(SEED)
            keep = subsample(labels, rng)
            if len(keep) == 0:
                continue
            full = np.zeros((len(genes), len(keep)), dtype=np.float32)
            for i, g in enumerate(genes):
                if g in expr_rows:
                    full[i] = expr_rows[g][keep]
            store_unit(bucket, "GSE189357", unit, quartile, labels, keep, full, lib[keep], genes)
            print(f"  189357 {unit} kept={len(keep)}", flush=True)


def extract_131(tab: pd.DataFrame, genes: list[str], bucket: list) -> None:
    sub = q_units(tab, "GSE131907")
    q_of = dict(zip(sub.unit_id, sub.cldn4_quartile))
    ann_path = GEO / "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
    chosen = []  # (col is filled after header), store index -> meta
    # cell index string -> assignment
    wanted_index = {}
    with gzip.open(ann_path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        by_sample = defaultdict(list)
        for line in f:
            p = line.rstrip("\n").split("\t")
            sample = p[idx["Sample"]]
            if sample not in q_of:
                continue
            ctype = p[idx["Cell_type"]]
            subtype = p[idx["Cell_subtype"]]
            if subtype == "Malignant cells":
                lab = "malignant"
            elif ctype == "T lymphocytes":
                lab = "T"
            elif ctype == "NK cells":
                lab = "NK"
            elif ctype == "Myeloid cells":
                lab = "myeloid"
            elif ctype == "B lymphocytes":
                lab = "B"
            else:
                continue
            by_sample[sample].append((p[idx["Index"]], lab))
    rng = np.random.default_rng(SEED)
    keep_index = {}
    for sample, cells in by_sample.items():
        labels_map = defaultdict(list)
        for index, lab in cells:
            labels_map[lab].append(index)
        for lab, ids in labels_map.items():
            if len(ids) < MIN_CELLS:
                continue
            take = ids if len(ids) <= PER_UNIT else list(rng.choice(ids, PER_UNIT, replace=False))
            for index in take:
                keep_index[index] = (sample, lab, q_of[sample])
    print(f"  131907 cells to pull {len(keep_index)}", flush=True)
    mat_path = GEO / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    gene_to_row = {g: i for i, g in enumerate(genes)}
    with gzip.open(mat_path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        barcodes = header[1:]
        col_of = {b: i for i, b in enumerate(barcodes)}
        ordered = [(col_of[b], b) for b in keep_index if b in col_of]
        ordered.sort()
        cols = [c for c, _ in ordered]
        ids = [b for _, b in ordered]
        colpos = {c: j for j, c in enumerate(cols)}
        expr = np.zeros((len(genes), len(cols)), dtype=np.float32)
        lib = np.zeros(len(cols), dtype=np.float64)
        # Fast column picker: split is acceptable (~40s for the whole file).
        target = set(cols)
        n_genes_seen = 0
        for line in f:
            tab_i = line.find("\t")
            gene = line[:tab_i].upper()
            n_genes_seen += 1
            parts = line.rstrip("\n").split("\t")
            # library size and panel
            for c in cols:
                val = parts[c + 1]
                if val == "0":
                    continue
                v = float(val)
                lib[colpos[c]] += v
                gi = gene_to_row.get(gene)
                if gi is not None:
                    expr[gi, colpos[c]] = v
            if n_genes_seen % 5000 == 0:
                print(f"    streamed {n_genes_seen} genes", flush=True)
    for j, index in enumerate(ids):
        sample, lab, quartile = keep_index[index]
        bucket.append(
            {
                "dataset": "GSE131907",
                "unit_id": sample,
                "quartile": quartile,
                "cell_type": lab,
                "lib_size": float(lib[j]),
                "expr": expr[:, j],
            }
        )
    print(f"  131907 stored {len(ids)}", flush=True)


def save_bucket(bucket: list, genes: list[str], path: Path) -> None:
    meta = pd.DataFrame(
        [
            {k: r[k] for k in ("dataset", "unit_id", "quartile", "cell_type", "lib_size")}
            for r in bucket
        ]
    )
    expr = np.vstack([r["expr"] for r in bucket])
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, expr=expr, genes=np.array(genes))
    meta.to_csv(path.with_suffix(".tsv"), sep="\t", index=False)
    print(f"wrote {path} cells={len(meta)} genes={len(genes)}", flush=True)


def main() -> None:
    tab = pd.read_csv(ROOT / "results" / "tables" / "composition_counts.tsv", sep="\t")
    genes = panel_genes()
    bucket: list = []
    print("GSE123902", flush=True)
    extract_123(tab, genes, bucket)
    print("GSE189357", flush=True)
    extract_189(tab, genes, bucket)
    print("GSE131907", flush=True)
    extract_131(tab, genes, bucket)
    expr_path = OUT / "gse205335_expr.tsv.gz"
    meta_path = OUT / "gse205335_panel.tsv"
    if expr_path.exists() and meta_path.exists():
        meta205 = pd.read_csv(meta_path, sep="\t")
        expr205 = pd.read_csv(expr_path, sep="\t")
        expr205.columns = [c.upper() for c in expr205.columns]
        if list(expr205.columns) != genes:
            missing = [g for g in genes if g not in expr205.columns]
            extra = [g for g in expr205.columns if g not in genes]
            raise SystemExit(f"GSE205335 panel mismatch missing={missing[:8]} extra={extra[:8]}")
        arr = expr205[genes].to_numpy(dtype=np.float32)
        for i, rec in meta205.iterrows():
            bucket.append(
                {
                    "dataset": "GSE205335",
                    "unit_id": rec.unit_id,
                    "quartile": rec.quartile,
                    "cell_type": rec.cell_type,
                    "lib_size": float(rec.lib_size),
                    "expr": arr[i],
                }
            )
        print(f"added GSE205335 cells={len(meta205)}", flush=True)
    else:
        print("GSE205335 panel not ready; Augur will be marked partial", flush=True)
    save_bucket(bucket, genes, OUT / "augur_panel.npz")


if __name__ == "__main__":
    main()
