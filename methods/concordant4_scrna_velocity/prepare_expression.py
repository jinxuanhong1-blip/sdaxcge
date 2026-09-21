#!/usr/bin/env python3
"""Build a capped epithelial count matrix for concordant-4 trajectory sweeps.

Tumor units are the locked 65 (GSE123902 13, GSE131907 21, GSE205335 22,
GSE189357 9). Uninvolved-lung / normal epithelial cells are kept only as a
root pool and are not inferential units.

Outputs /tmp/c4data/epi_counts.npz (not committed; matrices are large).
"""

from __future__ import annotations

import csv
import gzip
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import io as spio
from scipy import sparse

DATA = Path("/tmp/c4data")
OUT = DATA / "epi_counts.npz"
RNG = np.random.default_rng(4)
CAP_TUMOR = 160
CAP_ROOT = 80

TUMOR_131907 = {
    "BRONCHO_11", "BRONCHO_58", "EBUS_06", "EBUS_10", "EBUS_12", "EBUS_13",
    "EBUS_15", "EBUS_19", "EBUS_28", "EBUS_49", "EBUS_51", "NS_02", "NS_03",
    "NS_04", "NS_06", "NS_07", "NS_12", "NS_13", "NS_16", "NS_17", "NS_19",
}
TUMOR_123902 = {
    "LX255B", "LX653", "LX661", "LX666", "LX675", "LX676", "LX679", "LX680",
    "LX681", "LX682", "LX684", "LX699", "LX701",
}
TUMOR_189 = {f"TD{i}" for i in range(1, 10)}


def cap_indices(labels: np.ndarray, cap: int, rng: np.random.Generator) -> np.ndarray:
    keep = []
    for lab in pd.unique(labels):
        idx = np.flatnonzero(labels == lab)
        if len(idx) > cap:
            idx = rng.choice(idx, size=cap, replace=False)
        keep.append(idx)
    return np.sort(np.concatenate(keep)) if keep else np.array([], dtype=int)


def epi_mask_from_counts(genes: list[str], get_count) -> np.ndarray:
    """get_count(gene) -> 1d array over cells. Epithelial = PTPRC==0 and an epithelial gene >0."""
    upper = {g.upper(): g for g in genes}
    def vec(name):
        g = upper.get(name)
        if g is None:
            return None
        return get_count(g)
    ptprc = vec("PTPRC")
    if ptprc is None:
        raise RuntimeError("PTPRC missing")
    epi = np.zeros(len(ptprc), dtype=bool)
    for name in ("EPCAM", "KRT8", "KRT18", "KRT19"):
        v = vec(name)
        if v is not None:
            epi |= v > 0
    return epi & (ptprc == 0)


def store_block(blocks: list, genes: list[str], counts_cg: np.ndarray, dataset: str, unit: np.ndarray, pool: np.ndarray, tissue: np.ndarray) -> None:
    """counts_cg: genes x cells uint32."""
    blocks.append({
        "genes": [g.upper() for g in genes],
        "X": counts_cg.T.astype(np.float32, copy=False),  # cells x genes
        "dataset": np.full(counts_cg.shape[1], dataset),
        "unit": unit.astype(str),
        "pool": pool.astype(str),
        "tissue": tissue.astype(str),
    })


def load_123902(blocks: list) -> None:
    folder = DATA / "gse123"
    # filename pattern GSM*_MSK_<unit>_<TISSUE>_dense.csv.gz
    for path in sorted(folder.glob("*.csv.gz")):
        m = re.search(r"MSK_(.+?)_(PRIMARY_TUMOUR|METASTASIS|NORMAL)_dense", path.name)
        if not m:
            print("skip", path.name)
            continue
        unit, tissue = m.group(1), m.group(2)
        if tissue == "NORMAL":
            pool_name = "root"
        elif unit in TUMOR_123902:
            pool_name = "tumor"
        else:
            continue
        df = pd.read_csv(path, index_col=0)
        genes = [str(c) for c in df.columns]
        mat = df.to_numpy(dtype=np.float32).T  # genes x cells
        # SEQC dense matrices are counts; guard against a log matrix
        if np.nanmax(mat) < 30:
            counts = np.rint(np.expm1(mat)).astype(np.float32)
        else:
            counts = mat
        epi = epi_mask_from_counts(genes, lambda g: counts[genes.index(g)])
        sel = np.flatnonzero(epi)
        if sel.size == 0:
            print(path.name, "no epi")
            continue
        cap = CAP_ROOT if pool_name == "root" else CAP_TUMOR
        if sel.size > cap:
            sel = RNG.choice(sel, size=cap, replace=False)
            sel.sort()
        store_block(
            blocks, genes, counts[:, sel], "GSE123902",
            np.full(sel.size, unit), np.full(sel.size, pool_name), np.full(sel.size, tissue),
        )
        print(f"GSE123902 {unit} {tissue} {pool_name} {sel.size}", flush=True)


def load_189357(blocks: list) -> None:
    folder = DATA / "gse189"
    for td in [f"TD{i}" for i in range(1, 10)]:
        feat = next(folder.glob(f"*_{td}_features.tsv.gz"))
        bar = next(folder.glob(f"*_{td}_barcodes.tsv.gz"))
        mtx = next(folder.glob(f"*_{td}_matrix.mtx.gz"))
        genes = [ln.split("\t")[1].upper() for ln in gzip.open(feat, "rt")]
        _barcodes = [ln.strip() for ln in gzip.open(bar, "rt")]
        with gzip.open(mtx, "rb") as fh:
            X = spio.mmread(fh).tocsr()  # genes x cells
        if X.shape[0] != len(genes):
            X = X.T.tocsr()
        lookup = {g: i for i, g in enumerate(genes)}
        def get_count(g):
            i = lookup.get(g)
            if i is None:
                return None
            return np.asarray(X[i].todense()).ravel()
        epi = epi_mask_from_counts(genes, get_count)
        sel = np.flatnonzero(epi)
        if sel.size > CAP_TUMOR:
            sel = RNG.choice(sel, size=CAP_TUMOR, replace=False)
            sel.sort()
        sub = X[:, sel].toarray().astype(np.float32)
        store_block(
            blocks, genes, sub, "GSE189357",
            np.full(sel.size, td), np.full(sel.size, "tumor"), np.full(sel.size, "TUMOR"),
        )
        print(f"GSE189357 {td} epi_kept {sel.size} / {epi.sum()}", flush=True)
        del X


def load_131907(blocks: list) -> None:
    ann = list(csv.DictReader(gzip.open(DATA / "GSE131907_ann.txt.gz", "rt"), delimiter="\t"))
    want = {}
    for r in ann:
        if r["Cell_type"] != "Epithelial cells":
            continue
        sample = r["Sample"]
        origin = r["Sample_Origin"]
        if sample in TUMOR_131907:
            want[r["Index"]] = (sample, "tumor", origin)
        elif origin == "nLung":
            want[r["Index"]] = (sample, "root", origin)
    print("GSE131907 epithelial candidates", len(want), flush=True)
    path = DATA / "GSE131907_UMI.txt.gz"
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        name_to_j = {name: j for j, name in enumerate(header)}
        chosen = [n for n in want if n in name_to_j]
        cols = np.array([name_to_j[n] for n in chosen], dtype=np.int32)
        meta = [want[n] for n in chosen]
        # cap per sample before reading genes
        units = np.array([m[0] for m in meta])
        pools = np.array([m[1] for m in meta])
        tissues = np.array([m[2] for m in meta])
        keep_local = []
        for lab in pd.unique(units):
            idx = np.flatnonzero(units == lab)
            cap = CAP_ROOT if pools[idx[0]] == "root" else CAP_TUMOR
            if len(idx) > cap:
                idx = RNG.choice(idx, size=cap, replace=False)
            keep_local.append(idx)
        keep_local = np.sort(np.concatenate(keep_local))
        cols = cols[keep_local]
        units, pools, tissues = units[keep_local], pools[keep_local], tissues[keep_local]
        print("GSE131907 cells kept", len(cols), flush=True)
        colset = {int(c): k for k, c in enumerate(cols)}
        genes = []
        rows = []
        for line in fh:
            # first field gene, then counts. Extract only selected columns.
            parts = line.rstrip("\n").split("\t")
            genes.append(parts[0].upper())
            vals = np.zeros(len(cols), dtype=np.float32)
            for src, dst in colset.items():
                vals[dst] = float(parts[src]) if src < len(parts) and parts[src] else 0.0
            rows.append(vals)
    counts = np.vstack(rows)
    store_block(blocks, genes, counts, "GSE131907", units, pools, tissues)
    print("GSE131907 genes", len(genes), "cells", counts.shape[1], flush=True)


def parse_205335_map() -> dict[str, dict]:
    """orig.ident -> patient/tissue. orig.ident matches the identity table."""
    out = {}
    for p in (DATA / "gsm").glob("GSM*.txt"):
        t = p.read_text(errors="replace")
        def grab(k):
            m = re.search(rf"!Sample_characteristics_ch1 = {k}: (.*)", t)
            return m.group(1).strip() if m else ""
        desc = re.search(r"!Sample_description = (.*)", t)
        plat = grab("platform")
        end = "3" if "3'" in plat else ("5" if "5'" in plat else "")
        if not desc or not end:
            continue
        orig = desc.group(1).strip().replace("_", "-") + f"-{end}P"
        out[orig] = {"patient": grab("patient"), "tissue": grab("tissue"), "gsm": p.stem}
    return out


def load_205335(blocks: list) -> None:
    mtx_path = DATA / "GSE205335_epi.mtx"
    gene_path = DATA / "GSE205335_epi_genes.txt"
    cell_path = DATA / "GSE205335_epi_cells.txt"
    if not mtx_path.exists():
        raise SystemExit("GSE205335 epithelial mtx missing; run subset_gse205335.R first")
    X = spio.mmread(mtx_path).tocsr()
    genes = [ln.strip().upper() for ln in gene_path.read_text().splitlines() if ln.strip()]
    cells = [ln.strip() for ln in cell_path.read_text().splitlines() if ln.strip()]
    if X.shape[0] != len(genes):
        X = X.T.tocsr()
    meta = parse_205335_map()
    ident = list(csv.DictReader(gzip.open(DATA / "GSE205335_id.txt.gz", "rt"), delimiter="\t"))
    by_bc = {r["barcode"]: r for r in ident}
    tumor_patients = set()
    # locked 22 from patient_units
    units_path = Path(__file__).resolve().parent / "data" / "patient_units.tsv"
    for r in csv.DictReader(units_path.open(), delimiter="\t"):
        if r["dataset"] == "GSE205335":
            tumor_patients.add(r["unit_id"])
    unit, pool, tissue = [], [], []
    keep = []
    for j, bc in enumerate(cells):
        rec = by_bc.get(bc)
        if rec is None or rec["lineage.total"] != "Epithelial cells":
            continue
        info = meta.get(rec["orig.ident"])
        if info is None:
            continue
        patient, tis = info["patient"], info["tissue"]
        if tis.lower().startswith("normal"):
            unit.append(patient); pool.append("root"); tissue.append(tis); keep.append(j)
        elif patient in tumor_patients:
            unit.append(patient); pool.append("tumor"); tissue.append(tis); keep.append(j)
    keep = np.array(keep, dtype=int)
    unit = np.array(unit); pool = np.array(pool); tissue = np.array(tissue)
    local = cap_indices_by_pool(unit, pool)
    keep = keep[local]
    unit, pool, tissue = unit[local], pool[local], tissue[local]
    sub = X[:, keep].toarray().astype(np.float32)
    store_block(blocks, genes, sub, "GSE205335", unit, pool, tissue)
    print("GSE205335 kept", keep.size, flush=True)


def cap_indices_by_pool(unit: np.ndarray, pool: np.ndarray) -> np.ndarray:
    keep = []
    for i in range(len(unit)):
        pass
    labels = np.array([f"{p}|{u}" for p, u in zip(pool, unit)])
    for lab in pd.unique(labels):
        idx = np.flatnonzero(labels == lab)
        cap = CAP_ROOT if lab.startswith("root|") else CAP_TUMOR
        if len(idx) > cap:
            idx = RNG.choice(idx, size=cap, replace=False)
        keep.append(idx)
    return np.sort(np.concatenate(keep))


def intersect_and_save(blocks: list) -> None:
    shared = set(blocks[0]["genes"])
    for b in blocks[1:]:
        shared &= set(b["genes"])
    shared = sorted(shared)
    print("shared genes", len(shared), "blocks", len(blocks), flush=True)
    if "CLDN4" not in shared:
        raise SystemExit("CLDN4 not in shared genes")
    Xs, ds, us, ps, ts = [], [], [], [], []
    for b in blocks:
        idx = {g: i for i, g in enumerate(b["genes"])}
        take = [idx[g] for g in shared]
        Xs.append(b["X"][:, take])
        ds.append(b["dataset"]); us.append(b["unit"]); ps.append(b["pool"]); ts.append(b["tissue"])
    X = np.vstack(Xs).astype(np.float32)
    # library-size log1p CP10k
    lib = X.sum(1, keepdims=True)
    lib[lib == 0] = 1
    logx = np.log1p(X / lib * 1e4).astype(np.float32)
    np.savez_compressed(
        OUT,
        X=logx,
        genes=np.array(shared),
        dataset=np.concatenate(ds),
        unit=np.concatenate(us),
        pool=np.concatenate(ps),
        tissue=np.concatenate(ts),
    )
    print("wrote", OUT, "cells", logx.shape[0], "genes", logx.shape[1], flush=True)


def main() -> None:
    blocks: list = []
    load_123902(blocks)
    load_189357(blocks)
    load_131907(blocks)
    load_205335(blocks)
    intersect_and_save(blocks)


if __name__ == "__main__":
    main()
