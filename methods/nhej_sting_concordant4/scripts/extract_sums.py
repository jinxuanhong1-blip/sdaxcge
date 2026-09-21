#!/usr/bin/env python3
"""Malignant UMI sums for concordant-4 NHEJ / STING / IFN / MHC-I.

Datasets: GSE123902, GSE131907, GSE189357.
GSE205335 is exported by export_gse205335.R (RDS).
Honest unit list is the locked concordant-4 patient/donor/sample table.
"""
from __future__ import annotations

import gzip
import json
import tarfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import io as spio
from scipy import sparse

ROOT = Path("/workspace/methods/nhej_sting_concordant4")
GEO = Path("/tmp/geo_nhej")
OUT = ROOT / "results" / "tables"
GENE_SETS = ROOT / "data" / "gene_sets.json"
LOCKED = ROOT / "data" / "locked_patient_units_pr539.tsv"

MARKERS = ["EPCAM", "KRT8", "KRT18", "KRT19", "PTPRC", "CD3D", "CD3E", "CD8A", "NKG7", "GNLY", "KLRD1", "CLDN4"]


def load_sets():
    gs = json.loads(GENE_SETS.read_text())
    aliases = {k.upper(): v.upper() for k, v in gs["aliases"].items()}
    canonical = set()
    for key in ("nhej_kegg", "nhej_core", "sting_core", "ifn_hallmark", "mhc_i", "prolif"):
        canonical.update(g.upper() for g in gs[key])
    canonical.add("CLDN4")
    # symbols to pull from a matrix: canonical + aliases
    raw_needed = set(canonical) | set(aliases)
    return gs, aliases, canonical, raw_needed


def canon_symbol(symbol: str, aliases: dict[str, str]) -> str:
    s = symbol.upper()
    return aliases.get(s, s)


def locked_units(dataset: str) -> list[str]:
    df = pd.read_csv(LOCKED, sep="\t")
    return df.loc[df["dataset"] == dataset, "unit_id"].astype(str).tolist()


def marker_masks(symbol_to_vec: dict[str, np.ndarray], n: int):
    def g(name):
        v = symbol_to_vec.get(name)
        if v is None:
            return np.zeros(n, dtype=np.float64)
        return v

    mal = (
        ((g("EPCAM") > 0) | (g("KRT8") > 0) | (g("KRT18") > 0) | (g("KRT19") > 0))
        & (g("PTPRC") == 0)
    )
    tnk = (
        (g("CD3D") > 0)
        | (g("CD3E") > 0)
        | (g("CD8A") > 0)
        | (g("NKG7") > 0)
        | (g("GNLY") > 0)
        | (g("KLRD1") > 0)
    ) & ~mal
    return mal, tnk


def rows_for_unit(dataset, unit, gene_umi: dict[str, float], n_cells, n_mal, n_tnk, cldn4, lib):
    cldn4 = np.asarray(cldn4, dtype=np.float64)
    meta = {
        "dataset": dataset,
        "unit_id": unit,
        "n_cells": int(n_cells),
        "n_malignant": int(n_mal),
        "n_tnk": int(n_tnk),
        "mal_CLDN4_pct": float(100.0 * np.mean(cldn4 > 0)) if n_mal else np.nan,
        "mal_CLDN4_mean": float(np.mean(np.log1p(cldn4))) if n_mal else np.nan,
        "lib_umi": float(lib),
    }
    grows = [
        {"dataset": dataset, "unit_id": unit, "gene": gene, "umi": float(umi)}
        for gene, umi in sorted(gene_umi.items())
    ]
    return meta, grows


def collapse_sums(symbol_sums: dict[str, float], aliases: dict[str, str], keep: set[str]) -> dict[str, float]:
    out = defaultdict(float)
    for sym, val in symbol_sums.items():
        c = canon_symbol(sym, aliases)
        if c in keep:
            out[c] += float(val)
    return dict(out)


def extract_gse123902(aliases, keep_genes):
    print("GSE123902", flush=True)
    marker = pd.read_csv(ROOT / "data" / "gse123902_files.tsv", sep="\t")
    tar_path = GEO / "GSE123902_RAW.tar"
    tf = tarfile.open(tar_path)
    metas, grows = [], []
    for rec in marker.itertuples(index=False):
        print(f"  {rec.unit_id} {rec.file}", flush=True)
        member = tf.extractfile(rec.file)
        raw = pd.read_csv(gzip.GzipFile(fileobj=member), index_col=0)
        raw.columns = [str(c).upper() for c in raw.columns]
        # cells x genes
        mat = raw.to_numpy(dtype=np.float64)
        genes = list(raw.columns)
        # collapse duplicate symbols
        sym_vec = {}
        for i, g in enumerate(genes):
            c = canon_symbol(g, aliases)
            col = mat[:, i]
            if c in sym_vec:
                sym_vec[c] = sym_vec[c] + col
            else:
                sym_vec[c] = col
        mal, tnk = marker_masks(sym_vec, mat.shape[0])
        if mal.sum() == 0:
            raise SystemExit(f"no malignant cells in {rec.unit_id}")
        mal_mat = mat[mal]
        lib = float(mal_mat.sum())
        # gene sums on malignant cells
        sums = {}
        for i, g in enumerate(genes):
            c = canon_symbol(g, aliases)
            if c in keep_genes:
                sums[c] = sums.get(c, 0.0) + float(mal_mat[:, i].sum())
        cldn4 = sym_vec.get("CLDN4", np.zeros(mat.shape[0]))[mal]
        meta, grow = rows_for_unit(
            "GSE123902", rec.unit_id, sums, mat.shape[0], int(mal.sum()), int(tnk.sum()), cldn4, lib
        )
        metas.append(meta)
        grows.extend(grow)
        print(
            f"    n_mal={meta['n_malignant']} pct={meta['mal_CLDN4_pct']:.4f} "
            f"mean={meta['mal_CLDN4_mean']:.4f}",
            flush=True,
        )
    return metas, grows


def extract_gse189357(aliases, keep_genes):
    print("GSE189357", flush=True)
    units = locked_units("GSE189357")
    tf = tarfile.open(GEO / "GSE189357_RAW.tar")
    names = tf.getnames()
    metas, grows = [], []
    dest = GEO / "gse189357"
    dest.mkdir(exist_ok=True)
    for unit in units:
        need = [n for n in names if f"_{unit}_" in Path(n).name]
        print(f"  {unit} files={len(need)}", flush=True)
        paths = {}
        for n in need:
            out = dest / Path(n).name
            if not out.exists():
                src = tf.extractfile(n)
                out.write_bytes(src.read())
            if n.endswith("_matrix.mtx.gz"):
                paths["mtx"] = out
            elif n.endswith("_features.tsv.gz"):
                paths["feat"] = out
            elif n.endswith("_barcodes.tsv.gz"):
                paths["bc"] = out
        feat = pd.read_csv(paths["feat"], sep="\t", header=None, compression="gzip")
        symbols = feat.iloc[:, 1].astype(str).str.upper().tolist() if feat.shape[1] > 1 else feat.iloc[:, 0].astype(str).str.upper().tolist()
        mat = spio.mmread(paths["mtx"]).tocsr()
        if mat.shape[0] != len(symbols):
            raise SystemExit(f"{unit} features {len(symbols)} != mtx rows {mat.shape[0]}")
        # collapse symbols into a dense marker block and sparse sums
        canon = [canon_symbol(s, aliases) for s in symbols]
        n = mat.shape[1]
        sym_vec = {}
        for gene in MARKERS:
            rows = [i for i, c in enumerate(canon) if c == gene]
            if not rows:
                sym_vec[gene] = np.zeros(n, dtype=np.float64)
            else:
                sym_vec[gene] = np.asarray(mat[rows].sum(axis=0)).ravel()
        mal, tnk = marker_masks(sym_vec, n)
        mal_idx = np.flatnonzero(mal)
        if mal_idx.size == 0:
            raise SystemExit(f"no malignant cells in {unit}")
        sub = mat[:, mal_idx]
        lib = float(sub.sum())
        # sum rows by canonical gene
        acc = defaultdict(float)
        # only rows whose canonical id is needed
        wanted_rows = {}
        for i, c in enumerate(canon):
            if c in keep_genes:
                wanted_rows.setdefault(c, []).append(i)
        sums = {}
        for c, rows in wanted_rows.items():
            sums[c] = float(sub[rows].sum())
        cldn4 = sym_vec["CLDN4"][mal_idx]
        meta, grow = rows_for_unit(
            "GSE189357", unit, sums, n, int(mal.sum()), int(tnk.sum()), cldn4, lib
        )
        metas.append(meta)
        grows.extend(grow)
        print(
            f"    n_mal={meta['n_malignant']} pct={meta['mal_CLDN4_pct']:.4f} "
            f"mean={meta['mal_CLDN4_mean']:.4f} genes={len(sums)}",
            flush=True,
        )
        del mat, sub
    return metas, grows


def extract_gse131907(aliases, keep_genes):
    print("GSE131907", flush=True)
    units = locked_units("GSE131907")
    unit_set = set(units)
    ann = pd.read_csv(GEO / "GSE131907_Lung_Cancer_cell_annotation.txt.gz", sep="\t")
    ann["author_malignant"] = ann["Cell_subtype"] == "Malignant cells"
    ann["author_tnk"] = ann["Cell_type"].isin(["T lymphocytes", "NK cells"])
    ann = ann[ann["Sample"].isin(unit_set)]
    # header barcodes are the Index column
    mal = ann[ann["author_malignant"]]
    mal_by_sample = {s: set(g["Index"]) for s, g in mal.groupby("Sample")}
    n_cells = ann.groupby("Sample").size().to_dict()
    n_tnk = ann[ann["author_tnk"]].groupby("Sample").size().to_dict()
    n_mal = mal.groupby("Sample").size().to_dict()
    print(
        f"  locked samples={len(units)} malignant cells={len(mal)}",
        flush=True,
    )

    path = GEO / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    # Discover header and which columns we keep.
    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
    barcodes = header[1:]
    col_sample = []
    usecols = ["Index"]
    for b in barcodes:
        # sample id is the suffix after the first underscore of Index: BARCODE_SAMPLE
        # Index itself is BARCODE_SAMPLE and matches annotation Index.
        if b in mal_by_sample.get(b.split("_", 1)[-1] if False else "", ()):
            pass
    # Build membership from the annotation Index set, not a parse of the barcode.
    mal_index = set(mal["Index"])
    sample_of = dict(zip(mal["Index"], mal["Sample"]))
    keep_bc = [b for b in barcodes if b in mal_index]
    missing = len(mal_index) - len(keep_bc)
    print(f"  matrix malignant columns={len(keep_bc)} missing={missing}", flush=True)
    if len(keep_bc) < 1000:
        raise SystemExit("GSE131907 malignant barcode match failed")

    # column groups in the frame we will read (Index + keep_bc)
    bc_to_sample = {b: sample_of[b] for b in keep_bc}
    lib = {s: 0.0 for s in units}
    cldn4_pos = {s: 0 for s in units}
    cldn4_log = {s: 0.0 for s in units}
    gene_sums = {s: defaultdict(float) for s in units}
    genes_seen = set()

    reader = pd.read_csv(
        path,
        sep="\t",
        usecols=["Index"] + keep_bc,
        chunksize=400,
        engine="c",
        low_memory=False,
    )
    n_genes = 0
    for chunk in reader:
        genes_idx = chunk["Index"].astype(str)
        chunk = chunk.drop(columns=["Index"])
        n_genes += chunk.shape[0]
        chunk_np = np.nan_to_num(chunk.to_numpy(dtype=np.float64), nan=0.0, posinf=0.0, neginf=0.0)
        # library size
        # group columns
        if n_genes == chunk.shape[0]:
            # cache column index lists once
            cols = list(chunk.columns)
            groups = {s: [] for s in units}
            for i, b in enumerate(cols):
                groups[bc_to_sample[b]].append(i)
            group_idx = {s: np.asarray(ix, dtype=np.int32) for s, ix in groups.items()}
        for s, ix in group_idx.items():
            if ix.size == 0:
                continue
            block = chunk_np[:, ix]
            lib[s] += float(block.sum())
        symbols = [canon_symbol(g, aliases) for g in genes_idx]
        for r, sym in enumerate(symbols):
            genes_seen.add(sym)
            if sym not in keep_genes and sym != "CLDN4":
                continue
            for s, ix in group_idx.items():
                if ix.size == 0:
                    continue
                vec = chunk_np[r, ix]
                if sym in keep_genes:
                    gene_sums[s][sym] += float(vec.sum())
                if sym == "CLDN4":
                    cldn4_pos[s] += int(np.count_nonzero(vec > 0))
                    cldn4_log[s] += float(np.log1p(vec).sum())
        if n_genes % 4000 < 400:
            print(f"  streamed {n_genes} genes", flush=True)
    print(f"  done genes={n_genes}", flush=True)

    metas, grows = [], []
    for s in units:
        nmal = int(n_mal.get(s, 0))
        meta = {
            "dataset": "GSE131907",
            "unit_id": s,
            "n_cells": int(n_cells.get(s, 0)),
            "n_malignant": nmal,
            "n_tnk": int(n_tnk.get(s, 0)),
            "mal_CLDN4_pct": float(100.0 * cldn4_pos[s] / nmal) if nmal else np.nan,
            "mal_CLDN4_mean": float(cldn4_log[s] / nmal) if nmal else np.nan,
            "lib_umi": float(lib[s]),
        }
        # write every keep gene that existed in the matrix, including zeros
        present = [g for g in sorted(keep_genes) if g in genes_seen or g in gene_sums[s]]
        # genes_seen is canonical of ALL genes, so membership is correct
        sums = {g: float(gene_sums[s].get(g, 0.0)) for g in sorted(keep_genes) if g in genes_seen}
        metas.append(meta)
        grows.extend(
            {"dataset": "GSE131907", "unit_id": s, "gene": g, "umi": u} for g, u in sums.items()
        )
        print(
            f"  {s} n_mal={nmal} pct={meta['mal_CLDN4_pct']:.4f} mean={meta['mal_CLDN4_mean']:.4f}",
            flush=True,
        )
    return metas, grows


def write_wanted(raw_needed: set[str]):
    path = ROOT / "data" / "wanted_genes.txt"
    path.write_text("\n".join(sorted(raw_needed)) + "\n")
    return path


def main():
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="123902,189357,131907", help="comma list of GSE suffixes")
    args = ap.parse_args()
    only = {x.strip() for x in args.only.split(",") if x.strip()}
    gs, aliases, canonical, raw_needed = load_sets()
    write_wanted(raw_needed | set(MARKERS))
    # file map for GSE123902: tumor rows only, one per donor, same rule as PR 539
    marker = pd.read_csv(ROOT / "data" / "GSE123902_marker_units.tsv", sep="\t")
    tumor = marker[marker["tissue"].isin(["PRIMARY", "METASTASIS"])].copy()
    tumor = tumor.sort_values(["patient", "tissue"])
    tumor = tumor.drop_duplicates("patient", keep="first")
    tumor = tumor[(tumor["n_malignant"] >= 20) & (tumor["n_tnk"] >= 20)]
    locked = set(locked_units("GSE123902"))
    if set(tumor["patient"]) != locked:
        raise SystemExit(f"GSE123902 unit mismatch {set(tumor['patient']) ^ locked}")
    file_map = tumor[["patient", "file", "tissue"]].rename(columns={"patient": "unit_id"})
    file_map.to_csv(ROOT / "data" / "gse123902_files.tsv", sep="\t", index=False)

    metas, grows = [], []
    if "123902" in only:
        m, g = extract_gse123902(aliases, canonical)
        metas.extend(m)
        grows.extend(g)
    if "189357" in only:
        m, g = extract_gse189357(aliases, canonical)
        metas.extend(m)
        grows.extend(g)
    if "131907" in only:
        m, g = extract_gse131907(aliases, canonical)
        metas.extend(m)
        grows.extend(g)

    OUT.mkdir(parents=True, exist_ok=True)
    tag = "python" if only == {"123902", "189357", "131907"} else "py_" + "_".join(sorted(only))
    pd.DataFrame(metas).to_csv(OUT / f"unit_meta_{tag}.tsv", sep="\t", index=False)
    pd.DataFrame(grows).to_csv(OUT / f"gene_umi_{tag}.tsv", sep="\t", index=False)
    print(f"wrote {tag} units={len(metas)} gene_rows={len(grows)}", flush=True)


if __name__ == "__main__":
    main()
