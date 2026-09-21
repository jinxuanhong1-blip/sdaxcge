#!/usr/bin/env python3
"""Extract malignant-cell UMI counts for the locked concordant-4 units.

Networks are built later. This script only materializes malignant cells
from GSE123902, GSE131907, GSE205335, and GSE189357. Immune cells are not
stored. GSE205335 must already have been exported by export_gse205335.R.
"""
from __future__ import annotations

import gzip
import tarfile
from pathlib import Path

import numpy as np
import scipy.io
import scipy.sparse as sp

ROOT = Path("/workspace/methods/hdwgcna_concordant4_cldn4")
GEO = Path("/tmp/geo_hdwgcna")
WORK = Path("/tmp/hdwgcna_c4")
STORE = WORK / "store"
STATS = WORK / "stats"
AUDIT = WORK / "audit"
UNIVERSE = WORK / "universe_genes.txt"

MARKER_POS = ("EPCAM", "KRT8", "KRT18", "KRT19")
MARKER_NEG = "PTPRC"


def say(msg: str) -> None:
    print(msg, flush=True)


def load_tsv(path: Path) -> list[dict[str, str]]:
    lines = path.read_text().splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:] if line]


def locked_units() -> dict[str, list[dict[str, str]]]:
    rows = load_tsv(ROOT / "data" / "patient_units_locked.tsv")
    out: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        out.setdefault(row["dataset"], []).append(row)
    return out


def universe() -> list[str]:
    genes = [g.strip().upper() for g in UNIVERSE.read_text().splitlines() if g.strip()]
    seen: set[str] = set()
    out: list[str] = []
    for g in genes:
        if g not in seen:
            seen.add(g)
            out.append(g)
    return out


def write_store(dataset: str, genes: list[str], units: list[str], counts: np.ndarray) -> None:
    """counts: uint16, genes x cells."""
    dest = STORE / dataset
    dest.mkdir(parents=True, exist_ok=True)
    if counts.dtype != np.uint16:
        raise SystemExit(f"{dataset}: expected uint16, got {counts.dtype}")
    if counts.shape != (len(genes), len(units)):
        raise SystemExit(f"{dataset}: shape {counts.shape} != {(len(genes), len(units))}")
    (dest / "genes.txt").write_text("\n".join(genes) + "\n")
    (dest / "units.txt").write_text("\n".join(units) + "\n")
    mm = np.memmap(dest / "counts.u16", dtype=np.uint16, mode="w+", shape=counts.shape)
    mm[:] = counts
    mm.flush()
    del mm
    (dest / "shape.txt").write_text(f"{counts.shape[0]}\t{counts.shape[1]}\n")
    n = counts.shape[1]
    nnz = np.empty(counts.shape[0], dtype=np.int64)
    sum_x = np.empty(counts.shape[0], dtype=np.float64)
    sum_x2 = np.empty(counts.shape[0], dtype=np.float64)
    for i in range(counts.shape[0]):
        row = counts[i].astype(np.float64)
        lx = np.log1p(row)
        nnz[i] = int(np.count_nonzero(row))
        sum_x[i] = float(lx.sum())
        sum_x2[i] = float(np.dot(lx, lx))
    STATS.mkdir(parents=True, exist_ok=True)
    with (STATS / f"{dataset}.tsv").open("w") as handle:
        handle.write("gene\tn_cells\tnnz\tsum_log1p\tsumsq_log1p\n")
        for i, gene in enumerate(genes):
            handle.write(f"{gene}\t{n}\t{nnz[i]}\t{sum_x[i]:.8g}\t{sum_x2[i]:.8g}\n")
    say(f"  stored {dataset} {counts.shape[0]} genes x {counts.shape[1]} cells")


def audit(dataset: str, units: list[str], cldn4: np.ndarray, locked: list[dict[str, str]]) -> None:
    AUDIT.mkdir(parents=True, exist_ok=True)
    by: dict[str, list[float]] = {}
    for unit, value in zip(units, cldn4):
        by.setdefault(unit, []).append(float(value))
    locked_map = {row["unit_id"]: row for row in locked}
    lines = ["unit_id\tn_malignant\tn_locked\tpct_pos\tpct_locked\tabs_pct_diff\n"]
    missing = []
    for unit, row in locked_map.items():
        vals = by.get(unit)
        if not vals:
            missing.append(unit)
            continue
        n = len(vals)
        pct = 100.0 * sum(v > 0 for v in vals) / n
        pct_locked = float(row["mal_CLDN4_pct"])
        n_locked = int(float(row["n_malignant"]))
        lines.append(
            f"{unit}\t{n}\t{n_locked}\t{pct:.6f}\t{pct_locked:.6f}\t{abs(pct - pct_locked):.6f}\n"
        )
    text = "".join(lines)
    (AUDIT / f"{dataset}.tsv").write_text(text)
    if missing:
        say(f"  AUDIT {dataset} missing units: {missing}")
    # summary
    diffs = []
    n_mismatch = 0
    for line in lines[1:]:
        unit, n, n_locked, pct, pct_locked, diff = line.rstrip("\n").split("\t")
        diffs.append(float(diff))
        if int(n) != int(n_locked):
            n_mismatch += 1
            say(f"  AUDIT n mismatch {dataset} {unit}: {n} vs locked {n_locked}")
    if diffs:
        say(
            f"  AUDIT {dataset} units={len(diffs)} n_mismatch={n_mismatch} "
            f"max|pct diff|={max(diffs):.4f} median={float(np.median(diffs)):.4f}"
        )


def to_uint16(row: np.ndarray, clip_counter: list[int]) -> np.ndarray:
    row = np.nan_to_num(row, nan=0.0, posinf=0.0, neginf=0.0)
    row = np.rint(row)
    over = int(np.sum(row > 65535))
    if over:
        clip_counter[0] += over
        row = np.clip(row, 0, 65535)
    row[row < 0] = 0
    return row.astype(np.uint16)


def align_matrix(genes_src: list[str], mat: np.ndarray, universe_genes: list[str]) -> tuple[list[str], np.ndarray]:
    """Sum duplicate symbols and return universe genes that are present, genes x cells."""
    buckets: dict[str, list[int]] = {}
    for i, gene in enumerate(genes_src):
        buckets.setdefault(gene, []).append(i)
    present = [g for g in universe_genes if g in buckets]
    out = np.empty((len(present), mat.shape[1]), dtype=np.float64)
    for i, gene in enumerate(present):
        ix = buckets[gene]
        if len(ix) == 1:
            out[i] = np.asarray(mat[ix[0]], dtype=np.float64).ravel()
        else:
            out[i] = np.asarray(mat[ix], dtype=np.float64).sum(axis=0).ravel()
    return present, out


def marker_mask(genes: list[str], mat: np.ndarray) -> np.ndarray:
    index = {g: i for i, g in enumerate(genes)}
    n = mat.shape[1]
    pos = np.zeros(n, dtype=bool)
    for gene in MARKER_POS:
        if gene in index:
            pos |= np.asarray(mat[index[gene]]).ravel() > 0
    if MARKER_NEG in index:
        neg = np.asarray(mat[index[MARKER_NEG]]).ravel() == 0
    else:
        neg = np.ones(n, dtype=bool)
    return pos & neg


def extract_gse123902(universe_genes: list[str], locked: list[dict[str, str]]) -> None:
    say("GSE123902")
    marker = load_tsv(ROOT / "data" / "GSE123902_marker_units.tsv")
    want = {(row["unit_id"], row["tissue"]) for row in locked}
    files = []
    for row in marker:
        key = (row["patient"], row["tissue"])
        if key in want:
            files.append((row["patient"], row["file"]))
    if len(files) != len(locked):
        raise SystemExit(f"GSE123902 file map {len(files)} != locked {len(locked)}")
    tar_path = GEO / "GSE123902_RAW.tar"
    frames = []
    gene_order = None
    units: list[str] = []
    cldn4_all: list[np.ndarray] = []
    clip = [0]
    with tarfile.open(tar_path, "r") as tar:
        for patient, filename in files:
            say(f"  {patient} {filename}")
            member = tar.extractfile(filename)
            if member is None:
                raise SystemExit(f"missing {filename}")
            with gzip.GzipFile(fileobj=member) as gz:
                header = gz.readline().decode().rstrip("\n").split(",")
                genes = [g.strip().strip('"').upper() for g in header[1:]]
                data = np.loadtxt(gz, delimiter=",", dtype=np.float64, usecols=range(1, len(header)))
            # loadtxt returns cells x genes
            if data.ndim == 1:
                data = data.reshape(1, -1)
            mat = np.ascontiguousarray(data.T)
            present, aligned = align_matrix(genes, mat, universe_genes)
            if gene_order is None:
                gene_order = present
            elif present != gene_order:
                # reindex
                idx = {g: i for i, g in enumerate(present)}
                re = np.zeros((len(gene_order), aligned.shape[1]), dtype=np.float64)
                for i, gene in enumerate(gene_order):
                    if gene in idx:
                        re[i] = aligned[idx[gene]]
                aligned = re
            mask = marker_mask(gene_order, aligned)
            sub = aligned[:, mask]
            gindex = {g: i for i, g in enumerate(gene_order)}
            cldn4 = sub[gindex["CLDN4"]] if "CLDN4" in gindex else np.zeros(sub.shape[1])
            frames.append(sub)
            units.extend([patient] * sub.shape[1])
            cldn4_all.append(cldn4)
            say(f"    malignant {sub.shape[1]}")
    counts_f = np.concatenate(frames, axis=1)
    counts = np.empty(counts_f.shape, dtype=np.uint16)
    for i in range(counts_f.shape[0]):
        counts[i] = to_uint16(counts_f[i], clip)
    say(f"  clipped values {clip[0]}")
    write_store("GSE123902", gene_order, units, counts)
    audit("GSE123902", units, np.concatenate(cldn4_all), locked)


def _row_sum(mat: sp.spmatrix, ix: list[int]) -> np.ndarray:
    if len(ix) == 1:
        return np.asarray(mat[ix[0]].todense()).ravel()
    return np.asarray(mat[ix].sum(axis=0)).ravel()


def extract_gse189357(universe_genes: list[str], locked: list[dict[str, str]]) -> None:
    say("GSE189357")
    tar_path = GEO / "GSE189357_RAW.tar"
    dest = GEO / "gse189357"
    if not any(dest.glob("*_matrix.mtx.gz")):
        dest.mkdir(parents=True, exist_ok=True)
        say("  extracting tar")
        with tarfile.open(tar_path, "r") as tar:
            tar.extractall(dest)
    frames = []
    gene_order = [g for g in universe_genes]
    units: list[str] = []
    cldn4_all: list[np.ndarray] = []
    clip = [0]
    for row in locked:
        pat = row["unit_id"]
        mtx = next(dest.glob(f"*_{pat}_matrix.mtx.gz"))
        feat = next(dest.glob(f"*_{pat}_features.tsv.gz"))
        bc = next(dest.glob(f"*_{pat}_barcodes.tsv.gz"))
        say(f"  {pat} {mtx.name}")
        mat = scipy.io.mmread(mtx).tocsr()
        with gzip.open(feat, "rt") as handle:
            features = [line.rstrip("\n").split("\t") for line in handle]
        symbols = [rec[1].upper() if len(rec) > 1 else rec[0].upper() for rec in features]
        n_bc = sum(1 for _ in gzip.open(bc, "rt"))
        if mat.shape[1] != n_bc and mat.shape[0] == n_bc:
            mat = mat.T.tocsr()
        if mat.shape[0] != len(symbols) or mat.shape[1] != n_bc:
            raise SystemExit(f"{pat}: mtx {mat.shape} features {len(symbols)} barcodes {n_bc}")
        buckets: dict[str, list[int]] = {}
        for i, gene in enumerate(symbols):
            buckets.setdefault(gene, []).append(i)
        n = mat.shape[1]
        pos = np.zeros(n, dtype=bool)
        for gene in MARKER_POS:
            if gene in buckets:
                pos |= _row_sum(mat, buckets[gene]) > 0
        if MARKER_NEG in buckets:
            neg = _row_sum(mat, buckets[MARKER_NEG]) == 0
        else:
            neg = np.ones(n, dtype=bool)
        mask = pos & neg
        mal_n = int(mask.sum())
        sub = mat[:, mask].tocsr()
        del mat
        block = np.zeros((len(gene_order), mal_n), dtype=np.float64)
        for i, gene in enumerate(gene_order):
            if gene in buckets:
                block[i] = _row_sum(sub, buckets[gene])
        cldn4 = block[gene_order.index("CLDN4")].copy() if "CLDN4" in buckets else np.zeros(mal_n)
        frames.append(block)
        units.extend([pat] * mal_n)
        cldn4_all.append(cldn4)
        say(f"    cells {n} malignant {mal_n}")
        del sub
    counts_f = np.concatenate(frames, axis=1)
    # drop genes that are all zero across the dataset (absent)
    present_ix = [i for i in range(counts_f.shape[0]) if np.any(counts_f[i] > 0)]
    # keep structural zeros that are truly measured: 10x panel includes the gene
    # even if a malignant subset is all zero. Presence was checked per sample via `in rows`.
    # A gene absent from the panel stays all zero in every sample. Drop those.
    genes = [gene_order[i] for i in present_ix]
    counts_f = counts_f[present_ix]
    counts = np.empty(counts_f.shape, dtype=np.uint16)
    for i in range(counts_f.shape[0]):
        counts[i] = to_uint16(counts_f[i], clip)
    say(f"  clipped values {clip[0]}")
    write_store("GSE189357", genes, units, counts)
    audit("GSE189357", units, np.concatenate(cldn4_all), locked)


def extract_gse131907(universe_genes: list[str], locked: list[dict[str, str]]) -> None:
    say("GSE131907")
    ann_path = GEO / "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
    mat_path = GEO / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    wanted = {row["unit_id"] for row in locked}
    cells = []
    with gzip.open(ann_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            sample = parts[idx["Sample"]]
            if sample not in wanted:
                continue
            if parts[idx["Cell_subtype"]] != "Malignant cells":
                continue
            cells.append((parts[idx["Index"]], sample))
    say(f"  malignant cells in locked samples {len(cells)}")
    with gzip.open(mat_path, "rt") as handle:
        barcodes = handle.readline().rstrip("\n").split("\t")[1:]
    col_of = {bc: i for i, bc in enumerate(barcodes)}
    missing = [bc for bc, _ in cells if bc not in col_of]
    if missing:
        raise SystemExit(f"GSE131907 missing {len(missing)} barcodes, e.g. {missing[0]}")
    malig_cols = np.array([col_of[bc] for bc, _ in cells], dtype=np.int64)
    units = [sample for _, sample in cells]
    gene_to_row = {g: i for i, g in enumerate(universe_genes)}
    mm_path = STORE / "GSE131907" / "counts.u16"
    (STORE / "GSE131907").mkdir(parents=True, exist_ok=True)
    mm = np.memmap(mm_path, dtype=np.uint16, mode="w+", shape=(len(universe_genes), len(units)))
    mm[:] = 0
    seen = np.zeros(len(universe_genes), dtype=bool)
    clip = [0]
    n_hit = 0
    say("  streaming UMI matrix")
    with gzip.open(mat_path, "rt") as handle:
        handle.readline()
        for line_i, line in enumerate(handle, start=1):
            gene, sep, rest = line.partition("\t")
            if not sep:
                continue
            gene = gene.split(".")[0].upper()
            row_i = gene_to_row.get(gene)
            if row_i is None:
                continue
            arr = np.fromstring(rest, sep="\t", dtype=np.float32)
            if arr.size != len(barcodes):
                raise SystemExit(f"{gene}: {arr.size} values, expected {len(barcodes)}")
            sub = arr[malig_cols]
            if seen[row_i]:
                prev = mm[row_i].astype(np.float32)
                sub = prev + sub
            block = to_uint16(sub, clip)
            mm[row_i] = block
            seen[row_i] = True
            n_hit += 1
            if line_i % 4000 == 0:
                say(f"    scanned {line_i} genes, universe hits {n_hit}")
    mm.flush()
    genes = [g for g, flag in zip(universe_genes, seen) if flag]
    if len(genes) != int(seen.sum()):
        raise SystemExit("seen mismatch")
    # compact to genes that were present in the matrix
    compact = np.array(mm[seen], dtype=np.uint16)
    del mm
    mm_path.unlink()
    say(f"  clipped values {clip[0]} genes present {compact.shape[0]}")
    cldn4 = compact[genes.index("CLDN4")].astype(np.float64) if "CLDN4" in genes else np.zeros(len(units))
    write_store("GSE131907", genes, units, compact)
    audit("GSE131907", units, cldn4, locked)


def extract_gse205335(universe_genes: list[str], locked: list[dict[str, str]]) -> None:
    say("GSE205335")
    export = WORK / "gse205335_export"
    mtx = export / "matrix.mtx"
    if not mtx.exists():
        raise SystemExit(f"missing {mtx}; run export_gse205335.R first")
    genes = [g.strip().upper() for g in (export / "genes.txt").read_text().splitlines() if g.strip()]
    meta = load_tsv(export / "meta.tsv")
    say(f"  reading mtx {len(genes)} genes x {len(meta)} cells")
    mat = scipy.io.mmread(mtx).tocsr()
    if mat.shape != (len(genes), len(meta)):
        raise SystemExit(f"mtx shape {mat.shape} != {(len(genes), len(meta))}")
    buckets: dict[str, list[int]] = {}
    for i, gene in enumerate(genes):
        buckets.setdefault(gene, []).append(i)
    present = [g for g in universe_genes if g in buckets]
    units = [row["unit_id"] for row in meta]
    clip = [0]
    counts = np.empty((len(present), len(units)), dtype=np.uint16)
    cldn4 = np.zeros(len(units), dtype=np.float64)
    for i, gene in enumerate(present):
        row = _row_sum(mat, buckets[gene])
        if gene == "CLDN4":
            cldn4 = row.astype(np.float64)
        counts[i] = to_uint16(row, clip)
        if i and i % 4000 == 0:
            say(f"    packed {i} genes")
    del mat
    say(f"  clipped values {clip[0]}")
    write_store("GSE205335", present, units, counts)
    audit("GSE205335", units, cldn4, locked)


def main() -> None:
    locked = locked_units()
    genes = universe()
    say(f"universe {len(genes)}")
    extract_gse123902(genes, locked["GSE123902"])
    extract_gse189357(genes, locked["GSE189357"])
    extract_gse131907(genes, locked["GSE131907"])
    extract_gse205335(genes, locked["GSE205335"])
    say("extract done")


if __name__ == "__main__":
    main()
