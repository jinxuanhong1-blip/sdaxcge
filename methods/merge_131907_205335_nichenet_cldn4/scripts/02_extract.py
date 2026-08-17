#!/usr/bin/env python3
"""Extract NicheNet ligands/receptors + T/NK programs from both UMI matrices."""
from __future__ import annotations

import gzip
import shutil
import sys
import tempfile
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
CACHE = Path("/tmp/merge_131907_205335_nichenet_cldn4")
sys.path.insert(0, str(ROOT / "scripts"))
from gene_sets import CYTOTOXICITY, EXHAUSTION, EXTRA_TNK, IFN, LINEAGE  # noqa: E402

LR = DATA / "lr_network.tsv"


def panel_genes() -> set[str]:
    lr = pd.read_csv(LR, sep="\t")
    genes = set(lr["from"].astype(str)) | set(lr["to"].astype(str))
    for block in (CYTOTOXICITY, IFN, EXHAUSTION, EXTRA_TNK, LINEAGE):
        genes.update(block)
    return genes


def stream_umi_txt(path: Path, keep: set[str]) -> tuple[np.ndarray, dict[str, np.ndarray], np.ndarray, int]:
    print(f"[stream] {path} keep={len(keep)}", flush=True)
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        cells = np.array(header[1:], dtype=object)
        n = len(cells)
        totals = np.zeros(n, dtype=np.float64)
        store: dict[str, np.ndarray] = {}
        n_genes = 0
        for line in fh:
            gene, sep, rest = line.partition("\t")
            if not sep:
                continue
            gene = gene.split(".")[0]
            n_genes += 1
            vals = np.fromstring(rest, sep="\t", dtype=np.float32)
            if vals.size != n:
                raise ValueError(f"column mismatch for {gene}: {vals.size} != {n}")
            totals += vals
            if gene in keep and gene not in store:
                store[gene] = vals
            if n_genes % 5000 == 0:
                print(f"[stream] genes_seen={n_genes} kept={len(store)}", flush=True)
    print(f"[stream] cells={n} genes_in_file={n_genes} genes_kept={len(store)}", flush=True)
    return cells, store, totals, n_genes


def extract_gse205335(path: Path, keep: set[str]) -> tuple[np.ndarray, dict[str, np.ndarray], np.ndarray, list[str]]:
    import rdata
    from scipy import sparse

    with tempfile.TemporaryDirectory(prefix="gse205335-") as tmp:
        matrix_path = path
        if path.suffix == ".gz":
            matrix_path = Path(tmp) / path.stem
            print(f"decompress {path.name}", flush=True)
            with gzip.open(path, "rb") as source, matrix_path.open("wb") as dest:
                shutil.copyfileobj(source, dest, 16 * 1024 * 1024)
        print("read RDS", flush=True)
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", message='Missing constructor for R class "dgCMatrix"'
            )
            obj = rdata.read_rds(matrix_path)
    expected = {"i", "p", "Dim", "Dimnames", "x"}
    if not expected.issubset(vars(obj)):
        raise TypeError("RDS is not the expected Matrix::dgCMatrix object")
    genes = np.asarray(obj.Dimnames[0], dtype=str)
    barcodes = np.asarray(obj.Dimnames[1], dtype=str)
    print(f"build CSC {tuple(obj.Dim)}", flush=True)
    matrix = sparse.csc_matrix((obj.x, obj.i, obj.p), shape=tuple(obj.Dim), copy=False)
    library_umi = np.asarray(matrix.sum(axis=0)).ravel()
    print("convert CSR for row extract", flush=True)
    matrix = matrix.tocsr()
    name_to_row = {g: i for i, g in enumerate(genes)}
    extracted: dict[str, np.ndarray] = {}
    for gene in sorted(keep):
        row = name_to_row.get(gene)
        if row is None:
            continue
        extracted[gene] = np.asarray(matrix.getrow(row).toarray()).ravel().astype(np.float32)
    print(f"extracted {len(extracted)} / {len(keep)} genes", flush=True)
    return barcodes, extracted, library_umi, genes.tolist()


def write_panel(barcodes, store, totals, dest: Path) -> None:
    df = pd.DataFrame(store, index=pd.Index(barcodes, name="barcode"))
    df["nCount_RNA"] = totals
    dest.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(dest)
    print(f"wrote {dest} shape={df.shape} bytes={dest.stat().st_size}", flush=True)


def main() -> int:
    want = panel_genes()
    print(f"panel requested: {len(want)}", flush=True)

    p131 = CACHE / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    cells, store, totals, n_genes = stream_umi_txt(p131, want)
    write_panel(cells, store, totals, CACHE / "gse131907_panel.parquet")
    (CACHE / "gse131907_extract_n.txt").write_text(f"n_genes={n_genes}\n")

    p205 = CACHE / "GSE205335_Lung_IO_UMI_matrix.rds.gz"
    barcodes, store2, totals2, all_genes = extract_gse205335(p205, want)
    write_panel(barcodes, store2, totals2, CACHE / "gse205335_panel.parquet")
    (CACHE / "gse205335_n_genes.txt").write_text(f"n_genes={len(all_genes)}\n")
    print("ok", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
