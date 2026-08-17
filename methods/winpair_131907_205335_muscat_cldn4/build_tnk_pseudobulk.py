#!/usr/bin/env python3
"""Patient-level T/NK UMI-sum matrices for the winning pair only.

GSE131907: author Cell_type in {T lymphocytes, NK cells}.
GSE205335: author lineage.total == T/NK cells; drop normal tissues.
Same patients/samples as PR #320 author-malignant CLDN4 %pos vs T/NK.
"""
from __future__ import annotations

import gzip
import shutil
import tempfile
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
GEO = Path("/tmp/winpair_geo")
TUMOR_ORIGINS = {"tLung", "tL/B", "mLN", "PE", "mBrain"}
NORMAL_205335 = {"Normal Brain", "Normal LN", "Normal Lung"}
MIN_TNK = 20


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


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


def accumulate_from_txt(path: Path, cell_to_patient: dict[str, str]) -> pd.DataFrame:
    cells, fh = stream_txt_header(path)
    n = len(cells)
    patients = sorted(set(cell_to_patient.values()))
    pidx = {p: i for i, p in enumerate(patients)}
    col_patient = np.full(n, -1, dtype=np.int32)
    for i, c in enumerate(cells):
        pt = cell_to_patient.get(c)
        if pt is not None:
            col_patient[i] = pidx[pt]
    use = col_patient >= 0
    use_i = np.where(use)[0]
    assign = col_patient[use]
    log(f"  stream {path.name}: cells={n} kept={int(use.sum())} units={len(patients)}")
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
                log(f"    genes={n_genes} nonzero={len(genes)}")
    finally:
        fh.close()
    return pd.DataFrame(np.vstack(rows) if rows else np.zeros((0, len(patients))), index=genes, columns=patients)


def save_counts(stem: str, counts: pd.DataFrame, meta: pd.DataFrame) -> None:
    counts = counts.groupby(counts.index.astype(str)).sum()
    counts = counts.loc[counts.sum(axis=1) > 0]
    path = DATA / f"{stem}.tsv.gz"
    counts.to_csv(path, sep="\t", compression="gzip")
    meta.to_csv(DATA / f"{stem.replace('_counts', '_meta')}.tsv", sep="\t", index=False)
    log(f"  wrote {path.name} genes={counts.shape[0]} units={counts.shape[1]}")


def build_gse131907() -> None:
    log("=== GSE131907 T/NK ===")
    scores = pd.read_csv(DATA / "GSE131907_samples.tsv", sep="\t")
    keep = scores.loc[
        scores["origin"].isin(TUMOR_ORIGINS) & (scores["n_malignant"] >= 20),
        "sample",
    ].astype(str)
    annot = pd.read_csv(GEO / "GSE131907" / "cell_annotation.txt.gz", sep="\t")
    tnk = annot["Cell_type"].isin(["T lymphocytes", "NK cells"])
    tumor = annot["Sample_Origin"].isin(TUMOR_ORIGINS)
    annot = annot.loc[tnk & tumor & annot["Sample"].astype(str).isin(keep)].copy()
    n_tnk = annot.groupby("Sample").size()
    keep_pt = n_tnk[n_tnk >= MIN_TNK].index.astype(str)
    annot = annot.loc[annot["Sample"].astype(str).isin(keep_pt)]
    cell_to_pt = dict(zip(annot["Index"].astype(str), annot["Sample"].astype(str)))
    log(f"  T/NK cells={len(annot)} samples={len(keep_pt)} (min {MIN_TNK})")
    counts = accumulate_from_txt(GEO / "GSE131907" / "raw_UMI_matrix.txt.gz", cell_to_pt)
    counts = counts.loc[:, [c for c in keep_pt if c in counts.columns]]
    pmeta = (
        annot.groupby("Sample")
        .agg(
            n_tnk=("Index", "size"),
            Sample_Origin=("Sample_Origin", "first"),
            n_T=("Cell_type", lambda s: int((s == "T lymphocytes").sum())),
            n_NK=("Cell_type", lambda s: int((s == "NK cells").sum())),
        )
        .reset_index()
        .rename(columns={"Sample": "patient"})
    )
    pmeta["tnk_rule"] = "author Cell_type in {T lymphocytes, NK cells}; tumor origins; n_mal>=20"
    pmeta["cohort"] = "GSE131907"
    pmeta["unit"] = "sample"
    save_counts("GSE131907_tnk_counts", counts, pmeta)


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


def build_gse205335() -> None:
    log("=== GSE205335 T/NK ===")
    import rdata
    from scipy import sparse

    scores = pd.read_csv(DATA / "GSE205335_patients.tsv", sep="\t")
    keep_patients = set(scores["patient"].astype(str))
    gsm = pd.read_csv(DATA / "GSE205335_gsm_sample_metadata.csv")
    ident = pd.read_csv(GEO / "GSE205335_Lung_IO_CellIdentity.txt.gz", sep="\t")
    ident = ident.merge(gsm[["orig.ident", "patient", "tissue"]], on="orig.ident", how="left")
    tnk = (
        ident["lineage.total"].eq("T/NK cells")
        & ~ident["tissue"].isin(NORMAL_205335)
        & ident["patient"].astype(str).isin(keep_patients)
    )
    ident_t = ident.loc[tnk].copy()
    n_tnk = ident_t.groupby("patient").size()
    keep_pt = n_tnk[n_tnk >= MIN_TNK].index.astype(str)
    ident_t = ident_t.loc[ident_t["patient"].astype(str).isin(keep_pt)]
    barcode_to_pt = dict(zip(ident_t["barcode"].astype(str), ident_t["patient"].astype(str)))
    log(f"  T/NK cells={len(ident_t)} patients={len(keep_pt)} (min {MIN_TNK})")

    with tempfile.TemporaryDirectory(prefix="gse205335-tnk-") as tmp:
        rds_path = Path(tmp) / "matrix.rds"
        log("  decompressing RDS...")
        _ungzip_until_rds(GEO / "GSE205335_Lung_IO_UMI_matrix.rds.gz", rds_path)
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
    log(f"  summed T/NK columns={n_used}")
    counts = pd.DataFrame(acc, index=genes, columns=patients)
    pmeta = (
        ident_t.groupby("patient")
        .agg(
            n_tnk=("barcode", "size"),
            tissue=("tissue", lambda s: ",".join(sorted(set(map(str, s))))),
            n_CD4=("lineage.sub", lambda s: int((s == "CD4+ T cells").sum())),
            n_CD8=("lineage.sub", lambda s: int((s == "CD8+ T cells").sum())),
            n_NK=("lineage.sub", lambda s: int((s == "NK cells").sum())),
        )
        .reset_index()
    )
    pmeta["tnk_rule"] = "author lineage.total==T/NK cells; drop normal tissues"
    pmeta["cohort"] = "GSE205335"
    pmeta["unit"] = "patient"
    save_counts("GSE205335_tnk_counts", counts, pmeta)


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    build_gse131907()
    build_gse205335()
    log("done")


if __name__ == "__main__":
    main()
