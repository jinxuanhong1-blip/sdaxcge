#!/usr/bin/env python3
"""Build concordant-4 cell-class pseudobulk CPM profiles.

See PROTOCOL.md. Raw matrices are downloaded under data/raw and are not
committed. Outputs are results/reference/*.tsv.
"""
from __future__ import annotations

import gzip
import struct
import subprocess
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import io as spio
from scipy import sparse

ROOT = Path(__file__).resolve().parent
RAW = ROOT / "data" / "raw"
REF = ROOT / "results" / "reference"
CLASSES = [
    "malignant", "cd8", "cd4", "nk", "b", "myeloid", "endothelial", "fibroblast",
]
MIN_CELLS = 30

GSE131907_CD8 = {
    "Exhausted CD8+ T", "CD8 low T", "Cytotoxic CD8+ T", "Naive CD8+ T",
}
GSE131907_CD4 = {"CD4+ Th", "Naive CD4+ T", "Treg", "Exhausted Tfh"}


def _download(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        print(f"have {dest.name}", flush=True)
        return dest
    print(f"GET {url}", flush=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    urllib.request.urlretrieve(url, tmp)
    tmp.rename(dest)
    return dest


def _symbol_upper(x: pd.Index | pd.Series) -> pd.Index:
    return pd.Index(pd.Series(x).astype(str).str.upper().str.strip())


def class_cpm(sums: pd.DataFrame, ncells: pd.Series, libs: pd.Series) -> pd.DataFrame:
    sums = sums.copy()
    sums.index = _symbol_upper(sums.index)
    sums = sums.groupby(level=0).sum()
    cpm = pd.DataFrame(index=sums.index, columns=CLASSES, dtype=float)
    for k in CLASSES:
        lib = float(libs.get(k, 0.0))
        n = int(ncells.get(k, 0))
        if n >= MIN_CELLS and lib > 0 and k in sums.columns:
            cpm[k] = sums[k].to_numpy(dtype=float) / lib * 1e6
        else:
            cpm[k] = np.nan
    return cpm


def marker_labels(expr: pd.DataFrame) -> pd.Series:
    """expr: cells x genes, raw counts, gene symbols already upper-case."""

    def col(name: str) -> np.ndarray:
        if name not in expr.columns:
            return np.zeros(expr.shape[0], dtype=float)
        return expr[name].to_numpy(dtype=float)

    cd3 = (col("CD3D") > 0) | (col("CD3E") > 0)
    cd8 = col("CD8A") + col("CD8B")
    cd4 = col("CD4")
    nk = (col("NKG7") > 0) | (col("GNLY") > 0)
    b = (col("MS4A1") > 0) | (col("CD79A") > 0)
    my = (col("LYZ") > 0) | (col("CD68") > 0) | (col("CD14") > 0)
    en = (col("PECAM1") > 0) | (col("VWF") > 0)
    fi = (col("COL1A1") > 0) | (col("DCN") > 0)
    epi = ((col("EPCAM") > 0) | (col("KRT19") > 0)) & (col("PTPRC") == 0)
    lab = np.array([""] * expr.shape[0], dtype=object)
    both = cd3 & (cd8 > 0) & (cd4 > 0)
    only8 = cd3 & (cd8 > 0) & (cd4 == 0)
    only4 = cd3 & (cd4 > 0) & (cd8 == 0)
    lab[only8] = "cd8"
    lab[only4] = "cd4"
    lab[both & (cd8 > cd4)] = "cd8"
    lab[both & (cd4 > cd8)] = "cd4"
    rest = lab == ""
    lab[rest & nk] = "nk"
    rest = lab == ""
    lab[rest & b] = "b"
    rest = lab == ""
    lab[rest & my] = "myeloid"
    rest = lab == ""
    lab[rest & en] = "endothelial"
    rest = lab == ""
    lab[rest & fi] = "fibroblast"
    rest = lab == ""
    lab[rest & epi] = "malignant"
    return pd.Series(lab, index=expr.index)


def _pool_labeled(counts_cells_by_genes: pd.DataFrame, labels: pd.Series,
                  sums: pd.DataFrame, ncells: pd.Series, libs: pd.Series) -> None:
    counts_cells_by_genes = counts_cells_by_genes.groupby(level=0, axis=1).sum()
    lib = counts_cells_by_genes.sum(axis=1)
    for k in CLASSES:
        m = labels == k
        if not m.any():
            continue
        block = counts_cells_by_genes.loc[m]
        summed = block.sum(axis=0)
        for g, v in summed.items():
            sums.at[g, k] = float(sums.at[g, k]) + float(v) if g in sums.index else float(v)
            if g not in sums.index:
                sums.loc[g, k] = float(v)
        ncells[k] = int(ncells[k]) + int(m.sum())
        libs[k] = float(libs[k]) + float(lib.loc[m].sum())


def accumulate_dense(df_cells_genes: pd.DataFrame, sums: pd.DataFrame,
                     ncells: pd.Series, libs: pd.Series) -> int:
    df_cells_genes = df_cells_genes.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    df_cells_genes.columns = _symbol_upper(df_cells_genes.columns)
    grouped = df_cells_genes.T.groupby(level=0).sum().T
    labels = marker_labels(grouped)
    before = {k: int(ncells[k]) for k in CLASSES}
    lib = grouped.sum(axis=1)
    for k in CLASSES:
        m = labels.to_numpy() == k
        if not np.any(m):
            continue
        summed = grouped.iloc[np.flatnonzero(m)].sum(axis=0)
        if k not in sums.columns or sums.empty:
            sums[k] = summed
        else:
            sums[k] = sums[k].add(summed, fill_value=0.0)
        ncells[k] = int(ncells[k]) + int(np.sum(m))
        libs[k] = float(libs[k]) + float(lib.iloc[np.flatnonzero(m)].sum())
    return int(sum(int(ncells[k]) - before[k] for k in CLASSES))


def build_gse123902() -> tuple[pd.DataFrame, pd.Series]:
    out_cpm = REF / "cpm_GSE123902.tsv"
    out_n = REF / "ncells_GSE123902.tsv"
    if out_cpm.exists() and out_n.exists():
        return pd.read_csv(out_cpm, sep="\t", index_col=0), pd.read_csv(out_n, sep="\t", index_col=0).squeeze()
    rows = []
    # GSM3516662 .. GSM3516678 from the GEO file list; skip NORMAL.
    names = {
        "GSM3516662": "GSM3516662_MSK_LX653_PRIMARY_TUMOUR_dense.csv.gz",
        "GSM3516663": "GSM3516663_MSK_LX661_PRIMARY_TUMOUR_dense.csv.gz",
        "GSM3516664": "GSM3516664_MSK_LX666_METASTASIS_dense.csv.gz",
        "GSM3516665": "GSM3516665_MSK_LX675_PRIMARY_TUMOUR_dense.csv.gz",
        "GSM3516667": "GSM3516667_MSK_LX676_PRIMARY_TUMOUR_dense.csv.gz",
        "GSM3516668": "GSM3516668_MSK_LX255B_METASTASIS_dense.csv.gz",
        "GSM3516669": "GSM3516669_MSK_LX679_PRIMARY_TUMOUR_dense.csv.gz",
        "GSM3516670": "GSM3516670_MSK_LX680_PRIMARY_TUMOUR_dense.csv.gz",
        "GSM3516671": "GSM3516671_MSK_LX681_METASTASIS_dense.csv.gz",
        "GSM3516672": "GSM3516672_MSK_LX682_PRIMARY_TUMOUR_dense.csv.gz",
        "GSM3516674": "GSM3516674_MSK_LX684_PRIMARY_TUMOUR_dense.csv.gz",
        "GSM3516677": "GSM3516677_MSK_LX699_METASTASIS_dense.csv.gz",
        "GSM3516678": "GSM3516678_MSK_LX701_METASTASIS_dense.csv.gz",
    }
    sums = pd.DataFrame(columns=CLASSES, dtype=float)
    ncells = pd.Series(0, index=CLASSES, dtype=int)
    libs = pd.Series(0.0, index=CLASSES, dtype=float)
    inventory = []
    for gsm, name in names.items():
        url = f"https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3516nnn/{gsm}/suppl/{name}"
        path = _download(url, RAW / name)
        df = pd.read_csv(path, index_col=0)
        n_kept = accumulate_dense(df, sums, ncells, libs)
        inventory.append({"sample": name, "n_cells_matrix": df.shape[0], "n_cells_assigned": n_kept})
        print(f"GSE123902 {name} cells={df.shape[0]} assigned_new={n_kept}", flush=True)
        del df
    sums = sums.fillna(0.0)
    cpm = class_cpm(sums, ncells, libs)
    REF.mkdir(parents=True, exist_ok=True)
    cpm.to_csv(out_cpm, sep="\t")
    ncells.to_csv(out_n, sep="\t", header=["n_cells"])
    pd.DataFrame(inventory).to_csv(REF / "inventory_GSE123902.tsv", sep="\t", index=False)
    print("GSE123902 n", ncells.to_dict(), flush=True)
    return cpm, ncells


def _gse131907_labels(ann: pd.DataFrame) -> pd.Series:
    ct = ann["Cell_type"]
    st = ann["Cell_subtype"]
    lab = pd.Series("", index=ann.index)
    lab[(ct == "T lymphocytes") & st.isin(GSE131907_CD8)] = "cd8"
    lab[(ct == "T lymphocytes") & st.isin(GSE131907_CD4)] = "cd4"
    lab[(ct == "NK cells") & (st == "NK")] = "nk"
    lab[st.isin(["tS1", "tS2", "tS3"])] = "malignant"
    lab[(lab == "") & (ct == "B lymphocytes")] = "b"
    lab[(lab == "") & (ct == "Myeloid cells")] = "myeloid"
    lab[(lab == "") & (ct == "Endothelial cells")] = "endothelial"
    lab[(lab == "") & (ct == "Fibroblasts")] = "fibroblast"
    return lab


def build_gse131907() -> tuple[pd.DataFrame, pd.Series]:
    out_cpm = REF / "cpm_GSE131907.tsv"
    out_n = REF / "ncells_GSE131907.tsv"
    if out_cpm.exists() and out_n.exists():
        return pd.read_csv(out_cpm, sep="\t", index_col=0), pd.read_csv(out_n, sep="\t", index_col=0).squeeze()
    ann_path = _download(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_cell_annotation.txt.gz",
        RAW / "GSE131907_annotation.txt.gz",
    )
    mat_path = _download(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz",
        RAW / "GSE131907_raw_UMI.txt.gz",
    )
    ann = pd.read_csv(ann_path, sep="\t")
    ann = ann.loc[ann["Sample_Origin"] == "tLung"].copy()
    ann["label"] = _gse131907_labels(ann)
    ann = ann.loc[ann["label"] != ""]
    print("reading GSE131907 header", flush=True)
    with gzip.open(mat_path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
    index_of = {c: i for i, c in enumerate(header)}
    cls = np.full(len(header), -1, dtype=np.int16)
    n_match = 0
    for _, row in ann.iterrows():
        i = index_of.get(row["Index"])
        if i is None:
            continue
        cls[i] = np.int16(CLASSES.index(row["label"]))
        n_match += 1
    print(f"GSE131907 labeled tLung cells matched to matrix: {n_match}", flush=True)
    if n_match < 1000:
        raise SystemExit("GSE131907 barcode match failed")
    cls_path = RAW / "GSE131907_class.bin"
    with cls_path.open("wb") as fh:
        fh.write(struct.pack("<i", len(header)))
        cls.astype("<i2").tofile(fh)
    bin_path = ROOT / "extract_sums"
    if not bin_path.exists():
        subprocess.check_call(["gcc", "-O3", "-o", str(bin_path), str(ROOT / "extract_sums.c"), "-lz"])
    sums_path = RAW / "GSE131907_sums.tsv"
    log_path = RAW / "GSE131907_extract.log"
    print("pooling GSE131907 (full matrix scan)", flush=True)
    with log_path.open("w") as log:
        subprocess.check_call([str(bin_path), str(mat_path), str(cls_path), str(sums_path)], stderr=log)
    lib_line = ""
    for line in log_path.read_text().splitlines():
        if line.startswith("LIB"):
            lib_line = line
    libs_v = [float(x) for x in lib_line.split()[1:]]
    libs = pd.Series(libs_v, index=CLASSES)
    ncells = pd.Series([(cls == i).sum() for i in range(len(CLASSES))], index=CLASSES)
    sums = pd.read_csv(sums_path, sep="\t", header=None, names=["gene", *CLASSES])
    sums = sums.set_index("gene")
    cpm = class_cpm(sums, ncells, libs)
    cpm.to_csv(out_cpm, sep="\t")
    ncells.to_csv(out_n, sep="\t", header=["n_cells"])
    print("GSE131907 n", ncells.to_dict(), "lib", libs.to_dict(), flush=True)
    return cpm, ncells


def build_gse205335() -> tuple[pd.DataFrame, pd.Series]:
    out_cpm = REF / "cpm_GSE205335.tsv"
    out_n = REF / "ncells_GSE205335.tsv"
    if out_cpm.exists() and out_n.exists():
        return pd.read_csv(out_cpm, sep="\t", index_col=0), pd.read_csv(out_n, sep="\t", index_col=0).squeeze()
    mat = _download(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_UMI_matrix.rds.gz",
        RAW / "GSE205335_UMI.rds.gz",
    )
    ids = _download(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_CellIdentity.txt.gz",
        RAW / "GSE205335_CellIdentity.txt.gz",
    )
    plain = RAW / "GSE205335_UMI.rds"
    if not plain.exists():
        print("gunzip GSE205335 RDS", flush=True)
        subprocess.check_call(["gunzip", "-c", str(mat)], stdout=plain.open("wb"))
    id_txt = RAW / "GSE205335_CellIdentity.txt"
    if not id_txt.exists():
        subprocess.check_call(["gunzip", "-c", str(ids)], stdout=id_txt.open("wb"))
    subprocess.check_call(
        ["Rscript", str(ROOT / "aggregate_gse205335.R"), str(plain), str(id_txt), str(RAW)]
    )
    sums = pd.read_csv(RAW / "gse205335_sums.tsv", sep="\t").set_index("gene")
    meta = pd.read_csv(RAW / "gse205335_ncells.tsv", sep="\t")
    ncells = meta.set_index("class")["n_cells"].reindex(CLASSES).fillna(0).astype(int)
    libs = meta.set_index("class")["library_umi"].reindex(CLASSES).fillna(0.0)
    cpm = class_cpm(sums[CLASSES], ncells, libs)
    cpm.to_csv(out_cpm, sep="\t")
    ncells.to_csv(out_n, sep="\t", header=["n_cells"])
    print("GSE205335 n", ncells.to_dict(), flush=True)
    return cpm, ncells


def _load_10x(gsm: str, td: str) -> tuple[sparse.csr_matrix, pd.Index]:
    base = f"https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5699nnn/{gsm}/suppl"
    feat = _download(f"{base}/{gsm}_{td}_features.tsv.gz", RAW / f"{gsm}_{td}_features.tsv.gz")
    bar = _download(f"{base}/{gsm}_{td}_barcodes.tsv.gz", RAW / f"{gsm}_{td}_barcodes.tsv.gz")
    mtx = _download(f"{base}/{gsm}_{td}_matrix.mtx.gz", RAW / f"{gsm}_{td}_matrix.mtx.gz")
    features = pd.read_csv(feat, sep="\t", header=None)
    symbols = features.iloc[:, 1].astype(str).str.upper()
    barcodes = pd.read_csv(bar, sep="\t", header=None).iloc[:, 0]
    X = spio.mmread(mtx).tocsr()
    if X.shape[0] != len(symbols) and X.shape[1] == len(symbols):
        X = X.T.tocsr()
    if X.shape[1] != len(barcodes):
        raise SystemExit(f"{td} barcode/matrix mismatch {X.shape} vs {len(barcodes)}")
    # collapse duplicate symbols
    codes, uniques = pd.factorize(symbols, sort=False)
    acc = sparse.csr_matrix(
        (np.ones(len(codes), dtype=np.float64), (codes, np.arange(len(codes)))),
        shape=(len(uniques), len(codes)),
    )
    Xg = acc @ X
    return Xg.tocsr(), pd.Index(uniques)


def build_gse189357() -> tuple[pd.DataFrame, pd.Series]:
    out_cpm = REF / "cpm_GSE189357.tsv"
    out_n = REF / "ncells_GSE189357.tsv"
    if out_cpm.exists() and out_n.exists():
        return pd.read_csv(out_cpm, sep="\t", index_col=0), pd.read_csv(out_n, sep="\t", index_col=0).squeeze()
    sums = None
    ncells = pd.Series(0, index=CLASSES, dtype=int)
    libs = pd.Series(0.0, index=CLASSES, dtype=float)
    genes = None
    inventory = []
    for i in range(1, 10):
        gsm = f"GSM{5699776 + i}"
        td = f"TD{i}"
        X, symbols = _load_10x(gsm, td)
        marker_names = [
            "CD3D", "CD3E", "CD8A", "CD8B", "CD4", "NKG7", "GNLY", "MS4A1", "CD79A",
            "LYZ", "CD68", "CD14", "PECAM1", "VWF", "COL1A1", "DCN", "EPCAM", "KRT19", "PTPRC",
        ]
        pos = {g: int(np.flatnonzero(symbols == g)[0]) for g in marker_names if (symbols == g).any()}
        sub = pd.DataFrame(
            {g: np.asarray(X[i].todense()).ravel() for g, i in pos.items()}
        )
        labels = marker_labels(sub)
        lab = labels.to_numpy()
        ind = sparse.csc_matrix(
            (np.ones(X.shape[1]), (np.arange(X.shape[1]),
                                   np.array([CLASSES.index(x) if x in CLASSES else 0 for x in lab]))),
            shape=(X.shape[1], len(CLASSES)),
        )
        # zero-out unassigned by rebuilding indicator properly
        cols = []
        data_cols = []
        for ki, k in enumerate(CLASSES):
            m = np.flatnonzero(lab == k)
            ncells[k] = int(ncells[k]) + int(m.size)
            cols.append(m)
        # class sums
        class_sums = np.zeros((X.shape[0], len(CLASSES)), dtype=np.float64)
        col_lib = np.asarray(X.sum(axis=0)).ravel()
        for ki, k in enumerate(CLASSES):
            m = np.flatnonzero(lab == k)
            if m.size == 0:
                continue
            class_sums[:, ki] = np.asarray(X[:, m].sum(axis=1)).ravel()
            libs[k] = float(libs[k]) + float(col_lib[m].sum())
        if sums is None:
            sums = class_sums
            genes = symbols
        else:
            if not genes.equals(symbols):
                # align
                cur = pd.DataFrame(class_sums, index=symbols)
                base = pd.DataFrame(sums, index=genes)
                both = base.add(cur, fill_value=0.0)
                genes = both.index
                sums = both.to_numpy()
                class_sums = None
            else:
                sums += class_sums
        n_asg = int((lab != "").sum())
        inventory.append({"sample": td, "n_cells_matrix": int(X.shape[1]), "n_cells_assigned": n_asg})
        print(f"GSE189357 {td} cells={X.shape[1]} assigned={n_asg}", flush=True)
        # free the mtx download to keep disk reasonable
        mtx = RAW / f"{gsm}_{td}_matrix.mtx.gz"
        if mtx.exists():
            mtx.unlink()
        del X
    sums_df = pd.DataFrame(sums, index=genes, columns=CLASSES)
    cpm = class_cpm(sums_df, ncells, libs)
    cpm.to_csv(out_cpm, sep="\t")
    ncells.to_csv(out_n, sep="\t", header=["n_cells"])
    pd.DataFrame(inventory).to_csv(REF / "inventory_GSE189357.tsv", sep="\t", index=False)
    print("GSE189357 n", ncells.to_dict(), flush=True)
    return cpm, ncells


def combine(profiles: dict[str, pd.DataFrame], malignant_datasets: list[str] | None,
            tag: str) -> pd.DataFrame:
    pieces = []
    for ds, cpm in profiles.items():
        block = cpm.copy()
        if malignant_datasets is not None:
            if ds not in malignant_datasets:
                block["malignant"] = np.nan
        pieces.append(block)
    union = pieces[0].index
    for p in pieces[1:]:
        union = union.union(p.index)
    union = union.sort_values() if hasattr(union, "sort_values") else union
    arrs = [p.reindex(union) for p in pieces]
    stacked = np.stack([a.to_numpy(dtype=float) for a in arrs], axis=0)
    with np.errstate(all="ignore"):
        mean = np.nanmean(stacked, axis=0)
    present_ds = []
    for ds, cpm in profiles.items():
        present_ds.append(pd.Series(1, index=cpm.index).reindex(union).fillna(0).to_numpy())
    n_ds = np.sum(np.vstack(present_ds), axis=0)
    out = pd.DataFrame(mean, index=union, columns=CLASSES)
    out = out.loc[n_ds >= 2]
    out = out.dropna(how="all")
    out.to_csv(REF / f"cpm_concordant4_{tag}.tsv", sep="\t")
    return out


def main() -> None:
    REF.mkdir(parents=True, exist_ok=True)
    RAW.mkdir(parents=True, exist_ok=True)
    profiles = {
        "GSE123902": build_gse123902()[0],
        "GSE131907": build_gse131907()[0],
        "GSE205335": build_gse205335()[0],
        "GSE189357": build_gse189357()[0],
    }
    # reload n for the inventory table
    rows = []
    for ds in profiles:
        n = pd.read_csv(REF / f"ncells_{ds}.tsv", sep="\t", index_col=0).squeeze()
        for k, v in n.items():
            rows.append({"dataset": ds, "class": k, "n_cells": int(v)})
    pd.DataFrame(rows).to_csv(REF / "ncells_by_dataset.tsv", sep="\t", index=False)
    primary = combine(profiles, None, "primary")
    author_mal = combine(profiles, ["GSE131907", "GSE205335"], "author_malignant")
    print("primary genes", primary.shape, "author-malignant genes", author_mal.shape, flush=True)
    markers = ["CD8A", "CD8B", "CD4", "NKG7", "MS4A1", "CD79A", "LYZ", "EPCAM",
               "KRT19", "KRT8", "PECAM1", "COL1A1", "CLDN4", "PTPRC", "CD3D"]
    qc = primary.reindex([m for m in markers if m in primary.index])
    qc.to_csv(REF / "marker_cpm_primary.tsv", sep="\t")
    print(qc.round(1).to_string(), flush=True)


if __name__ == "__main__":
    main()
