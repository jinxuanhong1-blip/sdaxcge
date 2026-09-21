#!/usr/bin/env python3
"""Public mouse K/KP/KL lung scRNA: epi vs T labels and program scores.

Reads processed GEO matrices only. Does not download SRA/FASTQ and does not
touch private 8-KL matrices. Does not pool mice across accessions.

Epithelial and T labels are marker rules (plus author labels when GEO
deposits them). Program scores are means of log-normalized expression.
The Cldn4 threshold sweep is computed inside each accession.
"""

from __future__ import annotations

import gzip
import io
import json
import tarfile
from pathlib import Path

import anndata as ad
import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
from scipy import stats

EPI_GENES = ["Epcam", "Krt8", "Krt18", "Krt19", "Sftpc", "Sftpb", "Nkx2-1", "Cdh1"]
T_GENES = ["Cd3d", "Cd3e", "Cd3g", "Cd2", "Trac", "Cd8a", "Cd4"]
NK_GENES = ["Ncr1", "Klrb1c", "Nkg7"]
NHEJ_CORE = ["Xrcc6", "Xrcc5", "Prkdc", "Dclre1c", "Nhej1", "Xrcc4", "Lig4", "Poll", "Polm"]
# Paxx was 1110057K04Rik on older gene builds. Use the symbol that is present.
PAXX_ALIASES = ["Paxx", "1110057K04Rik"]
STING_ALIASES = ["Sting1", "Tmem173"]
CGAS_ALIASES = ["Cgas", "Mb21d1"]
IFN_GENES = [
    "Stat1",
    "Irf7",
    "Isg15",
    "Ifit1",
    "Ifit2",
    "Ifit3",
    "Mx1",
    "Oasl2",
    "Rsad2",
    "Cxcl9",
    "Cxcl10",
    "Cxcl11",
    "Ccl5",
    "Ifng",
    "Ifnb1",
    "Gbp4",
]
# Pre-specified grid. q=0.75 is the primary threshold. gt0 is Cldn4 > 0.
SWEEP_QS = [0.50, 0.60, 0.70, 0.75, 0.80, 0.90, 0.95]
PRIMARY_Q = 0.75
MIN_EPI = 20
MIN_GENES_DROPLET = 200
MAX_PCT_MITO = 20.0

# GEO characteristics retrieved for this run (sample pages, not inferred).
GEO_179501 = {
    "CM2260": ("KL_XTR", "KrasLSL-G12D;Lkb1XTR/XTR;Rosa26LSL-tdTomato/FLPo-ERT2", "Restored", "Female"),
    "CM2319": ("KL_XTR", "KrasLSL-G12D;Lkb1XTR/XTR;Rosa26LSL-tdTomato/FLPo-ERT2", "Restored", "Female"),
    "CM2324": ("KL_XTR", "KrasLSL-G12D;Lkb1XTR/XTR;Rosa26LSL-tdTomato", "Non-Restored", "Female"),
    "CM2328": ("KL_XTR", "KrasLSL-G12D;Lkb1XTR/XTR;Rosa26LSL-tdTomato", "Non-Restored", "Female"),
}
GEO_179502 = {
    "CM0875": ("KL_XTR", "KrasLSL-G12D;Lkb1XTR/XTR;Rosa26LSL-tdTomato", "NonRestored", "Tamoxifen"),
    "CM0879": ("KL_XTR", "KrasLSL-G12D;Lkb1XTR/XTR;Rosa26LSL-tdTomato/FLPo-ERT2", "Restored", "Tamoxifen"),
    "CM0884": ("KL_XTR", "KrasLSL-G12D;Lkb1XTR/XTR;Rosa26LSL-tdTomato/FLPo-ERT2", "Restored", "Tamoxifen"),
    "ZR1932": ("KL_XTR", "KrasLSL-G12D;Lkb1XTR/XTR;Rosa26LSL-tdTomato", "NonRestored", "Vehicle"),
    "ZR1966": ("KL_XTR", "KrasLSL-G12D;Lkb1XTR/XTR;Rosa26LSL-tdTomato/FLPo-ERT2", "NonRestored", "Vehicle"),
    "ZR1969": ("KL_XTR", "KrasLSL-G12D;Lkb1XTR/XTR;Rosa26LSL-tdTomato/FLPo-ERT2", "Restored", "Tamoxifen"),
}


def open_text(path: Path):
    raw = path.open("rb")
    magic = raw.read(2)
    raw.seek(0)
    if magic == b"\x1f\x8b":
        return gzip.open(raw, "rt")
    return io.TextIOWrapper(raw, encoding="utf-8")


def read_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(open_text(path))


def unique_symbols(symbols: list[str]) -> list[str]:
    seen = {}
    out = []
    for i, sym in enumerate(symbols):
        if sym not in seen:
            seen[sym] = i
            out.append(sym)
        else:
            out.append(f"{sym}__dup{i}")
    return out


def csr_from_coo(data, rows, cols, shape) -> sp.csr_matrix:
    mat = sp.coo_matrix((data, (rows, cols)), shape=shape, dtype=np.float32)
    return mat.tocsr()


def triplet_h5_to_csr(path: Path, n_genes: int, n_cells: int) -> sp.csr_matrix:
    with h5py.File(path, "r") as handle:
        gene_i = np.asarray(handle["i"][0], dtype=np.int64) - 1
        cell_j = np.asarray(handle["j"][0], dtype=np.int64) - 1
        values = np.asarray(handle["v"][0], dtype=np.float32)
    if gene_i.min() < 0 or cell_j.min() < 0:
        raise SystemExit(f"{path.name}: negative index after 1-based conversion")
    if gene_i.max() >= n_genes or cell_j.max() >= n_cells:
        raise SystemExit(
            f"{path.name}: index exceeds tables "
            f"gene {gene_i.max()} vs {n_genes}, cell {cell_j.max()} vs {n_cells}"
        )
    return csr_from_coo(values, cell_j, gene_i, (n_cells, n_genes))


def mtx_to_cells_by_genes(path: Path, n_genes: int, n_cells: int) -> sp.csr_matrix:
    matrix = scipy.io.mmread(path).tocsr().astype(np.float32)
    if matrix.shape == (n_genes, n_cells):
        return matrix.T.tocsr()
    if matrix.shape == (n_cells, n_genes):
        return matrix
    raise SystemExit(f"{path.name}: shape {matrix.shape} != {(n_genes, n_cells)} or {(n_cells, n_genes)}")


def read_10x_dir(directory: Path) -> ad.AnnData:
    mtx = directory / "matrix.mtx"
    if not mtx.exists():
        mtx = directory / "matrix.mtx.gz"
    feat = directory / "features.tsv"
    if not feat.exists():
        feat = directory / "features.tsv.gz"
    if not feat.exists():
        feat = directory / "genes.tsv"
    if not feat.exists():
        feat = directory / "genes.tsv.gz"
    barcodes = directory / "barcodes.tsv"
    if not barcodes.exists():
        barcodes = directory / "barcodes.tsv.gz"
    genes = pd.read_csv(open_text(feat), sep="\t", header=None)
    bars = pd.read_csv(open_text(barcodes), sep="\t", header=None)
    symbols = unique_symbols(genes.iloc[:, 1].astype(str).tolist() if genes.shape[1] > 1 else genes.iloc[:, 0].astype(str).tolist())
    matrix = mtx_to_cells_by_genes(mtx, n_genes=len(symbols), n_cells=len(bars))
    obs = pd.DataFrame(index=bars.iloc[:, 0].astype(str).tolist())
    var = pd.DataFrame(index=symbols)
    return ad.AnnData(X=matrix, obs=obs, var=var)


def extract_tar_member(tar_path: Path, dest: Path) -> None:
    dest.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tar_path) as handle:
        handle.extractall(dest)


def ensure_10x_layouts(data: Path) -> None:
    root = data / "x10"
    root.mkdir(exist_ok=True)
    if not (root / "K" / "matrix.mtx").exists():
        extract_tar_member(data / "GSE180963_RAW.tar", root)
        for name in ("GSM5481386_K.tar.gz", "GSM5481387_KL.tar.gz"):
            with tarfile.open(root / name) as handle:
                handle.extractall(root)
    if not (root / "KL1_count" / "filtered_feature_bc_matrix" / "matrix.mtx.gz").exists():
        extract_tar_member(data / "GSE165641_RAW.tar", root)
        for name in ("GSM5047302_KL1_count.tar.gz", "GSM5047303_KL2_count.tar.gz"):
            with tarfile.open(root / name) as handle:
                handle.extractall(root)
    g149 = data / "GSE149813"
    if not (g149 / "GSM4513594_AM1_matrix.mtx.gz").exists():
        g149.mkdir(exist_ok=True)
        extract_tar_member(data / "GSE149813_RAW.tar", g149)
    counts = data / "extract" / "00_normalized_combined_counts.txt"
    if not counts.exists():
        (data / "extract").mkdir(exist_ok=True)
        with tarfile.open(data / "GSE136246_mouse_normalized_gene_counts.txt.gz") as handle:
            handle.extractall(data / "extract")


def dense_symbol_matrix(path: Path, sep: str, header_cells: bool) -> ad.AnnData:
    """Gene x cell text matrix -> cells x genes sparse. One pass."""
    with open_text(path) as handle:
        if header_cells:
            header = handle.readline().rstrip("\n")
            cells = [c.strip() for c in header.split(sep) if c.strip()]
        else:
            cells = None
        gene_names = []
        chunks_r = []
        chunks_c = []
        chunks_v = []
        n_cells = None
        for line in handle:
            if sep == " ":
                sym, rest = line.rstrip("\n").split(" ", 1)
                values = np.fromstring(rest, sep=" ", dtype=np.float32)
            else:
                sym, rest = line.rstrip("\n").split(sep, 1)
                values = np.fromstring(rest, sep=sep, dtype=np.float32)
            if n_cells is None:
                n_cells = int(values.size)
                if cells is not None and len(cells) != n_cells:
                    raise SystemExit(f"{path.name}: header cells {len(cells)} != row width {n_cells}")
            elif values.size != n_cells:
                raise SystemExit(f"{path.name}: ragged row for {sym}")
            nz = np.flatnonzero(values)
            gi = len(gene_names)
            gene_names.append(sym)
            if nz.size:
                chunks_r.append(np.full(nz.size, gi, dtype=np.int32))
                chunks_c.append(nz.astype(np.int32))
                chunks_v.append(values[nz])
    if cells is None:
        cells = [f"cell{i}" for i in range(n_cells)]
    # File is genes x cells. COO below is built gene-major then transposed.
    gene_by_cell = sp.coo_matrix(
        (np.concatenate(chunks_v), (np.concatenate(chunks_r), np.concatenate(chunks_c))),
        shape=(len(gene_names), n_cells),
        dtype=np.float32,
    ).T.tocsr()
    symbols = unique_symbols(gene_names)
    return ad.AnnData(X=gene_by_cell, obs=pd.DataFrame(index=cells), var=pd.DataFrame(index=symbols))


def load_gse154989(data: Path) -> ad.AnnData:
    genes = read_table(data / "GSE154989_mmLungPlate_fQC_geneTable.csv.gz")
    samples = read_table(data / "GSE154989_mmLungPlate_fQC_smpTable.csv.gz")
    annot = read_table(data / "GSE154989_mmLungPlate_fQC_dZ_annot_smpTable.csv.gz")
    symbols = unique_symbols(genes["geneSymbol"].astype(str).tolist())
    matrix = triplet_h5_to_csr(
        data / "GSE154989_mmLungPlate_fQC_dSp_rawCount.h5",
        n_genes=len(symbols),
        n_cells=len(samples),
    )
    obs = samples.merge(annot, on="sampleID", how="left")
    obs.index = obs["sampleID"].astype(str)
    mouse = obs["mouseID"].astype(str).str.replace(r"_T\d+$", "", regex=True)
    obs["mouse_id"] = mouse
    obs["library_id"] = obs["mouseID"].astype(str)
    obs["genotype_group"] = np.where(mouse.str.startswith("KP"), "KP", np.where(mouse.str.startswith("K_"), "K", "unparsed"))
    obs["genotype_geo"] = obs["genotype_group"]
    obs["cohort"] = obs["timesimple"].astype(str)
    obs["tissue"] = "lung_tumor_plate"
    obs["author_label"] = obs["clusterK12"].astype(str)
    adata = ad.AnnData(X=matrix, obs=obs, var=pd.DataFrame(index=symbols))
    adata.uns["value_kind"] = "author_rawCount_h5"
    adata.uns["qc_policy"] = "author_fQC_keep_all"
    adata.uns["design"] = "sorted_or_tumor_cell_plate"
    adata.uns["sweep"] = "ineligible_tumor_cell_plate"
    return adata


def load_gse154977(data: Path) -> ad.AnnData:
    genes = read_table(data / "GSE154977_mmLung10x_cis_geneTable.csv.gz")
    samples = read_table(data / "GSE154977_mmLung10x_cis_smpTable.csv.gz")
    annot = read_table(data / "GSE154977_mmLung10x_cis_dZ_annot_annot_smpTable.csv.gz")
    symbol_col = "geneSymbol" if "geneSymbol" in genes.columns else "geneID"
    symbols = unique_symbols(genes[symbol_col].astype(str).tolist())
    matrix = triplet_h5_to_csr(
        data / "GSE154977_mmLung10x_cis_dSp_rawCount.h5",
        n_genes=len(symbols),
        n_cells=len(samples),
    )
    obs = samples.merge(annot, on="sampleID", how="left")
    obs.index = obs["sampleID"].astype(str)
    lib = obs["sampleID"].astype(str).str.replace(r"_id-.*$", "", regex=True)
    mouse = lib.str.replace(r"_PT$", "", regex=True)
    obs["mouse_id"] = mouse
    obs["library_id"] = lib
    obs["genotype_group"] = "KP"
    obs["genotype_geo"] = "KP"
    obs["cohort"] = np.where(mouse.str.contains("Cis"), "Cis72", "ND")
    obs["tissue"] = "lung_tumor_10x"
    obs["author_label"] = obs["timecourse_pred_cluster"].astype(str)
    adata = ad.AnnData(X=matrix, obs=obs, var=pd.DataFrame(index=symbols))
    adata.uns["value_kind"] = "raw_counts"
    adata.uns["qc_policy"] = "droplet"
    adata.uns["design"] = "kp_tumor_10x"
    adata.uns["sweep"] = "ineligible_until_T_seen"
    return adata


def load_xtr(data: Path, prefix: str, geo: dict, design: str, sweep: str) -> ad.AnnData:
    feat_name = f"{prefix}_features.tsv.gz"
    bar_name = f"{prefix}_barcodes.tsv.gz"
    mtx_name = f"{prefix}_matrix.mtx.gz"
    genes = pd.read_csv(open_text(data / feat_name), sep="\t", header=None)
    bars = pd.read_csv(open_text(data / bar_name), sep="\t", header=None)
    symbols = unique_symbols(genes.iloc[:, 1].astype(str).tolist())
    matrix = mtx_to_cells_by_genes(data / mtx_name, n_genes=len(symbols), n_cells=len(bars))
    cell_ids = bars.iloc[:, 0].astype(str)
    mouse = cell_ids.str.split("_").str[0]
    unknown = sorted(set(mouse) - set(geo))
    if unknown:
        raise SystemExit(f"{prefix}: barcode prefixes not in GEO table: {unknown}")
    obs = pd.DataFrame(index=cell_ids.tolist())
    obs["mouse_id"] = mouse.to_numpy()
    obs["library_id"] = obs["mouse_id"]
    obs["genotype_group"] = obs["mouse_id"].map(lambda m: geo[m][0])
    obs["genotype_geo"] = obs["mouse_id"].map(lambda m: geo[m][1])
    obs["cohort"] = obs["mouse_id"].map(lambda m: geo[m][2] + "|" + geo[m][3])
    obs["tissue"] = "lung"
    obs["author_label"] = ""
    adata = ad.AnnData(X=matrix, obs=obs, var=pd.DataFrame(index=symbols))
    adata.uns["value_kind"] = "raw_counts"
    adata.uns["qc_policy"] = "droplet"
    adata.uns["design"] = design
    adata.uns["sweep"] = sweep
    return adata


def load_gse165641(data: Path) -> ad.AnnData:
    pieces = []
    for folder, mouse, geo in (
        ("KL1_count/filtered_feature_bc_matrix", "KL1", "KrasLSL-G12D/+ Lkb1fl/fl"),
        ("KL2_count/filtered_feature_bc_matrix", "KL2", "KrasLSL-G12D/+ Lkb2fl/fl"),
    ):
        part = read_10x_dir(data / "x10" / folder)
        part.obs["mouse_id"] = mouse
        part.obs["library_id"] = mouse
        part.obs["genotype_group"] = "KL"
        part.obs["genotype_geo"] = geo
        part.obs["cohort"] = "AdCre_10w"
        part.obs["tissue"] = "lung"
        part.obs["author_label"] = ""
        part.obs_names = [f"{mouse}_{bc}" for bc in part.obs_names]
        pieces.append(part)
    adata = ad.concat(pieces, join="outer", fill_value=0)
    adata.uns["value_kind"] = "raw_counts"
    adata.uns["qc_policy"] = "droplet"
    adata.uns["design"] = "mixed_lung"
    adata.uns["sweep"] = "descriptive_n2"
    return adata


def load_gse180963(data: Path) -> ad.AnnData:
    pieces = []
    meta = {
        "K": ("K", "KrasG12D/+", "lenti-Tomato_or_K"),
        "KL": ("KL", "KrasG12D/+Lkb1fl/fl", "lenti-Lkb1_or_KL"),
    }
    for folder, (group, geo, cohort) in meta.items():
        part = read_10x_dir(data / "x10" / folder)
        part.obs["mouse_id"] = folder
        part.obs["library_id"] = folder
        part.obs["genotype_group"] = group
        part.obs["genotype_geo"] = geo
        part.obs["cohort"] = cohort
        part.obs["tissue"] = "lung"
        part.obs["author_label"] = ""
        part.obs_names = [f"{folder}_{bc}" for bc in part.obs_names]
        pieces.append(part)
    adata = ad.concat(pieces, join="outer", fill_value=0)
    adata.uns["value_kind"] = "raw_counts"
    adata.uns["qc_policy"] = "droplet"
    adata.uns["design"] = "mixed_lung"
    adata.uns["sweep"] = "descriptive_n2"
    return adata


def load_gse149813(data: Path) -> ad.AnnData:
    root = data / "GSE149813"
    spec = [
        ("GSM4513594_AM1", "mouse1", "YFP_negative"),
        ("GSM4513595_AM2", "mouse1", "YFP_positive"),
        ("GSM4513596_AM4", "mouse2", "YFP_negative"),
        ("GSM4513597_AM5", "mouse2", "YFP_positive"),
    ]
    pieces = []
    for stem, mouse, yfp in spec:
        folder = root / stem
        folder.mkdir(exist_ok=True)
        for kind in ("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz"):
            src = root / f"{stem}_{kind}"
            dest = folder / kind
            if not dest.exists():
                dest.symlink_to(src)
        part = read_10x_dir(folder)
        part.obs["mouse_id"] = mouse
        part.obs["library_id"] = stem
        part.obs["genotype_group"] = "K"
        part.obs["genotype_geo"] = "Kras LSL-G12D/WT; Rosa26 LSL-YFP/LSL-YFP"
        part.obs["cohort"] = yfp
        part.obs["tissue"] = "sorted_lung_epithelium"
        part.obs["author_label"] = yfp
        part.obs_names = [f"{stem}_{bc}" for bc in part.obs_names]
        pieces.append(part)
    adata = ad.concat(pieces, join="inner")
    adata.uns["value_kind"] = "raw_counts"
    adata.uns["qc_policy"] = "droplet"
    adata.uns["design"] = "sorted_epithelium"
    adata.uns["sweep"] = "ineligible_sorted_epithelium"
    return adata


def load_gse127465(data: Path) -> ad.AnnData:
    genes = pd.read_csv(open_text(data / "GSE127465_gene_names_mouse_28205.tsv.gz"), header=None)
    meta = pd.read_csv(data / "GSE127465_mouse_cell_metadata_15939x12.tsv.gz", sep="\t")
    symbols = unique_symbols(genes.iloc[:, 0].astype(str).tolist())
    matrix = mtx_to_cells_by_genes(
        data / "GSE127465_mouse_counts_normalized_15939x28205.mtx.gz",
        n_genes=len(symbols),
        n_cells=len(meta),
    )
    obs = meta.copy()
    obs.index = [
        f"{rep}_{lib}_{bc}"
        for rep, lib, bc in zip(obs["Biological replicate"], obs["Library"], obs["Barcode"])
    ]
    obs["mouse_id"] = obs["Tumor or healthy"].astype(str) + "|" + obs["Biological replicate"].astype(str)
    obs["library_id"] = obs["Library"].astype(str)
    obs["genotype_group"] = obs["Tumor or healthy"].astype(str)
    obs["genotype_geo"] = "not in deposited mouse metadata columns"
    obs["cohort"] = obs["Tumor or healthy"].astype(str)
    obs["tissue"] = "lung_CD45"
    obs["author_label"] = obs["Major cell type"].astype(str)
    adata = ad.AnnData(X=matrix, obs=obs, var=pd.DataFrame(index=symbols))
    adata.uns["value_kind"] = "author_normalized"
    adata.uns["qc_policy"] = "author_normalized"
    adata.uns["design"] = "cd45_immune"
    adata.uns["sweep"] = "ineligible_cd45"
    return adata


def load_gse267321(data: Path) -> ad.AnnData:
    adata = dense_symbol_matrix(
        data / "GSE267321_Normalized_expression_matrix_02122026.csv.gz",
        sep=",",
        header_cells=True,
    )
    cells = pd.Index(adata.obs_names.astype(str))
    # K_BARCODE_LKR13.K.1
    prefix = cells.str.split("_").str[0]
    suffix = cells.str.replace(r"^[^_]+_[^_]+_", "", regex=True)
    if not suffix.str.startswith("LKR13.").all():
        raise SystemExit("GSE267321: unexpected barcode suffix")
    adata.obs["mouse_id"] = suffix.to_numpy()
    adata.obs["library_id"] = suffix.to_numpy()
    adata.obs["genotype_group"] = prefix.to_numpy()
    adata.obs["genotype_geo"] = prefix.map(
        {"K": "KrasG12D", "KK": "KrasG12D, KEAP1 knockout", "KLK": "KrasG12D, KEAP1/LKB1 knockout"}
    ).to_numpy()
    adata.obs["cohort"] = "subcutaneous_no_treatment"
    adata.obs["tissue"] = "subcutaneous"
    adata.obs["author_label"] = ""
    adata.uns["value_kind"] = "author_normalized"
    adata.uns["qc_policy"] = "author_normalized"
    adata.uns["design"] = "subcutaneous_nonmalignant"
    adata.uns["sweep"] = "ineligible_subcutaneous_nonmalignant"
    return adata


def load_gse136246(data: Path) -> ad.AnnData:
    adata = dense_symbol_matrix(
        data / "extract" / "00_normalized_combined_counts.txt",
        sep=" ",
        header_cells=True,
    )
    cells = pd.Index(adata.obs_names.astype(str))
    # bcXXXX_KP_NegCtrl_VHL1_HD -> library token after the barcode
    library = cells.str.replace(r"^bc[^_]+_", "", regex=True)
    if library.str.contains(" ").any() or (library == "").any():
        raise SystemExit("GSE136246: failed to parse library token")
    groupings = pd.read_csv(data / "GSE136246_mouse_cell_groupings.csv.gz", header=None)
    row = groupings.iloc[0].astype(str).tolist()
    if row and row[0] == "graph_partition":
        row = row[1:]
    if len(row) != adata.n_obs:
        raise SystemExit(f"GSE136246 groupings {len(row)} != cells {adata.n_obs}")
    adata.obs["mouse_id"] = library.to_numpy()
    adata.obs["library_id"] = library.to_numpy()
    adata.obs["genotype_group"] = "KP"
    adata.obs["genotype_geo"] = "LSL-K-RasG12D, p53 flox"
    adata.obs["cohort"] = library.to_numpy()
    adata.obs["tissue"] = "lung"
    adata.obs["author_label"] = row
    adata.uns["value_kind"] = "author_normalized"
    adata.uns["qc_policy"] = "author_normalized"
    adata.uns["design"] = "mixed_lung"
    adata.uns["sweep"] = "eligible_if_both_compartments"
    return adata


def symbol_map(var_names: pd.Index) -> dict[str, int]:
    mapping = {}
    for i, name in enumerate(var_names.astype(str)):
        base = name.split("__dup")[0]
        mapping.setdefault(base, i)
    return mapping


def first_index(mapping: dict[str, int], aliases: list[str]) -> tuple[str, int] | None:
    for name in aliases:
        if name in mapping:
            return name, mapping[name]
    return None


def present_indices(mapping: dict[str, int], names: list[str]) -> tuple[list[str], list[int], list[str]]:
    found_names, found_idx, missing = [], [], []
    for name in names:
        if name in mapping:
            found_names.append(name)
            found_idx.append(mapping[name])
        else:
            missing.append(name)
    return found_names, found_idx, missing


def column(matrix: sp.spmatrix, idx: int) -> np.ndarray:
    return np.asarray(matrix[:, idx].todense()).ravel()


def row_means(matrix: sp.spmatrix, idxs: list[int]) -> np.ndarray:
    if not idxs:
        return np.full(matrix.shape[0], np.nan, dtype=np.float32)
    return np.asarray(matrix[:, idxs].mean(axis=1)).ravel().astype(np.float32)


def choose_scale(matrix: sp.csr_matrix, declared: str) -> str:
    """Decide whether to log-normalize. Declared kind is checked against the matrix."""
    rng = np.random.default_rng(0)
    n = matrix.shape[0]
    take = rng.choice(n, size=min(n, 2000), replace=False)
    sample = matrix[take]
    sums = np.asarray(sample.sum(axis=1)).ravel()
    data = sample.data
    if data.size == 0:
        raise SystemExit("empty matrix")
    frac_int = float(np.mean(np.isclose(data, np.round(data))))
    p99 = float(np.quantile(data, 0.99))
    med = float(np.median(sums))
    if declared == "raw_counts" or declared == "author_rawCount_h5":
        if p99 < 5 and frac_int < 0.5:
            raise SystemExit(f"declared counts but p99={p99} frac_int={frac_int}")
        return "log1p_cp10k"
    if declared == "author_normalized":
        if 5_000 <= med <= 20_000 and p99 > 20:
            return "log1p_of_author_linear"
        return "author_values_as_stored"
    raise SystemExit(f"unknown declared kind {declared}")


def to_lognorm(matrix: sp.csr_matrix, mode: str) -> sp.csr_matrix:
    if mode == "author_values_as_stored":
        out = matrix.astype(np.float32)
        out.eliminate_zeros()
        return out
    sums = np.asarray(matrix.sum(axis=1)).ravel()
    scale = np.zeros_like(sums, dtype=np.float32)
    ok = sums > 0
    scale[ok] = np.float32(1e4) / sums[ok].astype(np.float32)
    scaled = matrix.multiply(scale[:, None]).tocsr()
    scaled.data = np.log1p(scaled.data)
    return scaled


def apply_qc(adata: ad.AnnData, lognorm: sp.csr_matrix) -> np.ndarray:
    counts = adata.X.tocsr()
    n_genes = np.asarray((counts > 0).sum(axis=1)).ravel()
    mapping = symbol_map(adata.var_names)
    mito_idx = [i for name, i in mapping.items() if name.startswith("mt-")]
    if mito_idx and adata.uns["qc_policy"] == "droplet":
        mito_sum = np.asarray(counts[:, mito_idx].sum(axis=1)).ravel()
        total = np.asarray(counts.sum(axis=1)).ravel()
        pct = np.divide(mito_sum, total, out=np.zeros_like(mito_sum), where=total > 0) * 100
    elif mito_idx:
        mito_sum = np.asarray(lognorm[:, mito_idx].sum(axis=1)).ravel()
        total = np.asarray(lognorm.sum(axis=1)).ravel()
        pct = np.divide(mito_sum, total, out=np.zeros_like(mito_sum, dtype=np.float64), where=total > 0) * 100
    else:
        pct = np.full(adata.n_obs, np.nan)
    adata.obs["n_genes"] = n_genes
    adata.obs["pct_mito"] = pct
    if adata.uns["qc_policy"] == "droplet":
        keep = (n_genes >= MIN_GENES_DROPLET) & (pct <= MAX_PCT_MITO)
    elif adata.uns["qc_policy"] == "author_normalized":
        keep = n_genes >= MIN_GENES_DROPLET
    elif adata.uns["qc_policy"] == "author_fQC_keep_all":
        keep = np.ones(adata.n_obs, dtype=bool)
    else:
        raise SystemExit(adata.uns["qc_policy"])
    return keep


def marker_labels(lognorm: sp.csr_matrix, mapping: dict[str, int]) -> tuple[np.ndarray, dict]:
    def vec(name: str) -> np.ndarray:
        hit = first_index(mapping, [name])
        if hit is None:
            return np.zeros(lognorm.shape[0], dtype=np.float32)
        return column(lognorm, hit[1]).astype(np.float32)

    epi_names, epi_idx, epi_missing = present_indices(mapping, EPI_GENES)
    t_names, t_idx, t_missing = present_indices(mapping, T_GENES)
    nk_names, nk_idx, nk_missing = present_indices(mapping, NK_GENES)
    epi_score = row_means(lognorm, epi_idx)
    t_score = row_means(lognorm, t_idx)
    nk_score = row_means(lognorm, nk_idx)
    epi_pos = np.zeros(lognorm.shape[0], dtype=np.int16)
    for name in ("Epcam", "Krt8", "Krt18", "Sftpc", "Sftpb", "Nkx2-1"):
        if name in mapping:
            epi_pos += (column(lognorm, mapping[name]) > 0).astype(np.int16)
    cd3 = (vec("Cd3d") > 0) | (vec("Cd3e") > 0)
    # Sftpc alone is not a doublet call: lung T cells often carry ambient Sftpc.
    # A doublet needs Cd3 plus Epcam and a second epithelial marker.
    epcam = vec("Epcam") > 0
    labels = np.full(lognorm.shape[0], "other", dtype=object)
    doublet = cd3 & epcam & (epi_pos >= 2) & (epi_score > 0) & (t_score > 0)
    labels[doublet] = "doublet"
    is_t = (~doublet) & cd3 & (t_score >= epi_score)
    labels[is_t] = "T"
    nk_pos = np.zeros(lognorm.shape[0], dtype=bool)
    for name in ("Ncr1", "Klrb1c"):
        if name in mapping:
            nk_pos |= column(lognorm, mapping[name]) > 0
    is_nk = (~doublet) & (~is_t) & nk_pos & (~cd3) & (nk_score > epi_score)
    labels[is_nk] = "NK"
    is_epi = (~doublet) & (~is_t) & (~is_nk) & (epi_pos >= 2) & (~cd3) & (epi_score > t_score)
    labels[is_epi] = "epi"
    info = {
        "epi_genes": epi_names,
        "epi_missing": epi_missing,
        "t_genes": t_names,
        "t_missing": t_missing,
        "nk_genes": nk_names,
        "nk_missing": nk_missing,
        "epi_score": epi_score,
        "t_score": t_score,
    }
    return labels, info


def program_scores(lognorm: sp.csr_matrix, mapping: dict[str, int]) -> dict[str, tuple[np.ndarray, list[str]]]:
    out = {}
    for key, names in (
        ("Cldn4", ["Cldn4"]),
        ("Tacstd2", ["Tacstd2"]),
        ("IFN", IFN_GENES),
        ("NHEJ", NHEJ_CORE),
    ):
        found, idxs, _missing = present_indices(mapping, names)
        out[key] = (row_means(lognorm, idxs), found)
    paxx = first_index(mapping, PAXX_ALIASES)
    nhej_names = list(out["NHEJ"][1])
    nhej_idx = [mapping[n] for n in nhej_names]
    if paxx is not None:
        nhej_names.append(paxx[0])
        nhej_idx.append(paxx[1])
    out["NHEJ"] = (row_means(lognorm, nhej_idx), nhej_names)
    sting_names = []
    sting_idx = []
    for aliases in (STING_ALIASES, CGAS_ALIASES):
        hit = first_index(mapping, aliases)
        if hit is not None:
            sting_names.append(hit[0])
            sting_idx.append(hit[1])
    out["STING"] = (row_means(lognorm, sting_idx), sting_names)
    if not out["Cldn4"][1]:
        raise SystemExit("Cldn4 is absent")
    if not out["Tacstd2"][1]:
        raise SystemExit("Tacstd2 is absent")
    if len(out["NHEJ"][1]) < 5:
        raise SystemExit(f"NHEJ panel too small: {out['NHEJ'][1]}")
    if len(out["STING"][1]) < 1:
        raise SystemExit("STING panel empty")
    if len(out["IFN"][1]) < 8:
        raise SystemExit(f"IFN panel too small: {out['IFN'][1]}")
    return out


def final_labels(marker: np.ndarray, author: pd.Series) -> np.ndarray:
    author = author.fillna("").astype(str).to_numpy()
    # Cluster ids and library tokens are not cell-type calls. Use the author
    # field only when it actually says T cells or NK cells (GSE127465).
    if not np.any(np.isin(author, ["T cells", "NK cells"])):
        return marker
    labels = np.full(len(marker), "other", dtype=object)
    labels[author == "T cells"] = "T"
    labels[author == "NK cells"] = "NK"
    labels[marker == "epi"] = "epi"
    both = (author == "T cells") & (marker == "epi")
    labels[both] = "doublet"
    both_nk = (author == "NK cells") & (marker == "epi")
    labels[both_nk] = "doublet"
    return labels


def process(adata: ad.AnnData, dataset: str) -> tuple[ad.AnnData, dict]:
    if not sp.isspmatrix_csr(adata.X):
        adata.X = sp.csr_matrix(adata.X)
    mode = choose_scale(adata.X, adata.uns["value_kind"])
    lognorm = to_lognorm(adata.X, mode)
    keep = apply_qc(adata, lognorm)
    n_before = int(adata.n_obs)
    adata = adata[keep].copy()
    lognorm = lognorm[keep]
    mapping = symbol_map(adata.var_names)
    marker, info = marker_labels(lognorm, mapping)
    labels = final_labels(marker, adata.obs["author_label"])
    scores = program_scores(lognorm, mapping)
    adata.obs["marker_label"] = marker
    adata.obs["label"] = labels
    adata.obs["dataset"] = dataset
    adata.obs["score_epi"] = info["epi_score"]
    adata.obs["score_t"] = info["t_score"]
    for key, (values, _names) in scores.items():
        adata.obs[f"score_{key}"] = values
    adata.layers["lognorm"] = lognorm
    if mode != "author_values_as_stored":
        adata.layers["counts"] = adata.X.copy()
        adata.X = lognorm
    else:
        adata.X = lognorm
    meta = {
        "dataset": dataset,
        "scale_mode": mode,
        "n_before_qc": n_before,
        "n_after_qc": int(adata.n_obs),
        "design": adata.uns["design"],
        "sweep_rule": adata.uns["sweep"],
        "genes": {key: names for key, (_v, names) in scores.items()},
        "epi_genes": info["epi_genes"],
        "epi_missing": info["epi_missing"],
        "t_genes": info["t_genes"],
        "t_missing": info["t_missing"],
    }
    # Label check when both compartments are large enough to compare.
    check_label_separation(adata, dataset)
    return adata, meta


def check_label_separation(adata: ad.AnnData, dataset: str) -> None:
    epi = adata.obs["label"].to_numpy() == "epi"
    tcell = adata.obs["label"].to_numpy() == "T"
    if epi.sum() < MIN_EPI or tcell.sum() < MIN_EPI:
        return
    mapping = symbol_map(adata.var_names)
    lognorm = adata.layers["lognorm"]
    def mean(mask, gene):
        if gene not in mapping:
            return np.nan
        return float(column(lognorm, mapping[gene])[mask].mean())
    epcam_epi, epcam_t = mean(epi, "Epcam"), mean(tcell, "Epcam")
    krt_epi, krt_t = mean(epi, "Krt8"), mean(tcell, "Krt8")
    cd3_epi, cd3_t = mean(epi, "Cd3e"), mean(tcell, "Cd3e")
    epi_marker_ok = (epcam_epi > epcam_t) or (krt_epi > krt_t)
    t_marker_ok = cd3_t > cd3_epi
    if not (epi_marker_ok and t_marker_ok):
        raise SystemExit(
            f"{dataset}: label check failed "
            f"Epcam epi/T {epcam_epi}/{epcam_t} Krt8 {krt_epi}/{krt_t} Cd3e {cd3_epi}/{cd3_t}"
        )


def mouse_table(adata: ad.AnnData, meta: dict) -> pd.DataFrame:
    rows = []
    obs = adata.obs
    for mouse, sub in obs.groupby("mouse_id", sort=True):
        epi = sub[sub["label"] == "epi"]
        n_epi = int(len(epi))
        n_t = int((sub["label"] == "T").sum())
        n_nk = int((sub["label"] == "NK").sum())
        n_cells = int(len(sub))
        rec = {
            "dataset": meta["dataset"],
            "mouse_id": mouse,
            "genotype_group": sub["genotype_group"].iloc[0],
            "genotype_geo": sub["genotype_geo"].iloc[0],
            "cohort": sub["cohort"].iloc[0],
            "tissue": sub["tissue"].iloc[0],
            "n_libraries": int(sub["library_id"].nunique()),
            "n_cells": n_cells,
            "n_epi": n_epi,
            "n_T": n_t,
            "n_NK": n_nk,
            "n_doublet": int((sub["label"] == "doublet").sum()),
            "n_other": int((sub["label"] == "other").sum()),
            "frac_T": n_t / n_cells,
            "frac_T_of_epi_plus_T": (n_t / (n_t + n_epi)) if (n_t + n_epi) else np.nan,
            "epi_score_ok": n_epi >= MIN_EPI,
        }
        for key in ("Cldn4", "Tacstd2", "NHEJ", "STING", "IFN"):
            col = f"score_{key}"
            if n_epi >= MIN_EPI:
                rec[f"epi_mean_{key}"] = float(epi[col].mean())
                rec[f"epi_frac_pos_{key}"] = float((epi[col] > 0).mean()) if key in ("Cldn4", "Tacstd2") else np.nan
            else:
                rec[f"epi_mean_{key}"] = np.nan
                rec[f"epi_frac_pos_{key}"] = np.nan
        if n_t >= MIN_EPI:
            rec["T_mean_IFN"] = float(sub.loc[sub["label"] == "T", "score_IFN"].mean())
        else:
            rec["T_mean_IFN"] = np.nan
        # YFP-positive epithelium, only when that cohort exists.
        if (sub["cohort"] == "YFP_positive").any():
            yfp = sub[(sub["cohort"] == "YFP_positive") & (sub["label"] == "epi")]
            rec["n_epi_YFPpos"] = int(len(yfp))
            if len(yfp) >= MIN_EPI:
                rec["epi_YFPpos_mean_Cldn4"] = float(yfp["score_Cldn4"].mean())
                rec["epi_YFPpos_mean_Tacstd2"] = float(yfp["score_Tacstd2"].mean())
            else:
                rec["epi_YFPpos_mean_Cldn4"] = np.nan
                rec["epi_YFPpos_mean_Tacstd2"] = np.nan
        rows.append(rec)
    out = pd.DataFrame(rows)
    for key, genes in meta["genes"].items():
        out[f"genes_{key}"] = ",".join(genes)
    out["scale_mode"] = meta["scale_mode"]
    out["design"] = meta["design"]
    return out


def sweep_dataset(dataset: str, adata: ad.AnnData, mice: pd.DataFrame, rule: str) -> list[dict]:
    epi = adata.obs[adata.obs["label"] == "epi"]
    usable = mice[mice["epi_score_ok"]].copy()
    base = {
        "dataset": dataset,
        "design": adata.uns["design"],
        "sweep_rule": rule,
        "n_mice_with_epi": int(len(usable)),
        "n_epi_cells": int(len(epi)),
    }
    rows = []
    if len(epi) == 0 or len(usable) == 0:
        rows.append({**base, "threshold_name": "none", "threshold_value": np.nan, "in_sweep": "no", "rho_frac_T": np.nan, "p_frac_T": np.nan, "n_corr": 0, "reason": "no epithelium"})
        return rows
    values = epi["score_Cldn4"].to_numpy()
    # Per-mouse Cldn4 vector aligned to usable mice, built once per threshold.
    thresholds = [("gt0", 0.0, "yes")]
    for q in SWEEP_QS:
        thresholds.append((f"q{q:.2f}", float(np.quantile(values, q)), "yes"))
    # A single stray T cell can create a rank correlation. Require a real T compartment.
    t_ok = int((usable["n_T"] > 0).sum()) >= 3 and int(usable["n_T"].sum()) >= 20
    eligible = rule == "eligible_if_both_compartments" and len(usable) >= 4 and t_ok
    for name, thr, in_sweep in thresholds:
        high_frac = []
        for mouse in usable["mouse_id"]:
            block = epi[epi["mouse_id"] == mouse]
            high_frac.append(float((block["score_Cldn4"] > thr).mean()))
        usable = usable.copy()
        usable["cldn4_high_frac"] = high_frac
        sp_t = spearman(usable["cldn4_high_frac"], usable["frac_T"])
        sp_ratio = spearman(usable["cldn4_high_frac"], usable["frac_T_of_epi_plus_T"])
        if len(usable) >= 4 and not t_ok:
            sp_t = {"n": int(len(usable)), "rho": np.nan, "p": np.nan, "reason": "T compartment too thin"}
            sp_ratio = sp_t
        rows.append(
            {
                **base,
                "threshold_name": name,
                "threshold_value": thr,
                "primary_prespecified": name == f"q{PRIMARY_Q:.2f}",
                "in_sweep": in_sweep,
                "eligible_for_strongest": bool(eligible and in_sweep == "yes"),
                "rho_frac_T": sp_t["rho"],
                "p_frac_T": sp_t["p"],
                "rho_frac_T_of_epi_plus_T": sp_ratio["rho"],
                "p_frac_T_of_epi_plus_T": sp_ratio["p"],
                "n_corr": sp_t["n"],
                "reason": sp_t["reason"],
                "mean_high_frac": float(np.mean(high_frac)),
                "mean_frac_T": float(usable["frac_T"].mean()),
            }
        )
    # Continuous reference, not part of the min search.
    sp_mean = spearman(usable["epi_mean_Cldn4"], usable["frac_T"])
    if len(usable) >= 4 and not t_ok:
        sp_mean = {"n": int(len(usable)), "rho": np.nan, "p": np.nan, "reason": "T compartment too thin"}
    rows.append(
        {
            **base,
            "threshold_name": "mouse_mean",
            "threshold_value": np.nan,
            "primary_prespecified": False,
            "in_sweep": "no",
            "eligible_for_strongest": False,
            "rho_frac_T": sp_mean["rho"],
            "p_frac_T": sp_mean["p"],
            "rho_frac_T_of_epi_plus_T": np.nan,
            "p_frac_T_of_epi_plus_T": np.nan,
            "n_corr": sp_mean["n"],
            "reason": sp_mean["reason"] or "continuous mean, not a threshold",
            "mean_high_frac": np.nan,
            "mean_frac_T": float(usable["frac_T"].mean()),
        }
    )
    return rows


def spearman(x, y) -> dict:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    if len(x) < 4:
        return {"n": int(len(x)), "rho": np.nan, "p": np.nan, "reason": "n<4"}
    if np.unique(x).size < 2 or np.unique(y).size < 2:
        return {"n": int(len(x)), "rho": np.nan, "p": np.nan, "reason": "zero variance"}
    rho, p = stats.spearmanr(x, y)
    return {"n": int(len(x)), "rho": float(rho), "p": float(p), "reason": ""}


def genotype_rows(mice: pd.DataFrame) -> list[dict]:
    rows = []
    dataset = mice["dataset"].iloc[0]
    if dataset == "GSE154989":
        groups = [("K", "KP")]
        subsets = [("all_stages", mice)]
        for stage, sub in mice.groupby("cohort"):
            subsets.append((str(stage), sub))
    elif dataset == "GSE267321":
        groups = [("K", "KK"), ("K", "KLK"), ("KK", "KLK")]
        subsets = [("all", mice)]
    elif dataset == "GSE180963":
        groups = [("K", "KL")]
        subsets = [("all", mice)]
    elif dataset == "GSE179501":
        groups = []
        subsets = []
        # cohort string is "Restored|Female" or "Non-Restored|Female"
        restored = mice[mice["cohort"].str.startswith("Restored")]
        non = mice[mice["cohort"].str.startswith("Non-Restored")]
        if len(restored) and len(non):
            groups = [("Non-Restored", "Restored")]
            # Rebuild a temporary frame with genotype_group overwritten by cohort arm.
            tmp = mice.copy()
            tmp["genotype_group"] = np.where(tmp["cohort"].str.startswith("Restored"), "Restored", "Non-Restored")
            subsets = [("all_mice", tmp)]
    else:
        return rows
    scores = ["epi_mean_Cldn4", "epi_mean_Tacstd2", "epi_mean_NHEJ", "epi_mean_STING", "epi_mean_IFN"]
    for subset_name, frame in subsets:
        for left, right in groups:
            a = frame[frame["genotype_group"] == left]
            b = frame[frame["genotype_group"] == right]
            for score in scores:
                xa = a[score].dropna().to_numpy()
                xb = b[score].dropna().to_numpy()
                rec = {
                    "dataset": dataset,
                    "subset": subset_name,
                    "contrast": f"{right}_minus_{left}",
                    "score": score,
                    "n_left": int(len(xa)),
                    "n_right": int(len(xb)),
                    "median_left": float(np.median(xa)) if len(xa) else np.nan,
                    "median_right": float(np.median(xb)) if len(xb) else np.nan,
                }
                if len(xa) >= 3 and len(xb) >= 3:
                    stat = stats.mannwhitneyu(xb, xa, alternative="two-sided")
                    rec["U"] = float(stat.statistic)
                    rec["p"] = float(stat.pvalue)
                else:
                    rec["U"] = np.nan
                    rec["p"] = np.nan
                rec["delta_median"] = rec["median_right"] - rec["median_left"]
                rows.append(rec)
    return rows


def label_summary(adata: ad.AnnData, dataset: str) -> pd.DataFrame:
    mapping = symbol_map(adata.var_names)
    lognorm = adata.layers["lognorm"]
    rows = []
    for label, sub_idx in adata.obs.groupby("label").groups.items():
        mask = np.asarray(adata.obs_names.isin(sub_idx))
        rec = {"dataset": dataset, "label": label, "n_cells": int(mask.sum())}
        for gene in ("Epcam", "Krt8", "Sftpc", "Cd3e", "Cd3d", "Ncr1", "Cldn4", "Tacstd2"):
            if gene in mapping:
                rec[f"mean_{gene}"] = float(column(lognorm, mapping[gene])[mask].mean())
            else:
                rec[f"mean_{gene}"] = np.nan
        rows.append(rec)
    return pd.DataFrame(rows)


def write_finding(path: Path, mice: pd.DataFrame, sweep: pd.DataFrame, geno: pd.DataFrame, qc: pd.DataFrame) -> None:
    eligible = sweep[sweep["eligible_for_strongest"] == True].copy()  # noqa: E712
    lines = []
    lines.append("# Public mouse K/KP/KL lung scRNA: epithelium vs T, with NHEJ/STING/IFN")
    lines.append("")
    lines.append("Processed GEO matrices only. Private 8-KL matrices were not used and were not merged.")
    lines.append("Mice are not pooled across accessions. Scores are means on the log-normalized matrix described in `tables/qc_scale.tsv`.")
    lines.append("Epithelial means require at least 20 epithelial cells in that mouse. T fraction is T cells / QC cells in that mouse.")
    lines.append("")
    lines.append("## Cell compartments")
    lines.append("")
    lines.append("Marker rule: T if Cd3d or Cd3e is detected and the T-marker mean is at least the epithelial-marker mean; epithelium if at least two of Epcam/Krt8/Krt18/Sftpc/Sftpb/Nkx2-1 are detected, Cd3 is absent, and the epithelial mean exceeds the T mean. A doublet requires Cd3, Epcam, and a second epithelial marker. Ambient Sftpc alone does not make a doublet. GSE127465 uses deposited Major cell type for T and NK, and the marker rule for epithelium.")
    lines.append("")
    counts = (
        mice.groupby("dataset")
        .agg(n_mice=("mouse_id", "nunique"), n_cells=("n_cells", "sum"), n_epi=("n_epi", "sum"), n_T=("n_T", "sum"))
        .reset_index()
    )
    lines.append(counts.to_markdown(index=False))
    lines.append("")
    lines.append("## Threshold sweep")
    lines.append("")
    lines.append("Inside each accession, epithelial Cldn4 is cut at 0 and at the pooled epithelial quantiles 0.50, 0.60, 0.70, 0.75, 0.80, 0.90, and 0.95. The mouse feature is the fraction of that mouse's epithelial cells above the cut. The outcome is the T-cell fraction. Spearman is reported only for n>=4 mice with an epithelial score. The pre-specified primary cut is the 75th percentile. The strongest grid point is the minimum Spearman among eligible cuts. Eligibility requires a mixed lung design (not CD45-sorted, not epithelium-sorted, not a tumor-cell plate, not the subcutaneous non-malignant series), at least four mice, and a non-zero T count. The minimum over the grid is a grid result, not a locked single test. When Cldn4 is zero in most epithelial cells, several nominal quantiles collapse to the same cut (Cldn4 > 0).")
    lines.append("")
    if eligible.empty or not np.isfinite(eligible["rho_frac_T"]).any():
        lines.append("No accession produced an eligible Spearman on this grid.")
    else:
        finite = eligible[np.isfinite(eligible["rho_frac_T"])].sort_values(
            ["rho_frac_T", "threshold_name", "dataset"]
        )
        best = finite.iloc[0]
        tied = finite[
            (finite["dataset"] == best["dataset"])
            & np.isclose(finite["rho_frac_T"], best["rho_frac_T"], atol=1e-8)
        ]
        lines.append(
            f"Strongest grid point: {best['dataset']} threshold {best['threshold_name']} "
            f"(value {best['threshold_value']:.6g}), n={int(best['n_corr'])}, "
            f"Spearman rho={best['rho_frac_T']:.4f}, two-sided p={best['p_frac_T']:.4g}, "
            f"mean T fraction {best['mean_frac_T']:.4f}."
        )
        lines.append(
            "Same rho in that accession for: "
            + ", ".join(
                f"{r.threshold_name}={r.threshold_value:.6g}" for r in tied.itertuples()
            )
            + "."
        )
        primary = sweep[(sweep["primary_prespecified"] == True) & (sweep["dataset"] == best["dataset"])]  # noqa: E712
        if len(primary) == 1:
            row = primary.iloc[0]
            lines.append(
                f"Primary q=0.75 in that accession: value {row['threshold_value']:.6g}, "
                f"rho={row['rho_frac_T']:.4f}, p={row['p_frac_T']:.4g}."
            )
        others = (
            finite[finite["dataset"] != best["dataset"]]
            .sort_values(["dataset", "threshold_name"])
            .groupby("dataset", as_index=False)
            .first()
        )
        for _, other in others.iterrows():
            lines.append(
                f"Other eligible accession {other['dataset']} {other['threshold_name']}: "
                f"n={int(other['n_corr'])}, rho={other['rho_frac_T']:.4f}, p={other['p_frac_T']:.4g}."
            )
        lines.append(
            "Read the grid minimum together with the other eligible rows. "
            "A single negative rho on a zero-inflated Cldn4 cut is not a stable Cldn4-high / T-low result."
        )
    lines.append("")
    show = sweep[sweep["in_sweep"] == "yes"][
        ["dataset", "threshold_name", "threshold_value", "n_corr", "rho_frac_T", "p_frac_T", "eligible_for_strongest", "reason"]
    ]
    lines.append(show.to_markdown(index=False))
    lines.append("")
    lines.append("## Genotype contrasts on epithelial means")
    lines.append("")
    lines.append("Mann-Whitney p is filled only when both sides have at least 3 mice. GSE154989 contrasts K vs KP on the plate tumor-cell series (no T-cell fraction). The all-stage contrast mixes stages that are not balanced across genotypes; within-stage rows are the ones to read, and none of those rows has n>=3 on both sides. GSE180963 is one K library and one KL library.")
    lines.append("")
    if len(geno):
        lines.append(geno.to_markdown(index=False))
    lines.append("")
    lines.append("## Design notes")
    lines.append("")
    lines.append("- GSE165641: two KL lung libraries, mixed cells. GSM5047303 characteristics literally say Lkb2fl/fl; the series title is Lkb1. Genotype group is KL for both. n=2, so the sweep is descriptive.")
    lines.append("- GSE180963: one KrasG12D/+ library and one KrasG12D/+;Lkb1fl/fl library. The GEO growth protocol also describes lenti targeting. n=2.")
    lines.append("- GSE154977: four KP 30-week 10x libraries (two ND, two Cis72). Sweep eligibility is revoked if T cells are absent after the marker rule.")
    lines.append("- GSE179502: six Lkb1-XTR mice, sorted neoplastic cells. Epithelial scores only.")
    lines.append("- GSE179501: four Lkb1-XTR mice, total viable cells (the unsorted sister series). This is the KL mixed-lung series with mouse n=4.")
    lines.append("- GSE154989: K and KP plate tumor cells. Mouse id strips the tumor suffix (`_T#`). T-cell sweep is not eligible.")
    lines.append("- GSE267321: LKR13 subcutaneous K / KK / KLK, deposited as non-malignant cells. Not lung, not eligible for the lung T-cell sweep.")
    lines.append("- GSE127465: mouse arm is CD45-positive cells from two healthy and two tumor-bearing lungs. Author T/NK labels are used. Epithelial Cldn4 is not an epithelial compartment here.")
    lines.append("- GSE136246: Laughney KP lung, RBC-depleted, seven library tokens parsed from cell names. Treatment is the library token (NegCtrl/VHL/PTC), not relabeled.")
    lines.append("- GSE149813: two K mice, YFP-sorted epithelium at 7 weeks. Extra YFP-positive epithelial means are in the mouse table.")
    lines.append("- Found and not scored: GSE277777 (KP tumor-cell objects; combined h5ad is 8.0 GB) and GSE319598 (KcP eGFP+ tumor cells; stroma is pooled across mice, so per-mouse T is not in the deposit).")
    lines.append("")
    lines.append("QC scale modes are in `tables/qc_scale.tsv`. Label means are in `tables/label_means.tsv`.")
    path.write_text("\n".join(lines) + "\n")


def df_to_markdown(frame: pd.DataFrame, index: bool = False) -> str:
    if frame.empty:
        return "(empty)"
    cols = list(frame.columns)
    header = "| " + " | ".join(cols) + " |"
    sep = "| " + " | ".join("---" for _ in cols) + " |"
    body = []
    for _, row in frame.iterrows():
        cells = []
        for col in cols:
            val = row[col]
            if isinstance(val, float):
                cells.append("" if not np.isfinite(val) else f"{val:.6g}")
            else:
                cells.append(str(val))
        body.append("| " + " | ".join(cells) + " |")
    return "\n".join([header, sep, *body])


pd.DataFrame.to_markdown = df_to_markdown  # type: ignore[method-assign]


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=Path, default=Path("/tmp/public_mouse_scrna"))
    parser.add_argument("--out", type=Path, default=Path("methods/public_mouse_klkp_epi_t"))
    args = parser.parse_args()
    data = args.data
    out = args.out
    tables = out / "tables"
    figures = out / "figures"
    objects = out / "objects"
    for folder in (tables, figures, objects):
        folder.mkdir(parents=True, exist_ok=True)
    ensure_10x_layouts(data)

    loaders = [
        ("GSE154989", load_gse154989),
        ("GSE154977", load_gse154977),
        ("GSE179501", lambda p: load_xtr(p, "GSE179501_XTR_scRNAseq", GEO_179501, "mixed_lung", "eligible_if_both_compartments")),
        ("GSE179502", lambda p: load_xtr(p, "GSE179502_XTR_sorted_scRNAseq", GEO_179502, "sorted_neoplastic", "ineligible_sorted_neoplastic")),
        ("GSE165641", load_gse165641),
        ("GSE180963", load_gse180963),
        ("GSE149813", load_gse149813),
        ("GSE127465", load_gse127465),
        ("GSE267321", load_gse267321),
        ("GSE136246", load_gse136246),
    ]
    mouse_frames = []
    sweep_rows = []
    geno_rows = []
    label_frames = []
    qc_rows = []
    for dataset, loader in loaders:
        print(f"=== {dataset} ===", flush=True)
        adata = loader(data)
        adata, meta = process(adata, dataset)
        # Revoke sweep if the marker rule finds no T cells.
        if adata.uns["sweep"] == "ineligible_until_T_seen":
            if int((adata.obs["label"] == "T").sum()) == 0:
                adata.uns["sweep"] = "ineligible_no_T_cells"
            else:
                adata.uns["sweep"] = "eligible_if_both_compartments"
        mice = mouse_table(adata, meta)
        mouse_frames.append(mice)
        sweep_rows.extend(sweep_dataset(dataset, adata, mice, adata.uns["sweep"]))
        geno_rows.extend(genotype_rows(mice))
        label_frames.append(label_summary(adata, dataset))
        qc_rows.append(
            {
                "dataset": dataset,
                "scale_mode": meta["scale_mode"],
                "n_before_qc": meta["n_before_qc"],
                "n_after_qc": meta["n_after_qc"],
                "design": meta["design"],
                "sweep_rule": adata.uns["sweep"],
                "genes_NHEJ": ",".join(meta["genes"]["NHEJ"]),
                "genes_STING": ",".join(meta["genes"]["STING"]),
                "genes_IFN": ",".join(meta["genes"]["IFN"]),
                "epi_missing": ",".join(meta["epi_missing"]),
                "t_missing": ",".join(meta["t_missing"]),
            }
        )
        # Slim object: lognorm of scored genes plus obs. Full gene matrix stays on disk as h5ad too.
        keep_genes = sorted(set(EPI_GENES + T_GENES + NK_GENES + NHEJ_CORE + PAXX_ALIASES + STING_ALIASES + CGAS_ALIASES + IFN_GENES + ["Cldn4", "Tacstd2"]))
        present = [g for g in keep_genes if g in set(adata.var_names.astype(str).str.split("__dup").str[0])]
        # var names may include __dup; select by base name first occurrence via symbol map
        mapping = symbol_map(adata.var_names)
        idxs = [mapping[g] for g in present if g in mapping]
        slim = adata[:, idxs].copy()
        slim.write_h5ad(objects / f"{dataset}.h5ad")
        meta_path = tables / f"{dataset}_cell_metadata.tsv.gz"
        adata.obs.to_csv(meta_path, sep="\t")
        del adata
        print(f"wrote {dataset} cells={meta['n_after_qc']} scale={meta['scale_mode']}", flush=True)

    mice_all = pd.concat(mouse_frames, ignore_index=True)
    sweep = pd.DataFrame(sweep_rows)
    geno = pd.DataFrame(geno_rows)
    labels = pd.concat(label_frames, ignore_index=True)
    qc = pd.DataFrame(qc_rows)
    mice_all.to_csv(tables / "mouse_level_scores.tsv", sep="\t", index=False)
    sweep.to_csv(tables / "threshold_sweep.tsv", sep="\t", index=False)
    geno.to_csv(tables / "genotype_contrasts.tsv", sep="\t", index=False)
    labels.to_csv(tables / "label_means.tsv", sep="\t", index=False)
    qc.to_csv(tables / "qc_scale.tsv", sep="\t", index=False)

    # Recompute the reported strongest rho from the mouse table so the sentence cannot drift.
    eligible = sweep[sweep["eligible_for_strongest"] == True].copy()  # noqa: E712
    check = {"strongest": None}
    if len(eligible) and np.isfinite(eligible["rho_frac_T"]).any():
        finite = eligible[np.isfinite(eligible["rho_frac_T"])].sort_values(["rho_frac_T", "threshold_name", "dataset"])
        best = finite.iloc[0]
        block = mice_all[(mice_all["dataset"] == best["dataset"]) & (mice_all["epi_score_ok"])]
        # Reconstruct high fraction from stored cell metadata.
        cells = pd.read_csv(tables / f"{best['dataset']}_cell_metadata.tsv.gz", sep="\t")
        epi = cells[cells["label"] == "epi"]
        thr = float(best["threshold_value"])
        xs, ys = [], []
        for mouse, sub in block.groupby("mouse_id"):
            ee = epi[epi["mouse_id"] == mouse]
            xs.append(float((ee["score_Cldn4"] > thr).mean()))
            ys.append(float(sub["frac_T"].iloc[0]))
        again = spearman(xs, ys)
        if not np.isclose(again["rho"], best["rho_frac_T"], atol=1e-8):
            raise SystemExit(f"sweep rho drift {again['rho']} vs {best['rho_frac_T']}")
        check["strongest"] = {
            "dataset": best["dataset"],
            "threshold_name": best["threshold_name"],
            "rho": again["rho"],
            "p": again["p"],
            "n": again["n"],
        }
    (tables / "sweep_recompute.json").write_text(json.dumps(check, indent=2) + "\n")
    write_finding(out / "FINDING.md", mice_all, sweep, geno, qc)
    plot_sweep(mice_all, sweep, figures / "cldn4_high_vs_T.png")
    print(json.dumps(check, indent=2))


def plot_sweep(mice: pd.DataFrame, sweep: pd.DataFrame, path: Path) -> None:
    eligible_datasets = sorted(sweep.loc[sweep["eligible_for_strongest"] == True, "dataset"].unique())  # noqa: E712
    if not len(eligible_datasets):
        return
    fig, axes = plt.subplots(1, len(eligible_datasets), figsize=(4.2 * len(eligible_datasets), 3.6), squeeze=False)
    for ax, dataset in zip(axes[0], eligible_datasets):
        row = sweep[(sweep["dataset"] == dataset) & (sweep["threshold_name"] == f"q{PRIMARY_Q:.2f}")]
        block = mice[(mice["dataset"] == dataset) & (mice["epi_score_ok"])]
        ax.scatter(block["epi_mean_Cldn4"], block["frac_T"], s=28)
        for _, mouse in block.iterrows():
            ax.annotate(str(mouse["mouse_id"])[:18], (mouse["epi_mean_Cldn4"], mouse["frac_T"]), fontsize=6)
        title = dataset
        if len(row):
            title += f"\nq0.75 rho={row.iloc[0]['rho_frac_T']:.3f} n={int(row.iloc[0]['n_corr'])}"
        ax.set_title(title, fontsize=9)
        ax.set_xlabel("epithelial mean Cldn4 (continuous)")
        ax.set_ylabel("T fraction of QC cells")
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


if __name__ == "__main__":
    main()
