#!/usr/bin/env python3
"""Triple-merge CLDN4-only AUCell on GSE131907 + GSE148071 + GSE205335 malignant cells.

ADDITIVE. Patient is the unit. No TACSTD2 gate / no dual-high.
A10 ELF3–CLDN4 is taken as given. pySCENIC is not run; AUCell uses public priors
on a fixed gene universe (regulons + background). Scores are computed
within each dataset, then pooled by within-dataset ranks.
"""
from __future__ import annotations

import argparse
import gzip
import json
import shutil
import tempfile
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse, stats

HERE = Path(__file__).resolve().parent
TABLES = HERE / "tables"
FIGS = HERE / "figures"
MIN_MAL = 20
AUC_THR = 0.05
TUMOR_131907 = {"tLung", "tL/B", "mLN", "mBrain", "PE"}
MAL_SUB_131907 = {"Malignant cells", "tS1", "tS2", "tS3"}

# CLDN4 is the predictor and is held out of TJ / ELF3 sets.
IFN_GENES = [
    "STAT1", "STAT2", "IRF1", "IRF2", "IRF7", "IRF9", "JAK2", "IFNGR1", "IFNGR2",
    "CXCL9", "CXCL10", "CXCL11", "IDO1", "GBP1", "GBP2", "GBP4", "GBP5",
    "IFI27", "IFI35", "IFI44", "IFI44L", "IFIT1", "IFIT2", "IFIT3",
    "IFITM1", "IFITM3", "ISG15", "ISG20", "MX1", "MX2",
    "OAS1", "OAS2", "OAS3", "OASL", "RSAD2", "SAMD9", "SAMD9L",
    "IFIH1", "DDX58", "SOCS1", "XAF1", "BST2", "LY6E", "EPSTI1",
    "CMPK2", "HERC5", "HERC6", "USP18", "TRIM21", "TRIM22",
    "PARP9", "PARP14", "APOL6", "BATF2", "SP110", "SP100", "NMI", "PML",
    "WARS1", "WARS",
]
MHC_GENES = [
    "NLRC5", "CIITA", "RFX5", "RFXANK", "RFXAP",
    "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G",
    "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1",
    "HLA-DMA", "HLA-DMB", "B2M", "TAP1", "TAP2", "TAPBP", "TAPBPL",
    "PSMB8", "PSMB9", "PSMB10", "PSME1", "PSME2", "ERAP1", "ERAP2",
    "CALR", "CANX", "PDIA3", "CD74",
]
TJ_GENES = [
    "ELF3", "GRHL1", "GRHL2", "GRHL3", "KLF4",
    "CLDN1", "CLDN3", "CLDN7", "CLDN8", "OCLN",
    "TJP1", "TJP2", "TJP3", "F11R", "JAM2", "JAM3",
    "CGN", "MARVELD2", "MARVELD3", "CRB3", "PARD3", "PARD6B", "PRKCI",
    "CDH1", "CDH3", "EPCAM", "DSP", "PKP3", "DSG2", "DSC2",
    "NECTIN1", "NECTIN2", "NECTIN4", "CEACAM1", "CEACAM5", "CEACAM6",
    "LSR", "MAGI1", "CGNL1", "LLGL2", "SCRIB", "RAB25", "RAB11A",
    "CTNNA1", "CTNNB1", "CTNND1",
]
ELF3_PRIOR = [
    "ANGPT1", "AR", "CCL20", "CCND3", "CLDN7", "COL2A1", "EGR1", "EHF",
    "ERBB2", "FOLH1", "GADD45GIP1", "HP", "IL1B", "ITGB6", "KLK3", "KRT4",
    "KRT8", "LYZ", "MED23", "MMP13", "MMP9", "NOS2", "NOTCH3", "POU5F1",
    "PRR9", "PTGS2", "SPMAP2", "SPRR1A", "SPRR1B", "SPRR2A", "SPRR3",
    "TGFB1", "TGFBR2", "TIMP3",
]
# Background ranks only. Not a tested regulon except HK_CTRL.
BACKGROUND = [
    "ACTB", "ACTG1", "GAPDH", "RPLP0", "RPL13A", "RPS18", "RPS6", "EEF1A1",
    "UBC", "YWHAZ", "HPRT1", "TBP", "PPIA", "GUSB", "SDHA", "PGK1", "LDHA",
    "ENO1", "PKM", "ALDOA", "TUBA1B", "TUBB", "VIM", "KRT8", "KRT18", "KRT19",
    "KRT7", "KRT5", "NKX2-1", "NAPSA", "SFTPB", "MUC1", "SOX2", "SOX4",
    "MYC", "JUN", "FOS", "TP53", "EGFR", "KRAS", "CDKN1A", "CCND1", "CDK4",
    "MDM2", "BCL2", "BAX", "HSP90AA1", "HSPA8", "ATP5F1A", "ATP5F1B",
    "COX4I1", "NDUFA4", "VDAC1", "SLC25A5", "EIF4A1", "EIF4G1", "PABPC1",
]
HK_CTRL = [
    "ACTB", "GAPDH", "RPLP0", "RPL13A", "RPS18", "EEF1A1", "UBC", "PPIA",
    "PGK1", "LDHA", "ENO1", "TUBA1B", "HSP90AA1", "HSPA8", "PABPC1",
]
ALWAYS = ["CLDN4", "TACSTD2", "EPCAM", "PTPRC", "ELF3"]

REGULONS = {
    "IFN_STAT1": IFN_GENES,
    "MHC_NLRC5": MHC_GENES,
    "TJ_barrier": TJ_GENES,          # CLDN4 held out
    "ELF3_CollecTRI": ELF3_PRIOR,    # A10 given; CLDN4 not in prior
    "HK_CTRL": HK_CTRL,
}


def wanted_genes() -> list[str]:
    s: set[str] = set(ALWAYS)
    for genes in REGULONS.values():
        s.update(genes)
    s.update(BACKGROUND)
    return sorted(s)


def fmt_p(value: float) -> str:
    if not np.isfinite(value):
        return "NA"
    return f"{value:.2e}" if value < 0.001 else f"{value:.3g}"


def fmt_num(value: float, digits: int = 3) -> str:
    if not np.isfinite(value):
        return "NA"
    return f"{value:+.{digits}f}"


def spearman(x, y) -> dict:
    xa = np.asarray(x, float)
    ya = np.asarray(y, float)
    m = np.isfinite(xa) & np.isfinite(ya)
    n = int(m.sum())
    if n < 4 or np.unique(xa[m]).size < 2 or np.unique(ya[m]).size < 2:
        return {"n": n, "rho": np.nan, "p": np.nan}
    rho, p = stats.spearmanr(xa[m], ya[m])
    return {"n": n, "rho": float(rho), "p": float(p)}


def q4_vs_q1(predictor, endpoint) -> dict:
    s = pd.DataFrame(
        {"c": np.asarray(predictor, float), "i": np.asarray(endpoint, float)}
    )
    s = s[np.isfinite(s["c"]) & np.isfinite(s["i"])].copy()
    rec = {
        "n": int(len(s)),
        "n_q1": np.nan,
        "n_q4": np.nan,
        "n_compared": np.nan,
        "median_q1": np.nan,
        "median_q4": np.nan,
        "delta_median": np.nan,
        "r_rb": np.nan,
        "p": np.nan,
        "thin": True,
    }
    if len(s) < 6:
        return rec
    ranks = s["c"].rank(method="average")
    try:
        qs = pd.qcut(ranks, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")
    except ValueError:
        return rec
    if qs.nunique() < 4:
        return rec
    q1 = s.loc[qs == "Q1", "i"]
    q4 = s.loc[qs == "Q4", "i"]
    n1, n4 = int(len(q1)), int(len(q4))
    if n1 < 2 or n4 < 2:
        return rec
    u, p = stats.mannwhitneyu(q4.to_numpy(), q1.to_numpy(), alternative="two-sided")
    rec.update(
        {
            "n_q1": n1,
            "n_q4": n4,
            "n_compared": n1 + n4,
            "median_q1": float(q1.median()),
            "median_q4": float(q4.median()),
            "delta_median": float(q4.median() - q1.median()),
            "r_rb": float((2.0 * float(u)) / (n4 * n1) - 1.0),
            "p": float(p),
            "thin": n1 < 3 or n4 < 3 or (n1 + n4) < 8,
        }
    )
    return rec


def assign_quartiles(values: pd.Series) -> pd.Series:
    ranks = values.rank(method="average")
    return pd.qcut(ranks, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")


def aucell(expr: np.ndarray, member: np.ndarray, auc_threshold: float = AUC_THR) -> np.ndarray:
    """Aibar-style AUCell on the supplied gene universe (not full transcriptome)."""
    n_cells, n_genes = expr.shape
    max_rank = max(int(np.ceil(auc_threshold * n_genes)), 1)
    idx = np.where(member)[0]
    n_set = int(idx.size)
    if n_set == 0:
        return np.full(n_cells, np.nan)
    order = np.argsort(-expr, axis=1, kind="mergesort")
    ranks = np.empty_like(order)
    row = np.arange(n_cells)[:, None]
    ranks[row, order] = np.arange(n_genes)[None, :]
    r = ranks[:, idx].astype(np.float64)
    contrib = np.clip(max_rank - r, 0, None)
    return contrib.sum(axis=1) / (n_set * max_rank)


def module_score(expr: np.ndarray, member: np.ndarray) -> np.ndarray:
    if not member.any():
        return np.full(expr.shape[0], np.nan)
    return expr[:, member].mean(axis=1)


def parse_series_matrix(path: Path) -> pd.DataFrame:
    fields: dict[str, list[str]] = {}
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if not line.startswith("!Sample_"):
                continue
            key, _, rest = line.rstrip("\n").partition("\t")
            key = key.removeprefix("!Sample_")
            values = [v.strip('"') for v in rest.split("\t")]
            if key == "characteristics_ch1":
                prefix, _, _ = values[0].partition(": ")
                key = prefix.strip().lower().replace(" ", "_")
                values = [v.split(": ", 1)[-1] if ": " in v else v for v in values]
            fields[key] = values
    n = len(fields["title"])
    return pd.DataFrame({k: v for k, v in fields.items() if len(v) == n})


def parse_geo_soft(path: Path) -> pd.DataFrame:
    records: list[dict[str, str]] = []
    current: dict[str, str] | None = None
    descriptions: list[str] = []
    titles: list[str] = []
    with gzip.open(path, "rt", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current is not None:
                    current["description"] = descriptions[0] if descriptions else ""
                    current["title"] = titles[0] if titles else ""
                    records.append(current)
                current = {"gsm": line.split(" = ", 1)[1]}
                descriptions, titles = [], []
            elif current is not None and line.startswith("!Sample_title = "):
                titles.append(line.split(" = ", 1)[1])
            elif current is not None and line.startswith("!Sample_description = "):
                descriptions.append(line.split(" = ", 1)[1])
            elif current is not None and line.startswith("!Sample_characteristics_ch1 = "):
                value = line.split(" = ", 1)[1]
                if ": " in value:
                    key, item = value.split(": ", 1)
                    current[key] = item
        if current is not None:
            current["description"] = descriptions[0] if descriptions else ""
            current["title"] = titles[0] if titles else ""
            records.append(current)
    return pd.DataFrame(records)


def stream_gse131907(matrix: Path, wanted: set[str]) -> tuple[list[str], dict[str, np.ndarray], np.ndarray]:
    found: dict[str, np.ndarray] = {}
    with gzip.open(matrix, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        cells = header[1:]
        n = len(cells)
        total = np.zeros(n, dtype=np.float64)
        print(f"[GSE131907] cells={n}", flush=True)
        n_genes = 0
        for line in handle:
            n_genes += 1
            gene, sep, rest = line.partition("\t")
            if not sep:
                continue
            gene = gene.split(".")[0]
            arr = np.fromstring(rest, sep="\t", dtype=np.float32)
            if arr.size != n:
                raise ValueError(f"{gene}: {arr.size} != {n}")
            total += arr
            if gene in wanted:
                found[gene] = arr.copy()
            if n_genes % 5000 == 0:
                print(f"  genes={n_genes} kept={len(found)}", flush=True)
    print(f"[GSE131907] done genes={n_genes} kept={len(found)}", flush=True)
    return cells, found, total


def load_gse131907(root: Path, wanted: set[str]) -> pd.DataFrame:
    ann = pd.read_csv(root / "GSE131907_Lung_Cancer_cell_annotation.txt.gz", sep="\t", dtype=str)
    meta = parse_series_matrix(root / "GSE131907_series_matrix.txt.gz")
    sample_meta = meta.rename(columns={"title": "Sample"})[
        ["Sample", "patient_id", "tumor_stage", "tissue_origin_abbrevation"]
    ]
    cells, expr, total = stream_gse131907(
        root / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz", wanted
    )
    per = ann.set_index("Index").reindex(cells).reset_index()
    if per["Sample"].isna().any():
        raise SystemExit("GSE131907 matrix barcodes do not match annotation Index")
    per = per.merge(sample_meta, on="Sample", how="left")
    per["total_umi"] = total
    tumor = per["Sample_Origin"].isin(TUMOR_131907)
    mal = (
        (per["Cell_type"] == "Epithelial cells")
        & per["Cell_subtype"].isin(MAL_SUB_131907)
        & tumor
    )
    per = per.loc[mal].copy()
    lib = np.maximum(per["total_umi"].to_numpy(), 1.0)
    out = pd.DataFrame(
        {
            "dataset": "GSE131907",
            "patient": per["patient_id"].to_numpy(),
            "sample": per["Sample"].to_numpy(),
            "site": per["Sample_Origin"].to_numpy(),
        }
    )
    for gene, arr in expr.items():
        umi = arr[mal.to_numpy()]
        out[gene] = np.log1p(umi / lib * 1e4).astype(np.float32)
        out[f"{gene}__umi"] = umi
    print(f"[GSE131907] malignant cells kept={len(out)} patients={out.patient.nunique()}", flush=True)
    return out


def load_gse205335_rds(path: Path, wanted: set[str]):
    with tempfile.TemporaryDirectory(prefix="gse205335-") as tmp:
        matrix_path = path
        if path.suffix == ".gz":
            matrix_path = Path(tmp) / path.stem
            print(f"[GSE205335] decompress {path.name}", flush=True)
            with gzip.open(path, "rb") as src, matrix_path.open("wb") as dest:
                shutil.copyfileobj(src, dest, 16 * 1024 * 1024)
        print("[GSE205335] read RDS", flush=True)
        import rdata

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message='Missing constructor for R class "dgCMatrix"')
            obj = rdata.read_rds(matrix_path)
    if not {"i", "p", "Dim", "Dimnames", "x"}.issubset(vars(obj)):
        raise TypeError("GSE205335 RDS is not Matrix::dgCMatrix")
    genes = np.asarray(obj.Dimnames[0], dtype=str)
    barcodes = np.asarray(obj.Dimnames[1], dtype=str)
    matrix = sparse.csc_matrix((obj.x, obj.i, obj.p), shape=tuple(obj.Dim), copy=False)
    library = np.asarray(matrix.sum(axis=0)).ravel()
    matrix = matrix.tocsr()
    name_to_row = {g: i for i, g in enumerate(genes)}
    extracted = {}
    for gene in sorted(wanted):
        row = name_to_row.get(gene)
        if row is None:
            continue
        extracted[gene] = np.asarray(matrix.getrow(row).toarray()).ravel()
    print(f"[GSE205335] extracted {len(extracted)}/{len(wanted)} genes", flush=True)
    return extracted, library, barcodes


def load_gse205335(root: Path, wanted: set[str]) -> pd.DataFrame:
    identities = pd.read_csv(root / "GSE205335_Lung_IO_CellIdentity.txt.gz", sep="\t")
    soft = parse_geo_soft(root / "GSE205335_family.soft.gz")
    soft["desc_key"] = soft["description"].str.replace("_", "-", regex=False)
    extracted, library, barcodes = load_gse205335_rds(
        root / "GSE205335_Lung_IO_UMI_matrix.rds.gz", wanted
    )
    cells = identities.set_index("barcode").loc[barcodes].reset_index()
    cells["total_umi"] = library
    cells["desc_key"] = cells["orig.ident"].str.replace(r"-[35]P$", "", regex=True)
    cells = cells.merge(
        soft[["desc_key", "patient", "tissue", "cancer subtype", "recist"]],
        on="desc_key",
        how="left",
        validate="many_to_one",
    )
    if cells["patient"].isna().any():
        missing = sorted(cells.loc[cells["patient"].isna(), "orig.ident"].unique())
        raise SystemExit(f"GSE205335 orig.ident unmatched: {missing}")
    normal = cells["tissue"].astype(str).str.contains("Normal", case=False, na=False)
    mal = cells["lineage.sub"].eq("Malignant cells") & ~normal
    keep = cells.loc[mal].copy()
    lib = np.maximum(keep["total_umi"].to_numpy(), 1.0)
    out = pd.DataFrame(
        {
            "dataset": "GSE205335",
            "patient": keep["patient"].to_numpy(),
            "sample": keep["orig.ident"].to_numpy(),
            "site": keep["tissue"].to_numpy(),
        }
    )
    mal_idx = np.flatnonzero(mal.to_numpy())
    for gene, arr in extracted.items():
        umi = arr[mal_idx]
        out[gene] = np.log1p(umi / lib * 1e4).astype(np.float32)
        out[f"{gene}__umi"] = umi
    print(f"[GSE205335] malignant cells kept={len(out)} patients={out.patient.nunique()}", flush=True)
    return out


def _as_str(arr) -> np.ndarray:
    out = []
    for x in arr:
        out.append(x.decode() if isinstance(x, bytes) else str(x))
    return np.array(out, dtype=object)


def load_tisch_genes(path: Path, wanted: list[str]) -> tuple[pd.DataFrame, list[str]]:
    import h5py

    with h5py.File(path, "r") as f:
        grp = None
        if all(k in f for k in ("data", "indices", "indptr")):
            grp = f
        else:
            for key in f.keys():
                g = f[key]
                if isinstance(g, h5py.Group) and all(k in g for k in ("data", "indices", "indptr")):
                    grp = g
                    break
        if grp is None:
            raise ValueError(f"no sparse matrix in {path}")

        def names(cands):
            for c in cands:
                node = grp.get(c) if hasattr(grp, "get") else None
                if node is None:
                    node = f.get(c)
                if node is None:
                    continue
                if isinstance(node, h5py.Dataset):
                    return _as_str(node[:])
                if isinstance(node, h5py.Group):
                    for sub in ("name", "id", "gene_names"):
                        if sub in node:
                            return _as_str(node[sub][:])
            return None

        genes = names(["features", "gene_names", "genes", "rownames"])
        barcodes = names(["barcodes", "cell_names", "colnames"])
        shape = tuple(int(x) for x in grp["shape"][:]) if "shape" in grp else None
        data = grp["data"][:]
        indices = grp["indices"][:]
        indptr = grp["indptr"][:]
    if shape is None:
        raise ValueError("missing shape")
    n0, n1 = int(shape[0]), int(shape[1])
    if genes is not None and barcodes is not None and len(genes) == n1 and len(barcodes) == n0:
        n_genes, n_cells = n1, n0
        csc = sparse.csc_matrix((data, indices, indptr), shape=(n_cells, n_genes))

        def col(i: int) -> np.ndarray:
            return np.asarray(csc[:, i].todense()).ravel()
    else:
        n_genes, n_cells = n0, n1
        csc = sparse.csc_matrix((data, indices, indptr), shape=(n_genes, n_cells))

        def col(i: int) -> np.ndarray:
            return np.asarray(csc.getrow(i).todense()).ravel()
    present = [g for g in wanted if g in set(genes)]
    gene_index = {g: int(np.where(genes == g)[0][0]) for g in present}
    table = {"barcode": barcodes}
    for g, i in gene_index.items():
        table[g] = col(i)
    return pd.DataFrame(table).set_index("barcode"), present


def load_gse148071(root: Path, wanted: set[str]) -> pd.DataFrame:
    meta = pd.read_csv(root / "NSCLC_GSE148071_CellMetainfo_table.tsv", sep="\t")
    expr, present = load_tisch_genes(root / "NSCLC_GSE148071_expression.h5", sorted(wanted))
    print(f"[GSE148071] TISCH genes present {len(present)}/{len(wanted)}", flush=True)
    meta = meta.copy()
    meta["barcode"] = meta["Cell"].astype(str)
    # TISCH barcodes may be Cell, or the token after @
    if not meta["barcode"].isin(expr.index).mean() > 0.8:
        alt = meta["Cell"].astype(str).str.split("@").str[-1]
        if alt.isin(expr.index).mean() > 0.8:
            meta["barcode"] = alt
        else:
            # try Cell as-is vs index containing Cell
            idx_map = {str(i): i for i in expr.index}
            hit = meta["Cell"].astype(str).map(idx_map)
            if hit.notna().mean() > 0.8:
                meta["barcode"] = hit
            else:
                raise SystemExit(
                    f"GSE148071 barcode overlap too low: "
                    f"{meta['barcode'].isin(expr.index).mean():.3f}"
                )
    mal = meta["Celltype (major-lineage)"].eq("Malignant")
    meta = meta.loc[mal].copy()
    joined = meta.set_index("barcode").join(expr, how="inner")
    if len(joined) < 100:
        raise SystemExit(f"GSE148071 malignant join too small: {len(joined)}")
    out = pd.DataFrame(
        {
            "dataset": "GSE148071",
            "patient": joined["Patient"].to_numpy(),
            "sample": joined["Sample"].to_numpy(),
            "site": joined["Tissue"].to_numpy(),
        }
    )
    for gene in present:
        out[gene] = joined[gene].to_numpy(dtype=np.float32)
    print(f"[GSE148071] malignant cells kept={len(out)} patients={out.patient.nunique()}", flush=True)
    return out


def score_cells(df: pd.DataFrame, universe: list[str]) -> pd.DataFrame:
    genes = [g for g in universe if g in df.columns and not g.endswith("__umi")]
    expr = df[genes].to_numpy(dtype=np.float32)
    expr = np.nan_to_num(expr, nan=0.0)
    scored = df[["dataset", "patient", "sample", "site"]].copy()
    if "CLDN4" in df.columns:
        scored["CLDN4"] = df["CLDN4"].to_numpy()
        umi_col = "CLDN4__umi"
        if umi_col in df.columns:
            scored["CLDN4_pos"] = (df[umi_col].to_numpy() > 0).astype(np.float32)
        else:
            scored["CLDN4_pos"] = (df["CLDN4"].to_numpy() > 0).astype(np.float32)
    for name, members in REGULONS.items():
        member = np.array([g in set(members) for g in genes])
        scored[f"auc_{name}"] = aucell(expr, member)
        scored[f"mod_{name}"] = module_score(expr, member)
    scored["n_universe"] = len(genes)
    return scored


def patient_table(cells: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (dataset, patient), sub in cells.groupby(["dataset", "patient"], sort=True):
        rec = {
            "dataset": dataset,
            "patient": patient,
            "n_malignant": int(len(sub)),
            "n_samples": int(sub["sample"].nunique()),
            "sites": ",".join(sorted(sub["site"].astype(str).unique())),
            "eligible": int(len(sub) >= MIN_MAL),
        }
        rec["CLDN4_mean"] = float(sub["CLDN4"].mean())
        rec["CLDN4_pctpos"] = float(100.0 * sub["CLDN4_pos"].mean())
        for name in REGULONS:
            rec[f"auc_{name}"] = float(sub[f"auc_{name}"].mean())
            rec[f"mod_{name}"] = float(sub[f"mod_{name}"].mean())
        rows.append(rec)
    out = pd.DataFrame(rows)
    out["cldn4_rank_in_dataset"] = out.groupby("dataset")["CLDN4_mean"].rank(method="average", pct=True)
    out["cldn4_pct_rank_in_dataset"] = out.groupby("dataset")["CLDN4_pctpos"].rank(method="average", pct=True)
    for name in REGULONS:
        out[f"auc_{name}_rank"] = out.groupby("dataset")[f"auc_{name}"].rank(method="average", pct=True)
        out[f"mod_{name}_rank"] = out.groupby("dataset")[f"mod_{name}"].rank(method="average", pct=True)
    # Quartiles inside each dataset on eligible patients only, then mapped back.
    out["q_mean"] = pd.NA
    out["q_pct"] = pd.NA
    for ds, sub in out.groupby("dataset"):
        elig = sub["eligible"] == 1
        if elig.sum() >= 6:
            qm = assign_quartiles(sub.loc[elig, "CLDN4_mean"]).astype(str)
            qp = assign_quartiles(sub.loc[elig, "CLDN4_pctpos"]).astype(str)
            out.loc[qm.index, "q_mean"] = qm
            out.loc[qp.index, "q_pct"] = qp
    return out


def contrast_rows(patients: pd.DataFrame) -> pd.DataFrame:
    elig = patients[patients["eligible"] == 1].copy()
    rows = []

    def add(scope, frame, predictor, pred_name, endpoint, end_name):
        sp = spearman(frame[predictor], frame[endpoint])
        q = q4_vs_q1(frame[predictor], frame[endpoint])
        rows.append(
            {
                "scope": scope,
                "n_patients": sp["n"],
                "predictor": pred_name,
                "regulon": end_name,
                "spearman_rho": sp["rho"],
                "spearman_p": sp["p"],
                "n_q1": q["n_q1"],
                "n_q4": q["n_q4"],
                "n_compared": q["n_compared"],
                "median_q1": q["median_q1"],
                "median_q4": q["median_q4"],
                "delta_median": q["delta_median"],
                "q4q1_r": q["r_rb"],
                "q4q1_p": q["p"],
                "thin": q["thin"],
            }
        )

    scopes = [("triple", elig)]
    for ds in ["GSE131907", "GSE148071", "GSE205335"]:
        scopes.append((ds, elig[elig.dataset == ds]))
    pair = elig[elig.dataset.isin(["GSE131907", "GSE205335"])]
    scopes.append(("pair_131907_205335_sensitivity", pair))

    endpoints = []
    for name in REGULONS:
        endpoints.append((f"auc_{name}", f"AUCell:{name}"))
        endpoints.append((f"mod_{name}", f"module:{name}"))
    # ranked-within-dataset versions for pooled scopes
    for scope, frame in scopes:
        if len(frame) < 4:
            continue
        pred_mean = "cldn4_rank_in_dataset" if scope.startswith(("triple", "pair")) else "CLDN4_mean"
        pred_pct = "cldn4_pct_rank_in_dataset" if scope.startswith(("triple", "pair")) else "CLDN4_pctpos"
        for end_col, end_name in endpoints:
            end_use = f"{end_col}_rank" if scope.startswith(("triple", "pair")) and f"{end_col}_rank" in frame else end_col
            add(scope, frame, pred_mean, "CLDN4_mean" if not scope.startswith(("triple", "pair")) else "CLDN4_mean_rank_in_dataset", end_use, end_name)
            add(scope, frame, pred_pct, "CLDN4_pctpos" if not scope.startswith(("triple", "pair")) else "CLDN4_pctpos_rank_in_dataset", end_use, end_name)
        # Stratified Q4 vs Q1: within-dataset quartiles, then pool tails
        if scope.startswith(("triple", "pair")):
            for qcol, pred_lab in (("q_mean", "CLDN4_mean_stratified_Q"), ("q_pct", "CLDN4_pctpos_stratified_Q")):
                tails = frame[frame[qcol].isin(["Q1", "Q4"])]
                for end_col, end_name in endpoints:
                    q1 = tails.loc[tails[qcol] == "Q1", end_col]
                    q4 = tails.loc[tails[qcol] == "Q4", end_col]
                    n1, n4 = int(len(q1)), int(len(q4))
                    if n1 < 2 or n4 < 2:
                        continue
                    u, p = stats.mannwhitneyu(q4.to_numpy(), q1.to_numpy(), alternative="two-sided")
                    rows.append(
                        {
                            "scope": scope,
                            "n_patients": int(len(frame)),
                            "predictor": pred_lab,
                            "regulon": end_name,
                            "spearman_rho": np.nan,
                            "spearman_p": np.nan,
                            "n_q1": n1,
                            "n_q4": n4,
                            "n_compared": n1 + n4,
                            "median_q1": float(q1.median()),
                            "median_q4": float(q4.median()),
                            "delta_median": float(q4.median() - q1.median()),
                            "q4q1_r": float((2.0 * float(u)) / (n4 * n1) - 1.0),
                            "q4q1_p": float(p),
                            "thin": n1 < 3 or n4 < 3,
                        }
                    )
    return pd.DataFrame(rows)


def honest_n(patients: pd.DataFrame, cells: pd.DataFrame) -> pd.DataFrame:
    rows = []
    rows.append({"item": "GEO patients GSE131907 (Kim 2020)", "n": 44, "note": "58 samples; not the test n"})
    rows.append({"item": "GEO patients GSE148071 (Wu 2021)", "n": 42, "note": "one sample each; not the test n"})
    rows.append({"item": "GEO patients GSE205335 (Hu 2022/2024)", "n": 26, "note": "33 samples; not the test n"})
    rows.append({"item": "GEO patients summed", "n": 44 + 42 + 26, "note": "do not write this as the test n"})
    for ds in ["GSE131907", "GSE148071", "GSE205335"]:
        sub = patients[patients.dataset == ds]
        elig = sub[sub.eligible == 1]
        rows.append({"item": f"{ds} patients with any scored malignant cell", "n": int(len(sub)), "note": ""})
        rows.append({"item": f"{ds} eligible (≥{MIN_MAL} malignant)", "n": int(len(elig)), "note": "test n for this dataset"})
        rows.append({"item": f"{ds} malignant cells scored", "n": int((cells.dataset == ds).sum()), "note": "cell n; not the test n"})
        if len(elig):
            q = elig["q_mean"].value_counts()
            rows.append({"item": f"{ds} Q1/Q4 tails (CLDN4 mean)", "n": int(q.get("Q1", 0) + q.get("Q4", 0)), "note": f"Q1={int(q.get('Q1', 0))} Q4={int(q.get('Q4', 0))}"})
    elig = patients[patients.eligible == 1]
    rows.append({"item": "triple eligible patients", "n": int(len(elig)), "note": "primary pooled n"})
    q = elig["q_mean"].value_counts()
    rows.append({"item": "triple stratified Q1+Q4 (CLDN4 mean)", "n": int(q.get("Q1", 0) + q.get("Q4", 0)), "note": f"Q1={int(q.get('Q1', 0))} Q4={int(q.get('Q4', 0))}; n_compared not 112"})
    pair = elig[elig.dataset.isin(["GSE131907", "GSE205335"])]
    rows.append({"item": "pair 131907+205335 eligible (sensitivity only)", "n": int(len(pair)), "note": "not the primary; other agent owns this pair"})
    return pd.DataFrame(rows)


def plot_boxes(patients: pd.DataFrame, out: Path) -> None:
    elig = patients[patients.eligible == 1].copy()
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.6), sharey=False)
    for ax, name, title in zip(
        axes,
        ["IFN_STAT1", "MHC_NLRC5", "TJ_barrier"],
        ["IFN (STAT1/IRF)", "MHC / antigen presentation", "TJ / barrier (CLDN4 held out)"],
    ):
        data = []
        labels = []
        for ds in ["GSE131907", "GSE148071", "GSE205335"]:
            sub = elig[elig.dataset == ds]
            for q, color in (("Q1", "#4C78A8"), ("Q4", "#F58518")):
                vals = sub.loc[sub["q_mean"] == q, f"auc_{name}"].dropna()
                data.append(vals.to_numpy())
                labels.append(f"{ds.replace('GSE','')}\n{q}")
        ax.boxplot(data, labels=labels, showfliers=False)
        ax.set_title(title, fontsize=9)
        ax.set_ylabel("AUCell (patient mean)" if ax is axes[0] else "")
        ax.tick_params(axis="x", labelsize=7)
    fig.suptitle("CLDN4-high vs low malignant AUCell — triple merge (within-dataset Q)", fontsize=11)
    fig.tight_layout()
    fig.savefig(out / "fig_q4q1_aucell.png", dpi=140)
    fig.savefig(out / "fig_q4q1_aucell.pdf")
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.6))
    for ax, name, title in zip(
        axes,
        ["IFN_STAT1", "MHC_NLRC5", "TJ_barrier"],
        ["IFN AUCell", "MHC AUCell", "TJ AUCell"],
    ):
        for ds, c in zip(["GSE131907", "GSE148071", "GSE205335"], ["#4C78A8", "#54A24B", "#E45756"]):
            sub = elig[elig.dataset == ds]
            ax.scatter(sub["CLDN4_mean"], sub[f"auc_{name}"], s=18, alpha=0.75, label=ds, c=c)
        ax.set_xlabel("malignant CLDN4 mean")
        ax.set_ylabel(title)
        ax.legend(fontsize=7, frameon=False)
    fig.suptitle("Patient-level CLDN4 vs regulon AUCell (raw scores; tests use within-dataset ranks)", fontsize=10)
    fig.tight_layout()
    fig.savefig(out / "fig_scatter_cldn4_aucell.png", dpi=140)
    fig.savefig(out / "fig_scatter_cldn4_aucell.pdf")
    plt.close(fig)


def md_table(df: pd.DataFrame, cols: list[str]) -> str:
    lines = ["| " + " | ".join(cols) + " |", "|" + "|".join(["---"] * len(cols)) + "|"]
    for rec in df.itertuples(index=False):
        row = []
        for c in cols:
            val = getattr(rec, c)
            if isinstance(val, float):
                if not np.isfinite(val):
                    row.append("NA")
                elif c.endswith("_p"):
                    row.append(fmt_p(val))
                elif abs(val) >= 100:
                    row.append(f"{val:.1f}")
                else:
                    row.append(f"{val:.3f}" if c.startswith(("n_", "n")) is False else f"{val:.0f}")
            else:
                row.append(str(val))
        lines.append("| " + " | ".join(row) + " |")
    return "\n".join(lines)


def write_finding(patients: pd.DataFrame, contrasts: pd.DataFrame, ntab: pd.DataFrame, out: Path) -> None:
    elig = patients[patients.eligible == 1]
    n_trip = int(len(elig))
    n_q1 = int((elig.q_mean == "Q1").sum())
    n_q4 = int((elig.q_mean == "Q4").sum())
    counts = elig.groupby("dataset").size().to_dict()

    def pick(scope, pred, regulon):
        hit = contrasts[
            (contrasts.scope == scope)
            & (contrasts.predictor == pred)
            & (contrasts.regulon == regulon)
        ]
        return hit.iloc[0] if len(hit) else None

    primary_regs = [
        ("AUCell:IFN_STAT1", "IFN (STAT1/IRF + ISGs)"),
        ("AUCell:MHC_NLRC5", "MHC / antigen presentation"),
        ("AUCell:TJ_barrier", "TJ / barrier (CLDN4 held out)"),
        ("AUCell:ELF3_CollecTRI", "ELF3 CollecTRI (A10 given)"),
        ("AUCell:HK_CTRL", "housekeeping control"),
    ]

    lines = [
        "# Finding — triple-merge CLDN4-only AUCell (GSE131907 + GSE148071 + GSE205335)",
        "",
        "ADDITIVE. **CLDN4 only.** No TACSTD2 gate and no dual-high score.",
        "Patient is the unit. A10 ELF3–CLDN4 is **given** and is not re-proved.",
        "This slice is the **triple**. The 131907+205335-only SCENIC agent is not",
        "redone; that pair appears only as a sensitivity row.",
        "pySCENIC / cisTarget were not run. AUCell is Aibar recovery on a fixed",
        "prior + background gene universe. p-values are descriptive.",
        "",
        "## Verdict",
        "",
        f"Eligible malignant patients: **n={n_trip}** "
        f"(GSE131907 {counts.get('GSE131907', 0)}, "
        f"GSE148071 {counts.get('GSE148071', 0)}, "
        f"GSE205335 {counts.get('GSE205335', 0)}). "
        f"Stratified Q4 vs Q1 (CLDN4 mean, quartiles **inside each dataset**) "
        f"is **{n_q4} vs {n_q1}**. Do not write n=112 (GEO sum) or cell n.",
        "",
        "Primary tests use **within-dataset ranks** of malignant CLDN4 mean vs",
        "patient-mean AUCell, then a stratified Q4 vs Q1 on those same ranks.",
        "Positive r / ρ = CLDN4-high patients have **higher** regulon activity.",
        "",
    ]

    # triple primary table
    rows = []
    for reg, lab in primary_regs:
        sp = pick("triple", "CLDN4_mean_rank_in_dataset", reg)
        q = pick("triple", "CLDN4_mean_stratified_Q", reg)
        if sp is None:
            continue
        q_r = fmt_num(q.q4q1_r) if q is not None else "NA"
        q_p = fmt_p(q.q4q1_p) if q is not None else "NA"
        q_n = f"{int(q.n_q4)} vs {int(q.n_q1)}" if q is not None else "NA"
        rows.append(
            f"| {lab} | {int(sp.n_patients)} | {fmt_num(sp.spearman_rho)} ({fmt_p(sp.spearman_p)}) | "
            f"{q_r} ({q_p}; {q_n}) |"
        )
    lines += [
        "| Regulon | n | Spearman ρ (p) | stratified Q4 vs Q1 r (p; n_Q4/n_Q1) |",
        "|---|---:|---|---|",
        *rows,
        "",
    ]

    # interpret
    ifn = pick("triple", "CLDN4_mean_rank_in_dataset", "AUCell:IFN_STAT1")
    mhc = pick("triple", "CLDN4_mean_rank_in_dataset", "AUCell:MHC_NLRC5")
    tj = pick("triple", "CLDN4_mean_rank_in_dataset", "AUCell:TJ_barrier")
    elf = pick("triple", "CLDN4_mean_rank_in_dataset", "AUCell:ELF3_CollecTRI")
    hk = pick("triple", "CLDN4_mean_rank_in_dataset", "AUCell:HK_CTRL")

    def tone(rec, name):
        if rec is None or not np.isfinite(rec.spearman_rho):
            return f"{name}: not scored."
        direction = "higher" if rec.spearman_rho > 0 else "lower"
        hit = "tracks" if rec.spearman_p < 0.05 else "does not significantly track"
        return (
            f"{name} {hit} CLDN4 ({direction} in CLDN4-high; "
            f"ρ={fmt_num(rec.spearman_rho)}, p={fmt_p(rec.spearman_p)}, n={int(rec.n_patients)})."
        )

    lines += [
        "**What holds.** " + " ".join([tone(ifn, "IFN"), tone(mhc, "MHC"), tone(tj, "TJ")]),
        "",
        "**A10 given.** "
        + (
            tone(elf, "ELF3 CollecTRI (CLDN4 held out of the prior)")
            + " This is supporting, not a new ELF3–CLDN4 claim."
        ),
        "",
        "**Control.** " + tone(hk, "Housekeeping AUCell"),
        "",
        "Do not upgrade a cell-level pattern to this patient n. Do not read",
        "GSE205335 RECIST as MPR. GSE131907 has no ICI labels.",
        "",
        "## Honest n",
        "",
        "| item | n | note |",
        "|---|---:|---|",
    ]
    for rec in ntab.itertuples(index=False):
        lines.append(f"| {rec.item} | {int(rec.n)} | {rec.note} |")

    lines += [
        "",
        "Quartiles are assigned **inside each dataset** on eligible patients,",
        "then tails are stacked. n_compared = n_Q1 + n_Q4, not the GEO n and",
        "not the continuous n.",
        "",
        "## Per-dataset AUCell (CLDN4 mean, raw scores)",
        "",
        "| Dataset | Regulon | n | ρ (p) | Q4 vs Q1 r (p; n_Q4/n_Q1) |",
        "|---|---|---:|---|---|",
    ]
    for ds in ["GSE131907", "GSE148071", "GSE205335"]:
        for reg, lab in primary_regs[:3]:
            sp = pick(ds, "CLDN4_mean", reg)
            if sp is None:
                continue
            lines.append(
                f"| {ds} | {lab} | {int(sp.n_patients)} | "
                f"{fmt_num(sp.spearman_rho)} ({fmt_p(sp.spearman_p)}) | "
                f"{fmt_num(sp.q4q1_r)} ({fmt_p(sp.q4q1_p)}; "
                f"{int(sp.n_q4) if np.isfinite(sp.n_q4) else 'NA'} vs "
                f"{int(sp.n_q1) if np.isfinite(sp.n_q1) else 'NA'}) |"
            )

    pair_ifn = pick("pair_131907_205335_sensitivity", "CLDN4_mean_rank_in_dataset", "AUCell:IFN_STAT1")
    lines += [
        "",
        "## Pair 131907+205335 (sensitivity only — not this agent's claim)",
        "",
    ]
    if pair_ifn is not None:
        lines.append(
            f"Eligible pair n={int(pair_ifn.n_patients)}. IFN AUCell ρ="
            f"{fmt_num(pair_ifn.spearman_rho)} (p={fmt_p(pair_ifn.spearman_p)}). "
            "The primary is the triple that adds GSE148071."
        )
    else:
        lines.append("Pair sensitivity row was not scored.")

    lines += [
        "",
        "## Methods (this slice)",
        "",
        "- Predictor: malignant **CLDN4 only** (mean log1p(CP10k) on GEO UMI;",
        "  TISCH log-normalized values on GSE148071). TACSTD2 is recorded but",
        "  never a gate.",
        "- Malignant labels are **author / TISCH major-lineage**, not a new CNV call.",
        "  GSE131907: epithelial `Malignant cells` / tS1 / tS2 / tS3 in tumor",
        "  sites (tLung, tL/B, mLN, mBrain, PE). GSE148071: TISCH `Malignant`.",
        "  GSE205335: `lineage.sub == Malignant cells`, Normal* tissues dropped.",
        "- Eligible patient: ≥20 malignant cells after that filter. Multiple",
        "  tumor samples from one patient are pooled.",
        "- Regulons: IFN = STAT1/IRF + Hallmark-like ISGs; MHC = NLRC5/CIITA +",
        "  HLA/TAP/immunoproteasome; TJ = ELF3/GRHL + junction genes **minus CLDN4**;",
        "  ELF3 = CollecTRI prior (A10 given; CLDN4 was never in that prior);",
        "  HK = housekeeping control. Gene lists: `resources/regulons.json`.",
        "- AUCell: Aibar linear recovery, threshold 5% of the **fixed universe**",
        "  (regulon members + background). Not full-transcriptome pySCENIC.",
        "  Module score (mean of members) is written beside AUCell.",
        "- Pool: rank CLDN4 and each regulon **within dataset**, then one Spearman.",
        "  Q4 vs Q1: `pd.qcut` on average ranks inside each dataset; Mann–Whitney",
        "  on stacked tails. Rank-biserial r = 2U/(n4 n1) − 1.",
        "- GSE148071 uses TISCH2 `expression.h5` (already log-normalized).",
        "  Ranks are within-dataset, so the scale is not mixed into one AUCell.",
        "",
        "## What was not done",
        "",
        "- No dual-high TACSTD2×CLDN4 score.",
        "- No re-run of the 131907+205335-only SCENIC agent as the primary.",
        "- No pySCENIC GRNBoost2 / cisTarget motif ranking.",
        "- No ICI response / MPR test (GSE131907 and GSE148071 have none;",
        "  GSE205335 RECIST is not used as MPR).",
        "- No cell-level p-value as the claim (pseudoreplication).",
        "",
        "Tables: `tables/regulon_table.tsv` (primary), `tables/patient_scores.tsv`,",
        "`tables/honest_n.tsv`, `tables/contrasts_all.tsv`.",
        "",
        "Figures: `figures/fig_q4q1_aucell.png`, `figures/fig_scatter_cldn4_aucell.png`.",
        "",
        "Reproduce: `python3 methods/triple_scrna_scenic_cldn4/download.py && ",
        "python3 methods/triple_scrna_scenic_cldn4/analyze.py`",
        "",
    ]
    (out / "FINDING.md").write_text("\n".join(lines) + "\n")


def write_regulon_manifest(path: Path, present_by_ds: dict[str, list[str]]) -> None:
    payload = {
        "auc_threshold": AUC_THR,
        "cldn4_held_out_of": ["TJ_barrier", "ELF3_CollecTRI"],
        "regulons": {k: sorted(set(v)) for k, v in REGULONS.items()},
        "background": BACKGROUND,
        "genes_present": present_by_ds,
        "note": (
            "AUCell ranks genes inside this universe only. "
            "CLDN4 is the predictor and is excluded from TJ and ELF3 sets."
        ),
    }
    path.write_text(json.dumps(payload, indent=2) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gse131907", type=Path, default=Path("/tmp/gse131907"))
    ap.add_argument("--gse148071", type=Path, default=Path("/tmp/tisch_gse148071"))
    ap.add_argument("--gse205335", type=Path, default=Path("/tmp/gse205335"))
    args = ap.parse_args()
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)

    wanted = set(wanted_genes())
    parts = []
    present = {}
    d131 = load_gse131907(args.gse131907, wanted)
    parts.append(d131)
    present["GSE131907"] = sorted(g for g in wanted if g in d131.columns)
    d148 = load_gse148071(args.gse148071, wanted)
    parts.append(d148)
    present["GSE148071"] = sorted(g for g in wanted if g in d148.columns)
    d205 = load_gse205335(args.gse205335, wanted)
    parts.append(d205)
    present["GSE205335"] = sorted(g for g in wanted if g in d205.columns)

    universe = sorted(set.intersection(*(set(v) for v in present.values())))
    print(f"[universe] shared genes={len(universe)}", flush=True)
    scored_parts = [score_cells(p, universe) for p in parts]
    cells = pd.concat(scored_parts, ignore_index=True)
    patients = patient_table(cells)
    contrasts = contrast_rows(patients)
    ntab = honest_n(patients, cells)

    # Primary regulon table: triple AUCell, CLDN4 mean ranks + stratified Q
    prim = contrasts[
        (contrasts.scope == "triple")
        & (contrasts.regulon.str.startswith("AUCell:"))
        & (contrasts.predictor.isin(["CLDN4_mean_rank_in_dataset", "CLDN4_mean_stratified_Q"]))
    ].copy()
    prim.to_csv(TABLES / "regulon_table.tsv", sep="\t", index=False)
    contrasts.to_csv(TABLES / "contrasts_all.tsv", sep="\t", index=False)
    patients.to_csv(TABLES / "patient_scores.tsv", sep="\t", index=False)
    ntab.to_csv(TABLES / "honest_n.tsv", sep="\t", index=False)
    cells[["dataset", "patient", "sample", "site", "CLDN4", "CLDN4_pos"] + [c for c in cells.columns if c.startswith("auc_")]].to_csv(
        TABLES / "cell_scores.tsv.gz", sep="\t", index=False
    )
    write_regulon_manifest(HERE / "resources" / "regulons.json", present)
    plot_boxes(patients, FIGS)
    write_finding(patients, contrasts, ntab, HERE)

    summary = {
        "n_triple_eligible": int((patients.eligible == 1).sum()),
        "n_by_dataset": patients[patients.eligible == 1].groupby("dataset").size().to_dict(),
        "n_universe": len(universe),
        "primary_table": TABLES.joinpath("regulon_table.tsv").as_posix(),
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
