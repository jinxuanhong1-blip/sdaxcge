#!/usr/bin/env python3
"""Winning-pair DoRothEA wmean: malignant CLDN4-high vs low.

ADDITIVE. CLDN4 only. GSE131907 + GSE205335 only.
GSE207422 and GSE148071 are not loaded.
decoupleR is not imported; scoring is documented wmean.
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

ROOT = Path(__file__).resolve().parents[1]
TUMOR_ORIGINS = ("tLung", "tL/B", "mLN", "mBrain")
MALIG_SUBTYPES = {"Malignant cells", "tS1", "tS2", "tS3"}
MIN_MAL = 20
MIN_TAIL = 8
MIN_MEDIAN_SIDE = 10
MIN_TF_TARGETS = 5

# Pre-specified program TFs (DoRothEA A+B+C only; missing TFs are not imputed).
PROGRAM_TFS = {
    "IFN": ["STAT1", "STAT2", "IRF1", "IRF2", "IRF3", "IRF7", "IRF8", "IRF9", "STAT3"],
    "MHC": ["RFX5", "IRF1"],
    "TJ": ["GRHL2", "KLF4", "ELF3", "TFAP2A"],
    "keratin": ["TP63", "KLF5", "SOX2"],
}
FOCUS_ORDER = [
    ("IFN", "STAT1"),
    ("IFN", "STAT2"),
    ("IFN", "IRF1"),
    ("IFN", "IRF2"),
    ("IFN", "IRF3"),
    ("IFN", "IRF7"),
    ("IFN", "IRF8"),
    ("IFN", "IRF9"),
    ("IFN", "STAT3"),
    ("MHC", "RFX5"),
    ("TJ", "GRHL2"),
    ("TJ", "KLF4"),
    ("TJ", "ELF3"),
    ("TJ", "TFAP2A"),
    ("keratin", "TP63"),
    ("keratin", "KLF5"),
    ("keratin", "SOX2"),
]
MISSING_FROM_DOROTHEA = ["CIITA", "NLRC5", "GRHL1", "GRHL3", "OVOL1", "OVOL2", "RFXANK", "RFXAP"]

TJ_NO_CLDN4 = [
    "CLDN1", "CLDN3", "CLDN7", "OCLN", "TJP1", "TJP2", "F11R", "PARD3",
    "MARVELD2", "CGN", "CRB3", "JAM3",
]
KERATIN_ALL = ["KRT7", "KRT8", "KRT18", "KRT19", "KRT5", "KRT6A", "KRT6B", "KRT14", "KRT17"]
ISG_CORE = [
    "ISG15", "IFI6", "IFI27", "IFI44", "IFI44L", "IFIT1", "IFIT2", "IFIT3",
    "IFIT5", "IFITM1", "IFITM2", "IFITM3", "MX1", "MX2", "OAS1", "OAS2",
    "OAS3", "OASL", "RSAD2", "USP18", "BST2", "XAF1", "STAT1", "STAT2",
    "IRF7", "IRF9", "DDX58", "IFIH1", "SAMD9", "SAMD9L", "HERC5", "HERC6",
    "EPSTI1", "CMPK2", "PARP9", "DTX3L", "LY6E", "SP100", "SP110", "PLSCR1",
]
MHC1_APM = [
    "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "B2M", "TAP1", "TAP2",
    "TAPBP", "PSMB8", "PSMB9", "PSMB10", "PSME1", "PSME2", "NLRC5",
    "ERAP1", "ERAP2", "CALR", "PDIA3", "CANX", "IRF1",
]
CTRL_OXPHOS = [
    "NDUFA1", "NDUFB3", "COX5A", "COX7A2", "UQCRC1", "SDHA",
    "ATP5F1A", "ATP5F1B", "CYCS", "VDAC1",
]
EXTRA_GENES = ["CLDN4", "TACSTD2", "EPCAM", "PTPRC"]
MODULES = {
    "ifn_isg": ISG_CORE,
    "mhc1_apm": MHC1_APM,
    "tj_no_cldn4": TJ_NO_CLDN4,
    "keratin": KERATIN_ALL,
    "ctrl_oxphos": CTRL_OXPHOS,
}


def fmt_p(value: float) -> str:
    if not np.isfinite(value):
        return "NA"
    return f"{value:.2e}" if value < 0.001 else f"{value:.3g}"


def fmt_num(value: float, digits: int = 3) -> str:
    if not np.isfinite(value):
        return "NA"
    return f"{value:+.{digits}f}" if value != 0 else f"{value:.{digits}f}"


def all_wanted_genes(net: pd.DataFrame) -> set[str]:
    genes = set(net["target"]) | set(EXTRA_GENES)
    for block in MODULES.values():
        genes.update(block)
    genes.update(tf for _, tf in FOCUS_ORDER)
    return genes


def load_dorothea(path: Path) -> pd.DataFrame:
    net = pd.read_csv(path, sep="\t")
    need = {"source", "target", "mor", "confidence"}
    if not need.issubset(net.columns):
        raise SystemExit(f"DoRothEA table missing {need - set(net.columns)}")
    net = net[net["confidence"].isin(["A", "B", "C"])].copy()
    net["mor"] = net["mor"].astype(float)
    return net.drop_duplicates(["source", "target"], keep="first")


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
            if key not in fields:
                fields[key] = values
    n = len(fields["title"])
    return pd.DataFrame({k: v for k, v in fields.items() if len(v) == n})


def parse_geo_soft(path: Path) -> pd.DataFrame:
    opener = gzip.open if path.suffix == ".gz" else open
    records: list[dict[str, str]] = []
    current: dict[str, str] | None = None
    descriptions: list[str] = []
    titles: list[str] = []
    with opener(path, "rt", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
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
    metadata = pd.DataFrame(records)
    read_end = metadata["platform"].str.extract(r"Single Cell ([35])'")[0]
    metadata["orig.ident"] = (
        metadata["description"].str.replace("_", "-", regex=False) + "-" + read_end + "P"
    )
    return metadata.rename(
        columns={"tumor stage": "tumor_stage", "cancer subtype": "cancer_subtype"}
    )


def stream_matrix(matrix_path: Path, wanted: set[str], keep: np.ndarray):
    keep_idx = np.flatnonzero(keep)
    found: dict[str, np.ndarray] = {}
    with gzip.open(matrix_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        cell_ids = header[1:]
        n = len(cell_ids)
        if keep.size != n:
            raise ValueError(f"keep mask {keep.size} != matrix cells {n}")
        n_umi = np.zeros(keep_idx.size, dtype=np.float64)
        n_streamed = 0
        for line in handle:
            gene, sep, rest = line.partition("\t")
            if not sep:
                continue
            gene = gene.split(".")[0]
            arr = np.fromstring(rest, sep="\t", dtype=np.float32)
            if arr.size != n:
                raise ValueError(f"{gene}: {arr.size} != {n}")
            kept = arr[keep_idx]
            n_umi += kept
            if gene in wanted:
                found[gene] = np.asarray(kept, dtype=np.float32)
            n_streamed += 1
            if n_streamed % 4000 == 0:
                print(f"  GSE131907 stream genes={n_streamed} stored={len(found)}", flush=True)
    print(
        f"GSE131907 stream done genes={n_streamed} kept_cells={keep_idx.size} stored={len(found)}",
        flush=True,
    )
    return [cell_ids[i] for i in keep_idx], found, n_umi, n_streamed


def load_rds_genes(path: Path, wanted: set[str], keep_barcodes: set[str] | None):
    with tempfile.TemporaryDirectory(prefix="gse205335-") as tmp:
        matrix_path = path
        if path.suffix == ".gz":
            matrix_path = Path(tmp) / path.stem
            print(f"decompress {path.name}", flush=True)
            with gzip.open(path, "rb") as source, matrix_path.open("wb") as dest:
                shutil.copyfileobj(source, dest, 16 * 1024 * 1024)
        print("read GSE205335 RDS", flush=True)
        import rdata

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
    matrix = sparse.csc_matrix((obj.x, obj.i, obj.p), shape=tuple(obj.Dim), copy=False)
    if keep_barcodes is not None:
        keep = np.fromiter((b in keep_barcodes for b in barcodes), dtype=bool, count=len(barcodes))
        matrix = matrix[:, keep]
        barcodes = barcodes[keep]
    library_umi = np.asarray(matrix.sum(axis=0)).ravel()
    matrix = matrix.tocsr()
    name_to_row = {g: i for i, g in enumerate(genes)}
    extracted: dict[str, np.ndarray] = {}
    for gene in sorted(wanted):
        row = name_to_row.get(gene)
        if row is None:
            continue
        extracted[gene] = np.asarray(matrix.getrow(row).toarray(), dtype=np.float32).ravel()
    print(f"GSE205335 extracted {len(extracted)} / {len(wanted)} genes; cells={len(barcodes)}", flush=True)
    return extracted, library_umi, barcodes, set(genes.tolist())


def log1p_cp10k(umi: dict[str, np.ndarray], lib: np.ndarray) -> dict[str, np.ndarray]:
    scale = np.maximum(lib, 1.0)
    return {g: np.log1p(umi[g] / scale * 1e4).astype(np.float32) for g in umi}


def load_gse131907(args, wanted: set[str]) -> dict:
    print("== GSE131907 ==", flush=True)
    ann = pd.read_csv(args.gse131907_ann, sep="\t", dtype=str)
    with gzip.open(args.gse131907_matrix, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
    cell_ids = header[1:]
    per = ann.set_index("Index").reindex(cell_ids)
    if per["Sample"].isna().any():
        raise SystemExit("GSE131907 matrix/annotation mismatch")
    epi = per["Cell_type"].eq("Epithelial cells")
    mal = epi & per["Cell_subtype"].isin(MALIG_SUBTYPES)
    tumor = per["Sample_Origin"].isin(TUMOR_ORIGINS)
    keep = (mal & tumor).to_numpy()
    print(f"GSE131907 malignant tumor cells: {int(keep.sum())} / {len(cell_ids)}", flush=True)
    kept_ids, found, n_umi, n_streamed = stream_matrix(args.gse131907_matrix, wanted, keep)
    if "CLDN4" not in found:
        raise SystemExit("CLDN4 missing from GSE131907 UMI")
    kept_ann = per.loc[keep].reset_index()
    series = parse_series_matrix(args.gse131907_series)
    sample_meta = series.rename(columns={"title": "Sample"})
    pmap = sample_meta.set_index("Sample")["patient_id"]
    patient = np.array([pmap[s] if s in pmap.index else s for s in kept_ann["Sample"]], dtype=object)
    log_cp = log1p_cp10k(found, n_umi)
    return {
        "cohort": "GSE131907",
        "log_cp": log_cp,
        "umi": found,
        "lib": n_umi,
        "patient": patient,
        "sample": kept_ann["Sample"].to_numpy(),
        "origin": kept_ann["Sample_Origin"].to_numpy(),
        "subtype": kept_ann["Cell_subtype"].to_numpy(),
        "genes": set(found),
        "n_streamed": n_streamed,
        "n_cells_matrix": int(len(cell_ids)),
        "n_malignant": int(keep.sum()),
    }


def load_gse205335(args, wanted: set[str]) -> dict:
    print("== GSE205335 ==", flush=True)
    identities = pd.read_csv(args.gse205335_identities, sep="\t")
    if identities["barcode"].duplicated().any():
        raise ValueError("GSE205335 identity barcodes are not unique")
    metadata = parse_geo_soft(args.gse205335_soft)
    mal_barcodes = set(identities.loc[identities["lineage.sub"].eq("Malignant cells"), "barcode"])
    extracted, library_umi, barcodes, genes = load_rds_genes(
        args.gse205335_matrix, wanted, mal_barcodes
    )
    if "CLDN4" not in extracted:
        raise SystemExit("CLDN4 missing from GSE205335 UMI")
    indexed = identities.set_index("barcode")
    cells = indexed.loc[barcodes].reset_index()
    cells["total_umi"] = library_umi
    cells = cells.merge(
        metadata[["orig.ident", "gsm", "patient", "tissue", "recist", "cancer_subtype"]],
        on="orig.ident",
        how="left",
        validate="many_to_one",
    )
    if cells["patient"].isna().any():
        raise ValueError("GSE205335 identity samples did not match GEO metadata")
    log_cp = log1p_cp10k(extracted, cells["total_umi"].to_numpy())
    return {
        "cohort": "GSE205335",
        "log_cp": log_cp,
        "umi": extracted,
        "lib": cells["total_umi"].to_numpy(),
        "patient": cells["patient"].to_numpy(),
        "sample": cells["orig.ident"].to_numpy(),
        "origin": cells["tissue"].fillna("").to_numpy(),
        "subtype": np.array(["Malignant cells"] * len(cells), dtype=object),
        "histology": cells["cancer_subtype"].to_numpy(),
        "recist": cells["recist"].to_numpy(),
        "genes": genes,
        "n_cells_matrix": int(len(identities)),
        "n_malignant": int(len(cells)),
    }


def _group_index(labels: np.ndarray) -> list[tuple[object, np.ndarray]]:
    order = pd.Series(labels).groupby(labels, sort=True).indices
    return [(key, np.asarray(idx)) for key, idx in order.items()]


def qcut_tails(values: np.ndarray) -> tuple[np.ndarray, np.ndarray] | None:
    if values.size < (MIN_TAIL * 2):
        return None
    ranks = pd.Series(values).rank(method="average")
    try:
        qs = pd.qcut(ranks, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")
    except ValueError:
        return None
    if qs.nunique() < 4:
        return None
    low = np.flatnonzero(qs.eq("Q1").to_numpy())
    high = np.flatnonzero(qs.eq("Q4").to_numpy())
    if low.size < MIN_TAIL or high.size < MIN_TAIL:
        return None
    return high, low


def median_split(values: np.ndarray) -> tuple[np.ndarray, np.ndarray] | None:
    if values.size < (MIN_MEDIAN_SIDE * 2):
        return None
    med = float(np.median(values))
    if med == 0:
        high = np.flatnonzero(values > 0)
        low = np.flatnonzero(values == 0)
    else:
        high = np.flatnonzero(values >= med)
        low = np.flatnonzero(values < med)
    if high.size < MIN_MEDIAN_SIDE or low.size < MIN_MEDIAN_SIDE:
        return None
    return high, low


def module_score(log_cp: dict[str, np.ndarray], genes: list[str], idx: np.ndarray) -> float:
    present = [g for g in genes if g in log_cp]
    if not present or idx.size == 0:
        return float("nan")
    stacked = np.vstack([log_cp[g][idx] for g in present])
    return float(stacked.mean())


def score_tfs(
    log_cp: dict[str, np.ndarray],
    net: pd.DataFrame,
    shared_targets: set[str],
) -> tuple[dict[str, np.ndarray], pd.DataFrame]:
    """Per-cell wmean. activity[tf] shape = n_cells."""
    n_cells = next(iter(log_cp.values())).size
    gene_index = {g: i for i, g in enumerate(sorted(g for g in log_cp if g in shared_targets))}
    if not gene_index:
        raise SystemExit("no shared DoRothEA targets in the expression slice")
    mat = np.vstack([log_cp[g] for g in gene_index]).astype(np.float32, copy=False)
    # mat is n_genes x n_cells
    rows = []
    activities: dict[str, np.ndarray] = {}
    for tf, sub in net.groupby("source", sort=True):
        sub = sub[sub["target"].isin(gene_index)]
        if len(sub) < MIN_TF_TARGETS:
            continue
        idx = np.array([gene_index[t] for t in sub["target"]], dtype=int)
        w = sub["mor"].to_numpy(dtype=np.float64)
        denom = float(np.abs(w).sum())
        if denom <= 0:
            continue
        act = (w @ mat[idx]) / denom
        activities[tf] = act.astype(np.float32)
        rows.append(
            {
                "tf": tf,
                "n_targets_used": int(len(sub)),
                "n_targets_abc": int((net["source"] == tf).sum()),
                "n_repressor": int((sub["mor"] < 0).sum()),
                "weight_l1": denom,
            }
        )
    return activities, pd.DataFrame(rows)


def wilcoxon_paired(high: np.ndarray, low: np.ndarray) -> dict:
    high = np.asarray(high, dtype=float)
    low = np.asarray(low, dtype=float)
    mask = np.isfinite(high) & np.isfinite(low)
    high, low = high[mask], low[mask]
    n = int(high.size)
    if n < 3:
        return {
            "n": n,
            "mean_high": float("nan"),
            "mean_low": float("nan"),
            "median_high": float("nan"),
            "median_low": float("nan"),
            "delta_mean": float("nan"),
            "delta_median": float("nan"),
            "n_high_gt_low": 0,
            "W": float("nan"),
            "p": float("nan"),
        }
    diff = high - low
    try:
        w_stat, p = stats.wilcoxon(high, low, alternative="two-sided", zero_method="wilcox", method="auto")
    except ValueError:
        w_stat, p = float("nan"), float("nan")
    return {
        "n": n,
        "mean_high": float(high.mean()),
        "mean_low": float(low.mean()),
        "median_high": float(np.median(high)),
        "median_low": float(np.median(low)),
        "delta_mean": float(diff.mean()),
        "delta_median": float(np.median(diff)),
        "n_high_gt_low": int((high > low).sum()),
        "W": float(w_stat) if np.isfinite(w_stat) else float("nan"),
        "p": float(p) if np.isfinite(p) else float("nan"),
    }


def spearman(x, y) -> tuple[float, float, int]:
    xa = np.asarray(x, dtype=float)
    ya = np.asarray(y, dtype=float)
    mask = np.isfinite(xa) & np.isfinite(ya)
    xa, ya = xa[mask], ya[mask]
    n = int(xa.size)
    if n < 3:
        return float("nan"), float("nan"), n
    rho, p = stats.spearmanr(xa, ya)
    return float(rho), float(p), n


def patient_rows(ds: dict, activities: dict[str, np.ndarray], tf_meta: pd.DataFrame) -> pd.DataFrame:
    rows = []
    cldn = ds["log_cp"]["CLDN4"]
    for patient, idx in _group_index(ds["patient"]):
        n_mal = int(idx.size)
        cldn_p = cldn[idx]
        q = qcut_tails(cldn_p)
        m = median_split(cldn_p) if n_mal >= MIN_MAL else None
        q_ok = bool(n_mal >= MIN_MAL and q is not None)
        med_ok = bool(n_mal >= MIN_MAL and m is not None)
        rec = {
            "cohort": ds["cohort"],
            "patient": patient,
            "n_malignant": n_mal,
            "n_samples": int(len(set(ds["sample"][idx].tolist()))),
            "origins": ",".join(sorted(set(str(x) for x in ds["origin"][idx].tolist()))),
            "cldn4_mean": float(cldn_p.mean()),
            "cldn4_pct_pos": float(100.0 * (ds["umi"]["CLDN4"][idx] > 0).mean()),
            "eligible_q4q1": q_ok,
            "eligible_median": med_ok,
            "n_q4": int(q[0].size) if q_ok else 0,
            "n_q1": int(q[1].size) if q_ok else 0,
            "n_med_high": int(m[0].size) if med_ok else 0,
            "n_med_low": int(m[1].size) if med_ok else 0,
        }
        if "histology" in ds:
            rec["histology"] = ",".join(sorted(set(str(x) for x in ds["histology"][idx].tolist())))
            rec["recist"] = ",".join(sorted(set(str(x) for x in ds["recist"][idx].tolist())))
        for name, genes in MODULES.items():
            rec[f"all_{name}"] = module_score(ds["log_cp"], genes, idx)
            if q_ok:
                rec[f"q4_{name}"] = module_score(ds["log_cp"], genes, idx[q[0]])
                rec[f"q1_{name}"] = module_score(ds["log_cp"], genes, idx[q[1]])
            if med_ok:
                rec[f"medh_{name}"] = module_score(ds["log_cp"], genes, idx[m[0]])
                rec[f"medl_{name}"] = module_score(ds["log_cp"], genes, idx[m[1]])
        for tf, act in activities.items():
            rec[f"all_{tf}"] = float(act[idx].mean())
            if q_ok:
                rec[f"q4_{tf}"] = float(act[idx[q[0]]].mean())
                rec[f"q1_{tf}"] = float(act[idx[q[1]]].mean())
            if med_ok:
                rec[f"medh_{tf}"] = float(act[idx[m[0]]].mean())
                rec[f"medl_{tf}"] = float(act[idx[m[1]]].mean())
        rows.append(rec)
    frame = pd.DataFrame(rows)
    frame.attrs["tf_meta"] = tf_meta
    return frame


def test_table(patients: pd.DataFrame, items: list[tuple[str, str]], split: str) -> pd.DataFrame:
    high_key = "q4" if split == "q4q1" else "medh"
    low_key = "q1" if split == "q4q1" else "medl"
    elig = "eligible_q4q1" if split == "q4q1" else "eligible_median"
    rows = []
    for cohort, sub in [("GSE131907", patients[patients.cohort.eq("GSE131907")]),
                        ("GSE205335", patients[patients.cohort.eq("GSE205335")]),
                        ("pooled", patients)]:
        use = sub[sub[elig]].copy()
        for program, name in items:
            hcol, lcol = f"{high_key}_{name}", f"{low_key}_{name}"
            if hcol not in use.columns:
                continue
            stats_row = wilcoxon_paired(use[hcol].to_numpy(), use[lcol].to_numpy())
            rho, rp, rn = spearman(use["cldn4_mean"], use[f"all_{name}"]) if f"all_{name}" in use else (np.nan, np.nan, 0)
            rows.append(
                {
                    "split": split,
                    "cohort": cohort,
                    "program": program,
                    "name": name,
                    "n_attempted": int(len(sub)),
                    "n": stats_row["n"],
                    "n_high_gt_low": stats_row["n_high_gt_low"],
                    "mean_high": stats_row["mean_high"],
                    "mean_low": stats_row["mean_low"],
                    "delta_mean": stats_row["delta_mean"],
                    "delta_median": stats_row["delta_median"],
                    "W": stats_row["W"],
                    "p": stats_row["p"],
                    "spearman_rho": rho,
                    "spearman_p": rp,
                    "spearman_n": rn,
                    "thin": stats_row["n"] < 8,
                }
            )
    return pd.DataFrame(rows)


def md_tf_table(tf_tests: pd.DataFrame, tf_meta: pd.DataFrame, cohort: str = "pooled") -> str:
    sub = tf_tests[(tf_tests.cohort.eq(cohort)) & (tf_tests.split.eq("q4q1"))].copy()
    meta = tf_meta.set_index("tf")
    lines = [
        "| program | TF | n_targets | n | high vs low | Δ mean | high>low | W | p |",
        "|---|---|---:|---:|---|---:|---:|---:|---|",
    ]
    for program, tf in FOCUS_ORDER:
        row = sub[sub.name.eq(tf)]
        if row.empty:
            n_t = int(meta.loc[tf, "n_targets_used"]) if tf in meta.index else 0
            lines.append(f"| {program} | {tf} | {n_t} | 0 | NA | NA | NA | NA | not scored |")
            continue
        r = row.iloc[0]
        n_t = int(meta.loc[tf, "n_targets_used"]) if tf in meta.index else 0
        lines.append(
            f"| {program} | {tf} | {n_t} | {int(r.n)} | "
            f"{r.mean_high:.3f} vs {r.mean_low:.3f} | {fmt_num(r.delta_mean)} | "
            f"{int(r.n_high_gt_low)}/{int(r.n)} | {r.W:.0f} | **{fmt_p(r.p)}** |"
            if r.p < 0.05
            else f"| {program} | {tf} | {n_t} | {int(r.n)} | "
            f"{r.mean_high:.3f} vs {r.mean_low:.3f} | {fmt_num(r.delta_mean)} | "
            f"{int(r.n_high_gt_low)}/{int(r.n)} | {r.W:.0f} | {fmt_p(r.p)} |"
        )
    return "\n".join(lines)


def write_figures(patients: pd.DataFrame, tf_tests: pd.DataFrame, figdir: Path) -> list[str]:
    figdir.mkdir(parents=True, exist_ok=True)
    written = []
    use = patients[patients.eligible_q4q1].copy()
    key_tfs = ["STAT1", "IRF1", "RFX5", "GRHL2", "TP63"]

    # Honest n
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.0), sharey=True)
    for ax, cohort in zip(axes, ["GSE131907", "GSE205335"]):
        sub = patients[patients.cohort.eq(cohort)].sort_values("n_malignant")
        colors = ["#1b4f72" if ok else "#d0d0d0" for ok in sub.eligible_q4q1]
        ax.bar(np.arange(len(sub)), sub.n_malignant, color=colors, width=0.85)
        ax.axhline(MIN_MAL, color="#a04000", ls="--", lw=1, label=f"floor n={MIN_MAL}")
        ax.set_title(f"{cohort}  eligible {int(sub.eligible_q4q1.sum())}/{len(sub)}")
        ax.set_xlabel("patient (sorted)")
        ax.set_xticks([])
    axes[0].set_ylabel("malignant cells")
    axes[0].legend(frameon=False, loc="upper left")
    fig.suptitle("Honest n — malignant occupancy (blue = Q4 vs Q1 eligible)", y=1.02)
    fig.tight_layout()
    path = figdir / "fig_extra_honest_n.png"
    fig.savefig(path, dpi=140, bbox_inches="tight")
    fig.savefig(figdir / "fig_extra_honest_n.pdf", bbox_inches="tight")
    plt.close(fig)
    written.append(str(path.relative_to(ROOT)))

    # Paired key TFs
    fig, axes = plt.subplots(1, len(key_tfs), figsize=(12.5, 3.6), sharey=False)
    for ax, tf in zip(axes, key_tfs):
        if f"q4_{tf}" not in use.columns:
            ax.set_title(f"{tf} missing")
            continue
        for cohort, color in [("GSE131907", "#1b4f72"), ("GSE205335", "#b85c38")]:
            sub = use[use.cohort.eq(cohort)]
            for rec in sub.itertuples(index=False):
                ax.plot([0, 1], [getattr(rec, f"q1_{tf}"), getattr(rec, f"q4_{tf}")],
                        color=color, alpha=0.35, lw=1)
            ax.scatter(np.zeros(len(sub)), sub[f"q1_{tf}"], s=18, color=color, zorder=3)
            ax.scatter(np.ones(len(sub)), sub[f"q4_{tf}"], s=18, color=color, zorder=3, label=cohort)
        ax.set_xticks([0, 1], ["Q1", "Q4"])
        ax.set_title(tf)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    axes[0].set_ylabel("DoRothEA wmean")
    axes[-1].legend(frameon=False, fontsize=7, loc="best")
    fig.suptitle("Extra — paired malignant CLDN4 Q1 vs Q4 TF activity", y=1.03)
    fig.tight_layout()
    path = figdir / "fig_extra_paired_tf.png"
    fig.savefig(path, dpi=140, bbox_inches="tight")
    fig.savefig(figdir / "fig_extra_paired_tf.pdf", bbox_inches="tight")
    plt.close(fig)
    written.append(str(path.relative_to(ROOT)))

    # Forest
    forest = tf_tests[(tf_tests.cohort.eq("pooled")) & (tf_tests.split.eq("q4q1"))].copy()
    forest = forest[forest.name.isin([t for _, t in FOCUS_ORDER])]
    forest["label"] = forest.program + " / " + forest.name
    forest = forest.iloc[::-1]
    fig, ax = plt.subplots(figsize=(7.2, 5.8))
    y = np.arange(len(forest))
    colors = ["#1b4f72" if p < 0.05 else "#7f8c8d" for p in forest.p]
    ax.axvline(0, color="#444", lw=0.8)
    ax.scatter(forest.delta_mean, y, c=colors, s=36, zorder=3)
    ax.set_yticks(y, forest.label)
    ax.set_xlabel("Δ wmean (CLDN4 Q4 − Q1), pooled patients")
    ax.set_title(f"Extra — focused TF deltas (n={int(forest.n.max()) if len(forest) else 0})")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    path = figdir / "fig_extra_forest.png"
    fig.savefig(path, dpi=140, bbox_inches="tight")
    fig.savefig(figdir / "fig_extra_forest.pdf", bbox_inches="tight")
    plt.close(fig)
    written.append(str(path.relative_to(ROOT)))

    # Heatmap of patient deltas
    tfs = [t for _, t in FOCUS_ORDER if f"q4_{t}" in use.columns]
    if tfs and len(use):
        mat = np.vstack([(use[f"q4_{t}"] - use[f"q1_{t}"]).to_numpy() for t in tfs])
        fig, ax = plt.subplots(figsize=(max(8, 0.28 * len(use)), 5.6))
        vmax = np.nanpercentile(np.abs(mat), 95) or 1.0
        im = ax.imshow(mat, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        ax.set_yticks(np.arange(len(tfs)), tfs)
        ax.set_xticks(np.arange(len(use)), [f"{c[3:]}:{p}" for c, p in zip(use.cohort, use.patient)],
                      rotation=90, fontsize=6)
        fig.colorbar(im, ax=ax, shrink=0.7, label="Δ wmean Q4−Q1")
        ax.set_title("Extra — per-patient TF Δ (Q4 − Q1)")
        fig.tight_layout()
        path = figdir / "fig_extra_heatmap_delta.png"
        fig.savefig(path, dpi=140, bbox_inches="tight")
        fig.savefig(figdir / "fig_extra_heatmap_delta.pdf", bbox_inches="tight")
        plt.close(fig)
        written.append(str(path.relative_to(ROOT)))

    # Module companion
    fig, axes = plt.subplots(1, 4, figsize=(11.5, 3.5))
    for ax, name, title in zip(
        axes,
        ["ifn_isg", "mhc1_apm", "tj_no_cldn4", "keratin"],
        ["IFN ISG", "MHC-I APM", "TJ (no CLDN4)", "keratin"],
    ):
        for cohort, color in [("GSE131907", "#1b4f72"), ("GSE205335", "#b85c38")]:
            sub = use[use.cohort.eq(cohort)]
            ax.scatter(sub[f"q1_{name}"], sub[f"q4_{name}"], s=22, color=color, label=cohort, alpha=0.85)
        lims = [
            min(ax.get_xlim()[0], ax.get_ylim()[0]),
            max(ax.get_xlim()[1], ax.get_ylim()[1]),
        ]
        ax.plot(lims, lims, color="#888", lw=0.7)
        ax.set_xlabel("Q1")
        ax.set_ylabel("Q4")
        ax.set_title(title)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    axes[-1].legend(frameon=False, fontsize=7)
    fig.suptitle("Extra — companion gene-set modules (not TF activity)", y=1.03)
    fig.tight_layout()
    path = figdir / "fig_extra_modules.png"
    fig.savefig(path, dpi=140, bbox_inches="tight")
    fig.savefig(figdir / "fig_extra_modules.pdf", bbox_inches="tight")
    plt.close(fig)
    written.append(str(path.relative_to(ROOT)))
    return written


def write_finding(
    patients: pd.DataFrame,
    tf_tests: pd.DataFrame,
    mod_tests: pd.DataFrame,
    tf_meta: pd.DataFrame,
    inventory: dict,
    figures: list[str],
    path: Path,
) -> None:
    q = patients[patients.eligible_q4q1]
    n131 = int(q.cohort.eq("GSE131907").sum())
    n205 = int(q.cohort.eq("GSE205335").sum())
    n_pool = int(len(q))
    att131 = int(patients.cohort.eq("GSE131907").sum())
    att205 = int(patients.cohort.eq("GSE205335").sum())
    pooled = tf_tests[(tf_tests.cohort.eq("pooled")) & (tf_tests.split.eq("q4q1"))]

    def pick(tf: str) -> pd.Series | None:
        hit = pooled[pooled.name.eq(tf)]
        return None if hit.empty else hit.iloc[0]

    stat1 = pick("STAT1")
    grhl2 = pick("GRHL2")
    tp63 = pick("TP63")
    rfx5 = pick("RFX5")

    def one_line(row, label):
        if row is None:
            return f"- **{label}:** not scored."
        return (
            f"- **{label}:** n={int(row.n)} paired patients, "
            f"Δ={fmt_num(row.delta_mean)}, high>low {int(row.n_high_gt_low)}/{int(row.n)}, "
            f"p={fmt_p(row.p)}."
        )

    body = f"""# FINDING — winning-pair DoRothEA TF activity (CLDN4-only)

ADDITIVE **CLDN4 only** on the PR #320 winning pair **GSE131907 + GSE205335**.
No TACSTD2∩CLDN4 dual-high. **GSE207422 and GSE148071 are not merged in.**
Patient is the unit. p-values are descriptive.

decoupleR / dorothea R were **not installed**. Activity is documented
DoRothEA **wmean** (Badia-i-Mompel et al. 2022):
wmean = sum(w_g * x_g) / sum(|w_g|) on `log1p(CP10k)` targets,
w = +1 / -1 from OmniPath A+B+C. No VIPER NES.

## Honest n

| Item | n | Note |
|---|---:|---|
| GSE131907 patients with any tumor-origin malignant cells | {att131} | author `Malignant cells` + tS1/tS2/tS3; PE/nLung/nLN out |
| GSE205335 patients with any author-malignant cells | {att205} | `lineage.sub == Malignant cells` |
| **Paired Q4 vs Q1 (floor ≥20 malignant, ≥8/tail)** | **{n_pool}** | **{n131} + {n205}** |
| GSE131907 eligible | {n131} / {att131} | dropped if a tail <8 |
| GSE205335 eligible | {n205} / {att205} | dropped if a tail <8 |
| Cells as n | 0 | not used |

Do **not** cite the attempted header n as the test n. Thin tails stay out.

GSE131907 matrix: {inventory.get("gse131907_cells", "NA")} barcodes streamed,
{inventory.get("gse131907_malignant", "NA")} malignant tumor cells kept.
GSE205335 matrix: {inventory.get("gse205335_cells", "NA")} barcodes,
{inventory.get("gse205335_malignant", "NA")} author-malignant kept.
Shared DoRothEA A+B+C targets in both matrices: {inventory.get("n_shared_targets", "NA")}.
TFs scored (≥5 shared targets): {inventory.get("n_tfs_scored", "NA")}.

Absent from DoRothEA A+B+C (not imputed): {", ".join(MISSING_FROM_DOROTHEA)}.

## Verdict (pooled Q4 vs Q1)

{one_line(stat1, "IFN / STAT1")}
{one_line(pick("IRF1"), "IFN∩MHC / IRF1")}
{one_line(rfx5, "MHC / RFX5")}
{one_line(grhl2, "TJ / GRHL2")}
{one_line(tp63, "keratin / TP63")}

This is TF activity in the **same malignant cells** that define the CLDN4
split, not a T/NK test. The PR #320 T/NK association on this pair is
given and is not re-ranked (Q4 vs Q1 n=23 r=−0.705 p=0.0003; continuous
n=43 ρ=−0.479).

## Primary TF table — pooled patients, Q4 vs Q1

{md_tf_table(tf_tests, tf_meta, "pooled")}

High = within-patient malignant CLDN4 Q4; low = Q1. Mean is the mean of
patient-bin wmean values. **Bold p** is p<0.05 (descriptive).

### Same table by cohort

**GSE131907 (n={n131})**

{md_tf_table(tf_tests, tf_meta, "GSE131907")}

**GSE205335 (n={n205})**

{md_tf_table(tf_tests, tf_meta, "GSE205335")}

## Companion gene-set modules (not the TF table)

Mean `log1p(CP10k)` of the listed genes. CLDN4 is excluded from TJ.

{md_module_table(mod_tests, n_pool)}

## Median-split companion

Same patients are not required. Eligible median-split n =
{int(patients.eligible_median.sum())}
({int(patients[patients.cohort.eq("GSE131907")].eligible_median.sum())} +
{int(patients[patients.cohort.eq("GSE205335")].eligible_median.sum())}).
See `results/tf_tests.tsv` (`split=median`).

## Honest limits

1. **Winning pair only.** GSE207422 and GSE148071 were not downloaded.
2. **No dual-high.** TACSTD2 is not a gate.
3. **wmean, not decoupleR NES.** No permutation / VIPER normalization.
4. **Author labels.** CopyKAT was not re-run. GSE131907 tLung uses tS1/tS2/tS3.
5. **GSE205335 histology mix** (ADC/SQ/SCLC) is part of the honest n.
6. IRF7 / IRF8 regulons are small in A+B+C; they are flagged by `n_targets`.
7. p-values are descriptive. Exact Wilcoxon cannot go below 0.031 at n=6.

## Extra figures

{chr(10).join(f"- `{p}`" for p in figures)}

## Files

- `results/tf_table.tsv` — compact primary TF table (pooled Q4 vs Q1)
- `results/tf_tests.tsv` — all cohort × split tests
- `results/per_patient.tsv` — occupancy + bin activities
- `results/tf_meta.tsv` — targets used per TF
- `results/module_tests.tsv` — companion gene-set tests
- `resources/dorothea_hs_ABC.tsv` — OmniPath A+B+C weights
- `METHODS.md` — wmean formula, floors, labels

## Reproduce

```bash
python3 methods/winpair_131907_205335_dorothea_cldn4/scripts/download.py
python3 methods/winpair_131907_205335_dorothea_cldn4/scripts/analyze.py
```
"""
    path.write_text(body)
    print(f"wrote {path}", flush=True)


def md_module_table(mod_tests: pd.DataFrame, n_pool: int) -> str:
    sub = mod_tests[(mod_tests.cohort.eq("pooled")) & (mod_tests.split.eq("q4q1"))]
    lines = [
        "| module | n | high vs low | Δ mean | high>low | p |",
        "|---|---:|---|---:|---:|---|",
    ]
    labels = {
        "ifn_isg": "IFN ISG core",
        "mhc1_apm": "MHC-I APM",
        "tj_no_cldn4": "TJ (no CLDN4)",
        "keratin": "keratin",
        "ctrl_oxphos": "OXPHOS control",
    }
    for key, label in labels.items():
        row = sub[sub.name.eq(key)]
        if row.empty:
            lines.append(f"| {label} | {n_pool} | NA | NA | NA | NA |")
            continue
        r = row.iloc[0]
        lines.append(
            f"| {label} | {int(r.n)} | {r.mean_high:.3f} vs {r.mean_low:.3f} | "
            f"{fmt_num(r.delta_mean)} | {int(r.n_high_gt_low)}/{int(r.n)} | {fmt_p(r.p)} |"
        )
    return "\n".join(lines)


def compact_tf_table(tf_tests: pd.DataFrame, tf_meta: pd.DataFrame) -> pd.DataFrame:
    sub = tf_tests[(tf_tests.cohort.eq("pooled")) & (tf_tests.split.eq("q4q1"))].copy()
    meta = tf_meta.set_index("tf")
    rows = []
    for program, tf in FOCUS_ORDER:
        hit = sub[sub.name.eq(tf)]
        rec = {"program": program, "tf": tf}
        if tf in meta.index:
            rec["n_targets_used"] = int(meta.loc[tf, "n_targets_used"])
            rec["n_targets_abc"] = int(meta.loc[tf, "n_targets_abc"])
        if hit.empty:
            rec["status"] = "not_scored"
        else:
            r = hit.iloc[0]
            rec.update(
                {
                    "n": int(r.n),
                    "mean_high": r.mean_high,
                    "mean_low": r.mean_low,
                    "delta_mean": r.delta_mean,
                    "n_high_gt_low": int(r.n_high_gt_low),
                    "W": r.W,
                    "p": r.p,
                    "thin": bool(r.thin),
                    "status": "ok",
                }
            )
        rows.append(rec)
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--geo", type=Path, default=Path("/tmp/geo_winpair"))
    parser.add_argument(
        "--dorothea",
        type=Path,
        default=ROOT / "resources" / "dorothea_hs_ABC.tsv",
    )
    args = parser.parse_args()
    args.gse131907_matrix = args.geo / "gse131907/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    args.gse131907_ann = args.geo / "gse131907/GSE131907_Lung_Cancer_cell_annotation.txt.gz"
    args.gse131907_series = args.geo / "gse131907/GSE131907_series_matrix.txt.gz"
    args.gse205335_matrix = args.geo / "gse205335/GSE205335_Lung_IO_UMI_matrix.rds.gz"
    args.gse205335_identities = args.geo / "gse205335/GSE205335_Lung_IO_CellIdentity.txt.gz"
    args.gse205335_soft = args.geo / "gse205335/GSE205335_family.soft.gz"
    for path in (
        args.gse131907_matrix,
        args.gse131907_ann,
        args.gse131907_series,
        args.gse205335_matrix,
        args.gse205335_identities,
        args.gse205335_soft,
        args.dorothea,
    ):
        if not path.exists():
            raise SystemExit(f"missing {path}; run scripts/download.py")

    net = load_dorothea(args.dorothea)
    wanted = all_wanted_genes(net)
    print(f"DoRothEA ABC edges={len(net)} TFs={net.source.nunique()} wanted_genes={len(wanted)}", flush=True)

    ds131 = load_gse131907(args, wanted)
    ds205 = load_gse205335(args, wanted)
    shared = (set(ds131["log_cp"]) & set(ds205["genes"]) & set(net["target"])) | (
        set(ds205["log_cp"]) & set(ds131["genes"] if "genes" in ds131 else set(ds131["log_cp"])) & set(net["target"])
    )
    # Require target present as extracted gene in BOTH expression dicts.
    shared = set(ds131["log_cp"]) & set(ds205["log_cp"]) & set(net["target"])
    print(f"shared DoRothEA targets: {len(shared)}", flush=True)
    net_shared = net[net["target"].isin(shared)].copy()

    act131, meta131 = score_tfs(ds131["log_cp"], net_shared, shared)
    act205, meta205 = score_tfs(ds205["log_cp"], net_shared, shared)
    common_tfs = sorted(set(act131) & set(act205))
    tf_meta = meta131[meta131.tf.isin(common_tfs)].copy()
    print(f"TFs scored in both cohorts: {len(common_tfs)}", flush=True)

    patients = pd.concat(
        [
            patient_rows(ds131, {t: act131[t] for t in common_tfs}, tf_meta),
            patient_rows(ds205, {t: act205[t] for t in common_tfs}, tf_meta),
        ],
        ignore_index=True,
    )
    focus_items = [(p, t) for p, t in FOCUS_ORDER if t in common_tfs]
    extra_focus = [("IFN", t) for t in common_tfs if t not in {x for _, x in FOCUS_ORDER}]
    tf_items = focus_items
    mod_items = [("module", k) for k in MODULES]
    tf_tests = test_table(patients, tf_items, "q4q1")
    tf_tests_med = test_table(patients, tf_items, "median")
    mod_tests = pd.concat(
        [test_table(patients, mod_items, "q4q1"), test_table(patients, mod_items, "median")],
        ignore_index=True,
    )
    all_tf_tests = pd.concat([tf_tests, tf_tests_med], ignore_index=True)

    out = ROOT / "results"
    figdir = ROOT / "figures"
    out.mkdir(parents=True, exist_ok=True)
    compact = compact_tf_table(all_tf_tests, tf_meta)
    compact.to_csv(out / "tf_table.tsv", sep="\t", index=False)
    all_tf_tests.to_csv(out / "tf_tests.tsv", sep="\t", index=False)
    mod_tests.to_csv(out / "module_tests.tsv", sep="\t", index=False)
    tf_meta.to_csv(out / "tf_meta.tsv", sep="\t", index=False)
    # per-patient: keep occupancy + focused activities, not every TF
    keep_cols = [
        c
        for c in patients.columns
        if c.split("_", 1)[0] in {"cohort", "patient", "n", "origins", "cldn4", "eligible", "histology", "recist"}
        or any(c.endswith(f"_{tf}") or c == f"q4_{tf}" for tf in [t for _, t in FOCUS_ORDER])
        or any(c.endswith(f"_{m}") for m in MODULES)
        or c.startswith(("n_", "eligible_", "cldn4_", "cohort", "patient", "origins", "histology", "recist"))
    ]
    # simpler: drop all_* for non-focus TFs
    drop = []
    focus_names = {t for _, t in FOCUS_ORDER} | set(MODULES)
    for c in patients.columns:
        tail = c.split("_", 1)[-1] if "_" in c else c
        prefix = c.split("_", 1)[0]
        if prefix in {"all", "q4", "q1", "medh", "medl"} and tail not in focus_names:
            drop.append(c)
    patients.drop(columns=drop).to_csv(out / "per_patient.tsv", sep="\t", index=False)

    inventory = {
        "gse131907_cells": ds131["n_cells_matrix"],
        "gse131907_malignant": ds131["n_malignant"],
        "gse205335_cells": ds205["n_cells_matrix"],
        "gse205335_malignant": ds205["n_malignant"],
        "n_shared_targets": len(shared),
        "n_tfs_scored": len(common_tfs),
        "n_eligible_q4q1": int(patients.eligible_q4q1.sum()),
        "excluded_accessions": ["GSE207422", "GSE148071"],
        "dual_high": False,
        "scorer": "dorothea_wmean_documented",
        "decoupler": False,
    }
    (out / "summary.json").write_text(json.dumps(inventory, indent=2) + "\n")
    figures = write_figures(patients, all_tf_tests, figdir)
    write_finding(patients, all_tf_tests, mod_tests, tf_meta, inventory, figures, ROOT / "FINDING.md")
    print(json.dumps(inventory, indent=2))
    print(compact.to_string(index=False))


if __name__ == "__main__":
    main()
