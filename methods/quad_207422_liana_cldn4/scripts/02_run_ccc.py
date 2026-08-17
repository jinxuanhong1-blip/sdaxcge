#!/usr/bin/env python3
"""QUAD merge CLDN4-only LIANA / CellPhoneDB-style LR.

Cohorts: GSE207422 + GSE131907 + GSE148071 + GSE205335.
CLDN4 only. TACSTD2 is not a gate. Dual-high is not run.
This is not a re-audit of the GSE207422-only LIANA (PR #344).

Primary: documented CellPhoneDB mean-of-means, patient-level paired Wilcoxon.
Secondary: LIANA mt.cellphonedb per cohort if importable.
CellChat is not run (no R).
"""

from __future__ import annotations

import argparse
import gzip
import json
import shutil
import tempfile
import traceback
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from scipy import sparse, stats
from statsmodels.stats.multitest import multipletests

HERE = Path(__file__).resolve().parents[1]
FOCUS = {"T_recruit", "IFN", "MHC_I"}
FOCUS_PAIRS = [
    ("CXCL9", "CXCR3"),
    ("CXCL10", "CXCR3"),
    ("CXCL11", "CXCR3"),
    ("CXCL16", "CXCR6"),
    ("CCL5", "CCR5"),
    ("CCL4", "CCR5"),
    ("CX3CL1", "CX3CR1"),
    ("HLA-A", "CD8A"),
    ("HLA-B", "CD8A"),
    ("HLA-C", "CD8A"),
    ("HLA-A", "CD8B"),
    ("HLA-B", "CD8B"),
    ("HLA-C", "CD8B"),
    ("HLA-E", "KLRD1"),
    ("HLA-E", "KLRC1"),
    ("IFNG", "IFNGR1"),
    ("IFNG", "IFNGR2"),
]

PAPER_GROUP_207422 = {
    "BD_immune01": "TN",
    "BD_immune02": "NMPR",
    "BD_immune03": "MPR",
    "BD_immune04": "NMPR",
    "BD_immune05": "TN",
    "BD_immune06": "MPR",
    "BD_immune07": "NMPR",
    "BD_immune08": "TN",
    "BD_immune09": "NMPR",
    "BD_immune10": "NMPR",
    "BD_immune11": "MPR",
    "BD_immune12": "NMPR",
    "BD_immune13": "NMPR",
    "BD_immune14": "MPR",
    "BD_immune15": "NMPR",
}
PATIENT_OF_207422 = {f"BD_immune{i:02d}": f"P{i:02d}" for i in range(1, 16)}
TUMOR_ORIGINS_131907 = {"tLung", "tL/B", "mLN", "mBrain", "PE"}
MALIGNANT_SUB_131907 = {"Malignant cells", "tS1", "tS2", "tS3"}
T_SUBS_205335 = {"CD4+ T cells", "CD8+ T cells"}
NK_SUBS_205335 = {"NK cells"}


def log(msg: str) -> None:
    print(msg, flush=True)


def log1p_cp10k(umi: np.ndarray, total: np.ndarray) -> np.ndarray:
    with np.errstate(divide="ignore", invalid="ignore"):
        cp = np.where(total > 0, umi / total * 1e4, 0.0)
    return np.log1p(cp).astype(np.float32)


def group_gene_stats(logx: dict[str, np.ndarray], mask: np.ndarray):
    n = int(mask.sum())
    means, fracs = {}, {}
    if n == 0:
        return means, fracs, 0
    for gene, arr in logx.items():
        v = arr[mask]
        means[gene] = float(np.mean(v))
        fracs[gene] = float(np.mean(v > 0))
    return means, fracs, n


def partner_from_stats(units: list[str], means: dict, fracs: dict):
    m, f = [], []
    for gene in units:
        if gene not in means:
            return np.nan, np.nan
        m.append(means[gene])
        f.append(fracs[gene])
    return float(np.min(m)), float(np.min(f))


def score_pairs(pairs: pd.DataFrame, logx: dict[str, np.ndarray], sender, receiver, expr_prop: float):
    s_mean, s_frac, n_s = group_gene_stats(logx, sender)
    r_mean, r_frac, n_r = group_gene_stats(logx, receiver)
    rows = []
    for rec in pairs.itertuples(index=False):
        lig_u = str(rec.ligand).split("+")
        rec_u = str(rec.receptor).split("+")
        l_mean, l_frac = partner_from_stats(lig_u, s_mean, s_frac)
        rec_m, rec_f = partner_from_stats(rec_u, r_mean, r_frac)
        if not np.isfinite(l_mean) or not np.isfinite(rec_m):
            continue
        rows.append(
            {
                "ligand": rec.ligand,
                "receptor": rec.receptor,
                "pathway": rec.pathway,
                "pair_origin": getattr(rec, "pair_origin", ""),
                "n_sender": n_s,
                "n_receiver": n_r,
                "ligand_mean": l_mean,
                "receptor_mean": rec_m,
                "ligand_frac": l_frac,
                "receptor_frac": rec_f,
                "cpdb_mean_score": 0.5 * (l_mean + rec_m),
                "pass_expr_prop": bool((l_frac >= expr_prop) and (rec_f >= expr_prop)),
            }
        )
    return pd.DataFrame(rows)


def paired_delta(high_df: pd.DataFrame, low_df: pd.DataFrame) -> pd.DataFrame:
    key = ["ligand", "receptor", "pathway"]
    cols = ["cpdb_mean_score", "pass_expr_prop", "ligand_frac", "receptor_frac", "n_sender", "n_receiver"]
    merged = high_df[key + cols].merge(low_df[key + cols], on=key, suffixes=("_high", "_low"))
    merged["delta_high_minus_low"] = merged["cpdb_mean_score_high"] - merged["cpdb_mean_score_low"]
    merged["pass_either"] = merged["pass_expr_prop_high"] | merged["pass_expr_prop_low"]
    merged["pass_both"] = merged["pass_expr_prop_high"] & merged["pass_expr_prop_low"]
    return merged


def wilcoxon_safe(a, b) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if len(a) < 3 or np.allclose(a, b):
        return np.nan
    try:
        return float(stats.wilcoxon(a, b, zero_method="wilcox", alternative="two-sided").pvalue)
    except ValueError:
        return np.nan


def fmt_p(x) -> str:
    if x is None or not np.isfinite(x):
        return "NA"
    return f"{x:.3g}"


def stream_umi(path: Path, keep: set[str], strip_version: bool = False):
    log(f"[stream] {path} keep={len(keep)}")
    with gzip.open(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        cells = np.array(header[1:], dtype=object)
        n = len(cells)
        totals = np.zeros(n, dtype=np.float64)
        store: dict[str, np.ndarray] = {}
        n_genes = 0
        for line in handle:
            gene, sep, rest = line.partition("\t")
            if not sep:
                continue
            if strip_version:
                gene = gene.split(".")[0]
            n_genes += 1
            vals = np.fromstring(rest, sep="\t", dtype=np.float32)
            if vals.size != n:
                raise ValueError(f"column mismatch for {gene}: {vals.size} != {n}")
            totals += vals
            if gene in keep and gene not in store:
                store[gene] = vals
            if n_genes % 5000 == 0:
                log(f"[stream] genes_seen={n_genes} kept={len(store)}")
    log(f"[stream] cells={n} genes={n_genes} kept={len(store)}")
    return cells, store, totals, n_genes


def score_lineage(expr: dict[str, np.ndarray], genes: list[str], n: int) -> np.ndarray:
    mats = [np.log1p(expr[g].astype(np.float32)) for g in genes if g in expr]
    if not mats:
        return np.zeros(n, dtype=np.float32)
    return np.mean(np.vstack(mats), axis=0)


def assign_lineage(expr: dict[str, np.ndarray], lineages: dict, n: int) -> np.ndarray:
    names = list(lineages)
    scores = np.vstack([score_lineage(expr, lineages[k], n) for k in names])
    best = scores.argmax(axis=0)
    best_val = scores.max(axis=0)
    second = np.partition(scores, -2, axis=0)[-2]
    assigned = np.array(names, dtype=object)[best]
    assigned[(best_val < 0.15) | ((best_val - second) < 0.05)] = "other"
    return assigned


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
    opener = gzip.open if str(path).endswith(".gz") else open
    records: list[dict[str, str]] = []
    current: dict[str, str] | None = None
    descriptions: list[str] = []
    titles: list[str] = []
    with opener(path, "rt", errors="replace") as handle:
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
    metadata = pd.DataFrame(records)
    read_end = metadata["platform"].str.extract(r"Single Cell ([35])'")[0]
    metadata["orig.ident"] = (
        metadata["description"].str.replace("_", "-", regex=False) + "-" + read_end + "P"
    )
    return metadata.rename(columns={"tumor stage": "tumor_stage", "cancer subtype": "cancer_subtype"})


def load_rds_genes(path: Path, wanted: set[str]):
    with tempfile.TemporaryDirectory(prefix="gse205335-") as tmp:
        matrix_path = path
        if str(path).endswith(".gz"):
            matrix_path = Path(tmp) / "GSE205335_Lung_IO_UMI_matrix.rds"
            log(f"decompress {path.name}")
            with gzip.open(path, "rb") as source, matrix_path.open("wb") as dest:
                shutil.copyfileobj(source, dest, 16 * 1024 * 1024)
        import rdata

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message='Missing constructor for R class "dgCMatrix"')
            obj = rdata.read_rds(matrix_path)
    expected = {"i", "p", "Dim", "Dimnames", "x"}
    if not expected.issubset(vars(obj)):
        raise TypeError("RDS is not the expected Matrix::dgCMatrix object")
    genes = np.asarray(obj.Dimnames[0], dtype=str)
    barcodes = np.asarray(obj.Dimnames[1], dtype=str)
    log(f"build CSC {tuple(obj.Dim)}")
    matrix = sparse.csc_matrix((obj.x, obj.i, obj.p), shape=tuple(obj.Dim), copy=False)
    library_umi = np.asarray(matrix.sum(axis=0)).ravel()
    matrix = matrix.tocsr()
    name_to_row = {g: i for i, g in enumerate(genes)}
    extracted = {}
    for gene in sorted(wanted):
        row = name_to_row.get(gene)
        if row is None:
            continue
        extracted[gene] = np.asarray(matrix.getrow(row).toarray()).ravel().astype(np.float32)
    log(f"extracted {len(extracted)} / {len(wanted)} genes")
    return extracted, library_umi, barcodes


def load_tisch_genes(path: Path, wanted: list[str]):
    import h5py

    def as_str(arr):
        return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in arr], dtype=object)

    with h5py.File(path, "r") as handle:
        grp = None
        if all(k in handle for k in ("data", "indices", "indptr")):
            grp = handle
        else:
            for key in handle.keys():
                node = handle[key]
                if hasattr(node, "keys") and all(k in node for k in ("data", "indices", "indptr")):
                    grp = node
                    break
        if grp is None:
            raise ValueError(f"no sparse matrix in {path}")

        def names(cands):
            for cand in cands:
                node = grp.get(cand) if hasattr(grp, "get") else None
                if node is None:
                    node = handle.get(cand)
                if node is None:
                    continue
                if hasattr(node, "shape") and not hasattr(node, "keys"):
                    return as_str(node[:])
                if hasattr(node, "keys"):
                    for sub in ("name", "id", "gene_names"):
                        if sub in node:
                            return as_str(node[sub][:])
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
        mat = sparse.csc_matrix((data, indices, indptr), shape=(n_cells, n_genes))

        def col(i: int):
            return np.asarray(mat[:, i].todense()).ravel()
    else:
        n_genes, n_cells = n0, n1
        mat = sparse.csc_matrix((data, indices, indptr), shape=(n_genes, n_cells))

        def col(i: int):
            return np.asarray(mat.getrow(i).todense()).ravel()

    gene_index = {g: int(np.where(genes == g)[0][0]) for g in wanted if g in set(genes.tolist())}
    store = {g: col(i).astype(np.float32) for g, i in gene_index.items()}
    log(f"[tisch] cells={n_cells} genes={n_genes} kept={len(store)}")
    return store, barcodes.astype(str), n_genes


def patient_chunks(pairs, logx, patient, cld_hi, cld_lo, is_tnk, paired_ids, expr_prop, extra_fn):
    chunks = []
    for pat in paired_ids:
        m = patient == pat
        h, l, t = m & cld_hi, m & cld_lo, m & is_tnk
        hdf = score_pairs(pairs, logx, h, t, expr_prop)
        ldf = score_pairs(pairs, logx, l, t, expr_prop)
        d = paired_delta(hdf, ldf)
        d["patient"] = pat
        d["n_high"] = int(h.sum())
        d["n_low"] = int(l.sum())
        d["n_tnk"] = int(t.sum())
        extra_fn(d, pat)
        chunks.append(d)
        log(f"  outgoing {pat} high={int(h.sum())} low={int(l.sum())} tnk={int(t.sum())}")
    return pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame()


def incoming_chunks(pairs, logx, patient, cld_hi, cld_lo, is_tnk, paired_ids, expr_prop):
    chunks = []
    for pat in paired_ids:
        m = patient == pat
        h, l, t = m & cld_hi, m & cld_lo, m & is_tnk
        hdf = score_pairs(pairs, logx, t, h, expr_prop)
        ldf = score_pairs(pairs, logx, t, l, expr_prop)
        d = paired_delta(hdf, ldf)
        d["patient"] = pat
        chunks.append(d)
    return pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame()


def rank_table(raw: pd.DataFrame) -> pd.DataFrame:
    if raw.empty:
        return pd.DataFrame(
            columns=[
                "ligand",
                "receptor",
                "pathway",
                "n_patients",
                "n_cohorts",
                "median_delta",
                "mean_delta",
                "pval",
                "padj",
            ]
        )
    rows = []
    for (lig, recp, path), sub in raw.groupby(["ligand", "receptor", "pathway"], observed=True):
        keep_sub = sub[sub["pass_either"]]
        if len(keep_sub) < 3:
            continue
        if path not in FOCUS and keep_sub["pass_both"].sum() < 3:
            continue
        p = wilcoxon_safe(keep_sub["cpdb_mean_score_high"], keep_sub["cpdb_mean_score_low"])
        rows.append(
            {
                "ligand": lig,
                "receptor": recp,
                "pathway": path,
                "n_patients": int(keep_sub["patient"].nunique()),
                "n_cohorts": int(keep_sub["cohort"].nunique()) if "cohort" in keep_sub else 1,
                "median_delta": float(np.median(keep_sub["delta_high_minus_low"])),
                "mean_delta": float(np.mean(keep_sub["delta_high_minus_low"])),
                "mean_score_high": float(keep_sub["cpdb_mean_score_high"].mean()),
                "mean_score_low": float(keep_sub["cpdb_mean_score_low"].mean()),
                "frac_pass_high": float(keep_sub["pass_expr_prop_high"].mean()),
                "frac_pass_low": float(keep_sub["pass_expr_prop_low"].mean()),
                "n_gse207422": int((keep_sub["cohort"] == "GSE207422").sum()) if "cohort" in keep_sub else 0,
                "n_gse131907": int((keep_sub["cohort"] == "GSE131907").sum()) if "cohort" in keep_sub else 0,
                "n_gse148071": int((keep_sub["cohort"] == "GSE148071").sum()) if "cohort" in keep_sub else 0,
                "n_gse205335": int((keep_sub["cohort"] == "GSE205335").sum()) if "cohort" in keep_sub else 0,
                "pval": p,
            }
        )
    tab = pd.DataFrame(rows)
    if len(tab) and tab["pval"].notna().any():
        mask = tab["pval"].notna()
        tab.loc[mask, "padj"] = multipletests(tab.loc[mask, "pval"], method="fdr_bh")[1]
    else:
        tab["padj"] = np.nan
    return tab.sort_values(["pathway", "median_delta"]) if len(tab) else tab


def run_liana_cohort(logx, cld_hi, cld_lo, is_t, is_nk, out_csv, n_perms, max_cells, seed) -> str:
    if n_perms <= 0:
        return "LIANA_SKIPPED"
    try:
        import anndata as ad
        import liana as li
    except Exception as exc:
        return f"LIANA_IMPORT_FAILED: {exc}"
    try:
        mask = cld_hi | cld_lo | is_t | is_nk
        genes = sorted(logx)
        x = np.vstack([logx[g][mask] for g in genes]).T
        obs = pd.DataFrame(index=np.arange(int(mask.sum())))
        grp = np.array(["other"] * int(mask.sum()), dtype=object)
        grp[cld_hi[mask]] = "Malig_CLDN4high"
        grp[cld_lo[mask]] = "Malig_CLDN4low"
        grp[is_t[mask]] = "T"
        grp[is_nk[mask]] = "NK"
        obs["cc_group"] = grp
        rng = np.random.default_rng(seed)
        keep_idx = []
        for _, idx in obs.groupby("cc_group", observed=True).indices.items():
            if obs.iloc[idx]["cc_group"].iloc[0] == "other":
                continue
            if len(idx) > max_cells:
                idx = rng.choice(idx, size=max_cells, replace=False)
            keep_idx.append(np.asarray(idx))
        keep_idx = np.sort(np.concatenate(keep_idx)) if keep_idx else np.array([], dtype=int)
        adata = ad.AnnData(X=x[keep_idx], obs=obs.iloc[keep_idx].copy(), var=pd.DataFrame(index=genes))
        adata.obs["cc_group"] = pd.Categorical(adata.obs["cc_group"])
        adata.uns["log1p"] = {"base": None}
        li.mt.cellphonedb(
            adata,
            groupby="cc_group",
            resource_name="cellphonedb",
            expr_prop=0.10,
            n_perms=n_perms,
            use_raw=False,
            verbose=True,
            key_added="liana_res",
        )
        res = adata.uns["liana_res"].copy()
        res.to_csv(out_csv, index=False)
        return f"LIANA_OK n_edges={len(res)}"
    except Exception as exc:
        return f"LIANA_RUN_FAILED: {exc}\n{traceback.format_exc()}"


def score_gse207422(datadir: Path, pairs: pd.DataFrame, keep: set[str], cfg: dict, outdir: Path):
    log("=== GSE207422 ===")
    P = cfg["params"]
    cells, expr, total, n_genes = stream_umi(
        datadir / "gse207422" / "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz", keep
    )
    n = len(cells)
    sample = np.array([str(c).rsplit("_", 1)[0] for c in cells])
    patient = np.array([PATIENT_OF_207422.get(s, s) for s in sample], dtype=object)
    lineage = assign_lineage(expr, cfg["lineage_markers"], n)
    normal_umi = np.zeros(n, dtype=np.float64)
    for g in cfg["normal_lung"]:
        if g in expr:
            normal_umi += expr[g]
    is_epi = lineage == "epithelial"
    is_malig = is_epi & (normal_umi == 0)
    is_t = lineage == "T"
    is_nk = lineage == "NK"
    is_tnk = is_t | is_nk
    logx = {g: log1p_cp10k(expr[g], total) for g in expr}
    thr = float(np.median(logx["CLDN4"][is_malig]))
    cld_hi = is_malig & (logx["CLDN4"] >= thr)
    cld_lo = is_malig & (logx["CLDN4"] < thr)
    rows = []
    for pat in sorted(set(patient)):
        m = patient == pat
        n_hi, n_lo, n_t = int((m & cld_hi).sum()), int((m & cld_lo).sum()), int((m & is_tnk).sum())
        samp = str(sample[m][0])
        rows.append(
            {
                "cohort": "GSE207422",
                "patient": pat,
                "note": f"{samp}/{PAPER_GROUP_207422.get(samp, '')}",
                "n_cells": int(m.sum()),
                "n_malignant": int((m & is_malig).sum()),
                "n_T": int((m & is_t).sum()),
                "n_NK": int((m & is_nk).sum()),
                "n_T_NK": n_t,
                "n_cldn4_high": n_hi,
                "n_cldn4_low": n_lo,
                "paired": n_hi >= P["min_malig_per_state"] and n_lo >= P["min_malig_per_state"] and n_t >= P["min_tnk"],
            }
        )
    ntab = pd.DataFrame(rows)
    paired = ntab.loc[ntab["paired"], "patient"].tolist()
    raw = patient_chunks(
        pairs, logx, patient, cld_hi, cld_lo, is_tnk, paired, P["expr_prop"], lambda d, p: None
    )
    inc = incoming_chunks(pairs, logx, patient, cld_hi, cld_lo, is_tnk, paired, P["expr_prop"])
    raw["cohort"] = "GSE207422"
    inc["cohort"] = "GSE207422"
    stats_d = {
        "n_cells": n,
        "n_genes": n_genes,
        "n_malignant": int(is_malig.sum()),
        "n_T": int(is_t.sum()),
        "n_NK": int(is_nk.sum()),
        "n_tnk": int(is_tnk.sum()),
        "n_cldn4_high": int(cld_hi.sum()),
        "n_cldn4_low": int(cld_lo.sum()),
        "cldn4_threshold": thr,
        "n_patients_listed": int(ntab["patient"].nunique()),
        "n_patients_paired": int(len(paired)),
        "scale": "log1p_cp10k",
        "malignant_rule": "A3 epithelial AND zero normal-lung UMI",
    }
    liana_note = run_liana_cohort(
        logx, cld_hi, cld_lo, is_t, is_nk,
        outdir / "liana_gse207422.csv",
        P["liana_n_perms"], P["liana_max_cells_per_group"], P["random_seed"],
    )
    return ntab, raw, inc, stats_d, liana_note


def score_gse131907(datadir: Path, pairs: pd.DataFrame, keep: set[str], cfg: dict, outdir: Path):
    log("=== GSE131907 ===")
    P = cfg["params"]
    cells, expr, total, n_genes = stream_umi(
        datadir / "gse131907" / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz", keep, strip_version=True
    )
    ann = pd.read_csv(
        datadir / "gse131907" / "GSE131907_Lung_Cancer_cell_annotation.txt.gz", sep="\t", dtype=str
    )
    idx_col = ann.columns[0]
    ann = ann.set_index(idx_col)
    meta = parse_series_matrix(datadir / "gse131907" / "GSE131907_series_matrix.txt.gz")
    sample_meta = meta.rename(columns={"title": "Sample", "tissue_origin_abbrevation": "Sample_Origin_geo"})
    keep_cols = [c for c in ["Sample", "geo_accession", "patient_id"] if c in sample_meta.columns]
    sample_meta = sample_meta[keep_cols].drop_duplicates("Sample")
    per = ann.reindex(cells)
    if per["Sample"].isna().any():
        raise SystemExit("GSE131907 matrix cell IDs do not align with annotation Index")
    sample = per["Sample"].to_numpy()
    origin = per["Sample_Origin"].to_numpy()
    ctype = per["Cell_type"].fillna("").to_numpy()
    csub = per["Cell_subtype"].fillna("").to_numpy()
    patient_map = sample_meta.set_index("Sample")["patient_id"].to_dict()
    patient = np.array([patient_map.get(s, s) for s in sample], dtype=object)
    is_tumor = np.isin(origin, list(TUMOR_ORIGINS_131907))
    is_malig = is_tumor & np.isin(csub, list(MALIGNANT_SUB_131907))
    is_t = is_tumor & (ctype == "T lymphocytes")
    is_nk = is_tumor & (ctype == "NK cells")
    is_tnk = is_t | is_nk
    logx = {g: log1p_cp10k(expr[g], total) for g in expr}
    thr = float(np.median(logx["CLDN4"][is_malig])) if is_malig.any() else float("nan")
    cld_hi = is_malig & (logx["CLDN4"] >= thr)
    cld_lo = is_malig & (logx["CLDN4"] < thr)
    rows = []
    for pat in sorted(set(patient)):
        m = (patient == pat) & is_tumor
        if not m.any():
            m = patient == pat
            n_hi = n_lo = n_t = 0
            n_mal = 0
            paired = False
            note = "no_tumor_origin"
        else:
            n_hi, n_lo, n_t = int((m & cld_hi).sum()), int((m & cld_lo).sum()), int((m & is_tnk).sum())
            n_mal = int((m & is_malig).sum())
            paired = n_hi >= P["min_malig_per_state"] and n_lo >= P["min_malig_per_state"] and n_t >= P["min_tnk"]
            note = ",".join(sorted(set(origin[m].astype(str))))
        rows.append(
            {
                "cohort": "GSE131907",
                "patient": pat,
                "note": note,
                "n_cells": int((patient == pat).sum()),
                "n_malignant": n_mal,
                "n_T": int((m & is_t).sum()) if m.any() else 0,
                "n_NK": int((m & is_nk).sum()) if m.any() else 0,
                "n_T_NK": n_t,
                "n_cldn4_high": n_hi,
                "n_cldn4_low": n_lo,
                "paired": paired,
            }
        )
    ntab = pd.DataFrame(rows)
    paired = ntab.loc[ntab["paired"], "patient"].tolist()
    raw = patient_chunks(
        pairs, logx, patient, cld_hi, cld_lo, is_tnk, paired, P["expr_prop"], lambda d, p: None
    )
    inc = incoming_chunks(pairs, logx, patient, cld_hi, cld_lo, is_tnk, paired, P["expr_prop"])
    raw["cohort"] = "GSE131907"
    inc["cohort"] = "GSE131907"
    stats_d = {
        "n_cells": int(len(cells)),
        "n_genes": n_genes,
        "n_malignant": int(is_malig.sum()),
        "n_T": int(is_t.sum()),
        "n_NK": int(is_nk.sum()),
        "n_tnk": int(is_tnk.sum()),
        "n_cldn4_high": int(cld_hi.sum()),
        "n_cldn4_low": int(cld_lo.sum()),
        "cldn4_threshold": thr,
        "n_patients_listed": int(sample_meta["patient_id"].nunique()),
        "n_patients_paired": int(len(paired)),
        "n_samples_geo": int(sample_meta["Sample"].nunique()),
        "scale": "log1p_cp10k",
        "malignant_rule": "tumor-origin AND Cell_subtype in {Malignant cells, tS1, tS2, tS3}",
    }
    liana_note = run_liana_cohort(
        logx, cld_hi, cld_lo, is_t, is_nk,
        outdir / "liana_gse131907.csv",
        P["liana_n_perms"], P["liana_max_cells_per_group"], P["random_seed"],
    )
    return ntab, raw, inc, stats_d, liana_note


def score_gse148071(datadir: Path, pairs: pd.DataFrame, keep: set[str], cfg: dict, outdir: Path):
    log("=== GSE148071 ===")
    P = cfg["params"]
    meta = pd.read_csv(datadir / "gse148071" / "NSCLC_GSE148071_CellMetainfo_table.tsv", sep="\t", low_memory=False)
    cell_col = meta.columns[0]
    meta[cell_col] = meta[cell_col].astype(str)
    meta = meta.set_index(cell_col)
    lin_col = next(c for c in meta.columns if "major-lineage" in c.lower())
    meta["_lineage"] = meta[lin_col].astype(str).str.strip()
    meta["_patient"] = meta["Patient"].astype(str)
    tnk_set = set(cfg["tnk_lineages"])
    malig_set = set(cfg["malignant_lineages"])
    store, barcodes, n_genes = load_tisch_genes(
        datadir / "gse148071" / "NSCLC_GSE148071_expression.h5", sorted(keep)
    )
    expr_index = pd.Index(barcodes)
    common = meta.index.intersection(expr_index)
    log(f"aligned cells={len(common)} / meta={len(meta)} / h5={len(expr_index)}")
    order = pd.Index(common)
    pos = expr_index.get_indexer(order)
    logx = {g: store[g][pos] for g in store}
    meta = meta.loc[order]
    is_malig = meta["_lineage"].isin(malig_set).to_numpy()
    is_tnk = meta["_lineage"].isin(tnk_set).to_numpy()
    is_t = meta["_lineage"].isin(tnk_set - {"NK", "NKT"}).to_numpy()
    is_nk = meta["_lineage"].isin({"NK", "NKT"}).to_numpy()
    patient = meta["_patient"].to_numpy()
    thr = float(np.median(logx["CLDN4"][is_malig]))
    cld_hi = is_malig & (logx["CLDN4"] >= thr)
    cld_lo = is_malig & (logx["CLDN4"] < thr)
    rows = []
    for pat, sub in meta.groupby("_patient", sort=True):
        loc = meta.index.get_indexer(sub.index)
        n_hi, n_lo = int(cld_hi[loc].sum()), int(cld_lo[loc].sum())
        n_t = int(is_tnk[loc].sum())
        rows.append(
            {
                "cohort": "GSE148071",
                "patient": pat,
                "note": "TISCH2",
                "n_cells": int(len(sub)),
                "n_malignant": int(is_malig[loc].sum()),
                "n_T": int(is_t[loc].sum()),
                "n_NK": int(is_nk[loc].sum()),
                "n_T_NK": n_t,
                "n_cldn4_high": n_hi,
                "n_cldn4_low": n_lo,
                "paired": n_hi >= P["min_malig_per_state"] and n_lo >= P["min_malig_per_state"] and n_t >= P["min_tnk"],
            }
        )
    ntab = pd.DataFrame(rows)
    paired = ntab.loc[ntab["paired"], "patient"].tolist()
    raw = patient_chunks(
        pairs, logx, patient, cld_hi, cld_lo, is_tnk, paired, P["expr_prop"], lambda d, p: None
    )
    inc = incoming_chunks(pairs, logx, patient, cld_hi, cld_lo, is_tnk, paired, P["expr_prop"])
    raw["cohort"] = "GSE148071"
    inc["cohort"] = "GSE148071"
    stats_d = {
        "n_cells": int(len(meta)),
        "n_genes": n_genes,
        "n_malignant": int(is_malig.sum()),
        "n_T": int(is_t.sum()),
        "n_NK": int(is_nk.sum()),
        "n_tnk": int(is_tnk.sum()),
        "n_cldn4_high": int(cld_hi.sum()),
        "n_cldn4_low": int(cld_lo.sum()),
        "cldn4_threshold": thr,
        "n_patients_listed": int(ntab["patient"].nunique()),
        "n_patients_paired": int(len(paired)),
        "scale": "tisch2_log2_tpm10p1",
        "malignant_rule": "TISCH2 major-lineage Malignant",
        "tnk_note": "TISCH2 T/NK labels; NK often 0 in this object",
    }
    liana_note = run_liana_cohort(
        logx, cld_hi, cld_lo, is_t, is_nk,
        outdir / "liana_gse148071.csv",
        P["liana_n_perms"], P["liana_max_cells_per_group"], P["random_seed"],
    )
    return ntab, raw, inc, stats_d, liana_note


def score_gse205335(datadir: Path, pairs: pd.DataFrame, keep: set[str], cfg: dict, outdir: Path):
    log("=== GSE205335 ===")
    P = cfg["params"]
    identities = pd.read_csv(datadir / "gse205335" / "GSE205335_Lung_IO_CellIdentity.txt.gz", sep="\t")
    metadata = parse_geo_soft(datadir / "gse205335" / "GSE205335_family.soft.gz")
    extracted, library_umi, barcodes = load_rds_genes(
        datadir / "gse205335" / "GSE205335_Lung_IO_UMI_matrix.rds.gz", keep
    )
    indexed = identities.set_index("barcode")
    missing = pd.Index(barcodes).difference(indexed.index)
    extra = indexed.index.difference(pd.Index(barcodes))
    if len(missing) or len(extra):
        raise SystemExit(f"GSE205335 matrix/identity mismatch: {len(missing)} missing, {len(extra)} extra")
    cells = indexed.loc[barcodes].reset_index()
    cells = cells.merge(
        metadata[["orig.ident", "gsm", "patient", "tissue", "recist", "cancer_subtype"]],
        on="orig.ident",
        how="left",
        validate="many_to_one",
    )
    if cells["patient"].isna().any():
        raise SystemExit("GSE205335 identity samples did not match GEO metadata")
    cells["is_tumor_sample"] = ~cells["tissue"].astype(str).str.startswith("Normal")
    cells["is_malig"] = cells["is_tumor_sample"] & cells["lineage.sub"].eq("Malignant cells")
    cells["is_t"] = cells["is_tumor_sample"] & cells["lineage.sub"].isin(T_SUBS_205335)
    cells["is_nk"] = cells["is_tumor_sample"] & cells["lineage.sub"].isin(NK_SUBS_205335)
    cells["is_tnk"] = cells["is_tumor_sample"] & cells["lineage.total"].eq("T/NK cells")
    is_malig = cells["is_malig"].to_numpy()
    is_t = cells["is_t"].to_numpy()
    is_nk = cells["is_nk"].to_numpy()
    is_tnk = cells["is_tnk"].to_numpy()
    patient = cells["patient"].to_numpy()
    logx = {g: log1p_cp10k(extracted[g], library_umi) for g in extracted}
    thr = float(np.median(logx["CLDN4"][is_malig])) if is_malig.any() else float("nan")
    cld_hi = is_malig & (logx["CLDN4"] >= thr)
    cld_lo = is_malig & (logx["CLDN4"] < thr)
    rows = []
    for pat, sub in cells.groupby("patient", observed=True):
        tumor = sub[sub["is_tumor_sample"]]
        idx = tumor.index.to_numpy()
        n_hi = int(cld_hi[idx].sum()) if len(tumor) else 0
        n_lo = int(cld_lo[idx].sum()) if len(tumor) else 0
        n_t = int(tumor["is_tnk"].sum()) if len(tumor) else 0
        rows.append(
            {
                "cohort": "GSE205335",
                "patient": pat,
                "note": f"{sub.iloc[0]['cancer_subtype']}/{sub.iloc[0]['recist']}",
                "n_cells": int(len(tumor)),
                "n_malignant": int(tumor["is_malig"].sum()) if len(tumor) else 0,
                "n_T": int(tumor["is_t"].sum()) if len(tumor) else 0,
                "n_NK": int(tumor["is_nk"].sum()) if len(tumor) else 0,
                "n_T_NK": n_t,
                "n_cldn4_high": n_hi,
                "n_cldn4_low": n_lo,
                "paired": bool(len(tumor))
                and n_hi >= P["min_malig_per_state"]
                and n_lo >= P["min_malig_per_state"]
                and n_t >= P["min_tnk"],
            }
        )
    ntab = pd.DataFrame(rows)
    paired = ntab.loc[ntab["paired"], "patient"].tolist()
    raw = patient_chunks(
        pairs, logx, patient, cld_hi, cld_lo, is_tnk, paired, P["expr_prop"], lambda d, p: None
    )
    inc = incoming_chunks(pairs, logx, patient, cld_hi, cld_lo, is_tnk, paired, P["expr_prop"])
    raw["cohort"] = "GSE205335"
    inc["cohort"] = "GSE205335"
    stats_d = {
        "n_cells": int(len(cells)),
        "n_genes": "rds_full",
        "n_malignant": int(is_malig.sum()),
        "n_T": int(is_t.sum()),
        "n_NK": int(is_nk.sum()),
        "n_tnk": int(is_tnk.sum()),
        "n_cldn4_high": int(cld_hi.sum()),
        "n_cldn4_low": int(cld_lo.sum()),
        "cldn4_threshold": thr,
        "n_patients_listed": int(ntab["patient"].nunique()),
        "n_patients_paired": int(len(paired)),
        "scale": "log1p_cp10k",
        "malignant_rule": "tumor sample AND lineage.sub == Malignant cells",
    }
    liana_note = run_liana_cohort(
        logx, cld_hi, cld_lo, is_t, is_nk,
        outdir / "liana_gse205335.csv",
        P["liana_n_perms"], P["liana_max_cells_per_group"], P["random_seed"],
    )
    return ntab, raw, inc, stats_d, liana_note


def make_figures(figdir: Path, ntab: pd.DataFrame, raw: pd.DataFrame, ranks: pd.DataFrame, inc_ranks: pd.DataFrame):
    figdir.mkdir(parents=True, exist_ok=True)
    colors = {"T_recruit": "#4C72B0", "IFN": "#55A868", "MHC_I": "#C44E52"}
    cohort_colors = {
        "GSE207422": "#4C72B0",
        "GSE131907": "#DD8452",
        "GSE148071": "#55A868",
        "GSE205335": "#C44E52",
    }

    show = ntab.sort_values(["cohort", "patient"]).reset_index(drop=True)
    fig, axes = plt.subplots(2, 2, figsize=(11.2, 7.2), sharey=False)
    for ax, (cohort, sub) in zip(axes.ravel(), show.groupby("cohort", sort=True)):
        x = np.arange(len(sub))
        ax.bar(x - 0.2, sub["n_malignant"], width=0.4, color="#4C72B0", label="malignant")
        ax.bar(x + 0.2, sub["n_T_NK"], width=0.4, color="#DD8452", label="T/NK")
        ax.set_xticks(x)
        ax.set_xticklabels(sub["patient"], rotation=80, ha="right", fontsize=6)
        ax.set_title(f"{cohort}  paired={int(sub['paired'].sum())}/{len(sub)}")
        ax.legend(frameon=False, fontsize=7)
    fig.suptitle("Honest n: malignant vs T/NK cells per patient (cells are not the test n)")
    fig.tight_layout()
    fig.savefig(figdir / "n_cells_by_patient.png", dpi=160)
    fig.savefig(figdir / "n_cells_by_patient.pdf")
    plt.close(fig)

    foc = ranks[ranks["pathway"].isin(FOCUS)].copy() if len(ranks) else ranks
    if len(foc):
        foc = foc.sort_values("median_delta")
        fig, ax = plt.subplots(figsize=(8.2, max(3.4, 0.28 * len(foc) + 1.2)))
        y = np.arange(len(foc))
        ax.barh(y, foc["median_delta"], color=[colors.get(p, "#999") for p in foc["pathway"]], edgecolor="none")
        ax.axvline(0, color="k", lw=0.8)
        ax.set_yticks(y)
        ax.set_yticklabels([f"{a}–{b} ({p})" for a, b, p in zip(foc["ligand"], foc["receptor"], foc["pathway"])], fontsize=7)
        ax.set_xlabel("median patient Δ (CLDN4-high − low)")
        ax.set_title(f"Outgoing malignant → T/NK  paired n={int(ntab['paired'].sum())}")
        fig.tight_layout()
        fig.savefig(figdir / "cldn4_outgoing.png", dpi=160)
        fig.savefig(figdir / "cldn4_outgoing.pdf")
        plt.close(fig)

    ifc = inc_ranks[inc_ranks["pathway"].isin(FOCUS)].copy() if len(inc_ranks) else inc_ranks
    if len(ifc):
        ifc = ifc.sort_values("median_delta")
        fig, ax = plt.subplots(figsize=(8.2, max(3.0, 0.28 * len(ifc) + 1.0)))
        y = np.arange(len(ifc))
        ax.barh(y, ifc["median_delta"], color=[colors.get(p, "#999") for p in ifc["pathway"]], edgecolor="none")
        ax.axvline(0, color="k", lw=0.8)
        ax.set_yticks(y)
        ax.set_yticklabels([f"{a}–{b}" for a, b in zip(ifc["ligand"], ifc["receptor"])], fontsize=7)
        ax.set_xlabel("median patient Δ (CLDN4-high − low)")
        ax.set_title("Incoming T/NK → malignant CLDN4 state")
        fig.tight_layout()
        fig.savefig(figdir / "cldn4_incoming.png", dpi=160)
        fig.savefig(figdir / "cldn4_incoming.pdf")
        plt.close(fig)

    # Extra 1: honest n by cohort
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    cohorts = ["GSE207422", "GSE131907", "GSE148071", "GSE205335"]
    listed = [int((ntab["cohort"] == c).sum()) for c in cohorts]
    paired = [int(((ntab["cohort"] == c) & ntab["paired"]).sum()) for c in cohorts]
    x = np.arange(len(cohorts))
    ax.bar(x - 0.18, listed, width=0.36, color="#9B9B9B", label="listed patients")
    ax.bar(x + 0.18, paired, width=0.36, color="#4C72B0", label="paired (test n)")
    ax.set_xticks(x)
    ax.set_xticklabels(cohorts, rotation=15)
    ax.set_ylabel("patients")
    ax.legend(frameon=False)
    ax.set_title("Extra: honest patient n by cohort")
    fig.tight_layout()
    fig.savefig(figdir / "extra_n_by_cohort.png", dpi=160)
    fig.savefig(figdir / "extra_n_by_cohort.pdf")
    plt.close(fig)

    # Extra 2: high vs low malignant counts
    fig, ax = plt.subplots(figsize=(10.5, 3.8))
    paired_tab = ntab[ntab["paired"]].sort_values(["cohort", "patient"])
    x = np.arange(len(paired_tab))
    ax.bar(x - 0.18, paired_tab["n_cldn4_high"], width=0.36, color="#C44E52", label="CLDN4-high")
    ax.bar(x + 0.18, paired_tab["n_cldn4_low"], width=0.36, color="#4C72B0", label="CLDN4-low")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{c}:{p}" for c, p in zip(paired_tab["cohort"], paired_tab["patient"])], rotation=80, ha="right", fontsize=6)
    ax.legend(frameon=False)
    ax.set_title("Extra: CLDN4-high vs low malignant cells in paired patients")
    fig.tight_layout()
    fig.savefig(figdir / "extra_high_low_counts.png", dpi=160)
    fig.savefig(figdir / "extra_high_low_counts.pdf")
    plt.close(fig)

    # Extra 3: cohort-stratified forest for focus pairs
    if len(raw):
        rows = []
        for lig, recp in FOCUS_PAIRS:
            sub = raw[(raw["ligand"] == lig) & (raw["receptor"] == recp) & raw["pass_either"]]
            for cohort, ss in sub.groupby("cohort"):
                if len(ss) < 2:
                    continue
                rows.append(
                    {
                        "pair": f"{lig}–{recp}",
                        "cohort": cohort,
                        "n": int(ss["patient"].nunique()),
                        "median_delta": float(ss["delta_high_minus_low"].median()),
                    }
                )
        forest = pd.DataFrame(rows)
        if len(forest):
            fig, ax = plt.subplots(figsize=(8.8, max(3.4, 0.22 * len(forest) + 1.4)))
            y = 0
            yticks, ylabels = [], []
            for pair, sub in forest.groupby("pair", sort=False):
                for rec in sub.itertuples(index=False):
                    ax.plot(rec.median_delta, y, "o", color=cohort_colors.get(rec.cohort, "#333"))
                    yticks.append(y)
                    ylabels.append(f"{pair}  {rec.cohort} n={rec.n}")
                    y += 1
                y += 0.3
            ax.axvline(0, color="k", lw=0.8)
            ax.set_yticks(yticks)
            ax.set_yticklabels(ylabels, fontsize=6)
            ax.set_xlabel("median Δ within cohort")
            ax.set_title("Extra: cohort-stratified outgoing Δ (not a meta-analysis weight)")
            fig.tight_layout()
            fig.savefig(figdir / "extra_cohort_forest.png", dpi=160)
            fig.savefig(figdir / "extra_cohort_forest.pdf")
            plt.close(fig)

    # Extra 4: per-patient heatmap of MHC-I / T-recruit Δ
    if len(raw):
        focus_raw = raw[raw["pathway"].isin(FOCUS) & raw["pass_either"]].copy()
        focus_raw["pair"] = focus_raw["ligand"] + "–" + focus_raw["receptor"]
        if len(focus_raw):
            keep_pairs = [f"{a}–{b}" for a, b in FOCUS_PAIRS if ((focus_raw["ligand"] == a) & (focus_raw["receptor"] == b)).any()]
            mat = focus_raw[focus_raw["pair"].isin(keep_pairs)].pivot_table(
                index=["cohort", "patient"], columns="pair", values="delta_high_minus_low", aggfunc="median"
            )
            if mat.shape[0] and mat.shape[1]:
                fig, ax = plt.subplots(figsize=(max(7.5, 0.55 * mat.shape[1] + 3), max(4.2, 0.22 * mat.shape[0] + 1.5)))
                im = ax.imshow(mat.to_numpy(), aspect="auto", cmap="coolwarm", vmin=-0.3, vmax=0.3)
                ax.set_xticks(np.arange(mat.shape[1]))
                ax.set_xticklabels(mat.columns, rotation=70, ha="right", fontsize=7)
                ax.set_yticks(np.arange(mat.shape[0]))
                ax.set_yticklabels([f"{a}:{b}" for a, b in mat.index], fontsize=6)
                fig.colorbar(im, ax=ax, shrink=0.6, label="Δ")
                ax.set_title("Extra: per-patient outgoing Δ")
                fig.tight_layout()
                fig.savefig(figdir / "extra_patient_focus_heatmap.png", dpi=160)
                fig.savefig(figdir / "extra_patient_focus_heatmap.pdf")
                plt.close(fig)

    # Extra 5: paired strips for HLA-B–CD8A and CXCL16–CXCR6
    if len(raw):
        fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.6))
        for ax, (lig, recp) in zip(axes, [("HLA-B", "CD8A"), ("CXCL16", "CXCR6")]):
            sub = raw[(raw["ligand"] == lig) & (raw["receptor"] == recp)]
            if sub.empty:
                ax.set_title(f"{lig}–{recp} (empty)")
                continue
            for i, (cohort, ss) in enumerate(sub.groupby("cohort")):
                ax.scatter(
                    np.zeros(len(ss)) + i * 0.15,
                    ss["delta_high_minus_low"],
                    s=22,
                    color=cohort_colors.get(cohort, "#333"),
                    label=f"{cohort} n={ss['patient'].nunique()}",
                    alpha=0.85,
                )
            ax.axhline(0, color="k", lw=0.7)
            ax.set_xticks([])
            ax.set_ylabel("patient Δ")
            ax.set_title(f"{lig}–{recp}")
            ax.legend(frameon=False, fontsize=6)
        fig.suptitle("Extra: paired patient Δ, selected outgoing pairs", fontsize=10)
        fig.tight_layout()
        fig.savefig(figdir / "extra_paired_strips.png", dpi=160)
        fig.savefig(figdir / "extra_paired_strips.pdf")
        plt.close(fig)


def write_finding(root: Path, outdir: Path, summary: dict, ntab: pd.DataFrame, ranks: pd.DataFrame, inc: pd.DataFrame):
    paired = ntab[ntab["paired"]].copy()
    dropped = ntab[~ntab["paired"]].copy()
    foc = ranks[ranks["pathway"].isin(FOCUS)].copy() if len(ranks) else ranks
    inc_f = inc[inc["pathway"].isin(FOCUS)].copy() if len(inc) else inc

    def table(df: pd.DataFrame) -> list[str]:
        if df is None or df.empty:
            return ["No T-recruit / IFN / MHC-I pairs with ≥3 paired patients.", ""]
        df = df.sort_values(["pathway", "median_delta"])
        lines = [
            "| Pathway | Pair | n patients | n cohorts | median Δ | Wilcoxon p | FDR |",
            "|---|---|---:|---:|---:|---:|---:|",
        ]
        for r in df.itertuples(index=False):
            lines.append(
                f"| {r.pathway} | {r.ligand}–{r.receptor} | {int(r.n_patients)} | {int(r.n_cohorts)} | "
                f"{r.median_delta:+.3f} | {fmt_p(r.pval)} | {fmt_p(r.padj)} |"
            )
        n_neg = int((df["median_delta"] < 0).sum())
        n_sig = int((df["padj"] < 0.05).sum()) if df["padj"].notna().any() else 0
        lines += [
            "",
            f"{n_neg}/{len(df)} focus pairs have median Δ < 0 (weaker from/to CLDN4-high). "
            f"**{n_sig}/{len(df)} reach FDR < 0.05**.",
            "",
        ]
        return lines

    paired_list = "; ".join(
        f"{r.cohort}:{r.patient} ({r.note}; high={int(r.n_cldn4_high)}/low={int(r.n_cldn4_low)}, T/NK={int(r.n_T_NK)})"
        for r in paired.itertuples(index=False)
    )
    dropped_list = "; ".join(
        f"{r.cohort}:{r.patient} (malig={int(r.n_malignant)}, high={int(r.n_cldn4_high)}, low={int(r.n_cldn4_low)}, T/NK={int(r.n_T_NK)})"
        for r in dropped.itertuples(index=False)
    )
    by = (
        foc.groupby("pathway")
        .apply(
            lambda s: f"{int((s.median_delta < 0).sum())}/{len(s)} Δ<0, median Δ={s.median_delta.median():+.3f}",
            include_groups=False,
        )
        .to_dict()
        if len(foc)
        else {}
    )
    n_neg = int((foc["median_delta"] < 0).sum()) if len(foc) else 0
    n_sig = int((foc["padj"] < 0.05).sum()) if len(foc) and foc["padj"].notna().any() else 0
    n_pair = int(summary["n_patients_paired"])
    if n_pair == 0:
        verdict = "No patient met the paired gate. The LR table is empty for inference."
    else:
        verdict = (
            f"On the patient-level CellPhoneDB-style score (paired n={n_pair} across four public objects), "
            f"CLDN4-high vs CLDN4-low outgoing T-recruit/IFN/MHC-I pairs: "
            f"{n_neg}/{len(foc)} have median Δ < 0; {n_sig}/{len(foc)} reach FDR < 0.05. "
            f"By axis: {by}. This is a ranked co-expression list, not a recruitment mechanism."
        )

    cstats = summary["cohorts"]
    lines = [
        "# QUAD merge — LIANA/CellPhoneDB-style LR from CLDN4-high malignant to T/NK",
        "",
        f"**Verdict (paired n={n_pair} honest):** {verdict}",
        "",
        "ADDITIVE. **CLDN4 only.** TACSTD2 is not a gate. Dual-high is not run. "
        "GSE207422 is included because the user asked for it. This is **not** a re-audit of "
        "the GSE207422-only LIANA (PR #344).",
        "",
        "## Question",
        "",
        "Do CLDN4-high malignant cells show weaker outgoing T-recruit / MHC-I communication "
        "toward same-patient T/NK than CLDN4-low cells in the four-cohort public merge?",
        "",
        "## What was run",
        "",
        "| Method | Status |",
        "|---|---|",
        f"| Documented CellPhoneDB-style score (mean of partner means; {summary['n_pairs_scored']:,} pairs) | **primary** — patient-level paired Wilcoxon |",
        f"| LIANA `mt.cellphonedb` (per cohort, ≤1,500 cells/group, 20 permutations) | {summary['liana_status']} |",
        "| CellChat | Not run. R unavailable. |",
        "| Joint Harmony / concatenated expression LIANA | Not run. Platforms differ. |",
        "",
        "## Honest n",
        "",
        "The test unit is the **patient**. Cells are not n. CLDN4-high/low is the "
        "within-cohort malignant median. A patient enters if it has ≥10 malignant cells "
        "in both bins and ≥20 T/NK after pooling that patient’s tumor samples.",
        "",
        "| Item | GSE207422 | GSE131907 | GSE148071 | GSE205335 | Merge |",
        "|---|---:|---:|---:|---:|---:|",
        f"| Cells | {cstats['GSE207422']['n_cells']:,} | {cstats['GSE131907']['n_cells']:,} | {cstats['GSE148071']['n_cells']:,} | {cstats['GSE205335']['n_cells']:,} | {summary['n_cells_total']:,} |",
        f"| Malignant | {cstats['GSE207422']['n_malignant']:,} | {cstats['GSE131907']['n_malignant']:,} | {cstats['GSE148071']['n_malignant']:,} | {cstats['GSE205335']['n_malignant']:,} | {summary['n_malignant_total']:,} |",
        f"| T / NK | {cstats['GSE207422']['n_T']:,} / {cstats['GSE207422']['n_NK']:,} | {cstats['GSE131907']['n_T']:,} / {cstats['GSE131907']['n_NK']:,} | {cstats['GSE148071']['n_T']:,} / {cstats['GSE148071']['n_NK']:,} | {cstats['GSE205335']['n_T']:,} / {cstats['GSE205335']['n_NK']:,} | — |",
        f"| CLDN4-high / low | {cstats['GSE207422']['n_cldn4_high']:,} / {cstats['GSE207422']['n_cldn4_low']:,} | {cstats['GSE131907']['n_cldn4_high']:,} / {cstats['GSE131907']['n_cldn4_low']:,} | {cstats['GSE148071']['n_cldn4_high']:,} / {cstats['GSE148071']['n_cldn4_low']:,} | {cstats['GSE205335']['n_cldn4_high']:,} / {cstats['GSE205335']['n_cldn4_low']:,} | — |",
        f"| Listed patients | {cstats['GSE207422']['n_patients_listed']} | {cstats['GSE131907']['n_patients_listed']} | {cstats['GSE148071']['n_patients_listed']} | {cstats['GSE205335']['n_patients_listed']} | {len(ntab)} |",
        f"| **Paired patients (test n)** | **{cstats['GSE207422']['n_patients_paired']}** | **{cstats['GSE131907']['n_patients_paired']}** | **{cstats['GSE148071']['n_patients_paired']}** | **{cstats['GSE205335']['n_patients_paired']}** | **{n_pair}** |",
        f"| Scale | log1p CP10k | log1p CP10k | TISCH2 log2(TPM/10+1) | log1p CP10k | deltas stacked |",
        "",
        f"**Paired patients (n={n_pair}):** {paired_list}.",
        "",
        f"**Dropped (not in the Wilcoxon header n):** {dropped_list}.",
        "",
        "Per-patient counts: `results/n_cells_patients.tsv`.",
        "",
        "## Primary LR table — outgoing CLDN4-high malignant → T/NK",
        "",
        "Median Δ = high − low. Negative = weaker from CLDN4-high. "
        "Full table: `results/lr_table_cldn4_outgoing_tnk.tsv`.",
        "",
    ]
    lines += table(foc)
    lines += [
        "Incoming T/NK → malignant (focus):",
        "",
    ]
    lines += table(inc_f)
    lines += [
        "## Method (short)",
        "",
        "Partner expression = min of subunit means. Pair score = mean of the two partner means "
        "(Efremova 2020 / Garcia-Alonso 2022). `pass_expr_prop` = both partners in ≥10% of cells. "
        "Wilcoxon is paired across patients. FDR is BH within the contrast. "
        "This is not a CellChat probability and does not observe secretion or spatial contact.",
        "",
        f"LIANA notes: {summary['liana_notes']}",
        "",
        "## What this is not",
        "",
        "- Not a re-audit of PR #344 (GSE207422-only LIANA).",
        "- Not dual-high / TACSTD2-gated.",
        "- Not CellChat. Not inferCNV/CopyKAT re-calls.",
        "- Not ICI/MPR as the test (GSE131907 is treatment-naive; GSE205335 has RECIST, not MPR).",
        "- Not a joint embedding. Expression was not concatenated.",
        "- Cells are not the sample size.",
        "",
        "## Extra figures",
        "",
        "| File | Content |",
        "|---|---|",
        "| `results/figures/n_cells_by_patient.png` | Malignant vs T/NK per patient, by cohort |",
        "| `results/figures/cldn4_outgoing.png` | Focus outgoing Δ |",
        "| `results/figures/cldn4_incoming.png` | Focus incoming Δ |",
        "| `results/figures/extra_n_by_cohort.png` | Listed vs paired patient n |",
        "| `results/figures/extra_high_low_counts.png` | High vs low malignant n |",
        "| `results/figures/extra_cohort_forest.png` | Cohort-stratified Δ |",
        "| `results/figures/extra_patient_focus_heatmap.png` | Per-patient focus Δ |",
        "| `results/figures/extra_paired_strips.png` | HLA-B–CD8A and CXCL16–CXCR6 |",
        "",
        "## Reproduce",
        "",
        "```bash",
        "cd methods/quad_207422_liana_cldn4",
        "python3 scripts/00_download.py",
        "python3 scripts/01_build_pairs.py",
        "python3 scripts/02_run_ccc.py",
        "```",
        "",
    ]
    text = "\n".join(lines) + "\n"
    (root / "FINDING.md").write_text(text)
    (outdir / "FINDING.md").write_text(text)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--datadir", type=Path, default=Path("/tmp/quad_207422_liana"))
    parser.add_argument("--outdir", type=Path, default=HERE / "results")
    parser.add_argument("--skip-liana", action="store_true")
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "figures").mkdir(parents=True, exist_ok=True)

    cfg = yaml.safe_load((HERE / "config" / "gene_sets.yaml").read_text())
    if args.skip_liana:
        cfg["params"]["liana_n_perms"] = 0
    pairs = pd.read_csv(HERE / "resources" / "cellphonedb_v5_lr_pairs.tsv", sep="\t")
    keep = set(Path(HERE / "resources" / "lr_genes.txt").read_text().split())
    keep.update(cfg["state_genes"])
    keep.update(cfg["normal_lung"])
    for block in cfg["lineage_markers"].values():
        keep.update(block)
    keep.update(["TACSTD2", "EPCAM", "PTPRC", "CD3D", "CD3E", "CD8A", "NKG7", "CLDN4"])

    ntabs, raws, incs, cohort_stats, liana_notes = [], [], [], {}, []
    for name, fn in [
        ("GSE207422", score_gse207422),
        ("GSE131907", score_gse131907),
        ("GSE148071", score_gse148071),
        ("GSE205335", score_gse205335),
    ]:
        ntab, raw, inc, stats_d, liana_note = fn(args.datadir, pairs, keep, cfg, args.outdir)
        ntabs.append(ntab)
        raws.append(raw)
        incs.append(inc)
        cohort_stats[name] = stats_d
        liana_notes.append(f"{name}: {liana_note.splitlines()[0]}")
        ntab.to_csv(args.outdir / f"n_cells_{name.lower()}.tsv", sep="\t", index=False)
        if len(raw):
            raw.to_csv(args.outdir / f"patient_outgoing_{name.lower()}.tsv", sep="\t", index=False)
        log(f"[{name}] paired={stats_d['n_patients_paired']} liana={liana_note.splitlines()[0]}")

    ntab = pd.concat(ntabs, ignore_index=True)
    raw = pd.concat([d for d in raws if len(d)], ignore_index=True) if any(len(d) for d in raws) else pd.DataFrame()
    inc = pd.concat([d for d in incs if len(d)], ignore_index=True) if any(len(d) for d in incs) else pd.DataFrame()
    ntab.to_csv(args.outdir / "n_cells_patients.tsv", sep="\t", index=False)
    if len(raw):
        raw.to_csv(args.outdir / "patient_cldn4_outgoing.tsv", sep="\t", index=False)
    if len(inc):
        inc.to_csv(args.outdir / "patient_cldn4_incoming.tsv", sep="\t", index=False)

    ranks = rank_table(raw)
    inc_ranks = rank_table(inc)
    out_tab = ranks.copy()
    out_tab.insert(0, "direction", "malignant_CLDN4high_to_TNK")
    out_tab.to_csv(args.outdir / "lr_table_cldn4_outgoing_tnk.tsv", sep="\t", index=False)
    out_tab.to_csv(args.outdir / "lr_table.tsv", sep="\t", index=False)
    in_tab = inc_ranks.copy()
    in_tab.insert(0, "direction", "TNK_to_malignant_CLDN4")
    in_tab.to_csv(args.outdir / "lr_table_cldn4_incoming_tnk.tsv", sep="\t", index=False)
    focus = out_tab[out_tab["pathway"].isin(FOCUS)].copy()
    focus.to_csv(args.outdir / "lr_table_focus_outgoing.tsv", sep="\t", index=False)
    ranks.to_csv(args.outdir / "ranks_cldn4_outgoing.tsv", sep="\t", index=False)
    inc_ranks.to_csv(args.outdir / "ranks_cldn4_incoming.tsv", sep="\t", index=False)

    for cohort, sub in (raw.groupby("cohort") if len(raw) else []):
        rank_table(sub).to_csv(args.outdir / f"ranks_outgoing_{cohort.lower()}.tsv", sep="\t", index=False)

    path_rows = []
    for key, tab in [("outgoing", ranks), ("incoming", inc_ranks)]:
        for path, sub in tab.groupby("pathway") if len(tab) else []:
            if path not in FOCUS:
                continue
            path_rows.append(
                {
                    "contrast": key,
                    "pathway": path,
                    "n_pairs": int(len(sub)),
                    "n_pairs_delta_neg": int((sub["median_delta"] < 0).sum()),
                    "median_of_pair_deltas": float(sub["median_delta"].median()),
                    "n_pairs_fdr05": int((sub["padj"] < 0.05).sum()) if sub["padj"].notna().any() else 0,
                }
            )
    pd.DataFrame(path_rows).to_csv(args.outdir / "pathway_summary.tsv", sep="\t", index=False)

    liana_status = "; ".join(liana_notes)
    if all("LIANA_IMPORT_FAILED" in n for n in liana_notes):
        liana_short = "not_run_import_failed"
    elif any(n.split(": ", 1)[-1].startswith("LIANA_OK") for n in liana_notes):
        liana_short = "ran_cellphonedb_method_per_cohort"
    else:
        liana_short = "attempted"

    summary = {
        "merge": "GSE207422+GSE131907+GSE148071+GSE205335",
        "state_gene": "CLDN4",
        "tacstd2_used_as_gate": False,
        "dual_high": False,
        "not_a_reaudit_of": "PR #344 GSE207422-only LIANA",
        "n_cells_total": int(sum(c["n_cells"] if isinstance(c["n_cells"], int) else 0 for c in cohort_stats.values())),
        "n_malignant_total": int(sum(c["n_malignant"] for c in cohort_stats.values())),
        "n_patients_paired": int(ntab["paired"].sum()),
        "n_pairs_scored": int(len(pairs)),
        "min_malig_per_state": cfg["params"]["min_malig_per_state"],
        "min_tnk": cfg["params"]["min_tnk"],
        "liana_status": liana_short,
        "liana_notes": liana_status,
        "cellchat_status": "not_run_R_unavailable",
        "method": "cellphonedb_mean_of_means_patient_merge",
        "unit_of_inference": "patient",
        "cohorts": cohort_stats,
    }
    (args.outdir / "summary.json").write_text(json.dumps(summary, indent=2))
    make_figures(args.outdir / "figures", ntab, raw, ranks, inc_ranks)
    write_finding(HERE, args.outdir, summary, ntab, ranks, inc_ranks)
    log(json.dumps({k: summary[k] for k in ("n_patients_paired", "liana_status")}, indent=2))
    log(f"[done] {args.outdir}")


if __name__ == "__main__":
    main()
