#!/usr/bin/env python3
"""CLDN4-only high-end LR on the Harmony n=73 patient table.

Given (not re-audited): malignant-like CLDN4 vs T/NK ρ=−0.27, n=73
(GSE131907 + GSE253013 + GSE148071 + GSE127465).

Joint Harmony embedding is not stored in the repo. The GSE253013 RDS
(~9 GB) is not downloaded. LR is run on each contributing cohort that
has a public processed matrix, then patient deltas are meta-analyzed.

Patient is the unit. CLDN4 only — TACSTD2 is not a gate.
"""
from __future__ import annotations

import gzip
import json
import tarfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
DATA = Path("/tmp/harmony73_data")
OUT = ROOT / "results"
FIG = OUT / "figures"

KH = 0.5
TRIM = 0.10
EXPR_PROP = 0.10
MIN_MAL_ARM = 8
MIN_TNK = 20
SEED = 1

GSE131907_TUMOR_ORIGINS = ("tLung", "tL/B", "mLN", "PE", "mBrain")

LINEAGES = {
    "epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "CDH1"],
    "T": ["CD3D", "CD3E", "CD2", "CD8A", "CD4"],
    "NK": ["NKG7", "GNLY", "FGFBP2", "KLRD1", "KLRF1"],
    "B": ["CD79A", "MS4A1", "CD19", "CD79B"],
    "myeloid": ["LYZ", "CD68", "CD14", "FCGR3A"],
    "fibroblast": ["COL1A1", "COL1A2", "DCN"],
    "endothelial": ["VWF", "PECAM1", "RAMP2"],
}
NORMAL_LUNG = ["SFTPA2", "SFTPC", "AGER", "SCGB1A1", "SCGB3A1", "TPPP3", "FOXJ1"]

# Curated outgoing Mal → T/NK pairs shown in the FINDING.
FOCUS = [
    ("NECTIN2", "TIGIT", "barrier"),
    ("PVR", "TIGIT", "barrier"),
    ("CDH1", "ITGAE+ITGB7", "barrier"),
    ("F11R", "ITGAL+ITGB2", "barrier"),
    ("CEACAM1", "CD8A", "barrier"),
    ("LGALS9", "PTPRC", "inhibitory"),
    ("LGALS9", "HAVCR2", "inhibitory"),
    ("HLA-E", "KLRD1", "inhibitory"),
    ("HLA-E", "CD8A", "inhibitory"),
    ("CD274", "PDCD1", "inhibitory"),
    ("HLA-A", "CD8A", "MHC_I"),
    ("HLA-B", "CD8A", "MHC_I"),
    ("HLA-C", "CD8A", "MHC_I"),
    ("CXCL16", "CXCR6", "recruit"),
    ("CXCL9", "CXCR3", "recruit"),
    ("CXCL10", "CXCR3", "recruit"),
    ("CXCL11", "CXCR3", "recruit"),
    ("CCL5", "CCR5", "recruit"),
    ("CX3CL1", "CX3CR1", "recruit"),
]


def log(msg: str) -> None:
    print(msg, flush=True)


def parse_symbols(val) -> list[str]:
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return []
    text = str(val).replace(";", ",").replace("&", ",").replace("|", ",")
    return [x.strip() for x in text.split(",") if x.strip() and x.strip().lower() != "nan"]


def load_cellchat_pairs() -> pd.DataFrame:
    inter = pd.read_csv(ROOT / "db" / "interaction_cellchatdb_v2_protein.csv")
    rows = []
    for rec in inter.itertuples(index=False):
        lig = parse_symbols(getattr(rec, "ligand.symbol", None)) or parse_symbols(rec.ligand)
        recp = parse_symbols(getattr(rec, "receptor.symbol", None)) or parse_symbols(rec.receptor)
        if not lig or not recp:
            continue
        rows.append(
            {
                "source": "cellchat",
                "interaction": rec.interaction_name,
                "pathway": rec.pathway_name,
                "ligand": "+".join(lig),
                "receptor": "+".join(recp),
                "ligand_genes": tuple(lig),
                "receptor_genes": tuple(recp),
            }
        )
    return pd.DataFrame(rows)


def load_liana_pairs() -> pd.DataFrame:
    df = pd.read_csv(ROOT / "db" / "cellphonedb_v5_lr_pairs.tsv", sep="\t")
    rows = []
    for rec in df.itertuples(index=False):
        lig = str(rec.ligand).split("+")
        recp = str(rec.receptor).split("+")
        rows.append(
            {
                "source": "liana",
                "interaction": f"{rec.ligand}_{rec.receptor}",
                "pathway": rec.pathway if pd.notna(rec.pathway) else "",
                "ligand": rec.ligand,
                "receptor": rec.receptor,
                "ligand_genes": tuple(g for g in lig if g),
                "receptor_genes": tuple(g for g in recp if g),
            }
        )
    return pd.DataFrame(rows)


def wanted_from_pairs(pairs: pd.DataFrame) -> set[str]:
    genes = set()
    for rec in pairs.itertuples(index=False):
        genes.update(rec.ligand_genes)
        genes.update(rec.receptor_genes)
    for vs in LINEAGES.values():
        genes.update(vs)
    genes.update(NORMAL_LUNG)
    genes.update(["CLDN4", "TACSTD2", "PTPRC", "CD3D", "CD3E", "CD8A", "NCAM1"])
    return genes


def load_given_n73() -> pd.DataFrame:
    per = pd.read_csv(ROOT / "data" / "given" / "per_donor_metrics.tsv", sep="\t")
    mal = pd.to_numeric(per["mal_CLDN4_mean_log1p"], errors="coerce")
    keep = per["eligible"].astype(str).isin(["True", "true", "1"]) & mal.notna()
    out = per.loc[keep].copy()
    out["donor"] = out["donor"].astype(str)
    if len(out) != 73:
        raise SystemExit(f"expected 73 given patients, got {len(out)}")
    return out


def marker_score(expr: dict[str, np.ndarray], genes: list[str], n: int) -> np.ndarray:
    mats = [np.log1p(np.clip(expr[g], 0, None)) for g in genes if g in expr]
    if not mats:
        return np.zeros(n, dtype=np.float32)
    return np.mean(np.vstack(mats), axis=0).astype(np.float32)


def assign_lineage(expr: dict[str, np.ndarray], n: int) -> np.ndarray:
    names = list(LINEAGES)
    scores = np.vstack([marker_score(expr, LINEAGES[k], n) for k in names])
    best = scores.argmax(axis=0)
    best_val = scores.max(axis=0)
    second = np.partition(scores, -2, axis=0)[-2]
    assigned = np.array(names, dtype=object)[best]
    assigned[(best_val < 0.12) | ((best_val - second) < 0.04)] = "other"
    return assigned


def log1p_cp10k(umi: np.ndarray, total: np.ndarray) -> np.ndarray:
    with np.errstate(divide="ignore", invalid="ignore"):
        cp = np.where(total > 0, umi / total * 1e4, 0.0)
    return np.log1p(cp).astype(np.float32)


def trim_mean(x: np.ndarray, proportiontocut: float = TRIM) -> float:
    x = np.asarray(x, dtype=float)
    n = x.size
    if n == 0:
        return 0.0
    if n == 1:
        return float(x[0])
    k = int(n * proportiontocut)
    if k == 0:
        return float(x.mean())
    s = np.sort(x)
    return float(s[k : n - k].mean())


def geom_mean(vals: list[float]) -> float:
    if any(v <= 0 for v in vals):
        return 0.0
    if len(vals) == 1:
        return float(vals[0])
    return float(np.exp(np.mean(np.log(np.clip(vals, 1e-12, None)))))


def hill_prob(lig: float, rec: float, kh: float = KH) -> float:
    if lig <= 0 or rec <= 0:
        return 0.0
    prod = lig * rec
    return float(prod / (kh + prod))


def stream_gene_rows(path: Path, wanted: set[str], keep_idx: np.ndarray | None = None):
    opener = gzip.open if str(path).endswith(".gz") else open
    found: dict[str, np.ndarray] = {}
    with opener(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        if header[0] in {"", "gene", "Gene", "index", "Index", "GENE"}:
            cell_ids = header[1:]
        else:
            cell_ids = header
        n = len(cell_ids)
        totals = np.zeros(n, dtype=np.float64)
        n_genes = 0
        for line in handle:
            gene, _, rest = line.partition("\t")
            gene = gene.strip().strip('"').split(".")[0]
            arr = np.fromstring(rest, sep="\t", dtype=np.float32)
            if arr.size != n:
                toks = line.rstrip("\n").split("\t")
                arr = np.asarray(toks[1:], dtype=np.float32)
            if arr.size != n:
                continue
            totals += arr
            n_genes += 1
            if gene in wanted:
                found[gene] = arr
            if n_genes % 4000 == 0:
                log(f"  stream {path.name} genes={n_genes} stored={len(found)}")
    log(f"stream done {path.name} cells={n} genes={n_genes} stored={len(found)}")
    if keep_idx is not None:
        cell_ids = [cell_ids[i] for i in np.where(keep_idx)[0]]
        totals = totals[keep_idx]
        found = {g: v[keep_idx] for g, v in found.items()}
    return np.array(cell_ids, dtype=object), found, totals


def load_gse131907(wanted: set[str], donors: set[str]):
    annot = pd.read_csv(DATA / "GSE131907_Lung_Cancer_cell_annotation.txt.gz", sep="\t")
    annot["donor"] = annot["Sample"].astype(str)
    keep_donors = annot["donor"].isin(donors) & annot["Sample_Origin"].isin(GSE131907_TUMOR_ORIGINS)
    log(f"GSE131907 annot rows={len(annot)} keep_donors={int(keep_donors.sum())}")
    cells, expr, totals = stream_gene_rows(
        DATA / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz", wanted
    )
    annot = annot.set_index("Index", drop=False)
    if cells[0] not in annot.index:
        annot["alt"] = annot["Barcode"].astype(str) + "_" + annot["Sample"].astype(str)
        annot = annot.set_index("alt", drop=False)
    order = pd.Index(cells)
    missing = (~order.isin(annot.index)).sum()
    if missing:
        raise SystemExit(f"GSE131907 barcode mismatch missing={missing}")
    annot = annot.loc[order].reset_index(drop=True)
    keep = annot["donor"].isin(donors) & annot["Sample_Origin"].isin(GSE131907_TUMOR_ORIGINS)
    keep_idx = keep.to_numpy()
    annot = annot.loc[keep_idx].copy()
    expr = {g: v[keep_idx] for g, v in expr.items()}
    totals = totals[keep_idx]
    logx = {g: log1p_cp10k(v, totals) for g, v in expr.items()}
    n = len(annot)
    lineage = assign_lineage(expr, n)
    normal = marker_score(expr, [g for g in NORMAL_LUNG if g in expr], n)
    mal = (lineage == "epithelial") & (normal <= 0.05)
    tnk = np.isin(lineage, ["T", "NK"])
    cldn4 = logx["CLDN4"] if "CLDN4" in logx else np.zeros(n, dtype=np.float32)
    return {
        "dataset": "GSE131907",
        "donor": annot["donor"].to_numpy(),
        "logx": logx,
        "mal": mal,
        "tnk": tnk,
        "cldn4": cldn4,
        "scale": "log1p_cp10k",
        "n_cells": n,
    }


def load_gse148071(wanted: set[str], donors: set[str]):
    files = sorted((DATA / "GSE148071_files").glob("*_exp.txt.gz"))
    if not files:
        raise SystemExit("GSE148071 exp files missing")
    expr_all: dict[str, list[np.ndarray]] = {}
    totals_all = []
    donor_all = []
    for fp in files:
        stem = fp.name.replace(".txt.gz", "")
        patient = stem.split("_")[1] if "_" in stem else stem
        if patient not in donors:
            continue
        cells, found, totals = stream_gene_rows(fp, wanted)
        n = len(cells)
        donor_all.append(np.array([patient] * n, dtype=object))
        totals_all.append(totals)
        for g in wanted:
            expr_all.setdefault(g, []).append(found[g] if g in found else np.zeros(n, dtype=np.float32))
        log(f"  GSE148071 {patient} n={n} genes={len(found)}")
    if not donor_all:
        raise SystemExit("GSE148071: no overlapping donors")
    donor = np.concatenate(donor_all)
    totals = np.concatenate(totals_all)
    expr = {g: np.concatenate(vs) for g, vs in expr_all.items() if any(v.sum() for v in vs)}
    logx = {g: log1p_cp10k(v, totals) for g, v in expr.items()}
    n = len(donor)
    lineage = assign_lineage(expr, n)
    normal = marker_score(expr, [g for g in NORMAL_LUNG if g in expr], n)
    mal = (lineage == "epithelial") & (normal <= 0.05)
    tnk = np.isin(lineage, ["T", "NK"])
    cldn4 = logx["CLDN4"] if "CLDN4" in logx else np.zeros(n, dtype=np.float32)
    return {
        "dataset": "GSE148071",
        "donor": donor,
        "logx": logx,
        "mal": mal,
        "tnk": tnk,
        "cldn4": cldn4,
        "scale": "log1p_cp10k",
        "n_cells": n,
    }


def load_gse127465(wanted: set[str], donors: set[str]):
    genes = pd.read_csv(DATA / "GSE127465_gene_names_human_41861.tsv.gz", sep="\t", header=None)[0].astype(str)
    want_idx = {i: g for i, g in enumerate(genes) if g in wanted}
    meta = pd.read_csv(DATA / "GSE127465_human_cell_metadata_54773x25.tsv.gz", sep="\t")
    n_cells = len(meta)
    expr = {g: np.zeros(n_cells, dtype=np.float32) for g in want_idx.values()}
    totals = meta["Total counts"].to_numpy(dtype=np.float64) if "Total counts" in meta.columns else np.zeros(n_cells)
    mtx = DATA / "GSE127465_human_counts_normalized_54773x41861.mtx.gz"
    with gzip.open(mtx, "rt") as handle:
        header = handle.readline()
        while header.startswith("%"):
            header = handle.readline()
        nrow, ncol, nnz = [int(x) for x in header.split()[:3]]
        log(f"GSE127465 MTX {nrow}x{ncol} nnz={nnz}")
        if nrow == n_cells and ncol == len(genes):
            cell_dim, gene_dim = 0, 1
        elif ncol == n_cells and nrow == len(genes):
            cell_dim, gene_dim = 1, 0
        else:
            raise SystemExit(f"MTX shape {nrow}x{ncol} vs cells={n_cells} genes={len(genes)}")
        for i, line in enumerate(handle):
            a, b, val = line.split()[:3]
            r, c = int(a) - 1, int(b) - 1
            gene_i = r if gene_dim == 0 else c
            cell_i = r if cell_dim == 0 else c
            g = want_idx.get(gene_i)
            if g is not None:
                expr[g][cell_i] = float(val)
            if i and i % 5_000_000 == 0:
                log(f"  MTX {i}/{nnz}")
    keep = meta["Tissue"].astype(str).str.lower().eq("tumor") & meta["Patient"].astype(str).isin(donors)
    keep_idx = keep.to_numpy()
    meta = meta.loc[keep_idx].copy()
    expr = {g: v[keep_idx] for g, v in expr.items()}
    totals = totals[keep_idx]
    # Deposited matrix is already library-normalized; do not re-CP10k.
    logx = {g: np.log1p(np.clip(v, 0, None)).astype(np.float32) for g, v in expr.items()}
    n = len(meta)
    lineage = assign_lineage(expr, n)
    normal = marker_score(expr, [g for g in NORMAL_LUNG if g in expr], n)
    mal = (lineage == "epithelial") & (normal <= 0.05)
    tnk = np.isin(lineage, ["T", "NK"])
    cldn4 = logx["CLDN4"] if "CLDN4" in logx else np.zeros(n, dtype=np.float32)
    return {
        "dataset": "GSE127465",
        "donor": meta["Patient"].astype(str).to_numpy(),
        "logx": logx,
        "mal": mal,
        "tnk": tnk,
        "cldn4": cldn4,
        "scale": "log1p_deposited_normalized",
        "n_cells": n,
    }


def group_stats(logx: dict[str, np.ndarray], mask: np.ndarray) -> tuple[dict[str, float], dict[str, float], dict[str, float], int]:
    n = int(mask.sum())
    means, trims, fracs = {}, {}, {}
    if n == 0:
        return means, trims, fracs, 0
    for g, arr in logx.items():
        v = arr[mask]
        means[g] = float(np.mean(v))
        trims[g] = trim_mean(v)
        fracs[g] = float(np.mean(v > 0))
    return means, trims, fracs, n


def score_pairs(pairs: pd.DataFrame, logx: dict[str, np.ndarray], sender: np.ndarray, receiver: np.ndarray) -> pd.DataFrame:
    s_mean, s_trim, s_frac, n_s = group_stats(logx, sender)
    r_mean, r_trim, r_frac, n_r = group_stats(logx, receiver)
    rows = []
    for rec in pairs.itertuples(index=False):
        lig = list(rec.ligand_genes)
        recp = list(rec.receptor_genes)
        if any(g not in s_mean for g in lig) or any(g not in r_mean for g in recp):
            continue
        l_mean = min(s_mean[g] for g in lig)
        rec_mean = min(r_mean[g] for g in recp)
        lig_trim = geom_mean([s_trim[g] for g in lig])
        rec_trim = geom_mean([r_trim[g] for g in recp])
        l_frac = min(s_frac[g] for g in lig)
        rec_frac = min(r_frac[g] for g in recp)
        rows.append(
            {
                "source": rec.source,
                "interaction": rec.interaction,
                "pathway": rec.pathway,
                "ligand": rec.ligand,
                "receptor": rec.receptor,
                "n_sender": n_s,
                "n_receiver": n_r,
                "ligand_mean": l_mean,
                "receptor_mean": rec_mean,
                "ligand_frac": l_frac,
                "receptor_frac": rec_frac,
                "liana_score": 0.5 * (l_mean + rec_mean),
                "cellchat_P": hill_prob(lig_trim, rec_trim),
                "pass_expr_prop": (l_frac >= EXPR_PROP) and (rec_frac >= EXPR_PROP),
            }
        )
    return pd.DataFrame(rows)


def patient_deltas(obj: dict, pairs: pd.DataFrame, given_donors: set[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    donors = obj["donor"]
    mal = obj["mal"]
    tnk = obj["tnk"]
    cldn4 = obj["cldn4"]
    logx = obj["logx"]
    ntab_rows = []
    delta_rows = []
    for donor in sorted(set(donors.astype(str))):
        if donor not in given_donors:
            continue
        dmask = donors == donor
        d_mal = dmask & mal
        d_tnk = dmask & tnk
        n_mal = int(d_mal.sum())
        n_tnk = int(d_tnk.sum())
        if n_mal == 0:
            thr = np.nan
            high = d_mal
            low = d_mal
        else:
            thr = float(np.median(cldn4[d_mal]))
            high = d_mal & (cldn4 >= thr)
            low = d_mal & (cldn4 < thr)
        n_hi, n_lo = int(high.sum()), int(low.sum())
        paired = n_hi >= MIN_MAL_ARM and n_lo >= MIN_MAL_ARM and n_tnk >= MIN_TNK
        ntab_rows.append(
            {
                "dataset": obj["dataset"],
                "donor": donor,
                "n_cells": int(dmask.sum()),
                "n_malignant_like": n_mal,
                "n_cldn4_high": n_hi,
                "n_cldn4_low": n_lo,
                "n_tnk": n_tnk,
                "cldn4_threshold": thr,
                "mean_cldn4_high": float(cldn4[high].mean()) if n_hi else np.nan,
                "mean_cldn4_low": float(cldn4[low].mean()) if n_lo else np.nan,
                "paired": paired,
                "scale": obj["scale"],
            }
        )
        if not paired:
            continue
        hi = score_pairs(pairs, logx, high, d_tnk)
        lo = score_pairs(pairs, logx, low, d_tnk)
        if hi.empty or lo.empty:
            continue
        key = ["source", "interaction", "pathway", "ligand", "receptor"]
        m = hi.merge(lo, on=key, suffixes=("_high", "_low"))
        m["dataset"] = obj["dataset"]
        m["donor"] = donor
        m["delta_liana"] = m["liana_score_high"] - m["liana_score_low"]
        m["delta_cellchat"] = m["cellchat_P_high"] - m["cellchat_P_low"]
        m["pass_either"] = m["pass_expr_prop_high"] | m["pass_expr_prop_low"]
        delta_rows.append(m)
    ntab = pd.DataFrame(ntab_rows)
    deltas = pd.concat(delta_rows, ignore_index=True) if delta_rows else pd.DataFrame()
    return ntab, deltas


def wilcoxon_p(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) < 4 or np.allclose(x, 0):
        return np.nan
    try:
        return float(stats.wilcoxon(x, zero_method="wilcox", alternative="two-sided").pvalue)
    except ValueError:
        return np.nan


def dl_random_effects(means: np.ndarray, ses: np.ndarray) -> dict:
    means = np.asarray(means, dtype=float)
    ses = np.asarray(ses, dtype=float)
    ok = np.isfinite(means) & np.isfinite(ses) & (ses > 0)
    means, ses = means[ok], ses[ok]
    if len(means) < 2:
        return {"k": int(len(means)), "mu": np.nan, "se": np.nan, "p": np.nan, "I2": np.nan}
    w = 1.0 / (ses**2)
    mu_fe = np.sum(w * means) / np.sum(w)
    q = np.sum(w * (means - mu_fe) ** 2)
    df = len(means) - 1
    c = np.sum(w) - np.sum(w**2) / np.sum(w)
    tau2 = max(0.0, (q - df) / c) if c > 0 else 0.0
    w_re = 1.0 / (ses**2 + tau2)
    mu = np.sum(w_re * means) / np.sum(w_re)
    se = np.sqrt(1.0 / np.sum(w_re))
    z = mu / se if se > 0 else np.nan
    p = float(2 * stats.norm.sf(abs(z))) if np.isfinite(z) else np.nan
    i2 = max(0.0, (q - df) / q) * 100 if q > 0 else 0.0
    return {"k": int(len(means)), "mu": float(mu), "se": float(se), "p": p, "I2": float(i2)}


def rank_pairs(deltas: pd.DataFrame, score_col: str) -> pd.DataFrame:
    rows = []
    for (src, inter, path, lig, rec), g in deltas.groupby(
        ["source", "interaction", "pathway", "ligand", "receptor"], sort=False
    ):
        g = g[g["pass_either"]].copy()
        if g.empty:
            continue
        x = g[score_col].to_numpy(dtype=float)
        x = x[np.isfinite(x)]
        if len(x) < 4:
            continue
        med = float(np.median(x))
        p = wilcoxon_p(x)
        n_pos = int((x > 0).sum())
        n_neg = int((x < 0).sum())
        rec_row = {
            "source": src,
            "interaction": inter,
            "pathway": path,
            "ligand": lig,
            "receptor": rec,
            "n_patients": int(len(x)),
            "n_cohorts": int(g["dataset"].nunique()),
            "median_delta": med,
            "mean_delta": float(np.mean(x)),
            "n_pos": n_pos,
            "n_neg": n_neg,
            "wilcoxon_p": p,
        }
        # cohort RE on mean delta
        cmeans, cses, cnames = [], [], []
        for ds, sub in g.groupby("dataset"):
            xx = sub[score_col].to_numpy(dtype=float)
            xx = xx[np.isfinite(xx)]
            if len(xx) < 4:
                continue
            cmeans.append(float(np.mean(xx)))
            cses.append(float(np.std(xx, ddof=1) / np.sqrt(len(xx))) if len(xx) > 1 else np.nan)
            cnames.append(ds)
        re = dl_random_effects(np.array(cmeans), np.array(cses))
        rec_row.update(
            {
                "re_k": re["k"],
                "re_mu": re["mu"],
                "re_p": re["p"],
                "re_I2": re["I2"],
                "cohorts": "+".join(cnames),
            }
        )
        rows.append(rec_row)
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    mask = out["wilcoxon_p"].notna()
    if mask.any():
        out.loc[mask, "fdr"] = multipletests(out.loc[mask, "wilcoxon_p"], method="fdr_bh")[1]
    else:
        out["fdr"] = np.nan
    return out.sort_values(["fdr", "wilcoxon_p", "median_delta"], na_position="last")


def make_figures(ntab: pd.DataFrame, ranks_liana: pd.DataFrame, ranks_cc: pd.DataFrame, focus_pat: pd.DataFrame) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    paired = ntab[ntab["paired"]].copy()
    if paired.empty:
        return

    fig, ax = plt.subplots(figsize=(max(8, 0.22 * len(paired)), 4.2))
    x = np.arange(len(paired))
    ax.bar(x - 0.2, paired["n_cldn4_high"], width=0.2, label="CLDN4-high mal", color="#4C72B0")
    ax.bar(x, paired["n_cldn4_low"], width=0.2, label="CLDN4-low mal", color="#55A868")
    ax.bar(x + 0.2, paired["n_tnk"], width=0.2, label="T/NK", color="#DD8452")
    ax.set_xticks(x)
    labels = [f"{r.dataset.replace('GSE','')}\n{r.donor}" for r in paired.itertuples(index=False)]
    ax.set_xticklabels(labels, fontsize=6)
    ax.set_ylabel("cells")
    ax.set_title(f"Paired LR patients (n={len(paired)}; min {MIN_MAL_ARM}/{MIN_MAL_ARM} mal + {MIN_TNK} T/NK)")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "fig_n_patients.png", dpi=160)
    fig.savefig(FIG / "fig_n_patients.pdf")
    plt.close(fig)

    if not focus_pat.empty:
        order = [f"{a}–{b}" for a, b, _ in FOCUS]
        focus_pat = focus_pat.copy()
        focus_pat["pair"] = focus_pat["ligand"] + "–" + focus_pat["receptor"]
        present = [p for p in order if p in set(focus_pat["pair"])]
        fig, ax = plt.subplots(figsize=(8.5, max(3.5, 0.38 * len(present))))
        for i, p in enumerate(present):
            sub = focus_pat[focus_pat["pair"] == p]
            y = sub["delta_liana"].to_numpy(dtype=float)
            y = y[np.isfinite(y)]
            jitter = np.random.default_rng(SEED).uniform(-0.12, 0.12, size=len(y))
            ax.scatter(y, np.full(len(y), i) + jitter, s=16, alpha=0.7, c="#4C72B0", linewidths=0)
            if len(y):
                ax.plot([np.median(y)], [i], "D", color="#C44E52", ms=6)
        ax.axvline(0, color="0.4", lw=0.8)
        ax.set_yticks(range(len(present)))
        ax.set_yticklabels(present, fontsize=8)
        ax.set_xlabel("patient Δ LIANA-style score (CLDN4-high − low)")
        ax.set_title("Outgoing CLDN4-high malignant → T/NK (patient unit)")
        fig.tight_layout()
        fig.savefig(FIG / "fig_focus_patient_deltas.png", dpi=160)
        fig.savefig(FIG / "fig_focus_patient_deltas.pdf")
        plt.close(fig)

    if not ranks_liana.empty:
        show = ranks_liana[ranks_liana["n_patients"] >= 8].head(20)
        if not show.empty:
            fig, ax = plt.subplots(figsize=(8.2, 5.6))
            y = np.arange(len(show))
            ax.barh(y, show["median_delta"], color=np.where(show["median_delta"] >= 0, "#4C72B0", "#C44E52"))
            ax.set_yticks(y)
            ax.set_yticklabels([f"{r.ligand}–{r.receptor}" for r in show.itertuples(index=False)], fontsize=7)
            ax.axvline(0, color="0.3", lw=0.8)
            ax.invert_yaxis()
            ax.set_xlabel("median patient Δ (LIANA-style)")
            ax.set_title("Top outgoing pairs by FDR (n_patients ≥ 8)")
            fig.tight_layout()
            fig.savefig(FIG / "fig_top_outgoing_liana.png", dpi=160)
            fig.savefig(FIG / "fig_top_outgoing_liana.pdf")
            plt.close(fig)

    # ligand-table image for focus ranks
    if not ranks_liana.empty:
        foc_keys = {(a, b) for a, b, _ in FOCUS}
        foc = ranks_liana[ranks_liana.apply(lambda r: (r.ligand, r.receptor) in foc_keys, axis=1)].copy()
        if not foc.empty:
            fig, ax = plt.subplots(figsize=(10.5, 0.42 * (len(foc) + 3)))
            ax.axis("off")
            cols = ["ligand", "receptor", "n_patients", "median_delta", "wilcoxon_p", "fdr", "n_pos", "n_neg"]
            cell = foc[cols].copy()
            cell["median_delta"] = cell["median_delta"].map(lambda v: f"{v:+.3f}")
            cell["wilcoxon_p"] = cell["wilcoxon_p"].map(lambda v: "NA" if not np.isfinite(v) else f"{v:.3g}")
            cell["fdr"] = cell["fdr"].map(lambda v: "NA" if not np.isfinite(v) else f"{v:.3g}")
            tbl = ax.table(
                cellText=cell.values,
                colLabels=cols,
                loc="center",
                cellLoc="center",
            )
            tbl.auto_set_font_size(False)
            tbl.set_fontsize(7)
            tbl.scale(1, 1.25)
            ax.set_title("Focus outgoing pairs — LIANA-style patient Δ (honest n)", pad=12)
            fig.tight_layout()
            fig.savefig(FIG / "fig_ligand_table.png", dpi=160)
            fig.savefig(FIG / "fig_ligand_table.pdf")
            plt.close(fig)


def write_finding(summary: dict, ntab: pd.DataFrame, ranks_liana: pd.DataFrame, ranks_cc: pd.DataFrame) -> None:
    paired = ntab[ntab["paired"]].copy()
    dropped = ntab[~ntab["paired"]].copy()
    given = summary["n_given"]
    n_lr = int(paired.shape[0])
    n_skip253 = summary["n_gse253013"]

    def fmt(v, nd=3):
        if v is None or not np.isfinite(v):
            return "NA"
        return f"{v:.{nd}g}" if nd == 3 and abs(v) < 0.01 else f"{v:.3f}" if isinstance(v, float) and abs(v) >= 0.01 else f"{v:.3g}"

    foc_keys = {(a, b): c for a, b, c in FOCUS}
    foc = ranks_liana[ranks_liana.apply(lambda r: (r.ligand, r.receptor) in foc_keys, axis=1)].copy()
    if not foc.empty:
        foc["pair_class"] = [foc_keys[(r.ligand, r.receptor)] for r in foc.itertuples(index=False)]
        foc = foc.sort_values(["pair_class", "median_delta"])

    n_neg = int((foc["median_delta"] < 0).sum()) if not foc.empty else 0
    n_sig = int((foc["fdr"] < 0.05).sum()) if (not foc.empty and foc["fdr"].notna().any()) else 0

    # one-sentence verdict from focus recruit vs barrier
    rec_sub = foc[foc["pair_class"] == "recruit"] if not foc.empty else foc
    bar_sub = foc[foc["pair_class"].isin(["barrier", "inhibitory"])] if not foc.empty else foc
    rec_med = float(rec_sub["median_delta"].median()) if len(rec_sub) else np.nan
    bar_med = float(bar_sub["median_delta"].median()) if len(bar_sub) else np.nan
    verdict = (
        f"On the Harmony n=73 table (given ρ=−0.27, not re-audited), outgoing LR was scored "
        f"in **{n_lr}** paired patients from GSE131907+GSE148071+GSE127465 "
        f"(GSE253013 n={n_skip253} skipped: 9 GB RDS not downloaded). "
        f"Focus pairs: {n_neg}/{len(foc) if not foc.empty else 0} have median Δ<0; "
        f"**{n_sig}/{len(foc) if not foc.empty else 0} reach FDR<0.05**. "
        f"Recruit median Δ={fmt(rec_med)}; barrier/inhibitory median Δ={fmt(bar_med)}."
    )

    lines = [
        "# FINDING — Harmony n=73 CLDN4-only high-end (CellChat-style + LIANA-style)",
        "",
        "**ADDITIVE. CLDN4 only. No dual-high.**",
        "",
        "The Harmony / multi-cohort patient table is taken as given and was **not** re-audited:",
        "malignant-like *CLDN4* vs T/NK **n=73, ρ=−0.27** (GSE131907 + GSE253013 + GSE148071 + GSE127465;",
        "`data/given/per_donor_metrics.tsv`, `association_statistics.tsv`).",
        "",
        verdict,
        "",
        "## What was run",
        "",
        "| Method | Status |",
        "|---|---|",
        "| Given Harmony malignant CLDN4 vs T/NK (n=73, ρ=−0.27) | **taken as given — not re-audited** |",
        "| CellChat-style Hill P (10% truncated means, K_h=0.5, CellChatDB v2 protein) | **run** — patient-level Δ, then meta |",
        "| LIANA / CellPhoneDB-style score (mean of min-subunit means on log1p) | **primary** — patient-level Wilcoxon + cohort DL |",
        "| LIANA R/Python `liana.mt.cellphonedb` | **not run** — no joint object; documented score used instead |",
        "| CellChat R | **not run** — R package not used |",
        "| Milo-style nhood on the Harmony joint embedding | **not run** — embedding is not stored; GSE253013 RDS not downloaded so it cannot be rebuilt |",
        "",
        "## Honest n",
        "",
        "| Set | n |",
        "|---|---:|",
        f"| Given Harmony patients (malignant CLDN4 vs T/NK) | **73** |",
        f"| GSE253013 in that table (no public processed UMI used here) | {n_skip253} |",
        f"| Patients with a public processed matrix (131907+148071+127465) | {summary['n_with_matrix']} |",
        f"| Patients scored for LR (seen in a matrix) | {int(len(ntab))} |",
        f"| **Paired LR patients** (≥{MIN_MAL_ARM} CLDN4-high **and** ≥{MIN_MAL_ARM} CLDN4-low malignant-like, ≥{MIN_TNK} T/NK) | **{n_lr}** |",
        f"| Paired GSE131907 / GSE148071 / GSE127465 | {summary['n_paired_131907']} / {summary['n_paired_148071']} / {summary['n_paired_127465']} |",
        f"| Dropped (one-sided / thin malignant or T/NK) | {int(len(dropped))} |",
        f"| CellChatDB v2 protein pairs scored | {summary['n_cellchat_pairs']} |",
        f"| CellPhoneDB v5 pairs scored | {summary['n_liana_pairs']} |",
        "",
        "The header n for the Wilcoxon is the paired n, **not 73**. GSE253013's 9 patients stay in the given ρ=−0.27 row and are absent from the LR table.",
        "",
        "CLDN4-high / low is the **patient-specific median** of malignant-like *CLDN4* (same Harmony malignant-like rule: epithelial marker-argmax and near-zero normal-lung score). TACSTD2 is not a gate.",
        "",
        "### Paired patients",
        "",
    ]
    if paired.empty:
        lines.append("None.")
    else:
        lines += [
            "| Dataset | Donor | mal high | mal low | T/NK |",
            "|---|---|---:|---:|---:|",
        ]
        for r in paired.itertuples(index=False):
            lines.append(
                f"| {r.dataset} | {r.donor} | {int(r.n_cldn4_high)} | {int(r.n_cldn4_low)} | {int(r.n_tnk)} |"
            )
    lines += [
        "",
        "Per-patient counts: `results/n_patients.tsv`.",
        "",
        "## Score",
        "",
        "- **LIANA-style (primary):** on log1p(CP10k) for GSE131907 / GSE148071, and log1p of the deposited normalized counts for GSE127465. Partner expression = minimum subunit mean. Pair score = mean of the two partner means (Efremova 2020; Garcia-Alonso 2022). `pass_expr_prop` = both partners in ≥10% of cells.",
        "- **CellChat-style (companion):** 10% truncated means, geometric mean of subunits, Hill P = LR / (0.5 + LR). Not a CellChat R communication probability.",
        "- Unit = **patient**. Receiver T/NK are that patient's own T/NK (so outgoing Δ is ligand-driven). Cells are not replicates.",
        "- Meta: Wilcoxon signed-rank on patient Δ (high−low); BH-FDR within method. Cohort DerSimonian–Laird on cohort mean Δ when ≥2 cohorts have ≥4 patients.",
        "",
        "## Primary ligand table — outgoing CLDN4-high malignant → T/NK (LIANA-style)",
        "",
        "Median patient Δ = high − low. Negative = weaker from the CLDN4-high state.",
        "",
    ]
    if foc.empty:
        lines.append("No focus pairs reached n_patients ≥ 4.")
    else:
        lines += [
            "| Class | Pair | n patients | n cohorts | median Δ | Wilcoxon p | FDR | n+ / n− |",
            "|---|---|---:|---:|---:|---:|---:|---:|",
        ]
        for r in foc.itertuples(index=False):
            lines.append(
                f"| {r.pair_class} | {r.ligand}–{r.receptor} | {int(r.n_patients)} | {int(r.n_cohorts)} | "
                f"{r.median_delta:+.3f} | {fmt(r.wilcoxon_p)} | {fmt(r.fdr)} | {int(r.n_pos)}/{int(r.n_neg)} |"
            )
        lines += [
            "",
            f"{n_neg}/{len(foc)} focus pairs have median Δ < 0. **{n_sig}/{len(foc)} reach FDR < 0.05**.",
            "",
        ]
    lines += [
        "Full ranked table: `results/ligand_table_liana_outgoing.tsv` (primary) and `results/ligand_table_cellchat_outgoing.tsv`.",
        "Patient-level deltas: `results/patient_deltas.tsv`.",
        "",
        "## CellChat-style companion (same patients)",
        "",
    ]
    foc_cc = ranks_cc[ranks_cc.apply(lambda r: (r.ligand, r.receptor) in foc_keys, axis=1)].copy() if not ranks_cc.empty else ranks_cc
    if foc_cc.empty:
        lines.append("No focus CellChat-style pairs with n_patients ≥ 4.")
    else:
        foc_cc["pair_class"] = [foc_keys[(r.ligand, r.receptor)] for r in foc_cc.itertuples(index=False)]
        foc_cc = foc_cc.sort_values(["pair_class", "median_delta"])
        lines += [
            "| Class | Pair | n | median ΔP | Wilcoxon p | FDR |",
            "|---|---|---:|---:|---:|---:|",
        ]
        for r in foc_cc.itertuples(index=False):
            lines.append(
                f"| {r.pair_class} | {r.ligand}–{r.receptor} | {int(r.n_patients)} | "
                f"{r.median_delta:+.3f} | {fmt(r.wilcoxon_p)} | {fmt(r.fdr)} |"
            )
        lines += [""]
    lines += [
        "## Milo / joint embedding",
        "",
        "The Harmony joint UMAP/PCA is **not** in the repository (only donor scores and PNG figures).",
        "Rebuilding it requires the GSE253013 Garnett RDS (~9 GB), which was not downloaded.",
        "No Milo-style neighborhood table is claimed. File: `results/nhood_skip.json`.",
        "",
        "## What is not claimed",
        "",
        "- The given n=73 ρ=−0.27 was not recomputed.",
        "- Dual-high (TACSTD2 **and** CLDN4) was not run.",
        "- Cell-level p-values are not reported.",
        "- GSE253013 LR is absent (honest skip).",
        "- A rebuilt Harmony embedding / Milo DA test.",
        "- ICI response (these four series are treatment-naïve / diagnostic atlases).",
        "",
        "## Files",
        "",
        "| File | Role |",
        "|---|---|",
        "| `FINDING.md` | This note |",
        "| `results/ligand_table_liana_outgoing.tsv` | Primary ligand table (honest n) |",
        "| `results/ligand_table_cellchat_outgoing.tsv` | CellChat-style companion |",
        "| `results/ligand_table_focus.tsv` | Curated outgoing pairs |",
        "| `results/n_patients.tsv` | Per-patient cell counts |",
        "| `results/patient_deltas.tsv` | Per-patient pair Δ |",
        "| `results/nhood_skip.json` | Why Milo was not run |",
        "| `results/summary.json` | Machine-readable n |",
        "| `results/figures/` | Extra figures |",
        "",
        "## Reproduce",
        "",
        "```bash",
        "cd methods/harmony73_cldn4_hiend",
        "python3 scripts/00_download.py   # does not fetch GSE253013",
        "python3 scripts/01_run_lr.py",
        "```",
        "",
    ]
    text = "\n".join(lines) + "\n"
    (ROOT / "FINDING.md").write_text(text)
    (OUT / "FINDING.md").write_text(text)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    given = load_given_n73()
    given.to_csv(OUT / "given_n73_patients.tsv", sep="\t", index=False)
    donors_by = {ds: set(sub["donor"].astype(str)) for ds, sub in given.groupby("dataset")}
    log(f"given n=73 by cohort: { {k: len(v) for k, v in donors_by.items()} }")

    cc = load_cellchat_pairs()
    li = load_liana_pairs()
    # Score the union but keep source tags. Dedup within source on ligand+receptor.
    pairs = pd.concat([cc, li], ignore_index=True)
    pairs = pairs.drop_duplicates(["source", "ligand", "receptor"])
    wanted = wanted_from_pairs(pairs)
    log(f"pairs cellchat={len(cc)} liana={len(li)} union_scored={len(pairs)} genes={len(wanted)}")

    ntabs = []
    deltas = []
    loaders = {
        "GSE131907": load_gse131907,
        "GSE148071": load_gse148071,
        "GSE127465": load_gse127465,
    }
    for ds, fn in loaders.items():
        dset = donors_by.get(ds, set())
        if not dset:
            continue
        log(f"==== {ds} donors={len(dset)} ====")
        obj = fn(wanted, dset)
        nt, de = patient_deltas(obj, pairs, dset)
        ntabs.append(nt)
        if not de.empty:
            deltas.append(de)
        log(f"{ds}: scored_patients={int(nt['paired'].sum()) if len(nt) else 0} cells={obj['n_cells']}")
        del obj

    ntab = pd.concat(ntabs, ignore_index=True) if ntabs else pd.DataFrame()
    all_delta = pd.concat(deltas, ignore_index=True) if deltas else pd.DataFrame()
    ntab.to_csv(OUT / "n_patients.tsv", sep="\t", index=False)
    if not all_delta.empty:
        all_delta.to_csv(OUT / "patient_deltas.tsv", sep="\t", index=False)

    ranks_liana = rank_pairs(all_delta[all_delta["source"] == "liana"], "delta_liana") if not all_delta.empty else pd.DataFrame()
    ranks_cc = rank_pairs(all_delta[all_delta["source"] == "cellchat"], "delta_cellchat") if not all_delta.empty else pd.DataFrame()
    if not ranks_liana.empty:
        ranks_liana.to_csv(OUT / "ligand_table_liana_outgoing.tsv", sep="\t", index=False)
    if not ranks_cc.empty:
        ranks_cc.to_csv(OUT / "ligand_table_cellchat_outgoing.tsv", sep="\t", index=False)

    foc_keys = {(a, b): c for a, b, c in FOCUS}
    if not ranks_liana.empty:
        foc = ranks_liana[ranks_liana.apply(lambda r: (r.ligand, r.receptor) in foc_keys, axis=1)].copy()
        if not foc.empty:
            foc["pair_class"] = [foc_keys[(r.ligand, r.receptor)] for r in foc.itertuples(index=False)]
            foc.to_csv(OUT / "ligand_table_focus.tsv", sep="\t", index=False)
    focus_pat = pd.DataFrame()
    if not all_delta.empty:
        focus_pat = all_delta[
            (all_delta["source"] == "liana")
            & all_delta.apply(lambda r: (r.ligand, r.receptor) in foc_keys, axis=1)
            & all_delta["pass_either"]
        ].copy()

    make_figures(ntab, ranks_liana, ranks_cc, focus_pat)

    paired = ntab[ntab["paired"]] if len(ntab) else ntab
    summary = {
        "n_given": 73,
        "given_rho_cldn4_tnk": -0.27,
        "given_note": "not_re_audited",
        "n_gse253013": int((given["dataset"] == "GSE253013").sum()),
        "n_with_matrix": int((given["dataset"] != "GSE253013").sum()),
        "n_scored_in_matrix": int(len(ntab)),
        "n_paired": int(len(paired)),
        "n_paired_131907": int((paired["dataset"] == "GSE131907").sum()) if len(paired) else 0,
        "n_paired_148071": int((paired["dataset"] == "GSE148071").sum()) if len(paired) else 0,
        "n_paired_127465": int((paired["dataset"] == "GSE127465").sum()) if len(paired) else 0,
        "n_cellchat_pairs": int(len(cc)),
        "n_liana_pairs": int(len(li)),
        "min_mal_arm": MIN_MAL_ARM,
        "min_tnk": MIN_TNK,
        "expr_prop": EXPR_PROP,
        "gse253013": "skipped_9gb_rds_not_downloaded",
        "milo": "skipped_no_stored_joint_embedding",
        "dual_high": False,
        "unit": "patient",
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (OUT / "nhood_skip.json").write_text(
        json.dumps(
            {
                "ran": False,
                "reason": (
                    "Harmony joint embedding (harmonypy PCA/UMAP) is not stored in the repo. "
                    "Rebuilding it requires GSE253013_all_luad_garnett_temp.rds.gz (~9 GB), "
                    "which this analysis does not download. No Milo-style nhood table."
                ),
                "n_given": 73,
            },
            indent=2,
        )
        + "\n"
    )
    write_finding(summary, ntab, ranks_liana, ranks_cc)
    log(f"DONE paired={summary['n_paired']} ligand_rows={0 if ranks_liana.empty else len(ranks_liana)}")


if __name__ == "__main__":
    main()
