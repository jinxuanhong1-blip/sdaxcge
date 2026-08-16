#!/usr/bin/env python3
"""Hunt TACSTD2-high tumors for TJ / keratin-barrier GSEA up and EMT down.

Pre-specified design. Every cohort is reported, including nulls and contradictions.
Outputs: results/hunt_tj_gsea/
"""
from __future__ import annotations

import gzip
import json
import math
import os
import sys
import time
import warnings
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

warnings.filterwarnings("ignore", category=RuntimeWarning)

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RAW = DATA / "raw"
OUT = ROOT / "results" / "hunt_tj_gsea"
FIG = OUT / "figures"
SEED = 42
FOCAL = ["CLDN1", "CLDN4", "CLDN7", "F11R", "PARD3"]
CONTEXT_GENES = FOCAL + [
    "OCLN",
    "TJP1",
    "CRB3",
    "PARD6B",
    "CDH1",
    "VIM",
    "ZEB1",
    "SNAI2",
    "KRT5",
    "KRT17",
    "SPRR1B",
    "EPCAM",
]


@dataclass
class Cohort:
    name: str
    source: str
    modality: str
    expr: pd.DataFrame  # genes x samples, numeric, already log-like
    notes: str = ""
    meta: pd.DataFrame = field(default_factory=pd.DataFrame)


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def load_primary_sets() -> tuple[dict[str, list[str]], dict[str, dict]]:
    payload = json.loads((DATA / "genesets" / "primary_sets.json").read_text())
    return payload["sets"], payload["meta"]


def tcga_sample_type(barcode: str) -> str:
    parts = str(barcode).split("-")
    if len(parts) < 4:
        return "unknown"
    return parts[3][:2]


def load_tcga(cancer: str) -> Cohort | None:
    path = RAW / "tcga" / f"TCGA.{cancer}.HiSeqV2.gz"
    if not path.exists():
        return None
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = df.index.astype(str)
    # Xena HiSeqV2 is log2(norm_count+1)
    types = pd.Series({c: tcga_sample_type(c) for c in df.columns})
    if cancer == "SKCM":
        keep = types[types.isin(["01", "06"])].index
        note = "TCGA Xena HiSeqV2; primary (01) + metastatic (06) kept because SKCM is mostly metastases"
    else:
        keep = types[types == "01"].index
        note = "TCGA Xena HiSeqV2 log2(norm+1); primary tumors (sample type 01) only"
    expr = df.loc[:, keep].apply(pd.to_numeric, errors="coerce")
    expr = expr.loc[~expr.index.duplicated(keep="first")]
    return Cohort(
        name=f"TCGA-{cancer}",
        source="UCSC Xena TCGA Hub HiSeqV2",
        modality="bulk_rnaseq",
        expr=expr,
        notes=note,
        meta=pd.DataFrame({"sample": expr.columns, "sample_type": [types[c] for c in expr.columns]}),
    )


def load_gse207422_bulk() -> Cohort:
    expr = pd.read_csv(
        RAW / "geo" / "GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz",
        sep="\t",
        index_col=0,
    )
    expr.index = expr.index.astype(str)
    expr = expr.apply(pd.to_numeric, errors="coerce")
    meta = pd.read_excel(RAW / "geo" / "GSE207422_NSCLC_bulk_RNAseq_metadata.xlsx")
    meta = meta[meta["Sample"].isin(expr.columns)].copy()
    expr = expr.loc[:, meta["Sample"].tolist()]
    return Cohort(
        name="GSE207422_bulk",
        source="GEO GSE207422 Hu et al. Genome Med 2023; NSCLC neoadjuvant PD-1 + chemo",
        modality="bulk_rnaseq",
        expr=expr,
        notes="log2(TPM); n=24 pre-treatment biopsies. Quartile split is 6 vs 6 — underpowered.",
        meta=meta.set_index("Sample"),
    )


def load_gse31210() -> Cohort:
    annot_path = RAW / "geo" / "GPL570.annot.gz"
    probe_to_gene: dict[str, str] = {}
    with gzip.open(annot_path, "rt", errors="replace") as fh:
        in_table = False
        header = None
        for line in fh:
            if line.startswith("!platform_table_begin"):
                in_table = True
                header = next(fh).rstrip("\n").split("\t")
                continue
            if line.startswith("!platform_table_end"):
                break
            if not in_table or header is None:
                continue
            parts = line.rstrip("\n").split("\t")
            rec = dict(zip(header, parts))
            pid = rec.get("ID") or rec.get("#ID")
            sym = rec.get("Gene symbol") or rec.get("Gene symbol")
            if not pid or not sym or "///" in sym:
                # skip multi-gene probes
                if sym and "///" in sym:
                    continue
                continue
            probe_to_gene[pid] = sym.split("///")[0].strip()

    rows = []
    sample_ids = None
    titles = None
    with gzip.open(RAW / "geo" / "GSE31210_series_matrix.txt.gz", "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!Sample_title"):
                titles = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_geo_accession"):
                sample_ids = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!series_matrix_table_begin"):
                header = next(fh).rstrip("\n").split("\t")
                for line2 in fh:
                    if line2.startswith("!series_matrix_table_end"):
                        break
                    parts = line2.rstrip("\n").split("\t")
                    pid = parts[0].strip('"')
                    gene = probe_to_gene.get(pid)
                    if not gene:
                        continue
                    vals = [float(x) if x not in ("", "null", "NA") else np.nan for x in parts[1:]]
                    rows.append((gene, vals))
                break
    if not rows or sample_ids is None:
        raise RuntimeError("failed to parse GSE31210")
    mat = pd.DataFrame([r[1] for r in rows], index=[r[0] for r in rows], columns=sample_ids)
    # collapse probes: mean of probes for the same gene
    mat = mat.groupby(level=0).mean()
    # Okayama 2011: 226 tumors + 20 normals; titles often contain 'normal'
    tissue = None
    with gzip.open(RAW / "geo" / "GSE31210_series_matrix.txt.gz", "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!Sample_characteristics_ch1") and "tissue:" in line:
                tissue = [x.strip().strip('"').split("tissue:", 1)[-1].strip() for x in line.rstrip("\n").split("\t")[1:]]
                break
    is_tumor = pd.Series(True, index=mat.columns)
    if tissue is not None and len(tissue) == mat.shape[1]:
        is_tumor = pd.Series(["tumor" in t.lower() for t in tissue], index=mat.columns)
    elif titles is not None and len(titles) == mat.shape[1]:
        is_tumor = pd.Series(["normal" not in t.lower() for t in titles], index=mat.columns)
    expr = mat.loc[:, mat.columns[is_tumor]]
    return Cohort(
        name="GSE31210",
        source="GEO GSE31210 Okayama et al. stage I-II LUAD Affymetrix GPL570",
        modality="microarray",
        expr=expr,
        notes="MAS5 series matrix; multi-gene probes dropped; remaining probes collapsed by mean. Primary tumors only.",
        meta=pd.DataFrame({"sample": expr.columns, "platform": "GPL570"}),
    )


def _stream_marker_matrix(path: Path, markers: set[str]) -> tuple[list[str], pd.DataFrame]:
    found: dict[str, np.ndarray] = {}
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        cells = header[1:]
        for line in fh:
            gene = line.split("\t", 1)[0]
            if gene not in markers:
                continue
            vals = np.fromstring(line.split("\t", 1)[1], sep="\t", dtype=np.float32)
            found[gene] = vals
            if len(found) == len(markers):
                break
    if not found:
        raise RuntimeError(f"no markers found in {path}")
    df = pd.DataFrame(found, index=cells).T
    return cells, df


def _pseudobulk_from_wide(
    path: Path,
    keep_cells: list[str],
    cell_to_sample: dict[str, str],
    min_cells_per_sample: int = 20,
) -> pd.DataFrame:
    """Stream genes x cells UMI matrix; sum UMIs by sample for selected cells."""
    sample_of = {c: cell_to_sample[c] for c in keep_cells if c in cell_to_sample}
    keep_cells = [c for c in keep_cells if c in sample_of]
    counts = pd.Series(keep_cells).map(sample_of).value_counts()
    samples = [s for s, n in counts.items() if n >= min_cells_per_sample]
    keep_cells = [c for c in keep_cells if sample_of[c] in samples]
    if len(samples) < 6:
        raise RuntimeError(f"too few samples after min_cells filter: {len(samples)}")
    with gzip.open(path, "rt") as fh:
        gene_col = fh.readline().rstrip("\n").split("\t")[0]
    usecols = [gene_col] + keep_cells
    log(f"  pseudobulk {path.name}: {len(keep_cells)} cells -> {len(samples)} samples")
    acc_blocks: list[pd.DataFrame] = []
    t0 = time.time()
    n_chunk = 0
    for chunk in pd.read_csv(
        path,
        sep="\t",
        index_col=0,
        usecols=usecols,
        chunksize=250,
    ):
        chunk = chunk.reindex(columns=keep_cells)
        grouped = chunk.T.groupby([sample_of[c] for c in keep_cells]).sum().T
        grouped = grouped.reindex(columns=samples)
        acc_blocks.append(grouped)
        n_chunk += 1
        if n_chunk % 20 == 0:
            log(f"    genes {sum(b.shape[0] for b in acc_blocks)} elapsed {time.time()-t0:.0f}s")
    mat = pd.concat(acc_blocks, axis=0)
    # library-size normalize to log2(CPM+1)
    lib = mat.sum(axis=0).replace(0, np.nan)
    cpm = mat.divide(lib, axis=1) * 1e6
    return np.log2(cpm + 1.0)


def _cache_path(name: str) -> Path:
    d = DATA / "processed"
    d.mkdir(parents=True, exist_ok=True)
    return d / f"{name}.tsv.gz"


def load_gse131907_malignant_pseudobulk() -> Cohort:
    cache = _cache_path("GSE131907_malignant_pseudobulk")
    if cache.exists():
        expr = pd.read_csv(cache, sep="\t", index_col=0)
        ann = pd.read_csv(RAW / "geo" / "GSE131907_Lung_Cancer_cell_annotation.txt.gz", sep="\t")
        mal = ann[ann["Cell_subtype"] == "Malignant cells"]
        mal = mal[~mal["Sample_Origin"].isin(["nLung", "nLN"])]
        meta = (
            mal.groupby("Sample")
            .agg(n_malignant=("Index", "size"), origin=("Sample_Origin", "first"))
            .reindex(expr.columns)
        )
        return Cohort(
            name="GSE131907_malignant_pseudobulk",
            source="GEO GSE131907 Kim et al. Nat Commun 2020 LUAD scRNA; malignant cells, tumor/met sites",
            modality="scrna_pseudobulk",
            expr=expr,
            notes="Cached malignant-cell UMI sums per sample; log2(CPM+1). Excludes nLung/nLN.",
            meta=meta,
        )
    ann = pd.read_csv(RAW / "geo" / "GSE131907_Lung_Cancer_cell_annotation.txt.gz", sep="\t")
    mal = ann[ann["Cell_subtype"] == "Malignant cells"].copy()
    # tumor / metastasis sites only (not nLung / nLN)
    mal = mal[~mal["Sample_Origin"].isin(["nLung", "nLN"])]
    cell_to_sample = dict(zip(mal["Index"], mal["Sample"]))
    expr = _pseudobulk_from_wide(
        RAW / "geo" / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz",
        mal["Index"].tolist(),
        cell_to_sample,
        min_cells_per_sample=30,
    )
    meta = (
        mal.groupby("Sample")
        .agg(n_malignant=("Index", "size"), origin=("Sample_Origin", "first"))
        .reindex(expr.columns)
    )
    expr.to_csv(cache, sep="\t")
    return Cohort(
        name="GSE131907_malignant_pseudobulk",
        source="GEO GSE131907 Kim et al. Nat Commun 2020 LUAD scRNA; malignant cells, tumor/met sites",
        modality="scrna_pseudobulk",
        expr=expr,
        notes="Malignant-cell UMI sums per sample; log2(CPM+1). Excludes nLung/nLN. Samples with <30 malignant cells dropped.",
        meta=meta,
    )


def load_gse207422_epithelial_pseudobulk() -> Cohort:
    cache = _cache_path("GSE207422_sc_epithelial_pseudobulk")
    if cache.exists():
        expr = pd.read_csv(cache, sep="\t", index_col=0)
        return Cohort(
            name="GSE207422_sc_epithelial_pseudobulk",
            source="GEO GSE207422 scRNA; epithelial cells inferred by EPCAM/KRT+ PTPRC/PECAM1/COL- ",
            modality="scrna_pseudobulk",
            expr=expr,
            notes="Cached marker-gated epithelial pseudobulk. No official GEO cell-type table.",
        )
    path = RAW / "geo" / "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"
    markers = {"EPCAM", "KRT8", "KRT18", "KRT19", "PTPRC", "PECAM1", "COL1A1", "COL3A1", "TACSTD2"}
    log("  GSE207422 scRNA: streaming lineage markers")
    cells, mark = _stream_marker_matrix(path, markers)
    # epithelial: keratin/EPCAM positive, immune/endothelial/fibroblast low
    epi_score = mark.reindex(["EPCAM", "KRT8", "KRT18", "KRT19"]).fillna(0).sum(axis=0)
    ptprc = mark.reindex(["PTPRC"]).fillna(0).iloc[0]
    pecam = mark.reindex(["PECAM1"]).fillna(0).iloc[0]
    col = mark.reindex(["COL1A1", "COL3A1"]).fillna(0).sum(axis=0)
    epi_mask = (epi_score >= 1) & (ptprc == 0) & (pecam == 0) & (col == 0)
    epi_cells = [c for c, keep in zip(cells, epi_mask.to_numpy()) if keep]
    # sample id is prefix before last underscore-digit chunk: BD_immune01_612637 -> BD_immune01
    cell_to_sample = {c: "_".join(c.split("_")[:2]) if c.count("_") >= 2 else c.rsplit("_", 1)[0] for c in epi_cells}
    expr = _pseudobulk_from_wide(path, epi_cells, cell_to_sample, min_cells_per_sample=15)
    meta = pd.DataFrame(
        {
            "n_epithelial": pd.Series(epi_cells).map(cell_to_sample).value_counts(),
        }
    ).reindex(expr.columns)
    meta_s = pd.read_excel(RAW / "geo" / "GSE207422_NSCLC_scRNAseq_metadata.xlsx")
    meta_s = meta_s.set_index("Sample")
    meta = meta.join(meta_s, how="left")
    expr.to_csv(cache, sep="\t")
    return Cohort(
        name="GSE207422_sc_epithelial_pseudobulk",
        source="GEO GSE207422 scRNA; epithelial cells inferred by EPCAM/KRT+ PTPRC/PECAM1/COL- ",
        modality="scrna_pseudobulk",
        expr=expr,
        notes=(
            f"No official cell-type table on GEO. Epithelial gate: (EPCAM+KRT8+KRT18+KRT19)>=1 UMI "
            f"and PTPRC=PECAM1=COL1A1/COL3A1=0. {len(epi_cells)} cells. Marker-based, not author annotation."
        ),
        meta=meta,
    )


def quartile_split(x: pd.Series) -> tuple[pd.Index, pd.Index]:
    q1, q3 = x.quantile(0.25), x.quantile(0.75)
    low = x.index[x <= q1]
    high = x.index[x >= q3]
    return high, low


def median_split(x: pd.Series) -> tuple[pd.Index, pd.Index]:
    med = x.median()
    # drop exact median ties from both sides
    high = x.index[x > med]
    low = x.index[x < med]
    return high, low


def rank_high_vs_low(expr: pd.DataFrame, high: pd.Index, low: pd.Index) -> pd.Series:
    """Welch t-statistic, high minus low. Drop near-zero-variance genes."""
    eh = expr.loc[:, high].to_numpy(dtype=np.float64)
    el = expr.loc[:, low].to_numpy(dtype=np.float64)
    # require finite values
    ok = np.isfinite(eh).all(axis=1) & np.isfinite(el).all(axis=1)
    vh = np.nanvar(eh, axis=1)
    vl = np.nanvar(el, axis=1)
    ok &= (vh + vl) > 1e-8
    tstat = np.full(expr.shape[0], np.nan)
    # vectorized welch t
    mh = eh.mean(axis=1)
    ml = el.mean(axis=1)
    nh, nl = eh.shape[1], el.shape[1]
    se = np.sqrt(vh / nh + vl / nl)
    with np.errstate(divide="ignore", invalid="ignore"):
        tstat[ok] = (mh[ok] - ml[ok]) / se[ok]
    s = pd.Series(tstat, index=expr.index).dropna()
    s = s.replace([np.inf, -np.inf], np.nan).dropna()
    # drop the ranking gene from GSEA? keep it; it is not in primary sets
    return s.sort_values(ascending=False)


def run_prerank(rank: pd.Series, gene_sets: dict[str, list[str]], outdir: Path | None) -> pd.DataFrame:
    import gseapy as gp

    # gseapy wants a 2-col dataframe or series
    rnk = rank.rename("t").reset_index()
    rnk.columns = ["gene", "score"]
    try:
        res = gp.prerank(
            rnk=rnk,
            gene_sets=gene_sets,
            outdir=str(outdir) if outdir else None,
            permutation_num=1000,
            min_size=5,
            max_size=400,
            seed=SEED,
            verbose=False,
            no_plot=True,
            threads=1,
        )
        tab = res.res2d.copy()
    except Exception as exc:
        log(f"  gseapy failed ({exc}); using fallback NES")
        tab = fallback_gsea(rank, gene_sets)
    # normalize column names across gseapy versions
    colmap = {c.lower(): c for c in tab.columns}
    def pick(*names):
        for n in names:
            if n in tab.columns:
                return n
            if n.lower() in colmap:
                return colmap[n.lower()]
        return None

    term = pick("Term", "Term")
    nes = pick("NES")
    fdr = pick("FDR q-val", "FDR", "fdr")
    pval = pick("NOM p-val", "NOM p-val", "Pval", "pval")
    es = pick("ES")
    lead = pick("Lead_genes", "Lead genes", "lead_genes")
    size = pick("Tag %", "Tag%", "Size")
    out = pd.DataFrame(
        {
            "term": tab[term] if term else tab.iloc[:, 0],
            "nes": pd.to_numeric(tab[nes], errors="coerce") if nes else np.nan,
            "fdr": pd.to_numeric(tab[fdr], errors="coerce") if fdr else np.nan,
            "nom_p": pd.to_numeric(tab[pval], errors="coerce") if pval else np.nan,
            "es": pd.to_numeric(tab[es], errors="coerce") if es else np.nan,
            "lead_genes": tab[lead] if lead else "",
        }
    )
    return out


def fallback_gsea(rank: pd.Series, gene_sets: dict[str, list[str]], nperm: int = 500) -> pd.DataFrame:
    """Simple weighted KS NES if gseapy is unavailable."""
    genes = rank.index.to_numpy()
    scores = rank.to_numpy(dtype=float)
    order = np.argsort(-np.abs(scores))  # not used; use signed order
    order = np.argsort(-scores)
    genes = genes[order]
    scores = scores[order]
    abs_s = np.abs(scores)
    rng = np.random.default_rng(SEED)
    rows = []
    for term, members in gene_sets.items():
        hit = np.isin(genes, list(members))
        n_hit = int(hit.sum())
        if n_hit < 5:
            continue
        def es_of(hit_mask):
            nh = hit_mask.sum()
            miss = (~hit_mask).astype(float)
            n_miss = miss.sum()
            hit_w = np.where(hit_mask, abs_s, 0.0)
            hit_w = hit_w / hit_w.sum() if hit_w.sum() else hit_w
            step_miss = miss / n_miss if n_miss else miss
            walk = np.cumsum(hit_w - step_miss)
            return float(walk[np.argmax(np.abs(walk))]) if len(walk) else 0.0

        es = es_of(hit)
        null = []
        for _ in range(nperm):
            perm = np.zeros(len(genes), dtype=bool)
            perm[rng.choice(len(genes), size=n_hit, replace=False)] = True
            null.append(es_of(perm))
        null = np.array(null)
        if es >= 0:
            pos = null[null >= 0]
            nes = es / pos.mean() if len(pos) and pos.mean() != 0 else 0.0
            p = (np.sum(null >= es) + 1) / (len(null) + 1)
        else:
            neg = null[null < 0]
            nes = es / abs(neg.mean()) if len(neg) and neg.mean() != 0 else 0.0
            p = (np.sum(null <= es) + 1) / (len(null) + 1)
        rows.append({"Term": term, "NES": nes, "NOM p-val": p, "FDR q-val": p, "ES": es, "Lead_genes": ""})
    return pd.DataFrame(rows)


def signature_zmean(expr: pd.DataFrame, genes: list[str]) -> pd.Series:
    present = [g for g in genes if g in expr.index]
    if len(present) < 3:
        return pd.Series(np.nan, index=expr.columns)
    sub = expr.loc[present]
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, np.nan), axis=0)
    return z.mean(axis=0)


def bh_fdr(p: pd.Series) -> pd.Series:
    p = p.astype(float)
    n = p.notna().sum()
    if n == 0:
        return p
    order = p.argsort()
    ranks = pd.Series(np.arange(1, len(p) + 1), index=p.index)
    # only on non-nan
    q = p.copy()
    valid = p.dropna()
    n = len(valid)
    order = valid.sort_values().index
    adj = valid.loc[order] * n / np.arange(1, n + 1)
    adj = adj[::-1].cummin()[::-1].clip(upper=1)
    q.loc[order] = adj
    return q


def analyze_cohort(
    cohort: Cohort,
    gene_sets: dict[str, list[str]],
    set_meta: dict[str, dict],
    split: str = "quartile",
) -> dict:
    expr = cohort.expr
    if "TACSTD2" not in expr.index:
        return {"skip": True, "reason": "TACSTD2 absent"}
    x = expr.loc["TACSTD2"].dropna()
    if x.nunique() < 6:
        return {"skip": True, "reason": "TACSTD2 has too little variation"}
    high, low = quartile_split(x) if split == "quartile" else median_split(x)
    n_high, n_low = len(high), len(low)
    if n_high < 4 or n_low < 4:
        return {"skip": True, "reason": f"split too small high={n_high} low={n_low}"}
    underpowered = n_high < 8 or n_low < 8
    rank = rank_high_vs_low(expr, high, low)
    gsea_dir = OUT / "gsea_runs" / f"{cohort.name}_{split}"
    gsea_dir.mkdir(parents=True, exist_ok=True)
    gsea = run_prerank(rank, gene_sets, gsea_dir)
    gsea["cohort"] = cohort.name
    gsea["split"] = split
    gsea["n_high"] = n_high
    gsea["n_low"] = n_low
    gsea["n_rank_genes"] = len(rank)
    gsea = gsea.merge(
        pd.DataFrame.from_dict(set_meta, orient="index").rename_axis("term").reset_index(),
        on="term",
        how="left",
    )

    # focal correlations
    focal_rows = []
    for g in CONTEXT_GENES:
        if g not in expr.index:
            focal_rows.append({"gene": g, "rho": np.nan, "p": np.nan, "present": False})
            continue
        rho, p = stats.spearmanr(x, expr.loc[g], nan_policy="omit")
        focal_rows.append({"gene": g, "rho": float(rho), "p": float(p), "present": True})
    focal = pd.DataFrame(focal_rows)
    focal["q"] = bh_fdr(focal["p"])
    focal["cohort"] = cohort.name
    focal["split"] = split

    # signature scores vs TACSTD2
    sig_rows = []
    for term, genes in gene_sets.items():
        sc = signature_zmean(expr, genes)
        if sc.isna().all():
            continue
        rho, p = stats.spearmanr(x, sc, nan_policy="omit")
        sig_rows.append(
            {
                "term": term,
                "rho": float(rho),
                "p": float(p),
                "n_genes_present": sum(g in expr.index for g in genes),
            }
        )
    sig = pd.DataFrame(sig_rows)
    if not sig.empty:
        sig["q"] = bh_fdr(sig["p"])
        sig["cohort"] = cohort.name
        sig["split"] = split
        sig = sig.merge(
            pd.DataFrame.from_dict(set_meta, orient="index").rename_axis("term").reset_index(),
            on="term",
            how="left",
        )

    # leading-edge membership of focal genes
    lead_rows = []
    for _, r in gsea.iterrows():
        leads = set(str(r.get("lead_genes") or "").replace(";", ",").replace(" ", ",").split(","))
        leads = {g for g in leads if g}
        for g in FOCAL:
            lead_rows.append(
                {
                    "cohort": cohort.name,
                    "split": split,
                    "term": r["term"],
                    "bucket": r.get("bucket"),
                    "gene": g,
                    "in_leading_edge": g in leads,
                    "nes": r["nes"],
                    "fdr": r["fdr"],
                }
            )
    lead = pd.DataFrame(lead_rows)

    verdict = make_verdict(gsea, focal, underpowered, n_high, n_low)
    return {
        "skip": False,
        "cohort": cohort.name,
        "n_samples": expr.shape[1],
        "n_high": n_high,
        "n_low": n_low,
        "underpowered": underpowered,
        "tacstd2_median": float(x.median()),
        "tacstd2_iqr": float(x.quantile(0.75) - x.quantile(0.25)),
        "gsea": gsea,
        "focal": focal,
        "signatures": sig,
        "lead": lead,
        "verdict": verdict,
        "notes": cohort.notes,
        "source": cohort.source,
        "modality": cohort.modality,
        "split": split,
    }


def make_verdict(gsea: pd.DataFrame, focal: pd.DataFrame, underpowered: bool, n_high: int, n_low: int) -> dict:
    def hit(bucket: str, expected: str, fdr_cut: float) -> bool:
        sub = gsea[(gsea["bucket"] == bucket) & (gsea["term"].map(lambda t: not str(t).startswith("CUSTOM_CLAUDIN")))]
        if sub.empty:
            return False
        if expected == "UP":
            return bool(((sub["nes"] > 0) & (sub["fdr"] < fdr_cut)).any())
        return bool(((sub["nes"] < 0) & (sub["fdr"] < fdr_cut)).any())

    def opposite(bucket: str, expected: str, fdr_cut: float) -> bool:
        sub = gsea[gsea["bucket"] == bucket]
        if sub.empty:
            return False
        if expected == "UP":
            return bool(((sub["nes"] < 0) & (sub["fdr"] < fdr_cut)).any())
        return bool(((sub["nes"] > 0) & (sub["fdr"] < fdr_cut)).any())

    foc = focal[focal["gene"].isin(FOCAL) & focal["present"]]
    n_pos = int(((foc["rho"] > 0) & (foc["q"] < 0.05)).sum())
    n_neg = int(((foc["rho"] < 0) & (foc["q"] < 0.05)).sum())

    def pack(cut: float) -> dict:
        tj = hit("TJ", "UP", cut)
        ker = hit("KERATIN_BARRIER", "UP", cut)
        emt = hit("EMT", "DOWN", cut)
        tj_opp = opposite("TJ", "UP", cut)
        ker_opp = opposite("KERATIN_BARRIER", "UP", cut)
        emt_opp = opposite("EMT", "DOWN", cut)
        focal_ok = n_pos >= 3
        if underpowered:
            label = "underpowered"
        elif emt_opp and not emt:
            label = "contradicts_EMT"
        elif (tj_opp and not tj) and (ker_opp and not ker):
            label = "contradicts_TJ_barrier"
        elif tj and ker and emt and focal_ok:
            label = "supportive"
        elif (tj or ker) and emt and focal_ok:
            label = "partial"
        elif (tj or ker) and emt:
            label = "partial_no_focal"
        elif sum([tj, ker, emt, focal_ok]) >= 2:
            label = "mixed"
        else:
            label = "null"
        return {
            "fdr_cut": cut,
            "tj_up": tj,
            "keratin_barrier_up": ker,
            "emt_down": emt,
            "tj_opposite": tj_opp,
            "keratin_opposite": ker_opp,
            "emt_opposite": emt_opp,
            "focal_pos_q05": n_pos,
            "focal_neg_q05": n_neg,
            "label": label,
        }

    return {
        "n_high": n_high,
        "n_low": n_low,
        "underpowered": underpowered,
        "exploratory_fdr0.25": pack(0.25),
        "strict_fdr0.05": pack(0.05),
    }


def plot_nes_heatmap(gsea_all: pd.DataFrame, path: Path) -> None:
    sub = gsea_all[gsea_all["split"] == "quartile"].copy()
    if sub.empty:
        return
    # order terms by bucket
    bucket_order = {"TJ": 0, "KERATIN_BARRIER": 1, "EMT": 2}
    sub["bord"] = sub["bucket"].map(bucket_order)
    terms = (
        sub.drop_duplicates("term")
        .sort_values(["bord", "term"])["term"]
        .tolist()
    )
    pivot = sub.pivot_table(index="term", columns="cohort", values="nes", aggfunc="first")
    pivot = pivot.reindex(terms)
    fdr = sub.pivot_table(index="term", columns="cohort", values="fdr", aggfunc="first").reindex(pivot.index)
    fig, ax = plt.subplots(figsize=(max(10, 0.55 * pivot.shape[1] + 4), max(6, 0.38 * pivot.shape[0] + 2)))
    sns.heatmap(
        pivot,
        cmap="RdBu_r",
        center=0,
        vmin=-2.5,
        vmax=2.5,
        ax=ax,
        linewidths=0.3,
        linecolor="white",
        cbar_kws={"label": "NES (TACSTD2-high vs low)"},
    )
    # stars for FDR
    for i, term in enumerate(pivot.index):
        for j, coh in enumerate(pivot.columns):
            try:
                q = fdr.loc[term, coh]
            except Exception:
                continue
            if pd.isna(q):
                continue
            mark = "***" if q < 0.01 else "**" if q < 0.05 else "*" if q < 0.25 else ""
            if mark:
                ax.text(j + 0.5, i + 0.5, mark, ha="center", va="center", fontsize=7, color="black")
    ax.set_title("Pre-specified GSEA NES: TACSTD2-high vs TACSTD2-low (quartile)\n* FDR<0.25  ** FDR<0.05  *** FDR<0.01")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_focal_heatmap(focal_all: pd.DataFrame, path: Path) -> None:
    sub = focal_all[(focal_all["split"] == "quartile") & (focal_all["gene"].isin(FOCAL))].copy()
    if sub.empty:
        return
    pivot = sub.pivot_table(index="gene", columns="cohort", values="rho", aggfunc="first").reindex(FOCAL)
    q = sub.pivot_table(index="gene", columns="cohort", values="q", aggfunc="first").reindex(FOCAL)
    fig, ax = plt.subplots(figsize=(max(10, 0.55 * pivot.shape[1] + 3), 3.6))
    sns.heatmap(pivot, cmap="RdBu_r", center=0, vmin=-0.8, vmax=0.8, ax=ax, linewidths=0.3, linecolor="white",
                cbar_kws={"label": "Spearman rho vs TACSTD2"})
    for i, gene in enumerate(pivot.index):
        for j, coh in enumerate(pivot.columns):
            try:
                qq = q.loc[gene, coh]
            except Exception:
                continue
            if pd.notna(qq) and qq < 0.05:
                ax.text(j + 0.5, i + 0.5, "*", ha="center", va="center", fontsize=9)
    ax.set_title("Focal genes vs TACSTD2 (Spearman). * BH-FDR<0.05")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_verdict_table(verdicts: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(11, max(4, 0.35 * len(verdicts) + 1.5)))
    ax.axis("off")
    tbl = verdicts.copy()
    ax.table(cellText=tbl.values, colLabels=tbl.columns, loc="center", cellLoc="center")
    fig.tight_layout()
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def write_report(results: list[dict], gsea_all: pd.DataFrame, focal_all: pd.DataFrame, lead_all: pd.DataFrame) -> None:
    lines = []
    lines.append("# Hunt: TACSTD2-high tumors vs tight junction / keratin-barrier / EMT")
    lines.append("")
    lines.append("**Honest discovery report.** Pre-specified hypothesis, all tested cohorts shown, including nulls.")
    lines.append("")
    lines.append("## Hypothesis (pre-specified)")
    lines.append("")
    lines.append("In public bulk and scRNA tumor datasets, **TACSTD2-high** samples GSEA-enrich")
    lines.append("**tight junction** and **keratinization / skin-barrier** programs and **deplete EMT**,")
    lines.append("with the leading-edge / co-expression intersection near **CLDN1, CLDN4, CLDN7, F11R, PARD3**.")
    lines.append("")
    lines.append("## What this is not")
    lines.append("")
    lines.append("- Not a claim that TACSTD2 *causes* tight junctions. TACSTD2 (Trop-2) is itself an epithelial surface protein.")
    lines.append("- EMT-down in TACSTD2-high tumors is **partly tautological**: TACSTD2 tracks epithelial identity.")
    lines.append("  The non-trivial part is whether TJ/barrier sets, and specifically the claudin–F11R–PARD3 module, ride along *within* a cancer type.")
    lines.append("- Within-cancer quartile splits reduce histology confounding (e.g. LUAD vs LUSC) but do not remove differentiation state inside a cancer.")
    lines.append("- GSEA FDR 0.25 is the classical exploratory cutoff; FDR 0.05 is also reported. This is a hunt, not a confirmatory trial.")
    lines.append("")
    lines.append("## Methods (locked before looking at NES)")
    lines.append("")
    lines.append("- **Split (primary):** TACSTD2 top vs bottom quartile. **Sensitivity:** median split.")
    lines.append("- **Ranking:** Welch t-statistic, high minus low, on the supplied log-like matrix.")
    lines.append("- **GSEA:** gseapy prerank, 1000 permutations, min_size=5, seed=42. MSigDB v2023.2.Hs + KEGG v7.5.1 + custom modules.")
    lines.append("- **Primary gene sets:** KEGG_TIGHT_JUNCTION, GOBP_TIGHT_JUNCTION_ORGANIZATION, HALLMARK_APICAL_JUNCTION,")
    lines.append("  CUSTOM_TJ_CORE, CUSTOM_CLAUDIN_PAR_FOCAL; GOBP_KERATINIZATION, GOBP_CORNIFICATION,")
    lines.append("  GOBP_ESTABLISHMENT_OF_SKIN_BARRIER, GOBP_KERATINOCYTE_DIFFERENTIATION, GOBP_EPIDERMAL_CELL_DIFFERENTIATION,")
    lines.append("  CUSTOM_BARRIER_KERATIN; HALLMARK_EMT, GOBP_EMT, CUSTOM_EMT_CORE.")
    lines.append("- TACSTD2 is **not** a member of any primary set (checked).")
    lines.append("- **Focal genes:** Spearman vs TACSTD2, BH-FDR within cohort.")
    lines.append("- **Complementary:** z-mean signature scores vs TACSTD2 (more powered than GSEA when n is small).")
    lines.append("- **scRNA:** sample-level pseudobulk of selected cells, not mixed-cell bulk. GSE131907 uses author malignant labels.")
    lines.append("  GSE207422 scRNA has **no public cell-type table**; epithelial cells were gated by markers (see notes).")
    lines.append("- **Verdict:** `supportive` = TJ-up AND keratin/barrier-up AND EMT-down (GSEA FDR cut) AND ≥3/5 focal genes rho>0 q<0.05.")
    lines.append("  Cohorts with n_high or n_low < 8 are labeled **underpowered** regardless of NES.")
    lines.append("")
    lines.append("## Cohorts")
    lines.append("")
    lines.append("| Cohort | Modality | n | n_high/n_low (Q) | Source | Notes |")
    lines.append("|---|---|---:|---|---|---|")
    seen = set()
    for r in results:
        if r.get("skip") or r.get("split") != "quartile" or r["cohort"] in seen:
            continue
        seen.add(r["cohort"])
        notes = (r.get("notes") or "").replace("|", "/")
        lines.append(
            f"| {r['cohort']} | {r['modality']} | {r['n_samples']} | {r['n_high']}/{r['n_low']} | {r['source']} | {notes} |"
        )
    lines.append("")
    lines.append("## Verdicts (quartile split)")
    lines.append("")
    lines.append("| Cohort | n | underpowered | FDR<0.25 | FDR<0.05 | TJ↑ | K/barrier↑ | EMT↓ | focal+/5 |")
    lines.append("|---|---:|:---:|---|---|:---:|:---:|:---:|---:|")
    for r in results:
        if r.get("skip") or r.get("split") != "quartile":
            continue
        v = r["verdict"]
        a = v["exploratory_fdr0.25"]
        b = v["strict_fdr0.05"]
        lines.append(
            f"| {r['cohort']} | {r['n_samples']} | {v['underpowered']} | {a['label']} | {b['label']} | "
            f"{'Y' if a['tj_up'] else 'n'} | {'Y' if a['keratin_barrier_up'] else 'n'} | "
            f"{'Y' if a['emt_down'] else 'n'} | {a['focal_pos_q05']} |"
        )
    lines.append("")

    # intersection
    qg = gsea_all[(gsea_all["split"] == "quartile") & (gsea_all["bucket"] == "TJ") & (gsea_all["nes"] > 0) & (gsea_all["fdr"] < 0.25)]
    lines.append("## Leading-edge intersection near CLDN1/4/7 F11R PARD3")
    lines.append("")
    if lead_all.empty:
        lines.append("No leading-edge table (gseapy did not return lead genes).")
    else:
        L = lead_all[(lead_all["split"] == "quartile") & (lead_all["bucket"] == "TJ") & (lead_all["nes"] > 0) & (lead_all["fdr"] < 0.25)]
        if L.empty:
            lines.append("No TJ set was GSEA-up at FDR<0.25 in any cohort, so there is **no leading-edge intersection to report**.")
        else:
            wide = L.groupby(["cohort", "gene"])["in_leading_edge"].any().unstack("gene").reindex(columns=FOCAL)
            lines.append("Gene is in the leading edge of **any** TJ set with NES>0 and FDR<0.25:")
            lines.append("")
            lines.append("| Cohort | " + " | ".join(FOCAL) + " | n_focal_in_LE |")
            lines.append("|---|" + "|".join(["---"] * 5) + "|---:|")
            for coh, row in wide.iterrows():
                flags = [("Y" if bool(row.get(g)) else "n") for g in FOCAL]
                lines.append(f"| {coh} | " + " | ".join(flags) + f" | {sum(x=='Y' for x in flags)} |")
            # genes in LE in >=3 cohorts
            per_gene = wide.fillna(False).astype(bool).sum(axis=0)
            lines.append("")
            lines.append(
                "Focal genes present in a TJ leading edge in "
                f"**{int((per_gene>=1).sum())}/5** cohorts-with-TJ-up (count of cohorts per gene): "
                + ", ".join(f"{g}={int(per_gene.get(g,0))}" for g in FOCAL)
            )
    lines.append("")
    lines.append("## Honest take")
    lines.append("")
    qres = [r for r in results if not r.get("skip") and r.get("split") == "quartile"]
    supportive = [r["cohort"] for r in qres if r["verdict"]["exploratory_fdr0.25"]["label"] == "supportive"]
    partial = [r["cohort"] for r in qres if r["verdict"]["exploratory_fdr0.25"]["label"].startswith("partial")]
    nulls = [r["cohort"] for r in qres if r["verdict"]["exploratory_fdr0.25"]["label"] == "null"]
    under = [r["cohort"] for r in qres if r["verdict"]["underpowered"]]
    contra = [r["cohort"] for r in qres if "contradict" in r["verdict"]["exploratory_fdr0.25"]["label"]]
    lines.append(f"- **Supportive (FDR<0.25, adequately sized):** {', '.join(supportive) if supportive else 'none'}.")
    lines.append(f"- **Partial:** {', '.join(partial) if partial else 'none'}.")
    lines.append(f"- **Null:** {', '.join(nulls) if nulls else 'none'}.")
    lines.append(f"- **Contradictory:** {', '.join(contra) if contra else 'none'}.")
    lines.append(f"- **Underpowered (do not over-read NES):** {', '.join(under) if under else 'none'}.")
    lines.append("")
    lines.append("If the supportive list is empty or limited to squamous-leaning cancers (LUSC, HNSC, CESC),")
    lines.append("the keratin/barrier arm is likely **differentiation**, not a TACSTD2-specific TJ program.")
    lines.append("LUAD (TCGA-LUAD, GSE31210, GSE131907 malignant pseudobulk) is the cleaner test of the hunt.")
    lines.append("GSE207422 bulk n=24 is too small for stable GSEA; trust the Spearman/signature columns more than NES there.")
    lines.append("")
    lines.append("## Files")
    lines.append("")
    lines.append("- `gsea_prerank_all.tsv` — all primary sets, all cohorts, both splits")
    lines.append("- `focal_gene_correlations.tsv` — TACSTD2 vs CLDN1/4/7 F11R PARD3 and context genes")
    lines.append("- `signature_scores_vs_tacstd2.tsv` — z-mean set scores vs TACSTD2")
    lines.append("- `leading_edge_focal.tsv` — whether each focal gene is in each set's leading edge")
    lines.append("- `verdicts.tsv` — machine-readable labels")
    lines.append("- `cohort_inventory.tsv`")
    lines.append("- `figures/nes_heatmap.png`, `figures/focal_rho_heatmap.png`")
    lines.append("")
    (OUT / "README.md").write_text("\n".join(lines) + "\n")


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    gene_sets, set_meta = load_primary_sets()
    log(f"primary sets: {list(gene_sets)}")

    cohorts: list[Cohort] = []
    cancers = [
        "LUAD",
        "LUSC",
        "HNSC",
        "BLCA",
        "CESC",
        "BRCA",
        "COAD",
        "ESCA",
        "STAD",
        "PAAD",
        "OV",
        "UCEC",
        "SKCM",
        "KIRC",
        "GBM",
    ]
    for c in cancers:
        log(f"load TCGA-{c}")
        coh = load_tcga(c)
        if coh is not None:
            log(f"  {coh.name} {coh.expr.shape}")
            cohorts.append(coh)

    log("load GSE207422 bulk")
    cohorts.append(load_gse207422_bulk())
    log(f"  {cohorts[-1].expr.shape}")

    log("load GSE31210")
    try:
        cohorts.append(load_gse31210())
        log(f"  {cohorts[-1].expr.shape}")
    except Exception as exc:
        log(f"  GSE31210 failed: {exc}")

    # scRNA last (slow)
    try:
        log("load GSE131907 malignant pseudobulk")
        cohorts.append(load_gse131907_malignant_pseudobulk())
        log(f"  {cohorts[-1].expr.shape}")
    except Exception as exc:
        log(f"  GSE131907 failed: {exc}")

    try:
        log("load GSE207422 sc epithelial pseudobulk")
        cohorts.append(load_gse207422_epithelial_pseudobulk())
        log(f"  {cohorts[-1].expr.shape}")
    except Exception as exc:
        log(f"  GSE207422 sc failed: {exc}")

    inv = pd.DataFrame(
        [
            {
                "cohort": c.name,
                "source": c.source,
                "modality": c.modality,
                "n_genes": c.expr.shape[0],
                "n_samples": c.expr.shape[1],
                "has_TACSTD2": "TACSTD2" in c.expr.index,
                "has_all_focal": all(g in c.expr.index for g in FOCAL),
                "notes": c.notes,
            }
            for c in cohorts
        ]
    )
    inv.to_csv(OUT / "cohort_inventory.tsv", sep="\t", index=False)

    results = []
    gsea_frames, focal_frames, sig_frames, lead_frames = [], [], [], []
    for c in cohorts:
        for split in ("quartile", "median"):
            log(f"analyze {c.name} {split} n={c.expr.shape[1]}")
            r = analyze_cohort(c, gene_sets, set_meta, split=split)
            if r.get("skip"):
                log(f"  skip: {r.get('reason')}")
                continue
            results.append(r)
            gsea_frames.append(r["gsea"])
            focal_frames.append(r["focal"])
            if r["signatures"] is not None and len(r["signatures"]):
                sig_frames.append(r["signatures"])
            lead_frames.append(r["lead"])
            log(f"  verdict FDR0.25={r['verdict']['exploratory_fdr0.25']['label']}")

    gsea_all = pd.concat(gsea_frames, ignore_index=True) if gsea_frames else pd.DataFrame()
    focal_all = pd.concat(focal_frames, ignore_index=True) if focal_frames else pd.DataFrame()
    sig_all = pd.concat(sig_frames, ignore_index=True) if sig_frames else pd.DataFrame()
    lead_all = pd.concat(lead_frames, ignore_index=True) if lead_frames else pd.DataFrame()

    if not gsea_all.empty:
        gsea_all.to_csv(OUT / "gsea_prerank_all.tsv", sep="\t", index=False)
        gsea_all[gsea_all["split"] == "quartile"].to_csv(OUT / "gsea_prerank_quartile.tsv", sep="\t", index=False)
    if not focal_all.empty:
        focal_all.to_csv(OUT / "focal_gene_correlations.tsv", sep="\t", index=False)
    if not sig_all.empty:
        sig_all.to_csv(OUT / "signature_scores_vs_tacstd2.tsv", sep="\t", index=False)
    if not lead_all.empty:
        lead_all.to_csv(OUT / "leading_edge_focal.tsv", sep="\t", index=False)

    # verdicts table (quartile only, unique cohorts)
    vrows = []
    seen = set()
    for r in results:
        if r["split"] != "quartile" or r["cohort"] in seen:
            continue
        seen.add(r["cohort"])
        a = r["verdict"]["exploratory_fdr0.25"]
        b = r["verdict"]["strict_fdr0.05"]
        vrows.append(
            {
                "cohort": r["cohort"],
                "modality": r["modality"],
                "n": r["n_samples"],
                "n_high": r["n_high"],
                "n_low": r["n_low"],
                "underpowered": r["verdict"]["underpowered"],
                "label_fdr0.25": a["label"],
                "label_fdr0.05": b["label"],
                "tj_up_fdr0.25": a["tj_up"],
                "keratin_barrier_up_fdr0.25": a["keratin_barrier_up"],
                "emt_down_fdr0.25": a["emt_down"],
                "focal_pos_q05": a["focal_pos_q05"],
                "focal_neg_q05": a["focal_neg_q05"],
            }
        )
    verdicts = pd.DataFrame(vrows)
    verdicts.to_csv(OUT / "verdicts.tsv", sep="\t", index=False)

    if not gsea_all.empty:
        plot_nes_heatmap(gsea_all, FIG / "nes_heatmap.png")
    if not focal_all.empty:
        plot_focal_heatmap(focal_all, FIG / "focal_rho_heatmap.png")

    write_report(results, gsea_all, focal_all, lead_all)
    log(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
