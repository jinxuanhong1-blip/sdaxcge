#!/usr/bin/env python3
"""CLDN4-high malignant GRN proxy on public GSE131907.

Additive SCENIC/pySCENIC proxy (Pearson + public priors + AUCell) on
author-labeled malignant cells, split CLDN4-high vs low.

pySCENIC cisTarget is not run. No ChIP peaks are invented. Regulons are:

  1. public curated TF–target priors (TRRUST / DoRothEA / CollecTRI)
  2. Pearson co-expression in CLDN4-high cells (top-N)
  3. high-specific edges: r_high >= 0.15 and r_low < 0.10

A10 ELF3–CLDN4 bulk RNA is taken as given and is not re-tested as a discovery.
CLDN4 is held out of every regulon (it is the split gene).
CLDN4-low edges are used only to mark high-specific targets.
"""
from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
REPO = ROOT.parent.parent
RES = REPO / "results" / "gse131907_scenic_cldn4"
RES.mkdir(parents=True, exist_ok=True)

GEO_CANDIDATES = [
    Path("/tmp/gse131907_scenic_cldn4/geo"),
    Path("/tmp/gse131907"),
    REPO / "data" / "gse131907_scenic_cldn4" / "GSE131907",
    Path("data/gse131907_scenic_cldn4/GSE131907"),
]
PRIORS = ROOT / "resources" / "tf_targets_public.tsv"
TF_LIST = ROOT / "resources" / "human_tfs_from_priors.txt"
FINDING = ROOT / "FINDING.md"

FOCUS = ["ELF3", "GRHL1", "GRHL2", "KLF4", "KLF5",
         "OVOL1", "OVOL2", "TFAP2A", "SP1", "TP63"]
CONTROLS = ["NKX2-1", "SOX2", "FOXA1"]
GIVEN_TFS = {"ELF3"}
HOLD_OUT = {"CLDN4"}
MALIGNANT_SUBTYPES = {"Malignant cells", "tS1", "tS2", "tS3"}

MIN_UMI = 200
MIN_CELLS = 20
MIN_SAMPLES = 6
HVG_N = 4000
MIN_DET = 0.02
TF_MIN_DET = 0.05
PEARSON_TOP = 50
PEARSON_MIN_R = 0.10
HIGH_SPEC_R = 0.15
LOW_SPEC_R = 0.10
AUC_THR = 0.05
TOP_SCREEN = 25


def find_geo() -> Path:
    for p in GEO_CANDIDATES:
        if (p / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz").exists() and (
            p / "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
        ).exists():
            return p
    raise SystemExit("GSE131907 UMI/annotation not found. Run scripts/download_gse131907.py")


def log1p_cp10k(counts: np.ndarray, total: np.ndarray) -> np.ndarray:
    tot = np.asarray(total, float)
    tot[tot <= 0] = np.nan
    return np.log1p(np.asarray(counts, float) / tot[:, None] * 1e4)


def spear(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    n = int(m.sum())
    if n < 4:
        return dict(rho=np.nan, p=np.nan, n=n)
    r, p = stats.spearmanr(x[m], y[m])
    return dict(rho=float(r), p=float(p), n=n)


def aucell(expr: np.ndarray, member: np.ndarray, auc_threshold: float = AUC_THR) -> np.ndarray:
    """Aibar-style AUCell linear recovery. Not R AUCell binary; not binding."""
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


def gene_name(raw: bytes) -> str:
    g = raw.decode("ascii", errors="replace")
    return g.split(".")[0]


def stream_pass1(matrix: Path, keep_idx: np.ndarray, always: set[str]):
    """Library sizes + gene stats on author-malignant columns; store always-keep genes."""
    keep_idx = np.asarray(keep_idx, dtype=int)
    n_keep = int(keep_idx.size)
    with gzip.open(matrix, "rb") as fh:
        header = fh.readline().decode("ascii").rstrip("\n")
        n = len(header.split("\t")) - 1
        total = np.zeros(n_keep, dtype=np.float64)
        kept: dict[str, np.ndarray] = {}
        stats_rows = []
        n_genes = 0
        for raw in fh:
            n_genes += 1
            tab = raw.find(b"\t")
            gene = gene_name(raw[:tab])
            arr = np.fromstring(raw[tab + 1 :].rstrip(b"\r\n"), sep="\t", dtype=np.float32)
            if arr.size != n:
                raise ValueError(f"{gene}: {arr.size} != {n}")
            sub = arr[keep_idx]
            total += sub
            det = float((sub > 0).mean())
            mu = float(sub.mean())
            var = float(sub.astype(np.float64).var())
            stats_rows.append((gene, mu, var, det))
            if gene in always:
                kept[gene] = sub.copy()
            if n_genes % 2000 == 0:
                print(f"  pass1 genes={n_genes} always_kept={len(kept)}", flush=True)
    info = pd.DataFrame(stats_rows, columns=["gene", "mean", "var", "detection"])
    print(f"[pass1] done genes={n_genes} always_kept={len(kept)} malig_cols={n_keep}", flush=True)
    return total, kept, info, n_genes


def stream_pass2(matrix: Path, keep_idx: np.ndarray, want: set[str]):
    keep_idx = np.asarray(keep_idx, dtype=int)
    n_keep = int(keep_idx.size)
    with gzip.open(matrix, "rb") as fh:
        header = fh.readline().decode("ascii").rstrip("\n")
        n = len(header.split("\t")) - 1
        genes: list[str] = []
        mats: list[np.ndarray] = []
        n_genes = 0
        for raw in fh:
            n_genes += 1
            tab = raw.find(b"\t")
            gene = gene_name(raw[:tab])
            if gene not in want:
                continue
            arr = np.fromstring(raw[tab + 1 :].rstrip(b"\r\n"), sep="\t", dtype=np.float32)
            if arr.size != n:
                raise ValueError(f"{gene}: {arr.size} != {n}")
            genes.append(gene)
            mats.append(arr[keep_idx].copy())
            if len(genes) % 500 == 0:
                print(f"  pass2 stored={len(genes)} scanned={n_genes}", flush=True)
    print(f"[pass2] stored {len(genes)} / {n_genes} genes for {n_keep} cells", flush=True)
    X = np.vstack(mats).T.astype(np.float32)
    return X, genes


def load_priors() -> pd.DataFrame:
    df = pd.read_csv(PRIORS, sep="\t")
    df["tf"] = df["tf"].astype(str)
    df["target"] = df["target"].astype(str)
    return df


def load_tf_list() -> set[str]:
    return {ln.strip() for ln in TF_LIST.read_text().splitlines() if ln.strip()}


def prior_targets(priors: pd.DataFrame, tf: str, genes: set[str]) -> list[str]:
    hits = priors.loc[priors.tf == tf, "target"].unique().tolist()
    return sorted(g for g in hits if g in genes and g != tf and g not in HOLD_OUT)


def pearson_block(logX: np.ndarray, genes: list[str], tfs: list[str]) -> pd.DataFrame:
    gi = {g: i for i, g in enumerate(genes)}
    present = [t for t in tfs if t in gi]
    if not present:
        return pd.DataFrame(columns=["tf", "target", "pearson_r"])
    X = np.asarray(logX, float)
    ok = np.isfinite(X).all(axis=1)
    X = X[ok]
    if X.shape[0] < 20:
        return pd.DataFrame(columns=["tf", "target", "pearson_r"])
    X = X - X.mean(axis=0, keepdims=True)
    sd = X.std(axis=0, ddof=1)
    sd[sd == 0] = np.nan
    Z = X / sd
    tf_idx = np.array([gi[t] for t in present], dtype=int)
    corr = (Z[:, tf_idx].T @ Z) / (Z.shape[0] - 1)
    rows = []
    gene_arr = np.array(genes)
    for i, tf in enumerate(present):
        r = corr[i]
        for j, g in enumerate(gene_arr):
            if g == tf:
                continue
            val = float(r[j])
            if not np.isfinite(val):
                continue
            rows.append((tf, g, val))
    return pd.DataFrame(rows, columns=["tf", "target", "pearson_r"])


def top_regulon(pr: pd.DataFrame, tf: str) -> pd.DataFrame:
    sub = pr.loc[(pr.tf == tf) & (~pr.target.isin(HOLD_OUT))].copy()
    if sub.empty:
        return sub
    sub = sub.sort_values("pearson_r", ascending=False)
    top = sub.loc[sub.pearson_r >= PEARSON_MIN_R].head(PEARSON_TOP)
    if len(top) < 10:
        top = sub.loc[sub.pearson_r > 0].head(PEARSON_TOP)
    top = top.reset_index(drop=True)
    top["rank"] = np.arange(1, len(top) + 1)
    return top


def wilcoxon_paired(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    m = np.isfinite(a) & np.isfinite(b)
    n = int(m.sum())
    if n < 4:
        return dict(stat=np.nan, p=np.nan, n=n, delta=np.nan)
    try:
        w = stats.wilcoxon(a[m], b[m], zero_method="wilcox", alternative="two-sided")
        return dict(stat=float(w.statistic), p=float(w.pvalue), n=n,
                    delta=float(np.median(a[m] - b[m])))
    except ValueError:
        return dict(stat=np.nan, p=np.nan, n=n, delta=float(np.median(a[m] - b[m])))


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
    n = len(fields.get("title", []))
    return pd.DataFrame({k: v for k, v in fields.items() if len(v) == n})


def md_table(df: pd.DataFrame, cols: list[str]) -> str:
    use = df.loc[:, [c for c in cols if c in df.columns]].copy()
    lines = ["| " + " | ".join(use.columns) + " |",
             "| " + " | ".join("---" for _ in use.columns) + " |"]
    for rec in use.itertuples(index=False):
        cells = []
        for v in rec:
            if v is None or (isinstance(v, float) and not np.isfinite(v)):
                cells.append("")
            else:
                cells.append(str(v).replace("|", "/"))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def write_finding(regulons: pd.DataFrame, n_tab: pd.DataFrame, paired_df: pd.DataFrame,
                  summary: dict) -> None:
    focus = regulons.loc[
        regulons.tf.isin(FOCUS)
        & regulons.kind.isin(["public_prior", "pearson_cldn4_high", "high_specific", "prior_AND_pearson_high"])
    ].copy()
    show = focus[["regulon", "tf", "kind", "n_targets", "n_high_specific", "mean_r", "elf3_given"]].copy()
    n_show = n_tab.copy()
    paired_show = paired_df.copy() if len(paired_df) else pd.DataFrame()
    if len(paired_show):
        keep = [c for c in ["regulon", "kind", "n_samples", "delta_median_high_minus_low", "p", "fdr"]
                if c in paired_show.columns]
        paired_show = paired_show[keep]
        if "delta_median_high_minus_low" in paired_show:
            paired_show["delta_median_high_minus_low"] = paired_show["delta_median_high_minus_low"].map(
                lambda x: f"{x:.4g}" if pd.notna(x) else ""
            )
        if "p" in paired_show:
            paired_show["p"] = paired_show["p"].map(lambda x: f"{x:.3g}" if pd.notna(x) else "")
        if "fdr" in paired_show:
            paired_show["fdr"] = paired_show["fdr"].map(lambda x: f"{x:.3g}" if pd.notna(x) else "")

    n_high = summary["n_cldn4_high"]
    n_low = summary["n_cldn4_low"]
    n_samp = summary["n_samples_primary"]
    n_paired = summary["n_samples_paired"]
    n_mal = summary["n_malignant_author"]
    body = f"""# FINDING — GSE131907 SCENIC/GRN proxy, malignant CLDN4-high vs low

Additive public scRNA GRN on **GSE131907** (Kim et al., *Nat Commun* 2020, PMID 32385277).
**A10 ELF3–CLDN4 bulk RNA is taken as given** and is not re-tested as a discovery.

This is a **lightweight SCENIC proxy** (Pearson + public priors + AUCell). Full pySCENIC cisTarget was **not** run. No ChIP peaks are invented. Pearson neighborhoods are co-expression, not binding.

**CLDN4-high regulons only.** CLDN4-low edges are used only to mark high-specific targets. CLDN4 is held out of every regulon (it is the split gene).

Primary table: [`results/gse131907_scenic_cldn4/regulons.tsv`](../../results/gse131907_scenic_cldn4/regulons.tsv).

## Honest n

Primary unit = **sample**. Cell n is labeled as cells and is not the claim (pseudoreplication). EBUS_28 and NS_07 can dominate pooled cell counts.

| item | n | note |
| --- | ---: | --- |
"""
    for rec in n_show.itertuples(index=False):
        body += f"| {rec.item} | {rec.n} | {rec.note} |\n"

    body += f"""
Author-labeled malignant cells = **{n_mal}**. Primary GRN uses samples with ≥{MIN_CELLS} such cells (**{n_samp} samples**): CLDN4-high **{n_high} cells**, CLDN4-low **{n_low} cells**. Paired prior-AUCell unit = **{n_paired} samples** with ≥{MIN_CELLS} high and low.

Unlabeled tumor epithelial cells (empty `Cell_subtype`) are **not** counted as malignant. PE / nLung / nLN contribute **0** author-malignant cells. LUNG_T09 has 5 malignant cells and is excluded from the primary sample set.

## CLDN4-high regulon table (focus TFs)

{md_table(show, ["regulon", "tf", "kind", "n_targets", "n_high_specific", "mean_r", "elf3_given"])}

Target lists (semicolon-separated) are in `regulons.tsv`. `elf3_given=True` means ELF3 is A10-given, not a discovery from this folder.

## Sample-paired public-prior AUCell (CLDN4-high − low)

"""
    if len(paired_show):
        body += md_table(paired_show, list(paired_show.columns)) + "\n\n"
    else:
        body += "_No paired samples met the ≥20 high and ≥20 low rule._\n\n"

    body += f"""## What is / is not claimed

- **Claimed:** a descriptive CLDN4-high regulon table on author-malignant GSE131907 cells, with honest sample n.
- **Not claimed:** TF binding, pySCENIC cisTarget, or a new ELF3–CLDN4 discovery (A10 is given).
- **Not claimed:** ICI / MPR association (this atlas is treatment-naive).
- **Not claimed:** cell-level p-values as the finding.

Method flags: pyscenic importable={summary["pyscenic"]["importable"]}; cisTarget_run={summary["pyscenic"]["cistarget_run"]}; chip_peaks_invented={summary["pyscenic"]["chip_peaks_invented"]}.
"""
    FINDING.write_text(body)
    print(f"[finding] wrote {FINDING}", flush=True)


def main() -> None:
    geo = find_geo()
    matrix = geo / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    ann_p = geo / "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
    series_p = geo / "GSE131907_series_matrix.txt.gz"
    print(f"[data] {matrix}", flush=True)

    try:
        import pyscenic
        pyscenic_ver = getattr(pyscenic, "__version__", "unknown")
        pyscenic_ok = True
    except Exception as exc:
        pyscenic_ver = str(exc)
        pyscenic_ok = False
    print(f"[pyscenic] importable={pyscenic_ok} version={pyscenic_ver}", flush=True)
    print("[pyscenic] cisTarget NOT run; no ChIP peaks invented.", flush=True)

    priors = load_priors()
    human_tfs = load_tf_list()
    focus_all = FOCUS + CONTROLS

    ann = pd.read_csv(ann_p, sep="\t", dtype=str)
    if "Index" not in ann.columns:
        raise SystemExit(f"annotation missing Index; columns={list(ann.columns)}")
    ann["author_malignant"] = ann["Cell_subtype"].isin(MALIGNANT_SUBTYPES)
    print(f"[ann] cells={len(ann)} author_malignant={int(ann.author_malignant.sum())}", flush=True)
    print(ann.loc[ann.author_malignant].groupby(["Sample_Origin", "Cell_subtype"]).size().to_string(),
          flush=True)

    sample_meta = pd.DataFrame()
    if series_p.exists():
        meta = parse_series_matrix(series_p)
        if "title" in meta.columns:
            sample_meta = meta.rename(columns={"title": "Sample"})
            keep_cols = [c for c in ["Sample", "geo_accession", "patient_id", "tumor_stage",
                                     "tissue_origin_abbrevation", "source_name_ch1"]
                         if c in sample_meta.columns]
            sample_meta = sample_meta[keep_cols].drop_duplicates("Sample")

    with gzip.open(matrix, "rb") as fh:
        header = fh.readline().decode("ascii").rstrip("\n")
    cells = header.split("\t")[1:]
    n_raw = len(cells)
    print(f"[matrix] cells_in_header={n_raw}", flush=True)

    ann_i = ann.set_index("Index")
    missing = [c for c in cells if c not in ann_i.index]
    if missing:
        raise SystemExit(f"{len(missing)} matrix barcodes missing from annotation (e.g. {missing[:3]})")
    obs = ann_i.reindex(cells).reset_index()
    obs = obs.rename(columns={"Index": "barcode", "Sample": "sample"})
    if sample_meta is not None and len(sample_meta):
        obs = obs.merge(sample_meta, left_on="sample", right_on="Sample", how="left")
        if "Sample" in obs.columns:
            obs = obs.drop(columns=["Sample"])

    mal_mask = obs.author_malignant.to_numpy()
    mal_idx = np.where(mal_mask)[0]
    if mal_idx.size < 50:
        raise SystemExit(f"too few author-malignant cells: {mal_idx.size}")

    always = set(focus_all) | {"CLDN4", "TACSTD2", "EPCAM"} | set(priors.target) | set(priors.tf)
    print("[pass1] library sizes + gene stats on malignant columns ...", flush=True)
    total, kept, gene_stats, n_genes_full = stream_pass1(matrix, mal_idx, always)
    qc = total >= MIN_UMI
    print(f"[qc] malignant UMI>={MIN_UMI}: {int(qc.sum())}/{mal_idx.size}", flush=True)

    mal = obs.iloc[mal_idx].reset_index(drop=True)
    mal["total_umi"] = total.astype(np.int64)
    mal["qc_umi"] = qc
    if "CLDN4" not in kept:
        raise SystemExit("CLDN4 absent from public UMI matrix")
    cldn4_counts = np.column_stack([kept["CLDN4"]])
    mal["CLDN4_log1p"] = log1p_cp10k(cldn4_counts, total)[:, 0]
    if "TACSTD2" in kept:
        mal["TACSTD2_log1p"] = log1p_cp10k(np.column_stack([kept["TACSTD2"]]), total)[:, 0]
    for tf in focus_all:
        if tf in kept:
            mal[f"{tf}_log1p"] = log1p_cp10k(np.column_stack([kept[tf]]), total)[:, 0]

    mal_qc = mal.loc[mal.qc_umi].copy()
    samp_n = mal_qc.groupby("sample").size()
    keep_samples = set(samp_n[samp_n >= MIN_CELLS].index)
    mal_qc = mal_qc.loc[mal_qc["sample"].isin(keep_samples)].copy()
    med_within = mal_qc.groupby("sample")["CLDN4_log1p"].transform("median")
    mal_qc["cldn4_high_within"] = mal_qc["CLDN4_log1p"] >= med_within
    global_med = float(mal.loc[mal.qc_umi, "CLDN4_log1p"].median())
    mal_qc["cldn4_high_global"] = mal_qc["CLDN4_log1p"] >= global_med

    mal = mal.merge(
        mal_qc[["barcode", "cldn4_high_within", "cldn4_high_global"]],
        on="barcode", how="left",
    )
    mal["in_primary_samples"] = mal["sample"].isin(keep_samples) & mal.qc_umi
    mal["cldn4_high_within"] = mal["cldn4_high_within"].fillna(False).astype(bool)
    mal["cldn4_high_global"] = mal["cldn4_high_global"].fillna(False).astype(bool)

    n_high = int((mal.in_primary_samples & mal.cldn4_high_within).sum())
    n_low = int((mal.in_primary_samples & ~mal.cldn4_high_within).sum())
    print(f"[split] within-sample CLDN4-high={n_high} low={n_low} "
          f"samples={len(keep_samples)} global_med={global_med:.3f}", flush=True)

    # HVG among malignant QC cells: use pass1 variance, require detection
    gs = gene_stats.drop_duplicates("gene", keep="first").copy()
    always_mask = gs.gene.isin(always | human_tfs)
    cand = gs.loc[(~always_mask) & (gs.detection >= MIN_DET)].sort_values("var", ascending=False)
    hvg = set(cand.head(HVG_N).gene.tolist())
    want = set(gs.loc[always_mask, "gene"]) | hvg | always
    print(f"[hvg] want={len(want)} (always/tfs + top {HVG_N} var, det>={MIN_DET})", flush=True)

    print("[pass2] extract HVG + always-keep on malignant columns ...", flush=True)
    X_counts, genes = stream_pass2(matrix, mal_idx, want)
    # drop duplicate gene names if any (keep first)
    seen = {}
    keep_g = []
    for i, g in enumerate(genes):
        if g not in seen:
            seen[g] = i
            keep_g.append(i)
    if len(keep_g) != len(genes):
        X_counts = X_counts[:, keep_g]
        genes = [genes[i] for i in keep_g]
    gene_set = set(genes)
    logX = log1p_cp10k(X_counts, mal.total_umi.to_numpy())
    print(f"[expr] author-malignant matrix {logX.shape}", flush=True)

    primary = mal.in_primary_samples.to_numpy()
    high = mal.cldn4_high_within.to_numpy() & primary
    low = (~mal.cldn4_high_within.to_numpy()) & primary
    print(f"[grn cells] high={int(high.sum())} low={int(low.sum())}", flush=True)

    det_high = (X_counts[high] > 0).mean(axis=0)
    tf_present = [
        g for g, d in zip(genes, det_high)
        if g in human_tfs and d >= TF_MIN_DET and g not in HOLD_OUT
    ]
    print(f"[tfs] detected>={TF_MIN_DET} in CLDN4-high: {len(tf_present)}", flush=True)

    print("[pearson] CLDN4-high ...", flush=True)
    pr_high = pearson_block(logX[high], genes, tf_present)
    print("[pearson] CLDN4-low (contrast only, not reported as regulons) ...", flush=True)
    pr_low = pearson_block(logX[low], genes, tf_present)
    pr_high.to_csv(RES / "pearson_cldn4_high.tsv.gz", sep="\t", index=False)
    pr_low_s = pr_low.rename(columns={"pearson_r": "pearson_r_low"})
    joined = pr_high.merge(pr_low_s, on=["tf", "target"], how="left")
    joined["high_specific"] = (
        (joined.pearson_r >= HIGH_SPEC_R)
        & (joined.pearson_r_low.fillna(0) < LOW_SPEC_R)
        & (~joined.target.isin(HOLD_OUT))
    )
    joined.to_csv(RES / "pearson_high_vs_low.tsv.gz", sep="\t", index=False)

    rows = []
    for tf in focus_all:
        prior = prior_targets(priors, tf, gene_set)
        rows.append(dict(
            regulon=f"{tf}_prior", tf=tf, kind="public_prior",
            n_targets=len(prior), n_high_specific="",
            mean_r="",
            targets=";".join(prior),
            elf3_given=tf in GIVEN_TFS,
            note="TRRUST+DoRothEA+CollecTRI union; CLDN4 held out; not ChIP",
        ))

    tf_stats = []
    for tf, sub in pr_high.groupby("tf"):
        sub = sub.loc[~sub.target.isin(HOLD_OUT)]
        top = sub.loc[sub.pearson_r >= PEARSON_MIN_R].nlargest(PEARSON_TOP, "pearson_r")
        if len(top) < 10:
            top = sub.loc[sub.pearson_r > 0].nlargest(PEARSON_TOP, "pearson_r")
        hs = joined.loc[(joined.tf == tf) & joined.high_specific]
        tf_stats.append(dict(
            tf=tf,
            n_r10=int((sub.pearson_r >= PEARSON_MIN_R).sum()),
            n_top=int(len(top)),
            mean_r_top=float(top.pearson_r.mean()) if len(top) else np.nan,
            n_high_specific=int(len(hs)),
            det_high=float(det_high[genes.index(tf)]) if tf in genes else np.nan,
        ))
    tf_stats = pd.DataFrame(tf_stats).sort_values(
        ["n_high_specific", "mean_r_top"], ascending=False
    )
    tf_stats.to_csv(RES / "tf_screen_cldn4_high.tsv", sep="\t", index=False)

    screen_tfs = tf_stats.head(TOP_SCREEN).tf.tolist() if len(tf_stats) else []
    report_tfs = list(dict.fromkeys(focus_all + screen_tfs))

    for tf in report_tfs:
        top = top_regulon(pr_high, tf)
        hs = joined.loc[(joined.tf == tf) & joined.high_specific, "target"].tolist()
        hs_in_top = [g for g in top.target.tolist() if g in set(hs)]
        rows.append(dict(
            regulon=f"{tf}_pearson_cldn4_high", tf=tf, kind="pearson_cldn4_high",
            n_targets=int(len(top)), n_high_specific=int(len(hs)),
            mean_r=f"{top.pearson_r.mean():.3f}" if len(top) else "",
            targets=";".join(top.target.tolist()),
            elf3_given=tf in GIVEN_TFS,
            note=(
                f"CLDN4-high only; top {len(top)} Pearson r>={PEARSON_MIN_R} "
                f"(or top positive); CLDN4 held out; not ChIP"
            ),
        ))
        rows.append(dict(
            regulon=f"{tf}_high_specific", tf=tf, kind="high_specific",
            n_targets=int(len(hs_in_top)), n_high_specific=int(len(hs)),
            mean_r="",
            targets=";".join(hs_in_top),
            elf3_given=tf in GIVEN_TFS,
            note=(
                f"subset of Pearson top with r_high>={HIGH_SPEC_R} and "
                f"r_low<{LOW_SPEC_R}; CLDN4-low regulons are not reported"
            ),
        ))
        prior = prior_targets(priors, tf, gene_set)
        if prior:
            inter = [g for g in prior if g in set(top.target)]
            rows.append(dict(
                regulon=f"{tf}_intersect", tf=tf, kind="prior_AND_pearson_high",
                n_targets=len(inter), n_high_specific="",
                mean_r="",
                targets=";".join(inter),
                elf3_given=tf in GIVEN_TFS,
                note="public prior ∩ CLDN4-high Pearson top; not ChIP",
            ))

    regulons = pd.DataFrame(rows)
    regulons.to_csv(RES / "regulons.tsv", sep="\t", index=False)
    print("\n=== CLDN4-high regulons (focus TFs) ===", flush=True)
    show = regulons.loc[
        regulons.tf.isin(FOCUS) & regulons.kind.isin(["public_prior", "pearson_cldn4_high"])
    ]
    print(show[["regulon", "n_targets", "n_high_specific", "mean_r", "elf3_given"]].to_string(index=False),
          flush=True)

    prior_regs = regulons.loc[regulons.kind == "public_prior"]
    auc = {}
    for rec in prior_regs.itertuples(index=False):
        members = [g for g in (rec.targets.split(";") if rec.targets else []) if g in gene_set]
        mask = np.array([g in set(members) for g in genes])
        auc[rec.regulon] = aucell(logX, mask) if mask.any() else np.full(logX.shape[0], np.nan)
        print(f"[aucell] {rec.regulon} n={int(mask.sum())} mean={np.nanmean(auc[rec.regulon]):.4f}",
              flush=True)
        mal[f"AUC_{rec.regulon}"] = auc[rec.regulon]

    paired_rows = []
    samp_rows = []
    use = mal.loc[mal.in_primary_samples].copy()
    for sample, g in use.groupby("sample"):
        hi = g.loc[g.cldn4_high_within]
        lo = g.loc[~g.cldn4_high_within]
        rec = dict(
            sample=sample,
            Sample_Origin=g.Sample_Origin.iloc[0],
            n_malignant=len(g),
            n_high=len(hi),
            n_low=len(lo),
            Cell_subtype_mode=g.Cell_subtype.mode().iloc[0] if len(g.Cell_subtype.mode()) else "",
            CLDN4_high=float(hi.CLDN4_log1p.mean()) if len(hi) else np.nan,
            CLDN4_low=float(lo.CLDN4_log1p.mean()) if len(lo) else np.nan,
        )
        if "patient_id" in g.columns:
            rec["patient_id"] = g.patient_id.iloc[0]
        if "tumor_stage" in g.columns:
            rec["tumor_stage"] = g.tumor_stage.iloc[0]
        for name in auc:
            rec[f"{name}_high"] = float(hi[f"AUC_{name}"].mean()) if len(hi) else np.nan
            rec[f"{name}_low"] = float(lo[f"AUC_{name}"].mean()) if len(lo) else np.nan
        for tf in focus_all:
            col = f"{tf}_log1p"
            if col in g:
                rec[f"{tf}_RNA_high"] = float(hi[col].mean()) if len(hi) else np.nan
                rec[f"{tf}_RNA_low"] = float(lo[col].mean()) if len(lo) else np.nan
        samp_rows.append(rec)
    samp = pd.DataFrame(samp_rows)
    samp.to_csv(RES / "sample_means_high_vs_low.tsv", sep="\t", index=False)

    paired = samp.loc[(samp.n_high >= MIN_CELLS) & (samp.n_low >= MIN_CELLS)].copy()
    print(f"[paired] samples with >={MIN_CELLS} high and low: {len(paired)}", flush=True)
    for name in auc:
        w = wilcoxon_paired(paired[f"{name}_high"], paired[f"{name}_low"])
        paired_rows.append(dict(
            regulon=name, kind="public_prior_AUCell",
            n_samples=w["n"], delta_median_high_minus_low=w["delta"],
            wilcoxon_stat=w["stat"], p=w["p"],
            note="paired sample means; public prior AUCell; not inferred on these cells",
        ))
    for tf in focus_all:
        hcol, lcol = f"{tf}_RNA_high", f"{tf}_RNA_low"
        if hcol in paired.columns:
            w = wilcoxon_paired(paired[hcol], paired[lcol])
            paired_rows.append(dict(
                regulon=f"{tf}_RNA", kind="tf_rna",
                n_samples=w["n"], delta_median_high_minus_low=w["delta"],
                wilcoxon_stat=w["stat"], p=w["p"],
                note="paired sample means of TF log1p; ELF3 RNA vs CLDN4 is A10-given, not a discovery",
            ))
    paired_df = pd.DataFrame(paired_rows)
    if len(paired_df) and paired_df.p.notna().any():
        mask = paired_df.p.notna()
        paired_df.loc[mask, "fdr"] = multipletests(paired_df.loc[mask, "p"], method="fdr_bh")[1]
    paired_df.to_csv(RES / "paired_high_vs_low.tsv", sep="\t", index=False)
    print("\n=== paired prior AUCell / TF RNA (sample) ===", flush=True)
    if len(paired_df):
        print(paired_df[["regulon", "n_samples", "delta_median_high_minus_low", "p"]].to_string(index=False),
              flush=True)

    tlung = (mal.Sample_Origin == "tLung") & primary
    tlung_high = high & (mal.Sample_Origin == "tLung")
    print(f"[tLung] malignant primary={int(tlung.sum())} CLDN4-high={int(tlung_high.sum())}",
          flush=True)
    if int(tlung_high.sum()) >= 50:
        pr_tlung = pearson_block(logX[tlung_high], genes, focus_all)
        pr_tlung.to_csv(RES / "pearson_cldn4_high_tLung.tsv.gz", sep="\t", index=False)
        tlung_rows = []
        for tf in focus_all:
            top = top_regulon(pr_tlung, tf)
            tlung_rows.append(dict(
                regulon=f"{tf}_pearson_cldn4_high_tLung", tf=tf,
                kind="pearson_cldn4_high_tLung",
                n_targets=int(len(top)),
                mean_r=float(top.pearson_r.mean()) if len(top) else np.nan,
                targets=";".join(top.target.tolist()),
                n_cells=int(tlung_high.sum()),
                n_samples=int(mal.loc[tlung_high, "sample"].nunique()),
                note="tLung-only sensitivity; CLDN4-high within-sample; CLDN4 held out",
            ))
        pd.DataFrame(tlung_rows).to_csv(RES / "regulons_tLung.tsv", sep="\t", index=False)

    origin_counts = (
        mal.loc[mal.in_primary_samples, "Sample_Origin"].value_counts(dropna=False).to_dict()
    )
    n_tab = pd.DataFrame([
        dict(item="cells_in_matrix", n=n_raw, note="GSE131907 public UMI header"),
        dict(item="author_malignant", n=int(mal_mask.sum()),
             note="Cell_subtype in {Malignant cells, tS1, tS2, tS3}"),
        dict(item="malignant_UMI>=200", n=int(qc.sum()), note=""),
        dict(item="samples_ge20_malignant", n=len(keep_samples),
             note=f"min {MIN_CELLS} author-malignant after UMI QC; primary GRN samples"),
        dict(item="CLDN4_high_within_sample", n=n_high,
             note=">= sample median CLDN4 log1p among author-malignant"),
        dict(item="CLDN4_low_within_sample", n=n_low, note="complement in the same samples"),
        dict(item="samples_paired_ge20_high_and_low", n=int(len(paired)),
             note="unit of prior-AUCell high vs low"),
        dict(item="tLung_samples_ge20", n=int(mal.loc[tlung, "sample"].nunique()),
             note="primary-tumor sensitivity"),
        dict(item="tLung_CLDN4_high", n=int(tlung_high.sum()), note="cells"),
        dict(item="TFs_screened_in_high", n=len(tf_present),
             note=f"human TF list ∩ det>={TF_MIN_DET} in CLDN4-high"),
    ])
    n_tab.to_csv(RES / "n_table.tsv", sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(8.2, 5.4))
    if len(tf_stats):
        foc = tf_stats.loc[tf_stats.tf.isin(FOCUS)].copy()
        extra = tf_stats.loc[~tf_stats.tf.isin(FOCUS + CONTROLS)].head(10)
        plot = pd.concat([foc, extra], ignore_index=True)
        plot = plot.sort_values("n_high_specific")
        colors = ["#c45c26" if t in GIVEN_TFS else ("#1f4e79" if t in FOCUS else "#7a7a7a")
                  for t in plot.tf]
        ax.barh(plot.tf, plot.n_high_specific, color=colors)
    ax.set_xlabel("n high-specific targets (r_high≥0.15 and r_low<0.10)")
    ax.set_title("GSE131907 CLDN4-high author-malignant — Pearson GRN proxy\n"
                 f"n_high={n_high} cells / {len(keep_samples)} samples; "
                 "orange=ELF3 given; navy=focus; grey=data-driven")
    fig.tight_layout()
    fig.savefig(RES / "fig_high_specific_counts.png", dpi=160)
    fig.savefig(RES / "fig_high_specific_counts.pdf")
    plt.close(fig)

    if len(paired):
        fig, ax = plt.subplots(figsize=(7.2, 4.6))
        regs = [f"{t}_prior" for t in FOCUS if f"{t}_prior_high" in paired.columns]
        ys, labs = [], []
        for i, name in enumerate(regs):
            h = paired[f"{name}_high"]
            l = paired[f"{name}_low"]
            w = wilcoxon_paired(h, l)
            ax.scatter(w["delta"], i, s=50, color="#1f4e79")
            ax.plot([0, w["delta"] if np.isfinite(w["delta"]) else 0], [i, i],
                    color="#1f4e79", lw=1.4)
            labs.append(f"{name} p={w['p']:.3g}" if np.isfinite(w["p"]) else name)
            ys.append(i)
        ax.axvline(0, color="0.5", lw=0.8)
        ax.set_yticks(ys)
        ax.set_yticklabels(labs, fontsize=8)
        ax.set_xlabel("median Δ AUCell (CLDN4-high − low), sample-paired")
        ax.set_title(f"Public-prior AUCell · n={len(paired)} samples with ≥{MIN_CELLS} high and low")
        fig.tight_layout()
        fig.savefig(RES / "fig_prior_aucell_paired.png", dpi=160)
        fig.savefig(RES / "fig_prior_aucell_paired.pdf")
        plt.close(fig)

    primary_regs = regulons.loc[
        regulons.kind.isin(["pearson_cldn4_high", "public_prior", "high_specific"])
        & regulons.tf.isin(FOCUS)
    ].copy()
    primary_regs.to_csv(RES / "regulons_focus.tsv", sep="\t", index=False)

    high_by_samp = (
        mal.loc[high].groupby("sample").size().to_dict() if high.any() else {}
    )
    summary = {
        "dataset": "GSE131907",
        "pmid": "32385277",
        "question": "Which TF regulons are recovered from malignant CLDN4-high vs low cells?",
        "a10_taken_as_given": "ELF3–CLDN4 bulk RNA; ELF3 is labeled given, not a discovery.",
        "method": "lightweight GRN proxy (Pearson + public priors + AUCell); not pySCENIC cisTarget; not ChIP",
        "pyscenic": {"importable": pyscenic_ok, "version": pyscenic_ver,
                     "cistarget_run": False, "chip_peaks_invented": False},
        "n_genes_in_matrix": n_genes_full,
        "n_cells_raw": n_raw,
        "n_malignant_author": int(mal_mask.sum()),
        "malignant_definition": "author Cell_subtype in {Malignant cells, tS1, tS2, tS3}",
        "n_cldn4_high": n_high,
        "n_cldn4_low": n_low,
        "n_samples_primary": len(keep_samples),
        "n_samples_paired": int(len(paired)),
        "cldn4_split": "within-sample median among author-malignant",
        "global_cldn4_median": global_med,
        "n_genes_scored": len(genes),
        "n_tfs_screened": len(tf_present),
        "origin_primary_cells": {str(k): int(v) for k, v in origin_counts.items()},
        "high_cells_by_sample": {str(k): int(v) for k, v in high_by_samp.items()},
        "primary_output": "CLDN4-high regulon table only",
        "primary_unit": "sample (paired high vs low for prior AUCell)",
        "cell_level": "exploratory only",
        "notes": [
            "Author malignant labels are on GEO (unlike GSE207422 CopyKAT).",
            "CLDN4 held out of every regulon (split gene).",
            "CLDN4-low regulons are computed only to mark high-specific edges and are not reported as a table.",
            "Pearson neighborhoods are co-expression, not binding.",
            "EBUS_28 and NS_07 can dominate pooled cell counts; sample n is the claim.",
            "log2TPM text (2.86 GB) and EGA FASTQ were not used.",
            "No ICI / MPR labels in this treatment-naive atlas.",
        ],
    }
    (RES / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")
    write_finding(regulons, n_tab, paired_df, summary)
    print(f"[done] {RES}", flush=True)


if __name__ == "__main__":
    sys.exit(main() or 0)
