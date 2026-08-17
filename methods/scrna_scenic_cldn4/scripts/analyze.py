#!/usr/bin/env python3
"""CLDN4-high malignant GRN proxy on public GSE207422.

Additive to methods/scrna_scenic (ELF3/GRHL1/KLF4 AUCell vs TACSTD2/CLDN4).
This script infers regulons in CLDN4-high malignant-like cells only.

pySCENIC cisTarget is not run. No ChIP peaks are invented. Regulons are:

  1. public curated TF–target priors (TRRUST / DoRothEA / CollecTRI)
  2. Pearson co-expression in CLDN4-high cells (top-N)
  3. high-specific edges: r_high >= 0.15 and r_low < 0.10

A10 ELF3–CLDN4 bulk RNA is taken as given and is not re-tested as a discovery.
CLDN4 is held out of every regulon (it is the split gene).
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
RES = REPO / "results" / "scrna_scenic_cldn4"
RES.mkdir(parents=True, exist_ok=True)

GEO_CANDIDATES = [
    Path("/tmp/scrna_scenic_cldn4/geo"),
    Path("/tmp/scrna_scenic/geo"),
    REPO / "data" / "scrna_scenic_cldn4" / "GSE207422",
    REPO / "data" / "scrna_scenic" / "GSE207422",
    Path("data/scrna_scenic_cldn4/GSE207422"),
    Path("data/scrna_scenic/GSE207422"),
]
PRIORS = ROOT / "resources" / "tf_targets_public.tsv"
TF_LIST = ROOT / "resources" / "human_tfs_from_priors.txt"

FOCUS = ["ELF3", "GRHL1", "GRHL2", "KLF4", "KLF5",
         "OVOL1", "OVOL2", "TFAP2A", "SP1", "TP63"]
CONTROLS = ["NKX2-1", "SOX2", "FOXA1"]
GIVEN_TFS = {"ELF3"}  # A10 ELF3–CLDN4 taken as given
HOLD_OUT = {"CLDN4"}

LINEAGE = {
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT5", "KRT7", "KRT17",
                   "ELF3", "CDH1", "MUC1"],
    "T_NK": ["CD3D", "CD3E", "CD3G", "TRAC", "CD2", "NKG7", "GNLY", "KLRD1"],
    "B_Plasma": ["CD79A", "CD79B", "MS4A1", "JCHAIN", "MZB1", "IGHM"],
    "Myeloid": ["LYZ", "CD68", "CD14", "C1QA", "C1QB", "FCN1", "ITGAX"],
    "Mast": ["TPSAB1", "TPSB2", "CPA3", "MS4A2"],
    "Endothelial": ["PECAM1", "VWF", "CLDN5", "CDH5", "RAMP2"],
    "Fibroblast": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "TAGLN"],
}
NORMAL_LUNG = ["SFTPA1", "SFTPA2", "SFTPB", "SFTPC", "SFTPD", "AGER", "NAPSA",
               "SCGB1A1", "SCGB3A2", "TPPP3", "FOXJ1", "CAPS"]

MIN_UMI = 200
NL_Q = 0.75
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
        if (p / "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz").exists():
            return p
    raise SystemExit("GSE207422 matrix not found. Run scripts/download_gse207422.py")


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


def stream_pass1(matrix: Path, keep_genes: set[str]):
    with gzip.open(matrix, "rb") as fh:
        header = fh.readline().decode("ascii").rstrip("\n")
        cells = header.split("\t")[1:]
        n = len(cells)
        print(f"[pass1] cells={n}", flush=True)
        total = np.zeros(n, dtype=np.float64)
        kept: dict[str, np.ndarray] = {}
        n_genes = 0
        for raw in fh:
            n_genes += 1
            tab = raw.find(b"\t")
            gene = raw[:tab].decode("ascii")
            arr = np.fromstring(raw[tab + 1 :].rstrip(b"\r\n"), sep="\t", dtype=np.float32)
            if arr.size != n:
                raise ValueError(f"{gene}: {arr.size} != {n}")
            total += arr
            if gene in keep_genes:
                kept[gene] = arr.copy()
            if n_genes % 2000 == 0:
                print(f"  pass1 genes={n_genes} kept={len(kept)}", flush=True)
    print(f"[pass1] done genes={n_genes} kept={len(kept)}", flush=True)
    return cells, total, kept, n_genes


def stream_pass2(matrix: Path, keep_idx: np.ndarray, always: set[str], tfs: set[str]):
    keep_idx = np.asarray(keep_idx, dtype=int)
    n_keep = int(keep_idx.size)
    with gzip.open(matrix, "rb") as fh:
        header = fh.readline().decode("ascii").rstrip("\n")
        n = len(header.split("\t")) - 1
        genes: list[str] = []
        mats: list[np.ndarray] = []
        means: list[float] = []
        dets: list[float] = []
        n_genes = 0
        for raw in fh:
            n_genes += 1
            tab = raw.find(b"\t")
            gene = raw[:tab].decode("ascii")
            arr = np.fromstring(raw[tab + 1 :].rstrip(b"\r\n"), sep="\t", dtype=np.float32)
            if arr.size != n:
                raise ValueError(f"{gene}: {arr.size} != {n}")
            sub = arr[keep_idx]
            det = float((sub > 0).mean())
            if gene in always or gene in tfs or det >= MIN_DET:
                genes.append(gene)
                mats.append(sub.copy())
                means.append(float(sub.mean()))
                dets.append(det)
            if n_genes % 2000 == 0:
                print(f"  pass2 genes={n_genes} stored={len(genes)}", flush=True)
    print(f"[pass2] stored {len(genes)} / {n_genes} genes for {n_keep} cells", flush=True)
    X = np.vstack(mats).T.astype(np.float32)
    info = pd.DataFrame({"gene": genes, "mean": means, "detection": dets})
    always_mask = info.gene.isin(always | tfs).to_numpy()
    if (~always_mask).sum() > HVG_N:
        var = X[:, ~always_mask].astype(np.float64).var(axis=0)
        keep_rel = np.argsort(-var)[:HVG_N]
        rel_idx = np.where(~always_mask)[0][keep_rel]
        keep_g = np.sort(np.unique(np.concatenate([np.where(always_mask)[0], rel_idx])))
        X = X[:, keep_g]
        info = info.iloc[keep_g].reset_index(drop=True)
        print(f"[pass2] HVG filter -> {X.shape[1]} genes", flush=True)
    return X, info


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
    """Vectorized Pearson of each TF column vs all genes. Returns long table."""
    gi = {g: i for i, g in enumerate(genes)}
    present = [t for t in tfs if t in gi]
    if not present:
        return pd.DataFrame(columns=["tf", "target", "pearson_r"])
    X = np.asarray(logX, float)
    # drop cells with non-finite library
    ok = np.isfinite(X).all(axis=1)
    X = X[ok]
    if X.shape[0] < 20:
        return pd.DataFrame(columns=["tf", "target", "pearson_r"])
    X = X - X.mean(axis=0, keepdims=True)
    sd = X.std(axis=0, ddof=1)
    sd[sd == 0] = np.nan
    Z = X / sd
    tf_idx = np.array([gi[t] for t in present], dtype=int)
    # (n_tf x n_cells) @ (n_cells x n_genes) / (n-1)
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


def main() -> None:
    geo = find_geo()
    matrix = geo / "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"
    meta_p = geo / "GSE207422_NSCLC_scRNAseq_metadata.xlsx"
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

    keep = set(focus_all) | {"CLDN4", "TACSTD2", "EPCAM"}
    for gs in LINEAGE.values():
        keep.update(gs)
    keep.update(NORMAL_LUNG)
    keep.update(priors.target.tolist())
    keep.update(priors.tf.tolist())

    cells, total, kept, n_genes_full = stream_pass1(matrix, keep)
    n = len(cells)
    qc = total >= MIN_UMI
    print(f"[qc] UMI>={MIN_UMI}: {int(qc.sum())}/{n}", flush=True)

    marker_genes = sorted(g for g in kept if g in keep)
    M = np.column_stack([kept[g] for g in marker_genes])
    ln_all = log1p_cp10k(M, total)
    gix = {g: i for i, g in enumerate(marker_genes)}

    def score(names):
        idx = [gix[g] for g in names if g in gix]
        if not idx:
            return np.zeros(n)
        return ln_all[:, idx].mean(axis=1)

    scores = {lin: score(gs) for lin, gs in LINEAGE.items()}
    S = np.column_stack([scores[k] for k in LINEAGE])
    lineage = np.array(list(LINEAGE.keys()))[S.argmax(axis=1)]
    nl = score(NORMAL_LUNG)
    epi = (lineage == "Epithelial") & qc
    nl_cut = float(np.quantile(nl[epi], NL_Q)) if epi.any() else 0.0
    mal = epi & (nl <= nl_cut)
    print(f"[gate] epithelial={int(epi.sum())} malignant_like={int(mal.sum())} nl_cut={nl_cut:.4f}",
          flush=True)

    if "CLDN4" not in gix:
        raise SystemExit("CLDN4 absent from public UMI matrix")
    cldn4 = ln_all[:, gix["CLDN4"]]

    obs = pd.DataFrame({
        "barcode": cells,
        "sample": pd.Series(cells).str.extract(r"(BD_immune\d+)", expand=False),
        "total_umi": total.astype(np.int64),
        "lineage": lineage,
        "normal_lung": nl,
        "epithelial": epi,
        "malignant_like": mal,
        "CLDN4_log1p": cldn4,
    })
    if "TACSTD2" in gix:
        obs["TACSTD2_log1p"] = ln_all[:, gix["TACSTD2"]]
    for tf in focus_all:
        if tf in gix:
            obs[f"{tf}_log1p"] = ln_all[:, gix[tf]]

    meta = pd.read_excel(meta_p)
    meta = meta[meta["Sample"].astype(str).str.startswith("BD_immune")].copy()
    meta = meta.rename(columns={"Sample": "sample"})
    meta["histology"] = meta["Pathology"].map({"Adeno": "LUAD", "Squamous": "LUSC"})
    meta["timing"] = meta["Resource"].map({
        "Pre-treatment biopsy": "pre",
        "Post-treatment surgery": "post",
    })
    obs = obs.merge(
        meta[["sample", "Patient", "histology", "timing", "Pathologic Response", "Pathology"]],
        on="sample", how="left",
    )

    # Primary split: within-sample median among malignant-like, samples with >=20
    mal_obs = obs.loc[obs.malignant_like].copy()
    samp_n = mal_obs.groupby("sample").size()
    keep_samples = set(samp_n[samp_n >= MIN_CELLS].index)
    mal_obs = mal_obs.loc[mal_obs.sample.isin(keep_samples)].copy()
    med_within = mal_obs.groupby("sample")["CLDN4_log1p"].transform("median")
    mal_obs["cldn4_high_within"] = mal_obs["CLDN4_log1p"] >= med_within
    global_med = float(obs.loc[obs.malignant_like, "CLDN4_log1p"].median())
    mal_obs["cldn4_high_global"] = mal_obs["CLDN4_log1p"] >= global_med
    # primary flag on full obs
    obs["in_primary_samples"] = obs.sample.isin(keep_samples) & obs.malignant_like
    obs = obs.merge(
        mal_obs[["barcode", "cldn4_high_within", "cldn4_high_global"]],
        on="barcode", how="left",
    )
    obs["cldn4_high_within"] = obs["cldn4_high_within"].fillna(False).astype(bool)
    obs["cldn4_high_global"] = obs["cldn4_high_global"].fillna(False).astype(bool)

    n_high = int(obs.cldn4_high_within.sum())
    n_low = int((obs.in_primary_samples & ~obs.cldn4_high_within).sum())
    print(f"[split] within-sample CLDN4-high={n_high} low={n_low} "
          f"samples={len(keep_samples)} global_med={global_med:.3f}", flush=True)

    always = set(focus_all) | {"CLDN4", "TACSTD2"} | set(priors.target) | set(NORMAL_LUNG)
    for gs in LINEAGE.values():
        always.update(gs)

    mal_idx = np.where(obs.malignant_like.to_numpy())[0]
    if mal_idx.size < 50:
        raise SystemExit(f"too few malignant-like cells: {mal_idx.size}")

    X_counts, gene_info = stream_pass2(matrix, mal_idx, always, human_tfs)
    genes = gene_info.gene.tolist()
    gene_set = set(genes)
    mal_full = obs.iloc[mal_idx].reset_index(drop=True)
    logX = log1p_cp10k(X_counts, mal_full.total_umi.to_numpy())
    print(f"[expr] malignant-like matrix {logX.shape}", flush=True)

    # restrict GRN to primary-sample cells
    primary = mal_full.in_primary_samples.to_numpy()
    high = mal_full.cldn4_high_within.to_numpy() & primary
    low = (~mal_full.cldn4_high_within.to_numpy()) & primary
    print(f"[grn cells] high={int(high.sum())} low={int(low.sum())}", flush=True)

    # TFs detected in CLDN4-high
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
    # keep a compact low table for the high-specific join only
    pr_low_s = pr_low.rename(columns={"pearson_r": "pearson_r_low"})
    joined = pr_high.merge(pr_low_s, on=["tf", "target"], how="left")
    joined["high_specific"] = (
        (joined.pearson_r >= HIGH_SPEC_R)
        & (joined.pearson_r_low.fillna(0) < LOW_SPEC_R)
        & (~joined.target.isin(HOLD_OUT))
    )
    joined.to_csv(RES / "pearson_high_vs_low.tsv.gz", sep="\t", index=False)

    # Build CLDN4-high regulon table (high only)
    rows = []
    # public priors for focus TFs
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

    # Pearson regulons for focus + top data-driven TFs
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

    screen_tfs = tf_stats.head(TOP_SCREEN).tf.tolist()
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

    # AUCell of public priors on all malignant-like (not circular: priors are external)
    prior_regs = regulons.loc[regulons.kind == "public_prior"]
    auc = {}
    for rec in prior_regs.itertuples(index=False):
        members = [g for g in (rec.targets.split(";") if rec.targets else []) if g in gene_set]
        mask = np.array([g in set(members) for g in genes])
        auc[rec.regulon] = aucell(logX, mask) if mask.any() else np.full(logX.shape[0], np.nan)
        print(f"[aucell] {rec.regulon} n={int(mask.sum())} mean={np.nanmean(auc[rec.regulon]):.4f}",
              flush=True)
        mal_full[f"AUC_{rec.regulon}"] = auc[rec.regulon]

    # sample-level paired high vs low (primary unit)
    paired_rows = []
    samp_rows = []
    use = mal_full.loc[mal_full.in_primary_samples].copy()
    for sample, g in use.groupby("sample"):
        hi = g.loc[g.cldn4_high_within]
        lo = g.loc[~g.cldn4_high_within]
        rec = dict(
            sample=sample,
            n_malignant=len(g),
            n_high=len(hi),
            n_low=len(lo),
            histology=g.histology.iloc[0],
            Patient=g.Patient.iloc[0],
            timing=g.timing.iloc[0],
            CLDN4_high=float(hi.CLDN4_log1p.mean()) if len(hi) else np.nan,
            CLDN4_low=float(lo.CLDN4_log1p.mean()) if len(lo) else np.nan,
        )
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

    # paired tests: samples with >=20 high AND >=20 low
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

    # LUAD-only sensitivity Pearson for focus TFs
    luad = (mal_full.histology == "LUAD") & primary
    luad_high = high & (mal_full.histology == "LUAD")
    print(f"[luad] malignant-like primary={int(luad.sum())} CLDN4-high={int(luad_high.sum())}",
          flush=True)
    if int(luad_high.sum()) >= 50:
        pr_luad = pearson_block(logX[luad_high], genes, focus_all)
        pr_luad.to_csv(RES / "pearson_cldn4_high_LUAD.tsv.gz", sep="\t", index=False)
        luad_rows = []
        for tf in focus_all:
            top = top_regulon(pr_luad, tf)
            luad_rows.append(dict(
                regulon=f"{tf}_pearson_cldn4_high_LUAD", tf=tf,
                kind="pearson_cldn4_high_LUAD",
                n_targets=int(len(top)),
                mean_r=float(top.pearson_r.mean()) if len(top) else np.nan,
                targets=";".join(top.target.tolist()),
                n_cells=int(luad_high.sum()),
                n_samples=int(mal_full.loc[luad_high, "sample"].nunique()),
                note="LUAD-only sensitivity; CLDN4-high within-sample; CLDN4 held out",
            ))
        pd.DataFrame(luad_rows).to_csv(RES / "regulons_LUAD.tsv", sep="\t", index=False)

    # n table
    n_tab = pd.DataFrame([
        dict(item="cells_raw", n=n, note="GSE207422 public UMI"),
        dict(item="cells_UMI>=200", n=int(qc.sum()), note=""),
        dict(item="epithelial", n=int(epi.sum()), note="marker argmax"),
        dict(item="malignant_like", n=int(mal.sum()),
             note=f"epithelial AND normal-lung <= p{int(NL_Q*100)}; not CopyKAT"),
        dict(item="samples_ge20_malignant", n=len(keep_samples),
             note=f"min {MIN_CELLS} malignant-like; primary GRN samples"),
        dict(item="CLDN4_high_within_sample", n=n_high,
             note=">= sample median CLDN4 log1p among malignant-like"),
        dict(item="CLDN4_low_within_sample", n=n_low, note="complement in the same samples"),
        dict(item="samples_paired_ge20_high_and_low", n=int(len(paired)),
             note="unit of prior-AUCell high vs low"),
        dict(item="LUAD_samples_ge20", n=int(mal_full.loc[luad, "sample"].nunique()),
             note="sensitivity"),
        dict(item="LUAD_CLDN4_high", n=int(luad_high.sum()), note="cells"),
        dict(item="TFs_screened_in_high", n=len(tf_present),
             note=f"human TF list ∩ det>={TF_MIN_DET} in CLDN4-high"),
    ])
    n_tab.to_csv(RES / "n_table.tsv", sep="\t", index=False)

    # figures
    fig, ax = plt.subplots(figsize=(8.2, 5.4))
    foc = tf_stats.loc[tf_stats.tf.isin(FOCUS)].copy()
    extra = tf_stats.loc[~tf_stats.tf.isin(FOCUS + CONTROLS)].head(10)
    plot = pd.concat([foc, extra], ignore_index=True)
    plot = plot.sort_values("n_high_specific")
    colors = ["#c45c26" if t in GIVEN_TFS else ("#1f4e79" if t in FOCUS else "#7a7a7a")
              for t in plot.tf]
    ax.barh(plot.tf, plot.n_high_specific, color=colors)
    ax.set_xlabel("n high-specific targets (r_high≥0.15 and r_low<0.10)")
    ax.set_title("GSE207422 CLDN4-high malignant-like — Pearson GRN proxy\n"
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

    # compact primary regulon markdown-friendly TSV
    primary = regulons.loc[
        regulons.kind.isin(["pearson_cldn4_high", "public_prior", "high_specific"])
        & regulons.tf.isin(FOCUS)
    ].copy()
    primary.to_csv(RES / "regulons_focus.tsv", sep="\t", index=False)

    hist_counts = mal_full.loc[mal_full.in_primary_samples, "histology"].value_counts(dropna=False).to_dict()
    high_by_samp = (
        mal_full.loc[high].groupby("sample").size().to_dict()
        if high.any() else {}
    )
    summary = {
        "dataset": "GSE207422",
        "pmid": "36869384",
        "question": "Which TF regulons are recovered from CLDN4-high malignant-like cells?",
        "a10_taken_as_given": "ELF3–CLDN4 bulk RNA; ELF3 is labeled given, not a discovery.",
        "method": "lightweight GRN proxy (Pearson + public priors + AUCell); not pySCENIC cisTarget; not ChIP",
        "pyscenic": {"importable": pyscenic_ok, "version": pyscenic_ver,
                     "cistarget_run": False, "chip_peaks_invented": False},
        "n_genes_in_matrix": n_genes_full,
        "n_cells_raw": n,
        "n_epithelial": int(epi.sum()),
        "n_malignant_like": int(mal.sum()),
        "nl_cut": nl_cut,
        "n_cldn4_high": n_high,
        "n_cldn4_low": n_low,
        "n_samples_primary": len(keep_samples),
        "n_samples_paired": int(len(paired)),
        "cldn4_split": "within-sample median among malignant-like",
        "global_cldn4_median": global_med,
        "n_genes_scored": len(genes),
        "n_tfs_screened": len(tf_present),
        "histology_primary_cells": hist_counts,
        "high_cells_by_sample": {str(k): int(v) for k, v in high_by_samp.items()},
        "primary_output": "CLDN4-high regulon table only",
        "primary_unit": "sample (paired high vs low for prior AUCell)",
        "cell_level": "exploratory only",
        "notes": [
            "Author CopyKAT malignant labels are not on GEO; malignant-like is a marker proxy.",
            "CLDN4 held out of every regulon (split gene).",
            "CLDN4-low regulons are computed only to mark high-specific edges and are not reported as a table.",
            "Pearson neighborhoods are co-expression, not binding.",
            "BD_immune07 (LUSC) can dominate pooled cell counts; sample n is the claim.",
            "GSE253013 and GSE131907 were not required.",
        ],
    }
    (RES / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")
    print(f"[done] {RES}", flush=True)


if __name__ == "__main__":
    sys.exit(main() or 0)
