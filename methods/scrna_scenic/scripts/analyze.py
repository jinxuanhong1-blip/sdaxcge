#!/usr/bin/env python3
"""ELF3/GRHL1/KLF4 regulons vs TACSTD2/CLDN4 on public GSE207422 epithelium.

pySCENIC is importable here, but cisTarget motif databases are not fetched and
no ChIP peaks are invented. Regulons are:

  1. public curated TF–target priors (TRRUST / DoRothEA / CollecTRI)
  2. Pearson co-expression in this matrix (top-N, TACSTD2/CLDN4 held out)
  3. optional GRNBoost2 (arboreto) for the three TFs only

AUCell is implemented from the Aibar 2017 recovery-curve formula.
Cell-level Spearman is exploratory (pseudoreplication). Primary n is samples
with >=20 malignant-like cells. LUAD-only = Pathology==Adeno; all = Adeno+Squamous.

A10 ELF3–CLDN4 bulk RNA is taken as given and is not re-tested.
"""
from __future__ import annotations

import gzip
import json
import sys
import traceback
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
RES = REPO / "results" / "scrna_scenic"
RES.mkdir(parents=True, exist_ok=True)

# Prefer a local GEO cache if the download script was pointed at /tmp.
GEO_CANDIDATES = [
    Path("/tmp/scrna_scenic/geo"),
    REPO / "data" / "scrna_scenic" / "GSE207422",
    Path("data/scrna_scenic/GSE207422"),
]
PRIORS = ROOT / "resources" / "tf_targets_public.tsv"

FOCUS = ["ELF3", "GRHL1", "KLF4"]
CONTROLS = ["TFAP2A", "NKX2-1"]
ALL_TFS = FOCUS + CONTROLS
TARGETS = ["TACSTD2", "CLDN4"]
HOLD_OUT = set(TARGETS)

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
PEARSON_TOP = 50
PEARSON_MIN_R = 0.10
GRNBOOST_TOP = 50
AUC_THR = 0.05


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
    """Aibar-style AUCell: normalized recovery of set members in the top ranks.

    For each cell, genes are ranked by expression (0 = highest). Members with
    rank < max_rank contribute (max_rank - rank). Score is divided by
    n_set * max_rank. This is the common linear-recovery form; it is not the
    R AUCell package binary, and it is not a binding assay.
    """
    n_cells, n_genes = expr.shape
    max_rank = max(int(np.ceil(auc_threshold * n_genes)), 1)
    idx = np.where(member)[0]
    n_set = int(idx.size)
    if n_set == 0:
        return np.full(n_cells, np.nan)
    # ranks: 0 = most expressed
    order = np.argsort(-expr, axis=1, kind="mergesort")
    ranks = np.empty_like(order)
    row = np.arange(n_cells)[:, None]
    ranks[row, order] = np.arange(n_genes)[None, :]
    r = ranks[:, idx].astype(np.float64)
    contrib = np.clip(max_rank - r, 0, None)
    return contrib.sum(axis=1) / (n_set * max_rank)


def stream_pass1(matrix: Path, keep_genes: set[str]):
    """Return barcodes, total UMI, and a dense cells x keep-genes count frame."""
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


def stream_pass2(matrix: Path, epi_idx: np.ndarray, always: set[str]):
    """Extract epithelial columns; keep HVGs + always-keep genes."""
    epi_idx = np.asarray(epi_idx, dtype=int)
    n_epi = int(epi_idx.size)
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
            sub = arr[epi_idx]
            det = float((sub > 0).mean())
            if gene in always or det >= MIN_DET:
                genes.append(gene)
                mats.append(sub.copy())
                means.append(float(sub.mean()))
                dets.append(det)
            if n_genes % 2000 == 0:
                print(f"  pass2 genes={n_genes} stored={len(genes)}", flush=True)
    print(f"[pass2] stored {len(genes)} / {n_genes} genes for {n_epi} epithelial cells", flush=True)
    X = np.vstack(mats).T.astype(np.float32)  # cells x genes
    info = pd.DataFrame({"gene": genes, "mean": means, "detection": dets})
    # HVG among non-always genes
    always_mask = info.gene.isin(always).to_numpy()
    if (~always_mask).sum() > HVG_N:
        # variance on log1p of raw counts (library-size later)
        var = X[:, ~always_mask].astype(np.float64).var(axis=0)
        keep_rel = np.argsort(-var)[:HVG_N]
        rel_idx = np.where(~always_mask)[0][keep_rel]
        keep_idx = np.sort(np.unique(np.concatenate([np.where(always_mask)[0], rel_idx])))
        X = X[:, keep_idx]
        info = info.iloc[keep_idx].reset_index(drop=True)
        print(f"[pass2] HVG filter -> {X.shape[1]} genes", flush=True)
    return X, info


def load_priors() -> pd.DataFrame:
    df = pd.read_csv(PRIORS, sep="\t")
    df["tf"] = df["tf"].astype(str)
    df["target"] = df["target"].astype(str)
    return df


def prior_targets(priors: pd.DataFrame, tf: str, genes: set[str]) -> list[str]:
    hits = priors.loc[priors.tf == tf, "target"].unique().tolist()
    return sorted(g for g in hits if g in genes and g != tf and g not in HOLD_OUT)


def pearson_targets(logX: np.ndarray, genes: list[str], tf: str) -> pd.DataFrame:
    if tf not in genes:
        return pd.DataFrame(columns=["tf", "target", "pearson_r", "p", "fdr", "rank"])
    gi = {g: i for i, g in enumerate(genes)}
    x = logX[:, gi[tf]]
    rows = []
    for g, j in gi.items():
        if g == tf:
            continue
        y = logX[:, j]
        m = np.isfinite(x) & np.isfinite(y)
        if m.sum() < 20 or np.nanstd(x[m]) == 0 or np.nanstd(y[m]) == 0:
            continue
        r, p = stats.pearsonr(x[m], y[m])
        rows.append((tf, g, float(r), float(p)))
    out = pd.DataFrame(rows, columns=["tf", "target", "pearson_r", "p"])
    if out.empty:
        return out
    out["fdr"] = multipletests(out["p"], method="fdr_bh")[1]
    out = out.sort_values("pearson_r", ascending=False).reset_index(drop=True)
    out["rank"] = np.arange(1, len(out) + 1)
    return out


def try_grnboost(logX: np.ndarray, genes: list[str], tfs: list[str]) -> pd.DataFrame:
    import os

    present = [t for t in tfs if t in genes]
    empty = pd.DataFrame(columns=["tf", "target", "importance"])
    if not present:
        return empty
    if os.environ.get("SKIP_GRNBOOST", "").strip() in {"1", "true", "yes"}:
        print("[grnboost] skipped (SKIP_GRNBOOST)", flush=True)
        return empty
    try:
        from arboreto.algo import grnboost2
    except Exception as exc:
        print(f"[grnboost] arboreto not usable: {exc}", flush=True)
        return empty
    X = logX
    gnames = list(genes)
    if X.shape[0] > 3000:
        rng = np.random.default_rng(7)
        take = rng.choice(X.shape[0], 3000, replace=False)
        X = X[take]
        print(f"[grnboost] subsampled to {X.shape[0]} cells", flush=True)
    if X.shape[1] > 2500:
        # keep TFs + highest-variance genes
        tf_idx = [i for i, g in enumerate(gnames) if g in present]
        var = X.astype(np.float64).var(axis=0)
        keep = set(tf_idx) | set(np.argsort(-var)[:2500].tolist())
        keep_idx = sorted(keep)
        X = X[:, keep_idx]
        gnames = [gnames[i] for i in keep_idx]
        print(f"[grnboost] gene cap -> {X.shape[1]}", flush=True)
    expr = pd.DataFrame(X, columns=gnames)
    print(f"[grnboost] {expr.shape[0]} cells x {expr.shape[1]} genes; TFs={present}", flush=True)
    try:
        adj = grnboost2(expression_data=expr, tf_names=present, verbose=True, seed=7)
    except Exception:
        traceback.print_exc()
        return empty
    if adj is None or len(adj) == 0:
        return empty
    adj = adj.rename(columns={"TF": "tf", "target": "target", "importance": "importance"})
    adj = adj[adj.tf.isin(present)].copy()
    return adj.sort_values(["tf", "importance"], ascending=[True, False])


def build_regulon_table(priors, pearson, grn, genes: set[str]) -> pd.DataFrame:
    rows = []
    for tf in ALL_TFS:
        prior = prior_targets(priors, tf, genes)
        rows.append(dict(regulon=f"{tf}_prior", tf=tf, kind="public_prior",
                         n_targets=len(prior), targets=";".join(prior),
                         note="TRRUST+DoRothEA+CollecTRI union; TACSTD2/CLDN4 held out; not ChIP"))
        pr = pearson.get(tf)
        if pr is not None and len(pr):
            top = pr.loc[(pr.pearson_r >= PEARSON_MIN_R) & (~pr.target.isin(HOLD_OUT))].head(PEARSON_TOP)
            # if too few pass the r cut, still take top-N positive
            if len(top) < 10:
                top = pr.loc[(pr.pearson_r > 0) & (~pr.target.isin(HOLD_OUT))].head(PEARSON_TOP)
            rows.append(dict(regulon=f"{tf}_pearson", tf=tf, kind="pearson_coexpression",
                             n_targets=len(top), targets=";".join(top.target.tolist()),
                             note=f"top {len(top)} Pearson r>={PEARSON_MIN_R} (or top positive); held out TACSTD2/CLDN4"))
            inter = [g for g in prior if g in set(top.target)]
            rows.append(dict(regulon=f"{tf}_intersect", tf=tf, kind="prior_AND_pearson",
                             n_targets=len(inter), targets=";".join(inter),
                             note="intersection of public prior and Pearson top; not ChIP"))
        if grn is not None and len(grn):
            gtop = grn.loc[(grn.tf == tf) & (~grn.target.isin(HOLD_OUT))].head(GRNBOOST_TOP)
            if len(gtop):
                rows.append(dict(regulon=f"{tf}_grnboost", tf=tf, kind="grnboost2",
                                 n_targets=len(gtop), targets=";".join(gtop.target.tolist()),
                                 note="arboreto GRNBoost2 top importance; no cisTarget; not ChIP"))
    # combinatorial union of the three focus priors / pearson
    for kind, suffix in (("public_prior", "prior"), ("pearson_coexpression", "pearson")):
        members = []
        for tf in FOCUS:
            rec = next((r for r in rows if r["regulon"] == f"{tf}_{suffix}"), None)
            if rec:
                members.extend(rec["targets"].split(";") if rec["targets"] else [])
        members = sorted({g for g in members if g})
        rows.append(dict(regulon=f"ELF3_GRHL1_KLF4_{suffix}_union", tf="ELF3+GRHL1+KLF4",
                         kind=f"combinatorial_{kind}", n_targets=len(members),
                         targets=";".join(members),
                         note="union of the three focus regulons; TACSTD2/CLDN4 held out"))
    return pd.DataFrame(rows)


def verdict(sample_rows: pd.DataFrame) -> dict:
    """Apply pre-specified rules to sample-level primary regulons."""
    primary_regs = [f"{tf}_prior" for tf in FOCUS] + [f"{tf}_pearson" for tf in FOCUS]
    primary_regs += ["ELF3_GRHL1_KLF4_prior_union", "ELF3_GRHL1_KLF4_pearson_union"]
    sub = sample_rows.loc[
        (sample_rows.unit == "sample")
        & (sample_rows.compartment == "malignant_like")
        & (sample_rows.regulon.isin(primary_regs))
    ].copy()

    def meet(split, gene):
        hit = []
        for reg in primary_regs:
            r = sub[(sub.split == split) & (sub.regulon == reg) & (sub.gene == gene)]
            if r.empty:
                continue
            n = int(r.n.iloc[0])
            rho = float(r.rho.iloc[0])
            p = float(r.p.iloc[0])
            if n < MIN_SAMPLES:
                hit.append(("underpowered", reg, n, rho, p))
            elif np.isfinite(rho) and rho >= 0.30 and np.isfinite(p) and p < 0.05:
                hit.append(("support", reg, n, rho, p))
            else:
                hit.append(("no", reg, n, rho, p))
        return hit

    summary = {}
    for split in ("all", "LUAD"):
        t2 = meet(split, "TACSTD2")
        c4 = meet(split, "CLDN4")
        n_s = next((h[2] for h in t2), 0)
        if n_s < MIN_SAMPLES:
            call = "UNDERPOWERED"
        else:
            any_t2 = any(h[0] == "support" for h in t2)
            any_c4 = any(h[0] == "support" for h in c4)
            both_same = False
            for a in t2:
                for b in c4:
                    if a[1] == b[1] and a[0] == "support" and b[0] == "support":
                        both_same = True
            if both_same:
                call = "SUPPORT"
            elif any_t2 or any_c4:
                call = "PARTIAL"
            else:
                call = "NOT_SUPPORTED"
        summary[split] = {"call": call, "TACSTD2": t2, "CLDN4": c4, "n_samples": n_s}
    # overall
    calls = {summary["all"]["call"], summary["LUAD"]["call"]}
    if calls == {"SUPPORT"}:
        overall = "SUPPORT"
    elif "SUPPORT" in calls or "PARTIAL" in calls:
        overall = "PARTIAL"
    elif calls == {"UNDERPOWERED"}:
        overall = "UNDERPOWERED"
    else:
        overall = "NOT_SUPPORTED"
    return {"overall": overall, "by_split": {
        k: {"call": v["call"], "n_samples": v["n_samples"]} for k, v in summary.items()
    }, "detail": summary}


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
    print("[pyscenic] cisTarget motif DBs NOT downloaded; ctx/aucell CLI not run; no ChIP peaks.", flush=True)

    priors = load_priors()
    prior_genes = set(priors.target) | set(ALL_TFS) | set(TARGETS)
    keep = set(prior_genes)
    for gs in LINEAGE.values():
        keep.update(gs)
    keep.update(NORMAL_LUNG)
    keep.update(["EPCAM"])

    cells, total, kept, n_genes_full = stream_pass1(matrix, keep)
    n = len(cells)
    qc = total >= MIN_UMI
    print(f"[qc] UMI>={MIN_UMI}: {int(qc.sum())}/{n}", flush=True)

    # lineage on all cells (markers only)
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
    print(f"[gate] epithelial={int(epi.sum())} malignant_like={int(mal.sum())} nl_cut={nl_cut:.4f}", flush=True)

    obs = pd.DataFrame({
        "barcode": cells,
        "sample": pd.Series(cells).str.extract(r"(BD_immune\d+)", expand=False),
        "total_umi": total.astype(np.int64),
        "lineage": lineage,
        "normal_lung": nl,
        "epithelial": epi,
        "malignant_like": mal,
    })
    for g in ALL_TFS + TARGETS + ["EPCAM"]:
        obs[g] = kept[g] if g in kept else np.nan
        obs[f"{g}_log1p"] = ln_all[:, gix[g]] if g in gix else np.nan

    meta = pd.read_excel(meta_p)
    meta = meta[meta["Sample"].astype(str).str.startswith("BD_immune")].copy()
    meta = meta.rename(columns={"Sample": "sample"})
    meta["histology"] = meta["Pathology"].map({"Adeno": "LUAD", "Squamous": "LUSC"})
    meta["timing"] = meta["Resource"].map({
        "Pre-treatment biopsy": "pre",
        "Post-treatment surgery": "post",
    })
    obs = obs.merge(meta[["sample", "Patient", "histology", "timing",
                          "Pathologic Response", "Pathology"]], on="sample", how="left")

    always = set(ALL_TFS) | set(TARGETS) | set(priors.target) | set(NORMAL_LUNG)
    for gs in LINEAGE.values():
        always.update(gs)

    epi_idx = np.where(obs.epithelial.to_numpy())[0]
    if epi_idx.size < 50:
        raise SystemExit(f"too few epithelial cells: {epi_idx.size}")

    X_counts, gene_info = stream_pass2(matrix, epi_idx, always)
    genes = gene_info.gene.tolist()
    gene_set = set(genes)
    epi_obs = obs.iloc[epi_idx].reset_index(drop=True)
    logX = log1p_cp10k(X_counts, epi_obs.total_umi.to_numpy())
    print(f"[expr] epithelial matrix {logX.shape}", flush=True)

    # Pearson GRN for each TF
    pearson = {}
    pearson_rows = []
    for tf in ALL_TFS:
        pr = pearson_targets(logX, genes, tf)
        pearson[tf] = pr
        if len(pr):
            pearson_rows.append(pr)
            print(f"[pearson] {tf} top5: " +
                  ", ".join(f"{a}={b:.3f}" for a, b in pr.head(5)[["target", "pearson_r"]].values),
                  flush=True)
            for gene in TARGETS:
                hit = pr.loc[pr.target == gene]
                if len(hit):
                    print(f"         {tf} vs {gene}: r={hit.pearson_r.iloc[0]:.3f} "
                          f"fdr={hit.fdr.iloc[0]:.2e} rank={int(hit['rank'].iloc[0])}", flush=True)
    pearson_all = pd.concat(pearson_rows, ignore_index=True) if pearson_rows else pd.DataFrame()
    pearson_all.to_csv(RES / "pearson_tf_targets.tsv", sep="\t", index=False)

    grn = try_grnboost(logX, genes, FOCUS)
    if len(grn):
        grn.to_csv(RES / "grnboost2_adjacencies.tsv", sep="\t", index=False)
        print(f"[grnboost] edges={len(grn)}", flush=True)
    else:
        print("[grnboost] skipped or empty — Pearson+priors remain", flush=True)

    regulons = build_regulon_table(priors, pearson, grn, gene_set)
    regulons.to_csv(RES / "regulons.tsv", sep="\t", index=False)
    print(regulons[["regulon", "kind", "n_targets"]].to_string(index=False), flush=True)

    # AUCell on epithelial cells
    auc = {}
    for rec in regulons.itertuples(index=False):
        members = [g for g in (rec.targets.split(";") if rec.targets else []) if g in gene_set]
        mask = np.array([g in set(members) for g in genes])
        auc[rec.regulon] = aucell(logX, mask) if mask.any() else np.full(logX.shape[0], np.nan)
        print(f"[aucell] {rec.regulon} n={int(mask.sum())} mean={np.nanmean(auc[rec.regulon]):.4f}", flush=True)

    # also score TF RNA itself as a comparator (not a regulon)
    gix_e = {g: i for i, g in enumerate(genes)}
    for tf in ALL_TFS:
        if tf in gix_e:
            epi_obs[f"{tf}_log1p"] = logX[:, gix_e[tf]]
    for gene in TARGETS:
        if gene in gix_e:
            epi_obs[f"{gene}_log1p"] = logX[:, gix_e[gene]]
    for name, vec in auc.items():
        epi_obs[f"AUC_{name}"] = vec

    epi_obs.to_csv(RES / "epithelial_cells.tsv.gz", sep="\t", index=False, compression="gzip")

    # --- correlations ---
    corr_rows = []

    def add_corr(unit, split, compartment, regulon, gene, x, y, note):
        s = spear(x, y)
        corr_rows.append(dict(
            unit=unit, split=split, compartment=compartment, regulon=regulon,
            gene=gene, rho=s["rho"], p=s["p"], n=s["n"], note=note,
        ))

    splits = {
        "all": np.ones(len(epi_obs), dtype=bool),
        "LUAD": (epi_obs.histology == "LUAD").to_numpy(),
        "LUSC": (epi_obs.histology == "LUSC").to_numpy(),
    }
    compartments = {
        "epithelial": np.ones(len(epi_obs), dtype=bool),
        "malignant_like": epi_obs.malignant_like.to_numpy(),
    }
    auc_names = list(auc.keys()) + [f"{tf}_RNA" for tf in ALL_TFS]

    for split, smask in splits.items():
        for comp, cmask in compartments.items():
            m = smask & cmask
            if m.sum() < 20:
                continue
            for reg in auc_names:
                if reg.endswith("_RNA"):
                    tf = reg.replace("_RNA", "")
                    x = epi_obs.loc[m, f"{tf}_log1p"] if f"{tf}_log1p" in epi_obs else None
                else:
                    x = epi_obs.loc[m, f"AUC_{reg}"] if f"AUC_{reg}" in epi_obs else None
                if x is None:
                    continue
                for gene in TARGETS:
                    y = epi_obs.loc[m, f"{gene}_log1p"]
                    add_corr("cell", split, comp, reg, gene, x, y,
                             "exploratory; cells are not independent (pseudoreplication)")

            # sample-level means
            sub = epi_obs.loc[m].copy()
            for sample, g in sub.groupby("sample"):
                if len(g) < MIN_CELLS:
                    continue
            samp_rows = []
            for sample, g in sub.groupby("sample"):
                if len(g) < MIN_CELLS:
                    continue
                rec = dict(sample=sample, n_cells=len(g),
                           histology=g.histology.iloc[0], Patient=g.Patient.iloc[0])
                for gene in TARGETS:
                    rec[gene] = float(g[f"{gene}_log1p"].mean())
                for tf in ALL_TFS:
                    if f"{tf}_log1p" in g:
                        rec[f"{tf}_RNA"] = float(g[f"{tf}_log1p"].mean())
                for reg in auc:
                    rec[reg] = float(g[f"AUC_{reg}"].mean())
                samp_rows.append(rec)
            samp = pd.DataFrame(samp_rows)
            if samp.empty:
                continue
            samp["split"] = split
            samp["compartment"] = comp
            samp.to_csv(RES / f"sample_means_{split}_{comp}.tsv", sep="\t", index=False)
            for reg in auc_names:
                col = reg
                if col not in samp.columns:
                    continue
                for gene in TARGETS:
                    add_corr("sample", split, comp, reg, gene, samp[col], samp[gene],
                             f"primary unit=sample; min {MIN_CELLS} {comp} cells")

    corr = pd.DataFrame(corr_rows)
    corr.to_csv(RES / "correlations.tsv", sep="\t", index=False)

    # print the honest primary table
    prim = corr.loc[
        (corr.unit == "sample")
        & (corr.compartment == "malignant_like")
        & (corr.split.isin(["all", "LUAD"]))
        & (corr.regulon.isin(
            [f"{t}_prior" for t in FOCUS]
            + [f"{t}_pearson" for t in FOCUS]
            + [f"{t}_RNA" for t in FOCUS]
            + ["ELF3_GRHL1_KLF4_prior_union", "ELF3_GRHL1_KLF4_pearson_union"]
        ))
    ].copy()
    print("\n=== PRIMARY sample-level malignant-like ===", flush=True)
    if len(prim):
        show = prim.copy()
        show["rho"] = show["rho"].map(lambda x: f"{x:.3f}" if np.isfinite(x) else "NA")
        show["p"] = show["p"].map(lambda x: f"{x:.2e}" if np.isfinite(x) and x < 0.001 else (f"{x:.4f}" if np.isfinite(x) else "NA"))
        print(show[["split", "regulon", "gene", "rho", "p", "n"]].to_string(index=False), flush=True)
    else:
        print("no primary rows", flush=True)

    v = verdict(corr)
    print(f"\n[verdict] overall={v['overall']}  all={v['by_split']['all']}  LUAD={v['by_split']['LUAD']}", flush=True)

    # --- figures ---
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 5.2), sharey=True)
    regs_plot = [f"{t}_prior" for t in FOCUS] + [f"{t}_pearson" for t in FOCUS] + [
        "ELF3_GRHL1_KLF4_prior_union", "ELF3_GRHL1_KLF4_pearson_union"
    ]
    for ax, gene in zip(axes, TARGETS):
        ylabels = []
        ys = []
        for i, split in enumerate(["all", "LUAD"]):
            sub = corr[(corr.unit == "sample") & (corr.compartment == "malignant_like")
                       & (corr.split == split) & (corr.gene == gene)
                       & (corr.regulon.isin(regs_plot))]
            for j, reg in enumerate(regs_plot):
                r = sub[sub.regulon == reg]
                if r.empty:
                    continue
                y = j + (0 if i == 0 else 0)
                # plot all and LUAD as paired points on same y
        # rebuild with offset
        ytick = []
        yticklab = []
        for j, reg in enumerate(regs_plot):
            ytick.append(j)
            yticklab.append(reg)
            for i, split, col, dx in ((0, "all", "#1f4e79", -0.12), (1, "LUAD", "#c45c26", 0.12)):
                r = corr[(corr.unit == "sample") & (corr.compartment == "malignant_like")
                         & (corr.split == split) & (corr.gene == gene) & (corr.regulon == reg)]
                if r.empty or not np.isfinite(r.rho.iloc[0]):
                    continue
                ax.scatter(r.rho.iloc[0], j + dx, color=col, s=42, zorder=3,
                           label=f"{split} n={int(r.n.iloc[0])}" if j == 0 else None)
                ax.plot([0, r.rho.iloc[0]], [j + dx, j + dx], color=col, lw=1.4, alpha=0.8)
        ax.axvline(0, color="0.5", lw=0.8)
        ax.axvline(0.30, color="0.75", ls="--", lw=0.8)
        ax.set_xlabel(f"sample Spearman ρ  (AUCell vs {gene})")
        ax.set_yticks(ytick)
        ax.set_yticklabels(yticklab, fontsize=8)
        ax.set_title(gene)
        ax.set_xlim(-0.6, 1.0)
        handles, labels = ax.get_legend_handles_labels()
        if handles:
            ax.legend(handles[:2], labels[:2], fontsize=8, loc="lower right")
    fig.suptitle("GSE207422 malignant-like epithelium — regulon AUCell vs TACSTD2/CLDN4\n"
                 "sample means (n = samples with ≥20 cells). Dashed line = pre-specified ρ=0.30.",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(RES / "fig_sample_rho_forest.png", dpi=160)
    fig.savefig(RES / "fig_sample_rho_forest.pdf")
    plt.close(fig)

    # scatter: ELF3 prior AUC vs CLDN4, all vs LUAD
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.2))
    for ax, split in zip(axes, ("all", "LUAD")):
        p = RES / f"sample_means_{split}_malignant_like.tsv"
        if not p.exists():
            ax.set_title(f"{split}: no table")
            continue
        samp = pd.read_csv(p, sep="\t")
        if "ELF3_prior" not in samp.columns:
            ax.set_title(f"{split}: no ELF3_prior")
            continue
        c = np.where(samp.histology == "LUAD", "#c45c26", "#1f4e79")
        ax.scatter(samp["ELF3_prior"], samp["CLDN4"], c=c, s=46, edgecolor="k", lw=0.4)
        s = spear(samp["ELF3_prior"], samp["CLDN4"])
        ax.set_xlabel("mean ELF3_prior AUCell")
        ax.set_ylabel("mean CLDN4 log1p CP10k")
        ax.set_title(f"{split}  ρ={s['rho']:.3f}  p={s['p']:.3g}  n={s['n']}")
    fig.suptitle("ELF3 public-prior regulon vs CLDN4 (sample means, malignant-like)", fontsize=10)
    fig.tight_layout()
    fig.savefig(RES / "fig_elf3_prior_vs_cldn4.png", dpi=160)
    fig.savefig(RES / "fig_elf3_prior_vs_cldn4.pdf")
    plt.close(fig)

    # cell vs sample honesty panel
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    sub = corr[(corr.compartment == "malignant_like") & (corr.split == "all")
               & (corr.regulon.isin(regs_plot)) & (corr.gene == "CLDN4")]
    cell = sub[sub.unit == "cell"]
    samp = sub[sub.unit == "sample"]
    merged = cell.merge(samp, on="regulon", suffixes=("_cell", "_sample"))
    if len(merged):
        ax.scatter(merged.rho_cell, merged.rho_sample, s=40)
        for rec in merged.itertuples():
            ax.annotate(rec.regulon, (rec.rho_cell, rec.rho_sample), fontsize=6, xytext=(4, 2),
                        textcoords="offset points")
        ax.axhline(0, color="0.5", lw=0.6)
        ax.axvline(0, color="0.5", lw=0.6)
        ax.set_xlabel("cell-level ρ (exploratory, large n)")
        ax.set_ylabel("sample-level ρ (honest n)")
        ax.set_title("CLDN4 · all histologies · malignant-like")
    fig.tight_layout()
    fig.savefig(RES / "fig_cell_vs_sample_rho.png", dpi=160)
    fig.savefig(RES / "fig_cell_vs_sample_rho.pdf")
    plt.close(fig)

    summary = {
        "dataset": "GSE207422",
        "pmid": "36869384",
        "question": "Do ELF3/GRHL1/KLF4 regulons track TACSTD2/CLDN4 in public lung epithelium?",
        "a10_taken_as_given": "ELF3–CLDN4 bulk RNA couple; this test is regulon AUCell, not bulk coexpression and not ChIP.",
        "pyscenic": {"importable": pyscenic_ok, "version": pyscenic_ver,
                     "cistarget_run": False, "chip_peaks_invented": False},
        "n_genes_in_matrix": n_genes_full,
        "n_cells_raw": n,
        "n_epithelial": int(epi.sum()),
        "n_malignant_like": int(mal.sum()),
        "nl_cut": nl_cut,
        "n_genes_scored": len(genes),
        "histology_cells": epi_obs.histology.value_counts(dropna=False).to_dict(),
        "verdict": v["overall"],
        "verdict_by_split": v["by_split"],
        "primary_unit": "sample (malignant-like, min 20 cells)",
        "cell_level": "exploratory only",
        "priors_list_TACSTD2_or_CLDN4": False,
        "notes": [
            "Author CopyKAT malignant labels are not on GEO; malignant-like is a marker proxy.",
            "Public priors do not list TACSTD2 or CLDN4 as ELF3/GRHL1/KLF4 targets.",
            "Regulons hold out TACSTD2/CLDN4 so AUCell vs those genes is not circular.",
            "GRHL1 public prior is tiny (2–3 genes); Pearson regulon is the informative GRHL1 test.",
            "GSE253013 (9.3 GB) and GSE131907 (multi-GB) were not required for the LUAD vs all split.",
        ],
    }
    (RES / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")
    (RES / "verdict.json").write_text(json.dumps(v, indent=2, default=str) + "\n")
    print(f"[done] {RES}", flush=True)


if __name__ == "__main__":
    sys.exit(main() or 0)
