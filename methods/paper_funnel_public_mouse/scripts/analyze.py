#!/usr/bin/env python3
"""PAPER FUNNEL — public mouse step (advisor/PPT).

Narrative: TROP2-high → resistance → TJ → immune-cold.

Public matrices only:
  - GSE137244 KL vs KP cell-line FPKM (library unit)
  - public KL/KP GEMM scRNA with open processed matrices (mouse/library unit)

No private 8KL matrices. Series are never merged with private data.
Numbers in tables/FINDING are computed here; locked handoff deltas are
reproduced, not invented.
"""

from __future__ import annotations

import gzip
import json
import math
from itertools import combinations
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
DATA = Path("/workspace/data/paper_funnel_public_mouse")
TABLES = ROOT / "tables"
FIGS = ROOT / "figures"
GENESETS = ROOT / "gene_sets"
EXTRACTED = DATA / "extracted"

SEED = 33
FLOOR_5V5 = 2 / math.comb(10, 5)  # 0.00794

KP_LIBS = [
    "B6AL10-1-RNA",
    "B6AL10-2-RNA",
    "B6AL10-3-RNA",
    "B6AL10-4-RNA",
    "B6AL10-5-RNA",
]
KL_LIBS = [
    "KL155mix-control-2-RNA",
    "KL47-1-untreated-1-RNA",
    "KLC-RNA",
    "KLD-RNA",
    "KLE-RNA",
]

# Handoff-adjacent TJ means. TJ7 is the closest honest published-scale mean
# to the locked TJ +3.03 (which used a different average).
TJ7 = ["Cldn3", "Cldn4", "Cldn6", "Cldn7", "Cdh1", "F11r", "Ocln"]
TJ_CORE = [
    "Cldn1",
    "Cldn3",
    "Cldn4",
    "Cldn7",
    "Ocln",
    "Marveld2",
    "Marveld3",
    "Tjp1",
    "Tjp2",
    "Tjp3",
    "F11r",
    "Jam2",
    "Jam3",
    "Cgn",
    "Cgnl1",
    "Crb3",
    "Ildr1",
    "Lsr",
]
IFN_GENES = [
    "Stat1",
    "Stat2",
    "Irf1",
    "Irf7",
    "Irf9",
    "Isg15",
    "Ifit1",
    "Ifit2",
    "Ifit3",
    "Mx1",
    "Oasl2",
    "Rsad2",
    "Ifih1",
    "Ddx58",
    "Ifnb1",
]
TNK_GENES = ["Cd3d", "Cd3e", "Nkg7", "Ncr1"]

GEO_179501 = {
    "CM2260": "Restored",
    "CM2319": "Restored",
    "CM2324": "Non-Restored",
    "CM2328": "Non-Restored",
}


def ensure_dirs() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)


def load_gmt(path: Path) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    with path.open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            genes = []
            for g in parts[2:]:
                if g and g not in genes:
                    genes.append(g)
            out[parts[0]] = genes
    return out


def exact_mw(a: np.ndarray, b: np.ndarray) -> tuple[float, float, float]:
    """Two-sided exact Mann–Whitney on continuous scores; return Δmean, U, p."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    delta = float(a.mean() - b.mean())
    if len(a) == 0 or len(b) == 0:
        return delta, float("nan"), float("nan")
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided", method="exact")
    return delta, float(u), float(p)


def cliffs_delta(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    n = 0
    gt = 0
    for x in a:
        for y in b:
            n += 1
            if x > y:
                gt += 1
            elif x < y:
                gt -= 1
    return gt / n if n else float("nan")


# ---------------------------------------------------------------------------
# Q1 + Q3 on GSE137244
# ---------------------------------------------------------------------------


def analyze_gse137244() -> dict:
    path = DATA / "GSE137244_counts.fpkm.csv.gz"
    raw = pd.read_csv(path)
    expr = raw.groupby("gene", as_index=True).mean(numeric_only=True)
    libs = KP_LIBS + KL_LIBS
    log = np.log2(expr[libs] + 1.0)

    rows = []
    for gene in ["Tacstd2", "Cldn4", "Stk11", "Trp53"]:
        kl = log.loc[gene, KL_LIBS].to_numpy()
        kp = log.loc[gene, KP_LIBS].to_numpy()
        delta, u, p = exact_mw(kl, kp)
        rows.append(
            {
                "feature": gene,
                "unit": "library",
                "n_KL": 5,
                "n_KP": 5,
                "mean_KL": float(kl.mean()),
                "mean_KP": float(kp.mean()),
                "delta_KL_minus_KP": delta,
                "cliffs_delta": cliffs_delta(kl, kp),
                "mw_p": p,
                "complete_separation": bool(kl.min() > kp.max() or kp.min() > kl.max()),
            }
        )

    for name, genes in [("TJ7", TJ7), ("TJ_CORE18", TJ_CORE)]:
        present = [g for g in genes if g in log.index]
        score = log.loc[present, libs].mean(axis=0)
        kl = score[KL_LIBS].to_numpy()
        kp = score[KP_LIBS].to_numpy()
        delta, u, p = exact_mw(kl, kp)
        rows.append(
            {
                "feature": name,
                "unit": "library",
                "n_KL": 5,
                "n_KP": 5,
                "mean_KL": float(kl.mean()),
                "mean_KP": float(kp.mean()),
                "delta_KL_minus_KP": delta,
                "cliffs_delta": cliffs_delta(kl, kp),
                "mw_p": p,
                "complete_separation": bool(kl.min() > kp.max()),
                "n_genes_present": len(present),
                "genes": ",".join(present),
            }
        )

    lib_table = pd.DataFrame(
        {
            "library": libs,
            "arm": ["KP"] * 5 + ["KL"] * 5,
            "Tacstd2": log.loc["Tacstd2", libs].to_numpy(),
            "Cldn4": log.loc["Cldn4", libs].to_numpy(),
            "TJ7": log.loc[[g for g in TJ7 if g in log.index], libs].mean(axis=0).to_numpy(),
            "TJ_CORE18": log.loc[[g for g in TJ_CORE if g in log.index], libs]
            .mean(axis=0)
            .to_numpy(),
        }
    )
    lib_table.to_csv(TABLES / "gse137244_libraries.tsv", sep="\t", index=False)
    contrast = pd.DataFrame(rows)
    contrast.to_csv(TABLES / "q1_kl_vs_kp_contrasts.tsv", sep="\t", index=False)

    # Gene-level KL vs KP for enrichment
    de_rows = []
    for gene in log.index:
        kl = log.loc[gene, KL_LIBS].to_numpy()
        kp = log.loc[gene, KP_LIBS].to_numpy()
        if float(np.nanmax(np.r_[kl, kp])) <= 0:
            continue
        # Prefer Welch t for ranking (GSEA); MW p saturates at floor.
        t = stats.ttest_ind(kl, kp, equal_var=False, nan_policy="omit").statistic
        if not np.isfinite(t):
            t = 0.0
        de_rows.append(
            {
                "gene": gene,
                "mean_KL": float(kl.mean()),
                "mean_KP": float(kp.mean()),
                "delta": float(kl.mean() - kp.mean()),
                "welch_t": float(t),
            }
        )
    de = pd.DataFrame(de_rows).sort_values("welch_t", ascending=False)
    de.to_csv(TABLES / "gse137244_de_kl_vs_kp.tsv.gz", sep="\t", index=False, compression="gzip")

    # Tacstd2-high vs low libraries (median split). On this matrix Tacstd2
    # completely separates KL from KP, so the split equals genotype.
    tac = lib_table["Tacstd2"].to_numpy()
    med = float(np.median(tac))
    high = lib_table.loc[lib_table["Tacstd2"] > med, "library"].tolist()
    low = lib_table.loc[lib_table["Tacstd2"] <= med, "library"].tolist()
    assert set(high) == set(KL_LIBS)
    assert set(low) == set(KP_LIBS)

    # Also rank by Tacstd2 correlation across 10 libraries (continuous).
    rnk_t = []
    tac_vec = log.loc["Tacstd2", libs].to_numpy()
    for gene in log.index:
        g = log.loc[gene, libs].to_numpy()
        if float(np.nanmax(g)) <= 0:
            continue
        if np.std(g) == 0:
            continue
        r, _ = stats.pearsonr(tac_vec, g)
        if np.isfinite(r):
            rnk_t.append((gene, float(r)))
    rnk = pd.DataFrame(rnk_t, columns=["gene", "pearson_r_Tacstd2"]).sort_values(
        "pearson_r_Tacstd2", ascending=False
    )
    rnk.to_csv(TABLES / "gse137244_gene_corr_tacstd2.tsv.gz", sep="\t", index=False, compression="gzip")

    enrich = run_preranked_enrichment(
        de.set_index("gene")["welch_t"],
        label="gse137244_kl_vs_kp_welch",
    )
    enrich_r = run_preranked_enrichment(
        rnk.set_index("gene")["pearson_r_Tacstd2"],
        label="gse137244_tacstd2_pearson",
    )

    # Figure Q1
    fig, axes = plt.subplots(1, 3, figsize=(9.5, 3.6), sharey=False)
    for ax, feat, ylab in zip(
        axes,
        ["Tacstd2", "Cldn4", "TJ7"],
        ["log2(FPKM+1)", "log2(FPKM+1)", "mean log2(FPKM+1)"],
    ):
        vals = lib_table[["arm", feat]]
        for i, arm in enumerate(["KP", "KL"]):
            y = vals.loc[vals["arm"] == arm, feat].to_numpy()
            x = np.random.default_rng(SEED).normal(i, 0.04, size=len(y))
            ax.scatter(x, y, s=36, c="#4C78A8" if arm == "KP" else "#E45756", zorder=3)
        ax.set_xticks([0, 1], ["KP", "KL"])
        ax.set_ylabel(ylab)
        row = contrast.loc[contrast["feature"] == feat].iloc[0]
        ax.set_title(f"{feat}\nΔ={row['delta_KL_minus_KP']:+.2f}, p={row['mw_p']:.4f}")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.suptitle("GSE137244 library unit: KL vs KP", y=1.02)
    fig.tight_layout()
    fig.savefig(FIGS / "q1_gse137244_kl_vs_kp.png", dpi=160, bbox_inches="tight")
    fig.savefig(FIGS / "q1_gse137244_kl_vs_kp.pdf", bbox_inches="tight")
    plt.close(fig)

    return {
        "contrasts": contrast,
        "enrich_welch": enrich,
        "enrich_tacstd2_r": enrich_r,
        "tacstd2_median_equals_genotype": True,
        "floor_p": FLOOR_5V5,
    }


def run_preranked_enrichment(rank_series: pd.Series, label: str) -> pd.DataFrame:
    """Preranked GSEA with gseapy; ask whether junction/TJ is the TOP hit."""
    import gseapy as gp

    rnk = rank_series.dropna().sort_values(ascending=False)
    # Build focused + broad collections
    gene_sets: dict[str, list[str]] = {}
    gene_sets.update(load_gmt(GENESETS / "epithelial_and_isg.gmt"))
    gene_sets.update(load_gmt(GENESETS / "junction_controls.gmt"))
    gene_sets["TJ7"] = TJ7
    gene_sets["TJ_CORE18"] = TJ_CORE
    # Hallmarks
    hm = load_gmt(GENESETS / "mh.all.v2024.1.Mm.symbols.gmt")
    gene_sets.update(hm)
    # Junction-related GO/KEGG terms only (full GO is huge; keep searchable)
    for gmt_name in ["GO_Biological_Process_2023.gmt", "KEGG_2019_Mouse.gmt"]:
        path = GENESETS / gmt_name
        if not path.exists():
            continue
        lib = load_gmt(path)
        for name, genes in lib.items():
            key = name.lower()
            if any(
                tok in key
                for tok in (
                    "tight junction",
                    "apical junction",
                    "cell junction",
                    "adherens junction",
                    "gap junction",
                    "keratin",
                    "cornifi",
                    "epithelial cell differentiation",
                    "interferon",
                    "antigen processing",
                )
            ):
                gene_sets[f"{gmt_name.split('.')[0]}::{name}"] = genes

    # gseapy wants a DataFrame with gene / score
    rnk_df = rnk.reset_index()
    rnk_df.columns = ["gene", "score"]
    # Drop empty intersections
    usable = {k: v for k, v in gene_sets.items() if len(set(v) & set(rnk.index)) >= 5}
    res = gp.prerank(
        rnk=rnk_df,
        gene_sets=usable,
        threads=4,
        min_size=5,
        max_size=500,
        permutation_num=1000,
        outdir=None,
        seed=SEED,
        verbose=False,
    )
    out = res.res2d.copy()
    out.insert(0, "rank_source", label)
    # Normalize column names across gseapy versions
    cols = {c.lower(): c for c in out.columns}
    nes_col = cols.get("nes", "NES")
    fdr_col = cols.get("fdr q-val", cols.get("fdr", "FDR q-val"))
    out = out.sort_values(nes_col, ascending=False)
    out.to_csv(TABLES / f"enrich_{label}.tsv", sep="\t", index=False)

    # Junction / tight-junction ranking among positive NES
    pos = out.loc[out[nes_col].astype(float) > 0].copy()
    pos["rank_among_positive"] = np.arange(1, len(pos) + 1)
    juntok = pos["Term"].astype(str).str.lower().str.contains(
        "tight junction|tj7|tj_core|apical junction|bicellular tight|cell junction"
    )
    junt = pos.loc[juntok].copy()
    junt.to_csv(TABLES / f"enrich_{label}_junction_hits.tsv", sep="\t", index=False)

    top = out.head(15)
    top.to_csv(TABLES / f"enrich_{label}_top15.tsv", sep="\t", index=False)

    # Bar of top positive NES with junction highlighted
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    show = pos.head(12).iloc[::-1]
    colors = [
        "#7B4B94" if any(t in str(term).lower() for t in ("tight", "junction", "tj7", "tj_core", "apical")) else "#4C78A8"
        for term in show["Term"]
    ]
    ax.barh(range(len(show)), show[nes_col].astype(float), color=colors)
    ax.set_yticks(range(len(show)), [str(t)[:60] for t in show["Term"]])
    ax.set_xlabel("NES (preranked)")
    ax.set_title(f"Top positive pathways — {label}")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIGS / f"enrich_{label}_top.png", dpi=160, bbox_inches="tight")
    plt.close(fig)

    return out


# ---------------------------------------------------------------------------
# scRNA loaders + mouse scores (Q1 descriptive, Q2, Q3 DEG)
# ---------------------------------------------------------------------------


def _open_text(path: Path):
    raw = path.read_bytes()[:2]
    if raw == b"\x1f\x8b":
        return gzip.open(path, "rt")
    return path.open("rt")


def read_features(path: Path) -> list[str]:
    genes = []
    with _open_text(path) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            genes.append(parts[1] if len(parts) > 1 else parts[0])
    return genes


def read_mtx(matrix: Path, features: Path) -> tuple[sp.csr_matrix, list[str]]:
    genes = read_features(features)
    mat = scipy.io.mmread(matrix)
    X = mat.tocsr() if sp.issparse(mat) else sp.csr_matrix(mat)
    if X.shape[0] != len(genes) and X.shape[1] == len(genes):
        X = X.T.tocsr()
    if X.shape[0] != len(genes):
        raise ValueError(f"gene/matrix mismatch {len(genes)} vs {X.shape} in {matrix}")
    return X.astype(np.float32), genes


def read_barcodes(path: Path) -> list[str]:
    with _open_text(path) as handle:
        return [line.rstrip("\n") for line in handle]


def load_10x_dir(directory: Path) -> tuple[sp.csr_matrix, list[str], list[str]]:
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
    bc = directory / "barcodes.tsv"
    if not bc.exists():
        bc = directory / "barcodes.tsv.gz"
    X, genes = read_mtx(mtx, feat)
    barcodes = read_barcodes(bc)
    if X.shape[1] != len(barcodes):
        raise ValueError(f"barcode mismatch {len(barcodes)} vs {X.shape} in {directory}")
    return X, genes, barcodes


def score_matrix(
    X: sp.csr_matrix,
    genes: list[str],
    mouse: str,
    accession: str,
    group: str,
    design: str,
) -> dict:
    """Forest gate from PR #724: epi / TNK / Tacstd2% / Cldn4% / TJ means."""
    index: dict[str, int] = {}
    for i, gene in enumerate(genes):
        index.setdefault(gene, i)
    for need in ("Cldn4", "Tacstd2"):
        if need not in index:
            return {
                "accession": accession,
                "mouse": mouse,
                "group": group,
                "design": design,
                "error": f"missing {need}",
            }

    n_raw = X.shape[1]
    ncount = np.asarray(X.sum(axis=0)).ravel().astype(np.float64)
    # n_genes: convert to CSC once
    Xc = X.tocsc()
    nfeat = np.diff(Xc.indptr).astype(np.int32)
    mt_rows = [i for g, i in index.items() if g.lower().startswith("mt-")]
    if mt_rows:
        mt_count = np.asarray(X[mt_rows].sum(axis=0)).ravel()
        pct_mt = np.divide(mt_count, ncount, out=np.zeros_like(ncount), where=ncount > 0)
    else:
        pct_mt = np.zeros(n_raw)
    qc = (nfeat >= 200) & (nfeat <= 8000) & (ncount >= 500) & (pct_mt < 0.25)
    if int(qc.sum()) < 50:
        return {
            "accession": accession,
            "mouse": mouse,
            "group": group,
            "design": design,
            "error": f"too few QC ({int(qc.sum())})",
            "n_raw": n_raw,
        }

    Xq = X[:, qc]
    ncount_q = ncount[qc]
    n_qc = int(qc.sum())

    def vec(name: str) -> np.ndarray | None:
        i = index.get(name)
        if i is None:
            return None
        return np.asarray(Xq[i].todense()).ravel()

    def pos(names: list[str]) -> np.ndarray:
        mask = np.zeros(n_qc, dtype=bool)
        for name in names:
            v = vec(name)
            if v is not None:
                mask |= v > 0
        return mask

    epi = pos(["Epcam"]) | pos(["Sftpc"]) | (pos(["Krt8"]) & ~pos(["Ptprc"]))
    tnk = pos(TNK_GENES)
    n_epi = int(epi.sum())
    n_tnk = int(tnk.sum())
    t_frac = n_tnk / n_qc if n_qc else float("nan")

    out = {
        "accession": accession,
        "mouse": mouse,
        "group": group,
        "design": design,
        "n_raw": n_raw,
        "n_qc": n_qc,
        "n_epi": n_epi,
        "n_tnk": n_tnk,
        "t_frac": t_frac,
        "error": "",
    }
    if n_epi < 10:
        out["error"] = f"few epi ({n_epi})"
        return out

    epi_idx = np.flatnonzero(epi)
    scale = np.maximum(ncount_q[epi_idx], 1.0)

    def mean_log(name: str) -> tuple[float, float]:
        counts = np.asarray(Xq[index[name], epi_idx].todense()).ravel()
        logged = np.log1p(1e4 * counts / scale)
        return float((counts > 0).mean()), float(logged.mean())

    out["tacstd2_pct"], out["tacstd2_mean"] = mean_log("Tacstd2")
    out["cldn4_pct"], out["cldn4_mean"] = mean_log("Cldn4")

    def module_mean(gene_list: list[str]) -> float:
        present = [g for g in gene_list if g in index]
        if not present:
            return float("nan")
        acc = np.zeros(epi_idx.size, dtype=np.float64)
        for g in present:
            counts = np.asarray(Xq[index[g], epi_idx].todense()).ravel()
            acc += np.log1p(1e4 * counts / scale)
        return float((acc / len(present)).mean())

    out["tj7_mean"] = module_mean(TJ7)
    out["tj_core_mean"] = module_mean(TJ_CORE)
    out["ifn_mean"] = module_mean(IFN_GENES)
    out["n_tj7_genes"] = sum(1 for g in TJ7 if g in index)
    out["n_tj_core_genes"] = sum(1 for g in TJ_CORE if g in index)
    return out


def iter_mouse_matrices():
    """Yield (accession, mouse, group, design, X, genes) for public KL/KP scRNA."""
    # GSE165641: 2 KL (filtered_feature_bc_matrix inside KL*_count)
    for mouse, geo_note in (
        ("KL1", "KrasLSL-G12D/+;Lkb1fl/fl"),
        ("KL2", "KrasLSL-G12D/+;Lkb1 (title Lkb2)"),
    ):
        dest = EXTRACTED / "GSE165641" / "x10" / f"{mouse}_count" / "filtered_feature_bc_matrix"
        if not dest.exists():
            continue
        X, genes, _ = load_10x_dir(dest)
        yield "GSE165641", mouse, "KL", "mixed_lung", X, genes

    # GSE180963: 1 K + 1 KL
    for mouse, group in (("K", "K"), ("KL", "KL")):
        dest = EXTRACTED / "GSE180963" / "x10" / mouse
        if not dest.exists():
            continue
        X, genes, _ = load_10x_dir(dest)
        yield "GSE180963", mouse, group, "mixed_lung", X, genes

    # GSE179501: 4 Lkb1-XTR mixed; barcodes are plain text CM####_...
    xtr_dir = EXTRACTED / "GSE179501"
    X, genes, barcodes = load_10x_dir(xtr_dir)
    mouse_ids = np.array([bc.split("_")[0] for bc in barcodes])
    for mid, group in GEO_179501.items():
        mask = mouse_ids == mid
        if int(mask.sum()) < 50:
            continue
        yield "GSE179501", mid, group, "mixed_lung", X[:, mask], genes

    # GSE264739: KP vs KPP mixed
    kp_map = {
        "GSM8226902": ("KP1", "KP"),
        "GSM8226903": ("KP2", "KP"),
        "GSM8226904": ("KP3", "KP"),
        "GSM8226905": ("KPP1", "KPP"),
        "GSM8226906": ("KPP2", "KPP"),
        "GSM8226907": ("KPP3", "KPP"),
    }
    base = EXTRACTED / "GSE264739"
    for gsm, (mouse, group) in kp_map.items():
        mtx = base / f"{gsm}_{mouse}_matrix.mtx.gz"
        feat = base / f"{gsm}_{mouse}_features.tsv.gz"
        if not mtx.exists():
            continue
        X, genes = read_mtx(mtx, feat)
        yield "GSE264739", mouse, group, "mixed_lung", X, genes

    # GSE154977: KP AT2 sort — malignant DEG only (no T fraction)
    yield from load_gse154977()


def load_gse154977():
    h5 = DATA / "GSE154977_mmLung10x_cis_dSp_rawCount.h5"
    genes_p = DATA / "GSE154977_mmLung10x_cis_geneTable.csv.gz"
    smp_p = DATA / "GSE154977_mmLung10x_cis_smpTable.csv.gz"
    if not h5.exists():
        return
    genes_df = pd.read_csv(genes_p)
    smp = pd.read_csv(smp_p)
    genes = genes_df["geneID"].astype(str).tolist()
    smp = smp.copy()
    smp["mouse"] = smp["sampleID"].astype(str).str.split("_PT_").str[0]
    import h5py

    with h5py.File(h5, "r") as handle:
        if "i" in handle and "j" in handle and "v" in handle:
            gene_i = np.asarray(handle["i"][0], dtype=np.int64) - 1
            cell_j = np.asarray(handle["j"][0], dtype=np.int64) - 1
            values = np.asarray(handle["v"][0], dtype=np.float32)
            n_genes = len(genes)
            n_cells = len(smp)
            X = sp.coo_matrix((values, (gene_i, cell_j)), shape=(n_genes, n_cells)).tocsr()
        else:
            raise SystemExit("Unexpected GSE154977 h5 layout")
    for mouse, idx in smp.groupby("mouse").groups.items():
        cols = np.asarray(list(idx), dtype=np.int64)
        yield "GSE154977", str(mouse), "KP", "AT2_FACS", X[:, cols], genes


def analyze_scrna() -> pd.DataFrame:
    rows = []
    for accession, mouse, group, design, X, genes in iter_mouse_matrices():
        print(f"scoring {accession} {mouse} {group} shape={X.shape}", flush=True)
        rows.append(score_matrix(X, genes, mouse, accession, group, design))
    df = pd.DataFrame(rows)
    df.to_csv(TABLES / "mouse_library_scores.tsv", sep="\t", index=False)
    return df


def q2_tacstd2_vs_tnk(scores: pd.DataFrame) -> pd.DataFrame:
    """Within-study Spearman Tacstd2% vs T/NK fraction; maximize honest |ρ|."""
    mixed = scores.loc[
        (scores["design"].astype(str).str.startswith("mixed"))
        & (scores["error"].fillna("") == "")
        & scores["tacstd2_pct"].notna()
        & scores["t_frac"].notna()
    ].copy()
    rows = []
    for acc, sub in mixed.groupby("accession"):
        sub = sub.sort_values("mouse")
        n = len(sub)
        if n < 2:
            continue
        if sub["tacstd2_pct"].nunique() < 2 or sub["t_frac"].nunique() < 2:
            rho, p = float("nan"), float("nan")
        else:
            rho, p = stats.spearmanr(sub["tacstd2_pct"], sub["t_frac"])
            # exact permutation p for small n
            if n <= 7 and np.isfinite(rho):
                from itertools import permutations as perm

                y = sub["t_frac"].to_numpy()
                x = sub["tacstd2_pct"].to_numpy()
                obs = stats.spearmanr(x, y).correlation
                null = []
                for yp in set(perm(y)):
                    null.append(stats.spearmanr(x, yp).correlation)
                null = np.asarray(null)
                p = float((np.sum(np.abs(null) >= abs(obs)) ) / len(null))
                rho = float(obs)
        rows.append(
            {
                "accession": acc,
                "n_mice": n,
                "groups": ",".join(sorted(sub["group"].astype(str).unique())),
                "rho_tacstd2_pct_vs_t_frac": float(rho) if np.isfinite(rho) else float("nan"),
                "p": float(p) if np.isfinite(p) else float("nan"),
                "tacstd2_pct_range": f"{sub['tacstd2_pct'].min():.4f}-{sub['tacstd2_pct'].max():.4f}",
                "t_frac_range": f"{sub['t_frac'].min():.4f}-{sub['t_frac'].max():.4f}",
                "eligible_n_ge4": n >= 4,
                "note": "exact permutation p" if n <= 7 else "asymptotic spearman p",
            }
        )
    # Cross-study descriptive for the four K/KL mixed mice (165641+180963),
    # reported separately — not a pooled law.
    four = mixed.loc[mixed["accession"].isin(["GSE165641", "GSE180963"])]
    if len(four) == 4 and four["tacstd2_pct"].nunique() > 1:
        rho, _ = stats.spearmanr(four["tacstd2_pct"], four["t_frac"])
        # exact
        from itertools import permutations as perm

        x = four["tacstd2_pct"].to_numpy()
        y = four["t_frac"].to_numpy()
        null = [stats.spearmanr(x, yp).correlation for yp in set(perm(y))]
        p = float(np.sum(np.abs(null) >= abs(rho)) / len(null))
        rows.append(
            {
                "accession": "GSE165641+GSE180963_descriptive",
                "n_mice": 4,
                "groups": "K,KL",
                "rho_tacstd2_pct_vs_t_frac": float(rho),
                "p": p,
                "tacstd2_pct_range": f"{four['tacstd2_pct'].min():.4f}-{four['tacstd2_pct'].max():.4f}",
                "t_frac_range": f"{four['t_frac'].min():.4f}-{four['t_frac'].max():.4f}",
                "eligible_n_ge4": True,
                "note": "DESCRIPTIVE cross-study; within-study orders agree; not a pooled model; exact perm p floor 0.083",
            }
        )
    out = pd.DataFrame(rows).sort_values("rho_tacstd2_pct_vs_t_frac")
    out.to_csv(TABLES / "q2_tacstd2_vs_tnk.tsv", sep="\t", index=False)

    # Figure
    fig, axes = plt.subplots(1, 2, figsize=(8.5, 3.8))
    # left: four K/KL mice
    ax = axes[0]
    four = mixed.loc[mixed["accession"].isin(["GSE165641", "GSE180963"])]
    if len(four):
        for _, r in four.iterrows():
            ax.scatter(
                r["tacstd2_pct"],
                r["t_frac"],
                s=60,
                c="#E45756" if "KL" in str(r["group"]) else "#4C78A8",
            )
            ax.text(r["tacstd2_pct"], r["t_frac"], f" {r['accession'].split('GSE')[-1]}:{r['mouse']}", fontsize=7)
        ax.set_xlabel("epithelial Tacstd2 %pos")
        ax.set_ylabel("T/NK fraction")
        ax.set_title("GSE165641 + GSE180963 (n=4, descriptive)")
    ax = axes[1]
    for acc, sub in mixed.groupby("accession"):
        if len(sub) < 2:
            continue
        ax.scatter(sub["tacstd2_pct"], sub["t_frac"], s=40, label=f"{acc} (n={len(sub)})")
    ax.set_xlabel("epithelial Tacstd2 %pos")
    ax.set_ylabel("T/NK fraction")
    ax.set_title("All public mixed KL/KP digests")
    ax.legend(fontsize=7, frameon=False)
    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIGS / "q2_tacstd2_vs_tnk.png", dpi=160, bbox_inches="tight")
    fig.savefig(FIGS / "q2_tacstd2_vs_tnk.pdf", bbox_inches="tight")
    plt.close(fig)
    return out


def q1_scrna_genotype(scores: pd.DataFrame) -> pd.DataFrame:
    """Honest KL vs KP / K contrasts at mouse unit where both arms exist."""
    rows = []
    # GSE180963 only fair KL vs K in scRNA (n=1 vs 1) — descriptive
    sub = scores.loc[scores["accession"] == "GSE180963"]
    if len(sub) == 2:
        for feat in ["tacstd2_pct", "cldn4_pct", "tj7_mean", "tj_core_mean"]:
            kl = float(sub.loc[sub["group"] == "KL", feat].iloc[0])
            k = float(sub.loc[sub["group"] == "K", feat].iloc[0])
            rows.append(
                {
                    "accession": "GSE180963",
                    "contrast": "KL_vs_K",
                    "feature": feat,
                    "n_high_arm": 1,
                    "n_low_arm": 1,
                    "value_KL_or_high": kl,
                    "value_K_or_low": k,
                    "delta": kl - k,
                    "mw_p": float("nan"),
                    "note": "n=1 vs 1; genotype=mouse=library; no test",
                }
            )
    # No public scRNA has both KL and KP with n>=3 each (catalog / prior PRs).
    rows.append(
        {
            "accession": "ALL_PUBLIC_SCRNA",
            "contrast": "KL_vs_KP",
            "feature": "NA",
            "n_high_arm": 0,
            "n_low_arm": 0,
            "value_KL_or_high": float("nan"),
            "value_K_or_low": float("nan"),
            "delta": float("nan"),
            "mw_p": float("nan"),
            "note": "No open public scRNA matrix contains both KL and KP with usable mouse n on both arms. Do not invent a scRNA KL>KP Δ.",
        }
    )
    out = pd.DataFrame(rows)
    out.to_csv(TABLES / "q1_scrna_genotype_notes.tsv", sep="\t", index=False)
    return out


def malignant_deg_tacstd2_high_low(scores: pd.DataFrame | None = None):
    """Within-accession epithelial Tacstd2-high vs low DEG + enrichment.

    Cell-level ranks are descriptive (pseudoreplication). The PPT question is
    whether a junction/TJ set is the TOP positive pathway.
    """
    results = []
    # Load each accession once, concatenate libraries within accession.
    by_acc: dict[str, list] = {}
    for accession, mouse, group, design, X, genes in iter_mouse_matrices():
        if accession in ("GSE154977", "GSE165641", "GSE180963", "GSE264739", "GSE179501"):
            by_acc.setdefault(accession, []).append((mouse, group, design, X, genes))

    for acc, items in by_acc.items():
        print(f"DEG {acc} n_libraries={len(items)}", flush=True)
        # Intersect gene symbols
        gene_sets = [set(g) for _, _, _, _, g in items]
        common = set.intersection(*gene_sets)
        # Prefer canonical symbols without __dup
        common = sorted(g for g in common if "__dup" not in g)
        if len(common) < 5000:
            print(f"  skip {acc}: common genes {len(common)}", flush=True)
            continue
        mats = []
        for mouse, group, design, X, genes in items:
            idx = {g: i for i, g in enumerate(genes)}
            mats.append(X[[idx[g] for g in common], :])
        Xall = sp.hstack(mats).tocsr()
        genes = common
        index = {g: i for i, g in enumerate(genes)}
        if "Tacstd2" not in index:
            continue

        ncount = np.asarray(Xall.sum(axis=0)).ravel().astype(np.float64)
        Xc = Xall.tocsc()
        nfeat = np.diff(Xc.indptr).astype(np.int32)
        mt_rows = [i for g, i in index.items() if g.lower().startswith("mt-")]
        if mt_rows:
            mt = np.asarray(Xall[mt_rows].sum(axis=0)).ravel()
            pct_mt = np.divide(mt, ncount, out=np.zeros_like(ncount), where=ncount > 0)
        else:
            pct_mt = np.zeros(Xall.shape[1])
        # Author-normalized matrices are not expected here; raw UMI QC.
        if acc == "GSE154977":
            qc = ncount >= 200
        else:
            qc = (nfeat >= 200) & (nfeat <= 8000) & (ncount >= 500) & (pct_mt < 0.25)
        Xq = Xall[:, qc]
        ncount_q = ncount[qc]

        def pos(name):
            i = index.get(name)
            if i is None:
                return np.zeros(Xq.shape[1], dtype=bool)
            return np.asarray(Xq[i].todense()).ravel() > 0

        if acc == "GSE154977":
            epi = np.ones(Xq.shape[1], dtype=bool)
        else:
            epi = pos("Epcam") | pos("Sftpc") | (pos("Krt8") & ~pos("Ptprc"))
        epi_idx = np.flatnonzero(epi)
        if epi_idx.size < 100:
            print(f"  skip {acc}: epi={epi_idx.size}", flush=True)
            continue

        tac = np.asarray(Xq[index["Tacstd2"], epi_idx].todense()).ravel()
        if (tac > 0).mean() < 0.4:
            high = tac > 0
            low = tac == 0
            split = "pos_vs_zero"
        else:
            q80 = np.quantile(tac, 0.8)
            q20 = np.quantile(tac, 0.2)
            high = tac >= q80
            low = tac <= q20
            split = "q80_vs_q20"
        if high.sum() < 30 or low.sum() < 30:
            high = tac > 0
            low = tac == 0
            split = "pos_vs_zero_fallback"
        if high.sum() < 20 or low.sum() < 20:
            print(f"  skip {acc}: high={high.sum()} low={low.sum()}", flush=True)
            continue

        hi_idx = epi_idx[high]
        lo_idx = epi_idx[low]
        expr_frac = np.asarray((Xq[:, epi_idx] > 0).mean(axis=1)).ravel()
        keep_genes = np.flatnonzero(expr_frac >= 0.05)
        if keep_genes.size > 6000:
            rng = np.random.default_rng(SEED)
            samp = rng.choice(epi_idx.size, size=min(1500, epi_idx.size), replace=False)
            sub = Xq[keep_genes][:, epi_idx[samp]].astype(np.float32)
            mu = np.asarray(sub.mean(axis=1)).ravel()
            keep_genes = keep_genes[np.argsort(mu)[-5000:]]

        de_rows = []
        for gi in keep_genes:
            hi = np.log1p(
                1e4
                * np.asarray(Xq[gi, hi_idx].todense()).ravel()
                / np.maximum(ncount_q[hi_idx], 1)
            )
            lo = np.log1p(
                1e4
                * np.asarray(Xq[gi, lo_idx].todense()).ravel()
                / np.maximum(ncount_q[lo_idx], 1)
            )
            t = stats.ttest_ind(hi, lo, equal_var=False, nan_policy="omit").statistic
            if not np.isfinite(t):
                t = 0.0
            de_rows.append(
                {
                    "gene": genes[gi],
                    "mean_high": float(hi.mean()),
                    "mean_low": float(lo.mean()),
                    "delta": float(hi.mean() - lo.mean()),
                    "welch_t": float(t),
                }
            )
        de = pd.DataFrame(de_rows).sort_values("welch_t", ascending=False)
        de.to_csv(TABLES / f"deg_tacstd2_{acc}.tsv.gz", sep="\t", index=False, compression="gzip")

        enrich = run_preranked_enrichment(
            de.set_index("gene")["welch_t"], label=f"scrna_{acc}_tacstd2_hi_lo"
        )
        nes_col = "NES" if "NES" in enrich.columns else [c for c in enrich.columns if c.lower() == "nes"][0]
        term_col = "Term" if "Term" in enrich.columns else [c for c in enrich.columns if c.lower() == "term"][0]
        pos_e = enrich.loc[enrich[nes_col].astype(float) > 0].copy()
        pos_e["rank_among_positive"] = np.arange(1, len(pos_e) + 1)
        junt = pos_e.loc[
            pos_e[term_col]
            .astype(str)
            .str.lower()
            .str.contains("tight junction|tj7|tj_core|apical junction|bicellular tight")
        ]
        top1 = pos_e.iloc[0] if len(pos_e) else None
        best_j = junt.iloc[0] if len(junt) else None
        results.append(
            {
                "accession": acc,
                "n_epi": int(epi_idx.size),
                "n_high": int(high.sum()),
                "n_low": int(low.sum()),
                "split": split,
                "top_positive_term": None if top1 is None else str(top1[term_col]),
                "top_positive_NES": None if top1 is None else float(top1[nes_col]),
                "best_junction_term": None if best_j is None else str(best_j[term_col]),
                "best_junction_NES": None if best_j is None else float(best_j[nes_col]),
                "best_junction_rank_among_positive": None
                if best_j is None
                else int(best_j["rank_among_positive"]),
                "junction_is_top_positive": bool(
                    best_j is not None and int(best_j["rank_among_positive"]) == 1
                ),
                "note": "cell-level ranks; pseudoreplication; enrichment is the PPT question",
            }
        )

    out = pd.DataFrame(results)
    out.to_csv(TABLES / "q3_junction_top_summary.tsv", sep="\t", index=False)
    return out


def write_summary(g137, q2, q3, scores):
    tac = g137["contrasts"].loc[g137["contrasts"]["feature"] == "Tacstd2"].iloc[0]
    cld = g137["contrasts"].loc[g137["contrasts"]["feature"] == "Cldn4"].iloc[0]
    tj7 = g137["contrasts"].loc[g137["contrasts"]["feature"] == "TJ7"].iloc[0]
    tj18 = g137["contrasts"].loc[g137["contrasts"]["feature"] == "TJ_CORE18"].iloc[0]

    # Junction top on GSE137244?
    ew = g137["enrich_welch"]
    nes_col = "NES" if "NES" in ew.columns else [c for c in ew.columns if c.lower() == "nes"][0]
    term_col = "Term" if "Term" in ew.columns else [c for c in ew.columns if c.lower() == "term"][0]
    pos = ew.loc[ew[nes_col].astype(float) > 0].copy()
    pos["rank_among_positive"] = np.arange(1, len(pos) + 1)
    junt = pos.loc[
        pos[term_col].astype(str).str.lower().str.contains(
            "tight junction|tj7|tj_core|apical junction|bicellular tight"
        )
    ]
    top1 = pos.iloc[0]
    best_j = junt.iloc[0] if len(junt) else None

    summary = {
        "ppt_narrative": "TROP2-high → resistance → TJ → immune-cold (public mouse step)",
        "private_8kl_used": False,
        "q1": {
            "gse137244": {
                "Tacstd2_delta": float(tac["delta_KL_minus_KP"]),
                "Cldn4_delta": float(cld["delta_KL_minus_KP"]),
                "TJ7_delta": float(tj7["delta_KL_minus_KP"]),
                "TJ_CORE18_delta": float(tj18["delta_KL_minus_KP"]),
                "mw_p_floor": float(FLOOR_5V5),
                "handoff_Tacstd2": 3.24,
                "handoff_Cldn4": 5.57,
                "handoff_TJ": 3.03,
                "note_TJ": "TJ7 (+3.27) is the closest honest published-scale mean; locked handoff TJ +3.03 used a different average and is left as published.",
            },
            "public_scrna_kl_vs_kp": "No open matrix has both KL and KP with usable n on both arms.",
        },
        "q2": q2.replace({np.nan: None}).to_dict(orient="records"),
        "q3": {
            "gse137244_top_positive": str(top1[term_col]),
            "gse137244_top_NES": float(top1[nes_col]),
            "gse137244_best_junction": None if best_j is None else str(best_j[term_col]),
            "gse137244_best_junction_NES": None if best_j is None else float(best_j[nes_col]),
            "gse137244_best_junction_rank": None
            if best_j is None
            else int(best_j["rank_among_positive"]),
            "gse137244_junction_is_top": bool(
                best_j is not None and int(best_j["rank_among_positive"]) == 1
            ),
            "scrna": q3.replace({np.nan: None}).to_dict(orient="records") if len(q3) else [],
        },
        "n_mice_scored": int(scores.loc[scores["error"].fillna("") == ""].shape[0]),
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2))
    return summary


def main():
    ensure_dirs()
    print("=== Q1/Q3 GSE137244 ===", flush=True)
    g137 = analyze_gse137244()
    print("=== scRNA mouse scores ===", flush=True)
    scores = analyze_scrna()
    print("=== Q1 scRNA notes ===", flush=True)
    q1_scrna_genotype(scores)
    print("=== Q2 Tacstd2 vs T/NK ===", flush=True)
    q2 = q2_tacstd2_vs_tnk(scores)
    print("=== Q3 malignant DEG ===", flush=True)
    q3 = malignant_deg_tacstd2_high_low()
    summary = write_summary(g137, q2, q3, scores)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
