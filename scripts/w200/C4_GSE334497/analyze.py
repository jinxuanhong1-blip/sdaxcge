#!/usr/bin/env python3
"""C4 analog: GSE334497 4T1 Trop2 KO vs WT — Cldn4, Cxcl9, IFN/APM.

Author-normalized counts (GEO supplementary). Contrast = KO minus WT.
Does not claim tumor-intrinsic IFN: bulk immunocompetent tumors.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import os
import sys
import urllib.request
from itertools import combinations
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "results" / "w200" / "C4_GSE334497"
RAW = OUT / "raw"
TABLES = OUT / "tables"
FIGS = OUT / "figures"
CACHE = Path(os.environ.get("C4_GSE334497_CACHE", "/tmp/gse334497"))

COUNTS_URL = (
    "https://www.ncbi.nlm.nih.gov/geo/download/"
    "?acc=GSE334497&format=file&file=GSE334497_normalized_counts.csv.gz"
)
GENE2ENS_URL = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz"
GENEINFO_URL = (
    "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/"
    "Mus_musculus.gene_info.gz"
)

# Library names from GEO series matrix Sample_description
KO = ["KO162", "KO164", "KO165", "KO172", "RESUB-KO163R"]
WT = ["control170", "RESUB-171R", "RESUB-170R", "RESUB-169R", "RESUB-168R"]

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gene_sets import (  # noqa: E402
    CTRL_MYC,
    CTRL_OXPHOS,
    CXCL_IFN,
    FOCAL,
    ISG_CORE,
    MHC1_APM,
    QC_GENES,
    REPORT_ORDER,
    SETS,
    T_CYT,
    USER_CORE6,
)


def md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    print(f"download {url} -> {dest}", flush=True)
    tmp = dest.with_suffix(dest.suffix + ".tmp")
    urllib.request.urlretrieve(url, tmp)
    tmp.replace(dest)
    return dest


def build_symbol_map(ensembl_ids: list[str]) -> dict[str, str]:
    """NCBI gene2ensembl + gene_info → Ensembl gene ID to official symbol."""
    shipped = RAW / "ensembl_to_symbol.tsv"
    cache_map = CACHE / "ensembl_to_symbol.tsv"
    for cand in (shipped, cache_map):
        if cand.exists():
            m = {}
            with open(cand) as f:
                next(f)
                for line in f:
                    ens, sym = line.rstrip("\n").split("\t")
                    m[ens] = sym
            if sum(1 for e in ensembl_ids if e in m) > 0.7 * len(ensembl_ids):
                return m

    g2e = download(GENE2ENS_URL, CACHE / "gene2ensembl.gz")
    ginfo = download(GENEINFO_URL, CACHE / "Mus_musculus.gene_info.gz")
    gid2sym: dict[str, str] = {}
    with gzip.open(ginfo, "rt") as f:
        next(f)
        for line in f:
            p = line.split("\t")
            gid2sym[p[1]] = p[2]
    ens2gid: dict[str, str] = {}
    with gzip.open(g2e, "rt") as f:
        next(f)
        for line in f:
            p = line.split("\t")
            if p[0] != "10090":
                continue
            ens2gid[p[2]] = p[1]
    m = {}
    for ens in ensembl_ids:
        gid = ens2gid.get(ens)
        if gid and gid in gid2sym:
            m[ens] = gid2sym[gid]
    cache_map.parent.mkdir(parents=True, exist_ok=True)
    with open(cache_map, "w") as f:
        f.write("ensembl\tsymbol\n")
        for ens, sym in sorted(m.items()):
            f.write(f"{ens}\t{sym}\n")
    return m


def welch_t(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    if np.allclose(a, a[0]) and np.allclose(b, b[0]) and a[0] == b[0]:
        return 0.0, 1.0
    t, p = stats.ttest_ind(a, b, equal_var=False)
    if np.isnan(t) or np.isnan(p):
        return 0.0, 1.0
    return float(t), float(p)


def mwu(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    try:
        u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    except ValueError:
        return float("nan"), 1.0
    return float(u), float(p)


def cohens_d(a: np.ndarray, b: np.ndarray) -> float:
    n1, n2 = len(a), len(b)
    s1, s2 = a.var(ddof=1), b.var(ddof=1)
    sp = math.sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2))
    if sp == 0:
        return 0.0
    return float((a.mean() - b.mean()) / sp)


def exact_perm_p(scores: np.ndarray, is_ko: np.ndarray, observed: float) -> dict:
    """Exact two-group mean-difference permutation (5v5 → 252 splits)."""
    idx = np.arange(len(scores))
    n_ko = int(is_ko.sum())
    n = 0
    n_ge = 0
    n_le = 0
    n_abs = 0
    for combo in combinations(idx, n_ko):
        mask = np.zeros(len(scores), dtype=bool)
        mask[list(combo)] = True
        delta = scores[mask].mean() - scores[~mask].mean()
        n += 1
        if delta >= observed - 1e-15:
            n_ge += 1
        if delta <= observed + 1e-15:
            n_le += 1
        if abs(delta) >= abs(observed) - 1e-15:
            n_abs += 1
    return {
        "n_perm": n,
        "perm_p_up": n_ge / n,
        "perm_p_down": n_le / n,
        "perm_p_two": n_abs / n,
    }


def auc_from_mwu(u: float, n1: int, n2: int) -> float:
    if n1 * n2 == 0:
        return float("nan")
    return float(u / (n1 * n2))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    RAW.mkdir(parents=True, exist_ok=True)
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(parents=True, exist_ok=True)

    counts_gz = download(COUNTS_URL, CACHE / "GSE334497_normalized_counts.csv.gz")
    dest_counts = RAW / "GSE334497_normalized_counts.csv.gz"
    if not dest_counts.exists():
        dest_counts.write_bytes(counts_gz.read_bytes())

    raw = pd.read_csv(counts_gz, index_col=0)
    missing = [c for c in KO + WT if c not in raw.columns]
    if missing:
        raise SystemExit(f"missing columns: {missing}; have {list(raw.columns)}")

    ens2sym = build_symbol_map(raw.index.astype(str).tolist())
    symbols = pd.Series(raw.index.astype(str).map(ens2sym), index=raw.index)
    n_mapped = int(symbols.notna().sum())
    print(f"mapped {n_mapped}/{len(raw)} Ensembl IDs to symbols", flush=True)

    # Collapse duplicate symbols by mean (rare)
    expr = raw.copy()
    expr["symbol"] = symbols
    expr = expr.dropna(subset=["symbol"])
    expr = expr.groupby("symbol", sort=False).mean(numeric_only=True)

    samples = KO + WT
    group = pd.Series(["KO"] * len(KO) + ["WT"] * len(WT), index=samples)
    mat = expr[samples].astype(float)
    log2 = np.log2(mat + 1.0)

    # Expression filter for genome-wide tests
    keep = mat.mean(axis=1) >= 1.0
    log2_f = log2.loc[keep]
    print(f"genes mean-count>=1: {int(keep.sum())}/{len(mat)}", flush=True)

    # Per-gene DE
    rows = []
    for gene, row in log2_f.iterrows():
        a = row[KO].to_numpy()
        b = row[WT].to_numpy()
        lfc = float(a.mean() - b.mean())
        t, p_t = welch_t(a, b)
        u, p_u = mwu(a, b)
        rows.append(
            {
                "gene": gene,
                "log2FC_KO_minus_WT": lfc,
                "mean_log2_KO": float(a.mean()),
                "mean_log2_WT": float(b.mean()),
                "mean_norm_KO": float(mat.loc[gene, KO].mean()),
                "mean_norm_WT": float(mat.loc[gene, WT].mean()),
                "welch_t": t,
                "welch_p": p_t,
                "mwu_U": u,
                "mwu_p": p_u,
                "cohens_d": cohens_d(a, b),
            }
        )
    de = pd.DataFrame(rows)
    de["welch_fdr"] = multipletests(de["welch_p"], method="fdr_bh")[1]
    de["mwu_fdr"] = multipletests(de["mwu_p"], method="fdr_bh")[1]
    de = de.sort_values("welch_p")
    de.to_csv(TABLES / "de_all_genes.csv", index=False)

    # Focal + QC gene table (include genes even if filtered, if present)
    want = []
    for g in QC_GENES + FOCAL + USER_CORE6 + CXCL_IFN + T_CYT:
        if g not in want:
            want.append(g)
    focal_rows = []
    for g in want:
        if g not in log2.index:
            focal_rows.append({"gene": g, "present": False})
            continue
        a = log2.loc[g, KO].to_numpy()
        b = log2.loc[g, WT].to_numpy()
        t, p_t = welch_t(a, b)
        u, p_u = mwu(a, b)
        rec = de[de["gene"] == g]
        fdr = float(rec["welch_fdr"].iloc[0]) if len(rec) else float("nan")
        focal_rows.append(
            {
                "gene": g,
                "present": True,
                "log2FC_KO_minus_WT": float(a.mean() - b.mean()),
                "mean_norm_KO": float(mat.loc[g, KO].mean()),
                "mean_norm_WT": float(mat.loc[g, WT].mean()),
                "welch_t": t,
                "welch_p": p_t,
                "welch_fdr_genomewide": fdr,
                "mwu_p": p_u,
                "cohens_d": cohens_d(a, b),
                "values_KO": ",".join(f"{x:.4f}" for x in a),
                "values_WT": ",".join(f"{x:.4f}" for x in b),
            }
        )
    focal = pd.DataFrame(focal_rows)
    focal.to_csv(TABLES / "focal_genes.csv", index=False)

    # Sample-level z-score means
    is_ko = group.loc[samples].eq("KO").to_numpy()
    set_rows = []
    sample_scores = {"sample": samples, "group": group.loc[samples].tolist()}
    for name in REPORT_ORDER:
        genes = [g for g in SETS[name] if g in log2_f.index]
        missing_g = [g for g in SETS[name] if g not in log2_f.index]
        if len(genes) < 2:
            set_rows.append(
                {
                    "set": name,
                    "n_in_set": len(SETS[name]),
                    "n_tested": len(genes),
                    "missing": ",".join(missing_g),
                    "note": "too few genes present",
                }
            )
            continue
        sub = log2_f.loc[genes]
        z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1, ddof=1).replace(0, np.nan), axis=0)
        score = z.mean(axis=0)
        sample_scores[name] = score.loc[samples].to_numpy()
        a = score.loc[KO].to_numpy()
        b = score.loc[WT].to_numpy()
        delta = float(a.mean() - b.mean())
        t, p_t = welch_t(a, b)
        _, p_u = mwu(a, b)
        perm = exact_perm_p(score.loc[samples].to_numpy(), is_ko, delta)

        # Competitive: set log2FC vs background
        set_lfc = de.set_index("gene").loc[genes, "log2FC_KO_minus_WT"].to_numpy()
        bg_lfc = de.loc[~de["gene"].isin(genes), "log2FC_KO_minus_WT"].to_numpy()
        u_comp, p_comp = mwu(set_lfc, bg_lfc)
        # one-sided greater (set up vs background)
        try:
            _, p_comp_up = stats.mannwhitneyu(set_lfc, bg_lfc, alternative="greater")
        except ValueError:
            p_comp_up = 1.0
        set_rows.append(
            {
                "set": name,
                "n_in_set": len(SETS[name]),
                "n_tested": len(genes),
                "missing": ",".join(missing_g),
                "median_lfc_set": float(np.median(set_lfc)),
                "median_lfc_bg": float(np.median(bg_lfc)),
                "delta_median_lfc": float(np.median(set_lfc) - np.median(bg_lfc)),
                "mean_lfc_set": float(np.mean(set_lfc)),
                "competitive_mwu_p_two": p_comp,
                "competitive_mwu_p_up": float(p_comp_up),
                "competitive_auc": auc_from_mwu(u_comp, len(set_lfc), len(bg_lfc)),
                "sample_score_delta_KO_minus_WT": delta,
                "sample_score_welch_t": t,
                "sample_score_welch_p": p_t,
                "sample_score_mwu_p": p_u,
                "sample_score_cohens_d": cohens_d(a, b),
                **perm,
                "genes_tested": ",".join(genes),
            }
        )
    set_df = pd.DataFrame(set_rows)
    set_df.to_csv(TABLES / "geneset_stats.csv", index=False)
    pd.DataFrame(sample_scores).to_csv(TABLES / "sample_set_scores.csv", index=False)

    # Per-sample raw values for key genes
    key = [g for g in ["Tacstd2", "Cldn4", "Cldn7", "Cxcl9", "Cxcl10", "Isg15", "Stat1", "B2m", "H2-K1", "Cd8a", "Ifng", "Gzmb"] if g in log2.index]
    key_tab = log2.loc[key, samples].T
    key_tab.insert(0, "group", group.loc[samples].values)
    key_tab.to_csv(TABLES / "key_genes_log2.csv")

    # Figures
    def box(ax, gene):
        if gene not in log2.index:
            ax.set_title(f"{gene} (absent)")
            ax.axis("off")
            return
        data = [log2.loc[gene, WT].to_numpy(), log2.loc[gene, KO].to_numpy()]
        bp = ax.boxplot(data, tick_labels=["WT", "KO"], patch_artist=True, widths=0.55)
        for patch, c in zip(bp["boxes"], ["#4C78A8", "#F58518"]):
            patch.set_facecolor(c)
            patch.set_alpha(0.7)
        for i, vals in enumerate(data, start=1):
            x = np.random.default_rng(0).normal(i, 0.04, size=len(vals))
            ax.scatter(x, vals, c="black", s=18, zorder=3)
        rec = focal[focal["gene"] == gene]
        if len(rec) and bool(rec["present"].iloc[0]):
            lfc = rec["log2FC_KO_minus_WT"].iloc[0]
            p = rec["welch_p"].iloc[0]
            ax.set_title(f"{gene}\nlog2FC={lfc:+.2f}  p={p:.3g}")
        else:
            ax.set_title(gene)
        ax.set_ylabel("log2(norm+1)")

    fig, axes = plt.subplots(2, 3, figsize=(10.5, 7.2))
    for ax, g in zip(axes.ravel(), ["Tacstd2", "Cldn4", "Cldn7", "Cxcl9", "Isg15", "Cd8a"]):
        box(ax, g)
    fig.suptitle("GSE334497 4T1 Trop2 KO vs WT (bulk tumor, n=5+5)", fontsize=12)
    fig.tight_layout()
    fig.savefig(FIGS / "focal_boxplots.png", dpi=160)
    plt.close(fig)

    # Heatmap of IFN/APM + Cxcl + T_CYT
    heat_genes = []
    for lst in (FOCAL, USER_CORE6, MHC1_APM[:12], CXCL_IFN, T_CYT[:8]):
        for g in lst:
            if g in log2.index and g not in heat_genes:
                heat_genes.append(g)
    if "Tacstd2" in log2.index:
        heat_genes = ["Tacstd2"] + [g for g in heat_genes if g != "Tacstd2"]
    H = log2.loc[heat_genes, WT + KO]
    Hz = H.sub(H.mean(axis=1), axis=0).div(H.std(axis=1, ddof=1).replace(0, np.nan), axis=0)
    fig, ax = plt.subplots(figsize=(8.2, max(4.5, 0.28 * len(heat_genes) + 1.4)))
    im = ax.imshow(Hz.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    ax.set_xticks(range(len(Hz.columns)))
    ax.set_xticklabels(Hz.columns, rotation=75, ha="right", fontsize=7)
    ax.set_yticks(range(len(heat_genes)))
    ax.set_yticklabels(heat_genes, fontsize=8)
    ax.axvline(4.5, color="black", lw=1)
    ax.set_title("Row-z log2(norm+1)  |  WT (left) vs Trop2 KO (right)")
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="z")
    fig.tight_layout()
    fig.savefig(FIGS / "ifn_apm_heatmap.png", dpi=160)
    plt.close(fig)

    # Set-score boxes
    fig, axes = plt.subplots(2, 4, figsize=(12, 6.4))
    axes = axes.ravel()
    for i, name in enumerate(REPORT_ORDER):
        ax = axes[i]
        if name not in sample_scores:
            ax.axis("off")
            continue
        sc = np.array(sample_scores[name])
        data = [sc[~is_ko], sc[is_ko]]
        bp = ax.boxplot(data, tick_labels=["WT", "KO"], patch_artist=True, widths=0.55)
        for patch, c in zip(bp["boxes"], ["#4C78A8", "#F58518"]):
            patch.set_facecolor(c)
            patch.set_alpha(0.7)
        rec = set_df[set_df["set"] == name].iloc[0]
        ax.set_title(
            f"{name} (n={int(rec['n_tested'])})\n"
            f"Δ={rec['sample_score_delta_KO_minus_WT']:+.2f}  "
            f"perm p={rec['perm_p_two']:.3f}"
        )
        ax.set_ylabel("mean z-score")
    fig.suptitle("Pre-specified set scores, GSE334497 Trop2 KO − WT")
    fig.tight_layout()
    fig.savefig(FIGS / "set_score_boxplots.png", dpi=160)
    plt.close(fig)

    # Manifest / machine summary
    def gene_rec(g: str) -> dict:
        r = focal[focal["gene"] == g]
        if not len(r):
            return {"gene": g, "present": False}
        return {k: (None if (isinstance(v, float) and math.isnan(v)) else v) for k, v in r.iloc[0].to_dict().items()}

    summary = {
        "accession": "GSE334497",
        "claim": "C4 analog (not CLDN4 KD): 4T1 Trop2 KO vs WT bulk tumors",
        "contrast": "KO minus WT",
        "n_KO": len(KO),
        "n_WT": len(WT),
        "KO": KO,
        "WT": WT,
        "n_genes_matrix": int(len(raw)),
        "n_genes_mapped": n_mapped,
        "n_genes_tested": int(keep.sum()),
        "counts_md5": md5(counts_gz),
        "counts_url": COUNTS_URL,
        "perturbation_Tacstd2": gene_rec("Tacstd2"),
        "Cldn4": gene_rec("Cldn4"),
        "Cldn7": gene_rec("Cldn7"),
        "Cxcl9": gene_rec("Cxcl9"),
        "sets": set_df.replace({np.nan: None}).to_dict(orient="records"),
        "honesty": {
            "is_cldn4_kd": False,
            "is_trop2_ko": True,
            "tissue": "mouse 4T1 mammary tumor, immunocompetent BALB/c, bulk frozen section",
            "paper_mediator": "Cldn7 (protein localization), not Cldn4",
            "same_paper_as_hypothesis": True,
            "tumor_intrinsic_IFN_identifiable": False,
        },
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    # Human-readable key lines
    print("=== FOCAL ===")
    print(focal[["gene", "present", "log2FC_KO_minus_WT", "welch_p", "welch_fdr_genomewide", "cohens_d"]].to_string(index=False))
    print("=== SETS ===")
    cols = [
        "set",
        "n_tested",
        "median_lfc_set",
        "competitive_mwu_p_two",
        "sample_score_delta_KO_minus_WT",
        "sample_score_welch_p",
        "perm_p_two",
    ]
    print(set_df[cols].to_string(index=False))


if __name__ == "__main__":
    main()
