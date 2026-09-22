#!/usr/bin/env python3
"""Tacstd2-high vs Tacstd2-low malignant/epithelial DEG → tight-junction enrichment.

Public integrate cohorts only: GSE154977, GSE180963, GSE165641.
Private 8 KL matrices are not read and are not merged.

Unit for inference is the mouse. Within each mouse, locked-gate epithelial
cells are split into Tacstd2 count>0 (high) vs Tacstd2 count==0 (low). Cell-level
Mann–Whitney on a frozen TJ module score is the primary within-mouse test.
Pseudobulk mean-log1p deltas are then ranked and tested for TJ set enrichment
by a one-sided hypergeometric on genes with mean delta > 0 and by a simple
enrichment score on the ranked list (mean rank of TJ genes vs all genes).

Tacstd2 is never used to call epithelium.
"""

from __future__ import annotations

import argparse
import gzip
import json
import platform
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import __version__ as scipy_version
from scipy.io import mmread
from scipy.stats import hypergeom, mannwhitneyu, rankdata

TJ_EPITHELIAL = [
    "Cldn1", "Cldn3", "Cldn4", "Cldn7", "Ocln", "Marveld2", "Marveld3",
    "Tjp1", "Tjp2", "Tjp3", "F11r", "Jam2", "Jam3", "Cgn", "Cgnl1", "Crb3",
    "Ildr1", "Lsr",
]
TJ_TISMO = ["Cldn3", "Cldn4", "Cldn6", "Cldn7", "Cdh1", "F11r", "Ocln"]
CLDN4_TJ_EDGE = [
    "Cldn4", "Cldn1", "Cldn7", "Cgnl1", "Marveld2", "Marveld3",
    "Tjp1", "Tjp2", "Ildr1", "Cldn3", "Ocln",
]
SETS = {
    "TJ_EPITHELIAL": TJ_EPITHELIAL,
    "TJ_TISMO": TJ_TISMO,
    "CLDN4_TJ_EDGE": CLDN4_TJ_EDGE,
}

EPI_CORE = ["Epcam", "Cdh1", "Krt8", "Krt18", "Krt19"]
STRUCT = ["Cdh1", "Krt8", "Krt18", "Krt19", "Cldn18"]
CALLER = sorted(set(EPI_CORE + STRUCT + ["Epcam", "Ptprc", "Tacstd2", "Sftpc"]))


def _read_symbols(path: Path) -> list[str]:
    opener = gzip.open if str(path).endswith(".gz") else open
    symbols = []
    with opener(path, "rt") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[1] not in {"", "Gene Expression"} and not parts[1].startswith("ENSMUS"):
                if parts[0].startswith("ENS") or parts[0].startswith("ENSMUS"):
                    symbols.append(parts[1])
                else:
                    symbols.append(parts[0])
            else:
                symbols.append(parts[0])
    return symbols


def _load_mtx(path: Path):
    if str(path).endswith(".gz"):
        with gzip.open(path, "rb") as handle:
            mat = mmread(handle)
    else:
        mat = mmread(path)
    return mat.tocsr()


def _qc(mat, symbols):
    ncount = np.asarray(mat.sum(axis=0)).ravel().astype(np.float64)
    nfeature = np.asarray((mat > 0).sum(axis=0)).ravel().astype(np.int32)
    mt_idx = [i for i, s in enumerate(symbols) if s.startswith("mt-")]
    if mt_idx:
        mt = np.asarray(mat[mt_idx].sum(axis=0)).ravel().astype(np.float64)
    else:
        mt = np.zeros(mat.shape[1], dtype=np.float64)
    pct_mt = np.where(ncount > 0, 100.0 * mt / ncount, 100.0)
    keep = (nfeature >= 200) & (ncount >= 500) & (pct_mt < 25)
    return keep, ncount


def _first_index(symbols):
    out = {}
    for i, s in enumerate(symbols):
        if s not in out:
            out[s] = i
    return out


def locked_epi(mat, symbols, digest: str, keep: np.ndarray) -> np.ndarray:
    idx = _first_index(symbols)
    def pos(genes):
        rows = [idx[g] for g in genes if g in idx]
        if not rows:
            return np.zeros(int(keep.sum()), dtype=bool)
        return np.asarray(mat[rows][:, keep].tocsr().sum(axis=0)).ravel() > 0

    epcam = pos(["Epcam"])
    struct = np.zeros(int(keep.sum()), dtype=bool)
    for g in STRUCT:
        struct |= pos([g])
    ptprc = pos(["Ptprc"])
    not_p = ~ptprc
    if digest == "AT2_lineage_FACS":
        return not_p & (epcam | struct)
    return (epcam & struct & not_p)


def library_specs(root: Path):
    specs = []
    # GSE154977: four mice in one h5
    specs.append(("GSE154977", None, "KP", "AT2_lineage_FACS", "h5"))
    for folder, mouse, geno in (("K", "GSE180963_K", "K"), ("KL", "GSE180963_KL", "KL")):
        specs.append(("GSE180963", mouse, geno, "mixed", root / "GSE180963" / folder))
    for hint, mouse in (("KL1", "GSE165641_KL1"), ("KL2", "GSE165641_KL2")):
        hits = [p for p in (root / "GSE165641").rglob("filtered_feature_bc_matrix") if hint in str(p)]
        if len(hits) != 1:
            raise RuntimeError(f"expected one matrix for {hint}, found {hits}")
        specs.append(("GSE165641", mouse, "KL", "mixed", hits[0]))
    return specs


def load_10x(directory: Path):
    mtx = next(directory.glob("matrix.mtx*"))
    genes = directory / "features.tsv.gz"
    if not genes.exists():
        genes = directory / "features.tsv"
    if not genes.exists():
        genes = directory / "genes.tsv.gz"
    if not genes.exists():
        genes = directory / "genes.tsv"
    symbols = _read_symbols(genes)
    mat = _load_mtx(mtx)
    if mat.shape[0] != len(symbols):
        if mat.shape[1] == len(symbols):
            mat = mat.T.tocsr()
        else:
            raise RuntimeError(f"shape {mat.shape} vs genes {len(symbols)}")
    return symbols, mat


def load_gse154977(root: Path):
    d = root / "GSE154977"
    genes = pd.read_csv(d / "GSE154977_mmLung10x_cis_geneTable.csv.gz")
    smp = pd.read_csv(d / "GSE154977_mmLung10x_cis_smpTable.csv.gz")
    symbols = genes["geneID"].tolist()
    with h5py.File(d / "GSE154977_mmLung10x_cis_dSp_rawCount.h5", "r") as handle:
        ii = np.array(handle["i"]).ravel()
        jj = np.array(handle["j"]).ravel()
        vv = np.array(handle["v"]).ravel()
    if ii.min() == 0 or jj.min() == 0:
        ii = ii + 1
        jj = jj + 1
    mat = __import__("scipy").sparse.coo_matrix(
        (vv, (ii.astype(np.int64) - 1, jj.astype(np.int64) - 1)),
        shape=(len(symbols), len(smp)),
    ).tocsr()
    return symbols, mat, smp


def score_module(lognorm: np.ndarray, idx: dict[str, int], genes: list[str]) -> np.ndarray:
    rows = [idx[g] for g in genes if g in idx]
    if not rows:
        return np.full(lognorm.shape[1], np.nan)
    return lognorm[rows].mean(axis=0)


def hypergeo_up(ranked_genes: list[str], set_genes: set[str], top_n: int) -> dict:
    top = set(ranked_genes[:top_n])
    k = len(top & set_genes)
    M = len(ranked_genes)
    n = len(set_genes & set(ranked_genes))
    N = top_n
    if n == 0 or N == 0:
        return dict(k=0, n_set=n, top_n=N, universe=M, p=1.0, enrichment=np.nan)
    # P(X >= k)
    p = float(hypergeom.sf(k - 1, M, n, N))
    expected = N * n / M
    return dict(k=k, n_set=n, top_n=N, universe=M, p=p, enrichment=k / expected if expected else np.nan)


def enrichment_score(deltas: pd.Series, set_genes: set[str]) -> dict:
    """Simple score: difference between mean rank of set genes and mean rank of all."""
    present = [g for g in deltas.index if g in set_genes]
    if len(present) < 3:
        return dict(n=len(present), mean_delta=np.nan, mean_rank=np.nan, es=np.nan, p_mw=np.nan)
    ranks = rankdata(-deltas.to_numpy(), method="average")  # high delta = rank 1
    rank_map = dict(zip(deltas.index.tolist(), ranks))
    set_ranks = np.array([rank_map[g] for g in present], dtype=float)
    other_genes = [g for g in deltas.index.tolist() if g not in set_genes]
    other_ranks = np.array([rank_map[g] for g in other_genes], dtype=float)
    mean_delta = float(deltas.loc[present].mean())
    # One-sided: set genes have higher delta than other genes
    p = float(mannwhitneyu(
        deltas.loc[present].to_numpy(),
        deltas.loc[other_genes].to_numpy(),
        alternative="greater",
    ).pvalue)
    es = float(other_ranks.mean() - set_ranks.mean())  # positive if set ranks better (lower rank number)
    return dict(
        n=len(present),
        mean_delta=mean_delta,
        mean_rank=float(set_ranks.mean()),
        es=es,
        p_mw=p,
        genes=",".join(present),
    )


def process_mouse(symbols, mat, keep, digest, mouse, dataset, genotype, min_high=20, min_low=20):
    idx = _first_index(symbols)
    if "Tacstd2" not in idx:
        return None, None, None
    epi_mask = locked_epi(mat, symbols, digest, keep)
    ncount = np.asarray(mat[:, keep].sum(axis=0)).ravel().astype(np.float64)
    # Restrict to epi cells
    epi_cells = np.where(epi_mask)[0]
    if len(epi_cells) == 0:
        meta = dict(
            mouse=mouse, dataset=dataset, genotype=genotype, digest=digest,
            n_epi=0, n_high=0, n_low=0, status="no_epi",
        )
        return meta, None, None
    tac = np.asarray(mat[idx["Tacstd2"], keep][:, epi_cells].todense()).ravel()
    high = tac > 0
    low = tac == 0
    n_high = int(high.sum())
    n_low = int(low.sum())
    if n_high < min_high or n_low < min_low or len(epi_cells) < (min_high + min_low):
        meta = dict(
            mouse=mouse, dataset=dataset, genotype=genotype, digest=digest,
            n_epi=int(len(epi_cells)), n_high=n_high, n_low=n_low,
            tacstd2_pct=100.0 * n_high / len(epi_cells),
            status="below_floor",
        )
        return meta, None, None

    # log1p CP10k for epi cells (dense only for needed rows would be better; use full sparse slice)
    sub = mat[:, keep][:, epi_cells].tocsr()
    # Convert to lognorm gene x epi
    # Memory: genes x n_epi float32 — for ~20k x 4k = 320MB OK
    dens = sub.toarray().astype(np.float32)
    denom = ncount[epi_cells]
    with np.errstate(divide="ignore", invalid="ignore"):
        logn = np.log1p(dens / denom * 10000.0)

    # Module scores
    score_rows = []
    for name, genes in SETS.items():
        rows = [idx[g] for g in genes if g in idx]
        if len(rows) < 3:
            continue
        sc = logn[rows].mean(axis=0)
        a = sc[high]
        b = sc[low]
        delta = float(a.mean() - b.mean())
        p = float(mannwhitneyu(a, b, alternative="greater").pvalue)
        score_rows.append(dict(
            mouse=mouse, dataset=dataset, genotype=genotype, set=name,
            n_genes=len(rows), n_high=n_high, n_low=n_low,
            mean_high=float(a.mean()), mean_low=float(b.mean()),
            delta=delta, p_greater=p,
        ))

    # Pseudobulk gene deltas (high - low)
    mean_h = logn[:, high].mean(axis=1)
    mean_l = logn[:, low].mean(axis=1)
    delta = mean_h - mean_l
    # Drop all-zero genes in this mouse epi
    expressed = (dens.sum(axis=1) > 0)
    genes = np.array(symbols)[expressed]
    deltas = pd.Series(delta[expressed], index=genes)
    # Collapse duplicate symbols by mean
    if deltas.index.duplicated().any():
        deltas = deltas.groupby(level=0).mean()
    meta = dict(
        mouse=mouse, dataset=dataset, genotype=genotype, digest=digest,
        n_epi=int(len(epi_cells)), n_high=n_high, n_low=n_low, status="ok",
        tacstd2_pct=100.0 * n_high / len(epi_cells),
    )
    return meta, pd.DataFrame(score_rows), deltas


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    out = args.out
    tab = out / "tables"
    fig = out / "figures"
    tab.mkdir(parents=True, exist_ok=True)
    fig.mkdir(parents=True, exist_ok=True)

    meta_rows = []
    score_parts = []
    delta_parts = []

    # GSE154977 mice
    print("DEG GSE154977", flush=True)
    symbols, mat, smp = load_gse154977(args.data)
    keep, _ = _qc(mat, symbols)
    libraries = []
    for sid in smp.loc[keep, "sampleID"]:
        library = sid.split("_id-")[0]
        mouse_raw = library[:-3] if library.endswith("_PT") else library
        libraries.append(f"GSE154977_{mouse_raw}")
    library_arr = np.array(libraries)
    for mouse in sorted(set(libraries)):
        cell_mask = (library_arr == mouse)
        # rebuild keep relative mask: keep is over all cells; cell_mask over kept cells
        # Map: positions in keep.true cells
        keep_idx = np.where(keep)[0]
        mouse_keep = np.zeros(mat.shape[1], dtype=bool)
        mouse_keep[keep_idx[cell_mask]] = True
        # Process using only this mouse's cells as "keep"
        meta, scores, deltas = process_mouse(
            symbols, mat, mouse_keep, "AT2_lineage_FACS", mouse, "GSE154977", "KP",
            min_high=30, min_low=30,
        )
        if meta:
            meta_rows.append(meta)
        if scores is not None:
            score_parts.append(scores)
            deltas.name = mouse
            delta_parts.append(deltas)
            print(f"  {mouse} high={meta['n_high']} low={meta['n_low']}", flush=True)

    # 10x libraries
    for dataset, mouse, geno, digest, path in library_specs(args.data):
        if dataset == "GSE154977":
            continue
        print(f"DEG {mouse}", flush=True)
        symbols, mat = load_10x(path)
        keep, _ = _qc(mat, symbols)
        meta, scores, deltas = process_mouse(
            symbols, mat, keep, digest, mouse, dataset, geno, min_high=20, min_low=20,
        )
        if meta:
            meta_rows.append(meta)
        if scores is not None:
            score_parts.append(scores)
            deltas.name = mouse
            delta_parts.append(deltas)
            print(f"  {mouse} high={meta['n_high']} low={meta['n_low']}", flush=True)
        else:
            print(f"  {mouse} skipped {meta}", flush=True)

    meta_df = pd.DataFrame(meta_rows)
    meta_df.to_csv(tab / "deg_mouse_inventory.tsv", sep="\t", index=False)

    if not score_parts:
        raise SystemExit("no mouse passed Tacstd2-high/low floor")

    scores = pd.concat(score_parts, ignore_index=True)
    scores.to_csv(tab / "tj_module_high_vs_low.tsv", sep="\t", index=False)

    # Across-mouse consistency for module scores
    cons = []
    for name, sub in scores.groupby("set"):
        n = len(sub)
        n_pos = int((sub["delta"] > 0).sum())
        n_sig = int(((sub["delta"] > 0) & (sub["p_greater"] < 0.05)).sum())
        # Exact binomial one-sided: P(X>=n_pos) under p=0.5 for direction
        from scipy.stats import binom
        p_dir = float(binom.sf(n_pos - 1, n, 0.5)) if n else np.nan
        cons.append(dict(
            set=name, n_mice=n, n_delta_pos=n_pos, n_p_lt_05_and_pos=n_sig,
            mean_delta=float(sub["delta"].mean()),
            median_delta=float(sub["delta"].median()),
            binomial_p_direction=p_dir,
        ))
    cons_df = pd.DataFrame(cons).sort_values("set")
    cons_df.to_csv(tab / "tj_module_consistency.tsv", sep="\t", index=False)

    # Combined gene deltas: mean across mice that have the gene
    delta_mat = pd.concat(delta_parts, axis=1)
    mean_delta = delta_mat.mean(axis=1, skipna=True)
    n_mice_gene = delta_mat.notna().sum(axis=1)
    mean_delta = mean_delta[n_mice_gene >= 2]
    mean_delta = mean_delta.sort_values(ascending=False)
    mean_delta.to_csv(tab / "gene_mean_delta_tacstd2_high_minus_low.tsv", sep="\t", header=["mean_delta"])

    enrich_rows = []
    for name, genes in SETS.items():
        es = enrichment_score(mean_delta, set(genes))
        ranked = mean_delta.index.tolist()
        for top_n in (50, 100, 200, 500):
            hg = hypergeo_up(ranked, set(genes), top_n)
            enrich_rows.append(dict(set=name, method="hypergeo_top", **hg, **{f"es_{k}": v for k, v in es.items()}))
        enrich_rows.append(dict(
            set=name, method="rank_mw", top_n=np.nan, k=np.nan, n_set=es["n"],
            universe=len(mean_delta), p=es["p_mw"], enrichment=es["es"],
            mean_delta=es["mean_delta"], mean_rank=es["mean_rank"],
        ))
    enrich = pd.DataFrame(enrich_rows)
    enrich.to_csv(tab / "tj_enrichment.tsv", sep="\t", index=False)

    # Top positive genes that are TJ members
    tj_all = sorted(set(TJ_EPITHELIAL + TJ_TISMO + CLDN4_TJ_EDGE))
    tj_hits = mean_delta.loc[[g for g in mean_delta.index if g in tj_all]].sort_values(ascending=False)
    tj_hits.to_csv(tab / "tj_gene_deltas.tsv", sep="\t", header=["mean_delta"])

    # Figure: module deltas per mouse
    fig1, ax = plt.subplots(figsize=(7.2, 4.2))
    sets_order = ["TJ_TISMO", "TJ_EPITHELIAL", "CLDN4_TJ_EDGE"]
    xpos = np.arange(len(sets_order))
    width = 0.1
    mice = sorted(scores["mouse"].unique())
    colors = plt.cm.tab10(np.linspace(0, 1, max(len(mice), 1)))
    for i, mouse in enumerate(mice):
        vals = []
        for s in sets_order:
            hit = scores[(scores.mouse == mouse) & (scores.set == s)]
            vals.append(float(hit["delta"].iloc[0]) if len(hit) else np.nan)
        ax.bar(xpos + (i - len(mice) / 2) * width, vals, width=width, label=mouse.replace("GSE", ""), color=colors[i])
    ax.axhline(0, color="0.4", lw=0.8)
    ax.set_xticks(xpos)
    ax.set_xticklabels(sets_order, rotation=15, ha="right")
    ax.set_ylabel("Tacstd2-high − low TJ module (log1p)")
    ax.set_title("Within-mouse TJ module: Tacstd2-high vs low epithelium")
    ax.legend(fontsize=7, ncol=2)
    fig1.tight_layout()
    fig1.savefig(fig / "tj_module_deltas.png", dpi=140)
    fig1.savefig(fig / "tj_module_deltas.pdf")
    plt.close(fig1)

    # Primary headline: TJ_TISMO consistency
    primary = cons_df[cons_df["set"] == "TJ_TISMO"].iloc[0]
    best_hg = enrich[(enrich.set == "TJ_TISMO") & (enrich.method == "hypergeo_top")].sort_values("p").iloc[0]
    rank_mw = enrich[(enrich.set == "TJ_TISMO") & (enrich.method == "rank_mw")].iloc[0]

    finding = f"""# Tacstd2-high epithelial DEG → tight junction (public integrate mice)

Public processed counts only: GSE154977, GSE180963, GSE165641. Private 8 KL matrices were not read and were not merged. Epithelium is the locked gate from the integrate (Epcam∧structural∧Ptprc−, or AT2-FACS Ptprc− with Epcam∨structural). **Tacstd2 is not an epithelial caller.**

Within each mouse, epithelial cells with Tacstd2 count > 0 are compared to epithelial cells with Tacstd2 count = 0. Mice below the high/low cell floor are inventoried and dropped from DEG.

## Module scores (primary)

Frozen TJ lists from the GSE137244 transfer signatures (`gene_sets/signatures.gmt`). Cell-level mean log1p module score; one-sided Mann–Whitney that Tacstd2-high > Tacstd2-low; binomial test on the sign of the mouse-level delta.

| set | n mice | Δ>0 | Δ>0 and p<0.05 | mean Δ | binomial p (direction) |
|---|---:|---:|---:|---:|---:|
"""
    for _, r in cons_df.iterrows():
        finding += (
            f"| {r['set']} | {int(r['n_mice'])} | {int(r['n_delta_pos'])} | "
            f"{int(r['n_p_lt_05_and_pos'])} | {r['mean_delta']:+.4f} | {r['binomial_p_direction']:.4g} |\n"
        )
    finding += f"""
**Headline (TJ_TISMO 7-gene):** {int(primary['n_delta_pos'])}/{int(primary['n_mice'])} mice have Tacstd2-high epithelium above Tacstd2-low (mean Δ={primary['mean_delta']:+.4f}; one-sided binomial p={primary['binomial_p_direction']:.4g}). {int(primary['n_p_lt_05_and_pos'])}/{int(primary['n_mice'])} also have within-mouse MW p<0.05.

## Ranked gene enrichment (secondary)

Mean across-mouse pseudobulk delta (Tacstd2-high − low) for genes present in ≥2 mice. TJ_TISMO rank Mann–Whitney (set deltas greater than background) p={rank_mw['p']:.4g}. Strongest hypergeometric among top-N cutoffs: top {int(best_hg['top_n'])}, k={int(best_hg['k'])}/{int(best_hg['n_set'])}, enrichment={best_hg['enrichment']:.2f}, p={best_hg['p']:.4g}.

Cell-level p-values inside a mouse are descriptive for calling the module direction; the mouse is the inferential unit for the binomial and for the across-mouse gene ranking.

Mice scored: {', '.join(meta_df.loc[meta_df.status=='ok','mouse'].tolist())}.
Private 8 KL mice used: 0.
"""
    (out / "FINDING_DEG_TJ.md").write_text(finding)

    summary = {
        "private_8kl_mice": 0,
        "merged_with_private_8kl": False,
        "n_mice_scored": int((meta_df.status == "ok").sum()),
        "n_mice_inventoried": int(len(meta_df)),
        "tj_tismo": cons_df[cons_df.set == "TJ_TISMO"].iloc[0].to_dict(),
        "python": platform.python_version(),
        "scipy": scipy_version,
    }
    (tab / "deg_tj_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print("wrote", out / "FINDING_DEG_TJ.md", flush=True)


if __name__ == "__main__":
    main()
