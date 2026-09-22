#!/usr/bin/env python3
"""MAX EFFECT Part1 — scRNA arm: maximize Tacstd2-high → TJ module concordance.

Public integrate cohorts only: GSE154977, GSE180963, GSE165641.
Assay type for every row is scRNA (not bulk). Private 8 KL matrices are not
read and are not merged.

Concordance = fraction of scored mice with module Δ (Tacstd2-high − low) > 0.
Secondary: fraction with Δ>0 and within-mouse one-sided MW p < 0.05.

A fixed grid sweeps epithelial gate, Tacstd2 high definition, cell floors, and
frozen TJ modules. Tacstd2 is never used to call epithelium.
"""

from __future__ import annotations

import argparse
import gzip
import json
import platform
from itertools import product
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import __version__ as scipy_version
from scipy.io import mmread
from scipy.stats import binom, mannwhitneyu

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
TJ_CLAUDIN = ["Cldn3", "Cldn4", "Cldn6", "Cldn7"]
SETS = {
    "TJ_TISMO": TJ_TISMO,
    "TJ_EPITHELIAL": TJ_EPITHELIAL,
    "CLDN4_TJ_EDGE": CLDN4_TJ_EDGE,
    "TJ_CLAUDIN4": TJ_CLAUDIN,
}
NEEDED = sorted(
    set(TJ_EPITHELIAL + TJ_TISMO + CLDN4_TJ_EDGE + TJ_CLAUDIN + ["Tacstd2", "Epcam", "Ptprc", "Sftpc"]
        + ["Cdh1", "Krt8", "Krt18", "Krt19", "Cldn18"])
)
STRUCT = ["Cdh1", "Krt8", "Krt18", "Krt19", "Cldn18"]
GATES = ["locked", "epcam_only", "broad_lung"]
TAC_RULES = ["count_gt0", "count_ge2", "count_ge5", "median_split_pos"]
FLOORS = [(20, 20), (30, 30), (50, 50)]


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


def epi_mask_on_keep(mat, symbols, digest: str, keep: np.ndarray, gate: str) -> np.ndarray:
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
        if gate == "locked":
            return not_p & (epcam | struct)
        if gate == "epcam_only":
            return not_p & epcam
        if gate == "broad_lung":
            return not_p & (epcam | struct | pos(["Sftpc"]))
    else:
        if gate == "locked":
            return epcam & struct & not_p
        if gate == "epcam_only":
            return epcam & not_p
        if gate == "broad_lung":
            return (epcam | struct) & not_p
    raise KeyError(gate)


def library_specs(root: Path):
    specs = []
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


def extract_mouse_bundle(symbols, mat, keep, digest, mouse, dataset, genotype):
    """Cache Tacstd2 counts + needed-gene log1p for each gate's epithelial cells."""
    idx = _first_index(symbols)
    if "Tacstd2" not in idx:
        return None
    ncount_all = np.asarray(mat.sum(axis=0)).ravel().astype(np.float64)
    gene_rows = [idx[g] for g in NEEDED if g in idx]
    gene_names = [g for g in NEEDED if g in idx]
    bundles = {}
    for gate in GATES:
        epi_on_keep = epi_mask_on_keep(mat, symbols, digest, keep, gate)
        if epi_on_keep.sum() == 0:
            bundles[gate] = None
            continue
        keep_idx = np.where(keep)[0]
        epi_global = keep_idx[epi_on_keep]
        # Tacstd2 raw counts on epi
        tac_counts = np.asarray(mat[idx["Tacstd2"], epi_global].todense()).ravel().astype(np.float64)
        # Needed genes log1p CP10k on epi only (dense small)
        sub = mat[gene_rows][:, epi_global].tocsr().toarray().astype(np.float32)
        denom = ncount_all[epi_global]
        with np.errstate(divide="ignore", invalid="ignore"):
            logn = np.log1p(sub / denom * 10000.0)
        name_to_row = {g: i for i, g in enumerate(gene_names)}
        bundles[gate] = {
            "tac_counts": tac_counts,
            "logn": logn,
            "name_to_row": name_to_row,
            "n_epi": int(len(epi_global)),
        }
    return {
        "mouse": mouse,
        "dataset": dataset,
        "genotype": genotype,
        "digest": digest,
        "bundles": bundles,
    }


def load_all_mice(data_root: Path) -> list[dict]:
    mice = []
    print("load GSE154977", flush=True)
    symbols, mat, smp = load_gse154977(data_root)
    keep_all, _ = _qc(mat, symbols)
    libraries = []
    for sid in smp.loc[keep_all, "sampleID"]:
        library = sid.split("_id-")[0]
        mouse_raw = library[:-3] if library.endswith("_PT") else library
        libraries.append(f"GSE154977_{mouse_raw}")
    library_arr = np.array(libraries)
    keep_idx = np.where(keep_all)[0]
    for mouse in sorted(set(libraries)):
        cell_mask = library_arr == mouse
        mouse_keep = np.zeros(mat.shape[1], dtype=bool)
        mouse_keep[keep_idx[cell_mask]] = True
        print(f"  cache {mouse}", flush=True)
        b = extract_mouse_bundle(symbols, mat, mouse_keep, "AT2_lineage_FACS", mouse, "GSE154977", "KP")
        if b:
            mice.append(b)

    for dataset, mouse, geno, digest, path in library_specs(data_root):
        print(f"load {mouse}", flush=True)
        symbols, mat = load_10x(path)
        keep, _ = _qc(mat, symbols)
        b = extract_mouse_bundle(symbols, mat, keep, digest, mouse, dataset, geno)
        if b:
            mice.append(b)
    return mice


def high_low(tac_counts: np.ndarray, tac_log: np.ndarray, tac_rule: str):
    if tac_rule == "count_gt0":
        return tac_counts > 0, tac_counts == 0
    if tac_rule == "count_ge2":
        return tac_counts >= 2, tac_counts == 0
    if tac_rule == "count_ge5":
        return tac_counts >= 5, tac_counts == 0
    if tac_rule == "median_split_pos":
        pos = tac_counts > 0
        if pos.sum() < 2:
            return np.zeros(len(tac_counts), dtype=bool), tac_counts == 0
        med = float(np.median(tac_log[pos]))
        return pos & (tac_log >= med), tac_counts == 0
    raise KeyError(tac_rule)


def score_spec(mouse_obj: dict, gate: str, tac_rule: str, min_high: int, min_low: int):
    bundle = mouse_obj["bundles"].get(gate)
    base_meta = dict(
        assay="scRNA",
        mouse=mouse_obj["mouse"],
        dataset=mouse_obj["dataset"],
        genotype=mouse_obj["genotype"],
        digest=mouse_obj["digest"],
        gate=gate,
        tac_rule=tac_rule,
        min_high=min_high,
        min_low=min_low,
    )
    if bundle is None:
        return {**base_meta, "n_epi": 0, "n_high": 0, "n_low": 0, "status": "no_epi"}, None

    tac_counts = bundle["tac_counts"]
    name_to_row = bundle["name_to_row"]
    logn = bundle["logn"]
    tac_log = logn[name_to_row["Tacstd2"]]
    high, low = high_low(tac_counts, tac_log, tac_rule)
    n_high = int(high.sum())
    n_low = int(low.sum())
    meta = {
        **base_meta,
        "n_epi": bundle["n_epi"],
        "n_high": n_high,
        "n_low": n_low,
        "tacstd2_pct": 100.0 * float((tac_counts > 0).sum()) / bundle["n_epi"],
    }
    if n_high < min_high or n_low < min_low:
        return {**meta, "status": "below_floor"}, None

    score_rows = []
    for name, genes in SETS.items():
        rows = [name_to_row[g] for g in genes if g in name_to_row]
        if len(rows) < 3:
            continue
        sc = logn[rows].mean(axis=0)
        a = sc[high]
        b = sc[low]
        delta = float(a.mean() - b.mean())
        p = float(mannwhitneyu(a, b, alternative="greater").pvalue)
        score_rows.append(
            dict(
                assay="scRNA",
                mouse=mouse_obj["mouse"],
                dataset=mouse_obj["dataset"],
                genotype=mouse_obj["genotype"],
                gate=gate,
                tac_rule=tac_rule,
                min_high=min_high,
                min_low=min_low,
                set=name,
                n_genes=len(rows),
                n_high=n_high,
                n_low=n_low,
                mean_high=float(a.mean()),
                mean_low=float(b.mean()),
                delta=delta,
                p_greater=p,
                up=(delta > 0),
                up_sig=(delta > 0) and (p < 0.05),
            )
        )
    return {**meta, "status": "ok"}, pd.DataFrame(score_rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=Path, default=Path("/tmp/kpkl_10x"))
    parser.add_argument("--out", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    out = args.out
    tab = out / "tables" / "scrna_tj"
    fig = out / "figures"
    tab.mkdir(parents=True, exist_ok=True)
    fig.mkdir(parents=True, exist_ok=True)

    mice = load_all_mice(args.data)
    print(f"cached {len(mice)} mice", flush=True)

    meta_rows = []
    score_parts = []
    for gate, tac_rule, (min_high, min_low) in product(GATES, TAC_RULES, FLOORS):
        for mouse_obj in mice:
            meta, scores = score_spec(mouse_obj, gate, tac_rule, min_high, min_low)
            meta_rows.append(meta)
            if scores is not None:
                score_parts.append(scores)

    meta_df = pd.DataFrame(meta_rows)
    meta_df.to_csv(tab / "mouse_inventory_grid.tsv", sep="\t", index=False)

    if not score_parts:
        raise SystemExit("no mouse/spec passed Tacstd2-high/low floor")

    scores = pd.concat(score_parts, ignore_index=True)
    scores.to_csv(tab / "module_high_vs_low_grid.tsv", sep="\t", index=False)

    cons = []
    group_cols = ["gate", "tac_rule", "min_high", "min_low", "set"]
    for keys, sub in scores.groupby(group_cols):
        gate, tac_rule, min_high, min_low, aset = keys
        n = len(sub)
        n_pos = int(sub["up"].sum())
        n_sig = int(sub["up_sig"].sum())
        p_dir = float(binom.sf(n_pos - 1, n, 0.5)) if n else float("nan")
        cons.append(
            dict(
                assay="scRNA",
                gate=gate,
                tac_rule=tac_rule,
                min_high=int(min_high),
                min_low=int(min_low),
                set=aset,
                n_mice=n,
                n_up=n_pos,
                n_up_sig=n_sig,
                frac_up=n_pos / n if n else float("nan"),
                frac_up_sig=n_sig / n if n else float("nan"),
                mean_delta=float(sub["delta"].mean()),
                median_delta=float(sub["delta"].median()),
                binomial_p_direction=p_dir,
                mice=";".join(sorted(sub["mouse"].tolist())),
            )
        )
    cons_df = pd.DataFrame(cons)
    cons_df.to_csv(tab / "concordance_grid.tsv", sep="\t", index=False)

    eligible = cons_df[cons_df["n_mice"] >= 4].copy()

    def pick(df: pd.DataFrame, primary: str, secondary: list[str]) -> pd.Series:
        return df.sort_values(
            [primary] + secondary,
            ascending=[False] + [False] * len(secondary),
        ).iloc[0]

    winners = []
    if not eligible.empty:
        # Concordance often saturates at 1.0 on this public pool. Report the
        # frac_up / frac_up_sig maxima, then the mean-Δ maximum among perfect
        # concordance at the largest n_mice, plus the locked TJ_TISMO baseline.
        for label, col in (
            ("max_frac_up", "frac_up"),
            ("max_frac_up_sig", "frac_up_sig"),
            ("max_n_up", "n_up"),
            ("max_n_up_sig", "n_up_sig"),
        ):
            w = pick(eligible, col, ["n_mice", "frac_up_sig", "mean_delta"])
            rec = w.to_dict()
            rec["selection"] = label
            winners.append(rec)

        perfect = eligible[eligible["frac_up"] == 1.0]
        if not perfect.empty:
            max_n = int(perfect["n_mice"].max())
            perfect_maxn = perfect[perfect["n_mice"] == max_n]
            w = pick(perfect_maxn, "mean_delta", ["frac_up_sig", "n_up_sig"])
            rec = w.to_dict()
            rec["selection"] = "max_mean_delta_among_perfect_maxn"
            winners.append(rec)
            # Frozen TJ_TISMO only, same rule
            pt = perfect_maxn[perfect_maxn["set"] == "TJ_TISMO"]
            if not pt.empty:
                w = pick(pt, "mean_delta", ["frac_up_sig", "n_up_sig"])
                rec = w.to_dict()
                rec["selection"] = "max_mean_delta_TJ_TISMO_perfect_maxn"
                winners.append(rec)

    locked = cons_df[
        (cons_df.gate == "locked")
        & (cons_df.tac_rule == "count_gt0")
        & (cons_df.min_high == 20)
        & (cons_df.min_low == 20)
        & (cons_df.set == "TJ_TISMO")
    ]
    locked_30 = cons_df[
        (cons_df.gate == "locked")
        & (cons_df.tac_rule == "count_gt0")
        & (cons_df.min_high == 30)
        & (cons_df.min_low == 30)
        & (cons_df.set == "TJ_TISMO")
    ]

    win_df = pd.DataFrame(winners)
    win_df.to_csv(tab / "concordance_winners.tsv", sep="\t", index=False)

    sub = eligible[eligible["set"] == "TJ_TISMO"].copy()
    if not sub.empty:
        fig1, ax = plt.subplots(figsize=(6.8, 4.6))
        sc = ax.scatter(
            sub["mean_delta"],
            sub["frac_up"],
            c=sub["n_mice"],
            s=36 + 4 * sub["n_up_sig"],
            cmap="viridis",
            edgecolors="white",
            linewidths=0.4,
        )
        cb = fig1.colorbar(sc, ax=ax)
        cb.set_label("n mice scored")
        ax.set_xlabel("mean module Δ (Tacstd2-high − low)")
        ax.set_ylabel("concordance (mice up / mice scored)")
        ax.set_title("scRNA MAX EFFECT Part1: Tacstd2-high → TJ_TISMO")
        ax.set_ylim(-0.05, 1.05)
        ax.axhline(1.0, color="0.6", lw=0.7, ls=":")
        fig1.tight_layout()
        fig1.savefig(fig / "scrna_tj_concordance.png", dpi=150)
        fig1.savefig(fig / "scrna_tj_concordance.pdf")
        plt.close(fig1)

    if winners:
        # Prefer the mean-Δ maximizer among perfect concordance for the bar plot.
        best = next(
            (w for w in winners if w["selection"] == "max_mean_delta_among_perfect_maxn"),
            next((w for w in winners if w["selection"] == "max_frac_up"), winners[0]),
        )
        best_scores = scores[
            (scores.gate == best["gate"])
            & (scores.tac_rule == best["tac_rule"])
            & (scores.min_high == best["min_high"])
            & (scores.min_low == best["min_low"])
            & (scores.set == best["set"])
        ].copy()
        best_scores.to_csv(tab / "winner_mouse_deltas.tsv", sep="\t", index=False)

        fig2, ax = plt.subplots(figsize=(7.0, 3.8))
        order = best_scores.sort_values("delta")
        colors = ["#54A24B" if d > 0 else "#E45756" for d in order["delta"]]
        ax.barh(order["mouse"], order["delta"], color=colors)
        ax.axvline(0, color="0.4", lw=0.8)
        ax.set_xlabel("TJ module Δ (Tacstd2-high − low)")
        ax.set_title(
            f"Winner: {best['set']} {best['gate']}/{best['tac_rule']} "
            f"floor {best['min_high']}/{best['min_low']} — {best['n_up']}/{best['n_mice']} up"
        )
        fig2.tight_layout()
        fig2.savefig(fig / "scrna_winner_mouse_deltas.png", dpi=150)
        fig2.savefig(fig / "scrna_winner_mouse_deltas.pdf")
        plt.close(fig2)

    summary = {
        "assay": "scRNA",
        "endpoint": "Tacstd2-high → TJ module concordance",
        "private_8kl": 0,
        "merged_with_private_8kl": False,
        "n_grid_specs": int(len(cons_df)),
        "n_eligible_ge4": int(len(eligible)),
        "locked_tj_tismo_floor20": locked.iloc[0].to_dict() if len(locked) else {},
        "locked_tj_tismo_floor30": locked_30.iloc[0].to_dict() if len(locked_30) else {},
        "winners": winners,
        "python": platform.python_version(),
        "scipy": scipy_version,
        "rule": (
            "Maximize mice-up / mice-scored on a fixed grid of gate × Tacstd2 rule × "
            "cell floor × frozen TJ module. Eligible specs need ≥4 mice. "
            "Tacstd2 is not an epithelial caller."
        ),
    }
    (tab / "scrna_summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")

    show = eligible.sort_values(["frac_up", "n_mice", "mean_delta"], ascending=False).head(15)
    print(show.to_string(index=False))
    print("--- winners ---")
    if not win_df.empty:
        cols = [
            "selection", "set", "gate", "tac_rule", "min_high", "min_low",
            "n_mice", "n_up", "frac_up", "n_up_sig", "mean_delta",
        ]
        print(win_df[cols].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
