#!/usr/bin/env python3
"""Ligand–receptor communication on GSE207422 (public processed UMI).

User A3 is taken as given: Hu lineage markers; malignant-like = epithelial
AND zero UMI of SFTPA2/AGER/SCGB1A1/SCGB3A1/TPPP3; T/NK = T or NK.

Primary executable method
    Documented CellPhoneDB-style score (Efremova 2020 / Garcia-Alonso 2022):
    partner expression = min(subunit means) on log1p(CP10k);
    score = mean of the two partner means.
    Patient-level paired Wilcoxon is the inferential unit.

Secondary method (if import succeeds)
    LIANA `cellphonedb` on a downsampled object. Labeled as LIANA, not CellChat.

CellChat is not run: R is not available in this environment. No CellChat
tables are written.
"""

from __future__ import annotations

import argparse
import gzip
import json
import traceback
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from scipy import stats
from statsmodels.stats.multitest import multipletests

HERE = Path(__file__).resolve().parents[1]

PAPER_GROUP = {
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


def log(msg: str) -> None:
    print(msg, flush=True)


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


def stream_matrix(path: Path, keep: set[str]) -> tuple[np.ndarray, list[str], dict[str, np.ndarray], np.ndarray]:
    log(f"[stream] {path} keep={len(keep)} symbols")
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        cells = np.array(header[1:], dtype=object)
        n = len(cells)
        totals = np.zeros(n, dtype=np.float64)
        store: dict[str, np.ndarray] = {}
        n_genes = 0
        for line in fh:
            gene, rest = line.split("\t", 1)
            n_genes += 1
            vals = np.fromstring(rest, sep="\t", dtype=np.float32)
            if vals.size != n:
                raise ValueError(f"column mismatch for {gene}: {vals.size} != {n}")
            totals += vals
            if gene in keep:
                store[gene] = vals
    log(f"[stream] cells={n} genes_in_file={n_genes} genes_kept={len(store)}")
    return cells, list(store), store, totals


def log1p_cp10k(umi: np.ndarray, total: np.ndarray) -> np.ndarray:
    with np.errstate(divide="ignore", invalid="ignore"):
        cp = np.where(total > 0, umi / total * 1e4, 0.0)
    return np.log1p(cp).astype(np.float32)


def partner_stats(units: list[str], logx: dict[str, np.ndarray], mask: np.ndarray) -> tuple[float, float, int]:
    """CellPhoneDB complex rule: min of subunit means; min of subunit %pos."""
    means, fracs = [], []
    n = int(mask.sum())
    if n == 0:
        return np.nan, np.nan, 0
    for g in units:
        if g not in logx:
            return np.nan, np.nan, n
        v = logx[g][mask]
        means.append(float(np.mean(v)))
        fracs.append(float(np.mean(v > 0)))
    return float(np.min(means)), float(np.min(fracs)), n


def score_pairs(
    pairs: pd.DataFrame,
    logx: dict[str, np.ndarray],
    sender: np.ndarray,
    receiver: np.ndarray,
    expr_prop: float,
) -> pd.DataFrame:
    rows = []
    for rec in pairs.itertuples(index=False):
        lig_u = str(rec.ligand).split("+")
        rec_u = str(rec.receptor).split("+")
        l_mean, l_frac, n_s = partner_stats(lig_u, logx, sender)
        r_mean, r_frac, n_r = partner_stats(rec_u, logx, receiver)
        if not np.isfinite(l_mean) or not np.isfinite(r_mean):
            continue
        pass_prop = (l_frac >= expr_prop) and (r_frac >= expr_prop)
        rows.append(
            {
                "ligand": rec.ligand,
                "receptor": rec.receptor,
                "pathway": rec.pathway,
                "pair_origin": rec.pair_origin,
                "directionality": rec.directionality,
                "n_sender": n_s,
                "n_receiver": n_r,
                "ligand_mean": l_mean,
                "receptor_mean": r_mean,
                "ligand_frac": l_frac,
                "receptor_frac": r_frac,
                "cpdb_mean_score": 0.5 * (l_mean + r_mean),
                "product_score": l_mean * r_mean,
                "pass_expr_prop": pass_prop,
            }
        )
    return pd.DataFrame(rows)


def high_low_masks(values: np.ndarray, malig: np.ndarray, rule: str) -> tuple[np.ndarray, np.ndarray, float]:
    v = values[malig]
    if rule == "tertile":
        lo_t, hi_t = np.quantile(v, [1 / 3, 2 / 3])
        high = malig & (values >= hi_t)
        low = malig & (values <= lo_t)
        return high, low, float(hi_t)
    thr = float(np.median(v))
    high = malig & (values >= thr)
    low = malig & (values < thr)
    return high, low, thr


def paired_delta(high_df: pd.DataFrame, low_df: pd.DataFrame, score_col: str) -> pd.DataFrame:
    key = ["ligand", "receptor", "pathway"]
    h = high_df[key + [score_col, "pass_expr_prop", "ligand_frac", "receptor_frac", "n_sender", "n_receiver"]].copy()
    l = low_df[key + [score_col, "pass_expr_prop", "ligand_frac", "receptor_frac", "n_sender", "n_receiver"]].copy()
    m = h.merge(l, on=key, suffixes=("_high", "_low"))
    m["delta_high_minus_low"] = m[f"{score_col}_high"] - m[f"{score_col}_low"]
    m["pass_either"] = m["pass_expr_prop_high"] | m["pass_expr_prop_low"]
    m["pass_both"] = m["pass_expr_prop_high"] & m["pass_expr_prop_low"]
    return m


def wilcoxon_safe(a, b) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if len(a) < 3 or np.allclose(a, b):
        return np.nan
    try:
        return float(stats.wilcoxon(a, b, zero_method="wilcox", alternative="two-sided").pvalue)
    except ValueError:
        return np.nan


def run_liana(adata, groupby: str, out_csv: Path, n_perms: int) -> str:
    try:
        import liana as li
    except Exception as exc:  # pragma: no cover
        return f"LIANA_IMPORT_FAILED: {exc}"
    try:
        li.mt.cellphonedb(
            adata,
            groupby=groupby,
            resource_name="cellphonedb",
            expr_prop=0.10,
            n_perms=n_perms,
            use_raw=False,
            verbose=True,
            key_added="liana_res",
        )
        res = adata.uns["liana_res"].copy()
        res.to_csv(out_csv, index=False)
        return f"LIANA_OK n_edges={len(res)} file={out_csv.name}"
    except Exception as exc:
        return f"LIANA_RUN_FAILED: {exc}\n{traceback.format_exc()}"


def make_figures(outdir: Path, ntab: pd.DataFrame, ranks: dict[str, pd.DataFrame]) -> None:
    figdir = outdir / "figures"
    figdir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    samples = ntab["Sample"].tolist()
    x = np.arange(len(samples))
    ax.bar(x - 0.2, ntab["n_malig_like"], width=0.4, label="malignant-like", color="#4C72B0")
    ax.bar(x + 0.2, ntab["n_T_NK"], width=0.4, label="T/NK", color="#DD8452")
    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=75, ha="right", fontsize=7)
    ax.set_ylabel("cells")
    ax.legend(frameon=False)
    ax.set_title("GSE207422 cells used (A3 malignant-like definition)")
    fig.tight_layout()
    fig.savefig(figdir / "n_cells_by_sample.png", dpi=160)
    fig.savefig(figdir / "n_cells_by_sample.pdf")
    plt.close(fig)

    for key, df in ranks.items():
        focus = df[df["pathway"].isin(["T_recruit", "IFN", "MHC_I"])].copy()
        if focus.empty:
            continue
        focus = focus.sort_values("median_delta")
        colors = {"T_recruit": "#4C72B0", "IFN": "#55A868", "MHC_I": "#C44E52"}
        fig, ax = plt.subplots(figsize=(8.0, max(3.2, 0.28 * len(focus) + 1.2)))
        y = np.arange(len(focus))
        ax.barh(
            y,
            focus["median_delta"],
            color=[colors.get(p, "#999") for p in focus["pathway"]],
            edgecolor="none",
        )
        ax.axvline(0, color="k", lw=0.8)
        ax.set_yticks(y)
        labels = [f"{a}–{b}  ({p})" for a, b, p in zip(focus["ligand"], focus["receptor"], focus["pathway"])]
        ax.set_yticklabels(labels, fontsize=7)
        ax.set_xlabel("median patient Δ score (high − low)")
        ax.set_title(key.replace("_", " "))
        fig.tight_layout()
        fig.savefig(figdir / f"{key}.png", dpi=160)
        fig.savefig(figdir / f"{key}.pdf")
        plt.close(fig)


def write_results_md(
    outdir: Path,
    summary: dict,
    ntab: pd.DataFrame,
    ranks: dict[str, pd.DataFrame],
    liana_notes: list[str],
) -> None:
    lines = [
        "# GSE207422 malignant TACSTD2 / CLDN4 ligand–receptor communication",
        "",
        "Public processed UMI from Hu et al., *Genome Med* 2023 (GSE207422; PMID 36869384).",
        "User A3 lineage and malignant-like rules are used as given. This note reports",
        "outgoing signals from malignant-like cells to T/NK and incoming signals from T/NK",
        "to malignant-like cells, split on TACSTD2 and (separately) CLDN4.",
        "",
        "## What was run",
        "",
        f"- **Documented CellPhoneDB-style score** (not CellChat): {summary['n_pairs_scored']} pairs from CellPhoneDB v5 + curated overlay.",
        f"- **LIANA:** {summary['liana_status']}",
        f"- **CellChat:** {summary['cellchat_status']}",
        "",
        "## n cells / patients",
        "",
        f"- Matrix: **{summary['n_cells']:,}** cells, **{summary['n_patients']}** patients (BD Rhapsody WTA).",
        f"- Lineage (A3): epithelial={summary['n_epithelial']:,}; malignant-like={summary['n_malig']:,}; T={summary['n_T']:,}; NK={summary['n_NK']:,}; T/NK={summary['n_tnk']:,}.",
        f"- TACSTD2 split (global median log1p CP10k = {summary['tacstd2_threshold']:.3f}): high={summary['n_tacstd2_high']:,} cells, low={summary['n_tacstd2_low']:,} cells.",
        f"- CLDN4 split (global median log1p CP10k = {summary['cldn4_threshold']:.3f}): high={summary['n_cldn4_high']:,} cells, low={summary['n_cldn4_low']:,} cells.",
        f"- Patients entering the paired TACSTD2 test (≥{summary['min_malig_per_state']} malignant cells per bin and ≥{summary['min_tnk']} T/NK): **{summary['n_patients_tacstd2_paired']}**.",
        f"- Patients entering the paired CLDN4 test: **{summary['n_patients_cldn4_paired']}**.",
        "",
        "Per-sample counts are in `n_cells_patients.tsv`.",
        "",
        "## Score",
        "",
        "On log1p(CP10k), each partner’s expression is the **minimum subunit mean**",
        "(CellPhoneDB complex rule). The pair score is the **mean of the two partner means**",
        "(Efremova et al. 2020 *Nat Protoc*; Garcia-Alonso et al. 2022 *Nat Protoc*).",
        "A pair is flagged `pass_expr_prop` when both partners are detected in ≥10% of",
        "cells in their group. Patient-level tests use that patient’s own T/NK and that",
        "patient’s high vs low malignant-like cells. Cells are not treated as replicates.",
        "",
        "## TACSTD2-high vs TACSTD2-low",
        "",
    ]

    def _block(title: str, key: str) -> list[str]:
        df = ranks[key]
        focus = df[df["pathway"].isin(["T_recruit", "IFN", "MHC_I"])].copy()
        if focus.empty:
            return [f"### {title}", "", "No T-recruit / IFN / MHC-I pairs with ≥3 paired patients.", ""]
        focus = focus.sort_values(["pathway", "median_delta"])
        out = [
            f"### {title}",
            "",
            "Median patient Δ = high − low. Negative = weaker from/to the high state.",
            "",
            "| Pathway | Pair | n patients | median Δ | Wilcoxon p | FDR |",
            "|---|---|---:|---:|---:|---:|",
        ]
        for r in focus.itertuples(index=False):
            p = f"{r.pval:.3g}" if np.isfinite(r.pval) else "NA"
            q = f"{r.padj:.3g}" if np.isfinite(r.padj) else "NA"
            out.append(
                f"| {r.pathway} | {r.ligand}–{r.receptor} | {int(r.n_patients)} | {r.median_delta:+.3f} | {p} | {q} |"
            )
        n_neg = int((focus["median_delta"] < 0).sum())
        n_sig = int(((focus["padj"] < 0.05) & (focus["median_delta"] < 0)).sum()) if focus["padj"].notna().any() else 0
        n_sig_any = int((focus["padj"] < 0.05).sum()) if focus["padj"].notna().any() else 0
        out += [
            "",
            f"{n_neg}/{len(focus)} focus pairs have median Δ < 0 (weaker in the high state). "
            f"{n_sig} of those are FDR < 0.05. {n_sig_any}/{len(focus)} focus pairs are FDR < 0.05 in either direction.",
            "",
        ]
        return out

    lines += _block("Outgoing malignant → T/NK", "tacstd2_outgoing")
    lines += _block("Incoming T/NK → malignant", "tacstd2_incoming")
    lines += ["## CLDN4-high vs CLDN4-low", ""]
    lines += _block("Outgoing malignant → T/NK", "cldn4_outgoing")
    lines += _block("Incoming T/NK → malignant", "cldn4_incoming")

    lines += [
        "## LIANA (secondary)",
        "",
    ]
    for note in liana_notes:
        lines.append(f"- {note}")
    lines += [
        "",
        "LIANA permutation p-values test cluster-label specificity inside one pooled object,",
        "not reproducibility across patients. Patient-level ranks above are the ones to cite.",
        "",
        "## CellChat",
        "",
        "R / CellChat were not available. No CellChat probability tables were written.",
        "",
        "## Readout vs a TROP2-high reduced T-recruit / IFN / MHC-I pattern",
        "",
        summary["trend_sentence"],
        "",
        "## Limits",
        "",
        "- BD Rhapsody, not 10x. GEO deposited no per-cell labels; A3 marker gates are used as given.",
        "- Malignant-like is a normal-lung-marker exclusion, not public CopyKAT calls.",
        "- Ambient RNA cannot be re-estimated from the processed matrix.",
        "- Chemokine dropout is high; `pass_expr_prop` should be read with the ranks.",
        f"- Paired n = {summary['n_patients_tacstd2_paired']} (TACSTD2) and {summary['n_patients_cldn4_paired']} (CLDN4).",
        "",
    ]
    (outdir / "RESULTS.md").write_text("\n".join(lines))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--datadir", type=Path, default=HERE / "data")
    ap.add_argument("--outdir", type=Path, default=HERE / "results")
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    cfg = yaml.safe_load((HERE / "config" / "gene_sets.yaml").read_text())
    P = cfg["params"]
    pairs = pd.read_csv(HERE / "resources" / "cellphonedb_v5_lr_pairs.tsv", sep="\t")

    keep = set()
    for block in cfg["lineage_markers"].values():
        keep.update(block)
    keep.update(cfg["normal_lung"])
    keep.update(cfg["state_genes"])
    keep.update(Path(HERE / "resources" / "lr_genes.txt").read_text().split())
    try:
        import liana as li

        res = li.rs.select_resource("cellphonedb")
        keep.update(res["ligand"].astype(str))
        keep.update(res["receptor"].astype(str))
        keep.update(g for s in res["ligand"].astype(str) for g in s.split("_"))
        keep.update(g for s in res["receptor"].astype(str) for g in s.split("_"))
        liana_import = "ok"
    except Exception as exc:
        liana_import = f"failed: {exc}"

    mtx = args.datadir / "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"
    cells, _genes, expr, total = stream_matrix(mtx, keep)
    n = len(cells)
    sample = np.array([str(c).rsplit("_", 1)[0] for c in cells])
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
    tac = logx["TACSTD2"] if "TACSTD2" in logx else np.zeros(n, dtype=np.float32)
    cld = logx["CLDN4"] if "CLDN4" in logx else np.zeros(n, dtype=np.float32)
    tac_hi, tac_lo, tac_thr = high_low_masks(tac, is_malig, P["high_rule"])
    cld_hi, cld_lo, cld_thr = high_low_masks(cld, is_malig, P["high_rule"])

    meta = pd.read_excel(args.datadir / "GSE207422_NSCLC_scRNAseq_metadata.xlsx")
    meta = meta.dropna(subset=["Sample"]).copy()
    meta = meta[meta["Sample"].astype(str).str.startswith("BD_immune")]
    meta["paper_group"] = meta["Sample"].map(PAPER_GROUP)

    ntab_rows = []
    for s in sorted(set(sample)):
        m = sample == s
        ntab_rows.append(
            {
                "Sample": s,
                "patient": meta.set_index("Sample")["Patient"].get(s, s),
                "paper_group": PAPER_GROUP.get(s, ""),
                "n_cells": int(m.sum()),
                "n_epithelial": int((m & is_epi).sum()),
                "n_malig_like": int((m & is_malig).sum()),
                "n_T": int((m & is_t).sum()),
                "n_NK": int((m & is_nk).sum()),
                "n_T_NK": int((m & is_tnk).sum()),
                "n_tacstd2_high": int((m & tac_hi).sum()),
                "n_tacstd2_low": int((m & tac_lo).sum()),
                "n_cldn4_high": int((m & cld_hi).sum()),
                "n_cldn4_low": int((m & cld_lo).sum()),
            }
        )
    ntab = pd.DataFrame(ntab_rows)
    ntab.to_csv(args.outdir / "n_cells_patients.tsv", sep="\t", index=False)

    # Pooled scores (descriptive).
    pooled = {}
    for name, hi, lo in (
        ("tacstd2", tac_hi, tac_lo),
        ("cldn4", cld_hi, cld_lo),
    ):
        out_hi = score_pairs(pairs, logx, hi, is_tnk, P["expr_prop"])
        out_lo = score_pairs(pairs, logx, lo, is_tnk, P["expr_prop"])
        in_hi = score_pairs(pairs, logx, is_tnk, hi, P["expr_prop"])
        in_lo = score_pairs(pairs, logx, is_tnk, lo, P["expr_prop"])
        out_d = paired_delta(out_hi, out_lo, "cpdb_mean_score")
        in_d = paired_delta(in_hi, in_lo, "cpdb_mean_score")
        out_d.to_csv(args.outdir / f"pooled_outgoing_{name}.tsv", sep="\t", index=False)
        in_d.to_csv(args.outdir / f"pooled_incoming_{name}.tsv", sep="\t", index=False)
        pooled[name] = (out_d, in_d)

    # Patient-level.
    def patient_table(hi: np.ndarray, lo: np.ndarray, direction: str) -> pd.DataFrame:
        chunks = []
        for s in sorted(set(sample)):
            m = sample == s
            h, l, t = m & hi, m & lo, m & is_tnk
            if h.sum() < P["min_malig_per_state"] or l.sum() < P["min_malig_per_state"] or t.sum() < P["min_tnk"]:
                continue
            if direction == "outgoing":
                hdf = score_pairs(pairs, logx, h, t, P["expr_prop"])
                ldf = score_pairs(pairs, logx, l, t, P["expr_prop"])
            else:
                hdf = score_pairs(pairs, logx, t, h, P["expr_prop"])
                ldf = score_pairs(pairs, logx, t, l, P["expr_prop"])
            d = paired_delta(hdf, ldf, "cpdb_mean_score")
            d["Sample"] = s
            d["n_high"] = int(h.sum())
            d["n_low"] = int(l.sum())
            d["n_tnk"] = int(t.sum())
            chunks.append(d)
        return pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame()

    ranks: dict[str, pd.DataFrame] = {}
    n_paired = {}
    for state, hi, lo in (("tacstd2", tac_hi, tac_lo), ("cldn4", cld_hi, cld_lo)):
        for direction in ("outgoing", "incoming"):
            raw = patient_table(hi, lo, direction)
            raw.to_csv(args.outdir / f"patient_{state}_{direction}.tsv", sep="\t", index=False)
            if raw.empty:
                ranks[f"{state}_{direction}"] = pd.DataFrame(
                    columns=["ligand", "receptor", "pathway", "n_patients", "median_delta", "pval", "padj"]
                )
                n_paired[f"{state}_{direction}"] = 0
                continue
            n_paired[f"{state}_{direction}"] = int(raw["Sample"].nunique())
            rows = []
            for (lig, recp, path), sub in raw.groupby(["ligand", "receptor", "pathway"], observed=True):
                # Prefer pairs that pass expr_prop in at least one state in ≥ half of patients.
                keep_sub = sub[sub["pass_either"]]
                if len(keep_sub) < 3:
                    continue
                # Focus pathways always kept if present; others need pass_both in ≥3 patients.
                if path not in {"T_recruit", "IFN", "MHC_I"} and keep_sub["pass_both"].sum() < 3:
                    continue
                delta = keep_sub["delta_high_minus_low"].to_numpy()
                p = wilcoxon_safe(
                    keep_sub["cpdb_mean_score_high"],
                    keep_sub["cpdb_mean_score_low"],
                )
                rows.append(
                    {
                        "ligand": lig,
                        "receptor": recp,
                        "pathway": path,
                        "n_patients": int(keep_sub["Sample"].nunique()),
                        "median_delta": float(np.median(delta)),
                        "mean_delta": float(np.mean(delta)),
                        "mean_score_high": float(keep_sub["cpdb_mean_score_high"].mean()),
                        "mean_score_low": float(keep_sub["cpdb_mean_score_low"].mean()),
                        "frac_pass_high": float(keep_sub["pass_expr_prop_high"].mean()),
                        "frac_pass_low": float(keep_sub["pass_expr_prop_low"].mean()),
                        "pval": p,
                    }
                )
            tab = pd.DataFrame(rows)
            if len(tab) and tab["pval"].notna().any():
                mask = tab["pval"].notna()
                tab.loc[mask, "padj"] = multipletests(tab.loc[mask, "pval"], method="fdr_bh")[1]
            else:
                tab["padj"] = np.nan
            tab = tab.sort_values(["pathway", "median_delta"])
            tab.to_csv(args.outdir / f"ranks_{state}_{direction}.tsv", sep="\t", index=False)
            ranks[f"{state}_{direction}"] = tab

    # Pathway-level summary (mean of pair median-deltas).
    path_rows = []
    for key, tab in ranks.items():
        for path, sub in tab.groupby("pathway"):
            if path not in {"T_recruit", "IFN", "MHC_I"}:
                continue
            path_rows.append(
                {
                    "contrast": key,
                    "pathway": path,
                    "n_pairs": int(len(sub)),
                    "n_pairs_delta_neg": int((sub["median_delta"] < 0).sum()),
                    "median_of_pair_deltas": float(sub["median_delta"].median()),
                    "n_pairs_fdr05_neg": int(((sub["padj"] < 0.05) & (sub["median_delta"] < 0)).sum())
                    if sub["padj"].notna().any()
                    else 0,
                }
            )
    path_df = pd.DataFrame(path_rows)
    path_df.to_csv(args.outdir / "pathway_summary.tsv", sep="\t", index=False)

    # LIANA secondary (downsampled).
    liana_notes = [f"import: {liana_import}"]
    liana_status = "not_run"
    try:
        import anndata as ad
        import scanpy as sc

        def build_adata(hi, lo, label):
            mask = hi | lo | is_t | is_nk
            genes = sorted(g for g in logx if g in expr)
            X = np.vstack([logx[g][mask] for g in genes]).T
            obs = pd.DataFrame(
                {
                    "sample": sample[mask],
                    "lineage": lineage[mask],
                },
                index=pd.Index(cells[mask], name="cell"),
            )
            grp = np.array(["other"] * int(mask.sum()), dtype=object)
            grp[hi[mask]] = f"Malig_{label}high"
            grp[lo[mask]] = f"Malig_{label}low"
            grp[is_t[mask]] = "T"
            grp[is_nk[mask]] = "NK"
            obs["cc_group"] = grp
            # downsample
            rng = np.random.default_rng(P["random_seed"])
            keep_idx = []
            for g, idx in obs.groupby("cc_group", observed=True).indices.items():
                if g == "other":
                    continue
                if len(idx) > P["liana_max_cells_per_group"]:
                    idx = rng.choice(idx, size=P["liana_max_cells_per_group"], replace=False)
                keep_idx.append(np.asarray(idx))
            keep_idx = np.sort(np.concatenate(keep_idx))
            a = ad.AnnData(X=X[keep_idx], obs=obs.iloc[keep_idx].copy(), var=pd.DataFrame(index=genes))
            a.obs["cc_group"] = pd.Categorical(a.obs["cc_group"])
            a.uns["log1p"] = {"base": None}
            return a

        for label, hi, lo in (("TACSTD2", tac_hi, tac_lo), ("CLDN4", cld_hi, cld_lo)):
            adata = build_adata(hi, lo, label)
            counts = adata.obs["cc_group"].value_counts().to_dict()
            log(f"[liana] {label} groups {counts}")
            note = run_liana(
                adata,
                "cc_group",
                args.outdir / f"liana_cellphonedb_{label.lower()}.csv",
                P["liana_n_perms"],
            )
            liana_notes.append(f"{label}: {note}; downsampled groups={counts}")
            if note.startswith("LIANA_OK"):
                liana_status = "ran_cellphonedb_method"
                res = pd.read_csv(args.outdir / f"liana_cellphonedb_{label.lower()}.csv")
                senders = [f"Malig_{label}high", f"Malig_{label}low"]
                focus = res[res["source"].isin(senders) & res["target"].isin(["T", "NK"])].copy()
                incoming = res[res["source"].isin(["T", "NK"]) & res["target"].isin(senders)].copy()
                focus.to_csv(args.outdir / f"liana_{label.lower()}_outgoing_tnk.csv", index=False)
                incoming.to_csv(args.outdir / f"liana_{label.lower()}_incoming_tnk.csv", index=False)
    except Exception as exc:
        liana_notes.append(f"LIANA wrapper failed: {exc}")
        liana_status = f"failed: {exc}"

    # Trend sentence from TACSTD2 outgoing focus pairs.
    foc = ranks.get("tacstd2_outgoing", pd.DataFrame())
    foc = foc[foc["pathway"].isin(["T_recruit", "IFN", "MHC_I"])] if len(foc) else foc
    if len(foc):
        n_neg = int((foc["median_delta"] < 0).sum())
        n_sig = int(((foc["padj"] < 0.05) & (foc["median_delta"] < 0)).sum()) if foc["padj"].notna().any() else 0
        by = (
            foc.groupby("pathway")
            .apply(lambda s: f"{int((s.median_delta < 0).sum())}/{len(s)} Δ<0, median Δ={s.median_delta.median():+.3f}", include_groups=False)
            .to_dict()
        )
        trend = (
            f"On the patient-level CellPhoneDB-style score, TACSTD2-high vs TACSTD2-low outgoing "
            f"T-recruit/IFN/MHC-I pairs: {n_neg}/{len(foc)} have median Δ < 0; {n_sig} reach FDR < 0.05. "
            f"By axis: {by}. "
            "This is the observed rank in this 15-patient public object, not a general rule."
        )
    else:
        trend = "Too few patient-level T-recruit/IFN/MHC-I pairs passed filters to rank a TROP2-high reduction."

    summary = {
        "dataset": "GSE207422",
        "n_cells": int(n),
        "n_patients": 15,
        "n_epithelial": int(is_epi.sum()),
        "n_malig": int(is_malig.sum()),
        "n_T": int(is_t.sum()),
        "n_NK": int(is_nk.sum()),
        "n_tnk": int(is_tnk.sum()),
        "n_tacstd2_high": int(tac_hi.sum()),
        "n_tacstd2_low": int(tac_lo.sum()),
        "n_cldn4_high": int(cld_hi.sum()),
        "n_cldn4_low": int(cld_lo.sum()),
        "tacstd2_threshold": tac_thr,
        "cldn4_threshold": cld_thr,
        "min_malig_per_state": P["min_malig_per_state"],
        "min_tnk": P["min_tnk"],
        "n_patients_tacstd2_paired": n_paired.get("tacstd2_outgoing", 0),
        "n_patients_cldn4_paired": n_paired.get("cldn4_outgoing", 0),
        "n_pairs_scored": int(len(pairs)),
        "liana_status": liana_status,
        "cellchat_status": "not_run_R_unavailable",
        "liana_notes": liana_notes,
        "trend_sentence": trend,
        "a3_taken_as_given": True,
        "method": "cellphonedb_mean_of_means_on_log1p_cp10k",
    }
    (args.outdir / "summary.json").write_text(json.dumps(summary, indent=2))
    make_figures(args.outdir, ntab, ranks)
    write_results_md(outdir=args.outdir, summary=summary, ntab=ntab, ranks=ranks, liana_notes=liana_notes)
    log(json.dumps({k: summary[k] for k in ("n_malig", "n_tnk", "n_patients_tacstd2_paired", "liana_status", "trend_sentence")}, indent=2))
    log(f"[done] {args.outdir}")


if __name__ == "__main__":
    main()
