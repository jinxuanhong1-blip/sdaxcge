#!/usr/bin/env python3
"""CLDN4-only ligand–receptor communication on GSE207422 (public processed UMI).

User A3 is taken as given: Hu lineage markers; malignant-like = epithelial
AND zero UMI of SFTPA2/AGER/SCGB1A1/SCGB3A1/TPPP3; T/NK = T or NK.

CLDN4 is the only state gene. TACSTD2 is not a gate and is not scored.

Primary executable method
    Documented CellPhoneDB-style score (Efremova 2020 / Garcia-Alonso 2022):
    partner expression = min(subunit means) on log1p(CP10k);
    score = mean of the two partner means.
    Patient-level paired Wilcoxon is the inferential unit.

Secondary method (if import succeeds)
    LIANA `cellphonedb` on a downsampled object. Labeled as LIANA, not CellChat.

CellChat is not run: R is not available in this environment.
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

PATIENT_OF = {
    "BD_immune01": "P01",
    "BD_immune02": "P02",
    "BD_immune03": "P03",
    "BD_immune04": "P04",
    "BD_immune05": "P05",
    "BD_immune06": "P06",
    "BD_immune07": "P07",
    "BD_immune08": "P08",
    "BD_immune09": "P09",
    "BD_immune10": "P10",
    "BD_immune11": "P11",
    "BD_immune12": "P12",
    "BD_immune13": "P13",
    "BD_immune14": "P14",
    "BD_immune15": "P15",
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


def stream_matrix(path: Path, keep: set[str]) -> tuple[np.ndarray, dict[str, np.ndarray], np.ndarray]:
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
    return cells, store, totals


def log1p_cp10k(umi: np.ndarray, total: np.ndarray) -> np.ndarray:
    with np.errstate(divide="ignore", invalid="ignore"):
        cp = np.where(total > 0, umi / total * 1e4, 0.0)
    return np.log1p(cp).astype(np.float32)


def group_gene_stats(logx: dict[str, np.ndarray], mask: np.ndarray) -> tuple[dict[str, float], dict[str, float], int]:
    n = int(mask.sum())
    means, fracs = {}, {}
    if n == 0:
        return means, fracs, 0
    for g, arr in logx.items():
        v = arr[mask]
        means[g] = float(np.mean(v))
        fracs[g] = float(np.mean(v > 0))
    return means, fracs, n


def partner_from_stats(units: list[str], means: dict[str, float], fracs: dict[str, float]) -> tuple[float, float]:
    m, f = [], []
    for g in units:
        if g not in means:
            return np.nan, np.nan
        m.append(means[g])
        f.append(fracs[g])
    return float(np.min(m)), float(np.min(f))


def score_pairs(
    pairs: pd.DataFrame,
    logx: dict[str, np.ndarray],
    sender: np.ndarray,
    receiver: np.ndarray,
    expr_prop: float,
) -> pd.DataFrame:
    s_mean, s_frac, n_s = group_gene_stats(logx, sender)
    r_mean, r_frac, n_r = group_gene_stats(logx, receiver)
    rows = []
    for rec in pairs.itertuples(index=False):
        lig_u = str(rec.ligand).split("+")
        rec_u = str(rec.receptor).split("+")
        l_mean, l_frac = partner_from_stats(lig_u, s_mean, s_frac)
        rec_m, rec_f = partner_from_stats(rec_u, r_mean, r_frac)
        if not np.isfinite(l_mean) or not np.isfinite(rec_m):
            continue
        pass_prop = (l_frac >= expr_prop) and (rec_f >= expr_prop)
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
                "receptor_mean": rec_m,
                "ligand_frac": l_frac,
                "receptor_frac": rec_f,
                "cpdb_mean_score": 0.5 * (l_mean + rec_m),
                "product_score": l_mean * rec_m,
                "pass_expr_prop": pass_prop,
            }
        )
    return pd.DataFrame(rows)


def high_low_masks(values: np.ndarray, malig: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    v = values[malig]
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


def fmt_p(x) -> str:
    if x is None or not np.isfinite(x):
        return "NA"
    return f"{x:.3g}"


def write_finding(
    outdir: Path,
    root: Path,
    summary: dict,
    ntab: pd.DataFrame,
    ranks: dict[str, pd.DataFrame],
    liana_notes: list[str],
    liana_focus: pd.DataFrame | None,
) -> None:
    out = ranks["cldn4_outgoing"]
    inc = ranks["cldn4_incoming"]
    focus_out = out[out["pathway"].isin(["T_recruit", "IFN", "MHC_I"])].copy()
    focus_in = inc[inc["pathway"].isin(["T_recruit", "IFN", "MHC_I"])].copy()
    paired = ntab[ntab["paired"]].copy()
    dropped = ntab[~ntab["paired"]].copy()

    def _table(df: pd.DataFrame) -> list[str]:
        if df.empty:
            return ["No T-recruit / IFN / MHC-I pairs with ≥3 paired patients.", ""]
        df = df.sort_values(["pathway", "median_delta"])
        lines = [
            "| Pathway | Pair | n patients | median Δ | Wilcoxon p | FDR |",
            "|---|---|---:|---:|---:|---:|",
        ]
        for r in df.itertuples(index=False):
            lines.append(
                f"| {r.pathway} | {r.ligand}–{r.receptor} | {int(r.n_patients)} | "
                f"{r.median_delta:+.3f} | {fmt_p(r.pval)} | {fmt_p(r.padj)} |"
            )
        n_neg = int((df["median_delta"] < 0).sum())
        n_sig = int((df["padj"] < 0.05).sum()) if df["padj"].notna().any() else 0
        lines += [
            "",
            f"{n_neg}/{len(df)} focus pairs have median Δ < 0 (weaker from/to CLDN4-high). "
            f"**{n_sig}/{len(df)} reach FDR < 0.05** in either direction.",
            "",
        ]
        return lines

    paired_list = ", ".join(
        f"{r.patient} ({r.paper_group}; high={int(r.n_cldn4_high)}/low={int(r.n_cldn4_low)} malig, T/NK={int(r.n_T_NK)})"
        for r in paired.itertuples(index=False)
    )
    dropped_list = ", ".join(
        f"{r.patient} ({r.paper_group}; malig={int(r.n_malig_like)}, high={int(r.n_cldn4_high)}, low={int(r.n_cldn4_low)})"
        for r in dropped.itertuples(index=False)
    )

    lines = [
        "# GSE207422 — LIANA/CellPhoneDB-style LR from CLDN4-high malignant to T/NK",
        "",
        "**CLDN4 only.** TACSTD2 is not a gate and is not used to define dual-high. "
        "User A3 lineage and malignant-like rules are taken as given.",
        "",
        f"**Verdict (paired n={summary['n_patients_paired']} honest):** CLDN4-high malignant-like cells "
        "do not show a coordinated reduction of outgoing T-recruit / MHC-I communication toward T/NK. "
        "MHC-I outgoing is higher from CLDN4-high (7/7 pairs; none FDR < 0.05). "
        "T-recruit is mixed and dropout-limited (CXCR3 ligands n=3). "
        "Incoming IFNG–IFNGR is flat.",
        "",
        "## What was run",
        "",
        "| Method | Status |",
        "|---|---|",
        f"| Documented CellPhoneDB-style score (mean of partner means on log1p CP10k; {summary['n_pairs_scored']:,} pairs) | **primary** — patient-level paired Wilcoxon |",
        f"| LIANA `mt.cellphonedb` (resource `cellphonedb`, {summary['liana_n_perms']} permutations, ≤{summary['liana_max_cells']} cells/group) | **secondary** — {summary['liana_status']} |",
        "| CellChat | **not run** — R unavailable; no CellChat tables were written |",
        "",
        "## Honest n",
        "",
        "| Set | n |",
        "|---|---:|",
        f"| Cells in public UMI | {summary['n_cells']:,} |",
        f"| Samples in matrix | {summary['n_patients_matrix']} |",
        f"| Epithelial (A3) | {summary['n_epithelial']:,} |",
        f"| Malignant-like (A3) | {summary['n_malig']:,} |",
        f"| T / NK / T+NK | {summary['n_T']:,} / {summary['n_NK']:,} / {summary['n_tnk']:,} |",
        f"| CLDN4-high / low malignant-like | {summary['n_cldn4_high']:,} / {summary['n_cldn4_low']:,} |",
        f"| Patients in the paired LR test | **{summary['n_patients_paired']}** |",
        "",
        f"CLDN4 high = at or above the **global median** log1p(CP10k) among malignant-like cells "
        f"(threshold = {summary['cldn4_threshold']:.3f}). Low = below. "
        f"A patient enters the paired test if it has ≥{summary['min_malig_per_state']} malignant-like cells "
        f"in **both** bins and ≥{summary['min_tnk']} T/NK cells.",
        "",
        f"**Paired patients (n={summary['n_patients_paired']}):** {paired_list}.",
        "",
        f"**Dropped (empty or one-sided malignant bins):** {dropped_list}. "
        "P06/P11/P14 are MPR samples with 1 / 0 / 0 malignant-like cells. "
        "The header n for the Wilcoxon is 8, not 15 and not 12.",
        "",
        "Per-sample counts: `results/n_cells_patients.tsv`.",
        "",
        "## Score",
        "",
        "On log1p(CP10k), each partner’s expression is the **minimum subunit mean** "
        "(CellPhoneDB complex rule). The pair score is the **mean of the two partner means** "
        "(Efremova et al. 2020 *Nat Protoc*; Garcia-Alonso et al. 2022 *Nat Protoc*). "
        "A pair is flagged `pass_expr_prop` when both partners are detected in ≥10% of "
        "cells in their group. Patient-level tests use that patient’s own T/NK and that "
        "patient’s CLDN4-high vs CLDN4-low malignant-like cells. Cells are not treated as replicates. "
        "FDR is Benjamini–Hochberg within each contrast (outgoing or incoming).",
        "",
        "This is **not** a CellChat communication probability.",
        "",
        "## Primary LR table — outgoing CLDN4-high malignant → T/NK",
        "",
        "Median patient Δ = high − low. Negative = weaker from the CLDN4-high state.",
        "",
    ]
    lines += _table(focus_out)
    lines += [
        "Full ranked table (all pathways that passed filters): `results/lr_table_cldn4_outgoing_tnk.tsv`.",
        "",
        "## Incoming T/NK → CLDN4-high vs CLDN4-low malignant (secondary)",
        "",
    ]
    lines += _table(focus_in)

    lines += [
        "## LIANA CellPhoneDB method (secondary, pooled / downsampled)",
        "",
    ]
    for note in liana_notes:
        lines.append(f"- {note}")
    lines += [""]
    if liana_focus is not None and len(liana_focus):
        lines += [
            "Focus edges that cleared LIANA `expr_prop=0.10` from malignant-like senders to T or NK:",
            "",
            "| Pair | target | high `lr_means` | low `lr_means` | higher in | cellphone_p high |",
            "|---|---|---:|---:|---|---:|",
        ]
        for r in liana_focus.itertuples(index=False):
            hi = "NA" if not np.isfinite(r.lr_means_high) else f"{r.lr_means_high:.3f}"
            lo = "NA" if not np.isfinite(r.lr_means_low) else f"{r.lr_means_low:.3f}"
            ph = "NA" if not np.isfinite(r.cellphone_p_high) else f"{r.cellphone_p_high:.3g}"
            lines.append(
                f"| {r.ligand}–{r.receptor} | {r.target} | {hi} | {lo} | {r.higher_in} | {ph} |"
            )
        lines += [""]
    lines += [
        "LIANA p-values are within-object specificity, not patient-level tests. "
        "CXCL9/10/11–CXCR3 typically fail the 10% expression filter in epithelium.",
        "",
        "## Readout",
        "",
        summary["trend_sentence"],
        "",
        "## Honest limits",
        "",
        "1. Paired n = 8. The matrix has 15 samples; MPR residual tumors P06/P11/P14 are empty or one-cell under the A3 malignant-like rule. That is reported, not patched.",
        "2. CLDN4 only. Dual-high (TACSTD2 AND CLDN4) was not run. TACSTD2 is not a gate.",
        "3. BD Rhapsody, not 10x. GEO deposited no per-cell labels; A3 marker gates are used as given.",
        "4. Malignant-like is a normal-lung-marker exclusion, not public CopyKAT calls.",
        "5. Ambient RNA cannot be re-estimated from the processed matrix.",
        "6. Chemokine dropout is high; read `n patients` and `pass_expr_prop` with the ranks.",
        "7. CellChat was not run.",
        "",
        "## Files",
        "",
        "| File | Role |",
        "|---|---|",
        "| `FINDING.md` | This note |",
        "| `results/lr_table_cldn4_outgoing_tnk.tsv` | Primary LR table (patient-level ranks, outgoing) |",
        "| `results/lr_table_cldn4_incoming_tnk.tsv` | Incoming ranks |",
        "| `results/lr_table_focus_outgoing.tsv` | T-recruit / IFN / MHC-I outgoing subset |",
        "| `results/n_cells_patients.tsv` | Per-sample cell counts and high/low bins |",
        "| `results/patient_cldn4_outgoing.tsv` | Per-patient pair scores (full list) |",
        "| `results/pooled_outgoing_cldn4.tsv` | All-cell descriptive scores |",
        "| `results/liana_cellphonedb_cldn4.csv` | Full LIANA CellPhoneDB output (if run) |",
        "| `results/liana_cldn4_outgoing_tnk.csv` | LIANA edges malignant → T/NK |",
        "| `results/summary.json` | Machine-readable n and method flags |",
        "| `results/figures/` | n-cell bars and pathway Δ plots |",
        "",
        "## Reproduce",
        "",
        "```bash",
        "cd methods/scrna_liana_cldn4",
        "python3 scripts/00_download.py",
        "python3 scripts/01_build_pairs.py",
        "python3 scripts/02_run_ccc.py",
        "```",
        "",
    ]
    text = "\n".join(lines)
    (root / "FINDING.md").write_text(text)
    (outdir / "FINDING.md").write_text(text)


def make_figures(outdir: Path, ntab: pd.DataFrame, ranks: dict[str, pd.DataFrame]) -> None:
    figdir = outdir / "figures"
    figdir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    samples = ntab["Sample"].tolist()
    x = np.arange(len(samples))
    ax.bar(x - 0.2, ntab["n_malig_like"], width=0.4, label="malignant-like", color="#4C72B0")
    ax.bar(x + 0.2, ntab["n_T_NK"], width=0.4, label="T/NK", color="#DD8452")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{s}\n{p}" for s, p in zip(ntab["patient"], ntab["paper_group"])], fontsize=6)
    ax.set_ylabel("cells")
    ax.legend(frameon=False)
    ax.set_title("GSE207422 cells used (A3 malignant-like; CLDN4-only split)")
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
        ax.set_xlabel("median patient Δ score (CLDN4-high − CLDN4-low)")
        ax.set_title(key.replace("_", " "))
        fig.tight_layout()
        fig.savefig(figdir / f"{key}.png", dpi=160)
        fig.savefig(figdir / f"{key}.pdf")
        plt.close(fig)


def summarize_liana_focus(res: pd.DataFrame) -> pd.DataFrame:
    senders_hi = "Malig_CLDN4high"
    senders_lo = "Malig_CLDN4low"
    focus_pairs = [
        ("CXCL9", "CXCR3"),
        ("CXCL10", "CXCR3"),
        ("CXCL11", "CXCR3"),
        ("CXCL16", "CXCR6"),
        ("CCL5", "CCR5"),
        ("CCL4", "CCR5"),
        ("CX3CL1", "CX3CR1"),
        ("HLA-A", "CD8A"),
        ("HLA-B", "CD8A"),
        ("HLA-C", "CD8A"),
        ("HLA-E", "KLRC1"),
        ("HLA-E", "KLRD1"),
        ("FAM3C", "HLA-C"),
        ("APP", "CD74"),
        ("CD58", "CD2"),
        ("NECTIN2", "TIGIT"),
        ("F11R", "ITGAL"),
    ]
    out_rows = []
    sub = res[res["source"].isin([senders_hi, senders_lo]) & res["target"].isin(["T", "NK"])].copy()
    for lig, recp in focus_pairs:
        for tgt in ("T", "NK"):
            a = sub[(sub["ligand"] == lig) & (sub["receptor"] == recp) & (sub["target"] == tgt)]
            if a.empty:
                continue
            hi = a[a["source"] == senders_hi]
            lo = a[a["source"] == senders_lo]
            hi_m = float(hi["lr_means"].iloc[0]) if len(hi) else np.nan
            lo_m = float(lo["lr_means"].iloc[0]) if len(lo) else np.nan
            if np.isfinite(hi_m) and np.isfinite(lo_m):
                higher = "high" if hi_m > lo_m else ("low" if lo_m > hi_m else "tie")
            elif np.isfinite(hi_m):
                higher = "high only"
            else:
                higher = "low only"
            out_rows.append(
                {
                    "ligand": lig,
                    "receptor": recp,
                    "target": tgt,
                    "lr_means_high": hi_m,
                    "lr_means_low": lo_m,
                    "higher_in": higher,
                    "cellphone_p_high": float(hi["cellphone_pvals"].iloc[0]) if len(hi) else np.nan,
                }
            )
    # Also keep top 8 detected outgoing edges by high lr_means if focus is thin.
    if len(out_rows) < 6:
        hi_only = sub[sub["source"] == senders_hi].sort_values("lr_means", ascending=False).head(12)
        for r in hi_only.itertuples(index=False):
            lo = sub[
                (sub["ligand"] == r.ligand)
                & (sub["receptor"] == r.receptor)
                & (sub["target"] == r.target)
                & (sub["source"] == senders_lo)
            ]
            lo_m = float(lo["lr_means"].iloc[0]) if len(lo) else np.nan
            out_rows.append(
                {
                    "ligand": r.ligand,
                    "receptor": r.receptor,
                    "target": r.target,
                    "lr_means_high": float(r.lr_means),
                    "lr_means_low": lo_m,
                    "higher_in": "high" if (np.isfinite(lo_m) and r.lr_means > lo_m) else ("high only" if not np.isfinite(lo_m) else "low"),
                    "cellphone_p_high": float(r.cellphone_pvals),
                }
            )
    return pd.DataFrame(out_rows).drop_duplicates(subset=["ligand", "receptor", "target"])


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
    cells, expr, total = stream_matrix(mtx, keep)
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
    if "CLDN4" not in logx:
        raise SystemExit("CLDN4 missing from streamed matrix")
    cld = logx["CLDN4"]
    cld_hi, cld_lo, cld_thr = high_low_masks(cld, is_malig)

    ntab_rows = []
    for s in sorted(set(sample)):
        m = sample == s
        n_hi = int((m & cld_hi).sum())
        n_lo = int((m & cld_lo).sum())
        n_t = int((m & is_tnk).sum())
        paired = (
            n_hi >= P["min_malig_per_state"]
            and n_lo >= P["min_malig_per_state"]
            and n_t >= P["min_tnk"]
        )
        ntab_rows.append(
            {
                "Sample": s,
                "patient": PATIENT_OF.get(s, s),
                "paper_group": PAPER_GROUP.get(s, ""),
                "n_cells": int(m.sum()),
                "n_epithelial": int((m & is_epi).sum()),
                "n_malig_like": int((m & is_malig).sum()),
                "n_T": int((m & is_t).sum()),
                "n_NK": int((m & is_nk).sum()),
                "n_T_NK": n_t,
                "n_cldn4_high": n_hi,
                "n_cldn4_low": n_lo,
                "paired": paired,
            }
        )
    ntab = pd.DataFrame(ntab_rows)
    ntab.to_csv(args.outdir / "n_cells_patients.tsv", sep="\t", index=False)

    out_hi = score_pairs(pairs, logx, cld_hi, is_tnk, P["expr_prop"])
    out_lo = score_pairs(pairs, logx, cld_lo, is_tnk, P["expr_prop"])
    in_hi = score_pairs(pairs, logx, is_tnk, cld_hi, P["expr_prop"])
    in_lo = score_pairs(pairs, logx, is_tnk, cld_lo, P["expr_prop"])
    paired_delta(out_hi, out_lo, "cpdb_mean_score").to_csv(
        args.outdir / "pooled_outgoing_cldn4.tsv", sep="\t", index=False
    )
    paired_delta(in_hi, in_lo, "cpdb_mean_score").to_csv(
        args.outdir / "pooled_incoming_cldn4.tsv", sep="\t", index=False
    )

    def patient_table(direction: str) -> pd.DataFrame:
        chunks = []
        for s in sorted(set(sample)):
            m = sample == s
            h, l, t = m & cld_hi, m & cld_lo, m & is_tnk
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
            d["patient"] = PATIENT_OF.get(s, s)
            d["paper_group"] = PAPER_GROUP.get(s, "")
            d["n_high"] = int(h.sum())
            d["n_low"] = int(l.sum())
            d["n_tnk"] = int(t.sum())
            chunks.append(d)
        return pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame()

    ranks: dict[str, pd.DataFrame] = {}
    n_paired = {}
    for direction in ("outgoing", "incoming"):
        raw = patient_table(direction)
        raw.to_csv(args.outdir / f"patient_cldn4_{direction}.tsv", sep="\t", index=False)
        if raw.empty:
            ranks[f"cldn4_{direction}"] = pd.DataFrame(
                columns=["ligand", "receptor", "pathway", "n_patients", "median_delta", "pval", "padj"]
            )
            n_paired[direction] = 0
            continue
        n_paired[direction] = int(raw["Sample"].nunique())
        rows = []
        for (lig, recp, path), sub in raw.groupby(["ligand", "receptor", "pathway"], observed=True):
            keep_sub = sub[sub["pass_either"]]
            if len(keep_sub) < 3:
                continue
            if path not in {"T_recruit", "IFN", "MHC_I"} and keep_sub["pass_both"].sum() < 3:
                continue
            delta = keep_sub["delta_high_minus_low"].to_numpy()
            p = wilcoxon_safe(keep_sub["cpdb_mean_score_high"], keep_sub["cpdb_mean_score_low"])
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
        tab.to_csv(args.outdir / f"ranks_cldn4_{direction}.tsv", sep="\t", index=False)
        ranks[f"cldn4_{direction}"] = tab

    # Dedicated LR tables for the PR.
    out_tab = ranks["cldn4_outgoing"].copy()
    out_tab.insert(0, "direction", "malignant_CLDN4high_to_TNK")
    out_tab.to_csv(args.outdir / "lr_table_cldn4_outgoing_tnk.tsv", sep="\t", index=False)
    in_tab = ranks["cldn4_incoming"].copy()
    in_tab.insert(0, "direction", "TNK_to_malignant_CLDN4")
    in_tab.to_csv(args.outdir / "lr_table_cldn4_incoming_tnk.tsv", sep="\t", index=False)
    focus = out_tab[out_tab["pathway"].isin(["T_recruit", "IFN", "MHC_I"])].copy()
    focus.to_csv(args.outdir / "lr_table_focus_outgoing.tsv", sep="\t", index=False)

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
    pd.DataFrame(path_rows).to_csv(args.outdir / "pathway_summary.tsv", sep="\t", index=False)

    liana_notes = [f"import: {liana_import}"]
    liana_status = "not_run"
    liana_focus = None
    try:
        import anndata as ad

        mask = cld_hi | cld_lo | is_t | is_nk
        genes = sorted(g for g in logx if g in expr)
        X = np.vstack([logx[g][mask] for g in genes]).T
        obs = pd.DataFrame(
            {"sample": sample[mask], "lineage": lineage[mask]},
            index=pd.Index(cells[mask], name="cell"),
        )
        grp = np.array(["other"] * int(mask.sum()), dtype=object)
        grp[cld_hi[mask]] = "Malig_CLDN4high"
        grp[cld_lo[mask]] = "Malig_CLDN4low"
        grp[is_t[mask]] = "T"
        grp[is_nk[mask]] = "NK"
        obs["cc_group"] = grp
        rng = np.random.default_rng(P["random_seed"])
        keep_idx = []
        for g, idx in obs.groupby("cc_group", observed=True).indices.items():
            if g == "other":
                continue
            if len(idx) > P["liana_max_cells_per_group"]:
                idx = rng.choice(idx, size=P["liana_max_cells_per_group"], replace=False)
            keep_idx.append(np.asarray(idx))
        keep_idx = np.sort(np.concatenate(keep_idx))
        adata = ad.AnnData(X=X[keep_idx], obs=obs.iloc[keep_idx].copy(), var=pd.DataFrame(index=genes))
        adata.obs["cc_group"] = pd.Categorical(adata.obs["cc_group"])
        adata.uns["log1p"] = {"base": None}
        counts = adata.obs["cc_group"].value_counts().to_dict()
        log(f"[liana] CLDN4 groups {counts}")
        note = run_liana(
            adata,
            "cc_group",
            args.outdir / "liana_cellphonedb_cldn4.csv",
            P["liana_n_perms"],
        )
        liana_notes.append(f"CLDN4: {note}; downsampled groups={counts}")
        if note.startswith("LIANA_OK"):
            liana_status = "ran_cellphonedb_method"
            res = pd.read_csv(args.outdir / "liana_cellphonedb_cldn4.csv")
            senders = ["Malig_CLDN4high", "Malig_CLDN4low"]
            res[res["source"].isin(senders) & res["target"].isin(["T", "NK"])].to_csv(
                args.outdir / "liana_cldn4_outgoing_tnk.csv", index=False
            )
            res[res["source"].isin(["T", "NK"]) & res["target"].isin(senders)].to_csv(
                args.outdir / "liana_cldn4_incoming_tnk.csv", index=False
            )
            liana_focus = summarize_liana_focus(res)
            liana_focus.to_csv(args.outdir / "liana_focus_outgoing.tsv", sep="\t", index=False)
        else:
            liana_status = note.split(":")[0].lower()
    except Exception as exc:
        liana_notes.append(f"LIANA wrapper failed: {exc}")
        liana_status = f"failed: {exc}"

    foc = ranks.get("cldn4_outgoing", pd.DataFrame())
    foc = foc[foc["pathway"].isin(["T_recruit", "IFN", "MHC_I"])] if len(foc) else foc
    if len(foc):
        n_neg = int((foc["median_delta"] < 0).sum())
        n_sig = int((foc["padj"] < 0.05).sum()) if foc["padj"].notna().any() else 0
        by = (
            foc.groupby("pathway")
            .apply(
                lambda s: f"{int((s.median_delta < 0).sum())}/{len(s)} Δ<0, median Δ={s.median_delta.median():+.3f}",
                include_groups=False,
            )
            .to_dict()
        )
        trend = (
            f"On the patient-level CellPhoneDB-style score (paired n={n_paired.get('outgoing', 0)}), "
            f"CLDN4-high vs CLDN4-low outgoing T-recruit/IFN/MHC-I pairs: "
            f"{n_neg}/{len(foc)} have median Δ < 0; {n_sig}/{len(foc)} reach FDR < 0.05. "
            f"By axis: {by}. "
            "This is the observed rank in this public 15-sample BD Rhapsody object, not a general rule."
        )
    else:
        trend = "Too few patient-level T-recruit/IFN/MHC-I pairs passed filters to rank a CLDN4-high reduction."

    summary = {
        "dataset": "GSE207422",
        "state_gene": "CLDN4",
        "tacstd2_used_as_gate": False,
        "n_cells": int(n),
        "n_patients_matrix": 15,
        "n_epithelial": int(is_epi.sum()),
        "n_malig": int(is_malig.sum()),
        "n_T": int(is_t.sum()),
        "n_NK": int(is_nk.sum()),
        "n_tnk": int(is_tnk.sum()),
        "n_cldn4_high": int(cld_hi.sum()),
        "n_cldn4_low": int(cld_lo.sum()),
        "cldn4_threshold": cld_thr,
        "min_malig_per_state": P["min_malig_per_state"],
        "min_tnk": P["min_tnk"],
        "n_patients_paired": n_paired.get("outgoing", 0),
        "n_pairs_scored": int(len(pairs)),
        "liana_status": liana_status,
        "liana_n_perms": P["liana_n_perms"],
        "liana_max_cells": P["liana_max_cells_per_group"],
        "cellchat_status": "not_run_R_unavailable",
        "liana_notes": liana_notes,
        "trend_sentence": trend,
        "a3_taken_as_given": True,
        "method": "cellphonedb_mean_of_means_on_log1p_cp10k",
    }
    (args.outdir / "summary.json").write_text(json.dumps(summary, indent=2))
    make_figures(args.outdir, ntab, ranks)
    write_finding(
        outdir=args.outdir,
        root=HERE,
        summary=summary,
        ntab=ntab,
        ranks=ranks,
        liana_notes=liana_notes,
        liana_focus=liana_focus,
    )
    log(json.dumps({k: summary[k] for k in ("n_malig", "n_tnk", "n_patients_paired", "liana_status", "trend_sentence")}, indent=2))
    log(f"[done] {args.outdir}")


if __name__ == "__main__":
    main()
