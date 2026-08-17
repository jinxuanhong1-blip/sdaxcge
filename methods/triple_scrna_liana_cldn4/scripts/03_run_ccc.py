#!/usr/bin/env python3
"""Patient-level CLDN4-only CellPhoneDB-style LR on the triple merge.

Primary: documented mean-of-means on log1p(CP10k), paired Wilcoxon.
Secondary: LIANA cellphonedb if importable. CellChat is not run.
"""

from __future__ import annotations

import argparse
import json
import traceback
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from scipy import stats
from statsmodels.stats.multitest import multipletests

HERE = Path(__file__).resolve().parents[1]


def log(msg: str) -> None:
    print(msg, flush=True)


def group_gene_stats(mat: np.ndarray, genes: list[str], mask: np.ndarray) -> tuple[dict[str, float], dict[str, float], int]:
    n = int(mask.sum())
    means, fracs = {}, {}
    if n == 0:
        return means, fracs, 0
    sub = mat[mask]
    for j, g in enumerate(genes):
        v = sub[:, j]
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


def score_pairs(pairs: pd.DataFrame, mat, genes, sender, receiver, expr_prop: float) -> pd.DataFrame:
    s_mean, s_frac, n_s = group_gene_stats(mat, genes, sender)
    r_mean, r_frac, n_r = group_gene_stats(mat, genes, receiver)
    rows = []
    for rec in pairs.itertuples(index=False):
        lig_u = str(rec.ligand).split("+")
        rec_u = str(rec.receptor).split("+")
        l_mean, l_frac = partner_from_stats(lig_u, s_mean, s_frac)
        rec_m, rec_f = partner_from_stats(rec_u, r_mean, r_frac)
        if not np.isfinite(l_mean) or not np.isfinite(rec_m):
            continue
        rows.append(
            {
                "ligand": rec.ligand,
                "receptor": rec.receptor,
                "pathway": rec.pathway,
                "pair_origin": rec.pair_origin,
                "n_sender": n_s,
                "n_receiver": n_r,
                "ligand_mean": l_mean,
                "receptor_mean": rec_m,
                "ligand_frac": l_frac,
                "receptor_frac": rec_f,
                "cpdb_mean_score": 0.5 * (l_mean + rec_m),
                "pass_expr_prop": (l_frac >= expr_prop) and (rec_f >= expr_prop),
            }
        )
    return pd.DataFrame(rows)


def paired_delta(high_df: pd.DataFrame, low_df: pd.DataFrame) -> pd.DataFrame:
    key = ["ligand", "receptor", "pathway"]
    cols = key + ["cpdb_mean_score", "pass_expr_prop", "ligand_frac", "receptor_frac", "n_sender", "n_receiver"]
    m = high_df[cols].merge(low_df[cols], on=key, suffixes=("_high", "_low"))
    m["delta_high_minus_low"] = m["cpdb_mean_score_high"] - m["cpdb_mean_score_low"]
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


def fmt_p(x) -> str:
    if x is None or not np.isfinite(x):
        return "NA"
    return f"{x:.3g}"


def run_liana(adata, groupby: str, out_csv: Path, n_perms: int) -> str:
    try:
        import liana as li
    except Exception as exc:
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


def summarize_liana_focus(res: pd.DataFrame) -> pd.DataFrame:
    senders = ["Malig_CLDN4high", "Malig_CLDN4low"]
    sub = res[res["source"].isin(senders) & res["target"].isin(["T", "NK", "TNK"])].copy()
    if sub.empty:
        return sub
    focus_lig = {
        "CXCL9",
        "CXCL10",
        "CXCL11",
        "CXCL16",
        "CX3CL1",
        "HLA-A",
        "HLA-B",
        "HLA-C",
        "HLA-E",
        "FAM3C",
        "APP",
        "CD58",
    }
    if "ligand_complex" in sub.columns:
        sub = sub[sub["ligand_complex"].isin(focus_lig) | sub["ligand_complex"].astype(str).str.startswith("HLA")]
    return sub


def make_figures(outdir: Path, ntab: pd.DataFrame, ranks: dict[str, pd.DataFrame], cells: pd.DataFrame, thr: float) -> None:
    figdir = outdir / "figures"
    figdir.mkdir(parents=True, exist_ok=True)
    colors = {"GSE131907": "#4c6a92", "GSE148071": "#b23a48", "GSE205335": "#2a9d8f"}

    # 1. Honest n by dataset
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ds = ["GSE131907", "GSE148071", "GSE205335"]
    deposited = [int((ntab.dataset == d).sum()) for d in ds]
    paired = [int(((ntab.dataset == d) & ntab.paired).sum()) for d in ds]
    x = np.arange(len(ds))
    ax.bar(x - 0.18, deposited, 0.36, label="deposited patients in extract", color="#c5c5c5")
    ax.bar(x + 0.18, paired, 0.36, label="paired (both bins + T/NK)", color="#4c6a92")
    ax.set_xticks(x)
    ax.set_xticklabels(ds)
    ax.set_ylabel("patients")
    ax.set_title("Honest paired n is not deposited n")
    ax.legend(frameon=False)
    for i, (a, b) in enumerate(zip(deposited, paired)):
        ax.text(i - 0.18, a + 0.3, str(a), ha="center", fontsize=8)
        ax.text(i + 0.18, b + 0.3, str(b), ha="center", fontsize=8)
    fig.tight_layout()
    fig.savefig(figdir / "fig_honest_n.png", dpi=160)
    fig.savefig(figdir / "fig_honest_n.pdf")
    plt.close(fig)

    # 2. Per-patient cell counts
    show = ntab.sort_values(["dataset", "patient"]).reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(max(8, 0.22 * len(show) + 2), 4.8))
    x = np.arange(len(show))
    ax.bar(x, show["n_cldn4_high"], color="#b23a48", label="CLDN4-high malig")
    ax.bar(x, show["n_cldn4_low"], bottom=show["n_cldn4_high"], color="#e9c46a", label="CLDN4-low malig")
    ax.plot(x, show["n_T_NK"], color="#264653", marker="o", ms=2, lw=0.8, label="T/NK")
    ax.set_xticks(x)
    ax.set_xticklabels([p.split("|", 1)[-1] for p in show["patient"]], rotation=90, fontsize=6)
    ax.set_ylabel("cells")
    ax.set_title("Malignant CLDN4 bins and T/NK per patient")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(figdir / "fig_n_per_patient.png", dpi=160)
    fig.savefig(figdir / "fig_n_per_patient.pdf")
    plt.close(fig)

    # 3. CLDN4 density by dataset
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    mal = cells[cells["is_malig"]]
    for d, sub in mal.groupby("dataset"):
        ax.hist(sub["cldn4_log1p"], bins=40, density=True, histtype="step", lw=1.6, color=colors.get(d, "grey"), label=f"{d} n={len(sub)}")
    ax.axvline(thr, color="black", ls="--", lw=1, label=f"global median {thr:.3f}")
    ax.set_xlabel("malignant CLDN4 log1p(CP10k)")
    ax.set_ylabel("density")
    ax.set_title("CLDN4 among malignant cells — one global cut")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(figdir / "fig_cldn4_density.png", dpi=160)
    fig.savefig(figdir / "fig_cldn4_density.pdf")
    plt.close(fig)

    # 4–5. Focus outgoing / incoming
    for key, title, fname in (
        ("cldn4_outgoing", "Outgoing CLDN4-high malignant → T/NK", "fig_outgoing_focus_delta"),
        ("cldn4_incoming", "Incoming T/NK → CLDN4-high vs low malignant", "fig_incoming_focus_delta"),
    ):
        tab = ranks.get(key, pd.DataFrame())
        foc = tab[tab["pathway"].isin(["T_recruit", "IFN", "MHC_I"])].copy() if len(tab) else tab
        fig, ax = plt.subplots(figsize=(8.2, max(3.5, 0.32 * max(len(foc), 1) + 1.5)))
        if foc.empty:
            ax.text(0.5, 0.5, "no focus pairs", ha="center")
        else:
            foc = foc.sort_values(["pathway", "median_delta"])
            y = np.arange(len(foc))
            cols = [{"T_recruit": "#2a9d8f", "IFN": "#e76f51", "MHC_I": "#4c6a92"}.get(p, "grey") for p in foc["pathway"]]
            ax.barh(y, foc["median_delta"], color=cols)
            ax.axvline(0, color="black", lw=0.8)
            ax.set_yticks(y)
            ax.set_yticklabels([f"{a}–{b} (n={int(n)})" for a, b, n in zip(foc["ligand"], foc["receptor"], foc["n_patients"])], fontsize=8)
            ax.set_xlabel("median patient Δ (high − low)")
        ax.set_title(title)
        fig.tight_layout()
        fig.savefig(figdir / f"{fname}.png", dpi=160)
        fig.savefig(figdir / f"{fname}.pdf")
        plt.close(fig)

    # 6. Extra ligand table figure
    out = ranks.get("cldn4_outgoing", pd.DataFrame())
    foc = out[out["pathway"].isin(["T_recruit", "IFN", "MHC_I"])].copy() if len(out) else out
    fig, ax = plt.subplots(figsize=(10.5, max(3.2, 0.38 * max(len(foc), 1) + 1.8)))
    ax.axis("off")
    ax.set_title("Primary LR table — outgoing CLDN4-high malignant → T/NK", pad=12)
    if foc.empty:
        ax.text(0.5, 0.5, "no pairs", ha="center")
    else:
        foc = foc.sort_values(["pathway", "median_delta"])
        tbl = [
            [
                r.pathway,
                f"{r.ligand}–{r.receptor}",
                str(int(r.n_patients)),
                f"{r.median_delta:+.3f}",
                fmt_p(r.pval),
                fmt_p(r.padj),
            ]
            for r in foc.itertuples(index=False)
        ]
        table = ax.table(
            cellText=tbl,
            colLabels=["Pathway", "Pair", "n", "median Δ", "p", "FDR"],
            loc="center",
            cellLoc="center",
        )
        table.auto_set_font_size(False)
        table.set_fontsize(8)
        table.scale(1, 1.25)
    fig.tight_layout()
    fig.savefig(figdir / "fig_extra_ligand_table.png", dpi=160)
    fig.savefig(figdir / "fig_extra_ligand_table.pdf")
    plt.close(fig)

    # 7. Dataset contribution among paired patients
    fig, ax = plt.subplots(figsize=(5.6, 4.2))
    paired = ntab[ntab["paired"]]
    counts = paired["dataset"].value_counts().reindex(ds).fillna(0)
    ax.bar(ds, counts.to_numpy(), color=[colors[d] for d in ds])
    ax.set_ylabel("paired patients")
    ax.set_title(f"Paired n={int(paired.shape[0])} by dataset")
    for i, v in enumerate(counts.to_numpy()):
        ax.text(i, v + 0.15, str(int(v)), ha="center")
    fig.tight_layout()
    fig.savefig(figdir / "fig_dataset_contribution.png", dpi=160)
    fig.savefig(figdir / "fig_dataset_contribution.pdf")
    plt.close(fig)

    # 8. Top outgoing |delta| among all pathways that passed
    tab = ranks.get("cldn4_outgoing", pd.DataFrame())
    fig, ax = plt.subplots(figsize=(8.0, 5.2))
    if tab.empty:
        ax.text(0.5, 0.5, "no pairs", ha="center")
    else:
        top = tab.reindex(tab["median_delta"].abs().sort_values(ascending=False).index).head(20)
        top = top.sort_values("median_delta")
        y = np.arange(len(top))
        ax.barh(y, top["median_delta"], color=["#b23a48" if v < 0 else "#4c6a92" for v in top["median_delta"]])
        ax.axvline(0, color="black", lw=0.8)
        ax.set_yticks(y)
        ax.set_yticklabels([f"{a}–{b}" for a, b in zip(top["ligand"], top["receptor"])], fontsize=8)
        ax.set_xlabel("median patient Δ (high − low)")
        ax.set_title("Top |Δ| outgoing pairs (any pathway)")
    fig.tight_layout()
    fig.savefig(figdir / "fig_top_outgoing_absdelta.png", dpi=160)
    fig.savefig(figdir / "fig_top_outgoing_absdelta.pdf")
    plt.close(fig)


def write_finding(root: Path, outdir: Path, summary: dict, ntab: pd.DataFrame, ranks: dict, liana_notes: list[str], liana_focus: pd.DataFrame | None) -> None:
    out = ranks["cldn4_outgoing"]
    inc = ranks["cldn4_incoming"]
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

    def _plist(df: pd.DataFrame) -> str:
        bits = []
        for r in df.itertuples(index=False):
            bits.append(
                f"{r.patient} ({r.dataset}; high={int(r.n_cldn4_high)}/low={int(r.n_cldn4_low)} malig, T/NK={int(r.n_T_NK)})"
            )
        return ", ".join(bits) if bits else "(none)"

    by_ds = (
        ntab.groupby("dataset")
        .agg(deposited=("patient", "nunique"), paired=("paired", "sum"), malig=("n_malig", "sum"), tnk=("n_T_NK", "sum"))
        .reset_index()
    )
    ds_lines = ["| Dataset | Deposited in extract | Paired | Malignant cells | T/NK |", "|---|---:|---:|---:|---:|"]
    for r in by_ds.itertuples(index=False):
        ds_lines.append(f"| {r.dataset} | {int(r.deposited)} | {int(r.paired)} | {int(r.malig)} | {int(r.tnk)} |")

    focus_out = out[out["pathway"].isin(["T_recruit", "IFN", "MHC_I"])].copy() if len(out) else out
    focus_in = inc[inc["pathway"].isin(["T_recruit", "IFN", "MHC_I"])].copy() if len(inc) else inc
    n_paired = int(summary["n_patients_paired"])

    liana_block = ["## LIANA CellPhoneDB method (secondary, pooled / downsampled)", ""]
    for note in liana_notes:
        liana_block.append(f"- {note}")
    liana_block.append("")
    if liana_focus is not None and len(liana_focus):
        liana_block += [
            "Focus edges that cleared LIANA `expr_prop=0.10` from malignant senders to T/NK:",
            "",
            "| Pair | target | high `lr_means` | low `lr_means` | higher in |",
            "|---|---|---:|---:|---|",
        ]
        # try to pair high/low
        if {"source", "target", "ligand_complex", "receptor_complex", "lr_means"}.issubset(liana_focus.columns):
            hi = liana_focus[liana_focus["source"] == "Malig_CLDN4high"]
            lo = liana_focus[liana_focus["source"] == "Malig_CLDN4low"]
            key = ["ligand_complex", "receptor_complex", "target"]
            m = hi.merge(lo, on=key, suffixes=("_high", "_low"), how="outer")
            for r in m.head(12).itertuples(index=False):
                pair = f"{getattr(r, 'ligand_complex')}–{getattr(r, 'receptor_complex')}"
                hv = getattr(r, "lr_means_high", np.nan)
                lv = getattr(r, "lr_means_low", np.nan)
                winner = "high" if (np.isfinite(hv) and np.isfinite(lv) and hv >= lv) else "low"
                liana_block.append(
                    f"| {pair} | {r.target} | {hv if np.isfinite(hv) else float('nan'):.3f} | "
                    f"{lv if np.isfinite(lv) else float('nan'):.3f} | {winner} |"
                )
        liana_block.append("")
        liana_block.append("LIANA p-values are within-object specificity, not patient-level tests.")
        liana_block.append("")

    lines = [
        "# Triple merge GSE131907+GSE148071+GSE205335 — CLDN4-only LIANA/CellPhoneDB",
        "",
        "**CLDN4 only.** TACSTD2 is not a gate and is not used to define dual-high.",
        "**Not GSE207422.** **Not** the 131907+205335-only pair merge.",
        "",
        f"**Verdict (paired n={n_paired} honest):** {summary['trend_sentence']}",
        "",
        "## What was run",
        "",
        "| Method | Status |",
        "|---|---|",
        f"| Documented CellPhoneDB-style score (mean of partner means on log1p CP10k; {summary['n_pairs_scored']} pairs) | **primary** — patient-level paired Wilcoxon |",
        f"| LIANA `mt.cellphonedb` (resource `cellphonedb`, 50 permutations, ≤2,000 cells/group) | **secondary** — {summary['liana_status']} |",
        "| CellChat | **not run** — R unavailable; no CellChat tables were written |",
        "",
        "## Honest n",
        "",
        "| Set | n |",
        "|---|---:|",
        f"| Cells kept (malignant or T/NK) | {summary['n_cells_kept']} |",
        f"| Malignant / malignant-like | {summary['n_malig']} |",
        f"| T / NK / T+NK | {summary['n_T']} / {summary['n_NK']} / {summary['n_tnk']} |",
        f"| CLDN4-high / low malignant | {summary['n_cldn4_high']} / {summary['n_cldn4_low']} |",
        f"| Patients in the extract | {summary['n_patients_extract']} |",
        f"| Patients in the paired LR test | **{n_paired}** |",
        "",
        *ds_lines,
        "",
        f"CLDN4 high = at or above the **global median** log1p(CP10k) among malignant cells in the merge (threshold = {summary['cldn4_threshold']:.3f}). "
        "Low = below. A patient enters the paired test if it has ≥10 malignant cells in **both** bins and ≥20 T/NK cells.",
        "",
        f"**Paired patients (n={n_paired}):** {_plist(paired)}.",
        "",
        f"**Dropped (empty or one-sided malignant bins, or thin T/NK):** {_plist(dropped)}.",
        "",
        "The header n for the Wilcoxon is the paired count, not deposited n and not the number of GEO samples.",
        "Per-sample counts: `results/n_cells_patients.tsv`.",
        "",
        "## Score",
        "",
        "On log1p(CP10k), each partner’s expression is the **minimum subunit mean** (CellPhoneDB complex rule). "
        "The pair score is the **mean of the two partner means** (Efremova et al. 2020 *Nat Protoc*; Garcia-Alonso et al. 2022 *Nat Protoc*). "
        "A pair is flagged `pass_expr_prop` when both partners are detected in ≥10% of cells in their group. "
        "Patient-level tests use that patient’s own T/NK and that patient’s CLDN4-high vs CLDN4-low malignant cells. "
        "Cells are not treated as replicates. FDR is Benjamini–Hochberg within each contrast.",
        "",
        "This is **not** a CellChat communication probability.",
        "",
        "## Primary LR table — outgoing CLDN4-high malignant → T/NK",
        "",
        "Median patient Δ = high − low. Negative = weaker from the CLDN4-high state.",
        "",
        *_table(focus_out),
        "Full ranked table (all pathways that passed filters): `results/lr_table_cldn4_outgoing_tnk.tsv`.",
        "",
        "## Incoming T/NK → CLDN4-high vs CLDN4-low malignant (secondary)",
        "",
        *_table(focus_in),
        *liana_block,
        "## Readout",
        "",
        summary["trend_sentence"],
        "",
        "## Honest limits",
        "",
        f"1. Paired n = {n_paired}. Deposited patients in the extract = {summary['n_patients_extract']}. One-sided or empty malignant bins are reported, not patched.",
        "2. CLDN4 only. Dual-high (TACSTD2 AND CLDN4) was not run. TACSTD2 is not a gate.",
        "3. GSE207422 is excluded. This is not the 131907+205335-only pair.",
        "4. Labels are mixed by design: author malignant/T/NK in GSE131907 and GSE205335; A3 marker-argmax malignant-like in GSE148071 (no GEO labels).",
        "5. GSE131907 uses tumor-site cells only. Normal-lung epithelium is not counted as malignant.",
        "6. Ambient RNA cannot be re-estimated from the processed matrices.",
        "7. Chemokine dropout is high; read `n patients` and `pass_expr_prop` with the ranks.",
        "8. CellChat was not run. LIANA p-values are not patient-level tests.",
        "9. The three series are different platforms and clinical settings. The merge is additive, not batch-corrected Harmony/scVI.",
        "",
        "## Files",
        "",
        "| File | Role |",
        "|---|---|",
        "| `FINDING.md` | This note |",
        "| `results/lr_table_cldn4_outgoing_tnk.tsv` | Primary LR table (patient-level ranks, outgoing) |",
        "| `results/lr_table_cldn4_incoming_tnk.tsv` | Incoming ranks |",
        "| `results/lr_table_focus_outgoing.tsv` | T-recruit / IFN / MHC-I outgoing subset |",
        "| `results/n_cells_patients.tsv` | Per-patient cell counts and high/low bins |",
        "| `results/patient_cldn4_outgoing.tsv` | Per-patient pair scores |",
        "| `results/summary.json` | Machine-readable n and method flags |",
        "| `results/figures/` | Honest-n, density, focus Δ, extra ligand table |",
        "",
        "## Reproduce",
        "",
        "```bash",
        "cd methods/triple_scrna_liana_cldn4",
        "python3 scripts/00_download.py",
        "python3 scripts/01_build_pairs.py",
        "python3 scripts/02_extract.py",
        "python3 scripts/03_run_ccc.py",
        "```",
        "",
    ]
    text = "\n".join(lines)
    (outdir / "FINDING.md").write_text(text)
    (root / "FINDING.md").write_text(text)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--extract", type=Path, default=HERE / "data" / "extract")
    ap.add_argument("--outdir", type=Path, default=HERE / "results")
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "figures").mkdir(exist_ok=True)

    cfg = yaml.safe_load((HERE / "config" / "gene_sets.yaml").read_text())
    P = cfg["params"]
    pairs = pd.read_csv(HERE / "resources" / "cellphonedb_v5_lr_pairs.tsv", sep="\t")
    cells = pd.read_parquet(args.extract / "cells.parquet")
    blob = np.load(args.extract / "log1p_cp10k.npz", allow_pickle=True)
    mat = blob["log1p"]
    genes = [str(g) for g in blob["genes"]]
    if len(cells) != mat.shape[0]:
        raise SystemExit(f"extract mismatch {len(cells)} cells vs {mat.shape}")

    is_malig = cells["is_malig"].to_numpy()
    is_t = cells["is_t"].to_numpy()
    is_nk = cells["is_nk"].to_numpy()
    is_tnk = cells["is_tnk"].to_numpy()
    cldn = cells["cldn4_log1p"].to_numpy()
    patient = cells["patient"].to_numpy()
    dataset = cells["dataset"].to_numpy()

    thr = float(np.median(cldn[is_malig])) if is_malig.any() else float("nan")
    cld_hi = is_malig & (cldn >= thr)
    cld_lo = is_malig & (cldn < thr)
    log(f"CLDN4 global median={thr:.4f} high={int(cld_hi.sum())} low={int(cld_lo.sum())}")

    rows = []
    for p, sub in cells.groupby("patient", observed=True):
        idx = sub.index.to_numpy()
        hi = int(cld_hi[idx].sum())
        lo = int(cld_lo[idx].sum())
        tnk = int(is_tnk[idx].sum())
        paired = hi >= P["min_malig_per_state"] and lo >= P["min_malig_per_state"] and tnk >= P["min_tnk"]
        rows.append(
            {
                "patient": p,
                "dataset": str(sub["dataset"].iloc[0]),
                "n_cells_kept": int(len(sub)),
                "n_malig": int(sub["is_malig"].sum()),
                "n_cldn4_high": hi,
                "n_cldn4_low": lo,
                "n_T": int(sub["is_t"].sum()),
                "n_NK": int(sub["is_nk"].sum()),
                "n_T_NK": tnk,
                "paired": paired,
                "label_source": str(sub["label_source"].iloc[0]),
            }
        )
    ntab = pd.DataFrame(rows).sort_values(["dataset", "patient"])
    ntab.to_csv(args.outdir / "n_cells_patients.tsv", sep="\t", index=False)

    def patient_table(direction: str) -> pd.DataFrame:
        chunks = []
        for p in ntab.loc[ntab["paired"], "patient"]:
            m = patient == p
            h, l, t = m & cld_hi, m & cld_lo, m & is_tnk
            if direction == "outgoing":
                hdf = score_pairs(pairs, mat, genes, h, t, P["expr_prop"])
                ldf = score_pairs(pairs, mat, genes, l, t, P["expr_prop"])
            else:
                hdf = score_pairs(pairs, mat, genes, t, h, P["expr_prop"])
                ldf = score_pairs(pairs, mat, genes, t, l, P["expr_prop"])
            d = paired_delta(hdf, ldf)
            d["patient"] = p
            d["dataset"] = str(ntab.loc[ntab["patient"] == p, "dataset"].iloc[0])
            d["n_high"] = int(h.sum())
            d["n_low"] = int(l.sum())
            d["n_tnk"] = int(t.sum())
            chunks.append(d)
            log(f"[{direction}] {p} high={int(h.sum())} low={int(l.sum())} tnk={int(t.sum())} pairs={len(d)}")
        return pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame()

    ranks: dict[str, pd.DataFrame] = {}
    n_paired: dict[str, int] = {}
    for direction in ("outgoing", "incoming"):
        raw = patient_table(direction)
        raw.to_csv(args.outdir / f"patient_cldn4_{direction}.tsv", sep="\t", index=False)
        if raw.empty:
            ranks[f"cldn4_{direction}"] = pd.DataFrame(
                columns=["ligand", "receptor", "pathway", "n_patients", "median_delta", "pval", "padj"]
            )
            n_paired[direction] = 0
            continue
        n_paired[direction] = int(raw["patient"].nunique())
        rows = []
        for (lig, recp, path), sub in raw.groupby(["ligand", "receptor", "pathway"], observed=True):
            keep_sub = sub[sub["pass_either"]]
            if len(keep_sub) < 3:
                continue
            if path not in {"T_recruit", "IFN", "MHC_I"} and keep_sub["pass_both"].sum() < 3:
                continue
            delta = keep_sub["delta_high_minus_low"].to_numpy()
            pval = wilcoxon_safe(keep_sub["cpdb_mean_score_high"], keep_sub["cpdb_mean_score_low"])
            rows.append(
                {
                    "ligand": lig,
                    "receptor": recp,
                    "pathway": path,
                    "n_patients": int(keep_sub["patient"].nunique()),
                    "n_GSE131907": int(keep_sub.loc[keep_sub["dataset"] == "GSE131907", "patient"].nunique()),
                    "n_GSE148071": int(keep_sub.loc[keep_sub["dataset"] == "GSE148071", "patient"].nunique()),
                    "n_GSE205335": int(keep_sub.loc[keep_sub["dataset"] == "GSE205335", "patient"].nunique()),
                    "median_delta": float(np.median(delta)),
                    "mean_delta": float(np.mean(delta)),
                    "mean_score_high": float(keep_sub["cpdb_mean_score_high"].mean()),
                    "mean_score_low": float(keep_sub["cpdb_mean_score_low"].mean()),
                    "frac_pass_high": float(keep_sub["pass_expr_prop_high"].mean()),
                    "frac_pass_low": float(keep_sub["pass_expr_prop_low"].mean()),
                    "pval": pval,
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

    out_tab = ranks["cldn4_outgoing"].copy()
    out_tab.insert(0, "direction", "malignant_CLDN4high_to_TNK")
    out_tab.to_csv(args.outdir / "lr_table_cldn4_outgoing_tnk.tsv", sep="\t", index=False)
    in_tab = ranks["cldn4_incoming"].copy()
    in_tab.insert(0, "direction", "TNK_to_malignant_CLDN4")
    in_tab.to_csv(args.outdir / "lr_table_cldn4_incoming_tnk.tsv", sep="\t", index=False)
    focus = out_tab[out_tab["pathway"].isin(["T_recruit", "IFN", "MHC_I"])].copy()
    focus.to_csv(args.outdir / "lr_table_focus_outgoing.tsv", sep="\t", index=False)

    # Pooled descriptive scores
    hdf = score_pairs(pairs, mat, genes, cld_hi, is_tnk, P["expr_prop"])
    ldf = score_pairs(pairs, mat, genes, cld_lo, is_tnk, P["expr_prop"])
    pooled = paired_delta(hdf, ldf)
    pooled.to_csv(args.outdir / "pooled_outgoing_cldn4.tsv", sep="\t", index=False)

    liana_notes = []
    liana_status = "not_run"
    liana_focus = None
    try:
        import anndata as ad

        mask = cld_hi | cld_lo | is_tnk
        X = mat[mask]
        obs = cells.loc[mask, ["dataset", "patient"]].copy()
        grp = np.array(["other"] * int(mask.sum()), dtype=object)
        grp[cld_hi[mask]] = "Malig_CLDN4high"
        grp[cld_lo[mask]] = "Malig_CLDN4low"
        grp[is_t[mask]] = "T"
        grp[is_nk[mask] & ~is_t[mask]] = "NK"
        # leftover T/NK without T/NK split
        leftover = is_tnk[mask] & (grp == "other")
        grp[leftover] = "T"
        obs["cc_group"] = grp
        rng = np.random.default_rng(P["random_seed"])
        keep_idx = []
        for g, idx in obs.groupby("cc_group", observed=True).indices.items():
            if g == "other":
                continue
            if len(idx) > P["liana_max_cells_per_group"]:
                idx = rng.choice(idx, size=P["liana_max_cells_per_group"], replace=False)
            keep_idx.append(np.asarray(idx))
        keep_idx = np.sort(np.concatenate(keep_idx)) if keep_idx else np.array([], dtype=int)
        adata = ad.AnnData(X=X[keep_idx], obs=obs.iloc[keep_idx].copy(), var=pd.DataFrame(index=genes))
        adata.obs["cc_group"] = pd.Categorical(adata.obs["cc_group"])
        adata.uns["log1p"] = {"base": None}
        counts = adata.obs["cc_group"].value_counts().to_dict()
        log(f"[liana] groups {counts}")
        note = run_liana(adata, "cc_group", args.outdir / "liana_cellphonedb_cldn4.csv", P["liana_n_perms"])
        liana_notes.append(f"CLDN4: {note}; downsampled groups={counts}")
        if note.startswith("LIANA_OK"):
            liana_status = "ran_cellphonedb_method"
            res = pd.read_csv(args.outdir / "liana_cellphonedb_cldn4.csv")
            senders = ["Malig_CLDN4high", "Malig_CLDN4low"]
            res[res["source"].isin(senders) & res["target"].isin(["T", "NK"])].to_csv(
                args.outdir / "liana_cldn4_outgoing_tnk.csv", index=False
            )
            liana_focus = summarize_liana_focus(res)
            if len(liana_focus):
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
            f"CLDN4-high vs CLDN4-low outgoing T-recruit / IFN / MHC-I pairs: "
            f"{n_neg}/{len(foc)} have median Δ < 0; {n_sig}/{len(foc)} reach FDR < 0.05. "
            f"By axis: {by}. "
            "This is the observed rank in this public triple merge, not a general rule."
        )
    else:
        trend = "Too few patient-level T-recruit/IFN/MHC-I pairs passed filters to rank a CLDN4-high reduction."

    summary = {
        "datasets": ["GSE131907", "GSE148071", "GSE205335"],
        "excluded": ["GSE207422"],
        "state_gene": "CLDN4",
        "tacstd2_used_as_gate": False,
        "dual_high": False,
        "n_cells_kept": int(len(cells)),
        "n_malig": int(is_malig.sum()),
        "n_T": int(is_t.sum()),
        "n_NK": int(is_nk.sum()),
        "n_tnk": int(is_tnk.sum()),
        "n_cldn4_high": int(cld_hi.sum()),
        "n_cldn4_low": int(cld_lo.sum()),
        "cldn4_threshold": thr,
        "min_malig_per_state": P["min_malig_per_state"],
        "min_tnk": P["min_tnk"],
        "n_patients_extract": int(ntab["patient"].nunique()),
        "n_patients_paired": n_paired.get("outgoing", 0),
        "n_pairs_scored": int(len(pairs)),
        "liana_status": liana_status,
        "cellchat_status": "not_run_R_unavailable",
        "liana_notes": liana_notes,
        "trend_sentence": trend,
        "method": "cellphonedb_mean_of_means_on_log1p_cp10k",
        "not_the_131907_205335_only_pair": True,
    }
    (args.outdir / "summary.json").write_text(json.dumps(summary, indent=2))
    make_figures(args.outdir, ntab, ranks, cells, thr)
    write_finding(HERE, args.outdir, summary, ntab, ranks, liana_notes, liana_focus)
    log(json.dumps({k: summary[k] for k in ("n_malig", "n_tnk", "n_patients_paired", "liana_status", "trend_sentence")}, indent=2))
    log(f"[done] {args.outdir}")


if __name__ == "__main__":
    main()
