#!/usr/bin/env python3
"""Patient-level Spearman + random-effects meta for TLS/B vs TACSTD2/CLDN4.

Self-contained public-data analysis.
Filters are not tuned to a user-claimed number.
A6/B6 immune-cold is taken as given and is not re-tested here.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from genes import TLS12

ROOT = Path(__file__).resolve().parents[2]
EXTRACTED = ROOT / "data" / "scrna_tls_meta" / "extracted"
OUT = ROOT / "results" / "scrna_tls_meta"

MIN_MAL = 20
MIN_B_INNER = 10
MIN_N_SPEARMAN = 5

PRIMARY = [
    ("tacstd2_mal_mean_log1p_cp10k", "frac_B", "TACSTD2 vs B fraction"),
    ("tacstd2_mal_mean_log1p_cp10k", "tls12_z", "TACSTD2 vs TLS12"),
    ("cldn4_mal_mean_log1p_cp10k", "frac_B", "CLDN4 vs B fraction"),
    ("cldn4_mal_mean_log1p_cp10k", "tls12_z", "CLDN4 vs TLS12"),
]


def log1p_cp10k(umi: np.ndarray, lib: np.ndarray) -> np.ndarray:
    lib = np.asarray(lib, dtype=float)
    umi = np.asarray(umi, dtype=float)
    out = np.full(len(umi), np.nan, dtype=float)
    ok = lib > 0
    out[ok] = np.log1p(umi[ok] * 1e4 / lib[ok])
    return out


def partial_spearman(x, y, z, contrast: str, cohort: str) -> dict:
    """Spearman of residuals after OLS of ranks(x), ranks(y) on ranks(z)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    x, y, z = x[m], y[m], z[m]
    rec = {
        "cohort": cohort,
        "contrast": contrast,
        "n": int(len(x)),
        "spearman_rho": np.nan,
        "spearman_p": np.nan,
        "note": "partial_on_epi_fraction",
    }
    if len(x) < MIN_N_SPEARMAN:
        rec["note"] = "underpowered; inconclusive"
        return rec
    rx, ry, rz = stats.rankdata(x), stats.rankdata(y), stats.rankdata(z)
    A = np.column_stack([np.ones(len(x)), rz])
    bx, *_ = np.linalg.lstsq(A, rx, rcond=None)
    by, *_ = np.linalg.lstsq(A, ry, rcond=None)
    ex, ey = rx - A @ bx, ry - A @ by
    if np.unique(ex).size < 2 or np.unique(ey).size < 2:
        rec["note"] = "no_variance"
        return rec
    rho, p = stats.pearsonr(ex, ey)
    rec["spearman_rho"] = float(rho)
    rec["spearman_p"] = float(p)
    return rec


def spearman(x, y, contrast: str, cohort: str) -> dict:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    n = int(len(x))
    rec = {
        "cohort": cohort,
        "contrast": contrast,
        "n": n,
        "spearman_rho": np.nan,
        "spearman_p": np.nan,
        "note": "",
    }
    if n < MIN_N_SPEARMAN:
        rec["note"] = "underpowered; inconclusive"
        return rec
    if np.unique(x).size < 2 or np.unique(y).size < 2:
        rec["note"] = "no_variance"
        return rec
    rho, p = stats.spearmanr(x, y)
    rec["spearman_rho"] = float(rho)
    rec["spearman_p"] = float(p)
    return rec


def fisher_z(rho: float) -> float:
    r = float(np.clip(rho, -0.999999, 0.999999))
    return float(np.arctanh(r))


def random_effects_meta(rows: list[dict]) -> dict:
    usable = [r for r in rows if r["n"] >= MIN_N_SPEARMAN and np.isfinite(r.get("spearman_rho", np.nan))]
    if len(usable) < 2:
        return {
            "k": len(usable),
            "n_total": int(sum(r["n"] for r in usable)),
            "rho": np.nan,
            "p": np.nan,
            "ci_lo": np.nan,
            "ci_hi": np.nan,
            "tau2": np.nan,
            "note": "fewer than 2 cohorts",
        }
    z = np.array([fisher_z(r["spearman_rho"]) for r in usable])
    n = np.array([r["n"] for r in usable], dtype=float)
    se = 1.0 / np.sqrt(n - 3.0)
    w = 1.0 / se**2
    zbar = np.sum(w * z) / np.sum(w)
    q = np.sum(w * (z - zbar) ** 2)
    k = len(usable)
    c = np.sum(w) - np.sum(w**2) / np.sum(w)
    tau2 = max(0.0, (q - (k - 1.0)) / c) if c > 0 else 0.0
    wstar = 1.0 / (se**2 + tau2)
    zmeta = np.sum(wstar * z) / np.sum(wstar)
    se_meta = math.sqrt(1.0 / np.sum(wstar))
    p = float(2 * stats.norm.sf(abs(zmeta / se_meta)))
    lo = math.tanh(zmeta - 1.96 * se_meta)
    hi = math.tanh(zmeta + 1.96 * se_meta)
    return {
        "k": k,
        "n_total": int(n.sum()),
        "rho": float(math.tanh(zmeta)),
        "p": p,
        "ci_lo": float(lo),
        "ci_hi": float(hi),
        "tau2": float(tau2),
        "Q": float(q),
        "note": "",
        "cohorts": [r["cohort"] for r in usable],
    }


def bh_fdr(pvals: list[float]) -> list[float]:
    p = np.asarray(pvals, dtype=float)
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return q.tolist()
    idx = np.where(ok)[0]
    order = idx[np.argsort(p[idx])]
    m = len(order)
    ranked = p[order]
    adj = ranked * m / np.arange(1, m + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    q[order] = adj
    return q.tolist()


def gene_col(df: pd.DataFrame, gene: str) -> np.ndarray | None:
    if gene in df.columns:
        return df[gene].to_numpy(dtype=float)
    return None


def patient_table(df: pd.DataFrame, cohort: str) -> pd.DataFrame:
    tumor = df[df["tissue"].astype(str).str.lower().eq("tumor")].copy()
    if tumor.empty:
        tumor = df.copy()
        tumor["_note"] = "no_tissue_flag_used_all"
    rows = []
    for pid, g in tumor.groupby("patient"):
        n = len(g)
        lib = g["n_umi"].to_numpy(dtype=float)
        mal = g[g["is_malignant"].astype(bool)]
        epi = g[g["lineage"].astype(str).eq("epithelial") | g["lineage"].astype(str).eq("gate")]
        use = mal if len(mal) >= MIN_MAL else epi
        if len(use) < MIN_MAL:
            compartment = "insufficient"
        elif (use["lineage"].astype(str) == "gate").mean() >= 0.5:
            compartment = "gate"
        elif len(mal) >= MIN_MAL:
            compartment = "malignant"
        else:
            compartment = "epithelial"
        n_b = int((g["lineage"].astype(str) == "B").sum())
        n_plasma = int((g["lineage"].astype(str) == "plasma").sum())
        rec = {
            "cohort": cohort,
            "patient": str(pid),
            "n_cells": n,
            "n_malignant": int(len(mal)),
            "n_epithelial_or_gate": int(len(epi)),
            "n_B": n_b,
            "n_plasma": n_plasma,
            "frac_B": n_b / n if n else np.nan,
            "frac_B_plasma": (n_b + n_plasma) / n if n else np.nan,
            "frac_epithelial": float(((g["lineage"].astype(str) == "epithelial") | (g["lineage"].astype(str) == "gate")).mean()) if n else np.nan,
            "compartment": compartment,
            "eligible": compartment != "insufficient",
        }
        for gene, prefix in (("TACSTD2", "tacstd2"), ("CLDN4", "cldn4")):
            arr = gene_col(use, gene) if len(use) else None
            lib_u = use["n_umi"].to_numpy(dtype=float) if len(use) else np.array([])
            if arr is None or len(use) == 0:
                rec[f"{prefix}_mal_mean_log1p_cp10k"] = np.nan
                rec[f"{prefix}_mal_pct_pos"] = np.nan
                rec[f"{prefix}_mal_pb_cpm"] = np.nan
            else:
                rec[f"{prefix}_mal_mean_log1p_cp10k"] = float(np.nanmean(log1p_cp10k(arr, lib_u)))
                rec[f"{prefix}_mal_pct_pos"] = float(np.mean(arr > 0) * 100.0)
                rec[f"{prefix}_mal_pb_cpm"] = float(arr.sum() / lib_u.sum() * 1e6) if lib_u.sum() > 0 else np.nan
        # TLS / CXCL13 / MS4A1 on all tumor cells
        for gene, key in (("CXCL13", "cxcl13"), ("MS4A1", "ms4a1")):
            arr = gene_col(g, gene)
            if arr is None:
                rec[f"{key}_mean_log1p_cp10k"] = np.nan
                rec[f"{key}_pct_pos"] = np.nan
            else:
                rec[f"{key}_mean_log1p_cp10k"] = float(np.nanmean(log1p_cp10k(arr, lib)))
                rec[f"{key}_pct_pos"] = float(np.mean(arr > 0) * 100.0)
        tcells = g[g["lineage"].astype(str) == "T"]
        if "CXCL13" in g.columns and len(tcells):
            rec["frac_CXCL13pos_T"] = float(np.mean(tcells["CXCL13"].to_numpy(dtype=float) > 0))
        else:
            rec["frac_CXCL13pos_T"] = np.nan
        tls_means = {}
        for gene in TLS12:
            arr = gene_col(g, gene)
            if arr is None:
                continue
            tls_means[gene] = float(np.nanmean(log1p_cp10k(arr, lib)))
        rec["tls12_k"] = len(tls_means)
        rec["_tls12_means"] = tls_means
        rows.append(rec)
    out = pd.DataFrame(rows)
    # within-cohort z of TLS12 gene means, then average
    if out.empty:
        return out
    zmat = []
    genes_used = []
    for gene in TLS12:
        vals = []
        for rec in rows:
            vals.append(rec["_tls12_means"].get(gene, np.nan))
        vals = np.asarray(vals, dtype=float)
        if np.isfinite(vals).sum() >= MIN_N_SPEARMAN and np.nanstd(vals) > 0:
            z = (vals - np.nanmean(vals)) / np.nanstd(vals, ddof=1)
            zmat.append(z)
            genes_used.append(gene)
    if zmat:
        out["tls12_z"] = np.nanmean(np.vstack(zmat), axis=0)
    else:
        out["tls12_z"] = np.nan
    out["tls12_genes"] = ",".join(genes_used)
    out = out.drop(columns=["_tls12_means"])
    return out


def tls_single_cell_table(df: pd.DataFrame, cohort: str) -> pd.DataFrame:
    tumor = df[df["tissue"].astype(str).str.lower().eq("tumor")]
    if tumor.empty:
        tumor = df
    rows = []
    lineages = sorted(tumor["lineage"].astype(str).unique())
    for lin in lineages + ["ALL"]:
        sub = tumor if lin == "ALL" else tumor[tumor["lineage"].astype(str) == lin]
        rec = {
            "cohort": cohort,
            "lineage": lin,
            "n_cells": int(len(sub)),
            "n_patients": int(sub["patient"].nunique()),
        }
        for gene in ["CXCL13", "MS4A1", "CD79A", "CCL19", "CCL21"]:
            if gene in sub.columns:
                rec[f"{gene}_pct_pos"] = float(np.mean(sub[gene].to_numpy(dtype=float) > 0) * 100.0)
            else:
                rec[f"{gene}_pct_pos"] = np.nan
        rec["tls12_genes_present"] = sum(1 for g in TLS12 if g in sub.columns)
        rows.append(rec)
    return pd.DataFrame(rows)


def forest_plot(stats_df: pd.DataFrame, contrast: str, meta: dict, dest: Path) -> None:
    sub = stats_df[stats_df["contrast"] == contrast].copy()
    sub = sub[np.isfinite(sub["spearman_rho"])]
    if sub.empty:
        return
    labels = list(sub["cohort"]) + [f"RE meta (k={meta.get('k', 0)})"]
    rhos = list(sub["spearman_rho"]) + [meta.get("rho", np.nan)]
    ns = list(sub["n"]) + [meta.get("n_total", 0)]
    fig, ax = plt.subplots(figsize=(8.2, 0.55 * len(labels) + 1.6))
    y = np.arange(len(labels))[::-1]
    for i, (lab, rho, n) in enumerate(zip(labels, rhos, ns)):
        color = "#1f4e79" if i < len(sub) else "#8b1e1e"
        if np.isfinite(rho):
            # Fisher-z 95% CI for cohort rows; meta CI if present
            if i < len(sub) and n >= 4:
                se = 1.0 / math.sqrt(n - 3.0)
                lo, hi = math.tanh(fisher_z(rho) - 1.96 * se), math.tanh(fisher_z(rho) + 1.96 * se)
            else:
                lo, hi = meta.get("ci_lo", rho), meta.get("ci_hi", rho)
            if np.isfinite(lo) and np.isfinite(hi):
                ax.plot([lo, hi], [y[i], y[i]], color=color, lw=1.4)
            ax.plot(rho, y[i], "o", color=color, ms=7)
            ax.text(1.02, y[i], f"n={int(n)}  ρ={rho:.3f}", va="center", fontsize=8, transform=ax.get_yaxis_transform())
        else:
            ax.text(0, y[i], "NA", va="center", ha="center", fontsize=8)
    ax.axvline(0, color="0.5", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlim(-1.05, 1.05)
    ax.set_xlabel("Spearman ρ (patient)")
    ax.set_title(contrast)
    fig.tight_layout()
    dest.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(dest, dpi=160)
    fig.savefig(dest.with_suffix(".pdf"))
    plt.close(fig)


def fmt_rho(row: dict) -> str:
    if not np.isfinite(row.get("spearman_rho", np.nan)):
        return f"n={row['n']} NA ({row.get('note') or 'no ρ'})"
    return f"n={row['n']} ρ={row['spearman_rho']:+.3f} p={row['spearman_p']:.3g}"


def writeup(text_path: Path, payload: dict) -> None:
    lines = []
    lines.append("# Public lung scRNA: B-cell / TLS-like vs malignant TACSTD2 / CLDN4")
    lines.append("")
    lines.append("**Slice:** `methods|scripts|results/scrna_tls_meta/` only.")
    lines.append("**Additive.** User A6 / B6 immune-cold is taken as given and is not re-cut.")
    lines.append("**ICI labels:** not the estimand. This slice cannot test ICI benefit or histologic TLS.")
    lines.append("")
    lines.append("### TL;DR (per-feature, not a pooled overclaim)")
    lines.append("")
    t2b = payload["meta"]["TACSTD2 vs B fraction"]
    t2t = payload["meta"]["TACSTD2 vs TLS12"]
    c4b = payload["meta"]["CLDN4 vs B fraction"]
    c4t = payload["meta"]["CLDN4 vs TLS12"]

    def _meta_s(m):
        if not np.isfinite(m.get("rho", np.nan)):
            return "NA"
        return f"k={m['k']} n={m['n_total']} ρ={m['rho']:+.3f} (95% CI {m['ci_lo']:+.3f} to {m['ci_hi']:+.3f}) p={m['p']:.3g} τ²={m['tau2']:.3f}"

    lines.append(f"- TACSTD2 vs B fraction RE-meta: {_meta_s(t2b)}. Heterogeneous; not a general inverse.")
    lines.append(f"- TACSTD2 vs TLS12 RE-meta: {_meta_s(t2t)}. Not a general inverse.")
    lines.append(f"- CLDN4 vs B fraction RE-meta: {_meta_s(c4b)}.")
    lines.append(f"- CLDN4 vs TLS12 RE-meta: {_meta_s(c4t)}.")
    lines.append("- A6/B6 T-cell immune-cold is taken as given and is not tested here.")
    lines.append("")
    lines.append("### What was tested")
    lines.append("")
    lines.append("Patient-level Spearman of malignant (or epithelial / author gate) TACSTD2 and CLDN4 versus B-cell fraction and the 12-chemokine TLS score, then DerSimonian–Laird random-effects meta. Public processed matrices only.")
    lines.append("")
    lines.append("### Coverage (honest missingness)")
    lines.append("")
    lines.append("| Cohort | Status | Patients eligible | Cells (tumor) | TACSTD2 | CLDN4 | TLS12 k/12 | B cells | Note |")
    lines.append("|---|---|---:|---:|---|---|---:|---:|---|")
    for c in payload["coverage"]:
        lines.append(
            f"| {c['cohort']} | {c['status']} | {c.get('n_eligible','')} | {c.get('n_tumor_cells','')} | "
            f"{c.get('tacstd2','')} | {c.get('cldn4','')} | {c.get('tls12_k','')} | {c.get('n_B','')} | {c.get('note','')} |"
        )
    lines.append("")
    lines.append("### Primary patient-level Spearman (honest n / ρ / p)")
    lines.append("")
    for x, y, title in PRIMARY:
        key = f"{x}__{y}"
        lines.append(f"**{title}**")
        lines.append("")
        lines.append("| Cohort | n | ρ | p | note |")
        lines.append("|---|---:|---:|---:|---|")
        for r in payload["stats"]:
            if r["contrast"] != title:
                continue
            rho = r["spearman_rho"]
            p = r["spearman_p"]
            rho_s = "NA" if not np.isfinite(rho) else f"{rho:+.3f}"
            p_s = "NA" if not np.isfinite(p) else f"{p:.3g}"
            lines.append(f"| {r['cohort']} | {r['n']} | {rho_s} | {p_s} | {r.get('note','')} |")
        m = payload["meta"][title]
        rho_s = "NA" if not np.isfinite(m.get("rho", np.nan)) else f"{m['rho']:+.3f}"
        p_s = "NA" if not np.isfinite(m.get("p", np.nan)) else f"{m['p']:.3g}"
        ci = ""
        if np.isfinite(m.get("ci_lo", np.nan)):
            ci = f" (95% CI {m['ci_lo']:+.3f} to {m['ci_hi']:+.3f})"
        tau = m.get("tau2")
        tau_s = "" if not np.isfinite(tau) else f"τ²={tau:.3f}"
        lines.append(f"| RE meta k={m.get('k')} n={m.get('n_total')} | {m.get('n_total')} | {rho_s}{ci} | {p_s} | {tau_s} {m.get('note','')} |")
        lines.append("")
    lines.append("### Sensitivity: partial Spearman residualizing epithelial fraction")
    lines.append("")
    lines.append("Composition check. `frac_B` can fall when more epithelial cells are captured. Ranks of x and y are residualized on ranks of epithelial/gate fraction.")
    lines.append("")
    lines.append("| Contrast | Cohort | n | ρ | p |")
    lines.append("|---|---|---:|---:|---:|")
    for r in payload["stats"]:
        if "given epi fraction" not in r["contrast"]:
            continue
        rho = r["spearman_rho"]
        p = r["spearman_p"]
        rho_s = "NA" if not np.isfinite(rho) else f"{rho:+.3f}"
        p_s = "NA" if not np.isfinite(p) else f"{p:.3g}"
        lines.append(f"| {r['contrast']} | {r['cohort']} | {r['n']} | {rho_s} | {p_s} |")
    lines.append("")
    lines.append("### Extra: TLS at single-cell")
    lines.append("")
    lines.append("CXCL13+ / MS4A1+ detection by lineage. This is not a follicle or spatial TLS call.")
    lines.append("")
    lines.append("See `tables/tls_single_cell.tsv`. CXCL13+ is T-enriched in mixed TME objects; MS4A1+ is B-restricted. GSE154826 TLS12 is scored inside a CD45-bead library (immune-internal), not a dissociated whole-tumor fraction.")
    lines.append("")
    lines.append("| Cohort | lineage | n_cells | CXCL13 %pos | MS4A1 %pos |")
    lines.append("|---|---|---:|---:|---:|")
    for r in payload.get("tls_sc", []):
        if r.get("lineage") not in {"T", "B", "ALL", "gate"}:
            continue
        c13 = r.get("CXCL13_pct_pos")
        ms = r.get("MS4A1_pct_pos")
        c13s = "NA" if c13 is None or (isinstance(c13, float) and (c13 != c13)) else f"{float(c13):.2f}"
        mss = "NA" if ms is None or (isinstance(ms, float) and (ms != ms)) else f"{float(ms):.2f}"
        lines.append(f"| {r['cohort']} | {r['lineage']} | {r['n_cells']} | {c13s} | {mss} |")
    lines.append("")
    lines.append("### Caveats")
    lines.append("")
    lines.append("1. Dissociated 10x / Rhapsody is not histologic TLS.")
    lines.append("2. GSE154826 is CD45-bead CITE-seq; epithelium sits in the author gate. Sort bias is recorded.")
    lines.append("3. GSE253013 and GSE148071 malignant calls are marker proxies, not CopyKAT.")
    lines.append("4. GSE241934 IIT and RWC are non-overlapping patients and are two strata of one paper.")
    lines.append("5. Underpowered cohorts are labeled inconclusive. A non-significant ρ is not evidence of no association.")
    lines.append("6. A6/B6 T-cell immune-cold is not re-tested.")
    lines.append("")
    lines.append("### Reproduce")
    lines.append("")
    lines.append("```bash")
    lines.append("python3 methods/scrna_tls_meta/download.py")
    lines.append("python3 methods/scrna_tls_meta/extract.py")
    lines.append("python3 methods/scrna_tls_meta/analyze.py")
    lines.append("```")
    lines.append("")
    text_path.write_text("\n".join(lines) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--extracted", type=Path, default=EXTRACTED)
    ap.add_argument("--outdir", type=Path, default=OUT)
    args = ap.parse_args()
    out = args.outdir
    (out / "tables").mkdir(parents=True, exist_ok=True)
    (out / "figures").mkdir(parents=True, exist_ok=True)

    parquets = sorted(args.extracted.glob("*.parquet"))
    coverage = []
    patients = []
    sc_rows = []
    skipped = []

    expected = [
        "GSE131907",
        "GSE253013",
        "GSE207422",
        "GSE241934_IIT",
        "GSE241934_RWC",
        "GSE154826",
        "GSE148071",
    ]
    have = {p.stem for p in parquets}
    for name in expected:
        if name not in have:
            coverage.append(
                {
                    "cohort": name,
                    "status": "skipped",
                    "note": "no extracted parquet (download/extract failed or not run)",
                    "n_eligible": 0,
                    "n_tumor_cells": 0,
                    "tacstd2": "NA",
                    "cldn4": "NA",
                    "tls12_k": 0,
                    "n_B": 0,
                }
            )
            skipped.append(name)

    for p in parquets:
        df = pd.read_parquet(p)
        cohort = str(df["cohort"].iloc[0]) if "cohort" in df.columns else p.stem
        tumor = df[df["tissue"].astype(str).str.lower().eq("tumor")]
        if tumor.empty:
            tumor = df
        pt = patient_table(df, cohort)
        patients.append(pt)
        sc_rows.append(tls_single_cell_table(df, cohort))
        elig = pt[pt["eligible"]]
        coverage.append(
            {
                "cohort": cohort,
                "status": "included" if len(elig) >= MIN_N_SPEARMAN else "underpowered",
                "n_eligible": int(len(elig)),
                "n_patients_total": int(len(pt)),
                "n_tumor_cells": int(len(tumor)),
                "tacstd2": "yes" if "TACSTD2" in df.columns else "missing",
                "cldn4": "yes" if "CLDN4" in df.columns else "missing",
                "tls12_k": int(sum(1 for g in TLS12 if g in df.columns)),
                "n_B": int((tumor["lineage"].astype(str) == "B").sum()),
                "note": "; ".join(sorted(elig["compartment"].unique())) if len(elig) else "n<5 after floors",
            }
        )

    if patients:
        pat = pd.concat(patients, ignore_index=True)
    else:
        pat = pd.DataFrame()
    sc = pd.concat(sc_rows, ignore_index=True) if sc_rows else pd.DataFrame()

    stats_rows = []
    if not pat.empty:
        for cohort, g in pat.groupby("cohort"):
            use = g[g["eligible"]]
            for x, y, title in PRIMARY:
                stats_rows.append(spearman(use[x], use[y], title, cohort))
            # extra honest rows
            stats_rows.append(spearman(use["tacstd2_mal_mean_log1p_cp10k"], use["cxcl13_mean_log1p_cp10k"], "TACSTD2 vs CXCL13 mean", cohort))
            stats_rows.append(spearman(use["tacstd2_mal_mean_log1p_cp10k"], use["ms4a1_mean_log1p_cp10k"], "TACSTD2 vs MS4A1 mean", cohort))
            stats_rows.append(spearman(use["cldn4_mal_mean_log1p_cp10k"], use["cxcl13_mean_log1p_cp10k"], "CLDN4 vs CXCL13 mean", cohort))
            stats_rows.append(spearman(use["cldn4_mal_mean_log1p_cp10k"], use["ms4a1_mean_log1p_cp10k"], "CLDN4 vs MS4A1 mean", cohort))
            stats_rows.append(spearman(use["tacstd2_mal_mean_log1p_cp10k"], use["frac_B_plasma"], "TACSTD2 vs B+plasma fraction", cohort))
            stats_rows.append(
                partial_spearman(
                    use["tacstd2_mal_mean_log1p_cp10k"],
                    use["frac_B"],
                    use["frac_epithelial"],
                    "TACSTD2 vs B fraction given epi fraction",
                    cohort,
                )
            )
            stats_rows.append(
                partial_spearman(
                    use["cldn4_mal_mean_log1p_cp10k"],
                    use["tls12_z"],
                    use["frac_epithelial"],
                    "CLDN4 vs TLS12 given epi fraction",
                    cohort,
                )
            )

    stats_df = pd.DataFrame(stats_rows)
    # FDR on primary only
    prim_titles = {t for _, _, t in PRIMARY}
    if not stats_df.empty:
        mask = stats_df["contrast"].isin(prim_titles)
        q = np.full(len(stats_df), np.nan)
        q[mask.to_numpy()] = bh_fdr(stats_df.loc[mask, "spearman_p"].tolist())
        stats_df["bh_q_primary"] = q

    meta = {}
    for x, y, title in PRIMARY:
        sub = [r for r in stats_rows if r["contrast"] == title]
        meta[title] = random_effects_meta(sub)
        forest_plot(stats_df, title, meta[title], out / "figures" / f"forest_{x}_{y}.png")

    # scatter grid for primary TACSTD2 vs frac_B
    if not pat.empty:
        elig = pat[pat["eligible"]]
        cohorts = list(elig["cohort"].unique())
        ncol = min(3, max(1, len(cohorts)))
        nrow = int(math.ceil(len(cohorts) / ncol)) if cohorts else 1
        fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.4 * nrow), squeeze=False)
        for i, cohort in enumerate(cohorts):
            ax = axes[i // ncol][i % ncol]
            g = elig[elig["cohort"] == cohort]
            ax.scatter(g["tacstd2_mal_mean_log1p_cp10k"], g["frac_B"], s=22, c="#1f4e79", alpha=0.85)
            ax.set_title(f"{cohort} n={len(g)}", fontsize=9)
            ax.set_xlabel("malignant TACSTD2 mean log1p(CP10K)")
            ax.set_ylabel("B fraction")
        for j in range(len(cohorts), nrow * ncol):
            axes[j // ncol][j % ncol].axis("off")
        fig.tight_layout()
        fig.savefig(out / "figures" / "scatter_tacstd2_vs_fracB.png", dpi=160)
        plt.close(fig)

    pat.to_csv(out / "tables" / "patient_scores.tsv", sep="\t", index=False)
    stats_df.to_csv(out / "tables" / "spearman_patient.tsv", sep="\t", index=False)
    sc.to_csv(out / "tables" / "tls_single_cell.tsv", sep="\t", index=False)
    pd.DataFrame(coverage).to_csv(out / "tables" / "coverage.tsv", sep="\t", index=False)
    meta_df = []
    for title, m in meta.items():
        rec = {"contrast": title}
        rec.update(m)
        rec["cohorts"] = ",".join(m.get("cohorts", [])) if isinstance(m.get("cohorts"), list) else m.get("cohorts", "")
        meta_df.append(rec)
    pd.DataFrame(meta_df).to_csv(out / "tables" / "meta_random_effects.tsv", sep="\t", index=False)

    payload = {
        "slice": "scrna_tls_meta",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "coverage": coverage,
        "stats": stats_rows,
        "meta": meta,
        "n_patients_eligible": int(pat["eligible"].sum()) if not pat.empty else 0,
    }
    # json-safe
    def _clean(o):
        if isinstance(o, dict):
            return {k: _clean(v) for k, v in o.items()}
        if isinstance(o, list):
            return [_clean(v) for v in o]
        if isinstance(o, float) and (math.isnan(o) or math.isinf(o)):
            return None
        if isinstance(o, (np.floating,)):
            v = float(o)
            return None if (math.isnan(v) or math.isinf(v)) else v
        if isinstance(o, (np.integer,)):
            return int(o)
        return o

    (out / "summary.json").write_text(json.dumps(_clean(payload), indent=2) + "\n")
    payload["tls_sc"] = sc.to_dict(orient="records") if not sc.empty else []
    writeup(out / "WRITEUP.md", payload)
    # also drop a pointer in methods/
    (ROOT / "methods" / "scrna_tls_meta" / "RESULTS_POINTER.md").write_text(
        "Results live in `results/scrna_tls_meta/` (WRITEUP.md, tables/, figures/).\n"
        "This methods folder is playbook + runnable scripts only.\n"
    )
    print("wrote", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
