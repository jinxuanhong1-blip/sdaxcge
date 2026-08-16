#!/usr/bin/env python3
"""A9: F11R and PARD3 versus TACSTD2-high public lung tumors.

Pre-specified claim genes: F11R (JAM-A) and PARD3.
Primary public cohorts: TCGA-LUAD and TCGA-LUSC primary tumors, analyzed
separately. Pooled NSCLC is a sensitivity only (TACSTD2-high is LUSC-enriched).
Independent LUAD: OncoSG 2020 via the public cBioPortal API.
Stroma-free sensitivity: DepMap 24Q4 lung cell lines, if the extract is present.

Call rule (pre-specified, same spirit as the sibling A8/A9 panel):
  up_in_tacstd2_high iff Spearman ρ > 0 AND tertile log2FC > 0
  AND Welch t-test BH-FDR < 0.05.
Primary FDR is within-cohort across the two claim genes only.
Comparators (CLDN1/4/7, EPCAM, KRT8, PTPRC, CD8A) are context, not the claim.

Does not tune splits, cohorts, or methods to manufacture a positive call.
"""

from __future__ import annotations

import argparse
import json
import math
import time
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from scipy import stats
from statsmodels.stats.multitest import multipletests

CLAIM_GENES = ["F11R", "PARD3"]
TJ_CONTEXT = ["CLDN1", "CLDN4", "CLDN7"]
EPITHELIAL_CONTROLS = ["EPCAM", "KRT8"]
IMMUNE_CONTROLS = ["PTPRC", "CD8A"]
ALL_PANEL = CLAIM_GENES + TJ_CONTEXT + EPITHELIAL_CONTROLS + IMMUNE_CONTROLS
TARGET = "TACSTD2"

LOWQ, HIGHQ = 1 / 3, 2 / 3

CBIO = "https://www.cbioportal.org/api"
ONCOSG_STUDY = "luad_oncosg_2020"
ONCOSG_PROFILE = f"{ONCOSG_STUDY}_rna_seq_v2_mrna_median_all_sample_Zscores"
ONCOSG_SAMPLES = f"{ONCOSG_STUDY}_rna_seq_v2_mrna"

# Magnitude labels are descriptive, not a second call rule.
RHO_WEAK = 0.20
RHO_MODERATE = 0.40


def _http(method: str, url: str, **kw):
    last = None
    for attempt in range(4):
        try:
            r = requests.request(method, url, timeout=120, **kw)
            r.raise_for_status()
            return r
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(2**attempt)
    raise last


def read_xena(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = df.index.astype(str)
    if df.index.duplicated().any():
        df = df.groupby(level=0).mean()
    return df


def primary_tumors(expr: pd.DataFrame) -> pd.DataFrame:
    cols = [c for c in expr.columns if str(c).endswith("-01")]
    return expr[cols]


def load_purity(path: Path) -> pd.DataFrame:
    raw = pd.read_excel(path, header=None)
    hdr = raw.index[raw.iloc[:, 0].astype(str).str.strip() == "Sample ID"][0]
    names = [str(v).strip() for v in raw.iloc[hdr].tolist()]
    df = raw.iloc[hdr + 1 :].copy()
    df.columns = names
    df = df.rename(columns={"Sample ID": "SampleID", "Cancer type": "CancerType"})
    df = df[df["SampleID"].astype(str).str.startswith("TCGA")].copy()
    df["barcode15"] = df["SampleID"].astype(str).str[:15]
    for col in ["ESTIMATE", "ABSOLUTE", "CPE"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.drop_duplicates("barcode15")
    return df.set_index("barcode15")[["CancerType", "ESTIMATE", "ABSOLUTE", "CPE"]]


def spearman(x: np.ndarray, y: np.ndarray) -> tuple[float, float, int]:
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    n = int(len(x))
    if n < 4:
        return float("nan"), float("nan"), n
    r, p = stats.spearmanr(x, y)
    return float(r), float(p), n


def pearson(x: np.ndarray, y: np.ndarray) -> tuple[float, float, int]:
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    n = int(len(x))
    if n < 4:
        return float("nan"), float("nan"), n
    r, p = stats.pearsonr(x, y)
    return float(r), float(p), n


def partial_spearman(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple[float, float, int]:
    """Partial Spearman = Pearson on ranks after residualizing z."""
    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    x, y, z = x[mask], y[mask], z[mask]
    n = int(len(x))
    if n < 5:
        return float("nan"), float("nan"), n
    rx, ry, rz = stats.rankdata(x), stats.rankdata(y), stats.rankdata(z)
    rxy = np.corrcoef(rx, ry)[0, 1]
    rxz = np.corrcoef(rx, rz)[0, 1]
    ryz = np.corrcoef(ry, rz)[0, 1]
    denom = math.sqrt((1 - rxz**2) * (1 - ryz**2))
    if denom == 0 or not np.isfinite(denom):
        return float("nan"), float("nan"), n
    pr = (rxy - rxz * ryz) / denom
    pr = max(min(float(pr), 1.0), -1.0)
    df = n - 3
    if df <= 0 or abs(pr) >= 1:
        return pr, float("nan"), n
    t = pr * math.sqrt(df / (1 - pr**2))
    p = 2 * stats.t.sf(abs(t), df)
    return pr, float(p), n


def magnitude(rho: float | None) -> str:
    if rho is None or not np.isfinite(rho):
        return "NA"
    a = abs(rho)
    if a < RHO_WEAK:
        return "weak"
    if a < RHO_MODERATE:
        return "moderate"
    return "strong"


def tertile_groups(trop2: pd.Series) -> tuple[pd.Index, pd.Index, float, float]:
    lo_thr = float(trop2.quantile(LOWQ))
    hi_thr = float(trop2.quantile(HIGHQ))
    high = trop2[trop2 >= hi_thr].index
    low = trop2[trop2 <= lo_thr].index
    return high, low, lo_thr, hi_thr


def median_groups(trop2: pd.Series) -> tuple[pd.Index, pd.Index]:
    med = float(trop2.median())
    # drop exact median ties so groups stay disjoint
    high = trop2[trop2 > med].index
    low = trop2[trop2 < med].index
    return high, low


def gene_block(
    trop2: pd.Series,
    gene: pd.Series,
    high: pd.Index,
    low: pd.Index,
    purity: pd.DataFrame | None,
) -> dict:
    x = trop2.to_numpy(float)
    y = gene.reindex(trop2.index).to_numpy(float)
    rho, rp, n = spearman(x, y)
    r, pp, _ = pearson(x, y)
    yh = gene.reindex(high).to_numpy(float)
    yl = gene.reindex(low).to_numpy(float)
    yh = yh[np.isfinite(yh)]
    yl = yl[np.isfinite(yl)]
    if len(yh) >= 3 and len(yl) >= 3:
        t, tp = stats.ttest_ind(yh, yl, equal_var=False)
        u, up = stats.mannwhitneyu(yh, yl, alternative="two-sided")
        log2fc = float(np.mean(yh) - np.mean(yl))
    else:
        t = tp = u = up = log2fc = float("nan")
    rec = {
        "n": n,
        "n_high": int(len(yh)),
        "n_low": int(len(yl)),
        "spearman_rho": rho,
        "spearman_p": rp,
        "pearson_r": r,
        "pearson_p": pp,
        "log2FC_high_minus_low": log2fc,
        "welch_t": float(t) if np.isfinite(t) else float("nan"),
        "welch_p": float(tp) if np.isfinite(tp) else float("nan"),
        "mwu_U": float(u) if np.isfinite(u) else float("nan"),
        "mwu_p": float(up) if np.isfinite(up) else float("nan"),
        "mean_high": float(np.mean(yh)) if len(yh) else float("nan"),
        "mean_low": float(np.mean(yl)) if len(yl) else float("nan"),
        "magnitude": magnitude(rho),
    }
    if purity is not None:
        for pm in ["ESTIMATE", "ABSOLUTE", "CPE"]:
            if pm not in purity.columns:
                continue
            z = purity.reindex(trop2.index)[pm].to_numpy(float)
            pr, ppr, pn = partial_spearman(x, y, z)
            rec[f"partial_rho_{pm}"] = pr
            rec[f"partial_p_{pm}"] = ppr
            rec[f"partial_n_{pm}"] = pn
    return rec


def apply_calls(rows: list[dict], p_key: str, q_key: str, claim_only: bool) -> None:
    idxs = [
        i
        for i, r in enumerate(rows)
        if (not claim_only or r["role"] == "claim") and np.isfinite(r.get(p_key, float("nan")))
    ]
    if not idxs:
        for r in rows:
            r[q_key] = float("nan")
        return
    pvals = [rows[i][p_key] for i in idxs]
    qvals = multipletests(pvals, method="fdr_bh")[1]
    for i, q in zip(idxs, qvals):
        rows[i][q_key] = float(q)
    for i, r in enumerate(rows):
        if q_key not in r:
            r[q_key] = float("nan")


def call_up(row: dict, q_key: str) -> bool:
    rho = row.get("spearman_rho")
    fc = row.get("log2FC_high_minus_low")
    q = row.get(q_key)
    return bool(
        np.isfinite(rho)
        and rho > 0
        and np.isfinite(fc)
        and fc > 0
        and np.isfinite(q)
        and q < 0.05
    )


def role_of(gene: str) -> str:
    if gene in CLAIM_GENES:
        return "claim"
    if gene in TJ_CONTEXT:
        return "tj_context"
    if gene in EPITHELIAL_CONTROLS:
        return "epithelial_control"
    if gene in IMMUNE_CONTROLS:
        return "immune_control"
    return "other"


def cohort_frame(expr: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    present = [g for g in genes if g in expr.index]
    return expr.loc[present].T.astype(float)


def run_expr_cohort(
    name: str,
    expr: pd.DataFrame,
    purity: pd.DataFrame | None,
    note: str,
    split: str = "tertile",
) -> tuple[list[dict], dict, pd.DataFrame]:
    if TARGET not in expr.index:
        raise SystemExit(f"{TARGET} missing in {name}")
    trop2 = expr.loc[TARGET].astype(float)
    if split == "tertile":
        high, low, lo_thr, hi_thr = tertile_groups(trop2)
    elif split == "median":
        high, low = median_groups(trop2)
        lo_thr = hi_thr = float(trop2.median())
    else:
        raise ValueError(split)

    rows = []
    for g in ALL_PANEL:
        if g not in expr.index:
            rows.append(
                {
                    "cohort": name,
                    "split": split,
                    "gene": g,
                    "role": role_of(g),
                    "present": False,
                    "note": note,
                }
            )
            continue
        rec = gene_block(trop2, expr.loc[g].astype(float), high, low, purity)
        rec.update(
            {
                "cohort": name,
                "split": split,
                "gene": g,
                "role": role_of(g),
                "present": True,
                "note": note,
            }
        )
        rows.append(rec)

    apply_calls(rows, "welch_p", "welch_fdr_claim", claim_only=True)
    apply_calls(rows, "spearman_p", "spearman_fdr_claim", claim_only=True)
    apply_calls(rows, "welch_p", "welch_fdr_panel", claim_only=False)
    apply_calls(rows, "spearman_p", "spearman_fdr_panel", claim_only=False)
    for r in rows:
        if r.get("present"):
            r["up_in_tacstd2_high"] = call_up(r, "welch_fdr_claim") if r["role"] == "claim" else call_up(r, "welch_fdr_panel")
        else:
            r["up_in_tacstd2_high"] = False

    meta = {
        "cohort": name,
        "split": split,
        "n_samples": int(trop2.notna().sum()),
        "n_high": int(len(high)),
        "n_low": int(len(low)),
        "tacstd2_low_threshold": lo_thr,
        "tacstd2_high_threshold": hi_thr,
        "mean_TACSTD2_high": float(trop2.loc[high].mean()) if len(high) else None,
        "mean_TACSTD2_low": float(trop2.loc[low].mean()) if len(low) else None,
        "note": note,
    }
    sample = pd.DataFrame(
        {
            "sample_id": trop2.index,
            "cohort": name,
            "TACSTD2": trop2.values,
            "tacstd2_group": [
                "high" if s in high else "low" if s in low else "mid" for s in trop2.index
            ],
        }
    )
    for g in ALL_PANEL:
        if g in expr.index:
            sample[g] = expr.loc[g].reindex(trop2.index).to_numpy(float)
    if purity is not None:
        for pm in ["ESTIMATE", "ABSOLUTE", "CPE"]:
            if pm in purity.columns:
                sample[pm] = purity.reindex(trop2.index)[pm].to_numpy(float)
    return rows, meta, sample


def fetch_oncosg() -> tuple[pd.DataFrame, pd.Series]:
    genes = [TARGET] + ALL_PANEL
    g = _http(
        "POST",
        f"{CBIO}/genes/fetch",
        params={"geneIdType": "HUGO_GENE_SYMBOL"},
        json=genes,
    ).json()
    entrez2hugo = {x["entrezGeneId"]: x["hugoGeneSymbol"] for x in g}
    data = _http(
        "POST",
        f"{CBIO}/molecular-profiles/{ONCOSG_PROFILE}/molecular-data/fetch",
        params={"projection": "SUMMARY"},
        json={"entrezGeneIds": list(entrez2hugo.keys()), "sampleListId": ONCOSG_SAMPLES},
    ).json()
    recs = [
        (entrez2hugo[d["entrezGeneId"]], d["sampleId"], float(d["value"]))
        for d in data
        if d.get("value") is not None and d["entrezGeneId"] in entrez2hugo
    ]
    df = pd.DataFrame(recs, columns=["gene", "sample", "value"])
    expr = df.pivot_table(index="gene", columns="sample", values="value")
    clin = _http(
        "GET",
        f"{CBIO}/studies/{ONCOSG_STUDY}/clinical-data",
        params={"clinicalDataType": "SAMPLE", "projection": "SUMMARY"},
    ).json()
    purity = pd.Series(
        {
            d["sampleId"]: float(d["value"])
            for d in clin
            if d.get("clinicalAttributeId") == "PURITY" and d.get("value") not in (None, "", "NA")
        },
        name="PURITY",
    )
    return expr, purity


def load_depmap(data_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame] | None:
    expr_path = data_dir / "depmap24q4_panel_all_models.csv"
    model_path = data_dir / "Model.csv"
    if not expr_path.exists() or not model_path.exists():
        return None
    expr = pd.read_csv(expr_path)
    model = pd.read_csv(model_path, low_memory=False)
    return expr, model


def run_depmap(expr: pd.DataFrame, model: pd.DataFrame) -> tuple[list[dict], dict, pd.DataFrame]:
    keep = [
        c
        for c in [
            "ModelID",
            "CellLineName",
            "StrippedCellLineName",
            "ModelType",
            "OncotreeLineage",
            "OncotreePrimaryDisease",
            "OncotreeSubtype",
        ]
        if c in model.columns
    ]
    df = expr.merge(model[keep], on="ModelID", how="left")
    for g in [TARGET] + ALL_PANEL:
        if g in df.columns:
            df[g] = pd.to_numeric(df[g], errors="coerce")
    lung = df[(df["OncotreeLineage"] == "Lung") & (df["ModelType"].fillna("") == "Cell Line")].copy()
    lung = lung.dropna(subset=[TARGET]).copy()
    # Build a fake genes x samples frame for reuse of run_expr_cohort
    genes_present = [g for g in [TARGET] + ALL_PANEL if g in lung.columns]
    mat = lung.set_index("ModelID")[genes_present].T
    rows, meta, sample = run_expr_cohort(
        "DepMap24Q4_lung_cell_lines",
        mat,
        purity=None,
        note="SENSITIVITY: DepMap Public 24Q4 lung cell lines only (no stroma). log2(TPM+1).",
        split="tertile",
    )
    extra = lung[
        [c for c in ["ModelID", "CellLineName", "OncotreePrimaryDisease", "OncotreeSubtype"] if c in lung.columns]
    ].copy()
    sample = sample.merge(extra, left_on="sample_id", right_on="ModelID", how="left")
    return rows, meta, sample


def _finite(v) -> bool:
    return v is not None and isinstance(v, (int, float)) and np.isfinite(v)


def honest_for_gene(gene: str, rows: list[dict], primary_cohorts: list[str]) -> dict:
    prim = [r for r in rows if r["gene"] == gene and r["cohort"] in primary_cohorts and r.get("split", "tertile") == "tertile"]
    calls = [bool(r.get("up_in_tacstd2_high")) for r in prim if r.get("present")]
    rhos = [r["spearman_rho"] for r in prim if r.get("present") and _finite(r.get("spearman_rho"))]
    if not prim or not calls:
        label = "NOT_COMPUTED"
        statement = f"{gene}: not computed in a primary public lung cohort."
    elif all(calls):
        label = "SUPPORTED_BOTH_HISTOLOGIES"
        statement = (
            f"{gene} meets the pre-specified up_in_tacstd2_high call in both TCGA-LUAD and "
            f"TCGA-LUSC (tertile split). Spearman ρ = "
            + ", ".join(f"{r['cohort']} {r['spearman_rho']:.3f}" for r in prim if r.get("present"))
            + ". This is a bulk RNA association, not proof of a protein complex or regulation."
        )
    elif any(calls):
        label = "PARTIAL_ONE_HISTOLOGY"
        hit = [r["cohort"] for r in prim if r.get("up_in_tacstd2_high")]
        miss = [r["cohort"] for r in prim if r.get("present") and not r.get("up_in_tacstd2_high")]
        statement = (
            f"{gene} meets the call in {', '.join(hit) or 'none'} and fails it in "
            f"{', '.join(miss) or 'none'}. Do not treat this as a histology-independent partner."
        )
    else:
        label = "NOT_SUPPORTED"
        statement = (
            f"{gene} does not meet the pre-specified up_in_tacstd2_high call in either "
            f"TCGA-LUAD or TCGA-LUSC. "
            + (
                f"ρ values are {', '.join(f'{r:.3f}' for r in rhos)} — "
                "statistically detectable weak correlation is not the same as a robust partner."
                if rhos
                else ""
            )
        )
    # magnitude honesty
    if rhos and all(abs(r) < RHO_WEAK for r in rhos) and label == "SUPPORTED_BOTH_HISTOLOGIES":
        label = "SUPPORTED_BUT_WEAK"
        statement += " Effect size is weak (|ρ|<0.20) in both histologies; do not over-interpret."
    return {
        "gene": gene,
        "label": label,
        "primary_calls": {r["cohort"]: bool(r.get("up_in_tacstd2_high")) for r in prim},
        "primary_rho": {r["cohort"]: r.get("spearman_rho") for r in prim},
        "statement": statement,
    }


def write_figures(out: Path, sample_tcga: pd.DataFrame, stats_df: pd.DataFrame) -> list[str]:
    out.mkdir(parents=True, exist_ok=True)
    written = []
    # Scatter per claim gene, LUAD vs LUSC
    fig, axes = plt.subplots(2, 2, figsize=(9.2, 8.4), sharex=False, sharey=False)
    palette = {"LUAD": "#1f77b4", "LUSC": "#d62728"}
    for i, gene in enumerate(CLAIM_GENES):
        for j, hist in enumerate(["LUAD", "LUSC"]):
            ax = axes[i, j]
            sub = sample_tcga[sample_tcga["histology"] == hist]
            ax.scatter(sub[TARGET], sub[gene], s=10, alpha=0.45, c=palette[hist], edgecolors="none")
            row = stats_df[
                (stats_df["cohort"] == f"TCGA_{hist}")
                & (stats_df["gene"] == gene)
                & (stats_df["split"] == "tertile")
            ]
            rho = float(row["spearman_rho"].iloc[0]) if len(row) else float("nan")
            n = int(row["n"].iloc[0]) if len(row) else 0
            ax.set_title(f"TCGA-{hist}  {gene} vs TACSTD2\nρ={rho:.3f}  n={n}")
            ax.set_xlabel("TACSTD2  log2(norm+1)")
            ax.set_ylabel(f"{gene}  log2(norm+1)")
            ax.grid(True, alpha=0.25)
    fig.tight_layout()
    p = out / "fig_scatter_claim_tcga.png"
    fig.savefig(p, dpi=160)
    fig.savefig(out / "fig_scatter_claim_tcga.pdf")
    plt.close(fig)
    written.append(str(p))

    # Forest of Spearman ρ for claim + context in primary cohorts
    show_genes = CLAIM_GENES + TJ_CONTEXT
    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    y = 0
    yticks = []
    ylabels = []
    for gene in show_genes[::-1]:
        for hist, color in [("LUAD", "#1f77b4"), ("LUSC", "#d62728")]:
            row = stats_df[
                (stats_df["cohort"] == f"TCGA_{hist}")
                & (stats_df["gene"] == gene)
                & (stats_df["split"] == "tertile")
            ]
            if row.empty:
                continue
            rho = float(row["spearman_rho"].iloc[0])
            ax.plot([rho], [y], "o", color=color, markersize=7)
            ax.hlines(y, 0, rho, color=color, linewidth=2)
            yticks.append(y)
            ylabels.append(f"{gene}  {hist}")
            y += 1
        y += 0.35
    ax.axvline(0, color="black", linewidth=0.8)
    ax.axvline(RHO_WEAK, color="gray", linewidth=0.6, linestyle="--")
    ax.axvline(RHO_MODERATE, color="gray", linewidth=0.6, linestyle=":")
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels)
    ax.set_xlabel("Spearman ρ vs TACSTD2")
    ax.set_title("TCGA primary tumors, tertile analysis (context CLDNs shown)")
    ax.grid(True, axis="x", alpha=0.25)
    fig.tight_layout()
    p = out / "fig_rho_forest_tcga.png"
    fig.savefig(p, dpi=160)
    plt.close(fig)
    written.append(str(p))

    # High vs low boxes
    fig, axes = plt.subplots(2, 2, figsize=(8.8, 7.6))
    for i, gene in enumerate(CLAIM_GENES):
        for j, hist in enumerate(["LUAD", "LUSC"]):
            ax = axes[i, j]
            sub = sample_tcga[
                (sample_tcga["histology"] == hist) & (sample_tcga["tacstd2_group"].isin(["high", "low"]))
            ]
            data = [sub.loc[sub["tacstd2_group"] == g, gene].dropna().values for g in ("low", "high")]
            bp = ax.boxplot(data, labels=["TACSTD2-low", "TACSTD2-high"], patch_artist=True, widths=0.55)
            for patch, c in zip(bp["boxes"], ["#9ecae1", "#fc9272"]):
                patch.set_facecolor(c)
            ax.set_title(f"TCGA-{hist}  {gene}")
            ax.set_ylabel(f"{gene}  log2(norm+1)")
            ax.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    p = out / "fig_box_high_low_tcga.png"
    fig.savefig(p, dpi=160)
    plt.close(fig)
    written.append(str(p))
    return written


def _json_safe(obj):
    if isinstance(obj, dict):
        return {k: _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_json_safe(v) for v in obj]
    if isinstance(obj, (np.floating, float)):
        v = float(obj)
        return v if np.isfinite(v) else None
    if isinstance(obj, (np.integer, int)):
        return int(obj)
    if isinstance(obj, (np.bool_, bool)):
        return bool(obj)
    return obj


def write_report(out: Path, summary: dict, stats_df: pd.DataFrame) -> None:
    hv = summary["honest_verdict"]
    f11 = hv["F11R"]
    pard = hv["PARD3"]

    def line(gene: str, cohort: str) -> str:
        row = stats_df[(stats_df["cohort"] == cohort) & (stats_df["gene"] == gene) & (stats_df["split"] == "tertile")]
        if row.empty:
            return f"{cohort}: not computed"
        r = row.iloc[0]
        q = r.get("welch_fdr_claim") if gene in CLAIM_GENES else r.get("welch_fdr_panel")
        qv = f"{float(q):.2g}" if _finite(q) else "NA"
        call = "CALL" if bool(r.get("up_in_tacstd2_high")) else "no-call"
        pr = r.get("partial_rho_CPE")
        prs = f", partial-CPE ρ={float(pr):.3f}" if _finite(pr) else ""
        return (
            f"{cohort}: ρ={float(r['spearman_rho']):.3f} (n={int(r['n'])}), "
            f"tertile log2FC={float(r['log2FC_high_minus_low']):+.3f}, "
            f"Welch FDR={qv}, {call}{prs}"
        )

    md = f"""# A9 — F11R and PARD3 vs TACSTD2-high public lung

## Honest result

**F11R:** {f11['label']}. {f11['statement']}

**PARD3:** {pard['label']}. {pard['statement']}

These are bulk RNA associations in public tumors. They do not show that TROP2
binds F11R or PARD3, that either gene is required for a TROP2-high state, or
that the private cohort's effect sizes transfer.

### Primary public lung (TCGA, within-histology tertiles)

- {line('F11R', 'TCGA_LUAD')}
- {line('F11R', 'TCGA_LUSC')}
- {line('PARD3', 'TCGA_LUAD')}
- {line('PARD3', 'TCGA_LUSC')}

### Independent LUAD (OncoSG 2020)

- {line('F11R', 'OncoSG_LUAD')}
- {line('PARD3', 'OncoSG_LUAD')}

### What this is not

- Not a protein / IHC / spatial result.
- Not a malignant-cell-only result (except the optional DepMap cell-line slice).
- Not evidence that F11R or PARD3 belong to a TROP2-high *program* on their
  own: CLDN1/CLDN4 are stronger public correlates; PARD3 is in the weak tail.
- Pooled LUAD+LUSC is **not** the primary test. TACSTD2-high tumors are
  LUSC-enriched, so an unstratified pool can inflate a junction-gene signal.

## Pre-specified methods

| Item | Choice |
|------|--------|
| Claim genes | F11R, PARD3 |
| Primary cohorts | TCGA-LUAD and TCGA-LUSC primary tumors (`-01`), separate |
| TACSTD2-high | Within-cohort upper tertile; lower tertile = low |
| Primary association | Spearman ρ vs continuous TACSTD2 |
| High-vs-low test | Welch t-test on log2 expression |
| Call | ρ>0 AND log2FC>0 AND Welch BH-FDR<0.05, FDR across the 2 claim genes |
| Purity | Partial Spearman vs Aran CPE / ESTIMATE / ABSOLUTE |
| Independent | OncoSG LUAD 2020 (cBioPortal z-scored RNA; Spearman is invariant) |
| Stroma-free | DepMap 24Q4 lung cell lines, if downloaded |
| Comparators | CLDN1/4/7 (other user TJ genes), EPCAM/KRT8, PTPRC/CD8A |

Median split, pooled NSCLC, and DepMap are labeled sensitivity. They were not
used to flip a primary call.

## Sources

- UCSC Xena `TCGA.LUAD.sampleMap/HiSeqV2` and `TCGA.LUSC.sampleMap/HiSeqV2` (log2(norm_count+1))
- Aran, Sirota & Butte, *Nat Commun* 2015, Supplementary Data 1 (purity)
- OncoSG LUAD (Chen et al. 2020) via [cBioPortal](https://www.cbioportal.org/study/summary?id=luad_oncosg_2020)
- DepMap Public 24Q4 (optional): https://doi.org/10.25452/figshare.plus.27993248.v1

## Reproduce

```bash
python3 scripts/w200/A9_F11R_PARD3/download.py
python3 scripts/w200/A9_F11R_PARD3/analyze.py
```
"""
    (out / "README.md").write_text(md)


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--data-dir", default="results/w200/A9_F11R_PARD3/data")
    p.add_argument("--out-dir", default="results/w200/A9_F11R_PARD3")
    p.add_argument("--skip-oncosg", action="store_true")
    args = p.parse_args()
    data = Path(args.data_dir)
    out = Path(args.out_dir)
    figdir = out / "figures"
    out.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(parents=True, exist_ok=True)

    luad = primary_tumors(read_xena(data / "TCGA.LUAD.HiSeqV2.gz"))
    lusc = primary_tumors(read_xena(data / "TCGA.LUSC.HiSeqV2.gz"))
    purity = load_purity(data / "Aran_CPE_purity.xlsx")

    all_rows: list[dict] = []
    metas: list[dict] = []
    samples: list[pd.DataFrame] = []

    for name, mat, note in [
        ("TCGA_LUAD", luad, "PRIMARY: TCGA-LUAD primary solid tumors (barcode -01)."),
        ("TCGA_LUSC", lusc, "PRIMARY: TCGA-LUSC primary solid tumors (barcode -01)."),
    ]:
        rows, meta, sample = run_expr_cohort(name, mat, purity, note, split="tertile")
        all_rows.extend(rows)
        metas.append(meta)
        sample = sample.copy()
        sample["histology"] = name.replace("TCGA_", "")
        samples.append(sample)
        # median-split sensitivity
        rows_m, meta_m, _ = run_expr_cohort(
            name, mat, purity, "SENSITIVITY: median split (not the primary call).", split="median"
        )
        all_rows.extend(rows_m)
        metas.append(meta_m)

    # Pooled NSCLC sensitivity
    genes = luad.index.intersection(lusc.index)
    pooled = pd.concat([luad.loc[genes], lusc.loc[genes]], axis=1)
    rows_p, meta_p, sample_p = run_expr_cohort(
        "TCGA_NSCLC_pooled",
        pooled,
        purity,
        "SENSITIVITY ONLY: pooled LUAD+LUSC. TACSTD2-high is LUSC-enriched; not the primary test.",
        split="tertile",
    )
    all_rows.extend(rows_p)
    metas.append(meta_p)
    sample_p = sample_p.copy()
    sample_p["histology"] = sample_p["sample_id"].map(
        lambda s: "LUAD" if s in luad.columns else "LUSC"
    )
    samples.append(sample_p)

    oncosg_ok = False
    if not args.skip_oncosg:
        try:
            oexpr, opurity = fetch_oncosg()
            op = pd.DataFrame({"PURITY": opurity})
            # map PURITY into the ESTIMATE-style slot so gene_block can use it
            op = op.rename(columns={"PURITY": "CPE"})
            rows_o, meta_o, sample_o = run_expr_cohort(
                "OncoSG_LUAD",
                oexpr,
                op,
                "INDEPENDENT: OncoSG LUAD 2020, cBioPortal z-scored RNA. Spearman is invariant to the z-transform. Partial uses OncoSG PURITY (stored as CPE column).",
                split="tertile",
            )
            all_rows.extend(rows_o)
            metas.append(meta_o)
            sample_o = sample_o.copy()
            sample_o["histology"] = "LUAD"
            samples.append(sample_o)
            oncosg_ok = True
        except Exception as exc:  # noqa: BLE001
            metas.append({"cohort": "OncoSG_LUAD", "error": str(exc)})

    depmap = load_depmap(data)
    depmap_ok = False
    if depmap is not None:
        rows_d, meta_d, sample_d = run_depmap(*depmap)
        all_rows.extend(rows_d)
        metas.append(meta_d)
        samples.append(sample_d)
        depmap_ok = True

    stats_df = pd.DataFrame(all_rows)
    # stable column order
    front = [
        "cohort",
        "split",
        "gene",
        "role",
        "present",
        "n",
        "n_high",
        "n_low",
        "spearman_rho",
        "spearman_p",
        "spearman_fdr_claim",
        "spearman_fdr_panel",
        "pearson_r",
        "pearson_p",
        "log2FC_high_minus_low",
        "welch_t",
        "welch_p",
        "welch_fdr_claim",
        "welch_fdr_panel",
        "mwu_p",
        "up_in_tacstd2_high",
        "magnitude",
        "partial_rho_CPE",
        "partial_p_CPE",
        "partial_n_CPE",
        "partial_rho_ESTIMATE",
        "partial_p_ESTIMATE",
        "partial_n_ESTIMATE",
        "partial_rho_ABSOLUTE",
        "partial_p_ABSOLUTE",
        "partial_n_ABSOLUTE",
        "mean_high",
        "mean_low",
        "note",
    ]
    cols = [c for c in front if c in stats_df.columns] + [c for c in stats_df.columns if c not in front]
    stats_df = stats_df[cols]
    stats_df.to_csv(out / "stats.tsv", sep="\t", index=False)

    sample_all = pd.concat(samples, ignore_index=True)
    sample_all.to_csv(out / "sample_expression.tsv", sep="\t", index=False)
    pd.DataFrame(metas).to_csv(out / "cohort_summary.tsv", sep="\t", index=False)

    primary = ["TCGA_LUAD", "TCGA_LUSC"]
    honest = {
        "F11R": honest_for_gene("F11R", all_rows, primary),
        "PARD3": honest_for_gene("PARD3", all_rows, primary),
    }

    # Purity attenuation note
    purity_note = {}
    for gene in CLAIM_GENES:
        for cohort in primary:
            recs = [
                r
                for r in all_rows
                if r["gene"] == gene and r["cohort"] == cohort and r.get("split") == "tertile" and r.get("present")
            ]
            if not recs:
                continue
            r = recs[0]
            purity_note[f"{gene}:{cohort}"] = {
                "rho": r.get("spearman_rho"),
                "partial_CPE": r.get("partial_rho_CPE"),
                "partial_ESTIMATE": r.get("partial_rho_ESTIMATE"),
            }

    sample_tcga = pd.concat(
        [s for s in samples if s["cohort"].iloc[0] in ("TCGA_LUAD", "TCGA_LUSC")],
        ignore_index=True,
    )
    figures = write_figures(figdir, sample_tcga, stats_df)

    summary = {
        "task": "A9_F11R_PARD3",
        "question": "Are F11R and PARD3 higher / co-expressed in TACSTD2-high public lung tumors?",
        "claim_genes": CLAIM_GENES,
        "primary_cohorts": primary,
        "call_rule": "spearman_rho>0 AND tertile_log2FC>0 AND welch_fdr_claim<0.05",
        "did_we_tune_to_a_positive_call": False,
        "honest_verdict": honest,
        "purity_partial": purity_note,
        "cohorts": metas,
        "oncosg_ok": oncosg_ok,
        "depmap_ok": depmap_ok,
        "figures": figures,
        "methods": {
            "expression": "UCSC Xena HiSeqV2 log2(norm_count+1); primary tumors only",
            "tacstd2_high": "within-cohort tertiles (low<=Q1/3, high>=Q2/3)",
            "purity": "Aran et al. 2015 CPE/ESTIMATE/ABSOLUTE; first-order partial Spearman on ranks",
            "oncosg": f"{ONCOSG_STUDY} {ONCOSG_PROFILE}" if oncosg_ok else "failed_or_skipped",
            "comparators": {
                "tj_context": TJ_CONTEXT,
                "epithelial": EPITHELIAL_CONTROLS,
                "immune": IMMUNE_CONTROLS,
            },
        },
        "generated_utc": datetime.now(timezone.utc).isoformat(),
    }
    (out / "summary.json").write_text(json.dumps(_json_safe(summary), indent=2) + "\n")
    write_report(out, summary, stats_df)

    print(json.dumps(_json_safe(honest), indent=2))
    show = stats_df[
        (stats_df["split"] == "tertile")
        & (stats_df["gene"].isin(CLAIM_GENES))
        & (stats_df["cohort"].isin(primary + ["OncoSG_LUAD", "DepMap24Q4_lung_cell_lines", "TCGA_NSCLC_pooled"]))
    ][
        [
            "cohort",
            "gene",
            "n",
            "spearman_rho",
            "spearman_p",
            "log2FC_high_minus_low",
            "welch_fdr_claim",
            "up_in_tacstd2_high",
            "partial_rho_CPE",
            "magnitude",
        ]
    ]
    print(show.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
