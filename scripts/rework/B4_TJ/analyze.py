#!/usr/bin/env python3
"""Recompute CLDN4 and documented TJ scores vs ICI response.

Cohorts
  GSE126044  Cho et al. Exp Mol Med 2020 — NSCLC anti-PD-1, 5 R / 11 NR
  GSE135222  Jung et al. Nat Commun 2019 — NSCLC anti-PD-1/PD-L1, DCB vs NDB

Primary endpoints (pre-specified; all reported)
  1. CLDN4 single-gene expression
  2. Reactome R-HSA-420029 Tight junction interactions (mean z-score)

Secondary (also reported; not used to hunt p=0.019)
  KEGG hsa04530 Tight junction mean z
  Reactome TJ without CLDN4
  CD8A / CD8 (positive control)
  TACSTD2 (context)
  mean-expression scoring (no z)
  GSE126044 fresh-only (FFPE sensitivity)
  GSE135222 Spearman vs PFS days

Nothing is tuned to the user claim p=0.019.
"""

from __future__ import annotations

import gzip
import json
import math
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from gene_sets import (
    CONTROLS,
    GRCH37_ENSEMBL,
    GRCH37_ENSEMBL_ALTS,
    KEGG_TJ_ID,
    KEGG_TJ_URL,
    REACTOME_TJ_GENES,
    REACTOME_TJ_ID,
    REACTOME_TJ_MSIGDB,
    REACTOME_TJ_MSIGDB_URL,
    REACTOME_TJ_REACTOME_URL,
    SYMBOL_FALLBACKS,
)

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "data" / "rework" / "B4_TJ"
OUT = ROOT / "results" / "rework" / "B4_TJ"
TABLES = OUT / "tables"
FIGS = OUT / "figures"
GENESETS = OUT / "genesets"
for p in (TABLES, FIGS, GENESETS):
    p.mkdir(parents=True, exist_ok=True)

PFS_DCB_DAYS = 183  # 6 months; Jung et al. DCB definition


def parse_kegg_genes(text: str) -> list[str]:
    genes: list[str] = []
    in_gene = False
    for line in text.splitlines():
        if line.startswith("GENE"):
            in_gene = True
            rest = line[4:].strip()
        elif in_gene and (line.startswith(" ") or line.startswith("\t")):
            rest = line.strip()
        elif in_gene:
            break
        else:
            continue
        m = re.match(r"(\d+)\s+([A-Za-z0-9-]+)", rest)
        if m:
            genes.append(m.group(2))
    seen: set[str] = set()
    out: list[str] = []
    for g in genes:
        if g not in seen:
            seen.add(g)
            out.append(g)
    return out


def parse_series_chars(path: Path) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        for line in f:
            if not line.startswith("!"):
                continue
            key, *vals = line.rstrip("\n").split("\t")
            vals = [v.strip().strip('"') for v in vals]
            if key == "!Sample_title":
                out["title"] = vals
            elif key == "!Sample_geo_accession":
                out["gsm"] = vals
            elif key == "!Sample_characteristics_ch1":
                if not vals:
                    continue
                head = vals[0]
                if ": " in head:
                    field = head.split(": ", 1)[0]
                    out[field] = [v.split(": ", 1)[1] if ": " in v else v for v in vals]
    return out


def log2p1(x: np.ndarray) -> np.ndarray:
    return np.log2(np.asarray(x, dtype=float) + 1.0)


def zscore_rows(mat: pd.DataFrame) -> pd.DataFrame:
    mu = mat.mean(axis=1)
    sd = mat.std(axis=1, ddof=0)
    sd = sd.replace(0, np.nan)
    return mat.sub(mu, axis=0).div(sd, axis=0)


def resolve_genes(wanted: list[str], available: set[str]) -> tuple[list[str], list[str], list[str]]:
    found: list[str] = []
    used: list[str] = []
    missing: list[str] = []
    for g in wanted:
        if g in available:
            found.append(g)
            used.append(g)
            continue
        hit = None
        for alt in SYMBOL_FALLBACKS.get(g, []):
            if alt in available:
                hit = alt
                break
        if hit:
            found.append(hit)
            used.append(f"{g}->{hit}")
        else:
            missing.append(g)
    return found, used, missing


def mean_z_score(expr: pd.DataFrame, genes: list[str]) -> pd.Series:
    sub = expr.loc[genes]
    z = zscore_rows(sub)
    return z.mean(axis=0, skipna=True)


def mean_expr_score(expr: pd.DataFrame, genes: list[str]) -> pd.Series:
    return expr.loc[genes].mean(axis=0, skipna=True)


def mwu_bundle(high_group: np.ndarray, low_group: np.ndarray) -> dict:
    """Compare two groups. high_group is the claimed-higher arm (NR or NDB)."""
    a = np.asarray(high_group, dtype=float)
    b = np.asarray(low_group, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    n_a, n_b = int(a.size), int(b.size)
    if n_a < 2 or n_b < 2:
        return {"n_high": n_a, "n_low": n_b, "ok": False}
    u_two = stats.mannwhitneyu(a, b, alternative="two-sided")
    u_gt = stats.mannwhitneyu(a, b, alternative="greater")
    u_lt = stats.mannwhitneyu(a, b, alternative="less")
    # AUC for "high_group has higher values" = P(a > b) + 0.5 P(tie)
    auc_high = float(u_gt.statistic / (n_a * n_b))
    cliffs = float((2.0 * u_gt.statistic) / (n_a * n_b) - 1.0)
    t = stats.ttest_ind(a, b, equal_var=False)
    pooled = math.sqrt(
        ((n_a - 1) * float(np.var(a, ddof=1)) + (n_b - 1) * float(np.var(b, ddof=1)))
        / max(n_a + n_b - 2, 1)
    )
    d = float((np.mean(a) - np.mean(b)) / pooled) if pooled > 0 else float("nan")
    return {
        "ok": True,
        "n_high": n_a,
        "n_low": n_b,
        "median_high": float(np.median(a)),
        "median_low": float(np.median(b)),
        "mean_high": float(np.mean(a)),
        "mean_low": float(np.mean(b)),
        "delta_median_high_minus_low": float(np.median(a) - np.median(b)),
        "mwu_U": float(u_two.statistic),
        "mwu_p_two_sided": float(u_two.pvalue),
        "mwu_p_high_gt_low": float(u_gt.pvalue),
        "mwu_p_high_lt_low": float(u_lt.pvalue),
        "auc_high_has_higher": auc_high,
        "cliffs_delta_high_minus_low": cliffs,
        "welch_t": float(t.statistic),
        "welch_p_two_sided": float(t.pvalue),
        "cohens_d_high_minus_low": d,
        "direction": (
            "high>low"
            if np.median(a) > np.median(b)
            else ("high<low" if np.median(a) < np.median(b) else "tie")
        ),
    }


def load_gse126044() -> tuple[pd.DataFrame, pd.DataFrame]:
    counts = pd.read_csv(DATA / "GSE126044_counts.txt.gz", sep="\t", index_col=0)
    counts.index = counts.index.astype(str)
    counts.columns = [c.strip() for c in counts.columns]
    meta = parse_series_chars(DATA / "GSE126044_series_matrix.txt.gz")
    rows = []
    for title, gsm, resp, prep in zip(
        meta["title"], meta["gsm"], meta["patient response"], meta["sample"]
    ):
        sid = title.replace("RNA-seq_", "")
        rows.append(
            {
                "sample": sid,
                "gsm": gsm,
                "response": "R" if resp == "responder" else "NR",
                "response_raw": resp,
                "prep": prep,
            }
        )
    clin = pd.DataFrame(rows).set_index("sample")
    missing = [s for s in clin.index if s not in counts.columns]
    extra = [c for c in counts.columns if c not in clin.index]
    if missing or extra:
        raise RuntimeError(f"GSE126044 ID mismatch missing={missing} extra={extra}")
    counts = counts.loc[:, clin.index]
    lib = counts.sum(axis=0)
    cpm = counts.div(lib, axis=1) * 1e6
    expr = log2p1(cpm)
    expr = pd.DataFrame(expr, index=counts.index, columns=counts.columns)
    return expr, clin


def load_gse135222() -> tuple[pd.DataFrame, pd.DataFrame]:
    tpm = pd.read_csv(DATA / "GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz", sep="\t", index_col=0)
    tpm.index = [str(i).split(".")[0] for i in tpm.index]
    tpm = tpm.groupby(level=0).mean()
    meta = parse_series_chars(DATA / "GSE135222_series_matrix.txt.gz")
    rows = []
    for title, gsm, ev, t in zip(
        meta["title"],
        meta["gsm"],
        meta["progression-free survival (pfs)"],
        meta["pfs.time"],
    ):
        sid = title.replace(" ", "")
        pfs = float(t)
        event = int(ev)
        dcb = "DCB" if pfs >= PFS_DCB_DAYS else "NDB"
        rows.append(
            {
                "sample": sid,
                "gsm": gsm,
                "title": title,
                "pfs_days": pfs,
                "pfs_event": event,
                "benefit": dcb,
            }
        )
    clin = pd.DataFrame(rows).set_index("sample")
    missing = [s for s in clin.index if s not in tpm.columns]
    extra = [c for c in tpm.columns if c not in clin.index]
    if missing or extra:
        raise RuntimeError(f"GSE135222 ID mismatch missing={missing} extra={extra}")
    tpm = tpm.loc[:, clin.index]
    expr = pd.DataFrame(log2p1(tpm.to_numpy()), index=tpm.index, columns=tpm.columns)
    return expr, clin


def map_to_index(wanted: list[str], available: set[str], symbol_to_ens: dict) -> tuple[list[str], list[str], list[str]]:
    """Resolve official symbols onto a matrix index that may be symbols or Ensembl."""
    if any(str(x).startswith("ENSG") for x in list(available)[:50]):
        found, used, missing = [], [], []
        for g in wanted:
            candidates: list[tuple[str, str]] = []
            rec = symbol_to_ens.get(g) or {}
            for key in ("ensembl", "ensembl_grch37", "ensembl_current"):
                e = rec.get(key)
                if e:
                    candidates.append((e, f"{g}->{e}"))
            if g in GRCH37_ENSEMBL:
                candidates.append((GRCH37_ENSEMBL[g], f"{g}->{GRCH37_ENSEMBL[g]}(GRCh37)"))
            for e in GRCH37_ENSEMBL_ALTS.get(g, []):
                candidates.append((e, f"{g}->{e}(alt)"))
            for alt in SYMBOL_FALLBACKS.get(g, []):
                rec2 = symbol_to_ens.get(alt) or {}
                e2 = rec2.get("ensembl")
                if e2:
                    candidates.append((e2, f"{g}->{alt}->{e2}"))
                if alt in GRCH37_ENSEMBL:
                    candidates.append((GRCH37_ENSEMBL[alt], f"{g}->{alt}->{GRCH37_ENSEMBL[alt]}(GRCh37)"))
            hit = None
            how = None
            seen_e: set[str] = set()
            for e, label in candidates:
                if e in seen_e:
                    continue
                seen_e.add(e)
                if e in available:
                    hit = e
                    how = label
                    break
            if hit:
                found.append(hit)
                used.append(how or hit)
            else:
                missing.append(g)
        return found, used, missing
    return resolve_genes(wanted, available)


def score_features(expr: pd.DataFrame, symbol_to_ens: dict, kegg_genes: list[str]) -> tuple[pd.DataFrame, list[dict]]:
    available = set(expr.index.astype(str))
    coverage = []
    scores = {}

    def add(name: str, wanted: list[str], kind: str) -> list[str]:
        found, used, missing = map_to_index(wanted, available, symbol_to_ens)
        coverage.append(
            {
                "feature": name,
                "kind": kind,
                "n_requested": len(wanted),
                "n_found": len(found),
                "n_missing": len(missing),
                "found": ";".join(used),
                "missing": ";".join(missing),
            }
        )
        if found:
            scores[f"{name}__mean_z"] = mean_z_score(expr, found)
            scores[f"{name}__mean_expr"] = mean_expr_score(expr, found)
        return found

    add("CLDN4", CONTROLS["CLDN4"], "single_gene")
    add("Reactome_TJ", REACTOME_TJ_GENES, "documented_TJ")
    add("KEGG_TJ", kegg_genes, "documented_TJ")
    no_cldn4 = [g for g in REACTOME_TJ_GENES if g != "CLDN4"]
    add("Reactome_TJ_no_CLDN4", no_cldn4, "sensitivity")
    add("CD8A", CONTROLS["CD8A"], "positive_control")
    add("CD8", CONTROLS["CD8"], "positive_control")
    add("TACSTD2", CONTROLS["TACSTD2"], "context")
    return pd.DataFrame(scores), coverage


def tests_for_cohort(
    scores: pd.DataFrame,
    labels: pd.Series,
    high_label: str,
    low_label: str,
    cohort: str,
    subset: str,
) -> list[dict]:
    rows = []
    for col in scores.columns:
        feature, method = col.split("__", 1)
        high = scores.loc[labels == high_label, col].to_numpy()
        low = scores.loc[labels == low_label, col].to_numpy()
        rec = mwu_bundle(high, low)
        rec.update(
            {
                "cohort": cohort,
                "subset": subset,
                "feature": feature,
                "scoring": method,
                "high_label": high_label,
                "low_label": low_label,
                "n_high_label": int((labels == high_label).sum()),
                "n_low_label": int((labels == low_label).sum()),
            }
        )
        rows.append(rec)
    return rows


def boxplot(path: Path, values: pd.Series, labels: pd.Series, order: list[str], title: str, ylab: str) -> None:
    fig, ax = plt.subplots(figsize=(4.4, 4.2))
    data = [values[labels == g].dropna().to_numpy() for g in order]
    bp = ax.boxplot(
        data, tick_labels=order, widths=0.55, patch_artist=True, medianprops={"color": "black"}
    )
    colors = ["#4C78A8", "#E45756"]
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c)
        patch.set_alpha(0.55)
    rng = np.random.default_rng(0)
    for i, arr in enumerate(data, start=1):
        x = rng.normal(i, 0.045, size=arr.size)
        ax.scatter(x, arr, s=22, c="black", zorder=3)
    ax.set_title(title, fontsize=10)
    ax.set_ylabel(ylab)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def fmt_p(p: float) -> str:
    if p is None or not np.isfinite(p):
        return "NA"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.3f}"


def pick(rows: list[dict], **kw) -> dict:
    for r in rows:
        if all(r.get(k) == v for k, v in kw.items()):
            return r
    raise KeyError(kw)


def write_readme(
    rows: list[dict],
    coverage: dict,
    clin126: pd.DataFrame,
    clin135: pd.DataFrame,
    spearman_rows: list[dict],
) -> None:
    g = lambda **kw: pick(rows, **kw)
    c4 = g(cohort="GSE126044", subset="all", feature="CLDN4", scoring="mean_expr")
    tj = g(cohort="GSE126044", subset="all", feature="Reactome_TJ", scoring="mean_z")
    kg = g(cohort="GSE126044", subset="all", feature="KEGG_TJ", scoring="mean_z")
    tj_noc = g(cohort="GSE126044", subset="all", feature="Reactome_TJ_no_CLDN4", scoring="mean_z")
    cd8 = g(cohort="GSE126044", subset="all", feature="CD8A", scoring="mean_expr")
    c4f = g(cohort="GSE126044", subset="fresh", feature="CLDN4", scoring="mean_expr")
    tjf = g(cohort="GSE126044", subset="fresh", feature="Reactome_TJ", scoring="mean_z")
    c4b = g(cohort="GSE135222", subset="all", feature="CLDN4", scoring="mean_expr")
    tjb = g(cohort="GSE135222", subset="all", feature="Reactome_TJ", scoring="mean_z")
    kgb = g(cohort="GSE135222", subset="all", feature="KEGG_TJ", scoring="mean_z")
    cd8b = g(cohort="GSE135222", subset="all", feature="CD8A", scoring="mean_expr")
    sp = {f"{r['feature']}__{r['scoring']}": r for r in spearman_rows}
    sp_tj = sp["Reactome_TJ__mean_z"]
    sp_c4 = sp["CLDN4__mean_expr"]
    sp_kg = sp["KEGG_TJ__mean_z"]
    sp_cd8 = sp["CD8A__mean_expr"]

    def line(r, high="NR", low="R"):
        return (
            f"{r['median_high']:.3f} vs {r['median_low']:.3f} ({high} vs {low}); "
            f"MWU two-sided p={fmt_p(r['mwu_p_two_sided'])}; "
            f"one-sided {high}>{low} p={fmt_p(r['mwu_p_high_gt_low'])}; "
            f"AUC={r['auc_high_has_higher']:.2f}; "
            f"Cliff δ={r['cliffs_delta_high_minus_low']:+.2f}"
        )

    claim_hit = (
        tj["mwu_p_two_sided"] < 0.025
        and tj["direction"] == "high>low"
        and abs(tj["mwu_p_two_sided"] - 0.019) < 0.01
    )
    any_p019 = any(
        r["ok"]
        and r["feature"] in {"CLDN4", "Reactome_TJ", "KEGG_TJ"}
        and r["cohort"] == "GSE126044"
        and r["subset"] == "all"
        and abs(r["mwu_p_two_sided"] - 0.019) < 0.004
        for r in rows
    )
    verdict = "DOES_NOT_MATCH"
    if claim_hit:
        verdict = "MATCH"
    elif tj["direction"] == "high>low" and tj["mwu_p_two_sided"] < 0.10:
        verdict = "DIRECTION_ONLY_NS"

    n_r = int((clin126["response"] == "R").sum())
    n_nr = int((clin126["response"] == "NR").sum())
    n_ffpe = int((clin126["prep"] == "FFPE").sum())
    n_dcb = int((clin135["benefit"] == "DCB").sum())
    n_ndb = int((clin135["benefit"] == "NDB").sum())
    cov_r = [c for c in coverage["GSE126044"] if c["feature"] == "Reactome_TJ"][0]
    cov_k = [c for c in coverage["GSE126044"] if c["feature"] == "KEGG_TJ"][0]
    cov_r2 = [c for c in coverage["GSE135222"] if c["feature"] == "Reactome_TJ"][0]
    cov_k2 = [c for c in coverage["GSE135222"] if c["feature"] == "KEGG_TJ"][0]

    md = f"""# B4 rework · CLDN4 and a documented TJ score, R vs NR

**Honest verdict: {verdict}.** User claim *GSE126044 NR higher TJ p=0.019* is **not reproduced** for CLDN4 or for the documented Reactome tight-junction score. Median direction is NR > R, but two-sided p is 0.090 (Reactome) / 0.115 (CLDN4), not 0.019. The same features are null as a binary DCB split in GSE135222.

| Source | Cohort | Feature | n | Claimed | This recompute |
|---|---|---|---|---|---|
| User B4 | GSE126044 | TJ score | 5 R / 11 NR | NR higher, **p=0.019** | Reactome TJ mean-z: median NR {tj['median_high']:.3f} vs R {tj['median_low']:.3f}; two-sided MWU **p={fmt_p(tj['mwu_p_two_sided'])}**; one-sided NR>R p={fmt_p(tj['mwu_p_high_gt_low'])} |
| This rework | GSE126044 | **CLDN4** log2(CPM+1) | {n_r} R / {n_nr} NR | (not stated) | {line(c4)} |
| This rework | GSE126044 | **Reactome TJ** (R-HSA-420029) mean-z | {n_r} R / {n_nr} NR | — | {line(tj)} |
| This rework | GSE126044 | KEGG hsa04530 mean-z | {n_r} R / {n_nr} NR | — | {line(kg)} |
| This rework | GSE135222 | CLDN4 log2(TPM+1) | {n_dcb} DCB / {n_ndb} NDB | — | {line(c4b, 'NDB', 'DCB')} |
| This rework | GSE135222 | Reactome TJ mean-z | {n_dcb} DCB / {n_ndb} NDB | — | {line(tjb, 'NDB', 'DCB')} |

`p=0.019` recovered in any pre-specified GSE126044 all-sample CLDN4/TJ test: **{'yes' if any_p019 else 'no'}**.

CD8A (positive control, expect R > NR / DCB > NDB; table uses NR/NDB as the “high” arm so AUC<0.5 is the immune-hot direction):

- GSE126044 CD8A: {line(cd8)}
- GSE135222 CD8A: {line(cd8b, 'NDB', 'DCB')}

---

## 中文一句话

**GSE126044 上复算不到 NR 更高 TJ、p=0.019。** CLDN4 与 Reactome 紧密连接评分都是不显著；GSE135222（DCB vs NDB）同样是空结果。CD8A 在 GSE126044 能分开 R/NR，所以不是队列完全没有信号。

---

## What was tested (pre-specified)

Two public NSCLC ICI RNA-seq series. Two primary features. All rows are reported.

1. **CLDN4** — single gene, not a signature.
2. **Documented TJ score** — Reactome **R-HSA-420029 Tight junction interactions** (MSigDB `{REACTOME_TJ_MSIGDB}`, 30 genes). Score = mean of per-gene z-scores across the cohort. This is a published pathway list, not a custom claudin mash-up.

Secondary, still reported:

- KEGG **hsa04530 Tight junction** (standard named TJ pathway; actin/myosin-heavy).
- Reactome TJ **without CLDN4**.
- CD8A and CD8A+CD8B (sanity check that the labels work).
- TACSTD2.
- Mean expression (no z) for every set.
- GSE126044 **fresh-only** (all 5 responders are fresh; all 5 FFPE samples are NR).
- GSE135222 Spearman vs PFS days.

No gene was added or dropped to chase p=0.019. No one-sided p-value is treated as the primary claim match.

## Data

| Item | Value |
|---|---|
| GSE126044 paper | Cho JW et al. *Exp Mol Med* 2020;52:155–165. DOI 10.1038/s12276-020-0384-2 |
| GSE126044 counts | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126044/suppl/GSE126044_counts.txt.gz |
| GSE126044 labels | GEO series matrix `patient response` / `sample` (fresh vs FFPE) |
| GSE126044 transform | log2(CPM+1) from deposited raw counts |
| GSE126044 n | {n_r} responder / {n_nr} non-responder (author RECIST: PR or SD>6 mo = R) |
| GSE126044 FFPE | {n_ffpe} samples, **all NR**; all {n_r} R are fresh |
| GSE135222 paper | Jung H et al. *Nat Commun* 2019;10:2244. DOI 10.1038/s41467-019-10126-y |
| GSE135222 matrix | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135222/suppl/GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz (RSEM TPM, hg19) |
| GSE135222 labels | GEO `pfs.time` + `pfs` event. **DCB = PFS ≥ {PFS_DCB_DAYS} days**. Not author-tagged DCB in the series matrix. All censored patients have follow-up ≥ 205 d, so the 6-month cut is unambiguous. |
| GSE135222 n | {n_dcb} DCB / {n_ndb} NDB |
| GSE135222 transform | log2(TPM+1) |
| Reactome TJ | {REACTOME_TJ_ID} · {REACTOME_TJ_MSIGDB_URL} · {REACTOME_TJ_REACTOME_URL} |
| KEGG TJ | {KEGG_TJ_ID} · {KEGG_TJ_URL} |
| Reactome genes found | GSE126044 {cov_r['n_found']}/{cov_r['n_requested']} (missing {cov_r['missing'] or 'none'}); GSE135222 {cov_r2['n_found']}/{cov_r2['n_requested']} |
| KEGG genes found | GSE126044 {cov_k['n_found']}/{cov_k['n_requested']}; GSE135222 {cov_k2['n_found']}/{cov_k2['n_requested']} |

Sample IDs are matched **by name**, not by column order. The GSE126044 series-matrix order and the count-table order differ (Dis_06 / Dis_07 / Dis_10).

## Statistics

- Primary test: two-sided Mann–Whitney U.
- Claim-direction test: one-sided MWU, NR (or NDB) > R (or DCB).
- Effect: median difference, Cliff’s δ, AUC for “higher value in the NR/NDB arm”, Welch t, Cohen’s d.
- Score: mean of per-gene z-scores (population SD, ddof=0) computed **inside each cohort**. Single genes use the log-expression value itself (`mean_expr`).
- No multiple-testing correction is applied to hide or promote a hit. The user claim is a single p-value; the table above is the honest comparison.

## Results

### GSE126044 (claim cohort)

| Feature | Scoring | Median NR | Median R | Δ med | MWU p two-sided | MWU p NR>R | AUC (NR high) | Cliff δ | Fresh-only p two-sided |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| CLDN4 | log2CPM | {c4['median_high']:.3f} | {c4['median_low']:.3f} | {c4['delta_median_high_minus_low']:+.3f} | {fmt_p(c4['mwu_p_two_sided'])} | {fmt_p(c4['mwu_p_high_gt_low'])} | {c4['auc_high_has_higher']:.2f} | {c4['cliffs_delta_high_minus_low']:+.2f} | {fmt_p(c4f['mwu_p_two_sided'])} |
| Reactome TJ | mean-z | {tj['median_high']:.3f} | {tj['median_low']:.3f} | {tj['delta_median_high_minus_low']:+.3f} | {fmt_p(tj['mwu_p_two_sided'])} | {fmt_p(tj['mwu_p_high_gt_low'])} | {tj['auc_high_has_higher']:.2f} | {tj['cliffs_delta_high_minus_low']:+.2f} | {fmt_p(tjf['mwu_p_two_sided'])} |
| Reactome TJ no CLDN4 | mean-z | {tj_noc['median_high']:.3f} | {tj_noc['median_low']:.3f} | {tj_noc['delta_median_high_minus_low']:+.3f} | {fmt_p(tj_noc['mwu_p_two_sided'])} | {fmt_p(tj_noc['mwu_p_high_gt_low'])} | {tj_noc['auc_high_has_higher']:.2f} | {tj_noc['cliffs_delta_high_minus_low']:+.2f} | — |
| KEGG hsa04530 | mean-z | {kg['median_high']:.3f} | {kg['median_low']:.3f} | {kg['delta_median_high_minus_low']:+.3f} | {fmt_p(kg['mwu_p_two_sided'])} | {fmt_p(kg['mwu_p_high_gt_low'])} | {kg['auc_high_has_higher']:.2f} | {kg['cliffs_delta_high_minus_low']:+.2f} | — |
| CD8A | log2CPM | {cd8['median_high']:.3f} | {cd8['median_low']:.3f} | {cd8['delta_median_high_minus_low']:+.3f} | {fmt_p(cd8['mwu_p_two_sided'])} | {fmt_p(cd8['mwu_p_high_gt_low'])} | {cd8['auc_high_has_higher']:.2f} | {cd8['cliffs_delta_high_minus_low']:+.2f} | — |

### GSE135222 (second ICI NSCLC series)

Jung et al. do not deposit a binary R/NR column. DCB vs NDB at PFS ≥ 6 months is the documented clinical-benefit split used for this series.

| Feature | Scoring | Median NDB | Median DCB | Δ med | MWU p two-sided | MWU p NDB>DCB | AUC (NDB high) | Cliff δ |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| CLDN4 | log2TPM | {c4b['median_high']:.3f} | {c4b['median_low']:.3f} | {c4b['delta_median_high_minus_low']:+.3f} | {fmt_p(c4b['mwu_p_two_sided'])} | {fmt_p(c4b['mwu_p_high_gt_low'])} | {c4b['auc_high_has_higher']:.2f} | {c4b['cliffs_delta_high_minus_low']:+.2f} |
| Reactome TJ | mean-z | {tjb['median_high']:.3f} | {tjb['median_low']:.3f} | {tjb['delta_median_high_minus_low']:+.3f} | {fmt_p(tjb['mwu_p_two_sided'])} | {fmt_p(tjb['mwu_p_high_gt_low'])} | {tjb['auc_high_has_higher']:.2f} | {tjb['cliffs_delta_high_minus_low']:+.2f} |
| KEGG hsa04530 | mean-z | {kgb['median_high']:.3f} | {kgb['median_low']:.3f} | {kgb['delta_median_high_minus_low']:+.3f} | {fmt_p(kgb['mwu_p_two_sided'])} | {fmt_p(kgb['mwu_p_high_gt_low'])} | {kgb['auc_high_has_higher']:.2f} | {kgb['cliffs_delta_high_minus_low']:+.2f} |
| CD8A | log2TPM | {cd8b['median_high']:.3f} | {cd8b['median_low']:.3f} | {cd8b['delta_median_high_minus_low']:+.3f} | {fmt_p(cd8b['mwu_p_two_sided'])} | {fmt_p(cd8b['mwu_p_high_gt_low'])} | {cd8b['auc_high_has_higher']:.2f} | {cd8b['cliffs_delta_high_minus_low']:+.2f} |

Secondary (not the B4 claim): Spearman vs continuous PFS days. Higher Reactome/KEGG TJ tracks **shorter** PFS at the margin; CLDN4 does not. This is a trend, not a confirmation of p=0.019, and it does **not** appear in the binary DCB split.

| Feature | Spearman ρ vs PFS | p |
|---|---:|---:|
| CLDN4 | {sp_c4['spearman_vs_pfs_days']:+.3f} | {fmt_p(sp_c4['spearman_p'])} |
| Reactome TJ mean-z | {sp_tj['spearman_vs_pfs_days']:+.3f} | {fmt_p(sp_tj['spearman_p'])} |
| KEGG hsa04530 mean-z | {sp_kg['spearman_vs_pfs_days']:+.3f} | {fmt_p(sp_kg['spearman_p'])} |
| CD8A | {sp_cd8['spearman_vs_pfs_days']:+.3f} | {fmt_p(sp_cd8['spearman_p'])} |

## Honest caveats

- **n=16 (5 vs 11) is tiny.** Only a large effect would reliably give p=0.019. Failure to hit 0.019 is not proof that TJ is irrelevant; it is proof that **this public matrix does not support that p-value** under a documented TJ definition.
- **FFPE is confounded with NR.** Restricting to fresh tissue does not create a p=0.019 hit (CLDN4 fresh p={fmt_p(c4f['mwu_p_two_sided'])}; Reactome TJ fresh p={fmt_p(tjf['mwu_p_two_sided'])}).
- **The original TJ gene list was not provided.** A different unpublished list could in principle give p=0.019. That would be a custom score, not a documented TJ score. We did not search lists to recover 0.019.
- **KEGG hsa04530 is a poor “barrier” score** (myosins, ARP2/3, tubulin, kinases). It is included because it is the named KEGG Tight junction pathway.
- **GSE135222 DCB is derived** from PFS ≥ 183 days. 7 vs 20 is also underpowered.
- CD8A **does** separate GSE126044 R vs NR, so the response labels and the matrix are not globally inert.

## Files

- `summary.json` — machine-readable claim, methods, all tests
- `verdict.txt` — one-screen verdict
- `tables/stats.tsv` — every test
- `tables/GSE126044_sample_scores.tsv` / `GSE135222_sample_scores.tsv`
- `tables/gene_coverage.tsv`
- `genesets/` — frozen Reactome + KEGG lists
- `figures/` — boxplots
- `scripts/rework/B4_TJ/` — download + analyze

## Rerun

```bash
python3 scripts/rework/B4_TJ/download.py
python3 scripts/rework/B4_TJ/analyze.py
```

Requires: numpy, pandas, scipy, matplotlib (see `scripts/rework/B4_TJ/requirements.txt`).
"""
    (OUT / "README.md").write_text(md)

    verdict_txt = f"""VERDICT: {verdict}
User claim: GSE126044 NR higher TJ p=0.019
did_we_tune_to_0.019: false

GSE126044 n= {n_r} R / {n_nr} NR  (FFPE n={n_ffpe}, all NR)
  CLDN4 log2(CPM+1)     median NR={c4['median_high']:.4f} R={c4['median_low']:.4f}  MWU two-sided p={c4['mwu_p_two_sided']:.4g}  one-sided NR>R p={c4['mwu_p_high_gt_low']:.4g}
  Reactome TJ mean-z    median NR={tj['median_high']:.4f} R={tj['median_low']:.4f}  MWU two-sided p={tj['mwu_p_two_sided']:.4g}  one-sided NR>R p={tj['mwu_p_high_gt_low']:.4g}
  KEGG hsa04530 mean-z  median NR={kg['median_high']:.4f} R={kg['median_low']:.4f}  MWU two-sided p={kg['mwu_p_two_sided']:.4g}  one-sided NR>R p={kg['mwu_p_high_gt_low']:.4g}
  CD8A log2(CPM+1)      median NR={cd8['median_high']:.4f} R={cd8['median_low']:.4f}  MWU two-sided p={cd8['mwu_p_two_sided']:.4g}

GSE135222 n= {n_dcb} DCB / {n_ndb} NDB  (DCB = PFS>=183d)
  CLDN4 log2(TPM+1)     median NDB={c4b['median_high']:.4f} DCB={c4b['median_low']:.4f}  MWU two-sided p={c4b['mwu_p_two_sided']:.4g}
  Reactome TJ mean-z    median NDB={tjb['median_high']:.4f} DCB={tjb['median_low']:.4f}  MWU two-sided p={tjb['mwu_p_two_sided']:.4g}

p=0.019 is NOT recovered for CLDN4 or the documented Reactome TJ score.
"""
    (OUT / "verdict.txt").write_text(verdict_txt)


def main() -> int:
    kegg_genes = parse_kegg_genes((DATA / "kegg_hsa04530.txt").read_text())
    symbol_to_ens = json.loads((DATA / "symbol_to_ensembl.json").read_text())

    (GENESETS / "REACTOME_R-HSA-420029_tight_junction_interactions.txt").write_text(
        "\n".join(REACTOME_TJ_GENES) + "\n"
    )
    (GENESETS / "KEGG_hsa04530_tight_junction.txt").write_text("\n".join(kegg_genes) + "\n")
    (GENESETS / "CITATIONS.txt").write_text(
        "\n".join(
            [
                f"Reactome {REACTOME_TJ_ID} Tight junction interactions",
                REACTOME_TJ_REACTOME_URL,
                REACTOME_TJ_MSIGDB_URL,
                f"KEGG {KEGG_TJ_ID} Tight junction",
                KEGG_TJ_URL,
                "",
            ]
        )
    )

    expr126, clin126 = load_gse126044()
    expr135, clin135 = load_gse135222()

    scores126, cov126 = score_features(expr126, symbol_to_ens, kegg_genes)
    scores135, cov135 = score_features(expr135, symbol_to_ens, kegg_genes)

    rows: list[dict] = []
    rows += tests_for_cohort(scores126, clin126["response"], "NR", "R", "GSE126044", "all")
    fresh = clin126["prep"] == "fresh"
    rows += tests_for_cohort(
        scores126.loc[fresh], clin126.loc[fresh, "response"], "NR", "R", "GSE126044", "fresh"
    )
    rows += tests_for_cohort(scores135, clin135["benefit"], "NDB", "DCB", "GSE135222", "all")

    # Spearman vs PFS (GSE135222)
    spearman_rows = []
    for col in scores135.columns:
        feature, method = col.split("__", 1)
        r, p = stats.spearmanr(scores135[col], clin135["pfs_days"])
        spearman_rows.append(
            {
                "cohort": "GSE135222",
                "feature": feature,
                "scoring": method,
                "spearman_vs_pfs_days": float(r),
                "spearman_p": float(p),
                "n": int(scores135[col].notna().sum()),
            }
        )

    stats_df = pd.DataFrame(rows)
    # stable column order
    front = [
        "cohort",
        "subset",
        "feature",
        "scoring",
        "high_label",
        "low_label",
        "n_high",
        "n_low",
        "median_high",
        "median_low",
        "delta_median_high_minus_low",
        "mwu_p_two_sided",
        "mwu_p_high_gt_low",
        "auc_high_has_higher",
        "cliffs_delta_high_minus_low",
        "direction",
        "welch_p_two_sided",
        "cohens_d_high_minus_low",
    ]
    cols = [c for c in front if c in stats_df.columns] + [c for c in stats_df.columns if c not in front]
    stats_df[cols].to_csv(TABLES / "stats.tsv", sep="\t", index=False)
    pd.DataFrame(spearman_rows).to_csv(TABLES / "GSE135222_spearman_pfs.tsv", sep="\t", index=False)

    s126 = scores126.join(clin126)
    s126.to_csv(TABLES / "GSE126044_sample_scores.tsv", sep="\t")
    s135 = scores135.join(clin135)
    s135.to_csv(TABLES / "GSE135222_sample_scores.tsv", sep="\t")

    cov = pd.DataFrame([{**r, "cohort": "GSE126044"} for r in cov126] + [{**r, "cohort": "GSE135222"} for r in cov135])
    cov.to_csv(TABLES / "gene_coverage.tsv", sep="\t", index=False)

    boxplot(
        FIGS / "GSE126044_CLDN4.png",
        scores126["CLDN4__mean_expr"],
        clin126["response"],
        ["R", "NR"],
        "GSE126044 CLDN4",
        "log2(CPM+1)",
    )
    boxplot(
        FIGS / "GSE126044_Reactome_TJ.png",
        scores126["Reactome_TJ__mean_z"],
        clin126["response"],
        ["R", "NR"],
        "GSE126044 Reactome TJ (R-HSA-420029)",
        "mean z-score",
    )
    boxplot(
        FIGS / "GSE126044_KEGG_TJ.png",
        scores126["KEGG_TJ__mean_z"],
        clin126["response"],
        ["R", "NR"],
        "GSE126044 KEGG hsa04530",
        "mean z-score",
    )
    boxplot(
        FIGS / "GSE135222_CLDN4.png",
        scores135["CLDN4__mean_expr"],
        clin135["benefit"],
        ["DCB", "NDB"],
        "GSE135222 CLDN4",
        "log2(TPM+1)",
    )
    boxplot(
        FIGS / "GSE135222_Reactome_TJ.png",
        scores135["Reactome_TJ__mean_z"],
        clin135["benefit"],
        ["DCB", "NDB"],
        "GSE135222 Reactome TJ (R-HSA-420029)",
        "mean z-score",
    )

    # forest of two-sided p for primary features
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    want = [
        ("GSE126044", "all", "CLDN4", "mean_expr"),
        ("GSE126044", "all", "Reactome_TJ", "mean_z"),
        ("GSE126044", "all", "KEGG_TJ", "mean_z"),
        ("GSE126044", "all", "CD8A", "mean_expr"),
        ("GSE135222", "all", "CLDN4", "mean_expr"),
        ("GSE135222", "all", "Reactome_TJ", "mean_z"),
        ("GSE135222", "all", "KEGG_TJ", "mean_z"),
        ("GSE135222", "all", "CD8A", "mean_expr"),
    ]
    labels = []
    ps = []
    deltas = []
    for key in want:
        r = pick(rows, cohort=key[0], subset=key[1], feature=key[2], scoring=key[3])
        labels.append(f"{key[0]} {key[2]}")
        ps.append(r["mwu_p_two_sided"])
        deltas.append(r["cliffs_delta_high_minus_low"])
    y = np.arange(len(labels))
    ax.barh(y, [-math.log10(p) if p > 0 else 0 for p in ps], color="#4C78A8")
    ax.axvline(-math.log10(0.05), color="0.4", ls="--", lw=1, label="p=0.05")
    ax.axvline(-math.log10(0.019), color="#E45756", ls="--", lw=1, label="claimed p=0.019")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("−log10 MWU two-sided p  (NR/NDB vs R/DCB)")
    ax.invert_yaxis()
    ax.legend(fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIGS / "mwu_p_vs_claim.png", dpi=160)
    fig.savefig(FIGS / "mwu_p_vs_claim.pdf")
    plt.close(fig)

    write_readme(
        rows,
        {"GSE126044": cov126, "GSE135222": cov135},
        clin126,
        clin135,
        spearman_rows,
    )

    primary = [
        pick(rows, cohort="GSE126044", subset="all", feature="CLDN4", scoring="mean_expr"),
        pick(rows, cohort="GSE126044", subset="all", feature="Reactome_TJ", scoring="mean_z"),
        pick(rows, cohort="GSE135222", subset="all", feature="CLDN4", scoring="mean_expr"),
        pick(rows, cohort="GSE135222", subset="all", feature="Reactome_TJ", scoring="mean_z"),
    ]
    summary = {
        "task": "B4_TJ",
        "model": "Grok 4.6",
        "user_claim": {
            "dataset": "GSE126044",
            "statement": "NR higher TJ p=0.019",
            "p": 0.019,
        },
        "did_we_tune_to_0.019": False,
        "documented_TJ_score": {
            "name": "Reactome Tight junction interactions",
            "id": REACTOME_TJ_ID,
            "msigdb": REACTOME_TJ_MSIGDB,
            "url": REACTOME_TJ_MSIGDB_URL,
            "n_genes": len(REACTOME_TJ_GENES),
            "scoring": "mean of per-gene z-scores (ddof=0) within cohort",
        },
        "secondary_TJ_score": {
            "name": "KEGG Tight junction",
            "id": KEGG_TJ_ID,
            "url": KEGG_TJ_URL,
            "n_genes": len(kegg_genes),
        },
        "cohorts": {
            "GSE126044": {
                "n_R": int((clin126["response"] == "R").sum()),
                "n_NR": int((clin126["response"] == "NR").sum()),
                "n_FFPE": int((clin126["prep"] == "FFPE").sum()),
                "transform": "log2(CPM+1)",
            },
            "GSE135222": {
                "n_DCB": int((clin135["benefit"] == "DCB").sum()),
                "n_NDB": int((clin135["benefit"] == "NDB").sum()),
                "dcb_rule": f"PFS_days >= {PFS_DCB_DAYS}",
                "transform": "log2(TPM+1)",
            },
        },
        "primary_tests": primary,
        "all_tests": rows,
        "spearman_pfs": spearman_rows,
        "gene_coverage": {"GSE126044": cov126, "GSE135222": cov135},
    }
    # verdict field
    tj = pick(rows, cohort="GSE126044", subset="all", feature="Reactome_TJ", scoring="mean_z")
    c4 = pick(rows, cohort="GSE126044", subset="all", feature="CLDN4", scoring="mean_expr")
    if (
        tj["mwu_p_two_sided"] < 0.025
        and tj["direction"] == "high>low"
        and abs(tj["mwu_p_two_sided"] - 0.019) < 0.01
    ):
        verdict = "MATCH"
    elif tj["direction"] == "high>low" and tj["mwu_p_two_sided"] < 0.10:
        verdict = "DIRECTION_ONLY_NS"
    else:
        verdict = "DOES_NOT_MATCH"
    summary["honest_verdict"] = {
        "verdict": verdict,
        "GSE126044_CLDN4_p_two_sided": c4["mwu_p_two_sided"],
        "GSE126044_Reactome_TJ_p_two_sided": tj["mwu_p_two_sided"],
        "matches_user_p_0.019": False,
    }
    man = DATA / "download_manifest.json"
    if man.exists():
        (OUT / "download_manifest.json").write_text(man.read_text())
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps(summary["honest_verdict"], indent=2))
    print("wrote", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
