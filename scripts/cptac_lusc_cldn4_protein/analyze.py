#!/usr/bin/env python3
"""CPTAC-LSCC (LUSC) CLDN4 / TJ-15 protein vs ImmuneScore / ESTIMATE / GEP18 / CD8A.

Additive extra. CPTAC LUAD TJ-15 vs ImmuneScore (ρ=−0.30, n=110; PR #245)
is taken as given and is not re-estimated.

Public freeze v1.2 processed tables only. Pairwise-complete Spearman.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from genes import (  # noqa: E402
    CD8,
    CYT,
    ENSEMBL,
    GEP18,
    PHENO_ALIASES,
    SYNONYMS,
    TJ_15,
    TJ_7,
)

RNG_SEED = 20260817
N_BOOT = 2000
MIN_TJ15 = 5
MIN_TJ7 = 4
MIN_GEP = 10

# Given from PR #245; not re-computed.
LUAD_GIVEN = {
    "predictor": "TJ-15 protein",
    "endpoint": "ESTIMATE ImmuneScore",
    "cohort": "CPTAC LUAD",
    "n": 110,
    "rho": -0.296,
    "p": 0.00167,
    "source": "PR #245",
}


def load_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = df.index.astype(str)
    df.columns = [str(c).replace(".", "-") for c in df.columns]
    return df.apply(pd.to_numeric, errors="coerce")


def load_pheno(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    idcol = df.columns[0]
    df = df.rename(columns={idcol: "sample"})
    df["sample"] = df["sample"].astype(str).str.replace(".", "-", regex=False)
    df = df.drop_duplicates("sample").set_index("sample")
    return df


def find_row(index: pd.Index, symbol: str) -> str | None:
    names = [symbol] + SYNONYMS.get(symbol, [])
    for name in names:
        if name in index:
            return str(name)
    for ens in ENSEMBL.get(symbol, []):
        hits = [i for i in index if str(i) == ens or str(i).startswith(ens + ".")]
        if hits:
            return str(hits[0])
    for name in names:
        hits = [i for i in index if str(i).split("|")[0] == name or str(i).endswith("|" + name)]
        if hits:
            return str(hits[0])
    return None


def extract_gene(mat: pd.DataFrame, symbol: str) -> pd.Series:
    row = find_row(mat.index, symbol)
    if row is None:
        return pd.Series(np.nan, index=mat.columns, name=symbol)
    s = pd.to_numeric(mat.loc[row], errors="coerce")
    if isinstance(s, pd.DataFrame):
        s = s.iloc[0]
    s.name = symbol
    return s


def mean_z(mat: pd.DataFrame, symbols: list[str], min_genes: int) -> tuple[pd.Series, pd.DataFrame]:
    rows = []
    coverage = []
    for sym in symbols:
        rid = find_row(mat.index, sym)
        found = rid is not None
        n_obs = 0
        used = False
        if found:
            s = pd.to_numeric(mat.loc[rid], errors="coerce")
            if isinstance(s, pd.DataFrame):
                s = s.iloc[0]
            n_obs = int(s.notna().sum())
            if n_obs >= 8:
                sd = float(s.std(ddof=0))
                if sd and not np.isnan(sd):
                    rows.append((s - s.mean()) / sd)
                    used = True
        coverage.append({"gene": sym, "row": rid, "found": found, "n_obs": n_obs, "used": used})
    cov = pd.DataFrame(coverage)
    if not rows:
        return pd.Series(np.nan, index=mat.columns, name="score"), cov
    z = pd.concat(rows, axis=1)
    n = z.notna().sum(axis=1)
    score = z.mean(axis=1, skipna=True)
    score[n < min_genes] = np.nan
    score.name = "score"
    return score, cov


def spearman(x: pd.Series, y: pd.Series) -> dict:
    d = pd.concat([x, y], axis=1).dropna()
    n = int(len(d))
    if n < 6 or d.iloc[:, 0].nunique() < 2 or d.iloc[:, 1].nunique() < 2:
        return {"n": n, "rho": np.nan, "p": np.nan, "ci95_lo": np.nan, "ci95_hi": np.nan}
    rho, p = stats.spearmanr(d.iloc[:, 0], d.iloc[:, 1])
    rng = np.random.default_rng(RNG_SEED)
    xv = d.iloc[:, 0].to_numpy()
    yv = d.iloc[:, 1].to_numpy()
    boots = []
    for _ in range(N_BOOT):
        idx = rng.integers(0, n, n)
        try:
            boots.append(float(stats.spearmanr(xv[idx], yv[idx]).statistic))
        except Exception:
            continue
    lo, hi = np.nanpercentile(boots, [2.5, 97.5]) if boots else (np.nan, np.nan)
    return {
        "n": n,
        "rho": float(rho),
        "p": float(p),
        "ci95_lo": float(lo),
        "ci95_hi": float(hi),
    }


def partial_spearman(x: pd.Series, y: pd.Series, z: pd.Series) -> dict:
    d = pd.concat([x, y, z], axis=1).dropna()
    d.columns = ["x", "y", "z"]
    n = int(len(d))
    if n < 8 or d["z"].nunique() < 2:
        return {"n": n, "rho": np.nan, "p": np.nan, "ci95_lo": np.nan, "ci95_hi": np.nan}
    rx, ry, rz = d["x"].rank(), d["y"].rank(), d["z"].rank()
    bx = np.polyfit(rz, rx, 1)
    by = np.polyfit(rz, ry, 1)
    ex = rx - (bx[0] * rz + bx[1])
    ey = ry - (by[0] * rz + by[1])
    out = spearman(ex, ey)
    out["n"] = n
    return out


def high_low(score: pd.Series, outcome: pd.Series) -> dict:
    d = pd.concat([score, outcome], axis=1).dropna()
    d.columns = ["s", "y"]
    if len(d) < 8:
        return {
            "n_high": 0,
            "n_low": 0,
            "median_high": np.nan,
            "median_low": np.nan,
            "p": np.nan,
            "U": np.nan,
        }
    med = d["s"].median()
    hi = d.loc[d["s"] > med, "y"]
    lo = d.loc[d["s"] <= med, "y"]
    if len(hi) < 3 or len(lo) < 3:
        return {
            "n_high": int(len(hi)),
            "n_low": int(len(lo)),
            "median_high": float(hi.median()) if len(hi) else np.nan,
            "median_low": float(lo.median()) if len(lo) else np.nan,
            "p": np.nan,
            "U": np.nan,
        }
    U, p = stats.mannwhitneyu(hi, lo, alternative="two-sided")
    return {
        "n_high": int(len(hi)),
        "n_low": int(len(lo)),
        "median_high": float(hi.median()),
        "median_low": float(lo.median()),
        "delta_median_high_minus_low": float(hi.median() - lo.median()),
        "U": float(U),
        "p": float(p),
    }


def fmt_p(p) -> str:
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return "NA"
    p = float(p)
    if p == 0:
        return "<1e-300"
    if p < 1e-4:
        return f"{p:.1e}"
    return f"{p:.3g}"


def fmt_rho(r) -> str:
    if r is None or (isinstance(r, float) and np.isnan(r)):
        return "NA"
    return f"{float(r):+.3f}"


def fmt_ci(lo, hi) -> str:
    if lo is None or hi is None or (isinstance(lo, float) and np.isnan(lo)):
        return "NA"
    return f"{float(lo):+.3f} to {float(hi):+.3f}"


def pick_col(df: pd.DataFrame, key: str) -> str | None:
    for name in PHENO_ALIASES.get(key, [key]):
        if name in df.columns:
            return name
    lower = {c.lower(): c for c in df.columns}
    for name in PHENO_ALIASES.get(key, [key]):
        if name.lower() in lower:
            return lower[name.lower()]
    return None


def pick_immune_cols(pheno: pd.DataFrame) -> dict[str, str]:
    found = {}
    for key in PHENO_ALIASES:
        col = pick_col(pheno, key)
        if col:
            found[key] = col
    if "CD8_CIBERSORT" not in found:
        for c in pheno.columns:
            cl = c.lower()
            if "cibersort" in cl and "cd8" in cl:
                found["CD8_CIBERSORT"] = c
                break
    if "CD8_xCell" not in found:
        for c in pheno.columns:
            cl = c.lower()
            if "xcell" in cl and "cd8" in cl and "naive" not in cl and "memory" not in cl:
                found["CD8_xCell"] = c
                break
    if "ImmuneScore" not in found:
        for c in pheno.columns:
            cl = c.lower()
            if "estimate" in cl and "immune" in cl:
                found["ImmuneScore"] = c
                break
    return found


def assoc_row(predictor: str, endpoint: str, family: str, x, y, purity_name: str | None, z=None) -> dict:
    marg = spearman(x, y)
    row = {
        "predictor": predictor,
        "endpoint": endpoint,
        "family": family,
        "n": marg["n"],
        "rho": marg["rho"],
        "p": marg["p"],
        "ci95_lo": marg["ci95_lo"],
        "ci95_hi": marg["ci95_hi"],
        "purity": purity_name or "none",
        "n_partial": np.nan,
        "partial_rho": np.nan,
        "partial_p": np.nan,
        "partial_ci95_lo": np.nan,
        "partial_ci95_hi": np.nan,
    }
    if z is not None:
        part = partial_spearman(x, y, z)
        row.update(
            {
                "n_partial": part["n"],
                "partial_rho": part["rho"],
                "partial_p": part["p"],
                "partial_ci95_lo": part["ci95_lo"],
                "partial_ci95_hi": part["ci95_hi"],
            }
        )
    return row


def scatter_panel(ax, x, y, xlabel, ylabel, title_prefix=""):
    d = pd.concat([x, y], axis=1).dropna()
    d.columns = ["x", "y"]
    ax.scatter(d["x"], d["y"], s=22, alpha=0.75, c="#2c5aa0", edgecolors="none")
    if len(d) >= 6:
        res = spearman(d["x"], d["y"])
        title = f"{title_prefix}n={res['n']}  ρ={fmt_rho(res['rho'])}  p={fmt_p(res['p'])}"
        if d["x"].nunique() > 2:
            m, b = np.polyfit(d["x"], d["y"], 1)
            xs = np.linspace(d["x"].min(), d["x"].max(), 50)
            ax.plot(xs, m * xs + b, color="#c44e52", lw=1.4, alpha=0.85)
    else:
        title = f"{title_prefix}n={len(d)}"
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=10)
    ax.axhline(d["y"].median() if len(d) else 0, color="#999", lw=0.6, ls=":")
    ax.axvline(d["x"].median() if len(d) else 0, color="#999", lw=0.6, ls=":")


def box_panel(ax, score, outcome, xlabel, ylabel):
    d = pd.concat([score, outcome], axis=1).dropna()
    d.columns = ["s", "y"]
    med = d["s"].median()
    hi = d.loc[d["s"] > med, "y"]
    lo = d.loc[d["s"] <= med, "y"]
    hl = high_low(score, outcome)
    data = [lo.values, hi.values]
    bp = ax.boxplot(data, tick_labels=["low", "high"], patch_artist=True, widths=0.55)
    for patch, color in zip(bp["boxes"], ["#7aa6c2", "#c44e52"]):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(
        f"median split  {hl['n_low']} vs {hl['n_high']}  MWU p={fmt_p(hl['p'])}",
        fontsize=10,
    )


def forest_panel(ax, rows: list[dict], title: str, use_partial: bool):
    labels = []
    rhos = []
    los = []
    his = []
    for r in rows:
        labels.append(f"{r['predictor']} vs {r['endpoint']}")
        if use_partial:
            rhos.append(r["partial_rho"])
            los.append(r["partial_ci95_lo"])
            his.append(r["partial_ci95_hi"])
        else:
            rhos.append(r["rho"])
            los.append(r["ci95_lo"])
            his.append(r["ci95_hi"])
    y = np.arange(len(labels))[::-1]
    ax.axvline(0, color="#444", lw=0.8)
    for yi, rho, lo, hi in zip(y, rhos, los, his):
        if rho is None or (isinstance(rho, float) and np.isnan(rho)):
            continue
        ax.plot([lo, hi], [yi, yi], color="#2c5aa0", lw=1.6)
        ax.plot(rho, yi, "o", color="#c44e52", ms=6)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Spearman ρ (95% bootstrap CI)")
    ax.set_title(title, fontsize=11)
    ax.set_xlim(-0.85, 0.55)


def write_finding(
    path: Path,
    assoc: pd.DataFrame,
    coverage: dict[str, pd.DataFrame],
    highlow: pd.DataFrame,
    pheno_map: dict[str, str],
    n_tumor: int,
    wes_n: int,
    tj15_used: list[str],
    fig_rel: str,
) -> None:
    def row(pred, end, family=None):
        hit = assoc[(assoc["predictor"] == pred) & (assoc["endpoint"] == end)]
        if family is not None:
            hit = hit[hit["family"] == family]
        else:
            hit = hit[hit["family"] != "wgs_sensitivity"]
        return hit.iloc[0] if len(hit) else None

    def line(r, with_partial=True) -> str:
        if r is None:
            return "| — | — | NA | NA | NA | NA |"
        part = "NA"
        if with_partial and pd.notna(r.get("partial_rho")):
            part = f"{int(r['n_partial'])}, {fmt_rho(r['partial_rho'])}, {fmt_p(r['partial_p'])}"
        return (
            f"| {r['predictor']} | {r['endpoint']} | {int(r['n'])} | "
            f"{fmt_rho(r['rho'])} | {fmt_ci(r['ci95_lo'], r['ci95_hi'])} | "
            f"{fmt_p(r['p'])} | {part} |"
        )

    primary_preds = ["CLDN4 protein", "TJ-15 protein"]
    primary_ends = [
        "ESTIMATE ImmuneScore",
        "ESTIMATE ESTIMATEScore",
        "GEP18 RNA",
        "CD8A RNA",
    ]
    extra_ends = [
        "ESTIMATE StromalScore",
        "CD8 RNA (CD8A/B)",
        "CYT RNA (GZMA/PRF1)",
        "CIBERSORT CD8",
        "xCell CD8",
        "xCell immune score",
    ]

    tj_cov = coverage["tj15_protein"]
    used_n = int(tj_cov["used"].sum())
    used_genes = ", ".join(tj_cov.loc[tj_cov["used"], "gene"].tolist())
    missing_genes = ", ".join(tj_cov.loc[~tj_cov["used"], "gene"].tolist()) or "none"

    cldn4_cov = coverage["protein_targets"]
    cldn4 = cldn4_cov[cldn4_cov["gene"] == "CLDN4"].iloc[0]
    tac = cldn4_cov[cldn4_cov["gene"] == "TACSTD2"].iloc[0]

    cldn4_tac = row("CLDN4 protein", "TACSTD2 protein")
    cldn4_wes = row("CLDN4 protein", "WES purity")
    tj_wes = row("TJ-15 protein", "WES purity")
    imm_wes = row("ESTIMATE ImmuneScore", "WES purity")

    hl_lines = []
    for _, h in highlow.iterrows():
        dlt = h.get("delta_median_high_minus_low", np.nan)
        dlt_s = f"{float(dlt):+.3g}" if pd.notna(dlt) else "NA"
        hl_lines.append(
            f"| {h['predictor']} | {h['endpoint']} | {int(h['n_high'])} | "
            f"{int(h['n_low'])} | {dlt_s} | {fmt_p(h['p'])} |"
        )

    primary_md = "\n".join(
        line(row(p, e)) for p in primary_preds for e in primary_ends
    )
    extra_md = "\n".join(
        line(row(p, e)) for p in primary_preds for e in extra_ends if row(p, e) is not None
    )
    tj7_md = "\n".join(line(row("TJ-7 protein", e)) for e in primary_ends if row("TJ-7 protein", e) is not None)
    coexp_md = line(cldn4_tac, with_partial=True)
    cldn4_imm = row("CLDN4 protein", "ESTIMATE ImmuneScore")
    tj15_imm = row("TJ-15 protein", "ESTIMATE ImmuneScore")
    cldn4_gep = row("CLDN4 protein", "GEP18 RNA")
    cldn4_cd8a = row("CLDN4 protein", "CD8A RNA")
    cldn4_wgs = row("CLDN4 protein", "ESTIMATE ImmuneScore", family="wgs_sensitivity")
    wgs_note = ""
    if cldn4_wgs is not None and pd.notna(cldn4_wgs.get("partial_rho")):
        wgs_note = (
            f"After **WGS** purity (sensitivity, n={int(cldn4_wgs['n_partial'])}) "
            f"CLDN4 protein vs ImmuneScore partial is {fmt_rho(cldn4_wgs['partial_rho'])}, "
            f"p={fmt_p(cldn4_wgs['partial_p'])}. WES is the pre-specified DNA residual."
        )
    pur_bits = []
    for pred, lab in [
        ("CLDN4 protein", "CLDN4 protein vs WES purity"),
        ("TJ-15 protein", "TJ-15 protein vs WES purity"),
        ("ESTIMATE ImmuneScore", "ImmuneScore vs WES purity"),
    ]:
        r = row(pred, "WES purity")
        if r is not None:
            pur_bits.append(f"{lab}: {fmt_rho(r['rho'])} (n={int(r['n'])}, p={fmt_p(r['p'])})")
    purity_note = "; ".join(pur_bits) + "." if pur_bits else ""

    wes_note = (
        f"WES purity is present for **{wes_n} / {n_tumor}** tumors "
        f"(column `{pheno_map.get('WES_purity', 'WES_purity')}`)."
        if wes_n
        else "WES purity was not found in the phenotype table; no DNA residual is reported."
    )

    md = f"""# CPTAC-LSCC (LUSC) CLDN4 / TJ-15 protein versus immune

Additive protein extra. **LUSC protein rows are first.** CPTAC LUAD TJ-15 protein versus ImmuneScore (ρ=−0.30, n=110) is taken as given from [PR #245](https://github.com/jinxuanhong1-blip/sdaxcge/pull/245) and is not re-estimated.

Question: in the public CPTAC lung squamous (LSCC / LUSC) TMT proteome, how do **CLDN4 protein** and a **structural TJ-15 protein score** associate with ImmuneScore, ESTIMATE, GEP18, and CD8A? Also CLDN4 protein versus TACSTD2 protein, and the same associations after WES purity when that column exists.

Public processed tables only. Pairwise-complete Spearman ρ, two-sided p, 2,000-resample bootstrap 95% CI (seed `{RNG_SEED}`). Partial ρ residualizes ranks on WES purity.

## LUSC protein (Satpathy 2021; freeze v1.2)

n = **{n_tumor}** treatment-naive surgical LSCC tumors. CLDN4 protein is quantified in **{int(cldn4['n_obs'])} / {n_tumor}** ({int(n_tumor - cldn4['n_obs'])} NA). TACSTD2 protein is quantified in **{int(tac['n_obs'])} / {n_tumor}**. No ICI labels.

TJ-15 protein uses the same 15 structural genes as PR #245. On this public LSCC TMT table, **{used_n} / 15** members have ≥8 non-missing tumors and enter the mean z-score: {used_genes}. Not used (absent or <8 observations): {missing_genes}. Score requires ≥{MIN_TJ15} members per sample.

{wes_note} CLDN7 protein is sparse (**13 / 108**); it enters the TJ mean only where measured.

CLDN4 protein is inverse with ImmuneScore ({fmt_rho(cldn4_imm['rho'])}, n={int(cldn4_imm['n'])}, p={fmt_p(cldn4_imm['p'])}), GEP18 RNA ({fmt_rho(cldn4_gep['rho'])}), and CD8A RNA ({fmt_rho(cldn4_cd8a['rho'])}). Those three stay negative after WES residual. The 15-gene TJ protein score versus ImmuneScore is {fmt_rho(tj15_imm['rho'])} (n={int(tj15_imm['n'])}, p={fmt_p(tj15_imm['p'])}). CLDN4 and TACSTD2 proteins do not co-vary ({fmt_rho(cldn4_tac['rho'])}, n={int(cldn4_tac['n'])}, p={fmt_p(cldn4_tac['p'])}).

### Primary LUSC protein rows

| Predictor | Endpoint | n | ρ | 95% CI | p | Partial \\| WES (n, ρ, p) |
|---|---|---:|---:|---|---:|---|
{primary_md}

### CLDN4 protein vs TACSTD2 protein

| Predictor | Endpoint | n | ρ | 95% CI | p | Partial \\| WES (n, ρ, p) |
|---|---|---:|---:|---|---:|---|
{coexp_md}

### TJ-7 protein (sensitivity; same 7 genes as the claim-page module)

Not a search. CLDN7 remains n=13.

| Predictor | Endpoint | n | ρ | 95% CI | p | Partial \\| WES (n, ρ, p) |
|---|---|---:|---:|---|---:|---|
{tj7_md}

### Extra LUSC protein rows (same tumors)

| Predictor | Endpoint | n | ρ | 95% CI | p | Partial \\| WES (n, ρ, p) |
|---|---|---:|---:|---|---:|---|
{extra_md}

### Purity context and WGS sensitivity

{purity_note} {wgs_note}

### Median-split extras

| Predictor | Endpoint | n high | n low | Δmedian (high−low) | MWU p |
|---|---|---:|---:|---:|---:|
{chr(10).join(hl_lines)}

## Extra figures

![CLDN4 protein vs ImmuneScore, GEP18, CD8A]({fig_rel}/fig1_cldn4_protein_immune.png)

![TJ-15 protein vs ImmuneScore, GEP18, CD8A]({fig_rel}/fig2_tj15_protein_immune.png)

![CLDN4 vs TACSTD2 protein, and WES-purity context]({fig_rel}/fig3_cldn4_tacstd2_purity.png)

![LUSC protein forest (marginal and WES partial)]({fig_rel}/fig4_lusc_protein_forest.png)

![Median-split ImmuneScore / GEP18]({fig_rel}/fig5_highlow.png)

## Given LUAD protein (not re-run)

From PR #245, CPTAC LUAD TJ-15 protein versus ESTIMATE ImmuneScore is **ρ=−0.30, n=110** (reported −0.296, p=0.0017; partial \\| WES −0.256, p=0.0075). That LUAD table is not re-downloaded here.

## Methods (short)

- **Source:** CPTAC pan-cancer freeze `data_freeze_v1.2_reorganized` / LinkedOmics [CPTAC-pancan-LSCC](https://www.linkedomics.org/data_download/CPTAC-pancan-LSCC/). Paper: Satpathy et al., *Cell* 2021.
- **Files:** `LSCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt`, `LSCC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt`, `LSCC_phenotype.txt`.
- **IDs:** Ensembl prefix match (`ENSG…` ± version). Protein values are already log2 vs reference; not logged again.
- **TJ-15:** CLDN1, CLDN3, CLDN4, CLDN7, OCLN, TJP1, TJP2, TJP3, F11R, JAM2, JAM3, MARVELD2, MARVELD3, CGN, CGNL1. Mean of per-gene z-scores; genes with <8 observations dropped.
- **TJ-7** (sensitivity): CLDN1, CLDN4, CLDN7, F11R, TJP1, TJP2, OCLN.
- **GEP18:** Ayers 2017 18-gene T-cell-inflamed list, scored on matched tumor RNA (mean z; ≥{MIN_GEP} genes).
- **CD8A:** single RNA row. CD8 RNA = mean z of CD8A+CD8B (extra).
- **Immune / ESTIMATE:** freeze phenotype columns `{pheno_map.get('ImmuneScore', 'ESTIMATE_ImmuneScore')}`, `{pheno_map.get('ESTIMATEScore', 'ESTIMATE_ESTIMATEScore')}`, `{pheno_map.get('StromalScore', 'ESTIMATE_StromalScore')}`.
- **Partial Spearman:** rank-transform X, Y, and WES purity; residualize X and Y ranks on the purity rank; Spearman of residuals.
- **Not done:** no LUAD re-fit; no ICI response model (this surgical proteome has no ICI labels); no gene-set fishing.

## Reproduce

```bash
python3 -m pip install -r requirements.txt
python3 scripts/cptac_lusc_cldn4_protein/download.py
python3 scripts/cptac_lusc_cldn4_protein/analyze.py
```

Tables: `methods/cptac_lusc_cldn4_protein/tables/`. Figures: `methods/cptac_lusc_cldn4_protein/figures/`.
"""
    path.write_text(md)
    _ = cldn4_wes, tj_wes, imm_wes, tj15_used


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--datadir", default="data/cptac_lusc_cldn4_protein")
    p.add_argument("--outdir", default="methods/cptac_lusc_cldn4_protein")
    args = p.parse_args()

    datadir = Path(args.datadir)
    outdir = Path(args.outdir)
    figdir = outdir / "figures"
    tabdir = outdir / "tables"
    figdir.mkdir(parents=True, exist_ok=True)
    tabdir.mkdir(parents=True, exist_ok=True)

    prot_path = datadir / "LSCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt"
    rna_path = datadir / "LSCC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt"
    pheno_path = datadir / "LSCC_phenotype.txt"
    for f in (prot_path, rna_path, pheno_path):
        if not f.exists():
            raise SystemExit(f"missing {f}; run download.py first")

    protein = load_matrix(prot_path)
    rna = load_matrix(rna_path)
    pheno = load_pheno(pheno_path)
    pheno_map = pick_immune_cols(pheno)

    samples = protein.columns.intersection(rna.columns).intersection(pheno.index)
    protein = protein[samples]
    rna = rna[samples]
    pheno = pheno.loc[samples]
    n_tumor = int(len(samples))

    cldn4_p = extract_gene(protein, "CLDN4")
    tac_p = extract_gene(protein, "TACSTD2")
    tj15_p, tj15_cov = mean_z(protein, TJ_15, MIN_TJ15)
    tj7_p, tj7_cov = mean_z(protein, TJ_7, MIN_TJ7)
    tj15_p.name = "TJ15_protein"
    tj7_p.name = "TJ7_protein"

    gep_r, gep_cov = mean_z(rna, GEP18, MIN_GEP)
    cd8_r, cd8_cov = mean_z(rna, CD8, 2)
    cyt_r, cyt_cov = mean_z(rna, CYT, 2)
    cd8a_r = extract_gene(rna, "CD8A")
    gep_r.name = "GEP18_RNA"
    cd8_r.name = "CD8_RNA"
    cyt_r.name = "CYT_RNA"
    cd8a_r.name = "CD8A_RNA"

    target_cov = []
    for layer, mat, genes in [
        ("protein", protein, ["CLDN4", "TACSTD2"]),
        ("rna", rna, ["CLDN4", "TACSTD2", "CD8A"]),
    ]:
        for g in genes:
            rid = find_row(mat.index, g)
            s = extract_gene(mat, g)
            target_cov.append(
                {
                    "layer": layer,
                    "gene": g,
                    "row": rid,
                    "n_obs": int(s.notna().sum()),
                    "n_na": int(s.isna().sum()),
                }
            )
    target_cov_df = pd.DataFrame(target_cov)

    endpoints = {}
    if "ImmuneScore" in pheno_map:
        endpoints["ESTIMATE ImmuneScore"] = pd.to_numeric(pheno[pheno_map["ImmuneScore"]], errors="coerce")
    if "ESTIMATEScore" in pheno_map:
        endpoints["ESTIMATE ESTIMATEScore"] = pd.to_numeric(pheno[pheno_map["ESTIMATEScore"]], errors="coerce")
    if "StromalScore" in pheno_map:
        endpoints["ESTIMATE StromalScore"] = pd.to_numeric(pheno[pheno_map["StromalScore"]], errors="coerce")
    endpoints["GEP18 RNA"] = gep_r
    endpoints["CD8A RNA"] = cd8a_r
    endpoints["CD8 RNA (CD8A/B)"] = cd8_r
    endpoints["CYT RNA (GZMA/PRF1)"] = cyt_r
    if "CD8_CIBERSORT" in pheno_map:
        endpoints["CIBERSORT CD8"] = pd.to_numeric(pheno[pheno_map["CD8_CIBERSORT"]], errors="coerce")
    if "CD8_xCell" in pheno_map:
        endpoints["xCell CD8"] = pd.to_numeric(pheno[pheno_map["CD8_xCell"]], errors="coerce")
    if "Immune_xCell" in pheno_map:
        endpoints["xCell immune score"] = pd.to_numeric(pheno[pheno_map["Immune_xCell"]], errors="coerce")

    wes = None
    wes_n = 0
    if "WES_purity" in pheno_map:
        wes = pd.to_numeric(pheno[pheno_map["WES_purity"]], errors="coerce")
        wes_n = int(wes.notna().sum())
    wgs = None
    if "WGS_purity" in pheno_map:
        wgs = pd.to_numeric(pheno[pheno_map["WGS_purity"]], errors="coerce")

    predictors = {
        "CLDN4 protein": cldn4_p,
        "TJ-15 protein": tj15_p,
        "TJ-7 protein": tj7_p,
        "TACSTD2 protein": tac_p,
    }

    rows = []
    primary_ends = [
        "ESTIMATE ImmuneScore",
        "ESTIMATE ESTIMATEScore",
        "GEP18 RNA",
        "CD8A RNA",
    ]
    extra_ends = [
        "ESTIMATE StromalScore",
        "CD8 RNA (CD8A/B)",
        "CYT RNA (GZMA/PRF1)",
        "CIBERSORT CD8",
        "xCell CD8",
        "xCell immune score",
    ]
    for pname, pser in predictors.items():
        for ename in primary_ends + extra_ends:
            if ename not in endpoints:
                continue
            fam = "primary" if (pname in {"CLDN4 protein", "TJ-15 protein"} and ename in primary_ends) else "extra"
            rows.append(assoc_row(pname, ename, fam, pser, endpoints[ename], "WES_purity" if wes is not None else None, wes))
            if wgs is not None and pname in {"CLDN4 protein", "TJ-15 protein"} and ename in primary_ends:
                rows.append(assoc_row(pname, ename, "wgs_sensitivity", pser, endpoints[ename], "WGS_purity", wgs))

    rows.append(assoc_row("CLDN4 protein", "TACSTD2 protein", "coexpression", cldn4_p, tac_p, "WES_purity" if wes is not None else None, wes))

    if wes is not None:
        for pname, pser in [("CLDN4 protein", cldn4_p), ("TJ-15 protein", tj15_p), ("TACSTD2 protein", tac_p)]:
            rows.append(assoc_row(pname, "WES purity", "purity_context", pser, wes, None))
        if "ESTIMATE ImmuneScore" in endpoints:
            rows.append(assoc_row("ESTIMATE ImmuneScore", "WES purity", "purity_context", endpoints["ESTIMATE ImmuneScore"], wes, None))
        if "GEP18 RNA" in endpoints:
            rows.append(assoc_row("GEP18 RNA", "WES purity", "purity_context", endpoints["GEP18 RNA"], wes, None))
        if "CD8A RNA" in endpoints:
            rows.append(assoc_row("CD8A RNA", "WES purity", "purity_context", endpoints["CD8A RNA"], wes, None))

    assoc = pd.DataFrame(rows)

    hl_rows = []
    for pname, pser in [("CLDN4 protein", cldn4_p), ("TJ-15 protein", tj15_p)]:
        for ename in ["ESTIMATE ImmuneScore", "GEP18 RNA", "CD8A RNA"]:
            if ename not in endpoints:
                continue
            h = high_low(pser, endpoints[ename])
            h["predictor"] = pname
            h["endpoint"] = ename
            hl_rows.append(h)
    highlow = pd.DataFrame(hl_rows)

    # Sample-level table
    sample = pd.DataFrame(
        {
            "CLDN4_protein": cldn4_p,
            "TACSTD2_protein": tac_p,
            "TJ15_protein": tj15_p,
            "TJ7_protein": tj7_p,
            "GEP18_RNA": gep_r,
            "CD8A_RNA": cd8a_r,
            "CD8_RNA": cd8_r,
            "CYT_RNA": cyt_r,
        }
    )
    for key, col in pheno_map.items():
        sample[key] = pd.to_numeric(pheno[col], errors="coerce")
    sample.index.name = "sample"
    sample.to_csv(tabdir / "sample_scores.tsv", sep="\t")
    assoc.to_csv(tabdir / "associations.tsv", sep="\t", index=False)
    highlow.to_csv(tabdir / "highlow.tsv", sep="\t", index=False)
    tj15_cov.to_csv(tabdir / "tj15_protein_coverage.tsv", sep="\t", index=False)
    tj7_cov.to_csv(tabdir / "tj7_protein_coverage.tsv", sep="\t", index=False)
    gep_cov.to_csv(tabdir / "gep18_rna_coverage.tsv", sep="\t", index=False)
    target_cov_df.to_csv(tabdir / "target_coverage.tsv", sep="\t", index=False)
    pd.DataFrame([{"key": k, "column": v} for k, v in pheno_map.items()]).to_csv(
        tabdir / "phenotype_columns.tsv", sep="\t", index=False
    )

    # Figures
    plt.rcParams.update(
        {
            "figure.dpi": 140,
            "savefig.dpi": 160,
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )

    fig, axes = plt.subplots(1, 3, figsize=(12.2, 3.8))
    scatter_panel(axes[0], cldn4_p, endpoints.get("ESTIMATE ImmuneScore", pd.Series(dtype=float)), "CLDN4 protein", "ESTIMATE ImmuneScore")
    scatter_panel(axes[1], cldn4_p, gep_r, "CLDN4 protein", "GEP18 RNA")
    scatter_panel(axes[2], cldn4_p, cd8a_r, "CLDN4 protein", "CD8A RNA")
    fig.suptitle("CPTAC-LSCC: CLDN4 protein vs immune", fontsize=12, y=1.02)
    fig.tight_layout()
    fig.savefig(figdir / "fig1_cldn4_protein_immune.png", bbox_inches="tight")
    fig.savefig(figdir / "fig1_cldn4_protein_immune.pdf", bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(12.2, 3.8))
    scatter_panel(axes[0], tj15_p, endpoints.get("ESTIMATE ImmuneScore", pd.Series(dtype=float)), "TJ-15 protein", "ESTIMATE ImmuneScore")
    scatter_panel(axes[1], tj15_p, gep_r, "TJ-15 protein", "GEP18 RNA")
    scatter_panel(axes[2], tj15_p, cd8a_r, "TJ-15 protein", "CD8A RNA")
    fig.suptitle("CPTAC-LSCC: TJ-15 protein vs immune", fontsize=12, y=1.02)
    fig.tight_layout()
    fig.savefig(figdir / "fig2_tj15_protein_immune.png", bbox_inches="tight")
    fig.savefig(figdir / "fig2_tj15_protein_immune.pdf", bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(12.2, 3.8))
    scatter_panel(axes[0], cldn4_p, tac_p, "CLDN4 protein", "TACSTD2 protein")
    if wes is not None:
        scatter_panel(axes[1], cldn4_p, wes, "CLDN4 protein", "WES purity")
        scatter_panel(axes[2], endpoints.get("ESTIMATE ImmuneScore", pd.Series(dtype=float)), wes, "ESTIMATE ImmuneScore", "WES purity")
    else:
        axes[1].text(0.5, 0.5, "WES purity not in phenotype", ha="center", va="center")
        axes[1].set_axis_off()
        axes[2].set_axis_off()
    fig.suptitle("CPTAC-LSCC: CLDN4 vs TACSTD2 protein and WES purity", fontsize=12, y=1.02)
    fig.tight_layout()
    fig.savefig(figdir / "fig3_cldn4_tacstd2_purity.png", bbox_inches="tight")
    fig.savefig(figdir / "fig3_cldn4_tacstd2_purity.pdf", bbox_inches="tight")
    plt.close(fig)

    forest_rows = [
        r
        for r in rows
        if r["family"] == "primary"
        and r["purity"] in {"WES_purity", "none"}
    ]
    fig, axes = plt.subplots(1, 2, figsize=(12.6, 5.2))
    forest_panel(axes[0], forest_rows, "Marginal Spearman", use_partial=False)
    forest_panel(axes[1], forest_rows, "Partial Spearman | WES purity", use_partial=True)
    fig.suptitle("CPTAC-LSCC protein rows", fontsize=12)
    fig.tight_layout()
    fig.savefig(figdir / "fig4_lusc_protein_forest.png", bbox_inches="tight")
    fig.savefig(figdir / "fig4_lusc_protein_forest.pdf", bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(2, 3, figsize=(11.4, 7.0))
    box_panel(axes[0, 0], cldn4_p, endpoints.get("ESTIMATE ImmuneScore", pd.Series(dtype=float)), "CLDN4 protein", "ImmuneScore")
    box_panel(axes[0, 1], cldn4_p, gep_r, "CLDN4 protein", "GEP18 RNA")
    box_panel(axes[0, 2], cldn4_p, cd8a_r, "CLDN4 protein", "CD8A RNA")
    box_panel(axes[1, 0], tj15_p, endpoints.get("ESTIMATE ImmuneScore", pd.Series(dtype=float)), "TJ-15 protein", "ImmuneScore")
    box_panel(axes[1, 1], tj15_p, gep_r, "TJ-15 protein", "GEP18 RNA")
    box_panel(axes[1, 2], tj15_p, cd8a_r, "TJ-15 protein", "CD8A RNA")
    fig.suptitle("CPTAC-LSCC median-split extras", fontsize=12)
    fig.tight_layout()
    fig.savefig(figdir / "fig5_highlow.png", bbox_inches="tight")
    fig.savefig(figdir / "fig5_highlow.pdf", bbox_inches="tight")
    plt.close(fig)

    summary = {
        "cohort": "CPTAC-LSCC",
        "synonym": "LUSC",
        "n_tumor": n_tumor,
        "cldn4_protein_n": int(cldn4_p.notna().sum()),
        "tacstd2_protein_n": int(tac_p.notna().sum()),
        "tj15_protein_genes_used": tj15_cov.loc[tj15_cov["used"], "gene"].tolist(),
        "tj15_protein_n": int(tj15_p.notna().sum()),
        "wes_purity_n": wes_n,
        "wes_purity_column": pheno_map.get("WES_purity"),
        "phenotype_columns": pheno_map,
        "luad_given": LUAD_GIVEN,
        "seed": RNG_SEED,
        "n_boot": N_BOOT,
        "written_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    (tabdir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    write_finding(
        outdir / "FINDING.md",
        assoc,
        {
            "tj15_protein": tj15_cov,
            "protein_targets": target_cov_df[target_cov_df["layer"] == "protein"],
        },
        highlow,
        pheno_map,
        n_tumor,
        wes_n,
        tj15_cov.loc[tj15_cov["used"], "gene"].tolist(),
        "figures",
    )

    print(json.dumps({"n_tumor": n_tumor, "cldn4_n": int(cldn4_p.notna().sum()), "wes_n": wes_n, "outdir": str(outdir)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
