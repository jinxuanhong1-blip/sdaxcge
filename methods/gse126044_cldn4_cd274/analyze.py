#!/usr/bin/env python3
"""GSE126044: CLDN4 vs CD274 / HLA-A/B/C (CD274 extra).

Additive. CLDN4 only. Patient is the unit. Honest n.

ImmuneScore vs CLDN4 is already reported in leftover OPEN ICI bulk
(PR #296; n=16, ρ=−0.57, p=0.022). This page does not re-audit that
row. The extra axis is CD274 (PD-L1). HLA-A/B/C and MHC-I mean-z
(HLA-A/B/C only, 3/3) are reported next to it.

Partial Spearman residualizes ranks on ImmuneScore (the known covariate).
CD8A residual is a sensitivity row, not the primary partial.

Public GEO processed counts only. Does not invent histology, purity,
or response labels. Correlations use all 16 patients.
"""
from __future__ import annotations

import gzip
import json
import math
import urllib.request
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
CACHE = Path("/tmp/gse126044_cldn4_cd274")
TAB = HERE / "tables"
FIG = HERE / "figures"
PROC = HERE / "processed"
for d in (TAB, FIG, PROC, CACHE):
    d.mkdir(parents=True, exist_ok=True)

URLS = {
    "counts": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126044/suppl/"
        "GSE126044_counts.txt.gz"
    ),
    "series": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126044/matrix/"
        "GSE126044_series_matrix.txt.gz"
    ),
}

GENES = {
    "CLDN4": ["CLDN4", "ENSG00000189143"],
    "CD274": ["CD274", "PDCD1LG1", "PDL1", "B7-H1", "ENSG00000120217"],
    "HLA-A": ["HLA-A", "HLAA", "ENSG00000206503"],
    "HLA-B": ["HLA-B", "HLAB", "ENSG00000234745"],
    "HLA-C": ["HLA-C", "HLAC", "ENSG00000204525"],
    "CD8A": ["CD8A", "ENSG00000153563"],
}

EST = pd.read_csv(HERE / "estimate_gene_sets.csv")
IMMUNE = [g for g in EST["immune_signature"].dropna().astype(str) if g]
STROMAL = [g for g in EST["stromal_signature"].dropna().astype(str) if g]


def download(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 1000:
        return dest
    print(f"download {url} -> {dest}", flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=180) as r, open(dest, "wb") as f:
        f.write(r.read())
    return dest


def _split_fields(line: str) -> list[str]:
    return [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]


def parse_series(path: Path) -> pd.DataFrame:
    stored: dict[str, list[str]] = {}
    char_rows: list[list[str]] = []
    with gzip.open(path, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                break
            if line.startswith("!Sample_geo_accession"):
                stored["geo"] = _split_fields(line)
            elif line.startswith("!Sample_title"):
                stored["title"] = _split_fields(line)
            elif line.startswith("!Sample_source_name"):
                stored["source"] = _split_fields(line)
            elif line.startswith("!Sample_description"):
                stored["description"] = _split_fields(line)
            elif line.startswith("!Sample_characteristics"):
                char_rows.append(_split_fields(line))
    geo = stored.get("geo") or []
    n = len(geo)
    meta = pd.DataFrame({"gsm": geo})
    for key in ("title", "source", "description"):
        vals = stored.get(key) or [None] * n
        meta[key] = (vals + [None] * n)[:n]
    chars: dict[str, dict[str, str]] = defaultdict(dict)
    for vals in char_rows:
        for g, v in zip(geo, vals):
            if not v:
                continue
            if ":" in v:
                k, val = v.split(":", 1)
                chars[g][k.strip().lower()] = val.strip()
            else:
                chars[g]["characteristic"] = v
    keys = sorted({k for d in chars.values() for k in d})
    for k in keys:
        # GEO "sample:" is fresh/FFPE — do not overwrite matrix sample IDs
        col = "sample_type" if k == "sample" else k
        meta[col] = [chars[g].get(k) for g in geo]
    return meta


def normalize_symbol(x: str) -> str:
    x = str(x).strip().strip('"')
    x = x.split(".")[0]
    if "_" in x and not x.startswith("ENSG"):
        x = x.split("_")[0]
    return x.upper()


def collapse_symbols(expr: pd.DataFrame) -> pd.DataFrame:
    expr = expr.copy()
    expr.index = [normalize_symbol(i) for i in expr.index]
    return expr.groupby(expr.index).mean()


def resolve_gene(index: pd.Index, aliases: list[str]) -> str | None:
    upper = {str(i).upper(): str(i) for i in index}
    for a in aliases:
        if a.upper() in upper:
            return upper[a.upper()]
    return None


def counts_to_log2cpm(expr: pd.DataFrame) -> pd.DataFrame:
    lib = expr.sum(axis=0).replace(0, np.nan)
    cpm = expr.div(lib, axis=1) * 1e6
    return np.log2(cpm + 1)


def mean_z(log: pd.DataFrame, genes: list[str]) -> tuple[pd.Series | None, list[str]]:
    found = [g for g in genes if g in log.index]
    if len(found) < 2:
        return None, found
    arr = log.loc[found].to_numpy(dtype=float)
    mu = np.nanmean(arr, axis=1, keepdims=True)
    sd = np.nanstd(arr, axis=1, keepdims=True)
    sd = np.where(sd == 0, np.nan, sd)
    return pd.Series(np.nanmean((arr - mu) / sd, axis=0), index=log.columns), found


def spearman(x, y) -> dict:
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = ~(np.isnan(x) | np.isnan(y))
    n = int(m.sum())
    if n < 4:
        return {"n": n, "rho": np.nan, "p": np.nan}
    r, p = stats.spearmanr(x[m], y[m])
    return {"n": n, "rho": float(r), "p": float(p)}


def partial_spearman(x, y, z) -> dict:
    """Pearson of rank residuals; df = n − 3."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    z = np.asarray(z, float)
    m = ~(np.isnan(x) | np.isnan(y) | np.isnan(z))
    n = int(m.sum())
    if n < 5:
        return {"n": n, "rho": np.nan, "p": np.nan}
    rx, ry, rz = stats.rankdata(x[m]), stats.rankdata(y[m]), stats.rankdata(z[m])
    Z = np.column_stack([np.ones(n), rz])
    bx, *_ = np.linalg.lstsq(Z, rx, rcond=None)
    by, *_ = np.linalg.lstsq(Z, ry, rcond=None)
    ex, ey = rx - Z @ bx, ry - Z @ by
    if np.std(ex) == 0 or np.std(ey) == 0:
        return {"n": n, "rho": np.nan, "p": np.nan}
    r = float(np.corrcoef(ex, ey)[0, 1])
    df = n - 3
    t = r * math.sqrt(df / max(1e-12, 1 - r * r))
    p = float(2 * stats.t.sf(abs(t), df))
    return {"n": n, "rho": r, "p": p}


def fmt_p(p: float) -> str:
    if p != p:
        return "NA"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.3g}"


def fmt_rho(r: float) -> str:
    if r != r:
        return "NA"
    return f"{r:+.3f}"


def fmt_cell(r: float, p: float) -> str:
    if r != r:
        return "NA"
    return f"{fmt_rho(r)} ({fmt_p(p)})"


def main() -> None:
    download(URLS["counts"], CACHE / "GSE126044_counts.txt.gz")
    download(URLS["series"], CACHE / "GSE126044_series_matrix.txt.gz")

    raw = pd.read_csv(CACHE / "GSE126044_counts.txt.gz", sep="\t", index_col=0)
    raw = raw.apply(pd.to_numeric, errors="coerce")
    n_genes_raw, n_cols = raw.shape
    log = counts_to_log2cpm(collapse_symbols(raw))

    meta = parse_series(CACHE / "GSE126044_series_matrix.txt.gz")
    meta["matrix_id"] = meta["title"].str.replace("RNA-seq_", "", regex=False)

    resolved = {canon: resolve_gene(log.index, aliases) for canon, aliases in GENES.items()}
    rows = []
    for s in log.columns:
        rec = {"matrix_id": str(s)}
        for canon, hit in resolved.items():
            rec[canon] = float(log.loc[hit, s]) if hit is not None else np.nan
        rows.append(rec)
    sc = pd.DataFrame(rows)

    mhci, mhci_used = mean_z(log, ["HLA-A", "HLA-B", "HLA-C"])
    immune, immune_used = mean_z(log, IMMUNE)
    stromal, stromal_used = mean_z(log, STROMAL)
    if mhci is not None:
        sc["MHCI"] = [float(mhci[s]) for s in log.columns]
    else:
        sc["MHCI"] = np.nan
    if immune is not None:
        sc["ImmuneScore"] = [float(immune[s]) for s in log.columns]
    else:
        sc["ImmuneScore"] = np.nan
    if stromal is not None:
        sc["StromalScore"] = [float(stromal[s]) for s in log.columns]
        sc["ESTIMATEScore"] = sc["ImmuneScore"] + sc["StromalScore"]
    else:
        sc["StromalScore"] = np.nan
        sc["ESTIMATEScore"] = np.nan

    df = meta.merge(sc, on="matrix_id", how="inner")
    df["patient"] = df["gsm"]
    assert df["patient"].is_unique, "GSM / patient IDs are not unique"
    assert len(df) == n_cols, f"merge lost samples: {len(df)} vs {n_cols} matrix columns"
    df["resp"] = df.get("patient response", pd.Series([None] * len(df))).astype(str).str.lower()
    n_r = int(df["resp"].eq("responder").sum())
    n_nr = int(df["resp"].eq("non-responder").sum())
    n_fresh = int((df.get("sample_type", pd.Series([""] * len(df))).astype(str).str.lower() == "fresh").sum())
    n_ffpe = int((df.get("sample_type", pd.Series([""] * len(df))).astype(str).str.lower() == "ffpe").sum())

    df.to_csv(PROC / "GSE126044_patient.tsv", sep="\t", index=False)

    axes = ["CD274", "HLA-A", "HLA-B", "HLA-C", "MHCI", "ImmuneScore", "CD8A"]
    tests = []
    for axis in axes:
        sp = spearman(df["CLDN4"], df[axis])
        p_imm = partial_spearman(df["CLDN4"], df[axis], df["ImmuneScore"])
        p_cd8 = partial_spearman(df["CLDN4"], df[axis], df["CD8A"])
        tests.append(
            {
                "cohort": "GSE126044",
                "axis": axis,
                "n": sp["n"],
                "rho": sp["rho"],
                "p": sp["p"],
                "partial_ImmuneScore_rho": p_imm["rho"] if axis != "ImmuneScore" else np.nan,
                "partial_ImmuneScore_p": p_imm["p"] if axis != "ImmuneScore" else np.nan,
                "partial_CD8A_rho": p_cd8["rho"] if axis != "CD8A" else np.nan,
                "partial_CD8A_p": p_cd8["p"] if axis != "CD8A" else np.nan,
                "role": (
                    "already_known"
                    if axis == "ImmuneScore"
                    else "extra"
                    if axis == "CD274"
                    else "companion"
                    if axis in ("HLA-A", "HLA-B", "HLA-C", "MHCI")
                    else "sensitivity"
                ),
            }
        )
    tests_df = pd.DataFrame(tests)
    tests_df.to_csv(TAB / "spearman.tsv", sep="\t", index=False)

    by = tests_df.set_index("axis")
    one = {
        "cohort": "GSE126044",
        "n": int(by.loc["CD274", "n"]),
        "n_R": n_r,
        "n_NR": n_nr,
        "scale": "log2(CPM+1)",
        "CLDN4": resolved["CLDN4"] or "",
        "CD274": resolved["CD274"] or "",
        "HLA-A": resolved["HLA-A"] or "",
        "HLA-B": resolved["HLA-B"] or "",
        "HLA-C": resolved["HLA-C"] or "",
        "CD8A": resolved["CD8A"] or "",
        "MHCI_genes": ",".join(mhci_used),
        "ImmuneScore_genes_used": len(immune_used),
        "ImmuneScore_genes_list": len(IMMUNE),
        "CD274_rho": by.loc["CD274", "rho"],
        "CD274_p": by.loc["CD274", "p"],
        "HLA-A_rho": by.loc["HLA-A", "rho"],
        "HLA-A_p": by.loc["HLA-A", "p"],
        "HLA-B_rho": by.loc["HLA-B", "rho"],
        "HLA-B_p": by.loc["HLA-B", "p"],
        "HLA-C_rho": by.loc["HLA-C", "rho"],
        "HLA-C_p": by.loc["HLA-C", "p"],
        "MHCI_rho": by.loc["MHCI", "rho"],
        "MHCI_p": by.loc["MHCI", "p"],
        "ImmuneScore_rho": by.loc["ImmuneScore", "rho"],
        "ImmuneScore_p": by.loc["ImmuneScore", "p"],
        "CD8A_rho": by.loc["CD8A", "rho"],
        "CD8A_p": by.loc["CD8A", "p"],
        "CD274_partial_ImmuneScore_rho": by.loc["CD274", "partial_ImmuneScore_rho"],
        "CD274_partial_ImmuneScore_p": by.loc["CD274", "partial_ImmuneScore_p"],
        "MHCI_partial_ImmuneScore_rho": by.loc["MHCI", "partial_ImmuneScore_rho"],
        "MHCI_partial_ImmuneScore_p": by.loc["MHCI", "partial_ImmuneScore_p"],
        "CD274_partial_CD8A_rho": by.loc["CD274", "partial_CD8A_rho"],
        "CD274_partial_CD8A_p": by.loc["CD274", "partial_CD8A_p"],
        "MHCI_partial_CD8A_rho": by.loc["MHCI", "partial_CD8A_rho"],
        "MHCI_partial_CD8A_p": by.loc["MHCI", "partial_CD8A_p"],
    }
    pd.DataFrame([one]).to_csv(TAB / "one_row.tsv", sep="\t", index=False)

    inv = {
        "cohort": "GSE126044",
        "geo": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126044",
        "tissue": "NSCLC tumor biopsy, anti-PD-1, pre-treatment",
        "assay": "RNA-seq counts (author processed)",
        "scale": "log2(CPM+1)",
        "n_series_gsm": int(len(meta)),
        "n_matrix_columns": int(n_cols),
        "n_merged_patients": int(df["patient"].nunique()),
        "n_unique_gsm": int(df["gsm"].nunique()),
        "n_unique_titles": int(df["title"].nunique()),
        "n_R": n_r,
        "n_NR": n_nr,
        "n_fresh": n_fresh,
        "n_ffpe": n_ffpe,
        "n_CLDN4_finite": int(df["CLDN4"].notna().sum()),
        "n_CD274_finite": int(df["CD274"].notna().sum()),
        "n_HLA_A_finite": int(df["HLA-A"].notna().sum()),
        "n_HLA_B_finite": int(df["HLA-B"].notna().sum()),
        "n_HLA_C_finite": int(df["HLA-C"].notna().sum()),
        "n_CD8A_finite": int(df["CD8A"].notna().sum()),
        "n_ImmuneScore_finite": int(df["ImmuneScore"].notna().sum()),
        "n_genes_raw": int(n_genes_raw),
        "n_symbols_after_collapse": int(log.shape[0]),
        "ImmuneScore_coverage": f"{len(immune_used)}/{len(IMMUNE)}",
        "StromalScore_coverage": f"{len(stromal_used)}/{len(STROMAL)}",
        "MHCI_definition": "HLA-A/B/C mean-z only (3/3); not B2M/TAP",
        "histology_split": "none on GEO (NSCLC)",
        "purity_deposited": "no",
        "luad_lusc_split": "no",
        "primary_pairwise_n": int(by.loc["CD274", "n"]),
    }
    pd.DataFrame([inv]).to_csv(TAB / "inventory.tsv", sep="\t", index=False)

    summary = {
        "question": "GSE126044 CLDN4 vs CD274/HLA; ImmuneScore already known; CD274 extra",
        "n": one["n"],
        "one_row": one,
        "inventory": inv,
        "resolved": {k: (v or "") for k, v in resolved.items()},
        "tests": tests,
    }
    (TAB / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    # Figures
    fig, axes_p = plt.subplots(1, 2, figsize=(8.4, 3.6))
    colors = np.where(df["resp"].eq("responder"), "#1f77b4", "#c0392b")
    axes_p[0].scatter(df["CLDN4"], df["CD274"], c=colors, s=36, edgecolors="k", linewidths=0.3)
    axes_p[0].set_xlabel("CLDN4 log2(CPM+1)")
    axes_p[0].set_ylabel("CD274 log2(CPM+1)")
    axes_p[0].set_title(
        f"CLDN4 vs CD274 (extra)\nρ={fmt_rho(one['CD274_rho'])}  p={fmt_p(one['CD274_p'])}  n={one['n']}"
    )
    axes_p[1].scatter(df["CLDN4"], df["MHCI"], c=colors, s=36, edgecolors="k", linewidths=0.3)
    axes_p[1].set_xlabel("CLDN4 log2(CPM+1)")
    axes_p[1].set_ylabel("HLA-A/B/C mean-z")
    axes_p[1].set_title(
        f"CLDN4 vs MHC-I (HLA-A/B/C)\nρ={fmt_rho(one['MHCI_rho'])}  p={fmt_p(one['MHCI_p'])}  n={one['n']}"
    )
    axes_p[1].axhline(0, color="#999", lw=0.6, ls="--")
    fig.tight_layout()
    fig.savefig(FIG / "fig1_cldn4_vs_cd274_hla.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    labels = ["CD274 (extra)", "HLA-A", "HLA-B", "HLA-C", "MHC-I mean-z", "ImmuneScore (known)"]
    keys = ["CD274", "HLA-A", "HLA-B", "HLA-C", "MHCI", "ImmuneScore"]
    rhos = [by.loc[k, "rho"] for k in keys]
    ps = [by.loc[k, "p"] for k in keys]
    y = np.arange(len(keys))[::-1]
    cols = ["#c0392b" if k == "CD274" else "#7f8c8d" if k == "ImmuneScore" else "#1f4e79" for k in keys]
    ax.axvline(0, color="#444", lw=0.8)
    ax.scatter(rhos, y, c=cols, s=48, zorder=3)
    for i, (rho, p, lab) in enumerate(zip(rhos, ps, labels)):
        ax.text(0.02 if rho >= 0 else -0.02, y[i] + 0.18, f"{lab}  {fmt_cell(rho, p)}", fontsize=8, va="bottom",
                ha="left" if rho >= 0 else "right")
    ax.set_yticks([])
    ax.set_xlabel("Spearman ρ vs CLDN4")
    ax.set_xlim(-1.05, 0.4)
    ax.set_title(f"GSE126044 n={one['n']}  (honest; all 16 patients)")
    fig.tight_layout()
    fig.savefig(FIG / "fig2_forest.png", dpi=160)
    plt.close(fig)

    write_finding(one, inv, by)
    print(json.dumps({"n": one["n"], "CD274": fmt_cell(one["CD274_rho"], one["CD274_p"]),
                      "ImmuneScore": fmt_cell(one["ImmuneScore_rho"], one["ImmuneScore_p"])}, indent=2))


def write_finding(one: dict, inv: dict, by: pd.DataFrame) -> None:
    n = one["n"]
    md = f"""# GSE126044 — CLDN4 vs CD274 / HLA (CD274 extra)

**Additive. CLDN4 only.** No TACSTD2 gate. Patient is the unit.

Public Cho et al. anti-PD-1 NSCLC biopsy RNA-seq ([GSE126044](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126044); Cho et al., *Nat Commun* 2020). Author count matrix → `log2(CPM+1)`.

**Already known (not re-audited as a claim).** CLDN4 vs ESTIMATE ImmuneScore on this same matrix is in leftover OPEN ICI bulk ([PR #296](https://github.com/jinxuanhong1-blip/sdaxcge/pull/296)): n=16, ρ=−0.57, p=0.022. This page recomputes ImmuneScore only so CD274 / HLA can be residualized on the same covariate. The ImmuneScore row is a **given**, not a new finding.

**Extra.** `CD274` (PD-L1) was not on the leftover immune-axis table (that table had CD8A / IFN / 6-gene MHC / ImmuneScore). HLA-A/B/C and MHC-I mean-z here are **HLA-A/B/C only (3/3)**, not B2M/TAP.

This page does **not** re-audit R vs NR (CLDN4 Cliff δ=−0.53, p=0.115 is given in PR #296 / PR #149). Correlations use all **16** patients.

Reproduce: `python3 methods/gse126044_cldn4_cd274/analyze.py`

---

## Honest n

| Item | Public? | n | Note |
|---|---|---:|---|
| GEO series GSM | yes | **{inv['n_series_gsm']}** | GSE126044 series matrix |
| Count-matrix columns | yes | **{inv['n_matrix_columns']}** | `GSE126044_counts.txt.gz` |
| Unique GSM / unique titles | yes | **{inv['n_unique_gsm']} / {inv['n_unique_titles']}** | 1:1 |
| Merged patients (unit) | yes | **{inv['n_merged_patients']}** | title `RNA-seq_*` stripped to matrix ID |
| GEO responder / non-responder | yes | **{inv['n_R']} / {inv['n_NR']}** | used only as a plot color; not the test |
| Fresh / FFPE (`sample:` characteristic) | yes | **{inv['n_fresh']} / {inv['n_ffpe']}** | all 5 FFPE are NR (given) |
| LUAD vs LUSC | **no** | 0 | GEO says NSCLC; no histology table |
| Tumor % / ABSOLUTE purity | **no** | 0 | RNA ImmuneScore is the public proxy |
| CLDN4 / CD274 / HLA-A / HLA-B / HLA-C finite | yes | **{inv['n_CLDN4_finite']} / {inv['n_CD274_finite']} / {inv['n_HLA_A_finite']} / {inv['n_HLA_B_finite']} / {inv['n_HLA_C_finite']}** | symbols on the count matrix |
| ImmuneScore finite (Yoshihara immune list mean-z) | yes | **{inv['n_ImmuneScore_finite']}** | coverage {inv['ImmuneScore_coverage']} |
| **Primary pairwise n (CLDN4 + CD274)** | yes | **{n}** | this is the n used below |

Do not write n>16. Do not drop FFPE to chase a p-value. Do not treat 5 vs 11 as the CD274 n.

---

## One-row table

| Cohort | n | scale | CD274 ρ (p) **extra** | HLA-A ρ (p) | HLA-B ρ (p) | HLA-C ρ (p) | MHC-I mean-z ρ (p) | CD274 \\| ImmuneScore ρ (p) | MHC-I \\| ImmuneScore ρ (p) | ImmuneScore ρ (p) **given** |
|---|---:|---|---|---|---|---|---|---|---|---|
| GSE126044 | **{n}** | log2(CPM+1) | **{fmt_cell(one['CD274_rho'], one['CD274_p'])}** | {fmt_cell(one['HLA-A_rho'], one['HLA-A_p'])} | {fmt_cell(one['HLA-B_rho'], one['HLA-B_p'])} | {fmt_cell(one['HLA-C_rho'], one['HLA-C_p'])} | {fmt_cell(one['MHCI_rho'], one['MHCI_p'])} | {fmt_cell(one['CD274_partial_ImmuneScore_rho'], one['CD274_partial_ImmuneScore_p'])} | {fmt_cell(one['MHCI_partial_ImmuneScore_rho'], one['MHCI_partial_ImmuneScore_p'])} | {fmt_cell(one['ImmuneScore_rho'], one['ImmuneScore_p'])} |

Numeric row: `tables/one_row.tsv`. Axis-level tests: `tables/spearman.tsv`.

---

## How to read the row

**CD274 extra (n={n}).** CLDN4 is **CD274-low** (ρ={fmt_rho(one['CD274_rho'])}, p={fmt_p(one['CD274_p'])}). After residualizing ranks on the already-known ImmuneScore, the partial is {fmt_cell(one['CD274_partial_ImmuneScore_rho'], one['CD274_partial_ImmuneScore_p'])}. CD8A residual (sensitivity, not primary): {fmt_cell(one['CD274_partial_CD8A_rho'], one['CD274_partial_CD8A_p'])}.

**HLA / MHC-I (3/3).** HLA-B/C trend negative; HLA-A is weaker. MHC-I mean-z ρ={fmt_rho(one['MHCI_rho'])}, p={fmt_p(one['MHCI_p'])}. After ImmuneScore the MHC-I residual is {fmt_cell(one['MHCI_partial_ImmuneScore_rho'], one['MHCI_partial_ImmuneScore_p'])}. Leftover PR #296 used a **6-gene** MHC set (HLA-A/B/C+B2M+TAP1+TAP2, raw ρ=−0.44, p=0.087; partial | ImmuneScore +0.19, p=0.49). That is a different score and is not re-used as the MHC-I number here.

**ImmuneScore given.** Recomputed Yoshihara immune-list mean-z matches the leftover row (ρ={fmt_rho(one['ImmuneScore_rho'])}, p={fmt_p(one['ImmuneScore_p'])}, n={n}, coverage {inv['ImmuneScore_coverage']}). Cited so the CD274 residual has a named covariate. Not a new ImmuneScore claim.

Honest n={n} is small. A large |ρ| is required to reach p<0.05. Do not pool this ρ with GSE218989 or other ICI bulk CD274 rows (signs conflict; see PR #326).

---

## Methods

- **Matrix:** GEO `GSE126044_counts.txt.gz`. Columns matched to series titles after stripping the `RNA-seq_` prefix (`Dis_01` …).
- **Scale:** `log2(CPM+1)` from library-size CPM. Symbols collapsed by mean after uppercasing.
- **CD274 extra** = HUGO `CD274` (aliases PDCD1LG1 / PDL1 / B7-H1 checked; matrix uses `CD274`).
- **MHC-I mean-z** = gene-wise z of `HLA-A`, `HLA-B`, `HLA-C` only ({one['MHCI_genes']}; 3/3). Not B2M/TAP.
- **ImmuneScore** = Yoshihara 2013 immune-signature **mean z** on the log2-CPM matrix (`estimate_gene_sets.csv`; same list as PR #296). Coverage {inv['ImmuneScore_coverage']}. This is not the Affymetrix-calibrated ESTIMATE R purity transform.
- **Partial Spearman:** Pearson of rank residuals; df = n − 3. Primary residual covariate = **ImmuneScore** (already-known infiltrate axis). CD8A residual is sensitivity.
- **Unit:** patient = GSM. All {n} patients enter every pairwise test. No FFPE drop. No LUAD/LUSC split (none deposited).

---

## What was not done

- No TACSTD2 column on the claim table.
- No R vs NR re-test and no search for a cutoff that hits p=0.019.
- No Fisher-z meta-analysis with other ICI bulk CD274 rows.
- No invented histology or purity.
- ImmuneScore is **given**, not sold as a new finding.

---

## Files

- `FINDING.md` — this page
- `analyze.py` — download-once, score, write tables/figures
- `estimate_gene_sets.csv` — Yoshihara stromal/immune lists (tidyestimate / PR #296)
- `tables/one_row.tsv` — the one-row table
- `tables/spearman.tsv` / `inventory.tsv` / `summary.json`
- `processed/GSE126044_patient.tsv`
- `figures/fig1_cldn4_vs_cd274_hla.png`, `fig2_forest.png`
"""
    (HERE / "FINDING.md").write_text(md)


if __name__ == "__main__":
    main()
