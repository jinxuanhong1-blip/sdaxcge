#!/usr/bin/env python3
"""CLDN4 vs CD274 and HLA-A/B/C on four public lung ICI bulk series.

Additive. CLDN4 only. Patient is the unit. Honest n / Spearman ρ / p.
Partial Spearman residualizes on CD8A ranks when CD8A is on the matrix.

Cohorts: GSE218989, GSE126044, GSE166449, GSE182328.
Public GEO processed matrices only. Does not invent ICI labels.
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
CACHE = Path("/tmp/cldn4_cd274_ici")
TAB = HERE / "tables"
FIG = HERE / "figures"
PROC = HERE / "processed"
for d in (TAB, FIG, PROC, CACHE):
    d.mkdir(parents=True, exist_ok=True)

URLS = {
    "GSE218989_TPM": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE218nnn/GSE218989/suppl/"
        "GSE218989_SMC_KAIST_TPM_matrix_pc_rmdupli.txt.gz"
    ),
    "GSE218989_series": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE218nnn/GSE218989/matrix/"
        "GSE218989_series_matrix.txt.gz"
    ),
    "GSE126044_counts": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126044/suppl/"
        "GSE126044_counts.txt.gz"
    ),
    "GSE126044_series": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126044/matrix/"
        "GSE126044_series_matrix.txt.gz"
    ),
    "GSE166449_TPM": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166449/suppl/"
        "GSE166449_Raw_gene_TPM_matrix.txt.gz"
    ),
    "GSE166449_series": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166449/matrix/"
        "GSE166449_series_matrix.txt.gz"
    ),
    "GSE182328_counts": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE182nnn/GSE182328/suppl/"
        "GSE182328_Gene_counts_matrix.txt.gz"
    ),
    "GSE182328_series": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE182nnn/GSE182328/matrix/"
        "GSE182328_series_matrix.txt.gz"
    ),
}

# Locked symbols. Ensembl used only as a fallback if the matrix is Ensembl-keyed.
GENES = {
    "CLDN4": ["CLDN4", "ENSG00000189143"],
    "CD274": ["CD274", "PDCD1LG1", "PDL1", "B7-H1", "ENSG00000120217"],
    "HLA-A": ["HLA-A", "HLAA", "ENSG00000206503"],
    "HLA-B": ["HLA-B", "HLAB", "ENSG00000234745"],
    "HLA-C": ["HLA-C", "HLAC", "ENSG00000204525"],
    "CD8A": ["CD8A", "ENSG00000153563"],
}
AXES = ["CD274", "HLA-A", "HLA-B", "HLA-C", "MHCI"]


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
        meta[k] = [chars[g].get(k) for g in geo]
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
    expr = expr.groupby(expr.index).mean()
    return expr


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


def to_log(expr: pd.DataFrame, is_counts: bool) -> pd.DataFrame:
    expr = collapse_symbols(expr)
    expr = expr.apply(pd.to_numeric, errors="coerce")
    if is_counts:
        return counts_to_log2cpm(expr)
    mx = float(np.nanmax(expr.to_numpy()))
    if mx <= 20:
        return expr
    return np.log2(expr.clip(lower=0) + 1)


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
    """Rank residual Spearman of x vs y after z (CD8A)."""
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


def mean_z(log: pd.DataFrame, genes: list[str]) -> pd.Series | None:
    found = [g for g in genes if g in log.index]
    if len(found) < 2:
        return None
    arr = log.loc[found].to_numpy(dtype=float)
    mu = np.nanmean(arr, axis=1, keepdims=True)
    sd = np.nanstd(arr, axis=1, keepdims=True)
    sd = np.where(sd == 0, np.nan, sd)
    return pd.Series(np.nanmean((arr - mu) / sd, axis=0), index=log.columns)


def extract_patient_table(log: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    resolved = {canon: resolve_gene(log.index, aliases) for canon, aliases in GENES.items()}
    rows = []
    for s in log.columns:
        rec = {"sample": str(s)}
        for canon, hit in resolved.items():
            rec[canon] = float(log.loc[hit, s]) if hit is not None else np.nan
        rows.append(rec)
    out = pd.DataFrame(rows)
    hla_hits = [resolved[g] for g in ("HLA-A", "HLA-B", "HLA-C") if resolved[g] is not None]
    mz = mean_z(log, hla_hits) if hla_hits else None
    if mz is not None:
        out["MHCI"] = [float(mz[s]) for s in log.columns]
    else:
        out["MHCI"] = np.nan
    cov = {k: (v if v is not None else "") for k, v in resolved.items()}
    cov["MHCI_genes"] = ",".join(hla_hits)
    return out, cov


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


def load_expr(path: Path, index_col: int = 0) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", index_col=index_col, compression="infer")


def score_cohort(
    cohort: str,
    expr: pd.DataFrame,
    meta: pd.DataFrame,
    sample_key: str,
    is_counts: bool,
    scale: str,
    tissue: str,
    assay: str,
    ici_note: str,
) -> tuple[pd.DataFrame, list[dict], dict]:
    log = to_log(expr, is_counts=is_counts)
    sc, cov = extract_patient_table(log)
    df = meta.merge(sc, left_on=sample_key, right_on="sample", how="inner")
    if df.empty:
        # try GSM columns if titles fail
        if "gsm" in meta.columns and set(meta["gsm"]).intersection(set(sc["sample"])):
            df = meta.merge(sc, left_on="gsm", right_on="sample", how="inner")
    if "patient" not in df.columns:
        df["patient"] = df["gsm"] if "gsm" in df.columns else df["sample"]
    assert df["patient"].is_unique, f"{cohort} patient IDs not unique after merge"
    df.to_csv(PROC / f"{cohort}_patient.tsv", sep="\t", index=False)

    tests = []
    for axis in AXES:
        if axis not in df or df[axis].notna().sum() < 4:
            tests.append(
                {
                    "cohort": cohort,
                    "axis": axis,
                    "n": int(df[axis].notna().sum()) if axis in df else 0,
                    "rho": np.nan,
                    "p": np.nan,
                    "partial_cd8a_rho": np.nan,
                    "partial_cd8a_p": np.nan,
                    "cd8a_present": bool(df["CD8A"].notna().sum() >= 5) if "CD8A" in df else False,
                    "note": "axis missing" if axis not in df or df[axis].notna().sum() < 4 else "",
                }
            )
            continue
        sp = spearman(df["CLDN4"], df[axis])
        if df["CD8A"].notna().sum() >= 5:
            psp = partial_spearman(df["CLDN4"], df[axis], df["CD8A"])
        else:
            psp = {"n": int(df["CD8A"].notna().sum()), "rho": np.nan, "p": np.nan}
        tests.append(
            {
                "cohort": cohort,
                "axis": axis,
                "n": sp["n"],
                "rho": sp["rho"],
                "p": sp["p"],
                "partial_cd8a_rho": psp["rho"],
                "partial_cd8a_p": psp["p"],
                "cd8a_present": df["CD8A"].notna().sum() >= 5,
                "note": "",
            }
        )
    inv = {
        "cohort": cohort,
        "tissue": tissue,
        "assay": assay,
        "scale": scale,
        "n_patients": int(df["patient"].nunique()),
        "n_merged": int(len(df)),
        "CLDN4": cov.get("CLDN4", ""),
        "CD274": cov.get("CD274", ""),
        "HLA-A": cov.get("HLA-A", ""),
        "HLA-B": cov.get("HLA-B", ""),
        "HLA-C": cov.get("HLA-C", ""),
        "CD8A": cov.get("CD8A", ""),
        "MHCI_genes": cov.get("MHCI_genes", ""),
        "ici_note": ici_note,
    }
    return df, tests, inv


def main() -> None:
    for key, url in URLS.items():
        download(url, CACHE / Path(url).name)

    all_tests: list[dict] = []
    inventory: list[dict] = []
    patients: dict[str, pd.DataFrame] = {}

    # GSE218989 — SMC-KAIST NSCLC ICI TPM, patient columns
    tpm = load_expr(CACHE / "GSE218989_SMC_KAIST_TPM_matrix_pc_rmdupli.txt.gz")
    meta218 = parse_series(CACHE / "GSE218989_series_matrix.txt.gz")
    # TPM columns are patient IDs (SMC__Pat*); series titles match
    meta218["sample_key"] = meta218["title"]
    # keep only columns that are in the TPM matrix
    keep = [c for c in tpm.columns if str(c) in set(meta218["title"]) or str(c).startswith("SMC")]
    if not keep:
        keep = list(tpm.columns)
    tpm = tpm[keep]
    meta218 = meta218[meta218["title"].isin(tpm.columns)].copy()
    if meta218.empty:
        meta218 = pd.DataFrame({"gsm": tpm.columns, "title": tpm.columns, "sample_key": tpm.columns})
    df, tests, inv = score_cohort(
        "GSE218989",
        tpm,
        meta218,
        "sample_key",
        is_counts=False,
        scale="log2(TPM+1)",
        tissue="NSCLC tumor (SMC-KAIST PD-1/PD-L1)",
        assay="RNA-seq TPM, protein-coding",
        ici_note="GEO treatment outcome on 355 patients; this page is CLDN4 vs CD274/HLA, not response",
    )
    patients["GSE218989"] = df
    all_tests.extend(tests)
    inventory.append(inv)

    # GSE126044 — counts; sample IDs are title without RNA-seq_
    counts = load_expr(CACHE / "GSE126044_counts.txt.gz")
    meta126 = parse_series(CACHE / "GSE126044_series_matrix.txt.gz")
    meta126["sample_key"] = meta126["title"].str.replace("RNA-seq_", "", regex=False)
    df, tests, inv = score_cohort(
        "GSE126044",
        counts,
        meta126,
        "sample_key",
        is_counts=True,
        scale="log2(CPM+1)",
        tissue="NSCLC tumor biopsy, anti-PD-1",
        assay="RNA-seq counts",
        ici_note="GEO responder / non-responder (5 / 11); correlations use all 16 patients",
    )
    patients["GSE126044"] = df
    all_tests.extend(tests)
    inventory.append(inv)

    # GSE166449 — TPM; sample IDs are GEO description
    tpm166 = load_expr(CACHE / "GSE166449_Raw_gene_TPM_matrix.txt.gz")
    meta166 = parse_series(CACHE / "GSE166449_series_matrix.txt.gz")
    meta166["sample_key"] = meta166["description"]
    df, tests, inv = score_cohort(
        "GSE166449",
        tpm166,
        meta166,
        "sample_key",
        is_counts=False,
        scale="log2(TPM+1)",
        tissue="advanced LUAD tumor, pembrolizumab",
        assay="RNA-seq TPM",
        ici_note="GEO titles Responder 7 / nonResponder 15; correlations use all 22 patients",
    )
    patients["GSE166449"] = df
    all_tests.extend(tests)
    inventory.append(inv)

    # GSE182328 — counts; sample IDs are titles
    counts182 = load_expr(CACHE / "GSE182328_Gene_counts_matrix.txt.gz")
    meta182 = parse_series(CACHE / "GSE182328_series_matrix.txt.gz")
    meta182["sample_key"] = meta182["title"]
    df, tests, inv = score_cohort(
        "GSE182328",
        counts182,
        meta182,
        "sample_key",
        is_counts=True,
        scale="log2(CPM+1)",
        tissue="NSCLC tumor, ICI-treated",
        assay="RNA-seq counts",
        ici_note="GEO has Akkermansia only (not RECIST/PFS/MPR); correlations use all 44 patients",
    )
    patients["GSE182328"] = df
    all_tests.extend(tests)
    inventory.append(inv)

    tests_df = pd.DataFrame(all_tests)
    inv_df = pd.DataFrame(inventory)
    tests_df.to_csv(TAB / "spearman.tsv", sep="\t", index=False)
    inv_df.to_csv(TAB / "inventory.tsv", sep="\t", index=False)

    # One row per cohort
    rows = []
    for cohort in ["GSE218989", "GSE126044", "GSE166449", "GSE182328"]:
        sub = tests_df[tests_df["cohort"] == cohort].set_index("axis")
        inv = inv_df[inv_df["cohort"] == cohort].iloc[0]
        rec = {
            "cohort": cohort,
            "n": int(sub.loc["CD274", "n"]) if "CD274" in sub.index else int(inv["n_patients"]),
            "scale": inv["scale"],
            "tissue": inv["tissue"],
            "CD274_on_matrix": bool(inv["CD274"]),
            "HLA_A_on_matrix": bool(inv["HLA-A"]),
            "HLA_B_on_matrix": bool(inv["HLA-B"]),
            "HLA_C_on_matrix": bool(inv["HLA-C"]),
            "CD8A_on_matrix": bool(inv["CD8A"]),
        }
        for axis in AXES:
            if axis in sub.index:
                rec[f"{axis}_rho"] = sub.loc[axis, "rho"]
                rec[f"{axis}_p"] = sub.loc[axis, "p"]
                rec[f"{axis}_partial_cd8a_rho"] = sub.loc[axis, "partial_cd8a_rho"]
                rec[f"{axis}_partial_cd8a_p"] = sub.loc[axis, "partial_cd8a_p"]
            else:
                rec[f"{axis}_rho"] = np.nan
                rec[f"{axis}_p"] = np.nan
                rec[f"{axis}_partial_cd8a_rho"] = np.nan
                rec[f"{axis}_partial_cd8a_p"] = np.nan
        rec["ici_note"] = inv["ici_note"]
        rows.append(rec)
    one = pd.DataFrame(rows)
    one.to_csv(TAB / "one_row.tsv", sep="\t", index=False)

    # Forest of raw Spearman
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.2), sharey=True)
    cohort_order = ["GSE218989", "GSE126044", "GSE166449", "GSE182328"]
    y = np.arange(len(cohort_order))
    colors = {"CD274": "#c0392b", "MHCI": "#1f4e79"}
    for ax, axis, title in (
        (axes[0], "CD274", "CLDN4 vs CD274"),
        (axes[1], "MHCI", "CLDN4 vs HLA-A/B/C mean-z"),
    ):
        sub = tests_df[tests_df["axis"] == axis].set_index("cohort").loc[cohort_order]
        ax.axvline(0, color="#888", lw=0.8)
        ax.errorbar(
            sub["rho"],
            y,
            fmt="o",
            color=colors[axis],
            capsize=0,
        )
        for i, (rho, p, n) in enumerate(zip(sub["rho"], sub["p"], sub["n"])):
            ax.text(0.02, i + 0.22, f"n={int(n)}  ρ={rho:+.2f}  p={fmt_p(p)}", fontsize=8, color="#333")
        ax.set_yticks(y)
        ax.set_yticklabels(cohort_order)
        ax.set_xlabel("Spearman ρ")
        ax.set_title(title)
        ax.set_xlim(-1.05, 1.05)
    fig.suptitle("CLDN4-only, patient unit, four ICI bulk series (not pooled)", fontsize=11)
    fig.tight_layout()
    fig.savefig(FIG / "fig1_forest_raw.png", dpi=160)
    plt.close(fig)

    # Partial vs raw for CD274 and MHCI
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.2), sharey=True)
    for ax, axis, title in (
        (axes[0], "CD274", "CLDN4 vs CD274 | CD8A"),
        (axes[1], "MHCI", "CLDN4 vs HLA-A/B/C mean-z | CD8A"),
    ):
        sub = tests_df[tests_df["axis"] == axis].set_index("cohort").loc[cohort_order]
        ax.axvline(0, color="#888", lw=0.8)
        ax.scatter(sub["rho"], y - 0.12, s=36, c="#7f8c8d", label="raw", zorder=3)
        ax.scatter(sub["partial_cd8a_rho"], y + 0.12, s=36, c="#16a085", label="partial | CD8A", zorder=3)
        for i, (rho, p) in enumerate(zip(sub["partial_cd8a_rho"], sub["partial_cd8a_p"])):
            ax.text(0.02, i + 0.28, f"partial ρ={fmt_rho(rho)}  p={fmt_p(p)}", fontsize=8)
        ax.set_yticks(y)
        ax.set_yticklabels(cohort_order)
        ax.set_xlabel("Spearman ρ")
        ax.set_title(title)
        ax.set_xlim(-1.05, 1.05)
        ax.legend(fontsize=8, loc="lower right")
    fig.suptitle("Partial Spearman residualizes ranks on CD8A", fontsize=11)
    fig.tight_layout()
    fig.savefig(FIG / "fig2_partial_cd8a.png", dpi=160)
    plt.close(fig)

    # Per-cohort scatters
    fig, axes = plt.subplots(2, 4, figsize=(12.4, 6.4))
    for j, cohort in enumerate(cohort_order):
        df = patients[cohort]
        for i, (axis, ylab) in enumerate((("CD274", "CD274"), ("MHCI", "HLA-A/B/C mean-z"))):
            ax = axes[i, j]
            x = np.asarray(df["CLDN4"], float)
            yv = np.asarray(df[axis], float)
            m = ~(np.isnan(x) | np.isnan(yv))
            ax.scatter(x[m], yv[m], s=18, c="#3d5a80", alpha=0.8)
            if m.sum() >= 3:
                lr = stats.linregress(x[m], yv[m])
                xs = np.linspace(np.nanmin(x[m]), np.nanmax(x[m]), 40)
                ax.plot(xs, lr.intercept + lr.slope * xs, color="#e07a5f", lw=1.2)
            sp = spearman(df["CLDN4"], df[axis])
            ax.set_title(f"{cohort} n={sp['n']}\nρ={fmt_rho(sp['rho'])} p={fmt_p(sp['p'])}", fontsize=8)
            if i == 1:
                ax.set_xlabel("CLDN4")
            if j == 0:
                ax.set_ylabel(ylab)
    fig.suptitle("CLDN4 vs CD274 / MHC-I (raw)", fontsize=11)
    fig.tight_layout()
    fig.savefig(FIG / "fig3_scatters.png", dpi=160)
    plt.close(fig)

    summary = {
        "question": "CLDN4 vs CD274 and HLA-A/B/C on four lung ICI bulk series; partial | CD8A",
        "unit": "patient",
        "gene": "CLDN4 only",
        "n_cohorts": 4,
        "one_row": one.to_dict(orient="records"),
        "inventory": inv_df.to_dict(orient="records"),
    }
    with open(TAB / "summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)

    print(one.to_string(index=False))
    print("\nWrote", TAB / "one_row.tsv")


if __name__ == "__main__":
    main()
