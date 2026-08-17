#!/usr/bin/env python3
"""GSE121810 leftover: CLDN4 vs immune / CD274 on the public HUGO counts.

Cloughesy et al., Nat Med 2019 (PMID 30742122). Recurrent GBM, pembrolizumab
neoadjuvant (arm A, on-treatment RNA) vs adjuvant (arm B, pretreatment RNA).

The GEO series matrix is public metadata only (0 expression rows). Tests use
the deposited supplementary HUGO count matrix
GSE121810_Prins.PD1NeoAdjv.Jul2018.HUGO.PtID.xlsx (29 tumors).

OS / PFS / RECIST are not scored here (prior B5 analog PR already tested
CLDN4 vs survival). This slice is CLDN4 vs immune axes and CD274 only.

Downloads stay under $GSE121810_CLDN4_DATA (default /tmp/gse121810) and are
not committed.
"""

from __future__ import annotations

import gzip
import json
import os
import re
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
TABLES = HERE / "tables"
FIGURES = HERE / "figures"
DATA = Path(os.environ.get("GSE121810_CLDN4_DATA", "/tmp/gse121810"))

MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE121nnn/GSE121810/"
    "matrix/GSE121810_series_matrix.txt.gz"
)
COUNTS_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE121nnn/GSE121810/"
    "suppl/GSE121810_Prins.PD1NeoAdjv.Jul2018.HUGO.PtID.xlsx"
)

# Locked bulk T/NK list used on leftover lung ICI FPKM (GSE248378).
# All 16 are rows on this HUGO matrix. KLRF1 is present but was never locked.
TNK_GENES = [
    "CD8A",
    "CD8B",
    "CD2",
    "CD3D",
    "CD3E",
    "CD3G",
    "NKG7",
    "GNLY",
    "PRF1",
    "GZMA",
    "GZMB",
    "GZMK",
    "KLRD1",
    "KLRK1",
    "KLRB1",
    "NCR1",
]

# Ayers IFN-γ 6-gene. IFNG is a row but at the floor (honest 6/6 present,
# 0/29 with CPM ≥ 1). The score still uses the six genes; FINDING.md flags
# that IFNG itself is empty.
IFN_AYERS = ["IFNG", "STAT1", "CXCL9", "CXCL10", "IDO1", "HLA-DRA"]
MHC_I = ["HLA-A", "HLA-B", "HLA-C"]
SINGLE = [
    "CLDN4",
    "TACSTD2",
    "CD274",
    "CD8A",
    "NKG7",
    "IFNG",
    "STAT1",
    "CXCL9",
    "CXCL10",
    "PTPRC",
    "EPCAM",
    "GFAP",
    "PDCD1",
    "MKI67",
]


def dl(dest: Path, url: str) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    print(f"[get] {url}", flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": "gse121810-cldn4/1.0"})
    with urllib.request.urlopen(req, timeout=180) as r, open(dest, "wb") as out:
        out.write(r.read())
    return dest


def parse_series_meta(path: Path) -> tuple[pd.DataFrame, dict, int]:
    series: dict[str, list[str]] = {}
    sample_fields: dict[str, list[str]] = {}
    n_expr_rows = 0
    in_table = False
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if in_table:
                if not line.startswith('"ID_REF"') and not line.startswith("ID_REF"):
                    n_expr_rows += 1
                continue
            if line.startswith("!Series_"):
                k = line.split("\t", 1)[0][8:]
                v = line.split("\t", 1)[1].strip().strip('"') if "\t" in line else ""
                series.setdefault(k, []).append(v)
            elif line.startswith("!Sample_"):
                key = line.split("\t", 1)[0][8:]
                vals = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
                if key in sample_fields:
                    sample_fields[key] = [
                        (a + " | " + b) if a and b else (a or b)
                        for a, b in zip(sample_fields[key], vals)
                    ]
                else:
                    sample_fields[key] = vals
    meta = pd.DataFrame(sample_fields)
    if "geo_accession" in meta.columns:
        meta = meta.set_index("geo_accession")
    return meta, {k: " ".join(v) if isinstance(v, list) else v for k, v in series.items()}, n_expr_rows


def mean_z(log: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, list[str]]:
    found = [g for g in genes if g in log.index]
    if not found:
        return pd.Series(np.nan, index=log.columns), []
    z = log.loc[found].astype(float).apply(
        lambda r: (r - r.mean()) / (r.std(ddof=0) or np.nan), axis=1
    )
    return z.mean(axis=0), found


def spearman(x, y) -> dict:
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    n = int(x.size)
    if n < 5:
        return {"n": n, "rho": np.nan, "p": np.nan}
    rho, p = stats.spearmanr(x, y)
    return {"n": n, "rho": float(rho), "p": float(p)}


def partial_spearman(x, y, z) -> dict:
    x = pd.Series(x, dtype=float)
    y = pd.Series(y, dtype=float)
    z = pd.Series(z, dtype=float)
    ok = x.notna() & y.notna() & z.notna()
    n = int(ok.sum())
    if n < 6:
        return {"n": n, "rho": np.nan, "p": np.nan}
    rx = x[ok].rank()
    ry = y[ok].rank()
    rz = z[ok].rank()
    rx = rx - np.polyval(np.polyfit(rz, rx, 1), rz)
    ry = ry - np.polyval(np.polyfit(rz, ry, 1), rz)
    rho, p = stats.pearsonr(rx, ry)
    return {"n": n, "rho": float(rho), "p": float(p)}


def cliffs_delta(a, b) -> float:
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if a.size == 0 or b.size == 0:
        return float("nan")
    gt = sum((ai > b).sum() for ai in a)
    lt = sum((ai < b).sum() for ai in a)
    return float((gt - lt) / (a.size * b.size))


def mwu(a, b) -> dict:
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if a.size < 2 or b.size < 2:
        return {
            "n_high": int(a.size),
            "n_low": int(b.size),
            "median_high": float(np.median(a)) if a.size else np.nan,
            "median_low": float(np.median(b)) if b.size else np.nan,
            "U": np.nan,
            "p": np.nan,
            "delta": np.nan,
        }
    U, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    return {
        "n_high": int(a.size),
        "n_low": int(b.size),
        "median_high": float(np.median(a)),
        "median_low": float(np.median(b)),
        "U": float(U),
        "p": float(p),
        "delta": cliffs_delta(a, b),
    }


def detectability(counts: pd.DataFrame, gene: str) -> dict:
    if gene not in counts.index:
        return {
            "gene": gene,
            "present_in_matrix": False,
            "n": int(counts.shape[1]),
            "nonzero_n": 0,
            "raw_median": np.nan,
            "cpm_median": np.nan,
            "cpm_q1": np.nan,
            "cpm_q3": np.nan,
            "cpm_max": np.nan,
            "n_cpm_ge_1": 0,
            "n_cpm_ge_5": 0,
        }
    raw = counts.loc[gene].astype(float)
    lib = counts.sum(axis=0).replace(0, np.nan)
    cpm = raw / lib * 1e6
    return {
        "gene": gene,
        "present_in_matrix": True,
        "n": int(raw.size),
        "nonzero_n": int((raw > 0).sum()),
        "raw_median": float(raw.median()),
        "cpm_median": float(cpm.median()),
        "cpm_q1": float(cpm.quantile(0.25)),
        "cpm_q3": float(cpm.quantile(0.75)),
        "cpm_max": float(cpm.max()),
        "n_cpm_ge_1": int((cpm >= 1).sum()),
        "n_cpm_ge_5": int((cpm >= 5).sum()),
    }


def scatter(path, x, y, xlab, ylab, title, rho, p, n):
    fig, ax = plt.subplots(figsize=(4.4, 3.8))
    ax.scatter(x, y, s=28, c="#1f4e79", alpha=0.85, edgecolors="none")
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() >= 3:
        coef = np.polyfit(x[ok], y[ok], 1)
        xs = np.linspace(np.nanmin(x[ok]), np.nanmax(x[ok]), 50)
        ax.plot(xs, np.polyval(coef, xs), color="#c0392b", lw=1.2)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(title, fontsize=10)
    ax.text(
        0.04,
        0.96,
        f"ρ={rho:.3f}\np={p:.3g}\nn={n}",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=8,
        family="monospace",
    )
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def forest(path, rows):
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    ys = np.arange(len(rows))[::-1]
    rhos = [r["rho"] for r in rows]
    ax.axvline(0, color="#888", lw=0.8)
    ax.scatter(rhos, ys, s=36, c="#1f4e79", zorder=3)
    ax.set_yticks(ys)
    ax.set_yticklabels(
        [f"{r['pair']}  n={r['n']}  p={r['p']:.3g}" for r in rows],
        fontsize=8,
        family="monospace",
    )
    ax.set_xlabel("Spearman ρ (CLDN4 vs axis)")
    ax.set_xlim(-1.05, 1.05)
    ax.set_title("GSE121810 CLDN4 vs immune / CD274 (n=29)", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def detect_bar(path, det_df: pd.DataFrame, genes: list[str]):
    sub = det_df.set_index("gene").loc[genes]
    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    xs = np.arange(len(genes))
    ax.bar(xs, sub["cpm_median"].to_numpy(), color="#1f4e79")
    ax.set_xticks(xs)
    ax.set_xticklabels(genes, rotation=40, ha="right", fontsize=8)
    ax.set_ylabel("median CPM")
    ax.set_title("GSE121810 detectability (median CPM)", fontsize=10)
    ax.set_yscale("log")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def patient_num(title: str) -> int:
    m = re.search(r"Patient\s+(\d+)", str(title))
    if not m:
        raise ValueError(f"cannot parse patient from {title!r}")
    return int(m.group(1))


def xlsx_patient_num(col: str) -> int:
    m = re.match(r"Pt(\d+)_", str(col))
    if not m:
        raise ValueError(f"cannot parse xlsx column {col!r}")
    return int(m.group(1))


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)

    matrix_path = dl(DATA / "GSE121810_series_matrix.txt.gz", MATRIX_URL)
    counts_path = dl(
        DATA / "GSE121810_Prins.PD1NeoAdjv.Jul2018.HUGO.PtID.xlsx", COUNTS_URL
    )

    meta, series, n_expr_rows = parse_series_meta(matrix_path)
    counts = pd.read_excel(counts_path, index_col=0)
    counts.index = counts.index.astype(str)
    counts = counts[~counts.index.str.startswith("__")]
    counts = counts.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    if counts.index.duplicated().any():
        counts = counts.groupby(level=0).sum()

    meta = meta.copy()
    meta["patient"] = meta["title"].map(patient_num)
    xlsx_pts = [xlsx_patient_num(c) for c in counts.columns]
    geo_pts = meta["patient"].tolist()
    if xlsx_pts != geo_pts:
        raise SystemExit(f"patient order mismatch xlsx={xlsx_pts} geo={geo_pts}")

    counts.columns = meta.index
    lib = counts.sum(axis=0)
    cpm = counts.div(lib, axis=1) * 1e6
    log = np.log2(cpm + 1.0)

    tnk, tnk_used = mean_z(log, TNK_GENES)
    ifn, ifn_used = mean_z(log, IFN_AYERS)
    mhc, mhc_used = mean_z(log, MHC_I)

    def extract_therapy(s: str) -> str:
        m = re.search(r"therapy:\s*([^|]+)", str(s), flags=re.I)
        return m.group(1).strip() if m else str(s).strip()

    therapy = meta["characteristics_ch1"].map(extract_therapy)
    arm = np.where(
        therapy.str.contains("neoadjuvant", case=False),
        "A_neoadjuvant",
        "B_adjuvant",
    )

    per = pd.DataFrame(index=meta.index)
    per["title"] = meta["title"].to_numpy()
    per["xlsx_id"] = [f"Pt{p}_{'A' if a.startswith('A') else 'B'}" for p, a in zip(meta["patient"], arm)]
    per["patient"] = meta["patient"].to_numpy()
    per["therapy"] = therapy.to_numpy()
    per["arm"] = arm
    per["libsize"] = lib.to_numpy()
    for g in ["CLDN4", "TACSTD2", "CD274", "CD8A", "NKG7", "IFNG", "PTPRC", "EPCAM", "GFAP", "PDCD1"]:
        per[g] = log.loc[g].to_numpy()
        per[f"{g}_CPM"] = cpm.loc[g].to_numpy()
        per[f"{g}_raw"] = counts.loc[g].to_numpy()
    per["TNK"] = tnk.to_numpy()
    per["IFN"] = ifn.to_numpy()
    per["MHC"] = mhc.to_numpy()
    per.to_csv(TABLES / "per_sample.tsv", sep="\t")

    det_rows = [detectability(counts, g) for g in SINGLE + TNK_GENES + IFN_AYERS + MHC_I]
    # unique genes, keep first
    seen = set()
    uniq = []
    for row in det_rows:
        if row["gene"] in seen:
            continue
        seen.add(row["gene"])
        uniq.append(row)
    det = pd.DataFrame(uniq)
    det.to_csv(TABLES / "detectability.tsv", sep="\t", index=False)

    coverage_rows = []
    for name, locked, used in [
        ("TNK", TNK_GENES, tnk_used),
        ("IFN_Ayers", IFN_AYERS, ifn_used),
        ("MHC_I", MHC_I, mhc_used),
    ]:
        coverage_rows.append(
            {
                "signature": name,
                "locked_n": len(locked),
                "present_n": len(used),
                "present": ",".join(used),
                "missing": ",".join([g for g in locked if g not in used]),
            }
        )
    for g in SINGLE:
        coverage_rows.append(
            {
                "signature": g,
                "locked_n": 1,
                "present_n": int(g in log.index),
                "present": g if g in log.index else "",
                "missing": "" if g in log.index else g,
            }
        )
    pd.DataFrame(coverage_rows).to_csv(TABLES / "gene_coverage.tsv", sep="\t", index=False)

    arm_counts = per["arm"].value_counts().to_dict()
    cldn4_cpm = per["CLDN4_CPM"]
    label_inv = pd.DataFrame(
        [
            {
                "field": "series_matrix_expression_rows",
                "public": "yes (empty table)",
                "n": n_expr_rows,
                "note": "RNA-seq series matrix; 0 features",
            },
            {
                "field": "supplementary_HUGO_samples",
                "public": "yes",
                "n": int(counts.shape[1]),
                "note": "GSE121810_Prins.PD1NeoAdjv.Jul2018.HUGO.PtID.xlsx",
            },
            {
                "field": "supplementary_HUGO_genes",
                "public": "yes",
                "n": int(counts.shape[0]),
                "note": "HUGO symbols; htseq-count GRCh38",
            },
            {
                "field": "GEO_GSM",
                "public": "yes",
                "n": int(len(meta)),
                "note": "GSM3447008–GSM3447036",
            },
            {
                "field": "therapy_arm",
                "public": "yes",
                "n": int(len(per)),
                "note": (
                    f"neoadjuvant(A)={arm_counts.get('A_neoadjuvant', 0)}; "
                    f"adjuvant(B)={arm_counts.get('B_adjuvant', 0)}"
                ),
            },
            {
                "field": "CLDN4_row",
                "public": "yes",
                "n": int(len(per)),
                "note": "present on HUGO matrix",
            },
            {
                "field": "CLDN4_CPM_ge_1",
                "public": "yes",
                "n": int((cldn4_cpm >= 1).sum()),
                "note": f"median CPM {float(cldn4_cpm.median()):.3f}; IQR {float(cldn4_cpm.quantile(0.25)):.3f}–{float(cldn4_cpm.quantile(0.75)):.3f}",
            },
            {
                "field": "CD274_row",
                "public": "yes",
                "n": int(len(per)),
                "note": "present; 27/29 CPM ≥ 1",
            },
            {
                "field": "IFNG_usable",
                "public": "row present, floor",
                "n": int((per["IFNG_CPM"] >= 1).sum()),
                "note": f"nonzero {int((per['IFNG_raw'] > 0).sum())}/29; max CPM {float(per['IFNG_CPM'].max()):.3f}",
            },
            {
                "field": "RECIST_ORR",
                "public": "no",
                "n": 0,
                "note": "not a GEO characteristic; not recovered here",
            },
            {
                "field": "OS_PFS",
                "public": "not used",
                "n": 0,
                "note": "prior B5 analog scored survival; this leftover does not",
            },
        ]
    )
    label_inv.to_csv(TABLES / "label_inventory.tsv", sep="\t", index=False)

    tests = []

    def add_spearman(pair, x, y, note="", subset="all"):
        sp = spearman(x, y)
        tests.append({"test": "spearman", "subset": subset, "pair": pair, "note": note, **sp})

    add_spearman("CLDN4 vs CD274", per["CLDN4"], per["CD274"], "primary leftover")
    add_spearman("CLDN4 vs TNK", per["CLDN4"], per["TNK"], f"mean-z {len(tnk_used)}/{len(TNK_GENES)}")
    add_spearman("CLDN4 vs IFN", per["CLDN4"], per["IFN"], f"Ayers mean-z {len(ifn_used)}/{len(IFN_AYERS)}; IFNG floor")
    add_spearman("CLDN4 vs CD8A", per["CLDN4"], per["CD8A"])
    add_spearman("CLDN4 vs NKG7", per["CLDN4"], per["NKG7"])
    add_spearman("CLDN4 vs PTPRC", per["CLDN4"], per["PTPRC"], "CD45 / immune content")
    add_spearman("CLDN4 vs MHC", per["CLDN4"], per["MHC"], f"HLA-A/B/C {len(mhc_used)}/3")
    add_spearman("CLDN4 vs PDCD1", per["CLDN4"], per["PDCD1"], "supporting")
    add_spearman("CLDN4 vs TACSTD2", per["CLDN4"], per["TACSTD2"], "companion; also floor")
    add_spearman("CLDN4 vs EPCAM", per["CLDN4"], per["EPCAM"], "epithelial; floor in GBM")
    add_spearman("CLDN4 vs GFAP", per["CLDN4"], per["GFAP"], "glial content")
    add_spearman("CD274 vs TNK", per["CD274"], per["TNK"], "immune-axis sanity")
    add_spearman("CD8A vs TNK", per["CD8A"], per["TNK"], "positive-control; CD8A is in TNK")
    add_spearman("TNK vs IFN", per["TNK"], per["IFN"], "immune-axis sanity")
    add_spearman("TNK vs PTPRC", per["TNK"], per["PTPRC"], "immune-axis sanity")

    for cov_name in ["PTPRC", "GFAP"]:
        for axis in ["CD274", "TNK", "IFN", "CD8A", "NKG7", "MHC"]:
            psp = partial_spearman(per["CLDN4"], per[axis], per[cov_name])
            tests.append(
                {
                    "test": f"partial_spearman_|{cov_name}",
                    "subset": "all",
                    "pair": f"CLDN4 vs {axis}",
                    "note": f"rank residual on {cov_name}",
                    **psp,
                }
            )

    med = float(per["CLDN4"].median())
    hi = per["CLDN4"] >= med
    lo = ~hi
    for axis in ["CD274", "TNK", "IFN", "CD8A", "NKG7", "MHC", "PTPRC"]:
        blk = mwu(per.loc[hi, axis], per.loc[lo, axis])
        tests.append(
            {
                "test": "median_split_MWU",
                "subset": "all",
                "pair": f"CLDN4-high vs low, {axis}",
                "note": f"high n={int(hi.sum())} (>= median log2(CPM+1) {med:.3f}); low n={int(lo.sum())}",
                "n": int(hi.sum() + lo.sum()),
                "rho": blk["delta"],
                "p": blk["p"],
                **{f"mwu_{k}": v for k, v in blk.items()},
            }
        )

    floor = per["CLDN4_CPM"] < 1
    above = ~floor
    if int(floor.sum()) >= 5:
        for axis in ["CD274", "TNK", "IFN", "CD8A"]:
            add_spearman(
                f"CLDN4 vs {axis}",
                per.loc[floor, "CLDN4"],
                per.loc[floor, axis],
                "CLDN4 CPM < 1 (floor-only)",
                subset="CLDN4_CPM_lt_1",
            )
    # n=3 above floor cannot support Spearman (n<5); record the empty test.
    tests.append(
        {
            "test": "spearman",
            "subset": "CLDN4_CPM_ge_1",
            "pair": "CLDN4 vs CD274",
            "note": "only 3 tumors have CLDN4 CPM ≥ 1; Spearman not run",
            "n": int(above.sum()),
            "rho": np.nan,
            "p": np.nan,
        }
    )

    for name, mask in [
        ("A_neoadjuvant", per["arm"] == "A_neoadjuvant"),
        ("B_adjuvant", per["arm"] == "B_adjuvant"),
    ]:
        sub = per.loc[mask]
        for axis in ["CD274", "TNK", "IFN", "CD8A"]:
            sp = spearman(sub["CLDN4"], sub[axis])
            tests.append(
                {
                    "test": "spearman_arm",
                    "subset": name,
                    "pair": f"CLDN4 vs {axis}",
                    "note": "GEO therapy; A=on-treatment after 1 dose, B=pretreatment",
                    **sp,
                }
            )

    testdf = pd.DataFrame(tests)
    testdf.to_csv(TABLES / "stats.tsv", sep="\t", index=False)

    primary = testdf[(testdf["test"] == "spearman") & (testdf["subset"] == "all")].copy()

    def grab(pair: str, col: str) -> float:
        return float(primary.loc[primary["pair"] == pair, col].iloc[0])

    one = {
        "dataset": "GSE121810",
        "n": 29,
        "n_CLDN4_CPM_ge_1": int((cldn4_cpm >= 1).sum()),
        "CLDN4_CPM_median": float(cldn4_cpm.median()),
        "matrix": "supplementary HUGO counts (series matrix has 0 expression rows)",
        "CLDN4_vs_CD274_rho": grab("CLDN4 vs CD274", "rho"),
        "CLDN4_vs_CD274_p": grab("CLDN4 vs CD274", "p"),
        "CLDN4_vs_TNK_rho": grab("CLDN4 vs TNK", "rho"),
        "CLDN4_vs_TNK_p": grab("CLDN4 vs TNK", "p"),
        "CLDN4_vs_IFN_rho": grab("CLDN4 vs IFN", "rho"),
        "CLDN4_vs_IFN_p": grab("CLDN4 vs IFN", "p"),
        "CLDN4_vs_CD8A_rho": grab("CLDN4 vs CD8A", "rho"),
        "CLDN4_vs_CD8A_p": grab("CLDN4 vs CD8A", "p"),
        "TNK_genes": f"{len(tnk_used)}/{len(TNK_GENES)}",
        "IFN_genes": f"{len(ifn_used)}/{len(IFN_AYERS)}",
        "IFNG": "row present; 0/29 CPM≥1",
        "RECIST_n": 0,
        "OS_n_this_slice": 0,
    }
    pd.DataFrame([one]).to_csv(TABLES / "one_row.tsv", sep="\t", index=False)

    scatter(
        FIGURES / "fig1_cldn4_vs_cd274.png",
        per["CLDN4"].to_numpy(),
        per["CD274"].to_numpy(),
        "CLDN4 log2(CPM+1)",
        "CD274 log2(CPM+1)",
        "GSE121810 CLDN4 vs CD274",
        one["CLDN4_vs_CD274_rho"],
        one["CLDN4_vs_CD274_p"],
        29,
    )
    scatter(
        FIGURES / "fig2_cldn4_vs_tnk.png",
        per["CLDN4"].to_numpy(),
        per["TNK"].to_numpy(),
        "CLDN4 log2(CPM+1)",
        f"T/NK mean-z ({len(tnk_used)}/{len(TNK_GENES)})",
        "GSE121810 CLDN4 vs T/NK",
        one["CLDN4_vs_TNK_rho"],
        one["CLDN4_vs_TNK_p"],
        29,
    )
    scatter(
        FIGURES / "fig3_cldn4_vs_cd8a.png",
        per["CLDN4"].to_numpy(),
        per["CD8A"].to_numpy(),
        "CLDN4 log2(CPM+1)",
        "CD8A log2(CPM+1)",
        "GSE121810 CLDN4 vs CD8A",
        one["CLDN4_vs_CD8A_rho"],
        one["CLDN4_vs_CD8A_p"],
        29,
    )
    forest_rows = [
        {"pair": r["pair"].replace("CLDN4 vs ", ""), "n": int(r["n"]), "rho": float(r["rho"]), "p": float(r["p"])}
        for _, r in primary.iterrows()
        if r["pair"].startswith("CLDN4 vs ")
        and r["pair"]
        in {
            "CLDN4 vs CD274",
            "CLDN4 vs TNK",
            "CLDN4 vs IFN",
            "CLDN4 vs CD8A",
            "CLDN4 vs NKG7",
            "CLDN4 vs PTPRC",
            "CLDN4 vs MHC",
        }
    ]
    forest(FIGURES / "fig4_spearman_forest.png", forest_rows)
    detect_bar(
        FIGURES / "fig5_detectability.png",
        det,
        ["CLDN4", "TACSTD2", "EPCAM", "IFNG", "CD274", "CD8A", "PTPRC", "GFAP"],
    )

    summary = {
        "accession": "GSE121810",
        "pmid": "30742122",
        "n_samples": 29,
        "n_genes": int(counts.shape[0]),
        "series_matrix_expression_rows": n_expr_rows,
        "arm": arm_counts,
        "tnk_used": tnk_used,
        "ifn_used": ifn_used,
        "mhc_used": mhc_used,
        "cldn4_detectability": detectability(counts, "CLDN4"),
        "cd274_detectability": detectability(counts, "CD274"),
        "ifng_detectability": detectability(counts, "IFNG"),
        "one_row": one,
        "series_title": series.get("title", ""),
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(one, indent=2))
    print("wrote", TABLES, FIGURES)


if __name__ == "__main__":
    main()
