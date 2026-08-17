#!/usr/bin/env python3
"""GSE179351 leftover: CLDN4 vs immune / CD274 on public DESeq2 counts.

Parikh et al., Nat Cancer 2021 (PMID 35122060). Phase II radiation + ipilimumab
+ nivolumab in metastatic MSS CRC and PDAC. GEO GSE179351 is public
(2021-07-06) with series-matrix metadata and supplementary DESeq2-normalized
and raw counts.

RECIST / DCR / ORR are not GEO characteristics (n=0) and are not recovered
from the paper. Multiple timepoints exist per patient; the primary unit is
one baseline biopsy per patient (Pre-Tx, else Pre-xRT).

Downloads stay under $GSE179351_CLDN4_DATA (default /tmp/gse179351) and are
not committed.
"""

from __future__ import annotations

import gzip
import json
import os
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
DATA = Path(os.environ.get("GSE179351_CLDN4_DATA", "/tmp/gse179351"))

MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179351/"
    "matrix/GSE179351_series_matrix.txt.gz"
)
NORM_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179351/"
    "suppl/GSE179351_DESeq2_NormalizedCountsForAllSamples.txt.gz"
)

ANN_COLS = [
    "gene.id",
    "ensembl.gene.id",
    "chr",
    "gene.symbol",
    "gene.type",
    "length",
    "description",
    "entrez.gene.id",
    "uniprot.id",
    "go.gene.id",
]

# Locked bulk T/NK list used on leftover lung ICI FPKM pages.
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
IFN_AYERS = ["IFNG", "STAT1", "CXCL9", "CXCL10", "IDO1", "HLA-DRA"]
MHC_I = ["HLA-A", "HLA-B", "HLA-C"]
SINGLE = [
    "CLDN4",
    "TACSTD2",
    "CD274",
    "CD8A",
    "NKG7",
    "IFNG",
    "PDCD1",
    "CTLA4",
    "EPCAM",
    "MKI67",
]


def dl(dest: Path, url: str) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    print(f"[get] {url}", flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": "gse179351-cldn4/1.0"})
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
                        a + " | " + b if a else b
                        for a, b in zip(sample_fields[key], vals)
                    ]
                else:
                    sample_fields[key] = vals
    meta = pd.DataFrame(sample_fields)
    meta.index = meta["geo_accession"].astype(str)
    meta.index.name = "gsm"
    chars = meta["characteristics_ch1"].astype(str)
    meta["tissue"] = chars.str.extract(r"tissue:\s*([^|]+)", expand=False).str.strip()
    meta["cancer"] = chars.str.extract(r"cancer type:\s*([^|]+)", expand=False).str.strip()
    meta["timepoint"] = chars.str.extract(r"time point:\s*([^|]+)", expand=False).str.strip()
    meta["title"] = meta["title"].astype(str)
    meta["patient"] = meta["title"].str.extract(r"^(\d+)", expand=False)
    return meta, {k: " | ".join(v) for k, v in series.items()}, n_expr_rows


def count_title(col: str) -> str:
    s = col[1:] if col.startswith("X") else col
    return s.replace(".", "-")


def pick_symbol_matrix(norm: pd.DataFrame) -> pd.DataFrame:
    """One protein-coding row per HGNC symbol; skip empty symbols."""
    df = norm.copy()
    df["gene.symbol"] = df["gene.symbol"].astype(str)
    df = df[df["gene.symbol"].notna() & (df["gene.symbol"] != "") & (df["gene.symbol"] != "nan")]
    pc = df[df["gene.type"].astype(str) == "protein_coding"]
    # If a symbol has several protein-coding rows, keep the highest mean.
    sample_cols = [c for c in df.columns if c not in ANN_COLS]
    pc = pc.copy()
    pc["_mean"] = pc[sample_cols].astype(float).mean(axis=1)
    pc = pc.sort_values("_mean", ascending=False).drop_duplicates("gene.symbol", keep="first")
    mat = pc.set_index("gene.symbol")[sample_cols].astype(float)
    return mat


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


def pick_baseline(g: pd.DataFrame) -> pd.Series:
    order = {"Pre-Tx": 0, "Pre-xRT": 1, "Post-xRT": 2}
    g = g.copy()
    g["_ord"] = g["timepoint"].map(order).fillna(9)
    g = g.sort_values(["_ord", "title"])
    return g.iloc[0]


def scatter(path, x, y, xlab, ylab, title, rho, p, n):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
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


def forest(path, rows, title):
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
    ax.set_title(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def add_block(tests: list, per: pd.DataFrame, cohort: str, note: str) -> None:
    axes = [
        ("CD274", "CLDN4 vs CD274"),
        ("TNK", "CLDN4 vs TNK"),
        ("IFN", "CLDN4 vs IFN"),
        ("CD8A", "CLDN4 vs CD8A"),
        ("NKG7", "CLDN4 vs NKG7"),
        ("MHC", "CLDN4 vs MHC"),
        ("TACSTD2", "CLDN4 vs TACSTD2"),
        ("EPCAM", "CLDN4 vs EPCAM"),
        ("PDCD1", "CLDN4 vs PDCD1"),
    ]
    for col, pair in axes:
        if col not in per.columns:
            continue
        sp = spearman(per["CLDN4"], per[col])
        tests.append({"test": "spearman", "cohort": cohort, "pair": pair, "note": note, **sp})

    if "EPCAM" in per.columns:
        for col, pair in [
            ("CD274", "CLDN4 vs CD274"),
            ("TNK", "CLDN4 vs TNK"),
            ("IFN", "CLDN4 vs IFN"),
            ("CD8A", "CLDN4 vs CD8A"),
            ("NKG7", "CLDN4 vs NKG7"),
            ("MHC", "CLDN4 vs MHC"),
        ]:
            psp = partial_spearman(per["CLDN4"], per[col], per["EPCAM"])
            tests.append(
                {
                    "test": "partial_spearman_|EPCAM",
                    "cohort": cohort,
                    "pair": pair,
                    "note": note + "; rank residual on EPCAM",
                    **psp,
                }
            )

    med = per["CLDN4"].median()
    hi = per["CLDN4"] >= med
    lo = ~hi
    for col in ["CD274", "TNK", "IFN", "CD8A", "NKG7", "MHC"]:
        if col not in per.columns:
            continue
        blk = mwu(per.loc[hi, col], per.loc[lo, col])
        tests.append(
            {
                "test": "median_split_MWU",
                "cohort": cohort,
                "pair": f"CLDN4-high vs low, {col}",
                "note": f"{note}; high n={int(hi.sum())} (>= median {med:.3f}); low n={int(lo.sum())}",
                "n": int(hi.sum() + lo.sum()),
                "rho": blk["delta"],
                "p": blk["p"],
                **{f"mwu_{k}": v for k, v in blk.items()},
            }
        )


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)

    matrix_path = dl(DATA / "GSE179351_series_matrix.txt.gz", MATRIX_URL)
    norm_path = dl(DATA / "GSE179351_DESeq2_NormalizedCountsForAllSamples.txt.gz", NORM_URL)

    meta, series, n_expr_rows = parse_series_meta(matrix_path)
    norm = pd.read_csv(norm_path, sep="\t", low_memory=False)
    mat = pick_symbol_matrix(norm)
    mat.columns = [count_title(c) for c in mat.columns]

    series_titles = set(meta["title"])
    count_titles = set(mat.columns)
    geo_only = sorted(series_titles - count_titles)
    count_only = sorted(count_titles - series_titles)
    joined_titles = [t for t in meta["title"] if t in count_titles]

    log = np.log2(mat.astype(float) + 1.0)
    tnk, tnk_used = mean_z(log, TNK_GENES)
    ifn, ifn_used = mean_z(log, IFN_AYERS)
    mhc, mhc_used = mean_z(log, MHC_I)

    per = meta.loc[meta["title"].isin(joined_titles)].copy()
    per = per.set_index("title", drop=False)
    per.index.name = "sample_title"
    for gene in SINGLE:
        per[gene] = log.loc[gene, per.index].to_numpy() if gene in log.index else np.nan
    per["TNK"] = tnk.loc[per.index].to_numpy()
    per["IFN"] = ifn.loc[per.index].to_numpy()
    per["MHC"] = mhc.loc[per.index].to_numpy()
    per["CLDN4_norm"] = mat.loc["CLDN4", per.index].to_numpy()
    per["CD274_norm"] = mat.loc["CD274", per.index].to_numpy()

    # Patient 52 is MSS CRC on Pre-Tx and PDAC on later GEO rows.
    conflict = (
        per.groupby("patient")["cancer"]
        .nunique()
        .loc[lambda s: s > 1]
        .index.astype(str)
        .tolist()
    )

    baseline = pd.DataFrame(
        [pick_baseline(g) for _, g in per.groupby("patient", sort=True)]
    ).reset_index(drop=True)
    num_cols = [
        c
        for c in [
            "CLDN4",
            "TACSTD2",
            "CD274",
            "CD8A",
            "NKG7",
            "IFNG",
            "PDCD1",
            "CTLA4",
            "EPCAM",
            "MKI67",
            "TNK",
            "IFN",
            "MHC",
            "CLDN4_norm",
            "CD274_norm",
        ]
        if c in baseline.columns
    ]
    baseline[num_cols] = baseline[num_cols].apply(pd.to_numeric, errors="coerce")

    pretx = per[per["timepoint"] == "Pre-Tx"].copy()
    # One Pre-Tx row per patient (should already be unique).
    pretx = pretx.drop_duplicates("patient")

    keep = [
        "geo_accession",
        "title",
        "patient",
        "cancer",
        "timepoint",
        "tissue",
        "CLDN4",
        "CD274",
        "TACSTD2",
        "CD8A",
        "NKG7",
        "IFNG",
        "PDCD1",
        "CTLA4",
        "EPCAM",
        "TNK",
        "IFN",
        "MHC",
    ]
    per[keep].to_csv(TABLES / "per_sample.tsv", sep="\t")
    baseline[keep].to_csv(TABLES / "per_patient_baseline.tsv", sep="\t", index=False)

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

    cancer_pt = baseline["cancer"].value_counts().to_dict()
    time_all = per["timepoint"].value_counts().to_dict()
    cancer_all = per["cancer"].value_counts().to_dict()

    label_inv = pd.DataFrame(
        [
            {
                "field": "GEO_status",
                "public": "yes",
                "n": 54,
                "note": "Public on Jul 06 2021; PMID 35122060",
            },
            {
                "field": "series_matrix_GSM",
                "public": "yes",
                "n": int(len(meta)),
                "note": "GSM5416019–GSM5416072",
            },
            {
                "field": "series_matrix_expression_rows",
                "public": "yes (empty table)",
                "n": n_expr_rows,
                "note": "RNA-seq series matrix; expression is in supplementary counts",
            },
            {
                "field": "DESeq2_normalized_columns",
                "public": "yes",
                "n": int(mat.shape[1]),
                "note": "GSE179351_DESeq2_NormalizedCountsForAllSamples.txt.gz",
            },
            {
                "field": "joined_title_samples",
                "public": "yes",
                "n": int(len(per)),
                "note": f"GEO-only titles={geo_only}; count-only titles={count_only}",
            },
            {
                "field": "unique_patients_joined",
                "public": "yes",
                "n": int(per["patient"].nunique()),
                "note": "patient ID parsed from GEO title; 20-Pre-xRT does not join",
            },
            {
                "field": "baseline_one_per_patient",
                "public": "yes",
                "n": int(len(baseline)),
                "note": "Pre-Tx preferred, else Pre-xRT; primary n",
            },
            {
                "field": "Pre-Tx_only",
                "public": "yes",
                "n": int(len(pretx)),
                "note": "sensitivity; excludes Pre-xRT-only patients 7 and 51",
            },
            {
                "field": "cancer_type",
                "public": "yes",
                "n": int(len(per)),
                "note": f"samples {cancer_all}; baseline patients {cancer_pt}; patient 52 cancer label conflicts: {conflict}",
            },
            {
                "field": "timepoint",
                "public": "yes",
                "n": int(len(per)),
                "note": "; ".join(f"{k}={v}" for k, v in sorted(time_all.items())),
            },
            {
                "field": "RECIST_DCR_ORR",
                "public": "no",
                "n": 0,
                "note": "not a GEO characteristic; not recovered here",
            },
            {
                "field": "PFS_OS",
                "public": "no",
                "n": 0,
                "note": "not a GEO characteristic",
            },
            {
                "field": "CLDN4",
                "public": "yes",
                "n": int(len(per)),
                "note": "protein_coding ENSG00000189143",
            },
            {
                "field": "CD274",
                "public": "yes",
                "n": int(len(per)),
                "note": "protein_coding ENSG00000120217",
            },
            {
                "field": "IFNG",
                "public": "yes",
                "n": int(len(per)),
                "note": "present on DESeq2 matrix (Ayers 6/6)",
            },
        ]
    )
    label_inv.to_csv(TABLES / "label_inventory.tsv", sep="\t", index=False)

    tests: list[dict] = []
    add_block(
        tests,
        baseline,
        "baseline_patient",
        f"one biopsy/patient; n={len(baseline)}; Pre-Tx else Pre-xRT",
    )
    add_block(
        tests,
        pretx,
        "pre_tx",
        f"Pre-Tx only; n={len(pretx)}; still one/patient",
    )
    add_block(
        tests,
        per,
        "all_joined_samples",
        f"all joined biopsies n={len(per)}; pseudoreplication (not primary)",
    )

    for name, mask in [
        ("MSS CRC", baseline["cancer"] == "MSS CRC"),
        ("PDAC", baseline["cancer"] == "PDAC"),
    ]:
        sub = baseline.loc[mask]
        for col, pair in [
            ("CD274", "CLDN4 vs CD274"),
            ("TNK", "CLDN4 vs TNK"),
            ("IFN", "CLDN4 vs IFN"),
            ("CD8A", "CLDN4 vs CD8A"),
        ]:
            sp = spearman(sub["CLDN4"], sub[col])
            tests.append(
                {
                    "test": "spearman_cancer",
                    "cohort": f"baseline_{name.replace(' ', '_')}",
                    "pair": pair,
                    "note": f"GEO cancer type on baseline row; n={int(mask.sum())}",
                    **sp,
                }
            )

    testdf = pd.DataFrame(tests)
    testdf.to_csv(TABLES / "stats.tsv", sep="\t", index=False)

    def grab(cohort: str, pair: str, test: str = "spearman") -> dict:
        hit = testdf[(testdf["test"] == test) & (testdf["cohort"] == cohort) & (testdf["pair"] == pair)]
        if hit.empty:
            return {"n": 0, "rho": np.nan, "p": np.nan}
        r = hit.iloc[0]
        return {"n": int(r["n"]), "rho": float(r["rho"]) if pd.notna(r["rho"]) else np.nan, "p": float(r["p"]) if pd.notna(r["p"]) else np.nan}

    cd274 = grab("baseline_patient", "CLDN4 vs CD274")
    tnk_s = grab("baseline_patient", "CLDN4 vs TNK")
    ifn_s = grab("baseline_patient", "CLDN4 vs IFN")
    cd8_s = grab("baseline_patient", "CLDN4 vs CD8A")
    one = {
        "dataset": "GSE179351",
        "n_primary": int(len(baseline)),
        "n_joined_samples": int(len(per)),
        "n_geo_gsm": int(len(meta)),
        "n_pre_tx": int(len(pretx)),
        "n_response_geo": 0,
        "matrix": "suppl. DESeq2 normalized counts (series matrix empty)",
        "title_mismatch": f"GEO-only={geo_only}; count-only={count_only}",
        "CLDN4_vs_CD274_rho": cd274["rho"],
        "CLDN4_vs_CD274_p": cd274["p"],
        "CLDN4_vs_TNK_rho": tnk_s["rho"],
        "CLDN4_vs_TNK_p": tnk_s["p"],
        "CLDN4_vs_IFN_rho": ifn_s["rho"],
        "CLDN4_vs_IFN_p": ifn_s["p"],
        "CLDN4_vs_CD8A_rho": cd8_s["rho"],
        "CLDN4_vs_CD8A_p": cd8_s["p"],
        "TNK_genes": f"{len(tnk_used)}/{len(TNK_GENES)}",
        "IFN_genes": f"{len(ifn_used)}/{len(IFN_AYERS)}",
        "IFNG": "present",
        "baseline_CRC": int((baseline["cancer"] == "MSS CRC").sum()),
        "baseline_PDAC": int((baseline["cancer"] == "PDAC").sum()),
        "cancer_label_conflicts": ",".join(conflict) if conflict else "",
    }
    pd.DataFrame([one]).to_csv(TABLES / "one_row.tsv", sep="\t", index=False)

    xlab = "CLDN4 log2(DESeq2+1)"
    scatter(
        FIGURES / "fig1_cldn4_vs_cd274.png",
        baseline["CLDN4"].to_numpy(),
        baseline["CD274"].to_numpy(),
        xlab,
        "CD274 log2(DESeq2+1)",
        "GSE179351 baseline CLDN4 vs CD274",
        one["CLDN4_vs_CD274_rho"],
        one["CLDN4_vs_CD274_p"],
        one["n_primary"],
    )
    scatter(
        FIGURES / "fig2_cldn4_vs_tnk.png",
        baseline["CLDN4"].to_numpy(),
        baseline["TNK"].to_numpy(),
        xlab,
        f"T/NK mean-z ({len(tnk_used)}/{len(TNK_GENES)})",
        "GSE179351 baseline CLDN4 vs T/NK",
        one["CLDN4_vs_TNK_rho"],
        one["CLDN4_vs_TNK_p"],
        one["n_primary"],
    )
    scatter(
        FIGURES / "fig3_cldn4_vs_ifn.png",
        baseline["CLDN4"].to_numpy(),
        baseline["IFN"].to_numpy(),
        xlab,
        f"IFN Ayers mean-z ({len(ifn_used)}/{len(IFN_AYERS)})",
        "GSE179351 baseline CLDN4 vs IFN",
        one["CLDN4_vs_IFN_rho"],
        one["CLDN4_vs_IFN_p"],
        one["n_primary"],
    )
    forest_rows = []
    for pair in [
        "CLDN4 vs CD274",
        "CLDN4 vs TNK",
        "CLDN4 vs IFN",
        "CLDN4 vs CD8A",
        "CLDN4 vs NKG7",
        "CLDN4 vs MHC",
        "CLDN4 vs TACSTD2",
        "CLDN4 vs EPCAM",
    ]:
        r = grab("baseline_patient", pair)
        forest_rows.append(
            {"pair": pair.replace("CLDN4 vs ", ""), "n": r["n"], "rho": r["rho"], "p": r["p"]}
        )
    forest(
        FIGURES / "fig4_spearman_forest.png",
        forest_rows,
        f"GSE179351 CLDN4 vs immune/CD274 (baseline n={one['n_primary']})",
    )

    summary = {
        "accession": "GSE179351",
        "pmid": "35122060",
        "public": True,
        "n_geo_gsm": int(len(meta)),
        "n_count_columns": int(mat.shape[1]),
        "n_joined_samples": int(len(per)),
        "n_patients_joined": int(per["patient"].nunique()),
        "n_baseline": int(len(baseline)),
        "n_pre_tx": int(len(pretx)),
        "n_response_geo": 0,
        "geo_only_titles": geo_only,
        "count_only_titles": count_only,
        "cancer_label_conflicts": conflict,
        "baseline_cancer": cancer_pt,
        "timepoints_joined": time_all,
        "tnk_used": tnk_used,
        "ifn_used": ifn_used,
        "mhc_used": mhc_used,
        "one_row": one,
        "series_title": series.get("title", ""),
        "series_status": series.get("status", ""),
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(one, indent=2))
    print("wrote", TABLES, FIGURES)


if __name__ == "__main__":
    main()
