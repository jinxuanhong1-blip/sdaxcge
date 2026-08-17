#!/usr/bin/env python3
"""GSE248378 leftover: CLDN4 vs T/NK or IFN on the public post-durva FPKM.

The GEO series matrix is public metadata only (0 expression rows). Tests use
the deposited supplementary FPKM (GSE248378_Durva_Post_FPKMs.txt.gz).

Recurrence / MPR / PFS are not GEO characteristics and are not recovered here
(prior leftover PRs 175 / 224 already scored CLDN4 vs recurrence). This slice
is CLDN4 vs T/NK and IFN only.

Downloads stay under $GSE248378_CLDN4_DATA (default /tmp/gse248378) and are
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
DATA = Path(os.environ.get("GSE248378_CLDN4_DATA", "/tmp/gse248378"))

MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE248nnn/GSE248378/"
    "matrix/GSE248378_series_matrix.txt.gz"
)
FPKM_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE248nnn/GSE248378/"
    "suppl/GSE248378_Durva_Post_FPKMs.txt.gz"
)

# Locked bulk T/NK list: CD8 / CD3 plus cytotoxic NK genes present on this
# Twist-exome FPKM. KLRF1 is absent and is not added after the fact.
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

# Ayers IFN-γ 6-gene. IFNG is absent from the deposited FPKM (honest 5/6).
IFN_AYERS = ["IFNG", "STAT1", "CXCL9", "CXCL10", "IDO1", "HLA-DRA"]
MHC_I = ["HLA-A", "HLA-B", "HLA-C"]
SINGLE = ["CLDN4", "TACSTD2", "CD8A", "NKG7", "IFNG", "STAT1", "CXCL9", "CXCL10", "EPCAM", "MKI67"]


def dl(dest: Path, url: str) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    print(f"[get] {url}", flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": "gse248378-cldn4/1.0"})
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
                    # repeated characteristics_ch1 / data_processing lines
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
    meta["histology"] = chars.str.extract(r"cell type:\s*([^|]+)", expand=False).str.strip()
    meta["arm"] = chars.str.extract(r"treatment:\s*([^|]+)", expand=False).str.strip()
    meta["title"] = meta["title"].astype(str)
    return meta, {k: " | ".join(v) for k, v in series.items()}, n_expr_rows


def mean_z(log: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, list[str]]:
    found = [g for g in genes if g in log.index]
    if not found:
        return pd.Series(np.nan, index=log.columns), []
    z = log.loc[found].astype(float).apply(lambda r: (r - r.mean()) / (r.std(ddof=0) or np.nan), axis=1)
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


def fmt(x, nd=3):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "—"
    return f"{x:.{nd}f}"


def scatter(path, x, y, xlab, ylab, title, rho, p, n):
    fig, ax = plt.subplots(figsize=(4.4, 3.8))
    ax.scatter(x, y, s=28, c="#1f4e79", alpha=0.85, edgecolors="none")
    if np.isfinite(x).sum() >= 3:
        coef = np.polyfit(x, y, 1)
        xs = np.linspace(np.nanmin(x), np.nanmax(x), 50)
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
    fig, ax = plt.subplots(figsize=(6.2, 3.6))
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
    ax.set_title("GSE248378 CLDN4 vs T/NK / IFN (n=29)", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)

    matrix_path = dl(DATA / "GSE248378_series_matrix.txt.gz", MATRIX_URL)
    fpkm_path = dl(DATA / "GSE248378_Durva_Post_FPKMs.txt.gz", FPKM_URL)

    meta, series, n_expr_rows = parse_series_meta(matrix_path)
    fpkm = pd.read_csv(fpkm_path, sep="\t", index_col=0)
    fpkm.columns = fpkm.columns.astype(str)

    # Join FPKM titles to GEO sample titles (one column per GSM).
    title_to_gsm = {t: g for g, t in zip(meta.index, meta["title"])}
    missing_titles = [c for c in fpkm.columns if c not in title_to_gsm]
    extra_titles = [t for t in meta["title"] if t not in set(fpkm.columns)]
    if missing_titles or extra_titles:
        raise SystemExit(f"title mismatch FPKM-only={missing_titles} GEO-only={extra_titles}")

    fpkm = fpkm.loc[:, meta["title"].tolist()]
    fpkm.columns = meta.index
    log = np.log2(fpkm.astype(float) + 1.0)

    tnk, tnk_used = mean_z(log, TNK_GENES)
    ifn, ifn_used = mean_z(log, IFN_AYERS)
    mhc, mhc_used = mean_z(log, MHC_I)

    per = meta[["title", "histology", "arm", "tissue"]].copy()
    per["CLDN4"] = log.loc["CLDN4"].to_numpy()
    per["TACSTD2"] = log.loc["TACSTD2"].to_numpy()
    per["CD8A"] = log.loc["CD8A"].to_numpy()
    per["NKG7"] = log.loc["NKG7"].to_numpy()
    per["TNK"] = tnk.to_numpy()
    per["IFN"] = ifn.to_numpy()
    per["MHC"] = mhc.to_numpy()
    per["EPCAM"] = log.loc["EPCAM"].to_numpy()
    per["MKI67"] = log.loc["MKI67"].to_numpy()
    per["CLDN4_FPKM"] = fpkm.loc["CLDN4"].to_numpy()
    per.to_csv(TABLES / "per_sample.tsv", sep="\t")

    hist_counts = per["histology"].value_counts().to_dict()
    arm_counts = per["arm"].value_counts().to_dict()

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
    cov = pd.DataFrame(coverage_rows)
    cov.to_csv(TABLES / "gene_coverage.tsv", sep="\t", index=False)

    label_inv = pd.DataFrame(
        [
            {"field": "series_matrix_expression_rows", "public": "yes (empty table)", "n": n_expr_rows, "note": "RNA-seq series matrix; 0 features"},
            {"field": "supplementary_FPKM_samples", "public": "yes", "n": int(fpkm.shape[1]), "note": "GSE248378_Durva_Post_FPKMs.txt.gz"},
            {"field": "supplementary_FPKM_genes", "public": "yes", "n": int(fpkm.shape[0]), "note": "HGNC symbols; Cufflinks FPKM, Gencode v19"},
            {"field": "GEO_GSM", "public": "yes", "n": int(len(meta)), "note": "GSM7912309–GSM7912337"},
            {"field": "treatment_arm", "public": "yes", "n": int(len(per)), "note": f"Arm1={arm_counts.get('Arm1', 0)}; Arm2={arm_counts.get('Arm2', 0)}"},
            {"field": "histology", "public": "yes", "n": int(len(per)), "note": "; ".join(f"{k}={v}" for k, v in sorted(hist_counts.items()))},
            {"field": "MPR", "public": "no", "n": 0, "note": "not a GEO characteristic"},
            {"field": "recurrence", "public": "no", "n": 0, "note": "not a GEO characteristic; not recovered here"},
            {"field": "PFS_OS", "public": "no", "n": 0, "note": "not a GEO characteristic"},
            {"field": "CLDN4", "public": "yes", "n": int(len(per)), "note": "present on FPKM"},
            {"field": "IFNG", "public": "no", "n": 0, "note": "absent from deposited FPKM; IFNAR/IFNGR present"},
        ]
    )
    label_inv.to_csv(TABLES / "label_inventory.tsv", sep="\t", index=False)

    tests = []

    def add_spearman(pair, x, y, note=""):
        sp = spearman(x, y)
        tests.append({"test": "spearman", "pair": pair, "note": note, **sp})

    add_spearman("CLDN4 vs TNK", per["CLDN4"], per["TNK"], f"mean-z {len(tnk_used)}/{len(TNK_GENES)}")
    add_spearman("CLDN4 vs IFN", per["CLDN4"], per["IFN"], f"Ayers mean-z {len(ifn_used)}/{len(IFN_AYERS)}; IFNG absent")
    add_spearman("CLDN4 vs CD8A", per["CLDN4"], per["CD8A"])
    add_spearman("CLDN4 vs NKG7", per["CLDN4"], per["NKG7"])
    add_spearman("CLDN4 vs MHC", per["CLDN4"], per["MHC"], f"HLA-A/B/C {len(mhc_used)}/3")
    add_spearman("CLDN4 vs TACSTD2", per["CLDN4"], per["TACSTD2"], "companion")
    add_spearman("CLDN4 vs EPCAM", per["CLDN4"], per["EPCAM"], "epithelial content")
    add_spearman("TACSTD2 vs TNK", per["TACSTD2"], per["TNK"], "companion")
    add_spearman("TACSTD2 vs IFN", per["TACSTD2"], per["IFN"], "companion")
    add_spearman("CD8A vs TNK", per["CD8A"], per["TNK"], "positive-control; CD8A is in TNK")
    add_spearman("TNK vs IFN", per["TNK"], per["IFN"], "immune-axis sanity")

    for axis in ["TNK", "IFN", "CD8A", "NKG7", "MHC"]:
        psp = partial_spearman(per["CLDN4"], per[axis], per["EPCAM"])
        tests.append({"test": "partial_spearman_|EPCAM", "pair": f"CLDN4 vs {axis}", "note": "rank residual on EPCAM", **psp})

    med = per["CLDN4"].median()
    hi = per["CLDN4"] >= med
    lo = ~hi
    for axis in ["TNK", "IFN", "CD8A", "NKG7", "MHC"]:
        blk = mwu(per.loc[hi, axis], per.loc[lo, axis])
        tests.append(
            {
                "test": "median_split_MWU",
                "pair": f"CLDN4-high vs low, {axis}",
                "note": f"high n={int(hi.sum())} (>= median {med:.3f}); low n={int(lo.sum())}",
                "n": int(hi.sum() + lo.sum()),
                "rho": blk["delta"],
                "p": blk["p"],
                **{f"mwu_{k}": v for k, v in blk.items()},
            }
        )

    # Histology split only where n is honest and not a singleton.
    adeno = per["histology"] == "Adenocarcinoma"
    squam = per["histology"] == "Squamous"
    for name, mask in [("Adenocarcinoma", adeno), ("Squamous", squam)]:
        sub = per.loc[mask]
        for axis in ["TNK", "IFN", "CD8A"]:
            sp = spearman(sub["CLDN4"], sub[axis])
            tests.append(
                {
                    "test": "spearman_histology",
                    "pair": f"CLDN4 vs {axis} | {name}",
                    "note": f"GEO cell type; n={int(mask.sum())}",
                    **sp,
                }
            )

    for name, mask in [("Arm1", per["arm"] == "Arm1"), ("Arm2", per["arm"] == "Arm2")]:
        sub = per.loc[mask]
        for axis in ["TNK", "IFN"]:
            sp = spearman(sub["CLDN4"], sub[axis])
            tests.append(
                {
                    "test": "spearman_arm",
                    "pair": f"CLDN4 vs {axis} | {name}",
                    "note": "GEO treatment; Arm1=durva, Arm2=durva+SBRT (paper)",
                    **sp,
                }
            )

    testdf = pd.DataFrame(tests)
    testdf.to_csv(TABLES / "stats.tsv", sep="\t", index=False)

    primary = testdf[testdf["test"] == "spearman"].copy()
    one = {
        "dataset": "GSE248378",
        "n": 29,
        "matrix": "supplementary FPKM (series matrix has 0 expression rows)",
        "CLDN4_vs_TNK_rho": float(primary.loc[primary["pair"] == "CLDN4 vs TNK", "rho"].iloc[0]),
        "CLDN4_vs_TNK_p": float(primary.loc[primary["pair"] == "CLDN4 vs TNK", "p"].iloc[0]),
        "CLDN4_vs_IFN_rho": float(primary.loc[primary["pair"] == "CLDN4 vs IFN", "rho"].iloc[0]),
        "CLDN4_vs_IFN_p": float(primary.loc[primary["pair"] == "CLDN4 vs IFN", "p"].iloc[0]),
        "CLDN4_vs_CD8A_rho": float(primary.loc[primary["pair"] == "CLDN4 vs CD8A", "rho"].iloc[0]),
        "CLDN4_vs_CD8A_p": float(primary.loc[primary["pair"] == "CLDN4 vs CD8A", "p"].iloc[0]),
        "TNK_genes": f"{len(tnk_used)}/{len(TNK_GENES)}",
        "IFN_genes": f"{len(ifn_used)}/{len(IFN_AYERS)}",
        "IFNG": "absent",
        "MPR_n": 0,
        "recurrence_n": 0,
    }
    pd.DataFrame([one]).to_csv(TABLES / "one_row.tsv", sep="\t", index=False)

    scatter(
        FIGURES / "fig1_cldn4_vs_tnk.png",
        per["CLDN4"].to_numpy(),
        per["TNK"].to_numpy(),
        "CLDN4 log2(FPKM+1)",
        f"T/NK mean-z ({len(tnk_used)}/{len(TNK_GENES)})",
        "GSE248378 CLDN4 vs T/NK",
        one["CLDN4_vs_TNK_rho"],
        one["CLDN4_vs_TNK_p"],
        29,
    )
    scatter(
        FIGURES / "fig2_cldn4_vs_ifn.png",
        per["CLDN4"].to_numpy(),
        per["IFN"].to_numpy(),
        "CLDN4 log2(FPKM+1)",
        f"IFN Ayers mean-z ({len(ifn_used)}/{len(IFN_AYERS)}; no IFNG)",
        "GSE248378 CLDN4 vs IFN",
        one["CLDN4_vs_IFN_rho"],
        one["CLDN4_vs_IFN_p"],
        29,
    )
    scatter(
        FIGURES / "fig3_cldn4_vs_cd8a.png",
        per["CLDN4"].to_numpy(),
        per["CD8A"].to_numpy(),
        "CLDN4 log2(FPKM+1)",
        "CD8A log2(FPKM+1)",
        "GSE248378 CLDN4 vs CD8A",
        one["CLDN4_vs_CD8A_rho"],
        one["CLDN4_vs_CD8A_p"],
        29,
    )
    forest_rows = [
        {"pair": r["pair"].replace("CLDN4 vs ", ""), "n": int(r["n"]), "rho": float(r["rho"]), "p": float(r["p"])}
        for _, r in primary.iterrows()
        if r["pair"].startswith("CLDN4 vs ")
    ]
    forest(FIGURES / "fig4_spearman_forest.png", forest_rows)

    summary = {
        "accession": "GSE248378",
        "pmid": ["38114518", "38401548"],
        "n_samples": 29,
        "n_genes_fpkm": int(fpkm.shape[0]),
        "series_matrix_expression_rows": n_expr_rows,
        "histology": hist_counts,
        "arm": arm_counts,
        "tnk_used": tnk_used,
        "ifn_used": ifn_used,
        "mhc_used": mhc_used,
        "one_row": one,
        "series_title": series.get("title", ""),
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(one, indent=2))
    print("wrote", TABLES, FIGURES)


if __name__ == "__main__":
    main()
