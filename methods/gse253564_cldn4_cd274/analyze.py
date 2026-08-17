#!/usr/bin/env python3
"""ADDITIVE leftover — GSE253564 CLDN4 vs CD274 / HLA continuous.

Public pre-treatment FPKM only. No FASTQ / SRA. No GSEA (that leftover is
methods/gse253564_cldn4_gsea). Does not re-run MPR/PFS (opus_geo_leftover).

CLDN4 only. Patient/tumour is the unit. Honest pairwise-complete n.
Spearman vs CD274 and HLA-A/B/C; MHC-I mean-z is HLA-A/B/C only (3/3),
not B2M/TAP. Partial residualizes ranks on CD8A.
"""
from __future__ import annotations

import hashlib
import json
import math
import time
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
FIG = HERE / "figures"
TAB = HERE / "tables"
for d in (DATA, FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

ACCESSION = "GSE253564"
FPKM_NAME = "GSE253564_Pre-treatment_Samples_Pubs_FPKMs.txt.gz"
FPKM_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE253nnn/GSE253564/suppl/"
    "GSE253564_Pre-treatment_Samples_Pubs_FPKMs.txt.gz"
)

# Locked symbols. Ensembl used only as a fallback if the matrix is Ensembl-keyed.
GENES = {
    "CLDN4": ["CLDN4", "ENSG00000189143"],
    "CD274": ["CD274", "PDCD1LG1", "PDL1", "B7-H1", "ENSG00000120217"],
    "HLA-A": ["HLA-A", "HLAA", "ENSG00000206503"],
    "HLA-B": ["HLA-B", "HLAB", "ENSG00000234745"],
    "HLA-C": ["HLA-C", "HLAC", "ENSG00000204525"],
    "CD8A": ["CD8A", "ENSG00000153563"],
    "B2M": ["B2M", "ENSG00000166710"],
    "TACSTD2": ["TACSTD2", "TROP2", "ENSG00000184292"],
}
HLA_I = ["HLA-A", "HLA-B", "HLA-C"]
AXES = ["CD274", "HLA-A", "HLA-B", "HLA-C", "MHCI"]


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download_fpkm() -> Path:
    dest = DATA / FPKM_NAME
    if dest.exists() and dest.stat().st_size > 1000:
        log(f"using cached {dest.name} ({dest.stat().st_size} bytes)")
        return dest
    log(f"download {FPKM_URL}")
    req = urllib.request.Request(FPKM_URL, headers={"User-Agent": "gse253564-cldn4-cd274/1.0"})
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(req, timeout=600) as fh, tmp.open("wb") as out:
        while True:
            chunk = fh.read(1 << 20)
            if not chunk:
                break
            out.write(chunk)
    tmp.rename(dest)
    log(f"saved {dest.name} ({dest.stat().st_size} bytes) sha256={sha256(dest)}")
    return dest


def collapse_symbols(expr: pd.DataFrame) -> pd.DataFrame:
    expr = expr.copy()
    expr.index = expr.index.astype(str).str.replace(r"\.\d+$", "", regex=True)
    expr.index = expr.index.str.strip().str.strip('"')
    expr = expr.apply(pd.to_numeric, errors="coerce")
    if expr.index.duplicated().any():
        means = expr.mean(axis=1)
        keep_idx = means.groupby(level=0).idxmax()
        expr = expr.loc[keep_idx]
    return expr.loc[~expr.index.duplicated(keep="first")]


def load_expr(path: Path) -> pd.DataFrame:
    raw = pd.read_csv(path, sep="\t")
    gene_col = raw.columns[0]
    drop = [c for c in raw.columns if c.lower().replace(".", "_") in {"entrez_id", "entrez"}]
    expr = raw.set_index(gene_col).drop(columns=drop, errors="ignore")
    expr = collapse_symbols(expr)
    expr = np.log2(expr.clip(lower=0) + 1.0)
    return expr


def resolve_gene(index: pd.Index, aliases: list[str]) -> str | None:
    upper = {str(i).upper(): str(i) for i in index}
    for a in aliases:
        if a.upper() in upper:
            return upper[a.upper()]
    return None


def mean_z(expr: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, list[str]]:
    present = [g for g in genes if g in expr.index]
    if not present:
        return pd.Series(np.nan, index=expr.columns), []
    sub = expr.loc[present]
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1, ddof=1).replace(0, np.nan), axis=0)
    return z.mean(axis=0), present


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


def fmt_p(p) -> str:
    if p is None or (isinstance(p, float) and (p != p)):
        return "NA"
    if abs(p) < 1e-4:
        return f"{p:.2e}"
    return f"{p:.3f}"


def fmt_rho(x) -> str:
    if x is None or (isinstance(x, float) and (x != x)):
        return "NA"
    return f"{x:+.3f}"


def save_scatter(path: Path, x, y, title: str, xlab: str, ylab: str) -> None:
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = ~(np.isnan(x) | np.isnan(y))
    fig, ax = plt.subplots(figsize=(4.2, 4.0))
    ax.scatter(x[m], y[m], s=32, c="#2c3e50", edgecolors="white", linewidths=0.4, alpha=0.9)
    if m.sum() >= 3:
        lr = stats.linregress(x[m], y[m])
        xs = np.linspace(float(np.nanmin(x[m])), float(np.nanmax(x[m])), 50)
        ax.plot(xs, lr.intercept + lr.slope * xs, color="#c0392b", lw=1.2)
    sp = spearman(x, y)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(f"{title}\nρ={fmt_rho(sp['rho'])}  p={fmt_p(sp['p'])}  n={sp['n']}")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def write_finding(one: dict, n_rows: list[dict], tests: pd.DataFrame) -> None:
    def pick(axis: str) -> dict:
        hit = tests[tests["axis"] == axis]
        return hit.iloc[0].to_dict()

    cd = pick("CD274")
    ha = pick("HLA-A")
    hb = pick("HLA-B")
    hc = pick("HLA-C")
    mi = pick("MHCI")
    cd8 = pick("CD8A")
    b2m = pick("B2M")
    trop = pick("TACSTD2")

    n = int(one["n"])
    n_geo = int(one["n_geo_columns"])
    n_complete = int(one["n_complete_cldn4_cd274_hla"])

    md = f"""# FINDING — GSE253564 leftover CLDN4 vs CD274 / HLA continuous

**Additive only. CLDN4 only.** Public pre-treatment FPKM from the leftover neoadjuvant **durvalumab ± SBRT** series ([GSE253564](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253564); Altorki et al., *Cell Reports Medicine* 2024, PMID 38401548). This folder does **not** run GSEA (that leftover is `methods/gse253564_cldn4_gsea`) and does **not** re-audit MPR/PFS (`opus_geo_leftover`). TACSTD2 is a companion column only.

**Question.** Does continuous CLDN4 track PD-L1 (`CD274`) or classical MHC-I (`HLA-A`, `HLA-B`, `HLA-C`) on this leftover ICI-lung matrix? Partial Spearman residualizes ranks on `CD8A`.

**Scale.** Author FPKM → `log2(FPKM+1)`. MHC-I mean-z = gene-wise z of HLA-A/B/C only (**3/3**), not B2M/TAP. Patient/tumour is the unit (one pre-treatment column per tumour).

Reproduce: `python3 methods/gse253564_cldn4_cd274/analyze.py` (downloads the public FPKM if missing).

---

## Honest n

GEO deposits **{n_geo}** pre-treatment tumours. Continuous Spearman uses every complete pair. Do **not** write the GSEA quartile n (8 vs 8) as the n for these ρ values.

| Item | n |
|---|---|
| GEO pre-treatment FPKM columns | **{n_geo}** |
| Unique tumours (1 column = 1 tumour) | **{n}** |
| CLDN4 + CD274 + HLA-A/B/C non-NA | **{n_complete} / {n_geo}** |
| CD8A non-NA (partial covariate) | **{int(one["n_cd8a"])} / {n_geo}** |
| HLA-I genes in the mean-z | {one["mhci_genes"]} ({one["n_hla_in_score"]}/3) |
| GSEA Q4 vs Q1 (other leftover; not used here) | 8 vs 8 |
| RECIST / DCB on this page | **0** (not an ICI-response test) |

---

## One row (CLDN4, n={n})

| Cohort | n | scale | CD274 ρ (p) | HLA-A ρ (p) | HLA-B ρ (p) | HLA-C ρ (p) | MHC-I mean-z ρ (p) | CD274 \\| CD8A ρ (p) | MHC-I \\| CD8A ρ (p) |
|---|---:|---|---|---|---|---|---|---|---|
| **GSE253564** | **{n}** | log2(FPKM+1) | {fmt_rho(cd["rho"])} ({fmt_p(cd["p"])}) | {fmt_rho(ha["rho"])} ({fmt_p(ha["p"])}) | {fmt_rho(hb["rho"])} ({fmt_p(hb["p"])}) | {fmt_rho(hc["rho"])} ({fmt_p(hc["p"])}) | {fmt_rho(mi["rho"])} ({fmt_p(mi["p"])}) | {fmt_rho(cd["partial_cd8a_rho"])} ({fmt_p(cd["partial_cd8a_p"])}) | {fmt_rho(mi["partial_cd8a_rho"])} ({fmt_p(mi["partial_cd8a_p"])}) |

Full numeric row: `tables/one_row.tsv`. Axis-level tests: `tables/spearman.tsv`.

---

## How to read the row

Continuous CLDN4 vs **CD274** is **null** at n={n} (ρ={fmt_rho(cd["rho"])}, p={fmt_p(cd["p"])}). Classical HLA-A/B/C and the HLA-A/B/C mean-z are also **null** (|ρ|≤{max(abs(ha["rho"]), abs(hb["rho"]), abs(hc["rho"]), abs(mi["rho"])):.3f}; MHC-I ρ={fmt_rho(mi["rho"])}, p={fmt_p(mi["p"])}). Residualising ranks on CD8A does not create a CD274 or MHC-I hit (partial CD274 ρ={fmt_rho(cd["partial_cd8a_rho"])}, p={fmt_p(cd["partial_cd8a_p"])}; partial MHC-I ρ={fmt_rho(mi["partial_cd8a_rho"])}, p={fmt_p(mi["partial_cd8a_p"])}).

The sibling GSEA leftover (`gse253564_cldn4_gsea`) already printed CLDN4 vs CD274 ρ=−0.228, p=0.209, n=32 on the same matrix. This page keeps that CD274 number and adds the **gene-level HLA-A/B/C** row plus the **CD8A partial**. Do not upgrade the GSEA MHC-I/APM set (which mixed B2M/TAP) into an HLA-gene claim. HLA-A/B/C mean-z here is {fmt_rho(mi["rho"])}.

Companion (not a CD274/HLA claim): CLDN4 vs CD8A ρ={fmt_rho(cd8["rho"])}, p={fmt_p(cd8["p"])}; vs B2M ρ={fmt_rho(b2m["rho"])}, p={fmt_p(b2m["p"])}; vs TACSTD2 ρ={fmt_rho(trop["rho"])}, p={fmt_p(trop["p"])}.

This is pre-treatment leftover durvalumab ± SBRT bulk, not an MPR test. CLDN4 vs MPR was already null in `opus_geo_leftover`.

---

## Locked design (before ρ)

| Item | Choice |
|---|---|
| Matrix | Public author FPKM only. No FASTQ / SRA. |
| File | `{FPKM_NAME}` |
| Scale | log2(FPKM+1); duplicate symbols collapsed to the highest-mean locus |
| Unit | one GEO pre-treatment column = one tumour |
| Honest n | complete cases with CLDN4 + CD274 + HLA-A/B/C |
| Predictor | CLDN4. TACSTD2 is a companion only. |
| HLA-I score | unweighted mean of per-gene z-scores for HLA-A, HLA-B, HLA-C (3/3). Not B2M/TAP. |
| Continuous test | two-sided Spearman |
| Partial | first-order rank residual on CD8A |
| Out of scope | GSEA; Q4 vs Q1 NES; GSE248378 post-treatment; MPR/PFS re-test; TACSTD2 splits |

---

## What this does not claim

- Not GSEA. The IFN/MHC/TJ NES leftover is a different agent and a different folder.
- Not a TACSTD2 analysis and not a dual-high gate.
- Not an MPR / PFS / recurrence test.
- Not CLDN4 vs PD-L1 protein or HLA protein. These are RNA rows.
- Not a claim that CLDN4 *induces* PD-L1. Association only. Bulk FFPE mixes epithelium and infiltrate.
- Do not pool this ρ with the four ICI-bulk CD274 rows in `methods/cldn4_cd274_ici`.

---

## Files

- `FINDING.md` — this page
- `analyze.py` — download-once, score, write tables/figures
- `tables/one_row.tsv` — the one-row table
- `tables/spearman.tsv` / `n_table.tsv` / `samples.tsv` / `summary.json`
- `figures/cldn4_vs_cd274.png`, `cldn4_vs_hla_i.png`
"""
    (HERE / "FINDING.md").write_text(md)
    log(f"wrote {HERE / 'FINDING.md'}")


def main() -> None:
    path = download_fpkm()
    expr = load_expr(path)
    n_geo = int(expr.shape[1])
    log(f"matrix genes={expr.shape[0]} tumours={n_geo}")

    resolved = {canon: resolve_gene(expr.index, aliases) for canon, aliases in GENES.items()}
    missing_core = [g for g in ["CLDN4", "CD274", *HLA_I] if resolved[g] is None]
    if missing_core:
        raise SystemExit(f"missing core genes on GSE253564 FPKM: {missing_core}")

    scores = pd.DataFrame({canon: expr.loc[hit].astype(float) for canon, hit in resolved.items() if hit})
    scores.index.name = "sample"
    mhci, hla_used = mean_z(expr, [resolved[g] for g in HLA_I if resolved[g]])
    scores["MHCI"] = mhci.reindex(scores.index)

    complete = scores[["CLDN4", "CD274", *HLA_I]].dropna()
    n_complete = int(len(complete))
    n = n_complete
    if n != n_geo:
        log(f"WARNING complete-case dropout: {n} of {n_geo}")

    tests = []
    for axis in AXES + ["CD8A", "B2M", "TACSTD2"]:
        if axis not in scores:
            tests.append(
                {
                    "cohort": ACCESSION,
                    "axis": axis,
                    "n": 0,
                    "rho": np.nan,
                    "p": np.nan,
                    "partial_cd8a_rho": np.nan,
                    "partial_cd8a_p": np.nan,
                    "cd8a_present": False,
                    "note": "axis missing",
                }
            )
            continue
        sp = spearman(scores["CLDN4"], scores[axis])
        if "CD8A" in scores and scores["CD8A"].notna().sum() >= 5:
            psp = partial_spearman(scores["CLDN4"], scores[axis], scores["CD8A"])
        else:
            psp = {"n": 0, "rho": np.nan, "p": np.nan}
        tests.append(
            {
                "cohort": ACCESSION,
                "axis": axis,
                "n": sp["n"],
                "rho": sp["rho"],
                "p": sp["p"],
                "partial_cd8a_rho": psp["rho"],
                "partial_cd8a_p": psp["p"],
                "cd8a_present": bool("CD8A" in scores and scores["CD8A"].notna().sum() >= 5),
                "note": "companion" if axis in {"CD8A", "B2M", "TACSTD2"} else "primary",
            }
        )
    tests_df = pd.DataFrame(tests)
    tests_df.to_csv(TAB / "spearman.tsv", sep="\t", index=False)

    by_axis = tests_df.set_index("axis")
    one = {
        "cohort": ACCESSION,
        "n": n,
        "n_geo_columns": n_geo,
        "n_complete_cldn4_cd274_hla": n_complete,
        "n_cd8a": int(scores["CD8A"].notna().sum()) if "CD8A" in scores else 0,
        "scale": "log2(FPKM+1)",
        "tissue": "NSCLC tumor, neoadjuvant durvalumab ± SBRT, pre-treatment FFPE",
        "assay": "RNA-seq FPKM (author processed)",
        "mhci_genes": ",".join(HLA_I),
        "n_hla_in_score": len(hla_used),
        "CD274_rho": float(by_axis.loc["CD274", "rho"]),
        "CD274_p": float(by_axis.loc["CD274", "p"]),
        "HLA-A_rho": float(by_axis.loc["HLA-A", "rho"]),
        "HLA-A_p": float(by_axis.loc["HLA-A", "p"]),
        "HLA-B_rho": float(by_axis.loc["HLA-B", "rho"]),
        "HLA-B_p": float(by_axis.loc["HLA-B", "p"]),
        "HLA-C_rho": float(by_axis.loc["HLA-C", "rho"]),
        "HLA-C_p": float(by_axis.loc["HLA-C", "p"]),
        "MHCI_rho": float(by_axis.loc["MHCI", "rho"]),
        "MHCI_p": float(by_axis.loc["MHCI", "p"]),
        "CD274_partial_cd8a_rho": float(by_axis.loc["CD274", "partial_cd8a_rho"]),
        "CD274_partial_cd8a_p": float(by_axis.loc["CD274", "partial_cd8a_p"]),
        "MHCI_partial_cd8a_rho": float(by_axis.loc["MHCI", "partial_cd8a_rho"]),
        "MHCI_partial_cd8a_p": float(by_axis.loc["MHCI", "partial_cd8a_p"]),
        "resolved": {k: (v if v is not None else "") for k, v in resolved.items()},
        "file": FPKM_NAME,
        "url": FPKM_URL,
        "sha256": sha256(path),
        "bytes": int(path.stat().st_size),
        "note": "CLDN4 only; leftover Durva pre-treatment bulk; continuous CD274/HLA; no GSEA",
    }
    pd.DataFrame([one]).to_csv(TAB / "one_row.tsv", sep="\t", index=False)

    n_rows = [
        {"item": "GEO pre-treatment FPKM columns", "n": n_geo},
        {"item": "Unique tumours (1 column = 1 tumour)", "n": n},
        {"item": "CLDN4 + CD274 + HLA-A/B/C non-NA", "n": n_complete},
        {"item": "CD8A non-NA (partial covariate)", "n": one["n_cd8a"]},
        {"item": "HLA-I genes in the mean-z", "n": ",".join(HLA_I)},
        {"item": "GSEA Q4 vs Q1 (other leftover; not used here)", "n": "8 vs 8"},
        {"item": "RECIST / DCB on this page", "n": 0},
    ]
    pd.DataFrame(n_rows).to_csv(TAB / "n_table.tsv", sep="\t", index=False)

    sample_tbl = scores.reset_index()
    sample_tbl.to_csv(TAB / "samples.tsv", sep="\t", index=False)
    (TAB / "summary.json").write_text(json.dumps(one, indent=2, default=str) + "\n")

    save_scatter(
        FIG / "cldn4_vs_cd274.png",
        scores["CLDN4"],
        scores["CD274"],
        "GSE253564 leftover — CLDN4 vs CD274",
        "CLDN4 log2(FPKM+1)",
        "CD274 log2(FPKM+1)",
    )
    save_scatter(
        FIG / "cldn4_vs_hla_i.png",
        scores["CLDN4"],
        scores["MHCI"],
        "GSE253564 leftover — CLDN4 vs HLA-A/B/C mean-z",
        "CLDN4 log2(FPKM+1)",
        "HLA-I mean-z",
    )

    write_finding(one, n_rows, tests_df)
    log("done")
    print(tests_df.to_string(index=False))
    print(pd.DataFrame([one])[["cohort", "n", "CD274_rho", "CD274_p", "MHCI_rho", "MHCI_p"]].to_string(index=False))


if __name__ == "__main__":
    main()
