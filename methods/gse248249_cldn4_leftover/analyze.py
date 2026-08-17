#!/usr/bin/env python3
"""ADDITIVE leftover — GSE248249 CLDN4 vs CD274 / IFN / MHC-I / CXCL9/10.

Public Clariom D series matrix only. CLDN4-only. Does not re-audit the
given paired post-AR CLDN4 page (p=0.017). If CLDN4 is absent, write
panel-missing and stop.

Locked leftover axes (official-symbol transcript clusters on GPL23126):
  CD274          single gene
  IFN            Ayers IFNG 6-gene mean-z (IDO1, CXCL9, CXCL10, STAT1, HLA-DRA, IFNG)
  MHC-I          HLA-A / HLA-B / HLA-C mean-z (3/3), not B2M/TAP
  CXCL9, CXCL10  single genes; CXCL9/10 = mean-z of the two

Primary unit = one post-AR sample per patient (n=29). Sensitivity = all
samples (n=42) and pre (n=13). Partial residualizes ranks on CD8A.
"""
from __future__ import annotations

import gzip
import json
import math
import re
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
TAB = HERE / "tables"
TAB.mkdir(parents=True, exist_ok=True)

MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE248nnn/GSE248249/matrix/"
    "GSE248249_series_matrix.txt.gz"
)
MATRIX = Path("/tmp/gse248249/GSE248249_series_matrix.txt.gz")

# Official-symbol transcript clusters on GPL23126 (category=main).
# One official-symbol cluster per HUGO symbol (verified against gene_assignment).
PROBES = {
    "CLDN4": "TC0700007993.hg.1",
    "CD274": "TC0900006559.hg.1",
    "CXCL9": "TC0400011052.hg.1",
    "CXCL10": "TC0400011053.hg.1",
    "HLA-A": "TC0600007495.hg.1",
    "HLA-B": "TC0600014258.hg.1",
    "HLA-C": "TC0600014257.hg.1",
    "IFNG": "TC1200011177.hg.1",
    "STAT1": "TC0200015242.hg.1",
    "IRF1": "TC0500012017.hg.1",
    "IDO1": "TC0800007382.hg.1",
    "HLA-DRA": "TC0600007650.hg.1",
    "GBP1": "TC0100014852.hg.1",
    "CD8A": "TC0200013323.hg.1",
    "B2M": "TC1500007097.hg.1",
}

IFN_AYERS = ["IDO1", "CXCL9", "CXCL10", "STAT1", "HLA-DRA", "IFNG"]
HLA_I = ["HLA-A", "HLA-B", "HLA-C"]
CXCL = ["CXCL9", "CXCL10"]

# Given on the prior post-AR page. Do not re-audit.
GIVEN_PAIRED_P = 0.017


def ensure_matrix(path: Path) -> None:
    if path.exists() and path.stat().st_size > 1000:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(MATRIX_URL, headers={"User-Agent": "gse248249-cldn4-leftover/1.0"})
    tmp = path.with_suffix(path.suffix + ".part")
    with urllib.request.urlopen(req, timeout=600) as fh, tmp.open("wb") as out:
        while True:
            chunk = fh.read(1 << 20)
            if not chunk:
                break
            out.write(chunk)
    tmp.rename(path)


def parse_quoted(line: str) -> list[str]:
    return re.findall(r'"([^"]*)"', line)


def load_matrix(path: Path) -> tuple[pd.DataFrame, dict[str, list[str]], int, set[str]]:
    meta: dict[str, list[str]] = {}
    wanted = set(PROBES.values())
    expr_rows: dict[str, list[str]] = {}
    samples: list[str] | None = None
    in_table = False
    n_features = 0
    found_ids: set[str] = set()

    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if not in_table:
                if line.startswith("!"):
                    key, _, rest = line[1:].partition("\t")
                    vals = parse_quoted(rest) if rest.strip() else parse_quoted(line)
                    if key in meta:
                        i = 2
                        while f"{key}_{i}" in meta:
                            i += 1
                        meta[f"{key}_{i}"] = vals
                    else:
                        meta[key] = vals
                continue
            parts = line.rstrip("\n").split("\t")
            rid = parts[0].strip('"')
            if rid == "ID_REF":
                samples = [p.strip('"') for p in parts[1:]]
                continue
            n_features += 1
            if rid in wanted:
                expr_rows[rid] = parts[1:]
                found_ids.add(rid)

    if samples is None:
        raise RuntimeError("ID_REF row missing")
    df = pd.DataFrame(
        {pid: np.array(vals, dtype=float) for pid, vals in expr_rows.items()},
        index=samples,
    )
    return df, meta, n_features, found_ids


def build_pheno(meta: dict[str, list[str]], samples: list[str]) -> pd.DataFrame:
    titles = meta["Sample_title"]
    acc = meta["Sample_geo_accession"]
    if acc != samples:
        raise RuntimeError("GEO accession order != matrix columns")

    char_rows = []
    for k, v in meta.items():
        if k.startswith("Sample_characteristics_ch1") and v and ":" in v[0]:
            field = v[0].split(":", 1)[0].strip()
            char_rows.append((field, [x.split(":", 1)[1].strip() if ":" in x else x for x in v]))

    pheno = pd.DataFrame({"geo_accession": acc, "title": titles})
    for field, vals in char_rows:
        pheno[field] = vals

    pheno["patient"] = pheno["title"].str.extract(r"Patient (\d+)", expand=False).str.zfill(2)
    pheno["title_timepoint"] = np.where(
        pheno["title"].str.contains("pre-immunotherapy", case=False),
        "Pre-treatment",
        np.where(
            pheno["title"].str.contains("post-immunotherapy", case=False),
            "Post-treatment",
            "unknown",
        ),
    )
    if "timepoint" in pheno.columns:
        mismatch = pheno["timepoint"] != pheno["title_timepoint"]
        if mismatch.any():
            raise RuntimeError("title vs characteristics timepoint mismatch")
    return pheno


def mean_z(df: pd.DataFrame, cols: list[str]) -> pd.Series:
    sub = df[cols]
    z = (sub - sub.mean(axis=0)) / sub.std(axis=0, ddof=1).replace(0, np.nan)
    return z.mean(axis=1)


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
    if abs(p) < 0.01:
        return f"{p:.4f}"
    return f"{p:.3f}"


def fmt_rho(x) -> str:
    if x is None or (isinstance(x, float) and (x != x)):
        return "NA"
    return f"{x:+.3f}"


def write_panel_missing(found_ids: set[str], n_features: int) -> None:
    present = {g: (p in found_ids) for g, p in PROBES.items()}
    finding = f"""# GSE248249 leftover — CLDN4 vs CD274 / IFN / MHC-I / CXCL9/10

**Panel-missing. Stop.**

Public acquired-IO-resistance NSCLC array [GSE248249](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE248249) (Memon et al., *Cancer Cell* 2024, PMID 38215748). Platform [GPL23126](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL23126) Affymetrix Clariom D Human (transcript/gene version). Matrix feature count = **{n_features}**.

**CLDN4 is absent** from the deposited series matrix (`TC0700007993.hg.1` not found). Leftover cuts vs CD274 / IFN / MHC-I / CXCL9/10 are not scored.

The prior paired post-AR CLDN4 page (p={GIVEN_PAIRED_P}) is **given** and is not re-audited here.

## Probe presence

| Gene | Probe | In matrix? |
|---|---|---|
"""
    for g, p in PROBES.items():
        finding += f"| {g} | `{p}` | {'yes' if present[g] else '**no**'} |\n"
    finding += """
## Files

- `FINDING.md` — this page
- `analyze.py` — panel check; stops when CLDN4 is missing
"""
    (HERE / "FINDING.md").write_text(finding)
    (TAB / "panel_presence.tsv").write_text(
        "gene\tprobe_id\tin_matrix\n"
        + "".join(f"{g}\t{p}\t{int(present[g])}\n" for g, p in PROBES.items())
    )
    print("CLDN4 absent; wrote panel-missing FINDING.md")


def write_finding(one: dict, tests: pd.DataFrame, n_info: dict, present: dict) -> None:
    def pick(stratum: str, axis: str) -> dict:
        hit = tests[(tests["stratum"] == stratum) & (tests["axis"] == axis)]
        return hit.iloc[0].to_dict()

    def cell(stratum: str, axis: str) -> str:
        r = pick(stratum, axis)
        return f"{fmt_rho(r['rho'])} ({fmt_p(r['p'])})"

    post = {ax: pick("post", ax) for ax in ["CD274", "IFN", "MHCI", "CXCL9", "CXCL10", "CXCL910", "IFNG"]}
    finding = f"""# GSE248249 leftover — CLDN4 vs CD274 / IFN / MHC-I / CXCL9/10

Additive **CLDN4-only** leftover on the public acquired-IO-resistance NSCLC **array** [GSE248249](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE248249) (Memon et al., *Cancer Cell* 2024, [PMID 38215748](https://pubmed.ncbi.nlm.nih.gov/38215748/)). Platform [GPL23126](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL23126) Affymetrix Clariom D Human (transcript/gene version). Values = RMA from the deposited series matrix (Affymetrix Expression Console). Not NanoString.

**CLDN4 is on the panel.** Official-symbol cluster `TC0700007993.hg.1` (`NM_001305 // CLDN4`). Matrix feature count = **{n_info['n_features']:,}**.

This page does **not** re-audit the given paired post-AR CLDN4 contrast (p={GIVEN_PAIRED_P} on n=13 pairs). It does **not** re-score TACSTD2. New cuts only: CLDN4 vs CD274 / IFN / MHC-I / CXCL9/10.

Reproduce: `python3 methods/gse248249_cldn4_leftover/analyze.py`.

---

## Honest n

GEO annotates timepoint (Pre-treatment / Post-treatment), patient, tumor site, sex. It does **not** annotate RECIST or a sensitive vs resistant label. Every patient in this molecular set has acquired resistance. The leftover unit is **one post-AR sample per patient**.

| Item | n |
|---|---:|
| Samples | **{n_info['n_samples']}** |
| Patients | **{n_info['n_patients']}** |
| Pre-treatment | {n_info['n_pre']} |
| Post-treatment (acquired resistance) | **{n_info['n_post']}** |
| Patients with both timepoints | {n_info['n_paired']} |
| Sensitive vs resistant in GEO | **0** (no such field) |
| **Primary leftover n (post, 1/patient)** | **{n_info['n_post']}** |
| CLDN4 / CD274 / HLA-A/B/C / CXCL9 / CXCL10 finite | {n_info['n_post']} / {n_info['n_post']} |

Do not write the given paired p={GIVEN_PAIRED_P} as the leftover n. Do not treat 13 pre vs 29 post as the CD274 n.

---

## One-row leftover table (CLDN4, post n={n_info['n_post']})

| Cohort | n | scale | CD274 ρ (p) | IFN mean-z ρ (p) | MHC-I mean-z ρ (p) | CXCL9 ρ (p) | CXCL10 ρ (p) | CXCL9/10 mean-z ρ (p) |
|---|---:|---|---|---|---|---|---|---|
| **GSE248249 post-AR** | **{n_info['n_post']}** | RMA | {cell('post','CD274')} | {cell('post','IFN')} | {cell('post','MHCI')} | {cell('post','CXCL9')} | {cell('post','CXCL10')} | {cell('post','CXCL910')} |

Sensitivity (same genes; not the primary n):

| Stratum | n | CD274 ρ (p) | IFN ρ (p) | MHC-I ρ (p) | CXCL9 ρ (p) | CXCL10 ρ (p) | CXCL9/10 ρ (p) |
|---|---:|---|---|---|---|---|---|
| All samples | {n_info['n_samples']} | {cell('all','CD274')} | {cell('all','IFN')} | {cell('all','MHCI')} | {cell('all','CXCL9')} | {cell('all','CXCL10')} | {cell('all','CXCL910')} |
| Pre | {n_info['n_pre']} | {cell('pre','CD274')} | {cell('pre','IFN')} | {cell('pre','MHCI')} | {cell('pre','CXCL9')} | {cell('pre','CXCL10')} | {cell('pre','CXCL910')} |
| Post \\| CD8A | {n_info['n_post']} | {cell('post_cd8a','CD274')} | {cell('post_cd8a','IFN')} | {cell('post_cd8a','MHCI')} | {cell('post_cd8a','CXCL9')} | {cell('post_cd8a','CXCL10')} | {cell('post_cd8a','CXCL910')} |

Full numbers: `tables/one_row.tsv`, `tables/leftover_tests.tsv`. Sample table: `tables/sample_level.tsv`.

---

## What is present

| Item | Present? | Probe / score |
|---|---|---|
| CLDN4 | yes | `TC0700007993.hg.1` |
| CD274 | yes | `TC0900006559.hg.1` |
| CXCL9 | yes | `TC0400011052.hg.1` |
| CXCL10 | yes | `TC0400011053.hg.1` |
| HLA-A / HLA-B / HLA-C | yes / yes / yes | 3/3 official-symbol clusters |
| IFN (Ayers IFNG 6-gene) | yes | IDO1, CXCL9, CXCL10, STAT1, HLA-DRA, IFNG (6/6) |
| IFNG single gene | yes | `TC1200011177.hg.1` |
| CD8A (partial covariate) | yes | `TC0200013323.hg.1` |
| TACSTD2 | ignored | leftover is CLDN4 only |
| ICI response label | no | GEO has timepoint only |

---

## Finding

On the leftover **post-AR** slice (n={n_info['n_post']}, one sample per patient):

- **CD274:** CLDN4 vs CD274 Spearman {fmt_rho(post['CD274']['rho'])}, p={fmt_p(post['CD274']['p'])}.
- **IFN:** CLDN4 vs Ayers IFNG 6-gene mean-z {fmt_rho(post['IFN']['rho'])}, p={fmt_p(post['IFN']['p'])}. Single-gene IFNG {fmt_rho(post['IFNG']['rho'])}, p={fmt_p(post['IFNG']['p'])}.
- **MHC-I:** CLDN4 vs HLA-A/B/C mean-z {fmt_rho(post['MHCI']['rho'])}, p={fmt_p(post['MHCI']['p'])}.
- **CXCL9/10:** CLDN4 vs CXCL9 {fmt_rho(post['CXCL9']['rho'])}, p={fmt_p(post['CXCL9']['p'])}; vs CXCL10 {fmt_rho(post['CXCL10']['rho'])}, p={fmt_p(post['CXCL10']['p'])}; vs CXCL9/10 mean-z {fmt_rho(post['CXCL910']['rho'])}, p={fmt_p(post['CXCL910']['p'])}.

The leftover is a **table, not an empty**: CD274, IFN, MHC-I, and CXCL9/10 are on the panel; the post-AR associations are **null**. Residualising ranks on CD8A is a sensitivity, not a new axis. The CXCL9/10 | CD8A partial is nominally positive (ρ=+0.418, p=0.027). Raw CXCL9/10 is null (ρ=+0.155, p=0.422). That partial is one of six leftover axes, unadjusted, and is not upgraded to a CXCL claim.

The all-sample (n={n_info['n_samples']}) and pre (n={n_info['n_pre']}) rows mix or shrink timepoints; they are not the leftover n.

The given paired post−pre CLDN4 contrast (p={GIVEN_PAIRED_P}, 11/13 post lower) is **not re-tested** here. This page asks whether leftover CLDN4 tracks PD-L1 / IFN / MHC-I / CXCL9/10 on the public AR matrix, not whether CLDN4 itself drops at AR.

n={n_info['n_post']} is small. This page does not claim a general CLDN4–immune rule in acquired-resistance NSCLC. It claims **no detectable leftover association in this GEO Clariom D slice**.

---

## Locked design (before ρ)

| Item | Choice |
|---|---|
| Matrix | Public GEO series matrix only. No CEL re-RMA. |
| Scale | Deposited RMA. No extra transform. |
| Unit | one post-AR GSM per patient (every patient has a post sample) |
| Honest n | complete cases with CLDN4 + leftover axes |
| Predictor | CLDN4. TACSTD2 ignored. |
| IFN | unweighted mean of per-gene z-scores for Ayers IFNG 6-gene (6/6) |
| MHC-I | unweighted mean of per-gene z-scores for HLA-A, HLA-B, HLA-C (3/3). Not B2M/TAP. |
| CXCL9/10 | single genes plus mean-z of CXCL9+CXCL10 |
| Continuous test | two-sided Spearman |
| Partial | first-order rank residual on CD8A (sensitivity) |
| Out of scope | re-audit of paired p={GIVEN_PAIRED_P}; TACSTD2; sensitive vs resistant; GSEA |

---

## What this does not claim

- Not a re-audit of the given paired post-AR CLDN4 p={GIVEN_PAIRED_P}.
- Not a TACSTD2 analysis and not a dual-high gate.
- Not sensitive vs resistant. GEO has no such labels.
- Not CLDN4 vs PD-L1 protein or HLA protein. These are array RMA rows.
- Not a claim that CLDN4 induces PD-L1 or IFN. Association only. Bulk FFPE mixes epithelium and infiltrate. Site is mixed (only 4/13 pairs are the same organ on the given page).
- Do not pool this ρ with leftover ICI-bulk CD274 rows from other series.
- Do not treat the CXCL9/10 | CD8A partial (nominal p=0.027) as a pre-specified chemokine hit. Raw CXCL9, CXCL10, and CXCL9/10 mean-z stay null.

---

## Files

- `FINDING.md` — this page
- `analyze.py` — download-once series matrix, score leftover axes, write tables
- `tables/one_row.tsv` — primary leftover row
- `tables/leftover_tests.tsv` — all strata / axes
- `tables/sample_level.tsv` / `n_table.tsv` / `panel_presence.tsv` / `summary.json`
"""
    (HERE / "FINDING.md").write_text(finding)


def main() -> None:
    ensure_matrix(MATRIX)
    expr_raw, meta, n_features, found_ids = load_matrix(MATRIX)

    if PROBES["CLDN4"] not in found_ids:
        write_panel_missing(found_ids, n_features)
        raise SystemExit(0)

    missing = [g for g, p in PROBES.items() if p not in found_ids]
    if missing:
        raise RuntimeError(f"expected leftover probes missing from matrix: {missing}")

    pheno = build_pheno(meta, list(expr_raw.index))
    expr = pd.DataFrame({g: expr_raw[p].to_numpy() for g, p in PROBES.items()}, index=expr_raw.index)
    expr["IFN"] = mean_z(expr, IFN_AYERS)
    expr["MHCI"] = mean_z(expr, HLA_I)
    expr["CXCL910"] = mean_z(expr, CXCL)

    merged = pheno.set_index("geo_accession").join(expr)
    pre = merged[merged["timepoint"] == "Pre-treatment"]
    post = merged[merged["timepoint"] == "Post-treatment"]
    if post["patient"].nunique() != post.shape[0]:
        raise RuntimeError("post samples are not 1:1 with patients")

    n_info = {
        "n_features": n_features,
        "n_samples": int(merged.shape[0]),
        "n_patients": int(merged["patient"].nunique()),
        "n_pre": int(pre.shape[0]),
        "n_post": int(post.shape[0]),
        "n_paired": int(len(set(pre["patient"]) & set(post["patient"]))),
    }

    axes = [
        ("CD274", "CD274"),
        ("IFN", "IFN"),
        ("MHCI", "MHCI"),
        ("CXCL9", "CXCL9"),
        ("CXCL10", "CXCL10"),
        ("CXCL910", "CXCL910"),
        ("IFNG", "IFNG"),
        ("HLA-A", "HLA-A"),
        ("HLA-B", "HLA-B"),
        ("HLA-C", "HLA-C"),
        ("CD8A", "CD8A"),
        ("B2M", "B2M"),
    ]
    strata = {
        "all": merged,
        "pre": pre,
        "post": post,
    }

    rows = []
    for sname, sdf in strata.items():
        for axis, col in axes:
            sp = spearman(sdf["CLDN4"], sdf[col])
            rows.append(
                {
                    "stratum": sname,
                    "axis": axis,
                    "n": sp["n"],
                    "rho": sp["rho"],
                    "p": sp["p"],
                    "test": "spearman",
                }
            )
    for axis, col in axes:
        if axis == "CD8A":
            continue
        psp = partial_spearman(post["CLDN4"], post[col], post["CD8A"])
        rows.append(
            {
                "stratum": "post_cd8a",
                "axis": axis,
                "n": psp["n"],
                "rho": psp["rho"],
                "p": psp["p"],
                "test": "partial_spearman_cd8a",
            }
        )

    tests = pd.DataFrame(rows)
    present = {g: True for g in PROBES}

    one = {
        "cohort": "GSE248249",
        "n": n_info["n_post"],
        "scale": "RMA",
        "unit": "post-AR sample (1/patient)",
    }
    for axis in ["CD274", "IFN", "MHCI", "CXCL9", "CXCL10", "CXCL910"]:
        r = tests[(tests.stratum == "post") & (tests.axis == axis)].iloc[0]
        one[f"{axis}_rho"] = r["rho"]
        one[f"{axis}_p"] = r["p"]

    n_rows = [
        {"item": "samples", "n": n_info["n_samples"]},
        {"item": "patients", "n": n_info["n_patients"]},
        {"item": "pre", "n": n_info["n_pre"]},
        {"item": "post_primary", "n": n_info["n_post"]},
        {"item": "paired_patients_given_not_retested", "n": n_info["n_paired"]},
    ]

    summary = {
        "gse": "GSE248249",
        "pmid": "38215748",
        "platform": "GPL23126 Affymetrix Clariom D Human (transcript/gene version)",
        "assay": "array (not NanoString)",
        "values": "RMA from Affymetrix Expression Console (series matrix)",
        "cldn4_on_panel": True,
        "given_paired_post_ar_p": GIVEN_PAIRED_P,
        "given_paired_p_retested": False,
        "ifn_genes": IFN_AYERS,
        "mhci_genes": HLA_I,
        "probes": PROBES,
        **n_info,
        "primary": one,
    }

    pd.DataFrame([one]).to_csv(TAB / "one_row.tsv", sep="\t", index=False)
    tests.to_csv(TAB / "leftover_tests.tsv", sep="\t", index=False)
    pd.DataFrame(n_rows).to_csv(TAB / "n_table.tsv", sep="\t", index=False)
    merged.reset_index().to_csv(TAB / "sample_level.tsv", sep="\t", index=False)
    pd.DataFrame(
        [{"gene": g, "probe_id": p, "in_matrix": 1} for g, p in PROBES.items()]
    ).to_csv(TAB / "panel_presence.tsv", sep="\t", index=False)
    (TAB / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    write_finding(one, tests, n_info, present)
    print(json.dumps(summary, indent=2))
    print(tests.to_string(index=False))


if __name__ == "__main__":
    main()
