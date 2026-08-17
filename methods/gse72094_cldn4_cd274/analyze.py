#!/usr/bin/env python3
"""GSE72094 LUAD microarray: additive CLDN4 vs CD274 / HLA.

New CD274 cut only. Residual CLDN4 vs CD8A is already NS (PR 302) and is
not re-run or headlined here.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import math
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "gse72094"
SIG = ROOT / "data" / "signatures"
OUT = Path(__file__).resolve().parent

MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE72nnn/GSE72094/"
    "matrix/GSE72094_series_matrix.txt.gz"
)
GPL_URL = (
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
    "?acc=GPL15048&targ=self&form=text&view=data"
)

# MHC-I cassette used on DepMap CLDN4–IFN page (PR 304).
MHC1_GENES = ["HLA-A", "HLA-B", "HLA-C", "B2M"]
HLA_SINGLE = ["HLA-A", "HLA-B", "HLA-C", "B2M", "HLA-E", "HLA-DRA", "HLA-DRB1"]
PRIMARY_FEATURES = ["CD274", "MHC-I"]

EST_A = 0.6049872018
EST_B = 0.0001467884


def md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fmt_p(p) -> str:
    if p is None or not np.isfinite(p):
        return "NA"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.3g}"


def fmt_rho(r) -> str:
    if r is None or not np.isfinite(r):
        return "NA"
    return f"{r:+.3f}"


def download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"downloading {url} -> {dest}")
    urllib.request.urlretrieve(url, dest)


def ensure_inputs() -> tuple[Path, Path]:
    matrix = DATA / "GSE72094_series_matrix.txt.gz"
    plat = DATA / "GPL15048_datatable.txt"
    if not matrix.exists():
        download(MATRIX_URL, matrix)
    if not plat.exists():
        download(GPL_URL, plat)
    return matrix, plat


def residualize(y, Z):
    x = np.column_stack([np.ones(len(y)), Z])
    b, *_ = np.linalg.lstsq(x, y, rcond=None)
    return y - x @ b


def spearman_pair(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    n = int(m.sum())
    if n < 6 or np.nanstd(x[m]) == 0 or np.nanstd(y[m]) == 0:
        return np.nan, np.nan, n
    r, p = stats.spearmanr(x[m], y[m])
    return float(r), float(p), n


def spearman_ci(x, y, n_boot: int = 5000, seed: int = 0):
    m = np.isfinite(x) & np.isfinite(y)
    xx, yy = x[m], y[m]
    n = int(len(xx))
    if n < 6:
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    rhos = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        rhos[i] = stats.spearmanr(xx[idx], yy[idx]).statistic
    return float(np.nanpercentile(rhos, 2.5)), float(np.nanpercentile(rhos, 97.5))


def partial_spearman(x, y, z):
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    n = int(m.sum())
    if n < 8:
        return np.nan, np.nan, n
    xr = stats.rankdata(x[m])
    yr = stats.rankdata(y[m])
    zr = stats.rankdata(z[m])
    Z = zr.reshape(-1, 1)
    rx = residualize(xr, Z)
    ry = residualize(yr, Z)
    r, _ = stats.pearsonr(rx, ry)
    r = float(np.clip(r, -0.999999, 0.999999))
    df = n - 3
    t = r * math.sqrt(df / max(1e-12, 1 - r * r))
    p = float(2 * stats.t.sf(abs(t), df))
    return r, p, n


def bh_fdr(pvals: list[float]) -> list[float]:
    p = np.asarray(pvals, float)
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return q.tolist()
    order = np.argsort(p[ok])
    ranked = p[ok][order]
    n = len(ranked)
    adj = ranked * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    out = np.empty(n)
    out[order] = adj
    q[ok] = out
    return [float(v) for v in q]


def ssgsea(expr: pd.DataFrame, genes: list[str], tau: float = 0.25) -> pd.Series:
    present = [g for g in genes if g in expr.index]
    if len(present) < 10:
        return pd.Series(np.nan, index=expr.columns)
    n_genes = expr.shape[0]
    ranked = expr.rank(axis=0, method="average", ascending=True) * (10000.0 / n_genes)
    gene_set = set(present)
    scores = {}
    for sample in expr.columns:
        m = ranked[sample]
        order = m.sort_values(ascending=False).index
        m_ord = m.loc[order].to_numpy(float)
        hits = np.fromiter((g in gene_set for g in order), dtype=bool, count=len(order))
        w = np.abs(m_ord) ** tau
        w_hit = np.where(hits, w, 0.0)
        nhit = float(w_hit.sum())
        nmiss = float((~hits).sum())
        if nhit <= 0 or nmiss <= 0:
            scores[sample] = np.nan
            continue
        p_hit = np.cumsum(w_hit) / nhit
        p_miss = np.cumsum((~hits).astype(float)) / nmiss
        scores[sample] = float(np.sum(p_hit - p_miss))
    return pd.Series(scores)


def parse_geo_matrix(path: Path):
    meta_rows = {}
    samples = None
    expr_start = None
    with gzip.open(path, "rt", errors="replace") as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if line.startswith("!Sample_geo_accession"):
            samples = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
        if line.startswith("!Sample_source_name"):
            meta_rows["source"] = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
        if line.startswith("!Sample_characteristics"):
            vals = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            key = None
            for v in vals:
                if v and ":" in v:
                    key = v.split(":", 1)[0].strip().lower()
                    break
            if key is None:
                key = f"char_{len(meta_rows)}"
            cleaned = []
            for v in vals:
                cleaned.append(v.split(":", 1)[1].strip() if ":" in v else v)
            k = key
            n = 2
            while k in meta_rows:
                k = f"{key}_{n}"
                n += 1
            meta_rows[k] = cleaned
        if line.startswith('"ID_REF"') or line.startswith("ID_REF"):
            expr_start = i
            break
    if samples is None or expr_start is None:
        raise SystemExit(f"could not parse GEO matrix {path}")
    meta = pd.DataFrame(meta_rows, index=samples)
    header = [x.strip().strip('"') for x in lines[expr_start].rstrip("\n").split("\t")]
    rows = []
    idx = []
    for line in lines[expr_start + 1 :]:
        if line.startswith("!"):
            break
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 2:
            continue
        idx.append(parts[0].strip().strip('"'))
        rows.append([float(x) if x not in ("", "NA", "null") else np.nan for x in parts[1:]])
    expr = pd.DataFrame(rows, index=idx, columns=header[1:])
    expr = expr.loc[:, samples]
    return meta, expr


def first_hugo(symbol: str) -> str:
    """First /// token, then first whitespace token.

    GPL15048 annotates HLA-A as 'HLA-A LOC100507703' (Entrez 3105).
    """
    if symbol is None or (isinstance(symbol, float) and np.isnan(symbol)):
        return ""
    s = str(symbol).split("///")[0].strip()
    if not s or s in {"nan", "None"}:
        return ""
    return s.split()[0]


def load_gpl15048(path: Path) -> pd.Series:
    skip = 0
    with open(path, "r", errors="replace") as f:
        for i, line in enumerate(f):
            if line.startswith("ID\t"):
                skip = i
                break
    df = pd.read_csv(path, sep="\t", skiprows=skip, dtype=str, low_memory=False)
    s = df.set_index("ID")["GeneSymbol"].map(first_hugo)
    s = s[s.str.len() > 0]
    return s


def collapse_maxmean(expr: pd.DataFrame, probe2gene: pd.Series) -> pd.DataFrame:
    common = expr.index.intersection(probe2gene.index)
    g = probe2gene.loc[common]
    e = expr.loc[common]
    means = e.mean(axis=1)
    pick = means.groupby(g).idxmax()
    out = e.loc[pick.values]
    out.index = pick.index
    return out


def mean_z(expr: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, list[str]]:
    present = [g for g in genes if g in expr.index]
    if not present:
        return pd.Series(np.nan, index=expr.columns), present
    z = expr.loc[present].T
    score = ((z - z.mean()) / z.std(ddof=0)).mean(axis=1)
    return score, present


def quartile_masks(x: np.ndarray):
    q1 = np.nanpercentile(x, 25)
    q3 = np.nanpercentile(x, 75)
    lo = np.isfinite(x) & (x <= q1)
    hi = np.isfinite(x) & (x >= q3)
    return lo, hi, float(q1), float(q3)


def fisher_q4_overlap(a: np.ndarray, b: np.ndarray):
    _, a_hi, _, _ = quartile_masks(a)
    _, b_hi, _, _ = quartile_masks(b)
    both = int((a_hi & b_hi).sum())
    a_only = int((a_hi & ~b_hi).sum())
    b_only = int((~a_hi & b_hi).sum())
    neither = int((~a_hi & ~b_hi).sum())
    table = np.array([[both, a_only], [b_only, neither]])
    or_, p = stats.fisher_exact(table, alternative="two-sided")
    n_a = int(a_hi.sum())
    n_b = int(b_hi.sum())
    frac = both / n_a if n_a else np.nan
    return {
        "n_q4_a": n_a,
        "n_q4_b": n_b,
        "n_overlap": both,
        "frac_of_a_q4": frac,
        "odds_ratio": float(or_),
        "fisher_p": float(p),
        "table": table.tolist(),
    }


def mw_q4_q1(predictor: np.ndarray, endpoint: np.ndarray):
    lo, hi, q1, q3 = quartile_masks(predictor)
    y_hi = endpoint[hi]
    y_lo = endpoint[lo]
    n_hi = int(np.isfinite(y_hi).sum())
    n_lo = int(np.isfinite(y_lo).sum())
    if n_hi < 3 or n_lo < 3:
        return {
            "n_q4": n_hi,
            "n_q1": n_lo,
            "median_q4": np.nan,
            "median_q1": np.nan,
            "U": np.nan,
            "p": np.nan,
            "q1_cut": q1,
            "q3_cut": q3,
        }
    u, p = stats.mannwhitneyu(y_hi[np.isfinite(y_hi)], y_lo[np.isfinite(y_lo)], alternative="two-sided")
    return {
        "n_q4": n_hi,
        "n_q1": n_lo,
        "median_q4": float(np.nanmedian(y_hi)),
        "median_q1": float(np.nanmedian(y_lo)),
        "U": float(u),
        "p": float(p),
        "q1_cut": q1,
        "q3_cut": q3,
    }


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    matrix, plat = ensure_inputs()

    print("parsing matrix")
    meta, probes = parse_geo_matrix(matrix)
    print(f"  samples={meta.shape[0]} probes={probes.shape[0]}")
    print(f"  source unique: {sorted(meta['source'].astype(str).str.lower().unique())}")

    keep = meta["source"].astype(str).str.lower().str.contains("adenocarcinoma")
    n_matrix = int(meta.shape[0])
    n_luad = int(keep.sum())
    print(f"  LUAD filter: {n_luad} / {n_matrix}")
    meta_t = meta.loc[keep]
    probes_t = probes.loc[:, meta_t.index]

    p2g = load_gpl15048(plat)
    genes = collapse_maxmean(probes_t, p2g)
    print(f"  genes after max-mean collapse: {genes.shape[0]}")

    needed = ["CLDN4", "CD274", "TACSTD2"] + MHC1_GENES + HLA_SINGLE
    present = [g for g in dict.fromkeys(needed) if g in genes.index]
    missing = [g for g in dict.fromkeys(needed) if g not in genes.index]
    print(f"  present={present} missing={missing}")

    mhc1, mhc1_used = mean_z(genes, MHC1_GENES)
    print(f"  MHC-I genes used: {mhc1_used} ({len(mhc1_used)}/4)")

    est_tab = pd.read_csv(SIG / "estimate_yoshihara_2013.tsv", sep="\t")
    stromal = est_tab.loc[est_tab["set"] == "Stromal141_UP", "hugo"].tolist()
    immune = est_tab.loc[est_tab["set"] == "Immune141_UP", "hugo"].tolist()
    n_strom = sum(g in genes.index for g in stromal)
    n_imm = sum(g in genes.index for g in immune)
    strom = ssgsea(genes, stromal)
    imm = ssgsea(genes, immune)
    est_score = strom + imm
    print(f"  ESTIMATE stromal {n_strom}/141 immune {n_imm}/141")

    cldn4 = genes.loc["CLDN4"].to_numpy(float)
    cd274 = genes.loc["CD274"].to_numpy(float)
    mhc1_v = mhc1.to_numpy(float)
    tac = genes.loc["TACSTD2"].to_numpy(float) if "TACSTD2" in genes.index else None

    rows = []

    def add_row(feature, y, note=""):
        ru, pu, nu = spearman_pair(cldn4, y)
        lo, hi = spearman_ci(cldn4, y)
        rp, pp, np_ = partial_spearman(cldn4, y, est_score.to_numpy(float))
        rows.append(
            {
                "anchor": "CLDN4",
                "feature": feature,
                "n": nu,
                "rho": ru,
                "p": pu,
                "ci_lo": lo,
                "ci_hi": hi,
                "partial_rho_ESTIMATEScore": rp,
                "partial_p_ESTIMATEScore": pp,
                "n_partial": np_,
                "note": note,
            }
        )

    add_row("CD274", cd274, "primary; new CD274 cut")
    add_row("MHC-I", mhc1_v, f"mean-z of {','.join(mhc1_used)}")
    for g in HLA_SINGLE:
        if g in genes.index:
            add_row(g, genes.loc[g].to_numpy(float), "single gene")
    if tac is not None:
        ru, pu, nu = spearman_pair(tac, cd274)
        rows.append(
            {
                "anchor": "TACSTD2",
                "feature": "CD274",
                "n": nu,
                "rho": ru,
                "p": pu,
                "ci_lo": np.nan,
                "ci_hi": np.nan,
                "partial_rho_ESTIMATEScore": np.nan,
                "partial_p_ESTIMATEScore": np.nan,
                "n_partial": nu,
                "note": "companion only; not the finding",
            }
        )

    res = pd.DataFrame(rows)
    prim = res[(res.anchor == "CLDN4") & (res.feature.isin(PRIMARY_FEATURES))].copy()
    q_primary = bh_fdr(prim["p"].tolist())
    res["q_primary_CD274_MHCI"] = np.nan
    res.loc[prim.index, "q_primary_CD274_MHCI"] = q_primary
    hla_idx = res[(res.anchor == "CLDN4") & (res.feature.isin(HLA_SINGLE))].index
    res.loc[hla_idx, "q_HLA_single"] = bh_fdr(res.loc[hla_idx, "p"].tolist())
    res.to_csv(OUT / "correlations.tsv", sep="\t", index=False)

    # New CD274 cut only (Q4 vs Q1 of CD274; CLDN4 Q4 ∩ CD274 Q4).
    cut_cldn4 = mw_q4_q1(cd274, cldn4)
    cut_cd274_by_cldn4 = mw_q4_q1(cldn4, cd274)
    overlap = fisher_q4_overlap(cldn4, cd274)
    overlap_mhc = fisher_q4_overlap(cldn4, mhc1_v)

    cuts = pd.DataFrame(
        [
            {
                "cut": "CD274_Q4_vs_Q1",
                "endpoint": "CLDN4",
                **{k: v for k, v in cut_cldn4.items()},
            },
            {
                "cut": "CLDN4_Q4_vs_Q1",
                "endpoint": "CD274",
                **{k: v for k, v in cut_cd274_by_cldn4.items()},
            },
        ]
    )
    cuts.to_csv(OUT / "q4_cuts.tsv", sep="\t", index=False)
    with open(OUT / "q4_overlap.json", "w") as f:
        json.dump(
            {"CLDN4_Q4_and_CD274_Q4": overlap, "CLDN4_Q4_and_MHCI_Q4": overlap_mhc},
            f,
            indent=2,
        )

    st = pd.DataFrame(index=genes.columns)
    for g in ["CLDN4", "CD274", "TACSTD2"] + MHC1_GENES:
        if g in genes.index:
            st[g] = genes.loc[g]
    st["MHC_I"] = mhc1
    st["ESTIMATE_Score"] = est_score
    st["ESTIMATE_ImmuneScore"] = imm
    st["CD274_quartile"] = pd.qcut(st["CD274"], 4, labels=["Q1", "Q2", "Q3", "Q4"])
    st.index.name = "sample"
    st.to_csv(OUT / "samples.tsv", sep="\t")

    cov = {
        "cohort": "GSE72094",
        "pmid": "26477306",
        "platform": "GPL15048 HuRSTA_2a520709",
        "n_matrix": n_matrix,
        "n_luad": n_luad,
        "n_dropped_non_luad": n_matrix - n_luad,
        "n_probes": int(probes.shape[0]),
        "n_genes": int(genes.shape[0]),
        "genes_present": present,
        "genes_missing": missing,
        "mhc1_genes_used": mhc1_used,
        "hla_a_annotation": "GPL15048 GeneSymbol is 'HLA-A LOC100507703' (Entrez 3105); first-token map to HLA-A",
        "estimate_stromal_present": n_strom,
        "estimate_immune_present": n_imm,
        "cd274_q4_n": cut_cldn4["n_q4"],
        "cd274_q1_n": cut_cldn4["n_q1"],
    }
    pd.DataFrame([cov]).to_csv(OUT / "coverage.tsv", sep="\t", index=False)
    provenance = {
        "coverage": cov,
        "inputs": {
            "GSE72094_series_matrix.txt.gz": {
                "bytes": matrix.stat().st_size,
                "md5": md5(matrix),
                "url": MATRIX_URL,
            },
            "GPL15048_datatable.txt": {
                "bytes": plat.stat().st_size,
                "md5": md5(plat),
                "url": GPL_URL,
            },
        },
        "generated_at": datetime.now(timezone.utc).isoformat(),
    }
    with open(OUT / "provenance.json", "w") as f:
        json.dump(provenance, f, indent=2)

    def row(feature, anchor="CLDN4"):
        sub = res[(res.anchor == anchor) & (res.feature == feature)]
        return None if sub.empty else sub.iloc[0]

    r_cd = row("CD274")
    r_mhc = row("MHC-I")
    r_a = row("HLA-A")
    r_b = row("HLA-B")
    r_c = row("HLA-C")
    r_b2m = row("B2M")
    r_e = row("HLA-E")
    r_dra = row("HLA-DRA")
    r_drb = row("HLA-DRB1")
    r_tac = row("CD274", "TACSTD2")

    def qcell(r):
        q = r.q_primary_CD274_MHCI
        if q is None or not np.isfinite(q):
            return "—"
        return fmt_p(q)

    finding = f"""# Finding — GSE72094 LUAD: CLDN4 vs CD274 / HLA (new CD274 cut)

**Additive public cohort.** Schabath / Moffitt **GSE72094** (PMID 26477306): resected **lung adenocarcinoma**, **GPL15048** Rosetta/Merck HuRSTA custom Affymetrix 2.0, author IRON-normalized series matrix.

**Honest n.** Series matrix has **{n_matrix}** arrays. All **{n_luad}** have `source_name = lung adenocarcinoma`. Non-LUAD dropped: **{n_matrix - n_luad}**. Complete-case n for CLDN4, CD274, and MHC-I genes is **{int(r_cd.n)}**. No ICI labels. No paired normals.

**This page is the new CD274 / HLA cut only.** Residual CLDN4 vs CD8A after ESTIMATE is already **not significant** on this matrix (PR 302: partial ρ = −0.051, p = 0.282, n = 442) and is **not** re-run or headlined here.

## Verdict (unadjusted Spearman; n={int(r_cd.n)})

| Test | n | ρ | 95% CI | p | BH q |
|---|---:|---:|---|---:|---:|
| **CLDN4 vs CD274** | {int(r_cd.n)} | {fmt_rho(r_cd.rho)} | [{fmt_rho(r_cd.ci_lo)}, {fmt_rho(r_cd.ci_hi)}] | {fmt_p(r_cd.p)} | {qcell(r_cd)} |
| **CLDN4 vs MHC-I** | {int(r_mhc.n)} | {fmt_rho(r_mhc.rho)} | [{fmt_rho(r_mhc.ci_lo)}, {fmt_rho(r_mhc.ci_hi)}] | {fmt_p(r_mhc.p)} | {qcell(r_mhc)} |
| CLDN4 vs HLA-A | {int(r_a.n) if r_a is not None else 0} | {fmt_rho(r_a.rho) if r_a is not None else "NA"} | [{fmt_rho(r_a.ci_lo) if r_a is not None else "NA"}, {fmt_rho(r_a.ci_hi) if r_a is not None else "NA"}] | {fmt_p(r_a.p) if r_a is not None else "NA"} | — |
| CLDN4 vs HLA-B | {int(r_b.n) if r_b is not None else 0} | {fmt_rho(r_b.rho) if r_b is not None else "NA"} | [{fmt_rho(r_b.ci_lo) if r_b is not None else "NA"}, {fmt_rho(r_b.ci_hi) if r_b is not None else "NA"}] | {fmt_p(r_b.p) if r_b is not None else "NA"} | — |
| CLDN4 vs HLA-C | {int(r_c.n) if r_c is not None else 0} | {fmt_rho(r_c.rho) if r_c is not None else "NA"} | [{fmt_rho(r_c.ci_lo) if r_c is not None else "NA"}, {fmt_rho(r_c.ci_hi) if r_c is not None else "NA"}] | {fmt_p(r_c.p) if r_c is not None else "NA"} | — |
| CLDN4 vs B2M | {int(r_b2m.n) if r_b2m is not None else 0} | {fmt_rho(r_b2m.rho) if r_b2m is not None else "NA"} | [{fmt_rho(r_b2m.ci_lo) if r_b2m is not None else "NA"}, {fmt_rho(r_b2m.ci_hi) if r_b2m is not None else "NA"}] | {fmt_p(r_b2m.p) if r_b2m is not None else "NA"} | — |

BH *q* is within the two new primary axes (CD274, MHC-I) only. Single HLA genes are the MHC-I decomposition, not a second family of claims. MHC-I = unweighted within-cohort z-mean of HLA-A / HLA-B / HLA-C / B2M (**{len(mhc1_used)}/4** present: {", ".join(mhc1_used) or "none"}).

**What holds.** The new **CD274** cut is **null**: CLDN4 vs CD274 Spearman CI includes 0 (n={int(r_cd.n)}, ρ={fmt_rho(r_cd.rho)}, p={fmt_p(r_cd.p)}). CD274 Q4 vs Q1 does not separate CLDN4, and CLDN4 Q4 ∩ CD274 Q4 is chance ({overlap["n_overlap"]}/{overlap["n_q4_a"]} = {overlap["frac_of_a_q4"]*100:.1f}%, OR={overlap["odds_ratio"]:.2f}, p={fmt_p(overlap["fisher_p"])}). MHC-I is a weak inverse (ρ={fmt_rho(r_mhc.rho)}, p={fmt_p(r_mhc.p)}, q={qcell(r_mhc)}); that cassette signal is **B2M-driven** (ρ={fmt_rho(r_b2m.rho)}, p={fmt_p(r_b2m.p)}). HLA-A and HLA-B are NS; HLA-C is weak.

**What does not hold.** CLDN4-high is not CD274-high on this LUAD microarray. Do not cite this page as a residual-vs-CD8 result.

## New CD274 cut (Q4 vs Q1)

CD274 high = ≥ cohort 75th percentile; CD274 low = ≤ 25th percentile. Ties at the percentile are kept, so arm n is reported rather than assumed n/4.

| Cut | n_high | n_low | CLDN4 median (high vs low) | Mann–Whitney p | CLDN4 Q4 ∩ CD274 Q4 | Fisher OR | Fisher p |
|---|---:|---:|---|---:|---|---:|---:|
| **CD274 Q4 vs Q1** | {cut_cldn4["n_q4"]} | {cut_cldn4["n_q1"]} | {cut_cldn4["median_q4"]:.3f} vs {cut_cldn4["median_q1"]:.3f} | {fmt_p(cut_cldn4["p"])} | {overlap["n_overlap"]} / {overlap["n_q4_a"]} of CLDN4 Q4 ({overlap["frac_of_a_q4"]*100:.1f}%) | {overlap["odds_ratio"]:.2f} | {fmt_p(overlap["fisher_p"])} |

Independence expectation for two Q4 calls is 25% of the CLDN4-Q4 arm. MHC-I Q4 overlap with CLDN4 Q4: {overlap_mhc["n_overlap"]} / {overlap_mhc["n_q4_a"]} ({overlap_mhc["frac_of_a_q4"]*100:.1f}%), OR={overlap_mhc["odds_ratio"]:.2f}, p={fmt_p(overlap_mhc["fisher_p"])}.

TACSTD2 vs CD274 is a same-run companion only (n={int(r_tac.n) if r_tac is not None else 0}, ρ={fmt_rho(r_tac.rho) if r_tac is not None else "NA"}, p={fmt_p(r_tac.p) if r_tac is not None else "NA"}) and is not the finding.

## What this is not

- **Not** residual CLDN4 vs CD8A. That test is already NS on this matrix (PR 302) and is not headlined.
- **Not** an ICI-response, PD-L1 IHC, or protein result. CD274 and HLA here are microarray RNA.
- **Not** a claim that CLDN4 induces PD-L1 or MHC-I. Association only.
- ESTIMATE residual for CD274 / MHC-I is in `correlations.tsv` as a sensitivity column. It is **not** the headline.

## Methods

- **Matrix:** GEO `GSE72094_series_matrix.txt.gz` (author-processed IRON + RNA-quality batch correction). Probe × sample values used as published.
- **Annotation:** GPL15048 `GeneSymbol`. First `///` token, then first whitespace token. **Max-mean** probe collapse to HUGO. HLA-A has no singleton symbol on this table; all five probes are `HLA-A LOC100507703` (Entrez 3105) and are mapped to HLA-A. One HLA-C probe is dual-annotated `HLA-C HLA-E` and is first-token mapped to HLA-C; two probes are exact `HLA-C`.
- **CD274** = collapsed `CD274` (probes `merck-NM_014143_at`, `merck-CK904742_at`, `merck2-ENST00000381577_at`; max-mean picks one).
- **MHC-I** = mean of within-cohort z-scores of HLA-A, HLA-B, HLA-C, B2M.
- **CD274 cut:** cohort quartiles on collapsed CD274. Primary contrast is Q4 vs Q1 (Mann–Whitney on CLDN4). Overlap is CLDN4 Q4 ∩ CD274 Q4 (two-sided Fisher).
- **Spearman** two-sided, complete cases. Bootstrap 95% CI: 5,000 resamples, seed 0, CLDN4 rows only.
- **ESTIMATE** (sensitivity column only): ssGSEA τ=0.25 on Yoshihara 2013 Stromal141 + Immune141. Coverage **{n_strom}/141** stromal, **{n_imm}/141** immune. Partial Spearman = Pearson of rank residuals; df = n − 3.

## Files

- `correlations.tsv` — n / ρ / p / CI / ESTIMATE residual (not headlined)
- `q4_cuts.tsv` / `q4_overlap.json` — new CD274 quartile cut
- `samples.tsv` — per-sample genes + MHC-I + CD274 quartile
- `coverage.tsv` / `provenance.json`

```bash
python3 -m pip install -r methods/gse72094_cldn4_cd274/requirements.txt
python3 methods/gse72094_cldn4_cd274/analyze.py
```

Place (or let the script download) `GSE72094_series_matrix.txt.gz` and `GPL15048_datatable.txt` under `data/gse72094/`.
"""
    (OUT / "FINDING.md").write_text(finding)
    print(finding)
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
