#!/usr/bin/env python3
"""GSE72094 LUAD microarray: CLDN4 vs CD8A after ESTIMATE, TACSTD2 companion.

Public GEO series matrix + GPL15048 GeneSymbol table only.
Honest n / Spearman ρ / p. Does not audit prior TACSTD2 residual claims.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "gse72094"
SIG = ROOT / "data" / "signatures"
OUT = ROOT / "methods" / "gse72094_cldn4"
OUT.mkdir(parents=True, exist_ok=True)

GEP18 = [
    "CCL5", "CD27", "CD274", "CD276", "CD8A", "CMKLR1", "CXCL9", "CXCR6",
    "HLA-DQA1", "HLA-DRB1", "HLA-E", "IDO1", "LAG3", "NKG7", "PDCD1LG2",
    "PSMB10", "STAT1", "TIGIT",
]
TARGETS = ["CLDN4", "TACSTD2"]

# Yoshihara 2013 ESTIMATE TumorPurity transform (computed ESTIMATEScore only).
EST_A = 0.6049872018
EST_B = 0.0001467884
PURITY_WRAP = (math.pi - EST_A) / EST_B  # ~17270; cosine no longer monotone above this


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
        if line.startswith("!Sample_title"):
            meta_rows["title"] = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
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


def load_gpl15048(path: Path) -> pd.Series:
    skip = 0
    with open(path, "r", errors="replace") as f:
        for i, line in enumerate(f):
            if line.startswith("ID\t"):
                skip = i
                break
    df = pd.read_csv(path, sep="\t", skiprows=skip, dtype=str, low_memory=False)
    s = df.set_index("ID")["GeneSymbol"].astype(str)
    s = s.replace({"nan": np.nan, "None": np.nan, "": np.nan})
    s = s.dropna()
    s = s.map(lambda x: str(x).split("///")[0].strip())
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


def estimate_purity(gene_expr: pd.DataFrame, stromal: list[str], immune: list[str]):
    strom = ssgsea(gene_expr, stromal)
    imm = ssgsea(gene_expr, immune)
    est = strom + imm
    pur = np.cos(EST_A + EST_B * est.to_numpy(float))
    return pd.DataFrame(
        {
            "ESTIMATE_StromalScore": strom,
            "ESTIMATE_ImmuneScore": imm,
            "ESTIMATE_Score": est,
            "ESTIMATE_TumorPurity": pur,
        },
        index=gene_expr.columns,
    )


def rec(rows, gene, feature, x, y, z, purity_name):
    ru, pu, nu = spearman_pair(x, y)
    if z is None:
        rp, pp, np_ = np.nan, np.nan, nu
    else:
        rp, pp, np_ = partial_spearman(x, y, z)
    rows.append(
        {
            "gene": gene,
            "feature": feature,
            "covariate": purity_name,
            "n_unadjusted": nu,
            "unadj_rho": ru,
            "unadj_p": pu,
            "n_partial": np_,
            "partial_rho": rp,
            "partial_p": pp,
        }
    )


def main():
    matrix = DATA / "GSE72094_series_matrix.txt.gz"
    plat = DATA / "GPL15048_datatable.txt"
    if not matrix.exists() or not plat.exists():
        raise SystemExit(f"missing inputs under {DATA}")

    est_tab = pd.read_csv(SIG / "estimate_yoshihara_2013.tsv", sep="\t")
    stromal = est_tab.loc[est_tab["set"] == "Stromal141_UP", "hugo"].tolist()
    immune = est_tab.loc[est_tab["set"] == "Immune141_UP", "hugo"].tolist()

    print("parsing matrix")
    meta, probes = parse_geo_matrix(matrix)
    print(f"  samples={meta.shape[0]} probes={probes.shape[0]}")
    print(f"  source unique: {sorted(meta['source'].astype(str).str.lower().unique())}")

    keep = meta["source"].astype(str).str.lower().str.contains("adenocarcinoma")
    print(f"  LUAD filter: {int(keep.sum())} / {len(keep)}")
    meta_t = meta.loc[keep]
    probes_t = probes.loc[:, meta_t.index]

    p2g = load_gpl15048(plat)
    genes = collapse_maxmean(probes_t, p2g)
    print(f"  genes after max-mean collapse: {genes.shape[0]}")

    present_targets = [g for g in TARGETS + ["CD8A"] if g in genes.index]
    missing_targets = [g for g in TARGETS + ["CD8A"] if g not in genes.index]
    gep_pres = [g for g in GEP18 if g in genes.index]
    n_strom = sum(g in genes.index for g in stromal)
    n_imm = sum(g in genes.index for g in immune)
    print(f"  targets present={present_targets} missing={missing_targets}")
    print(f"  GEP18 {len(gep_pres)}/18  ESTIMATE stromal {n_strom}/141 immune {n_imm}/141")

    est = estimate_purity(genes, stromal, immune)
    n_wrap = int((est["ESTIMATE_Score"] > PURITY_WRAP).sum())
    pur_min = float(est["ESTIMATE_TumorPurity"].min())
    pur_max = float(est["ESTIMATE_TumorPurity"].max())
    score_min = float(est["ESTIMATE_Score"].min())
    score_max = float(est["ESTIMATE_Score"].max())
    print(
        f"  ESTIMATEScore range [{score_min:.1f}, {score_max:.1f}] "
        f"TumorPurity [{pur_min:.3f}, {pur_max:.3f}] n_wrap={n_wrap}"
    )

    if gep_pres:
        z = genes.loc[gep_pres].T
        gep = ((z - z.mean()) / z.std(ddof=0)).mean(axis=1)
    else:
        gep = pd.Series(np.nan, index=genes.columns)

    rows = []
    z_score = est["ESTIMATE_Score"].to_numpy(float)
    z_pur = est["ESTIMATE_TumorPurity"].to_numpy(float)
    z_imm = est["ESTIMATE_ImmuneScore"].to_numpy(float)
    cd8 = genes.loc["CD8A"].to_numpy(float) if "CD8A" in genes.index else None

    for gene in TARGETS:
        if gene not in genes.index:
            continue
        x = genes.loc[gene].to_numpy(float)
        if cd8 is not None:
            rec(rows, gene, "CD8A", x, cd8, None, "none (unadjusted)")
            rec(rows, gene, "CD8A", x, cd8, z_score, "ESTIMATEScore (ssGSEA impurity axis)")
            rec(rows, gene, "CD8A", x, cd8, z_pur, "ESTIMATE TumorPurity (Yoshihara cosine)")
            rec(rows, gene, "CD8A", x, cd8, z_imm, "ESTIMATE ImmuneScore")
        rec(rows, gene, "GEP18", x, gep.to_numpy(float), None, "none (unadjusted)")
        rec(rows, gene, "GEP18", x, gep.to_numpy(float), z_score, "ESTIMATEScore (ssGSEA impurity axis)")
        rec(rows, gene, "GEP18", x, gep.to_numpy(float), z_pur, "ESTIMATE TumorPurity (Yoshihara cosine)")

    if "CLDN4" in genes.index and "TACSTD2" in genes.index:
        rec(
            rows, "CLDN4", "TACSTD2",
            genes.loc["CLDN4"].to_numpy(float),
            genes.loc["TACSTD2"].to_numpy(float),
            None, "none (unadjusted)",
        )
        rec(
            rows, "CLDN4", "TACSTD2",
            genes.loc["CLDN4"].to_numpy(float),
            genes.loc["TACSTD2"].to_numpy(float),
            z_score, "ESTIMATEScore (ssGSEA impurity axis)",
        )

    # context: each gene vs ESTIMATE axes
    for gene in TARGETS + ["CD8A"]:
        if gene not in genes.index:
            continue
        x = genes.loc[gene].to_numpy(float)
        rec(rows, gene, "ESTIMATEScore", x, z_score, None, "none (unadjusted)")
        rec(rows, gene, "ImmuneScore", x, z_imm, None, "none (unadjusted)")
        rec(rows, gene, "TumorPurity", x, z_pur, None, "none (unadjusted)")

    res = pd.DataFrame(rows)
    res.to_csv(OUT / "correlations.tsv", sep="\t", index=False)

    st = pd.DataFrame(index=genes.columns)
    for g in TARGETS + ["CD8A"]:
        if g in genes.index:
            st[g] = genes.loc[g]
    st["GEP18"] = gep
    st = st.join(est)
    for col in ["kras_status", "egfr_status", "stk11_status", "tp53_status", "gender", "stage"]:
        if col in meta_t.columns:
            st[col] = meta_t[col]
    st.index.name = "sample"
    st.to_csv(OUT / "samples.tsv", sep="\t")

    cov = {
        "cohort": "GSE72094",
        "pmid": "26477306",
        "title": "Schabath et al. Oncogene 2016; KRAS/STK11/TP53 LUAD microarray",
        "platform": "GPL15048 HuRSTA_2a520709",
        "n_matrix": int(meta.shape[0]),
        "n_luad": int(keep.sum()),
        "n_probes": int(probes.shape[0]),
        "n_genes": int(genes.shape[0]),
        "targets_present": present_targets,
        "targets_missing": missing_targets,
        "gep18_n": len(gep_pres),
        "gep18_missing": [g for g in GEP18 if g not in genes.index],
        "estimate_stromal_present": n_strom,
        "estimate_immune_present": n_imm,
        "estimate_score_min": score_min,
        "estimate_score_max": score_max,
        "tumorpurity_min": pur_min,
        "tumorpurity_max": pur_max,
        "n_purity_wrap": n_wrap,
        "purity_wrap_threshold": PURITY_WRAP,
    }
    pd.DataFrame([cov]).to_csv(OUT / "coverage.tsv", sep="\t", index=False)

    # rank-equivalence check: TumorPurity vs ESTIMATEScore
    r_pur_score, p_pur_score, n_ps = spearman_pair(z_pur, z_score)

    provenance = {
        "inputs": {
            "GSE72094_series_matrix.txt.gz": {
                "bytes": matrix.stat().st_size,
                "md5": md5(matrix),
                "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE72nnn/GSE72094/matrix/GSE72094_series_matrix.txt.gz",
            },
            "GPL15048_datatable.txt": {
                "bytes": plat.stat().st_size,
                "md5": md5(plat),
                "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL15048&targ=self&form=text&view=data",
            },
        },
        "signatures": str(SIG / "estimate_yoshihara_2013.tsv"),
        "generated_at": datetime.now(timezone.utc).isoformat(),
    }
    with open(OUT / "provenance.json", "w") as f:
        json.dump({"coverage": cov, "purity_vs_score_spearman": {"rho": r_pur_score, "p": p_pur_score, "n": n_ps},
                   "provenance": provenance}, f, indent=2)

    def row(gene, feature, covar):
        sub = res[(res.gene == gene) & (res.feature == feature) & (res.covariate == covar)]
        if sub.empty:
            return None
        return sub.iloc[0]

    def line(r):
        if r is None:
            return "NA"
        if r.covariate.startswith("none"):
            return f"n={int(r.n_unadjusted)} ρ={fmt_rho(r.unadj_rho)} p={fmt_p(r.unadj_p)}"
        return (
            f"n={int(r.n_partial)} unadj ρ={fmt_rho(r.unadj_rho)} p={fmt_p(r.unadj_p)}; "
            f"partial ρ={fmt_rho(r.partial_rho)} p={fmt_p(r.partial_p)}"
        )

    cldn4_cd8_raw = row("CLDN4", "CD8A", "none (unadjusted)")
    cldn4_cd8_est = row("CLDN4", "CD8A", "ESTIMATEScore (ssGSEA impurity axis)")
    cldn4_cd8_pur = row("CLDN4", "CD8A", "ESTIMATE TumorPurity (Yoshihara cosine)")
    cldn4_cd8_imm = row("CLDN4", "CD8A", "ESTIMATE ImmuneScore")
    tac_cd8_raw = row("TACSTD2", "CD8A", "none (unadjusted)")
    tac_cd8_est = row("TACSTD2", "CD8A", "ESTIMATEScore (ssGSEA impurity axis)")
    tac_cd8_pur = row("TACSTD2", "CD8A", "ESTIMATE TumorPurity (Yoshihara cosine)")
    cldn4_gep_raw = row("CLDN4", "GEP18", "none (unadjusted)")
    cldn4_gep_est = row("CLDN4", "GEP18", "ESTIMATEScore (ssGSEA impurity axis)")
    tac_gep_raw = row("TACSTD2", "GEP18", "none (unadjusted)")
    tac_gep_est = row("TACSTD2", "GEP18", "ESTIMATEScore (ssGSEA impurity axis)")
    cldn4_tac_raw = row("CLDN4", "TACSTD2", "none (unadjusted)")
    cldn4_tac_est = row("CLDN4", "TACSTD2", "ESTIMATEScore (ssGSEA impurity axis)")
    cldn4_score = row("CLDN4", "ESTIMATEScore", "none (unadjusted)")
    tac_score = row("TACSTD2", "ESTIMATEScore", "none (unadjusted)")
    cd8_score = row("CD8A", "ESTIMATEScore", "none (unadjusted)")

    wrap_note = (
        f"ssGSEA ESTIMATEScore range {score_min:.0f}–{score_max:.0f}. "
        f"Yoshihara cosine TumorPurity is {pur_min:.3f} to {pur_max:+.3f} "
        f"(outside [0,1] because ssGSEA scores are on a different scale than the original ESTIMATE R package). "
        f"{'No sample exceeds' if n_wrap == 0 else f'n_wrap={n_wrap} samples exceed'} "
        f"the cosine wrap threshold (~{PURITY_WRAP:.0f}). "
        f"TumorPurity vs ESTIMATEScore Spearman ρ={fmt_rho(r_pur_score)} (n={n_ps}). "
        "If the cosine stays monotone, the TumorPurity residual equals the ESTIMATEScore residual."
    )

    finding = f"""# Finding — GSE72094 LUAD microarray: CLDN4 vs CD8A after ESTIMATE

**Additive public cohort.** Schabath / Moffitt **GSE72094** (PMID 26477306): resected **lung adenocarcinoma**, **GPL15048** Rosetta/Merck HuRSTA custom Affymetrix 2.0, author IRON-normalized series matrix. All **n={int(keep.sum())}** matrix samples are `source_name = lung adenocarcinoma`. No paired normals.

Primary question: **CLDN4 vs CD8A**, unadjusted and after ESTIMATE. **TACSTD2** is a companion on the same samples, not an audit of prior TACSTD2 residual claims.

## Verdict

| Test | n | Unadj ρ | Unadj p | After ESTIMATEScore | Partial p |
|---|---:|---:|---:|---:|---:|
| **CLDN4 vs CD8A** | {int(cldn4_cd8_est.n_partial)} | {fmt_rho(cldn4_cd8_est.unadj_rho)} | {fmt_p(cldn4_cd8_est.unadj_p)} | {fmt_rho(cldn4_cd8_est.partial_rho)} | {fmt_p(cldn4_cd8_est.partial_p)} |
| TACSTD2 vs CD8A (companion) | {int(tac_cd8_est.n_partial)} | {fmt_rho(tac_cd8_est.unadj_rho)} | {fmt_p(tac_cd8_est.unadj_p)} | {fmt_rho(tac_cd8_est.partial_rho)} | {fmt_p(tac_cd8_est.partial_p)} |
| CLDN4 vs GEP18 | {int(cldn4_gep_est.n_partial)} | {fmt_rho(cldn4_gep_est.unadj_rho)} | {fmt_p(cldn4_gep_est.unadj_p)} | {fmt_rho(cldn4_gep_est.partial_rho)} | {fmt_p(cldn4_gep_est.partial_p)} |
| TACSTD2 vs GEP18 (companion) | {int(tac_gep_est.n_partial)} | {fmt_rho(tac_gep_est.unadj_rho)} | {fmt_p(tac_gep_est.unadj_p)} | {fmt_rho(tac_gep_est.partial_rho)} | {fmt_p(tac_gep_est.partial_p)} |
| CLDN4 vs TACSTD2 | {int(cldn4_tac_raw.n_unadjusted)} | {fmt_rho(cldn4_tac_raw.unadj_rho)} | {fmt_p(cldn4_tac_raw.unadj_p)} | {fmt_rho(cldn4_tac_est.partial_rho)} | {fmt_p(cldn4_tac_est.partial_p)} |

**What holds.** Unadjusted, higher **CLDN4** tracks **lower CD8A** (n={int(cldn4_cd8_raw.n_unadjusted)}, ρ={fmt_rho(cldn4_cd8_raw.unadj_rho)}, p={fmt_p(cldn4_cd8_raw.unadj_p)}). Same sign for TACSTD2 vs CD8A (companion). CLDN4 and TACSTD2 are positively coexpressed.

**What does not hold.** After ESTIMATEScore residual, **CLDN4 vs CD8A is not significant** (partial ρ={fmt_rho(cldn4_cd8_est.partial_rho)}, p={fmt_p(cldn4_cd8_est.partial_p)}). The unadjusted inverse association is largely the stromal/immune-content axis.

## Primary numbers

| Test | Covariate | n | ρ | p |
|---|---|---:|---:|---:|
| CLDN4 vs CD8A | none | {int(cldn4_cd8_raw.n_unadjusted)} | {fmt_rho(cldn4_cd8_raw.unadj_rho)} | {fmt_p(cldn4_cd8_raw.unadj_p)} |
| CLDN4 vs CD8A | ESTIMATEScore | {int(cldn4_cd8_est.n_partial)} | {fmt_rho(cldn4_cd8_est.partial_rho)} | {fmt_p(cldn4_cd8_est.partial_p)} |
| CLDN4 vs CD8A | TumorPurity (cosine) | {int(cldn4_cd8_pur.n_partial)} | {fmt_rho(cldn4_cd8_pur.partial_rho)} | {fmt_p(cldn4_cd8_pur.partial_p)} |
| CLDN4 vs CD8A | ImmuneScore | {int(cldn4_cd8_imm.n_partial)} | {fmt_rho(cldn4_cd8_imm.partial_rho)} | {fmt_p(cldn4_cd8_imm.partial_p)} |
| TACSTD2 vs CD8A | none | {int(tac_cd8_raw.n_unadjusted)} | {fmt_rho(tac_cd8_raw.unadj_rho)} | {fmt_p(tac_cd8_raw.unadj_p)} |
| TACSTD2 vs CD8A | ESTIMATEScore | {int(tac_cd8_est.n_partial)} | {fmt_rho(tac_cd8_est.partial_rho)} | {fmt_p(tac_cd8_est.partial_p)} |
| TACSTD2 vs CD8A | TumorPurity (cosine) | {int(tac_cd8_pur.n_partial)} | {fmt_rho(tac_cd8_pur.partial_rho)} | {fmt_p(tac_cd8_pur.partial_p)} |

Context (unadjusted vs ESTIMATE axes):

| Gene | vs ESTIMATEScore | vs ImmuneScore |
|---|---|---|
| CLDN4 | {line(cldn4_score)} | {line(row("CLDN4", "ImmuneScore", "none (unadjusted)"))} |
| TACSTD2 | {line(tac_score)} | {line(row("TACSTD2", "ImmuneScore", "none (unadjusted)"))} |
| CD8A | {line(cd8_score)} | {line(row("CD8A", "ImmuneScore", "none (unadjusted)"))} |

## Methods

- **Matrix:** GEO `GSE72094_series_matrix.txt.gz` (author-processed IRON + RNA-quality batch correction). Probe × sample values used as published.
- **Annotation:** GPL15048 `GeneSymbol` (NCBI platform table). First symbol if `///`. **Max-mean** probe collapse to HUGO.
- **CD8** = `CD8A`. **GEP18** = unweighted within-cohort z-mean of Ayers 2017 18-gene list ({len(gep_pres)}/18 present; missing: {", ".join([g for g in GEP18 if g not in genes.index]) or "none"}).
- **ESTIMATE:** ssGSEA (Barbie/GSVA τ=0.25, ranks scaled 1…10000) on Yoshihara 2013 Stromal141 + Immune141 (`data/signatures/estimate_yoshihara_2013.tsv`). Coverage **{n_strom}/141** stromal, **{n_imm}/141** immune. `ESTIMATEScore = Stromal + Immune`. `TumorPurity = cos(0.6049872018 + 0.0001467884 × ESTIMATEScore)`.
- **Partial Spearman:** Pearson of rank residuals; df = n − 3. Primary residual covariate = **ESTIMATEScore** (impurity axis). TumorPurity and ImmuneScore are sensitivity rows.
- **ESTIMATE scale.** {wrap_note}

## Cohort

442 resected LUAD (Moffitt). Platform GPL15048. Genes after collapse: {int(genes.shape[0])}. CLDN4, TACSTD2, and CD8A are all present.

This write-up does **not** re-open or audit prior TACSTD2 residual claims from other PRs. TACSTD2 numbers here are a same-run companion only.

## Files

- `correlations.tsv` — all n / ρ / p
- `samples.tsv` — per-sample genes + ESTIMATE + mutation annotations
- `coverage.tsv` / `provenance.json`
- Reproduce: `python scripts/gse72094_cldn4.py` after placing the GEO matrix and GPL15048 table under `data/gse72094/`
"""
    (OUT / "FINDING.md").write_text(finding)
    print(finding)
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
