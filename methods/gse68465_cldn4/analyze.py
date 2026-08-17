#!/usr/bin/env python3
"""GSE68465 LUAD: CLDN4 vs CD8A after ESTIMATE TumorPurity.

Additive CLDN4 slice. TACSTD2 partial ρ = −0.17, p = 3e-4 is taken as given
from PR #235 and is not recomputed here.

Public GEO series matrix + GPL96 annotation + Yoshihara 2013 ESTIMATE lists.
Honest n / ρ / p only.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
DATA = ROOT / "data" / "gse68465"
SIG = ROOT / "data" / "signatures"
OUT = HERE
TABLES = OUT / "tables"
FIG = OUT / "figures"
TABLES.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

# Yoshihara 2013 ESTIMATE TumorPurity transform
EST_A = 0.6049872018
EST_B = 0.0001467884
BOOT_SEED = 20260817
N_BOOT = 2000

# Named Affy U133A probes (official GEO GPL96 annotation)
NAMED_PROBES = {
    "CLDN4": "201428_at",
    "CD8A": "205758_at",
}


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


def bootstrap_spearman(x, y, n_boot=N_BOOT, seed=BOOT_SEED):
    m = np.isfinite(x) & np.isfinite(y)
    xx, yy = x[m], y[m]
    n = int(len(xx))
    rng = np.random.default_rng(seed)
    rhos = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        r, _ = stats.spearmanr(xx[idx], yy[idx])
        rhos[i] = r
    lo, hi = np.quantile(rhos, [0.025, 0.975])
    return float(lo), float(hi)


def bootstrap_partial(x, y, z, n_boot=N_BOOT, seed=BOOT_SEED):
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    xx, yy, zz = x[m], y[m], z[m]
    n = int(len(xx))
    rng = np.random.default_rng(seed)
    rhos = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        r, _, _ = partial_spearman(xx[idx], yy[idx], zz[idx])
        rhos[i] = r
    lo, hi = np.quantile(rhos, [0.025, 0.975])
    return float(lo), float(hi)


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
    expr_start = None
    with gzip.open(path, "rt", errors="replace") as f:
        lines = f.readlines()
    samples = None
    for i, line in enumerate(lines):
        if line.startswith("!Sample_geo_accession"):
            samples = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
        if line.startswith("!Sample_title"):
            vals = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            meta_rows["title"] = vals
        if line.startswith("!Sample_source_name"):
            vals = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            meta_rows["source"] = vals
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


def load_gpl_annot(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt", errors="replace") as f:
        skip = 0
        for i, line in enumerate(f):
            if line.startswith("ID\t") or line.startswith("ID "):
                skip = i
                break
    df = pd.read_csv(path, sep="\t", skiprows=skip, dtype=str, low_memory=False)
    return df


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


def rec_assoc(predictor, endpoint, x, y, z=None, purity_name=""):
    ru, pu, nu = spearman_pair(x, y)
    lo, hi = bootstrap_spearman(x, y)
    if z is None:
        rp, pp, np_ = np.nan, np.nan, nu
        plo, phi = np.nan, np.nan
    else:
        rp, pp, np_ = partial_spearman(x, y, z)
        plo, phi = bootstrap_partial(x, y, z)
    return {
        "predictor": predictor,
        "endpoint": endpoint,
        "n": nu,
        "rho": ru,
        "rho_ci_lo": lo,
        "rho_ci_hi": hi,
        "p": pu,
        "n_partial": np_,
        "partial_rho": rp,
        "partial_ci_lo": plo,
        "partial_ci_hi": phi,
        "partial_p": pp,
        "purity": purity_name,
    }


def main():
    matrix = DATA / "GSE68465_series_matrix.txt.gz"
    annot = DATA / "GPL96.annot.gz"
    est_path = SIG / "estimate_yoshihara_2013.tsv"
    if not matrix.exists() or not annot.exists():
        raise SystemExit(
            "Missing downloads. See FINDING.md Reproduce. "
            f"Need {matrix} and {annot}."
        )

    est_tab = pd.read_csv(est_path, sep="\t")
    stromal = est_tab.loc[est_tab["set"] == "Stromal141_UP", "hugo"].tolist()
    immune = est_tab.loc[est_tab["set"] == "Immune141_UP", "hugo"].tolist()

    gpl = load_gpl_annot(annot)
    id_col = "ID" if "ID" in gpl.columns else gpl.columns[0]
    sym_col = "Gene symbol" if "Gene symbol" in gpl.columns else [c for c in gpl.columns if "symbol" in c.lower()][0]
    p2g = gpl.set_index(id_col)[sym_col].astype(str)
    p2g = p2g.replace({"nan": np.nan, "None": np.nan, "": np.nan}).dropna()
    p2g = p2g.map(lambda x: str(x).split("///")[0].strip())
    p2g = p2g[p2g.str.len() > 0]

    # platform confirmation for named probes
    probe_rows = []
    for gene, probe in NAMED_PROBES.items():
        hit = gpl[gpl[id_col] == probe]
        if hit.empty:
            probe_rows.append({"gene": gene, "probe": probe, "gpl96_symbol": "ABSENT", "in_matrix": False})
            continue
        probe_rows.append(
            {
                "gene": gene,
                "probe": probe,
                "gpl96_symbol": str(hit.iloc[0][sym_col]).split("///")[0].strip(),
                "entrez": str(hit.iloc[0].get("Gene ID", "")),
                "title": str(hit.iloc[0].get("Gene title", "")),
            }
        )
    probe_df = pd.DataFrame(probe_rows)

    print("Parsing GSE68465 series matrix…")
    meta, probes = parse_geo_matrix(matrix)
    print(f"  matrix samples={meta.shape[0]} probes={probes.shape[0]}")
    disease = meta["disease_state"].astype(str)
    keep = disease.str.contains("Adenocarcinoma", case=False, na=False)
    print(f"  disease_state counts:\n{disease.value_counts().to_string()}")
    meta_t = meta.loc[keep]
    probes_t = probes.loc[:, meta_t.index]
    n_tumor = int(keep.sum())
    print(f"  LUAD tumors n={n_tumor}")

    for rec in probe_rows:
        rec["in_matrix"] = rec["probe"] in probes.index
    probe_df = pd.DataFrame(probe_rows)
    probe_df.to_csv(TABLES / "platform_probe_confirm.tsv", sep="\t", index=False)

    genes = collapse_maxmean(probes_t, p2g)
    print(f"  genes after max-mean collapse: {genes.shape[0]}")
    for g in ("CLDN4", "CD8A"):
        if g not in genes.index:
            raise SystemExit(f"{g} missing after collapse")

    # which probe won max-mean
    used_probes = {}
    for g in ("CLDN4", "CD8A"):
        # invert: find probe whose values match the collapsed gene
        gvals = genes.loc[g]
        cands = p2g[p2g == g].index.intersection(probes_t.index)
        for pr in cands:
            if np.allclose(probes_t.loc[pr].to_numpy(float), gvals.to_numpy(float), equal_nan=True):
                used_probes[g] = pr
                break
        print(f"  {g} max-mean probe={used_probes.get(g, 'NA')} named={NAMED_PROBES[g]}")

    est = estimate_purity(genes, stromal, immune)
    strom_n = sum(g in genes.index for g in stromal)
    imm_n = sum(g in genes.index for g in immune)
    print(f"  ESTIMATE coverage stromal={strom_n}/141 immune={imm_n}/141")

    cldn4 = genes.loc["CLDN4"].to_numpy(float)
    cd8a = genes.loc["CD8A"].to_numpy(float)
    pur = est["ESTIMATE_TumorPurity"].to_numpy(float)
    imm = est["ESTIMATE_ImmuneScore"].to_numpy(float)
    strom = est["ESTIMATE_StromalScore"].to_numpy(float)
    estsc = est["ESTIMATE_Score"].to_numpy(float)

    purity_name = "ESTIMATE TumorPurity (ssGSEA Yoshihara 141+141)"
    rows = [
        rec_assoc("CLDN4", "CD8A", cldn4, cd8a, pur, purity_name),
        rec_assoc("CLDN4", "ESTIMATE_ImmuneScore", cldn4, imm, pur, purity_name),
        rec_assoc("CLDN4", "ESTIMATE_StromalScore", cldn4, strom, pur, purity_name),
        rec_assoc("CLDN4", "ESTIMATE_TumorPurity", cldn4, pur, None, ""),
        rec_assoc("CD8A", "ESTIMATE_TumorPurity", cd8a, pur, None, ""),
        rec_assoc("CD8A", "ESTIMATE_ImmuneScore", cd8a, imm, None, ""),
        # sensitivity: residualize on ImmuneScore instead of TumorPurity
        rec_assoc("CLDN4", "CD8A", cldn4, cd8a, imm, "ESTIMATE ImmuneScore (sensitivity residual)"),
        rec_assoc("CLDN4", "CD8A", cldn4, cd8a, estsc, "ESTIMATE Score (sensitivity residual)"),
    ]
    assoc = pd.DataFrame(rows)
    assoc.to_csv(TABLES / "associations.tsv", sep="\t", index=False)

    # named-probe sensitivity (no collapse)
    named_rows = []
    if NAMED_PROBES["CLDN4"] in probes_t.index and NAMED_PROBES["CD8A"] in probes_t.index:
        nx = probes_t.loc[NAMED_PROBES["CLDN4"]].to_numpy(float)
        ny = probes_t.loc[NAMED_PROBES["CD8A"]].to_numpy(float)
        named_rows.append(rec_assoc("CLDN4 named 201428_at", "CD8A named 205758_at", nx, ny, pur, purity_name))
    named = pd.DataFrame(named_rows)
    if not named.empty:
        named.to_csv(TABLES / "named_probe_sensitivity.tsv", sep="\t", index=False)

    st = pd.DataFrame(
        {
            "CLDN4": genes.loc["CLDN4"],
            "CD8A": genes.loc["CD8A"],
            "CLDN4_probe": used_probes.get("CLDN4", ""),
            "CD8A_probe": used_probes.get("CD8A", ""),
        },
        index=genes.columns,
    )
    st = st.join(est)
    for col in ("sex", "age", "disease_stage", "histologic_grade", "smoking_history"):
        if col in meta_t.columns:
            st[col] = meta_t[col]
    st.index.name = "sample"
    st.to_csv(TABLES / "sample_scores.tsv", sep="\t")

    cov = pd.DataFrame(
        [
            {
                "cohort": "GSE68465",
                "title": "Director's Challenge LUAD (Jacob / Shedden multi-site)",
                "platform": "GPL96 U133A",
                "n_matrix": int(meta.shape[0]),
                "n_normal": int((~keep).sum()),
                "n_tumor": n_tumor,
                "n_genes_collapsed": int(genes.shape[0]),
                "cldn4_present": True,
                "cd8a_present": True,
                "cldn4_maxmean_probe": used_probes.get("CLDN4", ""),
                "cd8a_maxmean_probe": used_probes.get("CD8A", ""),
                "estimate_stromal_present": strom_n,
                "estimate_immune_present": imm_n,
                "purity": purity_name,
            }
        ]
    )
    cov.to_csv(TABLES / "coverage.tsv", sep="\t", index=False)

    # figures
    fig, axes = plt.subplots(1, 3, figsize=(12.2, 3.8))
    axes[0].scatter(cldn4, cd8a, s=10, alpha=0.45, c="#1f4e79", edgecolors="none")
    axes[0].set_xlabel("CLDN4 (max-mean, log2)")
    axes[0].set_ylabel("CD8A (max-mean, log2)")
    r0 = assoc.loc[0]
    axes[0].set_title(f"Unadj ρ={r0['rho']:+.3f} p={fmt_p(r0['p'])}\nn={int(r0['n'])}")

    xr = residualize(stats.rankdata(cldn4), stats.rankdata(pur).reshape(-1, 1))
    yr = residualize(stats.rankdata(cd8a), stats.rankdata(pur).reshape(-1, 1))
    axes[1].scatter(xr, yr, s=10, alpha=0.45, c="#922b21", edgecolors="none")
    axes[1].axhline(0, color="k", lw=0.6)
    axes[1].axvline(0, color="k", lw=0.6)
    axes[1].set_xlabel("CLDN4 rank residual | purity")
    axes[1].set_ylabel("CD8A rank residual | purity")
    axes[1].set_title(f"Partial ρ={r0['partial_rho']:+.3f} p={fmt_p(r0['partial_p'])}\nn={int(r0['n_partial'])}")

    axes[2].scatter(pur, cldn4, s=10, alpha=0.45, c="#1e8449", edgecolors="none")
    axes[2].set_xlabel("ESTIMATE TumorPurity")
    axes[2].set_ylabel("CLDN4")
    rpur = assoc.loc[3]
    axes[2].set_title(f"CLDN4 vs purity ρ={rpur['rho']:+.3f}\np={fmt_p(rpur['p'])} n={int(rpur['n'])}")
    fig.suptitle("GSE68465 LUAD n=443 — CLDN4 vs CD8A after ESTIMATE purity", fontsize=11)
    fig.tight_layout()
    fig.savefig(FIG / "fig1_cldn4_cd8a_purity.png", dpi=150)
    fig.savefig(FIG / "fig1_cldn4_cd8a_purity.pdf")
    plt.close(fig)

    summary = {
        "task": "gse68465_cldn4",
        "note": "TACSTD2 from PR #235 taken as given; not recomputed",
        "tacstd2_given_pr235": {
            "n": 443,
            "partial_rho": -0.17,
            "partial_p": 3e-4,
            "source": "https://github.com/jinxuanhong1-blip/sdaxcge/pull/235",
        },
        "n_tumor": n_tumor,
        "cldn4_vs_cd8a": assoc.iloc[0].to_dict(),
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "bootstrap": {"n": N_BOOT, "seed": BOOT_SEED},
        "estimate_coverage": {"stromal": strom_n, "immune": imm_n},
        "used_probes": used_probes,
    }
    with open(TABLES / "summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=float)

    prov = {
        "inputs": {},
        "signatures": str(est_path),
        "methods": {
            "collapse": "max-mean probe to HUGO",
            "ssgsea": "Barbie/GSVA tau=0.25, ranks 1..10000",
            "purity": "TumorPurity = cos(0.6049872018 + 0.0001467884 × ESTIMATEScore)",
            "partial": "Pearson of rank residuals; df = n − 3",
        },
    }
    for fn in (matrix, annot, est_path):
        if fn.exists():
            prov["inputs"][fn.name] = {"bytes": fn.stat().st_size, "md5": md5(fn)}
    with open(TABLES / "provenance.json", "w") as f:
        json.dump(prov, f, indent=2)

    print("\n=== CLDN4 associations ===")
    print(
        assoc[["predictor", "endpoint", "n", "rho", "p", "partial_rho", "partial_p"]].to_string(
            index=False, float_format=lambda v: f"{v:.4g}"
        )
    )
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
