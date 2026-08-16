#!/usr/bin/env python3
"""B3 extra: CPTAC LUAD (and remaining open freeze proteomes) CLDN4/TJ protein
vs ImmuneScore/CD8/ESTIMATE, plus TCGA-LUSC TJ vs CD8/GEP.

B3 TCGA-LUAD RNA (TJ-high → CD8/GEP low) is taken as given and is not re-run.
CPTAC LSCC CLDN4 protein vs ImmuneScore is already public and is not re-audited.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from genes import (  # noqa: E402
    CD8,
    CYT,
    ENSEMBL,
    GEP18,
    PHENO_ALIASES,
    REMAINING_PROTEOME,
    SYNONYMS,
    TJ_15,
    TJ_7,
)

REPO = HERE.parents[1]
RNG_SEED = 20260816
N_BOOT = 2000
MIN_TJ_GENES = 5
MIN_TJ7_GENES = 4


def load_matrix(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", index_col=0)


def find_row(index: pd.Index, symbol: str) -> str | None:
    """Match a gene in an Ensembl or symbol index."""
    names = [symbol] + SYNONYMS.get(symbol, [])
    for name in names:
        if name in index:
            return str(name)
    for ens in ENSEMBL.get(symbol, []):
        hits = [i for i in index if str(i) == ens or str(i).startswith(ens + ".")]
        if hits:
            return str(hits[0])
    # last resort: symbol after a pipe / trailing annotation
    for name in names:
        hits = [i for i in index if str(i).split("|")[0] == name or str(i).endswith("|" + name)]
        if hits:
            return str(hits[0])
    return None


def extract_gene(mat: pd.DataFrame, symbol: str) -> pd.Series:
    row = find_row(mat.index, symbol)
    if row is None:
        return pd.Series(np.nan, index=mat.columns, name=symbol)
    s = pd.to_numeric(mat.loc[row], errors="coerce")
    s.name = symbol
    return s


def mean_z(mat: pd.DataFrame, symbols: list[str], min_genes: int) -> tuple[pd.Series, pd.DataFrame]:
    rows = []
    used = []
    coverage = []
    for sym in symbols:
        rid = find_row(mat.index, sym)
        found = rid is not None
        n_obs = 0
        if found:
            s = pd.to_numeric(mat.loc[rid], errors="coerce")
            n_obs = int(s.notna().sum())
            if n_obs >= 8:
                sd = s.std(ddof=0)
                if sd and not np.isnan(sd):
                    rows.append((s - s.mean()) / sd)
                    used.append(sym)
        coverage.append({"gene": sym, "row": rid, "found": found, "n_obs": n_obs})
    cov = pd.DataFrame(coverage)
    if not rows:
        return pd.Series(np.nan, index=mat.columns), cov
    z = pd.concat(rows, axis=1)
    n = z.notna().sum(axis=1)
    score = z.mean(axis=1, skipna=True)
    score[n < min_genes] = np.nan
    score.name = "score"
    cov.attrs["used"] = used
    return score, cov


def spearman(x: pd.Series, y: pd.Series) -> dict:
    d = pd.concat([x, y], axis=1).dropna()
    n = int(len(d))
    if n < 6 or d.iloc[:, 0].nunique() < 2 or d.iloc[:, 1].nunique() < 2:
        return {"n": n, "rho": np.nan, "p": np.nan, "ci95_lo": np.nan, "ci95_hi": np.nan}
    rho, p = stats.spearmanr(d.iloc[:, 0], d.iloc[:, 1])
    rng = np.random.default_rng(RNG_SEED)
    boots = []
    xv = d.iloc[:, 0].to_numpy()
    yv = d.iloc[:, 1].to_numpy()
    for _ in range(N_BOOT):
        idx = rng.integers(0, n, n)
        try:
            boots.append(stats.spearmanr(xv[idx], yv[idx]).statistic)
        except Exception:
            continue
    lo, hi = np.nanpercentile(boots, [2.5, 97.5]) if boots else (np.nan, np.nan)
    return {
        "n": n,
        "rho": float(rho),
        "p": float(p),
        "ci95_lo": float(lo),
        "ci95_hi": float(hi),
    }


def partial_spearman(x: pd.Series, y: pd.Series, z: pd.Series) -> dict:
    d = pd.concat([x, y, z], axis=1).dropna()
    d.columns = ["x", "y", "z"]
    n = int(len(d))
    if n < 8 or d["z"].nunique() < 2:
        return {"n": n, "rho": np.nan, "p": np.nan, "ci95_lo": np.nan, "ci95_hi": np.nan}
    rx, ry, rz = d["x"].rank(), d["y"].rank(), d["z"].rank()
    bx = np.polyfit(rz, rx, 1)
    by = np.polyfit(rz, ry, 1)
    ex = rx - (bx[0] * rz + bx[1])
    ey = ry - (by[0] * rz + by[1])
    out = spearman(ex, ey)
    out["n"] = n
    return out


def high_low(score: pd.Series, outcome: pd.Series) -> dict:
    d = pd.concat([score, outcome], axis=1).dropna()
    d.columns = ["s", "y"]
    if len(d) < 8:
        return {"n_high": 0, "n_low": 0, "median_high": np.nan, "median_low": np.nan, "p": np.nan, "U": np.nan}
    med = d["s"].median()
    hi = d.loc[d["s"] > med, "y"]
    lo = d.loc[d["s"] <= med, "y"]
    if len(hi) < 3 or len(lo) < 3:
        return {
            "n_high": int(len(hi)),
            "n_low": int(len(lo)),
            "median_high": float(hi.median()) if len(hi) else np.nan,
            "median_low": float(lo.median()) if len(lo) else np.nan,
            "p": np.nan,
            "U": np.nan,
        }
    U, p = stats.mannwhitneyu(hi, lo, alternative="two-sided")
    return {
        "n_high": int(len(hi)),
        "n_low": int(len(lo)),
        "median_high": float(hi.median()),
        "median_low": float(lo.median()),
        "delta_median_high_minus_low": float(hi.median() - lo.median()),
        "U": float(U),
        "p": float(p),
    }


def fmt_p(p) -> str:
    if p is None or (isinstance(p, float) and (np.isnan(p))):
        return "NA"
    p = float(p)
    if p == 0:
        return "<1e-300"
    if p < 1e-4:
        return f"{p:.1e}"
    return f"{p:.3g}"


def fmt_rho(r) -> str:
    if r is None or (isinstance(r, float) and np.isnan(r)):
        return "NA"
    return f"{float(r):+.3f}"


def pick_col(df: pd.DataFrame, key: str) -> str | None:
    for name in PHENO_ALIASES.get(key, [key]):
        if name in df.columns:
            return name
    # fuzzy: case-insensitive exact
    lower = {c.lower(): c for c in df.columns}
    for name in PHENO_ALIASES.get(key, [key]):
        if name.lower() in lower:
            return lower[name.lower()]
    # contains
    tokens = {
        "ImmuneScore": ["immunescore", "immune_score"],
        "CD8_CIBERSORT": ["cibersort"] ,
        "CD8_xCell": ["xcell"],
        "WES_purity": ["wes_purity"],
    }
    return None


def pick_immune_cols(pheno: pd.DataFrame) -> dict[str, str]:
    found = {}
    for key in PHENO_ALIASES:
        col = pick_col(pheno, key)
        if col:
            found[key] = col
    # extra CD8 hunt if aliases missed
    if "CD8_CIBERSORT" not in found:
        for c in pheno.columns:
            cl = c.lower()
            if "cibersort" in cl and "cd8" in cl:
                found["CD8_CIBERSORT"] = c
                break
    if "CD8_xCell" not in found:
        for c in pheno.columns:
            cl = c.lower()
            if "xcell" in cl and "cd8" in cl:
                found["CD8_xCell"] = c
                break
    if "ImmuneScore" not in found:
        for c in pheno.columns:
            cl = c.lower()
            if "estimate" in cl and "immune" in cl:
                found["ImmuneScore"] = c
                break
    return found


def harmonize_ids(idx) -> pd.Index:
    return pd.Index([str(s).replace(".", "-") for s in idx])


def load_pheno(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    # first column is usually sample id
    idcol = df.columns[0]
    df = df.rename(columns={idcol: "sample"})
    df["sample"] = df["sample"].astype(str).str.replace(".", "-", regex=False)
    df = df.drop_duplicates("sample").set_index("sample")
    return df


def join_protein_pheno(protein: pd.DataFrame, pheno: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    protein = protein.copy()
    protein.columns = harmonize_ids(protein.columns)
    common = [s for s in protein.columns if s in pheno.index]
    return protein[common], pheno.loc[common]


def estimate_purity(estimate_score: pd.Series) -> pd.Series:
    arg = 0.6049872018 + 0.0001467884 * estimate_score
    purity = np.cos(arg)
    purity = purity.where((arg >= 0) & (arg <= np.pi))
    return purity


# ---------------------------------------------------------------------------
# CPTAC LUAD
# ---------------------------------------------------------------------------
def analyze_luad(data: Path, out: Path) -> dict:
    prot = load_matrix(
        data / "cptac/LUAD/LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt"
    )
    rna = load_matrix(data / "cptac/LUAD/LUAD_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt")
    pheno = load_pheno(data / "cptac/LUAD/LUAD_phenotype.txt")
    prot, pheno_m = join_protein_pheno(prot, pheno)
    rna.columns = harmonize_ids(rna.columns)
    rna = rna[[c for c in prot.columns if c in rna.columns]]

    cols = pick_immune_cols(pheno_m)
    cldn4_p = extract_gene(prot, "CLDN4")
    tj15_p, cov15_p = mean_z(prot, TJ_15, MIN_TJ_GENES)
    tj7_p, cov7_p = mean_z(prot, TJ_7, MIN_TJ7_GENES)
    tj15_r, cov15_r = mean_z(rna, TJ_15, MIN_TJ_GENES)
    tj7_r, _ = mean_z(rna, TJ_7, MIN_TJ7_GENES)
    cd8_r, _ = mean_z(rna, CD8, 2)
    cyt_r, _ = mean_z(rna, CYT, 2)
    gep_r, cov_gep = mean_z(rna, GEP18, 10)

    immune = {}
    for key, col in cols.items():
        immune[key] = pd.to_numeric(pheno_m[col], errors="coerce")

    scores = pd.DataFrame({
        "sample": prot.columns,
        "CLDN4_protein": cldn4_p.reindex(prot.columns).values,
        "TJ15_protein": tj15_p.reindex(prot.columns).values,
        "TJ7_protein": tj7_p.reindex(prot.columns).values,
        "TJ15_RNA": tj15_r.reindex(prot.columns).values,
        "TJ7_RNA": tj7_r.reindex(prot.columns).values,
        "CD8_RNA": cd8_r.reindex(prot.columns).values,
        "CYT_RNA": cyt_r.reindex(prot.columns).values,
        "GEP18_RNA": gep_r.reindex(prot.columns).values,
    }).set_index("sample")
    for key, s in immune.items():
        scores[key] = s.reindex(scores.index).values

    # associations
    predictors = {
        "CLDN4_protein": scores["CLDN4_protein"],
        "TJ15_protein": scores["TJ15_protein"],
        "TJ7_protein": scores["TJ7_protein"],
        "TJ15_RNA": scores["TJ15_RNA"],
    }
    endpoints = {}
    for name in ["ImmuneScore", "CD8_CIBERSORT", "CD8_xCell", "Immune_xCell", "StromalScore"]:
        if name in scores.columns:
            endpoints[name] = scores[name]
    endpoints["CD8_RNA"] = scores["CD8_RNA"]
    endpoints["GEP18_RNA"] = scores["GEP18_RNA"]
    endpoints["CYT_RNA"] = scores["CYT_RNA"]

    rows = []
    for pname, px in predictors.items():
        for ename, ey in endpoints.items():
            sp = spearman(px, ey)
            rec = {"cohort": "CPTAC_LUAD", "predictor": pname, "endpoint": ename, "test": "spearman", **sp}
            if "WES_purity" in scores.columns:
                rec.update({f"partial_{k}": v for k, v in partial_spearman(px, ey, scores["WES_purity"]).items()})
            rows.append(rec)

    hl_rows = []
    for pname in ["CLDN4_protein", "TJ15_protein", "TJ15_RNA"]:
        for ename in endpoints:
            rec = {"cohort": "CPTAC_LUAD", "predictor": pname, "endpoint": ename, **high_low(scores[pname], endpoints[ename])}
            hl_rows.append(rec)

    assoc = pd.DataFrame(rows)
    hl = pd.DataFrame(hl_rows)
    scores.to_csv(out / "tables/luad_sample_scores.tsv", sep="\t")
    cov15_p.to_csv(out / "tables/luad_tj15_protein_coverage.tsv", sep="\t", index=False)
    cov7_p.to_csv(out / "tables/luad_tj7_protein_coverage.tsv", sep="\t", index=False)
    cov15_r.to_csv(out / "tables/luad_tj15_rna_coverage.tsv", sep="\t", index=False)
    cov_gep.to_csv(out / "tables/luad_gep_rna_coverage.tsv", sep="\t", index=False)
    assoc.to_csv(out / "tables/luad_associations.tsv", sep="\t", index=False)
    hl.to_csv(out / "tables/luad_highlow.tsv", sep="\t", index=False)

    n_cldn4 = int(scores["CLDN4_protein"].notna().sum())
    n_tum = int(len(scores))
    summary = {
        "cohort": "CPTAC_LUAD",
        "citation": "Gillette et al. Cell 2020; CPTAC freeze v1.2",
        "n_tumors": n_tum,
        "n_CLDN4_protein": n_cldn4,
        "n_CLDN4_protein_NA": n_tum - n_cldn4,
        "tj15_protein_genes_used": cov15_p.loc[cov15_p["n_obs"] >= 8, "gene"].tolist(),
        "tj15_protein_n_genes": int((cov15_p["n_obs"] >= 8).sum()),
        "phenotype_columns": cols,
        "headline": {},
    }
    for pred, end in [
        ("CLDN4_protein", "ImmuneScore"),
        ("CLDN4_protein", "CD8_CIBERSORT"),
        ("CLDN4_protein", "CD8_xCell"),
        ("CLDN4_protein", "CD8_RNA"),
        ("CLDN4_protein", "GEP18_RNA"),
        ("TJ15_protein", "ImmuneScore"),
        ("TJ15_protein", "CD8_CIBERSORT"),
        ("TJ15_protein", "CD8_xCell"),
        ("TJ15_protein", "CD8_RNA"),
        ("TJ15_protein", "GEP18_RNA"),
        ("TJ15_RNA", "ImmuneScore"),
        ("TJ15_RNA", "CD8_RNA"),
        ("TJ15_RNA", "GEP18_RNA"),
    ]:
        hit = assoc[(assoc.predictor == pred) & (assoc.endpoint == end)]
        if len(hit):
            summary["headline"][f"{pred}__{end}"] = hit.iloc[0].to_dict()
    return {"summary": summary, "scores": scores, "assoc": assoc, "hl": hl, "cols": cols}


# ---------------------------------------------------------------------------
# Remaining open proteomes
# ---------------------------------------------------------------------------
def analyze_remaining(data: Path, out: Path) -> pd.DataFrame:
    rows = []
    coverage_rows = []
    for cohort in REMAINING_PROTEOME:
        pdir = data / "cptac" / cohort
        prot_path = pdir / f"{cohort}_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt"
        pheno_path = pdir / f"{cohort}_phenotype.txt"
        if not prot_path.exists() or not pheno_path.exists():
            rows.append({"cohort": cohort, "status": "missing_files"})
            continue
        prot = load_matrix(prot_path)
        pheno = load_pheno(pheno_path)
        prot, pheno_m = join_protein_pheno(prot, pheno)
        cols = pick_immune_cols(pheno_m)
        cldn4 = extract_gene(prot, "CLDN4")
        tj15, cov = mean_z(prot, TJ_15, MIN_TJ_GENES)
        cov = cov.copy()
        cov.insert(0, "cohort", cohort)
        coverage_rows.append(cov)
        n_cldn4 = int(cldn4.notna().sum())
        ends = {}
        for key in ["ImmuneScore", "CD8_CIBERSORT", "CD8_xCell", "Immune_xCell"]:
            if key in cols:
                ends[key] = pd.to_numeric(pheno_m[cols[key]], errors="coerce")
        for pred_name, pred in [("CLDN4_protein", cldn4), ("TJ15_protein", tj15)]:
            for ename, ey in ends.items():
                sp = spearman(pred, ey)
                rec = {
                    "cohort": cohort,
                    "n_tumors": int(prot.shape[1]),
                    "n_CLDN4_protein": n_cldn4,
                    "tj15_n_genes": int((cov["n_obs"] >= 8).sum()),
                    "predictor": pred_name,
                    "endpoint": ename,
                    "pheno_column": cols.get(ename),
                    **sp,
                }
                rows.append(rec)
    assoc = pd.DataFrame(rows)
    assoc.to_csv(out / "tables/remaining_proteome_associations.tsv", sep="\t", index=False)
    if coverage_rows:
        pd.concat(coverage_rows, ignore_index=True).to_csv(
            out / "tables/remaining_proteome_tj_coverage.tsv", sep="\t", index=False
        )
    return assoc


# ---------------------------------------------------------------------------
# TCGA-LUSC
# ---------------------------------------------------------------------------
def analyze_lusc(data: Path, out: Path) -> dict:
    expr = load_matrix(data / "tcga/TCGA.LUSC.HiSeqV2.gz")
    est = pd.read_csv(data / "tcga/ESTIMATE_LUSC_RNAseqV2.txt", sep="\t")
    # MD Anderson tables use ID / Name-like first column
    idcol = est.columns[0]
    est = est.rename(columns={idcol: "sample"}).set_index("sample")
    # normalize score column names
    rename = {}
    for c in est.columns:
        cl = c.lower().replace(" ", "_")
        if "stromal" in cl:
            rename[c] = "Stromal_score"
        elif "immune" in cl:
            rename[c] = "Immune_score"
        elif "estimate" in cl:
            rename[c] = "ESTIMATE_score"
    est = est.rename(columns=rename)

    expr_cols = {c[:15]: c for c in expr.columns}
    common = sorted({s for s in est.index if s in expr_cols and str(s).endswith("-01")})
    expr_m = expr[[expr_cols[s] for s in common]]
    expr_m.columns = common
    est_m = est.loc[common]

    tj15, cov15 = mean_z(expr_m, TJ_15, MIN_TJ_GENES)
    tj7, cov7 = mean_z(expr_m, TJ_7, MIN_TJ7_GENES)
    cd8, _ = mean_z(expr_m, CD8, 2)
    cyt, _ = mean_z(expr_m, CYT, 2)
    gep, covg = mean_z(expr_m, GEP18, 10)
    purity = estimate_purity(pd.to_numeric(est_m["ESTIMATE_score"], errors="coerce"))
    immune = pd.to_numeric(est_m["Immune_score"], errors="coerce")

    scores = pd.DataFrame({
        "sample": common,
        "TJ15": tj15.reindex(common).values,
        "TJ7": tj7.reindex(common).values,
        "CD8": cd8.reindex(common).values,
        "CYT": cyt.reindex(common).values,
        "GEP18": gep.reindex(common).values,
        "Immune_score": immune.reindex(common).values,
        "ESTIMATE_score": pd.to_numeric(est_m["ESTIMATE_score"], errors="coerce").reindex(common).values,
        "purity": purity.reindex(common).values,
    }).set_index("sample")

    rows = []
    for pred in ["TJ15", "TJ7"]:
        for end in ["CD8", "CYT", "GEP18", "Immune_score"]:
            sp = spearman(scores[pred], scores[end])
            part = partial_spearman(scores[pred], scores[end], scores["purity"])
            rec = {"cohort": "TCGA_LUSC", "predictor": pred, "endpoint": end, "test": "spearman", **sp}
            rec.update({f"partial_{k}": v for k, v in part.items()})
            rows.append(rec)
    assoc = pd.DataFrame(rows)
    hl_rows = []
    for pred in ["TJ15", "TJ7"]:
        for end in ["CD8", "CYT", "GEP18"]:
            rec = {"cohort": "TCGA_LUSC", "predictor": pred, "endpoint": end, **high_low(scores[pred], scores[end])}
            hl_rows.append(rec)
    hl = pd.DataFrame(hl_rows)

    scores.to_csv(out / "tables/lusc_sample_scores.tsv", sep="\t")
    cov15.to_csv(out / "tables/lusc_tj15_coverage.tsv", sep="\t", index=False)
    cov7.to_csv(out / "tables/lusc_tj7_coverage.tsv", sep="\t", index=False)
    covg.to_csv(out / "tables/lusc_gep_coverage.tsv", sep="\t", index=False)
    assoc.to_csv(out / "tables/lusc_associations.tsv", sep="\t", index=False)
    hl.to_csv(out / "tables/lusc_highlow.tsv", sep="\t", index=False)

    summary = {
        "cohort": "TCGA_LUSC",
        "n_primary": int(len(scores)),
        "n_with_purity": int(scores["purity"].notna().sum()),
        "tj15_genes": cov15.loc[cov15["found"], "gene"].tolist(),
        "headline": {},
        "highlow": {},
    }
    for pred, end in [("TJ15", "CD8"), ("TJ15", "GEP18"), ("TJ15", "CYT"), ("TJ7", "CD8"), ("TJ7", "GEP18")]:
        hit = assoc[(assoc.predictor == pred) & (assoc.endpoint == end)]
        if len(hit):
            summary["headline"][f"{pred}__{end}"] = hit.iloc[0].to_dict()
        hhit = hl[(hl.predictor == pred) & (hl.endpoint == end)]
        if len(hhit):
            summary["highlow"][f"{pred}__{end}"] = hhit.iloc[0].to_dict()
    return {"summary": summary, "scores": scores, "assoc": assoc, "hl": hl}


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
def _scatter(ax, x, y, title, xlabel, ylabel, color=None, color_label=None):
    d = pd.concat([x, y], axis=1).dropna()
    d.columns = ["x", "y"]
    sp = spearman(d["x"], d["y"])
    c = color.reindex(d.index) if color is not None else None
    if c is not None and c.notna().sum() > 10:
        sc = ax.scatter(d["x"], d["y"], c=c, cmap="viridis", s=18, alpha=0.75, edgecolors="none")
        plt.colorbar(sc, ax=ax, fraction=0.046, pad=0.04).set_label(color_label or "covariate")
    else:
        ax.scatter(d["x"], d["y"], s=18, alpha=0.7, c="#2c5aa0", edgecolors="none")
    if len(d) >= 6:
        m, b = np.polyfit(d["x"].rank(), d["y"].rank(), 1)
        # visual aid: least-squares on raw values
        mr, br = np.polyfit(d["x"], d["y"], 1)
        xs = np.linspace(d["x"].min(), d["x"].max(), 50)
        ax.plot(xs, mr * xs + br, color="#c0392b", lw=1.2, alpha=0.8)
    ax.set_title(title, fontsize=10)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.text(
        0.03,
        0.97,
        f"n={sp['n']}   ρ={fmt_rho(sp['rho'])}   p={fmt_p(sp['p'])}",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=8,
        bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.7", alpha=0.9),
    )


def fig_luad(luad: dict, figdir: Path) -> None:
    s = luad["scores"]
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.0), constrained_layout=True)
    _scatter(
        axes[0], s["CLDN4_protein"], s.get("ImmuneScore"),
        "CPTAC LUAD  ·  CLDN4 protein vs ImmuneScore",
        "CLDN4 protein (log2 TMT)", "ESTIMATE ImmuneScore",
        s.get("WES_purity"),
        "WES purity",
    )
    cd8 = s["CD8_xCell"] if "CD8_xCell" in s.columns else s.get("CD8_CIBERSORT")
    _scatter(
        axes[1], s["CLDN4_protein"], cd8,
        "CPTAC LUAD  ·  CLDN4 protein vs CD8",
        "CLDN4 protein (log2 TMT)", "CD8 (xCell)",
        s.get("WES_purity"),
        "WES purity",
    )
    fig.savefig(figdir / "fig1_luad_cldn4_protein.png", dpi=160)
    fig.savefig(figdir / "fig1_luad_cldn4_protein.pdf")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.0), constrained_layout=True)
    _scatter(
        axes[0], s["TJ15_protein"], s.get("ImmuneScore"),
        "CPTAC LUAD  ·  TJ-15 protein vs ImmuneScore",
        "TJ-15 protein (mean z)", "ESTIMATE ImmuneScore",
        s.get("WES_purity"),
        "WES purity",
    )
    _scatter(
        axes[1], s["TJ15_protein"], s["CD8_RNA"],
        "CPTAC LUAD  ·  TJ-15 protein vs CD8 RNA",
        "TJ-15 protein (mean z)", "CD8A/CD8B RNA (mean z)",
        s.get("WES_purity"),
        "WES purity",
    )
    fig.savefig(figdir / "fig2_luad_tj15_protein.png", dpi=160)
    fig.savefig(figdir / "fig2_luad_tj15_protein.pdf")
    plt.close(fig)

    # forest of LUAD headline associations
    want = [
        ("CLDN4_protein", "ImmuneScore"),
        ("CLDN4_protein", "CD8_CIBERSORT"),
        ("CLDN4_protein", "CD8_xCell"),
        ("CLDN4_protein", "CD8_RNA"),
        ("CLDN4_protein", "GEP18_RNA"),
        ("TJ15_protein", "ImmuneScore"),
        ("TJ15_protein", "CD8_CIBERSORT"),
        ("TJ15_protein", "CD8_xCell"),
        ("TJ15_protein", "CD8_RNA"),
        ("TJ15_protein", "GEP18_RNA"),
        ("TJ15_RNA", "ImmuneScore"),
        ("TJ15_RNA", "CD8_RNA"),
        ("TJ15_RNA", "GEP18_RNA"),
    ]
    a = luad["assoc"]
    labels, rhos, los, his = [], [], [], []
    for pred, end in want:
        hit = a[(a.predictor == pred) & (a.endpoint == end)]
        if not len(hit) or pd.isna(hit.iloc[0]["rho"]):
            continue
        r = hit.iloc[0]
        labels.append(f"{pred} vs {end}  (n={int(r['n'])})")
        rhos.append(r["rho"])
        los.append(r["ci95_lo"])
        his.append(r["ci95_hi"])
    fig, ax = plt.subplots(figsize=(8.4, 0.42 * max(len(labels), 1) + 1.2), constrained_layout=True)
    y = np.arange(len(labels))[::-1]
    ax.axvline(0, color="0.5", lw=0.8)
    ax.errorbar(rhos, y, xerr=[np.array(rhos) - np.array(los), np.array(his) - np.array(rhos)],
                fmt="o", color="#1f4e79", ecolor="#1f4e79", capsize=2, ms=5)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Spearman ρ (95% bootstrap CI)")
    ax.set_title("CPTAC LUAD — CLDN4 / TJ protein and RNA vs immune readouts")
    fig.savefig(figdir / "fig3_luad_forest.png", dpi=160)
    fig.savefig(figdir / "fig3_luad_forest.pdf")
    plt.close(fig)


def fig_remaining(assoc: pd.DataFrame, figdir: Path) -> None:
    a = assoc.dropna(subset=["rho"]).copy()
    if a.empty:
        return
    # primary extra: CLDN4 protein vs ImmuneScore, then TJ15 vs ImmuneScore
    sub = a[a.endpoint == "ImmuneScore"].copy()
    if sub.empty:
        sub = a.copy()
    sub = sub.sort_values(["predictor", "cohort"])
    labels = [f"{r.cohort}  {r.predictor}  n={int(r.n)}" for r in sub.itertuples()]
    fig, ax = plt.subplots(figsize=(8.2, 0.38 * max(len(sub), 1) + 1.3), constrained_layout=True)
    y = np.arange(len(sub))[::-1]
    ax.axvline(0, color="0.5", lw=0.8)
    colors = ["#1f4e79" if p == "CLDN4_protein" else "#b35c00" for p in sub["predictor"]]
    ax.errorbar(
        sub["rho"], y,
        xerr=[sub["rho"] - sub["ci95_lo"], sub["ci95_hi"] - sub["rho"]],
        fmt="none", ecolor="0.45", capsize=2,
    )
    ax.scatter(sub["rho"], y, c=colors, s=28, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Spearman ρ vs ImmuneScore (95% CI)")
    ax.set_title("Remaining open CPTAC freeze proteomes (LSCC not re-run)")
    fig.savefig(figdir / "fig4_remaining_proteome_forest.png", dpi=160)
    fig.savefig(figdir / "fig4_remaining_proteome_forest.pdf")
    plt.close(fig)


def fig_lusc(lusc: dict, figdir: Path) -> None:
    s = lusc["scores"]
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.0), constrained_layout=True)
    _scatter(
        axes[0], s["TJ15"], s["CD8"],
        "TCGA-LUSC  ·  TJ-15 vs CD8",
        "TJ-15 RNA (mean z)", "CD8A/CD8B (mean z)",
        s["purity"],
        "ESTIMATE purity",
    )
    _scatter(
        axes[1], s["TJ15"], s["GEP18"],
        "TCGA-LUSC  ·  TJ-15 vs GEP18",
        "TJ-15 RNA (mean z)", "Ayers GEP18 (mean z)",
        s["purity"],
        "ESTIMATE purity",
    )
    fig.savefig(figdir / "fig5_lusc_tj_cd8_gep.png", dpi=160)
    fig.savefig(figdir / "fig5_lusc_tj_cd8_gep.pdf")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.8), constrained_layout=True)
    for ax, end, title in [
        (axes[0], "CD8", "TCGA-LUSC  ·  CD8 by TJ-15 median split"),
        (axes[1], "GEP18", "TCGA-LUSC  ·  GEP18 by TJ-15 median split"),
    ]:
        d = pd.concat([s["TJ15"], s[end]], axis=1).dropna()
        d.columns = ["tj", "y"]
        med = d["tj"].median()
        hi = d.loc[d["tj"] > med, "y"]
        lo = d.loc[d["tj"] <= med, "y"]
        hl = high_low(s["TJ15"], s[end])
        ax.boxplot([lo.values, hi.values], tick_labels=["TJ-low", "TJ-high"], widths=0.55,
                   patch_artist=True,
                   boxprops=dict(facecolor="#d6e4f0", edgecolor="#1f4e79"),
                   medianprops=dict(color="#c0392b", lw=1.6))
        ax.set_ylabel(end)
        ax.set_title(title, fontsize=10)
        ax.text(
            0.5, 0.97,
            f"n_low={hl['n_low']}  n_high={hl['n_high']}  MWU p={fmt_p(hl['p'])}",
            transform=ax.transAxes, ha="center", va="top", fontsize=8,
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.7", alpha=0.9),
        )
    fig.savefig(figdir / "fig6_lusc_tj_highlow.png", dpi=160)
    fig.savefig(figdir / "fig6_lusc_tj_highlow.pdf")
    plt.close(fig)


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--data", default="data/b3_extra_proteome_lusc")
    p.add_argument("--out", default="results/B3_extra_proteome_lusc")
    args = p.parse_args()
    data = Path(args.data)
    out = Path(args.out)
    (out / "tables").mkdir(parents=True, exist_ok=True)
    (out / "figures").mkdir(parents=True, exist_ok=True)

    luad = analyze_luad(data, out)
    remaining = analyze_remaining(data, out)
    lusc = analyze_lusc(data, out)

    fig_luad(luad, out / "figures")
    fig_remaining(remaining, out / "figures")
    fig_lusc(lusc, out / "figures")

    summary = {
        "luad": luad["summary"],
        "lusc": lusc["summary"],
        "remaining_n_rows": int(len(remaining)),
        "notes": {
            "b3_tcga_luad": "Taken as given; not re-run.",
            "cptac_lscc": "CLDN4 protein vs ImmuneScore already public; not re-audited.",
        },
    }
    (out / "tables/summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")
    print(json.dumps(summary, indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
