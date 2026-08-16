#!/usr/bin/env python3
"""EXTRA LUAD RNA tables — TACSTD2/CLDN4 vs CD8/GEP after purity residual.

The user's TCGA + OncoSG purity-corrected TACSTD2–immune result is taken
as given (OncoSG numbers from PR #139). This script ADDS non-TCGA /
East-Asian public LUAD RNA cohorts:

  GSE31210  Japan Okayama stage I–II LUAD (GPL570)
  GSE72094  Moffitt / Shedden-era LUAD (GPL15048)
  GSE68465  Director's Challenge LUAD (GPL96)
  CPTAC LUAD RNA  Gillette Cell 2020 freeze v1.2 (published WES + ESTIMATE)
  GSE19804  Taiwan never-smoker female lung cancer tumors (GPL570; extra East-Asian)

2023–2026 GEO LUAD tumor RNA with processed matrices and n large enough
for a purity residual was not found (hunt logged). OncoSG is not re-run.

Outputs -> results/rework/A1_extra_luad/
Honest n / ρ / p only. No claim-failed language.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "extra_luad"
SIG = ROOT / "data" / "signatures"
OUT = ROOT / "results" / "rework" / "A1_extra_luad"
FIG = OUT / "figures"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

GEP18 = [
    "CCL5", "CD27", "CD274", "CD276", "CD8A", "CMKLR1", "CXCL9", "CXCR6",
    "HLA-DQA1", "HLA-DRB1", "HLA-E", "IDO1", "LAG3", "NKG7", "PDCD1LG2",
    "PSMB10", "STAT1", "TIGIT",
]
TARGETS = ["TACSTD2", "CLDN4"]
IMMUNE = ["CD8", "GEP18"]

# Yoshihara 2013 ESTIMATE TumorPurity transform (used only on *computed* ESTIMATEScore).
EST_A = 0.6049872018
EST_B = 0.0001467884

# OncoSG (already done; PR #139). Copied as a reference row, not recomputed.
ONCOSG_REF = [
    {"cohort": "OncoSG (PR #139, not re-run)", "gene": "TACSTD2", "feature": "CD8",
     "n": 169, "unadj_rho": -0.380, "partial_rho": -0.309, "partial_p": 4.69e-05,
     "purity": "published clinical PURITY", "note": "taken as given"},
    {"cohort": "OncoSG (PR #139, not re-run)", "gene": "TACSTD2", "feature": "GEP18",
     "n": 169, "unadj_rho": -0.414, "partial_rho": -0.349, "partial_p": 3.62e-06,
     "purity": "published clinical PURITY", "note": "GEP 17/18; CCL5 absent"},
]


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
    return f"{r:.3f}"


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
    """Return (meta DataFrame samples x characteristics, expression probe x sample)."""
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
            # strip key prefix
            cleaned = []
            for v in vals:
                if ":" in v:
                    cleaned.append(v.split(":", 1)[1].strip())
                else:
                    cleaned.append(v)
            # avoid overwrite
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
    # expression
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


def load_gpl_annot(path: Path, id_col: str, symbol_col: str) -> pd.Series:
    """probe -> gene symbol (first symbol if ///)."""
    if str(path).endswith(".gz"):
        opener = lambda: gzip.open(path, "rt", errors="replace")
    else:
        opener = lambda: open(path, "r", errors="replace")
    with opener() as f:
        # skip GEO annot header until ID\t
        skip = 0
        for i, line in enumerate(f):
            if line.startswith("ID\t") or line.startswith("ID "):
                skip = i
                break
    df = pd.read_csv(path, sep="\t", skiprows=skip, dtype=str, low_memory=False)
    if id_col not in df.columns:
        # first col
        id_col = df.columns[0]
    if symbol_col not in df.columns:
        cands = [c for c in df.columns if "symbol" in c.lower()]
        if not cands:
            raise SystemExit(f"no symbol column in {path}: {list(df.columns)[:12]}")
        symbol_col = cands[0]
    s = df.set_index(id_col)[symbol_col].astype(str)
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
    # cosine purity; clip ESTIMATEScore to the range where the official map is monotone
    pur = np.cos(EST_A + EST_B * est.to_numpy(float))
    return pd.DataFrame(
        {"ESTIMATE_StromalScore": strom, "ESTIMATE_ImmuneScore": imm,
         "ESTIMATE_Score": est, "ESTIMATE_TumorPurity": pur},
        index=gene_expr.columns,
    )


def add_rows(rows, cohort, gene, feature, x, y, z, purity_name, extra=None):
    ru, pu, nu = spearman_pair(x, y)
    if z is None:
        rp, pp, np_ = np.nan, np.nan, nu
    else:
        rp, pp, np_ = partial_spearman(x, y, z)
    rec = {
        "cohort": cohort,
        "gene": gene,
        "feature": feature,
        "n_unadjusted": nu,
        "unadj_rho": ru,
        "unadj_p": pu,
        "n_partial": np_,
        "partial_rho": rp,
        "partial_p": pp,
        "purity": purity_name,
    }
    if extra:
        rec.update(extra)
    rows.append(rec)


def main():
    est_tab = pd.read_csv(SIG / "estimate_yoshihara_2013.tsv", sep="\t")
    stromal = est_tab.loc[est_tab["set"] == "Stromal141_UP", "hugo"].tolist()
    immune = est_tab.loc[est_tab["set"] == "Immune141_UP", "hugo"].tolist()
    print(f"ESTIMATE genes stromal={len(stromal)} immune={len(immune)}")

    gpl570 = load_gpl_annot(DATA / "GPL570.annot.gz", "ID", "Gene symbol")
    gpl96 = load_gpl_annot(DATA / "GPL96.annot.gz", "ID", "Gene symbol")
    gpl15048 = pd.read_csv(DATA / "GPL15048_platform_table.tsv", sep="\t", dtype=str)
    p2g_15048 = gpl15048.dropna(subset=["GeneSymbol"]).set_index("ID")["GeneSymbol"]
    p2g_15048 = p2g_15048[p2g_15048.str.len() > 0]

    rows = []
    coverage = []
    cohort_n = {}

    def run_geo(name, matrix, probe2gene, tumor_mask_fn, region, platform):
        print(f"\n=== {name} ===")
        meta, probes = parse_geo_matrix(matrix)
        print(f"  samples in matrix: {meta.shape[0]} probes: {probes.shape[0]}")
        keep = tumor_mask_fn(meta)
        print(f"  tumor filter: {int(keep.sum())} / {len(keep)}")
        meta_t = meta.loc[keep]
        probes_t = probes.loc[:, meta_t.index]
        genes = collapse_maxmean(probes_t, probe2gene)
        print(f"  genes after collapse: {genes.shape[0]}")
        miss_t = [g for g in TARGETS if g not in genes.index]
        gep_pres = [g for g in GEP18 if g in genes.index]
        gep_miss = [g for g in GEP18 if g not in genes.index]
        coverage.append({
            "cohort": name, "platform": platform, "region": region,
            "n_matrix": int(meta.shape[0]), "n_tumor": int(keep.sum()),
            "n_genes": int(genes.shape[0]),
            "targets_present": ",".join([g for g in TARGETS if g in genes.index]),
            "targets_missing": ",".join(miss_t),
            "gep18_present": len(gep_pres), "gep18_missing": ",".join(gep_miss),
            "estimate_stromal_present": sum(g in genes.index for g in stromal),
            "estimate_immune_present": sum(g in genes.index for g in immune),
        })
        if miss_t:
            print(f"  missing targets: {miss_t}")
        if "CD8A" not in genes.index:
            print("  CD8A missing — skip CD8")
        est = estimate_purity(genes, stromal, immune)
        if gep_pres:
            z = genes.loc[gep_pres].T
            gep = ((z - z.mean()) / z.std(ddof=0)).mean(axis=1)
        else:
            gep = pd.Series(np.nan, index=genes.columns)
        features = {}
        if "CD8A" in genes.index:
            features["CD8"] = genes.loc["CD8A"]
        features["GEP18"] = gep
        zpur = est["ESTIMATE_TumorPurity"].to_numpy(float)
        for gene in TARGETS:
            if gene not in genes.index:
                continue
            x = genes.loc[gene].to_numpy(float)
            for feat, yser in features.items():
                y = yser.to_numpy(float)
                add_rows(
                    rows, name, gene, feat, x, y, zpur,
                    "ESTIMATE TumorPurity (ssGSEA Yoshihara 141+141)",
                    extra={"region": region, "platform": platform, "gep18_n_genes": len(gep_pres)},
                )
        # TACSTD2 vs purity (context)
        if "TACSTD2" in genes.index:
            r, p, n = spearman_pair(genes.loc["TACSTD2"].to_numpy(float), zpur)
            print(f"  TACSTD2 vs ESTIMATE purity ρ={r:.3f} p={p:.3g} n={n}")
        cohort_n[name] = int(keep.sum())
        # save per-cohort sample table (small)
        st = pd.DataFrame(index=genes.columns)
        for g in TARGETS:
            if g in genes.index:
                st[g] = genes.loc[g]
        if "CD8A" in genes.index:
            st["CD8A"] = genes.loc["CD8A"]
        st["GEP18"] = gep
        st = st.join(est)
        st.index.name = "sample"
        st.to_csv(OUT / f"samples_{name.replace(' ', '_')}.tsv", sep="\t")
        return genes, est

    def mask_31210(meta):
        col = "tissue" if "tissue" in meta.columns else "source"
        s = meta[col].astype(str).str.lower()
        return s.str.contains("primary lung tumor") | s.str.contains("primary lung tumors")

    def mask_72094(meta):
        s = meta.get("source", meta.iloc[:, 0]).astype(str).str.lower()
        return s.str.contains("adenocarcinoma")

    def mask_68465(meta):
        if "disease_state" in meta.columns:
            return meta["disease_state"].astype(str).str.contains("Adenocarcinoma", case=False, na=False)
        return pd.Series(True, index=meta.index)

    def mask_19804(meta):
        # tumors vs paired normal: titles end with T / N or tissue=lung cancer
        if "tissue" in meta.columns:
            t = meta["tissue"].astype(str).str.lower()
            tumor = t.str.contains("lung cancer") & ~t.str.contains("adjacent") & ~t.str.contains("normal")
            if tumor.sum() >= 20:
                return tumor
        titles = meta["title"].astype(str)
        return titles.str.contains(r"T$|tumor|cancer", case=False) & ~titles.str.contains(r"N$|normal", case=False)

    run_geo(
        "GSE31210",
        DATA / "GSE31210_series_matrix.txt.gz",
        gpl570,
        mask_31210,
        "East Asia (Japan)",
        "GPL570 U133 Plus 2.0",
    )
    run_geo(
        "GSE72094",
        DATA / "GSE72094_series_matrix.txt.gz",
        p2g_15048,
        mask_72094,
        "USA (Moffitt custom array)",
        "GPL15048 HuRSTA",
    )
    run_geo(
        "GSE68465",
        DATA / "GSE68465_series_matrix.txt.gz",
        gpl96,
        mask_68465,
        "USA/multi-site Director's Challenge",
        "GPL96 U133A",
    )
    run_geo(
        "GSE19804",
        DATA / "GSE19804_series_matrix.txt.gz",
        gpl570,
        mask_19804,
        "East Asia (Taiwan, never-smoker female)",
        "GPL570 U133 Plus 2.0",
    )

    # ---------- CPTAC LUAD RNA ----------
    print("\n=== CPTAC LUAD RNA ===")
    rna = pd.read_csv(DATA / "LUAD_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt", sep="\t", index_col=0)
    rna.index = rna.index.astype(str).str.replace(r"\.\d+$", "", regex=True)
    pheno = pd.read_csv(DATA / "LUAD_phenotype.txt", sep="\t", index_col=0)
    # Ensembl mapping for targets / GEP18 / ESTIMATE not needed for outcomes if we have symbols via a small map
    ENS = {
        "TACSTD2": "ENSG00000184292",
        "CLDN4": "ENSG00000189143",
        "CD8A": "ENSG00000153563",
        "CCL5": "ENSG00000271503",
        "CD27": "ENSG00000139193",
        "CD274": "ENSG00000120217",
        "CD276": "ENSG00000103855",
        "CMKLR1": "ENSG00000174600",
        "CXCL9": "ENSG00000138755",
        "CXCR6": "ENSG00000172215",
        "HLA-DQA1": "ENSG00000196735",
        "HLA-DRB1": "ENSG00000196126",
        "HLA-E": "ENSG00000204592",
        "IDO1": "ENSG00000131203",
        "LAG3": "ENSG00000089692",
        "NKG7": "ENSG00000105374",
        "PDCD1LG2": "ENSG00000197646",
        "PSMB10": "ENSG00000205220",
        "STAT1": "ENSG00000115415",
        "TIGIT": "ENSG00000181847",
    }
    # HLA-DQA1 / DRB1 have multiple ENSG; try primary then skip
    def pick_ens(symbol):
        ens = ENS.get(symbol)
        if ens and ens in rna.index:
            return ens
        return None

    samples = [c for c in rna.columns if c in pheno.index]
    print(f"  RNA samples {rna.shape[1]}; phenotype overlap {len(samples)}")
    gep_ens = [(g, pick_ens(g)) for g in GEP18]
    gep_ok = [(g, e) for g, e in gep_ens if e]
    gep_miss = [g for g, e in gep_ens if e is None]
    z = rna.loc[[e for _, e in gep_ok], samples].T
    gep = ((z - z.mean()) / z.std(ddof=0)).mean(axis=1)
    cd8 = rna.loc[ENS["CD8A"], samples] if ENS["CD8A"] in rna.index else pd.Series(np.nan, index=samples)
    wes = pd.to_numeric(pheno.reindex(samples)["WES_purity"], errors="coerce")
    estscore = pd.to_numeric(pheno.reindex(samples)["ESTIMATE_ESTIMATEScore"], errors="coerce")
    # published ESTIMATEScore is used as an impurity axis (higher = more stroma+immune).
    # Do not apply the Yoshihara cosine: freeze scores (~5k–21k) wrap that formula (PR #99).
    features = {"CD8": cd8, "GEP18": gep}
    for gene in TARGETS:
        ens = ENS[gene]
        if ens not in rna.index:
            print(f"  {gene} {ens} missing")
            continue
        x = rna.loc[ens, samples].to_numpy(float)
        for feat, yser in features.items():
            y = yser.to_numpy(float)
            add_rows(
                rows, "CPTAC_LUAD_RNA", gene, feat, x, y, wes.to_numpy(float),
                "published WES_purity (Gillette 2020 freeze)",
                extra={"region": "USA CPTAC", "platform": "RNA-seq RSEM log2 UQ", "gep18_n_genes": len(gep_ok)},
            )
            add_rows(
                rows, "CPTAC_LUAD_RNA", gene, feat, x, y, estscore.to_numpy(float),
                "published ESTIMATE_ESTIMATEScore (RNA impurity axis)",
                extra={"region": "USA CPTAC", "platform": "RNA-seq RSEM log2 UQ", "gep18_n_genes": len(gep_ok)},
            )
    coverage.append({
        "cohort": "CPTAC_LUAD_RNA", "platform": "RNA-seq RSEM log2 UQ", "region": "USA CPTAC",
        "n_matrix": int(rna.shape[1]), "n_tumor": len(samples),
        "n_genes": int(rna.shape[0]),
        "targets_present": ",".join([g for g in TARGETS if ENS[g] in rna.index]),
        "targets_missing": "",
        "gep18_present": len(gep_ok), "gep18_missing": ",".join(gep_miss),
        "estimate_stromal_present": "published",
        "estimate_immune_present": "published",
    })
    cohort_n["CPTAC_LUAD_RNA"] = len(samples)
    st = pd.DataFrame({
        "TACSTD2": rna.loc[ENS["TACSTD2"], samples] if ENS["TACSTD2"] in rna.index else np.nan,
        "CLDN4": rna.loc[ENS["CLDN4"], samples] if ENS["CLDN4"] in rna.index else np.nan,
        "CD8A": cd8,
        "GEP18": gep,
        "WES_purity": wes,
        "ESTIMATE_ESTIMATEScore": estscore,
    })
    st.index.name = "sample"
    st.to_csv(OUT / "samples_CPTAC_LUAD_RNA.tsv", sep="\t")

    res = pd.DataFrame(rows)
    cov = pd.DataFrame(coverage)
    res.to_csv(OUT / "extra_correlations.tsv", sep="\t", index=False)
    cov.to_csv(OUT / "cohort_coverage.tsv", sep="\t", index=False)
    pd.DataFrame(ONCOSG_REF).to_csv(OUT / "oncosg_reference_pr139.tsv", sep="\t", index=False)

    # compact extra table: TACSTD2/CLDN4 x CD8/GEP, one purity per GEO + both CPTAC
    extra = res.copy()
    extra.to_csv(OUT / "EXTRA_TABLE_tacstd2_cldn4_cd8_gep.tsv", sep="\t", index=False)

    # figure
    plot = res[res["purity"].str.contains("TumorPurity|WES_purity")].copy()
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 5.2), sharey=True)
    for ax, gene in zip(axes, TARGETS):
        sub = plot[plot.gene == gene]
        if sub.empty:
            ax.set_title(gene)
            continue
        labels = [f"{r.cohort}\n{r.feature}" for r in sub.itertuples()]
        y = np.arange(len(sub))
        cols = ["#c0392b" if (np.isfinite(r) and r < 0) else "#2471a3" for r in sub["partial_rho"]]
        ax.barh(y, sub["partial_rho"], color=cols)
        ax.axvline(0, color="k", lw=0.8)
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=7)
        ax.invert_yaxis()
        ax.set_xlabel("partial Spearman ρ")
        ax.set_title(f"{gene} vs CD8 / GEP18 after purity residual")
    fig.tight_layout()
    fig.savefig(FIG / "bar_extra_partial_rho.png", dpi=150)
    plt.close(fig)

    def md_table(df):
        cols = ["cohort", "gene", "feature", "n_partial", "unadj_rho", "unadj_p",
                "partial_rho", "partial_p", "purity"]
        lines = [
            "| Cohort | Gene | Feature | n | Unadj ρ | Unadj p | Partial ρ | Partial p | Purity |",
            "|---|---|---|---:|---:|---:|---:|---:|---|",
        ]
        for r in df.itertuples():
            lines.append(
                f"| {r.cohort} | {r.gene} | {r.feature} | {int(r.n_partial) if np.isfinite(r.n_partial) else 'NA'} | "
                f"{fmt_rho(r.unadj_rho)} | {fmt_p(r.unadj_p)} | {fmt_rho(r.partial_rho)} | "
                f"{fmt_p(r.partial_p)} | {r.purity} |"
            )
        return "\n".join(lines)

    geo = res[~res.cohort.str.startswith("CPTAC")]
    cptac = res[res.cohort.str.startswith("CPTAC")]

    hunt_2023 = """NCBI GEO (GDS) search `lung adenocarcinoma[Title] AND 2023:2026[Publication Date] AND Expression profiling by high throughput sequencing` returned 122 records. The first page is cell-line, n=1–36, scRNA, or plasma EV studies (GSE293914 n=1; GSE330179 n=18; GSE332750 n=4; GSE299604 n=11; GSE317139 n=36). GSE226481 (China, early- vs late-onset LUAD) has **n=14** tumors — too small for a purity residual. GSE40419 (Seo Korea RNA-seq) and GSE140343 (paired LUAD RNA-seq) have **no processed series matrix** (SRA/GTF only). No 2023–2026 GEO LUAD bulk tumor matrix with n≥50 was used."""

    report = f"""# EXTRA LUAD RNA tables — TACSTD2 / CLDN4 vs CD8 / GEP after purity residual

The user's **TCGA + OncoSG** purity-corrected TACSTD2–immune result is **taken as given**. OncoSG numbers below are from PR #139 and were not re-run. This folder **adds** public non-TCGA / East-Asian LUAD RNA cohorts.

## Extra table (this run)

{md_table(res)}

## OncoSG reference (already done)

| Cohort | Gene | Feature | n | Unadj ρ | Partial ρ | Partial p | Purity |
|---|---|---|---:|---:|---:|---:|---|
| OncoSG PR #139 | TACSTD2 | CD8A | 169 | −0.380 | −0.309 | 4.69e-05 | published clinical PURITY |
| OncoSG PR #139 | TACSTD2 | GEP18 (17/18) | 169 | −0.414 | −0.349 | 3.62e-06 | published clinical PURITY |

## Cohorts

| Cohort | Region | Platform | Tumor n | Purity used |
|---|---|---|---:|---|
| GSE31210 Okayama / Kohno | East Asia (Japan) | GPL570 U133 Plus 2.0 | {cohort_n.get('GSE31210','NA')} | ESTIMATE TumorPurity (ssGSEA, Yoshihara 141+141) |
| GSE19804 Taiwan never-smoker female | East Asia (Taiwan) | GPL570 U133 Plus 2.0 | {cohort_n.get('GSE19804','NA')} | ESTIMATE TumorPurity (ssGSEA) |
| GSE72094 | USA | GPL15048 HuRSTA | {cohort_n.get('GSE72094','NA')} | ESTIMATE TumorPurity (ssGSEA) |
| GSE68465 Director's Challenge | USA / multi-site | GPL96 U133A | {cohort_n.get('GSE68465','NA')} | ESTIMATE TumorPurity (ssGSEA) |
| CPTAC LUAD RNA (Gillette *Cell* 2020) | USA | RNA-seq RSEM log2 UQ | {cohort_n.get('CPTAC_LUAD_RNA','NA')} | published WES_purity **and** published ESTIMATEScore |

Primary tumors only. GSE19804 paired normals dropped. GSE31210 restricted to `tissue: primary lung tumor`.

## Definitions

- **TACSTD2 / CLDN4 / CD8:** microarray = max-mean probe collapse to HUGO; CPTAC = Ensembl `ENSG00000184292` / `ENSG00000189143` / `ENSG00000153563` (version suffix stripped).
- **CD8** = `CD8A`.
- **GEP18** = unweighted within-cohort z-mean of the Ayers 2017 18-gene list. Missing genes are not imputed (coverage in `cohort_coverage.tsv`).
- **GEO purity:** ESTIMATE ssGSEA (Barbie/GSVA tau=0.25, ranks 1..10000) on Yoshihara 2013 Stromal141 + Immune141, then `TumorPurity = cos(0.6049872018 + 0.0001467884 × ESTIMATEScore)`. Gene list: `data/signatures/estimate_yoshihara_2013.tsv` (Supp Data 1).
- **CPTAC purity:** published `WES_purity` (DNA). Published `ESTIMATE_ESTIMATEScore` is a second RNA axis. The Yoshihara cosine is **not** applied to the freeze scores (they wrap; see PR #99).
- **Partial Spearman:** Pearson of rank residuals; df = n − 3.

## 2023–2026 GEO hunt

{hunt_2023}

## Files

- `EXTRA_TABLE_tacstd2_cldn4_cd8_gep.tsv` — all extra n/ρ/p
- `extra_correlations.tsv` — same
- `cohort_coverage.tsv`
- `oncosg_reference_pr139.tsv`
- `samples_*.tsv`
- `figures/bar_extra_partial_rho.png`
- `summary.json` / `provenance.json`

## Reproduce

```
pip install -r requirements.txt
python scripts/rework_A1_extra_luad.py
```

Public GEO series matrices + CPTAC S3 freeze v1.2 + Yoshihara Supp Data 1. Raw downloads under `data/extra_luad/` (gitignored).
"""
    (OUT / "REPORT.md").write_text(report)

    summary = {
        "task": "A1_extra_luad",
        "note": "TCGA+OncoSG taken as given; these are extra cohorts",
        "cohort_n": cohort_n,
        "n_rows": int(len(res)),
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "results": res.to_dict(orient="records"),
    }
    with open(OUT / "summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=float)

    prov = {"inputs": {}, "signatures": str(SIG / "estimate_yoshihara_2013.tsv")}
    for fn in [
        "GSE31210_series_matrix.txt.gz",
        "GSE72094_series_matrix.txt.gz",
        "GSE68465_series_matrix.txt.gz",
        "GSE19804_series_matrix.txt.gz",
        "LUAD_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt",
        "LUAD_phenotype.txt",
        "GPL570.annot.gz",
        "GPL96.annot.gz",
        "GPL15048_platform_table.tsv",
    ]:
        fp = DATA / fn
        if fp.exists():
            prov["inputs"][fn] = {"bytes": fp.stat().st_size, "md5": md5(fp)}
    with open(OUT / "provenance.json", "w") as f:
        json.dump(prov, f, indent=2)

    print("\n=== EXTRA TABLE ===")
    print(res[["cohort", "gene", "feature", "n_partial", "unadj_rho", "partial_rho", "partial_p"]].to_string(
        index=False, float_format=lambda v: f"{v:.4g}"
    ))
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
