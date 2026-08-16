#!/usr/bin/env python3
"""
Claims B3 + B4 — independent, honest reproduction attempt.

B3 (TCGA-LUAD): Is a tight-junction (TJ) gene-expression signature associated
    with a CD8 T-cell signature and with the T-cell-inflamed GEP signature?
    User claim: p < 1e-6.

B4 (GSE126044, NSCLC pre-treatment anti-PD-1): Do non-responders (NR) show a
    HIGHER TJ signature score than responders (R)?
    User claim: p = 0.019.

Design notes / honesty caveats
------------------------------
* The exact gene membership of the "TJ signature" was NOT provided with the
  task. We therefore anchor it to citable, reproducible sources:
    - TJ_core   : curated core *structural* tight-junction components
                  (claudins, occludin, ZO-1/2/3, JAMs, cingulin, PATJ, ...).
                  This is the PRIMARY signature.
    - TJ_kegg   : the full KEGG "Tight junction" pathway (hsa04530), fetched
                  live from the KEGG REST API. Reported as a SENSITIVITY check
                  (note: it contains many actin/myosin/tubulin genes and is
                  therefore less TJ-specific).
* CD8 signature = mean of CD8A, CD8B (canonical cytotoxic T-cell markers).
* GEP = the 18-gene T-cell-inflamed Gene Expression Profile (Ayers et al.,
  J Clin Invest 2017). We use the unweighted z-score mean as a transparent
  proxy for the published (housekeeping-normalized, weighted) score.
* Signature score = mean across signature genes of the per-gene z-score
  (gene standardized across samples). Simple, transparent, widely used.
* We report BOTH Spearman and Pearson for B3, and BOTH two-sided and the
  directional (NR>R) one-sided Mann-Whitney U for B4. We report the ACTUAL
  numbers, whether or not they match the user's stated p-values.
"""

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

DATA = Path("/workspace/data")
OUT = Path("/workspace/results/claim_B3B4")
OUT.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Signature definitions
# ---------------------------------------------------------------------------
TJ_CORE = [
    # Claudins (barrier-forming TJ strand proteins)
    "CLDN1", "CLDN2", "CLDN3", "CLDN4", "CLDN5", "CLDN7", "CLDN8",
    "CLDN10", "CLDN11", "CLDN15", "CLDN18", "CLDN23",
    # Occludin + TJ-associated MARVEL proteins
    "OCLN", "MARVELD2", "MARVELD3",
    # Zonula occludens (scaffolding)
    "TJP1", "TJP2", "TJP3",
    # Junctional adhesion molecules
    "F11R", "JAM2", "JAM3",
    # Cingulin / polarity / scaffolds
    "CGN", "CGNL1", "CRB3", "PARD3", "PARD6A", "MPDZ", "PATJ", "SYMPK",
]

CD8 = ["CD8A", "CD8B"]

# Ayers et al. 2017 18-gene T-cell-inflamed GEP
GEP = [
    "CCL5", "CD27", "CD274", "CD276", "CD8A", "CMKLR1", "CXCL9", "CXCR6",
    "HLA-DQA1", "HLA-DRB1", "HLA-E", "IDO1", "LAG3", "NKG7", "PDCD1LG2",
    "PSMB10", "STAT1", "TIGIT",
]

# Common symbol synonyms seen in older TCGA/GEO annotations
SYNONYMS = {
    "PATJ": ["INADL"],
    "PDCD1LG2": ["PDL2", "PDCD1LG2"],
    "CD276": ["B7H3"],
    "F11R": ["JAM1", "JAMA"],
    "HLA-DQA1": ["HLA.DQA1"],
    "HLA-DRB1": ["HLA.DRB1"],
    "HLA-E": ["HLA.E"],
}


def parse_kegg_symbols(path: Path) -> list:
    """Extract HGNC symbols from a KEGG `get` flat file GENE section."""
    syms, in_gene = [], False
    for line in path.read_text().splitlines():
        if line.startswith("GENE"):
            in_gene = True
        elif in_gene and line[:1].strip() and not line.startswith(" "):
            break  # next top-level section
        if in_gene:
            m = re.search(r"\d+\s+([A-Za-z0-9\-]+);", line)
            if m:
                syms.append(m.group(1))
    return sorted(set(syms))


TJ_KEGG = parse_kegg_symbols(DATA / "kegg_hsa04530.txt")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def resolve_genes(wanted, available):
    """Map wanted symbols to those available, using synonyms. Returns
    (matched_available_symbols, missing_wanted_symbols)."""
    avail = set(available)
    matched, missing = [], []
    for g in wanted:
        cand = [g] + SYNONYMS.get(g, [])
        hit = next((c for c in cand if c in avail), None)
        if hit is not None:
            matched.append(hit)
        else:
            missing.append(g)
    return matched, missing


def zscore_rows(df: pd.DataFrame) -> pd.DataFrame:
    """Z-score each gene (row) across samples (columns)."""
    mu = df.mean(axis=1)
    sd = df.std(axis=1, ddof=0).replace(0, np.nan)
    return df.sub(mu, axis=0).div(sd, axis=0)


def signature_score(expr: pd.DataFrame, genes) -> pd.Series:
    """Mean per-gene z-score across signature genes. expr = genes x samples."""
    matched, missing = resolve_genes(genes, expr.index)
    z = zscore_rows(expr.loc[matched])
    return z.mean(axis=0), matched, missing


# ---------------------------------------------------------------------------
# Load TCGA-LUAD (genes x samples, log2(norm_count+1))
# ---------------------------------------------------------------------------
tcga = pd.read_csv(DATA / "TCGA_LUAD_HiSeqV2", sep="\t", index_col=0)
tcga.index = tcga.index.astype(str)
# keep primary tumor samples only: TCGA barcode sample-type code == '01'
sample_type = pd.Series(tcga.columns).str.split("-").str[3].str[:2].values
tumor_cols = tcga.columns[sample_type == "01"]
tcga_t = tcga[tumor_cols]

# ---------------------------------------------------------------------------
# Load GSE126044 (raw counts, genes x samples) -> log2 CPM
# ---------------------------------------------------------------------------
gse = pd.read_csv(DATA / "GSE126044_counts.txt", sep="\t", index_col=0)
gse.index = gse.index.astype(str)
lib = gse.sum(axis=0)
cpm = gse.div(lib, axis=1) * 1e6
logcpm = np.log2(cpm + 1.0)

# response labels from series matrix
RESPONDERS = {"Dis_02", "Dis_15", "Dis_04", "Dis_17", "Dis_10"}
labels = pd.Series(
    ["R" if c in RESPONDERS else "NR" for c in gse.columns], index=gse.columns
)

# ---------------------------------------------------------------------------
# Run analyses
# ---------------------------------------------------------------------------
results = {"datasets": {}, "signatures": {}, "B3": {}, "B4": {}}

results["datasets"]["TCGA_LUAD"] = {
    "source": "UCSC Xena TCGA.LUAD.sampleMap/HiSeqV2 (log2 norm_count+1, RSEM)",
    "n_genes": int(tcga.shape[0]),
    "n_samples_total": int(tcga.shape[1]),
    "n_primary_tumor": int(tcga_t.shape[1]),
}
results["datasets"]["GSE126044"] = {
    "source": "GEO GSE126044 raw counts (Cho et al., NSCLC pre-treatment anti-PD-1)",
    "n_genes": int(gse.shape[0]),
    "n_samples": int(gse.shape[1]),
    "n_responder": int((labels == "R").sum()),
    "n_nonresponder": int((labels == "NR").sum()),
    "normalization": "log2(CPM+1)",
}


def sig_report(name, genes, expr):
    score, matched, missing = signature_score(expr, genes)
    return score, {"n_requested": len(genes), "n_matched": len(matched),
                   "matched": matched, "missing": missing}


# ---- B3: TCGA-LUAD correlations ----
tj_tcga, tj_meta = sig_report("TJ_core", TJ_CORE, tcga_t)
tjk_tcga, tjk_meta = sig_report("TJ_kegg", TJ_KEGG, tcga_t)
cd8_tcga, cd8_meta = sig_report("CD8", CD8, tcga_t)
gep_tcga, gep_meta = sig_report("GEP", GEP, tcga_t)

results["signatures"]["TCGA"] = {
    "TJ_core": tj_meta, "TJ_kegg": tjk_meta, "CD8": cd8_meta, "GEP": gep_meta,
}


def corr_block(x, y):
    r_s, p_s = stats.spearmanr(x, y)
    r_p, p_p = stats.pearsonr(x, y)
    return {
        "n": int(len(x)),
        "spearman_rho": float(r_s), "spearman_p": float(p_s),
        "pearson_r": float(r_p), "pearson_p": float(p_p),
    }


for tjname, tjscore in [("TJ_core", tj_tcga), ("TJ_kegg", tjk_tcga)]:
    results["B3"][tjname] = {
        "vs_CD8": corr_block(tjscore, cd8_tcga),
        "vs_GEP": corr_block(tjscore, gep_tcga),
    }

# ---- B4: GSE126044 NR vs R ----
tj_gse, tj_gse_meta = sig_report("TJ_core", TJ_CORE, logcpm)
tjk_gse, tjk_gse_meta = sig_report("TJ_kegg", TJ_KEGG, logcpm)
results["signatures"]["GSE126044"] = {"TJ_core": tj_gse_meta, "TJ_kegg": tjk_gse_meta}


def group_test(score):
    nr = score[labels == "NR"].values
    r = score[labels == "R"].values
    u_two, p_two = stats.mannwhitneyu(nr, r, alternative="two-sided")
    u_gr, p_gr = stats.mannwhitneyu(nr, r, alternative="greater")  # NR > R
    t, p_t = stats.ttest_ind(nr, r, equal_var=False)
    return {
        "n_NR": int(len(nr)), "n_R": int(len(r)),
        "mean_NR": float(np.mean(nr)), "mean_R": float(np.mean(r)),
        "median_NR": float(np.median(nr)), "median_R": float(np.median(r)),
        "direction": "NR>R" if np.median(nr) > np.median(r) else "NR<=R",
        "mannwhitney_U_two_sided": float(u_two), "mannwhitney_p_two_sided": float(p_two),
        "mannwhitney_U_NRgtR": float(u_gr), "mannwhitney_p_NRgtR_one_sided": float(p_gr),
        "welch_t": float(t), "welch_p_two_sided": float(p_t),
    }


for tjname, tjscore in [("TJ_core", tj_gse), ("TJ_kegg", tjk_gse)]:
    results["B4"][tjname] = group_test(tjscore)

# ---------------------------------------------------------------------------
# Persist
# ---------------------------------------------------------------------------
with open(OUT / "results_B3B4.json", "w") as f:
    json.dump(results, f, indent=2)

# per-sample score tables (for transparency / plotting)
pd.DataFrame({
    "TJ_core": tj_tcga, "TJ_kegg": tjk_tcga, "CD8": cd8_tcga, "GEP": gep_tcga,
}).to_csv(OUT / "TCGA_LUAD_signature_scores.csv")

pd.DataFrame({
    "response": labels, "TJ_core": tj_gse, "TJ_kegg": tjk_gse,
}).to_csv(OUT / "GSE126044_signature_scores.csv")

print(json.dumps(results, indent=2))
