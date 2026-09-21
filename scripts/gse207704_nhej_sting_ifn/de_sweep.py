#!/usr/bin/env python3
"""Replicate-level DE sweep for GSE207704 CLDN4 knockout.

Counts are kallisto transcript estimates (Ensembl 110 cDNA, SE fragment length
200 ± 30) summed to gene symbols and rounded. The eight GEO runs are 2 vs 2
per line:

  T47D WT  GSM6310640 SRR20029125, GSM6310641 SRR20029124
  T47D KO  GSM6310642 SRR20029123, GSM6310643 SRR20029122
  MCF7 WT  GSM6310644 SRR20029121, GSM6310645 SRR20029120
  MCF7 KO  GSM6310646 SRR20029119, GSM6310647 SRR20029118

Methods, all KO minus WT (positive = higher after CLDN4 loss):
  DESeq2 Wald, fitType parametric with mean-dispersion fallback (inmoose 0.9.1)
  edgeR quasi-likelihood F-test on a common dispersion, TMM norm factors
  limma-trend: lmFit on log2 CPM, eBayes(trend=True). Robust squeezeVar is not in inmoose 0.9.1.

The thesis used to rank rows is pre-declared: NHEJ down; STING, IFN, and APM up.
Every specification is written to the sweep table. Nothing is dropped because
its sign disagrees with that thesis.
"""
from __future__ import annotations

import gzip
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import patsy
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gsea_core import NPERM, SEED, bh_fdr, gsea_prerank  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "methods" / "gse207704_nhej_sting_ifn"
TABLES = OUT / "tables"
KAL = Path("/tmp/gse207704_quant/kal")
FASTA = Path("/tmp/gse207704_quant/Homo_sapiens.GRCh38.cdna.all.fa.gz")
GMT = OUT / "reactome_sets.gmt"
HALLMARK = Path(
    "/home/ubuntu/.cursor/projects/workspace/agent-tools/fb4ac5cf-0521-4111-a8fd-1159ae4ef70b.txt"
)

SAMPLES = [
    ("SRR20029125", "GSM6310640", "T47D", "WT"),
    ("SRR20029124", "GSM6310641", "T47D", "WT"),
    ("SRR20029123", "GSM6310642", "T47D", "KO"),
    ("SRR20029122", "GSM6310643", "T47D", "KO"),
    ("SRR20029121", "GSM6310644", "MCF7", "WT"),
    ("SRR20029120", "GSM6310645", "MCF7", "WT"),
    ("SRR20029119", "GSM6310646", "MCF7", "KO"),
    ("SRR20029118", "GSM6310647", "MCF7", "KO"),
]

NHEJ_CORE = ["PRKDC", "LIG4", "XRCC4", "XRCC5", "XRCC6", "NHEJ1"]
STING_AXIS = ["CGAS", "STING1", "TBK1", "IRF3"]
IFN_GENES = [
    "IFI27", "OAS2", "IFIT1", "MX1", "ISG15", "IFIT2", "IFIT3", "IFIT5", "MX2",
    "OAS1", "OAS3", "OASL", "RSAD2", "USP18", "IFI6", "IFI35", "IFI44", "IFI44L",
    "IFI16", "BST2", "IFITM1", "IFITM2", "IFITM3", "DDX58", "IFIH1", "DDX60",
    "HERC5", "XAF1", "EIF2AK2", "LY6E", "CMPK2", "EPSTI1", "SAMD9L", "TRIM22",
    "ZBP1", "PLSCR1", "STAT1", "STAT2", "IRF1", "IRF7", "IRF9", "JAK1", "JAK2",
    "IFNAR1", "IFNAR2", "IFNGR1", "IFNGR2", "SOCS1", "SOCS3", "GBP1", "GBP2",
    "GBP4", "GBP5", "CXCL9", "CXCL10", "CXCL11", "CIITA", "IDO1", "IFNG", "IFNB1",
]
APM_GENES = [
    "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G", "B2M", "NLRC5",
    "TAP1", "TAP2", "TAPBP", "PSMB8", "PSMB9", "PSMB10", "ERAP1", "ERAP2",
    "CALR", "CANX", "PDIA3", "SEC61A1", "PSME1", "PSME2",
]
ALIAS = {"CGAS": "MB21D1", "STING1": "TMEM173", "H2AX": "H2AFX", "RIGI": "DDX58"}
# Pre-declared thesis. Used only to flag rows, never to edit statistics.
THESIS = {}
CUSTOM = {
    "NHEJ_CORE": NHEJ_CORE,
    "STING_AXIS": STING_AXIS,
    "IFN": IFN_GENES,
    "APM": APM_GENES,
}
for name in (
    "NHEJ_CORE",
    "REACTOME_NONHOMOLOGOUS_END_JOINING",
):
    THESIS[name] = "down"
for name in (
    "STING_AXIS",
    "REACTOME_STING_MEDIATED_INDUCTION",
    "IFN",
    "APM",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING",
    "REACTOME_INTERFERON_GAMMA_SIGNALING",
    "REACTOME_CLASS_I_MHC_PEPTIDE_LOADING",
):
    THESIS[name] = "up"

FLOORS = (0.0, 0.5, 1.0)
NAMED = ["CLDN4", "TACSTD2", *NHEJ_CORE, *STING_AXIS]


def tx_to_gene(fasta: Path) -> dict[str, str]:
    mapping = {}
    opener = gzip.open if str(fasta).endswith(".gz") else open
    with opener(fasta, "rt") as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            parts = line[1:].strip().split()
            tid = parts[0]
            symbol = ""
            for tok in parts[1:]:
                if tok.startswith("gene_symbol:"):
                    symbol = tok.split(":", 1)[1]
            if symbol and symbol != ".":
                mapping[tid] = symbol
    return mapping


def gene_counts_from_kallisto() -> pd.DataFrame:
    print("building gene counts from kallisto")
    txmap = tx_to_gene(FASTA)
    frames = []
    for srr, gsm, line, geno in SAMPLES:
        path = KAL / srr / "abundance.tsv"
        if not path.exists():
            raise FileNotFoundError(path)
        ab = pd.read_csv(path, sep="\t")
        ab["gene"] = ab["target_id"].map(txmap)
        ab = ab.dropna(subset=["gene"])
        summed = ab.groupby("gene")["est_counts"].sum()
        summed.name = srr
        frames.append(summed)
        info = (KAL / srr / "run_info.json").read_text()
        print(srr, gsm, line, geno, "genes", summed.shape[0], info.strip().replace("\n", " ")[:180])
    mat = pd.concat(frames, axis=1).fillna(0.0)
    mat = mat.round().astype(int)
    mat.index.name = "gene"
    return mat


def load_counts() -> pd.DataFrame:
    path = TABLES / "kallisto_gene_counts.tsv.gz"
    if path.exists():
        mat = pd.read_csv(path, sep="\t", index_col=0)
        print("loaded", path, mat.shape)
        return mat
    mat = gene_counts_from_kallisto()
    TABLES.mkdir(parents=True, exist_ok=True)
    mat.to_csv(path, sep="\t")
    print("wrote", path, mat.shape)
    return mat


def tmm_factors(counts: np.ndarray) -> np.ndarray:
    """Robinson-Oshlack TMM norm factors. counts is genes x samples."""
    lib = counts.sum(axis=0).astype(float)
    q = np.quantile(counts, 0.75, axis=0)
    ref = int(np.argmin(np.abs(q - q.mean())))
    factors = np.ones(counts.shape[1], dtype=float)
    for j in range(counts.shape[1]):
        if j == ref:
            continue
        y = counts[:, j].astype(float)
        r = counts[:, ref].astype(float)
        ok = (y > 0) & (r > 0)
        if int(ok.sum()) < 50:
            continue
        m = np.log2((y[ok] / lib[j]) / (r[ok] / lib[ref]))
        a = 0.5 * np.log2((y[ok] / lib[j]) * (r[ok] / lib[ref]))
        m_lo, m_hi = np.quantile(m, [0.30, 0.70])
        a_lo, a_hi = np.quantile(a, [0.05, 0.95])
        keep = (m >= m_lo) & (m <= m_hi) & (a >= a_lo) & (a <= a_hi)
        if int(keep.sum()) < 20:
            keep = np.ones(m.shape[0], dtype=bool)
        n1 = y[ok][keep]
        n2 = r[ok][keep]
        var = (lib[j] - n1) / lib[j] / n1 + (lib[ref] - n2) / lib[ref] / n2
        w = np.zeros_like(var)
        good = np.isfinite(var) & (var > 0)
        w[good] = 1.0 / var[good]
        mm = m[keep]
        factors[j] = 2 ** (np.sum(w * mm) / w.sum()) if w.sum() else float(np.mean(mm))
    factors /= math.exp(np.mean(np.log(factors)))
    return factors


def log_cpm(counts: pd.DataFrame, norm: np.ndarray, prior: float = 3.0) -> pd.DataFrame:
    lib = counts.sum(axis=0).to_numpy(dtype=float) * norm
    prior_scaled = prior * lib / lib.mean()
    num = counts.to_numpy(dtype=float) + prior_scaled
    den = lib + 2 * prior_scaled
    logged = np.log2(num / den * 1e6)
    return pd.DataFrame(logged, index=counts.index, columns=counts.columns)


def filter_counts(counts: pd.DataFrame, min_count: int = 10, min_samples: int = 2) -> pd.DataFrame:
    keep = (counts >= min_count).sum(axis=1) >= min_samples
    return counts.loc[keep]


def design_matrix(meta: pd.DataFrame):
    """Patsy DesignMatrix aligned to meta's row order, plus the KO-vs-WT column."""
    obs = meta.copy()
    obs["line"] = pd.Categorical(obs["line"], categories=["T47D", "MCF7"])
    obs["genotype"] = pd.Categorical(obs["genotype"], categories=["WT", "KO"])
    formula = "~ genotype" if obs["line"].nunique() == 1 else "~ line + genotype"
    design = patsy.dmatrix(formula, obs)
    coef = [c for c in design.design_info.column_names if "genotype" in c and "KO" in c]
    if len(coef) != 1:
        raise RuntimeError(f"KO coefficient not unique: {design.design_info.column_names}")
    return design, coef[0]


def run_deseq2(counts: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    from inmoose.deseq2 import DESeq, DESeqDataSet

    obs = meta.copy()
    obs["line"] = pd.Categorical(obs["line"], categories=["T47D", "MCF7"])
    obs["genotype"] = pd.Categorical(obs["genotype"], categories=["WT", "KO"])
    formula = "~ genotype" if obs["line"].nunique() == 1 else "~ line + genotype"
    # inmoose wants samples x genes
    cdf = counts.T.loc[obs.index]
    dds = DESeqDataSet(cdf, clinicalData=obs, design=formula)
    try:
        dds = DESeq(dds, fitType="parametric", quiet=True)
        fit = "parametric"
    except Exception as exc:  # parametric curve is not always implemented
        print("DESeq2 parametric failed, using fitType=mean:", type(exc).__name__, exc)
        dds = DESeqDataSet(cdf, clinicalData=obs, design=formula)
        dds = DESeq(dds, fitType="mean", quiet=True)
        fit = "mean"
    res = dds.results(name="genotype_KO_vs_WT", independentFiltering=False)
    out = pd.DataFrame(
        {
            "log2FC": res["log2FoldChange"].astype(float),
            "stat": res["stat"].astype(float),
            "pvalue": res["pvalue"].astype(float),
            "baseMean": res["baseMean"].astype(float),
        },
        index=res.index,
    )
    out["fdr"] = bh_fdr(out["pvalue"])
    out["fit"] = fit
    return out


def run_edger(counts: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    from inmoose.edgepy import DGEList, glmQLFTest

    meta = meta.loc[list(counts.columns)]
    design, coef_name = design_matrix(meta)
    y = DGEList(counts)
    y.samples["norm_factors"] = tmm_factors(counts.to_numpy())
    y = y.estimateGLMCommonDisp(design)
    fit = y.glmQLFit(design, robust=False)
    coef_idx = list(design.design_info.column_names).index(coef_name)
    qlf = glmQLFTest(fit, coef=coef_idx)
    logfc_col = "log2FoldChange" if "log2FoldChange" in qlf.columns else "logFC"
    signed = qlf["stat"].astype(float).to_numpy() * np.sign(qlf[logfc_col].astype(float).to_numpy())
    out = pd.DataFrame(
        {
            "log2FC": qlf[logfc_col].astype(float).to_numpy(),
            "stat": signed,
            "pvalue": qlf["pvalue"].astype(float).to_numpy(),
        },
        index=qlf.index,
    )
    out["fdr"] = bh_fdr(out["pvalue"])
    out["fit"] = f"QLF_common_disp_{y.common_dispersion:.4g}"
    return out


def run_limma(counts: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    from inmoose.limma import eBayes, lmFit, topTable

    meta = meta.loc[list(counts.columns)]
    design, coef_name = design_matrix(meta)
    norm = tmm_factors(counts.to_numpy())
    logged = log_cpm(counts, norm)
    fit = lmFit(logged, design)
    fit = eBayes(fit, trend=True, robust=False)
    # inmoose 0.9.1 rejects sort_by="none" even though the docstring lists it.
    tt = topTable(fit, coef=coef_name, number=logged.shape[0], adjust_method="fdr_bh", sort_by="P")
    # topTable may reorder. Reindex to genes.
    cols = {c.lower(): c for c in tt.columns}
    log_col = cols.get("logfc", cols.get("log2fc", cols.get("log2foldchange")))
    t_col = cols.get("t", cols.get("stat"))
    p_col = cols.get("p_value", cols.get("pvalue", cols.get("p.value")))
    if log_col is None or t_col is None or p_col is None:
        raise RuntimeError(f"unexpected limma columns {list(tt.columns)}")
    out = pd.DataFrame(
        {
            "log2FC": tt[log_col].astype(float),
            "stat": tt[t_col].astype(float),
            "pvalue": tt[p_col].astype(float),
        },
        index=tt.index,
    )
    out = out.reindex(logged.index)
    out["fdr"] = bh_fdr(out["pvalue"])
    out["fit"] = "limma_trend"
    return out


def read_sets(symbols: pd.Index) -> dict[str, list[str]]:
    sets = {k: list(v) for k, v in CUSTOM.items()}
    for line in GMT.read_text().splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        sets[parts[0]] = parts[2:]
    if HALLMARK.exists():
        for line in HALLMARK.read_text().splitlines():
            if line.startswith("Interferon Alpha Response "):
                sets["HALLMARK_INTERFERON_ALPHA_RESPONSE"] = line.split()[3:]
            elif line.startswith("Interferon Gamma Response "):
                sets["HALLMARK_INTERFERON_GAMMA_RESPONSE"] = line.split()[3:]
    resolved = {}
    index = set(symbols)
    for name, genes in sets.items():
        hits = []
        for g in genes:
            if g in index and g not in hits:
                hits.append(g)
            else:
                alt = ALIAS.get(g)
                if alt and alt in index and alt not in hits:
                    hits.append(alt)
        resolved[name] = hits
    return resolved


def set_tests(de: pd.DataFrame, sets: dict[str, list[str]], floor: float) -> list[dict]:
    base = de.replace([np.inf, -np.inf], np.nan).dropna(subset=["log2FC", "stat"])
    if floor > 0:
        base = base.loc[base["log2FC"].abs() >= floor]
    rows = []
    if base.empty:
        return rows
    rank = base["stat"].sort_values(ascending=False)
    print(f"  set tests floor {floor}: {len(base)} genes", flush=True)
    # GSEA once per floor, all sets.
    min_floor = 3
    gsea = gsea_prerank(rank, sets, nperm=NPERM, seed=SEED, min_size=min_floor, max_size=500)
    gsea_map = {r.term: r for r in gsea.itertuples()} if len(gsea) else {}
    rng = np.random.default_rng(SEED)
    values = rank.to_numpy()
    names = rank.index.to_numpy()
    pos = {g: i for i, g in enumerate(names)}
    n = len(names)
    for term, members in sets.items():
        idx = [pos[g] for g in members if g in pos]
        idx = np.unique(np.array(idx, dtype=int))
        mean_lfc = float(base.loc[[names[i] for i in idx], "log2FC"].mean()) if len(idx) else np.nan
        rec = {
            "term": term,
            "n_set_in_test": int(len(idx)),
            "mean_log2FC": mean_lfc,
            "lfc_floor": floor,
        }
        if term in gsea_map:
            g = gsea_map[term]
            rec.update(
                {
                    "gsea_nes": float(g.nes),
                    "gsea_p": float(g.nom_p),
                    "gsea_lead": g.lead_genes,
                }
            )
        else:
            rec.update({"gsea_nes": np.nan, "gsea_p": np.nan, "gsea_lead": ""})
        if len(idx) >= 3:
            obs = float(values[idx].mean())
            null = np.empty(NPERM)
            for i in range(NPERM):
                pick = rng.choice(n, size=len(idx), replace=False)
                null[i] = values[pick].mean()
            rec["camera_stat"] = obs
            rec["camera_p"] = float((np.sum(np.abs(null) >= abs(obs)) + 1) / (NPERM + 1))
            lfc = base.loc[[names[i] for i in idx], "log2FC"].to_numpy()
            # Wilcoxon needs non-zero differences. Skip if all ~0.
            if np.allclose(lfc, 0):
                rec["wilcox_p"] = np.nan
                rec["wilcox_thesis_p"] = np.nan
            else:
                rec["wilcox_p"] = float(stats.wilcoxon(lfc, alternative="two-sided").pvalue)
                alt = "less" if THESIS.get(term) == "down" else "greater"
                rec["wilcox_thesis_p"] = float(stats.wilcoxon(lfc, alternative=alt).pvalue)
        else:
            rec["camera_stat"] = np.nan
            rec["camera_p"] = np.nan
            rec["wilcox_p"] = np.nan
            rec["wilcox_thesis_p"] = np.nan
        rows.append(rec)
    return rows


def long_sweep(raw_rows: list[dict]) -> pd.DataFrame:
    """One row per test so FDR vs nominal can be compared without hiding signs."""
    long_rows = []
    for r in raw_rows:
        direction = "up" if r["mean_log2FC"] > 0 else "down" if r["mean_log2FC"] < 0 else "zero"
        thesis = THESIS.get(r["term"], "")
        match = bool(thesis) and direction == thesis
        specs = [
            ("gsea", r["gsea_nes"], r["gsea_p"], "two-sided gene-set permutation"),
            ("mean_stat_perm", r["camera_stat"], r["camera_p"], "two-sided gene-label permutation of mean statistic"),
            ("wilcoxon_log2FC", r["mean_log2FC"], r["wilcox_p"], "two-sided Wilcoxon of gene log2FC vs 0"),
            ("wilcoxon_thesis_sided", r["mean_log2FC"], r["wilcox_thesis_p"], "one-sided Wilcoxon in the thesis direction"),
        ]
        for test, effect, p, note in specs:
            long_rows.append(
                {
                    "method": r["method"],
                    "contrast": r["contrast"],
                    "term": r["term"],
                    "thesis": thesis,
                    "direction": direction,
                    "thesis_sign_match": match,
                    "test": test,
                    "effect": effect,
                    "mean_log2FC": r["mean_log2FC"],
                    "pvalue": p,
                    "n_set_in_test": r["n_set_in_test"],
                    "lfc_floor": r["lfc_floor"],
                    "family": "primary" if r["lfc_floor"] == 0 and test != "wilcoxon_thesis_sided" else "sensitivity",
                    "note": note,
                    "gsea_lead": r["gsea_lead"] if test == "gsea" else "",
                    "de_fit": r["fit"],
                }
            )
    out = pd.DataFrame(long_rows)
    out["fdr_within_family"] = np.nan
    for _, idx in out.groupby(["method", "contrast", "test", "lfc_floor"]).groups.items():
        out.loc[idx, "fdr_within_family"] = bh_fdr(out.loc[idx, "pvalue"]).to_numpy()
    out["fdr_global"] = bh_fdr(out["pvalue"])
    return out


def plot(gene_de: pd.DataFrame, sweep: pd.DataFrame) -> None:
    fig_dir = OUT / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    # Named-gene log2FC from the three methods, two lines.
    sub = gene_de[gene_de["contrast"].isin(["T47D", "MCF7"]) & gene_de["gene"].isin(NAMED)].copy()
    methods = ["DESeq2", "edgeR_QL", "limma_trend"]
    genes = [g for g in NAMED if g in set(sub["gene"])]
    mats = []
    for contrast in ["T47D", "MCF7"]:
        chunk = sub[sub["contrast"].eq(contrast)]
        mat = np.full((len(genes), len(methods)), np.nan)
        for j, method in enumerate(methods):
            m = chunk[chunk["method"].eq(method)].set_index("gene")
            for i, gene in enumerate(genes):
                if gene in m.index:
                    mat[i, j] = m.loc[gene, "log2FC"]
        mats.append(mat)
    finite = np.concatenate([m[np.isfinite(m)] for m in mats])
    vmax = max(1.0, float(np.max(np.abs(finite))) if finite.size else 1.0)
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 7.4), sharey=True, constrained_layout=True)
    for ax, contrast, mat in zip(axes, ["T47D", "MCF7"], mats):
        im = ax.imshow(mat, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
        ax.set_xticks(range(len(methods)), ["DESeq2", "edgeR", "limma"], fontsize=9)
        ax.set_yticks(range(len(genes)), genes, fontsize=9)
        ax.set_title(f"{contrast} KO vs WT")
        ax.tick_params(length=0)
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                val = mat[i, j]
                if np.isfinite(val):
                    ax.text(j, i, f"{val:+.2f}", ha="center", va="center", fontsize=8, color="black")
    fig.colorbar(im, ax=axes, fraction=0.03, pad=0.02, label="log2FC")
    fig.suptitle("GSE207704 replicate-level log2FC (kallisto counts)")
    fig.savefig(fig_dir / "fig_replicate_de_named_genes.png", dpi=160)
    fig.savefig(fig_dir / "fig_replicate_de_named_genes.pdf")
    plt.close(fig)

    # Primary GSEA sweep: lfc floor 0, NES. Shared color scale across methods.
    gsea = sweep[(sweep["test"].eq("gsea")) & (sweep["lfc_floor"].eq(0))].copy()
    terms = [t for t in THESIS if t in set(gsea["term"])]
    contrasts = ["T47D", "MCF7", "pooled"]
    short = [
        t.replace("REACTOME_", "R ")
        .replace("HALLMARK_", "H ")
        .replace("INTERFERON_", "IFN ")
        .replace("NONHOMOLOGOUS_END_JOINING", "NHEJ")
        .replace("CLASS_I_MHC_PEPTIDE_LOADING", "MHC-I loading")
        .replace("STING_MEDIATED_INDUCTION", "STING induction")
        .replace("_", " ")
        for t in terms
    ]
    mats = []
    for method in methods:
        mat = np.full((len(terms), len(contrasts)), np.nan)
        chunk = gsea[gsea["method"].eq(method)]
        for j, contrast in enumerate(contrasts):
            m = chunk[chunk["contrast"].eq(contrast)].set_index("term")
            for i, term in enumerate(terms):
                if term in m.index and np.isfinite(m.loc[term, "effect"]):
                    mat[i, j] = m.loc[term, "effect"]
        mats.append(mat)
    finite = np.concatenate([m[np.isfinite(m)] for m in mats])
    vmax = float(np.max(np.abs(finite))) if finite.size else 1.0
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 7.6), sharey=True, constrained_layout=True)
    for ax, method, mat in zip(axes, methods, mats):
        im = ax.imshow(mat, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
        ax.set_xticks(range(len(contrasts)), contrasts, fontsize=9)
        ax.set_title(method.replace("_", " "))
        ax.tick_params(length=0)
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                val = mat[i, j]
                if np.isfinite(val):
                    ax.text(j, i, f"{val:+.2f}", ha="center", va="center", fontsize=8, color="black")
    axes[0].set_yticks(range(len(terms)), short, fontsize=8)
    fig.colorbar(im, ax=axes, fraction=0.03, pad=0.02, label="GSEA NES (positive = up after KO)")
    fig.suptitle("Primary prerank GSEA, no effect-size floor")
    fig.savefig(fig_dir / "fig_de_sweep_gsea_nes.png", dpi=160)
    fig.savefig(fig_dir / "fig_de_sweep_gsea_nes.pdf")
    plt.close(fig)

    # Mean log2FC of each pre-declared set. One point per DE method.
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 7.2), sharey=True, constrained_layout=True)
    colors = {"DESeq2": "#1b4f72", "edgeR_QL": "#b9770e", "limma_trend": "#1e8449"}
    offsets = {"DESeq2": -0.18, "edgeR_QL": 0.0, "limma_trend": 0.18}
    for ax, contrast in zip(axes, contrasts):
        chunk = gsea[gsea["contrast"].eq(contrast)]
        ax.axvline(0, color="#888888", lw=0.8, zorder=0)
        for method in methods:
            m = chunk[chunk["method"].eq(method)].set_index("term")
            ys = np.arange(len(terms)) + offsets[method]
            xs = [m.loc[t, "mean_log2FC"] if t in m.index else np.nan for t in terms]
            ax.scatter(xs, ys, s=28, color=colors[method], label=method.replace("_", " "), zorder=2)
        ax.set_yticks(range(len(terms)), short, fontsize=8)
        ax.set_title(contrast)
        ax.set_xlabel("mean log2FC")
        ax.set_xlim(-0.55, 0.55)
    axes[2].legend(frameon=False, fontsize=8, loc="lower right")
    fig.suptitle("Set mean log2FC after CLDN4 knockout (floor 0)")
    fig.savefig(fig_dir / "fig_de_sweep_mean_log2fc.png", dpi=160)
    fig.savefig(fig_dir / "fig_de_sweep_mean_log2fc.pdf")
    plt.close(fig)


def write_report(gene_de: pd.DataFrame, sweep: pd.DataFrame, counts: pd.DataFrame) -> None:
    primary = sweep[(sweep["family"].eq("primary")) & (sweep["test"].eq("gsea"))].copy()
    matched = primary[primary["thesis_sign_match"].eq(True) & primary["pvalue"].notna()]
    opposed = primary[primary["thesis_sign_match"].eq(False) & primary["pvalue"].notna()]

    def best(frame: pd.DataFrame) -> pd.Series | None:
        if frame.empty:
            return None
        return frame.sort_values(["pvalue", "fdr_within_family"]).iloc[0]

    best_match = best(matched)
    best_any = best(primary)
    best_opp = best(opposed)
    # Sensitivity: smallest nominal p among thesis-matching rows anywhere in the sweep.
    sens = sweep[sweep["thesis_sign_match"].eq(True) & sweep["pvalue"].notna()]
    best_sens = best(sens)

    def fmt_row(r: pd.Series | None) -> str:
        if r is None:
            return "none"
        return (
            f"{r['method']} / {r['contrast']} / {r['term']} / {r['test']} / floor {r['lfc_floor']}: "
            f"mean log2FC {r['mean_log2FC']:+.3f}, effect {r['effect']:+.3f}, "
            f"nominal p {r['pvalue']:.4g}, within-family FDR {r['fdr_within_family']:.4g}, "
            f"global FDR {r['fdr_global']:.4g}, n={int(r['n_set_in_test'])}"
        )

    n_match_fdr = int((sweep["thesis_sign_match"].eq(True) & (sweep["fdr_within_family"] < 0.05)).sum())
    n_match_global = int((sweep["thesis_sign_match"].eq(True) & (sweep["fdr_global"] < 0.05)).sum())
    n_rows = int(len(sweep))
    dir_lines = [
        "| Set | thesis | fits with that sign | median mean log2FC | smallest nominal p (any sign) | best thesis-match p | smallest within-family FDR (any sign) |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for term, thesis in THESIS.items():
        sub = primary[primary["term"].eq(term)]
        if sub.empty:
            continue
        n_fit = int(len(sub))
        n_ok = int(sub["thesis_sign_match"].sum())
        med = float(sub["mean_log2FC"].median())
        best_p = float(sub["pvalue"].min())
        matched_term = sub[sub["thesis_sign_match"].eq(True)]
        best_match_p = float(matched_term["pvalue"].min()) if len(matched_term) else float("nan")
        best_fdr = float(sub["fdr_within_family"].min())
        match_txt = f"{best_match_p:.4g}" if np.isfinite(best_match_p) else "no row"
        dir_lines.append(
            f"| {term} | {thesis} | {n_ok}/{n_fit} | {med:+.3f} | {best_p:.4g} | {match_txt} | {best_fdr:.4g} |"
        )
    direction_table = "\n".join(dir_lines)

    def gene_md(contrast: str, method: str) -> str:
        sub = gene_de[(gene_de["contrast"].eq(contrast)) & (gene_de["method"].eq(method))]
        sub = sub.set_index("gene")
        lines = ["| Gene | log2FC | nominal p | FDR |", "|---|---:|---:|---:|"]
        for gene in NAMED:
            if gene not in sub.index:
                lines.append(f"| {gene} | absent |  |  |")
                continue
            r = sub.loc[gene]
            if isinstance(r, pd.DataFrame):
                r = r.iloc[0]
            lines.append(
                f"| {gene} | {r['log2FC']:+.3f} | {r['pvalue']:.3g} | {r['fdr']:.3g} |"
            )
        return "\n".join(lines)

    cldn = gene_de[gene_de["gene"].eq("CLDN4")]
    cldn_bits = ", ".join(
        f"{r.method} {r.contrast} {r.log2FC:+.3f} (p {r.pvalue:.3g})"
        for r in cldn.itertuples()
    )
    text = f"""# FINDING — GSE207704 replicate-level DE sweep

Additive counts analysis. The earlier FPKM scores stand. This file uses the eight deposited runs, not the collapsed FPKM table.

## Design that was actually fit

Kallisto 0.51.1, Ensembl 110 cDNA, single-end, fragment length 200 (sd 30). Transcript estimates were summed to gene symbols and rounded. n = 2 per genotype per line. Positive log2FC = higher in CLDN4−/− than in parental WT.

Methods: DESeq2 Wald (inmoose 0.9.1), edgeR quasi-likelihood F on a common dispersion with TMM factors, limma-trend (`eBayes(trend=True)` on log2 CPM). Robust empirical-Bayes squeezing is not implemented in this inmoose build, so both edgeR QL and limma use the standard squeeze. Genes with fewer than 2 samples at count ≥ 10 were dropped before each fit. FDR within a method is Benjamini-Hochberg on the genes tested in that contrast.

CLDN4 across fits: {cldn_bits}.

Library sizes are in `tables/kallisto_sample_qc.tsv`. Gene counts: `tables/kallisto_gene_counts.tsv.gz` ({counts.shape[0]} genes × {counts.shape[1]} samples).

## Pre-declared thesis and the rule for “best”

Thesis, fixed before looking at these p-values: NHEJ sets go **down** after CLDN4 loss; STING, IFN, and APM sets go **up**. A row matches when the set’s mean log2FC has that sign.

Primary family: GSEA on the signed DE statistic, effect-size floor 0, two-sided gene-set permutation, three methods × three contrasts (T47D, MCF7, pooled `line + genotype`). The best primary match is the smallest nominal p inside that family among sign-matching rows. Sensitivity rows (floors 0.5 and 1, one-sided Wilcoxon, mean-statistic permutation) are in the same table and are not substituted for the primary result.

Global BH-FDR is computed across every test in the sweep, because the sweep itself is many looks.

## What is most significant

The strongest primary result opposes the thesis. IFN sets go down after CLDN4 loss. No specification in this sweep — method, contrast, test, or effect-size floor — puts a thesis-matching set at within-family FDR < 0.05. Of {n_rows} tests, {n_match_fdr} thesis-matching rows have within-family FDR < 0.05 and {n_match_global} have global FDR < 0.05.

Gene-set permutation p-values cannot be smaller than 1/1001 with 1000 permutations.

- Strongest primary GSEA row, any direction: {fmt_row(best_any)}
- Strongest primary GSEA row whose sign matches the thesis: {fmt_row(best_match)}
- Strongest primary GSEA row whose sign opposes the thesis: {fmt_row(best_opp)}
- Strongest thesis-matching row anywhere in the sensitivity sweep: {fmt_row(best_sens)}

The sensitivity winner is an effect-size floor that keeps only genes with |log2FC| ≥ 0.5. It is not the primary result.

## Primary GSEA direction (floor 0)

Nine fits = DESeq2, edgeR QL, and limma-trend, each on T47D, MCF7, and the line-adjusted pool. The smallest nominal p and the smallest within-family FDR are taken over those nine whatever the sign. A small FDR on a down IFN set is evidence against the thesis.

{direction_table}

NHEJ_CORE mean log2FC is positive in every fit, so the 6-gene panel does not go down. Reactome NHEJ is slightly negative in every fit and is not significant. Reactome STING-mediated induction contains PRKDC, XRCC5, XRCC6, and MRE11, so it is not a pure CGAS–STING1–TBK1–IRF3 test. STING1 itself is present in the cDNA counts but the counts are small (see the QC table); the positive log2FC is not a significant gene-level result.

IFN, Reactome IFN-α/β, and Hallmark IFN-γ are negative in every primary fit. Hallmark IFN-α is negative in the fits that reach nominal p < 0.05.

## Named genes, DESeq2

### T47D

{gene_md("T47D", "DESeq2")}

### MCF7

{gene_md("MCF7", "DESeq2")}

edgeR and limma log2FC values are in `tables/gene_de_replicate.tsv` and in `figures/fig_replicate_de_named_genes.png`. Primary GSEA NES values are in `figures/fig_de_sweep_gsea_nes.png`. Set mean log2FC is in `figures/fig_de_sweep_mean_log2fc.png`. The full grid is `tables/de_sweep.tsv`.
"""
    (OUT / "SWEEP.md").write_text(text)
    print(text)


def main() -> None:
    counts = load_counts()
    counts = counts.loc[counts.sum(axis=1) > 0]
    meta = pd.DataFrame(SAMPLES, columns=["srr", "gsm", "line", "genotype"]).set_index("srr")
    meta = meta.loc[counts.columns]
    sets = read_sets(counts.index)
    print({k: len(v) for k, v in sets.items()})

    qc_rows = []
    for srr, gsm, line, geno in SAMPLES:
        info_path = KAL / srr / "run_info.json"
        n_proc = n_pseudo = ""
        if info_path.exists():
            import json
            info = json.loads(info_path.read_text())
            n_proc = info.get("n_processed", "")
            n_pseudo = info.get("n_pseudoaligned", "")
        qc_rows.append(
            {
                "srr": srr,
                "gsm": gsm,
                "line": line,
                "genotype": geno,
                "n_processed": n_proc,
                "n_pseudoaligned": n_pseudo,
                "gene_count_sum": int(counts[srr].sum()),
                "CLDN4": int(counts.loc["CLDN4", srr]) if "CLDN4" in counts.index else 0,
                "STING1": int(counts.loc["STING1", srr]) if "STING1" in counts.index else 0,
                "CGAS": int(counts.loc["CGAS", srr]) if "CGAS" in counts.index else 0,
            }
        )
    TABLES.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(qc_rows).to_csv(TABLES / "kallisto_sample_qc.tsv", sep="\t", index=False)

    runners = {
        "DESeq2": run_deseq2,
        "edgeR_QL": run_edger,
        "limma_trend": run_limma,
    }
    contrasts = {}
    for name, lines in (
        ("T47D", ["T47D"]),
        ("MCF7", ["MCF7"]),
        ("pooled", ["T47D", "MCF7"]),
    ):
        cols = [s for s in counts.columns if meta.loc[s, "line"] in lines]
        contrasts[name] = cols

    gene_rows = []
    sweep_raw = []
    for contrast, cols in contrasts.items():
        sub_counts = filter_counts(counts[cols])
        sub_meta = meta.loc[cols]
        print(f"contrast {contrast}: {sub_counts.shape[0]} genes, samples {cols}")
        for method, fn in runners.items():
            print(" running", method, contrast)
            de = fn(sub_counts, sub_meta)
            de = de.replace([np.inf, -np.inf], np.nan).dropna(subset=["log2FC", "stat", "pvalue"])
            for gene in NAMED:
                if gene not in de.index and gene in ALIAS and ALIAS[gene] in de.index:
                    src = ALIAS[gene]
                elif gene in de.index:
                    src = gene
                else:
                    continue
                r = de.loc[src]
                gene_rows.append(
                    {
                        "gene": gene,
                        "deposited_symbol": src,
                        "method": method,
                        "contrast": contrast,
                        "log2FC": float(r["log2FC"]),
                        "stat": float(r["stat"]),
                        "pvalue": float(r["pvalue"]),
                        "fdr": float(r["fdr"]),
                        "fit": r["fit"],
                    }
                )
            for floor in FLOORS:
                for rec in set_tests(de, sets, floor):
                    rec["method"] = method
                    rec["contrast"] = contrast
                    rec["fit"] = str(de["fit"].iloc[0])
                    sweep_raw.append(rec)

    gene_de = pd.DataFrame(gene_rows)
    sweep = long_sweep(sweep_raw)
    gene_de.to_csv(TABLES / "gene_de_replicate.tsv", sep="\t", index=False, float_format="%.6g")
    sweep.to_csv(TABLES / "de_sweep.tsv", sep="\t", index=False, float_format="%.6g")
    # CLDN4 must fall in both lines. A positive CLDN4 logFC means the count matrix is mislabeled.
    cldn = gene_de[gene_de["gene"].eq("CLDN4") & gene_de["contrast"].isin(["T47D", "MCF7"])]
    if (cldn["log2FC"] > 0).any():
        raise SystemExit(f"CLDN4 did not fall after knockout:\n{cldn}")
    plot(gene_de, sweep)
    write_report(gene_de, sweep, counts)


if __name__ == "__main__":
    main()
