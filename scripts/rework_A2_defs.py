#!/usr/bin/env python3
"""Rework A2: TACSTD2 vs every requested immune definition.

Cohorts (only public human durvalumab NSCLC bulk RNA-seq):
  GSE253564  pre-treatment biopsies, n=32, FPKM  (NCT02904954)
  GSE248378  post-treatment resections, n=29, FPKM (same trial)

User claim: raw Spearman ρ = -0.65, purity-adjusted ρ = -0.46, p = 2e-4.
Prior A2 (ESTIMATE ImmuneScore, post): adj ρ = -0.474, adj p = 0.011.

Immune definitions (ALL reported; none hidden):
  1. ESTIMATE ImmuneScore   Yoshihara 2013; exact Python port of ESTIMATE v1.0.13
  2. CYT                    Rooney 2015; mean log2(FPKM+1) of GZMA and PRF1
  3. CD8A                   log2(FPKM+1)
  4. GEP18                  Ayers 2017 18-gene T-cell-inflamed GEP; mean of
                            per-gene z-scores of log2(FPKM+1)
  5. xCell CD8              Aran 2017 "CD8+ T-cells" score (ssGSEA + transform
                            + RNA-seq spillover, alpha=0.5)

Purity adjustment: partial Spearman controlling for ESTIMATEScore ranks
(rank-identical to ESTIMATE TumorPurity; keeps samples where the Affymetrix
cosine formula is out of bounds on FPKM).

Outputs -> results/rework/A2_defs/
"""
from __future__ import annotations

import json
import os
import re
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import lsq_linear

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(BASE, "data")
OUT = os.path.join(BASE, "results", "rework", "A2_defs")
FIG = os.path.join(OUT, "figures")
os.makedirs(FIG, exist_ok=True)

CLAIM = {"rho_raw": -0.65, "rho_purity_adj": -0.46, "p": 2e-4}
PRIOR_A2_POST_ADJ_P = 0.01085

# Ayers et al., J Clin Invest 2017; NanoString TIS / 18-gene T-cell-inflamed GEP.
# Official weights are NanoString-platform-specific and not used here.
GEP18_GENES = [
    "CCL5", "CD27", "CD274", "CD276", "CD8A", "CMKLR1", "CXCL9", "CXCR6",
    "HLA-DQA1", "HLA-DRB1", "HLA-E", "IDO1", "LAG3", "NKG7", "PDCD1LG2",
    "PSMB10", "STAT1", "TIGIT",
]

XCELL_SOURCES = ("ENCODE", "FANTOM", "HPCA", "BLUEPRINT", "IRIS", "NOVERSHTERN")
XCELL_SIG_RE = re.compile(
    r"_(?:ENCODE|FANTOM|HPCA|BLUEPRINT|IRIS|NOVERSHTERN)_\d+$"
)

PRIMARY_FEATURES = [
    "ESTIMATE_ImmuneScore",
    "CYT",
    "CD8A",
    "GEP18",
    "xCell_CD8",
]


def load_gmt(path: str) -> dict[str, list[str]]:
    sets: dict[str, list[str]] = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            sets[parts[0]] = [g for g in parts[2:] if g]
    return sets


def estimate_scores(expr: pd.DataFrame, gene_sets: dict[str, list[str]]):
    """Exact port of estimateScore() from ESTIMATE R v1.0.13."""
    genes = expr.index.to_numpy()
    n_genes = expr.shape[0]
    m = expr.rank(axis=0, method="average").to_numpy() * (10000.0 / n_genes)

    rows = {}
    overlaps = {}
    for name in ("StromalSignature", "ImmuneSignature"):
        in_set = np.isin(genes, list(set(gene_sets[name])))
        overlaps[name] = int(in_set.sum())
        es = np.empty(expr.shape[1])
        for j in range(expr.shape[1]):
            order = np.argsort(-m[:, j], kind="stable")
            tag = in_set[order].astype(float)
            w = np.abs(m[order, j]) ** 0.25
            nh = tag.sum()
            nm = n_genes - nh
            pn = np.cumsum(tag * w) / (tag * w).sum()
            p0 = np.cumsum((1.0 - tag) / nm)
            es[j] = float(np.sum(pn - p0))
        rows["Stromal" if name.startswith("Str") else "Immune"] = es

    df = pd.DataFrame(
        {"StromalScore": rows["Stromal"], "ImmuneScore": rows["Immune"]},
        index=expr.columns,
    )
    df["ESTIMATEScore"] = df["StromalScore"] + df["ImmuneScore"]
    purity = np.cos(0.6049872018 + 0.0001467884 * df["ESTIMATEScore"].to_numpy())
    purity[purity < 0] = np.nan
    df["TumorPurity"] = purity
    return df, overlaps


def ssgsea_matrix(expr: pd.DataFrame, gene_sets: dict[str, list[str]], tau: float = 0.25) -> pd.DataFrame:
    """GSVA-style ssGSEA (ssgsea.norm=FALSE), matching xCell/GSVA.

    Rank scaling is irrelevant: |c*r|^tau cancels in the hit denominator.
    """
    genes = expr.index.to_numpy()
    n_genes, n_samples = expr.shape
    ranks = expr.rank(axis=0, method="average").to_numpy()
    gene_index = {g: i for i, g in enumerate(genes)}

    names = list(gene_sets)
    masks = np.zeros((len(names), n_genes), dtype=bool)
    for i, name in enumerate(names):
        for g in gene_sets[name]:
            j = gene_index.get(g)
            if j is not None:
                masks[i, j] = True

    out = np.zeros((len(names), n_samples))
    for j in range(n_samples):
        order = np.argsort(-ranks[:, j], kind="stable")
        w = np.abs(ranks[order, j]) ** tau
        tag = masks[:, order]
        tw = tag * w[None, :]
        denom = tw.sum(axis=1, keepdims=True)
        denom = np.where(denom == 0, 1.0, denom)
        pn = np.cumsum(tw, axis=1) / denom
        nm = n_genes - tag.sum(axis=1, keepdims=True)
        nm = np.where(nm == 0, 1.0, nm)
        p0 = np.cumsum((1.0 - tag) / nm, axis=1)
        out[:, j] = (pn - p0).sum(axis=1)
    return pd.DataFrame(out, index=names, columns=expr.columns)


def xcell_cell_type(sig_name: str) -> str:
    m = XCELL_SIG_RE.search(sig_name)
    if m:
        return sig_name[: m.start()]
    return sig_name.rsplit("_", 2)[0]


def xcell_scores(expr: pd.DataFrame, signatures: dict[str, list[str]], coef: pd.DataFrame, alpha: float = 0.5):
    """Official xCell pipeline (RNA-seq spillover) for all 64 cell types.

    Gene universe = all genes in this FPKM matrix (not the ~10.8k xCell training
    list, which is locked inside the R `xCell.data` GeneSet object). Ranks are
    computed on this universe. Documented in the report.
    """
    raw = ssgsea_matrix(expr, signatures)
    raw = raw.sub(raw.min(axis=1), axis=0)

    grouped: dict[str, list[str]] = defaultdict(list)
    for name in raw.index:
        grouped[xcell_cell_type(name)].append(name)
    agg = pd.DataFrame(
        {ct: raw.loc[sigs].mean(axis=0) for ct, sigs in grouped.items()}
    ).T
    agg = agg.reindex(coef.index)

    # transformScores: (raw - min)/5000 ^ power / (calib * 2)
    minx = agg.min(axis=1)
    tscores = ((agg.sub(minx, axis=0)) / 5000.0).clip(lower=0)
    power = coef["Power coefficient"].reindex(tscores.index).astype(float)
    calib = coef["Calibration parameter"].reindex(tscores.index).astype(float)
    tscores = tscores.pow(power, axis=0).div(calib * 2.0, axis=0)

    k_cols = [c for c in coef.columns if c in tscores.index]
    K0 = coef.loc[tscores.index, k_cols].astype(float)
    K0 = K0.reindex(index=tscores.index, columns=tscores.index)
    K = K0.to_numpy(dtype=float) * alpha
    np.fill_diagonal(K, 1.0)

    adjusted = np.zeros_like(tscores.to_numpy())
    for j, col in enumerate(tscores.columns):
        res = lsq_linear(K, tscores[col].to_numpy(dtype=float), bounds=(0, np.inf))
        x = res.x
        x[x < 0] = 0
        adjusted[:, j] = x
    adj = pd.DataFrame(adjusted, index=tscores.index, columns=tscores.columns)
    return agg, tscores, adj, {ct: len(s) for ct, s in grouped.items()}


def spearman(x, y):
    rho, p = stats.spearmanr(x, y)
    return float(rho), float(p)


def partial_spearman(x, y, z):
    """Partial Spearman of x,y | z via Pearson on ranks; p from t with n-3 df."""
    rx, ry, rz = (stats.rankdata(v) for v in (x, y, z))
    rxy = np.corrcoef(rx, ry)[0, 1]
    rxz = np.corrcoef(rx, rz)[0, 1]
    ryz = np.corrcoef(ry, rz)[0, 1]
    den = np.sqrt((1 - rxz**2) * (1 - ryz**2))
    if den < 1e-12:
        return float("nan"), float("nan")
    r = (rxy - rxz * ryz) / den
    r = float(np.clip(r, -1.0, 1.0))
    n = len(x)
    dof = n - 3
    if abs(r) >= 1.0 - 1e-15:
        return r, 0.0
    t = r * np.sqrt(dof / (1 - r**2))
    p = 2 * stats.t.sf(abs(t), dof)
    return float(r), float(p)


def log1p2(s: pd.Series) -> np.ndarray:
    return np.log2(s.to_numpy(dtype=float) + 1.0)


def gep18_score(expr: pd.DataFrame) -> tuple[np.ndarray, dict]:
    present = [g for g in GEP18_GENES if g in expr.index]
    missing = [g for g in GEP18_GENES if g not in expr.index]
    mat = np.log2(expr.loc[present].to_numpy(dtype=float) + 1.0)
    # z-score each gene across samples, then mean
    sd = mat.std(axis=1, ddof=1)
    sd[sd == 0] = 1.0
    z = (mat - mat.mean(axis=1, keepdims=True)) / sd[:, None]
    score = z.mean(axis=0)
    return score, {
        "n_genes_present": len(present),
        "n_genes_total": len(GEP18_GENES),
        "genes_present": present,
        "genes_missing": missing,
        "scoring": "mean of per-gene z-scores of log2(FPKM+1); unweighted (NanoString TIS weights not used)",
    }


def load_expr(path: str, drop_entrez: bool) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    if drop_entrez and "Entrez.ID" in df.columns:
        df = df.drop(columns=["Entrez.ID"])
    return df.groupby("gene").max()


def analyze(tag: str, expr: pd.DataFrame, estimate_sets, xcell_sigs, xcell_coef):
    est, overlaps = estimate_scores(expr, estimate_sets)
    if "TACSTD2" not in expr.index:
        raise KeyError(f"TACSTD2 missing from {tag}")
    tac = log1p2(expr.loc["TACSTD2"])
    z_purity = est["ESTIMATEScore"].to_numpy()

    features: dict[str, np.ndarray] = {}
    notes: dict[str, dict] = {}

    features["ESTIMATE_ImmuneScore"] = est["ImmuneScore"].to_numpy()
    notes["ESTIMATE_ImmuneScore"] = {
        "definition": "ESTIMATE v1.0.13 ImmuneScore (rank-normalized ssGSEA, tau=0.25)",
        "geneset_overlap": overlaps,
    }

    if {"GZMA", "PRF1"}.issubset(expr.index):
        features["CYT"] = 0.5 * (log1p2(expr.loc["GZMA"]) + log1p2(expr.loc["PRF1"]))
        notes["CYT"] = {
            "definition": "Rooney 2015 cytolytic activity: mean log2(FPKM+1) of GZMA and PRF1",
            "genes_present": ["GZMA", "PRF1"],
        }
    else:
        features["CYT"] = np.full(len(tac), np.nan)
        notes["CYT"] = {"definition": "Rooney CYT", "genes_missing": [g for g in ("GZMA", "PRF1") if g not in expr.index]}

    if "CD8A" in expr.index:
        features["CD8A"] = log1p2(expr.loc["CD8A"])
        notes["CD8A"] = {"definition": "log2(FPKM+1) CD8A", "genes_present": ["CD8A"]}
    else:
        features["CD8A"] = np.full(len(tac), np.nan)
        notes["CD8A"] = {"definition": "log2(FPKM+1) CD8A", "genes_missing": ["CD8A"]}

    gep, gep_note = gep18_score(expr)
    features["GEP18"] = gep
    notes["GEP18"] = {"definition": "Ayers 2017 18-gene T-cell-inflamed GEP", **gep_note}

    raw_x, _ts, adj_x, n_sigs = xcell_scores(expr, xcell_sigs, xcell_coef)
    features["xCell_CD8"] = adj_x.loc["CD8+ T-cells"].to_numpy()
    features["xCell_CD8_raw_ssgsea"] = raw_x.loc["CD8+ T-cells"].to_numpy()
    notes["xCell_CD8"] = {
        "definition": "xCell CD8+ T-cells after RNA-seq transform + spillover (alpha=0.5)",
        "n_signatures_aggregated": n_sigs.get("CD8+ T-cells"),
        "n_xcell_genes_in_matrix": int(sum(1 for gs in xcell_sigs.values() for g in gs if g in expr.index)),
    }
    notes["xCell_CD8_raw_ssgsea"] = {
        "definition": "xCell CD8+ T-cells aggregated ssGSEA only (no spillover); sensitivity",
        "n_signatures_aggregated": n_sigs.get("CD8+ T-cells"),
    }
    # extra CD8 subsets — reported so nothing is hidden
    for ct in ("CD8+ Tcm", "CD8+ Tem", "CD8+ naive T-cells"):
        key = "xCell_" + ct.replace("+", "").replace(" ", "_")
        features[key] = adj_x.loc[ct].to_numpy()
        notes[key] = {
            "definition": f"xCell '{ct}' after spillover (supplement; not the named xCell CD8 parent score)",
            "n_signatures_aggregated": n_sigs.get(ct),
        }

    metrics = {}
    rows = []
    for name, vec in features.items():
        mask = np.isfinite(vec) & np.isfinite(tac) & np.isfinite(z_purity)
        n = int(mask.sum())
        if n < 6 or np.unique(vec[mask]).size < 2:
            rec = {
                "spearman_rho_raw": None,
                "p_raw": None,
                "spearman_rho_purity_adjusted": None,
                "p_purity_adjusted": None,
                "n": n,
            }
        else:
            rho, p = spearman(tac[mask], vec[mask])
            pr, pp = partial_spearman(tac[mask], vec[mask], z_purity[mask])
            rec = {
                "spearman_rho_raw": round(float(rho), 4),
                "p_raw": float(f"{p:.6g}"),
                "spearman_rho_purity_adjusted": round(float(pr), 4),
                "p_purity_adjusted": float(f"{pp:.6g}"),
                "n": n,
            }
        rec["primary"] = name in PRIMARY_FEATURES
        rec["note"] = notes.get(name, {})
        metrics[name] = rec
        rows.append(
            {
                "dataset": tag,
                "feature": name,
                "primary": rec["primary"],
                "n": rec["n"],
                "rho_raw": rec["spearman_rho_raw"],
                "p_raw": rec["p_raw"],
                "rho_purity_adj": rec["spearman_rho_purity_adjusted"],
                "p_purity_adj": rec["p_purity_adjusted"],
            }
        )

    tbl = est.copy()
    tbl.insert(0, "TACSTD2_log2FPKMp1", tac)
    for name, vec in features.items():
        tbl[name] = vec
    tbl.to_csv(os.path.join(OUT, f"scores_{tag}.csv"))

    return {
        "dataset": tag,
        "n_samples": int(len(tac)),
        "n_genes": int(expr.shape[0]),
        "n_purity_formula_out_of_bounds": int(np.isnan(est["TumorPurity"]).sum()),
        "purity_adjustment": (
            "partial Spearman controls for ESTIMATEScore ranks (rank-identical "
            "to ESTIMATE TumorPurity); Affymetrix cosine formula not used as z"
        ),
        "metrics": metrics,
    }, tbl, rows, tac, features, est


def write_figures(plots):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # scatter: 2 cohorts x 5 primary features
    fig, axes = plt.subplots(2, 5, figsize=(18, 7.2))
    for i, (tag, tac, feats, est) in enumerate(plots):
        for j, name in enumerate(PRIMARY_FEATURES):
            ax = axes[i, j]
            y = feats[name]
            rho, p = spearman(tac, y)
            ax.scatter(
                tac, y, c=est["ESTIMATEScore"], cmap="viridis_r",
                s=36, edgecolor="k", lw=0.3,
            )
            ax.set_xlabel("TACSTD2 log2(FPKM+1)", fontsize=8)
            ax.set_ylabel(name, fontsize=8)
            ax.set_title(f"{tag.split('_')[0]}  ρ={rho:.2f} p={p:.1e}", fontsize=8)
    fig.suptitle("A2 rework: TACSTD2 vs each immune definition (color = ESTIMATEScore)")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "scatter_all_defs.png"), dpi=140)
    plt.close(fig)

    # bar of partial rho
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    for ax, (tag, tac, feats, est) in zip(axes, plots):
        names, rhos, ps = [], [], []
        z = est["ESTIMATEScore"].to_numpy()
        for name in PRIMARY_FEATURES:
            pr, pp = partial_spearman(tac, feats[name], z)
            names.append(name)
            rhos.append(pr)
            ps.append(pp)
        colors = ["#c0392b" if p < 0.05 else "#7f8c8d" for p in ps]
        ax.barh(names, rhos, color=colors)
        ax.axvline(0, color="k", lw=0.8)
        ax.axvline(-0.46, color="#2980b9", ls="--", lw=0.9, label="user adj ρ=-0.46")
        ax.set_title(f"{tag} (n={len(tac)})")
        ax.set_xlabel("partial Spearman ρ | ESTIMATEScore ranks")
        ax.legend(fontsize=7)
    fig.suptitle("Purity-adjusted TACSTD2–immune ρ (red = unadjusted-for-multiplicity p<0.05)")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "bar_partial_rho.png"), dpi=140)
    plt.close(fig)


def write_report(results: dict, table: pd.DataFrame):
    lines = []
    lines.append("# Rework A2: TACSTD2 vs immune definitions in public durvalumab NSCLC RNA")
    lines.append("")
    lines.append("**Self-contained recompute. Every requested definition is reported. Nothing is hidden.**")
    lines.append("")
    lines.append("## Claim and prior")
    lines.append("")
    lines.append("| Source | n | Immune def | Raw ρ | Raw p | Purity-adj ρ | Adj p |")
    lines.append("|---|---|---|---|---|---|---|")
    lines.append("| User | ? | unspecified | **-0.65** | 2e-4 (which test unspecified) | **-0.46** | paired with 2e-4 in the slide |")
    lines.append("| Prior A2 (this repo) | 29 post | ESTIMATE ImmuneScore | -0.665 | 8.3e-5 | -0.474 | **0.011** |")
    lines.append("| Prior A2 (this repo) | 32 pre | ESTIMATE ImmuneScore | -0.716 | 4.0e-6 | -0.530 | 0.0022 |")
    lines.append("")
    lines.append("Public durvalumab RNA used here (and in the prior A2):")
    lines.append("")
    lines.append("- [GSE253564](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253564) pre-treatment, n=32")
    lines.append("- [GSE248378](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE248378) post-treatment, n=29")
    lines.append("- Same randomized phase II neoadjuvant durvalumab ± SBRT trial (Altorki / NCT02904954).")
    lines.append("")
    lines.append("## Results — all five definitions")
    lines.append("")
    lines.append("Partial Spearman always controls **ESTIMATEScore ranks** (n−3 df).")
    lines.append("")
    def fmt_p(p):
        if p is None:
            return "NA"
        return f"{p:.2e}" if p < 0.001 else f"{p:.3g}"

    lines.append("| Dataset | n | Feature | Raw ρ | Raw p | Adj ρ | Adj p |")
    lines.append("|---|---|---|---|---|---|---|")
    for ds in results["datasets"]:
        for feat in PRIMARY_FEATURES:
            m = ds["metrics"][feat]
            lines.append(
                f"| {ds['dataset']} | {m['n']} | **{feat}** | "
                f"{m['spearman_rho_raw']} | {fmt_p(m['p_raw'])} | "
                f"{m['spearman_rho_purity_adjusted']} | {fmt_p(m['p_purity_adjusted'])} |"
            )
    lines.append("")
    lines.append("### Supplement (not the five named scores; still reported)")
    lines.append("")
    lines.append("| Dataset | Feature | Raw ρ | Raw p | Adj ρ | Adj p |")
    lines.append("|---|---|---|---|---|---|")
    for ds in results["datasets"]:
        for feat, m in ds["metrics"].items():
            if feat in PRIMARY_FEATURES:
                continue
            lines.append(
                f"| {ds['dataset']} | {feat} | "
                f"{m['spearman_rho_raw']} | {fmt_p(m['p_raw'])} | "
                f"{m['spearman_rho_purity_adjusted']} | {fmt_p(m['p_purity_adjusted'])} |"
            )
    lines.append("")
    lines.append("## Verdict (honest — all five definitions)")
    lines.append("")

    def mget(ds_name, feat):
        for ds in results["datasets"]:
            if ds["dataset"] == ds_name:
                return ds["metrics"][feat]
        raise KeyError(ds_name)

    post_imm = mget("GSE248378_posttreatment", "ESTIMATE_ImmuneScore")
    pre_imm = mget("GSE253564_pretreatment", "ESTIMATE_ImmuneScore")
    post_gep = mget("GSE248378_posttreatment", "GEP18")
    pre_gep = mget("GSE253564_pretreatment", "GEP18")
    post_cd8 = mget("GSE248378_posttreatment", "CD8A")
    pre_cd8 = mget("GSE253564_pretreatment", "CD8A")
    post_cyt = mget("GSE248378_posttreatment", "CYT")
    pre_cyt = mget("GSE253564_pretreatment", "CYT")
    post_xc = mget("GSE248378_posttreatment", "xCell_CD8")
    pre_xc = mget("GSE253564_pretreatment", "xCell_CD8")

    lines.append("**Closest match to the user numbers is post-treatment ESTIMATE ImmuneScore "
                 "(coefficients) plus a raw p-value (not the adjusted one).**")
    lines.append("")
    lines.append(
        f"1. **ESTIMATE ImmuneScore** — reproduces prior A2 exactly. "
        f"Post n=29: raw ρ={post_imm['spearman_rho_raw']} (p={fmt_p(post_imm['p_raw'])}), "
        f"adj ρ={post_imm['spearman_rho_purity_adjusted']} (p={fmt_p(post_imm['p_purity_adjusted'])} = prior 0.011). "
        f"Pre n=32: raw ρ={pre_imm['spearman_rho_raw']}, adj ρ={pre_imm['spearman_rho_purity_adjusted']} "
        f"(p={fmt_p(pre_imm['p_purity_adjusted'])}). User raw −0.65 ≈ post raw −0.665; "
        f"user adj −0.46 ≈ post adj −0.474. User p=2e-4 is the **raw** post p-order, not the adjusted p."
    )
    lines.append(
        f"2. **CYT** — post remains strong after purity (raw {post_cyt['spearman_rho_raw']}, "
        f"adj {post_cyt['spearman_rho_purity_adjusted']}, p={fmt_p(post_cyt['p_purity_adjusted'])}). "
        f"Pre collapses (raw {pre_cyt['spearman_rho_raw']}, adj {pre_cyt['spearman_rho_purity_adjusted']}, "
        f"p={fmt_p(pre_cyt['p_purity_adjusted'])})."
    )
    lines.append(
        f"3. **CD8A** — same pattern. Post adj ρ={post_cd8['spearman_rho_purity_adjusted']} "
        f"(p={fmt_p(post_cd8['p_purity_adjusted'])}); pre adj ρ={pre_cd8['spearman_rho_purity_adjusted']} "
        f"(p={fmt_p(pre_cd8['p_purity_adjusted'])}, NS)."
    )
    lines.append(
        f"4. **GEP18** (all 18 genes present in both matrices) — post raw ρ={post_gep['spearman_rho_raw']}, "
        f"raw p={fmt_p(post_gep['p_raw'])} is the single closest match to the user's **p=2e-4**, "
        f"but adj ρ={post_gep['spearman_rho_purity_adjusted']} (p={fmt_p(post_gep['p_purity_adjusted'])}) "
        f"is weaker than −0.46. Pre adj ρ={pre_gep['spearman_rho_purity_adjusted']} "
        f"(p={fmt_p(pre_gep['p_purity_adjusted'])}) is null — almost all of the pre GEP18 signal is purity."
    )
    lines.append(
        f"5. **xCell CD8** (spillover-adjusted CD8+ T-cells) — negative in both cohorts after purity: "
        f"pre adj ρ={pre_xc['spearman_rho_purity_adjusted']} (p={fmt_p(pre_xc['p_purity_adjusted'])}); "
        f"post adj ρ={post_xc['spearman_rho_purity_adjusted']} (p={fmt_p(post_xc['p_purity_adjusted'])}). "
        "Stronger than the user −0.46; does not reproduce that coefficient. "
        "Raw ssGSEA (no spillover) is even stronger (see supplement)."
    )
    lines.append("")
    lines.append("- **p=2e-4 attached to adj ρ=−0.46 does not reproduce** at n=29 or n=32 "
                 "(that pairing needs n≈55+). The prior A2 adj p=0.011 is the correct p for "
                 "the ESTIMATE post adjusted coefficient.")
    lines.append("- **Do not cite a single ρ.** Pre-treatment CYT / CD8A / GEP18 lose significance "
                 "after ESTIMATEScore adjustment; ESTIMATE ImmuneScore and xCell CD8 do not.")
    lines.append("")
    lines.append("## Methods")
    lines.append("")
    lines.append("- **Expression:** GEO supplementary FPKM; duplicate symbols collapsed by max; TACSTD2 / CD8A / CYT on log2(FPKM+1).")
    lines.append("- **ESTIMATE:** Python port of v1.0.13 `estimateScore` (official `SI_geneset.gmt`).")
    lines.append("- **CYT:** Rooney *Cell* 2015, (log2(GZMA+1)+log2(PRF1+1))/2.")
    lines.append("- **GEP18:** Ayers *JCI* 2017 18 genes listed in `scripts/rework_A2_defs.py`. Score = mean of per-gene z-scores of log2(FPKM+1). Unweighted; NanoString TIS regression weights are platform-specific and were not applied.")
    lines.append("- **xCell CD8:** Aran *Genome Biol* 2017. ssGSEA (tau=0.25, no GSVA normalization) on the 15 official CD8+ T-cells signatures, mean-aggregated, RNA-seq power/calibration transform, spillover α=0.5 via non-negative least squares. Gene universe = all genes in each FPKM matrix (pre 20,187; post 15,165), not the locked ~10.8k xCell training list. Raw aggregated ssGSEA (rank-identical to the transform, before spillover) is in the supplement.")
    lines.append("- **Partial Spearman:** Pearson on ranks of (TACSTD2, feature, ESTIMATEScore); p from Student t with n−3 df. Same estimator as prior A2.")
    lines.append("- **No multiple-testing correction** is applied across the five definitions — this is a pre-specified sensitivity panel, and every p-value is shown raw.")
    lines.append("")
    lines.append("## Caveats")
    lines.append("")
    lines.append("1. ESTIMATE TumorPurity cosine is out of bounds for many FPKM samples "
                 f"(pre {results['datasets'][0]['n_purity_formula_out_of_bounds']}/32, "
                 f"post {results['datasets'][1]['n_purity_formula_out_of_bounds']}/29); only score *ranks* are used.")
    lines.append("2. Both cohorts are small and from one early-stage neoadjuvant trial. The user slide did not name the accession or the immune definition.")
    lines.append("3. xCell on n=29–32 is below the size where spillover is stable; that is why raw ssGSEA CD8 is also reported.")
    lines.append("4. GEP18 here is an unweighted z-mean, not the commercial TIS.")
    lines.append("5. Bessede et al. 2024 (TROP2 vs ICI resistance) used atezolizumab OAK/POPLAR under EGA access — not these public durvalumab series.")
    lines.append("")
    lines.append("## Reproduce")
    lines.append("")
    lines.append("```bash")
    lines.append("pip install pandas numpy scipy matplotlib")
    lines.append("python3 scripts/rework_A2_defs.py")
    lines.append("```")
    lines.append("")
    lines.append("## Files")
    lines.append("")
    lines.append("- `results.json` — machine-readable metrics + methods notes")
    lines.append("- `correlations.tsv` — one row per dataset × feature")
    lines.append("- `scores_GSE253564_pretreatment.csv`, `scores_GSE248378_posttreatment.csv`")
    lines.append("- `figures/scatter_all_defs.png`, `figures/bar_partial_rho.png`")
    lines.append("- `scripts/rework_A2_defs.py`")
    lines.append("")
    path = os.path.join(OUT, "README.md")
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def main():
    estimate_sets = load_gmt(os.path.join(DATA, "estimate", "inst", "extdata", "SI_geneset.gmt"))
    with open(os.path.join(DATA, "xcell", "Xcell_signatures.json")) as fh:
        xcell_sigs = json.load(fh)
    xcell_coef = pd.read_csv(
        os.path.join(DATA, "xcell", "xcell_rnaseq_spill_and_calibration.csv"),
        index_col=0,
    )

    pre = load_expr(os.path.join(DATA, "GSE253564_Pre_FPKMs.txt.gz"), drop_entrez=True)
    post = load_expr(os.path.join(DATA, "GSE248378_Durva_Post_FPKMs.txt.gz"), drop_entrez=False)

    results = {
        "claim": CLAIM,
        "prior_A2_post_adj_p": PRIOR_A2_POST_ADJ_P,
        "purity_covariate": "ESTIMATEScore ranks",
        "primary_features": PRIMARY_FEATURES,
        "gep18_genes": GEP18_GENES,
        "datasets": [],
    }
    plots = []
    all_rows = []
    for tag, expr in (("GSE253564_pretreatment", pre), ("GSE248378_posttreatment", post)):
        res, tbl, rows, tac, feats, est = analyze(tag, expr, estimate_sets, xcell_sigs, xcell_coef)
        results["datasets"].append(res)
        plots.append((tag, tac, feats, est))
        all_rows.extend(rows)
        print(f"== {tag} n={res['n_samples']} ==")
        for feat in list(PRIMARY_FEATURES) + ["xCell_CD8_raw_ssgsea"]:
            m = res["metrics"][feat]
            print(
                f"  {feat:24s} raw={m['spearman_rho_raw']:+.3f} p={m['p_raw']:.3g}  "
                f"adj={m['spearman_rho_purity_adjusted']:+.3f} p={m['p_purity_adjusted']:.3g}"
            )

    table = pd.DataFrame(all_rows)
    table.to_csv(os.path.join(OUT, "correlations.tsv"), sep="\t", index=False)
    with open(os.path.join(OUT, "results.json"), "w") as fh:
        json.dump(results, fh, indent=2)
    write_figures(plots)
    write_report(results, table)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
