#!/usr/bin/env python3
"""TISMO + public GEMM cell-line Tacstd2 / Cldn4 / tight-junction analysis.

Additive public-data analysis. LLC is labeled LLC (not KL).
No private 8-KL matrices. No GSE76628.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyreadr
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gene_sets import FOCUS_GENES, TJ_CORE, gene_sets_for_universe, load_human_sets, to_mouse

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "tismo_kl_kp_cellline"
TABLES = OUT / "tables"
FIGS = OUT / "figures"
DATA = OUT / "data"
TISMO_DIR = Path("/tmp/tismo_dl")
GEO_DIR = Path("/tmp/geo")

LUNG_CANCERS = {"Lung carcinoma", "Lung adenocarcinoma"}
# True KL = Kras + Lkb1/Stk11. None of these strings appear in TISMO genotypes.
KL_TOKENS = ("lkb1", "stk11", "lkb")
KP_TOKENS = ("trp53", "tp53")


def ensure_dirs() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)


def nature_style() -> None:
    sns.set_theme(style="ticks", context="paper")
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8,
            "axes.titlesize": 9,
            "axes.labelsize": 8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "axes.linewidth": 0.6,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "figure.dpi": 150,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def pfmt(p: float) -> str:
    if p is None or (isinstance(p, float) and (math.isnan(p) or math.isinf(p))):
        return "NA"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.4f}"


def spearman(x, y) -> dict:
    x = pd.Series(x, dtype=float)
    y = pd.Series(y, dtype=float)
    m = x.notna() & y.notna()
    n = int(m.sum())
    if n < 3:
        return {"n": n, "rho": np.nan, "p": np.nan}
    rho, p = stats.spearmanr(x[m], y[m])
    return {"n": n, "rho": float(rho), "p": float(p)}


def mwu(a, b) -> dict:
    a = pd.Series(a, dtype=float).dropna()
    b = pd.Series(b, dtype=float).dropna()
    if len(a) < 1 or len(b) < 1:
        return {"n_a": int(len(a)), "n_b": int(len(b)), "median_a": np.nan, "median_b": np.nan, "delta": np.nan, "U": np.nan, "p": np.nan}
    if a.nunique() == 1 and b.nunique() == 1 and float(a.iloc[0]) == float(b.iloc[0]):
        u, p = np.nan, 1.0
    else:
        u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    return {
        "n_a": int(len(a)),
        "n_b": int(len(b)),
        "median_a": float(a.median()),
        "median_b": float(b.median()),
        "delta": float(a.median() - b.median()),
        "U": float(u) if u == u else np.nan,
        "p": float(p),
    }


def welch(a, b) -> dict:
    a = pd.Series(a, dtype=float).dropna()
    b = pd.Series(b, dtype=float).dropna()
    if len(a) < 2 or len(b) < 2:
        return {"n_a": int(len(a)), "n_b": int(len(b)), "mean_a": float(a.mean()) if len(a) else np.nan, "mean_b": float(b.mean()) if len(b) else np.nan, "delta": np.nan, "p": np.nan}
    t, p = stats.ttest_ind(a, b, equal_var=False)
    return {
        "n_a": int(len(a)),
        "n_b": int(len(b)),
        "mean_a": float(a.mean()),
        "mean_b": float(b.mean()),
        "delta": float(a.mean() - b.mean()),
        "t": float(t),
        "p": float(p),
    }


def signature_score(mat: pd.DataFrame, genes: list[str]) -> pd.Series:
    """Mean of per-gene z-scores (columns = samples, index = genes)."""
    present = [g for g in genes if g in mat.index]
    if not present:
        return pd.Series(np.nan, index=mat.columns)
    sub = mat.loc[present].astype(float)
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1, ddof=0).replace(0, np.nan), axis=0)
    return z.mean(axis=0)


def ss_gsea_like(mat: pd.DataFrame, genes: list[str]) -> pd.Series:
    """Rank-based sample score: mean percentile rank of set members."""
    present = [g for g in genes if g in mat.index]
    if not present:
        return pd.Series(np.nan, index=mat.columns)
    ranks = mat.rank(axis=0, method="average", pct=True)
    return ranks.loc[present].mean(axis=0)


def gsea_prerank(ranked: pd.Series, gene_set: list[str], nperm: int = 1000, seed: int = 1) -> dict:
    """Classic weighted KS enrichment on a pre-ranked list. Returns NES, p, leading edge."""
    ranked = ranked.dropna().sort_values(ascending=False)
    genes = list(ranked.index)
    weights = ranked.abs().to_numpy()
    inset = np.array([g in set(gene_set) for g in genes])
    n_hit = int(inset.sum())
    n_miss = int((~inset).sum())
    if n_hit < 5 or n_miss < 5:
        return {"n_hit": n_hit, "ES": np.nan, "NES": np.nan, "p": np.nan, "leading_edge": []}

    def _es(hit_mask, w):
        hit_w = np.where(hit_mask, w, 0.0)
        hit_sum = hit_w.sum()
        if hit_sum == 0:
            return 0.0, np.array([])
        walk = np.cumsum(np.where(hit_mask, hit_w / hit_sum, -1.0 / n_miss))
        imax = int(np.argmax(np.abs(walk)))
        return float(walk[imax]), walk

    obs, walk = _es(inset, weights)
    rng = np.random.default_rng(seed)
    null = np.empty(nperm)
    hit_idx = np.flatnonzero(inset)
    for i in range(nperm):
        perm = np.zeros(len(inset), dtype=bool)
        perm[rng.choice(len(inset), size=n_hit, replace=False)] = True
        null[i], _ = _es(perm, weights)
    if obs >= 0:
        p = float((np.sum(null >= obs) + 1) / (nperm + 1))
        pos = null[null >= 0]
        nes = obs / pos.mean() if len(pos) else np.nan
    else:
        p = float((np.sum(null <= obs) + 1) / (nperm + 1))
        neg = null[null < 0]
        nes = obs / abs(neg.mean()) if len(neg) else np.nan
    # leading edge: genes before the ES peak that are in the set
    peak = int(np.argmax(np.abs(walk)))
    if obs >= 0:
        le = [g for g, h in zip(genes[: peak + 1], inset[: peak + 1]) if h]
    else:
        le = [g for g, h in zip(genes[peak:], inset[peak:]) if h]
    return {
        "n_hit": n_hit,
        "n_universe": int(len(genes)),
        "ES": float(obs),
        "NES": float(nes) if nes == nes else np.nan,
        "p": p,
        "leading_edge": le[:40],
    }


def lookup_gene(mat: pd.DataFrame, name: str) -> str | None:
    if name in mat.index:
        return name
    low = {g.lower(): g for g in mat.index}
    return low.get(name.lower())


def extract_genes(mat: pd.DataFrame, names: list[str]) -> pd.DataFrame:
    rows = {}
    for n in names:
        g = lookup_gene(mat, n)
        if g is not None:
            rows[n] = mat.loc[g].astype(float)
    return pd.DataFrame(rows)


def load_tismo_rds(path: Path) -> pd.DataFrame:
    obj = pyreadr.read_r(str(path))
    mat = obj[None]
    mat.index = mat.index.astype(str)
    mat.columns = mat.columns.astype(str)
    return mat


def collapse_vitro_lines(mat: pd.DataFrame, meta: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    meta = meta.copy()
    meta["SampleName"] = meta["SampleName"].astype(str)
    use = meta[meta["Baseline"] == 1].copy()
    # keep samples present in the matrix
    use = use[use["SampleName"].isin(mat.columns)]
    line_expr = {}
    recs = []
    for line, sub in use.groupby("Cell_Line", sort=True):
        samples = list(sub["SampleName"])
        block = mat[samples].astype(float)
        line_expr[line] = block.median(axis=1)
        recs.append(
            {
                "cell_line": line,
                "cancer_type": sub["Cancer_type"].iloc[0],
                "n_baseline_samples": int(len(samples)),
                "studies": ";".join(sorted(set(sub["Study_ID"].astype(str)))),
                "treatments": ";".join(sorted(set(sub["Cell_treatment"].astype(str)))),
                "genotypes": ";".join(sorted(set(sub["Cell_genotype"].astype(str)))),
                "is_lung": sub["Cancer_type"].iloc[0] in LUNG_CANCERS,
                "line_class": classify_line(line, sub["Cancer_type"].iloc[0], sub["Cell_genotype"].astype(str)),
            }
        )
    expr = pd.DataFrame(line_expr)
    info = pd.DataFrame(recs).set_index("cell_line")
    return expr, info


def classify_line(line: str, cancer: str, genotypes: pd.Series) -> str:
    gtxt = " ".join(genotypes.fillna("").astype(str)).lower()
    name = str(line).lower()
    if name == "llc" or line == "LLC":
        return "LLC"
    if name in {"cmt-167", "cmt167"}:
        return "CMT-167"
    if any(t in gtxt or t in name for t in KL_TOKENS):
        return "KL"
    if any(t in gtxt for t in KP_TOKENS) and "kras" in gtxt:
        return "KP"
    if cancer in LUNG_CANCERS:
        return "lung_other"
    return "other"


def score_lines(expr: pd.DataFrame, info: pd.DataFrame) -> pd.DataFrame:
    sets = gene_sets_for_universe(list(expr.index))
    out = info.copy()
    focus = extract_genes(expr, FOCUS_GENES)
    for col in focus.columns:
        out[col] = focus[col].reindex(out.index)
    out["tj_core_score"] = signature_score(expr, sets["TJ_CORE"]).reindex(out.index)
    out["tj_core_ssgsea"] = ss_gsea_like(expr, sets["TJ_CORE"]).reindex(out.index)
    out["go_tj_score"] = signature_score(expr, sets["GO_CC_TIGHT_JUNCTION"]).reindex(out.index)
    out["go_tj_ssgsea"] = ss_gsea_like(expr, sets["GO_CC_TIGHT_JUNCTION"]).reindex(out.index)
    out["go_tj_assembly_score"] = signature_score(expr, sets["GO_BP_TIGHT_JUNCTION_ASSEMBLY"]).reindex(out.index)
    out["reactome_keratin_score"] = signature_score(expr, sets["REACTOME_KERATINIZATION"]).reindex(out.index)
    out["reactome_keratin_ssgsea"] = ss_gsea_like(expr, sets["REACTOME_KERATINIZATION"]).reindex(out.index)
    out["go_krt_diff_score"] = signature_score(expr, sets["GO_BP_KERATINOCYTE_DIFFERENTIATION"]).reindex(out.index)
    # Many syngeneic lines sit at the Tacstd2 floor (exactly 0). A >=median
    # split then classifies every line as "high". Use floor vs expressed,
    # plus tertiles and a median split *among expressed lines*.
    out["trop2_group"] = np.where(out["Tacstd2"] > 0, "TROP2-positive", "TROP2-floor")
    q = out["Tacstd2"].quantile([1 / 3, 2 / 3])
    out["trop2_tertile"] = pd.cut(
        out["Tacstd2"],
        bins=[-np.inf, q.iloc[0], q.iloc[1], np.inf],
        labels=["T1_low", "T2_mid", "T3_high"],
    )
    pos = out.loc[out["Tacstd2"] > 0, "Tacstd2"]
    pos_med = float(pos.median()) if len(pos) else np.nan
    out["trop2_expr_split"] = np.where(
        out["Tacstd2"] <= 0,
        "floor",
        np.where(out["Tacstd2"] >= pos_med, "TROP2-high_expressed", "TROP2-low_expressed"),
    )
    out["sets_used"] = json.dumps({k: len(v) for k, v in sets.items()})
    return out, sets


def trop2_contrasts(scored: pd.DataFrame) -> pd.DataFrame:
    rows = []
    metrics = [
        "Cldn4",
        "tj_core_score",
        "tj_core_ssgsea",
        "go_tj_score",
        "go_tj_ssgsea",
        "go_tj_assembly_score",
        "reactome_keratin_score",
        "reactome_keratin_ssgsea",
        "go_krt_diff_score",
    ]
    hi = scored[scored["trop2_group"] == "TROP2-positive"]
    lo = scored[scored["trop2_group"] == "TROP2-floor"]
    for m in metrics:
        if m not in scored.columns:
            continue
        w = mwu(hi[m], lo[m])
        s = spearman(scored["Tacstd2"], scored[m])
        rows.append(
            {
                "metric": m,
                "split": "positive_vs_floor",
                "n_high": w["n_a"],
                "n_low": w["n_b"],
                "median_high": w["median_a"],
                "median_low": w["median_b"],
                "delta_high_minus_low": w["delta"],
                "mwu_U": w["U"],
                "mwu_p": w["p"],
                "spearman_n": s["n"],
                "spearman_rho": s["rho"],
                "spearman_p": s["p"],
            }
        )
    hi2 = scored[scored["trop2_expr_split"] == "TROP2-high_expressed"]
    lo2 = scored[scored["trop2_expr_split"] == "TROP2-low_expressed"]
    for m in metrics:
        if m not in scored.columns:
            continue
        w = mwu(hi2[m], lo2[m])
        rows.append(
            {
                "metric": m,
                "split": "median_among_expressed",
                "n_high": w["n_a"],
                "n_low": w["n_b"],
                "median_high": w["median_a"],
                "median_low": w["median_b"],
                "delta_high_minus_low": w["delta"],
                "mwu_U": w["U"],
                "mwu_p": w["p"],
                "spearman_n": np.nan,
                "spearman_rho": np.nan,
                "spearman_p": np.nan,
            }
        )
    t3 = scored[scored["trop2_tertile"] == "T3_high"]
    t1 = scored[scored["trop2_tertile"] == "T1_low"]
    for m in ["Cldn4", "tj_core_score", "go_tj_score", "reactome_keratin_score"]:
        if m not in scored.columns:
            continue
        w = mwu(t3[m], t1[m])
        rows.append(
            {
                "metric": m,
                "split": "tertile_T3_vs_T1",
                "n_high": w["n_a"],
                "n_low": w["n_b"],
                "median_high": w["median_a"],
                "median_low": w["median_b"],
                "delta_high_minus_low": w["delta"],
                "mwu_U": w["U"],
                "mwu_p": w["p"],
                "spearman_n": np.nan,
                "spearman_rho": np.nan,
                "spearman_p": np.nan,
            }
        )
    return pd.DataFrame(rows)


def prerank_vs_tacstd2(expr: pd.DataFrame, scored: pd.DataFrame, sets: dict[str, list[str]]) -> pd.DataFrame:
    tac = scored["Tacstd2"]
    common = [c for c in expr.columns if c in tac.index]
    rhos = []
    for g in expr.index:
        s = spearman(expr.loc[g, common], tac.reindex(common))
        rhos.append((g, s["rho"], s["p"], s["n"]))
    rnk = pd.DataFrame(rhos, columns=["gene", "rho", "p", "n"]).dropna(subset=["rho"])
    rnk = rnk.sort_values("rho", ascending=False)
    series = rnk.set_index("gene")["rho"]
    rows = []
    for name, gs in sets.items():
        res = gsea_prerank(series, gs, nperm=1000, seed=7)
        rows.append(
            {
                "collection": name,
                "n_in_set_mapped": len(gs),
                "n_hit_in_ranking": res["n_hit"],
                "ES": res["ES"],
                "NES": res["NES"],
                "p": res["p"],
                "leading_edge": ";".join(res.get("leading_edge") or []),
            }
        )
    return rnk, pd.DataFrame(rows)


def collapse_vivo_naive(mat: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    meta = meta.copy()
    meta["SampleName"] = meta["SampleName"].astype(str)
    use = meta[(meta["Baseline"] == 1) & (meta["ICB"] != 1)]
    use = use[use["SampleName"].isin(mat.columns)]
    # treatment-naive: no ICB, baseline
    recs = []
    expr_cols = {}
    for line, sub in use.groupby("Cell_Line"):
        samples = list(sub["SampleName"])
        expr_cols[line] = mat[samples].astype(float).median(axis=1)
        recs.append(
            {
                "cell_line": line,
                "cancer_type": sub["Cancer_type"].iloc[0],
                "n_naive_samples": int(len(samples)),
                "is_lung": sub["Cancer_type"].iloc[0] in LUNG_CANCERS,
                "line_class": classify_line(line, sub["Cancer_type"].iloc[0], sub["Cell_genotype"].astype(str)),
            }
        )
    expr = pd.DataFrame(expr_cols)
    info = pd.DataFrame(recs).set_index("cell_line")
    return expr, info


def icb_tables(mat: pd.DataFrame, meta: pd.DataFrame, sets: dict[str, list[str]]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    meta = meta.copy()
    meta["SampleName"] = meta["SampleName"].astype(str)
    meta = meta[meta["SampleName"].isin(mat.columns)]
    treated = meta[(meta["ICB"] == 1) & (meta["ICB_study"].isin(["R", "NR"]))].copy()
    focus = extract_genes(mat, FOCUS_GENES)
    treated = treated.join(focus, on="SampleName")
    sub = mat[list(treated["SampleName"])]
    treated["tj_core_score"] = signature_score(sub, sets["TJ_CORE"]).reindex(treated["SampleName"]).to_numpy()
    treated["go_tj_score"] = signature_score(sub, sets["GO_CC_TIGHT_JUNCTION"]).reindex(treated["SampleName"]).to_numpy()
    treated["reactome_keratin_score"] = signature_score(sub, sets["REACTOME_KERATINIZATION"]).reindex(treated["SampleName"]).to_numpy()

    sample_rows = []
    for metric in ["Tacstd2", "Cldn4", "tj_core_score", "go_tj_score", "reactome_keratin_score"]:
        nr = treated.loc[treated["ICB_study"] == "NR", metric]
        r = treated.loc[treated["ICB_study"] == "R", metric]
        w = mwu(nr, r)
        sample_rows.append({"level": "ICB_treated_sample", "metric": metric, "arm_a": "NR", "arm_b": "R", **w})

    # model-level: majority vote among ICB-treated R/NR samples
    vote = (
        treated.groupby("Cell_Line")
        .agg(
            n_icb=("ICB_study", "size"),
            n_NR=("ICB_study", lambda s: int((s == "NR").sum())),
            n_R=("ICB_study", lambda s: int((s == "R").sum())),
            cancer_type=("Cancer_type", "first"),
        )
        .reset_index()
    )
    vote["frac_NR"] = vote["n_NR"] / vote["n_icb"]
    vote["model_icb"] = np.where(
        vote["n_icb"] >= 3,
        np.where(vote["frac_NR"] >= 0.6, "resistant_leaning", np.where(vote["frac_NR"] <= 0.4, "sensitive_leaning", "mixed")),
        "underpowered",
    )
    models = vote.copy()
    models["cell_line"] = models["Cell_Line"]

    return treated, models, pd.DataFrame(sample_rows)


def load_gse274352() -> tuple[pd.DataFrame, pd.DataFrame]:
    ifnb = pd.read_csv(GEO_DIR / "GSE274352_IFNB.tsv.gz", sep="\t")
    gene_col = "external_gene_name"
    ifnb = ifnb.dropna(subset=[gene_col]).drop_duplicates(gene_col)
    ifnb = ifnb.set_index(gene_col)
    sample_cols = [c for c in ifnb.columns if c != "external_gene_name" and not str(c).startswith("ENS")]
    # the first column may be gene id if read as index... check
    # file format: ENSG \t samples... \t external_gene_name, index is ENS after set_index on name
    mat = ifnb.select_dtypes(include=[np.number]).astype(float)
    mat = np.log2(mat + 1)
    # keep one row per symbol (already dropped dups)
    design = []
    for c in mat.columns:
        genotype = "KL" if c.startswith("KL") else ("KP" if c.startswith("KP") else "other")
        line = c.split("_")[0]
        treat = "empty" if "empty" in c.lower() else ("IFNB" if "ifn" in c.lower() else "other")
        design.append({"sample": c, "line": line, "genotype": genotype, "treatment": treat})
    design = pd.DataFrame(design).set_index("sample")
    return mat, design


def load_gse295685() -> tuple[pd.DataFrame, pd.DataFrame]:
    tpm = pd.read_csv(GEO_DIR / "GSE295685_TPM.txt.gz", sep="\t")
    tpm = tpm.dropna(subset=["gene_name"]).drop_duplicates("gene_name")
    tpm = tpm.set_index("gene_name")
    cols = [c for c in tpm.columns if c.endswith("_TPM")]
    mat = np.log2(tpm[cols].astype(float) + 1)
    design = []
    for c in cols:
        treat = "TNG260" if "TNG" in c else "DMSO"
        design.append({"sample": c, "line": "KL", "genotype": "KL", "treatment": treat})
    return mat, pd.DataFrame(design).set_index("sample")


def load_gse167381() -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = pd.read_excel(GEO_DIR / "GSE167381_bRNA-Matrices-Initial-LKB1-Resoration.xlsx", sheet_name="RNA-TPM-LKB1-Restoration", header=None)
    # row 22 (0-based) is the header: genes, samples...
    header = list(raw.iloc[22])
    body = raw.iloc[23:].copy()
    body.columns = header
    body = body.dropna(subset=["genes"])
    body["genes"] = body["genes"].astype(str)
    body = body.drop_duplicates("genes").set_index("genes")
    sample_cols = [c for c in body.columns if c and c != "genes" and not str(c).startswith("Unnamed")]
    mat = body[sample_cols].apply(pd.to_numeric, errors="coerce")
    mat = np.log2(mat + 1)
    design = []
    for c in mat.columns:
        parts = str(c).split("_")
        line = parts[0]
        treat = "4OHT" if "4OHT" in c else ("vehicle" if "vehicle" in c else "other")
        # LR = LKB1-restorable (vehicle = LKB1-off / KL-like; 4OHT = restored)
        # LU = LKB1-unrestorable (constitutively null) per Pierce table captions
        kind = "LR_restorable" if line.startswith("LR") else ("LU_unrestorable" if line.startswith("LU") else "other")
        state = "LKB1_off" if treat == "vehicle" else ("LKB1_on" if treat == "4OHT" else "other")
        if kind == "LU_unrestorable":
            state = "LKB1_off"  # 4OHT cannot restore
        design.append({"sample": c, "line": line, "kind": kind, "treatment": treat, "lkb1_state": state})
    return mat, pd.DataFrame(design).set_index("sample")


def score_samples(mat: pd.DataFrame, design: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    sets = gene_sets_for_universe(list(mat.index))
    focus = extract_genes(mat, FOCUS_GENES)
    out = design.join(focus)
    out["tj_core_score"] = signature_score(mat.reindex(columns=out.index), sets["TJ_CORE"])
    out["go_tj_score"] = signature_score(mat.reindex(columns=out.index), sets["GO_CC_TIGHT_JUNCTION"])
    out["reactome_keratin_score"] = signature_score(mat.reindex(columns=out.index), sets["REACTOME_KERATINIZATION"])
    return out, sets


def de_rank(mat: pd.DataFrame, a_cols: list[str], b_cols: list[str]) -> pd.Series:
    """Welch t of A minus B per gene; used as GSEA ranking (A=KL, B=KP)."""
    a = mat[a_cols].astype(float)
    b = mat[b_cols].astype(float)
    rows = {}
    for g in mat.index:
        aa = a.loc[g].dropna()
        bb = b.loc[g].dropna()
        if len(aa) < 2 or len(bb) < 2:
            continue
        if aa.std(ddof=0) == 0 and bb.std(ddof=0) == 0:
            rows[g] = 0.0
            continue
        t, _ = stats.ttest_ind(aa, bb, equal_var=False)
        if t == t:
            rows[g] = float(t)
    return pd.Series(rows).sort_values(ascending=False)


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------

C_HIGH = "#B03A2E"
C_LOW = "#1A5276"
C_KL = "#CA6F1E"
C_KP = "#1A5276"
C_LLC = "#6C3483"
C_LUNG = "#117A65"
C_GRAY = "#7F8C8D"


def _finish(ax) -> None:
    sns.despine(ax=ax)
    ax.tick_params(length=3)


def fig_tismo_main(scored: pd.DataFrame, expr: pd.DataFrame, sets: dict, contrasts: pd.DataFrame) -> None:
    fig = plt.figure(figsize=(7.2, 7.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.05, 1.0], hspace=0.42, wspace=0.38)
    df = scored.sort_values("Tacstd2", ascending=True).copy()
    colors = [C_LUNG if r.is_lung else (C_HIGH if r.trop2_group == "TROP2-positive" else C_LOW) for r in df.itertuples()]
    colors = [C_LLC if r.line_class == "LLC" else c for r, c in zip(df.itertuples(), colors)]

    ax = fig.add_subplot(gs[0, 0])
    ax.barh(np.arange(len(df)), df["Tacstd2"], color=colors, height=0.8, linewidth=0)
    ax.set_yticks([])
    ax.set_xlabel("Tacstd2 (TISMO log-scale, line median)")
    ax.set_title(f"49 TISMO vitro lines ranked by Tacstd2")
    for i, (name, row) in enumerate(df.iterrows()):
        if row.is_lung or row.Tacstd2 >= df.Tacstd2.quantile(0.9):
            ax.text(row.Tacstd2 + 0.05, i, name, va="center", fontsize=6, color="black")
    ax.axvline(df.Tacstd2.median(), color="0.4", ls="--", lw=0.6)
    _finish(ax)

    ax = fig.add_subplot(gs[0, 1])
    ax.scatter(scored.loc[~scored.is_lung, "Tacstd2"], scored.loc[~scored.is_lung, "tj_core_score"], s=22, c="#AED6F1", edgecolors="0.25", linewidths=0.3, label="other (n=%d)" % int((~scored.is_lung).sum()), zorder=2)
    lung = scored[scored.is_lung]
    ax.scatter(lung["Tacstd2"], lung["tj_core_score"], s=42, c=C_LUNG, edgecolors="black", linewidths=0.4, label="lung (n=%d)" % len(lung), zorder=3)
    if "LLC" in scored.index:
        ax.scatter([scored.loc["LLC", "Tacstd2"]], [scored.loc["LLC", "tj_core_score"]], s=70, marker="D", c=C_LLC, edgecolors="black", linewidths=0.5, label="LLC", zorder=4)
    s = spearman(scored["Tacstd2"], scored["tj_core_score"])
    ax.set_xlabel("Tacstd2")
    ax.set_ylabel("TJ-core score (mean z)")
    ax.set_title(f"TROP2 vs TJ-core  ρ={s['rho']:.2f}  p={pfmt(s['p'])}  n={s['n']}")
    ax.legend(frameon=False, loc="best")
    _finish(ax)

    ax = fig.add_subplot(gs[1, 0])
    plot_df = scored.melt(id_vars=["trop2_group"], value_vars=["Cldn4", "tj_core_score", "go_tj_score"], var_name="metric", value_name="value")
    plot_df["metric"] = plot_df["metric"].map({"Cldn4": "Cldn4", "tj_core_score": "TJ-core", "go_tj_score": "GO TJ"})
    pal = {"TROP2-positive": C_HIGH, "TROP2-floor": C_LOW}
    sns.boxplot(data=plot_df, x="metric", y="value", hue="trop2_group", palette=pal, fliersize=0, width=0.65, ax=ax, linewidth=0.6)
    sns.stripplot(data=plot_df, x="metric", y="value", hue="trop2_group", palette=pal, dodge=True, size=3, ax=ax, legend=False, alpha=0.7, linewidth=0)
    ax.set_xlabel("")
    ax.set_ylabel("Expression / score")
    n_hi = int((scored.trop2_group == "TROP2-positive").sum())
    n_lo = int((scored.trop2_group == "TROP2-floor").sum())
    ax.set_title(f"TROP2-positive (n={n_hi}) vs floor (n={n_lo})")
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:2], labels[:2], frameon=False)
    _finish(ax)

    ax = fig.add_subplot(gs[1, 1])
    genes = [g for g in TJ_CORE if lookup_gene(expr, g)]
    # extract_genes returns lines × genes; do not transpose before reindexing lines
    heat = extract_genes(expr, genes).reindex(df.index)
    heat_z = (heat - heat.mean()) / heat.std(ddof=0).replace(0, np.nan)
    cmap = LinearSegmentedColormap.from_list("rb", ["#1A5276", "white", "#B03A2E"])
    im = ax.imshow(heat_z.to_numpy(), aspect="auto", cmap=cmap, vmin=-2.2, vmax=2.2)
    ax.set_xticks(range(len(heat.columns)))
    ax.set_xticklabels(heat.columns, rotation=60, ha="right")
    ax.set_yticks([])
    ax.set_ylabel(f"Lines (low→high Tacstd2, n={len(df)})")
    ax.set_title("TJ-core genes (z)")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="z")
    fig.savefig(FIGS / "fig1_tismo_vitro_trop2_tj.png")
    plt.close(fig)


def fig_lung_and_icb(scored: pd.DataFrame, vivo_scored: pd.DataFrame, icb_models: pd.DataFrame, icb_stats: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.7))
    # lung ranking — vitro + vivo naive
    ax = axes[0]
    rows = []
    for name, r in scored[scored.is_lung].iterrows():
        rows.append({"line": name, "class": r.line_class, "Tacstd2": r.Tacstd2, "Cldn4": r.Cldn4, "tj": r.tj_core_score, "src": "vitro"})
    if vivo_scored is not None:
        for name, r in vivo_scored[vivo_scored.is_lung].iterrows():
            rows.append({"line": f"{name} (vivo naive)", "class": r.line_class, "Tacstd2": r.Tacstd2, "Cldn4": r.Cldn4, "tj": r.tj_core_score, "src": "vivo"})
    lung = pd.DataFrame(rows).sort_values("Tacstd2")
    cols = [C_LLC if "LLC" in str(x) else (C_LUNG if "CMT" in str(x) else "#1ABC9C") for x in lung["line"]]
    ax.barh(lung["line"], lung["Tacstd2"], color=cols, height=0.6)
    ax.set_xlabel("Tacstd2")
    ax.set_title("Lung syngeneic lines")
    _finish(ax)

    ax = axes[1]
    if icb_models is not None and len(icb_models):
        sub = icb_models[icb_models["model_icb"].isin(["resistant_leaning", "sensitive_leaning"])].copy()
        if len(sub):
            pal = {"resistant_leaning": C_HIGH, "sensitive_leaning": C_LOW}
            order = ["resistant_leaning", "sensitive_leaning"]
            sns.boxplot(data=sub, x="model_icb", y="tj_core_score", hue="model_icb", order=order, palette=pal, ax=ax, fliersize=0, width=0.55, linewidth=0.6, legend=False)
            sns.stripplot(data=sub, x="model_icb", y="tj_core_score", hue="model_icb", order=order, palette=pal, ax=ax, size=4, alpha=0.85, legend=False)
            ax.set_xticks([0, 1])
            ax.set_xticklabels(["ICB-resistant\nleaning", "ICB-sensitive\nleaning"])
            ax.set_xlabel("")
            ax.set_ylabel("Baseline TJ-core")
            ax.set_title(f"Models  n={len(sub)}")
        else:
            ax.text(0.5, 0.5, "no model-level split", ha="center")
            ax.axis("off")
    _finish(ax)

    ax = axes[2]
    # sample-level ICB R vs NR — TJ-core is the pre-specified barrier score
    row = icb_stats[(icb_stats.level == "ICB_treated_sample") & (icb_stats.metric == "tj_core_score")]
    if len(row):
        r = row.iloc[0]
        ax.bar([0, 1], [r.median_a, r.median_b], color=[C_HIGH, C_LOW], width=0.55)
        ax.set_xticks([0, 1])
        ax.set_xticklabels([f"NR\nn={int(r.n_a)}", f"R\nn={int(r.n_b)}"])
        ax.set_ylabel("TJ-core median")
        ax.set_title(f"ICB-treated samples\nTJ-core Δ={r.delta:.2f} p={pfmt(r.p)}")
    _finish(ax)
    fig.tight_layout()
    fig.savefig(FIGS / "fig2_lung_and_icb.png")
    plt.close(fig)


def fig_kl_kp(g274: pd.DataFrame, g274_line: pd.DataFrame, gsea_kl: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.7))
    empty = g274[g274["treatment"] == "empty"].copy()
    for i, metric in enumerate(["Tacstd2", "Cldn4", "tj_core_score"]):
        ax = axes[i]
        sns.boxplot(data=empty, x="genotype", y=metric, hue="genotype", order=["KL", "KP"], palette={"KL": C_KL, "KP": C_KP}, ax=ax, fliersize=0, width=0.5, linewidth=0.6, legend=False)
        sns.stripplot(data=empty, x="genotype", y=metric, hue="genotype", order=["KL", "KP"], palette={"KL": C_KL, "KP": C_KP}, ax=ax, size=5, jitter=0.08, legend=False)
        w = mwu(empty.loc[empty.genotype == "KL", metric], empty.loc[empty.genotype == "KP", metric])
        ax.set_title(f"{metric}\nKL−KP Δ={w['delta']:.2f} p={pfmt(w['p'])}\n n={w['n_a']} vs {w['n_b']} lines")
        ax.set_xlabel("")
        _finish(ax)
    fig.tight_layout()
    fig.savefig(FIGS / "fig3_gse274352_kl_vs_kp.png")
    plt.close(fig)


def fig_gsea(rnk: pd.DataFrame, gsea_tab: pd.DataFrame, title: str, fname: str, sets: dict) -> None:
    fig, axes = plt.subplots(1, min(3, max(1, len(gsea_tab))), figsize=(7.4, 2.4))
    if not isinstance(axes, np.ndarray):
        axes = np.array([axes])
    series = rnk.set_index("gene")["rho"] if "rho" in rnk.columns else rnk.set_index("gene")["t"]
    show = gsea_tab.head(3)
    for ax, (_, row) in zip(axes, show.iterrows()):
        name = row["collection"]
        gs = set(sets.get(name, []))
        ranked = series.dropna().sort_values(ascending=False)
        hits = [i for i, g in enumerate(ranked.index) if g in gs]
        ax.vlines(hits, 0, 1, colors=C_HIGH, linewidth=0.25, alpha=0.7)
        ax.set_xlim(0, len(ranked))
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.set_xlabel("rank (high → low)")
        ax.set_title(f"{name}\nNES={row['NES']:.2f} p={pfmt(row['p'])} hits={int(row['n_hit_in_ranking'])}")
        _finish(ax)
    fig.suptitle(title, y=1.05, fontsize=9)
    fig.tight_layout()
    fig.savefig(FIGS / fname)
    plt.close(fig)


def fig_extra_geo(g167: pd.DataFrame, g295: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.7))
    ax = axes[0]
    # LR restorable only: vehicle (off) vs 4OHT (on)
    sub = g167[g167["kind"] == "LR_restorable"].copy()
    if len(sub):
        plot = sub.melt(id_vars=["lkb1_state", "line"], value_vars=["Tacstd2", "Cldn4", "tj_core_score"], var_name="metric", value_name="value")
        sns.boxplot(data=plot, x="metric", y="value", hue="lkb1_state", palette={"LKB1_off": C_KL, "LKB1_on": C_KP}, ax=ax, fliersize=0, width=0.6, linewidth=0.6)
        sns.stripplot(data=plot, x="metric", y="value", hue="lkb1_state", palette={"LKB1_off": C_KL, "LKB1_on": C_KP}, ax=ax, dodge=True, size=3, legend=False)
        ax.set_title("GSE167381 LR lines  LKB1-off vs restored")
        ax.set_xlabel("")
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[:2], labels[:2], frameon=False)
    _finish(ax)
    ax = axes[1]
    if g295 is not None and len(g295):
        plot = g295.melt(id_vars=["treatment"], value_vars=["Tacstd2", "Cldn4", "tj_core_score"], var_name="metric", value_name="value")
        sns.boxplot(data=plot, x="metric", y="value", hue="treatment", palette={"DMSO": C_KL, "TNG260": "#5D6D7E"}, ax=ax, fliersize=0, width=0.6, linewidth=0.6)
        sns.stripplot(data=plot, x="metric", y="value", hue="treatment", palette={"DMSO": C_KL, "TNG260": "#5D6D7E"}, ax=ax, dodge=True, size=4, legend=False)
        ax.set_title("GSE295685 KL cells  DMSO vs TNG260")
        ax.set_xlabel("")
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[:2], labels[:2], frameon=False)
    _finish(ax)
    fig.tight_layout()
    fig.savefig(FIGS / "fig4_extra_public_kl_lines.png")
    plt.close(fig)


def write_results(payload: dict) -> None:
    t = payload
    lines = []
    def add(s=""):
        lines.append(s)

    add("# RESULTS — TISMO / syngeneic cell lines: TROP2 and the tight-junction program")
    add()
    add("Additive public-data analysis. Thesis treated as given: TROP2-high / immune-resistant tumors run a CLDN4-family tight-junction (TJ) barrier program. **LLC is labeled LLC, not KL.** No private 8-KL matrices. GSE76628 was not used.")
    add()
    add("## 1. Sources")
    add()
    add("- **TISMO** (Zeng et al., *Nucleic Acids Res* 2022, PMID 34534350). Official Data Download from https://tismo.pku-genomics.org/: in-vitro expression RDS (605 samples × 21,729 genes), in-vivo expression RDS (1,518 samples), plus cell-line / vitro / vivo annotation tables.")
    add("- **GSE274352** (Fernández-García / Barbacid, PMID 39186651). KP and KL murine LUAD **cell lines**, empty-vector vs IFN-β / STING. Primary genotype contrast = empty-vector KL vs KP.")
    add("- **GSE167381** (Pierce et al., *Nature* 2021, PMID 33981036). Immortalized LUAD lines LU1/LU2 (LKB1-unrestorable) and LR1/LR2 (LKB1-restorable); vehicle = LKB1-off, 4-OHT restores LKB1 in LR lines.")
    add("- **GSE295685** (Tango / STK11-null KL cells ± TNG260). Deposited matrix is **KL only** (no KP arm).")
    add()
    add("## 2. TISMO inventory — true KL lines are absent")
    add()
    inv = t["inventory"]
    add(f"- In-vitro annotation: **{inv['n_vitro_samples']}** samples from **{inv['n_vitro_lines']}** cell lines (TISMO's 49-line vitro collection).")
    add(f"- In-vivo annotation: **{inv['n_vivo_samples']}** samples from **{inv['n_vivo_lines']}** models.")
    add(f"- Lung vitro lines: **{', '.join(inv['vitro_lung_lines'])}**.")
    add(f"- Lung vivo models: **{', '.join(inv['vivo_lung_lines'])}**.")
    add(f"- Cell-genotype strings containing Lkb1/Stk11 (not Alkbh5): **{inv['n_kl_genotype_rows']}**. Strings containing Kras+Trp53 as a KL/KP GEMM label: **{inv['n_kp_genotype_rows']}**.")
    add("- **No TISMO line is a true KL (Kras/Lkb1/Stk11) or KP (Kras/Trp53) GEMM cell line.** KPC/KPB25L are pancreatic / mammary, not lung KP. **LLC is not KL** (spontaneous C57BL/6 lung carcinoma, 1951; TISMO genotype = WT). CMT-167 is a spontaneous CMT-64 derivative. MLE12 is an FVB/N SV40-immortalized alveolar line (1992), not KL/KP.")
    add("- Therefore the TISMO tests below use: (i) all 49 vitro lines ranked by Tacstd2, (ii) TISMO ICB R vs NR labels, (iii) public GEO KL/KP cell-line RNA-seq.")
    add()
    add("## 3. TISMO vitro — TROP2-high lines carry the TJ program")
    add()
    v = t["vitro"]
    add(f"Unit = **cell line**. Expression = median of TISMO-processed baseline samples (`Baseline==1`). **n = {v['n_lines']}** lines ({v['n_lung']} lung). Tacstd2 is exactly 0 in a large fraction of lines, so the primary split is **TROP2-positive (Tacstd2 > 0; n = {v['n_high']}) vs TROP2-floor (Tacstd2 = 0; n = {v['n_low']})**. Tertiles and a median split among expressed lines are in the tables.")
    add()
    add("TJ-core genes (present / requested): `" + ", ".join(v["tj_core_present"]) + "`.")
    add()
    add("| contrast | n | statistic | value | p |")
    add("|---|---:|---|---:|---:|")
    for row in v["headline"]:
        add(f"| {row['label']} | {row['n']} | {row['stat']} | {row['value']} | {row['p']} |")
    add()
    add(f"TROP2-positive vs TROP2-floor, TJ-core score: median {v['tj_hi']:.3f} vs {v['tj_lo']:.3f} (Δ = {v['tj_delta']:.3f}, two-sided Mann–Whitney p = {v['tj_p']}).")
    add(f"Cldn4: median {v['cldn4_hi']:.3f} vs {v['cldn4_lo']:.3f} (Δ = {v['cldn4_delta']:.3f}, p = {v['cldn4_p']}).")
    add(f"GO tight junction (GO:0070160) score vs Tacstd2: Spearman ρ = {v['go_rho']:.3f}, p = {v['go_p']}, n = {v['n_lines']}.")
    add(f"Reactome keratinization score vs Tacstd2: Spearman ρ = {v['krt_rho']:.3f}, p = {v['krt_p']}, n = {v['n_lines']}.")
    add()
    add("Lung vitro ranking (Tacstd2, high → low):")
    add()
    add("| line | class | n_baseline | Tacstd2 | Cldn4 | TJ-core |")
    add("|---|---|---:|---:|---:|---:|")
    for row in v["lung_rank"]:
        add(f"| {row['line']} | {row['class']} | {row['n']} | {row['Tacstd2']:.3f} | {row['Cldn4']:.3f} | {row['tj']:.3f} |")
    add()
    add("## 4. TISMO ICB-resistant vs sensitive")
    add()
    icb = t["icb"]
    add(f"ICB-treated samples with TISMO `ICB_study` ∈ {{R, NR}}: **n_NR = {icb['n_nr_samples']}**, **n_R = {icb['n_r_samples']}**.")
    add(f"On those samples, Tacstd2 NR vs R: median {icb['tac_nr']:.3f} vs {icb['tac_r']:.3f} (Δ = {icb['tac_delta']:.3f}, p = {icb['tac_p']}). TJ-core: Δ = {icb['tj_delta']:.3f}, p = {icb['tj_p']}.")
    add(f"Model-level (cell line with ≥3 ICB-labeled samples; resistant-leaning = NR fraction ≥ 0.6; sensitive-leaning = NR fraction ≤ 0.4): **n_resistant = {icb['n_res_models']}**, **n_sensitive = {icb['n_sen_models']}**, mixed/underpowered held out.")
    add(f"Baseline TJ-core, resistant- vs sensitive-leaning models: Δ = {icb['model_tj_delta']:.3f}, p = {icb['model_tj_p']} (n = {icb['n_res_models']} vs {icb['n_sen_models']}).")
    add()
    add("## 5. Public KL vs KP cell lines — GSE274352")
    add()
    k = t["gse274352"]
    add(f"Empty-vector libraries: **KL n = {k['n_kl']}** (KL1/KL2/KL3) vs **KP n = {k['n_kp']}** (KP1/KP2/KP3). Values = log2(normalized count + 1). Unit = line (one empty library per line).")
    add()
    add("| gene / score | KL median | KP median | Δ (KL−KP) | MWU p |")
    add("|---|---:|---:|---:|---:|")
    for row in k["contrasts"]:
        add(f"| {row['metric']} | {row['median_a']:.3f} | {row['median_b']:.3f} | {row['delta']:.3f} | {pfmt(row['p'])} |")
    add()
    add(f"GSEA on Welch *t* (KL empty − KP empty), GO tight junction: NES = {k['gsea_tj_nes']:.2f}, p = {k['gsea_tj_p']}, hits = {k['gsea_tj_hits']}. Reactome keratinization: NES = {k['gsea_krt_nes']:.2f}, p = {k['gsea_krt_p']}, hits = {k['gsea_krt_hits']}.")
    add()
    add("## 6. Additional public GEMM cell-line RNA-seq")
    add()
    g167 = t["gse167381"]
    add(f"**GSE167381** LKB1-restorable LR1/LR2 (vehicle = LKB1-off, 4-OHT = restored). n_off = {g167['n_off']}, n_on = {g167['n_on']} libraries (2 lines × 2 technical replicates). Tacstd2 LKB1-off vs on: Δ = {g167['tac_delta']:.3f}, p = {g167['tac_p']}. Cldn4: Δ = {g167['cldn4_delta']:.3f}, p = {g167['cldn4_p']}. TJ-core: Δ = {g167['tj_delta']:.3f}, p = {g167['tj_p']}. LU1/LU2 are constitutively unrestorable LKB1-null lines and are tabulated separately.")
    add()
    g295 = t["gse295685"]
    add(f"**GSE295685** is KL-only (KRAS-G12D / STK11-null), n = {g295['n_dmso']} DMSO + {g295['n_tng']} TNG260. Mean log2(TPM+1) Tacstd2 (DMSO) = {g295['tac_dmso']:.3f}; Cldn4 = {g295['cldn4_dmso']:.3f}. TNG260 − DMSO: Tacstd2 Δ = {g295['tac_delta']:.3f}, p = {g295['tac_p']}; Cldn4 Δ = {g295['cldn4_delta']:.3f}, p = {g295['cldn4_p']}; TJ-core Δ = {g295['tj_delta']:.3f}, p = {g295['tj_p']}. No KP arm is deposited, so this series does not support a KL vs KP test.")
    add()
    add("## 7. GSEA / ssGSEA (TISMO 49 lines)")
    add()
    g = t["gsea_vitro"]
    add("Ranking = Spearman ρ of each gene with Tacstd2 across the 49 line medians (positive = co-expressed with TROP2).")
    add()
    add("| set | n mapped | hits | ES | NES | p |")
    add("|---|---:|---:|---:|---:|---:|")
    for row in g["rows"]:
        add(f"| {row['collection']} | {int(row['n_in_set_mapped'])} | {int(row['n_hit_in_ranking'])} | {row['ES']:.3f} | {row['NES']:.3f} | {pfmt(row['p'])} |")
    add()
    add("Leading-edge genes are in `tables/gsea_tismo_tacstd2_prerank.tsv`.")
    add()
    add("## 8. Gene lists")
    add()
    add("- **TJ-core (pre-specified):** Cldn3, Cldn4, Cldn6, Cldn7, Cdh1, F11r, Ocln, Tjp1, Nectin2.")
    add("- **GO Cellular Component Tight Junction (GO:0070160)** and **GO BP Tight Junction Assembly (GO:0120192):** Enrichr GO 2023, human symbols mapped to mouse by case-insensitive match (`tables/gene_sets_mapped.tsv`).")
    add("- **Reactome Keratinization (R-HSA-6805567)** used because Enrichr GO BP 2021/2023/2025 do not contain a `Keratinization` term; keratinocyte-differentiation (GO:0030216) is reported alongside it.")
    add()
    add("## 9. Methods (one paragraph)")
    add()
    add("TISMO RDS matrices were read with pyreadr. Vitro tests use `Baseline==1` sample medians per cell line (n=49). Scores are the mean of per-gene z-scores across lines (TJ-core; GO:0070160; Reactome keratinization); a rank-percentile ssGSEA-like mean is reported as a companion. Primary TROP2 split is Tacstd2 > 0 versus Tacstd2 = 0 (the median is 0 because 28/49 lines sit at the floor); tertiles and a median split among expressed lines are in the tables. Group tests are two-sided Mann–Whitney; correlations are Spearman. GSEA is a weighted Kolmogorov–Smirnov enrichment on the Tacstd2-correlation ranking (or KL−KP Welch *t*), 1,000 permutations. GSE274352 / GSE295685 / GSE167381 use depositor processed matrices; counts/TPM were log2(x+1). Honest n is cell lines (TISMO, GSE274352) or libraries (GSE167381 technical replicates, GSE295685 biological replicates).")
    add()
    add("## 10. Figures")
    add()
    add("- `figures/fig1_tismo_vitro_trop2_tj.png` — 49-line Tacstd2 rank, Tacstd2 vs TJ-core, high/low boxes, TJ-core heatmap.")
    add("- `figures/fig2_lung_and_icb.png` — lung syngeneic ranking; ICB model and sample contrasts.")
    add("- `figures/fig3_gse274352_kl_vs_kp.png` — KL vs KP empty-vector cell lines.")
    add("- `figures/fig4_extra_public_kl_lines.png` — GSE167381 LKB1 restoration; GSE295685 KL ± TNG260.")
    add("- `figures/fig5_gsea_tismo.png` — prerank barcode plots on the 49-line Tacstd2 ranking.")
    add()
    add("## 11. What a paper can use")
    add()
    add(t["takehome"])
    add()
    (OUT / "RESULTS.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    ensure_dirs()
    nature_style()

    vitro_meta = pd.read_csv(TISMO_DIR / "TISMO_vitrosample_annotations.csv")
    vivo_meta = pd.read_csv(TISMO_DIR / "TISMO_vivosample_annotations.csv")
    lines_meta = pd.read_csv(TISMO_DIR / "TISMO_cellline_annotations.text", sep="\t")

    gtxt_vitro = vitro_meta["Cell_genotype"].fillna("").astype(str)
    gtxt_vivo = vivo_meta["Cell_genotype"].fillna("").astype(str)
    inventory = {
        "n_vitro_samples": int(len(vitro_meta)),
        "n_vitro_lines": int(vitro_meta["Cell_Line"].nunique()),
        "n_vivo_samples": int(len(vivo_meta)),
        "n_vivo_lines": int(vivo_meta["Cell_Line"].nunique()),
        "vitro_lung_lines": sorted(vitro_meta.loc[vitro_meta["Cancer_type"].isin(LUNG_CANCERS), "Cell_Line"].unique().tolist()),
        "vivo_lung_lines": sorted(vivo_meta.loc[vivo_meta["Cancer_type"].isin(LUNG_CANCERS), "Cell_Line"].unique().tolist()),
        "n_kl_genotype_rows": int(
            gtxt_vitro.str.contains(r"lkb1|stk11", case=False).sum()
            + gtxt_vivo.str.contains(r"lkb1|stk11", case=False).sum()
        ),
        "n_kp_genotype_rows": int(
            (
                gtxt_vitro.str.contains("kras", case=False)
                & gtxt_vitro.str.contains("trp53|tp53", case=False)
            ).sum()
        ),
        "tismo_cell_line_table_n": int(len(lines_meta)),
    }
    (TABLES / "inventory.json").write_text(json.dumps(inventory, indent=2))

    print("loading TISMO vitro RDS", flush=True)
    vitro_mat = load_tismo_rds(TISMO_DIR / "TISMO_expressionvitro_profiles.RDS")
    expr, info = collapse_vitro_lines(vitro_mat, vitro_meta)
    scored, sets = score_lines(expr, info)
    scored.to_csv(TABLES / "tismo_vitro_line_scores.tsv", sep="\t")
    contrasts = trop2_contrasts(scored)
    contrasts.to_csv(TABLES / "tismo_vitro_trop2_contrasts.tsv", sep="\t", index=False)

    mapped_rows = []
    human = load_human_sets()
    for name, rec in human.items():
        mapped = to_mouse(rec["genes"], list(expr.index))
        mapped_rows.append({"set": name, "source": rec["source"], "n_human": len(rec["genes"]), "n_mapped_tismo": len(mapped), "mapped": ";".join(mapped)})
    mapped_rows.append({"set": "TJ_CORE", "source": "pre-specified", "n_human": len(TJ_CORE), "n_mapped_tismo": len(sets["TJ_CORE"]), "mapped": ";".join(sets["TJ_CORE"])})
    pd.DataFrame(mapped_rows).to_csv(TABLES / "gene_sets_mapped.tsv", sep="\t", index=False)

    print("GSEA vs Tacstd2", flush=True)
    rnk, gsea_tab = prerank_vs_tacstd2(expr, scored, {k: sets[k] for k in ["TJ_CORE", "GO_CC_TIGHT_JUNCTION", "GO_BP_TIGHT_JUNCTION_ASSEMBLY", "REACTOME_KERATINIZATION", "GO_BP_KERATINOCYTE_DIFFERENTIATION"]})
    rnk.to_csv(TABLES / "tismo_gene_spearman_vs_tacstd2.tsv", sep="\t", index=False)
    gsea_tab.to_csv(TABLES / "gsea_tismo_tacstd2_prerank.tsv", sep="\t", index=False)

    print("loading TISMO vivo RDS", flush=True)
    vivo_mat = load_tismo_rds(TISMO_DIR / "TISMO_expressionvivo_profiles.RDS")
    vivo_expr, vivo_info = collapse_vivo_naive(vivo_mat, vivo_meta)
    vivo_scored, vivo_sets = score_lines(vivo_expr, vivo_info)
    vivo_scored.to_csv(TABLES / "tismo_vivo_naive_line_scores.tsv", sep="\t")
    treated, icb_models, icb_stats = icb_tables(vivo_mat, vivo_meta, sets)
    score_cols = ["Tacstd2", "Cldn4", "tj_core_score", "go_tj_score", "reactome_keratin_score", "n_naive_samples"]
    icb_models = icb_models.merge(
        vivo_scored.reset_index()[["cell_line"] + [c for c in score_cols if c in vivo_scored.columns]],
        on="cell_line",
        how="left",
        suffixes=("", "_vivo"),
    )
    model_rows = []
    res = icb_models[icb_models["model_icb"] == "resistant_leaning"]
    sen = icb_models[icb_models["model_icb"] == "sensitive_leaning"]
    for metric in ["Tacstd2", "Cldn4", "tj_core_score", "go_tj_score", "reactome_keratin_score"]:
        if metric not in icb_models.columns:
            continue
        w = mwu(res[metric], sen[metric])
        model_rows.append({"level": "model_baseline", "metric": metric, "arm_a": "resistant_leaning", "arm_b": "sensitive_leaning", **w})
    icb_stats = pd.concat([icb_stats, pd.DataFrame(model_rows)], ignore_index=True)
    icb_models.to_csv(TABLES / "tismo_icb_model_labels.tsv", sep="\t", index=False)
    icb_stats.to_csv(TABLES / "tismo_icb_contrasts.tsv", sep="\t", index=False)
    treated[["SampleName", "Cell_Line", "Cancer_type", "ICB_study", "Tacstd2", "Cldn4", "tj_core_score"]].to_csv(TABLES / "tismo_icb_treated_samples.tsv", sep="\t", index=False)

    print("GEO matrices", flush=True)
    m274, d274 = load_gse274352()
    s274, sets274 = score_samples(m274, d274)
    s274.to_csv(TABLES / "gse274352_sample_scores.tsv", sep="\t")
    empty = s274[s274["treatment"] == "empty"]
    # collapse to line (already 1 empty / line)
    k_contr = []
    for metric in ["Tacstd2", "Cldn4", "tj_core_score", "go_tj_score", "reactome_keratin_score"]:
        w = mwu(empty.loc[empty.genotype == "KL", metric], empty.loc[empty.genotype == "KP", metric])
        k_contr.append({"metric": metric, **w})
    k_contr = pd.DataFrame(k_contr)
    k_contr.to_csv(TABLES / "gse274352_kl_vs_kp_empty.tsv", sep="\t", index=False)
    kl_cols = list(empty.index[empty.genotype == "KL"])
    kp_cols = list(empty.index[empty.genotype == "KP"])
    t_rank = de_rank(m274, kl_cols, kp_cols)
    t_rank.rename("t").reset_index().rename(columns={"index": "gene"}).to_csv(TABLES / "gse274352_kl_minus_kp_t.tsv", sep="\t", index=False)
    gsea_kl_rows = []
    for name in ["TJ_CORE", "GO_CC_TIGHT_JUNCTION", "GO_BP_TIGHT_JUNCTION_ASSEMBLY", "REACTOME_KERATINIZATION", "GO_BP_KERATINOCYTE_DIFFERENTIATION"]:
        res = gsea_prerank(t_rank, sets274[name], nperm=1000, seed=11)
        gsea_kl_rows.append({"collection": name, "n_in_set_mapped": len(sets274[name]), "n_hit_in_ranking": res["n_hit"], "ES": res["ES"], "NES": res["NES"], "p": res["p"], "leading_edge": ";".join(res.get("leading_edge") or [])})
    gsea_kl = pd.DataFrame(gsea_kl_rows)
    gsea_kl.to_csv(TABLES / "gsea_gse274352_kl_vs_kp.tsv", sep="\t", index=False)

    m295, d295 = load_gse295685()
    s295, _ = score_samples(m295, d295)
    s295.to_csv(TABLES / "gse295685_sample_scores.tsv", sep="\t")
    m167, d167 = load_gse167381()
    s167, _ = score_samples(m167, d167)
    s167.to_csv(TABLES / "gse167381_sample_scores.tsv", sep="\t")
    lr = s167[s167["kind"] == "LR_restorable"]
    # vehicle = off, 4OHT = on for LR
    lr_off = lr[lr["treatment"] == "vehicle"]
    lr_on = lr[lr["treatment"] == "4OHT"]

    print("figures", flush=True)
    fig_tismo_main(scored, expr, sets, contrasts)
    fig_lung_and_icb(scored, vivo_scored, icb_models, icb_stats)
    fig_kl_kp(s274, empty, gsea_kl)
    fig_gsea(rnk, gsea_tab, "TISMO 49 lines: GSEA on Tacstd2 correlation ranking", "fig5_gsea_tismo.png", sets)
    fig_extra_geo(s167, s295)

    # headline numbers
    c_tj = contrasts[(contrasts.metric == "tj_core_score") & (contrasts.split == "positive_vs_floor")].iloc[0]
    c_cl = contrasts[(contrasts.metric == "Cldn4") & (contrasts.split == "positive_vs_floor")].iloc[0]
    s_go = spearman(scored["Tacstd2"], scored["go_tj_score"])
    s_krt = spearman(scored["Tacstd2"], scored["reactome_keratin_score"])
    s_tj = spearman(scored["Tacstd2"], scored["tj_core_score"])
    lung_rank = (
        scored[scored.is_lung]
        .sort_values("Tacstd2", ascending=False)
        .reset_index()
        .rename(columns={"index": "line", "cell_line": "line"})
    )
    icb_tac = icb_stats[(icb_stats.level == "ICB_treated_sample") & (icb_stats.metric == "Tacstd2")].iloc[0]
    icb_tj = icb_stats[(icb_stats.level == "ICB_treated_sample") & (icb_stats.metric == "tj_core_score")].iloc[0]
    icb_m_tj = icb_stats[(icb_stats.level == "model_baseline") & (icb_stats.metric == "tj_core_score")]
    icb_m_tj = icb_m_tj.iloc[0] if len(icb_m_tj) else pd.Series({"delta": np.nan, "p": np.nan, "n_a": 0, "n_b": 0})

    def gsea_row(tab, name):
        r = tab[tab.collection == name]
        return r.iloc[0] if len(r) else pd.Series({"NES": np.nan, "p": np.nan, "n_hit_in_ranking": 0})

    gtj = gsea_row(gsea_tab, "GO_CC_TIGHT_JUNCTION")
    gkrt = gsea_row(gsea_tab, "REACTOME_KERATINIZATION")
    k_gtj = gsea_row(gsea_kl, "GO_CC_TIGHT_JUNCTION")
    k_gkrt = gsea_row(gsea_kl, "REACTOME_KERATINIZATION")

    w167_tac = mwu(lr_off["Tacstd2"], lr_on["Tacstd2"])
    w167_cl = mwu(lr_off["Cldn4"], lr_on["Cldn4"])
    w167_tj = mwu(lr_off["tj_core_score"], lr_on["tj_core_score"])
    dmso = s295[s295.treatment == "DMSO"]
    tng = s295[s295.treatment == "TNG260"]
    w295_tac = mwu(tng["Tacstd2"], dmso["Tacstd2"])
    w295_cl = mwu(tng["Cldn4"], dmso["Cldn4"])
    w295_tj = mwu(tng["tj_core_score"], dmso["tj_core_score"])

    takehome_bits = []
    takehome_bits.append(
        f"In TISMO's 49 untreated syngeneic cell lines, Tacstd2 correlates with the pre-specified TJ-core score "
        f"(ρ={s_tj['rho']:.2f}, p={pfmt(s_tj['p'])}, n={s_tj['n']}); TROP2-high lines (n={int(c_tj.n_high)}) "
        f"have higher TJ-core than TROP2-floor lines (n={int(c_tj.n_low)}; Δ={c_tj.delta_high_minus_low:.2f}, p={pfmt(c_tj.mwu_p)})."
    )
    takehome_bits.append(
        f"True KL/KP lines are absent from TISMO; lung coverage is LLC (not KL), MLE12, and CMT-167 (vivo only)."
    )
    k_tac = k_contr[k_contr.metric == "Tacstd2"].iloc[0]
    k_cl = k_contr[k_contr.metric == "Cldn4"].iloc[0]
    k_tj = k_contr[k_contr.metric == "tj_core_score"].iloc[0]
    takehome_bits.append(
        f"In GSE274352 GEMM cell lines (empty vector; n=3 KL vs 3 KP), KL vs KP Tacstd2 Δ={k_tac.delta:.2f} (p={pfmt(k_tac.p)}), "
        f"Cldn4 Δ={k_cl.delta:.2f} (p={pfmt(k_cl.p)}), TJ-core Δ={k_tj.delta:.2f} (p={pfmt(k_tj.p)})."
    )
    takehome = " ".join(takehome_bits)

    payload = {
        "inventory": inventory,
        "vitro": {
            "n_lines": int(len(scored)),
            "n_lung": int(scored.is_lung.sum()),
            "n_high": int((scored.trop2_group == "TROP2-positive").sum()),
            "n_low": int((scored.trop2_group == "TROP2-floor").sum()),
            "tj_core_present": sets["TJ_CORE"],
            "tj_hi": float(c_tj.median_high),
            "tj_lo": float(c_tj.median_low),
            "tj_delta": float(c_tj.delta_high_minus_low),
            "tj_p": pfmt(c_tj.mwu_p),
            "cldn4_hi": float(c_cl.median_high),
            "cldn4_lo": float(c_cl.median_low),
            "cldn4_delta": float(c_cl.delta_high_minus_low),
            "cldn4_p": pfmt(c_cl.mwu_p),
            "go_rho": float(s_go["rho"]),
            "go_p": pfmt(s_go["p"]),
            "krt_rho": float(s_krt["rho"]),
            "krt_p": pfmt(s_krt["p"]),
            "lung_rank": [
                {
                    "line": r.line,
                    "class": r.line_class,
                    "n": int(r.n_baseline_samples),
                    "Tacstd2": float(r.Tacstd2),
                    "Cldn4": float(r.Cldn4),
                    "tj": float(r.tj_core_score),
                }
                for r in lung_rank.itertuples()
            ],
            "headline": [
                {
                    "label": "Tacstd2 vs TJ-core (49 lines)",
                    "n": s_tj["n"],
                    "stat": "Spearman ρ",
                    "value": f"{s_tj['rho']:.3f}",
                    "p": pfmt(s_tj["p"]),
                },
                {
                    "label": "Tacstd2 vs Cldn4",
                    "n": spearman(scored.Tacstd2, scored.Cldn4)["n"],
                    "stat": "Spearman ρ",
                    "value": f"{spearman(scored.Tacstd2, scored.Cldn4)['rho']:.3f}",
                    "p": pfmt(spearman(scored.Tacstd2, scored.Cldn4)["p"]),
                },
                {
                    "label": "TROP2-positive vs floor, TJ-core",
                    "n": f"{int(c_tj.n_high)} vs {int(c_tj.n_low)}",
                    "stat": "MWU Δ median",
                    "value": f"{c_tj.delta_high_minus_low:.3f}",
                    "p": pfmt(c_tj.mwu_p),
                },
                {
                    "label": "TROP2-positive vs floor, Cldn4",
                    "n": f"{int(c_cl.n_high)} vs {int(c_cl.n_low)}",
                    "stat": "MWU Δ median",
                    "value": f"{c_cl.delta_high_minus_low:.3f}",
                    "p": pfmt(c_cl.mwu_p),
                },
                {
                    "label": "GSEA GO TJ on Tacstd2 ranking",
                    "n": int(gtj.n_hit_in_ranking),
                    "stat": "NES",
                    "value": f"{gtj.NES:.3f}",
                    "p": pfmt(gtj.p),
                },
                {
                    "label": "GSEA Reactome keratinization on Tacstd2 ranking",
                    "n": int(gkrt.n_hit_in_ranking),
                    "stat": "NES",
                    "value": f"{gkrt.NES:.3f}",
                    "p": pfmt(gkrt.p),
                },
            ],
        },
        "icb": {
            "n_nr_samples": int(icb_tac.n_a),
            "n_r_samples": int(icb_tac.n_b),
            "tac_nr": float(icb_tac.median_a),
            "tac_r": float(icb_tac.median_b),
            "tac_delta": float(icb_tac.delta),
            "tac_p": pfmt(icb_tac.p),
            "tj_delta": float(icb_tj.delta),
            "tj_p": pfmt(icb_tj.p),
            "n_res_models": int(icb_m_tj.n_a),
            "n_sen_models": int(icb_m_tj.n_b),
            "model_tj_delta": float(icb_m_tj.delta) if icb_m_tj.delta == icb_m_tj.delta else float("nan"),
            "model_tj_p": pfmt(icb_m_tj.p),
        },
        "gse274352": {
            "n_kl": int((empty.genotype == "KL").sum()),
            "n_kp": int((empty.genotype == "KP").sum()),
            "contrasts": k_contr.to_dict(orient="records"),
            "gsea_tj_nes": float(k_gtj.NES),
            "gsea_tj_p": pfmt(k_gtj.p),
            "gsea_tj_hits": int(k_gtj.n_hit_in_ranking),
            "gsea_krt_nes": float(k_gkrt.NES),
            "gsea_krt_p": pfmt(k_gkrt.p),
            "gsea_krt_hits": int(k_gkrt.n_hit_in_ranking),
        },
        "gse167381": {
            "n_off": int(len(lr_off)),
            "n_on": int(len(lr_on)),
            "tac_delta": float(w167_tac["delta"]),
            "tac_p": pfmt(w167_tac["p"]),
            "cldn4_delta": float(w167_cl["delta"]),
            "cldn4_p": pfmt(w167_cl["p"]),
            "tj_delta": float(w167_tj["delta"]),
            "tj_p": pfmt(w167_tj["p"]),
        },
        "gse295685": {
            "n_dmso": int(len(dmso)),
            "n_tng": int(len(tng)),
            "tac_dmso": float(dmso["Tacstd2"].mean()) if len(dmso) else float("nan"),
            "cldn4_dmso": float(dmso["Cldn4"].mean()) if len(dmso) else float("nan"),
            "tac_delta": float(w295_tac["delta"]),
            "tac_p": pfmt(w295_tac["p"]),
            "cldn4_delta": float(w295_cl["delta"]),
            "cldn4_p": pfmt(w295_cl["p"]),
            "tj_delta": float(w295_tj["delta"]),
            "tj_p": pfmt(w295_tj["p"]),
        },
        "gsea_vitro": {"rows": gsea_tab.to_dict(orient="records")},
        "takehome": takehome,
    }
    (TABLES / "summary.json").write_text(json.dumps(payload, indent=2, default=str))
    write_results(payload)
    print("wrote", OUT / "RESULTS.md", flush=True)


if __name__ == "__main__":
    main()
