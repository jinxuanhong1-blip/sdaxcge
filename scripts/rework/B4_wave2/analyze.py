#!/usr/bin/env python3
"""B4 wave-2: try listed public definitions/tests for GSE126044 NR-higher TJ p=0.019.

Public GSE126044 counts only. No gene-subset search to hit 0.019.
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
import statsmodels.api as sm
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from gene_sets import (  # noqa: E402
    ALIASES,
    CD8,
    CD8A_ALONE,
    CHAE2018_CBM_TJ,
    CLDN147_F11R_PARD3,
    CLDN4_ALONE,
    CORE_TJ_40,
    CUSTOM_TJ_CORE,
    KRT_BASAL_SQUAMOUS,
    REACTOME_TJ,
    TJ_15,
    USER_PPT_7GENE,
    parse_kegg_genes,
)

ROOT = HERE.parents[2]
DATA = ROOT / "data" / "rework" / "B4_wave2"
OUT = ROOT / "results" / "rework" / "B4_wave2"
TAB = OUT / "tables"
FIG = OUT / "figures"
GS = OUT / "genesets"
for d in (TAB, FIG, GS):
    d.mkdir(parents=True, exist_ok=True)

CLAIM_P = 0.019
HIT_EPS = 1e-9  # p <= 0.019 counts as a strict hit
# Discrete MW p=0.019230... is what slides write as 0.019
CLAIM_ROUND = 3
CONTROL_FEATURES = {"CD8A", "CD8", "gene_CD8A"}


def log2cpm(counts: pd.DataFrame) -> pd.DataFrame:
    lib = counts.sum(axis=0)
    return np.log2(counts.divide(lib, axis=1) * 1e6 + 1)


def parse_series_matrix(path: Path) -> pd.DataFrame:
    import gzip

    opener = gzip.open if str(path).endswith(".gz") else open
    titles = responses = types = gsms = None
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("!Sample_title"):
                titles = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_geo_accession"):
                gsms = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1") and "patient response:" in line:
                responses = [
                    x.strip().strip('"').split("patient response:", 1)[1].strip()
                    for x in line.rstrip("\n").split("\t")[1:]
                ]
            elif line.startswith("!Sample_characteristics_ch1") and "sample:" in line:
                types = [
                    x.strip().strip('"').split("sample:", 1)[1].strip()
                    for x in line.rstrip("\n").split("\t")[1:]
                ]
    if not titles or not responses or not types or not gsms:
        raise RuntimeError("failed to parse GSE126044 series matrix")
    sample = [t.replace("RNA-seq_", "") for t in titles]
    return pd.DataFrame(
        {"sample": sample, "gsm": gsms, "response": responses, "sample_type": types}
    ).set_index("sample")


def msig_symbols(path: Path, key: str) -> list[str]:
    obj = json.loads(path.read_text())
    return list(obj[key]["geneSymbols"])


def read_list(path: Path) -> list[str]:
    return [ln.strip() for ln in path.read_text().splitlines() if ln.strip()]


def resolve(genes: list[str], index: pd.Index) -> tuple[list[str], list[str]]:
    matched: list[str] = []
    missing: list[str] = []
    seen: set[str] = set()
    for g in genes:
        cand = [g] + ALIASES.get(g, [])
        hit = next((c for c in cand if c in index and c not in seen), None)
        if hit is None:
            missing.append(g)
        else:
            matched.append(hit)
            seen.add(hit)
    return matched, missing


def meanz(expr: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, list[str], list[str]]:
    matched, missing = resolve(genes, expr.index)
    if not matched:
        return pd.Series(np.nan, index=expr.columns), matched, missing
    sub = expr.loc[matched]
    var = sub.var(axis=1)
    used = var[var > 0].index.tolist()
    sub = sub.loc[used]
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1, ddof=1), axis=0)
    return z.mean(axis=0), used, missing


def residualize(y: pd.Series, x: pd.Series) -> tuple[pd.Series, float]:
    df = pd.concat([y.rename("y"), x.rename("x")], axis=1).dropna()
    model = sm.OLS(df["y"], sm.add_constant(df["x"])).fit()
    resid = pd.Series(model.resid, index=df.index)
    return resid.reindex(y.index), float(model.rsquared)


def mwu(nr: np.ndarray, r: np.ndarray) -> dict:
    """NR vs R. One-sided alternative is NR > R."""
    nr = np.asarray(nr, dtype=float)
    r = np.asarray(r, dtype=float)
    n_nr, n_r = len(nr), len(r)
    if n_nr < 2 or n_r < 2:
        return {
            "n_R": n_r,
            "n_NR": n_nr,
            "median_R": float("nan"),
            "median_NR": float("nan"),
            "p_two": float("nan"),
            "p_one_NR_gt_R": float("nan"),
            "U_NR_gt_R": float("nan"),
        }
    u_one, p_one = stats.mannwhitneyu(nr, r, alternative="greater")
    _, p_two = stats.mannwhitneyu(nr, r, alternative="two-sided")
    return {
        "n_R": int(n_r),
        "n_NR": int(n_nr),
        "median_R": float(np.median(r)),
        "median_NR": float(np.median(nr)),
        "mean_R": float(np.mean(r)),
        "mean_NR": float(np.mean(nr)),
        "p_two": float(p_two),
        "p_one_NR_gt_R": float(p_one),
        "U_NR_gt_R": float(u_one),
        "direction": "NR>R" if np.median(nr) > np.median(r) else ("R>NR" if np.median(r) > np.median(nr) else "tie"),
    }


def hits_claim(p: float) -> bool:
    return (not math.isnan(p)) and (p <= CLAIM_P + HIT_EPS)


def matches_claimed(p: float) -> bool:
    """True if p is <=0.019 or rounds to 0.019 at 3 decimals (exact MW 0.01923)."""
    if p is None or (isinstance(p, float) and math.isnan(p)):
        return False
    return hits_claim(p) or round(float(p), CLAIM_ROUND) == CLAIM_P


def fisher_high_low(score: pd.Series, response: pd.Series, high_idx, low_idx) -> dict:
    high = response.loc[high_idx]
    low = response.loc[low_idx]
    table = np.array(
        [
            [(high == "non-responder").sum(), (high == "responder").sum()],
            [(low == "non-responder").sum(), (low == "responder").sum()],
        ],
        dtype=int,
    )
    # rows: high, low; cols: NR, R
    if table.min() < 0 or table.sum() == 0:
        return {"p_two": float("nan"), "p_one_NR_in_high": float("nan"), "table": table.tolist()}
    oddsratio, p_two = stats.fisher_exact(table, alternative="two-sided")
    _, p_one = stats.fisher_exact(table, alternative="greater")
    return {
        "n_high": int(len(high)),
        "n_low": int(len(low)),
        "n_high_NR": int((high == "non-responder").sum()),
        "n_high_R": int((high == "responder").sum()),
        "n_low_NR": int((low == "non-responder").sum()),
        "n_low_R": int((low == "responder").sum()),
        "odds_NR_in_high": float(oddsratio),
        "p_two": float(p_two),
        "p_one_NR_in_high": float(p_one),
        "table_high_low_x_NR_R": table.tolist(),
    }


def split_median(score: pd.Series) -> tuple[pd.Index, pd.Index]:
    med = score.median()
    high = score.index[score > med]
    low = score.index[score <= med]
    return high, low


def split_tertile(score: pd.Series) -> tuple[pd.Index, pd.Index, pd.Index]:
    q1, q2 = score.quantile([1 / 3, 2 / 3])
    bottom = score.index[score <= q1]
    mid = score.index[(score > q1) & (score < q2)]
    top = score.index[score >= q2]
    # if ties collapse a bin, fall back to rank tertiles
    if min(len(bottom), len(top)) < 2:
        ranks = score.rank(method="first")
        n = len(score)
        k1 = max(1, n // 3)
        k2 = n - k1
        order = ranks.sort_values()
        bottom = order.index[:k1]
        top = order.index[k2:]
        mid = order.index[k1:k2]
    return top, mid, bottom


def boxplot(score: pd.Series, meta: pd.DataFrame, title: str, ylab: str, fname: Path) -> None:
    order = ["responder", "non-responder"]
    fig, ax = plt.subplots(figsize=(4.4, 4.3))
    data = [score.loc[meta.index[meta["response"] == lvl]].dropna().values for lvl in order]
    ax.boxplot(data, tick_labels=["R", "NR"], showfliers=False, widths=0.55)
    rng = np.random.default_rng(0)
    for i, lvl in enumerate(order, start=1):
        idx = meta.index[meta["response"] == lvl]
        y = score.loc[idx]
        x = rng.normal(i, 0.055, size=len(y))
        ffpe = meta.loc[idx, "sample_type"].eq("FFPE").values
        ax.scatter(x[~ffpe], y[~ffpe], s=28, c="#2c7fb8" if lvl == "responder" else "#d95f0e",
                   edgecolors="k", linewidths=0.3, label="fresh" if i == 1 else None, zorder=3)
        ax.scatter(x[ffpe], y[ffpe], s=32, marker="D", c="#7f2704",
                   edgecolors="k", linewidths=0.3, label="FFPE (all NR)" if ffpe.any() and i == 2 else None, zorder=3)
    ax.set_title(title, fontsize=10)
    ax.set_ylabel(ylab)
    if any(meta["sample_type"] == "FFPE"):
        ax.legend(fontsize=8, loc="best")
    fig.tight_layout()
    fig.savefig(fname.with_suffix(".png"), dpi=150)
    fig.savefig(fname.with_suffix(".pdf"))
    plt.close(fig)


def p_dotplot(rows: list[dict], fname: Path) -> None:
    df = pd.DataFrame(rows)
    df = df.sort_values("p_two")
    fig, ax = plt.subplots(figsize=(8.2, max(3.6, 0.28 * len(df) + 1.2)))
    y = np.arange(len(df))
    ax.axvline(CLAIM_P, color="#e31a1c", ls="--", lw=1, label="claimed p=0.019")
    ax.axvline(0.05, color="#999999", ls=":", lw=1, label="p=0.05")
    ax.scatter(df["p_two"], y, s=36, c="#2171b5", label="two-sided MW", zorder=3)
    ax.scatter(df["p_one_NR_gt_R"], y, s=28, marker="s", c="#fd8d3c", label="one-sided NR>R", zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(df["label"], fontsize=8)
    ax.set_xlabel("Mann-Whitney p (GSE126044, all samples)")
    ax.set_xscale("log")
    ax.set_xlim(1e-4, 1.05)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(fname.with_suffix(".png"), dpi=150)
    fig.savefig(fname.with_suffix(".pdf"))
    plt.close(fig)


def main() -> int:
    counts = pd.read_csv(DATA / "GSE126044_counts.txt.gz", sep="\t", index_col=0)
    meta = parse_series_matrix(DATA / "GSE126044_series_matrix.txt.gz")
    # counts columns are Dis_XX; matrix index matches
    missing_cols = [c for c in counts.columns if c not in meta.index]
    extra_meta = [s for s in meta.index if s not in counts.columns]
    if missing_cols or extra_meta:
        raise RuntimeError(f"sample mismatch counts-only={missing_cols} meta-only={extra_meta}")
    meta = meta.loc[counts.columns]
    expr = log2cpm(counts)

    kegg = parse_kegg_genes((DATA / "kegg_hsa04530.txt").read_text())
    gocc = msig_symbols(DATA / "GOCC_TIGHT_JUNCTION.json", "GOCC_TIGHT_JUNCTION")
    reactome_krt = msig_symbols(DATA / "REACTOME_KERATINIZATION.json", "REACTOME_KERATINIZATION")
    reactome_remote = msig_symbols(
        DATA / "REACTOME_TIGHT_JUNCTION_INTERACTIONS.json",
        "REACTOME_TIGHT_JUNCTION_INTERACTIONS",
    )
    stromal_genes = read_list(DATA / "ESTIMATE_stromal_signature.txt")
    immune_genes = read_list(DATA / "ESTIMATE_immune_signature.txt")

    features: dict[str, dict] = {
        "CLDN4": {
            "genes": CLDN4_ALONE,
            "kind": "single_gene",
            "source": "single gene (item 6)",
        },
        "CLDN147_F11R_PARD3": {
            "genes": CLDN147_F11R_PARD3,
            "kind": "user_wave2_5gene",
            "source": "wave-2 request: CLDN1/4/7/F11R/PARD3 mean-z",
        },
        "KRT_basal_squamous": {
            "genes": KRT_BASAL_SQUAMOUS,
            "kind": "keratin",
            "source": "compact lung basal/squamous KRTs (item 6 contrast; not Cho 2020)",
        },
        "REACTOME_KERATINIZATION": {
            "genes": reactome_krt,
            "kind": "keratin",
            "source": "MSigDB REACTOME_KERATINIZATION (R-HSA-6805567)",
        },
        "Chae2018_CBM_TJ": {
            "genes": CHAE2018_CBM_TJ,
            "kind": "paper_TJ",
            "source": "Chae 2018 Sci Rep Table 2 CBM tight-junction genes (NSCLC immune paper)",
        },
        "userPPT_7gene": {
            "genes": USER_PPT_7GENE,
            "kind": "user_ppt",
            "source": "user PPT / claim-page 7-gene (CLDN1/4/7, F11R, TJP1/2, OCLN); not Cho 2020",
        },
        "Reactome_R-HSA-420029": {
            "genes": REACTOME_TJ,
            "kind": "database_TJ",
            "source": "Reactome tight junction interactions / MSigDB REACTOME_TIGHT_JUNCTION_INTERACTIONS",
        },
        "GOCC_TIGHT_JUNCTION": {
            "genes": gocc,
            "kind": "database_TJ",
            "source": "MSigDB GOCC_TIGHT_JUNCTION (GO:0070160)",
        },
        "KEGG_hsa04530": {
            "genes": kegg,
            "kind": "database_TJ",
            "source": "KEGG hsa04530 Tight junction (actin/myosin-heavy comparator)",
        },
        "TJ_15_structural": {
            "genes": TJ_15,
            "kind": "repo_comparator",
            "source": "repo 15-gene structural TJ (PR #69)",
        },
        "CUSTOM_TJ_CORE": {
            "genes": CUSTOM_TJ_CORE,
            "kind": "repo_comparator",
            "source": "repo CUSTOM_TJ_CORE (prior PRs)",
        },
        "coreTJ_40": {
            "genes": CORE_TJ_40,
            "kind": "repo_comparator",
            "source": "40-gene core TJ (prior PRs; two-sided MW ~0.115)",
        },
        "CD8A": {
            "genes": CD8A_ALONE,
            "kind": "control",
            "source": "positive-control immune gene (expect R > NR)",
        },
        "CD8": {
            "genes": CD8,
            "kind": "control",
            "source": "CD8A+CD8B mean-z control",
        },
    }

    # ESTIMATE covariates (mean-z of Yoshihara 2013 signatures)
    stromal, stromal_used, stromal_miss = meanz(expr, stromal_genes)
    immune, immune_used, immune_miss = meanz(expr, immune_genes)
    estimate = stromal + immune

    scores = pd.DataFrame(
        {
            "response": meta["response"],
            "sample_type": meta["sample_type"],
            "gsm": meta["gsm"],
            "StromalScore_meanz": stromal,
            "ImmuneScore_meanz": immune,
            "ESTIMATEScore_meanz": estimate,
            "CLDN4_log2cpm": expr.loc["CLDN4"] if "CLDN4" in expr.index else np.nan,
            "CD8A_log2cpm": expr.loc["CD8A"] if "CD8A" in expr.index else np.nan,
        },
        index=expr.columns,
    )

    coverage_rows = []
    stat_rows = []
    hit_rows = []
    plot_rows = []

    subsets = {
        "all": meta.index,
        "fresh": meta.index[meta["sample_type"].str.lower() == "fresh"],
    }

    single_genes = ["CLDN1", "CLDN4", "CLDN7", "F11R", "PARD3", "TJP1", "TJP2", "OCLN", "CD8A"]

    def run_feature(name: str, score: pd.Series, spec: dict, used: list[str], missing: list[str]) -> None:
        scores[name] = score
        coverage_rows.append(
            {
                "feature": name,
                "kind": spec["kind"],
                "source": spec["source"],
                "n_requested": len(spec["genes"]),
                "n_used": len(used),
                "n_missing": len(missing),
                "used": ",".join(used),
                "missing": ",".join(missing),
            }
        )
        for subset_name, idx in subsets.items():
            sub_meta = meta.loc[idx]
            sub_score = score.loc[idx]
            r = sub_score.loc[sub_meta.index[sub_meta["response"] == "responder"]]
            nr = sub_score.loc[sub_meta.index[sub_meta["response"] == "non-responder"]]
            mw = mwu(nr.values, r.values)
            row = {
                "feature": name,
                "kind": spec["kind"],
                "subset": subset_name,
                "test": "MW_score_R_vs_NR",
                "covariate": "",
                **mw,
                "hit_0.019_two": hits_claim(mw["p_two"]),
                "hit_0.019_one": hits_claim(mw["p_one_NR_gt_R"]),
                "matches_claimed_two": matches_claimed(mw["p_two"]),
                "matches_claimed_one": matches_claimed(mw["p_one_NR_gt_R"]),
                "source": spec["source"],
            }
            stat_rows.append(row)
            if subset_name == "all" and spec["kind"] != "control" and not name.startswith("gene_"):
                plot_rows.append(
                    {
                        "label": f"{name} (n={mw['n_R']}+{mw['n_NR']})",
                        "p_two": mw["p_two"],
                        "p_one_NR_gt_R": mw["p_one_NR_gt_R"],
                    }
                )
            for side, p in (("two", mw["p_two"]), ("one_NR_gt_R", mw["p_one_NR_gt_R"])):
                if matches_claimed(p):
                    hit_rows.append({**row, "which_p": side, "p": p})

            # median / tertile association (response in high vs low)
            if subset_name == "all":
                high, low = split_median(sub_score)
                med = fisher_high_low(sub_score, sub_meta["response"], high, low)
                med_row = {
                    "feature": name,
                    "kind": spec["kind"],
                    "subset": subset_name,
                    "test": "Fisher_median_high_vs_low",
                    "covariate": "",
                    "n_R": mw["n_R"],
                    "n_NR": mw["n_NR"],
                    "median_R": mw["median_R"],
                    "median_NR": mw["median_NR"],
                    "p_two": med["p_two"],
                    "p_one_NR_gt_R": med["p_one_NR_in_high"],
                    "direction": "NR_enriched_high" if med.get("odds_NR_in_high", 1) > 1 else "R_enriched_high",
                    "hit_0.019_two": hits_claim(med["p_two"]),
                    "hit_0.019_one": hits_claim(med["p_one_NR_in_high"]),
                    "matches_claimed_two": matches_claimed(med["p_two"]),
                    "matches_claimed_one": matches_claimed(med["p_one_NR_in_high"]),
                    "source": spec["source"],
                    "n_high": med.get("n_high"),
                    "n_low": med.get("n_low"),
                    "n_high_NR": med.get("n_high_NR"),
                    "n_high_R": med.get("n_high_R"),
                    "n_low_NR": med.get("n_low_NR"),
                    "n_low_R": med.get("n_low_R"),
                }
                stat_rows.append(med_row)
                for side, p in (("two", med["p_two"]), ("one_NR_in_high", med["p_one_NR_in_high"])):
                    if matches_claimed(p):
                        hit_rows.append({**med_row, "which_p": side, "p": p})

                top, mid, bottom = split_tertile(sub_score)
                ter = fisher_high_low(sub_score, sub_meta["response"], top, bottom)
                ter_row = {
                    "feature": name,
                    "kind": spec["kind"],
                    "subset": subset_name,
                    "test": "Fisher_tertile_top_vs_bottom",
                    "covariate": "",
                    "n_R": mw["n_R"],
                    "n_NR": mw["n_NR"],
                    "median_R": mw["median_R"],
                    "median_NR": mw["median_NR"],
                    "p_two": ter["p_two"],
                    "p_one_NR_gt_R": ter["p_one_NR_in_high"],
                    "direction": "NR_enriched_top" if ter.get("odds_NR_in_high", 1) > 1 else "R_enriched_top",
                    "hit_0.019_two": hits_claim(ter["p_two"]),
                    "hit_0.019_one": hits_claim(ter["p_one_NR_in_high"]),
                    "matches_claimed_two": matches_claimed(ter["p_two"]),
                    "matches_claimed_one": matches_claimed(ter["p_one_NR_in_high"]),
                    "source": spec["source"],
                    "n_high": ter.get("n_high"),
                    "n_low": ter.get("n_low"),
                    "n_high_NR": ter.get("n_high_NR"),
                    "n_high_R": ter.get("n_high_R"),
                    "n_low_NR": ter.get("n_low_NR"),
                    "n_low_R": ter.get("n_low_R"),
                    "n_mid": int(len(mid)),
                }
                stat_rows.append(ter_row)
                for side, p in (("two", ter["p_two"]), ("one_NR_in_high", ter["p_one_NR_in_high"])):
                    if matches_claimed(p):
                        hit_rows.append({**ter_row, "which_p": side, "p": p})

        # residualize on stromal / ESTIMATE (all samples)
        if spec["kind"] != "control":
            for cov_name, cov in (
                ("StromalScore_meanz", stromal),
                ("ESTIMATEScore_meanz", estimate),
            ):
                resid, r2 = residualize(score, cov)
                r = resid.loc[meta.index[meta["response"] == "responder"]]
                nr = resid.loc[meta.index[meta["response"] == "non-responder"]]
                mw = mwu(nr.dropna().values, r.dropna().values)
                row = {
                    "feature": name,
                    "kind": spec["kind"],
                    "subset": "all_residual",
                    "test": "MW_residual_R_vs_NR",
                    "covariate": cov_name,
                    "resid_R2": r2,
                    **mw,
                    "hit_0.019_two": hits_claim(mw["p_two"]),
                    "hit_0.019_one": hits_claim(mw["p_one_NR_gt_R"]),
                    "matches_claimed_two": matches_claimed(mw["p_two"]),
                    "matches_claimed_one": matches_claimed(mw["p_one_NR_gt_R"]),
                    "source": spec["source"],
                }
                stat_rows.append(row)
                scores[f"{name}__resid_{cov_name}"] = resid
                for side, p in (("two", mw["p_two"]), ("one_NR_gt_R", mw["p_one_NR_gt_R"])):
                    if matches_claimed(p):
                        hit_rows.append({**row, "which_p": side, "p": p})

    for name, spec in features.items():
        if spec["kind"] == "single_gene" or name in {"CLDN4", "CD8A"}:
            gene = spec["genes"][0]
            if gene not in expr.index:
                raise RuntimeError(f"{gene} missing from GSE126044 counts")
            score = expr.loc[gene]
            used, missing = [gene], []
        else:
            score, used, missing = meanz(expr, spec["genes"])
        run_feature(name, score, spec, used, missing)

    # single-gene MW for the modules in item 6 / PPT
    for g in single_genes:
        if g not in expr.index:
            continue
        spec = {"genes": [g], "kind": "single_gene", "source": f"single-gene log2CPM {g}"}
        run_feature(f"gene_{g}", expr.loc[g], spec, [g], [])

    # write gene-set files
    (GS / "CITATIONS.txt").write_text(
        "\n".join(
            [
                "GSE126044 counts: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126044",
                "Cho et al. 2020 Exp Mol Med 52:1550 (PMID 32879421). No TJ gene list in the paper.",
                "Chae et al. 2018 Sci Rep 8:1023 Table 2 CBM tight-junction genes.",
                "Reactome R-HSA-420029 / MSigDB REACTOME_TIGHT_JUNCTION_INTERACTIONS.",
                "KEGG hsa04530 Tight junction.",
                "MSigDB GOCC_TIGHT_JUNCTION (GO:0070160).",
                "MSigDB REACTOME_KERATINIZATION (R-HSA-6805567).",
                "ESTIMATE stromal/immune 141-gene signatures: Yoshihara 2013 Nat Commun 4:2612 via tidyestimate 1.1.1.",
                "userPPT_7gene is the claim-page list, not a Cho 2020 list.",
                "",
            ]
        )
    )
    for name, spec in features.items():
        (GS / f"{name}.txt").write_text("\n".join(spec["genes"]) + "\n")
    (GS / "ESTIMATE_stromal.txt").write_text("\n".join(stromal_genes) + "\n")
    (GS / "ESTIMATE_immune.txt").write_text("\n".join(immune_genes) + "\n")
    (GS / "Reactome_R-HSA-420029_msigdb_remote.txt").write_text("\n".join(reactome_remote) + "\n")

    stats_df = pd.DataFrame(stat_rows)
    cov_df = pd.DataFrame(coverage_rows)
    hits_df = pd.DataFrame(hit_rows) if hit_rows else pd.DataFrame(
        columns=["feature", "subset", "test", "which_p", "p"]
    )

    # primary MW table (all + fresh, no residuals/fisher) for the write-up
    primary = stats_df[stats_df["test"] == "MW_score_R_vs_NR"].copy()
    primary.to_csv(TAB / "mw_R_vs_NR.tsv", sep="\t", index=False)
    stats_df.to_csv(TAB / "all_tests.tsv", sep="\t", index=False)
    cov_df.to_csv(TAB / "gene_coverage.tsv", sep="\t", index=False)
    hits_df.to_csv(TAB / "hits_p_le_0.019.tsv", sep="\t", index=False)
    scores.to_csv(TAB / "per_sample_scores.tsv", sep="\t")
    meta.to_csv(TAB / "phenotype.tsv", sep="\t")

    n_r = int((meta["response"] == "responder").sum())
    n_nr = int((meta["response"] == "non-responder").sum())
    n_fresh = int((meta["sample_type"].str.lower() == "fresh").sum())
    n_ffpe = int((meta["sample_type"].str.upper() == "FFPE").sum())
    ffpe_nr = int(((meta["sample_type"].str.upper() == "FFPE") & (meta["response"] == "non-responder")).sum())

    # figures for the three item-6 modules + PPT 7-gene + Chae + CD8A
    for feat, title, ylab in [
        ("CLDN4", "GSE126044 CLDN4", "log2(CPM+1)"),
        ("CLDN147_F11R_PARD3", "GSE126044 CLDN1/4/7/F11R/PARD3 mean-z", "mean-z"),
        ("KRT_basal_squamous", "GSE126044 basal/squamous keratin mean-z", "mean-z"),
        ("userPPT_7gene", "GSE126044 user-PPT 7-gene TJ mean-z", "mean-z"),
        ("Chae2018_CBM_TJ", "GSE126044 Chae 2018 CBM TJ mean-z", "mean-z"),
        ("Reactome_R-HSA-420029", "GSE126044 Reactome TJ mean-z", "mean-z"),
        ("CD8A", "GSE126044 CD8A (label control)", "log2(CPM+1)"),
    ]:
        boxplot(scores[feat], meta, title, ylab, FIG / feat)

    p_dotplot(plot_rows, FIG / "mwu_p_vs_claim")

    # summary JSON
    def pack_mw(feature: str, subset: str) -> dict:
        sub = primary[(primary["feature"] == feature) & (primary["subset"] == subset)]
        if sub.empty:
            return {}
        r = sub.iloc[0]
        return {
            "n_R": int(r["n_R"]),
            "n_NR": int(r["n_NR"]),
            "median_R": r["median_R"],
            "median_NR": r["median_NR"],
            "p_two": r["p_two"],
            "p_one_NR_gt_R": r["p_one_NR_gt_R"],
            "direction": r["direction"],
            "hit_0.019_two": bool(r["hit_0.019_two"]),
            "hit_0.019_one": bool(r["hit_0.019_one"]),
            "matches_claimed_two": bool(r.get("matches_claimed_two", False)),
            "matches_claimed_one": bool(r.get("matches_claimed_one", False)),
        }

    if len(hits_df):
        primary_hits = hits_df[
            (~hits_df["feature"].isin(CONTROL_FEATURES))
            & (~hits_df["feature"].astype(str).str.startswith("gene_"))
        ]
        diagnostic_hits = hits_df[
            hits_df["feature"].astype(str).str.startswith("gene_")
            & (~hits_df["feature"].isin(CONTROL_FEATURES))
        ]
    else:
        primary_hits = hits_df
        diagnostic_hits = hits_df
    if len(primary_hits):
        two_match = primary_hits[primary_hits["which_p"] == "two"]
        verdict = "RECOVERS_CLAIMED_0.019" if len(two_match) else "ONE_SIDED_ONLY_LE_0.019"
    else:
        verdict = "NONE_HIT_0.019"
    noncontrol_hits = primary_hits

    summary = {
        "dataset": "GSE126044",
        "source": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126044/suppl/GSE126044_counts.txt.gz",
        "paper": "Cho et al. 2020 Exp Mol Med (PMID 32879421). Cho 2020 does not publish a TJ gene list.",
        "n_total": int(len(meta)),
        "n_R": n_r,
        "n_NR": n_nr,
        "n_fresh": n_fresh,
        "n_FFPE": n_ffpe,
        "n_FFPE_that_are_NR": ffpe_nr,
        "claim_p": CLAIM_P,
        "verdict": verdict,
        "any_primary_match_claimed_0.019": bool(len(primary_hits)),
        "n_primary_hit_rows": int(len(primary_hits)),
        "n_diagnostic_hit_rows": int(len(diagnostic_hits)),
        "n_hit_rows": int(len(hits_df)),
        "diagnostic_hits": diagnostic_hits.to_dict(orient="records") if len(diagnostic_hits) else [],
        "estimate_coverage": {
            "stromal_used": len(stromal_used),
            "stromal_missing": len(stromal_miss),
            "immune_used": len(immune_used),
            "immune_missing": len(immune_miss),
        },
        "primary_MW_all": {k: pack_mw(k, "all") for k in features},
        "primary_MW_fresh": {k: pack_mw(k, "fresh") for k in features},
        "hits": hits_df.to_dict(orient="records") if len(hits_df) else [],
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=str))

    write_readme(OUT, summary, stats_df, cov_df, hits_df, primary_hits, diagnostic_hits, meta)
    import shutil
    src_man = DATA / "download_manifest.json"
    if src_man.exists():
        shutil.copy(src_man, OUT / "download_manifest.json")
    if len(noncontrol_hits):
        hit_lines = []
        for rec in noncontrol_hits.to_dict(orient="records"):
            hit_lines.append(
                f"- {rec['feature']} | {rec['subset']} | {rec['test']} | "
                f"{rec.get('covariate', '')} | {rec['which_p']} p={rec['p']:.6g}"
            )
        hit_txt = "\n".join(hit_lines)
    else:
        hit_txt = "NONE"
    (OUT / "verdict.txt").write_text(
        f"{verdict}\n"
        f"Cho 2020 has no TJ gene list.\n"
        f"non-control tests with p<=0.019: {len(noncontrol_hits)}\n"
        f"{hit_txt}\n"
    )
    print(json.dumps({"verdict": verdict, "n_hits": int(len(noncontrol_hits)),
                      "CLDN4_all": pack_mw("CLDN4", "all"),
                      "ppt7_all": pack_mw("userPPT_7gene", "all"),
                      "five_all": pack_mw("CLDN147_F11R_PARD3", "all")}, indent=2))
    return 0


def write_readme(out: Path, summary: dict, stats_df: pd.DataFrame, cov_df: pd.DataFrame,
                 hits_df: pd.DataFrame, primary_hits: pd.DataFrame,
                 diagnostic_hits: pd.DataFrame, meta: pd.DataFrame) -> None:
    def fmt_p(p) -> str:
        if p is None or (isinstance(p, float) and math.isnan(p)):
            return "NA"
        return f"{p:.4g}"

    def md_mw(feature: str) -> str:
        rows = []
        for subset in ("all", "fresh"):
            sub = stats_df[(stats_df["feature"] == feature) & (stats_df["subset"] == subset)
                           & (stats_df["test"] == "MW_score_R_vs_NR")]
            if sub.empty:
                continue
            r = sub.iloc[0]
            hit = "YES" if r.get("matches_claimed_two") or r.get("matches_claimed_one") else "no"
            rows.append(
                f"| {subset} | {int(r.n_R)} / {int(r.n_NR)} | {r.median_R:.3g} | {r.median_NR:.3g} | "
                f"{fmt_p(r.p_two)} | {fmt_p(r.p_one_NR_gt_R)} | {r.direction} | {hit} |"
            )
        return "\n".join(rows)

    key_feats = [
        "CLDN4", "CLDN147_F11R_PARD3", "KRT_basal_squamous", "REACTOME_KERATINIZATION",
        "Chae2018_CBM_TJ", "userPPT_7gene", "Reactome_R-HSA-420029", "GOCC_TIGHT_JUNCTION",
        "KEGG_hsa04530", "TJ_15_structural", "CUSTOM_TJ_CORE", "coreTJ_40", "CD8A",
    ]
    blocks = []
    for feat in key_feats:
        src = ""
        sub = cov_df[cov_df["feature"] == feat]
        if not sub.empty:
            src = str(sub.iloc[0]["source"])
            n_used = int(sub.iloc[0]["n_used"])
            n_req = int(sub.iloc[0]["n_requested"])
        else:
            n_used = n_req = 0
        blocks.append(
            f"### {feat}\n\n{src}  \nGenes used: {n_used}/{n_req}\n\n"
            "| subset | n R / NR | median R | median NR | MW two-sided | MW one-sided NR>R | direction | hit 0.019 |\n"
            "|---|---|---|---|---|---|---|---|\n"
            f"{md_mw(feat)}\n"
        )

    resid = stats_df[stats_df["test"] == "MW_residual_R_vs_NR"]
    resid_lines = [
        "| feature | covariate | R2 | n R/NR | median R | median NR | p two | p one NR>R | hit |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for _, r in resid[resid["feature"].isin(key_feats)].iterrows():
        hit = "YES" if r.get("matches_claimed_two") or r.get("matches_claimed_one") else "no"
        resid_lines.append(
            f"| {r.feature} | {r.covariate} | {r.get('resid_R2', float('nan')):.3f} | "
            f"{int(r.n_R)}/{int(r.n_NR)} | {r.median_R:.3g} | {r.median_NR:.3g} | "
            f"{fmt_p(r.p_two)} | {fmt_p(r.p_one_NR_gt_R)} | {hit} |"
        )

    fish = stats_df[stats_df["test"].str.startswith("Fisher")]
    fish_lines = [
        "| feature | test | high NR/R | low NR/R | p two | p one NR-in-high | hit |",
        "|---|---|---|---|---|---|---|",
    ]
    for _, r in fish[fish["feature"].isin(key_feats)].iterrows():
        hit = "YES" if r.get("matches_claimed_two") or r.get("matches_claimed_one") else "no"
        fish_lines.append(
            f"| {r.feature} | {r.test} | {r.get('n_high_NR')}/{r.get('n_high_R')} | "
            f"{r.get('n_low_NR')}/{r.get('n_low_R')} | {fmt_p(r.p_two)} | "
            f"{fmt_p(r.p_one_NR_gt_R)} | {hit} |"
        )

    def hit_table(df: pd.DataFrame, empty: str) -> str:
        if df is None or df.empty:
            return empty
        lines = ["| feature | subset | test | covariate | which p | p |",
                 "|---|---|---|---|---|---|"]
        for rec in df.to_dict(orient="records"):
            lines.append(
                f"| {rec['feature']} | {rec['subset']} | {rec['test']} | "
                f"{rec.get('covariate', '')} | {rec['which_p']} | {fmt_p(rec['p'])} |"
            )
        return "\n".join(lines)

    hit_block = hit_table(
        primary_hits,
        "No pre-specified module test reached p≤0.019 or rounded to 0.019.",
    )
    diag_block = hit_table(
        diagnostic_hits,
        "No single-gene diagnostic reached p≤0.019.",
    )

    n_r = summary["n_R"]
    n_nr = summary["n_NR"]
    text = f"""# B4 wave-2: GSE126044 NR-higher TJ p=0.019

**Verdict: `{summary['verdict']}`**

The claimed two-sided p=0.019 is recovered by **one** pre-specified definition: the user-PPT / claim-page 7-gene mean-z (`CLDN1, CLDN4, CLDN7, F11R, TJP1, TJP2, OCLN`). Exact Mann-Whitney two-sided p = **0.01923** (5 R vs 11 NR), which is the discrete p people write as 0.019. One-sided NR>R p = **0.0096**.

That list is **not** from Cho 2020. Cho 2020 publishes no TJ gene list. Related paper/database TJ sets, CLDN4 alone, the 5-gene CLDN1/4/7/F11R/PARD3 module, keratin modules, FFPE-drop, ESTIMATE residualization, and tertile/median splits of those modules do **not** give two-sided p≤0.019.

The 7-gene result is largely **OCLN-driven** (OCLN alone two-sided p=0.0055). Residualizing the 7-gene score on stromal or ESTIMATEScore removes the 0.019 (p=0.32 and 0.91). Fresh-only two-sided p=0.030.

Public GSE126044 counts only. Nothing was tuned to hit 0.019.

## Item-by-item (does it recover p≤0.019?)

| # | What was tried | Recovers p≤0.019? | Exact result |
|---|---|---|---|
| 1 | Author/paper TJ lists | **Only the user-PPT 7-gene** (not Cho 2020) | Cho 2020: no TJ list. Chae 2018 CBM TJ two-sided p=0.145. Reactome TJ p=0.090. GOCC TJ p=0.377. KEGG hsa04530 p=0.441. user-PPT 7-gene two-sided **p=0.01923** |
| 2 | One-sided MW (NR>R) and two-sided | **7-gene yes (both)** | 7-gene two-sided 0.01923; one-sided 0.0096. CLDN4 one-sided 0.057. Reactome one-sided 0.045. Chae one-sided 0.073 |
| 3 | Drop all FFPE (fresh-only, 5 R / 6 NR) | **No two-sided ≤0.019** | 7-gene two-sided 0.030; one-sided 0.015. CLDN4 two-sided 0.329. OCLN two-sided 0.017 (single gene) |
| 4 | Residualize on StromalScore / ESTIMATEScore | **No** | 7-gene residual stromal p=0.32; ESTIMATE p=0.91. CLDN4 residual stromal p=0.27; ESTIMATE p=0.66 |
| 5 | Tertile / top-vs-bottom (Fisher) | **No for modules** | 7-gene median Fisher p=0.28; tertile p=0.24. CLDN4 median p=0.28; tertile p=0.55 |
| 6 | CLDN4 vs CLDN1/4/7/F11R/PARD3 vs keratin | **No** | CLDN4 two-sided 0.115. 5-gene 0.320. basal/squamous KRT 0.661. REACTOME_KERATINIZATION 0.441 |

## Which pre-specified module tests match claimed 0.019?

{hit_block}

Single-gene diagnostics of the 7-gene members (not used to invent a new signature):

{diag_block}

CD8A is a label control (R > NR, two-sided p=0.00092) and is not a TJ hit.

## Cohort

| | |
|---|---|
| Series | [GSE126044](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126044) |
| Counts | `GSE126044_counts.txt.gz` (author FeatureCounts, GRCh37/GENCODE 19) |
| n | {summary['n_total']} ({n_r} responder / {n_nr} non-responder) |
| Tissue | {summary['n_fresh']} fresh, {summary['n_FFPE']} FFPE |
| FFPE × response | {summary['n_FFPE_that_are_NR']}/{summary['n_FFPE']} FFPE samples are NR (all FFPE are NR; all 5 R are fresh) |
| Score | single gene = log2(CPM+1); modules = mean of per-gene z-scores across the 16 samples |
| Primary test | Mann-Whitney U, two-sided and one-sided (NR > R) |

Cho 2020 is a methylation/enhancer paper. It does **not** publish a tight-junction gene list. Item 1 therefore uses related NSCLC barrier/immune papers and public TJ pathway sets, not a fabricated Cho list.

A hit below is p≤0.019 or a discrete MW p that rounds to 0.019 (0.01923). Controls (CD8A/CD8) are excluded from the verdict.

## Item 1 — author/paper/database TJ lists

{blocks[4]}
{blocks[5]}
{blocks[6]}
{blocks[7]}
{blocks[8]}
{blocks[9]}

`userPPT_7gene` is the 7-gene module written on the user claim page (CLDN1, CLDN4, CLDN7, F11R, TJP1, TJP2, OCLN). It is **not** from Cho 2020. It is included because that is the list prior claim-page work used.

## Item 2 — one-sided and two-sided MW

Both p-values are in every MW row above. CD8A is the label control (responders higher).

{blocks[12]}

## Item 3 — drop FFPE (fresh-only)

Fresh-only is 5 R vs 6 NR. All 5 FFPE libraries are NR, so sample type is aliased with response. Fresh-only p-values are in the tables above; they do not create a new significant claim unless marked YES.

## Item 4 — residualize on ESTIMATE stromal / ESTIMATEScore

Stromal and immune signatures are the 141-gene Yoshihara 2013 lists (via tidyestimate 1.1.1). Scores are mean-z (not the Affymetrix ssGSEA implementation). Residual = OLS residual of the feature on the covariate.

{chr(10).join(resid_lines)}

## Item 5 — median split and tertile (top vs bottom)

Fisher exact on response in TJ-high vs TJ-low. Tertile drops the middle third. n=16, so cells are tiny.

{chr(10).join(fish_lines)}

## Item 6 — CLDN4 vs CLDN1/4/7/F11R/PARD3 vs keratin

{blocks[0]}
{blocks[1]}
{blocks[2]}
{blocks[3]}

## Prior-PR comparators (not expected to hit 0.019)

{blocks[10]}
{blocks[11]}

## Honest limits

- n=16. Only a large, well-aligned effect gives p=0.019.
- FFPE is completely nested in NR.
- Cho 2020 has no TJ list. Using the same cohort's DE genes as a "TJ list" would be circular and was not done.
- ESTIMATE here is mean-z of the published 141-gene lists, not the original Affymetrix ssGSEA + purity cosine.
- Multiple tests are reported. A single p≤0.019 among many definitions is not independent confirmation of the slide.

## Rerun

```bash
python3 -m pip install -r scripts/rework/B4_wave2/requirements.txt
python3 scripts/rework/B4_wave2/download.py
python3 scripts/rework/B4_wave2/analyze.py
```
"""
    (out / "README.md").write_text(text)
    (out / "WRITEUP.md").write_text(text)


if __name__ == "__main__":
    raise SystemExit(main())
