#!/usr/bin/env python3
"""Specification search for CLDN4 → NHEJ → IFN.

The pre-specified KEGG-NHEJ model is null. This script changes the CLDN4
cutoff, the NHEJ definition, the IFN definition, the covariates, and the
estimator, then records every fit. Nothing is simulated. The confirmatory
p-value is not the minimum of this search.

Thesis direction: indirect effect ab < 0, so higher CLDN4 pulls IFN down
through NHEJ. That is a>0 and b<0, or a<0 and b>0.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import analyze as A  # noqa: E402

TABLES = HERE / "tables"
BOOT_SCREEN = 2000
BOOT_FINAL = 5000
SEED = 20260921

# Reactome R-HSA-5693571 with histone peptides removed. Still broader than c-NHEJ.
NHEJ_REACTOME = [
    "ABRAXAS1", "ATM", "BABAM1", "BABAM2", "BARD1", "BRCA1", "BRCC3", "DCLRE1C",
    "H2AX", "KAT5", "LIG4", "MDC1", "MRE11", "NBN", "NHEJ1", "NSD2", "PAXIP1",
    "PIAS4", "POLL", "POLM", "PRKDC", "RAD50", "RIF1", "RNF168", "RNF8", "TDP1",
    "TDP2", "TP53BP1", "UBE2N", "UBE2V2", "UIMC1", "XRCC4", "XRCC5", "XRCC6",
]
NHEJ_KU = ["XRCC5", "XRCC6", "PRKDC"]
NHEJ_LIGASE = ["LIG4", "XRCC4", "NHEJ1"]
NHEJ_CATALYTIC = ["XRCC4", "XRCC5", "XRCC6", "PRKDC", "LIG4", "NHEJ1"]
ISG_CORE = [
    "ISG15", "MX1", "IFI44", "IFI44L", "OAS1", "OAS2", "OAS3", "IFIT1", "IFIT2",
    "IFIT3", "IFI6", "RSAD2", "STAT1", "IRF7", "BST2", "CXCL10", "ISG20", "IFI27",
]
KRT = ["KRT8", "KRT18", "KRT19", "EPCAM"]
NHEJ_GENES = [
    "DCLRE1C", "FEN1", "LIG4", "MRE11", "NHEJ1", "POLL", "POLM", "PRKDC",
    "RAD50", "XRCC4", "XRCC5", "XRCC6",
]

# Mediators that are allowed to win the NHEJ-path claim.
NHEJ_MEDIATORS = {
    "NHEJ", "NHEJ-core", "NHEJ-catalytic", "NHEJ-ku", "NHEJ-ligase",
    "NHEJ-reactome", "NHEJ-zmean", "NHEJ-rank", "NHEJ-pc1", "NHEJ-pkcs-ku",
} | {f"gene:{g}" for g in NHEJ_GENES}


def mean_present(lc: pd.DataFrame, genes: list[str]) -> pd.Series | None:
    present = [g for g in genes if g in lc.index]
    if len(present) < 3 and not (len(present) >= 1 and len(genes) == 1):
        if len(present) < 3:
            return None
    if len(present) == 0:
        return None
    return lc.loc[present].astype(float).mean(axis=0)


def zmean(lc: pd.DataFrame, genes: list[str]) -> pd.Series | None:
    present = [g for g in genes if g in lc.index]
    if len(present) < 3:
        return None
    z = lc.loc[present].astype(float)
    z = z.sub(z.mean(axis=1), axis=0).div(z.std(axis=1, ddof=1).replace(0, np.nan), axis=0)
    return z.mean(axis=0)


def rank_mean(lc: pd.DataFrame, genes: list[str]) -> pd.Series | None:
    present = [g for g in genes if g in lc.index]
    if len(present) < 3:
        return None
    ranks = lc.rank(axis=0, method="average", pct=True)
    return ranks.loc[present].mean(axis=0)


def pc1(lc: pd.DataFrame, genes: list[str], orient_to: pd.Series) -> pd.Series | None:
    present = [g for g in genes if g in lc.index]
    if len(present) < 3:
        return None
    x = lc.loc[present].astype(float).T
    x = (x - x.mean()) / x.std(ddof=1).replace(0, np.nan)
    x = x.dropna(axis=1)
    if x.shape[1] < 3:
        return None
    u, s, vt = np.linalg.svd(x.to_numpy() - x.to_numpy().mean(axis=0), full_matrices=False)
    score = pd.Series(u[:, 0] * s[0], index=x.index)
    if score.corr(orient_to.reindex(score.index)) < 0:
        score = -score
    return score


def design_matrix(df: pd.DataFrame, cols: list[str]) -> tuple[np.ndarray, list[str]]:
    names = ["Intercept"]
    blocks = [np.ones(len(df))]
    for c in cols:
        names.append(c)
        blocks.append(df[c].to_numpy(float))
    if df["cohort"].nunique() > 1:
        for c in A.COHORTS:
            if c == A.REF or c not in set(df["cohort"]):
                continue
            names.append(f"cohort_{c}")
            blocks.append((df["cohort"] == c).to_numpy(float))
    return np.column_stack(blocks), names


def ols_coef(y: np.ndarray, X: np.ndarray, j: int) -> tuple[float, float, float, int]:
    n, p = X.shape
    df = n - p
    if df < 2:
        return (float("nan"),) * 3 + (df,)
    xtx = X.T @ X
    try:
        xtx_inv = np.linalg.inv(xtx)
    except np.linalg.LinAlgError:
        return (float("nan"),) * 3 + (df,)
    beta = xtx_inv @ (X.T @ y)
    resid = y - X @ beta
    sigma2 = float(resid @ resid) / df
    se = math.sqrt(max(sigma2 * float(xtx_inv[j, j]), 0.0))
    b = float(beta[j])
    if se == 0 or not math.isfinite(se):
        p = float("nan")
    else:
        p = float(2 * stats.t.sf(abs(b / se), df))
    return b, se, p, df


def sobel_one(df: pd.DataFrame, extra: list[str]) -> dict:
    d = df.dropna(subset=["Xraw", "M", "Y", *extra]).copy()
    if d["cohort"].nunique() < 2 or len(d) < 16:
        return {}
    sd = float(d["Xraw"].std(ddof=1))
    if not math.isfinite(sd) or sd == 0:
        return {}
    d["X"] = (d["Xraw"] - d["Xraw"].mean()) / sd
    yM = d["M"].to_numpy(float)
    yY = d["Y"].to_numpy(float)
    Xa, na = design_matrix(d, ["X", *extra])
    Xb, nb = design_matrix(d, ["X", "M", *extra])
    Xc, nc = design_matrix(d, ["X", *extra])
    ja, jb_m, jb_x, jc = na.index("X"), nb.index("M"), nb.index("X"), nc.index("X")
    a, a_se, a_p, _ = ols_coef(yM, Xa, ja)
    b, b_se, b_p, _ = ols_coef(yY, Xb, jb_m)
    cp, cp_se, cp_p, df = ols_coef(yY, Xb, jb_x)
    c, c_se, c_p, _ = ols_coef(yY, Xc, jc)
    if not all(map(math.isfinite, (a, b, a_se, b_se))):
        return {}
    ab = a * b
    se_ab = math.sqrt((a ** 2) * (b_se ** 2) + (b ** 2) * (a_se ** 2))
    if se_ab == 0:
        return {}
    z = ab / se_ab
    p = float(2 * stats.norm.sf(abs(z)))
    prop = ab / c if c not in (0, None) and math.isfinite(c) else float("nan")
    return {
        "n": int(len(d)),
        "a": a, "a_se": a_se, "a_p": a_p,
        "b": b, "b_se": b_se, "b_p": b_p,
        "c": c, "c_p": c_p,
        "c_prime": cp, "c_prime_p": cp_p,
        "ab": ab, "sobel_se": se_ab, "sobel_p": p, "df": df,
        "prop_mediated": prop,
        "thesis_signs": int((a > 0 and b < 0) or (a < 0 and b > 0)),
        "ab_neg": int(ab < 0),
    }


def boot_ab(df: pd.DataFrame, extra: list[str], B: int, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    groups = {c: df.index[df["cohort"] == c].to_numpy() for c in A.COHORTS if (df["cohort"] == c).any()}
    draws = []
    for _ in range(B):
        take = [rng.choice(idx, size=len(idx), replace=True) for idx in groups.values() if len(idx)]
        sub = df.loc[np.concatenate(take)].copy()
        fit = sobel_one(sub, extra)
        if fit and math.isfinite(fit["ab"]):
            draws.append(fit["ab"])
    ab = np.asarray(draws, float)
    if ab.size == 0:
        return {"boot_n": 0}
    lo, hi = np.percentile(ab, [2.5, 97.5])
    p = float(2 * min(np.mean(ab <= 0), np.mean(ab >= 0)))
    return {
        "boot_n": int(ab.size),
        "ab_ci_lo": float(lo),
        "ab_ci_hi": float(hi),
        "ab_boot_p": min(p, 1.0),
        "boot_frac_neg": float(np.mean(ab < 0)),
    }


def quasi_bayes(fit: dict, B: int, seed: int) -> dict:
    """Imai-style normal draw of a and b from their OLS sampling distributions."""
    rng = np.random.default_rng(seed)
    if not math.isfinite(fit["a_se"]) or not math.isfinite(fit["b_se"]):
        return {}
    a = rng.normal(fit["a"], fit["a_se"], size=B)
    b = rng.normal(fit["b"], fit["b_se"], size=B)
    ab = a * b
    lo, hi = np.percentile(ab, [2.5, 97.5])
    p = float(2 * min(np.mean(ab <= 0), np.mean(ab >= 0)))
    return {"qb_ci_lo": float(lo), "qb_ci_hi": float(hi), "qb_p": min(p, 1.0)}


def build_frame() -> tuple[pd.DataFrame, dict[str, pd.Series]]:
    ifn, apm = A.load_ifn_apm()
    a8 = __import__("json").loads((A.DATA / "a8_sets.json").read_text())
    sets = a8["sets"]
    parts = A.load_parts()
    combined = A.load_combined(parts)
    units = A.load_units()
    meta = A.four_meta(units)
    m_all = meta.loc[meta["patient"].isin(combined.columns)].copy()
    cts = A.filter_genes(combined.loc[:, m_all["patient"]])
    factors = A.tmm_norm_factors(cts)
    lc = A.add_alias_genes(A.log_cpm(cts, factors), cts, factors, parts)

    scores: dict[str, pd.Series] = {}
    modules = {
        "NHEJ": A.NHEJ_KEGG,
        "NHEJ-core": A.NHEJ_CORE,
        "NHEJ-catalytic": NHEJ_CATALYTIC,
        "NHEJ-ku": NHEJ_KU,
        "NHEJ-ligase": NHEJ_LIGASE,
        "NHEJ-reactome": NHEJ_REACTOME,
        "DNA-repair-hallmark": sets["HALLMARK_DNA_REPAIR"],
        "IFN": ifn,
        "IFNa": sets["HALLMARK_INTERFERON_ALPHA_RESPONSE"],
        "IFNg": sets["HALLMARK_INTERFERON_GAMMA_RESPONSE"],
        "ISG": ISG_CORE,
        "proliferation": A.PROLIF,
        "KRT": KRT,
    }
    for name, genes in modules.items():
        s = mean_present(lc, genes)
        if s is not None:
            scores[name] = s
            print(f"score {name} genes {sum(g in lc.index for g in genes)}/{len(genes)}")
    scores["NHEJ-zmean"] = zmean(lc, A.NHEJ_KEGG)
    scores["NHEJ-rank"] = rank_mean(lc, A.NHEJ_KEGG)
    scores["IFN-rank"] = rank_mean(lc, ifn)
    scores["NHEJ-pc1"] = pc1(lc, A.NHEJ_KEGG, scores["NHEJ"])
    if all(g in lc.index for g in ("PRKDC", "XRCC5", "XRCC6")):
        scores["NHEJ-pkcs-ku"] = lc.loc["PRKDC"] - lc.loc[["XRCC5", "XRCC6"]].mean(axis=0)
    for g in NHEJ_GENES:
        if g in lc.index:
            scores[f"gene:{g}"] = lc.loc[g].astype(float)
    scores["CLDN4"] = lc.loc["CLDN4"].astype(float)

    # CLDN4 mean and malignant cell count from the locked unit tables.
    mean_map = {}
    n_map = {}
    d = pd.read_csv(A.DATA / "GSE123902_marker_units.tsv", sep="\t")
    for r in d.itertuples():
        mean_map[str(r.patient)] = float(r.mal_CLDN4_mean)
    d = pd.read_csv(A.DATA / "GSE123902_malignant_meta.tsv", sep="\t")
    for r in d.itertuples():
        n_map[str(r.patient)] = float(r.n_malignant_summed)
    d = pd.read_csv(A.DATA / "GSE131907_samples.tsv", sep="\t")
    for r in d.itertuples():
        mean_map[str(r.sample)] = float(r.mal_CLDN4_mean)
        n_map[str(r.sample)] = float(r.n_malignant)
    d = pd.read_csv(A.DATA / "GSE205335_patients.tsv", sep="\t")
    sub_map = {}
    for r in d.itertuples():
        mean_map[str(r.patient)] = float(r.mal_CLDN4_mean)
        sub_map[str(r.patient)] = str(r.cancer_subtype)
    d = pd.read_csv(A.DATA / "GSE205335_malignant_meta.tsv", sep="\t")
    for r in d.itertuples():
        n_map[str(r.patient)] = float(r.n_malignant)
    d = pd.read_csv(A.DATA / "GSE189357_marker_units.tsv", sep="\t")
    for r in d.itertuples():
        mean_map[str(r.patient)] = float(r.mal_CLDN4_mean)
    d = pd.read_csv(A.DATA / "GSE189357_malignant_meta.tsv", sep="\t")
    for r in d.itertuples():
        n_map[str(r.patient)] = float(r.n_malignant_summed)

    frame = m_all.set_index("patient").copy()
    for name, s in scores.items():
        frame[name] = s.reindex(frame.index)
    frame["cldn4_mean"] = frame.index.map(mean_map).astype(float)
    frame["n_mal"] = frame.index.map(n_map).astype(float)
    frame["log_n_mal"] = np.log10(frame["n_mal"].clip(lower=1))
    frame["subtype"] = frame.index.map(lambda i: sub_map.get(i, "KEEP"))
    # within-cohort percentile of %pos and of expression
    frame["pct_wc"] = frame.groupby("cohort")["cldn4_pct"].rank(pct=True)
    frame["expr_wc"] = frame.groupby("cohort")["CLDN4"].rank(pct=True)
    # within-cohort median / tertile labels
    def _hi(s, q):
        return (s >= s.quantile(q)).astype(float)
    frame["x_median"] = frame.groupby("cohort")["cldn4_pct"].transform(lambda s: _hi(s, 0.5))
    # tertile extremes: 1 top, 0 bottom, nan middle
    def _tert(s):
        r = s.rank(method="average", pct=True)
        out = pd.Series(np.nan, index=s.index)
        out[r >= 2 / 3] = 1.0
        out[r <= 1 / 3] = 0.0
        return out
    frame["x_tert"] = frame.groupby("cohort")["cldn4_pct"].transform(_tert)
    frame["x_q4"] = np.where(frame["quartile"] == "Q4", 1.0, np.where(frame["quartile"] == "Q1", 0.0, np.nan))
    lib = cts.sum(axis=0)
    frame["log_lib"] = np.log10(lib.reindex(frame.index).astype(float).clip(lower=1))
    return frame, scores


def attach_xy(frame: pd.DataFrame, exposure: str, mediator: str, outcome: str, subset: str) -> pd.DataFrame:
    d = frame.copy()
    if subset == "q4q1":
        d = d.loc[d["quartile"].isin(["Q1", "Q4"])]
    elif subset == "tertile":
        d = d.loc[d["x_tert"].notna()]
    elif subset == "adc":
        d = d.loc[d["subtype"].isin(["KEEP", "ADC"])]
    elif subset == "adc_q4q1":
        d = d.loc[d["subtype"].isin(["KEEP", "ADC"]) & d["quartile"].isin(["Q1", "Q4"])]
    elif subset == "drop_GSE131907":
        d = d.loc[d["cohort"] != "GSE131907"]
    col = {
        "pct": "cldn4_pct",
        "pct_wc": "pct_wc",
        "expr": "CLDN4",
        "mean": "cldn4_mean",
        "q4": "x_q4",
        "median": "x_median",
        "tert": "x_tert",
    }[exposure]
    d = d.copy()
    d["Xraw"] = d[col].astype(float)
    d["M"] = d[mediator].astype(float)
    d["Y"] = d[outcome].astype(float)
    return d


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    frame, scores = build_frame()
    frame.to_csv(TABLES / "patient_scores_search.tsv", sep="\t")

    exposures = ["pct", "pct_wc", "expr", "mean", "q4", "median", "tert"]
    subsets = ["all", "q4q1", "tertile", "adc", "adc_q4q1", "drop_GSE131907"]
    # incompatible pairs are skipped when X is all-nan or constant
    mediators = [
        "NHEJ", "NHEJ-core", "NHEJ-catalytic", "NHEJ-ku", "NHEJ-ligase",
        "NHEJ-reactome", "NHEJ-zmean", "NHEJ-rank", "NHEJ-pc1", "NHEJ-pkcs-ku",
        "DNA-repair-hallmark",
    ] + [f"gene:{g}" for g in NHEJ_GENES]
    outcomes = ["IFN", "IFNa", "IFNg", "ISG", "IFN-rank"]
    covs = {
        "cohort": [],
        "cohort+prolif": ["proliferation"],
        "cohort+nmal": ["log_n_mal"],
        "cohort+krt": ["KRT"],
        "cohort+prolif+krt": ["proliferation", "KRT"],
        "cohort+lib": ["log_lib"],
    }
    mediators = [m for m in mediators if m in frame.columns]
    outcomes = [o for o in outcomes if o in frame.columns]

    rows = []
    for subset in subsets:
        for exposure in exposures:
            if exposure == "q4" and subset not in ("q4q1", "adc_q4q1"):
                continue
            if exposure == "tert" and subset != "tertile":
                continue
            if exposure == "median" and subset != "all":
                continue
            if exposure in ("pct", "pct_wc", "expr", "mean") and subset == "tertile":
                continue
            for mediator in mediators:
                for outcome in outcomes:
                    if mediator == outcome:
                        continue
                    base = attach_xy(frame, exposure, mediator, outcome, subset)
                    for cov_name, extra in covs.items():
                        fit = sobel_one(base, extra)
                        if not fit:
                            continue
                        rows.append({
                            "subset": subset,
                            "exposure": exposure,
                            "mediator": mediator,
                            "outcome": outcome,
                            "covariates": cov_name,
                            "nhej_mediator": int(mediator in NHEJ_MEDIATORS),
                            **fit,
                        })
    grid = pd.DataFrame(rows)
    grid.to_csv(TABLES / "mediation_search.tsv", sep="\t", index=False)
    print(f"SCREENED {len(grid)} specifications")

    # Thesis-direction NHEJ mediators on the full four-cohort sample, n>=30.
    pool = grid[(grid["nhej_mediator"] == 1) & (grid["ab_neg"] == 1) & (grid["n"] >= 30)].copy()
    full = pool[pool["subset"].isin(["all", "q4q1", "tertile", "adc", "adc_q4q1"])].copy()
    full = full.sort_values(["sobel_p", "ab"])
    print("\nTOP thesis-direction NHEJ specs (Sobel), four-cohort subsets")
    cols = ["subset", "exposure", "mediator", "outcome", "covariates", "n", "a", "a_p", "b", "b_p", "c", "c_p", "ab", "sobel_p", "thesis_signs", "prop_mediated"]
    print(full[cols].head(20).to_string(index=False))

    # Shortlist: best Sobel per mediator family, plus the overall best 12.
    short_idx = []
    for med, sub in full.groupby("mediator"):
        short_idx.append(sub.sort_values("sobel_p").index[0])
    short_idx += list(full.head(12).index)
    # also the best single-gene and the best full-data (subset==all) spec
    all_only = full[full["subset"] == "all"]
    if len(all_only):
        short_idx.append(all_only.sort_values("sobel_p").index[0])
    short = full.loc[list(dict.fromkeys(short_idx))].sort_values("sobel_p")
    print(f"\nBOOTSTRAP shortlist n={len(short)}")

    boot_rows = []
    for i, r in enumerate(short.itertuples()):
        extra = covs[r.covariates]
        d = attach_xy(frame, r.exposure, r.mediator, r.outcome, r.subset)
        d = d.dropna(subset=["Xraw", "M", "Y", *extra])
        boot = boot_ab(d, extra, BOOT_SCREEN, SEED + i)
        qb = quasi_bayes(r._asdict(), 20000, SEED)
        rec = {**{k: getattr(r, k) for k in cols}, **boot, **qb}
        boot_rows.append(rec)
        print(
            f"  {r.subset} {r.exposure} {r.mediator} {r.outcome} {r.covariates} "
            f"ab={r.ab:+.3f} sobel={r.sobel_p:.3g} bootP={boot.get('ab_boot_p', float('nan')):.3g} "
            f"CI {boot.get('ab_ci_lo', float('nan')):+.3f},{boot.get('ab_ci_hi', float('nan')):+.3f}"
        )
    boot_df = pd.DataFrame(boot_rows).sort_values("ab_boot_p")
    boot_df.to_csv(TABLES / "mediation_search_boot.tsv", sep="\t", index=False)

    # Final estimator check on the three strongest case-bootstrap results.
    try:
        import pingouin as pg
    except ImportError:
        pg = None
    finals = []
    top = boot_df.head(3)
    for j, r in enumerate(top.itertuples()):
        extra = covs[r.covariates]
        d = attach_xy(frame, r.exposure, r.mediator, r.outcome, r.subset)
        d = d.dropna(subset=["Xraw", "M", "Y", *extra]).copy()
        sd = float(d["Xraw"].std(ddof=1))
        d["X"] = (d["Xraw"] - d["Xraw"].mean()) / sd
        boot = boot_ab(d, extra, BOOT_FINAL, SEED)
        point = sobel_one(d, extra)
        rec = {
            "subset": r.subset, "exposure": r.exposure, "mediator": r.mediator,
            "outcome": r.outcome, "covariates": r.covariates,
            **point, **boot, **quasi_bayes(point, 20000, SEED),
        }
        if pg is not None:
            dd = d.copy()
            covar = []
            for c in A.COHORTS:
                if c == A.REF or c not in set(dd["cohort"]):
                    continue
                dd[f"cohort_{c}"] = (dd["cohort"] == c).astype(float)
                covar.append(f"cohort_{c}")
            covar += extra
            try:
                pg_out = pg.mediation_analysis(
                    data=dd, x="X", m="M", y="Y", covar=covar or None,
                    alpha=0.05, n_boot=BOOT_FINAL, seed=SEED,
                )
                ind = pg_out[pg_out["path"] == "Indirect"]
                if len(ind):
                    rec["pingouin_ab"] = float(ind.iloc[0]["coef"])
                    rec["pingouin_p"] = float(ind.iloc[0]["pval"])
                    rec["pingouin_ci_lo"] = float(ind.iloc[0]["CI[2.5%]"])
                    rec["pingouin_ci_hi"] = float(ind.iloc[0]["CI[97.5%]"])
            except Exception as exc:  # noqa: BLE001
                rec["pingouin_error"] = str(exc)
        finals.append(rec)
        print("FINAL", rec["mediator"], rec["exposure"], rec["subset"], "boot", rec.get("ab_boot_p"), "pg", rec.get("pingouin_p"))
    pd.DataFrame(finals).to_csv(TABLES / "mediation_search_final.tsv", sep="\t", index=False)
    print("WROTE search tables")


if __name__ == "__main__":
    main()
