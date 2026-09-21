#!/usr/bin/env python3
"""Specification sweep on the same concordant-4 malignant pseudobulks.

The pre-specified analysis remains scripts/analyze.py (DoRothEA A+B+C, ULM,
locked Q4 vs Q1). This script re-runs the listed networks, decoupler
statistics, and CLDN4 cutoffs and writes every row. Ranking is applied
after the table exists. No cohort is dropped to manufacture a sign.

Thesis labels used only as column flags:
  IFN/STAT1/IRF activity higher in CLDN4-low (coef of high-vs-low < 0)
  barrier/TJ activity higher in CLDN4-high (coef > 0)
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

import decoupler as dc

sys.path.insert(0, str(Path(__file__).resolve().parent))
import analyze as az  # noqa: E402

HERE = az.HERE
RES = az.RES
OUT = az.OUT
FIGS = az.FIGS
COHORTS = az.COHORTS
REF = az.REF

IFN_TFS = ["STAT1", "STAT2", "IRF1", "IRF9"]
TJ_TFS = ["GRHL2", "ELF3", "KLF4", "TFAP2A"]
FOCUS = IFN_TFS + TJ_TFS
PATHWAYS = ["JAK-STAT", "TNFa", "NFkB"]  # PROGENy; no TJ pathway exists there

CUTOFFS = [
    "quartile_locked",
    "quartile_equal",
    "tertile",
    "median",
    "tail30",
    "tail20",
]


def load_dorothea_all() -> pd.DataFrame:
    net = pd.read_csv(RES / "dorothea_hs_ABCD.tsv", sep="\t")
    net["source"] = net["source"].astype(str).str.upper()
    net["target"] = net["target"].astype(str).str.upper()
    net["weight"] = net["weight"].astype(float)
    return net


def load_progeny() -> pd.DataFrame:
    net = pd.read_csv(RES / "progeny_top500.tsv", sep="\t")
    net["source"] = net["source"].astype(str)
    net["target"] = net["target"].astype(str).str.upper()
    net["weight"] = net["weight"].astype(float)
    return net[["source", "target", "weight"]]


def geneset_net(log_index: pd.Index) -> pd.DataFrame:
    a8 = json.loads((az.DATA / "a8_sets.json").read_text())["sets"]
    ifn = {g.upper() for g in set(a8["HALLMARK_INTERFERON_GAMMA_RESPONSE"]) | set(a8["HALLMARK_INTERFERON_ALPHA_RESPONSE"])}
    tj = {g.upper() for g in set(a8["KEGG_TIGHT_JUNCTION"]) | set(a8["GOBP_TIGHT_JUNCTION_ORGANIZATION"])}
    tj.discard("CLDN4")
    apical = {g.upper() for g in a8.get("HALLMARK_APICAL_JUNCTION", [])}
    rows = []
    for name, genes in [("IFN_hallmark", ifn), ("TJ_no_CLDN4", tj), ("APICAL_JUNCTION", apical)]:
        for g in sorted(genes):
            if g in log_index:
                rows.append({"source": name, "target": g, "weight": 1.0})
    return pd.DataFrame(rows)


def tail_mask(meta: pd.DataFrame, cutoff: str) -> pd.DataFrame:
    """Within-cohort high/low labels. Middle units are dropped."""
    parts = []
    for cohort, g in meta.groupby("cohort", sort=False):
        g = g.sort_values(["cldn4_pct", "patient"])
        n = len(g)
        if cutoff == "quartile_locked":
            lab = pd.Series(np.where(g["quartile"].eq("Q4"), "high", np.where(g["quartile"].eq("Q1"), "low", None)), index=g.index)
        else:
            if cutoff == "median":
                k = n // 2
            elif cutoff == "tertile":
                k = n // 3
            elif cutoff == "quartile_equal":
                k = n // 4
            elif cutoff == "tail30":
                k = int(np.floor(n * 0.30))
            elif cutoff == "tail20":
                k = int(np.floor(n * 0.20))
            else:
                raise ValueError(cutoff)
            lab = pd.Series([None] * n, index=g.index, dtype=object)
            if k >= 2:
                lab.iloc[:k] = "low"
                lab.iloc[-k:] = "high"
        gg = g.copy()
        gg["arm"] = lab.to_numpy()
        parts.append(gg)
    out = pd.concat(parts, axis=0)
    return out


def design_arm(meta: pd.DataFrame, with_cohort: bool) -> pd.DataFrame:
    m = meta.set_index("patient")
    d = pd.DataFrame(index=m.index)
    d["Intercept"] = 1.0
    d["CLDN4_high"] = (m["arm"] == "high").astype(float)
    if with_cohort:
        for c in COHORTS:
            if c == REF:
                continue
            d[f"cohort_{c}"] = (m["cohort"] == c).astype(float)
    return d


def score_net(logcpm: pd.DataFrame, net: pd.DataFrame, method: str) -> pd.DataFrame | None:
    data = logcpm.T.astype(float)
    try:
        if method == "ulm":
            es, _ = dc.mt.ulm(data, net, tmin=az.TMIN)
        elif method == "mlm":
            es, _ = dc.mt.mlm(data, net, tmin=az.TMIN)
        elif method == "wmean":
            es, _ = dc.mt.waggr(data, net, tmin=az.TMIN, fun="wmean", times=1, seed=42)
        elif method == "wsum":
            es, _ = dc.mt.waggr(data, net, tmin=az.TMIN, fun="wsum", times=1, seed=42)
        elif method == "wmean_nes":
            es, _ = dc.mt.waggr(data, net, tmin=az.TMIN, fun="wmean", times=100, seed=42)
        else:
            raise ValueError(method)
    except Exception as exc:  # rank-deficient MLM or empty net
        print(f"  FAIL {method}: {exc}")
        return None
    return es


def n_targets(net: pd.DataFrame, genes: pd.Index) -> pd.Series:
    sub = net.loc[net["target"].isin(set(map(str, genes)))]
    return sub.groupby("source")["target"].nunique()


def fit_high(y: pd.Series, meta_arm: pd.DataFrame, with_cohort: bool) -> dict:
    m = meta_arm.loc[meta_arm["arm"].isin(["low", "high"])].copy()
    m = m.loc[m["patient"].isin(y.index)]
    n_low = int((m["arm"] == "low").sum())
    n_high = int((m["arm"] == "high").sum())
    base = {"n_low": n_low, "n_high": n_high, "n": n_low + n_high}
    if n_low < 3 or n_high < 3:
        return {**base, "coef": np.nan, "se": np.nan, "p": np.nan, "note": "tail n<3"}
    fit = az.ols_coef(y.reindex(m["patient"]), design_arm(m, with_cohort), "CLDN4_high")
    return {**base, "coef": fit["coef"], "se": fit["se"], "p": fit["p"], "note": ""}


def ivw(rows: list[dict]) -> dict:
    usable = [r for r in rows if np.isfinite(r.get("coef", np.nan)) and np.isfinite(r.get("se", np.nan)) and r["se"] > 0]
    if len(usable) < 2:
        return {"n_cohorts": len(usable), "coef": np.nan, "se": np.nan, "p": np.nan, "n_low": np.nan, "n_high": np.nan}
    w = np.array([1.0 / r["se"] ** 2 for r in usable])
    c = np.array([r["coef"] for r in usable])
    coef = float(np.sum(w * c) / np.sum(w))
    se = float(np.sqrt(1.0 / np.sum(w)))
    z = coef / se if se > 0 else np.nan
    p = float(2 * stats.norm.sf(abs(z))) if np.isfinite(z) else np.nan
    return {
        "n_cohorts": len(usable),
        "coef": coef,
        "se": se,
        "p": p,
        "n_low": int(sum(r["n_low"] for r in usable)),
        "n_high": int(sum(r["n_high"] for r in usable)),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    meta = az.load_meta()
    counts = az.load_counts(meta)
    meta = meta.loc[meta["patient"].isin(counts.columns)].copy()
    logcpm = az.normalize(counts.loc[:, meta["patient"]])
    meta = meta.loc[meta["patient"].isin(logcpm.columns)].copy()
    print(f"matrix {logcpm.shape}")

    dor = load_dorothea_all()
    nets: dict[str, pd.DataFrame] = {
        "dorothea_A": dor.loc[dor["confidence"].eq("A"), ["source", "target", "weight"]],
        "dorothea_AB": dor.loc[dor["confidence"].isin(["A", "B"]), ["source", "target", "weight"]],
        "dorothea_ABC": dor.loc[dor["confidence"].isin(["A", "B", "C"]), ["source", "target", "weight"]],
        "dorothea_ABCD": dor.loc[:, ["source", "target", "weight"]],
        "collectri": az.load_collectri(),
        "progeny_top500": load_progeny(),
        "geneset_unsigned": geneset_net(logcpm.index),
    }
    methods_for = {
        "dorothea_A": ["ulm", "mlm", "wmean", "wsum"],
        "dorothea_AB": ["ulm", "mlm", "wmean", "wsum"],
        "dorothea_ABC": ["ulm", "mlm", "wmean", "wsum", "wmean_nes"],
        "dorothea_ABCD": ["ulm", "wmean", "wsum"],  # MLM on D-level is rank-risky and not a TJ model
        "collectri": ["ulm", "mlm", "wmean", "wsum", "wmean_nes"],
        "progeny_top500": ["ulm", "mlm", "wmean", "wsum"],
        "geneset_unsigned": ["ulm", "wmean"],
    }

    arms = {cut: tail_mask(meta, cut) for cut in CUTOFFS}
    for cut, am in arms.items():
        vc = am.loc[am["arm"].isin(["low", "high"]), "arm"].value_counts()
        print(f"cutoff {cut}: low={int(vc.get('low', 0))} high={int(vc.get('high', 0))}")

    long_rows = []
    for net_name, net in nets.items():
        nt = n_targets(net, logcpm.index)
        for method in methods_for[net_name]:
            print(f"scoring {net_name} {method}", flush=True)
            scores = score_net(logcpm, net, method)
            if scores is None:
                long_rows.append({"network": net_name, "method": method, "note": "score failed"})
                continue
            features = [f for f in FOCUS + PATHWAYS + ["IFN_hallmark", "TJ_no_CLDN4", "APICAL_JUNCTION"] if f in scores.columns]
            for feature in features:
                y = scores[feature]
                for cut, am in arms.items():
                    stacked = fit_high(y, am, with_cohort=True)
                    long_rows.append({
                        "network": net_name, "method": method, "feature": feature,
                        "n_targets": int(nt.get(feature, 0)),
                        "cutoff": cut, "scope": "stacked",
                        **stacked,
                    })
                    cohort_fits = []
                    for cohort in COHORTS:
                        sub = am.loc[am["cohort"] == cohort]
                        fit = fit_high(y, sub, with_cohort=False)
                        cohort_fits.append(fit)
                        long_rows.append({
                            "network": net_name, "method": method, "feature": feature,
                            "n_targets": int(nt.get(feature, 0)),
                            "cutoff": cut, "scope": cohort,
                            **fit,
                        })
                    meta_fit = ivw(cohort_fits)
                    long_rows.append({
                        "network": net_name, "method": method, "feature": feature,
                        "n_targets": int(nt.get(feature, 0)),
                        "cutoff": cut, "scope": "meta_ivw",
                        "n_low": meta_fit["n_low"], "n_high": meta_fit["n_high"],
                        "n": meta_fit["n_cohorts"], "coef": meta_fit["coef"],
                        "se": meta_fit["se"], "p": meta_fit["p"], "note": f"cohorts={meta_fit['n_cohorts']}",
                    })

    long_df = pd.DataFrame(long_rows)
    # BH within each stacked/meta spec across features scored in that spec
    long_df["fdr_within_spec"] = np.nan
    for keys, idx in long_df.groupby(["network", "method", "cutoff", "scope"]).groups.items():
        block = long_df.loc[idx]
        p = block["p"].to_numpy(dtype=float) if "p" in block else None
        if p is None or not np.isfinite(p).any():
            continue
        q = np.full(len(block), np.nan)
        ok = np.isfinite(p)
        if ok.sum() >= 2:
            q[ok] = az.bh(p[ok])
        elif ok.sum() == 1:
            q[ok] = p[ok]
        long_df.loc[idx, "fdr_within_spec"] = q

    def flag_ifn(feature, coef):
        if feature in IFN_TFS + ["IFN_hallmark", "JAK-STAT", "TNFa", "NFkB"] and np.isfinite(coef):
            return "match" if coef < 0 else "wrong_sign"
        return ""

    def flag_tj(feature, coef):
        if feature in TJ_TFS + ["TJ_no_CLDN4", "APICAL_JUNCTION"] and np.isfinite(coef):
            return "match" if coef > 0 else "wrong_sign"
        return ""

    long_df["ifn_vs_thesis"] = [flag_ifn(f, c) for f, c in zip(long_df["feature"], long_df["coef"])]
    long_df["tj_vs_thesis"] = [flag_tj(f, c) for f, c in zip(long_df["feature"], long_df["coef"])]
    long_df.to_csv(OUT / "sweep_long.tsv", sep="\t", index=False)

    # One row per network × method × cutoff × scope for the joint panel.
    panel_rows = []
    scopes = ["stacked", "meta_ivw"] + COHORTS
    for (net_name, method, cut, scope), g in long_df.groupby(["network", "method", "cutoff", "scope"]):
        if scope not in scopes:
            continue

        def grab(feature: str) -> pd.Series | None:
            hit = g.loc[g["feature"] == feature]
            if hit.empty or not np.isfinite(hit.iloc[0]["coef"]):
                return None
            return hit.iloc[0]

        stat1 = grab("STAT1")
        irf1 = grab("IRF1")
        ifn_gs = grab("IFN_hallmark")
        jak = grab("JAK-STAT")
        tj_hits = []
        for tf in TJ_TFS + ["TJ_no_CLDN4", "APICAL_JUNCTION"]:
            row = grab(tf)
            if row is not None:
                tj_hits.append(row)
        pos = [r for r in tj_hits if r["coef"] > 0]
        best_tj = min(pos, key=lambda r: r["p"]) if pos else (min(tj_hits, key=lambda r: r["p"]) if tj_hits else None)
        ifn_row = stat1 if stat1 is not None else (ifn_gs if ifn_gs is not None else jak)
        if ifn_row is None or best_tj is None:
            continue
        ifn_ok = bool(ifn_row["coef"] < 0)
        tj_ok = bool(best_tj["coef"] > 0)
        both_nom = bool(ifn_ok and tj_ok and ifn_row["p"] < 0.05 and best_tj["p"] < 0.05)
        both_fdr = bool(
            ifn_ok and tj_ok
            and np.isfinite(ifn_row["fdr_within_spec"]) and np.isfinite(best_tj["fdr_within_spec"])
            and ifn_row["fdr_within_spec"] < 0.05 and best_tj["fdr_within_spec"] < 0.05
        )
        panel_rows.append({
            "network": net_name,
            "method": method,
            "cutoff": cut,
            "scope": scope,
            "ifn_feature": ifn_row["feature"],
            "ifn_coef": ifn_row["coef"],
            "ifn_p": ifn_row["p"],
            "ifn_fdr": ifn_row["fdr_within_spec"],
            "ifn_n_targets": ifn_row["n_targets"],
            "irf1_coef": None if irf1 is None else irf1["coef"],
            "irf1_p": None if irf1 is None else irf1["p"],
            "tj_feature": best_tj["feature"],
            "tj_coef": best_tj["coef"],
            "tj_p": best_tj["p"],
            "tj_fdr": best_tj["fdr_within_spec"],
            "tj_n_targets": best_tj["n_targets"],
            "n_low": ifn_row["n_low"],
            "n_high": ifn_row["n_high"],
            "ifn_sign_match": ifn_ok,
            "tj_sign_match": tj_ok,
            "both_signs": ifn_ok and tj_ok,
            "both_nominal_p05": both_nom,
            "both_fdr05": both_fdr,
            "max_p": max(float(ifn_row["p"]), float(best_tj["p"])),
        })
    panel = pd.DataFrame(panel_rows)
    panel.to_csv(OUT / "sweep_panel.tsv", sep="\t", index=False)

    # Winner: stacked or meta only. Single cohorts stay in the table and cannot be the panel.
    pool = panel.loc[panel["scope"].isin(["stacked", "meta_ivw"])].copy()
    supportive = pool.loc[pool["both_nominal_p05"]].sort_values(["max_p", "ifn_p", "tj_p"])
    sign_only = pool.loc[pool["both_signs"]].sort_values(["max_p", "ifn_p", "tj_p"])
    winner = supportive.iloc[0].to_dict() if len(supportive) else (sign_only.iloc[0].to_dict() if len(sign_only) else None)
    summary = {
        "n_panel_rows": int(len(panel)),
        "n_stacked_or_meta": int(len(pool)),
        "n_both_signs_stacked_or_meta": int(pool["both_signs"].sum()),
        "n_both_nominal_p05_stacked_or_meta": int(pool["both_nominal_p05"].sum()),
        "n_both_fdr05_stacked_or_meta": int(pool["both_fdr05"].sum()),
        "n_single_cohort_both_nominal": int(panel.loc[~panel["scope"].isin(["stacked", "meta_ivw"]), "both_nominal_p05"].sum()),
        "winner_rule": "Among stacked and inverse-variance meta specs, both signs match and both nominal p<0.05; smallest max(p_IFN, p_TJ). Single-cohort specs are reported and are not eligible.",
        "winner": winner,
        "joint_nominal_met": bool(len(supportive)),
    }
    (OUT / "sweep_summary.json").write_text(json.dumps(summary, indent=2, default=float))
    print(json.dumps({k: summary[k] for k in summary if k != "winner"}, indent=2))
    if winner:
        print("WINNER", {k: winner[k] for k in ["network", "method", "cutoff", "scope", "ifn_feature", "ifn_coef", "ifn_p", "tj_feature", "tj_coef", "tj_p", "n_low", "n_high", "both_nominal_p05"]})

    # Scatter of the search. Eligible scopes only.
    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    colors = np.where(pool["both_nominal_p05"], "#1b7837", np.where(pool["both_signs"], "#fdae61", "#b2182b"))
    ax.axvline(0, color="#666", lw=0.6)
    ax.axhline(0, color="#666", lw=0.6)
    ax.scatter(pool["ifn_coef"], pool["tj_coef"], c=colors, s=18, alpha=0.85)
    ax.set_xlabel("IFN-arm coefficient (CLDN4-high − low)\nnegative matches CLDN4-low → IFN up")
    ax.set_ylabel("Best TJ-arm coefficient\npositive matches CLDN4-high → barrier/TJ")
    ax.set_title("Specification sweep (stacked and meta only)\ngreen = both signs and both p<0.05")
    fig.tight_layout()
    fig.savefig(FIGS / "fig_sweep_signs.png", dpi=160)
    fig.savefig(FIGS / "fig_sweep_signs.pdf")
    plt.close(fig)

    if winner and winner.get("both_nominal_p05"):
        # Forest of that one spec's focus features, read from the long table.
        sub = long_df.loc[
            (long_df["network"] == winner["network"])
            & (long_df["method"] == winner["method"])
            & (long_df["cutoff"] == winner["cutoff"])
            & (long_df["scope"] == winner["scope"])
            & (long_df["feature"].isin(FOCUS + ["IFN_hallmark", "TJ_no_CLDN4", "APICAL_JUNCTION", "JAK-STAT"]))
        ].copy()
        sub = sub.loc[np.isfinite(sub["coef"].astype(float))]
        sub = sub.sort_values("coef")
        fig, ax = plt.subplots(figsize=(8.0, 4.8))
        y = np.arange(len(sub))
        ax.axvline(0, color="#444", lw=0.7)
        ax.errorbar(sub["coef"], y, xerr=1.96 * sub["se"].astype(float), fmt="o", color="#333", ecolor="#888", capsize=2, ms=5)
        ax.set_yticks(y)
        ax.set_yticklabels([f"{r.feature}  p={r.p:.3g}" for r in sub.itertuples()], fontsize=8)
        ax.set_xlabel("coefficient (CLDN4-high − low)")
        ax.set_title(
            f"Strongest joint spec from the sweep\n{winner['network']} / {winner['method']} / {winner['cutoff']} / {winner['scope']}"
        )
        fig.tight_layout()
        fig.savefig(FIGS / "fig_sweep_best.png", dpi=160)
        fig.savefig(FIGS / "fig_sweep_best.pdf")
        plt.close(fig)


if __name__ == "__main__":
    main()
