#!/usr/bin/env python3
"""Score sweep: malignant ELF3 and a 4-gene module vs locked T/NK.

The endpoint is fixed: frac_tnk on the locked 65 units, DerSimonian–Laird
Spearman. The search is over malignant score definitions only. p-values on
the searched maximum are descriptive.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
LOCKED_PATH = ROOT / "data" / "locked_units.tsv"
TABLES = ROOT / "results" / "tables"
FIGURES = ROOT / "results" / "figures"
COHORTS = ["GSE123902", "GSE131907", "GSE205335", "GSE189357"]
COLORS = {
    "GSE123902": "#1b4f72",
    "GSE131907": "#b85c38",
    "GSE205335": "#1e7f4f",
    "GSE189357": "#6c3483",
}

UMI_K = [1, 2, 3, 4, 5, 6, 8, 10, 12, 15]
CP_T = [0.25, 0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 10.0, 15.0]
QUANTILES = [0.50, 0.60, 0.70, 0.75, 0.80, 0.90]
CORE = ["ELF3", "TACSTD2", "CLDN4"]
# Fourth gene is swapped inside this list. CLDN7 is the pre-specified member
# from the concordant-4 ELF3 coexpression module. KRT5 is not a candidate.
FOURTHS = [
    "CLDN7", "CDH1", "EPCAM", "KRT8", "KRT18", "KRT19", "GRHL2",
    "CLDN3", "CLDN1", "MUC1", "CDH3", "KLF5", "TJP1", "F11R", "KRT7", "SPDEF", "OCLN", "OVOL2", "GRHL1",
]
PRESPEC_FOURTH = "CLDN7"

# Locked pipeline check (malignant CLDN4 UMI>0 vs frac_tnk).
LOCKED_CLDN4 = {
    "GSE123902": -0.6593406593406593,
    "GSE131907": -0.5220779220779221,
    "GSE205335": -0.43534726143421804,
    "GSE189357": -0.6,
}
LOCKED_DL = -0.5311678045689989


def say(msg: str) -> None:
    print(msg, flush=True)


def spearman_pair(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    n = int(x.size)
    if n < 4 or np.unique(x).size < 2 or np.unique(y).size < 2:
        return np.nan, np.nan, n
    r, p = stats.spearmanr(x, y)
    return float(r), float(p), n


def dl_spearman(rhos, ns):
    rhos = np.asarray(rhos, dtype=float)
    ns = np.asarray(ns, dtype=float)
    ok = np.isfinite(rhos) & np.isfinite(ns) & (ns > 3)
    rhos, ns = rhos[ok], ns[ok]
    if rhos.size == 0:
        return dict(rho=np.nan, p=np.nan, I2=np.nan, ci_lo=np.nan, ci_hi=np.nan, k=0, N=0)
    z = np.arctanh(np.clip(rhos, -0.999999, 0.999999))
    var = 1.0 / (ns - 3.0)
    w = 1.0 / var
    zbar = np.sum(w * z) / np.sum(w)
    q = np.sum(w * (z - zbar) ** 2)
    k = int(rhos.size)
    dfree = k - 1
    cdenom = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    tau2 = max(0.0, (q - dfree) / cdenom) if dfree > 0 and cdenom > 0 else 0.0
    wstar = 1.0 / (var + tau2)
    zre = np.sum(wstar * z) / np.sum(wstar)
    se = float(np.sqrt(1.0 / np.sum(wstar)))
    p = float(2 * stats.norm.sf(abs(zre / se))) if se > 0 else np.nan
    i2 = float(max(0.0, (q - dfree) / q)) if q > 0 else 0.0
    ci = np.tanh(zre + np.array([-1.0, 1.0]) * 1.96 * se)
    return dict(
        rho=float(np.tanh(zre)), p=p, I2=i2, ci_lo=float(ci[0]), ci_hi=float(ci[1]),
        k=k, N=int(ns.sum()),
    )


def within_quartile(x: pd.Series) -> pd.Series:
    r = x.rank(method="average")
    breaks = np.quantile(r.to_numpy(dtype=float), [0, 0.25, 0.5, 0.75, 1.0])
    # identical edges if the score is nearly constant
    if np.unique(breaks).size < 2:
        return pd.Series(["NA"] * len(x), index=x.index)
    breaks = np.maximum.accumulate(breaks)
    for i in range(1, len(breaks)):
        if breaks[i] <= breaks[i - 1]:
            breaks[i] = breaks[i - 1] + 1e-9
    return pd.cut(r, bins=breaks, include_lowest=True, labels=["Q1", "Q2", "Q3", "Q4"]).astype(str)


def rank_biserial_stacked(patient: pd.DataFrame, score: pd.Series) -> dict:
    q_parts = []
    y_parts = []
    for ds, idx in patient.groupby("dataset").groups.items():
        qq = within_quartile(score.loc[idx])
        q_parts.append(qq)
        y_parts.append(patient.loc[idx, "frac_tnk"])
    q = pd.concat(q_parts)
    y = pd.concat(y_parts)
    x4 = y[q.eq("Q4")].to_numpy(dtype=float)
    x1 = y[q.eq("Q1")].to_numpy(dtype=float)
    x4, x1 = x4[np.isfinite(x4)], x1[np.isfinite(x1)]
    if x4.size == 0 or x1.size == 0:
        return dict(r=np.nan, p=np.nan, n_q1=int(x1.size), n_q4=int(x4.size))
    res = stats.mannwhitneyu(x4, x1, alternative="two-sided", method="asymptotic")
    r = 2.0 * float(res.statistic) / (x4.size * x1.size) - 1.0
    return dict(r=r, p=float(res.pvalue), n_q1=int(x1.size), n_q4=int(x4.size))


def load_panels(panel_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    frames = []
    measured_rows = []
    for ds in COHORTS:
        z = np.load(panel_dir / f"{ds}.npz", allow_pickle=True)
        genes = [str(g) for g in z["genes"]]
        unit_id = np.asarray(z["unit_id"]).astype(str)
        counts = z["counts"].astype(np.float64)
        lib = z["lib"].astype(np.float64)
        cp = counts / np.maximum(lib, 1.0)[:, None] * 1e4
        part = pd.DataFrame(counts, columns=[f"umi:{g}" for g in genes])
        part_cp = pd.DataFrame(cp, columns=[f"cp:{g}" for g in genes])
        part.insert(0, "dataset", ds)
        part.insert(1, "unit_id", unit_id)
        part = pd.concat([part, part_cp], axis=1)
        frames.append(part)
        mu = np.asarray(z["measured_units"]).astype(str)
        md = np.asarray(z["measured"]).astype(bool)
        mdf = pd.DataFrame(md, columns=genes)
        mdf.insert(0, "dataset", ds)
        mdf.insert(1, "unit_id", mu)
        measured_rows.append(mdf)
    cells = pd.concat(frames, ignore_index=True)
    measured = pd.concat(measured_rows, ignore_index=True)
    return cells, measured


def gene_usable(measured: pd.DataFrame, gene: str) -> bool:
    if gene not in measured.columns:
        return False
    return bool(measured[gene].all())


def patient_gene_scores(cells: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    """One row per locked unit. Scores that do not use frac_tnk."""
    locked = pd.read_csv(LOCKED_PATH, sep="\t", dtype={"unit_id": str})
    base = locked[["dataset", "unit_id", "n_cells", "n_malignant", "n_tnk", "frac_tnk"]].copy()
    base["unit_id"] = base["unit_id"].astype(str)
    rows = []
    for (ds, unit), sub in cells.groupby(["dataset", "unit_id"], sort=False):
        rec = {"dataset": ds, "unit_id": str(unit)}
        for g in genes:
            umi = sub[f"umi:{g}"].to_numpy(dtype=float)
            cp = sub[f"cp:{g}"].to_numpy(dtype=float)
            for k in UMI_K:
                rec[f"{g}|pct_umi>={k:g}"] = 100.0 * float(np.mean(umi >= k))
            rec[f"{g}|mean_log1p"] = float(np.mean(np.log1p(umi)))
            rec[f"{g}|mean_cp"] = float(np.mean(np.log1p(cp)))
            for t in CP_T:
                rec[f"{g}|pct_cp>={t:g}"] = 100.0 * float(np.mean(cp >= t))
            pos = umi >= 1
            rec[f"{g}|mean_log1p_detected"] = float(np.mean(np.log1p(umi[pos]))) if pos.any() else np.nan
            if umi.size >= 4:
                thr = np.quantile(umi, 0.75)
                tail = umi >= thr
                rec[f"{g}|tail75_mean_log1p"] = float(np.mean(np.log1p(umi[tail]))) if tail.any() else np.nan
            else:
                rec[f"{g}|tail75_mean_log1p"] = np.nan
        rows.append(rec)
    scores = pd.DataFrame(rows)
    # cohort-relative CP10k quantiles (threshold from malignant cells, not from T/NK)
    qthr = {}
    for ds, sub in cells.groupby("dataset", sort=False):
        for g in genes:
            cp = sub[f"cp:{g}"].to_numpy(dtype=float)
            qthr[(ds, g)] = {q: float(np.quantile(cp, q)) for q in QUANTILES}
    qrows = []
    for (ds, unit), sub in cells.groupby(["dataset", "unit_id"], sort=False):
        rec = {"dataset": ds, "unit_id": str(unit)}
        for g in genes:
            cp = sub[f"cp:{g}"].to_numpy(dtype=float)
            for q in QUANTILES:
                rec[f"{g}|pct_cp_q{int(q*100)}"] = 100.0 * float(np.mean(cp >= qthr[(ds, g)][q]))
        qrows.append(rec)
    scores = scores.merge(pd.DataFrame(qrows), on=["dataset", "unit_id"], how="left")
    out = base.merge(scores, on=["dataset", "unit_id"], how="left")
    if out[f"CLDN4|pct_umi>=1"].isna().any():
        raise SystemExit("a locked unit has no malignant counts")
    return out


def meta_of(patient: pd.DataFrame, score: pd.Series) -> dict:
    rhos, ns = [], []
    per = {}
    tmp = patient[["dataset", "frac_tnk"]].copy()
    tmp["_s"] = score.to_numpy()
    for ds in COHORTS:
        g = tmp[tmp.dataset == ds]
        r, p, n = spearman_pair(g["_s"], g["frac_tnk"])
        rhos.append(r)
        ns.append(n)
        per[ds] = r
    dl = dl_spearman(rhos, ns)
    dl["per"] = per
    dl["n_neg"] = int(sum(np.isfinite(v) and v < 0 for v in per.values()))
    dl["abs_rho"] = abs(dl["rho"]) if np.isfinite(dl["rho"]) else np.nan
    return dl


def members_vary(patient: pd.DataFrame, cols: list[str]) -> bool:
    """Every member must vary inside every cohort. A constant member is dropped
    by z-scoring and would make the label '4-gene' false in that cohort."""
    for col in cols:
        for _, idx in patient.groupby("dataset").groups.items():
            x = patient.loc[idx, col].astype(float)
            if int(np.unique(x[np.isfinite(x)]).size) < 2:
                return False
    return True


def z_mean(patient: pd.DataFrame, cols: list[str]) -> pd.Series:
    parts = []
    for col in cols:
        z = pd.Series(np.nan, index=patient.index, dtype=float)
        for _, idx in patient.groupby("dataset").groups.items():
            x = patient.loc[idx, col].astype(float)
            sd = float(x.std(ddof=1))
            if not np.isfinite(sd) or sd == 0 or x.notna().sum() < 4:
                continue
            z.loc[idx] = (x - float(x.mean())) / sd
        parts.append(z)
    return pd.concat(parts, axis=1).mean(axis=1)


def rank_mean(patient: pd.DataFrame, cols: list[str]) -> pd.Series:
    parts = []
    for col in cols:
        r = pd.Series(np.nan, index=patient.index, dtype=float)
        for _, idx in patient.groupby("dataset").groups.items():
            r.loc[idx] = patient.loc[idx, col].rank(method="average")
        parts.append(r)
    return pd.concat(parts, axis=1).mean(axis=1)


def pc1_mean(patient: pd.DataFrame, cols: list[str], align_col: str) -> pd.Series:
    out = pd.Series(np.nan, index=patient.index, dtype=float)
    for _, idx in patient.groupby("dataset").groups.items():
        X = patient.loc[idx, cols].to_numpy(dtype=float)
        if np.isnan(X).any() or X.shape[0] < 4:
            continue
        sd = X.std(axis=0, ddof=1)
        if np.any(sd == 0) or not np.all(np.isfinite(sd)):
            continue
        X = (X - X.mean(axis=0)) / sd
        _u, s, vt = np.linalg.svd(X, full_matrices=False)
        scores = X @ vt[0]
        align = patient.loc[idx, align_col].to_numpy(dtype=float)
        if np.corrcoef(scores, align)[0, 1] < 0:
            scores = -scores
        out.loc[idx] = scores
    return out


def joint_pct(cells: pd.DataFrame, genes: list[str], kind: str, cut: float, how: str) -> pd.Series:
    """how = 'all' or 'any3'. kind = 'umi' or 'cp'."""
    prefix = "umi" if kind == "umi" else "cp"
    recs = []
    for (ds, unit), sub in cells.groupby(["dataset", "unit_id"], sort=False):
        mat = np.column_stack([sub[f"{prefix}:{g}"].to_numpy(dtype=float) >= cut for g in genes])
        if how == "all":
            hit = mat.all(axis=1)
        elif how == "any3":
            hit = mat.sum(axis=1) >= min(3, mat.shape[1])
        else:
            raise ValueError(how)
        recs.append({"dataset": ds, "unit_id": str(unit), "score": 100.0 * float(hit.mean())})
    return pd.DataFrame(recs)


def attach(patient: pd.DataFrame, spec: pd.DataFrame) -> pd.Series:
    key = patient[["dataset", "unit_id"]].merge(spec, on=["dataset", "unit_id"], how="left")
    if key["score"].isna().any():
        raise SystemExit("joint score missed a unit")
    return pd.Series(key["score"].to_numpy(), index=patient.index)


def row_from_meta(family, name, genes, score_kind, detail, meta) -> dict:
    rec = {
        "family": family,
        "name": name,
        "genes": ",".join(genes),
        "fourth": genes[3] if len(genes) == 4 else "",
        "prespec_module": int(genes == CORE + [PRESPEC_FOURTH]),
        "score_kind": score_kind,
        "detail": detail,
        "rho": meta["rho"],
        "p": meta["p"],
        "I2": meta["I2"],
        "ci_lo": meta["ci_lo"],
        "ci_hi": meta["ci_hi"],
        "k": meta["k"],
        "N": meta["N"],
        "n_neg": meta["n_neg"],
        "abs_rho": meta["abs_rho"],
    }
    for ds in COHORTS:
        rec[f"rho_{ds}"] = meta["per"][ds]
    return rec


def main() -> None:
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--panel", type=Path, default=Path("/tmp/geo_c4/panel"))
    args = p.parse_args()
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)

    say("load panels")
    cells, measured = load_panels(args.panel)
    genes = [g for g in ["ELF3", "TACSTD2", "CLDN4", *FOURTHS, "KRT5", "GRHL3"] if g in measured.columns]
    usable = {g: gene_usable(measured, g) for g in genes}
    say("usable " + ", ".join(f"{g}={int(v)}" for g, v in usable.items()))
    missing_primary = [g for g in CORE + [PRESPEC_FOURTH] if not usable.get(g, False)]
    if missing_primary:
        raise SystemExit(f"primary genes not measured in every unit: {missing_primary}")

    say("patient gene scores")
    patient = patient_gene_scores(cells, [g for g in genes if usable[g]])
    # pipeline check
    cldn = patient["CLDN4|pct_umi>=1"]
    meta = meta_of(patient, cldn)
    for ds, exp in LOCKED_CLDN4.items():
        got = meta["per"][ds]
        if abs(got - exp) > 1e-6:
            raise SystemExit(f"CLDN4 cohort rho mismatch {ds}: {got} != {exp}")
    if abs(meta["rho"] - LOCKED_DL) > 1e-6 or meta["N"] != 65:
        raise SystemExit(f"CLDN4 DL mismatch {meta}")
    locked_pct = pd.read_csv(LOCKED_PATH, sep="\t", dtype={"unit_id": str})
    cmp = patient.merge(locked_pct[["dataset", "unit_id", "mal_CLDN4_pct"]], on=["dataset", "unit_id"])
    delta = np.nanmax(np.abs(cmp["CLDN4|pct_umi>=1"] - cmp["mal_CLDN4_pct"]))
    if delta > 1e-4:
        raise SystemExit(f"CLDN4 %pos differs from locked table by {delta}")
    say(f"pipeline check ok  CLDN4 %pos DL rho={meta['rho']:.6f}  max|pct-locked|={delta:.3g}")

    records = []
    cache = {}

    def add(family, name, gene_list, kind, detail, series):
        meta = meta_of(patient, series)
        if meta["N"] != 65 or meta["k"] != 4:
            return
        records.append(row_from_meta(family, name, gene_list, kind, detail, meta))
        cache[name] = series

    # ELF3 marginal
    g = "ELF3"
    for k in UMI_K:
        add("ELF3", f"ELF3 pct UMI>={k:g}", [g], "pct_umi", f"k={k:g}", patient[f"{g}|pct_umi>={k:g}"])
    for t in CP_T:
        add("ELF3", f"ELF3 pct CP10k>={t:g}", [g], "pct_cp", f"t={t:g}", patient[f"{g}|pct_cp>={t:g}"])
    for q in QUANTILES:
        add("ELF3", f"ELF3 pct above cohort CP q{int(q*100)}", [g], "pct_cohort_q", f"q={q:.2f}", patient[f"{g}|pct_cp_q{int(q*100)}"])
    add("ELF3", "ELF3 mean log1p UMI", [g], "mean_log1p", "mean", patient[f"{g}|mean_log1p"])
    add("ELF3", "ELF3 mean log1p CP10k", [g], "mean_cp", "mean", patient[f"{g}|mean_cp"])
    add("ELF3", "ELF3 mean log1p among detected", [g], "mean_detected", "detected", patient[f"{g}|mean_log1p_detected"])
    add("ELF3", "ELF3 tail75 mean log1p", [g], "tail75", "within-unit q75", patient[f"{g}|tail75_mean_log1p"])

    memberships = []
    for fourth in FOURTHS:
        if not usable.get(fourth, False):
            say(f"skip fourth {fourth}: not measured in every unit")
            continue
        memberships.append(CORE + [fourth])

    say(f"module memberships {len(memberships)}")
    for genes4 in memberships:
        tag = "+".join(genes4)
        fourth = genes4[3]
        # patient-level aggregates
        for k in UMI_K:
            cols = [f"{g}|pct_umi>={k:g}" for g in genes4]
            if members_vary(patient, cols):
                add("module", f"zmean pct UMI>={k:g} {tag}", genes4, "zmean_pct_umi", f"k={k:g}", z_mean(patient, cols))
                add("module", f"rmean pct UMI>={k:g} {tag}", genes4, "rmean_pct_umi", f"k={k:g}", rank_mean(patient, cols))
        for t in CP_T:
            cols = [f"{g}|pct_cp>={t:g}" for g in genes4]
            if members_vary(patient, cols):
                add("module", f"zmean pct CP10k>={t:g} {tag}", genes4, "zmean_pct_cp", f"t={t:g}", z_mean(patient, cols))
        for q in QUANTILES:
            cols = [f"{g}|pct_cp_q{int(q*100)}" for g in genes4]
            if members_vary(patient, cols):
                add("module", f"zmean pct cohort-q{int(q*100)} {tag}", genes4, "zmean_pct_q", f"q={q:.2f}", z_mean(patient, cols))
        cols_mean = [f"{g}|mean_log1p" for g in genes4]
        cols_cp = [f"{g}|mean_cp" for g in genes4]
        if members_vary(patient, cols_mean):
            add("module", f"zmean mean-log1p {tag}", genes4, "zmean_mean_log1p", "mean", z_mean(patient, cols_mean))
            add("module", f"rmean mean-log1p {tag}", genes4, "rmean_mean_log1p", "mean", rank_mean(patient, cols_mean))
            add("module", f"PC1 mean-log1p {tag}", genes4, "pc1_mean_log1p", "PC1", pc1_mean(patient, cols_mean, "ELF3|mean_log1p"))
        if members_vary(patient, cols_cp):
            add("module", f"zmean mean-cp {tag}", genes4, "zmean_mean_cp", "mean", z_mean(patient, cols_cp))
        pct1_cols = [f"{g}|pct_umi>=1" for g in genes4]
        if members_vary(patient, pct1_cols):
            add("module", f"PC1 pct UMI>=1 {tag}", genes4, "pc1_pct_umi", "k=1", pc1_mean(patient, pct1_cols, "ELF3|pct_umi>=1"))

        # cell-level co-detection. Same cut on every member.
        for k in UMI_K:
            add("module", f"joint UMI>={k:g} {tag}", genes4, "joint_umi", f"k={k:g}", attach(patient, joint_pct(cells, genes4, "umi", k, "all")))
            add("module", f"any3 UMI>={k:g} {tag}", genes4, "any3_umi", f"k={k:g}", attach(patient, joint_pct(cells, genes4, "umi", k, "any3")))
        for t in CP_T:
            add("module", f"joint CP10k>={t:g} {tag}", genes4, "joint_cp", f"t={t:g}", attach(patient, joint_pct(cells, genes4, "cp", t, "all")))

        say(f"  finished {fourth}")

    grid = pd.DataFrame(records)
    grid.to_csv(TABLES / "sweep_grid.tsv", sep="\t", index=False)
    say(f"grid rows {len(grid)}")

    def pick(frame: pd.DataFrame, concordant: bool) -> pd.Series:
        d = frame[np.isfinite(frame.abs_rho)].copy()
        if concordant:
            d = d[d.n_neg == 4]
        if d.empty:
            raise SystemExit("no eligible row")
        d = d.sort_values(["abs_rho", "I2", "p"], ascending=[False, True, True])
        return d.iloc[0]

    elf = grid[grid.family == "ELF3"]
    prespec = grid[(grid.family == "module") & (grid.prespec_module == 1)]
    modules = grid[grid.family == "module"]
    winners = {
        "ELF3_concordant": pick(elf, True),
        "ELF3_unconstrained": pick(elf, False),
        "module_prespec_concordant": pick(prespec, True),
        "module_prespec_unconstrained": pick(prespec, False),
        "module_any_concordant": pick(modules, True),
        "module_any_unconstrained": pick(modules, False),
    }
    # baselines
    base_names = {
        "ELF3_pct_umi_ge1": "ELF3 pct UMI>=1",
        "module_prespec_zmean_pct_umi_ge1": f"zmean pct UMI>=1 {'+'.join(CORE + [PRESPEC_FOURTH])}",
    }
    summary_rows = []
    for key, rec in winners.items():
        summary_rows.append({"slot": key, **rec.to_dict()})
    for key, name in base_names.items():
        hit = grid[grid.name == name]
        if hit.empty:
            raise SystemExit(f"missing baseline {name}")
        summary_rows.append({"slot": key, **hit.iloc[0].to_dict()})
    # CLDN4 check row
    summary_rows.append({
        "slot": "CLDN4_pipeline_check",
        "family": "check",
        "name": "CLDN4 pct UMI>=1",
        "genes": "CLDN4",
        "fourth": "",
        "prespec_module": 0,
        "score_kind": "pct_umi",
        "detail": "k=1",
        "rho": meta["rho"],
        "p": meta["p"],
        "I2": meta["I2"],
        "ci_lo": meta["ci_lo"],
        "ci_hi": meta["ci_hi"],
        "k": meta["k"],
        "N": meta["N"],
        "n_neg": meta["n_neg"],
        "abs_rho": abs(meta["rho"]),
        **{f"rho_{ds}": meta["per"][ds] for ds in COHORTS},
    })
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(TABLES / "winners.tsv", sep="\t", index=False)

    # patient vectors for headline scores + baselines + check
    keep_names = [winners["ELF3_concordant"]["name"], winners["module_prespec_concordant"]["name"], winners["module_any_concordant"]["name"]]
    keep_names += list(base_names.values())
    # unconstrained too if different
    for slot in winners:
        keep_names.append(winners[slot]["name"])
    vectors = patient[["dataset", "unit_id", "n_malignant", "frac_tnk"]].copy()
    vectors["CLDN4_pct_umi_ge1"] = patient["CLDN4|pct_umi>=1"]
    for name in dict.fromkeys(keep_names):
        if name not in cache:
            continue
        vectors[name] = cache[name].to_numpy()
    vectors.to_csv(TABLES / "patient_winners.tsv", sep="\t", index=False)

    qrows = []
    for name in dict.fromkeys(keep_names + ["CLDN4 pct UMI>=1"]):
        if name == "CLDN4 pct UMI>=1":
            series = patient["CLDN4|pct_umi>=1"]
        elif name not in cache:
            continue
        else:
            series = cache[name]
        rb = rank_biserial_stacked(patient, series)
        qrows.append({"name": name, **rb})
    pd.DataFrame(qrows).to_csv(TABLES / "winners_q4q1.tsv", sep="\t", index=False)

    # where the pre-specified concordant winner sits in the full module list
    ranked = modules[modules.n_neg == 4].sort_values(["abs_rho", "I2"], ascending=[False, True]).reset_index(drop=True)
    ranked.head(15).to_csv(TABLES / "top_concordant_modules.tsv", sep="\t", index=False)
    pres_name = winners["module_prespec_concordant"]["name"]
    rank_pos = int(ranked.index[ranked.name == pres_name][0]) + 1
    meta_out = {
        "n_specs": int(len(grid)),
        "n_elf3": int(len(elf)),
        "n_module": int(len(modules)),
        "n_prespec": int(len(prespec)),
        "n_module_concordant": int((modules.n_neg == 4).sum()),
        "prespec_concordant_rank_among_concordant_modules": rank_pos,
        "n_fourth_genes_tested": int(modules.fourth.nunique()),
        "fourths_tested": sorted(modules.fourth.unique()),
    }
    (TABLES / "sweep_meta.json").write_text(json.dumps(meta_out, indent=2) + "\n")
    say(json.dumps(meta_out, indent=2))
    for slot, rec in winners.items():
        say(
            f"{slot}: rho={rec['rho']:+.4f} p={rec['p']:.3g} I2={rec['I2']:.3f} "
            f"n_neg={int(rec['n_neg'])} {rec['name']}"
        )

    plot(patient, grid, winners, cache)
    say("done")


def plot(patient, grid, winners, cache) -> None:
    # 1. ELF3 threshold curves
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.8), sharey=True)
    elf = grid[grid.family == "ELF3"]
    for ax, kind, xlabel in (
        (axes[0], "pct_umi", "malignant cells with ELF3 UMI ≥ k (%)"),
        (axes[1], "pct_cp", "malignant cells with ELF3 CP10k ≥ t (%)"),
    ):
        sub = elf[elf.score_kind == kind].copy()
        # parse numeric from detail
        sub["x"] = sub["detail"].str.split("=").str[1].astype(float)
        sub = sub.sort_values("x")
        ax.plot(sub["x"], sub["rho"], color="#1b4f72", marker="o", ms=4)
        ax.axhline(LOCKED_DL, color="#888888", lw=0.8, ls="--")
        ax.axhline(0, color="#cccccc", lw=0.6)
        ax.set_xlabel(xlabel.split("with ")[-1].replace(" (%)", ""))
        ax.set_title(kind)
    axes[0].set_ylabel("DL Spearman vs T/NK fraction")
    fig.suptitle("ELF3 score vs threshold (dashed = locked CLDN4 %pos)", fontsize=11)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig_elf3_threshold.png", dpi=160)
    fig.savefig(FIGURES / "fig_elf3_threshold.pdf")
    plt.close()

    # 2. prespec module vs k
    fig, ax = plt.subplots(figsize=(6.6, 3.9))
    pres = grid[(grid.prespec_module == 1) & (grid.score_kind.isin(["zmean_pct_umi", "rmean_pct_umi", "joint_umi", "any3_umi"]))].copy()
    pres["x"] = pres["detail"].str.split("=").str[1].astype(float)
    for kind, color, label in (
        ("rmean_pct_umi", "#6c3483", "mean of ranks of %pos"),
        ("zmean_pct_umi", "#1b4f72", "mean of z(%pos)"),
        ("joint_umi", "#b85c38", "all four genes ≥ k"),
        ("any3_umi", "#1e7f4f", "at least three ≥ k"),
    ):
        sub = pres[pres.score_kind == kind].sort_values("x")
        if sub.empty:
            continue
        ax.plot(sub["x"], sub["rho"], marker="o", ms=4, color=color, label=label)
    ax.axhline(LOCKED_DL, color="#888888", lw=0.8, ls="--", label="CLDN4 %pos")
    ax.axhline(0, color="#cccccc", lw=0.6)
    ax.set_xlabel("UMI threshold k")
    ax.set_ylabel("DL Spearman vs T/NK fraction")
    ax.set_title("ELF3 + TACSTD2 + CLDN4 + CLDN7")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig_module_threshold.png", dpi=160)
    fig.savefig(FIGURES / "fig_module_threshold.pdf")
    plt.close()

    # 3. forest of headline + baselines
    prespec_w = winners["module_prespec_concordant"]
    search_w = winners["module_any_concordant"]
    slots = [
        ("CLDN4 %pos, pipeline check", "CLDN4 pct UMI>=1"),
        ("ELF3 %pos, max |ρ|", winners["ELF3_concordant"]["name"]),
        ("CLDN7 module, z(%pos)", f"zmean pct UMI>=1 {'+'.join(CORE + [PRESPEC_FOURTH])}"),
        (f"{prespec_w['fourth']} module, concordant max", prespec_w["name"]),
        (f"{search_w['fourth']} module, search max", search_w["name"]),
    ]
    # build rho table from grid / check
    fig, ax = plt.subplots(figsize=(8.2, 4.4))
    ypos = np.arange(len(slots))[::-1]
    for y, (label, name) in zip(ypos, slots):
        if name == "CLDN4 pct UMI>=1":
            rhos = [LOCKED_CLDN4[ds] for ds in COHORTS]
            # DL from grid isn't stored; use locked
            dl = LOCKED_DL
        else:
            hit = grid[grid.name == name].iloc[0]
            rhos = [hit[f"rho_{ds}"] for ds in COHORTS]
            dl = hit["rho"]
        ax.scatter(rhos, [y] * 4, c=[COLORS[ds] for ds in COHORTS], s=28, zorder=3)
        ax.plot([dl], [y], marker="D", color="black", ms=5, zorder=4)
    ax.axvline(0, color="#cccccc", lw=0.8)
    ax.set_yticks(ypos)
    ax.set_yticklabels([s[0] for s in slots])
    ax.set_xlabel("Spearman vs T/NK fraction")
    ax.set_title("Cohort ρ (dots) and DL pool (diamond)")
    handles = [plt.Line2D([0], [0], marker="o", color=c, ls="", label=ds) for ds, c in COLORS.items()]
    handles.append(plt.Line2D([0], [0], marker="D", color="black", ls="", label="DL"))
    ax.legend(handles=handles, frameon=False, fontsize=8, loc="lower left")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig_forest.png", dpi=160)
    fig.savefig(FIGURES / "fig_forest.pdf")
    plt.close()

    # 4. scatter of the two headlines
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.2), sharey=True)
    panels = (
        (axes[0], winners["ELF3_concordant"], "ELF3 % of malignant cells, UMI ≥ 1"),
        (axes[1], winners["module_any_concordant"], "Rank mean of % positive\nELF3, TACSTD2, CLDN4, " + winners["module_any_concordant"]["fourth"]),
    )
    for ax, rec, xlabel in panels:
        series = cache[rec["name"]]
        for ds, idx in patient.groupby("dataset").groups.items():
            ax.scatter(series.loc[idx], patient.loc[idx, "frac_tnk"], s=22, color=COLORS[ds], label=ds)
        ax.set_xlabel(xlabel)
        ax.set_title(f"DL ρ = {rec['rho']:+.3f}")
    axes[0].set_ylabel("T/NK fraction")
    axes[1].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig_scatter.png", dpi=160)
    fig.savefig(FIGURES / "fig_scatter.pdf")
    plt.close()


if __name__ == "__main__":
    main()
