#!/usr/bin/env python3
"""Additive B4 extras: public lung ICI cohorts (not only GSE126044) + TCGA TJ vs CD8/GEP.

Paper-facing figures/tables. Honest n / medians / ρ / p. Public data only.
"""

from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.duration.hazard_regression import PHReg

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from gene_sets import (  # noqa: E402
    ALIASES,
    CD8,
    CLDN147_F11R_PARD3,
    CLDN4_ALONE,
    USER_PPT_7GENE,
)

ROOT = HERE.parents[2]
DATA = ROOT / "data" / "rework" / "B4_wave2"
EXTRA_DATA = DATA / "extra"
OUT = ROOT / "results" / "rework" / "B4_wave2" / "extra"
TAB = OUT / "tables"
FIG = OUT / "figures"
for d in (TAB, FIG):
    d.mkdir(parents=True, exist_ok=True)

# Ayers et al. 2017 GEP (18 genes)
GEP = [
    "CCL5", "CD27", "CD274", "CD276", "CD8A", "CMKLR1", "CXCL9", "CXCR6",
    "HLA-DQA1", "HLA-DRB1", "HLA-E", "IDO1", "LAG3", "NKG7", "PDCD1LG2",
    "PSMB10", "STAT1", "TIGIT",
]
CYT = ["GZMA", "PRF1"]
OCLN_ALONE = ["OCLN"]

# Frozen Ensembl (strip version at match time). Used only for GSE135222.
ENSEMBL = {
    "CLDN1": "ENSG00000163347",
    "CLDN4": "ENSG00000189143",
    "CLDN7": "ENSG00000181885",
    "F11R": "ENSG00000158769",
    "PARD3": "ENSG00000148498",
    "TJP1": "ENSG00000104067",
    "TJP2": "ENSG00000119139",
    "OCLN": "ENSG00000197826",
    "CD8A": "ENSG00000153563",
    "CD8B": "ENSG00000172116",
    "GZMA": "ENSG00000145649",
    "PRF1": "ENSG00000180644",
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
# GRCh37 CCL5 often ENSG00000161570
ENSEMBL_ALT = {"CCL5": ["ENSG00000161570"]}


def log2p1(df: pd.DataFrame) -> pd.DataFrame:
    return np.log2(df.clip(lower=0) + 1)


def resolve(genes: list[str], index: pd.Index) -> list[str]:
    idx = pd.Index(index.astype(str))
    # Ensembl matrix?
    stripped = pd.Index([x.split(".")[0] for x in idx])
    is_ens = stripped.str.startswith("ENSG").mean() > 0.5
    matched = []
    seen = set()
    for g in genes:
        if is_ens:
            ids = [ENSEMBL[g]] + ENSEMBL_ALT.get(g, []) if g in ENSEMBL else []
            hit = None
            for eid in ids:
                if eid in set(stripped):
                    hit = idx[stripped == eid][0]
                    break
            if hit and hit not in seen:
                matched.append(hit)
                seen.add(hit)
            continue
        cands = [g] + ALIASES.get(g, [])
        hit = next((c for c in cands if c in idx and c not in seen), None)
        if hit:
            matched.append(hit)
            seen.add(hit)
    return matched


def meanz(expr: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, list[str]]:
    used = resolve(genes, expr.index)
    if not used:
        return pd.Series(np.nan, index=expr.columns), []
    sub = expr.loc[used]
    var = sub.var(axis=1)
    used = var[var > 0].index.tolist()
    sub = sub.loc[used]
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1, ddof=1), axis=0)
    return z.mean(axis=0), used


def mwu_high_in_b(a: pd.Series, b: pd.Series) -> dict:
    """Compare group B vs group A. Report whether B has higher scores."""
    a = a.dropna().astype(float)
    b = b.dropna().astype(float)
    if len(a) < 2 or len(b) < 2:
        return {"n_A": len(a), "n_B": len(b), "p_two": np.nan}
    _, p = stats.mannwhitneyu(b.values, a.values, alternative="two-sided")
    return {
        "n_A": int(len(a)),
        "n_B": int(len(b)),
        "median_A": float(a.median()),
        "median_B": float(b.median()),
        "p_two": float(p),
        "direction": "B>A" if b.median() > a.median() else ("A>B" if a.median() > b.median() else "tie"),
    }


def spearman(x: pd.Series, y: pd.Series) -> dict:
    d = pd.concat([x.rename("x"), y.rename("y")], axis=1).dropna()
    if len(d) < 5:
        return {"n": int(len(d)), "rho": np.nan, "p": np.nan}
    rho, p = stats.spearmanr(d["x"], d["y"])
    return {"n": int(len(d)), "rho": float(rho), "p": float(p)}


def cox_per_sd(time: pd.Series, event: pd.Series, score: pd.Series) -> dict:
    d = pd.concat([time.rename("t"), event.rename("e"), score.rename("x")], axis=1).dropna()
    d["t"] = pd.to_numeric(d["t"], errors="coerce")
    d["e"] = pd.to_numeric(d["e"], errors="coerce")
    d = d.dropna()
    d = d[d["t"] > 0]
    if len(d) < 8 or d["e"].sum() < 3 or d["x"].std(ddof=1) == 0:
        return {"n": int(len(d)), "n_event": int(d["e"].sum()) if len(d) else 0, "hr_per_sd": np.nan, "p": np.nan}
    x = ((d["x"] - d["x"].mean()) / d["x"].std(ddof=1)).values
    res = PHReg(d["t"].values, x, status=d["e"].values).fit()
    return {
        "n": int(len(d)),
        "n_event": int(d["e"].sum()),
        "hr_per_sd": float(np.exp(res.params[0])),
        "loghr": float(res.params[0]),
        "se": float(res.bse[0]),
        "p": float(res.pvalues[0]),
    }


def logrank_median(time: pd.Series, event: pd.Series, score: pd.Series) -> dict:
    d = pd.concat([time.rename("t"), event.rename("e"), score.rename("x")], axis=1).dropna()
    d["t"] = pd.to_numeric(d["t"], errors="coerce")
    d["e"] = pd.to_numeric(d["e"], errors="coerce")
    d = d.dropna()
    d = d[d["t"] > 0]
    if len(d) < 8 or d["e"].sum() < 3:
        return {"n": int(len(d)), "p": np.nan}
    med = float(d["x"].median())
    d["hi"] = (d["x"] > med).astype(int)
    # Mantel–Haenszel log-rank
    times = np.sort(d.loc[d["e"] == 1, "t"].unique())
    o1 = e1 = v = 0.0
    for t in times:
        at = d[d["t"] >= t]
        n = len(at)
        n1 = int((at["hi"] == 1).sum())
        dth = d[(d["t"] == t) & (d["e"] == 1)]
        di = len(dth)
        d1 = int((dth["hi"] == 1).sum())
        if n <= 1 or n1 == 0 or n1 == n:
            continue
        ei = di * n1 / n
        vi = (di * (n - di) * n1 * (n - n1)) / (n ** 2 * (n - 1)) if n > 1 else 0
        o1 += d1
        e1 += ei
        v += vi
    if v <= 0:
        return {"n": int(len(d)), "n_high": int(d["hi"].sum()), "median": med, "p": np.nan}
    z = (o1 - e1) / np.sqrt(v)
    p = float(2 * stats.norm.sf(abs(z)))
    return {
        "n": int(len(d)),
        "n_high": int(d["hi"].sum()),
        "n_low": int((d["hi"] == 0).sum()),
        "median": med,
        "p": p,
        "z": float(z),
    }


def km_curve(time: np.ndarray, event: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    order = np.argsort(time)
    t = time[order]
    e = event[order]
    n = len(t)
    surv = 1.0
    xs = [0.0]
    ys = [1.0]
    i = 0
    while i < n:
        ti = t[i]
        at_risk = n - i
        died = 0
        while i < n and t[i] == ti:
            died += int(e[i] == 1)
            i += 1
        if died and at_risk:
            surv *= 1 - died / at_risk
        xs.extend([ti, ti])
        ys.extend([ys[-1], surv])
    return np.array(xs), np.array(ys)


def km_plot(time: pd.Series, event: pd.Series, score: pd.Series, title: str, fname: Path) -> None:
    d = pd.concat([time.rename("t"), event.rename("e"), score.rename("x")], axis=1).dropna()
    d["t"] = pd.to_numeric(d["t"], errors="coerce")
    d["e"] = pd.to_numeric(d["e"], errors="coerce")
    d = d.dropna()
    d = d[d["t"] > 0]
    if len(d) < 8:
        return
    med = float(d["x"].median())
    fig, ax = plt.subplots(figsize=(4.6, 4.0))
    for lab, mask, c in [
        (f"≤ median (n={(d['x'] <= med).sum()})", d["x"] <= med, "#2c7fb8"),
        (f"> median (n={(d['x'] > med).sum()})", d["x"] > med, "#d95f0e"),
    ]:
        xs, ys = km_curve(d.loc[mask, "t"].values, d.loc[mask, "e"].values)
        ax.step(xs, ys, where="post", color=c, label=lab, lw=1.6)
    lr = logrank_median(d["t"], d["e"], d["x"])
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("time")
    ax.set_ylabel("PFS")
    ax.set_title(f"{title}\nlog-rank p={lr.get('p', np.nan):.3g}", fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(fname.with_suffix(".png"), dpi=150)
    fig.savefig(fname.with_suffix(".pdf"))
    plt.close(fig)


def parse_chars(path: Path) -> pd.DataFrame:
    titles = gsms = None
    rows: dict[str, list[str]] = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("!Sample_title"):
                titles = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_geo_accession"):
                gsms = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                vals = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
                keys = []
                parsed = []
                for v in vals:
                    if ": " in v:
                        k, rest = v.split(": ", 1)
                    else:
                        k, rest = "value", v
                    keys.append(k)
                    parsed.append(rest)
                key = max(set(keys), key=keys.count)
                rows[key] = parsed
    meta = pd.DataFrame(rows)
    meta["title"] = titles
    meta["gsm"] = gsms
    return meta


def load_gse135222() -> tuple[pd.DataFrame, pd.DataFrame]:
    expr = pd.read_csv(EXTRA_DATA / "GSE135222_exp.tsv.gz", sep="\t", index_col=0)
    expr = log2p1(expr)
    meta = parse_chars(EXTRA_DATA / "GSE135222_series_matrix.txt.gz")
    meta["sample"] = meta["title"].str.replace(" ", "", regex=False)
    meta = meta.set_index("sample")
    meta["pfs_days"] = pd.to_numeric(meta["pfs.time"], errors="coerce")
    meta["pfs_event"] = pd.to_numeric(meta["progression-free survival (pfs)"], errors="coerce")
    # Public GEO DCB: PFS ≥ 183 days (6 months). All censored patients have long follow-up.
    meta["group"] = np.where(meta["pfs_days"] >= 183, "DCB", "NDB")
    meta["group_high"] = "NDB"  # for MW: is NDB higher? (same direction as B4 NR-higher TJ)
    common = [c for c in expr.columns if c in meta.index]
    return expr[common], meta.loc[common]


def load_gse207422() -> tuple[pd.DataFrame, pd.DataFrame]:
    expr = pd.read_csv(EXTRA_DATA / "GSE207422_log2TPM.txt.gz", sep="\t", index_col=0)
    md = pd.read_excel(EXTRA_DATA / "GSE207422_metadata.xlsx")
    md = md[md["Sample"].isin(expr.columns)].copy()
    md["group"] = md["Pathologic Response"].map(
        lambda x: "MPR" if isinstance(x, str) and str(x).startswith("MPR") else (
            "NMPR" if x == "NMPR" else np.nan
        )
    )
    md = md.dropna(subset=["group"]).set_index("Sample")
    return expr[md.index], md


def load_gse190265() -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = pd.read_csv(EXTRA_DATA / "GSE190265_TPM.csv.gz", sep=None, engine="python", index_col=0)
    expr = log2p1(raw.T)  # genes x samples
    samp = pd.read_csv(EXTRA_DATA / "GSE190265_samples.csv.gz", sep=";").set_index("sample")
    samp["pfs_months"] = pd.to_numeric(samp["time_PFS"], errors="coerce")
    samp["pfs_event"] = pd.to_numeric(samp["evtPFS"], errors="coerce")
    samp["group"] = np.where(samp["pfs_months"] >= 6, "DCB", "NDB")
    common = [c for c in expr.columns if c in samp.index]
    return expr[common], samp.loc[common]


def load_gse166449() -> tuple[pd.DataFrame, pd.DataFrame]:
    expr = pd.read_csv(EXTRA_DATA / "GSE166449_TPM.txt.gz", sep="\t", index_col=0)
    expr = log2p1(expr)
    meta = parse_chars(EXTRA_DATA / "GSE166449_series_matrix.txt.gz")
    # Supplementary columns are in GEO sample order (C001 = first GSM).
    cols = list(expr.columns)
    if len(cols) != len(meta):
        raise RuntimeError(f"GSE166449 column/meta mismatch {len(cols)} vs {len(meta)}")
    meta = meta.copy()
    meta["sample"] = cols
    meta["group"] = np.where(meta["title"].str.contains("nonResponder", case=False), "NR", "R")
    # Confirm from title: Responder vs nonResponder
    meta.loc[meta["title"].str.contains("Immunotherapy_Responder") & ~meta["title"].str.contains("non"), "group"] = "R"
    meta.loc[meta["title"].str.contains("nonResponder"), "group"] = "NR"
    meta = meta.set_index("sample")
    return expr[meta.index], meta


def load_gse190266() -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = pd.read_csv(EXTRA_DATA / "GSE190266_TPM.csv.gz", sep=";", decimal=",", index_col=0)
    expr = log2p1(raw.T)  # genes x samples
    meta = parse_chars(EXTRA_DATA / "GSE190266_series_matrix.txt.gz")
    meta = meta.set_index("title")
    meta["pfs_months"] = pd.to_numeric(meta["pfs_time (6 months)"], errors="coerce")
    meta["pfs_event"] = pd.to_numeric(meta["pfs_evt (6 months)"], errors="coerce")
    # GEO field is PFS truncated at 6 months; DCB = still progression-free at 6 mo.
    meta["group"] = np.where(meta["pfs_months"] >= 6, "DCB", "NDB")
    common = [c for c in expr.columns if c in meta.index]
    return expr[common], meta.loc[common]


def load_gse161537() -> tuple[pd.DataFrame, pd.DataFrame]:
    # Targeted ~2.5k-gene log2CPM (HTG / similar). No CLDN4 / CLDN1 / CLDN7 / TJP.
    expr = pd.read_csv(EXTRA_DATA / "GSE161537_log2cpm.csv.gz", sep=";", decimal=",", index_col=0)
    meta = parse_chars(EXTRA_DATA / "GSE161537_series_matrix.txt.gz")
    meta["sample"] = meta["patient id"].astype(str)
    meta = meta.set_index("sample")
    recist = meta["best response on immunotherapy (recist)"].astype(str)
    meta["recist"] = recist
    meta["group"] = recist.map({"CR": "R", "PR": "R", "PD": "NR"})
    meta["pfs_months"] = pd.to_numeric(meta["pfs (month)"], errors="coerce")
    expr.columns = expr.columns.astype(str)
    common = [c for c in expr.columns if c in meta.index]
    return expr[common], meta.loc[common]


def load_gse126044() -> tuple[pd.DataFrame, pd.DataFrame]:
    counts = pd.read_csv(DATA / "GSE126044_counts.txt.gz", sep="\t", index_col=0)
    lib = counts.sum(axis=0)
    expr = np.log2(counts.divide(lib, axis=1) * 1e6 + 1)
    meta = parse_chars(DATA / "GSE126044_series_matrix.txt.gz")
    meta["sample"] = meta["title"].str.replace("RNA-seq_", "", regex=False)
    meta = meta.set_index("sample")
    meta["group"] = meta["patient response"].map({"responder": "R", "non-responder": "NR"})
    return expr[meta.index], meta


def load_tcga(path: Path) -> pd.DataFrame:
    expr = pd.read_csv(path, sep="\t", index_col=0)
    tumors = [c for c in expr.columns if str(c).split("-")[-1].startswith("01")]
    return expr[tumors]


def score_features(expr: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    feats = {
        "CLDN4": CLDN4_ALONE,
        "OCLN": OCLN_ALONE,
        "CLDN147_F11R_PARD3": CLDN147_F11R_PARD3,
        "TJ_7gene": USER_PPT_7GENE,
        "CD8": CD8,
        "CYT": CYT,
        "GEP": GEP,
    }
    scores = pd.DataFrame(index=expr.columns)
    used = {}
    for name, genes in feats.items():
        if name in ("CLDN4", "OCLN"):
            hit = resolve(genes, expr.index)
            scores[name] = expr.loc[hit[0]] if hit else np.nan
            used[name] = hit
        else:
            s, u = meanz(expr, genes)
            scores[name] = s
            used[name] = u
    return scores, used


def box_by_group(score: pd.Series, group: pd.Series, title: str, ylab: str, fname: Path, order: list[str]) -> None:
    fig, ax = plt.subplots(figsize=(4.2, 4.2))
    data = [score.loc[group.index[group == g]].dropna().values for g in order]
    ax.boxplot(data, tick_labels=order, showfliers=False, widths=0.55)
    rng = np.random.default_rng(0)
    colors = ["#2c7fb8", "#d95f0e"]
    for i, (vals, c) in enumerate(zip(data, colors), start=1):
        x = rng.normal(i, 0.06, size=len(vals))
        ax.scatter(x, vals, s=26, c=c, edgecolors="k", linewidths=0.3, zorder=3)
    ax.set_title(title, fontsize=10)
    ax.set_ylabel(ylab)
    fig.tight_layout()
    fig.savefig(fname.with_suffix(".png"), dpi=150)
    fig.savefig(fname.with_suffix(".pdf"))
    plt.close(fig)


def scatter_xy(x: pd.Series, y: pd.Series, title: str, xlab: str, ylab: str, fname: Path) -> None:
    d = pd.concat([x.rename("x"), y.rename("y")], axis=1).dropna()
    fig, ax = plt.subplots(figsize=(4.3, 4.3))
    ax.scatter(d["x"], d["y"], s=8, alpha=0.35, c="#225ea8", edgecolors="none")
    if len(d) >= 5:
        rho, p = stats.spearmanr(d["x"], d["y"])
        ax.set_title(f"{title}\nρ={rho:.2f} p={p:.2g} n={len(d)}", fontsize=10)
    else:
        ax.set_title(title)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    fig.tight_layout()
    fig.savefig(fname.with_suffix(".png"), dpi=140)
    fig.savefig(fname.with_suffix(".pdf"))
    plt.close(fig)


def main() -> int:
    ici_specs = [
        {
            "id": "GSE126044",
            "loader": load_gse126044,
            "note": "Index B4 cohort (Cho 2020). Pre-treatment NSCLC anti-PD-1. GEO responder/non-responder.",
            "low": "R",
            "high": "NR",
            "endpoint": "ORR (GEO responder vs non-responder)",
            "role": "index",
        },
        {
            "id": "GSE135222",
            "loader": load_gse135222,
            "note": "Jung/Kim 2020. Advanced NSCLC anti-PD-1/PD-L1. DCB = PFS ≥ 183 days from GEO pfs.time.",
            "low": "DCB",
            "high": "NDB",
            "endpoint": "DCB (PFS≥183d) / continuous PFS",
            "role": "extra",
        },
        {
            "id": "GSE207422",
            "loader": load_gse207422,
            "note": "Hu 2023. Neoadjuvant PD-1 + chemo. Pre-treatment bulk. MPR vs NMPR.",
            "low": "MPR",
            "high": "NMPR",
            "endpoint": "pathologic MPR vs NMPR",
            "role": "extra",
        },
        {
            "id": "GSE190265",
            "loader": load_gse190265,
            "note": "France3 / Leduc 2022-linked. NSCLC chemo±ICI biopsies. DCB = PFS ≥ 6 months.",
            "low": "DCB",
            "high": "NDB",
            "endpoint": "DCB (PFS≥6 mo) / continuous PFS",
            "role": "extra",
        },
        {
            "id": "GSE166449",
            "loader": load_gse166449,
            "note": "Lee/SMC immunotherapy lung RNA. GEO titles Responder vs nonResponder (7 vs 15).",
            "low": "R",
            "high": "NR",
            "endpoint": "GEO immunotherapy responder vs non-responder",
            "role": "extra",
        },
        {
            "id": "GSE190266",
            "loader": load_gse190266,
            "note": "France4 / Leduc 2022-linked. NSCLC ICI biopsies. Public TPM; PFS truncated at 6 months in GEO. Public matrix lacks TJP1/TJP2/OCLN/PARD3.",
            "low": "DCB",
            "high": "NDB",
            "endpoint": "DCB (PFS≥6 mo, GEO-capped) / continuous PFS",
            "role": "extra",
        },
        {
            "id": "GSE161537",
            "loader": load_gse161537,
            "note": "NivoBio targeted RNA (~2.5k genes). Pre-tx NSCLC PD-1/PD-L1. RECIST CR/PR vs PD. Panel has OCLN/F11R; no CLDN4/CLDN1/CLDN7/TJP.",
            "low": "R",
            "high": "NR",
            "endpoint": "ORR (GEO RECIST CR/PR vs PD)",
            "role": "extra",
        },
    ]

    ici_rows = []
    coverage = []
    sample_frames = []
    held: dict[str, tuple[pd.DataFrame, pd.DataFrame, dict]] = {}
    for spec in ici_specs:
        expr, meta = spec["loader"]()
        scores, used = score_features(expr)
        held[spec["id"]] = (scores, meta, used)
        scores = scores.join(meta[["group"]], how="inner")
        coverage.append({
            "cohort": spec["id"],
            "n": int(len(scores)),
            "endpoint": spec["endpoint"],
            "n_low": int((scores["group"] == spec["low"]).sum()),
            "n_high": int((scores["group"] == spec["high"]).sum()),
            "genes_CLDN4": ",".join(map(str, used["CLDN4"])),
            "genes_OCLN": ",".join(map(str, used["OCLN"])),
            "genes_5": len(used["CLDN147_F11R_PARD3"]),
            "genes_7": len(used["TJ_7gene"]),
            "note": spec["note"],
        })
        out_scores = scores.copy()
        out_scores.insert(0, "cohort", spec["id"])
        sample_frames.append(out_scores.reset_index().rename(columns={"index": "sample"}))

        for feat in ["CLDN4", "OCLN", "CLDN147_F11R_PARD3", "TJ_7gene", "CD8", "GEP"]:
            if scores[feat].notna().sum() < 4:
                continue
            a = scores.loc[scores["group"] == spec["low"], feat]
            b = scores.loc[scores["group"] == spec["high"], feat]
            mw = mwu_high_in_b(a, b)
            row = {
                "cohort": spec["id"],
                "role": spec["role"],
                "feature": feat,
                "n_genes": len(used[feat]),
                "endpoint": spec["endpoint"],
                "group_A": spec["low"],
                "group_B": spec["high"],
                "n_A": mw.get("n_A"),
                "n_B": mw.get("n_B"),
                "median_A": mw.get("median_A"),
                "median_B": mw.get("median_B"),
                "p_two": mw.get("p_two"),
                "direction": mw.get("direction"),
                "note": spec["note"],
            }
            # continuous PFS where available
            if "pfs_days" in meta.columns:
                sp = spearman(scores[feat], meta.loc[scores.index, "pfs_days"])
                row["spearman_pfs_rho"] = sp["rho"]
                row["spearman_pfs_p"] = sp["p"]
                row["spearman_pfs_n"] = sp["n"]
                row["pfs_unit"] = "days"
            if "pfs_months" in meta.columns:
                sp = spearman(scores[feat], meta.loc[scores.index, "pfs_months"])
                row["spearman_pfs_rho"] = sp["rho"]
                row["spearman_pfs_p"] = sp["p"]
                row["spearman_pfs_n"] = sp["n"]
                row["pfs_unit"] = "months"
            ici_rows.append(row)

        if spec["role"] == "extra":
            ep = spec["endpoint"].split("/")[0].strip()
            if len(used["TJ_7gene"]) >= 6:
                box_by_group(
                    scores["TJ_7gene"], scores["group"],
                    f"{spec['id']} TJ 7-gene vs {ep}",
                    "TJ 7-gene mean-z", FIG / f"{spec['id']}_TJ7",
                    [spec["low"], spec["high"]],
                )
            if used["CLDN4"]:
                box_by_group(
                    scores["CLDN4"], scores["group"],
                    f"{spec['id']} CLDN4 vs {ep}",
                    "CLDN4", FIG / f"{spec['id']}_CLDN4",
                    [spec["low"], spec["high"]],
                )
            elif used["OCLN"]:
                box_by_group(
                    scores["OCLN"], scores["group"],
                    f"{spec['id']} OCLN vs {ep}",
                    "OCLN", FIG / f"{spec['id']}_OCLN",
                    [spec["low"], spec["high"]],
                )
            if spec["id"] == "GSE190266" and len(used["CLDN147_F11R_PARD3"]) >= 3:
                box_by_group(
                    scores["CLDN147_F11R_PARD3"], scores["group"],
                    f"{spec['id']} CLDN1/4/7/F11R vs {ep}",
                    "CLDN1/4/7/F11R mean-z (PARD3 absent)", FIG / f"{spec['id']}_CLDN147",
                    [spec["low"], spec["high"]],
                )

    ici_df = pd.DataFrame(ici_rows)
    cov_df = pd.DataFrame(coverage)
    ici_df.to_csv(TAB / "lung_ICI_TJ_vs_response.tsv", sep="\t", index=False)
    cov_df.to_csv(TAB / "cohort_coverage.tsv", sep="\t", index=False)
    pd.concat(sample_frames, ignore_index=True).to_csv(TAB / "per_sample_scores.tsv", sep="\t", index=False)

    # PFS Cox / log-rank (public time + event only)
    surv_rows = []
    immune_rows = []
    for spec in ici_specs:
        scores, meta, used = held[spec["id"]]
        joined = scores.join(meta, how="inner", rsuffix="_meta")
        time_col = "pfs_days" if "pfs_days" in joined.columns else ("pfs_months" if "pfs_months" in joined.columns else None)
        event_col = "pfs_event" if "pfs_event" in joined.columns else None
        unit = "days" if time_col == "pfs_days" else ("months" if time_col == "pfs_months" else "")
        for feat in ["CLDN4", "OCLN", "CLDN147_F11R_PARD3", "TJ_7gene", "CD8", "GEP"]:
            if scores[feat].notna().sum() < 8:
                continue
            if time_col and event_col:
                cx = cox_per_sd(joined[time_col], joined[event_col], joined[feat])
                lr = logrank_median(joined[time_col], joined[event_col], joined[feat])
                surv_rows.append({
                    "cohort": spec["id"],
                    "feature": feat,
                    "n_genes": len(used[feat]),
                    "n": cx["n"],
                    "n_event": cx["n_event"],
                    "time_unit": unit,
                    "cox_hr_per_sd": cx.get("hr_per_sd"),
                    "cox_p": cx.get("p"),
                    "logrank_median_p": lr.get("p"),
                    "note": spec["note"],
                })
                if spec["role"] == "extra" and feat in ("TJ_7gene", "CLDN4", "OCLN") and (
                    (feat == "TJ_7gene" and len(used[feat]) >= 6)
                    or (feat == "CLDN4" and used[feat])
                    or (feat == "OCLN" and spec["id"] == "GSE161537")
                ):
                    km_plot(
                        joined[time_col], joined[event_col], joined[feat],
                        f"{spec['id']} {feat} median split PFS",
                        FIG / f"{spec['id']}_{feat}_KM",
                    )
            if feat in ("CLDN4", "OCLN", "CLDN147_F11R_PARD3", "TJ_7gene"):
                for imm in ["CD8", "GEP"]:
                    if scores[imm].notna().sum() < 8:
                        continue
                    sp = spearman(scores[feat], scores[imm])
                    immune_rows.append({
                        "cohort": spec["id"],
                        "n": sp["n"],
                        "TJ": feat,
                        "n_TJ_genes": len(used[feat]),
                        "immune": imm,
                        "spearman_rho": sp["rho"],
                        "spearman_p": sp["p"],
                    })
                if spec["role"] == "extra" and feat == "TJ_7gene" and len(used[feat]) >= 6:
                    scatter_xy(scores[feat], scores["CD8"], f"{spec['id']} TJ 7-gene vs CD8",
                               "TJ 7-gene mean-z", "CD8 mean-z", FIG / f"{spec['id']}_TJ7_vs_CD8")
                if spec["role"] == "extra" and feat == "CLDN4" and used[feat]:
                    scatter_xy(scores[feat], scores["CD8"], f"{spec['id']} CLDN4 vs CD8",
                               "CLDN4", "CD8 mean-z", FIG / f"{spec['id']}_CLDN4_vs_CD8")

    surv_df = pd.DataFrame(surv_rows)
    imm_df = pd.DataFrame(immune_rows)
    if len(surv_df):
        surv_df.to_csv(TAB / "lung_ICI_TJ_vs_PFS_cox.tsv", sep="\t", index=False)
    if len(imm_df):
        imm_df.to_csv(TAB / "lung_ICI_TJ_vs_CD8_GEP.tsv", sep="\t", index=False)

    # Forest of extra cohorts (plus index) for TJ_7gene and CLDN4
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 5.2), sharey=True)
    extras = ici_df[ici_df["feature"].isin(["TJ_7gene", "CLDN4"])].copy()
    cohorts = ["GSE135222", "GSE207422", "GSE190265", "GSE190266", "GSE166449", "GSE126044"]
    for ax, feat, title in zip(axes, ["TJ_7gene", "CLDN4"], ["TJ 7-gene (≥6 genes)", "CLDN4"]):
        sub = extras[extras["feature"] == feat].copy()
        if "n_genes" in sub.columns and feat == "TJ_7gene":
            sub = sub[sub["n_genes"] >= 6]
        sub = sub.set_index("cohort")
        keep = [c for c in cohorts if c in sub.index]
        sub = sub.loc[keep]
        y = np.arange(len(sub))
        ax.axvline(0.05, color="#999", ls=":", lw=1)
        ax.scatter(sub["p_two"], y, s=42, c=["#d95f0e" if c != "GSE126044" else "#2c7fb8" for c in sub.index])
        ax.set_yticks(y)
        labels = [f"{c} ({int(sub.loc[c,'n_A'])}/{int(sub.loc[c,'n_B'])})" for c in sub.index]
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_xscale("log")
        ax.set_xlabel("two-sided MW p")
        ax.set_title(title, fontsize=11)
        ax.set_xlim(1e-3, 1.05)
    fig.suptitle("Public lung ICI: TJ / CLDN4 vs response (extra cohorts + B4 index)", fontsize=11)
    fig.tight_layout()
    fig.savefig(FIG / "forest_ICI_p.png", dpi=150)
    fig.savefig(FIG / "forest_ICI_p.pdf")
    plt.close(fig)

    # TCGA
    tcga_rows = []
    for cohort, fname in [("TCGA-LUAD", "TCGA_LUAD_HiSeqV2.gz"), ("TCGA-LUSC", "TCGA_LUSC_HiSeqV2.gz")]:
        expr = load_tcga(EXTRA_DATA / fname)
        scores, used = score_features(expr)
        for tj in ["CLDN4", "CLDN147_F11R_PARD3", "TJ_7gene"]:
            for imm in ["CD8", "CYT", "GEP"]:
                sp = spearman(scores[tj], scores[imm])
                tcga_rows.append({
                    "cohort": cohort,
                    "n_tumors": int(len(scores)),
                    "TJ": tj,
                    "immune": imm,
                    "n": sp["n"],
                    "spearman_rho": sp["rho"],
                    "spearman_p": sp["p"],
                    "n_TJ_genes": len(used[tj]) if tj != "CLDN4" else 1,
                    "n_immune_genes": len(used[imm]),
                })
        scatter_xy(scores["TJ_7gene"], scores["CD8"], f"{cohort} TJ 7-gene vs CD8",
                   "TJ 7-gene mean-z", "CD8 mean-z", FIG / f"{cohort}_TJ7_vs_CD8")
        scatter_xy(scores["TJ_7gene"], scores["GEP"], f"{cohort} TJ 7-gene vs GEP",
                   "TJ 7-gene mean-z", "Ayers GEP mean-z", FIG / f"{cohort}_TJ7_vs_GEP")
        scatter_xy(scores["CLDN4"], scores["CD8"], f"{cohort} CLDN4 vs CD8",
                   "CLDN4", "CD8 mean-z", FIG / f"{cohort}_CLDN4_vs_CD8")

    tcga_df = pd.DataFrame(tcga_rows)
    tcga_df.to_csv(TAB / "TCGA_TJ_vs_CD8_GEP.tsv", sep="\t", index=False)

    write_extra(OUT, ici_df, cov_df, tcga_df, surv_df, imm_df)
    (OUT / "summary.json").write_text(json.dumps({
        "ici": ici_df.to_dict(orient="records"),
        "coverage": cov_df.to_dict(orient="records"),
        "tcga": tcga_df.to_dict(orient="records"),
        "pfs_cox": surv_df.to_dict(orient="records") if len(surv_df) else [],
        "ici_tj_vs_immune": imm_df.to_dict(orient="records") if len(imm_df) else [],
        "skipped": [
            {"id": "GSE136961", "reason": "Oncomine Immune Response 395-gene panel; no CLDN/TJ genes."},
            {"id": "GSE253564", "reason": "Pre-treatment FPKM public, but GEO has arm only (no MPR/ORR/PFS)."},
            {"id": "GSE182328", "reason": "Public counts, but GEO has Akkermansia detectability only (no ORR/PFS)."},
            {"id": "GSE248378", "reason": "Post-treatment FPKM; GEO has treatment arm only (no MPR/ORR/PFS)."},
            {"id": "GSE93157", "reason": "nCounter PanCancer Immune 730-gene; no CLDN4/TJ genes. Lung subset n=35."},
            {"id": "GSE110390", "reason": "Durvalumab 21-gene IFN panel; no CLDN/TJ genes."},
            {"id": "GSE162520", "reason": "Same targeted panel as GSE161537; TUMADOR early-stage surgical cohort, not ICI response."},
        ],
    }, indent=2, default=str))
    print(ici_df[ici_df["feature"].isin(["TJ_7gene", "CLDN4"])][
        ["cohort", "feature", "n_A", "n_B", "median_A", "median_B", "p_two", "direction"]
    ].to_string(index=False))
    print(tcga_df.to_string(index=False))
    return 0


def write_extra(
    out: Path,
    ici: pd.DataFrame,
    cov: pd.DataFrame,
    tcga: pd.DataFrame,
    surv: pd.DataFrame,
    imm: pd.DataFrame,
) -> None:
    def fmt(p) -> str:
        if p is None or (isinstance(p, float) and (np.isnan(p))):
            return "NA"
        return f"{p:.3g}"

    def ici_md(feat: str, min_genes: int = 1) -> str:
        lines = [
            "| cohort | endpoint | n A / B | genes | median A | median B | MW p | direction | PFS ρ (p) |",
            "|---|---|---|---|---|---|---|---|---|",
        ]
        sub = ici[ici["feature"] == feat]
        for _, r in sub.iterrows():
            ng = int(r["n_genes"]) if "n_genes" in r and pd.notna(r.get("n_genes")) else 0
            if ng and ng < min_genes:
                continue
            if pd.isna(r.get("p_two")):
                continue
            pfs = ""
            if pd.notna(r.get("spearman_pfs_rho")):
                pfs = f"ρ={r.spearman_pfs_rho:.2f} (p={fmt(r.spearman_pfs_p)}, n={int(r.spearman_pfs_n)})"
            lines.append(
                f"| {r.cohort} | {r.endpoint} | {int(r.n_A)} / {int(r.n_B)} | {ng or ''} | "
                f"{r.median_A:.3g} | {r.median_B:.3g} | {fmt(r.p_two)} | {r.direction} | {pfs} |"
            )
        return "\n".join(lines)

    def surv_md() -> str:
        if surv is None or len(surv) == 0:
            return "_No public PFS time+event in these matrices._"
        lines = [
            "| cohort | feature | genes | n (events) | Cox HR / SD | Cox p | log-rank median p |",
            "|---|---|---|---|---|---|---|",
        ]
        for _, r in surv.iterrows():
            if int(r.n_genes) < 1:
                continue
            if r.feature == "TJ_7gene" and int(r.n_genes) < 6:
                continue
            lines.append(
                f"| {r.cohort} | {r.feature} | {int(r.n_genes)} | {int(r.n)} ({int(r.n_event)}) | "
                f"{r.cox_hr_per_sd:.2f} | {fmt(r.cox_p)} | {fmt(r.logrank_median_p)} |"
            )
        return "\n".join(lines)

    def imm_md() -> str:
        if imm is None or len(imm) == 0:
            return ""
        lines = [
            "| cohort | n | TJ | immune | Spearman ρ | p |",
            "|---|---|---|---|---|---|",
        ]
        keep_tj = {"CLDN4", "TJ_7gene", "CLDN147_F11R_PARD3", "OCLN"}
        for _, r in imm.iterrows():
            if r.TJ not in keep_tj:
                continue
            if r.TJ == "TJ_7gene" and int(r.n_TJ_genes) < 6:
                continue
            if r.TJ == "CLDN147_F11R_PARD3" and int(r.n_TJ_genes) < 4:
                continue
            lines.append(
                f"| {r.cohort} | {int(r.n)} | {r.TJ} | {r.immune} | {r.spearman_rho:.2f} | {fmt(r.spearman_p)} |"
            )
        return "\n".join(lines)

    def tcga_md() -> str:
        lines = [
            "| cohort | n tumors | TJ | immune | Spearman ρ | p |",
            "|---|---|---|---|---|---|",
        ]
        for _, r in tcga.iterrows():
            lines.append(
                f"| {r.cohort} | {int(r.n_tumors)} | {r.TJ} | {r.immune} | {r.spearman_rho:.2f} | {fmt(r.spearman_p)} |"
            )
        return "\n".join(lines)

    cov_lines = [
        "| cohort | n | endpoint | note |",
        "|---|---|---|---|",
    ]
    for _, r in cov.iterrows():
        cov_lines.append(f"| {r.cohort} | {int(r.n)} | {r.endpoint} | {r.note} |")

    text = f"""# Extra paper figures: TJ / CLDN4 in additional public lung ICI cohorts + TCGA

Additive analyses on top of the GSE126044 B4 index cohort (that TJ/NR result is taken as given). Public data only. Honest n / medians / ρ / p / HR.

Signatures (same definitions as the B4 7-gene / 5-gene modules):

- **CLDN4** — single gene
- **OCLN** — single gene (used when a targeted panel lacks CLDN4)
- **CLDN1/4/7/F11R/PARD3** — 5-gene mean-z (score the genes present)
- **TJ 7-gene** — CLDN1, CLDN4, CLDN7, F11R, TJP1, TJP2, OCLN mean-z (tables below require ≥6 genes)
- **CD8** — CD8A+CD8B mean-z; **GEP** — Ayers 18-gene mean-z (immune context)

Group B is the poorer-outcome class (NR / NDB / NMPR) so direction **B>A** means higher TJ in non-responders, matching the B4 index direction.

## Extra lung ICI cohorts (not only GSE126044)

{chr(10).join(cov_lines)}

Not scored vs response: **GSE136961** (Oncomine 395-gene immune panel; no CLDN/TJ). **GSE253564** / **GSE248378** (GEO deposits treatment arm only). **GSE182328** (public counts; GEO has Akkermansia detectability, no ORR/PFS). **GSE93157** (nCounter immune 730-gene; no CLDN/TJ). **GSE110390** (21-gene IFN panel). **GSE162520** (same targeted panel as NivoBio, but TUMADOR early-stage surgical, not ICI).

### TJ 7-gene vs response (≥6 genes present)

{ici_md("TJ_7gene", min_genes=6)}

### CLDN4 vs response

{ici_md("CLDN4")}

### CLDN1/4/7/F11R/PARD3 vs response

{ici_md("CLDN147_F11R_PARD3", min_genes=4)}

### OCLN vs response (targeted-panel extra)

{ici_md("OCLN")}

### Immune context in the same ICI matrices

CD8 (CD8A+CD8B):

{ici_md("CD8")}

Ayers GEP:

{ici_md("GEP")}

## PFS Cox / log-rank (public time + event)

Cox HR is per +1 SD of the score (higher TJ → HR>1 means shorter PFS). Log-rank is a median split.

{surv_md()}

## TJ vs CD8 / GEP inside the ICI matrices

{imm_md()}

## TCGA supporting context: TJ vs CD8 / GEP

Xena `HiSeqV2` primary tumors (`-01`). Spearman, two-sided.

{tcga_md()}

## Figures

- Extra-cohort boxplots: `figures/GSE135222_TJ7.png`, `GSE207422_TJ7.png`, `GSE190265_TJ7.png`, `GSE166449_TJ7.png`, `GSE190266_CLDN4.png`, `GSE161537_OCLN.png` (and matching CLDN4 where the gene is present)
- `figures/forest_ICI_p.png` — two-sided MW p across extra cohorts + B4 index
- KM median-split PFS where GEO has time+event
- ICI and TCGA TJ vs CD8/GEP scatters

## Rerun

```bash
python3 scripts/rework/B4_wave2/download.py
python3 scripts/rework/B4_wave2/extra_cohorts.py
```
"""
    (out / "EXTRA_FOR_PAPER.md").write_text(text)
    (out / "README.md").write_text(text)


if __name__ == "__main__":
    raise SystemExit(main())
