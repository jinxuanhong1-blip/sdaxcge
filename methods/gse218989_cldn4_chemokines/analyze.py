#!/usr/bin/env python3
"""GSE218989 additive: CLDN4 vs CXCL9 / CXCL10 / CXCL13 and GEP-like.

Public SMC–KAIST PD-1/PD-L1 NSCLC bulk TPM only (Kang et al., Nat Commun 2024).
Patient unit. Does not headline CD274 / HLA (separate launch) and does not
re-audit methods/gse218989_cldn4_ici.
"""
from __future__ import annotations

import gzip
import json
import math
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
TABLES = HERE / "tables"
FIGURES = HERE / "figures"
CACHE = Path("/tmp/gse218989")

TPM_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE218nnn/GSE218989/suppl/"
    "GSE218989_SMC_KAIST_TPM_matrix_pc_rmdupli.txt.gz"
)
MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE218nnn/GSE218989/matrix/"
    "GSE218989_series_matrix.txt.gz"
)

# Ayers JCI 2017 public unweighted 18-gene T-cell-inflamed GEP
# (same list as data/signatures/gep18_genes.tsv on sibling branches).
GEP18 = [
    "CCL5",
    "CD27",
    "CD274",
    "CD276",
    "CD8A",
    "CMKLR1",
    "CXCL9",
    "CXCR6",
    "HLA-DQA1",
    "HLA-DRB1",
    "HLA-E",
    "IDO1",
    "LAG3",
    "NKG7",
    "PDCD1LG2",
    "PSMB10",
    "STAT1",
    "TIGIT",
]
CHEMOKINES = ["CXCL9", "CXCL10", "CXCL13"]
# CXCL10 is in Ayers IFN-γ 6-gene, not in GEP18. Keep that explicit.
AYERS6 = ["IDO1", "CXCL10", "CXCL9", "HLA-DRA", "STAT1", "IFNG"]


def download(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 1000:
        return dest
    print(f"download {url} -> {dest}", flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=180) as r, open(dest, "wb") as f:
        f.write(r.read())
    return dest


def load_tpm(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0, compression="gzip")
    df.index = df.index.astype(str)
    df = df[~df.index.duplicated(keep="first")]
    return df.astype(float)


def load_geo_meta(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt") as f:
        lines = f.read().splitlines()
    titles = None
    cols: dict[str, list[str]] = {}
    for line in lines:
        if line.startswith("!Sample_title"):
            titles = [x.strip('"') for x in line.split("\t")[1:]]
        elif line.startswith("!Sample_geo_accession"):
            cols["gsm"] = [x.strip('"') for x in line.split("\t")[1:]]
        elif line.startswith("!Sample_characteristics_ch1"):
            vals = [x.strip('"') for x in line.split("\t")[1:]]
            key = vals[0].split(":", 1)[0].strip()
            cols[key] = [v.split(":", 1)[1].strip() if ":" in v else v for v in vals]
    return pd.DataFrame({"patient": titles, **cols})


def mean_z(logx: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, list[str]]:
    present = [g for g in genes if g in logx.index]
    if not present:
        return pd.Series(np.nan, index=logx.columns), present
    z = logx.loc[present].T.apply(
        lambda s: (s - s.mean()) / (s.std(ddof=0) or 1.0), axis=0
    )
    return z.mean(axis=1), present


def spearman(a: pd.Series, b: pd.Series) -> dict:
    m = pd.concat([a, b], axis=1).dropna()
    if len(m) < 5:
        return {"n": int(len(m)), "rho": np.nan, "p": np.nan}
    r, p = stats.spearmanr(m.iloc[:, 0], m.iloc[:, 1])
    return {"n": int(len(m)), "rho": float(r), "p": float(p)}


def mwu_auc(score: pd.Series, y: pd.Series) -> dict:
    m = pd.concat([score, y], axis=1).dropna()
    m.columns = ["s", "y"]
    a = m.loc[m.y == 1, "s"]
    b = m.loc[m.y == 0, "s"]
    if len(a) < 2 or len(b) < 2:
        return {
            "n_pos": int(len(a)),
            "n_neg": int(len(b)),
            "median_pos": np.nan,
            "median_neg": np.nan,
            "U": np.nan,
            "p": np.nan,
            "auc": np.nan,
        }
    U, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    return {
        "n_pos": int(len(a)),
        "n_neg": int(len(b)),
        "median_pos": float(a.median()),
        "median_neg": float(b.median()),
        "U": float(U),
        "p": float(p),
        "auc": float(U / (len(a) * len(b))),
    }


def fmt_p(p: float) -> str:
    if p is None or (isinstance(p, float) and math.isnan(p)):
        return "NA"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.4g}"


def fmt_num(x: float, nd: int = 3) -> str:
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "NA"
    return f"{x:.{nd}f}"


def fmt_rho_p(sp: dict) -> str:
    return f"{fmt_num(sp.get('rho', np.nan))} ({fmt_p(sp.get('p', np.nan))})"


def bh_q(pvals: list[float]) -> list[float]:
    n = len(pvals)
    order = np.argsort(pvals)
    q = np.empty(n, dtype=float)
    prev = 1.0
    for rank, i in enumerate(order[::-1], start=0):
        adj = pvals[i] * n / (n - rank)
        prev = min(prev, adj)
        q[i] = prev
    return [float(min(1.0, x)) for x in q]


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(parents=True, exist_ok=True)

    tpm_path = download(TPM_URL, CACHE / "GSE218989_SMC_KAIST_TPM_matrix_pc_rmdupli.txt.gz")
    mtx_path = download(MATRIX_URL, CACHE / "GSE218989_series_matrix.txt.gz")

    tpm = load_tpm(tpm_path)
    geo = load_geo_meta(mtx_path).set_index("patient")
    patients = list(tpm.columns)
    if len(patients) != len(set(patients)):
        raise SystemExit("duplicate patient columns")
    if set(patients) != set(geo.index):
        raise SystemExit(
            f"TPM vs GEO ID mismatch: tpm={len(patients)} geo={len(geo)} "
            f"overlap={len(set(patients) & set(geo.index))}"
        )

    logx = np.log2(tpm + 1.0)
    needed = ["CLDN4", "TACSTD2", *CHEMOKINES, *GEP18, *AYERS6]
    presence = pd.DataFrame(
        [
            {
                "gene": g,
                "present": g in tpm.index,
                "n_patients_if_present": int(tpm.shape[1]) if g in tpm.index else 0,
                "median_tpm": float(tpm.loc[g].median()) if g in tpm.index else np.nan,
                "n_zero_tpm": int((tpm.loc[g] == 0).sum()) if g in tpm.index else np.nan,
            }
            for g in dict.fromkeys(needed)
        ]
    )
    presence.to_csv(TABLES / "gene_presence.tsv", sep="\t", index=False)

    missing = presence.loc[~presence["present"], "gene"].tolist()
    if "CLDN4" in missing:
        raise SystemExit("CLDN4 missing from TPM matrix")

    cldn4 = logx.loc["CLDN4"]
    tacstd2 = logx.loc["TACSTD2"] if "TACSTD2" in logx.index else pd.Series(np.nan, index=patients)
    chemokine_scores = {g: logx.loc[g] for g in CHEMOKINES if g in logx.index}
    gep, gep_used = mean_z(logx, GEP18)
    cx_trio, cx_used = mean_z(logx, CHEMOKINES)
    ayers, ayers_used = mean_z(logx, AYERS6)
    gep_no_cxcl9, gep_no_cxcl9_used = mean_z(logx, [g for g in GEP18 if g != "CXCL9"])

    geo_r = (
        geo.loc[patients, "treatment outcome"]
        .map({"Responder": 1, "Non-responder": 0})
        .astype(int)
    )
    n_r = int((geo_r == 1).sum())
    n_nr = int((geo_r == 0).sum())

    rows: list[dict] = []

    def add_spearman(feature: str, axis: str, a: pd.Series, b: pd.Series, family: str) -> dict:
        sp = spearman(a, b)
        rec = {
            "family": family,
            "feature": feature,
            "axis": axis,
            "n": sp["n"],
            "rho": sp["rho"],
            "p": sp["p"],
        }
        rows.append(rec)
        return rec

    primary_pairs = []
    for g in CHEMOKINES:
        if g in chemokine_scores:
            rec = add_spearman("CLDN4", g, cldn4, chemokine_scores[g], "primary")
            primary_pairs.append(rec)
        else:
            rows.append(
                {
                    "family": "primary",
                    "feature": "CLDN4",
                    "axis": g,
                    "n": 0,
                    "rho": np.nan,
                    "p": np.nan,
                }
            )
    rec_gep = add_spearman("CLDN4", "GEP18_meanz", cldn4, gep, "primary")
    primary_pairs.append(rec_gep)
    add_spearman("CLDN4", "CXCL9_10_13_meanz", cldn4, cx_trio, "extra")

    # companions — not headline
    if "TACSTD2" in logx.index:
        for g, s in chemokine_scores.items():
            add_spearman("TACSTD2", g, tacstd2, s, "companion")
        add_spearman("TACSTD2", "GEP18_meanz", tacstd2, gep, "companion")
        add_spearman("TACSTD2", "CXCL9_10_13_meanz", tacstd2, cx_trio, "companion")
        add_spearman("CLDN4", "TACSTD2", cldn4, tacstd2, "companion")
    add_spearman("CLDN4", "Ayers6_meanz", cldn4, ayers, "companion")
    add_spearman("CLDN4", "GEP18_noCXCL9_meanz", cldn4, gep_no_cxcl9, "sensitivity")
    if "CXCL9" in chemokine_scores and "CXCL10" in chemokine_scores:
        add_spearman("CXCL9", "CXCL10", chemokine_scores["CXCL9"], chemokine_scores["CXCL10"], "control")
    if "CXCL9" in chemokine_scores:
        add_spearman("CXCL9", "GEP18_meanz", chemokine_scores["CXCL9"], gep, "control")

    primary_p = [r["p"] for r in rows if r["family"] == "primary" and pd.notna(r["p"])]
    primary_idx = [i for i, r in enumerate(rows) if r["family"] == "primary" and pd.notna(r["p"])]
    qvals = bh_q(primary_p) if primary_p else []
    for r in rows:
        r["q_bh_primary"] = np.nan
    for i, q in zip(primary_idx, qvals):
        rows[i]["q_bh_primary"] = q

    stats_df = pd.DataFrame(rows)
    stats_df.to_csv(TABLES / "stats.tsv", sep="\t", index=False)

    mwu_rows = []
    for name, s in [
        ("CLDN4", cldn4),
        ("TACSTD2", tacstd2),
        *[(g, chemokine_scores[g]) for g in CHEMOKINES if g in chemokine_scores],
        ("GEP18_meanz", gep),
        ("CXCL9_10_13_meanz", cx_trio),
        ("Ayers6_meanz", ayers),
    ]:
        mw = mwu_auc(s, geo_r)
        mwu_rows.append({"feature": name, **mw, "n_total": mw["n_pos"] + mw["n_neg"]})
    mwu_df = pd.DataFrame(mwu_rows)
    mwu_df.to_csv(TABLES / "response_mwu_note.tsv", sep="\t", index=False)

    q = pd.qcut(cldn4, 4, labels=["Q1", "Q2", "Q3", "Q4"])
    q_n = q.value_counts().sort_index()
    qtab = pd.DataFrame({"quartile": q_n.index.astype(str), "n": q_n.values})
    qtab.to_csv(TABLES / "quartile_n.tsv", sep="\t", index=False)
    q4 = q == "Q4"
    q1 = q == "Q1"
    q_rows = []
    for name, s in [
        *[(g, chemokine_scores[g]) for g in CHEMOKINES if g in chemokine_scores],
        ("GEP18_meanz", gep),
        ("CXCL9_10_13_meanz", cx_trio),
    ]:
        y = pd.Series(np.where(q4, 1, np.where(q1, 0, np.nan)), index=cldn4.index)
        mw = mwu_auc(s, y)
        q_rows.append({"feature": name, "contrast": "CLDN4_Q4_vs_Q1", **mw})
    pd.DataFrame(q_rows).to_csv(TABLES / "quartile_mwu.tsv", sep="\t", index=False)

    def pick(feature: str, axis: str) -> dict:
        hit = stats_df[(stats_df.feature == feature) & (stats_df.axis == axis)]
        return hit.iloc[0].to_dict() if len(hit) else {"n": 0, "rho": np.nan, "p": np.nan}

    c4_cxcl9 = pick("CLDN4", "CXCL9")
    c4_cxcl10 = pick("CLDN4", "CXCL10")
    c4_cxcl13 = pick("CLDN4", "CXCL13")
    c4_gep = pick("CLDN4", "GEP18_meanz")
    c4_trio = pick("CLDN4", "CXCL9_10_13_meanz")

    one = pd.DataFrame(
        [
            {
                "dataset": "GSE218989 SMC-KAIST",
                "n_matrix_patients": int(tpm.shape[1]),
                "n_genes": int(tpm.shape[0]),
                "n_spearman": int(c4_gep.get("n", 0)),
                "CXCL9_present": "CXCL9" in tpm.index,
                "CXCL10_present": "CXCL10" in tpm.index,
                "CXCL13_present": "CXCL13" in tpm.index,
                "GEP18_n_present": len(gep_used),
                "GEP18_n_listed": len(GEP18),
                "GEP18_missing": ",".join([g for g in GEP18 if g not in gep_used]) or "none",
                "CLDN4_vs_CXCL9_n": int(c4_cxcl9.get("n", 0)),
                "CLDN4_vs_CXCL9_rho": c4_cxcl9.get("rho", np.nan),
                "CLDN4_vs_CXCL9_p": c4_cxcl9.get("p", np.nan),
                "CLDN4_vs_CXCL10_n": int(c4_cxcl10.get("n", 0)),
                "CLDN4_vs_CXCL10_rho": c4_cxcl10.get("rho", np.nan),
                "CLDN4_vs_CXCL10_p": c4_cxcl10.get("p", np.nan),
                "CLDN4_vs_CXCL13_n": int(c4_cxcl13.get("n", 0)),
                "CLDN4_vs_CXCL13_rho": c4_cxcl13.get("rho", np.nan),
                "CLDN4_vs_CXCL13_p": c4_cxcl13.get("p", np.nan),
                "CLDN4_vs_GEP18_n": int(c4_gep.get("n", 0)),
                "CLDN4_vs_GEP18_rho": c4_gep.get("rho", np.nan),
                "CLDN4_vs_GEP18_p": c4_gep.get("p", np.nan),
                "CLDN4_vs_CXCL9_10_13_rho": c4_trio.get("rho", np.nan),
                "CLDN4_vs_CXCL9_10_13_p": c4_trio.get("p", np.nan),
                "geo_R": n_r,
                "geo_NR": n_nr,
                "histology": "not_deposited",
            }
        ]
    )
    one.to_csv(TABLES / "one_row.tsv", sep="\t", index=False)

    coverage = pd.DataFrame(
        [
            {
                "set": "CXCL9_CXCL10_CXCL13",
                "n_listed": 3,
                "n_present": len(cx_used),
                "genes": ",".join(cx_used),
                "missing": ",".join([g for g in CHEMOKINES if g not in cx_used]) or "none",
            },
            {
                "set": "GEP18_Ayers2017",
                "n_listed": len(GEP18),
                "n_present": len(gep_used),
                "genes": ",".join(gep_used),
                "missing": ",".join([g for g in GEP18 if g not in gep_used]) or "none",
            },
            {
                "set": "GEP18_noCXCL9",
                "n_listed": len(GEP18) - 1,
                "n_present": len(gep_no_cxcl9_used),
                "genes": ",".join(gep_no_cxcl9_used),
                "missing": ",".join(
                    [g for g in GEP18 if g != "CXCL9" and g not in gep_no_cxcl9_used]
                )
                or "none",
            },
            {
                "set": "Ayers6_companion",
                "n_listed": len(AYERS6),
                "n_present": len(ayers_used),
                "genes": ",".join(ayers_used),
                "missing": ",".join([g for g in AYERS6 if g not in ayers_used]) or "none",
            },
        ]
    )
    coverage.to_csv(TABLES / "signature_coverage.tsv", sep="\t", index=False)

    inventory = pd.DataFrame(
        [
            {
                "field": "TPM_patients",
                "present": True,
                "n": int(tpm.shape[1]),
                "notes": "GSE218989_SMC_KAIST_TPM_matrix_pc_rmdupli.txt.gz; paper text 497, public matrix 355",
            },
            {
                "field": "GEO_treatment_outcome",
                "present": True,
                "n": int(tpm.shape[1]),
                "notes": f"Responder {n_r} / Non-responder {n_nr}; series matrix",
            },
            {
                "field": "GEO_histology",
                "present": False,
                "n": 0,
                "notes": "not a series-matrix characteristic; LUAD vs LUSC not split",
            },
            {
                "field": "CXCL9",
                "present": "CXCL9" in tpm.index,
                "n": int(tpm.shape[1]) if "CXCL9" in tpm.index else 0,
                "notes": "single-gene log2(TPM+1)",
            },
            {
                "field": "CXCL10",
                "present": "CXCL10" in tpm.index,
                "n": int(tpm.shape[1]) if "CXCL10" in tpm.index else 0,
                "notes": "in Ayers6; not in GEP18",
            },
            {
                "field": "CXCL13",
                "present": "CXCL13" in tpm.index,
                "n": int(tpm.shape[1]) if "CXCL13" in tpm.index else 0,
                "notes": "TLS-associated chemokine; not in GEP18 or Ayers6",
            },
            {
                "field": "GEP18",
                "present": len(gep_used) > 0,
                "n": int(tpm.shape[1]) if gep_used else 0,
                "notes": f"{len(gep_used)}/{len(GEP18)} genes; unweighted mean-z; includes CD274 as 1/18, not isolated here",
            },
        ]
    )
    inventory.to_csv(TABLES / "label_inventory.tsv", sep="\t", index=False)

    per = pd.DataFrame(
        {
            "patient": patients,
            "CLDN4_log2tpm1": cldn4.values,
            "TACSTD2_log2tpm1": tacstd2.values,
            "CXCL9_log2tpm1": chemokine_scores["CXCL9"].values if "CXCL9" in chemokine_scores else np.nan,
            "CXCL10_log2tpm1": chemokine_scores["CXCL10"].values if "CXCL10" in chemokine_scores else np.nan,
            "CXCL13_log2tpm1": chemokine_scores["CXCL13"].values if "CXCL13" in chemokine_scores else np.nan,
            "GEP18_meanz": gep.values,
            "CXCL9_10_13_meanz": cx_trio.values,
            "Ayers6_meanz": ayers.values,
            "CLDN4_quartile": q.astype(str).values,
            "GEO_responder": geo_r.values,
        }
    )
    per.to_csv(TABLES / "per_patient.tsv", sep="\t", index=False)

    summary = {
        "accession": "GSE218989",
        "matrix": "GSE218989_SMC_KAIST_TPM_matrix_pc_rmdupli.txt.gz",
        "n_patients": int(len(patients)),
        "n_genes": int(tpm.shape[0]),
        "missing_requested_genes": missing,
        "GEP18_coverage": f"{len(gep_used)}/{len(GEP18)}",
        "chemokine_coverage": f"{len(cx_used)}/3",
        "histology": "not_deposited",
        "geo_R_NR": [n_r, n_nr],
        "one_row": one.iloc[0].to_dict(),
        "primary_spearman": stats_df[stats_df.family == "primary"].to_dict(orient="records"),
        "note": "Additive chemokine/GEP page. CD274/HLA isolated tests are out of scope.",
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2, default=str))

    plt.rcParams.update(
        {
            "figure.dpi": 140,
            "savefig.dpi": 160,
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )

    fig, axes = plt.subplots(2, 2, figsize=(8.6, 7.6))
    panels = [
        (chemokine_scores.get("CXCL9"), "CXCL9 log2(TPM+1)", c4_cxcl9),
        (chemokine_scores.get("CXCL10"), "CXCL10 log2(TPM+1)", c4_cxcl10),
        (chemokine_scores.get("CXCL13"), "CXCL13 log2(TPM+1)", c4_cxcl13),
        (gep, "GEP18 mean-z (Ayers 2017)", c4_gep),
    ]
    for ax, (y, ylab, sp) in zip(axes.ravel(), panels):
        if y is None:
            ax.text(0.5, 0.5, f"{ylab}\nnot on matrix", ha="center", va="center")
            ax.set_axis_off()
            continue
        ax.scatter(cldn4, y, s=10, alpha=0.4, c="#4C78A8", linewidths=0)
        ax.set_xlabel("CLDN4 log2(TPM+1)")
        ax.set_ylabel(ylab)
        ax.set_title(
            f"ρ={fmt_num(sp.get('rho', np.nan))}  p={fmt_p(sp.get('p', np.nan))}  n={int(sp.get('n', 0))}",
            fontsize=9,
        )
    fig.suptitle("GSE218989 CLDN4 vs chemokines / GEP-like (patient unit)", y=1.01)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig1_cldn4_chemokine_scatter.png", bbox_inches="tight")
    plt.close(fig)

    prim = stats_df[stats_df.family == "primary"].copy()
    fig, ax = plt.subplots(figsize=(6.8, 3.6))
    y = np.arange(len(prim))
    ax.axvline(0, color="0.5", lw=1)
    colors = ["#4C78A8" if (isinstance(r, float) and r < 0) else "#E45756" for r in prim["rho"]]
    for i, (_, r) in enumerate(prim.iterrows()):
        ax.plot(r["rho"], i, "o", color=colors[i], ms=8)
        ax.plot([0, r["rho"]], [i, i], color=colors[i], lw=2)
    ax.set_yticks(y, [f"CLDN4 vs {a}" for a in prim["axis"]], fontsize=8)
    ax.set_xlabel(f"Spearman ρ (n={int(c4_gep.get('n', 0))})")
    ax.set_title("Primary chemokine / GEP-like axes")
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(FIGURES / "fig2_spearman_forest.png", bbox_inches="tight")
    plt.close(fig)

    print("wrote", TABLES)
    print("wrote", FIGURES)
    print(one.T.to_string())
    print(stats_df[stats_df.family == "primary"].to_string(index=False))


if __name__ == "__main__":
    main()
