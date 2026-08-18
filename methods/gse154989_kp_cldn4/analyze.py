#!/usr/bin/env python3
"""ADDITIVE public mouse: GSE154989 KP GEMM Smart-seq2, Cldn4-only.

Public processed matrix only. Epithelial / malignant Cldn4 vs IFN/MHC/TJ
at mouse level. T/NK fraction is a design no-go (FACS CD45−).

Thesis is not rewritten. No dual-high TACSTD2×CLDN4.
"""
from __future__ import annotations

import json
import math
import re
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
TABLES = HERE / "tables"
FIGS = HERE / "figures"

GEO_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154989/suppl"
FILES = {
    "h5": "GSE154989_mmLungPlate_fQC_dSp_normTPM.h5",
    "smp": "GSE154989_mmLungPlate_fQC_smpTable.csv.gz",
    "annot": "GSE154989_mmLungPlate_fQC_dZ_annot_smpTable.csv.gz",
    "gene": "GSE154989_mmLungPlate_fQC_geneTable.csv.gz",
}

MIN_CELLS_PRIMARY = 20
MIN_N_SPEARMAN = 4
MIN_N_Q4 = 8
MOUSE_RE = re.compile(r"^([A-Za-z]+)_(\d+w)_ND_([mf]\d+)(?:_T(\d+))?$")

# Curated mouse MHC-I / APM (classical + nonclassical + machinery).
# Human HLA-* genes in A8 have no 1:1 title-case mouse name.
MOUSE_MHC_I_APM = [
    "H2-K1",
    "H2-D1",
    "H2-Q1",
    "H2-Q2",
    "H2-Q4",
    "H2-Q6",
    "H2-Q7",
    "H2-Q10",
    "H2-T23",
    "H2-M3",
    "B2m",
    "Tap1",
    "Tap2",
    "Tapbp",
    "Tapbpl",
    "Nlrc5",
    "Psmb8",
    "Psmb9",
    "Psmb10",
    "Erap1",
    "Calr",
    "Canx",
    "Pdia3",
    "Irf1",
]

# Human A8 leftovers that title-case misses.
HUMAN_TO_MOUSE = {
    "WARS1": "Wars",
    "C1R": "C1ra",
    "C1S": "C1s1",
    "FCGR1A": "Fcgr1",
    "HLA-A": "H2-K1",
    "HLA-B": "H2-D1",
    "HLA-C": "H2-Q4",
    "HLA-E": "H2-T23",
    "HLA-F": "H2-Q10",
    "HLA-G": "H2-Q6",
    "HLA-DMA": "H2-DMa",
    "HLA-DQA1": "H2-Aa",
    "HLA-DRB1": "H2-Eb1",
    "SECTM1": "Sectm1a",
    "TENT5A": "Tent5a",
    "MARCHF1": "Marchf1",
}

TNK_AUDIT = ["Ptprc", "Cd3d", "Cd3e", "Cd3g", "Cd2", "Cd8a", "Cd8b1", "Nkg7", "Klrd1", "Gzma"]
EPI_AUDIT = ["Cldn4", "Epcam", "Krt8", "Krt18", "Sftpc", "Nkx2-1"]


def _fmt(x, nd=3):
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return "—"
    return f"{x:.{nd}f}"


def download(name: str) -> Path:
    dest = DATA / name
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 1000:
        return dest
    gz = dest if name.endswith(".gz") else None
    url = f"{GEO_BASE}/{name}"
    print("DOWNLOAD", url)
    urllib.request.urlretrieve(url, dest)
    return dest


def human_to_mouse(sym: str) -> str:
    if sym in HUMAN_TO_MOUSE:
        return HUMAN_TO_MOUSE[sym]
    if sym.startswith("HLA-"):
        return ""
    if sym.isupper() or (sym.replace("-", "").replace("_", "").isupper()):
        parts = []
        for p in re.split(r"(-)", sym):
            if p == "-":
                parts.append(p)
            elif p:
                parts.append(p[0] + p[1:].lower())
        return "".join(parts)
    return sym


def load_families(available: set[str]) -> dict[str, list[str]]:
    a8 = json.loads((DATA / "a8_families.json").read_text())
    sets = a8["sets"]
    ifn_h = set(sets["HALLMARK_INTERFERON_GAMMA_RESPONSE"]) | set(
        sets["HALLMARK_INTERFERON_ALPHA_RESPONSE"]
    )
    ifn = []
    seen = set()
    for g in sorted(ifn_h):
        m = human_to_mouse(g)
        if m and m in available and m not in seen:
            ifn.append(m)
            seen.add(m)

    mhc = [g for g in MOUSE_MHC_I_APM if g in available]

    tj_h = set(sets["KEGG_TIGHT_JUNCTION"]) | set(sets["GOBP_TIGHT_JUNCTION_ORGANIZATION"])
    for g in a8.get("focal_genes", []):
        if g not in sets.get("KRT_EPITHELIAL", []):
            tj_h.add(g)
    tj_h.discard("CLDN4")
    tj = []
    seen = set()
    for g in sorted(tj_h):
        m = human_to_mouse(g)
        if m == "Cldn4":
            continue
        if m and m in available and m not in seen:
            tj.append(m)
            seen.add(m)

    return {"IFN": ifn, "MHC-I/APM": mhc, "TJ": tj}


def parse_mouse(mouse_id: str) -> dict:
    m = MOUSE_RE.match(str(mouse_id))
    if not m:
        return {
            "mouse_id": mouse_id,
            "animal": mouse_id,
            "genotype": "?",
            "week": "",
            "tumor": "",
        }
    geno, week, animal_n, tumor = m.group(1), m.group(2), m.group(3), m.group(4)
    animal = f"{geno}_{week}_ND_{animal_n}"
    return {
        "mouse_id": mouse_id,
        "animal": animal,
        "genotype": geno,
        "week": week,
        "tumor": f"T{tumor}" if tumor else "T0",
    }


def q4_vs_q1(cldn4, y, min_n: int = MIN_N_Q4) -> dict:
    s = pd.DataFrame({"c": np.asarray(cldn4, float), "y": np.asarray(y, float)})
    s = s[np.isfinite(s["c"]) & np.isfinite(s["y"])].copy()
    n = int(len(s))
    out = {
        "n": n,
        "n_q1": np.nan,
        "n_q4": np.nan,
        "median_q1": np.nan,
        "median_q4": np.nan,
        "delta_median": np.nan,
        "r_rb": np.nan,
        "p": np.nan,
        "usable": False,
    }
    if n < min_n:
        return out
    ranks = s["c"].rank(method="average")
    try:
        qs = pd.qcut(ranks, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")
    except ValueError:
        return out
    if qs.nunique() < 4:
        return out
    q1 = s.loc[qs == "Q1", "y"]
    q4 = s.loc[qs == "Q4", "y"]
    n1, n4 = int(len(q1)), int(len(q4))
    if n1 < 2 or n4 < 2:
        return out
    u, p = stats.mannwhitneyu(q4.values, q1.values, alternative="two-sided")
    r_rb = (2.0 * float(u)) / (n4 * n1) - 1.0
    out.update(
        {
            "n_q1": n1,
            "n_q4": n4,
            "median_q1": float(q1.median()),
            "median_q4": float(q4.median()),
            "delta_median": float(q4.median() - q1.median()),
            "r_rb": float(r_rb),
            "p": float(p),
            "usable": True,
        }
    )
    return out


def spearman(x, y) -> dict:
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    n = int(m.sum())
    if n < MIN_N_SPEARMAN:
        return {"n": n, "rho": np.nan, "p": np.nan, "usable": False}
    rho, p = stats.spearmanr(x[m], y[m])
    return {"n": n, "rho": float(rho), "p": float(p), "usable": True}


def load_sparse_matrix():
    import gzip

    import h5py

    gene_path = download(FILES["gene"])
    if str(gene_path).endswith(".gz"):
        with gzip.open(gene_path, "rt") as fh:
            gene_df = pd.read_csv(fh)
    else:
        gene_df = pd.read_csv(gene_path)
    symbols = gene_df["geneSymbol"].astype(str).tolist()
    assert len(symbols) == 52638

    smp_path = download(FILES["smp"])
    if str(smp_path).endswith(".gz"):
        with gzip.open(smp_path, "rt") as fh:
            smp = pd.read_csv(fh)
    else:
        smp = pd.read_csv(smp_path)
    assert len(smp) == 3891

    annot_path = download(FILES["annot"])
    if str(annot_path).endswith(".gz"):
        with gzip.open(annot_path, "rt") as fh:
            annot = pd.read_csv(fh)
    else:
        annot = pd.read_csv(annot_path)

    h5_path = download(FILES["h5"])
    with h5py.File(h5_path, "r") as f:
        ii = f["i"][0].astype(np.int32)
        jj = f["j"][0].astype(np.int32)
        vv = f["v"][0].astype(np.float64)
    return symbols, smp, annot, ii, jj, vv


def extract_needed(symbols, needed, ii, jj, vv, n_cells: int) -> dict[str, np.ndarray]:
    want_idx = {i + 1: g for i, g in enumerate(symbols) if g in needed}
    out = {g: np.zeros(n_cells, dtype=np.float64) for g in needed if g in set(symbols)}
    for k in range(len(ii)):
        g = want_idx.get(int(ii[k]))
        if g is None:
            continue
        out[g][int(jj[k]) - 1] = vv[k]
    return out


def mean_score(vecs: dict[str, np.ndarray], genes: list[str], n_cells: int) -> np.ndarray:
    present = [g for g in genes if g in vecs]
    if not present:
        return np.full(n_cells, np.nan)
    return np.mean(np.vstack([vecs[g] for g in present]), axis=0)


def cohort_mask(meta: pd.DataFrame, name: str) -> pd.Series:
    n = meta["n_cells"] >= MIN_CELLS_PRIMARY
    if name == "KP_n20":
        return n & (meta["genotype"] == "KP")
    if name == "KplusKP_n20":
        return n & meta["genotype"].isin(["K", "KP"])
    if name == "all_n20":
        return n
    if name == "KP_all":
        return meta["genotype"] == "KP"
    raise KeyError(name)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)

    symbols, smp, annot, ii, jj, vv = load_sparse_matrix()
    available = set(symbols)
    families = load_families(available)

    needed = {"Cldn4"}
    needed.update(EPI_AUDIT)
    needed.update(TNK_AUDIT)
    for genes in families.values():
        needed.update(genes)

    vecs = extract_needed(symbols, needed, ii, jj, vv, n_cells=len(smp))
    missing_audit = [g for g in (EPI_AUDIT + TNK_AUDIT) if g not in vecs]
    print("missing audit genes", missing_audit)
    for fam, genes in families.items():
        print(f"family {fam}: {len(genes)} present")

    parsed = smp["mouseID"].map(parse_mouse).apply(pd.Series)
    meta_cells = pd.concat([smp.reset_index(drop=True), parsed], axis=1)
    meta_cells["clusterK12"] = annot.set_index("sampleID").loc[smp["sampleID"], "clusterK12"].to_numpy()
    meta_cells["Cldn4"] = vecs["Cldn4"]
    meta_cells["Cldn4_pos"] = (vecs["Cldn4"] > 0).astype(int)
    for fam, genes in families.items():
        meta_cells[fam] = mean_score(vecs, genes, len(smp))
    for g in TNK_AUDIT:
        if g in vecs:
            meta_cells[g] = vecs[g]
    for g in EPI_AUDIT:
        if g in vecs:
            meta_cells[g] = vecs[g]

    # Cell-level immune leak audit (not a T/NK fraction).
    tnk_present = [g for g in TNK_AUDIT if g in vecs]
    leak = {
        "n_cells": int(len(meta_cells)),
        "Cldn4_pos": int((meta_cells["Cldn4"] > 0).sum()),
        "Epcam_pos": int((meta_cells["Epcam"] > 0).sum()) if "Epcam" in meta_cells else None,
        "Ptprc_pos": int((meta_cells["Ptprc"] > 0).sum()) if "Ptprc" in meta_cells else None,
        "any_tnk_marker_pos": int((meta_cells[tnk_present] > 0).any(axis=1).sum()),
        "Cd3d_pos": int((meta_cells["Cd3d"] > 0).sum()) if "Cd3d" in meta_cells else None,
        "Nkg7_pos": int((meta_cells["Nkg7"] > 0).sum()) if "Nkg7" in meta_cells else None,
        "design": "FACS tdTomato+/CD45−/CD11b−/TER119−/CD31− Smart-seq2",
        "tnk_fraction_scored": False,
        "reason": "Immune cells were FACS-excluded. Residual Cd3d/Nkg7/Ptprc is leak, not a T/NK compartment.",
    }

    def agg_units(df: pd.DataFrame, unit_col: str) -> pd.DataFrame:
        rows = []
        for uid, sub in df.groupby(unit_col, sort=True):
            rec = {
                "unit": uid,
                "unit_col": unit_col,
                "animal": sub["animal"].iloc[0],
                "mouse_id": sub["mouseID"].nunique(),
                "genotype": sub["genotype"].iloc[0],
                "week": sub["week"].iloc[0],
                "timesimple": ",".join(sorted(sub["timesimple"].astype(str).unique())),
                "n_cells": int(len(sub)),
                "n_tumors": int(sub["mouseID"].nunique()) if unit_col == "animal" else 1,
                "Cldn4_mean": float(sub["Cldn4"].mean()),
                "Cldn4_pct_pos": float(sub["Cldn4_pos"].mean()),
                "IFN_mean": float(sub["IFN"].mean()),
                "MHC-I/APM_mean": float(sub["MHC-I/APM"].mean()),
                "TJ_mean": float(sub["TJ"].mean()),
                "Epcam_mean": float(sub["Epcam"].mean()) if "Epcam" in sub else np.nan,
                "Ptprc_pct_pos": float((sub["Ptprc"] > 0).mean()) if "Ptprc" in sub else np.nan,
                "any_tnk_marker_pct": float((sub[tnk_present] > 0).any(axis=1).mean()),
            }
            rows.append(rec)
        return pd.DataFrame(rows)

    animals = agg_units(meta_cells, "animal")
    mouse_ids = agg_units(meta_cells, "mouseID")
    animals.to_csv(TABLES / "mouse_units.tsv", sep="\t", index=False)
    mouse_ids.to_csv(TABLES / "deposited_mouseID_units.tsv", sep="\t", index=False)

    cohorts = {
        "KP_n20": animals[cohort_mask(animals, "KP_n20")].copy(),
        "KplusKP_n20": animals[cohort_mask(animals, "KplusKP_n20")].copy(),
        "all_n20": animals[cohort_mask(animals, "all_n20")].copy(),
        "KP_all": animals[cohort_mask(animals, "KP_all")].copy(),
        "deposited_mouseID_KP30w_tumors": mouse_ids[mouse_ids["unit"].astype(str).str.startswith("KP_30w")].copy(),
    }

    family_rows = []
    for cname, d in cohorts.items():
        for fam, col in [("IFN", "IFN_mean"), ("MHC-I/APM", "MHC-I/APM_mean"), ("TJ", "TJ_mean")]:
            sp = spearman(d["Cldn4_mean"], d[col])
            q = q4_vs_q1(d["Cldn4_mean"], d[col])
            family_rows.append(
                {
                    "cohort": cname,
                    "family": fam,
                    "n_units": int(len(d)),
                    "n_cells": int(d["n_cells"].sum()) if len(d) else 0,
                    "spearman_rho": sp["rho"],
                    "spearman_p": sp["p"],
                    "spearman_usable": sp["usable"],
                    **{f"q4_{k}": v for k, v in q.items()},
                    "ifn_mhc_down": (
                        (q["usable"] and q["delta_median"] < 0)
                        or (sp["usable"] and sp["rho"] < 0)
                    )
                    if fam != "TJ"
                    else None,
                }
            )
    fam_df = pd.DataFrame(family_rows)
    fam_df.to_csv(TABLES / "family_q4q1.tsv", sep="\t", index=False)

    # Honest n table
    honest = pd.DataFrame(
        [
            {"item": "GEO cells (QC table)", "n": 3891, "note": "do not quote as analysis n"},
            {"item": "deposited mouseID values", "n": int(smp["mouseID"].nunique()), "note": "paper '39 mice'; KP 30w tumors split"},
            {"item": "biological animals (strip _T#)", "n": int(animals.shape[0]), "note": "true mouse unit"},
            {"item": "T / normal AT2 animals", "n": int((animals["genotype"] == "T").sum()), "note": "not KP"},
            {"item": "K-only animals", "n": int((animals["genotype"] == "K").sum()), "note": "Kras; Trp53 WT"},
            {"item": "KP animals", "n": int((animals["genotype"] == "KP").sum()), "note": "named GEMM"},
            {
                "item": "KP animals with ≥20 cells (PRIMARY)",
                "n": int(len(cohorts["KP_n20"])),
                "note": "Cldn4 vs IFN/MHC/TJ unit",
            },
            {
                "item": "K+KP animals with ≥20 cells",
                "n": int(len(cohorts["KplusKP_n20"])),
                "note": "sensitivity",
            },
            {
                "item": "T/NK fraction mice",
                "n": 0,
                "note": "FACS CD45−; no immune compartment",
            },
        ]
    )
    honest.to_csv(TABLES / "honest_n.tsv", sep="\t", index=False)

    # Gene-set coverage
    cov = []
    a8 = json.loads((DATA / "a8_families.json").read_text())
    ifn_h = set(a8["sets"]["HALLMARK_INTERFERON_GAMMA_RESPONSE"]) | set(
        a8["sets"]["HALLMARK_INTERFERON_ALPHA_RESPONSE"]
    )
    cov.append({"set": "IFN (Hallmark IFNα∪IFNγ → mouse)", "n_human": len(ifn_h), "n_mouse_present": len(families["IFN"])})
    cov.append({"set": "MHC-I/APM (curated mouse)", "n_human": 21, "n_mouse_present": len(families["MHC-I/APM"])})
    tj_h = set(a8["sets"]["KEGG_TIGHT_JUNCTION"]) | set(a8["sets"]["GOBP_TIGHT_JUNCTION_ORGANIZATION"])
    cov.append({"set": "TJ (KEGG∪GOBP, Cldn4 held out)", "n_human": len(tj_h), "n_mouse_present": len(families["TJ"])})
    pd.DataFrame(cov).to_csv(TABLES / "geneset_coverage.tsv", sep="\t", index=False)
    json.dump(families, (TABLES / "mouse_families_used.json").open("w"), indent=2)

    # Cluster audit (all should be epithelial states)
    clus = (
        meta_cells.groupby("clusterK12")
        .agg(
            n_cells=("sampleID", "size"),
            Cldn4_mean=("Cldn4", "mean"),
            Cldn4_pct_pos=("Cldn4_pos", "mean"),
            Epcam_mean=("Epcam", "mean") if "Epcam" in meta_cells else ("Cldn4", "mean"),
            Ptprc_pct=("Ptprc", lambda s: float((s > 0).mean())) if "Ptprc" in meta_cells else ("Cldn4", "size"),
        )
        .reset_index()
    )
    clus.to_csv(TABLES / "cluster_audit.tsv", sep="\t", index=False)

    primary = cohorts["KP_n20"]
    primary.to_csv(TABLES / "primary_KP_n20.tsv", sep="\t", index=False)

    # Figures
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.4), constrained_layout=True)
    for ax, fam, col in zip(
        axes, ["IFN", "MHC-I/APM", "TJ (Cldn4 held out)"], ["IFN_mean", "MHC-I/APM_mean", "TJ_mean"]
    ):
        d = primary
        ax.scatter(d["Cldn4_mean"], d[col], c="#1f4e79", s=36, alpha=0.85, edgecolors="none")
        if len(d) >= 2:
            z = np.polyfit(d["Cldn4_mean"], d[col], 1)
            xs = np.linspace(d["Cldn4_mean"].min(), d["Cldn4_mean"].max(), 50)
            ax.plot(xs, np.polyval(z, xs), color="#888", lw=1)
        sp = spearman(d["Cldn4_mean"], d[col])
        ax.set_title(f"{fam}\nρ={_fmt(sp['rho'])} p={_fmt(sp['p'])} n={sp['n']}")
        ax.set_xlabel("KP mouse mean Cldn4")
        ax.set_ylabel(f"{fam} mean")
    fig.savefig(FIGS / "kp_cldn4_vs_ifn_mhc_tj.png", dpi=150)
    fig.savefig(FIGS / "kp_cldn4_vs_ifn_mhc_tj.pdf")
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.4), constrained_layout=True)
    ranks = primary["Cldn4_mean"].rank(method="average")
    qs = pd.qcut(ranks, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")
    for ax, fam, col in zip(axes, ["IFN", "MHC-I/APM", "TJ"], ["IFN_mean", "MHC-I/APM_mean", "TJ_mean"]):
        q1 = primary.loc[qs == "Q1", col]
        q4 = primary.loc[qs == "Q4", col]
        ax.boxplot([q1.values, q4.values], widths=0.55)
        ax.set_xticklabels(["Cldn4 Q1", "Cldn4 Q4"])
        ax.scatter(np.repeat(1, len(q1)), q1, c="#4c78a8", s=22, zorder=3)
        ax.scatter(np.repeat(2, len(q4)), q4, c="#c44e52", s=22, zorder=3)
        q = q4_vs_q1(primary["Cldn4_mean"], primary[col])
        ax.set_title(f"{fam} Δ={_fmt(q['delta_median'])} p={_fmt(q['p'])}")
        ax.set_ylabel(f"{fam} mean")
    fig.savefig(FIGS / "kp_q4q1_ifn_mhc_tj.png", dpi=150)
    fig.savefig(FIGS / "kp_q4q1_ifn_mhc_tj.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 3.6), constrained_layout=True)
    labels = honest["item"].tolist()
    vals = honest["n"].tolist()
    ax.barh(range(len(labels)), vals, color="#1f4e79")
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("n")
    ax.set_title("Honest n (do not quote 3891 cells)")
    for i, v in enumerate(vals):
        ax.text(v + 0.3, i, str(v), va="center", fontsize=8)
    fig.savefig(FIGS / "honest_n.png", dpi=150)
    fig.savefig(FIGS / "honest_n.pdf")
    plt.close(fig)

    summary = {
        "accession": "GSE154989",
        "species": "Mus musculus",
        "model": "KP GEMM lung (also K and T AT2)",
        "assay": "Smart-seq2 plate; FACS epithelial/tumor",
        "matrix": FILES["h5"],
        "tnk_fraction": "NO-GO",
        "tnk_reason": leak["reason"],
        "immune_leak": leak,
        "primary_cohort": "KP animals with ≥20 epithelial cells",
        "primary_n": int(len(primary)),
        "families_n_genes": {k: len(v) for k, v in families.items()},
        "primary_family": fam_df[fam_df["cohort"] == "KP_n20"].to_dict(orient="records"),
        "sensitivity_KplusKP": fam_df[fam_df["cohort"] == "KplusKP_n20"].to_dict(orient="records"),
        "join_human_concordant4": False,
        "join_reason": "T/NK fraction structurally absent; mouse additive only. IFN/MHC is the remaining arm.",
    }
    json.dump(summary, (TABLES / "summary.json").open("w"), indent=2)
    print(json.dumps(summary, indent=2))
    print("WROTE", TABLES)
    print("WROTE", FIGS)


if __name__ == "__main__":
    main()
