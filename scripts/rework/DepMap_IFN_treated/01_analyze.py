#!/usr/bin/env python3
"""Conceptual rework after the unstimulated RNA mismatch.

Prior hunt (results/hunt_depmap_ifn): NSCLC TACSTD2 vs Hallmark IFN-α RNA
was Spearman ρ = +0.377 (positive), not negative. That test is closed.

This slice does **not** recompute unstimulated Hallmark IFN RNA. It asks
the four questions that test still left open:

1. PRISM IFN-pathway drugs (not recombinant IFN — none exist in PRISM).
2. CRISPR Chronos gene effect of IFN / APM genes vs TACSTD2 RNA.
3. Protein-level ISG / IFN scores (CCLE MS proteomics).
4. MHC-I **protein** (HLA-A/B/C, B2M), not RNA.

Primary cohort: DepMap 24Q4 NSCLC cell lines (same type rules as the
unstimulated hunt). SCLC is sensitivity only.
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

NSCLC_TYPES = {
    "LUAD", "LUSC", "LCLC", "NSCLC", "LUAS", "GCLC", "SMARCA4-UT",
    "NUTCL", "LUMEC", "NSCLCPD",
}
SCLC_TYPES = {"SCLC"}
EXCLUDE_FROM_MALIGNANT = {"ZIMMEPCL", "ZIMMLUNG", "ZIMMMPLC"}

# CRISPR / protein gene panels
IFN_TYPE1 = ["IFNAR1", "IFNAR2", "TYK2", "JAK1", "STAT1", "STAT2", "IRF9"]
IFN_TYPE2 = ["IFNGR1", "IFNGR2", "JAK2", "STAT1", "IRF1"]
IFN_NEG = ["SOCS1", "PTPN2", "ADAR", "USP18"]
APM = ["B2M", "TAP1", "TAP2", "TAPBP", "NLRC5", "HLA-A", "HLA-B", "HLA-C"]
MHC1 = ["HLA-A", "HLA-B", "HLA-C", "B2M"]
PROTEIN_ISG = [
    "STAT1", "STAT2", "IRF9", "IRF1", "MX1", "ISG15", "IFIT1", "IFIT3",
    "OAS1", "OAS2", "OAS3", "EIF2AK2", "IFI35", "BST2", "SAMHD1",
    "TRIM25", "ADAR", "HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2",
    "PSMB8", "PSMB9", "PSMB10",
]
CRISPR_EXTRA = ["TACSTD2", "CLDN4", "CD274", "JAK1", "JAK2", "TYK2"]

# PRISM compound classes (Drug.Name, uppercased in matcher)
IFN_INDUCERS = {
    "IMIQUIMOD", "TILORONE", "BROPIRIMINE", "MIW-815",
    "PIDOTIMOD", "PROPAGERMANIUM",
}
JAK_CORE = {
    "RUXOLITINIB", "TOFACITINIB", "BARICITINIB", "FEDRATINIB",
    "ITACITINIB", "UPADACITINIB", "AZD1480", "CYT387",
    "PEFICITINIB", "NVP-BSK805", "TG-101348", "WHI-P154",
    "CEP-33779", "LY2784544", "AZ960", "XL019", "NS-018",
    "SOLCITINIB", "DECERNOTINIB", "BMS-911543", "DELGOCITINIB",
}

EXPR_GENES = ["TACSTD2", "CLDN4"]


def subtype_group(depmap_type: str) -> str:
    if depmap_type == "LUAD":
        return "LUAD"
    if depmap_type == "LUSC":
        return "LUSC"
    if depmap_type in SCLC_TYPES:
        return "SCLC"
    if depmap_type == "LUCA":
        return "carcinoid"
    if depmap_type in EXCLUDE_FROM_MALIGNANT:
        return "noncancer"
    if depmap_type in NSCLC_TYPES:
        return "other_NSCLC"
    return "other"


def symbol_map(columns: list[str]) -> dict[str, str]:
    m: dict[str, str] = {}
    for c in columns:
        if not c:
            continue
        sym = c.split(" (")[0].strip()
        if sym and sym not in m:
            m[sym] = c
    return m


def load_wide_symbols(path: Path, wanted: list[str]) -> tuple[pd.DataFrame, dict[str, str]]:
    header = pd.read_csv(path, nrows=0)
    cols = list(header.columns)
    smap = symbol_map(cols)
    id_col = cols[0]
    use = [id_col]
    found: dict[str, str] = {}
    for s in wanted:
        if s in smap and smap[s] not in use:
            use.append(smap[s])
            found[s] = smap[s]
    df = pd.read_csv(path, usecols=use)
    df = df.rename(columns={id_col: "ModelID"})
    df = df.rename(columns={v: k for k, v in found.items()})
    return df.set_index("ModelID"), found


def zscore(s: pd.Series) -> pd.Series:
    mu = s.mean()
    sd = s.std(ddof=0)
    if sd == 0 or not np.isfinite(sd):
        return s * 0.0
    return (s - mu) / sd


def signature(df: pd.DataFrame, cols: list[str]) -> pd.Series:
    use = [c for c in cols if c in df.columns]
    if not use:
        return pd.Series(np.nan, index=df.index)
    return df[use].apply(zscore, axis=0).mean(axis=1)


def spearman(x: pd.Series, y: pd.Series) -> dict:
    a = pd.concat([x, y], axis=1).dropna()
    n = len(a)
    if n < 8:
        return {"n": n, "rho": np.nan, "p": np.nan}
    rho, p = stats.spearmanr(a.iloc[:, 0], a.iloc[:, 1])
    return {"n": int(n), "rho": float(rho), "p": float(p)}


def mw(a: pd.Series, b: pd.Series) -> dict:
    a = a.dropna()
    b = b.dropna()
    if len(a) < 5 or len(b) < 5:
        return {
            "n_a": int(len(a)), "n_b": int(len(b)),
            "median_a": np.nan, "median_b": np.nan, "p": np.nan, "rbc": np.nan,
        }
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    rbc = 1 - (2 * u) / (len(a) * len(b))
    return {
        "n_a": int(len(a)), "n_b": int(len(b)),
        "median_a": float(a.median()), "median_b": float(b.median()),
        "p": float(p), "rbc": float(rbc),
    }


def bh(pvals: list[float]) -> list[float]:
    arr = np.asarray(pvals, dtype=float)
    out = np.full(arr.shape, np.nan)
    ok = np.isfinite(arr)
    if ok.sum() == 0:
        return out.tolist()
    _, q, _, _ = multipletests(arr[ok], method="fdr_bh")
    out[ok] = q
    return out.tolist()


def parse_proteomics(path: Path) -> tuple[pd.DataFrame, list[str]]:
    """Return proteins x ModelID matrix (gene symbols as index, mean if dup)."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        header = f.readline().rstrip("\n").split(",")
    # find gene-symbol column and ACH- / CCLE sample columns
    gene_col = None
    for cand in ("Gene_Symbol", "gene_symbol", "GeneSymbol", "Gene", "gene"):
        if cand in header:
            gene_col = cand
            break
    if gene_col is None:
        # Nusinow file often has Gene_Symbol as col 1
        for i, h in enumerate(header):
            if "gene" in h.lower() and "symbol" in h.lower():
                gene_col = h
                break
    if gene_col is None:
        raise RuntimeError(f"no gene symbol column in proteomics header: {header[:20]}")

    sample_cols = [c for c in header if c.startswith("ACH-")]
    ccle_like = [c for c in header if c.endswith("_TenTen") or "_BREAST" in c or c.count("_") >= 1]
    # Prefer ACH- IDs. If absent, keep columns that look like CCLE names (contain underscore + tissue).
    if not sample_cols:
        skip = {
            gene_col, "Protein_Id", "Description", "Group_ID", "Uniprot",
            "First.TenTen", "Protein", "Protein.Id",
        }
        sample_cols = [c for c in header if c not in skip and c and "TenTen" not in c and not c.startswith("Uniprot")]
        # too greedy — restrict to columns with a tissue suffix
        sample_cols = [c for c in sample_cols if any(x in c for x in (
            "_LUNG", "_BREAST", "_SKIN", "_OVARY", "_KIDNEY", "_COLON",
            "_PANCREAS", "_LIVER", "_STOMACH", "_PROSTATE", "_THYROID",
            "_ESOPHAGUS", "_SOFT_TISSUE", "_HAEMATOPOIETIC", "_CENTRAL",
            "_UPPER_AERODIGESTIVE", "_URINARY", "_ENDOMETRIUM", "_BONE",
            "_PLEURA", "_AUTONOMIC", "_BILIARY", "_CERVIX",
        )) or c.startswith("ACH-")]

    usecols = [gene_col] + sample_cols
    df = pd.read_csv(path, usecols=lambda c: c in set(usecols))
    df = df.rename(columns={gene_col: "gene"})
    df["gene"] = df["gene"].astype(str).str.split(" ").str[0]
    df = df[df["gene"].isin(set(PROTEIN_ISG + MHC1 + ["TACSTD2", "CLDN4"]))]
    # mean duplicate genes
    num = df.drop(columns=["gene"]).apply(pd.to_numeric, errors="coerce")
    num["gene"] = df["gene"].values
    prot = num.groupby("gene").mean(numeric_only=True)
    return prot, list(prot.columns)


def map_proteomics_columns(sample_cols: list[str], model: pd.DataFrame) -> dict[str, str]:
    """Map proteomics sample column -> ModelID."""
    out: dict[str, str] = {}
    ccle = model["CCLEName"].dropna()
    ccle_map = {str(v): i for i, v in ccle.items()}
    stripped = model["StrippedCellLineName"].dropna()
    stripped_map = {str(v).upper(): i for i, v in stripped.items()}
    for c in sample_cols:
        if c.startswith("ACH-"):
            out[c] = c.split("_")[0]
            continue
        if c in ccle_map:
            out[c] = ccle_map[c]
            continue
        # strip TenTen / extra suffixes
        base = c.replace("_TenTen", "")
        if base in ccle_map:
            out[c] = ccle_map[base]
            continue
        token = base.split("_")[0].upper()
        if token in stripped_map:
            out[c] = stripped_map[token]
    return out


def forest_rows(rows: list[dict], title: str, out: Path, xlim=None) -> None:
    if not rows:
        return
    fig, ax = plt.subplots(figsize=(7.2, max(2.8, 0.32 * len(rows) + 1.2)))
    y = np.arange(len(rows))
    rhos = [r["rho"] for r in rows]
    ax.axvline(0, color="#444", lw=0.8)
    ax.scatter(rhos, y, s=28, c="#1f4e79", zorder=3)
    labels = []
    for r in rows:
        star = ""
        q = r.get("q", r.get("p"))
        if q is not None and np.isfinite(q) and q < 0.05:
            star = " *"
        labels.append(f"{r['label']}  n={r['n']}{star}")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Spearman ρ vs TACSTD2 RNA")
    ax.set_title(title, fontsize=10)
    if xlim:
        ax.set_xlim(*xlim)
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


def scatter_plot(x, y, xlabel, ylabel, title, out, hue=None) -> None:
    d = pd.concat([x.rename("x"), y.rename("y")], axis=1)
    if hue is not None:
        d = pd.concat([d, hue.rename("hue")], axis=1)
    d = d.dropna(subset=["x", "y"])
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    if "hue" in d.columns:
        for k, g in d.groupby("hue"):
            ax.scatter(g["x"], g["y"], s=16, alpha=0.75, label=str(k))
        ax.legend(fontsize=7, frameon=False)
    else:
        ax.scatter(d["x"], d["y"], s=16, alpha=0.7, c="#1f4e79")
    if len(d) >= 8:
        rho, p = stats.spearmanr(d["x"], d["y"])
        ax.set_title(f"{title}\nρ={rho:.3f}  p={p:.2e}  n={len(d)}", fontsize=9)
    else:
        ax.set_title(f"{title}\nn={len(d)}", fontsize=9)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data/rework/DepMap_IFN_treated")
    ap.add_argument("--outdir", default="results/rework/DepMap_IFN_treated")
    args = ap.parse_args()
    data = Path(args.data)
    out = Path(args.outdir)
    figdir = out / "figures"
    tabdir = out / "tables"
    figdir.mkdir(parents=True, exist_ok=True)
    tabdir.mkdir(parents=True, exist_ok=True)

    model = pd.read_csv(data / "Model.csv").set_index("ModelID")
    model["group"] = model["DepmapModelType"].astype(str).map(subtype_group)
    lung = model[model["OncotreeLineage"] == "Lung"].copy()
    nsclc = lung[lung["group"].isin(["LUAD", "LUSC", "other_NSCLC"])].copy()
    sclc = lung[lung["group"] == "SCLC"].copy()

    expr, expr_found = load_wide_symbols(
        data / "OmicsExpressionProteinCodingGenesTPMLogp1.csv",
        EXPR_GENES,
    )
    crispr_wanted = sorted(set(
        IFN_TYPE1 + IFN_TYPE2 + IFN_NEG + APM + CRISPR_EXTRA
    ))
    crispr, crispr_found = load_wide_symbols(
        data / "CRISPRGeneEffect.csv",
        crispr_wanted,
    )

    # --- protein ---
    prot_path = data / "protein_quant_current_normalized.csv.gz"
    protein_ok = prot_path.exists() and prot_path.stat().st_size > 1_000_000
    prot_long = None
    prot_found: list[str] = []
    if protein_ok:
        try:
            prot_mat, prot_cols = parse_proteomics(prot_path)
            colmap = map_proteomics_columns(prot_cols, model)
            keep = [c for c in prot_mat.columns if c in colmap]
            prot_mat = prot_mat[keep].rename(columns=colmap)
            # if duplicate ModelIDs, mean
            prot_mat = prot_mat.T.groupby(level=0).mean().T
            prot_long = prot_mat.T  # ModelID x gene
            prot_found = [g for g in prot_long.columns]
        except Exception as e:  # noqa: BLE001
            protein_ok = False
            prot_err = str(e)
    else:
        prot_err = "proteomics file missing or too small"

    # --- PRISM ---
    prism_mat = pd.read_csv(
        data / "Repurposing_Public_23Q2_Extended_Primary_Data_Matrix.csv",
        index_col=0,
    )
    compounds = pd.read_csv(
        data / "Repurposing_Public_23Q2_Extended_Primary_Compound_List.csv"
    )
    # row index of matrix should match Drug.Name or IDs
    prism_mat.index = prism_mat.index.astype(str)
    compounds["Drug.Name"] = compounds["Drug.Name"].astype(str)
    compounds["name_up"] = compounds["Drug.Name"].str.upper()
    # try to align matrix rows to drug names
    if set(prism_mat.index) & set(compounds["Drug.Name"]):
        name_of = {n: n for n in prism_mat.index}
    elif set(prism_mat.index) & set(compounds["IDs"].astype(str)):
        name_of = dict(zip(compounds["IDs"].astype(str), compounds["Drug.Name"]))
    else:
        # first column may already be names with extra quotes stripped
        name_of = {i: i for i in prism_mat.index}

    def classify_drug(name: str) -> str:
        u = str(name).upper()
        if u in IFN_INDUCERS:
            return "IFN_inducer_or_IFNR_agonist"
        if u in JAK_CORE:
            return "JAK_inhibitor_core"
        return "other"

    # collapse duplicate drug names (multiple IDs / screens) by median LFC
    prism_named = prism_mat.copy()
    prism_named["_name"] = [name_of.get(i, i) for i in prism_named.index]
    prism_named["_class"] = prism_named["_name"].map(classify_drug)
    # keep a per-row table for background
    lfc = prism_named.drop(columns=["_name", "_class"])

    # NSCLC models with TACSTD2
    nsclc_ids = nsclc.index.intersection(expr.index)
    sample = expr.loc[nsclc_ids, ["TACSTD2"]].copy()
    if "CLDN4" in expr.columns:
        sample["CLDN4"] = expr.loc[nsclc_ids, "CLDN4"]
    sample["group"] = nsclc.loc[sample.index, "group"]
    sample["DepmapModelType"] = nsclc.loc[sample.index, "DepmapModelType"]
    sample["CellLineName"] = nsclc.loc[sample.index, "CellLineName"]

    sclc_ids = sclc.index.intersection(expr.index)
    sclc_tac = expr.loc[sclc_ids, "TACSTD2"] if len(sclc_ids) else pd.Series(dtype=float)

    # ========== 1. PRISM ==========
    prism_rows = []
    # restrict to cell lines present in both
    prism_cells = [c for c in lfc.columns if c in sample.index]
    tac = sample.loc[prism_cells, "TACSTD2"]
    for idx, row in prism_named.iterrows():
        name = row["_name"]
        klass = row["_class"]
        y = pd.to_numeric(row[prism_cells], errors="coerce")
        y.index = prism_cells
        st = spearman(tac, y)
        st.update({
            "compound_id": str(idx),
            "name": str(name),
            "class": klass,
            "screen_guess": "PRISM_23Q2_extended_primary",
        })
        prism_rows.append(st)
    prism_df = pd.DataFrame(prism_rows)
    # collapse to unique drug name (median rho if multiple rows)
    def collapse_prism(g: pd.DataFrame) -> pd.Series:
        return pd.Series({
            "n_rows": len(g),
            "n": int(g["n"].median()) if g["n"].notna().any() else 0,
            "rho": float(g["rho"].median()) if g["rho"].notna().any() else np.nan,
            "p": float(g["p"].median()) if g["p"].notna().any() else np.nan,
            "class": g["class"].iloc[0],
        })

    parts = []
    for name, g in prism_df.dropna(subset=["rho"]).groupby("name"):
        rec = collapse_prism(g)
        rec["name"] = name
        parts.append(rec)
    prism_by_name = pd.DataFrame(parts) if parts else pd.DataFrame(
        columns=["name", "n_rows", "n", "rho", "p", "class"]
    )
    # FDR within all compounds (exploratory background)
    if len(prism_by_name):
        prism_by_name["q_all"] = bh(prism_by_name["p"].tolist())
    prism_by_name = prism_by_name.sort_values("rho")
    prism_by_name.to_csv(tabdir / "prism_all_compounds_vs_TACSTD2.tsv", sep="\t", index=False)

    focus = prism_by_name[prism_by_name["class"] != "other"].copy()
    if len(focus):
        focus["q_focus"] = bh(focus["p"].tolist())
    focus.to_csv(tabdir / "prism_IFN_pathway_compounds.tsv", sep="\t", index=False)

    # background percentile of focus drugs
    bg = prism_by_name["rho"].dropna()
    focus_summary = []
    for _, r in focus.iterrows():
        if bg.size and np.isfinite(r["rho"]):
            pct = float((bg < r["rho"]).mean())
        else:
            pct = np.nan
        focus_summary.append({**r.to_dict(), "background_percentile": pct})

    # class-level: mean LFC of JAK core / IFN inducers vs TACSTD2
    def class_lfc(names: set[str]) -> pd.Series:
        mask = prism_named["_name"].str.upper().isin(names)
        sub = lfc.loc[mask, prism_cells]
        if sub.empty:
            return pd.Series(dtype=float)
        return sub.apply(pd.to_numeric, errors="coerce").median(axis=0)

    jak_lfc = class_lfc(JAK_CORE)
    ifn_lfc = class_lfc(IFN_INDUCERS)
    class_stats = {
        "JAK_inhibitor_core_medianLFC": spearman(tac.reindex(jak_lfc.index), jak_lfc),
        "IFN_inducer_or_IFNR_agonist_medianLFC": spearman(tac.reindex(ifn_lfc.index), ifn_lfc),
    }

    # ========== 2. CRISPR ==========
    crispr_nsclc = crispr.reindex(sample.index)
    crispr_rows = []
    for g in crispr_wanted:
        if g not in crispr_nsclc.columns:
            continue
        panel = (
            "IFN_type1" if g in IFN_TYPE1 else
            "IFN_type2" if g in IFN_TYPE2 else
            "IFN_negative_reg" if g in IFN_NEG else
            "APM_MHC" if g in APM else
            "other"
        )
        st = spearman(sample["TACSTD2"], crispr_nsclc[g])
        st.update({"gene": g, "panel": panel, "present": True})
        crispr_rows.append(st)
    crispr_df = pd.DataFrame(crispr_rows)
    if len(crispr_df):
        # primary FDR: IFN type I/II + APM (not TACSTD2/CLDN4 self)
        prim_mask = crispr_df["panel"].isin(["IFN_type1", "IFN_type2", "APM_MHC", "IFN_negative_reg"])
        q = np.full(len(crispr_df), np.nan)
        if prim_mask.any():
            q[prim_mask.to_numpy()] = bh(crispr_df.loc[prim_mask, "p"].tolist())
        crispr_df["q_primary"] = q
    crispr_df.to_csv(tabdir / "crispr_gene_effect_vs_TACSTD2.tsv", sep="\t", index=False)

    # signature of IFN-sensing Chronos (mean of type I receptors + JAK1/STAT1/IRF9)
    type1_cols = [g for g in IFN_TYPE1 if g in crispr_nsclc.columns]
    apm_cols = [g for g in APM if g in crispr_nsclc.columns]
    crispr_nsclc = crispr_nsclc.copy()
    crispr_nsclc["sig_IFN_type1_effect"] = signature(crispr_nsclc, type1_cols)
    crispr_nsclc["sig_APM_effect"] = signature(crispr_nsclc, apm_cols)
    crispr_sig = {
        "IFN_type1_Chronos_mean": spearman(sample["TACSTD2"], crispr_nsclc["sig_IFN_type1_effect"]),
        "APM_Chronos_mean": spearman(sample["TACSTD2"], crispr_nsclc["sig_APM_effect"]),
    }
    if "TACSTD2" in crispr_nsclc.columns:
        crispr_sig["TACSTD2_dependency"] = {
            "n": int(crispr_nsclc["TACSTD2"].notna().sum()),
            "median_Chronos": float(crispr_nsclc["TACSTD2"].median()),
            "note": "Chronos ~0 means TACSTD2 is not a lung-lineage essential",
        }

    # SCLC sensitivity for STAT1 / IFNAR1
    sclc_crispr = {}
    if len(sclc_tac):
        for g in ["STAT1", "IFNAR1", "IFNGR1", "B2M", "JAK1"]:
            if g in crispr.columns:
                sclc_crispr[g] = spearman(sclc_tac, crispr.reindex(sclc_tac.index)[g])

    # ========== 3+4. Protein ==========
    protein_rows = []
    protein_sig = {}
    n_prot_nsclc = 0
    if protein_ok and prot_long is not None:
        prot_nsclc = prot_long.reindex(sample.index)
        n_prot_nsclc = int(prot_nsclc.dropna(how="all").shape[0])
        for g in prot_found:
            st = spearman(sample["TACSTD2"], prot_nsclc[g])
            st.update({
                "protein": g,
                "panel": (
                    "MHC1" if g in MHC1 else
                    "ISG" if g in PROTEIN_ISG else
                    "self"
                ),
            })
            protein_rows.append(st)
        prot_df = pd.DataFrame(protein_rows)
        if len(prot_df):
            mhc_isg = prot_df["panel"].isin(["MHC1", "ISG"])
            q = np.full(len(prot_df), np.nan)
            if mhc_isg.any():
                q[mhc_isg.to_numpy()] = bh(prot_df.loc[mhc_isg, "p"].tolist())
            prot_df["q"] = q
        prot_df.to_csv(tabdir / "protein_vs_TACSTD2_RNA.tsv", sep="\t", index=False)

        isg_cols = [g for g in PROTEIN_ISG if g in prot_nsclc.columns]
        mhc_cols = [g for g in MHC1 if g in prot_nsclc.columns]
        prot_nsclc = prot_nsclc.copy()
        prot_nsclc["sig_ISG_protein"] = signature(prot_nsclc, isg_cols)
        prot_nsclc["sig_MHC1_protein"] = signature(prot_nsclc, mhc_cols)
        protein_sig = {
            "ISG_protein_score": spearman(sample["TACSTD2"], prot_nsclc["sig_ISG_protein"]),
            "MHC1_protein_score": spearman(sample["TACSTD2"], prot_nsclc["sig_MHC1_protein"]),
            "genes_in_ISG_score": isg_cols,
            "genes_in_MHC1_score": mhc_cols,
            "n_NSCLC_with_any_protein": n_prot_nsclc,
        }
        if "TACSTD2" in prot_nsclc.columns:
            protein_sig["TACSTD2_protein_vs_MHC1_protein"] = spearman(
                prot_nsclc["TACSTD2"], prot_nsclc["sig_MHC1_protein"]
            )
            protein_sig["TACSTD2_protein_vs_ISG_protein"] = spearman(
                prot_nsclc["TACSTD2"], prot_nsclc["sig_ISG_protein"]
            )
            protein_sig["TACSTD2_RNA_vs_protein"] = spearman(
                sample["TACSTD2"], prot_nsclc["TACSTD2"]
            )
        # median-split MHC-I protein
        if prot_nsclc["sig_MHC1_protein"].notna().sum() >= 10:
            med = sample["TACSTD2"].median()
            hi = prot_nsclc.loc[sample["TACSTD2"] >= med, "sig_MHC1_protein"]
            lo = prot_nsclc.loc[sample["TACSTD2"] < med, "sig_MHC1_protein"]
            protein_sig["MHC1_protein_median_split"] = mw(hi, lo)
    else:
        prot_df = pd.DataFrame()
        protein_sig = {"error": prot_err if 'prot_err' in dir() else "unavailable"}

    # ========== figures ==========
    scatter_plot(
        sample["TACSTD2"],
        crispr_nsclc["sig_IFN_type1_effect"] if "sig_IFN_type1_effect" in crispr_nsclc.columns else pd.Series(dtype=float),
        "TACSTD2 RNA  log2(TPM+1)",
        "IFN type-I CRISPR effect (z-mean Chronos)",
        "NSCLC: TACSTD2 vs IFN type-I gene effect",
        figdir / "fig1_tacstd2_vs_IFN_type1_Chronos.png",
        hue=sample["group"],
    )
    if "B2M" in crispr_nsclc.columns:
        scatter_plot(
            sample["TACSTD2"], crispr_nsclc["B2M"],
            "TACSTD2 RNA  log2(TPM+1)", "B2M Chronos (more negative = more dependent)",
            "NSCLC: TACSTD2 vs B2M gene effect",
            figdir / "fig2_tacstd2_vs_B2M_Chronos.png",
            hue=sample["group"],
        )
    if protein_ok and prot_long is not None and "sig_MHC1_protein" in prot_nsclc.columns:
        scatter_plot(
            sample["TACSTD2"], prot_nsclc["sig_MHC1_protein"],
            "TACSTD2 RNA  log2(TPM+1)", "MHC-I protein score (z-mean HLA-A/B/C/B2M)",
            "NSCLC: TACSTD2 RNA vs MHC-I protein",
            figdir / "fig3_tacstd2_vs_MHC1_protein.png",
            hue=sample["group"],
        )
        scatter_plot(
            sample["TACSTD2"], prot_nsclc["sig_ISG_protein"],
            "TACSTD2 RNA  log2(TPM+1)", "ISG protein score (z-mean detected ISGs)",
            "NSCLC: TACSTD2 RNA vs ISG protein",
            figdir / "fig4_tacstd2_vs_ISG_protein.png",
            hue=sample["group"],
        )
    if len(jak_lfc):
        scatter_plot(
            tac.reindex(jak_lfc.index), jak_lfc,
            "TACSTD2 RNA  log2(TPM+1)", "Median LFC across core JAK inhibitors",
            "NSCLC PRISM: TACSTD2 vs JAK-inhibitor viability",
            figdir / "fig5_tacstd2_vs_JAK_PRISM.png",
            hue=sample.reindex(jak_lfc.index)["group"],
        )
    if len(ifn_lfc):
        scatter_plot(
            tac.reindex(ifn_lfc.index), ifn_lfc,
            "TACSTD2 RNA  log2(TPM+1)", "Median LFC across IFN inducers / IFNR agonists",
            "NSCLC PRISM: TACSTD2 vs IFN-inducer viability",
            figdir / "fig6_tacstd2_vs_IFNinducer_PRISM.png",
            hue=sample.reindex(ifn_lfc.index)["group"],
        )

    # forests
    forest_rows(
        [
            {"label": r["gene"], "n": r["n"], "rho": r["rho"], "q": r.get("q_primary"), "p": r["p"]}
            for _, r in crispr_df.sort_values("rho").iterrows()
            if r.get("panel") in {"IFN_type1", "IFN_type2", "APM_MHC", "IFN_negative_reg"}
            and np.isfinite(r.get("rho", np.nan))
        ],
        "CRISPR Chronos vs TACSTD2 RNA (NSCLC)",
        figdir / "fig7_crispr_forest.png",
        xlim=(-0.55, 0.55),
    )
    if len(prot_df):
        forest_rows(
            [
                {"label": r["protein"], "n": r["n"], "rho": r["rho"], "q": r.get("q"), "p": r["p"]}
                for _, r in prot_df.sort_values("rho").iterrows()
                if r.get("panel") in {"MHC1", "ISG"} and np.isfinite(r.get("rho", np.nan))
            ],
            "Protein abundance vs TACSTD2 RNA (NSCLC, CCLE MS)",
            figdir / "fig8_protein_forest.png",
            xlim=(-0.7, 0.7),
        )

    # sample table
    sample_out = sample.copy()
    if "sig_IFN_type1_effect" in crispr_nsclc.columns:
        sample_out["sig_IFN_type1_Chronos"] = crispr_nsclc["sig_IFN_type1_effect"]
    if protein_ok and prot_long is not None:
        if "sig_MHC1_protein" in prot_nsclc.columns:
            sample_out["sig_MHC1_protein"] = prot_nsclc["sig_MHC1_protein"]
        if "sig_ISG_protein" in prot_nsclc.columns:
            sample_out["sig_ISG_protein"] = prot_nsclc["sig_ISG_protein"]
        if "TACSTD2" in prot_nsclc.columns:
            sample_out["TACSTD2_protein"] = prot_nsclc["TACSTD2"]
    if len(jak_lfc):
        sample_out["PRISM_JAK_core_medianLFC"] = jak_lfc
    if len(ifn_lfc):
        sample_out["PRISM_IFNinducer_medianLFC"] = ifn_lfc
    sample_out.to_csv(tabdir / "sample_table.tsv", sep="\t")

    cohort = {
        "n_lung_annot": int((model["OncotreeLineage"] == "Lung").sum()),
        "n_NSCLC_annot": int(len(nsclc)),
        "n_NSCLC_expr": int(len(sample)),
        "n_NSCLC_by_group": sample["group"].value_counts().to_dict(),
        "n_SCLC_expr": int(len(sclc_tac)),
        "n_NSCLC_PRISM": int(len(prism_cells)),
        "n_NSCLC_CRISPR": int(crispr_nsclc.dropna(how="all").shape[0]),
        "n_NSCLC_protein": n_prot_nsclc,
    }
    pd.Series(cohort).to_csv(tabdir / "cohort_counts.tsv", sep="\t", header=["value"])

    coverage = pd.DataFrame([
        {"layer": "RNA", "item": k, "present": True} for k in expr_found
    ] + [
        {"layer": "CRISPR", "item": k, "present": True} for k in crispr_found
    ] + [
        {"layer": "protein", "item": k, "present": True} for k in prot_found
    ])
    coverage.to_csv(tabdir / "coverage.tsv", sep="\t", index=False)

    # honest missing-data flags
    missing = {
        "recombinant_IFN_treated_DepMap_RNAseq": False,
        "recombinant_IFN_biologic_in_PRISM": False,
        "note": (
            "DepMap 24Q4 expression is unstimulated culture RNA. "
            "PRISM 23Q2 / 19Q4 are small-molecule libraries; interferon alfa/gamma "
            "proteins are not compounds in the extended primary list. "
            "IFN inducers and JAK inhibitors are the public PRISM proxies."
        ),
    }

    key = {
        "slice": "results/rework/DepMap_IFN_treated",
        "why_rework": {
            "prior_slice": "results/hunt_depmap_ifn",
            "unstimulated_NSCLC_TACSTD2_vs_Hallmark_IFNA_rho": 0.377,
            "sign": "positive, not negative",
            "closed": True,
        },
        "release": {
            "depmap": "DepMap Public 24Q4",
            "depmap_doi": "10.25452/figshare.plus.27993248.v1",
            "prism": "PRISM Repurposing Public 23Q2",
            "prism_doi": "10.6084/m9.figshare.23600310.v4",
            "proteomics": "Nusinow et al. Cell 2020 CCLE MS (Gygi)",
        },
        "cohort": cohort,
        "missing_public_data": missing,
        "expr_genes_found": expr_found,
        "crispr_genes_found": list(crispr_found),
        "protein_genes_found": prot_found,
        "prism_class_stats": class_stats,
        "prism_focus": focus_summary,
        "crispr_vs_TACSTD2": crispr_df.to_dict(orient="records") if len(crispr_df) else [],
        "crispr_signatures": crispr_sig,
        "crispr_SCLC_sensitivity": sclc_crispr,
        "protein_vs_TACSTD2": prot_df.to_dict(orient="records") if len(prot_df) else [],
        "protein_signatures": protein_sig,
        "primary_tests": {
            "CRISPR_IFN_type1_mean_vs_TACSTD2": crispr_sig.get("IFN_type1_Chronos_mean"),
            "CRISPR_APM_mean_vs_TACSTD2": crispr_sig.get("APM_Chronos_mean"),
            "PRISM_JAK_core_medianLFC_vs_TACSTD2": class_stats.get("JAK_inhibitor_core_medianLFC"),
            "PRISM_IFNinducer_medianLFC_vs_TACSTD2": class_stats.get("IFN_inducer_or_IFNR_agonist_medianLFC"),
            "protein_ISG_vs_TACSTD2_RNA": protein_sig.get("ISG_protein_score"),
            "protein_MHC1_vs_TACSTD2_RNA": protein_sig.get("MHC1_protein_score"),
        },
    }
    (tabdir / "key_stats.json").write_text(json.dumps(key, indent=2, default=str))
    print(json.dumps({"cohort": cohort, "primary_tests": key["primary_tests"]}, indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
