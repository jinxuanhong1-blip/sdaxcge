#!/usr/bin/env python3
"""CLDN4 vs DNA-damage signatures, STING, and HLA in DepMap lung lines.

γH2AX protein is not in CCLE RPPA or MCLP. DNA-damage readouts are:
  - MSigDB Hallmark DNA repair (mean-z of RNA)
  - a pre-specified DSB-response RNA signature
  - MCLP Valid phospho marks ATM pS1981 and Rad17 pS645
STING is the core machinery (CGAS, STING1, TBK1, IKBKE, IRF3), not Hallmark IFN-γ.
HLA primary score is classical HLA-A/B/C. MHC-I (+B2M) is a replication of PR #304.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

# H2AX transcript is scored on its own and is not part of this signature.
# Transcript abundance is not γH2AX.
DSB_GENES = [
    "ATM",
    "ATR",
    "CHEK1",
    "CHEK2",
    "BRCA1",
    "BRCA2",
    "RAD51",
    "MDC1",
    "TP53BP1",
    "PARP1",
    "PRKDC",
    "XRCC5",
    "XRCC6",
    "RPA1",
    "RPA2",
    "MRE11",
    "RAD50",
    "NBN",
    "RAD17",
]
STING_CORE = ["CGAS", "STING1", "TBK1", "IKBKE", "IRF3"]
STING_ISG = ["CXCL10", "CCL5", "ISG15", "IFIT1", "MX1", "OAS1", "IRF7"]
HLA_ABC = ["HLA-A", "HLA-B", "HLA-C"]
MHC1 = ["HLA-A", "HLA-B", "HLA-C", "B2M"]
HLA_II = ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1"]

# role: primary endpoints enter BH. replication / secondary / control do not.
RNA_TARGETS = [
    ("sig_DNA_REPAIR", "Hallmark DNA repair", "rna_primary"),
    ("sig_DSB", "DSB-response RNA", "rna_primary"),
    ("sig_STING_CORE", "STING core", "rna_primary"),
    ("sig_HLA_ABC", "HLA-A/B/C", "rna_primary"),
    ("HLA-A", "HLA-A", "hla_gene"),
    ("HLA-B", "HLA-B", "hla_gene"),
    ("HLA-C", "HLA-C", "hla_gene"),
    ("sig_MHC1", "MHC-I (HLA-A/B/C+B2M)", "replication"),
    ("sig_IFNG", "Hallmark IFN-γ", "replication"),
    ("sig_STING_ISG", "STING-adjacent ISG", "secondary"),
    ("sig_HLA_II", "HLA class II", "secondary"),
    ("STING1", "STING1", "secondary"),
    ("CGAS", "CGAS", "secondary"),
    ("TBK1", "TBK1", "secondary"),
    ("IKBKE", "IKBKE", "secondary"),
    ("IRF3", "IRF3 RNA", "secondary"),
    ("H2AX", "H2AX transcript (not γH2AX)", "control"),
    ("TACSTD2", "TACSTD2", "companion"),
    ("B2M", "B2M", "secondary"),
    ("HLA-E", "HLA-E", "secondary"),
]

# Valid antibodies are the protein primary family. Caution antibodies are reported
# with their MCLP validation flag and are not in the primary FDR.
PROTEIN_TARGETS = [
    ("ATMPS1981", "ATM pS1981", "protein_primary", "Valid"),
    ("RAD17PS645", "Rad17 pS645", "protein_primary", "Valid"),
    ("CGAS", "cGAS protein", "protein_primary", "Valid"),
    ("IRF3", "IRF3 protein", "protein_primary", "Valid"),
    ("HLADQA1", "HLA-DQA1 protein", "protein_primary", "Valid"),
    ("CHK1PS345", "Chk1 pS345", "caution", "Caution"),
    ("CHK2PT68", "Chk2 pT68", "caution", "Caution"),
    ("HLADRDPDQDX", "HLA-DR/DP/DQ/DX protein", "caution", "Caution"),
    ("PDL1", "PD-L1 protein", "caution", "Caution"),
    ("IRF1", "IRF1 protein", "caution", "Caution"),
    ("ATM", "ATM total", "secondary", "Valid"),
    ("X53BP1", "53BP1 total", "secondary", "Valid"),
    ("RAD51", "RAD51 total", "secondary", "Valid"),
    ("RAD50", "RAD50 total", "secondary", "Valid"),
    ("HISTONEH3PS10", "Histone H3 pS10 (mitotic, not γH2AX)", "control", "Valid"),
    ("CLAUDIN7", "Claudin-7 protein", "companion", "Valid"),
]


def zscore(s: pd.Series) -> pd.Series:
    mu = s.mean()
    sd = s.std(ddof=0)
    if sd == 0 or not np.isfinite(sd):
        return s * 0.0
    return (s - mu) / sd


def signature(df: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, list[str]]:
    use = [g for g in genes if g in df.columns]
    if not use:
        return pd.Series(np.nan, index=df.index), []
    z = df[use].apply(lambda col: zscore(pd.to_numeric(col, errors="coerce")), axis=0)
    return z.mean(axis=1, skipna=True), use


def spearman(x: pd.Series, y: pd.Series) -> dict:
    a = pd.concat([pd.to_numeric(x, errors="coerce"), pd.to_numeric(y, errors="coerce")], axis=1).dropna()
    n = int(len(a))
    if n < 8:
        return {"n": n, "rho": None, "p": None}
    rho, p = stats.spearmanr(a.iloc[:, 0], a.iloc[:, 1])
    return {"n": n, "rho": float(rho), "p": float(p)}


def bootstrap_spearman_ci(x: np.ndarray, y: np.ndarray, n_boot: int = 5000, seed: int = 0) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    n = len(x)
    boots = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        boots[i] = stats.spearmanr(x[idx], y[idx]).statistic
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return float(lo), float(hi)


def bh_fdr(pvals: list[float]) -> list[float]:
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for rank, i in enumerate(order[::-1], start=1):
        k = n - rank + 1
        val = min(prev, p[i] * n / k)
        q[i] = val
        prev = val
    return [float(min(1.0, v)) for v in q]


def cohort_tag(row: pd.Series) -> str:
    disease = str(row.get("OncotreePrimaryDisease") or "")
    subtype = str(row.get("OncotreeSubtype") or "")
    if disease == "Non-Small Cell Lung Cancer":
        if subtype == "Lung Adenocarcinoma":
            return "LUAD"
        if subtype == "Lung Squamous Cell Carcinoma":
            return "LUSC"
        return "other_NSCLC"
    if disease == "Lung Neuroendocrine Tumor" or subtype == "Small Cell Lung Cancer":
        return "SCLC_NET"
    return "other_lung"


def norm_id(value):
    if pd.isna(value):
        return np.nan
    text = str(value).strip()
    if text == "" or text.lower() == "nan":
        return np.nan
    return "".join(ch for ch in text.upper() if ch.isalnum())


def partial_spearman(x: pd.Series, y: pd.Series, groups: pd.Series) -> dict:
    """Partial Spearman: Pearson correlation of rank-residuals after group indicators."""
    a = pd.DataFrame({"x": pd.to_numeric(x, errors="coerce"), "y": pd.to_numeric(y, errors="coerce"), "g": groups}).dropna()
    n = int(len(a))
    if n < 8 or a["g"].nunique() < 2:
        return {"n": n, "rho": None, "p": None, "n_groups": int(a["g"].nunique())}
    rx = stats.rankdata(a["x"].to_numpy())
    ry = stats.rankdata(a["y"].to_numpy())
    dummies = pd.get_dummies(a["g"], drop_first=True)
    design = np.column_stack([np.ones(n), dummies.to_numpy(dtype=float)])
    bx, *_ = np.linalg.lstsq(design, rx, rcond=None)
    by, *_ = np.linalg.lstsq(design, ry, rcond=None)
    rho, p = stats.pearsonr(rx - design @ bx, ry - design @ by)
    return {"n": n, "rho": float(rho), "p": float(p), "n_groups": int(a["g"].nunique())}


def resolve_aliases(columns: list[str]) -> dict[str, str]:
    """Map canonical symbols onto whatever DepMap called them."""
    have = set(columns)
    mapping = {}
    pairs = [("H2AX", "H2AFX"), ("CGAS", "MB21D1"), ("STING1", "TMEM173")]
    for preferred, alt in pairs:
        if preferred in have:
            mapping[preferred] = preferred
        elif alt in have:
            mapping[preferred] = alt
    return mapping


def add_scores(frame: pd.DataFrame, ifng_genes: list[str], dna_genes: list[str], alias: dict[str, str]) -> tuple[pd.DataFrame, dict]:
    out = frame.copy()
    dsb = [alias.get(g, g) for g in DSB_GENES]
    sting = [alias.get(g, g) for g in STING_CORE]
    # Drop duplicate columns if an alias resolved onto a gene already listed.
    dsb = list(dict.fromkeys(dsb))
    sting = list(dict.fromkeys(sting))
    sig_dna, used_dna = signature(out, dna_genes)
    sig_dsb, used_dsb = signature(out, dsb)
    sig_sting, used_sting = signature(out, sting)
    sig_isg, used_isg = signature(out, STING_ISG)
    sig_hla, used_hla = signature(out, HLA_ABC)
    sig_mhc, used_mhc = signature(out, MHC1)
    sig_ii, used_ii = signature(out, HLA_II)
    sig_ifng, used_ifng = signature(out, ifng_genes)
    out["sig_DNA_REPAIR"] = sig_dna
    out["sig_DSB"] = sig_dsb
    out["sig_STING_CORE"] = sig_sting
    out["sig_STING_ISG"] = sig_isg
    out["sig_HLA_ABC"] = sig_hla
    out["sig_MHC1"] = sig_mhc
    out["sig_HLA_II"] = sig_ii
    out["sig_IFNG"] = sig_ifng
    if "H2AX" not in out.columns and alias.get("H2AX") in out.columns:
        out["H2AX"] = out[alias["H2AX"]]
    if "CGAS" not in out.columns and alias.get("CGAS") in out.columns:
        out["CGAS"] = out[alias["CGAS"]]
    if "STING1" not in out.columns and alias.get("STING1") in out.columns:
        out["STING1"] = out[alias["STING1"]]
    cov = {
        "dna_repair_n_used": len(used_dna),
        "dna_repair_n_gmt": len(dna_genes),
        "dna_repair_missing": sorted(set(dna_genes) - set(used_dna)),
        "dsb_used": used_dsb,
        "sting_core_used": used_sting,
        "sting_isg_used": used_isg,
        "hla_abc_used": used_hla,
        "mhc1_used": used_mhc,
        "hla_ii_used": used_ii,
        "ifng_n_used": len(used_ifng),
        "ifng_n_gmt": len(ifng_genes),
        "alias": alias,
    }
    return out, cov


def correlate_table(frame: pd.DataFrame, targets: list[tuple], cohort: str, predictor: str = "CLDN4") -> list[dict]:
    rows = []
    for col, label, role, *rest in targets:
        if col not in frame.columns:
            rows.append({
                "cohort": cohort,
                "predictor": predictor,
                "target": col,
                "label": label,
                "role": role,
                "n": 0,
                "rho": None,
                "p": None,
                "ci95_low": None,
                "ci95_high": None,
                "validation": rest[0] if rest else "",
            })
            continue
        rec = spearman(frame[predictor], frame[col])
        rec.update({
            "cohort": cohort,
            "predictor": predictor,
            "target": col,
            "label": label,
            "role": role,
            "validation": rest[0] if rest else "",
            "ci95_low": None,
            "ci95_high": None,
        })
        if cohort in {"lung_cell_lines", "NSCLC"} and rec["n"] and rec["n"] >= 8 and rec["rho"] is not None:
            a = pd.concat([
                pd.to_numeric(frame[predictor], errors="coerce"),
                pd.to_numeric(frame[col], errors="coerce"),
            ], axis=1).dropna()
            lo, hi = bootstrap_spearman_ci(a.iloc[:, 0].to_numpy(), a.iloc[:, 1].to_numpy())
            rec["ci95_low"] = lo
            rec["ci95_high"] = hi
        rows.append(rec)
    return rows


def apply_fdr(corr: pd.DataFrame, cohort: str, role: str) -> None:
    mask = (corr["cohort"] == cohort) & (corr["role"] == role) & corr["p"].notna()
    if mask.sum() == 0:
        return
    corr.loc[mask, "q_bh"] = bh_fdr(corr.loc[mask, "p"].astype(float).tolist())


def scatter(ax, frame: pd.DataFrame, ycol: str, ylabel: str, colors: dict) -> None:
    for group, sub in frame.groupby("group"):
        ax.scatter(
            sub["CLDN4"],
            sub[ycol],
            s=18,
            alpha=0.85,
            c=colors.get(group, "#333333"),
            label=f"{group} (n={len(sub)})",
            linewidths=0,
        )
    a = pd.concat([
        pd.to_numeric(frame["CLDN4"], errors="coerce"),
        pd.to_numeric(frame[ycol], errors="coerce"),
    ], axis=1).dropna()
    if len(a) >= 8:
        slope, intercept, *_ = stats.linregress(a.iloc[:, 0], a.iloc[:, 1])
        xs = np.linspace(a.iloc[:, 0].min(), a.iloc[:, 0].max(), 50)
        ax.plot(xs, intercept + slope * xs, color="#222222", lw=1.0, zorder=0)
        rho, p = stats.spearmanr(a.iloc[:, 0], a.iloc[:, 1])
        ax.set_title(f"ρ={rho:+.3f}  p={p:.2g}  n={len(a)}", fontsize=10)
    ax.set_xlabel("CLDN4 log2(TPM+1)")
    ax.set_ylabel(ylabel)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data/ccle_cldn4_ddr_sting_hla")
    ap.add_argument("--outdir", default="methods/ccle_cldn4_ddr_sting_hla")
    args = ap.parse_args()
    data = Path(args.data)
    outdir = Path(args.outdir)
    figdir = outdir / "figures"
    tabdir = outdir / "tables"
    figdir.mkdir(parents=True, exist_ok=True)
    tabdir.mkdir(parents=True, exist_ok=True)

    manifest = json.loads((data / "download_manifest.json").read_text())
    model = pd.read_csv(data / "Model.csv")
    expr = pd.read_csv(data / "expression_panel.csv")
    dna_genes = [g.strip() for g in (data / "HALLMARK_DNA_REPAIR.txt").read_text().splitlines() if g.strip()]
    ifng_genes = [g.strip() for g in (data / "HALLMARK_INTERFERON_GAMMA_RESPONSE.txt").read_text().splitlines() if g.strip()]
    alias = resolve_aliases(list(expr.columns))

    df = expr.merge(model, on="ModelID", how="left")
    lung = df[(df["OncotreeLineage"] == "Lung") & (df["ModelType"] == "Cell Line")].copy()
    lung = lung.dropna(subset=["CLDN4"]).copy()
    lung["group"] = lung.apply(cohort_tag, axis=1)
    lung["is_NSCLC"] = lung["OncotreePrimaryDisease"] == "Non-Small Cell Lung Cancer"
    lung["is_SCLC_NET"] = lung["group"] == "SCLC_NET"
    lung["ccle_key"] = lung["CCLEName"].map(norm_id)

    lung, cov_all = add_scores(lung, ifng_genes, dna_genes, alias)
    cohorts = {
        "lung_cell_lines": lung,
        "NSCLC": add_scores(lung[lung["is_NSCLC"]].copy(), ifng_genes, dna_genes, alias)[0],
        "SCLC_NET": add_scores(lung[lung["is_SCLC_NET"]].copy(), ifng_genes, dna_genes, alias)[0],
        "LUAD": add_scores(lung[lung["group"] == "LUAD"].copy(), ifng_genes, dna_genes, alias)[0],
        "LUSC": add_scores(lung[lung["group"] == "LUSC"].copy(), ifng_genes, dna_genes, alias)[0],
        "other_NSCLC": add_scores(lung[lung["group"] == "other_NSCLC"].copy(), ifng_genes, dna_genes, alias)[0],
    }
    cov_nsclc = add_scores(lung[lung["is_NSCLC"]].copy(), ifng_genes, dna_genes, alias)[1]

    rows = []
    for cname, frame in cohorts.items():
        rows.extend(correlate_table(frame, RNA_TARGETS, cname))
    corr = pd.DataFrame(rows)
    for cname in ["lung_cell_lines", "NSCLC"]:
        apply_fdr(corr, cname, "rna_primary")
        apply_fdr(corr, cname, "hla_gene")

    # --- MCLP protein join (CLDN4 RNA vs RPPA). γH2AX column does not exist. ---
    rppa = pd.read_csv(data / "mclp_rppa_selected.csv")
    rppa["ccle_key"] = rppa["sample"].map(norm_id)
    n_missing_ccle_name = int(lung["ccle_key"].isna().sum())
    lung_keys = lung.dropna(subset=["ccle_key"]).copy()
    # One DepMap model per MCLP sample. True duplicate names are dropped, not averaged.
    key_counts = lung_keys["ccle_key"].value_counts()
    ambiguous = set(key_counts[key_counts > 1].index)
    lung_keys = lung_keys[~lung_keys["ccle_key"].isin(ambiguous)]
    prot = rppa.merge(
        lung_keys,
        on="ccle_key",
        how="inner",
        suffixes=("_rppa", ""),
    )
    # Protein values arrived from the RPPA side. After the merge the protein
    # columns keep their MCLP ids unless CLDN4-side names collided (they don't).
    prot["group"] = prot.apply(cohort_tag, axis=1)
    prot["is_NSCLC"] = prot["OncotreePrimaryDisease"] == "Non-Small Cell Lung Cancer"

    prot_cohorts = {
        "lung_cell_lines": prot,
        "NSCLC": prot[prot["is_NSCLC"]].copy(),
        "LUAD": prot[prot["group"] == "LUAD"].copy(),
        "LUSC": prot[prot["group"] == "LUSC"].copy(),
        "SCLC_NET": prot[prot["group"] == "SCLC_NET"].copy(),
    }
    prot_rows = []
    for cname, frame in prot_cohorts.items():
        prot_rows.extend(correlate_table(frame, PROTEIN_TARGETS, cname))
    pcorr = pd.DataFrame(prot_rows)
    for cname in ["lung_cell_lines", "NSCLC"]:
        apply_fdr(pcorr, cname, "protein_primary")

    # STING core vs Hallmark IFN-γ: are they the same score?
    sting_vs_ifn = spearman(lung["sig_STING_CORE"], lung["sig_IFNG"])

    # Subtype-adjusted partial Spearman. All-lung uses the five Oncotree groups;
    # NSCLC uses LUAD / LUSC / other NSCLC only.
    partial_targets = ["sig_DNA_REPAIR", "sig_DSB", "sig_STING_CORE", "sig_HLA_ABC", "HLA-A", "HLA-B", "HLA-C"]
    partial_rows = []
    for label, frame in (
        ("lung_cell_lines", lung),
        ("NSCLC", lung[lung["is_NSCLC"]].copy()),
    ):
        for target in partial_targets:
            rec = partial_spearman(frame["CLDN4"], frame[target], frame["group"])
            rec.update({"cohort": label, "predictor": "CLDN4", "target": target, "adjustment": "Oncotree group"})
            partial_rows.append(rec)
    partial = pd.DataFrame(partial_rows)

    corr.to_csv(tabdir / "correlations_rna.tsv", sep="\t", index=False)
    pcorr.to_csv(tabdir / "correlations_protein.tsv", sep="\t", index=False)
    partial.to_csv(tabdir / "partial_spearman_subtype.tsv", sep="\t", index=False)

    sample_cols = [
        "ModelID",
        "CellLineName",
        "CCLEName",
        "OncotreePrimaryDisease",
        "OncotreeSubtype",
        "group",
        "CLDN4",
        "TACSTD2",
        "sig_DNA_REPAIR",
        "sig_DSB",
        "sig_STING_CORE",
        "sig_STING_ISG",
        "sig_HLA_ABC",
        "sig_MHC1",
        "sig_HLA_II",
        "sig_IFNG",
        "HLA-A",
        "HLA-B",
        "HLA-C",
        "B2M",
        "STING1",
        "CGAS",
        "H2AX",
    ]
    keep = [c for c in sample_cols if c in lung.columns]
    lung[keep].sort_values(["group", "CellLineName"]).to_csv(tabdir / "lung_cell_lines.tsv", sep="\t", index=False)

    prot_keep = [
        "ModelID",
        "CellLineName",
        "sample",
        "group",
        "OncotreePrimaryDisease",
        "OncotreeSubtype",
        "CLDN4",
    ] + [t[0] for t in PROTEIN_TARGETS if t[0] in prot.columns]
    prot[prot_keep].sort_values(["group", "CellLineName"]).to_csv(
        tabdir / "lung_rppa_joined.tsv", sep="\t", index=False
    )

    counts = lung.groupby("group").size().rename("n_rna").reset_index()
    pcounts = prot.groupby("group").size().rename("n_rppa").reset_index()
    counts.merge(pcounts, on="group", how="outer").to_csv(tabdir / "cohort_counts.tsv", sep="\t", index=False)

    inventory = []
    for item in manifest["gamma_h2ax_inventory"]:
        inventory.append({
            "source": item.get("source"),
            "n_samples": item.get("n_samples"),
            "n_antibodies": item.get("n_antibodies"),
            "gamma_h2ax_present": item.get("gamma_h2ax_present"),
            "histone_antibodies": "; ".join(item.get("histone_antibodies") or []),
        })
    pd.DataFrame(inventory).to_csv(tabdir / "gammah2ax_inventory.tsv", sep="\t", index=False)

    def pack(table: pd.DataFrame, cohort: str, target: str) -> dict:
        hit = table[(table["cohort"] == cohort) & (table["target"] == target)]
        if hit.empty:
            return {}
        r = hit.iloc[0]
        def f(key):
            if key not in r or pd.isna(r[key]):
                return None
            if key == "n":
                return int(r[key])
            return float(r[key]) if key != "label" else r[key]
        return {
            "label": r["label"],
            "role": r["role"],
            "n": f("n"),
            "rho": None if pd.isna(r["rho"]) else float(r["rho"]),
            "p": None if pd.isna(r["p"]) else float(r["p"]),
            "q_bh": None if "q_bh" not in r or pd.isna(r["q_bh"]) else float(r["q_bh"]),
            "ci95_low": None if pd.isna(r["ci95_low"]) else float(r["ci95_low"]),
            "ci95_high": None if pd.isna(r["ci95_high"]) else float(r["ci95_high"]),
            "validation": r.get("validation") or "",
        }

    key = {
        "release": manifest.get("release"),
        "doi": manifest.get("doi"),
        "expression": "OmicsExpressionProteinCodingGenesTPMLogp1 = log2(TPM+1)",
        "lineage": "OncotreeLineage == Lung AND ModelType == Cell Line",
        "rppa": "MCLP CCLE RPPA release 20221116, norm level 4, joined on CCLE name",
        "gamma_h2ax": "absent from CCLE RPPA 20180123, CCLE RPPA 20181003, and MCLP 20221116",
        "n_lung_rna": int(len(lung)),
        "n_NSCLC_rna": int(lung["is_NSCLC"].sum()),
        "n_SCLC_NET_rna": int(lung["is_SCLC_NET"].sum()),
        "n_LUAD_rna": int((lung["group"] == "LUAD").sum()),
        "n_LUSC_rna": int((lung["group"] == "LUSC").sum()),
        "n_lung_rppa_joined": int(len(prot)),
        "n_NSCLC_rppa_joined": int(prot["is_NSCLC"].sum()),
        "n_missing_ccle_name": n_missing_ccle_name,
        "n_ambiguous_ccle_keys_dropped": int(len(ambiguous)),
        "ambiguous_ccle_keys": sorted(ambiguous),
        "partial_spearman_subtype": partial.to_dict(orient="records"),
        "n_mclp_lung_suffix": int(rppa["sample"].str.upper().str.endswith("_LUNG").sum()),
        "signature_coverage_all_lung": cov_all,
        "signature_coverage_NSCLC": cov_nsclc,
        "sting_core_vs_hallmark_ifng": sting_vs_ifn,
        "rna": {c: {t[0]: pack(corr, c, t[0]) for t in RNA_TARGETS} for c in cohorts},
        "protein": {c: {t[0]: pack(pcorr, c, t[0]) for t in PROTEIN_TARGETS} for c in prot_cohorts},
    }
    (tabdir / "key_stats.json").write_text(json.dumps(key, indent=2) + "\n")

    colors = {
        "LUAD": "#1f77b4",
        "LUSC": "#ff7f0e",
        "other_NSCLC": "#2ca02c",
        "SCLC_NET": "#d62728",
        "other_lung": "#7f7f7f",
    }

    fig, axes = plt.subplots(1, 3, figsize=(12.4, 4.15), constrained_layout=True)
    scatter(axes[0], lung, "sig_DNA_REPAIR", "Hallmark DNA repair (mean-z)", colors)
    scatter(axes[1], lung, "sig_STING_CORE", "STING core (mean-z)", colors)
    scatter(axes[2], lung, "sig_HLA_ABC", "HLA-A/B/C (mean-z)", colors)
    handles, labels = axes[2].get_legend_handles_labels()
    axes[2].legend_ = None
    fig.legend(handles, labels, frameon=False, fontsize=7, loc="lower center", ncol=5, bbox_to_anchor=(0.5, -0.08))
    fig.suptitle("DepMap 24Q4 lung cell lines: CLDN4 RNA", fontsize=12)
    fig.savefig(figdir / "fig_cldn4_rna_ddr_sting_hla.png", dpi=160, bbox_inches="tight")
    fig.savefig(figdir / "fig_cldn4_rna_ddr_sting_hla.pdf", bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(12.4, 4.15), constrained_layout=True)
    scatter(axes[0], prot, "ATMPS1981", "ATM pS1981 (MCLP, Valid)", colors)
    scatter(axes[1], prot, "RAD17PS645", "Rad17 pS645 (MCLP, Valid)", colors)
    scatter(axes[2], prot, "IRF3", "IRF3 protein (MCLP, Valid)", colors)
    handles, labels = axes[2].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, fontsize=7, loc="lower center", ncol=5, bbox_to_anchor=(0.5, -0.08))
    fig.suptitle("CLDN4 RNA vs MCLP RPPA (γH2AX antibody absent)", fontsize=12)
    fig.savefig(figdir / "fig_cldn4_rppa_ddr_sting.png", dpi=160, bbox_inches="tight")
    fig.savefig(figdir / "fig_cldn4_rppa_ddr_sting.pdf", bbox_inches="tight")
    plt.close(fig)

    # Forest: primary RNA + primary Valid protein on all-lung.
    forest_rows = []
    for target, label, role in [
        ("sig_DNA_REPAIR", "Hallmark DNA repair", "RNA"),
        ("sig_DSB", "DSB-response RNA", "RNA"),
        ("sig_STING_CORE", "STING core", "RNA"),
        ("sig_HLA_ABC", "HLA-A/B/C", "RNA"),
        ("HLA-A", "HLA-A", "RNA"),
        ("HLA-B", "HLA-B", "RNA"),
        ("HLA-C", "HLA-C", "RNA"),
    ]:
        rec = pack(corr, "lung_cell_lines", target)
        rec["name"] = label
        rec["panel"] = role
        forest_rows.append(rec)
    for target, label, _role, _val in PROTEIN_TARGETS:
        if _role != "protein_primary":
            continue
        rec = pack(pcorr, "lung_cell_lines", target)
        rec["name"] = label
        rec["panel"] = "RPPA"
        forest_rows.append(rec)
    forest = pd.DataFrame(forest_rows).iloc[::-1]
    fig, ax = plt.subplots(figsize=(7.4, 5.6), constrained_layout=True)
    ypos = np.arange(len(forest))
    ax.axvline(0, color="#888888", lw=0.8)
    for i, r in enumerate(forest.itertuples()):
        color = "#1f77b4" if r.panel == "RNA" else "#d62728"
        ax.plot([r.ci95_low, r.ci95_high], [i, i], color=color, lw=1.6)
        ax.plot(r.rho, i, "o", color=color, ms=5)
    ax.set_yticks(ypos)
    ax.set_yticklabels([f"{r.name}  n={r.n}" for r in forest.itertuples()], fontsize=8)
    ax.set_xlabel("Spearman ρ vs CLDN4 (95% bootstrap CI)")
    ax.set_title("All lung lines")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.savefig(figdir / "fig_forest_cldn4.png", dpi=160)
    fig.savefig(figdir / "fig_forest_cldn4.pdf")
    plt.close(fig)

    # Console summary of primary rows.
    show = corr[(corr["cohort"] == "lung_cell_lines") & (corr["role"].isin(["rna_primary", "hla_gene", "replication"]))]
    print(show[["target", "n", "rho", "p", "q_bh"]].to_string(index=False))
    print("--- protein ---")
    showp = pcorr[(pcorr["cohort"] == "lung_cell_lines") & (pcorr["role"] == "protein_primary")]
    print(showp[["target", "n", "rho", "p", "q_bh"]].to_string(index=False))
    print("n_lung", len(lung), "n_rppa", len(prot), "sting_vs_ifn", sting_vs_ifn)
    print(partial.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
