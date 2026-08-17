#!/usr/bin/env python3
"""CCLE/DepMap lung protein extra: CLDN4 vs CD274 / IFN proteins.

RNA Spearman ρ is already reported in methods/depmap_cldn4_ifn (PR #304).
This script does not re-estimate that RNA table. It inventories public
protein matrices and, where both proteins exist, scores CLDN4 protein
against CD274 and IFN / MHC-I proteins with pairwise-complete n.
"""
from __future__ import annotations

import argparse
import gzip
import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

MHC1 = ["HLA-A", "HLA-B", "HLA-C", "B2M"]
# Compact IFN/ISG cassette used when Hallmark coverage is sparse.
# Hallmark IFN-γ members that are actually present are preferred for the score.
IFN_CORE = [
    "STAT1",
    "STAT2",
    "IRF1",
    "IRF9",
    "MX1",
    "ISG15",
    "IFIT1",
    "IFIT3",
    "OAS1",
    "OAS2",
    "OAS3",
    "EIF2AK2",
    "IFI35",
    "BST2",
    "SAMHD1",
    "TRIM25",
    "ADAR",
    "TAP1",
    "TAP2",
    "PSMB8",
    "PSMB9",
    "PSMB10",
]
CD274_ALIASES = {"CD274", "PDCD1LG1", "PD-L1", "PDL1", "CD274_HUMAN", "PD1L1_HUMAN"}
CLDN4_ALIASES = {"CLDN4", "CLD4_HUMAN", "CLAUDIN-4", "CLAUDIN4"}
# Do not substring-match "PDL1" (that also hits SPDL1).
HALLMARK_IFNG = "HALLMARK_INTERFERON_GAMMA_RESPONSE"
TENPX = re.compile(r"_TenPx\d+$")


def zscore(s: pd.Series) -> pd.Series:
    mu = s.mean()
    sd = s.std(ddof=0)
    if sd == 0 or not np.isfinite(sd):
        return s * 0.0
    return (s - mu) / sd


def signature(df: pd.DataFrame, genes: list[str], min_members: int = 5) -> tuple[pd.Series, list[str]]:
    use = [g for g in genes if g in df.columns]
    if not use:
        return pd.Series(np.nan, index=df.index), []
    z = df[use].apply(zscore, axis=0)
    n_ok = z.notna().sum(axis=1)
    score = z.mean(axis=1)
    score = score.where(n_ok >= min_members)
    return score, use


def spearman(x: pd.Series, y: pd.Series) -> dict:
    a = pd.concat([x, y], axis=1).dropna()
    n = int(len(a))
    if n < 8:
        return {"n": n, "rho": None, "p": None, "ci95_low": None, "ci95_high": None}
    rho, p = stats.spearmanr(a.iloc[:, 0], a.iloc[:, 1])
    return {
        "n": n,
        "rho": float(rho),
        "p": float(p),
        "ci95_low": None,
        "ci95_high": None,
    }


def bootstrap_ci(x: pd.Series, y: pd.Series, n_boot: int = 5000, seed: int = 0) -> tuple[float | None, float | None]:
    a = pd.concat([x, y], axis=1).dropna()
    if len(a) < 8:
        return None, None
    rng = np.random.default_rng(seed)
    xv = a.iloc[:, 0].to_numpy()
    yv = a.iloc[:, 1].to_numpy()
    n = len(xv)
    boots = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        boots[i] = stats.spearmanr(xv[idx], yv[idx]).statistic
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


def parse_gmt_ifng(path: Path) -> list[str]:
    for line in path.read_text().splitlines():
        parts = line.rstrip("\n").split("\t")
        if parts and parts[0] == HALLMARK_IFNG and len(parts) >= 3:
            return [g for g in parts[2:] if g]
    raise SystemExit(f"{HALLMARK_IFNG} not found in {path}")


def map_gygi_columns(sample_cols: list[str], model: pd.DataFrame) -> dict[str, str]:
    out: dict[str, str] = {}
    ccle_map = {str(v): i for i, v in model["CCLEName"].dropna().items()}
    stripped_map = {str(v).upper(): i for i, v in model["StrippedCellLineName"].dropna().items()}
    for c in sample_cols:
        if c.startswith("ACH-"):
            out[c] = c.split("_")[0]
            continue
        base = TENPX.sub("", c)
        if base in ccle_map:
            out[c] = ccle_map[base]
            continue
        token = base.split("_")[0].upper()
        if token in stripped_map:
            out[c] = stripped_map[token]
    return out


def load_gygi_panel(path: Path, wanted: set[str]) -> tuple[pd.DataFrame, dict]:
    """Return samples × genes (symbols) plus inventory metadata."""
    import csv

    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", newline="") as fh:
        header = next(csv.reader(fh))
    if "Gene_Symbol" not in header:
        raise SystemExit(f"no Gene_Symbol in {path}; first cols={header[:8]}")
    meta_cols = {
        "Protein_Id",
        "Gene_Symbol",
        "Description",
        "Group_ID",
        "Uniprot",
        "Uniprot_Acc",
    }
    sample_cols = [
        c
        for c in header
        if c not in meta_cols and not c.endswith("_Peptides") and not c.startswith("TenPx")
    ]
    usecols = ["Gene_Symbol", "Protein_Id"] + sample_cols
    df = pd.read_csv(path, usecols=usecols)
    df["gene"] = df["Gene_Symbol"].astype(str).str.split(" ").str[0]
    all_genes = sorted(set(df["gene"].dropna()))
    cd274_hits = [g for g in all_genes if g.upper() in CD274_ALIASES or g.upper() in {"PD-L1", "PDL1"}]
    cldn4_hits = [g for g in all_genes if g.upper() in CLDN4_ALIASES or g.upper() == "CLDN4"]
    keep = set(wanted) | set(cd274_hits) | set(cldn4_hits)
    sub = df[df["gene"].isin(keep)].copy()
    num = sub[sample_cols].apply(pd.to_numeric, errors="coerce")
    num["gene"] = sub["gene"].values
    prot = num.groupby("gene").mean(numeric_only=True)
    inventory = {
        "n_proteins_in_file": int(df["gene"].nunique()),
        "n_quant_columns": len(sample_cols),
        "CLDN4_symbols": cldn4_hits,
        "CD274_symbols": cd274_hits,
        "wanted_present": sorted(g for g in wanted if g in prot.index),
        "wanted_absent": sorted(g for g in wanted if g not in prot.index),
        "protein_ids": (
            sub.groupby("gene")["Protein_Id"]
            .first()
            .to_dict()
        ),
    }
    return prot, sample_cols, inventory


def scan_rppa(ab_path: Path, mat_path: Path) -> dict:
    ab = pd.read_csv(ab_path)
    mat = pd.read_csv(mat_path)
    id_col = mat.columns[0]
    mat = mat.rename(columns={id_col: "CCLE_ID"})
    ab_str = ab.astype(str)
    def hits(pat: str) -> pd.DataFrame:
        mask = ab_str.apply(lambda c: c.str.contains(pat, case=False, na=False)).any(axis=1)
        return ab.loc[mask]
    cldn = hits(r"CLDN4|Claudin-4|Claudin4|Claudin 4")
    pdl1 = hits(r"CD274|PD-L1|PDL1|PD-L1|PDCD1LG1")
    ifn = hits(r"STAT1|STAT2|IRF1|IRF9|MX1|ISG15|IFIT|HLA-|B2M|PD-L1|PDL1|CD274|IFN")
    lung = mat[mat["CCLE_ID"].astype(str).str.contains("_LUNG", na=False)]
    return {
        "n_samples": int(len(mat)),
        "n_lung_suffix": int(len(lung)),
        "n_antibodies": int(len(ab)),
        "has_CLDN4": bool(len(cldn)),
        "has_CD274_or_PDL1": bool(len(pdl1)),
        "CLDN4_rows": cldn.astype(str).to_dict(orient="records"),
        "CD274_rows": pdl1.astype(str).to_dict(orient="records"),
        "IFN_related_antibodies": ifn.astype(str).to_dict(orient="records"),
        "note": "RPPA n is a headcount. No CLDN4 antibody means no CLDN4–CD274 protein ρ on this table.",
    }


def scan_procan(mat_path: Path, map_path: Path) -> dict:
    # First column is protein id; look at index/names only.
    header = pd.read_csv(mat_path, sep="\t", nrows=0)
    proteins = pd.read_csv(mat_path, sep="\t", usecols=[0])
    prot_col = proteins.columns[0]
    names = proteins[prot_col].astype(str)
    def has(pat: str) -> list[str]:
        return names[names.str.contains(pat, case=False, na=False)].tolist()
    mapping = pd.read_csv(map_path, sep="\t")
    lung_n = None
    nsclc_n = None
    for col in mapping.columns:
        s = mapping[col].astype(str)
        if s.str.contains("NSCLC", case=False, na=False).any():
            nsclc_n = int(s.str.contains("NSCLC", case=False, na=False).sum())
        if (s.str.lower() == "lung").any() or s.str.contains(r"\bLung\b", case=False, na=False).any():
            lung_n = int(s.str.contains(r"lung", case=False, na=False).sum())
    return {
        "n_proteins": int(len(names)),
        "n_samples": int(len(header.columns) - 1),
        "CLDN4_hits": has(r"CLDN4|CLD4_HUMAN"),
        "CD274_hits": has(r"CD274|PDCD1LG1|PD1L1_HUMAN|PD-L1"),
        "has_CLDN4": bool(has(r"CLDN4|CLD4_HUMAN")),
        "has_CD274": bool(has(r"CD274|PDCD1LG1|PD1L1_HUMAN")),
        "n_mapping_Lung_guess": lung_n,
        "n_mapping_NSCLC_guess": nsclc_n,
        "note": "Prior B2 audit: TACSTD2 present, CLDN4 absent. Re-checked here.",
    }


def fmt_rho(rec: dict) -> str:
    if rec.get("rho") is None:
        return f"n={rec.get('n', 0)}; not tested (n<8)"
    qbit = ""
    if rec.get("q_bh") is not None:
        qbit = f", q={rec['q_bh']:.2g}"
    ci = ""
    if rec.get("ci95_low") is not None:
        ci = f"; CI [{rec['ci95_low']:+.2f}, {rec['ci95_high']:+.2f}]"
    return f"{rec['rho']:+.3f} (n={rec['n']}, p={rec['p']:.2g}{qbit}){ci}"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data/ccle_cldn4_protein_ifn")
    ap.add_argument("--outdir", default="methods/ccle_cldn4_protein_ifn")
    args = ap.parse_args()
    data = Path(args.data)
    outdir = Path(args.outdir)
    figdir = outdir / "figures"
    tabdir = outdir / "tables"
    figdir.mkdir(parents=True, exist_ok=True)
    tabdir.mkdir(parents=True, exist_ok=True)

    ifng_genes = parse_gmt_ifng(data / "h.all.v2024.1.Hs.symbols.gmt")
    (data / "hallmark_ifng_genes.txt").write_text("\n".join(ifng_genes) + "\n")
    wanted = set(MHC1 + IFN_CORE + ["CLDN4", "TACSTD2", "CD274"] + ifng_genes)

    model = pd.read_csv(data / "Model.csv")
    model = model.set_index("ModelID", drop=False)

    gygi_path = data / "protein_quant_current_normalized.csv.gz"
    prot, sample_cols, gygi_inv = load_gygi_panel(gygi_path, wanted)
    col_to_model = map_gygi_columns(sample_cols, model)
    gygi_inv["n_columns_mapped_to_ModelID"] = len(col_to_model)
    gygi_inv["n_columns_unmapped"] = int(len(sample_cols) - len(col_to_model))

    # samples × genes
    wide = prot.T.copy()
    wide.index.name = "gygi_col"
    wide["ModelID"] = wide.index.map(col_to_model)
    wide = wide.dropna(subset=["ModelID"])
    # One ModelID can have replicate TenPx columns; mean them.
    num_genes = [c for c in wide.columns if c != "ModelID"]
    wide[num_genes] = wide[num_genes].apply(pd.to_numeric, errors="coerce")
    by_model = wide.groupby("ModelID")[num_genes].mean()

    ann = by_model.join(model[["CellLineName", "CCLEName", "OncotreeLineage", "OncotreePrimaryDisease", "OncotreeSubtype", "ModelType"]], how="left")
    lung = ann[(ann["OncotreeLineage"] == "Lung") & (ann["ModelType"] == "Cell Line")].copy()
    lung["group"] = lung.apply(cohort_tag, axis=1)
    lung["is_NSCLC"] = lung["OncotreePrimaryDisease"] == "Non-Small Cell Lung Cancer"
    lung["is_SCLC_NET"] = lung["group"] == "SCLC_NET"

    # Protein presence in lung
    present_lung = {
        g: int(lung[g].notna().sum()) if g in lung.columns else 0
        for g in ["CLDN4", "TACSTD2", "CD274"] + MHC1 + IFN_CORE
    }
    hallmark_in_gygi = [g for g in ifng_genes if g in lung.columns]
    hallmark_lung_n = {g: int(lung[g].notna().sum()) for g in hallmark_in_gygi}

    rppa_inv = scan_rppa(data / "CCLE_RPPA_Ab_info_20180123.csv", data / "CCLE_RPPA_20180123.csv")
    procan_inv = scan_procan(
        data / "ProCan-DepMapSanger_protein_matrix_8498_averaged.txt",
        data / "ProCan-DepMapSanger_mapping_file_averaged.txt",
    )

    def add_scores(frame: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
        out = frame.copy()
        # Hallmark IFN-γ proteins present in this frame (any non-NA).
        hm = [g for g in ifng_genes if g in out.columns and out[g].notna().any()]
        # Require the protein to have ≥8 values in this cohort before entering the score.
        hm = [g for g in hm if int(out[g].notna().sum()) >= 8]
        sig_ifn, used_ifn = signature(out, hm, min_members=5)
        sig_mhc, used_mhc = signature(out, MHC1, min_members=2)
        isg = [g for g in IFN_CORE if g in out.columns and int(out[g].notna().sum()) >= 8]
        sig_isg, used_isg = signature(out, isg, min_members=5)
        out["sig_IFNG_protein"] = sig_ifn
        out["sig_MHC1_protein"] = sig_mhc
        out["sig_ISG_protein"] = sig_isg
        cov = {
            "ifng_proteins_used": used_ifn,
            "ifng_n_used": len(used_ifn),
            "ifng_n_gmt": len(ifng_genes),
            "mhc1_proteins_used": used_mhc,
            "mhc1_n_used": len(used_mhc),
            "isg_proteins_used": used_isg,
            "isg_n_used": len(used_isg),
        }
        return out, cov

    lung_scored, cov_all = add_scores(lung)
    # Restrict primary protein tests to lines with CLDN4 protein.
    lung_cldn4 = lung_scored.dropna(subset=["CLDN4"]).copy() if "CLDN4" in lung_scored.columns else lung_scored.iloc[0:0].copy()

    def subset_scored(mask: pd.Series) -> pd.DataFrame:
        """Recompute z-scores only when the subset is large enough (n≥8)."""
        sub = lung_cldn4.loc[mask].copy()
        if len(sub) >= 8:
            return add_scores(sub)[0]
        return sub

    nsclc_scored, cov_nsclc = (
        add_scores(lung_cldn4[lung_cldn4["is_NSCLC"]].copy()) if len(lung_cldn4) else (lung_cldn4, {})
    )
    cohorts = {
        "lung_CLDN4_protein": lung_cldn4,
        "NSCLC_CLDN4_protein": nsclc_scored,
        "LUAD_CLDN4_protein": subset_scored(lung_cldn4["group"] == "LUAD"),
        "LUSC_CLDN4_protein": subset_scored(lung_cldn4["group"] == "LUSC"),
        "SCLC_NET_CLDN4_protein": subset_scored(lung_cldn4["is_SCLC_NET"]),
    }

    score_targets = ["sig_IFNG_protein", "sig_MHC1_protein", "sig_ISG_protein"]
    if "CD274" in lung.columns:
        score_targets.append("CD274")
    gene_targets = [g for g in MHC1 + IFN_CORE + ["TACSTD2"] if g in lung.columns and g != "CLDN4"]

    rows = []
    for cname, frame in cohorts.items():
        if "CLDN4" not in frame.columns or frame["CLDN4"].notna().sum() == 0:
            continue
        for target in score_targets + gene_targets:
            if target not in frame.columns:
                rec = {"n": 0, "rho": None, "p": None, "ci95_low": None, "ci95_high": None}
            else:
                rec = spearman(frame["CLDN4"], frame[target])
                if cname in {"lung_CLDN4_protein", "NSCLC_CLDN4_protein"} and rec["n"] >= 8:
                    lo, hi = bootstrap_ci(frame["CLDN4"], frame[target])
                    rec["ci95_low"] = lo
                    rec["ci95_high"] = hi
            rec.update(
                {
                    "dataset": "Gygi_CCLE_MS",
                    "cohort": cname,
                    "predictor": "CLDN4_protein",
                    "target": target,
                    "target_class": (
                        "score"
                        if target.startswith("sig_")
                        else ("self" if target == "TACSTD2" else "protein")
                    ),
                }
            )
            rows.append(rec)

    corr = pd.DataFrame(rows)
    # BH within the three primary axes on all-lung and NSCLC (Hallmark IFN score, MHC-I score, CD274).
    # Compact ISG score is a sensitivity and is excluded from this FDR.
    primary_axes = [t for t in ["sig_IFNG_protein", "sig_MHC1_protein", "CD274"] if t in set(corr["target"])]
    for cname in ["lung_CLDN4_protein", "NSCLC_CLDN4_protein"]:
        mask = (corr["cohort"] == cname) & (corr["target"].isin(primary_axes))
        if mask.sum() == 0:
            continue
        p = corr.loc[mask, "p"]
        ok = p.notna()
        if ok.sum() == 0:
            continue
        corr.loc[mask & ok, "q_bh"] = bh_fdr(p[ok].astype(float).tolist())

    # Per-protein FDR among gene targets (not scores, not TACSTD2) on all-lung.
    mask = (corr["cohort"] == "lung_CLDN4_protein") & (corr["target_class"] == "protein") & corr["p"].notna()
    if mask.sum():
        corr.loc[mask, "q_bh_proteins"] = bh_fdr(corr.loc[mask, "p"].astype(float).tolist())

    corr.to_csv(tabdir / "correlations.tsv", sep="\t", index=False)

    # Sample table
    keep_cols = [
        c
        for c in [
            "CellLineName",
            "CCLEName",
            "OncotreePrimaryDisease",
            "OncotreeSubtype",
            "group",
            "CLDN4",
            "TACSTD2",
            "CD274",
            "sig_IFNG_protein",
            "sig_MHC1_protein",
            "sig_ISG_protein",
        ]
        + MHC1
        if c in lung_scored.columns
    ]
    lung_scored[keep_cols].sort_values(["group", "CellLineName"]).to_csv(
        tabdir / "lung_protein_lines.tsv", sep="\t"
    )

    counts = (
        lung_cldn4.groupby("group", dropna=False).size().rename("n_CLDN4_protein").reset_index()
        if len(lung_cldn4)
        else pd.DataFrame(columns=["group", "n_CLDN4_protein"])
    )
    counts_all = lung.groupby("group", dropna=False).size().rename("n_any_protein").reset_index()
    counts = counts_all.merge(counts, on="group", how="left")
    counts.to_csv(tabdir / "cohort_counts.tsv", sep="\t", index=False)

    coverage_rows = []
    for g, n in sorted(present_lung.items(), key=lambda kv: (-kv[1], kv[0])):
        coverage_rows.append(
            {
                "dataset": "Gygi_CCLE_MS",
                "gene": g,
                "n_lung_cell_lines_quantified": n,
                "n_lung_with_CLDN4_and_gene": (
                    int(lung.dropna(subset=["CLDN4", g]).shape[0]) if g in lung.columns and "CLDN4" in lung.columns else 0
                ),
                "in_hallmark_IFNG": g in set(ifng_genes),
            }
        )
    pd.DataFrame(coverage_rows).to_csv(tabdir / "protein_coverage.tsv", sep="\t", index=False)

    def pick(cname: str, target: str) -> dict:
        hit = corr[(corr["cohort"] == cname) & (corr["target"] == target)]
        if hit.empty:
            return {"n": 0, "rho": None, "p": None, "q_bh": None, "ci95_low": None, "ci95_high": None}
        r = hit.iloc[0]
        return {
            "n": int(r["n"]) if pd.notna(r["n"]) else 0,
            "rho": None if pd.isna(r.get("rho")) else float(r["rho"]),
            "p": None if pd.isna(r.get("p")) else float(r["p"]),
            "q_bh": None if "q_bh" not in r or pd.isna(r["q_bh"]) else float(r["q_bh"]),
            "ci95_low": None if pd.isna(r.get("ci95_low")) else float(r["ci95_low"]),
            "ci95_high": None if pd.isna(r.get("ci95_high")) else float(r["ci95_high"]),
        }

    inventory = {
        "gygi": gygi_inv,
        "rppa": rppa_inv,
        "procan": procan_inv,
        "hallmark_IFNG_proteins_in_gygi_lung": hallmark_lung_n,
        "n_lung_models_with_any_gygi_protein": int(len(lung)),
        "n_lung_models_with_CLDN4_protein": int(len(lung_cldn4)),
        "n_NSCLC_with_CLDN4_protein": int(lung_cldn4["is_NSCLC"].sum()) if len(lung_cldn4) else 0,
        "CD274_protein_in_gygi": bool(gygi_inv["CD274_symbols"]),
        "CLDN4_protein_in_gygi": bool(gygi_inv["CLDN4_symbols"]),
        "CD274_protein_in_rppa": rppa_inv["has_CD274_or_PDL1"],
        "CLDN4_protein_in_rppa": rppa_inv["has_CLDN4"],
        "CD274_protein_in_procan": procan_inv["has_CD274"],
        "CLDN4_protein_in_procan": procan_inv["has_CLDN4"],
        "same_matrix_CLDN4_and_CD274": bool(gygi_inv["CD274_symbols"]) and bool(gygi_inv["CLDN4_symbols"]),
    }
    (tabdir / "inventory.json").write_text(json.dumps(inventory, indent=2) + "\n")

    key = {
        "task": "ccle_cldn4_protein_ifn",
        "additive_to": "methods/depmap_cldn4_ifn RNA ρ (PR #304)",
        "release_model": "DepMap Public 24Q4 Model.csv",
        "primary_protein_matrix": "Nusinow/Gygi CCLE TMT-MS protein_quant_current_normalized.csv.gz",
        "lineage": "OncotreeLineage == Lung AND ModelType == Cell Line",
        "n_lung_any_protein": int(len(lung)),
        "n_lung_CLDN4_protein": int(len(lung_cldn4)),
        "n_NSCLC_CLDN4_protein": int(lung_cldn4["is_NSCLC"].sum()) if len(lung_cldn4) else 0,
        "n_LUAD_CLDN4_protein": int((lung_cldn4["group"] == "LUAD").sum()) if len(lung_cldn4) else 0,
        "n_LUSC_CLDN4_protein": int((lung_cldn4["group"] == "LUSC").sum()) if len(lung_cldn4) else 0,
        "n_SCLC_NET_CLDN4_protein": int(lung_cldn4["is_SCLC_NET"].sum()) if len(lung_cldn4) else 0,
        "signature_coverage_all_lung": cov_all,
        "signature_coverage_NSCLC": cov_nsclc,
        "CLDN4_protein_vs": {
            cname: {t: pick(cname, t) for t in score_targets} for cname in cohorts if cname in set(corr["cohort"])
        },
        "inventory_flags": {
            "CD274_protein_in_gygi": inventory["CD274_protein_in_gygi"],
            "same_matrix_CLDN4_and_CD274": inventory["same_matrix_CLDN4_and_CD274"],
            "CLDN4_in_rppa": inventory["CLDN4_protein_in_rppa"],
            "CD274_in_rppa": inventory["CD274_protein_in_rppa"],
            "CLDN4_in_procan": inventory["CLDN4_protein_in_procan"],
            "CD274_in_procan": inventory["CD274_protein_in_procan"],
        },
        "notes": [
            "RNA ρ is given (PR #304). This folder is protein extra.",
            "n is pairwise-complete protein, not the RNA n=214.",
            "Hallmark IFN-γ protein score = mean of within-cohort z-scores; sample needs ≥5 members.",
            "MHC-I protein score = HLA-A/B/C + B2M; sample needs ≥2 members.",
            "No T-cell infiltrate. No IFN treatment. TMT abundance, not surface FACS.",
        ],
    }
    (tabdir / "key_stats.json").write_text(json.dumps(key, indent=2) + "\n")

    # --- figures ---
    color = {
        "LUAD": "#1f77b4",
        "LUSC": "#ff7f0e",
        "other_NSCLC": "#2ca02c",
        "SCLC_NET": "#d62728",
        "other_lung": "#7f7f7f",
    }
    panels = [
        ("sig_IFNG_protein", "Hallmark IFN-γ protein score (mean-z)"),
        ("sig_MHC1_protein", "MHC-I protein score (HLA-A/B/C + B2M)"),
    ]
    if "CD274" in lung_cldn4.columns and lung_cldn4["CD274"].notna().sum() >= 8:
        panels.append(("CD274", "CD274 protein (Gygi TMT)"))
    fig, axes = plt.subplots(1, len(panels), figsize=(4.1 * len(panels), 4.0), constrained_layout=True)
    if len(panels) == 1:
        axes = [axes]
    for ax, (col, ylab) in zip(axes, panels):
        for g, sub in lung_cldn4.groupby("group"):
            ax.scatter(
                sub["CLDN4"],
                sub[col],
                s=22,
                alpha=0.8,
                c=color.get(g, "#333333"),
                label=f"{g} n={sub[col].notna().sum()}",
                edgecolors="none",
            )
        rec = pick("lung_CLDN4_protein", col)
        ax.set_xlabel("CLDN4 protein (Gygi TMT)")
        ax.set_ylabel(ylab)
        rho_s = "NA" if rec["rho"] is None else f"{rec['rho']:.2f}"
        ax.set_title(f"all lung CLDN4+ n={rec['n']}\nρ={rho_s}", fontsize=9)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=5, frameon=False, fontsize=8, bbox_to_anchor=(0.5, -0.10))
    fig.suptitle("CCLE MS lung lines — CLDN4 protein vs IFN / MHC-I protein", fontsize=11, y=1.08)
    fig.savefig(figdir / "fig_cldn4_vs_ifn_mhc_protein.png", dpi=160, bbox_inches="tight")
    fig.savefig(figdir / "fig_cldn4_vs_ifn_mhc_protein.pdf", bbox_inches="tight")
    plt.close(fig)

    # Forest of per-protein associations (all-lung, n≥8)
    forest = corr[
        (corr["cohort"] == "lung_CLDN4_protein")
        & (corr["target_class"].isin(["protein", "score"]))
        & corr["rho"].notna()
    ].copy()
    forest = forest.sort_values("rho")
    if len(forest):
        fig, ax = plt.subplots(figsize=(7.4, max(3.0, 0.28 * len(forest) + 1.2)))
        y = np.arange(len(forest))
        ax.axvline(0, color="#444", lw=0.8)
        ax.scatter(forest["rho"], y, s=28, c="#1f4e79", zorder=3)
        labels = []
        for _, r in forest.iterrows():
            star = " *" if pd.notna(r.get("p")) and r["p"] < 0.05 else ""
            labels.append(f"{r['target']}  n={int(r['n'])}{star}")
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_xlabel("Spearman ρ vs CLDN4 protein")
        ax.set_title("Gygi CCLE MS lung — CLDN4 protein vs IFN/MHC proteins")
        ax.invert_yaxis()
        fig.tight_layout()
        fig.savefig(figdir / "fig_cldn4_protein_forest.png", dpi=160)
        fig.savefig(figdir / "fig_cldn4_protein_forest.pdf")
        plt.close(fig)

    print(
        json.dumps(
            {
                "n_lung_any_protein": key["n_lung_any_protein"],
                "n_lung_CLDN4_protein": key["n_lung_CLDN4_protein"],
                "n_NSCLC_CLDN4_protein": key["n_NSCLC_CLDN4_protein"],
                "CD274_in_gygi": inventory["CD274_protein_in_gygi"],
                "CD274_in_rppa": inventory["CD274_protein_in_rppa"],
                "CLDN4_in_rppa": inventory["CLDN4_protein_in_rppa"],
                "CLDN4_in_procan": inventory["CLDN4_protein_in_procan"],
                "primary": key["CLDN4_protein_vs"].get("lung_CLDN4_protein", {}),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
