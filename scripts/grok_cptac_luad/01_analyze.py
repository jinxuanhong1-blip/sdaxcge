#!/usr/bin/env python3
"""CPTAC LUAD TACSTD2 / CLDN4 vs protein, RNA, and immune scores.

Treatment-naive surgical cohort (Gillette et al. Cell 2020; pan-cancer freeze
v1.2). No ICI response labels exist — do not invent them.

Inputs: data/grok_cptac_luad/ (from 00_download.py)
Outputs: results/grok_cptac_luad/{tables,figures}/
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.duration.hazard_regression import PHReg
from statsmodels.stats.multitest import multipletests

TACSTD2 = "ENSG00000184292"
CLDN4 = "ENSG00000189143"

# RNA / protein signature members (GENCODE-style Ensembl, verified in freeze).
SIGS = {
    "CYT_cytolytic": ["ENSG00000145649", "ENSG00000180644"],  # GZMA, PRF1
    "CD8_Tcell": ["ENSG00000153563", "ENSG00000198851", "ENSG00000167286"],  # CD8A, CD3E, CD3D
    "IFNG_6gene": [
        "ENSG00000111537",  # IFNG
        "ENSG00000115415",  # STAT1
        "ENSG00000138755",  # CXCL9
        "ENSG00000169245",  # CXCL10
        "ENSG00000131203",  # IDO1
        "ENSG00000204287",  # HLA-DRA
    ],
    "GEP_Tcell_inflamed": [
        "ENSG00000271503",  # CCL5
        "ENSG00000139193",  # CD27
        "ENSG00000120217",  # CD274
        "ENSG00000103855",  # CD276
        "ENSG00000153563",  # CD8A
        "ENSG00000174600",  # CMKLR1
        "ENSG00000138755",  # CXCL9
        "ENSG00000172215",  # CXCR6
        "ENSG00000196735",  # HLA-DQA1
        "ENSG00000196126",  # HLA-DRB1
        "ENSG00000204592",  # HLA-E
        "ENSG00000131203",  # IDO1
        "ENSG00000089692",  # LAG3
        "ENSG00000105374",  # NKG7
        "ENSG00000197646",  # PDCD1LG2
        "ENSG00000205220",  # PSMB10
        "ENSG00000115415",  # STAT1
        "ENSG00000181847",  # TIGIT
    ],
    "TLS_chemokine": [
        "ENSG00000271503",  # CCL5
        "ENSG00000172724",  # CCL19
        "ENSG00000137077",  # CCL21
        "ENSG00000138755",  # CXCL9
        "ENSG00000169245",  # CXCL10
        "ENSG00000169248",  # CXCL11
        "ENSG00000156234",  # CXCL13
    ],
    "TGFB_axis": ["ENSG00000105329", "ENSG00000163513"],  # TGFB1, TGFBR2
}

# Pre-specified immune phenotype columns (Bessede-style infiltration / IFN / exclusion).
PRIMARY_IMMUNE = [
    "ESTIMATE_ImmuneScore",
    "ESTIMATE_StromalScore",
    "CIBERSORT_T_cell_CD8+",
    "CIBERSORT_T_cell_regulatory_(Tregs)",
    "CIBERSORT_Macrophage_M2",
    "CIBERSORT_Neutrophil",
    "xCell_T_cell_CD8+",
    "xCell_immune_score",
    "xCell_T_cell_regulatory_(Tregs)",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_TGF_BETA_SIGNALING",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "PROGENy_JAK-STAT",
    "PROGENy_TGFb",
    "PROGENy_NFkB",
    "TMB",
    "WES_purity",
]


def find_row(index: pd.Index, prefix: str) -> str | None:
    hits = [i for i in index if str(i) == prefix or str(i).startswith(prefix + ".")]
    if not hits:
        return None
    if len(hits) > 1:
        raise ValueError(f"multiple rows for {prefix}: {hits}")
    return hits[0]


def load_matrix(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", index_col=0)


def extract_gene(mat: pd.DataFrame, prefix: str) -> pd.Series:
    row = find_row(mat.index, prefix)
    if row is None:
        return pd.Series(np.nan, index=mat.columns, name=prefix)
    s = pd.to_numeric(mat.loc[row], errors="coerce")
    s.name = prefix
    return s


def mean_z_score(mat: pd.DataFrame, prefixes: list[str], min_genes: int = 2) -> pd.Series:
    rows = []
    used = []
    for p in prefixes:
        rid = find_row(mat.index, p)
        if rid is None:
            continue
        s = pd.to_numeric(mat.loc[rid], errors="coerce")
        if s.notna().sum() < 10:
            continue
        mu, sd = s.mean(), s.std(ddof=0)
        if sd == 0 or np.isnan(sd):
            continue
        rows.append((s - mu) / sd)
        used.append(rid)
    if not rows:
        return pd.Series(np.nan, index=mat.columns), used
    z = pd.concat(rows, axis=1)
    n = z.notna().sum(axis=1)
    score = z.mean(axis=1, skipna=True)
    score[n < min_genes] = np.nan
    return score, used


def spearman(x: pd.Series, y: pd.Series) -> dict:
    d = pd.concat([x, y], axis=1).dropna()
    n = len(d)
    if n < 4:
        return {"n": n, "rho": np.nan, "p": np.nan}
    if d.iloc[:, 0].nunique() < 2 or d.iloc[:, 1].nunique() < 2:
        return {"n": n, "rho": np.nan, "p": np.nan}
    rho, p = stats.spearmanr(d.iloc[:, 0], d.iloc[:, 1])
    return {"n": int(n), "rho": float(rho), "p": float(p)}


def mannwhitney(a: pd.Series, b: pd.Series) -> dict:
    a = pd.to_numeric(a, errors="coerce").dropna()
    b = pd.to_numeric(b, errors="coerce").dropna()
    if len(a) < 3 or len(b) < 3:
        return {
            "n_a": int(len(a)),
            "n_b": int(len(b)),
            "median_a": float(a.median()) if len(a) else np.nan,
            "median_b": float(b.median()) if len(b) else np.nan,
            "U": np.nan,
            "p": np.nan,
            "rank_biserial": np.nan,
        }
    U, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    # Rank-biserial r = 1 - 2U/(n1 n2); sign: positive if a > b on ranks.
    r = 1.0 - (2.0 * U) / (len(a) * len(b))
    # scipy U is for a vs b; convert so positive r means a stochastically greater
    r = -r
    return {
        "n_a": int(len(a)),
        "n_b": int(len(b)),
        "median_a": float(a.median()),
        "median_b": float(b.median()),
        "U": float(U),
        "p": float(p),
        "rank_biserial": float(r),
    }


def wilcoxon_paired(a: pd.Series, b: pd.Series) -> dict:
    d = pd.concat([a, b], axis=1).dropna()
    d.columns = ["a", "b"]
    if len(d) < 5:
        return {"n": int(len(d)), "median_a": np.nan, "median_b": np.nan, "W": np.nan, "p": np.nan}
    try:
        W, p = stats.wilcoxon(d["a"], d["b"], alternative="two-sided", zero_method="wilcox")
    except ValueError:
        return {
            "n": int(len(d)),
            "median_a": float(d["a"].median()),
            "median_b": float(d["b"].median()),
            "W": np.nan,
            "p": np.nan,
        }
    return {
        "n": int(len(d)),
        "median_a": float(d["a"].median()),
        "median_b": float(d["b"].median()),
        "median_delta_a_minus_b": float((d["a"] - d["b"]).median()),
        "W": float(W),
        "p": float(p),
    }


def partial_spearman(x: pd.Series, y: pd.Series, z: pd.Series) -> dict:
    """Spearman of residuals after ranking + linear residualization on z."""
    d = pd.concat([x, y, z], axis=1).dropna()
    d.columns = ["x", "y", "z"]
    n = len(d)
    if n < 6 or d["x"].nunique() < 2 or d["y"].nunique() < 2 or d["z"].nunique() < 2:
        return {"n": n, "rho": np.nan, "p": np.nan}
    rx = d["x"].rank()
    ry = d["y"].rank()
    rz = d["z"].rank()
    # residualize ranks on z ranks
    bx = np.polyfit(rz, rx, 1)
    by = np.polyfit(rz, ry, 1)
    ex = rx - (bx[0] * rz + bx[1])
    ey = ry - (by[0] * rz + by[1])
    rho, p = stats.spearmanr(ex, ey)
    return {"n": int(n), "rho": float(rho), "p": float(p)}


def logrank(time: pd.Series, event: pd.Series, group: pd.Series) -> dict:
    d = pd.concat([time, event, group], axis=1).dropna()
    d.columns = ["t", "e", "g"]
    d["t"] = pd.to_numeric(d["t"], errors="coerce")
    d["e"] = pd.to_numeric(d["e"], errors="coerce")
    d = d.dropna()
    labs = sorted(d["g"].unique())
    if len(labs) != 2:
        return {"n": int(len(d)), "n_high": np.nan, "n_low": np.nan, "events": int(d["e"].sum()) if len(d) else 0, "stat": np.nan, "p": np.nan}
    g0, g1 = labs
    t0, e0 = d.loc[d.g == g0, "t"], d.loc[d.g == g0, "e"]
    t1, e1 = d.loc[d.g == g1, "t"], d.loc[d.g == g1, "e"]
    try:
        stat, p = stats.logrank(t0, t1, e0, e1)
        stat = float(np.asarray(stat).ravel()[0])
        p = float(np.asarray(p).ravel()[0])
    except Exception:
        # fallback: two-sample logrank via contingency at event times
        times = np.sort(d.loc[d.e == 1, "t"].unique())
        o1 = e1_exp = v = 0.0
        for t in times:
            at0 = ((t0 >= t).sum())
            at1 = ((t1 >= t).sum())
            at = at0 + at1
            dth = ((d.t == t) & (d.e == 1)).sum()
            d1 = ((t1 == t) & (e1 == 1)).sum()
            if at <= 1:
                continue
            exp1 = dth * at1 / at
            var = dth * (at - dth) / (at - 1) * (at1 / at) * (at0 / at)
            o1 += d1
            e1_exp += exp1
            v += var
        if v <= 0:
            stat, p = np.nan, np.nan
        else:
            stat = (o1 - e1_exp) ** 2 / v
            p = float(stats.chi2.sf(stat, 1))
    return {
        "n": int(len(d)),
        "n_high": int((d.g == "high").sum()) if "high" in set(d.g) else int((d.g == g1).sum()),
        "n_low": int((d.g == "low").sum()) if "low" in set(d.g) else int((d.g == g0).sum()),
        "events": int(d["e"].sum()),
        "stat": float(stat) if stat == stat else np.nan,
        "p": float(p) if p == p else np.nan,
    }


def km_steps(time: pd.Series, event: pd.Series) -> tuple[np.ndarray, np.ndarray]:
    d = pd.concat([time, event], axis=1).dropna()
    d.columns = ["t", "e"]
    d = d.sort_values("t")
    t = d["t"].to_numpy(float)
    e = d["e"].to_numpy(float)
    uniq = np.unique(t)
    surv = 1.0
    xs = [0.0]
    ys = [1.0]
    n = len(d)
    for ti in uniq:
        at = (t >= ti).sum()
        di = ((t == ti) & (e == 1)).sum()
        if at > 0 and di > 0:
            surv *= 1.0 - di / at
        xs.extend([ti, ti])
        ys.extend([ys[-1], surv])
    return np.array(xs), np.array(ys)


def fdr(pvals: list[float]) -> list[float]:
    arr = np.asarray(pvals, dtype=float)
    out = np.full_like(arr, np.nan, dtype=float)
    ok = np.isfinite(arr)
    if ok.sum() == 0:
        return out.tolist()
    out[ok] = multipletests(arr[ok], method="fdr_bh")[1]
    return [float(x) if x == x else np.nan for x in out]


def style_ax(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=8)


def savefig(fig, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data/grok_cptac_luad")
    ap.add_argument("--outdir", default="results/grok_cptac_luad")
    args = ap.parse_args()
    data = Path(args.data)
    out = Path(args.outdir)
    tables = out / "tables"
    figs = out / "figures"
    tables.mkdir(parents=True, exist_ok=True)
    figs.mkdir(parents=True, exist_ok=True)

    rna_t = load_matrix(data / "LUAD_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt")
    rna_n = load_matrix(data / "LUAD_RNAseq_gene_RSEM_coding_UQ_1500_log2_Normal.txt")
    prot_t = load_matrix(data / "LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt")
    prot_n = load_matrix(data / "LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt")
    pheno = pd.read_csv(data / "LUAD_phenotype.txt", sep="\t", index_col=0)
    meta = pd.read_csv(data / "LUAD_meta.txt", sep="\t")
    meta = meta[meta["case_id"] != "data_type"].set_index("case_id")
    surv = pd.read_csv(data / "LUAD_survival.txt", sep="\t").set_index("case_id")

    samples = sorted(set(rna_t.columns) & set(prot_t.columns) & set(pheno.index))
    nat_rna = sorted(set(rna_t.columns) & set(rna_n.columns))
    nat_prot = sorted(set(prot_t.columns) & set(prot_n.columns))

    # Per-sample core table
    core = pd.DataFrame(index=samples)
    core["TACSTD2_RNA"] = extract_gene(rna_t, TACSTD2).reindex(samples)
    core["CLDN4_RNA"] = extract_gene(rna_t, CLDN4).reindex(samples)
    core["TACSTD2_protein"] = extract_gene(prot_t, TACSTD2).reindex(samples)
    core["CLDN4_protein"] = extract_gene(prot_t, CLDN4).reindex(samples)
    core["TACSTD2_RNA_NAT"] = extract_gene(rna_n, TACSTD2).reindex(samples)
    core["CLDN4_RNA_NAT"] = extract_gene(rna_n, CLDN4).reindex(samples)
    core["TACSTD2_protein_NAT"] = extract_gene(prot_n, TACSTD2).reindex(samples)
    core["CLDN4_protein_NAT"] = extract_gene(prot_n, CLDN4).reindex(samples)

    coverage = []
    for layer, mat, gene, col in [
        ("RNA_tumor", rna_t, TACSTD2, "TACSTD2_RNA"),
        ("RNA_tumor", rna_t, CLDN4, "CLDN4_RNA"),
        ("protein_tumor", prot_t, TACSTD2, "TACSTD2_protein"),
        ("protein_tumor", prot_t, CLDN4, "CLDN4_protein"),
        ("RNA_NAT", rna_n, TACSTD2, "TACSTD2_RNA_NAT"),
        ("RNA_NAT", rna_n, CLDN4, "CLDN4_RNA_NAT"),
        ("protein_NAT", prot_n, TACSTD2, "TACSTD2_protein_NAT"),
        ("protein_NAT", prot_n, CLDN4, "CLDN4_protein_NAT"),
    ]:
        rid = find_row(mat.index, gene)
        s = core[col]
        coverage.append(
            {
                "layer": layer,
                "symbol": "TACSTD2" if gene == TACSTD2 else "CLDN4",
                "ensembl_query": gene,
                "row_id": rid,
                "n_samples_in_matrix": int(mat.shape[1]),
                "n_non_na": int(s.notna().sum()),
                "n_na": int(s.isna().sum()),
                "min": float(s.min()) if s.notna().any() else np.nan,
                "median": float(s.median()) if s.notna().any() else np.nan,
                "max": float(s.max()) if s.notna().any() else np.nan,
            }
        )
    cov_df = pd.DataFrame(coverage)
    cov_df.to_csv(tables / "gene_coverage.tsv", sep="\t", index=False)

    # Signature scores
    sig_cov = []
    for name, genes in SIGS.items():
        rs, rused = mean_z_score(rna_t, genes)
        ps, pused = mean_z_score(prot_t, genes)
        core[f"{name}_RNA"] = rs.reindex(samples)
        core[f"{name}_protein"] = ps.reindex(samples)
        sig_cov.append(
            {
                "signature": name,
                "n_genes_requested": len(genes),
                "RNA_genes_used": ";".join(rused),
                "RNA_n_used": len(rused),
                "protein_genes_used": ";".join(pused),
                "protein_n_used": len(pused),
            }
        )
    pd.DataFrame(sig_cov).to_csv(tables / "signature_coverage.tsv", sep="\t", index=False)

    # Join phenotype / clinical / survival
    keep_pheno = [c for c in pheno.columns if c.startswith(("CIBERSORT_", "ESTIMATE_", "xCell_", "PROGENy_", "HALLMARK_", "TMB", "WES_purity", "WGS_purity"))]
    core = core.join(pheno.reindex(samples)[keep_pheno])
    clin_cols = [
        "Age",
        "Sex",
        "Tumor_Size_cm",
        "Histologic_Grade",
        "Stage",
        "Tobacco_smoking_history",
        "KEAP1_mutation",
        "STK11_mutation",
        "KRAS_mutation",
        "TP53_mutation",
        "EGFR_mutation",
    ]
    core = core.join(meta.reindex(samples)[clin_cols])
    core = core.join(surv.reindex(samples))
    core.index.name = "case_id"
    core.to_csv(tables / "sample_table.tsv", sep="\t")

    # ---- Tumor vs NAT ----
    tvn = []
    for gene, tcol, ncol, layer in [
        ("TACSTD2", "TACSTD2_RNA", "TACSTD2_RNA_NAT", "RNA"),
        ("CLDN4", "CLDN4_RNA", "CLDN4_RNA_NAT", "RNA"),
        ("TACSTD2", "TACSTD2_protein", "TACSTD2_protein_NAT", "protein"),
        ("CLDN4", "CLDN4_protein", "CLDN4_protein_NAT", "protein"),
    ]:
        r = wilcoxon_paired(core[tcol], core[ncol])
        r.update({"gene": gene, "layer": layer, "comparison": "tumor_vs_NAT_paired"})
        tvn.append(r)
    tvn_df = pd.DataFrame(tvn)
    tvn_df.to_csv(tables / "tumor_vs_nat.tsv", sep="\t", index=False)

    # ---- RNA-protein and coexpression ----
    rp = []
    for gene, rcol, pcol in [
        ("TACSTD2", "TACSTD2_RNA", "TACSTD2_protein"),
        ("CLDN4", "CLDN4_RNA", "CLDN4_protein"),
    ]:
        r = spearman(core[rcol], core[pcol])
        r.update({"gene": gene, "test": "RNA_vs_protein"})
        rp.append(r)
    for layer, a, b in [
        ("RNA", "TACSTD2_RNA", "CLDN4_RNA"),
        ("protein", "TACSTD2_protein", "CLDN4_protein"),
    ]:
        r = spearman(core[a], core[b])
        r.update({"gene": "TACSTD2_vs_CLDN4", "test": f"coexpression_{layer}"})
        rp.append(r)
    rp_df = pd.DataFrame(rp)
    rp_df.to_csv(tables / "rna_protein_coexpression.tsv", sep="\t", index=False)

    # ---- Immune Spearman (all keep_pheno + signatures) ----
    predictors = [
        "TACSTD2_RNA",
        "CLDN4_RNA",
        "TACSTD2_protein",
        "CLDN4_protein",
    ]
    immune_targets = [c for c in keep_pheno if not c.startswith(("WES_", "WGS_"))]
    immune_targets += [c for c in core.columns if c.endswith(("_RNA", "_protein")) and any(c.startswith(s) for s in SIGS)]
    # unique preserve order
    seen = set()
    immune_targets = [c for c in immune_targets if not (c in seen or seen.add(c))]

    rows = []
    for pred in predictors:
        for tgt in immune_targets:
            if tgt == pred:
                continue
            r = spearman(core[pred], core[tgt])
            r.update({"predictor": pred, "target": tgt, "adjust": "none"})
            rows.append(r)
            if "WES_purity" in core.columns:
                pr = partial_spearman(core[pred], core[tgt], core["WES_purity"])
                pr.update({"predictor": pred, "target": tgt, "adjust": "WES_purity"})
                rows.append(pr)
    imm_df = pd.DataFrame(rows)
    # FDR within predictor × adjust
    imm_df["q"] = np.nan
    for (pred, adj), idx in imm_df.groupby(["predictor", "adjust"]).groups.items():
        imm_df.loc[idx, "q"] = fdr(imm_df.loc[idx, "p"].tolist())
    imm_df["primary"] = imm_df["target"].isin(PRIMARY_IMMUNE)
    imm_df = imm_df.sort_values(["predictor", "adjust", "p"])
    imm_df.to_csv(tables / "immune_spearman.tsv", sep="\t", index=False)

    primary_df = imm_df[imm_df["primary"]].copy()
    primary_df.to_csv(tables / "immune_spearman_primary.tsv", sep="\t", index=False)

    # ---- Median-split MW vs ImmuneScore / CD8 / IFNG ----
    mw_rows = []
    endpoints = [
        "ESTIMATE_ImmuneScore",
        "CIBERSORT_T_cell_CD8+",
        "xCell_T_cell_CD8+",
        "HALLMARK_INTERFERON_GAMMA_RESPONSE",
        "CYT_cytolytic_RNA",
        "GEP_Tcell_inflamed_RNA",
    ]
    for pred in predictors:
        med = core[pred].median(skipna=True)
        grp = pd.Series(np.where(core[pred] >= med, "high", "low"), index=core.index)
        grp[core[pred].isna()] = np.nan
        for ep in endpoints:
            if ep not in core.columns:
                continue
            r = mannwhitney(core.loc[grp == "high", ep], core.loc[grp == "low", ep])
            r.update(
                {
                    "predictor": pred,
                    "endpoint": ep,
                    "split": "median",
                    "median_cut": float(med) if med == med else np.nan,
                    "n_high": int((grp == "high").sum()),
                    "n_low": int((grp == "low").sum()),
                }
            )
            # rename a/b
            r["median_high"] = r.pop("median_a")
            r["median_low"] = r.pop("median_b")
            r["n_high_used"] = r.pop("n_a")
            r["n_low_used"] = r.pop("n_b")
            mw_rows.append(r)
    mw_df = pd.DataFrame(mw_rows)
    mw_df["q"] = fdr(mw_df["p"].tolist())
    mw_df.to_csv(tables / "median_split_immune.tsv", sep="\t", index=False)

    # ---- Clinical associations ----
    clin_rows = []
    for pred in predictors:
        # smoking: never vs ever
        sm = core["Tobacco_smoking_history"].astype(str)
        never = core.loc[sm.str.contains("non-smoker", case=False, na=False), pred]
        ever = core.loc[sm.isin(["current smoker", "past smoker"]), pred]
        r = mannwhitney(ever, never)
        r.update({"predictor": pred, "variable": "smoking_ever_vs_never", "group_a": "ever", "group_b": "never"})
        clin_rows.append(r)
        for mut in ["STK11_mutation", "KEAP1_mutation", "EGFR_mutation", "KRAS_mutation", "TP53_mutation"]:
            a = core.loc[core[mut].astype(str) == "1", pred]
            b = core.loc[core[mut].astype(str) == "0", pred]
            r = mannwhitney(a, b)
            r.update({"predictor": pred, "variable": mut, "group_a": "mut", "group_b": "wt"})
            clin_rows.append(r)
        # sex
        r = mannwhitney(core.loc[core["Sex"] == "Male", pred], core.loc[core["Sex"] == "Female", pred])
        r.update({"predictor": pred, "variable": "Sex_male_vs_female", "group_a": "Male", "group_b": "Female"})
        clin_rows.append(r)
        # stage I vs II-IV
        st = core["Stage"].astype(str)
        r = mannwhitney(core.loc[st.str.contains("Stage I$|Stage I ", na=False) | (st == "Stage I"), pred], core.loc[st.isin(["Stage II", "Stage III", "Stage IV"]), pred])
        # fix stage I
        stage_i = st.eq("Stage I")
        later = st.isin(["Stage II", "Stage III", "Stage IV"])
        r = mannwhitney(core.loc[later, pred], core.loc[stage_i, pred])
        r.update({"predictor": pred, "variable": "stage_IIIV_vs_I", "group_a": "II-IV", "group_b": "I"})
        clin_rows.append(r)
    clin_df = pd.DataFrame(clin_rows)
    clin_df["q"] = fdr(clin_df["p"].tolist())
    clin_df.to_csv(tables / "clinical_associations.tsv", sep="\t", index=False)

    # ---- Survival ----
    surv_rows = []
    for pred in predictors:
        x = core[pred]
        med = x.median(skipna=True)
        grp = pd.Series(np.where(x >= med, "high", "low"), index=core.index)
        grp[x.isna()] = np.nan
        for endpoint, tcol, ecol in [("OS", "OS_days", "OS_event"), ("PFS", "PFS_days", "PFS_event")]:
            lr = logrank(core[tcol], core[ecol], grp)
            rec = {"predictor": pred, "endpoint": endpoint, "split": "median", "median_cut": float(med) if med == med else np.nan}
            rec.update(lr)
            # Cox continuous (per 1 SD)
            d = pd.concat([core[tcol], core[ecol], x], axis=1).dropna()
            d.columns = ["t", "e", "x"]
            rec["cox_n"] = int(len(d))
            rec["cox_events"] = int(d["e"].sum()) if len(d) else 0
            if len(d) >= 20 and d["e"].sum() >= 5 and d["x"].std() > 0:
                z = (d["x"] - d["x"].mean()) / d["x"].std(ddof=0)
                try:
                    model = PHReg(d["t"].to_numpy(float), z.to_numpy(float)[:, None], status=d["e"].to_numpy(float))
                    res = model.fit(disp=0)
                    rec["cox_hr_per_sd"] = float(np.exp(res.params[0]))
                    rec["cox_p"] = float(res.pvalues[0])
                    rec["cox_loghr"] = float(res.params[0])
                    rec["cox_loghr_se"] = float(res.bse[0])
                except Exception as exc:
                    rec["cox_hr_per_sd"] = np.nan
                    rec["cox_p"] = np.nan
                    rec["cox_error"] = str(exc)
            else:
                rec["cox_hr_per_sd"] = np.nan
                rec["cox_p"] = np.nan
            surv_rows.append(rec)
    surv_df = pd.DataFrame(surv_rows)
    surv_df.to_csv(tables / "survival.tsv", sep="\t", index=False)

    # ---- JSON key stats for writeup ----
    def row_get(df, **kw):
        q = df
        for k, v in kw.items():
            q = q[q[k] == v]
        return q.iloc[0].to_dict() if len(q) else {}

    key = {
        "cohort": "CPTAC LUAD pan-cancer freeze v1.2",
        "treatment": "treatment-naive surgical tumors; no ICI labels",
        "n_tumor": len(samples),
        "n_NAT_RNA": len(nat_rna),
        "n_NAT_protein": len(nat_prot),
        "coverage": coverage,
        "tumor_vs_nat": tvn,
        "rna_protein_coexpression": rp,
        "primary_immune_unadjusted": primary_df[primary_df.adjust == "none"].to_dict(orient="records"),
        "median_split": mw_df.to_dict(orient="records"),
        "clinical": clin_df.to_dict(orient="records"),
        "survival": surv_df.to_dict(orient="records"),
    }
    (tables / "key_stats.json").write_text(json.dumps(key, indent=2, default=str) + "\n")

    # ================= figures =================
    # 1. Tumor vs NAT
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.4))
    pairs = [
        (axes[0, 0], "TACSTD2_RNA", "TACSTD2_RNA_NAT", "TACSTD2 RNA (log2 RSEM-UQ)"),
        (axes[0, 1], "CLDN4_RNA", "CLDN4_RNA_NAT", "CLDN4 RNA (log2 RSEM-UQ)"),
        (axes[1, 0], "TACSTD2_protein", "TACSTD2_protein_NAT", "TACSTD2 protein (log2 TMT)"),
        (axes[1, 1], "CLDN4_protein", "CLDN4_protein_NAT", "CLDN4 protein (log2 TMT)"),
    ]
    for ax, tcol, ncol, title in pairs:
        d = core[[tcol, ncol]].dropna()
        data_box = [d[ncol].to_numpy(), d[tcol].to_numpy()]
        bp = ax.boxplot(data_box, tick_labels=["NAT", "Tumor"], widths=0.55, patch_artist=True, showfliers=False)
        for patch, color in zip(bp["boxes"], ["#9bb7d4", "#d4896a"]):
            patch.set_facecolor(color)
            patch.set_alpha(0.85)
        for i in range(len(d)):
            ax.plot([1, 2], [d.iloc[i, 1], d.iloc[i, 0]], color="0.7", lw=0.3, alpha=0.5)
        rec = [x for x in tvn if x["layer"] in title and x["gene"] in title]
        # annotate from tvn_df
        sub = tvn_df[(tvn_df.gene == title.split()[0]) & (tvn_df.layer == ("RNA" if "RNA" in title else "protein"))]
        if len(sub):
            p = sub.iloc[0]["p"]
            n = int(sub.iloc[0]["n"])
            ax.set_title(f"{title}\npaired Wilcoxon n={n} p={p:.2e}", fontsize=8)
        else:
            ax.set_title(title, fontsize=8)
        style_ax(ax)
    savefig(fig, figs / "fig1_tumor_vs_nat.png")

    # 2. RNA vs protein
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.4))
    for ax, gene, rcol, pcol in [
        (axes[0], "TACSTD2", "TACSTD2_RNA", "TACSTD2_protein"),
        (axes[1], "CLDN4", "CLDN4_RNA", "CLDN4_protein"),
    ]:
        d = core[[rcol, pcol]].dropna()
        ax.scatter(d[rcol], d[pcol], s=16, c="#3b6d9a", alpha=0.75, edgecolors="none")
        r = spearman(core[rcol], core[pcol])
        ax.set_xlabel(f"{gene} RNA")
        ax.set_ylabel(f"{gene} protein")
        ax.set_title(f"{gene} RNA–protein  ρ={r['rho']:.3f} p={r['p']:.2e} n={r['n']}", fontsize=8)
        style_ax(ax)
    savefig(fig, figs / "fig2_rna_vs_protein.png")

    # 3. Coexpression
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.4))
    for ax, layer, a, b in [
        (axes[0], "RNA", "TACSTD2_RNA", "CLDN4_RNA"),
        (axes[1], "protein", "TACSTD2_protein", "CLDN4_protein"),
    ]:
        d = core[[a, b]].dropna()
        ax.scatter(d[a], d[b], s=16, c="#6a4c93", alpha=0.75, edgecolors="none")
        r = spearman(core[a], core[b])
        ax.set_xlabel(f"TACSTD2 {layer}")
        ax.set_ylabel(f"CLDN4 {layer}")
        ax.set_title(f"TACSTD2 vs CLDN4 {layer}  ρ={r['rho']:.3f} p={r['p']:.2e} n={r['n']}", fontsize=8)
        style_ax(ax)
    savefig(fig, figs / "fig3_tacstd2_cldn4_coexpression.png")

    # 4. Primary immune forest (protein + RNA TACSTD2)
    def forest(ax, pred, title):
        sub = imm_df[(imm_df.predictor == pred) & (imm_df.adjust == "none") & (imm_df.primary)]
        sub = sub.sort_values("rho")
        y = np.arange(len(sub))
        colors = ["#b33" if r < 0 else "#2a7" for r in sub["rho"]]
        ax.barh(y, sub["rho"], color=colors, height=0.7, alpha=0.85)
        ax.set_yticks(y)
        ax.set_yticklabels(sub["target"], fontsize=6.5)
        ax.axvline(0, color="0.3", lw=0.6)
        ax.set_xlabel("Spearman ρ")
        ax.set_title(title, fontsize=8)
        style_ax(ax)

    fig, axes = plt.subplots(2, 2, figsize=(10.5, 8.2))
    forest(axes[0, 0], "TACSTD2_protein", "TACSTD2 protein vs primary immune")
    forest(axes[0, 1], "CLDN4_protein", "CLDN4 protein vs primary immune")
    forest(axes[1, 0], "TACSTD2_RNA", "TACSTD2 RNA vs primary immune")
    forest(axes[1, 1], "CLDN4_RNA", "CLDN4 RNA vs primary immune")
    savefig(fig, figs / "fig4_primary_immune_spearman.png")

    # 5. Heatmap: predictors × selected immune
    heat_targets = [
        "ESTIMATE_ImmuneScore",
        "xCell_immune_score",
        "CIBERSORT_T_cell_CD8+",
        "xCell_T_cell_CD8+",
        "CIBERSORT_T_cell_regulatory_(Tregs)",
        "CIBERSORT_Macrophage_M2",
        "CIBERSORT_Neutrophil",
        "HALLMARK_INTERFERON_GAMMA_RESPONSE",
        "HALLMARK_TGF_BETA_SIGNALING",
        "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
        "PROGENy_JAK-STAT",
        "CYT_cytolytic_RNA",
        "GEP_Tcell_inflamed_RNA",
        "TLS_chemokine_RNA",
        "CYT_cytolytic_protein",
        "WES_purity",
        "TMB",
    ]
    heat_targets = [t for t in heat_targets if t in core.columns]
    M = np.zeros((len(predictors), len(heat_targets)))
    P = np.zeros_like(M)
    for i, pred in enumerate(predictors):
        for j, tgt in enumerate(heat_targets):
            r = spearman(core[pred], core[tgt])
            M[i, j] = r["rho"]
            P[i, j] = r["p"]
    fig, ax = plt.subplots(figsize=(10.2, 3.6))
    im = ax.imshow(M, cmap="RdBu_r", vmin=-0.6, vmax=0.6, aspect="auto")
    ax.set_xticks(range(len(heat_targets)))
    ax.set_xticklabels(heat_targets, rotation=55, ha="right", fontsize=7)
    ax.set_yticks(range(len(predictors)))
    ax.set_yticklabels(predictors, fontsize=8)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            star = "***" if P[i, j] < 0.001 else "**" if P[i, j] < 0.01 else "*" if P[i, j] < 0.05 else ""
            ax.text(j, i, f"{M[i, j]:.2f}{star}", ha="center", va="center", fontsize=5.5, color="black")
    fig.colorbar(im, ax=ax, shrink=0.8, label="Spearman ρ")
    ax.set_title("CPTAC LUAD n=110 tumors (CLDN4 protein pairwise-complete ~79)")
    savefig(fig, figs / "fig5_immune_heatmap.png")

    # 6. KM
    fig, axes = plt.subplots(2, 2, figsize=(7.6, 6.4))
    km_plan = [
        (axes[0, 0], "TACSTD2_protein", "OS"),
        (axes[0, 1], "CLDN4_protein", "OS"),
        (axes[1, 0], "TACSTD2_protein", "PFS"),
        (axes[1, 1], "CLDN4_protein", "PFS"),
    ]
    for ax, pred, endpoint in km_plan:
        tcol, ecol = ("OS_days", "OS_event") if endpoint == "OS" else ("PFS_days", "PFS_event")
        x = core[pred]
        med = x.median(skipna=True)
        grp = pd.Series(np.where(x >= med, "high", "low"), index=core.index)
        grp[x.isna()] = np.nan
        for lab, color in [("low", "#3b6d9a"), ("high", "#c04b3a")]:
            m = grp == lab
            xs, ys = km_steps(core.loc[m, tcol], core.loc[m, ecol])
            n = int(m.sum())
            ev = int(pd.to_numeric(core.loc[m, ecol], errors="coerce").fillna(0).sum())
            ax.step(xs, ys, where="post", color=color, lw=1.5, label=f"{lab} n={n} ev={ev}")
        rec = surv_df[(surv_df.predictor == pred) & (surv_df.endpoint == endpoint)]
        p = rec.iloc[0]["p"] if len(rec) else np.nan
        ax.set_ylim(0, 1.02)
        ax.set_xlabel("Days")
        ax.set_ylabel(f"{endpoint} survival")
        ax.set_title(f"{pred} median split  logrank p={p:.3g}", fontsize=8)
        ax.legend(fontsize=6, frameon=False)
        style_ax(ax)
    savefig(fig, figs / "fig6_survival_km.png")

    # 7. Clinical boxes: STK11 and smoking for TACSTD2 protein
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4))
    # STK11
    ax = axes[0]
    a = core.loc[core["STK11_mutation"].astype(str) == "1", "TACSTD2_protein"].dropna()
    b = core.loc[core["STK11_mutation"].astype(str) == "0", "TACSTD2_protein"].dropna()
    bp = ax.boxplot([b, a], tick_labels=["STK11 WT", "STK11 mut"], widths=0.55, patch_artist=True, showfliers=True)
    for patch, color in zip(bp["boxes"], ["#9bb7d4", "#d4896a"]):
        patch.set_facecolor(color)
    r = mannwhitney(a, b)
    ax.set_ylabel("TACSTD2 protein")
    ax.set_title(f"STK11  MW p={r['p']:.3g}  n_mut={r['n_a']} n_wt={r['n_b']}", fontsize=8)
    style_ax(ax)
    ax = axes[1]
    never = core.loc[core["Tobacco_smoking_history"].astype(str).str.contains("non-smoker", case=False, na=False), "TACSTD2_protein"].dropna()
    ever = core.loc[core["Tobacco_smoking_history"].isin(["current smoker", "past smoker"]), "TACSTD2_protein"].dropna()
    bp = ax.boxplot([never, ever], tick_labels=["Never", "Ever"], widths=0.55, patch_artist=True, showfliers=True)
    for patch, color in zip(bp["boxes"], ["#9bb7d4", "#d4896a"]):
        patch.set_facecolor(color)
    r = mannwhitney(ever, never)
    ax.set_ylabel("TACSTD2 protein")
    ax.set_title(f"Smoking ever vs never  MW p={r['p']:.3g}", fontsize=8)
    style_ax(ax)
    savefig(fig, figs / "fig7_clinical_tacstd2_protein.png")

    # 8. Missingness bar for CLDN4 protein
    fig, ax = plt.subplots(figsize=(5.2, 3.2))
    labels = ["TACSTD2\nRNA", "CLDN4\nRNA", "TACSTD2\nprotein", "CLDN4\nprotein"]
    non = [
        int(core["TACSTD2_RNA"].notna().sum()),
        int(core["CLDN4_RNA"].notna().sum()),
        int(core["TACSTD2_protein"].notna().sum()),
        int(core["CLDN4_protein"].notna().sum()),
    ]
    na = [110 - x for x in non]
    ax.bar(labels, non, color="#3b6d9a", label="detected")
    ax.bar(labels, na, bottom=non, color="#d9d9d9", label="NA")
    ax.set_ylabel("Tumor samples (n=110)")
    ax.set_title("Feature missingness in CPTAC LUAD freeze")
    ax.legend(frameon=False, fontsize=8)
    style_ax(ax)
    savefig(fig, figs / "fig8_missingness.png")

    print(f"Wrote tables → {tables}")
    print(f"Wrote figures → {figs}")
    print("CLDN4 protein NA:", int(core["CLDN4_protein"].isna().sum()), "/110")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
