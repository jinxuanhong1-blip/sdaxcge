#!/usr/bin/env python3
"""NicheNet-style ligand activity: CLDN4-high malignant -> T/NK IFN / cytotoxicity.

Unit of inference is the patient (n=12 post-treatment). Cells are never
treated as independent replicates of MPR/NMPR.

Sender definition is CLDN4-only (malignant CLDN4 log1p(CP10k) >= malignant
median). TACSTD2 is scored as a companion column and is not used to call
the sender set.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import average_precision_score, roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
CACHE = Path("/tmp/scrna_nichenet_cldn4")
OUT = ROOT / "results"
OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(ROOT / "scripts"))
from gene_sets import (  # noqa: E402
    CD8_TYPE,
    CYTOTOXICITY,
    EXHAUSTION,
    IFN,
    MALIGNANT,
    MPR_LABELS,
    NK_TYPE,
    NMPR_LABELS,
    TNK_TYPES,
)

ANN = DATA / "drmref_cell_annotation.tsv.gz"
META = DATA / "geo_scRNAseq_sample_metadata.tsv"
PANEL = CACHE / "extracted_panel.parquet"
LR = DATA / "lr_network.tsv"
LT = CACHE / "prior_ligand_target.parquet"

DETECT_FRAC = 0.10  # NicheNet 10x-style expressed-gene rule
TOP_N = 15
EMPIRICAL_P = 0.15  # honest: n=4 vs 8; do not pretend FDR<0.05
PRIMARY_SETS = ("a_priori_ifn", "a_priori_cytotoxicity")


def mpr_group(label: str) -> str | None:
    s = str(label).strip()
    if s in MPR_LABELS or s.upper() == "PCR":
        return "MPR"
    if s in NMPR_LABELS:
        return "NMPR"
    return None


def log1p_cp10k(umi: np.ndarray, ncount: np.ndarray) -> np.ndarray:
    ncount = np.maximum(ncount.astype(float), 1.0)
    return np.log1p(1e4 * umi.astype(float) / ncount)


def mean_score(df: pd.DataFrame, genes: list[str], prefix: str) -> pd.Series:
    cols = [f"{prefix}{g}" for g in genes if f"{prefix}{g}" in df.columns]
    if not cols:
        return pd.Series(np.nan, index=df.index)
    return df[cols].mean(axis=1)


def wilcoxon_exact_groups(values: pd.Series, groups: pd.Series, a: str, b: str) -> dict:
    """Two-sided Mann-Whitney; exact enumeration when n is small."""
    xa = values[groups == a].dropna().to_numpy()
    xb = values[groups == b].dropna().to_numpy()
    na, nb = len(xa), len(xb)
    if na == 0 or nb == 0:
        return {"n_a": na, "n_b": nb, "mean_a": np.nan, "mean_b": np.nan, "U": np.nan, "p": np.nan}
    res = stats.mannwhitneyu(xb, xa, alternative="two-sided", method="auto")
    return {
        "n_a": na,
        "n_b": nb,
        "mean_a": float(np.mean(xa)),
        "mean_b": float(np.mean(xb)),
        "median_a": float(np.median(xa)),
        "median_b": float(np.median(xb)),
        "U": float(res.statistic),
        "p": float(res.pvalue),
    }


def spearman_safe(x: pd.Series, y: pd.Series) -> dict:
    m = x.notna() & y.notna()
    n = int(m.sum())
    if n < 4:
        return {"n": n, "rho": np.nan, "p": np.nan}
    rho, p = stats.spearmanr(x[m], y[m])
    return {"n": n, "rho": float(rho), "p": float(p)}


def ligand_activity(lt: pd.DataFrame, geneset: set[str], background: list[str], ligands: list[str]) -> pd.DataFrame:
    genes = [g for g in background if g in lt.index]
    if len(genes) < 20 or not ligands:
        return pd.DataFrame()
    y = np.array([1 if g in geneset else 0 for g in genes], dtype=int)
    if y.sum() < 2 or y.sum() > len(y) - 2:
        return pd.DataFrame()
    use = [L for L in ligands if L in lt.columns]
    X = lt.loc[genes, use].to_numpy(dtype=float)
    rows = []
    for j, L in enumerate(use):
        s = X[:, j]
        if not np.isfinite(s).all() or np.allclose(s, s[0]):
            continue
        pear = stats.pearsonr(s, y)
        try:
            auroc = float(roc_auc_score(y, s))
        except ValueError:
            auroc = np.nan
        try:
            aupr = float(average_precision_score(y, s))
        except ValueError:
            aupr = np.nan
        rows.append(
            {
                "ligand": L,
                "pearson": float(pear.statistic),
                "pearson_p": float(pear.pvalue),
                "auroc": auroc,
                "aupr": aupr,
                "n_background": len(genes),
                "n_geneset": int(y.sum()),
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.sort_values("pearson", ascending=False).reset_index(drop=True)


def expressed_genes(sub: pd.DataFrame, genes: list[str], frac: float = DETECT_FRAC) -> list[str]:
    keep = []
    n = len(sub)
    if n == 0:
        return keep
    for g in genes:
        if g not in sub.columns:
            continue
        if float((sub[g] > 0).mean()) >= frac:
            keep.append(g)
    return keep


def add_log_cols(cells: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    ncount = cells["nCount_RNA"].to_numpy()
    logs = {}
    for g in genes:
        if g in cells.columns:
            logs[f"log_{g}"] = log1p_cp10k(cells[g].to_numpy(), ncount)
    if logs:
        cells = pd.concat([cells, pd.DataFrame(logs, index=cells.index)], axis=1)
    return cells


def main() -> int:
    ann = pd.read_csv(ANN, sep="\t")
    meta = pd.read_csv(META, sep="\t")
    panel = pd.read_parquet(PANEL)
    lr = pd.read_csv(LR, sep="\t")
    lt = pd.read_parquet(LT)
    if lt.shape[0] < lt.shape[1]:
        print(f"warning: lt shape {lt.shape} (expected genes x ligands)", flush=True)

    meta = meta.rename(columns={"Sample": "orig.ident", "Patient": "patient_meta"})
    meta["mpr_group"] = meta["Pathologic_Response"].map(mpr_group)
    meta["timing"] = np.where(
        meta["Resource"].astype(str).str.contains("Post", case=False), "post", "pre"
    )

    cells = ann.merge(panel, left_on="cell_barcode", right_index=True, how="inner")
    cells = cells.merge(
        meta[["orig.ident", "mpr_group", "timing", "Pathology", "Pathologic_Response"]],
        on="orig.ident",
        how="left",
    )
    cells["mpr_group"] = cells["mpr_group"].fillna(cells["response"].map(mpr_group))
    post = cells[cells["timing"] == "post"].copy()
    print(f"joined cells {len(cells)}; post {len(post)}", flush=True)

    genes_in_panel = [c for c in panel.columns]
    post = add_log_cols(post, genes_in_panel)

    mal = post[post["celltype"] == MALIGNANT].copy()
    tnk = post[post["celltype"].isin(TNK_TYPES)].copy()
    print(f"malignant {len(mal)} patients={mal.patient.nunique()}; T/NK {len(tnk)}", flush=True)

    # CLDN4-only high sender: CLDN4 >= malignant-cell median (pre-registered)
    if "log_CLDN4" not in mal.columns:
        raise RuntimeError("CLDN4 missing from extracted panel")
    med_c = float(mal["log_CLDN4"].median())
    med_t = float(mal["log_TACSTD2"].median()) if "log_TACSTD2" in mal.columns else np.nan
    mal["cldn4_high"] = mal["log_CLDN4"] >= med_c
    # companion only — not used to define senders
    if "log_TACSTD2" in mal.columns:
        mal["dual_high_companion"] = (mal["log_TACSTD2"] >= med_t) & (mal["log_CLDN4"] >= med_c)
    n_high = int(mal["cldn4_high"].sum())
    n_dual = int(mal["dual_high_companion"].sum()) if "dual_high_companion" in mal.columns else np.nan
    print(
        f"CLDN4-high malignant {n_high}/{len(mal)} med_CLDN4={med_c:.3f} "
        f"(companion dual-high {n_dual}; not used as sender)",
        flush=True,
    )

    rows = []
    for pid, g in post.groupby("patient"):
        mg = mal[mal["patient"] == pid]
        tg = tnk[tnk["patient"] == pid]
        rec = {
            "patient": pid,
            "mpr_group": g["mpr_group"].iloc[0],
            "pathology": g["Pathology"].iloc[0] if "Pathology" in g else "",
            "n_annotated": int(len(g)),
            "n_malignant": int(len(mg)),
            "n_cldn4_high": int(mg["cldn4_high"].sum()) if len(mg) else 0,
            "n_tnk": int(len(tg)),
            "n_cd8": int((g["celltype"] == CD8_TYPE).sum()),
            "n_nk": int((g["celltype"] == NK_TYPE).sum()),
            "tnk_frac": float((g["celltype"].isin(TNK_TYPES)).mean()),
        }
        if len(mg):
            rec["mal_CLDN4"] = float(mg["log_CLDN4"].mean())
            rec["pct_cldn4_high"] = float(mg["cldn4_high"].mean())
            if "log_TACSTD2" in mg.columns:
                rec["mal_TACSTD2_companion"] = float(mg["log_TACSTD2"].mean())
        if len(tg):
            rec["tnk_cyto"] = float(mean_score(tg, CYTOTOXICITY, "log_").mean())
            rec["tnk_ifn"] = float(mean_score(tg, IFN, "log_").mean())
            rec["tnk_exh"] = float(mean_score(tg, EXHAUSTION, "log_").mean())
            rec["tnk_ifn_minus_cyto"] = rec["tnk_ifn"] - rec["tnk_cyto"]
        rows.append(rec)
    sample = pd.DataFrame(rows).sort_values("patient")
    sample.to_csv(OUT / "sample_metrics.tsv", sep="\t", index=False)

    tests = []
    post_s = sample[sample["mpr_group"].isin(["MPR", "NMPR"])].copy()
    for metric in [
        "mal_CLDN4",
        "pct_cldn4_high",
        "tnk_cyto",
        "tnk_ifn",
        "tnk_ifn_minus_cyto",
        "tnk_exh",
        "tnk_frac",
    ]:
        if metric not in post_s.columns:
            continue
        w = wilcoxon_exact_groups(post_s[metric], post_s["mpr_group"], "MPR", "NMPR")
        tests.append({"test": f"NMPR_vs_MPR_{metric}", "contrast": "NMPR vs MPR", **w})
    for x, y in [
        ("mal_CLDN4", "tnk_cyto"),
        ("mal_CLDN4", "tnk_ifn"),
        ("mal_CLDN4", "tnk_ifn_minus_cyto"),
        ("mal_CLDN4", "tnk_exh"),
        ("mal_CLDN4", "tnk_frac"),
        ("pct_cldn4_high", "tnk_cyto"),
        ("pct_cldn4_high", "tnk_ifn"),
    ]:
        if x not in post_s.columns or y not in post_s.columns:
            continue
        sp = spearman_safe(post_s[x], post_s[y])
        tests.append({"test": f"spearman_{x}__{y}", "contrast": f"{x} vs {y}", **sp})

    ligands_all = sorted(set(lr["from"].astype(str)) & set(lt.columns.astype(str)))
    receptors_all = sorted(set(lr["to"].astype(str)))
    high = mal[mal["cldn4_high"]]
    low = mal[~mal["cldn4_high"]]
    expr_lig_high = expressed_genes(high, ligands_all)
    expr_lig_low = expressed_genes(low, ligands_all)
    expr_rec_tnk = expressed_genes(tnk, receptors_all)
    potential = sorted(
        set(lr.loc[lr["from"].isin(expr_lig_high) & lr["to"].isin(expr_rec_tnk), "from"])
    )
    print(
        f"expressed ligands CLDN4-high={len(expr_lig_high)} low={len(expr_lig_low)} "
        f"T/NK receptors={len(expr_rec_tnk)} potential={len(potential)}",
        flush=True,
    )

    high_mpr = high[high["mpr_group"] == "MPR"]
    high_nmpr = high[high["mpr_group"] == "NMPR"]
    tnk_mpr = tnk[tnk["mpr_group"] == "MPR"]
    tnk_nmpr = tnk[tnk["mpr_group"] == "NMPR"]
    pot_mpr = sorted(
        set(
            lr.loc[
                lr["from"].isin(expressed_genes(high_mpr, ligands_all))
                & lr["to"].isin(expressed_genes(tnk_mpr, receptors_all)),
                "from",
            ]
        )
    )
    pot_nmpr = sorted(
        set(
            lr.loc[
                lr["from"].isin(expressed_genes(high_nmpr, ligands_all))
                & lr["to"].isin(expressed_genes(tnk_nmpr, receptors_all)),
                "from",
            ]
        )
    )

    tnk_expr = expressed_genes(tnk, list(lt.index.intersection(panel.columns)))
    background = sorted(set(tnk_expr) | set(CYTOTOXICITY) | set(IFN) | set(EXHAUSTION))
    background = [g for g in background if g in lt.index]
    print(f"background genes (panel-restricted) {len(background)}", flush=True)

    gene_de = []
    tnk_genes = [g for g in (CYTOTOXICITY + IFN + EXHAUSTION + list(set(tnk_expr))) if f"log_{g}" in tnk.columns]
    tnk_genes = sorted(set(tnk_genes))
    pat_means = []
    for pid, g in tnk.groupby("patient"):
        rec = {"patient": pid, "mpr_group": g["mpr_group"].iloc[0]}
        for gene in tnk_genes:
            rec[gene] = float(g[f"log_{gene}"].mean())
        pat_means.append(rec)
    pat_tnk = pd.DataFrame(pat_means)
    for gene in tnk_genes:
        w = wilcoxon_exact_groups(pat_tnk[gene], pat_tnk["mpr_group"], "MPR", "NMPR")
        gene_de.append(
            {
                "gene": gene,
                "in_cyto": gene in CYTOTOXICITY,
                "in_ifn": gene in IFN,
                "in_exh": gene in EXHAUSTION,
                "mean_MPR": w["mean_a"],
                "mean_NMPR": w["mean_b"],
                "delta_NMPR_minus_MPR": (w["mean_b"] - w["mean_a"]) if np.isfinite(w["mean_a"]) else np.nan,
                "p": w["p"],
                "n_MPR": w["n_a"],
                "n_NMPR": w["n_b"],
            }
        )
    de = pd.DataFrame(gene_de).sort_values("p")
    de.to_csv(OUT / "tnk_patient_gene_means_NMPR_vs_MPR.tsv", sep="\t", index=False)
    emp_up = set(de.loc[(de["p"] < EMPIRICAL_P) & (de["delta_NMPR_minus_MPR"] > 0), "gene"])
    emp_down = set(de.loc[(de["p"] < EMPIRICAL_P) & (de["delta_NMPR_minus_MPR"] < 0), "gene"])

    sets = {
        "a_priori_ifn": set(IFN) & set(lt.index),
        "a_priori_cytotoxicity": set(CYTOTOXICITY) & set(lt.index),
        "a_priori_exhaustion_sensitivity": set(EXHAUSTION) & set(lt.index),
        "empirical_NMPR_TNK_up": emp_up & set(lt.index),
        "empirical_NMPR_TNK_down": emp_down & set(lt.index),
    }

    activities = []
    for set_name, gs in sets.items():
        for setting, ligs in [
            ("cldn4_high_to_all_TNK", potential),
            ("cldn4_high_NMPR_to_NMPR_TNK", pot_nmpr),
            ("cldn4_high_MPR_to_MPR_TNK", pot_mpr),
        ]:
            act = ligand_activity(lt, gs, background, ligs)
            if act.empty:
                continue
            act.insert(0, "setting", setting)
            act.insert(1, "geneset", set_name)
            activities.append(act)
    act_all = pd.concat(activities, ignore_index=True) if activities else pd.DataFrame()
    if not act_all.empty:
        act_all.to_csv(OUT / "ligand_activity_all.tsv", sep="\t", index=False)

    tops = []
    if not act_all.empty:
        for (setting, gs), sub in act_all.groupby(["setting", "geneset"]):
            chunk = sub.head(TOP_N).copy()
            chunk["rank"] = np.arange(1, len(chunk) + 1)
            tops.append(chunk)
    top = pd.concat(tops, ignore_index=True) if tops else pd.DataFrame()
    if not top.empty:
        top.to_csv(OUT / "top_ligands.tsv", sep="\t", index=False)
        primary = top[
            (top["setting"] == "cldn4_high_to_all_TNK") & (top["geneset"].isin(PRIMARY_SETS))
        ].copy()
        primary.to_csv(OUT / "top_ligands_primary.tsv", sep="\t", index=False)

    lig_pat = []
    focus_ligands = sorted(set(potential) | set(pot_mpr) | set(pot_nmpr))
    if not top.empty:
        focus_ligands = sorted(set(focus_ligands) | set(top["ligand"]))
    for pid, g in high.groupby("patient"):
        rec = {"patient": pid, "mpr_group": g["mpr_group"].iloc[0], "n_cldn4_high": int(len(g))}
        for L in focus_ligands:
            if f"log_{L}" in g.columns:
                rec[L] = float(g[f"log_{L}"].mean())
        lig_pat.append(rec)
    lig_pat_df = pd.DataFrame(lig_pat)
    lig_pat_df.to_csv(OUT / "cldn4_high_ligand_patient_means.tsv", sep="\t", index=False)

    lig_tests = []
    if not lig_pat_df.empty:
        for L in focus_ligands:
            if L not in lig_pat_df.columns:
                continue
            w = wilcoxon_exact_groups(lig_pat_df[L], lig_pat_df["mpr_group"], "MPR", "NMPR")
            lig_tests.append({"ligand": L, **w})
    lig_test_df = pd.DataFrame(lig_tests)
    if not lig_test_df.empty:
        lig_test_df = lig_test_df.sort_values("p")
        lig_test_df.to_csv(OUT / "cldn4_high_ligand_NMPR_vs_MPR.tsv", sep="\t", index=False)

    lig_corr = []
    merged = lig_pat_df.merge(post_s[["patient", "tnk_cyto", "tnk_ifn", "tnk_ifn_minus_cyto"]], on="patient")
    for L in focus_ligands:
        if L not in merged.columns:
            continue
        for y in ["tnk_cyto", "tnk_ifn", "tnk_ifn_minus_cyto"]:
            sp = spearman_safe(merged[L], merged[y])
            lig_corr.append({"ligand": L, "vs": y, **sp})
    lig_corr_df = pd.DataFrame(lig_corr)
    if not lig_corr_df.empty:
        lig_corr_df.to_csv(OUT / "ligand_vs_tnk_program_spearman.tsv", sep="\t", index=False)

    n_tbl = pd.DataFrame(
        [
            {"item": "GEO_cells_in_UMI_matrix", "n": 92330, "note": "Hu et al. GSE207422"},
            {"item": "DRMref_annotated_cells", "n": int(len(ann)), "note": "Liu et al. NAR 2024; post-tx only"},
            {"item": "joined_DRMref_x_GEO", "n": int(len(cells)), "note": "barcode match"},
            {"item": "post_tx_patients", "n": int(post["patient"].nunique()), "note": "P02–P04, P06, P07, P09–P15"},
            {"item": "MPR_patients", "n": int((post_s["mpr_group"] == "MPR").sum()), "note": "pCR counted as MPR"},
            {"item": "NMPR_patients", "n": int((post_s["mpr_group"] == "NMPR").sum()), "note": ""},
            {"item": "malignant_cells_post", "n": int(len(mal)), "note": "DRMref Malignant cells"},
            {
                "item": "cldn4_high_malignant",
                "n": n_high,
                "note": "CLDN4-only: log1p(CP10k) >= malignant median",
            },
            {
                "item": "dual_high_companion_not_used",
                "n": int(n_dual) if np.isfinite(n_dual) else 0,
                "note": "TACSTD2 and CLDN4 >= median; companion count only",
            },
            {"item": "cldn4_high_MPR_cells", "n": int(len(high_mpr)), "note": "not a patient replicate"},
            {"item": "cldn4_high_NMPR_cells", "n": int(len(high_nmpr)), "note": "not a patient replicate"},
            {"item": "TNK_cells_post", "n": int(len(tnk)), "note": "CD4+T + CD8+T + NK"},
            {"item": "TNK_MPR_cells", "n": int(len(tnk_mpr)), "note": "not a patient replicate"},
            {"item": "TNK_NMPR_cells", "n": int(len(tnk_nmpr)), "note": "not a patient replicate"},
            {"item": "potential_ligands_pooled", "n": int(len(potential)), "note": "CLDN4-high ligand + T/NK receptor"},
            {"item": "potential_ligands_MPR", "n": int(len(pot_mpr)), "note": "combinatorial MPR senders/receivers"},
            {"item": "potential_ligands_NMPR", "n": int(len(pot_nmpr)), "note": "combinatorial NMPR senders/receivers"},
            {"item": "background_genes_restricted", "n": int(len(background)), "note": "panel ∩ prior ∩ T/NK-expressed"},
            {"item": "a_priori_ifn_in_prior", "n": int(len(sets["a_priori_ifn"])), "note": "IFN genes present in v2 prior"},
            {
                "item": "a_priori_cyto_in_prior",
                "n": int(len(sets["a_priori_cytotoxicity"])),
                "note": "cytotoxicity genes present in v2 prior",
            },
            {"item": "empirical_NMPR_up_genes", "n": int(len(emp_up)), "note": f"patient Wilcoxon p<{EMPIRICAL_P}"},
            {"item": "empirical_NMPR_down_genes", "n": int(len(emp_down)), "note": f"patient Wilcoxon p<{EMPIRICAL_P}"},
        ]
    )
    n_tbl.to_csv(OUT / "n_table.tsv", sep="\t", index=False)

    tests_df = pd.DataFrame(tests)
    tests_df.to_csv(OUT / "patient_level_tests.tsv", sep="\t", index=False)

    def top_list(setting: str, geneset: str, k: int = 10) -> list[dict]:
        if top.empty:
            return []
        sub = top[(top["setting"] == setting) & (top["geneset"] == geneset)].head(k)
        return [
            {
                "rank": int(r.rank),
                "ligand": r.ligand,
                "pearson": float(r.pearson),
                "pearson_p": float(r.pearson_p),
                "auroc": float(r.auroc),
            }
            for r in sub.itertuples()
        ]

    primary_tests = {
        row["test"]: {k: row[k] for k in row if k != "test"}
        for row in tests_df.to_dict(orient="records")
        if row["test"]
        in {
            "NMPR_vs_MPR_tnk_cyto",
            "NMPR_vs_MPR_tnk_ifn",
            "NMPR_vs_MPR_tnk_ifn_minus_cyto",
            "spearman_mal_CLDN4__tnk_cyto",
            "spearman_mal_CLDN4__tnk_ifn",
        }
    }

    summary = {
        "dataset": "GSE207422",
        "sender_definition": "CLDN4-only malignant high (log1p CP10k >= malignant median)",
        "skipped": {
            "GSE253013": "9.3 GB RDS; no MPR/NMPR labels; not used",
            "full_nichenetr": "R/nichenetr not installed; v2 ligand-target prior used in Python",
            "TACSTD2_in_sender_call": "CLDN4 only; TACSTD2 is a companion column",
        },
        "n": {
            "MPR_patients": int((post_s["mpr_group"] == "MPR").sum()),
            "NMPR_patients": int((post_s["mpr_group"] == "NMPR").sum()),
            "malignant_cells": int(len(mal)),
            "cldn4_high_malignant": n_high,
            "dual_high_companion_not_used": int(n_dual) if np.isfinite(n_dual) else None,
            "tnk_cells": int(len(tnk)),
            "potential_ligands": int(len(potential)),
            "background_genes": int(len(background)),
        },
        "cutoffs": {
            "median_log1p_cp10k_CLDN4": med_c,
            "median_log1p_cp10k_TACSTD2_companion": med_t,
        },
        "hypothesis_patient_tests": primary_tests,
        "empirical_sets": {
            "up": sorted(emp_up),
            "down": sorted(emp_down),
            "p_cutoff": EMPIRICAL_P,
        },
        "top_ligands": {
            "cldn4_high_to_all_TNK__ifn": top_list("cldn4_high_to_all_TNK", "a_priori_ifn"),
            "cldn4_high_to_all_TNK__cytotoxicity": top_list(
                "cldn4_high_to_all_TNK", "a_priori_cytotoxicity"
            ),
            "NMPR_to_NMPR__ifn": top_list("cldn4_high_NMPR_to_NMPR_TNK", "a_priori_ifn"),
            "MPR_to_MPR__ifn": top_list("cldn4_high_MPR_to_MPR_TNK", "a_priori_ifn"),
            "NMPR_to_NMPR__cytotoxicity": top_list(
                "cldn4_high_NMPR_to_NMPR_TNK", "a_priori_cytotoxicity"
            ),
            "MPR_to_MPR__cytotoxicity": top_list(
                "cldn4_high_MPR_to_MPR_TNK", "a_priori_cytotoxicity"
            ),
        },
        "caveats": [
            "n=4 MPR vs n=8 NMPR; all patient-level p-values are underpowered",
            "NicheNet prior is unsigned regulatory potential, not signed up/down",
            "Background genes are panel-restricted (not the full receiver transcriptome)",
            "DRMref labels are marker-based, not Hu CopyKAT malignant IDs",
            "Cell counts are descriptive; inference is patient-level",
            "Sender is CLDN4-only; this does not replace the dual-high TACSTD2∩CLDN4 folder",
        ],
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    colors = {"MPR": "#2c7bb6", "NMPR": "#d7191c"}
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.6))
    for ax, metric, title in [
        (axes[0], "mal_CLDN4", "Malignant CLDN4"),
        (axes[1], "tnk_ifn", "T/NK IFN"),
        (axes[2], "tnk_cyto", "T/NK cytotoxicity"),
    ]:
        rng = np.random.default_rng(0)
        for grp, sub in post_s.groupby("mpr_group"):
            ax.scatter(
                rng.uniform(-0.08, 0.08, len(sub)) + (0 if grp == "MPR" else 1),
                sub[metric],
                c=colors[grp],
                s=40,
                label=f"{grp} n={len(sub)}",
                zorder=3,
            )
        ax.set_xticks([0, 1], ["MPR", "NMPR"])
        ax.set_title(title)
        ax.set_ylabel("mean log1p(CP10k)")
        ax.legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "fig1_patient_programs.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(4.4, 4.0))
    for grp, sub in post_s.groupby("mpr_group"):
        ax.scatter(sub["mal_CLDN4"], sub["tnk_ifn"], c=colors[grp], s=50, label=f"{grp} n={len(sub)}")
        for r in sub.itertuples():
            ax.annotate(
                r.patient,
                (r.mal_CLDN4, r.tnk_ifn),
                fontsize=7,
                xytext=(3, 3),
                textcoords="offset points",
            )
    ax.set_xlabel("Malignant CLDN4 mean log1p(CP10k)")
    ax.set_ylabel("T/NK IFN mean log1p(CP10k)")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "fig2_cldn4_vs_tnk_ifn.png", dpi=160)
    plt.close(fig)

    if not top.empty:
        fig, axes = plt.subplots(1, 2, figsize=(11, 5.2))
        for ax, gs, title in [
            (axes[0], "a_priori_ifn", "Prior activity vs IFN set"),
            (axes[1], "a_priori_cytotoxicity", "Prior activity vs cytotoxicity set"),
        ]:
            sub = top[(top["setting"] == "cldn4_high_to_all_TNK") & (top["geneset"] == gs)].head(12)
            if sub.empty:
                ax.axis("off")
                continue
            ax.barh(sub["ligand"][::-1], sub["pearson"][::-1], color="#4C72B0")
            ax.set_xlabel("Pearson (prior vs gene-set membership)")
            ax.set_title(title)
        fig.tight_layout()
        fig.savefig(OUT / "fig3_ligand_activity_top.png", dpi=160)
        plt.close(fig)

        a = act_all[
            (act_all["geneset"] == "a_priori_ifn")
            & (act_all["setting"] == "cldn4_high_NMPR_to_NMPR_TNK")
        ][["ligand", "pearson"]].rename(columns={"pearson": "pearson_NMPR"})
        b = act_all[
            (act_all["geneset"] == "a_priori_ifn")
            & (act_all["setting"] == "cldn4_high_MPR_to_MPR_TNK")
        ][["ligand", "pearson"]].rename(columns={"pearson": "pearson_MPR"})
        both = a.merge(b, on="ligand", how="inner")
        if not both.empty:
            both.to_csv(OUT / "combinatorial_ifn_activity_MPR_vs_NMPR.tsv", sep="\t", index=False)
            fig, ax = plt.subplots(figsize=(5.2, 5.0))
            ax.scatter(both["pearson_MPR"], both["pearson_NMPR"], s=18, c="#555555")
            lim = [
                min(both["pearson_MPR"].min(), both["pearson_NMPR"].min()) - 0.02,
                max(both["pearson_MPR"].max(), both["pearson_NMPR"].max()) + 0.02,
            ]
            ax.plot(lim, lim, ls="--", c="0.7", lw=1)
            both = both.assign(delta=both["pearson_NMPR"] - both["pearson_MPR"])
            for r in both.reindex(both["delta"].abs().sort_values(ascending=False).index).head(8).itertuples():
                ax.annotate(r.ligand, (r.pearson_MPR, r.pearson_NMPR), fontsize=7)
            ax.set_xlabel("Ligand activity Pearson (MPR receivers)")
            ax.set_ylabel("Ligand activity Pearson (NMPR receivers)")
            ax.set_title("Combinatorial IFN: same prior, MPR vs NMPR T/NK")
            fig.tight_layout()
            fig.savefig(OUT / "fig4_combinatorial_MPR_vs_NMPR.png", dpi=160)
            plt.close(fig)

    print("wrote results to", OUT, flush=True)
    print(json.dumps({k: summary[k] for k in ("n", "hypothesis_patient_tests")}, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
