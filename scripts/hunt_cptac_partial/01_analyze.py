#!/usr/bin/env python3
"""CPTAC LUAD/LSCC: TACSTD2 protein/RNA vs xCell CD8/immune, partial Spearman.

Pre-specified question (user direction: negative):
  Does TACSTD2 protein or RNA anti-correlate with xCell CD8 / xCell immune
  score after controlling for any available tumor-purity proxy?

Cohorts are treatment-naive surgical CPTAC (Gillette Cell 2020 LUAD;
Satpathy Cell 2021 LSCC). No ICI labels exist — do not invent them.

Inputs:  data/hunt_cptac_partial/{LUAD,LSCC}/  (from 00_download.py)
Outputs: results/hunt_cptac_partial/{tables,figures}/
"""
from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

TACSTD2_ENSEMBL = "ENSG00000184292"
TACSTD2_SYMBOL = "TACSTD2"

PREDICTORS = ("TACSTD2_protein", "TACSTD2_RNA")
PRIMARY_TARGETS = ("xCell_T_cell_CD8+", "xCell_immune_score")
SECONDARY_TARGETS = ("CIBERSORT_T_cell_CD8+", "ESTIMATE_ImmuneScore")
ALL_TARGETS = PRIMARY_TARGETS + SECONDARY_TARGETS

# Yoshihara et al. 2013 Nat Commun (ESTIMATE tumor purity).
ESTIMATE_PURITY_A = 0.6049872018
ESTIMATE_PURITY_B = 0.0001467884

COHORTS = ("LUAD", "LSCC")


def find_row(index: pd.Index, *aliases: str) -> str | None:
    hits: list[str] = []
    for alias in aliases:
        for i in index:
            s = str(i)
            if s == alias or s.startswith(alias + ".") or s.upper() == alias.upper():
                hits.append(s)
    hits = list(dict.fromkeys(hits))
    if not hits:
        return None
    if len(hits) > 1:
        # Prefer exact symbol / exact Ensembl over versioned if both exist.
        exact = [h for h in hits if h in aliases]
        if len(exact) == 1:
            return exact[0]
        raise ValueError(f"multiple rows for {aliases}: {hits}")
    return hits[0]


def load_matrix(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", index_col=0)


def extract_gene(mat: pd.DataFrame, *aliases: str) -> tuple[pd.Series, str | None]:
    row = find_row(mat.index, *aliases)
    if row is None:
        return pd.Series(np.nan, index=mat.columns, name=aliases[0]), None
    s = pd.to_numeric(mat.loc[row], errors="coerce")
    s.name = aliases[0]
    return s, row


def estimate_tumor_purity(estimate_score: pd.Series) -> pd.Series:
    x = pd.to_numeric(estimate_score, errors="coerce")
    return pd.Series(np.cos(ESTIMATE_PURITY_A + ESTIMATE_PURITY_B * x), index=x.index)


def spearman_pair(x: pd.Series, y: pd.Series) -> dict:
    d = pd.concat([x, y], axis=1).dropna()
    d.columns = ["x", "y"]
    n = len(d)
    if n < 4 or d["x"].nunique() < 2 or d["y"].nunique() < 2:
        return {"n": int(n), "rho": np.nan, "p": np.nan}
    rho, p = stats.spearmanr(d["x"], d["y"])
    return {"n": int(n), "rho": float(rho), "p": float(p)}


def partial_spearman_formula(x: pd.Series, y: pd.Series, z: pd.Series) -> dict:
    """Closed-form partial Spearman: Pearson on ranks, then ρ_xy·z.

    p from Student-t with n-3 df. This is the pre-specified method.
    """
    d = pd.concat([x, y, z], axis=1).dropna()
    d.columns = ["x", "y", "z"]
    n = len(d)
    out = {
        "n": int(n),
        "rho": np.nan,
        "p": np.nan,
        "rho_xy": np.nan,
        "rho_xz": np.nan,
        "rho_yz": np.nan,
        "method": "partial_spearman_formula",
    }
    if n < 6 or d["x"].nunique() < 2 or d["y"].nunique() < 2 or d["z"].nunique() < 2:
        return out
    rxy, _ = stats.spearmanr(d["x"], d["y"])
    rxz, _ = stats.spearmanr(d["x"], d["z"])
    ryz, _ = stats.spearmanr(d["y"], d["z"])
    out["rho_xy"] = float(rxy)
    out["rho_xz"] = float(rxz)
    out["rho_yz"] = float(ryz)
    den = math.sqrt(max(0.0, (1.0 - rxz**2) * (1.0 - ryz**2)))
    if den < 1e-12:
        return out
    rho = (rxy - rxz * ryz) / den
    rho = float(np.clip(rho, -1.0, 1.0))
    out["rho"] = rho
    if abs(rho) >= 1.0 - 1e-15:
        out["p"] = 0.0
        return out
    t = rho * math.sqrt((n - 3) / (1.0 - rho**2))
    out["p"] = float(2.0 * stats.t.sf(abs(t), df=n - 3))
    return out


def partial_spearman_residual(x: pd.Series, y: pd.Series, z: pd.Series) -> dict:
    """Rank-residual Spearman (sensitivity; same estimator family as prior LUAD slice)."""
    d = pd.concat([x, y, z], axis=1).dropna()
    d.columns = ["x", "y", "z"]
    n = len(d)
    if n < 6 or d["x"].nunique() < 2 or d["y"].nunique() < 2 or d["z"].nunique() < 2:
        return {"n": int(n), "rho": np.nan, "p": np.nan, "method": "rank_residual"}
    rx, ry, rz = d["x"].rank(), d["y"].rank(), d["z"].rank()
    bx = np.polyfit(rz, rx, 1)
    by = np.polyfit(rz, ry, 1)
    ex = rx - (bx[0] * rz + bx[1])
    ey = ry - (by[0] * rz + by[1])
    rho, p = stats.spearmanr(ex, ey)
    return {"n": int(n), "rho": float(rho), "p": float(p), "method": "rank_residual"}


def fdr(pvals: list[float]) -> list[float]:
    arr = np.asarray(pvals, dtype=float)
    out = np.full(arr.shape, np.nan)
    mask = np.isfinite(arr)
    if mask.sum() == 0:
        return out.tolist()
    out[mask] = multipletests(arr[mask], method="fdr_bh")[1]
    return [float(v) if np.isfinite(v) else np.nan for v in out]


def fisher_z_meta(rows: list[dict]) -> dict:
    """Inverse-variance meta of Fisher-z Spearman (independent cohorts)."""
    zs, ws, ns = [], [], []
    for r in rows:
        n, rho = r.get("n"), r.get("rho")
        if n is None or rho is None or not np.isfinite(rho) or n < 6 or abs(rho) >= 1:
            continue
        z = math.atanh(float(np.clip(rho, -0.999999, 0.999999)))
        w = n - 3
        zs.append(z)
        ws.append(w)
        ns.append(n)
    if not zs:
        return {"n_cohorts": 0, "n_total": 0, "rho": np.nan, "p": np.nan}
    zbar = float(np.average(zs, weights=ws))
    se = math.sqrt(1.0 / sum(ws))
    p = float(2.0 * stats.norm.sf(abs(zbar / se)))
    return {
        "n_cohorts": int(len(zs)),
        "n_total": int(sum(ns)),
        "rho": float(math.tanh(zbar)),
        "p": p,
    }


def load_cohort(data: Path, cohort: str) -> dict:
    ddir = data / cohort
    rna_path = ddir / f"{cohort}_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt"
    prot_path = ddir / (
        f"{cohort}_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt"
    )
    pheno_path = ddir / f"{cohort}_phenotype.txt"
    rna = load_matrix(rna_path)
    prot = load_matrix(prot_path)
    pheno = pd.read_csv(pheno_path, sep="\t", index_col=0)

    tac_rna, rna_id = extract_gene(rna, TACSTD2_ENSEMBL, TACSTD2_SYMBOL)
    tac_prot, prot_id = extract_gene(prot, TACSTD2_ENSEMBL, TACSTD2_SYMBOL)

    samples = sorted(set(rna.columns) & set(prot.columns) & set(pheno.index))
    core = pd.DataFrame(index=samples)
    core["TACSTD2_RNA"] = tac_rna.reindex(samples)
    core["TACSTD2_protein"] = tac_prot.reindex(samples)

    target_map = {}
    for tgt in ALL_TARGETS:
        if tgt in pheno.columns:
            core[tgt] = pd.to_numeric(pheno.loc[samples, tgt], errors="coerce")
            target_map[tgt] = tgt
        else:
            core[tgt] = np.nan
            target_map[tgt] = None

    purity: dict[str, pd.Series] = {}
    notes: dict[str, str] = {}
    for col in ("WES_purity", "WGS_purity"):
        if col in pheno.columns:
            s = pd.to_numeric(pheno.loc[samples, col], errors="coerce")
            purity[col] = s
            core[col] = s
            notes[col] = "DNA-based purity column in freeze phenotype"

    estimate_col = next(
        (
            c
            for c in ("ESTIMATE_ESTIMATEScore", "ESTIMATEScore", "ESTIMATE_Score")
            if c in pheno.columns
        ),
        None,
    )
    yoshihara_ok = False
    if estimate_col is not None:
        es = pd.to_numeric(pheno.loc[samples, estimate_col], errors="coerce")
        core["ESTIMATEScore"] = es
        # RNA impurity axis. Partial Spearman is identical for ESTIMATEScore
        # and -ESTIMATEScore (both xz and yz flip). Higher score = more
        # immune+stroma = lower purity. Not independent of xCell immune.
        purity["ESTIMATEScore"] = es
        notes["ESTIMATEScore"] = (
            f"RNA ESTIMATE combined score ({estimate_col}). Impurity axis, "
            "not a 0–1 purity. Not independent of xCell/ESTIMATE immune."
        )
        ep = estimate_tumor_purity(es)
        wrap = (math.pi - ESTIMATE_PURITY_A) / ESTIMATE_PURITY_B
        n_wrap = int((es > wrap).sum())
        n_neg = int((ep < 0).sum())
        core["ESTIMATE_TumorPurity_Yoshihara"] = ep
        # Only treat cosine as a usable z if it stays in (0,1] and does not wrap.
        yoshihara_ok = n_wrap == 0 and n_neg == 0 and ep.notna().sum() >= 6
        notes["_yoshihara_qc"] = (
            f"Yoshihara cos(a+b*ESTIMATEScore) on {estimate_col}: "
            f"n_wrap(score>{wrap:.1f})={n_wrap}, n_negative_purity={n_neg}, "
            f"usable={yoshihara_ok}. This freeze's ESTIMATE scores are too "
            "large for the 2013 cosine (non-monotonic / unphysical)."
        )
        if yoshihara_ok:
            purity["ESTIMATE_TumorPurity"] = ep
            notes["ESTIMATE_TumorPurity"] = (
                "Yoshihara 2013 cosine of ESTIMATEScore; RNA-derived."
            )
    elif "ESTIMATE_TumorPurity" in pheno.columns:
        s = pd.to_numeric(pheno.loc[samples, "ESTIMATE_TumorPurity"], errors="coerce")
        purity["ESTIMATE_TumorPurity"] = s
        core["ESTIMATE_TumorPurity"] = s
        notes["ESTIMATE_TumorPurity"] = "Precomputed ESTIMATE_TumorPurity in phenotype"

    extra_purity = [
        c
        for c in pheno.columns
        if "purity" in str(c).lower() and c not in purity and c not in core.columns
    ]
    for col in extra_purity:
        s = pd.to_numeric(pheno.loc[samples, col], errors="coerce")
        if s.notna().sum() >= 6 and s.nunique(dropna=True) >= 2:
            purity[col] = s
            core[col] = s
            notes[col] = "Additional phenotype column matching *purity*"

    # Keep a few columns that help interpret confounding (not used as z).
    for extra in ("ESTIMATE_StromalScore", "xCell_stroma_score", "xCell_microenvironment_score"):
        if extra in pheno.columns:
            core[extra] = pd.to_numeric(pheno.loc[samples, extra], errors="coerce")

    return {
        "cohort": cohort,
        "core": core,
        "purity": purity,
        "purity_notes": notes,
        "target_map": target_map,
        "rna_id": rna_id,
        "prot_id": prot_id,
        "n_rna_genes": int(rna.shape[0]),
        "n_prot_genes": int(prot.shape[0]),
        "n_pheno_cols": int(pheno.shape[1]),
        "n_intersect": int(len(samples)),
        "yoshihara_qc": notes.get("_yoshihara_qc"),
    }


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--data", default="data/hunt_cptac_partial")
    p.add_argument("--outdir", default="results/hunt_cptac_partial")
    args = p.parse_args()
    data = Path(args.data)
    out = Path(args.outdir)
    tables = out / "tables"
    figs = out / "figures"
    tables.mkdir(parents=True, exist_ok=True)
    figs.mkdir(parents=True, exist_ok=True)

    loaded = [load_cohort(data, c) for c in COHORTS]

    coverage_rows = []
    purity_avail = []
    yoshihara_rows = []
    sample_frames = []
    unadj_rows = []
    partial_rows = []
    pred_pur_rows = []
    tgt_pur_rows = []

    for pack in loaded:
        cohort = pack["cohort"]
        core: pd.DataFrame = pack["core"]
        purity: dict[str, pd.Series] = pack["purity"]
        core = core.copy()
        core.insert(0, "cohort", cohort)
        core.insert(1, "sample_id", core.index)
        sample_frames.append(core.reset_index(drop=True))

        if pack.get("yoshihara_qc"):
            yoshihara_rows.append({"cohort": cohort, "note": pack["yoshihara_qc"]})
        coverage_rows.append(
            {
                "cohort": cohort,
                "n_intersect_rna_protein_pheno": pack["n_intersect"],
                "TACSTD2_RNA_id": pack["rna_id"],
                "TACSTD2_protein_id": pack["prot_id"],
                "n_TACSTD2_RNA": int(core["TACSTD2_RNA"].notna().sum()),
                "n_TACSTD2_protein": int(core["TACSTD2_protein"].notna().sum()),
                "n_rna_genes": pack["n_rna_genes"],
                "n_prot_genes": pack["n_prot_genes"],
                "n_pheno_cols": pack["n_pheno_cols"],
                **{f"n_{t}": int(core[t].notna().sum()) for t in ALL_TARGETS},
                "targets_missing": ",".join(
                    t for t, m in pack["target_map"].items() if m is None
                ),
            }
        )
        for name, s in purity.items():
            purity_avail.append(
                {
                    "cohort": cohort,
                    "proxy": name,
                    "n": int(s.notna().sum()),
                    "n_missing": int(s.isna().sum()),
                    "min": float(s.min()) if s.notna().any() else np.nan,
                    "median": float(s.median()) if s.notna().any() else np.nan,
                    "max": float(s.max()) if s.notna().any() else np.nan,
                    "independence_note": pack["purity_notes"].get(name, ""),
                }
            )

        for pred in PREDICTORS:
            for tgt in ALL_TARGETS:
                r = spearman_pair(core[pred], core[tgt])
                r.update(
                    {
                        "cohort": cohort,
                        "predictor": pred,
                        "target": tgt,
                        "adjust": "none",
                        "primary": tgt in PRIMARY_TARGETS,
                    }
                )
                unadj_rows.append(r)

            for pname, z in purity.items():
                pr = spearman_pair(core[pred], z)
                pr.update({"cohort": cohort, "predictor": pred, "purity_proxy": pname})
                pred_pur_rows.append(pr)

        for tgt in ALL_TARGETS:
            for pname, z in purity.items():
                tr = spearman_pair(core[tgt], z)
                tr.update({"cohort": cohort, "target": tgt, "purity_proxy": pname})
                tgt_pur_rows.append(tr)

        for pred in PREDICTORS:
            for tgt in ALL_TARGETS:
                for pname, z in purity.items():
                    form = partial_spearman_formula(core[pred], core[tgt], z)
                    form.update(
                        {
                            "cohort": cohort,
                            "predictor": pred,
                            "target": tgt,
                            "adjust": pname,
                            "primary": tgt in PRIMARY_TARGETS,
                        }
                    )
                    resid = partial_spearman_residual(core[pred], core[tgt], z)
                    form["rho_rank_residual"] = resid["rho"]
                    form["p_rank_residual"] = resid["p"]
                    partial_rows.append(form)

    cov_df = pd.DataFrame(coverage_rows)
    cov_df.to_csv(tables / "gene_coverage.tsv", sep="\t", index=False)
    pd.DataFrame(purity_avail).to_csv(tables / "purity_availability.tsv", sep="\t", index=False)
    if yoshihara_rows:
        pd.DataFrame(yoshihara_rows).to_csv(
            tables / "yoshihara_purity_qc.tsv", sep="\t", index=False
        )
    pd.concat(sample_frames, ignore_index=True).to_csv(
        tables / "sample_table.tsv", sep="\t", index=False
    )

    unadj = pd.DataFrame(unadj_rows)
    # FDR: pre-specified 8 primary unadjusted tests (2 cohort × 2 pred × 2 xCell).
    prim_mask = unadj["primary"]
    unadj["q_primary8"] = np.nan
    unadj.loc[prim_mask, "q_primary8"] = fdr(unadj.loc[prim_mask, "p"].tolist())
    # Also BH within cohort × predictor across the 4 immune targets (unadjusted).
    unadj["q_within_cohort_predictor"] = np.nan
    for _, idx in unadj.groupby(["cohort", "predictor"]).groups.items():
        unadj.loc[idx, "q_within_cohort_predictor"] = fdr(unadj.loc[idx, "p"].tolist())
    unadj = unadj.sort_values(["cohort", "predictor", "primary", "p"], ascending=[True, True, False, True])
    unadj.to_csv(tables / "spearman_unadjusted.tsv", sep="\t", index=False)

    pred_pur = pd.DataFrame(pred_pur_rows)
    pred_pur.to_csv(tables / "predictor_vs_purity.tsv", sep="\t", index=False)
    tgt_pur = pd.DataFrame(tgt_pur_rows)
    tgt_pur.to_csv(tables / "target_vs_purity.tsv", sep="\t", index=False)

    part = pd.DataFrame(partial_rows)
    part["q_primary_per_adjust"] = np.nan
    prim_part = part["primary"]
    for adj, idx in part.loc[prim_part].groupby("adjust").groups.items():
        part.loc[idx, "q_primary_per_adjust"] = fdr(part.loc[idx, "p"].tolist())
    part = part.sort_values(["cohort", "predictor", "target", "adjust"])
    part.to_csv(tables / "partial_spearman.tsv", sep="\t", index=False)

    # Focus table: primary pairs, unadjusted + each proxy.
    focus_un = unadj[unadj["primary"]].copy()
    focus_un["method"] = "unadjusted_spearman"
    focus_un["rho_xy"] = focus_un["rho"]
    focus_un["rho_xz"] = np.nan
    focus_un["rho_yz"] = np.nan
    focus_un["rho_rank_residual"] = np.nan
    focus_un["p_rank_residual"] = np.nan
    cols = [
        "cohort",
        "predictor",
        "target",
        "adjust",
        "method",
        "n",
        "rho",
        "p",
        "rho_xy",
        "rho_xz",
        "rho_yz",
        "rho_rank_residual",
        "p_rank_residual",
        "q_primary8",
        "q_primary_per_adjust",
    ]
    focus_part = part[part["primary"]].copy()
    focus_part["method"] = "partial_spearman_formula"
    focus_part["q_primary8"] = np.nan
    focus = pd.concat(
        [focus_un.reindex(columns=cols), focus_part.reindex(columns=cols)],
        ignore_index=True,
    )
    focus["direction_matches_user_negative"] = focus["rho"] < 0
    focus.to_csv(tables / "primary_focus.tsv", sep="\t", index=False)

    # Meta-analysis across LUAD+LSCC for primary pairs.
    meta_rows = []
    for pred in PREDICTORS:
        for tgt in PRIMARY_TARGETS:
            sub = unadj[(unadj.predictor == pred) & (unadj.target == tgt)]
            m = fisher_z_meta(sub.to_dict("records"))
            m.update(
                {
                    "predictor": pred,
                    "target": tgt,
                    "adjust": "none",
                    "cohorts": ",".join(sub["cohort"].tolist()),
                }
            )
            meta_rows.append(m)
            for adj in sorted(part["adjust"].dropna().unique()):
                subp = part[
                    (part.predictor == pred) & (part.target == tgt) & (part.adjust == adj)
                ]
                if subp.empty:
                    continue
                mp = fisher_z_meta(subp.to_dict("records"))
                mp.update(
                    {
                        "predictor": pred,
                        "target": tgt,
                        "adjust": adj,
                        "cohorts": ",".join(subp["cohort"].tolist()),
                    }
                )
                meta_rows.append(mp)
    meta = pd.DataFrame(meta_rows)
    meta.to_csv(tables / "meta_luad_lscc.tsv", sep="\t", index=False)

    # ---- Figures ----
    _plot_forest(focus, figs / "fig1_forest_primary_partial.png")
    for pack in loaded:
        _plot_scatter(
            pack["core"],
            pack["cohort"],
            figs / f"fig_scatter_{pack['cohort'].lower()}_protein_xcell_cd8.png",
        )
    _plot_purity_confound(loaded, figs / "fig2_purity_vs_xcell.png")

    key = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "question": (
            "TACSTD2 protein/RNA vs xCell CD8 / xCell immune after any purity proxy"
        ),
        "user_direction": "negative",
        "prior_claim": "LUAD TACSTD2 protein vs xCell CD8 Spearman rho ≈ -0.29",
        "cohorts": COHORTS,
        "predictors": list(PREDICTORS),
        "primary_targets": list(PRIMARY_TARGETS),
        "secondary_targets": list(SECONDARY_TARGETS),
        "partial_method": (
            "Closed-form partial Spearman (Pearson on ranks); p from t with n-3 df. "
            "Rank-residual Spearman stored as sensitivity."
        ),
        "purity_proxies_used": sorted({r["proxy"] for r in purity_avail}),
        "yoshihara_qc": yoshihara_rows,
        "coverage": coverage_rows,
        "purity_availability": purity_avail,
        "primary_unadjusted": focus_un[
            ["cohort", "predictor", "target", "n", "rho", "p", "q_primary8"]
        ].to_dict("records"),
        "primary_partial": focus_part[
            ["cohort", "predictor", "target", "adjust", "n", "rho", "p", "q_primary_per_adjust"]
        ].to_dict("records"),
        "meta": meta.to_dict("records"),
        "honesty": {
            "no_ici_labels": True,
            "estimate_purity_not_independent_of_xcell": True,
            "tacstd2_rna_vs_xcell_same_rna_layer": True,
            "wes_wgs_purity_are_preferred_covariates": True,
        },
    }
    (tables / "key_stats.json").write_text(json.dumps(key, indent=2, default=_json_default) + "\n")
    print(out)
    return 0


def _json_default(o):
    if isinstance(o, (np.floating,)):
        return float(o) if np.isfinite(o) else None
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, float) and not math.isfinite(o):
        return None
    raise TypeError(type(o))


def _plot_forest(focus: pd.DataFrame, path: Path) -> None:
    # One panel per primary pair; points = unadjusted + each adjust.
    pairs = [
        (c, pred, tgt)
        for c in COHORTS
        for pred in PREDICTORS
        for tgt in PRIMARY_TARGETS
    ]
    fig, axes = plt.subplots(len(COHORTS), 4, figsize=(14, 6.5), sharex=True)
    if axes.ndim == 1:
        axes = np.array([axes])
    short_tgt = {
        "xCell_T_cell_CD8+": "xCell CD8",
        "xCell_immune_score": "xCell immune",
    }
    short_pred = {"TACSTD2_protein": "protein", "TACSTD2_RNA": "RNA"}
    for i, cohort in enumerate(COHORTS):
        col = 0
        for pred in PREDICTORS:
            for tgt in PRIMARY_TARGETS:
                ax = axes[i, col]
                sub = focus[(focus.cohort == cohort) & (focus.predictor == pred) & (focus.target == tgt)]
                sub = sub.copy()
                # order: none first, then DNA proxies, then ESTIMATE
                order = ["none", "WES_purity", "WGS_purity", "ESTIMATEScore", "ESTIMATE_TumorPurity"]
                extra = [a for a in sub["adjust"].unique() if a not in order]
                order = [a for a in order if a in set(sub["adjust"])] + extra
                sub["ord"] = sub["adjust"].map({a: k for k, a in enumerate(order)})
                sub = sub.sort_values("ord")
                y = np.arange(len(sub))
                colors = ["#222222" if a == "none" else "#1f77b4" if "ESTIMATE" not in str(a) else "#d62728" for a in sub["adjust"]]
                ax.axvline(0, color="#888888", lw=0.8)
                ax.scatter(sub["rho"], y, c=colors, s=36, zorder=3)
                ax.set_yticks(y)
                ax.set_yticklabels([str(a) for a in sub["adjust"]], fontsize=7)
                ax.set_title(f"{cohort} {short_pred[pred]} vs {short_tgt[tgt]}", fontsize=8)
                ax.set_xlim(-0.7, 0.5)
                ax.grid(axis="x", ls=":", alpha=0.4)
                col += 1
        axes[i, 0].set_ylabel(cohort)
    axes[-1, 0].set_xlabel("Spearman ρ (unadjusted or partial)")
    fig.suptitle(
        "TACSTD2 vs xCell CD8/immune: unadjusted vs partial Spearman after purity",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _plot_scatter(core: pd.DataFrame, cohort: str, path: Path) -> None:
    d = core[["TACSTD2_protein", "xCell_T_cell_CD8+", "WES_purity"]].copy()
    if "WES_purity" not in core.columns:
        d["WES_purity"] = np.nan
    d = d.dropna(subset=["TACSTD2_protein", "xCell_T_cell_CD8+"])
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    c = d["WES_purity"] if d["WES_purity"].notna().any() else None
    sc = ax.scatter(
        d["TACSTD2_protein"],
        d["xCell_T_cell_CD8+"],
        c=c,
        cmap="viridis",
        s=28,
        alpha=0.85,
        edgecolors="none",
    )
    if c is not None and c.notna().any():
        cb = fig.colorbar(sc, ax=ax, shrink=0.8)
        cb.set_label("WES_purity")
    r = spearman_pair(d["TACSTD2_protein"], d["xCell_T_cell_CD8+"])
    ax.set_xlabel("TACSTD2 protein (TMT log2)")
    ax.set_ylabel("xCell CD8")
    ax.set_title(f"{cohort}  ρ={r['rho']:.3f}  p={r['p']:.3g}  n={r['n']}")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _plot_purity_confound(loaded: list[dict], path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.0))
    for ax, pack in zip(axes, loaded):
        core = pack["core"]
        if "WES_purity" not in core.columns:
            ax.set_title(f"{pack['cohort']}: no WES_purity")
            continue
        d = core[["WES_purity", "xCell_T_cell_CD8+", "xCell_immune_score"]].dropna()
        ax.scatter(d["WES_purity"], d["xCell_T_cell_CD8+"], s=18, alpha=0.7, label="xCell CD8")
        ax.scatter(
            d["WES_purity"],
            d["xCell_immune_score"],
            s=18,
            alpha=0.7,
            label="xCell immune",
        )
        r1 = spearman_pair(d["WES_purity"], d["xCell_T_cell_CD8+"])
        r2 = spearman_pair(d["WES_purity"], d["xCell_immune_score"])
        ax.set_xlabel("WES_purity")
        ax.set_ylabel("xCell score")
        ax.set_title(
            f"{pack['cohort']}  CD8 ρ={r1['rho']:.2f}  immune ρ={r2['rho']:.2f}"
        )
        ax.legend(fontsize=7)
    fig.suptitle("Purity vs xCell (why partial correlation is required)", fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    raise SystemExit(main())
