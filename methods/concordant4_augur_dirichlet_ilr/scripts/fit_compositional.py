#!/usr/bin/env python3
"""ILR, CLR, and Dirichlet-multinomial models on the locked 65 units."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parent))

from compositional import (  # noqa: E402
    BALANCE_NAMES_3,
    BALANCE_NAMES_4,
    SBP_TNK_MYE_REST,
    SBP_WITHIN_IMMUNE,
    design_matrix,
    dl_fisher_z,
    finite_difference_proportions,
    fit_dirichlet_multinomial,
    ilr_from_sbp,
    ols_hc1,
    permute_within_cohort,
    proportions_with_pseudocount,
    spearman_rho,
    student_p,
    within_cohort_z,
)

ROOT = Path(__file__).resolve().parents[1]
TAB = ROOT / "results" / "tables"
TAB.mkdir(parents=True, exist_ok=True)

REF = "GSE131907"
N_PERM_ILR = 1999
N_PERM_DM = 299
SEED = 1
PARTS4 = ["T", "NK", "myeloid", "Rest"]
COUNT4 = ["n_T", "n_NK", "n_myeloid", "n_rest"]
PARTS3 = ["T", "NK", "myeloid"]
COUNT3 = ["n_T", "n_NK", "n_myeloid"]


def rank_biserial(high: np.ndarray, low: np.ndarray) -> dict:
    high = np.asarray(high, dtype=float)
    low = np.asarray(low, dtype=float)
    res = mannwhitneyu(high, low, alternative="two-sided", method="asymptotic")
    u = float(res.statistic)
    n4, n1 = len(high), len(low)
    r = 2 * u / (n4 * n1) - 1
    return {"r": r, "p": float(res.pvalue), "n_q1": n1, "n_q4": n4}


def wilks_lambda(y: np.ndarray, x: np.ndarray, drop: int) -> float:
    beta = np.linalg.lstsq(x, y, rcond=None)[0]
    resid = y - x @ beta
    e = resid.T @ resid
    x_red = np.delete(x, drop, axis=1)
    beta_r = np.linalg.lstsq(x_red, y, rcond=None)[0]
    resid_r = y - x_red @ beta_r
    h_plus_e = resid_r.T @ resid_r
    sign_e, logdet_e = np.linalg.slogdet(e)
    sign_h, logdet_h = np.linalg.slogdet(h_plus_e)
    if sign_e <= 0 or sign_h <= 0:
        return 1.0
    return float(np.exp(logdet_e - logdet_h))


def fit_ilr_block(
    name: str,
    balances: np.ndarray,
    balance_names: list[str],
    cohorts: np.ndarray,
    z: np.ndarray,
    cldn4: np.ndarray,
    quartile: np.ndarray,
    rng: np.random.Generator,
) -> tuple[list[dict], list[dict], dict]:
    x, coef_names = design_matrix(cohorts, z, REF)
    cldn_i = coef_names.index("cldn4")
    rows = []
    for j, bname in enumerate(balance_names):
        fit = ols_hc1(balances[:, j], x)
        beta = float(fit["beta"][cldn_i])
        se = float(fit["se"][cldn_i])
        null = np.empty(N_PERM_ILR)
        for p_i in range(N_PERM_ILR):
            zp = permute_within_cohort(z, cohorts, rng)
            xp, _ = design_matrix(cohorts, zp, REF)
            null[p_i] = float(ols_hc1(balances[:, j], xp)["beta"][cldn_i])
        perm_p = (1 + np.sum(np.abs(null) >= abs(beta))) / (1 + N_PERM_ILR)
        rows.append(
            {
                "simplex": name,
                "balance": bname,
                "model": "continuous_z",
                "n": int(fit["n"]),
                "beta_per_sd": beta,
                "se_hc1": se,
                "p_hc1": student_p(beta, se, fit["dof"]),
                "p_perm": float(perm_p),
                "null_sd": float(null.std(ddof=1)),
            }
        )
    # joint Wilks on all balances
    obs = wilks_lambda(balances, x, cldn_i)
    null_w = np.empty(N_PERM_ILR)
    for p_i in range(N_PERM_ILR):
        zp = permute_within_cohort(z, cohorts, rng)
        xp, _ = design_matrix(cohorts, zp, REF)
        null_w[p_i] = wilks_lambda(balances, xp, cldn_i)
    joint = {
        "simplex": name,
        "stat": "wilks",
        "lambda": obs,
        "p_perm": float((1 + np.sum(null_w <= obs)) / (1 + N_PERM_ILR)),
        "n": int(len(z)),
        "n_balances": len(balance_names),
    }

    # cohort Spearmans of each balance vs raw CLDN4 %
    spearman_rows = []
    for j, bname in enumerate(balance_names):
        rhos, ns = [], []
        for coh in sorted(set(cohorts.tolist())):
            m = cohorts == coh
            rho, p = spearman_rho(cldn4[m], balances[m, j])
            spearman_rows.append(
                {
                    "simplex": name,
                    "balance": bname,
                    "dataset": coh,
                    "n": int(m.sum()),
                    "rho": rho,
                    "p": p,
                }
            )
            rhos.append(rho)
            ns.append(int(m.sum()))
        pooled = dl_fisher_z(np.array(rhos), np.array(ns, dtype=float))
        spearman_rows.append(
            {
                "simplex": name,
                "balance": bname,
                "dataset": "DL_meta",
                "n": pooled["N"],
                "rho": pooled["rho"],
                "p": pooled["p"],
                "I2": pooled["I2"],
                "ci_lo": pooled["ci_lo"],
                "ci_hi": pooled["ci_hi"],
                "k": pooled["k"],
            }
        )

    # Q4 vs Q1 on the stacked within-cohort labels
    q_rows = []
    q4 = quartile == "Q4"
    q1 = quartile == "Q1"
    for j, bname in enumerate(balance_names):
        rb = rank_biserial(balances[q4, j], balances[q1, j])
        # cohort-adjusted binary OLS on the Q1/Q4 subset
        sub_c = cohorts[q4 | q1]
        sub_y = balances[q4 | q1, j]
        sub_bin = (quartile[q4 | q1] == "Q4").astype(float)
        xb, names_b = design_matrix(sub_c, sub_bin, REF)
        fitb = ols_hc1(sub_y, xb)
        bi = names_b.index("cldn4")
        beta = float(fitb["beta"][bi])
        se = float(fitb["se"][bi])
        null = np.empty(N_PERM_ILR)
        for p_i in range(N_PERM_ILR):
            bp = permute_within_cohort(sub_bin, sub_c, rng)
            xp, _ = design_matrix(sub_c, bp, REF)
            null[p_i] = float(ols_hc1(sub_y, xp)["beta"][-1])
        q_rows.append(
            {
                "simplex": name,
                "balance": bname,
                "n_q1": rb["n_q1"],
                "n_q4": rb["n_q4"],
                "rank_biserial_q4_vs_q1": rb["r"],
                "p_mw": rb["p"],
                "beta_q4": beta,
                "se_hc1": se,
                "p_hc1": student_p(beta, se, fitb["dof"]),
                "p_perm": float((1 + np.sum(np.abs(null) >= abs(beta))) / (1 + N_PERM_ILR)),
            }
        )
    return rows, spearman_rows, {"joint": joint, "q4q1": q_rows}


def clr_rows(props: np.ndarray, names: list[str], cohorts, z, rng) -> list[dict]:
    logp = np.log(props)
    clr = logp - logp.mean(axis=1, keepdims=True)
    x, _ = design_matrix(cohorts, z, REF)
    out = []
    for j, name in enumerate(names):
        fit = ols_hc1(clr[:, j], x)
        beta = float(fit["beta"][-1])
        se = float(fit["se"][-1])
        null = np.empty(N_PERM_ILR)
        for p_i in range(N_PERM_ILR):
            zp = permute_within_cohort(z, cohorts, rng)
            xp, _ = design_matrix(cohorts, zp, REF)
            null[p_i] = float(ols_hc1(clr[:, j], xp)["beta"][-1])
        out.append(
            {
                "part": name,
                "beta_clr_per_sd": beta,
                "se_hc1": se,
                "p_hc1": student_p(beta, se, fit["dof"]),
                "p_perm": float((1 + np.sum(np.abs(null) >= abs(beta))) / (1 + N_PERM_ILR)),
                "n": int(fit["n"]),
            }
        )
    # Benjamini-Hochberg on permutation p
    order = np.argsort([r["p_perm"] for r in out])
    m = len(out)
    q = np.empty(m)
    prev = 1.0
    for rank, idx in enumerate(order[::-1], start=1):
        # rank from the largest
        orig_rank = m - rank + 1
        val = out[idx]["p_perm"] * m / orig_rank
        prev = min(prev, val)
        q[idx] = prev
    for i, r in enumerate(out):
        r["q_bh"] = float(min(q[i], 1.0))
    return out


def dm_block(counts, part_names, cohorts, z, rng, simplex: str) -> tuple[list[dict], dict]:
    x, coef_names = design_matrix(cohorts, z, REF)
    fit = fit_dirichlet_multinomial(counts, x, part_names, coef_names, n_starts=3, seed=SEED)
    delta = finite_difference_proportions(fit, x, coef_index=x.shape[1] - 1, delta=1.0)
    # permutation of the CLDN4 column coefficients (reference-coded, K-1 of them)
    null = np.zeros((N_PERM_DM, fit.beta.shape[0]))
    for p_i in range(N_PERM_DM):
        zp = permute_within_cohort(z, cohorts, rng)
        xp, _ = design_matrix(cohorts, zp, REF)
        fit_p = fit_dirichlet_multinomial(
            counts, xp, part_names, coef_names, n_starts=1, seed=SEED + p_i
        )
        null[p_i] = fit_p.beta[:, -1]
    rows = []
    for j, part in enumerate(part_names[:-1]):
        beta = float(fit.beta[j, -1])
        perm_p = float((1 + np.sum(np.abs(null[:, j]) >= abs(beta))) / (1 + N_PERM_DM))
        rows.append(
            {
                "simplex": simplex,
                "part": part,
                "reference": part_names[-1],
                "beta_log_vs_ref_per_sd": beta,
                "p_perm": perm_p,
                "delta_pi_plus_minus_0.5sd": float(delta[j]),
                "phi": float(np.exp(fit.log_phi)),
                "converged": fit.success,
                "n": int(counts.shape[0]),
                "n_perm": N_PERM_DM,
            }
        )
    rows.append(
        {
            "simplex": simplex,
            "part": part_names[-1],
            "reference": part_names[-1],
            "beta_log_vs_ref_per_sd": 0.0,
            "p_perm": float("nan"),
            "delta_pi_plus_minus_0.5sd": float(delta[-1]),
            "phi": float(np.exp(fit.log_phi)),
            "converged": fit.success,
            "n": int(counts.shape[0]),
            "n_perm": N_PERM_DM,
        }
    )
    return rows, {"phi": float(np.exp(fit.log_phi)), "nll": fit.nll, "success": fit.success}


def fraction_meta(tab: pd.DataFrame) -> list[dict]:
    rows = []
    for col, label in [
        ("frac_T", "T"),
        ("frac_NK", "NK"),
        ("frac_myeloid", "myeloid"),
        ("frac_rest", "Rest"),
        ("frac_tnk_split", "T_plus_NK"),
        ("locked_frac_tnk", "locked_T_NK_gate"),
    ]:
        rhos, ns = [], []
        for coh, sub in tab.groupby("dataset"):
            rho, p = spearman_rho(sub.mal_CLDN4_pct.to_numpy(), sub[col].to_numpy())
            rows.append({"fraction": label, "dataset": coh, "n": int(len(sub)), "rho": rho, "p": p})
            rhos.append(rho)
            ns.append(len(sub))
        pooled = dl_fisher_z(np.array(rhos), np.array(ns, dtype=float))
        rows.append(
            {
                "fraction": label,
                "dataset": "DL_meta",
                "n": pooled["N"],
                "rho": pooled["rho"],
                "p": pooled["p"],
                "I2": pooled["I2"],
                "ci_lo": pooled["ci_lo"],
                "ci_hi": pooled["ci_hi"],
                "k": pooled["k"],
            }
        )
        q4 = tab.cldn4_quartile == "Q4"
        q1 = tab.cldn4_quartile == "Q1"
        rb = rank_biserial(tab.loc[q4, col].to_numpy(), tab.loc[q1, col].to_numpy())
        rows.append(
            {
                "fraction": label,
                "dataset": "stacked_Q4_vs_Q1",
                "n": rb["n_q1"] + rb["n_q4"],
                "rho": rb["r"],
                "p": rb["p"],
                "n_q1": rb["n_q1"],
                "n_q4": rb["n_q4"],
            }
        )
    # malignant cellularity, not part of the 4-simplex numerator
    tab = tab.copy()
    tab["frac_malignant"] = tab["n_malignant"] / tab["n_cells"]
    rhos, ns = [], []
    for coh, sub in tab.groupby("dataset"):
        rho, p = spearman_rho(sub.mal_CLDN4_pct.to_numpy(), sub.frac_malignant.to_numpy())
        rows.append({"fraction": "malignant", "dataset": coh, "n": int(len(sub)), "rho": rho, "p": p})
        rhos.append(rho)
        ns.append(len(sub))
    pooled = dl_fisher_z(np.array(rhos), np.array(ns, dtype=float))
    rows.append(
        {
            "fraction": "malignant",
            "dataset": "DL_meta",
            "n": pooled["N"],
            "rho": pooled["rho"],
            "p": pooled["p"],
            "I2": pooled["I2"],
            "ci_lo": pooled["ci_lo"],
            "ci_hi": pooled["ci_hi"],
            "k": pooled["k"],
        }
    )
    return rows


def main() -> None:
    tab = pd.read_csv(TAB / "composition_counts.tsv", sep="\t")
    if len(tab) != 65:
        raise SystemExit(f"expected 65 units, got {len(tab)}")
    cohorts = tab.dataset.to_numpy()
    z = within_cohort_z(tab.mal_CLDN4_pct.to_numpy(), cohorts)
    cldn4 = tab.mal_CLDN4_pct.to_numpy()
    quartile = tab.cldn4_quartile.to_numpy()
    rng = np.random.default_rng(SEED)

    counts4 = tab[COUNT4].to_numpy(dtype=float)
    props4 = proportions_with_pseudocount(counts4, 0.5)
    bal4 = ilr_from_sbp(props4, SBP_TNK_MYE_REST)
    counts3 = tab[COUNT3].to_numpy(dtype=float)
    props3 = proportions_with_pseudocount(counts3, 0.5)
    bal3 = ilr_from_sbp(props3, SBP_WITHIN_IMMUNE)

    unit = tab[["dataset", "unit_id", "mal_CLDN4_pct", "cldn4_quartile"]].copy()
    unit["cldn4_z_within_cohort"] = z
    for j, name in enumerate(BALANCE_NAMES_4):
        unit[f"ilr4_{name}"] = bal4[:, j]
    for j, name in enumerate(BALANCE_NAMES_3):
        unit[f"ilr3_{name}"] = bal3[:, j]
    unit.to_csv(TAB / "ilr_coordinates.tsv", sep="\t", index=False)

    print("ILR 4-part", flush=True)
    rows4, sp4, extra4 = fit_ilr_block(
        "T_NK_myeloid_Rest", bal4, list(BALANCE_NAMES_4), cohorts, z, cldn4, quartile, rng
    )
    print("ILR 3-part", flush=True)
    rows3, sp3, extra3 = fit_ilr_block(
        "within_T_NK_myeloid", bal3, list(BALANCE_NAMES_3), cohorts, z, cldn4, quartile, rng
    )
    pd.DataFrame(rows4 + rows3).to_csv(TAB / "ilr_ols.tsv", sep="\t", index=False)
    pd.DataFrame(sp4 + sp3).to_csv(TAB / "ilr_spearman.tsv", sep="\t", index=False)
    pd.DataFrame(extra4["q4q1"] + extra3["q4q1"]).to_csv(TAB / "ilr_q4q1.tsv", sep="\t", index=False)

    print("CLR", flush=True)
    clr = clr_rows(props4, PARTS4, cohorts, z, rng)
    pd.DataFrame(clr).to_csv(TAB / "clr_ols.tsv", sep="\t", index=False)

    print("fraction meta", flush=True)
    frac_rows = fraction_meta(tab)
    pd.DataFrame(frac_rows).to_csv(TAB / "fraction_spearman.tsv", sep="\t", index=False)

    print("Dirichlet-multinomial 4-part", flush=True)
    dm4, info4 = dm_block(counts4, PARTS4, cohorts, z, rng, "T_NK_myeloid_Rest")
    print("Dirichlet-multinomial within immune", flush=True)
    dm3, info3 = dm_block(counts3, PARTS3, cohorts, z, rng, "within_T_NK_myeloid")
    pd.DataFrame(dm4 + dm3).to_csv(TAB / "dirichlet_multinomial.tsv", sep="\t", index=False)

    # reproduction of the locked T/NK gate
    locked_match = {
        "author_tnk_exact": bool(
            (
                tab.loc[tab.dataset.isin(["GSE131907", "GSE205335"]), "n_tnk_locked_gate"]
                == tab.loc[tab.dataset.isin(["GSE131907", "GSE205335"]), "locked_n_tnk"]
            ).all()
        ),
        "marker_gate_exact": bool(
            (
                tab.loc[tab.dataset.isin(["GSE123902", "GSE189357"]), "n_tnk_locked_gate"]
                == tab.loc[tab.dataset.isin(["GSE123902", "GSE189357"]), "locked_n_tnk"]
            ).all()
        ),
        "malignant_exact_all": bool((tab.n_malignant == tab.locked_n_malignant).all()),
        "cells_exact_all": bool((tab.n_cells == tab.locked_n_cells).all()),
        "extra_marker_cells_in_split": int(
            (
                tab.loc[tab.dataset.isin(["GSE123902", "GSE189357"]), "n_tnk_split"]
                - tab.loc[tab.dataset.isin(["GSE123902", "GSE189357"]), "locked_n_tnk"]
            ).sum()
        ),
    }
    summary = {
        "n_units": 65,
        "n_perm_ilr": N_PERM_ILR,
        "n_perm_dm": N_PERM_DM,
        "reference_cohort": REF,
        "covariate": "within-cohort z of malignant CLDN4 percent positive",
        "locked_match": locked_match,
        "ilr_joint": [extra4["joint"], extra3["joint"]],
        "dirichlet": {"four": info4, "three": info3},
        "quartile_n": {
            "Q1": int((tab.cldn4_quartile == "Q1").sum()),
            "Q4": int((tab.cldn4_quartile == "Q4").sum()),
        },
    }
    (TAB / "composition_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    print("DONE compositional", flush=True)


if __name__ == "__main__":
    main()
