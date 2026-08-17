#!/usr/bin/env python3
"""GSE148071 malignant CLDN4: given continuous ρ + extra Q4 vs Q1 vs T/NK and CXCL13+."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
TABLES = HERE / "tables"
FIGS = HERE / "figures"

# Locked continuous ρ from the TISCH NSCLC pool (eligible n=25).
GIVEN_TISCH = {
    "source": "TISCH NSCLC pool per_dataset_spearman.tsv",
    "contrast": "CLDN4_mean_vs_fracTNK",
    "n": 25,
    "rho": 0.13538461538461538,
    "p": 0.5187661950161317,
}

# Locked continuous ρ from PR #320 / #274 malignant TLS extract.
GIVEN_CXCL13 = {
    "source": "PR #274/#320 tls_malignant × cxcl13pos / mean",
    "contrast": "CLDN4_mal_mean vs frac_CXCL13pos_T",
    "n": 31,
    "rho": -0.4086020621109129,
    "p": 0.02248010028418731,
}


def spearman(x, y) -> dict:
    d = pd.DataFrame({"x": np.asarray(x, dtype=float), "y": np.asarray(y, dtype=float)})
    d = d[np.isfinite(d["x"]) & np.isfinite(d["y"])]
    n = int(len(d))
    if n < 4 or d["x"].nunique() < 2 or d["y"].nunique() < 2:
        return {"n": n, "rho": np.nan, "p": np.nan}
    rho, p = stats.spearmanr(d["x"], d["y"])
    return {"n": n, "rho": float(rho), "p": float(p)}


def q4_vs_q1(predictor, endpoint) -> dict:
    """CLDN4 quartiles on pairwise-complete rows; MWU on the endpoint.

    Rank-biserial r = 2U/(n4 n1) − 1. Positive r = Q4 stochastically higher.
    Quartiles use average ranks then qcut; ties that collapse four bins return
    n_q1/n_q4 as computed. n_compared is n_Q1 + n_Q4, not the full n.
    """
    s = pd.DataFrame(
        {"c": np.asarray(predictor, dtype=float), "i": np.asarray(endpoint, dtype=float)}
    )
    s = s[np.isfinite(s["c"]) & np.isfinite(s["i"])].copy()
    n = int(len(s))
    rec = {
        "n": n,
        "n_q1": np.nan,
        "n_q4": np.nan,
        "n_compared": np.nan,
        "median_q1": np.nan,
        "median_q4": np.nan,
        "delta_median": np.nan,
        "U": np.nan,
        "r_rb": np.nan,
        "p": np.nan,
        "thin": True,
        "q1_cut": np.nan,
        "q4_cut": np.nan,
    }
    if n < 6:
        return rec
    ranks = s["c"].rank(method="average")
    try:
        qs = pd.qcut(ranks, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")
    except ValueError:
        return rec
    if qs.nunique() < 4:
        return rec
    q1 = s.loc[qs == "Q1", "i"]
    q4 = s.loc[qs == "Q4", "i"]
    n1, n4 = int(len(q1)), int(len(q4))
    rec.update(
        {
            "n_q1": n1,
            "n_q4": n4,
            "n_compared": n1 + n4,
            "median_q1": float(q1.median()) if n1 else np.nan,
            "median_q4": float(q4.median()) if n4 else np.nan,
            "thin": n < 8 or n1 < 3 or n4 < 3,
            "q1_cut": float(s.loc[qs == "Q1", "c"].max()) if n1 else np.nan,
            "q4_cut": float(s.loc[qs == "Q4", "c"].min()) if n4 else np.nan,
        }
    )
    if n1 and n4:
        rec["delta_median"] = float(q4.median() - q1.median())
    if n1 < 2 or n4 < 2:
        return rec
    u, p = stats.mannwhitneyu(q4.values, q1.values, alternative="two-sided")
    rec.update(
        {
            "U": float(u),
            "p": float(p),
            "r_rb": (2.0 * float(u)) / (n4 * n1) - 1.0,
        }
    )
    return rec


def fmt_p(p) -> str:
    if p is None or not np.isfinite(p):
        return "—"
    if p < 1e-3:
        return f"{p:.2e}"
    return f"{p:.3g}"


def fmt_r(x, digits=3) -> str:
    if x is None or not np.isfinite(x):
        return "—"
    return f"{x:+.{digits}f}"


def row(analysis, endpoint, predictor, source, unit, spear, q, note="") -> dict:
    return {
        "analysis": analysis,
        "endpoint": endpoint,
        "predictor": predictor,
        "source": source,
        "unit": unit,
        "n": spear["n"],
        "rho": spear["rho"],
        "p_spearman": spear["p"],
        "n_q1": q["n_q1"],
        "n_q4": q["n_q4"],
        "n_compared": q["n_compared"],
        "median_q1": q["median_q1"],
        "median_q4": q["median_q4"],
        "delta_median": q["delta_median"],
        "U": q["U"],
        "r_rb": q["r_rb"],
        "p_q4q1": q["p"],
        "thin_q4q1": q["thin"],
        "q1_cut": q["q1_cut"],
        "q4_cut": q["q4_cut"],
        "note": note,
    }


def boxplot(q1, q4, ylabel, title, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(4.2, 3.8))
    bp = ax.boxplot(
        [q1, q4],
        tick_labels=[f"Q1\nn={len(q1)}", f"Q4\nn={len(q4)}"],
        patch_artist=True,
        widths=0.55,
    )
    for patch, color in zip(bp["boxes"], ["#6a8aaa", "#c44e52"]):
        patch.set_facecolor(color)
        patch.set_alpha(0.55)
    rng = np.random.default_rng(0)
    for i, vals in enumerate((q1, q4), start=1):
        ax.scatter(
            i + rng.uniform(-0.08, 0.08, size=len(vals)),
            vals,
            c="black",
            s=16,
            zorder=3,
        )
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def q_vectors(predictor, endpoint):
    s = pd.DataFrame(
        {"c": np.asarray(predictor, dtype=float), "i": np.asarray(endpoint, dtype=float)}
    )
    s = s[np.isfinite(s["c"]) & np.isfinite(s["i"])].copy()
    ranks = s["c"].rank(method="average")
    qs = pd.qcut(ranks, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")
    return s.loc[qs == "Q1", "i"].to_numpy(), s.loc[qs == "Q4", "i"].to_numpy()


def load() -> tuple[pd.DataFrame, pd.DataFrame]:
    tisch = pd.read_csv(DATA / "tisch_GSE148071_units.tsv", sep="\t")
    tls = pd.read_csv(DATA / "tls_GSE148071_patients.tsv", sep="\t")
    tisch["patient"] = tisch["patient"].astype(str)
    tls["patient"] = tls["patient"].astype(str)
    return tisch, tls


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    tisch, tls = load()

    n_geo = 42
    assert len(tisch) == n_geo, f"TISCH rows {len(tisch)} != 42"
    assert len(tls) == n_geo, f"TLS rows {len(tls)} != 42"

    tisch_elig = tisch[tisch["eligible"]].copy()
    tisch_malig = tisch_elig[tisch_elig["epi_definition"] == "Malignant"].copy()
    tisch_epi_like = tisch_elig[tisch_elig["epi_definition"] != "Malignant"].copy()

    tls_malig = tls[tls["compartment"] == "malignant"].copy()
    tls_cx = tls_malig[np.isfinite(tls_malig["frac_CXCL13pos_T"])].copy()
    tls_cx_mean = tls_malig[np.isfinite(tls_malig["cxcl13_mean_log1p_cp10k"])].copy()

    rows = []

    # --- TISCH given eligible set (locked continuous ρ) ---
    s_tnk_given = spearman(tisch_elig["CLDN4_epi_mean"], tisch_elig["frac_tnk"])
    q_tnk_given = q4_vs_q1(tisch_elig["CLDN4_epi_mean"], tisch_elig["frac_tnk"])
    if abs(s_tnk_given["rho"] - GIVEN_TISCH["rho"]) > 1e-9 or s_tnk_given["n"] != GIVEN_TISCH["n"]:
        raise SystemExit(
            f"TISCH given ρ mismatch: recomputed n={s_tnk_given['n']} "
            f"ρ={s_tnk_given['rho']} vs locked n={GIVEN_TISCH['n']} ρ={GIVEN_TISCH['rho']}"
        )
    rows.append(
        row(
            "tisch_eligible_given",
            "frac_TNK",
            "TISCH CLDN4 mean (Malignant if ≥20 else epithelial-like)",
            "TISCH NSCLC pool",
            "patient",
            s_tnk_given,
            q_tnk_given,
            "Locked continuous ρ taken as given. Q4 vs Q1 is extra. "
            f"Includes {len(tisch_epi_like)} epithelial-like fallbacks "
            f"({', '.join(sorted(tisch_epi_like['patient']))}).",
        )
    )

    s_tnk_pct = spearman(tisch_elig["CLDN4_epi_pctpos"], tisch_elig["frac_tnk"])
    q_tnk_pct = q4_vs_q1(tisch_elig["CLDN4_epi_pctpos"], tisch_elig["frac_tnk"])
    rows.append(
        row(
            "tisch_eligible_pct",
            "frac_TNK",
            "TISCH CLDN4 %pos",
            "TISCH NSCLC pool",
            "patient",
            s_tnk_pct,
            q_tnk_pct,
            "Secondary %pos on the same locked eligible n=25.",
        )
    )

    s_tnk_mal = spearman(tisch_malig["CLDN4_epi_mean"], tisch_malig["frac_tnk"])
    q_tnk_mal = q4_vs_q1(tisch_malig["CLDN4_epi_mean"], tisch_malig["frac_tnk"])
    rows.append(
        row(
            "tisch_malignant_only",
            "frac_TNK",
            "TISCH CLDN4 mean (Malignant only)",
            "TISCH NSCLC pool",
            "patient",
            s_tnk_mal,
            q_tnk_mal,
            "Drops the 3 epithelial-like eligible units. Honest malignant recut.",
        )
    )

    # --- TLS CXCL13+ given ---
    s_cx = spearman(tls_cx["cldn4_mal_mean_log1p_cp10k"], tls_cx["frac_CXCL13pos_T"])
    q_cx = q4_vs_q1(tls_cx["cldn4_mal_mean_log1p_cp10k"], tls_cx["frac_CXCL13pos_T"])
    if abs(s_cx["rho"] - GIVEN_CXCL13["rho"]) > 1e-9 or s_cx["n"] != GIVEN_CXCL13["n"]:
        raise SystemExit(
            f"CXCL13 given ρ mismatch: recomputed n={s_cx['n']} ρ={s_cx['rho']} "
            f"vs locked n={GIVEN_CXCL13['n']} ρ={GIVEN_CXCL13['rho']}"
        )
    rows.append(
        row(
            "tls_malignant_cxcl13pos_given",
            "frac_CXCL13pos_T",
            "TLS malignant CLDN4 mean log1p(CP10k)",
            "PR #274 TLS extract",
            "patient",
            s_cx,
            q_cx,
            "Locked continuous ρ taken as given. Q4 vs Q1 is extra. "
            "Requires malignant compartment and at least one T cell with CXCL13 call.",
        )
    )

    s_cx_pct = spearman(tls_cx["cldn4_mal_pct_pos"], tls_cx["frac_CXCL13pos_T"])
    q_cx_pct = q4_vs_q1(tls_cx["cldn4_mal_pct_pos"], tls_cx["frac_CXCL13pos_T"])
    rows.append(
        row(
            "tls_malignant_cxcl13pos_pct",
            "frac_CXCL13pos_T",
            "TLS malignant CLDN4 %pos",
            "PR #274 TLS extract",
            "patient",
            s_cx_pct,
            q_cx_pct,
            "Secondary %pos on the same CXCL13+ complete n=31.",
        )
    )

    s_cxm = spearman(tls_cx_mean["cldn4_mal_mean_log1p_cp10k"], tls_cx_mean["cxcl13_mean_log1p_cp10k"])
    q_cxm = q4_vs_q1(tls_cx_mean["cldn4_mal_mean_log1p_cp10k"], tls_cx_mean["cxcl13_mean_log1p_cp10k"])
    rows.append(
        row(
            "tls_malignant_cxcl13_mean",
            "CXCL13 mean log1p(CP10k)",
            "TLS malignant CLDN4 mean log1p(CP10k)",
            "PR #274 TLS extract",
            "patient",
            s_cxm,
            q_cxm,
            "Whole-sample CXCL13 mean (not restricted to T). n=35 malignant.",
        )
    )

    # --- Same-patient intersection: TLS malignant CLDN4 + CXCL13+ + TISCH T/NK ---
    both = tls_cx.merge(
        tisch_elig[["patient", "frac_tnk", "n_tnk", "n_malignant", "CLDN4_epi_mean", "epi_definition", "eligible"]],
        on="patient",
        how="inner",
        suffixes=("", "_tisch"),
    )
    both_mal = both[both["epi_definition"] == "Malignant"].copy()

    s_both_tnk = spearman(both["cldn4_mal_mean_log1p_cp10k"], both["frac_tnk"])
    q_both_tnk = q4_vs_q1(both["cldn4_mal_mean_log1p_cp10k"], both["frac_tnk"])
    rows.append(
        row(
            "same_patient_tnk",
            "frac_TNK",
            "TLS malignant CLDN4 mean (shared cut)",
            "TLS ∩ TISCH eligible",
            "patient",
            s_both_tnk,
            q_both_tnk,
            "Same patients as CXCL13+ row below. One TLS CLDN4 quartile cut. "
            "T/NK from TISCH. Continuous ρ here is not the locked TISCH-internal given.",
        )
    )

    s_both_cx = spearman(both["cldn4_mal_mean_log1p_cp10k"], both["frac_CXCL13pos_T"])
    q_both_cx = q4_vs_q1(both["cldn4_mal_mean_log1p_cp10k"], both["frac_CXCL13pos_T"])
    rows.append(
        row(
            "same_patient_cxcl13pos",
            "frac_CXCL13pos_T",
            "TLS malignant CLDN4 mean (shared cut)",
            "TLS ∩ TISCH eligible",
            "patient",
            s_both_cx,
            q_both_cx,
            "Same patients and same CLDN4 Q4/Q1 as the T/NK intersection row.",
        )
    )

    # Shared-cut Q4 vs Q1: assign quartiles once, test both endpoints.
    both = both.copy()
    ranks = both["cldn4_mal_mean_log1p_cp10k"].rank(method="average")
    both["quartile"] = pd.qcut(ranks, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")
    shared_patients = both[["patient", "quartile", "cldn4_mal_mean_log1p_cp10k", "frac_tnk", "frac_CXCL13pos_T", "n_malignant", "n_tnk"]].copy()

    out = pd.DataFrame(rows)
    out.to_csv(TABLES / "q4q1_table.tsv", sep="\t", index=False)
    shared_patients.to_csv(TABLES / "same_patient_quartiles.tsv", sep="\t", index=False)

    inventory = pd.DataFrame(
        [
            {"item": "GEO patients (Wu 2021)", "n": n_geo, "note": "one sample each; not the test n"},
            {"item": "TISCH units in extract", "n": len(tisch), "note": "all 42 samples present"},
            {
                "item": "TISCH eligible (≥20 scored epi + ≥20 T/NK)",
                "n": len(tisch_elig),
                "note": "locked given-ρ set",
            },
            {
                "item": "TISCH eligible but epithelial-like fallback",
                "n": len(tisch_epi_like),
                "note": ",".join(sorted(tisch_epi_like["patient"])),
            },
            {
                "item": "TISCH eligible and Malignant",
                "n": len(tisch_malig),
                "note": "honest malignant T/NK recut",
            },
            {
                "item": "TISCH ineligible",
                "n": int((~tisch["eligible"]).sum()),
                "note": "usually <20 T/NK",
            },
            {"item": "TLS eligible (author extract)", "n": int(tls["eligible"].sum()), "note": "41/42"},
            {
                "item": "TLS malignant compartment",
                "n": len(tls_malig),
                "note": "epithelial / insufficient dropped",
            },
            {
                "item": "TLS malignant + finite CXCL13+ among T",
                "n": len(tls_cx),
                "note": "locked given-ρ set; 4 malignant lack a T-cell CXCL13 call",
            },
            {
                "item": "TLS malignant + CXCL13 mean",
                "n": len(tls_cx_mean),
                "note": "whole-sample CXCL13; n=35",
            },
            {
                "item": "Same-patient TLS CXCL13+ ∩ TISCH eligible",
                "n": len(both),
                "note": "shared CLDN4 Q4/Q1",
            },
            {
                "item": "Same-patient ∩ TISCH Malignant only",
                "n": len(both_mal),
                "note": "epi-like TISCH units are not in the TLS malignant+CXCL13+ set",
            },
        ]
    )
    inventory.to_csv(TABLES / "honest_n.tsv", sep="\t", index=False)

    # Figures: given-set Q4 boxes + same-patient shared cut.
    q1, q4 = q_vectors(tisch_elig["CLDN4_epi_mean"], tisch_elig["frac_tnk"])
    boxplot(q1, q4, "T/NK fraction", "GSE148071 TISCH eligible\nCLDN4 Q4 vs Q1 vs T/NK", FIGS / "q4q1_tisch_tnk.png")
    q1, q4 = q_vectors(tls_cx["cldn4_mal_mean_log1p_cp10k"], tls_cx["frac_CXCL13pos_T"])
    boxplot(q1, q4, "CXCL13+ among T", "GSE148071 TLS malignant\nCLDN4 Q4 vs Q1 vs CXCL13+", FIGS / "q4q1_tls_cxcl13.png")
    q1 = both.loc[both["quartile"] == "Q1", "frac_tnk"].to_numpy()
    q4 = both.loc[both["quartile"] == "Q4", "frac_tnk"].to_numpy()
    boxplot(q1, q4, "T/NK fraction", "GSE148071 same-patient\nshared CLDN4 Q4 vs Q1 vs T/NK", FIGS / "q4q1_same_tnk.png")
    q1 = both.loc[both["quartile"] == "Q1", "frac_CXCL13pos_T"].to_numpy()
    q4 = both.loc[both["quartile"] == "Q4", "frac_CXCL13pos_T"].to_numpy()
    boxplot(q1, q4, "CXCL13+ among T", "GSE148071 same-patient\nshared CLDN4 Q4 vs Q1 vs CXCL13+", FIGS / "q4q1_same_cxcl13.png")

    write_finding(out, inventory, both)
    summary = {
        "cohort": "GSE148071",
        "paper": "Wu et al. Nat Commun 2021",
        "given_tisch_tnk": GIVEN_TISCH,
        "given_cxcl13": GIVEN_CXCL13,
        "recomputed_match": True,
        "n_same_patient": int(len(both)),
        "q4q1_is_extra": True,
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(out.to_string(index=False))
    print("\nHonest n:")
    print(inventory.to_string(index=False))


def write_finding(out: pd.DataFrame, inventory: pd.DataFrame, both: pd.DataFrame) -> None:
    def get(name: str) -> pd.Series:
        hit = out[out["analysis"] == name]
        if hit.empty:
            raise KeyError(name)
        return hit.iloc[0]

    tnk = get("tisch_eligible_given")
    tnk_mal = get("tisch_malignant_only")
    cx = get("tls_malignant_cxcl13pos_given")
    same_t = get("same_patient_tnk")
    same_c = get("same_patient_cxcl13pos")
    cxm = get("tls_malignant_cxcl13_mean")
    tnk_pct = get("tisch_eligible_pct")
    cx_pct = get("tls_malignant_cxcl13pos_pct")

    q_counts = both["quartile"].value_counts().to_dict()
    n_q1 = int(q_counts.get("Q1", 0))
    n_q4 = int(q_counts.get("Q4", 0))

    def qline(r) -> str:
        thin = " thin" if bool(r.thin_q4q1) else ""
        return (
            f"{int(r.n_q4)} vs {int(r.n_q1)}{thin} | "
            f"{r.median_q4:.4g} | {r.median_q1:.4g} | {r.delta_median:+.4g} | "
            f"{fmt_r(r.r_rb)} | {fmt_p(r.p_q4q1)}"
        )

    lines = [
        "# Finding — GSE148071 malignant CLDN4 Q4 vs Q1 vs same-patient T/NK and CXCL13+",
        "",
        "ADDITIVE. **CLDN4 only.** No TACSTD2 gate. Patient is the unit (42 advanced",
        "NSCLC tumors, Wu *et al.* *Nat Commun* 2021; one sample each). Continuous",
        "Spearman ρ is **taken as given** from the existing TISCH T/NK extract and",
        "the PR #274 / #320 malignant TLS extract. **Q4 vs Q1 is extra**",
        "(Mann–Whitney, rank-biserial *r*). p-values are descriptive.",
        "",
        "Existing processed scores only. The 10x matrices were not re-downloaded.",
        "No ICI response / RECIST / MPR labels are in these extracts — do not read",
        "T/NK or CXCL13+ as immunotherapy outcome.",
        "",
        "## Verdict",
        "",
        "| Contrast | Given continuous ρ | Extra Q4 vs Q1 r | Honest n |",
        "|---|---|---|---|",
        (
            f"| TISCH T/NK (locked eligible) | "
            f"n={int(tnk.n)} ρ={fmt_r(tnk.rho)} p={fmt_p(tnk.p_spearman)} | "
            f"n_Q1={int(tnk.n_q1)} n_Q4={int(tnk.n_q4)} r={fmt_r(tnk.r_rb)} p={fmt_p(tnk.p_q4q1)} | "
            f"25 eligible / 42; 3 epithelial-like inside the 25 |"
        ),
        (
            f"| TISCH T/NK (malignant only) | "
            f"n={int(tnk_mal.n)} ρ={fmt_r(tnk_mal.rho)} p={fmt_p(tnk_mal.p_spearman)} | "
            f"n_Q1={int(tnk_mal.n_q1)} n_Q4={int(tnk_mal.n_q4)} r={fmt_r(tnk_mal.r_rb)} p={fmt_p(tnk_mal.p_q4q1)} | "
            f"22; drops P5/P35/P39 epithelial-like |"
        ),
        (
            f"| TLS CXCL13+ among T (locked malignant) | "
            f"n={int(cx.n)} ρ={fmt_r(cx.rho)} p={fmt_p(cx.p_spearman)} | "
            f"n_Q1={int(cx.n_q1)} n_Q4={int(cx.n_q4)} r={fmt_r(cx.r_rb)} p={fmt_p(cx.p_q4q1)} | "
            f"31 malignant with a T-cell CXCL13 call |"
        ),
        (
            f"| Same-patient T/NK (shared TLS CLDN4 cut) | "
            f"n={int(same_t.n)} ρ={fmt_r(same_t.rho)} p={fmt_p(same_t.p_spearman)} | "
            f"n_Q1={int(same_t.n_q1)} n_Q4={int(same_t.n_q4)} r={fmt_r(same_t.r_rb)} p={fmt_p(same_t.p_q4q1)} | "
            f"{int(same_t.n)} patients in both extracts |"
        ),
        (
            f"| Same-patient CXCL13+ (shared TLS CLDN4 cut) | "
            f"n={int(same_c.n)} ρ={fmt_r(same_c.rho)} p={fmt_p(same_c.p_spearman)} | "
            f"n_Q1={int(same_c.n_q1)} n_Q4={int(same_c.n_q4)} r={fmt_r(same_c.r_rb)} p={fmt_p(same_c.p_q4q1)} | "
            f"same {int(same_c.n)} patients as the T/NK intersection |"
        ),
        "",
        "**What holds:** On the locked TLS malignant extract, higher malignant CLDN4",
        f"tracks **lower CXCL13+ among T** (given ρ={fmt_r(cx.rho)}, p={fmt_p(cx.p_spearman)},",
        f"n={int(cx.n)}). The extra Q4 vs Q1 cut stays negative",
        f"(r={fmt_r(cx.r_rb)}, p={fmt_p(cx.p_q4q1)}, {int(cx.n_q4)} vs {int(cx.n_q1)}).",
        "",
        "**What does not hold:** Same-patient T/NK is **not** a CLDN4-low-TME story.",
        f"The locked TISCH continuous ρ is weakly **positive** (n=25, ρ={fmt_r(tnk.rho)},",
        f"p={fmt_p(tnk.p_spearman)}). Extra Q4 vs Q1 on that set is also not a negative",
        f"hit (r={fmt_r(tnk.r_rb)}, p={fmt_p(tnk.p_q4q1)}). Joining TLS malignant CLDN4",
        f"to TISCH T/NK on the {int(same_t.n)} overlapping patients does not flip T/NK",
        "to a significant negative quartile effect. Do not write n=42. Do not pool",
        "T/NK and CXCL13+ into one immune score.",
        "",
        "## Honest n",
        "",
        "| item | n | note |",
        "|---|---:|---|",
    ]
    for r in inventory.itertuples(index=False):
        lines.append(f"| {r.item} | {int(r.n)} | {r.note} |")
    lines += [
        "",
        "Quartiles are assigned **inside** each analysis set. n_compared = n_Q1 + n_Q4",
        "(the tails), not the GEO n=42 and not the continuous n. Same-patient shared",
        f"cut on n={len(both)}: Q1={n_q1}, Q2={int(q_counts.get('Q2', 0))}, "
        f"Q3={int(q_counts.get('Q3', 0))}, Q4={n_q4}.",
        "",
        "## Methods (this slice)",
        "",
        "- Predictor: malignant CLDN4 only. TACSTD2 is not a gate.",
        "- T/NK (TISCH): `n_TNK / n_cells`. T/NK labels = CD8T, CD8Tex, CD4Tconv,",
        "  Treg, Tprolif, TMKI67, NK, Tcell, NKT, ILC. Eligible = ≥20 scored",
        "  epithelial/malignant cells and ≥20 T/NK (locked given set).",
        "- CXCL13+ (TLS PR #274): `frac_CXCL13pos_T` among T in the same patient;",
        "  CLDN4 = mean log1p(CP10k) in the **malignant** compartment. Epithelial",
        "  and insufficient rows are dropped.",
        "- Given continuous: Spearman on the locked complete sets above. Recomputed",
        "  here only to confirm the file still matches; the ρ is not a new audit.",
        "- Extra cut: `pd.qcut(rank(method='average'), 4)` on pairwise-complete",
        "  CLDN4. Two-sided Mann–Whitney U, Q4 vs Q1. Rank-biserial",
        "  `r = 2U/(n4 n1) − 1` (positive = Q4 higher).",
        "- Same-patient extra: inner join on `P1…P42`. One TLS malignant CLDN4",
        "  quartile cut; both endpoints tested on those tails.",
        "- Thin flag: n<8 or either tail <3. p-values are descriptive.",
        "",
        "## Extra Q4 vs Q1 (primary table)",
        "",
        "| Contrast | n Q4 / Q1 | median Q4 | median Q1 | Δ (Q4−Q1) | r | p |",
        "|---|---|---:|---:|---:|---:|---:|",
        f"| TISCH T/NK, locked eligible | {qline(tnk)} |",
        f"| TISCH T/NK, malignant only | {qline(tnk_mal)} |",
        f"| TLS CXCL13+ among T, locked malignant | {qline(cx)} |",
        f"| Same-patient T/NK, shared TLS CLDN4 cut | {qline(same_t)} |",
        f"| Same-patient CXCL13+, shared TLS CLDN4 cut | {qline(same_c)} |",
        f"| TLS CXCL13 mean (secondary) | {qline(cxm)} |",
        "",
        "## Given continuous Spearman (supporting; not re-audited)",
        "",
        "| Contrast | n | ρ | p | source |",
        "|---|---:|---:|---:|---|",
        (
            f"| TISCH CLDN4 mean vs T/NK | {int(tnk.n)} | {fmt_r(tnk.rho)} | "
            f"{fmt_p(tnk.p_spearman)} | TISCH pool, eligible units |"
        ),
        (
            f"| TISCH CLDN4 %pos vs T/NK | {int(tnk_pct.n)} | {fmt_r(tnk_pct.rho)} | "
            f"{fmt_p(tnk_pct.p_spearman)} | same n=25 |"
        ),
        (
            f"| TISCH CLDN4 mean vs T/NK, malignant only | {int(tnk_mal.n)} | "
            f"{fmt_r(tnk_mal.rho)} | {fmt_p(tnk_mal.p_spearman)} | drop 3 epi-like |"
        ),
        (
            f"| TLS CLDN4 mean vs CXCL13+ among T | {int(cx.n)} | {fmt_r(cx.rho)} | "
            f"{fmt_p(cx.p_spearman)} | PR #274/#320 malignant |"
        ),
        (
            f"| TLS CLDN4 %pos vs CXCL13+ among T | {int(cx_pct.n)} | {fmt_r(cx_pct.rho)} | "
            f"{fmt_p(cx_pct.p_spearman)} | same n=31 |"
        ),
        (
            f"| TLS CLDN4 mean vs CXCL13 mean | {int(cxm.n)} | {fmt_r(cxm.rho)} | "
            f"{fmt_p(cxm.p_spearman)} | malignant n=35 |"
        ),
        (
            f"| Same-patient TLS CLDN4 vs TISCH T/NK | {int(same_t.n)} | "
            f"{fmt_r(same_t.rho)} | {fmt_p(same_t.p_spearman)} | join; not the locked TISCH ρ |"
        ),
        (
            f"| Same-patient TLS CLDN4 vs CXCL13+ | {int(same_c.n)} | "
            f"{fmt_r(same_c.rho)} | {fmt_p(same_c.p_spearman)} | join subset of the n=31 |"
        ),
        "",
        "Q4 vs Q1 throws away Q2+Q3. That is why a continuous CXCL13+ hit can weaken",
        "on the intersection tails. Do not upgrade the n=31 Spearman to an n=42 claim,",
        "and do not upgrade the T/NK null to a negative Q4 story.",
        "",
        "## What was not done",
        "",
        "- No dual-high TACSTD2×CLDN4 score.",
        "- No new 10x download or re-annotation of GSE148071.",
        "- No ICI response test (none in these extracts).",
        "- TISCH and TLS annotations are different; cell counts do not match 1:1.",
        "  Same-patient uses the patient ID join, not a cell-level merge.",
        "- T/NK was not invented on the TLS table (no T/NK column there).",
        "",
        "Tables: `tables/q4q1_table.tsv`, `tables/honest_n.tsv`,",
        "`tables/same_patient_quartiles.tsv`.",
        "",
        "Figures: `figures/q4q1_tisch_tnk.png`, `figures/q4q1_tls_cxcl13.png`,",
        "`figures/q4q1_same_tnk.png`, `figures/q4q1_same_cxcl13.png`.",
        "",
        "Reproduce: `python3 methods/gse148071_cldn4_q4/analyze.py`",
        "",
    ]
    (HERE / "FINDING.md").write_text("\n".join(lines))


if __name__ == "__main__":
    main()
