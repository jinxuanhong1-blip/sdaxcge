"""Verification / robustness layer for the TACSTD2-CLDN4 findings.

Independent re-derivation and stress-tests of the two headline claims:

  Claim 1 (driver): TACSTD2 and CLDN4 are higher in EGFR-mutant vs EGFR/ALK-WT LUAD.
  Claim 2 (immune): TACSTD2 and CLDN4 correlate negatively with cytotoxic/CD8
                    immune signatures (immune-excluded phenotype).

Checks:
  A. Bootstrap 95% CIs for the key Spearman correlations (LUAD & LUSC).
  B. Permutation-test p-values for those correlations.
  C. Cross-cohort sign concordance (LUAD vs LUSC) across the full immune panel.
  D. Bootstrap 95% CI of the EGFR-mutant vs WT median difference (Claim 1).
  E. Sensitivity of Claim 2 to alternative immune scoring (single-gene CD8A,
     GZMB; and an ESTIMATE-style broad immune panel).
  F. Positive control recap: CD8 signature separates ICI responders (GSE126044).

Writes tables/verification_report.csv and prints a PASS/FAIL summary.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

import common as C

RNG = np.random.default_rng(C.RANDOM_SEED)
N_BOOT = 5000
N_PERM = 10000
KEY_IMMUNE = "sig_CD8_effector"


def load_sig(cohort: str) -> pd.DataFrame:
    return pd.read_csv(f"{C.TABLES}/tcga_{cohort.lower()}_signature_scores.csv")


def boot_spearman_ci(x, y, n=N_BOOT):
    x = np.asarray(x, float); y = np.asarray(y, float)
    ok = ~(np.isnan(x) | np.isnan(y))
    x, y = x[ok], y[ok]
    obs = stats.spearmanr(x, y).statistic
    idx = np.arange(len(x))
    boots = np.empty(n)
    for i in range(n):
        b = RNG.choice(idx, size=len(idx), replace=True)
        boots[i] = stats.spearmanr(x[b], y[b]).statistic
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return obs, lo, hi


def perm_spearman_p(x, y, n=N_PERM):
    x = np.asarray(x, float); y = np.asarray(y, float)
    ok = ~(np.isnan(x) | np.isnan(y))
    x, y = x[ok], y[ok]
    obs = stats.spearmanr(x, y).statistic
    count = 0
    for _ in range(n):
        count += abs(stats.spearmanr(x, RNG.permutation(y)).statistic) >= abs(obs)
    return obs, (count + 1) / (n + 1)


def boot_median_diff_ci(a, b, n=N_BOOT):
    a = np.asarray(a, float); b = np.asarray(b, float)
    obs = np.median(a) - np.median(b)
    boots = np.empty(n)
    for i in range(n):
        ba = RNG.choice(a, size=len(a), replace=True)
        bb = RNG.choice(b, size=len(b), replace=True)
        boots[i] = np.median(ba) - np.median(bb)
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return obs, lo, hi


def main() -> None:
    rows = []

    # ---- A + B : bootstrap CI & permutation p for target vs CD8 effector ----
    for cohort in ("LUAD", "LUSC"):
        df = load_sig(cohort)
        for tgt in C.TARGETS:
            obs, lo, hi = boot_spearman_ci(df[tgt], df[KEY_IMMUNE])
            _, pperm = perm_spearman_p(df[tgt], df[KEY_IMMUNE])
            rows.append({
                "check": "A/B_corr_target_vs_CD8effector", "cohort": cohort,
                "target": tgt, "metric": "spearman_rho", "value": obs,
                "ci_lo": lo, "ci_hi": hi, "perm_p": pperm,
                "pass": bool(hi < 0),  # CI strictly below 0 => robust negative
            })

    # ---- C : cross-cohort sign concordance across full immune panel ----
    ov_luad = pd.read_csv(f"{C.TABLES}/tcga_luad_target_immune_correlations.csv")
    ov_lusc = pd.read_csv(f"{C.TABLES}/tcga_lusc_target_immune_correlations.csv")
    for tgt in C.TARGETS:
        a = ov_luad[(ov_luad.group == "ALL") & (ov_luad.target == tgt)][["feature", "spearman_rho"]]
        b = ov_lusc[(ov_lusc.group == "ALL") & (ov_lusc.target == tgt)][["feature", "spearman_rho"]]
        m = a.merge(b, on="feature", suffixes=("_luad", "_lusc"))
        concord = (np.sign(m["spearman_rho_luad"]) == np.sign(m["spearman_rho_lusc"])).sum()
        total = len(m)
        binom_p = stats.binomtest(int(concord), total, 0.5, alternative="greater").pvalue
        rows.append({
            "check": "C_cross_cohort_sign_concordance", "cohort": "LUAD_vs_LUSC",
            "target": tgt, "metric": "concordant/total", "value": concord / total,
            "ci_lo": concord, "ci_hi": total, "perm_p": binom_p,
            "pass": bool(binom_p < 0.05),
        })

    # ---- D : bootstrap CI of EGFR-mut vs WT median difference ----
    for cohort in ("LUAD",):  # EGFR-mut only well powered in LUAD
        df = load_sig(cohort)
        for tgt in C.TARGETS:
            a = df.loc[df.driver_group == "EGFR_mut", tgt].dropna().values
            b = df.loc[df.driver_group == "EGFR_ALK_wt", tgt].dropna().values
            obs, lo, hi = boot_median_diff_ci(a, b)
            u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
            rows.append({
                "check": "D_EGFRmut_vs_WT_median_diff", "cohort": cohort,
                "target": tgt, "metric": "median_diff_log2TPM", "value": obs,
                "ci_lo": lo, "ci_hi": hi, "perm_p": p,
                "pass": bool(lo > 0),  # EGFR-mut strictly higher
            })

    # ---- E : sensitivity to alternative immune scoring (LUAD & LUSC) ----
    for cohort in ("LUAD", "LUSC"):
        sig = load_sig(cohort)
        samp = pd.read_csv(f"{C.PROCESSED}/tcga_{cohort.lower()}_samples.csv")
        m = sig.merge(samp[["sample", "CD8A", "GZMB"]], on="sample", how="left")
        for tgt in C.TARGETS:
            for alt in ("CD8A", "GZMB", "sig_Broad_immune"):
                rho, p = stats.spearmanr(m[tgt], m[alt], nan_policy="omit")
                rows.append({
                    "check": "E_alt_immune_scoring", "cohort": cohort,
                    "target": tgt, "metric": f"spearman_vs_{alt}", "value": rho,
                    "ci_lo": np.nan, "ci_hi": np.nan, "perm_p": p,
                    "pass": bool(rho < 0 and p < 0.05),
                })

    # ---- F : positive control (ICI immune signature predicts response) ----
    g126 = pd.read_csv(f"{C.TABLES}/ici_gse126044_target_response.csv")
    cd8 = g126[g126.feature == "sig_CD8_effector"].iloc[0]
    rows.append({
        "check": "F_poscontrol_CD8sig_predicts_response", "cohort": "GSE126044",
        "target": "sig_CD8_effector", "metric": "AUC_resp_high", "value": cd8["auc_resp_high"],
        "ci_lo": np.nan, "ci_hi": np.nan, "perm_p": cd8["p"],
        "pass": bool(cd8["auc_resp_high"] > 0.7 and cd8["p"] < 0.05),
    })

    rep = pd.DataFrame(rows)
    rep.to_csv(f"{C.TABLES}/verification_report.csv", index=False)

    pd.set_option("display.width", 170)
    print("=== VERIFICATION REPORT ===")
    print(rep.to_string(index=False))
    n_pass = int(rep["pass"].sum())
    print(f"\nPASSED {n_pass}/{len(rep)} checks")

    # concise headline verdict
    print("\n--- Headline claim verdicts ---")
    cA = rep[rep.check.str.startswith("A/B")]
    print("Claim 2 (immune-excluded, CD8 effector), all robustly negative CIs:",
          bool(cA["pass"].all()))
    cD = rep[rep.check.str.startswith("D_")]
    print("Claim 1 (EGFR-mut higher target in LUAD), CI>0 both targets:",
          bool(cD["pass"].all()))


if __name__ == "__main__":
    main()
