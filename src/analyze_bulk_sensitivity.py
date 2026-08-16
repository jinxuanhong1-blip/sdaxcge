"""Post-hoc sensitivity (not pre-registered): CD47 vs TACSTD2 after
housekeeping-score adjustment, plus a compact CD47-focused summary table.

Triggered because A4 (negative controls) showed that in LUAD, several
housekeepers correlate with TACSTD2 at |rho| similar to CD47. That has to
be shown, not buried.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C
from stats_utils import spearman, partial_spearman, cliffs_delta, zscore_module

DERIVED = C.DATA / "derived"
TAB = C.OUT / "tables"
TAB.mkdir(parents=True, exist_ok=True)

FOCUS = ["CD47", "SIRPA", "THBS1", "NECTIN4", "LGALS3", "TGFB1"]


def main():
    rows = []
    for cohort in ["LUAD", "LUSC"]:
        mat = pd.read_parquet(DERIVED / f"bulk_{cohort}_panel.parquet")
        samp = pd.read_parquet(DERIVED / f"bulk_{cohort}_samples.parquet")
        keep = samp.sample_type == "PrimaryTumor"
        mat, samp = mat.loc[keep], samp.loc[keep]
        tac = mat[C.TARGET]
        hk, _ = zscore_module(mat, C.NEG_CONTROL)
        purity = samp["purity"]

        for gene in FOCUS:
            y = mat[gene]
            naive = spearman(tac, y)
            adj_p = partial_spearman(tac, y, [purity])
            adj_hk = partial_spearman(tac, y, [hk])
            adj_both = partial_spearman(tac, y, [purity, hk])
            rows.append(dict(
                cohort=cohort, gene=gene,
                n_tumor=len(mat), n_with_purity=int(purity.notna().sum()),
                rho_naive=naive["rho"], p_naive=naive["p"],
                rho_adj_purity=adj_p["rho"], p_adj_purity=adj_p["p"],
                rho_adj_housekeeping=adj_hk["rho"], p_adj_housekeeping=adj_hk["p"],
                rho_adj_purity_and_HK=adj_both["rho"], p_adj_both=adj_both["p"],
                lo_adj_both=adj_both["lo"], hi_adj_both=adj_both["hi"],
            ))

        # Confounder correlations themselves.
        for name, y in [("purity", purity), ("housekeeping_score", hk)]:
            r = spearman(tac, y)
            rows.append(dict(
                cohort=cohort, gene=f"TACSTD2_vs_{name}",
                n_tumor=r["n"], n_with_purity=int(purity.notna().sum()),
                rho_naive=r["rho"], p_naive=r["p"],
                rho_adj_purity=np.nan, p_adj_purity=np.nan,
                rho_adj_housekeeping=np.nan, p_adj_housekeeping=np.nan,
                rho_adj_purity_and_HK=np.nan, p_adj_both=np.nan,
                lo_adj_both=r["lo"], hi_adj_both=r["hi"],
            ))

    pd.DataFrame(rows).to_csv(TAB / "14_bulk_cd47_sensitivity.csv", index=False)
    print("wrote", TAB / "14_bulk_cd47_sensitivity.csv")


if __name__ == "__main__":
    main()
