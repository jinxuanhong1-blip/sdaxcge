#!/usr/bin/env python3
"""Per-series table + detected-cut stats from already-scored leftover files."""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from common import OUT, dump_json, high_low_cut, spearman
from importlib import import_module

spatial_mod = import_module("03_score_spatial") if False else None


def _neighbor_mean(x, y, values, k=6):
    from importlib.machinery import SourceFileLoader
    from pathlib import Path

    mod = SourceFileLoader("score_spatial", str(Path(__file__).with_name("03_score_spatial.py"))).load_module()
    return mod.neighbor_mean(x, y, values, k=k)


def scrna_cuts():
    cells = pd.read_csv(OUT / "scrna_cell_scores.tsv.gz", sep="\t")
    rows = []
    for (gse, sample), sub in cells.groupby(["gse", "sample"]):
        epi = sub[sub["label"] == "epithelial"]
        tnk = sub[sub["label"] == "tnk"]
        mask, meta = high_low_cut(epi["Cldn4"], epi["Tscore"], how="detected")
        rho = spearman(epi["Cldn4"], epi["Tscore"]) if len(epi) >= 8 else {"n": len(epi), "rho": np.nan, "p": np.nan}
        rows.append(
            {
                "gse": gse,
                "sample": sample,
                "n_epi": int(len(epi)),
                "n_tnk": int(len(tnk)),
                "frac_tnk": float(len(tnk) / len(sub)) if len(sub) else np.nan,
                "epi_Cldn4_mean": float(epi["Cldn4"].mean()) if len(epi) else np.nan,
                "tnk_Cldn4_mean": float(tnk["Cldn4"].mean()) if len(tnk) else np.nan,
                "epi_Cldn4_pct": float((epi["Cldn4"] > 0).mean() * 100) if len(epi) else np.nan,
                "epi_rho_Cldn4_T": rho["rho"],
                "epi_rho_p": rho["p"],
                **{f"cut_{k}": v for k, v in meta.items()},
            }
        )
    return pd.DataFrame(rows)


def spatial_cuts():
    spots = pd.read_csv(OUT / "spatial_spot_scores.tsv.gz", sep="\t")
    rows = []
    for sample, sub in spots.groupby("sample"):
        mask, meta = high_low_cut(sub["Cldn4"], sub["Tscore"], how="detected")
        rho = spearman(sub["Cldn4"], sub["Tscore"])
        nei = None
        if sub["x"].notna().sum() >= 50:
            nei_t = _neighbor_mean(sub["x"].to_numpy(), sub["y"].to_numpy(), sub["Tscore"].to_numpy(), k=6)
            nei = spearman(sub["Cldn4"], nei_t)
        rows.append(
            {
                "gse": "GSE261890",
                "sample": sample,
                "n_bins": int(len(sub)),
                "Cldn4_mean": float(sub["Cldn4"].mean()),
                "Cldn4_pct": float((sub["Cldn4"] > 0).mean() * 100),
                "Tscore_mean": float(sub["Tscore"].mean()),
                "same_rho": rho["rho"],
                "same_p": rho["p"],
                "nei_rho": None if nei is None else nei["rho"],
                "nei_p": None if nei is None else nei["p"],
                "nei_n": None if nei is None else nei["n"],
                **{f"cut_{k}": v for k, v in meta.items()},
            }
        )
    return pd.DataFrame(rows)


def per_series_table(scrna_inv, contrasts, scuts, spcuts):
    """One row per series for the PR / FINDING."""
    rows = []
    # GSE297632
    a = scrna_inv[scrna_inv["sample"] == "LLC_Control"].iloc[0]
    b = scrna_inv[scrna_inv["sample"] == "LLC_aPD1_residual"].iloc[0]
    rows.append(
        {
            "series": "GSE297632",
            "platform": "scRNA 10x (GEM-X Flex)",
            "model": "LLC subcutaneous; residual 7d after anti-PD-1 vs untreated",
            "contrast": "aPD-1 residual vs control",
            "n": "1 vs 1 library",
            "ICB_vs_true_control": "yes (n=1); subcutaneous not orthotopic lung",
            "epithelial_present": "yes",
            "epi_n_cells": f"{int(a.n_epithelial)} vs {int(b.n_epithelial)}",
            "Cldn4_epi_mean": f"{a.epithelial_Cldn4_mean:.4f} vs {b.epithelial_Cldn4_mean:.4f}",
            "Cldn4_epi_pct": f"{a.epithelial_Cldn4_pct:.1f}% vs {b.epithelial_Cldn4_pct:.1f}%",
            "TNK_frac": f"{a.frac_tnk:.3f} vs {b.frac_tnk:.3f}",
            "Cldn4_direction": "up in residual (descriptive)",
            "T_direction": "T/NK fraction down in residual (descriptive)",
            "exclusion_note": "same-library Cldn4 is higher in epi than T/NK; n=1 vs 1 so no p",
            "welch_p": "NA (n<2)",
            "cut_Cldn4pos_Tlow": "see FINDING epithelial detected-cut",
        }
    )
    # GSE261890
    v = spcuts[spcuts["sample"] == "Sen.Res.Veh"].iloc[0]
    i = spcuts[spcuts["sample"] == "Sen.Res.ICI"].iloc[0]
    rows.append(
        {
            "series": "GSE261890",
            "platform": "BMKMANU S1000 spatial, L7 bins",
            "model": "mouse NSCLC ICI; public titles Sen.Res.Veh vs Sen.Res.ICI",
            "contrast": "Veh vs ICI (not a 4-group sen/res split)",
            "n": "1 vs 1 section",
            "ICB_vs_true_control": "Veh vs ICI; overall-design lists 4 groups, only 2 GSM deposited",
            "epithelial_present": "spatial bins (Epcam/Krt scored)",
            "epi_n_cells": f"{int(v.n_bins)} vs {int(i.n_bins)} L7 bins",
            "Cldn4_epi_mean": f"{v.Cldn4_mean:.4f} vs {i.Cldn4_mean:.4f} (all bins)",
            "Cldn4_epi_pct": f"{v.Cldn4_pct:.1f}% vs {i.Cldn4_pct:.1f}%",
            "TNK_frac": f"Tscore {v.Tscore_mean:.4f} vs {i.Tscore_mean:.4f}",
            "Cldn4_direction": "higher on ICI section (descriptive)",
            "T_direction": "T score also higher on ICI section (descriptive)",
            "exclusion_note": f"same-bin ρ(Cldn4,T) Veh {v.same_rho:.3f} / ICI {i.same_rho:.3f} (positive, not exclusion)",
            "welch_p": "NA (n<2)",
            "cut_Cldn4pos_Tlow": f"Veh {int(v.cut_n_hi_lo)}/{int(v.cut_n_cldn4_pos)} Cldn4+; ICI {int(i.cut_n_hi_lo)}/{int(i.cut_n_cldn4_pos)}",
        }
    )
    # GSE285606
    wt = scrna_inv[scrna_inv.arm == "WT_IgG"]
    r = scrna_inv[scrna_inv.arm == "PD1R1_IgG"]
    c = contrasts[(contrasts.gse == "GSE285606") & (contrasts.field == "epithelial_Cldn4_mean")].iloc[0]
    ct = contrasts[(contrasts.gse == "GSE285606") & (contrasts.field == "frac_tnk")].iloc[0]
    rows.append(
        {
            "series": "GSE285606",
            "platform": "scRNA 10x",
            "model": "Kras/p53 344SQ parental vs PD1R1 (140P) acquired aPD-1-resistant",
            "contrast": "PD1R1 vs WT; both IgG",
            "n": "2 vs 2 libraries (T1=wk4, T2=wk6)",
            "ICB_vs_true_control": "no — both arms IgG",
            "epithelial_present": "yes",
            "epi_n_cells": f"WT {int(wt.n_epithelial.sum())} / PD1R1 {int(r.n_epithelial.sum())} cells; n_lib=2 vs 2",
            "Cldn4_epi_mean": f"WT {wt.epithelial_Cldn4_mean.mean():.4f} vs PD1R1 {r.epithelial_Cldn4_mean.mean():.4f}",
            "Cldn4_epi_pct": f"WT {wt.epithelial_Cldn4_pct.mean():.1f}% vs PD1R1 {r.epithelial_Cldn4_pct.mean():.1f}%",
            "TNK_frac": f"WT {wt.frac_tnk.mean():.3f} vs PD1R1 {r.frac_tnk.mean():.3f}",
            "Cldn4_direction": "down in PD1R1 (resistant line)",
            "T_direction": "T/NK fraction up in PD1R1",
            "exclusion_note": "opposite a Cldn4-high / T-low pattern; genotype on IgG, not ICB vs control",
            "welch_p": f"Cldn4 {c.welch_p:.2f}; T/NK frac {ct.welch_p:.2f} (n=2 vs 2)",
            "cut_Cldn4pos_Tlow": "see FINDING epithelial detected-cut",
        }
    )
    # skip
    rows.append(
        {
            "series": "GSE303943",
            "platform": "scRNA 10x",
            "model": "CMT167R subcutaneous",
            "contrast": "PKCi vs solvent",
            "n": "skip",
            "ICB_vs_true_control": "no — not ICB",
            "epithelial_present": "—",
            "epi_n_cells": "—",
            "Cldn4_epi_mean": "—",
            "Cldn4_epi_pct": "—",
            "TNK_frac": "—",
            "Cldn4_direction": "skipped",
            "T_direction": "skipped",
            "exclusion_note": "subcutaneous PKCi, not ICB (given with PR #285)",
            "welch_p": "—",
            "cut_Cldn4pos_Tlow": "—",
        }
    )
    return pd.DataFrame(rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    inv = pd.read_csv(OUT / "scrna_sample_inventory.tsv", sep="\t")
    contrasts = pd.read_csv(OUT / "scrna_contrasts.tsv", sep="\t")
    scuts = scrna_cuts()
    spcuts = spatial_cuts()
    scuts.to_csv(OUT / "scrna_detected_cuts.tsv", sep="\t", index=False)
    spcuts.to_csv(OUT / "spatial_detected_cuts.tsv", sep="\t", index=False)
    table = per_series_table(inv, contrasts, scuts, spcuts)
    table.to_csv(OUT / "per_series_table.tsv", sep="\t", index=False)
    dump_json(OUT / "per_series_table.json", table.to_dict(orient="records"))
    print(table[["series", "n", "Cldn4_direction", "T_direction", "welch_p"]].to_string(index=False))
    print("spatial cuts:\n", spcuts.to_string(index=False))
    print("scrna cuts:\n", scuts.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
