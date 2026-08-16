"""Aggregate per-section tables into a machine-readable summary.json."""
import json
import numpy as np
import pandas as pd
from scipy import stats
import common

OUT = "results/claim_B6"


def wilcox(x):
    x = np.asarray([v for v in x if np.isfinite(v)], float)
    if len(x) < 3:
        return dict(n=len(x), median=float(np.median(x)) if len(x) else None,
                    n_neg=int(np.sum(x < 0)) if len(x) else 0,
                    n_pos=int(np.sum(x > 0)) if len(x) else 0,
                    wilcoxon_p=None, sign_p=None)
    try:
        p_w = float(stats.wilcoxon(x).pvalue)
    except ValueError:
        p_w = None
    s = common.sign_test_summary(x)
    return dict(n=int(len(x)), median=float(np.median(x)), mean=float(np.mean(x)),
                n_neg=s["n_neg"], n_pos=s["n_pos"],
                wilcoxon_p=p_w, sign_p=s["sign_p"])


def main():
    em = pd.read_csv(f"{OUT}/tables/emtab13530_per_section.csv")
    g189 = pd.read_csv(f"{OUT}/tables/gse189487_per_section.csv")
    vis = pd.concat([em, g189], ignore_index=True)
    vis.to_csv(f"{OUT}/tables/visium_samespot_all.csv", index=False)

    nb = pd.read_csv(f"{OUT}/tables/visium_neighborhood_rings.csv")
    like = pd.read_csv(f"{OUT}/tables/visium_pr110like_rings.csv")
    gx = pd.read_csv(f"{OUT}/tables/gse271689_geomx_correlations.csv")
    cx_nb = pd.read_csv(f"{OUT}/tables/gse287472_cosmx_neighborhood.csv")
    cx_same = pd.read_csv(f"{OUT}/tables/gse287472_cosmx_samecell.csv")

    summary = {
        "claim": "B6: CLDN4/TACSTD2-high spots anti-colocalize with T/B cells",
        "same_spot_visium": {},
        "neighborhood_visium": {},
        "pr110_like": {},
        "geomx_gse271689": gx.to_dict(orient="records"),
        "cosmx_gse287472": {
            "note": "1 patient, 1000-plex CosMx; CLDN4 detected in 3.6% of cells",
            "same_cell": cx_same.to_dict(orient="records"),
            "neighborhood": cx_nb.to_dict(orient="records"),
        },
    }

    for (ds, grp), sub in vis.groupby(["dataset", "group"]):
        key = f"{ds}|{grp}"
        summary["same_spot_visium"][key] = {
            "n_sections": int(len(sub)),
            "rho_epi_vs_TB": wilcox(sub["rho_epi_vs_TB"]),
            "rho_epi_vs_T": wilcox(sub["rho_epi_vs_T"]),
            "rho_epi_vs_B": wilcox(sub["rho_epi_vs_B"]),
            "rb_TB_in_epihigh": wilcox(sub["rb_TB_in_epihigh"]),
        }

    for (ds, grp, ring), sub in nb.groupby(["dataset", "group", "ring"]):
        key = f"{ds}|{grp}|ring{int(ring)}"
        summary["neighborhood_visium"][key] = {
            "n_sections": int(len(sub)),
            "rho_CLDN4_vs_TBnb": wilcox(sub["rho_CLDN4_vs_TBnb"]),
            "rho_CLDN4_vs_Tnb": wilcox(sub["rho_CLDN4_vs_Tnb"]),
            "rho_CLDN4_vs_Bnb": wilcox(sub["rho_CLDN4_vs_Bnb"]),
            "rho_partialEpi": wilcox(sub["rho_CLDN4_vs_TBnb_partialEpi"]),
        }

    for (ds, ring), sub in like.groupby(["dataset", "ring"]):
        key = f"{ds}|ring{int(ring)}"
        summary["pr110_like"][key] = {
            "n_sections": int(len(sub)),
            "rho_raw": wilcox(sub["rho_raw"]),
            "rho_partialEpi": wilcox(sub["rho_partialEpi"]),
        }

    # overall tumor visium same-spot
    tumor = vis[vis["group"] == "tumor"]
    summary["same_spot_visium_tumor_pooled"] = {
        "n_sections": int(len(tumor)),
        "rho_epi_vs_TB": wilcox(tumor["rho_epi_vs_TB"]),
        "rho_epi_vs_T": wilcox(tumor["rho_epi_vs_T"]),
        "rho_epi_vs_B": wilcox(tumor["rho_epi_vs_B"]),
    }

    with open(f"{OUT}/summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(summary["same_spot_visium"], indent=2))
    print("--- neighborhood tumor ring1 ---")
    for k, v in summary["neighborhood_visium"].items():
        if "tumor|ring1" in k:
            print(k, json.dumps(v, indent=2))
    print("--- pr110-like ---")
    print(json.dumps(summary["pr110_like"], indent=2))


if __name__ == "__main__":
    main()
