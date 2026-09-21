#!/usr/bin/env python3
"""Second-pass search: smoothed epithelial scores and epithelial-bin fields."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

from cosmx_hotspot_hmrf_modules import (
    ADHESION,
    EPI_MARKERS,
    IMM_MARKERS,
    KERATIN,
    MIN_CD8,
    MIN_EPI,
    PATIENT,
    PATIENT_ORDER,
    PX_TO_UM,
    SAMPLES,
    STR_MARKERS,
    library_lognorm,
    load_sample,
    mean_genes,
)

OUT = Path(__file__).resolve().parents[1] / "results" / "cosmx_hotspot_hmrf"
PROGRAM = ["CLDN4", "EPCAM", "TACSTD2", "CDH1", "KRT19", "KRT8", "KRT7", "EZR"]


def spearman(a, b, min_n=8) -> float:
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    m = np.isfinite(a) & np.isfinite(b)
    a, b = a[m], b[m]
    if a.size < min_n or np.unique(a).size < 2 or np.unique(b).size < 2:
        return np.nan
    return float(spearmanr(a, b).statistic)


def local_mean(xy, values, radius_px):
    tree = cKDTree(xy)
    groups = tree.query_ball_point(xy, r=radius_px)
    out = np.empty(len(values), dtype=np.float64)
    for i, ix in enumerate(groups):
        out[i] = float(np.mean(values[ix])) if ix else np.nan
    return out


def main():
    rows = []
    for sample in SAMPLES:
        pack = load_sample(sample)
        gmap = pack["gmap"]
        logx = library_lognorm(pack["counts"])
        stacks = np.vstack([
            mean_genes(logx, gmap, EPI_MARKERS),
            mean_genes(logx, gmap, IMM_MARKERS),
            mean_genes(logx, gmap, STR_MARKERS),
        ])
        comp = np.array(["epithelial", "immune", "stromal"])[stacks.argmax(0)]
        comp[stacks.max(0) <= 0] = "unassigned"
        is_epi = comp == "epithelial"
        cd8 = np.zeros(len(comp))
        cd3 = np.zeros(len(comp))
        for g in ("CD8A", "CD8B"):
            cd8 += pack["counts"][:, gmap[g]]
        for g in ("CD3D", "CD3E", "CD3G"):
            cd3 += pack["counts"][:, gmap[g]]
        is_cd8 = (cd8 > 0) & (cd3 > 0) & (~is_epi)
        cldn4 = logx[:, gmap["CLDN4"]]
        program = mean_genes(logx, gmap, PROGRAM)
        keratin = mean_genes(logx, gmap, KERATIN)
        adhesion = mean_genes(logx, gmap, ADHESION)
        patient = PATIENT[sample]
        for fov in np.unique(pack["fov"]):
            m = pack["fov"] == fov
            if int((m & is_epi).sum()) < MIN_EPI or int((m & is_cd8).sum()) < MIN_CD8:
                continue
            xy = pack["xy"][m]
            epi = is_epi[m]
            cd8m = is_cd8[m]
            epi_xy = xy[epi]
            cl = cldn4[m][epi]
            prog = program[m][epi]
            ker = keratin[m][epi]
            ad = adhesion[m][epi]
            tree_cd8 = cKDTree(xy[cd8m])
            smooth_r = 25.0 / PX_TO_UM
            cl_s = local_mean(epi_xy, cl, smooth_r)
            prog_s = local_mean(epi_xy, prog, smooth_r)
            ad_s = local_mean(epi_xy, ad, smooth_r)
            preds = {
                "cldn4": cl,
                "cldn4_smooth25": cl_s,
                "program": prog,
                "program_smooth25": prog_s,
                "adhesion_smooth25": ad_s,
                "keratin": ker,
            }
            # Binary high/low on smoothed CLDN4 and on smoothed program.
            for key in ("cldn4_smooth25", "program_smooth25"):
                x = preds[key]
                med = np.median(x)
                code = np.where(x > med, 1.0, np.where(x < med, -1.0, 0.0))
                preds[key + "_highlow"] = code
            q1, q2 = np.quantile(cl, [1 / 3, 2 / 3])
            extremes = (cl <= q1) | (cl >= q2)
            for um in (20.0, 30.0, 50.0):
                count = np.asarray(
                    tree_cd8.query_ball_point(epi_xy, r=um / PX_TO_UM, return_length=True),
                    dtype=np.float64,
                )
                for pname, x in preds.items():
                    rows.append((patient, f"{pname}|all|count_{int(um)}", spearman(x, count, 25)))
                    if pname in ("cldn4", "program", "cldn4_smooth25", "program_smooth25"):
                        rows.append((patient, f"{pname}|tertile|count_{int(um)}", spearman(x[extremes], count[extremes], 25)))
            # Epithelial-majority bins.
            for bin_um in (40.0, 60.0, 80.0):
                bin_px = bin_um / PX_TO_UM
                origin = xy.min(0)
                ijk = np.floor((xy - origin) / bin_px).astype(np.int32)
                keys = ijk[:, 0] * 100000 + ijk[:, 1]
                xs, ys, xs_p, ys_k = [], [], [], []
                for key in np.unique(keys):
                    sel = keys == key
                    n_e = int((sel & epi).sum())
                    n = int(sel.sum())
                    if n_e < 8 or n_e / n < 0.5:
                        continue
                    xs.append(float(cldn4[m][sel & epi].mean()))
                    ys.append(float((sel & cd8m).sum()))
                    xs_p.append(float(program[m][sel & epi].mean()))
                    ys_k.append(float(keratin[m][sel & epi].mean()))
                rows.append((patient, f"bin_cldn4|epi50|count_{int(bin_um)}", spearman(xs, ys)))
                rows.append((patient, f"bin_program|epi50|count_{int(bin_um)}", spearman(xs_p, ys)))
                rows.append((patient, f"bin_keratin|epi50|count_{int(bin_um)}", spearman(ys_k, ys)))
        print(sample, flush=True)

    df = pd.DataFrame(rows, columns=["patient", "spec", "rho"])
    summary = []
    for spec, g in df.groupby("spec"):
        means = []
        for patient in PATIENT_ORDER:
            v = g.loc[g.patient == patient, "rho"].to_numpy(float)
            v = v[np.isfinite(v)]
            means.append(float(v.mean()) if v.size else np.nan)
        if not np.all(np.isfinite(means)):
            continue
        summary.append({
            "spec": spec,
            "patient_mean_rho": float(np.mean(means)),
            "n_patients_neg": int(sum(x < 0 for x in means)),
            "n_fov": int(np.isfinite(g["rho"]).sum()),
            **{f"rho_{p}": means[i] for i, p in enumerate(PATIENT_ORDER)},
        })
    sm = pd.DataFrame(summary).sort_values("patient_mean_rho")
    sm.to_csv(OUT / "tables" / "spec_search2_summary.tsv", sep="\t", index=False)
    print(sm.head(30).to_string(index=False))
    print("--- 5/5 ---")
    print(sm[sm.n_patients_neg == 5].head(15).to_string(index=False))


if __name__ == "__main__":
    main()
