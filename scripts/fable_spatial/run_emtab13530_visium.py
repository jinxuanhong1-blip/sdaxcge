"""E-MTAB-13530 (10x Visium, human NSCLC + non-involved lung + healthy donors).

Per section: TACSTD2/CLDN4 in epithelial spots vs immune score of the
surrounding spot neighborhood (self excluded). Cross-section inference:
Wilcoxon signed-rank on per-section rho + signed Stouffer combination,
separately for tumour (T), non-involved background (B) and donor (D) sections.
Also compares epithelial-spot TACSTD2/CLDN4 between tumour and matched
non-involved sections per patient.
"""

import glob
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(__file__))
from common import (TARGETS, analyze_section, collapse_duplicate_genes, ensure_results,
                    load_10x_h5, load_positions_from_tar, score, stouffer, wilcoxon_rhos,
                    EPI_MARKERS)

DATA = "/workspace/data/fable_spatial/emtab13530"
OUT = ensure_results()


def tissue_class(section):
    if section.startswith("D"):
        return "donor"
    return {"T": "tumour", "B": "non_involved"}[re.match(r"P\d+_([TB])", section).group(1)]


def main():
    h5s = sorted(glob.glob(os.path.join(DATA, "*-filtered_feature_bc_matrix.h5")))
    all_rows, qcs, expr_rows = [], [], []
    for h5 in h5s:
        section = os.path.basename(h5).replace("-filtered_feature_bc_matrix.h5", "")
        tar = os.path.join(DATA, f"{section}-spatial.tar")
        if not os.path.exists(tar):
            print(f"skip {section}: no spatial tar")
            continue
        X, bcs, names = load_10x_h5(h5)
        pos = load_positions_from_tar(tar)
        rows, qc = analyze_section(X, bcs, names, pos, section)
        qc["tissue"] = tissue_class(section)
        qcs.append(qc)
        for r in rows:
            r["tissue"] = qc["tissue"]
        all_rows.extend(rows)

        # epithelial-spot mean expression of targets, for T-vs-B contrast
        if rows:
            Xc, namesc = collapse_duplicate_genes(X, names)
            posc = pos[pos.in_tissue == 1].set_index("barcode")
            keep = [i for i, b in enumerate(bcs) if b in posc.index]
            Xc = Xc[keep]
            tot = np.asarray(Xc.sum(1)).ravel()
            ok = tot >= 250
            Xc = Xc[ok]
            logX = Xc.multiply(1e4 / np.asarray(Xc.sum(1)).ravel()[:, None]).tocsr()
            logX.data = np.log1p(logX.data)
            idx = {g: i for i, g in enumerate(namesc)}
            epi, _ = score(logX, idx, EPI_MARKERS)
            mask = epi > np.median(epi)
            rec = dict(section=section, tissue=qc["tissue"],
                       patient=section.split("_")[0], n_epi=int(mask.sum()))
            for g in TARGETS:
                v = np.asarray(logX[:, idx[g]].todense()).ravel()[mask]
                rec[f"{g}_mean_epi"] = float(v.mean())
            expr_rows.append(rec)
        print(f"{section} [{qc['tissue']}]: {qc.get('n_spots')} spots, "
              f"{qc.get('n_epi_spots', 0)} epi spots, {len(rows)} result rows")

    res = pd.DataFrame(all_rows)
    res.to_csv(os.path.join(OUT, "emtab13530_per_section.csv"), index=False)
    pd.DataFrame(qcs).to_csv(os.path.join(OUT, "emtab13530_section_qc.csv"), index=False)
    ex = pd.DataFrame(expr_rows)
    ex.to_csv(os.path.join(OUT, "emtab13530_epi_expression.csv"), index=False)

    # ------------------------------------------------- cross-section meta --
    meta = []
    for (tissue, gene, nb), g in res.groupby(["tissue", "gene", "neighborhood"]):
        w, wp = wilcoxon_rhos(g.partial_rho)
        zc, zp = stouffer(g.partial_p, np.sign(g.partial_rho))
        meta.append(dict(dataset="E-MTAB-13530", tissue=tissue, gene=gene,
                         neighborhood=nb, n_sections=len(g),
                         median_partial_rho=round(float(g.partial_rho.median()), 4),
                         n_positive=int((g.partial_rho > 0).sum()),
                         wilcoxon_p=wp, stouffer_z=zc, stouffer_p=zp))
    meta = pd.DataFrame(meta)
    meta.to_csv(os.path.join(OUT, "emtab13530_meta.csv"), index=False)
    print("\n== meta across sections (partial Spearman) ==")
    print(meta.to_string(index=False))

    # -------------------------------------- tumour vs non-involved, paired --
    rows = []
    pt = ex.groupby(["patient", "tissue"])[[f"{g}_mean_epi" for g in TARGETS]].mean().reset_index()
    for g in TARGETS:
        wide = pt.pivot(index="patient", columns="tissue", values=f"{g}_mean_epi").dropna(subset=["tumour", "non_involved"])
        try:
            w, p = stats.wilcoxon(wide["tumour"], wide["non_involved"])
        except ValueError:
            w, p = np.nan, np.nan
        rows.append(dict(gene=g, n_patients=len(wide),
                         mean_tumour=round(float(wide["tumour"].mean()), 3),
                         mean_non_involved=round(float(wide["non_involved"].mean()), 3),
                         wilcoxon_p=float(p)))
    tb = pd.DataFrame(rows)
    tb.to_csv(os.path.join(OUT, "emtab13530_tumour_vs_noninvolved.csv"), index=False)
    print("\n== epithelial-spot expression: tumour vs matched non-involved ==")
    print(tb.to_string(index=False))


if __name__ == "__main__":
    main()
