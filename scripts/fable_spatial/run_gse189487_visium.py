"""GSE189487 (10x Visium, LUAD progression AIS -> MIA -> IAC, 6 sections).

Same neighborhood analysis as E-MTAB-13530, plus a stage trend test
(Jonckheere-type via Spearman of stage rank vs epithelial-spot expression,
patient-level, n=6 sections).
"""

import glob
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(__file__))
from common import (TARGETS, analyze_section, collapse_duplicate_genes, ensure_results,
                    load_mtx_triplet, load_positions_csv_gz, score, stouffer,
                    wilcoxon_rhos, EPI_MARKERS)

DATA = "/workspace/data/fable_spatial/gse189487"
OUT = ensure_results()

STAGE = {"TD1": "IAC", "TD2": "IAC", "TD3": "MIA", "TD5": "AIS", "TD6": "MIA", "TD8": "AIS"}
STAGE_RANK = {"AIS": 0, "MIA": 1, "IAC": 2}


def main():
    mtxs = sorted(glob.glob(os.path.join(DATA, "GSM*_matrix.mtx.gz")))
    all_rows, qcs, expr_rows = [], [], []
    for mtx in mtxs:
        prefix = mtx.replace("_matrix.mtx.gz", "")
        section = os.path.basename(prefix).split("_")[1]
        X, bcs, names = load_mtx_triplet(mtx, prefix + "_features.tsv.gz",
                                         prefix + "_barcodes.tsv.gz")
        pos = load_positions_csv_gz(prefix + "_tissue_positions_list.csv.gz")
        rows, qc = analyze_section(X, bcs, names, pos, section)
        qc["stage"] = STAGE[section]
        qcs.append(qc)
        for r in rows:
            r["stage"] = qc["stage"]
        all_rows.extend(rows)

        Xc, namesc = collapse_duplicate_genes(X, names)
        posc = pos[pos.in_tissue == 1].set_index("barcode")
        keep = [i for i, b in enumerate(bcs) if b in posc.index]
        Xc = Xc[keep]
        tot = np.asarray(Xc.sum(1)).ravel()
        Xc = Xc[tot >= 250]
        logX = Xc.multiply(1e4 / np.asarray(Xc.sum(1)).ravel()[:, None]).tocsr()
        logX.data = np.log1p(logX.data)
        idx = {g: i for i, g in enumerate(namesc)}
        epi, _ = score(logX, idx, EPI_MARKERS)
        mask = epi > np.median(epi)
        rec = dict(section=section, stage=qc["stage"], n_epi=int(mask.sum()))
        for g in TARGETS:
            v = np.asarray(logX[:, idx[g]].todense()).ravel()[mask]
            rec[f"{g}_mean_epi"] = float(v.mean())
        expr_rows.append(rec)
        print(f"{section} [{qc['stage']}]: {qc.get('n_spots')} spots, "
              f"{qc.get('n_epi_spots', 0)} epi spots")

    res = pd.DataFrame(all_rows)
    res.to_csv(os.path.join(OUT, "gse189487_per_section.csv"), index=False)
    pd.DataFrame(qcs).to_csv(os.path.join(OUT, "gse189487_section_qc.csv"), index=False)
    ex = pd.DataFrame(expr_rows)
    ex.to_csv(os.path.join(OUT, "gse189487_epi_expression.csv"), index=False)

    meta = []
    for (gene, nb), g in res.groupby(["gene", "neighborhood"]):
        w, wp = wilcoxon_rhos(g.partial_rho)
        zc, zp = stouffer(g.partial_p, np.sign(g.partial_rho))
        meta.append(dict(dataset="GSE189487", gene=gene, neighborhood=nb,
                         n_sections=len(g),
                         median_partial_rho=round(float(g.partial_rho.median()), 4),
                         n_positive=int((g.partial_rho > 0).sum()),
                         wilcoxon_p=wp, stouffer_z=zc, stouffer_p=zp))
    meta = pd.DataFrame(meta)
    meta.to_csv(os.path.join(OUT, "gse189487_meta.csv"), index=False)
    print("\n== meta across sections ==")
    print(meta.to_string(index=False))

    print("\n== stage trend (epithelial-spot mean, n=6 sections) ==")
    trend = []
    ex["stage_rank"] = ex["stage"].map(STAGE_RANK)
    for g in TARGETS:
        r, p = stats.spearmanr(ex["stage_rank"], ex[f"{g}_mean_epi"])
        trend.append(dict(gene=g, spearman_rho=round(float(r), 3), p=float(p), n=len(ex)))
    trend = pd.DataFrame(trend)
    trend.to_csv(os.path.join(OUT, "gse189487_stage_trend.csv"), index=False)
    print(trend.to_string(index=False))
    print(ex.to_string(index=False))


if __name__ == "__main__":
    main()
