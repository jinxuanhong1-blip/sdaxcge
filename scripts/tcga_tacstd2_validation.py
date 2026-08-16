#!/usr/bin/env python3
"""
Open-RNA validation: TACSTD2 (TROP2) vs immune signatures in TCGA LUAD/LUSC
(PanCanAtlas RSEM via cBioPortal API), purity-adjusted with ABSOLUTE purity
(TCGA_mastercalls.abs_tables_JSedit.fixed.txt from the official GDC
PanCanAtlas publication page, file UUID 4f277128-f793-4354-a13d-30cc7fe9f6b5).

ABSOLUTE purity is DNA-derived and independent of expression, so the partial
correlation here has no ESTIMATE-style circularity. TCGA is treatment-naive
resected NSCLC — a biological reference for the PACIFIC population, not
durvalumab-exposed tissue.
"""

import json
import os
import urllib.request

import numpy as np
import pandas as pd
from scipy import stats

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(BASE, "data", "hunt_pacific")
OUT = os.path.join(BASE, "results", "hunt_pacific")
API = "https://www.cbioportal.org/api"

GENES = {"CCL5": 6352, "CD27": 939, "CD274": 29126, "CD276": 80381,
         "CD8A": 925, "CMKLR1": 1240, "CXCL10": 3627, "CXCL9": 4283,
         "CXCR6": 10663, "GZMA": 3001, "HLA-DQA1": 3117, "HLA-DRA": 3122,
         "HLA-DRB1": 3123, "HLA-E": 3133, "IDO1": 3620, "IFNG": 3458,
         "LAG3": 3902, "NKG7": 4818, "PDCD1LG2": 80380, "PRF1": 5551,
         "PSMB10": 5699, "STAT1": 6772, "TACSTD2": 4070, "TIGIT": 201633}
E2S = {v: k for k, v in GENES.items()}

IFNG6 = ["IFNG", "STAT1", "IDO1", "CXCL10", "CXCL9", "HLA-DRA"]
GEP18 = ["CCL5", "CD27", "CD274", "CD276", "CD8A", "CMKLR1", "CXCL9",
         "CXCR6", "HLA-DQA1", "HLA-DRB1", "HLA-E", "IDO1", "LAG3", "NKG7",
         "PDCD1LG2", "PSMB10", "STAT1", "TIGIT"]
CYT = ["GZMA", "PRF1"]


def post(url, payload):
    req = urllib.request.Request(url, data=json.dumps(payload).encode(),
                                 headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.loads(r.read())


def fetch_expression(study):
    profile = f"{study}_rna_seq_v2_mrna"
    rows = post(f"{API}/molecular-profiles/{profile}/molecular-data/fetch"
                "?projection=SUMMARY",
                {"sampleListId": f"{study}_all",
                 "entrezGeneIds": sorted(GENES.values())})
    df = pd.DataFrame([(r["sampleId"], E2S[r["entrezGeneId"]], r["value"])
                       for r in rows], columns=["sample", "gene", "value"])
    return df.pivot_table(index="sample", columns="gene", values="value")


def partial_spearman(x, y, z):
    rx, ry, rz = (stats.rankdata(v) for v in (x, y, z))
    def resid(a, b):
        b1 = np.column_stack([np.ones_like(b), b])
        beta, *_ = np.linalg.lstsq(b1, a, rcond=None)
        return a - b1 @ beta
    ex, ey = resid(rx, rz), resid(ry, rz)
    r, _ = stats.pearsonr(ex, ey)
    dof = len(x) - 3
    t = r * np.sqrt(dof / (1 - r ** 2))
    return r, 2 * stats.t.sf(abs(t), dof)


def main():
    pur = pd.read_csv(os.path.join(
        DATA, "TCGA_mastercalls.abs_tables_JSedit.fixed.txt"), sep="\t")
    pur = pur.dropna(subset=["purity"]).set_index("array")["purity"]

    results, log = [], []
    for study, label in [("luad_tcga_pan_can_atlas_2018", "TCGA-LUAD"),
                         ("lusc_tcga_pan_can_atlas_2018", "TCGA-LUSC")]:
        expr = np.log2(fetch_expression(study) + 1)
        merged = expr.join(pur.rename("ABSOLUTE_purity"), how="inner").dropna(
            subset=["TACSTD2", "ABSOLUTE_purity"])
        log.append(f"== {label}: {len(expr)} expr samples, "
                   f"{len(merged)} with ABSOLUTE purity ==")
        merged["IFNG6_Ayers"] = merged[IFNG6].mean(axis=1)
        merged["GEP18_TcellInflamed"] = merged[GEP18].mean(axis=1)
        merged["CYT_Rooney"] = merged[CYT].mean(axis=1)
        merged.to_csv(os.path.join(OUT, f"scores_{label}.csv"))

        x = merged["TACSTD2"].values
        z = merged["ABSOLUTE_purity"].values
        for name in ["IFNG6_Ayers", "GEP18_TcellInflamed", "CYT_Rooney",
                     "CD8A"]:
            y = merged[name].values
            rho, p = stats.spearmanr(x, y)
            prho, pp = partial_spearman(x, y, z)
            results.append({"dataset": label, "immune_score": name,
                            "n": len(x), "spearman_rho": round(rho, 4),
                            "spearman_p": float(f"{p:.2e}"),
                            "purity_adj_rho": round(prho, 4),
                            "purity_adj_p": float(f"{pp:.2e}")})
            log.append(f"  TACSTD2 vs {name}: rho={rho:.4f} (p={p:.2e}), "
                       f"ABSOLUTE-purity-adj rho={prho:.4f} (p={pp:.2e})")
        rho_z, p_z = stats.spearmanr(x, z)
        log.append(f"  TACSTD2 vs ABSOLUTE purity: rho={rho_z:.4f} (p={p_z:.2e})")
        results.append({"dataset": label, "immune_score": "ABSOLUTE_purity",
                        "n": len(x), "spearman_rho": round(rho_z, 4),
                        "spearman_p": float(f"{p_z:.2e}"),
                        "purity_adj_rho": np.nan, "purity_adj_p": np.nan})

    pd.DataFrame(results).to_csv(
        os.path.join(OUT, "tacstd2_immune_correlations_tcga.csv"), index=False)
    with open(os.path.join(OUT, "analysis_log_tcga.txt"), "w") as fh:
        fh.write("\n".join(log) + "\n")
    print("\n".join(log))


if __name__ == "__main__":
    main()
