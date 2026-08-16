#!/usr/bin/env python3
"""Large-n test of the tumor-intrinsic hypothesis in TCGA-LUAD/LUSC and OncoSG.

These are not ICI cohorts, but they provide the statistical power to test the
mechanistic backbone of the user's direction:
  * does tumor TACSTD2 (TROP2) ANTI-correlate with CD8 / NK cytotoxic signatures
    after correcting for tumor purity (epithelial content)?

We fetch expression for the target + signature genes from cBioPortal, build the
same signatures as the ICI cohorts, and compute raw and purity-corrected
(partial Spearman, controlling epithelial content) correlations.
"""
import json
import os
import sys
import time
import urllib.request

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
import signatures as S

DATA = "/tmp/ici_bulk_data"
RES = "/workspace/results/fable_ici_bulk"
API = "https://www.cbioportal.org/api"

STUDIES = [
    ("TCGA_LUAD", "luad_tcga_pan_can_atlas_2018",
     "luad_tcga_pan_can_atlas_2018_rna_seq_v2_mrna", "log2p1"),
    ("TCGA_LUSC", "lusc_tcga_pan_can_atlas_2018",
     "lusc_tcga_pan_can_atlas_2018_rna_seq_v2_mrna", "log2p1"),
    ("OncoSG_LUAD", "luad_oncosg_2020",
     "luad_oncosg_2020_rna_seq_v2_mrna_median_Zscores", "asis"),
]

ALL_GENES = sorted(set(S.TARGET_GENES + sum(S.SIGNATURES.values(), [])))


def _post(url, body):
    for i in range(4):
        try:
            req = urllib.request.Request(
                url, data=json.dumps(body).encode(),
                headers={"Content-Type": "application/json",
                         "Accept": "application/json"}, method="POST")
            return json.load(urllib.request.urlopen(req, timeout=120))
        except Exception as e:  # noqa
            print("retry", i, e)
            time.sleep(2 ** i)
    raise RuntimeError("POST failed: " + url)


def fetch_expression(study, profile, entrez):
    url = f"{API}/molecular-profiles/{profile}/molecular-data/fetch?projection=SUMMARY"
    data = _post(url, {"entrezGeneIds": list(entrez.values()),
                       "sampleListId": f"{study}_all"})
    rev = {v: k for k, v in entrez.items()}
    rows = {}
    for d in data:
        g = rev.get(d["entrezGeneId"])
        if g is None:
            continue
        rows.setdefault(g, {})[d["sampleId"]] = d["value"]
    df = pd.DataFrame(rows).T  # genes x samples
    return df


def get_entrez_map():
    path = f"{DATA}/entrez_map.json"
    if os.path.exists(path):
        return json.load(open(path))
    url = f"{API}/genes/fetch?geneIdType=HUGO_GENE_SYMBOL"
    data = _post(url, ALL_GENES)
    m = {g["hugoGeneSymbol"]: g["entrezGeneId"] for g in data}
    os.makedirs(DATA, exist_ok=True)
    json.dump(m, open(path, "w"))
    return m


def main():
    entrez = get_entrez_map()
    corr_rows = []
    for name, study, profile, scale in STUDIES:
        print(f"Fetching {name} ...")
        expr = fetch_expression(study, profile, entrez)
        expr = expr.dropna(how="all")
        if scale == "log2p1":
            expr = np.log2(expr.clip(lower=0) + 1)
        # (z-score profiles used as-is; Spearman is rank-based)
        scores = {}
        found = {}
        for signame, genes in S.SIGNATURES.items():
            sc, fg = S.signature_score(expr, genes)
            scores[signame] = sc
            found[signame] = fg
        epi = scores["Epithelial"]
        leuk = scores["Leukocyte"]
        n_samples = expr.shape[1]
        for gene in S.TARGET_GENES:
            if gene not in expr.index:
                continue
            target = expr.loc[gene]
            for sig in ["CD8_Tcell", "NK_cell", "Cytotoxic", "CD8_NK"]:
                sv = scores[sig]
                r_raw, p_raw, n = S.spearman(target.values, sv.reindex(target.index).values)
                r_epi, p_epi, _ = S.partial_spearman(
                    target.values, sv.reindex(target.index).values,
                    [epi.reindex(target.index).values])
                r_both, p_both, _ = S.partial_spearman(
                    target.values, sv.reindex(target.index).values,
                    [epi.reindex(target.index).values,
                     leuk.reindex(target.index).values])
                corr_rows.append({
                    "dataset": name, "n_samples": n_samples, "gene": gene,
                    "immune_signature": sig, "n": n,
                    "spearman_raw": r_raw, "p_raw": p_raw,
                    "partial_epi": r_epi, "p_partial_epi": p_epi,
                    "partial_epi_leuk": r_both, "p_partial_epi_leuk": p_both,
                })
        print(f"  {name}: {n_samples} samples, {expr.shape[0]} genes fetched")

    cdf = pd.DataFrame(corr_rows)
    cdf.to_csv(f"{RES}/tables/cbioportal_correlations.csv", index=False)
    pd.set_option("display.width", 240)
    pd.set_option("display.max_columns", 40)
    print("\n=== TCGA/OncoSG: target vs immune, raw & purity-corrected ===")
    print(cdf.to_string(index=False))


if __name__ == "__main__":
    main()
