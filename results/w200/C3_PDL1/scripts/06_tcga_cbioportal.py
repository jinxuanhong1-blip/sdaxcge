#!/usr/bin/env python3
"""TCGA/cBioPortal CLDN4 vs CD274 Spearman — association only, not C3.

C3 is a causal claim (PD-L1 rises AFTER CLDN4 loss). Bulk co-expression is
confounded by purity (CLDN4 epithelial, CD274 immune + tumor). This table
exists so the association is measured and then explicitly not over-read.

Studies: PanCancer Atlas LUAD, LUSC, BRCA, OV, STAD, COADREAD.
"""

from __future__ import annotations

import json
import os
import urllib.request

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "omics", "tcga_cbioportal")
os.makedirs(OUT, exist_ok=True)

UA = "C3-PDL1-catalog/1.0"
API = "https://www.cbioportal.org/api"

# entrez
GENES = {"CLDN4": 1364, "CD274": 29126, "TACSTD2": 4070, "CD8A": 925, "EPCAM": 4072}

STUDIES = [
    ("luad_tcga_pan_can_atlas_2018", "LUAD"),
    ("lusc_tcga_pan_can_atlas_2018", "LUSC"),
    ("brca_tcga_pan_can_atlas_2018", "BRCA"),
    ("ov_tcga_pan_can_atlas_2018", "OV"),
    ("stad_tcga_pan_can_atlas_2018", "STAD"),
    ("coadread_tcga_pan_can_atlas_2018", "COADREAD"),
]


def get_json(url: str):
    req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=90) as fh:
        return json.load(fh)


def post_json(url: str, payload: dict):
    data = json.dumps(payload).encode()
    req = urllib.request.Request(
        url, data=data,
        headers={"User-Agent": UA, "Accept": "application/json", "Content-Type": "application/json"},
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=120) as fh:
        return json.load(fh)


def main():
    rows = []
    for study, label in STUDIES:
        print(study, flush=True)
        profiles = get_json(f"{API}/studies/{study}/molecular-profiles")
        rna = None
        for p in profiles:
            mid = p.get("molecularProfileId", "")
            if mid.endswith("rna_seq_v2_mrna") and "zscores" not in mid:
                rna = mid
                break
        if rna is None:
            # fallback: any mrna
            for p in profiles:
                if p.get("molecularAlterationType") == "MRNA_EXPRESSION" and "zscore" not in p.get("molecularProfileId", ""):
                    rna = p["molecularProfileId"]
                    break
        if rna is None:
            rows.append({"study": study, "label": label, "error": "no_mrna_profile"})
            continue
        samples = get_json(f"{API}/studies/{study}/samples?projection=ID&pageSize=20000")
        sample_ids = [s["sampleId"] for s in samples]
        # fetch gene expression
        body = {
            "entrezGeneIds": list(GENES.values()),
            "sampleIds": sample_ids,
        }
        data = post_json(f"{API}/molecular-profiles/{rna}/molecular-data/fetch", body)
        # pivot
        recs = {}
        for d in data:
            sid = d.get("sampleId")
            gid = d.get("entrezGeneId")
            val = d.get("value")
            recs.setdefault(sid, {})[gid] = val
        mat = pd.DataFrame.from_dict(recs, orient="index")
        mat.columns = [k for k, v in GENES.items() for e in [v] if e in mat.columns]  # may misorder
        # safer rename
        inv = {v: k for k, v in GENES.items()}
        mat = pd.DataFrame.from_dict(recs, orient="index").rename(columns=inv)
        mat = mat.apply(pd.to_numeric, errors="coerce")
        n = int(mat[["CLDN4", "CD274"]].dropna().shape[0])
        if n < 10:
            rows.append({"study": study, "label": label, "n": n, "error": "too_few"})
            continue
        pair = mat[["CLDN4", "CD274"]].dropna()
        rho, p = stats.spearmanr(pair["CLDN4"], pair["CD274"])
        # tertiles of CLDN4: CD274 in low vs high
        q = pair["CLDN4"].quantile([1/3, 2/3])
        low = pair.loc[pair["CLDN4"] <= q.iloc[0], "CD274"]
        high = pair.loc[pair["CLDN4"] >= q.iloc[1], "CD274"]
        lfc = float(np.log2((low.mean() + 1) / (high.mean() + 1)))
        _, p_t = stats.ttest_ind(low, high, equal_var=False)
        # also TACSTD2 vs CD274
        rho2 = p2 = None
        if "TACSTD2" in mat.columns:
            pair2 = mat[["TACSTD2", "CD274"]].dropna()
            rho2, p2 = stats.spearmanr(pair2["TACSTD2"], pair2["CD274"])
        rows.append({
            "study": study, "label": label, "profile": rna, "n": n,
            "spearman_CLDN4_CD274": float(rho),
            "spearman_p": float(p),
            "CD274_log2_CLDN4low_over_high": lfc,
            "CD274_low_vs_high_welch_p": float(p_t),
            "mean_CD274_CLDN4_low": float(low.mean()),
            "mean_CD274_CLDN4_high": float(high.mean()),
            "spearman_TACSTD2_CD274": float(rho2) if rho2 is not None else None,
            "spearman_TACSTD2_p": float(p2) if p2 is not None else None,
            "note": "bulk mRNA, no purity adjustment; not a loss-of-function experiment",
        })
        mat.to_csv(os.path.join(OUT, f"{label}_genes.tsv"), sep="\t")
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "spearman.tsv"), sep="\t", index=False)
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()
