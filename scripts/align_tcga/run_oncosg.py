"""
USER-ALIGN (OncoSG): repeat the TACSTD2 vs immune-signature analysis in the
public OncoSG lung adenocarcinoma cohort (Chen et al. 2020), which exposes both
mRNA expression and a per-sample tumour PURITY value via cBioPortal.

Expression profile used: `..._rna_seq_v2_mrna_median_all_sample_Zscores`
(per-gene z-scores of log RSEM). Spearman is invariant to the per-gene linear
z-transform, so single-gene results are identical to raw RSEM; signature scores
are the mean of the per-gene z-scores.

All data are public/open via the cBioPortal public API. Outputs to
results/align_tcga/.
"""
from __future__ import annotations

import os
import time
import numpy as np
import pandas as pd
import requests

import common as C

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "..", "results", "align_tcga")
API = "https://www.cbioportal.org/api"
STUDY = "luad_oncosg_2020"
PROFILE = f"{STUDY}_rna_seq_v2_mrna_median_all_sample_Zscores"
SAMPLE_LIST = f"{STUDY}_rna_seq_v2_mrna"


def _get(url, **kw):
    for attempt in range(4):
        try:
            r = requests.get(url, timeout=90, **kw)
            r.raise_for_status()
            return r
        except Exception as e:  # noqa: BLE001
            if attempt == 3:
                raise
            time.sleep(2 ** attempt)


def _post(url, **kw):
    for attempt in range(4):
        try:
            r = requests.post(url, timeout=120, **kw)
            r.raise_for_status()
            return r
        except Exception as e:  # noqa: BLE001
            if attempt == 3:
                raise
            time.sleep(2 ** attempt)


def all_query_genes() -> list[str]:
    genes = {C.TARGET_GENE}
    for gs in C.SIGNATURES.values():
        genes.update(gs)
    genes.update(C.SINGLE_IMMUNE_GENES)
    genes.update(C.TJ_GENES)
    # OncoSG/cBioPortal uses the current symbol VSIR (not C10orf54).
    genes.discard("C10orf54")
    genes.add("VSIR")
    return sorted(genes)


def fetch_expression() -> pd.DataFrame:
    genes = all_query_genes()
    g = _post(f"{API}/genes/fetch", json=genes,
              params={"geneIdType": "HUGO_GENE_SYMBOL"}).json()
    entrez2hugo = {x["entrezGeneId"]: x["hugoGeneSymbol"] for x in g}
    entrez_ids = list(entrez2hugo.keys())

    data = _post(
        f"{API}/molecular-profiles/{PROFILE}/molecular-data/fetch",
        params={"projection": "SUMMARY"},
        json={"entrezGeneIds": entrez_ids, "sampleListId": SAMPLE_LIST},
    ).json()

    recs = [(entrez2hugo[d["entrezGeneId"]], d["sampleId"], d["value"])
            for d in data if d.get("value") is not None]
    df = pd.DataFrame(recs, columns=["gene", "sample", "value"])
    expr = df.pivot_table(index="gene", columns="sample", values="value")
    # map legacy symbol back so common.py resolve works
    if "VSIR" in expr.index:
        expr = expr.rename(index={"VSIR": "C10orf54"})
    return expr


def fetch_purity() -> pd.Series:
    data = _get(
        f"{API}/studies/{STUDY}/clinical-data",
        params={"clinicalDataType": "SAMPLE", "projection": "SUMMARY"},
    ).json()
    recs = [(d["sampleId"], float(d["value"])) for d in data
            if d.get("clinicalAttributeId") == "PURITY"
            and d.get("value") not in (None, "", "NA")]
    s = pd.Series({sid: v for sid, v in recs}, name="PURITY")
    return s


def main() -> None:
    expr = fetch_expression()
    purity = fetch_purity()
    print(f"[OncoSG] expression samples: {expr.shape[1]}, genes: {expr.shape[0]}")
    print(f"[OncoSG] samples with PURITY: {purity.notna().sum()}")

    if C.TARGET_GENE not in expr.index:
        raise SystemExit("TACSTD2 not found in OncoSG expression")

    # build features (samples x features)
    feats = {C.TARGET_GENE: expr.loc[C.TARGET_GENE].astype(float)}
    for name, genes in C.SIGNATURES.items():
        score, _ = C.signature_score(expr, genes)
        feats[name] = score
    avail = set(expr.index)
    for gname in C.SINGLE_IMMUNE_GENES:
        r = C.resolve_symbol(gname, avail)
        if r is not None:
            feats["gene_" + gname] = expr.loc[r].astype(float)
    for gname in C.TJ_GENES:
        r = C.resolve_symbol(gname, avail)
        if r is not None:
            feats["TJ_" + gname] = expr.loc[r].astype(float)
    feat = pd.DataFrame(feats)
    feat = feat.join(purity, how="left")

    x = feat[C.TARGET_GENE].to_numpy(float)
    p = feat["PURITY"].to_numpy(float)
    rows = []
    for col in [c for c in feat.columns if c not in (C.TARGET_GENE, "PURITY")]:
        y = feat[col].to_numpy(float)
        r, pv, n = C.spearman(x, y)
        pr, pp, pn = C.partial_spearman(x, y, p)
        rows.append({"cohort": "OncoSG_LUAD", "feature": col,
                     "spearman_rho": r, "spearman_p": pv, "n": n,
                     "partial_rho_PURITY": pr, "partial_p_PURITY": pp,
                     "partial_n_PURITY": pn})
    res = pd.DataFrame(rows).round(4)
    out = os.path.join(OUT, "oncosg_tacstd2_correlations.csv")
    res.to_csv(out, index=False)
    print(f"[OncoSG] wrote {out} ({len(res)} rows)")

    with open(os.path.join(OUT, "oncosg_sample_counts.txt"), "w") as fh:
        fh.write(f"expression_samples\t{expr.shape[1]}\n")
        fh.write(f"with_PURITY\t{int(purity.notna().sum())}\n")
        both = feat.dropna(subset=[C.TARGET_GENE, "PURITY"]).shape[0]
        fh.write(f"expr_and_purity\t{both}\n")

    sigs = list(C.SIGNATURES.keys())
    print("\n=== OncoSG LUAD: TACSTD2 vs signatures ===")
    print(res[res.feature.isin(sigs)][
        ["feature", "spearman_rho", "spearman_p",
         "partial_rho_PURITY", "partial_p_PURITY", "partial_n_PURITY"]
    ].to_string(index=False))


if __name__ == "__main__":
    main()
