"""Extract the A11 gene panel from TCGA-LUAD / TCGA-LUSC and build a sample table.

Input : UCSC Xena GDC-hub STAR TPM matrices, log2(TPM+1), Ensembl-ID rows.
Output: data/derived/bulk_<cohort>_panel.parquet  (samples x genes)
        data/derived/bulk_<cohort>_samples.parquet (sample annotation + purity)
"""

import sys
import numpy as np
import pandas as pd

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import config as C

DERIVED = C.DATA / "derived"
DERIVED.mkdir(parents=True, exist_ok=True)


def load_probemap():
    pm = pd.read_csv(C.DATA / "bulk" / "gencode.v36.annotation.gtf.gene.probemap",
                     sep="\t")
    return dict(zip(pm["id"], pm["gene"]))


def load_purity():
    """ABSOLUTE purity/ploidy calls (GDC PanCanAtlas)."""
    ab = pd.read_csv(C.DATA / "bulk" / "absolute_purity.txt", sep="\t")
    ab = ab[["array", "purity", "ploidy", "Subclonal genome fraction"]].copy()
    ab.columns = ["sample15", "purity", "ploidy", "subclonal_frac"]
    ab = ab.dropna(subset=["purity"]).drop_duplicates("sample15")
    return ab.set_index("sample15")


SAMPLE_TYPE = {"01": "PrimaryTumor", "02": "RecurrentTumor",
               "06": "Metastatic", "11": "NormalAdjacent"}


def prep(cohort: str):
    print(f"[{cohort}] loading matrix ...", flush=True)
    df = pd.read_csv(C.DATA / "bulk" / f"TCGA-{cohort}.star_tpm.tsv.gz",
                     sep="\t", index_col=0)
    print(f"[{cohort}] matrix {df.shape[0]} genes x {df.shape[1]} samples", flush=True)

    sym = load_probemap()
    df.index = [sym.get(i, i) for i in df.index]

    wanted = set(C.all_genes())
    sub = df.loc[df.index.isin(wanted)]
    # Several symbols map to >1 Ensembl gene; keep the copy with the highest
    # mean expression, which is the expressed locus rather than a readthrough.
    sub = sub.assign(_m=sub.mean(axis=1)).sort_values("_m", ascending=False)
    sub = sub[~sub.index.duplicated(keep="first")].drop(columns="_m")

    missing = sorted(wanted - set(sub.index))
    if missing:
        print(f"[{cohort}] WARNING genes absent from annotation: {missing}", flush=True)

    mat = sub.T  # samples x genes
    mat.index.name = "sample"

    samp = pd.DataFrame(index=mat.index)
    samp["patient"] = [s[:12] for s in samp.index]
    samp["type_code"] = [s[13:15] for s in samp.index]
    samp["sample_type"] = samp["type_code"].map(SAMPLE_TYPE).fillna("Other")
    samp["cohort"] = cohort
    samp["sample15"] = [s[:15] for s in samp.index]

    pur = load_purity()
    samp = samp.join(pur, on="sample15")

    n_pur = samp["purity"].notna().sum()
    print(f"[{cohort}] {samp['sample_type'].value_counts().to_dict()}", flush=True)
    print(f"[{cohort}] ABSOLUTE purity available for {n_pur}/{len(samp)} samples", flush=True)

    mat.to_parquet(DERIVED / f"bulk_{cohort}_panel.parquet")
    samp.to_parquet(DERIVED / f"bulk_{cohort}_samples.parquet")
    print(f"[{cohort}] wrote panel {mat.shape}", flush=True)


if __name__ == "__main__":
    for coh in ["LUAD", "LUSC"]:
        prep(coh)
