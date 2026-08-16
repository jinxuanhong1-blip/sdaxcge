"""Harmonise TCGA-LUAD / TCGA-LUSC expression + covariates into a single cohort object.

Inputs (downloaded by 00_download.sh):
  UCSC Xena GDC hub  : TCGA-{LUAD,LUSC}.star_tpm.tsv.gz, clinical, survival, probemap
  GDC PanCanAtlas    : ABSOLUTE purity, leukocyte fraction, CIBERSORT, mutation load

Output: pickled dict of {expr, pheno} per cohort in $OPUS_TLS_DATA/processed/.
"""

from __future__ import annotations

import os
import pickle

import numpy as np
import pandas as pd

from opus_tls_lib import DATA_DIR, PROC_DIR

TCGA_DIR = os.path.join(DATA_DIR, "tcga")
META_DIR = os.path.join(DATA_DIR, "meta")


def load_probemap() -> pd.Series:
    pm = pd.read_csv(os.path.join(TCGA_DIR, "gencode.v36.annotation.gtf.gene.probemap"), sep="\t")
    return pm.set_index("id")["gene"]


def collapse_to_symbols(expr: pd.DataFrame, id2gene: pd.Series) -> pd.DataFrame:
    genes = expr.index.map(id2gene)
    expr = expr.loc[genes.notna()]
    expr.index = genes[genes.notna()]
    # For duplicated symbols keep the probe with the highest mean expression.
    order = expr.mean(axis=1).sort_values(ascending=False).index
    expr = expr.loc[order]
    expr = expr[~expr.index.duplicated(keep="first")]
    return expr.sort_index()


def short_barcode(x: pd.Series) -> pd.Series:
    return x.astype(str).str.replace(".", "-", regex=False).str.slice(0, 15)


def build_cohort(project: str, id2gene: pd.Series) -> dict:
    expr = pd.read_csv(os.path.join(TCGA_DIR, f"{project}.star_tpm.tsv.gz"), sep="\t", index_col=0)
    expr = collapse_to_symbols(expr, id2gene)

    sample_type = pd.Series(expr.columns, index=expr.columns).str.slice(13, 15)
    tumors = sample_type[sample_type.isin(["01"])].index  # primary solid tumour only
    normals = sample_type[sample_type == "11"].index

    pheno = pd.DataFrame(index=pd.Index(tumors, name="sample"))
    pheno["project"] = project
    pheno["patient"] = pd.Series(pheno.index, index=pheno.index).str.slice(0, 12)
    pheno["short"] = pd.Series(pheno.index, index=pheno.index).str.slice(0, 15)

    clin = pd.read_csv(os.path.join(TCGA_DIR, f"{project}.clinical.tsv.gz"), sep="\t", low_memory=False)
    clin = clin.drop_duplicates(subset="sample").set_index("sample")
    for src, dst in [
        ("gender.demographic", "sex"),
        ("age_at_index.demographic", "age"),
        ("ajcc_pathologic_stage.diagnoses", "stage_raw"),
        ("pack_years_smoked.exposures", "pack_years"),
        ("race.demographic", "race"),
    ]:
        if src in clin.columns:
            pheno[dst] = clin[src].reindex(pheno.index)

    stage_map = {"Stage I": 1, "Stage IA": 1, "Stage IB": 1, "Stage II": 2, "Stage IIA": 2,
                 "Stage IIB": 2, "Stage IIIA": 3, "Stage IIIB": 3, "Stage III": 3, "Stage IV": 4}
    pheno["stage"] = pheno.get("stage_raw", pd.Series(index=pheno.index)).map(stage_map)

    surv = pd.read_csv(os.path.join(TCGA_DIR, f"{project}.survival.tsv.gz"), sep="\t")
    surv = surv.drop_duplicates(subset="sample").set_index("sample")
    pheno["os_event"] = surv["OS"].reindex(pheno.index)
    pheno["os_time"] = surv["OS.time"].reindex(pheno.index)

    # ---- ABSOLUTE tumour purity (Aran et al. / PanCanAtlas mastercalls)
    pur = pd.read_csv(os.path.join(META_DIR, "absolute_purity.txt"), sep="\t")
    pur["short"] = short_barcode(pur["array"])
    pur = pur.drop_duplicates("short").set_index("short")
    pheno["purity_absolute"] = pheno["short"].map(pur["purity"])
    pheno["ploidy"] = pheno["short"].map(pur["ploidy"])

    # ---- Methylation-based leukocyte fraction (orthogonal to RNA)
    leuk = pd.read_csv(os.path.join(META_DIR, "leukocyte_fraction.tsv"), sep="\t", header=None,
                       names=["cancer", "barcode", "leuk_frac"])
    leuk["short"] = short_barcode(leuk["barcode"])
    leuk = leuk.drop_duplicates("short").set_index("short")
    pheno["leukocyte_fraction"] = pheno["short"].map(leuk["leuk_frac"])

    # ---- CIBERSORT relative fractions (orthogonal B/plasma deconvolution)
    cib = pd.read_csv(os.path.join(META_DIR, "cibersort_relative.tsv"), sep="\t")
    cib["short"] = short_barcode(cib["SampleID"])
    cib = cib.drop_duplicates("short").set_index("short")
    for col, dst in [("B.cells.naive", "cib_B_naive"), ("B.cells.memory", "cib_B_memory"),
                     ("Plasma.cells", "cib_plasma"), ("T.cells.follicular.helper", "cib_Tfh"),
                     ("T.cells.CD8", "cib_CD8")]:
        if col in cib.columns:
            pheno[dst] = pheno["short"].map(cib[col])
    if {"cib_B_naive", "cib_B_memory"}.issubset(pheno.columns):
        pheno["cib_B_total"] = pheno[["cib_B_naive", "cib_B_memory"]].sum(axis=1, min_count=1)

    # ---- TMB
    ml = pd.read_csv(os.path.join(META_DIR, "mutation_load.txt"), sep="\t")
    ml["short"] = ml["Tumor_Sample_ID"].astype(str).str.slice(0, 15)
    ml = ml.drop_duplicates("short").set_index("short")
    pheno["tmb_nonsilent"] = pheno["short"].map(ml["Non-silent per Mb"])

    expr_t = expr[tumors]
    expr_n = expr[normals] if len(normals) else None
    return {"expr": expr_t, "pheno": pheno, "expr_normal": expr_n}


def main() -> None:
    id2gene = load_probemap()
    out = {}
    for project in ["TCGA-LUAD", "TCGA-LUSC"]:
        print(f"[prep] {project} ...", flush=True)
        out[project] = build_cohort(project, id2gene)
        p = out[project]
        print(f"       genes={p['expr'].shape[0]} tumours={p['expr'].shape[1]} "
              f"normals={0 if p['expr_normal'] is None else p['expr_normal'].shape[1]}")
        print(f"       purity n={p['pheno']['purity_absolute'].notna().sum()}, "
              f"leukocyte n={p['pheno']['leukocyte_fraction'].notna().sum()}, "
              f"CIBERSORT n={p['pheno'].get('cib_plasma', pd.Series(dtype=float)).notna().sum()}")

    # Value sanity check: Xena GDC hub star_tpm is log2(TPM+1).
    e = out["TCGA-LUAD"]["expr"]
    print(f"[check] expression range {e.values.min():.3f} .. {e.values.max():.3f}; "
          f"median non-zero {np.median(e.values[e.values > 0]):.3f}")

    with open(os.path.join(PROC_DIR, "tcga.pkl"), "wb") as fh:
        pickle.dump(out, fh, protocol=4)
    print(f"[done] wrote {os.path.join(PROC_DIR, 'tcga.pkl')}")


if __name__ == "__main__":
    main()
