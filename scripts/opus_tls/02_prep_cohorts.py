"""Harmonise the open GEO lung cohorts (ICI-treated and treatment-naive atlases).

Every cohort ends up as {expr: genes x samples on a log2 scale, pheno: samples x vars}.
Ensembl IDs are mapped with the same GENCODE v36 probemap used for TCGA so that
gene identity is consistent across cohorts.
"""

from __future__ import annotations

import gzip
import io
import os
import pickle
import re

import numpy as np
import pandas as pd

from opus_tls_lib import DATA_DIR, PROC_DIR

GEO_DIR = os.path.join(DATA_DIR, "geo")
TCGA_DIR = os.path.join(DATA_DIR, "tcga")


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def ens2symbol() -> pd.Series:
    pm = pd.read_csv(os.path.join(TCGA_DIR, "gencode.v36.annotation.gtf.gene.probemap"), sep="\t")
    pm["ens"] = pm["id"].str.split(".").str[0]
    return pm.drop_duplicates("ens").set_index("ens")["gene"]


def collapse(expr: pd.DataFrame) -> pd.DataFrame:
    expr = expr[~expr.index.isna()]
    order = expr.mean(axis=1).sort_values(ascending=False).index
    expr = expr.loc[order]
    return expr[~expr.index.duplicated(keep="first")].sort_index()


def parse_series_matrix_meta(path: str) -> pd.DataFrame:
    """Return samples x characteristics from a GEO series matrix file."""
    rows, titles, accs = [], None, None
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
            if line.startswith("!Sample_title"):
                titles = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_geo_accession"):
                accs = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                rows.append([x.strip('"') for x in line.rstrip("\n").split("\t")[1:]])
    meta = pd.DataFrame(index=pd.Index(accs, name="gsm"))
    meta["title"] = titles
    for row in rows:
        keys = [v.split(":", 1)[0].strip().lower().replace(" ", "_") for v in row if ":" in v]
        if not keys:
            continue
        key = pd.Series(keys).mode()[0]
        vals = [v.split(":", 1)[1].strip() if ":" in v else np.nan for v in row]
        col = key
        i = 2
        while col in meta.columns:
            col = f"{key}_{i}"
            i += 1
        meta[col] = vals
    return meta


def read_series_matrix_table(path: str) -> pd.DataFrame:
    with gzip.open(path, "rt", errors="replace") as fh:
        lines = fh.read().split("\n")
    start = next(i for i, l in enumerate(lines) if l.startswith("!series_matrix_table_begin")) + 1
    end = next(i for i, l in enumerate(lines) if l.startswith("!series_matrix_table_end"))
    tbl = pd.read_csv(io.StringIO("\n".join(lines[start:end])), sep="\t", index_col=0)
    tbl.index = tbl.index.astype(str).str.strip('"')
    return tbl


def log2p1(df: pd.DataFrame) -> pd.DataFrame:
    return np.log2(df.astype(float) + 1)


# ---------------------------------------------------------------------------
# cohorts
# ---------------------------------------------------------------------------

def cohort_gse135222(e2s: pd.Series) -> dict:
    """NSCLC, anti-PD-1/PD-L1 monotherapy, pre-treatment RNA-seq, PFS. n=27."""
    expr = pd.read_csv(os.path.join(GEO_DIR, "GSE135222_exp.tsv.gz"), sep="\t", index_col=0)
    expr.index = expr.index.str.split(".").str[0].map(e2s)
    expr = collapse(log2p1(expr))
    meta = parse_series_matrix_meta(os.path.join(GEO_DIR, "GSE135222_matrix.txt.gz"))
    meta["sid"] = meta["title"].str.replace(" ", "", regex=False)
    meta = meta.set_index("sid")
    pheno = pd.DataFrame(index=expr.columns)
    pheno["sex"] = meta["gender"].reindex(pheno.index)
    pheno["age"] = pd.to_numeric(meta["age"].reindex(pheno.index), errors="coerce")
    pfs_col = [c for c in meta.columns if c.startswith("progression-free")][0]
    pheno["pfs_event"] = pd.to_numeric(meta[pfs_col].reindex(pheno.index), errors="coerce")
    pheno["pfs_time"] = pd.to_numeric(meta["pfs.time"].reindex(pheno.index), errors="coerce")
    # Durable clinical benefit proxy used in the original report: PFS >= 6 months.
    pheno["dcb_6mo"] = np.where(pheno["pfs_time"] >= 180, 1,
                                np.where(pheno["pfs_event"] == 1, 0, np.nan))
    return {"expr": expr, "pheno": pheno, "kind": "ici", "platform": "RNA-seq (TPM)",
            "note": "advanced NSCLC, anti-PD-1/PD-L1, pre-treatment biopsy"}


def cohort_gse126044(e2s: pd.Series) -> dict:
    """NSCLC, anti-PD-1, pre-treatment RNA-seq counts, RECIST responder labels. n=16."""
    cnt = pd.read_csv(os.path.join(GEO_DIR, "GSE126044_counts.txt.gz"), sep="\t", index_col=0)
    cpm = cnt.div(cnt.sum(axis=0), axis=1) * 1e6
    expr = collapse(log2p1(cpm))
    meta = parse_series_matrix_meta(os.path.join(GEO_DIR, "GSE126044_matrix.txt.gz"))
    meta["sid"] = meta["title"].str.replace("RNA-seq_", "", regex=False)
    meta = meta.set_index("sid")
    pheno = pd.DataFrame(index=expr.columns)
    resp = meta["patient_response"].reindex(pheno.index)
    pheno["response"] = resp
    pheno["responder"] = (resp == "responder").astype(float).where(resp.notna())
    if "sample" in meta.columns:
        pheno["preservation"] = meta["sample"].reindex(pheno.index)
    return {"expr": expr, "pheno": pheno, "kind": "ici", "platform": "RNA-seq (counts->logCPM)",
            "note": "NSCLC, anti-PD-1, pre-treatment biopsy"}


def cohort_gse207422() -> dict:
    """NSCLC, neoadjuvant anti-PD-1 + chemotherapy, bulk RNA-seq (log2 TPM), MPR labels."""
    expr = pd.read_csv(os.path.join(GEO_DIR, "GSE207422_bulk_log2TPM.txt.gz"), sep="\t", index_col=0)
    expr = collapse(expr)
    meta = pd.read_excel(os.path.join(GEO_DIR, "GSE207422_bulk_meta.xlsx")).set_index("Sample")
    pheno = pd.DataFrame(index=expr.columns)
    pheno["patient"] = meta["Patient"].reindex(pheno.index)
    pheno["timepoint"] = np.where(
        meta["Resource"].reindex(pheno.index).astype(str).str.lower().str.contains("pre"), "pre", "post")
    pheno["histology"] = meta["Pathology"].reindex(pheno.index)
    pheno["sex"] = meta["Sex"].reindex(pheno.index)
    pheno["age"] = pd.to_numeric(meta["Age"].reindex(pheno.index), errors="coerce")
    pr = meta["Pathologic Response"].reindex(pheno.index).astype(str)
    pheno["mpr"] = np.where(pr.str.startswith("MPR"), 1, np.where(pr == "NMPR", 0, np.nan))
    pheno["residual_tumor"] = pd.to_numeric(meta["Residual Tumor"].reindex(pheno.index), errors="coerce")
    pheno["recist"] = meta["RECIST"].reindex(pheno.index)
    return {"expr": expr, "pheno": pheno, "kind": "ici", "platform": "RNA-seq (log2 TPM)",
            "note": "resectable NSCLC, neoadjuvant anti-PD-1 + chemotherapy"}


def cohort_gse81089(e2s: pd.Series) -> dict:
    """Uppsala NSCLC surgical cohort, RNA-seq FPKM, OS. Treatment-naive atlas."""
    expr = pd.read_csv(os.path.join(GEO_DIR, "GSE81089_fpkm.tsv.gz"), sep="\t", index_col=0)
    expr.index = expr.index.astype(str).str.split(".").str[0].map(e2s)
    expr = collapse(log2p1(expr))
    tumor_cols = [c for c in expr.columns if c.endswith("T")]
    expr = expr[tumor_cols]
    meta = parse_series_matrix_meta(os.path.join(GEO_DIR, "GSE81089_matrix.txt.gz")).set_index("title")
    pheno = pd.DataFrame(index=expr.columns)
    pheno["sex"] = meta["gender"].reindex(pheno.index)
    pheno["age"] = pd.to_numeric(meta["age"].reindex(pheno.index), errors="coerce")
    stage = pd.to_numeric(meta["stage_tnm"].reindex(pheno.index), errors="coerce")
    pheno["stage"] = stage.replace({5: 3})  # coding in GEO: 1=IA,2=IB,3=IIA/B,4=IIIA,5=IIIB
    hist = pd.to_numeric(meta["histology"].reindex(pheno.index), errors="coerce")
    pheno["histology"] = hist.map({1: "AC", 2: "SqCC", 3: "LCC"})
    pheno["os_event"] = pd.to_numeric(meta["dead"].reindex(pheno.index), errors="coerce")
    d1 = pd.to_datetime(meta["surgery_date"].reindex(pheno.index), errors="coerce")
    d2 = pd.to_datetime(meta["vital_date"].reindex(pheno.index), errors="coerce")
    pheno["os_time"] = (d2 - d1).dt.days
    pheno["smoking"] = pd.to_numeric(meta["smoking"].reindex(pheno.index), errors="coerce")
    return {"expr": expr, "pheno": pheno, "kind": "atlas", "platform": "RNA-seq (FPKM)",
            "note": "Uppsala NSCLC, surgically resected, treatment-naive"}


def cohort_gse72094() -> dict:
    """Moffitt/Schabath LUAD array cohort (n>400), OS and driver mutation status."""
    tbl = read_series_matrix_table(os.path.join(GEO_DIR, "GSE72094_matrix.txt.gz"))
    ann = pd.read_csv(os.path.join(GEO_DIR, "GPL15048.txt"), sep="\t", skiprows=7, low_memory=False)
    ann = ann.dropna(subset=["GeneSymbol"]).drop_duplicates("ID").set_index("ID")
    tbl.index = tbl.index.map(ann["GeneSymbol"])
    expr = collapse(tbl)  # already log2 scale
    meta = parse_series_matrix_meta(os.path.join(GEO_DIR, "GSE72094_matrix.txt.gz"))
    meta = meta.reindex(expr.columns)
    pheno = pd.DataFrame(index=expr.columns)
    pheno["sex"] = meta["gender"]
    pheno["age"] = pd.to_numeric(meta["age_at_diagnosis"], errors="coerce")
    pheno["smoking"] = meta["smoking_status"]
    pheno["os_event"] = meta["vital_status"].map({"Dead": 1, "Alive": 0})
    pheno["os_time"] = pd.to_numeric(meta["survival_time_in_days"], errors="coerce")
    stage = meta["stage"].astype(str).str.extract(r"^(\d)")[0]
    pheno["stage"] = pd.to_numeric(stage, errors="coerce")
    for g in ["kras_status", "egfr_status", "stk11_status", "tp53_status"]:
        if g in meta.columns:
            pheno[g] = meta[g]
    return {"expr": expr, "pheno": pheno, "kind": "atlas", "platform": "Affymetrix HuRSTA (log2)",
            "note": "LUAD, surgically resected, treatment-naive"}


def cohort_gse190265() -> dict:
    """France3 NSCLC biopsies (TPM). Small independent bulk cohort."""
    df = pd.read_csv(os.path.join(GEO_DIR, "GSE190265_tpm.csv.gz"), sep=";", index_col=0, low_memory=False)
    expr = collapse(log2p1(df.T))
    meta = parse_series_matrix_meta(os.path.join(GEO_DIR, "GSE190265_matrix.txt.gz")).set_index("title")
    expr.columns = [str(c) for c in expr.columns]
    meta.index = [str(i) for i in meta.index]
    keep = [c for c in expr.columns if c in meta.index]
    expr = expr[keep]
    pheno = pd.DataFrame(index=expr.columns)
    pheno["histology"] = meta["disease_state"].reindex(pheno.index)
    return {"expr": expr, "pheno": pheno, "kind": "atlas", "platform": "RNA-seq (TPM)",
            "note": "NSCLC biopsies (France3 series of GSE190265)"}


def main() -> None:
    e2s = ens2symbol()
    builders = {
        "GSE135222": lambda: cohort_gse135222(e2s),
        "GSE126044": lambda: cohort_gse126044(e2s),
        "GSE207422": cohort_gse207422,
        "GSE81089": lambda: cohort_gse81089(e2s),
        "GSE72094": cohort_gse72094,
        "GSE190265": cohort_gse190265,
    }
    out = {}
    for name, fn in builders.items():
        print(f"[prep] {name} ...", flush=True)
        c = fn()
        out[name] = c
        miss = [g for g in ["TACSTD2", "CLDN4", "CXCL13", "MS4A1", "MZB1"] if g not in c["expr"].index]
        print(f"       genes={c['expr'].shape[0]} samples={c['expr'].shape[1]} "
              f"kind={c['kind']} missing_key_genes={miss}")
    with open(os.path.join(PROC_DIR, "geo_cohorts.pkl"), "wb") as fh:
        pickle.dump(out, fh, protocol=4)
    print(f"[done] wrote {os.path.join(PROC_DIR, 'geo_cohorts.pkl')}")


if __name__ == "__main__":
    main()
