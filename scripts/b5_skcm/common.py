"""Shared loaders and constants for the B5_SKCM analysis (CLDN4 vs ICI response in melanoma)."""

from __future__ import annotations

import gzip
import io
import re
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
RAW = REPO / "data" / "raw"
PROC = REPO / "data" / "processed"
OUT = REPO / "results" / "w200" / "B5_SKCM"

# log2 offset, fixed in the analysis plan
OFFSET = 0.1

GENE_OF_INTEREST = "CLDN4"

POSITIVE_CONTROLS = [
    "CD8A", "CD274", "IFNG", "GZMB", "PRF1", "CXCL9", "CXCL10", "STAT1", "HLA-DRA",
]
AYERS_IFNG_6 = ["IFNG", "STAT1", "CXCL9", "CXCL10", "IDO1", "HLA-DRA"]
EPITHELIAL_CONTROLS = [
    "CLDN3", "CLDN7", "EPCAM", "CDH1", "KRT5", "KRT14", "KRT1", "KRT10", "SFN",
]
# genes used to build the keratinocyte / epidermal-contamination score
KERATINOCYTE_SCORE_GENES = ["KRT1", "KRT5", "KRT10", "KRT14", "SFN", "DSP", "KRT6A"]
MELANOCYTE_GENES = ["MLANA", "PMEL", "TYR", "DCT", "SOX10", "MITF"]


def log2t(x: pd.DataFrame | pd.Series) -> pd.DataFrame | pd.Series:
    return np.log2(x + OFFSET)


def collapse_duplicate_symbols(df: pd.DataFrame) -> pd.DataFrame:
    """One row per gene symbol, keeping the row with the highest mean expression.

    Applied identically to every cohort. Recorded rather than silent: the caller logs how many
    symbols were collapsed.
    """
    if not df.index.has_duplicates:
        return df
    order = df.mean(axis=1).sort_values(ascending=False)
    df = df.loc[order.index]
    return df[~df.index.duplicated(keep="first")]


def entrez_to_symbol() -> pd.Series:
    """Map NCBI Entrez GeneID -> current HGNC-style symbol from Homo_sapiens.gene_info."""
    path = RAW / "Homo_sapiens.gene_info.gz"
    gi = pd.read_csv(path, sep="\t", usecols=["GeneID", "Symbol", "type_of_gene"], dtype={"GeneID": str})
    return gi.set_index("GeneID")["Symbol"]


def _parse_series_matrix(path: Path) -> pd.DataFrame:
    """Parse !Sample_* lines of a GEO series matrix into a sample x field table."""
    with gzip.open(path, "rt", errors="replace") as fh:
        lines = [ln.rstrip("\n") for ln in fh if ln.startswith("!Sample_")]
    rows: dict[str, list[str]] = {}
    char_rows: list[list[str]] = []
    for ln in lines:
        key, _, rest = ln.partition("\t")
        vals = [v.strip('"') for v in rest.split("\t")]
        key = key[1:]
        if key == "Sample_characteristics_ch1":
            char_rows.append(vals)
        elif key in ("Sample_title", "Sample_geo_accession"):
            rows[key] = vals
    meta = pd.DataFrame(rows)
    # characteristics are position-shuffled across samples in some series: parse "key: value"
    extra: dict[str, list[str | float]] = {}
    n = len(meta)
    for vals in char_rows:
        for i, v in enumerate(vals):
            if ":" not in v:
                continue
            k, _, val = v.partition(":")
            k = k.strip().lower()
            extra.setdefault(k, [np.nan] * n)[i] = val.strip()
    for k, v in extra.items():
        meta[k] = v
    return meta


# --------------------------------------------------------------------------------------------
# Cohort loaders. Each returns (log2 expression [genes x samples], clinical [samples x fields]).
# --------------------------------------------------------------------------------------------


def load_hugo() -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Hugo 2016, GSE78220, pembrolizumab. Author-provided FPKM keyed by gene symbol."""
    fpkm = pd.read_excel(RAW / "GSE78220_PatientFPKM.xlsx").set_index("Gene")
    n_dup = int(fpkm.index.duplicated().sum())
    fpkm = collapse_duplicate_symbols(fpkm)
    fpkm.columns = [c.replace(".baseline", "").replace(".OnTx", "") for c in fpkm.columns]

    meta = _parse_series_matrix(RAW / "GSE78220_series_matrix.txt.gz")
    meta = meta.rename(
        columns={
            "anti-pd-1 response": "response_raw",
            "biopsy time": "timepoint_raw",
            "patient id": "patient",
            "overall survival (days)": "os_days",
            "vital status": "vital_status",
            "study site": "site",
            "anatomical location": "biopsy_site",
            "stranded/unstranded rnaseq": "library",
            "age (yrs)": "age",
            "gender": "sex",
            "treatment": "treatment",
        }
    ).set_index("Sample_title")

    clin = pd.DataFrame(index=fpkm.columns)
    clin.index.name = "sample"
    for c in [
        "Sample_geo_accession", "patient", "response_raw", "timepoint_raw", "os_days",
        "vital_status", "site", "biopsy_site", "library", "age", "sex", "treatment",
    ]:
        clin[c] = meta[c].reindex(clin.index)
    clin["cohort"] = "Hugo_GSE78220"
    clin["pretreatment"] = clin["timepoint_raw"].eq("pre-treatment")
    clin["response"] = clin["response_raw"].map(
        {"Complete Response": "R", "Partial Response": "R", "Progressive Disease": "NR"}
    )
    clin["os_days"] = pd.to_numeric(clin["os_days"], errors="coerce")
    clin["os_event"] = clin["vital_status"].map({"Dead": 1, "Alive": 0})
    clin["age"] = pd.to_numeric(clin["age"], errors="coerce")
    clin["therapy"] = "anti-PD-1"
    clin["batch"] = clin["library"]
    clin["episode_order"] = 1
    return log2t(fpkm), clin, {"unit": "FPKM", "duplicate_symbols_collapsed": n_dup}


def load_riaz() -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Riaz 2017, GSE91061, nivolumab. Author FPKM on hg19 knownGene, rows are Entrez GeneIDs."""
    fpkm = pd.read_csv(RAW / "GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv.gz", index_col=0)
    fpkm.index = fpkm.index.astype(str)
    sym = entrez_to_symbol()
    mapped = fpkm.index.map(sym)
    n_unmapped = int(pd.isna(mapped).sum())
    fpkm = fpkm[~pd.isna(mapped)]
    fpkm.index = pd.Index([s for s in mapped if not pd.isna(s)], name="Gene")
    n_dup = int(fpkm.index.duplicated().sum())
    fpkm = collapse_duplicate_symbols(fpkm)

    meta = _parse_series_matrix(RAW / "GSE91061_series_matrix.txt.gz")
    meta = meta.rename(
        columns={"visit (pre or on treatment)": "timepoint_raw", "response": "response_raw"}
    ).set_index("Sample_title")

    clin = pd.DataFrame(index=fpkm.columns)
    clin.index.name = "sample"
    for c in ["Sample_geo_accession", "timepoint_raw", "response_raw"]:
        clin[c] = meta[c].reindex(clin.index)
    clin["patient"] = [s.split("_")[0] for s in clin.index]
    clin["cohort"] = "Riaz_GSE91061"
    clin["pretreatment"] = clin["timepoint_raw"].eq("Pre")
    clin["response"] = clin["response_raw"].map({"PRCR": "R", "PD": "NR"})  # SD/UNK -> NaN
    clin["response_sd_as_nr"] = clin["response_raw"].map({"PRCR": "R", "PD": "NR", "SD": "NR"})
    clin["therapy"] = "anti-PD-1"
    clin["batch"] = np.nan
    clin["episode_order"] = 1
    return log2t(fpkm), clin, {
        "unit": "FPKM",
        "entrez_ids_unmapped": n_unmapped,
        "duplicate_symbols_collapsed": n_dup,
    }


def load_mgh() -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Auslander 2018, GSE115821. Author raw counts + feature length -> TPM."""
    cnt = pd.read_csv(RAW / "GSE115821_MGH_counts.csv.gz")
    cnt = cnt.rename(columns={cnt.columns[0]: "Geneid"})
    length = pd.to_numeric(cnt["Length"], errors="coerce")
    mat = cnt.drop(columns=["Geneid", "Chr", "Start", "End", "Strand", "Length"]).apply(
        pd.to_numeric, errors="coerce"
    )
    mat.index = cnt["Geneid"].values
    length.index = cnt["Geneid"].values
    rate = mat.div(length / 1000.0, axis=0)
    tpm = rate.div(rate.sum(axis=0), axis=1) * 1e6
    tpm.columns = [c[:-4] if c.endswith(".bam") else c for c in tpm.columns]
    n_dup = int(tpm.index.duplicated().sum())
    tpm = collapse_duplicate_symbols(tpm)

    soft = gzip.open(RAW / "GSE115821_family.soft.gz", "rt", errors="replace").read()
    recs = []
    for block in soft.split("^SAMPLE = ")[1:]:
        gsm = block.split("\n", 1)[0].strip()
        title = re.search(r"!Sample_title = (.+)", block).group(1).strip()
        d = {"Sample_geo_accession": gsm, "title": title}
        for k, v in re.findall(r"!Sample_characteristics_ch1 = ([^:]+): (.+)", block):
            d[k.strip().lower()] = v.strip()
        recs.append(d)
    meta = pd.DataFrame(recs)
    meta["sample"] = [t[:-4] if t.endswith(".bam") else t for t in meta["title"]]
    meta = meta.set_index("sample")

    clin = pd.DataFrame(index=tpm.columns)
    clin.index.name = "sample"
    clin["Sample_geo_accession"] = meta["Sample_geo_accession"].reindex(clin.index)
    clin["patient"] = meta["patient id"].reindex(clin.index)
    clin["response_raw"] = meta["response"].reindex(clin.index)
    clin["timepoint_raw"] = meta["treatment state"].reindex(clin.index)
    clin["therapy"] = meta["antibody"].reindex(clin.index)
    clin["batch"] = meta["batch"].reindex(clin.index)
    clin["age"] = pd.to_numeric(meta["age at the baseline"].reindex(clin.index), errors="coerce")
    clin["timepoint_day"] = pd.to_numeric(meta["timepoint"].reindex(clin.index), errors="coerce")
    clin["cohort"] = "MGH_GSE115821"
    clin["pretreatment"] = clin["timepoint_raw"].str.startswith("PRE").fillna(False)
    clin["response"] = clin["response_raw"].map({"R": "R", "NR": "NR"})
    # a patient may have sequential treatment episodes (anti-CTLA-4 then anti-PD-1); rank the
    # baseline biopsies of each patient by biopsy day so the plan's "first episode" rule is explicit
    pre = clin[clin["pretreatment"]].copy()
    pre["episode_order"] = (
        pre.sort_values(["patient", "timepoint_day"]).groupby("patient").cumcount() + 1
    ).reindex(pre.index)
    clin["episode_order"] = pre["episode_order"].reindex(clin.index)
    return log2t(tpm), clin, {"unit": "TPM (recomputed from counts and feature length)",
                              "duplicate_symbols_collapsed": n_dup}


def load_tcga_skcm() -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """TCGA-SKCM from the UCSC Xena GDC hub (log2(TPM+1)). Non-ICI reference cohort."""
    tpm = pd.read_csv(RAW / "TCGA-SKCM.star_tpm.tsv.gz", sep="\t", index_col=0)
    probemap = pd.read_csv(RAW / "gencode.v36.annotation.gtf.gene.probemap", sep="\t", index_col=0)
    tpm.index = tpm.index.map(probemap["gene"])
    tpm = tpm[~pd.isna(tpm.index)]
    tpm = collapse_duplicate_symbols(tpm)
    # Xena distributes log2(tpm+1); convert to the analysis transform log2(tpm+OFFSET)
    lin = np.power(2.0, tpm) - 1.0
    lin[lin < 0] = 0.0
    clin = pd.DataFrame(index=tpm.columns)
    clin.index.name = "sample"
    clin["sample_type_code"] = [s.split("-")[3][:2] if len(s.split("-")) > 3 else None
                                for s in clin.index]
    clin["is_tumor"] = clin["sample_type_code"].isin(["01", "06"])
    clin["cohort"] = "TCGA_SKCM"
    return log2t(lin), clin, {"unit": "TPM (from Xena log2(TPM+1))"}
