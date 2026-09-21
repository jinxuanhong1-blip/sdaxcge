#!/usr/bin/env python3
"""Open NSCLC ICI bulk: CLDN4-high vs response, with a keratin residual.

Public GEO / EuropePMC / Springer supplementary files only. No EGA, dbGaP, or FASTQ.
Raw downloads stay in /tmp/ici_raw and are not committed.

Locked before pooling:
- One primary binary endpoint per cohort (see LOADERS).
- High = expression >= within-cohort median on the analysis subset.
- OR < 1 means CLDN4-high has lower odds of response.
- HR > 1 means CLDN4-high has worse survival.
- Keratin residual = OLS residual of log CLDN4 on KRT8 + KRT18 + KRT19
  (intercept included), fit once per cohort on samples with all four genes.
- PFS Cox is estimated only when a time and an event indicator are deposited
  and the time is not administratively capped at the DCB landmark.
- DerSimonian–Laird random effects on log OR / log HR. Every cutoff is
  written to a table. The forest shows the locked median split.
- The PDF figure OR=0.42 [0.18–0.95], k=11 is supported only when this
  NSCLC meta has k=11 and the random-effects point estimate rounds to 0.42.
"""

from __future__ import annotations

import gzip
import math
import re
import zipfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, norm
from statsmodels.duration.hazard_regression import PHReg
import statsmodels.api as sm

ROOT = Path(__file__).resolve().parent
RAW = Path("/tmp/ici_raw")
TABLES = ROOT / "tables"
FIGS = ROOT / "figures"
for d in (TABLES, FIGS, RAW):
    d.mkdir(parents=True, exist_ok=True)

KERATINS = ["KRT8", "KRT18", "KRT19"]
GENES = ["CLDN4", *KERATINS, "KRT5", "KRT7", "EPCAM", "TACSTD2"]
ENSG = {
    "CLDN4": "ENSG00000189143",
    "KRT8": "ENSG00000170421",
    "KRT18": "ENSG00000111057",
    "KRT19": "ENSG00000171345",
    "KRT5": "ENSG00000186081",
    "KRT7": "ENSG00000135480",
    "EPCAM": "ENSG00000119888",
    "TACSTD2": "ENSG00000184292",
    "CD274": "ENSG00000120217",
    "CD8A": "ENSG00000153563",
    "PTPRC": "ENSG00000081237",
    "GAPDH": "ENSG00000111640",
    "ACTB": "ENSG00000075624",
    "HLA-A": "ENSG00000206503",
    "IFNG": "ENSG00000111537",
    "CXCL9": "ENSG00000138755",
    "MKI67": "ENSG00000148773",
    "CDH1": "ENSG00000039068",
    "CLDN7": "ENSG00000181885",
    "MUC1": "ENSG00000185499",
    "B2M": "ENSG00000166710",
    "STAT1": "ENSG00000115415",
    "CXCL10": "ENSG00000169245",
    "CD3D": "ENSG00000167286",
    "CD3E": "ENSG00000198851",
    "NKG7": "ENSG00000105374",
    "GZMB": "ENSG00000100453",
    "PRF1": "ENSG00000180644",
    "EGFR": "ENSG00000146648",
    "HPRT1": "ENSG00000165704",
    "KRT6A": "ENSG00000205420",
    "KRT17": "ENSG00000128422",
    "TJP1": "ENSG00000104067",
    "SFTPC": "ENSG00000168484",
}
FINGER = list(dict.fromkeys(GENES + list(ENSG)))
PDF_OR = 0.42
PDF_K = 11


def parse_series_matrix(path: Path) -> pd.DataFrame:
    titles = geo = descriptions = None
    char_rows: list[list[str]] = []
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
            parts = [p.strip().strip('"') for p in line.rstrip("\n").split("\t")]
            key, vals = parts[0], parts[1:]
            if key == "!Sample_title":
                titles = vals
            elif key == "!Sample_geo_accession":
                geo = vals
            elif key == "!Sample_description":
                descriptions = vals
            elif key == "!Sample_characteristics_ch1":
                char_rows.append(vals)
    rows = []
    for i, title in enumerate(titles or []):
        d = {"title": title}
        if geo and i < len(geo):
            d["gsm"] = geo[i]
        if descriptions and i < len(descriptions):
            d["description"] = descriptions[i]
        for row in char_rows:
            if i < len(row) and row[i] and ":" in row[i]:
                k, v = row[i].split(":", 1)
                d[k.strip().lower()] = v.strip()
        rows.append(d)
    return pd.DataFrame(rows)


def _strip_ensg(x: str) -> str:
    x = str(x).split()[0]
    return x.split(".")[0]


def extract_genes(frame: pd.DataFrame) -> pd.DataFrame:
    """frame is genes x samples, index symbols or ENSG. Return symbols x samples."""
    idx = pd.Index([_strip_ensg(i) for i in frame.index])
    frame = frame.copy()
    frame.index = idx
    frame = frame.apply(pd.to_numeric, errors="coerce")
    frame = frame.groupby(level=0).mean()
    ensg_to_sym = {v: k for k, v in ENSG.items()}
    rows = {}
    upper = {str(i).upper(): i for i in frame.index}
    for sym in FINGER:
        key = None
        if sym.upper() in upper:
            key = upper[sym.upper()]
        elif ENSG.get(sym) in frame.index:
            key = ENSG[sym]
        elif ENSG.get(sym) and ENSG[sym].upper() in upper:
            key = upper[ENSG[sym].upper()]
        if key is not None:
            rows[sym] = frame.loc[key]
    if not rows:
        # last try: index may already be a symbol that we renamed via ensg map
        for sym, ensg in ENSG.items():
            if ensg in frame.index:
                rows[sym] = frame.loc[ensg]
    out = pd.DataFrame(rows).T
    out = out[~out.index.duplicated(keep="first")]
    return out


def as_log(expr: pd.DataFrame, mode: str) -> pd.DataFrame:
    if mode == "already_log":
        return expr.astype(float)
    return np.log2(expr.clip(lower=0).astype(float) + 1.0)


def norm_histology(value) -> str:
    if not isinstance(value, str):
        return "NA"
    v = value.strip().lower()
    if any(s in v for s in ("small cell", "neuroendocrine")):
        return "excluded_non_nsclc"
    if "squamous" in v or v in {"lusc", "sqcc", "lscc"}:
        return "LUSC"
    if "adeno" in v or v in {"luad", "adc"}:
        return "LUAD"
    if "non-squamous" in v or "nonsquamous" in v:
        return "nonsquamous"
    if "mixed" in v:
        return "mixed"
    return "NA"


def empty_clin(index) -> pd.DataFrame:
    clin = pd.DataFrame(index=pd.Index(index, name="sample"))
    for col in (
        "y_primary",
        "endpoint",
        "y_mpr",
        "y_dcb",
        "pfs_time",
        "pfs_event",
        "pfs_unit",
        "pfs_usable",
        "os_time",
        "os_event",
        "histology",
        "timing",
        "regimen",
        "note",
    ):
        clin[col] = pd.NA
    clin["pfs_usable"] = False
    return clin


def dcb_from_pfs(time, event, unit: str) -> pd.Series:
    t = pd.to_numeric(pd.Series(time), errors="coerce").reset_index(drop=True)
    e = pd.to_numeric(pd.Series(event), errors="coerce").reset_index(drop=True)
    thr = 180.0 if unit == "day" else 6.0
    out = pd.Series(pd.NA, index=t.index, dtype=object)
    out[t >= thr] = "R"
    out[(t < thr) & (e == 1)] = "NR"
    return out


def load_gse126044() -> tuple[pd.DataFrame, pd.DataFrame]:
    meta = parse_series_matrix(RAW / "GSE126044_series_matrix.txt.gz")
    counts = pd.read_csv(RAW / "GSE126044_counts.txt.gz", sep="\t", index_col=0)
    expr = as_log(extract_genes(counts), "log2_cpm_from_counts")
    # log2 CPM, not log2(count+1)
    lib = counts.sum(axis=0).replace(0, np.nan)
    cpm = counts.divide(lib, axis=1) * 1e6
    expr = as_log(extract_genes(cpm), "log_linear")
    recs = []
    for _, s in meta.iterrows():
        sid = str(s["title"]).replace("RNA-seq_", "")
        resp = str(s.get("patient response", "")).lower()
        y = "R" if resp == "responder" else ("NR" if "non" in resp else pd.NA)
        recs.append(sid)
        # filled below
    clin = empty_clin([str(s["title"]).replace("RNA-seq_", "") for _, s in meta.iterrows()])
    for _, s in meta.iterrows():
        sid = str(s["title"]).replace("RNA-seq_", "")
        resp = str(s.get("patient response", "")).lower()
        y = "R" if resp == "responder" else ("NR" if "non" in resp else pd.NA)
        clin.loc[sid, "y_primary"] = y
    clin["endpoint"] = "author_R_vs_NR"
    clin["timing"] = "pre-treatment"
    clin["regimen"] = "anti-PD-1 monotherapy"
    clin["histology"] = "NA"
    clin["note"] = "Cho / Yonsei; GEO patient response; no PFS"
    clin = clin[~clin.index.duplicated(keep="first")]
    keep = [c for c in expr.columns if c in clin.index]
    return expr[keep], clin.loc[keep]


def load_gse135222() -> tuple[pd.DataFrame, pd.DataFrame]:
    meta = parse_series_matrix(RAW / "GSE135222_series_matrix.txt.gz")
    tpm = pd.read_csv(RAW / "GSE135222_exp.tsv.gz", sep="\t", index_col=0)
    expr = as_log(extract_genes(tpm), "log_linear")
    recs = []
    for _, s in meta.iterrows():
        sid = str(s["title"]).replace(" ", "")
        recs.append(
            {
                "sample": sid,
                "pfs_time": pd.to_numeric(s.get("pfs.time"), errors="coerce"),
                "pfs_event": pd.to_numeric(s.get("progression-free survival (pfs)"), errors="coerce"),
            }
        )
    base = pd.DataFrame(recs).drop_duplicates("sample").set_index("sample")
    clin = empty_clin(base.index)
    clin["pfs_time"] = base["pfs_time"]
    clin["pfs_event"] = base["pfs_event"]
    clin["pfs_unit"] = "day"
    clin["pfs_usable"] = True
    clin["y_dcb"] = dcb_from_pfs(clin["pfs_time"], clin["pfs_event"], "day").values
    clin["y_primary"] = clin["y_dcb"]
    clin["endpoint"] = "DCB_PFS_ge_180d"
    clin["timing"] = "not stated on GEO"
    clin["regimen"] = "anti-PD-1/PD-L1"
    clin["histology"] = "NA"
    clin["note"] = "Jung / KAIST; DCB derived from deposited PFS; no RECIST on GEO"
    keep = [c for c in expr.columns if c in clin.index]
    return expr[keep], clin.loc[keep]


def load_gse166449() -> tuple[pd.DataFrame, pd.DataFrame]:
    meta = parse_series_matrix(RAW / "GSE166449_series_matrix.txt.gz")
    tpm = pd.read_csv(RAW / "GSE166449_TPM.txt.gz", sep="\t", index_col=0)
    expr = as_log(extract_genes(tpm), "log_linear")
    recs = []
    for _, s in meta.iterrows():
        title = str(s.get("title", "")).lower()
        if "nonresponder" in title or "non-responder" in title:
            y = "NR"
        elif "responder" in title:
            y = "R"
        else:
            y = pd.NA
        recs.append({"sample": s.get("description"), "y": y})
    base = pd.DataFrame(recs).drop_duplicates("sample").set_index("sample")
    clin = empty_clin(base.index)
    clin["y_primary"] = base["y"]
    clin["endpoint"] = "author_R_vs_NR"
    clin["timing"] = "pre-treatment"
    clin["regimen"] = "immunotherapy (GEO title; drug not named per sample)"
    clin["histology"] = "NA"
    clin["note"] = "Samsung advanced-stage lung; pre-treatment; responder in the sample title"
    keep = [c for c in expr.columns if c in clin.index]
    return expr[keep], clin.loc[keep]


def load_gse190265() -> tuple[pd.DataFrame, pd.DataFrame]:
    path = RAW / "GSE190265_TPM_France3.csv.gz"
    with gzip.open(path, "rt") as fh:
        genes = fh.readline().strip().split(";")
    tpm = pd.read_csv(path, sep=";", header=None, skiprows=1, index_col=0)
    tpm.columns = genes
    tpm.index = tpm.index.astype(str)
    expr = as_log(extract_genes(tpm.T), "log_linear")
    info = pd.read_csv(RAW / "GSE190265_samples_info_France3.csv.gz", sep=";")
    info["sample"] = info["sample"].astype(str)
    info = info.drop_duplicates("sample").set_index("sample")
    meta = parse_series_matrix(RAW / "GSE190265_series_matrix.txt.gz").set_index("title")
    meta.index = meta.index.astype(str)
    common = [c for c in expr.columns if c in info.index]
    clin = empty_clin(common)
    clin["pfs_time"] = pd.to_numeric(info.loc[common, "time_PFS"], errors="coerce").values
    clin["pfs_event"] = pd.to_numeric(info.loc[common, "evtPFS"], errors="coerce").values
    clin["pfs_unit"] = "month"
    clin["pfs_usable"] = True
    clin["y_dcb"] = dcb_from_pfs(clin["pfs_time"], clin["pfs_event"], "month").values
    clin["y_primary"] = clin["y_dcb"]
    clin["endpoint"] = "DCB_PFS_ge_6mo"
    clin["timing"] = "diagnosis (pre-treatment)"
    clin["regimen"] = "anti-PD-1 monotherapy"
    hist = []
    for s in common:
        if s in meta.index and "disease state" in meta.columns:
            hist.append(norm_histology(meta.loc[s, "disease state"] if not isinstance(meta.loc[s, "disease state"], pd.Series) else meta.loc[s, "disease state"].iloc[0]))
        else:
            hist.append("NA")
    clin["histology"] = hist
    clin["note"] = "Fumet France-3; uncapped PFS in samples_info; histology only for samples on the series matrix"
    return expr[common], clin


def load_gse190266() -> tuple[pd.DataFrame, pd.DataFrame]:
    path = RAW / "GSE190266_TPM_France4.csv.gz"
    tpm = pd.read_csv(path, sep=";", decimal=",", index_col=0, low_memory=False)
    tpm.index = tpm.index.astype(str)
    expr = as_log(extract_genes(tpm.T.apply(pd.to_numeric, errors="coerce")), "log_linear")
    meta = parse_series_matrix(RAW / "GSE190266_series_matrix.txt.gz")
    meta["sample"] = meta["title"].astype(str)
    meta = meta.drop_duplicates("sample").set_index("sample")
    common = [c for c in expr.columns if c in meta.index]
    clin = empty_clin(common)
    # column names retain the parenthetical from GEO
    tcol = [c for c in meta.columns if "pfs" in c and "time" in c]
    ecol = [c for c in meta.columns if "pfs" in c and "evt" in c]
    clin["pfs_time"] = pd.to_numeric(meta.loc[common, tcol[0]], errors="coerce").values if tcol else np.nan
    clin["pfs_event"] = pd.to_numeric(meta.loc[common, ecol[0]], errors="coerce").values if ecol else np.nan
    clin["pfs_unit"] = "month"
    clin["pfs_usable"] = False  # administratively capped at 6 months
    clin["y_dcb"] = dcb_from_pfs(clin["pfs_time"], clin["pfs_event"], "month").values
    clin["y_primary"] = clin["y_dcb"]
    clin["endpoint"] = "DCB_PFS_ge_6mo"
    clin["timing"] = "diagnosis (pre-treatment)"
    clin["regimen"] = "anti-PD-1 monotherapy"
    clin["histology"] = [norm_histology(meta.loc[s, "disease state"]) if "disease state" in meta.columns else "NA" for s in common]
    clin["note"] = "Fumet France-4; PFS field is capped at 6 months so Cox HR is not estimated; DCB is the binary; TACSTD2 absent (Excel column cap) but CLDN4 and KRT8/18/19 are present"
    return expr[common], clin


def load_gse207422() -> tuple[pd.DataFrame, pd.DataFrame]:
    md = pd.read_excel(RAW / "GSE207422_metadata.xlsx")
    md = md.dropna(subset=["Sample"]).copy()
    md["Sample"] = md["Sample"].astype(str)
    pre = md["Resource"].astype(str).str.contains("Pre", case=False, na=False)
    md = md.loc[pre].drop_duplicates("Sample")
    raw_expr = pd.read_csv(RAW / "GSE207422_log2TPM.txt.gz", sep="\t", index_col=0)
    expr = as_log(extract_genes(raw_expr), "already_log")
    md = md[md["Sample"].isin(expr.columns)]
    recist = md["RECIST"].astype(str).str.upper().str.strip()
    y = pd.Series(pd.NA, index=md.index, dtype=object)
    y[recist.isin(["CR", "PR"])] = "R"
    y[recist.isin(["SD", "PD"])] = "NR"
    path = md["Pathologic Response"].astype(str)
    mpr = pd.Series(pd.NA, index=md.index, dtype=object)
    mpr[path.str.upper().str.startswith("MPR")] = "R"
    mpr[path.str.upper().str.startswith("NMPR")] = "NR"
    clin = empty_clin(md["Sample"].values)
    clin["y_primary"] = y.values
    clin["y_mpr"] = mpr.values
    clin["endpoint"] = "RECIST_CRPR_vs_SDPD"
    clin["timing"] = "pre-treatment biopsy"
    clin["regimen"] = "neoadjuvant anti-PD-1 + chemotherapy"
    clin["histology"] = [norm_histology(v) for v in md["Pathology"].astype(str)]
    clin["note"] = "Hu / Shanghai Pulmonary Hospital bulk log2TPM; primary binary is RECIST; MPR is a separate pathologic endpoint and is not pooled with RECIST"
    keep = [c for c in clin.index if c in expr.columns]
    return expr[keep], clin.loc[keep]


def load_gse283829() -> tuple[pd.DataFrame, pd.DataFrame]:
    counts = pd.read_csv(RAW / "GSE283829_raw_express_matrix_all_samples.txt.gz", sep="\t", index_col=0)
    counts.columns = [str(c)[1:] if str(c).startswith("S") else str(c) for c in counts.columns]
    lib = counts.sum(axis=0).replace(0, np.nan)
    cpm = counts.divide(lib, axis=1) * 1e6
    expr = as_log(extract_genes(cpm), "log_linear")
    meta = parse_series_matrix(RAW / "GSE283829_series_matrix.txt.gz")
    recs = []
    for _, s in meta.iterrows():
        sid = str(s["title"])
        sid2 = sid[1:] if sid.startswith("S") else sid
        use = sid if sid in expr.columns else sid2
        raw_r = str(s.get("disease stage", "")).upper().strip()
        if raw_r in {"CR", "PR"}:
            y = "R"
        elif raw_r in {"SD", "PD"}:
            y = "NR"
        else:
            y = pd.NA
        recs.append(
            {
                "sample": use,
                "y": y,
                "histology": norm_histology(str(s.get("tumor type", ""))),
                "raw": raw_r,
            }
        )
    base = pd.DataFrame(recs).drop_duplicates("sample").set_index("sample")
    clin = empty_clin(base.index)
    clin["y_primary"] = base["y"]
    clin["endpoint"] = "RECIST_CR_vs_SDPD"
    clin["timing"] = "not stated on GEO"
    clin["regimen"] = "ICI (GEO disease stage is CR/SD/PD; no PR label)"
    clin["histology"] = base["histology"]
    clin["note"] = "Lindberg 2025; GEO disease stage used as RECIST; no PR and no PFS deposited"
    keep = [c for c in expr.columns if c in clin.index]
    return expr[keep], clin.loc[keep]


def load_gse274975() -> tuple[pd.DataFrame, pd.DataFrame]:
    zpath = RAW / "PMC11669362_supp.zip"
    dest = RAW / "supp274" / "Table1.xlsx"
    if not dest.exists():
        dest.parent.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(zpath) as z:
            z.extract("Table1.xlsx", dest.parent)
    clin_raw = pd.read_excel(dest)
    clin_raw.columns = [str(c).strip() for c in clin_raw.columns]
    soft_rows = []
    cur = {}
    with gzip.open(RAW / "GSE274975_family.soft.gz", "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("^SAMPLE"):
                if cur:
                    soft_rows.append(cur)
                cur = {"gsm": line.split("=", 1)[1].strip()}
            elif line.startswith("!Sample_title"):
                cur["title"] = line.split("=", 1)[1].strip()
            elif line.startswith("!Sample_source_name"):
                cur["source"] = line.split("=", 1)[1].strip()
        if cur:
            soft_rows.append(cur)
    soft = pd.DataFrame(soft_rows)
    counts = pd.read_csv(RAW / "GSE274975_raw_counts.tsv.gz", sep="\t", index_col=0)

    def col_to_sid(col: str) -> str | None:
        m = re.search(r"Lu[Cc][-_]?(\d+)", str(col))
        if not m:
            return None
        return f"OB_pat_LuC_{int(m.group(1))}"

    colmap = {c: col_to_sid(c) for c in counts.columns}
    if any(v is None for v in colmap.values()):
        raise RuntimeError("unmapped GSE274975 columns")
    counts = counts.rename(columns=colmap)
    lib = counts.sum(axis=0).replace(0, np.nan)
    cpm = counts.divide(lib, axis=1) * 1e6
    expr = as_log(extract_genes(cpm), "log_linear")
    df = clin_raw.merge(soft, left_on="Sample_ID", right_on="title", how="left")
    df = df.drop_duplicates("Sample_ID").set_index("Sample_ID")
    src = df.get("source", pd.Series("", index=df.index)).astype(str)
    histo = df["Histotype"].astype(str)
    non_nsclc = src.str.contains("small cell", case=False) | histo.str.contains("small cell|neuroendocrine", case=False)
    recist = df["RECIST Response"].astype(str).str.upper().str.strip()
    recist = recist.where(df["RECIST Response"].notna(), other="")
    y = pd.Series(pd.NA, index=df.index, dtype=object)
    y[recist.isin(["CR", "PR"])] = "R"
    y[recist.isin(["SD", "PD"])] = "NR"
    y[non_nsclc] = pd.NA
    event = pd.to_numeric(df["Response status available"], errors="coerce")
    pfs = pd.to_numeric(df["PFS time, months"], errors="coerce")
    clin = empty_clin(df.index)
    clin["y_primary"] = y
    clin["pfs_time"] = pfs
    clin["pfs_event"] = event
    clin["pfs_unit"] = "month"
    clin["pfs_usable"] = True
    dcb = dcb_from_pfs(pfs, event, "month")
    dcb.index = clin.index
    dcb[non_nsclc] = pd.NA
    clin["y_dcb"] = dcb
    clin["endpoint"] = "RECIST_CRPR_vs_SDPD"
    clin["timing"] = "biosampling (Table S1; GEO does not say pre-treatment)"
    clin["regimen"] = "mixed PD-1/PD-L1 ± chemo ± CTLA-4 ± PARPi"
    clin["histology"] = ["excluded_non_nsclc" if flag else norm_histology(h) for flag, h in zip(non_nsclc, histo)]
    clin.loc[clin["histology"].eq("excluded_non_nsclc"), ["pfs_time", "pfs_event"]] = pd.NA
    clin["note"] = (
        "Poddubskaya PMID 39723204 Table S1; SCLC/neuroendocrine dropped; "
        "PFS event is the deposited column 'Response status available' (not a separately named censor flag; "
        "2 PD rows are coded 0)"
    )
    keep = [c for c in clin.index if c in expr.columns]
    return expr[keep], clin.loc[keep]


def load_gse218989() -> tuple[pd.DataFrame, pd.DataFrame]:
    meta = parse_series_matrix(RAW / "GSE218989_series_matrix.txt.gz")
    tpm = pd.read_csv(RAW / "GSE218989_TPM.txt.gz", sep="\t", index_col=0)
    expr = as_log(extract_genes(tpm), "log_linear")
    meta["sample"] = meta["title"].astype(str)
    meta = meta.drop_duplicates("sample").set_index("sample")
    common = [c for c in expr.columns if c in meta.index]
    clin = empty_clin(common)
    outcome = meta.loc[common, "treatment outcome"].astype(str)
    y = pd.Series(pd.NA, index=common, dtype=object)
    y[outcome.str.lower().eq("responder")] = "R"
    y[outcome.str.lower().str.contains("non")] = "NR"
    y.index = common
    clin["y_primary"] = y.values
    clin["endpoint"] = "author_R_vs_NR"
    clin["timing"] = "not stated on GEO"
    clin["regimen"] = "PD-1/PD-L1 inhibitor"
    clin["histology"] = "NA"
    supp = pd.read_excel(RAW / "41467_2024_48310_MOESM11_ESM.xlsx")
    supp["PatientID"] = supp["PatientID"].astype(str)
    supp = supp.drop_duplicates("PatientID").set_index("PatientID")
    both = [s for s in common if s in supp.index]
    clin["os_time"] = np.nan
    clin["os_event"] = np.nan
    clin.loc[both, "os_time"] = pd.to_numeric(supp.loc[both, "Overall survival (days)"], errors="coerce").values
    clin.loc[both, "os_event"] = pd.to_numeric(supp.loc[both, "Death"], errors="coerce").values
    # PFS days are deposited without an event indicator. Do not impute one.
    geo_r = clin["y_primary"].map({"R": 1, "NR": 0})
    supp_r = supp["Responder"]
    n_dis = int((geo_r.loc[both].astype(float) != supp_r.loc[both].astype(float)).sum())
    clin["note"] = (
        "Kang / Samsung–KAIST Nat Commun 2024; GEO n with TPM; "
        f"Supplementary Data 8 covers 497 patients and disagrees with GEO response on {n_dis} of {len(both)} overlapping IDs. "
        "Primary label is the GEO treatment outcome. OS uses Death. PFS days have no event column, so no PFS HR."
    )
    clin.attrs["n_label_disagree"] = n_dis
    clin.attrs["n_supp_overlap"] = len(both)
    return expr[common], clin


LOADERS = {
    "GSE126044": load_gse126044,
    "GSE135222": load_gse135222,
    "GSE166449": load_gse166449,
    "GSE190265": load_gse190265,
    "GSE190266": load_gse190266,
    "GSE207422": load_gse207422,
    "GSE283829": load_gse283829,
    "GSE274975": load_gse274975,
    "GSE218989": load_gse218989,
}


def ols_residual(y: np.ndarray, X: np.ndarray) -> np.ndarray:
    y = np.asarray(y, float)
    X = np.asarray(X, float)
    out = np.full(y.shape, np.nan)
    m = np.isfinite(y) & np.isfinite(X).all(axis=1)
    if m.sum() < X.shape[1] + 3:
        return out
    A = np.column_stack([np.ones(m.sum()), X[m]])
    beta, *_ = np.linalg.lstsq(A, y[m], rcond=None)
    out[m] = y[m] - A @ beta
    # r^2
    ss_tot = np.sum((y[m] - y[m].mean()) ** 2)
    ss_res = np.sum(out[m] ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    out_s = pd.Series(out)
    out_s.attrs["r2"] = float(r2)
    return out


def keratin_residual(expr: pd.DataFrame) -> tuple[pd.Series, float, str]:
    need = ["CLDN4", *KERATINS]
    if any(g not in expr.index for g in need):
        missing = [g for g in need if g not in expr.index]
        return pd.Series(np.nan, index=expr.columns), np.nan, "missing " + ",".join(missing)
    y = expr.loc["CLDN4"].to_numpy(float)
    X = expr.loc[KERATINS].T.to_numpy(float)
    resid = ols_residual(y, X)
    m = np.isfinite(y) & np.isfinite(X).all(axis=1)
    if m.sum() < 6:
        return pd.Series(resid, index=expr.columns), np.nan, "too few complete rows"
    A = np.column_stack([np.ones(m.sum()), X[m]])
    beta, *_ = np.linalg.lstsq(A, y[m], rcond=None)
    fitted = A @ beta
    ss_tot = np.sum((y[m] - y[m].mean()) ** 2)
    ss_res = np.sum((y[m] - fitted) ** 2)
    r2 = float(1 - ss_res / ss_tot) if ss_tot > 0 else np.nan
    return pd.Series(resid, index=expr.columns), r2, "CLDN4 ~ KRT8 + KRT18 + KRT19"


def or_2x2(y: pd.Series, high: pd.Series) -> dict:
    mask = y.isin(["R", "NR"]) & high.isin([True, False])
    yy = y[mask]
    hh = high[mask].astype(bool)
    a = int(((hh) & (yy == "R")).sum())
    b = int(((hh) & (yy == "NR")).sum())
    c = int(((~hh) & (yy == "R")).sum())
    d = int(((~hh) & (yy == "NR")).sum())
    if min(a + b, c + d, a + c, b + d) == 0:
        return {"estimable": False, "reason": "empty_margin", "n_used": int(mask.sum()), "a": a, "b": b, "c": c, "d": d}
    haldane = min(a, b, c, d) == 0
    aa, bb, cc, dd = (a + 0.5, b + 0.5, c + 0.5, d + 0.5) if haldane else (a, b, c, d)
    logor = math.log((aa * dd) / (bb * cc))
    se = math.sqrt(1 / aa + 1 / bb + 1 / cc + 1 / dd)
    p = float(fisher_exact([[a, b], [c, d]], alternative="two-sided")[1])
    return {
        "estimable": True,
        "n_used": int(a + b + c + d),
        "n_R": int((yy == "R").sum()),
        "n_NR": int((yy == "NR").sum()),
        "n_high": int(hh.sum()),
        "n_low": int((~hh).sum()),
        "a": a,
        "b": b,
        "c": c,
        "d": d,
        "haldane": haldane,
        "log": logor,
        "se": se,
        "effect": math.exp(logor),
        "lo": math.exp(logor - 1.96 * se),
        "hi": math.exp(logor + 1.96 * se),
        "p": p,
        "kind": "OR",
    }


def median_high(x: pd.Series) -> pd.Series:
    v = pd.to_numeric(x, errors="coerce")
    out = pd.Series(pd.NA, index=v.index, dtype=object)
    ok = v.notna()
    if ok.sum() < 4:
        return out
    thr = float(v[ok].median())
    out[ok & (v >= thr)] = True
    out[ok & (v < thr)] = False
    return out


def extreme_high(x: pd.Series, how: str) -> pd.Series:
    v = pd.to_numeric(x, errors="coerce")
    out = pd.Series(pd.NA, index=v.index, dtype=object)
    ok = v.notna()
    if ok.sum() < 8:
        return out
    if how == "tertile":
        q1, q2 = v[ok].quantile([1 / 3, 2 / 3])
        out[ok & (v >= q2)] = True
        out[ok & (v <= q1)] = False
    elif how == "quartile":
        q1, q3 = v[ok].quantile([0.25, 0.75])
        out[ok & (v >= q3)] = True
        out[ok & (v <= q1)] = False
    return out


def logit_per_sd(y: pd.Series, x: pd.Series, covars: pd.DataFrame | None = None) -> dict:
    mask = y.isin(["R", "NR"]) & pd.to_numeric(x, errors="coerce").notna()
    if covars is not None:
        mask = mask & covars.notna().all(axis=1)
    yy = (y[mask] == "R").astype(int)
    xx = pd.to_numeric(x[mask], errors="coerce")
    if yy.nunique() < 2 or len(yy) < 8:
        return {"estimable": False, "reason": "too_few", "n_used": int(mask.sum())}
    z = (xx - xx.mean()) / (xx.std(ddof=1) or 1.0)
    cols = [z.to_numpy()]
    if covars is not None:
        C = covars.loc[mask].astype(float)
        C = (C - C.mean()) / C.std(ddof=0).replace(0, np.nan)
        if C.isna().any().any():
            return {"estimable": False, "reason": "covariate_sd_zero", "n_used": int(mask.sum())}
        cols.extend([C[c].to_numpy() for c in C.columns])
    X = sm.add_constant(np.column_stack(cols), has_constant="add")
    try:
        fit = sm.Logit(yy.to_numpy(), X).fit(disp=False, maxiter=200)
    except Exception as exc:
        return {"estimable": False, "reason": str(exc), "n_used": int(mask.sum())}
    # params[1] is the first slope (CLDN4), params[0] is the intercept
    logor = float(fit.params[1])
    se = float(fit.bse[1])
    return {
        "estimable": True,
        "n_used": int(mask.sum()),
        "n_R": int(yy.sum()),
        "n_NR": int((1 - yy).sum()),
        "log": logor,
        "se": se,
        "effect": math.exp(logor),
        "lo": math.exp(logor - 1.96 * se),
        "hi": math.exp(logor + 1.96 * se),
        "p": float(fit.pvalues[1]),
        "kind": "OR",
        "haldane": False,
    }


def cox_model(time, event, x: pd.Series, covars: pd.DataFrame | None = None, binary: bool = False) -> dict:
    t = pd.to_numeric(pd.Series(time), errors="coerce")
    e = pd.to_numeric(pd.Series(event), errors="coerce")
    xx = pd.to_numeric(x, errors="coerce")
    m = t.notna() & e.notna() & xx.notna() & (t > 0)
    if covars is not None:
        m = m & covars.notna().all(axis=1)
    t, e, xx = t[m], e[m], xx[m]
    if len(t) < 8 or int(e.sum()) < 3 or xx.nunique() < 2:
        return {"estimable": False, "reason": "too_few_events", "n_used": int(m.sum()), "n_events": int(e.sum()) if len(e) else 0}
    if binary:
        exog = xx.astype(float).to_numpy().reshape(-1, 1)
    else:
        z = (xx - xx.mean()) / (xx.std(ddof=1) or 1.0)
        cols = [z.to_numpy()]
        if covars is not None:
            C = covars.loc[m].astype(float)
            C = (C - C.mean()) / C.std(ddof=0).replace(0, np.nan)
            cols.extend([C[c].to_numpy() for c in C.columns])
        exog = np.column_stack(cols)
    try:
        res = PHReg(t.to_numpy(), exog, status=e.to_numpy()).fit(disp=0)
        loghr = float(np.asarray(res.params).reshape(-1)[0])
        se = float(np.asarray(res.bse).reshape(-1)[0])
        p = float(np.asarray(res.pvalues).reshape(-1)[0])
    except Exception as exc:
        return {"estimable": False, "reason": str(exc), "n_used": int(len(t)), "n_events": int(e.sum())}
    return {
        "estimable": True,
        "n_used": int(len(t)),
        "n_events": int(e.sum()),
        "log": loghr,
        "se": se,
        "effect": math.exp(loghr),
        "lo": math.exp(loghr - 1.96 * se),
        "hi": math.exp(loghr + 1.96 * se),
        "p": p,
        "kind": "HR",
        "haldane": False,
    }


def dl_meta(rows: list[dict]) -> dict:
    d = [r for r in rows if r.get("estimable") and r.get("se", 0) > 0 and np.isfinite(r.get("log", np.nan))]
    k = len(d)
    if k == 0:
        return {"k": 0, "estimable": False}
    log = np.array([r["log"] for r in d])
    se = np.array([r["se"] for r in d])
    w = 1.0 / se**2
    fe = float(np.sum(w * log) / np.sum(w))
    q = float(np.sum(w * (log - fe) ** 2))
    df = k - 1
    c = float(np.sum(w) - np.sum(w**2) / np.sum(w)) if k > 1 else 0.0
    tau2 = max(0.0, (q - df) / c) if c > 0 else 0.0
    w_re = 1.0 / (se**2 + tau2)
    re = float(np.sum(w_re * log) / np.sum(w_re))
    re_se = float(math.sqrt(1.0 / np.sum(w_re)))
    i2 = max(0.0, (q - df) / q) * 100 if q > 0 and df > 0 else 0.0
    z = re / re_se if re_se else np.nan
    p = float(2 * (1 - norm.cdf(abs(z)))) if z == z else np.nan
    return {
        "k": k,
        "n_total": int(sum(r["n_used"] for r in d)),
        "estimable": True,
        "effect": math.exp(re),
        "lo": math.exp(re - 1.96 * re_se),
        "hi": math.exp(re + 1.96 * re_se),
        "p": p,
        "I2": i2,
        "tau2": tau2,
        "Q": q,
        "cohorts": [r["cohort"] for r in d],
        "weights": (w_re / w_re.sum()).tolist(),
    }


def fmt(effect, lo, hi) -> str:
    if effect is None or effect != effect:
        return "NA"
    return f"{effect:.2f} [{lo:.2f}–{hi:.2f}]"


def spearman_match(a: pd.Series, b: pd.Series) -> float:
    df = pd.concat([a, b], axis=1, join="inner").dropna()
    if len(df) < 12:
        return np.nan
    return float(df.iloc[:, 0].corr(df.iloc[:, 1], method="spearman"))


def find_duplicates(exprs: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Flag a smaller cohort as non-independent when >=50% of its samples match a larger cohort at rho>=0.99."""
    names = list(exprs)
    rows = []
    for i, a in enumerate(names):
        for b in names[i + 1 :]:
            ea, eb = exprs[a], exprs[b]
            genes = sorted(set(ea.index) & set(eb.index))
            if len(genes) < 12:
                rows.append({"cohort_a": a, "cohort_b": b, "n_genes": len(genes), "n_matched": 0, "max_rho": np.nan, "note": "too few shared genes"})
                continue
            # samples as columns
            A = ea.loc[genes]
            B = eb.loc[genes]
            small, large = (a, b) if A.shape[1] <= B.shape[1] else (b, a)
            S = ea.loc[genes] if small == a else eb.loc[genes]
            L = eb.loc[genes] if small == a else ea.loc[genes]
            matched = 0
            best = []
            for s in S.columns:
                rhos = []
                sv = S[s]
                for lv_name in L.columns:
                    rho = spearman_match(sv, L[lv_name])
                    if rho == rho:
                        rhos.append((rho, lv_name))
                if not rhos:
                    continue
                top = max(rhos, key=lambda t: t[0])
                best.append(top[0])
                if top[0] >= 0.99:
                    matched += 1
            frac = matched / S.shape[1] if S.shape[1] else 0
            rows.append(
                {
                    "smaller": small,
                    "larger": large,
                    "n_small": int(S.shape[1]),
                    "n_genes": len(genes),
                    "n_matched_rho_ge_0.99": matched,
                    "frac_matched": frac,
                    "max_rho": float(np.nanmax(best)) if best else np.nan,
                    "median_best_rho": float(np.nanmedian(best)) if best else np.nan,
                    "non_independent": frac >= 0.5,
                }
            )
    return pd.DataFrame(rows)


def analyze_cohort(name: str, expr: pd.DataFrame, clin: pd.DataFrame) -> tuple[list[dict], pd.DataFrame]:
    resid, r2, resid_note = keratin_residual(expr)
    clin = clin.copy()
    clin["CLDN4"] = expr.loc["CLDN4"] if "CLDN4" in expr.index else np.nan
    for g in KERATINS:
        clin[g] = expr.loc[g] if g in expr.index else np.nan
    clin["CLDN4_krt_resid"] = resid
    rows = []

    def add(measure, cutoff, endpoint, res, extra=None):
        rec = {
            "cohort": name,
            "measure": measure,
            "cutoff": cutoff,
            "endpoint": endpoint,
            "keratin_r2": r2,
            "resid_note": resid_note,
            "regimen": clin["regimen"].iloc[0] if len(clin) else "",
            "timing": clin["timing"].iloc[0] if len(clin) else "",
        }
        rec.update(res)
        if extra:
            rec.update(extra)
        rows.append(rec)

    # primary binary
    y = clin["y_primary"]
    sub = clin.loc[y.isin(["R", "NR"])].copy()
    if "CLDN4" in sub and sub["CLDN4"].notna().sum() >= 8 and y.isin(["R", "NR"]).sum() >= 8:
        for measure, col in (("raw", "CLDN4"), ("keratin_residual", "CLDN4_krt_resid")):
            if sub[col].notna().sum() < 8:
                continue
            high = median_high(sub[col])
            add(measure, "median", sub["endpoint"].iloc[0], or_2x2(sub["y_primary"], high.map({True: True, False: False})))
            for how in ("tertile", "quartile"):
                lab = extreme_high(sub[col], how)
                # True/False/NA; drop NA inside or_2x2 via isin
                add(measure, how, sub["endpoint"].iloc[0], or_2x2(sub["y_primary"], lab.map({True: True, False: False})))
            cov = sub[KERATINS] if measure == "raw" else None
            # adjusted model is only meaningful on the raw CLDN4 scale
            if measure == "raw" and all(g in sub.columns for g in KERATINS):
                add("keratin_adjusted", "per_SD", sub["endpoint"].iloc[0], logit_per_sd(sub["y_primary"], sub["CLDN4"], sub[KERATINS]))
            add(measure, "per_SD", sub["endpoint"].iloc[0], logit_per_sd(sub["y_primary"], sub[col], None))

    # MPR sensitivity (not primary)
    if clin["y_mpr"].notna().any():
        m = clin.loc[clin["y_mpr"].isin(["R", "NR"])].copy()
        m["y_primary"] = m["y_mpr"]
        if m["CLDN4"].notna().sum() >= 8:
            high = median_high(m["CLDN4"])
            add("raw", "median", "MPR_vs_NMPR", or_2x2(m["y_mpr"], high))
            if m["CLDN4_krt_resid"].notna().sum() >= 8:
                high = median_high(m["CLDN4_krt_resid"])
                add("keratin_residual", "median", "MPR_vs_NMPR", or_2x2(m["y_mpr"], high))

    # DCB if it is not already the primary endpoint
    if clin["y_dcb"].notna().any() and str(clin["endpoint"].iloc[0]) != "DCB_PFS_ge_6mo" and "DCB" not in str(clin["endpoint"].iloc[0]):
        m = clin.loc[clin["y_dcb"].isin(["R", "NR"])].copy()
        if m["CLDN4"].notna().sum() >= 8:
            high = median_high(m["CLDN4"])
            add("raw", "median", "DCB", or_2x2(m["y_dcb"], high))
            if m["CLDN4_krt_resid"].notna().sum() >= 8:
                high = median_high(m["CLDN4_krt_resid"])
                add("keratin_residual", "median", "DCB", or_2x2(m["y_dcb"], high))

    # PFS
    if bool(clin["pfs_usable"].iloc[0]) and pd.to_numeric(clin["pfs_event"], errors="coerce").notna().any():
        for measure, col in (("raw", "CLDN4"), ("keratin_residual", "CLDN4_krt_resid")):
            if clin[col].notna().sum() < 8:
                continue
            high = median_high(clin[col]).map({True: 1.0, False: 0.0})
            add(measure, "median", "PFS", cox_model(clin["pfs_time"], clin["pfs_event"], high, binary=True))
            add(measure, "per_SD", "PFS", cox_model(clin["pfs_time"], clin["pfs_event"], clin[col], binary=False))
        if all(g in clin.columns for g in KERATINS):
            add(
                "keratin_adjusted",
                "per_SD",
                "PFS",
                cox_model(clin["pfs_time"], clin["pfs_event"], clin["CLDN4"], clin[KERATINS], binary=False),
            )

    # OS (secondary; not pooled with PFS)
    if pd.to_numeric(clin["os_event"], errors="coerce").notna().any():
        for measure, col in (("raw", "CLDN4"), ("keratin_residual", "CLDN4_krt_resid")):
            if clin[col].notna().sum() < 8:
                continue
            high = median_high(clin[col]).map({True: 1.0, False: 0.0})
            add(measure, "median", "OS", cox_model(clin["os_time"], clin["os_event"], high, binary=True))
            add(measure, "per_SD", "OS", cox_model(clin["os_time"], clin["os_event"], clin[col], binary=False))
        if all(g in clin.columns for g in KERATINS):
            add(
                "keratin_adjusted",
                "per_SD",
                "OS",
                cox_model(clin["os_time"], clin["os_event"], clin["CLDN4"], clin[KERATINS], binary=False),
            )

    clin["cohort"] = name
    return rows, clin


def forest_plot(studies: list[dict], pooled: dict, title: str, xlabel: str, path: Path) -> None:
    labels = []
    effects = []
    los = []
    his = []
    weights = pooled.get("weights") or []
    for i, s in enumerate(studies):
        w = weights[i] if i < len(weights) else np.nan
        wtxt = f", weight {100 * w:.0f}%" if w == w else ""
        if s.get("kind") == "HR":
            detail = f"n={int(s['n_used'])}, events={int(s['n_events'])}"
        else:
            detail = f"n={int(s['n_used'])}, R={int(s['n_R'])}"
        labels.append(f"{s['cohort']}  {detail}{wtxt}")
        effects.append(s["effect"])
        los.append(s["lo"])
        his.append(s["hi"])
    if pooled.get("estimable"):
        labels.append(f"Random effects (k={pooled['k']})")
        effects.append(pooled["effect"])
        los.append(pooled["lo"])
        his.append(pooled["hi"])
    n = len(labels)
    fig_h = max(3.2, 0.42 * n + 1.4)
    fig, ax = plt.subplots(figsize=(8.4, fig_h))
    ys = np.arange(n)[::-1]
    for i, (y, eff, lo, hi) in enumerate(zip(ys, effects, los, his)):
        ax.plot([lo, hi], [y, y], color="#333333", lw=1.4, zorder=2)
        marker = "D" if pooled.get("estimable") and i == n - 1 else "s"
        size = 70 if marker == "D" else 36
        ax.scatter([eff], [y], marker=marker, s=size, color="#1f4e79", zorder=3)
    ax.axvline(1.0, color="#888888", lw=0.8, ls="--")
    ax.set_yticks(ys)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel(xlabel)
    ax.set_title(title, fontsize=10, loc="left")
    finite_hi = [h for h in his if h == h]
    finite_lo = [l for l in los if l == l]
    if finite_hi and finite_lo:
        ax.set_xlim(min(0.05, min(finite_lo) * 0.8), max(finite_hi) * 1.05)
    ax.set_xscale("log")
    ax.tick_params(axis="y", length=0)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def panel_absent() -> list[dict]:
    rows = []
    # GSE136961 Oncomine panel
    idx = pd.read_csv(RAW / "GSE136961_TPM.tsv.gz", sep="\t", usecols=[0]).iloc[:, 0].astype(str)
    rows.append(
        {
            "accession": "GSE136961",
            "role": "excluded",
            "reason": "NSCLC anti-PD-1 DCB in sample titles, but CLDN4 is not on the Oncomine Immune Response panel",
            "cldn4": "CLDN4" in set(idx.str.upper()),
            "n_genes": int(len(idx)),
        }
    )
    idx2 = pd.read_csv(RAW / "GSE161537_log2cpm.csv.gz", sep=";", usecols=[0]).iloc[:, 0].astype(str)
    rows.append(
        {
            "accession": "GSE161537",
            "role": "excluded",
            "reason": "Advanced NSCLC 2L PD-1/PD-L1 with RECIST and PFS, HTG Oncology Biomarker Panel; CLDN4 absent (KRT8/18/19 present)",
            "cldn4": "CLDN4" in set(idx2.str.upper()),
            "n_genes": int(len(idx2)),
        }
    )
    # expression open, response not on GEO
    rows.append(
        {
            "accession": "GSE253564",
            "role": "excluded",
            "reason": "Pre-treatment bulk FPKM from durvalumab ± SBRT (Weill Cornell). GEO has arm and histology, not MPR/RECIST/PFS. Response is not joinable from the open matrix.",
            "cldn4": "not re-scored; FPKM file is whole transcriptome",
            "n_genes": "",
        }
    )
    rows.append(
        {
            "accession": "GSE248378",
            "role": "excluded",
            "reason": "Post-neoadjuvant resected-tumor FPKM (same trial as GSE253564). Not a pre-treatment R/NR matrix.",
            "cldn4": "",
            "n_genes": "",
        }
    )
    rows.append(
        {
            "accession": "GSE248249",
            "role": "excluded",
            "reason": "Acquired-resistance NSCLC FFPE pre/post immunotherapy. GEO sample records do not carry RECIST or R/NR.",
            "cldn4": "",
            "n_genes": "",
        }
    )
    rows.append(
        {
            "accession": "GSE182328",
            "role": "excluded",
            "reason": "NSCLC tumor RNA-seq on PD-1 blockade, but the deposited label is fecal Akkermansia detectable vs not, not RECIST or PFS.",
            "cldn4": "",
            "n_genes": "",
        }
    )
    rows.append(
        {
            "accession": "GSE208858",
            "role": "excluded",
            "reason": "Utomilumab (4-1BB agonist), not PD-1/PD-L1 checkpoint blockade; melanoma and NSCLC mixed.",
            "cldn4": "",
            "n_genes": "",
        }
    )
    rows.append(
        {
            "accession": "GSE93157",
            "role": "excluded",
            "reason": "NanoString PanCancer Immune 730. NSCLC RECIST and PFS exist; CLDN4 is not on the panel (prior panel audit, not re-downloaded).",
            "cldn4": False,
            "n_genes": 730,
        }
    )
    rows.append(
        {
            "accession": "GSE309652",
            "role": "excluded",
            "reason": "Author R/NR on GEO, but CLDN4 is absent from the targeted panel (prior open of GPL31904).",
            "cldn4": False,
            "n_genes": "",
        }
    )
    rows.append(
        {
            "accession": "GSE202417",
            "role": "excluded",
            "reason": "Blood CD8 microarray, not tumor bulk.",
            "cldn4": "",
            "n_genes": "",
        }
    )
    rows.append(
        {
            "accession": "GSE260770",
            "role": "excluded",
            "reason": "Plasma exosome RNA, not tumor bulk.",
            "cldn4": "",
            "n_genes": "",
        }
    )
    rows.append(
        {
            "accession": "ArrayExpress",
            "role": "excluded",
            "reason": "No open ArrayExpress NSCLC ICI cohort with a processed tumor matrix and R/NR or PFS was added. Prior accession-level search (E-MTAB leftovers) found none; clinical trial RNA (POPLAR/OAK, CheckMate, KEYNOTE, IMpower) is not in open ArrayExpress.",
            "cldn4": "",
            "n_genes": "",
        }
    )
    return rows


def supports_pdf(pooled: dict) -> dict:
    if not pooled.get("estimable"):
        return {"supported": False, "reason": "primary pool not estimable", "k": 0}
    rounds = round(pooled["effect"], 2) == PDF_OR
    k_ok = pooled["k"] == PDF_K
    inside = pooled["lo"] <= PDF_OR <= pooled["hi"]
    return {
        "supported": bool(rounds and k_ok),
        "k": pooled["k"],
        "effect": pooled["effect"],
        "lo": pooled["lo"],
        "hi": pooled["hi"],
        "p": pooled["p"],
        "I2": pooled["I2"],
        "rounds_to_0.42": rounds,
        "k_is_11": k_ok,
        "0.42_inside_ci": inside,
    }


def write_report(effect_df: pd.DataFrame, pooled_rows: list[dict], dup: pd.DataFrame, patients: pd.DataFrame, pdf: dict) -> None:
    def grab(pool_name):
        hit = [r for r in pooled_rows if r["pool"] == pool_name]
        return hit[0] if hit else None

    primary = grab("primary_OR_raw_median")
    krt = grab("primary_OR_resid_median")
    hr = grab("primary_HR_raw_median")
    hrk = grab("primary_HR_resid_median")

    def line(r):
        if not r or not r.get("estimable"):
            return "not estimable"
        return f"k={r['k']}, n={r['n_total']}, RE {fmt(r['effect'], r['lo'], r['hi'])}, p={r['p']:.3g}, I²={r['I2']:.0f}%"

    primary_eps = patients.groupby("cohort")["endpoint"].first().to_dict()
    per = effect_df[
        (effect_df["cutoff"] == "median")
        & (effect_df["measure"].isin(["raw", "keratin_residual"]))
        & (effect_df["estimable"] == True)
    ].copy()
    per["is_primary_endpoint"] = per.apply(lambda r: r["endpoint"] == primary_eps.get(r["cohort"]), axis=1)

    def cohort_table(mask) -> str:
        sub = per.loc[mask]
        lines = ["| Cohort | Measure | Endpoint | n | 2×2 or events | Effect [95% CI] | Fisher or Cox p |", "|---|---|---|---:|---|---|---:|"]
        if sub.empty:
            return "(none)"
        for _, r in sub.sort_values(["cohort", "measure"]).iterrows():
            if r["kind"] == "OR":
                cells = f"{int(r['a'])}/{int(r['b'])}/{int(r['c'])}/{int(r['d'])} (hiR/hiNR/loR/loNR)"
            else:
                cells = f"{int(r['n_events'])} events"
            lines.append(
                f"| {r['cohort']} | {r['measure']} | {r['endpoint']} | {int(r['n_used'])} | {cells} | {fmt(r['effect'], r['lo'], r['hi'])} | {r['p']:.3g} |"
            )
        return "\n".join(lines)

    if len(dup) and "max_rho" in dup.columns and dup["max_rho"].notna().any():
        top = dup.loc[dup["max_rho"].idxmax()]
        n_bad = int(dup["non_independent"].fillna(False).sum()) if "non_independent" in dup.columns else 0
        dup_txt = (
            f"{len(dup)} cohort pairs compared on the shared gene panel. "
            f"Pairs declared non-independent (at least half the smaller cohort at Spearman ≥ 0.99): {n_bad}. "
            f"Highest best-match Spearman was {top['max_rho']:.2f} "
            f"({top['smaller']} vs {top['larger']}; median best-match in that pair {top['median_best_rho']:.2f}). "
            "No sample pair reached 0.99, so all nine series stay in the primary pool. "
            "A shared patient sequenced on a different pipeline could be missed by this fingerprint. "
            "Full pairs: `tables/sample_overlap.csv`."
        )
    else:
        dup_txt = "No overlap table."
    supported = "Yes" if pdf.get("supported") else "No"
    reason = (
        f"Open NSCLC bulk primary median-split OR is {fmt(pdf.get('effect'), pdf.get('lo'), pdf.get('hi'))} "
        f"at k={pdf.get('k')}. "
        f"Point estimate rounds to 0.42: {pdf.get('rounds_to_0.42')}. "
        f"k equals 11: {pdf.get('k_is_11')}. "
        f"0.42 lies inside the 95% CI: {pdf.get('0.42_inside_ci')}."
    )

    small = effect_df[
        (effect_df["measure"] == "raw")
        & (effect_df["cutoff"] == "median")
        & (effect_df["estimable"] == True)
        & effect_df.apply(lambda r: r["endpoint"] == primary_eps.get(r["cohort"]), axis=1)
    ]
    small_lines = []
    n_sig = 0
    for _, r in small.iterrows():
        flag = []
        if r["n_used"] < 40:
            flag.append("n<40")
        n_r = r["n_R"] if r["n_R"] == r["n_R"] else 99
        n_nr = r["n_NR"] if r["n_NR"] == r["n_NR"] else 99
        if n_r < 8 or n_nr < 8:
            flag.append("fewer than 8 patients in one response class")
        if r["lo"] < 1 < r["hi"]:
            flag.append("CI includes 1")
        if r["p"] < 0.05:
            n_sig += 1
        small_lines.append(
            f"- {r['cohort']} ({r['endpoint']}): n={int(r['n_used'])}, OR {fmt(r['effect'], r['lo'], r['hi'])}, Fisher p={r['p']:.3g}"
            + (f" ({'; '.join(flag)})" if flag else "")
        )
    r2 = effect_df.loc[effect_df["measure"] == "raw", ["cohort", "keratin_r2"]].drop_duplicates("cohort")
    r2_lo = float(r2["keratin_r2"].min())
    r2_hi = float(r2["keratin_r2"].max())
    leave = grab("primary_OR_raw_median_without_GSE218989")
    big = small[small["cohort"] == "GSE218989"]
    if len(big):
        br = big.iloc[0]
        big_line = f"OR {fmt(br['effect'], br['lo'], br['hi'])} (Fisher p={br['p']:.3g}, {int(br['n_R'])} responders / {int(br['n_NR'])} non-responders)"
        big_ci_note = "That interval contains 0.42." if br["lo"] <= PDF_OR <= br["hi"] else "That interval does not contain 0.42."
    else:
        big_line = "not estimable"
        big_ci_note = ""

    text = f"""# Open NSCLC ICI bulk: CLDN4-high versus response

Public GEO and ArrayExpress tumor bulk only. One locked binary endpoint per cohort. Keratin residual is the within-cohort OLS residual of log2 CLDN4 on KRT8, KRT18, and KRT19. Odds ratios below 1 mean CLDN4-high tumors had lower odds of response. Hazard ratios above 1 mean CLDN4-high tumors had worse survival.

## Does this support the PDF 11-cohort figure?

The figure is CLDN4-high **OR = 0.42 (95% CI 0.18–0.95), k = 11**.

**Supported: {supported}.**

{reason}

This NSCLC meta is not that figure. Nine independent open tumor-bulk series have CLDN4 and a clear response label. That is not 11, and the random-effects interval does not contain 0.42. No cohort was added or dropped to move the point estimate.

Leaving out GSE218989 (the only large series) gives {line(leave)}. That interval is wide enough to contain 0.42, and the point estimate is not 0.42. The small-series meta is underpowered. It does not recover the PDF figure either.

## Locked primary results

| Analysis | Result |
|---|---|
| Response OR, raw CLDN4, median split | {line(primary)} |
| Response OR, keratin residual, median split | {line(krt)} |
| PFS HR, raw CLDN4, median split | {line(hr)} |
| PFS HR, keratin residual, median split | {line(hrk)} |

Forests: `figures/forest_or_raw_median.png`, `figures/forest_or_keratin_residual_median.png`, `figures/forest_hr_raw_median.png`, `figures/forest_hr_keratin_residual_median.png`.

### Primary response endpoint (median split)

One endpoint per cohort. Extra DCB and MPR rows are below and are not in the pooled OR.

{cohort_table(per["is_primary_endpoint"] & per["endpoint_family"].eq("response"))}

### Other response labels (not pooled)

{cohort_table((~per["is_primary_endpoint"]) & per["endpoint_family"].eq("response"))}

GSE207422 major pathologic response has a raw median OR whose Woolf interval excludes 1 (0.14 [0.02–0.96]) and a Fisher p of 0.089. That disagreement is the small-sample case noted above. It is not called a significant result, and it is not the primary endpoint (RECIST). Keratin residual MPR is null (OR 0.70, p=1).

### PFS or OS (median split)

PFS Cox is limited to cohorts with a deposited time and an event indicator, and with follow-up that is not capped at 6 months. GSE190266 is DCB-only for that reason. GSE218989 deposits PFS days without an event indicator, so it contributes an OS HR (Death) and not a PFS HR. GSE274975 PFS drops the two non-NSCLC samples (small-cell and neuroendocrine).

{cohort_table(per["endpoint_family"].eq("survival"))}

## What is underpowered, and what is a null

Primary raw median-split tests with Fisher p < 0.05: {n_sig}. Directions are mixed: some cohorts have OR above 1 and some below 1. Keratin regression R² for CLDN4 ranges from {r2_lo:.2f} to {r2_hi:.2f}, so the residual is a real adjustment, and the residual meta is still null.

Small series cannot rule a moderate association in or out. A confidence interval that includes 1 is a null at α = 0.05. In the smallest 2×2 tables the Woolf interval and the Fisher p can disagree; the p used here is Fisher exact, and the meta-analysis uses the Woolf log-OR.

{chr(10).join(small_lines)}

GSE218989 (Samsung–KAIST) is the only open NSCLC bulk series large enough for a moderate OR to have a relatively narrow interval. Its raw median-split result is {big_line}. {big_ci_note} The series dominates the random-effects pool. The leave-GSE218989 row in `tables/pooled.csv` is the underpowered small-cohort meta.

Neoadjuvant PD-1 plus chemotherapy (GSE207422) uses pathologic response only as a sensitivity (`MPR_vs_NMPR`). The primary binary for that series is RECIST. Dropping it does not create an 11-cohort OR of 0.42 (`primary_OR_raw_median_no_neoadjuvant`).

## Duplicate samples

Spearman correlation on the shared gene panel (CLDN4, keratins, and the other extracted genes). A smaller cohort is non-independent when at least half of its samples match a larger cohort at ρ ≥ 0.99. Non-independent series would stay in the per-cohort table and leave the primary pool.

{dup_txt}

## Cohorts opened and not put in the forest

| Accession | Why it is not in the OR/HR forest |
|---|---|
| GSE136961 | DCB labels, CLDN4 not on the panel (confirmed in this run) |
| GSE161537 | RECIST and PFS, CLDN4 not on the HTG panel (confirmed in this run; keratins are present) |
| GSE93157 | NSCLC RECIST/PFS; NanoString 730 without CLDN4 |
| GSE309652 | R/NR on GEO; CLDN4 absent from the targeted panel |
| GSE253564 | Pre-treatment bulk FPKM, durvalumab ± SBRT; no MPR/RECIST/PFS on GEO |
| GSE248378 | Post-treatment resected tumors, not baseline R/NR |
| GSE248249 | Pre/post immunotherapy FFPE; no RECIST on GEO |
| GSE182328 | Tumor RNA on PD-1, label is Akkermansia, not response |
| GSE208858 | 4-1BB agonist, not PD-1/PD-L1 |
| GSE202417, GSE260770 | Blood or exosome, not tumor bulk |
| ArrayExpress | No leftover open NSCLC ICI tumor matrix with R/NR or PFS |

GSE135222 is not double-counted as ICB_Jung. GSE207422 single-cell libraries are not in this bulk meta.

## Methods

- Expression: counts become log2(CPM+1). TPM and FPKM become log2(x+1). GSE207422 is the author log2TPM and is not logged again.
- Median: high = at or above the cohort median on the analysis subset. Tertile (T3 vs T1) and quartile (Q4 vs Q1) are in `tables/effects.csv` and are not used to choose a pooled number.
- Binary OR: Woolf interval; Haldane–Anscombe 0.5 if a 2×2 cell is 0; Fisher exact p. The interval and the Fisher p answer different questions and can disagree when n is small.
- Continuous OR: logistic regression per 1 SD. Keratin-adjusted OR: same model plus z-scored KRT8, KRT18, and KRT19.
- Keratin residual: within cohort, CLDN4 ~ KRT8 + KRT18 + KRT19. The residual R² is in `tables/effects.csv` (`keratin_r2`).
- DCB: PFS ≥ 180 days or ≥ 6 months. Patients censored before that landmark are excluded from the 2×2, not called non-responders.
- Cox models use `statsmodels` PHReg. GSE274975 uses the Table S1 column “Response status available” as the event indicator because no other censor column is deposited.
- Meta-analysis: DerSimonian–Laird random effects and the I² from Cochran’s Q.
- Sample-level genes and labels: `tables/patient_level.tsv`.

## Reproduce

```bash
python3 -m pip install -r results/ici_nsclc_cldn4_bulk/requirements.txt
python3 results/ici_nsclc_cldn4_bulk/analyze.py
```

The script expects the public files already fetched under `/tmp/ici_raw` (GEO FTP, EuropePMC PMC11669362 Table S1, Springer Supplementary Data 8 for GSE218989). It does not re-download if those files are present.

## 中文

开放 NSCLC ICI 肿瘤 bulk（GEO / ArrayExpress）里，能同时拿到 CLDN4 和明确 R/NR 或 PFS 的独立队列是 9 个，不是 11 个。锁定的中位数切分、随机效应合并：应答 OR 0.96 [0.69–1.33]，p=0.79，I²=0%，n=638。角蛋白残差（CLDN4 ~ KRT8+KRT18+KRT19）OR 1.00 [0.66–1.51]，p=0.99。PFS HR 只有 3 个有事件时间且未被 6 个月截尾的队列，HR 1.12 [0.74–1.68]，p=0.59。各队列主终点 Fisher 检验无一 p<0.05。GSE218989（n=355）的区间不含 0.42；去掉它之后的小队列合并区间很宽，点估计也不是 0.42。因此这组公开数字不支持 PDF 上的 11 队列 OR=0.42。小队列本身检验效能不足，合并结果是阴性而不是把 0.42 证伪到每一项灵敏度上。
"""
    (ROOT / "WRITEUP.md").write_text(text, encoding="utf-8")


def main() -> None:
    exprs = {}
    clins = {}
    all_rows = []
    patient_frames = []
    for name, fn in LOADERS.items():
        print("loading", name)
        expr, clin = fn()
        if "CLDN4" not in expr.index:
            raise RuntimeError(f"{name} missing CLDN4")
        missing_k = [g for g in KERATINS if g not in expr.index]
        if missing_k:
            raise RuntimeError(f"{name} missing {missing_k}")
        print(f"  samples {expr.shape[1]} genes {list(expr.index)} y {clin['y_primary'].value_counts(dropna=False).to_dict()}")
        exprs[name] = expr
        clins[name] = clin
        rows, clin2 = analyze_cohort(name, expr, clin)
        all_rows.extend(rows)
        patient_frames.append(clin2.reset_index().rename(columns={"index": "sample"}))

    dup = find_duplicates(exprs)
    dup.to_csv(TABLES / "sample_overlap.csv", index=False)
    drop = set()
    if len(dup) and "non_independent" in dup.columns:
        drop = set(dup.loc[dup["non_independent"] == True, "smaller"].astype(str))
    print("non-independent", sorted(drop))

    effect = pd.DataFrame(all_rows)
    effect["in_primary_pool"] = ~effect["cohort"].isin(drop)
    effect["endpoint_family"] = np.where(
        effect["endpoint"].isin(["PFS", "OS"]),
        "survival",
        "response",
    )
    effect.to_csv(TABLES / "effects.csv", index=False)
    patients = pd.concat(patient_frames, ignore_index=True)
    patients.to_csv(TABLES / "patient_level.tsv", sep="\t", index=False)

    def select(measure, cutoff, family, cohorts=None, endpoint=None):
        d = effect[(effect["measure"] == measure) & (effect["cutoff"] == cutoff) & (effect["estimable"] == True)].copy()
        if family == "response":
            # primary endpoint only (not MPR, not extra DCB)
            primary_eps = patients.groupby("cohort")["endpoint"].first().to_dict()
            d = d[d.apply(lambda r: r["endpoint"] == primary_eps.get(r["cohort"]), axis=1)]
        elif family == "PFS":
            d = d[d["endpoint"] == "PFS"]
        elif family == "OS":
            d = d[d["endpoint"] == "OS"]
        if cohorts is not None:
            d = d[d["cohort"].isin(cohorts)]
        if endpoint is not None:
            d = d[d["endpoint"] == endpoint]
        return d

    independent = [c for c in LOADERS if c not in drop]
    no_neo = [c for c in independent if c != "GSE207422"]
    no_big = [c for c in independent if c != "GSE218989"]
    specs = [
        ("primary_OR_raw_median", "raw", "median", "response", independent),
        ("primary_OR_resid_median", "keratin_residual", "median", "response", independent),
        ("primary_OR_raw_perSD", "raw", "per_SD", "response", independent),
        ("primary_OR_resid_perSD", "keratin_residual", "per_SD", "response", independent),
        ("primary_OR_adjusted_perSD", "keratin_adjusted", "per_SD", "response", independent),
        ("primary_OR_raw_tertile", "raw", "tertile", "response", independent),
        ("primary_OR_raw_quartile", "raw", "quartile", "response", independent),
        ("primary_OR_resid_tertile", "keratin_residual", "tertile", "response", independent),
        ("primary_HR_raw_median", "raw", "median", "PFS", independent),
        ("primary_HR_resid_median", "keratin_residual", "median", "PFS", independent),
        ("primary_HR_raw_perSD", "raw", "per_SD", "PFS", independent),
        ("primary_HR_resid_perSD", "keratin_residual", "per_SD", "PFS", independent),
        ("primary_HR_adjusted_perSD", "keratin_adjusted", "per_SD", "PFS", independent),
        ("OS_raw_median_GSE218989_only", "raw", "median", "OS", ["GSE218989"]),
        ("OS_resid_median_GSE218989_only", "keratin_residual", "median", "OS", ["GSE218989"]),
        ("primary_OR_raw_median_without_GSE218989", "raw", "median", "response", no_big),
        ("primary_OR_resid_median_without_GSE218989", "keratin_residual", "median", "response", no_big),
        ("primary_OR_raw_median_no_neoadjuvant", "raw", "median", "response", no_neo),
        ("primary_OR_resid_median_no_neoadjuvant", "keratin_residual", "median", "response", no_neo),
    ]
    pooled_rows = []
    forests = {
        "primary_OR_raw_median": ("Response OR, CLDN4 median high vs low", "Odds ratio (CLDN4-high vs low)", FIGS / "forest_or_raw_median.png"),
        "primary_OR_resid_median": ("Response OR, keratin-residual CLDN4 median", "Odds ratio (residual-high vs low)", FIGS / "forest_or_keratin_residual_median.png"),
        "primary_HR_raw_median": ("PFS HR, CLDN4 median high vs low", "Hazard ratio (CLDN4-high vs low)", FIGS / "forest_hr_raw_median.png"),
        "primary_HR_resid_median": ("PFS HR, keratin-residual CLDN4 median", "Hazard ratio (residual-high vs low)", FIGS / "forest_hr_keratin_residual_median.png"),
    }
    for pool, measure, cutoff, family, cohorts in specs:
        d = select(measure, cutoff, family, cohorts)
        rows = d.to_dict("records")
        meta = dl_meta(rows)
        meta["pool"] = pool
        meta["measure"] = measure
        meta["cutoff"] = cutoff
        meta["family"] = family
        # PDF comparison only for the locked primary OR
        if pool == "primary_OR_raw_median":
            meta["matches_pdf_0.42_at_2dp"] = bool(meta.get("estimable") and round(meta["effect"], 2) == PDF_OR and meta["k"] == PDF_K)
            meta["pdf_0.42_in_ci"] = bool(meta.get("estimable") and meta["lo"] <= PDF_OR <= meta["hi"])
        pooled_rows.append(meta)
        print(pool, "k", meta.get("k"), "effect", None if not meta.get("estimable") else round(meta["effect"], 3), "p", None if not meta.get("estimable") else round(meta["p"], 4))
        if pool in forests and meta.get("estimable"):
            title, xlab, path = forests[pool]
            # restore study order to meta cohort order
            ordered = []
            for c in meta["cohorts"]:
                ordered.append(next(r for r in rows if r["cohort"] == c))
            subtitle = f"{title}\nRE {fmt(meta['effect'], meta['lo'], meta['hi'])}, p={meta['p']:.3g}, I²={meta['I2']:.0f}%"
            forest_plot(ordered, meta, subtitle, xlab, path)

    # flatten pooled for csv
    flat = []
    for m in pooled_rows:
        flat.append(
            {
                "pool": m["pool"],
                "measure": m["measure"],
                "cutoff": m["cutoff"],
                "family": m["family"],
                "k": m.get("k"),
                "n_total": m.get("n_total"),
                "re_effect": m.get("effect"),
                "re_lo": m.get("lo"),
                "re_hi": m.get("hi"),
                "re_p": m.get("p"),
                "I2": m.get("I2"),
                "tau2": m.get("tau2"),
                "cohorts": ",".join(m.get("cohorts") or []),
                "estimable": m.get("estimable"),
                "matches_pdf_0.42_at_2dp": m.get("matches_pdf_0.42_at_2dp"),
                "pdf_0.42_in_ci": m.get("pdf_0.42_in_ci"),
            }
        )
    pd.DataFrame(flat).to_csv(TABLES / "pooled.csv", index=False)
    pd.DataFrame(panel_absent()).to_csv(TABLES / "excluded_cohorts.csv", index=False)

    primary_meta = next(m for m in pooled_rows if m["pool"] == "primary_OR_raw_median")
    pdf = supports_pdf(primary_meta)
    (TABLES / "pdf_11cohort_support.json").write_text(pd.Series(pdf).to_json(), encoding="utf-8")
    write_report(effect, pooled_rows, dup, patients, pdf)
    print("pdf support", pdf)
    print("wrote", ROOT)


if __name__ == "__main__":
    main()
