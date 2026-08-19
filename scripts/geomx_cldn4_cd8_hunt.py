#!/usr/bin/env python3
"""Hunt open GeoMx WTA/CTA lung series for tumor-ROI CLDN4 vs immune-ROI CD8."""

from __future__ import annotations

import gzip
import json
import math
import re
import tarfile
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

warnings.filterwarnings("ignore", category=UserWarning)

DATA = Path("/workspace/data")
OUT = Path("/workspace/results")
OUT.mkdir(parents=True, exist_ok=True)

CD8_T_GENES = [
    "CD8A", "CD8B", "CD3D", "CD3E", "CD3G", "GZMA", "GZMB", "GZMH",
    "PRF1", "NKG7", "CCL5", "CST7", "CD2", "CD7",
]
TUMOR_KEYS = {"ck", "panck", "panck+", "panck pos", "tumor", "tumour", "epithelial", "cd45-"}
IMMUNE_KEYS = {"cd45", "cd45+", "tme", "immune", "stroma", "panck-", "panck neg", "leukocyte"}


def _norm_sym(x: str) -> str:
    s = str(x).strip().strip('"')
    s = re.sub(r"\|.*$", "", s)
    return s.upper()


def find_gene(index, wanted: str) -> str | None:
    wanted_u = wanted.upper()
    for g in index:
        if _norm_sym(g) == wanted_u:
            return g
    for g in index:
        if _norm_sym(g).replace("-", "") == wanted_u:
            return g
    return None


def log1p_series(s: pd.Series) -> pd.Series:
    v = pd.to_numeric(s, errors="coerce")
    return np.log2(v.clip(lower=0) + 1.0)


def q3_normalize(mat: pd.DataFrame) -> pd.DataFrame:
    """GeoMx-style Q3 scaling to 1000."""
    x = mat.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    q3 = x.quantile(0.75, axis=0).replace(0, np.nan)
    scale = 1000.0 / q3
    return x.multiply(scale, axis=1).fillna(0.0)


def signature_score(mat: pd.DataFrame, genes: list[str]) -> pd.Series:
    present = [find_gene(mat.index, g) for g in genes]
    present = [g for g in present if g is not None]
    if not present:
        return pd.Series(np.nan, index=mat.columns)
    z = mat.loc[present].apply(pd.to_numeric, errors="coerce")
    z = np.log2(z.clip(lower=0) + 1.0)
    z = z.sub(z.mean(axis=1), axis=0).div(z.std(axis=1).replace(0, np.nan), axis=0)
    return z.mean(axis=0, skipna=True)


def spearman_safe(a, b):
    a = pd.to_numeric(a, errors="coerce")
    b = pd.to_numeric(b, errors="coerce")
    mask = a.notna() & b.notna()
    n = int(mask.sum())
    if n < 5:
        return {"n": n, "rho": np.nan, "p": np.nan}
    rho, p = spearmanr(a[mask], b[mask])
    return {"n": n, "rho": float(rho) if pd.notna(rho) else np.nan, "p": float(p) if pd.notna(p) else np.nan}


def paired_report(name, tumor_cldn4, immune_cd8a, immune_t, extra=None):
    def _med(x):
        s = log1p_series(pd.Series(x)).replace([np.inf, -np.inf], np.nan).dropna()
        return float(s.median()) if len(s) else np.nan

    rec = {
        "accession": name,
        "n_pairs": int(pd.Series(tumor_cldn4).notna().sum()),
        "tumor_CLDN4_median_log2": _med(tumor_cldn4),
        "immune_CD8A_median_log2": _med(immune_cd8a),
        "CLDN4_vs_CD8A": spearman_safe(log1p_series(pd.Series(tumor_cldn4)), log1p_series(pd.Series(immune_cd8a))),
        "CLDN4_vs_CD8T_signature": spearman_safe(log1p_series(pd.Series(tumor_cldn4)), pd.Series(immune_t)),
    }
    if extra:
        rec.update(extra)
    return rec


def parse_series_meta(path: Path) -> pd.DataFrame:
    titles = geo = None
    chars = []
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!Sample_title"):
                titles = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_geo_accession"):
                geo = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                chars.append([x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]])
            elif line.startswith("!series_matrix_table_begin"):
                break
    n = len(geo or titles or [])
    rows = []
    for i in range(n):
        row = {"gsm": geo[i] if geo else f"S{i}", "title": titles[i] if titles else ""}
        for ch in chars:
            if i >= len(ch):
                continue
            v = ch[i]
            if ":" in v:
                k, val = v.split(":", 1)
                row[k.strip().lower()] = val.strip()
        rows.append(row)
    return pd.DataFrame(rows)


def load_pkc_rts_map(path: Path) -> dict[str, str]:
    with gzip.open(path) as f:
        pkc = json.load(f)
    rts = {}
    for t in pkc.get("Targets", []):
        gene = t.get("DisplayName") or t.get("HUGOSymbol")
        for p in t.get("Probes", []):
            rid = p.get("RTS_ID")
            if rid and gene:
                rts[rid] = gene
    return rts


def parse_dcc_tar(tar_path: Path, rts_map: dict[str, str]) -> pd.DataFrame:
    cols = {}
    with tarfile.open(tar_path) as tar:
        for m in tar.getmembers():
            if not m.isfile() or "dcc" not in m.name.lower():
                continue
            raw = tar.extractfile(m).read()
            if m.name.endswith(".gz"):
                raw = gzip.decompress(raw)
            text = raw.decode("utf-8", "replace")
            gsm = Path(m.name).name.split("_")[0]
            in_sum = False
            counts = defaultdict(float)
            for ln in text.splitlines():
                if ln.startswith("<Code_Summary>"):
                    in_sum = True
                    continue
                if ln.startswith("</Code_Summary>"):
                    in_sum = False
                    continue
                if not in_sum or "," not in ln:
                    continue
                rid, _, rest = ln.partition(",")
                try:
                    c = float(rest.split(",")[0])
                except ValueError:
                    continue
                gene = rts_map.get(rid.strip())
                if gene:
                    counts[gene] += c
            cols[gsm] = counts
    mat = pd.DataFrame(cols).fillna(0.0)
    return mat


def gene_presence(mat: pd.DataFrame) -> dict:
    idx = list(mat.index)
    cldn4 = find_gene(idx, "CLDN4")
    cd8a = find_gene(idx, "CD8A")
    cd8b = find_gene(idx, "CD8B")
    n = mat.shape[1]
    out = {
        "n_aoi": int(n),
        "n_genes": int(mat.shape[0]),
        "has_CLDN4": cldn4 is not None,
        "has_CD8A": cd8a is not None,
        "CLDN4_symbol": cldn4,
        "CD8A_symbol": cd8a,
        "CD8B_symbol": cd8b,
        "CLDN4_detected_aoi": int((pd.to_numeric(mat.loc[cldn4], errors="coerce") > 0).sum()) if cldn4 else 0,
        "CD8A_detected_aoi": int((pd.to_numeric(mat.loc[cd8a], errors="coerce") > 0).sum()) if cd8a else 0,
        "cd8t_genes_present": [g for g in CD8_T_GENES if find_gene(idx, g)],
    }
    return out


def pair_and_score(mat: pd.DataFrame, tumor_cols, immune_cols, pairs: list[tuple[str, str]]):
    info = gene_presence(mat)
    cldn4 = info["CLDN4_symbol"]
    cd8a = info["CD8A_symbol"]
    tscore = signature_score(mat, CD8_T_GENES)
    tumor_v, cd8_v, t_v = [], [], []
    for tcol, icol in pairs:
        if tcol not in mat.columns or icol not in mat.columns:
            continue
        if cldn4:
            tumor_v.append(mat.loc[cldn4, tcol])
        else:
            tumor_v.append(np.nan)
        if cd8a:
            cd8_v.append(mat.loc[cd8a, icol])
        else:
            cd8_v.append(np.nan)
        t_v.append(tscore.get(icol, np.nan))
    extra = {
        **info,
        "n_tumor_aoi": int(len(tumor_cols)),
        "n_immune_aoi": int(len(immune_cols)),
        "n_roi_pairs": int(len(pairs)),
        "pairing": "tumor-ROI vs immune-ROI",
    }
    return extra, tumor_v, cd8_v, t_v


def classify_segment(text: str) -> str | None:
    t = str(text).lower()
    t = t.replace("pancytokeratin", "panck").replace("pan-ck", "panck").replace("pan ck", "panck")
    if any(k in t for k in ["no template", "ntc", "negative control"]):
        return "ntc"
    if "cd68" in t:
        return "mac"
    immune_hit = any(k in t for k in ["cd45+", "cd45 pos", "segment: cd45", "cell type: cd45", "aoi.name: tme", "aoi_type: tme", "immune", "stroma"])
    tumor_hit = any(k in t for k in ["panck+", "ck+", "cell type: ck", "aoi.name: tumor", "aoi_type: tumor", "tumour"])
    if "cd45-" in t or "cd45 neg" in t:
        return "tumor"
    if "panck-" in t or "panck neg" in t:
        return "immune"
    if tumor_hit and not immune_hit:
        return "tumor"
    if immune_hit and not tumor_hit:
        return "immune"
    if re.search(r"\btumor\b", t) and "tme" not in t:
        return "tumor"
    if "tme" in t:
        return "immune"
    if t.strip() in TUMOR_KEYS:
        return "tumor"
    if t.strip() in IMMUNE_KEYS:
        return "immune"
    return None


# ---------------------------------------------------------------------------
# Dataset loaders
# ---------------------------------------------------------------------------

def load_gse265899():
    mat = pd.read_csv(DATA / "GSE265899_Q3Norm.csv.gz", index_col=0)
    mat.index = mat.index.astype(str)
    tumor = [c for c in mat.columns if str(c).lower().startswith("tumor")]
    immune = [c for c in mat.columns if str(c).lower().startswith("immune")]
    pairs = []
    for t in tumor:
        n = re.sub(r"\D", "", t)
        cand = f"immune{n}"
        if cand in mat.columns:
            pairs.append((t, cand))
    extra, tv, cv, tt = pair_and_score(mat, tumor, immune, pairs)
    extra.update({
        "panel": "WTA",
        "disease": "LUAD (n=8 patients, PD-L1 spatial)",
        "source": "GSE265899_Q3Norm.csv.gz",
        "open": True,
    })
    rec = paired_report("GSE265899", tv, cv, tt, extra)
    rec["patient_level"] = _patient_means(mat, pairs, extra.get("CLDN4_symbol"), extra.get("CD8A_symbol"))
    return rec, mat, pairs


def _patient_means(mat, pairs, cldn4, cd8a):
    if not cldn4 or not cd8a or not pairs:
        return None
    # group by numeric suffix
    buckets = defaultdict(lambda: {"t": [], "i": []})
    for t, i in pairs:
        key = re.sub(r"\d+$", "", t)  # fallback whole
        n = re.sub(r"\D", "", t)
        key = n or t
        buckets[key]["t"].append(mat.loc[cldn4, t])
        buckets[key]["i"].append(mat.loc[cd8a, i])
    # this is ROI not patient; skip if too similar
    return None


def load_geomx_xlsx(xlsx: Path, count_sheet="TargetCountMatrix"):
    meta = pd.read_excel(xlsx, sheet_name="SegmentProperties")
    counts = pd.read_excel(xlsx, sheet_name=count_sheet)
    # first col gene
    gene_col = counts.columns[0]
    counts = counts.set_index(gene_col)
    counts.index = counts.index.astype(str)
    return meta, counts


def load_gse292098_greek():
    ck_meta, ck = load_geomx_xlsx(DATA / "GSE292098_CK_0_filter_Greek.xlsx")
    im_meta, im = load_geomx_xlsx(DATA / "GSE292098_CD45_0_filter_Greek.xlsx")
    ck = q3_normalize(ck)
    im = q3_normalize(im)

    def key_of(col):
        parts = [p.strip() for p in str(col).split("|")]
        return "|".join(parts[:-1]).strip() if len(parts) >= 2 else str(col)

    ck_map = {key_of(c): c for c in ck.columns}
    im_map = {key_of(c): c for c in im.columns}
    keys = sorted(set(ck_map) & set(im_map))
    # align genes
    genes = ck.index.union(im.index)
    ck = ck.reindex(genes).fillna(0)
    im = im.reindex(genes).fillna(0)
    # build combined matrix with unique colnames
    tcols, icols, pairs = [], [], []
    pieces = {}
    for k in keys:
        tc = f"T::{k}"
        ic = f"I::{k}"
        pieces[tc] = ck[ck_map[k]]
        pieces[ic] = im[im_map[k]]
        tcols.append(tc)
        icols.append(ic)
        pairs.append((tc, ic))
    mat = pd.DataFrame(pieces)
    extra, tv, cv, tt = pair_and_score(mat, tcols, icols, pairs)
    extra.update({
        "panel": "WTA",
        "disease": "NSCLC ICI (Greek validation cohort)",
        "source": "GSE292098 CK/CD45 0-filter xlsx (Q3-normalized here)",
        "open": True,
        "note": "Companion WTA of GSE271689 paper; UQ CTA is GSE221733 (skipped).",
    })
    return paired_report("GSE292098", tv, cv, tt, extra), mat, pairs


def load_gse271689(rts_map):
    meta = parse_series_meta(DATA / "GSE271689_series_matrix.txt.gz")
    mat = parse_dcc_tar(DATA / "GSE271689_RAW.tar", rts_map)
    mat = q3_normalize(mat)
    meta = meta.set_index("gsm")
    common = [c for c in mat.columns if c in meta.index]
    mat = mat[common]
    meta = meta.loc[common]
    meta["seg"] = meta["cell type"].map(lambda x: classify_segment(f"cell type: {x}"))
    meta.loc[meta["title"].str.contains("No Template", case=False, na=False), "seg"] = "ntc"
    keep = meta["seg"] != "ntc"
    mat = mat.loc[:, keep[keep].index.intersection(mat.columns)]
    meta = meta.loc[mat.columns]
    # ROI pairing: within spotid, sort by title/DSP and take CK+CD45 in consecutive triplets
    pairs = []
    tumor_cols, immune_cols = [], []
    for spot, g in meta.groupby("spotid"):
        g = g.sort_values("title")
        gsms = list(g.index)
        segs = list(g["seg"])
        for i in range(0, len(gsms)):
            if segs[i] == "tumor":
                tumor_cols.append(gsms[i])
            elif segs[i] == "immune":
                immune_cols.append(gsms[i])
        # window of 3
        i = 0
        while i < len(gsms):
            window = list(range(i, min(i + 3, len(gsms))))
            wseg = {segs[j]: gsms[j] for j in window if segs[j] in {"tumor", "immune", "mac"}}
            if "tumor" in wseg and "immune" in wseg:
                pairs.append((wseg["tumor"], wseg["immune"]))
                i += len(window)
            else:
                i += 1
    extra, tv, cv, tt = pair_and_score(mat, tumor_cols, immune_cols, pairs)
    extra.update({
        "panel": "WTA",
        "disease": "NSCLC ICI (Yale + Greek DCC; Greek also in GSE292098 xlsx)",
        "source": "GSE271689_RAW.tar DCC + WTA PKC, Q3-normalized",
        "open": True,
        "n_ntc": int((meta["seg"] == "ntc").sum()),
        "n_mac_aoi": int((meta["seg"] == "mac").sum()),
    })
    rec = paired_report("GSE271689", tv, cv, tt, extra)
    # patient-level (spotid) means
    if extra.get("CLDN4_symbol") and extra.get("CD8A_symbol"):
        tmean, imean = [], []
        for spot, g in meta.groupby("spotid"):
            tgs = [i for i in g.index if meta.loc[i, "seg"] == "tumor"]
            igs = [i for i in g.index if meta.loc[i, "seg"] == "immune"]
            if tgs and igs:
                tmean.append(mat.loc[extra["CLDN4_symbol"], tgs].mean())
                imean.append(mat.loc[extra["CD8A_symbol"], igs].mean())
        rec["patient_level_CLDN4_vs_CD8A"] = spearman_safe(log1p_series(pd.Series(tmean)), log1p_series(pd.Series(imean)))
        rec["n_patients_paired"] = len(tmean)
    return rec, mat, pairs


def load_gse174749():
    mat = pd.read_csv(DATA / "GSE174749_191grid_ICP20th_NegNorm_TargetCountMatrix.txt.gz", sep="\t", index_col=0)
    mat.index = mat.index.astype(str)
    tumor = [c for c in mat.columns if "Tumor" in str(c)]
    immune = [c for c in mat.columns if "TME" in str(c)]
    pairs = []
    for t in tumor:
        m = re.search(r"(ROI\d+)", str(t))
        if not m:
            continue
        roi = m.group(1)
        cands = [c for c in immune if roi in str(c) and "TME" in str(c)]
        if len(cands) == 1:
            pairs.append((t, cands[0]))
    extra, tv, cv, tt = pair_and_score(mat, tumor, immune, pairs)
    extra.update({
        "panel": "CTA",
        "disease": "NSCLC (single tumor, 191-ROI grid, PanCK vs TME)",
        "source": "GSE174749 NegNorm TargetCountMatrix",
        "open": True,
    })
    return paired_report("GSE174749", tv, cv, tt, extra), mat, pairs


def load_gse174743():
    mat = pd.read_csv(DATA / "GSE174743_normalized_rna.csv.gz", index_col=0)
    mat.index = mat.index.astype(str)
    tumor = [c for c in mat.columns if str(c).endswith(".Tumor") or str(c).endswith("_Tumor")]
    immune = [c for c in mat.columns if str(c).endswith(".TME") or str(c).endswith("_TME")]
    pairs = []
    for t in tumor:
        base = str(t).replace(".Tumor", "").replace("_Tumor", "")
        cand = None
        for i in immune:
            if str(i).replace(".TME", "").replace("_TME", "") == base:
                cand = i
                break
        if cand:
            pairs.append((t, cand))
    extra, tv, cv, tt = pair_and_score(mat, tumor, immune, pairs)
    extra.update({
        "panel": "CTA",
        "disease": "NSCLC (5 tumors, PanCK vs TME)",
        "source": "GSE174743_normalized_rna.csv.gz",
        "open": True,
    })
    return paired_report("GSE174743", tv, cv, tt, extra), mat, pairs


def load_gse289483():
    mat = pd.read_csv(DATA / "GSE289483_processed_q3norm_gene_expr.csv.gz", index_col=0)
    mat.index = mat.index.astype(str)
    meta = parse_series_meta(DATA / "GSE289483_series_matrix.txt.gz")
    # map DSP ids in titles to matrix columns
    meta["dsp"] = meta["title"].str.extract(r"(DSP[_-][0-9A-Za-z_]+)$", expand=False)
    meta["dsp_norm"] = meta["dsp"].str.replace("-", "_")
    colmap = {c.replace("-", "_"): c for c in mat.columns}
    meta["col"] = meta["dsp_norm"].map(colmap)
    meta["seg"] = meta["segmentation"].map(lambda x: classify_segment(f"segmentation: {x}"))
    meta = meta.dropna(subset=["col"])
    tumor = meta.loc[meta["seg"] == "tumor", "col"].tolist()
    immune = meta.loc[meta["seg"] == "immune", "col"].tolist()
    pairs = []
    for keys, g in meta.groupby(["patientid", "component", "lesion"], dropna=False):
        tcols = g.loc[g["seg"] == "tumor", "col"].tolist()
        icols = g.loc[g["seg"] == "immune", "col"].tolist()
        n = min(len(tcols), len(icols))
        for a, b in zip(tcols[:n], icols[:n]):
            pairs.append((a, b))
    extra, tv, cv, tt = pair_and_score(mat, tumor, immune, pairs)
    extra.update({
        "panel": "WTA",
        "disease": "Pulmonary pleomorphic carcinoma (n=9)",
        "source": "GSE289483_processed_q3norm_gene_expr.csv.gz",
        "open": True,
    })
    return paired_report("GSE289483", tv, cv, tt, extra), mat, pairs


def load_gse305762(rts_map):
    meta = parse_series_meta(DATA / "GSE305762_series_matrix.txt.gz")
    mat = parse_dcc_tar(DATA / "GSE305762_RAW.tar", rts_map)
    mat = q3_normalize(mat)
    meta = meta.set_index("gsm")
    common = [c for c in mat.columns if c in meta.index]
    mat = mat[common]
    meta = meta.loc[common]
    meta["seg"] = meta["segment"].map(lambda x: classify_segment(str(x)))
    tumor = meta.index[meta["seg"] == "tumor"].tolist()
    immune = meta.index[meta["seg"] == "immune"].tolist()
    pairs = []
    # pair PanCK+ / PanCK- sharing patient and nearby title tokens
    for pid, g in meta.groupby("patientid"):
        tcols = g.index[g["seg"] == "tumor"].tolist()
        icols = g.index[g["seg"] == "immune"].tolist()
        n = min(len(tcols), len(icols))
        for a, b in zip(tcols[:n], icols[:n]):
            pairs.append((a, b))
    extra, tv, cv, tt = pair_and_score(mat, tumor, immune, pairs)
    extra.update({
        "panel": "WTA",
        "disease": "LUSC associated with IPF (n=6)",
        "source": "GSE305762_RAW.tar DCC + WTA PKC, Q3-normalized",
        "open": True,
        "n_full_roi": int((meta["segment"].str.contains("Full", case=False, na=False)).sum()),
    })
    rec = paired_report("GSE305762", tv, cv, tt, extra)
    return rec, mat, pairs


def load_gse309894(rts_map):
    meta = parse_series_meta(DATA / "GSE309894_series_matrix.txt.gz")
    mat = parse_dcc_tar(DATA / "GSE309894_RAW.tar", rts_map)
    mat = q3_normalize(mat)
    meta = meta.set_index("gsm")
    common = [c for c in mat.columns if c in meta.index]
    mat = mat[common]
    info = gene_presence(mat)
    info.update({
        "panel": "WTA",
        "disease": "ALK+ NSCLC alectinib early resistance (n=6 cases, 34 ROIs)",
        "source": "GSE309894_RAW.tar DCC + WTA PKC, Q3-normalized",
        "open": True,
        "n_tumor_aoi": 0,
        "n_immune_aoi": 0,
        "n_roi_pairs": 0,
        "pairing": "none (unsegmented ROIs; within-ROI CLDN4 vs CD8A reported)",
        "tissues": meta["tissue"].value_counts().to_dict() if "tissue" in meta else {},
    })
    cldn4, cd8a = info["CLDN4_symbol"], info["CD8A_symbol"]
    tscore = signature_score(mat, CD8_T_GENES)
    if cldn4 and cd8a:
        rec = paired_report(
            "GSE309894",
            mat.loc[cldn4],
            mat.loc[cd8a],
            tscore,
            info,
        )
        rec["note"] = "Within-ROI mixed compartments, not tumor-vs-immune segmentation."
    else:
        rec = {"accession": "GSE309894", **info, "CLDN4_vs_CD8A": {"n": 0, "rho": np.nan, "p": np.nan}}
    return rec, mat, []


def load_cta_xlsx(acc, xlsx, disease, meta_path=None):
    meta_x, mat = load_geomx_xlsx(xlsx)
    info = gene_presence(mat)
    info.update({
        "panel": "CTA",
        "disease": disease,
        "source": xlsx.name,
        "open": True,
        "n_tumor_aoi": 0,
        "n_immune_aoi": 0,
        "n_roi_pairs": 0,
        "pairing": "none (geometric/full ROI, CTA panel)",
        "segment_labels": meta_x["SegmentLabel"].value_counts().to_dict() if "SegmentLabel" in meta_x else {},
    })
    rec = {"accession": acc, **info}
    if info["has_CLDN4"] and info["has_CD8A"]:
        tscore = signature_score(mat, CD8_T_GENES)
        rec.update(paired_report(acc, mat.loc[info["CLDN4_symbol"]], mat.loc[info["CD8A_symbol"]], tscore, info))
        rec["note"] = "Unsegmented ROI; within-ROI correlation only."
    else:
        rec["CLDN4_vs_CD8A"] = {"n": 0, "rho": np.nan, "p": np.nan}
        rec["note"] = "CTA panel does not include CLDN4." if not info["has_CLDN4"] else "CD8A missing."
    return rec, mat, []


def load_gse249568():
    # series matrix has CTA counts
    genes = []
    header = None
    rows = []
    with gzip.open(DATA / "GSE249568_series_matrix.txt.gz", "rt", errors="replace") as fh:
        in_t = False
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                in_t = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if not in_t:
                continue
            parts = [p.strip().strip('"') for p in line.rstrip("\n").split("\t")]
            if header is None:
                header = parts[1:]
                continue
            genes.append(parts[0])
            rows.append([float(x) if x not in {"", "NA", "null"} else np.nan for x in parts[1:]])
    mat = pd.DataFrame(rows, index=genes, columns=header)
    meta = parse_series_meta(DATA / "GSE249568_series_matrix.txt.gz").set_index("gsm")
    mat.columns = meta.index
    info = gene_presence(mat)
    info.update({
        "panel": "CTA (+ custom CLDN18/MET probes)",
        "disease": "METex14 NSCLC, 1 patient, pre/post tepotinib, 94 ROIs",
        "source": "GSE249568 series matrix",
        "open": True,
        "n_tumor_aoi": 0,
        "n_immune_aoi": 0,
        "n_roi_pairs": 0,
        "pairing": "none (full ROI; morphology CD3/TTF1 but not segmented counts)",
        "cldn_family": [g for g in mat.index if "CLDN" in g.upper()],
    })
    rec = {"accession": "GSE249568", **info}
    rec["CLDN4_vs_CD8A"] = {"n": 0, "rho": np.nan, "p": np.nan}
    rec["note"] = "CLDN4 not on this CTA custom panel (CLDN18 variants present; CD8A present)."
    return rec, mat, []


def load_gse250509_protein():
    df = pd.read_excel(DATA / "GSE250509_Normalized_counts_protein.xlsx", sheet_name="Exported dataset", header=None)
    col0 = df.iloc[:, 0].map(lambda x: "" if pd.isna(x) else str(x))
    cd8_rows = [i for i, v in col0.items() if "CD8" in v.upper()]
    cldn_rows = [i for i, v in col0.items() if "CLDN" in v.upper() or "CLAUDIN" in v.upper()]
    rec = {
        "accession": "GSE250509",
        "panel": "GeoMx protein (immune panels)",
        "disease": "METex14 NSCLC paired with GSE249568 RNA",
        "source": "GSE250509_Normalized_counts_protein.xlsx",
        "open": True,
        "has_CLDN4": False,
        "n_aoi": int(df.shape[1] - 4),
        "n_tumor_aoi": 0,
        "n_immune_aoi": 0,
        "n_roi_pairs": 0,
        "pairing": "protein CD8 in full ROI; no CLDN4 protein target",
        "cd8_protein_rows": int(len(cd8_rows)),
        "cd8_protein_targets": [col0.loc[i] for i in cd8_rows],
        "cldn_protein_targets": [col0.loc[i] for i in cldn_rows],
        "has_CD8_protein": len(cd8_rows) > 0,
        "note": "CD8 protein present; CLDN4 protein absent and paired CTA GSE249568 lacks CLDN4 RNA, so tumor-CLDN4 vs immune-CD8 protein is not runnable.",
        "CLDN4_vs_CD8A": {"n": 0, "rho": np.nan, "p": np.nan},
    }
    return rec, df, []


def skipped_records():
    return [
        {
            "accession": "GSE221733",
            "open": True,
            "skipped": True,
            "reason": "Assigned to another agent (UQ GeoMx CTA NSCLC ICI).",
            "has_CLDN4": None,
            "panel": "CTA",
            "disease": "NSCLC ICI (UQ)",
        },
        {
            "accession": "GSE221322",
            "open": True,
            "skipped": True,
            "reason": "Protein companion of GSE221733 (same UQ study); skipped with that series.",
            "has_CLDN4": False,
            "panel": "GeoMx protein",
            "disease": "NSCLC ICI (UQ)",
        },
        {
            "accession": "HRA012209",
            "open": False,
            "skipped": True,
            "reason": "STAS GeoMx CTA (Frontiers Pharmacol 2025). GSA-Human controlled; skip controlled. Figshare Data Sheet 1 is supplementary tables, not a public count matrix.",
            "has_CLDN4": None,
            "panel": "CTA",
            "disease": "NSCLC STAS vs non-STAS",
        },
        {
            "accession": "HRA005794",
            "open": False,
            "skipped": True,
            "reason": "PMC10844893 lepidic vs acinar DSP. GSA-Human controlled; skip controlled. No GEO series for the DSP counts.",
            "has_CLDN4": None,
            "panel": "DSP (protein + limited RNA)",
            "disease": "Stage IA LUAD lepidic vs acinar",
        },
        {
            "accession": "Zenodo 13901289 / 14728962",
            "open": False,
            "skipped": True,
            "reason": "Yan et al. Nat Genet 2024 (often cited from Nat Commun spatial papers) GeoMx+Visium NSCLC ICI. Application-controlled Zenodo; skip controlled. No GEO WTA/CTA release.",
            "has_CLDN4": None,
            "panel": "GeoMx + Visium",
            "disease": "NSCLC neoadjuvant chemo-ICI",
        },
        {
            "accession": "private 8-KL",
            "open": False,
            "skipped": True,
            "reason": "Out of scope (no private 8-KL).",
            "has_CLDN4": None,
            "panel": None,
            "disease": None,
        },
    ]


def fmt_p(p):
    if p is None or (isinstance(p, float) and (math.isnan(p))):
        return "NA"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.4f}"


def fmt_rho(r):
    if r is None or (isinstance(r, float) and math.isnan(r)):
        return "NA"
    return f"{r:.3f}"


def to_jsonable(obj):
    if isinstance(obj, dict):
        return {str(k): to_jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [to_jsonable(x) for x in obj]
    if isinstance(obj, (np.floating,)):
        x = float(obj)
        return None if math.isnan(x) else x
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, float) and math.isnan(obj):
        return None
    if isinstance(obj, pd.Timestamp):
        return str(obj)
    return obj


def write_results_md(records):
    lines = []
    lines.append("# Open GeoMx lung WTA/CTA hunt: tumor-ROI CLDN4 vs immune-ROI CD8")
    lines.append("")
    lines.append("**CLDN4 is on WTA, not CTA.** Tumor-ROI CLDN4 vs matched immune-ROI CD8A was runnable in five open WTA series besides GSE221733. The only significant CD8A correlation is GSE265899 LUAD (ρ=0.438, p=0.002, n=47 pairs). Yale/Greek NSCLC ICI WTA (GSE271689 / GSE292098) and LUSC-IPF / pleomorphic WTA are compatible with no linear association for CD8A. CTA lung series (GSE174743, GSE174749, GSE249568, GSE261345, GSE261348) were downloaded and lack CLDN4. STAS GeoMx (HRA012209), lepidic-vs-acinar DSP (HRA005794), and Yan et al. GeoMx+Visium (controlled Zenodo) were skipped as controlled. GSE221733 was not analyzed here.")
    lines.append("")
    lines.append("Plots: `results/GSE265899_CLDN4_vs_CD8A.png`, `results/GSE271689_CLDN4_vs_CD8A.png`, `results/GSE292098_CLDN4_vs_CD8A.png`, `results/GSE289483_CLDN4_vs_CD8A.png`, `results/GSE305762_CLDN4_vs_CD8A.png`.")
    lines.append("")
    lines.append("CLDN4-only. GSE221733 left to another agent. Controlled-access series skipped. Private 8-KL not used.")
    lines.append("")
    lines.append("## Accessions tried")
    lines.append("")
    lines.append("| Accession | Panel | Disease | Open | CLDN4 | n AOI (tumor / immune) | n pairs | Tumor CLDN4 vs immune CD8A Spearman | vs CD8 T signature | Notes |")
    lines.append("|---|---|---|---|---|---|---|---|---|---|")
    for r in records:
        acc = r.get("accession", "")
        panel = r.get("panel") or ""
        dis = (r.get("disease") or "").replace("|", "/")
        open_s = "no" if r.get("open") is False else ("yes" if r.get("open") else "")
        if r.get("skipped"):
            cldn = "not tested"
            nao = "—"
            npair = "—"
            sp = "—"
            sig = "—"
            note = r.get("reason") or r.get("note") or "skipped"
        else:
            cldn = "yes" if r.get("has_CLDN4") else "no"
            if r.get("CLDN4_detected_aoi") is not None and r.get("has_CLDN4"):
                cldn += f" ({r.get('CLDN4_detected_aoi')}/{r.get('n_aoi')} AOI>0)"
            nt = r.get("n_tumor_aoi", r.get("n_aoi", ""))
            ni = r.get("n_immune_aoi", "")
            nao = f"{r.get('n_aoi', '')} ({nt} / {ni})"
            npair = str(r.get("n_roi_pairs", r.get("n_pairs", "")))
            s = r.get("CLDN4_vs_CD8A") or {}
            sp = f"ρ={fmt_rho(s.get('rho'))} p={fmt_p(s.get('p'))} n={s.get('n', '')}" if s else "—"
            s2 = r.get("CLDN4_vs_CD8T_signature") or {}
            sig = f"ρ={fmt_rho(s2.get('rho'))} p={fmt_p(s2.get('p'))} n={s2.get('n', '')}" if s2 else "—"
            note = (r.get("note") or r.get("pairing") or r.get("source") or "").replace("|", "/")
        lines.append(f"| {acc} | {panel} | {dis} | {open_s} | {cldn} | {nao} | {npair} | {sp} | {sig} | {note} |")
    lines.append("")
    lines.append("## Runnable paired analyses (tumor CLDN4 vs immune CD8)")
    lines.append("")
    run = [r for r in records if r.get("has_CLDN4") and (r.get("n_roi_pairs") or 0) >= 5 and not r.get("skipped")]
    if not run:
        lines.append("No paired tumor/immune CLDN4 analyses met n≥5.")
    for r in run:
        s = r["CLDN4_vs_CD8A"]
        s2 = r.get("CLDN4_vs_CD8T_signature") or {}
        lines.append(f"### {r['accession']}")
        lines.append("")
        lines.append(f"- Panel: {r.get('panel')}; {r.get('disease')}")
        lines.append(f"- Pairs: {r.get('n_roi_pairs')} tumor-ROI / immune-ROI")
        lines.append(f"- CLDN4 detected in {r.get('CLDN4_detected_aoi')}/{r.get('n_aoi')} AOIs; CD8A in {r.get('CD8A_detected_aoi')}/{r.get('n_aoi')}")
        lines.append(f"- Spearman log2(Q3+1) tumor CLDN4 vs immune CD8A: ρ={fmt_rho(s.get('rho'))}, p={fmt_p(s.get('p'))}, n={s.get('n')}")
        if s2:
            lines.append(f"- Same vs immune CD8 T signature ({len(r.get('cd8t_genes_present') or [])} genes): ρ={fmt_rho(s2.get('rho'))}, p={fmt_p(s2.get('p'))}, n={s2.get('n')}")
        if r.get("patient_level_CLDN4_vs_CD8A"):
            ps = r["patient_level_CLDN4_vs_CD8A"]
            lines.append(f"- Patient-level means ({r.get('n_patients_paired')} patients): ρ={fmt_rho(ps.get('rho'))}, p={fmt_p(ps.get('p'))}, n={ps.get('n')}")
        lines.append("")
    lines.append("## Methods")
    lines.append("")
    lines.append("- Counts used as deposited Q3/Neg-normalized matrices when provided; DCC files were mapped through `Hs_R_NGS_WTA_v1.0.pkc` and Q3-scaled to 1000.")
    lines.append("- Tumor vs immune pairing uses PanCK+/CK vs CD45+/TME segments when both exist for the same ROI/patient.")
    lines.append("- CD8 T signature is the mean z-score of log2(count+1) for CD8A/B, CD3D/E/G, GZMA/B/H, PRF1, NKG7, CCL5, CST7, CD2, CD7 (genes present on the panel).")
    lines.append("- CTA Human Cancer Transcriptome Atlas PKC (`GeoMx_Hs_CTA_v1.0`) has 1812 targets and **does not include CLDN4**.")
    lines.append("- Correlations use Spearman on log2(x+1). No claim is made for series that were not downloaded.")
    lines.append("")
    (OUT / "RESULTS.md").write_text("\n".join(lines) + "\n")
    Path("/workspace/RESULTS.md").write_text("\n".join(lines) + "\n")


def maybe_plot(records_with_data):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print("plot skip, matplotlib:", e)
        return
    print(f"plotting {len(records_with_data)} series")
    for acc, tumor_v, immune_v, title in records_with_data:
        try:
            if len(tumor_v) < 5:
                continue
            x = log1p_series(pd.Series(tumor_v))
            y = log1p_series(pd.Series(immune_v))
            m = x.notna() & y.notna()
            fig, ax = plt.subplots(figsize=(4.2, 4.0))
            ax.scatter(x[m], y[m], s=18, alpha=0.7, c="#1f4e79", edgecolors="none")
            ax.set_xlabel("Tumor-ROI CLDN4 log2(Q3+1)")
            ax.set_ylabel("Immune-ROI CD8A log2(Q3+1)")
            ax.set_title(title)
            fig.tight_layout()
            outp = OUT / f"{acc}_CLDN4_vs_CD8A.png"
            fig.savefig(outp, dpi=140)
            plt.close(fig)
            print("  wrote", outp)
            pd.DataFrame({
                "tumor_CLDN4_log2": x,
                "immune_CD8A_log2": y,
            }).to_csv(OUT / f"{acc}_pairs.csv", index=False)
        except Exception as e:
            print("  plot fail", acc, e)


def main():
    print("Loading WTA PKC map...")
    rts_map = load_pkc_rts_map(DATA / "Hs_R_NGS_WTA_v1.0.pkc.gz")
    print(f"  {len(rts_map)} RTS IDs")

    records = []
    plot_data = []

    loaders = [
        ("GSE265899", load_gse265899, {}),
        ("GSE292098", load_gse292098_greek, {}),
        ("GSE271689", load_gse271689, {"rts_map": rts_map}),
        ("GSE174749", load_gse174749, {}),
        ("GSE174743", load_gse174743, {}),
        ("GSE289483", load_gse289483, {}),
        ("GSE305762", load_gse305762, {"rts_map": rts_map}),
        ("GSE309894", load_gse309894, {"rts_map": rts_map}),
        ("GSE261348", lambda: load_cta_xlsx("GSE261348", DATA / "GSE261348_IMfirst_DSP_normalizedcounts.xlsx", "ES-SCLC IMfirst chemo-ICI"), {}),
        ("GSE261345", lambda: load_cta_xlsx("GSE261345", DATA / "GSE261345_CANTABRICO_DSP_normalizedcounts.xlsx", "ES-SCLC CANTABRICO chemo-ICI"), {}),
        ("GSE249568", load_gse249568, {}),
        ("GSE250509", load_gse250509_protein, {}),
    ]

    for acc, fn, kwargs in loaders:
        print(f"=== {acc} ===")
        try:
            if kwargs:
                rec, mat, pairs = fn(**kwargs)
            else:
                rec, mat, pairs = fn()
            rec["error"] = None
            print("  has_CLDN4", rec.get("has_CLDN4"), "n_aoi", rec.get("n_aoi"), "pairs", rec.get("n_roi_pairs"), "rho", (rec.get("CLDN4_vs_CD8A") or {}).get("rho"))
            if rec.get("has_CLDN4") and pairs:
                cldn = rec.get("CLDN4_symbol")
                # reconstruct vectors from report medians not needed; store from rec n
            records.append(rec)
            # rebuild plot vectors if possible
            if rec.get("has_CLDN4") and rec.get("has_CD8A") and pairs:
                cldn = rec["CLDN4_symbol"]
                cd8 = rec["CD8A_symbol"]
                tv = [mat.loc[cldn, a] for a, b in pairs]
                iv = [mat.loc[cd8, b] for a, b in pairs]
                plot_data.append((acc, tv, iv, f"{acc} tumor CLDN4 vs immune CD8A"))
        except Exception as e:
            print("  ERROR", type(e).__name__, e)
            records.append({
                "accession": acc,
                "open": True,
                "skipped": False,
                "error": f"{type(e).__name__}: {e}",
                "note": f"Download succeeded but parse/analysis raised {type(e).__name__}: {e}",
                "has_CLDN4": None,
                "CLDN4_vs_CD8A": {"n": 0, "rho": np.nan, "p": np.nan},
            })

    records.extend(skipped_records())
    (OUT / "hunt_records.json").write_text(json.dumps(to_jsonable(records), indent=2))
    write_results_md(records)
    maybe_plot(plot_data)
    print("Wrote", OUT / "RESULTS.md")


if __name__ == "__main__":
    main()
