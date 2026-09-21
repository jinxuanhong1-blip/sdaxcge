#!/usr/bin/env python3
"""Open HRA/CNGB Chinese human lung scRNA: malignant CLDN4 vs STING/IFN/APM.

Catalog (no DAC requests):
  - GSA-Human finished.json: lung single-cell / single-nucleus titles, Open vs Controlled.
  - OMIX release list: open Homo sapiens lung single-cell processed tables.
  - CNGB: pages that were actually fetched (no invented access flags).

Analysis:
  - HRA009335 is the only Open human lung *tumor* snRNA accession. Public files are
    PE150 FASTQ. A read peek is recorded; per-nucleus CLDN4 is not scored unless a
    cell-barcode read is present.
  - OMIX002441 is the Open processed count matrix for the same BioProject as
    Controlled HRA001232 (SCLC, not LUAD). Malignant cells are the author label
    cell_type == tumor. Within each patient, CLDN4-positive vs CLDN4-zero malignant
    cells are compared on STING, IFN, and MHC-I/APM scores. Not merged with
    concordant-4 LUAD.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
from collections import Counter

import numpy as np
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "hra_cngb_lung_cldn4_sting")

LUNG_RE = re.compile(
    r"\b(lung|nsclc|sclc|luad|lusc|pulmonary|bronch|alveolar|pneumo|pleura)\b",
    re.I,
)
SC_RE = re.compile(
    r"(single[- ]cell|single[- ]nucleus|snRNA|scRNA|sc-RNA|10x|10X|"
    r"chromium|drop-seq|smart-seq|spatial transcript|visium|stereo-seq|"
    r"xenium|cosmx)",
    re.I,
)

# Tight programs. CLDN4 is the split, so it is not inside any score.
STING = ["CGAS", "STING1", "TMEM173", "TBK1", "IKBKE", "IRF3", "IRF7"]
IFN = [
    "STAT1", "STAT2", "IRF9", "IRF1", "MX1", "MX2", "OAS1", "OAS2", "OAS3",
    "OASL", "IFIT1", "IFIT2", "IFIT3", "ISG15", "ISG20", "IFI6", "IFI27",
    "IFI35", "IFI44", "IFI44L", "IFIH1", "DDX58", "RSAD2", "EIF2AK2", "BST2",
    "XAF1", "GBP1", "USP18", "HERC5", "IFITM1",
]
APM = [
    "HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "TAPBP", "PSMB8",
    "PSMB9", "PSMB10", "PSME1", "PSME2", "CALR", "PDIA3", "NLRC5", "ERAP1",
    "ERAP2", "HLA-E",
]
# Junction / epithelial control. CLDN4 held out.
TJ = ["CLDN3", "CLDN7", "TJP1", "OCLN", "CDH1", "EPCAM", "TACSTD2"]
KEY_GENES = [
    "CLDN4", "STING1", "CGAS", "TBK1", "IRF3", "STAT1", "ISG15", "MX1",
    "HLA-A", "B2M", "TAP1", "EPCAM", "PTPRC", "CD3D", "CD68",
]
FAMILIES = {"STING": STING, "IFN": IFN, "APM": APM, "TJ": TJ}

MIN_MALIGNANT = 40
MIN_GROUP = 15


def write_tsv(path: str, rows: list[dict], fields: list[str]) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def build_hra_catalog(finished_json: str) -> list[dict]:
    data = json.load(open(finished_json, encoding="utf-8"))
    rows = []
    for x in data:
        title = x.get("title") or ""
        if not LUNG_RE.search(title):
            continue
        if not SC_RE.search(title):
            continue
        controlled = str(x.get("isControlledAccess")) == "1"
        rows.append(
            {
                "accession": x.get("accession") or "",
                "bioproject": x.get("bioproject") or "",
                "access": "Controlled" if controlled else "Open",
                "is_controlled": "1" if controlled else "0",
                "title": title.replace("\t", " ").replace("\n", " "),
                "organization": (x.get("uorganization") or x.get("organization") or "").replace("\t", " "),
                "release_time_string": x.get("lastTimeString") or "",
                "study_type_id": x.get("studyTypeId") if x.get("studyTypeId") is not None else "",
                "url": f"https://ngdc.cncb.ac.cn/gsa-human/browse/{x.get('accession')}",
            }
        )
    rows.sort(key=lambda r: (r["access"], r["accession"]))
    return rows


def _strip(s: str) -> str:
    return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", "", s)).strip()


def build_omix_catalog(omix_html: str) -> list[dict]:
    raw = open(omix_html, encoding="utf-8", errors="replace").read()
    row_re = re.compile(
        r"<tr>\s*<td>(OMIX\d+)</td>\s*<td>(.*?)</td>\s*<td[^>]*>(.*?)</td>"
        r"\s*<td[^>]*>(.*?)</td>\s*<td>(.*?)</td>\s*<td>(.*?)</td>",
        re.S,
    )
    sc = re.compile(
        r"(single[- ]cell|scRNA|snRNA|single[- ]nucleus|10x|spatial transcript)",
        re.I,
    )
    lung = re.compile(r"(lung|nsclc|sclc|luad|lusc|pulmonary|bronch|pleural|alveolar)", re.I)
    rows = []
    for omix, prj, title, org, access, rel in row_re.findall(raw):
        title_s, org_s, access_s = _strip(title), _strip(org), access.strip()
        if not access_s.lower().startswith("open"):
            continue
        if "Homo sapiens" not in org_s:
            continue
        if not (sc.search(title_s) and lung.search(title_s)):
            continue
        rows.append(
            {
                "accession": omix,
                "bioproject": _strip(prj),
                "access": access_s,
                "organism": org_s,
                "release_date": rel.strip(),
                "title": title_s,
                "url": f"https://ngdc.cncb.ac.cn/omix/release/{omix}",
            }
        )
    return rows


def load_meta(path: str) -> list[dict]:
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    # first column name is often blank (cell id)
    id_key = [k for k in rows[0].keys() if k.strip() == ""][0]
    for r in rows:
        r["cell_id"] = r[id_key]
        r["nCount_RNA"] = float(r["nCount_RNA"])
    return rows


def stream_genes(matrix_path: str, wanted: set[str], cell_ids: list[str]) -> tuple[dict[str, np.ndarray], dict]:
    """One pass. Returns gene -> counts aligned to cell_ids, plus QC."""
    wanted_u = {g.upper(): g for g in wanted}
    found: dict[str, np.ndarray] = {}
    n_cells = len(cell_ids)
    qc = {"n_genes": 0, "n_header_cells": 0, "header_matches_meta": False, "colsum_cell0": None, "integerish": None}
    with open(matrix_path, "r", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split(",")
        # leading empty field, then cell ids
        header_cells = header[1:]
        qc["n_header_cells"] = len(header_cells)
        if header_cells == cell_ids:
            qc["header_matches_meta"] = True
            col_index = None  # already aligned
        else:
            pos = {c: i for i, c in enumerate(header_cells)}
            missing = [c for c in cell_ids if c not in pos]
            if missing:
                raise SystemExit(f"{len(missing)} metadata cells missing from matrix header, e.g. {missing[:3]}")
            col_index = np.array([pos[c] for c in cell_ids], dtype=np.int32)
        # column 0 of a mid-file gene, summed across ALL genes, to compare with nCount
        sum0 = 0.0
        integerish = True
        checked = 0
        for line in fh:
            if not line.strip():
                continue
            qc["n_genes"] += 1
            # gene name is the first field; counts follow. Gene symbols have no commas.
            comma = line.find(",")
            gene = line[:comma].strip().strip('"')
            key = wanted_u.get(gene.upper())
            vals = np.fromstring(line[comma + 1:], sep=",", dtype=np.float64)
            if vals.size != len(header_cells):
                # trailing newline already excluded by fromstring; size mismatch is real
                if vals.size == len(header_cells) + 1 and vals[-1] == 0 and line.rstrip("\n").endswith(","):
                    vals = vals[:-1]
                else:
                    raise SystemExit(f"{gene}: got {vals.size} values, header cells {len(header_cells)}")
            if col_index is not None:
                aligned = vals[col_index]
            else:
                aligned = vals
            if key is not None and key not in found:
                found[key] = aligned
            # QC on the first data column (header order), not a random cell
            sum0 += float(vals[0])
            if checked < 30:
                sample = vals[:20]
                if np.any(np.abs(sample - np.rint(sample)) > 1e-6):
                    integerish = False
                checked += 1
        qc["colsum_cell0"] = sum0
        qc["integerish"] = integerish
    return found, qc


def family_matrix(found: dict[str, np.ndarray], genes: list[str], logcpm: bool, ncount: np.ndarray) -> tuple[np.ndarray, list[str]]:
    present = []
    mats = []
    for g in genes:
        # prefer modern symbol; TMEM173 is an alias of STING1
        if g not in found:
            continue
        if g == "TMEM173" and "STING1" in found:
            continue
        present.append(g)
        mats.append(found[g])
    if not mats:
        return np.zeros((0, ncount.size)), []
    counts = np.vstack(mats)  # genes x cells
    if logcpm:
        scale = np.maximum(ncount, 1.0)
        expr = np.log1p(counts / scale * 1e4)
    else:
        expr = counts
    return expr.mean(axis=0), present


def summarize_patients(meta: list[dict], found: dict[str, np.ndarray], malignant_idx: np.ndarray) -> tuple[list[dict], list[dict], dict]:
    ncount = np.array([float(meta[i]["nCount_RNA"]) for i in range(len(meta))])
    cldn4 = found.get("CLDN4")
    if cldn4 is None:
        raise SystemExit("CLDN4 is not in the matrix")
    logcpm_cldn4 = np.log1p(cldn4 / np.maximum(ncount, 1.0) * 1e4)

    scores = {}
    present_map = {}
    for name, genes in FAMILIES.items():
        sc, present = family_matrix(found, genes, True, ncount)
        scores[name] = sc
        present_map[name] = present
    key_log = {}
    for g in KEY_GENES:
        if g in found:
            key_log[g] = np.log1p(found[g] / np.maximum(ncount, 1.0) * 1e4)

    ncount_arr = ncount
    nfeature = np.array([float(meta[i]["nFeature_RNA"]) for i in range(len(meta))])
    patients = sorted({meta[i]["patient"] for i in malignant_idx})
    per_patient = []
    for p in patients:
        idx = np.array([i for i in malignant_idx if meta[i]["patient"] == p], dtype=int)
        n = int(idx.size)
        c = cldn4[idx]
        n_pos = int(np.sum(c > 0))
        n_zero = n - n_pos
        rec = {
            "patient": p,
            "n_malignant": n,
            "n_CLDN4_pos": n_pos,
            "n_CLDN4_zero": n_zero,
            "pct_CLDN4_pos": n_pos / n if n else "",
            "median_CLDN4_count": float(np.median(c)) if n else "",
            "mean_CLDN4_logcpm": float(logcpm_cldn4[idx].mean()) if n else "",
            "spearman_included": "0",
            "poszero_included": "0",
            "reason": "",
        }
        if n < MIN_MALIGNANT:
            rec["reason"] = f"n_malignant<{MIN_MALIGNANT}"
            per_patient.append(rec)
            continue
        # Primary: continuous CLDN4 vs scores in every patient with enough malignant cells.
        rec["spearman_included"] = "1"
        rec["reason"] = f"spearman n>={MIN_MALIGNANT}"
        x = logcpm_cldn4[idx]
        for name, sc in scores.items():
            rho, pval = stats.spearmanr(x, sc[idx])
            rec[f"spearman_CLDN4_{name}"] = float(rho)
            rec[f"spearman_CLDN4_{name}_p"] = float(pval)
        rho_lib, p_lib = stats.spearmanr(x, np.log1p(ncount_arr[idx]))
        rec["spearman_CLDN4_log_nCount"] = float(rho_lib)
        rec["spearman_CLDN4_log_nCount_p"] = float(p_lib)
        rho_nf, p_nf = stats.spearmanr(x, nfeature[idx])
        rec["spearman_CLDN4_nFeature"] = float(rho_nf)
        rec["spearman_CLDN4_nFeature_p"] = float(p_nf)
        if n_pos < MIN_GROUP or n_zero < MIN_GROUP:
            rec["reason"] += f"; pos/zero group <{MIN_GROUP}"
            per_patient.append(rec)
            continue
        rec["poszero_included"] = "1"
        hi = idx[c > 0]
        lo = idx[c == 0]
        rec["n_high"] = int(hi.size)
        rec["n_low"] = int(lo.size)
        rec["median_nCount_high"] = float(np.median(ncount_arr[hi]))
        rec["median_nCount_low"] = float(np.median(ncount_arr[lo]))
        rec["median_nFeature_high"] = float(np.median(nfeature[hi]))
        rec["median_nFeature_low"] = float(np.median(nfeature[lo]))
        rec["mean_CLDN4_logcpm_high"] = float(logcpm_cldn4[hi].mean())
        rec["mean_CLDN4_logcpm_low"] = float(logcpm_cldn4[lo].mean())
        for name, sc in scores.items():
            rec[f"{name}_high"] = float(sc[hi].mean())
            rec[f"{name}_low"] = float(sc[lo].mean())
            rec[f"{name}_delta_high_minus_low"] = rec[f"{name}_high"] - rec[f"{name}_low"]
        for g, arr in key_log.items():
            rec[f"{g}_delta_high_minus_low"] = float(arr[hi].mean() - arr[lo].mean())
        per_patient.append(rec)

    spearman_set = [r for r in per_patient if r["spearman_included"] == "1"]
    included = [r for r in per_patient if r["poszero_included"] == "1"]
    summary = []
    for name in ["STING", "IFN", "APM", "TJ"]:
        vals = np.array([r[f"spearman_CLDN4_{name}"] for r in spearman_set], dtype=float)
        row = _delta_row(f"spearman_CLDN4_{name}", vals, spearman_set, "within_patient_spearman")
        row["n_genes_present"] = len(present_map[name])
        row["genes_present"] = ",".join(present_map[name])
        summary.append(row)
    for label in ("spearman_CLDN4_log_nCount", "spearman_CLDN4_nFeature"):
        vals = np.array([r[label] for r in spearman_set], dtype=float)
        summary.append(_delta_row(label, vals, spearman_set, "within_patient_spearman"))
    for name in ["STING", "IFN", "APM", "TJ"]:
        deltas = np.array([r[f"{name}_delta_high_minus_low"] for r in included], dtype=float)
        row = _delta_row(name, deltas, included, "pos_minus_zero_log1p_cpm")
        row["n_genes_present"] = len(present_map[name])
        row["genes_present"] = ",".join(present_map[name])
        summary.append(row)
    for g in KEY_GENES:
        if not included or f"{g}_delta_high_minus_low" not in included[0]:
            continue
        deltas = np.array([r[f"{g}_delta_high_minus_low"] for r in included], dtype=float)
        summary.append(_delta_row(g, deltas, included, "gene_pos_minus_zero"))

    # compartment detection (all cells, not the contrast)
    det = []
    for label, pred in [
        ("malignant_tumor", lambda r: r["cell_type"] == "tumor"),
        ("author_normal_epithelial", lambda r: r["cell_type"] == "normal"),
        ("Tcell", lambda r: r["cell_type"] == "Tcell"),
        ("myeloid", lambda r: r["cell_type"] == "Mye"),
        ("Bcell", lambda r: r["cell_type"] == "Bcell"),
    ]:
        idx = np.array([i for i, r in enumerate(meta) if pred(r)], dtype=int)
        if idx.size == 0 or "CLDN4" not in found:
            continue
        det.append(
            {
                "compartment": label,
                "n_cells": int(idx.size),
                "pct_CLDN4_pos": float(np.mean(found["CLDN4"][idx] > 0)),
                "mean_CLDN4_logcpm": float(logcpm_cldn4[idx].mean()),
                "mean_IFN": float(scores["IFN"][idx].mean()),
                "mean_APM": float(scores["APM"][idx].mean()),
                "mean_STING": float(scores["STING"][idx].mean()),
            }
        )
    extra = {"detection": det, "present": present_map, "n_patients_total": len(patients)}
    return per_patient, summary, extra


def _delta_row(name: str, deltas: np.ndarray, included: list[dict], kind: str) -> dict:
    n = int(deltas.size)
    row = {
        "feature": name,
        "kind": kind,
        "n_patients": n,
        "median_delta": "",
        "mean_delta": "",
        "n_negative": "",
        "n_positive": "",
        "n_zero": "",
        "wilcoxon_stat": "",
        "wilcoxon_p_two_sided": "",
    }
    if n == 0:
        return row
    row["median_delta"] = float(np.median(deltas))
    row["mean_delta"] = float(np.mean(deltas))
    row["n_negative"] = int(np.sum(deltas < 0))
    row["n_positive"] = int(np.sum(deltas > 0))
    row["n_zero"] = int(np.sum(deltas == 0))
    # signed-rank needs variation
    if n >= 5 and np.any(deltas != 0):
        stat, p = stats.wilcoxon(deltas, alternative="two-sided", zero_method="wilcox")
        row["wilcoxon_stat"] = float(stat)
        row["wilcoxon_p_two_sided"] = float(p)
    return row


def fmt(x, nd=4):
    if x == "" or x is None:
        return ""
    if isinstance(x, float):
        return f"{x:.6g}"
    return str(x)


def write_report(path: str, hra_rows: list[dict], omix_rows: list[dict], qc: dict, per_patient: list[dict], summary: list[dict], extra: dict, meta_n: int) -> None:
    open_rows = [r for r in hra_rows if r["access"] == "Open"]
    ctrl = [r for r in hra_rows if r["access"] == "Controlled"]
    spearman_pts = [r for r in per_patient if r.get("spearman_included") == "1"]
    included = [r for r in per_patient if r.get("poszero_included") == "1"]
    lines = []
    a = lines.append
    a("# Open HRA/CNGB Chinese human lung scRNA: malignant CLDN4 vs STING / IFN / APM")
    a("")
    a("Date of GSA-Human `finished.json`: pulled 2026-09-21 (7,312 records). OMIX release list pulled the same day (10,740 rows). Controlled HRA files were not requested.")
    a("")
    a("## Access")
    a("")
    a(f"Lung single-cell / single-nucleus / spatial **titles** in GSA-Human: **{len(hra_rows)}**. Open: **{len(open_rows)}**. Controlled: **{len(ctrl)}**.")
    a("")
    a("Title search misses a study whose title never says lung and never says single-cell. Every Open lung *disease* title in the same dump was read; the only Open single-cell/single-nucleus lung titles are the two below.")
    a("")
    a("### Open HRA (the whole open list)")
    a("")
    a("| Accession | What it is | Per-cell CLDN4 contrast |")
    a("|---|---|---|")
    a("| **HRA009335** (PRJCA031820) | Primary pulmonary lymphoepithelioma-like carcinoma, 4 snRNA runs (3 tumor, 1 adjacent). Sun Yat-sen. Open FASTQ only. | **Not scored.** PE150. In 4,000 read pairs from HRR2054880, R1 has a poly(T) stretch in the first 60 bp in 3,780/4,000 reads; R2 carries the Smart-seq TSO `AAGCAGTGGTATCAACGCAGAGTAC` in 215/4,000. Neither read is a 10x cell-barcode + UMI. One sample index (`TAGGACGT`) for the whole run. Processed matrix is not in the open deposit (paper: corresponding author on request). |")
    a("| **HRA011368** (PRJCA038849) | Single-nucleus RNA-seq of human **fetal** lung, twin-twin transfusion. Open FASTQ (HRR2360386 R1 ~10 GB, R2 ~25 GB). | Not a malignant-cell contrast. Not downloaded. |")
    a("")
    a("HRA009335 is the only Open human **lung tumor** single-nucleus accession. It cannot support a malignant CLDN4-high vs CLDN4-low score.")
    a("")
    a("### Controlled HRA lung sc/sn/spatial (not used)")
    a("")
    a(f"{len(ctrl)} titles. Raw files were not downloaded. Full table: `tables/hra_lung_scrna.tsv`. Examples that would have been the LUAD test if they were open: HRA001130 (lineage LUAD), HRA000154 / HRA000156 (subsolid / early LUAD), HRA002376 (preinvasive LUAD), HRA004391 (neoadjuvant immunotherapy LUAD), HRA012291 is not in the title-sc list under that accession; HRA008879 is the early-stage LUAD scRNA on PRJCA030984.")
    a("")
    a("### Open processed tables (OMIX), same Chinese archives")
    a("")
    a("These are Open OMIX releases, not Open HRA. The HRA raw run for the SCLC matrix is **Controlled**.")
    a("")
    a("| OMIX | BioProject | HRA of that BioProject | Used? |")
    a("|---|---|---|---|")
    a("| **OMIX002441** | PRJCA006026 | **HRA001232 Controlled** | **Yes.** 5,025 cells, author counts + cell info. SCLC primary tumors, matched adjacent, one relapse. |")
    a("| OMIX011746 | PRJCA045688 | no HRA row in finished.json | No. 1.67 GB Seurat RDS of pleural effusion and blood, TB vs LUAD **immunity**. Not a malignant-cell matrix in the file title. |")
    a("| OMIX003147 | PRJCA015245 | not a tumor HRA | No. Embryonic lung. |")
    a("| OMIX008248 | PRJCA033575 | COVID white lung | No. Neutrophil / COVID, not tumor CLDN4. |")
    a("| OMIX004145 | PRJCA017221 | tuberculous pleural effusion | No. Pleural immune comparison, not lung-tumor malignant cells. |")
    a("| OMIX007208 | PRJCA029532 | pediatric Mycoplasma BAL | No. Bronchoalveolar immune atlas, not tumor. |")
    a("")
    a("### CNGB / CNSA")
    a("")
    a("CNP0005129 (Tianjin Chest Hospital; LUAD scRNA + spatial, STAS; 14 samples, 535 GB) is listed on CNGBdb as **Apply for data**, with metadata download only. That is controlled. It was not downloaded.")
    a("")
    a("CNGBdb project search `lung single-cell` (235 hits) and `lung adenocarcinoma single cell` (25 hits) returned NCBI BioProject mirrors, not a public CNP tumor matrix. No other CNP lung-tumor scRNA page was confirmed Open in this pass. Absence from that search is not proof that no public CNP exists under a title that omits those words.")
    a("")
    a("## Analysis that was possible: OMIX002441 SCLC malignant cells")
    a("")
    a("Source: `OMIX002441-01.csv` (gene × cell counts) and `OMIX002441-02.csv` (cell info). Paper: Signal Transduction and Targeted Therapy, DOI 10.1038/s41392-022-01150-4 (Fuchou Tang / Peking University).")
    a("")
    a("This is **small-cell lung cancer**, not LUAD. It is not added to the concordant-4 LUAD T/NK correlation. A result here does not rewrite that LUAD result.")
    a("")
    a("**Malignant** = author `cell_type == tumor` (same cells have `NT == tumor`). n = 2,104 / 5,025. Immune and stromal cells are labeled `NT == normal` in this table even when the sample code is a tumor region, so they are not used as the low-CLDN4 malignant group.")
    a("")
    a(f"Matrix QC: {qc.get('n_genes')} genes, header cells {qc.get('n_header_cells')}, header matches metadata: {qc.get('header_matches_meta')}. First 30 genes look like integers: {qc.get('integerish')}. Sum of all genes in matrix column 0 = {fmt(qc.get('colsum_cell0'))}; metadata nCount_RNA for that cell = {fmt(qc.get('ncount_cell0'))}.")
    a("")
    a("Score = **mean of log1p(count / nCount_RNA × 10,000)** over genes present in the matrix. Library size is the author `nCount_RNA`, checked against the column sum above.")
    a("")
    a("**Primary contrast:** within each patient, Spearman correlation of malignant-cell CLDN4 log1p CPM with the score. Patients with ≥ 40 malignant cells are included. This keeps patients in whom almost every malignant cell is CLDN4-positive.")
    a("")
    a("**Secondary contrast:** CLDN4 count > 0 vs = 0, only when both groups have ≥ 15 cells. Quartiles are not used. In this matrix most malignant cells are CLDN4-positive, so the zero class is a small tail and the patients with the highest CLDN4 (too few zeros) drop out of the secondary contrast.")
    a("")
    a("A positive Spearman means higher CLDN4 goes with a higher score inside that patient's malignant cells. A positive pos-minus-zero delta means the CLDN4-positive malignant cells score higher than the CLDN4-zero malignant cells from the same patient. p is a two-sided Wilcoxon signed-rank test across patients (SciPy `zero_method='wilcox'`). It is descriptive. For 5 non-zero pairs that all share a sign, the smallest two-sided p this test can return is 0.0625. That floor is not a significance claim.")
    a("")
    a(f"Patients with any malignant cells: {extra['n_patients_total']}. Spearman set (n≥40): **{len(spearman_pts)}** ({', '.join(r['patient'] for r in spearman_pts)}). Pos-vs-zero set: **{len(included)}** ({', '.join(r['patient'] for r in included)}).")
    a("")
    a("### CLDN4 detection")
    a("")
    a("| compartment | n cells | % CLDN4 > 0 | mean log1p CPM |")
    a("|---|---:|---:|---:|")
    for d in extra["detection"]:
        a(f"| {d['compartment']} | {d['n_cells']} | {d['pct_CLDN4_pos']:.3f} | {d['mean_CLDN4_logcpm']:.3f} |")
    a("")
    a("### Patient-level result")
    a("")
    a("Delta = CLDN4-positive minus CLDN4-zero. For Spearman rows, the number is the within-patient Spearman of CLDN4 log1p CPM vs the score (not a high-minus-low delta).")
    a("")
    a("| feature | kind | n patients | median | mean | n down / up | Wilcoxon p |")
    a("|---|---|---:|---:|---:|---|---:|")
    show_genes = {"CLDN4", "STING1", "TBK1", "STAT1", "ISG15", "MX1", "HLA-A", "B2M", "TAP1", "EPCAM", "PTPRC", "CD3D", "CD68"}
    for r in summary:
        if r["kind"] == "gene_pos_minus_zero" and r["feature"] not in show_genes:
            continue
        nd = r["n_negative"]
        nu = r["n_positive"]
        a(
            f"| {r['feature']} | {r['kind']} | {r['n_patients']} | {fmt(r['median_delta'])} | {fmt(r['mean_delta'])} | {nd}/{nu} | {fmt(r['wilcoxon_p_two_sided'])} |"
        )
    a("")
    a("Genes present in each score are in `tables/sclc_family_summary.tsv`.")
    a("")
    a("### How to read it")
    a("")
    a("Author-malignant SCLC cells are mostly CLDN4-positive (about 90% of 2,104 cells; mean log1p CPM about 1.25). Adjacent epithelial cells labeled `normal` are lower (about 68%). Immune compartments are lower still (T cells about 25%). CLDN4 is measured in this matrix. It is not a LUAD cohort, and these patients are not added to concordant-4.")
    a("")
    a("Primary result, within-patient Spearman, 9 patients with ≥ 40 malignant cells (P2, P3, P4, P5, P7, P10, P11, P12, P13):")
    a("")
    a("- **STING** median ρ = −0.0004 (5/9 negative, 4/9 positive), Wilcoxon p = 0.57. Flat.")
    a("- **IFN** median ρ = +0.071 (7/9 positive), p = 0.13.")
    a("- **APM** median ρ = +0.097 (6/9 positive), p = 0.13.")
    a("- **TJ** (CLDN4 held out) median ρ = +0.30 (9/9 positive), p = 0.0039.")
    a("- CLDN4 vs log nCount median ρ = +0.053 (8/9 positive), p = 0.020. CLDN4 vs nFeature median ρ = +0.11 (8/9 positive), p = 0.012. The APM ρ is the same size as the nFeature ρ.")
    a("")
    a("The two patients with the highest CLDN4 detection (P2, 98.5% positive, median count 48; P7, 95% positive) have negative IFN Spearman (−0.15 and −0.01). The positive median is not those high-CLDN4 tumors.")
    a("")
    a("Secondary pos-vs-zero contrast is 5 patients (P4, P10, P11, P12, P13). CLDN4-zero malignant cells are the low-UMI tail: median nCount in the zero class is about half the positive class (P10 11,386 vs 32,438; P4 13,400 vs 23,516; P11 23,016 vs 42,187). IFN and APM means are slightly higher in the positive class (median deltas +0.027 and +0.066; 5/5 up; p = 0.0625, which is the floor of this test at n = 5). That contrast is library-size confounded. STING score delta median is +0.015 (4/5 up, p = 0.19). STING1, TBK1, and IRF3 do not move as one induced program (TBK1 median delta −0.045, 3/5 down).")
    a("")
    a("Immune-marker check on the same 5 patients, pos minus zero: PTPRC median delta −0.008 (3/5 down), CD3D +0.003, CD68 −0.018 (4/5 down). The CLDN4-positive class is not an immune doublet.")
    a("")
    a("TJ is a control family with **CLDN4 held out**. It is not the claim.")
    a("")
    a("## What was not done")
    a("")
    a("- No DAC application, no controlled FASTQ, no Cell Ranger on HRA009335 (no cell barcode in the reads).")
    a("- OMIX011746 RDS was not loaded.")
    a("- CNP0005129 was not downloaded.")
    a("- Concordant-4 GEO cohorts were not re-run and were not mixed with this SCLC table.")
    a("")
    a("## Files")
    a("")
    a("| file | |")
    a("|---|---|")
    a("| `tables/hra_lung_scrna.tsv` | every lung sc/sn/spatial HRA title, Open or Controlled |")
    a("| `tables/omix_open_human_lung_scrna.tsv` | Open human OMIX lung single-cell titles |")
    a("| `tables/sclc_patient_deltas.tsv` | one row per SCLC patient |")
    a("| `tables/sclc_family_summary.tsv` | paired deltas and Spearman |")
    a("| `tables/sclc_detection.tsv` | CLDN4 detection by author compartment |")
    a("| `tables/qc.json` | matrix QC |")
    a("")
    a(f"Cells in the metadata table: {meta_n}.")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--finished-json", default="/tmp/hra/finished.json")
    ap.add_argument("--omix-html", default="/tmp/hra/omix.html")
    ap.add_argument("--meta", default="/tmp/hra/omix/SCLC_cell_info.csv")
    ap.add_argument("--matrix", default="/tmp/hra/omix/SCLC_cell_count_matrix.csv")
    args = ap.parse_args()
    os.makedirs(os.path.join(OUT, "tables"), exist_ok=True)

    hra_rows = build_hra_catalog(args.finished_json)
    write_tsv(
        os.path.join(OUT, "tables", "hra_lung_scrna.tsv"),
        hra_rows,
        ["accession", "bioproject", "access", "is_controlled", "title", "organization", "release_time_string", "study_type_id", "url"],
    )
    omix_rows = build_omix_catalog(args.omix_html)
    write_tsv(
        os.path.join(OUT, "tables", "omix_open_human_lung_scrna.tsv"),
        omix_rows,
        ["accession", "bioproject", "access", "organism", "release_date", "title", "url"],
    )

    meta = load_meta(args.meta)
    cell_ids = [r["cell_id"] for r in meta]
    wanted = set(STING + IFN + APM + TJ + KEY_GENES + ["CLDN4", "TMEM173", "MB21D1"])
    found, qc = stream_genes(args.matrix, wanted, cell_ids)
    # cGAS is CGAS in current symbols; older matrices use MB21D1.
    if "CGAS" not in found and "MB21D1" in found:
        found["CGAS"] = found["MB21D1"]
    if "STING1" not in found and "TMEM173" in found:
        found["STING1"] = found["TMEM173"]
    by_id = {r["cell_id"]: r for r in meta}
    # column-0 sum is compared to the metadata cell that occupies matrix column 0
    # stream_genes records that sum; the id is the first metadata cell only when
    # the header is in metadata order. Re-read the first header field.
    with open(args.matrix, "r", encoding="utf-8", errors="replace") as fh:
        header0 = fh.readline().rstrip("\n").split(",")[1]
    qc["matrix_cell0"] = header0
    qc["ncount_cell0"] = by_id[header0]["nCount_RNA"]
    qc["genes_found"] = sorted(found)
    malignant_idx = np.array([i for i, r in enumerate(meta) if r["cell_type"] == "tumor"], dtype=int)
    per_patient, summary, extra = summarize_patients(meta, found, malignant_idx)

    # stable field order
    pat_fields = []
    for r in per_patient:
        for k in r:
            if k not in pat_fields:
                pat_fields.append(k)
    write_tsv(os.path.join(OUT, "tables", "sclc_patient_deltas.tsv"), per_patient, pat_fields)
    sum_fields = ["feature", "kind", "n_patients", "n_genes_present", "genes_present", "median_delta", "mean_delta", "n_negative", "n_positive", "n_zero", "wilcoxon_stat", "wilcoxon_p_two_sided"]
    write_tsv(os.path.join(OUT, "tables", "sclc_family_summary.tsv"), summary, sum_fields)
    write_tsv(
        os.path.join(OUT, "tables", "sclc_detection.tsv"),
        extra["detection"],
        ["compartment", "n_cells", "pct_CLDN4_pos", "mean_CLDN4_logcpm", "mean_IFN", "mean_APM", "mean_STING"],
    )
    with open(os.path.join(OUT, "tables", "qc.json"), "w", encoding="utf-8") as fh:
        json.dump({k: v for k, v in qc.items() if k != "genes_found"} | {"genes_found": qc["genes_found"]}, fh, indent=2)
        fh.write("\n")

    write_report(
        os.path.join(OUT, "REPORT.md"),
        hra_rows,
        omix_rows,
        qc,
        per_patient,
        summary,
        extra,
        len(meta),
    )
    print("Wrote", OUT)
    print("Open HRA", sum(1 for r in hra_rows if r["access"] == "Open"), "of", len(hra_rows))
    print("QC", {k: qc[k] for k in ("n_genes", "n_header_cells", "header_matches_meta", "integerish", "colsum_cell0", "ncount_cell0")})
    for r in summary:
        if r["kind"] == "within_patient_spearman":
            print(r["feature"], "n", r["n_patients"], "median", r["median_delta"], "neg", r["n_negative"], "pos", r["n_positive"], "p", r["wilcoxon_p_two_sided"])


if __name__ == "__main__":
    main()
