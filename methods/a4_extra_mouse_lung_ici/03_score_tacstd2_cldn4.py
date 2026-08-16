#!/usr/bin/env python3
"""Score Tacstd2 and Cldn4 in extra (non-TISMO) public mouse lung ICI / KL RNA."""
from __future__ import annotations

import csv
import gzip
import json
import math
import re
from pathlib import Path

import numpy as np
from scipy import stats

DATA = Path("notes/a4_extra_mouse_lung_ici/data")
OUT = Path("results/a4_extra_mouse_lung_ici")
OUT.mkdir(parents=True, exist_ok=True)

# ENSMUSG00000051397 = Tacstd2 (confirmed in GSE244452 / GSE338923 tables).
# ENSMUSG00000030798 is Cd37 — do not use it.
TACSTD2 = {
    "tacstd2",
    "trop2",
    "ensmusg00000051397",
    "ensmusg00000051397.1",
    "ensmusg00000051397.2",
    "ensmusg00000051397.3",
    "ensmusg00000051397.4",
    "ensmusg00000051397.5",
}
CLDN4 = {
    "cldn4",
    "ensmusg00000047501",
    "ensmusg00000047501.1",
    "ensmusg00000047501.2",
    "ensmusg00000047501.3",
}
# immune controls
CD8A = {"cd8a", "ensmusg00000053044"}
CD274 = {"cd274", "pdl1", "pd-l1", "ensmusg00000016496"}
STK11 = {"stk11", "lkb1", "ensmusg00000003032"}


def gzopen(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", errors="replace")
    return open(path, "rt", errors="replace")


def parse_series_matrix(path: Path) -> dict:
    meta = {"title": "", "summary": "", "design": "", "samples": []}
    titles, accs, chars = [], [], []
    if not path.exists() or path.stat().st_size < 200:
        return meta
    with gzopen(path) as fh:
        for line in fh:
            if line.startswith("!Series_title"):
                meta["title"] = line.split("\t", 1)[-1].strip().strip('"')
            elif line.startswith("!Series_summary"):
                meta["summary"] += " " + line.split("\t", 1)[-1].strip().strip('"')
            elif line.startswith("!Series_overall_design"):
                meta["design"] += " " + line.split("\t", 1)[-1].strip().strip('"')
            elif line.startswith("!Sample_title"):
                titles = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_geo_accession"):
                accs = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1") or line.startswith("!Sample_source_name_ch1"):
                chars.append([x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]])
            elif line.startswith("!series_matrix_table_begin"):
                break
    n = max(len(titles), len(accs), max((len(c) for c in chars), default=0))
    samples = []
    for i in range(n):
        s = {
            "title": titles[i] if i < len(titles) else "",
            "gsm": accs[i] if i < len(accs) else "",
            "chars": [c[i] for c in chars if i < len(c)],
        }
        samples.append(s)
    meta["samples"] = samples
    return meta


def norm_key(x: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(x).lower())


def gene_match(token: str, targets: set[str]) -> bool:
    t = str(token).strip().strip('"').lower()
    if t in targets:
        return True
    t0 = t.split(".")[0]
    if t0 in targets:
        return True
    # symbol in "id|symbol" or "ensembl_symbol"
    for part in re.split(r"[|;,\s/]+", t):
        if part in targets or part.split(".")[0] in targets:
            return True
    return False


def find_gene_row(header, rows, targets: set[str], id_cols: list[int] | None = None):
    """Return first matching numeric row as list[float|nan] aligned to sample columns."""
    if id_cols is None:
        id_cols = [0]
        if header:
            for i, h in enumerate(header[:8]):
                hl = str(h).lower()
                if any(k in hl for k in ("gene", "symbol", "id", "name")):
                    id_cols.append(i)
            id_cols = sorted(set(id_cols))
    for row in rows:
        hit = False
        for i in id_cols:
            if i < len(row) and gene_match(row[i], targets):
                hit = True
                break
        if not hit:
            continue
        vals = []
        start = max(id_cols) + 1
        for x in row[start:]:
            try:
                vals.append(float(str(x).replace(",", "")))
            except Exception:
                vals.append(np.nan)
        return row, vals
    return None, None


def read_delim(path: Path, sep=None, max_rows=None):
    with gzopen(path) as fh:
        first = fh.readline()
        if sep is None:
            if "\t" in first and first.count("\t") >= first.count(","):
                sep = "\t"
            elif ";" in first and first.count(";") > first.count(","):
                sep = ";"
            else:
                sep = ","
        header = [x.strip().strip('"') for x in first.rstrip("\n").split(sep)]
        rows = []
        for i, line in enumerate(fh):
            rows.append([x.strip().strip('"') for x in line.rstrip("\n").split(sep)])
            if max_rows and i >= max_rows:
                break
    return header, rows, sep


def welch_mwu(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    out = {
        "n_a": int(a.size),
        "n_b": int(b.size),
        "mean_a": float(np.mean(a)) if a.size else None,
        "mean_b": float(np.mean(b)) if b.size else None,
        "median_a": float(np.median(a)) if a.size else None,
        "median_b": float(np.median(b)) if b.size else None,
        "log2fc_mean": None,
        "welch_p": None,
        "mwu_p": None,
        "direction": None,
    }
    if a.size == 0 or b.size == 0:
        return out
    # log2FC on raw-ish values: if already log-like (range small / negatives), use mean diff
    if np.nanmin(np.concatenate([a, b])) < 0 or (np.nanmax(np.concatenate([a, b])) < 30 and np.nanmin(np.concatenate([a, b])) >= 0 and np.nanmedian(np.concatenate([a, b])) < 15):
        # treat as log-space if negatives present; else still compute log2((b+eps)/(a+eps))
        if np.nanmin(np.concatenate([a, b])) < 0:
            out["log2fc_mean"] = float(np.mean(b) - np.mean(a))
            out["space"] = "already_log_mean_diff"
        else:
            out["log2fc_mean"] = float(np.log2((np.mean(b) + 1e-6) / (np.mean(a) + 1e-6)))
            out["space"] = "log2_mean_plus_eps"
    else:
        out["log2fc_mean"] = float(np.log2((np.mean(b) + 1.0) / (np.mean(a) + 1.0)))
        out["space"] = "log2_mean_plus1"
    if a.size >= 2 and b.size >= 2:
        try:
            out["welch_p"] = float(stats.ttest_ind(a, b, equal_var=False, nan_policy="omit").pvalue)
        except Exception:
            out["welch_p"] = None
        try:
            out["mwu_p"] = float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)
        except Exception:
            out["mwu_p"] = None
    elif a.size == 1 or b.size == 1:
        out["welch_p"] = None
        out["mwu_p"] = None
        out["note"] = "n=1 arm; p not computed"
    if out["log2fc_mean"] is not None:
        if out["log2fc_mean"] > 0:
            out["direction"] = "up_in_b"
        elif out["log2fc_mean"] < 0:
            out["direction"] = "down_in_b"
        else:
            out["direction"] = "tie"
    return out


def bh_fdr(pvals):
    p = np.array([np.nan if v is None else v for v in pvals], dtype=float)
    out = [None] * len(p)
    mask = np.isfinite(p)
    if not mask.any():
        return out
    idx = np.where(mask)[0]
    pv = p[idx]
    order = np.argsort(pv)
    ranked = pv[order]
    m = len(ranked)
    q = ranked * m / (np.arange(1, m + 1))
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    for i, qi in zip(idx[order], q):
        out[i] = float(qi)
    return out


def extract_from_table(path: Path, sample_names: list[str] | None = None, gene_col_hints=None):
    header, rows, sep = read_delim(path)
    # drop empty trailing
    while header and header[-1] == "":
        header = header[:-1]
    genes = {}
    for name, targets in (("Tacstd2", TACSTD2), ("Cldn4", CLDN4), ("Cd8a", CD8A), ("Cd274", CD274), ("Stk11", STK11)):
        raw, vals = find_gene_row(header, rows, targets)
        if vals is None:
            genes[name] = None
        else:
            genes[name] = {"raw_id": raw[0] if raw else name, "values": vals, "n_cols": len(vals)}
    return {"header": header, "genes": genes, "n_rows": len(rows), "sep": sep}


def map_samples_to_values(header, values, sample_keys):
    """Map sample identifiers to values using fuzzy header match."""
    hnorm = [norm_key(h) for h in header]
    # values align to header[offset:]
    offset = len(header) - len(values)
    if offset < 0:
        offset = 0
    out = {}
    for key in sample_keys:
        nk = norm_key(key)
        best = None
        for i, h in enumerate(hnorm):
            if not h:
                continue
            if nk == h or nk in h or h in nk:
                vi = i - offset
                if 0 <= vi < len(values):
                    best = values[vi]
                    break
        out[key] = best
    return out


def xlsx_extract(path: Path):
    import openpyxl

    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    sheet = wb[wb.sheetnames[0]]
    rows = []
    for i, row in enumerate(sheet.iter_rows(values_only=True)):
        rows.append(["" if c is None else str(c) for c in row])
        if i > 200000:
            break
    wb.close()
    header = rows[0]
    body = rows[1:]
    genes = {}
    for name, targets in (("Tacstd2", TACSTD2), ("Cldn4", CLDN4), ("Cd8a", CD8A), ("Cd274", CD274)):
        raw, vals = find_gene_row(header, body, targets)
        genes[name] = None if vals is None else {"raw_id": raw[0] if raw else name, "values": vals}
    return {"header": header, "genes": genes, "n_rows": len(body), "sheet": sheet}


def values_by_titles(header, values, titles, pred):
    """Return values for samples whose title matches pred(title)."""
    # assume values correspond to last N header cols
    offset = len(header) - len(values)
    picked = []
    used = []
    for i, t in enumerate(titles):
        if not pred(t):
            continue
        # try title tokens against header
        nk = norm_key(t)
        found = None
        for j, h in enumerate(header):
            hn = norm_key(h)
            if not hn:
                continue
            if hn in nk or nk in hn or hn[:8] and hn[:8] in nk:
                vi = j - offset
                if 0 <= vi < len(values):
                    found = values[vi]
                    used.append((t, header[j]))
                    break
        if found is not None:
            picked.append(found)
    return picked, used


def main():
    catalog = []
    contrasts = []
    per_sample = []

    def add_contrast(acc, model, design, gene, group_a, group_b, a_vals, b_vals, note="", icb=False, response=False, genotype=False):
        st = welch_mwu(a_vals, b_vals)
        contrasts.append(
            {
                "accession": acc,
                "model": model,
                "design": design,
                "gene": gene,
                "group_a": group_a,
                "group_b": group_b,
                "icb_contrast": icb,
                "response_contrast": response,
                "genotype_contrast": genotype,
                "note": note,
                **st,
            }
        )

    def dump_samples(acc, model, gene, pairs):
        for sample, group, val in pairs:
            per_sample.append(
                {
                    "accession": acc,
                    "model": model,
                    "gene": gene,
                    "sample": sample,
                    "group": group,
                    "value": None if val is None or (isinstance(val, float) and not math.isfinite(val)) else float(val),
                }
            )

    # ---------- GSE114601 KP GEMM ICB ----------
    p = DATA / "GSE114601_counts.normalized.csv.gz"
    h, rows, _ = read_delim(p)
    meta = parse_series_matrix(DATA / "GSE114601_series_matrix.txt.gz")
    # header sample ids s1795...
    title_by_sid = {}
    treat_by_sid = {}
    for s in meta["samples"]:
        sid = s["title"].split()[0]
        title_by_sid[sid] = s["title"]
        treat = next((c.split(":", 1)[-1].strip() for c in s["chars"] if c.lower().startswith("treatment")), "")
        treat_by_sid[sid] = treat
    genes = {}
    for name, targets in (("Tacstd2", TACSTD2), ("Cldn4", CLDN4), ("Cd8a", CD8A), ("Cd274", CD274)):
        raw, vals = find_gene_row(h, rows, targets)
        genes[name] = dict(zip(h[1:], vals)) if vals is not None else None
    catalog.append({"accession": "GSE114601", "model": "KP GEMM nodules", "genes_found": {k: v is not None for k, v in genes.items()}})
    for gene, d in genes.items():
        if not d:
            continue
        veh = [d[s] for s, t in treat_by_sid.items() if t.lower() == "vehicle" and s in d]
        pd1 = [d[s] for s, t in treat_by_sid.items() if t.lower() == "anti-pd1" and s in d]
        jq1 = [d[s] for s, t in treat_by_sid.items() if t.lower() == "jq1" and s in d]
        combo = [d[s] for s, t in treat_by_sid.items() if "anti-pd1-jq1" in t.lower() and s in d]
        add_contrast("GSE114601", "KP GEMM lung nodules", "anti-PD1 vs vehicle", gene, "vehicle", "anti-PD1", veh, pd1, "n=2 vs 2; orthotopic GEMM", icb=True)
        add_contrast("GSE114601", "KP GEMM lung nodules", "anti-PD1+JQ1 vs vehicle", gene, "vehicle", "anti-PD1+JQ1", veh, combo, "BET+ICB combo", icb=True)
        for s, t in treat_by_sid.items():
            if s in d:
                dump_samples("GSE114601", "KP GEMM", gene, [(s, t, d[s])])

    # ---------- GSE157880 HKP1 RT ± PD-1 ----------
    p = DATA / "GSE157880_Bulk048.txt.gz"
    h, rows, _ = read_delim(p)
    meta = parse_series_matrix(DATA / "GSE157880_series_matrix.txt.gz")
    # columns after annotation
    sample_cols = h[8:]
    title_map = {}
    for s in meta["samples"]:
        # "0-1_S13 (IgG_0Gy_1)"
        left = s["title"].split()[0]
        title_map[left] = s["title"]
    genes = {}
    for name, targets in (("Tacstd2", TACSTD2), ("Cldn4", CLDN4), ("Cd8a", CD8A), ("Cd274", CD274)):
        raw, vals = find_gene_row(h, rows, targets, id_cols=[4, 6])
        if vals is None:
            genes[name] = None
        else:
            genes[name] = dict(zip(sample_cols, vals[: len(sample_cols)]))
    catalog.append({"accession": "GSE157880", "model": "HKP1 orthotopic lung", "genes_found": {k: v is not None for k, v in genes.items()}})

    def arm_from_title(t):
        m = re.search(r"\(([^)]+)\)", t)
        return m.group(1) if m else t

    for gene, d in genes.items():
        if not d:
            continue
        groups = {}
        for sid, title in title_map.items():
            arm = arm_from_title(title)
            if sid in d:
                groups.setdefault(arm.split("_")[0] + "_" + re.search(r"(\d+Gy)", arm).group(1) if re.search(r"(\d+Gy)", arm) else arm, [])
                key = re.sub(r"_\d+$", "", arm.replace(" ", "_"))
                groups.setdefault(key, []).append(d[sid])
                dump_samples("GSE157880", "HKP1 lung", gene, [(sid, key, d[sid])])
        igg0 = groups.get("IgG_0Gy", [])
        pd0 = groups.get("PD-1_0Gy", [])
        igg4 = groups.get("IgG_4Gy", [])
        pd4 = groups.get("PD-1_4Gy", [])
        add_contrast("GSE157880", "HKP1 orthotopic lung", "anti-PD-1 0 Gy vs IgG 0 Gy", gene, "IgG_0Gy", "PD-1_0Gy", igg0, pd0, "primary ICB without radiation", icb=True)
        add_contrast("GSE157880", "HKP1 orthotopic lung", "anti-PD-1 4 Gy vs IgG 0 Gy", gene, "IgG_0Gy", "PD-1_4Gy", igg0, pd4, "radiation confounder", icb=True)
        add_contrast("GSE157880", "HKP1 orthotopic lung", "IgG 4 Gy vs IgG 0 Gy", gene, "IgG_0Gy", "IgG_4Gy", igg0, igg4, "radiation-only control")

    # ---------- GSE239485 LLC polyI:C + aPD-1 ----------
    import openpyxl

    p = DATA / "GSE239485_Processed_data.xlsx"
    wb = openpyxl.load_workbook(p, read_only=True, data_only=True)
    sh = wb[wb.sheetnames[0]]
    rows = [["" if c is None else str(c) for c in r] for r in sh.iter_rows(values_only=True)]
    wb.close()
    h = rows[0]
    body = rows[1:]
    meta = parse_series_matrix(DATA / "GSE239485_series_matrix.txt.gz")
    treat = {}
    for s in meta["samples"]:
        tr = next((c.split(":", 1)[-1].strip() for c in s["chars"] if c.lower().startswith("treatment")), "")
        treat[s["title"]] = tr
        treat[s["title"].split(",")[0].strip()] = tr
    genes = {}
    for name, targets in (("Tacstd2", TACSTD2), ("Cldn4", CLDN4), ("Cd8a", CD8A), ("Cd274", CD274)):
        raw, vals = find_gene_row(h, body, targets)
        genes[name] = (h, vals, raw[0] if raw else None)
    catalog.append({"accession": "GSE239485", "model": "LLC tumors", "genes_found": {k: genes[k][1] is not None for k in genes}})
    for gene, pack in genes.items():
        header, vals, _ = pack
        if vals is None:
            continue
        offset = len(header) - len(vals)
        veh, combo, triple = [], [], []
        for j, col in enumerate(header[offset:]):
            if str(col).startswith("C_"):
                tr = "Vehicle"
                veh.append(vals[j])
            elif str(col).startswith("D_"):
                tr = "Poly I:C + anti-PD-1"
                combo.append(vals[j])
            elif str(col).startswith("T_"):
                tr = "Poly I:C + anti-PD-1 + C5aR1i"
                triple.append(vals[j])
            else:
                tr = "other"
            dump_samples("GSE239485", "LLC", gene, [(col, tr, vals[j])])
        add_contrast(
            "GSE239485",
            "LLC subcutaneous",
            "Poly I:C + anti-PD-1 vs vehicle",
            gene,
            "vehicle",
            "polyI:C+aPD-1",
            veh,
            combo,
            "NO monotherapy aPD-1 arm; Poly I:C confounder. C_=vehicle, D_=doublet, T_=triplet",
            icb=True,
        )
        add_contrast(
            "GSE239485",
            "LLC subcutaneous",
            "Poly I:C + anti-PD-1 + C5aR1i vs vehicle",
            gene,
            "vehicle",
            "triplet",
            veh,
            triple,
            "triple combo; still no aPD-1 monotherapy",
            icb=True,
        )

    # ---------- GSE262305 LLC aPD-L1 ----------
    p = DATA / "GSE262305_gene_expression.txt.gz"
    h, rows, _ = read_delim(p)
    meta = parse_series_matrix(DATA / "GSE262305_series_matrix.txt.gz")
    treat_by_title = {s["title"]: next((c.split(":", 1)[-1].strip() for c in s["chars"] if c.lower().startswith("treatment")), s["title"]) for s in meta["samples"]}
    genes = {}
    for name, targets in (("Tacstd2", TACSTD2), ("Cldn4", CLDN4), ("Cd8a", CD8A), ("Cd274", CD274)):
        raw, vals = find_gene_row(h, rows, targets)
        genes[name] = (h, vals)
    catalog.append({"accession": "GSE262305", "model": "LLC", "genes_found": {k: genes[k][1] is not None for k in genes}})
    for gene, (header, vals) in genes.items():
        if vals is None:
            continue
        offset = len(header) - len(vals)
        groups = {"isotype": [], "aPDL1": [], "BTZ+aPDL1": []}
        for j, col in enumerate(header[offset:]):
            if not str(col).endswith("_FPKM"):
                continue
            val = vals[j]
            cl = col.lower()
            if "blank" in cl:
                g = "isotype"
            elif "btz" in cl:
                g = "BTZ+aPDL1"
            elif "pdl1" in cl:
                g = "aPDL1"
            else:
                continue
            groups[g].append(val)
            dump_samples("GSE262305", "LLC", gene, [(col, g, val)])
        add_contrast("GSE262305", "LLC subcutaneous", "anti-PD-L1 vs isotype", gene, "isotype", "aPDL1", groups["isotype"], groups["aPDL1"], "monotherapy ICB; FPKM", icb=True)
        add_contrast("GSE262305", "LLC subcutaneous", "BTZ+anti-PD-L1 vs isotype", gene, "isotype", "BTZ+aPDL1", groups["isotype"], groups["BTZ+aPDL1"], "combo; FPKM", icb=True)

    # ---------- GSE274960 LL/2 aPD-1 ----------
    p = DATA / "GSE274960_gene_count_matrix.csv.gz"
    h, rows, _ = read_delim(p)
    meta = parse_series_matrix(DATA / "GSE274960_series_matrix.txt.gz")
    treat = {s["title"]: next((c.split(":", 1)[-1].strip() for c in s["chars"] if c.lower().startswith("treatment")), "") for s in meta["samples"]}
    genes = {}
    for name, targets in (("Tacstd2", TACSTD2), ("Cldn4", CLDN4), ("Cd8a", CD8A), ("Cd274", CD274)):
        raw, vals = find_gene_row(h, rows, targets)
        genes[name] = (h, vals)
    catalog.append({"accession": "GSE274960", "model": "LL/2 (LLC)", "genes_found": {k: genes[k][1] is not None for k in genes}})
    for gene, (header, vals) in genes.items():
        if vals is None:
            continue
        offset = len(header) - len(vals)
        groups = {}
        for j, col in enumerate(header[offset:]):
            # LL2_IgG_1 etc
            if "IgG" in col:
                g = "IgG"
            elif "E_PD" in col or "E+PD" in col:
                g = "Entrectinib+aPD1"
            elif re.search(r"PD_", col) or col.endswith("_PD") or "_PD_" in col:
                g = "aPD1"
            elif re.search(r"_E_\d", col):
                g = "Entrectinib"
            elif "shNTRK" in col:
                g = "shNTRK1"
            elif "CNTRL" in col or "CTRL" in col:
                g = "shCTRL"
            else:
                g = col
            groups.setdefault(g, []).append(vals[j])
            dump_samples("GSE274960", "LL/2", gene, [(col, g, vals[j])])
        add_contrast("GSE274960", "LL/2 LLC", "anti-PD-1 vs IgG", gene, "IgG", "aPD1", groups.get("IgG", []), groups.get("aPD1", []), "monotherapy ICB", icb=True)
        add_contrast("GSE274960", "LL/2 LLC", "Entrectinib+anti-PD-1 vs IgG", gene, "IgG", "Entrectinib+aPD1", groups.get("IgG", []), groups.get("Entrectinib+aPD1", []), "combo", icb=True)

    # ---------- GSE297630 LLC aPD-1 tolerant microarray ----------
    p = DATA / "GSE297630_processed_data.xlsx"
    wb = openpyxl.load_workbook(p, read_only=True, data_only=True)
    sh = wb[wb.sheetnames[0]]
    rows = [["" if c is None else str(c) for c in r] for r in sh.iter_rows(values_only=True)]
    wb.close()
    header_idx = next(i for i, r in enumerate(rows) if r and r[0] == "ID")
    h = rows[header_idx]
    body = rows[header_idx + 1 :]
    sym_i = h.index("Gene Symbol") if "Gene Symbol" in h else 12
    c_cols = [i for i, x in enumerate(h) if str(x).startswith("C-") and "Signal" in str(x)]
    p_cols = [i for i, x in enumerate(h) if str(x).startswith("P-") and "Signal" in str(x)]
    found297 = {}
    for row in body:
        sym = row[sym_i] if sym_i < len(row) else ""
        if gene_match(sym, TACSTD2):
            found297["Tacstd2"] = row
        elif gene_match(sym, CLDN4):
            found297["Cldn4"] = row
    catalog.append({"accession": "GSE297630", "model": "LLC aPD-1-tolerant microarray", "genes_found": {k: k in found297 for k in ("Tacstd2", "Cldn4")}})
    for gene, row in found297.items():
        ctrl = [float(row[i]) for i in c_cols]
        tol = [float(row[i]) for i in p_cols]
        for i, v in zip(c_cols, ctrl):
            dump_samples("GSE297630", "LLC", gene, [(h[i], "control", v)])
        for i, v in zip(p_cols, tol):
            dump_samples("GSE297630", "LLC", gene, [(h[i], "aPD1_tolerant", v)])
        add_contrast(
            "GSE297630",
            "LLC subcutaneous survivors",
            "anti-PD-1-tolerant remaining cells vs untreated control",
            gene,
            "control",
            "aPD1_tolerant",
            ctrl,
            tol,
            "Clariom S log2 signal; not R vs NR mice. Author C-vs-P DE also deposited.",
            icb=True,
            response=True,
        )

    # ---------- GSE169194 KPM total tumor A2V±PD1 ----------
    p = DATA / "GSE169194_annotated_log2_norm_count.txt.gz"
    h, rows, _ = read_delim(p)
    meta = parse_series_matrix(DATA / "GSE169194_series_matrix.txt.gz")
    # map 159_xx from titles
    samp = {}
    for s in meta["samples"]:
        m = re.search(r"\[(\d+_\d+)\]", s["title"])
        key = m.group(1) if m else None
        cell = next((c.split(":", 1)[-1].strip() for c in s["chars"] if c.lower().startswith("cell type")), "")
        tr = next((c.split(":", 1)[-1].strip() for c in s["chars"] if c.lower().startswith("treatment")), "")
        grp = next((c.split(":", 1)[-1].strip() for c in s["chars"] if c.lower().startswith("group")), "")
        samp[key] = {"cell": cell, "treatment": tr, "group": grp, "title": s["title"]}
    genes = {}
    for name, targets in (("Tacstd2", TACSTD2), ("Cldn4", CLDN4), ("Cd8a", CD8A), ("Cd274", CD274)):
        raw, vals = find_gene_row(h, rows, targets)
        genes[name] = (h, vals)
    catalog.append({"accession": "GSE169194", "model": "KPM GEMM", "genes_found": {k: genes[k][1] is not None for k in genes}})
    for gene, (header, vals) in genes.items():
        if vals is None:
            continue
        offset = len(header) - len(vals)
        groups = {}
        for j, col in enumerate(header[offset:]):
            m = re.search(r"(\d+_\d+)", col)
            key = m.group(1) if m else None
            info = samp.get(key, {})
            if info.get("cell", "").lower() != "total":
                continue
            g = info.get("group") or info.get("treatment") or col
            groups.setdefault(g, []).append(vals[j])
            dump_samples("GSE169194", "KPM total tumor", gene, [(col, g, vals[j])])
        # groups like total_IgG, total_A2V, total_A2V_PD1
        igg = groups.get("total_IgG", []) or groups.get("IgG", [])
        a2v = [v for k, vs in groups.items() if "A2V" in k and "PD" not in k.upper() for v in vs]
        combo = [v for k, vs in groups.items() if "PD" in k.upper() for v in vs]
        add_contrast("GSE169194", "KPM GEMM total tumor cells", "A2V+anti-PD1 vs IgG", gene, "IgG", "A2V+aPD1", igg, combo, "no aPD-1 monotherapy; anti-angiogenic combo", icb=True)
        add_contrast("GSE169194", "KPM GEMM total tumor cells", "A2V+anti-PD1 vs A2V", gene, "A2V", "A2V+aPD1", a2v, combo, "ICB on anti-angiogenic background", icb=True)

    # ---------- GSE246922 ICB-relapse LLC/KP CD45- ----------
    for acc_file, model in (("GSE246922_LLC1_RNA_counts_vst.csv.gz", "LLC1 CD45-"), ("GSE246922_KP_RNA_counts_vst.csv.gz", "KP CD45-")):
        p = DATA / acc_file
        h, rows, _ = read_delim(p)
        genes = {}
        for name, targets in (("Tacstd2", TACSTD2), ("Cldn4", CLDN4), ("Cd8a", CD8A), ("Cd274", CD274)):
            raw, vals = find_gene_row(h, rows, targets)
            genes[name] = (h, vals)
        catalog.append({"accession": "GSE246922", "model": model, "file": acc_file, "genes_found": {k: genes[k][1] is not None for k in genes}})
        for gene, (header, vals) in genes.items():
            if vals is None:
                continue
            offset = len(header) - len(vals)
            groups = {}
            for j, col in enumerate(header[offset:]):
                cl = col.lower()
                if cl.startswith("llc1y") or cl.startswith("kpy"):
                    g = "chronic_IFNG"
                elif cl.startswith("res"):
                    g = "ICB_relapse"
                elif cl.startswith("llc1") or cl.startswith("kp_"):
                    g = "parental"
                else:
                    g = col
                groups.setdefault(g, []).append(vals[j])
                dump_samples("GSE246922", model, gene, [(col, g, vals[j])])
            add_contrast(
                "GSE246922",
                model,
                "ICB-relapsed vs parental CD45- tumor cells",
                gene,
                "parental",
                "ICB_relapse",
                groups.get("parental", []),
                groups.get("ICB_relapse", []),
                "acquired-resistance model; not labeled R vs NR mice",
                icb=True,
                response=True,
            )

    # ---------- GSE179500 LKB1 restore (cite) ----------
    p = DATA / "GSE179500_20210620_XTR_Bulk_RNAseq_DESeq2norm_Counts.txt.gz"
    h, rows, _ = read_delim(p)
    meta = parse_series_matrix(DATA / "GSE179500_series_matrix.txt.gz")
    cohort = {s["title"]: next((c.split(":", 1)[-1].strip() for c in s["chars"] if c.lower().startswith("cohort")), "") for s in meta["samples"]}
    genes = {}
    for name, targets in (("Tacstd2", TACSTD2), ("Cldn4", CLDN4), ("Stk11", STK11), ("Cd8a", CD8A)):
        raw, vals = find_gene_row(h, rows, targets)
        genes[name] = dict(zip(h, vals)) if vals is not None else None
    catalog.append({"accession": "GSE179500", "model": "KL XTR GEMM neoplastic cells", "genes_found": {k: v is not None for k, v in genes.items()}, "note": "known; LKB1 off vs restore; not ICB"})
    for gene, d in genes.items():
        if not d:
            continue
        off, on, wt = [], [], []
        for title, coh in cohort.items():
            # columns are sample titles
            val = None
            for k, v in d.items():
                if norm_key(k) == norm_key(title) or norm_key(title) in norm_key(k):
                    val = v
                    break
            if val is None:
                continue
            dump_samples("GSE179500", "KL XTR", gene, [(title, coh, val)])
            cl = coh.lower()
            if "non-restored" in cl or "nonrestored" in cl:
                off.append(val)
            elif cl.startswith("restored"):
                on.append(val)
            elif "wild" in cl:
                wt.append(val)
        add_contrast("GSE179500", "KL XTR FACS tumor cells", "Lkb1 non-restored vs restored", gene, "Lkb1_restored", "Lkb1_off", on, off, "KNOWN series; genotype not ICB. Tacstd2 expected higher when LKB1 off", genotype=True)
        add_contrast("GSE179500", "KL XTR FACS tumor cells", "Lkb1 non-restored vs KT wild-type", gene, "KT_WT", "Lkb1_off", wt, off, "genotype", genotype=True)

    # author DESeq2 row
    p = DATA / "GSE179500_20210211_XTR_Bulk_RNAseq_DESeq2_results_RestoredvNon-Restored.csv.gz"
    h, rows, _ = read_delim(p)
    for name, targets in (("Tacstd2", TACSTD2), ("Cldn4", CLDN4), ("Stk11", STK11)):
        raw, _ = find_gene_row(h, rows, targets)
        if raw:
            # columns: baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, Gene
            rec = {h[i] if i < len(h) else str(i): raw[i] if i < len(raw) else "" for i in range(len(raw))}
            Path(OUT / f"GSE179500_deseq2_{name}.json").write_text(json.dumps({"gene": name, "row": raw[:12], "header": h[:12]}, indent=2))

    # ---------- GSE137396 KL vs KP nodules ----------
    p = DATA / "GSE137396_Normalized_logtransformed_medcentred_genetable_GEMMnodule.txt.gz"
    h, rows, _ = read_delim(p)
    genes = {}
    for name, targets in (("Tacstd2", TACSTD2), ("Cldn4", CLDN4), ("Stk11", STK11), ("Cd8a", CD8A)):
        raw, vals = find_gene_row(h, rows, targets)
        genes[name] = (h, vals)
    catalog.append({"accession": "GSE137396", "model": "KL vs KP GEMM nodules", "genes_found": {k: genes[k][1] is not None for k in genes}})
    for gene, (header, vals) in genes.items():
        if vals is None:
            continue
        offset = len(header) - len(vals)
        kl, kp = [], []
        for j, col in enumerate(header[offset:]):
            if col.startswith("KL") or "Lkb1" in col or "LKB1" in col:
                kl.append(vals[j])
                g = "KL"
            else:
                kp.append(vals[j])
                g = "KP"
            dump_samples("GSE137396", "GEMM nodule", gene, [(col, g, vals[j])])
        add_contrast("GSE137396", "GEMM lung nodules", "KL vs KP", gene, "KP", "KL", kp, kl, "already-log median-centered; not ICB", genotype=True)

    # ---------- GSE244452 KP vs KL DEG ----------
    p = DATA / "GSE244452_KPvsKL_deg_all.txt.gz"
    h, rows, _ = read_delim(p)
    for name, targets in (("Tacstd2", TACSTD2), ("Cldn4", CLDN4), ("Stk11", STK11), ("Cd8a", CD8A)):
        raw, _ = find_gene_row(h, rows, targets)
        if raw:
            Path(OUT / f"GSE244452_deg_{name}.json").write_text(json.dumps({"gene": name, "header": h[:15], "row": raw[:15]}, indent=2))
            # try to interpret logFC KP vs KL
            rec = {str(h[i]).lower(): raw[i] for i in range(min(len(h), len(raw)))}
            logfc = None
            pval = None
            for k, v in rec.items():
                if "log" in k and "fc" in k.replace("fold", "fc"):
                    try:
                        logfc = float(v)
                    except Exception:
                        pass
                if k in ("pvalue", "p-value", "pval", "p.value") or k == "pvalue":
                    try:
                        pval = float(v)
                    except Exception:
                        pass
            contrasts.append(
                {
                    "accession": "GSE244452",
                    "model": "KL vs KP subcutaneous",
                    "design": "author DEG KP vs KL (no sample matrix)",
                    "gene": name,
                    "group_a": "KP",
                    "group_b": "KL",
                    "icb_contrast": False,
                    "response_contrast": False,
                    "genotype_contrast": True,
                    "note": "DEG table only; RNA is genotype not ICB",
                    "n_a": 3,
                    "n_b": 3,
                    "log2fc_mean": logfc,
                    "welch_p": pval,
                    "mwu_p": None,
                    "author_row": rec,
                }
            )
    catalog.append({"accession": "GSE244452", "model": "KL vs KP tumors", "note": "DEG only; no ICB RNA"})

    # ---------- GSE274351 K/KP/KL LCM ----------
    p = DATA / "GSE274351_expredata_TPM_gene.txt.gz"
    h, rows, _ = read_delim(p)
    meta = parse_series_matrix(DATA / "GSE274351_series_matrix.txt.gz")
    genes = {}
    for name, targets in (("Tacstd2", TACSTD2), ("Cldn4", CLDN4), ("Stk11", STK11), ("Cd8a", CD8A)):
        raw, vals = find_gene_row(h, rows, targets)
        genes[name] = (h, vals)
    catalog.append({"accession": "GSE274351", "model": "K/KP/KL LCM adenomas", "genes_found": {k: genes[k][1] is not None for k in genes}})
    for gene, (header, vals) in genes.items():
        if vals is None:
            continue
        offset = len(header) - len(vals)
        groups = {"K": [], "KP": [], "KL": [], "NL": []}
        for j, col in enumerate(header[offset:]):
            cl = col.lower()
            if cl.startswith("kl"):
                g = "KL"
            elif cl.startswith("kp"):
                g = "KP"
            elif cl.startswith("nl"):
                g = "NL"
            elif cl.startswith("k") and not cl.startswith("kl") and not cl.startswith("kp"):
                g = "K"
            else:
                g = "K"
            groups[g].append(vals[j])
            dump_samples("GSE274351", "LCM adenoma", gene, [(col, g, vals[j])])
        add_contrast("GSE274351", "KRAS GEMM LCM adenomas", "KL vs KP", gene, "KP", "KL", groups["KP"], groups["KL"], "no ICB on these samples", genotype=True)
        add_contrast("GSE274351", "KRAS GEMM LCM adenomas", "KL vs K", gene, "K", "KL", groups["K"], groups["KL"], "no ICB", genotype=True)

    # ---------- GSE338923 STK11 KO Lacun3 ----------
    p = DATA / "GSE338923_Lacun3_STK11_RNAseq_voom_normalized_counts.tsv.gz"
    h, rows, _ = read_delim(p)
    genes = {}
    for name, targets in (("Tacstd2", TACSTD2), ("Cldn4", CLDN4), ("Stk11", STK11)):
        raw, vals = find_gene_row(h, rows, targets)
        genes[name] = (h, vals)
    catalog.append({"accession": "GSE338923", "model": "Lacun3 STK11 KO vitro", "genes_found": {k: genes[k][1] is not None for k in genes}})
    for gene, (header, vals) in genes.items():
        if vals is None:
            continue
        offset = len(header) - len(vals)
        wt, ko = [], []
        for j, col in enumerate(header[offset:]):
            cl = str(col)
            if cl.startswith("C"):
                wt.append(vals[j])
                g = "STK11_WT"
            elif cl.startswith("S"):
                ko.append(vals[j])
                g = "STK11_KO"
            else:
                g = col
            dump_samples("GSE338923", "Lacun3 vitro", gene, [(col, g, vals[j])])
        add_contrast("GSE338923", "Lacun3 LUAD cells", "STK11 KO vs WT", gene, "STK11_WT", "STK11_KO", wt, ko, "in vitro; no ICB", genotype=True)

    # ---------- GSE330941 Ago2KO LLC ICI-sensitized proxy ----------
    p = DATA / "GSE330941_filtered_tablecounts_tpm.csv.gz"
    h, rows, _ = read_delim(p)
    genes = {}
    for name, targets in (("Tacstd2", TACSTD2), ("Cldn4", CLDN4), ("Cd8a", CD8A), ("Cd274", CD274)):
        raw, vals = find_gene_row(h, rows, targets)
        genes[name] = (h, vals)
    catalog.append({"accession": "GSE330941", "model": "LLC Ago2KO vs WT (no ICI on RNA)", "genes_found": {k: genes[k][1] is not None for k in genes}})
    for gene, (header, vals) in genes.items():
        if vals is None:
            continue
        offset = len(header) - len(vals)
        # series matrix order: first 4 WT, last 4 Ago2KO (D1727T149-152 vs T157-160)
        wt_cols = {"D1727T149", "D1727T150", "D1727T151", "D1727T152"}
        ko_cols = {"D1727T157", "D1727T158", "D1727T159", "D1727T160"}
        wt, ko = [], []
        for j, col in enumerate(header[offset:]):
            c = str(col).strip('"')
            if c in wt_cols:
                wt.append(vals[j])
                g = "WT"
            elif c in ko_cols:
                ko.append(vals[j])
                g = "Ago2KO"
            else:
                g = c
            dump_samples("GSE330941", "LLC", gene, [(c, g, vals[j])])
        add_contrast("GSE330941", "LLC subcutaneous d12", "Ago2KO (ICI-sensitized) vs WT (ICI-refractory)", gene, "WT", "Ago2KO", wt, ko, "RNA without ICI; model-level sensitivity proxy", response=True)

    # ---------- GSE260596 344SQ all aPD-1 ----------
    p = DATA / "GSE260596_AllSamples_Genes_ReadCounts.txt.gz"
    h, rows, _ = read_delim(p)
    genes = {}
    for name, targets in (("Tacstd2", TACSTD2), ("Cldn4", CLDN4), ("Cd8a", CD8A), ("Cd274", CD274)):
        raw, vals = find_gene_row(h, rows, targets)
        genes[name] = (h, vals)
    catalog.append({"accession": "GSE260596", "model": "344SQ + aPD-1 ± LAIR1", "genes_found": {k: genes[k][1] is not None for k in genes}, "note": "all arms have aPD-1; no untreated"})
    for gene, (header, vals) in genes.items():
        if vals is None:
            continue
        offset = len(header) - len(vals)
        klh, lair = [], []
        for j, col in enumerate(header[offset:]):
            cl = col.lower()
            if "lair" in cl:
                lair.append(vals[j])
                g = "aPD1+aLAIR1"
            else:
                klh.append(vals[j])
                g = "aPD1+aKLH"
            dump_samples("GSE260596", "344SQ", gene, [(col, g, vals[j])])
        add_contrast("GSE260596", "344SQ KP-derived", "aPD-1+aLAIR1 vs aPD-1+control", gene, "aPD1+aKLH", "aPD1+aLAIR1", klh, lair, "no ICB-untreated arm", icb=True)

    # ---------- GSE197260 EGFR n=1 ----------
    p = DATA / "GSE197260_RNAseqTPM_MM_EGFR-TKI-CD8.txt.gz"
    h, rows, _ = read_delim(p)
    genes = {}
    for name, targets in (("Tacstd2", TACSTD2), ("Cldn4", CLDN4), ("Cd8a", CD8A), ("Cd274", CD274)):
        raw, vals = find_gene_row(h, rows, targets)
        genes[name] = (h, vals)
    catalog.append({"accession": "GSE197260", "model": "EGFR-mutant NSCLC", "genes_found": {k: genes[k][1] is not None for k in genes}, "note": "n=1/arm"})
    for gene, (header, vals) in genes.items():
        if vals is None:
            continue
        d = dict(zip(header[1:], vals))
        for k, v in d.items():
            dump_samples("GSE197260", "EGFR NSCLC", gene, [(k, k, v)])
        a = [d[k] for k in d if "gef_vehicle" in k]
        b = [d[k] for k in d if "gef_4h2" in k or "4h2" in k]
        add_contrast("GSE197260", "EGFR-mutant syngeneic", "gefitinib then aPD-1 (4H2) vs gefitinib then vehicle", gene, "gef_vehicle_d21", "gef_4h2_d21", a, b, "n=1 per arm; p not computed", icb=True)

    # ---------- GSE256071 SOX2 KP cells ----------
    p = DATA / "GSE256071_210602GerA_l2tpm.txt.gz"
    h, rows, _ = read_delim(p)
    genes = {}
    for name, targets in (("Tacstd2", TACSTD2), ("Cldn4", CLDN4)):
        raw, vals = find_gene_row(h, rows, targets)
        genes[name] = (h, vals)
    catalog.append({"accession": "GSE256071", "model": "KP ± SOX2 cells", "genes_found": {k: genes[k][1] is not None for k in genes}})
    for gene, (header, vals) in genes.items():
        if vals is None:
            continue
        offset = len(header) - len(vals)
        ct, sox = [], []
        for j, col in enumerate(header[offset:]):
            if "KPS2" in col or "KPS2" in str(col).upper():
                sox.append(vals[j])
                g = "SOX2"
            elif "KPCt" in col:
                ct.append(vals[j])
                g = "KP_ctrl"
            else:
                g = col
            dump_samples("GSE256071", "KP cells", gene, [(col, g, vals[j])])
        add_contrast("GSE256071", "KP cells vitro", "SOX2-high vs control", gene, "KP_ctrl", "SOX2", ct, sox, "checkpoint-resistance model; no ICB on RNA")

    # ---------- GSE182228 missing matrix ----------
    catalog.append(
        {
            "accession": "GSE182228",
            "model": "LKB1-deficient LUAD + aPD-1 ± palbociclib",
            "n_samples": 12,
            "note": "ICB design exists (vehicle / palbo / aPD-1 / combo, n=3) but no processed gene matrix on GEO (series matrix is metadata-only, 2.8 KB). Not scored.",
        }
    )
    catalog.append(
        {
            "accession": "GSE155972",
            "model": "TISMO LLC ICB",
            "note": "EXCLUDED: TISMO source. Taken as given (49/64 Tacstd2 up, p=5.8e-5) and not re-audited.",
        }
    )
    catalog.append(
        {
            "accession": "GSE298051",
            "model": "KL/KP TNG260 ± aPD-1",
            "note": "Mostly H3K9ac ChIP; last 8 samples are subcutaneous tumors including TNG260+aPD1. No processed RNA matrix deposited.",
        }
    )
    catalog.append(
        {
            "accession": "GSE114300",
            "model": "Kras GEMM TIL CD4/CD8 ± aPD-1",
            "note": "Sorted lymphocytes only; Tacstd2/Cldn4 not expected. No tumor epithelial matrix.",
        }
    )

    # FDR within primary ICB Tacstd2/Cldn4 welch tests
    primary = [c for c in contrasts if c.get("icb_contrast") and c.get("welch_p") is not None]
    qs = bh_fdr([c["welch_p"] for c in primary])
    for c, q in zip(primary, qs):
        c["welch_fdr_among_icb"] = q

    # write tables
    def write_tsv(path, rows):
        if not rows:
            path.write_text("")
            return
        keys = []
        for r in rows:
            for k in r:
                if k not in keys:
                    keys.append(k)
        with path.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore", delimiter="\t")
            w.writeheader()
            for r in rows:
                w.writerow({k: ("" if r.get(k) is None else r.get(k)) for k in keys})

    write_tsv(OUT / "contrasts.tsv", contrasts)
    write_tsv(OUT / "per_sample.tsv", per_sample)
    (OUT / "catalog.json").write_text(json.dumps(catalog, indent=2))
    (OUT / "contrasts.json").write_text(json.dumps(contrasts, indent=2, default=str))

    # summary counts
    icb_t = [c for c in contrasts if c.get("icb_contrast") and c["gene"] == "Tacstd2" and c.get("n_a", 0) >= 2 and c.get("n_b", 0) >= 2]
    print("ICB Tacstd2 contrasts with n>=2:", len(icb_t))
    for c in contrasts:
        if c["gene"] in ("Tacstd2", "Cldn4") and (c.get("icb_contrast") or c.get("genotype_contrast") or c.get("response_contrast")):
            print(
                f"{c['accession']}\t{c['gene']}\t{c['design']}\tn={c.get('n_a')}/{c.get('n_b')}\t"
                f"log2FC={c.get('log2fc_mean')}\twelch={c.get('welch_p')}\tmwu={c.get('mwu_p')}"
            )


if __name__ == "__main__":
    main()
