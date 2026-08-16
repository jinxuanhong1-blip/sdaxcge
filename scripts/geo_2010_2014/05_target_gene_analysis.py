"""TACSTD2 / CLDN4 / CD274 / PDCD1 audit + leftover analyses.

Honest scope of 2010-2014 leftover:
  * No human lung PD-1/PD-L1/CTLA-4 ICI treatment-response series exist in
    this window (nivolumab NSCLC approval is 2015).
  * Human lung immunotherapy leftovers are cancer-vaccine / epitope studies
    (GSE35640 MAGE-A3; GSE27556 HLA ligands). Analyze TACSTD2/CLDN4 there
    if the platform actually measures them.
  * Mouse GSE54351/2/3 is Lkb1/Pten lung SCC with elevated Pdl1 — biology,
    not ICI treatment. Cd274 is reported if present; TACSTD2/CLDN4 only if
    the mouse array annotation carries them.

Nothing is imputed. Missing genes are recorded as not measured.
"""
import glob
import gzip
import json
import math
import os
import re

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
import sys

sys.path.insert(0, HERE)
from eutils import _get  # noqa: E402

ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
OUT = os.path.join(ROOT, "results", "w200", "GEO_2010_2014")
DATA = os.path.join(OUT, "data")
CLIN = os.path.join(OUT, "clinical")
ANA = os.path.join(OUT, "analysis")
os.makedirs(ANA, exist_ok=True)

TARGET = {
    "TACSTD2": ["TACSTD2", "TROP2", "TROP-2", "GA733-1", "M1S1"],
    "CLDN4": ["CLDN4", "Cldn4"],
    "CD274": ["CD274", "PD-L1", "PDL1", "B7-H1", "Cd274", "Pdl1"],
    "PDCD1": ["PDCD1", "PD-1", "PD1", "Pdcd1"],
}


def _unq(c):
    return c.strip().strip('"')


def _norm(tok):
    return re.sub(r"[^A-Za-z0-9]", "", tok).upper()


ALIAS = {}
for g, als in TARGET.items():
    for a in [g, *als]:
        ALIAS[_norm(a)] = g


def fetch_platform_symbol_map(gpl):
    """Real SOFT platform table -> probe -> gene for TARGET genes only."""
    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gpl}&targ=self&form=text&view=full"
    txt = _get(url, timeout=180)
    header = []
    rows = []
    in_table = False
    for ln in txt.splitlines():
        if ln.startswith("!platform_table_begin"):
            in_table = True
            continue
        if ln.startswith("!platform_table_end"):
            break
        if in_table:
            cells = ln.split("\t")
            if not header:
                header = cells
            else:
                rows.append(cells)
    if not header:
        return {}, None, 0
    cand = [
        "Gene Symbol",
        "Gene symbol",
        "GENE_SYMBOL",
        "Symbol",
        "SYMBOL",
        "ILMN_Gene",
        "GeneSymbol",
        "gene_assignment",
        "gene",
    ]
    sym_idx = None
    sym_name = None
    for c in cand:
        for i, h in enumerate(header):
            if h.strip() == c:
                sym_idx, sym_name = i, h.strip()
                break
        if sym_idx is not None:
            break
    gene_probes = {g: [] for g in TARGET}
    if sym_idx is None:
        return gene_probes, None, len(rows)
    for cells in rows:
        if not cells:
            continue
        probe = cells[0].strip()
        if not probe or len(cells) <= sym_idx:
            continue
        raw_sym = cells[sym_idx]
        # Affy gene_assignment is "ENST // Symbol // description // ..."
        tokens = re.split(r"[ /,;|]+", raw_sym.replace(" // ", " "))
        for tok in tokens:
            g = ALIAS.get(_norm(tok))
            if g:
                gene_probes[g].append(probe)
                break
    return gene_probes, sym_name, len(rows)


def parse_series_matrix(path):
    opener = gzip.open if path.endswith(".gz") else open
    samples = []
    platform = ""
    in_table = False
    header = []
    values = {}
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for ln in fh:
            if ln.startswith("!Series_platform_id"):
                parts = ln.split("\t", 1)
                if len(parts) > 1:
                    platform = _unq(parts[1])
            elif ln.startswith("!Sample_geo_accession"):
                samples = [_unq(c) for c in ln.split("\t")[1:]]
            elif ln.startswith("!series_matrix_table_begin"):
                in_table = True
            elif ln.startswith("!series_matrix_table_end"):
                break
            elif in_table:
                cells = ln.rstrip("\n").split("\t")
                if not header:
                    header = [_unq(c) for c in cells]
                else:
                    values[_unq(cells[0])] = cells[1:]
    return {"platform": platform, "samples": samples, "values": values}


def to_float(c):
    c = _unq(c)
    if c == "" or c.upper() in {"NA", "NAN", "NULL"}:
        return math.nan
    try:
        return float(c)
    except ValueError:
        return math.nan


def collapse(parsed, gene_probes):
    n = len(parsed["samples"])
    out = {}
    used = {}
    for g, probes in gene_probes.items():
        present = [p for p in probes if p in parsed["values"]]
        used[g] = present
        if not present:
            out[g] = np.full(n, math.nan)
            continue
        stack = []
        for p in present:
            vals = parsed["values"][p]
            arr = np.array(
                [to_float(v) for v in vals[:n]] + [math.nan] * max(0, n - len(vals))
            )
            stack.append(arr[:n])
        with np.errstate(all="ignore"):
            out[g] = np.nanmean(np.vstack(stack), axis=0)
    return out, used


def summ(arr):
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        return {"n": 0, "mean": None, "median": None, "sd": None, "min": None, "max": None}
    return {
        "n": int(finite.size),
        "mean": round(float(np.mean(finite)), 4),
        "median": round(float(np.median(finite)), 4),
        "sd": round(float(np.std(finite, ddof=1)) if finite.size > 1 else 0.0, 4),
        "min": round(float(np.min(finite)), 4),
        "max": round(float(np.max(finite)), 4),
    }


def pearson(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    xs, ys = x[m], y[m]
    if xs.size < 3 or np.std(xs) == 0 or np.std(ys) == 0:
        return None, int(xs.size)
    return round(float(np.corrcoef(xs, ys)[0, 1]), 4), int(xs.size)


def spearman(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    xs, ys = x[m], y[m]
    if xs.size < 3 or np.std(xs) == 0 or np.std(ys) == 0:
        return None, int(xs.size)
    rx = np.argsort(np.argsort(xs)).astype(float)
    ry = np.argsort(np.argsort(ys)).astype(float)
    return round(float(np.corrcoef(rx, ry)[0, 1]), 4), int(xs.size)


def mannwhitney(a, b):
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if a.size < 3 or b.size < 3:
        return {"n_a": int(a.size), "n_b": int(b.size), "u": None, "p": None, "median_a": None, "median_b": None}
    # SciPy-free Mann-Whitney U (two-sided, no tie correction beyond midranks).
    comb = np.concatenate([a, b])
    ranks = np.argsort(np.argsort(comb)).astype(float) + 1.0
    # average ties
    order = np.argsort(comb)
    sorted_v = comb[order]
    i = 0
    avg = ranks.copy()
    while i < len(sorted_v):
        j = i
        while j < len(sorted_v) and sorted_v[j] == sorted_v[i]:
            j += 1
        if j - i > 1:
            mid = 0.5 * (i + 1 + j)
            avg[order[i:j]] = mid
        i = j
    u = float(avg[: a.size].sum() - a.size * (a.size + 1) / 2.0)
    mu = a.size * b.size / 2.0
    sigma = math.sqrt(a.size * b.size * (a.size + b.size + 1) / 12.0)
    if sigma == 0:
        p = 1.0
    else:
        z = (u - mu) / sigma
        # two-sided from erfc
        p = math.erfc(abs(z) / math.sqrt(2.0))
    return {
        "n_a": int(a.size),
        "n_b": int(b.size),
        "u": round(u, 3),
        "p": round(p, 4),
        "median_a": round(float(np.median(a)), 4),
        "median_b": round(float(np.median(b)), 4),
    }


def read_clinical(acc):
    path = os.path.join(CLIN, f"{acc}_clinical.tsv")
    if not os.path.exists(path):
        return []
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        return [dict(zip(header, ln.rstrip("\n").split("\t"))) for ln in f]


def main():
    annot_cache = {}
    audit = []
    analyses = {}
    matrices = sorted(glob.glob(os.path.join(DATA, "*", "*series_matrix*.txt.gz")))
    matrices += sorted(glob.glob(os.path.join(DATA, "*", "*series_matrix*.txt")))
    for path in matrices:
        acc = os.path.basename(os.path.dirname(path))
        fname = os.path.basename(path)
        print(f"analyze {acc}/{fname}", flush=True)
        try:
            parsed = parse_series_matrix(path)
        except Exception as exc:  # noqa: BLE001
            audit.append(
                {
                    "accession": acc,
                    "matrix_file": fname,
                    "status": f"parse_failed:{exc}",
                }
            )
            continue
        gpl = parsed["platform"]
        m = re.search(r"-(GPL\d+)_series_matrix", fname)
        if m:
            # Multi-platform series matrices: filename platform is authoritative.
            gpl = m.group(1)
        if not gpl:
            m = re.search(r"(GPL\d+)", fname)
            gpl = m.group(1) if m else ""
        if not gpl:
            audit.append(
                {
                    "accession": acc,
                    "matrix_file": fname,
                    "status": "no_platform",
                    "n_samples": len(parsed["samples"]),
                }
            )
            continue
        if gpl not in annot_cache:
            print(f"  platform annot {gpl}", flush=True)
            annot_cache[gpl] = fetch_platform_symbol_map(gpl)
        gene_probes, sym_col, n_probes = annot_cache[gpl]
        gm, used = collapse(parsed, gene_probes)
        summaries = {g: summ(gm[g]) for g in TARGET}
        pairs = {}
        for a, b in [
            ("TACSTD2", "CD274"),
            ("TACSTD2", "PDCD1"),
            ("CLDN4", "CD274"),
            ("CLDN4", "PDCD1"),
            ("TACSTD2", "CLDN4"),
            ("CD274", "PDCD1"),
        ]:
            pr, npr = pearson(gm[a], gm[b])
            sr, _ = spearman(gm[a], gm[b])
            pairs[f"{a}~{b}"] = {"pearson": pr, "spearman": sr, "n": npr}
        rec = {
            "accession": acc,
            "matrix_file": fname,
            "status": "ok",
            "platform": gpl,
            "symbol_col": sym_col,
            "n_platform_rows": n_probes,
            "n_samples": len(parsed["samples"]),
            "gene_probes": used,
            "gene_summary": summaries,
            "coexpression": pairs,
            "all_four_measured": all(summaries[g]["n"] > 0 for g in TARGET),
        }
        audit.append(
            {
                "accession": acc,
                "matrix_file": fname,
                "platform": gpl,
                "symbol_col": sym_col or "",
                "n_samples": len(parsed["samples"]),
                "TACSTD2_n": summaries["TACSTD2"]["n"],
                "CLDN4_n": summaries["CLDN4"]["n"],
                "CD274_n": summaries["CD274"]["n"],
                "PDCD1_n": summaries["PDCD1"]["n"],
                "TACSTD2_probes": ";".join(used["TACSTD2"]),
                "CLDN4_probes": ";".join(used["CLDN4"]),
                "CD274_probes": ";".join(used["CD274"]),
                "PDCD1_probes": ";".join(used["PDCD1"]),
            }
        )
        analyses[f"{acc}:{fname}"] = rec

        # Optional response contrast if a binary-ish clinical field exists.
        clin = read_clinical(acc)
        if clin and summaries["TACSTD2"]["n"] > 0:
            by_gsm = {r.get("geo_accession", ""): r for r in clin}
            # Look for a field that looks like response / histology / tissue.
            contrast_fields = []
            for key in clin[0]:
                if key in ("geo_accession", "title", "source_name", "matrix_file", "platform_from_filename", "other_characteristics"):
                    continue
                vals = [r.get(key, "") for r in clin]
                uniq = sorted({v for v in vals if v})
                if 2 <= len(uniq) <= 8:
                    contrast_fields.append((key, uniq))
            contrasts = {}
            for key, uniq in contrast_fields:
                # Only report TACSTD2/CLDN4 medians per level (no invented grouping).
                per = {}
                for level in uniq:
                    mask = np.array(
                        [
                            by_gsm.get(s, {}).get(key, "") == level
                            for s in parsed["samples"]
                        ]
                    )
                    per[level] = {
                        "n": int(mask.sum()),
                        "TACSTD2_median": None
                        if gm["TACSTD2"][mask].size == 0 or not np.isfinite(gm["TACSTD2"][mask]).any()
                        else round(float(np.nanmedian(gm["TACSTD2"][mask])), 4),
                        "CLDN4_median": None
                        if gm["CLDN4"][mask].size == 0 or not np.isfinite(gm["CLDN4"][mask]).any()
                        else round(float(np.nanmedian(gm["CLDN4"][mask])), 4),
                        "CD274_median": None
                        if gm["CD274"][mask].size == 0 or not np.isfinite(gm["CD274"][mask]).any()
                        else round(float(np.nanmedian(gm["CD274"][mask])), 4),
                    }
                contrasts[key] = per
            rec["clinical_contrasts"] = contrasts

            # If a field is clearly responder vs not, run MW.
            for key, uniq in contrast_fields:
                low = {u.lower() for u in uniq}
                pos = [
                    u
                    for u in uniq
                    if re.search(r"^(r|responder|cr|pr|yes|sensitive)$", u.strip(), re.I)
                ]
                neg = [
                    u
                    for u in uniq
                    if re.search(r"^(nr|non[- ]?responder|pd|no|resistant)$", u.strip(), re.I)
                ]
                if pos and neg:
                    mask_a = np.array(
                        [by_gsm.get(s, {}).get(key, "") in pos for s in parsed["samples"]]
                    )
                    mask_b = np.array(
                        [by_gsm.get(s, {}).get(key, "") in neg for s in parsed["samples"]]
                    )
                    rec.setdefault("response_tests", {})[key] = {
                        "positive_levels": pos,
                        "negative_levels": neg,
                        "TACSTD2": mannwhitney(gm["TACSTD2"][mask_a], gm["TACSTD2"][mask_b]),
                        "CLDN4": mannwhitney(gm["CLDN4"][mask_a], gm["CLDN4"][mask_b]),
                        "CD274": mannwhitney(gm["CD274"][mask_a], gm["CD274"][mask_b]),
                    }

    acols = [
        "accession",
        "matrix_file",
        "platform",
        "symbol_col",
        "n_samples",
        "TACSTD2_n",
        "CLDN4_n",
        "CD274_n",
        "PDCD1_n",
        "TACSTD2_probes",
        "CLDN4_probes",
        "CD274_probes",
        "PDCD1_probes",
    ]
    with open(os.path.join(ANA, "gene_availability_audit.tsv"), "w") as f:
        f.write("\t".join(acols) + "\n")
        for r in audit:
            if "TACSTD2_n" not in r:
                continue
            f.write("\t".join(str(r.get(c, "")) for c in acols) + "\n")

    # JSON-safe
    def _clean(o):
        if isinstance(o, dict):
            return {k: _clean(v) for k, v in o.items()}
        if isinstance(o, list):
            return [_clean(v) for v in o]
        if isinstance(o, float) and (math.isnan(o) or math.isinf(o)):
            return None
        return o

    with open(os.path.join(ANA, "analysis_results.json"), "w") as f:
        json.dump(_clean(analyses), f, indent=1)
    print(f"analyzed {len(analyses)} matrices -> results/w200/GEO_2010_2014/analysis/")


if __name__ == "__main__":
    main()
