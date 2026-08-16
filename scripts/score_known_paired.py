#!/usr/bin/env python3
"""Score the GEO series that actually have usable paired ICB + TACSTD2 tables.

This is the confirmatory arm: each dataset is handled from its real files (not a
one-regex-fits-all parser). Failures of the generic hunter are classified separately.
"""
import gzip
import json
import os
import re

import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(ROOT, "data", "geo", "cache")
OUT = os.path.join(ROOT, "results", "hunt_paired_up", "geo")


def wilcox_and_sign(pairs):
    n = len(pairs)
    if n == 0:
        return {
            "n_paired_subjects": 0,
            "n_up": 0,
            "n_down": 0,
            "n_flat": 0,
            "frac_up": None,
            "median_delta": None,
            "mean_pre": None,
            "mean_post": None,
            "wilcoxon_p": None,
            "sign_test_p": None,
        }
    up = int((pairs.delta > 0).sum())
    rec = {
        "n_paired_subjects": n,
        "n_up": up,
        "n_down": int((pairs.delta < 0).sum()),
        "n_flat": int((pairs.delta == 0).sum()),
        "frac_up": round(up / n, 3) if n else None,
        "median_delta": float(pairs.delta.median()) if n else None,
        "mean_pre": float(pairs.pre.mean()) if n else None,
        "mean_post": float(pairs.post.mean()) if n else None,
        "wilcoxon_p": (
            float(stats.wilcoxon(pairs.post, pairs.pre, alternative="two-sided", zero_method="wilcox").pvalue)
            if n >= 6 and (pairs.delta != 0).sum() >= 6
            else None
        ),
        "sign_test_p": float(stats.binomtest(up, n, 0.5).pvalue) if n else None,
    }
    return rec


def parse_matrix_chars(path):
    opener = gzip.open if path.endswith(".gz") else open
    blocks = {}
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.startswith("!"):
                continue
            key, _, rest = line[1:].partition("\t")
            key = key.strip()
            vals = [v.strip().strip('"') for v in rest.rstrip("\n").split("\t")]
            blocks.setdefault(key, []).append(vals)
    samples = blocks.get("Sample_geo_accession", [[]])[0]
    titles = dict(zip(samples, blocks.get("Sample_title", [[]])[0])) if samples else {}
    chars = {s: {} for s in samples}
    for vals in blocks.get("Sample_characteristics_ch1", []):
        for gsm, val in zip(samples, vals):
            if ":" in val:
                k, v = val.split(":", 1)
                chars[gsm][k.strip().lower()] = v.strip()
    return samples, titles, chars


def score_gse91061():
    fpkm = pd.read_csv(os.path.join(CACHE, "GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv.gz"), index_col=0)
    rld = pd.read_csv(os.path.join(CACHE, "GSE91061_BMS038109Sample.hg19KnownGene.rld.csv.gz"), index_col=0)
    out = {}
    for label, df in (("fpkm", fpkm), ("rld", rld)):
        row = df.loc[4070]
        pre, post = {}, {}
        for c in df.columns:
            parts = str(c).split("_")
            if len(parts) < 2:
                continue
            pt, tp = parts[0].lower(), parts[1].lower()
            if tp == "pre":
                pre.setdefault(pt, []).append(c)
            elif tp == "on":
                post.setdefault(pt, []).append(c)
        rows = []
        for pt in sorted(set(pre) & set(post)):
            a = float(row[pre[pt]].mean())
            b = float(row[post[pt]].mean())
            rows.append({"subject": pt, "pre": a, "post": b, "delta": b - a})
        pairs = pd.DataFrame(rows)
        pairs.to_csv(os.path.join(OUT, f"GSE91061_{label}_pairs.csv"), index=False)
        rec = wilcox_and_sign(pairs)
        rec.update(
            {
                "gse": "GSE91061",
                "unit": label,
                "title": "Riaz 2017 melanoma nivolumab pre/on-treatment (CheckMate 038 / BMS038)",
                "taxon": "Homo sapiens",
                "cancer": "melanoma",
                "treatment": "anti-PD-1 (nivolumab)",
                "pairing": "same patient, pre vs on-treatment tumor RNA-seq",
                "note": "Entrez 4070 in hg19KnownGene table. Primary unit is FPKM.",
            }
        )
        out[label] = rec
    return out["fpkm"], out


def score_gse115821():
    samples, titles, chars = parse_matrix_chars(os.path.join(CACHE, "GSE115821-GPL11154_series_matrix.txt.gz"))
    counts = pd.read_csv(os.path.join(CACHE, "GSE115821_MGH_counts.csv.gz"))
    gene = counts[counts.Geneid.astype(str).str.upper() == "TACSTD2"].iloc[0]
    # map title -> count column (titles lack .bam; columns have it)
    colmap = {}
    for c in counts.columns:
        if c in ("Geneid", "Chr", "Start", "End", "Strand", "Length"):
            continue
        key = re.sub(r"\.bam$", "", c)
        colmap[key] = c
        colmap[c] = c
    rows = []
    by = {}
    for gsm in samples:
        title = titles.get(gsm, "")
        col = colmap.get(title) or colmap.get(title + ".bam")
        if col is None:
            # fuzzy
            hits = [c for k, c in colmap.items() if title and title in k]
            col = hits[0] if hits else None
        if col is None:
            continue
        pt = chars[gsm].get("patient id") or chars[gsm].get("patient")
        state = (chars[gsm].get("treatment state") or "").lower()
        if "pre" in state:
            tp = "pre"
        elif re.search(r"\bon\b", state):
            tp = "post"
        else:
            continue
        by.setdefault(pt, {"pre": [], "post": []})[tp].append(float(gene[col]))
    for pt, tps in by.items():
        if tps["pre"] and tps["post"]:
            a, b = float(pd.Series(tps["pre"]).mean()), float(pd.Series(tps["post"]).mean())
            rows.append({"subject": pt, "pre": a, "post": b, "delta": b - a})
    pairs = pd.DataFrame(rows)
    pairs.to_csv(os.path.join(OUT, "GSE115821_pairs.csv"), index=False)
    rec = wilcox_and_sign(pairs)
    rec.update(
        {
            "gse": "GSE115821",
            "unit": "raw_counts",
            "title": "Auslander 2018 MGH metastatic melanoma ICB (subset with serial biopsies)",
            "taxon": "Homo sapiens",
            "cancer": "melanoma",
            "treatment": "anti-PD-1 and/or anti-CTLA-4",
            "pairing": "same patient, PRE vs On ICB tumor RNA-seq",
            "note": "Small n. Counts are not library-size normalized; direction is still usable.",
        }
    )
    return rec


def score_gse168204():
    samples, titles, chars = parse_matrix_chars(os.path.join(CACHE, "GSE168204-GPL11154_series_matrix.txt.gz"))
    path = os.path.join(CACHE, "GSE168204_MGH_counts.csv.gz")
    if not os.path.exists(path):
        return {
            "gse": "GSE168204",
            "status": "fail",
            "reason": "counts file missing",
        }
    counts = pd.read_csv(path)
    gene_col = counts.columns[0]
    hits = counts[counts[gene_col].astype(str).str.upper() == "TACSTD2"]
    if hits.empty:
        return {"gse": "GSE168204", "status": "fail", "reason": "TACSTD2 absent"}
    gene = hits.iloc[0]
    colmap = {re.sub(r"\.bam$", "", c): c for c in counts.columns}
    by = {}
    for gsm in samples:
        title = titles.get(gsm, "")
        col = colmap.get(title) or colmap.get(title + ".bam")
        if col is None:
            continue
        pt = chars[gsm].get("patient id") or chars[gsm].get("patient") or title.split("_")[0]
        state = (chars[gsm].get("treatment state") or "").lower()
        tp = "pre" if "pre" in state else ("post" if re.search(r"\bon\b", state) else None)
        if not tp:
            continue
        by.setdefault(pt, {"pre": [], "post": []})[tp].append(float(gene[col]))
    rows = []
    for pt, tps in by.items():
        if tps["pre"] and tps["post"]:
            a, b = float(pd.Series(tps["pre"]).mean()), float(pd.Series(tps["post"]).mean())
            rows.append({"subject": pt, "pre": a, "post": b, "delta": b - a})
    pairs = pd.DataFrame(rows)
    pairs.to_csv(os.path.join(OUT, "GSE168204_pairs.csv"), index=False)
    rec = wilcox_and_sign(pairs)
    rec.update(
        {
            "gse": "GSE168204",
            "unit": "raw_counts",
            "title": "Du 2021 MGH melanoma anti-PD-1/PD-L1 pre/on (companion to GSE115821)",
            "taxon": "Homo sapiens",
            "cancer": "melanoma",
            "treatment": "anti-PD-1 or anti-PD-L1",
            "pairing": "same patient, PRE vs ON tumor RNA-seq",
            "note": "Very small n; not a standalone confirmatory cohort.",
        }
    )
    return rec


def score_gse179351():
    df = pd.read_csv(os.path.join(CACHE, "GSE179351_DESeq2_NormalizedCountsForAllSamples.txt.gz"), sep="\t")
    row = df[df["gene.symbol"].astype(str).str.upper() == "TACSTD2"].iloc[0]
    by = {}
    for c in df.columns:
        m = re.search(r"(\d+).*(Pre|Post).*(Tx|xRT)", c)
        if not m:
            continue
        pt, tp, kind = m.group(1), m.group(2), m.group(3)
        by.setdefault(pt, {})[f"{tp}_{kind}"] = float(row[c])
    rows = []
    for pt, d in by.items():
        if "Pre_Tx" in d and "Post_xRT" in d:
            a, b = d["Pre_Tx"], d["Post_xRT"]
            rows.append({"subject": pt, "pre": a, "post": b, "delta": b - a})
    pairs = pd.DataFrame(rows)
    pairs.to_csv(os.path.join(OUT, "GSE179351_pairs.csv"), index=False)
    rec = wilcox_and_sign(pairs)
    rec.update(
        {
            "gse": "GSE179351",
            "unit": "DESeq2_normalized_counts",
            "title": "Phase II nivo+ipi plus radiation in MSS CRC / PDAC (Parikh et al.)",
            "taxon": "Homo sapiens",
            "cancer": "MSS colorectal and pancreatic adenocarcinoma",
            "treatment": "nivolumab + ipilimumab + radiation (not ICB alone)",
            "pairing": "same patient, Pre-Tx vs Post-xRT metastatic tumor",
            "note": "Confounded by radiation. Not a clean ICB-only confirmatory set.",
        }
    )
    return rec


def score_gse207422_status():
    return {
        "gse": "GSE207422",
        "status": "fail",
        "title": "Hu 2023 NSCLC neoadjuvant PD-1 + chemo (bulk + scRNA)",
        "reason": (
            "Bulk log2TPM (24 samples) is pre-treatment biopsy only. "
            "scRNA metadata has 3 pre and 12 post from different patients (P01–P15 each once). "
            "No within-patient paired TACSTD2 table is deposited."
        ),
        "cancer": "NSCLC",
        "n_paired_subjects": 0,
    }


def main():
    os.makedirs(OUT, exist_ok=True)
    fpkm, both = score_gse91061()
    recs = [
        {**fpkm, "status": "ok"},
        {**score_gse115821(), "status": "ok"},
        {**score_gse168204(), "status": "ok"},
        {**score_gse179351(), "status": "ok"},
        score_gse207422_status(),
    ]
    # drop failed 168204 if n<3
    cleaned = []
    for r in recs:
        if r.get("gse") == "GSE168204" and (r.get("n_paired_subjects") or 0) < 3:
            r["status"] = "fail"
            r["reason"] = f"only {r.get('n_paired_subjects')} paired subjects"
        cleaned.append(r)
    with open(os.path.join(OUT, "known_paired_scores.json"), "w") as fh:
        json.dump({"primary": cleaned, "gse91061_rld_sensitivity": both["rld"]}, fh, indent=2)
    print(json.dumps(cleaned, indent=2))


if __name__ == "__main__":
    main()
