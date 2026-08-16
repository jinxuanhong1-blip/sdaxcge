#!/usr/bin/env python3
"""Download GEO series matrices / supplementary tables and score TACSTD2 in paired ICB.

A series is scored only when:
  1. TACSTD2 / Tacstd2 is present in a usable expression table, and
  2. at least 3 subjects have both a pre and an on/post sample.

Anything else is recorded as a failure with a reason, not silently dropped.
"""
import gzip
import io
import json
import os
import re
import sys
import time
import urllib.error
import urllib.request
from collections import defaultdict

import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(ROOT, "data", "geo", "cache")
OUT = os.path.join(ROOT, "results", "hunt_paired_up", "geo")

# Curated seed: published ICB cohorts that are widely cited as having (or claiming)
# longitudinal tumour expression. Blood-only / scRNA-only / methylation-only series
# are included so the hunt can fail them honestly rather than skip them.
SEED = [
    # human tumour, classic paired ICB
    "GSE91061",  # Riaz 2017 melanoma nivo pre/on
    "GSE78220",  # Hugo 2016 melanoma anti-PD-1 (mostly pre)
    "GSE93157",  # Prat 2017 nanostring melanoma/NSCLC
    "GSE115821",
    "GSE123814",  # Yost 2019 BCC/SCC scRNA
    "GSE126044",  # Duruisseaux NSCLC methylation
    "GSE135222",  # Jung NSCLC ICB (mostly baseline)
    "GSE136961",
    "GSE145996",
    "GSE148673",
    "GSE162454",
    "GSE165252",
    "GSE168204",
    "GSE173839",  # I-SPY2
    "GSE176307",  # Rose bladder
    "GSE179351",
    "GSE181815",
    "GSE183924",
    "GSE189466",
    "GSE194040",  # I-SPY2
    "GSE199545",
    "GSE201425",
    "GSE205152",
    "GSE207422",  # Hu 2023 NSCLC neoadjuvant chemo-IO, HAS bulk TPM
    "GSE208652",
    "GSE215120",
    "GSE215868",
    "GSE217206",
    "GSE218989",
    "GSE221601",
    "GSE222219",
    "GSE223454",
    "GSE227666",  # NeoPembrOV
    "GSE231535",
    "GSE236581",
    "GSE240423",
    "GSE241176",
    "GSE247765",
    "GSE248378",  # neoadjuvant durvalumab NSCLC
    "GSE260770",  # sintilimab GGO
    "GSE262376",
    "GSE264434",
    "GSE268628",
    "GSE271757",
    "GSE272436",
    "GSE273160",
    "GSE279881",
    "GSE288199",
    "GSE298296",
    "GSE318645",
    "GSE319641",  # breast pre/on
    "GSE330465",
    # mouse / mixed
    "GSE168846",  # syngeneic anti-PD-1
    "GSE210547",  # NCI-H460 anti-PD1 in vivo
    "GSE155972",  # LLC in TISMO
]

PRE_RE = re.compile(
    r"(pre[-_ ]?(treat|therapy|tx|ici|icb|io|dose|nivo|pembro)?|baseline|screening|"
    r"untreated|naive|na[iï]ve|day\s*0|d0|c1d1|week\s*0|w0|before|preop|pre-op|"
    r"timepoint\s*0|\bt0\b)",
    re.I,
)
POST_RE = re.compile(
    r"(on[-_ ]?(treat|therapy|tx)|post[-_ ]?(treat|therapy|tx|ici|icb|io|surg)|"
    r"during|after|progress|relapse|resist|eot|end of treat|resect|surgery|"
    r"cycle\s*[2-9]|c[2-9]d|week\s*[1-9]|w[1-9]|day\s*[1-9]|postop|post-op|"
    r"on-treatment|ontreatment|post-treatment|posttreatment)",
    re.I,
)
SUBJ_KEYS = re.compile(
    r"^(patient|patient[\s_]*id|subject|subject[\s_]*id|case|donor|individual|"
    r"pt|pt[\s_]*id|sample[\s_]*id|donor[\s_]*id|participant)$",
    re.I,
)
GENE_ALIASES = ("TACSTD2", "Tacstd2", "tacstd2", "TROP2", "Trop2", "EGP1", "GA733-1", "M1S1")


def ftp_series_dir(gse):
    n = int(gse[3:])
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE{n // 1000}nnn/{gse}"


def fetch(url, dest, tries=4):
    if dest and os.path.exists(dest) and os.path.getsize(dest) > 200:
        return dest
    last = None
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "tacstd2-hunt/1.0"})
            with urllib.request.urlopen(req, timeout=180) as resp, open(dest, "wb") as fh:
                while True:
                    chunk = resp.read(1 << 20)
                    if not chunk:
                        break
                    fh.write(chunk)
            return dest
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(1.5 * (i + 1))
    raise RuntimeError(f"{url}: {last}")


def try_fetch(url, dest):
    try:
        return fetch(url, dest)
    except Exception as exc:  # noqa: BLE001
        return None, str(exc)


def parse_series_matrix(path):
    opener = gzip.open if path.endswith(".gz") else open
    header, samples, characteristics, titles, sources = {}, [], defaultdict(dict), {}, {}
    table_lines = []
    in_table = False
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if in_table:
                table_lines.append(line)
                continue
            if not line.startswith("!"):
                continue
            key, _, rest = line[1:].partition("\t")
            key = key.strip()
            vals = [v.strip().strip('"') for v in rest.rstrip("\n").split("\t")]
            if key == "Sample_geo_accession":
                samples = vals
            elif key == "Sample_title":
                titles = dict(zip(samples, vals))
            elif key == "Sample_source_name_ch1":
                sources = dict(zip(samples, vals))
            elif key.startswith("Sample_characteristics"):
                for gsm, val in zip(samples, vals):
                    if ":" in val:
                        k, v = val.split(":", 1)
                        characteristics[gsm][k.strip().lower()] = v.strip()
                    else:
                        characteristics[gsm][f"raw_{key}"] = val
            elif key.startswith("Series_"):
                header[key] = " ".join(vals)[:1500]
    expr = None
    if table_lines:
        buf = io.StringIO("".join(table_lines))
        expr = pd.read_csv(buf, sep="\t", index_col=0)
        expr.columns = [c.strip('"') for c in expr.columns]
        expr.index = [str(i).strip('"') for i in expr.index]
    return {
        "header": header,
        "samples": samples,
        "titles": titles,
        "sources": sources,
        "characteristics": dict(characteristics),
        "expr": expr,
    }


def subject_of(gsm, titles, sources, chars):
    for k, v in chars.get(gsm, {}).items():
        if SUBJ_KEYS.match(k) and v and v.lower() not in ("na", "n/a", "none"):
            return v.strip().lower()
    blob = " ".join([titles.get(gsm, ""), sources.get(gsm, "")])
    m = re.search(
        r"(patient[\s_-]*[a-z]?\d+|pt[\s_-]*\d+|p\d{1,4}\b|subject[\s_-]*\d+|"
        r"case[\s_-]*\d+|[a-z]{2,6}[-_]\d{2,4})",
        blob,
        re.I,
    )
    return m.group(1).lower().replace(" ", "") if m else None


def timepoint_of(gsm, titles, sources, chars):
    blob = " | ".join(
        [titles.get(gsm, ""), sources.get(gsm, "")]
        + [f"{k}:{v}" for k, v in chars.get(gsm, {}).items()]
    )
    pre, post = bool(PRE_RE.search(blob)), bool(POST_RE.search(blob))
    if pre and not post:
        return "pre"
    if post and not pre:
        return "post"
    if pre and post:
        # titles like "pre-treatment vs on-treatment" are series-level noise; prefer
        # the more specific token if only one side is in the sample title.
        t = titles.get(gsm, "")
        if PRE_RE.search(t) and not POST_RE.search(t):
            return "pre"
        if POST_RE.search(t) and not PRE_RE.search(t):
            return "post"
        return "ambiguous"
    return None


def pick_gene_row(expr):
    if expr is None or expr.empty:
        return None, None
    idx = pd.Index(expr.index.astype(str))
    for alias in GENE_ALIASES:
        hits = idx[idx.str.upper() == alias.upper()]
        if len(hits):
            return hits[0], expr.loc[hits[0]]
        hits = idx[idx.str.upper().str.startswith(alias.upper() + "|") | idx.str.upper().str.endswith("|" + alias.upper())]
        if len(hits):
            return hits[0], expr.loc[hits[0]]
    # Ensembl-style rows sometimes carry the symbol in a second column already used as index
    return None, None


def pair_table(meta_rows, gene_series):
    by = defaultdict(lambda: {"pre": [], "post": []})
    for r in meta_rows:
        if r["subject"] and r["timepoint"] in ("pre", "post") and r["gsm"] in gene_series.index:
            by[r["subject"]][r["timepoint"]].append(r["gsm"])
    pairs = []
    for subj, tps in by.items():
        if not tps["pre"] or not tps["post"]:
            continue
        pre = float(pd.to_numeric(gene_series[tps["pre"]], errors="coerce").mean())
        post = float(pd.to_numeric(gene_series[tps["post"]], errors="coerce").mean())
        pairs.append(
            {
                "subject": subj,
                "n_pre": len(tps["pre"]),
                "n_post": len(tps["post"]),
                "pre": pre,
                "post": post,
                "delta": post - pre,
                "pre_gsms": ";".join(tps["pre"]),
                "post_gsms": ";".join(tps["post"]),
            }
        )
    return pd.DataFrame(pairs)


def summarise_pairs(pairs, gse, extra):
    n = len(pairs)
    if n == 0:
        return None
    up = int((pairs.delta > 0).sum())
    rec = {
        "gse": gse,
        "n_paired_subjects": n,
        "n_up": up,
        "n_down": int((pairs.delta < 0).sum()),
        "n_flat": int((pairs.delta == 0).sum()),
        "frac_up": round(up / n, 3),
        "median_delta": float(pairs.delta.median()),
        "mean_pre": float(pairs.pre.mean()),
        "mean_post": float(pairs.post.mean()),
        "wilcoxon_p": (
            float(stats.wilcoxon(pairs.post, pairs.pre, alternative="two-sided", zero_method="wilcox").pvalue)
            if n >= 6 and (pairs.delta != 0).sum() >= 6
            else float("nan")
        ),
        "sign_test_p": float(stats.binomtest(up, n, 0.5).pvalue) if n else float("nan"),
    }
    rec.update(extra)
    return rec


def process_gse207422():
    """Special-case: bulk log2TPM lives in a supplementary file, not the series matrix."""
    dest_x = os.path.join(CACHE, "GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz")
    dest_m = os.path.join(CACHE, "GSE207422_NSCLC_bulk_RNAseq_metadata.xlsx")
    base = ftp_series_dir("GSE207422") + "/suppl/"
    fetch(base + "GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz", dest_x)
    fetch(base + "GSE207422_NSCLC_bulk_RNAseq_metadata.xlsx", dest_m)
    expr = pd.read_csv(dest_x, sep="\t", index_col=0)
    meta = pd.read_excel(dest_m)
    gene = None
    for alias in GENE_ALIASES:
        if alias in expr.index:
            gene = alias
            break
    if gene is None:
        return {"gse": "GSE207422", "status": "fail", "reason": "TACSTD2 absent from bulk TPM"}, None, meta
    # metadata columns vary; keep whatever is there
    meta = meta.copy()
    meta.columns = [str(c).strip() for c in meta.columns]
    sample_col = next((c for c in meta.columns if meta[c].astype(str).isin(expr.columns).sum() > 3), None)
    if sample_col is None:
        # try using the first column as sample id
        sample_col = meta.columns[0]
    meta["sample"] = meta[sample_col].astype(str)
    meta = meta[meta["sample"].isin(expr.columns)].copy()
    meta["TACSTD2"] = expr.loc[gene, meta["sample"]].to_numpy()
    # infer subject + timepoint from any column
    blob_cols = [c for c in meta.columns if c != "TACSTD2"]
    def row_blob(row):
        return " | ".join(str(row[c]) for c in blob_cols)
    meta["subject"] = None
    meta["timepoint"] = None
    for i, row in meta.iterrows():
        b = row_blob(row)
        m = re.search(r"(P\d+|patient[\s_-]*\d+|pt[\s_-]*\d+|[A-Z]{1,4}\d{1,3})", b, re.I)
        meta.at[i, "subject"] = m.group(1).lower() if m else None
        pre, post = bool(PRE_RE.search(b)), bool(POST_RE.search(b))
        meta.at[i, "timepoint"] = "pre" if pre and not post else ("post" if post and not pre else None)
    # if timepoint still missing, look at obvious column names
    for c in meta.columns:
        cl = c.lower()
        if "time" in cl or "treat" in cl or "group" in cl or "response" in cl:
            pass
    pairs_src = []
    for _, row in meta.iterrows():
        if row.subject and row.timepoint in ("pre", "post"):
            pairs_src.append({"gsm": row["sample"], "subject": row.subject, "timepoint": row.timepoint})
    gene_s = pd.Series(meta.set_index("sample")["TACSTD2"])
    pairs = pair_table(pairs_src, gene_s)
    extra = {
        "status": "ok" if len(pairs) >= 3 else "fail",
        "reason": None if len(pairs) >= 3 else f"only {len(pairs)} paired subjects in bulk TPM",
        "title": "NSCLC neoadjuvant PD-1 + chemo (Hu et al.; bulk TPM in GSE207422)",
        "taxon": "Homo sapiens",
        "tissue": "lung_tumor_bulk",
        "source_file": "GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz",
        "gene_row": gene,
        "n_samples_in_matrix": int(expr.shape[1]),
        "metadata_columns": list(meta.columns),
    }
    rec = summarise_pairs(pairs, "GSE207422", extra) or extra
    rec.update(extra)
    rec["n_paired_subjects"] = int(len(pairs))
    return rec, pairs, meta


def process_series_matrix(gse, title="", taxon=""):
    dest = os.path.join(CACHE, f"{gse}_series_matrix.txt.gz")
    url = f"{ftp_series_dir(gse)}/matrix/{gse}_series_matrix.txt.gz"
    try:
        fetch(url, dest)
    except Exception as exc:  # noqa: BLE001
        return {"gse": gse, "status": "fail", "reason": f"series matrix missing: {exc}", "title": title, "taxon": taxon}, None, None
    parsed = parse_series_matrix(dest)
    expr = parsed["expr"]
    gene_name, gene_s = pick_gene_row(expr)
    meta_rows = []
    for gsm in parsed["samples"]:
        meta_rows.append(
            {
                "gsm": gsm,
                "title": parsed["titles"].get(gsm, ""),
                "source": parsed["sources"].get(gsm, ""),
                "subject": subject_of(gsm, parsed["titles"], parsed["sources"], parsed["characteristics"]),
                "timepoint": timepoint_of(gsm, parsed["titles"], parsed["sources"], parsed["characteristics"]),
                "chars": parsed["characteristics"].get(gsm, {}),
            }
        )
    extra = {
        "title": title or parsed["header"].get("Series_title", ""),
        "taxon": taxon or parsed["header"].get("Series_sample_organism", ""),
        "source_file": os.path.basename(dest),
        "gene_row": gene_name,
        "n_samples_in_matrix": int(expr.shape[1]) if expr is not None else 0,
        "n_subjects_detected": len({r["subject"] for r in meta_rows if r["subject"]}),
        "n_pre_labelled": sum(r["timepoint"] == "pre" for r in meta_rows),
        "n_post_labelled": sum(r["timepoint"] == "post" for r in meta_rows),
    }
    if gene_s is None:
        extra.update({"gse": gse, "status": "fail", "reason": "TACSTD2/Tacstd2 not in series matrix"})
        return extra, None, pd.DataFrame(meta_rows)
    gene_s.index = [str(i) for i in gene_s.index]
    pairs = pair_table(meta_rows, gene_s)
    extra["n_paired_subjects"] = int(len(pairs))
    if len(pairs) < 3:
        extra.update({"gse": gse, "status": "fail", "reason": f"TACSTD2 present but only {len(pairs)} paired subjects"})
        return extra, pairs, pd.DataFrame(meta_rows)
    rec = summarise_pairs(pairs, gse, extra)
    rec["status"] = "ok"
    rec["reason"] = None
    rec["gse"] = gse
    return rec, pairs, pd.DataFrame(meta_rows)


def main():
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(CACHE, exist_ok=True)
    extra = [a for a in sys.argv[1:] if a.startswith("GSE")]
    wanted = list(dict.fromkeys(SEED + extra))
    rows, failures = [], []
    for gse in wanted:
        print(f"== {gse}", flush=True)
        try:
            if gse == "GSE207422":
                rec, pairs, meta = process_gse207422()
            else:
                rec, pairs, meta = process_series_matrix(gse)
        except Exception as exc:  # noqa: BLE001
            rec, pairs, meta = {"gse": gse, "status": "fail", "reason": f"exception: {exc}"}, None, None
        if rec.get("status") == "ok":
            rows.append(rec)
            if pairs is not None:
                pairs.assign(gse=gse).to_csv(os.path.join(OUT, f"{gse}_pairs.csv"), index=False)
            print(f"   OK paired={rec['n_paired_subjects']} up={rec['n_up']} p={rec.get('wilcoxon_p')}")
        else:
            failures.append(rec)
            print(f"   FAIL {rec.get('reason')}")
        if meta is not None and len(meta):
            meta.to_csv(os.path.join(OUT, f"{gse}_sample_meta.csv"), index=False)
        time.sleep(0.15)

    ok = pd.DataFrame(rows)
    fail = pd.DataFrame(failures)
    if len(ok):
        ok.sort_values("n_paired_subjects", ascending=False).to_csv(
            os.path.join(OUT, "geo_paired_tacstd2_ok.csv"), index=False
        )
    fail.to_csv(os.path.join(OUT, "geo_paired_tacstd2_failures.csv"), index=False)
    summary = {
        "n_series_attempted": len(wanted),
        "n_ok_paired_with_tacstd2": int(len(ok)),
        "n_failed": int(len(fail)),
        "ok": rows,
        "failures_by_reason": fail.reason.value_counts().to_dict() if len(fail) else {},
    }
    with open(os.path.join(OUT, "geo_extract_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps({k: summary[k] for k in ("n_series_attempted", "n_ok_paired_with_tacstd2", "n_failed")}, indent=2))


if __name__ == "__main__":
    main()
