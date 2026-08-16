#!/usr/bin/env python3
"""Score TACSTD2/Tacstd2 in public paired ICB GEO series.

A series is confirmatory only when the gene is quantified AND >=3 subjects have both
a pre and an on/post sample. Empty series matrices (common for RNA-seq) fall back to
supplementary expression tables. Gene IDs accepted: symbol, Entrez, Ensembl.
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
from html.parser import HTMLParser

import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(ROOT, "data", "geo", "cache")
OUT = os.path.join(ROOT, "results", "hunt_paired_up", "geo")

HUMAN_IDS = {"TACSTD2", "TROP2", "EGP1", "GA733-1", "M1S1", "4070", "ENSG00000184292"}
MOUSE_IDS = {"TACSTD2", "TACSTD2", "TROP2", "56753", "ENSMUSG00000050371"}

SEED = [
    "GSE91061", "GSE78220", "GSE93157", "GSE115821", "GSE123814", "GSE126044",
    "GSE135222", "GSE136961", "GSE145996", "GSE148673", "GSE162454", "GSE165252",
    "GSE168204", "GSE173839", "GSE176307", "GSE179351", "GSE181815", "GSE183924",
    "GSE189466", "GSE194040", "GSE199545", "GSE201425", "GSE205152", "GSE207422",
    "GSE208652", "GSE215120", "GSE215868", "GSE217206", "GSE218989", "GSE221601",
    "GSE222219", "GSE223454", "GSE227666", "GSE231535", "GSE236581", "GSE240423",
    "GSE241176", "GSE247765", "GSE248378", "GSE260770", "GSE262376", "GSE264434",
    "GSE268628", "GSE271757", "GSE272436", "GSE273160", "GSE279881", "GSE288199",
    "GSE298296", "GSE318645", "GSE319641", "GSE330465", "GSE168846", "GSE210547",
    "GSE155972",
]

PRE_RE = re.compile(
    r"(pre[-_ ]?(treat|therapy|tx|ici|icb|io|dose|nivo|pembro)?|baseline|screening|"
    r"untreated|naive|na[iï]ve|day\s*0|\bd0\b|c1d1|week\s*0|\bw0\b|before|preop|"
    r"pre-op|timepoint\s*0|\bt0\b|_pre_|\bpre\b)",
    re.I,
)
POST_RE = re.compile(
    r"(on[-_ ]?(treat|therapy|tx)|post[-_ ]?(treat|therapy|tx|ici|icb|io|surg)|"
    r"during|after|progress|relapse|resist|eot|end of treat|resect|surgery|"
    r"cycle\s*[2-9]|c[2-9]d|week\s*[1-9]|\bw[1-9]|day\s*[1-9]|postop|post-op|"
    r"on-treatment|ontreatment|post-treatment|posttreatment|_on_|\bon\b)",
    re.I,
)
SUBJ_KEYS = re.compile(
    r"^(patient|patient[\s_]*id|subject|subject[\s_]*id|case|donor|individual|"
    r"pt|pt[\s_]*id|participant)$",
    re.I,
)
EXPR_NAME = re.compile(
    r"(fpkm|tpm|counts?|expression|rld|rlog|vst|cpm|rpkm|norm|matrix|log2)",
    re.I,
)
SKIP_NAME = re.compile(
    r"(raw\.fastq|fastq\.gz|bam|bai|h5ad|rds|h5|mtx|barcodes|peaks|bed|vcf|bw|bigwig|hic)",
    re.I,
)


class HrefParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.hrefs = []

    def handle_starttag(self, tag, attrs):
        if tag == "a":
            href = dict(attrs).get("href")
            if href:
                self.hrefs.append(href)


def ftp_dir(gse):
    n = int(gse[3:])
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE{n // 1000}nnn/{gse}"


def fetch_bytes(url, dest, tries=4, min_size=200):
    if dest and os.path.exists(dest) and os.path.getsize(dest) >= min_size:
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
            time.sleep(1.2 * (i + 1))
    raise RuntimeError(f"{url}: {last}")


def list_ftp(url):
    dest = os.path.join(CACHE, "listings", re.sub(r"[^A-Za-z0-9]+", "_", url) + ".html")
    try:
        fetch_bytes(url if url.endswith("/") else url + "/", dest, min_size=20)
    except Exception:
        return []
    parser = HrefParser()
    parser.feed(open(dest, encoding="utf-8", errors="replace").read())
    out = []
    for h in parser.hrefs:
        if h in ("../", "./") or h.startswith("?") or h.startswith("http"):
            continue
        out.append(h.rstrip("/"))
    return out


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
                header[key] = " ".join(vals)[:2000]
    expr = None
    if table_lines:
        buf = io.StringIO("".join(table_lines))
        try:
            expr = pd.read_csv(buf, sep="\t", index_col=0)
            expr.columns = [str(c).strip().strip('"') for c in expr.columns]
            expr.index = [str(i).strip().strip('"') for i in expr.index]
            if expr.shape[0] < 2:
                expr = None
        except Exception:
            expr = None
    return {
        "header": header,
        "samples": samples,
        "titles": titles,
        "sources": sources,
        "characteristics": dict(characteristics),
        "expr": expr,
    }


def read_table(path):
    last = None
    for sep in ("\t", ",", None):
        try:
            df = pd.read_csv(path, sep=sep, index_col=0, compression="infer")
            if df.shape[1] >= 2:
                df.index = df.index.astype(str)
                df.columns = df.columns.astype(str)
                return df
        except Exception as exc:  # noqa: BLE001
            last = exc
    raise RuntimeError(f"cannot parse {path}: {last}")


def pick_gene_row(expr):
    if expr is None or expr.empty:
        return None, None
    idx = pd.Index(expr.index.astype(str))
    upper = idx.str.upper()
    aliases = HUMAN_IDS | {x.upper() for x in MOUSE_IDS} | {"TACSTD2", "TROP2"}
    for alias in aliases:
        hits = idx[upper == alias.upper()]
        if len(hits):
            return str(hits[0]), pd.to_numeric(expr.loc[hits[0]], errors="coerce")
        hits = idx[upper.str.startswith(alias.upper() + "|") | upper.str.endswith("|" + alias.upper())]
        if len(hits):
            return str(hits[0]), pd.to_numeric(expr.loc[hits[0]], errors="coerce")
        hits = idx[upper.str.contains(r"(?:^|[^A-Z0-9])" + re.escape(alias.upper()) + r"(?:[^A-Z0-9]|$)")]
        if len(hits) == 1:
            return str(hits[0]), pd.to_numeric(expr.loc[hits[0]], errors="coerce")
    return None, None


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
        t = titles.get(gsm, "")
        if PRE_RE.search(t) and not POST_RE.search(t):
            return "pre"
        if POST_RE.search(t) and not PRE_RE.search(t):
            return "post"
        return "ambiguous"
    return None


def pair_from_columns(gene_s):
    """Pair when sample names themselves encode PtX_Pre / PtX_On."""
    by = defaultdict(lambda: {"pre": [], "post": []})
    for col in gene_s.index.astype(str):
        m = re.match(r"(pt\d+|patient\d+|p\d+)[_\-](pre|on|post|baseline)", col, re.I)
        if not m:
            continue
        subj, tp = m.group(1).lower(), m.group(2).lower()
        by[subj]["pre" if tp in ("pre", "baseline") else "post"].append(col)
    return _pairs_from_map(by, gene_s)


def pair_from_meta(meta_rows, gene_s):
    by = defaultdict(lambda: {"pre": [], "post": []})
    for r in meta_rows:
        if r["subject"] and r["timepoint"] in ("pre", "post") and r["gsm"] in gene_s.index:
            by[r["subject"]][r["timepoint"]].append(r["gsm"])
    return _pairs_from_map(by, gene_s)


def _pairs_from_map(by, gene_s):
    rows = []
    for subj, tps in by.items():
        if not tps["pre"] or not tps["post"]:
            continue
        pre = float(pd.to_numeric(gene_s[tps["pre"]], errors="coerce").mean())
        post = float(pd.to_numeric(gene_s[tps["post"]], errors="coerce").mean())
        rows.append(
            {
                "subject": subj,
                "n_pre": len(tps["pre"]),
                "n_post": len(tps["post"]),
                "pre": pre,
                "post": post,
                "delta": post - pre,
                "pre_samples": ";".join(tps["pre"]),
                "post_samples": ";".join(tps["post"]),
            }
        )
    return pd.DataFrame(rows)


def summarise(pairs, gse, extra):
    n = len(pairs)
    rec = dict(extra)
    rec.update({"gse": gse, "n_paired_subjects": int(n)})
    if n < 3:
        rec["status"] = "fail"
        rec["reason"] = rec.get("reason") or f"TACSTD2 present but only {n} paired subjects"
        return rec
    up = int((pairs.delta > 0).sum())
    rec.update(
        {
            "status": "ok",
            "reason": None,
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
            "sign_test_p": float(stats.binomtest(up, n, 0.5).pvalue),
        }
    )
    return rec


def load_matrix_and_meta(gse):
    names = list_ftp(ftp_dir(gse) + "/matrix")
    matrices = [n for n in names if "series_matrix" in n]
    if not matrices:
        # default name even if listing failed
        matrices = [f"{gse}_series_matrix.txt.gz"]
    parsed = None
    last_err = None
    for name in matrices:
        dest = os.path.join(CACHE, name)
        try:
            fetch_bytes(f"{ftp_dir(gse)}/matrix/{name}", dest)
            parsed = parse_series_matrix(dest)
            if parsed["samples"]:
                return parsed
        except Exception as exc:  # noqa: BLE001
            last_err = exc
    if parsed is None:
        raise RuntimeError(f"no series matrix: {last_err}")
    return parsed


def candidate_suppl(gse):
    names = list_ftp(ftp_dir(gse) + "/suppl")
    keep = []
    for n in names:
        if SKIP_NAME.search(n):
            continue
        if n.endswith((".gz", ".txt", ".csv", ".tsv", ".xlsx", ".xls")) and (
            EXPR_NAME.search(n) or n.lower().endswith((".txt.gz", ".csv.gz", ".tsv.gz", ".txt", ".csv"))
        ):
            # skip huge scRNA UMI dumps
            if re.search(r"(scRNA|UMI|barcodes|filtered_feature)", n, re.I) and "bulk" not in n.lower():
                continue
            keep.append(n)
    # prefer smaller, gene-level tables
    keep.sort(key=lambda n: (0 if re.search(r"bulk|fpkm|tpm|rld", n, re.I) else 1, len(n)))
    return keep[:8]


def process(gse):
    reasons = []
    parsed = None
    try:
        parsed = load_matrix_and_meta(gse)
    except Exception as exc:  # noqa: BLE001
        reasons.append(f"matrix: {exc}")

    title = parsed["header"].get("Series_title", "") if parsed else ""
    taxon = parsed["header"].get("Series_sample_organism", "") if parsed else ""
    extra = {"title": title, "taxon": taxon}

    meta_rows = []
    if parsed:
        for gsm in parsed["samples"]:
            meta_rows.append(
                {
                    "gsm": gsm,
                    "title": parsed["titles"].get(gsm, ""),
                    "source": parsed["sources"].get(gsm, ""),
                    "subject": subject_of(gsm, parsed["titles"], parsed["sources"], parsed["characteristics"]),
                    "timepoint": timepoint_of(gsm, parsed["titles"], parsed["sources"], parsed["characteristics"]),
                }
            )
        extra.update(
            {
                "n_samples_meta": len(meta_rows),
                "n_subjects_detected": len({r["subject"] for r in meta_rows if r["subject"]}),
                "n_pre_labelled": sum(r["timepoint"] == "pre" for r in meta_rows),
                "n_post_labelled": sum(r["timepoint"] == "post" for r in meta_rows),
            }
        )

    tables = []
    if parsed and parsed["expr"] is not None:
        tables.append(("series_matrix", parsed["expr"]))
    for name in candidate_suppl(gse):
        dest = os.path.join(CACHE, name)
        url = f"{ftp_dir(gse)}/suppl/{name}"
        try:
            fetch_bytes(url, dest)
            if name.endswith((".xlsx", ".xls")):
                df = pd.read_excel(dest, index_col=0)
                df.index = df.index.astype(str)
                df.columns = df.columns.astype(str)
            else:
                df = read_table(dest)
            tables.append((name, df))
        except Exception as exc:  # noqa: BLE001
            reasons.append(f"suppl {name}: {exc}")

    best_pairs, best_rec, best_gene, best_src = None, None, None, None
    for src, expr in tables:
        gene, gene_s = pick_gene_row(expr)
        if gene_s is None:
            continue
        gene_s.index = gene_s.index.astype(str)
        pairs = pair_from_columns(gene_s)
        if len(pairs) < 3 and meta_rows:
            # map GSM columns if the table uses GSM ids
            pairs_m = pair_from_meta(meta_rows, gene_s)
            if len(pairs_m) > len(pairs):
                pairs = pairs_m
        extra_i = dict(extra)
        extra_i.update({"source_file": src, "gene_row": gene, "n_samples_in_matrix": int(expr.shape[1])})
        rec = summarise(pairs, gse, extra_i)
        if rec["status"] == "ok":
            return rec, pairs, pd.DataFrame(meta_rows)
        if best_rec is None or rec.get("n_paired_subjects", 0) > best_rec.get("n_paired_subjects", 0):
            best_rec, best_pairs, best_gene, best_src = rec, pairs, gene, src

    if best_rec is not None:
        best_rec["reason"] = best_rec.get("reason") or "gene found but pairing failed"
        best_rec["gene_row"] = best_gene
        best_rec["source_file"] = best_src
        return best_rec, best_pairs, pd.DataFrame(meta_rows)

    extra.update(
        {
            "gse": gse,
            "status": "fail",
            "reason": "; ".join(reasons) if reasons else "TACSTD2/Tacstd2 not found in matrix or supplementary tables",
            "n_tables_tried": len(tables),
        }
    )
    return extra, None, pd.DataFrame(meta_rows)


def main():
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(CACHE, exist_ok=True)
    extra = [a for a in sys.argv[1:] if a.startswith("GSE")]
    wanted = list(dict.fromkeys((extra or SEED)))
    rows, failures = [], []
    for gse in wanted:
        print(f"== {gse}", flush=True)
        try:
            rec, pairs, meta = process(gse)
        except Exception as exc:  # noqa: BLE001
            rec, pairs, meta = {"gse": gse, "status": "fail", "reason": f"exception: {exc}"}, None, None
        if rec.get("status") == "ok":
            rows.append(rec)
            pairs.assign(gse=gse).to_csv(os.path.join(OUT, f"{gse}_pairs.csv"), index=False)
            print(
                f"   OK n={rec['n_paired_subjects']} up={rec['n_up']} "
                f"frac={rec['frac_up']} wilcox={rec.get('wilcoxon_p')} src={rec.get('source_file')}",
                flush=True,
            )
        else:
            failures.append(rec)
            print(f"   FAIL {rec.get('reason')}", flush=True)
        if meta is not None and len(meta):
            meta.to_csv(os.path.join(OUT, f"{gse}_sample_meta.csv"), index=False)
        time.sleep(0.1)

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
        "failures": failures,
    }
    with open(os.path.join(OUT, "geo_extract_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2, default=str)
    print(json.dumps({k: summary[k] for k in ("n_series_attempted", "n_ok_paired_with_tacstd2", "n_failed")}, indent=2))


if __name__ == "__main__":
    main()
