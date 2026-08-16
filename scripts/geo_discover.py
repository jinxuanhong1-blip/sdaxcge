#!/usr/bin/env python3
"""Systematic hunt for GEO series that contain longitudinal (pre vs on/post) ICB samples.

Stage 1 of the GEO arm: enumerate candidate series with E-utilities, pull the lightweight
per-sample metadata for each, and detect subjects that have both a pre-treatment and an
on/post-treatment sample. Nothing here decides whether TACSTD2 is measurable; that happens
in geo_extract.py once a candidate has a usable expression matrix.
"""
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(ROOT, "data", "geo", "cache")
OUT = os.path.join(ROOT, "results", "hunt_paired_up", "geo")

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ACC_CGI = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"

ICB_TERMS = [
    "pembrolizumab", "nivolumab", "atezolizumab", "durvalumab", "avelumab", "cemiplimab",
    "ipilimumab", "tremelimumab", "toripalimab", "sintilimab", "camrelizumab", "tislelizumab",
    "anti-PD-1", "anti-PD1", "anti-PD-L1", "anti-PDL1", "anti-CTLA-4", "anti-CTLA4",
    "immune checkpoint blockade", "immune checkpoint inhibitor", "checkpoint blockade",
    "PD-1 blockade", "PD-L1 blockade", "immunotherapy",
]
LONGITUDINAL_TERMS = [
    "on-treatment", "on treatment", "post-treatment", "post treatment", "pre-treatment",
    "pre treatment", "paired", "longitudinal", "serial", "matched", "before and after",
    "neoadjuvant", "baseline",
]

PRE_PAT = re.compile(
    r"(?:^|[^a-z])(pre[-_ ]?(?:treat|therapy|tx|ici|icb|io|dose|nivo|pembro|anti)?\w*"
    r"|baseline|screening|untreated|naive|na[iï]ve|day\s*0|d0|c1d1|cycle\s*1\s*day\s*1|week\s*0|w0"
    r"|timepoint\s*0|t0|before)(?:[^a-z]|$)",
    re.I,
)
POST_PAT = re.compile(
    r"(?:^|[^a-z])(on[-_ ]?(?:treat|therapy|tx)\w*|post[-_ ]?(?:treat|therapy|tx|ici|icb|io|surg)?\w*"
    r"|during|edt|early\s*during|after|progress\w*|relapse\w*|resist\w*|eot|end\s*of\s*treat\w*"
    r"|resect\w*|surgery|cycle\s*[2-9]|c[2-9]d\d|week\s*[1-9]\d*|w[1-9]\d*|day\s*[1-9]\d*"
    r"|t[1-9]|timepoint\s*[1-9])(?:[^a-z]|$)",
    re.I,
)
SUBJECT_PAT = re.compile(
    r"(patient[\s_-]*[a-z]?\d+|pt[\s_-]*\d+|p\d{1,4}\b|subject[\s_-]*\d+|case[\s_-]*\d+"
    r"|donor[\s_-]*\d+|mrn\d+|[a-z]{2,4}[-_]\d{2,4})",
    re.I,
)


def fetch(url, dest=None, tries=4):
    if dest and os.path.exists(dest) and os.path.getsize(dest) > 0:
        with open(dest, "rb") as fh:
            return fh.read().decode("utf-8", "replace")
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "tacstd2-hunt/1.0"})
            with urllib.request.urlopen(req, timeout=180) as resp:
                body = resp.read()
            if dest:
                os.makedirs(os.path.dirname(dest), exist_ok=True)
                with open(dest, "wb") as fh:
                    fh.write(body)
            return body.decode("utf-8", "replace")
        except Exception as exc:  # noqa: BLE001 - network flakiness is expected; back off
            last = exc
            time.sleep(2 ** i)
    raise RuntimeError(f"failed {url}: {last}")


def esearch_series(term, retmax=400):
    q = urllib.parse.urlencode(
        {"db": "gds", "term": term, "retmax": retmax, "retmode": "json"}
    )
    data = json.loads(fetch(f"{EUTILS}/esearch.fcgi?{q}"))
    return data["esearchresult"].get("idlist", [])


def esummary(uids):
    out = {}
    for i in range(0, len(uids), 200):
        chunk = uids[i : i + 200]
        q = urllib.parse.urlencode(
            {"db": "gds", "id": ",".join(chunk), "retmode": "json"}
        )
        data = json.loads(fetch(f"{EUTILS}/esummary.fcgi?{q}"))
        for uid, rec in data.get("result", {}).items():
            if uid == "uids":
                continue
            out[uid] = rec
        time.sleep(0.34)
    return out


def gsm_metadata(gse):
    """Brief per-sample metadata for a series (much lighter than the full SOFT family file)."""
    url = f"{ACC_CGI}?acc={gse}&targ=gsm&form=text&view=brief"
    return fetch(url, os.path.join(CACHE, f"{gse}_gsm_brief.txt"))


def parse_samples(text):
    samples, cur = [], None
    for line in text.splitlines():
        line = line.strip()
        if line.startswith("^SAMPLE"):
            if cur:
                samples.append(cur)
            cur = {"gsm": line.split("=", 1)[1].strip(), "title": "", "fields": []}
        elif cur is not None and "=" in line:
            key, val = line.split("=", 1)
            key, val = key.strip("! ").strip(), val.strip()
            if key == "Sample_title":
                cur["title"] = val
            elif key in ("Sample_characteristics_ch1", "Sample_source_name_ch1", "Sample_description"):
                cur["fields"].append(val)
    if cur:
        samples.append(cur)
    return samples


def subject_of(sample):
    for f in sample["fields"]:
        low = f.lower()
        if re.match(r"^\s*(patient|patient[\s_]*id|subject|case|donor|individual|pt)\b", low):
            val = f.split(":", 1)[1].strip() if ":" in f else f
            if val:
                return val.lower()
    blob = sample["title"]
    m = SUBJECT_PAT.search(blob)
    return m.group(1).lower().replace(" ", "").replace("-", "_") if m else None


def timepoint_of(sample):
    blob = " | ".join([sample["title"]] + sample["fields"])
    pre, post = bool(PRE_PAT.search(blob)), bool(POST_PAT.search(blob))
    if pre and not post:
        return "pre"
    if post and not pre:
        return "post"
    if pre and post:
        return "ambiguous"
    return None


def profile_series(gse):
    text = gsm_metadata(gse)
    samples = parse_samples(text)
    subjects = {}
    for s in samples:
        subj, tp = subject_of(s), timepoint_of(s)
        if not subj or tp not in ("pre", "post"):
            continue
        subjects.setdefault(subj, set()).add(tp)
    paired = [s for s, tps in subjects.items() if {"pre", "post"} <= tps]
    return {
        "gse": gse,
        "n_samples": len(samples),
        "n_subjects_detected": len(subjects),
        "n_subjects_paired_pre_post": len(paired),
        "paired_subjects_example": sorted(paired)[:8],
    }


def main():
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(CACHE, exist_ok=True)

    icb_q = " OR ".join(f'"{t}"' for t in ICB_TERMS)
    long_q = " OR ".join(f'"{t}"' for t in LONGITUDINAL_TERMS)
    base = f'(({icb_q}) AND ({long_q})) AND "gse"[Filter]'
    queries = {
        "human_rnaseq": f'{base} AND "Homo sapiens"[Organism] AND ("expression profiling by high throughput sequencing"[DataSet Type] OR "expression profiling by array"[DataSet Type])',
        "mouse_rnaseq": f'{base} AND "Mus musculus"[Organism] AND "expression profiling by high throughput sequencing"[DataSet Type]',
        "lung_focus": f'{base} AND ("lung"[All Fields] OR "NSCLC"[All Fields] OR "non-small cell"[All Fields])',
    }

    uids = {}
    for name, term in queries.items():
        ids = esearch_series(term)
        print(f"[esearch] {name}: {len(ids)} series uids")
        uids[name] = ids

    all_uids = sorted({u for ids in uids.values() for u in ids})
    print(f"[esearch] union: {len(all_uids)} uids")
    summ = esummary(all_uids)

    series = {}
    for uid, rec in summ.items():
        acc = rec.get("accession", "")
        if not acc.startswith("GSE"):
            continue
        series[acc] = {
            "gse": acc,
            "title": rec.get("title", ""),
            "summary": (rec.get("summary", "") or "")[:1200],
            "taxon": rec.get("taxon", ""),
            "n_samples": int(rec.get("n_samples", 0) or 0),
            "gdstype": rec.get("gdstype", ""),
            "pdat": rec.get("PDAT", ""),
            "found_in_queries": [q for q, ids in uids.items() if uid in ids],
        }
    print(f"[esummary] {len(series)} distinct GSE")
    with open(os.path.join(OUT, "geo_candidate_series.json"), "w") as fh:
        json.dump(series, fh, indent=2)

    profiles = {}
    skipped = {}
    todo = []
    for acc, rec in sorted(series.items()):
        if rec["n_samples"] > 1500 or rec["n_samples"] < 4:
            skipped[acc] = f"n_samples={rec['n_samples']} outside 4..1500"
        else:
            todo.append((acc, rec))

    def work(item):
        acc, rec = item
        try:
            prof = profile_series(acc)
            prof.update({k: rec[k] for k in ("title", "taxon", "gdstype", "found_in_queries")})
            return acc, prof, None
        except Exception as exc:  # noqa: BLE001 - keep going, record the failure
            return acc, None, f"metadata fetch/parse failed: {exc}"

    with ThreadPoolExecutor(max_workers=6) as pool:
        for i, (acc, prof, err) in enumerate(pool.map(work, todo), 1):
            if prof is not None:
                profiles[acc] = prof
            else:
                skipped[acc] = err
            if i % 25 == 0:
                print(f"  profiled {i}/{len(todo)}", flush=True)

    with open(os.path.join(OUT, "geo_series_pairing_profile.json"), "w") as fh:
        json.dump({"profiles": profiles, "skipped": skipped}, fh, indent=2)

    ranked = sorted(
        profiles.values(), key=lambda r: r["n_subjects_paired_pre_post"], reverse=True
    )
    print(f"\nseries profiled: {len(profiles)}, skipped: {len(skipped)}")
    print("top candidates by detected pre/post paired subjects:")
    for r in ranked[:40]:
        print(
            f"  {r['gse']:>12}  paired={r['n_subjects_paired_pre_post']:>3}  "
            f"n={r['n_samples']:>4}  {r['taxon'][:20]:20}  {r['title'][:80]}"
        )


if __name__ == "__main__":
    sys.exit(main())
