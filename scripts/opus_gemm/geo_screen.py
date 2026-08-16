#!/usr/bin/env python3
"""Systematic screen of GEO for mouse lung-cancer studies with a real immune
checkpoint inhibitor (ICI) arm, in models other than plain LLC.

Steps
  1. esearch a set of broad queries against GEO DataSets (db=gds), mouse only.
  2. For every hit series, download the (small) MINiML family XML and read the
     overall design + per-sample characteristics.
  3. Score each series for: ICI evidence, model identity (KP / KL / EGFR /
     SCLC-GEMM / CMT167 / 344SQ / ...), in-vivo tissue, and LLC-only status.
  4. Emit a triage table plus a JSON cache of the parsed metadata.

The MINiML payloads are cached under --cache so re-runs are cheap and the
screen is reproducible.

Usage:
    python3 scripts/opus_gemm/geo_screen.py \
        --out results/opus_gemm/geo_screen.tsv \
        --json results/opus_gemm/geo_screen.json
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import os
import re
import sys
import tarfile
import time
import urllib.error
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"

MOUSE = '"Mus musculus"[Organism]'
SEQ = '("expression profiling by high throughput sequencing"[DataSet Type] OR "expression profiling by array"[DataSet Type])'

QUERIES = [
    f'lung AND (anti-PD-1 OR "PD-1 blockade" OR anti-PD-L1 OR anti-CTLA-4 OR "checkpoint blockade" OR "checkpoint inhibitor" OR immunotherapy) AND {MOUSE} AND {SEQ}',
    f'(KrasLSL-G12D OR "Kras G12D" OR KrasG12D) AND (Trp53 OR p53) AND lung AND {MOUSE} AND {SEQ}',
    f'(SCLC OR "small cell lung") AND (Rb1 OR Trp53 OR MycT58A OR RPM OR RPP) AND {MOUSE} AND {SEQ}',
    f'EGFR AND lung AND (immun* OR "T cell" OR PD-1) AND {MOUSE} AND {SEQ}',
    f'(344SQ OR 344P OR 393P OR "KrasLA1" OR "Kras;p53 metastatic") AND {MOUSE}',
    f'CMT167 AND {MOUSE}',
    f'(Tacstd2 OR Trop2 OR Cldn4 OR claudin-4) AND lung AND {MOUSE}',
    f'lung adenocarcinoma AND (orthotopic OR autochthonous OR GEMM OR syngeneic) AND {MOUSE} AND {SEQ}',
]

ICI_PAT = re.compile(
    r"anti-?\s?pd-?\s?1|anti-?\s?pd-?l1|\bapd-?1\b|\bapdl1\b|anti-?\s?ctla-?4|"
    r"\bactla4?\b|pd-?1 blockade|pd-?l1 blockade|checkpoint blockade|checkpoint inhibit|"
    r"\bici\b|\bicb\b|pembrolizumab|nivolumab|atezolizumab|durvalumab|serplulimab|"
    r"ipilimumab|immunotherapy|\bitc\b",
    re.I,
)

MODEL_PATS = {
    "344SQ_line": r"\b344sq\b|\b344p\b|\b393p\b|\b531ln\b|\b344lm\b",
    "CMT167": r"\bcmt-?167\b|\bcmt167\b",
    "KP_kras_p53": r"kras.{0,25}(trp53|p53)|(trp53|p53).{0,25}kras|\bkpc?\b.{0,30}lung|kras.{0,10}g12d.{0,30}(trp53|p53)|\bkp\b",
    "KL_kras_lkb1": r"kras.{0,25}(lkb1|stk11)|\bkl\b.{0,20}(cell|tumor|line)|stk11",
    "SCLC_GEMM": r"\brpm\b|\brpp\b|\brpr2\b|rb1.{0,20}(trp53|p53)|small cell lung|\bsclc\b|myct58a",
    "EGFR": r"egfr.{0,30}(del19|l858r|l860r|exon 19|t790m|mutant)|\begfr\b",
    "LLC": r"lewis lung|\bllc1?\b",
}

INVIVO_PAT = re.compile(
    r"tumou?r tissue|tumou?r nodule|lung tumou?r|subcutaneous|orthotopic|flank|"
    r"autochthonous|tumou?r-bearing|in vivo|allograft|xenograft|implant", re.I
)


def _get(url: str, timeout: int = 120, tries: int = 5) -> bytes:
    last: Exception | None = None
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "opus-gemm/1.0"})
            with urllib.request.urlopen(req, timeout=timeout) as fh:
                return fh.read()
        except urllib.error.HTTPError as exc:
            if exc.code in (403, 404):
                raise
            last = exc
        except Exception as exc:
            last = exc
        time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"GET failed {url}: {last}")


def esearch(term: str, retmax: int = 500) -> list[str]:
    url = f"{EUTILS}/esearch.fcgi?" + urllib.parse.urlencode(
        {"db": "gds", "term": term, "retmax": retmax, "retmode": "json"}
    )
    return json.loads(_get(url, timeout=90).decode())["esearchresult"].get("idlist", [])


def esummary(uids: list[str]) -> list[dict]:
    out: list[dict] = []
    for i in range(0, len(uids), 100):
        chunk = uids[i : i + 100]
        url = f"{EUTILS}/esummary.fcgi?" + urllib.parse.urlencode(
            {"db": "gds", "id": ",".join(chunk), "retmode": "json"}
        )
        res = json.loads(_get(url, timeout=120).decode()).get("result", {})
        out.extend(res[u] for u in res.get("uids", []))
        time.sleep(0.4)
    return out


def series_dir(acc: str) -> str:
    stub = acc[:-3] + "nnn"
    return f"{FTP}/{stub}/{acc}"


def fetch_miniml(acc: str, cache: str) -> str | None:
    os.makedirs(cache, exist_ok=True)
    path = os.path.join(cache, f"{acc}_family.xml")
    if os.path.exists(path):
        with open(path, encoding="utf-8", errors="replace") as fh:
            return fh.read()
    try:
        blob = _get(f"{series_dir(acc)}/miniml/{acc}_family.xml.tgz", timeout=300, tries=3)
    except Exception as exc:
        print(f"  [{acc}] miniml unavailable: {exc}", file=sys.stderr)
        return None
    try:
        with tarfile.open(fileobj=io.BytesIO(blob), mode="r:gz") as tf:
            member = next((m for m in tf.getmembers() if m.name.endswith("_family.xml")), None)
            if member is None:
                return None
            xml = tf.extractfile(member).read().decode("utf-8", "replace")
    except Exception as exc:
        print(f"  [{acc}] miniml unreadable: {exc}", file=sys.stderr)
        return None
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(xml)
    return xml


def _first(xml: str, tag: str) -> str:
    m = re.search(rf"<{tag}[^>]*>(.*?)</{tag}>", xml, re.S)
    return re.sub(r"\s+", " ", m.group(1)).strip() if m else ""


def _series_block(xml: str) -> str:
    m = re.search(r"<Series iid=\"GSE\d+\">(.*?)</Series>", xml, re.S)
    return m.group(1) if m else xml


def parse_series(acc: str, xml: str) -> dict:
    samples = []
    for gsm, body in re.findall(r'<Sample iid="(GSM\d+)">(.*?)</Sample>', xml, re.S):
        chars = {
            (tag or "attr").strip(): re.sub(r"\s+", " ", val).strip()
            for tag, val in re.findall(
                r'<Characteristics(?: tag="([^"]*)")?>(.*?)</Characteristics>', body, re.S
            )
        }
        samples.append(
            {
                "gsm": gsm,
                "title": re.sub(r"\s+", " ", _first(body, "Title")),
                "source": re.sub(r"\s+", " ", _first(body, "Source")),
                "strategy": _first(body, "Library-Strategy"),
                "characteristics": chars,
                "organism": _first(body, "Organism"),
            }
        )
    series = _series_block(xml)
    return {
        "accession": acc,
        "series_title": _first(series, "Title"),
        "summary": _first(series, "Summary"),
        "design": _first(series, "Overall-Design"),
        "pubmed": re.findall(r"<Pubmed-ID>(\d+)</Pubmed-ID>", xml),
        "samples": samples,
    }


def score(rec: dict) -> dict:
    sample_text = " ".join(
        f"{s['title']} {s['source']} " + " ".join(f"{k} {v}" for k, v in s["characteristics"].items())
        for s in rec["samples"]
    )
    design_text = f"{rec['series_title']} {rec['summary']} {rec['design']}"
    all_text = design_text + " " + sample_text

    ici_in_samples = bool(ICI_PAT.search(sample_text))
    models = sorted(k for k, pat in MODEL_PATS.items() if re.search(pat, all_text, re.I))
    strategies = sorted({s["strategy"] for s in rec["samples"] if s["strategy"]})
    return {
        "ici_any": bool(ICI_PAT.search(all_text)),
        "ici_in_sample_metadata": ici_in_samples,
        "models": ",".join(models),
        "llc_only": models == ["LLC"],
        "in_vivo": bool(INVIVO_PAT.search(all_text)),
        "strategies": ",".join(strategies),
        "n_samples": len(rec["samples"]),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--json", required=True)
    ap.add_argument("--cache", default="/tmp/opus_gemm_miniml")
    ap.add_argument("--max-series", type=int, default=600)
    args = ap.parse_args()

    uids: list[str] = []
    for q in QUERIES:
        hits = esearch(q)
        print(f"[esearch] {len(hits):4d}  {q[:90]}", file=sys.stderr)
        uids.extend(hits)
        time.sleep(0.4)
    uids = list(dict.fromkeys(uids))

    summaries = esummary(uids)
    series = {}
    for s in summaries:
        acc = s.get("accession", "")
        if acc.startswith("GSE") and "Mus musculus" in s.get("taxon", ""):
            series[acc] = s
    print(f"[esummary] {len(series)} unique mouse GSE", file=sys.stderr)

    accs = sorted(series)[: args.max_series]
    records = {}
    for i, acc in enumerate(accs, 1):
        xml = fetch_miniml(acc, args.cache)
        if xml is None:
            continue
        rec = parse_series(acc, xml)
        rec["gds_summary_type"] = series[acc].get("gdsType", "")
        rec["suppfile"] = series[acc].get("suppFile", "")
        rec["score"] = score(rec)
        records[acc] = rec
        if i % 25 == 0:
            print(f"  parsed {i}/{len(accs)}", file=sys.stderr)

    with open(args.json, "w") as fh:
        json.dump(records, fh, indent=1)

    fields = [
        "accession", "n_samples", "ici_any", "ici_in_sample_metadata", "models",
        "llc_only", "in_vivo", "strategies", "gds_summary_type", "suppfile",
        "pubmed", "series_title", "design",
    ]
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for acc, rec in sorted(records.items()):
            row = {"accession": acc}
            row.update(rec["score"])
            row["gds_summary_type"] = rec["gds_summary_type"]
            row["suppfile"] = rec["suppfile"]
            row["pubmed"] = ",".join(rec["pubmed"])
            row["series_title"] = rec["series_title"][:200]
            row["design"] = rec["design"][:400]
            w.writerow({k: row.get(k, "") for k in fields})
    print(f"wrote {len(records)} series -> {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
