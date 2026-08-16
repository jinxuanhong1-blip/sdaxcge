#!/usr/bin/env python3
"""Second sweep: lung/ICI series whose response labels live in a supplementary file, not in GEO metadata.

04/06 only see `!Sample_characteristics_ch1`. A number of submitters instead ship a small
"clinical"/"metadata"/"sample info" table alongside the expression matrix. This script lists those
candidates and, for small text payloads, downloads and greps them for response vocabulary so no
usable cohort is missed for a purely bureaucratic reason.
"""
from __future__ import annotations

import gzip
import io
import json
import re
import tarfile
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "opus_geo_leftover"
FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series"
SCAN_LIMIT = 30 * 1024**2  # only peek inside files small enough to read cheaply

LUNG = re.compile(r"\blung\b|NSCLC|non[- ]small[- ]cell|\bSCLC\b|LUAD|LUSC", re.I)
ICI = re.compile(
    r"nivolumab|pembrolizumab|atezolizumab|durvalumab|avelumab|ipilimumab|cemiplimab|"
    r"camrelizumab|sintilimab|tislelizumab|toripalimab|anti-?PD-?1|anti-?PD-?L1|anti-?PD-?\(L\)1|"
    r"immunotherap|immune checkpoint|checkpoint (inhibit|block)|\bICI\b|\bICB\b|PD-?1/PD-?L1",
    re.I,
)
CLINICAL_NAME = re.compile(
    r"clinic|metadata|meta_?data|sample.?(info|sheet|key|annot)|pheno|patient|cohort|"
    r"response|annotation|design|characteristic|coldata", re.I
)
RESPONSE_VALUES = re.compile(
    r"\b(responder|non-?responder|\bPR\b|\bCR\b|\bSD\b|\bPD\b|RECIST|MPR|pCR|DCB|NDB|"
    r"partial response|complete response|stable disease|progressive disease|durable)\b"
)
TEXTY = re.compile(r"\.(txt|tsv|csv|xlsx?|json)(\.gz)?$", re.I)


def fetch(url: str, limit: int | None = None) -> bytes | None:
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "opus-geo-leftover/1.0"})
        with urllib.request.urlopen(req, timeout=300) as fh:
            return fh.read(limit) if limit else fh.read()
    except Exception:  # noqa: BLE001
        return None


def peek(blob: bytes, name: str) -> str:
    if name.endswith(".gz"):
        try:
            blob = gzip.decompress(blob)
        except OSError:
            try:
                blob = gzip.GzipFile(fileobj=io.BytesIO(blob)).read(4 << 20)
            except Exception:  # noqa: BLE001
                return ""
    if name.endswith(".tar"):
        try:
            tf = tarfile.open(fileobj=io.BytesIO(blob))
            return " ".join(m.name for m in tf.getmembers())
        except Exception:  # noqa: BLE001
            return ""
    return blob[: 4 << 20].decode("utf-8", "replace")


def main() -> int:
    probe = json.loads((OUT / "geo_probe.json").read_text())
    cands = {c["accession"]: c for c in json.loads((OUT / "geo_triage_scored.json").read_text())}

    findings = []
    for r in probe:
        acc = r["accession"]
        cand = cands.get(acc, {})
        blob = " ".join([r["title"], cand.get("summary", ""), " ".join(r["characteristics"])])
        if not (LUNG.search(blob) and ICI.search(blob)):
            continue
        interesting = [
            f
            for f in r["supp_files"]
            if CLINICAL_NAME.search(f["name"]) and TEXTY.search(f["name"])
        ]
        if not interesting:
            continue
        for f in interesting:
            size = f["bytes"] or 0
            hit_text = ""
            if size and size <= SCAN_LIMIT and not f["name"].endswith((".xlsx", ".xls")):
                url = f"{FTP_BASE}/{stub(acc)}/{acc}/suppl/{urllib.parse.quote(f['name'])}"
                raw = fetch(url, limit=SCAN_LIMIT)
                if raw:
                    hit_text = peek(raw, f["name"])
            matches = sorted(set(RESPONSE_VALUES.findall(hit_text)))
            findings.append(
                {
                    "accession": acc,
                    "title": r["title"],
                    "file": f["name"],
                    "bytes": size,
                    "response_tokens_found": matches,
                    "preview": hit_text[:400].replace("\n", " | ") if hit_text else "",
                }
            )
            flag = "*" if matches else " "
            print(f"{flag} {acc:<11} {f['name']:<52} {size:>10} {matches[:8]}")

    (OUT / "supp_clinical_scan.json").write_text(json.dumps(findings, indent=2))
    print(f"\n[scan] {len(findings)} candidate clinical/annotation files")
    return 0


def stub(acc: str) -> str:
    n = acc[3:]
    return f"GSE{n[:-3]}nnn" if len(n) > 3 else "GSEnnn"


if __name__ == "__main__":
    raise SystemExit(main())
