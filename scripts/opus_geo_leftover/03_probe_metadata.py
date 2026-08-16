#!/usr/bin/env python3
"""Probe GEO FTP for each candidate series: supplementary file inventory + sample characteristics.

For every candidate we record
  * the supplementary files published at the series level (name + byte size),
  * whether a processed expression matrix appears to be present,
  * whether the series is FASTQ/raw-only (i.e. must be skipped per the >2GB / FASTQ rule),
  * the free-text sample characteristics, so response labels can be detected.

Series metadata comes from the small `*_series_matrix.txt.gz` header, which carries the
`!Sample_characteristics_ch1` lines for both array and RNA-seq series.
"""
from __future__ import annotations

import gzip
import io
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from html.parser import HTMLParser
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "opus_geo_leftover"
CACHE = Path("/tmp/geo_dl/meta")
CACHE.mkdir(parents=True, exist_ok=True)

FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series"
MAX_BYTES = 2 * 1024**3  # 2 GB hard ceiling from the task brief

RESPONSE_PAT = re.compile(
    r"\bresponse\b|\bresponder\b|non-?responder|\bNR\b|\bDCB\b|\bNDB\b|\bPD\b|\bSD\b|\bPR\b|\bCR\b|"
    r"RECIST|best overall|\bBOR\b|clinical benefit|durable|\bMPR\b|pathologic|\bpCR\b|\bPFS\b|"
    r"progression.free|overall survival|\bOS\b|efficacy|resistan|treatment.?outcome",
    re.I,
)
TREAT_PAT = re.compile(
    r"nivolumab|pembrolizumab|atezolizumab|durvalumab|avelumab|ipilimumab|cemiplimab|"
    r"camrelizumab|sintilimab|tislelizumab|toripalimab|anti-?PD-?1|anti-?PD-?L1|"
    r"immunotherap|checkpoint|\bICI\b|\bICB\b",
    re.I,
)
RAW_ONLY_PAT = re.compile(r"\.(fastq|fq|bam|sra|cram|bcl|cel|idat)(\.gz|\.bz2|\.tar)?$", re.I)
MATRIX_PAT = re.compile(
    r"(count|cpm|tpm|fpkm|rpkm|rsem|expr|expression|matrix|normali[sz]ed|rma|gene|rnaseq|"
    r"nanostring|log2|abundance|quant)",
    re.I,
)


class LinkParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.links: list[str] = []

    def handle_starttag(self, tag: str, attrs) -> None:
        if tag == "a":
            for k, v in attrs:
                if k == "href" and v:
                    self.links.append(v)


def fetch(url: str, binary: bool = False, retries: int = 4):
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "opus-geo-leftover/1.0"})
            with urllib.request.urlopen(req, timeout=180) as fh:
                data = fh.read()
            return data if binary else data.decode("utf-8", "replace")
        except urllib.error.HTTPError as exc:
            if exc.code in (403, 404):
                return None
            last = exc
        except Exception as exc:  # noqa: BLE001
            last = exc
        time.sleep(2 ** attempt)
    print(f"  ! fetch failed {url}: {last}", file=sys.stderr)
    return None


def stub(acc: str) -> str:
    n = acc[3:]
    return f"GSE{n[:-3]}nnn" if len(n) > 3 else "GSEnnn"


def head_size(url: str) -> int | None:
    """Content-Length for a GEO supplementary file, so the >2 GB rule can be enforced exactly."""
    for attempt in range(3):
        try:
            req = urllib.request.Request(
                url, method="HEAD", headers={"User-Agent": "opus-geo-leftover/1.0"}
            )
            with urllib.request.urlopen(req, timeout=90) as fh:
                cl = fh.headers.get("Content-Length")
                return int(cl) if cl else None
        except urllib.error.HTTPError:
            return None
        except Exception:  # noqa: BLE001
            time.sleep(1 + attempt)
    return None


def list_dir(url: str, with_sizes: bool = False) -> list[tuple[str, int | None]]:
    """Return (filename, size_bytes_or_None) for an NCBI FTP-over-HTTPS index page."""
    html = fetch(url)
    if html is None:
        return []
    p = LinkParser()
    p.feed(html)
    names: list[str] = []
    for href in p.links:
        if href.endswith("/") or href.startswith(("?", "/", "http://", "https://", "#")):
            continue
        name = urllib.parse.unquote(href)
        if name not in names:
            names.append(name)
    if not with_sizes:
        return [(n, None) for n in names]
    return [(n, head_size(url + urllib.parse.quote(n))) for n in names]


def series_matrix_meta(acc: str) -> dict:
    base = f"{FTP_BASE}/{stub(acc)}/{acc}/matrix/"
    files = [f for f, _ in list_dir(base) if f.endswith("series_matrix.txt.gz")]
    chars: list[str] = []
    header: dict[str, list[str]] = {}
    n_data_rows = 0
    for fname in files[:3]:
        cached = CACHE / fname
        if cached.exists():
            blob = cached.read_bytes()
        else:
            blob = fetch(base + fname, binary=True)
            if blob is None:
                continue
            cached.write_bytes(blob)
        try:
            text = gzip.decompress(blob).decode("utf-8", "replace")
        except OSError:
            continue
        in_table = False
        for line in io.StringIO(text):
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                in_table = False
                continue
            if in_table:
                n_data_rows += 1
                continue
            if line.startswith("!Sample_characteristics_ch"):
                chars.append(line.strip())
            elif line.startswith("!"):
                key, _, val = line[1:].partition("\t")
                header.setdefault(key.strip(), []).append(val.strip().strip('"'))
    return {
        "matrix_files": files,
        "characteristics": chars,
        "header": header,
        "matrix_data_rows": max(0, n_data_rows - 1),
    }


def main() -> int:
    cands = json.loads((OUT / "geo_triage_scored.json").read_text())
    accs = [c["accession"] for c in cands]

    probe_path = OUT / "geo_probe.json"
    results = json.loads(probe_path.read_text()) if probe_path.exists() else []
    done = {r["accession"] for r in results}
    accs = [a for a in accs if a not in done]
    print(f"[probe] {len(done)} cached, probing {len(accs)} new series", file=sys.stderr)

    for i, acc in enumerate(accs, 1):
        suppl = list_dir(f"{FTP_BASE}/{stub(acc)}/{acc}/suppl/", with_sizes=True)
        meta = series_matrix_meta(acc)
        chars_blob = " ".join(meta["characteristics"])
        header = meta["header"]
        title = (header.get("Series_title") or [""])[0]
        summary = " ".join(header.get("Series_summary") or [])
        overall = " ".join(header.get("Series_overall_design") or [])
        blob = f"{title} {summary} {overall} {chars_blob}"

        processed = [
            (f, s)
            for f, s in suppl
            if not RAW_ONLY_PAT.search(f) and (MATRIX_PAT.search(f) or f.endswith((".txt.gz", ".csv.gz", ".tsv.gz", ".xlsx", ".txt", ".csv", ".tsv")))
        ]
        oversized = [(f, s) for f, s in processed if s is not None and s > MAX_BYTES]
        usable = [(f, s) for f, s in processed if s is None or s <= MAX_BYTES]

        rec = {
            "accession": acc,
            "title": title,
            "n_supp_files": len(suppl),
            "supp_files": [{"name": f, "bytes": s} for f, s in suppl],
            "processed_candidates": [{"name": f, "bytes": s} for f, s in usable],
            "oversized_files": [{"name": f, "bytes": s} for f, s in oversized],
            "raw_only": bool(suppl) and all(RAW_ONLY_PAT.search(f) for f, _ in suppl),
            "has_series_matrix_values": meta["matrix_data_rows"] > 100,
            "matrix_data_rows": meta["matrix_data_rows"],
            "matrix_files": meta["matrix_files"],
            "n_characteristics_lines": len(meta["characteristics"]),
            "characteristics": meta["characteristics"],
            "has_response_annotation": bool(RESPONSE_PAT.search(chars_blob)),
            "mentions_ici": bool(TREAT_PAT.search(blob)),
            "platform_ids": header.get("Series_platform_id", []),
            "sample_count": len(header.get("Sample_geo_accession", [])),
            "pubmed": header.get("Series_pubmed_id", []),
        }
        results.append(rec)
        flag = "*" if rec["has_response_annotation"] and rec["mentions_ici"] else " "
        print(
            f"[{i}/{len(accs)}]{flag} {acc} supp={len(suppl)} proc={len(usable)} "
            f"resp={int(rec['has_response_annotation'])} ici={int(rec['mentions_ici'])} "
            f"smx_rows={rec['matrix_data_rows']}",
            file=sys.stderr,
        )
        probe_path.write_text(json.dumps(results, indent=2))
    probe_path.write_text(json.dumps(results, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
