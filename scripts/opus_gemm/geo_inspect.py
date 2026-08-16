#!/usr/bin/env python3
"""Inspect specific GEO series: summary, overall design, sample titles and the
supplementary files (name + byte size) that ship processed data.

Used to verify model identity (KP / EGFR / SCLC / CMT167 / 344SQ ...), the
presence of a real immune-checkpoint-inhibitor arm, and whether processed
gene-level matrices are small enough to vendor into this slice.

Usage:
    python3 scripts/opus_gemm/geo_inspect.py GSE217405 GSE236258 ...
    python3 scripts/opus_gemm/geo_inspect.py --json out.json GSE...
"""
from __future__ import annotations

import argparse
import gzip
import io
import json
import re
import sys
import time
import urllib.error
import urllib.request

FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"
ACC_URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"


def _get(url: str, timeout: int = 120, tries: int = 5) -> bytes:
    last: Exception | None = None
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "opus-gemm/1.0"})
            with urllib.request.urlopen(req, timeout=timeout) as fh:
                return fh.read()
        except urllib.error.HTTPError as exc:
            if exc.code in (404, 403):
                raise
            last = exc
        except Exception as exc:
            last = exc
        time.sleep(2 ** attempt)
    raise RuntimeError(f"GET failed {url}: {last}")


def series_dir(acc: str) -> str:
    stub = acc[:-3] + "nnn" if len(acc) > 6 else "GSEnnn"
    return f"{FTP}/{stub}/{acc}"


def suppl_listing(acc: str) -> list[dict]:
    """Parse the Apache-style HTML index of the series suppl/ directory."""
    out: list[dict] = []
    for sub in ("suppl", "matrix"):
        url = f"{series_dir(acc)}/{sub}/"
        try:
            html = _get(url, timeout=60).decode("utf-8", "replace")
        except Exception:
            continue
        for m in re.finditer(
            r'<a href="([^"]+)">[^<]+</a>\s*(\d{4}-\d\d-\d\d \d\d:\d\d)\s+([0-9.]+[KMG]?)',
            html,
        ):
            name, _date, size = m.groups()
            if name.startswith("/") or name.endswith("/"):
                continue
            out.append({"dir": sub, "name": name, "size": size, "url": url + name})
    return out


def miniml(acc: str) -> dict:
    """Pull the family MINiML XML (small) for sample-level characteristics."""
    url = f"{series_dir(acc)}/miniml/{acc}_family.xml.tgz"
    try:
        blob = _get(url, timeout=300)
    except Exception as exc:
        return {"error": str(exc)}
    import tarfile

    info: dict = {"samples": []}
    with tarfile.open(fileobj=io.BytesIO(blob), mode="r:gz") as tf:
        member = next((m for m in tf.getmembers() if m.name.endswith("_family.xml")), None)
        if member is None:
            return {"error": "no family.xml"}
        xml = tf.extractfile(member).read().decode("utf-8", "replace")
    info["title"] = _tag(xml, "Title")
    info["summary"] = _tag(xml, "Summary")
    info["design"] = _tag(xml, "Overall-Design")
    for block in re.findall(r'<Sample iid="(GSM\d+)">(.*?)</Sample>', xml, re.S):
        gsm, body = block
        chars = {
            (tag or "attr").strip(): re.sub(r"\s+", " ", val).strip()
            for tag, val in re.findall(
                r'<Characteristics(?: tag="([^"]*)")?>(.*?)</Characteristics>', body, re.S
            )
        }
        info["samples"].append(
            {
                "gsm": gsm,
                "title": re.sub(r"\s+", " ", _tag(body, "Title") or "").strip(),
                "source": re.sub(r"\s+", " ", _tag(body, "Source") or "").strip(),
                "library_strategy": _tag(body, "Library-Strategy"),
                "characteristics": chars,
                "supplementary": re.findall(r"<Supplementary-Data[^>]*>\s*([^<]+?)\s*</Supplementary-Data>", body),
            }
        )
    return info


def _tag(xml: str, tag: str) -> str | None:
    m = re.search(rf"<{tag}[^>]*>(.*?)</{tag}>", xml, re.S)
    return re.sub(r"\s+", " ", m.group(1)).strip() if m else None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("accessions", nargs="+")
    ap.add_argument("--json", help="write full records as JSON")
    ap.add_argument("--no-samples", action="store_true")
    args = ap.parse_args()

    records = {}
    for acc in args.accessions:
        print(f"\n{'=' * 100}\n{acc}", flush=True)
        rec: dict = {"accession": acc}
        meta = miniml(acc)
        rec["meta"] = meta
        if "error" in meta:
            print(f"  MINiML error: {meta['error']}")
        else:
            print(f"  TITLE   : {meta.get('title')}")
            print(f"  DESIGN  : {(meta.get('design') or '')[:600]}")
            print(f"  SUMMARY : {(meta.get('summary') or '')[:700]}")
            print(f"  N SAMPLE: {len(meta['samples'])}")
            if not args.no_samples:
                for s in meta["samples"]:
                    ch = "; ".join(f"{k}={v}" for k, v in s["characteristics"].items())
                    print(f"    {s['gsm']} [{s.get('library_strategy')}] {s['title']} || {ch}")
        sup = suppl_listing(acc)
        rec["suppl"] = sup
        print("  SUPPL:")
        for f in sup:
            print(f"    {f['size']:>8}  {f['dir']}/{f['name']}")
        records[acc] = rec
        time.sleep(0.5)

    if args.json:
        with open(args.json, "w") as fh:
            json.dump(records, fh, indent=1)
        print(f"\nwrote {args.json}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
