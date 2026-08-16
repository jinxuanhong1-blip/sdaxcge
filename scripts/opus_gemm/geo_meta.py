#!/usr/bin/env python3
"""GEO MINiML metadata: fetch, parse, and map matrix columns to GSM samples.

Every expression matrix used in this slice is keyed by GSM accession, and the
mapping from a vendor-specific column name to a GSM is resolved by one of three
auditable strategies (recorded per sample as ``column_source``):

  explicit_description : GEO sample description literally states the column name
                         ("Column name in <file>: <col>")
  supplementary_file   : the per-sample supplementary file name carries the GSM
  title_match          : normalised sample title matches the column name

Anything that cannot be resolved raises, so a silent mis-assignment cannot slip
into the statistics.
"""
from __future__ import annotations

import io
import os
import re
import tarfile
import time
import urllib.error
import urllib.request

FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"
CACHE = os.environ.get("OPUS_GEMM_MINIML", "results/opus_gemm/data/miniml")


def _get(url: str, timeout: int = 300, tries: int = 5) -> bytes:
    last: Exception | None = None
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "opus-gemm/1.0"})
            with urllib.request.urlopen(req, timeout=timeout) as fh:
                return fh.read()
        except urllib.error.HTTPError:
            raise
        except Exception as exc:
            last = exc
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"GET failed {url}: {last}")


def family_xml(acc: str, cache: str = CACHE) -> str:
    os.makedirs(cache, exist_ok=True)
    path = os.path.join(cache, f"{acc}_family.xml")
    if os.path.exists(path):
        return open(path, encoding="utf-8", errors="replace").read()
    blob = _get(f"{FTP}/{acc[:-3]}nnn/{acc}/miniml/{acc}_family.xml.tgz")
    with tarfile.open(fileobj=io.BytesIO(blob), mode="r:gz") as tf:
        member = next(m for m in tf.getmembers() if m.name.endswith("_family.xml"))
        xml = tf.extractfile(member).read().decode("utf-8", "replace")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(xml)
    return xml


def _tag(text: str, tag: str) -> str:
    m = re.search(rf"<{tag}[^>]*>(.*?)</{tag}>", text, re.S)
    return re.sub(r"\s+", " ", m.group(1)).strip() if m else ""


def series_info(acc: str) -> dict:
    xml = family_xml(acc)
    block = re.search(r"<Series iid=\"GSE\d+\">(.*?)</Series>", xml, re.S)
    body = block.group(1) if block else xml
    return {
        "accession": acc,
        "title": _tag(body, "Title"),
        "summary": _tag(body, "Summary"),
        "design": _tag(body, "Overall-Design"),
        "pubmed": re.findall(r"<Pubmed-ID>(\d+)</Pubmed-ID>", body),
    }


def samples(acc: str) -> list[dict]:
    xml = family_xml(acc)
    out: list[dict] = []
    for gsm, body in re.findall(r'<Sample iid="(GSM\d+)">(.*?)</Sample>', xml, re.S):
        chars = {
            (tag or "attr").strip(): re.sub(r"\s+", " ", val).strip()
            for tag, val in re.findall(
                r'<Characteristics(?: tag="([^"]*)")?>(.*?)</Characteristics>', body, re.S
            )
        }
        descriptions = [
            d.strip() for d in re.findall(r"<Description>(.*?)</Description>", body, re.S)
        ]
        out.append(
            {
                "gsm": gsm,
                "title": _tag(body, "Title"),
                "source": _tag(body, "Source"),
                "strategy": _tag(body, "Library-Strategy"),
                "characteristics": chars,
                "descriptions": descriptions,
                "supplementary": [
                    s.strip().split("/")[-1]
                    for s in re.findall(r"<Supplementary-Data[^>]*>(.*?)</Supplementary-Data>", body, re.S)
                    if s.strip() not in ("", "NONE")
                ],
            }
        )
    return out


def _norm(s: str) -> str:
    return re.sub(r"[^a-z0-9]", "", s.lower())


def declared_column(sample: dict) -> str | None:
    """Column name if the GEO description states it explicitly."""
    for d in sample["descriptions"]:
        m = re.search(r"[Cc]olumn name in [^:]*:\s*([^\s;]+)", d)
        if m:
            return m.group(1).strip()
    return None


def resolve_columns(acc: str, columns: list[str], strict: bool = True) -> dict[str, dict]:
    """Map matrix column name -> sample record (with ``column_source`` set).

    ``columns`` are the data columns of the matrix in file order.
    """
    smps = samples(acc)
    remaining = {c: _norm(c) for c in columns}
    mapping: dict[str, dict] = {}

    # 1. explicit column names from the GEO description
    for s in smps:
        col = declared_column(s)
        if col and col in remaining:
            mapping[col] = dict(s, column_source="explicit_description")
            remaining.pop(col)

    # 2. per-sample supplementary file names carrying the GSM id
    for s in smps:
        if s["gsm"] in [m["gsm"] for m in mapping.values()]:
            continue
        for col in list(remaining):
            if s["gsm"] in col:
                mapping[col] = dict(s, column_source="supplementary_file")
                remaining.pop(col)
                break

    # 3. normalised title / supplementary-stem matching
    for s in smps:
        if s["gsm"] in [m["gsm"] for m in mapping.values()]:
            continue
        cands = [s["title"]] + [re.sub(r"\.[^.]*$", "", f) for f in s["supplementary"]]
        keys = {_norm(c) for c in cands if c}
        hit = None
        for col, ncol in remaining.items():
            if ncol in keys or any(ncol == k or k.endswith(ncol) or ncol.endswith(k) for k in keys):
                hit = col
                break
        if hit:
            mapping[hit] = dict(s, column_source="title_match")
            remaining.pop(hit)

    if strict and remaining:
        raise RuntimeError(
            f"{acc}: could not map columns {sorted(remaining)}; "
            f"unmapped samples: {[s['gsm'] + '/' + s['title'] for s in smps if s['gsm'] not in {m['gsm'] for m in mapping.values()}]}"
        )
    return mapping
