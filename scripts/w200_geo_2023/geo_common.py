#!/usr/bin/env python3
"""Shared NCBI/GEO helpers for the w200 GEO-2023 leftover sweep."""
import gzip
import io
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
GEO_FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"
TOOL = "w200_geo_2023"
EMAIL = "cursoragent@cursor.com"

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "w200" / "GEO_2023"
NOTES = ROOT / "notes" / "w200_geo_2023"
CACHE = Path("/tmp/w200_cache")
for _d in (OUT, NOTES, CACHE):
    _d.mkdir(parents=True, exist_ok=True)

GENES = ["TACSTD2", "CLDN4"]
# Aliases as they appear in author-supplied matrices (symbol or Ensembl id).
GENE_ALIASES = {
    "TACSTD2": {"TACSTD2", "TROP2", "TROP-2", "GA733-1", "M1S1", "EGP-1",
                "ENSG00000184292"},
    "CLDN4": {"CLDN4", "CPE-R", "CPER", "WBSCR8", "CLDN-4",
              "ENSG00000189143"},
}


def _request(url, retries=5, timeout=120):
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(
                url, headers={"User-Agent": f"{TOOL} ({EMAIL})"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code in (400, 404):
                raise
            last = e
        except Exception as e:  # noqa: BLE001
            last = e
        time.sleep(2 ** attempt)
    raise RuntimeError(f"request failed: {url}: {last}")


def eutils(endpoint, **params):
    params = {**params, "tool": TOOL, "email": EMAIL}
    url = f"{EUTILS}/{endpoint}?" + urllib.parse.urlencode(params)
    time.sleep(0.35)  # stay under the 3 req/s anonymous limit
    return _request(url).decode("utf-8", "replace")


def esearch(term, db="gds", retmax=2000):
    txt = eutils("esearch.fcgi", db=db, term=term, retmax=retmax,
                 retmode="json")
    return json.loads(txt)["esearchresult"].get("idlist", [])


def esummary(uids, db="gds"):
    out = {}
    uids = list(uids)
    for i in range(0, len(uids), 200):
        chunk = uids[i:i + 200]
        txt = eutils("esummary.fcgi", db=db, id=",".join(chunk),
                     retmode="json")
        res = json.loads(txt).get("result", {})
        for uid in res.get("uids", []):
            out[uid] = res[uid]
    return out


def series_dir(gse):
    """GEO FTP directory for a series, e.g. GSE207422 -> GSE207nnn/GSE207422."""
    stub = gse[:-3] + "nnn" if len(gse) > 6 else gse[:3] + "nnn"
    return f"{GEO_FTP}/{stub}/{gse}"


def cached_get(url, name, max_bytes=None):
    """Download to /tmp cache once; return path or None on failure."""
    dest = CACHE / name
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    try:
        if max_bytes is not None:
            head = head_size(url)
            if head is not None and head > max_bytes:
                return None
        data = _request(url, timeout=600)
    except Exception:
        return None
    dest.write_bytes(data)
    return dest


def head_size(url):
    try:
        req = urllib.request.Request(
            url, method="HEAD",
            headers={"User-Agent": f"{TOOL} ({EMAIL})"})
        with urllib.request.urlopen(req, timeout=60) as r:
            n = r.headers.get("Content-Length")
            return int(n) if n else None
    except Exception:
        return None


def list_suppl(gse):
    """List supplementary files (name, bytes) from the GEO FTP HTML index."""
    url = f"{series_dir(gse)}/suppl/"
    try:
        html = _request(url, retries=3, timeout=90).decode("utf-8", "replace")
    except Exception:
        return []
    import re
    files = []
    for m in re.finditer(r'href="([^"?/][^"]*)"', html):
        fn = m.group(1)
        if fn in ("filelist.txt",) or fn.startswith("http"):
            continue
        files.append(fn)
    # sizes from the pre-formatted listing
    sizes = {}
    for m in re.finditer(
            r'href="([^"?/][^"]*)"[^\n]*?\s(\d{4}-\d\d-\d\d \d\d:\d\d)\s+([0-9.]+[KMG]?)',
            html):
        sizes[m.group(1)] = m.group(3)
    return [(f, sizes.get(f)) for f in sorted(set(files))]


def open_maybe_gz(path):
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def fetch_series_matrix(gse):
    """Download and return the series_matrix text (header + table) or None."""
    for suffix in ("", "-1"):
        url = (f"{series_dir(gse)}/matrix/{gse}{suffix}_series_matrix.txt.gz")
        p = cached_get(url, f"{gse}{suffix}_series_matrix.txt.gz")
        if p:
            try:
                with gzip.open(p, "rt", encoding="utf-8", errors="replace") as fh:
                    return fh.read()
            except Exception:
                continue
    return None


def parse_series_matrix_header(text):
    """Return dict of !Series_* and !Sample_* header lines."""
    meta = {}
    for line in text.splitlines():
        if not line.startswith("!"):
            if line.startswith('"ID_REF"'):
                break
            continue
        parts = line.rstrip("\n").split("\t")
        key = parts[0].lstrip("!")
        vals = [p.strip('"') for p in parts[1:]]
        meta.setdefault(key, []).append(vals if len(vals) > 1 else
                                        (vals[0] if vals else ""))
    return meta
