#!/usr/bin/env python3
"""Shared helpers: fetch GEO platform annotations and build probe->symbol maps."""
import gzip
import io
import re
import time
import urllib.request
from pathlib import Path

ANNOT_DIR = Path("results/geo_sweep/annot")
ANNOT_DIR.mkdir(parents=True, exist_ok=True)


def _get(url, tries=4):
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "geo-sweep/1.0"})
            with urllib.request.urlopen(req, timeout=180) as r:
                return r.read()
        except Exception as e:  # noqa: BLE001
            last = e
            code = getattr(e, "code", None)
            if code == 404:
                return None
            time.sleep(2 * (i + 1))
    raise RuntimeError(f"GET failed {url}: {last}")


def _annot_gz_url(gpl):
    stub = gpl[:-3] + "nnn"
    return f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/{stub}/{gpl}/annot/{gpl}.annot.gz"


def _self_full_url(gpl):
    return (f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gpl}"
            f"&targ=self&form=text&view=full")


def _download_platform_text(gpl):
    """Return the platform annotation text, cached locally."""
    cache = ANNOT_DIR / f"{gpl}.txt"
    if cache.exists():
        return cache.read_text(encoding="utf-8", errors="replace")
    # try curated annot.gz first
    raw = _get(_annot_gz_url(gpl))
    if raw is not None:
        text = gzip.decompress(raw).decode("utf-8", "replace")
    else:
        raw = _get(_self_full_url(gpl))
        text = raw.decode("utf-8", "replace") if raw else ""
    cache.write_text(text)
    return text


def _split_table(text):
    """Yield rows (as lists) of the platform_table / annotation table block."""
    lines = text.splitlines()
    start = None
    for i, ln in enumerate(lines):
        low = ln.lower()
        if low.startswith("!platform_table_begin"):
            start = i + 1
            break
    if start is None:
        # curated .annot files have no marker; header begins at first non-comment
        for i, ln in enumerate(lines):
            if ln and not ln.startswith(("^", "!", "#")):
                start = i
                break
    if start is None:
        return
    for ln in lines[start:]:
        if ln.lower().startswith("!platform_table_end"):
            break
        if not ln.strip():
            continue
        yield ln.split("\t")


SYMBOL_COLS = [
    "gene symbol", "gene_symbol", "symbol", "genesymbol",
    "gene", "ilmn_gene",
]


def build_symbol_map(gpl):
    """Return dict probe_id(str) -> uppercase gene symbol for the platform."""
    text = _download_platform_text(gpl)
    rows = list(_split_table(text))
    if not rows:
        return {}
    header = [h.strip().strip('"') for h in rows[0]]
    lower = [h.lower() for h in header]

    def find_col(cands):
        for c in cands:
            if c in lower:
                return lower.index(c)
        return None

    id_idx = 0
    sym_idx = find_col(SYMBOL_COLS)
    ga_idx = None
    for i, h in enumerate(lower):
        if h in ("gene_assignment", "gene assignment"):
            ga_idx = i
            break
    pname_idx = lower.index("probe name") if "probe name" in lower else None

    mp = {}
    for row in rows[1:]:
        if not row:
            continue
        pid = row[id_idx].strip().strip('"')
        sym = None
        if sym_idx is not None and sym_idx < len(row):
            sym = row[sym_idx].strip().strip('"')
        if (not sym or sym in ("", "---", "N/A")) and ga_idx is not None and ga_idx < len(row):
            # Clariom gene_assignment: "acc // SYMBOL // desc // ..."
            ga = row[ga_idx]
            m = re.search(r"//\s*([A-Za-z0-9\-\.]+)\s*//", ga)
            if m:
                sym = m.group(1)
        if (not sym or sym in ("", "---", "N/A")) and pname_idx is not None and pname_idx < len(row):
            sym = row[pname_idx].strip().strip('"')
        if sym and sym not in ("", "---", "N/A"):
            mp[pid] = sym.upper()
    return mp
