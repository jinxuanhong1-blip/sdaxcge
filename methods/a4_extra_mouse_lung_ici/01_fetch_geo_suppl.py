#!/usr/bin/env python3
"""Fetch GEO series pages + supplementary file lists for extra mouse-lung ICI candidates."""
from __future__ import annotations

import json
import re
import time
import urllib.request
from pathlib import Path

OUT = Path("notes/a4_extra_mouse_lung_ici/raw")
OUT.mkdir(parents=True, exist_ok=True)

# Extra (non-TISMO) mouse lung / NSCLC models with PD-1/PD-L1/CTLA4 or ICB-adjacent RNA.
# GSE155972 is TISMO LLC ICB and is excluded.
CANDIDATES = [
    # LLC / CMT-167 / 344SQ syngeneic
    "GSE239485",
    "GSE297630",
    "GSE297632",
    "GSE330941",
    "GSE333285",
    "GSE241978",
    "GSE241979",
    "GSE241980",
    "GSE190264",
    "GSE209766",
    "GSE260596",
    "GSE317011",
    "GSE304155",
    # KP / KL / STK11 / LKB1 / GEMM / orthotopic NSCLC
    "GSE179500",  # known KL genotype; cite, not ICB
    "GSE179501",
    "GSE180963",
    "GSE137396",
    "GSE137244",
    "GSE182228",
    "GSE244452",
    "GSE114601",
    "GSE114300",
    "GSE169194",
    "GSE169196",
    "GSE246922",
    "GSE236258",
    "GSE298051",
    "GSE295685",
    "GSE256071",
    "GSE275877",
    "GSE262305",
    "GSE227534",
    "GSE233862",
    "GSE261898",
    "GSE274960",
    "GSE268555",
    "GSE285341",
    "GSE285342",
    "GSE157880",
    "GSE157883",
    "GSE178177",
    "GSE197260",
    "GSE133604",
    "GSE274351",
    "GSE274352",
    "GSE338923",
    "GSE330733",
    "GSE214613",
    "GSE214485",
    "GSE210547",
    "GSE165517",
    "GSE193895",
    "GSE180964",
    "GSE271713",
    "GSE266364",
    "GSE277610",
    "GSE300288",
    "GSE291687",
    "GSE303940",
    "GSE303943",
    "GSE288231",
    "GSE276724",
    "GSE285606",
    "GSE209588",
]


def fetch(url: str, dest: Path, timeout: int = 90) -> bool:
    dest.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": "a4-extra-mouse-lung-ici/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            dest.write_bytes(r.read())
        return True
    except Exception as e:
        dest.write_text(f"ERROR {e}\nURL {url}\n")
        return False


def parse_suppl(html: str) -> list[dict]:
    files = []
    # GEO supplementary table rows
    for m in re.finditer(
        r'href="(ftp://ftp\.ncbi\.nlm\.nih\.gov/geo/series/[^"]+|https://www\.ncbi\.nlm\.nih\.gov/geo/download/\?acc=[^"]+)"[^>]*>([^<]+)</a>',
        html,
    ):
        files.append({"url": m.group(1).replace("&amp;", "&"), "name": m.group(2).strip()})
    # also file names in (txt|tsv|csv|xlsx) listings
    for m in re.finditer(r"(GSE\d+[A-Za-z0-9_\-\.]+?\.(txt|tsv|csv|xlsx|xls|gz|zip))", html):
        files.append({"name": m.group(1), "ext": m.group(2)})
    # unique by name
    seen = set()
    out = []
    for f in files:
        n = f.get("name")
        if n and n not in seen:
            seen.add(n)
            out.append(f)
    return out


def main() -> None:
    catalog = []
    for acc in CANDIDATES:
        html_path = OUT / f"{acc}.html"
        url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=self&form=html&view=full"
        ok = True
        if not html_path.exists() or html_path.stat().st_size < 200:
            ok = fetch(url, html_path)
            time.sleep(0.35)
        html = html_path.read_text(errors="replace") if html_path.exists() else ""
        title = ""
        m = re.search(r"<td[^>]*>Title</td>\s*<td[^>]*>(.*?)</td>", html, re.S | re.I)
        if m:
            title = re.sub("<[^>]+>", "", m.group(1)).strip()
        n_samples = None
        m = re.search(r"Samples\s*\((\d+)\)", html)
        if m:
            n_samples = int(m.group(1))
        organism = ""
        m = re.search(r"<td[^>]*>Organism</td>\s*<td[^>]*>(.*?)</td>", html, re.S | re.I)
        if m:
            organism = re.sub("<[^>]+>", " ", m.group(1))
            organism = re.sub(r"\s+", " ", organism).strip()
        suppl = parse_suppl(html)
        catalog.append(
            {
                "accession": acc,
                "fetch_ok": ok and "ERROR" not in html[:40],
                "title": title[:240],
                "n_samples": n_samples,
                "organism": organism[:200],
                "html_bytes": html_path.stat().st_size if html_path.exists() else 0,
                "suppl": suppl,
            }
        )
        print(f"{acc}\tn={n_samples}\tsuppl={len(suppl)}\t{title[:90]}")

    (OUT / "suppl_catalog.json").write_text(json.dumps(catalog, indent=2))
    print("wrote", OUT / "suppl_catalog.json", "n=", len(catalog))


if __name__ == "__main__":
    main()
