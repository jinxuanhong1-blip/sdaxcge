#!/usr/bin/env python3
"""Download GEO series_matrix files (metadata only, small) for the shortlist
and extract per-sample characteristics so we can see which series carry ICI
response labels. Writes one *.characteristics.txt per GSE and a combined
results/w200/GEO_2026/label_scan.txt summary.
"""
import gzip
import io
import time
import urllib.request
from pathlib import Path

RES = Path(__file__).resolve().parents[2] / "results" / "w200" / "GEO_2026"
META = RES / "series_matrix"
META.mkdir(exist_ok=True)

SHORTLIST = [
    "GSE253564", "GSE285029", "GSE292421", "GSE261345", "GSE261348",
    "GSE309652", "GSE265899", "GSE280232", "GSE300685", "GSE311200",
    "GSE270711", "GSE248249", "GSE241934", "GSE243013", "GSE266219",
    "GSE283829", "GSE225620", "GSE237087", "GSE271689", "GSE292098",
]

RESP_KEYS = ("respon", "recist", "pfs", "overall surviv", " os ", "dcb",
             "benefit", "pcr", "mpr", "pathologic", "progression",
             "sensitiv", "resistan", "outcome", "efficacy", "cr/pr", "pd/sd")


def fetch(gse: str) -> str | None:
    stub = gse[:-3] + "nnn"
    url = (f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/{gse}/matrix/"
           f"{gse}_series_matrix.txt.gz")
    for attempt in range(5):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "geo"})
            with urllib.request.urlopen(req, timeout=90) as r:
                raw = r.read()
            return gzip.GzipFile(fileobj=io.BytesIO(raw)).read().decode(
                "utf-8", "replace")
        except Exception as e:  # noqa: BLE001
            if attempt == 4:
                print(f"  {gse} ERROR {e}")
                return None
            time.sleep(2 ** attempt)
    return None


def main() -> None:
    summary = []
    for gse in SHORTLIST:
        txt = fetch(gse)
        if txt is None:
            summary.append(f"{gse}\tNO_MATRIX (may be multi-platform/superseries)")
            continue
        char_lines = [ln for ln in txt.splitlines()
                      if ln.startswith("!Sample_characteristics")
                      or ln.startswith("!Sample_title")
                      or ln.startswith("!Sample_source_name")
                      or ln.startswith("!Sample_geo_accession")]
        (META / f"{gse}.characteristics.txt").write_text("\n".join(char_lines))
        low = "\n".join(char_lines).lower()
        hits = sorted({k.strip() for k in RESP_KEYS if k in low})
        n_samp = txt.count("!Sample_geo_accession") and \
            len([c for c in char_lines if c.startswith("!Sample_geo_accession")])
        summary.append(f"{gse}\tchar_lines={len(char_lines)}\t"
                       f"resp_hits={hits}")
        print(f"{gse}: {len(char_lines)} char lines, resp_hits={hits}")
        time.sleep(0.3)
    (RES / "label_scan.txt").write_text("\n".join(summary))
    print("wrote", RES / "label_scan.txt")


if __name__ == "__main__":
    main()
