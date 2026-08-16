#!/usr/bin/env python3
"""Download open processed supplementary files (<2GB each) for the datasets
selected for TACSTD2/CLDN4 vs response analysis. All URLs are the live GEO
FTP paths recorded in supp_listing.tsv (no invented accessions).
Files land in results/fable_geo_2026/data/<GSE>/.
"""
import time
import urllib.request
from pathlib import Path

RES = Path(__file__).resolve().parents[2] / "results" / "fable_geo_2026"
DATA = RES / "data"

FILES = {
    "GSE261345": [
        "GSE261345_CANTABRICO_DSP_normalizedcounts.xlsx",
        "GSE261345_CANTABRICO_DSP_rawcounts.xlsx",
    ],
    "GSE261348": [
        "GSE261348_IMfirst_DSP_normalizedcounts.xlsx",
        "GSE261348_IMfirst_DSP_rawcounts.xlsx",
    ],
    "GSE309652": ["GSE309652_RAW.tar", "filelist.txt"],
    "GSE253564": ["GSE253564_Pre-treatment_Samples_Pubs_FPKMs.txt.gz"],
    "GSE292421": ["GSE292421_FPKM.csv.gz", "GSE292421_COUNTS.csv.gz"],
}


def url_for(gse: str, fname: str) -> str:
    stub = gse[:-3] + "nnn"
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/{gse}/suppl/{fname}"


def download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(5):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "geo"})
            with urllib.request.urlopen(req, timeout=300) as r, \
                    dest.open("wb") as fh:
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    fh.write(chunk)
            print(f"  OK {dest.name} ({dest.stat().st_size/1e6:.1f} MB)")
            return
        except Exception as e:  # noqa: BLE001
            if attempt == 4:
                print(f"  FAIL {url}: {e}")
                return
            time.sleep(2 ** attempt)


def main() -> None:
    for gse, files in FILES.items():
        print(gse)
        for f in files:
            download(url_for(gse, f), DATA / gse / f)


if __name__ == "__main__":
    main()
