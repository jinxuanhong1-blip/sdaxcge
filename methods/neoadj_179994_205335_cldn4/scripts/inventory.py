#!/usr/bin/env python3
"""Record GEO file sizes and GSE179994 T-cell metadata. Does not download the 421 MiB RDS."""
from __future__ import annotations

import json
import urllib.request
from pathlib import Path

FTP = {
    "GSE179994_all.Tcell.rawCounts.rds.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179994/suppl/GSE179994_all.Tcell.rawCounts.rds.gz",
    "GSE179994_Tcell.metadata.tsv.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179994/suppl/GSE179994_Tcell.metadata.tsv.gz",
    "GSE179994_RAW.tar": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179994/suppl/GSE179994_RAW.tar",
    "GSE205335_Lung_IO_UMI_matrix.rds.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_UMI_matrix.rds.gz",
    "GSE205335_Lung_IO_CellIdentity.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_CellIdentity.txt.gz",
}


def head_length(url: str) -> int:
    req = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(req, timeout=60) as resp:
        return int(resp.headers.get("Content-Length", "0"))


def main() -> None:
    sizes = {name: head_length(url) for name, url in FTP.items()}
    tcell_ok = sizes["GSE179994_all.Tcell.rawCounts.rds.gz"] < 2 * 1024**3
    tme_ok = False  # no TME matrix on GSE179994
    print(json.dumps({"bytes": sizes, "gse179994_tcell_under_2gb": tcell_ok, "gse179994_tme_matrix": tme_ok}, indent=2))


if __name__ == "__main__":
    main()
