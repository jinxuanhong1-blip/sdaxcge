#!/usr/bin/env python3
"""Record GEO file sizes for the 229353+205335 pair gate. Does not download the UMI RDS."""
from __future__ import annotations

import json
import urllib.request

FTP = {
    "GSE229353_RAW.tar": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE229nnn/GSE229353/suppl/GSE229353_RAW.tar",
    "GSE205335_Lung_IO_UMI_matrix.rds.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_UMI_matrix.rds.gz",
    "GSE205335_Lung_IO_CellIdentity.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_CellIdentity.txt.gz",
}


def head_length(url: str) -> int:
    req = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(req, timeout=60) as resp:
        return int(resp.headers.get("Content-Length", "0"))


def main() -> None:
    sizes = {name: head_length(url) for name, url in FTP.items()}
    print(
        json.dumps(
            {
                "bytes": sizes,
                "gse229353_tar_under_2gb": sizes["GSE229353_RAW.tar"] < 2 * 1024**3,
                "gse229353_malignant_matrix": False,
                "gse205335_umi_under_2gb": sizes["GSE205335_Lung_IO_UMI_matrix.rds.gz"]
                < 2 * 1024**3,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
