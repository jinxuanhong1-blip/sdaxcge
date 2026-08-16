#!/usr/bin/env python3
"""Download open-access TCGA-PAAD matrices for the B1_PAAD co-expression slice.

All files are fully open (no dbGaP / controlled access):
  * UCSC Xena GDC hub (https://gdc.xenahubs.net): STAR TPM expression
    (log2(TPM+1), GENCODE v36), clinical phenotype, gene probemap.
  * ABSOLUTE tumor purity/ploidy from the open PanCanAtlas supplement
    (TCGA_mastercalls.abs_tables_JSedit.fixed.txt) via the GDC API.

Data land in DATA_DIR (default /tmp/b1_paad_data, override with the
B1_PAAD_DATA env var). Raw matrices are intentionally NOT committed.

Usage:  python3 scripts/w200/B1_PAAD/download_data.py
"""

import hashlib
import os
import sys
import urllib.request

DATA_DIR = os.environ.get("B1_PAAD_DATA", "/tmp/b1_paad_data")

GDC_HUB = "https://gdc-hub.s3.us-east-1.amazonaws.com/download"

FILES = {
    # log2(TPM+1), Ensembl gene IDs, GENCODE v36
    "TCGA-PAAD.star_tpm.tsv.gz": f"{GDC_HUB}/TCGA-PAAD.star_tpm.tsv.gz",
    # clinical phenotype (histology, for the adenocarcinoma sensitivity subset)
    "TCGA-PAAD.clinical.tsv.gz": f"{GDC_HUB}/TCGA-PAAD.clinical.tsv.gz",
    # Ensembl -> HGNC symbol mapping
    "gencode.v36.probemap": f"{GDC_HUB}/gencode.v36.annotation.gtf.gene.probemap",
    # ABSOLUTE tumor purity/ploidy, PanCanAtlas open supplement
    "tcga_absolute_purity.txt": (
        "https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5"
    ),
}


def fetch(name: str, url: str) -> None:
    dest = os.path.join(DATA_DIR, name)
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        print(f"[skip] {name} already present")
        return
    print(f"[get ] {url}")
    tmp = dest + ".part"
    urllib.request.urlretrieve(url, tmp)
    os.replace(tmp, dest)
    md5 = hashlib.md5(open(dest, "rb").read()).hexdigest()
    print(f"[done] {name}  size={os.path.getsize(dest):,}  md5={md5}")


def main() -> int:
    os.makedirs(DATA_DIR, exist_ok=True)
    for name, url in FILES.items():
        fetch(name, url)
    print(f"\nAll files in {DATA_DIR}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
