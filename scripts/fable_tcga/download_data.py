#!/usr/bin/env python3
"""Download open-access TCGA LUAD / LUSC matrices used by the fable_tcga slice.

All files are fully open (no dbGaP / controlled access):
  * UCSC Xena GDC hub  (https://gdc.xenahubs.net): STAR TPM expression,
    OS survival, clinical phenotype, GENCODE v36 gene probemap.
  * UCSC Xena Pan-Cancer Atlas hub: TCGA-CDR curated survival endpoints
    (Liu et al., Cell 2018) providing OS / DSS / DFI / PFI.
  * MCP-counter marker-gene signatures (Becht et al., Genome Biology 2016)
    from the authors' public GitHub repository.

Data land in DATA_DIR (default /tmp/fable_tcga_data, override with the
FABLE_TCGA_DATA env var). Raw matrices are intentionally NOT committed to
the repository.

Usage:  python3 scripts/fable_tcga/download_data.py
"""

import hashlib
import os
import sys
import urllib.request

DATA_DIR = os.environ.get("FABLE_TCGA_DATA", "/tmp/fable_tcga_data")

GDC_HUB = "https://gdc-hub.s3.us-east-1.amazonaws.com/download"
PANCAN_HUB = "https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download"

FILES = {
    # cohort-level matrices (log2(TPM+1); Ensembl gene IDs)
    "TCGA-LUAD.star_tpm.tsv.gz": f"{GDC_HUB}/TCGA-LUAD.star_tpm.tsv.gz",
    "TCGA-LUSC.star_tpm.tsv.gz": f"{GDC_HUB}/TCGA-LUSC.star_tpm.tsv.gz",
    # OS survival + clinical covariates from the GDC hub
    "TCGA-LUAD.survival.tsv.gz": f"{GDC_HUB}/TCGA-LUAD.survival.tsv.gz",
    "TCGA-LUSC.survival.tsv.gz": f"{GDC_HUB}/TCGA-LUSC.survival.tsv.gz",
    "TCGA-LUAD.clinical.tsv.gz": f"{GDC_HUB}/TCGA-LUAD.clinical.tsv.gz",
    "TCGA-LUSC.clinical.tsv.gz": f"{GDC_HUB}/TCGA-LUSC.clinical.tsv.gz",
    # Ensembl -> HGNC symbol mapping
    "gencode.v36.probemap": f"{GDC_HUB}/gencode.v36.annotation.gtf.gene.probemap",
    # TCGA-CDR curated endpoints (OS/DSS/DFI/PFI), pan-cancer
    "tcga_cdr_survival.tsv": f"{PANCAN_HUB}/Survival_SupplementalTable_S1_20171025_xena_sp",
    # MCP-counter marker genes
    "mcpcounter_genes.txt": (
        "https://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt"
    ),
    # ABSOLUTE tumor purity/ploidy, PanCanAtlas supplement
    # (TCGA_mastercalls.abs_tables_JSedit.fixed.txt; open-access GDC file UUID)
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
