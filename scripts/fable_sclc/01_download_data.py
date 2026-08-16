"""Download processed (non-controlled) SCLC datasets used in this slice.

All accessions were manually verified on 2026-08-16 (see
notes/fable_sclc/DATASETS.md). Total processed download is < 2 GB:

  1. George et al. 2015 (Nature; EGA EGAS00001000925 is controlled raw data;
     processed expression is public via cBioPortal study sclc_ucologne_2015).
     Fetched here as the full FPKM matrix through the cBioPortal REST API.
  2. GSE60052 (Jiang et al. 2016): 79 SCLC tumors + 7 normal lungs,
     normalized log2 RNA-seq matrix (~10 MB).
  3. GSE261345 (CANTABRICO trial, durvalumab + platinum/etoposide) GeoMx DSP
     normalized counts + GEO series matrix (clinical annotations).
  4. GSE261348 (IMfirst trial, atezolizumab + platinum/etoposide) GeoMx DSP
     normalized counts + GEO series matrix.
  5. Chan et al. 2021 SCLC single-cell atlas (HTAN MSK / dbGaP phs002371),
     public cellxgene "Combined samples" h5ad (147,137 cells, ~1.46 GB).
     Skipped when --skip-scrna is passed.

Usage: python 01_download_data.py [--skip-scrna]
"""

import argparse
import json
import sys
import urllib.request

import pandas as pd

from common import DATA_DIR

CBIO = "https://www.cbioportal.org/api"

GEO_FILES = {
    "GSE60052_79tumor.7normal.normalized.log2.data.Rda.tsv.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/GSE60052/suppl/"
        "GSE60052_79tumor.7normal.normalized.log2.data.Rda.tsv.gz",
    "GSE261345_CANTABRICO_DSP_normalizedcounts.xlsx":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE261nnn/GSE261345/suppl/"
        "GSE261345_CANTABRICO_DSP_normalizedcounts.xlsx",
    "GSE261348_IMfirst_DSP_normalizedcounts.xlsx":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE261nnn/GSE261348/suppl/"
        "GSE261348_IMfirst_DSP_normalizedcounts.xlsx",
    "GSE261345_series_matrix.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE261nnn/GSE261345/matrix/"
        "GSE261345_series_matrix.txt.gz",
    "GSE261348_series_matrix.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE261nnn/GSE261348/matrix/"
        "GSE261348_series_matrix.txt.gz",
}

# Chan et al. 2021, cellxgene collection 62e8f058-9c37-48bc-9200-e767f318a8ec,
# dataset "Combined samples" (147,137 cells).
CHAN_H5AD_URL = ("https://datasets.cellxgene.cziscience.com/"
                 "a9d92e38-9a6e-401b-8484-74bb15122341.h5ad")


def download(url, dest):
    if dest.exists() and dest.stat().st_size > 0:
        print(f"[skip] {dest.name} already present")
        return
    print(f"[get ] {url}")
    urllib.request.urlretrieve(url, dest)
    print(f"[done] {dest.name}: {dest.stat().st_size / 1e6:.1f} MB")


def cbio_post(path, payload):
    req = urllib.request.Request(
        CBIO + path, data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req) as r:
        return json.load(r)


def fetch_george_full_matrix():
    """Fetch the full George 2015 RNA-seq FPKM matrix via the cBioPortal API.

    Genes are fetched in chunks of 2000 Entrez IDs; only genes with data in
    the sclc_ucologne_2015_rna_seq_mrna profile are returned (~20k genes,
    81 samples; ~30 MB on disk as csv.gz).
    """
    dest = DATA_DIR / "george2015_fpkm_full.csv.gz"
    if dest.exists():
        print(f"[skip] {dest.name} already present")
        return
    with urllib.request.urlopen(
            CBIO + "/genes?pageSize=100000&pageNumber=0&projection=SUMMARY") as r:
        catalog = json.load(r)
    entrez = {g["entrezGeneId"]: g["hugoGeneSymbol"] for g in catalog
              if g["entrezGeneId"] > 0}
    ids = sorted(entrez)
    frames = []
    for i in range(0, len(ids), 2000):
        chunk = ids[i:i + 2000]
        data = cbio_post(
            "/molecular-profiles/sclc_ucologne_2015_rna_seq_mrna/"
            "molecular-data/fetch?projection=SUMMARY",
            {"entrezGeneIds": chunk, "sampleListId": "sclc_ucologne_2015_all"})
        if data:
            frames.append(pd.DataFrame(data)[
                ["entrezGeneId", "sampleId", "value"]])
        print(f"[cbio] genes {i}-{i + len(chunk)}: "
              f"{len(data)} datapoints", flush=True)
    df = pd.concat(frames)
    df["gene"] = df["entrezGeneId"].map(entrez)
    mat = df.pivot_table(index="sampleId", columns="gene", values="value")
    mat.to_csv(dest)
    print(f"[done] {dest.name}: {mat.shape[0]} samples x {mat.shape[1]} genes")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-scrna", action="store_true",
                    help="skip the 1.46 GB Chan 2021 h5ad download")
    args = ap.parse_args()

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    for name, url in GEO_FILES.items():
        download(url, DATA_DIR / name)
    fetch_george_full_matrix()
    if not args.skip_scrna:
        download(CHAN_H5AD_URL, DATA_DIR / "chan2021_combined.h5ad")
    print("All downloads complete.")


if __name__ == "__main__":
    sys.exit(main())
