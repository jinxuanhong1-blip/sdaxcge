#!/usr/bin/env python3
"""Download public GEO processed files used by 03_score.py (idempotent)."""
from __future__ import annotations

import subprocess
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "notes/mouse_scrna_ici_pool/raw/data"
TMP = Path("/tmp/geo_dl")

SERIES_TAR = {
    "GSE157881": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157881/suppl/GSE157881_RAW.tar",
    "GSE157882": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157882/suppl/GSE157882_RAW.tar",
    "GSE268525": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE268nnn/GSE268525/suppl/GSE268525_RAW.tar",
    "GSE283827": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE283nnn/GSE283827/suppl/GSE283827_RAW.tar",
    "GSE303943": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE303nnn/GSE303943/suppl/GSE303943_RAW.tar",
    "GSE133604": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE133nnn/GSE133604/suppl/GSE133604_RAW.tar",
    "GSE129297": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE129nnn/GSE129297/suppl/GSE129297_RAW.tar",
    "GSE232730": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE232nnn/GSE232730/suppl/GSE232730_RAW.tar",
    "GSE222158": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE222nnn/GSE222158/suppl/GSE222158_RAW.tar",
}

SAMPLE_FTP = [
    # GSE176091 CD45 control vs CA170 only (full series RAW is 794 MB)
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5354nnn/GSM5354872/suppl/GSM5354872_CD45_Control_CD45_1_barcodes.tsv.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5354nnn/GSM5354872/suppl/GSM5354872_CD45_Control_CD45_1_features.tsv.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5354nnn/GSM5354872/suppl/GSM5354872_CD45_Control_CD45_1_matrix.mtx.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5354nnn/GSM5354873/suppl/GSM5354873_CD45_Control_CD45_2_barcodes.tsv.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5354nnn/GSM5354873/suppl/GSM5354873_CD45_Control_CD45_2_features.tsv.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5354nnn/GSM5354873/suppl/GSM5354873_CD45_Control_CD45_2_matrix.mtx.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5354nnn/GSM5354876/suppl/GSM5354876_CD45_CA170_1_barcodes.tsv.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5354nnn/GSM5354876/suppl/GSM5354876_CD45_CA170_1_features.tsv.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5354nnn/GSM5354876/suppl/GSM5354876_CD45_CA170_1_matrix.mtx.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5354nnn/GSM5354877/suppl/GSM5354877_CD45_CA170_2_barcodes.tsv.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5354nnn/GSM5354877/suppl/GSM5354877_CD45_CA170_2_features.tsv.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5354nnn/GSM5354877/suppl/GSM5354877_CD45_CA170_2_matrix.mtx.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE275nnn/GSE275877/suppl/GSE275877_features.tsv.gz",
]


def curl(url: str, dest: Path):
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        print("have", dest.name)
        return
    print("GET", dest.name)
    subprocess.check_call(["curl", "-fL", "--retry", "4", "--retry-delay", "4", "-o", str(dest), url])


def main():
    DATA.mkdir(parents=True, exist_ok=True)
    TMP.mkdir(parents=True, exist_ok=True)
    for acc, url in SERIES_TAR.items():
        tar_path = TMP / f"{acc}_RAW.tar"
        curl(url, tar_path)
        print("extract", acc)
        with tarfile.open(tar_path) as tf:
            tf.extractall(DATA)
    for url in SAMPLE_FTP:
        curl(url, DATA / url.rsplit("/", 1)[-1])
    print("done", DATA)


if __name__ == "__main__":
    main()
