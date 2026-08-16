#!/usr/bin/env python3
"""Download processed matrices only for viable mouse lung ICI datasets.

Rules enforced by the task:
  - processed matrices only (no FASTQ, no raw MS .raw/.mgf/.msf)
  - skip files > 2 GB
  - do not invent accessions (URLs come from real GEO/AE/PRIDE listings)
"""
import os
import sys
import time
import urllib.request

DATA = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "notes", "mouse", "data"))
os.makedirs(DATA, exist_ok=True)

MAX_BYTES = 2 * 1024 ** 3  # 2 GB

UA = {"User-Agent": "Mozilla/5.0 (mouse-ici-catalog)"}

# (accession, url, dest filename)
DOWNLOADS = [
    # --- bulk / processed tables ---
    ("E-MTAB-13704", "https://www.ebi.ac.uk/biostudies/files/E-MTAB-13704/GEMMS_raw_counts.csv", "E-MTAB-13704_GEMMS_raw_counts.csv"),
    ("E-MTAB-13704", "https://www.ebi.ac.uk/biostudies/files/E-MTAB-13704/E-MTAB-13704.sdrf.txt", "E-MTAB-13704.sdrf.txt"),
    ("GSE197260", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE197nnn/GSE197260/suppl/GSE197260_RNAseqTPM_MM_EGFR-TKI-CD8.txt.gz", "GSE197260_RNAseqTPM.txt.gz"),
    ("GSE239485", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE239nnn/GSE239485/suppl/GSE239485_Processed_data.xlsx", "GSE239485_Processed_data.xlsx"),
    ("GSE297630", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE297nnn/GSE297630/suppl/GSE297630_processed_data.xlsx", "GSE297630_processed_data.xlsx"),
    ("GSE241978", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241978/suppl/GSE241978_2020-07-21_Sherr_analysis_CMT_KO_vs_Cas9Ctrl.xlsx", "GSE241978_CMT_KO_vs_Cas9Ctrl.xlsx"),
    # bulk RNA-seq per-sample processed tables inside RAW.tar (xlsx, not FASTQ)
    ("GSE330658", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE330nnn/GSE330658/suppl/GSE330658_RAW.tar", "GSE330658_RAW.tar"),
    # --- scRNA-seq processed matrices (10x mtx) ---
    ("GSE129297", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE129nnn/GSE129297/suppl/GSE129297_barcodes.tsv.gz", "GSE129297_barcodes.tsv.gz"),
    ("GSE129297", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE129nnn/GSE129297/suppl/GSE129297_features.tsv.gz", "GSE129297_features.tsv.gz"),
    ("GSE129297", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE129nnn/GSE129297/suppl/GSE129297_RAW.tar", "GSE129297_RAW.tar"),
    ("GSE133604", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE133nnn/GSE133604/suppl/GSE133604_barcodes.tsv.gz", "GSE133604_barcodes.tsv.gz"),
    ("GSE133604", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE133nnn/GSE133604/suppl/GSE133604_genes.tsv.gz", "GSE133604_genes.tsv.gz"),
    ("GSE133604", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE133nnn/GSE133604/suppl/GSE133604_RAW.tar", "GSE133604_RAW.tar"),
    ("GSE222158", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE222nnn/GSE222158/suppl/GSE222158_RAW.tar", "GSE222158_RAW.tar"),
    ("GSE297632", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE297nnn/GSE297632/suppl/GSE297632_RAW.tar", "GSE297632_RAW.tar"),
    # --- PRIDE processed protein result (only sub-2GB table available) ---
    ("PXD059688", "https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/01/PXD059688/20230904_Tumor_AA.mzTab.gz", "PXD059688_Tumor_AA.mzTab.gz"),
]


def head_size(url):
    try:
        req = urllib.request.Request(url, headers=UA, method="HEAD")
        with urllib.request.urlopen(req, timeout=60) as r:
            cl = r.headers.get("Content-Length")
            return int(cl) if cl else None
    except Exception:
        return None


def download(url, dest):
    tmp = dest + ".part"
    req = urllib.request.Request(url, headers=UA)
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=300) as r, open(tmp, "wb") as f:
                total = 0
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    total += len(chunk)
                    if total > MAX_BYTES:
                        f.close()
                        os.remove(tmp)
                        return None, "exceeds 2GB during stream"
                    f.write(chunk)
            os.replace(tmp, dest)
            return total, None
        except Exception as e:  # noqa: BLE001
            sys.stderr.write(f"  attempt {attempt+1} failed: {e}\n")
            time.sleep(2 ** attempt)
    return None, "download failed"


def main():
    only = sys.argv[1:] if len(sys.argv) > 1 else None
    log = []
    for acc, url, name in DOWNLOADS:
        if only and acc not in only and name not in only:
            continue
        dest = os.path.join(DATA, name)
        if os.path.exists(dest) and os.path.getsize(dest) > 0:
            print(f"[skip exists] {name} ({os.path.getsize(dest)} B)")
            log.append((acc, name, os.path.getsize(dest), "exists"))
            continue
        sz = head_size(url)
        if sz is not None and sz > MAX_BYTES:
            print(f"[SKIP >2GB] {name} ({sz} B)  {url}")
            log.append((acc, name, sz, "skipped >2GB"))
            continue
        print(f"[get] {name}  (reported {sz} B)")
        got, err = download(url, dest)
        if err:
            print(f"   ERROR {name}: {err}")
            log.append((acc, name, sz, f"error: {err}"))
        else:
            print(f"   ok {got} B")
            log.append((acc, name, got, "ok"))
    print("\n== summary ==")
    for row in log:
        print(row)


if __name__ == "__main__":
    main()
