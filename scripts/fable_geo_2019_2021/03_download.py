#!/usr/bin/env python3
"""Download OPEN, PROCESSED GEO data + series-matrix metadata for verified series.

Only processed matrices (counts/TPM) and the series_matrix metadata files are
downloaded. Each file is size-checked (HTTP Content-Length) and the hard cap of
2 GB per file is enforced before download. Files are stored under
results/fable_geo_2019_2021/downloads/<GSE>/.
"""
import sys
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DL = ROOT / "results" / "fable_geo_2019_2021" / "downloads"
DL.mkdir(parents=True, exist_ok=True)
MAX_BYTES = 2 * 1024 ** 3  # 2 GB hard cap

FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"

# (GSE, [supplementary processed filenames]) - accessions & filenames are REAL,
# taken verbatim from the GEO suppl/ directory listings.
TARGETS = {
    "GSE126044": ["GSE126044_counts.txt.gz"],
    "GSE135222": ["GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz"],
    "GSE136961": ["GSE136961_TPM.tsv.gz", "GSE136961_raw_count.tsv.gz"],
    "GSE111414": ["GSE111414_gene_counts.csv.gz"],
    "GSE182328": ["GSE182328_Gene_counts_matrix.txt.gz"],
}


def nnn(gse):
    num = gse[3:]
    return f"GSE{num[:-3]}nnn" if len(num) > 3 else "GSEnnn"


def head_size(url):
    req = urllib.request.Request(url, method="HEAD")
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            cl = r.headers.get("Content-Length")
            return int(cl) if cl else None
    except Exception as e:  # noqa: BLE001
        print(f"  HEAD failed: {e}", file=sys.stderr)
        return None


def download(url, dest):
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=300) as r, open(dest, "wb") as f:
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    f.write(chunk)
            return True
        except Exception as e:  # noqa: BLE001
            wait = 2 ** attempt
            print(f"  download retry {attempt} after {wait}s ({e})", file=sys.stderr)
            time.sleep(wait)
    return False


def main():
    for gse, files in TARGETS.items():
        outdir = DL / gse
        outdir.mkdir(exist_ok=True)
        base = f"{FTP}/{nnn(gse)}/{gse}"
        # series matrix (metadata: titles, characteristics -> outcomes)
        files_to_get = [("matrix", f"{gse}_series_matrix.txt.gz")]
        files_to_get += [("suppl", f) for f in files]
        for sub, fname in files_to_get:
            url = f"{base}/{sub}/{fname}"
            dest = outdir / fname
            size = head_size(url)
            human = f"{size/1e6:.2f} MB" if size else "unknown"
            print(f"{gse}/{fname}: size={human}")
            if size and size > MAX_BYTES:
                print(f"  SKIP (>2GB): {fname}")
                continue
            if dest.exists() and dest.stat().st_size > 0:
                print("  already downloaded")
                continue
            ok = download(url, dest)
            print("  ->", "ok" if ok else "FAILED", dest.stat().st_size if dest.exists() else 0, "bytes")
            time.sleep(0.3)


if __name__ == "__main__":
    main()
