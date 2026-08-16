#!/usr/bin/env python3
"""Download leftover 2019-2021 processed matrices + series_matrix metadata.

No skip for tissue (blood/PBMC) or size under 2 GB. Files larger than 2 GB are
still attempted if they are processed expression matrices (logged, not silently
dropped). RAW.tar-only series are recorded as 'no processed matrix deposited'.
"""
import json
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "noskip" / "GEO_2019_2021"
DL = OUT / "downloads"
DL.mkdir(parents=True, exist_ok=True)
FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"
MAX_SOFT = 2 * 1024 ** 3

# Leftover processed expression-like files (real GEO filenames from inventory).
# Includes blood, large, cell-line, and non-outcome series. We do not skip them
# here; gene presence + outcome computability are decided after download.
TARGETS = {
    "GSE181820": ["GSE181820_RNAseq_NSCLC.txt.gz"],
    "GSE190265": ["GSE190265_TPM_France3.csv.gz"],
    "GSE190266": ["GSE190266_TPM_France4.csv.gz"],
    "GSE152590": ["GSE152590_TPM_matrix.xlsx"],
    "GSE145896": ["GSE145896_countData.csv.gz"],
    "GSE146100": ["GSE146100_NormData.txt.gz"],
    "GSE131933": [
        "GSE131933_T1T3D0_gene_count.txt.gz",
        "GSE131933_T3D11_gene_count.txt.gz",
        "GSE131933_T1D7_gene_count.txt.gz",
    ],
    "GSE184053": [
        "GSE184053_raw_gene_counts_1.txt.gz",
        "GSE184053_raw_gene_counts_2.txt.gz",
        "GSE184053_raw_gene_counts_3.txt.gz",
    ],
    "GSE128822": ["GSE128822_ICOSpCCRp_ICOSmCCR8m_raw_counts.csv.gz"],
    "GSE179994": ["GSE179994_all.Tcell.rawCounts.rds.gz"],
    "GSE179934": ["GSE179934_Raw_gene_counts_matrix.txt.gz"],
    "GSE129381": ["GSE129381_counts_filtered_anno.txt.gz", "GSE129381_quantseq_counts.txt.gz"],
    "GSE129968": ["GSE129968_PDL1.totalnorm.txt.gz"],
    "GSE171650": ["GSE171650_PC9ER_cells-_PD-L1_expression_in_wt_vs_ox_lung_cancer_cells_DEoutput_All.xlsx"],
    "GSE118933": ["GSE118933_IPF_invasive_non-invasive_RPKM.txt.gz"],
    "GSE173896": ["GSE173896_Arzamasov_raw_count_matrix.txt.gz"],
    "GSE162154": ["GSE162154_processed_data.xls.gz"],
    "GSE168707": ["GSE168707_non-normalized_data.txt.gz"],
    "GSE159787": ["GSE159787_analyzed_counts.txt.gz"],
    "GSE159785": ["GSE159785_analyzed_count_matrix.txt.gz"],
    "GSE139327": ["GSE139327_bulk_cells_gene_exp.diff.gz"],
    "GSE124885": ["GSE124885_T_5prime_raw_gene_expression.csv.gz"],
}

# Also pull series_matrix for leftover ICI-relevant / leftover-compute series
# (including RAW-only ones, so we can document outcomes honestly).
MATRIX_ONLY = [
    "GSE141479", "GSE180347", "GSE173351", "GSE154286", "GSE144945",
    "GSE120028", "GSE99995", "GSE117570", "GSE176021", "GSE176022",
    "GSE181820", "GSE190265", "GSE190266", "GSE152590", "GSE145896",
    "GSE146100", "GSE131933", "GSE184053", "GSE128822", "GSE179994",
    "GSE179934", "GSE129381", "GSE129968", "GSE171650",
]


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
        print(f"  HEAD fail {url}: {e}")
        return None


def download(url, dest):
    dest.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=600) as r, open(dest, "wb") as f:
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    f.write(chunk)
            return True
        except Exception as e:  # noqa: BLE001
            wait = 2 ** attempt
            print(f"  retry {attempt} {wait}s ({e})")
            time.sleep(wait)
    return False


def get_one(gse, sub, fname, log):
    url = f"{FTP}/{nnn(gse)}/{gse}/{sub}/{fname}"
    dest = DL / gse / fname
    size = head_size(url)
    human = f"{size/1e6:.2f} MB" if size else "unknown"
    rec = {"gse": gse, "file": fname, "url": url, "size_bytes": size, "status": ""}
    print(f"{gse}/{fname}: {human}")
    if dest.exists() and dest.stat().st_size > 0:
        rec["status"] = "already"
        rec["local_bytes"] = dest.stat().st_size
        print("  already")
        log.append(rec)
        return
    if size and size > MAX_SOFT:
        print(f"  NOTE >2GB ({human}); still downloading (no-skip-for-size)")
    ok = download(url, dest)
    rec["status"] = "ok" if ok else "FAILED"
    rec["local_bytes"] = dest.stat().st_size if dest.exists() else 0
    print("  ->", rec["status"], rec["local_bytes"])
    log.append(rec)
    time.sleep(0.2)


def main():
    log = []
    seen_matrix = set()
    for gse, files in TARGETS.items():
        get_one(gse, "matrix", f"{gse}_series_matrix.txt.gz", log)
        seen_matrix.add(gse)
        for fname in files:
            get_one(gse, "suppl", fname, log)
    for gse in MATRIX_ONLY:
        if gse in seen_matrix:
            continue
        get_one(gse, "matrix", f"{gse}_series_matrix.txt.gz", log)
    (OUT / "download_log.json").write_text(json.dumps(log, indent=2))
    print("Wrote", OUT / "download_log.json")


if __name__ == "__main__":
    main()
