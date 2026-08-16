"""Download GSE207422 processed scRNA files from NCBI GEO (<2 GB)."""
from pathlib import Path
import urllib.request

OUT = Path("data/gse207422")
OUT.mkdir(parents=True, exist_ok=True)
BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl"
FILES = {
    "scRNAseq_UMI_matrix.txt.gz": f"{BASE}/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz",
    "scRNAseq_metadata.xlsx": f"{BASE}/GSE207422_NSCLC_scRNAseq_metadata.xlsx",
}

for name, url in FILES.items():
    dest = OUT / name
    if dest.exists() and dest.stat().st_size > 1000:
        print("exists", dest, dest.stat().st_size)
        continue
    print("downloading", url)
    urllib.request.urlretrieve(url, dest)
    print("wrote", dest, dest.stat().st_size)
