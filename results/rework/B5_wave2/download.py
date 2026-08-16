#!/usr/bin/env python3
"""Download public processed ICI matrices for B5 wave 2.

Lung GEO first, then Zenodo BHK-lab ICB TSVs (CC-BY-4.0, 10.5281/zenodo.7058399).
No FASTQ. No EGA. Liu/Braun raw are controlled; this script uses the public
Zenodo processed expression+response tables only.

Raw archives stay under --raw-dir (default /tmp/b5_wave2_raw).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path
from urllib.request import Request, urlopen

HERE = Path(__file__).resolve().parent
FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"
ZENODO = "https://zenodo.org/api/records/7058399/files"

GEO_FILES = {
    "GSE126044_series_matrix.txt.gz": f"{FTP}/GSE126nnn/GSE126044/matrix/GSE126044_series_matrix.txt.gz",
    "GSE126044_counts.txt.gz": f"{FTP}/GSE126nnn/GSE126044/suppl/GSE126044_counts.txt.gz",
    "GSE135222_series_matrix.txt.gz": f"{FTP}/GSE135nnn/GSE135222/matrix/GSE135222_series_matrix.txt.gz",
    "GSE135222_exp.tsv.gz": f"{FTP}/GSE135nnn/GSE135222/suppl/GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz",
    "GSE166449_series_matrix.txt.gz": f"{FTP}/GSE166nnn/GSE166449/matrix/GSE166449_series_matrix.txt.gz",
    "GSE166449_TPM.txt.gz": f"{FTP}/GSE166nnn/GSE166449/suppl/GSE166449_Raw_gene_TPM_matrix.txt.gz",
    "GSE190265_series_matrix.txt.gz": f"{FTP}/GSE190nnn/GSE190265/matrix/GSE190265_series_matrix.txt.gz",
    "GSE190265_TPM_France3.csv.gz": f"{FTP}/GSE190nnn/GSE190265/suppl/GSE190265_TPM_France3.csv.gz",
    "GSE190265_samples_info_France3.csv.gz": f"{FTP}/GSE190nnn/GSE190265/suppl/GSE190265_samples_info_France3.csv.gz",
    "GSE190266_series_matrix.txt.gz": f"{FTP}/GSE190nnn/GSE190266/matrix/GSE190266_series_matrix.txt.gz",
    "GSE190266_TPM_France4.csv.gz": f"{FTP}/GSE190nnn/GSE190266/suppl/GSE190266_TPM_France4.csv.gz",
    "GSE207422_series_matrix.txt.gz": f"{FTP}/GSE207nnn/GSE207422/matrix/GSE207422_series_matrix.txt.gz",
    "GSE207422_log2TPM.txt.gz": f"{FTP}/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz",
    "GSE207422_metadata.xlsx": f"{FTP}/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_bulk_RNAseq_metadata.xlsx",
    "GSE283829_series_matrix.txt.gz": f"{FTP}/GSE283nnn/GSE283829/matrix/GSE283829_series_matrix.txt.gz",
    "GSE283829_raw_express_matrix_all_samples.txt.gz": f"{FTP}/GSE283nnn/GSE283829/suppl/GSE283829_raw_express_matrix_all_samples.txt.gz",
    "GSE253564_series_matrix.txt.gz": f"{FTP}/GSE253nnn/GSE253564/matrix/GSE253564_series_matrix.txt.gz",
    "GSE253564_Pre-treatment_Samples_Pubs_FPKMs.txt.gz": f"{FTP}/GSE253nnn/GSE253564/suppl/GSE253564_Pre-treatment_Samples_Pubs_FPKMs.txt.gz",
}

# Open processed ICB objects. Hwang/Jerby_Arnon/Roh are panels (expected gene-absent).
ZENODO_FILES = {
    name: f"{ZENODO}/{name}/content"
    for name in [
        "ICB_Mariathasan.zip",
        "ICB_Gide.zip",
        "ICB_Liu.zip",
        "ICB_Riaz.zip",
        "ICB_Braun.zip",
        "ICB_Jung.zip",
        "ICB_Hugo.zip",
        "ICB_Snyder.zip",
        "ICB_Shiuan.zip",
        "ICB_Puch.zip",
        "ICB_Miao1.zip",
        "ICB_Kim.zip",
        "ICB_Van_Allen.zip",
        "ICB_Nathanson.zip",
        "ICB_Hwang.zip",
        "ICB_Jerby_Arnon.zip",
        "ICB_Roh.zip",
    ]
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, dest: Path, retries: int = 5) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        print(f"exists {dest.name} ({dest.stat().st_size} bytes)")
        return
    tmp = dest.with_suffix(dest.suffix + ".part")
    last_err = None
    for i in range(retries):
        try:
            req = Request(url, headers={"User-Agent": "B5_wave2/1.0"})
            with urlopen(req, timeout=300) as r, open(tmp, "wb") as out:
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    out.write(chunk)
            if tmp.stat().st_size == 0:
                raise RuntimeError("empty download")
            tmp.replace(dest)
            print(f"OK {dest.name} ({dest.stat().st_size} bytes)")
            return
        except Exception as e:
            last_err = e
            print(f"retry {i + 1} {dest.name}: {e}")
            if tmp.exists():
                tmp.unlink()
            time.sleep(min(60, 2 ** (i + 2)))
    raise RuntimeError(f"failed {url}: {last_err}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-dir", default="/tmp/b5_wave2_raw")
    args = ap.parse_args()
    raw = Path(args.raw_dir)
    raw.mkdir(parents=True, exist_ok=True)

    manifest = []
    files = {**GEO_FILES, **ZENODO_FILES}
    print(f"downloading {len(files)} files to {raw}")
    for name, url in files.items():
        dest = raw / name
        download(url, dest)
        manifest.append(
            {
                "file": name,
                "bytes": dest.stat().st_size,
                "sha256": sha256_file(dest),
                "url": url,
            }
        )
    out = HERE / "tables"
    out.mkdir(parents=True, exist_ok=True)
    (out / "download_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"wrote {out / 'download_manifest.json'} n={len(manifest)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
