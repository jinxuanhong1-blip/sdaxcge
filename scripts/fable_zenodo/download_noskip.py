#!/usr/bin/env python3
"""Download open processed matrices that can contain TACSTD2/CLDN4. No size cap.

Restricted records are still skipped. Files land in results/noskip/zenodo/downloads/.
"""
import hashlib
import json
import sys
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "noskip" / "zenodo"
DL = OUT / "downloads"
UA = "fable-noskip/1.0 (research)"

# (record_id, filename) — open, processed or archives of processed matrices
SELECT = [
    # Nanostring nCounter, advanced NSCLC anti-PD1 (bulk-ish gene panel)
    ("2635194", "learn-data.csv"),
    ("2635194", "learn-metadata.csv"),
    ("2635194", "validate-data.csv"),
    ("2635194", "validate-metadata.csv"),
    # Early LUAD Visium 10x h5 + processed malignant/immune tables
    ("8417887", "P01_filtered_feature_bc_matrix.h5"),
    ("8417887", "P02_filtered_feature_bc_matrix.h5"),
    ("8417887", "P03_filtered_feature_bc_matrix.h5"),
    ("8417887", "P04_filtered_feature_bc_matrix.h5"),
    ("8417887", "P11_filtered_feature_bc_matrix.h5"),
    ("8417887", "P13_filtered_feature_bc_matrix.h5"),
    ("8417887", "P14_filtered_feature_bc_matrix.h5"),
    ("8417887", "Malignant_cell_data.txt"),
    ("8417887", "Malignant_cell_meta_data.txt"),
    ("8417887", "Immune_cell_meta_data.txt"),
    ("8417887", "Immune_cell_data.txt"),  # 2.67 GB processed table — do not skip
    # Paired normal-LUAD 10x h5 archive
    ("11205626", "24samples.h5.tar.gz"),
    ("11205626", "24samples_metadata.csv"),
    # Mouse LUAD companion matrix (Tacstd2 / Cldn4)
    ("10731914", "mouse_scRNAseq_adata_raw_counts.csv.gz"),
    ("10731914", "mouse_scRNAseq_adata_metadata.csv"),
    # I3LUNG immunotherapy clinical tables (small zip first)
    ("17535424", "I3LUNG_DATA.zip"),
    # Processed LUAD scRNA (1.5 GB zip)
    ("13947395", "Fig6_processed_scRNAseq.zip"),
]


def zenodo_file_url(rec_id, fname):
    return f"https://zenodo.org/api/records/{rec_id}/files/{urllib.request.quote(fname)}/content"


def download(url, dest: Path):
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    last = None
    for attempt in range(5):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=180) as r, open(tmp, "wb") as fh:
                h = hashlib.md5()
                total = 0
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    fh.write(chunk)
                    h.update(chunk)
                    total += len(chunk)
            tmp.rename(dest)
            return total, h.hexdigest()
        except Exception as e:  # noqa
            last = e
            print(f"   retry {attempt} ({e})", file=sys.stderr)
            time.sleep(2 ** attempt)
    raise RuntimeError(f"failed {url}: {last}")


def main():
    import urllib.parse
    urllib.request.quote = urllib.parse.quote
    manifest = []
    for rec_id, fname in SELECT:
        dest = DL / rec_id / fname
        url = zenodo_file_url(rec_id, fname)
        if dest.exists() and dest.stat().st_size > 0:
            print(f"= exists {dest} ({dest.stat().st_size})")
            got, md5 = dest.stat().st_size, None
        else:
            print(f"> {rec_id}/{fname}")
            got, md5 = download(url, dest)
            print(f"  done {got} bytes md5={md5}")
        manifest.append({
            "record": rec_id,
            "file": fname,
            "path": str(dest.relative_to(ROOT)),
            "downloaded_size": got,
            "md5": md5,
            "source_url": url,
        })
    (OUT / "download_manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"Wrote {OUT / 'download_manifest.json'} ({len(manifest)} files)")


if __name__ == "__main__":
    main()
