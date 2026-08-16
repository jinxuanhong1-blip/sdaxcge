#!/usr/bin/env python3
"""Download selected OPEN processed matrix files (<2GB) from verified Zenodo records.

Selection rationale (all CC-BY-4.0, access_right=open, verified DOIs):
  - 10731914: human 8-LUAD scRNA-seq processed count matrix + cell annotation (CSV)
              -> primary analysis target for TACSTD2 / CLDN4 expression.
  - 8041882 : anti-PD1 NSCLC scRNA SingleCellExperiment (RDS) + sample list
              -> most directly on-topic "lung ICI" dataset.

Restricted / >2GB files are intentionally skipped.
Downloads land in results/fable_zenodo/downloads/<record>/ and are checksummed.
"""
import hashlib
import json
import sys
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results" / "fable_zenodo"
DL = RESULTS / "downloads"
UA = "fable-zenodo-download/1.0 (research; contact: cloud-agent)"

MAX_BYTES = 2 * 1024 ** 3

# (record_id, filename) selections
SELECT = [
    ("10731914", "cell_matrix_sce_8LUADs.csv.gz"),
    ("10731914", "cell_annotation_sce_8LUADs.csv"),
    ("8041882", "sce-clustering-allsamples-original_final.RDS"),
    ("8041882", "AK_List.csv"),
]


def download(url, dest: Path):
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    for attempt in range(4):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=120) as r, open(tmp, "wb") as fh:
                total = 0
                h = hashlib.md5()
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
            print(f"   retry {attempt} ({e})", file=sys.stderr)
            time.sleep(2 ** attempt)
    raise RuntimeError(f"failed to download {url}")


def main():
    scored = json.loads((RESULTS / "shortlist.json").read_text())
    by_id = {str(r["id"]): r for r in scored}
    manifest = []
    for rec_id, fname in SELECT:
        rec = by_id.get(rec_id)
        if not rec:
            print(f"! record {rec_id} not in shortlist", file=sys.stderr)
            continue
        fmatch = next((f for f in rec["files"] if f.get("key") == fname), None)
        if not fmatch:
            print(f"! file {fname} not in record {rec_id}", file=sys.stderr)
            continue
        size = fmatch.get("size")
        if size and size > MAX_BYTES:
            print(f"! SKIP {fname} > 2GB", file=sys.stderr)
            continue
        dest = DL / rec_id / fname
        if dest.exists():
            print(f"= exists {dest}")
            got, md5 = dest.stat().st_size, None
        else:
            print(f"> {rec_id}/{fname} ({size} bytes) <- {fmatch['link']}")
            got, md5 = download(fmatch["link"], dest)
            print(f"  done {got} bytes md5={md5}")
        manifest.append({
            "record": rec_id,
            "doi": rec.get("doi"),
            "doi_url": rec.get("doi_url"),
            "license": rec.get("license"),
            "access_right": rec.get("access_right"),
            "file": fname,
            "expected_size": size,
            "downloaded_size": got,
            "md5": md5,
            "path": str(dest.relative_to(ROOT)),
            "source_url": fmatch["link"],
        })
    (RESULTS / "download_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False))
    print(f"\nManifest: {RESULTS / 'download_manifest.json'} ({len(manifest)} files)")


if __name__ == "__main__":
    main()
