#!/usr/bin/env python3
"""Download public DepMap portal files for lung TROP2 vs CLDN1/4/7.

The portal HTML is Cloudflare-gated. The no-captcha catalog is the public
file list:

  https://depmap.org/portal/api/no-captcha/download/files

DepMap Public 26Q1 is the newest release in that catalog, but its expression
URLs are blank (captcha-only). The newest release that still publishes direct
URLs is DepMap Public 24Q4. This script uses those URLs.

CCLE proteomics is the catalog release named "Proteomics"
(protein_quant_current_normalized.csv), whose portal URL is the Gygi lab
page. The matrix downloaded here is the gzip the page serves.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import sys
import urllib.request
from pathlib import Path

UA = "sdaxcge-depmap-lung-trop2-cldn147/1.0 (+https://github.com/jinxuanhong1-blip/sdaxcge)"
CATALOG_URL = "https://depmap.org/portal/api/no-captcha/download/files"
GYGI_GZ_URL = "https://gygi.hms.harvard.edu/data/ccle/protein_quant_current_normalized.csv.gz"
RNA_RELEASE = "DepMap Public 24Q4"
RNA_FILES = {
    "Model.csv": "Model.csv",
    "OmicsExpressionProteinCodingGenesTPMLogp1.csv": "OmicsExpressionProteinCodingGenesTPMLogp1.csv",
}
GENES = {
    "TACSTD2": "TACSTD2 (4070)",
    "CLDN1": "CLDN1 (9076)",
    "CLDN3": "CLDN3 (1365)",
    "CLDN4": "CLDN4 (1364)",
    "CLDN5": "CLDN5 (7122)",
    "CLDN7": "CLDN7 (1366)",
    "CLDN18": "CLDN18 (51208)",
    "EPCAM": "EPCAM (4072)",
    "KRT8": "KRT8 (3856)",
    "KRT18": "KRT18 (3875)",
    "KRT19": "KRT19 (3880)",
}
NEWER_RELEASES = ("DepMap Public 25Q2", "DepMap Public 25Q3", "DepMap Public 26Q1")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(url: str, dest: Path | None = None) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    print(f"GET {url}", flush=True)
    with urllib.request.urlopen(req, timeout=600) as resp:
        if dest is None:
            return resp.read()
        dest.parent.mkdir(parents=True, exist_ok=True)
        h = hashlib.sha256()
        n = 0
        with dest.open("wb") as out:
            while True:
                chunk = resp.read(1024 * 1024)
                if not chunk:
                    break
                out.write(chunk)
                h.update(chunk)
                n += len(chunk)
        print(f"  wrote {dest} ({n} bytes)", flush=True)
        return h.hexdigest().encode()


def download(url: str, dest: Path) -> str:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 1000:
        print(f"cache hit {dest} ({dest.stat().st_size} bytes)", flush=True)
        return sha256_file(dest)
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    print(f"GET {url} -> {dest}", flush=True)
    h = hashlib.sha256()
    with urllib.request.urlopen(req, timeout=900) as resp, dest.open("wb") as out:
        while True:
            chunk = resp.read(1024 * 1024)
            if not chunk:
                break
            out.write(chunk)
            h.update(chunk)
    return h.hexdigest()


def load_catalog() -> list[dict]:
    raw = fetch(CATALOG_URL)
    return list(csv.DictReader(io.StringIO(raw.decode("utf-8"))))


def catalog_row(rows: list[dict], release: str, filename: str) -> dict:
    hits = [r for r in rows if r.get("release") == release and r.get("filename") == filename]
    if len(hits) != 1:
        raise SystemExit(f"Expected one catalog row for {release} / {filename}, found {len(hits)}")
    return hits[0]


def extract_rna(expr_path: Path, out_path: Path) -> dict:
    print(f"Extracting selected genes from {expr_path}", flush=True)
    with expr_path.open("r", newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        missing = [col for col in GENES.values() if col not in header]
        if missing:
            raise SystemExit(f"Expression header missing {missing}")
        idx = {sym: header.index(col) for sym, col in GENES.items()}
        n = 0
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with out_path.open("w", newline="") as out:
            w = csv.writer(out)
            w.writerow(["ModelID", *GENES.keys()])
            for row in reader:
                if not row:
                    continue
                w.writerow([row[0], *[row[idx[sym]] for sym in GENES]])
                n += 1
    return {"n_models": n, "columns": GENES, "id_header": header[0]}


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--cache-dir", default="/tmp/depmap_trop2_cldn147")
    p.add_argument("--out-dir", default="results/depmap_lung_trop2_cldn147")
    args = p.parse_args()
    cache = Path(args.cache_dir)
    out = Path(args.out_dir)
    cache.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)

    rows = load_catalog()
    used = []
    files_meta = {}

    for filename in RNA_FILES:
        row = catalog_row(rows, RNA_RELEASE, filename)
        url = (row.get("url") or "").strip()
        if not url.startswith("http"):
            raise SystemExit(f"{RNA_RELEASE} {filename} has no public URL in the portal catalog")
        dest = cache / filename
        digest = download(url, dest)
        used.append(row)
        files_meta[filename] = {
            "release": RNA_RELEASE,
            "url": url,
            "catalog_md5": row.get("md5_hash") or "",
            "sha256": digest,
            "bytes": dest.stat().st_size,
        }

    prot_row = catalog_row(rows, "Proteomics", "protein_quant_current_normalized.csv")
    used.append(prot_row)
    page_url = (prot_row.get("url") or "").strip()
    page_html = fetch(page_url).decode("utf-8", errors="replace") if page_url.startswith("http") else ""
    if "protein_quant_current_normalized.csv" not in page_html:
        raise SystemExit(
            "Portal Proteomics URL did not mention protein_quant_current_normalized.csv: " + page_url
        )
    gz = cache / "protein_quant_current_normalized.csv.gz"
    gz_sha = download(GYGI_GZ_URL, gz)
    files_meta["protein_quant_current_normalized.csv.gz"] = {
        "release": "Proteomics",
        "portal_filename": "protein_quant_current_normalized.csv",
        "portal_url": page_url,
        "file_url": GYGI_GZ_URL,
        "sha256": gz_sha,
        "bytes": gz.stat().st_size,
        "note": "Portal catalog links the Gygi publications page. This is the gzip served from that site.",
    }

    extract_path = out / "depmap24q4_selected_genes_all_models.csv"
    extract_meta = extract_rna(cache / "OmicsExpressionProteinCodingGenesTPMLogp1.csv", extract_path)
    newer = []
    for row in rows:
        if row.get("release") not in NEWER_RELEASES:
            continue
        name = row.get("filename") or ""
        if name == "Model.csv" or ("TPMLogp1" in name and "ProteinCoding" in name):
            url = (row.get("url") or "").strip()
            newer.append(
                {
                    "release": row.get("release"),
                    "filename": name,
                    "url_published": bool(url.startswith("http")),
                }
            )

    catalog_out = out / "portal_catalog_rows_used.csv"
    with catalog_out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["release", "release_date", "filename", "url", "md5_hash"])
        w.writeheader()
        for row in used:
            w.writerow({k: row.get(k, "") for k in w.fieldnames})

    releases = sorted({r["release"] for r in rows if r["release"].startswith("DepMap Public")})
    manifest = {
        "portal_catalog": CATALOG_URL,
        "newest_depmap_public_release_in_catalog": releases[-1] if releases else None,
        "rna_release_used": RNA_RELEASE,
        "why_not_newer_rna": (
            "DepMap Public 25Q2, 25Q3, and 26Q1 expression rows are in the catalog, "
            "but the url column is blank. 24Q4 is the newest public release with direct file URLs."
        ),
        "newer_release_expression_urls_published": newer,
        "any_25q_26q_expression_url": any(item["url_published"] for item in newer),
        "expression_scale": "log2(TPM+1)",
        "files": files_meta,
        "rna_extract": extract_meta,
    }
    (out / "download_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
