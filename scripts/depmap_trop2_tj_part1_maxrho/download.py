#!/usr/bin/env python3
"""Download Gygi CCLE protein + DepMap 24Q4 Model/RNA for Part1 bridge max ρ.

Keeps the full Gygi matrix (peptide columns needed for QC). Extracts a slim
RNA panel of locked TJ / epithelial control genes, then deletes the full TPM.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import time
import urllib.request
from pathlib import Path

UA = "sdaxcge-depmap-trop2-tj-part1-maxrho/1.0 (+https://github.com/jinxuanhong1-blip/sdaxcge)"

MODEL_URL = "https://ndownloader.figshare.com/files/51065297"
EXPR_URL = "https://ndownloader.figshare.com/files/51065489"
GYGI_URLS = [
    "https://gygi.hms.harvard.edu/data/ccle/protein_quant_current_normalized.csv.gz",
    "https://gygi.med.harvard.edu/sites/gygi.med.harvard.edu/files/documents/protein_quant_current_normalized.csv.gz",
]
S1_URL = "https://gygi.hms.harvard.edu/data/ccle/Table_S1_Sample_Information.xlsx"

RELEASE = "DepMap Public 24Q4"
DOI = "10.25452/figshare.plus.27993248.v1"

# Locked TJ lists (same as PPT / GSE137244 transfer + TISMO).
TJ_EPITHELIAL = [
    "CLDN1",
    "CLDN3",
    "CLDN4",
    "CLDN7",
    "OCLN",
    "MARVELD2",
    "MARVELD3",
    "TJP1",
    "TJP2",
    "TJP3",
    "F11R",
    "JAM2",
    "JAM3",
    "CGN",
    "CGNL1",
    "CRB3",
    "ILDR1",
    "LSR",
]
TJ_TISMO = ["CLDN3", "CLDN4", "CLDN6", "CLDN7", "CDH1", "F11R", "OCLN"]
CLDN4_TJ_EDGE = [
    "CLDN4",
    "CLDN1",
    "CLDN7",
    "CGNL1",
    "MARVELD2",
    "MARVELD3",
    "TJP1",
    "TJP2",
    "ILDR1",
    "CLDN3",
    "OCLN",
]
CONTROLS = ["TACSTD2", "EPCAM", "CDH1", "KRT8", "KRT18", "KRT19", "VIM"]
FOCAL_SINGLE = ["CLDN1", "CLDN4", "CLDN7", "OCLN", "F11R", "TJP1"]

CORE_GENES = sorted(
    set(CONTROLS) | set(TJ_EPITHELIAL) | set(TJ_TISMO) | set(CLDN4_TJ_EDGE) | set(FOCAL_SINGLE)
)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, dest: Path, min_bytes: int) -> str:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    print(f"GET {url}", flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    h = hashlib.sha256()
    n = 0
    with urllib.request.urlopen(req, timeout=600) as resp, tmp.open("wb") as out:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            out.write(chunk)
            h.update(chunk)
            n += len(chunk)
    if n < min_bytes:
        tmp.unlink(missing_ok=True)
        raise RuntimeError(f"{url} wrote only {n} bytes")
    tmp.replace(dest)
    print(f"  {dest.name} {n} bytes sha256={h.hexdigest()[:12]}", flush=True)
    return h.hexdigest()


def fetch(urls: list[str], dest: Path, min_bytes: int) -> tuple[str, str]:
    last: Exception | None = None
    for url in urls:
        try:
            digest = download(url, dest, min_bytes)
            return url, digest
        except Exception as exc:
            last = exc
            print(f"WARN {url}: {exc}", flush=True)
            dest.unlink(missing_ok=True)
            time.sleep(2)
    raise SystemExit(f"failed to download {dest.name}: {last}")


def symbol_of(col: str) -> str:
    return col.split(" (")[0].strip()


def extract_model_by_gene(src: Path, wanted: set[str], dest: Path) -> dict:
    with src.open("r", newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        found: dict[str, int] = {}
        for i, col in enumerate(header[1:], start=1):
            sym = symbol_of(col)
            if sym in wanted and sym not in found:
                found[sym] = i
        keep = [0] + [found[s] for s in sorted(found)]
        missing = sorted(wanted - set(found))
        print(
            f"extract {src.name}: kept {len(found)}/{len(wanted)}; missing={missing}",
            flush=True,
        )
        dest.parent.mkdir(parents=True, exist_ok=True)
        n_rows = 0
        with dest.open("w", newline="") as out:
            writer = csv.writer(out)
            writer.writerow(["ModelID"] + [symbol_of(header[i]) for i in keep[1:]])
            for row in reader:
                if not row:
                    continue
                writer.writerow([row[i] if i < len(row) else "" for i in keep])
                n_rows += 1
    return {
        "n_models": n_rows,
        "n_genes_kept": len(found),
        "n_genes_wanted": len(wanted),
        "genes_kept": sorted(found),
        "missing_genes": missing,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="data/depmap_trop2_tj_part1_maxrho")
    args = ap.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    gene_sets = {
        "TJ_EPITHELIAL": TJ_EPITHELIAL,
        "TJ_TISMO": TJ_TISMO,
        "CLDN4_TJ_EDGE": CLDN4_TJ_EDGE,
        "CONTROLS": CONTROLS,
        "FOCAL_SINGLE": FOCAL_SINGLE,
    }
    (outdir / "gene_sets.json").write_text(json.dumps(gene_sets, indent=2) + "\n")

    manifest: dict = {
        "release": RELEASE,
        "doi": DOI,
        "core_genes": CORE_GENES,
        "files": {},
    }

    model = outdir / "Model.csv"
    if not model.exists() or model.stat().st_size < 100_000:
        url, digest = fetch([MODEL_URL], model, 100_000)
    else:
        url, digest = MODEL_URL, sha256_file(model)
        print(f"OK exists {model}", flush=True)
    manifest["files"]["Model.csv"] = {
        "url": url,
        "sha256": digest,
        "bytes": model.stat().st_size,
        "role": "DepMap Public 24Q4 Model.csv",
        "citation": f"DepMap Public 24Q4, Figshare+ {DOI} file 51065297",
    }

    s1 = outdir / "Table_S1_Sample_Information.xlsx"
    if not s1.exists() or s1.stat().st_size < 10_000:
        url, digest = fetch([S1_URL], s1, 10_000)
    else:
        url, digest = S1_URL, sha256_file(s1)
        print(f"OK exists {s1}", flush=True)
    manifest["files"]["Table_S1_Sample_Information.xlsx"] = {
        "url": url,
        "sha256": digest,
        "bytes": s1.stat().st_size,
        "role": "Nusinow 2020 Table S1 sample information",
        "citation": "Nusinow et al. Cell 2020;180:387-402.e16",
    }

    gygi = outdir / "protein_quant_current_normalized.csv.gz"
    if not gygi.exists() or gygi.stat().st_size < 5_000_000:
        url, digest = fetch(GYGI_URLS, gygi, 5_000_000)
    else:
        url, digest = GYGI_URLS[0], sha256_file(gygi)
        print(f"OK exists {gygi}", flush=True)
    manifest["files"]["protein_quant_current_normalized.csv.gz"] = {
        "url": url,
        "sha256": digest,
        "bytes": gygi.stat().st_size,
        "role": "Nusinow 2020 / Gygi normalized protein quantitation",
        "citation": "Nusinow et al. Cell 2020;180:387-402.e16",
    }

    expr_full = outdir / "OmicsExpressionProteinCodingGenesTPMLogp1.csv"
    expr_slim = outdir / "expression_panel.csv"
    if not expr_slim.exists() or expr_slim.stat().st_size < 50_000:
        if not expr_full.exists() or expr_full.stat().st_size < 400_000_000:
            url, digest = fetch([EXPR_URL], expr_full, 400_000_000)
        else:
            url, digest = EXPR_URL, sha256_file(expr_full)
        meta = extract_model_by_gene(expr_full, set(CORE_GENES), expr_slim)
        expr_full.unlink(missing_ok=True)
    else:
        url, digest = EXPR_URL, None
        with expr_slim.open() as fh:
            header = next(csv.reader(fh))
        meta = {
            "n_models": sum(1 for _ in expr_slim.open()) - 1,
            "genes_kept": header[1:],
            "missing_genes": sorted(set(CORE_GENES) - set(header[1:])),
            "note": "slim extract already present",
        }
        print(f"OK exists {expr_slim}", flush=True)
    manifest["files"]["expression"] = {
        "url": url,
        "sha256": digest,
        "extract": meta,
        "role": "DepMap 24Q4 OmicsExpressionProteinCodingGenesTPMLogp1 slim extract",
    }

    (outdir / "download_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print("manifest written", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
