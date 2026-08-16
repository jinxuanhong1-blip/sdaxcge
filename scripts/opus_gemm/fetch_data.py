#!/usr/bin/env python3
"""Download the processed GEO supplementary files used by this slice and record
a verifiable manifest (URL, bytes, md5, download time).

Only processed / summarised files are fetched (count, TPM, FPKM, CPM or VST
matrices, or tars of per-sample matrices) - no raw reads. The registry below is
the single source of truth for what enters the analysis.

Usage:
    python3 scripts/opus_gemm/fetch_data.py --outdir results/opus_gemm/data/raw \
        --manifest results/opus_gemm/data_manifest.tsv
    python3 scripts/opus_gemm/fetch_data.py ... --only GSE309199 GSE297817
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import os
import sys
import time
import urllib.error
import urllib.request

FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"

# accession -> supplementary file names (under <series>/suppl/)
REGISTRY: dict[str, list[str]] = {
    # --- in vivo tumours with an ICI arm -------------------------------------
    "GSE309199": ["GSE309199_Mouse_Azam.TPMcalculator.raw_counts.tsv.gz"],
    "GSE297817": ["GSE297817_RAW.tar"],
    "GSE297818": ["GSE297818_RAW.tar"],
    "GSE208614": [
        "GSE208614_201102_rawCounts.txt.gz",
        "GSE208614_210513_rawCounts.txt.gz",
        "GSE208614_220427_rawCounts.txt.gz",
    ],
    "GSE262975": ["GSE262975_CountsTable.txt.gz"],
    "GSE114601": ["GSE114601_counts.raw.csv.gz"],
    "GSE182228": ["GSE182228_RAW.tar"],
    "GSE260596": ["GSE260596_AllSamples_Genes_ReadCounts.txt.gz"],
    "GSE246922": ["GSE246922_KP_RNA_counts_vst.csv.gz", "GSE246922_LLC1_RNA_counts_vst.csv.gz"],
    "GSE197260": ["GSE197260_RNAseqTPM_MM_EGFR-TKI-CD8.txt.gz"],
    "GSE330658": ["GSE330658_RAW.tar"],
    # --- baseline / mechanism ------------------------------------------------
    "GSE309192": ["GSE309192_Mouse_Azam.TPMcalculator.raw_counts.tsv.gz"],
    "GSE274351": ["GSE274351_expredata_TPM_gene.txt.gz"],
    "GSE274352": [
        "GSE274352_normalizedcounts_genes_IFNB_vs_empty.tsv.gz",
        "GSE274352_normalizedcounts_genes_STING_vs_emtpy.tsv.gz",
    ],
    "GSE236258": ["GSE236258_RAW.tar"],
    "GSE241978": ["GSE241978_2020-07-21_Sherr_analysis_CMT_KO_vs_Cas9Ctrl.xlsx"],
    "GSE217405": ["GSE217405_mEGFR_del19.1_del19.2_L860R.1_CPM.xlsx"],
    "GSE303940": ["GSE303940_FPKMs_allSamples.txt.gz"],
    # --- single cell (344SQ parental vs anti-PD-1-resistant derivative) ------
    "GSE285606": ["GSE285606_RAW.tar"],
    # --- additional KRAS/LKB1 and KP NSCLC ICI / genotype series -------------
    "GSE169194": [
        "GSE169194_annotated_raw_count.txt.gz",
        "GSE169194_annotated_log2_norm_count.txt.gz",
    ],
    "GSE295685": ["GSE295685_TPM.txt.gz", "GSE295685_expected_count.txt.gz"],
    "GSE137396": [
        "GSE137396_Raw_genetable_GEMMnodule.txt.gz",
        "GSE137396_Normalized_logtransformed_medcentred_genetable_GEMMnodule.txt.gz",
    ],
    "GSE175479": ["GSE175479_Raw_gene_count_matrix.txt.gz"],
}


def url_for(acc: str, name: str) -> str:
    return f"{FTP}/{acc[:-3]}nnn/{acc}/suppl/{name}"


def download(url: str, dest: str, tries: int = 5) -> None:
    tmp = dest + ".part"
    last: Exception | None = None
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "opus-gemm/1.0"})
            with urllib.request.urlopen(req, timeout=600) as fh, open(tmp, "wb") as out:
                while chunk := fh.read(1 << 20):
                    out.write(chunk)
            os.replace(tmp, dest)
            return
        except urllib.error.HTTPError as exc:
            raise RuntimeError(f"HTTP {exc.code} for {url}") from exc
        except Exception as exc:
            last = exc
            print(f"  retry {attempt + 1}: {exc}", file=sys.stderr)
            time.sleep(4 * 2 ** attempt)
    raise RuntimeError(f"download failed {url}: {last}")


def md5(path: str) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        while chunk := fh.read(1 << 20):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="results/opus_gemm/data/raw")
    ap.add_argument("--manifest", default="results/opus_gemm/data_manifest.tsv")
    ap.add_argument("--only", nargs="*", default=None)
    args = ap.parse_args()

    accs = args.only or sorted(REGISTRY)
    rows: list[dict] = []
    if os.path.exists(args.manifest):
        with open(args.manifest) as fh:
            rows = [r for r in csv.DictReader(fh, delimiter="\t")]

    by_key = {(r["accession"], r["file"]): r for r in rows}
    total = 0
    for acc in accs:
        for name in REGISTRY[acc]:
            d = os.path.join(args.outdir, acc)
            os.makedirs(d, exist_ok=True)
            dest = os.path.join(d, name)
            url = url_for(acc, name)
            if not os.path.exists(dest):
                print(f"[get] {acc}/{name}", flush=True)
                download(url, dest)
            size = os.path.getsize(dest)
            digest = md5(dest)
            total += size
            by_key[(acc, name)] = {
                "accession": acc,
                "file": name,
                "url": url,
                "bytes": str(size),
                "md5": digest,
                "retrieved_utc": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            }
            print(f"      {size / 1e6:9.2f} MB  md5={digest}  {name}")

    out_rows = [by_key[k] for k in sorted(by_key)]
    os.makedirs(os.path.dirname(args.manifest) or ".", exist_ok=True)
    with open(args.manifest, "w", newline="") as fh:
        w = csv.DictWriter(
            fh, fieldnames=["accession", "file", "url", "bytes", "md5", "retrieved_utc"], delimiter="\t"
        )
        w.writeheader()
        w.writerows(out_rows)
    grand = sum(int(r["bytes"]) for r in out_rows)
    print(f"\nmanifest: {args.manifest}  files={len(out_rows)}  total={grand / 1e9:.3f} GB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
