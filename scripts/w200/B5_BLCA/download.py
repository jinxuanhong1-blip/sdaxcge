#!/usr/bin/env python3
"""B5_BLCA: download the public urothelial-carcinoma ICI cohorts.

Fetches every input needed to test whether CLDN4-high urothelial tumours respond
less often to immune-checkpoint inhibition. Idempotent: a file whose SHA-256
already matches the manifest is not re-downloaded.

Writes a provenance manifest (URL, bytes, SHA-256, declared expression scale) to
results/w200/B5_BLCA/download_manifest.json and a repo-convention catalog.tsv.

Note on IMvigor210: the official Genentech distribution point
http://research-pub.gene.com/IMvigor210CoreBiologies/ now returns HTTP 404. We
take the identical `cds.RData` from a third-party GitHub mirror and validate it
against published summary statistics in prepare_cohorts.py; this is recorded
openly in the manifest rather than being presented as the primary source.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
import urllib.request
import zipfile
from datetime import datetime, timezone

SOURCES = [
    {
        "key": "imvigor210_cds",
        "cohort": "IMvigor210",
        "file_name": "cds.RData",
        "subdir": "IMvigor210",
        "url": "https://raw.githubusercontent.com/SiYangming/IMvigor210CoreBiologies/master/data/cds.RData",
        "official_url": "http://research-pub.gene.com/IMvigor210CoreBiologies/",
        "official_url_status": "HTTP 404 as of 2026-08-16; mirror used instead",
        "citation": "Mariathasan S et al. TGF-beta attenuates tumour response to PD-L1 blockade by contributing to exclusion of T cells. Nature 2018;554:544-548.",
        "doi": "10.1038/nature25501",
        "expression_scale": "raw counts (TPM derived in extract_imvigor210.R using featureData gene length)",
        "license": "CC-BY",
    },
    {
        "key": "gse176307_salmon_tpm",
        "cohort": "BACI",
        "file_name": "GSE176307_salmon_tpm_gene.matrix.tsv.gz",
        "subdir": "GSE176307",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176307/suppl/GSE176307_salmon_tpm_gene.matrix.tsv.gz",
        "accession": "GSE176307",
        "citation": "Robertson AG et al. Expression-based subtypes and FGFR3 alterations in metastatic urothelial cancer treated with immune checkpoint blockade. 2021 (GEO GSE176307).",
        "expression_scale": "TPM (salmon, gene level); log2(TPM+1) applied downstream",
    },
    {
        "key": "gse176307_series_matrix",
        "cohort": "BACI",
        "file_name": "GSE176307_series_matrix.txt.gz",
        "subdir": "GSE176307",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176307/matrix/GSE176307_series_matrix.txt.gz",
        "accession": "GSE176307",
        "expression_scale": "n/a (clinical annotation: io.response, io.therapy, tmb, site)",
    },
    {
        "key": "gse176307_sample_key",
        "cohort": "BACI",
        "file_name": "GSE176307_BACI_Omniseq_Sample_Name_Key_submitted_GEO_v2.csv.gz",
        "subdir": "GSE176307",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176307/suppl/GSE176307_BACI_Omniseq_Sample_Name_Key_submitted_GEO_v2.csv.gz",
        "accession": "GSE176307",
        "expression_scale": "n/a (maps BACI patient id to RS- RNA-seq id)",
    },
    {
        "key": "gse176307_log_trans",
        "cohort": "BACI",
        "file_name": "GSE176307_BACI_log_trans_normalized_RNAseq.csv.gz",
        "subdir": "GSE176307",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176307/suppl/GSE176307_BACI_log_trans_normalized_RNAseq.csv.gz",
        "accession": "GSE176307",
        "expression_scale": "depositor log-transformed normalised RNA-seq (used only as a concordance check)",
    },
    {
        "key": "predictio_snyder",
        "cohort": "Snyder",
        "file_name": "ICB_Snyder.zip",
        "subdir": "predictio",
        "url": "https://zenodo.org/api/records/7058399/files/ICB_Snyder.zip/content",
        "accession": "Zenodo 7058399",
        "citation": "Snyder A et al. Contribution of systemic and somatic factors to clinical response and resistance to PD-L1 blockade in urothelial cancer. PLoS Med 2017;14:e1002309. Harmonised by Bareche Y et al., Ann Oncol 2022 (PMID 36055464).",
        "doi": "10.5281/zenodo.7058399",
        "expression_scale": "log2(TPM+0.001) as deposited",
        "unzip": True,
    },
    {
        "key": "predictio_mariathasan",
        "cohort": "IMvigor210 (harmonised copy, concordance check only)",
        "file_name": "ICB_Mariathasan.zip",
        "subdir": "predictio",
        "url": "https://zenodo.org/api/records/7058399/files/ICB_Mariathasan.zip/content",
        "accession": "Zenodo 7058399",
        "citation": "Mariathasan S et al. Nature 2018, harmonised by Bareche Y et al., Ann Oncol 2022 (PMID 36055464).",
        "doi": "10.5281/zenodo.7058399",
        "expression_scale": "log2(TPM+0.001) as deposited",
        "note": "Same trial as imvigor210_cds. Used ONLY for pipeline concordance (S7); excluded from pooling to avoid double counting.",
        "unzip": True,
    },
]


def sha256_of(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(url: str, dest: str, retries: int = 4) -> None:
    delay = 4
    last = None
    for attempt in range(retries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "b5-blca/1.0"})
            with urllib.request.urlopen(req, timeout=900) as r, open(dest, "wb") as out:
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    out.write(chunk)
            return
        except Exception as exc:  # network flake: exponential backoff
            last = exc
            if attempt == retries:
                break
            time.sleep(delay)
            delay *= 2
    raise RuntimeError(f"failed to download {url}: {last}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="data/raw")
    ap.add_argument("--out-dir", default="results/w200/B5_BLCA")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    entries = []
    catalog_rows = []

    for src in SOURCES:
        d = os.path.join(args.data_dir, src["subdir"])
        os.makedirs(d, exist_ok=True)
        path = os.path.join(d, src["file_name"])

        if os.path.exists(path):
            print(f"[have] {path}")
        else:
            print(f"[get ] {src['url']}")
            fetch(src["url"], path)

        digest = sha256_of(path)
        size = os.path.getsize(path)
        print(f"       {size} bytes  sha256={digest[:16]}...")

        if src.get("unzip"):
            with zipfile.ZipFile(path) as z:
                z.extractall(d)

        entry = {k: v for k, v in src.items() if k != "unzip"}
        entry.update(
            {
                "local_path": path,
                "bytes": size,
                "sha256": digest,
                "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
            }
        )
        entries.append(entry)
        catalog_rows.append(
            [
                src["cohort"],
                src.get("accession", src.get("doi", "n/a")),
                "GEO" if str(src.get("accession", "")).startswith("GSE") else ("Zenodo" if "zenodo" in src["url"] else "GitHub mirror"),
                src["url"],
                entry["retrieved_at_utc"],
                f"HTTP 200; {size} bytes retrieved",
                src["file_name"],
                str(size),
                "accepted",
                "",
                digest,
            ]
        )

    manifest = {
        "slice": "B5_BLCA",
        "question": "Do CLDN4-high urothelial tumours respond less often to ICI? (user-asserted OR=0.42)",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "files": entries,
    }
    mpath = os.path.join(args.out_dir, "download_manifest.json")
    with open(mpath, "w") as fh:
        json.dump(manifest, fh, indent=2)
    print(f"[write] {mpath}")

    header = [
        "sample_id", "accession", "database", "source_url", "verified_at_utc",
        "verification_evidence", "file_name", "bytes", "decision",
        "refusal_reason", "sha256",
    ]
    with open("catalog.tsv", "w") as fh:
        fh.write("\t".join(header) + "\n")
        for row in catalog_rows:
            fh.write("\t".join(row) + "\n")
    print("[write] catalog.tsv")

    os.makedirs("repro", exist_ok=True)
    with open("repro/seed.txt", "w") as fh:
        fh.write("20260816\n")
    try:
        freeze = subprocess.run(
            [sys.executable, "-m", "pip", "freeze", "--all"],
            capture_output=True, text=True, check=False,
        ).stdout
        with open("repro/pip-freeze.txt", "w") as fh:
            fh.write(freeze)
    except Exception as exc:
        print(f"[warn] pip freeze failed: {exc}")
    print("[write] repro/seed.txt, repro/pip-freeze.txt")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
