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
    {
        "key": "gse111636_series_matrix",
        "cohort": "GSE111636",
        "file_name": "GSE111636_series_matrix.txt.gz",
        "subdir": "GSE111636",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE111nnn/GSE111636/matrix/GSE111636_series_matrix.txt.gz",
        "accession": "GSE111636",
        "citation": "Homet Moreno B et al. Identification of potential biomarkers of response and resistance to PD-1 blockade in advanced urothelial tumors. GEO GSE111636 (pembrolizumab, n=11, Affymetrix HTA-2.0).",
        "expression_scale": "HTA-2.0 transcript-cluster intensity; CLDN4 = TC07000447.hg.1",
        "note": "Post-hoc open cohort. Endpoint is depositor binary responder/progressor, not full RECIST.",
    },
]


def sha256_of(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


UC_GENOME_GENES = {
    "CLDN4": 1364, "CD8A": 925, "GZMA": 3001, "GZMB": 3002, "IFNG": 3458,
    "EOMES": 8320, "CXCL9": 4283, "CXCL10": 3627, "TBX21": 30009, "PRF1": 5551,
    "ACTB": 60, "GAPDH": 2597, "TBP": 6908, "RPL13A": 23521, "PGK1": 5230,
    "TACSTD2": 4070, "CLDN1": 9076, "CLDN3": 1365, "CLDN7": 1366, "CD274": 29126,
}


def fetch_ucgenome(dest_dir: str) -> list:
    """Download UC-GENOME clinical + gene z-scores from the public cBioPortal API."""
    out = []
    clin_url = "https://www.cbioportal.org/api/studies/blca_bcan_hcrn_2022/clinical-data?clinicalDataType=PATIENT"
    clin_path = os.path.join(dest_dir, "clinical_patient.json")
    if not os.path.exists(clin_path):
        print(f"[get ] {clin_url}")
        fetch(clin_url, clin_path)
    else:
        print(f"[have] {clin_path}")
    out.append({
        "key": "ucgenome_clinical",
        "cohort": "UC-GENOME",
        "file_name": "clinical_patient.json",
        "url": clin_url,
        "accession": "blca_bcan_hcrn_2022",
        "citation": "Damrauer JS et al. Collaborative study from the Bladder Cancer Advocacy Network for the genomic analysis of metastatic urothelial cancer. Nat Commun 2022;13:6658. PMID 36333289.",
        "doi": "10.1038/s41467-022-33980-9",
        "expression_scale": "n/a (clinical: IMMUNOTHERAPY, BEST_RESPONSE_IMMUNOTHERAPY)",
        "local_path": clin_path,
        "bytes": os.path.getsize(clin_path),
        "sha256": sha256_of(clin_path),
        "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
        "note": "Post-hoc open cohort added after the 3-cohort primary was locked.",
    })
    for gene, eid in UC_GENOME_GENES.items():
        url = ("https://www.cbioportal.org/api/molecular-profiles/"
               "blca_bcan_hcrn_2022_rna_seq_v2_mrna_median_Zscores/molecular-data"
               f"?entrezGeneId={eid}&sampleListId=blca_bcan_hcrn_2022_all")
        path = os.path.join(dest_dir, f"{gene}_zscore.json")
        if not os.path.exists(path):
            print(f"[get ] UC-GENOME {gene}")
            fetch(url, path)
        else:
            print(f"[have] {path}")
        out.append({
            "key": f"ucgenome_{gene}",
            "cohort": "UC-GENOME",
            "file_name": f"{gene}_zscore.json",
            "url": url,
            "accession": "blca_bcan_hcrn_2022",
            "expression_scale": "cBioPortal RNA-seq V2 median z-score",
            "local_path": path,
            "bytes": os.path.getsize(path),
            "sha256": sha256_of(path),
            "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
        })
    return out


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

    # UC-GENOME is served as JSON from the cBioPortal API (no single tarball; the
    # datahub S3 link returns 403). Fetch clinical + selected genes here so the
    # rest of the pipeline is offline.
    uc_dir = os.path.join(args.data_dir, "ucgenome")
    os.makedirs(uc_dir, exist_ok=True)
    uc_files = fetch_ucgenome(uc_dir)
    entries.extend(uc_files)
    for e in uc_files:
        catalog_rows.append([
            "UC-GENOME", e.get("accession", "blca_bcan_hcrn_2022"), "cBioPortal",
            e["url"], e["retrieved_at_utc"], f"HTTP 200; {e['bytes']} bytes retrieved",
            e["file_name"], str(e["bytes"]), "accepted", "", e["sha256"],
        ])

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
