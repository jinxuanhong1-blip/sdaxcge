#!/usr/bin/env python3
"""Download public GSE126044 counts/phenotype and documented gene-set files."""

from __future__ import annotations

import hashlib
import json
import tarfile
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "data" / "rework" / "B4_wave2"
DATA.mkdir(parents=True, exist_ok=True)

UA = "B4_wave2/1.0 (public GEO recompute)"

URLS = {
    "GSE126044_counts.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126044/"
        "suppl/GSE126044_counts.txt.gz"
    ),
    "GSE126044_series_matrix.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126044/"
        "matrix/GSE126044_series_matrix.txt.gz"
    ),
    "kegg_hsa04530.txt": "https://rest.kegg.jp/get/hsa04530",
    "REACTOME_TIGHT_JUNCTION_INTERACTIONS.json": (
        "https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp"
        "?geneSetName=REACTOME_TIGHT_JUNCTION_INTERACTIONS&fileType=json"
    ),
    "GOCC_TIGHT_JUNCTION.json": (
        "https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp"
        "?geneSetName=GOCC_TIGHT_JUNCTION&fileType=json"
    ),
    "REACTOME_KERATINIZATION.json": (
        "https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp"
        "?geneSetName=REACTOME_KERATINIZATION&fileType=json"
    ),
}

EXTRA = DATA / "extra"
EXTRA_URLS = {
    "GSE135222_exp.tsv.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135222/"
        "suppl/GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz"
    ),
    "GSE135222_series_matrix.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135222/"
        "matrix/GSE135222_series_matrix.txt.gz"
    ),
    "GSE207422_log2TPM.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/"
        "suppl/GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz"
    ),
    "GSE207422_metadata.xlsx": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/"
        "suppl/GSE207422_NSCLC_bulk_RNAseq_metadata.xlsx"
    ),
    "GSE190265_TPM.csv.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190265/"
        "suppl/GSE190265_TPM_France3.csv.gz"
    ),
    "GSE190265_samples.csv.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190265/"
        "suppl/GSE190265_samples_info_France3.csv.gz"
    ),
    "GSE166449_TPM.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166449/"
        "suppl/GSE166449_Raw_gene_TPM_matrix.txt.gz"
    ),
    "GSE166449_series_matrix.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166449/"
        "matrix/GSE166449_series_matrix.txt.gz"
    ),
    "GSE190266_TPM.csv.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190266/"
        "suppl/GSE190266_TPM_France4.csv.gz"
    ),
    "GSE190266_series_matrix.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190266/"
        "matrix/GSE190266_series_matrix.txt.gz"
    ),
    "GSE161537_log2cpm.csv.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE161nnn/GSE161537/"
        "suppl/GSE161537_nivobio_log2cpm.csv.gz"
    ),
    "GSE161537_series_matrix.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE161nnn/GSE161537/"
        "matrix/GSE161537_series_matrix.txt.gz"
    ),
    "TCGA_LUAD_HiSeqV2.gz": (
        "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/"
        "TCGA.LUAD.sampleMap%2FHiSeqV2.gz"
    ),
    "TCGA_LUSC_HiSeqV2.gz": (
        "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/"
        "TCGA.LUSC.sampleMap%2FHiSeqV2.gz"
    ),
}

TIDYESTIMATE_URL = (
    "https://cran.r-project.org/src/contrib/tidyestimate_1.1.1.tar.gz"
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(url: str, dest: Path, timeout: int = 180) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=timeout) as r, dest.open("wb") as out:
        out.write(r.read())


def extract_estimate_from_tidyestimate(tarball: Path) -> None:
    import pyreadr

    with tarfile.open(tarball, "r:gz") as tf:
        member = tf.getmember("tidyestimate/data/gene_sets.rda")
        extracted = DATA / "tidyestimate_gene_sets.rda"
        with tf.extractfile(member) as src, extracted.open("wb") as dst:
            dst.write(src.read())
    df = pyreadr.read_r(str(extracted))["gene_sets"]
    (DATA / "ESTIMATE_stromal_signature.txt").write_text(
        "\n".join(df["stromal_signature"].tolist()) + "\n"
    )
    (DATA / "ESTIMATE_immune_signature.txt").write_text(
        "\n".join(df["immune_signature"].tolist()) + "\n"
    )


def main() -> int:
    manifest: dict = {"files": {}, "urls": URLS}
    for name, url in URLS.items():
        dest = DATA / name
        if dest.exists() and dest.stat().st_size > 0:
            print(f"KEEP {dest}", flush=True)
        else:
            print(f"GET {url} -> {dest}", flush=True)
            fetch(url, dest)
        manifest["files"][name] = {
            "url": url,
            "sha256": sha256(dest),
            "bytes": dest.stat().st_size,
        }

    EXTRA.mkdir(parents=True, exist_ok=True)
    for name, url in EXTRA_URLS.items():
        dest = EXTRA / name
        if dest.exists() and dest.stat().st_size > 0:
            print(f"KEEP {dest}", flush=True)
        else:
            print(f"GET {url} -> {dest}", flush=True)
            fetch(url, dest)
        manifest["files"][f"extra/{name}"] = {
            "url": url,
            "sha256": sha256(dest),
            "bytes": dest.stat().st_size,
        }

    stromal = DATA / "ESTIMATE_stromal_signature.txt"
    immune = DATA / "ESTIMATE_immune_signature.txt"
    if not (stromal.exists() and immune.exists()):
        tarball = DATA / "tidyestimate_1.1.1.tar.gz"
        print(f"GET {TIDYESTIMATE_URL}", flush=True)
        fetch(TIDYESTIMATE_URL, tarball)
        extract_estimate_from_tidyestimate(tarball)
        manifest["files"]["tidyestimate_1.1.1.tar.gz"] = {
            "url": TIDYESTIMATE_URL,
            "sha256": sha256(tarball),
            "bytes": tarball.stat().st_size,
            "note": "Yoshihara 2013 ESTIMATE stromal/immune 141-gene signatures via tidyestimate 1.1.1",
        }
    manifest["estimate"] = {
        "source": "Yoshihara et al. 2013 Nat Commun 4:2612; gene lists from tidyestimate 1.1.1",
        "n_stromal": len((DATA / "ESTIMATE_stromal_signature.txt").read_text().split()),
        "n_immune": len((DATA / "ESTIMATE_immune_signature.txt").read_text().split()),
    }
    (DATA / "download_manifest.json").write_text(json.dumps(manifest, indent=2))
    print(json.dumps({"n_files": len(manifest["files"]), "estimate": manifest["estimate"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
