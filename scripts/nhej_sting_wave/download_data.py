#!/usr/bin/env python3
"""Download open TCGA inputs for the NHEJ–STING logic wave.

Nothing here is controlled-access. Raw matrices stay outside the repo
(default /tmp/nhej_sting_data, override with NHEJ_STING_DATA).

Sources
-------
* UCSC Xena GDC hub: TCGA-OV and TCGA-LUAD STAR log2(TPM+1), GENCODE v36
  probemap.
* cBioPortal Firehose Legacy (ov_tcga, luad_tcga): MUTATION_COUNT, CLDN4
  RNA-seq V2 RSEM, and (OV only) CPTAC CLDN4 protein. This is the study
  Yamamoto et al., Mol Cancer Ther 2022 used for the mutational-burden claim.
* MCP-counter marker file (Becht et al., Genome Biol 2016).
* Enrichr copies of MSigDB Hallmark 2020 and KEGG 2021 Human, used only to
  freeze the IFN and KEGG NHEJ gene sets.
* PanCanAtlas ABSOLUTE purity (open GDC file).
"""

import json
import os
import sys
import urllib.request

DATA_DIR = os.environ.get("NHEJ_STING_DATA", "/tmp/nhej_sting_data")
GDC_HUB = "https://gdc-hub.s3.us-east-1.amazonaws.com/download"
CBIOPORTAL = "https://www.cbioportal.org/api"

FILES = {
    "TCGA-OV.star_tpm.tsv.gz": f"{GDC_HUB}/TCGA-OV.star_tpm.tsv.gz",
    "TCGA-LUAD.star_tpm.tsv.gz": f"{GDC_HUB}/TCGA-LUAD.star_tpm.tsv.gz",
    "gencode.v36.probemap": f"{GDC_HUB}/gencode.v36.annotation.gtf.gene.probemap",
    "mcpcounter_genes.txt": (
        "https://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt"
    ),
    "msigdb_hallmark_2020.txt": (
        "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=MSigDB_Hallmark_2020"
    ),
    "kegg_2021_human.txt": (
        "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2021_Human"
    ),
    "tcga_absolute_purity.txt": (
        "https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5"
    ),
}

# cBioPortal molecular / clinical pulls. CLDN4 Entrez 1364.
PULLS = [
    # Today's Firehose Legacy MUTATION_COUNT. In this extract the OV maximum
    # is far below the 1,899 cited in Yamamoto 2022 Fig. 4F.
    ("ov_mutation_count_firehose.tsv", "ov_tcga", "clinical", None),
    ("luad_mutation_count_firehose.tsv", "luad_tcga", "clinical", None),
    # Pan-Cancer Atlas MUTATION_COUNT. OV maximum is 1,893, which is the
    # public count that still matches the published 2–1,899 range.
    ("ov_mutation_count_pancan.tsv", "ov_tcga_pan_can_atlas_2018", "clinical", None),
    ("luad_mutation_count_pancan.tsv", "luad_tcga_pan_can_atlas_2018", "clinical", None),
    ("ov_cldn4_rsem.tsv", "ov_tcga_rna_seq_v2_mrna", "molecular", 1364),
    ("luad_cldn4_rsem.tsv", "luad_tcga_rna_seq_v2_mrna", "molecular", 1364),
    ("ov_cldn4_rsem_pancan.tsv", "ov_tcga_pan_can_atlas_2018_rna_seq_v2_mrna", "molecular", 1364),
    ("luad_cldn4_rsem_pancan.tsv", "luad_tcga_pan_can_atlas_2018_rna_seq_v2_mrna", "molecular", 1364),
    ("ov_cldn4_protein.tsv", "ov_tcga_protein_quantification", "molecular", 1364),
]


def fetch(name, url):
    dest = os.path.join(DATA_DIR, name)
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        print(f"[skip] {name}")
        return
    print(f"[get ] {url}")
    tmp = dest + ".part"
    urllib.request.urlretrieve(url, tmp)
    os.replace(tmp, dest)
    print(f"[done] {name}  size={os.path.getsize(dest):,}")


def _get_json(url):
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=180) as resp:
        return json.loads(resp.read().decode())


def pull_clinical(study):
    url = (
        f"{CBIOPORTAL}/studies/{study}/clinical-data"
        "?attributeId=MUTATION_COUNT&clinicalDataType=SAMPLE"
        "&projection=SUMMARY&pageSize=100000&pageNumber=0"
    )
    rows = _get_json(url)
    return [(r["sampleId"], r["patientId"], r["value"]) for r in rows]


def pull_molecular(profile, entrez):
    url = (
        f"{CBIOPORTAL}/molecular-profiles/{profile}/molecular-data"
        f"?entrezGeneId={entrez}&sampleListId={profile}"
        "&projection=SUMMARY&pageSize=100000&pageNumber=0"
    )
    rows = _get_json(url)
    return [(r["sampleId"], r["patientId"], r["value"]) for r in rows]


def write_tsv(name, rows):
    dest = os.path.join(DATA_DIR, name)
    with open(dest, "w") as fh:
        fh.write("sample_id\tpatient_id\tvalue\n")
        for sample_id, patient_id, value in rows:
            fh.write(f"{sample_id}\t{patient_id}\t{value}\n")
    print(f"[done] {name}  n={len(rows):,}")


def main():
    os.makedirs(DATA_DIR, exist_ok=True)
    for name, url in FILES.items():
        fetch(name, url)
    for name, key, kind, entrez in PULLS:
        dest = os.path.join(DATA_DIR, name)
        if os.path.exists(dest) and os.path.getsize(dest) > 0:
            print(f"[skip] {name}")
            continue
        if kind == "clinical":
            rows = pull_clinical(key)
        else:
            rows = pull_molecular(key, entrez)
        write_tsv(name, rows)
    print(f"\nAll files in {DATA_DIR}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
