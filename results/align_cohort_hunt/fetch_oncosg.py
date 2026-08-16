#!/usr/bin/env python3
"""Download a documented, analysis-ready OncoSG expression subset from cBioPortal."""

import csv
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "data"
PROV = ROOT / "provenance"
BASE = "https://www.cbioportal.org/api"
STUDY = "luad_oncosg_2020"
PROFILE = f"{STUDY}_rna_seq_v2_mrna_median_all_sample_Zscores"
SAMPLE_LIST = f"{STUDY}_rna_seq_v2_mrna"

SIGNATURES = {
    "target": ["TACSTD2"],
    "cd8": ["CD8A", "CD8B", "CCL5", "GZMK", "LCK", "TRAC"],
    "nk": ["NKG7", "KLRD1", "GNLY", "PRF1", "CTSW", "NCR1"],
    "pan_immune": ["PTPRC", "CD3D", "CD3E", "LST1", "CD74", "HLA-DRA"],
    "epithelial_proxy": ["EPCAM", "KRT7", "KRT8", "KRT18", "KRT19", "MUC1"],
}


def get_json(url):
    request = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(request, timeout=60) as response:
        return json.load(response)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    PROV.mkdir(parents=True, exist_ok=True)
    symbols = list(dict.fromkeys(g for genes in SIGNATURES.values() for g in genes))
    records = {}
    urls = []
    missing = []
    for symbol in symbols:
        gene_url = f"{BASE}/genes/{urllib.parse.quote(symbol)}"
        gene = get_json(gene_url)
        data_url = (
            f"{BASE}/molecular-profiles/{PROFILE}/molecular-data"
            f"?entrezGeneId={gene['entrezGeneId']}"
            f"&sampleListId={SAMPLE_LIST}&projection=DETAILED"
        )
        values = get_json(data_url)
        urls.extend([gene_url, data_url])
        if not values:
            missing.append(symbol)
        for row in values:
            records.setdefault(row["sampleId"], {})[symbol] = row["value"]
        time.sleep(0.03)

    clinical_url = (
        f"{BASE}/studies/{STUDY}/clinical-data"
        "?clinicalDataType=SAMPLE&projection=DETAILED"
    )
    clinical = get_json(clinical_url)
    urls.append(clinical_url)
    purity = {
        row["sampleId"]: row["value"]
        for row in clinical
        if row["clinicalAttributeId"] == "PURITY"
    }

    expression_path = OUT / "OncoSG_GIS031_selected_expression_zscores.tsv"
    with expression_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample_id"] + symbols, delimiter="\t")
        writer.writeheader()
        for sample_id in sorted(records):
            writer.writerow({"sample_id": sample_id, **records[sample_id]})

    purity_path = OUT / "OncoSG_GIS031_sample_purity.tsv"
    with purity_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample_id", "purity"])
        writer.writerows(sorted(purity.items()))

    manifest = {
        "study": STUDY,
        "canonical_OncoSG_id": "GIS031",
        "publication": "Chen et al., Nature Genetics 2020; PMID 32015526",
        "raw_controlled_access": ["EGAS00001002941", "EGAD00001004421"],
        "profile": PROFILE,
        "sample_list": SAMPLE_LIST,
        "expression_scale": "cBioPortal RNA-seq V2 RSEM z-scores versus all samples",
        "signatures": SIGNATURES,
        "missing_requested_genes": missing,
        "urls": urls,
    }
    (PROV / "oncosg_api_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Wrote {len(records)} expression samples; {len(purity)} purity values")
    if missing:
        print("Missing:", ", ".join(missing))


if __name__ == "__main__":
    main()
