#!/usr/bin/env python3
"""Write compact presence / inventory tables from already-computed results."""

from __future__ import annotations

import csv
import json
from pathlib import Path

OUT = Path("results/grok_pride")

inventory = [
    {
        "accession": "PXD042091",
        "host": "PRIDE",
        "verified_lung_ICI": True,
        "design": "Human plasma SWATH-MS; metastatic NSCLC on pembrolizumab; pre- and on-treatment (n=64 patients, 171 samples)",
        "paper": "Mondelo-Macía et al. Mol Cell Proteomics 2024; 23:100834; PMID 39216661",
        "processed_used": "2_dat_cs_all.txt; 6_data_tf_response.txt (PRIDE SEARCH, no raw MS)",
        "TACSTD2": "absent",
        "CLDN4": "absent",
        "note": "7115 quantified proteins; 323 baseline R vs NR DEPs at p<0.05 (paper reported 324). Seven-protein signature reproduced.",
    },
    {
        "accession": "PXD059688",
        "host": "PRIDE",
        "verified_lung_ICI": True,
        "design": "Mouse LLC1 syngeneic tumors: control / AA / anti-PD1 / AA+anti-PD1; nano-LC-ESI-MS/MS. Paper reports 6737 proteins.",
        "paper": "Kim et al. Front Immunol 2025; 15:1512605; DOI 10.3389/fimmu.2024.1512605",
        "processed_used": "20230904_Tumor_AA.mzTab.gz only (mgf 7.9G and msf 35G not downloaded)",
        "TACSTD2": "not_in_deposited_mztab",
        "CLDN4": "not_in_deposited_mztab",
        "note": "Deposited mzTab is a 93-protein methyltransferase subset, not the 6737-protein matrix. Full table is paper Supplementary Table S1; Frontiers file URLs were not openly fetchable.",
    },
    {
        "accession": "PXD019061",
        "host": "PanoramaPublic (+ MassIVE MSV000085049 for global MS)",
        "verified_lung_ICI": True,
        "design": "Human FFPE NSCLC; targeted PRM of checkpoint proteins (n=46) + global proteome spectral counts (n=26). ICI-naive tissue; checkpoint-focused.",
        "paper": "Liebler et al. Sci Rep 2020; 10:9805; PMID 32555232",
        "processed_used": "Sci Rep MOESM1–5 (Table S6 global spectral counts). No raw MS.",
        "TACSTD2": "present_global",
        "CLDN4": "absent_global",
        "note": "TACSTD2 in 24/26 tumors (median 10 spectral counts). CLDN4 not among 6819 gene symbols. Targeted panel does not include either protein.",
    },
    {
        "accession": "PXD039141",
        "host": "iProX (IPX0003162000)",
        "verified_lung_ICI": True,
        "design": "Human serum intact glycoproteomics; advanced lung cancer before/during anti-PD-1/PD-L1 (n=12 discovery)",
        "paper": "Acta Biochim Biophys Sin 2024; DOI 10.3724/abbs.2024110",
        "processed_used": "Glycopeptide search txt/xlsx (no raw)",
        "TACSTD2": "absent",
        "CLDN4": "absent",
        "note": "119 serum glycoprotein accessions; Ig/plasma-dominant. Membrane TACSTD2/CLDN4 not expected and not found.",
    },
    {
        "accession": "PXD019573",
        "host": "PRIDE",
        "verified_lung_ICI": True,
        "design": "Human NSCLC proteasome/degradome (MAPP); PA200/PSME4 vs immunoproteasome; Durvalumab response association in paper",
        "paper": "Javitt et al. Nat Cancer 2023; 4:629-647; also PXD028364 / PXD037365 (A549 immunopeptidome)",
        "processed_used": "Nature Cancer source xlsx (proteasome-subunit focused). No raw MS.",
        "TACSTD2": "not_in_paper_tables",
        "CLDN4": "not_in_paper_tables",
        "note": "Open ICI-related NSCLC proteome, but deposited/paper tables are proteasome-centric, not a TACSTD2/CLDN4 matrix.",
    },
]

# related but NOT verified lung ICI patient/model proteomes
excluded = [
    {
        "accession": "PXD020191",
        "reason": "Lehtiö NSCLC proteogenomics (immune phenotypes) but not ICI-treated cohort",
    },
    {
        "accession": "PXD034772",
        "reason": "Lung immunopeptidome vs T-cell inflammation; not an ICI-treatment proteome",
    },
    {
        "accession": "PXD019774",
        "reason": "HLA immunopeptidome of EGFR-mutant LUAD cell lines + melanoma; not ICI-treated tissue",
    },
    {
        "accession": "PXD044740",
        "reason": "IL-2/IL-12 T-cell phosphoproteome; murine lung tumors but not tumor TACSTD2/CLDN4 proteome",
    },
    {
        "accession": "PXD035347",
        "reason": "ULK1 / ICI resistance; melanoma, not lung",
    },
]

(OUT / "dataset_inventory.json").write_text(json.dumps({"included": inventory, "excluded": excluded}, indent=2))

with (OUT / "presence_table.tsv").open("w") as f:
    f.write("accession\thost\tverified_lung_ICI\tTACSTD2\tCLDN4\tmatrix\n")
    for d in inventory:
        f.write(
            f"{d['accession']}\t{d['host']}\t{d['verified_lung_ICI']}\t{d['TACSTD2']}\t{d['CLDN4']}\t{d['processed_used']}\n"
        )

# compact related-protein table from tf stats
rel_ids = {
    "P16422_EPCAM",
    "P15941_MUC1",
    "O00501_CLD5",
    "Q9NY35_CLDN1",
    "Q16625_OCLN",
    "P56856_CLD18",
    "P43121_MUC18",
    "Q9H3R2_MUC13",
    "Q9UKN1_MUC12",
    "Q8WXI7_MUC16",
}
sig_ids = {
    "Q7Z3C6_ATG9A",
    "Q9UHG0_DCDC2",
    "Q9UPZ3_HPS5",
    "Q4L180_FIL1L",
    "Q9NQ48_LZTL1",
    "Q92696_PGTA",
    "O15020_SPTN2",
}
tf = OUT / "PXD042091_tf_stats.tsv"
if tf.exists():
    with tf.open() as inf, (OUT / "PXD042091_signature_and_junction.tsv").open("w") as outf:
        r = csv.DictReader(inf, delimiter="\t")
        outf.write("class\tprotein\tn_R\tn_NR\tmedian_R\tmedian_NR\tp_val\tfoldchange\n")
        for rec in r:
            cls = None
            if rec["protein"] in sig_ids:
                cls = "signature"
            elif rec["protein"] in rel_ids:
                cls = "junction_or_epithelial"
            if cls:
                outf.write(
                    f"{cls}\t{rec['protein']}\t{rec['n_R']}\t{rec['n_NR']}\t{rec['median_R']}\t{rec['median_NR']}\t{rec['p_val']}\t{rec['foldchange']}\n"
                )

print("wrote inventory + presence tables")
