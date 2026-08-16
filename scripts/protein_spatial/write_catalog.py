#!/usr/bin/env python3
"""Write notes/protein_spatial_catalog.tsv from local downloads + analysis outputs."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from config import DATA, NOTES, RESULTS

ROWS = [
    {
        "dataset": "CPTAC_LUAD_v1.2",
        "accession": "CPTAC freeze v1.2 LUAD",
        "modality": "TMT_proteomics_gene",
        "species": "Homo sapiens",
        "cancer": "LUAD",
        "ici_related": "no_treatment_naive",
        "file_name": "LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt",
        "url": "https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/LUAD/LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt",
        "local_rel": "cptac/LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt",
        "downloaded": "yes",
        "skipped_reason": "",
        "notes": "Open processed tumor protein matrix. ENSG00000184292 and ENSG00000189143 present.",
    },
    {
        "dataset": "CPTAC_LSCC_v1.2",
        "accession": "CPTAC freeze v1.2 LSCC",
        "modality": "TMT_proteomics_gene",
        "species": "Homo sapiens",
        "cancer": "LSCC",
        "ici_related": "no_treatment_naive",
        "file_name": "LSCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt",
        "url": "https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/LSCC/LSCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt",
        "local_rel": "cptac/LSCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt",
        "downloaded": "yes",
        "skipped_reason": "",
        "notes": "LSCC sibling of the specified LUAD protein matrix.",
    },
    {
        "dataset": "CPTAC_LUAD_v1.2",
        "accession": "CPTAC freeze v1.2 LUAD",
        "modality": "phenotype_immune_scores",
        "species": "Homo sapiens",
        "cancer": "LUAD",
        "ici_related": "no_treatment_naive",
        "file_name": "LUAD_phenotype.txt",
        "url": "https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/LUAD/LUAD_phenotype.txt",
        "local_rel": "cptac/LUAD_phenotype.txt",
        "downloaded": "yes",
        "skipped_reason": "",
        "notes": "Same freeze. ESTIMATE/xCell/CIBERSORT (RNA-derived). No ICI labels.",
    },
    {
        "dataset": "CPTAC_LSCC_v1.2",
        "accession": "CPTAC freeze v1.2 LSCC",
        "modality": "phenotype_immune_scores",
        "species": "Homo sapiens",
        "cancer": "LSCC",
        "ici_related": "no_treatment_naive",
        "file_name": "LSCC_phenotype.txt",
        "url": "https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/LSCC/LSCC_phenotype.txt",
        "local_rel": "cptac/LSCC_phenotype.txt",
        "downloaded": "yes",
        "skipped_reason": "",
        "notes": "Same freeze phenotype sibling.",
    },
    {
        "dataset": "CPTAC_LUAD_v1.2",
        "accession": "CPTAC freeze v1.2 LUAD",
        "modality": "survival",
        "species": "Homo sapiens",
        "cancer": "LUAD",
        "ici_related": "no_treatment_naive",
        "file_name": "LUAD_survival.txt",
        "url": "https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/LUAD/LUAD_survival.txt",
        "local_rel": "cptac/LUAD_survival.txt",
        "downloaded": "yes",
        "skipped_reason": "",
        "notes": "OS/PFS after resection; not ICI survival.",
    },
    {
        "dataset": "CPTAC_LSCC_v1.2",
        "accession": "CPTAC freeze v1.2 LSCC",
        "modality": "survival",
        "species": "Homo sapiens",
        "cancer": "LSCC",
        "ici_related": "no_treatment_naive",
        "file_name": "LSCC_survival.txt",
        "url": "https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/LSCC/LSCC_survival.txt",
        "local_rel": "cptac/LSCC_survival.txt",
        "downloaded": "yes",
        "skipped_reason": "",
        "notes": "OS/PFS after resection; not ICI survival.",
    },
    {
        "dataset": "CPTAC_LUAD_v1.2",
        "accession": "CPTAC freeze v1.2 LUAD",
        "modality": "clinical_meta",
        "species": "Homo sapiens",
        "cancer": "LUAD",
        "ici_related": "no_treatment_naive",
        "file_name": "LUAD_meta.txt",
        "url": "https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/LUAD/LUAD_meta.txt",
        "local_rel": "cptac/LUAD_meta.txt",
        "downloaded": "yes",
        "skipped_reason": "",
        "notes": "Age/Sex/Stage/smoking/driver mutations. No ICI field.",
    },
    {
        "dataset": "CPTAC_LSCC_v1.2",
        "accession": "CPTAC freeze v1.2 LSCC",
        "modality": "clinical_meta",
        "species": "Homo sapiens",
        "cancer": "LSCC",
        "ici_related": "no_treatment_naive",
        "file_name": "LSCC_meta.txt",
        "url": "https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/LSCC/LSCC_meta.txt",
        "local_rel": "cptac/LSCC_meta.txt",
        "downloaded": "yes",
        "skipped_reason": "",
        "notes": "Same freeze clinical sibling.",
    },
    {
        "dataset": "PXD042091",
        "accession": "PXD042091",
        "modality": "SWATH_plasma_quant",
        "species": "Homo sapiens",
        "cancer": "metastatic_NSCLC",
        "ici_related": "yes_immunotherapy_plasma",
        "file_name": "6_data_tf_response.txt",
        "url": "https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/02/PXD042091/6_data_tf_response.txt",
        "local_rel": "pxd042091/6_data_tf_response.txt",
        "downloaded": "yes",
        "skipped_reason": "",
        "notes": "Processed SWATH table (responder vs nonresponder). TACSTD2/CLDN4 not quantified.",
    },
    {
        "dataset": "PXD042091",
        "accession": "PXD042091",
        "modality": "SWATH_plasma_quant",
        "species": "Homo sapiens",
        "cancer": "metastatic_NSCLC",
        "ici_related": "yes_immunotherapy_plasma",
        "file_name": "2_dat_cs_all.txt",
        "url": "https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/02/PXD042091/2_dat_cs_all.txt",
        "local_rel": "pxd042091/2_dat_cs_all.txt",
        "downloaded": "yes",
        "skipped_reason": "",
        "notes": "All-sample SWATH matrix. TACSTD2/CLDN4 not quantified.",
    },
    {
        "dataset": "PXD042091",
        "accession": "PXD042091",
        "modality": "spectral_library",
        "species": "Homo sapiens",
        "cancer": "NSCLC_library",
        "ici_related": "yes_immunotherapy_plasma",
        "file_name": "NSClibrary.txt",
        "url": "https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/02/PXD042091/NSClibrary.txt",
        "local_rel": "pxd042091/NSClibrary.txt",
        "downloaded": "yes",
        "skipped_reason": "",
        "notes": "Library includes TACD2_HUMAN and CLD4_HUMAN peptides; not sample quantification.",
    },
    {
        "dataset": "PXD042091",
        "accession": "PXD042091",
        "modality": "raw_SWATH",
        "species": "Homo sapiens",
        "cancer": "metastatic_NSCLC",
        "ici_related": "yes_immunotherapy_plasma",
        "file_name": "*.wiff / *.wiff.scan",
        "url": "https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/02/PXD042091/",
        "local_rel": "",
        "downloaded": "no",
        "skipped_reason": "raw_MS",
        "notes": "Skipped raw vendor files as required.",
    },
    {
        "dataset": "PXD059688",
        "accession": "PXD059688",
        "modality": "mzTab_quant",
        "species": "Mus musculus",
        "cancer": "LLC_syngeneic",
        "ici_related": "yes_antiPD1_plus_ascorbate",
        "file_name": "20230904_Tumor_AA.mzTab.gz",
        "url": "https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/01/PXD059688/20230904_Tumor_AA.mzTab.gz",
        "local_rel": "pxd059688/20230904_Tumor_AA.mzTab.gz",
        "downloaded": "yes",
        "skipped_reason": "",
        "notes": "Only processed table <2GB. 93 PRT rows, methyltransferase-focused. Tacstd2/Cldn4 absent.",
    },
    {
        "dataset": "PXD059688",
        "accession": "PXD059688",
        "modality": "peaklist",
        "species": "Mus musculus",
        "cancer": "LLC_syngeneic",
        "ici_related": "yes_antiPD1_plus_ascorbate",
        "file_name": "20230904_Tumor_AA.mgf",
        "url": "https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/01/PXD059688/20230904_Tumor_AA.mgf",
        "local_rel": "",
        "downloaded": "no",
        "skipped_reason": "file_gt_2GB_8492466377",
        "notes": "Skipped >2GB.",
    },
    {
        "dataset": "PXD059688",
        "accession": "PXD059688",
        "modality": "PD_msf",
        "species": "Mus musculus",
        "cancer": "LLC_syngeneic",
        "ici_related": "yes_antiPD1_plus_ascorbate",
        "file_name": "20230904_Tumor_AA.msf",
        "url": "https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/01/PXD059688/20230904_Tumor_AA.msf",
        "local_rel": "",
        "downloaded": "no",
        "skipped_reason": "file_gt_2GB_37928718336",
        "notes": "Skipped >2GB. Would be the full Proteome Discoverer result.",
    },
    {
        "dataset": "PXD059688",
        "accession": "PXD059688",
        "modality": "raw_MS",
        "species": "Mus musculus",
        "cancer": "LLC_syngeneic",
        "ici_related": "yes_antiPD1_plus_ascorbate",
        "file_name": "20230828_Tumor_*.raw",
        "url": "https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/01/PXD059688/",
        "local_rel": "",
        "downloaded": "no",
        "skipped_reason": "raw_MS",
        "notes": "Skipped raw Thermo files.",
    },
    {
        "dataset": "GSE271689",
        "accession": "GSE271689",
        "modality": "GeoMx_WTA_DCC",
        "species": "Homo sapiens",
        "cancer": "NSCLC",
        "ici_related": "yes_first_line_immunotherapy",
        "file_name": "GSE271689_RAW.tar",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE271nnn/GSE271689/suppl/GSE271689_RAW.tar",
        "local_rel": "gse271689/GSE271689_RAW.tar",
        "downloaded": "yes",
        "skipped_reason": "",
        "notes": "586 DCC.gz files. Series matrix has no OS/PFS. WTA panel includes TACSTD2 and CLDN4.",
    },
    {
        "dataset": "GSE271689",
        "accession": "GSE271689",
        "modality": "GEO_metadata",
        "species": "Homo sapiens",
        "cancer": "NSCLC",
        "ici_related": "yes_first_line_immunotherapy",
        "file_name": "GSE271689_series_matrix.txt.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE271nnn/GSE271689/matrix/GSE271689_series_matrix.txt.gz",
        "local_rel": "gse271689/GSE271689_series_matrix.txt.gz",
        "downloaded": "yes",
        "skipped_reason": "",
        "notes": "Empty expression table. Characteristics: spotid, tissue, cell type (CK/CD45/CD68), treatment=immunotherapy. No OS.",
    },
    {
        "dataset": "GSE271689",
        "accession": "GSE271689",
        "modality": "SRA_FASTQ",
        "species": "Homo sapiens",
        "cancer": "NSCLC",
        "ici_related": "yes_first_line_immunotherapy",
        "file_name": "SRA runs",
        "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE271689",
        "local_rel": "",
        "downloaded": "no",
        "skipped_reason": "raw_FASTQ",
        "notes": "Skipped SRA FASTQ as required.",
    },
    {
        "dataset": "NanoString_Hs_WTA_v1.0",
        "accession": "Hs_R_NGS_WTA_v1.0.pkc (Zenodo 12752405; official WTA panel map)",
        "modality": "GeoMx_probe_map",
        "species": "Homo sapiens",
        "cancer": "panel_not_cohort",
        "ici_related": "used_to_interpret_GSE271689_DCC",
        "file_name": "Hs_R_NGS_WTA_v1.0.pkc",
        "url": "https://zenodo.org/records/12752405/files/Hs_R_NGS_WTA_v1.0.pkc?download=1",
        "local_rel": "geomx_pkc/Hs_R_NGS_WTA_v1.0.pkc",
        "downloaded": "yes",
        "skipped_reason": "",
        "notes": "Not a lung-cancer accession. Official Hs WTA v1.0 PKC used only to map RTS IDs. TACSTD2=RTS0027086, CLDN4=RTS0026256.",
    },
    {
        "dataset": "E-MTAB-13530",
        "accession": "E-MTAB-13530",
        "modality": "Visium_WTA_h5",
        "species": "Homo sapiens",
        "cancer": "NSCLC",
        "ici_related": "no_in_SDRF",
        "file_name": "*-filtered_feature_bc_matrix.h5 (40 sections)",
        "url": "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/530/E-MTAB-13530/Files/",
        "local_rel": "emtab13530/",
        "downloaded": "yes",
        "skipped_reason": "",
        "notes": "Processed Space Ranger matrices. TACSTD2 and CLDN4 present. FASTQ/spatial.tar/html skipped.",
    },
    {
        "dataset": "E-MTAB-13530",
        "accession": "E-MTAB-13530",
        "modality": "SDRF",
        "species": "Homo sapiens",
        "cancer": "NSCLC",
        "ici_related": "no_in_SDRF",
        "file_name": "E-MTAB-13530.sdrf.txt",
        "url": "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/530/E-MTAB-13530/Files/E-MTAB-13530.sdrf.txt",
        "local_rel": "emtab13530/E-MTAB-13530.sdrf.txt",
        "downloaded": "yes",
        "skipped_reason": "",
        "notes": "Disease/sampling site/age/sex. No ICI or survival fields.",
    },
    {
        "dataset": "targeted_Xenium-IO_CosMx-1K",
        "accession": "not_used",
        "modality": "targeted_spatial",
        "species": "Homo sapiens",
        "cancer": "lung_cancer_panels",
        "ici_related": "skipped_unless_both_genes",
        "file_name": "",
        "url": "",
        "local_rel": "",
        "downloaded": "no",
        "skipped_reason": "targeted_panel_not_WTA",
        "notes": "Skipped as instructed. Did not invent CosMx/Xenium accessions.",
    },
]


def file_size(rel: str) -> str:
    if not rel:
        return ""
    p = DATA / rel
    if p.is_file():
        return str(p.stat().st_size)
    if p.is_dir():
        return str(sum(x.stat().st_size for x in p.rglob("*") if x.is_file()))
    return ""


def main() -> None:
    NOTES.mkdir(parents=True, exist_ok=True)
    presence = {}
    pride = RESULTS / "pride_gene_presence.tsv"
    if pride.exists():
        pr = pd.read_csv(pride, sep="\t")
        presence["pride"] = pr.to_dict(orient="records")
    miss = RESULTS / "cptac_missingness.tsv"
    if miss.exists():
        presence["cptac_missingness"] = pd.read_csv(miss, sep="\t").to_dict(orient="records")
    spatial = RESULTS / "spatial_run_summary.json"
    if spatial.exists():
        presence["spatial"] = json.loads(spatial.read_text())

    rows = []
    for r in ROWS:
        rec = dict(r)
        rec["file_size_bytes"] = file_size(r.get("local_rel", ""))
        # fill gene presence from analyses when available
        rec["contains_TACSTD2"] = ""
        rec["contains_CLDN4"] = ""
        if r["dataset"].startswith("CPTAC") and r["modality"] == "TMT_proteomics_gene":
            rec["contains_TACSTD2"] = "yes_complete"
            rec["contains_CLDN4"] = "yes_with_NAs"
        if r["dataset"] == "PXD042091" and r["file_name"] == "NSClibrary.txt":
            rec["contains_TACSTD2"] = "yes_library_only"
            rec["contains_CLDN4"] = "yes_library_only"
        if r["dataset"] == "PXD042091" and r["file_name"].endswith(".txt") and "library" not in r["file_name"]:
            rec["contains_TACSTD2"] = "no_not_quantified"
            rec["contains_CLDN4"] = "no_not_quantified"
        if r["dataset"] == "PXD059688" and r["file_name"].endswith("mzTab.gz"):
            rec["contains_TACSTD2"] = "no_in_processed_mzTab"
            rec["contains_CLDN4"] = "no_in_processed_mzTab"
        if r["dataset"] == "GSE271689" and r["modality"] == "GeoMx_WTA_DCC":
            rec["contains_TACSTD2"] = "yes_WTA_RTS0027086"
            rec["contains_CLDN4"] = "yes_WTA_RTS0026256"
        if r["dataset"] == "E-MTAB-13530" and r["modality"] == "Visium_WTA_h5":
            rec["contains_TACSTD2"] = "yes"
            rec["contains_CLDN4"] = "yes"
        rec.pop("local_rel", None)
        rows.append(rec)

    cols = [
        "dataset",
        "accession",
        "modality",
        "species",
        "cancer",
        "ici_related",
        "file_name",
        "url",
        "file_size_bytes",
        "downloaded",
        "skipped_reason",
        "contains_TACSTD2",
        "contains_CLDN4",
        "notes",
    ]
    df = pd.DataFrame(rows)[cols]
    out = NOTES / "protein_spatial_catalog.tsv"
    df.to_csv(out, sep="\t", index=False)
    print("wrote", out, "n=", len(df))


if __name__ == "__main__":
    main()
