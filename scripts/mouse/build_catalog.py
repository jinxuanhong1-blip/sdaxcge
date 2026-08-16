#!/usr/bin/env python3
"""Emit the mouse lung ICI dataset catalog (JSON + Markdown) from vetted metadata.

Every field is grounded in the real records fetched into notes/mouse/raw_meta and
the file listings verified during triage. Accessions are only those provided in the
task brief; none are invented.
"""
import json
import os

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
NOTES = os.path.join(ROOT, "notes", "mouse")

CATALOG = [
    {
        "accession": "GSE239485", "source": "GEO", "species": "Mus musculus",
        "model": "Lewis Lung Carcinoma (LLC), subcutaneous", "tissue": "lung tumor",
        "assay": "bulk RNA-seq", "ici": "anti-PD-1 (+ poly I:C +/- anti-C5aR1)",
        "groups": "Control vehicle (n=8); Poly I:C+anti-PD-1 (n=8); Poly I:C+anti-PD-1+anti-C5aR1 (n=8)",
        "processed_file": "GSE239485_Processed_data.xlsx (7.1 MB, normalized log2, gene_name keyed)",
        "size_ok": True, "decision": "ANALYZED (core, replicated)",
    },
    {
        "accession": "GSE133604", "source": "GEO", "species": "Mus musculus",
        "model": "KrasG12D;p53-/- (KP) lung tumor", "tissue": "lung tumor",
        "assay": "scRNA-seq (10x)", "ici": "anti-PD-1 (+ Asf1a KO)",
        "groups": "Ctrl; Ctrl+anti-PD-1; Asf1a-KO; Asf1a-KO+anti-PD-1 (1 sample/arm, ~6k cells each)",
        "processed_file": "GSE133604_RAW.tar -> per-sample matrix.mtx.gz + shared genes/barcodes (78 MB)",
        "size_ok": True, "decision": "ANALYZED (scRNA pseudobulk + single-cell)",
    },
    {
        "accession": "GSE222158", "source": "GEO", "species": "Mus musculus",
        "model": "Oncogene-driven NSCLC (murine)", "tissue": "lung tumor, FACS-sorted CD45+/CD3+ immune cells",
        "assay": "scRNA-seq (10x)", "ici": "anti-PD-1 (+ DC-CCL21 in situ vaccine)",
        "groups": "CD45: Ctrl/PD1/DCvax/Combo; CD3: Ctrl/Combo (1 sample/arm)",
        "processed_file": "GSE222158_RAW.tar -> per-sample barcodes/features/matrix (162 MB)",
        "size_ok": True, "decision": "ANALYZED (immune-compartment control; tumor markers ~absent)",
    },
    {
        "accession": "GSE129297", "source": "GEO", "species": "Mus musculus",
        "model": "SCLC (RPM) lung tumor", "tissue": "lung tumor",
        "assay": "scRNA-seq (10x)", "ici": "anti-PD-1 (+ CDK7i YKL-5-124)",
        "groups": "Ctrl; anti-PD-1; CDK7i(YKL); anti-PD-1+CDK7i (1 sample/arm, unfiltered droplet matrices)",
        "processed_file": "GSE129297_RAW.tar -> per-sample matrix.mtx.gz + shared features/barcodes (156 MB)",
        "size_ok": True, "decision": "ANALYZED (scRNA pseudobulk + single-cell)",
    },
    {
        "accession": "GSE241978", "source": "GEO", "species": "Mus musculus",
        "model": "CMT167 LUAD (C57BL/6 syngeneic) cell lines", "tissue": "lung adenocarcinoma cells",
        "assay": "bulk RNA-seq (DE summary table)", "ici": "immune-checkpoint pathway (AhR KO -> PD-L1/IDO1); no ICI drug",
        "groups": "AhR-KO vs Cas9 control (DE statistics table, no per-sample matrix)",
        "processed_file": "GSE241978_...CMT_KO_vs_Cas9Ctrl.xlsx (14 MB, DE table)",
        "size_ok": True, "decision": "SUPPORTING (genetic checkpoint-pathway perturbation, not ICI treatment)",
    },
    {
        "accession": "GSE330658", "source": "GEO", "species": "Mus musculus",
        "model": "Egfr-mutant lung cancer (C57BL/6J subcutaneous)", "tissue": "lung tumor",
        "assay": "bulk RNA-seq (per-sample TPM tables)", "ici": "anti-PD-L1 context; deposited arms = PTX / anti-VEGF",
        "groups": "Control; PTX; anti-VEGF; PTX+anti-VEGF (n=2 each). anti-PD-L1 arm NOT in processed deposit",
        "processed_file": "GSE330658_RAW.tar -> 8 per-sample xlsx with TPM (41 MB)",
        "size_ok": True, "decision": "ANALYZED (descriptive, n=2/arm)",
    },
    {
        "accession": "GSE197260", "source": "GEO", "species": "Mus musculus",
        "model": "Egfr-mutant NSCLC (syngeneic C57BL/6J)", "tissue": "subcutaneous lung tumor",
        "assay": "bulk RNA-seq (TPM)", "ici": "anti-PD-1 (4H2) + anti-VEGFR2 (DC101) after EGFR-TKI (gefitinib)",
        "groups": "vehicle_d3, gef_d3, gef_d14, gef_vehicle_d21, gef_dc101_d21, gef_4h2_d21, gef_comb_d21 (n=1/arm)",
        "processed_file": "GSE197260_RNAseqTPM_MM_EGFR-TKI-CD8.txt.gz (0.3 MB, symbol-keyed TPM)",
        "size_ok": True, "decision": "ANALYZED (descriptive trajectory, n=1/arm)",
    },
    {
        "accession": "E-MTAB-13704", "source": "ArrayExpress/BioStudies", "species": "Mus musculus",
        "model": "Lung GEMM (LUGEMM14)", "tissue": "lung tumor in situ",
        "assay": "bulk RNA-seq", "ici": "anti-PD-L1 mono + combos (ATRi, VEGFRi, Cisplatin/aPD-L1/aCTLA4)",
        "groups": "vehicle(5); aPD-L1(5); ATRi(4); ATRi/aPD-L1(5); VEGFRi/aPD-L1(4); Cisplatin/aPD-L1/aCTLA4(4)",
        "processed_file": "GEMMS_raw_counts.csv (14 MB, Ensembl-keyed counts) + sdrf",
        "size_ok": True, "decision": "ANALYZED (core, replicated, richest ICI design)",
    },
    {
        "accession": "GSE297630", "source": "GEO", "species": "Mus musculus",
        "model": "LLC subcutaneous", "tissue": "lung tumor cells (recovered post-treatment)",
        "assay": "microarray (Clariom S Mouse)", "ici": "anti-PD-1 (tolerant/surviving cells vs control)",
        "groups": "Control (C, n=3); anti-PD-1 tolerant (P, n=3)",
        "processed_file": "GSE297630_processed_data.xlsx (5.3 MB, per-sample log2 RMA + author stats)",
        "size_ok": True, "decision": "ANALYZED (core, replicated)",
    },
    {
        "accession": "GSE297632", "source": "GEO", "species": "Mus musculus",
        "model": "LLC subcutaneous", "tissue": "lung tumor (whole)",
        "assay": "scRNA-seq (10x)", "ici": "anti-PD-1 (tolerant vs control)",
        "groups": "Control (~8.4k cells); anti-PD-1 (~9.4k cells)",
        "processed_file": "GSE297632_RAW.tar -> per-sample barcodes/features/matrix (232 MB)",
        "size_ok": True, "decision": "ANALYZED (scRNA pseudobulk + single-cell; companion of GSE297630)",
    },
    {
        "accession": "PXD059688", "source": "PRIDE", "species": "Mus musculus",
        "model": "NSCLC syngeneic (C57BL/6J)", "tissue": "lung tumor",
        "assay": "LC-MS/MS proteomics", "ici": "anti-PD-1 +/- high-dose ascorbic acid (Con/P/AA/AP)",
        "groups": "Con, P(aPD1), AA(ascorbic acid), AP(AA+aPD1) x ~3 (raw MS only)",
        "processed_file": "Only raw .raw/.mgf/.msf (>2GB each) + one single-condition 93-protein mzTab (0.3 MB, no Tacstd2/Cldn4)",
        "size_ok": False, "decision": "CATALOGED, NOT ANALYZED (no <2GB cross-condition protein-quant table; raw MS skipped per rules)",
    },
]


def main():
    with open(os.path.join(NOTES, "catalog.json"), "w") as f:
        json.dump(CATALOG, f, indent=2)

    lines = ["# Mouse Lung ICI Dataset Catalog / 小鼠肺癌免疫检查点抑制 (ICI) 数据集目录",
             "",
             "Focus genes: **Tacstd2** (TROP2, ENSMUSG00000051397) and **Cldn4** (ENSMUSG00000047501).",
             "Scope: mouse (Mus musculus) lung tumor datasets with ICI / immune context.",
             "Rule: processed matrices only; skip FASTQ / raw MS; skip files > 2 GB; no invented accessions.",
             "",
             "| Accession | Source | Model / Tissue | Assay | ICI context | Groups | Processed file(s) | Decision |",
             "|---|---|---|---|---|---|---|---|"]
    for c in CATALOG:
        lines.append("| {accession} | {source} | {model}; {tissue} | {assay} | {ici} | {groups} | {processed_file} | {decision} |".format(**c))
    lines += ["",
              "## Notes",
              "- **ANALYZED (core, replicated)**: E-MTAB-13704, GSE239485, GSE297630 — replicated treated-vs-control designs used for formal statistics.",
              "- **ANALYZED (scRNA)**: GSE129297, GSE133604, GSE297632, GSE222158 — pseudobulk + single-cell correlations.",
              "- **ANALYZED (descriptive)**: GSE330658 (n=2/arm), GSE197260 (n=1/arm) — reported as fold changes/trends, underpowered for p-values.",
              "- **SUPPORTING**: GSE241978 — AhR-KO (genetic checkpoint-pathway perturbation), DE table only.",
              "- **NOT ANALYZED**: PXD059688 — all quantitative MS files are raw and >2 GB; the only <2 GB processed file is a single-condition 93-protein identification report that does not contain Tacstd2/Cldn4.",
              "",
              "Provenance: metadata in `notes/mouse/raw_meta/`; download URLs + sizes + checksums in `notes/mouse/data/MANIFEST.tsv`; gene ID map in `notes/mouse/gene_map.json`."]
    with open(os.path.join(NOTES, "catalog.md"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print("wrote catalog.json and catalog.md;", len(CATALOG), "datasets")


if __name__ == "__main__":
    main()
