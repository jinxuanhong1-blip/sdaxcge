"""Curated inclusion/exclusion decisions for this slice.

Every series listed here is either (a) named in the task, (b) an older
array/panel ICI cohort that survived sample-level screening, or (c) a
near-miss that must be documented so absence is auditable.

`decision` values:
  analyze          - genome-wide / DASL array, target genes present, response labels present
  document_absent  - ICI + response (or named series) but platform cannot measure the targets
  document_no_resp - target genes present but GEO has no per-sample response label
  exclude          - out of scope (RNA-seq, xenograft, no ICI, non-human, etc.)
"""

from __future__ import annotations

CATALOG = [
    {
        "gse": "GSE93157",
        "decision": "document_absent",
        "reason": "Named series. NanoString PanCancer Immune 730-gene panel (GPL19965). Per-sample RECIST (best.resp) and PFS are present for NSCLC / HNSCC / melanoma anti-PD-1, but TACSTD2 and CLDN4 are not on the panel.",
        "tissue": "tumor (NSCLC / HNSCC / melanoma)",
        "platform": "GPL19965",
        "analyze_controls": True,
        "response_fields": ["best.resp"],
        "or_labels": {"CR": 1, "PR": 1, "SD": 0, "PD": 0},
        "dcb_labels": {"CR": 1, "PR": 1, "SD": 1, "PD": 0},
        "subset_lung": True,
        "lung_field": "source_name",
        "lung_pattern": r"LUNG",
    },
    {
        "gse": "GSE202417",
        "decision": "analyze",
        "reason": "Affymetrix Clariom D (GPL23126) measures TACSTD2 and CLDN4. Peripheral-blood CD8+ T cells from advanced NSCLC treated with nivolumab + bezafibrate; phenotype R / NR is annotated. Primary contrast is pre-treatment R vs NR.",
        "tissue": "blood CD8+ T cells (NSCLC)",
        "platform": "GPL23126",
        "analyze_controls": True,
        "response_fields": ["phenotype"],
        "or_labels": {"R": 1, "NR": 0},
        "dcb_labels": {"R": 1, "NR": 0},
        "require": {"time": r"(?i)pre"},
    },
    {
        "gse": "GSE67501",
        "decision": "analyze",
        "reason": "Classic older ICI microarray (Illumina HT-12 WG-DASL). RCC (not lung) anti-PD-1 with RECIST + binary response. Included because it is one of the few public array cohorts that both measures TACSTD2/CLDN4 and reports ICI response. n=11; interpret as exploratory.",
        "tissue": "RCC tumor (primary / metastasis)",
        "platform": "GPL14951",
        "analyze_controls": True,
        "response_fields": [
            "response to anti-pd-1 (nivolumab) immunotherapy (response or no-response)",
        ],
        "or_labels": {"response": 1, "no_response": 0},
        "dcb_labels": {"response": 1, "no_response": 0},
    },
    {
        "gse": "GSE140901",
        "decision": "document_absent",
        "reason": "Same NanoString 730-gene immune panel as GSE93157 (HCC, not lung). Response labels exist; TACSTD2/CLDN4 are absent from the panel.",
        "tissue": "HCC tumor",
        "platform": "GPL19965",
        "analyze_controls": True,
        "response_fields": ["best_response"],
        "or_labels": {"CR": 1, "PR": 1, "SD": 0, "PD": 0},
        "dcb_labels": {"CR": 1, "PR": 1, "SD": 1, "PD": 0},
    },
    {
        "gse": "GSE248249",
        "decision": "document_no_resp",
        "reason": "Clariom D NSCLC FFPE pre/post immunotherapy; TACSTD2/CLDN4 are on the array, but GEO sample records do not carry RECIST or responder labels.",
        "tissue": "NSCLC FFPE tumor",
        "platform": "GPL23126",
        "analyze_controls": False,
    },
    {
        "gse": "GSE141479",
        "decision": "document_no_resp",
        "reason": "Clariom D blood CD8+ from NSCLC patients on nivolumab; TACSTD2/CLDN4 are on the array, but GEO has no response annotation (Havel et al. 2020; labels live in the paper, not the SOFT).",
        "tissue": "blood CD8+ T cells (NSCLC)",
        "platform": "GPL23126",
        "analyze_controls": False,
    },
    {
        "gse": "GSE305086",
        "decision": "document_no_resp",
        "reason": "HG-U133 Plus 2.0 whole-blood under immunotherapy (1LIO / 2LIO / CHTIO). TACSTD2/CLDN4 probes exist, but GEO annotates treatment line / timepoint only — no RECIST.",
        "tissue": "whole blood",
        "platform": "GPL570",
        "analyze_controls": False,
    },
    {
        "gse": "GSE99070",
        "decision": "document_no_resp",
        "reason": "Illumina HT-12 v4 malignant pleural mesothelioma; some samples treated with nivolumab/pembrolizumab vs treatment-naive. TACSTD2/CLDN4 present. No RECIST. Not lung NSCLC.",
        "tissue": "mesothelioma tumor",
        "platform": "GPL10558",
        "analyze_controls": False,
    },
    {
        "gse": "GSE180347",
        "decision": "exclude",
        "reason": "NanoString Cancer Immune panel on PD-L1+/- 'hot' LUAD. No ICI treatment, no response. Panel does not list TACSTD2/CLDN4.",
        "tissue": "LUAD tumor",
        "platform": "GPL29738",
        "analyze_controls": False,
    },
    {
        "gse": "GSE190731",
        "decision": "exclude",
        "reason": "Durvalumab/oleclumab-treated NSCLC xenografts (cell-line), not patient tumors.",
        "tissue": "xenograft",
        "platform": "GPL570",
        "analyze_controls": False,
    },
    {
        "gse": "GSE261345",
        "decision": "exclude",
        "reason": "ES-SCLC chemo-immunotherapy with RECIST, but platform GPL30173 is NextSeq 2000 (GeoMx DSP sequencing), not a microarray.",
        "tissue": "ES-SCLC (spatial)",
        "platform": "GPL30173",
        "analyze_controls": False,
    },
    {
        "gse": "GSE261348",
        "decision": "exclude",
        "reason": "Sister series to GSE261345 (atezolizumab + platinum/etoposide). Sequencing, not microarray.",
        "tissue": "ES-SCLC (spatial)",
        "platform": "GPL30173",
        "analyze_controls": False,
    },
    {
        "gse": "GSE135222",
        "decision": "exclude",
        "reason": "NSCLC anti-PD-1 with PFS, but RNA-seq (GPL16791), out of this microarray slice.",
        "tissue": "NSCLC tumor",
        "platform": "GPL16791",
        "analyze_controls": False,
    },
    {
        "gse": "GSE126044",
        "decision": "exclude",
        "reason": "NSCLC anti-PD-1, RNA-seq. Out of slice.",
        "tissue": "NSCLC tumor",
        "platform": "GPL16791",
        "analyze_controls": False,
    },
    {
        "gse": "GSE136961",
        "decision": "exclude",
        "reason": "NSCLC anti-PD-1 durable-benefit signatures, RNA-seq. Out of slice.",
        "tissue": "NSCLC tumor",
        "platform": "GPL24014",
        "analyze_controls": False,
    },
]
