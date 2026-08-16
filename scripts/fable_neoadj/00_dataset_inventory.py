"""Write the dataset inventory (all GEO series considered + include/exclude reasons)."""
import os
import pandas as pd
from common import TAB

ROWS = [
    ["GSE207422", "bulk RNA-seq + scRNA", "resectable NSCLC, neoadjuvant anti-PD-1 + platinum chemo",
     "pre-treatment biopsy (bulk, n=24) + post-tx scRNA (n=15)", "MPR/pCR/NMPR + residual%", "INCLUDED (primary/baseline + scRNA leftover)",
     "Bulk log2TPM 5.6MB; scRNA UMI 184MB; scRNA mostly post-tx, EPCAM+/PTPRC- gate"],
    ["GSE205335", "scRNA (10x 3'/5')", "26 lung cancer pts on ICI (mostly stage IV; ADC/SQ/SCLC/NUT)",
     "metastatic LN / lung / liver / effusion (n=33 samples)", "RECIST PR/SD/PD/NE (NOT MPR)", "INCLUDED as leftover (honest mismatch)",
     "RDS 524MB; 28,512 published malignant cells; advanced ICI, not neoadjuvant"],
    ["GSE241934", "scRNA/TCR (10x)", "NEOTIDE/CTONG2104 EGFR-mut sintilimab + real-world PD-1 +/- chemo",
     "post-treatment resected tumor (n=45)", "MPR/pCR/non-MPR + response rate", "INCLUDED (validation)",
     "IIT 466MB + Real 1.25GB (<2GB each); Epi compartment annotated"],
    ["GSE146100", "scRNA (10x)", "1 patient, 3 synchronous LUAD nodules, induction pembrolizumab",
     "surgical nodules (n=3)", "responded / non-responded (per nodule)", "INCLUDED (descriptive)",
     "229MB; n=1 patient -> anecdotal, EPCAM+ gating"],
    ["GSE243013", "scRNA/TCR (10x)", "234 NSCLC, post neoadjuvant chemo-immunotherapy",
     "post-treatment", "response labels", "EXCLUDED",
     "Counts matrix 7.1GB > 2GB cap; CD45/immune-heterogeneity atlas"],
    ["GSE229353", "scRNA (10x)", "7 NSCLC, NAC vs neoadjuvant pembrolizumab+chemo",
     "post-treatment", "response", "EXCLUDED",
     "CD45+ immune-sorted -> epithelial ADC targets not captured"],
    ["GSE280232", "scRNA/TCR", "13 KRAS-mutant resectable NSCLC, neoadjuvant ICB",
     "post-treatment", "recurrence / STK11 status", "EXCLUDED",
     "TIL/CD8-focused; KRASmut-only; no epithelial MPR readout"],
    ["GSE248378", "bulk RNA-seq (FFPE)", "phase II durvalumab +/- SBRT, stages I-III NSCLC",
     "post-treatment resected tumor (n=29)", "recurrence (no MPR/pCR)", "EXCLUDED",
     "Endpoint is freedom-from-recurrence; no MPR/pCR label in GEO"],
    ["GSE225620", "bulk RNA-seq", "31 NSCLC, neoadjuvant tislelizumab + chemo",
     "whole blood", "responder/non-responder", "EXCLUDED",
     "Whole blood -> epithelial TACSTD2/CLDN4 not informative"],
    ["GSE176021/176022", "scRNA + TCR-seq", "neoadjuvant anti-PD-1 NSCLC (neoantigen-specific T cells)",
     "T cells", "response", "EXCLUDED",
     "Neoantigen-specific T-cell focus; not tumor epithelium"],
]
COLS = ["accession", "assay", "cohort", "sample_type", "response_annotation", "decision", "notes"]

df = pd.DataFrame(ROWS, columns=COLS)
os.makedirs(TAB, exist_ok=True)
df.to_csv(os.path.join(TAB, "dataset_inventory.csv"), index=False)
print(df.to_string())
