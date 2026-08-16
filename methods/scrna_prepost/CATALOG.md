# Catalog: public human lung-tumor scRNA with pre and post ICI / chemo-IO

Inclusion: human lung **tumor** single-cell RNA with **both** a pre-ICI (or chemo-IO) timepoint and a post timepoint, even if patients are unmatched. Malignant TACSTD2 / CLDN4 is the endpoint. Blood-only, T-sorted-only, post-only, pre-only, and controlled-access matrices are catalogued and excluded from the pool.

GEO searches (2026-08-16): lung/NSCLC × scRNA × ICI/PD-1/neoadjuvant (104 series); lung × scRNA × pre/post/paired (112 series); title Pre/Post (41 series). Must-try accessions below were opened on GEO (SOFT + supplementary file lists).

## Testable for malignant pre vs post

| Accession | Year | Design | Pairing | Tumor cells? | Pre n | Post n | Decision |
|---|---|---|---|---|---|---|---|
| **GSE207422** | 2023 | Neoadjuvant PD-1 + chemo, Hu et al. *Genome Med* | Unmatched (different patients) | Yes (whole digest) | 3 biopsies (P01, P05, P08) | 12 resections (4 MPR/pCR, 8 NMPR) | **Analyze.** Sample-unit MWU. |

## Must-try: opened, not a pre+post tumor-cell test

| Accession | What GEO actually deposits | Why not in the pool |
|---|---|---|
| **GSE337519** | Series text: one patient, paired pre + post after 3 cycles chemo-IO. Supplementary: **one** 10x library, GSM9856929 titled “luad1 / one sample from non small cell lung cancer”. 8,612 barcodes, no HTO / sample tag, no timepoint column. | Cannot assign cells to pre vs post. Wilcoxon not defined. |
| **GSE205335** | 33 samples / 26 patients receiving ICI (eLife 2024). GEO characteristics: patient, tissue, stage, histology, RECIST. **No timepoint field.** Multi-sample patients are same-time different sites (EBUS + effusion / liver / neck), not serial biopsies. EGA EGAD00001008703 is titled pre-treatment. Paper excluded some “acquired after ICI” samples from the core set; those IDs are not in GEO SOFT. | No public pre vs post labels. Same-patient pairs in GEO are multi-site, not serial. |
| **GSE241934** | NEOTIDE/CTONG2104 + real-world neoadjuvant IO+chemo. IIT 11 patients, Real 34 patients. Metadata has cycles, PD-1 agent, MPR — **all resected after neoadjuvant therapy**. No pre column. | Post-only. |
| **GSE146100** | Three surgically resected nodules from **one** patient after pembrolizumab (W1/W3 NR, W2 R). | Post-only; no pre scRNA. |
| **GSE291670** | 6 post-resection samples (3 MPR, 3 NMPR) after anlotinib + camrelizumab. Paper’s 3 pre biopsies are HRA001033 (= GSE207422 TN), not in this GEO. | Post-only in the public accession. |

## 2023–2026 leftovers with Pre/Post or ICI in the record (not both tumor timepoints)

| Accession | Why excluded |
|---|---|
| GSE243013 | 234 **post**-neoadjuvant chemo-IO tumors. No pre. |
| GSE233203 | 7 pleural-effusion PDCs **before** ABCP. Pre-only. |
| GSE229353 | 7 **post**-neoadjuvant CD45+ sorts (1 chemo, 6 PD-1+chemo). No pre; no malignant cells. |
| GSE179994 | Liu et al. *Nat Cancer*: 33 TN + 9 post-R + 5 post-NR biopsies, **T cells only**. No malignant compartment. |
| GSE176021 / GSE185204 / GSE185206 | Neoadjuvant or post-ICI **T-sorted** scRNA/TCR. No malignant TACSTD2/CLDN4. |
| GSE280232 | Neoadjuvant nivo ± ipi / chemo; GEX is sorted T cells or mixed CD45± **after** neoadjuvant. No pre tumor arm. |
| GSE253013 | LUAD stroma atlas. No ICI timepoints. |
| GSE189357 / GSE150938 | Early LUAD / GGN atlases. Not ICI. |
| GSE286228 | PBMC, anti-synthetase ILD, not lung tumor ICI. |
| HRA006493 / PRJCA022724 | Paired pre/post NSCLC chemo-IO or anti-VEGFA. **GSA-Human controlled access.** Not public. |
| PRJNA1068179 | Molecular Cancer 2025: 6 unmatched TN + 6 post chemo-IO. **SRA raw only**; no processed gene-barcode matrix on GEO. |

## Pool

Only **GSE207422** contributes an independent public malignant-cell pre vs post contrast. Combined estimate = that one study (direction + Stouffer on the sample-unit MWU). GSE337519 is catalogued but not pooled.
