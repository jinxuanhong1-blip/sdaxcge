# Ku70/Ku80 KO/KD RNA-seq inventory

Queried 2026-09-21 via NCBI E-utilities (`db=gds`) and the ENCODE portal.
No accession in this note was typed in from memory without a live hit.

## Search queries

1. `("XRCC5" OR "XRCC6" OR "Ku70" OR "Ku80" OR "Ku86") AND (knockout OR knockdown OR CRISPR OR shRNA OR siRNA OR deficient OR deletion) AND "gse"[Entry Type]` — 82 GDS hits
2. Title query for XRCC5/XRCC6/Ku70/Ku80/Ku86 — 65 hits
3. Same gene names AND `Expression profiling by high throughput sequencing` — 78 hits
4. Tighter follow-ups (XRCC5 or XRCC6 with shRNA/CRISPR/knockout/knockdown; Ku70 or Ku80 depletion plus RNA-seq; mouse Xrcc6 knockout) — 33 unique series after summary

Unique series across the broad queries: 125. Titles and summaries were filtered for an actual Ku70/Ku80 perturbation rather than a passing mention.

## Included (processed counts, mammalian, real KO/KD)

| Accession | What it is | Why it is in the meta |
|---|---|---|
| GSE294709 | HCT116 Ku80-AID, HCT116 Ku86-flox CreERT2, HEK293 Ku70-KO with dox-inducible Ku70 rescue. PMID 40373806 | Author gene counts. Several pre-specified contrasts |
| GSE180581 | HEK293T Ku70+/−, Ku80+/− or DNA-PKcs+/− plus siRNA. PMID 34849385 | Author count matrix. DNA-PKcs arm is a comparator, not a Ku contrast |
| GSE247031 | Mouse Treg-specific Xrcc6 knockout in an LLC lung-tumor model | Author raw counts, n=3 vs 3 |
| GSE88488 | ENCODE K562 shXRCC5 (ENCSR715XZS) | Matched to non-targeting ENCSR815CVQ |
| GSE88126 | ENCODE K562 shXRCC6 (ENCSR232CPD) | Same K562 non-targeting control |
| GSE80921 | ENCODE HepG2 shXRCC5 (ENCSR732IYM) | Matched to ENCSR491FOC |
| GSE80894 | ENCODE HepG2 shXRCC6 (ENCSR500WHE) | Same HepG2 non-targeting control |
| GSE176960 | ENCODE K562 CRISPR XRCC5 (ENCSR276GMG) | Matched to ENCSR292PXV |
| GSE177139 | ENCODE K562 CRISPR XRCC6 (ENCSR516EPT) | Matched to ENCSR404OHQ (K562, despite a mismatched free-text description) |
| GSE177121 | ENCODE HepG2 CRISPR XRCC6 (ENCSR312VLS) | Matched to ENCSR964HKT |

ENCODE GEO series contain only the knockdown samples. Controls were taken from the ENCODE `possible_controls` field and quantified from GRCh38 V29 gene quantifications (`expected_count`).

## Reviewed and excluded

| Accession | Why excluded |
|---|---|
| GSE278950 | Foxp3/Ku70 ChIP-seq, not RNA-seq |
| GSE219220 | Ku70/Ku80 PAR-CLIP of bound RNA, not a perturbation transcriptome |
| GSE109026 | Ku-dependent rRNA / DNA-PKcs CLIP-style and related assays, not a Ku KO expression contrast used here |
| GSE26944 | Ku70-as-DNA-sensor study (type III IFN after DNA stimulation), not a Ku KO/KD |
| GSE2498 | Old Ku86/telomerase ablation microarray, outside the RNA-seq scope |
| GSE162453 | Ku70/Lig4 end-joining outcomes, not an expression matrix of Ku loss |
| GSE247031's companion GSE278950 | ChIP, see above |
| Fungal Ku70 deletions (GSE12893, GSE159974, GSE217909 and other Penicillium/Aspergillus hits) | No mammalian IFN/STING/APM program |
| GSE201844, GSE271695, GSE287014, GSE265857, GSE301051, GSE346664 | Ku70/Ku80 mentioned as a pathway or binding partner; the perturbation is a different gene |
| GSE71489 | Bisphenol A and DNA repair, not a Ku knockout |
| GSE120105 / GSE120110, GSE77634, GSE327675, GSE327677 | Binding or GRO-seq, not KO/KD RNA-seq |
| GSE186384 / GSE80196 and related SPAR-seq superseries | Splicing screens; Ku is not an isolated KO/KD expression contrast with a deposited count matrix used here |

GSE294709 also contains IAA day-4 and day-6 samples whose only DMSO baseline is the day-2 DMSO arm. Those time-unmatched comparisons were not tested. Matched contrasts (day-2 IAA vs day-2 DMSO, 24 h, CDKi 48 h, CDKi day 4) were tested. Day-4 CDKi is the primary AID unit; the other AID time points are the same clones and are not independent meta units.

## On-target rule

- siRNA, shRNA, CRISPR, flox excision, dox withdrawal: the perturbed mRNA must fall (log2FC < −0.4, or log2FC < 0 with FDR < 0.05). Failures stay in the table and leave the meta.
- AID degron: protein is degraded. An mRNA drop is reported but not required.
