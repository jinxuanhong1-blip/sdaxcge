# ENCODE / Roadmap: CLDN4 chromatin versus DNA-repair genes

Exploratory public-data check. It does not say that CLDN4 regulates DNA repair, and it does not touch the locked CosMx, concordant-4, GSE137244, TCGA, or TISMO results.

The only KEGG DNA-repair gene near CLDN4 is **RFC2**, 390,540 bp downstream on 7q11.23 (TSS-to-TSS 454,916 bp). Baseline epithelial Hi-C loop calls do not connect the two genes. Contact domains usually keep them apart. A subset of CTCF and POLR2A ChIA-PET experiments do link the 3' end of CLDN4 to the RFC2 gene body. Roadmap chromatin state marks the CLDN4 TSS as an enhancer or quiescent, and the RFC2 TSS as an active promoter.

## PR inventory

Searched this repository on 2026-09-21 (open PRs through #682). No PR runs ENCODE or Roadmap Hi-C, ChIA-PET, ChIP, or chromHMM at CLDN4. Parallel PRs ask a different question, RNA rather than chromatin, and are not re-cut here: #664 (CLDN4 NHEJ paper omics hunt), #665 (DepMap lung CLDN4 versus DNA-damage signatures), #675 (concordant-4 CLDN4 versus NHEJ / STING / MHC-I), and #678–#681 (Ku, DNA-PKcs, and NHEJ RNA sweeps). The closest older ChIP PR is #191, which is an ELF3 binding hunt, not this locus.

## Definitions

- Assembly: GRCh38. Ensembl gene coordinates are converted to 0-based half-open intervals. CLDN4 is chr7:73,799,541–73,832,690, TSS 73,799,541 (`ENSG00000189143`).
- DNA-repair set: KEGG hsa03410, hsa03420, hsa03430, hsa03440, hsa03450, and hsa03460. 176 symbols; 175 placed on a primary chromosome.
- Near, linear: gap between gene bodies ≤ 1 Mb. A 2 Mb gap is also reported.
- Hi-C / ChIA-PET loop hit: one anchor overlaps CLDN4 ± 5 kb and the other overlaps a repair gene ± 5 kb. Promoter–promoter is the stricter TSS ± 2 kb test.
- Same contact domain: one released domain overlaps both gene bodies and is ≤ 2 Mb.
- Epithelial Hi-C / ChIA-PET cohort: released ENCODE experiments whose biosample organ slim contains epithelium, plus A549 (lung carcinoma; ENCODE does not file it under the epithelium slim). Endothelial, mesothelial, melanoma, SK-N-MC, and HEK293 rows were dropped.
- Baseline: no auxin, degron, IAA treatment, knockout, or dexamethasone. A dexamethasone 0-hour label stays in baseline. Tagged-but-untreated degron lines are a sensitivity set, not the primary count.
- ChIP: Roadmap core 15-state chromHMM, hg38 lift, at the TSS. ENCODE adds one baseline peak file per biosample and target (TSS ± 2 kb). A shared peak has to overlap both promoters and be ≤ 50 kb.

Reproduce with `python3 methods/encode_roadmap_cldn4_dnarepair/analyze.py`.

## Linear neighborhood

RFC2 is the only KEGG repair gene within 2 Mb. It is in base-excision repair, nucleotide-excision repair, and mismatch repair. The next gene on chromosome 7 is SEM1, 22.6 Mb away. Eleven repair genes sit on chromosome 7; the other ten are 22.6–78.8 Mb from CLDN4.

Protein-coding genes between CLDN4 and RFC2 are METTL27, TMEM270, ELN, LIMK1, EIF4H, and LAT2. None of those are in the KEGG repair set. RFC2 is inside the Williams–Beuren interval with CLDN3 and CLDN4. Being in that interval is not the same as sharing a loop or a chromatin state.

## Hi-C

Thirty baseline epithelial Hi-C experiments have GRCh38 contact-domain calls (8 intact, 20 in situ, 2 dilution). **One** places CLDN4 and RFC2 in the same domain: HCT116 in situ Hi-C `ENCSR637QCS`, two calls of 1.01 Mb and 1.25 Mb (chr7:73,495,000–74,505,000 and chr7:73,420,000–74,670,000). The other 29 do not. Other HCT116 maps keep CLDN4 in a ~120 kb domain around 73.74–73.86 Mb and RFC2 in a separate domain around 74.20–74.42 Mb.

Lung maps do not put them together:

- A549 parental in situ Hi-C `ENCSR444WCZ` has no called domain over the CLDN4 gene body. The neighboring domain ends at 73,785,670, about 14 kb upstream.
- A549 dexamethasone 0-hour in situ Hi-C `ENCSR662QKG` puts CLDN4 in chr7:73,740,000–73,860,000 and RFC2 in chr7:74,205,000–74,400,000.
- NCI-H460 in situ Hi-C `ENCSR489OCU` has no called domain over CLDN4. That experiment has no released loop file.

Twenty-seven baseline epithelial Hi-C experiments have loop calls. **None** link CLDN4 ± 5 kb to RFC2 or to any other KEGG repair gene. Seventeen of those experiments have at least one loop that touches CLDN4. In every case the other anchor is still at least 353,230 bp from RFC2 (median closest gap 367,230 bp). Those partners sit just downstream of CLDN4, around 73.85–73.86 Mb. They are local CLDN4 loops, not loops that reach RFC2. Zero of the 17 are within 200 kb of RFC2.

The A549 0-hour loop file (`ENCFF803ZOW`) contains no loop with an anchor on CLDN4 ± 5 kb. Parental A549 Hi-C has contact matrices and domains, not a released loop file. Tagged-untreated HCT116 degron Hi-C (8 experiments) also has no repair-gene loop hit. That sensitivity set is not part of the 0/27 count.

## ChIA-PET

ChIA-PET is the ChIP-tethered contact assay, reported separately from Hi-C. Eighteen baseline epithelial experiments were scanned. **Five** link CLDN4 ± 5 kb to RFC2, and none of the five is promoter-to-promoter (TSS ± 2 kb at both ends). Anchor coordinates below are 0-based.

| Biosample | Target | Experiment | CLDN4-side anchor | RFC2-side anchor | Overlaps CLDN4 gene body |
|---|---|---|---|---|---|
| A549 | CTCF | ENCSR911ZMB | 73,814,320–73,814,863 | 74,252,034–74,252,626 | yes |
| A549 | POLR2A | ENCSR138NSW | 73,833,921–73,834,572 | 74,249,921–74,250,526 | no; starts 1,231 bp downstream |
| MCF-7 | POLR2A | ENCSR059HDE | 73,829,526–73,830,127 | 74,258,137–74,258,777 | yes |
| Panc1 | CTCF | ENCSR145PYF | 73,821,850–73,822,376 | 74,221,851–74,222,503 | yes |
| RWPE1 | CTCF | ENCSR030KAB | 73,804,640–73,805,279 | 74,259,032–74,259,656 | yes |

By target, that is CTCF in 3/8 baseline experiments, POLR2A in 2/8, and RAD21 in 0/2. Anchors are about 0.5–0.7 kb, so these are focal contacts across the ~0.39 Mb gap, not a single broad peak. Four additional HCT116 tagged-untreated ChIA-PET experiments (CTCF or POLR2A) also link the same pair. Those are outside the 5/18 baseline count. HepG2, MCF10A, HMEC, and parental HCT116 baseline ChIA-PET do not.

## ChIP and chromHMM

Roadmap epithelial **cell and carcinoma lines** (9 epigenomes with an hg38-lifted 15-state model):

| EID | Sample | CLDN4 TSS | RFC2 TSS |
|---|---|---|---|
| E027 | Breast myoepithelial | Quies | TssA |
| E028 | vHMEC | Enh | TssA |
| E057 | Foreskin keratinocyte skin02 | Enh | TssBiv |
| E058 | Foreskin keratinocyte skin03 | Enh | TssA |
| E114 | A549 lung carcinoma | Enh | TssA |
| E117 | HeLa-S3 | Quies | TssA |
| E118 | HepG2 | Enh | TssA |
| E119 | HMEC | Enh | TssA |
| E127 | NHEK | Enh | TssA |

CLDN4 TSS: Enh in 7/9, Quies in 2/9 (E027, E117). RFC2 TSS: TssA in 8/9, TssBiv in 1/9 (E057). Mucosal tissues (8) keep RFC2 at TssA. CLDN4 there is Quies (colon, esophagus, gastric), Enh (rectal mucosa, small intestine), TxWk (sigmoid colon), or ReprPC (stomach mucosa). Bulk lung E096, which is not an epithelial isolate, is Quies at CLDN4 and TssA at RFC2.

ENCODE peak files, one baseline experiment per biosample and target where one exists (50 files scored): H3K4me3 overlaps the RFC2 promoter in 8/8 files and the CLDN4 promoter in 1/8 (HCT116 only). H3K27ac overlaps RFC2 in 8/8 and CLDN4 in 2/8 (HCT116, Panc1). A549 matches the chromHMM call: H3K4me3 and H3K27ac at RFC2, H3K4me1 but not H3K4me3/H3K27ac at the CLDN4 TSS. Across those 50 files, no peak of ≤ 50 kb overlaps both promoters.

NCI-H460 has Hi-C and no released histone or CTCF ChIP in this baseline search. MCF10A has CTCF ChIP and Hi-C, and no baseline histone ChIP under the names used here.

## What this does not say

RFC2 sharing the Williams–Beuren neighborhood with CLDN4 is not evidence of a DNA-repair program at CLDN4. Hi-C loops do not join them in the baseline epithelial maps surveyed here. The ChIA-PET links that do exist are 3'-end contacts in some cell lines, not a promoter connection, and they are absent in other epithelial ChIA-PET experiments. Chromatin state at the two TSS positions does not match.

## Files

- `results/encode_roadmap_cldn4_dnarepair/kegg_dna_repair_distance_to_cldn4.tsv`
- `results/encode_roadmap_cldn4_dnarepair/experiment_summary.tsv`
- `results/encode_roadmap_cldn4_dnarepair/loop_links_to_dna_repair.tsv`
- `results/encode_roadmap_cldn4_dnarepair/domain_overlaps_dna_repair.tsv`
- `results/encode_roadmap_cldn4_dnarepair/hic_cldn4_loop_distance_to_rfc2.tsv`
- `results/encode_roadmap_cldn4_dnarepair/roadmap_chromhmm_tss.tsv`
- `results/encode_roadmap_cldn4_dnarepair/encode_chip_tss_peaks.tsv`
- `results/encode_roadmap_cldn4_dnarepair/fig_locus_domains.png`
- `results/encode_roadmap_cldn4_dnarepair/fig_chromhmm_tss.png`
