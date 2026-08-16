# C2 public literature — CLDN4 lysosomal degradation after ADC

**Scope.** Catalog only. C2 remains private. This folder lists public papers that do or do not document **CLDN4 protein lysosomal degradation after an ADC**. Private claim text is not copied or reconstructed.

**Date.** 2026-08-16. Nothing invented: every row has a DOI/PMID that was opened or API-confirmed, or is explicitly marked NONE.

## Found vs none

| # | Analog class | Verdict | Public analog (if any) |
| --- | --- | --- | --- |
| 1 | CLDN4 protein is lysosomally degraded after an ADC | **NONE** | No paper has ADC on + CLDN4 protein down + lysosome dependence |
| 2 | TROP2-ADC (SG / Dato-DXd / SKB264) causes CLDN4 lysosomal loss | **NONE** | Those ADCs traffic **TROP2** to lysosomes; CLDN4 is not read out |
| 3 | Published CLDN4-targeted ADC with CLDN4 protein-fate data | **NONE** | No peer-reviewed CLDN4-ADC. CLDN6/18.2 ADCs are the wrong claudin |

Do not treat near-misses as hits. They fail the analog as written.

## What would have counted as FOUND

A public paper that does all of:

1. Treats cells or tumors with an ADC (TROP2-ADC or anti-CLDN4 ADC).
2. Shows CLDN4 protein decrease, or CLDN4 in LAMP1+ / lysosomal compartments.
3. Shows lysosome dependence (chloroquine, bafilomycin, or equivalent rescue / colocalization).

No such paper was found.

## Near-misses (catalogued, not hits)

### A. TROP2 / claudin lysosome — wrong claudin, not ADC

**Wu et al., *Cells* 2020** — DOI `10.3390/cells9041027`, PMID 32326212, PMC7226414.

Matriptase cleaves TROP2 and EpCAM in keratinocytes. Associated claudins fall. Chloroquine rescues **claudin-1 and claudin-7** and accumulates cleaved TROP2/EpCAM. This is the closest public TROP2 → lysosomal claudin-loss experiment. It is not CLDN4 and not an ADC.

**Wu et al., *J Clin Invest* 2017** — DOI `10.1172/jci88428`, PMID 28094766.

EpCAM cleavage by matriptase destabilizes intestinal claudins (mainly claudin-7). Same pathway family; EpCAM not TROP2-ADC; not CLDN4.

**Zhao / Wu et al., *JITC* 2026** — DOI `10.1136/jitc-2025-012265`, PMID 41932810.

TROP2 KO and naked **hRS7** (the SG antibody, not the ADC) lower **claudin-7** / occludin by IF/IHC in TNBC and increase T-cell infiltration. Wrong claudin; not ADC; no lysosome assay.

### B. TROP2-ADC lysosome — right drug class, wrong antigen readout

**Okajima et al., *Mol Cancer Ther* 2021** — DOI `10.1158/1535-7163.mct-21-0206`, PMID 34413126.

Dato-DXd binds TROP2, internalizes, and releases DXd after lysosomal protease cleavage. No CLDN4 blot or IF.

**Cardillo 2011** (`10.1158/1078-0432.ccr-10-2939`) and **Goldenberg 2015** (`10.18632/oncotarget.4318`): SG/IMMU-132 Trop-2 binding, internalization, SN-38 delivery. No CLDN4 protein fate.

**ASCO 2026 e13006** (`10.1200/jco.2026.44.16_suppl.e13006`): paired pre/post-SG MS of Trop2, TOPO1, HER2/3, Nectin4, SLFN11. No TJ/claudin panel.

### C. TROP2–CLDN4 proximity — right pair, no degradation, no ADC

**Nakatsukasa et al., *Am J Pathol* 2010** — DOI `10.2353/ajpath.2010.100149`, PMID 20651236.

PLA: TACSTD2 is next to CLDN4 (and CLDN1/7/TJP1) on cornea. Co-IP of CLDN4 is **negative**.

**Van Itallie / Fredriksson et al., *PLOS ONE* 2015** — DOI `10.1371/journal.pone.0117074`, PMID 25789658.

CLDN4-BioID in MDCK II lists TROP-2 among Cldn4-neighborhood proteins (Table 2). Proximity MS, not treatment.

### D. CLDN4 can be endocytosed and degraded — right protein, no ADC

These show that CLDN4 is not a static junction protein. They do not show ADC-triggered lysosomal CLDN4 loss.

| Paper | DOI | PMID | Stimulus | Fate evidence |
| --- | --- | --- | --- | --- |
| Cong 2015 *J Cell Sci* | 10.1242/jcs.165878 | 25948584 | carbachol / mAChR, ERK, CLDN4-S195, β-arrestin2, clathrin | internalized CLDN4 then ubiquitin-dependent loss |
| Li 2021 *FASEB J* | 10.1096/fj.202002098r | 33484473 | nutrient starvation, IPEC-J2 | CLDN1/3/4 loss is lysosome-dependent; CLDN3/4 use dynamin |
| Li 2021 *Arch Toxicol* | 10.1007/s00204-021-03044-w | 33847777 | deoxynivalenol | endocytosis + lysosomal degradation of CLDN4; dynasore / JNK |
| Ogawa 2012 *Histochem Cell Biol* | 10.1007/s00418-012-0956-x | 22544349 | EGF, ovarian serous lines | post-translational CLDN4 decrease via MEK/ERK |
| Ikari 2019 *Sci Rep* | 10.1038/s41598-019-46250-4 | 31273276 | constitutive MDCK turnover | chloroquine (not lactacystin) raises CLDN4 |

PMTPV (*Int J Mol Sci* 2020, PMID 32824620) is lysosome-dependent **CLDN1** loss after a CLDN1 peptide. Wrong claudin.

### E. CLDN4 ligands that enter endosome/lysosome — right receptor, cargo not CLDN4 protein

**Yuan / Saeki et al., *BMC Cancer* 2011** — DOI `10.1186/1471-2407-11-61`, PMID 21303546.

CPE peptide targeting CLDN3/4 delivers gelonin into vesicles. The **toxin** is trapped in the endosomal/lysosomal compartment unless an R9 escape motif is added. This is cargo trafficking, not a CLDN4-degradation-after-ADC experiment.

**Hashimoto et al., *Int J Oncol* 2013** — DOI `10.3892/ijo.2013.1881`, PMID 23563899.

CPE-ETA' binds CLDN4 and is drawn as the ETA endosome → Golgi/ER route. Immunotoxin, not ADC lysosomal payload release, and not CLDN4 protein fate.

**Patent WO2016175551A2** describes a CLDN4 antibody–dendron–doxorubicin conjugate with claimed acidic endosome/lysosome release. Not a paper; not counted as FOUND.

CLDN6-23-ADC (McDermott 2023, PMID 36884217) and CLDN18.2 ADCs internalize the **wrong claudin**.

## Files

| File | Contents |
| --- | --- |
| `VERDICT.txt` | One-screen found vs none |
| `catalog.tsv` | All inspected items with status |
| `found.tsv` | Qualifying hits (header only; 0 rows) |
| `none.tsv` | Direct analog lanes |
| `near_misses.tsv` | Related papers that fail the analog as written |
| `excluded_private.tsv` | C2 marked private |
| `search_log.md` | Query-by-query log |
