# FINDING — PXD031094 CLDN4 bait CoIP partners (MaxQuant SEARCH)

**Verdict: ranked CLDN4 CoIP partners exist in the public MaxQuant table; IFN and TJ genes are in the primary list; classical MHC is detected but not primary.** Human CLDN4 bait (O14493) is recovered in 4/4 Cldn4 IPs and 0/4 GFP controls (log2FC 10.38). After removing reverse/contaminant/site-only rows and requiring ≥2 unique peptides and LFQ in ≥2 Cldn4 replicates, **81 non-bait primary partners** remain (stringent FDR_mod < 0.05 and log2FC > 1, or exclusive: ≥2 Cldn4 / 0 GFP).

This is **not** abundance vs ICI. No proteins were added that are absent from `proteinGroups_Cldn4.txt`.

---

## English

### Data

| Item | Value |
| --- | --- |
| Accession | [PXD031094](https://www.ebi.ac.uk/pride/archive/projects/PXD031094) |
| Paper | Suarez-Artiles et al. *Cell Rep.* 2022; 41:111588 (PMID 36351382) |
| File used | PRIDE SEARCH `proteinGroups_Cldn4.txt` (3.8 MB; SHA-1 `4dbde6c0d827105412ef1f34cbb969b8dd8acb36`) |
| Skipped | all 610 RAW files; PRISMA `proteinGroups_PRISMA.txt` (C-terminal peptide screen, not full-length bait CoIP) |
| Search | MaxQuant 1.5.2.8 vs *Canis lupus familiaris* UniProt 2018; Fast LFQ; MBR on |
| Design | YFP/CFP-CLDN4 GFP-Trap CoIP in MDCK-C7 vs cytosolic eGFP; **n=4 vs 4** |
| Protein groups in file | 2788 |
| After reverse / contaminant / site-only / unique peptides ≥2 | 2188 |
| LFQ in ≥2 Cldn4 replicates (ranked universe) | **1778** |
| Primary (stringent or exclusive) | **82** (1 bait + 81 partners) |
| Author text | CLDN4 had **24** significant interactions (their moderated t + GFP-derived second cutoff). We cannot replay that exact call without Table S1. |

Gene symbols are MaxQuant `Gene names` or UniProt `GN=` from FASTA headers. One primary row (E2RRA6) has **no** `GN=` and is left blank — not named.

### Ranking

1. Drop Reverse, Potential contaminant, Only identified by site.  
2. Keep unique peptides ≥2 and LFQ > 0 in ≥2 of 4 Cldn4 samples.  
3. log2 LFQ; zeros imputed at global min-positive / 2 (`1.8321e6`).  
4. Two-sample comparison Cldn4 vs GFP: Welch t and a simple variance-moderated t (df_prior=10; **not** limma). BH-FDR over the 2188 filtered groups.  
5. **Primary partner** = (FDR_mod < 0.05 and log2FC > 1) **or** exclusive (≥2 Cldn4, 0 GFP).  
6. Rank primary by bait first, then stringent, then exclusive, then log2FC.

Category labels are annotation only (MSigDB Hallmark/Reactome/KEGG/GO + a compact MHC/TJ core). MHC uses the 29-gene peptide-loading set plus DLA88, not the 381-gene “class I processing” dump.

### Primary list — IFN / MHC / TJ / trafficking / immune

Non-bait primary proteins that map to those sets (from `results/cldn4_coip_category_hits_primary.tsv`):

| Rank in category table | Gene | UniProt | log2FC | FDR_mod | Cldn4 / GFP LFQ n | Call | Categories |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 2 | **OAS1** | F1PLW6 | 5.76 | 8.9e-13 | 4 / 0 | exclusive + stringent | IFN; immune |
| 3 | **CLDN3** | Q95KM5 | 2.49 | 0.014 | 4 / 2 | stringent | TJ |
| 4 | **MPP7** | F1P8X0 | 1.50 | 2.6e-4 | 4 / 4 | stringent | TJ |
| 5 | **ISG15** | E2R7R1 | 1.27 | 0.0053 | 4 / 4 | stringent | IFN; immune |
| 6 | **OAS2** | F1PTR7 | 1.26 | 0.0021 | 4 / 4 | stringent | IFN; immune |
| 7 | KIF22 | E2RRZ0 | 1.20 | 0.0036 | 4 / 4 | stringent | trafficking |
| 8 | CDC42 | P60952 | 1.00 | 0.0032 | 4 / 4 | stringent | TJ |
| 9 | SYNPO | J9P7S2 | 1.00 | 0.0060 | 4 / 4 | stringent | TJ |
| 10 | **PARD3B** | F1PUY2 | 2.66 | 0.056 | 2 / 0 | exclusive | TJ |
| 11 | **EPN2** | F1P9H0 | 2.47 | 0.056 | 2 / 0 | exclusive | trafficking |
| 12 | **ICAM1** | P33729 | 2.21 | 0.056 | 2 / 0 | exclusive | IFN; TJ; immune |
| 13 | **RABL3** | E2R3T1 | 2.16 | 0.056 | 2 / 0 | exclusive | trafficking |

Primary non-bait category counts: IFN 4, TJ 6, trafficking 3, immune 4, **MHC 0**.

### MHC (detected, not primary)

No peptide-loading / DLA class I gene passes the primary cutoff. These MHC-set proteins **are in the MaxQuant table** (not invented) but fail FDR or log2FC:

| Gene | UniProt | log2FC | FDR_mod | Cldn4 / GFP n | Note |
| --- | --- | --- | --- | --- | --- |
| TAPBP | Q5TJE4 | 1.33 | 0.32 | 3 / 2 | tapasin; peptide loading |
| PDIA3 | E2RD86 | 0.44 | 0.13 | 4 / 4 | MHC folding set |
| DLA88 | P18466 | 0.27 | 0.35 | 4 / 4 | canine MHC class I (HA19) |
| CANX | P24643 | −3.81 | 0.013 | 2 / 4 | higher in GFP |

### Exclusive non-bait partners (n=27)

GFP-absent, ≥2 Cldn4 LFQ, unique peptides ≥2. Includes the strongest IFN hit (**OAS1**) and TJ/trafficking exclusives (PARD3B, EPN2, ICAM1, RABL3, SDCBP). Full list: `results/cldn4_coip_partners_primary.tsv` rows with `exclusive=True` and `is_bait=False`.

Top exclusive by log2FC (bait omitted): OAS1, INCENP, EHMT1, SDCBP, PTPRA, TIMP3, DDX49, LOC488767, MGAT4B, then n=2 exclusives (LRRC8A, CDCA8, CERS3, PARD3B, …). **LOC488767** (J9P1F7) stays LOC488767 — UniProt has no approved symbol.

### Other recovered junction / trafficking proteins (not primary)

Present in ≥2 Cldn4 IPs but not primary: OCLN (log2FC 2.08, FDR 0.15), CGN (1.53, 0.27), PATJ (0.41, 0.12), PARD3 (0.27, 0.29), TJP2 (0.17, 0.52), SLC2A1 / GLUT-1 (1.11, 0.0018 — **is** stringent, not in the five sets). **TJP1** has 0/4 Cldn4 LFQ and 1/4 GFP — not a CLDN4 partner in this table.

### What this is not

- Not ICI response, not tumor vs NAT abundance, not CPTAC.  
- Not a re-search of RAW files.  
- Not the authors’ Table S1 call of exactly 24.  
- Not PRISMA C-terminal SLiM partners.  
- Not human HLA-A/B/C (MDCK-C7 is dog; class I gene in the table is **DLA88**).

### Files

`results/cldn4_coip_partners_primary.tsv` — ranked primary partner list (the deliverable).  
`results/cldn4_coip_partners_ranked.tsv` — all 1778 bait-detected groups.  
`results/cldn4_coip_category_hits_primary.tsv` — primary rows in IFN/MHC/TJ/trafficking/immune.  
`results/cldn4_coip_category_mapped.tsv` — set membership in the 1778.  
`results/summary.json` — counts.  
`data/proteinGroups_Cldn4.txt` — public SEARCH table.

---

## 中文

**结论：公开 MaxQuant 表里有可排名的 CLDN4 CoIP 伴侣；IFN 和紧密连接基因在主名单里；经典 MHC 被鉴定到但未过主阈值。** 人源 bait CLDN4（O14493）4/4 出现、GFP 对照 0/4。过滤后主名单 **81 个非 bait 伴侣**（FDR_mod<0.05 且 log2FC>1，或 Cldn4≥2 / GFP=0）。

IFN：OAS1（仅 bait）、ISG15、OAS2、ICAM1。TJ：CLDN3、MPP7、CDC42、SYNPO、PARD3B、ICAM1。转运：KIF22、EPN2、RABL3。MHC 肽装载/犬类 I 类（TAPBP、DLA88）在表中，但不是主伴侣。没有编造蛋白。这不是 ICI 丰度比较。未用 RAW。
