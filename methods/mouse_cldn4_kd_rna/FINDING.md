# Mouse Cldn4 KO/KD RNA-seq hunt (public only)

**Hunt date:** 2026-08-18  
**Question:** After Cldn4 loss, do IFN/MHC go up and tight-junction (TJ) genes go down?  
**Scope:** Public *Mus musculus* Cldn4 knockout or knockdown RNA-seq. Lung preferred; any epithelial/tumor if no lung. Human lung CLDN4-KD RNA-seq is already known absent. User private 8-KL data were not used.

## Verdict

One public mouse Cldn4-loss RNA-seq series exists: **GSE50927** (whole lung, germline Cldn4 KO ± ventilator-induced lung injury; Kage/Flodby/Borok 2014).

**Scored contrast (KD-direction analog):** naive Cldn4 KO vs naive WT lung, author-deposited edgeR table `GSE50927_Cldn4lungWTvsKOgenes.csv`.

| Module | After Cldn4 loss | Score |
|---|---|---|
| Cldn4 itself | Down (logFC −6.06, FDR 4.1e−26) | Confirmed KO |
| IFN / ISG | Up (19/57 panel genes FDR<0.05; all 19 up) | **Supported** |
| MHC / antigen presentation | Mixed. B2m and Psmb9 up at FDR<0.05. Several H2-Q/T/M/Bl genes are strongly up, but the line is mixed 129/B6/BALB and H2-K1 is unchanged — haplotype/mapping risk. H2-Oa/Ob down. | **Partial, haplotype-caveated** |
| TJ program (other than Cldn4) | Not down at FDR<0.05. Cldn3/5/7/18 and Tjp1/Ocln/Cdh1 unchanged, matching the paper’s AT2 qPCR/Western. | **Not supported** |

**Overall analog score:** IFN-up **yes**; MHC-up **partial**; TJ-down **no** beyond the deleted gene. This is a usable public KD-direction analog for IFN, not for a general TJ-collapse program.

No other public mouse Cldn4 KO/KD RNA-seq (intestine, organoid, tumor, or additional lung) was found. Human CLDN4-loss expression sets exist (breast RNA-seq GSE207704; ovarian microarray GSE22493) and were catalogued only as exclusions.

Companion table: `GSE50927_naiveKO_vs_WT_IFN_MHC_TJ.tsv`.

---

## 1. Hunt methods (no fabrication)

Sources queried 2026-08-18:

1. **NCBI GEO / GDS** via E-utilities (`esearch` + `esummary` on `gds`).
2. **Literature** (PubMed/PMC/Google) for Cldn4/CLDN4 knockout, knockdown, siRNA, shRNA, CRISPR + RNA-seq / transcriptome / GSE / E-MTAB.
3. **OmicsDI / ArrayExpress / BioStudies** cross-checks of GEO hits (E-GEOD-50927 = GSE50927).
4. **SRA/BioProject** for GSE50927: PRJNA219385, SRP029974.

Query families (representative):

- `(Cldn4 OR CLDN4 OR "claudin 4" OR "claudin-4") AND (knockout OR knock-out OR knockdown OR siRNA OR shRNA OR CRISPR) AND "Mus musculus"[Organism]`
- `Cldn4 AND "Mus musculus"[Organism] AND ("expression profiling by high throughput sequencing" OR RNA-seq)`
- `(CLDN4 OR Cldn4) AND (knockdown OR knockout OR siRNA OR shRNA OR CRISPR) AND ("expression profiling by high throughput sequencing")` (any organism)
- `Cldn4-/- OR "Cldn4 deficient" OR "claudin-4-deficient"`
- Literature: intestine / organoid / urothelium / tumor / lung injury Cldn4 KO RNA-seq

Inclusion for scoring: *Mus musculus*, Cldn4 is the perturbed gene (KO or KD), transcriptome is RNA-seq, and a KO/KD vs control contrast exists.

Exclusion: Cldn4 used only as a marker (DATP/Krt8+ AT2, CLDN4+ sorting); other claudin KOs (Cldn7, Cldn18); human CLDN4 KD/KO; microarray-only unless no RNA-seq exists (none did for mouse Cldn4-loss); private 8-KL.

Scoring used the **author-deposited edgeR tables** (not a re-alignment). Sign convention: Cldn4 logFC is negative, so logFC = KO − WT.

---

## 2. Catalog

### 2.1 Scored (include)

| Accession | Species / tissue | Perturbation | Assay | n (GEO) | Why scored |
|---|---|---|---|---|---|
| **GSE50927** (E-GEOD-50927; PRJNA219385; SRP029974) | Mouse whole lung | Germline Cldn4 KO vs WT; naive and VILI | RNA-seq, HiSeq 2000, SR 100 nt, mm9/edgeR | 5 GEO samples (paper: 2 mice/condition, “all in duplicates”) | Only public mouse Cldn4-loss RNA-seq. Lung matches the preferred tissue. PMID 25106430 |

Samples (GSM):

| GSM | Title | Role in analog |
|---|---|---|
| GSM1232580 | WT no VILI | Naive WT |
| GSM1232581 | Cldn4 KO no VILI | Naive KO (**primary**) |
| GSM1232582 | WT VILI | Injury confounder |
| GSM1232583 | Cldn4 KO VILIlow | Injury + genotype |
| GSM1232584 | Cldn4 KO VILIhigh | Injury + genotype |

Deposited DE files:

- `GSE50927_Cldn4lungWTvsKOgenes.csv` — naive KO vs WT (**scored**)
- `GSE50927_VILIwtkohiGenes.csv` — KOhigh VILI vs WT VILI (injury-confounded)
- `GSE50927_VILIwtkoloGenes.csv` — KOlow VILI vs WT VILI (injury-confounded)
- `GSE50927_VILIwtGenes.csv` — WT VILI vs naive WT (not a Cldn4-loss contrast; Cldn4 **up** 16-fold, as reported)

### 2.2 Related mouse, not Cldn4-loss RNA-seq (exclude from scoring)

| Accession | What it is | Why excluded |
|---|---|---|
| GSE48443 / GDS4961 | Cldn**18** KO lung microarray | Wrong gene |
| GSE106233 | Cldn18 KO AT2 microarray | Wrong gene |
| GSE84742 | WT colon crypt LCM microarray; Hopx/Klf4 → Cldn4 | Cldn4 is a readout, not KO/KD |
| GSE11759 / GDS5284 | Hnf4a colon KO microarray | Not Cldn4 perturbation |
| GSE310539 / GSE309751 | AP-1 multiome; CLDN4+ transitional AT2 | Marker, not Cldn4 KO/KD |
| GSE247130 / GSE247271 / GSE264098 | Cebpa AT2 RNA/ChIP/ATAC; KRT8/CLDN4+ states | Marker, not Cldn4 KO/KD |
| GSE244820 | Coronavirus AT2 ER-stress RNA-seq | Mentions Cldn4 biology; not Cldn4 KO/KD |
| GSE108417 | Mouse asthma airway array | Cldn4 as barrier gene, not KO/KD |

### 2.3 Human CLDN4-loss expression (exclude; lung RNA-seq still absent)

| Accession | What it is | Why excluded |
|---|---|---|
| GSE207704 | Human T47D/MCF-7 CLDN4−/− breast RNA-seq (PMID 37069688) | Human tumor, not mouse; not lung |
| GSE22493 | SKOV-3 CLDN4 siRNA lentivirus microarray (n=3) | Human, microarray |
| GSE22421 | SKOV-3 C-CPE (CLDN4-directed) microarray | Pharmacologic, human, microarray |
| GSE99417 / GSE99415 / GSE99416 | CLDN4 ceRNA network arrays (gastric) | Not a Cldn4 KO/KD RNA-seq of mouse epithelium |

No public **human lung** CLDN4-KD/KO RNA-seq was found (NHBE/AEC CLDN4 siRNA papers are qPCR/TEER only: e.g. Wray 2009 AJP Lung; Lee 2018 AAIR).

### 2.4 Mouse Cldn4 KO papers without deposited RNA-seq

| Paper | Model | Transcriptome? |
|---|---|---|
| Fujita et al. 2012 PLoS One 7:e52272 | Germline Cldn4−/−; urothelial hyperplasia / hydronephrosis | qPCR of other Cldns in kidney/bladder; **no RNA-seq accession** |
| Kage / Flodby / Borok 2014 AJP Lung (this GEO) | Germline Cldn4 KO lung ± VILI/hyperoxia | **GSE50927** |
| Horng et al. 2017 JCI (astrocyte mGfap-Cre Cldn4fl/fl) | CNS conditional Cldn4 CKO | Lesion histology/IHC; **no RNA-seq accession found** |
| Deane / Lyons et al. 2025 Shock (PMID-linked “Claudin 4 Deletion Improves Gut Permeability…”) | Germline Cldn4−/− abdominal sepsis | Permeability, cytokines, flow cytometry; **no GEO/RNA-seq found** |

Intestine/organoid Cldn4 KO RNA-seq was specifically searched and **not found**. Colon organoid RNA-seq hits were Smad4 / Cldn7 / other genes, with Cldn4 as a downstream TJ readout.

---

## 3. Scored accession: GSE50927

**Citation:** Kage H, Flodby P, Gao D, Kim YH, Marconett CN, DeMaio L, Kim KJ, Crandall ED, Borok Z. *Claudin 4 knockout mice: normal physiological phenotype with increased susceptibility to lung injury.* Am J Physiol Lung Cell Mol Physiol 307:L524–L536 (2014). PMID 25106430. GEO GSE50927.

**Model:** Floxed Cldn4 from 129S6 BAC in W4 ES cells; CMV-Cre deletion; cre selected away on 129S6/SvEvTac. Paper: “parallel WT line … same mixed genetic background” (129S6/SvEvTac, C57BL/6, and BALB/c). Later MMRRC stock B6.Cg-Cldn4tm1.2Zboro was backcrossed to B6; **the 2013/2014 RNA-seq is the mixed-background line**.

**Why naive KO vs WT is the analog:** VILI induces Cldn4 16-fold in WT and a large injury transcriptome (Egr1/Tnf/Il1b/Il6). The paper used VILI KOhigh vs KOlow subtraction to find injury mediators, not the Cldn4-loss program. Naive KO vs WT is the only contrast that isolates genotype without mechanical injury.

**Limitations (do not over-read):**

- n = 2 mice per condition; GEO lists 5 samples (duplicates likely pooled or deposited as one sample/condition).
- Whole-lung bulk, not sorted epithelium.
- Mixed 129 / B6 / BALB background: BALB/c is H2^d vs H2^b on B6/129. Extreme H2-K2 / H2-M2 / H2-Bl fold-changes are **not** treated as clean MHC induction.
- Author alignment: Bowtie to UCSC mm9 RefSeq; edgeR after EDAseq. MHC gene annotation on mm9 is noisy.
- Paper itself: baseline physiology “normal”; AT2 Cldn3/5/7/18 mRNA and Cldn3/18 protein unchanged.

---

## 4. KD IFN table (naive KO vs WT)

Source: `GSE50927_Cldn4lungWTvsKOgenes.csv` (downloaded from GEO 2026-08-18). logFC is KO − WT.

Pre-specified IFN/ISG panel (mouse symbols). Rows with FDR < 0.05 **or** (|logFC| ≥ 0.5 and logCPM > −1) are shown. Full panel is in the TSV.

| Gene | Module | logFC | logCPM | P | FDR | Direction |
|---|---|---:|---:|---:|---:|---|
| Isg15 | IFN | +1.075 | 2.79 | 3.9e-05 | 0.0031 | up |
| Oas1g | IFN | +1.192 | 1.84 | 7.2e-04 | 0.028 | up |
| Oas2 | IFN | +1.489 | 2.67 | 1.7e-05 | 0.0015 | up |
| Oas3 | IFN | +1.816 | 1.42 | 1.9e-05 | 0.0016 | up |
| Oasl2 | IFN | +0.542 | 5.24 | 3.7e-04 | 0.018 | up |
| Ifit1 | IFN | +0.576 | 4.62 | 7.2e-04 | 0.028 | up |
| Ifitm1 | IFN | +0.528 | 5.62 | 2.8e-04 | 0.014 | up |
| Ifitm3 | IFN | +0.471 | 8.09 | 8.6e-04 | 0.032 | up |
| Gbp3 | IFN | +0.718 | 5.27 | 5.1e-04 | 0.022 | up |
| Gbp4 | IFN | +2.385 | 6.97 | 1.6e-09 | 4.8e-07 | up |
| Gbp5 | IFN | +1.067 | 4.54 | 4.0e-04 | 0.019 | up |
| Gbp6 | IFN | +1.371 | 3.27 | 4.3e-06 | 4.9e-04 | up |
| Ifi44 | IFN | +0.749 | 3.36 | 8.4e-04 | 0.031 | up |
| Ifi47 | IFN | +1.029 | 5.46 | 8.3e-07 | 1.3e-04 | up |
| Irgm2 | IFN | +1.221 | 6.27 | 1.9e-05 | 0.0017 | up |
| Igtp | IFN | +1.384 | 5.46 | 1.3e-05 | 0.0012 | up |
| Tgtp2 | IFN | +1.153 | 6.39 | 6.4e-04 | 0.026 | up |
| Ubd | IFN | +5.615 | −0.75 | 8.3e-04 | 0.031 | up (low CPM) |
| Irg1 (Acod1) | IFN/inflam | +5.492 | −0.82 | 0.0014 | 0.044 | up (low CPM) |
| Stat1 | IFN | +0.391 | 6.95 | 0.026 | 0.29 | up, ns FDR |
| Stat2 | IFN | +0.446 | 5.91 | 0.0028 | 0.071 | up, ns FDR |
| Irf1 | IFN | +0.382 | 7.04 | 0.015 | 0.21 | up, ns FDR |
| Irf7 | IFN | +0.227 | 3.98 | 0.33 | 0.98 | ns |
| Cxcl9 | IFN | +2.883 | 1.34 | 0.0023 | 0.062 | up, ns FDR |
| Cxcl10 | IFN | +2.543 | 0.76 | 0.0043 | 0.094 | up, ns FDR |
| Nlrc5 | IFN/MHC | +0.748 | 4.59 | 0.0042 | 0.093 | up, ns FDR |
| Ifit2 | IFN | −0.369 | 5.38 | 0.013 | 0.19 | down, ns FDR |

**IFN call:** After Cldn4 loss, a type I/II-like ISG cassette is up (OAS, GBP, IFITM, Isg15, Irgm2/Igtp/Tgtp2). Core TFs (Stat1/2, Irf7) move in the same direction but are not FDR-significant at n=2.

---

## 5. MHC and TJ (same contrast)

### MHC / antigen presentation

| Gene | logFC | logCPM | P | FDR | Note |
|---|---:|---:|---:|---:|---|
| B2m | +0.655 | 9.73 | 9.3e-04 | 0.034 | MHC-I light chain; less haplotype-sensitive |
| Psmb9 | +0.902 | 5.01 | 0.0011 | 0.037 | immunoproteasome |
| Psmb8 | +0.432 | 6.02 | 0.026 | 0.29 | ns FDR |
| H2-K1 | −0.071 | 9.66 | 0.72 | 1 | classical MHC-I **unchanged** |
| H2-D1 | +0.355 | 9.94 | 0.067 | 0.48 | ns |
| H2-K2 | +4.262 | 6.61 | 1.9e-69 | 1.4e-65 | **haplotype/mapping suspect** |
| H2-M2 | +5.182 | 4.00 | 6.3e-58 | 3.4e-54 | **haplotype/mapping suspect** |
| H2-Bl | +3.901 | 3.55 | 4.9e-28 | 9.6e-25 | **haplotype/mapping suspect** |
| H2-Q1 | +1.071 | 6.46 | 1.2e-07 | 2.4e-05 | Q-family; interpret with strain caution |
| H2-Q2 | +1.329 | 4.85 | 3.1e-09 | 8.6e-07 | same |
| H2-Q10 | +2.173 | 6.05 | 4.4e-26 | 7.9e-23 | same |
| H2-T23 | +1.045 | 7.36 | 2.9e-06 | 3.7e-04 | Qa-1 |
| H2-M3 | +1.191 | 4.23 | 3.7e-08 | 8.0e-06 | |
| H2-Oa | −1.400 | 2.20 | 1.0e-05 | 0.0010 | MHC-II DO down |
| H2-Ob | −2.288 | 4.50 | 6.8e-26 | 1.1e-22 | MHC-II DO down |
| Tap1 | +0.184 | 6.19 | 0.44 | 1 | ns |
| Cd74 | +0.486 | 10.10 | 0.0037 | 0.085 | ns FDR |

**MHC call:** Antigen-presentation machinery that is not a polymorphic class I gene (B2m, Psmb9) is up. Classical H2-K1 is not. Extreme H2-K2/M2/Bl changes are not scored as MHC induction.

### Tight junction (predicted: down)

| Gene | logFC | logCPM | P | FDR |
|---|---:|---:|---:|---:|
| **Cldn4** | **−6.061** | 2.25 | 1.8e-29 | **4.1e-26** |
| Cldn1 | −0.080 | 3.83 | 0.72 | 1 |
| Cldn3 | +0.091 | 7.13 | 0.52 | 1 |
| Cldn5 | +0.208 | 8.81 | 0.14 | 1 |
| Cldn7 | +0.113 | 5.42 | 0.47 | 1 |
| Cldn8 | −0.549 | 3.52 | 0.012 | 0.19 |
| Cldn15 | −0.608 | 3.19 | 0.059 | 0.45 |
| Cldn18 | −0.452 | 10.00 | 0.041 | 0.38 |
| Tjp1 | −0.104 | 9.29 | 0.46 | 1 |
| Tjp2 | −0.013 | 7.44 | 0.94 | 1 |
| Ocln | −0.046 | 6.31 | 0.76 | 1 |
| F11r | −0.117 | 7.79 | 0.46 | 1 |
| Cdh1 | −0.140 | 8.92 | 0.37 | 1 |
| Epcam | −0.088 | 5.99 | 0.53 | 1 |

**TJ call:** Only Cldn4 is FDR-significant. This matches the paper: “expression of all other claudin genes … unchanged.” A general TJ-down program is **not** observed in naive Cldn4 KO lung.

---

## 6. Secondary contrasts (not scored as the analog)

VILI WT vs naive WT: Cldn4 logFC **+3.96** (FDR 7e−91); IFN/MHC mostly **down**; injury cytokines (Il6, Cxcl1/2, Tnf, Il1b) up. Injury is the opposite of the Cldn4-loss IFN direction.

VILI KOhigh vs WT VILI: Cldn4 logFC −10.2; residual H2-K2/M2/Bl up (same haplotype pattern); Irf1, Gbp4, Rsad2 up at FDR<0.05; Egr1/Tnf/Il1b up as the paper reported. Confounded by greater leak/injury.

---

## 7. What this does and does not claim

- Public mouse lung Cldn4 KO RNA-seq **exists** (GSE50927). It is the only scored accession.
- After Cldn4 loss at baseline, **IFN/ISG transcripts are up**.
- **MHC is not a clean yes:** B2m/Psmb9 up; classical H2-K1 flat; several H2 genes look like strain/mapping.
- **TJ-down besides Cldn4 is not supported** in this lung KO.
- No public mouse intestine/organoid/tumor Cldn4 KO RNA-seq was found.
- Human lung CLDN4-KD RNA-seq remains absent.
- Private 8-KL samples were not used.

Done criterion: catalog exists; KD IFN table exists for the one scored accession.
