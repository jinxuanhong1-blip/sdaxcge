# Claim A10 — honest verdict

**Claim (operationalised in `claims/claim_A10.md` before these numbers were computed):**
ELF3, GRHL1, KLF4 and TFAP2A form a coherent transcription-factor module that is
positively coupled to TACSTD2 and CLDN4, in a lung-epithelial state where NKX2-1
is down. Expected in human and mouse public data.

**Overall: NOT SUPPORTED.**  
The four TFs are not a module. ELF3 does travel with TACSTD2/CLDN4. NKX2-1 is
**not** systematically down in that program; in normal human lung it is
*positively* correlated with ELF3, TACSTD2 and CLDN4. Experimental NKX2-1 loss
does not raise the four-TF set in human LUAD lines. Mouse *Nkx2-1*-null LUAD
and AT1 deletion raise *Elf3* and *Cldn4* — that is a real, narrower finding,
not the claim as stated.

Decision rules were pre-registered: SUPPORT requires |ρ| ≥ 0.3, FDR < 0.05,
extreme background percentile, and reproduction in ≥ 2 datasets per species,
plus P6 (not just bulk composition).

Machine-readable labels: `tables/verdicts.json`. Number dump: `TABLES_DUMP.md`.

---

## Per-prediction verdicts

| id | prediction | verdict |
|---|---|---|
| **P1** | ELF3/GRHL1/KLF4/TFAP2A are a coherent module | **NOT SUPPORTED** |
| **P2** | Module / each TF correlates + with TACSTD2 and CLDN4 | **PARTIALLY SUPPORTED** (ELF3 and TACSTD2–CLDN4; not the four-TF set) |
| **P3** | NKX2-1 anti-correlates with the module and the targets | **NOT SUPPORTED** (opposite in GTEx; null / mixed in TCGA) |
| **P4** | Experimental NKX2-1 loss raises MODULE and TARGET | **NOT SUPPORTED** in human lines; **partial** in mouse LUAD/AT1 for Elf3/Cldn4 only |
| **P5** | Same directions in mouse | **PARTIALLY SUPPORTED** for Elf3/Cldn4 after Nkx2-1 loss; observational mouse lung is underpowered for airway types |
| **P6** | Associations survive within an epithelial type | **PARTIALLY SUPPORTED** for ELF3–TACSTD2/CLDN4 (positive, modest); **fails** for NKX2-1-down |

---

## Data actually used (all public)

| dataset | n | role |
|---|---|---|
| GTEx v8 lung, recount3 G026 | 655 | human bulk, normal |
| TCGA-LUAD, recount3 G026 | 601 (542 tumour / 59 adjacent) | human bulk, tumour |
| CELLxGENE Census 2025-11-08, human lung epithelium, primary, normal + LUAD | 895,386 cells | P6 / cell-type means |
| CELLxGENE Census, mouse lung primary | 220,706 cells | P5/P6 |
| GSE129340 H441 / H209 TTF-1 siRNA | 2 KD + 1 control each | P4 human (underpowered) |
| GSE229541 NKX2-1 sg/sh/CRISPRi/OE in LUAD lines | 11 named contrasts | P4 human |
| GSE129583 *Nkx2-1* mutant AT1 (P5) and AT2 (P8/P9) | 3 vs 3 each | P4 mouse alveolus |
| GSE115899 Kras LUAD *Nkx2-1*− vs + | 3 vs 3 | P4/P5 mouse tumour |
| GSE145152 Braf/p53 LUAD *Nkx2-1* f/f vs f/+, control chow | 4 vs 5 | P4/P5 mouse tumour |

GSE188435 was downloaded and **excluded**: it deletes FoxA1/2 in an NKX2-1-positive
background and does not manipulate NKX2-1.

AGER and PTPRC resolved to all-zero rows in the recount3 G026 symbol collapse
(annotation mismatch). Composition adjustment therefore uses EPCAM, SFTPC,
SCGB1A1, KRT5, TP63, FOXJ1, MUC5B, COL1A1, PECAM1 — not AGER/PTPRC.

---

## P1 — not a four-TF module

Mean pairwise Spearman among ELF3/GRHL1/KLF4/TFAP2A:

| dataset | mean pairwise ρ | matched-null mean | matched-null p95 | empirical p | random-pair percentile |
|---|---:|---:|---:|---:|---:|
| GTEx-lung | 0.157 | 0.014 | 0.200 | 0.084 | 68 |
| TCGA-LUAD | 0.099 | 0.031 | 0.172 | 0.171 | 64 |

Neither dataset clears |ρ| ≥ 0.3 or the matched-null 95th percentile. GRHL1 and
TFAP2A are weakly expressed in GTEx lung (TFAP2A mean log2CPM 0.86; detected
>1 in only 28% of samples) and do not co-vary with KLF4 (GTEx KLF4–TFAP2A
ρ = −0.01; TCGA KLF4–GRHL1 ρ = −0.15). Calling these four genes a “module”
is not justified by the data.

---

## P2 — ELF3 (not the four-TF set) couples to TACSTD2/CLDN4

| dataset | pair | ρ | 95% CI | FDR | bg percentile |
|---|---|---:|---|---:|---:|
| GTEx-lung | ELF3–TACSTD2 | 0.409 | [0.34, 0.47] | 3.3e-27 | 93 |
| GTEx-lung | ELF3–CLDN4 | 0.651 | [0.60, 0.70] | 2.7e-79 | 99 |
| GTEx-lung | TACSTD2–CLDN4 | 0.756 | [0.72, 0.79] | 6.9e-121 | 100 |
| GTEx-lung | GRHL1–TACSTD2 | −0.111 | | 6.4e-3 | 26 |
| GTEx-lung | TFAP2A–TACSTD2 | −0.017 | | 0.68 | 40 |
| TCGA-LUAD | ELF3–TACSTD2 | 0.376 | [0.30, 0.45] | 6.6e-21 | 96 |
| TCGA-LUAD | ELF3–CLDN4 | 0.567 | [0.51, 0.62] | 5.6e-51 | 100 |
| TCGA-LUAD | TACSTD2–CLDN4 | 0.520 | [0.46, 0.58] | 9.3e-42 | 99 |
| TCGA-LUAD | GRHL1–TACSTD2 | 0.409 | | 1.1e-24 | 97 |
| TCGA-LUAD | KLF4–CLDN4 | −0.219 | | 1.4e-7 | 8 |

Module-score vs TACSTD2/CLDN4 is positive (GTEx 0.28 / 0.44; TCGA 0.51 / 0.48)
but that score is dominated by ELF3. GRHL1 and TFAP2A do not reproduce the
target coupling in normal lung.

Partial correlation: ELF3–CLDN4 stays large after residualising on EPCAM +
SFTPC + SCGB1A1 + all lineage markers (GTEx 0.65 → 0.60; TCGA 0.57 → 0.45).
ELF3–TACSTD2 shrinks in GTEx after EPCAM (0.41 → 0.19) but remains positive
and FDR-significant. So P2 for **ELF3–CLDN4** is not only “how much epithelium
is in the biopsy.”

---

## P3 — NKX2-1 is not down

This is the part of the claim that fails most clearly.

| dataset | pair | ρ | FDR | bg percentile | claim direction? |
|---|---|---:|---:|---:|---|
| GTEx-lung | NKX2-1–ELF3 | **+0.351** | 6.7e-20 | 89 | no (opposite) |
| GTEx-lung | NKX2-1–TACSTD2 | **+0.558** | 3.5e-54 | 98 | no (opposite) |
| GTEx-lung | NKX2-1–CLDN4 | **+0.715** | 3.1e-102 | 100 | no (opposite) |
| TCGA-LUAD | NKX2-1–ELF3 | **+0.201** | 1.5e-6 | 82 | no |
| TCGA-LUAD | NKX2-1–TACSTD2 | −0.034 | 0.45 | 35 | null |
| TCGA-LUAD | NKX2-1–CLDN4 | **+0.242** | 4.7e-9 | 87 | no (opposite) |
| TCGA-LUAD | NKX2-1–KLF4 | −0.383 | 1.4e-21 | 0.9 | yes, KLF4 only |
| TCGA-LUAD | NKX2-1–TFAP2A | −0.218 | 1.6e-7 | 7.8 | yes, TFAP2A only |

Module score vs NKX2-1: GTEx **+0.16**, TCGA **−0.13**. Neither is a
reproducible anti-correlation of a module.

TCGA NKX2-1-low vs NKX2-1-high tumours (n=136 vs 136, quartile split):

| gene | log2FC (low − high) | Cliff's δ | FDR | claim? |
|---|---:|---:|---:|---|
| ELF3 | **−0.61** | −0.28 | 1.6e-4 | opposite |
| GRHL1 | **−0.35** | −0.20 | 7.3e-3 | opposite |
| KLF4 | +1.47 | +0.68 | 1.8e-21 | yes |
| TFAP2A | +0.96 | +0.38 | 3.5e-7 | yes |
| TACSTD2 | −0.08 | +0.10 | 0.19 | null |
| CLDN4 | **−0.73** | −0.36 | 8.0e-7 | opposite |

NKX2-1-low LUAD is a KRT5-high / SFTPB-low state (KRT5 +1.33; SFTPB −4.03),
not an ELF3/TACSTD2/CLDN4-high state. CLDN4 and ELF3 are *lower* where NKX2-1
is lower.

This agrees with the independent TCGA-LUAD analysis in PR #103 (NKX2-1 vs
TACSTD2 null; vs CLDN4 positive). That PR did not test the ELF3 set; this one
does, and the inverse still does not appear.

After composition adjustment, GTEx NKX2-1–TACSTD2 collapses from 0.56 to ~0
(EPCAM residual), i.e. it was epithelial content. NKX2-1–CLDN4 stays positive
(0.71 → 0.16–0.38). Nothing becomes a significant inverse.

---

## P4 — experimental NKX2-1 loss

Positive log2FC = higher in the NKX-low / KD arm (claim direction).
NKX2-1 itself is the positive control and does go down in every loss-of-function
contrast that has power.

**Human LUAD lines (GSE229541).** ELF3 falls or is unchanged after NKX2-1
sg/sh/CRISPRi (H358 sg −0.61; H441 sg −0.29; H441 sh −0.41; H2087 sh −0.94;
HCC78 CRISPRi −0.25). TACSTD2 and CLDN4 move <0.5 and flip sign across lines.
GSE129340 H441 TTF-1 siRNA (n=2 vs 1): TACSTD2/CLDN4 ≈ 0; ELF3 +0.31. H209 is
SCLC and is not a LUAD test. Overexpression in PC9/H1975 (NKX2-1-low parents)
lowers TACSTD2 (PC9 +0.90 and H1975 +1.46 on the GFP-minus-OE scale) but
CLDN4 goes the other way. No human contrast meets |log2FC| ≥ 0.5 with
reproducible direction for the four-TF set plus both targets.

**Mouse.** This is the only place a piece of P4 is real.

| dataset | Elf3 | Cldn4 | Tacstd2 | Klf4 | Grhl1 | Tfap2a |
|---|---:|---:|---:|---:|---:|---:|
| GSE129583 AT1 P5 mutant vs control | **+2.32** | **+2.76** | +0.52 | −0.82 | +0.84 | +0.07 |
| GSE129583 AT2 P8/P9 | −0.04 | −0.58 | +0.23 | +0.17 | 0 | −1.15 |
| GSE115899 Kras Nkx2-1− vs + | **+3.02** | **+2.62** | +0.79 | +0.95 | −0.52 | +0.49 |
| GSE145152 Braf/p53 f/f vs f/+ | **+2.09** | +0.84 | −0.41 | **+1.76** | +0.08 | +0.19 |

AT2 in GSE129583 is near floor for these genes (including *Nkx2-1* itself in
the processed table) and is not interpretable. AT1 deletion and two independent
mouse LUAD *Nkx2-1*-null models raise *Elf3* and *Cldn4*. *Tacstd2* is
inconsistent. *Grhl1* and *Tfap2a* are low and uninformative. That is **not**
“the ELF3/GRHL1/KLF4/TFAP2A module vs TACSTD2/CLDN4.”

---

## P5 / P6 — mouse observational and within-type human scRNA

**Cell-type means (Census, normal human lung, full-transcriptome `raw_sum` CPM).**
Airway types (club, basal, secretory, ciliated, goblet) have higher ELF3,
TACSTD2 and CLDN4 than AT2. NKX2-1 is present in AT2 and also detectable in
several airway labels; it is not a clean “off in the TACSTD2-high type”
switch. GRHL1/KLF4/TFAP2A means are low in almost every type.

**Within-type Spearman (8,000-cell cap per type).** ELF3–CLDN4 and
TACSTD2–CLDN4 stay positive inside basal, club, ciliated, AT1, AT2 and LUAD
malignant cells (typical ρ 0.25–0.50). ELF3–TACSTD2 is positive but smaller
(0.09–0.34). NKX2-1 vs TACSTD2/CLDN4/ELF3 is **weakly positive** in almost
every type (e.g. AT2 NKX2-1–CLDN4 0.16; basal 0.16; malignant 0.06–0.12).
The one inverse (epithelial cell of lung, NKX2-1–TACSTD2 −0.24) is not
reproduced for CLDN4 in the same label (+0.35) and is not treated as evidence.

So P6 says: the ELF3–CLDN4/TACSTD2 coupling is not *only* airway-vs-AT2
mixing. The “NKX2-1 down” half of the claim is not a within-cell program.

**Mouse Census lung** is endothelium / fibroblast / AT2-heavy. Labelled basal
cells n=34; club n=970. *Nkx2-1* is high in AT1/AT2/progenitor; *Tacstd2* is
higher in AT1 than AT2. There is no powered mouse airway atlas here to test
the four-TF module. That is a gap, not a negative result.

---

## What would have counted as support, and what we actually have

The claim needs all of: a four-TF module, coupling to both targets, NKX2-1
down, human + mouse, not just composition. We have:

- a **two-gene / three-gene epithelial program** (ELF3 with TACSTD2 and
  especially CLDN4) that is real in bulk and inside cell types;
- **no** four-TF module;
- **no** NKX2-1 inverse in human bulk, human scRNA, or human LUAD lines;
- a **mouse-specific** *Nkx2-1* loss → *Elf3*/*Cldn4* up in AT1 and in
  *Nkx2-1*-null LUAD, which is interesting and should not be papered over,
  but it does not rescue the human claim or the four-TF framing.

---

## What this is not

- Not ChIP, motif, or reporter evidence that ELF3 binds TACSTD2 or CLDN4.
- Not TTF-1 IHC. NKX2-1 here is RNA.
- Not ICI or TROP2-ADC outcome.
- Not a test of any claim other than A10.
- Not a disagreement with the narrower PR #103 result (NKX2-1 vs TACSTD2/CLDN4
  in TCGA-LUAD). This analysis reproduces that null/positive pattern and adds
  the ELF3 set, GTEx, Census, and perturbation.

---

## Limitations

1. recount3 gene sums are coverage, not strict counts. Rank statistics are
   invariant to per-gene length; absolute CPM is not used for inference.
2. AGER and PTPRC failed to map (all-zero). Composition covariates are the
   remaining lineage markers.
3. Human perturbation n is small (2–4 / arm). We therefore require sign
   consistency, not per-contrast FDR.
4. GSE129340 has n=1 control; MWU is undefined. Log2FC is descriptive.
5. Census within-type labels still mix subtypes and datasets; positive
   within-type ρ can include residual composition and depth. The *absence*
   of an NKX2-1 inverse is the more robust P6 result.
6. Mouse Census has almost no labelled airway basal/club compartment.
7. Cell lines are not primary lung.

---

## Files

- Spec: `claims/claim_A10.md`
- Tables: `results/claim_A10/tables/`
- Figures: `results/claim_A10/figures/fig1`–`fig7`
- Logs: `results/claim_A10/logs/perturbation_notes.txt`
- Scripts: `scripts/01_human_bulk.py` … `scripts/05_write_report.py`
