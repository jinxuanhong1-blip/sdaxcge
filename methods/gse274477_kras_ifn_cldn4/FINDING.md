# FINDING: Cldn4 vs impaired type-I IFN in KRAS / KL epithelium

**Question.** Does Cldn4 track the impaired type-I IFN state in KRAS / KL mouse lung epithelium (Cldn4-only)? High Cldn4 vs low: IFN down, MHC down? KL vs K/KP Cldn4 if GSE274351 has a public gene matrix?

**Answer.** No. On the public processed matrices, high Cldn4 does **not** mark the IFN-low / MHC-low state. KRAS-mutant AT2 cells as a class are IFN-down and MHC-down vs WT (paper recovered). Cldn4 is rare. Among KRAS-mutant AT2, the Cldn4+ minority has **higher** IFN_core than Cldn4− cells (p=0.0045), not lower. MHC is not lower. LCM adenomas (n=14) show the same direction and are not significant. KL Cldn4 is not higher than K/KP (n_KL=5 vs n_K+KP=9, p=0.69).

Thesis unchanged (Cldn4-only). This is additive public mouse evidence against a simple Cldn4 = IFN-impaired rule in this Fernández-García KRAS / KL epithelium.

## Data (public processed only)

| accession | deposited file | assay | labels | honest n |
| --- | --- | --- | --- | --- |
| [GSE274477](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274477) | `GSE274477_Feature_counts_Matrix_14Aug.tsv.gz` | Fluidigm C1 scRNA-seq feature counts. AT2, 4 wk after Ad-Cre. KRAS WT vs KRAS(G12V) | `CONTROL` = WT, `MUT` = KRAS(G12V); chips `FCNIO1` / `FCNIO2` | Matrix: 1640 wells (820 WT, 820 MUT, 40 named Undetermined). Paper after author QC: 139 WT + 176 MUT (n=315). **This QC: 189 WT + 268 MUT (n=457)** |
| [GSE274351](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274351) | `GSE274351_expredata_TPM_gene.txt.gz` | LCM small adenomas + adjacent alveolus, QuantSeq 3′ TPM | **TPM column labels** (not GEO titles) | **K=4, KP=5, KL=5, NL=4** (18 columns) |

PMID 39186651 (Fernández-García et al., *PNAS* 2024). *Mus musculus*. No SRA / FASTQ / realignment.

**GSE274351 label mismatch.** GEO sample titles: 5 K, 5 KP, 4 KL, 4 NL. Deposited TPM columns: `K1–K4`, `KP1–KP5`, `KL1–KL5`, `NL1–NL4`. Genotype n and contrasts use the **file**. That is why n is 4/5/5/4, not the GEO title counts.

## Genes

Cldn4-only. GRCm38 `ENSMUSG00000047501` → `Cldn4`. Present on both matrices.

**IFN_core** (paper Fig. 2B / 3A): Stat1, Irf7, Ifih1, Ddx58 (Rigi), Oasl2, Bst2, Ifi27l2a, Ifitm3. All 8 present on both matrices.

**MHC / antigen presentation** (paper Fig. 2B): B2m, Tap1, Tapbp, H2-K1, H2-D1, H2-T23, H2-Q6, H2-Q7, H2-Aa. All 9 present on both matrices.

**H_IFNA** (secondary): MSigDB mouse `HALLMARK_INTERFERON_ALPHA_RESPONSE` (2024.1.Mm), 94 genes; 93 recovered on each matrix (`Wars1` absent).

Done condition: genes and labels are present; the Cldn4–IFN table is Table 2 (also `cldn4_ifn_table.tsv`).

## Methods

1. GEO supplementary processed files only (FTP `GSE274nnn`).
2. `ENSMUSG` IDs version-stripped and mapped with Ensembl GRCm38.102 GTF.
3. **GSE274477 QC.** Author well list is not deposited. Drop `Undetermined`. Keep library size ≥ 50,000, 1,000–8,000 detected genes, mitochondrial fraction ≤ 10% (MT genes from the same GTF). This is a depth filter near the paper’s reported averages (~100,000 reads, ~2,354 genes), not a reproduction of n=315.
4. scRNA: log1p(10,000 × count / library size). Gene-set score = mean of per-gene z-scores across QC cells.
5. High vs low Cldn4 = median split. In KRAS-mutant AT2 the median is 0, so the split is **Cldn4+ (counts > 0) vs Cldn4−**. Primary: within MUT. Also all QC AT2.
6. **GSE274351:** deposited TPM. Scores = mean z of log2(TPM+1) across the 18 samples. Adenoma-only median split (NL out), n=14. KL vs K/KP uses matrix labels.
7. Two-sided Mann–Whitney U; Spearman. Rank-biserial r > 0 means the first named group is higher. No imputation. No dropped labeled samples.

## Table 1. KRAS MUT vs WT AT2 (GSE274477 QC)

Sanity check: deposited counts recover the paper’s IFN / MHC impairment.

| score | n WT | n MUT | median WT | median MUT | Δ MUT−WT | MWU p | rank-biserial (MUT>WT) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Cldn4 | 189 | 268 | 0.000 | 0.000 | 0.000 | 0.0458 | 0.051 |
| IFN_core | 189 | 268 | 0.066 | −0.310 | −0.375 | 5.54e-17 | −0.460 |
| MHC | 189 | 268 | 0.169 | −0.175 | −0.344 | 8.31e-11 | −0.357 |
| H_IFNA | 189 | 268 | −0.004 | −0.090 | −0.086 | 2.51e-07 | −0.283 |

Cldn4 detection (counts > 0): WT 9/189 (4.8%); MUT 26/268 (9.7%). Both medians are 0. The small positive rank-biserial is higher MUT detection, not a shift in the median.

## Table 2. Cldn4–IFN / MHC

High = Cldn4 above median. In MUT AT2 that is Cldn4+ n=26 vs Cldn4− n=242.

| comparison | score | n high | n low | median high | median low | Δ high−low | MWU p | rank-biserial (high>low) | Spearman Cldn4 vs score |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE274477 MUT AT2 | IFN_core | 26 | 242 | −0.136 | −0.324 | **+0.187** | **0.0045** | +0.339 | +0.179 (p=0.0033) |
| GSE274477 MUT AT2 | MHC | 26 | 242 | −0.132 | −0.176 | +0.045 | 0.7524 | +0.038 | +0.019 (p=0.7617) |
| GSE274477 MUT AT2 | H_IFNA | 26 | 242 | −0.043 | −0.099 | +0.056 | 0.0782 | +0.210 | +0.115 (p=0.0607) |
| GSE274477 all QC AT2 | IFN_core | 35 | 422 | −0.088 | −0.192 | +0.104 | 0.0768 | +0.180 | +0.085 (p=0.0705) |
| GSE274477 all QC AT2 | MHC | 35 | 422 | 0.079 | −0.023 | +0.102 | 0.5213 | +0.065 | +0.029 (p=0.5347) |
| GSE274351 LCM adenomas | IFN_core | 7 | 7 | −0.183 | −0.386 | +0.203 | 0.9015 | +0.061 | +0.082 (p=0.7799) |
| GSE274351 LCM adenomas | MHC | 7 | 7 | −0.161 | −0.209 | +0.048 | 0.7104 | +0.143 | +0.149 (p=0.6114) |
| GSE274351 LCM adenomas | H_IFNA | 7 | 7 | −0.257 | −0.348 | +0.091 | 0.7104 | +0.143 | +0.200 (p=0.4930) |

Loose QC (lib ≥ 20k, genes ≥ 800, mito ≤ 15%; MUT n=404, Cldn4+ n=34): IFN_core Δ high−low = +0.196, p=0.0008; MHC Δ = +0.072, p=0.57. Same direction.

LCM high-Cldn4 adenomas: K1, K2, KL1, KL4, KP2, KP3, KP5. Low: K3, K4, KL2, KL3, KL5, KP1, KP4.

## Table 3. LCM Cldn4 by genotype (matrix labels)

| genotype | n | samples | Cldn4 TPM median (IQR) | IFN_core z median | MHC z median |
| --- | --- | --- | --- | --- | --- |
| NL | 4 | NL1–NL4 | 0.00 (0.00–0.00) | 0.868 | 0.571 |
| K | 4 | K1–K4 | 5.06 (0.08–14.66) | −0.366 | −0.023 |
| KP | 5 | KP1–KP5 | 22.48 (2.68–25.44) | −0.152 | 0.474 |
| KL | 5 | KL1–KL5 | 5.06 (0.00–8.32) | −0.525 | −0.452 |

NL IFN / MHC are high and Cldn4 is 0 in all four adjacent alveolus samples. Adenomas are IFN-lower than NL as a class (see per-sample table).

## Table 4. KL vs K / KP Cldn4 (GSE274351)

| contrast | n a | n b | Cldn4 TPM median a | Cldn4 TPM median b | MWU p | rank-biserial (a>b) |
| --- | --- | --- | --- | --- | --- | --- |
| KL vs K | 5 | 4 | 5.06 | 5.06 | 0.9009 | −0.100 |
| KL vs KP | 5 | 5 | 5.06 | 22.48 | 0.6723 | −0.200 |
| KL vs K+KP | 5 | 9 | 5.06 | 10.00 | 0.6859 | −0.156 |
| KL vs NL | 5 | 4 | 5.06 | 0.00 | 0.1094 | +0.600 |
| K vs NL | 4 | 4 | 5.06 | 0.00 | 0.0689 | +0.750 |
| KP vs NL | 5 | 4 | 22.48 | 0.00 | 0.0442 | +0.800 |

KL is not Cldn4-high vs K or KP. Point estimate is lower than KP (5 vs 22 TPM) and not significant at this n.

## Table 5. Per-sample LCM

| sample | genotype | Cldn4 TPM | IFN_core z | MHC z | Hallmark IFNα z |
| --- | --- | --- | --- | --- | --- |
| K1_S8_L00 | K | 10.00 | −0.183 | 0.229 | −0.209 |
| K2_S28_L00 | K | 28.62 | −0.859 | −0.607 | −0.282 |
| K3_S9_L00 | K | 0.00 | −0.386 | −0.209 | −0.348 |
| K4_S10_L00 | K | 0.11 | −0.347 | 0.164 | −0.001 |
| KL1_S16_L00 | KL | 8.32 | 0.334 | −0.161 | −0.011 |
| KL2_S17_L00 | KL | 5.06 | 0.136 | 0.032 | 0.235 |
| KL3_S29_L00 | KL | 0.00 | −0.525 | −0.452 | −0.379 |
| KL4_S30_L00 | KL | 38.53 | −0.620 | −0.525 | −0.500 |
| KL5_S31_L00 | KL | 0.00 | −0.628 | −1.477 | −0.365 |
| KP1_S11_L00 | KP | 2.68 | 0.820 | 0.975 | 0.388 |
| KP2_S12_L00 | KP | 22.48 | 0.698 | 0.474 | 0.442 |
| KP3_S13_L00 | KP | 25.44 | −0.617 | 0.792 | −0.316 |
| KP4_S14_L00 | KP | 0.00 | −0.850 | −1.184 | −0.761 |
| KP5_S15_L00 | KP | 34.40 | −0.152 | −0.558 | −0.257 |
| NL1_S18_L00 | NL | 0.00 | 0.360 | 0.616 | 0.504 |
| NL2_S19_L00 | NL | 0.00 | 0.877 | 0.527 | 0.600 |
| NL3_S20_L00 | NL | 0.00 | 0.859 | 0.891 | 0.669 |
| NL4_S21_L00 | NL | 0.00 | 1.083 | 0.475 | 0.590 |

## Interpretation

1. **Paper state is in the file.** KRAS(G12V) AT2: IFN_core down, MHC down, Hallmark IFNα down vs WT. LCM adenomas sit below adjacent NL on IFN / MHC. The impaired type-I IFN / MHC program is real in these deposits.
2. **Cldn4 does not track that impaired state inside KRAS epithelium.** Predicted if it did: high Cldn4 → IFN down, MHC down. Observed: Cldn4+ MUT AT2 → IFN_core **up** vs Cldn4− (p=0.0045); MHC not down. LCM high-Cldn4 adenomas are not IFN/MHC-low (n=7 vs 7, n.s.).
3. **Cldn4 is a sparse AT2 transcript here.** 26/268 MUT and 9/189 WT QC cells. High/low is detect vs not-detect. Cldn4+ MUT wells have slightly higher library size (median 194k vs 171k, p=0.060) and more genes (3826 vs 3277). Depth may contribute; it does not flip the sign to IFN-down.
4. **KL is not a Cldn4-high genotype on this matrix.** KL median TPM 5.06 vs K 5.06 vs KP 22.48. KL vs K+KP p=0.69. n is 5/4/5. Honest: no KL Cldn4 elevation is shown.

## Caveats

- Public processed only. Author C1 barcodes after QC are not in GEO; n here is 457, not 315.
- GSE274351 GEO titles ≠ TPM column labels.
- LCM is bulk adenoma, not pure AT2. Stroma, immune cells can move IFN / MHC.
- Fluidigm C1 dropout: Cldn4 absence is not proof of no protein.
- Hallmark IFNα is secondary. Paper-core lists are primary.
- No other claudins. No human data. No raw-read reprocessing.
- LCM n=18. p-values are descriptive.

## Files

- `analyze.py` — download-once, public-matrix analysis
- `results.json` — numbers
- `cldn4_ifn_table.tsv` — Table 2
- `FINDING.md` — this note
