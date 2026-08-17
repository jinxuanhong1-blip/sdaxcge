# Playbook — NicheNet prior, CLDN4-only malignant senders → T/NK IFN / cytotoxicity

> **Public processed files only.** No EGA/GSA FASTQ. No invented statistics.
> English first; 中文 at [中文版](#中文版).

---

## 0. Why this module exists

The claim to test is mechanistic and **CLDN4-only**:

> Ligands from **CLDN4-high malignant** cells predict **T/NK** gene changes (**IFN** or **cytotoxicity**), and that sender→receiver map can be ranked with a public NicheNet-v2 prior.

That is a **NicheNet-class** estimand (Browaeys et al., *Nat Methods* 2020): rank ligands by how well a **public ligand–target prior** explains a receiver gene set.

This folder is **additive**. It does not replace `methods/scrna_nichenet/` (TACSTD2∩CLDN4 dual-high) and does not use TACSTD2 to call senders.

---

## 1. Which cohort (and which one was skipped)

| Series | Why it qualifies | Why it may fail | Decision here |
| --- | --- | --- | --- |
| [GSE207422](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207422) | Neoadjuvant PD-1 + chemo NSCLC. **MPR n=4** (pCR counted as MPR) vs **NMPR n=8**. GEO UMI matrix 175 MB + sample sheet. Hu et al., *Genome Med* 2023, PMID 36869384. | Author CopyKAT malignant IDs were never deposited. | **Used.** Malignant / T / NK from **DRMref** public labels (Liu et al., *NAR* 2024). |
| [GSE253013](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253013) | Public LUAD scRNA RDS. | **9.3 GB.** No MPR/NMPR / ICI response column in the release used by sibling audits. | **Not downloaded.** Wrong phenotype for combinatorial MPR vs NMPR. |

Do not silently swap RECIST for MPR. P06 `pCR` is MPR.

---

## 2. Estimands (write these down before touching the matrix)

| ID | Estimand | Unit | Contrast |
| --- | --- | --- | --- |
| **E1** | T/NK **IFN** and **cytotoxicity** scores | patient | NMPR vs MPR (n=8 vs 4) |
| **E2** | Same scores | patient | vs malignant **CLDN4** (n=12) |
| **E3** | NicheNet **ligand activity** (Pearson / AUROC of prior vs gene-set membership) for ligands expressed in **CLDN4-high** malignant cells whose receptors are expressed on T/NK | ligand rank | a priori IFN; a priori cytotoxicity |
| **E4** | Combinatorial | ligand rank | **NMPR** CLDN4-high senders → **NMPR** T/NK vs **MPR** CLDN4-high senders → **MPR** T/NK |

**Stop rules**

- Cells are **not** replicates of MPR or of “CLDN4-high tumor”.
- NicheNet’s prior is **unsigned regulatory potential**. It does **not** encode “this ligand represses GZMB” or “this ligand induces ISGs”. Direction of IFN / cytotoxicity is **E1/E2 only**.
- If E1/E2 are null, E3/E4 are still reportable as a prior ranking, but they do **not** prove a signed claim.
- Do not AND TACSTD2 into the sender call in this folder.

---

## 3. Public files

| Role | File | Use |
| --- | --- | --- |
| UMI | `GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz` | gene × 92,330 barcodes |
| Sample sheet | `GSE207422_NSCLC_scRNAseq_metadata.xlsx` / committed TSV | patient, MPR/NMPR, timing |
| Cell types | DRMref `GSE207422_Tor` + `GSE207422_Sin` `meta.data` (30,877 post-tx barcodes) | `Malignant cells`, `CD8+ T cells`, `CD4+ T cells`, `NK cells` |
| Prior | Zenodo [7074291](https://doi.org/10.5281/zenodo.7074291) `ligand_target_matrix_nsga2r_final.rds` (targets × ligands) | ligand activity |
| LR | same record, `lr_network_human_21122021.rds` | potential ligands = expressed ligand ∩ expressed receptor |

Full `nichenetr` / `nichenetpy` is optional. This repo implements the **geneset ligand-activity** step in Python on the published v2 matrix (R is not required).

---

## 4. Definitions (pre-registered)

**CLDN4-high malignant cell:** DRMref `Malignant cells` with `log1p(CP10k)` **CLDN4 ≥ malignant median**. Library size = DRMref `nCount_RNA`. TACSTD2 is ignored for this call.

**T/NK:** DRMref `CD8+ T cells` + `CD4+ T cells` + `NK cells`.

**Cytotoxicity (a priori):** GZMB GZMA GZMH GZMK GZMM PRF1 GNLY NKG7 IFNG FASLG TNF CST7 FGFBP2 KLRK1.

**IFN (a priori):** IFNG STAT1 IRF1 IRF7 ISG15 MX1 IFIT1 IFIT2 IFIT3 OAS1 IFI6 RSAD2 CXCL9 CXCL10 CXCL11 GBP1.

**Expressed gene (10x-style):** detected in ≥ 10% of cells in that compartment (NicheNet vignette default for 10x).

**Potential ligand:** ligand expressed in the sender set, with at least one receptor expressed in the receiver set, and present as a column in the v2 matrix.

**Background genes:** intersection of (T/NK-expressed panel genes) with prior target rows. This is **narrower** than a full-transcriptome NicheNet run. Report `n_background`.

**Empirical T/NK gene sets:** patient-mean `log1p(CP10k)` per gene; Wilcoxon NMPR vs MPR; retain genes with p < 0.15 (honest: n=4 vs 8 will not yield FDR < 0.05). Split by sign of NMPR − MPR.

**Ligand activity:** for each ligand, Pearson / AUROC / AUPR between the prior target scores and gene-set membership on the background. Rank by Pearson (NicheNet default).

---

## 5. What not to do

- Download GSE253013 “to have a second cohort” when it has no MPR label.
- Rank ligands on cell-level DE (pseudoreplication).
- Drop P06 (15 malignant cells) after seeing p-values.
- Treat a high Pearson on the IFN list as evidence that IFN is up (or down) in NMPR.
- Require TACSTD2-high in this folder.
- Commit the 175 MB matrix or the 262 MB RDS.

---

## 6. Run

```bash
pip install -r methods/scrna_nichenet_cldn4/env/requirements.txt
python3 methods/scrna_nichenet_cldn4/scripts/00_download.py
python3 methods/scrna_nichenet_cldn4/scripts/01_convert_prior.py
python3 methods/scrna_nichenet_cldn4/scripts/02_extract_panel.py
python3 methods/scrna_nichenet_cldn4/scripts/03_analyze.py
```

Outputs: `results/n_table.tsv`, `top_ligands.tsv`, `top_ligands_primary.tsv`, `patient_level_tests.tsv`, `FINDING.md`.

---

## 中文版

要检验的是：**CLDN4 高恶性细胞**（只用 CLDN4，不用 TACSTD2 定义发送端）上的配体，能否用公开 **NicheNet-v2 配体–靶基因先验**解释 **T/NK IFN 或杀伤**基因集，并按 **MPR vs NMPR** 接收端拆开。

- 队列用 **GSE207422**（术后 MPR 4 例含 pCR，NMPR 8 例）。**GSE253013** 不下载：9.3 GB，且无 MPR/NMPR。
- 恶性 / T / NK 用 **DRMref 公开标签**，不是作者未上传的 CopyKAT。
- **推断单位是患者**，不是细胞。
- NicheNet 先验是**无符号**调控潜能，不能单独证明 IFN / 杀伤升降；方向只看患者水平 E1/E2。
- CLDN4 高恶性细胞：CLDN4 的 `log1p(CP10k)` **≥ 恶性细胞中位数**。
- 配体活性：先验分数 vs 基因集成员的 Pearson / AUROC；按 Pearson 排序。
- 经验基因集用患者均值 Wilcoxon，p < 0.15（n=4 vs 8，不要假装 FDR<0.05）。
