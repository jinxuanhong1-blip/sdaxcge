# TROP2 / TACSTD2 loss-of-function: CLDN4, junctions and IFN/immune programs
# TROP2 / TACSTD2 功能缺失：CLDN4、细胞连接与干扰素/免疫程序

**Question / 问题:** Across public TACSTD2/TROP2 knockdown / knockout / shRNA / siRNA /
CRISPR datasets, does TROP2 loss change **CLDN4**, **cell junctions**, and
**IFN / immune genes**?
在公开的 TACSTD2/TROP2 敲低/敲除/shRNA/siRNA/CRISPR 数据集中，TROP2 缺失是否改变
**CLDN4**、**细胞连接** 与 **干扰素/免疫基因**？

**Short answer / 简答:** Yes, reproducibly, but with **context-dependent direction**.
The single most reproducible effect is **cell-intrinsic de-repression of
interferon-stimulated genes (ISGs)** upon TROP2 (or its desmosomal partner DSG2) loss.
Junctions (tight-junction / desmosome) tend to **weaken**, and CLDN4 goes **down** where
it changes significantly.
是的，且可重复，但**方向依语境而定**。最可重复的效应是 TROP2（或其桥粒伙伴 DSG2）
缺失后**肿瘤细胞内在的干扰素刺激基因（ISG）去抑制（上调）**；细胞连接（紧密连接/桥粒）
倾向于**减弱**；CLDN4 在显著变化时呈**下调**。

---

## 1. Datasets (exhaustive, verified) / 数据集（穷尽检索、已核实）

Search: NCBI GEO DataSets (E-utilities) + EBI ArrayExpress/BioStudies, 2026-08-16.
Field-restricted title/summary queries + sample-tag probes (`shTACSTD2`, `siTACSTD2`,
`Trop2 knockout`, …). Bare `TACSTD2`/`TROP2` (~33k/~2k hits) and broad `TACSTD2 shRNA`
(129) hits were CRISPR/RNAi library screens or gene-list matches, **not** perturbation
experiments (verified by reading each esummary). ArrayExpress added only TROP2⁺
**cell-sorting** (marker) studies. Details: `notes/fable_tacstd2_kdko/01_dataset_inventory.md`.

All four true loss-of-function transcriptomic datasets have **processed data <2 GB**
(raw FASTQ / >2 GB skipped):

| Accession | Verified | Organism / model | Perturbation | Design | Processed file |
|---|---|---|---|---|---|
| [GSE334497](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE334497) | esummary + FTP | mouse 4T1 TNBC tumour | Trop2 **CRISPR KO** vs WT | 5 vs 5 | `GSE334497_normalized_counts.csv.gz` (1.2 MB) |
| [GSE289287](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289287) | esummary + FTP | human T-47D breast, xenograft | Trop-2 **KO** vs WT (also DSG2 KO) | 4 vs 3 | author DESeq2 tables (~2 MB) |
| [GSE245459](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245459) | esummary + FTP | human ovarian cancer cells | sh**TACSTD2** vs shNC | 3 vs 3 | `GSE245459_fpkm.anno.txt.gz` (14 MB) |
| [GSE15212](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15212) | esummary + GSM tables | human colorectal RKO cells | **TACSTD2** siRNA vs neg-ctrl | 6 vs 12 | per-sample Agilent GPL4133 |

GSE289287 also ships **author DESeq2** results for **DSG2 KO** (the desmosomal partner of
TROP2); we analyse these as a mechanistic comparison.

## 2. Methods / 方法

- Contrast is always **perturbation vs control**; `log2FC > 0` = **up on TROP2/TACSTD2 loss**.
- GSE289287: use the **author DESeq2** (log2FC, p, padj) directly — the most rigorous.
- GSE334497 / GSE245459 / GSE15212: **Welch t-test + Benjamini-Hochberg FDR** and
  Mann-Whitney U on log-scale expression (log2 norm-counts+0.5 / log2 FPKM+1 /
  quantile-normalized log2 signal), with Cohen's d.
- Gene sets: a simplified **CAMERA-style competitive** test (Mann-Whitney of the per-gene
  signed statistic, set vs rest of genome) **plus** a **self-contained** Wilcoxon that the
  set's log2FC is shifted from zero. Panels in `scripts/fable_tacstd2_kdko/gene_panels.py`.
- Perturbation QC: TACSTD2/TROP2 itself must drop.
- Full pipeline: `scripts/fable_tacstd2_kdko/01…07`. Reproduce with the commands in §6.

## 3. Perturbation worked / 扰动确实生效

TACSTD2 is strongly reduced in every loss dataset (log2FC −2.6 to −3.9; significant in
all), and unchanged in the DSG2-KO specificity controls. 每个缺失数据集中 TACSTD2 均显著
下调（log2FC −2.6 ~ −3.9），DSG2 敲除对照中 TACSTD2 不变。
(`SYNTHESIS_perturbation_qc.tsv`)

## 4. Results by axis / 分轴结果

Competitive-test direction; `***`p<1e-3, `**`p<1e-2, `*`p<0.05 (`SYNTHESIS_direction_grid.tsv`):

| Contrast | junctions | tight_jn | desmosome | ifn_immune | ifn_isg | antigen_pres | tcell_cytotox |
|---|---|---|---|---|---|---|---|
| mouse 4T1 Trop2 KO | down*** | down*** | down*** | **up*** | up** | up*** | up*** |
| human T-47D Trop-2 KO | up ns | up ns | up ns | **up*** | up*** | up ns | down ns |
| human ovarian shTACSTD2 | down*** | down*** | down*** | **down*** | down*** | down*** | up ns |
| human CRC siTACSTD2 | down ns | down ns | down (sc*) | **up** | up*** | down ns | down ns |
| DSG2 KO tumour (partner) | down ns | down* | up ns | **up*** | up*** | up*** | up ns |
| DSG2 KO cells (partner) | up ns | down ns | down ns | **up*** | up*** | up** | down ns |

**CLDN4** (`SYNTHESIS_cldn4.tsv`): significantly **down** in ovarian KD
(log2FC −1.92, padj 3.7e-3) and with DSG2 loss (−0.91, padj 2.1e-4); down trend in mouse
4T1 (−0.82, ns); slight ns up in the T-47D and CRC contrasts. So CLDN4 is **not**
uniformly regulated, but when it changes significantly it **decreases** — consistent with
TROP2 supporting, not repressing, the claudin program.

**Junctions / 连接:** weakened (down) in mouse breast and human ovarian (both p<1e-7),
down-trending in colorectal; the desmosome arm is the most consistently reduced. The
T-47D Trop-2 KO xenograft is the outlier (slight ns up).

**IFN / immune / 干扰素-免疫:** **up** in 3/4 loss datasets (mouse 4T1 p=5e-12; human
T-47D p=1.9e-8; human CRC p=6e-3) and in **both** DSG2-KO controls (ISG induction up to
p=1e-27). Because the T-47D xenografts grow in immunodeficient hosts, this ISG induction
is **tumour-cell-intrinsic**. In the immunocompetent 4T1 model it broadens to antigen
presentation + T-cell/cytotoxicity signatures — the transcriptomic correlate of the
"TROP2 → immune exclusion" phenotype. **Ovarian is the sole exception** (ISGs down).

## 5. Conclusion / 结论

TROP2/TACSTD2 loss reproducibly perturbs all three axes:
1. **IFN/immune** — the strongest, most reproducible effect: cell-intrinsic ISG
   de-repression on TROP2 (and DSG2) loss, extending to immune infiltration in an
   immunocompetent model; ovarian is the exception.
2. **Junctions** — generally weakened, especially the desmosome/claudin arm.
3. **CLDN4** — context-dependent; significantly down in ovarian KD and with DSG2 loss.

TROP2 缺失可重复地扰动全部三个轴：最强且最可重复的是**细胞内在 ISG 去抑制（上调）**；
连接（尤其桥粒/claudin）总体**减弱**；CLDN4 依语境显著**下调**（卵巢 KD 与 DSG2 缺失）。

## 6. Reproduce / 复现

```bash
pip install --user pandas numpy scipy statsmodels GEOparse requests   # + network to NCBI/mygene
python scripts/fable_tacstd2_kdko/01_search_geo.py       # exhaustive search
bash   scripts/fable_tacstd2_kdko/02_download.sh         # processed files (<2GB)
python scripts/fable_tacstd2_kdko/03_analyze_gse334497.py
python scripts/fable_tacstd2_kdko/04_analyze_gse289287.py
python scripts/fable_tacstd2_kdko/05_analyze_gse245459.py
python scripts/fable_tacstd2_kdko/06_analyze_gse15212.py
python scripts/fable_tacstd2_kdko/07_synthesize.py       # SYNTHESIS_* tables
```

## 7. Limitations / 局限

Small n in several contrasts (3/3, 4/3); GSE334497 has only normalized (not raw) counts
so its test is Welch/BH rather than DESeq2; competitive gene-set p-values use a simplified
Mann-Whitney that ignores inter-gene correlation (treat as direction evidence,
corroborated by the self-contained test and per-gene tables); mouse↔human orthology is by
symbol + mygene. See `notes/fable_tacstd2_kdko/02_findings.md`.
