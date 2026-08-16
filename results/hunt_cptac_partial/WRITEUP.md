# CPTAC LUAD/LSCC: TACSTD2 vs xCell CD8/immune after purity (partial Spearman)

**Slice:** `notes|scripts|results/hunt_cptac_partial/` only.
**Question:** TACSTD2 protein/RNA vs xCell CD8 and xCell immune score, unadjusted and after **any available purity proxy**. User direction: **negative**. Prior LUAD slice: protein vs xCell CD8 **ρ ≈ −0.29**.
**ICI labels:** none. Gillette *Cell* 2020 (LUAD) and Satpathy *Cell* 2021 (LSCC) are treatment-naive surgical cohorts. This slice cannot test immunotherapy resistance.

---

## English

### Verdict (honest)

| Claim | LUAD | LSCC |
|---|---|---|
| TACSTD2 **protein** vs xCell CD8 / immune is negative | **Yes, and it survives DNA purity** | Same sign, **not significant** |
| TACSTD2 **RNA** vs xCell CD8 / immune is negative | Direction only; **not significant** | Direction only; **not significant** |
| Signal is just tumor-purity confounding | **No** for LUAD protein (protein–purity ρ ≈ +0.07) | Cannot claim a signal to confound |

The prior **LUAD protein vs xCell CD8 ρ = −0.289** (n=110, p=0.00219) **replicates exactly** on the same freeze. After `WES_purity` it is **ρ = −0.261** (n=108, p=0.00670); after `WGS_purity` **ρ = −0.291** (n=104, p=0.00285). LSCC does **not** independently confirm (protein vs CD8 ρ = −0.080, p=0.41).

Pooled Fisher-z meta (LUAD+LSCC) stays negative for protein, but that average is **LUAD-driven**. Do not read the meta as a second histology replication.

### Data (open S3 freeze v1.2, HEAD HTTP 200)

Prefix: `https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/{LUAD,LSCC}/`

Index pages: [CPTAC-pancan-LUAD](https://www.linkedomics.org/data_download/CPTAC-pancan-LUAD/), [CPTAC-pancan-LSCC](https://www.linkedomics.org/data_download/CPTAC-pancan-LSCC/).

| Cohort | RNA tumor | Protein tumor | Phenotype | Intersect n |
|---|---|---|---|---|
| LUAD | `LUAD_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt` | `…_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt` | `LUAD_phenotype.txt` | **110** |
| LSCC | same name with `LSCC_` | same | `LSCC_phenotype.txt` | **108** |

TACSTD2 row: `ENSG00000184292.7` in both RNA and protein, both cohorts (110/110 and 108/108 quantified). xCell CD8, xCell immune, CIBERSORT CD8, and ESTIMATE ImmuneScore are complete in phenotype. NAT files were not used.

xCell / CIBERSORT / ESTIMATE were **not recomputed**. They are freeze phenotype scores (Li et al. *Cell Syst.* 2023; immunedeconv). Protein-vs-xCell is cross-layer. RNA-vs-xCell is the **same RNA matrix** and is a weaker test.

### Purity proxies that actually exist

| Proxy | Layer | LUAD n | LSCC n | Use |
|---|---|---|---|---|
| `WES_purity` | DNA | 108 (2 NA) | 107 (1 NA) | **Preferred** |
| `WGS_purity` | DNA | 104 (6 NA) | 104 (4 NA) | Preferred alternate |
| `ESTIMATE_ESTIMATEScore` | RNA impurity | 110 | 108 | Sensitivity only; **not independent of xCell** |
| Yoshihara 2013 `cos(0.605 + 1.47e-4 × ESTIMATEScore)` | derived | — | — | **Rejected** |

Yoshihara cosine is unusable on this freeze: ESTIMATE scores are ~5k–21k, 24/110 LUAD and 30/108 LSCC wrap past π, and **110/110 LUAD “purities” are negative**. QC: `tables/yoshihara_purity_qc.tsv`. No ABSOLUTE / CPE column exists in phenotype. Partial Spearman on `ESTIMATEScore` vs `−ESTIMATEScore` is identical (both pairwise ρ with z flip).

xCell scores **are** purity-associated (this is why partial correlation was required):

| Target vs `WES_purity` | LUAD ρ | LSCC ρ |
|---|---|---|
| xCell CD8 | −0.231 (p=0.016) | −0.396 (p=2.4×10⁻⁵) |
| xCell immune | −0.402 (p=1.6×10⁻⁵) | −0.724 (p=1.3×10⁻¹⁸) |

TACSTD2 itself is **not** a purity surrogate (protein vs WES: LUAD ρ=+0.074, p=0.45; LSCC ρ=+0.066, p=0.50).

### Methods

1. Samples = RNA tumor ∩ protein tumor ∩ phenotype IDs. Pairwise-complete per test.
2. Unadjusted Spearman (`scipy.stats.spearmanr`).
3. Partial Spearman: Pearson on ranks, then  
   ρ<sub>xy·z</sub> = (ρ<sub>xy</sub> − ρ<sub>xz</sub>ρ<sub>yz</sub>) / √[(1−ρ<sub>xz</sub>²)(1−ρ<sub>yz</sub>²)],  
   p from Student-t with n−3 df.
4. Rank-residual Spearman stored as sensitivity (agrees with the formula to ~0.01).
5. BH-FDR on the **8 pre-specified unadjusted tests** (2 cohorts × protein/RNA × CD8/immune). A second FDR is within those 8 pairs **per purity adjust**.
6. Inverse-variance Fisher-z meta of the two histologies — descriptive only.

No ICI labels, no fabricated accessions, no Yoshihara purity used as a covariate after it failed QC.

### Primary results

**Unadjusted (8 pre-specified tests)**

| Cohort | Predictor | Target | n | ρ | p | q (BH among 8) | User (−) |
|---|---|---|---|---|---|---|---|
| LUAD | protein | xCell immune | 110 | **−0.309** | **0.00100** | **0.0080** | yes |
| LUAD | protein | xCell CD8 | 110 | **−0.289** | **0.00219** | **0.0088** | yes |
| LSCC | protein | xCell immune | 108 | −0.133 | 0.169 | 0.450 | yes, ns |
| LSCC | protein | xCell CD8 | 108 | −0.080 | 0.413 | 0.472 | yes, ns |
| LSCC | RNA | xCell immune | 108 | −0.118 | 0.225 | 0.450 | yes, ns |
| LSCC | RNA | xCell CD8 | 108 | −0.100 | 0.302 | 0.472 | yes, ns |
| LUAD | RNA | xCell CD8 | 110 | −0.086 | 0.370 | 0.472 | yes, ns |
| LUAD | RNA | xCell immune | 110 | −0.007 | 0.943 | 0.943 | ~0 |

Only the two LUAD **protein** tests pass FDR. RNA is not a TROP2-protein substitute here.

**LUAD protein after purity (the extension that was asked for)**

| Target | Adjust | n | ρ_partial | p | q (among 8 primary / adjust) |
|---|---|---|---|---|---|
| xCell CD8 | WES_purity | 108 | **−0.261** | **0.00670** | **0.027** |
| xCell CD8 | WGS_purity | 104 | **−0.291** | **0.00285** | **0.021** |
| xCell CD8 | ESTIMATEScore | 110 | −0.254 | 0.00772 | 0.031 |
| xCell immune | WES_purity | 108 | **−0.272** | **0.00466** | **0.027** |
| xCell immune | WGS_purity | 104 | **−0.273** | **0.00521** | **0.021** |
| xCell immune | ESTIMATEScore | 110 | −0.317 | 0.00078 | 0.006 |

DNA purity moves LUAD protein–CD8 from −0.289 to about **−0.26 to −0.29**. Direction and p<0.01 hold. ESTIMATEScore is shown because the prompt said “any” proxy; it is RNA immune+stroma and **over-adjusts** an immune endpoint.

**LSCC protein after purity:** still |ρ|<0.13 and p>0.18 after WES/WGS. After ESTIMATEScore the CD8/immune partials flip to **+0.023 / +0.047** (p>0.63). That flip is over-adjustment, not a positive biology claim.

**Secondary (not in the primary 8):** LUAD protein vs CIBERSORT CD8 ρ=−0.133, p=0.16; vs ESTIMATE ImmuneScore ρ=−0.185, p=0.052. Same direction, weaker / method-dependent. Do not treat CIBERSORT as confirming the xCell CD8 hit.

**Pooled meta (do not over-read):** protein vs CD8 unadjusted ρ=−0.187, p=0.0058; after WES ρ=−0.163, p=0.017. Heterogeneous: LUAD is the signal, LSCC is a null with the same sign.

### What this is not

- Not ICI response. No atezolizumab / pembrolizumab / durvalumab labels in this freeze.
- Not a causal “TROP2 excludes T cells” proof. Residual confounding (histology subtype, STK11, smoking, stroma) remains.
- Not a license to use TACSTD2 **RNA** as the protein phenotype. LUAD RNA–xCell is null.
- Not a successful ESTIMATE tumor-purity transform. The 2013 cosine failed QC here.

### Reproducibility

```bash
python3 scripts/hunt_cptac_partial/00_download.py --outdir data/hunt_cptac_partial
python3 scripts/hunt_cptac_partial/01_analyze.py --data data/hunt_cptac_partial --outdir results/hunt_cptac_partial
```

Raw matrices stay in `data/` (gitignored). Tables and figures are under `results/hunt_cptac_partial/`.

---

## 中文

### 结论（如实）

用户方向是**负相关**。LUAD 里 **TACSTD2 蛋白** vs xCell CD8 **ρ = −0.289**（n=110，p=0.00219）与此前切片一致；对 `WES_purity` 做偏相关后为 **ρ = −0.261**（n=108，p=0.00670），对 `WGS_purity` 为 **ρ = −0.291**（n=104，p=0.00285）。**不是纯度假象**：蛋白与 WES 纯度只有 ρ=+0.074。

**LSCC 不能独立重复**：蛋白 vs CD8 ρ=−0.080，p=0.41；vs immune ρ=−0.133，p=0.17。符号同向，但都不显著。两癌种 Fisher-z 合并仍为负，**主要由 LUAD 拉动**，不能当成鳞癌复现。

**TACSTD2 RNA** 在两个癌种对 xCell CD8/immune 都接近零（|ρ|≤0.12，p>0.22）。RNA 不能代替蛋白。

Yoshihara 2013 的 ESTIMATE 纯度余弦在本 freeze **不可用**（分数约 5k–21k，LUAD 110/110 算出负纯度，且部分样本越过 π）。RNA 侧只把 `ESTIMATEScore` 当作杂质轴做敏感性；它与 xCell 不独立，LSCC 上偏相关甚至会翻到接近 0 / 弱正，**不能解释为生物学反转**。

本队列是手术、治疗前标本，**没有 ICI 疗效标签**，不能检验免疫治疗耐药。
