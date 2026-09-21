# GSE245459 SKOV3 TACSTD2 shRNA — IFN/APM direction

**This is TACSTD2/TROP2 shRNA in SKOV3, not a CLDN4 knockdown.** Ovarian line, author FPKM, n = 3. The thesis wants IFN/APM up after barrier-gene loss. The search below asks where that increase exists.

---

## Search for an IFN/APM increase

Script: `scripts/gse245459_specificity/sweep_thesis_up.py`.
Eight transforms (log2(FPKM+1), pseudocount 0.1 and 0.01, positive-only log2, upper-quartile log2, CLR, housekeeping-centered log2, within-sample rank). Three thesis estimands (untreated knockdown, knockdown on cisplatin, their interaction). Thirteen IFN/APM sets (C4 panel, APM and its HLA / peptide-loading / immunoproteasome splits, a custom ISG list, Hallmark IFN-α, Hallmark IFN-γ, their overlap, IFN-α with APM genes removed, CXCL9/10/11). That is 312 IFN/APM tests. Wilcoxon signed-rank p-values for an increase were BH-adjusted across the 264 tests with at least 5 genes. The search was run after the untreated decrease was already known, so a small p here is the best hit inside that search, not a pre-registered confirmation.

Context sets (E2F, G2/M, EMT, apoptosis, TNFα, IL-6, inflammatory, allograft rejection, epithelial identity, housekeeping) were scored on the same contrasts and were not entered in the IFN/APM FDR.

### Untreated knockdown — no IFN/APM increase

**0 of 104** untreated IFN/APM tests have a positive median effect. On the reference transform log2(FPKM+1):

| Set | n | median log2FC | up/down | one-sided Wilcoxon p for an increase |
|---|---:|---:|---|---:|
| C4 panel | 6 | −0.571 | 0/6 | 1 |
| APM | 18 | −0.448 | 3/15 | 0.998 |
| APM without HLA-A | 17 | −0.402 | 3/14 | 0.997 |
| Custom IFN/ISG | 44 | −0.531 | 2/42 | 1 |
| Hallmark IFN-α | 83 | −0.542 | 11/72 | 1 |
| Hallmark IFN-γ | 146 | −0.465 | 33/113 | 1 |

Sample-level one-sided p-values for an increase are > 0.998. BH q for an increase is 1. Leave-one-out of each untreated library keeps every core median negative. Dropping, together, the scramble sample with the highest IFN score and the knockdown sample with the lowest IFN score still leaves Hallmark IFN-α at −0.52 and the custom IFN set at −0.50. Within-group Pearson r on genes with mean FPKM ≥ 1 is 0.986–0.990. Replicates sit together on PC1/PC2.

Held-out co-regulation does not rescue the sign. IFN/APM genes with Spearman ρ ≤ −0.4 or ≥ 0.4 versus CLDN4, or versus TACSTD2, on the six cisplatin libraries are still down when the knockdown is tested on the untreated libraries (anti-CLDN4 median −0.53, n = 110; co-CLDN4 median −0.25, n = 58).

CD274 is down (−0.10). STAT1 −0.49, IRF1 −0.29, HLA-A −4.17, IFIT1 −1.89, ISG15 −0.62, B2M −0.51.

**Untreated call.** There is no TACSTD2-shRNA normalization, gene set, sample drop, or CLDN4/TACSTD2 co-regulation subset in this search that puts IFN/APM up. The untreated decrease stands. What that decrease is, relative to the rest of the transcriptome, is in the section below.

### Knockdown on cisplatin — the IFN/APM increase

Estimand: shTACSTD2+DDP versus shNC+DDP (cisplatin 10 µg/ml, 48 h). Both sides are TACSTD2 shRNA versus scramble. CLDN4 is not the perturbed gene, and on this arm it is not further reduced (PR #97: log2FC −0.14, p = 0.099).

Every one of the 104 knockdown-on-cisplatin IFN/APM tests has a positive median, in all eight transforms. Reference transform:

| Set | n | median log2FC | up/down | Wilcoxon p (up) | sample p (up) | BH q inside the 264-test search |
|---|---:|---:|---|---:|---:|---:|
| Hallmark IFN-α | 79 | **+0.406** | 60/19 | 3.7×10⁻⁷ | 0.0031 | **2.3×10⁻⁶** |
| Hallmark IFN-α without APM genes | 72 | +0.393 | 54/18 | 1.6×10⁻⁶ | 0.0034 | 8.7×10⁻⁶ |
| Custom IFN/ISG | 37 | **+0.504** | 30/7 | 4.9×10⁻⁶ | 0.0053 | 2.4×10⁻⁵ |
| Hallmark IFN-α ∩ IFN-γ | 57 | +0.356 | 41/16 | 9.6×10⁻⁵ | 0.0029 | 3.8×10⁻⁴ |
| Hallmark IFN-γ | 134 | +0.188 | 80/54 | 7.1×10⁻⁴ | 0.0044 | 0.0023 |
| APM without HLA-A | 17 | +0.397 | 14/3 | 0.0064 | 0.0070 | 0.017 |
| APM | 18 | +0.397 | 14/4 | 0.033 | 0.017 | 0.060 |
| C4 six-gene panel | 4 | +0.410 | 3/1 | fewer than 5 genes | 0.071 | — |

The C4 genes on this arm: IFIT1 +2.29, ISG15 +0.50, IFI27 +0.32, MX1 +0.24, OAS2 +0.009, HLA-A **−2.11**. HLA-A stays down, which is why the six-gene score is not a sample-level increase (p = 0.071) and why full APM loses the search FDR (q = 0.060) until HLA-A is removed. CD274 is down (−0.22). STAT1 is +0.09. IRF1 is −0.03. NLRC5 is +0.95. B2M is +0.44.

On the same contrast the context sets do not rise: E2F median −0.30, G2/M −0.21, TNFα −0.45, epithelial identity −0.16, housekeeping −0.14, EMT +0.15 (Wilcoxon p for an increase 0.19). Leaving out any one cisplatin library keeps Hallmark IFN-α, the custom IFN set, and APM positive (IFN-α medians +0.36 to +0.45).

**Cisplatin-arm call.** The strongest thesis-aligned IFN/APM increase in the reference normalization is Hallmark IFN-α, median log2FC +0.41, 60/79 genes up, search q = 2.3×10⁻⁶. It is a TACSTD2-shRNA effect on a cisplatin background. It is not an untreated effect and it is not a CLDN4 knockdown.

### Interaction

(KD effect on cisplatin) − (KD effect off cisplatin) is positive for the IFN sets (Hallmark IFN-α median +0.90 on log2(FPKM+1)). That number is large because the untreated effect is negative. EMT’s interaction median is +0.87 and epithelial identity’s is +1.03, while E2F’s is −1.20. The interaction is the arithmetic consequence of a decrease off drug and an increase on drug. It is not an IFN-specific level increase, and it is not reported as one.

### Co-regulation, held out

Genes were split by Spearman correlation with CLDN4 or TACSTD2 on one arm and tested on the other. Anti-correlated genes (ρ ≤ −0.4) are down on both tests. Genes that move with CLDN4 on the untreated arm (ρ ≥ 0.4, 130 genes after the expression floor) are up on the cisplatin knockdown (median +0.37, 96/130 up). That subset is most of the IFN/APM universe, because most of those genes fall together with CLDN4 off drug. It restates the cisplatin-arm increase. It does not create an untreated increase.

### 中文（这次搜索）

扰动是 SKOV3 的 **TACSTD2/TROP2 shRNA**，不是 CLDN4 敲低。未加药的 104 个 IFN/APM 检验中位数全部为负（0/104 上升）。留一法、连最有利于“上升”的双样本剔除、以及用顺铂臂定义的 CLDN4/TACSTD2 共调控子集，都不能把未加药方向改成上升。

上升只出现在两边都加顺铂时：Hallmark IFN-α 中位 log2FC **+0.41**（79 个基因，60 升 / 19 降，搜索内 BH q = 2.3×10⁻⁶）。自定义 IFN 集 +0.50。去掉 HLA-A 的 APM +0.40（q = 0.017）。C4 六基因里 IFIT1 +2.29、ISG15 +0.50，但 HLA-A 仍为 −2.11，六基因样本检验 p = 0.071。同一对比里 E2F、G2/M、上皮身份基因不升。交互项很大，是因为未加药降幅被减掉了，EMT 的交互同样为正，不把它当成 IFN 水平上升。

---

# Earlier slice — is the untreated IFN/APM decrease specific?

**Slice:** `scripts/gse245459_specificity/` and `results/gse245459_specificity/`
**Dataset:** [GSE245459](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245459) (Han et al., *J Cancer* 2024). Author FPKM, SKOV3, n = 3 per group.
**Prior slice:** PR #97 (`results/w200/C4_GSE245459/`). Untreated shTACSTD2 lowers TACSTD2 and CLDN4 and **closes** IFN/APM. That direction is kept here. The four published no-drug medians (C4 panel −0.571, APM −0.402, IFN −0.445, junction −0.229) reproduce with absolute difference 0.

This slice asks a different question. The identity genes that fall with CLDN4 (claudins, EPCAM, KRT8/KRT19, CD24, ELF3) can be read as an epithelial collapse, and a collapse can drag IFN/APM with it or make a relative IFN rise look real. The quantification below separates those readings.

FPKM columns are relative abundance. Column sums sit in a narrow band (2.87–3.00×10⁵) in every library, so this matrix cannot show an absolute transcriptional shutdown. A global crash in these units is a broad relative shift: most expressed genes down, housekeeping down, and the lost mass absorbed by a few survivors.

---

## English

### Test

Primary contrast: **shTACSTD2 vs shNC, no cisplatin** (the C4 analog). Secondary: the same knockdown on cisplatin (shTACSTD2+DDP vs shNC+DDP), which PR #97 already reported as a drug × knockdown interaction.

log2FC is the difference of mean `log2(FPKM+1)`. Genes enter the specificity tests when mean FPKM is at least 0.3 in the treatment group or the control group, so a gene that is present in scramble and then falls is kept, and genes that are absent in both groups are not. Background for that contrast is every other gene passing the same floor (no-drug n = 17,525; cisplatin n = 18,032).

For each program: median log2FC, Mann–Whitney against that background, an unmatched permutation of the median (5,000 draws), and an expression-decile matched permutation (10 deciles of mean log2(FPKM+1) in the six libraries of the contrast). The matched test is the one used for the call. Unmatched tests on housekeeping and ribosome look “up” because, in the no-drug contrast, log2FC rises with expression (Spearman ρ = +0.25). Matching removes that.

The all-gene median, with no expression floor, is what PR #97 published. It is reported as `median_log2FC_all` and is not the specificity call.

Programs, fixed before this comparison: the prior C4 panel, APM, and IFN/ISG lists; junction; an epithelial-identity list that **excludes TACSTD2 and CLDN4** so the knockdown gene and the C4 gene are not the evidence of collapse; keratins; housekeeping; ribosome; cell cycle; mesenchymal; apoptosis.

### No-drug knockdown — IFN/APM goes down, and the genome does not crash

| Program | n | median log2FC | up/down | vs genome | expression-matched p |
|---|---:|---:|---|---:|---:|
| Genome (eligible) | 17,525 | **−0.152** | 63% down; 6.0% below −1 | — | — |
| Housekeeping | 15 | +0.255 | 11/4 | +0.41 | **0.42** |
| Ribosome | 23 | +0.339 | 21/2 | +0.49 | **0.15** |
| Cell cycle | 20 | **+1.505** | **20/0** | +1.66 | **0.0002** |
| Junction | 19 | −0.415 | 3/16 | −0.26 | 0.055 |
| Epithelial identity, no TACSTD2/CLDN4 | 18 | **−1.282** | 2/16 | −1.13 | **0.0002** |
| Keratin | 4 | −2.602 | 1/3 | −2.45 | **0.0002** |
| IFN/ISG | 44 | **−0.531** | 2/42 | −0.38 | **0.0002** |
| APM | 18 | **−0.448** | 3/15 | −0.30 | **0.0072** |
| APM without HLA-A | 17 | −0.402 | 3/14 | −0.25 | **0.020** |
| C4 panel (6 genes) | 6 | −0.571 | 0/6 | −0.42 | 0.059 |

Genome median −0.15 is a mild drift. Six percent of eligible genes fall more than twofold. Housekeeping and ribosome are null once expression is matched. The cell-cycle set moves the other way: 20/20 up, median +1.50, and every gene in that set has q < 0.05 in the prior genome-wide FDR. shTACSTD2 libraries still sum to the same FPKM total. They do detect fewer genes at FPKM ≥ 1 (11,326–11,639 vs 12,569–12,871 in shNC), about a tenth fewer, which is a real drop in complexity and is already why the matched test exists.

IFN/ISG (median −0.53, 42 of 44 down, matched p = 0.0002) and APM (median −0.45, matched p = 0.007) sit below that genome. Dropping HLA-A leaves APM down (median −0.40, matched p = 0.020). The six C4 genes are all down (median −0.57, identical to PR #97). With only six genes the matched permutation is p = 0.059; without HLA-A the median is −0.52 (5/5 down) and the matched p is 0.16. The powered call is the 44-gene IFN set and the 18-gene APM set, not the six-gene panel alone.

Direct comparisons, same eligible genes:

| Comparison | median difference (IFN or APM minus reference) | MW p |
|---|---:|---:|
| IFN vs housekeeping | −0.79 | 6.5×10⁻⁸ |
| IFN vs ribosome | −0.87 | 6.6×10⁻¹¹ |
| IFN vs cell cycle | −2.04 | 1.9×10⁻¹⁰ |
| IFN vs junction | −0.12 | 0.47 |
| IFN vs epithelial identity | **+0.75** | 0.073 |
| APM vs epithelial identity | **+0.83** | 0.033 |
| APM without HLA-A vs epithelial identity | +0.88 | 0.014 |

IFN/APM is below housekeeping and opposite the cell-cycle rise. It is not below the epithelial-identity set. The identity genes fall further.

Those identity genes, no-drug log2FC, excluding TACSTD2 and CLDN4: CD24 −4.90, KRT19 −4.39, CLDN1 −3.86, KRT80 −2.82, KRT8 −2.39, SPINT2 −1.52, CLDN3 −1.48, CLDN7 −1.39, F11R −1.31, EPCAM −1.25, SPINT1 −1.19, ELF3 −1.02 (all q < 0.05), then MUC1 −0.42 (q = 0.035). Spared: TJP1 −0.13 (q = 0.39), OCLN +0.10 (q = 0.35), and the abundant keratin KRT18 +0.25 (q = 0.14, contrast-mean FPKM 20). The keratin median of −2.60 is KRT8/KRT19/KRT80. It is not a KRT18 collapse, and it is not a tight-junction scaffold collapse. CLDN4 itself remains the PR #97 result (log2FC −1.92, q = 0.0037) and was left out of the reference set on purpose.

Fold-change and FPKM mass are different. IFN genes account for 0.93% of the negative FPKM mass among eligible genes; the epithelial-identity set accounts for 0.22%. Cell cycle contributes none of the negative mass and gains +969 FPKM. The largest gainers are ACTB, mitochondrial transcripts, HSP90AB1, HSPA8, HSPA5, MMP10. The top 25 gainers are 19% of positive mass, so the compensation is spread out. The identity and IFN shifts are large in log2FC and small in transcriptome mass.

**No-drug call.** IFN/APM decreases, and that decrease is specific relative to the genome, housekeeping, and ribosome. The same libraries show a partial epithelial-identity loss that is larger than the IFN/APM decrease, with TJP1, OCLN, and KRT18 spared, plus a coherent cell-cycle increase. There is no IFN/APM increase in this arm.

### Cisplatin background — the IFN/APM increase is specific

This is the only arm in which IFN/APM goes up. It is the interaction from PR #97, reported again only to answer the specificity question.

| Program | n | median log2FC | up/down | expression-matched p |
|---|---:|---:|---|---:|
| Genome | 18,032 | **+0.042** | — | — |
| Housekeeping | 16 | −0.136 | 4/12 | **0.93** |
| Junction | 16 | +0.074 | 8/8 | **0.85** |
| Epithelial identity, no TACSTD2/CLDN4 | 16 | −0.160 | 4/12 | 0.084 |
| Keratin | 4 | **−1.333** | 0/4 | **0.0002** |
| Cell cycle | 20 | −0.354 | 3/17 | 0.12 |
| IFN/ISG | 37 | **+0.504** | 30/7 | **0.0002** |
| APM | 18 | **+0.397** | 14/4 | **0.0020** |
| APM without HLA-A | 17 | +0.397 | 14/3 | **0.0032** |
| C4 panel | 4 eligible (6 in the all-gene median) | +0.410 (all-gene **+0.278**) | 3/1 | 0.15 |

Genome median is +0.04. IFN (matched p = 0.0002) and APM (matched p = 0.002, unchanged without HLA-A) sit above it. Keratin keeps falling (median −1.33, 4/4 down). IFN vs keratin MW p = 2.0×10⁻⁵; IFN vs epithelial identity p = 2.6×10⁻⁴; IFN vs housekeeping p = 2.6×10⁻⁴. The six-gene C4 panel is underpowered here (matched p = 0.15); its all-gene median (+0.278) matches PR #97, where the panel-level MW was null.

**Cisplatin call.** The IFN/APM increase on this background is specific relative to a flat genome and relative to keratin, which is still down. It is not an epithelial rebound and it is not a global shift. CLDN4 is not further reduced in this arm (PR #97: log2FC −0.14, p = 0.099). This interaction is not the untreated C4 test.

### Bottom line

| Question | Call |
|---|---|
| Does untreated shTACSTD2 raise IFN/APM? | No. IFN median −0.53, APM −0.45. Same direction as PR #97. |
| Is that decrease a global crash? | No. Genome median −0.15, housekeeping and ribosome null after expression matching, cell cycle +1.50 (20/20). |
| Is there an epithelial-identity loss? | Yes, partial. Claudins, EPCAM, ELF3, CD24, KRT8/KRT19/KRT80 fall. TJP1, OCLN, and KRT18 do not. It is larger than the IFN/APM decrease. |
| Where IFN/APM does rise (cisplatin background), is the rise specific? | Yes. Genome median +0.04, IFN +0.50, APM +0.40, keratin −1.33. |

### Caveats

Ovarian SKOV3, not lung. The perturbation is TACSTD2 shRNA, not CLDN4 KD. Author FPKM, n = 3, Welch tests, no limma/voom. The keratin set has four genes above the floor, and KRT18 disagrees with the other three. HLA-A (log2FC −4.17) is extreme; APM without it is still down in the no-drug arm and still up on cisplatin. No protein, surface MHC, or killing assay in this deposit.

### Reproduce

```
python3 scripts/gse245459_specificity/analyze.py
python3 scripts/gse245459_specificity/sweep_thesis_up.py
```

The script downloads `GSE245459_fpkm.anno.txt.gz` from GEO if `results/gse245459_specificity/raw/` is empty. Tables are in `results/gse245459_specificity/tables/`. Figure: `results/gse245459_specificity/figures/ifn_vs_global_epithelial.png`.

---

## 中文

主对比是未加顺铂的 shTACSTD2 vs shNC。PR #97 的四个中位数原样复现（C4 −0.571，APM −0.402，IFN −0.445，连接 −0.229）。

未加药：全基因组中位 log2FC 为 −0.15（17,525 个基因里只有 6% 低于 −1）。看家基因和核糖体在表达量匹配后为 null。细胞周期 20/20 上调，中位 +1.50。IFN/ISG 中位 −0.53（44 个基因，42 个下调，匹配 p = 0.0002），APM 中位 −0.45（去掉 HLA-A 后 −0.40，匹配 p = 0.020）。上皮身份基因（不含 TACSTD2/CLDN4）中位 −1.28，降幅大于 IFN/APM（APM 对上皮 MW p = 0.033）。KRT8/19/80、CLDN1/3/7、EPCAM、ELF3、CD24 明显下降；TJP1、OCLN、KRT18 不降。所以未加药时 IFN/APM 是下降，而且相对全基因组和看家基因是特异的；同时存在一部分上皮身份丢失，幅度更大。这不是全转录组崩溃。

顺铂背景是 IFN/APM 上升的唯一臂：基因组中位 +0.04，IFN +0.50（匹配 p = 0.0002），APM +0.40（去掉 HLA-A 仍 p = 0.003），角蛋白仍为 −1.33。上升是特异的，不是角蛋白回升，也不是全基因组平移。这是药×敲低交互，不是未加药的 C4 检验。
