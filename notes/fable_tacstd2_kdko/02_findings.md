# Findings — TROP2/TACSTD2 loss and CLDN4 / junctions / IFN-immune programs

All numbers come from `results/fable_tacstd2_kdko/` (per-dataset `*_pergene_panel.tsv`,
`*_setlevel.tsv`, and the `SYNTHESIS_*` tables). Contrast is always
**perturbation (KD/KO) vs control**, so `log2FC > 0` = up upon TROP2/TACSTD2 loss.

## 1. The knockdown/knockout worked (perturbation QC)

| Dataset | Contrast | TACSTD2 log2FC | p | padj/FDR |
|---|---|---|---|---|
| GSE334497 | mouse 4T1 Trop2 KO | **-3.94** | 1.8e-3 | 0.71† |
| GSE289287 | human T-47D Trop-2 KO xenograft | **-3.26** | 1.6e-65 | 1.9e-61 |
| GSE245459 | human ovarian shTACSTD2 | **-2.63** | 5.3e-5 | 4.1e-3 |
| GSE15212 | human colorectal siTACSTD2 | **-3.02** | 5.3e-15 | 1.0e-10 |
| GSE289287_DSG2tumor | DSG2 KO (partner control) | -0.15 | 0.28 | 0.54 |
| GSE289287_DSG2cell | DSG2 KO (partner control) | +0.21 | 0.15 | 0.35 |

† FDR high only because a plain per-gene Welch t-test on n=5/5 dilutes single-gene FDR;
the effect size (−3.9 log2, ~15-fold) and nominal p are unambiguous. The two DSG2-KO
sets correctly leave TACSTD2 unchanged (specificity control).

## 2. CLDN4 (headline gene)

| Dataset | log2FC | significance |
|---|---|---|
| GSE334497 (mouse 4T1 KO) | -0.82 | ns (p=0.25), **down trend** |
| GSE289287 (human T-47D KO) | +0.28 | ns (padj=0.47) |
| GSE245459 (human ovarian KD) | **-1.92** | **padj=3.7e-3 (down)** |
| GSE15212 (human CRC KD) | +0.27 | ns |
| GSE289287_DSG2tumor (partner) | **-0.91** | **padj=2.1e-4 (down)** |

CLDN4 is **not uniformly regulated** by TROP2 loss: significantly **down** only in the
ovarian knockdown, a down trend in mouse 4T1, and slight (ns) up in the two human
breast/CRC contrasts. Loss of the desmosomal partner DSG2 significantly lowers CLDN4.
Net: when CLDN4 does move significantly, it goes **down** — consistent with TROP2
supporting (rather than repressing) the claudin program — but the effect is
context-dependent, not universal.

## 3. Junctions (tight junction + adherens + desmosome)

Competitive gene-set test direction (`***`p<1e-3, `**`p<1e-2, `*`p<0.05):

| Contrast | junctions | tight_junction | desmosome |
|---|---|---|---|
| mouse 4T1 Trop2 KO | **down***** | **down***** | **down***** |
| human T-47D Trop-2 KO | up ns | up ns | up ns |
| human ovarian shTACSTD2 | **down***** | **down***** | **down***** |
| human CRC siTACSTD2 | down ns | down ns | down (sc p=0.02) |
| DSG2 KO tumour (partner) | down ns | down * | up ns |

TROP2 loss **weakens the junction/claudin/desmosome program**, strongly in mouse breast
and human ovarian (both p<1e-7), directionally down in colorectal. The human T-47D
Trop-2 KO xenograft is the outlier (slight, non-significant up). Desmosome genes are the
most consistently reduced arm (down in 4T1 p=1.6e-7, ovarian p=6.7e-5, CRC self-contained
p=0.02).

## 4. IFN / immune genes

| Contrast | ifn_immune | ifn_isg | antigen_presentation | tcell_cytotox |
|---|---|---|---|---|
| mouse 4T1 Trop2 KO | **up***** | up ** | **up***** | **up***** |
| human T-47D Trop-2 KO (xenograft) | **up***** | **up***** | up ns | down ns |
| human ovarian shTACSTD2 | **down***** | **down***** | **down***** | up ns |
| human CRC siTACSTD2 | **up*** | **up***** | down ns | down ns |
| DSG2 KO tumour (partner) | **up***** | **up***** | **up***** | up ns |
| DSG2 KO cells (partner) | **up***** | **up***** | **up*** | down ns |

The most reproducible signal in the whole analysis: **TROP2 loss de-represses
interferon-stimulated genes**. ISGs go **up** in mouse 4T1 (p=5e-12), human T-47D
Trop-2 KO (p=1.9e-8), and human colorectal (p=6e-3); loss of the desmosomal partner DSG2
reproduces a very large ISG induction (up to p=1e-27). Because the T-47D xenografts grow
in immunodeficient hosts, this ISG induction is **tumour-cell-intrinsic**, not merely
immune-cell infiltrate. In the immunocompetent 4T1 model the immune response is broader:
antigen presentation and T-cell/cytotoxicity signatures also rise — the transcriptomic
correlate of the "immune exclusion reversed by TROP2 loss" phenotype reported for that
series. The human **ovarian** knockdown is the sole exception (ISGs and antigen
presentation go down).

## 5. One-line answer

Yes — TROP2/TACSTD2 loss reproducibly changes all three axes, but with
**context-dependent direction**:
- **IFN/immune**: up (de-repressed) in 3/4 loss datasets + both DSG2-KO controls; the
  single exception is ovarian. Cell-intrinsic ISG induction is the strongest, most
  reproducible effect.
- **Junctions**: weakened (down), strongly in mouse breast and human ovarian; desmosome
  arm most consistently down.
- **CLDN4**: not universal — significantly down in the ovarian KD (and with DSG2 loss),
  down trend in mouse 4T1, ns elsewhere.

## Caveats
- GSE245459 and the T-47D contrasts have small n (3/3 and 4/3); GSE334497 gives only
  normalized (not raw) counts, so its per-gene test is Welch/BH rather than DESeq2.
- Competitive gene-set p-values use a simplified CAMERA-style Mann-Whitney that ignores
  inter-gene correlation, so treat them as effect-direction evidence, corroborated by the
  self-contained Wilcoxon and per-gene tables, rather than exact FWER-controlled claims.
- Mouse↔human ortholog matching is by symbol (case-insensitive) + mygene mapping.
