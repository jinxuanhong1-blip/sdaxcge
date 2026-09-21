# DepMap lung lines: CLDN4 vs DNA-damage signatures, STING, and HLA

Additive public analysis. DepMap Public 24Q4 RNA (`log2(TPM+1)`), lung cell lines with expression: **n = 214** (260 lung models in `Model.csv`; 46 have no 24Q4 RNA). Same lineage filter as PR #304: `OncotreeLineage == Lung` and `ModelType == Cell Line`. NSCLC **n = 143** (LUAD 80 / LUSC 27 / other NSCLC 36). SCLC+NET n = 60; other lung n = 11. Cultured lines: basal transcription, no infiltrate, no IFN treatment.

MHC-I (HLA-A/B/C + B2M) ρ = **+0.217** (p = 0.0014) and Hallmark IFN-γ ρ = **+0.284** (p = 2.5×10⁻⁵) reproduce PR #304 on this n. They are not re-claimed.

## Paper sentence

CLDN4 RNA is positively associated with a STING-core score (ρ = **+0.277**, p = 4.1×10⁻⁵, q = 1.6×10⁻⁴) and with classical HLA-A/B/C (ρ = **+0.228**, p = 7.8×10⁻⁴, q = 0.0010). Both remain after Oncotree-group adjustment (partial ρ = +0.243 and +0.192). A 19-gene DSB-response RNA score is inverse (ρ = **−0.260**, q = 2.4×10⁻⁴; partial ρ = −0.161, p = 0.018). That inverse is absent in LUAD (n = 80, ρ = +0.054, p = 0.64). **γH2AX (H2AX-pS139) was not measured.** Valid MCLP phospho-DDR antibodies ATM pS1981 and Rad17 pS645 are null versus CLDN4 (n = 159; ρ = −0.101 and +0.027).

## γH2AX is absent

| source | samples | antibodies | H2AX-pS139 |
|---|---:|---:|---|
| CCLE RPPA 20180123 | 899 | 214 | absent |
| CCLE RPPA 20181003 | 899 | 214 | absent |
| MCLP CCLE RPPA 20221116 | 878 | 447 | absent |

MCLP histones on that release are total H2A, H2B, H3, H3-pS10, and a few H3 methylation marks. None is γH2AX. H2AX **transcript** is a different measurement and runs the other way (ρ = +0.172, p = 0.012). It is not in the DSB score.

## RNA (n = 214)

Signatures are the mean of gene-wise z-scores computed inside the cohort being tested. Spearman, two-sided. BH q is within the four primary scores, and separately within HLA-A/B/C. 95% CIs are 5,000 bootstrap resamples (seed 0).

| endpoint | ρ | p | q | 95% CI |
|---|---:|---:|---:|---|
| Hallmark DNA repair (150/150) | −0.160 | 0.019 | 0.019 | [−0.29, −0.04] |
| DSB-response RNA (19 genes) | −0.260 | 1.2×10⁻⁴ | 2.4×10⁻⁴ | [−0.38, −0.14] |
| STING core (CGAS, STING1, TBK1, IKBKE, IRF3) | +0.277 | 4.1×10⁻⁵ | 1.6×10⁻⁴ | [+0.15, +0.40] |
| HLA-A/B/C | +0.228 | 7.8×10⁻⁴ | 0.0010 | [+0.09, +0.35] |
| HLA-A | +0.186 | 0.0065 | 0.0086 | [+0.05, +0.32] |
| HLA-B | +0.179 | 0.0086 | 0.0086 | [+0.04, +0.31] |
| HLA-C | +0.280 | 3.2×10⁻⁵ | 9.5×10⁻⁵ | [+0.15, +0.40] |

NSCLC (n = 143): DSB ρ = −0.256 (q = 0.0080); STING core ρ = +0.197 (q = 0.037); HLA-A/B/C ρ = +0.182 (q = 0.040). Hallmark DNA repair ρ = −0.159 (q = 0.057; 95% CI [−0.31, −0.01]). Inside NSCLC, HLA-C holds (ρ = +0.247, q = 0.0089); HLA-A (ρ = +0.119, p = 0.16) and HLA-B (ρ = +0.154, p = 0.067) do not.

Partial Spearman (ranks residualized on Oncotree group):

| endpoint | all lung, 5 groups | NSCLC, 3 groups |
|---|---|---|
| Hallmark DNA repair | −0.064 (p = 0.35) | −0.141 (p = 0.094) |
| DSB-response RNA | −0.161 (p = 0.018) | −0.254 (p = 0.0022) |
| STING core | +0.243 (p = 3.3×10⁻⁴) | +0.183 (p = 0.028) |
| HLA-A/B/C | +0.192 (p = 0.0049) | +0.170 (p = 0.043) |

The DNA-repair inverse association does not survive subtype adjustment. The DSB inverse association does, but it is not a LUAD result: LUAD ρ = +0.054 (p = 0.64, n = 80). It is present in LUSC (ρ = −0.570, p = 0.0019, n = 27) and other NSCLC (ρ = −0.572, p = 2.6×10⁻⁴, n = 36). SCLC+NET DSB is null (ρ = +0.034, p = 0.80, n = 60).

STING-core is not CGAS. Gene-level, all lung: STING1 ρ = +0.309 (p = 4.1×10⁻⁶), IKBKE ρ = +0.260 (p = 1.2×10⁻⁴), IRF3 ρ = +0.251 (p = 2.1×10⁻⁴), TBK1 ρ = −0.101 (p = 0.14), CGAS ρ = +0.050 (p = 0.46). STING core correlates with Hallmark IFN-γ at ρ = 0.719, so this is the same direction as PR #304, not an independent IFN story. LUAD STING1 remains positive (ρ = +0.372, p = 6.7×10⁻⁴, n = 80); the five-gene STING-core score in LUAD does not (ρ = +0.201, p = 0.074).

## Protein (MCLP CCLE RPPA 20221116, norm level 4)

Joined to DepMap lung RNA on CCLE name: **n = 159** (NSCLC 111). Five lung models have no CCLE name and are not joined. There were no duplicate CCLE names. RNA and RPPA are not the same lysate. BH q is within the five Valid primary antibodies.

| antibody | validation | ρ | p | q | 95% CI |
|---|---|---:|---:|---:|---|
| ATM pS1981 | Valid | −0.101 | 0.21 | 0.52 | [−0.26, +0.06] |
| Rad17 pS645 | Valid | +0.027 | 0.74 | 0.82 | [−0.14, +0.19] |
| cGAS | Valid | +0.052 | 0.52 | 0.82 | [−0.10, +0.21] |
| IRF3 | Valid | +0.247 | 0.0017 | 0.0086 | [+0.09, +0.39] |
| HLA-DQA1 | Valid | +0.019 | 0.82 | 0.82 | [−0.14, +0.18] |

IRF3 protein agrees with the RNA direction. cGAS protein does not, matching CGAS RNA. Classical HLA-A/B/C protein is not on this array; HLA-DQA1 is class II and is null. In NSCLC (n = 111) none of the five Valid primaries has q < 0.05 (IRF3 ρ = +0.211, q = 0.10).

Reported, not primary: Chk1 pS345 and Chk2 pT68 are MCLP **Caution** antibodies (ρ = −0.222, p = 0.0049, and ρ = −0.091, p = 0.25). Total ATM protein is inverse (ρ = −0.330, p = 2.1×10⁻⁵) and is not the phospho mark. Histone H3 pS10, a mitotic mark, is also inverse (ρ = −0.242, p = 0.0021), so an inverse RNA DDR score is not specific evidence of fewer double-strand breaks.

## What this is

Cell-intrinsic basal RNA and RPPA. Positive STING and HLA means CLDN4-high lung lines are not STING-low or HLA-low in this matrix. That is the same direction as the IFN-γ / MHC-I result in PR #304. It is not a statement about immune cells next to a tumor cell.

## Reproduce

```bash
python3 methods/ccle_cldn4_ddr_sting_hla/download.py data/ccle_cldn4_ddr_sting_hla
python3 methods/ccle_cldn4_ddr_sting_hla/analyze.py \
  --data data/ccle_cldn4_ddr_sting_hla \
  --outdir methods/ccle_cldn4_ddr_sting_hla
```

Sources: DepMap 24Q4, doi:10.25452/figshare.plus.27993248.v1; MSigDB 2024.1 Hs Hallmark; MCLP/TCPA CCLE RPPA release 20221116 (https://tcpa.drbioright.org/rppa500mclp/); CCLE RPPA headers from https://data.broadinstitute.org/ccle/.
