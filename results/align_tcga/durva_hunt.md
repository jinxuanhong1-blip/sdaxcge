# Durvalumab NSCLC RNA hunt

Question: which **public** durvalumab NSCLC bulk-RNA matrix produced the
user-cited Spearman values ρ ≈ −0.65 / −0.46 (TACSTD2 vs immune / cytotoxic)?

## What was searched

Open catalogues only (GEO, cBioPortal public studies, Xena, papers that
deposit a processed matrix). Controlled-access / pharma-only files were
not used (PACIFIC, MYSTIC, NEPTUNE, COAST, NeoCOAST raw RNA, OAK/POPLAR
Genentech transcriptome).

| Candidate | Why considered | Public processed matrix? | Used? |
|-----------|----------------|--------------------------|-------|
| GSE253564 | Altorki *Nat Commun* 2024 neoadjuvant durvalumab ± SBRT, **pre-treatment** FPKM | Yes (`…_Pre-treatment_Samples_Pubs_FPKMs.txt.gz`, n=32) | Yes |
| GSE248378 | Same trial, **post-treatment** FPKM | Yes (`GSE248378_Durva_Post_FPKMs.txt.gz`, n=29) | Yes |
| GSE131933 | Durvalumab ± tremelimumab scRNA (not bulk NSCLC) | scRNA, not a bulk TACSTD2-vs-signature matrix | No |
| Bessede *CCR* / OAK+POPLAR | TACSTD2 vs ICI outcome, n=891 | Atezolizumab (not durvalumab); transcriptome not open | No |
| PACIFIC / MYSTIC / COAST / NeoCOAST | Landmark durvalumab NSCLC trials | No open gene-level matrix found | No |

No other GEO series with (i) durvalumab-treated NSCLC, (ii) a downloadable
gene × sample expression matrix, and (iii) enough samples for Spearman
was identified as of this run.

## What the open matrices actually give

Computed with the same signatures and Spearman / ESTIMATE-partial method
as the TCGA analysis (`durva_tacstd2_correlations.csv`). ESTIMATE scores
come from a Python port of Yoshihara 2013 that matches the official
`sample_estimate.gct` (max \|Δ\| < 1e-9).

**GSE253564 pre-treatment (n=32)** — closest match to the cited pair:

| Feature | Spearman ρ | p | partial ρ (ESTIMATE) | partial p |
|---------|-----------:|---:|---------------------:|----------:|
| Exhaustion | −0.642 | 7.4e-5 | −0.357 | 0.049 |
| Cytolytic_CYT | −0.461 | 0.0079 | −0.081 | 0.66 |
| Tcell_inflamed_GEP | −0.470 | 0.0067 | +0.021 | 0.91 |
| CD8_Tcell | −0.552 | 0.0011 | −0.173 | 0.35 |

The pair **Exhaustion −0.64 / CYT −0.46** (or Exhaustion / GEP −0.47) is
the same direction and the same magnitude as the user-cited −0.65 / −0.46.
We cannot prove this is their exact signature definition.

**GSE248378 post-treatment (n=29)** — same direction, larger |ρ|
(CYT −0.81, CD8 −0.76, GEP −0.63, IFNγ-6 −0.51). IFNG itself is absent
from this matrix; the 6-gene IFNγ score uses the other 5 genes. Most
signatures remain negative and significant after ESTIMATE (IFNγ-6 does not).

Neither GEO matrix ships a published purity column; the partial uses
`purity_proxy = −ESTIMATEScore`.
