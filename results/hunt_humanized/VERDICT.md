# Verdict

**No public CD34- or PBMC-humanized NSCLC IO RNA/scRNA dataset measures CD8 nest entry, and none can test TACSTD2/CLDN4 versus nest entry.**

Nest entry is a spatial phenotype (CD8 inside epithelial tumor nests vs stroma). Dissociated bulk and scRNA destroy that axis. No public HIS-NSCLC spatial transcriptomic or nest-scored multiplex dataset was found.

What exists and was actually opened:

| Dataset | What it is | TACSTD2 / CLDN4 vs CD8 nest? |
|---|---|---|
| **GSE260575** | PBMC → MHC I/II dKO NSG + H460 + pembro ± BJIKT, QuantSeq 3', n=12 | **No.** TACSTD2 = 0/12. CD8A/CD8B at the floor. CLDN4 present; vs CD8A ρ ≈ −0.22, p ≈ 0.48. |
| **GSE276724** | Paper: CD34-NSG + NSCLC PDX + denfivontinib/pembro. Deposit: 5 PDX × 3 unlabeled bulk counts | **No.** Cross-PDX bulk CD8 ≠ nest. TACSTD2 vs CD8A ρ = −0.26 (n=15) → **−0.09 after dropping one epithelial-dropout sample**. Effective n is 5 models, not 15 mice. |
| **GSE293914** | PBMC + H1975 scRNA, 449 MB, 1 pooled GEO sample | **Not nest.** No cell metadata. Paper IF is not a TACSTD2/CLDN4–CD8 nest score. |
| Yu 2023 HuNCG H226 10x | Real CD34 HIS LUSC scRNA on paper | **Not public** (“available upon request”). |
| SITC 2024 #807 | CD34 SGM3 + H358 + TROP2 ADC + sabestomig, bulk RNA + CellDive | **Not deposited.** Closest TROP2+HIS IO study; closed. |
| Meraz eLife NPRL2 | CD34 NSG + A549 mets, NanoString immune panel n=12 | **No GEO.** Panel is immune-centric; not a TROP2/CLDN4 nest assay. |

Do not treat patient PBMC scRNA (e.g. GSE285888) or hHGF-NSG “humanized” mice (GSE296572) as HIS tumor models.

Honest residual: GSE276724 shows a weak, purity-tangled, PDX-clustered anti-correlation of epithelial TACSTD2/CLDN4 with bulk CD8B. That is not nest entry and should not be cited as one.
