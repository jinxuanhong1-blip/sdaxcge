# Claims Index / 主张索引

TROP2 (TACSTD2) × immune-cold × tight-junction (CLDN4) hypothesis tracker.
TROP2（TACSTD2）× 免疫冷 × 紧密连接（CLDN4）假说的主张追踪表。

- **Source / 来源:** User PPT 2026-08-17, 金炫宏
- **Rule / 规则:** Do **not** invent results. User-reported statistics are labeled as such; every page
  records whether the public data/analog is **found / not found** and whether independent
  verification was run (none has been run in this repo yet).
  **不得编造结果。** 用户报告的统计量已明确标注；每页记录公共数据/替代数据是否 **已找到 / 未找到**，
  以及是否已独立验证（本仓库尚未运行任何分析）。

## Page format / 页面格式
Each page contains / 每页包含:
1. **Claim / 主张** — as stated by the user / 用户原始陈述
2. **User number / 用户编号** — e.g., A1, B3, C7
3. **Public test plan / 公共验证方案**
4. **Database query / 数据库查询** — GEO / Xena / TISMO / CCLE(DepMap) / cBioPortal / PubMed
5. **Done = found / not found / 完成状态 = 已找到 / 未找到**

## Legend / 图例
- **FOUND / 已找到** — a public dataset or analog is identified and accessible.
- **NOT FOUND (to hunt) / 未找到（待检索）** — public analog still needs locating/confirming.
- **PRIVATE / 私有** — internal data; not publicly reproducible (list public analog to hunt).
- **NOT DONE / 未执行** — analysis not yet run in this repo.

---

## A — TROP2 immune-cold association / TROP2 免疫冷关联 (A1–A11)

| # | Claim (short) / 主张（简） | Key dataset / 关键数据 | Public data / 公共数据 |
|---|---|---|---|
| [A1](A1.md) | TACSTD2 vs immune negative after purity (TCGA+OncoSG) | TCGA-LUAD, OncoSG | FOUND |
| [A2](A2.md) | Durvalumab ρ=-0.65; purity-adj ρ=-0.46 (p=2e-4) | durvalumab cohort | NOT FOUND (accession unconfirmed) |
| [A3](A3.md) | scRNA malignant TACSTD2 NMPR>MPR; vs T/NK ρ -0.40~-0.50 | GSE207422 | FOUND |
| [A4](A4.md) | Tacstd2 up 49/64 models (p=5.8e-5) | TISMO | FOUND |
| [A5](A5.md) | TROP2-high = immune-resistant subset | GSE76628 | FOUND |
| [A6](A6.md) | TROP2-high tumor frac up, CD8/NK down | TCGA/OncoSG/GSE207422 | FOUND |
| [A7](A7.md) | Zhejiang 25-pair IHC (PRIVATE); hunt paired pre/post | PRIVATE | PRIVATE / hunt |
| [A8](A8.md) | Keratin/TJ up, EMT down (GSEA) | TCGA-LUAD + MSigDB | FOUND |
| [A9](A9.md) | Intersection CLDN1/4/7, F11R, PARD3 | TCGA/OncoSG | FOUND |
| [A10](A10.md) | TFs ELF3/GRHL1/KLF4/TFAP2A up; NKX2-1 down | TCGA/OncoSG | FOUND |
| [A11](A11.md) | Galectin/Nectin/TGF-β/CD47 axes | TCGA/OncoSG + MSigDB | FOUND |

## B — CLDN4 tight-junction barrier / CLDN4 紧密连接屏障 (B1–B7)

| # | Claim (short) / 主张（简） | Key dataset / 关键数据 | Public data / 公共数据 |
|---|---|---|---|
| [B1](B1.md) | CLDN4 top surface coexpression with TROP2 | TCGA pan-cancer | FOUND |
| [B2](B2.md) | CLDN4–TROP2 protein ρ=0.69 | CCLE/DepMap proteomics | FOUND |
| [B3](B3.md) | TJ-high → low CD8/GEP (p<1e-6) | TCGA-LUAD | FOUND |
| [B4](B4.md) | NR higher TJ (p=0.019) | GSE126044 | FOUND |
| [B5](B5.md) | 11-cohort ICI meta CLDN4 OR=0.42 | 11 ICI cohorts | PARTIAL / confirm |
| [B6](B6.md) | Spatial: CLDN4 avoids immune niches | spatial NSCLC | NOT FOUND (to hunt) |
| [B7](B7.md) | CLDN4 TEER / barrier literature | PubMed literature | to curate (real PMIDs) |

## C — Private validation / 私有验证 (C1–C10)
All C pages are **PRIVATE** (SKB264 / PDX / mIF / Co-IP / AF2) and list the public analog to hunt.
Per-claim themes below are **inferred** and should be confirmed against the PPT.
所有 C 页面均为 **私有**（SKB264 / PDX / mIF / Co-IP / AF2），并列出待检索的公共替代数据。
以下每条主题为 **推断**，需与 PPT 核对确认。

| # | Inferred theme / 推断主题 | Type / 类型 | Status / 状态 |
|---|---|---|---|
| [C1](C1.md) | SKB264 efficacy in TROP2-high | SKB264 ADC | PRIVATE / hunt |
| [C2](C2.md) | SKB264 + ICI combination | SKB264 + ICI | PRIVATE / hunt |
| [C3](C3.md) | PDX TROP2-high response | PDX | PRIVATE / hunt |
| [C4](C4.md) | Humanized/immune PDX | PDX (immune) | PRIVATE / partial (TISMO) |
| [C5](C5.md) | mIF TROP2/CLDN4 vs CD8 exclusion | mIF | PRIVATE / hunt |
| [C6](C6.md) | mIF pre/post dynamics | mIF | PRIVATE / hunt |
| [C7](C7.md) | Co-IP TROP2–CLDN4 | Co-IP | PRIVATE / hunt |
| [C8](C8.md) | Co-IP TROP2–TJ complex | Co-IP | PRIVATE / hunt |
| [C9](C9.md) | AF2 TROP2–CLDN4 interface | AF2 | PRIVATE / reproducible via ColabFold |
| [C10](C10.md) | AF2 TROP2–TJ complex | AF2 | PRIVATE / reproducible via ColabFold |

---

## Notes / 说明
- Statistics such as ρ, OR, and p-values are **user-reported** and are recorded verbatim, not
  independently confirmed here. 统计量（ρ、OR、p 值等）均为 **用户报告**，仅照录，未在此独立确认。
- Confirm exact accessions for A2 (durvalumab) and the B5 11-cohort list from the source PPT.
  A2（durvalumab）与 B5（11 队列）的确切数据编号需从 PPT 核对确认。
