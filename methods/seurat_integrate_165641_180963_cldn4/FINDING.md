# FINDING — integrate GSE165641 + GSE180963 (Cldn4-only, Seurat + Harmony)

**ADDITIVE. Public mouse. Cldn4-only. No dual-high.** Thesis is already correct and is not rewritten. This folder **integrates** two public KL/K GEMM 10x series — [GSE165641](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165641) (Wang / Zhong, *Adv Sci* 2021, [PMID 34369094](https://pubmed.ncbi.nlm.nih.gov/34369094/)) and [GSE180963](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE180963) (Bai / Guo / Zhang / Long / Dong, [CAN-22-1740](https://doi.org/10.1158/0008-5472.can-22-1740)) — into **one** Seurat object. The point is mouse-level n **after merge**, not two separate catalogs. Primary engine is **R + Seurat 5.5.1** (`CreateSeuratObject` per dataset → `IntegrateLayers` / Seurat_IntegrateLayers_HarmonyIntegration). No Python-only primary. Private 8-KL matrices were not opened. FACS-only epithelium (GSE179502/GSE154989), CD45-only (GSE127465), subcutaneous Cldn4-floor (GSE267321), and injury KO (GSE50927) were **not** merged.

**Verdict.** Both series have processed 10x MTX. An integrated object **was built**. Honest n after merge = **4 mice** (3 KL + 1 K). Cldn4 vs T/NK is scored at that unit. Epithelial IFN/MHC/TJ Q4 vs Q1 is **no-go (n_mice<6)**. Do not write n = 20953 cells. Thesis unchanged.

---

## Decision

| Question | Answer |
|---|---|
| GSE165641 matrix | **yes** — Cell Ranger filtered MTX, 2 KL mice |
| GSE180963 matrix | **yes** — author-QC MTX, 1 K + 1 KL |
| CreateSeuratObject per dataset | **yes** — Seurat 5.5.1 |
| Integration | **Seurat_IntegrateLayers_HarmonyIntegration** (Harmony package 2.0.5) |
| Mice after merge (unit) | **4** |
| Marker epithelium (tight) | **569** cells |
| Cldn4 row present | **yes** — **323** cells > 0 |
| Cldn4 vs T/NK at mouse unit | **scored** — n=4 mice, ρ=-1.000, p=0 |
| Epithelial IFN/MHC/TJ Q4 vs Q1 | **no-go — n_mice<6** |
| Unit | **mouse** |
| Dual-high TACSTD2 × Cldn4 | **not defined** |
| Private 8 KL | **not used** |
| ICI / PD-1 | **no** |

---

## Honest n

| item | n | note |
|---|---:|---|
| GEO series with matrix | **2** | GSE165641 + GSE180963 |
| Mice (unit, after merge) | **4** | not cells |
| K mice | **1** | GSE180963 only |
| KL mice | **3** | GSE165641 KL1/KL2 + GSE180963 KL |
| Cells after QC | **20953** | not the unit |
| Genes after merge | **31053** | symbol union |
| Marker epithelial cells | **569** | Epcam+ structural+ Ptprc- |
| T/NK cells | **7802** | marker call |
| Cldn4-positive cells | **323** | count > 0 |
| Dual-high | **0** | not defined |
| Private 8-KL mice | **0** | public GEO only |
| ICI arms | **0** | untreated GEMM |

Do not write n = 20953. Do not treat two GEO series as two mice.

---

## Per-mouse table (the actual unit)

| mouse | dataset | GSM | genotype | strain | n cells | n epi | n T/NK | frac T/NK | Cldn4 all | Cldn4 %pos | Cldn4 epi | Cldn4 epi %pos | IFN epi | MHC epi | TJ epi |
|---|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| GSE165641_KL1 | GSE165641 | GSM5047302 | KL | C57BL/6 | 4199 | 207 |   97 | 0.023 | 0.0756 | 4.74 | 1.0804 | 64.73 | 0.136 | 0.376 | 0.591 |
| GSE165641_KL2 | GSE165641 | GSM5047303 | KL | C57BL/6 | 2494 | 160 |  253 | 0.101 | 0.0414 | 3.93 | 0.5637 | 51.25 | 0.075 | 0.611 | 0.773 |
| GSE180963_K | GSE180963 | GSM5481386 | K | FVB | 6696 |  28 | 3634 | 0.543 | 0.0001 | 0.01 | 0.0270 | 3.57 | 0.115 | 0.822 | 0.277 |
| GSE180963_KL | GSE180963 | GSM5481387 | KL | FVB | 7564 | 174 | 3818 | 0.505 | 0.0031 | 0.33 | 0.0970 | 9.77 | 0.242 | 1.687 | 0.623 |

---

## Cldn4 vs T/NK (mouse unit)

n=4 mice, ρ=-1.000 (perfect rank anti-correlation; `cor.test` prints p=0) for epithelial Cldn4 **mean** vs T/NK fraction. Same for **%pos**. That rank order **is the two-lab split**: both GSE165641 KL mice sit high-Cldn4 / low-T/NK, both GSE180963 mice sit low-Cldn4 / high-T/NK. It is not a within-lab slope and is **not powered**. Strain / dataset is a confounder (C57BL/6 vs FVB). Harmony corrects the embedding, not the mouse-level scores. Do not cite p=0. Do not write a Cldn4–T/NK law from this object.

Cldn4 vs epithelial IFN / MHC / TJ (descriptive Spearman, same n): IFN n=4 mice, ρ=0.000, p=1; MHC n=4 mice, ρ=-0.800, p=0.2; TJ n=4 mice, ρ=0.400, p=0.6.

## Epithelial IFN / MHC / TJ Q4 vs Q1

**No-go.** User lock: Q4 vs Q1 only if n_mice≥6. Here n_mice = **4**. Quartile split on 4 mice is 1 vs 1 and is not a test. Mouse-level IFN / MHC / TJ means are in the table above; they are not a high-vs-low claim.

---

## Epithelium / cell class (honest)

GEO deposits **no** author `Malignant` / `Epithelial` column on either series. Calling is markers on the **integrated** object:

- **Tight epithelium (primary)** = Epcam+ **and** (Cdh1 or Krt8/18/19 or Cldn18)+ **and** Ptprc−.
- **T/NK** = Cd3d / Cd3e / Cd3g / Cd8a / Nkg7 / Ncr1 / Klrb1c. Tight epithelium wins if both fire.
- **Cldn4 is not a caller.** Tacstd2 is inventory only. No dual-high gate.
- GSE180963 has heavy Sftpc ambient in the sibling digest; Sftpc is an ambient flag, not host AT2.

Lineage genes present: Epcam, Cdh1, Krt8, Krt18, Krt19, Sftpc, Scgb1a1, Ager, Cd3d, Cd3e, Cd3g, Cd8a, Nkg7, Ncr1, Klrb1c, Cldn4, Ptprc, Stk11, Tacstd2, Cldn18, Nkx2-1, Cldn1, Cldn3, Cldn7, Cldn18, Tjp1, Tjp2, Tjp3, Ocln, F11r, Marveld2, Marveld3, Cgn, Crb3.
Missing: none of the core set.

---

## Methods (short)

1. Public only. Downloaded per-GSM 10x MTX from NCBI GEO. No SRA / FASTQ. No private 8-KL object.
2. **GSE165641:** Cell Ranger filtered MTX, 2 KL mice (C57BL/6), Ad-Cre 10 weeks. GSM5047303 SOFT `Lkb2fl/fl` treated as a typo for Lkb1 (series title + PMID).
3. **GSE180963:** author-QC MTX (Seurat 3.1.5: 500–6000 features, mito < 20%). 1 K + 1 KL (FVB). Sibling digest: the two samples were mixed in one 10x library and demultiplexed by label.
4. `CreateSeuratObject` **per dataset** (mice as metadata), then `merge`. QC: nFeature≥200, nCount≥500, percent.mt<25.
5. Integration: split RNA layers by `dataset`, NormalizeData / VST / ScaleData / PCA, then `Seurat::IntegrateLayers(method = HarmonyIntegration)` (harmony 2.0.5). UMAP on the Harmony reduction. Scores stay on log-normalized counts (not Harmony-corrected expression).
6. Cldn4-only. Tacstd2 is inventory. No dual-high. TJ module **excludes** Cldn4.
7. Positive = raw count > 0. Module scores = mean lognorm of present genes.
8. Unit = mouse after merge. Q4 vs Q1 only if n_mice≥6.
9. Thesis is not rewritten.

```bash
bash methods/seurat_integrate_165641_180963_cldn4/scripts/install_r.sh
bash methods/seurat_integrate_165641_180963_cldn4/scripts/download.sh /tmp/geo/work
Rscript methods/seurat_integrate_165641_180963_cldn4/scripts/analyze.R --data /tmp/geo/work --out methods/seurat_integrate_165641_180963_cldn4
```

---

## How to read this

- **This is an integration job.** Separate GSE165641 and GSE180963 catalogs already exist; they are not repeated here as the result.
- **Honest n is 4 mice after merge.** Not 20953 cells. Not 2 series.
- **Cldn4 vs T/NK is four (or fewer) points.** Report the table. Do not write a law.
- **Q4 vs Q1 IFN/MHC/TJ did not fire** unless n_mice≥6.
- **Strain / lab batch remains.** Harmony is for the UMAP, not a license to ignore C57BL/6 vs FVB.
- **No dual-high. No private 8 KL. No FACS-only / CD45-only / subQ / injury-KO merge. Thesis unchanged.**

## Files

- `tables/per_mouse.tsv` / `tables/mouse_level_scores.tsv` — unit-level table
- `tables/mouse_level_associations.tsv` — Cldn4 vs T/NK and IFN/MHC/TJ Spearman
- `tables/q4q1_ifn_mhc_tj.tsv` — no-go or results
- `tables/cldn4_positive_cells.tsv`, `tables/compartment_by_mouse.tsv`
- `tables/gene_inventory.tsv`, `tables/honest_n.tsv`, `tables/series_status.tsv`, `tables/summary.json`
- `figures/fig_umap_dataset.png`, `fig_umap_cellclass.png` — required DimPlots
- `figures/fig_cldn4_vs_tnk.png`, `fig_cldn4_epi.png`, `fig_epi_ifn_mhc_tj.png`, `fig_honest_n.png`
- `objects/integrated_gse165641_gse180963.rds` — integrated Seurat object (gitignored if large)
- `scripts/analyze.R` — R + Seurat primary

## 结论

GSE165641 与 GSE180963 都有公开 10x MTX，已用 **R + Seurat `CreateSeuratObject`（每个 dataset）+ Harmony `IntegrateLayers`** 做成 **一个** 整合对象。诚实 n = **4 只鼠**（合并后），不是 20953 个细胞。Cldn4 对 T/NK 在鼠单位上已计分；上皮 IFN/MHC/TJ 的 Q4 vs Q1 因 n<6 未做。未打开私有 8 只 KL。无 dual-high，不改 thesis。

