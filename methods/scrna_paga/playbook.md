# Methods: PAGA / diffusion on public LUAD epithelium (TACSTD2 / CLDN4)

中文与英文并列。This playbook is **additive**. It does not replace `methods/scrna/` (ICI scRNA stack) or `methods/trajectory/` (Palantir / CellRank 2 / velocity gates). Those documents stay frozen.

**Question.** On a **public** LUAD epithelial object, after a go/no-go gate, where do **`TACSTD2`** and **`CLDN4`** sit on **PAGA** connectivities and **diffusion** pseudotime along AT2 / club / basal / malignant-like / barrier-keratin programs?

**Not a lineage proof.** PAGA is a coarse-grained graph of cluster connectivities (Wolf et al., *Genome Biol* 2019). Diffusion maps / DPT are an ordering (Haghverdi et al., 2015; 2016). Neither proves that AT2 becomes tumor. `TACSTD2` and `CLDN4` are **readouts on the graph**, not roots, not terminals, and not a both-high gate.

---

## 范围 / Scope

**In scope**

- One public human LUAD epithelial slice. Preferred accession: **GSE131907** (Kim et al., *Nat Commun* 2020) **tLung** epithelium, with **nLung** epithelium from the same atlas used only as an AT2 / club / AT1 **prior** and as an external root (not as a second cancer type).
- Alternatives named by the task, used only if GSE131907 is unavailable: GSE253013 (skip if the processed object exceeds the public-file size cap) or GSE207422 **LUAD malignant epithelium only** (do not pool LUSC).
- Author cell-type labels when deposited. Marker **scores** for AT2, club, basal, ciliated, transitional/barrier-keratin, and malignant-like programs — computed **after** restriction, not used to redefine the download.
- PAGA + diffusion maps + DPT on the **epithelial** neighbor graph. Gene-on-axis tests for `TACSTD2` and `CLDN4` separately, plus their Spearman co-expression.
- An **extra** barrier/keratin figure **only if** the pre-specified sample-level test says `TACSTD2` tracks that program (see Locked design). A negative test is written as negative; the figure is then omitted.

**Out of scope**

- Palantir, CellRank 2, scVelo / veloVI streamlines (see `methods/trajectory/`).
- Whole-tissue UMAPs that mix immune / stroma / epithelium.
- GSE253013 9 GB RDS, EGA/dbGaP FASTQ, or any file over the ~2 GB processed-public cap.
- Pooling LUAD with LUSC. GSE207422 is NSCLC; if used, histology must be LUAD or the slice stops.
- A TACSTD2∩CLDN4 “both-high” responder or lineage gate.
- ICI / MPR / RECIST claims from GSE131907 (treatment-naive atlas).
- Tuning Leiden resolution, root cells, or gene lists after seeing `TACSTD2` trends.

---

## 锁定设计 / Locked design (before looking at DPT)

| 项目 Item | 选择 Choice | 诚实限制 Limitation |
| --- | --- | --- |
| Cohort | GSE131907 LUAD; epithelium with `Sample_Origin` ∈ {tLung, nLung} | tLung is the tumor object. nLung is the alveolar/airway prior. Mets / PE / LN are not in this slice. |
| Counts | Author **raw UMI** text matrix (public, under the size cap). Not the 2.9 GB log2TPM text. | No re-alignment. No CellBender on the deposited matrix. |
| Restriction | `Cell_type == Epithelial cells` before neighbors | Immune leak is a QC check (`PTPRC`), not a second trajectory. |
| Graph | Seurat-v3 HVG on counts → PCA → k-NN (`n_neighbors=30`, `n_pcs=30`) **recomputed on epithelium** | Do not inherit the whole-atlas graph. |
| Clusters | Leiden on that graph; PAGA on Leiden | Resolution is pre-set (0.6). Do not retune to make a pretty path. |
| Root (DPT) | nLung cells in the Leiden cluster with highest **AT2 score**, or nLung `Cell_subtype == AT2` if that label exists | Never root on `TACSTD2` / `CLDN4` high. |
| Terminals | Author tLung tumor states (tS1/tS2/tS3) and/or highest malignant-like score | `TACSTD2` high is **not** malignancy. |
| Airway | Club / basal / ciliated are **scored and drawn on PAGA**. They are **not** stitched onto AT2 if the k-NN graph is disconnected | A disconnected airway component is reported as discrete states, not forced DPT. |
| Genes | `TACSTD2` and `CLDN4` tested **separately**. Combinatorial = Spearman co-expression along DPT / by cluster, not a quadrant gate | No high-high filter. |
| Inferential n | **Sample** (and patient, if it differs) is the unit for Spearman of mean gene vs mean DPT / mean module | Cell-level ρ is descriptive only. Huge n_cells must not be sold as a p-value. |
| Extra figure | Emit `fig_extra_barrier_keratin` iff sample-level Spearman(`TACSTD2`, barrier-keratin score) has ρ > 0 and two-sided p < 0.05 | If not, write the numbers and skip the figure. |
| Histology | LUAD only | Do not add LUSC from another accession to “complete” the axis. |

---

## 状态语法 / State grammar (pre-specified scores)

Scores are the mean of log-normalized expression of **detected** genes in the set (genes absent from the matrix are recorded as absent, not replaced).

| State | Positive RNA | Role |
| --- | --- | --- |
| AT2 | `SFTPC`, `SFTPB`, `SFTPA1`, `NAPSA`, `LAMP3`, `ABCA3` | Default **source** / early |
| AT1 | `AGER`, `PDPN`, `CAV1` | Repair terminal (not tumor) |
| Club | `SCGB1A1`, `SCGB3A2`, `SCGB3A1` | Airway secretory; often `TACSTD2` high at baseline |
| Basal | `KRT5`, `KRT15`, `TP63`, `NGFR` | Airway; often `TACSTD2` high at baseline |
| Ciliated | `FOXJ1`, `TPPP3`, `PIFO` | Off-axis check |
| Barrier / keratin (transitional, KAC-like) | `KRT8`, `KRT18`, `KRT19`, `KRT7`, `CLDN4`, `CDKN1A`, `PLAUR` | Waypoint / injury-keratin program. `CLDN4` is **inside** this set; `TACSTD2` is **not** |
| Malignant-like | Author tS1/tS2/tS3 when present; else `CEACAM5`, `CEACAM6`, `MKI67` as a weak expression proxy | Expression proxy **is not CNV**. Say so. |
| Lineage QC | `EPCAM` high; `PTPRC`, `PECAM1`, `COL1A1` low | Fail the slice if epithelium is not epithelial |

`TACSTD2` is airway-high at baseline (club / basal) and can rise in injured or neoplastic alveolar epithelium. A `TACSTD2` gradient on a mixed epithelial graph is often **airway vs alveolar**, not AT2 → tumor. Always pair it with AT2 genes, `CLDN4` + `KRT8`, and author tumor states.

---

## PAGA / 扩散 / Diffusion procedure

1. QC on the epithelial subset: genes in ≥ 10 cells; cells with a minimum UMI/gene floor written in the run log (do not tune after seeing `TACSTD2`).
2. `normalize_total` to 10⁴, `log1p`. Store raw counts in `layers["counts"]`.
3. Highly variable genes (Seurat v3, 3 000) on counts. Scale, PCA (50), neighbors (30/30).
4. Leiden (resolution 0.6). `sc.tl.paga`. `sc.tl.diffmap` (15 components). `sc.tl.dpt` with the locked AT2 root.
5. UMAP is **display only**. Neighbors are never computed on `X_umap`. PAGA connectivities, not UMAP edges, are the graph claim.
6. If PAGA shows two or more components that separate nLung AT2 from tLung tumor states **and** from club/basal, **do not** glue them. Report discrete PAGA vertices and skip a single global DPT interpretation (still compute DPT inside the AT2-containing component).

Wolf et al. 2019: PAGA threshold is a visualization parameter. Export the **full** connectivity matrix. Do not pick a threshold that creates an AT2→tumor path after seeing the genes.

---

## 基因在轴上 / Genes on the axis

For each of `TACSTD2` and `CLDN4`:

1. Cluster-level: mean log1p, percent UMI>0, n_cells, n_samples contributing ≥ 10 cells.
2. DPT: cell-level Spearman (descriptive) and **sample-level** Spearman of sample-mean gene vs sample-mean DPT (primary). Report n_samples and n_patients.
3. Module Spearman (sample-level): gene vs AT2, club, basal, barrier-keratin, malignant-like scores.
4. Co-expression: sample-level Spearman(`TACSTD2`, `CLDN4`). Not a both-high gate.

Control panel (same tests, same FDR family): `SFTPC`, `KRT8`, `SCGB1A1`, `KRT5`, `EPCAM`, `PTPRC`.

BH-FDR is applied only inside this pre-specified list.

---

## 功效与措辞 / Power and wording

- GSE131907 tLung / nLung epithelial **sample** counts are small. Write **n** next to every ρ and every PAGA vertex.
- Underpowered sample-level tests are **inconclusive**, not “no trajectory.”
- Do not write “AT2 differentiates into LUAD because PAGA connects tS2.”
- Do not write “TROP2 marks the malignant terminal.”
- Do not interpret GSE131907 as ICI biology.

---

## 复现 / Reproduce

```bash
pip install -r requirements.txt
bash scripts/scrna_paga/download.sh /tmp/scrna_paga_data
python3 scripts/scrna_paga/extract_epithelium.py --data /tmp/scrna_paga_data --out /tmp/scrna_paga_data/epithelium.h5ad
python3 scripts/scrna_paga/analyze_paga.py --input /tmp/scrna_paga_data/epithelium.h5ad --outdir results/scrna_paga
```

Public files only (NCBI GEO FTP). No FASTQ.

---

## 中文版

本手册是**加法**：不改写 `methods/scrna/` 与 `methods/trajectory/`。对象是公开 **LUAD** 上皮（首选 GSE131907 tLung，nLung 仅作 AT2/club 先验与 DPT 根）。工具是 **PAGA + 扩散伪时间**，不是 Palantir / CellRank / velocity。`TACSTD2` 与 `CLDN4` **分别**检验，共表达用 Spearman，不用双高门控。根细胞必须是 nLung AT2，不能是 TROP2 高细胞。若 k-NN 上 club/basal 与 AT2 不连通，如实写离散状态，不强行缝合。样本为推断单位；细胞数很大时的 p 值不当作主结论。仅当样本水平 `TACSTD2` 与屏障/角蛋白程序正相关（ρ>0 且 p<0.05）时才出额外图。GSE131907 不能检验 ICI。
