# Methods — GSE189357 REAL Slingshot/PAGA, CLDN4 only

ADDITIVE. Does **not** rewrite winning-pair Slingshot (`methods/winpair_131907_205335_slingshot_cldn4/`) or GSE131907-only PAGA.
TACSTD2 does not define groups. No dual-high gate. Spatial GSE189487 is not used.

**Question.** On GSE189357 **tumor epithelium only**, where do **CLDN4**, a CLDN4-excluded barrier/keratin program, and Hallmark IFN sit on a Slingshot lineage whose root is **not** CLDN4-high?

## Dataset

| Item | Value |
| --- | --- |
| Accession | GSE189357 (Zhu et al., *Exp Mol Med* 2022, DOI 10.1038/s12276-022-00896-9) |
| Design | 9 treatment-naïve resected LUAD lesions; 3 AIS + 3 MIA + 3 IAC |
| Patients | TD1–TD9. **n may be 9 — say so.** Patient = unit (one 10x sample each). |
| Stages | TD5/TD7/TD8 AIS; TD3/TD4/TD6 MIA; TD1/TD2/TD9 IAC (GEO `histolgical type`) |
| Matrix | public GEO 10x MTX/TSV in `GSE189357_RAW.tar` (no FASTQ, no SRA) |
| Epithelium | marker-malignant: `(EPCAM\|KRT8\|KRT18\|KRT19)>0` and `PTPRC==0` |
| Unused | GSE189487 ST; dual-high TACSTD2∩CLDN4; nLung (none exists) |

There is **no uninvolved lung** in this accession. Rooting uses AT2-high **tumor** epithelium.

## Trajectory clock

**REAL Slingshot** (Street et al. 2018; Bioconductor `slingshot`) is required. Tools are installed by `scripts/install_tools.sh`. PAGA (Wolf et al. 2019) is geometry only. scanpy DPT is a companion ordering, not the primary clock.

| Item | Choice |
| --- | --- |
| Graph | Seurat-v3 HVG 3000 → PCA → Harmony (`batch = patient`) → k-NN 30 |
| Clusters | Leiden 0.6; PAGA on Leiden |
| Root | Leiden with highest AT2 among clusters that are **not** the CLDN4-highest cluster. Root cell = median AT2 in that cluster, **excluding CLDN4-high tertile**. Never CLDN4-high. |
| Slingshot | `slingshot(rd, clusterLabels, start.clus=root)` on Harmony PCs |
| Score | CLDN4 continuous + tertiles; barrier/keratin **without CLDN4**; Hallmark IFNα ∪ IFNγ |
| Inferential n | **patient**. Cell-level ρ is descriptive. n may be 9. |
| Extra figures | CLDN4 + barrier + IFN along Slingshot; within-patient CLDN4-high vs low |
| Unused | dual-high TACSTD2∩CLDN4; stage as a powered contrast (n=3/stage) |

Harmony on patient removes between-patient / between-stage mean shifts (each patient is one stage). The remaining axis is within-epithelium programs. Do not read Slingshot as AIS→IAC invasion time.

## Tests

Primary: patient-level Spearman of mean CLDN4 / barrier / IFN vs mean Slingshot PT. BH inside that list only.

Sensitivity (not BH): drop IAC, drop AIS, DPT companion, IFN core, CLDN4-ward lineage.

Paired extra: within-patient CLDN4-high vs low program scores (min 8 cells/arm). n may be 9.

## Done criterion

`results/tables/slingshot_lineages.tsv` exists after a successful REAL Slingshot run.
