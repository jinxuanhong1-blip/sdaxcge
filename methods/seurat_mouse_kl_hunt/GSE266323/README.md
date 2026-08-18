# GSE266323 Seurat Cldn4-only scores

KP LUAD TME leftover. Malat1 CRISPRa n=2 vs Tomato n=2.

Unit: one row per biological mouse.
Cldn4-only. Processed GEO matrix. No FASTQ.
Seurat NormalizeData + AddModuleScore. Epithelial = Epcam>0 or (Krt8>0 and Ptprc==0).
T/NK cells = Cd3e/Cd3d/Nkg7/Ncr1 > 0. Mouse-level Cldn4 is mean log-normalized Cldn4 in epithelial cells (require >=10 epi cells).
