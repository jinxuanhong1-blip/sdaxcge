# GSE201247 Seurat Cldn4-only scores

Kras vs ATTAC whole-lung leftover. Honest n=6 mice (2 WT ATTAC, 2 Kras, 2 ATTAC;Kras).

Unit: one row per biological mouse.
Cldn4-only. Processed GEO matrix. No FASTQ.
Seurat NormalizeData + AddModuleScore. Epithelial = Epcam>0 or (Krt8>0 and Ptprc==0).
T/NK cells = Cd3e/Cd3d/Nkg7/Ncr1 > 0. Mouse-level Cldn4 is mean log-normalized Cldn4 in epithelial cells (require >=10 epi cells).
