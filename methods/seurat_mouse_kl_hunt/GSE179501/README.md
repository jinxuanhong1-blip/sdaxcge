# GSE179501 Seurat Cldn4-only scores

Lkb1-XTR total-viable leftover (sister of assigned GSE179502). Restored n=2 vs Non-Restored n=2.

Unit: one row per biological mouse.
Cldn4-only. Processed GEO matrix. No FASTQ.
Seurat NormalizeData + AddModuleScore. Epithelial = Epcam>0 or (Krt8>0 and Ptprc==0).
T/NK cells = Cd3e/Cd3d/Nkg7/Ncr1 > 0. Mouse-level Cldn4 is mean log-normalized Cldn4 in epithelial cells (require >=10 epi cells).
