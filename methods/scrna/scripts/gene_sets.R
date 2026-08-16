## Curated gene sets for the lung-cancer ICI scRNA playbook (R).
## Starting points, not gospel; cite the source when used. See playbook.md 6-7.

gene_sets <- list(
  TACSTD2_CLDN4_junction = c("TACSTD2","CLDN4","CLDN3","CLDN7","CLDN18",
                             "CDH1","TJP1","OCLN","F11R","EPCAM","ELF3"),
  epithelial_lineage     = c("EPCAM","KRT8","KRT18","KRT19","CDH1","SFTPC","SCGB1A1"),
  CD8_cytotoxic          = c("CD8A","CD8B","GZMB","GZMK","GZMH","PRF1","IFNG","NKG7","KLRG1"),
  CD8_exhaustion         = c("PDCD1","HAVCR2","LAG3","TIGIT","CTLA4","TOX","ENTPD1"),
  TRM                    = c("ITGAE","ZNF683","CXCR6","ITGA1"),
  TLS_chemokine          = c("CCL2","CCL3","CCL4","CCL5","CCL8","CCL18","CCL19",
                             "CCL21","CXCL9","CXCL10","CXCL11","CXCL13"),
  TLS_core               = c("CXCL13","CCL19","CCL21","CR2","CXCR5","LTB","SELL","MS4A1"),
  CXCL13                 = c("CXCL13")
)

lineage_markers <- list(
  `T/NK`       = c("CD3D","CD3E","CD2","TRAC","NKG7","GNLY","KLRD1"),
  `B/Plasma`   = c("MS4A1","CD79A","CD79B","BANK1","MZB1","IGHG1","XBP1"),
  Myeloid      = c("LYZ","CD68","CD14","FCGR3A","C1QA","C1QB","SPP1"),
  DC           = c("CLEC9A","CD1C","LAMP3","LILRA4"),
  Mast         = c("TPSAB1","TPSB2","CPA3","MS4A2"),
  Epithelial   = c("EPCAM","KRT8","KRT18","KRT19","CDH1","SFTPC","SCGB1A1"),
  Endothelial  = c("PECAM1","VWF","CLDN5","CLEC14A"),
  Fibroblast   = c("COL1A1","COL1A2","DCN","LUM","PDGFRB","ACTA2")
)
