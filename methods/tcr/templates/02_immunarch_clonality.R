#!/usr/bin/env Rscript
# immunarch >= 0.10: load 10x filtered_contig_annotations and write
# clonality / diversity tables. Put metadata.txt next to the contig files.
#
#   Rscript 02_immunarch_clonality.R data/public/gse243013/tcr workdir/immunarch

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("usage: 02_immunarch_clonality.R <tcr_folder> <out_dir>")
}
tcr_dir <- args[[1]]
out_dir <- args[[2]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("immunarch", quietly = TRUE)) {
  stop("install immunarch (>= 0.10) from CRAN before running this template")
}

imm <- immunarch::repLoad(tcr_dir, .mode = "paired")

vol <- immunarch::repExplore(imm$data, .method = "volume")
homeo <- immunarch::repClonality(imm$data, .method = "homeo")
cprop <- immunarch::repClonality(imm$data, .method = "clonal.prop", .perc = 10)
chao1 <- immunarch::repDiversity(imm$data, .method = "chao1")
invsimp <- immunarch::repDiversity(imm$data, .method = "inv.simp")
d50 <- immunarch::repDiversity(imm$data, .method = "d50")

write.table(as.data.frame(vol), file.path(out_dir, "explore_volume.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(as.data.frame(homeo), file.path(out_dir, "clonality_homeo.tsv"),
            sep = "\t", quote = FALSE, row.names = TRUE)
write.table(as.data.frame(cprop), file.path(out_dir, "clonality_prop.tsv"),
            sep = "\t", quote = FALSE, row.names = TRUE)
write.table(as.data.frame(chao1), file.path(out_dir, "diversity_chao1.tsv"),
            sep = "\t", quote = FALSE, row.names = TRUE)
write.table(as.data.frame(invsimp), file.path(out_dir, "diversity_invsimp.tsv"),
            sep = "\t", quote = FALSE, row.names = TRUE)
write.table(as.data.frame(d50), file.path(out_dir, "diversity_d50.tsv"),
            sep = "\t", quote = FALSE, row.names = TRUE)

if (length(imm$data) >= 2) {
  ov <- immunarch::repOverlap(imm$data, .method = "jaccard", .col = "aa")
  write.table(as.data.frame(as.matrix(ov)), file.path(out_dir, "overlap_jaccard.tsv"),
              sep = "\t", quote = FALSE, row.names = TRUE)
}

if (!is.null(imm$meta)) {
  write.table(imm$meta, file.path(out_dir, "metadata_used.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

message("wrote ", out_dir)
