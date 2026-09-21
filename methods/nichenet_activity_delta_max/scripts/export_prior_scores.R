#!/usr/bin/env Rscript
# Export NicheNet v2 target shifts on the locked concordant-4 T/NK CLDN4 axis.
# Matrix: Zenodo 7074291 ligand_target_matrix_nsga2r_final.rds
# Gene axis and sender deltas are the author-annotated summaries from the
# concordant-4 quartile extract (GSE131907 + GSE205335 for the receiver axis).

args <- commandArgs(trailingOnly = TRUE)
opt <- function(flag, default) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^", flag, "="), "", hit[[1]]) else default
}
PRIOR <- opt("--prior", "/tmp/nichenet_prior/ligand_target_matrix_nsga2r_final.rds")
INDIR <- opt("--in", "/tmp/nichenet_inputs")
OUT <- opt("--out", "/tmp/nichenet_work/prior_scores")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

message("load matrix")
lt <- readRDS(PRIOR)
stopifnot(is.matrix(lt))

bg <- read.delim(file.path(INDIR, "background_genes.tsv"), stringsAsFactors = FALSE)
bg <- bg$gene[bg$frac >= 0.1]
zdf <- read.delim(file.path(INDIR, "tnk_gene_vs_cldn4.tsv"), stringsAsFactors = FALSE)
zdf$z <- as.numeric(zdf$z)
zdf$epithelial_leak <- as.logical(zdf$epithelial_leak)
zdf <- zdf[is.finite(zdf$z) & !is.na(zdf$epithelial_leak), ]

genes <- intersect(rownames(lt), intersect(bg, zdf$gene[!zdf$epithelial_leak]))
zdf <- zdf[match(genes, zdf$gene), ]
stopifnot(all(zdf$gene == genes))
zz <- zdf$z
message("universe ", length(genes), " mean z ", mean(zz))

sub <- lt[genes, , drop = FALSE]
storage.mode(sub) <- "double"
ligands <- colnames(sub)

# Binary export for the Python AUROC / AUPR grid (column-major).
con <- file(file.path(OUT, "rp_universe.bin"), "wb")
writeBin(as.vector(sub), con)
close(con)
writeLines(genes, file.path(OUT, "genes.txt"))
writeLines(ligands, file.path(OUT, "ligands.txt"))
write.table(
  zdf[, c("gene", "z", "rho", "p", "k", "n_pos", "n_neg")],
  file.path(OUT, "universe_z.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Top-K mean z. K is applied after dropping genes with non-finite z (already done).
Ks <- c(25L, 50L, 100L, 200L)
# rank within each column: a matrix of ranks would be huge. Loop ligands.
shift <- matrix(NA_real_, nrow = length(ligands), ncol = length(Ks),
                dimnames = list(ligands, paste0("mean_z_top", Ks)))
mean_rp <- shift
for (j in seq_along(ligands)) {
  rp <- sub[, j]
  ord <- order(rp, decreasing = TRUE)
  for (ki in seq_along(Ks)) {
    k <- Ks[[ki]]
    top <- ord[seq_len(k)]
    shift[j, ki] <- mean(zz[top])
    mean_rp[j, ki] <- mean(rp[top])
  }
  if (j %% 200 == 0) message("  ligand ", j, "/", length(ligands))
}
out <- data.frame(ligand = ligands, stringsAsFactors = FALSE)
out$mean_rp_all <- colMeans(sub)
for (ki in seq_along(Ks)) {
  out[[paste0("mean_z_top", Ks[[ki]])]] <- shift[, ki]
  out[[paste0("mean_rp_top", Ks[[ki]])]] <- mean_rp[, ki]
}
out$bg_mean_z <- mean(zz)
out$bg_sd_z <- sd(zz)
write.table(out, file.path(OUT, "ligand_target_shift.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("wrote ", nrow(out), " ligands")
