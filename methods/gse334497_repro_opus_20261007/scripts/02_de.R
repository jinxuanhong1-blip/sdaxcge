#!/usr/bin/env Rscript
# Differential expression for GSE334497 (4T1 Trop2-KO vs Trop2-WT tumours, BALB/c).
# Contrast everywhere: KO vs WT (positive log2FC = higher in Trop2-KO).
suppressPackageStartupMessages({
  library(DESeq2); library(edgeR); library(limma)
})

args <- commandArgs(trailingOnly = FALSE)
here <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
raw_dir <- file.path(here, "raw"); tab <- file.path(here, "tables")
dir.create(tab, showWarnings = FALSE)

counts <- as.matrix(read.csv(gzfile(file.path(raw_dir, "GSE334497_raw_counts_recovered.csv.gz")),
                             row.names = 1, check.names = FALSE))
sheet <- read.csv(file.path(tab, "sample_sheet.csv"), stringsAsFactors = FALSE)
rownames(sheet) <- sheet$library_name
sheet <- sheet[colnames(counts), ]
ann <- read.delim(gzfile(file.path(raw_dir, "gene_annotation_ensembl102.tsv.gz")), row.names = 1)

source_panels <- function() {
  py <- readLines(file.path(here, "scripts", "gene_sets.py"))
  txt <- paste(py, collapse = "\n")
  get <- function(name) {
    m <- regmatches(txt, regexpr(paste0("\\n", name, " = \\[[^]]*\\]"), txt))
    unlist(regmatches(m, gregexpr('"[^"]+"', m))) |> gsub(pattern = '"', replacement = "")
  }
  names_map <- c(APM_core = "APM_CORE", IFN_core = "IFN_CORE", Immune_core = "IMMUNE_CORE",
                 Claudin_1_3_4_7 = "CLAUDIN_FOCAL", Paper_TJ_trio_Cldn1_Cldn7_Ocln = "PAPER_TJ_TRIO",
                 TJ_claudin_program = "TJ_CLAUDIN_PROGRAM", CD8_T = "CD8_T", Effector = "EFFECTOR",
                 Leukocyte = "LEUKOCYTE", NK = "NK", MHC_II = "MHC_II", Checkpoint = "CHECKPOINT")
  lapply(names_map, get)
}
panels <- source_panels()

designs <- list(
  M1_all10_genotype          = list(keep = rownames(sheet), formula = ~ genotype),
  M2_all10_wave_plus_genotype = list(keep = rownames(sheet), formula = ~ library_wave + genotype),
  M3_drop_RESUB170R_outlier  = list(keep = setdiff(rownames(sheet), "RESUB-170R"), formula = ~ genotype),
  M4_drop_control170         = list(keep = setdiff(rownames(sheet), "control170"), formula = ~ genotype),
  # Post hoc: drop the three WT tumours carrying epidermis (SKIN_FLAG_GENES in gene_sets.py;
  # see 07_tissue_contamination.py). Leaves 2 WT (both RESUB wave) vs 5 KO.
  M5_drop_skin_positive_WT   = list(keep = setdiff(rownames(sheet), c("control170", "RESUB-170R", "RESUB-169R")),
                                    formula = ~ genotype)
)

run_deseq <- function(keep, formula) {
  cd <- sheet[keep, ]
  cd$genotype <- factor(cd$genotype, levels = c("WT", "KO"))
  cd$library_wave <- factor(cd$library_wave, levels = c("original", "RESUB"))
  m <- counts[, keep]
  m <- m[rowSums(m) > 0, ]
  dds <- DESeqDataSetFromMatrix(m, cd, formula)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, name = "genotype_KO_vs_WT", alpha = 0.05)
  shr <- lfcShrink(dds, coef = "genotype_KO_vs_WT", type = "normal", quiet = TRUE)
  out <- data.frame(ensembl_gene_id = rownames(res), gene_name = ann[rownames(res), "gene_name"],
                    baseMean = res$baseMean, log2FC_KO_vs_WT = res$log2FoldChange, lfcSE = res$lfcSE,
                    log2FC_KO_vs_WT_shrunk_normal = shr$log2FoldChange,
                    wald_stat = res$stat, pvalue = res$pvalue, padj_BH = res$padj,
                    n_WT = sum(cd$genotype == "WT"), n_KO = sum(cd$genotype == "KO"),
                    stringsAsFactors = FALSE)
  out[order(out$pvalue), ]
}

all_res <- list()
for (nm in names(designs)) {
  d <- designs[[nm]]
  r <- run_deseq(d$keep, d$formula)
  r$model <- nm
  all_res[[nm]] <- r
  write.csv(r, gzfile(file.path(tab, sprintf("de_deseq2_%s.csv.gz", nm))), row.names = FALSE)
  cat(sprintf("%s: n=%d genes tested; padj<0.05: %d (up in KO %d, down in KO %d)\n", nm, nrow(r),
              sum(r$padj_BH < 0.05, na.rm = TRUE),
              sum(r$padj_BH < 0.05 & r$log2FC_KO_vs_WT > 0, na.rm = TRUE),
              sum(r$padj_BH < 0.05 & r$log2FC_KO_vs_WT < 0, na.rm = TRUE)))
}

# edgeR quasi-likelihood and limma-voom on the primary model (all 10, ~genotype)
genotype <- factor(sheet$genotype, levels = c("WT", "KO"))
design <- model.matrix(~ genotype)
y <- DGEList(counts, genes = data.frame(gene_name = ann[rownames(counts), "gene_name"]))
keep_expr <- filterByExpr(y, design)
y_f <- calcNormFactors(y[keep_expr, , keep.lib.sizes = FALSE])
y_f <- estimateDisp(y_f, design)
qfit <- glmQLFit(y_f, design, robust = TRUE)
qres <- topTags(glmQLFTest(qfit, coef = "genotypeKO"), n = Inf, sort.by = "none")$table
v <- voom(y_f, design)
vfit <- eBayes(lmFit(v, design), robust = TRUE)
vres <- topTable(vfit, coef = "genotypeKO", n = Inf, sort.by = "none")

# Low-expression genes (Tacstd2 included) are dropped by filterByExpr; report them from an
# unfiltered fit too so that every focal gene has an edgeR/limma estimate.
y_u <- calcNormFactors(y[rowSums(counts) > 0, , keep.lib.sizes = FALSE])
y_u <- estimateDisp(y_u, design)
qres_u <- topTags(glmQLFTest(glmQLFit(y_u, design, robust = TRUE), coef = "genotypeKO"),
                  n = Inf, sort.by = "none")$table
v_u <- voom(y_u, design)
vres_u <- topTable(eBayes(lmFit(v_u, design), robust = TRUE), coef = "genotypeKO", n = Inf, sort.by = "none")

alt <- data.frame(ensembl_gene_id = rownames(y_u), gene_name = y_u$genes$gene_name,
                  passes_filterByExpr = rownames(y_u) %in% rownames(y_f),
                  edgeR_QL_log2FC_unfiltered = qres_u$logFC, edgeR_QL_p_unfiltered = qres_u$PValue,
                  limma_voom_log2FC_unfiltered = vres_u$logFC, limma_voom_p_unfiltered = vres_u$P.Value,
                  stringsAsFactors = FALSE)
alt$edgeR_QL_log2FC <- qres$logFC[match(alt$ensembl_gene_id, rownames(qres))]
alt$edgeR_QL_p <- qres$PValue[match(alt$ensembl_gene_id, rownames(qres))]
alt$edgeR_QL_FDR <- qres$FDR[match(alt$ensembl_gene_id, rownames(qres))]
alt$limma_voom_log2FC <- vres$logFC[match(alt$ensembl_gene_id, rownames(vres))]
alt$limma_voom_p <- vres$P.Value[match(alt$ensembl_gene_id, rownames(vres))]
alt$limma_voom_FDR <- vres$adj.P.Val[match(alt$ensembl_gene_id, rownames(vres))]
write.csv(alt, gzfile(file.path(tab, "de_edgeR_limma_M1.csv.gz")), row.names = FALSE)
cat(sprintf("edgeR QL (filterByExpr, %d genes): FDR<0.05 = %d; limma-voom: FDR<0.05 = %d\n",
            nrow(qres), sum(qres$FDR < 0.05), sum(vres$adj.P.Val < 0.05)))

# Panel tests on voom (filterByExpr genes). fry = self-contained rotation test.
# camera = competitive test; limma's default fixes inter-gene correlation at 0.01,
# which is anti-conservative for co-regulated panels, so both are reported.
sym <- y_f$genes$gene_name
idx <- lapply(panels, function(g) which(sym %in% g))
fr <- fry(v, index = idx, design = design, contrast = "genotypeKO", sort = "none")
cm <- camera(v, index = idx, design = design, contrast = "genotypeKO", sort = FALSE, inter.gene.cor = NA)
cm_fixed <- camera(v, index = idx, design = design, contrast = "genotypeKO", sort = FALSE)
panel_tab <- data.frame(panel = names(panels),
                        n_genes_defined = sapply(panels, length),
                        n_genes_tested = sapply(idx, length),
                        genes_tested = sapply(idx, function(i) paste(sym[i], collapse = ";")),
                        mean_voom_log2FC_KO_vs_WT = sapply(idx, function(i) mean(vres$logFC[i])),
                        fry_direction = fr$Direction, fry_p = fr$PValue, fry_p_mixed = fr$PValue.Mixed,
                        camera_direction = cm$Direction, camera_p_estimated_cor = cm$PValue,
                        camera_estimated_intergene_cor = cm$Correlation,
                        camera_p_default_cor0.01 = cm_fixed$PValue,
                        stringsAsFactors = FALSE)
panel_tab$fry_FDR <- p.adjust(panel_tab$fry_p, "BH")
panel_tab$camera_FDR_estimated_cor <- p.adjust(panel_tab$camera_p_estimated_cor, "BH")
write.csv(panel_tab, file.path(tab, "panel_tests_limma_fry_camera.csv"), row.names = FALSE)
print(panel_tab[, c("panel", "n_genes_tested", "mean_voom_log2FC_KO_vs_WT", "fry_direction", "fry_p",
                    "camera_direction", "camera_p_estimated_cor", "camera_estimated_intergene_cor",
                    "camera_p_default_cor0.01")], row.names = FALSE, digits = 3)

# Focal genes across all models/methods
focal_src <- readLines(file.path(here, "scripts", "gene_sets.py")) |> paste(collapse = "\n")
fm <- regmatches(focal_src, regexpr("FOCAL_GENES = \\[[^]]*\\]", focal_src))
focal <- gsub('"', "", unlist(regmatches(fm, gregexpr('"[^"]+"', fm))))
all_cldn <- sort(unique(ann$gene_name[grepl("^Cldn[0-9]", ann$gene_name)]))
focal <- unique(c(focal, all_cldn))
rows <- list()
for (g in focal) {
  ids <- rownames(ann)[which(ann$gene_name == g)]
  if (length(ids) == 0) {
    rows[[length(rows) + 1]] <- data.frame(gene_name = g, ensembl_gene_id = NA, in_matrix = FALSE)
    next
  }
  for (id in ids) {
    row <- data.frame(gene_name = g, ensembl_gene_id = id, in_matrix = id %in% rownames(counts))
    for (nm in names(all_res)) {
      r <- all_res[[nm]][all_res[[nm]]$ensembl_gene_id == id, ]
      row[[paste0(nm, "_baseMean")]] <- if (nrow(r)) r$baseMean else NA
      row[[paste0(nm, "_log2FC")]] <- if (nrow(r)) r$log2FC_KO_vs_WT else NA
      row[[paste0(nm, "_lfcSE")]] <- if (nrow(r)) r$lfcSE else NA
      row[[paste0(nm, "_p")]] <- if (nrow(r)) r$pvalue else NA
      row[[paste0(nm, "_padj")]] <- if (nrow(r)) r$padj_BH else NA
    }
    a <- alt[alt$ensembl_gene_id == id, ]
    for (cc in c("passes_filterByExpr", "edgeR_QL_log2FC_unfiltered", "edgeR_QL_p_unfiltered",
                 "limma_voom_log2FC_unfiltered", "limma_voom_p_unfiltered", "edgeR_QL_FDR", "limma_voom_FDR"))
      row[[cc]] <- if (nrow(a)) a[[cc]] else NA
    rows[[length(rows) + 1]] <- row
  }
}
focal_tab <- do.call(rbind, lapply(rows, function(r) {
  full <- setNames(rep(NA, length(names(rows[[which.max(sapply(rows, ncol))]]))),
                   names(rows[[which.max(sapply(rows, ncol))]]))
  full[names(r)] <- unlist(r[1, ]); as.data.frame(as.list(full), stringsAsFactors = FALSE, check.names = FALSE)
}))
write.csv(focal_tab, file.path(tab, "focal_genes_all_models.csv"), row.names = FALSE)

vsd <- vst(DESeqDataSetFromMatrix(counts[rowSums(counts) > 0, ], transform(sheet, genotype = factor(genotype)), ~ genotype), blind = TRUE)
write.csv(assay(vsd), gzfile(file.path(tab, "vst_blind.csv.gz")))
writeLines(capture.output(sessionInfo()), file.path(tab, "sessionInfo_R.txt"))
