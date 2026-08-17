# Extract IMvigor210 sample-level CLDN4 / CD274 / CD8A / CXCL9 + clinical
# from the official IMvigor210CoreBiologies 1.0.0 CountDataSet.
#
# Does not write the full 31k-gene matrix. TPM uses all genes for the
# per-sample scaling factor, then only the requested rows are kept.
#
# Usage: Rscript extract_cds.R <cds.RData> <out_tsv>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript extract_cds.R <cds.RData> <out_tsv>")
}
rdata <- args[[1]]
out_path <- args[[2]]

load(rdata)
if (!exists("cds")) stop("cds object not found in ", rdata)

ad <- attr(cds, "assayData")
counts <- get("counts", envir = ad)
pd <- attr(attr(cds, "phenoData"), "data")
fd <- attr(attr(cds, "featureData"), "data")

if (!identical(colnames(counts), rownames(pd))) {
  stop("Sample names in counts and phenoData do not match")
}
if (!identical(rownames(counts), rownames(fd))) {
  stop("Feature names in counts and featureData do not match")
}

wanted <- c(CLDN4 = "CLDN4", CD274 = "CD274", CD8A = "CD8A", CXCL9 = "CXCL9")
idx <- integer(length(wanted))
names(idx) <- names(wanted)
for (sym in names(wanted)) {
  hits <- which(fd$symbol == wanted[[sym]])
  if (length(hits) != 1) {
    stop("Expected exactly one ", wanted[[sym]], " row, found ", length(hits))
  }
  idx[[sym]] <- hits
}

len <- as.numeric(fd$length)
if (any(!is.finite(len)) || any(len <= 0)) stop("Invalid gene lengths")
rpk <- counts / len
scale <- colSums(rpk)
if (any(!is.finite(scale)) || any(scale <= 0)) stop("Invalid TPM scale")
tpm <- sweep(rpk, 2, scale, "/") * 1e6

sf <- as.numeric(pd$sizeFactor)
if (any(!is.finite(sf)) || any(sf <= 0)) stop("Invalid sizeFactor")
libsize <- colSums(counts)

pick <- function(sym) {
  i <- idx[[sym]]
  raw <- as.numeric(counts[i, ])
  list(
    entrez = as.character(fd$entrez_id[i]),
    symbol = as.character(fd$symbol[i]),
    length = as.numeric(fd$length[i]),
    raw = raw,
    tpm = as.numeric(tpm[i, ]),
    log2tpm1 = log2(as.numeric(tpm[i, ]) + 1),
    log2sf = log2(raw / sf + 1)
  )
}

g <- lapply(names(wanted), pick)
names(g) <- names(wanted)

na_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == "NA"] <- NA_character_
  x
}

out <- data.frame(
  sample_id = rownames(pd),
  ANONPT_ID = as.character(pd$ANONPT_ID),
  sizeFactor = sf,
  library_size = as.numeric(libsize),
  Tissue = na_chr(pd$Tissue),
  Best_Confirmed_Overall_Response = na_chr(pd[["Best Confirmed Overall Response"]]),
  binaryResponse = na_chr(pd$binaryResponse),
  IC_Level = na_chr(pd[["IC Level"]]),
  TC_Level = na_chr(pd[["TC Level"]]),
  Enrollment_IC = na_chr(pd[["Enrollment IC"]]),
  Immune_phenotype = na_chr(pd[["Immune phenotype"]]),
  Sex = na_chr(pd$Sex),
  Received_platinum = na_chr(pd[["Received platinum"]]),
  os_months = as.numeric(pd$os),
  censOS = as.integer(as.character(pd$censOS)),
  TCGA_Subtype = na_chr(pd[["TCGA Subtype"]]),
  Lund2 = na_chr(pd$Lund2),
  stringsAsFactors = FALSE
)

for (sym in names(wanted)) {
  out[[paste0(sym, "_entrez")]] <- g[[sym]]$entrez
  out[[paste0(sym, "_length")]] <- g[[sym]]$length
  out[[paste0(sym, "_raw")]] <- g[[sym]]$raw
  out[[paste0(sym, "_TPM")]] <- g[[sym]]$tpm
  out[[paste0(sym, "_log2TPM1")]] <- g[[sym]]$log2tpm1
  out[[paste0(sym, "_log2SF")]] <- g[[sym]]$log2sf
}

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
write.table(out, out_path, sep = "\t", quote = FALSE, row.names = FALSE)

cat("Wrote", nrow(out), "rows to", out_path, "\n")
cat("Unique ANONPT_ID:", length(unique(out$ANONPT_ID)), "\n")
for (sym in names(wanted)) {
  cat(
    sprintf(
      "%s entrez=%s length=%s TPM range=%.4g-%.4g\n",
      sym, g[[sym]]$entrez, g[[sym]]$length,
      min(g[[sym]]$tpm), max(g[[sym]]$tpm)
    )
  )
}
