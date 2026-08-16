#!/usr/bin/env Rscript
# Inspect the IMvigor210 cds.RData object structure without assuming Biobase/DESeq
# class definitions are loadable.
suppressWarnings({
  e <- new.env()
  loaded <- load("/workspace/data/raw/IMvigor210/cds.RData", envir = e)
})
cat("objects loaded:", paste(loaded, collapse = ", "), "\n")
obj <- get(loaded[1], envir = e)
cat("class:", paste(class(obj), collapse = "/"), "\n")
cat("isS4:", isS4(obj), "\n")
cat("slots/attrs:", paste(names(attributes(obj)), collapse = ", "), "\n\n")

ad <- attr(obj, "assayData")
cat("assayData class:", paste(class(ad), collapse = "/"), "\n")
if (is.environment(ad)) {
  cat("assayData contents:", paste(ls(ad), collapse = ", "), "\n")
  for (nm in ls(ad)) {
    m <- get(nm, envir = ad)
    cat("  ", nm, ":", paste(class(m), collapse = "/"), "dim=", paste(dim(m), collapse = "x"), "\n")
  }
} else if (is.list(ad)) {
  cat("assayData names:", paste(names(ad), collapse = ", "), "\n")
}

pd <- attr(obj, "phenoData")
pdd <- attr(pd, "data")
cat("\nphenoData dim:", paste(dim(pdd), collapse = "x"), "\n")
cat("phenoData columns:\n")
print(colnames(pdd))

fd <- attr(obj, "featureData")
fdd <- attr(fd, "data")
cat("\nfeatureData dim:", paste(dim(fdd), collapse = "x"), "\n")
cat("featureData columns:\n")
print(colnames(fdd))
cat("\nfeatureData head:\n")
print(utils::head(fdd, 3))
