#!/usr/bin/env python3
"""Download the public signature resources the pipeline needs.

Nothing signature-sized is vendored in the repository: every file here is
fetched from its canonical URL and written under methods/bulk_immune/resources/
(or --dest). Re-run this script to refresh versions; the log records what
arrived.

LM22 is *not* fetched. Download it yourself from https://cibersortx.stanford.edu
after accepting the Stanford academic licence and pass --lm22 PATH, or drop
the file at resources/LM22.txt. Without it the pipeline uses the GPL-licensed
quanTIseq TIL10 matrix as a smoke-test substitute and says so in the output.
"""

from __future__ import annotations

import argparse
import io
import sys
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from bulkimmune.geo import download  # noqa: E402
from bulkimmune.signatures import write_gmt  # noqa: E402

URLS = {
    "mcp_counter_genes.txt": (
        "https://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt"
    ),
    "tip_signature_annotation.txt": (
        "http://biocc.hrbmu.edu.cn/TIP/download/signature%20annotation.txt"
    ),
    "TIL10_signature.txt": (
        "https://raw.githubusercontent.com/icbi-lab/quanTIseq/master/quantiseq/deconvolution/TIL10_signature.txt"
    ),
    "TIL10_mRNA_scaling.txt": (
        "https://raw.githubusercontent.com/icbi-lab/quanTIseq/master/quantiseq/deconvolution/TIL10_mRNA_scaling.txt"
    ),
    "TIL10_rmgenes.txt": (
        "https://raw.githubusercontent.com/icbi-lab/quanTIseq/master/quantiseq/deconvolution/TIL10_rmgenes.txt"
    ),
    "h.all.v2025.1.Hs.symbols.gmt": (
        "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Hs/h.all.v2025.1.Hs.symbols.gmt"
    ),
    "hgnc_complete_set.txt": (
        "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
    ),
}

ESTIMATE_TAR = "https://r-forge.r-project.org/src/contrib/estimate_1.0.13.tar.gz"
XCELL_RDA = "https://github.com/dviraran/xCell/raw/master/data/xCell.data.rda"


def _extract_estimate(tar_path: Path, dest: Path) -> None:
    """Pull SI_geneset + common_genes out of the GPL-2 estimate tarball."""
    try:
        import pyreadr
    except ImportError:
        pyreadr = None

    with tarfile.open(tar_path, "r:gz") as tf:
        tf.extractall(tar_path.parent)

    pkg = tar_path.parent / "estimate"
    si_rdata = pkg / "data" / "SI_geneset.RData"
    common_rdata = pkg / "data" / "common_genes.RData"
    if pyreadr is None:
        # Fall back to R if present
        import subprocess

        rscript = f"""
        dir.create("{dest}", showWarnings=FALSE, recursive=TRUE)
        load("{si_rdata}")
        # SI_geneset: 2 x 142, first column is 'estimate'
        sets <- list()
        for (i in seq_len(nrow(SI_geneset))) {{
          genes <- as.character(SI_geneset[i, -1])
          genes <- genes[genes != "" & !is.na(genes)]
          sets[[rownames(SI_geneset)[i]]] <- genes
        }}
        lines <- vapply(names(sets), function(nm) paste(c(nm, "estimate", sets[[nm]]), collapse="\\t"), "")
        writeLines(lines, "{dest / 'estimate_signatures.gmt'}")
        load("{common_rdata}")
        write.table(common_genes, "{dest / 'estimate_common_genes.tsv'}", sep="\\t", quote=FALSE, row.names=FALSE)
        """
        subprocess.run(["Rscript", "-e", rscript], check=True)
        return

    si = pyreadr.read_r(str(si_rdata))["SI_geneset"]
    sets = {}
    for name, row in si.iterrows():
        genes = [g for g in row.tolist()[1:] if isinstance(g, str) and g]
        sets[str(name)] = genes
    write_gmt(sets, dest / "estimate_signatures.gmt", description="estimate")
    common = pyreadr.read_r(str(common_rdata))["common_genes"]
    common.to_csv(dest / "estimate_common_genes.tsv", sep="\t", index=False)


def _extract_xcell(rda: Path, dest: Path) -> None:
    import subprocess
    import tempfile

    script = rda.parent / "_extract_xcell.R"
    script.write_text(
        """
args <- commandArgs(trailingOnly=TRUE)
rda <- args[1]; outdir <- args[2]
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
env <- new.env(); load(rda, envir=env)
obj <- get(ls(env)[1], envir=env)
sigs <- unclass(obj$signatures)
lines <- vapply(seq_along(sigs), function(i) {
  gs <- sigs[[i]]
  paste(c(attr(gs, "setName"), "na", attr(gs, "geneIds")), collapse="\\t")
}, character(1))
writeLines(lines, file.path(outdir, "xCell_signatures.gmt"))
writeLines(as.character(obj$genes), file.path(outdir, "xCell_genes.txt"))
for (spname in c("spill", "spill.array")) {
  sp <- obj[[spname]]
  tag <- ifelse(spname == "spill", "rnaseq", "array")
  write.table(sp$K, file.path(outdir, paste0("xCell_spillK_", tag, ".tsv")),
              sep="\\t", quote=FALSE, col.names=NA)
  write.table(sp$fv, file.path(outdir, paste0("xCell_fv_", tag, ".tsv")),
              sep="\\t", quote=FALSE, col.names=NA)
}
cat("xCell extracted\\n")
"""
    )
    subprocess.run(["Rscript", str(script), str(rda), str(dest)], check=True)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--dest", type=Path, default=ROOT / "resources")
    p.add_argument("--lm22", type=Path, default=None, help="Optional local LM22.txt to copy in")
    p.add_argument("--skip-xcell", action="store_true")
    args = p.parse_args()
    args.dest.mkdir(parents=True, exist_ok=True)

    for name, url in URLS.items():
        dest = args.dest / name
        print(f"[fetch] {name}")
        try:
            download(url, dest)
        except Exception as exc:
            print(f"  FAILED {url}: {exc}")

    print("[fetch] ESTIMATE 1.0.13 (R-Forge, GPL-2)")
    tar = args.dest / "_estimate_1.0.13.tar.gz"
    try:
        download(ESTIMATE_TAR, tar)
        _extract_estimate(tar, args.dest)
    except Exception as exc:
        print(f"  FAILED estimate: {exc}")

    if not args.skip_xcell:
        print("[fetch] xCell.data.rda (requires R to unserialise S4 signatures)")
        rda = args.dest / "_xCell.data.rda"
        try:
            download(XCELL_RDA, rda)
            _extract_xcell(rda, args.dest)
        except Exception as exc:
            print(f"  FAILED xCell extract (is Rscript on PATH?): {exc}")

    if args.lm22 is not None:
        dest = args.dest / "LM22.txt"
        dest.write_bytes(args.lm22.read_bytes())
        print(f"[fetch] copied LM22 -> {dest}")
    elif not (args.dest / "LM22.txt").exists():
        print(
            "[fetch] LM22.txt not present. CIBERSORTx/LM22 will be skipped; "
            "quanTIseq TIL10 will be used as the licence-free substitute."
        )

    print("done. resources in", args.dest)


if __name__ == "__main__":
    main()
