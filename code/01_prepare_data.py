#!/usr/bin/env python3
"""Convert the merged ICB compendium TSV into fast-loading numpy arrays.

Input (data/raw/, downloaded from Zenodo record 7459023 "ICB Data TSV Files"):
  merged_expr.tsv      genes x samples, log2 expression (see units check below)
  merged_metadata.tsv  one row per patient/sample
  gene_metadata.tsv    Ensembl gene id -> HUGO symbol (Gencode v40)

Output (data/derived/):
  expr.npy         float32 matrix, genes x samples
  expr_genes.txt   Ensembl ids (row order)
  expr_samples.txt sample ids (column order)
  gene_symbols.tsv ensembl_no_version -> symbol
"""
import csv
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
DER = ROOT / "data" / "derived"
DER.mkdir(parents=True, exist_ok=True)


def main() -> None:
    expr_path = RAW / "merged_expr.tsv"
    genes: list[str] = []
    rows: list[np.ndarray] = []

    import pandas as pd

    with expr_path.open() as fh:
        samples = next(csv.reader(fh, delimiter="\t", quotechar='"'))

    for chunk in pd.read_csv(
        expr_path,
        sep="\t",
        header=0,
        names=["gene_id"] + samples,
        skiprows=1,
        na_values=["NA", "NaN", ""],
        chunksize=2000,
        quotechar='"',
    ):
        genes.extend(chunk["gene_id"].tolist())
        rows.append(chunk.drop(columns=["gene_id"]).to_numpy(dtype=np.float32))
        print(f"  {len(genes)} genes parsed", file=sys.stderr)

    mat = np.vstack(rows)
    print(f"matrix {mat.shape}, samples {len(samples)}", file=sys.stderr)
    assert mat.shape[1] == len(samples), (mat.shape, len(samples))

    np.save(DER / "expr.npy", mat)
    (DER / "expr_genes.txt").write_text("\n".join(genes) + "\n")
    (DER / "expr_samples.txt").write_text("\n".join(samples) + "\n")

    with (RAW / "gene_metadata.tsv").open() as fh, (DER / "gene_symbols.tsv").open("w") as out:
        reader = csv.reader(fh, delimiter="\t", quotechar='"')
        header = next(reader)
        # The published file carries a duplicated header line and the row names
        # are shifted by one column relative to the header.
        idx_id = header.index("gene_id")
        idx_sym = header.index("gene_name")
        out.write("gene_id\tgene_id_no_ver\tsymbol\n")
        seen = 0
        for rec in reader:
            if rec[0] == "seqnames":
                continue
            gid = rec[0]
            sym = rec[idx_sym + 1] if len(rec) > idx_sym + 1 else ""
            gid_nov = gid.split(".")[0]
            out.write(f"{gid}\t{gid_nov}\t{sym}\n")
            seen += 1
    print(f"gene symbol rows: {seen}", file=sys.stderr)


if __name__ == "__main__":
    main()
