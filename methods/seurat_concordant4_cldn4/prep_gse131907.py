#!/usr/bin/env python3
"""I/O only: stream GSE131907 UMI text into a capped 10x atlas + malignant UMI sums.

Primary statistics and Harmony integration are in analyze.R (Seurat).
"""
from __future__ import annotations

import gzip
import math
import random
from collections import defaultdict
from pathlib import Path

ANN = Path("/tmp/geo_seurat/GSE131907_Lung_Cancer_cell_annotation.txt.gz")
MAT = Path("/tmp/geo_seurat/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz")
OUT = Path("/tmp/geo_seurat/gse131907_extract")
CAP = 350
SEED = 1
ELIG_ORIGIN = {"tLung", "tL/B", "mLN", "PE", "mBrain"}
MIN_MAL = 20


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    cells = []
    with gzip.open(ANN, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            rec = {
                "index": p[idx["Index"]],
                "barcode": p[idx["Barcode"]],
                "sample": p[idx["Sample"]],
                "origin": p[idx["Sample_Origin"]],
                "cell_type": p[idx["Cell_type"]],
                "cell_subtype": p[idx["Cell_subtype"]],
            }
            rec["author_malignant"] = rec["cell_subtype"] == "Malignant cells"
            rec["author_tnk"] = rec["cell_type"] in {"T lymphocytes", "NK cells"}
            rec["author_nk"] = rec["cell_type"] == "NK cells"
            rec["author_t"] = rec["cell_type"] == "T lymphocytes"
            cells.append(rec)

    by_sample: dict[str, list[dict]] = defaultdict(list)
    for rec in cells:
        by_sample[rec["sample"]].append(rec)

    units = []
    atlas_barcodes = set()
    malig_barcodes = set()
    rng = random.Random(SEED)
    for sample, recs in sorted(by_sample.items()):
        origin = recs[0]["origin"]
        n = len(recs)
        n_mal = sum(r["author_malignant"] for r in recs)
        n_tnk = sum(r["author_tnk"] for r in recs)
        eligible = origin in ELIG_ORIGIN and n_mal >= MIN_MAL
        units.append(
            {
                "dataset": "GSE131907",
                "unit_id": sample,
                "unit_type": "sample",
                "origin": origin,
                "n_cells": n,
                "n_malignant": n_mal,
                "n_tnk": n_tnk,
                "frac_tnk": n_tnk / n if n else 0.0,
                "eligible": eligible,
                "malig_def": "author_malig",
            }
        )
        if not eligible:
            continue
        for r in recs:
            if r["author_malignant"]:
                malig_barcodes.add(r["index"])
        chosen = recs if len(recs) <= CAP else rng.sample(recs, CAP)
        for r in chosen:
            atlas_barcodes.add(r["index"])

    with open(OUT / "units.tsv", "w") as f:
        cols = [
            "dataset",
            "unit_id",
            "unit_type",
            "origin",
            "n_cells",
            "n_malignant",
            "n_tnk",
            "frac_tnk",
            "eligible",
            "malig_def",
        ]
        f.write("\t".join(cols) + "\n")
        for u in units:
            f.write("\t".join(str(u[c]) for c in cols) + "\n")

    keep = atlas_barcodes | malig_barcodes
    print(
        f"eligible units={sum(u['eligible'] for u in units)} "
        f"atlas={len(atlas_barcodes)} malig={len(malig_barcodes)} keep={len(keep)}",
        flush=True,
    )

    # Stream genes x cells. Header barcodes are the Index column (barcode_sample).
    with gzip.open(MAT, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        barcodes = header[1:]
        col_of = {b: i for i, b in enumerate(barcodes)}
        keep_cols = sorted(col_of[b] for b in keep if b in col_of)
        atlas_cols = [c for c in keep_cols if barcodes[c] in atlas_barcodes]
        malig_cols = [c for c in keep_cols if barcodes[c] in malig_barcodes]
        missing = len(keep) - len(keep_cols)
        print(
            f"matrix cells={len(barcodes)} keep_cols={len(keep_cols)} "
            f"atlas_cols={len(atlas_cols)} malig_cols={len(malig_cols)} missing={missing}",
            flush=True,
        )

        atlas_bc = [barcodes[c] for c in atlas_cols]
        malig_bc = [barcodes[c] for c in malig_cols]
        sample_of = {}
        for rec in cells:
            sample_of[rec["index"]] = rec["sample"]
        elig_samples = [u["unit_id"] for u in units if u["eligible"]]
        malig_by_sample = {s: [i for i, b in enumerate(malig_bc) if sample_of.get(b) == s] for s in elig_samples}

        genes = []
        # triplet lists for atlas MTX
        rows_i = []
        cols_j = []
        data_x = []
        malig_sums = {s: [] for s in elig_samples}
        cldn4_pos = {s: [0, 0] for s in elig_samples}  # pos, n
        cldn4_sum_log1p = {s: 0.0 for s in elig_samples}

        gene_i = 0
        for line in f:
            parts = line.rstrip("\n").split("\t")
            gene = parts[0]
            genes.append(gene)
            vals = parts[1:]
            for j, c in enumerate(atlas_cols):
                v = vals[c]
                if v not in {"0", "0.0", ""}:
                    fv = float(v)
                    if fv != 0:
                        rows_i.append(gene_i + 1)
                        cols_j.append(j + 1)
                        data_x.append(int(fv) if fv.is_integer() else fv)
            is_cldn4 = gene.upper() == "CLDN4"
            for s, idxs in malig_by_sample.items():
                total = 0.0
                npos = 0
                sm = 0.0
                for i in idxs:
                    v = vals[malig_cols[i]]
                    fv = 0.0 if v in {"0", "0.0", ""} else float(v)
                    total += fv
                    if is_cldn4:
                        sm += math.log1p(fv)
                        if fv > 0:
                            npos += 1
                malig_sums[s].append(total)
                if is_cldn4:
                    cldn4_pos[s] = [npos, len(idxs)]
                    cldn4_sum_log1p[s] = sm
            gene_i += 1
            if gene_i % 2000 == 0:
                print(f"  streamed {gene_i} genes, nnz={len(data_x)}", flush=True)

    n_genes = len(genes)
    n_atlas = len(atlas_bc)
    print(f"done stream genes={n_genes} atlas_nnz={len(data_x)}", flush=True)

    # write 10x-like MTX for atlas
    mtx_dir = OUT / "atlas_10x"
    mtx_dir.mkdir(exist_ok=True)
    with gzip.open(mtx_dir / "barcodes.tsv.gz", "wt") as f:
        for b in atlas_bc:
            f.write(b + "\n")
    with gzip.open(mtx_dir / "features.tsv.gz", "wt") as f:
        for g in genes:
            f.write(f"{g}\t{g}\tGene Expression\n")
    with gzip.open(mtx_dir / "matrix.mtx.gz", "wt") as f:
        f.write("%%MatrixMarket matrix coordinate real general\n")
        f.write(f"{n_genes} {n_atlas} {len(data_x)}\n")
        for r, c, x in zip(rows_i, cols_j, data_x):
            f.write(f"{r} {c} {x}\n")

    rec_by_index = {r["index"]: r for r in cells}
    with open(OUT / "atlas_meta.tsv", "w") as f:
        f.write(
            "barcode\tunit_id\torigin\tcell_type\tcell_subtype\t"
            "author_malignant\tauthor_tnk\tauthor_t\tauthor_nk\n"
        )
        for b in atlas_bc:
            r = rec_by_index[b]
            f.write(
                f"{b}\t{r['sample']}\t{r['origin']}\t{r['cell_type']}\t{r['cell_subtype']}\t"
                f"{r['author_malignant']}\t{r['author_tnk']}\t{r['author_t']}\t{r['author_nk']}\n"
            )

    with gzip.open(OUT / "malig_pseudobulk.tsv.gz", "wt") as f:
        f.write("gene\t" + "\t".join(elig_samples) + "\n")
        for i, g in enumerate(genes):
            f.write(g + "\t" + "\t".join(str(malig_sums[s][i]) for s in elig_samples) + "\n")

    # attach CLDN4 to units
    unit_path = OUT / "units.tsv"
    with open(unit_path) as f:
        rows = [line.rstrip("\n").split("\t") for line in f]
    header = rows[0] + ["mal_CLDN4_pct", "mal_CLDN4_mean"]
    out_rows = [header]
    for row in rows[1:]:
        rec = dict(zip(rows[0], row))
        s = rec["unit_id"]
        if rec["eligible"] == "True" and s in cldn4_pos and cldn4_pos[s][1] > 0:
            npos, n = cldn4_pos[s]
            pct = 100.0 * npos / n
            mean = cldn4_sum_log1p[s] / n
        else:
            pct = ""
            mean = ""
        out_rows.append(row + [str(pct), str(mean)])
    with open(unit_path, "w") as f:
        for row in out_rows:
            f.write("\t".join(row) + "\n")
    print("wrote", OUT, flush=True)


if __name__ == "__main__":
    main()
