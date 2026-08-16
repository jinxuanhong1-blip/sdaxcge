#!/usr/bin/env python3
"""End-to-end demo: public LUAD/LSCC tables -> analyses -> methods/proteomics/demo/."""

from __future__ import annotations

import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from fetch_public import fetch_linkedomics, fetch_star_file, gdc_query_lung_star  # noqa: E402
from lib_io import (  # noqa: E402
    ENSEMBL_TO_SYMBOL,
    extract_star_targets,
    normalize_sample_id,
    read_expression_matrix,
    read_gdc_star_gene_counts,
)


def run(cmd: list[str]) -> None:
    print("+", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def write_tumor_nat_group(tumor: Path, normal: Path, dest: Path) -> Path:
    t = read_expression_matrix(tumor)
    n = read_expression_matrix(normal)
    rows = [{"sample_id": normalize_sample_id(c), "group": "Tumor"} for c in t.columns]
    rows += [{"sample_id": normalize_sample_id(c), "group": "NAT"} for c in n.columns]
    dest.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).drop_duplicates("sample_id").to_csv(dest, sep="\t", index=False)
    return dest


def merge_tumor_nat(tumor: Path, normal: Path, dest: Path) -> Path:
    t = read_expression_matrix(tumor)
    n = read_expression_matrix(normal)
    t.columns = [normalize_sample_id(c) for c in t.columns]
    n.columns = [normalize_sample_id(c) for c in n.columns]
    # overlapping genes only
    genes = t.index.intersection(n.index)
    merged = pd.concat([t.loc[genes], n.loc[genes]], axis=1)
    # drop duplicated sample IDs keeping tumor
    merged = merged.loc[:, ~merged.columns.duplicated()]
    dest.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(dest, sep="\t")
    return dest


def extract_star_bundle(cache: Path, dest: Path, n_each: int = 2) -> pd.DataFrame:
    frames = []
    hits = []
    for histo in ("LUAD", "LSCC"):
        for h in gdc_query_lung_star(histo, size=n_each):
            p = fetch_star_file(h, cache)
            star = read_gdc_star_gene_counts(p)
            df = extract_star_targets(star, h["submitter_id"])
            df["histology"] = histo
            df["file_id"] = h["file_id"]
            df["s3_uri"] = h["s3_uri"]
            frames.append(df)
            hits.append(h)
    out = pd.concat(frames, ignore_index=True)
    dest.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(dest, sep="\t", index=False)
    (dest.parent / "gdc_star_hits.json").write_text(json.dumps(hits, indent=2), encoding="utf-8")
    return out


def main() -> int:
    cache = ROOT / "demo" / "cache"
    tables = ROOT / "demo" / "tables"
    figs = ROOT / "demo" / "figures"
    tables.mkdir(parents=True, exist_ok=True)
    figs.mkdir(parents=True, exist_ok=True)

    paths = fetch_linkedomics(cache)
    star_tbl = extract_star_bundle(cache, tables / "gdc_s3_star_targets.tsv", n_each=2)

    # missingness on protein + one STAR file
    star_example = next(cache.joinpath("star").glob("*.star.tsv"))
    run(
        [
            sys.executable,
            str(HERE / "log2_missingness.py"),
            str(paths["LUAD_protein_tumor"]),
            str(paths["LSCC_protein_tumor"]),
            str(paths["LUAD_rna_tumor"]),
            str(paths["LSCC_rna_tumor"]),
            str(star_example),
            "--out-json",
            str(tables / "missingness_log2.json"),
            "--out-tsv",
            str(tables / "target_presence.tsv"),
        ]
    )

    run(
        [
            sys.executable,
            str(HERE / "protein_rna_concordance.py"),
            "--protein",
            "LUAD",
            str(paths["LUAD_protein_tumor"]),
            "--protein",
            "LSCC",
            str(paths["LSCC_protein_tumor"]),
            "--rna",
            "LUAD",
            str(paths["LUAD_rna_tumor"]),
            "--rna",
            "LSCC",
            str(paths["LSCC_rna_tumor"]),
            "--genes",
            "TACSTD2",
            "CLDN4",
            "--out-stats",
            str(tables / "protein_rna_concordance.tsv"),
            "--out-pairs",
            str(tables / "protein_rna_pairs.tsv"),
        ]
    )

    merged = merge_tumor_nat(paths["LUAD_protein_tumor"], paths["LUAD_protein_normal"], cache / "LUAD_protein_tumor_nat.tsv")
    groups = write_tumor_nat_group(paths["LUAD_protein_tumor"], paths["LUAD_protein_normal"], tables / "LUAD_tumor_nat_groups.tsv")
    run(
        [
            sys.executable,
            str(HERE / "limma_trend.py"),
            "--matrix",
            str(merged),
            "--group-tsv",
            str(groups),
            "--method",
            "limma",
            "--contrast",
            "Tumor",
            "NAT",
            "--out",
            str(tables / "limma_LUAD_protein_tumor_vs_nat.tsv"),
        ]
    )

    run(
        [
            sys.executable,
            str(HERE / "join_deconvolution.py"),
            "--protein",
            str(paths["LSCC_protein_tumor"]),
            "--rna",
            str(paths["LSCC_rna_tumor"]),
            "--phenotype",
            str(paths["LSCC_molecular"]),
            "--cohort",
            "LSCC",
            "--genes",
            "TACSTD2",
            "CLDN4",
            "--out",
            str(tables / "LSCC_protein_rna_immune_join.tsv"),
        ]
    )

    run(
        [
            sys.executable,
            str(HERE / "make_demo_figures.py"),
            "--protein-luad",
            str(paths["LUAD_protein_tumor"]),
            "--protein-lscc",
            str(paths["LSCC_protein_tumor"]),
            "--rna-luad",
            str(paths["LUAD_rna_tumor"]),
            "--rna-lscc",
            str(paths["LSCC_rna_tumor"]),
            "--pairs",
            str(tables / "protein_rna_pairs.tsv"),
            "--concord-stats",
            str(tables / "protein_rna_concordance.tsv"),
            "--limma",
            str(tables / "limma_LUAD_protein_tumor_vs_nat.tsv"),
            "--joined",
            str(tables / "LSCC_protein_rna_immune_join.tsv"),
            "--outdir",
            str(figs),
        ]
    )

    # compact target extract (safe to commit; not the 60 MB matrices)
    compact_rows = []
    for cohort, pkey, rkey in (
        ("LUAD", "LUAD_protein_tumor", "LUAD_rna_tumor"),
        ("LSCC", "LSCC_protein_tumor", "LSCC_rna_tumor"),
    ):
        prot = read_expression_matrix(paths[pkey])
        rna = read_expression_matrix(paths[rkey])
        for gene in ("TACSTD2", "CLDN4"):
            for layer, mat in (("protein", prot), ("rna", rna)):
                if gene not in mat.index:
                    compact_rows.append({"cohort": cohort, "layer": layer, "gene": gene, "status": "absent"})
                    continue
                s = mat.loc[gene]
                compact_rows.append(
                    {
                        "cohort": cohort,
                        "layer": layer,
                        "gene": gene,
                        "status": "present",
                        "n": int(s.notna().sum()),
                        "na_frac": float(s.isna().mean()),
                        "median": float(s.median()),
                        "min": float(s.min()),
                        "max": float(s.max()),
                    }
                )
    pd.DataFrame(compact_rows).to_csv(tables / "target_gene_summary.tsv", sep="\t", index=False)

    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "ensembl": ENSEMBL_TO_SYMBOL,
        "star_targets_n": int(star_tbl.shape[0]),
        "figures": sorted(p.name for p in figs.glob("*.png")),
        "note": (
            "Full CCT matrices stay in demo/cache (gitignored). "
            "CLDN4 protein is absent from NArm TMT tables; RNA and GDC STAR are not."
        ),
    }
    (ROOT / "demo" / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print("demo complete", json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
