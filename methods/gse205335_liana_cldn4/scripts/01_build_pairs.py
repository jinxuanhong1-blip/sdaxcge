#!/usr/bin/env python3
"""Build a documented CellPhoneDB v5 ligand–receptor pair table.

Parses the official `interactors` column. Hyphenated gene symbols (HLA-A)
are split only when both sides match the CellPhoneDB gene list.
Complexes use '+' subunits. This is the pair list used by the lightweight
score. It is not CellChat output.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import yaml

HERE = Path(__file__).resolve().parents[1]


def parse_interactors(text: str, genes: set[str]) -> tuple[str, str] | None:
    if not isinstance(text, str) or "-" not in text:
        return None
    for i, ch in enumerate(text):
        if ch != "-":
            continue
        left, right = text[:i], text[i + 1 :]
        if not left or not right:
            continue
        left_ok = all(part in genes for part in left.split("+"))
        right_ok = all(part in genes for part in right.split("+"))
        if left_ok and right_ok:
            return left, right
    return None


def pathway_of(ligand: str, receptor: str, cfg: dict) -> str:
    lig_units = set(ligand.split("+"))
    rec_units = set(receptor.split("+"))
    for name, block in cfg["pathways"].items():
        if lig_units & set(block["ligands"]) and rec_units & set(block["receptors"]):
            return name
    for name, block in cfg["pathways"].items():
        if lig_units & set(block["ligands"]) or rec_units & set(block["receptors"]):
            return name + "_partial"
    return "other"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--datadir", type=Path, default=Path("/tmp/gse205335"))
    parser.add_argument("--outdir", type=Path, default=HERE / "resources")
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    cfg = yaml.safe_load((HERE / "config" / "gene_sets.yaml").read_text())
    genes = pd.read_csv(args.datadir / "cellphonedb" / "gene_input.csv")
    valid = set(genes["gene_name"].dropna().astype(str)) | set(
        genes["hgnc_symbol"].dropna().astype(str)
    )
    valid.update(
        {
            "HLA-A",
            "HLA-B",
            "HLA-C",
            "HLA-E",
            "HLA-F",
            "HLA-G",
            "HLA-DRA",
            "HLA-DRB1",
            "B2M",
            "IFNG",
            "IFNGR1",
            "IFNGR2",
        }
    )

    inter = pd.read_csv(args.datadir / "cellphonedb" / "interaction_input.csv")
    rows = []
    n_fail = 0
    for rec in inter.itertuples(index=False):
        parsed = parse_interactors(getattr(rec, "interactors", None), valid)
        if parsed is None:
            n_fail += 1
            continue
        lig, recp = parsed
        rows.append(
            {
                "ligand": lig,
                "receptor": recp,
                "classification": getattr(rec, "classification", ""),
                "directionality": getattr(rec, "directionality", ""),
                "is_ppi": getattr(rec, "is_ppi", True),
                "source": getattr(rec, "source", ""),
                "version": getattr(rec, "version", ""),
                "pair_origin": "cellphonedb_v5_interactors",
            }
        )
    pairs = pd.DataFrame(rows).drop_duplicates(subset=["ligand", "receptor"])

    overlay = []
    for item in cfg["curated_overlay"]:
        overlay.append(
            {
                "ligand": item["ligand"],
                "receptor": item["receptor"],
                "classification": "curated_overlay",
                "directionality": "Ligand-Receptor",
                "is_ppi": True,
                "source": "config/gene_sets.yaml",
                "version": "overlay",
                "pair_origin": "curated_overlay",
            }
        )
    pairs = pd.concat([pairs, pd.DataFrame(overlay)], ignore_index=True)
    pairs = pairs.drop_duplicates(subset=["ligand", "receptor"], keep="first")
    pairs["pathway"] = [
        pathway_of(a, b, cfg) for a, b in zip(pairs["ligand"], pairs["receptor"])
    ]
    pairs = pairs.sort_values(["pathway", "ligand", "receptor"]).reset_index(drop=True)

    out = args.outdir / "cellphonedb_v5_lr_pairs.tsv"
    pairs.to_csv(out, sep="\t", index=False)
    genes_needed = sorted(
        {u for s in pairs["ligand"] for u in str(s).split("+")}
        | {u for s in pairs["receptor"] for u in str(s).split("+")}
    )
    (args.outdir / "lr_genes.txt").write_text("\n".join(genes_needed) + "\n")
    print(
        f"pairs={len(pairs)} parse_fail={n_fail} unique_genes={len(genes_needed)} -> {out}",
        flush=True,
    )
    print(pairs["pathway"].value_counts().to_string(), flush=True)


if __name__ == "__main__":
    main()
