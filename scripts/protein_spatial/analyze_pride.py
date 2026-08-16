#!/usr/bin/env python3
"""Confirm TACSTD2/CLDN4 presence in PXD042091 (human ICI plasma) and PXD059688 (mouse LLC)."""

from __future__ import annotations

import gzip
import json
import re
from pathlib import Path

import pandas as pd

from config import (
    CLDN4_MOUSE,
    CLDN4_MOUSE_UNIPROT,
    CLDN4_SYMBOL,
    CLDN4_UNIPROT,
    DATA,
    RESULTS,
    TACSTD2_MOUSE,
    TACSTD2_MOUSE_UNIPROT,
    TACSTD2_SYMBOL,
    TACSTD2_UNIPROT,
)

HUMAN_PAT = re.compile(
    r"TACSTD2|TACD2_HUMAN|TROP-?2|P09758|CLDN4|CLD4_HUMAN|O14493",
    re.I,
)
MOUSE_PAT = re.compile(
    r"Tacstd2|TACSTD2|TROP-?2|Q8BGV3|Cldn4|CLDN4|O35114",
    re.I,
)


def scan_text(path: Path, pattern: re.Pattern, max_hits: int = 8) -> dict:
    n = 0
    examples = []
    with path.open(errors="ignore") as f:
        for i, line in enumerate(f, 1):
            if pattern.search(line):
                n += 1
                if len(examples) < max_hits:
                    examples.append({"line": i, "text": line.strip()[:240]})
    return {"n_matching_lines": n, "examples": examples}


def unique_library_proteins(path: Path) -> tuple[set[str], set[str]]:
    names = set()
    hits = set()
    with path.open(errors="ignore") as f:
        header = f.readline()
        for line in f:
            parts = line.split("\t")
            if len(parts) > 3:
                names.add(parts[3])
                if HUMAN_PAT.search(parts[3]) or (len(parts) > 13 and HUMAN_PAT.search(parts[13])):
                    hits.add(parts[3])
    return names, hits


def quantified_ids(path: Path, col: str) -> pd.Series:
    df = pd.read_csv(path, sep="\t")
    return df[col].astype(str)


def mztab_proteins(path: Path) -> list[dict]:
    rows = []
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", errors="ignore") as f:
        for line in f:
            if line.startswith("PRT\t"):
                parts = line.rstrip("\n").split("\t")
                rows.append(
                    {
                        "accession": parts[1] if len(parts) > 1 else "",
                        "description": parts[2] if len(parts) > 2 else "",
                    }
                )
    return rows


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    recs = []

    # PXD042091
    lib = DATA / "pxd042091/NSClibrary.txt"
    q2 = DATA / "pxd042091/2_dat_cs_all.txt"
    q6 = DATA / "pxd042091/6_data_tf_response.txt"
    lib_names, lib_hits = unique_library_proteins(lib)
    acc2 = quantified_ids(q2, "accessions")
    acc6 = quantified_ids(q6, "protein")

    def in_series(s: pd.Series, *needles: str) -> bool:
        up = s.str.upper()
        return any(up.str.contains(n.upper(), regex=False).any() for n in needles)

    recs.append(
        {
            "dataset": "PXD042091",
            "file": lib.name,
            "role": "spectral_library",
            "n_unique_protein_name": len(lib_names),
            "TACSTD2_present": any("TACD2_HUMAN" in x or "GN=TACSTD2" in x for x in lib_hits)
            or any("TACD2_HUMAN" in x or "GN=TACSTD2" in x for x in lib_names),
            "CLDN4_present": any("CLD4_HUMAN" in x or "GN=CLDN4" in x for x in lib_hits)
            or any("CLD4_HUMAN" in x or "GN=CLDN4" in x for x in lib_names),
            "note": "Library contains TACD2_HUMAN and CLD4_HUMAN peptide entries; this is not quantification.",
        }
    )
    recs.append(
        {
            "dataset": "PXD042091",
            "file": q6.name,
            "role": "SWATH_quant_response_table",
            "n_proteins": int(acc6.nunique()),
            "TACSTD2_present": bool(in_series(acc6, TACSTD2_UNIPROT, "TACSTD2", "TACD2")),
            "CLDN4_present": bool(in_series(acc6, CLDN4_UNIPROT, "CLDN4", "CLD4")),
            "note": "protein column is UniProt_Gene (e.g. P0CG40_SP9). Neither P09758 nor O14493 is present.",
        }
    )
    recs.append(
        {
            "dataset": "PXD042091",
            "file": q2.name,
            "role": "SWATH_quant_all_samples",
            "n_proteins": int(acc2.nunique()),
            "TACSTD2_present": bool(in_series(acc2, TACSTD2_UNIPROT, "TACSTD2", "TACD2")),
            "CLDN4_present": bool(in_series(acc2, CLDN4_UNIPROT, "CLDN4", "CLD4")),
            "note": "Same accessions scheme as 6_data_tf_response.txt; TACSTD2/CLDN4 absent.",
        }
    )

    # PXD059688
    mz = DATA / "pxd059688/20230904_Tumor_AA.mzTab.gz"
    prt = mztab_proteins(mz)
    genes = []
    for r in prt:
        m = re.search(r"GN=([A-Za-z0-9]+)", r["description"])
        genes.append(m.group(1) if m else "")
    recs.append(
        {
            "dataset": "PXD059688",
            "file": mz.name,
            "role": "processed_mzTab_quantification",
            "n_proteins": len(prt),
            "TACSTD2_present": any(
                TACSTD2_MOUSE.lower() == g.lower() or r["accession"] == TACSTD2_MOUSE_UNIPROT
                for r, g in zip(prt, genes)
            ),
            "CLDN4_present": any(
                CLDN4_MOUSE.lower() == g.lower() or r["accession"] == CLDN4_MOUSE_UNIPROT
                for r, g in zip(prt, genes)
            ),
            "note": (
                "Only open processed table <2GB. 93 PRT rows are almost entirely "
                "methyltransferases (export of 20230904_Tumor_AA.pdResult). "
                "Tacstd2/Cldn4 are absent from this table. Full proteome is in "
                "skipped 8.5GB mgf / 37.9GB msf — presence in the full LLC proteome "
                "cannot be confirmed from allowed files."
            ),
        }
    )
    pd.DataFrame(prt).assign(gene=genes).to_csv(
        RESULTS / "pxd059688_mztab_proteins.tsv", sep="\t", index=False
    )

    out = pd.DataFrame(recs)
    out.to_csv(RESULTS / "pride_gene_presence.tsv", sep="\t", index=False)

    detail = {
        "PXD042091_library_target_hits": sorted(lib_hits),
        "PXD042091_quant_TACSTD2": False,
        "PXD042091_quant_CLDN4": False,
        "PXD059688_example_genes": sorted({g for g in genes if g})[:40],
        "skipped_raw": [
            "PXD042091 wiff/wiff.scan",
            "PXD059688 mgf 8492466377 bytes",
            "PXD059688 msf 37928718336 bytes",
        ],
    }
    (RESULTS / "pride_gene_presence_detail.json").write_text(json.dumps(detail, indent=2))
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
