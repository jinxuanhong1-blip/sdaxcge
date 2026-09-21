#!/usr/bin/env python3
"""Exploratory inventory: cell painting and DNA-damage morphology vs CLDN4.

Downloads small public metadata (not image archives) and writes:
  results/tables/cldn4_checks.tsv
  results/tables/inventory_cellpainting.tsv

Cell Painting DNA/Hoechst/DAPI features are nuclear morphology. They are not
a gamma-H2AX or comet DNA-damage assay.
"""

from __future__ import annotations

import csv
import gzip
import urllib.request
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "cache"
OUT = ROOT / "results" / "tables"

ORF_ID = "JCP2022_905637"
CRISPR_ID = "JCP2022_801379"

URLS = {
    "orf": "https://raw.githubusercontent.com/jump-cellpainting/datasets/main/metadata/orf.csv.gz",
    "crispr": "https://raw.githubusercontent.com/jump-cellpainting/datasets/main/metadata/crispr.csv.gz",
    "well": "https://raw.githubusercontent.com/jump-cellpainting/datasets/main/metadata/well.csv.gz",
    "jump_target_orf": "https://raw.githubusercontent.com/jump-cellpainting/JUMP-Target/master/JUMP-Target-1_orf_metadata.tsv",
    "geneinfo": "https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/geneinfo_beta.txt",
    "cell_health": "https://raw.githubusercontent.com/broadinstitute/cell-health/master/1.generate-profiles/tables/supplementary_table_1_perturbation_details.tsv",
    "ccle": "https://raw.githubusercontent.com/broadinstitute/2022_PERISCOPE/main/Hit_Calling/inputs/CCLE_expression_HeLa_A549.csv",
    "a549_ngs": "https://raw.githubusercontent.com/broadinstitute/2022_PERISCOPE/main/Extended_Data_2_and_8/inputs/A549_NGS_Counts.csv",
    "a549_whole": "https://raw.githubusercontent.com/broadinstitute/2022_PERISCOPE/main/4_A549_Screen_Summary/outputs/A549_plate_level_median_per_feat_sig_genes_1_FDR_whole_cell_hits.csv",
    "a549_comp": "https://raw.githubusercontent.com/broadinstitute/2022_PERISCOPE/main/4_A549_Screen_Summary/outputs/A549_plate_level_median_per_feat_sig_genes_1_FDR_compartment_specific_hits.csv",
    "hela_dmem_whole": "https://raw.githubusercontent.com/broadinstitute/2022_PERISCOPE/main/2_HeLa_Screens_Summary/outputs/HeLa_DMEM_plate_level_median_per_feat_sig_genes_1_FDR_whole_cell_hits.csv",
    "hela_hplm_whole": "https://raw.githubusercontent.com/broadinstitute/2022_PERISCOPE/main/2_HeLa_Screens_Summary/outputs/HeLa_HPLM_plate_level_median_per_feat_sig_genes_1_FDR_whole_cell_hits.csv",
    "hela_dmem_comp": "https://raw.githubusercontent.com/broadinstitute/2022_PERISCOPE/main/2_HeLa_Screens_Summary/outputs/HeLa_DMEM_plate_level_median_per_feat_sig_genes_1_FDR_compartment_specific_hits.csv",
    "hela_hplm_comp": "https://raw.githubusercontent.com/broadinstitute/2022_PERISCOPE/main/2_HeLa_Screens_Summary/outputs/HeLa_HPLM_plate_level_median_per_feat_sig_genes_1_FDR_compartment_specific_hits.csv",
    "hela_dmem_sig": "https://raw.githubusercontent.com/broadinstitute/2022_PERISCOPE/main/2_HeLa_Screens_Summary/outputs/HeLa_DMEM_significant_features_expressed_genes.csv",
    "hela_hplm_sig": "https://raw.githubusercontent.com/broadinstitute/2022_PERISCOPE/main/2_HeLa_Screens_Summary/outputs/HeLa_HPLM_significant_features_expressed_genes.csv",
    "hela_dmem_zero": "https://raw.githubusercontent.com/broadinstitute/2022_PERISCOPE/main/2_HeLa_Screens_Summary/outputs/HeLa_DMEM_significant_features_0TPM_genes.csv",
    "rohban": "https://raw.githubusercontent.com/carpenterlab/2017_rohban_elife/master/results/master/ORFs_sequence_matching_transcripts_percentage/matching.per.table.csv",
    "luad_genes": "https://raw.githubusercontent.com/broadinstitute/luad-cell-painting/main/inputs/metadata/l1k/gene_funcs.csv",
    "rxrx3": "https://huggingface.co/datasets/recursionpharma/rxrx3-core/resolve/main/metadata_rxrx3_core.csv",
}


def fetch(name: str) -> Path:
    CACHE.mkdir(parents=True, exist_ok=True)
    dest = CACHE / name
    # keep original suffix when the url has one
    url = URLS[name]
    suffix = Path(url.split("?")[0]).suffix
    if suffix and dest.suffix != suffix:
        dest = CACHE / f"{name}{suffix}"
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    print(f"download {name}")
    urllib.request.urlretrieve(url, dest)
    return dest


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", newline="")
    return path.open("rt", newline="")


def read_csv(path: Path, delimiter: str = ","):
    with open_text(path) as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def percentile_leq(values: list[float], value: float) -> float:
    return sum(1 for item in values if item <= value) / len(values)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    checks: list[dict] = []

    orf_rows = read_csv(fetch("orf"))
    orf_hits = [row for row in orf_rows if row.get("Metadata_Symbol") == "CLDN4"]
    checks.append(
        {
            "check": "jump_orf_metadata",
            "n_rows": len(orf_rows),
            "cldn4_rows": len(orf_hits),
            "detail": ""
            if not orf_hits
            else (
                f"{orf_hits[0]['Metadata_JCP2022']}; {orf_hits[0]['Metadata_Transcript']}; "
                f"insert {orf_hits[0]['Metadata_Insert_Length']} nt; "
                f"prot match {orf_hits[0]['Metadata_Prot_Match']}; "
                f"vector {orf_hits[0]['Metadata_Vector']}"
            ),
        }
    )

    crispr_rows = read_csv(fetch("crispr"))
    crispr_hits = [row for row in crispr_rows if row.get("Metadata_Symbol") == "CLDN4"]
    checks.append(
        {
            "check": "jump_crispr_metadata",
            "n_rows": len(crispr_rows),
            "cldn4_rows": len(crispr_hits),
            "detail": "" if not crispr_hits else crispr_hits[0]["Metadata_JCP2022"],
        }
    )

    well_counts: Counter[str] = Counter()
    well_examples: dict[str, list[str]] = {ORF_ID: [], CRISPR_ID: []}
    with gzip.open(fetch("well"), "rt") as handle:
        reader = csv.DictReader(handle)
        n_wells = 0
        for row in reader:
            n_wells += 1
            jcp = row["Metadata_JCP2022"]
            if jcp in well_examples:
                well_counts[jcp] += 1
                well_examples[jcp].append(
                    f"{row['Metadata_Source']}:{row['Metadata_Plate']}:{row['Metadata_Well']}"
                )
    checks.append(
        {
            "check": "jump_wells",
            "n_rows": n_wells,
            "cldn4_rows": well_counts[ORF_ID] + well_counts[CRISPR_ID],
            "detail": (
                f"ORF {ORF_ID} n={well_counts[ORF_ID]} ({'; '.join(well_examples[ORF_ID])}); "
                f"CRISPR {CRISPR_ID} n={well_counts[CRISPR_ID]} ({'; '.join(well_examples[CRISPR_ID])})"
            ),
        }
    )

    target_rows = read_csv(fetch("jump_target_orf"), delimiter="\t")
    target_hits = [row for row in target_rows if row.get("gene") == "CLDN4"]
    checks.append(
        {
            "check": "jump_target1_orf_plate",
            "n_rows": len(target_rows),
            "cldn4_rows": len(target_hits),
            "detail": "384-well target-subset plate, not the genome-scale JUMP ORF library",
        }
    )

    l1000 = read_csv(fetch("geneinfo"), delimiter="\t")
    l1000_hits = [row for row in l1000 if row.get("gene_symbol") == "CLDN4"]
    checks.append(
        {
            "check": "l1000_geneinfo_beta",
            "n_rows": len(l1000),
            "cldn4_rows": len(l1000_hits),
            "detail": ""
            if not l1000_hits
            else (
                f"gene_id {l1000_hits[0]['gene_id']}; "
                f"feature_space={l1000_hits[0]['feature_space']}; "
                f"{l1000_hits[0]['ensembl_id']}"
            ),
        }
    )

    health = read_csv(fetch("cell_health"), delimiter="\t")
    health_genes = sorted({row["gene_name"] for row in health})
    checks.append(
        {
            "check": "cell_health_supp_table1",
            "n_rows": len(health),
            "cldn4_rows": sum(1 for gene in health_genes if gene.upper() == "CLDN4"),
            "detail": f"unique gene_name field n={len(health_genes)}; CLDN4 absent",
        }
    )

    ccle = read_csv(fetch("ccle"))
    ccle_hits = [row for row in ccle if row.get("Gene") == "CLDN4"]
    if len(ccle_hits) != 1:
        raise SystemExit(f"expected one CCLE CLDN4 row, found {len(ccle_hits)}")
    a549_log = float(ccle_hits[0]["A549"])
    hela_log = float(ccle_hits[0]["HeLa"])
    checks.append(
        {
            "check": "periscope_depmap21q4_ccle",
            "n_rows": len(ccle),
            "cldn4_rows": 1,
            "detail": (
                f"A549 log2(TPM+1)={a549_log:.6f} (TPM={2 ** a549_log - 1:.2f}); "
                f"HeLa log2(TPM+1)={hela_log:.6f} (TPM={2 ** hela_log - 1:.2f}). "
                "Baseline cell-line RNA from DepMap Public 21Q4 CCLE_expression.csv, "
                "not RNA measured in the imaging wells."
            ),
        }
    )

    ngs = read_csv(fetch("a549_ngs"))
    ngs_hits = [row for row in ngs if row.get("gene_symbol") == "CLDN4"]
    checks.append(
        {
            "check": "periscope_a549_ngs_barcodes",
            "n_rows": len(ngs),
            "cldn4_rows": len(ngs_hits),
            "detail": "; ".join(
                f"{row['Barcode']} count={row['Count']}" for row in ngs_hits
            ),
        }
    )

    hit_files = {
        "a549_whole": "A549 1% FDR whole-cell",
        "a549_comp": "A549 1% FDR compartment",
        "hela_dmem_whole": "HeLa DMEM 1% FDR whole-cell",
        "hela_hplm_whole": "HeLa HPLM 1% FDR whole-cell",
        "hela_dmem_comp": "HeLa DMEM 1% FDR compartment",
        "hela_hplm_comp": "HeLa HPLM 1% FDR compartment",
    }
    for key, label in hit_files.items():
        rows = read_csv(fetch(key))
        hits = [row for row in rows if row.get("Gene") == "CLDN4"]
        other = sorted({row["Gene"] for row in rows if row.get("Gene", "").startswith("CLDN")})
        checks.append(
            {
                "check": f"periscope_hit_{key}",
                "n_rows": len(rows),
                "cldn4_rows": len(hits),
                "detail": f"{label}; other CLDN* on this list: {', '.join(other) if other else 'none'}",
            }
        )

    for key, label in (
        ("hela_dmem_sig", "HeLa DMEM expressed-gene feature counts"),
        ("hela_hplm_sig", "HeLa HPLM expressed-gene feature counts"),
    ):
        rows = read_csv(fetch(key))
        sums = [float(row["Sum"]) for row in rows]
        hit = next(row for row in rows if row["Gene"] == "CLDN4")
        total = float(hit["Sum"])
        checks.append(
            {
                "check": f"periscope_{key}",
                "n_rows": len(rows),
                "cldn4_rows": 1,
                "detail": (
                    f"{label}; Mito={hit['Mito']} ConA={hit['ConA']} DAPI={hit['DAPI']} "
                    f"WGA={hit['WGA']} Phalloidin={hit['Phalloidin']} Sum={hit['Sum']}; "
                    f"Sum percentile<={percentile_leq(sums, total):.3f}; "
                    f"median Sum={sorted(sums)[len(sums)//2]:.0f}"
                ),
            }
        )

    zero = read_csv(fetch("hela_dmem_zero"))
    checks.append(
        {
            "check": "periscope_hela_dmem_0tpm_table",
            "n_rows": len(zero),
            "cldn4_rows": sum(1 for row in zero if row.get("Gene") == "CLDN4"),
            "detail": "CLDN4 is not in the 0 TPM feature-count table",
        }
    )

    rohban = read_csv(fetch("rohban"))
    rohban_hits = [
        row for row in rohban if row.get("Symbol") == "CLDN4" or "CLDN4" in "".join(row.values())
    ]
    checks.append(
        {
            "check": "rohban_orf_sequence_match_table",
            "n_rows": len(rohban),
            "cldn4_rows": len(rohban_hits),
            "detail": "carpenterlab 2017_rohban_elife matching.per.table.csv; Symbol column",
        }
    )

    luad = read_csv(fetch("luad_genes"))
    luad_genes = set()
    for row in luad:
        name = row["x_mutation_status"]
        if name.endswith("_ctl") or name in {"-666", ""}:
            continue
        luad_genes.add(name.split("_")[0])
    checks.append(
        {
            "check": "luad_allele_gene_funcs",
            "n_rows": len(luad),
            "cldn4_rows": int("CLDN4" in luad_genes),
            "detail": f"unique gene prefixes n={len(luad_genes)}; CLDN4 absent",
        }
    )

    rxrx_n = 0
    rxrx_rows = []
    with fetch("rxrx3").open() as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            rxrx_n += 1
            if row.get("gene") == "CLDN4":
                rxrx_rows.append(row)
    guides = Counter(row["treatment"] for row in rxrx_rows)
    experiments = Counter(row["experiment_name"] for row in rxrx_rows)
    checks.append(
        {
            "check": "rxrx3_core_metadata",
            "n_rows": rxrx_n,
            "cldn4_rows": len(rxrx_rows),
            "detail": (
                f"HUVEC CRISPR query guides; experiments {dict(experiments)}; "
                f"guides {dict(guides)}; expression column absent"
            ),
        }
    )

    check_path = OUT / "cldn4_checks.tsv"
    with check_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["check", "n_rows", "cldn4_rows", "detail"]
        )
        writer.writeheader()
        writer.writerows(checks)

    inventory = [
        {
            "dataset": "cpg0016 JUMP ORF",
            "class": "cell_painting_genetic",
            "system": "U2OS arrayed ORF",
            "cldn4_status": "overexpressed",
            "cldn4_expression": "no same-well RNA; DepMap U2OS TPM is an external gene-level annotation in the Morphmap paper and was not recomputed here",
            "n": f"{well_counts[ORF_ID]} wells",
            "evidence": f"{ORF_ID}; {orf_hits[0]['Metadata_Transcript']}; pLX_304; insert {orf_hits[0]['Metadata_Insert_Length']} nt",
            "url": "https://github.com/jump-cellpainting/datasets",
        },
        {
            "dataset": "cpg0016 JUMP CRISPR",
            "class": "cell_painting_genetic",
            "system": "U2OS arrayed CRISPR",
            "cldn4_status": "knockout reagent",
            "cldn4_expression": "no same-well RNA",
            "n": f"{well_counts[CRISPR_ID]} wells",
            "evidence": f"{CRISPR_ID}; source_13 plates CP-CC9-R1-06 through R5-06 well G18",
            "url": "https://github.com/jump-cellpainting/datasets",
        },
        {
            "dataset": "RxRx3-core",
            "class": "cell_painting_genetic",
            "system": "HUVEC CRISPR, six-channel modified Cell Painting",
            "cldn4_status": "unblinded knockout",
            "cldn4_expression": "no expression column in metadata",
            "n": f"{len(rxrx_rows)} wells, 5 guides",
            "evidence": "metadata_rxrx3_core.csv gene=CLDN4; experiments gene-003 and gene-176",
            "url": "https://huggingface.co/datasets/recursionpharma/rxrx3-core",
        },
        {
            "dataset": "cpg0021 PERISCOPE",
            "class": "cell_painting_genetic",
            "system": "A549 and HeLa pooled CRISPR, Cell Painting variant",
            "cldn4_status": "in the screened genome; not a 1% FDR hit",
            "cldn4_expression": (
                f"baseline DepMap 21Q4 only: A549 log2(TPM+1)={a549_log:.3f} (TPM={2 ** a549_log - 1:.2f}); "
                f"HeLa log2(TPM+1)={hela_log:.3f} (TPM={2 ** hela_log - 1:.2f})"
            ),
            "n": "genome-wide; A549 NGS lists 4 CLDN4 barcodes",
            "evidence": "absent from six published 1% FDR hit lists; HeLa feature-count sums near the expressed-gene median; DAPI channel is a DNA stain, not gamma-H2AX",
            "url": "https://github.com/broadinstitute/2022_PERISCOPE",
        },
        {
            "dataset": "cpg0003 Rosetta / Haghighi 2022",
            "class": "cell_painting_plus_l1000",
            "system": "five matched Cell Painting + L1000 sets (CDRP, CDRP-bio, LUAD alleles, TA-ORF, LINCS pilot)",
            "cldn4_status": "not a queried perturbation in Rohban or LUAD allele lists",
            "cldn4_expression": "LINCS geneinfo_beta feature_space=best inferred, not a 978-gene landmark. The 8.5 GB cpg0003 matrices were not opened, so a BING-inferred CLDN4 column was not confirmed or excluded by a file read",
            "n": "12,328 BING genes include CLDN4 as inferred",
            "evidence": "gene_id 1364; Haghighi et al. describe the released L1000 profiles as the 978 landmark genes",
            "url": "https://doi.org/10.1038/s41592-022-01667-0",
        },
        {
            "dataset": "cpg0017 Rohban TA-ORF / BBBC037",
            "class": "cell_painting_genetic",
            "system": "U2OS ORF overexpression",
            "cldn4_status": "absent",
            "cldn4_expression": "no",
            "n": f"{len(rohban)} constructs in the sequence-match table",
            "evidence": "matching.per.table.csv has no CLDN4; eLife 24060 text has no CLDN4",
            "url": "https://doi.org/10.7554/eLife.24060",
        },
        {
            "dataset": "cpg0031 Caicedo LUAD alleles",
            "class": "cell_painting_genetic",
            "system": "A549 ORF alleles, paired L1000 in the Rosetta set",
            "cldn4_status": "absent",
            "cldn4_expression": "L1000 may infer CLDN4 as a response gene; CLDN4 is not one of the perturbed genes",
            "n": f"{len(luad_genes)} gene prefixes in gene_funcs.csv",
            "evidence": "inputs/metadata/l1k/gene_funcs.csv",
            "url": "https://github.com/broadinstitute/luad-cell-painting",
        },
        {
            "dataset": "JUMP-Target-1 ORF plate",
            "class": "cell_painting_genetic",
            "system": "compound-target ORF subset",
            "cldn4_status": "absent",
            "cldn4_expression": "no",
            "n": f"{len(target_rows)} ORF metadata rows",
            "evidence": "JUMP-Target-1_orf_metadata.tsv",
            "url": "https://github.com/jump-cellpainting/JUMP-Target",
        },
        {
            "dataset": "idr0080 Cell Health",
            "class": "dna_damage_morphology",
            "system": "A549, ES2, HCC44 CRISPR; Cell Painting plus a separate cell-cycle panel",
            "cldn4_status": "absent from the perturbation table",
            "cldn4_expression": "no expression matrix",
            "n": f"{len(health)} perturbation rows; paper describes 59 genes / 119 guides",
            "evidence": "gamma-H2AX spot counts are measured; supplementary_table_1 has no CLDN4",
            "url": "https://idr.openmicroscopy.org/study/idr0080/",
        },
        {
            "dataset": "cpg0012 Bray/Wawer compounds",
            "class": "cell_painting_compound",
            "system": "U2OS, about 30,000 compounds",
            "cldn4_status": "not a genetic perturbation",
            "cldn4_expression": "no transcriptome in the gallery release",
            "n": "compound screen",
            "evidence": "Later DDR-prediction papers use these profiles as morphology only. DNA channel is not gamma-H2AX.",
            "url": "https://gigadb.org/dataset/100351",
        },
        {
            "dataset": "cpg0004 LINCS Cell Painting",
            "class": "cell_painting_compound",
            "system": "A549, 1,571 compounds, 6 doses",
            "cldn4_status": "not a genetic perturbation",
            "cldn4_expression": "not in this image release; any matched L1000 CLDN4 value would be best-inferred, not landmark",
            "n": "compound screen",
            "evidence": "Way et al. 2022 Cell Systems predict cell-health readouts, including gamma-H2AX, from these images. Predictions are not measurements.",
            "url": "https://github.com/broadinstitute/lincs-cell-painting",
        },
        {
            "dataset": "idr0133 Dahlin reference injury",
            "class": "cell_painting_compound",
            "system": "U2OS cytotoxic and nuisance compounds",
            "cldn4_status": "not a genetic perturbation",
            "cldn4_expression": "no",
            "n": "218 concentration-response plus 283 single-dose compounds",
            "evidence": "Cell Painting injury morphology, not a gamma-H2AX assay",
            "url": "https://idr.openmicroscopy.org/study/idr0133/",
        },
        {
            "dataset": "cpg0022 cmQTL",
            "class": "cell_painting_donor",
            "system": "297 iPSC donors",
            "cldn4_status": "not a perturbation",
            "cldn4_expression": "no public matched RNA-seq; WGS is dbGaP phs002032.v1.p1",
            "n": "297 lines",
            "evidence": "Tegtmeyer et al. 2024 Nat Commun. Morphology QTLs, not a CLDN4 expression matrix.",
            "url": "https://doi.org/10.1038/s41467-023-44045-w",
        },
        {
            "dataset": "Zenodo 7673199",
            "class": "dna_damage_morphology",
            "system": "gamma-H2AX microscopy, valinomycin pilot",
            "cldn4_status": "absent",
            "cldn4_expression": "no",
            "n": "small microscopy set",
            "evidence": "Deposit describes Hoechst plus pH2AX images for a valinomycin pilot. Archive was not unpacked; the record does not list a CLDN4 expression matrix",
            "url": "https://zenodo.org/records/7673199",
        },
        {
            "dataset": "Zenodo 13683162",
            "class": "dna_damage_morphology",
            "system": "nanomaterial HTS including gamma-H2AX and 8-oxoG",
            "cldn4_status": "absent",
            "cldn4_expression": "no",
            "n": "two assay collections in the deposit",
            "evidence": "Deposit describes CellTiter-Glo, DAPI, gamma-H2AX, 8-oxoG, and caspase assays for nanomaterials. Archive was not unpacked; the record does not list CLDN4",
            "url": "https://zenodo.org/records/13683162",
        },
    ]

    inv_path = OUT / "inventory_cellpainting.tsv"
    with inv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(inventory[0].keys()))
        writer.writeheader()
        writer.writerows(inventory)

    print(f"wrote {check_path}")
    print(f"wrote {inv_path}")
    for row in checks:
        if row["cldn4_rows"]:
            print(f"HIT {row['check']}: {row['detail']}")


if __name__ == "__main__":
    main()
