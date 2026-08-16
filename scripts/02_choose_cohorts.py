#!/usr/bin/env python3
"""Turn the lung-study screen into a download list.

Rules are fixed before any expression matrix is opened:

  * Bulk-ish tissue: n_tissue >= 12, single_cell_frac < 0.25, median_spots >= 5e5,
    and either lung_run_frac >= 0.40 or selection_basis == run_level.
  * Culture / organoid / ALI: n_culture_organoid >= 12, same sc/depth filters.
    These are kept as composition-controlled tests, not as bulk-tissue evidence.
  * Cell line: n_cell_line >= 12, same filters. Same caveat.
  * Sorted cells / fluid are never downloaded for the primary hunt.
  * Titles matching blood / platelet / isolated immune cells are dropped even if
    the regex labelled them tissue.

GTEx lung and TCGA LUAD/LUSC are added as non-SRA anchors. They are not chosen
because they look promising; they are the largest uniformly processed public
human lung RNA-seq matrices in recount3.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

DROP_TITLE = re.compile(
    r"\b(platelet|pbmc|whole blood|peripheral blood|t[- ]cell|b[- ]cell|"
    r"tumor-educated|lymph node t cell|dendritic cell|"
    r"bodyMap|promoter level mammalian expression atlas|"
    r"reference epigenomes of human gastrointestinal|"
    r"interferon signature induced by IFN|"
    r"HPA RNA-seq normal tissues)",
    re.I,
)

RECOUNT_ROOT = "https://recount-opendata.s3.amazonaws.com/recount3/release"


def gene_sums_url(organism: str, source: str, project: str) -> str:
    suffix = "G026" if organism == "human" else "M023"
    last2 = project[-2:]
    return (
        f"{RECOUNT_ROOT}/{organism}/data_sources/{source}/gene_sums/"
        f"{last2}/{project}/{source}.gene_sums.{project}.{suffix}.gz"
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--screen", default="results/hunt_tf/tables/studies_screened.tsv")
    ap.add_argument("--out", default="results/hunt_tf/tables/cohorts_to_download.tsv")
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.screen), delimiter="\t"))
    chosen: list[dict] = []

    for r in rows:
        if DROP_TITLE.search(r["study_title"] or ""):
            continue
        if float(r["single_cell_frac"]) >= 0.25:
            continue
        if int(r["median_spots"]) < 500_000:
            continue
        lungish = float(r["lung_run_frac"]) >= 0.40 or r["selection_basis"] == "run_level"
        if not lungish:
            continue

        n_tis = int(r["n_tissue"])
        n_cul = int(r["n_culture_organoid"])
        n_cl = int(r["n_cell_line"])
        if n_tis >= 12:
            material = "tissue"
            n = n_tis
        elif n_cul >= 12:
            material = "culture_or_organoid"
            n = n_cul
        elif n_cl >= 12:
            material = "cell_line"
            n = n_cl
        else:
            continue

        chosen.append(
            {
                "cohort_id": r["project"],
                "source": "sra",
                "organism": r["organism"],
                "material": material,
                "n_nominated": n,
                "n_runs_lung_rnaseq": r["n_runs_lung_rnaseq"],
                "single_cell_frac": r["single_cell_frac"],
                "median_spots": r["median_spots"],
                "study_title": r["study_title"],
                "url": gene_sums_url(r["organism"], "sra", r["project"]),
                "selected_runs": r["selected_runs"],
            }
        )

    anchors = [
        {
            "cohort_id": "GTEx_LUNG",
            "source": "gtex",
            "organism": "human",
            "material": "tissue",
            "n_nominated": 578,
            "n_runs_lung_rnaseq": 578,
            "single_cell_frac": "0.0",
            "median_spots": "0",
            "study_title": "GTEx v8 lung (recount3 G026 gene sums)",
            "url": gene_sums_url("human", "gtex", "LUNG"),
            "selected_runs": "",
        },
        {
            "cohort_id": "TCGA_LUAD",
            "source": "tcga",
            "organism": "human",
            "material": "tumor",
            "n_nominated": 585,
            "n_runs_lung_rnaseq": 585,
            "single_cell_frac": "0.0",
            "median_spots": "0",
            "study_title": "TCGA LUAD (recount3 G026 gene sums; tumor vs normal split later)",
            "url": gene_sums_url("human", "tcga", "LUAD"),
            "selected_runs": "",
        },
        {
            "cohort_id": "TCGA_LUSC",
            "source": "tcga",
            "organism": "human",
            "material": "tumor",
            "n_nominated": 550,
            "n_runs_lung_rnaseq": 550,
            "single_cell_frac": "0.0",
            "median_spots": "0",
            "study_title": "TCGA LUSC (recount3 G026 gene sums; tumor vs normal split later)",
            "url": gene_sums_url("human", "tcga", "LUSC"),
            "selected_runs": "",
        },
    ]
    chosen = anchors + chosen

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    cols = [
        "cohort_id", "source", "organism", "material", "n_nominated",
        "n_runs_lung_rnaseq", "single_cell_frac", "median_spots",
        "study_title", "url", "selected_runs",
    ]
    with out.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for rec in chosen:
            fh.write("\t".join(str(rec[c]).replace("\t", " ") for c in cols) + "\n")

    summary = {
        "n_cohorts": len(chosen),
        "by_organism": {
            org: sum(1 for c in chosen if c["organism"] == org)
            for org in ("human", "mouse")
        },
        "by_material": {},
    }
    for rec in chosen:
        summary["by_material"][rec["material"]] = summary["by_material"].get(rec["material"], 0) + 1
    Path(out.parent, "cohorts_to_download_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
