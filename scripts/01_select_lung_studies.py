#!/usr/bin/env python3
"""Select candidate lung RNA-seq studies from the recount3 SRA metadata dump.

Reads the per-study SRA metadata files fetched by scripts/00_fetch_metadata.sh and
writes one row per study describing how many runs look like lung material, what kind
of material it is (tissue / cell line / culture-organoid), and whether the study looks
like single-cell data. Nothing here decides that a study supports the hypothesis; this
step only defines the search space, and every study that is skipped is written out with
the reason so the selection can be audited.
"""
from __future__ import annotations

import argparse
import gzip
import json
import re
from pathlib import Path

# Free-text evidence that a sample is lung / respiratory material. Deliberately broad:
# precision is recovered later by the material-type regexes and by manual review of the
# shortlist.
LUNG_RE = re.compile(
    r"\b("
    r"lung|lungs|pulmonary|pulmo|alveol|airway|bronch|"
    r"nsclc|luad|lusc|lung adenocarcinoma|lung squamous|small cell lung|"
    r"idiopathic pulmonary fibrosis|\bipf\b|copd|emphysema|bleomycin|"
    r"pneumocyte|pneumonia|pneumonitis|asthma|respiratory epitheli"
    r")",
    re.I,
)

# Tissues that, when they dominate a study, mean the "lung" hit came from a comparison
# arm rather than from the samples themselves.
OTHER_TISSUE_RE = re.compile(
    r"\b(liver|kidney|brain|cortex|heart|muscle|spleen|thymus|testis|ovary|"
    r"pancreas|colon|intestine|stomach|skin|breast|prostate|bladder|"
    r"blood|pbmc|bone marrow|adipose|placenta|retina|esophag)",
    re.I,
)

CELL_LINE_RE = re.compile(
    r"\b(a549|nci-?h\d{2,4}|calu-?\d|beas-?2b|pc-?9|hcc827|h1299|h1975|h358|h441|h226|"
    r"16hbe|hbec3|nl20|mle-?12|mle15|la-?4|lle|mrc-?5|imr-?90|wi-?38|"
    r"hela|hek293|cell line|cell-line|celline)",
    re.I,
)

# Applied to per-run text only. A tissue study that mentions macrophages in the
# abstract must not be discarded.
SORTED_RUN_RE = re.compile(
    r"\b(cd4\+|cd8\+|cd45\+|cd3\+|treg|th2|th17|"
    r"pbmc|whole blood|peripheral blood|plasma|serum|exosome|"
    r"bronchoalveolar lavage|\bbal fluid\b|sputum|"
    r"facs[- ]sorted|flow[- ]sorted|magnetic bead)",
    re.I,
)

# Applied to the study title only. These designs are not bulk lung parenchyma.
SORTED_TITLE_RE = re.compile(
    r"\b(t[- ]cell|b[- ]cell|treg|th2|th17|ilc\d?|nk cell|"
    r"macrophage|monocyte|neutrophil|eosinophil|platelet|"
    r"pbmc|whole blood|peripheral blood|lymph node|"
    r"bronchoalveolar lavage|\bbal\b|sputum|tumor-educated platelet)",
    re.I,
)

CULTURE_RE = re.compile(
    r"\b(organoid|air-?liquid interface|\bali\b|nhbe|hbec|primary .{0,25}epithelial cell|"
    r"ipsc|ips cell|embryonic stem|esc-derived|spheroid|explant|3d culture|"
    r"differentiat\w* (cell|culture)|passage \d)",
    re.I,
)

SINGLE_CELL_RE = re.compile(
    r"(single[- ]cell|single[- ]nucle|scrna|snrna|sc-?rna-?seq|smart-?seq|"
    r"10x genomics|10x chromium|drop-?seq|cel-?seq|mars-?seq|inDrop|"
    r"single cell rna)",
    re.I,
)

RUN_TEXT_FIELDS = [
    "sample_attributes",
    "sample_title",
    "sample_description",
    "experiment_title",
    "sample_name",
    "library_name",
    "design_description",
]
STUDY_TEXT_FIELDS = ["study_title", "study_abstract", "study_description"]


def read_study(path: Path) -> tuple[list[str], list[list[str]]]:
    with gzip.open(path, "rt", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = [line.rstrip("\n").split("\t") for line in fh]
    return header, rows


def classify_study(project: str, organism: str, path: Path) -> dict | None:
    header, rows = read_study(path)
    if not rows:
        return None
    idx = {name: i for i, name in enumerate(header)}

    def cell(row: list[str], name: str) -> str:
        i = idx.get(name)
        if i is None or i >= len(row):
            return ""
        return row[i]

    study_text = " ".join(cell(rows[0], f) for f in STUDY_TEXT_FIELDS)
    study_title = cell(rows[0], "study_title")

    runs: list[dict] = []
    for row in rows:
        run_text = " ".join(cell(row, f) for f in RUN_TEXT_FIELDS)
        full_text = f"{run_text} {study_text}"
        if CELL_LINE_RE.search(run_text):
            material = "cell_line"
        elif SORTED_RUN_RE.search(run_text) or SORTED_TITLE_RE.search(study_title):
            material = "sorted_cells_or_fluid"
        elif CULTURE_RE.search(run_text) or CULTURE_RE.search(study_text):
            material = "culture_or_organoid"
        else:
            material = "tissue"
        runs.append(
            {
                "run": cell(row, "external_id"),
                "lung_run": bool(LUNG_RE.search(run_text)),
                "other_tissue_run": bool(OTHER_TISSUE_RE.search(run_text)),
                "material": material,
                "strategy": cell(row, "library_strategy"),
                "source": cell(row, "library_source"),
                "single_cell": bool(SINGLE_CELL_RE.search(full_text)),
                "spots": cell(row, "num_spots") or cell(row, "sample_spots") or "0",
            }
        )

    n_runs = len(runs)
    lung_run_frac = sum(r["lung_run"] for r in runs) / n_runs
    study_lung = bool(LUNG_RE.search(study_text))

    # A run is analysed if its own text says lung; if no run-level text mentions lung but
    # the study does and no run mentions a different tissue, the whole study is treated as
    # lung (common for studies whose per-sample attributes are uninformative).
    if lung_run_frac > 0:
        selected = [r for r in runs if r["lung_run"]]
        basis = "run_level"
    elif study_lung and not any(r["other_tissue_run"] for r in runs):
        selected = runs
        basis = "study_level"
    else:
        selected = []
        basis = "none"

    rnaseq = [r for r in selected if r["strategy"] == "RNA-Seq" and r["source"] == "TRANSCRIPTOMIC"]
    materials = {
        m: sum(1 for r in rnaseq if r["material"] == m)
        for m in ("tissue", "cell_line", "culture_or_organoid", "sorted_cells_or_fluid")
    }
    single_cell_frac = (sum(r["single_cell"] for r in rnaseq) / len(rnaseq)) if rnaseq else 0.0
    median_spots = 0
    if rnaseq:
        spots = sorted(int(r["spots"] or 0) for r in rnaseq)
        median_spots = spots[len(spots) // 2]

    return {
        "project": project,
        "organism": organism,
        "n_runs_total": n_runs,
        "n_runs_lung_selected": len(selected),
        "n_runs_lung_rnaseq": len(rnaseq),
        "selection_basis": basis,
        "lung_run_frac": round(lung_run_frac, 3),
        "study_level_lung": study_lung,
        "n_tissue": materials["tissue"],
        "n_cell_line": materials["cell_line"],
        "n_culture_organoid": materials["culture_or_organoid"],
        "n_sorted_or_fluid": materials["sorted_cells_or_fluid"],
        "single_cell_frac": round(single_cell_frac, 3),
        "median_spots": median_spots,
        "study_title": (rows[0][idx["study_title"]] if "study_title" in idx else "")[:300],
        "selected_runs": ",".join(r["run"] for r in rnaseq),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--meta-dir", default=str(Path.home() / "hunt_data" / "meta"))
    ap.add_argument("--out", default="results/hunt_tf/tables/studies_screened.tsv")
    args = ap.parse_args()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    cols = [
        "project", "organism", "n_runs_total", "n_runs_lung_selected", "n_runs_lung_rnaseq",
        "selection_basis", "lung_run_frac", "study_level_lung", "n_tissue", "n_cell_line",
        "n_culture_organoid", "n_sorted_or_fluid", "single_cell_frac", "median_spots", "study_title", "selected_runs",
    ]
    summary = {}
    with out_path.open("w") as out:
        out.write("\t".join(cols) + "\n")
        for organism in ("human", "mouse"):
            files = sorted(Path(args.meta_dir, organism).glob("*.MD.gz"))
            kept = 0
            for path in files:
                rec = classify_study(path.stem.replace(".MD", ""), organism, path)
                if rec is None or rec["n_runs_lung_rnaseq"] == 0:
                    continue
                kept += 1
                out.write("\t".join(str(rec[c]).replace("\t", " ") for c in cols) + "\n")
            summary[organism] = {"studies_scanned": len(files), "studies_with_lung_rnaseq_runs": kept}
            print(f"{organism}: scanned {len(files)} studies, {kept} have >=1 lung RNA-seq run")

    Path(out_path.parent, "studies_screened_summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
