#!/usr/bin/env python3
"""Audit leftover Zenodo/figshare lung ICI records for TACSTD2/CLDN4."""

from __future__ import annotations

import csv
import hashlib
import json
import statistics
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "source_data"
TARGETS = ("TACSTD2", "CLDN4")
CONTEXT_GENES = ("EPCAM", "CDH1", "MUC1", "CEACAM6")
EXPECTED_MD5 = {
    "learn-data.csv": "90dcff5eba025880a7b55a925165ff96",
    "learn-metadata.csv": "57f426b6e54f47827fd88b5481e761b9",
    "validate-data.csv": "69d0f9ca7472ff459eebdc467e166764",
    "validate-metadata.csv": "4a9016ad18ce0cda0cad352182ddc382",
}
GEOMX_EXTRACT = ROOT / "geomx_tacstd2_rois.csv"
ROI_ORDER = ("immunity hub", "hybrid hub", "bystander TLS", "non-hub", "exclude")


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        return list(reader.fieldnames or []), list(reader)


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def count_text(rows: list[dict[str, str]], field: str) -> str:
    counts = Counter((row.get(field) or "missing").strip() or "missing" for row in rows)
    return json.dumps(dict(sorted(counts.items())), sort_keys=True)


def audit_nanostring() -> dict[str, object]:
    manifest_rows: list[dict[str, object]] = []
    for name, expected in EXPECTED_MD5.items():
        observed = md5(DATA / name)
        if observed != expected:
            raise RuntimeError(f"Checksum mismatch for {name}: {observed} != {expected}")
        manifest_rows.append(
            {
                "file": name,
                "md5": observed,
                "bytes": (DATA / name).stat().st_size,
                "verified": True,
            }
        )

    coverage_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    for cohort, stem in (("discovery", "learn"), ("validation", "validate")):
        expression_fields, expression = read_csv(DATA / f"{stem}-data.csv")
        metadata_fields, metadata = read_csv(DATA / f"{stem}-metadata.csv")
        if len(expression) != len(metadata):
            raise RuntimeError(
                f"{cohort}: expression rows ({len(expression)}) != metadata rows ({len(metadata)})"
            )
        if [row[""] for row in expression] != [row[""] for row in metadata]:
            raise RuntimeError(f"{cohort}: expression and metadata row identifiers do not align")

        genes = set(expression_fields)
        for gene in TARGETS + CONTEXT_GENES:
            coverage_rows.append(
                {
                    "cohort": cohort,
                    "gene": gene,
                    "requested_target": gene in TARGETS,
                    "present_exact_symbol": gene in genes,
                    "n_samples": len(expression),
                }
            )

        summary_rows.append(
            {
                "cohort": cohort,
                "n_samples": len(metadata),
                "n_expression_features_including_controls": len(expression_fields) - 1,
                "best_response_binary_counts": count_text(metadata, "best response binary"),
                "best_response_raw_counts": count_text(metadata, "best response IO"),
                "histology_counts": count_text(metadata, "entity"),
                "metadata_columns": len(metadata_fields) - 1,
            }
        )

    write_csv(ROOT / "source_manifest.csv", manifest_rows, ["file", "md5", "bytes", "verified"])
    write_csv(
        ROOT / "target_coverage.csv",
        coverage_rows,
        ["cohort", "gene", "requested_target", "present_exact_symbol", "n_samples"],
    )
    write_csv(
        ROOT / "cohort_summary.csv",
        summary_rows,
        [
            "cohort",
            "n_samples",
            "n_expression_features_including_controls",
            "best_response_binary_counts",
            "best_response_raw_counts",
            "histology_counts",
            "metadata_columns",
        ],
    )
    return {
        "record": "10.5281/zenodo.2635194",
        "assay": "NanoString nCounter PanCancer Immune Profiling panel",
        "population": "advanced NSCLC treated with anti-PD-1 immunotherapy",
        "n_discovery": summary_rows[0]["n_samples"],
        "n_validation": summary_rows[1]["n_samples"],
        "targets": {
            target: {
                "present_discovery": any(
                    row["cohort"] == "discovery"
                    and row["gene"] == target
                    and row["present_exact_symbol"]
                    for row in coverage_rows
                ),
                "present_validation": any(
                    row["cohort"] == "validation"
                    and row["gene"] == target
                    and row["present_exact_symbol"]
                    for row in coverage_rows
                ),
            }
            for target in TARGETS
        },
        "valid_target_analysis": False,
        "reason": "Neither requested gene is measured by the deposited targeted panel.",
        "surrogate_analysis_performed": False,
    }


def group_stats(rows: list[dict[str, str]], key: str) -> list[dict[str, object]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row[key]].append(row)
    out: list[dict[str, object]] = []
    for name in sorted(grouped, key=lambda item: (ROI_ORDER.index(item) if item in ROI_ORDER else 99, item)):
        subset = grouped[name]
        values = [float(row["q_norm"]) for row in subset]
        loqs = [float(row["LOQ.Hs_R_NGS_CTA_v1.0"]) for row in subset]
        above = sum(row["above_LOQ"] == "TRUE" for row in subset)
        out.append(
            {
                key: name,
                "n_rois": len(subset),
                "n_above_loq": above,
                "frac_above_loq": round(above / len(subset), 6),
                "median_q_norm": statistics.median(values),
                "median_loq": statistics.median(loqs),
                "q_norm_used_as_quantitative_endpoint": False,
            }
        )
    return out


def write_detection_svg(by_roi: list[dict[str, object]], by_slide: list[dict[str, object]]) -> None:
    def bars(rows: list[dict[str, object]], label_key: str, x0: float, title: str) -> str:
        width = 360
        height = 220
        max_frac = max((float(row["frac_above_loq"]) for row in rows), default=0) or 1
        bar_w = 48
        gap = 18
        parts = [
            f'<text x="{x0}" y="18" font-family="DejaVu Sans, Arial, sans-serif" font-size="13">{title}</text>'
        ]
        for idx, row in enumerate(rows):
            frac = float(row["frac_above_loq"])
            h = 140 * frac / max_frac
            x = x0 + 20 + idx * (bar_w + gap)
            y = 180 - h
            parts.append(
                f'<rect x="{x}" y="{y:.1f}" width="{bar_w}" height="{h:.1f}" fill="#4C78A8"/>'
            )
            parts.append(
                f'<text x="{x + bar_w / 2:.1f}" y="198" text-anchor="middle" font-size="9" '
                f'font-family="DejaVu Sans, Arial, sans-serif">{row[label_key]}</text>'
            )
            parts.append(
                f'<text x="{x + bar_w / 2:.1f}" y="{y - 6:.1f}" text-anchor="middle" font-size="10" '
                f'font-family="DejaVu Sans, Arial, sans-serif">{row["n_above_loq"]}/{row["n_rois"]}</text>'
            )
        return "\n".join(parts)

    roi_rows = [row for row in by_roi if row["ROI_Type"] != "exclude"]
    svg = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="920" height="260" viewBox="0 0 920 260">',
        '<rect width="920" height="260" fill="white"/>',
        '<text x="16" y="20" font-family="DejaVu Sans, Arial, sans-serif" font-size="14">'
        "Leftover GeoMx TACSTD2 above LOQ (not an ICI-response test)</text>",
        bars(roi_rows, "ROI_Type", 16, "By ROI type"),
        bars(by_slide, "Slide_ID", 470, "By slide"),
        "</svg>",
        "",
    ]
    (ROOT / "geomx_tacstd2_above_loq.svg").write_text("\n".join(svg), encoding="utf-8")


def analyze_geomx() -> dict[str, object]:
    fields, rows = read_csv(GEOMX_EXTRACT)
    required = {
        "SampleName",
        "Slide_ID",
        "roi",
        "ROI_Type",
        "q_norm",
        "LOQ.Hs_R_NGS_CTA_v1.0",
        "above_LOQ",
    }
    missing = required.difference(fields)
    if missing:
        raise RuntimeError(f"GeoMx extract missing columns: {sorted(missing)}")
    if len(rows) != 285:
        raise RuntimeError(f"Expected 285 TACSTD2 ROI rows, found {len(rows)}")

    by_roi = group_stats(rows, "ROI_Type")
    by_slide = group_stats(rows, "Slide_ID")
    write_csv(
        ROOT / "geomx_detection_by_roi_type.csv",
        by_roi,
        [
            "ROI_Type",
            "n_rois",
            "n_above_loq",
            "frac_above_loq",
            "median_q_norm",
            "median_loq",
            "q_norm_used_as_quantitative_endpoint",
        ],
    )
    write_csv(
        ROOT / "geomx_detection_by_slide.csv",
        by_slide,
        [
            "Slide_ID",
            "n_rois",
            "n_above_loq",
            "frac_above_loq",
            "median_q_norm",
            "median_loq",
            "q_norm_used_as_quantitative_endpoint",
        ],
    )
    write_detection_svg(by_roi, by_slide)

    n_above = sum(row["above_LOQ"] == "TRUE" for row in rows)
    slide_above = Counter(row["Slide_ID"] for row in rows if row["above_LOQ"] == "TRUE")
    dominant_slide, dominant_n = slide_above.most_common(1)[0]
    return {
        "record": "10.5281/zenodo.11198494",
        "file": "geomx.csv",
        "source_md5": "8421f72c459e8aa6eb06fa48ba974501",
        "assay": "GeoMx Cancer Transcriptome Atlas plus TCR module",
        "population": "immunotherapy-naive primary NSCLC spatial ROIs from a PD-1-response spatial study",
        "n_features_in_long_table": 2018,
        "n_rois": len(rows),
        "n_slides": len({row["Slide_ID"] for row in rows}),
        "targets": {
            "TACSTD2": {"present": True, "n_rois_above_loq": n_above},
            "CLDN4": {"present": False, "n_rois_above_loq": None},
        },
        "valid_ici_outcome_analysis": False,
        "valid_quantitative_expression_analysis": False,
        "reason": (
            "TACSTD2 is measured but above LOQ in only "
            f"{n_above}/{len(rows)} ROIs, and {dominant_n}/{n_above} detections "
            f"come from one slide ({dominant_slide}). The GeoMx table has no ICI "
            "response labels. CLDN4 is absent."
        ),
        "surrogate_analysis_performed": False,
    }


def main() -> None:
    nanostring = audit_nanostring()
    geomx = analyze_geomx()
    result = {
        "question": "leftover processed Zenodo/figshare lung ICI matrices for TACSTD2/CLDN4 not analyzed in PR42",
        "valid_target_vs_ici_result": False,
        "nanostring_zenodo_2635194": nanostring,
        "geomx_zenodo_11198494": geomx,
        "honest_conclusion": (
            "No leftover open compact matrix supports a TACSTD2 or CLDN4 comparison "
            "against ICI outcome. The NanoString ICI cohort lacks both genes. The "
            "leftover GeoMx table measures TACSTD2 but not CLDN4, is mostly below "
            "LOQ, is slide-confounded, and has no response labels."
        ),
    }
    (ROOT / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
