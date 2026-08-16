#!/usr/bin/env python3
"""Turn the FTP probe into a reviewable shortlist of usable lung-ICI expression cohorts.

Applies the task's hard constraints:
  * open access at the series level (a processed supplementary file or a populated series matrix),
  * no FASTQ/raw-only series, nothing above 2 GB,
  * human lung tumour material (blood/platelet/stool compartments are flagged, not silently used),
  * bulk rather than single-cell/spatial,
and records the reason each series was kept or dropped so the audit trail is explicit.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "opus_geo_leftover"
MAX_BYTES = 2 * 1024**3

LUNG = re.compile(r"\blung\b|NSCLC|non[- ]small[- ]cell|\bSCLC\b|small[- ]cell lung|LUAD|LUSC|pulmonary", re.I)
ICI = re.compile(
    r"nivolumab|pembrolizumab|atezolizumab|durvalumab|avelumab|ipilimumab|cemiplimab|"
    r"camrelizumab|sintilimab|tislelizumab|toripalimab|anti-?PD-?1|anti-?PD-?L1|anti-?PD-?\(L\)1|"
    r"immunotherap|immune checkpoint|checkpoint (inhibit|block)|\bICI\b|\bICB\b|PD-?1/PD-?L1",
    re.I,
)
RESPONSE_FIELD = re.compile(
    r"(best )?response|responder|\bRECIST\b|\bBOR\b|clinical benefit|\bDCB\b|\bMPR\b|"
    r"pathologic(al)? response|\bpCR\b|recurrence|relapse|progression|efficacy|outcome|"
    r"os \(|pfs \(|overall survival|progression.free",
    re.I,
)
NON_TUMOUR = re.compile(
    r"thrombocyte|platelet|\bPBMC\b|peripheral blood|whole blood|\bserum\b|plasma|\bstool\b|"
    r"faecal|fecal|saliva|\bBALF\b|bronchoalveolar|urine|buffy",
    re.I,
)
CELL_LINE = re.compile(
    r"\bA549\b|\bH1975\b|H1299|HCC827|PC-?9|\bH460\b|\bH358\b|cell line|xenograft|organoid|spheroid|"
    r"\bshRNA\b|\bsiRNA\b|knockdown|knockout|overexpress|CRISPR|transfect|treated with .* for \d+ ?h",
    re.I,
)
SINGLE_CELL = re.compile(
    r"single[- ]cell|scRNA|snRNA|10x|smart-seq|CITE-seq|single[- ]nucleus|spatial|Visium|GeoMx|"
    r"digital spatial|CosMx|\bDSP\b|\bROI\b",
    re.I,
)


def main() -> int:
    probe = json.loads((OUT / "geo_probe.json").read_text())
    cands = {c["accession"]: c for c in json.loads((OUT / "geo_triage_scored.json").read_text())}

    kept, dropped = [], []
    for r in probe:
        acc = r["accession"]
        cand = cands.get(acc, {})
        chars = " ".join(r["characteristics"])
        blob = " ".join([r["title"], cand.get("summary", ""), chars, cand.get("sample_titles", "")])
        n_chars_fields = max(
            (len(c.split("\t")) - 1 for c in r["characteristics"]), default=0
        )

        reasons = []
        if not LUNG.search(blob):
            reasons.append("not lung")
        if not ICI.search(blob):
            reasons.append("no ICI mention")
        if not RESPONSE_FIELD.search(chars):
            reasons.append("no response/outcome field in sample characteristics")
        if SINGLE_CELL.search(blob):
            reasons.append("single-cell or spatial")
        if CELL_LINE.search(blob) and n_chars_fields < 15:
            reasons.append("cell line / model system")
        if NON_TUMOUR.search(chars) and not re.search(r"tissue: (tumor|tumour|lung|NSCLC)", chars, re.I):
            reasons.append("non-tumour compartment")
        if n_chars_fields < 10:
            reasons.append(f"too few samples with metadata (n={n_chars_fields})")

        usable = [f for f in r["processed_candidates"] if (f["bytes"] or 0) <= MAX_BYTES]
        has_open_values = bool(usable) or r["has_series_matrix_values"]
        if not has_open_values:
            reasons.append("no open processed matrix (raw/FASTQ only or controlled access)")
        if r["raw_only"]:
            reasons.append("supplementary files are raw/FASTQ only")
        if r["oversized_files"] and not usable:
            reasons.append(">2GB only")

        rec = {
            "accession": acc,
            "title": r["title"],
            "n_samples_with_metadata": n_chars_fields,
            "pubmed": r["pubmed"],
            "platform_ids": r["platform_ids"],
            "usable_files": usable,
            "oversized_files": r["oversized_files"],
            "series_matrix_rows": r["matrix_data_rows"],
            "reasons_excluded": reasons,
        }
        (kept if not reasons else dropped).append(rec)

    kept.sort(key=lambda x: -x["n_samples_with_metadata"])
    (OUT / "dataset_shortlist.json").write_text(json.dumps(kept, indent=2))
    (OUT / "dataset_excluded.json").write_text(json.dumps(dropped, indent=2))

    print(f"{'ACCESSION':<12} {'N':>4}  {'SMXROWS':>8}  FILES / TITLE")
    for k in kept:
        files = ",".join(f["name"] for f in k["usable_files"])[:70]
        print(f"{k['accession']:<12} {k['n_samples_with_metadata']:>4}  {k['series_matrix_rows']:>8}  {files}")
        print(f"{'':<12} {'':>4}  {'':>8}  {k['title'][:100]}")
    print(f"\n[select] kept {len(kept)}, dropped {len(dropped)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
