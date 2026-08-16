#!/usr/bin/env python3
"""Step 3 - platform feasibility: can each platform even measure TACSTD2 / CLDN4?

For every platform used by a shortlisted series we download the GEO platform
table (annot when available, otherwise the platform family SOFT) and search the
annotation columns for the target genes and for control genes, matching on
official symbol *and* known aliases (TROP2/M1S1/EGP-1 for TACSTD2,
CPE-R/CPETR1/WBSCR8 for CLDN4, ...).

This is the gate for the whole slice: a series whose platform lacks a probe for
a gene cannot inform that gene, no matter how good its clinical annotation is.

Outputs
  results/opus_microarray/platform_summary.tsv     technology + probe counts per platform
  results/opus_microarray/platform_probes.tsv      every matched probe (id, symbol, column)
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (  # noqa: E402
    CACHE_DIR,
    CONTROL_GENES,
    RESULTS_DIR,
    TARGET_GENES,
    alias_set,
    ensure_dirs,
    fetch,
    gpl_annot_url,
    gpl_soft_url,
    norm_gene,
    open_maybe_gzip,
)

ACC_CGI = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"

# Platforms of the shortlisted / documented series.
PLATFORMS = {
    "GPL19965": "GSE93157, GSE140901 (NanoString nCounter PanCancer Immune 730-gene panel)",
    "GPL30173": "GSE261345, GSE261348 (ES-SCLC chemo-immunotherapy multi-region transcriptomics)",
    "GPL23126": "GSE202417, GSE248249, GSE141479 (Affymetrix Clariom D Human)",
    "GPL570":   "GSE305086 (Affymetrix HG-U133 Plus 2.0, whole blood under ICI)",
    "GPL29738": "GSE180347 (PD-L1 +/- 'hot' lung adenocarcinoma, no ICI treatment)",
    "GPL14951": "GSE67501 (Illumina HumanHT-12 WG-DASL, renal cell carcinoma anti-PD-1)",
    "GPL10558": "GSE99070 (Illumina HumanHT-12 v4, malignant pleural mesothelioma)",
}

GENES = TARGET_GENES + CONTROL_GENES
SPLIT_RE = re.compile(r"[;,/|]| // ")


def platform_header(gpl: str) -> dict[str, str]:
    dest = CACHE_DIR / "platform" / f"{gpl}_self_brief.txt"
    fetch(f"{ACC_CGI}?acc={gpl}&targ=self&form=text&view=brief", dest)
    info: dict[str, str] = {}
    for line in dest.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("!Platform_") and "=" in line:
            k, v = line.split("=", 1)
            k = k.replace("!Platform_", "").strip()
            info.setdefault(k, v.strip())
    return info


def platform_table(gpl: str) -> Path | None:
    """Prefer the curated .annot.gz; fall back to the submitter SOFT table."""
    for url_fn, name in ((gpl_annot_url, "annot"), (gpl_soft_url, "soft")):
        dest = CACHE_DIR / "platform" / f"{gpl}_{name}.gz"
        try:
            return fetch(url_fn(gpl), dest)
        except Exception as exc:  # noqa: BLE001
            print(f"  [miss] {gpl} {name}: {exc}")
    # Last resort: full platform table straight from acc.cgi
    dest = CACHE_DIR / "platform" / f"{gpl}_table.txt"
    try:
        return fetch(f"{ACC_CGI}?acc={gpl}&targ=self&form=text&view=data", dest)
    except Exception as exc:  # noqa: BLE001
        print(f"  [miss] {gpl} acc.cgi table: {exc}")
    return None


def scan_table(gpl: str, path: Path) -> tuple[list[dict], int, list[str]]:
    """Walk the platform table; return probe hits, row count and column names."""
    hits: list[dict] = []
    header: list[str] | None = None
    n_rows = 0
    want = {g: alias_set(g) for g in GENES}
    with open_maybe_gzip(path) as fh:
        in_table = False
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(("^", "!", "#")):
                if line.startswith("!platform_table_begin"):
                    in_table = True
                elif line.startswith("!platform_table_end"):
                    break
                continue
            if header is None:
                header = line.split("\t")
                continue
            if not line:
                continue
            n_rows += 1
            fields = line.split("\t")
            probe_id = fields[0]
            for gene, aliases in want.items():
                # Include the ID column: NanoString / custom panels often use the
                # official symbol as the probe identifier and leave Alias empty.
                for ci, cell in enumerate(fields):
                    if not cell or len(cell) > 4000:
                        continue
                    toks = {norm_gene(t) for t in SPLIT_RE.split(cell)}
                    if toks & aliases:
                        col = header[ci] if ci < len(header) else f"col{ci}"
                        hits.append(
                            {
                                "platform": gpl,
                                "gene": gene,
                                "probe_id": probe_id,
                                "matched_column": col,
                                "matched_value": cell[:160],
                            }
                        )
                        break
    _ = in_table
    return hits, n_rows, header or []


def main() -> int:
    ensure_dirs()
    all_hits: list[dict] = []
    summary: list[dict] = []
    for gpl, used_by in PLATFORMS.items():
        print(f"[platform] {gpl} - {used_by}", flush=True)
        info = platform_header(gpl)
        path = platform_table(gpl)
        if path is None:
            summary.append({"platform": gpl, "used_by": used_by, "status": "no_table"})
            continue
        hits, n_rows, header = scan_table(gpl, path)
        all_hits.extend(hits)
        per_gene = {g: sum(1 for h in hits if h["gene"] == g) for g in GENES}
        summary.append(
            {
                "platform": gpl,
                "used_by": used_by,
                "platform_title": info.get("title", ""),
                "technology": info.get("technology", ""),
                "organism": info.get("organism", ""),
                "table_rows": n_rows,
                "table_columns": ";".join(header),
                "source_file": path.name,
                **{f"n_probes_{g}": per_gene[g] for g in GENES},
                "targets_measurable": int(all(per_gene[g] > 0 for g in TARGET_GENES)),
            }
        )
        print(
            f"    rows={n_rows}  "
            + "  ".join(f"{g}={per_gene[g]}" for g in TARGET_GENES + ["CD274", "CD8A"])
        )

    sm = pd.DataFrame(summary)
    sm.to_csv(RESULTS_DIR / "platform_summary.tsv", sep="\t", index=False)
    ph = pd.DataFrame(all_hits)
    ph.to_csv(RESULTS_DIR / "platform_probes.tsv", sep="\t", index=False)
    print(f"\n[write] platform_summary.tsv ({len(sm)} platforms), platform_probes.tsv ({len(ph)} probe hits)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
