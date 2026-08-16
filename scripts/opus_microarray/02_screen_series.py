#!/usr/bin/env python3
"""Step 2 - screen candidate series on *sample-level* metadata.

Series-level titles/summaries are unreliable for deciding whether a GEO Series
actually contains per-patient checkpoint-inhibitor response labels, so we pull
the brief GSM metadata for every plausible hit
(`acc.cgi?acc=GSExxx&targ=gsm&form=text&view=brief`, no data tables => small)
and score each series on:

  * human tumour/blood patient material (vs cell line / xenograft / mouse)
  * per-sample ICI drug or ICI-treatment annotation
  * per-sample response annotation (RECIST, responder/non-responder, DCB/NDB, ...)
  * assay technology of the platform (array vs sequencer vs targeted panel)

Outputs
  results/opus_microarray/series_screen.tsv         one row per screened series
  results/opus_microarray/series_screen_evidence/   per-series matched metadata lines
"""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import CACHE_DIR, RESULTS_DIR, ensure_dirs, fetch  # noqa: E402

ACC_CGI = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"

RESPONSE_KEY_RE = re.compile(
    r"(best[\s._-]*resp|^resp|response|recist|responder|clinical[\s._-]*benefit|dcb|ndb|"
    r"objective[\s._-]*resp|treatment[\s._-]*outcome|outcome|efficacy|progression|pfs|"
    r"os[\s._-]*month|survival|durable)",
    re.I,
)
RESPONSE_VAL_RE = re.compile(
    r"\b(CR|PR|SD|PD|complete response|partial response|stable disease|progressive disease|"
    r"responder|non[- ]?responder|nonresponder|DCB|NDB|sensitive|refractory|resistant)\b"
)
ICI_RE = re.compile(
    r"nivolumab|pembrolizumab|atezolizumab|durvalumab|avelumab|ipilimumab|cemiplimab|"
    r"tislelizumab|camrelizumab|sintilimab|toripalimab|anti[-\s]?pd[-\s]?1|anti[-\s]?pd[-\s]?l1|"
    r"pd-?1 blockade|checkpoint (inhibitor|blockade)|immunotherapy|immune checkpoint",
    re.I,
)
LUNG_RE = re.compile(r"\blung\b|nsclc|non[-\s]?small[-\s]?cell|adenocarcinom|squamous|luad|lusc|pulmonary", re.I)
CELLLINE_RE = re.compile(
    r"cell line|A549|H1975|H1299|PC9|PC-9|HCC827|Calu-?[0-9]|xenograft|organoid|"
    r"mus musculus|mouse|murine|in vitro|transfect|knockdown|overexpress|siRNA|shRNA",
    re.I,
)


def brief_metadata(gse: str) -> str:
    dest = CACHE_DIR / "brief" / f"{gse}_gsm_brief.txt"
    url = f"{ACC_CGI}?acc={gse}&targ=gsm&form=text&view=brief"
    try:
        fetch(url, dest)
    except Exception as exc:  # noqa: BLE001
        return f"__FETCH_ERROR__ {exc}"
    return dest.read_text(encoding="utf-8", errors="replace")


def screen_one(gse: str, text: str) -> dict:
    lines = text.splitlines()
    n_gsm = sum(1 for l in lines if l.startswith("^SAMPLE"))
    char_lines = [l for l in lines if l.startswith("!Sample_characteristics_ch")]
    title_lines = [l for l in lines if l.startswith("!Sample_title")]
    source_lines = [l for l in lines if l.startswith("!Sample_source_name_ch")]
    org_lines = [l for l in lines if l.startswith("!Sample_organism_ch")]
    plat_lines = [l for l in lines if l.startswith("!Sample_platform_id")]
    platforms = sorted({l.split("=", 1)[1].strip() for l in plat_lines if "=" in l})
    organisms = sorted({l.split("=", 1)[1].strip() for l in org_lines if "=" in l})

    # characteristic keys, e.g. "!Sample_characteristics_ch1 = best.resp: SD" -> "best.resp"
    keys: set[str] = set()
    for l in char_lines:
        if "=" not in l:
            continue
        val = l.split("=", 1)[1].strip()
        if ":" in val:
            keys.add(val.split(":", 1)[0].strip().lower())
    resp_keys = sorted(k for k in keys if RESPONSE_KEY_RE.search(k))

    blob = "\n".join(char_lines + title_lines + source_lines)
    resp_val_hits = sorted(set(RESPONSE_VAL_RE.findall(blob)))
    ici_hits = sorted({m.group(0).lower() for m in ICI_RE.finditer(blob)})

    return {
        "gse": gse,
        "n_gsm_meta": n_gsm,
        "platforms_meta": ";".join(platforms),
        "organisms": ";".join(organisms),
        "n_char_fields": len(keys),
        "response_keys": ";".join(resp_keys),
        "response_value_tokens": ";".join(resp_val_hits),
        "ici_tokens": ";".join(ici_hits),
        "has_response_annot": int(bool(resp_keys) or len(resp_val_hits) >= 2),
        "has_ici_annot": int(bool(ici_hits)),
        "lung_sample_mention": int(bool(LUNG_RE.search(blob))),
        "cellline_mention": int(bool(CELLLINE_RE.search(blob))),
        "fetch_error": int(text.startswith("__FETCH_ERROR__")),
    }


def main() -> int:
    ensure_dirs()
    hits = pd.read_csv(RESULTS_DIR / "geo_search_hits.tsv", sep="\t")

    to_screen = hits[
        (hits.is_array_type == 1)
        | ((hits.lung_mention == 1) & (hits.ici_mention == 1))
        | hits.found_by.fillna("").str.contains("seed|panel")
    ].copy()
    accs = sorted(to_screen.gse.unique(), key=lambda g: int(g[3:]))
    print(f"[screen] {len(accs)} series", flush=True)

    ev_dir = RESULTS_DIR / "series_screen_evidence"
    ev_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for i, gse in enumerate(accs, 1):
        text = brief_metadata(gse)
        row = screen_one(gse, text)
        rows.append(row)
        if row["has_response_annot"] and row["has_ici_annot"]:
            keep = [
                l
                for l in text.splitlines()
                if l.startswith(("!Sample_characteristics_ch", "!Sample_source_name_ch", "!Sample_title", "!Sample_platform_id"))
            ]
            (ev_dir / f"{gse}_evidence.txt").write_text("\n".join(keep) + "\n")
        if i % 25 == 0:
            print(f"  .. {i}/{len(accs)}", flush=True)

    scr = pd.DataFrame(rows).merge(
        hits[["gse", "n_samples", "pdat", "gds_type", "is_array_type", "platforms", "platform_titles", "pubmed", "title"]],
        on="gse",
        how="left",
    )
    scr["priority"] = (
        scr.has_response_annot * 4
        + scr.has_ici_annot * 4
        + scr.lung_sample_mention * 2
        + scr.is_array_type.fillna(0).astype(int) * 2
        - scr.cellline_mention
    )
    scr = scr.sort_values(["priority", "gse"], ascending=[False, True])
    out = RESULTS_DIR / "series_screen.tsv"
    scr.to_csv(out, sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)
    print(f"[write] {out} ({len(scr)} rows)")

    top = scr[(scr.has_response_annot == 1) & (scr.has_ici_annot == 1) & (scr.lung_sample_mention == 1)]
    print(f"\n[shortlist] response + ICI + lung annotated at sample level: {len(top)}")
    for _, r in top.iterrows():
        print(
            f"  {r.gse}  n={r.n_gsm_meta:>4}  {str(r.pdat)[:10]}  {str(r.gds_type)[:40]:40s}  "
            f"{str(r.platforms_meta)[:24]:24s} resp_keys={str(r.response_keys)[:60]}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
