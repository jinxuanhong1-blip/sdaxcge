#!/usr/bin/env python3
"""
Fetch GEO series-matrix *metadata* for the relevant candidates and extract
per-sample characteristics so we can tell which series carry usable, per-patient
ICI outcome annotations (response / survival / treatment).

For each relevant GSE:
  * download all GSE*_series_matrix.txt.gz files (metadata header only),
  * parse !Sample_* lines,
  * record the distinct characteristic keys and a sample of values,
  * detect whether the matrix embeds an expression table (typical for arrays)
    vs. metadata-only (typical for RNA-seq, where counts live in suppl/).

Outputs:
  results/fable_geo_2022_2023/series_matrix_meta/<GSE>.json   (per-series)
  results/fable_geo_2022_2023/series_characteristics_summary.csv
"""
import gzip
import io
import json
import re
import sys
import time
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUTDIR = ROOT / "results" / "fable_geo_2022_2023"
METADIR = OUTDIR / "series_matrix_meta"
METADIR.mkdir(parents=True, exist_ok=True)
FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"


def matrix_prefix(acc):
    num = re.match(r"GSE(\d+)", acc).group(1)
    return "GSE" + (num[:-3] + "nnn" if len(num) > 3 else "nnn")


def list_matrix_files(acc):
    url = f"{FTP}/{matrix_prefix(acc)}/{acc}/matrix/"
    try:
        with urllib.request.urlopen(url, timeout=60) as r:
            html = r.read().decode("utf-8", "replace")
    except Exception as e:  # noqa: BLE001
        return [], str(e)
    files = re.findall(r'href="([^"]*series_matrix\.txt\.gz)"', html)
    return sorted(set(files)), ""


def fetch_matrix_meta(acc, fname):
    url = f"{FTP}/{matrix_prefix(acc)}/{acc}/matrix/{fname}"
    with urllib.request.urlopen(url, timeout=120) as r:
        raw = r.read()
    meta_lines = []
    has_table = False
    n_table_rows = 0
    with gzip.open(io.BytesIO(raw), "rt", errors="replace") as fh:
        in_table = False
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                has_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                in_table = False
                continue
            if in_table:
                n_table_rows += 1
                continue
            if line.startswith("!"):
                meta_lines.append(line.rstrip("\n"))
    return meta_lines, has_table, n_table_rows


def parse_meta(meta_lines):
    d = {}
    char_keys = {}
    for ln in meta_lines:
        parts = ln.split("\t")
        key = parts[0].lstrip("!")
        vals = [p.strip().strip('"') for p in parts[1:]]
        if key.startswith("Sample_characteristics_ch"):
            for v in vals:
                if ":" in v:
                    k2, val = v.split(":", 1)
                    char_keys.setdefault(k2.strip().lower(), set()).add(val.strip())
        elif key in ("Sample_title", "Sample_source_name_ch1",
                     "Sample_treatment_protocol_ch1", "Sample_description"):
            d.setdefault(key, [])
            for v in vals:
                if v and v not in d[key]:
                    d[key].append(v)
    return d, char_keys


OUTCOME_KEY_RE = re.compile(
    r"(response|responder|recist|benefit|resist|sensitiv|pfs|"
    r"progression|overall survival|\bos\b|survival|outcome|efficac|"
    r"dcb|\borr\b|relapse|event|status|treatment|therapy|drug|regimen|"
    r"pd-?1|pd-?l1|checkpoint|immunother)", re.I)


def main():
    ver = pd.read_csv(OUTDIR / "geo_verified.csv")
    rel = ver[ver.relevant].sort_values("accession")
    summary_rows = []
    for _, r in rel.iterrows():
        acc = r["accession"]
        files, err = list_matrix_files(acc)
        series_meta = {"accession": acc, "matrix_files": files,
                       "n_samples": int(r["n_samples"]) if str(r["n_samples"]).isdigit() else None,
                       "characteristics": {}, "has_expression_table": False,
                       "n_table_rows": 0, "error": err}
        all_char = {}
        for fn in files:
            try:
                meta_lines, has_table, n_rows = fetch_matrix_meta(acc, fn)
            except Exception as e:  # noqa: BLE001
                series_meta["error"] += f" | {fn}: {e}"
                continue
            _, char_keys = parse_meta(meta_lines)
            for k, vs in char_keys.items():
                all_char.setdefault(k, set()).update(vs)
            series_meta["has_expression_table"] = series_meta["has_expression_table"] or has_table
            series_meta["n_table_rows"] = max(series_meta["n_table_rows"], n_rows)
            time.sleep(0.15)
        # serialize characteristics (cap values)
        series_meta["characteristics"] = {
            k: sorted(list(v))[:25] for k, v in sorted(all_char.items())}
        (METADIR / f"{acc}.json").write_text(
            json.dumps(series_meta, indent=2, ensure_ascii=False))

        outcome_keys = [k for k in all_char if OUTCOME_KEY_RE.search(k)]
        summary_rows.append({
            "accession": acc,
            "n_samples": series_meta["n_samples"],
            "n_matrix_files": len(files),
            "has_expression_table": series_meta["has_expression_table"],
            "n_table_rows": series_meta["n_table_rows"],
            "n_char_keys": len(all_char),
            "char_keys": ";".join(sorted(all_char)),
            "outcome_keys": ";".join(sorted(outcome_keys)),
            "has_outcome_annotation": bool(outcome_keys),
            "error": series_meta["error"],
        })
        print(f"{acc}: char_keys={len(all_char)} expr_table={series_meta['has_expression_table']} "
              f"table_rows={series_meta['n_table_rows']} outcome_keys={outcome_keys}",
              file=sys.stderr)

    pd.DataFrame(summary_rows).to_csv(
        OUTDIR / "series_characteristics_summary.csv", index=False)
    print("\nDone.", file=sys.stderr)


if __name__ == "__main__":
    main()
