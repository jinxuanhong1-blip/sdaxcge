#!/usr/bin/env python3
"""
Step 4: inspect each downloaded series matrix.

For every GSE matrix we record:
  * key !Series_* header fields (title, summary, type, overall design),
  * the sample titles and characteristic lines (used to detect
    ICI-response / KD-KO experimental labels at the sample level),
  * whether an expression data table is embedded, how many data rows it has,
    the platform id, and the first few row identifiers (to infer id type).

Writes notes/geo_sweep/matrix_inspect.json.
"""
import gzip
import json
from pathlib import Path

NOTES = Path("notes/geo_sweep")
MATRIX_DIR = Path("results/geo_sweep/matrices")


def read_matrix(path):
    series = {}
    samples = {}
    n_data_rows = 0
    first_ids = []
    in_table = False
    platform_id = None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                in_table = False
                continue
            if in_table:
                if line.startswith('"ID_REF"') or line.startswith("ID_REF"):
                    continue  # header row of table
                if not line.strip():
                    continue
                n_data_rows += 1
                if len(first_ids) < 8:
                    first_ids.append(line.split("\t", 1)[0].strip('"'))
                continue
            if line.startswith("!Series_"):
                k = line.split("\t", 1)[0][1:]
                v = line.split("\t", 1)[1].strip('"') if "\t" in line else ""
                series.setdefault(k, []).append(v)
            elif line.startswith("!Sample_"):
                k = line.split("\t", 1)[0][1:]
                vals = [c.strip('"') for c in line.split("\t")[1:]]
                samples.setdefault(k, []).extend(vals) if k in samples else samples.update({k: vals})
            elif line.startswith("!Series_platform_id") or line.lower().startswith("!series_platform_id"):
                pass
    platform_id = None
    for k, v in series.items():
        if k.lower() == "series_platform_id" and v:
            platform_id = v[0]
    return series, samples, n_data_rows, first_ids, platform_id


def flatten_series(series):
    out = {}
    for k, v in series.items():
        out[k] = " | ".join(x for x in v if x)
    return out


def main():
    manifest = json.load(open(NOTES / "download_manifest.json"))
    results = []
    for m in manifest:
        acc = m["accession"]
        if m.get("status") != "ok":
            results.append({"accession": acc, "status": m.get("status")})
            continue
        for f in m["files"]:
            if "path" not in f:
                continue
            path = Path(f["path"])
            series, samples, nrows, first_ids, plat = read_matrix(path)
            sflat = flatten_series(series)
            rec = {
                "accession": acc,
                "file": f["name"],
                "platform_id": plat,
                "n_data_rows": nrows,
                "has_expression_table": nrows > 0,
                "first_row_ids": first_ids,
                "series_title": sflat.get("Series_title", ""),
                "series_type": sflat.get("Series_type", ""),
                "series_summary": sflat.get("Series_summary", ""),
                "series_overall_design": sflat.get("Series_overall_design", ""),
                "sample_titles": samples.get("Sample_title", []),
                "sample_characteristics": samples.get(
                    "Sample_characteristics_ch1", []),
                "sample_source": samples.get("Sample_source_name_ch1", []),
                "n_samples": len(samples.get("Sample_geo_accession", [])),
            }
            results.append(rec)
    (NOTES / "matrix_inspect.json").write_text(json.dumps(results, indent=2))
    n_expr = sum(1 for r in results if r.get("has_expression_table"))
    print(f"Inspected {len(results)} matrix files; "
          f"{n_expr} contain an embedded expression table.")
    for r in results:
        if r.get("has_expression_table"):
            print(f"  {r['accession']:12s} rows={r['n_data_rows']:7d} "
                  f"plat={r['platform_id']} ids={r['first_row_ids'][:3]}")


if __name__ == "__main__":
    main()
