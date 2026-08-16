"""Helpers to parse GEO series_matrix metadata for the ICI bulk analysis."""
import gzip
import re


def _split_line(line):
    # series_matrix rows are tab separated; values are wrapped in double quotes
    parts = line.rstrip("\n").split("\t")
    key = parts[0].strip().strip('"')
    vals = [p.strip().strip('"') for p in parts[1:]]
    return key, vals


def parse_series_matrix(path):
    """Return dict with 'samples' = list of per-sample dicts.

    Collects !Sample_title, !Sample_description, and every
    !Sample_characteristics_ch1 line (key: value) into aligned columns.
    """
    titles = None
    descriptions = None
    geo_acc = None
    char_rows = []  # list of (label, [values])
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!Sample_title"):
                _, titles = _split_line(line)
            elif line.startswith("!Sample_geo_accession"):
                _, geo_acc = _split_line(line)
            elif line.startswith("!Sample_description"):
                _, descriptions = _split_line(line)
            elif line.startswith("!Sample_characteristics_ch1"):
                _, vals = _split_line(line)
                char_rows.append(vals)
            elif line.startswith("!series_matrix_table_begin"):
                break

    n = len(titles) if titles else 0
    samples = []
    for i in range(n):
        d = {"title": titles[i]}
        if geo_acc and i < len(geo_acc):
            d["geo_accession"] = geo_acc[i]
        if descriptions and i < len(descriptions):
            d["description"] = descriptions[i]
        for row in char_rows:
            if i < len(row) and row[i]:
                cell = row[i]
                if ":" in cell:
                    k, v = cell.split(":", 1)
                    d[k.strip().lower()] = v.strip()
        samples.append(d)
    return {"samples": samples}
