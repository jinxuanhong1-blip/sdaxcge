"""Parse GEO series-matrix files without loading the full expression table.

We stream the header (sample characteristics) and then pull only the probe
rows we care about. Clariom D matrices are ~10^5 rows; this keeps memory flat.
"""

from __future__ import annotations

import gzip
import io
import re
from pathlib import Path

import pandas as pd

from common import open_maybe_gzip


CHAR_RE = re.compile(r"^!Sample_characteristics_ch(\d+)\t(.*)$")
TITLE_RE = re.compile(r"^!Sample_title\t(.*)$")
ACC_RE = re.compile(r"^!Sample_geo_accession\t(.*)$")
SRC_RE = re.compile(r"^!Sample_source_name_ch1\t(.*)$")
PLAT_RE = re.compile(r"^!Sample_platform_id\t(.*)$")
SERIES_TITLE_RE = re.compile(r"^!Series_title\t(.*)$")
SERIES_PLAT_RE = re.compile(r"^!Series_platform_id\t(.*)$")
SERIES_TYPE_RE = re.compile(r"^!Series_type\t(.*)$")
SERIES_SUM_RE = re.compile(r"^!Series_summary\t(.*)$")
SERIES_PMID_RE = re.compile(r"^!Series_pubmed_id\t(.*)$")


def _split_quoted(line_rest: str) -> list[str]:
    """Split a GEO matrix sample row: tab-separated, typically quoted fields."""
    out = []
    for tok in line_rest.split("\t"):
        tok = tok.strip()
        if len(tok) >= 2 and tok[0] == '"' and tok[-1] == '"':
            tok = tok[1:-1]
        out.append(tok)
    return out


def parse_series_matrix_header(path: Path) -> dict:
    """Return series metadata + per-sample phenotype table (no expression)."""
    series: dict[str, str] = {"source_file": str(path)}
    samples: dict[str, list[str]] = {}
    extra_keys: list[str] = []
    n_char = 0
    with open_maybe_gzip(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("!series_matrix_table_begin"):
                break
            m = SERIES_TITLE_RE.match(line)
            if m:
                series["title"] = _split_quoted(m.group(1))[0]
                continue
            m = SERIES_PLAT_RE.match(line)
            if m:
                series["platform"] = _split_quoted(m.group(1))[0]
                continue
            m = SERIES_TYPE_RE.match(line)
            if m:
                series["type"] = _split_quoted(m.group(1))[0]
                continue
            m = SERIES_SUM_RE.match(line)
            if m:
                series.setdefault("summary", "")
                series["summary"] += _split_quoted(m.group(1))[0] + " "
                continue
            m = SERIES_PMID_RE.match(line)
            if m:
                series["pubmed"] = _split_quoted(m.group(1))[0]
                continue
            m = ACC_RE.match(line)
            if m:
                samples["geo_accession"] = _split_quoted(m.group(1))
                continue
            m = TITLE_RE.match(line)
            if m:
                samples["title"] = _split_quoted(m.group(1))
                continue
            m = SRC_RE.match(line)
            if m:
                samples["source_name"] = _split_quoted(m.group(1))
                continue
            m = PLAT_RE.match(line)
            if m:
                samples["platform_id"] = _split_quoted(m.group(1))
                continue
            m = CHAR_RE.match(line)
            if m:
                vals = _split_quoted(m.group(2))
                keys = []
                cleaned = []
                for v in vals:
                    if ":" in v:
                        k, rest = v.split(":", 1)
                        keys.append(k.strip().lower())
                        cleaned.append(rest.strip())
                    else:
                        keys.append("")
                        cleaned.append(v)
                # GEO repeats !Sample_characteristics_ch1 for each field; the
                # key is usually constant across samples.
                key = next((k for k in keys if k), f"characteristic_{n_char}")
                # disambiguate duplicate keys
                base = key
                i = 2
                while key in samples:
                    key = f"{base}_{i}"
                    i += 1
                samples[key] = cleaned
                extra_keys.append(key)
                n_char += 1
                continue
    pheno = pd.DataFrame(samples)
    if "geo_accession" in pheno.columns:
        pheno = pheno.set_index("geo_accession", drop=False)
    series["n_samples"] = int(len(pheno))
    series["pheno_fields"] = extra_keys
    return {"series": series, "pheno": pheno}


def extract_probe_rows(path: Path, probe_ids: set[str]) -> pd.DataFrame:
    """Return a probes x samples DataFrame for the requested probe IDs.

    Index = probe ID (quotes stripped). Columns = GSM accessions.
    """
    want = {p.strip().strip('"') for p in probe_ids}
    header: list[str] | None = None
    rows: dict[str, list[float]] = {}
    with open_maybe_gzip(path) as fh:
        in_table = False
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if not in_table:
                continue
            if header is None:
                header = _split_quoted(line)
                continue
            if not line:
                continue
            fields = line.split("\t")
            pid = fields[0].strip().strip('"')
            if pid not in want:
                continue
            vals = []
            for tok in fields[1:]:
                tok = tok.strip().strip('"')
                if tok in {"", "null", "NA", "NaN"}:
                    vals.append(float("nan"))
                else:
                    vals.append(float(tok))
            rows[pid] = vals
    if header is None:
        raise RuntimeError(f"no expression table in {path}")
    cols = header[1:]
    if not rows:
        return pd.DataFrame(index=pd.Index([], name="probe_id"), columns=cols, dtype=float)
    return pd.DataFrame.from_dict(rows, orient="index", columns=cols)


def list_all_probe_ids(path: Path) -> list[str]:
    ids = []
    with open_maybe_gzip(path) as fh:
        in_table = False
        header_done = False
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if not in_table:
                continue
            if not header_done:
                header_done = True
                continue
            if not line:
                continue
            ids.append(line.split("\t", 1)[0].strip().strip('"'))
    return ids


def maybe_log2(expr: pd.DataFrame) -> tuple[pd.DataFrame, bool]:
    """Log2-transform if the matrix looks like raw / linear intensity."""
    vals = expr.to_numpy(dtype=float).ravel()
    vals = vals[~pd.isna(vals)]
    if vals.size == 0:
        return expr, False
    mx = float(vals.max())
    # NanoString / linear intensities typically exceed ~30; already-logged
    # Affymetrix MAS5/RMA sit roughly in 0–16 (sometimes up to ~20).
    if mx > 30:
        import numpy as np
        out = np.log2(expr.clip(lower=0) + 1.0)
        return out, True
    return expr, False
