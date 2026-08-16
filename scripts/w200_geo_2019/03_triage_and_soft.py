#!/usr/bin/env python3
"""Download series-matrix / SOFT sample headers for leftover 2019 candidates.

Skip the one 2019 series already fully analyzed in the fable 2019-2021 PR
(GSE135222). For every other 2019 hit, pull the series matrix (or SOFT) and
extract sample characteristics so triage is based on deposited metadata, not
title keywords alone.
"""
import gzip
import io
import json
import re
import sys
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SEARCH = ROOT / "results" / "w200" / "GEO_2019" / "search"
TRIAGE = ROOT / "results" / "w200" / "GEO_2019" / "triage"
DL = ROOT / "results" / "w200" / "GEO_2019" / "downloads"
TRIAGE.mkdir(parents=True, exist_ok=True)
DL.mkdir(parents=True, exist_ok=True)

# Already analyzed with TACSTD2/CLDN4 vs ICI outcomes in
# cursor/fable-geo-2019-2021-c71e. Do not re-analyze.
FABLE_ALREADY_COVERED = {"GSE135222", "GSE126044", "GSE111414", "GSE182328", "GSE136961"}

FTP_TMPL = "https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{acc}/matrix/{acc}_series_matrix.txt.gz"
SOFT_TMPL = "https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{acc}/soft/{acc}_family.soft.gz"
ACC_CGI = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=gsm&form=text&view=brief"


def prefix(acc):
    return acc[:-3] + "nnn"


def _get(url, timeout=180, tries=4):
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "w200-geo-2019/1.0"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except Exception as exc:  # noqa: BLE001
            print(f"  fail {url} ({exc})", file=sys.stderr)
            time.sleep(2 ** attempt)
    return None


def parse_series_matrix(raw: bytes):
    """Return (series_fields, sample_table) from a GEO series matrix."""
    if raw[:2] == b"\x1f\x8b":
        txt = gzip.decompress(raw).decode("utf-8", "replace")
    else:
        txt = raw.decode("utf-8", "replace")
    series = {}
    sample_rows = {}
    for line in txt.splitlines():
        if line.startswith("!Series_"):
            k, _, v = line.partition("\t")
            series.setdefault(k[1:], []).append(v.strip().strip('"'))
        elif line.startswith("!Sample_") or line.startswith('"ID_REF"') or line.startswith("ID_REF"):
            parts = next(csv_split(line))
            key = parts[0].strip().strip('"').lstrip("!")
            sample_rows[key] = [p.strip().strip('"') for p in parts[1:]]
        elif line.startswith("!series_matrix_table_begin"):
            break
    n = max((len(v) for v in sample_rows.values()), default=0)
    samples = []
    ids = sample_rows.get("Sample_geo_accession", sample_rows.get("Sample_title", [""] * n))
    for i in range(n):
        rec = {"geo_accession": ids[i] if i < len(ids) else f"col{i}"}
        for k, vals in sample_rows.items():
            rec[k] = vals[i] if i < len(vals) else ""
        samples.append(rec)
    return series, samples, txt


def csv_split(line):
    import csv
    yield from csv.reader(io.StringIO(line), delimiter="\t")


def characteristics_blob(sample):
    bits = []
    for k, v in sample.items():
        if "characteristics" in k.lower() or "title" in k.lower() or "source" in k.lower() or "description" in k.lower():
            bits.append(f"{k}={v}")
    return " | ".join(bits)


def main():
    recs = json.loads((SEARCH / "candidates_metadata.json").read_text())
    leftover = [r for r in recs if r["accession"] not in FABLE_ALREADY_COVERED]
    print(f"leftover candidates: {len(leftover)} (skipped covered={sorted(FABLE_ALREADY_COVERED & {r['accession'] for r in recs})})")

    all_samples = []
    matrix_notes = []
    for r in leftover:
        acc = r["accession"]
        url = FTP_TMPL.format(prefix=prefix(acc), acc=acc)
        print(f"== {acc} {url}")
        raw = _get(url)
        note = {"accession": acc, "matrix_url": url, "ok": False, "n_samples_parsed": 0,
                "has_expression_table": False, "n_char_keys": 0, "error": ""}
        if raw is None:
            # fallback: SOFT family (can be large) — try brief GSM via acc.cgi
            print(f"  no matrix, trying acc.cgi GSM brief")
            raw2 = _get(ACC_CGI.format(acc=acc))
            if raw2:
                (DL / acc).mkdir(exist_ok=True)
                (DL / acc / f"{acc}_gsm_brief.txt").write_bytes(raw2)
                note["error"] = "no series matrix; saved GSM brief"
            else:
                note["error"] = "no series matrix and no GSM brief"
            matrix_notes.append(note)
            time.sleep(0.3)
            continue
        (DL / acc).mkdir(exist_ok=True)
        (DL / acc / f"{acc}_series_matrix.txt.gz").write_bytes(raw)
        try:
            series, samples, txt = parse_series_matrix(raw)
        except Exception as exc:  # noqa: BLE001
            note["error"] = f"parse failed: {exc}"
            matrix_notes.append(note)
            continue
        note["ok"] = True
        note["n_samples_parsed"] = len(samples)
        note["has_expression_table"] = "!series_matrix_table_begin" in txt
        char_keys = [k for k in (samples[0] if samples else {}) if "characteristics" in k.lower()]
        note["n_char_keys"] = len(char_keys)
        note["sample_keys"] = sorted((samples[0] if samples else {}).keys())
        matrix_notes.append(note)
        for s in samples:
            all_samples.append({
                "accession": acc,
                "gsm": s.get("geo_accession") or s.get("Sample_geo_accession", ""),
                "title": s.get("Sample_title", s.get("title", "")),
                "source": s.get("Sample_source_name_ch1", s.get("source_name_ch1", "")),
                "organism": s.get("Sample_organism_ch1", ""),
                "characteristics": characteristics_blob(s),
            })
        print(f"  samples={len(samples)} expr_table={note['has_expression_table']} char_keys={char_keys}")
        time.sleep(0.25)

    (TRIAGE / "matrix_notes.json").write_text(json.dumps(matrix_notes, indent=2))
    (TRIAGE / "leftover_samples.json").write_text(json.dumps(all_samples, indent=2))
    print(f"wrote {TRIAGE / 'matrix_notes.json'} and leftover_samples.json ({len(all_samples)} samples)")


if __name__ == "__main__":
    main()
