#!/usr/bin/env python3
"""
Stage 03 - probe leftover *relevant* 2023 series.

For each leftover_relevant GSE:
  * list supplementary files and sizes
  * download the series_matrix header and parse per-sample characteristics
  * flag outcome-like fields
  * flag whether a processed expression matrix exists and is < 2 GB
  * if a small processed matrix exists, peek the first column for TACSTD2/CLDN4

Does not re-download already-analysed GSE207422 / GSE243238.

Outputs:
  results/w200/GEO_2023/leftover_probe.csv
  results/w200/GEO_2023/leftover_suppl_files.csv
  results/w200/GEO_2023/series_matrix_meta/<GSE>.json
"""
import csv
import json
import re
from pathlib import Path

from geo_common import (
    OUT, CACHE, GENE_ALIASES, list_suppl, fetch_series_matrix,
    parse_series_matrix_header, series_dir, cached_get, open_maybe_gz,
    head_size,
)

PROC_RE = re.compile(
    r"(counts?|fpkm|tpm|rpkm|cpm|expression|matrix|normalized|processed|"
    r"\.csv|\.tsv|\.txt|\.xlsx?|\.rds)",
    re.I)
RAW_RE = re.compile(r"(\.bam|\.fastq|\.fq|\.cel(\.gz)?$|\.idat|SRA)", re.I)
OUTCOME_KEY_RE = re.compile(
    r"(response|responder|recist|mpr|pcr|pfs|os\b|outcome|benefit|"
    r"residual|patholog|progression|survival|sensitive|resist)",
    re.I)


def parse_size_token(tok):
    if not tok or tok == "-":
        return None
    m = re.match(r"([0-9.]+)\s*([KMGT]?)", str(tok))
    if not m:
        return None
    val = float(m.group(1))
    return int(val * {"K": 1024, "M": 1024**2, "G": 1024**3,
                      "T": 1024**4}.get(m.group(2), 1))


def peek_genes(path, max_lines=80000):
    """Return which of TACSTD2/CLDN4 appear in the first column of a table."""
    found = {g: False for g in GENE_ALIASES}
    aliases = {a.upper(): g for g, als in GENE_ALIASES.items() for a in als}
    try:
        with open_maybe_gz(path) as fh:
            for i, line in enumerate(fh):
                if i >= max_lines:
                    break
                if i == 0:
                    continue
                token = line.split("\t", 1)[0].split(",", 1)[0].strip().strip('"')
                token = token.split(".")[0]  # ENSGxxxx.N
                if token.upper() in aliases:
                    found[aliases[token.upper()]] = True
                # also scan the whole first field for symbol-in-annotation
                up = token.upper()
                for g, als in GENE_ALIASES.items():
                    if any(a.upper() == up for a in als):
                        found[g] = True
    except Exception as e:  # noqa: BLE001
        return found, str(e)
    return found, ""


def sample_chars(meta):
    """Flatten Sample_characteristics_ch* into unique keys + example values."""
    keys = {}
    for k, blocks in meta.items():
        if not k.lower().startswith("sample_characteristics"):
            continue
        for block in blocks:
            vals = block if isinstance(block, list) else [block]
            for cell in vals:
                if ":" in str(cell):
                    ck, cv = str(cell).split(":", 1)
                    keys.setdefault(ck.strip().lower(), set()).add(cv.strip()[:80])
    return {k: sorted(v)[:8] for k, v in keys.items()}


def main():
    cls = list(csv.DictReader((OUT / "leftover_classification.csv").open()))
    targets = [r for r in cls if r["leftover_status"] == "leftover_relevant"]
    print(f"probing {len(targets)} leftover_relevant series")

    meta_dir = OUT / "series_matrix_meta"
    meta_dir.mkdir(parents=True, exist_ok=True)

    probe_rows = []
    file_rows = []

    for r in targets:
        acc = r["accession"]
        print(f"== {acc} ==")
        files = list_suppl(acc)
        proc_files = []
        total_bytes = 0
        for name, size_tok in files:
            size = parse_size_token(size_tok)
            is_proc = bool(PROC_RE.search(name)) and not RAW_RE.search(name)
            file_rows.append({
                "accession": acc, "file": name,
                "size_token": size_tok or "",
                "size_bytes": size if size is not None else "",
                "is_processed_candidate": is_proc,
            })
            if isinstance(size, int):
                total_bytes += size
            if is_proc:
                proc_files.append((name, size))

        text = fetch_series_matrix(acc)
        char_keys = {}
        outcome_keys = []
        n_gsm = 0
        has_table = False
        n_table_rows = 0
        if text:
            meta = parse_series_matrix_header(text)
            char_keys = sample_chars(meta)
            outcome_keys = [k for k in char_keys if OUTCOME_KEY_RE.search(k)]
            gsm = meta.get("Sample_geo_accession", [[]])
            if gsm:
                n_gsm = len(gsm[0]) if isinstance(gsm[0], list) else 1
            # crude table-row count
            in_table = False
            for line in text.splitlines():
                if line.startswith("!series_matrix_table_begin"):
                    in_table = True
                    has_table = True
                    continue
                if line.startswith("!series_matrix_table_end"):
                    in_table = False
                    continue
                if in_table:
                    n_table_rows += 1
            (meta_dir / f"{acc}.json").write_text(json.dumps({
                "accession": acc,
                "n_gsm": n_gsm,
                "char_keys": {k: v for k, v in char_keys.items()},
                "outcome_keys": outcome_keys,
                "has_expression_table": has_table and n_table_rows > 5,
                "n_table_rows": n_table_rows,
            }, indent=2))
        else:
            (meta_dir / f"{acc}.json").write_text(json.dumps({
                "accession": acc, "error": "no series_matrix"}, indent=2))

        # Peek genes in the smallest processed matrix under 80 MB if present.
        has_tacstd2 = ""
        has_cldn4 = ""
        peeked_file = ""
        peek_err = ""
        peekable = [(n, s) for n, s in proc_files
                    if s is not None and s < 80 * 1024 * 1024]
        peekable.sort(key=lambda x: x[1] or 0)
        if peekable:
            fname, _ = peekable[0]
            url = f"{series_dir(acc)}/suppl/{fname}"
            p = cached_get(url, f"{acc}__{fname}", max_bytes=80 * 1024 * 1024)
            if p:
                found, peek_err = peek_genes(p)
                has_tacstd2 = found.get("TACSTD2", False)
                has_cldn4 = found.get("CLDN4", False)
                peeked_file = fname
            else:
                peek_err = "download_failed_or_too_large"

        # Also peek series_matrix table if it has many rows (arrays).
        if n_table_rows > 100 and text:
            found_m = {g: False for g in GENE_ALIASES}
            aliases = {a.upper(): g for g, als in GENE_ALIASES.items() for a in als}
            in_table = False
            for line in text.splitlines():
                if line.startswith("!series_matrix_table_begin"):
                    in_table = True
                    continue
                if line.startswith("!series_matrix_table_end"):
                    break
                if not in_table:
                    continue
                token = line.split("\t", 1)[0].strip().strip('"').split(".")[0]
                if token.upper() in aliases:
                    found_m[aliases[token.upper()]] = True
            if found_m["TACSTD2"]:
                has_tacstd2 = True
            if found_m["CLDN4"]:
                has_cldn4 = True
            if not peeked_file:
                peeked_file = "series_matrix"

        under_2gb = total_bytes < 2 * 1024 ** 3
        has_outcome = bool(outcome_keys) or r["mentions_outcome"] == "True"

        # Honest usability: leftover lung ICI with open processed expression
        # containing both genes AND a per-sample ICI outcome. Stage 04 decides
        # whether the outcome is actually usable.
        usable_maybe = (
            r["is_human"] == "True" and r["is_lung"] == "True"
            and r["is_ici"] == "True"
            and (has_tacstd2 is True) and (has_cldn4 is True)
            and has_outcome
            and (bool(proc_files) or n_table_rows > 100)
        )

        probe_rows.append({
            "accession": acc,
            "title": r["title"],
            "n_samples": r["n_samples"],
            "study_category": r["study_category"],
            "n_suppl_files": len(files),
            "suppl_total_bytes": total_bytes,
            "under_2gb": under_2gb,
            "n_processed_files": len(proc_files),
            "processed_files": ";".join(n for n, _ in proc_files)[:400],
            "n_gsm_in_matrix": n_gsm,
            "has_expression_table": has_table and n_table_rows > 5,
            "n_table_rows": n_table_rows,
            "char_keys": ";".join(sorted(char_keys)),
            "outcome_keys": ";".join(outcome_keys),
            "has_outcome_annotation": bool(outcome_keys),
            "mentions_outcome_text": r["mentions_outcome"],
            "peeked_file": peeked_file,
            "has_TACSTD2": has_tacstd2,
            "has_CLDN4": has_cldn4,
            "peek_error": peek_err,
            "usable_maybe": usable_maybe,
        })

    with (OUT / "leftover_probe.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(probe_rows[0].keys()) if probe_rows else ["accession"])
        w.writeheader()
        w.writerows(probe_rows)
    with (OUT / "leftover_suppl_files.csv").open("w", newline="") as fh:
        fields = ["accession", "file", "size_token", "size_bytes",
                  "is_processed_candidate"]
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(file_rows)
    print(f"probed {len(probe_rows)}; usable_maybe="
          f"{sum(1 for x in probe_rows if x['usable_maybe'])}")


if __name__ == "__main__":
    main()
