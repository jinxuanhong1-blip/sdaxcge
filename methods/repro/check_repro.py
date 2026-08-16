#!/usr/bin/env python3
"""Check a project against methods/repro/playbook.md using only stdlib."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

MAX_FASTQ_BYTES = 2_000_000_000
REQUIRED_CATALOG_COLUMNS = (
    "sample_id",
    "accession",
    "database",
    "source_url",
    "verified_at_utc",
    "verification_evidence",
    "file_name",
    "bytes",
    "decision",
    "refusal_reason",
    "sha256",
)
ACCESSION_PATTERNS = (
    re.compile(r"^(?:SRR|ERR|DRR)\d+$"),
    re.compile(r"^(?:SRS|ERS|DRS)\d+$"),
    re.compile(r"^(?:SRX|ERX|DRX)\d+$"),
    re.compile(r"^(?:SRP|ERP|DRP)\d+$"),
    re.compile(r"^(?:SAMN|SAMEA|SAMD)\d+$"),
    re.compile(r"^(?:PRJNA|PRJEB|PRJDB)\d+$"),
    re.compile(r"^GSM\d+$"),
    re.compile(r"^GSE\d+$"),
)
SHA256_RE = re.compile(r"^[0-9a-fA-F]{64}$")
UTC_RE = re.compile(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z$")
FASTQ_SUFFIXES = (".fastq", ".fq", ".fastq.gz", ".fq.gz")


class Check:
    def __init__(self) -> None:
        self.errors: list[str] = []
        self.warnings: list[str] = []

    def error(self, message: str) -> None:
        self.errors.append(message)

    def warn(self, message: str) -> None:
        self.warnings.append(message)


def resolve_from_root(root: Path, value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else root / path


def check_catalog(root: Path, catalog: Path, result: Check) -> None:
    if not catalog.is_file():
        result.error(f"missing catalog: {catalog}")
        return

    try:
        with catalog.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fieldnames = reader.fieldnames or []
            columns = set(fieldnames)
            missing = set(REQUIRED_CATALOG_COLUMNS) - columns
            if missing:
                result.error(
                    "catalog.tsv missing columns: " + ", ".join(sorted(missing))
                )
                return
            if tuple(fieldnames[: len(REQUIRED_CATALOG_COLUMNS)]) != REQUIRED_CATALOG_COLUMNS:
                result.error("catalog.tsv required columns are not in the specified order")
            rows = list(reader)
    except (OSError, UnicodeError, csv.Error) as exc:
        result.error(f"cannot parse {catalog}: {exc}")
        return

    if not rows:
        result.error("catalog.tsv has no data rows")
        return

    seen_files: set[tuple[str, str]] = set()
    for line_number, row in enumerate(rows, start=2):
        prefix = f"catalog.tsv line {line_number}"
        cell = lambda column: (row.get(column) or "").strip()
        always_required = set(REQUIRED_CATALOG_COLUMNS) - {"refusal_reason", "sha256"}
        empty = [column for column in always_required if not cell(column)]
        if empty:
            result.error(f"{prefix}: empty required fields: {', '.join(sorted(empty))}")

        accession = cell("accession")
        if accession and not any(pattern.fullmatch(accession) for pattern in ACCESSION_PATTERNS):
            result.warn(
                f"{prefix}: unsupported accession syntax {accession!r}; "
                "verify manually because syntax is not proof of existence"
            )

        source_url = cell("source_url")
        if source_url and not source_url.startswith("https://"):
            result.error(f"{prefix}: source_url must use HTTPS")

        verified_at = cell("verified_at_utc")
        if verified_at and not UTC_RE.fullmatch(verified_at):
            result.error(f"{prefix}: verified_at_utc must use YYYY-MM-DDTHH:MM:SSZ")

        evidence_value = cell("verification_evidence")
        if evidence_value:
            if Path(evidence_value).is_absolute() or ".." in Path(evidence_value).parts:
                result.error(
                    f"{prefix}: verification_evidence must be a project-relative path"
                )
            evidence = resolve_from_root(root, evidence_value)
            if not evidence.is_file():
                result.error(f"{prefix}: verification evidence not found: {evidence}")
            elif accession and accession not in evidence.read_text(
                encoding="utf-8", errors="replace"
            ):
                result.error(
                    f"{prefix}: verification evidence does not contain {accession!r}"
                )

        file_name = cell("file_name")
        file_key = (accession, file_name)
        if file_key in seen_files:
            result.error(
                f"{prefix}: duplicate accession/file_name pair "
                f"{accession!r}/{file_name!r}"
            )
        seen_files.add(file_key)

        byte_value = cell("bytes")
        byte_count: int | None = None
        try:
            byte_count = int(byte_value)
            if byte_count < 0:
                raise ValueError
        except ValueError:
            if byte_value:
                result.error(f"{prefix}: bytes must be a non-negative integer")

        decision = cell("decision").lower()
        refusal_reason = cell("refusal_reason")
        digest = cell("sha256")
        if decision not in {"accepted", "refused"}:
            result.error(f"{prefix}: decision must be accepted or refused")
        elif decision == "accepted":
            if refusal_reason not in {"", "NA"}:
                result.error(f"{prefix}: accepted row must not have a refusal reason")
            if not SHA256_RE.fullmatch(digest):
                result.error(
                    f"{prefix}: accepted row sha256 must contain exactly "
                    "64 hexadecimal characters"
                )
            if (
                byte_count is not None
                and file_name.lower().endswith(FASTQ_SUFFIXES)
                and byte_count > MAX_FASTQ_BYTES
            ):
                result.error(
                    f"{prefix}: REFUSE FASTQ {file_name!r}: "
                    f"{byte_count} bytes exceeds {MAX_FASTQ_BYTES}"
                )
        elif decision == "refused":
            if not refusal_reason or refusal_reason == "NA":
                result.error(f"{prefix}: refused row requires a refusal_reason")
            if digest not in {"", "NA"}:
                result.error(f"{prefix}: refused row sha256 must be empty or NA")


def check_seed(root: Path, result: Check) -> None:
    seed_path = root / "repro" / "seed.txt"
    if not seed_path.is_file():
        result.error(f"missing seed record: {seed_path}")
        return
    text = seed_path.read_text(encoding="utf-8").strip()
    if not re.fullmatch(r"\d+", text):
        result.error("repro/seed.txt must contain exactly one non-negative integer")


def check_environment(root: Path, result: Check) -> None:
    conda_path = root / "repro" / "conda-explicit.txt"
    pip_path = root / "repro" / "pip-freeze.txt"
    if not conda_path.is_file():
        result.error(f"missing conda export: {conda_path}")
    elif "@EXPLICIT" not in conda_path.read_text(encoding="utf-8", errors="replace"):
        result.error("repro/conda-explicit.txt must be produced by `conda list --explicit`")
    if not pip_path.is_file():
        result.error(f"missing pip freeze: {pip_path}")
    else:
        pip_text = pip_path.read_text(encoding="utf-8", errors="replace").strip()
        if not pip_text:
            result.error("repro/pip-freeze.txt is empty")
        elif pip_text.startswith("NOT_USED"):
            result.error(
                "repro/pip-freeze.txt must be produced by "
                "`python3 -m pip freeze --all`"
            )


def check_writeup(root: Path, result: Check) -> None:
    writeup = root / "WRITEUP.md"
    if not writeup.is_file():
        result.error(f"missing bilingual report: {writeup}")
        return
    text = writeup.read_text(encoding="utf-8", errors="replace")
    required_markers = (
        "## English",
        "## 中文",
        "### Sample size / 样本量",
        "### Multiple testing / 多重检验",
    )
    for marker in required_markers:
        if marker not in text:
            result.error(f"WRITEUP.md missing required heading: {marker}")
    if re.search(r"\b(?:TBD|TODO|FIXME)\b|待补|待定", text, flags=re.IGNORECASE):
        result.error("WRITEUP.md contains an unresolved placeholder")


def check_local_fastq(root: Path, max_bytes: int, result: Check) -> None:
    ignored = {".git", ".conda", ".venv", "venv", "node_modules"}
    for path in root.rglob("*"):
        if any(part in ignored for part in path.parts):
            continue
        if not path.is_file() or not path.name.lower().endswith(FASTQ_SUFFIXES):
            continue
        size = path.stat().st_size
        if size > max_bytes:
            result.error(
                f"REFUSE local FASTQ {path}: {size} bytes exceeds {max_bytes}; "
                "remove it and record the refusal"
            )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate project-wide reproducibility artifacts."
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=Path.cwd(),
        help="project root (default: current directory)",
    )
    parser.add_argument(
        "--catalog",
        default="catalog.tsv",
        help="catalog path, absolute or relative to --root",
    )
    parser.add_argument(
        "--max-fastq-bytes",
        type=int,
        default=MAX_FASTQ_BYTES,
        help="local FASTQ refusal threshold (default: 2000000000)",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = args.root.resolve()
    if args.max_fastq_bytes < 0:
        print("ERROR: --max-fastq-bytes must be non-negative", file=sys.stderr)
        return 2
    if not root.is_dir():
        print(f"ERROR: project root is not a directory: {root}", file=sys.stderr)
        return 2

    result = Check()
    catalog = resolve_from_root(root, args.catalog)
    check_catalog(root, catalog, result)
    check_seed(root, result)
    check_environment(root, result)
    check_writeup(root, result)
    check_local_fastq(root, args.max_fastq_bytes, result)

    for warning in result.warnings:
        print(f"WARNING: {warning}")
    for error in result.errors:
        print(f"ERROR: {error}")
    if result.errors:
        print(f"FAIL: {len(result.errors)} error(s), {len(result.warnings)} warning(s)")
        return 1
    print(f"PASS: reproducibility checklist ({len(result.warnings)} warning(s))")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
