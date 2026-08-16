#!/usr/bin/env python3
"""End-to-end public-only scCODA / DM grid."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def run(script: str, extra: list[str] | None = None) -> None:
    cmd = [sys.executable, str(ROOT / "scripts" / script), "--root", str(ROOT)]
    if extra:
        cmd.extend(extra)
    print("+", " ".join(cmd), flush=True)
    subprocess.check_call(cmd)


def main() -> None:
    run("00_download.py")
    run("01_extract.py", ["--which", "all"])
    run("02_build_composition.py")
    run("03_fit_grid.py")
    run("04_figures.py")


if __name__ == "__main__":
    main()
