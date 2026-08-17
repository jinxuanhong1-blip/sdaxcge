#!/usr/bin/env python3
"""Download, extract, and score the pair trajectory."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = Path("/tmp/pair_123902_205335_traj")


def run(cmd: list[str]) -> None:
    print("+", " ".join(cmd), flush=True)
    subprocess.check_call(cmd)


def main() -> None:
    run([sys.executable, str(HERE / "download.py"), "--outdir", str(DATA)])
    run(
        [
            sys.executable,
            str(HERE / "extract_epithelium.py"),
            "--data",
            str(DATA),
            "--out",
            str(DATA / "epithelium.h5ad"),
        ]
    )
    run(
        [
            sys.executable,
            str(HERE / "analyze.py"),
            "--input",
            str(DATA / "epithelium.h5ad"),
            "--outdir",
            str(ROOT / "results"),
            "--finding",
            str(ROOT / "FINDING.md"),
        ]
    )


if __name__ == "__main__":
    main()
