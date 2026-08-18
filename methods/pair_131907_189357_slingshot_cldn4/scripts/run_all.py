#!/usr/bin/env python3
"""Download, extract, and run REAL Slingshot/PAGA on GSE131907+GSE189357."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = Path("/tmp/pair_131907_189357")
H5AD = DATA / "epithelium.h5ad"


def run(cmd: list[str]) -> None:
    print("+", " ".join(cmd), flush=True)
    subprocess.check_call(cmd)


def main() -> None:
    run([sys.executable, str(HERE / "download.py"), "--out", str(DATA)])
    run(
        [
            sys.executable,
            str(HERE / "extract_epithelium.py"),
            "--data",
            str(DATA),
            "--out",
            str(H5AD),
            "--scratch",
            str(DATA / "extract"),
        ]
    )
    run(
        [
            sys.executable,
            str(HERE / "analyze.py"),
            "--input",
            str(H5AD),
            "--outdir",
            str(ROOT / "results"),
            "--finding",
            str(ROOT / "FINDING.md"),
            "--scratch",
            str(DATA / "slingshot_scratch"),
        ]
    )


if __name__ == "__main__":
    main()
