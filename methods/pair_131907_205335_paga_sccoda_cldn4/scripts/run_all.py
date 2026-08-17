#!/usr/bin/env python3
"""scCODA first (tables), then download / extract / PAGA if requested."""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = Path(__file__).resolve().parent


def run(cmd: list[str]) -> None:
    print("+", " ".join(cmd), flush=True)
    subprocess.check_call(cmd)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--data", type=Path, default=Path("/tmp/pair_131907_205335"))
    p.add_argument("--skip-download", action="store_true")
    p.add_argument("--sccoda-only", action="store_true")
    args = p.parse_args()

    run([sys.executable, str(SCRIPTS / "run_sccoda.py"), "--outdir", str(ROOT / "results")])
    if args.sccoda_only:
        run([sys.executable, str(SCRIPTS / "write_finding.py")])
        return
    if not args.skip_download:
        run([sys.executable, str(SCRIPTS / "download.py"), "--out", str(args.data)])
    run(
        [
            sys.executable,
            str(SCRIPTS / "extract.py"),
            "--data",
            str(args.data),
            "--out",
            str(args.data / "pair_subsample.h5ad"),
        ]
    )
    run(
        [
            sys.executable,
            str(SCRIPTS / "run_paga.py"),
            "--input",
            str(args.data / "pair_subsample.h5ad"),
            "--outdir",
            str(ROOT / "results"),
        ]
    )
    run([sys.executable, str(SCRIPTS / "write_finding.py")])


if __name__ == "__main__":
    main()
