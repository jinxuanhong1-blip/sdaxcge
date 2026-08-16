#!/usr/bin/env python3
"""Launch the GSE166449 pembrolizumab TACSTD2/CLDN4 analysis."""

from pathlib import Path
import runpy


if __name__ == "__main__":
    runpy.run_path(
        str(Path(__file__).resolve().parents[2] / "results/w200/GSE166449/analyze.py"),
        run_name="__main__",
    )
