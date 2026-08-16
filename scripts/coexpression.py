"""Backward-compatible wrapper. Prefer scripts/w200/B1_BRCA/run_analysis.py."""
from pathlib import Path
import runpy
import sys

target = Path(__file__).resolve().parent / "w200" / "B1_BRCA" / "run_analysis.py"
sys.exit(runpy.run_path(str(target), run_name="__main__") or 0)
