#!/usr/bin/env python3
from pathlib import Path
import runpy

HERE = Path(__file__).resolve().parent
for name in ["01_build_composition.py", "02_fit.py", "03_figures.py"]:
    print("==>", name)
    runpy.run_path(str(HERE / name), run_name="__main__")
