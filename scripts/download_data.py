"""Backward-compatible wrapper. Prefer scripts/w200/B1_BRCA/download_data.py."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent / "w200" / "B1_BRCA"))
from download_data import main  # noqa: E402

if __name__ == "__main__":
    sys.exit(main())
