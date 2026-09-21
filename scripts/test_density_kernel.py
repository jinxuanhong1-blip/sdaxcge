#!/usr/bin/env python3
"""Kernel accuracy and planted CLDN4–CD8 association checks."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from cosmx_density_field_cldn4 import self_check


def main() -> None:
    self_check()
    print("test_density_kernel: PASS")


if __name__ == "__main__":
    main()
