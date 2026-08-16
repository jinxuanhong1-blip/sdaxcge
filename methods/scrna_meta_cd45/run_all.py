#!/usr/bin/env python3
"""Run the additive CD45+ immune-leak extra (GSE243013 + GSE229353 + leftover audit)."""

from analyze_gse229353 import main as g229
from analyze_gse243013 import main as g243
from audit_leftovers import main as leftover
from meta_forest import main as meta
from plot_extra import main as plot


if __name__ == "__main__":
    leftover()
    g229()
    g243()
    meta()
    plot()
