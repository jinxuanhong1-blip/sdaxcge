"""Run the public-data hunt (George + Chan + IMpower133 leftovers).

Does NOT contact EGA. Does NOT fabricate IMpower133 expression.
"""
from __future__ import annotations

from . import bulk_george
from . import plots_george
from . import plots_chan
from . import impower133_leftovers
from . import impower133_template


def main():
    print("=== IMpower133 leftovers (published numbers only) ===")
    impower133_leftovers.run()
    print("=== IMpower133 template guard (must refuse) ===")
    try:
        impower133_template.run()
    except impower133_template.Impower133NotPublic:
        print("template refused (expected; no EGA files)")
    print("=== George 2015 bulk + figures ===")
    plots_george.run()
    print("=== Chan 2021 atlas + figures ===")
    plots_chan.run()
    print("Done. See results/hunt_sclc_ici/REPORT.md")


if __name__ == "__main__":
    main()
