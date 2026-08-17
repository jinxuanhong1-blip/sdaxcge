#!/usr/bin/env python3
"""Download public processed GEO files for the leftover Cldn4 slice."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from common import DATA

BASE = "https://ftp.ncbi.nlm.nih.gov"

FILES = [
    # GSE297632 LLC ± aPD-1
    ("GSE297632", "geo/samples/GSM8995nnn/GSM8995911/suppl/GSM8995911_Control_GEX_barcodes.tsv.gz"),
    ("GSE297632", "geo/samples/GSM8995nnn/GSM8995911/suppl/GSM8995911_Control_GEX_features.tsv.gz"),
    ("GSE297632", "geo/samples/GSM8995nnn/GSM8995911/suppl/GSM8995911_Control_GEX_matrix.mtx.gz"),
    ("GSE297632", "geo/samples/GSM8995nnn/GSM8995912/suppl/GSM8995912_PD_1_treatment_GEX_barcodes.tsv.gz"),
    ("GSE297632", "geo/samples/GSM8995nnn/GSM8995912/suppl/GSM8995912_PD_1_treatment_GEX_features.tsv.gz"),
    ("GSE297632", "geo/samples/GSM8995nnn/GSM8995912/suppl/GSM8995912_PD_1_treatment_GEX_matrix.mtx.gz"),
    # GSE285606 KP 344SQ WT vs PD1R1, both IgG
    ("GSE285606", "geo/samples/GSM8705nnn/GSM8705589/suppl/GSM8705589_140P-IgG-T1.barcodes.tsv.gz"),
    ("GSE285606", "geo/samples/GSM8705nnn/GSM8705589/suppl/GSM8705589_140P-IgG-T1.features.tsv.gz"),
    ("GSE285606", "geo/samples/GSM8705nnn/GSM8705589/suppl/GSM8705589_140P-IgG-T1.matrix.mtx.gz"),
    ("GSE285606", "geo/samples/GSM8705nnn/GSM8705590/suppl/GSM8705590_140P-IgG-T2.barcodes.tsv.gz"),
    ("GSE285606", "geo/samples/GSM8705nnn/GSM8705590/suppl/GSM8705590_140P-IgG-T2.features.tsv.gz"),
    ("GSE285606", "geo/samples/GSM8705nnn/GSM8705590/suppl/GSM8705590_140P-IgG-T2.matrix.mtx.gz"),
    ("GSE285606", "geo/samples/GSM8705nnn/GSM8705591/suppl/GSM8705591_WT-IgG-T1.barcodes.tsv.gz"),
    ("GSE285606", "geo/samples/GSM8705nnn/GSM8705591/suppl/GSM8705591_WT-IgG-T1.features.tsv.gz"),
    ("GSE285606", "geo/samples/GSM8705nnn/GSM8705591/suppl/GSM8705591_WT-IgG-T1.matrix.mtx.gz"),
    ("GSE285606", "geo/samples/GSM8705nnn/GSM8705592/suppl/GSM8705592_WT-IgG-T2.barcodes.tsv.gz"),
    ("GSE285606", "geo/samples/GSM8705nnn/GSM8705592/suppl/GSM8705592_WT-IgG-T2.features.tsv.gz"),
    ("GSE285606", "geo/samples/GSM8705nnn/GSM8705592/suppl/GSM8705592_WT-IgG-T2.matrix.mtx.gz"),
    # GSE261890 spatial processed tars (not RAW 1.8 GB, not RDS)
    ("GSE261890", "geo/samples/GSM8153nnn/GSM8153570/suppl/GSM8153570_Sen.Res.Veh.Spatial.tar.gz"),
    ("GSE261890", "geo/samples/GSM8153nnn/GSM8153571/suppl/GSM8153571_Sen.Res.ICI.Spatial.tar.gz"),
]


def curl(url: str, dest: Path):
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        print("exists", dest.name, dest.stat().st_size, flush=True)
        return
    tmp = dest.with_suffix(dest.suffix + ".partial")
    cmd = ["curl", "-fL", "--retry", "5", "--retry-delay", "4", "-o", str(tmp), url]
    print("GET", url, flush=True)
    subprocess.check_call(cmd)
    tmp.rename(dest)
    print("OK", dest.name, dest.stat().st_size, flush=True)


def main():
    for gse, rel in FILES:
        dest = DATA / gse / Path(rel).name
        curl(f"{BASE}/{rel}", dest)
    return 0


if __name__ == "__main__":
    sys.exit(main())
