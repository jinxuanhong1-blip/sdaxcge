#!/usr/bin/env python3
"""GSE15212 - human colorectal (RKO) cells, TACSTD2 siRNA vs neg-ctrl siRNA (Agilent array).

Downloads per-sample quantile-normalized log2 signal (GPL4133) and the platform
probe->symbol annotation, assembles a matrix, and tests the TACSTD2-siRNA vs
negative-control (72 h) contrast.
Stats: Welch t-test + BH FDR + Mann-Whitney on log2 signal. log2FC>0 => up in siTACSTD2.
"""
import sys
import time
import urllib.request
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import analysis_common as AC

DATA = AC.DATA
GEO = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"

CASE = {  # siTACSTD2 (two siRNAs, 72h)
    "GSM379987": "TACSTD2.1.72h.4", "GSM379988": "TACSTD2.1.72h.5", "GSM379989": "TACSTD2.1.72h.6",
    "GSM379990": "TACSTD2.4.72h.4", "GSM379991": "TACSTD2.4.72h.5", "GSM379992": "TACSTD2.4.72h.6",
}
CTRL = {  # negative control siRNA, 72h
    "GSM379957": "Neg.0.72h.1", "GSM379958": "Neg.0.72h.2", "GSM379959": "Neg.0.72h.3",
    "GSM379960": "Neg.2.72h.1", "GSM379961": "Neg.2.72h.2", "GSM379962": "Neg.2.72h.3",
    "GSM379963": "Neg.2.72h.4", "GSM379964": "Neg.2.72h.5", "GSM379965": "Neg.2.72h.6",
    "GSM379966": "Neg.2.72h.7", "GSM379967": "Neg.2.72h.8", "GSM379968": "Neg.2.72h.9",
}


def fetch(url):
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=90) as r:
                return r.read().decode("utf-8", "replace")
        except Exception:
            if attempt == 3:
                raise
            time.sleep(2 ** attempt)


def gsm_values(gsm):
    cache = DATA / f"{gsm}.tsv"
    if cache.exists():
        txt = cache.read_text()
    else:
        txt = fetch(f"{GEO}?acc={gsm}&targ=self&form=text&view=full")
        cache.write_text(txt)
    rows = []
    in_tab = False
    for line in txt.splitlines():
        if line.startswith("!sample_table_begin"):
            in_tab = True; continue
        if line.startswith("!sample_table_end"):
            break
        if in_tab:
            parts = line.split("\t")
            if len(parts) >= 2 and parts[0] != "ID_REF":
                rows.append((parts[0], parts[1]))
    s = pd.Series({pid: (float(v) if v not in ("", "null", "NA") else np.nan) for pid, v in rows})
    return s


def gpl_annotation():
    cache = DATA / "GPL4133_annot.tsv"
    if not cache.exists():
        txt = fetch(f"{GEO}?acc=GPL4133&targ=self&form=text&view=full")
        # keep only the platform data table
        lines = txt.splitlines()
        try:
            b = next(i for i, l in enumerate(lines) if l.startswith("!platform_table_begin"))
            e = next(i for i, l in enumerate(lines) if l.startswith("!platform_table_end"))
        except StopIteration:
            b, e = 0, len(lines)
        cache.write_text("\n".join(lines[b + 1:e]))
    ann = pd.read_csv(cache, sep="\t", low_memory=False, dtype=str)
    col = "GENE_SYMBOL" if "GENE_SYMBOL" in ann.columns else \
          next((c for c in ann.columns if "SYMBOL" in c.upper()), None)
    ann = ann[["ID", col]].rename(columns={col: "symbol"})
    ann["ID"] = ann["ID"].astype(str)
    return ann.set_index("ID")["symbol"]


def main():
    ann = gpl_annotation()
    print("annotation probes:", len(ann), "with symbol:", ann.notna().sum())

    gsms = list(CASE) + list(CTRL)
    cols = {}
    for i, g in enumerate(gsms):
        cols[g] = gsm_values(g)
        time.sleep(0.2)
    mat_df = pd.DataFrame(cols)
    mat_df = mat_df.dropna(how="any")
    print("probe matrix:", mat_df.shape)

    symbols = mat_df.index.map(lambda pid: ann.get(pid, "")).astype(str).values
    logmat = mat_df.values.astype(float)  # already quantile-normalized log2 signal
    keep = np.array([s not in ("", "nan", "None") for s in symbols])
    symbols, logmat = symbols[keep], logmat[keep]
    symbols, logmat, _ = AC.collapse_to_gene(symbols, logmat)
    print("gene-level matrix:", logmat.shape)

    case_idx = list(range(len(CASE)))
    ctrl_idx = list(range(len(CASE), len(CASE) + len(CTRL)))
    summary, pergene, setdf = AC.run_expression_dataset(
        "GSE15212",
        "Human colorectal (RKO) cells, TACSTD2 siRNA vs neg-ctrl siRNA 72h (Agilent GPL4133)",
        symbols, logmat, case_idx, ctrl_idx,
        value_kind="quantile-normalized log2 signal",
        notes="6 siTACSTD2 (2 siRNAs x3) vs 12 neg-ctrl siRNA at 72h. log2FC>0 = up in siTACSTD2.",
    )
    print("QC perturbation:", summary["perturbation_qc"])
    print("CLDN4:", summary["cldn4"])
    print(setdf[["gene_set", "n_detected", "comp_direction", "comp_rank_biserial",
                 "comp_p", "sc_mean_log2FC", "sc_n_up", "sc_n_down", "sc_wilcoxon_p"]].to_string(index=False))


if __name__ == "__main__":
    main()
