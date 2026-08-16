#!/usr/bin/env python3
"""Extract Tacstd2 (Trop2) and Cldn4 by treatment across the mouse-lung ICI slice.

Loads processed data for each accession, extracts per-sample log2 expression for
the two target genes, assigns treatment groups from the authoritative GEO/
ArrayExpress metadata, and runs real statistics (Welch t-test + Mann-Whitney U,
with Benjamini-Hochberg FDR) for every treatment-vs-control contrast that has
biological replicates. Datasets with a single sample per condition (several
single-cell series and the GSE197260 time-course) are reported descriptively;
no p-values are fabricated for n=1.

Inputs live in /tmp/fable_ici_data (see download_data.sh). Outputs are written
under results/fable_mouse_ici/.
"""
import gzip
import io
import json
import tarfile
from pathlib import Path

import numpy as np
import openpyxl
from scipy import stats
from scipy.io import mmread

DATA = Path("/tmp/fable_ici_data")
RESULTS = Path(__file__).resolve().parents[2] / "results" / "fable_mouse_ici"
RESULTS.mkdir(parents=True, exist_ok=True)

GENES = ["Tacstd2", "Cldn4"]
ALIASES = {"Tacstd2": {"tacstd2", "trop2"}, "Cldn4": {"cldn4"}}
ENSEMBL = {"Tacstd2": "ENSMUSG00000051397", "Cldn4": "ENSMUSG00000041378"}

# long records: dataset, platform, gene, sample, group, is_control, value_log2, value_type
LONG = []
# author-reported precomputed stats
AUTHOR = []


def add(dataset, platform, gene, sample, group, is_control, value_log2, vtype):
    LONG.append(dict(dataset=dataset, platform=platform, gene=gene, sample=sample,
                     group=group, is_control=is_control,
                     value_log2=float(value_log2), value_type=vtype))


def match_gene(name):
    if name is None:
        return None
    n = str(name).strip().lower()
    for g, al in ALIASES.items():
        if n in al:
            return g
    return None


# ---------------------------------------------------------------------------
# GSE239485 -- bulk RNA-seq, LLC tumours, log2-normalised matrix (24 samples)
# groups from column prefix: C=Control vehicle, D=PolyI:C+aPD1, T=PolyI:C+aPD1+aC5aR1
# ---------------------------------------------------------------------------
def load_gse239485():
    ds, plat = "GSE239485", "bulk RNA-seq"
    grp = {"C": "Control vehicle", "D": "PolyIC+anti-PD1",
           "T": "PolyIC+anti-PD1+anti-C5aR1"}
    wb = openpyxl.load_workbook(DATA / "GSE239485_Processed_data.xlsx",
                                read_only=True, data_only=True)
    ws = wb["DataNorm"]
    rows = ws.iter_rows(values_only=True)
    hdr = next(rows)
    gname_i = hdr.index("gene_name")
    sample_cols = list(range(11, len(hdr)))  # C_/D_/T_ columns
    for row in rows:
        g = match_gene(row[gname_i])
        if not g:
            continue
        for ci in sample_cols:
            col = hdr[ci]
            if col is None:
                continue
            pref = str(col).split("_")[0]
            if pref not in grp:
                continue
            val = row[ci]
            if val is None:
                continue
            add(ds, plat, g, str(col), grp[pref], pref == "C", float(val),
                "log2_normalised")
    wb.close()


# ---------------------------------------------------------------------------
# GSE297630 -- Clariom S array, LLC tumours (3 control vs 3 anti-PD-1)
# per-sample log2 values + TC->symbol map both from the series family SOFT file
# ---------------------------------------------------------------------------
def load_gse297630():
    ds, plat = "GSE297630", "microarray (Clariom S)"
    path = DATA / "GSE297630_family.soft.gz"
    # 1) platform table: TC id -> gene symbol via mrna_assignment text
    tc2gene = {}
    sample_group = {}  # GSMxxx -> group
    sample_title = {}
    with gzip.open(path, "rt") as f:
        state = None
        cur_gsm = None
        for line in f:
            if line.startswith("^PLATFORM"):
                state = "platform_head"
            if line.startswith("^SAMPLE"):
                cur_gsm = line.split("=")[1].strip()
            if line.startswith("!Sample_title"):
                sample_title[cur_gsm] = line.split("=", 1)[1].strip()
            if line.startswith("!Sample_source_name_ch1"):
                src = line.split("=", 1)[1].strip().lower()
                sample_group[cur_gsm] = ("anti-PD-1" if "with anti-pd-1" in src
                                         else "Control")
            if line.startswith("!platform_table_begin"):
                phdr = next(f).rstrip("\n").split("\t")
                mi = phdr.index("mrna_assignment")
                for l in f:
                    if l.startswith("!platform_table_end"):
                        break
                    p = l.rstrip("\n").split("\t")
                    if len(p) <= mi:
                        continue
                    txt = p[mi]
                    for g in GENES:
                        # match "(Cldn4)," / "(Tacstd2)," in refseq description
                        if f"({g})" in txt:
                            tc2gene.setdefault(p[0], g)
    # 2) per-sample values from each sample table
    want = set(tc2gene)
    with gzip.open(path, "rt") as f:
        cur_gsm = None
        in_tab = False
        for line in f:
            if line.startswith("^SAMPLE"):
                cur_gsm = line.split("=")[1].strip()
            if line.startswith("!sample_table_begin"):
                in_tab = True
                next(f)  # header ID_REF VALUE
                continue
            if line.startswith("!sample_table_end"):
                in_tab = False
                continue
            if in_tab:
                p = line.rstrip("\n").split("\t")
                if p[0] in want and len(p) >= 2 and p[1] not in ("", "null"):
                    g = tc2gene[p[0]]
                    grp = sample_group.get(cur_gsm, "?")
                    add(ds, plat, g, sample_title.get(cur_gsm, cur_gsm), grp,
                        grp == "Control", float(p[1]), "log2_RMA")
    # author-reported summary (C vs P avg log2, fold change, p-val) from xlsx
    wb = openpyxl.load_workbook(DATA / "GSE297630_processed_data.xlsx",
                                read_only=True, data_only=True)
    ws = wb["Expression"]
    tc_for = {v: k for k, v in tc2gene.items()}
    rows = ws.iter_rows(values_only=True)
    hdr = None
    for row in rows:
        if row and row[0] == "ID":
            hdr = list(row)
            break
    for row in rows:
        if row[0] in tc2gene:
            g = tc2gene[row[0]]
            AUTHOR.append(dict(dataset=ds, gene=g, comparison="anti-PD-1 vs Control",
                               metric="author Clariom P vs C",
                               fold_change=row[hdr.index("Fold Change")],
                               p_value=row[hdr.index("P-val")]))
    wb.close()


# ---------------------------------------------------------------------------
# GSE241978 -- CMT167 lung line, AhR-KO vs Control (3 vs 3), normalised log block
# ---------------------------------------------------------------------------
def load_gse241978():
    ds, plat = "GSE241978", "bulk RNA-seq"
    wb = openpyxl.load_workbook(DATA / "GSE241978.xlsx", read_only=True,
                                data_only=True)
    ws = wb["2020-07-21_Sherr"]
    rows = ws.iter_rows(values_only=True)
    next(rows)  # merged title row
    hdr = next(rows)
    sym_i = hdr.index("Symbol")
    # normalised (log) per-sample block = columns 23..28
    block = {23: ("Control", True), 24: ("Control", True), 25: ("Control", True),
             26: ("AhR-KO", False), 27: ("AhR-KO", False), 28: ("AhR-KO", False)}
    p_i = 7  # authors' p (first block)
    fc_i = 10
    for row in rows:
        g = match_gene(row[sym_i])
        if not g:
            continue
        for ci, (grp, isc) in block.items():
            val = row[ci]
            if val is None:
                continue
            add(ds, plat, g, hdr[ci], grp, isc, float(val), "log2_normalised")
        AUTHOR.append(dict(dataset=ds, gene=g, comparison="AhR-KO vs Control",
                           metric="author DE", fold_change=row[fc_i],
                           p_value=row[p_i]))
    wb.close()


# ---------------------------------------------------------------------------
# GSE330658 -- bulk RNA-seq lung tumours, per-sample xlsx (TPM), 2 per group
# ---------------------------------------------------------------------------
def load_gse330658():
    ds, plat = "GSE330658", "bulk RNA-seq"

    def grp_of(fn):
        base = fn.split("_")[1] if "_" in fn else fn
        base = base.replace(".xlsx", "")
        for key in ("PTX-anti-VEGF", "anti-VEGF", "PTX", "Control"):
            if key.lower() in fn.lower():
                return key
        return base

    with tarfile.open(DATA / "GSE330658_RAW.tar") as t:
        for m in t.getmembers():
            if not m.name.endswith(".xlsx"):
                continue
            grp = grp_of(m.name)
            data = t.extractfile(m).read()
            wb = openpyxl.load_workbook(io.BytesIO(data), data_only=True)
            ws = wb[wb.sheetnames[0]]
            rows = ws.iter_rows(values_only=True)
            hdr = list(next(rows))
            ni = hdr.index("Name")
            ti = hdr.index("TPM")
            for row in rows:
                g = match_gene(row[ni])
                if not g:
                    continue
                tpm = row[ti]
                if tpm is None:
                    continue
                add(ds, plat, g, m.name.split("_")[0] + "_" + grp, grp,
                    grp == "Control", np.log2(float(tpm) + 1.0), "log2(TPM+1)")
            wb.close()


# ---------------------------------------------------------------------------
# E-MTAB-13704 -- GEMM lung, raw counts (27 samples, 6 arms). CPM -> log2.
# ---------------------------------------------------------------------------
def load_emtab13704():
    ds, plat = "E-MTAB-13704", "bulk RNA-seq"
    import csv
    # R# -> stimulus from sdrf scan names
    r2grp = {}
    with open(DATA / "E-MTAB-13704.sdrf.txt") as f:
        rd = list(csv.reader(f, delimiter="\t"))
    h = rd[0]
    scan_i = h.index("Scan Name")
    fac_i = [i for i, c in enumerate(h) if c.startswith("Factor Value")][0]
    for row in rd[1:]:
        scan = row[scan_i]
        rid = scan.split("_")[0]  # e.g. R13
        r2grp[rid] = row[fac_i].strip()
    ctrl = "vehicle"
    # read counts
    with open(DATA / "E-MTAB-13704_GEMMS_raw_counts.csv") as f:
        rd = csv.reader(f)
        cols = next(rd)[1:]  # R1..R27
        counts = {}
        libsize = np.zeros(len(cols))
        target_rows = {}
        for row in rd:
            ens = row[0].split(".")[0]
            vals = np.array([float(x) for x in row[1:]])
            libsize += vals
            for g, eid in ENSEMBL.items():
                if ens == eid:
                    target_rows[g] = vals
    for g, vals in target_rows.items():
        cpm = vals / libsize * 1e6
        log2 = np.log2(cpm + 1.0)
        for j, c in enumerate(cols):
            grp = r2grp.get(c, "?")
            add(ds, plat, g, c, grp, grp == ctrl, log2[j], "log2(CPM+1)")


# ---------------------------------------------------------------------------
# GSE197260 -- bulk TPM, 1 sample per condition (descriptive only)
# ---------------------------------------------------------------------------
def load_gse197260():
    ds, plat = "GSE197260", "bulk RNA-seq (n=1/arm)"
    labels = {
        "vehicle_d3": "vehicle", "gef_d3": "gefitinib d3", "gef_d14": "gefitinib d14",
        "gef_vehicle_d21": "gefitinib+vehicle", "gef_dc101_d21": "gefitinib+anti-VEGFR2",
        "gef_4h2_d21": "gefitinib+anti-PD-1", "gef_comb_d21": "gefitinib+anti-PD-1+anti-VEGFR2",
    }
    with gzip.open(DATA / "GSE197260_TPM.txt.gz", "rt") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        for line in f:
            p = line.rstrip("\n").split("\t")
            g = match_gene(p[0])
            if not g:
                continue
            for ci in range(1, len(hdr)):
                col = hdr[ci]
                add(ds, plat, g, col, labels.get(col, col),
                    col in ("vehicle_d3", "gef_vehicle_d21"),
                    np.log2(float(p[ci]) + 1.0), "log2(TPM+1)")


# ---------------------------------------------------------------------------
# single-cell pseudobulk helpers
# ---------------------------------------------------------------------------
def gene_index_from_tsv(path, symbol_col=1):
    idx = {}
    with gzip.open(path, "rt") as f:
        for i, line in enumerate(f):
            parts = line.rstrip("\n").split("\t")
            sym = parts[symbol_col] if len(parts) > symbol_col else parts[0]
            g = match_gene(sym)
            if g:
                idx[g] = i
    return idx


def pseudobulk_from_mtx(mtx_bytes, gene_rows):
    """Return {gene: summed_count} and total library from a gzipped mtx (genes x cells)."""
    with gzip.open(io.BytesIO(mtx_bytes), "rb") as fh:
        M = mmread(fh).tocsr()
    total = float(M.sum())
    out = {}
    for g, ri in gene_rows.items():
        out[g] = float(M[ri, :].sum())
    return out, total


def load_scrna(ds, tar_name, group_map, gene_tsv=None, per_sample_features=False,
               plat="scRNA-seq pseudobulk (n=1/arm)"):
    """Generic pseudobulk loader.

    group_map: dict prefix->group label keyed on GSM or filename token.
    gene_tsv: series-level features file (shared gene order) when not per-sample.
    per_sample_features: if True, read each sample's own *_features.tsv.gz.
    """
    shared_idx = None
    if gene_tsv is not None:
        shared_idx = gene_index_from_tsv(DATA / gene_tsv)
    with tarfile.open(DATA / tar_name) as t:
        members = t.getmembers()
        mtx = [m for m in members if "matrix.mtx" in m.name]
        for m in mtx:
            gsm = m.name.split("_")[0]
            token = m.name
            grp = None
            for k, v in group_map.items():
                if k in token:
                    grp = v
                    break
            if grp is None:
                continue
            if per_sample_features:
                prefix = m.name.rsplit("matrix.mtx", 1)[0]
                fmem = [x for x in members if x.name.startswith(prefix)
                        and "features" in x.name][0]
                gidx = gene_index_from_tsv(_tmp_tsv(t, fmem))
            else:
                gidx = shared_idx
            counts, total = pseudobulk_from_mtx(t.extractfile(m).read(), gidx)
            for g, c in counts.items():
                cpm = c / total * 1e6 if total else 0.0
                add(ds, plat, g, m.name.split(".")[0], grp, False,
                    np.log2(cpm + 1.0), "log2(pseudobulk CPM+1)")


def _tmp_tsv(t, fmem):
    """gene_index_from_tsv expects a path; write features to a temp gz path."""
    import tempfile
    data = t.extractfile(fmem).read()
    tf = tempfile.NamedTemporaryFile(suffix=".tsv.gz", delete=False)
    tf.write(data)
    tf.close()
    return tf.name


def main():
    load_gse239485()
    load_gse297630()
    load_gse241978()
    load_gse330658()
    load_emtab13704()
    load_gse197260()
    # single cell (all descriptive, 1 sample per condition)
    load_scrna("GSE133604", "GSE133604_RAW.tar",
               {"_Ctrl_": "Control", "_CtrlplusPD1_": "Control+anti-PD1",
                "_ko_": "Asf1a-KO", "_KOplusPD1_": "Asf1a-KO+anti-PD1"},
               gene_tsv="GSE133604_genes.tsv.gz")
    load_scrna("GSE129297", "GSE129297_RAW.tar",
               {"_Ctrl.": "Control", "_PD1.": "anti-PD-1", "_YKL.": "CDK7i (YKL)",
                "_combo.": "CDK7i+anti-PD-1"},
               gene_tsv="GSE129297_features.tsv.gz")
    load_scrna("GSE297632", "GSE297632_RAW.tar",
               {"_Control_": "Control", "_PD_1_treatment_": "anti-PD-1"},
               per_sample_features=True)
    load_scrna("GSE222158", "GSE222158_RAW.tar",
               {"_A_": "CD45 Control", "_B_": "CD45 anti-PD-1", "_C_": "CD45 21DC",
                "_D_": "CD45 Combo", "_AT_": "CD3 Control", "_DT_": "CD3 Combo"},
               per_sample_features=True)

    import pandas as pd
    df = pd.DataFrame(LONG)
    df.to_csv(RESULTS / "per_sample_expression.csv", index=False)
    print("per-sample rows:", len(df))
    print(df.groupby(["dataset", "gene"]).size())

    if AUTHOR:
        pd.DataFrame(AUTHOR).to_csv(RESULTS / "authors_reported_stats.csv",
                                    index=False)
    return df


if __name__ == "__main__":
    main()
