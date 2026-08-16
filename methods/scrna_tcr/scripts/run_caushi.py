#!/usr/bin/env python3
"""Caushi public processed only: GSE176021 VDJ + annotation RDS, GSE176022 bulk TCR.

Skip EGA EGAS00001005343 / EGAD00001007728 FASTQ.
GSE176021 GEX tarballs (4.1 GB RAW.tar) are not required for clone-expansion E1.
CXCL13 RNA is not in the public annotation RDS (no expression matrix).
Tfh cluster from the RDS is a proxy secondary, not a CXCL13 RNA gate.
"""

from __future__ import annotations

import argparse
import gzip
import io
import re
import sys
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
import pyreadr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from stats_util import mw_row, paired_mw, write_json, write_tsv


def parse_series_matrix(path: Path) -> pd.DataFrame:
    titles, patients, responses, tissues, gsms = [], [], [], [], []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("!Sample_geo_accession"):
                gsms = [x.strip().strip('"') for x in line.split("\t")[1:]]
            elif line.startswith("!Sample_title"):
                titles = [x.strip().strip('"') for x in line.split("\t")[1:]]
            elif line.startswith("!Sample_source_name"):
                tissues = [x.strip().strip('"') for x in line.split("\t")[1:]]
            elif line.startswith("!Sample_characteristics") and "patientid:" in line:
                patients = [
                    x.strip().strip('"').split("patientid:")[-1].strip()
                    for x in line.split("\t")[1:]
                ]
            elif line.startswith("!Sample_characteristics") and "response status:" in line:
                responses = [
                    x.strip().strip('"').split("response status:")[-1].strip()
                    for x in line.split("\t")[1:]
                ]
    n = max(len(gsms), len(titles))
    df = pd.DataFrame(
        {
            "gsm": gsms[:n] if gsms else [""] * n,
            "title": titles[:n] if titles else [""] * n,
            "tissue": tissues[:n] if tissues else [""] * n,
            "patient": patients[:n] if patients else [""] * n,
            "response": responses[:n] if responses else [""] * n,
        }
    )
    df["response_bin"] = df["response"].str.replace(" ", "", regex=False)
    df["response_bin"] = df["response_bin"].replace({"Non-MPR": "non-MPR", "MPR": "MPR"})
    return df


def _open_double_gzip_tar(path: Path) -> tarfile.TarFile:
    with gzip.open(path, "rb") as g:
        inner = g.read()
    if inner[:2] == b"\x1f\x8b":
        payload = gzip.GzipFile(fileobj=io.BytesIO(inner)).read()
    else:
        payload = inner
    return tarfile.open(fileobj=io.BytesIO(payload), mode="r:")


def load_one_vdj(path: Path) -> pd.DataFrame:
    with _open_double_gzip_tar(path) as tf:
        member = None
        for n in tf.getnames():
            if n.endswith("filtered_contig_annotations.csv") and "/._" not in n and not Path(n).name.startswith("._"):
                member = n
                break
        if member is None:
            raise FileNotFoundError(f"no contig table in {path.name}")
        raw = tf.extractfile(member).read()
    df = pd.read_csv(io.BytesIO(raw))
    stem = path.name
    # GSM5352886_MD01-024_tumor_1.vdj.tar.gz
    m = re.match(r"(GSM\d+)_(.+)\.vdj\.tar\.gz$", stem)
    if not m:
        m = re.match(r"(GSM\d+)_(.+)\.(vdt|tdj)\.tar\.gz$", stem)
    df["gsm"] = m.group(1) if m else ""
    df["capture"] = m.group(2) if m else stem
    return df


def productive_pairs(contigs: pd.DataFrame) -> pd.DataFrame:
    need = {"barcode", "chain", "cdr3", "productive"}
    if not need.issubset(contigs.columns):
        raise ValueError(f"unexpected contig columns: {list(contigs.columns)}")
    d = contigs.copy()
    prod = d["productive"].astype(str).str.lower().isin(["true", "t", "yes"])
    d = d[prod & d["cdr3"].notna() & (d["cdr3"].astype(str) != "None")]
    tra = d[d["chain"].eq("TRA")].sort_values("umis" if "umis" in d.columns else "reads", ascending=False)
    trb = d[d["chain"].eq("TRB")].sort_values("umis" if "umis" in d.columns else "reads", ascending=False)
    tra = tra.drop_duplicates(["gsm", "barcode"])
    trb = trb.drop_duplicates(["gsm", "barcode"])
    keys = ["gsm", "capture", "barcode"]
    a = tra[keys + ["cdr3"]].rename(columns={"cdr3": "TRA_cdr3"})
    b = trb[keys + ["cdr3"]].rename(columns={"cdr3": "TRB_cdr3"})
    pair = a.merge(b, on=keys, how="inner")
    pair["clone_id"] = pair["TRA_cdr3"] + "|" + pair["TRB_cdr3"]
    return pair


def capture_endpoints(pairs: pd.DataFrame, min_tcr: int) -> pd.DataFrame:
    pairs = pairs.copy()
    pairs["clone_size"] = pairs.groupby(["gsm", "clone_id"], observed=True)["barcode"].transform("size")
    pairs["expanded"] = pairs["clone_size"] >= 2
    rows = []
    for gsm, g in pairs.groupby("gsm", observed=True):
        n = int(len(g))
        rows.append(
            {
                "gsm": gsm,
                "capture": g["capture"].iloc[0],
                "n_tcr": n,
                "n_clones": int(g["clone_id"].nunique()),
                "n_exp_clones": int(g.loc[g["expanded"], "clone_id"].nunique()),
                "frac_exp_cells": float(g["expanded"].mean()) if n else np.nan,
                "pass_depth": n >= min_tcr,
            }
        )
    return pd.DataFrame(rows)


def join_tfh(pairs: pd.DataFrame, cd3: pd.DataFrame) -> pd.DataFrame:
    """Join author CellType on barcode+patient when possible. Not a CXCL13 RNA gate."""
    cd3 = cd3.copy()
    cd3["barcode"] = cd3.index.to_series().astype(str).str.split("_").str[-1]
    cd3["imid"] = cd3["imid"].astype(str)
    cd3["is_tfh"] = cd3["CellType"].astype(str).str.contains("Tfh", case=False, na=False)
    # capture like MD01-024_tumor_1 -> patient MD01-024
    pairs = pairs.copy()
    pairs["patient"] = pairs["capture"].str.split("_").str[0]
    j = pairs.merge(
        cd3[["imid", "barcode", "CellType", "is_tfh"]],
        left_on=["patient", "barcode"],
        right_on=["imid", "barcode"],
        how="left",
    )
    return j


def load_gse176022(tar_path: Path) -> pd.DataFrame:
    rows = []
    with tarfile.open(tar_path, "r") as tf:
        for m in tf.getmembers():
            if not m.name.endswith(".txt.gz"):
                continue
            raw = tf.extractfile(m).read()
            df = pd.read_csv(io.BytesIO(gzip.decompress(raw)), sep="\t")
            stem = Path(m.name).name  # GSM5266261_MD01-004_HIV.txt.gz
            parts = stem.replace(".txt.gz", "").split("_", 1)
            gsm = parts[0]
            label = parts[1] if len(parts) > 1 else stem
            patient = label.split("_")[0]
            n = int(len(df))
            top_freq = float(df["freq"].max()) if "freq" in df.columns and n else np.nan
            rows.append(
                {
                    "gsm": gsm,
                    "patient": patient,
                    "label": label,
                    "n_clonotypes": n,
                    "top_freq": top_freq,
                    "is_control": bool(re.search(r"HIV|NoPep|No_Pep|CEF", label, re.I)),
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--vdj-dir", type=Path, required=True)
    p.add_argument("--series", type=Path, required=True)
    p.add_argument("--cd3-rds", type=Path, required=True)
    p.add_argument("--gse176022-tar", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--min-tcr", type=int, default=50)
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    series = parse_series_matrix(args.series)
    write_tsv(series, args.out / "gse176021_sample_map.tsv")

    vdj_files = sorted(args.vdj_dir.glob("*.tar.gz"))
    contig_frames = []
    failed = []
    for f in vdj_files:
        try:
            contig_frames.append(load_one_vdj(f))
        except Exception as e:
            failed.append({"file": f.name, "error": str(e)})
    if failed:
        write_tsv(pd.DataFrame(failed), args.out / "vdj_read_failures.tsv")
    contigs = pd.concat(contig_frames, ignore_index=True) if contig_frames else pd.DataFrame()
    pairs = productive_pairs(contigs) if len(contigs) else pd.DataFrame()
    cap = capture_endpoints(pairs, args.min_tcr) if len(pairs) else pd.DataFrame()
    if len(cap):
        cap = cap.merge(series, on="gsm", how="left")
        write_tsv(cap, args.out / "capture_endpoints.tsv")

        # tumor captures only, then collapse to patient (sum cells / clones)
        tumor = cap[cap["tissue"].astype(str).str.lower().eq("tumor")].copy()
        # some titles use tumor_N; source_name is the reliable tissue
        pat_rows = []
        for pid, g in tumor.groupby("patient", observed=True):
            n_tcr = int(g["n_tcr"].sum())
            n_clones = int(g["n_clones"].sum())  # not unique across captures
            pat_rows.append(
                {
                    "patient": pid,
                    "response_bin": g["response_bin"].iloc[0],
                    "n_tumor_captures": int(len(g)),
                    "n_tcr": n_tcr,
                    "n_clones_sum_captures": n_clones,
                    "frac_exp_cells_mean": float(g["frac_exp_cells"].mean()),
                    "frac_exp_cells_median": float(g["frac_exp_cells"].median()),
                    "pass_depth": n_tcr >= args.min_tcr,
                }
            )
        pat = pd.DataFrame(pat_rows)
        write_tsv(pat, args.out / "patient_tumor_endpoints.tsv")

        stats_rows = []
        work = pat[pat["pass_depth"] & pat["response_bin"].isin(["MPR", "non-MPR"])]
        stats_rows.append(
            mw_row(
                work.loc[work["response_bin"] == "MPR", "frac_exp_cells_mean"],
                work.loc[work["response_bin"] == "non-MPR", "frac_exp_cells_mean"],
                "frac_exp_cells_mean",
                "tumor_MPR_vs_non-MPR",
            )
        )
        # also capture-level (secondary; still not cells-as-replicates, but captures nest in patients)
        twork = tumor[tumor["pass_depth"] & tumor["response_bin"].isin(["MPR", "non-MPR"])]
        stats_rows.append(
            mw_row(
                twork.loc[twork["response_bin"] == "MPR", "frac_exp_cells"],
                twork.loc[twork["response_bin"] == "non-MPR", "frac_exp_cells"],
                "frac_exp_cells",
                "tumor_capture_MPR_vs_non-MPR_PSEUDOREPLICATE_NOTE",
            )
        )
        write_tsv(pd.DataFrame(stats_rows), args.out / "stats_vdj_vs_mpr.tsv")

        # Tfh proxy join
        cd3 = pyreadr.read_r(str(args.cd3_rds))[None]
        joined = join_tfh(pairs, cd3)
        joined = joined.merge(series[["gsm", "tissue", "response_bin"]], on="gsm", how="left")
        if "patient" not in joined.columns:
            joined["patient"] = joined["capture"].astype(str).str.split("_").str[0]
        # combinatorial Tfh among expanded vs not, tumor only, patient-level
        tissue = joined["tissue"].astype(str).str.lower()
        jt = joined[tissue.eq("tumor")].copy()
        if jt.empty:
            jt = joined.copy()
        jt["clone_size"] = jt.groupby(["patient", "clone_id"], observed=True)["barcode"].transform("size")
        jt["expanded"] = jt["clone_size"] >= 2
        tfh_rows = []
        for pid, g in jt.groupby("patient", observed=True):
            if g["is_tfh"].isna().all():
                continue
            gg = g.dropna(subset=["is_tfh"])
            n_exp = int(gg["expanded"].sum())
            n_non = int((~gg["expanded"]).sum())
            tfh_rows.append(
                {
                    "patient": pid,
                    "response_bin": gg["response_bin"].iloc[0],
                    "n_joined": int(len(gg)),
                    "frac_tfh": float(gg["is_tfh"].mean()) if len(gg) else np.nan,
                    "frac_tfh_among_expanded": float(gg.loc[gg["expanded"], "is_tfh"].mean()) if n_exp else np.nan,
                    "frac_tfh_among_nonexpanded": float(gg.loc[~gg["expanded"], "is_tfh"].mean()) if n_non else np.nan,
                }
            )
        tfh = pd.DataFrame(tfh_rows)
        write_tsv(tfh, args.out / "patient_tfh_proxy.tsv")
        combo = []
        if len(tfh):
            combo.append(
                paired_mw(
                    tfh["frac_tfh_among_expanded"],
                    tfh["frac_tfh_among_nonexpanded"],
                    "Tfh_proxy_expanded_vs_nonexpanded",
                )
            )
            combo.append(
                mw_row(
                    tfh.loc[tfh["response_bin"] == "MPR", "frac_tfh"],
                    tfh.loc[tfh["response_bin"] == "non-MPR", "frac_tfh"],
                    "frac_tfh",
                    "tumor_Tfh_MPR_vs_non-MPR",
                )
            )
        write_tsv(pd.DataFrame(combo), args.out / "stats_tfh_proxy.tsv")

    bulk = load_gse176022(args.gse176022_tar)
    # attach MPR from series (patient-level)
    pat_resp = series[["patient", "response_bin"]].drop_duplicates("patient")
    bulk = bulk.merge(pat_resp, on="patient", how="left")
    write_tsv(bulk, args.out / "gse176022_bulk_tcr_summary.tsv")
    # GSE176022 is MANAFEST culture TCR, not TIL scTCR and not CXCL13
    (args.out / "gse176022_status.txt").write_text(
        "GSE176022 GEO processed exists (MiXCR culture tables, 10.9 MB). "
        "It is bulk TCR from MANAFEST/virus peptide cultures, not scRNA and not CXCL13. "
        "No TACSTD2/CLDN4. Used only as an open processed leftover TCR table. "
        "EGA raw FASTQ was not downloaded.\n"
    )

    spec = {
        "dataset": "GSE176021 + GSE176022 (Caushi Nature 2021)",
        "pmid": 34290408,
        "skipped": ["EGA EGAS00001005343", "EGAD00001007728 FASTQ", "GSE176021_RAW.tar 4.1GB GEX"],
        "used": [
            "GSE176021 per-sample .vdj.tar.gz (GEO processed)",
            "GSE176021_CD3_annotations.rds (CellType only; no CXCL13 RNA)",
            "GSE176022_RAW.tar (bulk culture TCR)",
        ],
        "cxcl13_rna": "not estimable from public annotation RDS",
        "e2_tumor_tacstd2": "not estimable: CD3-sorted T cells",
        "n_vdj_files": len(vdj_files),
        "n_vdj_failed": len(failed),
        "n_paired_tcr_cells": int(len(pairs)) if len(pairs) else 0,
        "n_patients_in_rds": int(pyreadr.read_r(str(args.cd3_rds))[None]["imid"].nunique()),
        "unit": "patient (tumor captures collapsed)",
    }
    write_json(spec, args.out / "run_info.json")
    print(f"Caushi wrote {args.out}  paired_tcr={spec['n_paired_tcr_cells']}")


if __name__ == "__main__":
    main()
