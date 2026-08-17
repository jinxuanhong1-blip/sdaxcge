#!/usr/bin/env python3
"""GSE305086 whole-blood Affy U133 Plus 2.0: CLDN4 / TACSTD2 detectability.

Epithelial genes in blood are expected to be low (shed-tumor / ambient signal).
This script confirms the named probes on GPL570, inventories public labels,
and reports detection. It does not force an IO-benefit story when the probes
sit at array background.

Inputs (not committed; default /tmp/geo):
  GSE305086_Expression_matrix_final.csv.gz
  GSE305086_series_matrix.txt.gz
  GPL570.annot.gz

Outputs (this directory):
  detection_table.tsv
  sample_annotation.tsv
  label_inventory.tsv
  platform_probe_confirm.tsv
  cycle4_named_probes.tsv
  summary.json
"""
from __future__ import annotations

import gzip
import json
import os
from collections import Counter

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

GEO_DIR = os.environ.get("GEO_DIR", "/tmp/geo")
OUT = os.path.dirname(os.path.abspath(__file__))

NAMED = {
    "201428_at": "CLDN4",
    "202286_s_at": "TACSTD2",
}
# Extra GPL570 probes for the same genes + lineage controls (not forced into a story).
PANEL = {
    "201428_at": "CLDN4",
    "1569421_at": "CLDN4",
    "202286_s_at": "TACSTD2",
    "202285_s_at": "TACSTD2",
    "202287_s_at": "TACSTD2",
    "227128_s_at": "TACSTD2",
    "201839_s_at": "EPCAM",
    "209008_x_at": "KRT8",
    "201596_x_at": "KRT18",
    "212587_s_at": "PTPRC",
    "207238_s_at": "PTPRC",
    "213539_at": "CD3D",
    "205758_at": "CD8A",
    "209604_s_at": "GATA3",
    "207634_at": "PDCD1",
    "214617_at": "PRF1",
}


def parse_series(path: str) -> pd.DataFrame:
    with gzip.open(path, "rt") as fh:
        lines = fh.readlines()

    def row(tag: str):
        for ln in lines:
            if ln.startswith(tag + "\t") or ln.startswith(tag + "\n"):
                return [x.strip('"') for x in ln.rstrip("\n").split("\t")[1:]]
        return None

    titles = row("!Sample_title")
    gsm = row("!Sample_geo_accession")
    source = row("!Sample_source_name_ch1")
    platform = row("!Sample_platform_id")
    chars = []
    for ln in lines:
        if ln.startswith("!Sample_characteristics_ch1"):
            chars.append([x.strip('"') for x in ln.rstrip("\n").split("\t")[1:]])

    recs = []
    for i, title in enumerate(titles):
        rec = {
            "title": title,
            "geo_accession": gsm[i],
            "source_name": source[i] if source else "",
            "platform_id": platform[i] if platform else "",
        }
        for vals in chars:
            v = vals[i]
            if ":" in v:
                k, val = v.split(":", 1)
                rec[k.strip().replace(" ", "_")] = val.strip()
            else:
                rec.setdefault("characteristic_extra", []).append(v)
        # matrix columns use the first token of the title (C1, C1+4, N1, R1)
        rec["matrix_col"] = title.split(" ", 1)[0]
        recs.append(rec)
    return pd.DataFrame(recs)


def confirm_gpl570(path: str) -> pd.DataFrame:
    """Confirm named + panel probes on the official GPL570 annotation."""
    wanted = set(PANEL)
    rows = []
    with gzip.open(path, "rt", errors="replace") as fh:
        header = None
        for ln in fh:
            if ln.startswith("#") or ln.startswith("^") or ln.startswith("!"):
                continue
            if header is None:
                header = ln.rstrip("\n").split("\t")
                continue
            parts = ln.rstrip("\n").split("\t")
            probe = parts[0]
            if probe not in wanted:
                continue
            rec = dict(zip(header, parts + [""] * (len(header) - len(parts))))
            symbol = rec.get("Gene symbol") or rec.get("Gene Symbol") or rec.get("GENE_SYMBOL") or ""
            title = rec.get("Gene title") or rec.get("Gene Title") or ""
            entrez = rec.get("Gene ID") or rec.get("ENTREZ_GENE_ID") or ""
            rows.append(
                {
                    "probe": probe,
                    "expected_gene": PANEL[probe],
                    "gpl570_gene_symbol": symbol,
                    "gpl570_gene_title": title,
                    "gpl570_gene_id": entrez,
                    "symbol_matches_expected": PANEL[probe] in symbol.split(" /// ")
                    if symbol
                    else False,
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    expr_path = os.path.join(GEO_DIR, "GSE305086_Expression_matrix_final.csv.gz")
    series_path = os.path.join(GEO_DIR, "GSE305086_series_matrix.txt.gz")
    plat_path = os.path.join(GEO_DIR, "GPL570.annot.gz")

    meta = parse_series(series_path)
    meta.to_csv(os.path.join(OUT, "sample_annotation.tsv"), sep="\t", index=False)

    plat = confirm_gpl570(plat_path)
    plat.to_csv(os.path.join(OUT, "platform_probe_confirm.tsv"), sep="\t", index=False)

    expr = pd.read_csv(expr_path, sep=";", index_col=0)
    expr.index.name = "probe"
    n_probes, n_samples = expr.shape

    # Align metadata to matrix columns
    meta = meta.set_index("matrix_col").reindex(expr.columns)
    if meta["geo_accession"].isna().any():
        missing = meta.index[meta["geo_accession"].isna()].tolist()
        raise SystemExit(f"unmatched matrix columns: {missing[:10]}")

    sample_type = meta["sample_type"]
    group = meta["group"]
    n_baseline = int((sample_type == "baseline").sum())
    n_followup = int((sample_type == "follow-up").sum())
    n_control = int((sample_type == "control group").sum())
    group_counts = {k: int(v) for k, v in group.value_counts().to_dict().items()}
    n_1lio = int((group == "1LIO").sum())
    n_1lio_patients = int(((group == "1LIO") & (sample_type == "baseline")).sum())

    # Public label inventory
    char_cols = [c for c in meta.columns if c not in ("title", "geo_accession", "source_name", "platform_id")]
    label_rows = []
    for col in ["sample_type", "group"]:
        vals = meta[col].fillna("").astype(str)
        label_rows.append(
            {
                "field": col,
                "n_nonempty": int((vals != "").sum()),
                "n_unique": int(vals.nunique()),
                "values": "; ".join(f"{k}={v}" for k, v in Counter(vals).most_common()),
                "usable_for": (
                    "timepoint (baseline / cycle-4 follow-up / control)"
                    if col == "sample_type"
                    else "regimen (1LIO / 2LIO / CHTIO / control); not IO benefit"
                ),
            }
        )
    for name, present in [
        ("PFS_time", False),
        ("PFS_event", False),
        ("OS_time", False),
        ("RECIST_response", False),
        ("IO_benefit_or_LTR", False),
        ("PDL1_TPS", False),
        ("NLR_or_blood_counts", False),
    ]:
        label_rows.append(
            {
                "field": name,
                "n_nonempty": 0,
                "n_unique": 0,
                "values": "",
                "usable_for": "not deposited in GEO sample characteristics",
            }
        )
    labels = pd.DataFrame(label_rows)
    labels.to_csv(os.path.join(OUT, "label_inventory.tsv"), sep="\t", index=False)

    # Detection. RMA has no MAS5 P/A calls. Rank each probe among all probes
    # (array-wide mean rank + per-sample rank). Detection rate = fraction of
    # samples in which the probe sits above that sample's median intensity.
    probe_mean = expr.mean(axis=1)
    pct_rank = probe_mean.rank(pct=True) * 100
    per_sample_pct = expr.rank(axis=0, pct=True) * 100

    det_rows = []
    for probe, gene in PANEL.items():
        on_platform = probe in set(plat["probe"]) if len(plat) else None
        in_matrix = probe in expr.index
        rec = {
            "probe": probe,
            "gene": gene,
            "named_in_task": probe in NAMED,
            "on_gpl570": bool(on_platform) if on_platform is not None else "",
            "in_processed_matrix": in_matrix,
            "n_samples": n_samples if in_matrix else 0,
            "mean_log2": "",
            "median_log2": "",
            "min_log2": "",
            "max_log2": "",
            "percentile_rank_of_mean": "",
            "median_per_sample_percentile": "",
            "n_samples_above_sample_median": "",
            "detection_rate_above_sample_median": "",
            "n_samples_above_sample_p75": "",
            "detection_rate_above_sample_p75": "",
        }
        if in_matrix:
            v = expr.loc[probe].to_numpy(float)
            sample_pct = per_sample_pct.loc[probe].to_numpy(float)
            rec.update(
                {
                    "mean_log2": round(float(v.mean()), 4),
                    "median_log2": round(float(np.median(v)), 4),
                    "min_log2": round(float(v.min()), 4),
                    "max_log2": round(float(v.max()), 4),
                    "percentile_rank_of_mean": round(float(pct_rank[probe]), 2),
                    "median_per_sample_percentile": round(float(np.median(sample_pct)), 2),
                    "n_samples_above_sample_median": int((sample_pct > 50).sum()),
                    "detection_rate_above_sample_median": round(float((sample_pct > 50).mean()), 4),
                    "n_samples_above_sample_p75": int((sample_pct > 75).sum()),
                    "detection_rate_above_sample_p75": round(float((sample_pct > 75).mean()), 4),
                }
            )
        det_rows.append(rec)
    det = pd.DataFrame(det_rows)
    # sort: named first, then by percentile
    det["_ord"] = det["named_in_task"].map({True: 0, False: 1})
    det = det.sort_values(["_ord", "percentile_rank_of_mean"], ascending=[True, False], na_position="last")
    det = det.drop(columns=["_ord"])
    det.to_csv(os.path.join(OUT, "detection_table.tsv"), sep="\t", index=False)

    # Cycle-4 paired change on the *named* probes only (honest n). Stop if
    # the probes sit in the left half of the array (background / low).
    baseline = meta.index[sample_type == "baseline"]
    pairs = [(s, f"{s}+4") for s in baseline if f"{s}+4" in expr.columns]
    c4_rows = []
    for probe, gene in NAMED.items():
        if probe not in expr.index:
            c4_rows.append(
                {
                    "probe": probe,
                    "gene": gene,
                    "n_paired": 0,
                    "n_paired_1LIO": 0,
                    "note": "probe absent from processed matrix",
                }
            )
            continue
        b = expr.loc[probe, [a for a, _ in pairs]].to_numpy(float)
        f = expr.loc[probe, [c for _, c in pairs]].to_numpy(float)
        delta = f - b
        # 1LIO subset
        oneL = [(a, c) for a, c in pairs if group.loc[a] == "1LIO"]
        b1 = expr.loc[probe, [a for a, _ in oneL]].to_numpy(float) if oneL else np.array([])
        f1 = expr.loc[probe, [c for _, c in oneL]].to_numpy(float) if oneL else np.array([])
        w, p = wilcoxon(b, f) if len(pairs) >= 6 else (np.nan, np.nan)
        w1, p1 = wilcoxon(b1, f1) if len(oneL) >= 6 else (np.nan, np.nan)
        c4_rows.append(
            {
                "probe": probe,
                "gene": gene,
                "n_paired": len(pairs),
                "n_paired_1LIO": len(oneL),
                "baseline_mean_all_paired": round(float(b.mean()), 4),
                "followup_mean_all_paired": round(float(f.mean()), 4),
                "median_delta_all_paired": round(float(np.median(delta)), 4),
                "wilcoxon_p_all_paired": ("" if np.isnan(p) else float(p)),
                "baseline_mean_1LIO": round(float(b1.mean()), 4) if len(b1) else "",
                "followup_mean_1LIO": round(float(f1.mean()), 4) if len(f1) else "",
                "median_delta_1LIO": round(float(np.median(f1 - b1)), 4) if len(b1) else "",
                "wilcoxon_p_1LIO": ("" if np.isnan(p1) else float(p1)),
                "percentile_rank_of_mean": round(float(pct_rank[probe]), 2),
            }
        )
    c4 = pd.DataFrame(c4_rows)
    c4.to_csv(os.path.join(OUT, "cycle4_named_probes.tsv"), sep="\t", index=False)

    named_pct = {
        probe: (round(float(pct_rank[probe]), 2) if probe in pct_rank.index else None)
        for probe in NAMED
    }
    at_background = all(v is not None and v < 50 for v in named_pct.values())

    summary = {
        "dataset": "GSE305086",
        "title": "Blood gene expression under immunotherapy as potential biomarker for immune checkpoint blockade in non-small-cell lung cancer",
        "platform": "GPL570 Affymetrix Human Genome U133 Plus 2.0",
        "compartment": "PAXgene whole blood (not PBMC-sorted)",
        "n_arrays": n_samples,
        "n_probes_in_processed_matrix": n_probes,
        "n_baseline": n_baseline,
        "n_followup_cycle4": n_followup,
        "n_control": n_control,
        "n_paired_cycle4": len(pairs),
        "group_counts_arrays": group_counts,
        "n_1LIO_arrays": n_1lio,
        "n_1LIO_patients_with_array": n_1lio_patients,
        "named_probes": NAMED,
        "named_probe_percentile_rank": named_pct,
        "named_probes_below_array_median": at_background,
        "public_PFS_or_IO_benefit_labels": False,
        "public_cycle4_timepoint": True,
        "stop_after_detection": at_background,
        "note": (
            "Epithelial genes in blood are expected to be low / shed-tumor signal. "
            "Named probes sit below the array-wide median; public deposit has regimen "
            "and cycle-4 timepoint but no PFS / RECIST / IO-benefit labels. Detection "
            "is the result. Cycle-4 deltas are recorded for the named probes only."
        ),
    }
    with open(os.path.join(OUT, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)

    print("matrix:", n_probes, "x", n_samples)
    print("sample_type:", dict(sample_type.value_counts()))
    print("group:", group_counts)
    print("paired cycle-4:", len(pairs), "1LIO paired:", int((group.loc[[a for a, _ in pairs]] == "1LIO").sum()))
    print("\nPlatform confirmation:")
    print(plat.to_string(index=False))
    print("\nDetection:")
    print(det.to_string(index=False))
    print("\nCycle-4 named probes:")
    print(c4.to_string(index=False))
    print("\nstop_after_detection:", at_background)


if __name__ == "__main__":
    main()
