#!/usr/bin/env python3
"""Patient-level TACSTD2/CLDN4 detection in GSE229353 CD45+ 10x MTX.

Libraries are CD45-bead sorted immune cells (Hui et al., npj Precis Oncol 2023).
This is residual/ambient/epithelial leak in an immune library, not malignant RNA.
MPR labels come from the public paper Supplementary Table 1 (not the GEO matrix).
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd

from config import GENES, GSE229353_DIR, GSE229353_RAW, GSE229353_TABLE_S1, RESULTS
from stats_lib import two_group


SAMPLES = [
    ("GSM7159183", "P01", "Chemo"),
    ("GSM7159184", "P02", "anti-PD1+Chemo"),
    ("GSM7159185", "P03", "anti-PD1+Chemo"),
    ("GSM7159186", "P04", "anti-PD1+Chemo"),
    ("GSM7159187", "P05", "anti-PD1+Chemo"),
    ("GSM7159188", "P06", "anti-PD1+Chemo"),
    ("GSM7159189", "P07", "anti-PD1+Chemo"),
]


def _read_features(path: Path) -> list[str]:
    names = []
    with gzip.open(path, "rt") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            names.append(parts[1] if len(parts) > 1 else parts[0])
    return names


def _open_mtx_text(mtx_path: Path):
    """Yield MTX lines. GEO P03 (GSM7159185) is a truncated gzip; recover what we can."""
    import zlib
    raw = Path(mtx_path).read_bytes()
    d = zlib.decompressobj(16 + zlib.MAX_WBITS)
    try:
        out = d.decompress(raw)
        complete = bool(d.eof)
    except Exception:
        d = zlib.decompressobj(16 + zlib.MAX_WBITS)
        out = d.decompress(raw)
        complete = False
    text = out.decode("ascii", errors="replace")
    if not text.endswith("\n"):
        text = text.rsplit("\n", 1)[0] + "\n"
    return text.splitlines(), complete


def _gene_counts(mtx_path: Path, gene_idx: dict[str, int]):
    """Return n_genes, n_cells, counts, lib, and a QC dict."""
    wanted = {idx: name for name, idx in gene_idx.items()}
    lines, gzip_complete = _open_mtx_text(mtx_path)
    header = None
    data_i = 0
    for i, line in enumerate(lines):
        if line.startswith("%") or not line.strip():
            continue
        header = line
        data_i = i + 1
        break
    if header is None:
        raise SystemExit(f"no MTX header in {mtx_path}")
    n_genes, n_cells, nnz = map(int, header.split())
    counts = {name: np.zeros(n_cells, dtype=np.int32) for name in gene_idx}
    lib = np.zeros(n_cells, dtype=np.int64)
    scanned = 0
    for line in lines[data_i:]:
        if not line.strip():
            continue
        parts = line.split()
        if len(parts) != 3:
            break
        g, c, v = parts
        gi = int(g) - 1
        ci = int(c) - 1
        val = int(v)
        if 0 <= ci < n_cells:
            lib[ci] += val
            if gi in wanted:
                counts[wanted[gi]][ci] = val
        scanned += 1
    qc = {
        "gzip_complete": gzip_complete,
        "declared_nnz": nnz,
        "scanned_nnz": scanned,
        "complete": bool(gzip_complete and scanned == nnz),
    }
    return n_genes, n_cells, counts, lib, qc


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    rows = []
    cell_rows = []
    for gsm, pid, geo_tx in SAMPLES:
        feat = GSE229353_RAW / f"{gsm}_{pid}_C_features.tsv.gz"
        mtx = GSE229353_RAW / f"{gsm}_{pid}_C_matrix.mtx.gz"
        bc = GSE229353_RAW / f"{gsm}_{pid}_C_barcodes.tsv.gz"
        names = _read_features(feat)
        gene_idx = {}
        for g in GENES:
            hits = [i for i, n in enumerate(names) if n == g]
            if len(hits) != 1:
                raise SystemExit(f"{pid}: {g} hits={hits}")
            gene_idx[g] = hits[0]
        n_barcodes = sum(1 for _ in gzip.open(bc, "rt"))
        n_genes, n_cells, counts, lib, qc = _gene_counts(mtx, gene_idx)
        s1 = GSE229353_TABLE_S1[pid]
        rec = {
            "dataset": "GSE229353",
            "patient": pid,
            "gsm": gsm,
            "geo_treatment": geo_tx,
            "therapy": s1["therapy"],
            "pathology": s1["pathology"],
            "rvt_pct": s1["rvt_pct"],
            "pathologic_response": s1["pathologic_response"],
            "mpr_any": "yes" if s1["pathologic_response"] == "MPR" else "no",
            "n_cells": int(n_cells),
            "n_barcodes": int(n_barcodes),
            "author_n_cells": s1["author_n_cells"],
            "sum_total_counts": int(lib.sum()),
            "mtx_complete": qc["complete"],
            "gzip_complete": qc["gzip_complete"],
            "declared_nnz": qc["declared_nnz"],
            "scanned_nnz": qc["scanned_nnz"],
            "compartment": "CD45+_bead_sorted_immune_library",
            "not_malignant_rna": True,
        }
        for g in GENES:
            vec = counts[g]
            rec[f"{g}_n_pos"] = int(np.sum(vec > 0))
            rec[f"{g}_frac_pos"] = float(np.mean(vec > 0))
            rec[f"{g}_sum"] = int(vec.sum())
            rec[f"{g}_mean_count"] = float(vec.mean())
            rec[f"{g}_mean_cpm"] = float(np.mean(np.where(lib > 0, 1e6 * vec / lib, 0.0)))
            pos = np.flatnonzero(vec > 0)
            for i in pos:
                cell_rows.append(
                    {"patient": pid, "cell_index_0based": int(i), "gene": g, "count": int(vec[i])}
                )
        rows.append(rec)

    patient = pd.DataFrame(rows)
    patient_path = RESULTS / "gse229353_patient_level.tsv"
    patient.to_csv(patient_path, sep="\t", index=False)
    pd.DataFrame(cell_rows).to_csv(RESULTS / "gse229353_positive_cells.tsv", sep="\t", index=False)

    tests = []
    usable = patient["mtx_complete"].astype(bool)
    for gene in GENES:
        for contrast, mask_note, a_mask, b_mask in [
            (
                "mpr_vs_nonMPR_complete_mtx",
                "paper Table S1; drop GEO-truncated P03 MTX",
                usable & (patient["mpr_any"] == "yes"),
                usable & (patient["mpr_any"] == "no"),
            ),
            (
                "mpr_vs_nonMPR_NAPC_complete",
                "NAPC only + complete MTX (drop P01 NAC and truncated P03)",
                usable & (patient["mpr_any"] == "yes") & (patient["therapy"] == "NAPC"),
                usable & (patient["mpr_any"] == "no") & (patient["therapy"] == "NAPC"),
            ),
        ]:
            a = patient.loc[a_mask, f"{gene}_frac_pos"].to_numpy()
            b = patient.loc[b_mask, f"{gene}_frac_pos"].to_numpy()
            t = two_group(a, b, "MPR", "non-MPR")
            t.update({"dataset": "GSE229353", "gene": gene, "metric": f"{gene}_frac_pos",
                      "contrast": contrast, "note": mask_note})
            tests.append(t)

    tests_df = pd.DataFrame(tests)
    tests_df.to_csv(RESULTS / "gse229353_stats.tsv", sep="\t", index=False)
    summary = {
        "dataset": "GSE229353",
        "compartment": "CD45+ bead-sorted immune scRNA (not malignant RNA)",
        "n_patients": int(len(patient)),
        "n_mpr": int((patient["mpr_any"] == "yes").sum()),
        "n_nonmpr": int((patient["mpr_any"] == "no").sum()),
        "n_cells_total": int(patient["n_cells"].sum()),
        "tacstd2_cells_pos": int(patient["TACSTD2_n_pos"].sum()),
        "cldn4_cells_pos": int(patient["CLDN4_n_pos"].sum()),
        "mpr_source": "Hui et al. npj Precis Oncol 2023 Supplementary Table 1 (public PDF)",
        "geo_has_mpr": False,
        "tests": tests_df.to_dict(orient="records"),
        "patients": patient.to_dict(orient="records"),
    }
    (RESULTS / "gse229353_summary.json").write_text(json.dumps(summary, indent=2))
    print(patient.to_string(index=False))
    print(tests_df[["gene", "contrast", "n_a", "n_b", "g_median_a", "g_median_b", "mw_p", "g_hedges_g"]].to_string(index=False))


if __name__ == "__main__":
    main()
