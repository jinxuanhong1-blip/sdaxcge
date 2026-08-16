#!/usr/bin/env python3
"""Cross-method concordance and cross-cohort sign table.

Reads existing results/<cohort>/scores.tsv and correlation_spearman.tsv.
Does not re-score.

    python scripts/07_concordance.py --root methods/bulk_immune
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from bulkimmune.concordance import AXES, pairwise_spearman, sign_concordance  # noqa: E402


COHORTS = ("GSE126044", "GSE135222", "TCGA_NSCLC")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, default=ROOT)
    args = p.parse_args()
    results = args.root / "results"
    demo = results / "demo"
    demo.mkdir(parents=True, exist_ok=True)

    pair_frames = []
    corr_tables = {}
    for cohort in COHORTS:
        scores = pd.read_csv(results / cohort / "scores.tsv", sep="\t", index_col=0)
        for axis, cols in AXES.items():
            pw = pairwise_spearman(scores, cols)
            pw.insert(0, "axis", axis)
            pw.insert(0, "cohort", cohort)
            pair_frames.append(pw)
        corr_path = results / cohort / "correlation_spearman.tsv"
        if corr_path.exists():
            corr_tables[cohort] = pd.read_csv(corr_path, sep="\t")

    pairs = pd.concat(pair_frames, ignore_index=True)
    pairs.to_csv(results / "method_concordance.tsv", sep="\t", index=False)
    pairs.to_csv(demo / "method_concordance.tsv", sep="\t", index=False)

    signs = sign_concordance(corr_tables)
    signs.to_csv(results / "sign_concordance.tsv", sep="\t", index=False)
    signs.to_csv(demo / "sign_concordance.tsv", sep="\t", index=False)

    # Overlap audit
    ov_frames = []
    for cohort in COHORTS:
        for kind in ("signature", "mcp"):
            path = results / cohort / f"report_{kind}_overlap.tsv"
            if path.exists():
                fr = pd.read_csv(path, sep="\t")
                fr.insert(0, "kind", kind)
                fr.insert(0, "cohort", cohort)
                ov_frames.append(fr)
    if ov_frames:
        ov = pd.concat(ov_frames, ignore_index=True)
        ov.to_csv(demo / "signature_overlap_all.tsv", sep="\t", index=False)

    # Compact JSON for the playbook
    summary = {"axes": {}, "signs": signs.to_dict(orient="records")}
    for axis in AXES:
        block = pairs[pairs["axis"] == axis]
        summary["axes"][axis] = {}
        for cohort, sub in block.groupby("cohort"):
            finite = sub.dropna(subset=["spearman_r"])
            summary["axes"][axis][cohort] = {
                "n_pairs": int(len(finite)),
                "median_r": float(finite["spearman_r"].median()) if len(finite) else None,
                "min_r": float(finite["spearman_r"].min()) if len(finite) else None,
                "max_r": float(finite["spearman_r"].max()) if len(finite) else None,
                "n_r_ge_0.5": int((finite["spearman_r"] >= 0.5).sum()),
                "n_r_le_0": int((finite["spearman_r"] <= 0).sum()),
            }
    (demo / "CONCORDANCE.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps({"n_method_pairs": int(len(pairs)), "n_sign_rows": int(len(signs))}, indent=2))
    print(signs.to_string(index=False, float_format=lambda x: f"{x:.3g}"))


if __name__ == "__main__":
    main()
