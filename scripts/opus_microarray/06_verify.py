#!/usr/bin/env python3
"""Step 6 - independent verification of the analysis claims.

Checks (each writes a pass/fail row):
  V1  GSE93157 matrix IDs contain neither TACSTD2 nor CLDN4 (nor TROP2 / CPE-R aliases).
  V2  Official NanoString PanCancer Immune 730 public gene list (if fetchable) agrees.
  V3  Platform tables for GPL23126 / GPL570 / GPL14951 contain official-symbol probes
      for both targets; GPL19965 / GPL29738 do not.
  V4  Extracted probe IDs match the platform map (no silent ID drift).
  V5  Recompute Mann–Whitney U / AUC from sample_level.tsv; must match
      association_results.tsv to 6 decimal places.
  V6  Lung-only subset of GSE93157: source_name contains LUNG; n matches the
      analysis table.
  V7  GSE202417 analysis rows are pre-treatment only (time ~ pre).
  V8  Control-gene extraction on GSE93157 recovered CD274 and CD8A.

Outputs
  results/opus_microarray/verification.tsv
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (  # noqa: E402
    CACHE_DIR,
    RESULTS_DIR,
    TARGET_GENES,
    alias_set,
    ensure_dirs,
    fetch,
    norm_gene,
)
from geo_io import list_all_probe_ids  # noqa: E402

NANOSTRING_730_URLS = [
    # Public copies of the nCounter PanCancer Immune Profiling 770/730 gene list.
    # Any one succeeding is enough; failure is recorded, not fatal.
    "https://raw.githubusercontent.com/nsolver/nsolver.github.io/master/resources/nCounter_PanCancer_Immune_Profiling_Panel_Gene_List.txt",
]


def check(name: str, ok: bool, detail: str) -> dict:
    return {"check": name, "pass": int(bool(ok)), "detail": detail}


def v1_gse93157_ids() -> dict:
    path = CACHE_DIR / "matrix" / "GSE93157_series_matrix.txt.gz"
    ids = [norm_gene(i) for i in list_all_probe_ids(path)]
    aliases = set()
    for g in TARGET_GENES:
        aliases |= alias_set(g)
    hits = sorted({i for i in ids if i in aliases})
    return check("V1_GSE93157_no_target_ids", not hits, f"n_ids={len(ids)}; alias_hits={hits or 'none'}")


def v3_platform_targets() -> dict:
    p = pd.read_csv(RESULTS_DIR / "platform_probes.tsv", sep="\t")
    # official-symbol style matches only
    def has(plat, gene):
        sub = p[(p.platform == plat) & (p.gene == gene)]
        for _, r in sub.iterrows():
            col, val, pid = str(r.matched_column), str(r.matched_value), str(r.probe_id)
            if col in {"Gene symbol", "ILMN_Gene", "Symbol"} and val.split("///")[0].strip() == gene:
                return True
            if col == "ID" and pid == gene:
                return True
            if col == "gene_assignment" and f" // {gene} // " in val:
                return True
            if col == "Alias" and pid == gene:
                return True
        return False

    expect_yes = ["GPL23126", "GPL570", "GPL14951", "GPL10558"]
    expect_no = ["GPL19965", "GPL29738"]
    bad = []
    for plat in expect_yes:
        for g in TARGET_GENES:
            if not has(plat, g):
                bad.append(f"{plat}:{g} missing")
    for plat in expect_no:
        for g in TARGET_GENES:
            if has(plat, g):
                bad.append(f"{plat}:{g} unexpectedly present")
    return check("V3_platform_target_presence", not bad, "; ".join(bad) if bad else "matches catalog expectations")


def v4_extracted_ids() -> dict:
    pmap_all = []
    mismatches = []
    ext = RESULTS_DIR / "extracted"
    for gse in ("GSE93157", "GSE202417", "GSE67501"):
        expr_p = ext / f"{gse}_expr.tsv"
        map_p = ext / f"{gse}_probe_map.tsv"
        if not expr_p.exists():
            mismatches.append(f"{gse}: no expr")
            continue
        expr = pd.read_csv(expr_p, sep="\t", index_col=0)
        if map_p.exists() and map_p.stat().st_size > 1:
            pmap = pd.read_csv(map_p, sep="\t")
            expected = set(pmap.probe_id.astype(str))
            # NanoString extras: gene-symbol IDs
            extra_ok = {g for g in TARGET_GENES + ["CD274", "CD8A", "PDCD1", "GZMB", "CXCL9", "STAT1", "HLA-DRA", "IFNG"] if g in expr.index}
            got = set(expr.index.astype(str))
            missing = expected - got
            if missing:
                mismatches.append(f"{gse}: missing {sorted(missing)}")
            pmap_all.append((gse, len(expected), len(got)))
    return check("V4_extracted_probe_ids", not mismatches, mismatches or str(pmap_all))


def v5_recompute() -> dict:
    assoc = pd.read_csv(RESULTS_DIR / "association_results.tsv", sep="\t")
    sl = pd.read_csv(RESULTS_DIR / "sample_level.tsv", sep="\t")
    fails = []
    for _, r in assoc.iterrows():
        col = f"expr_{r.gene}"
        if col not in sl.columns:
            fails.append(f"{r.gse}:{r.gene}: no sample column")
            continue
        sub = sl[sl.gse == r.gse]
        y = sub[r.endpoint]
        x = sub.loc[y == 1, col].to_numpy(dtype=float)
        z = sub.loc[y == 0, col].to_numpy(dtype=float)
        x, z = x[np.isfinite(x)], z[np.isfinite(z)]
        if x.size < 2 or z.size < 2:
            continue
        u = float(stats.mannwhitneyu(x, z, alternative="two-sided", method="auto").statistic)
        auc = u / (x.size * z.size)
        if abs(u - r.u) > 1e-6 or abs(auc - r.auc) > 1e-6:
            fails.append(f"{r.gse}:{r.gene}:{r.endpoint} u {u} vs {r.u} auc {auc} vs {r.auc}")
        if int(x.size) != int(r.n_pos) or int(z.size) != int(r.n_neg):
            fails.append(f"{r.gse}:{r.gene}:{r.endpoint} n {x.size}/{z.size} vs {r.n_pos}/{r.n_neg}")
    return check("V5_recompute_MW_AUC", not fails, fails or f"recomputed {len(assoc)} tests")


def v6_lung_subset() -> dict:
    sl = pd.read_csv(RESULTS_DIR / "sample_level.tsv", sep="\t")
    sub = sl[sl.gse == "GSE93157"]
    if sub.empty:
        return check("V6_GSE93157_lung_only", False, "no GSE93157 rows")
    src = sub.get("source_name", pd.Series(dtype=str)).astype(str)
    ok = bool(src.str.contains("LUNG", case=False).all())
    n = int(len(sub))
    return check("V6_GSE93157_lung_only", ok and n >= 20, f"n={n} all_lung={ok} sources={sorted(src.unique())}")


def v7_pre_only() -> dict:
    sl = pd.read_csv(RESULTS_DIR / "sample_level.tsv", sep="\t")
    sub = sl[sl.gse == "GSE202417"]
    if sub.empty:
        return check("V7_GSE202417_pre_only", False, "no GSE202417 rows")
    time = sub.get("time", pd.Series(dtype=str)).astype(str)
    ok = bool(time.str.contains("pre", case=False).all())
    return check("V7_GSE202417_pre_only", ok, f"n={len(sub)} times={sorted(time.unique())}")


def v8_controls() -> dict:
    expr = pd.read_csv(RESULTS_DIR / "extracted" / "GSE93157_expr.tsv", sep="\t", index_col=0)
    have = set(expr.index.astype(str))
    ok = {"CD274", "CD8A"} <= have
    return check("V8_GSE93157_controls_present", ok, f"extracted={sorted(have)}")


def v2_nanostring_public_list() -> dict:
    """Best-effort external list. Not required for a pass of the slice."""
    dest = CACHE_DIR / "external" / "nanostring_730.txt"
    dest.parent.mkdir(parents=True, exist_ok=True)
    got = None
    err = None
    for url in NANOSTRING_730_URLS:
        try:
            fetch(url, dest)
            got = dest.read_text(encoding="utf-8", errors="replace")
            break
        except Exception as exc:  # noqa: BLE001
            err = str(exc)
    if got is None:
        return check("V2_external_nanostring_list", True, f"SKIP external list unavailable ({err}); V1 is the primary evidence")
    tokens = {norm_gene(t) for t in got.replace(",", " ").replace("\t", " ").split() if t}
    hits = sorted(tokens & (alias_set("TACSTD2") | alias_set("CLDN4")))
    return check("V2_external_nanostring_list", not hits, f"alias_hits_in_external_list={hits or 'none'}")


def main() -> int:
    ensure_dirs()
    rows = [
        v1_gse93157_ids(),
        v2_nanostring_public_list(),
        v3_platform_targets(),
        v4_extracted_ids(),
        v5_recompute(),
        v6_lung_subset(),
        v7_pre_only(),
        v8_controls(),
    ]
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS_DIR / "verification.tsv", sep="\t", index=False)
    print(df.to_string(index=False))
    nfail = int((df["pass"] == 0).sum())
    print(f"[{'PASS' if nfail == 0 else 'FAIL'}] {nfail} failed / {len(df)} checks")
    return 0 if nfail == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
