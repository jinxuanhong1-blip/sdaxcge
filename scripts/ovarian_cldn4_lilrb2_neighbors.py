#!/usr/bin/env python3
"""CLDN4-only spot-neighbor proxy on the six ovarian Visium sections in Fig. 3H.

This is not SpaMAP and not the paper's enrichment score. Liu et al. pooled all
CLDNs. Here each GSE211956 section is scored from its own spot matrix.

Pre-specified:
- in-tissue spots only
- index spots: EPCAM > 0 (spot gate, not a cell-type call)
- CLDN4-high vs CLDN4-low: median split of CLDN4 on those spots; if the median
  is 0, high means CLDN4 > 0 and low means CLDN4 == 0
- LILRB2+ spot: LILRB2 UMI > 0 (no macrophage label is in the GEO matrix)
- outcome: mean count of in-tissue hex ring-1 neighbors that are LILRB2+
- self is not a neighbor
- section is the unit
"""

from __future__ import annotations

import csv
import gzip
import zipfile
from pathlib import Path

SAMPLES = [
    ("GSM6506110", "SP1", "Poor"),
    ("GSM6506111", "SP2", "Good"),
    ("GSM6506112", "SP3", "Good"),
    ("GSM6506114", "SP5", "Good"),
    ("GSM6506116", "SP7", "Poor"),
    ("GSM6506117", "SP8", "Poor"),
]

CLDN_PANEL = [
    "CLDN1", "CLDN2", "CLDN3", "CLDN4", "CLDN5", "CLDN6", "CLDN7", "CLDN8",
    "CLDN9", "CLDN10", "CLDN11", "CLDN12", "CLDN14", "CLDN15", "CLDN16",
    "CLDN17", "CLDN18", "CLDN19", "CLDN20", "CLDN22",
]


def hex_neighbors(row: int, col: int) -> list[tuple[int, int]]:
    # 10x Visium array: odd rows are shifted.
    if row % 2 == 0:
        deltas = [(0, -1), (0, 1), (-1, 0), (-1, -1), (1, 0), (1, -1)]
    else:
        deltas = [(0, -1), (0, 1), (-1, 0), (-1, 1), (1, 0), (1, 1)]
    return [(row + dr, col + dc) for dr, dc in deltas]


def load_gene_columns(mtx_dir: Path, gsm: str, symbols: list[str]) -> tuple[list[str], dict[str, list[float]]]:
    feat = mtx_dir / f"{gsm}_features.tsv.gz"
    bars = mtx_dir / f"{gsm}_barcodes.tsv.gz"
    mtx = mtx_dir / f"{gsm}_matrix.mtx.gz"
    wanted = {}
    with gzip.open(feat, "rt") as handle:
        for i, line in enumerate(handle, start=1):
            parts = line.rstrip("\n").split("\t")
            symbol = parts[1] if len(parts) > 1 else parts[0]
            if symbol in symbols:
                wanted[i] = symbol
    barcodes = [line.strip() for line in gzip.open(bars, "rt")]
    cols = {symbol: [0.0] * len(barcodes) for symbol in symbols}
    missing = [s for s in symbols if s not in set(wanted.values())]
    if missing:
        raise SystemExit(f"{gsm} missing genes: {missing}")
    with gzip.open(mtx, "rt") as handle:
        line = handle.readline()
        while line.startswith("%"):
            line = handle.readline()
        n_genes, n_cells, _nnz = map(int, line.split())
        if n_cells != len(barcodes):
            raise SystemExit(f"{gsm} barcode/matrix mismatch {len(barcodes)} vs {n_cells}")
        for line in handle:
            r_s, c_s, v_s = line.split()
            symbol = wanted.get(int(r_s))
            if symbol is None:
                continue
            cols[symbol][int(c_s) - 1] = float(v_s)
    return barcodes, cols


def load_positions(mtx_dir: Path, gsm: str) -> dict[str, tuple[int, int, int]]:
    zpath = mtx_dir / f"{gsm}_spatial.zip"
    with zipfile.ZipFile(zpath) as zf:
        raw = zf.read("spatial/tissue_positions_list.csv").decode()
    out = {}
    for line in raw.splitlines():
        if not line or line.startswith("barcode"):
            continue
        barcode, in_tissue, row, col, _y, _x = line.split(",")[:6]
        out[barcode] = (int(in_tissue), int(row), int(col))
    return out


def median(vals: list[float]) -> float:
    s = sorted(vals)
    n = len(s)
    if n == 0:
        return float("nan")
    mid = n // 2
    if n % 2:
        return s[mid]
    return 0.5 * (s[mid - 1] + s[mid])


def sign_flip_p(deltas: list[float]) -> float:
    """Two-sided exact sign-flip p for the sum of paired deltas. n<=12."""
    n = len(deltas)
    obs = abs(sum(deltas))
    extreme = 0
    total = 1 << n
    for mask in range(total):
        signed = 0.0
        for i, d in enumerate(deltas):
            signed += d if (mask >> i) & 1 else -d
        if abs(signed) >= obs - 1e-12:
            extreme += 1
    return extreme / total


def score_sample(mtx_dir: Path, gsm: str, outcome: str) -> dict:
    symbols = ["EPCAM", "CLDN4", "LILRB2", *CLDN_PANEL]
    # unique preserve order
    seen = []
    for s in symbols:
        if s not in seen:
            seen.append(s)
    barcodes, cols = load_gene_columns(mtx_dir, gsm, seen)
    pos = load_positions(mtx_dir, gsm)
    tissue_idx = []
    for i, bc in enumerate(barcodes):
        info = pos.get(bc)
        if info and info[0] == 1:
            tissue_idx.append(i)
    by_rc = {}
    for i in tissue_idx:
        _intissue, row, col = pos[barcodes[i]]
        by_rc[(row, col)] = i
    lilrb2_pos = {i for i in tissue_idx if cols["LILRB2"][i] > 0}
    epi = [i for i in tissue_idx if cols["EPCAM"][i] > 0]
    cldn4_vals = [cols["CLDN4"][i] for i in epi]
    med = median(cldn4_vals)
    if med == 0:
        high = [i for i in epi if cols["CLDN4"][i] > 0]
        low = [i for i in epi if cols["CLDN4"][i] == 0]
        split = "CLDN4>0 vs CLDN4==0 (median of EPCAM>0 spots is 0)"
    else:
        high = [i for i in epi if cols["CLDN4"][i] > med]
        low = [i for i in epi if cols["CLDN4"][i] <= med]
        split = "CLDN4 > median vs <= median among EPCAM>0 spots"

    def mean_lilrb2_neighbors(idxs: list[int]) -> float:
        if not idxs:
            return float("nan")
        total = 0
        for i in idxs:
            _t, row, col = pos[barcodes[i]]
            n = 0
            for rc in hex_neighbors(row, col):
                j = by_rc.get(rc)
                if j is not None and j in lilrb2_pos:
                    n += 1
            total += n
        return total / len(idxs)

    pooled_vals = []
    for i in epi:
        pooled_vals.append(sum(cols[g][i] for g in CLDN_PANEL))
    pooled_med = median(pooled_vals)
    if pooled_med == 0:
        phigh = [i for i in epi if sum(cols[g][i] for g in CLDN_PANEL) > 0]
        plow = [i for i in epi if sum(cols[g][i] for g in CLDN_PANEL) == 0]
    else:
        phigh = [i for i in epi if sum(cols[g][i] for g in CLDN_PANEL) > pooled_med]
        plow = [i for i in epi if sum(cols[g][i] for g in CLDN_PANEL) <= pooled_med]

    high_mean = mean_lilrb2_neighbors(high)
    low_mean = mean_lilrb2_neighbors(low)
    return {
        "gsm": gsm,
        "fig3h_outcome": outcome,
        "n_tissue_spots": len(tissue_idx),
        "n_epcam_pos": len(epi),
        "n_lilrb2_pos_tissue": len(lilrb2_pos),
        "cldn4_pos_among_epcam": sum(v > 0 for v in cldn4_vals),
        "cldn4_median_epcam_pos": med,
        "split": split,
        "n_cldn4_high": len(high),
        "n_cldn4_low": len(low),
        "mean_lilrb2_neighbors_cldn4_high": f"{high_mean:.6f}",
        "mean_lilrb2_neighbors_cldn4_low": f"{low_mean:.6f}",
        "delta_high_minus_low": f"{(high_mean - low_mean):.6f}",
        "mean_lilrb2_neighbors_pooled_high": f"{mean_lilrb2_neighbors(phigh):.6f}",
        "mean_lilrb2_neighbors_pooled_low": f"{mean_lilrb2_neighbors(plow):.6f}",
        "n_pooled_high": len(phigh),
        "n_pooled_low": len(plow),
    }


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--mtx-dir", type=Path, default=Path("/tmp/ovarian"))
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("results/gse316655_cldn_lilrb/ovarian_cldn4_lilrb2_neighbors.tsv"),
    )
    args = parser.parse_args()
    rows = [score_sample(args.mtx_dir, gsm, outcome) for gsm, _tag, outcome in SAMPLES]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    deltas = [float(r["delta_high_minus_low"]) for r in rows]
    n_neg = sum(d < 0 for d in deltas)
    n_pos = sum(d > 0 for d in deltas)
    print(f"sections {len(rows)} delta high-low neg {n_neg} pos {n_pos}")
    print(f"sign-flip p {sign_flip_p(deltas):.4f}")
    for r in rows:
        print(
            r["gsm"],
            r["fig3h_outcome"],
            "epi",
            r["n_epcam_pos"],
            "cldn4+",
            r["cldn4_pos_among_epcam"],
            "LILRB2+",
            r["n_lilrb2_pos_tissue"],
            "d",
            r["delta_high_minus_low"],
        )


if __name__ == "__main__":
    main()
