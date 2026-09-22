#!/usr/bin/env python3
"""MAX EFFECT Part1 — bulk arm: maximize Tacstd2 KL>KP Δ and |Welch t|.

Public mouse bulk / array / cell-line matrices only, plus the locked GSE137244
cell-line reference. Assay type for every row is bulk (not scRNA).

Does not read or merge private 8KL matrices. Does not label LLC as KL.
Leave-one-library and n=3 floors are stored with reduced-n flags and are not
substituted for fair n≥5-vs-5 (or the largest fair n available per series).
"""

from __future__ import annotations

import gzip
import json
import math
import platform
import urllib.request
from io import StringIO
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import __version__ as scipy_version
from scipy.stats import mannwhitneyu, ttest_ind

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "tables" / "bulk_kl_kp"
FIG = ROOT / "figures"
CACHE = ROOT / "cache"
UA = "max-effect-part1-public-mouse-tacstd2/1.0"

TJ7 = ["Cldn3", "Cldn4", "Cldn6", "Cldn7", "Cdh1", "F11r", "Ocln"]
FOCUS = ["Tacstd2", "Cldn4", "Epcam", "Krt8"] + TJ7

GSE6135_TUMORS = [
    ("GSM176566", "K_268", "K", "Ad", False),
    ("GSM176567", "K_287", "K", "Ad", False),
    ("GSM176568", "K_405", "K", "Ad", False),
    ("GSM176569", "K_484", "K", "Ad", False),
    ("GSM176571", "K_498", "K", "Ad", False),
    ("GSM176573", "KP_186", "KP", "Ad", True),
    ("GSM176574", "KP_196", "KP", "Ad", True),
    ("GSM176576", "KP_197", "KP", "Ad", True),
    ("GSM176578", "KP_498", "KP", "Ad", True),
    ("GSM176580", "KP_500", "KP", "Ad", True),
    ("GSM176585", "KL_861", "KL", "Ad", True),
    ("GSM176587", "KL_113", "KL", "Sq", True),
    ("GSM176591", "KL_452", "KL", "Ad-sq", True),
    ("GSM176593", "KL_452", "KL", "Ad-sq", True),
    ("GSM176595", "KL_459", "KL", "Ad", True),
    ("GSM176596", "KL_540", "KL", "Sq", True),
    ("GSM176597", "KL_540", "KL", "Sq", True),
    ("GSM176598", "KL_547", "KL", "Ad", True),
    ("GSM176599", "KL_592", "KL", "Ad", True),
    ("GSM176600", "KL_592", "KL", "Ad", True),
]
GSE6135_PROBES = {
    "Tacstd2": ["1423323_at"],
    "Cldn4": ["1418283_at"],
    "Cldn3": ["1426332_a_at", "1434651_a_at", "1451701_x_at", "1460569_x_at"],
    "Cldn6": ["1417845_at"],
    "Cldn7": ["1448393_at"],
    "Cdh1": ["1448261_at"],
    "F11r": ["1424595_at", "1436374_x_at"],
    "Ocln": ["1448873_at"],
    "Epcam": ["1424010_at", "1454936_a_at"],
    "Krt8": ["1423756_s_at", "1426464_at"],
}


def fetch(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=300) as resp:
        dest.write_bytes(resp.read())
    return dest


def log2p1(df: pd.DataFrame) -> pd.DataFrame:
    return np.log2(df.astype(float) + 1.0)


def residual_on(y: pd.Series, x: pd.Series, cols: list[str]) -> pd.Series:
    yy = y.loc[cols].to_numpy(dtype=float)
    xx = x.loc[cols].to_numpy(dtype=float)
    X = np.column_stack([np.ones(len(cols)), xx])
    beta, *_ = np.linalg.lstsq(X, yy, rcond=None)
    return pd.Series(yy - X @ beta, index=cols, name=y.name)


def residual_on_multi(y: pd.Series, xs: list[pd.Series], cols: list[str]) -> pd.Series:
    yy = y.loc[cols].to_numpy(dtype=float)
    mats = [np.ones(len(cols))]
    for x in xs:
        mats.append(x.loc[cols].to_numpy(dtype=float))
    X = np.column_stack(mats)
    beta, *_ = np.linalg.lstsq(X, yy, rcond=None)
    return pd.Series(yy - X @ beta, index=cols, name=y.name)


def contrast(y: pd.Series, kl: list[str], kp: list[str]) -> dict:
    a = y.loc[kl].to_numpy(dtype=float)
    b = y.loc[kp].to_numpy(dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    out = {
        "n_kl": int(len(a)),
        "n_kp": int(len(b)),
        "mean_kl": float(np.mean(a)) if len(a) else float("nan"),
        "mean_kp": float(np.mean(b)) if len(b) else float("nan"),
        "delta": float(np.mean(a) - np.mean(b)) if len(a) and len(b) else float("nan"),
        "welch_t": float("nan"),
        "abs_t": float("nan"),
        "welch_p": float("nan"),
        "mw_p": float("nan"),
        "complete_separation": False,
        "kl_min_gt_kp_max": False,
    }
    if len(a) < 2 or len(b) < 2:
        return out
    if np.std(a) == 0 and np.std(b) == 0:
        out["welch_t"] = 0.0 if float(a[0]) == float(b[0]) else float("nan")
        out["abs_t"] = abs(out["welch_t"]) if np.isfinite(out["welch_t"]) else float("nan")
        out["welch_p"] = 1.0 if float(a[0]) == float(b[0]) else float("nan")
    else:
        res = ttest_ind(a, b, equal_var=False)
        out["welch_t"] = float(res.statistic)
        out["abs_t"] = abs(float(res.statistic))
        out["welch_p"] = float(res.pvalue)
    if np.allclose(a, a[0]) and np.allclose(b, b[0]) and float(a[0]) == float(b[0]):
        out["mw_p"] = 1.0
    else:
        out["mw_p"] = float(mannwhitneyu(a, b, alternative="two-sided", method="exact").pvalue)
    out["complete_separation"] = bool(a.min() > b.max() or b.min() > a.max())
    out["kl_min_gt_kp_max"] = bool(a.min() > b.max())
    return out


def mw_floor(n: int, m: int) -> float:
    if n < 1 or m < 1:
        return float("nan")
    return 2.0 / math.comb(n + m, n)


def _unique(df: pd.DataFrame) -> pd.DataFrame:
    if df.index.duplicated().any():
        return df.groupby(level=0).mean(numeric_only=True)
    return df


def add_row(
    rows: list[dict],
    *,
    accession: str,
    setting: str,
    scale: str,
    filter_name: str,
    y: pd.Series,
    kl: list[str],
    kp: list[str],
    fair_n: bool,
    fair_genotype: bool,
    note: str,
    keep_n_rule: str,
) -> None:
    st = contrast(y, kl, kp)
    rows.append(
        {
            "assay": "bulk",
            "accession": accession,
            "setting": setting,
            "scale": scale,
            "endpoint": "Tacstd2",
            "filter": filter_name,
            "fair_n": bool(fair_n),
            "fair_genotype": bool(fair_genotype),
            "keep_n_rule": keep_n_rule,
            "mw_floor": mw_floor(st["n_kl"], st["n_kp"]),
            **st,
            "note": note,
        }
    )


def load_gse137244() -> tuple[pd.DataFrame, list[str], list[str], str, str]:
    path = fetch(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/GSE137244/suppl/GSE137244_counts.fpkm.csv.gz",
        CACHE / "GSE137244_counts.fpkm.csv.gz",
    )
    df = pd.read_csv(path, index_col=0)
    kl = [c for c in df.columns if str(c).startswith("KL")]
    kp = [c for c in df.columns if str(c).startswith("B6AL10")]
    genes = [g for g in FOCUS if g in df.index]
    log = log2p1(df.loc[genes])
    return log, kl, kp, "cultured_cells", "log2(FPKM+1)"


def load_gse137396() -> tuple[pd.DataFrame, list[str], list[str], str, str]:
    path = fetch(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/GSE137396/suppl/GSE137396_Raw_genetable_GEMMnodule.txt.gz",
        CACHE / "GSE137396_Raw_genetable_GEMMnodule.txt.gz",
    )
    df = pd.read_csv(path, sep="\t", index_col=0)
    kl = [c for c in df.columns if str(c).startswith("KL")]
    kp = [c for c in df.columns if "Kras-Trp53" in str(c) or str(c).startswith("KP")]
    genes = [g for g in FOCUS if g in df.index]
    return log2p1(df.loc[genes]), kl, kp, "in_vivo_nodule", "log2(abundance+1)"


def _homer_symbol_matrix(path: Path, symbols: set[str]) -> pd.DataFrame:
    with gzip.open(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        ann_i = header.index("Annotation/Divergence")
        samples = [h.replace(" FPKM", "").replace("Aligned.out.sam", "") for h in header[ann_i + 1 :]]
        buckets: dict[str, list[list[float]]] = {}
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= ann_i:
                continue
            symbol = parts[ann_i].split("|", 1)[0]
            if symbol not in symbols:
                continue
            vals = [float(x) if x not in ("", "NA") else math.nan for x in parts[ann_i + 1 :]]
            buckets.setdefault(symbol, []).append(vals)
    mat = {sym: np.nanmean(np.vstack(rows), axis=0) for sym, rows in buckets.items()}
    return pd.DataFrame(mat, index=samples).T


def load_gse164758() -> tuple[pd.DataFrame, list[str], list[str], str, str]:
    path = fetch(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE164nnn/GSE164758/suppl/GSE164758_primary_tumors_fpkm.txt.gz",
        CACHE / "GSE164758_primary_tumors_fpkm.txt.gz",
    )
    mat = _homer_symbol_matrix(path, set(FOCUS))
    kl, kp = [], []
    for col in mat.columns:
        if col.startswith("KL") and not col.startswith("KPL") and "H457" not in col:
            kl.append(col)
        elif col.startswith("KP") and not col.startswith("KPL") and not col.startswith("KPa"):
            kp.append(col)
    return log2p1(mat), kl, kp, "primary_bulk_tumor", "log2(FPKM+1)"


def load_gse244452() -> tuple[pd.DataFrame, list[str], list[str], str, str]:
    path = fetch(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE244nnn/GSE244452/suppl/GSE244452_KPvsKL_deg_all.txt.gz",
        CACHE / "GSE244452_KPvsKL_deg_all.txt.gz",
    )
    df = pd.read_csv(path, sep="\t")
    wide = df.set_index("gene_name")[["KP1", "KP2", "KP3", "KL1", "KL2", "KL3"]]
    wide = _unique(wide)
    genes = [g for g in FOCUS if g in wide.index]
    return log2p1(wide.loc[genes]), ["KL1", "KL2", "KL3"], ["KP1", "KP2", "KP3"], "syngeneic_bulk_tumor", "log2(count+1)"


def load_gse274352() -> tuple[pd.DataFrame, list[str], list[str], str, str]:
    path = fetch(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE274nnn/GSE274352/suppl/GSE274352_normalizedcounts_genes_IFNB_vs_empty.tsv.gz",
        CACHE / "GSE274352_normalizedcounts_genes_IFNB_vs_empty.tsv.gz",
    )
    df = pd.read_csv(path, sep="\t", index_col=0)
    df = df.set_index("external_gene_name")
    df = _unique(df)
    kl = [c for c in df.columns if str(c).startswith("KL") and "empty" in str(c)]
    kp = [c for c in df.columns if str(c).startswith("KP") and "empty" in str(c)]
    genes = [g for g in FOCUS if g in df.index]
    return log2p1(df.loc[genes]), kl, kp, "cultured_cells_empty_vector", "log2(normalized_count+1)"


def load_gse274351() -> tuple[pd.DataFrame, list[str], list[str], str, str]:
    path = fetch(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE274nnn/GSE274351/suppl/GSE274351_expredata_TPM_gene.txt.gz",
        CACHE / "GSE274351_expredata_TPM_gene.txt.gz",
    )
    df = pd.read_csv(path, sep="\t", index_col=0)
    ens = {
        "Tacstd2": "ENSMUSG00000051397",
        "Cldn4": "ENSMUSG00000047501",
        "Cldn3": "ENSMUSG00000070473",
        "Cldn6": "ENSMUSG00000023906",
        "Cldn7": "ENSMUSG00000018569",
        "Cdh1": "ENSMUSG00000000303",
        "F11r": "ENSMUSG00000038235",
        "Ocln": "ENSMUSG00000021638",
        "Epcam": "ENSMUSG00000045394",
        "Krt8": "ENSMUSG00000049382",
    }
    present = {g: eid for g, eid in ens.items() if eid in df.index}
    sub = df.loc[list(present.values())].copy()
    sub.index = list(present.keys())
    kl = [c for c in sub.columns if str(c).startswith("KL")]
    kp = [c for c in sub.columns if str(c).startswith("KP")]
    return log2p1(sub), kl, kp, "lcm_early_adenoma", "log2(TPM+1)"


def _series_matrix_table(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt", errors="replace") as handle:
        lines = []
        started = False
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                started = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if started:
                lines.append(line)
    df = pd.read_csv(StringIO("".join(lines)), sep="\t", index_col=0)
    df.index = [str(i).strip('"') for i in df.index]
    df.columns = [str(c).strip('"') for c in df.columns]
    return df.apply(pd.to_numeric, errors="coerce")


def load_gse6135_cohorts() -> list[tuple[str, pd.DataFrame, list[str], list[str], str, str, bool]]:
    path = fetch(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE6nnn/GSE6135/matrix/GSE6135-GPL8321_series_matrix.txt.gz",
        CACHE / "GSE6135-GPL8321_series_matrix.txt.gz",
    )
    probes = _series_matrix_table(path)
    gene_rows = {}
    for gene, ids in GSE6135_PROBES.items():
        present = [i for i in ids if i in probes.index]
        if present:
            gene_rows[gene] = probes.loc[present].mean(axis=0)
    log_mat = pd.DataFrame(gene_rows).T
    meta = pd.DataFrame(GSE6135_TUMORS, columns=["gsm", "mouse", "arm", "histology", "primary"])
    meta = meta[meta["gsm"].isin(log_mat.columns)]

    def mouse_means(sub: pd.DataFrame) -> pd.DataFrame:
        pieces = []
        for mouse, grp in sub.groupby("mouse"):
            vals = log_mat[grp["gsm"].tolist()].mean(axis=1)
            vals.name = mouse
            pieces.append(vals)
        return pd.concat(pieces, axis=1)

    out = []
    cohorts = {
        "mouse_all_histology": (meta[meta["primary"]], True),
        "mouse_adeno_only": (meta[meta["primary"] & meta["histology"].eq("Ad")], True),
        "mouse_adeno_plus_mixed": (
            meta[meta["primary"] & meta["histology"].isin(["Ad", "Ad-sq"])],
            True,
        ),
    }
    for name, (sub, fair) in cohorts.items():
        mm = mouse_means(sub)
        kl = [m for m in mm.columns if m.startswith("KL_")]
        kp = [m for m in mm.columns if m.startswith("KP_")]
        out.append((name, mm, kl, kp, name, "array_value", fair))
    return out


def endpoint_series(log: pd.DataFrame, filter_name: str, cols: list[str]) -> pd.Series | None:
    if "Tacstd2" not in log.index:
        return None
    y = log.loc["Tacstd2"]
    if filter_name == "none":
        return y.loc[cols]
    if filter_name == "minus_Epcam":
        if "Epcam" not in log.index:
            return None
        return (y.loc[cols] - log.loc["Epcam", cols]).rename("Tacstd2")
    if filter_name == "residual_Epcam":
        if "Epcam" not in log.index:
            return None
        return residual_on(y, log.loc["Epcam"], cols)
    if filter_name == "residual_Epcam_Krt8":
        if "Epcam" not in log.index or "Krt8" not in log.index:
            return None
        return residual_on_multi(y, [log.loc["Epcam"], log.loc["Krt8"]], cols)
    if filter_name == "minus_Krt8":
        if "Krt8" not in log.index:
            return None
        return (y.loc[cols] - log.loc["Krt8", cols]).rename("Tacstd2")
    raise KeyError(filter_name)


FILTERS = ["none", "minus_Epcam", "residual_Epcam", "residual_Epcam_Krt8", "minus_Krt8"]


def sweep_matrix(
    accession: str,
    setting: str,
    scale: str,
    log: pd.DataFrame,
    kl: list[str],
    kp: list[str],
    *,
    fair_genotype: bool = True,
    note: str = "",
) -> list[dict]:
    rows: list[dict] = []
    cols = list(kl) + list(kp)
    full_n = (len(kl), len(kp))
    for filt in FILTERS:
        y = endpoint_series(log, filt, cols)
        if y is None:
            continue
        add_row(
            rows,
            accession=accession,
            setting=setting,
            scale=scale,
            filter_name=filt,
            y=y,
            kl=kl,
            kp=kp,
            fair_n=True,
            fair_genotype=fair_genotype,
            note=note or "full arm membership",
            keep_n_rule=f"full {full_n[0]} vs {full_n[1]}",
        )
    # Leave-one library / mouse: recorded with fair_n=False when it changes n.
    for drop in cols:
        kl2 = [c for c in kl if c != drop]
        kp2 = [c for c in kp if c != drop]
        if len(kl2) < 2 or len(kp2) < 2:
            continue
        y = endpoint_series(log, "none", kl2 + kp2)
        if y is None:
            continue
        reduced = (len(kl2) != len(kl)) or (len(kp2) != len(kp))
        add_row(
            rows,
            accession=accession,
            setting=setting,
            scale=scale,
            filter_name=f"leaveone:{drop}",
            y=y,
            kl=kl2,
            kp=kp2,
            fair_n=not reduced,
            fair_genotype=fair_genotype,
            note="leave-one; n may drop. Not substituted for full-n primary.",
            keep_n_rule=f"leaveone {len(kl2)} vs {len(kp2)}",
        )
    return rows


def select_max(grid: pd.DataFrame, column: str, mask: pd.Series | None = None) -> pd.Series:
    sub = grid if mask is None else grid.loc[mask]
    sub = sub[np.isfinite(sub[column])].copy()
    if sub.empty:
        return pd.Series(dtype=object)
    ranked = sub.sort_values([column, "delta", "abs_t"], ascending=[False, False, False])
    return ranked.iloc[0]


def plot_pareto(fair: pd.DataFrame, path: Path) -> None:
    if fair.empty:
        return
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    colors = {
        "GSE137244": "#C44E52",
        "GSE164758": "#4C78A8",
        "GSE6135": "#F58518",
        "GSE137396": "#54A24B",
        "GSE244452": "#B279A2",
        "GSE274352": "#8C564B",
        "GSE274351": "#E45756",
    }
    for acc, sub in fair.groupby("accession"):
        ax.scatter(
            sub["delta"],
            sub["abs_t"],
            s=42,
            c=colors.get(acc, "0.4"),
            label=acc,
            alpha=0.85,
            edgecolors="white",
            linewidths=0.4,
        )
    # Pareto front among fair rows with KL>KP (delta>0)
    pos = fair[fair["delta"] > 0].copy()
    if not pos.empty:
        pts = pos.sort_values("delta")
        front = []
        best_t = -1.0
        for _, r in pts.iloc[::-1].iterrows():
            if r["abs_t"] >= best_t:
                front.append(r)
                best_t = float(r["abs_t"])
        if front:
            fx = [r["delta"] for r in front[::-1]]
            fy = [r["abs_t"] for r in front[::-1]]
            ax.plot(fx, fy, color="0.2", lw=1.0, ls="--", label="Pareto (Δ>0)")
    ax.axvline(0, color="0.55", lw=0.8)
    ax.set_xlabel("Tacstd2 Δ = mean(KL) − mean(KP)")
    ax.set_ylabel("Welch |t|")
    ax.set_title("Bulk MAX EFFECT Part1: Tacstd2 KL>KP (fair n)")
    ax.legend(frameon=False, fontsize=8, loc="best")
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def main() -> int:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(parents=True, exist_ok=True)

    rows: list[dict] = []

    print("bulk GSE137244", flush=True)
    log, kl, kp, setting, scale = load_gse137244()
    rows += sweep_matrix(
        "GSE137244",
        setting,
        scale,
        log,
        kl,
        kp,
        note="Locked cell-line reference. Normal lung excluded.",
    )

    print("bulk GSE137396", flush=True)
    log, kl, kp, setting, scale = load_gse137396()
    rows += sweep_matrix("GSE137396", setting, scale, log, kl, kp, note="Same paper, nodules.")

    print("bulk GSE164758", flush=True)
    log, kl, kp, setting, scale = load_gse164758()
    rows += sweep_matrix(
        "GSE164758",
        setting,
        scale,
        log,
        kl,
        kp,
        note="Primary bulk tumors; stroma included.",
    )

    print("bulk GSE244452", flush=True)
    log, kl, kp, setting, scale = load_gse244452()
    rows += sweep_matrix(
        "GSE244452",
        setting,
        scale,
        log,
        kl,
        kp,
        note="n=3 vs 3; MW floor 0.10.",
    )

    print("bulk GSE274352", flush=True)
    log, kl, kp, setting, scale = load_gse274352()
    rows += sweep_matrix(
        "GSE274352",
        setting,
        scale,
        log,
        kl,
        kp,
        note="Empty-vector lines; n=3 vs 3.",
    )

    print("bulk GSE274351", flush=True)
    log, kl, kp, setting, scale = load_gse274351()
    rows += sweep_matrix("GSE274351", setting, scale, log, kl, kp, note="LCM early adenomas.")

    print("bulk GSE6135", flush=True)
    for name, mm, kl, kp, setting, scale, fair in load_gse6135_cohorts():
        rows += sweep_matrix(
            "GSE6135",
            setting,
            scale,
            mm,
            kl,
            kp,
            fair_genotype=fair,
            note=f"Mouse-level array; cohort={name}.",
        )

    grid = pd.DataFrame(rows)
    grid.to_csv(TABLES / "tacstd2_filter_grid.tsv", sep="\t", index=False)

    fair = grid[grid["fair_n"] & grid["fair_genotype"]].copy()
    # Primary keepers: full-arm filters (not leaveone)
    fair_primary = fair[~fair["filter"].astype(str).str.startswith("leaveone:")].copy()
    # Prefer n such that MW can be <0.05 (need comb floor <0.05 → typically n≥4 vs 4 or better)
    fair_sigcapable = fair_primary[fair_primary["mw_floor"] < 0.05].copy()

    winners = []
    for label, subset, col in (
        ("max_delta_fair_primary", fair_primary, "delta"),
        ("max_abs_t_fair_primary", fair_primary, "abs_t"),
        ("max_delta_sigcapable", fair_sigcapable, "delta"),
        ("max_abs_t_sigcapable", fair_sigcapable, "abs_t"),
        ("max_delta_all_fair_incl_leaveone", fair, "delta"),
        ("max_abs_t_all_fair_incl_leaveone", fair, "abs_t"),
        ("max_delta_any_row", grid, "delta"),
        ("max_abs_t_any_row", grid, "abs_t"),
    ):
        w = select_max(subset, col)
        if w.empty:
            continue
        rec = w.to_dict()
        rec["selection"] = label
        winners.append(rec)
    win_df = pd.DataFrame(winners)
    win_df.to_csv(TABLES / "tacstd2_winners.tsv", sep="\t", index=False)

    # Locked unadjusted GSE137244 none filter
    locked = fair_primary[
        (fair_primary.accession == "GSE137244") & (fair_primary["filter"] == "none")
    ]
    locked_row = locked.iloc[0].to_dict() if len(locked) else {}

    plot_pareto(fair_primary, FIG / "bulk_tacstd2_delta_vs_abs_t.png")

    # Focus table for the writeup: one row per accession at filter=none, plus winners
    none_rows = fair_primary[fair_primary["filter"] == "none"].copy()
    none_rows = none_rows.sort_values("delta", ascending=False)
    none_rows.to_csv(TABLES / "tacstd2_unadjusted_by_accession.tsv", sep="\t", index=False)

    summary = {
        "assay": "bulk",
        "endpoint": "Tacstd2 KL>KP",
        "private_8kl": 0,
        "n_grid_rows": int(len(grid)),
        "n_fair_primary": int(len(fair_primary)),
        "locked_gse137244_none": locked_row,
        "winners": winners,
        "python": platform.python_version(),
        "scipy": scipy_version,
        "rule": (
            "Maximize Tacstd2 Δ and Welch |t| on a fixed filter grid. "
            "Primary reporting keeps full arm n and fair genotype. "
            "Leave-one and n=3 floors are inventoried, not primary."
        ),
    }
    (TABLES / "bulk_summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")

    print(none_rows[["accession", "setting", "n_kl", "n_kp", "delta", "abs_t", "mw_p"]].to_string(index=False))
    print("--- winners ---")
    print(win_df[["selection", "accession", "filter", "n_kl", "n_kp", "delta", "abs_t", "mw_p"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
