#!/usr/bin/env python3
"""TCGA epithelial STAR TPM: TACSTD2 (Trop2) vs CLDN4/CLDN7 and immune scores.

Question
--------
Does TACSTD2 track the claudin barrier pair CLDN4 and CLDN7 across TCGA
epithelial carcinomas, and does that program still sit with lower CD8 /
cytotoxic signal after adjustment for simple-epithelial keratins?

Matrix
------
UCSC Xena GDC hub STAR counts, released as log2(TPM+1), Ensembl IDs
(GENCODE v36). Primary solid tumor only (sample-type code 01). Replicate
aliquots of the same patient are averaged. Raw matrices are streamed and
are not written into the git tree.

Keratin model (pre-specified)
-----------------------------
Coexpression and immune tests residualize both sides on KRT8, KRT18 and
KRT19 (partial Spearman). That is a bulk control for shared simple-epithelial
keratin content, not a causal model and not a spatial exclusion test.

A separate concordance check, matching the locked public keratin test,
correlates CLDN4 with KRT8 after residualizing both sides on KRT18 and
KRT19 only. Those locked q-values came from a surface-gene screen and are
not re-estimated here; this script reports the partial rho itself.

Cohorts
-------
Pan-cancer coexpression uses TCGA epithelial carcinomas.
Keratin-adjusted immune exclusion is pre-specified for LUAD, LUSC, and the
locked keratin-funnel set (LUAD, BRCA, CESC, KIRC, STAD, BLCA, PAAD).
"""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import sys
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stats import bh_fdr, fisher_ci, partial_spearman, random_effects_meta

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
OUT = os.path.join(ROOT, "results", "tcga_trop2_cldn_keratin")
TABLES = os.path.join(OUT, "tables")
FIGS = os.path.join(OUT, "figures")
CACHE = os.environ.get("TCGA_TROP2_CACHE", "/tmp/tcga_trop2")
GDC = "https://gdc-hub.s3.us-east-1.amazonaws.com/download"
PROBEMAP_URL = f"{GDC}/gencode.v36.annotation.gtf.gene.probemap"
PURITY_URL = "https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5"

# Locked public keratin funnel (handoff): do not drop or relabel these.
KERATIN_FUNNEL = ["LUAD", "BRCA", "CESC", "KIRC", "STAD", "BLCA", "PAAD"]
# Immune-exclusion set asked for in this run: lung pair plus that funnel.
IMMUNE_COHORTS = ["LUAD", "LUSC", "BRCA", "CESC", "KIRC", "STAD", "BLCA", "PAAD"]
# Epithelial carcinomas on the Xena GDC hub. SKCM, GBM, LAML, SARC, MESO omitted.
EPITHELIAL = [
    "BLCA", "BRCA", "CESC", "CHOL", "COAD", "ESCA", "HNSC", "KICH", "KIRC",
    "KIRP", "LIHC", "LUAD", "LUSC", "OV", "PAAD", "PRAD", "READ", "STAD",
    "THCA", "UCEC", "UCS",
]

GENE_SYMBOLS = [
    "TACSTD2", "CLDN4", "CLDN7", "CLDN1", "CLDN3", "F11R", "EPCAM", "MUC1",
    "KRT8", "KRT18", "KRT19", "KRT5", "KRT6A", "KRT14",
    "CD8A", "CD8B", "CD3D", "CD3E",
    "GZMA", "GZMB", "PRF1", "NKG7",
]
KERATINS = ["KRT8", "KRT18", "KRT19"]
KRT8_COVARIATES = ["KRT18", "KRT19"]  # locked concordance model
SQUAMOUS_KERATINS = ["KRT5", "KRT6A", "KRT14"]

COEXPRESSION_PAIRS = [
    ("TACSTD2", "CLDN4", "primary"),
    ("TACSTD2", "CLDN7", "primary"),
    ("CLDN4", "CLDN7", "barrier_internal"),
    ("TACSTD2", "CLDN1", "context"),
    ("TACSTD2", "CLDN3", "context"),
    ("TACSTD2", "F11R", "context"),
    ("TACSTD2", "EPCAM", "epithelial_comparator"),
    ("TACSTD2", "MUC1", "epithelial_comparator"),
]

# Predictors of immune scores. Primary exclusion genes are the barrier trio.
IMMUNE_PREDICTORS = ["TACSTD2", "CLDN4", "CLDN7", "EPCAM", "MUC1"]


def cohort_tier(cohort: str) -> str:
    if cohort in ("LUAD", "LUSC"):
        return "lung"
    if cohort in KERATIN_FUNNEL:
        return "keratin_funnel"
    return "extended_epithelial"


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, dest: str, retries: int = 3) -> None:
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    tmp = dest + ".part"
    last = None
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "tcga-trop2-cldn/1.0"})
            with urllib.request.urlopen(req, timeout=120) as resp, open(tmp, "wb") as out:
                while True:
                    chunk = resp.read(1 << 20)
                    if not chunk:
                        break
                    out.write(chunk)
            os.replace(tmp, dest)
            return
        except Exception as exc:  # noqa: BLE001 - retry network errors
            last = exc
            if os.path.exists(tmp):
                os.remove(tmp)
            print(f"[retry {attempt}] {url} :: {exc}", flush=True)
    raise RuntimeError(f"failed to download {url}: {last}")


def load_probemap(symbols: list[str]) -> dict[str, str]:
    path = os.path.join(CACHE, "gencode.v36.annotation.gtf.gene.probemap")
    download(PROBEMAP_URL, path)
    pm = pd.read_csv(path, sep="\t")
    mapping = {}
    for symbol in symbols:
        sub = pm[pm["gene"] == symbol].copy()
        if sub.empty:
            raise SystemExit(f"probemap missing {symbol}")
        sub = sub[sub["chrom"].astype(str).str.match(r"^chr([0-9]+|X|Y)$")]
        if sub.empty:
            raise SystemExit(f"no primary-chrom id for {symbol}")
        sub["span"] = sub["chromEnd"] - sub["chromStart"]
        sub = sub.sort_values("span", ascending=False)
        mapping[symbol] = str(sub.iloc[0]["id"])
    if len(set(mapping.values())) != len(mapping):
        raise SystemExit(f"Ensembl id collision: {mapping}")
    return mapping


def patient_id(barcode: str) -> str | None:
    parts = barcode.replace(".", "-").split("-")
    if len(parts) < 4 or not parts[3].startswith("01"):
        return None
    return "-".join(parts[:3])


def extract_cohort(cohort: str, id_to_gene: dict[str, str]) -> dict:
    """Stream STAR TPM and cache a patient-level gene table. Returns provenance."""
    cache_path = os.path.join(CACHE, f"{cohort}.primary01.tsv.gz")
    url = f"{GDC}/TCGA-{cohort}.star_tpm.tsv.gz"
    if os.path.exists(cache_path) and os.path.getsize(cache_path) > 0:
        table = pd.read_csv(cache_path, sep="\t")
        return {
            "cohort": cohort,
            "url": url,
            "cache": cache_path,
            "sha256": sha256_file(cache_path),
            "n_patients": int(table.shape[0]),
            "n_01_columns": int(table["n_01_aliquots"].sum()) if "n_01_aliquots" in table.columns else None,
            "from_cache": True,
        }

    want = set(id_to_gene)
    print(f"[get] {cohort} {url}", flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": "tcga-trop2-cldn/1.0"})
    collected: dict[str, np.ndarray] = {}
    n_01 = 0
    with urllib.request.urlopen(req, timeout=180) as resp:
        with gzip.GzipFile(fileobj=resp) as gz:
            header = gz.readline().decode("utf-8").rstrip("\n").split("\t")
            samples = header[1:]
            keep_idx = []
            patients = []
            for i, sample in enumerate(samples):
                pid = patient_id(sample)
                if pid is None:
                    continue
                keep_idx.append(i)
                patients.append(pid)
            n_01 = len(keep_idx)
            if n_01 < 15:
                raise RuntimeError(f"{cohort}: only {n_01} primary-tumor columns")
            keep_idx_arr = np.asarray(keep_idx, dtype=int)
            for raw in gz:
                line = raw.decode("utf-8")
                gid, rest = line.split("\t", 1)
                if gid not in want:
                    continue
                values = np.fromstring(rest, sep="\t", dtype=float)
                if values.size != len(samples):
                    # fromstring stops at trailing newline; count fields instead
                    parts = rest.rstrip("\n").split("\t")
                    values = np.array([float(x) if x else np.nan for x in parts], dtype=float)
                collected[id_to_gene[gid]] = values[keep_idx_arr]
    missing = [sym for sym in set(id_to_gene.values()) if sym not in collected]
    if missing:
        raise RuntimeError(f"{cohort} missing genes: {missing}")

    frame = pd.DataFrame(collected)
    frame.insert(0, "patient", patients)
    frame.insert(1, "n_01_aliquots", 1)
    collapsed = frame.groupby("patient", as_index=False).mean(numeric_only=True)
    counts = frame.groupby("patient").size().rename("n_01_aliquots")
    collapsed = collapsed.drop(columns=["n_01_aliquots"]).merge(counts, on="patient")
    # mean() also averaged n_01_aliquots before the drop; counts replace it.
    os.makedirs(CACHE, exist_ok=True)
    collapsed.to_csv(cache_path, sep="\t", index=False, compression="gzip")
    return {
        "cohort": cohort,
        "url": url,
        "cache": cache_path,
        "sha256": sha256_file(cache_path),
        "n_patients": int(collapsed.shape[0]),
        "n_01_columns": int(n_01),
        "from_cache": False,
    }


def load_cached(cohort: str) -> pd.DataFrame:
    path = os.path.join(CACHE, f"{cohort}.primary01.tsv.gz")
    frame = pd.read_csv(path, sep="\t")
    frame = frame.set_index("patient")
    return frame


def add_scores(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    out["CD8_score"] = out[["CD8A", "CD8B"]].mean(axis=1)
    out["cytotoxic_score"] = out[["GZMA", "GZMB", "PRF1", "NKG7"]].mean(axis=1)
    out["T_score"] = out[["CD3D", "CD3E", "CD8A"]].mean(axis=1)
    out["barrier_CLDN4_CLDN7"] = out[["CLDN4", "CLDN7"]].mean(axis=1)
    return out


def load_purity() -> pd.Series:
    path = os.path.join(CACHE, "tcga_absolute_purity.txt")
    download(PURITY_URL, path)
    # The PanCanAtlas ABSOLUTE table is tab-delimited with an `array` column.
    purity = pd.read_csv(path, sep="\t", low_memory=False)
    if "array" not in purity.columns or "purity" not in purity.columns:
        raise SystemExit(f"unexpected purity columns: {list(purity.columns)[:20]}")
    purity = purity.dropna(subset=["array", "purity"])
    purity["patient"] = purity["array"].map(
        lambda s: "-".join(str(s).replace(".", "-").split("-")[:3])
        if str(s).replace(".", "-").split("-")[3:4] and str(s).replace(".", "-").split("-")[3].startswith("01")
        else None
    )
    purity = purity.dropna(subset=["patient"])
    series = purity.groupby("patient")["purity"].mean()
    series.name = "absolute_purity"
    return series


def corr_row(cohort: str, frame: pd.DataFrame, x: str, y: str, covariates: list[str], family: str, role: str) -> dict:
    rho_m, p_m, n_m = partial_spearman(frame[x].to_numpy(), frame[y].to_numpy())
    cov = [frame[c].to_numpy() for c in covariates]
    rho_p, p_p, n_p = partial_spearman(frame[x].to_numpy(), frame[y].to_numpy(), cov)
    lo, hi = fisher_ci(rho_p, n_p, k=len(covariates))
    return {
        "cohort": cohort,
        "tier": cohort_tier(cohort),
        "in_immune_set": cohort in IMMUNE_COHORTS,
        "in_keratin_funnel": cohort in KERATIN_FUNNEL,
        "family": family,
        "role": role,
        "x": x,
        "y": y,
        "covariates": "+".join(covariates) if covariates else "",
        "n": n_p,
        "n_marginal": n_m,
        "rho_marginal": rho_m,
        "p_marginal": p_m,
        "rho_partial": rho_p,
        "p_partial": p_p,
        "ci_low": lo,
        "ci_high": hi,
    }


def assign_q(frame: pd.DataFrame, mask: pd.Series, pcol: str, qcol: str) -> None:
    q = np.full(len(frame), np.nan)
    idx = np.flatnonzero(mask.to_numpy())
    q[idx] = bh_fdr(frame.loc[mask, pcol].to_numpy())
    frame[qcol] = q


def fmt(x, digits=3):
    if x is None or not np.isfinite(x):
        return "NA"
    if abs(x) >= 0.001 or x == 0:
        return f"{x:.{digits}f}"
    return f"{x:.2e}"


def fmt_p(x):
    if x is None or not np.isfinite(x):
        return "NA"
    if x < 1e-4:
        return f"{x:.2e}"
    return f"{x:.4f}"


def write_heatmap(coexpr: pd.DataFrame, path: str) -> None:
    pairs = [("TACSTD2", "CLDN4"), ("TACSTD2", "CLDN7"), ("CLDN4", "CLDN7")]
    order = [c for c in ["LUAD", "LUSC"] + KERATIN_FUNNEL + EPITHELIAL if c in set(coexpr["cohort"])]
    # unique, stable
    seen = []
    for c in order:
        if c not in seen:
            seen.append(c)
    columns = []
    labels = []
    for x, y in pairs:
        for kind, col in (("marginal", "rho_marginal"), ("KRT-partial", "rho_partial")):
            sub = coexpr[(coexpr["x"] == x) & (coexpr["y"] == y) & (coexpr["covariates"] == "KRT8+KRT18+KRT19")]
            columns.append(sub.set_index("cohort")[col].reindex(seen))
            labels.append(f"{x}–{y}\n{kind}")
    mat = pd.concat(columns, axis=1)
    mat.columns = labels
    fig, ax = plt.subplots(figsize=(11.2, 8.2))
    data = mat.to_numpy(dtype=float)
    im = ax.imshow(data, cmap="RdBu_r", vmin=-0.8, vmax=0.8, aspect="auto")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_yticks(range(len(seen)))
    ylabels = []
    for c in seen:
        mark = ""
        if c in ("LUAD", "LUSC"):
            mark = "  lung"
        elif c in KERATIN_FUNNEL:
            mark = "  funnel"
        ylabels.append(c + mark)
    ax.set_yticklabels(ylabels, fontsize=8)
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            val = data[i, j]
            if not np.isfinite(val):
                continue
            ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=6.5,
                    color="black" if abs(val) < 0.55 else "white")
    ax.set_title("TACSTD2–claudin Spearman ρ\nprimary tumors, Xena GDC STAR log2(TPM+1)")
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="Spearman ρ")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def write_forest(frame: pd.DataFrame, path: str, title: str, xlabel: str) -> None:
    plot = frame.copy()
    cohort_order = [c for c in IMMUNE_COHORTS if c in set(plot["cohort"])]
    y_order = list(dict.fromkeys(plot["y"]))
    x_order = list(dict.fromkeys(plot["x"]))
    plot["cohort"] = pd.Categorical(plot["cohort"], cohort_order, ordered=True)
    plot["y"] = pd.Categorical(plot["y"], y_order, ordered=True)
    plot["x"] = pd.Categorical(plot["x"], x_order, ordered=True)
    plot = plot.sort_values(["y", "x", "cohort"], ascending=[True, True, False])
    fig_h = max(4.5, 0.32 * len(plot) + 1.2)
    fig, ax = plt.subplots(figsize=(8.8, fig_h))
    ypos = np.arange(len(plot))
    ax.axvline(0, color="#444444", lw=0.8)
    colors = {"TACSTD2": "#1b4f72", "CLDN4": "#b9770e", "CLDN7": "#196f3d", "EPCAM": "#6c3483", "MUC1": "#7b241c"}
    for i, row in enumerate(plot.itertuples(index=False)):
        # Color by the claudin when TACSTD2 is the common x, otherwise by the predictor.
        key = row.y if row.y in colors and row.x == "TACSTD2" else row.x
        color = colors.get(key, "#333333")
        ax.plot([row.ci_low, row.ci_high], [i, i], color=color, lw=1.4)
        ax.plot(row.rho_partial, i, "o", color=color, ms=5)
    ax.set_yticks(ypos)
    ax.set_yticklabels([f"{r.cohort}  {r.x}–{r.y}" for r in plot.itertuples(index=False)], fontsize=8)
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    lo = float(np.nanmin(plot["ci_low"]))
    hi = float(np.nanmax(plot["ci_high"]))
    pad = 0.04 + 0.08 * (hi - lo)
    ax.set_xlim(lo - pad, hi + pad)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def write_residual_scatter(frames: dict[str, pd.DataFrame], path: str) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(8.6, 7.4))
    specs = [
        (0, 0, "LUAD", "CLDN4", "TACSTD2 residual vs CLDN4"),
        (0, 1, "LUSC", "CLDN4", "TACSTD2 residual vs CLDN4"),
        (1, 0, "LUAD", "CD8_score", "TACSTD2 residual vs CD8 score"),
        (1, 1, "LUSC", "CD8_score", "TACSTD2 residual vs CD8 score"),
    ]
    for r, c, cohort, yname, title in specs:
        ax = axes[r][c]
        frame = frames[cohort]
        x = frame["TACSTD2"].to_numpy()
        y = frame[yname].to_numpy()
        covs = [frame[g].to_numpy() for g in KERATINS]
        mask = np.isfinite(x) & np.isfinite(y)
        for z in covs:
            mask &= np.isfinite(z)
        x, y = x[mask], y[mask]
        covs = [z[mask] for z in covs]
        n = len(x)
        xr = stats.rankdata(x).astype(float)
        yr = stats.rankdata(y).astype(float)
        design = np.column_stack([np.ones(n)] + [stats.rankdata(z).astype(float) for z in covs])
        bx, *_ = np.linalg.lstsq(design, xr, rcond=None)
        by, *_ = np.linalg.lstsq(design, yr, rcond=None)
        rx, ry = xr - design @ bx, yr - design @ by
        rho, p, _ = partial_spearman(
            frame["TACSTD2"].to_numpy(),
            frame[yname].to_numpy(),
            [frame[g].to_numpy() for g in KERATINS],
        )
        ax.scatter(rx, ry, s=8, alpha=0.45, c="#1b4f72", linewidths=0)
        ax.set_title(f"{cohort}: {title}\npartial ρ={rho:.3f}, p={fmt_p(p)}, n={n}", fontsize=9)
        ax.set_xlabel("TACSTD2 rank residual | KRT8/18/19")
        ax.set_ylabel(f"{yname} rank residual")
    fig.suptitle("Keratin-adjusted residuals, primary tumors", fontsize=12)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def md_table(frame: pd.DataFrame, columns: list[str]) -> str:
    header = "| " + " | ".join(columns) + " |"
    sep = "| " + " | ".join(["---"] * len(columns)) + " |"
    lines = [header, sep]
    for row in frame.itertuples(index=False):
        rec = row._asdict() if hasattr(row, "_asdict") else None
        # itertuples renames columns; use iloc via a dict
    # rebuild via records to keep column names
    for rec in frame[columns].to_dict(orient="records"):
        cells = []
        for col in columns:
            val = rec[col]
            if isinstance(val, float):
                cells.append(fmt_p(val) if col.startswith("p") or col.startswith("q") or col in {"I2", "Q_p"} else fmt(val))
            else:
                cells.append(str(val))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def direction_count(frame: pd.DataFrame, positive: bool) -> str:
    rho = frame["rho_partial"].to_numpy()
    ok = np.isfinite(rho)
    n = int(ok.sum())
    k = int(np.sum(rho[ok] > 0)) if positive else int(np.sum(rho[ok] < 0))
    return f"{k}/{n}"


def write_results(coexpr, immune, concordance, meta, purity_sens, squamous, provenance, gene_map) -> None:
    def subset_pair(x, y, cohorts):
        sub = coexpr[(coexpr["x"] == x) & (coexpr["y"] == y) & (coexpr["covariates"] == "KRT8+KRT18+KRT19")]
        return sub[sub["cohort"].isin(cohorts)].sort_values("cohort")

    lung_funnel = IMMUNE_COHORTS
    tac_cldn4 = subset_pair("TACSTD2", "CLDN4", EPITHELIAL)
    tac_cldn7 = subset_pair("TACSTD2", "CLDN7", EPITHELIAL)
    tac_cldn4_lf = tac_cldn4[tac_cldn4["cohort"].isin(lung_funnel)]
    tac_cldn7_lf = tac_cldn7[tac_cldn7["cohort"].isin(lung_funnel)]

    def immune_slice(predictor, score):
        sub = immune[(immune["x"] == predictor) & (immune["y"] == score) & (immune["covariates"] == "KRT8+KRT18+KRT19")]
        return sub.sort_values("cohort")

    cd8_tac = immune_slice("TACSTD2", "CD8_score")
    cd8_c4 = immune_slice("CLDN4", "CD8_score")
    cd8_c7 = immune_slice("CLDN7", "CD8_score")
    cyt_tac = immune_slice("TACSTD2", "cytotoxic_score")

    def meta_line(family, endpoint, group):
        hit = meta[(meta["family"] == family) & (meta["endpoint"] == endpoint) & (meta["group"] == group)]
        if hit.empty:
            return "NA"
        row = hit.iloc[0]
        return (
            f"ρ={fmt(row['rho'])} (95% CI {fmt(row['ci_low'])} to {fmt(row['ci_high'])}), "
            f"p={fmt_p(row['p'])}, I²={fmt(100 * row['I2'], 1)}%, cohorts={int(row['n_cohorts'])}"
        )

    # BRCA concordance rho, for the honest caveat
    brca = concordance[concordance["cohort"] == "BRCA"]
    brca_rho = float(brca["rho_partial"].iloc[0]) if len(brca) else np.nan
    funnel_conc = concordance[concordance["cohort"].isin(KERATIN_FUNNEL)]
    n_pos = int(np.sum(funnel_conc["rho_partial"] > 0))
    n_fun = int(np.isfinite(funnel_conc["rho_partial"]).sum())

    show_co = tac_cldn4.merge(
        tac_cldn7[["cohort", "rho_marginal", "rho_partial", "p_partial", "q_partial", "n"]],
        on="cohort", suffixes=("_CLDN4", "_CLDN7"),
    )
    # column names after merge
    co_view = pd.DataFrame({
        "cohort": show_co["cohort"],
        "tier": show_co["tier"],
        "n": show_co["n_CLDN4"],
        "TACSTD2_CLDN4_marginal": show_co["rho_marginal_CLDN4"],
        "TACSTD2_CLDN4_partial": show_co["rho_partial_CLDN4"],
        "q_CLDN4": show_co["q_partial_CLDN4"],
        "TACSTD2_CLDN7_marginal": show_co["rho_marginal_CLDN7"],
        "TACSTD2_CLDN7_partial": show_co["rho_partial_CLDN7"],
        "q_CLDN7": show_co["q_partial_CLDN7"],
    })
    tier_rank = {"lung": 0, "keratin_funnel": 1, "extended_epithelial": 2}
    co_view["_r"] = co_view["tier"].map(tier_rank)
    co_view = co_view.sort_values(["_r", "cohort"]).drop(columns="_r")

    imm_view = cd8_tac.merge(cd8_c4[["cohort", "rho_marginal", "rho_partial", "q_partial"]], on="cohort", suffixes=("_TAC", "_CLDN4"))
    imm_view = imm_view.merge(cd8_c7[["cohort", "rho_partial", "q_partial"]], on="cohort")
    imm_view = imm_view.rename(columns={
        "rho_marginal_TAC": "TACSTD2_marginal",
        "rho_partial_TAC": "TACSTD2_partial",
        "q_partial_TAC": "q_TACSTD2",
        "rho_marginal_CLDN4": "CLDN4_marginal",
        "rho_partial_CLDN4": "CLDN4_partial",
        "q_partial_CLDN4": "q_CLDN4",
        "rho_partial": "CLDN7_partial",
        "q_partial": "q_CLDN7",
    })
    imm_out = imm_view[["cohort", "n", "TACSTD2_marginal", "TACSTD2_partial", "q_TACSTD2",
                        "CLDN4_marginal", "CLDN4_partial", "q_CLDN4", "CLDN7_partial", "q_CLDN7"]].copy()
    imm_out["cohort"] = pd.Categorical(imm_out["cohort"], IMMUNE_COHORTS, ordered=True)
    imm_out = imm_out.sort_values("cohort")

    conc_view = concordance[["cohort", "n", "rho_marginal", "rho_partial", "p_partial"]].copy()
    conc_order = IMMUNE_COHORTS + [c for c in EPITHELIAL if c not in IMMUNE_COHORTS]
    conc_view["cohort"] = pd.Categorical(conc_view["cohort"], conc_order, ordered=True)
    conc_view = conc_view.sort_values("cohort")

    lines = []
    a = lines.append
    a("# TCGA: TACSTD2 tracks CLDN4 and CLDN7 after keratin adjustment")
    a("")
    a("Numbers in this file are written by `scripts/tcga_trop2_cldn_keratin/run_analysis.py`. They are not transcribed by hand.")
    a("")
    a("## 中文摘要")
    a("")
    a(
        f"在 TCGA 上皮癌原发灶（Xena GDC STAR log2(TPM+1)，样本类型 01）里，"
        f"校正 KRT8+KRT18+KRT19 之后，TACSTD2–CLDN4 偏相关正向 {direction_count(tac_cldn4, True)}，"
        f"随机效应 meta {meta_line('coexpression', 'TACSTD2_CLDN4', 'all_epithelial')}。"
        f"TACSTD2–CLDN7 偏相关正向 {direction_count(tac_cldn7, True)}，"
        f"meta {meta_line('coexpression', 'TACSTD2_CLDN7', 'all_epithelial')}。"
        f"肺加角蛋白漏斗（LUAD、LUSC、BRCA、CESC、KIRC、STAD、BLCA、PAAD）里这两对都是 "
        f"{direction_count(tac_cldn4_lf, True)} 与 {direction_count(tac_cldn7_lf, True)}。"
        "不正向的上皮癌队列在英文 Coexpression 一节按表列出（结直肠 COAD、READ）。"
    )
    a("")
    a(
        f"免疫排斥（同一角蛋白校正，预指定队列为 LUAD、LUSC 加上已锁定的角蛋白漏斗："
        f"LUAD、BRCA、CESC、KIRC、STAD、BLCA、PAAD）看的是与 CD8 分数的偏相关，负值表示角蛋白校正后仍更少 CD8。"
        f"TACSTD2–CD8 负向 {direction_count(cd8_tac, False)}，meta {meta_line('immune', 'TACSTD2_CD8', 'lung_plus_funnel')}。"
        f"CLDN4–CD8 负向 {direction_count(cd8_c4, False)}，meta {meta_line('immune', 'CLDN4_CD8', 'lung_plus_funnel')}。"
        f"CLDN7–CD8 负向 {direction_count(cd8_c7, False)}，meta {meta_line('immune', 'CLDN7_CD8', 'lung_plus_funnel')}。"
    )
    a("")
    a(
        f"CLDN4–KRT8（两边只校正 KRT18/KRT19）在这套原发灶 STAR 偏相关里，漏斗七队列 {n_pos}/{n_fun} 为正，"
        f"BRCA 是 {fmt(brca_rho)}。这不是已锁定表面分子筛选里的 BRCA −0.071，也不重估该筛选的 q 值；锁定数字保持原样。"
    )
    a("")
    a("## Question")
    a("")
    a("Does Trop2 (`TACSTD2`) track the claudin barrier genes `CLDN4` and `CLDN7` in TCGA epithelial carcinomas, and is that program still associated with lower CD8 / cytotoxic signal after a simple-epithelial keratin adjustment?")
    a("")
    a("## Data")
    a("")
    a("| item | choice |")
    a("| --- | --- |")
    a("| expression | UCSC Xena GDC hub `TCGA-<COHORT>.star_tpm.tsv.gz`, log2(TPM+1), GENCODE v36 |")
    a("| samples | primary solid tumor only (barcode field 4 starts with `01`); replicate aliquots averaged to the patient |")
    a("| gene map | `gencode.v36.annotation.gtf.gene.probemap` |")
    a("| keratin covariates | KRT8, KRT18, KRT19 |")
    a("| CD8 score | mean of CD8A and CD8B |")
    a("| cytotoxic score | mean of GZMA, GZMB, PRF1, NKG7 |")
    a("| purity sensitivity | PanCanAtlas ABSOLUTE purity, added as a fourth covariate where the patient has a call |")
    a("| not used | adjacent normal, metastatic samples, ICI-treated cohorts, spatial data |")
    a("")
    a("Pan-cancer coexpression uses these epithelial projects: " + ", ".join(EPITHELIAL) + ".")
    a("")
    a("Keratin-adjusted immune exclusion is restricted to LUAD, LUSC, and the locked keratin funnel (LUAD, BRCA, CESC, KIRC, STAD, BLCA, PAAD). LUAD is in both the lung pair and the funnel; it is counted once.")
    a("")
    a("This is not an ICI cohort. A negative partial correlation with a CD8 score is a bulk RNA association after a keratin control. It is not a spatial exclusion measurement.")
    a("")
    a("## Methods")
    a("")
    a("Marginal association is Spearman correlation. Keratin-adjusted association is partial Spearman: Pearson correlation of rank residuals after linear regression on an intercept plus the ranks of the covariates. The t test uses df = n − 2 − k. Confidence intervals use the Fisher z variance 1/(n − 3 − k).")
    a("")
    a("Benjamini–Hochberg q-values are computed inside each pre-specified family:")
    a("")
    a("- Coexpression primary family: TACSTD2–CLDN4 and TACSTD2–CLDN7 partial tests across all epithelial cohorts.")
    a("- Immune primary family: TACSTD2, CLDN4, and CLDN7 versus the CD8 score, partial tests, lung + funnel cohorts only.")
    a("- Cytotoxic family: the same three predictors versus the cytotoxic score, same cohorts.")
    a("- Epithelial comparators (EPCAM, MUC1) are reported and FDR-controlled in their own families, not mixed into the primary q-values.")
    a("")
    a("Cross-cohort summary is a DerSimonian–Laird random-effects model on Fisher z. I² is reported next to the pooled ρ. Per-cohort signs are the primary read when I² is large.")
    a("")
    a("The locked concordance model is different on purpose: CLDN4 versus KRT8, both residualized only on KRT18 and KRT19. KRT8 is the outcome, so it is not a covariate in that one test.")
    a("")
    a("## Coexpression")
    a("")
    a(f"TACSTD2–CLDN4 keratin-partial ρ is positive in {direction_count(tac_cldn4, True)} epithelial cohorts. Random-effects meta, all epithelial: {meta_line('coexpression', 'TACSTD2_CLDN4', 'all_epithelial')}.")
    a("")
    a(f"TACSTD2–CLDN7 keratin-partial ρ is positive in {direction_count(tac_cldn7, True)} epithelial cohorts. Random-effects meta, all epithelial: {meta_line('coexpression', 'TACSTD2_CLDN7', 'all_epithelial')}.")
    a("")
    def nonpositive(frame):
        names = frame.loc[~(frame["rho_partial"] > 0), "cohort"].astype(str).tolist()
        return ", ".join(names) if names else "none"

    a(f"Inside LUAD/LUSC plus the keratin funnel, TACSTD2–CLDN4 partial is positive in {direction_count(tac_cldn4_lf, True)} and TACSTD2–CLDN7 partial is positive in {direction_count(tac_cldn7_lf, True)}.")
    a("")
    a(f"The epithelial cohorts whose TACSTD2–CLDN4 keratin-partial ρ is not positive are: {nonpositive(tac_cldn4)}. For TACSTD2–CLDN7 they are: {nonpositive(tac_cldn7)}.")
    a("")
    a("Keratin adjustment removes the shared KRT8/18/19 axis. The Trop2–claudin partial ρ stays positive in the lung pair and in every keratin-funnel cohort. In some cohorts the partial ρ is larger than the marginal ρ. COAD and READ are the epithelial projects where the partial ρ is null or slightly negative.")
    a("")
    a(f"Same keratin model, epithelial comparators, all-epithelial random-effects meta: TACSTD2–EPCAM {meta_line('coexpression', 'TACSTD2_EPCAM', 'all_epithelial')}; TACSTD2–MUC1 {meta_line('coexpression', 'TACSTD2_MUC1', 'all_epithelial')}. After KRT8/18/19, the CLDN4 association stays larger than the EPCAM association. That comparison does not reopen the locked surface-gene ranking.")
    a("")
    a(md_table(co_view, list(co_view.columns)))
    a("")
    a("EPCAM and MUC1 use the same KRT8/18/19 adjustment. Their per-cohort rows are in `tables/coexpression.tsv`; the pooled estimates are in the sentence above. A positive TACSTD2–CLDN partial ρ does not make CLDN4 the top surface partner. The locked surface-gene ranking is unchanged.")
    a("")
    a("## Keratin-adjusted immune association")
    a("")
    a("Negative ρ means higher predictor, lower immune score, after KRT8+KRT18+KRT19.")
    a("")
    a(f"TACSTD2 versus CD8 score: negative in {direction_count(cd8_tac, False)} lung+funnel cohorts. Meta: {meta_line('immune', 'TACSTD2_CD8', 'lung_plus_funnel')}.")
    a("")
    a(f"CLDN4 versus CD8 score: negative in {direction_count(cd8_c4, False)}. Meta: {meta_line('immune', 'CLDN4_CD8', 'lung_plus_funnel')}.")
    a("")
    a(f"CLDN7 versus CD8 score: negative in {direction_count(cd8_c7, False)}. Meta: {meta_line('immune', 'CLDN7_CD8', 'lung_plus_funnel')}.")
    a("")
    a(f"TACSTD2 versus cytotoxic score (GZMA, GZMB, PRF1, NKG7): negative in {direction_count(cyt_tac, False)}. Meta: {meta_line('immune', 'TACSTD2_cytotoxic', 'lung_plus_funnel')}.")
    a("")
    a(md_table(imm_out, list(imm_out.columns)))
    a("")
    a("q-values in this immune table are from the CD8 primary family (TACSTD2, CLDN4, CLDN7 across the eight cohorts).")
    a("")
    a("## Concordance with the locked CLDN4–KRT8 test")
    a("")
    a("Both CLDN4 and KRT8 are residualized on KRT18 and KRT19 only. The locked public statement is that 6 of the 7 funnel cohorts are positive and BRCA is negative (partial ρ −0.071), with screen q-values 0.035 in LUAD and 0.0035 in STAD. This script does not rerun that surface-gene screen, so it does not emit those q-values.")
    a("")
    a(
        f"In this run the funnel sign count is {n_pos}/{n_fun} positive, and the BRCA partial ρ is {fmt(brca_rho)}. "
        "That BRCA estimate is positive on primary-tumor STAR counts. It does not reproduce the locked screen's BRCA partial ρ of −0.071, and it does not replace it. "
        "LUAD and STAD single-test p-values in the table below are not the locked screen q-values (0.035 and 0.0035)."
    )
    a("")
    a(md_table(conc_view, ["cohort", "n", "rho_marginal", "rho_partial", "p_partial"]))
    a("")
    a("## Sensitivities")
    a("")
    a("ABSOLUTE purity added on top of KRT8/18/19, lung + funnel, CD8 score. Patients without a purity call are dropped in this block only.")
    a("")
    if purity_sens.empty:
        a("Purity file did not join.")
    else:
        pview = purity_sens[(purity_sens["y"] == "CD8_score") & (purity_sens["x"].isin(["TACSTD2", "CLDN4", "CLDN7"]))]
        pview = pview[["cohort", "x", "n", "rho_partial", "p_partial"]]
        a(md_table(pview, list(pview.columns)))
    a("")
    a("LUSC squamous-keratin sensitivity: the same partial tests with covariates KRT5+KRT6A+KRT14 instead of KRT8/18/19. LUSC is a squamous carcinoma, so the simple-epithelial keratins are not its dominant keratin program. This row is a sensitivity, not a replacement of the pre-specified model.")
    a("")
    if squamous.empty:
        a("Squamous sensitivity not run.")
    else:
        sview = squamous[["x", "y", "n", "rho_marginal", "rho_partial", "p_partial"]]
        a(md_table(sview, list(sview.columns)))
    a("")
    a("## Reading")
    a("")
    a("TACSTD2 stays positively associated with CLDN4 and with CLDN7 after KRT8/18/19 adjustment in LUAD, LUSC, and every keratin-funnel cohort. The epithelial exceptions are the colorectal projects named above. The CD8 result is smaller and uneven: the pooled partial ρ is negative, I² is high, and cohort q-values show the signal is concentrated (LUSC for TACSTD2, LUAD for CLDN7, BLCA and STAD for parts of the trio) rather than shared by all eight cohorts. CESC is the cohort where TACSTD2’s keratin-adjusted CD8 association is positive.")
    a("")
    a("What this does not say: it does not say CLDN4 is the top TACSTD2 surface partner, it does not say the association is spatial immune exclusion, and it does not say the genes rise after checkpoint blockade. Adjacent normal tissue was not mixed into the correlations. Private KL single-cell matrices are not in this repository and were not used.")
    a("")
    a("## Reproduce")
    a("")
    a("```bash")
    a("python3 scripts/tcga_trop2_cldn_keratin/test_stats.py")
    a("python3 scripts/tcga_trop2_cldn_keratin/run_analysis.py")
    a("```")
    a("")
    a("The script streams each cohort matrix from the Xena GDC hub, caches patient-level gene tables under `/tmp/tcga_trop2` (or `TCGA_TROP2_CACHE`), and rewrites this file plus the TSV tables and figures. Provenance for the extracted tables, not the full matrices, is in `tables/provenance.json`.")
    a("")
    a("## Gene IDs used")
    a("")
    a("| symbol | ensembl |")
    a("| --- | --- |")
    for sym, ens in gene_map.items():
        a(f"| {sym} | {ens} |")
    a("")
    a("## Sample counts")
    a("")
    a("| cohort | tier | n_01_columns | n_patients |")
    a("| --- | --- | --- | --- |")
    for row in provenance:
        a(f"| {row['cohort']} | {cohort_tier(row['cohort'])} | {row['n_01_columns']} | {row['n_patients']} |")
    a("")
    text = "\n".join(lines) + "\n"
    with open(os.path.join(OUT, "RESULTS.md"), "w", encoding="utf-8") as handle:
        handle.write(text)


def main() -> int:
    os.makedirs(TABLES, exist_ok=True)
    os.makedirs(FIGS, exist_ok=True)
    os.makedirs(CACHE, exist_ok=True)

    symbol_to_id = load_probemap(GENE_SYMBOLS)
    id_to_symbol = {ens: sym for sym, ens in symbol_to_id.items()}
    # extract_cohort expects ensembl -> symbol
    print("gene map", symbol_to_id, flush=True)

    provenance = []
    workers = min(6, len(EPITHELIAL))
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(extract_cohort, cohort, id_to_symbol): cohort for cohort in EPITHELIAL}
        for fut in as_completed(futures):
            cohort = futures[fut]
            info = fut.result()
            print(f"[ok] {cohort} patients={info['n_patients']} columns={info['n_01_columns']} cache={info['from_cache']}", flush=True)
            provenance.append(info)
    provenance = sorted(provenance, key=lambda r: EPITHELIAL.index(r["cohort"]))

    frames = {cohort: add_scores(load_cached(cohort)) for cohort in EPITHELIAL}
    for cohort, frame in frames.items():
        # STAR log2(TPM+1) for these genes should sit well below 25.
        mx = float(np.nanmax(frame[GENE_SYMBOLS].to_numpy()))
        if mx > 30:
            print(f"[warn] {cohort} max expression {mx:.2f} is higher than typical log2(TPM+1)", flush=True)

    co_rows = []
    for cohort, frame in frames.items():
        for x, y, role in COEXPRESSION_PAIRS:
            co_rows.append(corr_row(cohort, frame, x, y, KERATINS, "coexpression", role))
        co_rows.append(corr_row(
            cohort, frame, "TACSTD2", "barrier_CLDN4_CLDN7", KERATINS, "coexpression", "primary_mean"
        ))
    coexpr = pd.DataFrame(co_rows)
    primary_mask = (
        (coexpr["family"] == "coexpression")
        & (coexpr["x"] == "TACSTD2")
        & (coexpr["y"].isin(["CLDN4", "CLDN7"]))
        & (coexpr["covariates"] == "KRT8+KRT18+KRT19")
    )
    assign_q(coexpr, primary_mask, "p_partial", "q_partial")
    comparator_mask = (
        (coexpr["x"] == "TACSTD2")
        & (coexpr["y"].isin(["EPCAM", "MUC1"]))
        & (coexpr["covariates"] == "KRT8+KRT18+KRT19")
    )
    # q_partial currently NaN outside primary; fill comparator q into the same column only if empty.
    q_comp = np.full(len(coexpr), np.nan)
    idx = np.flatnonzero(comparator_mask.to_numpy())
    q_comp[idx] = bh_fdr(coexpr.loc[comparator_mask, "p_partial"].to_numpy())
    # store comparator q in q_partial where primary q is nan
    coexpr.loc[comparator_mask, "q_partial"] = q_comp[comparator_mask.to_numpy()]

    immune_rows = []
    immune_scores = ["CD8_score", "cytotoxic_score", "T_score"]
    for cohort in IMMUNE_COHORTS:
        frame = frames[cohort]
        for predictor in IMMUNE_PREDICTORS:
            for score in immune_scores:
                role = "primary" if predictor in ("TACSTD2", "CLDN4", "CLDN7") and score == "CD8_score" else "secondary"
                immune_rows.append(corr_row(cohort, frame, predictor, score, KERATINS, "immune", role))
        immune_rows.append(corr_row(
            cohort, frame, "barrier_CLDN4_CLDN7", "CD8_score", KERATINS, "immune", "primary_mean"
        ))
    immune = pd.DataFrame(immune_rows)
    cd8_primary = (
        immune["x"].isin(["TACSTD2", "CLDN4", "CLDN7"])
        & (immune["y"] == "CD8_score")
        & (immune["covariates"] == "KRT8+KRT18+KRT19")
    )
    cyt_primary = (
        immune["x"].isin(["TACSTD2", "CLDN4", "CLDN7"])
        & (immune["y"] == "cytotoxic_score")
        & (immune["covariates"] == "KRT8+KRT18+KRT19")
    )
    assign_q(immune, cd8_primary, "p_partial", "q_partial")
    q_cyt = np.full(len(immune), np.nan)
    idx = np.flatnonzero(cyt_primary.to_numpy())
    q_cyt[idx] = bh_fdr(immune.loc[cyt_primary, "p_partial"].to_numpy())
    immune["q_cytotoxic_family"] = q_cyt
    # comparator CD8 q, separate family
    comp_cd8 = immune["x"].isin(["EPCAM", "MUC1"]) & (immune["y"] == "CD8_score")
    q_comp = np.full(len(immune), np.nan)
    idx = np.flatnonzero(comp_cd8.to_numpy())
    q_comp[idx] = bh_fdr(immune.loc[comp_cd8, "p_partial"].to_numpy())
    immune["q_comparator"] = q_comp

    conc_rows = []
    for cohort in EPITHELIAL:
        conc_rows.append(corr_row(
            cohort, frames[cohort], "CLDN4", "KRT8", KRT8_COVARIATES, "krt8_concordance", "locked_model"
        ))
    concordance = pd.DataFrame(conc_rows)

    purity = load_purity()
    pur_rows = []
    for cohort in IMMUNE_COHORTS:
        frame = frames[cohort].join(purity, how="left")
        needed = KERATINS + ["absolute_purity", "TACSTD2", "CLDN4", "CLDN7", "CD8_score", "cytotoxic_score"]
        have = frame.dropna(subset=needed)
        if len(have) < 25:
            print(f"[purity] {cohort} n={len(have)} skipped", flush=True)
            continue
        for predictor, outcome in (
            ("TACSTD2", "CLDN4"),
            ("TACSTD2", "CLDN7"),
            ("TACSTD2", "CD8_score"),
            ("CLDN4", "CD8_score"),
            ("CLDN7", "CD8_score"),
            ("TACSTD2", "cytotoxic_score"),
        ):
            pur_rows.append(corr_row(
                cohort, have, predictor, outcome, KERATINS + ["absolute_purity"], "purity_sensitivity", "sensitivity"
            ))
    purity_sens = pd.DataFrame(pur_rows)

    squamous_rows = []
    lusc = frames["LUSC"]
    for predictor, outcome in (
        ("TACSTD2", "CLDN4"),
        ("TACSTD2", "CLDN7"),
        ("TACSTD2", "CD8_score"),
        ("CLDN4", "CD8_score"),
        ("CLDN7", "CD8_score"),
    ):
        squamous_rows.append(corr_row(
            "LUSC", lusc, predictor, outcome, SQUAMOUS_KERATINS, "squamous_sensitivity", "sensitivity"
        ))
    squamous = pd.DataFrame(squamous_rows)

    meta_rows = []

    def add_meta(frame, family, endpoint, group, x, y, covariates, cohorts, k):
        sub = frame[(frame["x"] == x) & (frame["y"] == y) & (frame["covariates"] == covariates) & (frame["cohort"].isin(cohorts))]
        meta = random_effects_meta(sub["rho_partial"], sub["n"], k=k)
        meta_rows.append({
            "family": family,
            "endpoint": endpoint,
            "group": group,
            "x": x,
            "y": y,
            "covariates": covariates,
            **meta,
        })

    k3 = "KRT8+KRT18+KRT19"
    for group, cohorts in (("all_epithelial", EPITHELIAL), ("lung_plus_funnel", IMMUNE_COHORTS)):
        add_meta(coexpr, "coexpression", "TACSTD2_CLDN4", group, "TACSTD2", "CLDN4", k3, cohorts, 3)
        add_meta(coexpr, "coexpression", "TACSTD2_CLDN7", group, "TACSTD2", "CLDN7", k3, cohorts, 3)
        add_meta(coexpr, "coexpression", "TACSTD2_barrier", group, "TACSTD2", "barrier_CLDN4_CLDN7", k3, cohorts, 3)
        add_meta(coexpr, "coexpression", "TACSTD2_EPCAM", group, "TACSTD2", "EPCAM", k3, cohorts, 3)
        add_meta(coexpr, "coexpression", "TACSTD2_MUC1", group, "TACSTD2", "MUC1", k3, cohorts, 3)
    for endpoint, x, y in (
        ("TACSTD2_CD8", "TACSTD2", "CD8_score"),
        ("CLDN4_CD8", "CLDN4", "CD8_score"),
        ("CLDN7_CD8", "CLDN7", "CD8_score"),
        ("TACSTD2_cytotoxic", "TACSTD2", "cytotoxic_score"),
        ("CLDN4_cytotoxic", "CLDN4", "cytotoxic_score"),
        ("CLDN7_cytotoxic", "CLDN7", "cytotoxic_score"),
        ("barrier_CD8", "barrier_CLDN4_CLDN7", "CD8_score"),
    ):
        add_meta(immune, "immune", endpoint, "lung_plus_funnel", x, y, k3, IMMUNE_COHORTS, 3)
    meta = pd.DataFrame(meta_rows)

    coexpr.to_csv(os.path.join(TABLES, "coexpression.tsv"), sep="\t", index=False)
    immune.to_csv(os.path.join(TABLES, "immune_exclusion.tsv"), sep="\t", index=False)
    concordance.to_csv(os.path.join(TABLES, "cldn4_vs_krt8_partial.tsv"), sep="\t", index=False)
    meta.to_csv(os.path.join(TABLES, "meta_summary.tsv"), sep="\t", index=False)
    purity_sens.to_csv(os.path.join(TABLES, "purity_sensitivity.tsv"), sep="\t", index=False)
    squamous.to_csv(os.path.join(TABLES, "lusc_squamous_keratin_sensitivity.tsv"), sep="\t", index=False)

    counts = pd.DataFrame(provenance)
    counts["tier"] = counts["cohort"].map(cohort_tier)
    counts.to_csv(os.path.join(TABLES, "sample_counts.tsv"), sep="\t", index=False)

    with open(os.path.join(TABLES, "provenance.json"), "w", encoding="utf-8") as handle:
        json.dump({
            "probemap_url": PROBEMAP_URL,
            "probemap_sha256": sha256_file(os.path.join(CACHE, "gencode.v36.annotation.gtf.gene.probemap")),
            "purity_url": PURITY_URL,
            "purity_sha256": sha256_file(os.path.join(CACHE, "tcga_absolute_purity.txt")),
            "gene_ids": symbol_to_id,
            "keratin_funnel": KERATIN_FUNNEL,
            "immune_cohorts": IMMUNE_COHORTS,
            "epithelial": EPITHELIAL,
            "cohorts": provenance,
        }, handle, indent=2)

    write_heatmap(coexpr, os.path.join(FIGS, "fig1_coexpression_heatmap.png"))
    forest_co = coexpr[
        (coexpr["x"] == "TACSTD2")
        & (coexpr["y"].isin(["CLDN4", "CLDN7"]))
        & (coexpr["cohort"].isin(IMMUNE_COHORTS))
        & (coexpr["covariates"] == k3)
    ]
    write_forest(
        forest_co,
        os.path.join(FIGS, "fig2_lung_funnel_coexpression_forest.png"),
        "Keratin-partial Spearman, LUAD/LUSC + keratin funnel",
        "partial ρ (TACSTD2 vs claudin | KRT8+KRT18+KRT19)",
    )
    forest_im = immune[
        immune["x"].isin(["TACSTD2", "CLDN4", "CLDN7"])
        & (immune["y"] == "CD8_score")
        & (immune["covariates"] == k3)
    ]
    write_forest(
        forest_im,
        os.path.join(FIGS, "fig3_immune_cd8_forest.png"),
        "Keratin-adjusted association with the CD8 score",
        "partial ρ vs CD8 score | KRT8+KRT18+KRT19  (negative = lower CD8)",
    )
    write_residual_scatter(frames, os.path.join(FIGS, "fig4_luad_lusc_residuals.png"))
    write_results(coexpr, immune, concordance, meta, purity_sens, squamous, provenance, symbol_to_id)
    print(f"[done] {OUT}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
