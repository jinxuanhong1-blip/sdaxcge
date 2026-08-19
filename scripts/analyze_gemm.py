#!/usr/bin/env python3
"""Public KP/KL anti-PD1 RNA-seq: Tacstd2, Cldn4, TJ score (ICB − naive)."""

from __future__ import annotations

import gzip
import io
import json
import tarfile
from pathlib import Path
from urllib.request import Request, urlopen

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import binomtest, wilcoxon

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "gemm"
OUT = ROOT / "results" / "gemm"
FIG = OUT / "figures"
TABLES = OUT / "tables"
UA = "tismo-icb-recompute/1.0"

FOCUS = ["Tacstd2", "Cldn4", "Cldn3", "Cldn6", "Cldn7", "Cdh1", "F11r", "Ocln"]
TJ = ["Cldn3", "Cldn4", "Cldn6", "Cldn7", "Cdh1", "F11r", "Ocln"]

# Ensembl 109/GRCm39 symbols (verified via REST below; fallback hardcoded).
ENSEMBL_FALLBACK = {
    "Tacstd2": "ENSMUSG00000051397",
    "Cldn4": "ENSMUSG00000047501",
    "Cldn3": "ENSMUSG00000070473",
    "Cldn6": "ENSMUSG00000043454",
    "Cldn7": "ENSMUSG00000018569",
    "Cdh1": "ENSMUSG00000000303",
    "F11r": "ENSMUSG00000038235",
    "Ocln": "ENSMUSG00000021638",
}


def fetch(url: str) -> bytes:
    req = Request(url, headers={"User-Agent": UA})
    with urlopen(req, timeout=300) as resp:
        return resp.read()


def fetch_text(url: str) -> str:
    blob = fetch(url)
    if url.endswith(".gz"):
        blob = gzip.decompress(blob)
    return blob.decode("utf-8", errors="replace")


def ensembl_map() -> dict[str, str]:
    out = dict(ENSEMBL_FALLBACK)
    for sym in FOCUS:
        try:
            raw = fetch(f"https://rest.ensembl.org/lookup/symbol/mus_musculus/{sym}?content-type=application/json")
            rec = json.loads(raw.decode())
            if rec.get("id"):
                out[sym] = rec["id"]
        except Exception:
            pass
    return out


def write_focus(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t")


def resolve_genes(df: pd.DataFrame, ens: dict[str, str]) -> pd.DataFrame:
    """Return gene x sample matrix with FOCUS index (symbols)."""
    idx = df.index.astype(str)
    # already symbols
    if set(FOCUS) & set(idx):
        keep = df.loc[[g for g in FOCUS if g in idx]]
        return keep
    # Ensembl
    rev = {v: k for k, v in ens.items()}
    mapped = []
    for i in idx:
        key = i.split(".")[0]
        if key in rev:
            mapped.append((rev[key], i))
    if mapped:
        sub = df.loc[[i for _, i in mapped]].copy()
        sub.index = [s for s, _ in mapped]
        return sub.groupby(level=0).mean()
    # gene symbol column
    for col in df.columns:
        if df[col].astype(str).isin(FOCUS).any():
            tmp = df[df[col].astype(str).isin(FOCUS)].copy()
            tmp = tmp.set_index(col)
            num = tmp.select_dtypes(include=[np.number])
            return num.groupby(level=0).mean().reindex([g for g in FOCUS if g in num.index.astype(str) or True])
    return pd.DataFrame()


def pair_stats(naive: np.ndarray, icb: np.ndarray) -> dict:
    naive = np.asarray(naive, dtype=float)
    icb = np.asarray(icb, dtype=float)
    d_mean = float(np.nanmean(icb) - np.nanmean(naive))
    d_med = float(np.nanmedian(icb) - np.nanmedian(naive))
    return {
        "n_naive": int(np.isfinite(naive).sum()),
        "n_icb": int(np.isfinite(icb).sum()),
        "mean_naive": float(np.nanmean(naive)),
        "mean_icb": float(np.nanmean(icb)),
        "median_naive": float(np.nanmedian(naive)),
        "median_icb": float(np.nanmedian(icb)),
        "delta_mean": d_mean,
        "delta_median": d_med,
        "direction_mean": int(np.sign(d_mean)) if d_mean != 0 else 0,
    }


def tj_series(mat: pd.DataFrame) -> pd.Series:
    present = [g for g in TJ if g in mat.index]
    return mat.loc[present].mean(axis=0)


def main() -> int:
    DATA.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    TABLES.mkdir(parents=True, exist_ok=True)

    ens = ensembl_map()
    (DATA / "ensembl_ids.json").write_text(json.dumps(ens, indent=2) + "\n")
    print("Ensembl", ens)

    rows = []
    sample_rows = []

    # ----- GSE182228 KL -----
    tar_path = DATA / "GSE182228_RAW.tar"
    if not tar_path.exists():
        print("download GSE182228")
        tar_path.write_bytes(fetch("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE182nnn/GSE182228/suppl/GSE182228_RAW.tar"))
    unpack = DATA / "GSE182228_RAW"
    unpack.mkdir(exist_ok=True)
    with tarfile.open(tar_path) as tar:
        tar.extractall(unpack)
    # FPKM files: typically gene / FPKM
    gsm_treat = {
        "GSM5525443": "combo",
        "GSM5525444": "combo",
        "GSM5525445": "combo",
        "GSM5525446": "aPD1",
        "GSM5525447": "aPD1",
        "GSM5525448": "aPD1",
        "GSM5525449": "palbo",
        "GSM5525450": "palbo",
        "GSM5525451": "palbo",
        "GSM5525452": "vehicle",
        "GSM5525453": "vehicle",
        "GSM5525454": "vehicle",
    }
    expr = {}
    for p in sorted(unpack.glob("*.txt.gz")) + sorted(unpack.glob("*.txt")):
        gsm = p.name.split("_")[0]
        df = pd.read_csv(p, sep="\t")
        # find gene + FPKM columns
        cols = {c.lower(): c for c in df.columns}
        gene_col = None
        val_col = None
        for cand in ("gene_name", "gene", "symbol", "gene_id", "id"):
            if cand in cols:
                gene_col = cols[cand]
                break
        if gene_col is None:
            gene_col = df.columns[0]
        for cand in ("fpkm", "tpm", "fpkm_uncorrected", "value"):
            if cand in cols:
                val_col = cols[cand]
                break
        if val_col is None:
            # last numeric
            num = df.select_dtypes(include=[np.number])
            val_col = num.columns[-1] if len(num.columns) else df.columns[-1]
        s = pd.to_numeric(df.set_index(gene_col)[val_col], errors="coerce")
        s = s.groupby(level=0).mean()
        # if Ensembl ids
        if str(s.index[0]).startswith("ENSMUSG"):
            s.index = [i.split(".")[0] for i in s.index.astype(str)]
            s = s.rename(index={v: k for k, v in ens.items()})
            s = s.groupby(level=0).mean()
        expr[gsm] = s
    mat182 = pd.DataFrame(expr)
    focus182 = mat182.reindex(FOCUS)
    write_focus(focus182, TABLES / "GSE182228_focus.tsv")
    treat = pd.Series(gsm_treat)
    log182 = np.log2(focus182.astype(float) + 1)
    tj182 = tj_series(log182)

    def add_contrast(acc, model, genotype, contrast, naive_ids, icb_ids, mat_log, tj_s):
        for gene in FOCUS + ["TJ"]:
            if gene == "TJ":
                naive = tj_s[naive_ids].to_numpy()
                icb = tj_s[icb_ids].to_numpy()
            else:
                if gene not in mat_log.index:
                    continue
                naive = mat_log.loc[gene, naive_ids].astype(float).to_numpy()
                icb = mat_log.loc[gene, icb_ids].astype(float).to_numpy()
            st = pair_stats(naive, icb)
            rows.append(
                {
                    "accession": acc,
                    "model": model,
                    "genotype": genotype,
                    "contrast": contrast,
                    "gene": gene,
                    **st,
                }
            )
        for sid, arm, valmap in (
            [(i, "naive", mat_log[i]) for i in naive_ids]
            + [(i, "icb", mat_log[i]) for i in icb_ids]
        ):
            rec = {
                "accession": acc,
                "sample": sid,
                "arm": arm,
                "contrast": contrast,
                "TJ": float(tj_s[sid]) if sid in tj_s.index else np.nan,
            }
            for g in FOCUS:
                rec[g] = float(valmap[g]) if g in valmap.index else np.nan
            sample_rows.append(rec)

    veh = treat[treat == "vehicle"].index.tolist()
    apd = treat[treat == "aPD1"].index.tolist()
    pal = treat[treat == "palbo"].index.tolist()
    comb = treat[treat == "combo"].index.tolist()
    add_contrast("GSE182228", "LKB1-deficient LUAD s.c.", "KL", "aPD1 vs vehicle", veh, apd, log182, tj182)
    add_contrast("GSE182228", "LKB1-deficient LUAD s.c.", "KL", "aPD1+palbo vs palbo", pal, comb, log182, tj182)

    # ----- GSE114601 KP -----
    print("download GSE114601")
    raw114 = fetch("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE114nnn/GSE114601/suppl/GSE114601_counts.normalized.csv.gz")
    df114 = pd.read_csv(io.BytesIO(gzip.decompress(raw114)), index_col=0)
    focus114 = resolve_genes(df114, ens)
    write_focus(focus114, TABLES / "GSE114601_focus.tsv")
    log114 = np.log2(focus114.astype(float) + 1)
    tj114 = tj_series(log114)
    # columns are s1795 etc
    veh114 = [c for c in log114.columns if c in ("s2596", "s2617")]
    apd114 = [c for c in log114.columns if c in ("s1795", "s2521")]
    add_contrast("GSE114601", "KP GEMM lung nodules", "KP", "aPD1 vs vehicle", veh114, apd114, log114, tj114)

    # ----- GSE157880 HKP1 0 Gy -----
    print("download GSE157880")
    raw157 = fetch("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157880/suppl/GSE157880_Bulk048.txt.gz")
    df157 = pd.read_csv(io.BytesIO(gzip.decompress(raw157)), sep="\t")
    gene_col = "Gene Symbol" if "Gene Symbol" in df157.columns else "gene_name"
    num_cols = [c for c in df157.columns if c not in (
        "Chromosome", "Start", "Stop", "Strand", "Gene Symbol", "gene_biotype", "gene_name", "gene_source"
    )]
    mat157 = df157.set_index(gene_col)[num_cols]
    mat157 = mat157.groupby(level=0).mean()
    focus157 = mat157.reindex(FOCUS)
    write_focus(focus157, TABLES / "GSE157880_focus.tsv")
    log157 = np.log2(focus157.astype(float) + 1)
    tj157 = tj_series(log157)
    igg0 = [c for c in log157.columns if c.startswith("0-") and int(c.split("_")[0].split("-")[1]) <= 3]
    # titles: 0-1,0-2,0-3 = IgG 0Gy; 0-4,0-5 = PD-1 0Gy
    igg0 = [c for c in log157.columns if c.split("_")[0] in ("0-1", "0-2", "0-3")]
    pd0 = [c for c in log157.columns if c.split("_")[0] in ("0-4", "0-5")]
    add_contrast("GSE157880", "HKP1 orthotopic lung (0 Gy)", "KP", "aPD1 vs IgG (0 Gy)", igg0, pd0, log157, tj157)

    # ----- GSE169194 KPM total viable -----
    print("download GSE169194")
    raw169 = fetch(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE169nnn/GSE169194/suppl/GSE169194_annotated_log2_norm_count.txt.gz"
    )
    df169 = pd.read_csv(io.BytesIO(gzip.decompress(raw169)), sep="\t", index_col=0)
    df169.index = [str(i).split(".")[0] for i in df169.index]
    focus169 = resolve_genes(df169, ens)
    write_focus(focus169, TABLES / "GSE169194_focus.tsv")
    # already log2
    log169 = focus169.astype(float)
    tj169 = tj_series(log169)
    # columns X1740R_159_NN ; total IgG = 04,08,12 ; total A2V = 16,20,24 ; total A2V_aPD1 = 32,36,40
    def col169(*nums):
        out = []
        for n in nums:
            key = f"159_{n:02d}"
            hits = [c for c in log169.columns if key in c]
            out.extend(hits)
        return out

    igg169 = col169(4, 8, 12)
    a2v169 = col169(16, 20, 24)
    combo169 = col169(32, 36, 40)
    add_contrast(
        "GSE169194",
        "KPM total viable cells",
        "KPM",
        "A2V+aPD1 vs A2V (no PD-1 mono)",
        a2v169,
        combo169,
        log169,
        tj169,
    )

    stats = pd.DataFrame(rows)
    samples = pd.DataFrame(sample_rows)
    stats.to_csv(TABLES / "gemm_contrasts.tsv", sep="\t", index=False)
    samples.to_csv(TABLES / "gemm_samples.tsv", sep="\t", index=False)

    # panel: one primary pair per series (mean delta sign)
    primary = stats[
        stats["contrast"].isin(
            ["aPD1 vs vehicle", "aPD1 vs IgG (0 Gy)", "A2V+aPD1 vs A2V (no PD-1 mono)"]
        )
        & stats["gene"].isin(["Tacstd2", "Cldn4", "TJ"])
    ].copy()
    # GSE182228 contributes two rows if we keep both; keep aPD1 vs vehicle as primary
    primary = primary[~((primary["accession"] == "GSE182228") & (primary["contrast"] != "aPD1 vs vehicle"))]

    panel = []
    for gene, sub in primary.groupby("gene"):
        v = sub["delta_mean"]
        n_up = int((v > 0).sum())
        n_down = int((v < 0).sum())
        n_tie = int((v == 0).sum())
        n = n_up + n_down
        p_bin = float(binomtest(n_up, n, 0.5, alternative="two-sided").pvalue) if n else np.nan
        v_nz = v[v != 0]
        if len(v_nz) >= 1:
            w = wilcoxon(v_nz, alternative="two-sided", zero_method="wilcox")
            wp = float(w.pvalue)
        else:
            wp = np.nan
        panel.append(
            {
                "gene": gene,
                "n_pairs": int(len(v)),
                "n_up": n_up,
                "n_down": n_down,
                "n_tie": n_tie,
                "binomial_p_two_sided": p_bin,
                "wilcoxon_p_two_sided": wp,
                "mean_delta": float(v.mean()),
                "pairs": "; ".join(f"{r.accession} {r.genotype} Δ={r.delta_mean:.3f}" for r in sub.itertuples()),
            }
        )
    panel_df = pd.DataFrame(panel)
    panel_df.to_csv(TABLES / "gemm_panel_sign.tsv", sep="\t", index=False)

    # figures: per-study paired dots for Tacstd2/Cldn4/TJ
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 4.0), sharey=False)
    for ax, gene in zip(axes, ["Tacstd2", "Cldn4", "TJ"]):
        sub = primary[primary["gene"] == gene].reset_index(drop=True)
        for i, r in sub.iterrows():
            color = {"KL": "#C44E52", "KP": "#4C78A8", "KPM": "#54A24B"}.get(r["genotype"], "0.4")
            ax.plot([0, 1], [r["mean_naive"], r["mean_icb"]], color=color, lw=1.4)
            ax.scatter([0, 1], [r["mean_naive"], r["mean_icb"]], color=color, s=28, zorder=3)
            ax.text(1.05, r["mean_icb"], f"{r['accession']} {r['genotype']}", fontsize=7, va="center", color=color)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["naive / control", "ICB"])
        ax.set_xlim(-0.2, 1.7)
        ax.set_title(gene)
        ax.set_ylabel("mean log2(expr+1) or depositor log2")
    fig.suptitle("Public KP / KL / KPM ICB pairs (not TISMO; LLC excluded here)", y=1.02)
    fig.tight_layout()
    fig.savefig(FIG / "paired_gemm_kp_kl.png", dpi=160, bbox_inches="tight")
    plt.close(fig)

    # waterfall of GEMM deltas
    fig, ax = plt.subplots(figsize=(8.5, 4.0))
    order = ["Tacstd2", "Cldn4", "TJ"]
    x = 0
    ticks = []
    labels = []
    for gene in order:
        sub = primary[primary["gene"] == gene]
        for r in sub.itertuples():
            color = "#4C78A8" if r.delta_mean >= 0 else "#E45756"
            ax.bar(x, r.delta_mean, color=color, width=0.8)
            ticks.append(x)
            labels.append(f"{r.accession}\n{r.genotype}\n{gene}")
            x += 1
        x += 0.4
    ax.axhline(0, color="0.4", lw=0.8)
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, fontsize=7)
    ax.set_ylabel("mean ICB − mean naive")
    ax.set_title("GEMM / GEMM-derived ICB deltas")
    fig.tight_layout()
    fig.savefig(FIG / "waterfall_gemm.png", dpi=160)
    plt.close(fig)

    summary = {
        "note": (
            "TISMO has no KL or KP ICB pairs. These are independent public GEO "
            "processed matrices. LLC is not included here (LLC is not KL; LLC ICB is GSE155972 in TISMO)."
        ),
        "ensembl": ens,
        "panel": panel_df.to_dict(orient="records"),
        "contrasts": stats.to_dict(orient="records"),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=float) + "\n")
    print(panel_df.to_string(index=False))
    print(stats[stats["gene"].isin(["Tacstd2", "Cldn4", "TJ"])].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
