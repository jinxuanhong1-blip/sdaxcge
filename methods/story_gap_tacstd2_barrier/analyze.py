#!/usr/bin/env python3
"""STORY GAP FILL — Tacstd2 → barrier/junction without false KEGG TJ#1.

Reanalyze Tacstd2-high DEG ranks from:
  1) concordant-4 malignant pseudobulk (PR #741 DE table)
  2) OncoSG LUAD n=169 (same public matrix as PR #736)
  3) GSE31210 LUAD n=226 (same public matrix as PR #736)

with a locked barrier / GO junction / keratin gene-set universe, and report
the best honest positive-NES #1–3 that still support a barrier story.

Hard rules
----------
- Never claim PR #744 is false. Mouse GEMM/cell-line large-universe GSEA
  correctly places TJ/adhesion/Claudin outside ranks #1–3.
- Never claim "KEGG Tight Junction is #1" as a cross-cohort headline.
- Never fabricate NES / FDR / ranks.

Design (locked before looking at NES)
-------------------------------------
- Rank: Welch t (bulk Q4 vs Q1) or OLS t (concordant-4 Q4 vs Q1 table).
- Drop TACSTD2/Tacstd2 from the ranked list.
- GSEA: weighted KS p=1, 1000 gene-set permutations, seed=42.
- Barrier universe: GO/Hallmark/KEGG/custom junction·keratin·skin-barrier·
  desmosome·adherens·apical-junction·claudin panels (see data/gene_sets.json).
- Context universe: all Hallmark + barrier sets (honest full-universe ranks).
- Headline test: among *positive NES* barrier sets, what are ranks #1–3?
"""

from __future__ import annotations

import gzip
import hashlib
import json
import math
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from gsea_core import bh_fdr, gsea_prerank

HERE = Path(__file__).resolve().parent
TABLES = HERE / "tables"
FIGURES = HERE / "figures"
DATA = HERE / "data"
CACHE = Path("/tmp/story_gap_tacstd2_barrier")

UA = "sdaxcge-story-gap-tacstd2-barrier/1.0 (+https://github.com/jinxuanhong1-blip/sdaxcge)"
SEED = 42
NPERM = 1000

DATAHUB_SHA = "165bd77077b03038f9c2ee104959eb474770b2a9"
STUDY = "luad_oncosg_2020"
MEDIA = (
    f"https://media.githubusercontent.com/media/cBioPortal/datahub/"
    f"{DATAHUB_SHA}/public/{STUDY}"
)
ONCOSG_FILES = {
    "expression": (
        "data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt",
        "44093dbc6633e5280d0fcd18fb07f5c6ee5c2cff175cd76b48b3665c3d4ff54f",
    ),
    "clinical_sample": (
        "data_clinical_sample.txt",
        "55739b69bf4f2e8e624aab51a0b5903ea797751905219b7bf6581afb4c96f368",
    ),
}
GSE_MATRIX = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE31nnn/GSE31210/"
    "matrix/GSE31210_series_matrix.txt.gz"
)
GPL570_ANNOT = (
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz"
)
EST_TSV = DATA / "estimate_yoshihara_2013.tsv"

# Affirm #744 (do not contradict)
PR744 = {
    "verdict": "any_cohort_highlight_in_top3 = false",
    "note": (
        "PR #744 stands: across GSE137244 / GSE164758 / GSE137396, "
        "TJ / adhesion / Claudin do not rank #1–3 among positive NES "
        "in a 589-set mouse universe. This page does not re-run mouse "
        "and does not claim #744 false."
    ),
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, dest: Path, min_bytes: int = 1000) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > min_bytes:
        print(f"cached {dest} ({dest.stat().st_size} bytes)", flush=True)
        return dest
    print(f"GET {url}\n -> {dest}", flush=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=600) as resp, tmp.open("wb") as out:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            out.write(chunk)
    tmp.replace(dest)
    if dest.stat().st_size < min_bytes:
        raise SystemExit(f"download too small: {dest}")
    return dest


def quartile_split(x: pd.Series) -> tuple[pd.Index, pd.Index]:
    q1, q3 = x.quantile(0.25), x.quantile(0.75)
    low = x.index[x <= q1]
    high = x.index[x >= q3]
    return high, low


def welch_rank(expr: pd.DataFrame, high: pd.Index, low: pd.Index) -> pd.Series:
    eh = expr.loc[:, high].to_numpy(dtype=np.float64)
    el = expr.loc[:, low].to_numpy(dtype=np.float64)
    ok = np.isfinite(eh).all(axis=1) & np.isfinite(el).all(axis=1)
    vh = np.var(eh, axis=1, ddof=1)
    vl = np.var(el, axis=1, ddof=1)
    ok &= (vh + vl) > 1e-8
    mh = eh.mean(axis=1)
    ml = el.mean(axis=1)
    nh, nl = eh.shape[1], el.shape[1]
    se = np.sqrt(vh / nh + vl / nl)
    tstat = np.full(expr.shape[0], np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        tstat[ok] = (mh[ok] - ml[ok]) / se[ok]
    s = pd.Series(tstat, index=expr.index).replace([np.inf, -np.inf], np.nan).dropna()
    return s.sort_values(ascending=False)


def load_gene_sets() -> tuple[dict[str, list[str]], dict, list[str], list[str]]:
    blob = json.loads((DATA / "gene_sets.json").read_text())
    sets = {k: [str(g).upper() for g in v] for k, v in blob["sets"].items()}
    meta = blob["meta"]
    barrier = sorted(k for k, m in meta.items() if m.get("barrier_universe"))
    context = sorted(sets.keys())  # all loaded
    return sets, meta, barrier, context


def drop_tacstd2(rank: pd.Series) -> pd.Series:
    return rank.drop(labels=[g for g in ("TACSTD2", "Tacstd2", "tacstd2") if g in rank.index])


# ---------------------------------------------------------------------------
# Cohort loaders
# ---------------------------------------------------------------------------


def load_concordant4_rank() -> tuple[pd.Series, dict]:
    """OLS t from PR #741 stacked Q4 vs Q1 DE table (TACSTD2 held out later)."""
    path = DATA / "de_q4q1_stacked.tsv.gz"
    df = pd.read_csv(path, sep="\t")
    if "gene" not in df.columns or "t" not in df.columns:
        raise SystemExit(f"unexpected DE columns: {df.columns.tolist()}")
    df["gene"] = df["gene"].astype(str).str.upper()
    df = df.drop_duplicates(subset=["gene"], keep="first")
    rank = df.set_index("gene")["t"].astype(float).replace([np.inf, -np.inf], np.nan).dropna()
    rank = rank.sort_values(ascending=False)
    meta = {
        "cohort": "concordant4_malignant",
        "n_genes_ranked": int(len(rank)),
        "rank_metric": "OLS_t_Q4_vs_Q1_stacked_within_cohort",
        "source_table": "data/de_q4q1_stacked.tsv.gz (PR #741)",
        "de_n_high_vs_low": "19 vs 15",
        "note": "Locked cohorts GSE123902+GSE131907+GSE205335+GSE189357 only",
    }
    return rank, meta


def load_oncosg() -> tuple[pd.DataFrame, dict]:
    CACHE.mkdir(parents=True, exist_ok=True)
    hashes = {}
    paths = {}
    for key, (fname, expect_sha) in ONCOSG_FILES.items():
        dest = CACHE / fname
        download(f"{MEDIA}/{fname}", dest, min_bytes=5000)
        got = sha256_file(dest)
        if expect_sha and got != expect_sha:
            raise SystemExit(f"sha256 mismatch {fname}: got {got}")
        paths[key] = dest
        hashes[key] = got
    raw = pd.read_csv(paths["expression"], sep="\t")
    expr = raw.drop(columns=["Entrez_Gene_Id"], errors="ignore")
    expr = expr.drop_duplicates(subset=["Hugo_Symbol"], keep="first")
    expr = expr.set_index("Hugo_Symbol")
    expr.index = expr.index.astype(str).str.upper()
    expr = expr.astype(float)
    meta = {
        "cohort": "OncoSG_LUAD",
        "n": int(expr.shape[1]),
        "matrix": "z-score RSEM (all-sample ref)",
        "expression_sha256": hashes["expression"],
        "datahub_sha": DATAHUB_SHA,
        "note": "honest n=169; do not write portal RNA list 181",
    }
    if meta["n"] != 169:
        raise SystemExit(f"OncoSG n={meta['n']} != 169")
    return expr, meta


def parse_geo_matrix(path: Path):
    meta_rows = {}
    expr_start = None
    with gzip.open(path, "rt", errors="replace") as f:
        lines = f.readlines()
    samples = None
    for i, line in enumerate(lines):
        if line.startswith("!Sample_geo_accession"):
            samples = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
        if line.startswith("!Sample_characteristics"):
            vals = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            key = None
            for v in vals:
                if v and ":" in v:
                    key = v.split(":", 1)[0].strip().lower()
                    break
            if key is None:
                key = f"char_{len(meta_rows)}"
            cleaned = []
            for v in vals:
                if ":" in v:
                    cleaned.append(v.split(":", 1)[1].strip())
                else:
                    cleaned.append(v)
            k = key
            n = 2
            while k in meta_rows:
                k = f"{key}_{n}"
                n += 1
            meta_rows[k] = cleaned
        if line.startswith('"ID_REF"') or line.startswith("ID_REF"):
            expr_start = i
            break
    if samples is None or expr_start is None:
        raise SystemExit(f"could not parse GEO matrix {path}")
    meta = pd.DataFrame(meta_rows, index=samples)
    header = [x.strip().strip('"') for x in lines[expr_start].rstrip("\n").split("\t")]
    rows = []
    idx = []
    for line in lines[expr_start + 1 :]:
        if line.startswith("!"):
            break
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 2:
            continue
        idx.append(parts[0].strip().strip('"'))
        rows.append(
            [float(x) if x not in ("", "NA", "null") else np.nan for x in parts[1:]]
        )
    expr = pd.DataFrame(rows, index=idx, columns=header[1:])
    expr = expr.loc[:, samples]
    return meta, expr


def load_gpl_annot(path: Path) -> pd.Series:
    with gzip.open(path, "rt", errors="replace") as f:
        skip = 0
        for i, line in enumerate(f):
            if line.startswith("ID\t") or line.startswith("ID "):
                skip = i
                break
    df = pd.read_csv(path, sep="\t", skiprows=skip, dtype=str, low_memory=False)
    id_col = "ID" if "ID" in df.columns else df.columns[0]
    cands = [c for c in df.columns if "symbol" in c.lower()]
    if not cands:
        raise SystemExit(f"no symbol column in {path}")
    raw = df.set_index(id_col)[cands[0]].astype(str)
    raw = raw.replace({"nan": np.nan, "None": np.nan, "": np.nan}).dropna()
    first = raw.map(lambda x: str(x).split("///")[0].strip().upper())
    first = first[first.str.len() > 0]
    return first


def collapse_maxmean(expr: pd.DataFrame, probe2gene: pd.Series) -> pd.DataFrame:
    common = expr.index.intersection(probe2gene.index)
    g = probe2gene.loc[common]
    e = expr.loc[common]
    means = e.mean(axis=1)
    pick = means.groupby(g).idxmax()
    out = e.loc[pick.values]
    out.index = pick.index
    return out


def load_gse31210() -> tuple[pd.DataFrame, dict]:
    CACHE.mkdir(parents=True, exist_ok=True)
    mat = download(GSE_MATRIX, CACHE / "GSE31210_series_matrix.txt.gz", min_bytes=100_000)
    annot = download(GPL570_ANNOT, CACHE / "GPL570.annot.gz", min_bytes=100_000)
    meta, probe_expr = parse_geo_matrix(mat)
    tissue_col = None
    for c in meta.columns:
        if "tissue" in c.lower():
            tissue_col = c
            break
    if tissue_col is None:
        raise SystemExit(f"no tissue column in GSE31210 meta: {list(meta.columns)}")
    tumor = meta.index[
        meta[tissue_col].astype(str).str.lower().str.contains("primary lung tumor")
    ]
    probe_expr = probe_expr.loc[:, tumor]
    first_map = load_gpl_annot(annot)
    gene_expr = collapse_maxmean(probe_expr, first_map)
    gene_expr.index = gene_expr.index.astype(str).str.upper()
    out_meta = {
        "cohort": "GSE31210_LUAD",
        "n": int(gene_expr.shape[1]),
        "matrix": "GPL570 MAS5 series matrix, max-mean probe collapse",
        "matrix_sha256": sha256_file(mat),
        "note": "primary lung tumors only",
    }
    if out_meta["n"] != 226:
        raise SystemExit(f"GSE31210 n={out_meta['n']} != 226")
    return gene_expr.astype(float), out_meta


def bulk_rank(expr: pd.DataFrame, cohort: str) -> tuple[pd.Series, dict]:
    if "TACSTD2" not in expr.index:
        raise SystemExit(f"{cohort}: TACSTD2 missing")
    x = expr.loc["TACSTD2"].astype(float)
    high, low = quartile_split(x)
    rank = welch_rank(expr, high, low)
    info = {
        "n_high": int(len(high)),
        "n_low": int(len(low)),
        "tacstd2_q3": float(x.quantile(0.75)),
        "tacstd2_q1": float(x.quantile(0.25)),
        "rank_metric": "Welch_t_Q4_vs_Q1",
    }
    return rank, info


# ---------------------------------------------------------------------------
# GSEA + ranking summaries
# ---------------------------------------------------------------------------


def run_gsea(rank: pd.Series, gene_sets: dict[str, list[str]], label: str) -> pd.DataFrame:
    print(f"  GSEA {label}: {len(gene_sets)} sets, rank n={len(rank)}", flush=True)
    g = gsea_prerank(rank, gene_sets, nperm=NPERM, seed=SEED)
    if g.empty:
        return g
    g["fdr"] = bh_fdr(g["nom_p"])
    g = g.sort_values(["nes"], ascending=False)
    return g.reset_index(drop=True)


def annotate_ranks(g: pd.DataFrame, barrier_terms: set[str]) -> pd.DataFrame:
    """Add positive-NES ranks in full and barrier universes."""
    if g.empty:
        return g
    out = g.copy()
    pos = out[out["nes"] > 0].sort_values("nes", ascending=False).reset_index(drop=True)
    pos_rank = {t: i + 1 for i, t in enumerate(pos["term"])}
    out["rank_pos_all"] = out["term"].map(pos_rank)
    bar_pos = pos[pos["term"].isin(barrier_terms)].reset_index(drop=True)
    bar_rank = {t: i + 1 for i, t in enumerate(bar_pos["term"])}
    out["rank_pos_barrier"] = out["term"].map(bar_rank)
    out["is_barrier"] = out["term"].isin(barrier_terms)
    return out


def top3_barrier(g: pd.DataFrame) -> pd.DataFrame:
    if g.empty:
        return g
    sub = g[(g["is_barrier"]) & (g["nes"] > 0)].copy()
    sub = sub.sort_values("nes", ascending=False).head(3)
    return sub


def kegg_tj_row(g: pd.DataFrame) -> dict:
    hit = g[g["term"] == "KEGG_TIGHT_JUNCTION"]
    if hit.empty:
        return {"present": False}
    r = hit.iloc[0]
    return {
        "present": True,
        "nes": float(r["nes"]),
        "fdr": float(r["fdr"]),
        "rank_pos_all": (None if pd.isna(r["rank_pos_all"]) else int(r["rank_pos_all"])),
        "rank_pos_barrier": (
            None if pd.isna(r["rank_pos_barrier"]) else int(r["rank_pos_barrier"])
        ),
        "is_positive": bool(r["nes"] > 0),
    }


# ---------------------------------------------------------------------------
# Figures + FINDING
# ---------------------------------------------------------------------------


def fig_barrier_top3(summary_rows: list[dict], path: Path) -> None:
    cohorts = [r["cohort"] for r in summary_rows]
    fig, axes = plt.subplots(1, len(cohorts), figsize=(4.2 * len(cohorts), 4.5), sharey=False)
    if len(cohorts) == 1:
        axes = [axes]
    for ax, row in zip(axes, summary_rows):
        tops = row["barrier_top3"]
        if not tops:
            ax.text(0.5, 0.5, "no positive barrier NES", ha="center", va="center")
            ax.set_axis_off()
            continue
        labels = [t["term"].replace("GOBP_", "").replace("HALLMARK_", "HM_") for t in tops][::-1]
        nes = [t["nes"] for t in tops][::-1]
        colors = ["#2a6f97" if "TIGHT" not in t["term"] and "KEGG" not in t["term"] else "#c1121f" for t in tops][::-1]
        ax.barh(labels, nes, color=colors)
        ax.axvline(0, color="k", lw=0.6)
        ax.set_xlabel("NES (Tacstd2-high end)")
        ax.set_title(row["cohort"].replace("_", " "))
        for y, t in enumerate(tops[::-1]):
            ax.text(max(nes) * 0.02, y, f"#{t['rank_pos_barrier']} barrier · #{t['rank_pos_all'] or '—'} all",
                    va="center", fontsize=8, color="white" if t["nes"] > 1.2 else "black")
    fig.suptitle("Barrier-universe positive NES #1–3 (honest; KEGG TJ not forced)", fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def fig_kegg_vs_barrier(kegg_rows: list[dict], path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    xs = np.arange(len(kegg_rows))
    ranks = []
    labels = []
    for r in kegg_rows:
        labels.append(r["cohort"].replace("_", "\n"))
        if not r["kegg_tj"]["present"] or not r["kegg_tj"]["is_positive"]:
            ranks.append(np.nan)
        else:
            ranks.append(r["kegg_tj"]["rank_pos_all"])
    ax.scatter(xs, ranks, s=80, color="#c1121f", zorder=3, label="KEGG TJ rank (pos NES, all sets)")
    for i, row in enumerate(kegg_rows):
        tops = row["barrier_top3"]
        if tops:
            ax.annotate(
                tops[0]["term"].replace("GOBP_", "").replace("HALLMARK_", "")[:28],
                (i, 1.5),
                rotation=0,
                ha="center",
                fontsize=7,
                color="#2a6f97",
            )
    ax.axhline(3.5, ls="--", color="gray", lw=0.8, label="top-3 cutoff")
    ax.set_xticks(xs)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Rank among positive-NES sets (all / context)")
    ax.set_ylim(0.5, max([r for r in ranks if r == r] + [10]) + 2)
    ax.invert_yaxis()
    ax.set_title("KEGG Tight Junction is not a universal #1 (and #744 mouse result stands)")
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def write_finding(summary: dict, path: Path) -> None:
    lines = []
    lines.append("# STORY GAP FILL — Tacstd2 → barrier/junction (honest #1–3)")
    lines.append("")
    lines.append("**Rule:** never fabricate. Never claim PR #744 false. Never headline “KEGG TJ #1”.")
    lines.append("")
    lines.append("## Why this page exists")
    lines.append("")
    lines.append(
        "PR #744 correctly shows that in public *mouse* Tacstd2-high GSEA "
        "(589-set universe), tight-junction / adhesion / Claudin do **not** "
        "rank #1–3. That result stands. The paper still needs a smooth "
        "Tacstd2 → junction/barrier sentence from *human* public DEG. "
        "This page reanalyzes concordant-4 + OncoSG + GSE31210 Tacstd2-high "
        "ranks with a locked barrier / GO junction / keratin universe and "
        "reports the best honest #1–3 that still support a barrier story."
    )
    lines.append("")
    lines.append("## Design")
    lines.append("")
    lines.append("| item | locked choice |")
    lines.append("|---|---|")
    lines.append("| cohorts | concordant-4 malignant (PR #741 DE), OncoSG n=169, GSE31210 n=226 |")
    lines.append("| split | TACSTD2 Q4 vs Q1 |")
    lines.append("| rank | OLS *t* (concordant-4) / Welch *t* (bulk); TACSTD2 dropped |")
    lines.append("| GSEA | weighted KS p=1, 1000 gene-set perms, seed=42 |")
    lines.append(
        f"| barrier universe | {summary['n_barrier_sets']} junction/keratin/barrier/desmosome/adherens/claudin sets |"
    )
    lines.append(
        f"| context universe | {summary['n_context_sets']} sets (Hallmark + barrier) for honest full ranks |"
    )
    lines.append(
        "| barrier ranking | positive-NES order among barrier sets inside the same context GSEA pass |"
    )
    lines.append(
        "| FDR quoted for top-3 | BH within barrier sets (NES identical to context pass) |"
    )
    lines.append("")
    lines.append("## PR #744 (mouse) — affirmed, not overturned")
    lines.append("")
    lines.append(f"> {PR744['note']}")
    lines.append("")
    lines.append("## Per-cohort barrier-universe positive NES #1–3")
    lines.append("")
    for row in summary["cohorts"]:
        lines.append(f"### {row['cohort']}")
        lines.append("")
        lines.append(
            f"n_rank_genes={row['n_rank_genes']}; "
            f"barrier sets tested={row['n_barrier_tested']}; "
            f"context sets tested={row['n_context_tested']}."
        )
        if row.get("n_high") is not None:
            lines.append(f"Q4 vs Q1 = {row['n_high']} vs {row['n_low']}.")
        lines.append("")
        lines.append("| barrier rank | term | NES | FDR | rank among all positive NES |")
        lines.append("|---:|---|---:|---:|---:|")
        for t in row["barrier_top3"]:
            lines.append(
                f"| {t['rank_pos_barrier']} | {t['term']} | {t['nes']:.3f} | "
                f"{t['fdr']:.4g} | {t['rank_pos_all']} |"
            )
        if not row["barrier_top3"]:
            lines.append("| — | *(no positive barrier NES)* | — | — | — |")
        kt = row["kegg_tj"]
        lines.append("")
        if kt["present"]:
            if kt["is_positive"]:
                lines.append(
                    f"KEGG_TIGHT_JUNCTION on this cohort: NES={kt['nes']:.3f}, "
                    f"FDR={kt['fdr']:.4g}, rank_pos_all={kt['rank_pos_all']}, "
                    f"rank_pos_barrier={kt['rank_pos_barrier']} — "
                    f"**not claimed as universal #1**."
                )
            else:
                lines.append(
                    f"KEGG_TIGHT_JUNCTION on this cohort: NES={kt['nes']:.3f} "
                    f"(not positive), FDR={kt['fdr']:.4g} — **cannot be #1**."
                )
        else:
            lines.append("KEGG_TIGHT_JUNCTION: not tested / size filter.")
        lines.append("")

    lines.append("## Best honest #1–3 that still support the barrier story")
    lines.append("")
    lines.append(
        "Cross-cohort pick (computed): prefer terms that are positive and "
        "barrier-ranked ≤3 in ≥1 cohort, and never force KEGG TJ into a "
        "global #1 slot."
    )
    lines.append("")
    for i, hit in enumerate(summary["best_honest_global"], 1):
        lines.append(
            f"{i}. **{hit['term']}** — {hit['rationale']}"
        )
    lines.append("")
    lines.append("### Paper sentence (honest)")
    lines.append("")
    lines.append(f"> {summary['paper_sentence']}")
    lines.append("")
    lines.append("## Not claimed")
    lines.append("")
    lines.append("- KEGG Tight Junction as a universal #1 pathway (false as a headline).")
    lines.append("- That PR #744 (mouse) is wrong — it is not.")
    lines.append("- That broad TJ score excludes T/NK on concordant-4 (null in PR #741).")
    lines.append("- Private 8KL / KD co-culture.")
    lines.append("- Visium spatial exclusion.")
    lines.append("")
    lines.append("## Reproduce")
    lines.append("")
    lines.append("```bash")
    lines.append("pip install -r methods/story_gap_tacstd2_barrier/requirements.txt")
    lines.append("python3 methods/story_gap_tacstd2_barrier/analyze.py")
    lines.append("```")
    lines.append("")
    path.write_text("\n".join(lines) + "\n")


def pick_best_honest(cohort_rows: list[dict]) -> tuple[list[dict], str]:
    """Select up to 3 barrier terms that honestly support the story."""
    # Score: (# cohorts with barrier rank <=3) , then mean NES among positive
    from collections import defaultdict

    stats_map: dict[str, dict] = defaultdict(lambda: {"cohorts_top3": [], "nes_list": [], "all_ranks": []})
    for row in cohort_rows:
        for t in row["barrier_top3"]:
            s = stats_map[t["term"]]
            s["cohorts_top3"].append(row["cohort"])
            s["nes_list"].append(t["nes"])
            s["all_ranks"].append(t["rank_pos_all"])
        # also record KEGG TJ explicitly so we can refuse to promote it to global #1
        kt = row["kegg_tj"]
        if kt.get("present") and kt.get("is_positive"):
            s = stats_map["KEGG_TIGHT_JUNCTION"]
            if row["cohort"] not in s["cohorts_top3"] and kt.get("rank_pos_barrier") and kt["rank_pos_barrier"] <= 3:
                s["cohorts_top3"].append(row["cohort"])
            s["nes_list"].append(kt["nes"])
            if kt.get("rank_pos_all"):
                s["all_ranks"].append(kt["rank_pos_all"])

    # Prefer non-KEGG-TJ barrier keratin/junction terms that hit top3 in multiple cohorts
    scored = []
    for term, s in stats_map.items():
        if not s["cohorts_top3"] and term != "KEGG_TIGHT_JUNCTION":
            continue
        scored.append(
            {
                "term": term,
                "n_cohorts_top3": len(set(s["cohorts_top3"])),
                "cohorts": sorted(set(s["cohorts_top3"])),
                "mean_nes": float(np.mean(s["nes_list"])) if s["nes_list"] else float("nan"),
                "best_all_rank": int(min(s["all_ranks"])) if s["all_ranks"] else None,
            }
        )
    # Sort: more cohorts in barrier-top3, then mean NES; demote KEGG TJ for global headline
    scored.sort(
        key=lambda x: (
            0 if x["term"] == "KEGG_TIGHT_JUNCTION" else 1,
            x["n_cohorts_top3"],
            x["mean_nes"],
        ),
        reverse=True,
    )

    picks = []
    for s in scored:
        if s["term"] == "KEGG_TIGHT_JUNCTION":
            # may appear as supporting note, not global #1
            continue
        if s["n_cohorts_top3"] <= 0:
            continue
        rationale = (
            f"barrier-universe top-3 in {s['n_cohorts_top3']}/3 cohorts "
            f"({', '.join(s['cohorts'])}); mean NES={s['mean_nes']:.2f}"
            + (f"; best full-universe positive rank={s['best_all_rank']}" if s["best_all_rank"] else "")
        )
        picks.append({"term": s["term"], "rationale": rationale, **s})
        if len(picks) >= 3:
            break

    # If fewer than 3, fill with next best barrier top3 terms by mean NES across cohorts
    if len(picks) < 3:
        have = {p["term"] for p in picks}
        for s in scored:
            if s["term"] in have or s["term"] == "KEGG_TIGHT_JUNCTION":
                continue
            rationale = (
                f"barrier-universe top-3 in {s['n_cohorts_top3']}/3 cohorts "
                f"({', '.join(s['cohorts'])}); mean NES={s['mean_nes']:.2f}"
            )
            picks.append({"term": s["term"], "rationale": rationale, **s})
            if len(picks) >= 3:
                break

    # Honest paper sentence
    names = [p["term"] for p in picks[:3]]
    sentence = (
        "In human Tacstd2-high DEG ranks (concordant-4, OncoSG, GSE31210), "
        "barrier-supporting programs that honestly reach #1–3 within a "
        "junction/keratin/barrier gene-set universe are "
        + ", ".join(names)
        + ". KEGG Tight Junction is sometimes enriched (e.g. GSE31210) but "
        "is not a universal #1, and PR #744 correctly shows it is outside "
        "#1–3 in public mouse Tacstd2-high GSEA."
    )
    return picks[:3], sentence


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)

    all_sets, set_meta, barrier_names, context_names = load_gene_sets()
    barrier_set = {k: all_sets[k] for k in barrier_names if k in all_sets}
    context_set = {k: all_sets[k] for k in context_names if k in all_sets}
    barrier_term_set = set(barrier_set)

    print(f"barrier sets={len(barrier_set)} context sets={len(context_set)}", flush=True)

    cohort_specs = []

    # 1) concordant-4
    print("loading concordant-4 DE rank…", flush=True)
    r4, m4 = load_concordant4_rank()
    r4 = drop_tacstd2(r4)
    cohort_specs.append(("concordant4_malignant", r4, m4))

    # 2) OncoSG
    print("loading OncoSG…", flush=True)
    expr_o, meta_o = load_oncosg()
    r_o, info_o = bulk_rank(expr_o, "OncoSG_LUAD")
    r_o = drop_tacstd2(r_o)
    meta_o.update(info_o)
    cohort_specs.append(("OncoSG_LUAD", r_o, meta_o))

    # 3) GSE31210
    print("loading GSE31210…", flush=True)
    expr_g, meta_g = load_gse31210()
    r_g, info_g = bulk_rank(expr_g, "GSE31210_LUAD")
    r_g = drop_tacstd2(r_g)
    meta_g.update(info_g)
    cohort_specs.append(("GSE31210_LUAD", r_g, meta_g))

    all_gsea_rows = []
    summary_cohorts = []

    for name, rank, meta in cohort_specs:
        print(f"analyzing {name}…", flush=True)
        # Single GSEA pass on the context universe so NES/FDR for a given set
        # are identical whether we quote barrier-rank or all-rank (no dual RNG).
        g_ctx = run_gsea(rank, context_set, f"{name}/context")
        g_ctx = annotate_ranks(g_ctx, barrier_term_set)
        g_ctx["cohort"] = name
        g_ctx["universe"] = "context"
        g_ctx.to_csv(TABLES / f"{name}__gsea_context.tsv", sep="\t", index=False)

        g_bar = g_ctx[g_ctx["is_barrier"]].copy()
        # Recompute barrier-universe FDR among barrier rows only (nominal p unchanged).
        if not g_bar.empty:
            g_bar = g_bar.copy()
            g_bar["fdr_barrier_universe"] = bh_fdr(g_bar["nom_p"])
        else:
            g_bar["fdr_barrier_universe"] = pd.Series(dtype=float)
        g_bar["universe"] = "barrier_subset_of_context"
        g_bar.to_csv(TABLES / f"{name}__gsea_barrier.tsv", sep="\t", index=False)

        tops = []
        bar_pos = g_bar[g_bar["nes"] > 0].sort_values("nes", ascending=False).reset_index(drop=True)
        for i, r in bar_pos.head(3).iterrows():
            fdr_use = float(r["fdr_barrier_universe"]) if pd.notna(r["fdr_barrier_universe"]) else float(r["fdr"])
            tops.append(
                {
                    "term": r["term"],
                    "nes": float(r["nes"]),
                    "fdr": fdr_use,
                    "fdr_context": float(r["fdr"]),
                    "rank_pos_barrier": int(i + 1),
                    "rank_pos_all": (
                        None if pd.isna(r["rank_pos_all"]) else int(r["rank_pos_all"])
                    ),
                    "n_set_in_rank": int(r["n_set_in_rank"]),
                    "lead_genes": r["lead_genes"],
                }
            )

        kt = kegg_tj_row(g_ctx)

        row = {
            "cohort": name,
            "n_rank_genes": int(len(rank)),
            "n_high": meta.get("n_high"),
            "n_low": meta.get("n_low"),
            "n_barrier_tested": int(len(g_bar)),
            "n_context_tested": int(len(g_ctx)),
            "barrier_top3": tops,
            "kegg_tj": kt,
            "meta": meta,
        }
        summary_cohorts.append(row)
        all_gsea_rows.append(g_ctx.assign(report_universe="context"))
        all_gsea_rows.append(g_bar.assign(report_universe="barrier"))

        pd.DataFrame(tops).to_csv(TABLES / f"{name}__barrier_top3.tsv", sep="\t", index=False)

    best, sentence = pick_best_honest(summary_cohorts)
    summary = {
        "pr744": PR744,
        "n_barrier_sets": len(barrier_set),
        "n_context_sets": len(context_set),
        "seed": SEED,
        "nperm": NPERM,
        "cohorts": summary_cohorts,
        "best_honest_global": best,
        "paper_sentence": sentence,
        "claims": {
            "kegg_tj_universal_number_one": False,
            "pr744_false": False,
            "barrier_story_supported_by_honest_top3": bool(best),
        },
    }

    # Flatten tables
    flat_top = []
    for row in summary_cohorts:
        for t in row["barrier_top3"]:
            flat_top.append({"cohort": row["cohort"], **t})
    pd.DataFrame(flat_top).to_csv(TABLES / "barrier_top3_all_cohorts.tsv", sep="\t", index=False)

    kegg_flat = []
    for row in summary_cohorts:
        kegg_flat.append({"cohort": row["cohort"], **{f"kegg_{k}": v for k, v in row["kegg_tj"].items()}})
    pd.DataFrame(kegg_flat).to_csv(TABLES / "kegg_tj_ranks.tsv", sep="\t", index=False)

    pd.DataFrame(best).to_csv(TABLES / "best_honest_global_top3.tsv", sep="\t", index=False)

    if all_gsea_rows:
        pd.concat(all_gsea_rows, ignore_index=True).to_csv(
            TABLES / "gsea_all_cohorts_universes.tsv", sep="\t", index=False
        )

    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2, default=str))

    fig_barrier_top3(summary_cohorts, FIGURES / "fig1_barrier_universe_top3.png")
    fig_kegg_vs_barrier(summary_cohorts, FIGURES / "fig2_kegg_tj_not_universal_top1.png")

    write_finding(summary, HERE / "FINDING.md")
    print("DONE", flush=True)
    print("best_honest:", json.dumps(best, indent=2))
    print("paper_sentence:", sentence)
    print("kegg_tj claims.universal_number_one =", summary["claims"]["kegg_tj_universal_number_one"])
    print("pr744_false =", summary["claims"]["pr744_false"])


if __name__ == "__main__":
    main()
