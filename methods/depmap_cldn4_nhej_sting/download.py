#!/usr/bin/env python3
"""Download DepMap Public 24Q4 slices for CLDN4 vs NHEJ / cGAS–STING.

Source (Figshare+, not the bot-gated portal):
  DepMap, Broad (2024). DepMap 24Q4 Public.
  https://doi.org/10.25452/figshare.plus.27993248.v1

Integrated Chronos (CRISPRGeneEffect) is the same gene-effect matrix used by
the earlier DepMap pages in this repo. PRKDC is absent from that matrix.
README.txt calls the CD screens Humagne-CD Cas12; ScreenGeneEffect keeps
those screens, and PRKDC is populated only there. This script keeps that
screen-level slice separate from the integrated Chronos extract.

Full matrices are deleted after the slim tables and the CLDN4-versus-all
Spearman files are written.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

RELEASE = "DepMap Public 24Q4"
DOI = "10.25452/figshare.plus.27993248.v1"
UA = (
    "sdaxcge-depmap-cldn4-nhej-sting/1.0 "
    "(public DepMap extract; +https://github.com/jinxuanhong1-blip/sdaxcge)"
)

FILES = {
    "Model.csv": "https://ndownloader.figshare.com/files/51065297",
    "README.txt": "https://ndownloader.figshare.com/files/51065795",
    "ScreenSequenceMap.csv": "https://ndownloader.figshare.com/files/51065828",
    "CRISPRGeneEffect.csv": "https://ndownloader.figshare.com/files/51064667",
    "OmicsExpressionProteinCodingGenesTPMLogp1.csv": (
        "https://ndownloader.figshare.com/files/51065489"
    ),
    "ScreenGeneEffect.csv": "https://ndownloader.figshare.com/files/51065804",
}

# Classical NHEJ. PRKDC is the DNA-PKcs gene; it is not in CRISPRGeneEffect.
NHEJ_GENES = [
    "PRKDC",
    "LIG4",
    "XRCC4",
    "XRCC5",
    "XRCC6",
    "NHEJ1",
    "DCLRE1C",
    "PAXX",
]
# cGAS–STING machinery. CGAS is the HGNC symbol for cGAS. Downstream IFN
# transcripts are intentionally not in this list.
STING_GENES = ["CGAS", "STING1", "TBK1", "IRF3"]
QUERY = "CLDN4"
EXTRA_RNA = ["MKI67", "EPCAM"]
SCREEN_GENES = [QUERY] + NHEJ_GENES + STING_GENES
MIN_GW_N = 100


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    print(f"GET {url} -> {dest}", flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    n = 0
    with urllib.request.urlopen(req, timeout=600) as resp, tmp.open("wb") as out:
        while True:
            chunk = resp.read(1024 * 1024)
            if not chunk:
                break
            out.write(chunk)
            n += len(chunk)
            if n % (50 * 1024 * 1024) < 1024 * 1024:
                print(f"  {dest.name}: {n / 1e6:.0f} MB", flush=True)
    tmp.replace(dest)
    print(f"  wrote {dest} ({dest.stat().st_size} bytes)", flush=True)


def ensure(url: str, dest: Path, min_bytes: int) -> None:
    if dest.exists() and dest.stat().st_size >= min_bytes:
        print(f"OK exists {dest} ({dest.stat().st_size} bytes)", flush=True)
        return
    download(url, dest)
    if dest.stat().st_size < min_bytes:
        raise SystemExit(f"{dest} too small: {dest.stat().st_size}")


def symbol_of(col: str) -> str:
    return col.split(" (")[0].strip()


def read_gene_matrix(path: Path, wanted: set[str] | None = None) -> pd.DataFrame:
    """Read a DepMap gene matrix. First column is the row id.

    Column labels are rewritten to HGNC symbols (text before ' (').
    If wanted is set, only those genes plus the id are read.
    """
    header = pd.read_csv(path, nrows=0)
    id_col = header.columns[0]
    usecols = None
    if wanted is not None:
        keep = [id_col]
        seen: set[str] = set()
        for col in header.columns[1:]:
            sym = symbol_of(str(col))
            if sym in wanted and sym not in seen:
                keep.append(col)
                seen.add(sym)
        missing = sorted(wanted - seen)
        print(f"{path.name}: slim {len(seen)}/{len(wanted)}; missing={missing}", flush=True)
        usecols = keep
    df = pd.read_csv(path, usecols=usecols, low_memory=False)
    df = df.rename(columns={id_col: "row_id"})
    rename = {}
    seen = set()
    for col in df.columns[1:]:
        sym = symbol_of(str(col))
        if sym in seen:
            continue
        rename[col] = sym
        seen.add(sym)
    df = df.rename(columns=rename)
    # Drop duplicate-symbol columns that were not renamed.
    df = df.loc[:, ~df.columns.duplicated()]
    for col in df.columns[1:]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["row_id"] = df["row_id"].astype(str)
    return df


def spearman_vs_query(df: pd.DataFrame, query: str, min_n: int) -> pd.DataFrame:
    """Pairwise Spearman of one query column against every other numeric column.

    Rows with a missing query value are dropped first. Partner columns still
    use their own pairwise-complete subset. CLDN4 against itself is omitted
    by the caller when the column name equals the query.
    """
    if query not in df.columns:
        raise SystemExit(f"{query} not in matrix")
    q = df[query].to_numpy(dtype=float)
    base = np.isfinite(q)
    qv = q[base]
    if len(qv) < min_n:
        raise SystemExit(f"{query} has only {len(qv)} finite rows")
    rq = rankdata(qv).astype(float)
    rq -= rq.mean()
    rq_ss = float(np.dot(rq, rq))
    genes = []
    ns = []
    rhos = []
    for gene in df.columns:
        if gene in ("row_id", query):
            continue
        y_all = df[gene].to_numpy(dtype=float)
        y = y_all[base]
        m = np.isfinite(y)
        n = int(m.sum())
        if n < min_n:
            continue
        yy = y[m]
        if float(np.std(yy)) == 0.0:
            continue
        ry = rankdata(yy).astype(float)
        ry -= ry.mean()
        if m.all():
            denom = np.sqrt(rq_ss * float(np.dot(ry, ry)))
            rho = float(np.dot(rq, ry) / denom) if denom else np.nan
        else:
            rx = rankdata(qv[m]).astype(float)
            rx -= rx.mean()
            denom = np.sqrt(float(np.dot(rx, rx)) * float(np.dot(ry, ry)))
            rho = float(np.dot(rx, ry) / denom) if denom else np.nan
        if not np.isfinite(rho):
            continue
        genes.append(gene)
        ns.append(n)
        rhos.append(rho)
    out = pd.DataFrame({"gene": genes, "n": ns, "rho": rhos})
    return out.sort_values("rho", ascending=False).reset_index(drop=True)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="data/depmap_cldn4_nhej_sting")
    args = ap.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    model_path = outdir / "Model.csv"
    readme_path = outdir / "README.txt"
    screen_map_path = outdir / "ScreenSequenceMap.csv"
    crispr_path = outdir / "CRISPRGeneEffect.csv"
    expr_path = outdir / "OmicsExpressionProteinCodingGenesTPMLogp1.csv"
    screen_path = outdir / "ScreenGeneEffect.csv"

    ensure(FILES["Model.csv"], model_path, 100_000)
    ensure(FILES["README.txt"], readme_path, 5_000)
    ensure(FILES["ScreenSequenceMap.csv"], screen_map_path, 50_000)

    model = pd.read_csv(model_path)
    if "ModelID" not in model.columns:
        raise SystemExit(f"Model.csv columns: {list(model.columns)[:12]}")
    lung_ids = set(
        model.loc[
            (model["OncotreeLineage"] == "Lung") & (model["ModelType"] == "Cell Line"),
            "ModelID",
        ].astype(str)
    )
    print(f"lung cell lines in Model.csv: {len(lung_ids)}", flush=True)

    seq = pd.read_csv(screen_map_path)
    libraries = sorted(seq["Library"].dropna().astype(str).unique()) if "Library" in seq.columns else []
    print(f"ScreenSequenceMap libraries: {libraries}", flush=True)
    cd = seq[seq["ScreenID"].astype(str).str.contains(r"\.CD", regex=True)] if "ScreenID" in seq.columns else seq.iloc[0:0]
    cd_libraries = sorted(cd["Library"].dropna().astype(str).unique()) if "Library" in cd.columns else []
    print(f"CD-screen libraries: {cd_libraries} n_rows={len(cd)}", flush=True)

    wanted = set(SCREEN_GENES + EXTRA_RNA)

    # --- integrated Chronos ---
    crispr_slim_path = outdir / "crispr_gene_effect_slim.csv"
    crispr_gw_path = outdir / "crispr_cldn4_vs_all_rho.tsv"
    if not crispr_slim_path.exists() or not crispr_gw_path.exists():
        ensure(FILES["CRISPRGeneEffect.csv"], crispr_path, 300_000_000)
        print("Reading CRISPRGeneEffect (integrated Chronos)...", flush=True)
        crispr_full = read_gene_matrix(crispr_path)
        print(f"  shape {crispr_full.shape}", flush=True)
        present = [g for g in SCREEN_GENES if g in crispr_full.columns]
        missing = [g for g in SCREEN_GENES if g not in crispr_full.columns]
        print(f"  integrated present={present} missing={missing}", flush=True)
        gw = spearman_vs_query(crispr_full, QUERY, MIN_GW_N)
        gw.to_csv(crispr_gw_path, sep="\t", index=False)
        print(f"  genome-wide Chronos rhos: {len(gw)} genes", flush=True)
        slim_cols = ["row_id"] + present
        crispr_full[slim_cols].rename(columns={"row_id": "ModelID"}).to_csv(
            crispr_slim_path, index=False
        )
        del crispr_full
        try:
            crispr_path.unlink()
            print("Removed full CRISPRGeneEffect.csv", flush=True)
        except OSError:
            pass
    else:
        print("OK exists integrated Chronos slim + genome-wide rho", flush=True)

    # --- expression, lung genome-wide + slim (all models, selected genes) ---
    expr_slim_path = outdir / "expression_slim.csv"
    expr_gw_path = outdir / "expr_lung_cldn4_vs_all_rho.tsv"
    if not expr_slim_path.exists() or not expr_gw_path.exists():
        ensure(FILES["OmicsExpressionProteinCodingGenesTPMLogp1.csv"], expr_path, 400_000_000)
        print("Reading expression matrix...", flush=True)
        expr_full = read_gene_matrix(expr_path)
        print(f"  shape {expr_full.shape}", flush=True)
        present = [g for g in wanted if g in expr_full.columns]
        missing = sorted(wanted - set(present))
        print(f"  expression missing={missing}", flush=True)
        lung_expr = expr_full[expr_full["row_id"].isin(lung_ids)].copy()
        print(f"  lung cell lines in expression: {len(lung_expr)}", flush=True)
        gw = spearman_vs_query(lung_expr, QUERY, min_n=80)
        gw.to_csv(expr_gw_path, sep="\t", index=False)
        print(f"  lung genome-wide expression rhos: {len(gw)} genes", flush=True)
        slim_cols = ["row_id"] + present
        expr_full[slim_cols].rename(columns={"row_id": "ModelID"}).to_csv(
            expr_slim_path, index=False
        )
        del expr_full, lung_expr
        try:
            expr_path.unlink()
            print("Removed full expression matrix", flush=True)
        except OSError:
            pass
    else:
        print("OK exists expression slim + lung genome-wide rho", flush=True)

    # --- screen-level Chronos (Humagne-CD holds PRKDC) ---
    screen_slim_path = outdir / "screen_gene_effect_slim.csv"
    if not screen_slim_path.exists():
        ensure(FILES["ScreenGeneEffect.csv"], screen_path, 300_000_000)
        print("Reading ScreenGeneEffect (slim)...", flush=True)
        screen = read_gene_matrix(screen_path, wanted=set(SCREEN_GENES))
        screen = screen.rename(columns={"row_id": "ScreenID"})
        screen.to_csv(screen_slim_path, index=False)
        print(f"  screens {len(screen)}", flush=True)
        try:
            screen_path.unlink()
            print("Removed full ScreenGeneEffect.csv", flush=True)
        except OSError:
            pass
    else:
        print(f"OK exists {screen_slim_path}", flush=True)

    manifest = {
        "release": RELEASE,
        "doi": DOI,
        "files": FILES,
        "nhej_genes": NHEJ_GENES,
        "sting_genes": STING_GENES,
        "note_prkdc": (
            "PRKDC is not a column in CRISPRGeneEffect.csv (integrated Chronos). "
            "It is populated only on Humagne-CD Cas12 screens in ScreenGeneEffect.csv."
        ),
        "cd_libraries": cd_libraries,
        "libraries": libraries,
        "n_lung_cell_lines_model_csv": len(lung_ids),
        "sha256": {
            "Model.csv": sha256_file(model_path),
            "crispr_gene_effect_slim.csv": sha256_file(crispr_slim_path),
            "crispr_cldn4_vs_all_rho.tsv": sha256_file(crispr_gw_path),
            "expression_slim.csv": sha256_file(expr_slim_path),
            "expr_lung_cldn4_vs_all_rho.tsv": sha256_file(expr_gw_path),
            "screen_gene_effect_slim.csv": sha256_file(screen_slim_path),
        },
    }
    (outdir / "download_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({k: manifest[k] for k in ("release", "cd_libraries", "n_lung_cell_lines_model_csv")}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
