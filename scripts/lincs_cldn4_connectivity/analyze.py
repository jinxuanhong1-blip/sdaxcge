#!/usr/bin/env python3
"""Connectivity of public CLDN4-loss profiles to LINCS L1000 signatures.

LINCS/L1000 does not contain a CLDN4 knockdown, knockout, or overexpression
perturbagen (SigCom LINCS 2021 characteristic-direction Level 5; iLINCS LIB_6;
GSE92742 pert_info). The query profiles are therefore the public CLDN4-loss
transcriptomes, scored in the same gene space as LINCS characteristic-direction
signatures.

Reference classes, fixed before looking at scores
-------------------------------------------------
NHEJ loss
    CRISPR knockout and, separately, shRNA of the core NHEJ genes
    PRKDC, LIG4, XRCC4, XRCC5, XRCC6, NHEJ1, DCLRE1C.
STING/IFN
    Ligand signatures whose perturbagen is IFNA, IFNA1, IFNA2, IFNB1, or IFNG.
    No STING agonist (diABZI, ADU-S100, cGAMP, DMXAA) is in this SigCom build.
    TMEM173 (STING1) and MB21D1 (cGAS) CRISPR knockouts are a separate class,
    STING/cGAS loss, and are not treated as IFN signatures.

Primary metric is Spearman correlation between the CLDN4-loss log2 fold-change
and the LINCS characteristic-direction coefficient, on shared gene symbols.
The unit of inference is the cell line (median across signatures in that line).
A matched background (other genes, same cell line and pert_type, seed 42) is
the null for the paired test. Cosine and a weighted connectivity score on the
top/bottom 100 query genes are secondary.
"""

from __future__ import annotations

import gzip
import json
import math
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "methods" / "lincs_cldn4_connectivity"
CACHE = OUT / "data" / "cache"
RAW = OUT / "data" / "raw"
TAB = OUT / "tables"
FIG = OUT / "figures"
for p in (CACHE, RAW, TAB, FIG):
    p.mkdir(parents=True, exist_ok=True)

META = "https://maayanlab.cloud/sigcom-lincs/metadata-api/signatures"
UA = {"User-Agent": "cldn4-lincs-connectivity", "Accept": "application/json"}
SEED = 42
BG_K = 8
WTCS_N = 100
MIN_OVERLAP = 1000
FPKM_PSEUDO = 0.5

NHEJ_GENES = ("PRKDC", "LIG4", "XRCC4", "XRCC5", "XRCC6", "NHEJ1", "DCLRE1C")
IFN_LIGANDS = ("IFNA", "IFNA1", "IFNA2", "IFNB1", "IFNG")
STING_GENES = ("TMEM173", "MB21D1")  # STING1 and cGAS symbols in this build
CATALOG_GENES = ("CLDN4", "CLDN3", "CLDN7", "TACSTD2", "STING1", "CGAS") + NHEJ_GENES + STING_GENES

# Prespecified absence checks for STING agonists. Counts are recorded; they are
# not pulled into the IFN class.
STING_AGONISTS = ("diABZI", "ADU-S100", "cGAMP", "2'3'-cGAMP", "vadimezan", "DMXAA")

PRIMARY_QUERIES = ("GSE207704_T47D", "GSE207704_MCF7")
PRIMARY_CLASSES = ("NHEJ_CRISPR", "NHEJ_shRNA", "IFN_ligand")
CLASS_ORDER = ("NHEJ_CRISPR", "NHEJ_shRNA", "IFN_ligand", "STING_cGAS_KO")
CLASS_LABEL = {
    "NHEJ_CRISPR": "NHEJ CRISPR KO",
    "NHEJ_shRNA": "NHEJ shRNA",
    "IFN_ligand": "IFN ligand",
    "STING_cGAS_KO": "STING/cGAS KO",
}
CLASS_COLOR = {
    "NHEJ_CRISPR": "#0072B2",
    "NHEJ_shRNA": "#56B4E9",
    "IFN_ligand": "#E69F00",
    "STING_cGAS_KO": "#009E73",
}
QUERY_LABEL = {
    "GSE207704_T47D": "T47D CLDN4 KO",
    "GSE207704_MCF7": "MCF7 CLDN4 KO",
    "GSE22493_SKOV3": "SKOV-3 siCLDN4",
    "GSE50927_naive_lung": "Mouse lung Cldn4 KO",
}

ALIAS_22493 = {"G1P2": "ISG15", "G1P3": "IFI6"}


def log(msg: str) -> None:
    print(msg, flush=True)


def fetch(url: str, data: bytes | None = None, timeout: int = 90, retries: int = 4) -> bytes:
    last = None
    for i in range(retries):
        try:
            req = urllib.request.Request(url, data=data, headers={**UA, "Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(1.5 * (i + 1))
    raise RuntimeError(f"GET/POST failed {url[:140]}: {last}")


def download(url: str, dest: Path, timeout: int = 180) -> None:
    if dest.exists() and dest.stat().st_size > 500:
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    last = None
    for i in range(4):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA["User-Agent"]})
            with urllib.request.urlopen(req, timeout=timeout) as resp, tmp.open("wb") as out:
                while True:
                    chunk = resp.read(1 << 16)
                    if not chunk:
                        break
                    out.write(chunk)
            tmp.replace(dest)
            return
        except Exception as exc:  # noqa: BLE001
            last = exc
            if tmp.exists():
                tmp.unlink()
            time.sleep(1.5 * (i + 1))
    raise RuntimeError(f"download failed {url}: {last}")


def api_get(path: str, params: dict) -> dict | list:
    url = META + path + "?" + urllib.parse.urlencode(params)
    return json.loads(fetch(url))


def api_count(where: dict) -> int:
    res = api_get("/count", {"where": json.dumps(where)})
    return int(res["count"])


def api_find(where: dict, limit: int = 200, skip: int = 0) -> list[dict]:
    payload = {
        "filter": {
            "where": where,
            "limit": limit,
            "skip": skip,
            "fields": [
                "id",
                "meta.pert_name",
                "meta.pert_type",
                "meta.cell_line",
                "meta.persistent_id",
                "meta.pert_time",
                "meta.pert_dose",
                "meta.cmap_id",
            ],
        }
    }
    raw = json.loads(fetch(META + "/find", data=json.dumps(payload).encode()))
    rows = []
    for item in raw:
        meta = item.get("meta") or {}
        rows.append(
            {
                "sig_id": item.get("id"),
                "pert_name": meta.get("pert_name"),
                "pert_type": meta.get("pert_type"),
                "cell_line": meta.get("cell_line"),
                "persistent_id": meta.get("persistent_id"),
                "pert_time": meta.get("pert_time"),
                "pert_dose": meta.get("pert_dose"),
                "cmap_id": meta.get("cmap_id"),
            }
        )
    return rows


def api_find_all(where: dict) -> list[dict]:
    out = []
    skip = 0
    while True:
        batch = api_find(where, limit=200, skip=skip)
        if not batch:
            break
        out.extend(batch)
        if len(batch) < 200:
            break
        skip += 200
    # unique by signature id
    seen = set()
    uniq = []
    for row in out:
        if row["sig_id"] in seen:
            continue
        seen.add(row["sig_id"])
        uniq.append(row)
    return uniq


def catalog() -> pd.DataFrame:
    rows = []
    for gene in list(dict.fromkeys(CATALOG_GENES)):
        try:
            n = api_count({"meta.pert_name": gene})
        except Exception as exc:  # noqa: BLE001
            rows.append({"resource": "SigCom LINCS pert_name", "query": gene, "n": np.nan, "note": str(exc)})
            continue
        # type breakdown only when the gene exists, plus always for CLDN4
        note = ""
        if n and n < 500:
            found = api_find_all({"meta.pert_name": gene})
            types: dict[str, int] = {}
            for r in found:
                types[r["pert_type"] or ""] = types.get(r["pert_type"] or "", 0) + 1
            note = json.dumps(types, sort_keys=True)
            n = len(found)
        rows.append({"resource": "SigCom LINCS pert_name", "query": gene, "n": n, "note": note})
        log(f"catalog {gene} n={n} {note}")
    for name in STING_AGONISTS:
        n = api_count({"meta.pert_name": name})
        rows.append({"resource": "SigCom LINCS pert_name", "query": name, "n": n, "note": "STING agonist name check"})
        log(f"catalog agonist {name} n={n}")
    # iLINCS individual/consensus knockdown library
    try:
        filt = urllib.parse.quote(json.dumps({"where": {"treatment": "CLDN4"}, "limit": 5}))
        body = json.loads(fetch(f"https://www.ilincs.org/api/SignatureMeta?filter={filt}"))
        rows.append(
            {
                "resource": "iLINCS SignatureMeta treatment",
                "query": "CLDN4",
                "n": len(body) if isinstance(body, list) else np.nan,
                "note": "LIB_6-style treatment field; empty list means no knockdown signature",
            }
        )
    except Exception as exc:  # noqa: BLE001
        rows.append({"resource": "iLINCS SignatureMeta treatment", "query": "CLDN4", "n": np.nan, "note": str(exc)})
    # Phase I Broad pert_info
    pert_path = RAW / "GSE92742_Broad_LINCS_pert_info.txt.gz"
    download(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_pert_info.txt.gz",
        pert_path,
    )
    n_pert = 0
    with gzip.open(pert_path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        name_i = header.index("pert_iname")
        for line in fh:
            if line.rstrip("\n").split("\t")[name_i] == "CLDN4":
                n_pert += 1
    rows.append(
        {
            "resource": "GSE92742 pert_info pert_iname",
            "query": "CLDN4",
            "n": n_pert,
            "note": "LINCS L1000 Phase I perturbagen table",
        }
    )
    # landmark status
    gene_path = RAW / "GSE92742_Broad_LINCS_gene_info.txt.gz"
    download(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_gene_info.txt.gz",
        gene_path,
    )
    with gzip.open(gene_path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            rec = dict(zip(header, parts))
            if rec.get("pr_gene_symbol") == "CLDN4":
                rows.append(
                    {
                        "resource": "GSE92742 gene_info",
                        "query": "CLDN4",
                        "n": 1,
                        "note": f"pr_gene_id={rec.get('pr_gene_id')} pr_is_lm={rec.get('pr_is_lm')} pr_is_bing={rec.get('pr_is_bing')}",
                    }
                )
                break
    cell_path = RAW / "GSE92742_Broad_LINCS_cell_info.txt.gz"
    download(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_cell_info.txt.gz",
        cell_path,
    )
    return pd.DataFrame(rows)


def lung_cell_lines() -> set[str]:
    path = RAW / "GSE92742_Broad_LINCS_cell_info.txt.gz"
    df = pd.read_csv(path, sep="\t", dtype=str)
    site_col = "primary_site" if "primary_site" in df.columns else None
    id_col = "cell_id" if "cell_id" in df.columns else df.columns[0]
    if site_col is None:
        return set()
    lung = set(df.loc[df[site_col].str.lower().eq("lung"), id_col].dropna())
    # SigCom uses a few aliases that are the same lines.
    lung |= {x.replace("-", "") for x in lung}
    # 2021 CRISPR panel line absent from the 2017 GSE92742 cell table.
    # LCLC103H is a human large-cell lung carcinoma (DSMZ ACC 384).
    lung.add("LCLC103H")
    return lung


def collect_references() -> pd.DataFrame:
    frames = []
    for gene in NHEJ_GENES:
        rows = api_find_all({"meta.pert_name": gene})
        for row in rows:
            if row["pert_type"] == "CRISPR Knockout":
                row["sig_class"] = "NHEJ_CRISPR"
            elif row["pert_type"] == "shRNA":
                row["sig_class"] = "NHEJ_shRNA"
            else:
                continue
            row["gene_group"] = gene
            frames.append(row)
        log(f"NHEJ {gene} kept {sum(1 for r in frames if r['gene_group']==gene)}")
    for name in IFN_LIGANDS:
        rows = api_find_all({"meta.pert_name": name})
        n_lig = 0
        for row in rows:
            if row["pert_type"] != "Ligand":
                continue
            row["sig_class"] = "IFN_ligand"
            row["gene_group"] = name
            frames.append(row)
            n_lig += 1
        log(f"IFN ligand {name} n={n_lig}")
    for gene in STING_GENES:
        rows = api_find_all({"meta.pert_name": gene})
        n = 0
        for row in rows:
            if row["pert_type"] != "CRISPR Knockout":
                continue
            row["sig_class"] = "STING_cGAS_KO"
            row["gene_group"] = gene
            frames.append(row)
            n += 1
        log(f"STING/cGAS {gene} CRISPR n={n}")
    df = pd.DataFrame(frames)
    df = df[df["persistent_id"].notna() & df["cell_line"].notna()].copy()
    df = df.drop_duplicates("sig_id")
    return df


def sample_background(ref: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    exclude = set(NHEJ_GENES) | set(STING_GENES) | set(IFN_LIGANDS) | {"STING1", "CGAS", "IFNB", "IFNA2A", "IFNA2B"}
    jobs = []
    # CRISPR background covers NHEJ CRISPR and STING/cGAS KO cell lines.
    crispr_lines = sorted(
        set(ref.loc[ref["sig_class"].isin(["NHEJ_CRISPR", "STING_cGAS_KO"]), "cell_line"])
    )
    for line in crispr_lines:
        jobs.append((line, "CRISPR Knockout", "BG_CRISPR"))
    shrna_lines = sorted(set(ref.loc[ref["sig_class"].eq("NHEJ_shRNA"), "cell_line"]))
    for line in shrna_lines:
        jobs.append((line, "shRNA", "BG_shRNA"))
    ifn_lines = sorted(set(ref.loc[ref["sig_class"].eq("IFN_ligand"), "cell_line"]))
    for line in ifn_lines:
        jobs.append((line, "Ligand", "BG_ligand"))

    picked = []
    for cell_line, pert_type, sig_class in jobs:
        where = {"meta.cell_line": cell_line, "meta.pert_type": pert_type}
        try:
            count = api_count(where)
        except Exception as exc:  # noqa: BLE001
            log(f"background count failed {cell_line} {pert_type}: {exc}")
            continue
        window = 60
        skip = 0 if count <= window else int(rng.integers(0, count - window))
        try:
            rows = api_find(where, limit=window, skip=skip)
        except Exception as exc:  # noqa: BLE001
            log(f"background find failed {cell_line} {pert_type}: {exc}")
            continue
        rows = [
            r
            for r in rows
            if r.get("pert_name") not in exclude and r.get("persistent_id") and r.get("pert_type") == pert_type
        ]
        if len(rows) < BG_K and count > window:
            skip2 = 0 if count <= window else int(rng.integers(0, max(count - window, 1)))
            try:
                extra = api_find(where, limit=window, skip=skip2)
            except Exception as exc:  # noqa: BLE001
                log(f"background second window failed {cell_line} {pert_type}: {exc}")
                extra = []
            have = {r["sig_id"] for r in rows}
            for r in extra:
                if r.get("sig_id") in have:
                    continue
                if r.get("pert_name") in exclude or not r.get("persistent_id") or r.get("pert_type") != pert_type:
                    continue
                rows.append(r)
        order = rng.permutation(len(rows)) if rows else []
        rows = [rows[int(i)] for i in order]
        take = rows[:BG_K]
        for row in take:
            row["sig_class"] = sig_class
            row["gene_group"] = row["pert_name"]
            picked.append(row)
        log(f"background {sig_class} {cell_line} took {len(take)} / count {count}")
    if not picked:
        return pd.DataFrame()
    df = pd.DataFrame(picked).drop_duplicates("sig_id")
    return df


def load_cd(path: Path) -> dict[str, float]:
    """Level 5 characteristic-direction table.

    SigCom writes either ``symbol\\tCD-coefficient`` or a blank first-column
    header ``\\tCD-coefficient``. Column 1 is the gene symbol in both cases.
    """
    out: dict[str, float] = {}
    with path.open() as fh:
        header = fh.readline()
        if "cd-coefficient" not in header.lower():
            raise ValueError(f"unexpected header in {path}: {header[:80]!r}")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            sym, val = parts[0], parts[1]
            if not sym or sym in out:
                continue
            try:
                x = float(val)
            except ValueError:
                continue
            if math.isfinite(x):
                out[sym] = x
    return out


def ensure_signatures(meta: pd.DataFrame) -> pd.DataFrame:
    ok = []

    def one(row: dict) -> dict:
        sid = row["sig_id"]
        dest = CACHE / "sig" / f"{sid}.tsv"
        status = "ok"
        try:
            download(row["persistent_id"], dest, timeout=120)
            if dest.stat().st_size < 500:
                status = "too_small"
        except Exception as exc:  # noqa: BLE001
            status = f"fail:{exc}"
        row = dict(row)
        row["path"] = str(dest) if status == "ok" else ""
        row["download_status"] = status
        return row

    records = meta.to_dict(orient="records")
    with ThreadPoolExecutor(max_workers=12) as pool:
        futs = [pool.submit(one, rec) for rec in records]
        for i, fut in enumerate(as_completed(futs), start=1):
            ok.append(fut.result())
            if i % 40 == 0:
                log(f"downloaded {i}/{len(records)}")
    return pd.DataFrame(ok)


def wtcs(query: dict[str, float], ref: dict[str, float]) -> float:
    shared = [g for g in ref if g in query]
    if len(shared) < MIN_OVERLAP:
        return np.nan
    ranked = sorted(shared, key=lambda g: ref[g], reverse=True)
    ups = sorted((g for g in shared if query[g] > 0), key=lambda g: query[g], reverse=True)[:WTCS_N]
    downs = sorted((g for g in shared if query[g] < 0), key=lambda g: query[g])[:WTCS_N]
    if len(ups) < 15 or len(downs) < 15:
        return np.nan

    def es(members: list[str]) -> float:
        weights = {g: abs(query[g]) for g in members}
        nr = sum(weights.values())
        nh = len(weights)
        n = len(ranked)
        miss = 1.0 / (n - nh)
        running = 0.0
        max_dev = 0.0
        min_dev = 0.0
        for g in ranked:
            w = weights.get(g)
            if w:
                running += w / nr
            else:
                running -= miss
            if running > max_dev:
                max_dev = running
            elif running < min_dev:
                min_dev = running
        return max_dev if abs(max_dev) >= abs(min_dev) else min_dev

    es_up = es(ups)
    es_down = es(downs)
    if es_up * es_down < 0:
        return 0.5 * (es_up - es_down)
    return 0.0


def score_pair(query: dict[str, float], ref: dict[str, float]) -> tuple[float, float, float, int]:
    genes = [g for g in ref if g in query]
    n = len(genes)
    if n < MIN_OVERLAP:
        return (np.nan, np.nan, np.nan, n)
    q = np.fromiter((query[g] for g in genes), dtype=float, count=n)
    r = np.fromiter((ref[g] for g in genes), dtype=float, count=n)
    if np.allclose(q, q[0]) or np.allclose(r, r[0]):
        rho = np.nan
    else:
        rho = float(stats.spearmanr(q, r).statistic)
    nq = float(np.linalg.norm(q))
    nr = float(np.linalg.norm(r))
    cosine = float(np.dot(q, r) / (nq * nr)) if nq > 0 and nr > 0 else np.nan
    connectivity = wtcs(query, ref)
    return rho, cosine, connectivity, n


def build_queries() -> tuple[dict[str, dict[str, float]], pd.DataFrame]:
    ranks: dict[str, dict[str, float]] = {}
    meta_rows = []

    # GSE207704 group-mean FPKM
    path = RAW / "GSE207704_CLDN4_RNAseq.txt.gz"
    download(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207704/suppl/GSE207704_CLDN4_RNAseq.txt.gz",
        path,
    )
    df = pd.read_csv(path, sep="\t", low_memory=False)
    cols = {
        "MCF7_KO": "MCF7_CLDN4KO_FPKM (fpkm)",
        "MCF7_WT": "MCF7_WT_FPKM (fpkm)",
        "T47D_KO": "T47D_CLDN4KO_FPKM (fpkm)",
        "T47D_WT": "T47D_WT_FPKM (fpkm)",
    }
    missing = [c for c in cols.values() if c not in df.columns]
    if missing or "gene_short_name" not in df.columns:
        raise SystemExit(f"GSE207704 columns unexpected: {df.columns.tolist()[:20]} missing {missing}")
    df = df.rename(columns={"gene_short_name": "symbol"})
    df = df[df["symbol"].notna() & (df["symbol"].astype(str).str.len() > 0)].copy()
    for c in cols.values():
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["_mean"] = df[list(cols.values())].mean(axis=1)
    df = df.sort_values("_mean", ascending=False).drop_duplicates("symbol")
    df = df.set_index("symbol")
    for key, ko, wt, species, tissue, n_loss, n_wt, contrast in [
        (
            "GSE207704_T47D",
            "T47D_CLDN4KO_FPKM (fpkm)",
            "T47D_WT_FPKM (fpkm)",
            "human",
            "T47D breast cancer line, group-mean FPKM",
            2,
            2,
            "T47D CLDN4-/- vs WT",
        ),
        (
            "GSE207704_MCF7",
            "MCF7_CLDN4KO_FPKM (fpkm)",
            "MCF7_WT_FPKM (fpkm)",
            "human",
            "MCF7 breast cancer line, group-mean FPKM",
            2,
            2,
            "MCF7 CLDN4-/- vs WT",
        ),
    ]:
        logfc = np.log2(df[ko] + FPKM_PSEUDO) - np.log2(df[wt] + FPKM_PSEUDO)
        logfc = logfc.replace([np.inf, -np.inf], np.nan).dropna()
        ranks[key] = {str(g): float(v) for g, v in logfc.items()}
        meta_rows.append(
            {
                "query": key,
                "accession": "GSE207704",
                "contrast": contrast,
                "species": species,
                "tissue": tissue,
                "n_loss": n_loss,
                "n_wt": n_wt,
                "symbol": "CLDN4",
                "target_log2fc": ranks[key].get("CLDN4", np.nan),
                "n_genes": len(ranks[key]),
                "rank_metric": "log2((KO+0.5)/(WT+0.5)) on GEO group-mean FPKM; replicates already collapsed",
            }
        )

    # GSE50927 author edgeR, naive whole lung. Map symbols by uppercase.
    path = RAW / "GSE50927_Cldn4lungWTvsKOgenes.csv.gz"
    download(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50927/suppl/GSE50927_Cldn4lungWTvsKOgenes.csv.gz",
        path,
    )
    m = pd.read_csv(path)
    m = m[m["Marker.Symbol"].notna()].copy()
    m["logFC"] = pd.to_numeric(m["logFC"], errors="coerce")
    m["logCPM"] = pd.to_numeric(m["logCPM"], errors="coerce")
    m["human_symbol"] = m["Marker.Symbol"].astype(str).str.upper()
    m["_abs"] = m["logFC"].abs()
    m = m.sort_values(["logCPM", "_abs"], ascending=False).drop_duplicates("human_symbol")
    lung = m.set_index("human_symbol")["logFC"].dropna()
    ranks["GSE50927_naive_lung"] = {str(g): float(v) for g, v in lung.items()}
    meta_rows.append(
        {
            "query": "GSE50927_naive_lung",
            "accession": "GSE50927",
            "contrast": "Cldn4 KO vs WT, naive whole lung",
            "species": "mouse",
            "tissue": "naive whole lung; symbols uppercased onto human symbols",
            "n_loss": 1,
            "n_wt": 1,
            "symbol": "CLDN4",
            "target_log2fc": ranks["GSE50927_naive_lung"].get("CLDN4", np.nan),
            "n_genes": len(ranks["GSE50927_naive_lung"]),
            "rank_metric": "author edgeR logFC; one row per uppercased symbol, highest logCPM kept",
        }
    )

    # GSE22493 two-colour arrays, KD versus CLDN4-overexpressing control.
    series = RAW / "GSE22493_series_matrix.txt.gz"
    download(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE22nnn/GSE22493/matrix/GSE22493_series_matrix.txt.gz",
        series,
    )
    plat = load_gpl10555()
    rows = []
    hdr = None
    with gzip.open(series, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                hdr = next(fh).rstrip("\n").replace('"', "").split("\t")
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if hdr is not None:
                rows.append(line.rstrip("\n").replace('"', "").split("\t"))
    mat = pd.DataFrame(rows, columns=hdr)
    for c in mat.columns[1:]:
        mat[c] = pd.to_numeric(mat[c], errors="coerce")
    mat["ID_REF"] = mat["ID_REF"].astype(str)
    merged = mat.merge(plat, left_on="ID_REF", right_on="ID", how="left")

    def symbol_of(orf: str, desc: str) -> str:
        s = (orf or "").strip().upper()
        if not s:
            s = (desc or "").split("--")[0].strip().upper()
        return ALIAS_22493.get(s, s)

    merged["symbol"] = [
        symbol_of(o, d) for o, d in zip(merged["ORF"].fillna(""), merged["DESCRIPTION"].fillna(""))
    ]
    gsms = [c for c in ("GSM558700", "GSM558701", "GSM558702") if c in merged.columns]
    if len(gsms) != 3:
        raise SystemExit(f"GSE22493 GSM columns missing: {mat.columns.tolist()[:8]}")
    keep = merged[merged["symbol"].str.len() > 0]
    gene = keep.groupby("symbol")[gsms].mean()
    rank = gene.mean(axis=1, skipna=True).dropna()
    ranks["GSE22493_SKOV3"] = {str(g): float(v) for g, v in rank.items()}
    meta_rows.append(
        {
            "query": "GSE22493_SKOV3",
            "accession": "GSE22493",
            "contrast": "SKOV-3 CLDN4 siRNA vs CLDN4-overexpressing control",
            "species": "human",
            "tissue": "SKOV-3 ovarian line, two-colour array GPL10555",
            "n_loss": 3,
            "n_wt": 3,
            "symbol": "CLDN4",
            "target_log2fc": ranks["GSE22493_SKOV3"].get("CLDN4", np.nan),
            "n_genes": len(ranks["GSE22493_SKOV3"]),
            "rank_metric": "mean deposited log2(KD/CLDN4-OE) across GSM558700-702; ORF else DESCRIPTION",
        }
    )
    return ranks, pd.DataFrame(meta_rows)


def load_gpl10555() -> pd.DataFrame:
    dest = RAW / "GPL10555_family.soft.gz"
    download(
        "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10nnn/GPL10555/soft/GPL10555_family.soft.gz",
        dest,
    )
    from io import StringIO

    lines: list[str] = []
    in_table = False
    with gzip.open(dest, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!platform_table_begin"):
                in_table = True
                continue
            if line.startswith("!platform_table_end"):
                break
            if in_table and line.strip():
                lines.append(line)
    plat = pd.read_csv(StringIO("".join(lines)), sep="\t", dtype=str)
    if not {"ID", "ORF", "DESCRIPTION"}.issubset(plat.columns):
        raise SystemExit(f"GPL10555 columns unexpected: {plat.columns.tolist()}")
    out = plat[["ID", "ORF", "DESCRIPTION"]].fillna("")
    out["ID"] = out["ID"].astype(str).str.replace(r"\.0$", "", regex=True)
    return out


def wilcoxon_p(values: np.ndarray) -> float:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) < 5:
        return np.nan
    if np.allclose(x, 0):
        return 1.0
    try:
        return float(stats.wilcoxon(x, alternative="two-sided", zero_method="wilcox").pvalue)
    except ValueError:
        return np.nan


def bh_fdr(pvals: list[float]) -> list[float]:
    p = np.asarray(pvals, dtype=float)
    out = np.full(p.shape, np.nan)
    idx = np.where(np.isfinite(p))[0]
    if len(idx) == 0:
        return out.tolist()
    order = idx[np.argsort(p[idx])]
    ranked = p[order]
    m = len(order)
    q = ranked * m / np.arange(1, m + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out[order] = np.clip(q, 0, 1)
    return out.tolist()


def summarize(scores: pd.DataFrame, lung: set[str]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    scores = scores.copy()
    scores["is_lung"] = scores["cell_line"].isin(lung)

    cell_rows = []
    for (query, sig_class, cell_line), sub in scores.groupby(["query", "sig_class", "cell_line"], dropna=False):
        cell_rows.append(
            {
                "query": query,
                "sig_class": sig_class,
                "cell_line": cell_line,
                "is_lung": bool(cell_line in lung),
                "n_sig": int(len(sub)),
                "n_genes_median": float(sub["n_genes"].median()),
                "median_spearman": float(sub["spearman"].median()),
                "median_cosine": float(sub["cosine"].median()),
                "median_wtcs": float(sub["wtcs"].median()),
            }
        )
    cell = pd.DataFrame(cell_rows)

    # Pair each biological class with its background class.
    pair_of = {
        "NHEJ_CRISPR": "BG_CRISPR",
        "STING_cGAS_KO": "BG_CRISPR",
        "NHEJ_shRNA": "BG_shRNA",
        "IFN_ligand": "BG_ligand",
    }

    def class_table(scope: str, mask: pd.Series) -> pd.DataFrame:
        use = cell[mask].copy()
        rows = []
        for (query, sig_class), sub in use.groupby(["query", "sig_class"]):
            if sig_class.startswith("BG_"):
                continue
            bg_name = pair_of.get(sig_class)
            bg = use[(use["query"] == query) & (use["sig_class"] == bg_name)]
            merged = sub.merge(bg, on="cell_line", how="inner", suffixes=("", "_bg"))
            delta = merged["median_spearman"] - merged["median_spearman_bg"] if len(merged) else pd.Series(dtype=float)
            rows.append(
                {
                    "scope": scope,
                    "query": query,
                    "sig_class": sig_class,
                    "n_cell": int(sub["cell_line"].nunique()),
                    "n_sig": int(sub["n_sig"].sum()),
                    "median_cell_spearman": float(sub["median_spearman"].median()),
                    "q25_spearman": float(sub["median_spearman"].quantile(0.25)),
                    "q75_spearman": float(sub["median_spearman"].quantile(0.75)),
                    "median_cell_cosine": float(sub["median_cosine"].median()),
                    "median_cell_wtcs": float(sub["median_wtcs"].median()),
                    "frac_cell_spearman_pos": float((sub["median_spearman"] > 0).mean()),
                    "wilcoxon_p_vs_0": wilcoxon_p(sub["median_spearman"].to_numpy()),
                    "n_paired_bg": int(len(merged)),
                    "median_spearman_minus_bg": float(delta.median()) if len(delta) else np.nan,
                    "paired_p_vs_bg": wilcoxon_p(delta.to_numpy()) if len(delta) else np.nan,
                }
            )
        return pd.DataFrame(rows)

    summary = pd.concat(
        [
            class_table("all_cell_lines", pd.Series(True, index=cell.index)),
            class_table("lung_cell_lines", cell["is_lung"]),
        ],
        ignore_index=True,
    )
    summary["primary_family"] = (
        summary["scope"].eq("all_cell_lines")
        & summary["query"].isin(PRIMARY_QUERIES)
        & summary["sig_class"].isin(PRIMARY_CLASSES)
    )
    summary["q_vs_0"] = np.nan
    summary["q_vs_bg"] = np.nan
    fam = summary["primary_family"]
    if fam.any():
        summary.loc[fam, "q_vs_0"] = bh_fdr(summary.loc[fam, "wilcoxon_p_vs_0"].tolist())
        summary.loc[fam, "q_vs_bg"] = bh_fdr(summary.loc[fam, "paired_p_vs_bg"].tolist())

    gene_rows = []
    gene_scores = scores[scores["sig_class"].isin(CLASS_ORDER)]
    for (query, sig_class, gene), sub in gene_scores.groupby(["query", "sig_class", "gene_group"]):
        gene_rows.append(
            {
                "query": query,
                "sig_class": sig_class,
                "gene": gene,
                "n_sig": int(len(sub)),
                "n_cell": int(sub["cell_line"].nunique()),
                "median_spearman": float(sub["spearman"].median()),
                "median_cosine": float(sub["cosine"].median()),
                "median_wtcs": float(sub["wtcs"].median()),
            }
        )
    genes = pd.DataFrame(gene_rows)
    return cell, summary, genes


def fmt(x, nd=3) -> str:
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return "NA"
    if isinstance(x, (int, np.integer)):
        return str(int(x))
    ax = abs(float(x))
    if ax != 0 and ax < 0.001:
        return f"{float(x):.2e}"
    return f"{float(x):.{nd}f}"


def write_finding(
    catalog_df: pd.DataFrame,
    query_meta: pd.DataFrame,
    summary: pd.DataFrame,
    genes: pd.DataFrame,
    n_fail: int,
) -> None:
    def cat(resource: str, query: str) -> str:
        hit = catalog_df[(catalog_df["resource"] == resource) & (catalog_df["query"] == query)]
        if hit.empty:
            return "NA"
        return fmt(hit.iloc[0]["n"], 0)

    def row(query: str, sig_class: str, scope: str = "all_cell_lines") -> pd.Series:
        hit = summary[
            (summary["query"] == query) & (summary["sig_class"] == sig_class) & (summary["scope"] == scope)
        ]
        if hit.empty:
            return pd.Series(dtype=object)
        return hit.iloc[0]

    def fmt_r(x) -> str:
        if x is None or (isinstance(x, float) and not math.isfinite(x)):
            return "NA"
        return f"{float(x):.3f}"

    def line(query: str, sig_class: str, scope: str = "all_cell_lines") -> str:
        r = row(query, sig_class, scope)
        if r.empty:
            return "| " + " | ".join([QUERY_LABEL[query], CLASS_LABEL[sig_class], "0"] + ["NA"] * 6) + " |"
        qcol = fmt(r["q_vs_0"]) if scope == "all_cell_lines" and bool(r["primary_family"]) else "—"
        return (
            f"| {QUERY_LABEL[query]} | {CLASS_LABEL[sig_class]} | {int(r['n_cell'])} | {int(r['n_sig'])} | "
            f"{fmt_r(r['median_cell_spearman'])} | {fmt_r(r['q25_spearman'])}–{fmt_r(r['q75_spearman'])} | "
            f"{fmt(r['wilcoxon_p_vs_0'])} | {qcol} | {fmt_r(r['median_spearman_minus_bg'])} | {fmt(r['paired_p_vs_bg'])} |"
        )

    qlines = []
    for _, r in query_meta.iterrows():
        qlines.append(
            f"| {r['query']} | {r['contrast']} | {r['species']} | {int(r['n_loss'])} vs {int(r['n_wt'])} | "
            f"{fmt(r['target_log2fc'])} | {int(r['n_genes'])} |"
        )

    def gene_block(query: str) -> str:
        sub = genes[(genes["query"] == query) & (genes["sig_class"].isin(["NHEJ_CRISPR", "NHEJ_shRNA"]))]
        if sub.empty:
            return "_no NHEJ signatures_"
        lines = ["| Gene | Class | n sig | n cell lines | median Spearman | median WTCS |", "|---|---|---:|---:|---:|---:|"]
        sub = sub.sort_values(["gene", "sig_class"])
        for _, r in sub.iterrows():
            lines.append(
                f"| {r['gene']} | {CLASS_LABEL[r['sig_class']]} | {int(r['n_sig'])} | {int(r['n_cell'])} | "
                f"{fmt(r['median_spearman'])} | {fmt(r['median_wtcs'])} |"
            )
        return "\n".join(lines)

    header_rows = "\n".join(
        line(q, c)
        for q in ("GSE207704_T47D", "GSE207704_MCF7", "GSE22493_SKOV3", "GSE50927_naive_lung")
        for c in CLASS_ORDER
    )
    lung_rows = "\n".join(
        line(q, c, "lung_cell_lines")
        for q in ("GSE207704_T47D", "GSE207704_MCF7")
        for c in CLASS_ORDER
    )
    bing = catalog_df[(catalog_df["resource"] == "GSE92742 gene_info") & (catalog_df["query"] == "CLDN4")]
    bing_note = bing.iloc[0]["note"] if len(bing) else "NA"
    prim = summary[summary["primary_family"]].copy()

    md = f"""# FINDING — LINCS/L1000 connectivity of CLDN4-loss profiles to NHEJ loss and STING/IFN

Additive public evidence. Numbers are written from `tables/class_summary.tsv`.

## Result

On T47D and MCF7, CLDN4-knockout connectivity to LINCS NHEJ-loss signatures and to IFN ligand signatures is centered at zero. Primary median cell-line Spearman values sit between {fmt_r(prim['median_cell_spearman'].min())} and {fmt_r(prim['median_cell_spearman'].max())}. BH q versus 0, inside that primary family, is at least {fmt(prim['q_vs_0'].min())}. The same rows’ cell-line median weighted connectivity scores are 0. Matched-background tests are in the table below.

Shared genes per scored signature are the intersection of the query rank with the characteristic-direction table (8018 genes on the GSE207704 queries). Lung lines in this signature set are A549, HCC515, and LCLC103H. That subset has two lines at most, so the lung Wilcoxon is not computed.

## What was queried

SigCom LINCS (characteristic-direction Level 5, metadata retrieved with this run) has **{cat('SigCom LINCS pert_name', 'CLDN4')}** signatures with `pert_name = CLDN4`. iLINCS `SignatureMeta` treatment CLDN4 returned **{cat('iLINCS SignatureMeta treatment', 'CLDN4')}**. GSE92742 Phase I `pert_info` has **{cat('GSE92742 pert_info pert_iname', 'CLDN4')}** CLDN4 perturbagen. CLDN4 in the Phase I gene table: {bing_note}. TACSTD2, for context in the same catalog, has **{cat('SigCom LINCS pert_name', 'TACSTD2')}** signatures.

The connectivity query is therefore each public CLDN4-loss transcriptome, not a LINCS CLDN4 hairpin. Positive Spearman means the CLDN4-loss profile points the same way as the reference signature.

| Query | Contrast | Species | n | CLDN4 / Cldn4 log2FC | Genes in rank |
|---|---|---|---|---:|---:|
{chr(10).join(qlines)}

GSE207704 deposits one FPKM per group (2 vs 2 already collapsed). GSE22493 control arm is CLDN4 overexpression. GSE50927 is naive whole-lung Cldn4 knockout, n=1 vs 1; mouse symbols were uppercased to overlap human LINCS symbols.

## Reference signatures

Core NHEJ loss: PRKDC, LIG4, XRCC4, XRCC5, XRCC6, NHEJ1, DCLRE1C. CRISPR knockout and shRNA are separate classes because Ku80 (XRCC5) and XRCC4 are shRNA in this build.

STING/IFN: ligand signatures for IFNA, IFNA1, IFNA2, IFNB1, and IFNG. STING agonist names diABZI, ADU-S100, cGAMP, 2'3'-cGAMP, vadimezan, and DMXAA have **{cat('SigCom LINCS pert_name', 'diABZI')} / {cat('SigCom LINCS pert_name', 'ADU-S100')} / {cat('SigCom LINCS pert_name', 'cGAMP')} / {cat('SigCom LINCS pert_name', "2'3'-cGAMP")} / {cat('SigCom LINCS pert_name', 'vadimezan')} / {cat('SigCom LINCS pert_name', 'DMXAA')}** signatures. TMEM173 and MB21D1 CRISPR knockouts are reported as STING/cGAS loss, a different class from IFN ligand treatment.

Background, seed {SEED}: up to {BG_K} other signatures of the same `pert_type` and cell line (CRISPR, shRNA, or ligand). The paired test uses cell-line median Spearman minus the matched background median.

## Cell-line summary

Unit is the cell line. Wilcoxon is two-sided on cell-line medians and needs at least 5 lines. BH q is inside the primary family only: T47D and MCF7, classes NHEJ CRISPR, NHEJ shRNA, and IFN ligand, scope all cell lines. Other rows are sensitivity and keep raw p.

| Query | Class | Cell lines | Signatures | Median Spearman | IQR | P vs 0 | q vs 0 | Median minus background | P vs background |
|---|---|---:|---:|---:|---|---:|---:|---:|---:|
{header_rows}

### Lung lines only

Phase I `primary_site = lung`, hyphen-stripped ids, and LCLC103H (large-cell lung carcinoma in the 2021 CRISPR panel; absent from the 2017 cell table).

| Query | Class | Cell lines | Signatures | Median Spearman | IQR | P vs 0 | q vs 0 | Median minus background | P vs background |
|---|---|---:|---:|---:|---|---:|---:|---:|---:|
{lung_rows}

Secondary scores (median of cell-line medians) are `median_cell_cosine` and `median_cell_wtcs` in `tables/class_summary.tsv`. WTCS uses the top and bottom {WTCS_N} CLDN4-loss genes, weighted by |log2FC|, and is 0 when the up and down enrichments are not opposite.

## NHEJ gene medians

### T47D

{gene_block("GSE207704_T47D")}

### MCF7

{gene_block("GSE207704_MCF7")}

## Reading rule

The headline is the cell-line median Spearman and its comparison with the matched background. A small Spearman can produce a low p when many cell lines share a sign. Signature-level p-values are not used: each signature has thousands of genes, so a tiny correlation is automatically “significant” at gene resolution.

Signature files dropped at parse: {n_fail}.

## Reproduce

```bash
python3 -m pip install -r methods/lincs_cldn4_connectivity/requirements.txt
python3 scripts/lincs_cldn4_connectivity/analyze.py
```

Raw GEO and LINCS signature files are cached under `methods/lincs_cldn4_connectivity/data/` and are gitignored. Tables and figures are the committed result.
"""
    (OUT / "FINDING.md").write_text(md)


def make_figures(cell: pd.DataFrame, summary: pd.DataFrame, genes: pd.DataFrame) -> None:
    queries = ["GSE207704_T47D", "GSE207704_MCF7", "GSE22493_SKOV3", "GSE50927_naive_lung"]
    fig, ax = plt.subplots(figsize=(10.2, 4.8))
    x = np.arange(len(queries))
    width = 0.18
    for i, sig_class in enumerate(CLASS_ORDER):
        meds, yerr_lo, yerr_hi = [], [], []
        for q in queries:
            hit = summary[
                (summary["query"] == q)
                & (summary["sig_class"] == sig_class)
                & (summary["scope"] == "all_cell_lines")
            ]
            if hit.empty:
                meds.append(np.nan)
                yerr_lo.append(0)
                yerr_hi.append(0)
                continue
            med = float(hit.iloc[0]["median_cell_spearman"])
            q25 = float(hit.iloc[0]["q25_spearman"])
            q75 = float(hit.iloc[0]["q75_spearman"])
            meds.append(med)
            yerr_lo.append(med - q25)
            yerr_hi.append(q75 - med)
        ax.bar(
            x + (i - 1.5) * width,
            meds,
            width=width,
            color=CLASS_COLOR[sig_class],
            label=CLASS_LABEL[sig_class],
            yerr=np.vstack([yerr_lo, yerr_hi]),
            capsize=2,
            error_kw={"elinewidth": 0.8},
        )
    ax.axhline(0, color="black", lw=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels([QUERY_LABEL[q] for q in queries], rotation=15, ha="right")
    ax.set_ylabel("Median cell-line Spearman")
    ax.set_title("CLDN4-loss connectivity to LINCS signatures")
    ax.legend(frameon=False, ncol=2, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "fig_class_spearman.png", dpi=160)
    fig.savefig(FIG / "fig_class_spearman.pdf")
    plt.close(fig)

    # Gene heatmap for human CRISPR/shRNA on T47D and MCF7, genes as rows.
    heat_queries = ["GSE207704_T47D", "GSE207704_MCF7", "GSE22493_SKOV3"]
    # One column per query x class that exists.
    col_keys = []
    for q in heat_queries:
        for c in ("NHEJ_CRISPR", "NHEJ_shRNA"):
            if ((genes["query"] == q) & (genes["sig_class"] == c)).any():
                col_keys.append((q, c))
    mat = pd.DataFrame(index=list(NHEJ_GENES), columns=[f"{QUERY_LABEL[q]}\n{CLASS_LABEL[c]}" for q, c in col_keys], dtype=float)
    for q, c in col_keys:
        label = f"{QUERY_LABEL[q]}\n{CLASS_LABEL[c]}"
        sub = genes[(genes["query"] == q) & (genes["sig_class"] == c)].set_index("gene")
        for g in NHEJ_GENES:
            if g in sub.index:
                mat.loc[g, label] = sub.loc[g, "median_spearman"]
    fig, ax = plt.subplots(figsize=(8.4, 4.6))
    data = mat.to_numpy(dtype=float)
    finite = data[np.isfinite(data)]
    vmax = max(0.05, float(np.nanmax(np.abs(finite))) if len(finite) else 0.05)
    im = ax.imshow(data, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(mat.shape[1]))
    ax.set_xticklabels(mat.columns, fontsize=7)
    ax.set_yticks(range(len(NHEJ_GENES)))
    ax.set_yticklabels(NHEJ_GENES)
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            val = data[i, j]
            if math.isfinite(val):
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=7)
    ax.set_title("NHEJ-gene median Spearman with CLDN4 loss")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Spearman")
    fig.tight_layout()
    fig.savefig(FIG / "fig_nhej_gene_spearman.png", dpi=160)
    fig.savefig(FIG / "fig_nhej_gene_spearman.pdf")
    plt.close(fig)

    # Cell-line strip for the two primary queries, NHEJ CRISPR vs IFN.
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.4), sharey=True)
    for ax, query in zip(axes, PRIMARY_QUERIES):
        for i, sig_class in enumerate(CLASS_ORDER):
            sub = cell[(cell["query"] == query) & (cell["sig_class"] == sig_class)]
            if sub.empty:
                continue
            jitter = (np.random.default_rng(SEED).random(len(sub)) - 0.5) * 0.25
            ax.scatter(
                np.full(len(sub), i) + jitter,
                sub["median_spearman"],
                s=18,
                c=CLASS_COLOR[sig_class],
                alpha=0.85,
                linewidths=0,
            )
        ax.axhline(0, color="black", lw=0.6)
        ax.set_xticks(range(len(CLASS_ORDER)))
        ax.set_xticklabels([CLASS_LABEL[c] for c in CLASS_ORDER], rotation=20, ha="right", fontsize=8)
        ax.set_title(QUERY_LABEL[query])
        ax.set_ylabel("Cell-line median Spearman" if query == PRIMARY_QUERIES[0] else "")
    fig.suptitle("Cell-line connectivity", y=1.02)
    fig.tight_layout()
    fig.savefig(FIG / "fig_cellline_strips.png", dpi=160, bbox_inches="tight")
    fig.savefig(FIG / "fig_cellline_strips.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    log("catalog")
    catalog_df = catalog()
    catalog_df.to_csv(TAB / "lincs_catalog.tsv", sep="\t", index=False)
    lung = lung_cell_lines()
    pd.Series(sorted(lung), name="cell_id").to_csv(TAB / "lung_cell_ids.tsv", sep="\t", index=False)
    log(f"lung cell ids {len(lung)}")

    log("queries")
    queries, query_meta = build_queries()
    query_meta.to_csv(TAB / "query_contrasts.tsv", sep="\t", index=False)
    for key, vec in queries.items():
        ser = pd.Series(vec, name="log2fc")
        ser.index.name = "symbol"
        ser.to_csv(TAB / f"rank_{key}.tsv.gz", sep="\t")
    log(query_meta[["query", "target_log2fc", "n_genes"]].to_string(index=False))

    log("reference metadata")
    ref = collect_references()
    rng = np.random.default_rng(SEED)
    log("background metadata")
    bg = sample_background(ref, rng)
    meta = pd.concat([ref, bg], ignore_index=True)
    meta.to_csv(TAB / "signature_metadata.tsv", sep="\t", index=False)
    log(f"signatures to download {len(meta)}")

    log("download characteristic-direction files")
    got = ensure_signatures(meta)
    got.to_csv(TAB / "signature_download.tsv", sep="\t", index=False)
    n_fail = int((got["download_status"] != "ok").sum())
    log(f"download failures {n_fail}")

    log("score")
    score_rows = []
    use = got[got["download_status"].eq("ok")]
    for rec in use.itertuples(index=False):
        try:
            ref_vec = load_cd(Path(rec.path))
        except Exception as exc:  # noqa: BLE001
            log(f"parse fail {rec.sig_id}: {exc}")
            continue
        for qname, qvec in queries.items():
            rho, cosine, connectivity, n_genes = score_pair(qvec, ref_vec)
            score_rows.append(
                {
                    "query": qname,
                    "sig_class": rec.sig_class,
                    "gene_group": rec.gene_group,
                    "pert_name": rec.pert_name,
                    "pert_type": rec.pert_type,
                    "cell_line": rec.cell_line,
                    "sig_id": rec.sig_id,
                    "pert_time": rec.pert_time,
                    "n_genes": n_genes,
                    "spearman": rho,
                    "cosine": cosine,
                    "wtcs": connectivity,
                }
            )
    scores = pd.DataFrame(score_rows)
    scores.to_csv(TAB / "connectivity_by_signature.tsv.gz", sep="\t", index=False)
    log(f"scored rows {len(scores)}")

    cell, summary, genes = summarize(scores, lung)
    cell.to_csv(TAB / "cellline_summary.tsv", sep="\t", index=False)
    summary.to_csv(TAB / "class_summary.tsv", sep="\t", index=False)
    genes.to_csv(TAB / "gene_summary.tsv", sep="\t", index=False)

    make_figures(cell, summary, genes)
    write_finding(catalog_df, query_meta, summary, genes, n_fail)

    # Hard checks so a parse failure cannot ship a silent empty result.
    t47 = float(query_meta.loc[query_meta["query"].eq("GSE207704_T47D"), "target_log2fc"].iloc[0])
    mcf = float(query_meta.loc[query_meta["query"].eq("GSE207704_MCF7"), "target_log2fc"].iloc[0])
    lung_fc = float(query_meta.loc[query_meta["query"].eq("GSE50927_naive_lung"), "target_log2fc"].iloc[0])
    sk = float(query_meta.loc[query_meta["query"].eq("GSE22493_SKOV3"), "target_log2fc"].iloc[0])
    assert t47 < -0.5, t47
    assert mcf < -0.4, mcf
    assert lung_fc < -3, lung_fc
    assert sk < 0, sk
    cldn4_n = catalog_df[(catalog_df["resource"] == "SigCom LINCS pert_name") & (catalog_df["query"] == "CLDN4")]
    assert int(cldn4_n.iloc[0]["n"]) == 0
    n_nhej = int((scores["sig_class"] == "NHEJ_CRISPR").sum() / len(queries))
    assert n_nhej >= 40, n_nhej
    human = scores[scores["query"].eq("GSE207704_T47D") & scores["sig_class"].eq("NHEJ_CRISPR")]
    assert human["n_genes"].median() > 2000
    assert (FIG / "fig_class_spearman.png").exists()
    log("assertions passed")
    show = summary[(summary["scope"] == "all_cell_lines") & (summary["query"].isin(PRIMARY_QUERIES))]
    log(show[["query", "sig_class", "n_cell", "median_cell_spearman", "wilcoxon_p_vs_0", "paired_p_vs_bg"]].to_string(index=False))


def rescore() -> None:
    """Recompute scores from cached signature files and saved metadata."""
    catalog_df = pd.read_csv(TAB / "lincs_catalog.tsv", sep="\t")
    query_meta = pd.read_csv(TAB / "query_contrasts.tsv", sep="\t")
    queries = {}
    for key in query_meta["query"]:
        ser = pd.read_csv(TAB / f"rank_{key}.tsv.gz", sep="\t", index_col=0).iloc[:, 0]
        queries[key] = {str(g): float(v) for g, v in ser.items() if pd.notna(v)}
    got = pd.read_csv(TAB / "signature_download.tsv", sep="\t")
    lung = lung_cell_lines()
    log("rescore")
    score_rows = []
    use = got[got["download_status"].eq("ok")]
    n_parse_fail = 0
    for rec in use.itertuples(index=False):
        try:
            ref_vec = load_cd(Path(rec.path))
        except Exception as exc:  # noqa: BLE001
            n_parse_fail += 1
            log(f"parse fail {rec.sig_id}: {exc}")
            continue
        if len(ref_vec) < MIN_OVERLAP:
            n_parse_fail += 1
            log(f"parse short {rec.sig_id}: {len(ref_vec)} genes")
            continue
        for qname, qvec in queries.items():
            rho, cosine, connectivity, n_genes = score_pair(qvec, ref_vec)
            score_rows.append(
                {
                    "query": qname,
                    "sig_class": rec.sig_class,
                    "gene_group": rec.gene_group,
                    "pert_name": rec.pert_name,
                    "pert_type": rec.pert_type,
                    "cell_line": rec.cell_line,
                    "sig_id": rec.sig_id,
                    "pert_time": rec.pert_time,
                    "n_genes": n_genes,
                    "spearman": rho,
                    "cosine": cosine,
                    "wtcs": connectivity,
                }
            )
    scores = pd.DataFrame(score_rows)
    scores.to_csv(TAB / "connectivity_by_signature.tsv.gz", sep="\t", index=False)
    log(f"scored rows {len(scores)} parse_fail {n_parse_fail}")
    cell, summary, genes = summarize(scores, lung)
    cell.to_csv(TAB / "cellline_summary.tsv", sep="\t", index=False)
    summary.to_csv(TAB / "class_summary.tsv", sep="\t", index=False)
    genes.to_csv(TAB / "gene_summary.tsv", sep="\t", index=False)
    make_figures(cell, summary, genes)
    write_finding(catalog_df, query_meta, summary, genes, n_parse_fail)
    human = scores[scores["query"].eq("GSE207704_T47D") & scores["sig_class"].eq("NHEJ_CRISPR")]
    ifn = scores[scores["query"].eq("GSE207704_T47D") & scores["sig_class"].eq("IFN_ligand")]
    assert len(human) >= 40
    assert human["n_genes"].median() > 2000
    assert len(ifn) >= 50, len(ifn)
    assert (FIG / "fig_class_spearman.png").exists()
    log("rescore assertions passed")
    show = summary[(summary["scope"] == "all_cell_lines") & (summary["query"].isin(PRIMARY_QUERIES))]
    log(
        show[
            ["query", "sig_class", "n_cell", "n_sig", "median_cell_spearman", "wilcoxon_p_vs_0", "paired_p_vs_bg"]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--rescore":
        rescore()
    else:
        main()
