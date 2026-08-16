"""I/O helpers for CPTAC/TMT gene-abundance, GDC STAR, LinkedOmics CCT, and PRIDE tables."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

TARGET_GENES = {
    "TACSTD2": {
        "ensembl": "ENSG00000184292",
        "symbol": "TACSTD2",
        "uniprot": "P09758",
        "aliases": ("TROP2", "EGP-1", "GA733-1", "M1S1", "TACSTD2"),
    },
    "CLDN4": {
        "ensembl": "ENSG00000189143",
        "symbol": "CLDN4",
        "uniprot": "O14493",
        "aliases": ("CLDN4", "CLD4", "CPE-R", "CPETR1", "WBSCR8"),
    },
}

ENSEMBL_TO_SYMBOL = {
    "ENSG00000184292": "TACSTD2",
    "ENSG00000189143": "CLDN4",
}
SYMBOL_TO_ENSEMBL = {v: k for k, v in ENSEMBL_TO_SYMBOL.items()}

_DOT_ID = re.compile(r"[.\-]")


def strip_ensembl_version(gene_id: str) -> str:
    s = str(gene_id).strip()
    if s.startswith("ENSG") and "." in s:
        return s.split(".", 1)[0]
    return s


def normalize_gene_id(label: str) -> str:
    """Map Ensembl / symbol / UniProt / alias to a canonical symbol when known."""
    raw = str(label).strip()
    if not raw or raw.lower() in {"nan", "na", "none"}:
        return raw
    ens = strip_ensembl_version(raw)
    if ens in ENSEMBL_TO_SYMBOL:
        return ENSEMBL_TO_SYMBOL[ens]
    up = raw.upper()
    for symbol, meta in TARGET_GENES.items():
        if up in {a.upper() for a in meta["aliases"]} or up == meta["uniprot"]:
            return symbol
    return raw


def normalize_sample_id(sample_id: str) -> str:
    """C3L.00001 / C3L-00001 / C3L_00001 -> C3L-00001 (CPTAC submitter style)."""
    s = str(sample_id).strip()
    s = s.replace(".rna_seq.augmented_star_gene_counts.tsv", "")
    s = re.sub(r"^X(?=\d)", "", s)  # LinkedOmics sometimes prefixes numeric IDs
    m = re.match(r"^(C3[LN])[.\-_]?(\d+)", s, flags=re.I)
    if m:
        return f"{m.group(1).upper()}-{m.group(2)}"
    return s.replace(".", "-")


def normalize_sample_columns(columns: Iterable) -> list[str]:
    return [normalize_sample_id(c) for c in columns]


def read_table(path: Path | str, index_col: int | str | None = 0) -> pd.DataFrame:
    path = Path(path)
    sep = "\t"
    if path.suffix.lower() in {".csv"}:
        sep = ","
    # LinkedOmics .tsi clinical: first column is attribute name, samples are columns
    df = pd.read_csv(path, sep=sep, index_col=index_col, comment=None, dtype=str)
    if df.index.name and str(df.index.name).lower() in {"attrib_name", "sample.id", "sample_id"}:
        pass
    return df


def _to_numeric_matrix(df: pd.DataFrame) -> pd.DataFrame:
    out = df.apply(pd.to_numeric, errors="coerce")
    out.index = [str(i) for i in out.index]
    out.columns = normalize_sample_columns(out.columns)
    # drop completely empty rows/cols created by headers
    out = out.dropna(axis=0, how="all").dropna(axis=1, how="all")
    return out


def read_expression_matrix(path: Path | str) -> pd.DataFrame:
    """Read a genes x samples numeric matrix (CCT, gene-abundance TSV, generic)."""
    path = Path(path)
    raw = pd.read_csv(path, sep="\t", index_col=0)
    # LinkedOmics clinical/phenotype .tsi stores data_type on first row of values
    if raw.index.name and "data_type" in {str(x).lower() for x in raw.index[:3]}:
        raw = raw.drop(index=[i for i in raw.index if str(i).lower() in {"data_type", "nan"}])
    if raw.shape[1] == 0:
        raw = pd.read_csv(path, sep=",", index_col=0)
    return _to_numeric_matrix(raw)


def read_gdc_star_gene_counts(path: Path | str) -> pd.DataFrame:
    """Parse one GDC STAR gene-counts TSV (public S3 / GDC data API).

    Returns a single-column frame (sample x metrics are not stacked):
    index = gene_id (version stripped), columns = [gene_name, unstranded, tpm, fpkm, fpkm_uq].
    """
    path = Path(path)
    df = pd.read_csv(path, sep="\t", comment="#")
    df = df[df["gene_id"].astype(str).str.startswith("ENSG")].copy()
    df["ensembl"] = df["gene_id"].map(strip_ensembl_version)
    df["gene_name"] = df["gene_name"].astype(str)
    keep = {
        "ensembl": "ensembl",
        "gene_name": "gene_name",
        "unstranded": "counts",
        "tpm_unstranded": "tpm",
        "fpkm_unstranded": "fpkm",
        "fpkm_uq_unstranded": "fpkm_uq",
    }
    out = df[list(keep)].rename(columns=keep)
    for c in ("counts", "tpm", "fpkm", "fpkm_uq"):
        out[c] = pd.to_numeric(out[c], errors="coerce")
    return out.set_index("ensembl")


def extract_star_targets(star_df: pd.DataFrame, sample_id: str) -> pd.DataFrame:
    rows = []
    for ens, symbol in ENSEMBL_TO_SYMBOL.items():
        if ens not in star_df.index:
            rows.append(
                {
                    "sample_id": normalize_sample_id(sample_id),
                    "ensembl": ens,
                    "symbol": symbol,
                    "counts": np.nan,
                    "tpm": np.nan,
                    "fpkm": np.nan,
                    "fpkm_uq": np.nan,
                    "found": False,
                }
            )
            continue
        r = star_df.loc[ens]
        if isinstance(r, pd.DataFrame):
            r = r.iloc[0]
        rows.append(
            {
                "sample_id": normalize_sample_id(sample_id),
                "ensembl": ens,
                "symbol": symbol,
                "counts": r.get("counts", np.nan),
                "tpm": r.get("tpm", np.nan),
                "fpkm": r.get("fpkm", np.nan),
                "fpkm_uq": r.get("fpkm_uq", np.nan),
                "found": True,
            }
        )
    return pd.DataFrame(rows)


def lookup_gene_rows(matrix: pd.DataFrame, query: str) -> pd.DataFrame:
    """Return rows matching symbol, Ensembl, or alias. Empty if the gene is absent (CLDN4 case)."""
    q = normalize_gene_id(query)
    ens = SYMBOL_TO_ENSEMBL.get(q, strip_ensembl_version(query))
    idx = pd.Index([str(i) for i in matrix.index])
    mask = idx.map(normalize_gene_id) == q
    mask = mask | (idx.map(strip_ensembl_version) == ens)
    mask = mask | idx.str.upper().isin({query.upper(), q.upper(), ens.upper()})
    return matrix.loc[mask]


def gene_presence(matrix: pd.DataFrame, query: str) -> dict:
    rows = lookup_gene_rows(matrix, query)
    if rows.empty:
        return {
            "query": query,
            "canonical": normalize_gene_id(query),
            "present": False,
            "n_rows": 0,
            "n_samples": matrix.shape[1],
            "n_observed": 0,
            "n_na": matrix.shape[1],
            "na_frac": 1.0,
            "note": "gene absent from table (not a row of NAs — never quantified / dropped by NArm)",
        }
    s = rows.iloc[0]
    n_obs = int(s.notna().sum())
    n = int(s.shape[0])
    return {
        "query": query,
        "canonical": normalize_gene_id(query),
        "present": True,
        "n_rows": int(rows.shape[0]),
        "n_samples": n,
        "n_observed": n_obs,
        "n_na": n - n_obs,
        "na_frac": float((n - n_obs) / n) if n else 1.0,
        "note": "row present",
    }


def detect_log_scale(values: pd.Series | np.ndarray, name: str = "matrix") -> dict:
    """Decide whether a table is already log2 (TMT ratio) or raw intensity/counts."""
    x = pd.to_numeric(pd.Series(np.asarray(values).ravel()), errors="coerce").dropna()
    if x.empty:
        return {"name": name, "already_log": False, "reason": "empty", "n": 0}
    xmin, xmax, med, frac_neg = float(x.min()), float(x.max()), float(x.median()), float((x < 0).mean())
    # TMT log2-ratio: negatives, |values| typically < 15
    if frac_neg > 0.01 and xmax < 40:
        return {
            "name": name,
            "already_log": True,
            "transform": "none",
            "reason": "negative values and compact range typical of TMT log2-ratio",
            "n": int(x.size),
            "min": xmin,
            "max": xmax,
            "median": med,
            "frac_neg": frac_neg,
        }
    # GDC / LinkedOmics RNA already log2(FPKM) or log2(UQ)
    if xmin >= -8 and xmax < 30 and med < 20 and frac_neg < 0.5:
        # could still be raw TPM (0-20). Use IQR and integer-ish counts.
        unique_ratio = x.nunique() / max(len(x), 1)
        if xmin >= 0 and med > 50:
            return {
                "name": name,
                "already_log": False,
                "transform": "log2(x+1)",
                "reason": "non-negative, large median — treat as raw abundance",
                "n": int(x.size),
                "min": xmin,
                "max": xmax,
                "median": med,
                "frac_neg": frac_neg,
            }
        return {
            "name": name,
            "already_log": True,
            "transform": "none",
            "reason": "compact signed/small-positive range; do not double-log",
            "n": int(x.size),
            "min": xmin,
            "max": xmax,
            "median": med,
            "frac_neg": frac_neg,
            "unique_ratio": unique_ratio,
        }
    if xmin >= 0 and xmax > 100:
        return {
            "name": name,
            "already_log": False,
            "transform": "log2(x+1)",
            "reason": "raw intensity/counts (wide non-negative range)",
            "n": int(x.size),
            "min": xmin,
            "max": xmax,
            "median": med,
            "frac_neg": frac_neg,
        }
    return {
        "name": name,
        "already_log": True,
        "transform": "none",
        "reason": "ambiguous; default to already-logged to avoid double-log",
        "n": int(x.size),
        "min": xmin,
        "max": xmax,
        "median": med,
        "frac_neg": frac_neg,
    }


def apply_log2_if_needed(matrix: pd.DataFrame, name: str = "matrix") -> tuple[pd.DataFrame, dict]:
    info = detect_log_scale(matrix.stack(dropna=True), name=name)
    if info.get("already_log", True):
        return matrix.copy(), info
    logged = np.log2(matrix.clip(lower=0) + 1.0)
    info["applied"] = "log2(x+1)"
    return logged, info


def read_pride_protein_table(path: Path | str) -> pd.DataFrame:
    """Best-effort parser for MaxQuant proteinGroups, DIA-NN pg_matrix, or generic PRIDE TSV."""
    path = Path(path)
    peek = pd.read_csv(path, sep="\t", nrows=0)
    cols = list(peek.columns)
    lower = {c.lower(): c for c in cols}

    if "majority protein ids" in lower or "gene names" in lower:
        df = pd.read_csv(path, sep="\t")
        if "Reverse" in df.columns:
            df = df[df["Reverse"].fillna("") != "+"]
        if "Potential contaminant" in df.columns:
            df = df[df["Potential contaminant"].fillna("") != "+"]
        gene_col = lower.get("gene names", lower.get("gene names", cols[0]))
        intensity_cols = [
            c
            for c in df.columns
            if c.startswith(("LFQ intensity", "Intensity ", "iBAQ ", "Reporter intensity corrected"))
            and not c.endswith("___")
        ]
        if not intensity_cols:
            intensity_cols = [c for c in df.columns if re.search(r"intensity|ibaq|lfq", c, re.I)]
        out = df.set_index(df[gene_col].astype(str).str.split(";").str[0])[intensity_cols]
        out.columns = [re.sub(r"^(LFQ intensity|Intensity|iBAQ|Reporter intensity corrected)\s*", "", c) for c in out.columns]
        return _to_numeric_matrix(out)

    if any("pg.quantity" in c.lower() or c.endswith(".raw") for c in cols) or "protein.group" in lower:
        df = pd.read_csv(path, sep="\t")
        id_col = lower.get("genes", lower.get("protein.group", cols[0]))
        value_cols = [c for c in df.columns if c not in {df.columns[0], id_col} and not c.lower().startswith("protein")]
        out = df.set_index(df[id_col].astype(str).str.split(";").str[0])[value_cols]
        return _to_numeric_matrix(out)

    return read_expression_matrix(path)


def read_cibersort(path: Path | str) -> pd.DataFrame:
    """CIBERSORT/CIBERSORTx: samples as rows, cell types as columns (Mixture first)."""
    path = Path(path)
    df = pd.read_csv(path, sep=None, engine="python")
    mix = None
    for c in df.columns:
        if str(c).lower() in {"mixture", "sample", "sample_id", "id"}:
            mix = c
            break
    if mix is None:
        mix = df.columns[0]
    df = df.set_index(mix)
    df.index = [normalize_sample_id(i) for i in df.index]
    numeric = df.apply(pd.to_numeric, errors="coerce")
    drop = [c for c in numeric.columns if re.search(r"p-value|rmse|correlation|absolute", str(c), re.I)]
    return numeric.drop(columns=drop, errors="ignore")


def read_xcell(path: Path | str) -> pd.DataFrame:
    """xCell: cell types as rows, samples as columns (or the transpose). Returns samples x cell types."""
    path = Path(path)
    df = pd.read_csv(path, sep="\t", index_col=0)
    # heuristic: more cell-type-like if first dim < 80 and second dim looks like samples
    if df.shape[0] < df.shape[1] and df.shape[0] <= 120:
        df = df.T
    df.index = [normalize_sample_id(i) for i in df.index]
    return df.apply(pd.to_numeric, errors="coerce")


def read_linkedomics_phenotype(path: Path | str) -> pd.DataFrame:
    """LinkedOmics .tsi: attributes as columns after transpose (samples x phenotypes)."""
    raw = pd.read_csv(path, sep="\t", index_col=0)
    raw = raw.loc[[i for i in raw.index if str(i).lower() not in {"nan", "data_type", "none"}]]
    # first row may be data_type codes (ORD/BIN/CAT) stored as a column named data_type
    if "data_type" in raw.columns:
        raw = raw.drop(columns=["data_type"])
    # file is attributes x samples OR samples x attributes
    if raw.shape[0] < raw.shape[1] and raw.shape[0] <= 80:
        pheno = raw.T
    else:
        # already samples x attributes (LSCC clinical uses sample IDs as index)
        pheno = raw
    pheno.index = [normalize_sample_id(i) for i in pheno.index]
    pheno.columns = [str(c).replace(".", "_") for c in pheno.columns]
    return pheno
