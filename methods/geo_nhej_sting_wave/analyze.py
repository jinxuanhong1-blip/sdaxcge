#!/usr/bin/env python3
"""GEO brute inventory: open matrices scored for NHEJ down and IFN up versus CLDN4."""

from __future__ import annotations

import argparse
import gzip
import io
import json
import math
import os
import re
import sys
import threading
import time
import traceback
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent
TABLES = ROOT / "tables"
FIGS = ROOT / "figures"
CACHE = Path("/tmp/geo_nhej_sting_cache")
GENE_SETS = json.loads((ROOT / "gene_sets.json").read_text())

NHEJ = list(GENE_SETS["nhej_core"])
IFN = list(GENE_SETS["ifn_isg"])
SENSORS = list(GENE_SETS["sting_sensors"])
ALIASES = {k.upper(): v for k, v in GENE_SETS["aliases"].items()}
MIN_NHEJ = int(GENE_SETS["min_nhej"])
MIN_IFN = int(GENE_SETS["min_ifn"])
MIN_N = int(GENE_SETS["min_samples_spearman"])
TARGETS = set(NHEJ + IFN + SENSORS + ["CLDN4"])
MAX_BYTES = 25 * 1024 * 1024
MIN_BYTES = 8 * 1024
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
FTP = "https://ftp.ncbi.nlm.nih.gov/geo"
UA = "sdaxcge-geo-nhej-sting/1.0"

# Download these query families. interferon_gse is count-only.
DOWNLOAD_QUERIES = {
    "brute_token",
    "cldn4_text",
    "nhej",
    "prkdc",
    "sting1",
    "cgas_title",
    "dna_repair",
}
COUNT_ONLY = {"strict_user", "interferon_gse", "sting_mesh", "dna_repair_title"}

QUERIES = [
    (
        "strict_user",
        '(CLDN4 OR Cldn4) AND ("DNA repair" OR NHEJ OR PRKDC OR STING OR cGAS OR interferon) AND gse[Entry Type]',
    ),
    (
        "brute_token",
        '(CLDN4 OR Cldn4 OR "claudin-4" OR "claudin 4" OR claudin OR CLDN) AND '
        '(NHEJ OR PRKDC OR STING1 OR TMEM173 OR "stimulator of interferon genes" OR '
        '"cGAS" OR "cyclic GMP-AMP" OR "cGAS-STING" OR "DNA repair" OR interferon) AND gse[Entry Type]',
    ),
    (
        "cldn4_text",
        '(CLDN4 OR Cldn4 OR "claudin-4" OR "claudin 4") AND gse[Entry Type]',
    ),
    ("nhej", "NHEJ AND gse[Entry Type]"),
    ("prkdc", "PRKDC AND gse[Entry Type]"),
    (
        "sting1",
        '(STING1 OR TMEM173 OR "stimulator of interferon genes" OR "cGAS-STING") AND gse[Entry Type]',
    ),
    (
        "cgas_title",
        '("cGAS"[Title] OR "cyclic GMP-AMP synthase"[Title] OR "cyclic GMP-AMP"[Title] OR MB21D1[Title]) AND gse[Entry Type]',
    ),
    ("dna_repair", '"DNA repair" AND gse[Entry Type]'),
    ("dna_repair_title", '"DNA repair"[Title] AND gse[Entry Type]'),
    ("interferon_gse", "interferon AND gse[Entry Type]"),
    ("sting_mesh", "STING AND gse[Entry Type]"),
]

LOW_RE = re.compile(
    r"(knock[\s_-]*down|knock[\s_-]*out|silenc\w*|crispri|"
    r"\bko\b|\bnull\b|heterozyg\w*)",
    re.I,
)
HIGH_RE = re.compile(
    r"(wild[\s_-]*type|\bwt\b|over[\s_-]*express\w*|\bparental\b|\bcontrol\b|\boe\b)",
    re.I,
)
SKIP_TITLE = re.compile(
    r"\b(chip-?seq|atac-?seq|dnase-?seq|bisulfite|wgbs|mirna|microrna|"
    r"small rna|ribo-seq|hi-c|scatac|snatac|cite-seq|methylation array)\b",
    re.I,
)
EXPR_NAME = re.compile(
    r"(count|fpkm|tpm|expr|expression|norm|rnaseq|rna-seq|matrix|htseq|"
    r"featurecount|rsem|salmon|kallisto|gene)",
    re.I,
)
SKIP_FILE = re.compile(
    r"(raw\.tar|\.cel|\.idat|\.bam|\.bw$|bigwig|\.bed\.|\.hic|fastq|\.sff|"
    r"\.gpr|barcodes\.tsv|features\.tsv|matrix\.mtx)",
    re.I,
)


class Limiter:
    def __init__(self, interval=0.2):
        self.interval = interval
        self.lock = threading.Lock()
        self.next_t = 0.0

    def wait(self):
        with self.lock:
            now = time.time()
            if now < self.next_t:
                time.sleep(self.next_t - now)
            self.next_t = time.time() + self.interval


LIMIT = Limiter(0.2)
WRITE_LOCK = threading.Lock()
GPL_LOCK = threading.Lock()


def fetch(url, timeout=90, retries=5):
    last = None
    for attempt in range(retries):
        LIMIT.wait()
        req = urllib.request.Request(url, headers={"User-Agent": UA})
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except urllib.error.HTTPError as exc:
            last = exc
            if exc.code in (429, 500, 502, 503, 504):
                time.sleep(2 ** attempt + 1)
                continue
            raise
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(1.5 ** attempt)
    raise last


def fetch_json(url):
    return json.loads(fetch(url).decode())


def esearch_count(term):
    url = (
        f"{EUTILS}/esearch.fcgi?db=gds&retmode=json&retmax=0&tool=sdaxcge&term="
        + urllib.parse.quote(term)
    )
    data = fetch_json(url)["esearchresult"]
    return int(data["count"])


def esearch_ids(term):
    ids = []
    retstart = 0
    count = None
    while True:
        url = (
            f"{EUTILS}/esearch.fcgi?db=gds&retmode=json&retmax=500&retstart={retstart}"
            f"&tool=sdaxcge&term=" + urllib.parse.quote(term)
        )
        data = fetch_json(url)["esearchresult"]
        count = int(data["count"])
        batch = data.get("idlist") or []
        ids.extend(batch)
        if not batch or len(ids) >= count:
            break
        retstart += len(batch)
    # preserve order, drop dupes
    seen = set()
    out = []
    for uid in ids:
        if uid not in seen:
            seen.add(uid)
            out.append(uid)
    return count, out


def esummary(uids):
    records = []
    for i in range(0, len(uids), 80):
        batch = uids[i : i + 80]
        url = f"{EUTILS}/esummary.fcgi?db=gds&retmode=json&tool=sdaxcge&id=" + ",".join(batch)
        try:
            data = fetch_json(url)["result"]
        except Exception:
            # split once
            mid = len(batch) // 2
            for part in (batch[:mid], batch[mid:]):
                if not part:
                    continue
                url = f"{EUTILS}/esummary.fcgi?db=gds&retmode=json&tool=sdaxcge&id=" + ",".join(part)
                data = fetch_json(url)["result"]
                for uid in data.get("uids", []):
                    records.append(data[uid])
            continue
        for uid in data.get("uids", []):
            records.append(data[uid])
    return records


def geo_stub(acc):
    prefix = "".join(ch for ch in acc if ch.isalpha())
    num = acc[len(prefix) :]
    if len(num) <= 3:
        return prefix + "nnn"
    return prefix + num[:-3] + "nnn"


def parse_size(token):
    token = token.strip()
    if token in {"-", ""}:
        return None
    mult = { "K": 1024, "M": 1024**2, "G": 1024**3}
    if token[-1] in mult:
        return int(float(token[:-1]) * mult[token[-1]])
    try:
        return int(token)
    except ValueError:
        return None


def list_ftp(url):
    try:
        html = fetch(url, timeout=40, retries=3).decode("utf-8", "replace")
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            return []
        raise
    files = []
    for name, size in re.findall(
        r'href="([^"/][^"]*)">.*?</a>\s+\d{4}-\d{2}-\d{2}\s+\d{2}:\d{2}\s+(\S+)',
        html,
    ):
        files.append({"name": name, "size": parse_size(size), "url": url + name})
    return files


def norm_symbol(raw):
    if raw is None or (isinstance(raw, float) and math.isnan(raw)):
        return None
    text = str(raw).strip()
    if not text or text.lower() in {"nan", "none", "---", "na"}:
        return None
    # Ensembl with version, or a bare symbol, or symbol///symbol
    parts = re.split(r"\s*(?://+|///|\||;|,)\s*", text)
    found = []
    for part in parts:
        token = part.strip()
        if not token:
            continue
        token = token.split()[0]
        if re.fullmatch(r"\d+\.0", token):
            token = token[:-2]
        if token.upper().startswith("ENS"):
            token = token.split(".")[0]
            canon = ID_TO_CANON.get(token)
            if canon and canon not in found:
                found.append(canon)
            continue
        token = token.upper()
        token = ALIASES.get(token, token)
        if token in TARGETS and token not in found:
            found.append(token)
        else:
            canon = ID_TO_CANON.get(token)
            if canon and canon not in found:
                found.append(canon)
    if len(found) == 1:
        return found[0]
    return None


def build_id_map():
    mapping = {}
    for canon, rec in GENE_SETS["ids"].items():
        mapping[canon] = canon
        for ens in rec.get("ensembl", []):
            mapping[ens] = canon
        for ent in rec.get("entrez", []):
            mapping[str(ent)] = canon
    for alias, canon in ALIASES.items():
        mapping[alias] = canon
    return mapping


ID_TO_CANON = build_id_map()


def arm_and_stem(name):
    text = str(name).lower()
    text = re.sub(r"\([^)]*\)", " ", text)
    text = text.replace("_", " ").replace("-", " ")
    # Split CLDN4 glued to an arm token (CLDN4KO) before the word-boundary tests.
    text = re.sub(r"cldn\s*4", " ", text)
    low = bool(LOW_RE.search(text))
    high = bool(HIGH_RE.search(text))
    stem = LOW_RE.sub(" ", text)
    stem = HIGH_RE.sub(" ", stem)
    stem = re.sub(r"\b(fpkm|tpm|counts|count|expression|norm|normalized)\b", " ", stem)
    stem = re.sub(r"[^a-z0-9]+", " ", stem)
    stem = " ".join(stem.split())
    if low and not high:
        arm = "low"
    elif high and not low:
        arm = "high"
    else:
        arm = None
    return arm, stem


def find_pairs(names):
    groups = defaultdict(lambda: {"high": [], "low": []})
    for i, name in enumerate(names):
        arm, stem = arm_and_stem(name)
        if arm:
            groups[stem][arm].append(i)
    pairs = []
    for stem, arms in groups.items():
        if arms["high"] and arms["low"]:
            pairs.append((stem, arms["high"], arms["low"]))
    if pairs:
        return pairs
    # drop trailing replicate numbers and try once
    groups = defaultdict(lambda: {"high": [], "low": []})
    for i, name in enumerate(names):
        arm, stem = arm_and_stem(name)
        if not arm:
            continue
        stem2 = " ".join(tok for tok in stem.split() if not tok.isdigit())
        groups[stem2][arm].append(i)
    pairs = []
    for stem, arms in groups.items():
        if arms["high"] and arms["low"]:
            pairs.append((stem or "all", arms["high"], arms["low"]))
    return pairs


def needs_log(matrix):
    values = matrix[np.isfinite(matrix)]
    if values.size == 0:
        return False
    if np.min(values) < -1e-6:
        return False
    if np.median(values) > 15:
        return True
    if np.percentile(values, 90) > 30:
        return True
    return False


def looks_like_log_ratio(matrix):
    values = matrix[np.isfinite(matrix)]
    if values.size < 30:
        return False
    return bool(
        np.mean(values < 0) > 0.15
        and abs(float(np.median(values))) < 1.0
        and float(np.percentile(values, 99)) < 30
    )


def collapse_rows(symbols, matrix):
    buckets = defaultdict(list)
    for i, symbol in enumerate(symbols):
        if symbol:
            buckets[symbol].append(matrix[i])
    out = {}
    for symbol, rows in buckets.items():
        stacked = np.vstack(rows)
        out[symbol] = np.nanmedian(stacked, axis=0)
    return out


def coverage_ok(expr):
    n_nhej = sum(g in expr for g in NHEJ)
    n_ifn = sum(g in expr for g in IFN)
    return n_nhej >= MIN_NHEJ and n_ifn >= MIN_IFN, n_nhej, n_ifn


def module_vector(expr, genes, zscore=True):
    rows = [expr[g] for g in genes if g in expr]
    if not rows:
        return None
    mat = np.vstack(rows)
    if not zscore:
        return np.nanmean(mat, axis=0)
    mu = np.nanmean(mat, axis=1, keepdims=True)
    sd = np.nanstd(mat, axis=1, ddof=1, keepdims=True)
    sd[~np.isfinite(sd) | (sd == 0)] = np.nan
    z = (mat - mu) / sd
    return np.nanmean(z, axis=0)


def safe_spearman(x, y):
    mask = np.isfinite(x) & np.isfinite(y)
    n = int(mask.sum())
    if n < MIN_N:
        return n, np.nan, np.nan
    if np.nanstd(x[mask]) == 0 or np.nanstd(y[mask]) == 0:
        return n, np.nan, np.nan
    rho, p = stats.spearmanr(x[mask], y[mask])
    return n, float(rho), float(p)


def paired_gene_deltas(expr, genes, pairs):
    deltas = []
    used = []
    for gene in genes:
        if gene not in expr:
            continue
        arr = expr[gene]
        ds = []
        for _stem, high, low in pairs:
            hi = np.nanmean(arr[high])
            lo = np.nanmean(arr[low])
            if np.isfinite(hi) and np.isfinite(lo):
                ds.append(hi - lo)
        if ds:
            deltas.append(float(np.mean(ds)))
            used.append(gene)
    if not deltas:
        return np.nan, 0, 0, []
    arr = np.array(deltas)
    return float(np.mean(arr)), int(np.sum(arr > 0)), int(np.sum(arr < 0)), used


def orient_ratio(expr):
    if "CLDN4" not in expr:
        return None
    mean = float(np.nanmean(expr["CLDN4"]))
    if mean < -0.05:
        return -1.0  # deposited values are low/high; flip to high minus low
    if mean > 0.05:
        return 1.0
    return None


def score_expr(expr, names, design):
    """expr values are already on the analysis scale (logged if needed)."""
    present_nhej = [g for g in NHEJ if g in expr]
    present_ifn = [g for g in IFN if g in expr]
    present_sens = [g for g in SENSORS if g in expr]
    ok, n_nhej, n_ifn = coverage_ok(expr)
    row = {
        "cldn4_present": "CLDN4" in expr,
        "n_nhej": n_nhej,
        "n_ifn": n_ifn,
        "n_sensors": len(present_sens),
        "genes_nhej": ",".join(present_nhej),
        "genes_ifn": ",".join(present_ifn),
        "genes_sensors": ",".join(present_sens),
        "coverage_ok": bool(ok and "CLDN4" in expr),
        "n_columns": int(len(names)),
        "design": design,
    }
    if "CLDN4" not in expr:
        row["skip_reason"] = "CLDN4_absent"
        return row
    cldn4 = expr["CLDN4"]
    row["cldn4_mean"] = float(np.nanmean(cldn4))
    row["cldn4_std"] = float(np.nanstd(cldn4))

    pairs = [] if design == "log_ratio" else find_pairs(names)
    row["n_pairs"] = len(pairs)
    row["pair_stems"] = ";".join(stem for stem, _h, _l in pairs)

    # A short log-ratio table is one contrast. A one-column logFC table is too,
    # even when the target-gene slice is too small for the ratio heuristic.
    one_logfc = len(names) == 1 and bool(re.search(r"log2?fc|logfc", str(names[0]), re.I))
    short_ratio = (design == "log_ratio" and len(names) < MIN_N) or one_logfc
    if one_logfc:
        row["design"] = "log_ratio"
        design = "log_ratio"
    if short_ratio:
        sign = orient_ratio(expr)
        row["ratio_orient"] = sign
        if sign is None:
            row["contrast_note"] = "log-ratio orientation ambiguous (CLDN4 mean near 0)"
        else:
            def module_mean(genes):
                vals = [float(np.nanmean(expr[g])) for g in genes if g in expr and np.isfinite(np.nanmean(expr[g]))]
                return float(sign * np.mean(vals)) if vals else np.nan

            row["delta_nhej"] = module_mean(NHEJ)
            row["delta_ifn"] = module_mean(IFN)
            row["delta_cgas"] = module_mean(["CGAS"]) if "CGAS" in expr else np.nan
            row["delta_sting"] = module_mean(["STING1"]) if "STING1" in expr else np.nan
            row["contrast_note"] = "short log-ratio or logFC table, oriented CLDN4-high minus CLDN4-low"
    elif pairs and ok:
        d_n, pos_n, neg_n, _ = paired_gene_deltas(expr, NHEJ, pairs)
        d_i, pos_i, neg_i, _ = paired_gene_deltas(expr, IFN, pairs)
        row["delta_nhej"] = d_n
        row["delta_ifn"] = d_i
        row["nhej_genes_up_in_high"] = pos_n
        row["nhej_genes_down_in_high"] = neg_n
        row["ifn_genes_up_in_high"] = pos_i
        row["ifn_genes_down_in_high"] = neg_i
        if "CGAS" in expr:
            d_c, _, _, _ = paired_gene_deltas(expr, ["CGAS"], pairs)
            row["delta_cgas"] = d_c
        if "STING1" in expr:
            d_s, _, _, _ = paired_gene_deltas(expr, ["STING1"], pairs)
            row["delta_sting"] = d_s
        d_c4, _, _, _ = paired_gene_deltas(expr, ["CLDN4"], pairs)
        row["delta_cldn4"] = d_c4
        row["contrast_note"] = "paired titles, high arm minus low arm, mean across genes of log values"

    # Spearman across samples, including multi-sample log-ratio series.
    if not short_ratio:
        nhej_score = module_vector(expr, NHEJ, zscore=True)
        ifn_score = module_vector(expr, IFN, zscore=True)
        if nhej_score is not None and ifn_score is not None and ok:
            composite = ifn_score - nhej_score
            n, rho, p = safe_spearman(cldn4, nhej_score)
            row["spearman_n"] = n
            row["rho_nhej"] = rho
            row["p_nhej"] = p
            _, rho_i, p_i = safe_spearman(cldn4, ifn_score)
            row["rho_ifn"] = rho_i
            row["p_ifn"] = p_i
            _, rho_c, p_c = safe_spearman(cldn4, composite)
            row["rho_composite"] = rho_c
            row["p_composite"] = p_c
            if "CGAS" in expr:
                _, rho_g, p_g = safe_spearman(cldn4, expr["CGAS"])
                row["rho_cgas"] = rho_g
                row["p_cgas"] = p_g
            if "STING1" in expr:
                _, rho_s, p_s = safe_spearman(cldn4, expr["STING1"])
                row["rho_sting"] = rho_s
                row["p_sting"] = p_s
    return row


def joint_from_signs(nhej_delta, ifn_delta, flat=0.0):
    if nhej_delta is None or ifn_delta is None:
        return ""
    try:
        nhej_delta = float(nhej_delta)
        ifn_delta = float(ifn_delta)
    except (TypeError, ValueError):
        return ""
    if not (math.isfinite(nhej_delta) and math.isfinite(ifn_delta)):
        return ""
    if abs(nhej_delta) <= flat and abs(ifn_delta) <= flat:
        return "flat"
    if abs(nhej_delta) <= flat or abs(ifn_delta) <= flat:
        return "mixed"
    nhej_down = nhej_delta < 0
    ifn_up = ifn_delta > 0
    if nhej_down and ifn_up:
        return "NHEJ_down_IFN_up"
    if (not nhej_down) and (not ifn_up):
        return "NHEJ_up_IFN_down"
    return "mixed"


def split_tsv(line):
    return [part.strip().strip('"') for part in line.rstrip("\n").split("\t")]


def parse_series_matrix(path):
    titles, accessions, chars = [], [], []
    platform = ""
    data_ids = []
    data_rows = []
    begun = False
    char_rows = []
    with gzip.open(path, "rt", errors="replace") as handle:
        for line in handle:
            if line.startswith("!Series_platform_id") and not platform:
                platform = split_tsv(line)[-1]
            elif line.startswith("!Sample_title\t"):
                titles = split_tsv(line)[1:]
            elif line.startswith("!Sample_geo_accession\t"):
                accessions = split_tsv(line)[1:]
            elif line.startswith("!Sample_characteristics"):
                char_rows.append(split_tsv(line)[1:])
            elif line.startswith("!series_matrix_table_begin"):
                begun = True
                header = split_tsv(next(handle))
                continue
            elif begun:
                if line.startswith("!series_matrix_table_end"):
                    break
                parts = split_tsv(line)
                if len(parts) < 2:
                    continue
                vals = []
                for item in parts[1:]:
                    try:
                        vals.append(float(item))
                    except ValueError:
                        vals.append(np.nan)
                data_ids.append(parts[0])
                data_rows.append(vals)
    if not data_rows:
        return None
    matrix = np.array(data_rows, dtype=float)
    names = titles or accessions or [f"s{i}" for i in range(matrix.shape[1])]
    # align characteristics
    meta = []
    width = matrix.shape[1]
    for i in range(width):
        bits = []
        if i < len(names):
            bits.append(names[i])
        for row in char_rows:
            if i < len(row) and row[i]:
                bits.append(row[i])
        meta.append(" | ".join(bits))
    symbols = [norm_symbol(i) for i in data_ids]
    if "CLDN4" not in symbols and not any(symbols):
        return {
            "names": names[:width],
            "meta": meta,
            "ids": data_ids,
            "matrix": matrix,
            "platform": platform.replace("GPL", "GPL") if platform else "",
            "expr": None,
        }
    expr = collapse_rows(symbols, matrix) if any(symbols) else None
    if expr and "CLDN4" in expr:
        return {
            "names": names[:width],
            "meta": meta,
            "platform": platform,
            "expr_raw": expr,
            "ids": None,
            "matrix": None,
        }
    return {
        "names": names[:width],
        "meta": meta,
        "ids": data_ids,
        "matrix": matrix,
        "platform": platform,
        "expr_raw": None,
    }


def gpl_target_map(gpl):
    gpl = gpl if str(gpl).upper().startswith("GPL") else "GPL" + str(gpl)
    cache = CACHE / "gpl" / f"{gpl}.json"
    with GPL_LOCK:
        if cache.exists():
            return json.loads(cache.read_text())
        return _gpl_target_map_unlocked(gpl, cache)


def _gpl_target_map_unlocked(gpl, cache):
    cache.parent.mkdir(parents=True, exist_ok=True)
    url = f"{FTP}/platforms/{geo_stub(gpl)}/{gpl}/annot/{gpl}.annot.gz"
    try:
        raw = fetch(url, timeout=120, retries=3)
    except Exception:
        cache.write_text("{}")
        return {}
    mapping = {}
    header = None
    sym_col = None
    with gzip.GzipFile(fileobj=io.BytesIO(raw)) as handle:
        for line in handle:
            text = line.decode("utf-8", "replace").rstrip("\n")
            if text.startswith("!") or text.startswith("^") or text.startswith("#") or not text:
                continue
            parts = text.split("\t")
            if header is None:
                header = [p.strip().lower() for p in parts]
                for i, col in enumerate(header):
                    if col in {"gene symbol", "gene_symbol", "symbol", "gene_assignment", "ilmn_gene"}:
                        sym_col = i
                        break
                if sym_col is None:
                    cache.write_text("{}")
                    return {}
                continue
            if sym_col >= len(parts):
                continue
            canon = norm_symbol(parts[sym_col])
            if canon:
                probe = parts[0].strip().strip('"')
                # one probe -> one gene; last write wins only if same
                prev = mapping.get(probe)
                if prev is None:
                    mapping[probe] = canon
                elif prev != canon:
                    mapping[probe] = ""  # ambiguous
    mapping = {k: v for k, v in mapping.items() if v}
    cache.write_text(json.dumps(mapping))
    return mapping


def apply_probe_map(ids, matrix, probe_map):
    symbols = []
    for probe in ids:
        symbols.append(probe_map.get(probe) or norm_symbol(probe))
    if not any(symbols):
        return None
    return collapse_rows(symbols, matrix)


def _symbolish(values):
    def ok(value):
        text = str(value).strip()
        if re.fullmatch(r"\d+(\.0)?", text):
            return True
        if text.upper().startswith("ENS"):
            return True
        return bool(re.match(r"^[A-Za-z][A-Za-z0-9_.-]*$", text))

    sample = list(values)[:40]
    if not sample:
        return 0.0
    return sum(ok(v) for v in sample) / len(sample)


def gene_column(df):
    lowered = {str(c).strip().lower(): c for c in df.columns}
    for key in (
        "gene_short_name",
        "gene_name",
        "gene_symbol",
        "hgnc_symbol",
        "external_gene_name",
        "marker.symbol",
        "marker_symbol",
        "symbol",
        "geneid",
        "gene_id",
        "entrezid",
        "entrez_id",
        "ensembl_gene_id",
        "ensembl",
        "ensg",
        "gene",
    ):
        if key in lowered:
            return lowered[key]
    first = df.columns[0]
    if str(first).startswith("Unnamed") or str(first).strip() == "":
        return first
    if _symbolish(df[first].tolist()) > 0.6:
        return first
    if "name" in lowered and _symbolish(df[lowered["name"]].tolist()) > 0.6:
        return lowered["name"]
    return None


def numeric_sample_columns(df, gene_col):
    drop_name = re.compile(
        r"\b(locus|refseq|tracking|tss|chromosome|chrom|start|end|strand|length|"
        r"gc|entrez|gene_id|description|biotype|uniprot|marker\.name)\b",
        re.I,
    )
    stat_name = re.compile(
        r"(log2?fc|logfc|logcpm|pvalue|p[\._-]?value|padj|adj[\._-]?p|fdr|qvalue|q[\._-]?value)",
        re.I,
    )
    numeric = []
    for col in df.columns:
        if col == gene_col:
            continue
        if drop_name.search(str(col)):
            continue
        series = pd.to_numeric(df[col], errors="coerce")
        if series.notna().mean() > 0.8:
            numeric.append(col)
    stat = [col for col in numeric if stat_name.search(str(col))]
    samples = [col for col in numeric if col not in stat]
    norm = [col for col in samples if re.search(r"(cpm|tpm|fpkm|rpkm|norm)", str(col), re.I)]
    raw = [col for col in samples if re.search(r"(rawcount|count)", str(col), re.I)]
    if len(norm) >= 2 and len(raw) >= 2:
        samples = [col for col in samples if col not in raw]
    if len(samples) >= 2:
        return samples
    logfc = [col for col in stat if re.search(r"log2?fc|logfc", str(col), re.I)]
    if len(logfc) == 1:
        return logfc
    return samples


def table_to_expr(df):
    gcol = gene_column(df)
    if gcol is None:
        return None
    scols = numeric_sample_columns(df, gcol)
    if len(scols) < 1:
        return None
    symbols = [norm_symbol(v) for v in df[gcol].tolist()]
    matrix = df[scols].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    expr = collapse_rows(symbols, matrix)
    if "CLDN4" not in expr:
        return None
    names = [str(c) for c in scols]
    return names, expr


def _peek_two(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", errors="replace") as handle:
        return handle.readline(), handle.readline()


def _read_shifted_header(path, sep):
    """First column is a gene id and the header row omitted its name."""
    opener = gzip.open if str(path).endswith(".gz") else open
    rows = []
    with opener(path, "rt", errors="replace") as handle:
        for line in handle:
            rows.append(line.rstrip("\n").split(sep))
    if len(rows) < 2:
        return None
    width = max(len(row) for row in rows[:50])
    header = rows[0]
    if len(header) != width - 1:
        return None
    header = ["gene_id"] + header
    body = [row + [""] * (width - len(row)) for row in rows[1:]]
    return pd.DataFrame([row[:width] for row in body], columns=header[:width])


def read_delimited(path):
    compression = "gzip" if str(path).endswith(".gz") else None
    line1, line2 = _peek_two(path)
    sep_guess = "\t" if line1.count("\t") >= line1.count(",") else ","
    if len(line2.rstrip("\n").split(sep_guess)) == len(line1.rstrip("\n").split(sep_guess)) + 1:
        df = _read_shifted_header(path, sep_guess)
        if df is not None and df.shape[1] >= 3:
            got = table_to_expr(df)
            if got:
                return got
    for sep in ("\t", ",", None):
        try:
            df = pd.read_csv(
                path,
                sep=sep,
                compression=compression,
                low_memory=False,
                dtype=str,
            )
        except Exception:
            continue
        if df.shape[1] < 2:
            continue
        got = table_to_expr(df)
        if got:
            return got
    return None


def read_xlsx(path):
    xl = pd.ExcelFile(path)
    for sheet in xl.sheet_names:
        df = xl.parse(sheet, dtype=str)
        got = table_to_expr(df)
        if got:
            return got
    return None


def prepare_expr(expr_raw):
    genes = list(expr_raw)
    matrix = np.vstack([expr_raw[g] for g in genes])
    if needs_log(matrix):
        matrix = np.log2(matrix + 1.0)
        scale = "log2(x+1)"
    else:
        scale = "as_deposited"
    ratio = looks_like_log_ratio(matrix)
    expr = {g: matrix[i] for i, g in enumerate(genes)}
    return expr, scale, ("log_ratio" if ratio else "single_channel")


def score_loaded(names, expr_raw, file_name):
    expr, scale, design = prepare_expr(expr_raw)
    row = score_expr(expr, names, design)
    design = row.get("design", design)
    row["scale"] = scale
    row["file_used"] = file_name
    row["joint_paired"] = joint_from_signs(row.get("delta_nhej"), row.get("delta_ifn"), flat=0.10)
    # |rho| <= 0.20 is too small to call that limb up or down.
    row["joint_spearman"] = joint_from_signs(row.get("rho_nhej"), row.get("rho_ifn"), flat=0.20)
    # Census uses Spearman once there are 6 samples. Paired and short-ratio
    # calls need both limbs to clear 0.10 log2 so a near-zero mean is not a hit.
    if row.get("spearman_n", 0) >= MIN_N and row.get("coverage_ok") and design != "log_ratio":
        row["joint_call"] = row["joint_spearman"]
        row["call_basis"] = "spearman"
    elif design == "log_ratio" and len(names) >= MIN_N and row.get("coverage_ok"):
        row["joint_call"] = row["joint_spearman"]
        row["call_basis"] = "spearman"
        row["contrast_note"] = "multi-sample log ratios; Spearman across samples, not one collapsed contrast"
    elif row.get("n_pairs"):
        row["joint_call"] = row["joint_paired"]
        row["call_basis"] = "paired_delta"
    elif design == "log_ratio" and row.get("joint_paired"):
        row["joint_call"] = row["joint_paired"]
        row["call_basis"] = "log_ratio"
    else:
        row["joint_call"] = ""
        row["call_basis"] = ""
    return row


def download_to(url, dest):
    dest.parent.mkdir(parents=True, exist_ok=True)
    raw = fetch(url, timeout=180, retries=3)
    dest.write_bytes(raw)
    return dest


def eligible_suppl(files):
    kept = []
    for item in files:
        name = item["name"]
        size = item["size"] or 0
        if size < MIN_BYTES or size > MAX_BYTES:
            continue
        if SKIP_FILE.search(name):
            continue
        lower = name.lower()
        if not re.search(r"\.(txt|tsv|csv|xlsx|xls)(\.gz)?$|\.zip$", lower):
            continue
        if not EXPR_NAME.search(name) and "series_matrix" not in lower:
            # still allow a single obvious table later; mark weak
            item = dict(item)
            item["weak"] = True
        else:
            item = dict(item)
            item["weak"] = False
        kept.append(item)
    def rank(item):
        name = item["name"].lower()
        if "vili" in name:
            base = 8
        elif "cldn" in name:
            base = 0
        elif "tpm" in name:
            base = 1
        elif "fpkm" in name:
            base = 2
        elif "norm" in name or "expression" in name:
            base = 3
        elif "count" in name:
            base = 4
        else:
            base = 5
        return (item["weak"], base, item["size"] or 0)

    kept.sort(key=rank)
    # drop weak files if a strong one exists
    strong = [i for i in kept if not i["weak"]]
    return strong or kept


def load_any_table(path):
    name = path.name.lower()
    if name.endswith(".xlsx") or name.endswith(".xls"):
        return read_xlsx(path)
    if name.endswith(".zip"):
        with zipfile.ZipFile(path) as zf:
            members = [n for n in zf.namelist() if re.search(r"\.(txt|tsv|csv)(\.gz)?$", n, re.I)]
            for member in members:
                dest = path.parent / Path(member).name
                dest.write_bytes(zf.read(member))
                got = read_delimited(dest)
                dest.unlink(missing_ok=True)
                if got:
                    return got
        return None
    return read_delimited(path)


def score_gse(rec):
    acc = rec["accession"]
    row = {
        "accession": acc,
        "queries": rec.get("queries", ""),
        "title": rec.get("title", ""),
        "taxon": rec.get("taxon", ""),
        "gdstype": rec.get("gdstype", ""),
        "n_samples_geo": rec.get("n_samples", ""),
        "gpl": rec.get("gpl", ""),
        "pubmed": rec.get("pubmed", ""),
    }
    try:
        stub = geo_stub(acc)
        base = f"{FTP}/series/{stub}/{acc}/"
        matrix_files = list_ftp(base + "matrix/")
        matrices = [
            f
            for f in matrix_files
            if f["name"].endswith("_series_matrix.txt.gz")
            and f["size"]
            and MIN_BYTES <= f["size"] <= MAX_BYTES
        ]
        loaded = None
        file_used = ""
        for item in matrices:
            dest = CACHE / acc / item["name"]
            download_to(item["url"], dest)
            parsed = parse_series_matrix(dest)
            dest.unlink(missing_ok=True)
            if not parsed:
                continue
            if parsed.get("expr_raw") and "CLDN4" in parsed["expr_raw"]:
                loaded = (parsed["names"], parsed["expr_raw"])
                file_used = item["name"]
                break
            if parsed.get("ids") is not None and parsed.get("platform"):
                probe_map = gpl_target_map(parsed["platform"])
                expr = apply_probe_map(parsed["ids"], parsed["matrix"], probe_map)
                if expr and "CLDN4" in expr:
                    loaded = (parsed["names"], expr)
                    file_used = item["name"]
                    break
        if loaded is None:
            suppl = eligible_suppl(list_ftp(base + "suppl/"))
            if not suppl and not matrices:
                # record why
                too_big = [
                    f["name"]
                    for f in matrix_files + list_ftp(base + "suppl/")
                    if (f["size"] or 0) > MAX_BYTES
                ]
                if too_big:
                    row["status"] = "skipped"
                    row["skip_reason"] = "over_25MB:" + ",".join(too_big[:4])
                elif any("RAW.tar" in f["name"] for f in list_ftp(base + "suppl/")):
                    row["status"] = "skipped"
                    row["skip_reason"] = "only_raw_tar_or_no_gene_table"
                else:
                    row["status"] = "skipped"
                    row["skip_reason"] = "no_open_matrix_under_25MB"
                return row
            for item in suppl:
                dest = CACHE / acc / item["name"]
                download_to(item["url"], dest)
                try:
                    got = load_any_table(dest)
                finally:
                    dest.unlink(missing_ok=True)
                if got:
                    loaded = got
                    file_used = item["name"]
                    break
        if loaded is None:
            row["status"] = "skipped"
            row["skip_reason"] = "open_file_has_no_CLDN4_or_unparsed"
            return row
        names, expr_raw = loaded
        scored = score_loaded(names, expr_raw, file_used)
        row.update(scored)
        if not row.get("cldn4_present"):
            row["status"] = "no_cldn4"
        elif not row.get("coverage_ok"):
            row["status"] = "low_coverage"
        else:
            row["status"] = "scored"
        return row
    except Exception as exc:  # noqa: BLE001
        row["status"] = "error"
        row["skip_reason"] = f"{type(exc).__name__}: {exc}"
        with WRITE_LOCK:
            with (TABLES / "errors.log").open("a") as handle:
                handle.write(f"\n== {acc} ==\n")
                handle.write(traceback.format_exc())
        return row


def is_expression(rec):
    gdstype = (rec.get("gdstype") or "").lower()
    title = rec.get("title") or ""
    if "non-coding rna" in gdstype:
        return False
    if re.search(r"\b(mirna|microrna|small rna)\b", title, re.I) and "mrna" not in title.lower():
        return False
    if "expression profiling" in gdstype:
        return True
    if any(token in gdstype for token in ("genome binding", "methylation", "genome variation", "snp")):
        return False
    if re.search(r"(rna-?seq|transcriptom|expression profil|mrna)", title, re.I):
        return True
    if not gdstype and not SKIP_TITLE.search(title):
        return True
    return False


def catalog():
    TABLES.mkdir(parents=True, exist_ok=True)
    membership = defaultdict(set)
    counts = []
    for name, term in QUERIES:
        print(f"search {name}", flush=True)
        if name in COUNT_ONLY and name != "dna_repair_title":
            count = esearch_count(term)
            counts.append({"query": name, "term": term, "count": count, "ids_fetched": 0})
            print(f"  count {count}", flush=True)
            continue
        count, ids = esearch_ids(term)
        counts.append({"query": name, "term": term, "count": count, "ids_fetched": len(ids)})
        print(f"  count {count} ids {len(ids)}", flush=True)
        for uid in ids:
            membership[uid].add(name)
    pd.DataFrame(counts).to_csv(TABLES / "query_counts.tsv", sep="\t", index=False)

    uids = list(membership)
    print(f"esummary {len(uids)}", flush=True)
    records = esummary(uids)
    rows = []
    for rec in records:
        if rec.get("entrytype") != "GSE":
            continue
        acc = rec.get("accession")
        if not acc:
            continue
        pubmed = rec.get("pubmedids") or []
        if isinstance(pubmed, list):
            pubmed = ",".join(str(x) for x in pubmed)
        gpl = rec.get("gpl")
        gpl = "" if gpl is None else str(gpl)
        if gpl and not gpl.upper().startswith("GPL"):
            gpl = "GPL" + gpl.split(";")[0]
        queries = sorted(membership.get(str(rec.get("uid")), []))
        rows.append(
            {
                "accession": acc,
                "uid": rec.get("uid"),
                "queries": ",".join(queries),
                "title": (rec.get("title") or "").replace("\t", " "),
                "taxon": rec.get("taxon") or "",
                "gdstype": rec.get("gdstype") or "",
                "n_samples": rec.get("n_samples") or "",
                "gpl": gpl,
                "pubmed": pubmed,
                "pdat": rec.get("pdat") or "",
            }
        )
    cat = pd.DataFrame(rows).drop_duplicates("accession")
    cat.to_csv(TABLES / "series_catalog.tsv", sep="\t", index=False)
    print(f"catalog {len(cat)}", flush=True)
    return cat


def priority(queries):
    q = set(queries.split(",")) if isinstance(queries, str) else set()
    if "cldn4_text" in q or "brute_token" in q:
        return 0
    if q & {"nhej", "prkdc", "sting1", "cgas_title"}:
        return 1
    return 2


def run_scores(cat, limit=None):
    TABLES.mkdir(parents=True, exist_ok=True)
    jsonl = TABLES / "scores.jsonl"
    done = set()
    if jsonl.exists():
        with jsonl.open() as handle:
            for line in handle:
                if line.strip():
                    done.add(json.loads(line)["accession"])
    queue = []
    for rec in cat.to_dict(orient="records"):
        if rec["accession"] in done:
            continue
        queries = set(str(rec["queries"]).split(","))
        if not (queries & DOWNLOAD_QUERIES):
            continue
        if not is_expression(rec):
            skipped = {
                "accession": rec["accession"],
                "queries": rec["queries"],
                "title": rec["title"],
                "taxon": rec["taxon"],
                "gdstype": rec["gdstype"],
                "n_samples_geo": rec["n_samples"],
                "status": "not_expression",
                "skip_reason": rec["gdstype"] or "title_filter",
            }
            with WRITE_LOCK:
                with jsonl.open("a") as handle:
                    handle.write(json.dumps(skipped) + "\n")
            continue
        try:
            n = int(rec["n_samples"])
        except (TypeError, ValueError):
            n = 0
        if n > 1500:
            skipped = {
                "accession": rec["accession"],
                "queries": rec["queries"],
                "title": rec["title"],
                "status": "skipped",
                "skip_reason": "n_samples_over_1500",
                "n_samples_geo": n,
            }
            with WRITE_LOCK:
                with jsonl.open("a") as handle:
                    handle.write(json.dumps(skipped) + "\n")
            continue
        queue.append(rec)
    queue.sort(key=lambda r: (priority(r["queries"]), r["accession"]))
    if limit:
        queue = queue[:limit]
    print(f"download queue {len(queue)} already {len(done)}", flush=True)

    from concurrent.futures import ThreadPoolExecutor, as_completed

    def _one(rec):
        result = score_gse(rec)
        with WRITE_LOCK:
            with jsonl.open("a") as handle:
                handle.write(json.dumps(result, default=str) + "\n")
        status = result.get("status")
        call = result.get("joint_call") or result.get("skip_reason") or ""
        print(f"{result['accession']} {status} {call}", flush=True)
        return result

    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = [pool.submit(_one, rec) for rec in queue]
        for fut in as_completed(futures):
            fut.result()
    rows = []
    with jsonl.open() as handle:
        for line in handle:
            if line.strip():
                rows.append(json.loads(line))
    inv = pd.DataFrame(rows)
    inv.to_csv(TABLES / "inventory.tsv", sep="\t", index=False)
    return inv


def add_fdr(inv):
    scored = inv.copy()
    mask = (
        (scored.get("status") == "scored")
        & (scored.get("call_basis") == "spearman")
        & scored["p_nhej"].apply(lambda v: isinstance(v, (int, float)) and math.isfinite(v))
    ) if "p_nhej" in scored.columns else pd.Series(False, index=scored.index)
    if "p_nhej" not in scored.columns:
        return scored
    from statsmodels.stats.multitest import multipletests

    scored["fdr_nhej"] = np.nan
    scored["fdr_ifn"] = np.nan
    scored["fdr_composite"] = np.nan
    idx = scored.index[mask.fillna(False)]
    if len(idx) >= 2:
        for src, dest in (("p_nhej", "fdr_nhej"), ("p_ifn", "fdr_ifn"), ("p_composite", "fdr_composite")):
            pvals = scored.loc[idx, src].astype(float).to_numpy()
            if np.isfinite(pvals).all():
                scored.loc[idx, dest] = multipletests(pvals, method="fdr_bh")[1]
    return scored


def summarize(inv):
    inv = add_fdr(inv)
    inv.to_csv(TABLES / "inventory.tsv", sep="\t", index=False)
    spearman = inv[inv.get("call_basis") == "spearman"] if "call_basis" in inv.columns else inv.iloc[0:0]
    paired = inv[inv.get("call_basis").isin(["paired_delta", "log_ratio"])] if "call_basis" in inv.columns else inv.iloc[0:0]

    def counts(frame, col):
        if frame.empty or col not in frame.columns:
            return {}
        return {str(k): int(v) for k, v in frame[col].fillna("").value_counts().to_dict().items()}

    status_counts = counts(inv, "status")
    summary = {
        "n_inventory_rows": int(len(inv)),
        "status_counts": status_counts,
        "spearman_n_series": int(len(spearman)),
        "spearman_joint": counts(spearman, "joint_call"),
        "paired_n_series": int(len(paired)),
        "paired_joint": counts(paired, "joint_call"),
        "spearman_fdr05_nhej_neg": int(
            (
                (spearman["fdr_nhej"] < 0.05) & (spearman["rho_nhej"] < 0)
            ).sum()
        )
        if len(spearman) and "fdr_nhej" in spearman.columns
        else 0,
        "spearman_fdr05_ifn_pos": int(
            (
                (spearman["fdr_ifn"] < 0.05) & (spearman["rho_ifn"] > 0)
            ).sum()
        )
        if len(spearman) and "fdr_ifn" in spearman.columns
        else 0,
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2), flush=True)
    write_figure(spearman)
    return summary


def write_figure(spearman):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    FIGS.mkdir(parents=True, exist_ok=True)
    if spearman.empty or "rho_nhej" not in spearman.columns:
        return
    df = spearman.dropna(subset=["rho_nhej", "rho_ifn"]).copy()
    if df.empty:
        return
    queries = df["queries"].fillna("") if "queries" in df.columns else ""
    focus = df[queries.str.contains("cldn4_text|brute_token")].copy()
    focus = focus.sort_values("rho_nhej")
    if focus.empty:
        focus = df.sort_values("rho_nhej").head(20)

    def color_for(call):
        if call == "NHEJ_down_IFN_up":
            return "#0072B2"
        if call == "NHEJ_up_IFN_down":
            return "#D55E00"
        return "#7F7F7F"

    colors = [color_for(c) for c in focus["joint_call"]]
    y = np.arange(len(focus))
    fig, axes = plt.subplots(
        1, 2, figsize=(9.4, max(3.8, 0.38 * len(focus) + 1.4)), sharey=True,
        gridspec_kw={"wspace": 0.08},
    )
    axes[0].scatter(focus["rho_nhej"], y, c=colors, s=36, zorder=3)
    axes[1].scatter(focus["rho_ifn"], y, c=colors, s=36, zorder=3)
    for ax, title in (
        (axes[0], "CLDN4 vs NHEJ"),
        (axes[1], "CLDN4 vs IFN"),
    ):
        ax.axvline(0, color="#444444", lw=0.8)
        ax.axvline(-0.2, color="#bbbbbb", lw=0.6, ls="--")
        ax.axvline(0.2, color="#bbbbbb", lw=0.6, ls="--")
        ax.set_xlabel("Spearman rho")
        ax.set_title(title)
        ax.set_xlim(-1.05, 1.05)
        ax.grid(axis="x", color="#eeeeee")
    axes[0].set_yticks(y)
    axes[0].set_yticklabels(focus["accession"].tolist(), fontsize=8)
    fig.suptitle(
        "CLDN4-text and brute-token series with n ≥ 6. Blue: NHEJ down and IFN up (|rho| > 0.20)",
        fontsize=10,
    )
    fig.tight_layout()
    fig.savefig(FIGS / "fig_spearman_nhej_ifn.png", dpi=140)
    fig.savefig(FIGS / "fig_spearman_nhej_ifn.pdf")
    plt.close(fig)


def self_test():
    url = f"{FTP}/series/GSE207nnn/GSE207704/suppl/GSE207704_CLDN4_RNAseq.txt.gz"
    dest = CACHE / "self" / "GSE207704_CLDN4_RNAseq.txt.gz"
    if not dest.exists():
        download_to(url, dest)
    names, expr_raw = read_delimited(dest)
    expr, scale, design = prepare_expr(expr_raw)
    assert design == "single_channel", design
    assert scale == "log2(x+1)", scale
    assert "PAXX" in expr and "CLDN4" in expr
    pairs = find_pairs(names)
    assert len(pairs) == 2, pairs
    d_n, _, _, _ = paired_gene_deltas(expr, NHEJ, pairs)
    d_i, _, _, _ = paired_gene_deltas(expr, IFN, pairs)
    d_c, _, _, _ = paired_gene_deltas(expr, ["CLDN4"], pairs)
    print("self-test", scale, design, "pairs", len(pairs), "CLDN4", d_c, "NHEJ", d_n, "IFN", d_i)
    assert abs(d_c - ((0.7059409658068487 + 1.0079322006889209) / 2)) < 1e-6
    assert abs(d_n - (-0.05125600234633107)) < 1e-6
    assert abs(d_i - 0.550370460757852) < 1e-6
    row = score_loaded(names, expr_raw, dest.name)
    assert row["joint_call"] == "mixed", row["joint_call"]
    assert row["n_nhej"] == 8
    assert row["n_ifn"] == 10
    print("self-test ok", row["joint_call"], "IFN genes", row["n_ifn"])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--scores-only", action="store_true")
    parser.add_argument("--summarize-only", action="store_true")
    args = parser.parse_args()
    TABLES.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(parents=True, exist_ok=True)
    if args.self_test:
        self_test()
        return
    if args.summarize_only:
        inv = pd.read_csv(TABLES / "inventory.tsv", sep="\t")
        summarize(inv)
        return
    if args.scores_only and (TABLES / "series_catalog.tsv").exists():
        cat = pd.read_csv(TABLES / "series_catalog.tsv", sep="\t", dtype=str).fillna("")
    else:
        cat = catalog()
    inv = run_scores(cat, limit=args.limit or None)
    summarize(inv)


if __name__ == "__main__":
    main()
