#!/usr/bin/env python3
"""IFN / NHEJ / STING scores inside Milo neighbourhoods.

Counts are read from the public GEO matrices that built the Harmony
object. A missing official symbol is taken from one HGNC previous
symbol when that alias is in the matrix. Genes still absent contribute
0. The score is the mean of log1p(count / nCount_RNA * 1e4) over the
fixed gene-set list.

Within each locked unit, delta = mean(score in the focal neighbourhood
compartment) - mean(score in other same-class neighbourhoods). The
signed-rank is across units. This is not a second Milo model and not
a spatial test.
"""

from __future__ import annotations

import csv
import gzip
import json
import subprocess
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.stats import wilcoxon
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from deepen_milo import CLASS_NAMES, build_graph, class_group_from_counts  # noqa: E402
from run_milo import class_fractions, count_by_sample, spatial_fdr_kdistance  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "results" / "cache"
DEEP = CACHE / "deepen"
TABLES = ROOT / "results" / "tables"
FIGS = ROOT / "figures"
GEO = Path("/tmp/geo_c4")
SETS_PATH = ROOT / "data" / "gene_sets.json"

# HGNC previous symbols. Used only when the official symbol is absent
# and this alias is present. One alias per gene.
ALIASES = {
    "STING1": "TMEM173",
    "CGAS": "MB21D1",
    "MRE11": "MRE11A",
    "MARCHF1": "MARCH1",
    "MARCHF5": "MARCH5",
    "WARS1": "WARS",
    "TENT5A": "FAM46A",
    "NSD2": "WHSC1",
    "CYREN": "C7orf49",
    "MRNIP": "C5orf45",
    "PAXX": "C9orf142",
    "SHLD1": "C20orf196",
    "SHLD2": "FAM35A",
}

SCORE_SETS = ["IFN", "NHEJ", "STING", "CGAS_STING_GO"]
SFDR = 0.05


def load_sets() -> dict[str, list[str]]:
    raw = json.loads(SETS_PATH.read_text())
    return {name: list(raw[name]["genes"]) for name in SCORE_SETS}


def gene_columns(universe_upper: dict[str, int], genes: list[str]) -> tuple[list[int], list[str]]:
    """Map each official gene to a column index in an upper-cased universe.

    Returns (index or -1, source label official/alias/absent).
    """
    idx = []
    src = []
    for g in genes:
        key = g.upper()
        if key in universe_upper:
            idx.append(universe_upper[key])
            src.append("official")
            continue
        alias = ALIASES.get(g)
        if alias and alias.upper() in universe_upper:
            idx.append(universe_upper[alias.upper()])
            src.append(f"alias:{alias}")
            continue
        idx.append(-1)
        src.append("absent")
    return idx, src


def first_upper_index(names: list[str]) -> dict[str, int]:
    out: dict[str, int] = {}
    for i, name in enumerate(names):
        key = name.upper()
        if key not in out:
            out[key] = i
    return out


def suffix_barcode(cell_id: str, dataset: str, unit_id: str) -> str:
    prefix = f"{dataset}_{unit_id}_"
    if not cell_id.startswith(prefix):
        raise ValueError(f"cell id {cell_id} does not start with {prefix}")
    return cell_id[len(prefix) :]


def extract_gse123(meta: pd.DataFrame, genes: list[str], counts: np.ndarray, row_of: dict[str, int]) -> list[dict]:
    sub = meta[meta["dataset"] == "GSE123902"]
    want = {}
    for cid, unit in zip(sub.index, sub["unit_id"].astype(str)):
        want.setdefault(unit, {})[suffix_barcode(cid, "GSE123902", unit)] = cid
    per_gene = [Counter() for _ in genes]
    n_donors = 0
    n_matched = 0
    # spot-check library size on the first 20 matched cells
    check_ids = list(sub.index[:20])
    check_sum = {c: 0.0 for c in check_ids}
    files = sorted(GEO.joinpath("gse123").glob("*_dense.csv.gz"))
    files = [p for p in files if "NORMAL" not in p.name]
    for path in files:
        # GSM*_MSK_{donor}_{tissue}_dense.csv.gz
        parts = path.name.split("_")
        # MSK donor is parts[2]
        donor = parts[2]
        if donor not in want:
            continue
        with gzip.open(path, "rt", newline="") as fh:
            reader = csv.reader(fh)
            header = next(reader)
            universe = first_upper_index(header[1:])
            col_idx, src = gene_columns(universe, genes)
            # Each dense file drops genes that are all zero in that donor.
            n_donors += 1
            for j, label in enumerate(src):
                per_gene[j][label] += 1
            wanted = want[donor]
            for row in reader:
                bc = row[0]
                cid = wanted.get(bc)
                if cid is None:
                    continue
                n_matched += 1
                r = row_of[cid]
                vals = row[1:]
                for j, c in enumerate(col_idx):
                    if c >= 0:
                        counts[r, j] = int(float(vals[c]))
                if cid in check_sum:
                    check_sum[cid] = float(sum(int(float(v)) for v in vals))
    print(f"GSE123902 matched {n_matched} / {len(sub)}", flush=True)
    if n_matched != len(sub):
        raise SystemExit(f"GSE123902 barcode match {n_matched} != {len(sub)}")
    # coverage recorded once from the first file; confirm later files share symbols
    # by re-checking only that every requested alias/official decision is stable.
    diffs = []
    for cid, total in check_sum.items():
        diffs.append(abs(total - float(meta.loc[cid, "nCount_RNA"])))
    print(
        f"GSE123902 nCount vs matrix sum: max abs {max(diffs):.4g} on {len(diffs)} cells",
        flush=True,
    )
    src = []
    for j in range(len(genes)):
        parts = [f"{k}:{v}/{n_donors}" for k, v in sorted(per_gene[j].items())]
        src.append(";".join(parts))
    return [{"dataset": "GSE123902", "genes": genes, "src": src}]


def extract_gse131(meta: pd.DataFrame, genes: list[str], counts: np.ndarray, row_of: dict[str, int]) -> list[dict]:
    sub = meta[meta["dataset"] == "GSE131907"]
    suffix_to_cid = {}
    for cid, unit in zip(sub.index, sub["unit_id"].astype(str)):
        suffix_to_cid[suffix_barcode(cid, "GSE131907", unit)] = cid
    path = GEO / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        name_to_col = {name: i for i, name in enumerate(header)}
        missing = [s for s in suffix_to_cid if s not in name_to_col]
        print(f"GSE131907 header cols {len(header)-1} unmatched cells {len(missing)}", flush=True)
        if missing:
            raise SystemExit(f"GSE131907 unmatched example {missing[:3]}")
        # column index in the data line (0 is the gene name, header[0] is Index)
        use_cols = []
        use_rows = []
        for suffix, cid in suffix_to_cid.items():
            use_cols.append(name_to_col[suffix])
            use_rows.append(row_of[cid])
        use_cols = np.array(use_cols, dtype=np.int32)
        use_rows = np.array(use_rows, dtype=np.int32)
        # gene universe: first pass is too expensive; resolve symbols on the fly
        # Build the set of matrix symbols we will accept.
        accept = {}
        for j, g in enumerate(genes):
            accept[g.upper()] = (j, "official")
            alias = ALIASES.get(g)
            if alias:
                accept.setdefault(alias.upper(), (j, f"alias:{alias}"))
        found_src = ["absent"] * len(genes)
        n_gene_rows = 0
        for line in fh:
            tab = line.find("\t")
            symbol = line[:tab]
            key = symbol.upper()
            hit = accept.get(key)
            if hit is None:
                continue
            j, src = hit
            if found_src[j] == "official" and src.startswith("alias"):
                continue
            if found_src[j] != "absent" and found_src[j] != src and found_src[j] == "official":
                continue
            # first official wins; alias only if still absent
            if found_src[j] == "absent" or (src == "official" and found_src[j].startswith("alias")):
                fields = line.rstrip("\n").split("\t")
                vals = np.fromiter((int(float(fields[c])) for c in use_cols), dtype=np.int32, count=len(use_cols))
                counts[use_rows, j] = vals
                found_src[j] = src
            n_gene_rows += 1
    print(
        f"GSE131907 matched {len(sub)} cells; genes official {(np.array(found_src)=='official').sum()} "
        f"alias {sum(s.startswith('alias') for s in found_src)} absent {found_src.count('absent')}",
        flush=True,
    )
    return [{"dataset": "GSE131907", "genes": genes, "src": found_src}]


def extract_gse189(meta: pd.DataFrame, genes: list[str], counts: np.ndarray, row_of: dict[str, int]) -> list[dict]:
    sub = meta[meta["dataset"] == "GSE189357"]
    want = {}
    for cid, unit in zip(sub.index, sub["unit_id"].astype(str)):
        bc = suffix_barcode(cid, "GSE189357", unit).replace(".", "-")
        want.setdefault(unit, {})[bc] = cid
    dest = GEO / "gse189"
    if not any(dest.glob("*_matrix.mtx.gz")):
        dest.mkdir(parents=True, exist_ok=True)
        subprocess.check_call(["tar", "-C", str(dest), "-xf", str(GEO / "GSE189357_RAW.tar")])
    per_gene = [Counter() for _ in genes]
    n_samples_seen = 0
    n_matched = 0
    check = {}
    for unit, bc_map in want.items():
        mtx = next(dest.glob(f"*_{unit}_matrix.mtx.gz"))
        feat = next(dest.glob(f"*_{unit}_features.tsv.gz"))
        bc_path = next(dest.glob(f"*_{unit}_barcodes.tsv.gz"))
        with gzip.open(feat, "rt") as fh:
            symbols = [line.rstrip("\n").split("\t")[1] for line in fh]
        universe = first_upper_index(symbols)
        col_idx, src = gene_columns(universe, genes)
        n_samples_seen += 1
        for j, label in enumerate(src):
            per_gene[j][label] += 1
        with gzip.open(bc_path, "rt") as fh:
            barcodes = [line.strip() for line in fh]
        keep_pos = []
        keep_rows = []
        for i, bc in enumerate(barcodes):
            cid = bc_map.get(bc)
            if cid is not None:
                keep_pos.append(i)
                keep_rows.append(row_of[cid])
        mat = mmread(gzip.open(mtx, "rb")).tocsr()
        if mat.shape[1] != len(barcodes) or mat.shape[0] != len(symbols):
            raise SystemExit(f"GSE189357 {unit} mtx shape {mat.shape} vs {len(symbols)} x {len(barcodes)}")
        gene_rows = [c for c in col_idx if c >= 0]
        gene_js = [j for j, c in enumerate(col_idx) if c >= 0]
        submat = mat[gene_rows][:, keep_pos].toarray()
        for local, j in enumerate(gene_js):
            counts[keep_rows, j] = submat[local].astype(np.int32)
        n_matched += len(keep_rows)
        # library-size check on up to 5 cells
        for pos, cid in list(zip(keep_pos, [meta.index[r] for r in keep_rows]))[:5]:
            check[cid] = float(mat[:, pos].sum())
        print(f"  GSE189357 {unit} matched {len(keep_rows)}", flush=True)
    print(f"GSE189357 matched {n_matched} / {len(sub)}", flush=True)
    if n_matched != len(sub):
        raise SystemExit(f"GSE189357 barcode match {n_matched} != {len(sub)}")
    diffs = [abs(check[c] - float(meta.loc[c, "nCount_RNA"])) for c in check]
    print(f"GSE189357 nCount vs matrix sum: max abs {max(diffs):.4g} on {len(diffs)} cells", flush=True)
    src = [";".join(f"{k}:{v}/{n_samples_seen}" for k, v in sorted(per_gene[j].items())) for j in range(len(genes))]
    return [{"dataset": "GSE189357", "genes": genes, "src": src}]


def extract_gse205(meta: pd.DataFrame, genes: list[str], counts: np.ndarray, row_of: dict[str, int]) -> list[dict]:
    sub = meta[meta["dataset"] == "GSE205335"]
    manifest = DEEP / "gse205_cells.tsv"
    DEEP.mkdir(parents=True, exist_ok=True)
    rows = []
    for cid, unit in zip(sub.index, sub["unit_id"].astype(str)):
        rows.append({"cell_id": cid, "barcode": suffix_barcode(cid, "GSE205335", unit), "nCount_RNA": meta.loc[cid, "nCount_RNA"]})
    pd.DataFrame(rows).to_csv(manifest, sep="\t", index=False)
    gene_file = DEEP / "gse205_genes.txt"
    gene_file.write_text("\n".join(genes) + "\n")
    alias_file = DEEP / "gse205_aliases.tsv"
    pd.DataFrame([{"official": k, "alias": v} for k, v in ALIASES.items()]).to_csv(alias_file, sep="\t", index=False)
    out_counts = DEEP / "gse205_counts.tsv.gz"
    out_src = DEEP / "gse205_src.tsv"
    rds = GEO / "work" / "gse205.rds"
    if not rds.exists():
        rds.parent.mkdir(parents=True, exist_ok=True)
        subprocess.check_call(
            f"gunzip -c {GEO / 'GSE205335_Lung_IO_UMI_matrix.rds.gz'} | gunzip > {rds}",
            shell=True,
        )
    r_code = f"""
    .libPaths(c(path.expand("~/R/library"), .libPaths()))
    suppressPackageStartupMessages(library(Matrix))
    mat <- readRDS({json.dumps(str(rds))})
    rn <- toupper(rownames(mat))
    keep <- !duplicated(rn)
    mat <- mat[keep, , drop = FALSE]
    rownames(mat) <- rn[keep]
    cells <- read.delim({json.dumps(str(manifest))}, stringsAsFactors = FALSE)
    genes <- readLines({json.dumps(str(gene_file))})
    genes <- genes[nzchar(genes)]
    alias <- read.delim({json.dumps(str(alias_file))}, stringsAsFactors = FALSE)
    alias_map <- setNames(alias$alias, alias$official)
    src <- rep("absent", length(genes))
    take <- integer(length(genes))
    for (i in seq_along(genes)) {{
      g <- toupper(genes[[i]])
      if (g %in% rownames(mat)) {{
        src[[i]] <- "official"
        take[[i]] <- match(g, rownames(mat))
      }} else {{
        ai <- match(genes[[i]], names(alias_map))
        if (!is.na(ai)) {{
          a <- toupper(alias_map[[ai]])
          if (a %in% rownames(mat)) {{
            src[[i]] <- paste0("alias:", alias_map[[ai]])
            take[[i]] <- match(a, rownames(mat))
          }}
        }}
      }}
    }}
    present <- which(take > 0)
    bc <- cells$barcode
    hit <- bc %in% colnames(mat)
    cat("matched", sum(hit), "of", length(bc), "\\n")
    if (sum(hit) != length(bc)) {{
      alt <- gsub("_", "-", colnames(mat))
      # second chance: embedding barcode uses the matrix name already
      stop("GSE205335 unmatched ", sum(!hit), " example ", bc[which(!hit)[1]])
    }}
    sub <- mat[take[present], bc, drop = FALSE]
    # write genes x cells
    out <- as.matrix(sub)
    storage.mode(out) <- "integer"
    rownames(out) <- genes[present]
    gz <- gzfile({json.dumps(str(out_counts))}, "w")
    write.table(out, file = gz, sep = "\\t", quote = FALSE, col.names = NA)
    close(gz)
    write.table(data.frame(gene = genes, src = src), {json.dumps(str(out_src))},
                sep = "\\t", quote = FALSE, row.names = FALSE)
    # library size check
    chk <- sample(bc, 20)
    lib <- Matrix::colSums(mat[, chk, drop = FALSE])
    ref <- setNames(cells$nCount_RNA, cells$barcode)
    cat("nCount max abs", max(abs(as.numeric(lib) - as.numeric(ref[chk]))), "\\n")
    """
    subprocess.check_call(["Rscript", "-e", r_code])
    src_df = pd.read_csv(out_src, sep="\t")
    mat = pd.read_csv(out_counts, sep="\t", index_col=0)
    # columns are barcodes
    for cid, bc in zip(sub.index, [suffix_barcode(c, "GSE205335", u) for c, u in zip(sub.index, sub["unit_id"].astype(str))]):
        if bc not in mat.columns:
            raise SystemExit(f"missing column {bc}")
        r = row_of[cid]
        for j, g in enumerate(genes):
            if g in mat.index:
                counts[r, j] = int(mat.at[g, bc])
    print(f"GSE205335 filled {len(sub)} cells", flush=True)
    return [{"dataset": "GSE205335", "genes": genes, "src": src_df["src"].tolist()}]


def coverage_rows(blocks: list[dict]) -> pd.DataFrame:
    rows = []
    for block in blocks:
        for g, src in zip(block["genes"], block["src"]):
            rows.append({"dataset": block["dataset"], "gene": g, "source": src})
    return pd.DataFrame(rows)


def score_matrix(counts: np.ndarray, ncount: np.ndarray, gene_index: dict[str, int], gene_list: list[str]) -> np.ndarray:
    idx = [gene_index[g] for g in gene_list]
    raw = counts[:, idx].astype(np.float64)
    lib = ncount.astype(np.float64)
    cp = np.log1p(raw / lib[:, None] * 1e4)
    return cp.mean(axis=1)


def verify_membership(members: list[np.ndarray], sub: pd.DataFrame, prefix: str) -> None:
    """The rebuilt graph must reproduce the cached neighbourhood counts."""
    cached = pd.read_csv(DEEP / f"counts_{prefix}.tsv", sep="\t", index_col=0)
    expected_ids = [f"{prefix}_nh{i}" for i in range(len(members))]
    if list(cached.index) != expected_ids:
        raise SystemExit(f"{prefix} nhood ids differ from the cached count matrix")
    sample_index = {s: i for i, s in enumerate(cached.columns)}
    codes = sub["unit_key"].map(sample_index)
    if codes.isna().any():
        raise SystemExit(f"{prefix} has cells with no sample column")
    got = count_by_sample(members, codes.to_numpy(dtype=np.int32), cached.shape[1])
    if not np.array_equal(got, cached.to_numpy()):
        raise SystemExit(f"{prefix} membership does not match the cached counts")
    print(f"  membership matches cached counts {cached.shape}", flush=True)


def assign_from_status(members: list[np.ndarray], status: np.ndarray, n_cells: int):
    """status per nhood: 2=up hit, 1=down hit, 0=other target, -1=ignore."""
    up = np.zeros(n_cells, dtype=bool)
    down = np.zeros(n_cells, dtype=bool)
    other = np.zeros(n_cells, dtype=bool)
    for i, cells in enumerate(members):
        st = int(status[i])
        if st < 0:
            continue
        if st == 2:
            up[cells] = True
        elif st == 1:
            down[cells] = True
        else:
            other[cells] = True
    both = int((up & down).sum())
    # priority up > down > other
    down = down & ~up
    other = other & ~up & ~down
    return up, down, other, both


def nhood_status(nhood_ids: list[str], groups: np.ndarray, da: pd.DataFrame, target: str) -> np.ndarray:
    da = da.set_index("nhood_id")
    # class labels in the DA table must match the rebuilt graph
    common = [i for i, nh in enumerate(nhood_ids) if nh in da.index]
    rebuilt = groups[common]
    saved = da.loc[[nhood_ids[i] for i in common], "class_group"].to_numpy()
    if not np.array_equal(rebuilt, saved):
        n_bad = int((rebuilt != saved).sum())
        raise SystemExit(f"rebuilt class_group mismatches DA on {n_bad} / {len(common)} nhoods")
    sfdr = da["SpatialFDR"]
    logfc = da["logFC"]
    status = np.full(len(nhood_ids), -1, dtype=np.int8)
    for i, nh in enumerate(nhood_ids):
        if groups[i] != target:
            continue
        if nh not in da.index:
            status[i] = 0  # untested target-class nhood counts as other
            continue
        if sfdr.loc[nh] < SFDR and logfc.loc[nh] > 0:
            status[i] = 2
        elif sfdr.loc[nh] < SFDR and logfc.loc[nh] < 0:
            status[i] = 1
        else:
            status[i] = 0
    return status


def paired_test(meta: pd.DataFrame, score: np.ndarray, focal: np.ndarray, background: np.ndarray, cell_ok: np.ndarray):
    use = cell_ok & (focal | background)
    rows = []
    deltas = []
    units = meta.loc[use, "unit_key"].to_numpy()
    datasets = meta.loc[use, "dataset"].to_numpy()
    # map unit -> dataset
    unit_ds = {}
    for u, d in zip(meta["unit_key"], meta["dataset"]):
        unit_ds[u] = d
    focal_u = meta.loc[cell_ok & focal, "unit_key"]
    back_u = meta.loc[cell_ok & background, "unit_key"]
    sc_f = pd.Series(score[cell_ok & focal], index=focal_u.index)
    sc_b = pd.Series(score[cell_ok & background], index=back_u.index)
    # groupby unit
    f_mean = sc_f.groupby(focal_u.to_numpy()).agg(["mean", "size"])
    b_mean = sc_b.groupby(back_u.to_numpy()).agg(["mean", "size"])
    both = f_mean.index.intersection(b_mean.index)
    for u in both:
        delta = float(f_mean.loc[u, "mean"] - b_mean.loc[u, "mean"])
        deltas.append(delta)
        rows.append(
            {
                "unit_key": u,
                "dataset": unit_ds[u],
                "n_focal": int(f_mean.loc[u, "size"]),
                "n_other": int(b_mean.loc[u, "size"]),
                "mean_focal": float(f_mean.loc[u, "mean"]),
                "mean_other": float(b_mean.loc[u, "mean"]),
                "delta": delta,
            }
        )
    arr = np.array(deltas, dtype=float)
    out = {
        "n_patients": int(len(arr)),
        "median_delta": float(np.median(arr)) if len(arr) else np.nan,
        "mean_delta": float(np.mean(arr)) if len(arr) else np.nan,
        "n_neg": int((arr < 0).sum()) if len(arr) else 0,
        "n_pos": int((arr > 0).sum()) if len(arr) else 0,
        "n_zero": int((arr == 0).sum()) if len(arr) else 0,
        "wilcoxon_stat": np.nan,
        "wilcoxon_p": np.nan,
    }
    nonzero = arr[arr != 0]
    if len(nonzero) >= 1 and len(arr) >= 5:
        res = wilcoxon(arr, zero_method="wilcox", alternative="two-sided", method="auto")
        out["wilcoxon_stat"] = float(res.statistic)
        out["wilcoxon_p"] = float(res.pvalue)
    return out, rows


def rebuild(meta: pd.DataFrame, X: np.ndarray, k: int, d: int, mask: np.ndarray, prefix: str):
    sub = meta.loc[mask]
    Xd = X[mask][:, :d]
    members, kdist = build_graph(Xd, k)
    class_map = {c: i for i, c in enumerate(CLASS_NAMES)}
    class_codes = sub["cell_class"].map(class_map).to_numpy(dtype=np.int32)
    cc = class_fractions(members, class_codes, len(CLASS_NAMES))
    groups = class_group_from_counts(cc)
    nhood_ids = [f"{prefix}_nh{i}" for i in range(len(members))]
    # members are indices into the masked cell list
    return sub, members, groups, nhood_ids, kdist


def global_to_local_mask(members, n_sub: int, n_all: int, mask: np.ndarray, local_bool: np.ndarray) -> np.ndarray:
    """Map a boolean over sub-cells to a boolean over all cells."""
    full = np.zeros(n_all, dtype=bool)
    full[np.flatnonzero(mask)] = local_bool
    return full


def main() -> None:
    DEEP.mkdir(parents=True, exist_ok=True)
    TABLES.mkdir(parents=True, exist_ok=True)
    sets = load_sets()
    # stable gene order: union in set order
    genes = []
    seen = set()
    for name in SCORE_SETS:
        for g in sets[name]:
            if g not in seen:
                genes.append(g)
                seen.add(g)
    gene_index = {g: i for i, g in enumerate(genes)}
    meta = pd.read_csv(CACHE / "meta.tsv", sep="\t", index_col=0)
    meta = meta.copy()
    meta["unit_key"] = meta["dataset"].astype(str) + "|" + meta["unit_id"].astype(str)
    counts_path = DEEP / "gene_counts.npy"
    genes_path = DEEP / "gene_counts.genes.txt"
    cov_path = TABLES / "gene_coverage.tsv"
    if counts_path.exists() and genes_path.exists() and genes_path.read_text().splitlines() == genes and cov_path.exists():
        counts = np.load(counts_path)
        print(f"loaded cached counts {counts.shape}", flush=True)
    else:
        counts = np.zeros((len(meta), len(genes)), dtype=np.int32)
        row_of = {cid: i for i, cid in enumerate(meta.index)}
        blocks = []
        blocks += extract_gse123(meta, genes, counts, row_of)
        blocks += extract_gse189(meta, genes, counts, row_of)
        blocks += extract_gse131(meta, genes, counts, row_of)
        blocks += extract_gse205(meta, genes, counts, row_of)
        np.save(counts_path, counts)
        genes_path.write_text("\n".join(genes) + "\n")
        cov = coverage_rows(blocks)
        # tag which set each gene belongs to (a gene may be in several)
        cov.to_csv(cov_path, sep="\t", index=False)
        print(cov.groupby(["dataset", "source"]).size() if False else cov.groupby("dataset")["source"].value_counts())

    ncount = meta["nCount_RNA"].to_numpy(dtype=np.float64)
    if np.any(ncount <= 0):
        raise SystemExit("nonpositive nCount_RNA")
    scores = {name: score_matrix(counts, ncount, gene_index, sets[name]) for name in SCORE_SETS}

    X = pd.read_csv(CACHE / "harmony.tsv.gz", sep="\t", index_col=0).loc[meta.index].to_numpy(dtype=np.float64)
    specs = [
        {
            "name": "primary_k30_d30_q34",
            "k": 30,
            "d": 30,
            "graph": "full",
            "prefix": "full_k30_d30",
            "da": TABLES / "da_grid_primary_k30_d30_q34.tsv",
            "mask": np.ones(len(meta), dtype=bool),
            "role": "primary",
        },
        {
            "name": "countmax_k15_d10_q4_vs_rest",
            "k": 15,
            "d": 10,
            "graph": "full",
            "prefix": "full_k15_d10",
            "da": TABLES / "da_max_tnk_sfdr05.tsv",
            "mask": np.ones(len(meta), dtype=bool),
            "role": "count_max_sensitivity",
        },
        {
            "name": "malignant_only_k30_d30_q34",
            "k": 30,
            "d": 30,
            "graph": "mal",
            "prefix": "mal_k30_d30",
            "da": DEEP / "edger" / "mal_k30_d30__q34.tsv",
            "mask": (meta["cell_class"] == "malignant").to_numpy(),
            "role": "malignant_intrinsic",
        },
    ]

    summary_rows = []
    delta_rows = []
    plot_payload = {}
    for spec in specs:
        print(f"graph {spec['name']}", flush=True)
        sub, members, groups, nhood_ids, kdist = rebuild(
            meta, X, spec["k"], spec["d"], spec["mask"], spec["prefix"]
        )
        verify_membership(members, sub, spec["prefix"])
        da = pd.read_csv(spec["da"], sep="\t")
        id_pos = {nh: i for i, nh in enumerate(nhood_ids)}
        if "SpatialFDR" not in da.columns:
            da = da.copy()
            da["SpatialFDR"] = spatial_fdr_kdistance(
                da["PValue"].to_numpy(),
                np.array([kdist[id_pos[i]] for i in da["nhood_id"]], dtype=float),
            )
            da["class_group"] = [groups[id_pos[i]] for i in da["nhood_id"]]
        else:
            recomputed = spatial_fdr_kdistance(
                da["PValue"].to_numpy(),
                np.array([kdist[id_pos[i]] for i in da["nhood_id"]], dtype=float),
            )
            gap = float(np.nanmax(np.abs(recomputed - da["SpatialFDR"].to_numpy())))
            print(f"  SpatialFDR rebuild max abs gap {gap:.3g}", flush=True)
            if gap > 1e-6:
                raise SystemExit(f"graph did not reproduce SpatialFDR for {spec['name']}")
        targets = ["malignant", "T/NK"] if spec["graph"] == "full" else ["malignant"]
        for target in targets:
            status = nhood_status(nhood_ids, groups, da, target)
            up_l, down_l, other_l, n_both = assign_from_status(members, status, len(sub))
            print(
                f"  {target} cells up {int(up_l.sum())} down {int(down_l.sum())} "
                f"other {int(other_l.sum())} in both up and down before priority {n_both}",
                flush=True,
            )
            # restrict to cells of that class
            if target == "T/NK":
                class_ok = np.isin(sub["cell_class"].to_numpy(), ["T", "NK"])
            else:
                class_ok = sub["cell_class"].to_numpy() == "malignant"
            up_l = up_l & class_ok
            down_l = down_l & class_ok
            other_l = other_l & class_ok
            comparisons = [("up_vs_other", up_l, other_l), ("down_vs_other", down_l, other_l)]
            if target == "T/NK":
                # thesis abundance claim is T/NK down; still report up_vs_other
                pass
            for set_name, score_all in scores.items():
                score_sub = score_all[spec["mask"]]
                for comp_name, focal, background in comparisons:
                    # malignant-intrinsic graph has no T/NK contrast; skip empty
                    if focal.sum() == 0 or background.sum() == 0:
                        summary_rows.append(
                            {
                                "spec": spec["name"],
                                "role": spec["role"],
                                "k": spec["k"],
                                "d": spec["d"],
                                "graph": spec["graph"],
                                "target": target,
                                "comparison": comp_name,
                                "gene_set": set_name,
                                "n_cells_focal": int(focal.sum()),
                                "n_cells_other": int(background.sum()),
                                "n_patients": 0,
                                "median_delta": np.nan,
                                "mean_delta": np.nan,
                                "n_neg": 0,
                                "n_pos": 0,
                                "n_zero": 0,
                                "wilcoxon_stat": np.nan,
                                "wilcoxon_p": np.nan,
                                "n_cells_in_both_hits": n_both,
                            }
                        )
                        continue
                    # paired_test expects arrays aligned to `meta` (all cells) OR to sub.
                    # Build a meta-like frame for sub cells.
                    stat, deltas = paired_test(sub, score_sub, focal, background, np.ones(len(sub), dtype=bool))
                    stat.update(
                        {
                            "spec": spec["name"],
                            "role": spec["role"],
                            "k": spec["k"],
                            "d": spec["d"],
                            "graph": spec["graph"],
                            "target": target,
                            "comparison": comp_name,
                            "gene_set": set_name,
                            "n_cells_focal": int(focal.sum()),
                            "n_cells_other": int(background.sum()),
                            "n_cells_in_both_hits": n_both,
                        }
                    )
                    summary_rows.append(stat)
                    for row in deltas:
                        row.update(
                            {
                                "spec": spec["name"],
                                "target": target,
                                "comparison": comp_name,
                                "gene_set": set_name,
                            }
                        )
                        delta_rows.append(row)
                    if (
                        spec["name"] == "primary_k30_d30_q34"
                        and set_name == "IFN"
                        and comp_name in ("up_vs_other", "down_vs_other")
                    ):
                        plot_payload[(target, comp_name)] = [r["delta"] for r in deltas]

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(TABLES / "nhood_gene_paired.tsv", sep="\t", index=False)
    pd.DataFrame(delta_rows).to_csv(TABLES / "nhood_gene_patient_deltas.tsv", sep="\t", index=False)
    print(summary[summary.gene_set.isin(["IFN", "NHEJ", "STING"])][
        ["spec", "target", "comparison", "gene_set", "n_patients", "median_delta", "n_neg", "n_pos", "wilcoxon_p"]
    ].to_string(index=False))

    # Figure: primary IFN deltas
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.6), sharey=False)
    panels = [
        ("malignant", "up_vs_other", "Malignant up − other\nIFN, primary k=30 d=30"),
        ("T/NK", "down_vs_other", "T/NK down − other\nIFN, primary k=30 d=30"),
    ]
    for ax, (target, comp, title) in zip(axes, panels):
        vals = np.sort(np.array(plot_payload.get((target, comp), []), dtype=float))
        if len(vals) == 0:
            ax.set_title(title)
            continue
        ax.axvline(0, color="0.5", lw=0.8)
        ax.scatter(vals, np.arange(len(vals)), s=12, c="#333333", zorder=3)
        ax.axvline(np.median(vals), color="#b22222", lw=1.2, label=f"median {np.median(vals):.3f}")
        ax.set_yticks([])
        ax.set_xlabel("Within-unit Δ log1p CP10k")
        ax.set_title(title, fontsize=10)
        ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGS / "nhood_ifn_delta_primary.png", dpi=160)
    fig.savefig(FIGS / "nhood_ifn_delta_primary.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
