#!/usr/bin/env python3
"""
Stage 04 - extract TACSTD2/CLDN4 and test vs ICI outcome for leftover series
that Stage 03 flagged as usable_maybe.

Honest rules
------------
* Do not re-analyse GSE207422 / GSE243238 (already done).
* Only leftover series with both genes *and* a per-sample ICI outcome.
* If no leftover series is usable, write an empty tests table and a
  documented skip list. That is a valid scientific result.

Outputs:
  results/w200/GEO_2023/analysis/marker_outcome_tests.csv
  results/w200/GEO_2023/analysis/<GSE>_per_sample.csv
  results/w200/GEO_2023/analysis/*.png  (only if a test ran)
"""
import csv
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from geo_common import (
    OUT, CACHE, GENE_ALIASES, series_dir, cached_get, open_maybe_gz,
    fetch_series_matrix,
)

ANA = OUT / "analysis"
ANA.mkdir(parents=True, exist_ok=True)

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAVE_MPL = True
except Exception:
    HAVE_MPL = False


def cliffs_delta(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) == 0 or len(b) == 0:
        return np.nan
    gt = sum(float(x > y) for x in a for y in b)
    lt = sum(float(x < y) for x in a for y in b)
    return (gt - lt) / (len(a) * len(b))


def parse_sample_meta(text):
    gsm = None
    fields = {}
    char_rows = []
    for line in text.splitlines():
        if line.startswith("!series_matrix_table_begin"):
            break
        if not line.startswith("!Sample_"):
            continue
        parts = line.rstrip("\n").split("\t")
        key = parts[0].lstrip("!")
        vals = [p.strip().strip('"') for p in parts[1:]]
        if key == "Sample_geo_accession":
            gsm = vals
        elif key == "Sample_title":
            fields["title"] = vals
        elif key == "Sample_source_name_ch1":
            fields["source_name"] = vals
        elif key.startswith("Sample_characteristics_ch"):
            char_rows.append(vals)
    if gsm is None:
        return pd.DataFrame()
    df = pd.DataFrame({"gsm": gsm})
    for k, v in fields.items():
        if len(v) == len(gsm):
            df[k] = v
    parsed = {}
    for row in char_rows:
        if len(row) != len(gsm):
            continue
        for i, cell in enumerate(row):
            if ":" in cell:
                k, val = cell.split(":", 1)
                parsed.setdefault(k.strip().lower(), [None] * len(gsm))
                parsed[k.strip().lower()][i] = val.strip()
    for k, col in parsed.items():
        df[k] = col
    return df


def load_expr_table(path):
    """Load a gene x sample or sample x gene table; return DataFrame."""
    # try pandas auto
    name = str(path).lower()
    if name.endswith(".xlsx") or name.endswith(".xls"):
        return pd.read_excel(path)
    sep = "\t"
    # sniff
    with open_maybe_gz(path) as fh:
        header = fh.readline()
    if header.count(",") > header.count("\t"):
        sep = ","
    df = pd.read_csv(path, sep=sep, index_col=0)
    return df


def extract_genes(expr):
    """Return {gene: Series indexed by sample} or empty."""
    # decide orientation: genes in index vs columns
    idx = [str(x).split(".")[0].upper() for x in expr.index]
    cols = [str(x).split(".")[0].upper() for x in expr.columns]
    alias_map = {a.upper(): g for g, als in GENE_ALIASES.items() for a in als}

    def pick(labels, axis):
        found = {}
        for i, lab in enumerate(labels):
            if lab in alias_map:
                g = alias_map[lab]
                if axis == 0:
                    found[g] = expr.iloc[i]
                else:
                    found[g] = expr.iloc[:, i]
        return found

    found = pick(idx, 0)
    if len(found) < 2:
        found2 = pick(cols, 1)
        if len(found2) > len(found):
            found = found2
    return found


def map_outcome(meta):
    """
    Try to build a binary responder / non-responder column.
    Returns (series_or_None, outcome_name, group1, group2, note).
    """
    cols = {c.lower(): c for c in meta.columns}
    # preferred keys
    candidates = []
    for key in cols:
        if re.search(r"(mpr|pcr|pathologic response|recist|response|"
                     r"clinical benefit|benefit|responder)", key):
            candidates.append(cols[key])
    if not candidates:
        return None, "", "", "", "no outcome column"

    col = candidates[0]
    raw = meta[col].astype(str)
    responder_pat = re.compile(
        r"\b(mpr|pcr|cpr|cr|pr|yes|responder|benefit|sensitive|"
        r"major pathologic|complete|partial response)\b", re.I)
    non_pat = re.compile(
        r"\b(nmpr|npr|sd|pd|no|non[- ]?responder|no benefit|resist|"
        r"progressive|stable)\b", re.I)
    labels = []
    for v in raw:
        vl = v.lower()
        if responder_pat.search(vl) and not non_pat.search(vl):
            labels.append("responder")
        elif non_pat.search(vl):
            labels.append("non_responder")
        else:
            labels.append("other")
    s = pd.Series(labels, index=meta.index)
    n1 = int((s == "responder").sum())
    n2 = int((s == "non_responder").sum())
    if n1 < 2 or n2 < 2:
        return s, col, "responder", "non_responder", \
            f"imbalanced/unusable groups n_resp={n1} n_non={n2}"
    return s, col, "responder", "non_responder", "ok"


def boxplot(groups, data, title, fname, ylabel):
    if not HAVE_MPL:
        return
    fig, ax = plt.subplots(figsize=(4.2, 4))
    vals = [data[g] for g in groups]
    bp = ax.boxplot(vals, tick_labels=groups, patch_artist=True, widths=0.55)
    colors = ["#4C72B0", "#C44E52"]
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c)
        patch.set_alpha(0.55)
    rng = np.random.default_rng(0)
    for i, g in enumerate(groups, 1):
        y = data[g]
        x = rng.normal(i, 0.05, len(y))
        ax.scatter(x, y, color="black", s=16, zorder=3, alpha=0.7)
    ax.set_title(title, fontsize=10)
    ax.set_ylabel(ylabel)
    fig.tight_layout()
    fig.savefig(ANA / fname, dpi=140)
    plt.close(fig)


def main():
    probe_path = OUT / "leftover_probe.csv"
    if not probe_path.exists():
        raise SystemExit("run 03_probe_leftovers.py first")
    probe = pd.read_csv(probe_path)
    usable = probe[probe["usable_maybe"] == True]  # noqa: E712
    tests = []
    skip = []

    if usable.empty:
        print("no leftover series flagged usable_maybe")
        pd.DataFrame(columns=[
            "cohort", "gene", "outcome", "group1", "n1", "median1",
            "group2", "n2", "median2", "mannwhitney_U", "p_value",
            "cliffs_delta", "note",
        ]).to_csv(ANA / "marker_outcome_tests.csv", index=False)
        (ANA / "SKIPPED.md").write_text(
            "No leftover 2023 lung ICI series passed usable_maybe "
            "(both genes present + per-sample ICI outcome + open processed "
            "expression). See leftover_probe.csv.\n")
        return

    for _, row in usable.iterrows():
        acc = row["accession"]
        print(f"analysing {acc}")
        text = fetch_series_matrix(acc)
        if not text:
            skip.append({"accession": acc, "reason": "no series_matrix"})
            continue
        meta = parse_sample_meta(text)
        if meta.empty:
            skip.append({"accession": acc, "reason": "empty sample meta"})
            continue

        expr = None
        peeked = str(row.get("peeked_file") or "")
        if peeked and peeked != "series_matrix":
            url = f"{series_dir(acc)}/suppl/{peeked}"
            p = cached_get(url, f"{acc}__{peeked}")
            if p:
                try:
                    expr = load_expr_table(p)
                except Exception as e:  # noqa: BLE001
                    skip.append({"accession": acc,
                                 "reason": f"expr parse failed: {e}"})
                    continue
        if expr is None and row.get("has_expression_table"):
            # parse table from series matrix
            rows = []
            in_table = False
            header = None
            for line in text.splitlines():
                if line.startswith("!series_matrix_table_begin"):
                    in_table = True
                    continue
                if line.startswith("!series_matrix_table_end"):
                    break
                if not in_table:
                    continue
                parts = [x.strip().strip('"') for x in line.split("\t")]
                if header is None:
                    header = parts
                    continue
                rows.append(parts)
            if header and rows:
                expr = pd.DataFrame(rows, columns=header).set_index(header[0])
                expr = expr.apply(pd.to_numeric, errors="coerce")

        if expr is None:
            skip.append({"accession": acc, "reason": "no expression table"})
            continue

        genes = extract_genes(expr)
        missing = [g for g in ("TACSTD2", "CLDN4") if g not in genes]
        if missing:
            skip.append({"accession": acc,
                         "reason": f"genes missing after load: {missing}"})
            continue

        # align samples
        # try GSM, then title, then any overlapping labels
        expr_cols = [str(c) for c in next(iter(genes.values())).index]
        meta = meta.copy()
        meta["_key"] = meta["gsm"].astype(str)
        overlap = set(meta["_key"]) & set(expr_cols)
        if len(overlap) < 4 and "title" in meta.columns:
            meta["_key"] = meta["title"].astype(str)
            overlap = set(meta["_key"]) & set(expr_cols)
        if len(overlap) < 4:
            skip.append({"accession": acc,
                         "reason": f"sample overlap <4 (n={len(overlap)})"})
            continue
        meta = meta[meta["_key"].isin(overlap)].copy()
        meta = meta.drop_duplicates("_key")
        meta = meta.set_index("_key")

        labels, oname, g1, g2, note = map_outcome(meta)
        if labels is None or note != "ok":
            skip.append({"accession": acc,
                         "reason": f"outcome unusable: {note} ({oname})"})
            # still write per-sample expression
            per = meta.copy()
            for g, s in genes.items():
                per[g] = pd.to_numeric(s.reindex(per.index), errors="coerce")
            per.to_csv(ANA / f"{acc}_per_sample.csv")
            continue

        per = meta.copy()
        per["outcome_group"] = labels.reindex(per.index)
        for g, s in genes.items():
            per[g] = pd.to_numeric(s.reindex(per.index), errors="coerce")
        per.to_csv(ANA / f"{acc}_per_sample.csv")

        for gene in ("TACSTD2", "CLDN4"):
            a = per.loc[per["outcome_group"] == g1, gene].dropna().astype(float)
            b = per.loc[per["outcome_group"] == g2, gene].dropna().astype(float)
            if len(a) < 2 or len(b) < 2:
                tests.append({
                    "cohort": acc, "gene": gene, "outcome": oname,
                    "group1": g1, "n1": len(a), "median1": "",
                    "group2": g2, "n2": len(b), "median2": "",
                    "mannwhitney_U": "", "p_value": "",
                    "cliffs_delta": "", "note": "n<2 in a group",
                })
                continue
            u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
            d = cliffs_delta(a, b)
            tests.append({
                "cohort": acc, "gene": gene, "outcome": oname,
                "group1": g1, "n1": int(len(a)),
                "median1": round(float(np.median(a)), 3),
                "group2": g2, "n2": int(len(b)),
                "median2": round(float(np.median(b)), 3),
                "mannwhitney_U": float(u),
                "p_value": round(float(p), 4),
                "cliffs_delta": round(float(d), 3),
                "note": "",
            })
            boxplot([g1, g2], {g1: a.values, g2: b.values},
                    f"{acc} {gene} vs {oname}",
                    f"{acc}_{gene}_outcome.png",
                    f"{gene} expression")

    pd.DataFrame(tests).to_csv(ANA / "marker_outcome_tests.csv", index=False)
    if skip:
        pd.DataFrame(skip).to_csv(ANA / "skipped_usable_maybe.csv", index=False)
    print(f"tests={len(tests)} skipped={len(skip)}")


if __name__ == "__main__":
    main()
