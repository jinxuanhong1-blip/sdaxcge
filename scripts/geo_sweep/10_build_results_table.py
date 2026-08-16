#!/usr/bin/env python3
"""
Build the single results table: results/geo_sweep/marker_stats.tsv

One row per (accession, gene) for every catalogued extra series × {TACSTD2, CLDN4}.
Additional rows are added only for GSE50927's four published DE contrasts (the
usable Cldn4-KO set). Numeric fields are left blank when a value was not
computed from a real table. No accessions or statistics are invented.
"""
import gzip
import io
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

NOTES = Path("notes/geo_sweep")
RESULTS = Path("results/geo_sweep")
SUPP = RESULTS / "supp"

GENES = ("TACSTD2", "CLDN4")

# columns that are annotation / DE stats, not sample expression
META_RE = re.compile(
    r"^(index|gene|symbol|entrez|ensembl|id|length|len|chr|chrom|"
    r"start|end|stop|strand|locus|description|genename|genetype|product|"
    r"dbxref|go_|go$|pathway|ko_entry|ec$|transcript|seqname|width|"
    r"basemean|log2foldchange|logfc|logcpm|lfcse|stat|pvalue|pval|padj|"
    r"fdr|qvalue|is_sig|comparison|gene_id|gene_name|geneid|sampleid|"
    r"marker\.|wts_psw|row\.names|supplementary)",
    re.IGNORECASE,
)

# Only emit a pairwise test when both group labels look biological.
# Blocks sample-ID prefixes (S105691) and assay-type splits (FPKM vs count).
BIO_GROUP_RE = re.compile(
    r"ko|kd|wt|wild|shgfp|shnik|shrna|ici|igg|pd-?1|pd-?l1|veh|dec|"
    r"dmso|control|ctrl|bom|tie|nik|vav|cldn|m1|m2|resistant|sensitive|"
    r"responder|treat",
    re.IGNORECASE,
)

DE_KEYS = [
    ("log2FoldChange", "de_logfc"),
    ("logFC", "de_logfc"),
    ("log2FC", "de_logfc"),
    ("pvalue", "de_pvalue"),
    ("PValue", "de_pvalue"),
    ("P.Value", "de_pvalue"),
    ("padj", "de_padj"),
    ("FDR", "de_padj"),
    ("adj.P.Val", "de_padj"),
]


def load_catalog():
    df = pd.read_csv(NOTES / "catalog.tsv", sep="\t", dtype=str).fillna("")
    return {r["accession"]: r for _, r in df.iterrows()}


def sample_values(row):
    """Numeric sample-like fields from a marker row dict."""
    out = {}
    for k, v in row.items():
        if v is None or (isinstance(v, float) and np.isnan(v)):
            continue
        key = str(k).strip()
        if key.lower().startswith("unnamed"):
            pass  # keep unnamed numeric sample columns
        elif META_RE.search(key):
            continue
        if isinstance(v, (int, float, np.integer, np.floating)):
            out[str(k)] = float(v)
    return out


def strip_rep(name):
    t = re.sub(r"[_.\- ]+(rep(licate)?|r)?[_.\- ]*\d+$", "", name, flags=re.I)
    t = re.sub(r"[_.\- ]+\d+$", "", t)
    return t.strip(" _-.") or name


def group_stats(values):
    """Group sample names by stripped replicate suffix. Return list of dicts."""
    buckets = {}
    for k, v in values.items():
        buckets.setdefault(strip_rep(k), []).append(v)
    rows = []
    for g, vals in buckets.items():
        arr = np.array(vals, dtype=float)
        rows.append({
            "group": g, "n": int(arr.size),
            "mean": float(np.mean(arr)),
            "median": float(np.median(arr)),
            "sd": float(np.std(arr, ddof=1)) if arr.size > 1 else None,
        })
    return rows


def pairwise(gst):
    """If exactly two groups with n>=2, Mann-Whitney + mean difference."""
    usable = [g for g in gst if g["n"] >= 2]
    if len(usable) != 2:
        return {}
    # rebuild arrays from means? we need raw values — caller passes values
    return None  # filled by pairwise_from_values


def pairwise_from_values(values, require_bio=True):
    buckets = {}
    for k, v in values.items():
        buckets.setdefault(strip_rep(k), []).append(v)
    items = [(g, np.array(v, dtype=float)) for g, v in buckets.items() if len(v) >= 2]
    if len(items) != 2:
        return {}
    (ga, va), (gb, vb) = items
    if require_bio and not (BIO_GROUP_RE.search(ga) and BIO_GROUP_RE.search(gb)):
        return {}
    try:
        _, p = stats.mannwhitneyu(va, vb, alternative="two-sided")
        p = float(p)
    except Exception:
        p = None
    return {
        "group_a": ga, "n_a": int(va.size), "mean_a": float(np.mean(va)),
        "group_b": gb, "n_b": int(vb.size), "mean_b": float(np.mean(vb)),
        "mean_diff_b_minus_a": float(np.mean(vb) - np.mean(va)),
        "mannwhitney_p": p,
    }


def collapse_sh_groups(values):
    """Map LM05_shGFP_* / LM05_shNIK_* onto two KD groups (unique keys)."""
    out = {}
    n_g = n_n = 0
    for k, v in values.items():
        kl = k.lower()
        if "shgfp" in kl:
            n_g += 1
            out[f"shGFP_{n_g}"] = v
        elif "shnik" in kl or "sh_nik" in kl:
            n_n += 1
            out[f"shNIK_{n_n}"] = v
        else:
            out[k] = v
    return out


def de_fields(row):
    out = {}
    lower = {str(k).lower(): (k, row[k]) for k in row}
    mapping = {
        "de_logfc": ["log2foldchange", "logfc", "log2fc"],
        "de_pvalue": ["pvalue", "p.value", "p_value"],
        "de_padj": ["padj", "fdr", "adj.p.val", "qvalue"],
    }
    for dest, cands in mapping.items():
        for c in cands:
            if c in lower:
                k, v = lower[c]
                if isinstance(v, (int, float, np.integer, np.floating)) and not pd.isna(v):
                    out[dest] = float(v)
                break
    return out


def read_gse50927_contrasts():
    files = [
        ("Cldn4lungWTvsKO", "GSE50927_Cldn4lungWTvsKOgenes.csv.gz"),
        ("VILI_WT", "GSE50927_VILIwtGenes.csv.gz"),
        ("VILI_WT_vs_KO_high", "GSE50927_VILIwtkohiGenes.csv.gz"),
        ("VILI_WT_vs_KO_low", "GSE50927_VILIwtkoloGenes.csv.gz"),
    ]
    out = []
    for contrast, name in files:
        path = SUPP / name
        if not path.exists():
            continue
        raw = gzip.decompress(path.read_bytes()).decode("utf-8", "replace")
        df = pd.read_csv(io.StringIO(raw))
        symcol = next((c for c in df.columns
                       if str(c).lower() in ("marker.symbol", "genesymbol", "symbol")),
                      None)
        if symcol is None:
            continue
        for gene, sym in (("TACSTD2", "Tacstd2"), ("CLDN4", "Cldn4")):
            hit = df[df[symcol].astype(str).str.lower() == sym.lower()]
            if hit.empty:
                continue
            r = hit.iloc[0]
            out.append({
                "accession": "GSE50927",
                "gene": gene,
                "contrast": contrast,
                "source_file": name,
                "de_logfc": float(r["logFC"]),
                "de_logcpm": float(r["logCPM"]),
                "de_pvalue": float(r["PValue"]),
                "de_padj": float(r["FDR"]),
            })
    return out


def blank_row(acc, gene, cat):
    return {
        "accession": acc,
        "gene": gene,
        "is_ici": cat.get("is_ici", ""),
        "is_kd_ko": cat.get("is_kd_ko", ""),
        "organism": cat.get("organism", ""),
        "category": cat.get("category", ""),
        "gds_type": cat.get("gds_type", ""),
        "n_geo_samples": cat.get("n_samples", ""),
        "title": cat.get("title", ""),
        "marker_status": "",
        "source": "",
        "source_file": "",
        "table_kind": "",
        "n_values": "",
        "value_mean": "",
        "value_median": "",
        "value_min": "",
        "value_max": "",
        "n_nonzero": "",
        "de_logfc": "",
        "de_pvalue": "",
        "de_padj": "",
        "contrast": "",
        "group_a": "",
        "n_a": "",
        "mean_a": "",
        "group_b": "",
        "n_b": "",
        "mean_b": "",
        "mean_diff_b_minus_a": "",
        "mannwhitney_p": "",
        "groups_json": "",
        "notes": "",
    }


def fill_numeric(row, **kw):
    for k, v in kw.items():
        if v is None:
            continue
        if isinstance(v, float):
            row[k] = f"{v:.6g}"
        else:
            row[k] = v


def main():
    catalog = load_catalog()
    embedded = {r["accession"]: r for r in json.loads((NOTES / "analysis_embedded.json").read_text())}
    supp = {r["accession"]: r for r in json.loads((NOTES / "analysis_supp.json").read_text())}
    inspect_titles = {}
    for ir in json.loads((NOTES / "matrix_inspect.json").read_text()):
        if ir.get("accession") and ir.get("sample_titles") and ir["accession"] not in inspect_titles:
            inspect_titles[ir["accession"]] = ir["sample_titles"]

    rows = []
    for acc, cat in catalog.items():
        for gene in GENES:
            rec = blank_row(acc, gene, cat)
            notes = []

            emb = embedded.get(acc, {})
            eg = (emb.get("genes") or {}).get(gene)
            if eg and eg.get("status") == "ok" and eg.get("per_sample"):
                raw_vals = {k: float(v) for k, v in eg["per_sample"].items() if v is not None}
                # prefer author sample titles as group keys when lengths match
                titles = emb.get("group_labels") or emb.get("sample_titles") or []
                if titles and len(titles) == len(raw_vals):
                    vals = {f"{t}__{i}": v for i, (v, t) in
                            enumerate(zip(raw_vals.values(), titles))}
                else:
                    vals = raw_vals
                arr = np.array(list(raw_vals.values()), dtype=float)
                rec["marker_status"] = "found"
                rec["source"] = "embedded_matrix"
                rec["source_file"] = f"{acc}_series_matrix.txt.gz"
                rec["table_kind"] = "expression"
                fill_numeric(rec,
                             n_values=int(arr.size),
                             value_mean=float(np.mean(arr)),
                             value_median=float(np.median(arr)),
                             value_min=float(np.min(arr)),
                             value_max=float(np.max(arr)),
                             n_nonzero=int(np.sum(arr != 0)))
                gst = group_stats(vals)
                rec["groups_json"] = json.dumps(
                    [{**g, "mean": round(g["mean"], 6),
                      "median": round(g["median"], 6),
                      "sd": None if g["sd"] is None else round(g["sd"], 6)}
                     for g in gst], separators=(",", ":"))
                pw = pairwise_from_values(vals, require_bio=False)
                if pw:
                    fill_numeric(rec, **pw)
                if emb.get("scale_guess"):
                    notes.append(f"scale={emb['scale_guess']}")
            elif eg and eg.get("status") == "gene_not_on_platform":
                rec["marker_status"] = "not_on_platform"
                rec["source"] = "embedded_matrix"
                rec["source_file"] = f"{acc}_series_matrix.txt.gz"
                notes.append("gene absent from embedded platform table")

            if rec["marker_status"] in ("", "not_on_platform"):
                sr = supp.get(acc, {})
                sm = (sr.get("markers") or {}).get(gene)
                if sm:
                    rec["marker_status"] = "found"
                    rec["source"] = "supplementary"
                    rec["source_file"] = sm.get("file", "")
                    row0 = (sm.get("rows") or [{}])[0]
                    if sm.get("is_de_table"):
                        rec["table_kind"] = "DE"
                        fill_numeric(rec, **de_fields(row0))
                        notes.append("DE table (author-computed contrast)")
                    else:
                        rec["table_kind"] = "scRNA" if sm.get("streamed") else "expression"
                        vals = sample_values(row0)
                        # GSE328294 mixes FPKM.* and count.* — keep FPKM only
                        fpkm = {k: v for k, v in vals.items()
                                if str(k).lower().startswith("fpkm")}
                        counts = {k: v for k, v in vals.items()
                                  if str(k).lower().startswith("count")}
                        if fpkm and counts:
                            vals = fpkm
                            notes.append("FPKM columns only (raw count columns ignored)")
                        if acc == "GSE182261":
                            vals = collapse_sh_groups(vals)
                        # when GEO sample titles line up 1:1 with numeric columns, use them
                        titles = inspect_titles.get(acc) or []
                        if titles and len(titles) == len(vals) and acc != "GSE182261":
                            vals = {f"{title}__{i}": v for i, (title, v) in
                                    enumerate(zip(titles, vals.values()))}
                        if vals:
                            arr = np.array(list(vals.values()), dtype=float)
                            fill_numeric(rec,
                                         n_values=int(arr.size),
                                         value_mean=float(np.mean(arr)),
                                         value_median=float(np.median(arr)),
                                         value_min=float(np.min(arr)),
                                         value_max=float(np.max(arr)),
                                         n_nonzero=int(np.sum(arr != 0)))
                            gst = group_stats(vals)
                            if 1 < len(gst) <= 12:
                                rec["groups_json"] = json.dumps(
                                    [{**g, "mean": round(g["mean"], 6),
                                      "median": round(g["median"], 6),
                                      "sd": None if g["sd"] is None else round(g["sd"], 6)}
                                     for g in gst], separators=(",", ":"))
                            if cat.get("is_ici") == "Y" or cat.get("is_kd_ko") == "Y":
                                pw = pairwise_from_values(vals, require_bio=True)
                                if pw:
                                    fill_numeric(rec, **pw)
                    if sm.get("streamed"):
                        notes.append("stream-extracted from >300MB supplementary matrix")
                elif rec["marker_status"] == "":
                    status = sr.get("status", "")
                    if status == "no_candidate_file":
                        rec["marker_status"] = "no_processed_table"
                        notes.append("no gene-level processed supplementary table listed")
                    elif status == "markers_absent":
                        rec["marker_status"] = "not_in_table"
                        tried = [f.get("name") for f in sr.get("files_tried", []) if f.get("name")]
                        if any(f.get("status") == "too_large" for f in sr.get("files_tried", [])):
                            rec["marker_status"] = "too_large_or_absent"
                            notes.append("candidate file too large or gene row not present")
                        else:
                            notes.append("processed table present; gene row not found")
                        if tried:
                            rec["source_file"] = tried[0]
                            rec["source"] = "supplementary"
                    elif acc in embedded and (embedded[acc].get("status") in
                                              ("analyzed", "no_expression_table")
                                              or embedded[acc].get("genes")):
                        rec["marker_status"] = rec["marker_status"] or "not_in_table"
                    else:
                        rec["marker_status"] = "no_processed_table"
                        notes.append("no usable processed expression for this gene")

            rec["notes"] = "; ".join(notes)
            rows.append(rec)

    # GSE50927: replace the single DE row pair with all four published contrasts
    rows = [r for r in rows if not (r["accession"] == "GSE50927" and r["marker_status"] == "found")]
    cat = catalog["GSE50927"]
    for extra in read_gse50927_contrasts():
        rec = blank_row("GSE50927", extra["gene"], cat)
        rec["marker_status"] = "found"
        rec["source"] = "supplementary"
        rec["source_file"] = extra["source_file"]
        rec["table_kind"] = "DE"
        rec["contrast"] = extra["contrast"]
        fill_numeric(rec,
                     de_logfc=extra["de_logfc"],
                     de_pvalue=extra["de_pvalue"],
                     de_padj=extra["de_padj"])
        rec["notes"] = (f"edgeR DE contrast {extra['contrast']}; "
                        f"logCPM={extra['de_logcpm']:.4g}")
        rows.append(rec)

    cols = list(blank_row("", "", {}).keys())
    out = RESULTS / "marker_stats.tsv"
    with open(out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    n = len(rows)
    n_found = sum(1 for r in rows if r["marker_status"] == "found")
    n_ici_found = sum(1 for r in rows if r["is_ici"] == "Y" and r["marker_status"] == "found")
    n_kd_found = sum(1 for r in rows if r["is_kd_ko"] == "Y" and r["marker_status"] == "found")
    n_cmp = sum(1 for r in rows if r["mannwhitney_p"] != "")
    n_de = sum(1 for r in rows if r["de_logfc"] != "")
    print(f"Wrote {out}")
    print(f"  rows={n}  found={n_found}  ICI-found={n_ici_found}  "
          f"KD-found={n_kd_found}  DE={n_de}  pairwise={n_cmp}")


if __name__ == "__main__":
    main()
