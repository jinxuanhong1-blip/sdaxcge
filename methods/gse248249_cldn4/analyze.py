#!/usr/bin/env python3
"""GSE248249: score CLDN4 and TACSTD2 on public Clariom D vs pre/post.

All numbers come from the GEO series matrix. No invented values.
"""
from __future__ import annotations

import gzip
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE248nnn/GSE248249/matrix/"
    "GSE248249_series_matrix.txt.gz"
)
MATRIX = Path("/tmp/gse248249/GSE248249_series_matrix.txt.gz")
OUT = Path(__file__).resolve().parent / "results"
OUT.mkdir(parents=True, exist_ok=True)


def ensure_matrix(path: Path) -> None:
    if path.exists() and path.stat().st_size > 0:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    import urllib.request

    urllib.request.urlretrieve(MATRIX_URL, path)

# Official-symbol transcript clusters on GPL23126 (Clariom D, gene version).
# Verified against GEO platform table gene_assignment (NM_002353 / NM_001305).
PROBES = {
    "TACSTD2": "TC0100014340.hg.1",
    "CLDN4": "TC0700007993.hg.1",
}
# Paper-highlighted immune genes used only as a processing sanity check.
# Official-symbol transcript clusters from GPL23126 gene_assignment (one each).
CONTROLS = {
    "CD8A": "TC0200013323.hg.1",
    "GZMA": "TC0500007413.hg.1",
    "CXCL9": "TC0400011052.hg.1",
    "CD274": "TC0900006559.hg.1",
    "IFNG": "TC1200011177.hg.1",
    "B2M": "TC1500007097.hg.1",
}


def parse_quoted(line: str) -> list[str]:
    return re.findall(r'"([^"]*)"', line)


def load_matrix(path: Path) -> tuple[pd.DataFrame, dict[str, list[str]], int, set[str]]:
    meta: dict[str, list[str]] = {}
    wanted = set(PROBES.values()) | set(CONTROLS.values())
    expr_rows: dict[str, list[str]] = {}
    samples: list[str] | None = None
    in_table = False
    n_features = 0
    found_ids: set[str] = set()

    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if not in_table:
                if line.startswith("!"):
                    key, _, rest = line[1:].partition("\t")
                    vals = parse_quoted(rest) if rest.strip() else parse_quoted(line)
                    if key in meta:
                        # characteristics fields repeat; keep last-seen lists
                        # by suffixing an index later if needed
                        i = 2
                        while f"{key}_{i}" in meta:
                            i += 1
                        meta[f"{key}_{i}"] = vals
                    else:
                        meta[key] = vals
                continue
            parts = line.rstrip("\n").split("\t")
            rid = parts[0].strip('"')
            if rid == "ID_REF":
                samples = [p.strip('"') for p in parts[1:]]
                continue
            n_features += 1
            if rid in wanted:
                expr_rows[rid] = parts[1:]
                found_ids.add(rid)

    if samples is None:
        raise RuntimeError("ID_REF row missing")
    if not expr_rows:
        raise RuntimeError("requested probes not found in series matrix")

    df = pd.DataFrame(
        {pid: np.array(vals, dtype=float) for pid, vals in expr_rows.items()},
        index=samples,
    )
    return df, meta, n_features, found_ids


def build_pheno(meta: dict[str, list[str]], samples: list[str]) -> pd.DataFrame:
    titles = meta["Sample_title"]
    acc = meta["Sample_geo_accession"]
    assert acc == samples

    char_rows = []
    for k, v in meta.items():
        if k.startswith("Sample_characteristics_ch1") and v and ":" in v[0]:
            field = v[0].split(":", 1)[0].strip()
            char_rows.append((field, [x.split(":", 1)[1].strip() if ":" in x else x for x in v]))

    pheno = pd.DataFrame({"geo_accession": acc, "title": titles})
    for field, vals in char_rows:
        pheno[field] = vals

    pheno["patient"] = pheno["title"].str.extract(r"Patient (\d+)", expand=False)
    pheno["patient"] = pheno["patient"].str.zfill(2)
    pheno["title_timepoint"] = np.where(
        pheno["title"].str.contains("pre-immunotherapy", case=False),
        "Pre-treatment",
        np.where(
            pheno["title"].str.contains("post-immunotherapy", case=False),
            "Post-treatment",
            "unknown",
        ),
    )
    if "timepoint" in pheno.columns:
        mismatch = pheno["timepoint"] != pheno["title_timepoint"]
        if mismatch.any():
            raise RuntimeError("title vs characteristics timepoint mismatch")
    return pheno


def summarize(x: np.ndarray) -> dict:
    x = np.asarray(x, dtype=float)
    return {
        "n": int(x.size),
        "mean": float(np.mean(x)),
        "median": float(np.median(x)),
        "min": float(np.min(x)),
        "max": float(np.max(x)),
    }


def cliffs_delta(a: np.ndarray, b: np.ndarray) -> float:
    """δ = P(a>b) - P(a<b); a=post, b=pre so positive = post higher."""
    a = np.asarray(a)
    b = np.asarray(b)
    gt = sum((ai > b).sum() for ai in a)
    lt = sum((ai < b).sum() for ai in a)
    n = a.size * b.size
    return float((gt - lt) / n)


def main() -> None:
    ensure_matrix(MATRIX)
    expr, meta, n_features, found_ids = load_matrix(MATRIX)
    pheno = build_pheno(meta, list(expr.index))

    missing_targets = [g for g, p in PROBES.items() if p not in found_ids]
    if "CLDN4" in missing_targets:
        raise SystemExit("CLDN4 absent from series matrix; stop.")

    gene_to_probe = {**PROBES, **{k: v for k, v in CONTROLS.items() if v in found_ids}}
    long_rows = []
    for gene, probe in gene_to_probe.items():
        tmp = pheno.copy()
        tmp["gene"] = gene
        tmp["probe_id"] = probe
        tmp["rma"] = expr[probe].to_numpy()
        long_rows.append(tmp)
    long = pd.concat(long_rows, ignore_index=True)

    pre_mask = pheno["timepoint"] == "Pre-treatment"
    post_mask = pheno["timepoint"] == "Post-treatment"
    n_pre = int(pre_mask.sum())
    n_post = int(post_mask.sum())
    n_patients = int(pheno["patient"].nunique())

    paired_patients = sorted(
        set(pheno.loc[pre_mask, "patient"]) & set(pheno.loc[post_mask, "patient"])
    )
    n_paired = len(paired_patients)

    # same-site pairs
    same_site = 0
    site_notes = []
    for p in paired_patients:
        pre_site = pheno.loc[(pheno.patient == p) & pre_mask, "tumor site"].iloc[0]
        post_site = pheno.loc[(pheno.patient == p) & post_mask, "tumor site"].iloc[0]
        same = pre_site == post_site
        same_site += int(same)
        site_notes.append({"patient": p, "pre_site": pre_site, "post_site": post_site, "same_site": same})

    tests = []
    for gene, probe in PROBES.items():
        y = expr[probe]
        pre = y.loc[pheno.loc[pre_mask, "geo_accession"]].to_numpy()
        post = y.loc[pheno.loc[post_mask, "geo_accession"]].to_numpy()
        u = stats.mannwhitneyu(post, pre, alternative="two-sided")
        # paired
        pre_p = []
        post_p = []
        for p in paired_patients:
            pre_p.append(y.loc[pheno.loc[(pheno.patient == p) & pre_mask, "geo_accession"]].iloc[0])
            post_p.append(y.loc[pheno.loc[(pheno.patient == p) & post_mask, "geo_accession"]].iloc[0])
        pre_p = np.array(pre_p)
        post_p = np.array(post_p)
        w = stats.wilcoxon(post_p, pre_p, alternative="two-sided", zero_method="wilcox")
        delta = post_p - pre_p
        tests.append(
            {
                "gene": gene,
                "probe_id": probe,
                "contrast": "unpaired_post_vs_pre",
                "n_pre": n_pre,
                "n_post": n_post,
                "n_patients": n_patients,
                "pre_mean": float(np.mean(pre)),
                "post_mean": float(np.mean(post)),
                "pre_median": float(np.median(pre)),
                "post_median": float(np.median(post)),
                "delta_mean": float(np.mean(post) - np.mean(pre)),
                "delta_median": float(np.median(post) - np.median(pre)),
                "U": float(u.statistic),
                "p": float(u.pvalue),
                "cliffs_delta_post_minus_pre": cliffs_delta(post, pre),
            }
        )
        tests.append(
            {
                "gene": gene,
                "probe_id": probe,
                "contrast": "paired_post_minus_pre",
                "n_pairs": n_paired,
                "n_same_site_pairs": same_site,
                "pre_mean": float(np.mean(pre_p)),
                "post_mean": float(np.mean(post_p)),
                "pre_median": float(np.median(pre_p)),
                "post_median": float(np.median(post_p)),
                "delta_mean": float(np.mean(delta)),
                "delta_median": float(np.median(delta)),
                "n_post_higher": int((delta > 0).sum()),
                "n_pre_higher": int((delta < 0).sum()),
                "W": float(w.statistic),
                "p": float(w.pvalue),
            }
        )
        same_idx = [i for i, p in enumerate(paired_patients) if site_notes[i]["same_site"]]
        if same_idx:
            pre_s = pre_p[same_idx]
            post_s = post_p[same_idx]
            d_s = post_s - pre_s
            w_s = stats.wilcoxon(post_s, pre_s, alternative="two-sided", zero_method="wilcox")
            tests.append(
                {
                    "gene": gene,
                    "probe_id": probe,
                    "contrast": "paired_same_site_only",
                    "n_pairs": int(len(same_idx)),
                    "n_same_site_pairs": int(len(same_idx)),
                    "pre_mean": float(np.mean(pre_s)),
                    "post_mean": float(np.mean(post_s)),
                    "pre_median": float(np.median(pre_s)),
                    "post_median": float(np.median(post_s)),
                    "delta_mean": float(np.mean(d_s)),
                    "delta_median": float(np.median(d_s)),
                    "n_post_higher": int((d_s > 0).sum()),
                    "n_pre_higher": int((d_s < 0).sum()),
                    "W": float(w_s.statistic),
                    "p": float(w_s.pvalue),
                }
            )

    # Spearman between the two genes (sample-level)
    rho_all = stats.spearmanr(expr[PROBES["CLDN4"]], expr[PROBES["TACSTD2"]])
    rho_pre = stats.spearmanr(
        expr.loc[pheno.loc[pre_mask, "geo_accession"], PROBES["CLDN4"]],
        expr.loc[pheno.loc[pre_mask, "geo_accession"], PROBES["TACSTD2"]],
    )
    rho_post = stats.spearmanr(
        expr.loc[pheno.loc[post_mask, "geo_accession"], PROBES["CLDN4"]],
        expr.loc[pheno.loc[post_mask, "geo_accession"], PROBES["TACSTD2"]],
    )

    control_tests = []
    for gene, probe in CONTROLS.items():
        if probe not in expr.columns:
            control_tests.append({"gene": gene, "probe_id": probe, "status": "probe_not_verified_skip"})
            continue
        y = expr[probe]
        pre = y.loc[pheno.loc[pre_mask, "geo_accession"]].to_numpy()
        post = y.loc[pheno.loc[post_mask, "geo_accession"]].to_numpy()
        u = stats.mannwhitneyu(post, pre, alternative="two-sided")
        control_tests.append(
            {
                "gene": gene,
                "probe_id": probe,
                "n_pre": n_pre,
                "n_post": n_post,
                "pre_median": float(np.median(pre)),
                "post_median": float(np.median(post)),
                "p": float(u.pvalue),
                "note": "sanity only; probe IDs from GPL23126 official-symbol scan if present in matrix",
            }
        )

    # GEO has no responder / sensitive label
    resp_fields = []
    for k, v in meta.items():
        blob = " ".join(v).lower()
        if any(w in blob for w in ("recist", "responder", "sensitive", "resistant", "cr/pr", "best response")):
            resp_fields.append(k)

    summary = {
        "gse": "GSE248249",
        "pmid": "38215748",
        "platform": "GPL23126 Affymetrix Clariom D Human (transcript/gene version)",
        "assay": "array (not NanoString)",
        "values": "RMA from Affymetrix Expression Console (as deposited in series matrix)",
        "n_features_in_matrix": n_features,
        "n_samples": int(pheno.shape[0]),
        "n_patients": n_patients,
        "n_pre": n_pre,
        "n_post": n_post,
        "n_paired_patients": n_paired,
        "n_same_site_pairs": same_site,
        "sensitive_vs_resistant_in_GEO": False,
        "response_fields_in_GEO": resp_fields,
        "cldn4_on_panel": "CLDN4" not in missing_targets,
        "tacstd2_on_panel": "TACSTD2" not in [g for g, p in PROBES.items() if p not in found_ids],
        "probes": PROBES,
        "cldn4_tacstd2_spearman_all": {"n": int(pheno.shape[0]), "rho": float(rho_all.statistic), "p": float(rho_all.pvalue)},
        "cldn4_tacstd2_spearman_pre": {"n": n_pre, "rho": float(rho_pre.statistic), "p": float(rho_pre.pvalue)},
        "cldn4_tacstd2_spearman_post": {"n": n_post, "rho": float(rho_post.statistic), "p": float(rho_post.pvalue)},
        "series_overall_design": meta.get("Series_overall_design", [""])[0],
    }

    pheno.to_csv(OUT / "pheno.tsv", sep="\t", index=False)
    long.to_csv(OUT / "sample_level.tsv", sep="\t", index=False)
    pd.DataFrame(tests).to_csv(OUT / "tests.tsv", sep="\t", index=False)
    pd.DataFrame(site_notes).to_csv(OUT / "paired_sites.tsv", sep="\t", index=False)
    pd.DataFrame(control_tests).to_csv(OUT / "control_sanity.tsv", sep="\t", index=False)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print("--- tests ---")
    print(pd.DataFrame(tests).to_string(index=False))


if __name__ == "__main__":
    main()
