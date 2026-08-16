#!/usr/bin/env python3
"""Compute TACSTD2 / CLDN4 vs ICI response in lung blood/PBMC/plasma/CTC series.

Do not drop a series for being blood. If the gene is in the public matrix, test it.
"""
from __future__ import annotations

import gzip
import json
import math
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

DATA = Path("/tmp/geo_blood")
OUT = Path("results/noskip/blood_all")
OUT.mkdir(parents=True, exist_ok=True)

PANEL = ["TACSTD2", "CLDN4", "EPCAM", "KRT8", "KRT18", "PTPRC", "CD3D", "CD8A", "GZMB", "PRF1"]
GPL570 = {
    "TACSTD2": ["202286_s_at", "202287_s_at"],
    "CLDN4": ["201428_at", "1569421_at"],
    "EPCAM": ["201839_s_at"],
    "KRT8": ["209008_x_at"],
    "KRT18": ["201596_x_at"],
    "PTPRC": ["212587_s_at"],
    "CD3D": ["213539_at"],
    "CD8A": ["205758_at"],
    "GZMB": ["210164_at"],
    "PRF1": ["214617_at"],
}
CLARIOM = {
    "TACSTD2": ["TC0100014340.hg.1"],
    "CLDN4": ["TC0700007993.hg.1"],
    "EPCAM": ["TC0200007506.hg.1"],
}


def mw(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return {"n1": int(len(a)), "n2": int(len(b)), "median1": None, "median2": None, "U": None, "p": None}
    U, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    return {
        "n1": int(len(a)),
        "n2": int(len(b)),
        "median1": float(np.median(a)),
        "median2": float(np.median(b)),
        "mean1": float(np.mean(a)),
        "mean2": float(np.mean(b)),
        "U": float(U),
        "p": float(p),
    }


def parse_series_matrix(path: Path) -> tuple[pd.DataFrame, pd.DataFrame | None]:
    """Return (sample_meta DataFrame indexed by GSM, optional expression DataFrame)."""
    meta_rows = {}
    expr_header = None
    expr_ids = []
    expr_data = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("!Sample_"):
                key, *vals = line.rstrip("\n").split("\t")
                key = key.replace("!Sample_", "")
                vals = [v.strip().strip('"') for v in vals]
                if key in meta_rows and key.startswith("characteristics"):
                    # GEO repeats the same key; keep every line
                    n = sum(1 for k in meta_rows if k.startswith("characteristics"))
                    key = f"characteristics_{n}"
                meta_rows[key] = vals
            elif line.startswith("ID_REF") or line.startswith('"ID_REF"'):
                parts = [p.strip().strip('"') for p in line.rstrip("\n").split("\t")]
                expr_header = parts[1:]
            elif expr_header is not None and line.strip() and not line.startswith("!"):
                parts = [p.strip().strip('"') for p in line.rstrip("\n").split("\t")]
                if len(parts) >= 2:
                    expr_ids.append(parts[0])
                    expr_data.append(parts[1 : 1 + len(expr_header)])
    # build meta
    gsms = meta_rows.get("geo_accession") or meta_rows.get("geo_accession\t") 
    if gsms is None:
        # key may be geo_accession
        for k in meta_rows:
            if k.startswith("geo_accession"):
                gsms = meta_rows[k]
                break
    if not gsms:
        raise RuntimeError(f"no GSM in {path}")
    n = len(gsms)
    meta = pd.DataFrame(index=gsms)
    char_i = 0
    for k, vals in meta_rows.items():
        vals = (vals + [""] * n)[:n]
        if k.startswith("characteristics"):
            # split field: value
            parsed = []
            field = None
            for v in vals:
                if ": " in v:
                    field, val = v.split(": ", 1)
                    parsed.append(val)
                else:
                    parsed.append(v)
            col = (field or f"char{char_i}").strip()
            char_i += 1
            meta[col] = parsed
        else:
            meta[k] = vals
    expr = None
    if expr_header and expr_data:
        arr = np.array(expr_data, dtype=object)
        # numeric where possible
        try:
            num = arr.astype(float)
        except Exception:
            num = pd.DataFrame(arr, index=expr_ids, columns=expr_header).apply(pd.to_numeric, errors="coerce").to_numpy()
        expr = pd.DataFrame(num, index=expr_ids, columns=expr_header)
    return meta, expr


def log2p1_cpm(counts: pd.Series, lib: pd.Series | None = None) -> pd.Series:
    if lib is None:
        # if already intensity-like, just return
        return counts
    cpm = counts / lib.replace(0, np.nan) * 1e6
    return np.log2(cpm + 1)


def boxplot(series_by_group: dict, title: str, ylabel: str, dest: Path):
    labels = list(series_by_group)
    data = [np.asarray(series_by_group[k], float) for k in labels]
    fig, ax = plt.subplots(figsize=(5.2, 4.0))
    ax.boxplot(data, tick_labels=labels, showfliers=True)
    for i, d in enumerate(data, 1):
        x = np.random.default_rng(0).normal(i, 0.06, size=len(d))
        ax.scatter(x, d, s=18, alpha=0.7, c="#333")
    ax.set_title(title, fontsize=10)
    ax.set_ylabel(ylabel)
    fig.tight_layout()
    fig.savefig(dest, dpi=140)
    plt.close(fig)


def record(rows, **kw):
    rows.append(kw)


# ---------------------------------------------------------------------------
# GSE111414  CD8 PBMC RNA-seq, responder vs non-responder
# ---------------------------------------------------------------------------
def analyze_gse111414(rows, per_sample):
    raw = pd.read_csv(DATA / "GSE111414_gene_counts.csv.gz")
    # ENTREZ, SYMBOL, GENENAME, LENGTH, then samples
    if "SYMBOL" in raw.columns:
        counts = raw.set_index("SYMBOL")
        sample_cols = [c for c in counts.columns if c not in ("GENENAME", "LENGTH") and not str(c).startswith("Unnamed")]
        counts = counts[sample_cols].apply(pd.to_numeric, errors="coerce").groupby(level=0).sum()
    else:
        counts = raw.set_index(raw.columns[0])
    meta, _ = parse_series_matrix(DATA / "GSE111414_series_matrix.txt.gz")
    meta = meta.copy()
    meta["sample"] = meta["description"]
    # source_name is responder_N1 / nonresponder_N2
    src = meta["source_name_ch1"].astype(str)
    meta["response"] = src.str.replace(r"_N\d+$", "", regex=True)
    meta["timepoint"] = src.str.extract(r"(N\d+)$")[0]
    resp_col = "response"
    time_col = "timepoint"
    lib = counts.sum(axis=0)
    genes_present = [g for g in PANEL if g in counts.index]
    # baseline N1 only for vs-response (one value per patient)
    base = meta[meta[time_col].astype(str).str.upper().eq("N1")].copy()
    for g in genes_present:
        logcpm = np.log2(counts.loc[g] / lib * 1e6 + 1)
        for _, r in meta.iterrows():
            sid = r["sample"]
            if sid not in logcpm.index:
                continue
            per_sample.append(
                {
                    "dataset": "GSE111414",
                    "compartment": "PBMC_CD8",
                    "sample": sid,
                    "patient": r.get("subject id", ""),
                    "timepoint": r.get(time_col, ""),
                    "response_raw": r.get(resp_col, ""),
                    "gene": g,
                    "value": float(logcpm[sid]),
                    "metric": "log2CPM",
                }
            )
        # patient-level baseline
        vals = {}
        for _, r in base.iterrows():
            sid = r["sample"]
            if sid in logcpm.index:
                raw = str(r[resp_col]).lower()
                grp = "responder" if raw.startswith("responder") else "nonresponder"
                vals.setdefault(grp, []).append(float(logcpm[sid]))
        st = mw(vals.get("responder", []), vals.get("nonresponder", []))
        record(
            rows,
            dataset="GSE111414",
            compartment="PBMC_CD8",
            gene=g,
            contrast="baseline_N1 responder vs nonresponder",
            metric="log2CPM",
            genes_in_matrix=True,
            n_responder=st["n1"],
            n_nonresponder=st["n2"],
            median_responder=st["median1"],
            median_nonresponder=st["median2"],
            p=st["p"],
            note="epithelial genes expected near-zero in sorted CD8 T cells",
        )
        if g in ("TACSTD2", "CLDN4") and vals:
            boxplot(
                {"responder": vals.get("responder", []), "nonresponder": vals.get("nonresponder", [])},
                f"GSE111414 {g} CD8 PBMC baseline",
                "log2 CPM",
                OUT / f"GSE111414_{g}_response.png",
            )
    return genes_present


# ---------------------------------------------------------------------------
# GSE202417  CD8 PBMC Clariom array, R vs NR, pre/post bezafibrate+nivo
# ---------------------------------------------------------------------------
def analyze_gse202417(rows, per_sample):
    meta, expr = parse_series_matrix(DATA / "GSE202417_series_matrix.txt.gz")
    if expr is None:
        record(rows, dataset="GSE202417", gene="TACSTD2", genes_in_matrix=False, note="no expression in series matrix")
        return []
    present = []
    # phenotype R/NR, time pre/post
    pheno = None
    for c in meta.columns:
        if c.lower() in ("phenotype", "response"):
            pheno = c
    time_c = None
    for c in meta.columns:
        if c.lower() in ("time", "timepoint"):
            time_c = c
    for gene, probes in CLARIOM.items():
        probes = [p for p in probes if p in expr.index]
        if not probes:
            record(rows, dataset="GSE202417", compartment="PBMC_CD8", gene=gene, genes_in_matrix=False, note="probe absent")
            continue
        present.append(gene)
        val = expr.loc[probes].astype(float).mean(axis=0)
        # align to GSM
        val.index = [str(x) for x in val.index]
        meta2 = meta.copy()
        meta2.index = [str(x) for x in meta2.index]
        for gsm in meta2.index:
            if gsm not in val.index:
                continue
            per_sample.append(
                {
                    "dataset": "GSE202417",
                    "compartment": "PBMC_CD8",
                    "sample": gsm,
                    "patient": meta2.loc[gsm].get("title", ""),
                    "timepoint": meta2.loc[gsm].get(time_c, "") if time_c else "",
                    "response_raw": meta2.loc[gsm].get(pheno, "") if pheno else "",
                    "gene": gene,
                    "value": float(val[gsm]),
                    "metric": "array_intensity",
                    "probe": ",".join(probes),
                }
            )
        # baseline pre-treatment, one per patient (title ComboN-R/NR-1 is pre)
        pre = meta2
        if time_c:
            pre = meta2[meta2[time_c].astype(str).str.contains("pre", case=False, na=False)]
        groups = {"R": [], "NR": []}
        for gsm, r in pre.iterrows():
            if gsm not in val.index:
                continue
            lab = str(r.get(pheno, "")).upper()
            if lab in groups:
                groups[lab].append(float(val[gsm]))
        st = mw(groups["R"], groups["NR"])
        record(
            rows,
            dataset="GSE202417",
            compartment="PBMC_CD8",
            gene=gene,
            contrast="pre-treatment R vs NR (bezafibrate+nivolumab trial)",
            metric="array_intensity",
            genes_in_matrix=True,
            n_responder=st["n1"],
            n_nonresponder=st["n2"],
            median_responder=st["median1"],
            median_nonresponder=st["median2"],
            p=st["p"],
            note="sorted CD8; probe " + ",".join(probes),
        )
        if gene in ("TACSTD2", "CLDN4"):
            boxplot(
                {"R": groups["R"], "NR": groups["NR"]},
                f"GSE202417 {gene} CD8 PBMC pre-tx",
                "array intensity",
                OUT / f"GSE202417_{gene}_response.png",
            )
    return present


# ---------------------------------------------------------------------------
# GSE249262  CTC microarray, progression vs stable, CRT + durvalumab
# ---------------------------------------------------------------------------
def analyze_gse249262(rows, per_sample):
    meta, expr = parse_series_matrix(DATA / "GSE249262_series_matrix.txt.gz")
    if expr is None:
        return []
    present = []
    status_c = None
    time_c = None
    for c in meta.columns:
        if c.lower() == "status":
            status_c = c
        if "timepoint" in c.lower() or c.lower() == "timepoint":
            time_c = c
    for gene, probes in CLARIOM.items():
        probes = [p for p in probes if p in expr.index]
        if not probes:
            record(rows, dataset="GSE249262", compartment="CTC", gene=gene, genes_in_matrix=False, note="probe absent")
            continue
        present.append(gene)
        val = expr.loc[probes].astype(float).mean(axis=0)
        val.index = [str(x) for x in val.index]
        meta2 = meta.copy()
        meta2.index = [str(x) for x in meta2.index]
        for gsm in meta2.index:
            if gsm not in val.index:
                continue
            per_sample.append(
                {
                    "dataset": "GSE249262",
                    "compartment": "CTC",
                    "sample": gsm,
                    "patient": meta2.loc[gsm].get("title", ""),
                    "timepoint": meta2.loc[gsm].get(time_c, "") if time_c else "",
                    "response_raw": meta2.loc[gsm].get(status_c, "") if status_c else "",
                    "gene": gene,
                    "value": float(val[gsm]),
                    "metric": "array_intensity",
                    "probe": ",".join(probes),
                }
            )
        # baseline samples, patient-level
        use = meta2
        if time_c:
            use = meta2[meta2[time_c].astype(str).str.contains("base", case=False, na=False)]
        groups = {"stable": [], "progression": []}
        for gsm, r in use.iterrows():
            if gsm not in val.index:
                continue
            lab = str(r.get(status_c, "")).lower()
            if "stable" in lab:
                groups["stable"].append(float(val[gsm]))
            elif "progress" in lab:
                groups["progression"].append(float(val[gsm]))
        st = mw(groups["stable"], groups["progression"])
        record(
            rows,
            dataset="GSE249262",
            compartment="CTC",
            gene=gene,
            contrast="baseline CTC stable vs progression (CRT + durvalumab, PMID 38261515)",
            metric="array_intensity",
            genes_in_matrix=True,
            n_responder=st["n1"],
            n_nonresponder=st["n2"],
            median_responder=st["median1"],
            median_nonresponder=st["median2"],
            p=st["p"],
            note="responder=stable; nonresponder=progression; CTC RNA is the biologically plausible blood compartment",
        )
        if gene in ("TACSTD2", "CLDN4"):
            boxplot(
                {"stable": groups["stable"], "progression": groups["progression"]},
                f"GSE249262 {gene} CTC baseline",
                "array intensity",
                OUT / f"GSE249262_{gene}_response.png",
            )
    return present


# ---------------------------------------------------------------------------
# GSE235048  PBMC RNA-seq TPM, ICI pre/post, NO response labels
# ---------------------------------------------------------------------------
def analyze_gse235048(rows, per_sample):
    expr = pd.read_csv(DATA / "GSE235048_24CFlow_TPM.txt.gz", sep="\t", index_col=0)
    meta, _ = parse_series_matrix(DATA / "GSE235048_series_matrix.txt.gz")
    # description is A2, B1...
    desc = meta["description"] if "description" in meta.columns else None
    meta = meta.copy()
    meta["well"] = desc.values if desc is not None else meta["title"].str.extract(r"\[([A-Z]\d+)\]")[0].values
    present = []
    for g in PANEL:
        if g not in expr.index:
            continue
        present.append(g)
        v = expr.loc[g].astype(float)
        for well, val in v.items():
            hit = meta[meta["well"].astype(str) == str(well)]
            rec = hit.iloc[0] if len(hit) else None
            per_sample.append(
                {
                    "dataset": "GSE235048",
                    "compartment": "PBMC",
                    "sample": str(well),
                    "patient": rec.get("title", "") if rec is not None else "",
                    "timepoint": rec.get("title", "") if rec is not None else "",
                    "response_raw": "",
                    "gene": g,
                    "value": float(val),
                    "metric": "TPM",
                }
            )
        # pre vs post from title
        pre, post = [], []
        for _, r in meta.iterrows():
            well = str(r["well"])
            if well not in v.index:
                continue
            t = str(r.get("title", "")).lower()
            if "pre-treatment" in t:
                pre.append(float(v[well]))
            elif "post-treatment" in t:
                post.append(float(v[well]))
        st = mw(pre, post)
        record(
            rows,
            dataset="GSE235048",
            compartment="PBMC",
            gene=g,
            contrast="pre vs post ICI (NO public RECIST/response label)",
            metric="TPM",
            genes_in_matrix=True,
            n_responder=st["n1"],
            n_nonresponder=st["n2"],
            median_responder=st["median1"],
            median_nonresponder=st["median2"],
            p=st["p"],
            note="groups are timepoint not response; n1=pre n2=post",
        )
    return present


# ---------------------------------------------------------------------------
# GSE225620  whole-blood featureCounts, neoadjuvant PD-1, NO public response
# ---------------------------------------------------------------------------
def analyze_gse225620(rows, per_sample):
    expr = pd.read_csv(DATA / "GSE225620_pre_vs_post_featureCounts.txt.gz", sep="\t", index_col=0)
    lib = expr.sum(axis=0)
    present = []
    for g in PANEL:
        if g not in expr.index:
            continue
        present.append(g)
        logcpm = np.log2(expr.loc[g] / lib * 1e6 + 1)
        pre = [c for c in logcpm.index if str(c).startswith("pre")]
        post = [c for c in logcpm.index if str(c).startswith("post")]
        for c, val in logcpm.items():
            per_sample.append(
                {
                    "dataset": "GSE225620",
                    "compartment": "whole_blood",
                    "sample": str(c),
                    "patient": "",
                    "timepoint": "pre" if str(c).startswith("pre") else ("post" if str(c).startswith("post") else "other"),
                    "response_raw": "",
                    "gene": g,
                    "value": float(val),
                    "metric": "log2CPM",
                }
            )
        st = mw(logcpm[pre], logcpm[post])
        record(
            rows,
            dataset="GSE225620",
            compartment="whole_blood",
            gene=g,
            contrast="pre vs post neoadjuvant PD-1 (NO public responder label in GEO)",
            metric="log2CPM",
            genes_in_matrix=True,
            n_responder=st["n1"],
            n_nonresponder=st["n2"],
            median_responder=st["median1"],
            median_nonresponder=st["median2"],
            p=st["p"],
            note="design text says responders exist but GEO characteristics have only timepoint/sex/age",
        )
        if g in ("TACSTD2", "CLDN4"):
            boxplot(
                {"pre": list(logcpm[pre]), "post": list(logcpm[post])},
                f"GSE225620 {g} whole blood",
                "log2 CPM",
                OUT / f"GSE225620_{g}_prepost.png",
            )
    return present


# ---------------------------------------------------------------------------
# GSE305086  whole-blood GPL570, baseline vs control / paired follow-up
# ---------------------------------------------------------------------------
def analyze_gse305086(rows, per_sample):
    # processed matrix is probe x sample, semicolon-separated
    expr = pd.read_csv(DATA / "GSE305086_Expression_matrix_final.csv.gz", sep=";", index_col=0)
    meta, _ = parse_series_matrix(DATA / "GSE305086_series_matrix.txt.gz")
    present = []
    # try to classify samples from columns / meta
    # columns of expr may be sample titles
    for gene, probes in GPL570.items():
        probes = [p for p in probes if p in expr.index]
        if not probes:
            # try symbol index
            if gene in expr.index:
                val = expr.loc[gene].astype(float)
                probes = [gene]
            else:
                continue
        else:
            val = expr.loc[probes].astype(float).mean(axis=0)
        present.append(gene)
        for sid, x in val.items():
            per_sample.append(
                {
                    "dataset": "GSE305086",
                    "compartment": "whole_blood",
                    "sample": str(sid),
                    "patient": "",
                    "timepoint": "",
                    "response_raw": "",
                    "gene": gene,
                    "value": float(x),
                    "metric": "log2_RMA",
                    "probe": ",".join(map(str, probes)),
                }
            )
        # detectability: percentile of mean among all probes
        if probes[0] in expr.index and expr.shape[0] > 100:
            means = expr.astype(float).mean(axis=1)
            pct = float(stats.percentileofscore(means, float(means.loc[probes[0]] if probes[0] in means.index else val.mean())))
        else:
            pct = None
        record(
            rows,
            dataset="GSE305086",
            compartment="whole_blood",
            gene=gene,
            contrast="detectability only — NO public ICI response/PFS label",
            metric="log2_RMA_mean_percentile",
            genes_in_matrix=True,
            n_responder=int(val.shape[0]),
            n_nonresponder=0,
            median_responder=float(np.median(val)),
            median_nonresponder=None,
            p=None,
            note=f"probe {','.join(map(str, probes))}; mean-intensity percentile among probes={pct}",
        )
    return present


# ---------------------------------------------------------------------------
# GSE152590  CD8 PBL TPM, 4 lung + 4 melanoma, NO response
# ---------------------------------------------------------------------------
def analyze_gse152590(rows, per_sample):
    expr = pd.read_excel(DATA / "GSE152590_TPM_matrix.xlsx", index_col=0)
    # index may be gene symbols or ensembl
    idx_map = {str(i).split(".")[0].upper(): i for i in expr.index}
    present = []
    for g in PANEL:
        key = None
        if g in expr.index:
            key = g
        else:
            for k, orig in idx_map.items():
                if k == g or k.endswith("|" + g) or g in k.split("|"):
                    key = orig
                    break
        if key is None:
            continue
        present.append(g)
        v = expr.loc[key].astype(float)
        for sid, x in v.items():
            per_sample.append(
                {
                    "dataset": "GSE152590",
                    "compartment": "PBMC_CD8",
                    "sample": str(sid),
                    "patient": str(sid),
                    "timepoint": "",
                    "response_raw": "",
                    "gene": g,
                    "value": float(x),
                    "metric": "TPM",
                }
            )
        record(
            rows,
            dataset="GSE152590",
            compartment="PBMC_CD8",
            gene=g,
            contrast="detectability only — 4 lung + 4 melanoma, NO response label",
            metric="TPM",
            genes_in_matrix=True,
            n_responder=int(v.shape[0]),
            n_nonresponder=0,
            median_responder=float(np.median(v)),
            median_nonresponder=None,
            p=None,
            note=f"mean TPM={float(v.mean()):.4g}; n_nonzero={int((v>0).sum())}",
        )
    return present


def main():
    rows = []
    per_sample = []
    presence = {}
    print("GSE111414", flush=True)
    presence["GSE111414"] = analyze_gse111414(rows, per_sample)
    print("GSE202417", flush=True)
    presence["GSE202417"] = analyze_gse202417(rows, per_sample)
    print("GSE249262", flush=True)
    presence["GSE249262"] = analyze_gse249262(rows, per_sample)
    print("GSE235048", flush=True)
    presence["GSE235048"] = analyze_gse235048(rows, per_sample)
    print("GSE225620", flush=True)
    presence["GSE225620"] = analyze_gse225620(rows, per_sample)
    print("GSE305086", flush=True)
    presence["GSE305086"] = analyze_gse305086(rows, per_sample)
    print("GSE152590", flush=True)
    presence["GSE152590"] = analyze_gse152590(rows, per_sample)

    pd.DataFrame(rows).to_csv(OUT / "response_stats.tsv", sep="\t", index=False)
    pd.DataFrame(per_sample).to_csv(OUT / "per_sample_expression.tsv", sep="\t", index=False)
    (OUT / "gene_presence.json").write_text(json.dumps(presence, indent=2))
    print("wrote", OUT / "response_stats.tsv", "n", len(rows))
    # print TACSTD2/CLDN4 rows
    for r in rows:
        if r.get("gene") in ("TACSTD2", "CLDN4"):
            print(f"  {r['dataset']} {r['gene']} {r.get('contrast')} p={r.get('p')} n={r.get('n_responder')}/{r.get('n_nonresponder')} med={r.get('median_responder')}/{r.get('median_nonresponder')}")


if __name__ == "__main__":
    main()
