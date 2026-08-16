#!/usr/bin/env python3
"""Analyze leftover 2026-only GEO lung ICI series for TACSTD2/CLDN4.

Usable leftover tumor cohorts:
  GSE329813  GeoMx DSP, 22 NSCLC patients, MPR/NMPR in sample titles
             TACSTD2 present; CLDN4 absent from the processed panel
  GSE292299  Visium, 16 NSCLC patients, Tx_Response in sample_metadata.csv
             per-sample H5 files are each < 2 GB (the 3.1 GB RAW.tar is not used)

Outputs land in results/w200/GEO_2026/.
"""
import gzip
import hashlib
import io
import json
import math
import re
import time
import urllib.request
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

RES = Path(__file__).resolve().parents[2] / "results" / "w200" / "GEO_2026"
DATA = RES / "data"
RESULTS = []


def fetch(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        print(f"  have {dest.name} ({dest.stat().st_size/1e6:.1f} MB)")
        return
    last = None
    for attempt in range(5):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "geo-2026"})
            with urllib.request.urlopen(req, timeout=300) as r, dest.open("wb") as fh:
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    fh.write(chunk)
            print(f"  OK {dest.name} ({dest.stat().st_size/1e6:.1f} MB)")
            return
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(2 ** attempt)
    raise RuntimeError(f"{url}: {last}")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def mannwhitney(a, b, method="auto"):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    if len(a) < 2 or len(b) < 2:
        return float("nan"), float("nan")
    return stats.mannwhitneyu(a, b, alternative="two-sided", method=method)


def cliffs_delta(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    if len(a) == 0 or len(b) == 0:
        return float("nan")
    gt = sum(1 for x in a for y in b if x > y)
    lt = sum(1 for x in a for y in b if x < y)
    return (gt - lt) / (len(a) * len(b))


def boxplot(groups, labels, title, ylab, path):
    fig, ax = plt.subplots(figsize=(4.2, 4))
    bp = ax.boxplot(groups, tick_labels=labels, showfliers=False,
                    patch_artist=True, widths=0.6)
    for patch, c in zip(bp["boxes"], ["#4C9F70", "#C0504D"]):
        patch.set_facecolor(c)
        patch.set_alpha(0.45)
    rng = np.random.default_rng(0)
    for i, g in enumerate(groups, 1):
        ax.scatter(rng.normal(i, 0.05, size=len(g)), g, s=18,
                   color="black", zorder=3, alpha=0.7)
    ax.set_title(title, fontsize=10)
    ax.set_ylabel(ylab)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def add_row(**kwargs):
    RESULTS.append(kwargs)
    print(" ", {k: kwargs[k] for k in (
        "dataset", "gene", "comparison", "n_group1", "n_group2",
        "p_value", "effect_value")})


def parse_series_titles(gse: str):
    stub = gse[:-3] + "nnn"
    url = (f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/{gse}/matrix/"
           f"{gse}_series_matrix.txt.gz")
    raw = urllib.request.urlopen(
        urllib.request.Request(url, headers={"User-Agent": "geo-2026"}),
        timeout=90).read()
    txt = gzip.GzipFile(fileobj=io.BytesIO(raw)).read().decode("utf-8", "replace")
    titles = accs = None
    for ln in txt.splitlines():
        if ln.startswith("!Sample_title"):
            titles = [x.strip().strip('"') for x in ln.split("\t")[1:]]
        if ln.startswith("!Sample_geo_accession"):
            accs = [x.strip().strip('"') for x in ln.split("\t")[1:]]
    return pd.DataFrame({"geo_accession": accs, "title": titles})


# ---------------------------------------------------------------------------
# GSE329813
# ---------------------------------------------------------------------------
def analyze_gse329813():
    print("\n" + "=" * 70 + "\n GSE329813")
    url = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE329nnn/GSE329813/suppl/"
           "GSE329813_processed_data_file_normalized_data.csv.gz")
    dest = DATA / "GSE329813" / "GSE329813_processed_data_file_normalized_data.csv.gz"
    fetch(url, dest)
    expr = pd.read_csv(dest)
    expr = expr.set_index(expr.columns[0])
    expr.index = expr.index.astype(str).str.upper()
    genes = {g: (g in expr.index) for g in ("TACSTD2", "CLDN4")}
    print(f"  genes present: {genes}; n_genes={expr.shape[0]} n_roi={expr.shape[1]}")

    meta = parse_series_titles("GSE329813")
    # titles: "ROI 1, Patient 1, Primary tumor bed, MPR"
    pat = []
    for t in meta["title"]:
        m = re.match(
            r"(ROI \d+),\s*Patient (\d+),\s*(Primary tumor bed|Lymph node),\s*(MPR|NMPR)$",
            t)
        if not m:
            raise ValueError(f"unparsed title: {t}")
        pat.append({"roi": m.group(1), "patient": int(m.group(2)),
                    "tissue": m.group(3), "mpr": m.group(4)})
    meta = pd.concat([meta, pd.DataFrame(pat)], axis=1)
    missing = [c for c in meta.roi if c not in expr.columns]
    print(f"  ROI columns missing from matrix: {missing}")

    rows = []
    for _, r in meta.iterrows():
        if r.roi not in expr.columns:
            continue
        rec = {"roi": r.roi, "patient": r.patient, "tissue": r.tissue,
               "mpr": r.mpr, "geo_accession": r.geo_accession}
        for g in ("TACSTD2", "CLDN4"):
            rec[g] = float(expr.loc[g, r.roi]) if g in expr.index else np.nan
        rows.append(rec)
    roi = pd.DataFrame(rows)
    roi.to_csv(RES / "GSE329813_roi_TACSTD2.tsv", sep="\t", index=False)

    # patient-level means, primary tumor bed first (pre-specified)
    for tissue, label in (
            ("Primary tumor bed", "primary tumor bed"),
            ("Lymph node", "draining lymph node")):
        sub = roi[roi.tissue == tissue]
        pat = sub.groupby(["patient", "mpr"], as_index=False).agg(
            TACSTD2=("TACSTD2", "mean"),
            n_roi=("roi", "count"),
        )
        pat.to_csv(RES / f"GSE329813_patient_{tissue.replace(' ', '_')}.tsv",
                   sep="\t", index=False)
        r = pat.loc[pat.mpr == "MPR", "TACSTD2"].dropna()
        nr = pat.loc[pat.mpr == "NMPR", "TACSTD2"].dropna()
        u, p = mannwhitney(r, nr)
        d = cliffs_delta(r, nr)
        print(f"  TACSTD2 {label}: MPR n={len(r)} med={r.median():.3f}; "
              f"NMPR n={len(nr)} med={nr.median():.3f}; p={p:.3f} δ={d:.3f}")
        add_row(
            dataset="GSE329813", cohort="NSCLC neoadjuvant chemo-IO (GeoMx DSP)",
            gene="TACSTD2",
            comparison=f"MPR vs NMPR, patient-mean, {label}",
            test="Mann-Whitney U (two-sided)",
            n_group1=len(r), n_group2=len(nr),
            median_group1=round(float(r.median()), 4) if len(r) else None,
            median_group2=round(float(nr.median()), 4) if len(nr) else None,
            statistic=round(float(u), 4) if not math.isnan(u) else None,
            p_value=round(float(p), 4) if not math.isnan(p) else None,
            effect_type="Cliff's delta (MPR minus NMPR)",
            effect_value=round(float(d), 4) if not math.isnan(d) else None,
        )
        if tissue == "Primary tumor bed":
            boxplot([r.values, nr.values],
                    [f"MPR\n(n={len(r)})", f"NMPR\n(n={len(nr)})"],
                    "GSE329813: TACSTD2 vs MPR (tumor bed)",
                    "TACSTD2 (normalized DSP)",
                    RES / "GSE329813_TACSTD2_MPR.png")
    add_row(
        dataset="GSE329813", cohort="NSCLC neoadjuvant chemo-IO (GeoMx DSP)",
        gene="CLDN4", comparison="panel presence",
        test="gene inventory of processed DSP matrix",
        n_group1=None, n_group2=None,
        median_group1=None, median_group2=None, statistic=None,
        p_value=None, effect_type="absent from processed panel",
        effect_value=None,
    )


# ---------------------------------------------------------------------------
# GSE292299
# ---------------------------------------------------------------------------
NSCLC_H5 = [
    ("GSM8855703", "NSCLC_P1"), ("GSM8855704", "NSCLC_P2"),
    ("GSM8855705", "NSCLC_P3"), ("GSM8855706", "NSCLC_P4"),
    ("GSM8855707", "NSCLC_P5"), ("GSM8855708", "NSCLC_P6"),
    ("GSM8855709", "NSCLC_P7"), ("GSM8855710", "NSCLC_P8"),
    ("GSM8855711", "NSCLC_P9"), ("GSM8855712", "NSCLC_P10"),
    ("GSM8855713", "NSCLC_P11"), ("GSM8855714", "NSCLC_P12"),
    ("GSM8855715", "NSCLC_P13"), ("GSM8855716", "NSCLC_P14"),
    ("GSM8855717", "NSCLC_P15"), ("GSM8855718", "NSCLC_P16"),
]


def read_10x_h5_gene_sums(path: Path, genes=("TACSTD2", "CLDN4")):
    with h5py.File(path, "r") as f:
        # 10x v3: matrix/{data,indices,indptr,shape,features/name}
        root = "matrix" if "matrix" in f else list(f.keys())[0]
        g = f[root]
        data = g["data"][:]
        indices = g["indices"][:]
        indptr = g["indptr"][:]
        shape = tuple(int(x) for x in g["shape"][:])
        if "features" in g and "name" in g["features"]:
            names = [x.decode() if isinstance(x, bytes) else str(x)
                     for x in g["features"]["name"][:]]
        elif "gene_names" in g:
            names = [x.decode() if isinstance(x, bytes) else str(x)
                     for x in g["gene_names"][:]]
        else:
            raise KeyError(f"no gene names in {path}")
        # CSC: genes x spots (10x stores genes as rows)
        n_genes, n_spots = shape
        from scipy.sparse import csc_matrix
        mat = csc_matrix((data, indices, indptr), shape=shape)
        total = float(mat.sum())
        name_u = [n.upper() for n in names]
        out = {"n_spots": n_spots, "total_counts": total}
        for gene in genes:
            if gene in name_u:
                i = name_u.index(gene)
                cpm = float(mat[i, :].sum()) / total * 1e6 if total else 0.0
                out[gene] = math.log2(cpm + 1)
            else:
                out[gene] = float("nan")
        return out


def analyze_gse292299():
    print("\n" + "=" * 70 + "\n GSE292299")
    meta_url = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE292nnn/GSE292299/"
                "suppl/GSE292299_sample_metadata.csv.gz")
    meta_path = DATA / "GSE292299" / "GSE292299_sample_metadata.csv.gz"
    fetch(meta_url, meta_path)
    meta = pd.read_csv(meta_path)
    nsclc = meta[meta.TumorType.eq("NSCLC")].copy()
    print("  NSCLC Tx_Response:", nsclc.Tx_Response.value_counts().to_dict())

    rows = []
    for gsm, pid in NSCLC_H5:
        stub = gsm[:-3] + "nnn"
        url = (f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{stub}/{gsm}/"
               f"suppl/{gsm}_{pid}_filtered_feature_bc_matrix.h5")
        dest = DATA / "GSE292299" / f"{gsm}_{pid}_filtered_feature_bc_matrix.h5"
        fetch(url, dest)
        vals = read_10x_h5_gene_sums(dest)
        rec = {"gsm": gsm, "patient": pid, **vals}
        m = nsclc[nsclc.ID_Patient.eq(pid)]
        if len(m) != 1:
            raise ValueError(f"metadata mismatch {pid}")
        rec.update({
            "response": m.Tx_Response.iloc[0],
            "tissue": m.Bx_Tissue.iloc[0],
            "pfs_months": float(m.PFS_Months.iloc[0]),
            "treatment": m.PostSample_Treatment.iloc[0],
        })
        rows.append(rec)
        print(f"  {pid} {rec['response']} {rec['tissue']}: "
              f"TACSTD2={rec['TACSTD2']:.2f} CLDN4={rec['CLDN4']:.2f} "
              f"spots={rec['n_spots']}")
    df = pd.DataFrame(rows)
    df.to_csv(RES / "GSE292299_NSCLC_pseudobulk.tsv", sep="\t", index=False)

    rho, pr = stats.spearmanr(df.TACSTD2, df.CLDN4)
    add_row(
        dataset="GSE292299", cohort="NSCLC Visium whole-sample pseudobulk",
        gene="TACSTD2~CLDN4", comparison="Marker concordance across patients",
        test="Spearman rank correlation",
        n_group1=len(df), n_group2=None,
        median_group1=None, median_group2=None,
        statistic=round(float(rho), 4), p_value=round(float(pr), 4),
        effect_type="Spearman rho", effect_value=round(float(rho), 4),
    )

    for g in ("TACSTD2", "CLDN4"):
        r = df.loc[df.response == "R", g].dropna()
        nr = df.loc[df.response == "NR", g].dropna()
        u, p = mannwhitney(r, nr, method="exact")
        d = cliffs_delta(r, nr)
        print(f"  {g} R n={len(r)} med={r.median():.3f}; "
              f"NR n={len(nr)} med={nr.median():.3f}; p={p:.3f} δ={d:.3f}")
        add_row(
            dataset="GSE292299", cohort="NSCLC Visium whole-sample pseudobulk",
            gene=g,
            comparison="Tx_Response R vs NR, one Visium slide per patient",
            test="Mann-Whitney U (two-sided, exact)",
            n_group1=len(r), n_group2=len(nr),
            median_group1=round(float(r.median()), 4),
            median_group2=round(float(nr.median()), 4),
            statistic=round(float(u), 4),
            p_value=round(float(p), 4),
            effect_type="Cliff's delta (R minus NR)",
            effect_value=round(float(d), 4),
        )
        boxplot([r.values, nr.values],
                [f"R\n(n={len(r)})", f"NR\n(n={len(nr)})"],
                f"GSE292299: {g} vs ICI response",
                "log2(whole-slide CPM+1)",
                RES / f"GSE292299_{g}_response.png")
        rho_p, pp = stats.spearmanr(df[g], df.pfs_months)
        add_row(
            dataset="GSE292299", cohort="NSCLC Visium whole-sample pseudobulk",
            gene=g, comparison="Spearman vs PFS months (censoring ignored)",
            test="Spearman rank correlation (censoring ignored)",
            n_group1=len(df), n_group2=None,
            median_group1=None, median_group2=None,
            statistic=round(float(rho_p), 4), p_value=round(float(pp), 4),
            effect_type="Spearman rho", effect_value=round(float(rho_p), 4),
        )


def write_manifest():
    rows = []
    for gse, files in {
        "GSE329813": ["GSE329813_processed_data_file_normalized_data.csv.gz"],
        "GSE292299": (
            ["GSE292299_sample_metadata.csv.gz"]
            + [f"{gsm}_{pid}_filtered_feature_bc_matrix.h5"
               for gsm, pid in NSCLC_H5]
        ),
    }.items():
        for name in files:
            path = DATA / gse / name
            if not path.exists():
                continue
            if name.startswith("GSM"):
                gsm = name.split("_", 1)[0]
                stub = gsm[:-3] + "nnn"
                url = (f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{stub}/"
                       f"{gsm}/suppl/{name}")
            else:
                stub = gse[:-3] + "nnn"
                url = (f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/{gse}/"
                       f"suppl/{name}")
            rows.append((gse, name, str(path.stat().st_size), sha256(path), url))
    with (RES / "leftover_2026_input_manifest.tsv").open("w") as fh:
        fh.write("gse\tfile\tbytes\tsha256\turl\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")


def main():
    analyze_gse329813()
    analyze_gse292299()
    write_manifest()
    out = pd.DataFrame(RESULTS)
    out.to_csv(RES / "leftover_2026_analysis_results.tsv", sep="\t", index=False)
    (RES / "leftover_2026_analysis_summary.json").write_text(
        json.dumps(RESULTS, indent=2))
    print("\nWrote", RES / "leftover_2026_analysis_results.tsv")
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
