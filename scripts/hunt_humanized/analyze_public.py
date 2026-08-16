#!/usr/bin/env python3
"""Reanalyze public HIS-NSCLC bulk RNA for TACSTD2/CLDN4 vs CD8 transcripts.

This does NOT measure nest entry. Dissociated bulk RNA has no spatial axis.
"""
from __future__ import annotations

import csv
import gzip
import math
import os
import zipfile
from collections import defaultdict
from xml.etree import ElementTree as ET

import matplotlib.pyplot as plt
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
OUT = os.path.join(ROOT, "results", "hunt_humanized")
FIG = os.path.join(OUT, "figures")
DATA = os.environ.get("HUNT_DATA", "/tmp/geo_data")
os.makedirs(FIG, exist_ok=True)

GENES = [
    "TACSTD2",
    "CLDN4",
    "CD8A",
    "CD8B",
    "CD3D",
    "CD3E",
    "CD3G",
    "PTPRC",
    "EPCAM",
    "KRT8",
    "KRT18",
    "GZMB",
    "PRF1",
    "ITGAE",
    "PDCD1",
    "LAG3",
    "CD4",
    "FOXP3",
    "CCL20",
    "CXCL9",
    "CXCL10",
    "IFNG",
    "CD274",
]


def ranks(xs):
    n = len(xs)
    order = sorted(range(n), key=lambda i: xs[i])
    out = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        avg = (i + j) / 2 + 1
        for k in range(i, j + 1):
            out[order[k]] = avg
        i = j + 1
    return out


def pearson(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    if a.std() == 0 or b.std() == 0:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def spearman(a, b):
    return pearson(ranks(list(a)), ranks(list(b)))


def spearman_p(rho, n):
    if n < 3 or not math.isfinite(rho) or abs(rho) >= 1:
        return float("nan")
    t = rho * math.sqrt((n - 2) / (1 - rho * rho))
    # two-sided Student-t via regularized incomplete beta
    x = (n - 2) / ((n - 2) + t * t)
    # I_x(a,b) with a=b=(n-2)/2 is symmetric; use continued fraction-ish numpy-free approx
    # Use scipy-free Wilson-Hilferty-ish: fall back to permutation if needed.
    return t_sf_two_sided(abs(t), n - 2)


def t_sf_two_sided(t, df):
    # Regularized incomplete beta I_x(df/2, 1/2) = 2*P(T>t) for two-sided? 
    # P(|T|>t) = I_{df/(df+t^2)}(df/2, 1/2)
    x = df / (df + t * t)
    a = df / 2.0
    b = 0.5
    return incbeta_reg(x, a, b)


def incbeta_reg(x, a, b):
    # Lentz continued fraction for regularized incomplete beta
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    # Use symmetry
    if x > (a + 1) / (a + b + 2):
        return 1 - incbeta_reg(1 - x, b, a)
    ln_beta = math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)
    front = math.exp(math.log(x) * a + math.log(1 - x) * b - ln_beta) / a
    # continued fraction (Numerical Recipes)
    maxiter = 200
    eps = 3e-12
    qab = a + b
    qap = a + 1
    qam = a - 1
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < 1e-30:
        d = 1e-30
    d = 1.0 / d
    h = d
    for m in range(1, maxiter + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < 1e-30:
            d = 1e-30
        c = 1.0 + aa / c
        if abs(c) < 1e-30:
            c = 1e-30
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < 1e-30:
            d = 1e-30
        c = 1.0 + aa / c
        if abs(c) < 1e-30:
            c = 1e-30
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < eps:
            break
    return front * h


def residualize(y, x):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    X = np.column_stack([np.ones(len(x)), x])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return y - X @ beta


def write_csv(path, header, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)


def analyze_gse260575():
    d = os.path.join(DATA, "gse260575")
    files = sorted(
        fn
        for fn in os.listdir(d)
        if fn.endswith(".txt.gz")
    )
    samples = []
    counts = {}
    lib = {}
    for fn in files:
        sample = fn.replace(".txt.gz", "").split("_", 1)[1]
        path = os.path.join(d, fn)
        gene_c = {}
        tot = 0.0
        with gzip.open(path, "rt") as fh:
            next(fh)
            for line in fh:
                g, c = line.rstrip("\n").split("\t")
                c = float(c)
                tot += c
                if g in GENES:
                    gene_c[g] = c
        counts[sample] = gene_c
        lib[sample] = tot
        if sample.startswith("Con"):
            grp = "Control"
        elif sample.startswith("Pem"):
            grp = "Pembro"
        elif sample.startswith("BJ"):
            grp = "BJIKT"
        elif sample.startswith("Combi"):
            grp = "Combo"
        else:
            grp = "?"
        samples.append((sample, grp))

    sample_names = [s for s, _ in samples]
    header = ["gene"] + [f"{s}|{g}|count" for s, g in samples] + [f"{s}|{g}|cpm" for s, g in samples]
    rows = []
    cpm = {}
    for g in GENES:
        row = [g]
        cpm[g] = []
        for s, _ in samples:
            row.append(counts[s].get(g, 0.0))
        for s, _ in samples:
            v = 1e6 * counts[s].get(g, 0.0) / lib[s]
            cpm[g].append(v)
            row.append(v)
        rows.append(row)
    write_csv(os.path.join(OUT, "gse260575_gene_table.csv"), header, rows)

    lib_rows = [[s, g, lib[s]] for s, g in samples]
    write_csv(os.path.join(OUT, "gse260575_library_sizes.csv"), ["sample", "group", "total_counts"], lib_rows)

    cor_rows = []
    pairs = [
        ("TACSTD2", "CD8A"),
        ("TACSTD2", "CD8B"),
        ("TACSTD2", "CD3D"),
        ("CLDN4", "CD8A"),
        ("CLDN4", "CD8B"),
        ("CLDN4", "CD3D"),
        ("CLDN4", "PTPRC"),
        ("KRT8", "CD8A"),
        ("EPCAM", "CD8A"),
    ]
    n = len(sample_names)
    for g1, g2 in pairs:
        rho = spearman(cpm[g1], cpm[g2])
        cor_rows.append(
            [
                "GSE260575",
                "all_n12",
                g1,
                g2,
                f"{rho:.4f}" if math.isfinite(rho) else "NA",
                f"{spearman_p(rho, n):.4g}" if math.isfinite(rho) else "NA",
                sum(v == 0 for v in cpm[g1]),
                sum(v == 0 for v in cpm[g2]),
                "bulk_cpm_not_nest",
            ]
        )
    return samples, cpm, cor_rows


def _xlsx_shared_strings(z):
    root = ET.fromstring(z.read("xl/sharedStrings.xml"))
    ns = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    ss = []
    for si in root.findall("m:si", ns):
        ss.append("".join(t.text or "" for t in si.findall(".//m:t", ns)))
    return ss, ns


def _colrow(ref):
    col = ""
    row = ""
    for ch in ref:
        if ch.isalpha():
            col += ch
        else:
            row += ch
    n = 0
    for ch in col:
        n = n * 26 + (ord(ch) - 64)
    return n, int(row)


def analyze_gse276724():
    path = os.path.join(DATA, "GSE276724_Bulk-RNA-Seq_PDX.xlsx")
    z = zipfile.ZipFile(path)
    ss, ns = _xlsx_shared_strings(z)
    root = ET.fromstring(z.read("xl/worksheets/sheet1.xml"))

    def cell_val(c):
        t = c.attrib.get("t")
        v = c.find("m:v", ns)
        if v is None or v.text is None:
            return None
        if t == "s":
            return ss[int(v.text)]
        try:
            return float(v.text)
        except ValueError:
            return v.text

    headers = {}
    found = {}
    want = set(GENES)
    for c in root.findall(".//m:c", ns):
        ref = c.attrib.get("r")
        if not ref:
            continue
        col, row = _colrow(ref)
        if row <= 2:
            headers.setdefault(row, {})[col] = cell_val(c)
        elif col == 3:
            g = cell_val(c)
            if isinstance(g, str) and g in want:
                found[g] = row

    models = []
    cur = None
    for col in range(4, 19):
        m = headers.get(1, {}).get(col)
        if m:
            cur = str(m).replace("YHIM-", "YHIM").replace("YHIM", "YHIM-")
            if cur == "YHIM--1062":
                cur = "YHIM-1062"
            cur = cur.replace("YHIM-1018", "YHIM-1018").replace("YHIM-2018", "YHIM-2018")
            cur = cur.replace("YHIM-2004", "YHIM-2004").replace("YHIM-1013", "YHIM-1013")
            if not cur.startswith("YHIM-"):
                cur = "YHIM-" + cur.replace("YHIM", "")
        sid = headers.get(2, {}).get(col)
        models.append((col, cur, sid))

    # normalize model names from first-row merged labels
    raw_models = []
    cur = None
    for col in range(4, 19):
        m = headers.get(1, {}).get(col)
        if m:
            cur = str(m).replace(" ", "")
        raw_models.append(cur)
    # canonicalize
    canon = {
        "YHIM-1062": "YHIM-1062",
        "YHIM1062": "YHIM-1062",
        "YHIM1018": "YHIM-1018",
        "YHIM-1018": "YHIM-1018",
        "YHIM2018": "YHIM-2018",
        "YHIM-2018": "YHIM-2018",
        "YHIM2004": "YHIM-2004",
        "YHIM-2004": "YHIM-2004",
        "YHIM1013": "YHIM-1013",
        "YHIM-1013": "YHIM-1013",
    }
    models = []
    for i, col in enumerate(range(4, 19)):
        models.append((col, canon[raw_models[i]], headers.get(2, {}).get(col)))

    gmat = defaultdict(dict)
    lib = defaultdict(float)
    for c in root.findall(".//m:c", ns):
        ref = c.attrib.get("r")
        if not ref:
            continue
        col, row = _colrow(ref)
        if row >= 3 and col >= 4:
            v = cell_val(c)
            if isinstance(v, (int, float)):
                lib[col] += float(v)
        if row in found.values():
            gmat[row][col] = cell_val(c)

    counts = {g: [float(gmat[found[g]].get(col, 0) or 0) for col, _, _ in models] for g in GENES if g in found}
    cpm = {
        g: [1e6 * counts[g][i] / lib[col] for i, (col, _, _) in enumerate(models)]
        for g in counts
    }

    header = ["gene"] + [f"{m}|{s}|count" for _, m, s in models] + [f"{m}|{s}|cpm" for _, m, s in models]
    rows = []
    for g in GENES:
        if g not in counts:
            continue
        rows.append([g] + counts[g] + cpm[g])
    write_csv(os.path.join(OUT, "gse276724_gene_table.csv"), header, rows)
    write_csv(
        os.path.join(OUT, "gse276724_samples.csv"),
        ["col", "pdx", "geo_sample_id", "library_size", "note"],
        [
            [
                col,
                m,
                s,
                lib[col],
                "possible_epithelial_dropout" if cpm["KRT8"][i] < 20 else "",
            ]
            for i, (col, m, s) in enumerate(models)
        ],
    )

    def cor_block(mask, label, extra_note=""):
        idx = [i for i, keep in enumerate(mask) if keep]
        n = len(idx)
        out = []
        pairs = [
            ("TACSTD2", "CD8A"),
            ("TACSTD2", "CD8B"),
            ("TACSTD2", "CD3D"),
            ("TACSTD2", "PTPRC"),
            ("CLDN4", "CD8A"),
            ("CLDN4", "CD8B"),
            ("CLDN4", "CD3D"),
            ("CLDN4", "PTPRC"),
            ("EPCAM", "CD8A"),
            ("KRT8", "CD8A"),
        ]
        for g1, g2 in pairs:
            a = [cpm[g1][i] for i in idx]
            b = [cpm[g2][i] for i in idx]
            rho = spearman(a, b)
            out.append(
                [
                    "GSE276724",
                    label,
                    g1,
                    g2,
                    f"{rho:.4f}" if math.isfinite(rho) else "NA",
                    f"{spearman_p(rho, n):.4g}" if math.isfinite(rho) else "NA",
                    sum(v == 0 for v in a),
                    sum(v == 0 for v in b),
                    "bulk_cpm_not_nest" + extra_note,
                ]
            )
        # residualize on log1p EPCAM
        ep = [math.log1p(cpm["EPCAM"][i]) for i in idx]
        for g1, g2 in [("TACSTD2", "CD8A"), ("CLDN4", "CD8A"), ("TACSTD2", "CD8B"), ("CLDN4", "CD8B")]:
            a = residualize([math.log1p(cpm[g1][i]) for i in idx], ep)
            b = residualize([math.log1p(cpm[g2][i]) for i in idx], ep)
            rho = spearman(a, b)
            out.append(
                [
                    "GSE276724",
                    label + "_resid_log1pEPCAM",
                    g1,
                    g2,
                    f"{rho:.4f}" if math.isfinite(rho) else "NA",
                    f"{spearman_p(rho, n):.4g}" if math.isfinite(rho) else "NA",
                    0,
                    0,
                    "purity_residual_still_not_nest",
                ]
            )
        return out

    mask_all = [True] * 15
    # YHIM-2004 sample SJ1RNAMFF004644 = col 14, index 10
    mask_drop = [models[i][2] != "SJ1RNAMFF004644" for i in range(15)]
    cor = cor_block(mask_all, "all_n15") + cor_block(mask_drop, "drop_epithelial_dropout_n14")
    return models, cpm, cor


def plot_figures(s260, cpm260, models276, cpm276):
    plt.rcParams.update({"font.size": 9, "figure.dpi": 140})

    # Fig 1: GSE260575 CPM bars
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.4))
    groups = ["Control", "Pembro", "BJIKT", "Combo"]
    grp_idx = {g: [i for i, (_, gg) in enumerate(s260) if gg == g] for g in groups}
    colors = ["#4C4C4C", "#2A6F97", "#B08968", "#9B2226"]
    for ax, gene in zip(axes, ["TACSTD2", "CLDN4", "CD8A"]):
        means = [np.mean([cpm260[gene][i] for i in grp_idx[g]]) for g in groups]
        sds = [np.std([cpm260[gene][i] for i in grp_idx[g]], ddof=1) for g in groups]
        ax.bar(groups, means, yerr=sds, color=colors, capsize=3, width=0.7)
        for g, col in zip(groups, colors):
            xs = np.linspace(-0.15, 0.15, len(grp_idx[g]))
            ax.scatter(
                [groups.index(g) + x for x in xs],
                [cpm260[gene][i] for i in grp_idx[g]],
                color="white",
                edgecolor="black",
                s=18,
                zorder=3,
            )
        ax.set_title(f"GSE260575 {gene} CPM")
        ax.set_ylabel("CPM")
        ax.tick_params(axis="x", rotation=20)
        if gene == "TACSTD2":
            ax.set_ylim(-0.05, 1.0)
            ax.text(0.5, 0.7, "all zeros", ha="center", transform=ax.transAxes, color="#9B2226")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig1_gse260575_cpm.png"))
    plt.close()

    # Fig 2: GSE276724 scatter
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.0))
    pdxs = sorted(set(m for _, m, _ in models276))
    cmap = {
        "YHIM-1062": "#264653",
        "YHIM-1018": "#2A9D8F",
        "YHIM-2018": "#E9C46A",
        "YHIM-2004": "#E76F51",
        "YHIM-1013": "#6D597A",
    }
    outlier = [s == "SJ1RNAMFF004644" for _, _, s in models276]
    for ax, xg, yg in zip(axes, ["TACSTD2", "CLDN4"], ["CD8A", "CD8A"]):
        for i, (_, m, s) in enumerate(models276):
            ax.scatter(
                cpm276[xg][i],
                cpm276[yg][i],
                c=cmap[m],
                s=70 if not outlier[i] else 110,
                marker="o" if not outlier[i] else "X",
                edgecolor="black",
                linewidth=0.5,
                label=m if (not outlier[i] and i == [j for j, t in enumerate(models276) if t[1] == m][0]) else None,
            )
        ax.set_xlabel(f"{xg} CPM")
        ax.set_ylabel(f"{yg} CPM")
        ax.set_title(f"GSE276724 {xg} vs {yg}\n(X = epithelial dropout)")
        ax.set_xscale("symlog", linthresh=10)
        ax.set_yscale("symlog", linthresh=1)
    handles = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=cmap[p], markeredgecolor="k", markersize=8, label=p)
        for p in pdxs
    ]
    axes[1].legend(handles=handles, fontsize=7, loc="upper right")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig2_gse276724_scatter.png"))
    plt.close()

    # Fig 3: honesty schematic
    fig, ax = plt.subplots(figsize=(8.8, 3.6))
    ax.axis("off")
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 4)
    boxes = [
        (0.2, 2.2, 3.0, 1.5, "#F4F1DE", "Question\nTACSTD2/CLDN4 vs\nCD8 nest entry"),
        (3.6, 2.2, 3.0, 1.5, "#E6E6E6", "Public HIS NSCLC\nRNA / scRNA"),
        (7.0, 2.2, 2.8, 1.5, "#E07A5F", "Nest entry\nNOT measurable"),
        (0.2, 0.25, 3.0, 1.6, "#81B29A", "GSE260575 PBMC\nTACSTD2 = 0\nCD8 ~ 0"),
        (3.6, 0.25, 3.0, 1.6, "#81B29A", "GSE276724 CD34 PDX\nbulk, no treatment\nlabels; purity-bound"),
        (7.0, 0.25, 2.8, 1.6, "#F2CC8F", "GSE293914 scRNA\nno annotation\nno spatial"),
    ]
    for x, y, w, h, c, t in boxes:
        ax.add_patch(plt.Rectangle((x, y), w, h, facecolor=c, edgecolor="#333", lw=1))
        ax.text(x + w / 2, y + h / 2, t, ha="center", va="center", fontsize=8)
    ax.annotate("", xy=(3.55, 2.95), xytext=(3.25, 2.95), arrowprops=dict(arrowstyle="->", color="#333"))
    ax.annotate("", xy=(6.95, 2.95), xytext=(6.65, 2.95), arrowprops=dict(arrowstyle="->", color="#333"))
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig0_verdict.png"))
    plt.close()


def main():
    s260, cpm260, cor260 = analyze_gse260575()
    models276, cpm276, cor276 = analyze_gse276724()
    write_csv(
        os.path.join(OUT, "correlations.csv"),
        ["dataset", "subset", "gene_x", "gene_y", "spearman_rho", "approx_p", "x_zeros", "y_zeros", "caveat"],
        cor260 + cor276,
    )
    plot_figures(s260, cpm260, models276, cpm276)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
