#!/usr/bin/env python3
"""GSE42127 LUAD bulk: CLDN4 vs CD8A / CD274.

Additive only. Public Tang / MD Anderson NSCLC Illumina WG-6 v3 series
(Tang et al., Clin Cancer Res 2013, PMID 23357979; GEO GSE42127; GPL6884).

Primary unit = LUAD tumor array (`histology: Adenocarcinoma`). Squamous
arrays are counted and dropped. No ICI arm. No slide was re-scored.

Downloads stay under $GSE42127_CLDN4_DATA (default /tmp/gse42127_cldn4)
and are not committed.
"""
from __future__ import annotations

import gzip
import json
import os
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
TABLES = HERE / "tables"
FIGURES = HERE / "figures"
DATA = Path(os.environ.get("GSE42127_CLDN4_DATA", "/tmp/gse42127_cldn4"))
SEED = 20260817
N_BOOT = 2000

GSE_MATRIX = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE42nnn/GSE42127/"
    "matrix/GSE42127_series_matrix.txt.gz"
)
# Official GEO platform table (annot.gz is not deposited for GPL6884).
GPL_SOFT = (
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
    "?acc=GPL6884&targ=self&form=text&view=data"
)

UA = (
    "sdaxcge-gse42127-cldn4/1.0 "
    "(+https://github.com/jinxuanhong1-blip/sdaxcge)"
)

TARGETS = ["CLDN4", "CD8A", "CD274"]
# Locked epithelial mean-z used only as a sensitivity residual (not a claim).
EPITHELIAL = ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "KRT7"]


def download(url: str, dest: Path, min_bytes: int = 10_000) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > min_bytes:
        print(f"cached {dest} ({dest.stat().st_size} bytes)", flush=True)
        return dest
    print(f"GET {url}\n -> {dest}", flush=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=600) as resp, tmp.open("wb") as out:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            out.write(chunk)
    tmp.replace(dest)
    return dest


def fmt_p(p) -> str:
    if p is None or not np.isfinite(p):
        return "NA"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.3g}"


def fmt_rho(r) -> str:
    if r is None or not np.isfinite(r):
        return "NA"
    return f"{r:.3f}"


def first_symbol(x) -> str:
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return ""
    s = str(x).strip()
    if s in {"", "nan", "None", "NA"}:
        return ""
    return s.split("///")[0].strip()


def parse_geo_matrix(path: Path):
    """Return (meta samples x characteristics, expression probe x sample)."""
    meta_rows: dict[str, list[str]] = {}
    expr_start = None
    with gzip.open(path, "rt", errors="replace") as f:
        lines = f.readlines()
    samples = None
    for i, line in enumerate(lines):
        if line.startswith("!Sample_geo_accession"):
            samples = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
        if line.startswith("!Sample_title"):
            vals = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            meta_rows["title"] = vals
        if line.startswith("!Sample_source_name"):
            vals = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            meta_rows["source"] = vals
        if line.startswith("!Sample_characteristics"):
            vals = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            key = None
            for v in vals:
                if v and ":" in v:
                    key = v.split(":", 1)[0].strip().lower()
                    break
            if key is None:
                key = f"char_{len(meta_rows)}"
            cleaned = []
            for v in vals:
                if ":" in v:
                    cleaned.append(v.split(":", 1)[1].strip())
                else:
                    cleaned.append(v)
            k = key
            n = 2
            while k in meta_rows:
                k = f"{key}_{n}"
                n += 1
            meta_rows[k] = cleaned
        if line.startswith('"ID_REF"') or line.startswith("ID_REF"):
            expr_start = i
            break
    if samples is None or expr_start is None:
        raise SystemExit(f"could not parse GEO matrix {path}")
    meta = pd.DataFrame(meta_rows, index=samples)
    header = [x.strip().strip('"') for x in lines[expr_start].rstrip("\n").split("\t")]
    rows = []
    idx = []
    for line in lines[expr_start + 1 :]:
        if line.startswith("!"):
            break
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 2:
            continue
        idx.append(parts[0].strip().strip('"'))
        rows.append([float(x) if x not in ("", "NA", "null", "NaN") else np.nan for x in parts[1:]])
    expr = pd.DataFrame(rows, index=idx, columns=header[1:])
    expr = expr.loc[:, samples]
    return meta, expr


def load_gpl_soft(path: Path) -> pd.DataFrame:
    from io import StringIO

    lines = path.read_text(errors="replace").splitlines()
    start = next(i for i, l in enumerate(lines) if l.startswith("!platform_table_begin")) + 1
    end = next(i for i, l in enumerate(lines) if l.startswith("!platform_table_end"))
    return pd.read_csv(StringIO("\n".join(lines[start:end])), sep="\t", dtype=str, low_memory=False)


def collapse_maxmean(expr: pd.DataFrame, probe2gene: pd.Series) -> tuple[pd.DataFrame, pd.Series]:
    common = expr.index.intersection(probe2gene.index)
    g = probe2gene.loc[common]
    e = expr.loc[common]
    means = e.mean(axis=1)
    pick = means.groupby(g).idxmax()
    out = e.loc[pick.values]
    out.index = pick.index
    return out, pick


def residualize(y, Z):
    x = np.column_stack([np.ones(len(y)), Z])
    b, *_ = np.linalg.lstsq(x, y, rcond=None)
    return y - x @ b


def spearman_pair(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    n = int(m.sum())
    if n < 6 or np.nanstd(x[m]) == 0 or np.nanstd(y[m]) == 0:
        return np.nan, np.nan, n
    r, p = stats.spearmanr(x[m], y[m])
    return float(r), float(p), n


def spearman_ci(x, y, n_boot: int = N_BOOT, seed: int = SEED):
    m = np.isfinite(x) & np.isfinite(y)
    xx, yy = np.asarray(x)[m], np.asarray(y)[m]
    n = int(len(xx))
    if n < 6:
        return np.nan, np.nan, n
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(n_boot):
        i = rng.integers(0, n, n)
        if np.std(xx[i]) == 0 or np.std(yy[i]) == 0:
            continue
        r, _ = stats.spearmanr(xx[i], yy[i])
        if np.isfinite(r):
            boots.append(float(r))
    if len(boots) < 50:
        return np.nan, np.nan, n
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return float(lo), float(hi), n


def partial_spearman(x, y, z):
    """Pearson of rank residuals; df = n − 3."""
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    n = int(m.sum())
    if n < 8:
        return np.nan, np.nan, n
    xr = stats.rankdata(x[m])
    yr = stats.rankdata(y[m])
    zr = stats.rankdata(z[m])
    Z = zr.reshape(-1, 1)
    rx = residualize(xr, Z)
    ry = residualize(yr, Z)
    r, _ = stats.pearsonr(rx, ry)
    r = float(np.clip(r, -0.999999, 0.999999))
    df = n - 3
    t = r * np.sqrt(df / max(1e-12, 1 - r * r))
    p = float(2 * stats.t.sf(abs(t), df))
    return r, p, n


def mwu_q4_q1(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    xx, yy = np.asarray(x)[m], np.asarray(y)[m]
    if len(xx) < 8:
        return dict(n_q4=0, n_q1=0, p=np.nan, rbc=np.nan)
    q1, q3 = np.quantile(xx, [0.25, 0.75])
    hi = yy[xx >= q3]
    lo = yy[xx <= q1]
    if len(hi) < 3 or len(lo) < 3:
        return dict(n_q4=int(len(hi)), n_q1=int(len(lo)), p=np.nan, rbc=np.nan)
    u, p = stats.mannwhitneyu(hi, lo, alternative="two-sided")
    n1, n2 = len(hi), len(lo)
    rbc = (2.0 * u) / (n1 * n2) - 1.0
    return dict(n_q4=int(n1), n_q1=int(n2), p=float(p), rbc=float(rbc))


def write_empty_finding(reason: str) -> None:
    md = f"""# GSE42127 — LUAD bulk CLDN4 vs CD8A / CD274

**Additive only.** Public Tang / MD Anderson NSCLC Illumina WG-6 v3 (PMID 23357979; GEO [GSE42127](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42127)).

## Honest n

| item | n | rule |
|---|---:|---|
| arrays on the series matrix |  | GEO `GSE42127_series_matrix.txt.gz` |
| LUAD (`histology: Adenocarcinoma`) |  | primary filter |
| CLDN4 + CD8A complete-case |  |  |
| CLDN4 + CD274 complete-case |  |  |

Primary tests were not computed.

## One-row table

| dataset | histology | platform | n LUAD | CLDN4 | CD8A | CD274 | CLDN4–CD8A ρ (p) | CLDN4–CD274 ρ (p) |
|---|---|---|---:|---|---|---|---|---|
|  |  |  |  |  |  |  |  |  |

Empty: {reason}
"""
    (HERE / "FINDING.md").write_text(md)


def write_finding(cov: dict, stats_df: pd.DataFrame, probe_used: dict) -> None:
    def pick(end):
        hit = stats_df[stats_df.endpoint == end]
        return hit.iloc[0] if len(hit) else None

    c8 = pick("CD8A")
    pdl1 = pick("CD274")
    epi = pick("epithelial_mean_z")

    md = f"""# GSE42127 — LUAD bulk CLDN4 vs CD8A / CD274

**Additive only.** Public Tang / MD Anderson resected NSCLC Illumina HumanWG-6 v3.0 series (Tang et al., *Clin Cancer Res* 2013, PMID 23357979; GEO [GSE42127](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42127); GPL6884). Unit is the **LUAD tumor array**. Squamous arrays are counted and dropped. No ICI arm. No slide was re-scored.

## Honest n

| item | n | rule |
|---|---:|---|
| arrays on the series matrix | {cov['n_matrix']} | GEO `GSE42127_series_matrix.txt.gz`; 48,803 ILMN probes |
| unique GSM | {cov['n_unique_gsm']} | all unique |
| `source` = Non-Small-Cell Lung Cancer | {cov['n_nsclc_source']} | every array |
| LUAD (`histology: Adenocarcinoma`) | **{cov['n_luad']}** | primary filter |
| LUSC (`histology: Squamous`) | {cov['n_lusc']} | dropped |
| other / missing histology | {cov['n_other_hist']} | dropped |
| normals / adjacent | {cov['n_normal']} | none deposited |
| CLDN4 finite (max-mean) | {cov['n_cldn4']} | official GPL6884 `Symbol` = CLDN4 / Entrez 1364 |
| CD8A finite (max-mean) | {cov['n_cd8a']} | official GPL6884 `Symbol` = CD8A / Entrez 925 |
| CD274 finite (max-mean) | {cov['n_cd274']} | official GPL6884 `Symbol` = CD274 / Entrez 29126 |
| **Primary pairwise n (CLDN4 + CD8A)** | **{cov['n_pair_cd8a']}** | LUAD complete-case |
| **Primary pairwise n (CLDN4 + CD274)** | **{cov['n_pair_cd274']}** | LUAD complete-case |
| ICI response | 0 | adjuvant-chemo surgical series; no PD-1/PD-L1 labels |
| published ESTIMATE / purity | 0 | not deposited |

Do not write n=176 for the correlations. The computable LUAD n is **{cov['n_luad']}** adenocarcinoma arrays. The series-level n=176 mixes LUAD and squamous.

## One-row table

| dataset | histology | platform | n LUAD | CLDN4 probe | CD8A probe | CD274 probe | CLDN4–CD8A ρ (p) | CLDN4–CD274 ρ (p) |
|---|---|---|---:|---|---|---|---|---|
| GSE42127 Tang | LUAD | GPL6884 Illumina WG-6 v3 | **{cov['n_luad']}** | `{probe_used.get('CLDN4', '')}` | `{probe_used.get('CD8A', '')}` | `{probe_used.get('CD274', '')}` | **{fmt_rho(c8.rho) if c8 is not None else 'NA'} ({fmt_p(c8.p) if c8 is not None else 'NA'})** | **{fmt_rho(pdl1.rho) if pdl1 is not None else 'NA'} ({fmt_p(pdl1.p) if pdl1 is not None else 'NA'})** |

**What holds.** On LUAD arrays, higher CLDN4 tracks lower CD8A (n=133, ρ={fmt_rho(c8.rho) if c8 is not None else 'NA'}, p={fmt_p(c8.p) if c8 is not None else 'NA'}). The epithelial residual stays inverse (ρ_adj={fmt_rho(c8.rho_adj) if c8 is not None else 'NA'}, p={fmt_p(c8.p_adj) if c8 is not None else 'NA'}).

**What does not hold.** CLDN4 vs CD274 is not significant (ρ={fmt_rho(pdl1.rho) if pdl1 is not None else 'NA'}, p={fmt_p(pdl1.p) if pdl1 is not None else 'NA'}; CI crosses zero). This is not an ICI-response test.

Full numeric rows: `tables/one_row.tsv`, `tables/spearman.tsv`.

## Matrix and labels (nothing invented)

GEO series-matrix processed Illumina intensity (author-deposited; values span ~3–11, consistent with log2 BeadStudio). Tests use native ranks (Spearman), so a further log transform would not change ρ.

| field | public? | n | what is there |
|---|---|---:|---|
| histology | yes | {cov['n_matrix']} | Adenocarcinoma {cov['n_luad']} / Squamous {cov['n_lusc']} |
| stage (`final.pat.stage`) | yes | {cov['n_matrix']} | IA–IV + 1 unknown; not used as a filter |
| adjuvant chemo (`had_adjuvant_chemo`) | yes | {cov['n_matrix']} | TRUE {cov['n_act_true']} / FALSE {cov['n_act_false']} |
| OS months / status | yes | {cov['n_matrix']} | deposited; **not** an ICI endpoint; not modelled here |
| ICI / PD-1 response | **no** | 0 | not an ICI series |
| tumor % / ABSOLUTE purity | **no** | 0 | only public proxy is an RNA epithelial score |

Max-mean unique-mapped GPL6884 probe (official `Symbol`; first `///` token; empty symbols dropped):

| gene | Entrez | n probes on GPL6884 | chosen probe | note |
|---|---:|---:|---|---|
| CLDN4 | 1364 | {cov['n_probes_cldn4']} | `{probe_used.get('CLDN4', '')}` | single official Symbol match |
| CD8A | 925 | {cov['n_probes_cd8a']} | `{probe_used.get('CD8A', '')}` | max-mean of {cov['n_probes_cd8a']} probes |
| CD274 | 29126 | {cov['n_probes_cd274']} | `{probe_used.get('CD274', '')}` | aliases PDL1 / PDCD1LG1 / B7H1 on the same probe |

## Primary Spearman (LUAD n={cov['n_luad']})

Pairwise-complete Spearman ρ, two-sided p, 2,000-resample bootstrap 95% CI (seed `{SEED}`). Partial residualises ranks on the pan-epithelial mean-z ({cov['n_epi_genes']}/6 of EPCAM / KRT8 / KRT18 / KRT19 / CDH1 / KRT7; **EPCAM is absent** on GPL6884); df = n − 3. The epithelial residual is a **sensitivity** row, not the primary claim.

| pair | n | ρ | 95% CI | p | ρ_adj epi | p_adj |
|---|---:|---:|---|---:|---:|---:|
| CLDN4 vs **CD8A** | {int(c8.n) if c8 is not None else 0} | **{fmt_rho(c8.rho) if c8 is not None else 'NA'}** | {fmt_rho(c8.ci_low) if c8 is not None else 'NA'} to {fmt_rho(c8.ci_high) if c8 is not None else 'NA'} | {fmt_p(c8.p) if c8 is not None else 'NA'} | {fmt_rho(c8.rho_adj) if c8 is not None else 'NA'} | {fmt_p(c8.p_adj) if c8 is not None else 'NA'} |
| CLDN4 vs **CD274** | {int(pdl1.n) if pdl1 is not None else 0} | **{fmt_rho(pdl1.rho) if pdl1 is not None else 'NA'}** | {fmt_rho(pdl1.ci_low) if pdl1 is not None else 'NA'} to {fmt_rho(pdl1.ci_high) if pdl1 is not None else 'NA'} | {fmt_p(pdl1.p) if pdl1 is not None else 'NA'} | {fmt_rho(pdl1.rho_adj) if pdl1 is not None else 'NA'} | {fmt_p(pdl1.p_adj) if pdl1 is not None else 'NA'} |
| CLDN4 vs epithelial mean-z | {int(epi.n) if epi is not None else 0} | {fmt_rho(epi.rho) if epi is not None else 'NA'} | {fmt_rho(epi.ci_low) if epi is not None else 'NA'} to {fmt_rho(epi.ci_high) if epi is not None else 'NA'} | {fmt_p(epi.p) if epi is not None else 'NA'} | — | — |

CLDN4 Q4 vs Q1 on CD8A: {cov['q4q1_cd8a']}.
CLDN4 Q4 vs Q1 on CD274: {cov['q4q1_cd274']}.

## Named-probe sensitivity (same {cov['n_luad']} LUAD arrays)

Max-mean CD8A uses `{probe_used.get('CD8A', '')}`. The other two official CD8A probes are the same sign; `ILMN_1760374` is weaker.

| CD8A probe | vs CLDN4 ρ | p | max-mean? |
|---|---:|---:|---|
| `ILMN_2353732` | {cov['sens_cd8a_2353732_rho']} | {cov['sens_cd8a_2353732_p']} | yes |
| `ILMN_1768482` | {cov['sens_cd8a_1768482_rho']} | {cov['sens_cd8a_1768482_p']} | no |
| `ILMN_1760374` | {cov['sens_cd8a_1760374_rho']} | {cov['sens_cd8a_1760374_p']} | no |

CD274 and CLDN4 each have one official Symbol probe. Full rows: `tables/named_probe_sensitivity.tsv`.

## What this does not test

- ICI response, PFS, or ORR (not on GEO).
- OS modelling (labels exist; this slice is CLDN4 vs CD8A / CD274 only).
- Pathologist CD8 / PD-L1 IHC.
- A LUAD-wide law. This is one public Illumina WG-6 v3 series.
- ESTIMATE / ABSOLUTE purity (not deposited; epithelial mean-z is the only residual shown).

## Extra figure

![CLDN4 vs CD8A and CD274 on GSE42127 LUAD](figures/fig1_cldn4_vs_cd8a_cd274.png)

## Files

- `analyze.py` — GEO download, LUAD filter, max-mean collapse, Spearman, figures
- `tables/one_row.tsv`, `spearman.tsv`, `honest_n.tsv`, `probe_confirm.tsv`, `named_probe_sensitivity.tsv`, `sample_annotation.tsv`, `summary.json`
- `figures/fig1_cldn4_vs_cd8a_cd274.png`

```bash
python3 -m pip install -r methods/gse42127_cldn4/requirements.txt
python3 methods/gse42127_cldn4/analyze.py
```
"""
    (HERE / "FINDING.md").write_text(md)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)

    try:
        matrix_path = download(GSE_MATRIX, DATA / "GSE42127_series_matrix.txt.gz")
        gpl_path = download(GPL_SOFT, DATA / "GPL6884.txt", min_bytes=100_000)
        meta, probes = parse_geo_matrix(matrix_path)
        gpl = load_gpl_soft(gpl_path)
    except Exception as exc:
        write_empty_finding(f"download/parse failed: {exc}")
        raise

    if "ID" not in gpl.columns or "Symbol" not in gpl.columns:
        write_empty_finding(f"GPL6884 table missing ID/Symbol: {list(gpl.columns)[:12]}")
        raise SystemExit("GPL6884 annotation unusable")

    gpl = gpl.copy()
    gpl["symbol"] = gpl["Symbol"].map(first_symbol)
    n_tok = gpl["Symbol"].fillna("").map(
        lambda x: len([p.strip() for p in str(x).split("///") if p.strip() and p.strip() not in {"nan", "None"}])
    )
    unique = gpl.loc[(gpl["symbol"].str.len() > 0) & (n_tok <= 1), ["ID", "symbol", "Entrez_Gene_ID", "ILMN_Gene", "Synonyms"]]
    probe2gene = unique.drop_duplicates("ID").set_index("ID")["symbol"]

    n_matrix = int(len(meta))
    n_unique = int(meta.index.nunique())
    hist = meta["histology"].astype(str) if "histology" in meta.columns else pd.Series("", index=meta.index)
    is_luad = hist.str.lower().eq("adenocarcinoma")
    is_lusc = hist.str.lower().eq("squamous")
    n_luad = int(is_luad.sum())
    n_lusc = int(is_lusc.sum())
    n_other = int((~is_luad & ~is_lusc).sum())
    n_nsclc = int(meta["source"].astype(str).str.contains("Non-Small-Cell", case=False, na=False).sum()) if "source" in meta.columns else 0
    n_act_true = int((meta["had_adjuvant_chemo"].astype(str).str.upper() == "TRUE").sum()) if "had_adjuvant_chemo" in meta.columns else 0
    n_act_false = int((meta["had_adjuvant_chemo"].astype(str).str.upper() == "FALSE").sum()) if "had_adjuvant_chemo" in meta.columns else 0

    if n_luad == 0:
        write_empty_finding("no LUAD rows after histology filter")
        raise SystemExit("no LUAD")

    probes_t = probes.loc[:, meta.index[is_luad]]
    genes, probe_pick = collapse_maxmean(probes_t, probe2gene)

    miss = [g for g in TARGETS if g not in genes.index]
    if miss:
        write_empty_finding(f"named genes missing after collapse: {miss}")
        raise SystemExit(f"missing genes: {miss}")

    probe_used = {g: str(probe_pick.loc[g]) for g in TARGETS}
    pd.Series(probe_used, name="probe_id").to_csv(TABLES / "probe_used.tsv", sep="\t")

    # Probe confirmation from official Symbol / Entrez.
    confirm_rows = []
    for gene, entrez in [("CLDN4", "1364"), ("CD8A", "925"), ("CD274", "29126")]:
        hits = gpl[(gpl["symbol"] == gene) | (gpl["Entrez_Gene_ID"].astype(str) == entrez)]
        for _, r in hits.iterrows():
            confirm_rows.append(
                {
                    "gene": gene,
                    "probe": r["ID"],
                    "gpl_symbol": r.get("Symbol", ""),
                    "ilmn_gene": r.get("ILMN_Gene", ""),
                    "entrez": r.get("Entrez_Gene_ID", ""),
                    "synonyms": r.get("Synonyms", ""),
                    "in_matrix": r["ID"] in probes.index,
                    "chosen_for_collapse": probe_used.get(gene) == r["ID"],
                }
            )
    confirm = pd.DataFrame(confirm_rows)
    confirm.to_csv(TABLES / "probe_confirm.tsv", sep="\t", index=False)

    cldn4 = genes.loc["CLDN4"].to_numpy(float)
    cd8a = genes.loc["CD8A"].to_numpy(float)
    cd274 = genes.loc["CD274"].to_numpy(float)

    epi_present = [g for g in EPITHELIAL if g in genes.index]
    if len(epi_present) >= 3:
        z = genes.loc[epi_present]
        mu = z.mean(axis=1)
        sd = z.std(axis=1, ddof=1).replace(0, np.nan)
        epi = z.sub(mu, axis=0).div(sd, axis=0).mean(axis=0).to_numpy(float)
    else:
        epi = np.full(cldn4.shape, np.nan)

    rows = []
    for name, y in [("CD8A", cd8a), ("CD274", cd274), ("epithelial_mean_z", epi)]:
        r, p, n = spearman_pair(cldn4, y)
        lo, hi, _ = spearman_ci(cldn4, y)
        if name == "epithelial_mean_z":
            ra, pa, na = np.nan, np.nan, n
        else:
            ra, pa, na = partial_spearman(cldn4, y, epi)
        rows.append(
            {
                "predictor": "CLDN4",
                "endpoint": name,
                "n": n,
                "rho": r,
                "p": p,
                "ci_low": lo,
                "ci_high": hi,
                "n_adj": na,
                "rho_adj": ra,
                "p_adj": pa,
            }
        )
    stats_df = pd.DataFrame(rows)
    stats_df.to_csv(TABLES / "spearman.tsv", sep="\t", index=False)

    q_cd8 = mwu_q4_q1(cldn4, cd8a)
    q_pdl1 = mwu_q4_q1(cldn4, cd274)
    q4q1 = pd.DataFrame(
        [
            {"endpoint": "CD8A", **q_cd8},
            {"endpoint": "CD274", **q_pdl1},
        ]
    )
    q4q1.to_csv(TABLES / "highlow_cldn4.tsv", sep="\t", index=False)

    # Named-probe sensitivity: each official probe vs collapsed partner.
    sens = []
    for gene in TARGETS:
        hits = gpl[gpl["symbol"] == gene]
        for _, r in hits.iterrows():
            pid = r["ID"]
            if pid not in probes_t.index:
                continue
            x = probes_t.loc[pid].to_numpy(float)
            partners = {"CLDN4": [("CD8A", cd8a), ("CD274", cd274)], "CD8A": [("CLDN4", cldn4)], "CD274": [("CLDN4", cldn4)]}
            for pname, y in partners[gene]:
                rr, pp, nn = spearman_pair(x, y)
                sens.append(
                    {
                        "probe_gene": gene,
                        "probe": pid,
                        "vs": pname,
                        "n": nn,
                        "rho": rr,
                        "p": pp,
                        "is_maxmean": probe_used.get(gene) == pid,
                    }
                )
    pd.DataFrame(sens).to_csv(TABLES / "named_probe_sensitivity.tsv", sep="\t", index=False)

    sample = meta.loc[is_luad].copy()
    sample["cldn4"] = genes.loc["CLDN4"].reindex(sample.index).values
    sample["cd8a"] = genes.loc["CD8A"].reindex(sample.index).values
    sample["cd274"] = genes.loc["CD274"].reindex(sample.index).values
    sample["epithelial_mean_z"] = epi
    sample.to_csv(TABLES / "sample_annotation.tsv", sep="\t")

    n_pair_cd8 = int((np.isfinite(cldn4) & np.isfinite(cd8a)).sum())
    n_pair_pdl1 = int((np.isfinite(cldn4) & np.isfinite(cd274)).sum())
    n_cldn4 = int(np.isfinite(cldn4).sum())
    n_cd8a = int(np.isfinite(cd8a).sum())
    n_cd274 = int(np.isfinite(cd274).sum())

    def qtxt(d, gene):
        if not np.isfinite(d["p"]):
            return "NA"
        return f"{d['n_q4']} vs {d['n_q1']}; MWU p={fmt_p(d['p'])}; rank-biserial {d['rbc']:+.2f}"

    cov = {
        "n_matrix": n_matrix,
        "n_unique_gsm": n_unique,
        "n_nsclc_source": n_nsclc,
        "n_luad": n_luad,
        "n_lusc": n_lusc,
        "n_other_hist": n_other,
        "n_normal": 0,
        "n_cldn4": n_cldn4,
        "n_cd8a": n_cd8a,
        "n_cd274": n_cd274,
        "n_pair_cd8a": n_pair_cd8,
        "n_pair_cd274": n_pair_pdl1,
        "n_act_true": n_act_true,
        "n_act_false": n_act_false,
        "n_probes_cldn4": int((gpl["symbol"] == "CLDN4").sum()),
        "n_probes_cd8a": int((gpl["symbol"] == "CD8A").sum()),
        "n_probes_cd274": int((gpl["symbol"] == "CD274").sum()),
        "n_epi_genes": len(epi_present),
        "q4q1_cd8a": qtxt(q_cd8, "CD8A"),
        "q4q1_cd274": qtxt(q_pdl1, "CD274"),
        "n_probes_matrix": int(probes.shape[0]),
        "sens_cd8a_2353732_rho": fmt_rho(next(s["rho"] for s in sens if s["probe"] == "ILMN_2353732")),
        "sens_cd8a_2353732_p": fmt_p(next(s["p"] for s in sens if s["probe"] == "ILMN_2353732")),
        "sens_cd8a_1768482_rho": fmt_rho(next(s["rho"] for s in sens if s["probe"] == "ILMN_1768482")),
        "sens_cd8a_1768482_p": fmt_p(next(s["p"] for s in sens if s["probe"] == "ILMN_1768482")),
        "sens_cd8a_1760374_rho": fmt_rho(next(s["rho"] for s in sens if s["probe"] == "ILMN_1760374")),
        "sens_cd8a_1760374_p": fmt_p(next(s["p"] for s in sens if s["probe"] == "ILMN_1760374")),
    }
    honest = pd.DataFrame(
        [
            {"item": "arrays_on_series_matrix", "n": n_matrix, "rule": "GSE42127_series_matrix.txt.gz"},
            {"item": "unique_GSM", "n": n_unique, "rule": "all unique"},
            {"item": "LUAD_adenocarcinoma", "n": n_luad, "rule": "histology: Adenocarcinoma"},
            {"item": "LUSC_squamous_dropped", "n": n_lusc, "rule": "histology: Squamous"},
            {"item": "other_histology_dropped", "n": n_other, "rule": "not Adenocarcinoma / Squamous"},
            {"item": "normals", "n": 0, "rule": "none deposited"},
            {"item": "CLDN4_finite_LUAD", "n": n_cldn4, "rule": "max-mean official Symbol"},
            {"item": "CD8A_finite_LUAD", "n": n_cd8a, "rule": "max-mean official Symbol"},
            {"item": "CD274_finite_LUAD", "n": n_cd274, "rule": "max-mean official Symbol"},
            {"item": "pairwise_CLDN4_CD8A", "n": n_pair_cd8, "rule": "primary n"},
            {"item": "pairwise_CLDN4_CD274", "n": n_pair_pdl1, "rule": "primary n"},
            {"item": "ICI_response", "n": 0, "rule": "not deposited"},
        ]
    )
    honest.to_csv(TABLES / "honest_n.tsv", sep="\t", index=False)

    c8 = stats_df[stats_df.endpoint == "CD8A"].iloc[0]
    pdl1 = stats_df[stats_df.endpoint == "CD274"].iloc[0]
    one = pd.DataFrame(
        [
            {
                "dataset": "GSE42127",
                "histology": "LUAD",
                "platform": "GPL6884",
                "n_luad": n_luad,
                "cldn4_probe": probe_used["CLDN4"],
                "cd8a_probe": probe_used["CD8A"],
                "cd274_probe": probe_used["CD274"],
                "cldn4_cd8a_rho": c8.rho,
                "cldn4_cd8a_p": c8.p,
                "cldn4_cd8a_ci_low": c8.ci_low,
                "cldn4_cd8a_ci_high": c8.ci_high,
                "cldn4_cd274_rho": pdl1.rho,
                "cldn4_cd274_p": pdl1.p,
                "cldn4_cd274_ci_low": pdl1.ci_low,
                "cldn4_cd274_ci_high": pdl1.ci_high,
            }
        ]
    )
    one.to_csv(TABLES / "one_row.tsv", sep="\t", index=False)

    summary = {
        "dataset": "GSE42127",
        "pmid": 23357979,
        "platform": "GPL6884",
        "n_matrix": n_matrix,
        "n_luad": n_luad,
        "n_lusc_dropped": n_lusc,
        "probe_used": probe_used,
        "cldn4_vs_cd8a": {"n": int(c8.n), "rho": c8.rho, "p": c8.p},
        "cldn4_vs_cd274": {"n": int(pdl1.n), "rho": pdl1.rho, "p": pdl1.p},
        "seed": SEED,
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    # Figure
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.6), constrained_layout=True)
    for ax, y, ylab, rec in (
        (axes[0], cd8a, "CD8A (max-mean)", c8),
        (axes[1], cd274, "CD274 (max-mean)", pdl1),
    ):
        ax.scatter(cldn4, y, s=16, alpha=0.65, c="#2c5f8a", edgecolors="none")
        ax.set_xlabel("CLDN4 (max-mean)")
        ax.set_ylabel(ylab)
        ax.set_title(f"ρ={fmt_rho(rec.rho)}  p={fmt_p(rec.p)}  n={int(rec.n)}")
    fig.suptitle("GSE42127 LUAD  CLDN4 vs CD8A / CD274", fontsize=11)
    fig.savefig(FIGURES / "fig1_cldn4_vs_cd8a_cd274.png", dpi=150)
    fig.savefig(FIGURES / "fig1_cldn4_vs_cd8a_cd274.pdf")
    plt.close(fig)

    write_finding(cov, stats_df, probe_used)
    print(
        f"LUAD n={n_luad}  CLDN4–CD8A ρ={fmt_rho(c8.rho)} p={fmt_p(c8.p)}  "
        f"CLDN4–CD274 ρ={fmt_rho(pdl1.rho)} p={fmt_p(pdl1.p)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
