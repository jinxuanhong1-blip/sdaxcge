"""Independent LUSC replication + TACSTD2×histology interaction + OS.

New public arrays/RNA (all <80 MB matrices):
  GSE4573   Raponi LUSC Affy U133A
  GSE17710  Wilkerson LUSC Agilent (has pathologist tumor_percent)
  GSE19188  Hou NSCLC Affy U133 Plus 2 — SCC / ADC split
  GSE50081  Der NSCLC Affy U133 Plus 2 — SCC / ADC split
  GSE103584 Zhou NSCLC RNA-seq — SCC / ADC split

Also: formal interaction in pooled TCGA, and OS Cox in TCGA LUAD/LUSC.
"""

from __future__ import annotations

import gzip
import io
import os
import pickle
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from opus_tls_lib import (
    FOCUS_GENES, PROC_DIR, RESULTS_DIR, SIGNATURES, bh_fdr, fmt_p,
    inverse_normal, ols_table, partial_spearman, signature_scores, spearman_ci,
)

TABLES = os.path.join(RESULTS_DIR, "tables")
GEO2 = os.path.join(os.environ.get("OPUS_TLS_DATA", "/tmp/opus_tls_data"), "geo2")
PRIMARY = ["TLS_Cabrita", "TLS_12chemokine", "TLS_imprint", "B_cell",
           "Plasma_cell", "Tfh", "T_cell_CD8", "IFNg_Ayers"]


def parse_meta(path: str) -> pd.DataFrame:
    """Per-sample dict of every 'key: value' characteristic (handles mixed rows)."""
    titles, accs = None, None
    bags: list[dict] | None = None
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
            if line.startswith("!Sample_title"):
                titles = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_geo_accession"):
                accs = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
                bags = [{} for _ in accs]
            elif line.startswith("!Sample_characteristics_ch1") and bags is not None:
                vals = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
                for i, v in enumerate(vals):
                    if i >= len(bags) or ":" not in v:
                        continue
                    k, val = v.split(":", 1)
                    bags[i][k.strip().lower().replace(" ", "_")] = val.strip()
    meta = pd.DataFrame(bags, index=pd.Index(accs, name="gsm"))
    meta["title"] = titles
    return meta


def read_matrix(path: str) -> pd.DataFrame:
    with gzip.open(path, "rt", errors="replace") as fh:
        lines = fh.read().split("\n")
    start = next(i for i, l in enumerate(lines) if l.startswith("!series_matrix_table_begin")) + 1
    end = next(i for i, l in enumerate(lines) if l.startswith("!series_matrix_table_end"))
    tbl = pd.read_csv(io.StringIO("\n".join(lines[start:end])), sep="\t", index_col=0)
    tbl.index = tbl.index.astype(str).str.strip('"')
    return tbl


def collapse(expr: pd.DataFrame) -> pd.DataFrame:
    expr = expr[~expr.index.isna() & (expr.index.astype(str) != "")]
    expr.index = expr.index.astype(str).str.split(" /// ").str[0].str.strip()
    order = expr.mean(axis=1).sort_values(ascending=False).index
    expr = expr.loc[order]
    return expr[~expr.index.duplicated(keep="first")].sort_index()


def load_gpl_symbols(path: str, symbol_col: str) -> pd.Series:
    with open(path) as fh:
        lines = fh.readlines()
    start = next(i for i, l in enumerate(lines) if l.startswith("!platform_table_begin")) + 1
    end = next(i for i, l in enumerate(lines) if l.startswith("!platform_table_end"))
    ann = pd.read_csv(io.StringIO("".join(lines[start:end])), sep="\t", dtype=str, low_memory=False)
    ann = ann.dropna(subset=["ID", symbol_col])
    return ann.drop_duplicates("ID").set_index("ID")[symbol_col]


def map_expr(tbl: pd.DataFrame, id2gene: pd.Series) -> pd.DataFrame:
    tbl = tbl.copy()
    tbl.index = tbl.index.astype(str).map(lambda x: id2gene.get(x, np.nan))
    return collapse(tbl.astype(float))


def assoc_block(name: str, expr: pd.DataFrame, cov: pd.Series | None,
                hist: str, platform: str) -> list[dict]:
    scores, used = signature_scores(expr)
    ph = pd.DataFrame({"epithelial_score": scores["Epithelial"]}, index=scores.index)
    if cov is not None:
        ph["purity_proxy"] = pd.to_numeric(cov.reindex(ph.index), errors="coerce")
    rows = []
    for g in FOCUS_GENES:
        if g not in expr.index:
            continue
        x = expr.loc[g, scores.index]
        for sig in PRIMARY:
            y = scores[sig]
            rho, p, ci, n = spearman_ci(x.values, y.values, n_boot=800)
            r_ep, p_ep, n_ep = partial_spearman(x.values, y.values, ph[["epithelial_score"]])
            rec = {"cohort": name, "histology": hist, "platform": platform, "gene": g,
                   "signature": sig, "n": n, "rho": rho, "p": p,
                   "ci_low": ci[0] if isinstance(ci, tuple) else np.nan,
                   "ci_high": ci[1] if isinstance(ci, tuple) else np.nan,
                   "rho_adj_epithelial": r_ep, "p_adj_epithelial": p_ep,
                   "n_genes_in_sig": len(used.get(sig, []))}
            if "purity_proxy" in ph.columns and ph["purity_proxy"].notna().sum() >= 20:
                r_pu, p_pu, n_pu = partial_spearman(x.values, y.values, ph[["purity_proxy"]])
                rec.update({"rho_adj_tumorpct": r_pu, "p_adj_tumorpct": p_pu, "n_adj_tumorpct": n_pu})
            n_use = n_ep if np.isfinite(n_ep) else n
            holds = bool(np.isfinite(r_ep) and r_ep < 0 and p_ep < 0.05 and n_use >= 40)
            opposite = bool(np.isfinite(r_ep) and r_ep > 0 and p_ep < 0.05 and n_use >= 40)
            rec["verdict"] = ("HOLDS" if holds else ("OPPOSITE" if opposite else
                              ("UNDERPOWERED" if n < 40 else "NO_EVIDENCE")))
            rows.append(rec)
    return rows


def build_new_cohorts() -> tuple[list[dict], dict]:
    gpl96 = load_gpl_symbols(os.path.join(GEO2, "GPL96.txt"), "Gene Symbol")
    gpl570 = load_gpl_symbols(os.path.join(GEO2, "GPL570.txt"), "Gene Symbol")
    gpl9053 = load_gpl_symbols(os.path.join(GEO2, "GPL9053.txt"), "ORF")

    rows = []
    scored = {}

    # --- GSE4573 all LUSC ---
    expr = map_expr(read_matrix(os.path.join(GEO2, "GSE4573_matrix.txt.gz")), gpl96)
    print(f"[GSE4573] genes={expr.shape[0]} n={expr.shape[1]} TACSTD2={'TACSTD2' in expr.index}")
    rows += assoc_block("GSE4573", expr, None, "LUSC", "Affy U133A")
    sc, _ = signature_scores(expr)
    scored["GSE4573"] = {"expr": expr, "scores": sc, "hist": "LUSC"}

    # --- GSE17710 all LUSC + tumor_percent ---
    expr = map_expr(read_matrix(os.path.join(GEO2, "GSE17710_matrix.txt.gz")), gpl9053)
    meta = parse_meta(os.path.join(GEO2, "GSE17710_matrix.txt.gz"))
    # matrix columns are GSM
    tp = pd.to_numeric(meta.reindex(expr.columns).get("tumor_percent"), errors="coerce")
    print(f"[GSE17710] genes={expr.shape[0]} n={expr.shape[1]} tumor_pct n={tp.notna().sum()}")
    rows += assoc_block("GSE17710", expr, tp, "LUSC", "Agilent GPL9053")
    sc, _ = signature_scores(expr)
    scored["GSE17710"] = {"expr": expr, "scores": sc, "hist": "LUSC", "tumor_percent": tp}

    # --- GSE19188 tumors, split histology ---
    expr = map_expr(read_matrix(os.path.join(GEO2, "GSE19188_matrix.txt.gz")), gpl570)
    meta = parse_meta(os.path.join(GEO2, "GSE19188_matrix.txt.gz")).reindex(expr.columns)
    tumor = meta["tissue_type"].astype(str).str.lower().eq("tumor")
    cell = meta["cell_type"].astype(str)
    print(f"[GSE19188] tumors={tumor.sum()} ADC={(tumor & cell.eq('ADC')).sum()} "
          f"SCC={(tumor & cell.eq('SCC')).sum()}")
    for lab, mask in [("ADC", tumor & cell.eq("ADC")), ("SCC", tumor & cell.eq("SCC")),
                      ("all_tumor", tumor)]:
        cols = expr.columns[mask.fillna(False)]
        if len(cols) < 12:
            continue
        rows += assoc_block(f"GSE19188_{lab}", expr[cols], None, lab, "Affy U133 Plus 2")
    scored["GSE19188"] = {"expr": expr, "scores": signature_scores(expr)[0],
                          "hist": cell, "tumor": tumor}

    # --- GSE50081 ---
    expr = map_expr(read_matrix(os.path.join(GEO2, "GSE50081_matrix.txt.gz")), gpl570)
    meta = parse_meta(os.path.join(GEO2, "GSE50081_matrix.txt.gz")).reindex(expr.columns)
    hist = meta["histology"].astype(str).str.lower()
    print(f"[GSE50081] n={expr.shape[1]} hist={hist.value_counts().to_dict()}")
    for lab, key in [("ADC", "adenocarcinoma"), ("SCC", "squamous cell carcinoma")]:
        cols = expr.columns[hist.eq(key)]
        if len(cols) < 12:
            continue
        rows += assoc_block(f"GSE50081_{lab}", expr[cols], None, lab, "Affy U133 Plus 2")
    scored["GSE50081"] = {"expr": expr, "scores": signature_scores(expr)[0], "hist": hist}

    # --- GSE103584 RNA-seq ---
    raw = pd.read_csv(os.path.join(GEO2, "GSE103584_rnaseq.txt.gz"), sep="\t", index_col=0)
    expr = collapse(np.log2(raw.astype(float).clip(lower=0) + 1))
    meta = parse_meta(os.path.join(GEO2, "GSE103584_matrix.txt.gz"))
    # expression columns are titles (R01-xxx)
    meta = meta.set_index("title")
    keep = [c for c in expr.columns if c in meta.index]
    expr = expr[keep]
    hist = meta.reindex(expr.columns)["histology"].astype(str)
    print(f"[GSE103584] n={expr.shape[1]} hist={hist.value_counts().to_dict()}")
    for lab, key in [("ADC", "Adenocarcinoma"), ("SCC", "Squamous cell carcinoma")]:
        cols = expr.columns[hist.eq(key)]
        if len(cols) < 12:
            continue
        rows += assoc_block(f"GSE103584_{lab}", expr[cols], None, lab, "RNA-seq")
    scored["GSE103584"] = {"expr": expr, "scores": signature_scores(expr)[0], "hist": hist}
    return rows, scored


def tcga_interaction() -> pd.DataFrame:
    """Does the TACSTD2–signature slope differ between LUAD and LUSC?"""
    import statsmodels.api as sm
    with open(os.path.join(PROC_DIR, "scored_cohorts.pkl"), "rb") as fh:
        S = pickle.load(fh)
    rows = []
    for sig in PRIMARY:
        parts = []
        for proj, hist in [("TCGA-LUAD", 0), ("TCGA-LUSC", 1)]:
            c = S[proj]
            df = pd.DataFrame({
                "sig": inverse_normal(c["scores"][sig]),
                "tac": inverse_normal(c["focus_expr"].loc["TACSTD2"]),
                "lusc": hist,
                "purity": c["pheno"].get("purity_absolute"),
            })
            parts.append(df)
        d = pd.concat(parts).dropna()
        X = sm.add_constant(pd.DataFrame({
            "tac": d["tac"], "lusc": d["lusc"], "tac_x_lusc": d["tac"] * d["lusc"],
            "purity": d["purity"],
        }))
        m = sm.OLS(d["sig"], X).fit(cov_type="HC3")
        rows.append({
            "signature": sig, "n": int(m.nobs), "r2": m.rsquared,
            "beta_tac_luad": m.params["tac"], "p_tac_luad": m.pvalues["tac"],
            "beta_interaction": m.params["tac_x_lusc"], "p_interaction": m.pvalues["tac_x_lusc"],
            "beta_lusc": m.params["lusc"], "p_lusc": m.pvalues["lusc"],
            "beta_purity": m.params["purity"], "p_purity": m.pvalues["purity"],
        })
    return pd.DataFrame(rows)


def tcga_os() -> pd.DataFrame:
    from lifelines import CoxPHFitter
    with open(os.path.join(PROC_DIR, "scored_cohorts.pkl"), "rb") as fh:
        S = pickle.load(fh)
    rows = []
    for proj in ["TCGA-LUAD", "TCGA-LUSC"]:
        c = S[proj]
        ph = c["pheno"]
        base = pd.DataFrame({
            "T": pd.to_numeric(ph.get("os_time"), errors="coerce"),
            "E": pd.to_numeric(ph.get("os_event"), errors="coerce"),
            "tac": inverse_normal(c["focus_expr"].loc["TACSTD2"].reindex(ph.index)),
            "purity": ph.get("purity_absolute"),
            "age": pd.to_numeric(ph.get("age"), errors="coerce"),
            "male": (ph.get("sex").astype(str).str.lower().str[0] == "m").astype(float)
                    if "sex" in ph.columns else np.nan,
            "stage": pd.to_numeric(ph.get("stage"), errors="coerce"),
        }, index=ph.index)
        for sig in ["TLS_Cabrita", "B_cell", "T_cell_CD8"]:
            df = base.copy()
            df["sig"] = inverse_normal(c["scores"][sig].reindex(ph.index))
            df["tac_x_sig"] = df["tac"] * df["sig"]
            df = df.dropna(subset=["T", "E", "tac"])
            df = df[df["T"] > 0]
            if df["E"].sum() < 20:
                continue
            for label, cols in [
                ("TACSTD2", ["tac", "T", "E"]),
                ("TACSTD2+purity", ["tac", "purity", "T", "E"]),
                ("TACSTD2+sig+purity", ["tac", "sig", "purity", "T", "E"]),
                ("interaction+purity", ["tac", "sig", "tac_x_sig", "purity", "T", "E"]),
                ("full_clinical", ["tac", "sig", "tac_x_sig", "purity", "age", "male", "stage", "T", "E"]),
            ]:
                use = df[cols].dropna()
                if use["E"].sum() < 15 or use.shape[0] < 40:
                    continue
                try:
                    m = CoxPHFitter()
                    m.fit(use, duration_col="T", event_col="E")
                    for term in [t for t in ["tac", "sig", "tac_x_sig"] if t in m.summary.index]:
                        s = m.summary.loc[term]
                        rows.append({"cohort": proj, "signature": sig, "model": label,
                                     "term": term, "n": len(use), "n_events": int(use["E"].sum()),
                                     "hr": float(s["exp(coef)"]),
                                     "hr_low": float(s["exp(coef) lower 95%"]),
                                     "hr_high": float(s["exp(coef) upper 95%"]),
                                     "p": float(s["p"])})
                except Exception as e:
                    rows.append({"cohort": proj, "signature": sig, "model": label,
                                 "note": str(e)[:120]})
    return pd.DataFrame(rows)


def main() -> None:
    print("[1] new LUSC / histology-split cohorts ...", flush=True)
    rows, _ = build_new_cohorts()
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(TABLES, "replication_lusc_histology.csv"), index=False)
    print("\n=== TACSTD2 after epithelial correction ===")
    t = df[(df.gene == "TACSTD2") & (df.signature.isin(
        ["TLS_Cabrita", "B_cell", "Plasma_cell", "T_cell_CD8"]))]
    print(t[["cohort", "histology", "n", "rho", "p", "rho_adj_epithelial",
             "p_adj_epithelial", "verdict"]].to_string(index=False))
    if "rho_adj_tumorpct" in df.columns:
        print("\n=== GSE17710 TACSTD2 adj tumor_percent ===")
        g = df[(df.cohort == "GSE17710") & (df.gene == "TACSTD2")]
        print(g[["signature", "n", "rho_adj_tumorpct", "p_adj_tumorpct"]].to_string(index=False))

    print("\n[2] TCGA TACSTD2 × LUSC interaction ...", flush=True)
    inter = tcga_interaction()
    inter.to_csv(os.path.join(TABLES, "tcga_tacstd2_histology_interaction.csv"), index=False)
    print(inter.to_string(index=False))

    print("\n[3] TCGA OS Cox ...", flush=True)
    os_df = tcga_os()
    os_df.to_csv(os.path.join(TABLES, "tcga_os_cox.csv"), index=False)
    print(os_df.to_string(index=False))
    print("[done]")


if __name__ == "__main__":
    main()
