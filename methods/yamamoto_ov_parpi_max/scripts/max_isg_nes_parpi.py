#!/usr/bin/env python3
"""Maximize Hallmark interferon-alpha NES on public ovarian PARPi RNA-seq.

Acute PARP-inhibitor versus matched vehicle only. Acquired-resistance
contrasts are not in this grid. Ranking statistics and the gene set are
fixed before looking at NES:

  log2FC     mean log2(CPM or TPM + 1), treated minus control
  welch_t    Welch t on that log matrix
  s2n        (mean difference) / (sd_treated + sd_control)

Gene set: Enrichr MSigDB Hallmark 2020 "Interferon Alpha Response"
(the ISG set). "Interferon Gamma Response" is scored in the same call and
reported, and it cannot replace the alpha set as the headline ISG NES.

GSEA is gseapy.prerank, weight 1, 1000 gene-set permutations, seed 123,
min size 15. A contrast can win only if each arm has at least 3 samples.
n=2 arms are still scored and written to the table.

Whole transcriptome means at least 8,000 genes entered the ranked list.
The AmpliSeq panel (GSE120500, 4,604 genes) is flagged and can win, but
the summary also records the best whole-transcriptome NES.
"""

from __future__ import annotations

import gzip
import io
import json
import tarfile
from pathlib import Path

import gseapy as gp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
GEO = Path("/tmp/ov_max/geo")
OUT_T = ROOT / "results" / "tables"
OUT_F = ROOT / "results" / "figures"
OUT_T.mkdir(parents=True, exist_ok=True)
OUT_F.mkdir(parents=True, exist_ok=True)

MIN_N = 3
MIN_SET = 15
WHOLE_TX_GENES = 8000
PERM = 1000
SEED = 123


def log2p(df: pd.DataFrame) -> pd.DataFrame:
    return np.log2(df.clip(lower=0) + 1.0)


def collapse(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.index = df.index.astype(str)
    df = df[~df.index.isin(["", "nan", "None"])]
    df = df.groupby(level=0).sum(numeric_only=True)
    return df


def counts_to_log2cpm(counts: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    counts = collapse(counts)
    lib = counts.sum(axis=0).replace(0, np.nan)
    cpm = counts.div(lib, axis=1) * 1e6
    return cpm, log2p(cpm)


def tpm_to_log(tpm: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    tpm = collapse(tpm)
    return tpm, log2p(tpm)


def symbol_map(path: Path) -> dict[str, str]:
    df = pd.read_csv(path, sep="\t")
    df = df.dropna()
    df = df[df.iloc[:, 1].astype(str).str.len() > 0]
    df = df.drop_duplicates(df.columns[0])
    return dict(zip(df.iloc[:, 0].astype(str), df.iloc[:, 1].astype(str)))


def map_ensembl(counts: pd.DataFrame, mapping: dict[str, str]) -> pd.DataFrame:
    ids = counts.index.astype(str).str.split(".").str[0]
    symbols = ids.map(mapping)
    out = counts.copy()
    out.index = symbols
    out = out[out.index.notna()]
    return collapse(out)


def load_hallmark() -> dict[str, list[str]]:
    sets: dict[str, list[str]] = {}
    with open(GEO / "hallmark.txt") as fh:
        for line in fh:
            parts = [p for p in line.rstrip("\n").split("\t") if p]
            if not parts:
                continue
            sets[parts[0]] = parts[1:]
    need = ["Interferon Alpha Response", "Interferon Gamma Response"]
    missing = [n for n in need if n not in sets]
    if missing:
        raise SystemExit(f"hallmark file missing {missing}")
    return {n: sets[n] for n in need}


def mouse_set(genes: list[str]) -> list[str]:
    return [g[:1].upper() + g[1:].lower() if g else g for g in genes]


def load_237361() -> pd.DataFrame:
    df = pd.read_csv(GEO / "GSE237361_counts.txt.gz", sep="\t", index_col=0)
    return df


def load_243208(human: dict[str, str]) -> pd.DataFrame:
    long = pd.read_csv(GEO / "GSE243208_counts.tsv.gz", sep="\t")
    long["Geneid"] = long["Geneid"].astype(str).str.split(".").str[0]
    wide = long.pivot_table(index="Geneid", columns="sample", values="count", aggfunc="sum")
    wide = map_ensembl(wide, human)
    rename = {
        "MP1237AA4P4": "ola_1",
        "MP1237AA5P5": "ola_2",
        "MP1237AA6P6": "ola_3",
        "MP1237AA28P28": "dmso_1",
        "MP1237AA29P29": "dmso_2",
        "MP1237AA30P30": "dmso_3",
    }
    missing = [c for c in rename if c not in wide.columns]
    if missing:
        raise SystemExit(f"GSE243208 missing samples {missing}")
    return wide[list(rename)].rename(columns=rename)


def load_285827() -> pd.DataFrame:
    frames = []
    with tarfile.open(GEO / "GSE285827_RAW.tar") as tar:
        for mem in tar.getmembers():
            if not mem.name.endswith(".gz"):
                continue
            raw = tar.extractfile(mem).read()
            one = pd.read_csv(io.BytesIO(raw), sep="\t", compression="gzip")
            gene_col, sample_col = one.columns[0], one.columns[1]
            ser = one.groupby(one[gene_col].astype(str))[sample_col].sum()
            ser.name = sample_col
            frames.append(ser)
    return pd.concat(frames, axis=1).fillna(0.0)


def load_120500() -> pd.DataFrame:
    xls_gz = GEO / "GSE120500_panel.xls.gz"
    xls_path = GEO / "GSE120500_panel.xls"
    if not xls_path.exists():
        with gzip.open(xls_gz, "rb") as src, open(xls_path, "wb") as dst:
            dst.write(src.read())
    amp = pd.read_excel(xls_path, engine="xlrd")
    amp = amp.rename(columns={amp.columns[0]: "gene"}).groupby("gene", as_index=True).sum(numeric_only=True)
    spots = pd.read_csv(ROOT / "data" / "gse120500_sra_spots.tsv", sep="\t")
    assign = {}
    for col, total in amp.sum(axis=0).items():
        rel = (spots["spots"] - float(total)).abs() / spots["spots"]
        j = int(rel.to_numpy().argmin())
        if float(rel.iloc[j]) > 0.02:
            raise SystemExit(f"GSE120500 column {col} does not match an SRA spot count")
        title = str(spots.iloc[j]["title"])
        group = "control" if title.startswith("Control") else "olaparib"
        other = spots.loc[~spots["title"].astype(str).str.startswith("Control" if group == "control" else "PARPi"), "spots"]
        other_rel = float(((other - float(total)).abs() / other).min())
        if other_rel < 0.03:
            raise SystemExit(f"GSE120500 group call ambiguous for {col}")
        assign[col] = f"{group}_{len(assign)}"
    # Stable names within group.
    counts = {"control": 0, "olaparib": 0}
    renamed = {}
    for col, label in assign.items():
        group = label.split("_", 1)[0]
        counts[group] += 1
        renamed[col] = f"{group}_{counts[group]}"
    return amp.rename(columns=renamed)


def load_191231() -> pd.DataFrame:
    df = pd.read_excel(GEO / "extra" / "GSE191231.xlsx", sheet_name=0)
    df = df.rename(columns={"id": "symbol"})
    tpm_cols = {
        "S_834C_1.fq.gz.genes.results.TPM": "vehicle_834",
        "S_846C_1.fq.gz.genes.results.TPM": "vehicle_846",
        "S_827_1.fq.gz.genes.results.TPM": "olaparib_827",
        "S_878_1.fq.gz.genes.results.TPM": "olaparib_878",
    }
    missing = [c for c in tpm_cols if c not in df.columns]
    if missing:
        raise SystemExit(f"GSE191231 missing {missing}")
    out = df.set_index("symbol")[list(tpm_cols)].rename(columns=tpm_cols)
    return out.apply(pd.to_numeric, errors="coerce").fillna(0.0)


def load_163854() -> pd.DataFrame:
    return pd.read_csv(GEO / "extra" / "GSE163854_counts.txt.gz", sep="\t", index_col=0)


def load_264286() -> pd.DataFrame:
    df = pd.read_csv(GEO / "extra" / "GSE264286_counts.txt.gz", sep="\t", index_col=0)
    # GEO Sample_description is the matrix column. Sample order is GSM8215511-518:
    # Control, Control, Niraparib, Niraparib, Lenvatinib, Lenvatinib, N+L, N+L.
    rename = {
        "MQ220701_001": "control_1",
        "MQ220701_002": "control_2",
        "MQ220701_003": "niraparib_1",
        "MQ220701_004": "niraparib_2",
        "MQ220701_005": "lenvatinib_1",
        "MQ220701_006": "lenvatinib_2",
        "MQ220701_007": "combo_1",
        "MQ220701_008": "combo_2",
    }
    missing = [c for c in rename if c not in df.columns]
    if missing:
        raise SystemExit(f"GSE264286 missing {missing}")
    return df[list(rename)].rename(columns=rename)


def load_309870() -> pd.DataFrame:
    df = pd.read_csv(GEO / "extra" / "GSE309870_TPM.txt.gz", sep="\t", encoding="utf-16")
    # DMSO columns are stored as shCoontrol (source spelling).
    df = df.set_index("symbol")
    df = df.drop(columns=[c for c in df.columns if c in {"entrezgene", "name"}], errors="ignore")
    return df.apply(pd.to_numeric, errors="coerce").fillna(0.0)


def load_246085() -> pd.DataFrame:
    df = pd.read_csv(GEO / "extra" / "GSE246085_Allgene.txt.gz", sep="\t")
    count_cols = [c for c in df.columns if c.startswith("Count_")]
    out = df.set_index("symbol")[count_cols]
    out.columns = [c.replace("Count_ ", "").replace("Count_", "") for c in out.columns]
    return out.apply(pd.to_numeric, errors="coerce").fillna(0.0)


def load_302064(mouse: dict[str, str]) -> pd.DataFrame:
    frames = []
    with tarfile.open(GEO / "extra" / "GSE302064_RAW.tar") as tar:
        for mem in tar.getmembers():
            if not mem.name.endswith(".gz"):
                continue
            raw = tar.extractfile(mem).read()
            one = pd.read_csv(io.BytesIO(raw), sep="\t", compression="gzip", usecols=["Feature", "Count"])
            name = Path(mem.name).name
            if "Scr-UT" in name:
                label = "scr_vehicle_" + name.split("Scr-UT_")[1][0]
            elif "Scr_ola" in name:
                label = "scr_olaparib_" + name.split("Scr_ola")[1][0]
            elif "KO1-UT" in name:
                label = "ko_vehicle_" + name.split("KO1-UT_")[1][0]
            elif "KO1-ola" in name:
                label = "ko_olaparib_" + name.split("KO1-ola")[1][0]
            else:
                continue
            one["symbol"] = one["Feature"].astype(str).str.split(".").str[0].map(mouse)
            one = one.dropna(subset=["symbol"]).groupby("symbol")["Count"].sum()
            one.name = label
            frames.append(one)
    return pd.concat(frames, axis=1).fillna(0.0)


def contrasts(mats: dict[str, pd.DataFrame]) -> list[dict]:
    m237 = mats["GSE237361"]
    m243 = mats["GSE243208"]
    m285 = mats["GSE285827"]
    m120 = mats["GSE120500"]
    m191 = mats["GSE191231"]
    m163 = mats["GSE163854"]
    m264 = mats["GSE264286"]
    m309 = mats["GSE309870"]
    m246 = mats["GSE246085"]
    m302 = mats["GSE302064"]
    rows = [
        dict(contrast_id="GSE237361_UWB_olaparib_96h", series="GSE237361", species="human",
             model="UWB1.289 BRCA1-mutant cell, 3.5 uM, 96 h", drug="olaparib", scale="counts",
             matrix=m237, treat=[c for c in m237.columns if "UWB_olaparib" in c],
             ctrl=[c for c in m237.columns if "UWB_DMSO" in c], panel=False),
        dict(contrast_id="GSE237361_OVCAR3_olaparib_96h", series="GSE237361", species="human",
             model="OVCAR3 BRCA1-WT cell, 7.5 uM, 96 h", drug="olaparib", scale="counts",
             matrix=m237, treat=[c for c in m237.columns if "OVC_olaparib" in c],
             ctrl=[c for c in m237.columns if "OVC_DMSO" in c], panel=False),
        dict(contrast_id="GSE243208_UWB_olaparib_24h", series="GSE243208", species="human",
             model="UWB1.289, 24 h (GSM 4 uM; series text 10 uM)", drug="olaparib", scale="counts",
             matrix=m243, treat=["ola_1", "ola_2", "ola_3"], ctrl=["dmso_1", "dmso_2", "dmso_3"], panel=False),
        dict(contrast_id="GSE285827_OVCAR3_talazoparib", series="GSE285827", species="human",
             model="OVCAR3 talazoparib vs DMSO, three dates", drug="talazoparib", scale="counts",
             matrix=m285, treat=[c for c in m285.columns if c.startswith("ovcar3_talap_")],
             ctrl=[c for c in m285.columns if c.startswith("ovcar3_dmso_")], panel=False),
        dict(contrast_id="GSE285827_OVCAR3_veliparib", series="GSE285827", species="human",
             model="OVCAR3 veliparib vs DMSO, three dates", drug="veliparib", scale="counts",
             matrix=m285, treat=[c for c in m285.columns if c.startswith("ovcar3_velip_")],
             ctrl=[c for c in m285.columns if c.startswith("ovcar3_dmso_")], panel=False),
        dict(contrast_id="GSE285827_CAOV3_talazoparib", series="GSE285827", species="human",
             model="CAOV3 talazoparib vs DMSO, three dates", drug="talazoparib", scale="counts",
             matrix=m285, treat=[c for c in m285.columns if c.startswith("caov3_talap_")],
             ctrl=[c for c in m285.columns if c.startswith("caov3_dmso_")], panel=False),
        dict(contrast_id="GSE285827_CAOV3_veliparib", series="GSE285827", species="human",
             model="CAOV3 veliparib vs DMSO, three dates", drug="veliparib", scale="counts",
             matrix=m285, treat=[c for c in m285.columns if c.startswith("caov3_velip_")],
             ctrl=[c for c in m285.columns if c.startswith("caov3_dmso_")], panel=False),
        dict(contrast_id="GSE120500_Brca1def_tumor_olaparib", series="GSE120500", species="mouse",
             model="Brca1-deficient mouse ovarian tumors, bulk AmpliSeq", drug="olaparib", scale="counts",
             matrix=m120, treat=[c for c in m120.columns if c.startswith("olaparib_")],
             ctrl=[c for c in m120.columns if c.startswith("control_")], panel=True),
        dict(contrast_id="GSE309870_A2780_shControl_olaparib", series="GSE309870", species="human",
             model="A2780 scramble, olaparib vs DMSO", drug="olaparib", scale="tpm",
             matrix=m309, treat=[c for c in m309.columns if c.startswith("A2780_shControl_Olaparib")],
             ctrl=[c for c in m309.columns if "shCoontrol_DMSO" in c], panel=False),
        dict(contrast_id="GSE309870_A2780_shA3B_olaparib", series="GSE309870", species="human",
             model="A2780 APOBEC3B knockdown, olaparib vs DMSO", drug="olaparib", scale="tpm",
             matrix=m309, treat=[c for c in m309.columns if c.startswith("A2780_shA3B_Olaparib")],
             ctrl=[c for c in m309.columns if c.startswith("A2780_shA3B_DMSO")], panel=False),
        dict(contrast_id="GSE246085_ES2_olaparib", series="GSE246085", species="human",
             model="ES2 ovarian clear-cell, olaparib vs control, n=6", drug="olaparib", scale="counts",
             matrix=m246, treat=[c for c in m246.columns if c.startswith("Ola-")],
             ctrl=[c for c in m246.columns if c.startswith("Control-")], panel=False),
        dict(contrast_id="GSE302064_ID8OR_scramble_olaparib", series="GSE302064", species="mouse",
             model="ID8 olaparib-resistant scramble, acute olaparib vs vehicle", drug="olaparib", scale="counts",
             matrix=m302, treat=[c for c in m302.columns if c.startswith("scr_olaparib_")],
             ctrl=[c for c in m302.columns if c.startswith("scr_vehicle_")], panel=False),
        dict(contrast_id="GSE302064_ID8OR_Rad52KO_olaparib", series="GSE302064", species="mouse",
             model="ID8 olaparib-resistant Rad52 KO, acute olaparib vs vehicle", drug="olaparib", scale="counts",
             matrix=m302, treat=[c for c in m302.columns if c.startswith("ko_olaparib_")],
             ctrl=[c for c in m302.columns if c.startswith("ko_vehicle_")], panel=False),
        dict(contrast_id="GSE163854_PH039_first_passage_niraparib", series="GSE163854", species="human",
             model="HGSOC PDX PH039 first-passage niraparib vs untreated source", drug="niraparib", scale="counts",
             matrix=m163, treat=["F1A3", "F1A4", "F1A5", "F251"], ctrl=["P2A1", "P2A2", "P2A3"], panel=False),
        dict(contrast_id="GSE191231_BRCAWT_PDX_olaparib", series="GSE191231", species="human",
             model="BRCA-WT ovarian PDX, olaparib vs vehicle", drug="olaparib", scale="tpm",
             matrix=m191, treat=["olaparib_827", "olaparib_878"], ctrl=["vehicle_834", "vehicle_846"], panel=False),
        dict(contrast_id="GSE264286_OCCC_PDX_niraparib", series="GSE264286", species="human",
             model="ovarian clear-cell PDX, niraparib vs control", drug="niraparib", scale="counts",
             matrix=m264, treat=["niraparib_1", "niraparib_2"], ctrl=["control_1", "control_2"], panel=False),
    ]
    for row in rows:
        missing = [c for c in row["treat"] + row["ctrl"] if c not in row["matrix"].columns]
        if missing or not row["treat"] or not row["ctrl"]:
            raise SystemExit(f"{row['contrast_id']} columns missing {missing}; treat={row['treat']} ctrl={row['ctrl']}")
    return rows


def rank_metrics(abundance: pd.DataFrame, logx: pd.DataFrame, treat: list[str], ctrl: list[str]) -> dict[str, pd.Series]:
    keep = (abundance[treat].mean(axis=1) >= 1.0) | (abundance[ctrl].mean(axis=1) >= 1.0)
    logx = logx.loc[keep]
    mt = logx[treat].mean(axis=1)
    mc = logx[ctrl].mean(axis=1)
    out = {"log2fc": (mt - mc).replace([np.inf, -np.inf], np.nan).dropna()}
    if len(treat) >= 2 and len(ctrl) >= 2:
        tt = stats.ttest_ind(logx[treat], logx[ctrl], axis=1, equal_var=False, nan_policy="omit")
        welch = pd.Series(tt.statistic, index=logx.index).replace([np.inf, -np.inf], np.nan).dropna()
        out["welch_t"] = welch
        st = logx[treat].std(axis=1, ddof=1)
        sc = logx[ctrl].std(axis=1, ddof=1)
        s2n = ((mt - mc) / (st + sc + 1e-6)).replace([np.inf, -np.inf], np.nan).dropna()
        out["s2n"] = s2n
    return out


def run_gsea(ranks: pd.Series, gene_sets: dict[str, list[str]]) -> pd.DataFrame:
    rnk = ranks.sort_values(ascending=False)
    rnk = rnk[~rnk.index.duplicated(keep="first")]
    pre = gp.prerank(
        rnk=rnk,
        gene_sets=gene_sets,
        permutation_num=PERM,
        weight=1.0,
        min_size=MIN_SET,
        max_size=500,
        seed=SEED,
        threads=4,
        no_plot=True,
        outdir=None,
        verbose=False,
    )
    res = pre.res2d.copy()
    res.columns = [str(c) for c in res.columns]
    return res


def nes_of(res: pd.DataFrame, term: str) -> dict:
    name_col = "Term" if "Term" in res.columns else "Name"
    hit = res[res[name_col].astype(str).eq(term)]
    if hit.empty:
        return {"nes": np.nan, "pval": np.nan, "fdr": np.nan, "n_set": 0, "lead": ""}
    row = hit.iloc[0]
    def num(*keys):
        for k in keys:
            if k in row.index and pd.notna(row[k]):
                return float(row[k])
        return np.nan
    lead = ""
    for k in ("Lead_genes", "Lead genes"):
        if k in row.index and pd.notna(row[k]):
            lead = str(row[k])
            break
    tag = str(row["Tag %"]) if "Tag %" in row.index else ""
    n_set = 0
    if tag and "/" in tag:
        try:
            n_set = int(tag.split("/")[-1])
        except ValueError:
            n_set = 0
    return {
        "nes": num("NES"),
        "pval": num("NOM p-val", "NOM p-value", "pval"),
        "fdr": num("FDR q-val", "FDR q-value", "FDR"),
        "n_set": n_set,
        "lead": lead.split(";")[0][:80] if lead else "",
    }


def main() -> None:
    print("[load] matrices")
    human = symbol_map(GEO / "human_ensembl_symbol.tsv")
    mouse = symbol_map(GEO / "mouse_ensembl_symbol.tsv")
    hallmark = load_hallmark()
    mats = {
        "GSE237361": load_237361(),
        "GSE243208": load_243208(human),
        "GSE285827": load_285827(),
        "GSE120500": load_120500(),
        "GSE191231": load_191231(),
        "GSE163854": load_163854(),
        "GSE264286": load_264286(),
        "GSE309870": load_309870(),
        "GSE246085": load_246085(),
        "GSE302064": load_302064(mouse),
    }
    specs = contrasts(mats)
    rows = []
    for spec in specs:
        mat = spec["matrix"][spec["treat"] + spec["ctrl"]]
        if spec["scale"] == "counts":
            abundance, logx = counts_to_log2cpm(mat)
        else:
            abundance, logx = tpm_to_log(mat)
        metrics = rank_metrics(abundance, logx, spec["treat"], spec["ctrl"])
        sets = hallmark if spec["species"] == "human" else {k: mouse_set(v) for k, v in hallmark.items()}
        n_treat, n_ctrl = len(spec["treat"]), len(spec["ctrl"])
        print(f"[gsea] {spec['contrast_id']} n={n_treat} vs {n_ctrl} genes={abundance.shape[0]}")
        for stat, ranks in metrics.items():
            present = {name: [g for g in genes if g in ranks.index] for name, genes in sets.items()}
            # Drop a set that would fall under min_size so the other set still runs.
            present = {k: v for k, v in present.items() if len(v) >= MIN_SET}
            if not present:
                for name in sets:
                    rows.append(_blank(spec, stat, name, len(ranks), n_treat, n_ctrl, 0))
                continue
            res = run_gsea(ranks, present)
            if rows == []:
                print("  gsea columns", list(res.columns))
            for name in sets:
                info = nes_of(res, name) if name in present else {"nes": np.nan, "pval": np.nan, "fdr": np.nan, "n_set": len(sets[name]), "lead": ""}
                if name not in present:
                    info["n_set"] = 0
                whole = int(len(ranks) >= WHOLE_TX_GENES and not spec["panel"])
                eligible = bool(n_treat >= MIN_N and n_ctrl >= MIN_N and name == "Interferon Alpha Response" and np.isfinite(info["nes"]))
                rows.append(
                    {
                        "contrast_id": spec["contrast_id"],
                        "series": spec["series"],
                        "species": spec["species"],
                        "model": spec["model"],
                        "drug": spec["drug"],
                        "n_treat": n_treat,
                        "n_ctrl": n_ctrl,
                        "scale": spec["scale"],
                        "panel_limited": bool(spec["panel"] or len(ranks) < WHOLE_TX_GENES),
                        "whole_transcriptome": bool(whole),
                        "n_genes_ranked": int(len(ranks)),
                        "rank_stat": stat,
                        "gene_set": name,
                        "n_set_in_list": int(info["n_set"] or len(present.get(name, []))),
                        "nes": info["nes"],
                        "nominal_p": info["pval"],
                        "fdr_q": info["fdr"],
                        "lead_gene": info["lead"],
                        "eligible": eligible,
                    }
                )
                print(f"  {stat:8} {name:28} NES={info['nes']}")
    table = pd.DataFrame(rows)
    alpha = table[table["gene_set"].eq("Interferon Alpha Response") & table["eligible"]]
    if alpha.empty or not np.isfinite(alpha["nes"]).any():
        raise SystemExit("no eligible interferon-alpha NES")
    winner = alpha.sort_values(["nes", "n_treat", "n_genes_ranked"], ascending=[False, False, False]).iloc[0]
    whole = alpha[alpha["whole_transcriptome"]]
    whole_winner = None
    if not whole.empty and np.isfinite(whole["nes"]).any():
        whole_winner = whole.sort_values(["nes", "n_treat"], ascending=[False, False]).iloc[0]
    gamma = table[table["gene_set"].eq("Interferon Gamma Response") & table["n_treat"].ge(MIN_N) & table["n_ctrl"].ge(MIN_N)]
    gamma_winner = None
    if not gamma.empty and np.isfinite(gamma["nes"]).any():
        gamma_winner = gamma.sort_values("nes", ascending=False).iloc[0]
    table.to_csv(OUT_T / "isg_nes_grid.tsv", sep="\t", index=False, float_format="%.6g")
    summary = {
        "gseapy": gp.__version__,
        "permutations": PERM,
        "weight": 1.0,
        "seed": SEED,
        "min_n_per_arm": MIN_N,
        "isg_set": "MSigDB Hallmark 2020 Interferon Alpha Response",
        "winner": winner.to_dict(),
        "whole_transcriptome_winner": None if whole_winner is None else whole_winner.to_dict(),
        "ifng_winner_n_at_least_3": None if gamma_winner is None else gamma_winner.to_dict(),
        "n_rows": int(len(table)),
        "n_eligible_alpha": int(alpha.shape[0]),
    }
    (OUT_T / "isg_nes_summary.json").write_text(json.dumps(summary, indent=2, default=float))
    plot_nes(table, winner["contrast_id"], winner["rank_stat"])
    print("[winner]", winner["contrast_id"], winner["rank_stat"], float(winner["nes"]))


def _blank(spec, stat, name, n_genes, n_treat, n_ctrl, n_set) -> dict:
    return {
        "contrast_id": spec["contrast_id"],
        "series": spec["series"],
        "species": spec["species"],
        "model": spec["model"],
        "drug": spec["drug"],
        "n_treat": n_treat,
        "n_ctrl": n_ctrl,
        "scale": spec["scale"],
        "panel_limited": True,
        "whole_transcriptome": False,
        "n_genes_ranked": int(n_genes),
        "rank_stat": stat,
        "gene_set": name,
        "n_set_in_list": int(n_set),
        "nes": np.nan,
        "nominal_p": np.nan,
        "fdr_q": np.nan,
        "lead_gene": "",
        "eligible": False,
    }


def plot_nes(table: pd.DataFrame, winner_id: str, winner_stat: str) -> None:
    sub = table[table["gene_set"].eq("Interferon Alpha Response") & table["rank_stat"].eq("log2fc")].copy()
    sub = sub.sort_values("nes")
    fig, ax = plt.subplots(figsize=(8.4, 6.4))
    colors = ["#b85c38" if (r.contrast_id == winner_id and winner_stat == "log2fc") else "#1f4e79" for r in sub.itertuples()]
    ax.barh(sub["contrast_id"], sub["nes"], color=colors)
    ax.axvline(0, color="0.4", lw=0.8)
    ax.set_xlabel("Hallmark interferon-alpha NES (log2FC rank, weight 1)")
    ax.set_title("Ovarian PARPi versus vehicle")
    fig.tight_layout()
    fig.savefig(OUT_F / "isg_nes_log2fc.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()
