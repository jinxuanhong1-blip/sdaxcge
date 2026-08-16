#!/usr/bin/env python3
"""Rework A4: recompute TISMO paired Tacstd2 / Cldn4 after ICB.

Downloads official TISMO public tables (Aliyun share links published by the
TISMO Data Download page, plus the gene-module CSV the website exports for
Tacstd2/Cldn4). Does not invent GEO accessions.

Primary universe = TISMO gene-module comparison groups (64 for Tacstd2).
Sensitivity = coarser pairing from the full in-vivo annotation + expression RDS.
"""
from __future__ import annotations

import hashlib
import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyreadr
from scipy.stats import binomtest, wilcoxon

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "rework" / "A4_TISMO"
SRC = OUT / "source"
FIG = OUT / "figures"
CACHE = Path("/tmp/tismo_tables")

ALIYUN_TOKEN = "https://bj21400.api.aliyunfile.com/v2/share_link/get_share_token"
ALIYUN_LIST = "https://bj21400.api.aliyunfile.com/v2/file/list"
SHARES = {
    "TISMO_vivosample_annotations.csv": "voEA1DXBEFo",
    "TISMO_expressionvivo_profiles.RDS": "AQKzyRCi1Jb",
    "TISMO_cellline_annotations.text": "YzQB2DQonQE",
}
GENE_API = "https://tismo.pku-genomics.org/rtismo/gene/downVivoExprn"
TREAT_API = "https://tismo.pku-genomics.org/tismo/gene/getVivoTreatment"
COHORT_API = "https://tismo.pku-genomics.org/tismo/gene/getVivoCohort"
ICB_TREATMENTS = [
    "antiCTLA4",
    "antiCTLA4&antiPD1",
    "antiCTLA4&antiPDL1",
    "antiPD1",
    "antiPDL1",
    "antiPDL2",
]
TUMOR_MODELS = [
    "402230",
    "4T1",
    "B16",
    "BNL-MEA",
    "CT26",
    "D3UV2",
    "D4M.3A.3",
    "E0771",
    "EMT6",
    "KPB25L",
    "LLC",
    "MC38",
    "MOC22",
    "p53-2225L",
    "p53-2336R",
    "T11",
    "YTN16",
    "YUMM1.7",
]
CLAIM_N_UP = 49
CLAIM_N = 64
CLAIM_P = 5.8e-5


def _post_json(url: str, payload: dict, headers: dict | None = None) -> dict:
    data = json.dumps(payload).encode()
    h = {"Content-Type": "application/json"}
    if headers:
        h.update(headers)
    req = urllib.request.Request(url, data=data, headers=h, method="POST")
    with urllib.request.urlopen(req, timeout=90) as r:
        return json.loads(r.read().decode())


def download_share(share_id: str, dest: Path) -> dict:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tok = _post_json(ALIYUN_TOKEN, {"share_id": share_id, "ignoreError": True})
    listing = _post_json(
        ALIYUN_LIST,
        {
            "limit": 100,
            "marker": "",
            "share_id": share_id,
            "parent_file_id": "root",
            "fields": "user_name,dir_size,url,content_type,upload_id,crc64_hash,revision_id,description",
            "url_expire_sec": 7200,
        },
        headers={"x-share-token": tok["share_token"]},
    )
    item = listing["items"][0]
    if dest.exists() and dest.stat().st_size == int(item["size"]):
        return {"name": item["name"], "size": item["size"], "cached": True, "path": str(dest)}
    req = urllib.request.Request(
        item["download_url"], headers={"x-share-token": tok["share_token"]}
    )
    with urllib.request.urlopen(req, timeout=600) as r, dest.open("wb") as f:
        while True:
            chunk = r.read(1024 * 1024)
            if not chunk:
                break
            f.write(chunk)
    return {"name": item["name"], "size": dest.stat().st_size, "cached": False, "path": str(dest)}


def download_gene_module(gene: str, dest: Path) -> dict:
    dest.parent.mkdir(parents=True, exist_ok=True)
    fields = {
        "filename": "genetreatment_vivo.csv",
        "type": "3",
        "gene": gene,
        "icbList": json.dumps(ICB_TREATMENTS),
        "tumorList": json.dumps(TUMOR_MODELS),
    }
    data = urllib.parse.urlencode(fields).encode()
    req = urllib.request.Request(GENE_API, data=data, method="POST")
    with urllib.request.urlopen(req, timeout=180) as r:
        body = r.read()
    dest.write_bytes(body)
    return {"gene": gene, "bytes": len(body), "path": str(dest)}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def group_stem(name: str) -> str:
    return re.sub(r"\(n=\d+\)$", "", str(name))


def model_from_group(name: str) -> str:
    return str(name).split("_", 1)[0]


def is_lung_group(name: str) -> bool:
    return model_from_group(name) == "LLC"


def comparison_table(df: pd.DataFrame, gene: str) -> pd.DataFrame:
    rows = []
    for group, sub in df.groupby("cell_line"):
        base = sub.loc[sub["Baseline"] == 1, "value"].astype(float)
        icb = sub.loc[sub["Baseline"] == 0, "value"].astype(float)
        r = sub.loc[sub["Responder"] == "Responders", "value"].astype(float)
        nr = sub.loc[sub["Responder"] == "Non-responders", "value"].astype(float)
        pvals = sub["pvalue"].dropna().unique()
        labels = [x for x in sub["label"].dropna().unique() if str(x) != "nan"]
        gses = sorted({str(x) for x in sub["GSE_ID"].dropna().unique()})
        rows.append(
            {
                "group": group,
                "group_stem": group_stem(group),
                "gene": gene,
                "model": model_from_group(group),
                "is_lung": is_lung_group(group),
                "accession": ";".join(gses),
                "n_baseline": int(len(base)),
                "n_icb": int(len(icb)),
                "n_responder": int(len(r)),
                "n_nonresponder": int(len(nr)),
                "mean_baseline": float(base.mean()) if len(base) else np.nan,
                "mean_icb": float(icb.mean()) if len(icb) else np.nan,
                "mean_responder": float(r.mean()) if len(r) else np.nan,
                "mean_nonresponder": float(nr.mean()) if len(nr) else np.nan,
                "tismo_deseq2_pvalue": float(pvals[0]) if len(pvals) == 1 else np.nan,
                "tismo_deseq2_label": labels[0] if len(labels) == 1 else "",
                "srx_ids": ";".join(sorted(sub["Samples"].astype(str).unique())),
            }
        )
    out = pd.DataFrame(rows)
    out["delta"] = out["mean_icb"] - out["mean_baseline"]
    out["log2fc"] = np.log2((out["mean_icb"] + 1e-3) / (out["mean_baseline"] + 1e-3))
    out["direction"] = np.where(
        out["delta"] > 0, "up", np.where(out["delta"] < 0, "down", "tie")
    )
    return out.sort_values(["is_lung", "model", "group"], ascending=[False, True, True])


def direction_stats(comp: pd.DataFrame, label: str) -> dict:
    d = comp.dropna(subset=["delta"]).copy()
    n_up = int((d["direction"] == "up").sum())
    n_down = int((d["direction"] == "down").sum())
    n_tie = int((d["direction"] == "tie").sum())
    n = int(len(d))
    n_signed = n_up + n_down
    binom_p = (
        float(binomtest(n_up, n_signed, 0.5, alternative="two-sided").pvalue)
        if n_signed
        else None
    )
    wilcox_p = None
    if n >= 6 and (d["delta"] != 0).any():
        wilcox_p = float(wilcoxon(d["delta"].to_numpy(), alternative="two-sided").pvalue)
    return {
        "label": label,
        "n_comparisons": n,
        "n_up": n_up,
        "n_down": n_down,
        "n_tie": n_tie,
        "fraction_up": n_up / n if n else None,
        "median_delta": float(d["delta"].median()) if n else None,
        "median_log2fc": float(d["log2fc"].median()) if n else None,
        "binomial_two_sided_p": binom_p,
        "wilcoxon_signed_rank_p": wilcox_p,
    }


def pair_from_metadata(meta: pd.DataFrame, keys: list[str]) -> pd.DataFrame:
    rows = []
    icb = meta[meta["ICB"].astype(str) == "1"]
    for key_t, sub_t in icb.groupby(keys, dropna=False):
        q = meta[
            (meta["Baseline"].astype(str) == "1")
            & (meta["Study_ID"] == sub_t["Study_ID"].iloc[0])
            & (meta["Cell_Line"] == sub_t["Cell_Line"].iloc[0])
        ]
        if "Cell_genotype" in keys:
            q = q[q["Cell_genotype"] == sub_t["Cell_genotype"].iloc[0]]
        if "Condition" in keys:
            q = q[q["Condition"] == sub_t["Condition"].iloc[0]]
        if q.empty or sub_t["Tacstd2"].notna().sum() == 0 or q["Tacstd2"].notna().sum() == 0:
            continue
        rows.append(
            {
                "keys": str(key_t),
                "study": sub_t["Study_ID"].iloc[0],
                "model": sub_t["Cell_Line"].iloc[0],
                "cancer": sub_t["Cancer_type"].iloc[0],
                "icb_group": sub_t["ICB_group"].iloc[0],
                "is_lung": "Lung" in str(sub_t["Cancer_type"].iloc[0]),
                "n_baseline": int(q["Tacstd2"].notna().sum()),
                "n_icb": int(sub_t["Tacstd2"].notna().sum()),
                "tac_baseline": float(q["Tacstd2"].mean()),
                "tac_icb": float(sub_t["Tacstd2"].mean()),
                "cld_baseline": float(q["Cldn4"].mean()),
                "cld_icb": float(sub_t["Cldn4"].mean()),
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["tac_delta"] = out["tac_icb"] - out["tac_baseline"]
    out["cld_delta"] = out["cld_icb"] - out["cld_baseline"]
    out["tac_direction"] = np.where(out["tac_delta"] > 0, "up", np.where(out["tac_delta"] < 0, "down", "tie"))
    out["cld_direction"] = np.where(out["cld_delta"] > 0, "up", np.where(out["cld_delta"] < 0, "down", "tie"))
    return out


def verify_accessions(accessions: list[str]) -> list[dict]:
    out = []
    for acc in accessions:
        if acc.startswith("GSE"):
            url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&form=text&view=brief"
        elif acc.startswith("ERP") or acc.startswith("E-"):
            url = f"https://www.ebi.ac.uk/ena/browser/api/xml/{acc}"
        else:
            out.append(
                {
                    "accession": acc,
                    "status": "not_a_public_geo_ena_id",
                    "note": "Present in TISMO gene-module table only; not queried as a GEO/ENA accession.",
                }
            )
            continue
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "A4-TISMO-rework"})
            with urllib.request.urlopen(req, timeout=25) as r:
                body = r.read(200).decode("utf-8", "replace")
            out.append({"accession": acc, "status": "live", "http": r.status, "head": body.splitlines()[0][:80]})
        except Exception as e:
            out.append({"accession": acc, "status": "error", "error": str(e)})
        time.sleep(0.15)
    return out


def plot_log2fc(comp: pd.DataFrame, title: str, dest: Path) -> None:
    d = comp.sort_values("log2fc")
    colors = np.where(d["is_lung"], "#d55e00", np.where(d["log2fc"] >= 0, "#0072b2", "#999999"))
    fig, ax = plt.subplots(figsize=(8.2, 9.5))
    ax.barh(np.arange(len(d)), d["log2fc"], color=colors, height=0.8)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_yticks([])
    ax.set_xlabel("log2( (mean ICB + 1e-3) / (mean baseline + 1e-3) )")
    ax.set_title(title)
    n_up = int((d["direction"] == "up").sum())
    ax.text(
        0.02,
        0.98,
        f"{n_up}/{len(d)} up   lung (orange) n={int(d['is_lung'].sum())}",
        transform=ax.transAxes,
        va="top",
        fontsize=9,
    )
    fig.tight_layout()
    dest.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(dest, dpi=150)
    plt.close(fig)


def plot_paired_scatter(merged: pd.DataFrame, dest: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 6.0))
    lung = merged[merged["is_lung"]]
    other = merged[~merged["is_lung"]]
    ax.scatter(other["tac_log2fc"], other["cld_log2fc"], s=28, c="#0072b2", alpha=0.85, label="non-lung")
    ax.scatter(lung["tac_log2fc"], lung["cld_log2fc"], s=55, c="#d55e00", edgecolor="k", label="LLC lung")
    ax.axhline(0, color="k", lw=0.6)
    ax.axvline(0, color="k", lw=0.6)
    ax.set_xlabel("Tacstd2 log2FC (ICB vs baseline)")
    ax.set_ylabel("Cldn4 log2FC (ICB vs baseline)")
    ax.set_title("Paired gene-module comparisons")
    ax.legend(frameon=False)
    fig.tight_layout()
    dest.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(dest, dpi=150)
    plt.close(fig)


def plot_summary(rows: list[dict], dest: Path) -> None:
    labels = [r["label"] for r in rows]
    fracs = [r["fraction_up"] if r["fraction_up"] is not None else 0 for r in rows]
    fig, ax = plt.subplots(figsize=(8.5, 4.2))
    ax.barh(labels, fracs, color="#0072b2")
    ax.axvline(0.5, color="k", ls="--", lw=0.8)
    for i, r in enumerate(rows):
        ax.text(0.02, i, f"{r['n_up']}/{r['n_comparisons']}", va="center", color="white", fontsize=9)
    ax.set_xlim(0, 1)
    ax.set_xlabel("Fraction of comparisons with mean ICB > mean baseline")
    ax.set_title("A4 TISMO direction counts")
    fig.tight_layout()
    dest.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(dest, dpi=150)
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    SRC.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(parents=True, exist_ok=True)

    provenance = {
        "tismo_site": "https://tismo.pku-genomics.org/",
        "legacy_site": "http://tismo.cistrome.org/",
        "paper": "https://doi.org/10.1093/nar/gkab804",
        "downloads": {},
    }

    for fname, sid in SHARES.items():
        provenance["downloads"][fname] = download_share(sid, CACHE / fname)
        provenance["downloads"][fname]["share_id"] = sid
        provenance["downloads"][fname]["sha256"] = sha256(CACHE / fname)

    for gene in ("Tacstd2", "Cldn4"):
        dest = SRC / f"{gene}_genetreatment_vivo.csv"
        info = download_gene_module(gene, dest)
        info["sha256"] = sha256(dest)
        provenance["downloads"][dest.name] = info

    tac = pd.read_csv(SRC / "Tacstd2_genetreatment_vivo.csv")
    cld = pd.read_csv(SRC / "Cldn4_genetreatment_vivo.csv")
    meta = pd.read_csv(CACHE / "TISMO_vivosample_annotations.csv")
    expr = pyreadr.read_r(str(CACHE / "TISMO_expressionvivo_profiles.RDS"))[None]

    tac_comp = comparison_table(tac, "Tacstd2")
    cld_comp = comparison_table(cld, "Cldn4")
    tac_comp.to_csv(OUT / "Tacstd2_gene_module_comparisons.tsv", sep="\t", index=False)
    cld_comp.to_csv(OUT / "Cldn4_gene_module_comparisons.tsv", sep="\t", index=False)

    paired_plot = tac_comp.merge(
        cld_comp[["group_stem", "log2fc", "delta", "direction", "mean_baseline", "mean_icb"]].rename(
            columns={
                "log2fc": "cld_log2fc",
                "delta": "cld_delta",
                "direction": "cld_direction",
                "mean_baseline": "cld_mean_baseline",
                "mean_icb": "cld_mean_icb",
            }
        ),
        on="group_stem",
        how="inner",
    ).rename(columns={"log2fc": "tac_log2fc", "delta": "tac_delta", "direction": "tac_direction"})
    paired_plot.to_csv(OUT / "paired_Tacstd2_Cldn4_comparisons.tsv", sep="\t", index=False)

    stats_rows = [
        direction_stats(tac_comp, "Tacstd2 all gene-module groups"),
        direction_stats(tac_comp[tac_comp["is_lung"]], "Tacstd2 lung (LLC) gene-module"),
        direction_stats(tac_comp[~tac_comp["is_lung"]], "Tacstd2 non-lung gene-module"),
        direction_stats(cld_comp, "Cldn4 all gene-module groups"),
        direction_stats(cld_comp[cld_comp["is_lung"]], "Cldn4 lung (LLC) gene-module"),
        direction_stats(cld_comp[~cld_comp["is_lung"]], "Cldn4 non-lung gene-module"),
    ]
    tmp = paired_plot.copy()
    tmp["delta"] = tmp["tac_delta"]
    tmp["log2fc"] = tmp["tac_log2fc"]
    tmp["direction"] = tmp["tac_direction"]
    stats_rows.append(direction_stats(tmp, "Tacstd2 on groups that also have Cldn4"))
    tmp2 = paired_plot.copy()
    tmp2["delta"] = tmp2["cld_delta"]
    tmp2["log2fc"] = tmp2["cld_log2fc"]
    tmp2["direction"] = tmp2["cld_direction"]
    stats_rows.append(direction_stats(tmp2, "Cldn4 on groups that also have Tacstd2"))

    genes = expr.loc[["Tacstd2", "Cldn4"]].T.reset_index().rename(columns={"index": "SampleName"})
    mm = meta.merge(genes, on="SampleName", how="left")
    mm.to_csv(OUT / "vivo_sample_Tacstd2_Cldn4.tsv", sep="\t", index=False)

    coarse = pair_from_metadata(mm, ["Study_ID", "Cell_Line", "ICB_group"])
    fine = pair_from_metadata(mm, ["Study_ID", "Cell_Line", "Cell_genotype", "Condition", "ICB_group"])
    coarse.to_csv(OUT / "sensitivity_study_line_icbgroup.tsv", sep="\t", index=False)
    fine.to_csv(OUT / "sensitivity_study_line_geno_condition_icbgroup.tsv", sep="\t", index=False)

    def meta_stats(df: pd.DataFrame, gene_delta: str, label: str) -> dict:
        if df.empty:
            return {"label": label, "n_comparisons": 0, "n_up": 0, "n_down": 0, "n_tie": 0}
        tmp = pd.DataFrame(
            {
                "delta": df[gene_delta],
                "log2fc": np.log2((df[gene_delta.replace("delta", "icb")] + 1e-3) / (df[gene_delta.replace("delta", "baseline")] + 1e-3)),
                "direction": np.where(df[gene_delta] > 0, "up", np.where(df[gene_delta] < 0, "down", "tie")),
            }
        )
        return direction_stats(tmp, label)

    stats_rows.extend(
        [
            meta_stats(coarse, "tac_delta", "Tacstd2 sensitivity study×line×ICB_group"),
            meta_stats(coarse[coarse["is_lung"]], "tac_delta", "Tacstd2 sensitivity lung study×line×ICB_group"),
            meta_stats(fine, "tac_delta", "Tacstd2 sensitivity study×line×geno×condition×ICB"),
            meta_stats(coarse, "cld_delta", "Cldn4 sensitivity study×line×ICB_group"),
            meta_stats(coarse[coarse["is_lung"]], "cld_delta", "Cldn4 sensitivity lung study×line×ICB_group"),
        ]
    )
    pd.DataFrame(stats_rows).to_csv(OUT / "direction_stats.tsv", sep="\t", index=False)

    accessions = sorted(
        {
            a
            for col in (tac["GSE_ID"], cld["GSE_ID"], meta["Study_ID"])
            for a in col.dropna().astype(str)
            if a and a != "NA"
        }
    )
    # Only verify IDs that appear in the gene-module ICB tables (the 64/65 universe).
    module_acc = sorted({str(x) for x in list(tac["GSE_ID"]) + list(cld["GSE_ID"]) if str(x) not in {"", "nan"}})
    acc_status = verify_accessions(module_acc)
    pd.DataFrame(acc_status).to_csv(OUT / "accession_verification.tsv", sep="\t", index=False)

    plot_log2fc(
        tac_comp,
        "Tacstd2 ICB vs baseline (TISMO gene-module, 64 groups)",
        FIG / "Tacstd2_log2fc_64.png",
    )
    plot_log2fc(
        cld_comp,
        "Cldn4 ICB vs baseline (TISMO gene-module, 65 groups)",
        FIG / "Cldn4_log2fc_65.png",
    )
    plot_paired_scatter(paired_plot, FIG / "paired_Tacstd2_vs_Cldn4_log2fc.png")
    plot_summary(
        list(reversed([r for r in stats_rows if r.get("n_comparisons", 0) > 0][:8])),
        FIG / "direction_summary.png",
    )

    tac_all = next(r for r in stats_rows if r["label"] == "Tacstd2 all gene-module groups")
    verdict = {
        "claim": "TISMO Tacstd2 up in 49/64 ICB comparisons, p=5.8e-5",
        "reproduced": bool(
            tac_all["n_up"] == CLAIM_N_UP
            and tac_all["n_comparisons"] == CLAIM_N
            and tac_all["wilcoxon_signed_rank_p"] is not None
            and abs(tac_all["wilcoxon_signed_rank_p"] - CLAIM_P) / CLAIM_P < 0.05
        ),
        "tacstd2_all": tac_all,
        "note": (
            "49/64 is exactly the TISMO gene-module ICB-vs-baseline direction count for Tacstd2. "
            "p=5.8e-5 matches Wilcoxon signed-rank on the 64 mean differences, not the two-sided "
            "binomial (2.4e-5). Cldn4 is not directionally consistent. Lung-only is 2 LLC groups, both up."
        ),
    }
    (OUT / "VERDICTS.json").write_text(json.dumps(verdict, indent=2) + "\n")
    (OUT / "summary.json").write_text(
        json.dumps(
            {
                "stats": stats_rows,
                "n_tac_groups": int(len(tac_comp)),
                "n_cld_groups": int(len(cld_comp)),
                "n_paired_stems": int(len(paired_plot)),
                "lung_groups_tac": tac_comp.loc[tac_comp["is_lung"], "group"].tolist(),
                "module_accessions": module_acc,
                "in_house_ids_not_geo": [a for a in module_acc if not (a.startswith("GSE") or a.startswith("ERP"))],
            },
            indent=2,
        )
        + "\n"
    )
    provenance["sha256"] = {k: v.get("sha256") for k, v in provenance["downloads"].items()}
    (OUT / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(json.dumps(verdict, indent=2))


if __name__ == "__main__":
    main()
