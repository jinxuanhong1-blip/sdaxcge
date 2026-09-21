#!/usr/bin/env python3
"""Deepen Yamamoto/Bitler MCT Table S1 for NHEJ-down and STING/IFN-up
under CLDN4-low, and record the Bitler-supplement RPPA scrape.

Table S1 is the TCGA ovarian transcriptome, CLDN4-high versus CLDN4-low.
log2 > 0 means higher in CLDN4-high, which is lower under CLDN4-low.

Predicted observational support:
  NHEJ down under CLDN4-low: q<0.05 and log2(high/low) > 0
  STING/IFN up under CLDN4-low: q<0.05 and log2(high/low) < 0

No RPPA fold is imputed. DNA-PKcs protein stays not-reported.
"""

from __future__ import annotations

import urllib.request
from pathlib import Path

import pandas as pd
from scipy.stats import fisher_exact

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "nhej_sting_rppa"
GENE_DIR = Path(__file__).resolve().parent / "genesets"
S1_URL = "https://ndownloader.figshare.com/files/39984958"
UA = "nhej-sting-rppa/1.0 (table-s1-direction)"


def load_genes(name: str) -> set[str]:
    genes = []
    for line in (GENE_DIR / name).read_text().splitlines():
        line = line.strip()
        if not line or line.startswith(">") or line.startswith("REACTOME") or line.startswith("HALLMARK"):
            continue
        genes.append(line.split()[0])
    return set(genes)


def is_histone(gene: str) -> bool:
    return gene.startswith("H2") or gene.startswith("H3") or gene.startswith("H4")


def ensure_s1(dest: Path) -> None:
    if dest.exists() and dest.stat().st_size > 100_000:
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(S1_URL, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=120) as resp:
        dest.write_bytes(resp.read())


def load_s1(path: Path) -> pd.DataFrame:
    df = pd.read_excel(path, header=1)
    df.columns = [str(c).strip() for c in df.columns]
    df["Gene"] = df["Gene"].astype(str).str.strip()
    df["log2"] = pd.to_numeric(df["Log2 of Ratio (unlogged CLDN4 High/CLDN4 Low)"], errors="coerce")
    df["p"] = pd.to_numeric(df["p-Value"], errors="coerce")
    df["q"] = pd.to_numeric(df["q-Value"], errors="coerce")
    return df


def direction_test(df: pd.DataFrame, genes: set[str], predicted: str) -> dict:
    """predicted is 'lower_in_cldn4_low' (log2>0) or 'higher_in_cldn4_low' (log2<0)."""
    sub = df[df["Gene"].isin(genes)]
    missing = sorted(genes - set(sub["Gene"]))
    sig = sub[sub["q"] < 0.05]
    if predicted == "lower_in_cldn4_low":
        hit = sig[sig["log2"] > 0]
        opp = sig[sig["log2"] < 0]
        pred_mask = df["log2"] > 0
    else:
        hit = sig[sig["log2"] < 0]
        opp = sig[sig["log2"] > 0]
        pred_mask = df["log2"] < 0
    allsig = df[df["q"] < 0.05]
    outside = allsig[~allsig["Gene"].isin(genes)]
    set_pred = int(len(hit))
    set_other = int(len(sig) - len(hit))
    out_pred = int(pred_mask.loc[outside.index].sum())
    out_other = int(len(outside) - out_pred)
    if (set_pred + set_other) == 0 or (out_pred + out_other) == 0:
        odds, pval = float("nan"), float("nan")
    else:
        odds, pval = fisher_exact([[set_pred, set_other], [out_pred, out_other]], alternative="two-sided")
    return {
        "n_in_table": int(len(sub)),
        "n_missing": len(missing),
        "missing_genes": ",".join(missing),
        "n_fdr05": int(len(sig)),
        "n_predicted": set_pred,
        "n_opposite": int(len(opp)),
        "fisher_odds_ratio": odds,
        "fisher_p": pval,
        "predicted_genes": ",".join(hit.sort_values("q")["Gene"].tolist()),
        "opposite_genes": ",".join(opp.sort_values("q")["Gene"].tolist()),
    }


def member_rows(df: pd.DataFrame, set_name: str, genes: set[str]) -> list[dict]:
    rows = []
    lookup = df.set_index("Gene")
    for gene in sorted(genes):
        if gene not in lookup.index:
            rows.append({"set": set_name, "gene": gene, "in_table": False})
            continue
        rec = lookup.loc[gene]
        log2 = float(rec["log2"])
        q = float(rec["q"])
        rows.append(
            {
                "set": set_name,
                "gene": gene,
                "in_table": True,
                "histone_in_reactome_nhej": is_histone(gene),
                "log2_high_over_low": log2,
                "p": float(rec["p"]),
                "q": q,
                "fdr05": bool(q < 0.05),
                "lower_in_cldn4_low": bool(log2 > 0),
                "higher_in_cldn4_low": bool(log2 < 0),
                "supports_nhej_down_q05": bool(q < 0.05 and log2 > 0),
                "supports_sting_ifn_up_q05": bool(q < 0.05 and log2 < 0),
            }
        )
    return rows


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    cache = OUT / "cache" / "table_S1.xlsx"
    ensure_s1(cache)
    df = load_s1(cache)

    nhej_full = load_genes("REACTOME_NONHOMOLOGOUS_END_JOINING_NHEJ.txt")
    nhej_core = {g for g in nhej_full if not is_histone(g)}
    nhej_histone = {g for g in nhej_full if is_histone(g)}
    sting = load_genes("REACTOME_STING_MEDIATED_INDUCTION_OF_HOST_IMMUNE_RESPONSES.txt")
    ifna = load_genes("HALLMARK_INTERFERON_ALPHA_RESPONSE.txt")
    ifng = load_genes("HALLMARK_INTERFERON_GAMMA_RESPONSE.txt")
    sting_ifn = sting | ifna | ifng

    sets = [
        ("nhej_core_no_histone", nhej_core, "lower_in_cldn4_low", "NHEJ down under CLDN4-low"),
        ("nhej_reactome_histones_only", nhej_histone, "lower_in_cldn4_low", "histone genes inside Reactome NHEJ, lower under CLDN4-low"),
        ("nhej_reactome_full", nhej_full, "lower_in_cldn4_low", "full Reactome NHEJ including histones, lower under CLDN4-low"),
        ("sting_reactome", sting, "higher_in_cldn4_low", "STING-pathway genes higher under CLDN4-low"),
        ("ifn_alpha_hallmark", ifna, "higher_in_cldn4_low", "Hallmark IFN-alpha higher under CLDN4-low"),
        ("ifn_gamma_hallmark", ifng, "higher_in_cldn4_low", "Hallmark IFN-gamma higher under CLDN4-low"),
        ("sting_union_ifn", sting_ifn, "higher_in_cldn4_low", "STING union IFN hallmarks higher under CLDN4-low"),
    ]
    tests = []
    members = []
    for name, genes, pred, label in sets:
        rec = direction_test(df, genes, pred)
        rec["set"] = name
        rec["prediction"] = label
        tests.append(rec)
        members.extend(member_rows(df, name, genes))

    tests_df = pd.DataFrame(tests).fillna("")
    tests_df.to_csv(OUT / "table_s1_direction_tests.tsv", sep="\t", index=False)
    pd.DataFrame(members).to_csv(OUT / "table_s1_set_members.tsv", sep="\t", index=False)

    # Supplement scrape inventory. Numeric RPPA folds from other experiments
    # are not copied onto the CLDN4-knockdown score.
    inventory = [
        {
            "paper": "Yamamoto 2022 Mol Cancer Ther",
            "pmcid": "PMC8988515",
            "file": "Figshare tables S1-S4 xlsx; figures S1-S4 png; supplement docx",
            "pdf_table": "no",
            "experiment": "OVCAR3 shCLDN4 RPPA described in the text",
            "numeric_rppa_matrix": "no",
            "used_for_cldn4_kd_score": "text direction only; no fold",
            "dna_pkcs": "not_reported",
            "note": "No PDF table. Data statement is request-from-the-author.",
        },
        {
            "paper": "Breed 2024 Cancer Gene Ther CASC4",
            "pmcid": "PMC10874890",
            "file": "41417_2023_703_MOESM1_ESM.pdf Table S1",
            "pdf_table": "yes",
            "experiment": "PEO1 shCASC4 vs shCTRL suspension RPPA",
            "numeric_rppa_matrix": "yes, per-replicate linear and log2 z-scores",
            "used_for_cldn4_kd_score": "no",
            "dna_pkcs": "row absent",
            "note": "53BP1-R-V, XRCC1-R-C, and STING-R-V rows are in this PDF. They are a CASC4 knockdown, not CLDN4. Values were not copied.",
        },
        {
            "paper": "Sottnik 2024 Cancer Res Commun WNT4",
            "pmcid": "PMC10793200",
            "file": "crc-23-0275-s11.xlsx Gene names + L4 matrices",
            "pdf_table": "no",
            "experiment": "103 gynecologic tumors, rs3820282, 484 antibodies",
            "numeric_rppa_matrix": "yes, different contrast",
            "used_for_cldn4_kd_score": "no",
            "dna_pkcs": "PRKDC/DNA-PKcs antibody absent from the 484-name list",
            "note": "The same core's antibody list includes 53BP1, XRCC1, STING, and LIG4. Absence of a DNA-PKcs antibody here is not a CLDN4-KD fold.",
        },
        {
            "paper": "Villagomez 2025 Sci Rep",
            "pmcid": "PMC12603150",
            "file": "MOESM5 docx, TCGA PanCancer Atlas RPPA significant proteins",
            "pdf_table": "no",
            "experiment": "TCGA OV CLDN4 mRNA high vs low, significant proteins only (18 rows)",
            "numeric_rppa_matrix": "significant subset only",
            "used_for_cldn4_kd_score": "no",
            "dna_pkcs": "not in the published significant list",
            "note": "53BP1 and XRCC1 are also absent from that 18-row significant list. This is not the knockdown matrix.",
        },
        {
            "paper": "Sanders 2022 Mol Cancer Ther DUSP",
            "pmcid": "PMC9357222",
            "file": "NIHMS1811007-supplement-1.docx",
            "pdf_table": "not retrieved; PMC package is closed",
            "experiment": "PEO1-OR plus DUSP inhibitor, 483 antibodies",
            "numeric_rppa_matrix": "not in the manuscript XML",
            "used_for_cldn4_kd_score": "no",
            "dna_pkcs": "not named",
            "note": "Manuscript text has no DNA-PKcs, 53BP1, or XRCC1 fold.",
        },
        {
            "paper": "McMellen 2025 Mol Carcinog VDX-111",
            "pmcid": "PMC12481673",
            "file": "NIHMS2111327-supplement-Supplemental_Material.docx",
            "pdf_table": "not retrieved; PMC package is closed",
            "experiment": "OVCAR3 plus VDX-111, 483 antibodies",
            "numeric_rppa_matrix": "manuscript says the dataset is in the supplement",
            "used_for_cldn4_kd_score": "no",
            "dna_pkcs": "not named in the manuscript XML",
            "note": "Not a CLDN4 knockdown. Binary supplement was not downloadable.",
        },
        {
            "paper": "Jordan 2020 Clin Cancer Res",
            "pmcid": "PMC7923250",
            "file": "NIHMS1630007 Table S1-S4 xlsx",
            "pdf_table": "not retrieved; PMC package is closed",
            "experiment": "HGSOC pre/post chemotherapy RPPA",
            "numeric_rppa_matrix": "text quotes other proteins (IL6 log2 0.43), not these three",
            "used_for_cldn4_kd_score": "no",
            "dna_pkcs": "not named",
            "note": "No CLDN4-knockdown RPPA fold in the manuscript text.",
        },
        {
            "paper": "Wheeler 2018 Oncogenesis CBX2",
            "pmcid": "PMC6255906",
            "file": "MOESM6 xlsx sheet T1 CBX2_RPPA",
            "pdf_table": "no",
            "experiment": "CBX2-high tumors, 7 protein rows",
            "numeric_rppa_matrix": "yes, short significant list",
            "used_for_cldn4_kd_score": "no",
            "dna_pkcs": "not in the sheet",
            "note": "53BP1 and XRCC1 are not in the sheet.",
        },
        {
            "paper": "Villagomez 2024 CRC autophagy; Villagomez 2025 CRC genome stability",
            "pmcid": "PMC11218812; PMC11705808",
            "file": "docx supplements",
            "pdf_table": "no",
            "experiment": "CLDN4 modulation, metabolomics and blots",
            "numeric_rppa_matrix": "no RPPA mention in the XML",
            "used_for_cldn4_kd_score": "no",
            "dna_pkcs": "not_reported",
            "note": "No numeric RPPA table.",
        },
    ]
    pd.DataFrame(inventory).to_csv(OUT / "bitler_rppa_supplement_inventory.tsv", sep="\t", index=False)
    print(pd.DataFrame(tests)[["set", "n_in_table", "n_fdr05", "n_predicted", "n_opposite", "fisher_p"]].to_string(index=False))
    print("fdr05 genes", int((df["q"] < 0.05).sum()))


if __name__ == "__main__":
    main()
