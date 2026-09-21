#!/usr/bin/env python3
"""Rank printed TCGA supplement rows for CLDN4-low vs the paper's ISG/STING list.

Inputs are the cell-for-cell reprints of MOESM2 and MOESM5. No TCGA matrix is
downloaded and no p-value or q-value is recomputed. A log2 odds printed as
``<-3`` stays censored.
"""

from __future__ import annotations

import csv
from pathlib import Path

HERE = Path(__file__).resolve().parent

# Named in the article before any supplement row is interpreted.
# ISGs: "ISG15, MX1, OAS1, IFIT1, EIF2AK2, IRF7, STAT1, RSAD2, BST2, and IFI44".
# cGAS-STING set named in the same Results section: cGAS, STING, LC3, TBK1, IRF3,
# plus Beclin-1 and Rab7.
ISG10 = {
    "ISG15",
    "MX1",
    "OAS1",
    "IFIT1",
    "EIF2AK2",
    "IRF7",
    "STAT1",
    "RSAD2",
    "BST2",
    "IFI44",
}
STING_NAMED = {
    "CGAS",
    "MB21D1",
    "STING1",
    "TMEM173",
    "MAP1LC3A",
    "MAP1LC3B",
    "MAP1LC3",
    "TBK1",
    "IRF3",
    "BECN1",
    "RAB7",
    "RAB7A",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def gene_token(symbol: str) -> str:
    """Map a cBioPortal protein id to a gene symbol without inventing a site.

    Phospho-site ids in this table are GENE_P<site> (STAT3_PY705, MTOR_PS2448,
    YAP1_PS127, RPS6_PS240_S244, PRKCB_PS660). The gene is the token before _P.
    Plain gene symbols are unchanged.
    """
    symbol = symbol.strip()
    if "_P" in symbol:
        return symbol.split("_P", 1)[0]
    return symbol


def membership(symbol: str) -> str:
    token = gene_token(symbol)
    if token in ISG10:
        return "isg10"
    if token in STING_NAMED:
        return "sting_named"
    return "none"


def parse_effect(printed: str) -> tuple[str, str]:
    text = printed.strip()
    if text.startswith("<") or text.startswith(">"):
        return "", "yes"
    return text, "no"


def abs_effect_sort_key(printed: str, censored: str) -> float:
    """Larger magnitude sorts first. Censored ``<-3`` is greater than 3, not equal to 3."""
    if censored == "yes":
        bound = printed.strip().lstrip("<>")
        return float(bound) + 0.01
    return abs(float(printed))


def tendency_class(tendency: str) -> str:
    low = tendency.strip().lower().replace("co-ocurrence", "co-occurrence")
    if low == "co-occurrence":
        return "co-occurrence"
    if low == "mutual exclusivity":
        return "mutual exclusivity"
    raise SystemExit(f"unrecognized tendency: {tendency!r}")


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row[key] for key in fields})


def main() -> None:
    isg_rows = read_tsv(HERE / "moesm2_tcga_isg_cbioportal.tsv")
    protein_rows = read_tsv(HERE / "moesm5_tcga_protein_vs_cldn4.tsv")
    if len(isg_rows) != 5:
        raise SystemExit(f"expected 5 ISG rows, got {len(isg_rows)}")
    if len(protein_rows) != 18:
        raise SystemExit(f"expected 18 protein rows, got {len(protein_rows)}")

    ranked: list[dict[str, str]] = []

    for row in isg_rows:
        gene1 = row["Gene 1 (level of expression)"]
        gene2 = row["Gene 2 (level of expression)"]
        g1 = gene_token(gene1.split("(")[0].strip())
        g2 = gene_token(gene2.split("(")[0].strip())
        self_row = g1 == "CLDN4" and g2 == "CLDN4"
        partner = "" if self_row else (g1 if g1 != "CLDN4" else g2)
        if self_row:
            cldn4_state = "low_and_high"
        elif g2 == "CLDN4" and "low" in gene2.lower():
            cldn4_state = "low"
        elif g2 == "CLDN4" and "high" in gene2.lower():
            cldn4_state = "high"
        else:
            raise SystemExit(f"CLDN4 state not on gene 2: {gene1!r} / {gene2!r}")
        tendency = tendency_class(row["Tendency"])
        printed = row["Log2 Odds Ratio"]
        numeric, censored = parse_effect(printed)
        member = "none" if self_row else membership(partner)
        positive_low = (
            cldn4_state == "low"
            and tendency == "co-occurrence"
            and member in {"isg10", "sting_named"}
        )
        ranked.append(
            {
                "table": "MOESM2_ISG",
                "feature": partner if partner else "CLDN4_self",
                "feature_printed_gene1": gene1,
                "feature_printed_gene2": gene2,
                "cldn4_state": cldn4_state,
                "direction_printed": row["Tendency"],
                "direction_class": tendency,
                "log2_effect_printed": printed,
                "log2_effect_numeric": numeric,
                "effect_censored": censored,
                "p_printed": row["p-Value"],
                "q_printed": row["q-Value"],
                "paper_set": member,
                "cldn4_low_positive_isg_sting": "yes" if positive_low else "no",
                "sign_matches_printed_direction": _isg_sign_ok(tendency, printed, censored),
            }
        )

    high_col = "Higher expression in"
    for row in protein_rows:
        symbol = row["Gene"]
        group = row[high_col]
        if group not in {"CLDN4: EXP>1", "CLDN4: EXP<-1"}:
            raise SystemExit(f"unexpected group: {group!r}")
        cldn4_state = "low" if group.endswith("<-1") else "high"
        printed = row["Log2 Ratio"]
        numeric, censored = parse_effect(printed)
        member = membership(symbol)
        positive_low = cldn4_state == "low" and member in {"isg10", "sting_named"}
        ranked.append(
            {
                "table": "MOESM5_protein",
                "feature": symbol,
                "feature_printed_gene1": symbol,
                "feature_printed_gene2": group,
                "cldn4_state": cldn4_state,
                "direction_printed": group,
                "direction_class": "higher_in_cldn4_low" if cldn4_state == "low" else "higher_in_cldn4_high",
                "log2_effect_printed": printed,
                "log2_effect_numeric": numeric,
                "effect_censored": censored,
                "p_printed": row["p-Value"],
                "q_printed": row["q-Value"],
                "paper_set": member,
                "cldn4_low_positive_isg_sting": "yes" if positive_low else "no",
                "sign_matches_printed_direction": _protein_sign_ok(cldn4_state, printed, censored),
            }
        )

    fields = list(ranked[0].keys())
    low_rows = [row for row in ranked if row["cldn4_state"] == "low"]
    low_rows.sort(key=_strength_key)
    counters: dict[str, int] = {}
    for row in low_rows:
        counters[row["table"]] = counters.get(row["table"], 0) + 1
        row["rank_within_table_by_q"] = str(counters[row["table"]])
    fields = list(ranked[0].keys()) + ["rank_within_table_by_q"]
    for row in ranked:
        row.setdefault("rank_within_table_by_q", "")
    write_tsv(HERE / "reanalysis_cldn4_low_ranked.tsv", low_rows, fields)

    aligned = [row for row in low_rows if row["cldn4_low_positive_isg_sting"] == "yes"]
    write_tsv(HERE / "reanalysis_cldn4_low_isg_sting_aligned.tsv", aligned, fields)

    isg_only = [row for row in ranked if row["table"] == "MOESM2_ISG" and row["feature"] != "CLDN4_self"]
    isg_only.sort(key=_strength_key)
    write_tsv(HERE / "reanalysis_isg_rows_ranked.tsv", isg_only, fields)

    present_genes = set()
    for row in ranked:
        if row["feature"] == "CLDN4_self":
            continue
        present_genes.add(gene_token(row["feature"]))
    sting_queries = [
        ("cGAS", ["CGAS", "MB21D1"]),
        ("STING", ["STING1", "TMEM173"]),
        ("LC3", ["MAP1LC3", "MAP1LC3A", "MAP1LC3B"]),
        ("TBK1", ["TBK1"]),
        ("IRF3", ["IRF3"]),
        ("Beclin-1", ["BECN1"]),
        ("Rab7", ["RAB7", "RAB7A"]),
    ]
    missing_isg = sorted(ISG10 - present_genes)
    with (HERE / "reanalysis_named_genes_absent.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["set", "name_in_paper", "symbols_searched", "present_in_moesm2_or_moesm5"])
        for gene in ["ISG15", "MX1", "OAS1", "IFIT1", "EIF2AK2", "IRF7", "STAT1", "RSAD2", "BST2", "IFI44"]:
            writer.writerow(["isg10", gene, gene, "yes" if gene in present_genes else "no"])
        for name, symbols in sting_queries:
            hit = "yes" if any(symbol in present_genes for symbol in symbols) else "no"
            writer.writerow(["sting_named", name, ",".join(symbols), hit])

    sign_fail = [row["feature"] for row in ranked if row["sign_matches_printed_direction"] != "yes"]
    print(f"isg_rows {len(isg_rows)} protein_rows {len(protein_rows)}")
    print(f"cldn4_low_rows {len(low_rows)} aligned_positive {len(aligned)}")
    print("aligned", [row["feature"] for row in aligned])
    print("strongest_cldn4_low", low_rows[0]["feature"], low_rows[0]["q_printed"], low_rows[0]["table"])
    print("isg_rank", [(row["feature"], row["cldn4_state"], row["q_printed"], row["direction_class"]) for row in isg_only])
    print("missing_isg", missing_isg)
    print("sting_present", [name for name, symbols in sting_queries if any(s in present_genes for s in symbols)])
    print("sign_fail", sign_fail)
    print("low_rank")
    for row in low_rows:
        print(
            row["table"],
            row["feature"],
            row["paper_set"],
            row["direction_class"],
            row["log2_effect_printed"],
            row["p_printed"],
            row["q_printed"],
            row["cldn4_low_positive_isg_sting"],
        )


def _isg_sign_ok(tendency: str, printed: str, censored: str) -> str:
    if censored == "yes":
        if tendency == "mutual exclusivity" and printed.strip().startswith("<"):
            return "yes"
        if tendency == "co-occurrence" and printed.strip().startswith(">"):
            return "yes"
        return "no"
    value = float(printed)
    if tendency == "co-occurrence" and value > 0:
        return "yes"
    if tendency == "mutual exclusivity" and value < 0:
        return "yes"
    return "no"


def _protein_sign_ok(cldn4_state: str, printed: str, censored: str) -> str:
    if censored == "yes":
        return "no"
    value = float(printed)
    if cldn4_state == "low" and value < 0:
        return "yes"
    if cldn4_state == "high" and value > 0:
        return "yes"
    return "no"


def _strength_key(row: dict[str, str]) -> tuple[float, float, float]:
    q = float(row["q_printed"])
    p = float(row["p_printed"])
    magnitude = abs_effect_sort_key(row["log2_effect_printed"], row["effect_censored"])
    return (q, p, -magnitude)


if __name__ == "__main__":
    main()
