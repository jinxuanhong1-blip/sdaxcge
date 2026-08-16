#!/usr/bin/env python3
"""
Stage 05 - assemble the leftover catalog + bilingual WRITEUP.

Reads the CSVs from stages 01-04 and writes:
  results/w200/GEO_2023/leftover_catalog.csv
  results/w200/GEO_2023/README.md
  notes/w200_geo_2023/WRITEUP.md
"""
import csv
from pathlib import Path

from geo_common import OUT, NOTES, ROOT


def read_csv(p):
    if not Path(p).exists():
        return []
    return list(csv.DictReader(Path(p).open()))


def main():
    cand = read_csv(OUT / "search_candidates_2023.csv")
    cls = read_csv(OUT / "leftover_classification.csv")
    probe = {r["accession"]: r for r in read_csv(OUT / "leftover_probe.csv")}
    tests = read_csv(OUT / "analysis" / "marker_outcome_tests.csv")
    skipped = read_csv(OUT / "analysis" / "skipped_usable_maybe.csv")

    catalog = []
    for r in cls:
        acc = r["accession"]
        pr = probe.get(acc, {})
        leftover = r["leftover"] in (True, "True", "true", "1")
        usable = pr.get("usable_maybe") in (True, "True", "true", "1")
        has_t = pr.get("has_TACSTD2", "")
        has_c = pr.get("has_CLDN4", "")
        verdict = _verdict(r, pr, leftover, usable)
        catalog.append({
            "accession": acc,
            "pdat": r.get("pdat", ""),
            "n_samples": r.get("n_samples", ""),
            "title": r.get("title", ""),
            "gdsType": r.get("gdsType", ""),
            "is_human": r.get("is_human", ""),
            "is_lung": r.get("is_lung", ""),
            "is_ici": r.get("is_ici", ""),
            "is_expression": r.get("is_expression", ""),
            "relevant_text": r.get("relevant_text", ""),
            "study_category": r.get("study_category", ""),
            "leftover_status": r.get("leftover_status", ""),
            "already_analysed_note": r.get("already_analysed_note", ""),
            "has_outcome_annotation": pr.get("has_outcome_annotation", ""),
            "outcome_keys": pr.get("outcome_keys", ""),
            "has_TACSTD2": has_t,
            "has_CLDN4": has_c,
            "peeked_file": pr.get("peeked_file", ""),
            "n_processed_files": pr.get("n_processed_files", ""),
            "under_2gb": pr.get("under_2gb", ""),
            "usable_maybe": usable,
            "verdict": verdict,
        })

    with (OUT / "leftover_catalog.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(catalog[0].keys()))
        w.writeheader()
        w.writerows(catalog)

    n = len(catalog)
    n_left_rel = sum(1 for x in catalog if x["leftover_status"] == "leftover_relevant")
    n_already = sum(1 for x in catalog if x["leftover_status"] == "already_analysed")
    n_usable = sum(1 for x in catalog if x["usable_maybe"] is True or x["usable_maybe"] == True)
    n_tests = len(tests)
    n_sig = sum(1 for t in tests
                if _is_sig(t.get("p_value"))
                and t.get("outcome") != "spearman_within_cohort")

    writeup = _writeup(
        n=n, n_left_rel=n_left_rel, n_already=n_already,
        n_usable=n_usable, n_tests=n_tests, n_sig=n_sig,
        catalog=catalog, tests=tests, skipped=skipped,
    )
    (NOTES / "WRITEUP.md").write_text(writeup)
    (OUT / "README.md").write_text(
        "# w200 / GEO_2023 leftover lung ICI TACSTD2/CLDN4\n\n"
        "See `notes/w200_geo_2023/WRITEUP.md` for the bilingual report.\n\n"
        f"- search candidates: {n}\n"
        f"- leftover relevant: {n_left_rel}\n"
        f"- already analysed (excluded): {n_already}\n"
        f"- leftover usable_maybe heuristic hits: {n_usable} "
        f"(cell-line / melanoma TIL; rejected)\n"
        f"- leftover actually tested: GSE248378 only\n"
        f"- outcome tests with p < 0.05: {n_sig}\n"
        f"- TACSTD2 vs recurrence (non-MPR, post-ICI): see analysis/\n"
    )
    print(f"catalog {n}; leftover_relevant {n_left_rel}; usable {n_usable}; tests {n_tests}")


def _is_sig(p):
    try:
        return float(p) < 0.05
    except (TypeError, ValueError):
        return False


def _verdict(r, pr, leftover, usable):
    if r.get("leftover_status") == "already_analysed":
        return r.get("already_analysed_note") or "already analysed"
    if r.get("leftover_status") != "leftover_relevant":
        return "not leftover-relevant (failed human+lung+ICI+expression text filter)"
    if r.get("accession") == "GSE248378":
        return ("leftover ANALYSED: post-durvalumab non-MPR bulk FPKM, both genes present; "
                "GEO has arm/histology only; recurrence recovered from Nat Commun source "
                "data Fig.5d via unique ITGAE FPKM (9 recur / 20 no recur)")
    if r.get("accession") == "GSE217451":
        return "leftover: H1650 si-hMENA cell-line RNA-seq; both genes present; not a patient ICI cohort"
    if r.get("accession") == "GSE224216":
        return "leftover: H2030 si-hMENA cell-line RNA-seq; both genes present; not a patient ICI cohort"
    if r.get("accession") == "GSE224099":
        return "leftover: melanoma TIL CD4 T-cell subsets (not lung epithelium); not TACSTD2/CLDN4 tumour ICI"
    if r.get("accession") == "GSE235603":
        return "leftover: SuperSeries of tumour-infiltrating Treg scRNA; no bulk tumour matrix"
    cat = r.get("study_category", "")
    if "single_cell" in cat:
        return "leftover: single-cell / sorted-immune; not bulk tumour epithelium for TACSTD2/CLDN4 vs ICI"
    if "platelet" in cat:
        return "leftover: platelet RNA, not tumour"
    if "blood_or_pbmc" in cat and "bulk_or_other" not in cat:
        return "leftover: blood/PBMC, not tumour epithelium"
    if "cell_line" in cat:
        return "leftover: cell-line / in-vitro, not a patient ICI outcome cohort"
    if "chromatin" in cat and "bulk_or_other" not in cat:
        return "leftover: chromatin assay (ATAC/ChIP), not gene expression of TACSTD2/CLDN4"
    if "flow_or_cytof" in cat:
        return "leftover: flow/CyTOF, not transcriptome"
    if "targeted_panel" in cat:
        return "leftover: targeted panel — gene presence must be checked; often lacks TACSTD2/CLDN4"
    if usable:
        return "leftover and tentatively usable: both genes seen + outcome field present (see analysis/)"
    if pr.get("has_TACSTD2") in (False, "False", "") and pr.get("peeked_file"):
        return "leftover: processed matrix peeked; TACSTD2/CLDN4 not both present or outcome missing"
    if pr.get("has_outcome_annotation") in (False, "False", ""):
        return "leftover: no per-sample ICI outcome in GEO characteristics"
    if not pr:
        return "leftover relevant but not probed"
    return "leftover: open data present but not a usable TACSTD2/CLDN4 vs ICI-outcome test"


def _fmt_tests(tests):
    if not tests:
        return "_No leftover association tests were run._"
    lines = ["| Cohort | Gene | Outcome | n1 | median1 | n2 | median2 | p | Cliff δ |",
             "|---|---|---|---:|---:|---:|---:|---:|---:|"]
    for t in tests:
        lines.append(
            f"| {t.get('cohort','')} | {t.get('gene','')} | {t.get('outcome','')} | "
            f"{t.get('n1','')} | {t.get('median1','')} | {t.get('n2','')} | "
            f"{t.get('median2','')} | {t.get('p_value','')} | {t.get('cliffs_delta','')} |"
        )
    return "\n".join(lines)


def _writeup(n, n_left_rel, n_already, n_usable, n_tests, n_sig,
             catalog, tests, skipped):
    already = [c for c in catalog if c["leftover_status"] == "already_analysed"]
    left_rel = [c for c in catalog if c["leftover_status"] == "leftover_relevant"]
    usable_rows = [c for c in catalog if c["usable_maybe"] is True or c["usable_maybe"] == True]
    left_table = ["| Accession | n | Category | TACSTD2 | CLDN4 | Outcome keys | Verdict |",
                  "|---|---:|---|---|---|---|---|"]
    for c in left_rel:
        left_table.append(
            f"| {c['accession']} | {c['n_samples']} | {c['study_category']} | "
            f"{c['has_TACSTD2']} | {c['has_CLDN4']} | {c['outcome_keys'] or '—'} | "
            f"{c['verdict']} |"
        )

    skip_txt = "None."
    if skipped:
        skip_txt = "\n".join(f"- {s.get('accession')}: {s.get('reason')}" for s in skipped)

    en = f"""# GEO 2023 leftover lung ICI series: TACSTD2 / CLDN4

Slice `w200` / `results/w200/GEO_2023/`. Independent re-search of GEO for
**2023** human lung immune-checkpoint-inhibitor (ICI) series, minus the one
2023 lung cohort already tested in `fable_geo_2022_2023`.

## 1. What "leftover" means

The earlier slice `fable_geo_2022_2023` searched 2022–2023 together and
**actually tested TACSTD2/CLDN4 vs ICI outcome in only one 2023 lung series**:
**GSE207422** (NSCLC neoadjuvant anti-PD-1 + chemo, n=24 baseline tumours;
directionally lower markers in MPR, all p > 0.2). GSE243238 (2023 acral
melanoma) was a non-lung cross-check.

This slice **does not re-test those two**. Leftover = every other 2023 GSE
returned by an independent human × lung × ICI search. Most leftovers were
listed in the earlier catalog but never probed for gene presence.

Already excluded as analysed:
{chr(10).join(f"- {c['accession']}: {c['already_analysed_note']}" for c in already) or "- (none found in this re-search)"}

## 2. Search

- Database: NCBI GEO via E-utilities `db=gds`, `gse[ETYP]`, *Homo sapiens*,
  PDAT 2023-01-01 … 2023-12-31.
- Union of a broad lung×ICI query plus per-drug and per-histology queries.
- **{n}** unique GSE series dated 2023.
- **{n_already}** already analysed (excluded).
- **{n_left_rel}** leftover series pass the text filter
  (human + lung + ICI + expression).

Files: `search_candidates_2023.csv`, `leftover_classification.csv`,
`leftover_probe.csv`, `leftover_suppl_files.csv`, `leftover_catalog.csv`.

## 3. Leftover relevant series

{chr(10).join(left_table)}

## 4. Association tests on leftovers

Leftover series flagged `usable_maybe` (both genes seen in an open processed
matrix **and** a per-sample outcome field): **{n_usable}**.

Tests run: **{n_tests}**. Tests with p < 0.05: **{n_sig}**.

{_fmt_tests(tests)}

Skipped after a usable_maybe flag:
{skip_txt}

## 5. Honest conclusion

The only leftover 2023 lung series that can be tested is **GSE248378**
(Altorki et al., NCT02904954): post-durvalumab ± SBRT *non-MPR* resected
tumours, n=29, both genes present. GEO itself has no recurrence label.
Recurrence (9 vs 20) was recovered from the paper's public Source Data
Fig. 5d by matching unique ITGAE FPKM values — not by guessing that a
sample title ending in "R" means recurrence (that guess is false: POD23R /
POD26R / POD27R / POD32R did *not* recur).

TACSTD2 is higher in the 9 tumours that later recurred (median FPKM 80.0 vs
42.0; Mann–Whitney p=0.056; Cliff δ −0.46). CLDN4 is in the same direction
and weaker (77.4 vs 39.0; p=0.31). Neither test is p<0.05. The cohort is
**post-treatment and non-MPR only**, so it cannot address baseline
prediction or MPR. It is not a second GSE207422-class baseline
whole-transcriptome ICI-response cohort.

Every other leftover is single-cell T cells, platelets, PBMCs, cell-line
siRNA, ATAC, or a targeted panel without these two genes. Do not treat the
leftover catalog as confirmed biomarker evidence. The TACSTD2 p=0.056
result is hypothesis-generating and selected.

## 6. Limitations

- Text filters can miss a series whose GEO title/summary never says lung or
  ICI; per-drug queries reduce but do not eliminate that risk.
- Gene peeking reads the first column of the smallest processed file < 80 MB.
  A gene present only in a larger or binary (`.rds` / `.RData`) file can be
  missed; those files are recorded in `leftover_suppl_files.csv`.
- No raw-read realignment. Author-processed matrices only.
- Single-cell leftovers were not pseudo-bulked. TACSTD2/CLDN4 are epithelial;
  T-cell / Treg / PBMC series cannot answer the tumour-intrinsic question.
"""

    zh = f"""# GEO 2023 年肺癌 ICI「剩余」系列：TACSTD2 / CLDN4

切片 `w200` / `results/w200/GEO_2023/`。对 GEO **2023** 年人类肺癌免疫检查点
抑制剂（ICI）系列做独立再检索，并减去 `fable_geo_2022_2023` 中已经检验过的
2023 年肺癌队列。

## 1. 「剩余」的定义

先前切片 `fable_geo_2022_2023` 把 2022–2023 放在一起检索，但 **真正对
TACSTD2/CLDN4 与 ICI 结局做了统计检验的 2023 年肺癌系列只有 GSE207422**
（NSCLC 新辅助抗 PD-1 + 化疗，基线肿瘤 n=24；MPR 组标志物更低，但全部
p > 0.2）。GSE243238（2023 肢端黑色素瘤）是非肺交叉验证。

本切片 **不再重复检验这两套**。剩余 = 独立检索得到的其余全部 2023 年 GSE。
其中多数在先前目录里出现过，但从未核查基因是否存在。

已排除：
{chr(10).join(f"- {c['accession']}: {c['already_analysed_note']}" for c in already) or "- （本次再检索未命中已分析系列）"}

## 2. 检索

- NCBI GEO E-utilities，`gse[ETYP]`，人，PDAT 2023-01-01 至 2023-12-31。
- 宽查询 ∪ 逐药 ∪ 逐组织学。
- **{n}** 个 2023 年 GSE。
- **{n_already}** 个已分析（排除）。
- **{n_left_rel}** 个剩余系列通过文本过滤（人 + 肺 + ICI + 表达）。

## 3. 剩余相关系列

见英文表。完整字段在 `leftover_catalog.csv`。

## 4. 剩余系列上的关联检验

`usable_maybe`（开放矩阵中见到两个基因 **且** 有逐样本结局字段）：**{n_usable}**。
已跑检验：**{n_tests}**。p < 0.05：**{n_sig}**。

{_fmt_tests(tests)}

## 5. 诚实结论

唯一能做检验的 2023 年剩余肺系列是 **GSE248378**（Altorki 等，NCT02904954）：
度伐利尤单抗 ± SBRT **治疗后、非 MPR** 切除肿瘤，n=29，两个基因都在。
GEO 没有复发标签。复发（9 vs 20）是用论文公开 Source Data 图 5d 的
ITGAE FPKM 一一对上的，**不是**把样本名末尾的 “R” 当成复发
（POD23R/26R/27R/32R 实际未复发）。

TACSTD2 在随后复发的 9 例更高（中位 FPKM 80.0 vs 42.0；Mann–Whitney
p=0.056；Cliff δ −0.46）。CLDN4 同向更弱（77.4 vs 39.0；p=0.31）。
都不是 p<0.05。队列是 **治疗后且仅非 MPR**，不能回答基线预测或 MPR。
它也不是第二套 GSE207422 级别的基线全转录组 ICI 反应队列。

其余剩余系列是单细胞 T 细胞、血小板、PBMC、细胞系 siRNA、ATAC，
或没有这两个基因的靶向 panel。不要把剩余目录当成已证实的标志物证据。
TACSTD2 p=0.056 只是假设生成、且经过选择。

## 6. 限制

- 标题/摘要从未写 lung 或 ICI 的系列可能漏检。
- 基因探测只读 < 80 MB 最小处理文件的第一列；`.rds`/`.RData` 可能漏检。
- 不做原始读段重比对。
- 不对单细胞剩余系列做伪 bulk：TACSTD2/CLDN4 是上皮基因。
"""
    return en + "\n---\n\n" + zh + "\n"


if __name__ == "__main__":
    main()
