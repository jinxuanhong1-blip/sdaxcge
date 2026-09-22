# Smooth-funnel elevator pitch (EN ∥ ZH)

**Story gap fill only.** Parallel ~120-word pitches for the paper/PPT smooth funnel.
No new matrices, no invented NES/FDR/ρ. Private 8KL / CD34 HIS / KD co-culture numbers are not quoted here.

**Must-mention checklist**

| Element | Present |
|---|---|
| SKB264 + ICI | yes |
| TROP2 label | yes |
| Junction program | yes |
| CLDN4 screen / pin | yes |
| No false TJ#1 public GSEA | yes |

**Provenance (locked public limbs, not re-run here):** paper-funnel PRs #735–#744, #737; CosMx / concordant-4 / GSE137244 / TISMO locks in the public handoff; GSEA honesty lock PR #744 (`any_cohort_highlight_in_top3 = false`).

---

## English (~120 words)

In STK11-cold lung cancer, the clinical hook is **SKB264 plus ICI**: the same **TROP2 label** that marks hard-to-treat disease is the ADC entry point for combination immunotherapy. Public cohorts place TACSTD2/Tacstd2-high tumors in immune-sparse niches with fewer T/NK neighbors, without claiming a sole causal antigen story. The smooth funnel does not jump from label to drug; it routes through an epithelial **junction program** that co-travels with Tacstd2 after KL>KP and Tacstd2-high DEG. Inside that program we **screen and pin CLDN4**, not a vague TJ mean, as the immune-inverse member across patient scRNA, CosMx exclusion, and protein. Critically, public Tacstd2-high GSEA does **not** rank tight-junction, adhesion, or Claudin sets at rank **#1**; we state that and refuse a false TJ#1 headline.

**Word count:** 120

---

## 中文（平行约 120 词；词间空格便于计数，口播可去空格）

在 STK11 冷 肺 癌 中 ， 临床 钩子 是 **SKB264 联合 ICI** ： 同一 **TROP2 标签** 既 标出 难治 亚群 ， 也是 ADC 进入 联合 免疫 入口 。 公共 队列 里 TACSTD2/Tacstd2 高 表达 肿瘤 落在 免疫 稀疏 微环境 ， T/NK 邻居 更 少 ， 但 不 写成 单一 抗原 决定论 。 顺滑 漏斗 不 从 标签 直 跳 到 药 ， 而是 经 与 Tacstd2 共走 的 上皮 **连接 junction 程序** ， 承接 KL>KP 与 Tacstd2-high DEG 。 程序 内 我们 **筛选 并 钉住 CLDN4** ， 而 不是 含糊 的 TJ 均值 ， 作为 免疫 负相关 成员 ， 覆盖 患者 scRNA 、 CosMx 排斥 与 蛋白 层 。 关键 是 ： 公共 Tacstd2-high **GSEA 并 未 把** 紧密 连接 、 黏附 或 Claudin 家族 排 到 **第 1** ； 我们 如实 写明 ， 拒绝 虚构 TJ#1 标题 。

**词数:** 120（空白分隔词；见 `count_words.py`）

---

## Parallel beat map (same order both languages)

1. Clinical hook → **SKB264 + ICI**
2. Disease marker → **TROP2 label**
3. Public phenotype → immune-sparse / fewer T/NK (no fabricated stats in the pitch)
4. Mechanism route → **junction program** with Tacstd2
5. Member priority → **CLDN4 screen / pin** (not mean TJ)
6. Honesty lock → **no false TJ#1** on public GSEA

---

## Explicit non-claims (do not upgrade in talk track)

- Do not say public GSEA puts KEGG TJ / Claudin at rank #1–3 (PR #744).
- Do not merge private 8KL with public mouse Harmony.
- Do not call Visium same-spot correlation “spatial exclusion.”
- Do not invent post-ICI CLDN4 rise from mixed human bulk.
- Surfaceome ranking alone does not nail CLDN4; public pin inside TJ leans on CLDN4 vs structural TJ-15 (CPTAC) plus immune-inverse limbs; private KD is out of this pitch text.

---

## Spoken Chinese (no counting spaces; same content)

在STK11冷肺癌中，临床钩子是SKB264联合ICI：同一TROP2标签既标出难治亚群，也是ADC进入联合免疫入口。公共队列里TACSTD2/Tacstd2高表达肿瘤落在免疫稀疏微环境，T/NK邻居更少，但不写成单一抗原决定论。顺滑漏斗不从标签直跳到药，而是经与Tacstd2共走的上皮连接（junction）程序，承接KL>KP与Tacstd2-high DEG。程序内我们筛选并钉住CLDN4，而不是含糊的TJ均值，作为免疫负相关成员，覆盖患者scRNA、CosMx排斥与蛋白层。关键是：公共Tacstd2-high GSEA并未把紧密连接、黏附或Claudin家族排到第1；我们如实写明，拒绝虚构TJ#1标题。

---

## Word-count check

```bash
python3 story/smooth_funnel/count_words.py
```
