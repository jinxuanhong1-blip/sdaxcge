# HPA lung TACSTD2 (TROP2) and CLDN4 — public IHC/protein tables

中 / EN. Official Human Protein Atlas downloads only. No invented scores.

## 结论 / Verdict

HPA 公开的是**序数 IHC**（not detected / low / medium / high）和少量患者级 TMA 记录，不是 H-score，也不是同一芯上的 TACSTD2+CLDN4 配对，更没有 ICI 结局。

HPA publishes **ordinal IHC** and small TMA patient records. It does not publish H-scores, paired TACSTD2+CLDN4 on the same cores, or immunotherapy labels.

## 正常肺实质 IHC / Normal lung parenchyma IHC

Source: v25.1 `normal_ihc_data.tsv` (`tables/02_…`). Reliability = Enhanced. v23 `normal_tissue.tsv` 相同 (`tables/04_…`).

| Gene | Tissue | Cell type | Level |
| --- | --- | --- | --- |
| TACSTD2 | Lung | alveolar cells | Not detected |
| TACSTD2 | Lung | macrophages | Not detected |
| CLDN4 | Lung | alveolar cells type I | Low |
| CLDN4 | Lung | alveolar cells type II | Not detected |
| CLDN4 | Lung | endothelial cells | Not detected |
| CLDN4 | Lung | macrophages | Not detected |

知识表把 **bronchus / nasopharynx respiratory epithelium** 标为 TACSTD2 Medium、CLDN4 Medium（CLDN4 无 nasopharynx 行）。抗体级 XML 与此一致：肺泡仍为 not detected；支气管纤毛/呼吸上皮为 medium 或 high，定位 cytoplasmic/membranous，quantity 多为 >75%（`tables/06_…`）。

The knowledge table also scores **bronchus (and TACSTD2 nasopharynx) respiratory epithelium as Medium**. Per-antibody XML agrees: alveoli remain not detected; airway epithelium is medium or high, cytoplasmic/membranous.

## 肺癌 TMA IHC / Lung-cancer TMA IHC

Official knowledge counts (`cancer_data.tsv`, `tables/03_…`; identical in v23 pathology):

| Gene | High | Medium | Low | Not detected | n |
| --- | --- | --- | --- | --- | --- |
| TACSTD2 | 0 | 2 | 3 | 5 | 10 |
| CLDN4 | 0 | 9 | 2 | 0 | 11 |

这些计数**只匹配一块抗体**：TACSTD2 = HPA043104；CLDN4 = CAB002610。另外两块 TACSTD2 抗体（HPA055067、CAB072852）在**另一组 12 名患者**上阳性更多（High 2 和 5），不能当成知识表的重复计数（`tables/14_…`）。

Those counts match **one antibody each**. Two other TACSTD2 antibodies stain a **different 12-patient set** more strongly and do not match the knowledge table.

患者级 XML（`tables/07_…`）：

- TACSTD2 与 CLDN4 的 `patient_id` **零重叠**。
- TACSTD2 HPA043104：5 ADC + 5 SQCC。
- CLDN4 CAB002610：6 ADC + 3 SQCC + 1 malignant neoplasm NOS + 1 malignant carcinoid。并非全部 NSCLC。
- 有公开 JPG/TIF URL，无数值 H-score，无 ICI 字段。

Patient IDs do not overlap between genes. CLDN4’s TMA is not NSCLC-only. Image URLs are public; H-scores and ICI labels are not.

## 蛋白（非 IHC）/ Protein that is not IHC

Bulk tissue MS (`ms_tissue.tsv`): lung intensity TACSTD2 = 22139708.2，CLDN4 = 440862.9。三份 lung 重复见 `tables/09_…`。

HPA-hosted CPTAC tumor-vs-normal (`cancer_cptac.tsv`):

- TACSTD2 Lung AC: logFC = −0.302233，adj p = 1.72E-4
- TACSTD2 Lung SQCC: logFC = 0.254722，adj p = 1.00E0
- **CLDN4 lung rows: absent**（文件里 CLDN4 只有 Colon AC / Liver HC / Ovary SC / Pancreatic DAC）

这是治疗初治 CPTAC，不是 ICI 队列。

## RNA 对照（标明非蛋白）/ RNA companion (not protein)

HPA consensus tissue RNA lung nTPM：TACSTD2 74.5，CLDN4 66.7（`tables/12_…`）。肺泡 IHC 对 TACSTD2 为 not detected，因此 **RNA ≠ 肺泡蛋白 IHC**。TCGA RNA–survival 在 HPA 肺癌行全部为 unprognostic（`tables/11_…`），也不是 ICI 结局。

Alveolar TACSTD2 IHC is not detected while lung RNA is present. Do not treat RNA as IHC.

## 未做 / Not done

未把序数改成 H-score，未合并不同抗体患者，未声称配对共定位，未下载 730 MB 的 `proteinatlas.xml.gz`。完整大表留在 gitignored `.cache/`。重跑：`python3 scripts/w200_hpa_lung/download_and_extract.py`。

## 引用 / Citation

Human Protein Atlas v25.1, https://www.proteinatlas.org (CC BY 4.0). Uhlén et al. Tissue-based map of the human proteome. *Science* (2015).
