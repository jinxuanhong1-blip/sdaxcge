# B7 — public CLDN4 TEER / barrier literature

Catalog of **public** papers on claudin-4 (CLDN4 / Cldn4) and
transepithelial electrical resistance (TEER / TER) or related barrier
readouts in **lung epithelium** and **cancer**.

Nothing here is a new experiment. Every number in `numbers.tsv` is
copied from a paper or its PMC HTML. If a paper only shows a figure
without printed Ω·cm² or a printed *P* value, the table says so.

## How to use

| File | What it is |
| --- | --- |
| `VERDICT.txt` | One-page finding |
| `papers.tsv` | Papers, PMIDs, systems, whether TEER was measured |
| `numbers.tsv` | Extracted public numbers |
| `search_log.md` | Queries, inclusion / exclusion |
| `extracted/` | Text from OA PDFs used for extraction |

## Inclusion

- CLDN4 (or Cldn4 / claudin-4) **and** a barrier endpoint: TEER/TER,
  Ussing *R*t / *I*sc, dye flux, alveolar fluid clearance, or ECIS *R*b.
- Lung epithelium (primary, airway, alveolar, mouse lung) **or** cancer
  (lung first; other epithelial cancers if they have public TER numbers).
- Foundational MDCK/LLC-PK1 CLDN4 TEER papers are tagged `context_not_lung`.

## Exclusion

- Invented *P* values or Ω values read off figures without a printed number.
- Private / user-only datasets.
- Papers that only mention CLDN4 IHC in tumors with no barrier assay.

## Do not invent *P* values

`p_value_as_published` is either the string printed in the paper
(e.g. `P < 0.05`, `P = 0.003`) or `not_stated`. No recalculation.
