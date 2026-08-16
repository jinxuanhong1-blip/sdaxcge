# Leftover 2024–2026 public human lung ICI / neoadjuvant scRNA

Additive hunt. **Not** the core set (GSE207422, GSE205335, GSE241934, GSE291670, GSE243013, GSE253013, GSE131907).

Public GEO / ArrayExpress only. EGA and dbGaP payloads were not used.

## Question

For each leftover series that has a public processed matrix **and** an epithelial / malignant compartment: patient-level epithelial TACSTD2 and CLDN4 versus T/NK fraction, and versus any public response / PD-L1 label. Honest n / ρ / p. Cell-level tests are not a patient-level claim.

## Reproduce

```bash
python3 methods/scrna_leftover_2426/scripts/01_fetch_geo_metadata.py
python3 methods/scrna_leftover_2426/scripts/02_fetch_more_meta.py
# download usable GEO matrices to /tmp/scrna_leftover_2426 (see 03)
python3 methods/scrna_leftover_2426/scripts/03_score_usable.py
python3 methods/scrna_leftover_2426/scripts/05_emtab_score.py   # remote stream of tumor h5ad
```

## Primary tables

- `tables/inventory.tsv` — one row per candidate series
- `tables/stats.tsv` — tests that exist
- per-series sample / patient TSVs in `tables/`

## Scoring

- Lineage = argmax of mean log1p(CP10k) marker scores (epithelial / T / NK / B / myeloid / fibroblast / endothelial).
- Endpoint is **epithelial** TACSTD2 / CLDN4. A normal-lung-marker malignant split was attempted; it collapsed when AT2/club markers were uniformly low, so malignant-like rows are not used as the claim.
- Spearman ρ at the patient (or sample-as-patient) unit. Eligible: ≥10 epithelial and ≥20 T/NK cells unless noted.
- Response contrasts only when a public label exists.

## Honest takeaway

No leftover series is a powered public lung-ICI scRNA cohort with malignant cells **and** MPR/RECIST labels outside the core set.

- **GSE267108** (n=8 LUAD, treatment-naïve, PD-L1 inferred from published cell counts): epi TACSTD2 vs T/NK ρ=−0.14, p=0.74; CLDN4 ρ=−0.36, p=0.39. PD-L1 pos/high vs neg is not significant (TACSTD2 p=0.11).
- **GSE274595** (n=8 resected NSCLC, HLA/NKEV paper, not ICI-treated): epi TACSTD2 vs T/NK ρ=+0.75, p=0.052 (n=7). Direction is **positive**, opposite a TROP2-high = immune-cold claim. Tissue-only n=2 after T/NK filter.
- **GSE337519**: n=1 paired neoadjuvant chemo-IO; descriptive only.
- **E-MTAB-13526**: Cvejic atlas, extra n, **not ICI**.
- **GSE176021**: GEO processed object is CD3/CD8 T cells (lymphocytes) with MPR labels; no epithelium.
- **GSE186446**: scRNA is T cells from 3 ICB patients. The only GEO expression matrix is bulk `fcount_aggr` (not scRNA).
- **GSE308745**: PBMC; skipped (no epithelium).
- **GSE302113**: scATAC + mtDNA, not RNA.
- **GSE274584 / GSE176022**: not scRNA expression (EV RNA / bulk TCR).
