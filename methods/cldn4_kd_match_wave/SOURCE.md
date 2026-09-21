# Source of the locked effects

The match wave does not re-estimate fold changes.

| File | Role |
|---|---|
| `locked/gene_effects_ifn_apm.tsv` | IFN_ISG and MHC_APM rows for the four contrasts below |
| `locked/geneset_summary_ifn_apm.tsv` | Locked medians, up-counts, and direction calls |

Copied from Claim C4, PR #88, branch `cursor/claim-c4-public-analogs-0b42`:

- `results/claim_C4/data/processed/gene_effects.tsv`
- `results/claim_C4/data/processed/geneset_summary.tsv`

SHA-256 of the copies in this directory:

- `gene_effects_ifn_apm.tsv` `5481331fb4a73f6aa9fe57f034f3e71f4d5aa1269436a6b522dfd834704e8a25`
- `geneset_summary_ifn_apm.tsv` `7731b66556ddbdba983e9f3569cac4a8c189dfff961ed0322a8e1c7781150f81`

Contrasts kept:

| Accession | Contrast | Class | Why it is in |
|---|---|---|---|
| GSE207704 | T47D KO vs WT | Cancer, breast | True CRISPR CLDN4 KO |
| GSE207704 | MCF-7 KO vs WT | Cancer, breast | True CRISPR CLDN4 KO |
| GSE22493 | SKOV-3 KD | Cancer, ovary | True CLDN4 siRNA; probes map to symbols |
| GSE50927 | baseline KO vs WT | Non-cancer, lung | True Cldn4 KO; uninjured baseline only |

Left out on purpose: GSE50927 VILI-low and VILI-high, and every TACSTD2/TROP2 knockdown (including GSE245459).

`build_concordance.py` drops GSE22493 if IFN/ISG or MHC-I/APM coverage is below 50% of the locked panel. Both panels clear that gate, so the ovarian contrast stays.
