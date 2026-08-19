# Additive public Visium, CLDN4-only

Two public LUAD Visium sources, run if downloadable:

- **A)** [GSE277206](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE277206) — CytAssist FFPE, never-smoker MIA (n=2)
- **B)** [Zenodo 13337961](https://zenodo.org/records/13337961) — lepidic vs solid ([10.3389/fimmu.2024.1430163](https://doi.org/10.3389/fimmu.2024.1430163))

No private 8-KL. Write-up: `results/visium_public_luad/RESULTS.md`.

```bash
python3 methods/visium_public_luad/download.py
python3 methods/visium_public_luad/analyze.py
```
