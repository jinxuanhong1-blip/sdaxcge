# CLDN4 vs DNA-damage signatures, STING, and HLA

DepMap Public 24Q4 lung cell lines, plus MCLP CCLE RPPA (20221116) for phospho-DDR and STING/HLA proteins. γH2AX (H2AX-pS139) is not in those RPPA panels.

Numbers: `FINDING.md`.

```bash
pip install -r methods/ccle_cldn4_ddr_sting_hla/requirements.txt
python3 methods/ccle_cldn4_ddr_sting_hla/download.py data/ccle_cldn4_ddr_sting_hla
python3 methods/ccle_cldn4_ddr_sting_hla/analyze.py \
  --data data/ccle_cldn4_ddr_sting_hla \
  --outdir methods/ccle_cldn4_ddr_sting_hla
```

`data/` is gitignored. The download script drops the full expression matrix after writing the gene panel.
