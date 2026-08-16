# Environment for `methods/tcr`

Python templates were checked against **scirpy 0.25.x** (`ir.io.read_10x_vdj`,
`ir.pp.index_chains`, `ir.tl.define_clonotypes`, `ir.tl.clonal_expansion`).

R templates target **immunarch ≥ 0.10** (`repLoad`, `repClonality`,
`repDiversity`, `select_barcodes`).

```bash
# Python (pip)
python3 -m pip install 'scirpy>=0.25,<0.26' 'scanpy>=1.10' 'mudata>=0.3' pyyaml

# R
Rscript -e 'install.packages("immunarch", repos="https://cloud.r-project.org")'
```

Do not put controlled FASTQ in this environment's working tree.
