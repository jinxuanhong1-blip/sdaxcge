Public ELF3 ChIP/CUT&Tag and ELF3-loss transcriptomes: classical NHEJ, cGAS–STING, and interferon.

```bash
pip install -r scripts/elf3_nhej_sting_ifn/requirements.txt
python scripts/elf3_nhej_sting_ifn/analyze.py
```

The script downloads GEO and annotation files into `results/elf3_nhej_sting_ifn/cache/` (gitignored) if they are not already there. Narrative and tables land in `results/elf3_nhej_sting_ifn/`.
