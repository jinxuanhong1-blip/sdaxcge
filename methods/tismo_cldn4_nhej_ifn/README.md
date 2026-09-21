# TISMO ICB: Cldn4 versus NHEJ and IFN scores

Scores Hallmark IFN-γ and GO:0006303 NHEJ on the 64 TISMO ICB slices locked in PR #542.

Tacstd2 stays 49/64. The analyzer exits if recomputed Tacstd2 or Cldn4 mean-deltas disagree with `locked_pairs_pr542.tsv`.

```bash
pip install -r methods/tismo_cldn4_nhej_ifn/requirements.txt
python3 methods/tismo_cldn4_nhej_ifn/download.py
python3 methods/tismo_cldn4_nhej_ifn/analyze.py
```

Narrative and numbers: `results/tismo_cldn4_nhej_ifn/RESULTS.md`.
