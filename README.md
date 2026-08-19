# TISMO ICB pairing: Tacstd2, Cldn4, TJ score

Public recompute of TISMO in-vivo ICB pairs (naive vs anti-PD1 / PD-L1 / CTLA4 / combo), plus Cldn4 and a tight-junction score (Cldn3/4/6/7, Cdh1, F11r, Ocln). Lung / KL / KP / LLC are called out separately. LLC is not KL.

See **[RESULTS.md](RESULTS.md)** for the paper sentence and exact numbers.

```bash
pip install -r scripts/requirements.txt
python3 scripts/download_tismo.py
python3 scripts/analyze_tismo.py
python3 scripts/analyze_gemm.py   # public KP/KL GEO ICB, because TISMO lung is LLC-only
```
