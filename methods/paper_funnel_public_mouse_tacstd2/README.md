# Paper funnel: public mouse Tacstd2

Public mouse integration only (GSE137244 and related KL/KP bulk; GSE154977 / GSE180963 / GSE165641 scRNA). **No private 8KL.**

## Arms

| Arm | Script | Finding |
|---|---|---|
| KL > KP Tacstd2 | `scripts/analyze_kl_kp.py` | `FINDING_KL_KP.md` |
| Tacstd2% vs T_frac (max \|ρ\| grid + locked) | `scripts/analyze_tacstd2_maxrho.py` | `FINDING_TACSTD2_MAXRHO.md` |
| Tacstd2-high epi DEG → TJ | `scripts/analyze_tacstd2_deg_tj.py` | `FINDING_DEG_TJ.md` |

Master write-up: `FINDING.md`.

## Reproduce

```bash
pip install -r methods/paper_funnel_public_mouse_tacstd2/requirements.txt
bash methods/paper_funnel_public_mouse_tacstd2/scripts/download_scrna.sh /tmp/kpkl_10x
python3 methods/paper_funnel_public_mouse_tacstd2/scripts/analyze_kl_kp.py
python3 methods/paper_funnel_public_mouse_tacstd2/scripts/analyze_tacstd2_maxrho.py \
  --data /tmp/kpkl_10x --out methods/paper_funnel_public_mouse_tacstd2
python3 methods/paper_funnel_public_mouse_tacstd2/scripts/analyze_tacstd2_deg_tj.py \
  --data /tmp/kpkl_10x --out methods/paper_funnel_public_mouse_tacstd2
```

Caches for bulk GEO downloads land in `cache/kl_kp/` (gitignored).
