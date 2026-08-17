# CLDN4-first combinatorial search (public lung scRNA)

Additive slice on top of PR #279 / #271 / #274. Those TACSTD2-primary
combinations are taken as given and are not re-audited.

Patient is the unit. The 8-cohort full-pool is **one row**. The search reports
combinations where **CLDN4 ρ < 0** vs T/NK or CXCL13+ (also B and T/NK
cytotoxicity), with honest n and p.

GSE207422 A3 TACSTD2 is taken as given. CLDN4 on the same 12-patient table is
one honest row. GSE253013 uses the existing extract; the 9 GB RDS is not
downloaded.

```bash
python3 methods/scrna_cldn4_combo/combinatorial_search.py
```

Writeup: [`FINDING.md`](FINDING.md). Tables: `tables/`. Extra figures: `figures/`.
