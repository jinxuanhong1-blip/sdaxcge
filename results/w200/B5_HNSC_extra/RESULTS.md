# B5 leftover HNSCC ICI CLDN4

Prat GSE93157 HNSCC (n=5, all nivolumab) cannot test CLDN4: NanoString Immune 730 has no CLDN family. Same for Foy GSE159067 (HTG) and ALPHA GSE190575. NIVACTOR train GSE212549 has Clariom D CLDN4 but no public ICI labels. EGA EGAD50000002506 is controlled FASTQ.

Two leftovers measure CLDN4 and have public ICI labels.

**GSE179730** (Liu 2021; neoadjuvant nivo OCSCC; n=11 pretreatment). Deposited matrix is linear CPM, not log2. Table S2: 3 pathologic responders, 3 stable, 5 progressors. Primary lock = clinical benefit vs progression (6 vs 5). CLDN4 is sparse (6/11 nonzero; median 0). Benefit vs progression: AUC 0.133, exact MWU p=0.045, BH q=0.091 versus locked 3-vs-rest sensitivity. Direction is opposite “high CLDN4 = benefit” and is driven by zeros. TACSTD2 is almost all zero (p=0.727). Detectable-vs-zero RFS/OS (last-observation text, not the contradictory Status column) are null (exact log-rank 1.0 and 0.455).

**GSE212550** (NIVACTOR test; R/M HNSCC ICI monotherapy; n=20 public LTS>18 mo vs STS<6 mo). Clariom D probe TC0700007993.hg.1. CLDN4 medians 3.51 vs 3.54; AUC 0.510; exact p=0.955.

Honest verdict: no reproducible leftover HNSCC ICI CLDN4 signal. Prat cannot be used. Liu is n=11 and sparse. NIVACTOR leftover that measures CLDN4 is null.
