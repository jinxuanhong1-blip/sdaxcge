# Trajectory / RNA-velocity Methods Playbook
# 轨迹与 RNA velocity 方法手册

> **METHODS ONLY / 仅方法学.** Protocol for *whether* and *how* to reconstruct directed epithelial state-change on **public lung** scRNA-seq. No results, no invented numbers, no decorative UMAP streamlines.
>
> **Hard gate / 硬门槛.** Run Palantir / CellRank 2 / velocity **only** if a biological transition is already real in epithelial state space (AT2 → injury/transitional → tumor). Do **not** run them to decorate a UMAP.
>
> **Axis genes / 轴基因.** `TACSTD2` (TROP2) and `CLDN4` are scored **along** that AT2-to-tumor axis — as gene trends, fate-associated drivers, and co-expression with the transitional program — never as a 2-D arrow on a layout.
>
> **Bilingual / 双语.** English first. Chinese mirror starts at [中文版](#中文版-轨迹与-rna-velocity-方法手册).

---

## 0. Scope and non-goals

**In scope**

- A **go / no-go gate** that decides whether any trajectory method is licensed on this object.
- **Palantir** (Setty et al., 2019) as the default static-snapshot fate model.
- **CellRank 2** (Weiler, Lange, Klein et al., 2024) as the 2024–2026 kernel/estimator framework.
- **RNA-velocity skepticism, 2024–2026**: when velocity is disallowed, when it is optional, and how not to circularly “validate” a UMAP.
- How to place **`TACSTD2` / `CLDN4`** on a real AT2-to-tumor epithelial axis in **public** human lung / LUAD scRNA-seq.

**Out of scope**

- Whole-tissue UMAPs that mix immune, stromal, and epithelial cells.
- Inferring direction from a 2-D embedding.
- Claiming that a velocity streamline *is* a lineage.
- Bulk RNA-seq, Visium-without-deconvolution, blood / BAL immune-only objects.
- Private or unpublished counts. This playbook is written for **public** accessions.

**Central premise.** Trajectory inference assumes gradual, memoryless motion on a phenotypic manifold (CellRank 2). That assumption is only as good as the manifold. In public LUAD, the manifold that can carry a *real* transition is the **alveolar epithelial continuum**: homeostatic AT2 ↔ injury-associated intermediates (DATP / PATS / KRT8<sup>+</sup> alveolar intermediate cells, KACs) ↔ malignant epithelium (Han et al., 2024; Moye et al., 2024; Strunz et al., 2020; Kobayashi et al., 2020; Choi et al., 2020). `CLDN4` is a recurrent marker of that intermediate / tumor-proximal state. `TACSTD2` is an epithelial surface gene that is high in airway epithelium at baseline and is often up along injury / tumor axes — it is **not** an AT2 identity gene. Both genes are *readouts on the axis*, not proof that the axis exists.

---

## 1. Go / no-go gate (read this before installing anything)

Do **not** start Palantir, CellRank, or scVelo until every box in **G1–G6** is true. If any box fails, **stop** and report the failed gate. A failed gate is a valid, publishable methods outcome.

### G1. The object is epithelium, not a tissue UMAP

- Input is a **pre-annotated epithelial subset** (`EPCAM` high, `PTPRC`/`PECAM1`/`COL1A1` low), or a malignant-epithelial subset defined by **CNV** (inferCNV / CopyKAT / SCEVAN) plus epithelial markers.
- Immune, endothelial, and fibroblast cells are **removed before** neighbor graphs used for diffusion / kernels.
- If the public object is immune-dominated (BAL, PBMC, TIL-only), the gate fails. Do not “find” an AT2-to-tumor axis in T cells.

### G2. The biological transition is independently real

At least **two** of the following must already exist *without* trajectory software:

1. **State grammar** recovered by clustering / annotation: AT2 (`SFTPC`, `SFTPB`, `NAPSA`, `LAMP3`) and a transitional / intermediate state (`KRT8`, `KRT19`, `CLDN4`, `CDKN1A`, often `PLAUR`) and a malignant or tumor-like state (CNV-altered, AT2-program down).
2. **Orthogonal prior** from the literature on *this* tissue: AT2 → DATP/PATS/KAC → AT1 *or* AT2 → KAC → KRAS-mutant malignant (Han et al., 2024; Moye et al., 2024).
3. **Composition shift** across a real covariate that is not the embedding: tumor vs uninvolved lung, lesion stage (AAH/AIS/MIA/LUAD), or a time / treatment label in the public metadata.
4. **Spatial or histology support** in the same paper / accession (IHC, CosMx, Visium with epithelial deconvolution) showing KRT8/CLDN4 intermediates at the AT2–tumor interface.

A pretty UMAP with a color gradient from AT2 to tumor is **not** evidence. Chari & Pachter (2023) document that 2-D embeddings systematically distort neighborhoods; Zheng et al. (2023) show that using RNA velocity to “confirm” a UMAP is circular.

### G3. Continuity, not a jump between unrelated epithelia

- AT2, intermediate, and tumor cells must share a **connected k-NN graph in PCA / scVI space** after restriction to epithelium.
- Do not stitch club / basal / ciliated cells onto AT2 just because `TACSTD2` is high in airways. That is a **lineage mix-up**, not a transition.
- If AT2 and tumor form two disconnected components, report **discrete states** (composition, DE, CNV). Do not force a pseudotime.

### G4. Sample size and coverage

- Enough epithelial cells in **each** state you intend to call a vertex (rule of thumb: hundreds per state after QC; do not invent a cutoff — compute it from the object and write it down).
- Multiple donors if the public accession has them. A one-patient “trajectory” is a case report.

### G5. Direction has an external arrow

Pseudotime is an **ordering**, not a clock. You need an external arrow before you interpret “early → late”:

| Arrow | Acceptable source on public lung | Not acceptable |
| --- | --- | --- |
| Palantir early cell | AT2 cells in uninvolved / normal-adjacent lung, or lowest CytoTRACE-2 potency among *non-malignant* AT2 | A cell picked because it sits on the left of the UMAP |
| Terminals | CNV-high malignant; or AT1 if the question is repair | “Tumor cluster” defined only by `TACSTD2` high |
| Experimental time | Rare in public LUAD; use `RealTimeKernel` only if the accession has real time points | Sample collection date, sequencing batch |
| RNA velocity | Only after Section 6 QC; never the sole arrow | Streamlines on UMAP |

### G6. The question is a gene-on-axis question

The scientific question this playbook is licensed to answer is:

> Conditional on a real AT2 → transitional → tumor (or AT1) axis, how do **`TACSTD2` and `CLDN4`** change as functions of Palantir / CellRank progress, and do they associate with fate probability toward the malignant terminal?

Questions this playbook is **not** licensed to answer: “Does the UMAP flow prove AT2 becomes tumor?”; “Is TROP2 a lineage driver because the arrow points at it?”

**If G1–G6 pass:** proceed to Palantir (Section 4), then CellRank 2 (Section 5). Velocity (Section 6) stays **off** unless its own QC passes.  
**If any gate fails:** write the failed gate into the methods, stop, and do not produce streamlines.

```text
                    public lung object
                            |
                     restrict to epithelium
                            |
              G1–G4 real states + continuity?
                     /                \
                   no                  yes
                   |                    |
              STOP (report)        external arrow (G5)?
                                     /            \
                                   no              yes
                                   |                |
                              STOP or          Palantir
                           composition-only         |
                                               CellRank 2
                                          (Pseudotime ± CytoTRACE)
                                                    |
                                      velocity QC (Section 6)?
                                         /                  \
                                       no                    yes
                                       |                      |
                                 do not use            optional VelocityKernel
                                 streamlines           (never UMAP-only)
```

---

## 2. Biological object: AT2-to-tumor epithelial axis

### 2.1 State grammar (use these labels; do not invent new ones without markers)

| State | Positive (RNA) | Negative / caution | Role on the axis |
| --- | --- | --- | --- |
| AT2 | `SFTPC`, `SFTPB`, `NAPSA`, `LAMP3`, `ABCA3` | Low `KRT5`, low CNV | Default **early** / source |
| AT1 | `AGER`, `PDPN`, `CAV1` | — | Repair **terminal** (not tumor) |
| Transitional / AIC / DATP / PATS | `KRT8`, `KRT19`, `CLDN4`, `CDKN1A`, often `PLAUR` | Not a dumpster for “stressed AT2” without these genes | **Waypoint** on both AT2→AT1 and AT2→tumor |
| KAC (KRT8-high AIC) | High `KRT8` plus the transitional set; Han et al. 2024 also note `CDKN2A`, `PLAUR`, `CLDN4` | Do not equate every `KRT8+` cell with KAC | Tumor-proximal intermediate |
| Malignant epithelium | CNV alteration + epithelial RNA; AT2 program down | `TACSTD2` high **alone** is not malignancy | Default **tumor terminal** |
| Club / secretory | `SCGB1A1`, `SCGB3A2`; `TACSTD2` often high | — | **Off-axis.** Do not use as AT2 source |
| Basal | `KRT5`, `KRT15`, `TP63`; `TACSTD2` often high | — | **Off-axis** (airway) |

`TACSTD2` is **airway-high at baseline** (club / basal) and can rise in injured or neoplastic alveolar epithelium. Therefore a `TACSTD2` gradient on a mixed epithelial UMAP is frequently an **airway-vs-alveolar** axis, not AT2-to-tumor. Always pair it with `CLDN4` + `KRT8` + AT2 genes + CNV.

### 2.2 What “along the axis” means for the two genes

Score each gene in **high-dimensional** space, then plot against Palantir pseudotime / CellRank absorption probabilities — not against UMAP-x.

1. **Raw and MAGIC-imputed** expression vs Palantir pseudotime, stratified by lineage (AT1 vs malignant).
2. **Branch-specific GAM trends** (Palantir gene trends; CellRank `gene_trends`).
3. **Lineage-driver correlation** (CellRank `compute_lineage_drivers`) of `TACSTD2` and `CLDN4` with fate probability toward the malignant terminal, next to a pre-specified control panel: AT2 (`SFTPC`), transitional (`KRT8`, `CDKN1A`), airway (`SCGB1A1`, `KRT5`), immune (`PTPRC`).
4. **Co-expression**, not just marginal means: fraction of cells that are `CLDN4`+`KRT8`+ inside the intermediate state; `TACSTD2` as a continuous covariate, not a cluster name.

Pre-specify the gene list. Do not hunt the transcriptome and then declare `TACSTD2` a driver.

### 2.3 Public lung object classes (inputs, not results)

Use an accession only if epithelial cells are actually present and labeled (or labelable):

| Class | Examples (public) | Typical use here | Usual gate failure |
| --- | --- | --- | --- |
| LUAD multi-region epithelium | Han et al., *Nature* 2024 atlas; Kim et al. GSE131907 (2020) | AT2 / AIC / KAC / malignant | None if subsetted to epithelium |
| Normal / HLCA alveolar | Travaglini et al., 2020; Sikkema et al. HLCA, 2023 | AT2 → AT1 prior; `TACSTD2` airway baseline | No tumor terminal |
| Injury / fibrosis intermediates | Strunz et al., 2020; Kobayashi et al., 2020; Choi et al., 2020 | DATP/PATS/Krt8 prior | Mouse; map orthologs (`Tacstd2`, `Cldn4`) |
| Immune-dominated ICI | BAL, blood, TIL-only GEO | **Out of scope** | G1 fails |

Record accession, genome build, gene namespace (symbol vs Ensembl), and whether spliced/unspliced matrices exist (they usually do **not** on processed GEO h5ad).

---

## 3. Inputs, restriction, and graphs

### 3.1 Required inputs

| Item | Why |
| --- | --- |
| Raw counts in `adata.layers["counts"]` | Palantir / CellRank / velocity all need a defined expression layer |
| Log-normalized `adata.X` + HVG + PCA or scVI latent | Neighbor graph |
| `adata.obs` with donor, sample, tissue (tumor vs nLung), histology if present | External arrow, batch |
| Epithelial / CNV labels | G1 |
| Optional: spliced / unspliced layers | Velocity only |
| Optional: real time points | `RealTimeKernel` / moscot |

Keep a frozen copy of the **full** object. All trajectory work happens on `adata_epi`.

### 3.2 Restrict, then recompute the graph

```python
# Schematic — fill paths / keys from the accession; do not invent labels.
import scanpy as sc

EPI_LABELS = ["AT2", "AT1", "AIC", "KAC", "transitional", "malignant", "tumor"]  # TODO: map to this object
adata_epi = adata[adata.obs["cell_class"].isin(EPI_LABELS)].copy()

# Drop off-axis airway if the question is alveolar → tumor
# adata_epi = adata_epi[~adata_epi.obs["cell_class"].isin(["club", "basal", "ciliated"])].copy()

sc.pp.highly_variable_genes(adata_epi, layer="counts", flavor="seurat_v3", n_top_genes=3000)
sc.pp.pca(adata_epi, n_comps=50)
sc.pp.neighbors(adata_epi, n_neighbors=30, n_pcs=30, use_rep="X_pca")  # or "X_scVI"
# UMAP is allowed for *display only*. Never compute neighbors on X_umap.
sc.tl.umap(adata_epi)
```

Recompute HVG / PCA / neighbors **on the epithelial subset**. A graph inherited from the whole-tissue object will leak immune neighborhoods into diffusion components.

### 3.3 Confounders to residualize or stratify *before* diffusion

- **Donor / protocol batch.** If you integrate, run Palantir on the **integrated latent** (scVI / Harmony) but score genes on **non-integrated** expression.
- **Cell cycle.** If a cycle ring dominates DC1, regress `S_score`/`G2M_score` or split cycling vs G0/G1 epithelium.
- **Ambient RNA.** Tumor droplets pick up surfactant and *vice versa*. Recheck `SFTPC` in CNV-high cells after SoupX/CellBender; do not treat residual `SFTPC` as proof of an AT2→tumor continuum.
- **Doublets.** Epithelial–immune doublets create fake intermediates. Run doublet detection **per sample** before subsetting.

### 3.4 Embedding policy (decorative UMAP ban)

| Allowed | Forbidden |
| --- | --- |
| UMAP/PHATE as a *legend* for already-computed states | Neighbor graph on 2-D UMAP |
| Palantir diffusion components / multi-scale space | Velocity streamlines as the primary figure |
| CellRank random-walk visualization on the kernel graph | “The arrow on the UMAP shows AT2 becoming tumor” |
| Gene trends vs pseudotime / fate probability | Using velocity to prove the UMAP is correct (Zheng et al., 2023) |

PHATE (Moon et al., 2019) is preferable to UMAP if you need a 2-D *illustration* of a continuous manifold. It is still not evidence.

---

## 4. Palantir (default static-snapshot method)

Palantir models differentiation as a Markov process on a diffusion-map nearest-neighbor graph, returns a **pseudotime**, **branch probabilities** to user- or automatically-chosen terminals, and an **entropy** (plasticity) (Setty et al., 2019). CellRank 2’s `PseudotimeKernel` is explicitly “inspired by Palantir.” Use Palantir first; feed its pseudotime into CellRank.

### 4.1 When Palantir is the right tool

- Snapshot scRNA-seq (almost all public LUAD).
- One early state (AT2 in nLung) and one or two terminals (malignant; optionally AT1).
- You want gene trends and a plasticity profile along a **pre-justified** axis (G2).

When it is the wrong tool: disconnected clusters; no credible early cell; terminals defined by the gene you intend to test (`TACSTD2`).

### 4.2 Procedure

```python
import palantir

# Diffusion maps on the epithelial PCA / kernel PCA representation
palantir.utils.run_diffusion_maps(adata_epi, n_components=10)          # TODO: inspect eigengap
palantir.utils.determine_multiscale_space(adata_epi)
palantir.utils.run_magic_imputation(adata_epi)                         # trends only; never DE

# Early cell: AT2 in uninvolved lung, nearest to the AT2 centroid in multi-scale space
# Do not pick by UMAP coordinates.
early = pick_early_cell(adata_epi, mask=at2_nlung_mask)                # TODO

# Terminals: pre-specified, not the TACSTD2-high quantile
terminals = {
    "malignant": pick_terminal(adata_epi, mask=cnv_malignant_mask),    # TODO
    # "AT1": pick_terminal(adata_epi, mask=at1_mask),                  # only if AT1 is in-object
}

pr_res = palantir.core.run_palantir(
    adata_epi,
    early_cell=early,
    terminal_states=terminals,
    num_waypoints=500,                                                 # TODO: sensitivity 250/500/1000
    knn=30,
)
# Writes: palantir_pseudotime, palantir_entropy, palantir_fate_probabilities
```

**Early-cell rule.** Among non-malignant AT2, take the cell with maximal mean distance *away* from malignant cells in multi-scale space, or the AT2 centroid’s nearest neighbor in nLung. Document the rule. Run a sensitivity: 10 alternative early cells in the same AT2 neighborhood; fate probabilities should not flip.

**Terminal rule.** Malignant terminal = CNV-high epithelial cells, or the public atlas’s malignant label. **Never** = top decile of `TACSTD2` or `CLDN4`. Those genes are the *estimands*.

### 4.3 Gene trends for `TACSTD2` and `CLDN4`

```python
genes = [
    "TACSTD2", "CLDN4",
    "SFTPC", "SFTPB", "NAPSA",          # AT2
    "KRT8", "KRT19", "CDKN1A", "PLAUR", # transitional / KAC
    "AGER", "PDPN",                     # AT1
    "SCGB1A1", "KRT5",                  # airway off-axis controls
    "EPCAM", "PTPRC",                   # lineage controls
]
trends = palantir.presults.compute_gene_trends(
    adata_epi, gene_list=genes, lineages=["malignant"]  # + "AT1" if present
)
```

Report for each gene: trend shape (monotonic up / down / transient peak), the pseudotime of the peak if transient, and whether the trend is **lineage-specific** (malignant vs AT1). A transient `CLDN4`/`KRT8` peak with later malignant rise is the expected *shape* under a KAC waypoint model; a `TACSTD2` trend that tracks `KRT5`/`SCGB1A1` is an airway contaminant.

MAGIC-imputed values are for visualization of trends only. Statistical tests use **unsmoothed** counts (pseudobulk by donor × state, or a GAM on log-normalized expression with a donor random effect).

### 4.4 Entropy / plasticity

Palantir entropy is high at branch points. On this axis, expect higher entropy in the transitional / KAC compartment than in mature AT2 or CNV-high tumor. Treat entropy as a **descriptive** waypoint score, not as a stemness claim. If you want a potency score, use **CytoTRACE 2** (Kang et al., 2025) as a separate, absolute-potency estimate — and still do not define terminals from it if the estimand is `TACSTD2`.

---

## 5. CellRank 2 (2024–2026 default fate framework)

CellRank 2 splits the problem into **kernels** (build a cell–cell transition matrix from a biological prior) and **estimators** (GPCCA: macrostates, initial/terminal states, fate probabilities, lineage drivers) (Weiler et al., 2024). Downstream analysis does not depend on which kernel built `T`. That is the point: you can swap priors without swapping the fate math.

### 5.1 Kernel choice for public lung epithelium

| Kernel | Prior | Use on this object? |
| --- | --- | --- |
| `PseudotimeKernel` | Palantir (or DPT) pseudotime biases the k-NN graph forward | **Default yes** after Section 4 |
| `CytoTRACEKernel` | Transcriptional diversity / CytoTRACE potency as the arrow | **Optional second view**; prefer CytoTRACE 2 scores if computed |
| `ConnectivityKernel` | Symmetric expression similarity (no direction) | Baseline / combination weight only |
| `RealTimeKernel` | Experimental time points (moscot OT under the hood; Klein, Palla, Lange et al., 2025) | **Only** if the accession has real time |
| `VelocityKernel` | RNA velocity | **Default no** — see Section 6 |
| `PrecomputedKernel` | A `T` you already trust | For sensitivity / published `T` |

**Default combination** for a passing G1–G6 LUAD epithelial object:

```text
T = 0.8 * PseudotimeKernel(Palantir) + 0.2 * ConnectivityKernel
```

Add `CytoTRACEKernel` as a *separate* run, not mixed into the same `T`, so you can see whether potency and Palantir arrows agree. Combine `VelocityKernel` only after Section 6 QC, at a modest weight (e.g. ≤ 0.3), and report the un-combined Palantir kernel as the primary result.

```python
import cellrank as cr
from cellrank.kernels import PseudotimeKernel, ConnectivityKernel, CytoTRACEKernel

pk = PseudotimeKernel(adata_epi, time_key="palantir_pseudotime")
pk.compute_transition_matrix()

ck = ConnectivityKernel(adata_epi)
ck.compute_transition_matrix()

vk = (0.8 * pk) + (0.2 * ck)   # CombinedKernel
# Optional second view:
# cy = CytoTRACEKernel(adata_epi).compute_cytotrace().compute_transition_matrix()
```

### 5.2 GPCCA estimator

```python
from cellrank.estimators import GPCCA

g = GPCCA(vk)
g.fit(cluster_key="cell_class", n_states=[3, 4, 5])   # TODO: choose by macrostate stability
g.predict_terminal_states(method="top_n", n_states=2) # expect malignant ± AT1
g.predict_initial_states()
g.compute_fate_probabilities()
```

**Sanity checks (methods, not results):**

- Initial macrostate overlaps AT2 / nLung AT2 — not club, not cycling immune doublets.
- At least one terminal overlaps CNV-high malignant.
- If GPCCA places a terminal on club/basal, you violated G3 (airway leaked in). Restrict and re-fit.
- Fate probabilities to malignant should increase from AT2 → transitional → malignant **as a rank trend**. If they do not, the kernel is not recovering the biology you claimed in G2; stop rather than retuning until the plot looks right.

### 5.3 Lineage drivers and `TACSTD2` / `CLDN4`

```python
drivers = g.compute_lineage_drivers(
    lineages=["malignant"],
    cluster_key="cell_class",
    clusters=["AT2", "AIC", "KAC", "transitional"],  # TODO: names in this object
)
# Inspect pre-specified rows: TACSTD2, CLDN4, SFTPC, KRT8, CDKN1A, SCGB1A1, KRT5, PTPRC
```

Drivers are **correlations** (or model scores) between gene expression and fate probability. They are not causal. A high `CLDN4`–malignant-fate correlation is consistent with a KAC waypoint; it does not prove CLDN4 pushes cells into tumor.

Fit CellRank gene trends (`cellrank.pl.gene_trends` / GAM models) for the same pre-specified list as Section 4.3. Concordance of Palantir trends and CellRank trends is a **robustness** check. Discordance is a stop sign, not a license to pick the prettier curve.

### 5.4 What CellRank 2 is not

- It does not make a bad manifold good.
- It does not rescue velocity that failed Section 6.
- Random-walk plots are qualitative overviews (Weiler et al., 2024), not statistics.
- Fate probabilities are not fractions of a clone. Public LUAD almost never has lineage tracing.

---

## 6. RNA velocity skepticism (2024–2026) — default OFF

### 6.1 Why the default is off on public LUAD

RNA velocity (La Manno et al., 2018; Bergen et al., 2020) infers a high-dimensional vector from spliced vs unspliced counts. The 2021–2026 literature is a sustained record of **assumption failure**, **workflow fragility**, and **benchmark trade-offs**:

- **Model assumptions are routinely false.** Steady-state / common-splicing assumptions fail; even dynamical scVelo assumes constant rates and well-shaped phase portraits (Bergen et al., 2021; Barile et al., 2021). Erythroid and other systems already showed inverted or empty portraits when rates change.
- **Quantification choices move the answer.** Spliced/unspliced pipelines change velocity (Soneson et al., 2021). Most **public processed** LUAD objects ship a single gene-count matrix — no `spliced`/`unspliced` layers. Without the FASTQ (or a velocyto/alevin-fry/STARsolo recount), velocity is **not computable**. That alone fails the velocity gate on most GEO h5ad files.
- **k-NN smoothing dominates the workflow.** Zheng, Stein-O’Brien, Boukas, Goff & Hansen (2023) show that direction and speed inherit the observed k-NN graph; when that graph is wrong, high- and low-dimensional velocity are wrong. Speed estimates are unreliable except at implausibly low noise. Mapping velocity onto an embedding is “effectively interpolating in the embedding space.”
- **Using velocity to validate a UMAP is circular** (Zheng et al., 2023). Low-dimensional embeddings themselves distort (Chari & Pachter, 2023). Projecting velocity onto UMAP blends orthogonal processes (e.g. cell cycle vs differentiation) (Aivazidis et al. / VeloCycle, 2024).
- **Negative controls produce arrows.** Mature, steady-state populations (e.g. PBMC) should not show coherent flow; many methods invent it (Zheng et al., 2023; Wu, Kong, Liao et al., *Genome Biology* 2026).
- **No method wins all tasks.** Wu et al. (2026) benchmarked 25 RNA-only methods across directional consistency, temporal precision, **negative-control robustness**, and depth stability. Directional performance is **anti-correlated** with negative-control robustness (Spearman ρ ≈ −0.57 in that study). UniTVelo (uni) / veloVI / Pyro-Velocity (m2) are their more balanced recommendations; LatentVelo (std) can look directionally excellent and fail negative controls. **No single method is consistently good.**
- **veloVI** (Gayoso, Weiler et al., 2024) is a better *probabilistic* dynamical model than point-estimate scVelo, not a license to skip QC.
- **Theoretical unraveling** (Gorin, Fang, Chari & Pachter, 2022) remains the right prior: the continuous transcription–splicing–degradation ODE is a severe reduction of real bursting, cell-specific rates, and gene–gene dependence.

For an AT2-to-tumor **cancer** axis the assumptions are worse, not better: rates are not developmental constants; malignant cells are often near a new steady state (exactly the negative-control regime); public LUAD is a snapshot, not a pulse-chase.

### 6.2 Velocity go / no-go (in addition to G1–G6)

Velocity may be computed only if **all** of V1–V5 hold. Otherwise write “RNA velocity not performed” and why.

| ID | Requirement |
| --- | --- |
| V1 | Spliced and unspliced counts exist (or you re-quantify public FASTQs with a documented tool: STARsolo + velocyto, alevin-fry, kallisto|bustools). |
| V2 | Phase portraits of a **pre-specified** gene set (AT2, transitional, cell-cycle) are inspected. Enough genes show induction/repression loops to take the model seriously. Genes with empty or blob portraits are removed, not averaged into a vector. |
| V3 | A **negative-control population** in the same object (e.g. mature T cells, or resting AT1 if they are truly static) does **not** yield a confident global arrow. If everything flows, the method is hallucinating. |
| V4 | High-dimensional velocity is analyzed with CellRank `VelocityKernel` (operates on the expression graph, not the UMAP). UMAP streamlines are supplementary or omitted. |
| V5 | Speed is **not** interpreted (Zheng et al., 2023). Direction is treated as a prior to be combined at low weight with Palantir, not as a discovery engine. |

If you do run a method after 2024, prefer **veloVI** or another method that ranked as *balanced* on direction **and** negative control in Wu et al. (2026). Still report scVelo dynamical as a sensitivity, not as ground truth.

```python
# Only after V1–V5. Schematic.
import scvelo as scv
# scv.pp.filter_and_normalize / moments on spliced+unspliced
# veloVI or scv.tl.recover_dynamics + scv.tl.velocity
from cellrank.kernels import VelocityKernel
velk = VelocityKernel(adata_epi).compute_transition_matrix()
# Optional: T = 0.6 * pk + 0.2 * ck + 0.2 * velk
# Primary text result remains the Palantir PseudotimeKernel.
```

### 6.3 Forbidden sentences

Do not write any of the following in a paper that followed this playbook:

- “RNA velocity confirmed the UMAP trajectory from AT2 to tumor.”
- “Cells are moving rapidly toward the TACSTD2-high state.” (speed + gene-defined terminal)
- “Streamlines demonstrate that KACs become LUAD.”
- “We used velocyto/scVelo with default parameters on the tissue UMAP.”

---

## 7. Putting `TACSTD2` and `CLDN4` on the axis

### 7.1 Estimands (pre-specify)

1. **Marginal trend.** Is `CLDN4` (and separately `TACSTD2`) a function of Palantir pseudotime on the malignant lineage? Transient, late, or flat?
2. **Fate association.** Spearman / CellRank driver score of each gene vs malignant fate probability inside **non-malignant** epithelium only (AT2 + transitional). Including the tumor terminal in the correlation double-counts the label.
3. **Waypoint specificity.** Does `CLDN4` peak with `KRT8`/`CDKN1A` (KAC-like) while `SFTPC` falls? Does `TACSTD2` follow that waypoint or the airway controls?
4. **Off-axis rejection.** Association of `TACSTD2` with `KRT5`/`SCGB1A1` after restricting to alveolar-labeled cells. If it survives, you are not on an AT2-to-tumor axis.

### 7.2 Statistics that do not use the UMAP

- **Donor-level:** for each donor, median gene expression in AT2 vs transitional vs malignant; Wilcoxon or mixed model across donors. This is the primary *test*. Trends in Section 4–5 are the primary *description*.
- **Within-donor rank:** correlation of gene vs Palantir pseudotime inside one donor, then meta-analyze (Stouffer / Fisher on ranks). Protects against one donor’s batch.
- **Permutation:** shuffle malignant fate labels within donor; recompute driver scores for the pre-specified genes.
- **Do not** test every gene and then highlight `TACSTD2`.

### 7.3 Protein and spatial (public, if present)

If the accession has CITE-seq, the TROP2 antibody is usually **absent** (immune ADT panels). Do not impute TROP2 protein. If a public spatial assay measures `TACSTD2`/`CLDN4` RNA or TROP2/CLDN4 protein, use it to test **interface localization** (KAC between AT2 and tumor), not to draw velocity.

### 7.4 Interpretation language

Allowed: “On a Palantir/CellRank axis whose early state is nLung AT2 and whose terminal is CNV-high malignant, `CLDN4` expression is associated with intermediate / malignant fate probability, consistent with published KAC / DATP markers.”

Not allowed: “Trajectory analysis showed that TACSTD2 drives AT2 transformation.”

---

## 8. Default pipeline (public lung, gates passed)

1. Download the public object; record accession, license, and gene IDs.
2. QC at the **sample** level (ambient, doublets, MT). Keep counts.
3. Annotate broadly; **restrict to epithelium**; recompute HVG/PCA/neighbors (Section 3).
4. Apply G1–G6. Write a one-paragraph gate report. Stop if failed.
5. Palantir with pre-specified early AT2 and CNV-malignant terminals (Section 4).
6. CellRank 2 `PseudotimeKernel` ± `ConnectivityKernel`; GPCCA; fate probabilities (Section 5).
7. Optional: CytoTRACE 2 as a second arrow; do not mix until compared.
8. Score pre-specified genes (`TACSTD2`, `CLDN4`, controls) as trends + donor-level tests (Section 7).
9. Velocity: skip unless V1–V5 pass (Section 6).
10. Figures: (i) state composition by tissue; (ii) diffusion / Palantir space, not a velocity-UMAP; (iii) gene trends vs pseudotime; (iv) fate-probability violin by state; (v) driver table for the pre-specified genes. UMAP, if shown, is a legend.

### 8.1 Software versions to pin (2024–2026 stack)

Pin exact versions in an `environment.yml`. Suggested lower bounds, not results:

| Package | Role | Pin ≥ |
| --- | --- | --- |
| `scanpy` | graphs, subset | 1.10 |
| `anndata` | I/O | 0.10 |
| `palantir` | pseudotime, trends | 1.3 |
| `cellrank` | kernels, GPCCA | 2.0 (CellRank 2) |
| `scvelo` | velocity I/O / dynamical (optional) | 0.3 |
| `velovi` | optional velocity model | 0.3 |
| `cytotrace2` / companion | optional potency | as published 2025 |
| `infercnvpy` or R `infercnv` | malignant call | current |

### 8.2 Reproducibility

- Seed neighbor graphs and Palantir waypoints; run two seeds.
- Sensitivity: `n_neighbors` ∈ {15, 30, 50}; Palantir `num_waypoints` ∈ {250, 500, 1000}; CellRank `n_states` around the eigengap.
- Leave-one-donor-out: does the malignant terminal still land on CNV-high cells?
- Do not tune these until `TACSTD2` “looks right.”

---

## 9. Methods-paragraph template (fill-in)

```text
Epithelial cells were isolated from the public scRNA-seq object <ACCESSION>
using <EPCAM / annotation / CNV rule>. Immune, endothelial, and mesenchymal
cells were excluded before graph construction. Highly variable genes, PCA
(or scVI), and a k-NN graph were recomputed on the epithelial subset.
We applied a pre-specified go/no-go gate requiring (i) recovery of AT2,
transitional/KRT8-CLDN4, and malignant states without trajectory software,
(ii) a connected epithelial graph, and (iii) an external direction
(uninvolved AT2 → CNV-high malignant). <If gate failed: trajectory
inference was not performed; we report composition and DE only.>

Palantir (Setty et al., 2019) was run on diffusion maps of the epithelial
graph with a documented AT2 early cell and pre-specified terminals that
were not defined by TACSTD2 or CLDN4. Gene trends for a pre-specified
panel (TACSTD2, CLDN4, AT2, transitional, airway, and lineage controls)
were fit on MAGIC-imputed expression for visualization and tested on
unsmoothed, donor-stratified expression.

CellRank 2 (Weiler et al., 2024) used a PseudotimeKernel on the Palantir
pseudotime, optionally combined with a ConnectivityKernel. GPCCA estimated
macrostates, fate probabilities, and lineage-driver scores. RNA velocity
was <not performed: no spliced/unspliced counts | not performed: failed
phase-portrait / negative-control QC | performed with veloVI and used
only as a low-weight kernel>. Velocity streamlines on UMAP were not used
as evidence of an AT2-to-tumor transition (Zheng et al., 2023;
Chari & Pachter, 2023; Wu et al., 2026).
```

---

## 10. References

1. Setty M, et al. Characterization of cell fate probabilities in single-cell data with Palantir. *Nat Biotechnol* 37, 451–460 (2019). https://doi.org/10.1038/s41587-019-0068-4
2. Weiler P, Lange M, Klein M, et al. CellRank 2: unified fate mapping in multiview single-cell data. *Nat Methods* 21, 1196–1205 (2024). https://doi.org/10.1038/s41592-024-02303-9
3. Lange M, et al. CellRank for directed single-cell fate mapping. *Nat Methods* 19, 159–170 (2022).
4. La Manno G, et al. RNA velocity of single cells. *Nature* 560, 494–498 (2018).
5. Bergen V, Lange M, Peidli S, Wolf FA, Theis FJ. Generalizing RNA velocity to transient cell states through dynamical modeling. *Nat Biotechnol* 38, 1408–1414 (2020).
6. Bergen V, Soldatov RA, Kharchenko PV, Theis FJ. RNA velocity—current challenges and future perspectives. *Mol Syst Biol* 17, e10282 (2021).
7. Gorin G, Fang M, Chari T, Pachter L. RNA velocity unraveled. *PLoS Comput Biol* 18, e1010492 (2022).
8. Zheng SC, Stein-O’Brien G, Boukas L, Goff LA, Hansen KD. Pumping the brakes on RNA velocity by understanding and interpreting RNA velocity estimates. *Genome Biol* 24, 246 (2023). https://doi.org/10.1186/s13059-023-03065-x
9. Soneson C, Srivastava A, Patro R, Stadler MB. Preprocessing choices affect RNA velocity results for droplet scRNA-seq data. *PLoS Comput Biol* 17, e1008585 (2021).
10. Gayoso A, Weiler P, et al. Deep generative modeling of transcriptional dynamics for RNA velocity analysis in single cells (veloVI). *Nat Methods* 21, 50–59 (2024).
11. Wu Y, Kong C, Liao X, et al. Comprehensive benchmarking of RNA velocity methods across single-cell datasets. *Genome Biol* 27, 242 (2026). https://doi.org/10.1186/s13059-026-04182-z
12. Luo Y, et al. Benchmarking RNA velocity methods across 17 independent studies. *Cell Rep Methods* 6, 101367 (2026).
13. Chari T, Pachter L. The specious art of single-cell genomics. *PLoS Comput Biol* 19, e1011288 (2023).
14. Kang M, Gulati GS, Brown EL, Qi Z, et al. Improved reconstruction of single-cell developmental potential with CytoTRACE 2. *Nat Methods* 22, 2258–2263 (2025).
15. Gulati GS, et al. Single-cell transcriptional diversity is a hallmark of developmental potential. *Science* 367, 405–411 (2020).
16. Klein D, Palla G, Lange M, et al. Mapping cells through time and space with moscot. *Nature* 638, 1065–1075 (2025).
17. Han G, Sinjab A, Rahal Z, et al. An atlas of epithelial cell states and plasticity in lung adenocarcinoma. *Nature* 627, 656–663 (2024).
18. Moye AL, et al. Early-stage lung cancer is driven by a transitional cell state dependent on a KRAS-ITGA3-SRC axis. *EMBO J* (2024). https://doi.org/10.1038/s44318-024-00113-5
19. Strunz M, et al. Alveolar regeneration through a Krt8+ transitional stem cell state. *Nat Commun* 11, 3559 (2020).
20. Kobayashi Y, et al. Persistence of a regeneration-associated, transitional alveolar epithelial cell state in pulmonary fibrosis. *Nat Cell Biol* 22, 934–946 (2020).
21. Choi J, et al. Inflammatory signals induce AT2 cell-derived damage-associated transient progenitors. *Cell Stem Cell* 27, 366–382.e7 (2020).
22. Travaglini KJ, et al. A molecular cell atlas of the human lung from single-cell RNA sequencing. *Nature* 587, 619–625 (2020).
23. Sikkema L, et al. An integrated cell atlas of the lung in health and disease. *Nat Med* 29, 1563–1577 (2023).
24. Kim N, et al. Single-cell RNA sequencing demonstrates the molecular and cellular reprogramming of metastatic lung adenocarcinoma. *Nat Commun* 11, 2285 (2020). (GSE131907)
25. Moon KR, et al. Visualizing structure and transitions in high-dimensional biological data (PHATE). *Nat Biotechnol* 37, 1482–1492 (2019).
26. van Dijk D, et al. Recovering gene interactions from single-cell data using data diffusion (MAGIC). *Cell* 174, 716–729.e27 (2018).
27. Barile M, et al. Coordinated changes in gene expression kinetics underlie both mouse and human erythroid maturation. *Genome Biol* 22, 197 (2021).
28. Aivazidis A, et al. Statistical inference with a manifold-constrained RNA velocity model uncovers cell cycle speed modulations. *Nat Methods* (2024). https://doi.org/10.1038/s41592-024-02471-8

---

# 中文版 · 轨迹与 RNA velocity 方法手册

> **仅方法学。** 规定在**公开肺**单细胞转录组上，*是否*以及*如何*重建有方向的上皮状态转变。无结果、无编造数字、无装饰性 UMAP 流线。
>
> **硬门槛。** 仅当上皮状态空间中已存在真实生物学转变（AT2 → 损伤/过渡态 → 肿瘤）时，才运行 Palantir / CellRank 2 / velocity。禁止用这些方法给 UMAP 做装饰。
>
> **轴基因。** `TACSTD2`（TROP2）与 `CLDN4` 沿着 AT2→肿瘤轴计分（基因趋势、命运相关 driver、与过渡程序共表达），而不是画成二维箭头。

---

## 0. 范围与非目标

**范围内：** 决定是否允许做轨迹的 **go/no-go 门槛**；**Palantir**（Setty 等, 2019）作为静态快照的默认命运模型；**CellRank 2**（Weiler、Lange、Klein 等, 2024）作为 2024–2026 的 kernel/estimator 框架；**2024–2026 年对 RNA velocity 的怀疑论**；以及如何在公开人肺 / LUAD 数据上把 `TACSTD2`/`CLDN4` 放到真实的 AT2→肿瘤轴上。

**范围外：** 免疫+间质+上皮混在一起的组织 UMAP；从二维嵌入推断方向；把 velocity 流线当成谱系；bulk、未反卷积的 Visium、血液/BAL 纯免疫对象；非公开数据。

**中心前提。** 轨迹推断假定细胞在表型流形上做渐进、无记忆的运动（CellRank 2）。流形不可靠，推断就不可靠。公开 LUAD 中能够承载*真实*转变的流形是**肺泡上皮连续体**：稳态 AT2 ↔ 损伤相关中间态（DATP / PATS / KRT8⁺ 肺泡中间细胞 KAC）↔ 恶性上皮（Han 等, 2024；Moye 等, 2024；Strunz / Kobayashi / Choi, 2020）。`CLDN4` 是该中间态 / 近肿瘤态的反复出现标记。`TACSTD2` 在气道上皮基线就高，损伤/肿瘤轴上常再升高——它**不是** AT2 身份基因。两基因是轴上的*读出*，不是轴存在的*证明*。

---

## 1. Go / no-go 门槛（装软件之前先读）

在 **G1–G6** 全部成立之前，不要启动 Palantir、CellRank 或 scVelo。任一失败即**停止**并报告失败门槛。失败本身是可发表的方法学结果。

### G1. 对象是上皮，不是组织 UMAP

输入必须是预注释的上皮子集（`EPCAM` 高，`PTPRC`/`PECAM1`/`COL1A1` 低），或由 CNV（inferCNV / CopyKAT / SCEVAN）加上皮标记定义的恶性上皮。免疫、内皮、成纤维细胞必须在用于扩散 / kernel 的邻居图**之前**去掉。若公开对象以免疫为主（BAL、PBMC、仅 TIL），门槛失败。不要在 T 细胞里“发现” AT2→肿瘤轴。

### G2. 生物学转变在轨迹软件之外已经成立

以下至少满足 **两项**，且不依赖轨迹软件：

1. **状态语法**已由聚类/注释恢复：AT2（`SFTPC`、`SFTPB`、`NAPSA`、`LAMP3`）、过渡/中间态（`KRT8`、`KRT19`、`CLDN4`、`CDKN1A`，常有 `PLAUR`）、恶性或肿瘤样态（CNV 改变、AT2 程序下降）。
2. **组织先验：** AT2 → DATP/PATS/KAC → AT1，或 AT2 → KAC → KRAS 突变恶性（Han 等, 2024；Moye 等, 2024）。
3. **组成随真实协变量变化**（肿瘤 vs 未受累肺、病变阶段、公开 metadata 中的时间/治疗），而不是随嵌入变化。
4. **同一 accession / 论文中的空间或组织学支持**（IHC、CosMx、经上皮反卷积的 Visium）显示 KRT8/CLDN4 中间态位于 AT2–肿瘤界面。

从 AT2 到肿瘤的 UMAP 颜色渐变**不是**证据。Chari & Pachter（2023）表明二维嵌入会系统扭曲邻域；Zheng 等（2023）表明用 RNA velocity“证实”UMAP 是循环论证。

### G3. 连续体，而不是把无关上皮缝在一起

AT2、中间态、肿瘤细胞在限制到上皮后的 PCA / scVI 空间中必须共享**连通的 k-NN 图**。不要因为气道里 `TACSTD2` 高，就把 club / 基底 / 纤毛细胞缝到 AT2 上——那是**谱系混杂**，不是转变。若 AT2 与肿瘤是两个连通分量，只报告离散状态（组成、DE、CNV），不要强行伪时间。

### G4. 样本量与覆盖

每个打算当作顶点的状态都要有足够上皮细胞（经验上 QC 后每态数百个；不要编造阈值，从对象计算并写明）。公开数据有多个供体就用多个供体。单病人“轨迹”只是病例报告。

### G5. 方向必须有外部箭头

伪时间是**排序**，不是时钟。解释“早 → 晚”之前需要外部箭头：

| 箭头 | 公开肺数据上可接受的来源 | 不可接受 |
| --- | --- | --- |
| Palantir 起点 | 未受累/癌旁 AT2，或非恶性 AT2 中 CytoTRACE 2 潜能最低者 | 因为在 UMAP 左边而挑的细胞 |
| 终点 | CNV 高的恶性；或修复问题中的 AT1 | 仅因 `TACSTD2` 高而叫“肿瘤簇” |
| 实验时间 | 公开 LUAD 中罕见；仅当 accession 真有时间点时用 `RealTimeKernel` | 取样日期、测序批次 |
| RNA velocity | 仅通过第 6 节 QC 后；绝不能当唯一箭头 | UMAP 上的流线 |

### G6. 问题必须是“轴上的基因”问题

本手册被授权回答的问题是：

> 在真实的 AT2 → 过渡态 → 肿瘤（或 AT1）轴成立的前提下，`TACSTD2` 与 `CLDN4` 如何随 Palantir / CellRank 进程变化，以及它们是否与朝向恶性终点的命运概率相关？

未被授权的问题：“UMAP 流是否证明 AT2 变成肿瘤？”；“箭头指向 TROP2，所以 TROP2 是谱系驱动因子？”

**G1–G6 通过：** 进入 Palantir（第 4 节），再进入 CellRank 2（第 5 节）。Velocity（第 6 节）默认关闭。  
**任一失败：** 把失败门槛写入方法，停止，不产流线。

---

## 2. 生物学对象：AT2→肿瘤上皮轴

### 2.1 状态语法（使用这些标签；没有标记就不要发明新名字）

| 状态 | RNA 阳性 | 阴性 / 注意 | 轴上角色 |
| --- | --- | --- | --- |
| AT2 | `SFTPC`、`SFTPB`、`NAPSA`、`LAMP3`、`ABCA3` | 低 `KRT5`、低 CNV | 默认**起点** |
| AT1 | `AGER`、`PDPN`、`CAV1` | — | 修复**终点**（不是肿瘤） |
| 过渡 / AIC / DATP / PATS | `KRT8`、`KRT19`、`CLDN4`、`CDKN1A`，常有 `PLAUR` | 不要把所有“应激 AT2”塞进来 | AT2→AT1 与 AT2→肿瘤的**途经点** |
| KAC（KRT8 高 AIC） | 高 `KRT8` + 过渡基因集；Han 等 2024 还提到 `CDKN2A`、`PLAUR`、`CLDN4` | 不是每个 `KRT8+` 都是 KAC | 近肿瘤中间态 |
| 恶性上皮 | CNV 改变 + 上皮 RNA；AT2 程序下降 | **单独**的 `TACSTD2` 高不是恶性 | 默认**肿瘤终点** |
| Club / 分泌 | `SCGB1A1`、`SCGB3A2`；`TACSTD2` 常高 | — | **轴外**。不要当 AT2 起点 |
| 基底 | `KRT5`、`KRT15`、`TP63`；`TACSTD2` 常高 | — | **轴外**（气道） |

`TACSTD2` **基线在气道高**（club / 基底），在损伤或肿瘤性肺泡上皮也可升高。因此混合上皮 UMAP 上的 `TACSTD2` 梯度常常是**气道 vs 肺泡**轴，不是 AT2→肿瘤。必须与 `CLDN4` + `KRT8` + AT2 基因 + CNV 联用。

### 2.2 两基因“沿轴”的含义

在**高维**空间计分，再对 Palantir 伪时间 / CellRank 吸收概率作图——不要对 UMAP-x 作图。

1. 原始与 MAGIC 插补表达对伪时间，按谱系分层（AT1 vs 恶性）。
2. 分支特异 GAM 趋势（Palantir gene trends；CellRank `gene_trends`）。
3. 与朝向恶性终点的命运概率的 **lineage-driver** 相关，对照组预先指定：AT2（`SFTPC`）、过渡（`KRT8`、`CDKN1A`）、气道（`SCGB1A1`、`KRT5`）、免疫（`PTPRC`）。
4. **共表达**，不只是边际均值：中间态内 `CLDN4`+`KRT8`+ 细胞比例；`TACSTD2` 作为连续协变量，不当簇名。

预先指定基因列表。禁止全转录组搜完再宣布 `TACSTD2` 是 driver。

### 2.3 公开肺对象类别（输入，不是结果）

仅当上皮细胞确实存在且可标注时才使用该 accession：

| 类别 | 公开例子 | 此处典型用途 | 常见门槛失败 |
| --- | --- | --- | --- |
| LUAD 多区域上皮 | Han 等 *Nature* 2024；Kim 等 GSE131907（2020） | AT2 / AIC / KAC / 恶性 | 限制到上皮后通常可通过 |
| 正常 / HLCA 肺泡 | Travaglini 等 2020；Sikkema 等 HLCA 2023 | AT2→AT1 先验；`TACSTD2` 气道基线 | 无肿瘤终点 |
| 损伤 / 纤维化中间态 | Strunz / Kobayashi / Choi, 2020 | DATP/PATS/Krt8 先验 | 小鼠；需直系同源（`Tacstd2`、`Cldn4`） |
| 免疫为主的 ICI | BAL、血液、仅 TIL 的 GEO | **范围外** | G1 失败 |

记录 accession、基因组版本、基因命名空间，以及是否存在 spliced/unspliced（GEO 上已处理的 h5ad **通常没有**）。

---

## 3. 输入、限制与图

### 3.1 必需输入

原始 counts 放在 `adata.layers["counts"]`；log 归一化的 `adata.X` + HVG + PCA 或 scVI；含供体/样本/组织（肿瘤 vs nLung）的 `obs`；上皮 / CNV 标签。可选：spliced/unspliced（仅 velocity）；真实时间点（`RealTimeKernel` / moscot）。完整对象冻结备份；轨迹只在 `adata_epi` 上做。

### 3.2 先限制，再重算图

在上皮子集上**重新**计算 HVG / PCA / neighbors。从全组织对象继承的图会把免疫邻域漏进扩散分量。UMAP 只用于展示，**禁止**在 `X_umap` 上算邻居。代码示意见英文第 3.2 节。

### 3.3 扩散之前要分层或残差掉的混杂

- **供体 / 实验批次。** 若做整合，Palantir 跑在整合潜空间（scVI / Harmony），基因计分必须用**未整合**表达。
- **细胞周期。** 若 DC1 被周期环主导，回归 `S_score`/`G2M_score` 或拆分周期 vs G0/G1 上皮。
- **环境 RNA。** 肿瘤液滴会沾表面活性蛋白，反之亦然。SoupX/CellBender 后再检查 CNV 高细胞里的 `SFTPC`；残留 `SFTPC` 不是 AT2→肿瘤连续体的证明。
- **双细胞。** 上皮–免疫双细胞会制造假中间态。子集化之前按**样本**做双细胞检测。

### 3.4 嵌入政策（禁止装饰性 UMAP）

允许：用 UMAP/PHATE 当*图例*；Palantir 扩散分量；CellRank 在 kernel 图上的随机游走示意；基因趋势对伪时间 / 命运概率。

禁止：在二维 UMAP 上建邻居图；把 velocity 流线当主图；写“UMAP 上的箭头表明 AT2 变成肿瘤”；用 velocity 证明 UMAP 正确（Zheng 等, 2023）。

若需要二维*示意*连续流形，PHATE（Moon 等, 2019）优于 UMAP。它仍然不是证据。

---

## 4. Palantir（静态快照默认方法）

Palantir 在扩散图近邻图上把分化建成 Markov 过程，输出**伪时间**、到终点的**分支概率**，以及作为可塑性的**熵**（Setty 等, 2019）。CellRank 2 的 `PseudotimeKernel` 明确“受 Palantir 启发”。先跑 Palantir，再把它的伪时间喂给 CellRank。

### 4.1 适用与不适用

适用：快照 scRNA-seq（几乎所有公开 LUAD）；一个起点（nLung 的 AT2）和一或两个终点（恶性；可选 AT1）；要在**事先成立**的轴上（G2）看基因趋势和可塑性。

不适用：不连通的簇；没有可信起点；用你打算检验的基因（`TACSTD2`）定义终点。

### 4.2 流程要点

在上皮 PCA 上算扩散图与多尺度空间；MAGIC 仅用于趋势可视化，**不做 DE**。起点：未受累肺 AT2、在多尺度空间中靠近 AT2 质心的细胞——**不要按 UMAP 坐标挑**。终点：CNV 高恶性上皮，或公开图谱的恶性标签。**禁止**用 `TACSTD2`/`CLDN4` 的高分位数当终点——它们是*待估量*。

起点规则写进方法；在同一 AT2 邻域换 10 个备选起点做敏感性，命运概率不应翻转。API 示意见英文第 4.2 节。

### 4.3 `TACSTD2` 与 `CLDN4` 的基因趋势

预指定列表：`TACSTD2`、`CLDN4`、AT2（`SFTPC`、`SFTPB`、`NAPSA`）、过渡/KAC（`KRT8`、`KRT19`、`CDKN1A`、`PLAUR`）、AT1（`AGER`、`PDPN`）、气道对照（`SCGB1A1`、`KRT5`）、谱系对照（`EPCAM`、`PTPRC`）。

对每个基因报告：趋势形状（单调升/降/瞬时峰）、若有峰则其伪时间、是否谱系特异（恶性 vs AT1）。KAC 途经点模型下，期望看到 `CLDN4`/`KRT8` 瞬时峰而后恶性上升；若 `TACSTD2` 跟着 `KRT5`/`SCGB1A1` 走，则是气道污染。

MAGIC 插补只用于看趋势。统计检验用**未平滑** counts（供体 × 状态的 pseudobulk，或带供体随机效应的 GAM）。

### 4.4 熵 / 可塑性

Palantir 熵在分支点高。此轴上，过渡/KAC 的熵应高于成熟 AT2 或 CNV 高肿瘤。把熵当作**描述性**途经点分数，不要当成干性声明。若要潜能分数，用 **CytoTRACE 2**（Kang 等, 2025）作为独立的绝对潜能估计——若待估量是 `TACSTD2`，仍不要用它定义终点。

---

## 5. CellRank 2（2024–2026 默认命运框架）

CellRank 2 把问题拆成 **kernel**（由生物学先验构建细胞–细胞转移矩阵）和 **estimator**（GPCCA：宏观态、起止态、命运概率、谱系 driver）（Weiler 等, 2024）。下游不依赖谁构建了 `T`：可以换先验而不换命运数学。

### 5.1 公开肺上皮的 kernel 选择

| Kernel | 先验 | 此对象上用不用 |
| --- | --- | --- |
| `PseudotimeKernel` | Palantir（或 DPT）伪时间把 k-NN 边向前偏置 | **默认用**（第 4 节之后） |
| `CytoTRACEKernel` | 转录多样性 / CytoTRACE 潜能当箭头 | **可选第二视角**；若已算 CytoTRACE 2 则优先用其分数 |
| `ConnectivityKernel` | 对称表达相似（无方向） | 仅作基线 / 组合权重 |
| `RealTimeKernel` | 实验时间点（底层 moscot OT；Klein、Palla、Lange 等, 2025） | **仅当** accession 有真实时间 |
| `VelocityKernel` | RNA velocity | **默认不用**——见第 6 节 |
| `PrecomputedKernel` | 你已信任的 `T` | 敏感性 / 已发表的 `T` |

**默认组合**（G1–G6 通过的 LUAD 上皮对象）：

```text
T = 0.8 * PseudotimeKernel(Palantir) + 0.2 * ConnectivityKernel
```

`CytoTRACEKernel` 作为*单独*一次运行，不要先混进同一个 `T`，以便比较潜能箭头与 Palantir 箭头是否一致。`VelocityKernel` 仅在第 6 节 QC 通过后以低权重（例如 ≤ 0.3）加入，正文主结果仍是未混合的 Palantir kernel。

### 5.2 GPCCA

拟合宏观态，预测终点 / 起点，计算命运概率。方法学核查（不是结果）：初始宏观态应落在 AT2 / nLung AT2，而不是 club 或周期免疫双细胞；至少一个终点落在 CNV 高恶性；若终点落在 club/基底，说明违反了 G3（气道漏进来），收紧子集后重拟合。朝向恶性的命运概率应沿 AT2 → 过渡 → 恶性呈**秩趋势**。若没有，说明 kernel 没有恢复你在 G2 声称的生物学——停止，而不是调参到图好看为止。

### 5.3 谱系 driver 与 `TACSTD2` / `CLDN4`

`compute_lineage_drivers` 给出表达与命运概率的**相关**（或模型分数），不是因果。`CLDN4` 与恶性命运高相关，与 KAC 途经点模型一致；并不证明 CLDN4 把细胞推进肿瘤。

对第 4.3 节同一预指定列表拟合 CellRank 基因趋势。Palantir 与 CellRank 趋势一致是**稳健性**检查；不一致是停车信号，不是挑选更好看曲线的许可证。

### 5.4 CellRank 2 不是什么

它不能把坏流形变好；不能拯救未过第 6 节的 velocity；随机游走图是定性概览（Weiler 等, 2024），不是统计量；命运概率不是克隆比例。公开 LUAD 几乎没有谱系示踪。

---

## 6. RNA velocity 怀疑论（2024–2026）——默认关闭

### 6.1 为什么公开 LUAD 上默认关闭

RNA velocity（La Manno 等, 2018；Bergen 等, 2020）从 spliced/unspliced 推断高维向量。2021–2026 文献持续记录**假设失败**、**流程脆弱**和**基准权衡**：

- **模型假设经常不成立。** 稳态 / 共用剪接率常被违反；即便动力学 scVelo 仍假定常数速率和形状良好的相位图（Bergen 等, 2021；Barile 等, 2021）。
- **定量选择会改变答案。** spliced/unspliced 流程影响 velocity（Soneson 等, 2021）。大多数**公开已处理**的 LUAD 对象只有一张基因计数矩阵，没有 `spliced`/`unspliced`。没有 FASTQ（或 velocyto/alevin-fry/STARsolo 重定量），velocity **算不了**。仅此一项，多数 GEO h5ad 的 velocity 门槛就失败。
- **k-NN 平滑支配整个流程。** Zheng、Stein-O’Brien、Boukas、Goff 与 Hansen（2023）表明方向和速度继承了观测到的 k-NN 图；图错了，高维和低维 velocity 都错。除极低噪声外，速度估计不可靠。把 velocity 映射到嵌入上“实际上是在嵌入空间里插值”。
- **用 velocity 验证 UMAP 是循环论证**（Zheng 等, 2023）。低维嵌入本身会扭曲（Chari & Pachter, 2023）。投影到 UMAP 会把正交过程（如周期 vs 分化）混在一起（VeloCycle, 2024）。
- **负对照也会画出箭头。** 成熟、稳态群体（如 PBMC）不应出现连贯流；许多方法会发明它（Zheng 等, 2023；Wu、Kong、Liao 等, *Genome Biology* 2026）。
- **没有方法赢得全部任务。** Wu 等（2026）在方向一致性、时间精度、**负对照稳健性**和测序深度稳定性上评了 25 个仅 RNA 方法。方向表现与负对照稳健性**负相关**（该研究中 Spearman ρ ≈ −0.57）。他们较均衡的推荐包括 UniTVelo (uni) / veloVI / Pyro-Velocity (m2)；LatentVelo (std) 可以方向极好而负对照很差。**没有单一方法始终好用。**
- **veloVI**（Gayoso、Weiler 等, 2024）是比点估计 scVelo 更好的*概率*动力学模型，不是跳过 QC 的许可证。
- **理论拆解**（Gorin、Fang、Chari 与 Pachter, 2022）仍是正确先验：连续转录–剪接–降解 ODE 严重简化了真实的爆发转录、细胞特异速率和基因依赖。

对 AT2→肿瘤这一**癌症**轴，假设不是更好而是更差：速率不是发育常数；恶性细胞往往接近新的稳态（恰是负对照区）；公开 LUAD 是快照，不是脉冲示踪。

### 6.2 Velocity 的 go / no-go（叠加在 G1–G6 之上）

仅当 **V1–V5 全部成立**才可计算 velocity。否则写明“未做 RNA velocity”及原因。

| ID | 要求 |
| --- | --- |
| V1 | 存在 spliced 与 unspliced（或用写明的工具对公开 FASTQ 重定量：STARsolo + velocyto、alevin-fry、kallisto\|bustools）。 |
| V2 | 检查**预指定**基因集（AT2、过渡、周期）的相位图。须有足够基因呈诱导/抑制环。空泡/团块相位图的基因删除，不平均进向量。 |
| V3 | 同一对象中的**负对照群体**（如成熟 T 细胞，或确实静态的静息 AT1）**不**产生自信的全局箭头。若一切都在流，方法在幻觉。 |
| V4 | 高维 velocity 用 CellRank `VelocityKernel` 分析（在表达图上，不在 UMAP 上）。UMAP 流线作附录或省略。 |
| V5 | **不解释速度**（Zheng 等, 2023）。方向只当作可与 Palantir 低权重组合的先验，不当发现引擎。 |

若 2024 年后仍要跑方法，优先 **veloVI** 或 Wu 等（2026）中方向**与**负对照都较均衡的方法。scVelo dynamical 只作敏感性，不当金标准。

### 6.3 禁止写的句子

按本手册做的论文里不要出现：

- “RNA velocity 证实了 UMAP 上从 AT2 到肿瘤的轨迹。”
- “细胞正迅速移向 TACSTD2 高的状态。”（速度 + 用基因定义的终点）
- “流线证明 KAC 变成 LUAD。”
- “我们在组织 UMAP 上用 velocyto/scVelo 默认参数。”

---

## 7. 把 `TACSTD2` 和 `CLDN4` 放到轴上

### 7.1 预先指定的待估量

1. **边际趋势。** 恶性谱系上，`CLDN4`（以及单独的 `TACSTD2`）是否为 Palantir 伪时间的函数？瞬时、晚期、还是平坦？
2. **命运关联。** 仅在**非恶性**上皮（AT2 + 过渡）内，各基因与恶性命运概率的 Spearman / CellRank driver 分数。把肿瘤终点放进相关会重复计算标签。
3. **途经点特异性。** `CLDN4` 是否与 `KRT8`/`CDKN1A` 一同达峰（KAC 样）而 `SFTPC` 下降？`TACSTD2` 是跟这个途经点走，还是跟气道对照走？
4. **轴外拒绝。** 限制到肺泡标注细胞后，`TACSTD2` 与 `KRT5`/`SCGB1A1` 是否仍相关。若仍相关，你不在 AT2→肿瘤轴上。

### 7.2 不用 UMAP 的统计

- **供体水平：** 每供体在 AT2 vs 过渡 vs 恶性的基因中位数；跨供体 Wilcoxon 或混合模型。这是主要*检验*。第 4–5 节趋势是主要*描述*。
- **供体内秩相关：** 单供体内基因 vs Palantir 伪时间，再做秩元分析。防止单个供体批次。
- **置换：** 在供体内打乱恶性命运标签，重算预指定基因的 driver 分数。
- **不要**检验所有基因再突出 `TACSTD2`。

### 7.3 蛋白与空间（若公开数据有）

CITE-seq 的 ADT 面板通常是免疫面板，**没有** TROP2 抗体。不要插补 TROP2 蛋白。若公开空间数据测了 `TACSTD2`/`CLDN4` RNA 或 TROP2/CLDN4 蛋白，用来检验**界面定位**（KAC 在 AT2 与肿瘤之间），不要画 velocity。

### 7.4 解释用语

允许：“在起点为 nLung AT2、终点为 CNV 高恶性的 Palantir/CellRank 轴上，`CLDN4` 表达与中间态 / 恶性命运概率相关，与已发表的 KAC / DATP 标记一致。”

不允许：“轨迹分析表明 TACSTD2 驱动 AT2 转化。”

---

## 8. 默认流程（公开肺，门槛已通过）

1. 下载公开对象；记录 accession、许可、基因 ID。
2. 按**样本** QC（环境 RNA、双细胞、线粒体）。保留 counts。
3. 粗注释；**限制到上皮**；重算 HVG/PCA/neighbors（第 3 节）。
4. 执行 G1–G6。写一段门槛报告。失败则停。
5. Palantir：预指定 AT2 起点与 CNV 恶性终点（第 4 节）。
6. CellRank 2 `PseudotimeKernel` ± `ConnectivityKernel`；GPCCA；命运概率（第 5 节）。
7. 可选：CytoTRACE 2 作为第二支箭头；比较后再考虑混合。
8. 对预指定基因（`TACSTD2`、`CLDN4`、对照）做趋势 + 供体水平检验（第 7 节）。
9. Velocity：除非 V1–V5 通过，否则跳过（第 6 节）。
10. 图：（i）按组织的状态组成；（ii）扩散 / Palantir 空间，而不是 velocity-UMAP；（iii）基因趋势对伪时间；（iv）按状态的命运概率；（v）预指定基因的 driver 表。若展示 UMAP，只当图例。

软件版本钉死在 `environment.yml`（见英文第 8.1 节）。敏感性：两种随机种子；`n_neighbors` ∈ {15, 30, 50}；Palantir `num_waypoints` ∈ {250, 500, 1000}；按特征间隙选择 CellRank `n_states`；留一供体：恶性终点是否仍落在 CNV 高细胞。**不要**为了让 `TACSTD2`“好看”而调这些旋钮。

