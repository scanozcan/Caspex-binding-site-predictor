# CASPEX TF-Binding Zone Predictor

A single-file R pipeline that turns **tiled-gRNA CasPEX / GLoPro / CAPLOCUS** proximity-labeling proteomics data into a spatial map of predicted transcription-factor (TF) binding zones on a gene's promoter.

## Background: locus-directed proximity proteomics

Proximity-labeling proteomics uses a promiscuous biotinylating enzyme — **APEX/APEX2** [^apex] or the engineered biotin ligase **TurboID** [^turboid] — to tag any protein sitting within a short (~10–300 nm) radius of the enzyme with biotin, so that downstream streptavidin enrichment and mass spectrometry report the local interactome of whatever the enzyme is fused to. Fusing the enzyme to **catalytically-dead Cas9 (dCas9)** and targeting it with a guide RNA turns the method into a *locus-directed* interactome reader: streptavidin-MS on cells expressing the fusion returns the set of proteins that were transiently near the targeted region of the genome during the biotinylation pulse.

The canonical embodiment is **CasPEX** — "dCas9-APEX2 at a defined locus" — introduced by Myers et al. [^caspex] in 2018, alongside the independently-reported **C-BERST** method [^cberst]. A family of variations has followed under names like **GLoPro**, **CAPLOCUS**, **CasID**, **enChIP-APEX**, and TurboID-based equivalents, all sharing the core idea of one gRNA → one localised biotin cloud → one MS readout of the local proteome. Reviews by Qin et al. [^qinreview] and Gingras et al. [^proxrev] survey the broader proximity-labeling field and the genomic-locus application.

The original CasPEX design uses a single gRNA per locus, which convolves all binders within the ~50–500 bp labeling radius into a single proteomic table. Tiling multiple gRNAs across a promoter samples the locus at several positions, and the resulting per-region differential-abundance tables carry *spatial* information about where each enriched protein was labeled — but that information is blurred by the labeling kernel and confounded by non-uniform gRNA tiling. This pipeline inverts the blur: it combines the per-region proteomics into a continuous labeling-weighted signal along the promoter, normalises that signal by the underlying labeling-opportunity map to remove the geometry of the tiling from the readout, and reports per-motif binding amplitudes that can be interpreted as coverage-corrected occupancy rates.

## What this pipeline does

Each gRNA tiles a ~20 bp window across the promoter. Proteomic logFC values from every region are re-projected onto the promoter coordinate system (relative to the TSS), combined into a Gaussian labeling-kernel signal, and decomposed into discrete binding events anchored to JASPAR (or HOCOMOCO) motif positions. Binding amplitudes are scored per motif against a **coverage-normalized signal**

```
β(x) = s(x) / max(C(x), cov_floor · max(C))
```

so that mid-gap binders — those sitting between two tiled gRNAs — are not down-weighted relative to binders sitting directly under a gRNA. An edge-guard mask (`edge_guard_frac`) prevents the coverage correction from inflating β at the flanks of the promoter window where labeling opportunity is thin.

- **Input**: a tab-separated manifest (`example_inputs/grnas.tsv`) pointing at one MS1-quantified table per gRNA-tiled promoter region, plus a gene symbol.
- **Output**: spatial model (CSV), binding-event table (CSV), 11 publication-quality PDFs.
- **Single file**: `caspex_analysis.R` sources cleanly, no package to install.
- **Alternative motif backend**: `caspex_analysis_hocomoco.R` runs the same pipeline with HOCOMOCO v11/v12 PWMs instead of JASPAR (offline after one-time download). See §11.

[^apex]: Rhee H-W, Zou P, Udeshi ND, Martell JD, Mootha VK, Carr SA, Ting AY. "Proteomic mapping of mitochondria in living cells via spatially restricted enzymatic tagging." *Science* 339, 1328–1331 (2013). doi:10.1126/science.1230593.

[^turboid]: Branon TC, Bosch JA, Sanchez AD, Udeshi ND, Svinkina T, Carr SA, Feldman JL, Perrimon N, Ting AY. "Efficient proximity labeling in living cells and organisms with TurboID." *Nat. Biotechnol.* 36, 880–887 (2018). doi:10.1038/nbt.4201.

[^caspex]: Myers SA, Wright J, Peckner R, Kalish BT, Zhang F, Carr SA. "Discovery of proteins associated with a predefined genomic locus via dCas9-APEX-mediated proximity labeling." *Nat. Methods* 15, 437–439 (2018). doi:10.1038/s41592-018-0007-1.

[^cberst]: Gao XD, Tu L-C, Mir A, Rodriguez T, Ding Y, Leszyk J, Dekker J, Shaffer SA, Zhu LJ, Wolfe SA, Sontheimer EJ. "C-BERST: defining subnuclear proteomic landscapes at genomic elements with dCas9-APEX2." *Nat. Methods* 15, 433–436 (2018). doi:10.1038/s41592-018-0006-2.

[^qinreview]: Qin W, Cho KF, Cavanagh PE, Ting AY. "Deciphering molecular interactions by proximity labeling." *Nat. Methods* 18, 133–143 (2021). doi:10.1038/s41592-020-01010-5.

[^proxrev]: Gingras A-C, Abe KT, Raught B. "Getting to know the neighborhood: using proximity-dependent biotinylation to characterize protein complexes and map organelles." *Curr. Opin. Chem. Biol.* 48, 44–54 (2019). doi:10.1016/j.cbpa.2018.10.017.

---

## 1. Quick start

The per-region gRNA sequences and logFC files are declared in a small manifest (`example_inputs/grnas.tsv`); the pipeline reads the manifest and wires everything up.

```r
source("caspex_analysis.R")

# Read example_inputs/grnas.tsv and resolve file paths under example_inputs/
inputs <- load_caspex_inputs("example_inputs")

result <- run_caspex(
  gene             = "gene_name",                # any HGNC symbol
  grnas            = inputs$grnas,
  data_files       = inputs$data_files,
  out_dir          = "caspex_output",
  coverage_correct = TRUE,                   # coverage-aware per-motif scoring
  cov_floor        = 0.05,
  edge_guard_frac  = 0.25
)

# Optional supplementary analyses (one-pagers, permutation null,
# σ sensitivity, jackknife, volcanoes, QC, coverage diagnostics).
source("caspex_extras.R")
extras <- run_caspex_extras(result, out_dir = "caspex_output/extras")
```

The thin wrappers `1-try.R` and `1-try-hocomoco.R` reproduce this end-to-end for the shipped example; run `source("1-try.R")` for the JASPAR-backed run or `source("1-try-hocomoco.R")` for the HOCOMOCO-backed run.

---

## 2. Input format

### 2.1 The `example_inputs/` folder

```
example_inputs/
├── grnas.tsv          # manifest: region × gRNA × data_file
├── Region1.txt        # per-region differential-abundance table
├── Region2.txt
├── ...
└── Region7.txt
```

### 2.2 `example_inputs/grnas.tsv` — per-region manifest

One row per tiled gRNA region. Lines starting with `#` are comments; blank lines are ignored. Columns (tab-separated):

| column | meaning |
|---|---|
| `region` | Region label used throughout the pipeline (`R1`, `R2`, ...). Must be unique. Case-sensitive. |
| `sequence` | Protospacer (17–23 bp). PAM (`NGG` / `CCN`) is optional; the loader strips it. Use `NA` (uppercase) if a region has no matched gRNA sequence — the pipeline will still use its logFC data but plot no cut-site tick for it. |
| `data_file` | Filename of the logFC / differential-expression table for that region, **relative to the `example_inputs/` directory**. |

Example (shipped with the repository):

```
region	sequence	data_file
R1	GGGCTCACCCACCGCCTGCG	Region1.txt
R2	CACTGCTGTCCCGCCCACGG	Region2.txt
R3	AAGTTGTTGAGGTGAGCAGC	Region3.txt
R4	TTTAAGTGACGTGTTAACAA	Region4.txt
R5	TACTGCTATTCAGGAAGCTA	Region5.txt
R6	TCATAAAAGCAACACTACAG	Region6.txt
R7	GGCTAACTTAGACATGTAGT	Region7.txt
```

### 2.3 Per-region differential-abundance tables

Each `Region*.txt` must be a tab-separated file with at least these columns:

| column | meaning |
|---|---|
| `Protein` or `protein` | Gene/protein symbol (HGNC) |
| `logFC` | log2 enrichment (labeling region vs. control) |
| `pval` or `p.value` | test p-value |
| `t` *(optional)* | moderated t-statistic from `limma::topTable`. If present, used as the default region weight (statistically preferred). If absent, a signed z-score is derived from the p-value. |
| `TFDatabase` *(optional)* | `"exist"` if the protein is in a curated TF list |

The loader (`load_region_data()`) normalises column casing and attaches a `region` column. The manifest loader (`load_caspex_inputs()`) resolves data-file paths relative to the `example_inputs/` directory, validates that the manifest has the three required columns, and returns both the gRNA vector and the data-file vector wired by region label.

---

## 3. Pipeline sections (inside `caspex_analysis.R`)

The file is structured into labelled sections; jump to each by searching for `SECTION N:`.

### Section 1 — IO & normalisation
- `load_region_data(files)` → long data.frame with `region, protein, lfc, pval, isTF`
- `load_caspex_inputs(inputs_dir, manifest)` — reads `example_inputs/grnas.tsv`, validates columns, resolves relative paths, returns `list(grnas, data_files, inputs_dir, manifest_path)`.
- Case-insensitive column matching; silently drops NA rows; prints per-region counts.

### Section 2 — Gene & promoter fetch (Ensembl REST)
- `lookup_gene(symbol, species = "human")` — HGNC → chromosome, TSS, strand; `species` is attached to the returned object and carries through to sequence fetch.
- `fetch_promoter_seq(gene_info, upstream, downstream, species = NULL)` — downloads the reference sequence via Ensembl Sequence API; returns a character vector plus coordinate bookkeeping. When `species` is not supplied it is read off `gene_info$species`, so the whole species-aware chain is consistent.

### Section 3 — gRNA mapping
- `match_all_grnas(grnas, promoter_info)` — Aho-Corasick-style substring scan in the promoter (both strands), returns a named numeric vector of TSS-relative positions (`pos_map`).
- Silently tolerates unmatched gRNAs (returned as `NA`).

### Section 4 — Spatial model
- `run_spatial_model(long_data, pos_map, ..., min_lfc = 0)` computes per-TF:
  - **Centroid**: `Σ(w⁺ · pos) / Σw⁺` with weights `w` chosen by `weight_mode`.
  - **Spread**: weighted SD around the centroid.
  - **Composite**: `mean(w⁺) × log1p(n_regions)`.
- Only `logFC > min_lfc` contributes to the spatial stages (`min_lfc = 0` by default; depletion is absence of evidence, not anti-binding).
- Output: `spatial_df` sorted by composite descending.

**Region weights (`weight_mode`).** Default is `"mod_t"`: use the moderated t-statistic from `limma::topTable` if the input files carry a `t` column; otherwise fall back to a signed z-score derived from the two-sided p-value:

```
z_i = sign(logFC_i) · Φ⁻¹(1 − p_i/2)
```

A z-score is the statistically principled "signal per standard error" — exactly the quantity an inverse-variance-weighted spatial mean should use. It also compresses extreme p-values (p = 1 × 10⁻¹⁰ → z ≈ 6.4, not 10), preventing one spectacular outlier from dominating a centroid estimate made from a handful of regions. The legacy `w = logFC × −log10(p)` is kept as `"lfc_x_negp"` for backward compatibility.

Other modes: `"z"` (always derive z from p, ignore t column), `"lfc_pos"` (effect size only), `"lfc_signed"` (raw logFC).

### Section 5 — JASPAR motif scanning (REST API v1)
- `fetch_jaspar_pwm(tf_name)` — queries `jaspar.elixir.no/api/v1/matrix/?tf_name=...&collection=CORE&tax_group=vertebrates`, paginates, applies three-tier fallback (`tf_name=` → `search=` → `name=`), and returns only PWMs whose matrix name **actually matches** the query (exact, case-insensitive match via `pick_match()`). Returns `NULL` if no real name match — never a first-hit fallback.
- `scan_pwm(pwm, seq, thresh)` — sliding-window log-odds score with reverse-complement scanning; returns all TSS-relative hit positions ≥ `thresh × max_score`.
- `run_motif_scan(tfs, promoter_info, thresh)` — wraps the two above for a TF set. Default threshold = 0.80.

### Section 5.5 — Binding-event decomposition (coverage-aware per-motif scoring)

Problem: proteomic CasPEX labeling is diffuse (~50–500 bp radius), so two binding events in adjacent regions can look like one, and one binding event can light up multiple regions. We disentangle them by scoring candidate motifs against a **coverage-normalized signal** that inverts the labeling-decay geometry.

**Signal and coverage maps.** On a 5-bp grid across the promoter window we build both the proteomic signal and an unweighted *labeling-opportunity* map:

```
s(x) = Σ_r w_r⁺ · G(x − pos_r; σ)       (positive-part region-weighted signal)
C(x) = Σ_r        G(x − pos_r; σ)       (unweighted — every gRNA contributes)
```

with `σ = 300 bp` (default `kernel_sigma`), the APEX labeling radius.

**Candidate positions.** Either JASPAR / HOCOMOCO motif hits for that TF (preferred) or local maxima of the raw signal `s(x)` inside the in-support mask (fallback for TFs without a PWM).

**Per-motif score.** For each candidate motif `m_k`:

```
β_k = s(m_k) / max( C(m_k), cov_floor · max(C) )
```

and keep `β_k > min_weight_frac · max_k β_k` (default `min_weight_frac = 0.15`). Close survivors are collapsed with the `merge_dist = 100 bp` rule by amplitude-weighted centroid.

*What β represents.* A per-unit-labeling occupancy rate: the effect size normalized by how much labeling opportunity was deposited at that position. In a 2-gRNA overlap both s and C roughly double so β returns the individual-region weight; in a gRNA-poor tail both approach zero and the `cov_floor` prevents unbounded amplification of tail noise.

**Edge guard.** The `cov_floor` clamp alone is not enough to stop the west/east edges from sprouting spurious calls. Just inside the clamp boundary — where `C(x)` is *barely above* `cov_floor · max(C)` — any residual kernel tail of `s(x)` gets divided by a tiny denominator, and the resulting inflated `β(x)` can form a contiguous above-threshold zone that touches the grid edge. The fix is a second, stricter threshold `edge_guard_frac` (default `0.25` = 5× default `cov_floor`) that defines the **in-support region** for β:

```
support_mask = C(x) > max(cov_floor, edge_guard_frac) · max(C)
β(x)[!support_mask] = 0     # zone / peak detection only inside the mask
```

Both the motif-based zone detection and the motif-less peak-finding branch respect this mask; the motif-less branch additionally peak-finds on the raw `s(x)` (not `β(x)`) so edge ramping cannot create phantom local maxima there. Raise `edge_guard_frac` toward `0.30`–`0.35` on sparsely-tiled genes if west-edge calls still appear; lower it toward `cov_floor` on densely-tiled regions where you trust β further into the transition band.

**Diagnostic columns.** `binding_events.csv` carries `distance_to_nearest_grna` and `local_coverage = C(event_pos)`. A high `β_k` with low `local_coverage` is a **coverage-rescued** call — a motif that a naive smoothed-signal fit would have clipped because the raw signal at its position is suppressed by labeling geometry. Inspect those first when validating results.

Public entry points:
- `build_caspex_signal(tf, long_data, pos_map, x_grid, kernel_sigma, weight_mode)`
- `predict_binding_events(tf, ..., motif_hits, coverage_correct = TRUE, cov_floor, edge_guard_frac, ...)` — per TF
- `predict_all_binding_events(tfs, ..., coverage_correct = TRUE, cov_floor, edge_guard_frac, ...)` — batch wrapper; returns a tidy data.frame with `tf, position, weight, motif_based, n_motifs_merged, distance_to_nearest_grna, local_coverage`.

### Section 6 — Plots
See the output deck below.

### Section 7 — `run_caspex()` main wrapper
The top-level function. Orchestrates load → spatial model → motif scan → decomposition → plots → CSV export. Returns `invisible(list(...))` so interactive users get a handle on every intermediate object.

### Section 8 — Convenience helpers
- `inspect_tf(result, tf)` — pretty console summary + per-region profile plot.
- `select_motif_tfs(...)` — TF-selection heuristic (see below).

---

## 4. TF selection heuristic (`select_motif_tfs`)

We scan JASPAR (or HOCOMOCO) for a curated subset rather than every modelled TF (API calls take ~1–2 s each under the JASPAR backend). Every motif-scanned TF is placed into **exactly one** of three disjoint buckets:

**(1) Common TFs** (default 10, `n_common`): top-N by composite score across all regions — broadly strong hits. Rendered in decks 06/07.

For the remaining buckets, we first build a long candidate table of every (TF, region) pair that passes `pval ≤ pval_thresh` **and** `logFC > 0` **and** is not already in common, scoring each row by

```
specificity(p, R) = lfc_R(p) − mean( lfc_{R'}(p) : R' ≠ R )
```

**(2) Shared-focal TFs** (default 10, `n_shared`): TFs that are significant in **≥ 2 regions**, ranked by their **summed specificity across the regions in which they are significant**. Each shared TF carries a tag like `R1+R5` listing *all* of its significant regions. Rendered in decks 11/12 with lane labels `POU3F2  [R1+R5]`.

**(3) Region-specific TFs** (default 10 per region, `n_specific`): TFs significant in **exactly one region**. Per-region top-N by specificity. Rendered in decks 08/09, one page per region.

Why three buckets and not two? A TF enriched strongly in R1 *and* R5 (but diffuse elsewhere) shouldn't be forced into a single region's "specific" lane (it isn't specific to R1 alone), nor should it inflate the "common" deck (it isn't broad). The shared-focal bucket gives it a home and preserves the *set* of regions where it's enriched.

**Disjointness guarantees** (enforced by post-condition assertions):
- `common ∩ shared        == ∅` — shared is built from non-common candidates only.
- `common ∩ region-specific == ∅`
- `shared ∩ region-specific == ∅` — region-specific is built from single-region candidates only.
- No TF appears on two pages of the region-specific deck.

The returned character vector carries attributes: `common`, `shared` (character vector), `shared_tag` (named: TF → `"R1+R5"`), `specific` (named list by region), `region_tag`, and `source` so downstream plotting can label lanes with provenance.

---

## 5. Output deck (`caspex_output/`)

| file | what it shows |
|---|---|
| `00_summary.pdf` | gRNA ruler + GC track + spatial ridge plot, stacked |
| `01_grna_positions.pdf` | gRNA cut sites on a promoter ruler |
| `02_spatial_track.pdf` | Top-25 TFs as Gaussian ridges on the promoter (centroid ± spread) |
| `03_heatmap.pdf` | logFC × region matrix, top-30 TFs |
| `04_gc_content.pdf` | 50-bp GC% track |
| `05_motif_overlay.pdf` | Legacy per-TF motif density track (kept for reference) |
| `06_track_common.pdf` | Common TFs (top-N composite) as mini genome-browser lanes, **no** motif overlay |
| `07_track_common_motif.pdf` | Same lanes + JASPAR motif ticks in a sub-lane **below** each ribbon |
| `08_track_region_specific.pdf` | **Multi-page**: one page per region, showing that region's specific TFs |
| `09_track_region_specific_motif.pdf` | Same, with JASPAR ticks |
| `10_binding_deconvolution.pdf` | **Multi-page**: one page per TF, showing signal + per-region logFCs + called events in a strip **below** the baseline, with `C(x)` overlaid as a dashed grey curve |
| `11_track_shared.pdf` | **Shared-focal deck**: TFs significant (p ≤ thr & logFC > 0) in ≥ 2 regions, not in common. Lane labels carry the full region-set tag (e.g. `POU3F2  [R1+R5]`). Emitted only when ≥ 1 shared TF exists. |
| `12_track_shared_motif.pdf` | Same lanes + JASPAR motif ticks. |
| `<GENE>_spatial_predictions.csv` | `spatial_df` |
| `<GENE>_binding_events.csv` | Tidy decomposition calls; columns `tf, position, weight, motif_based, n_motifs_merged, distance_to_nearest_grna, local_coverage, deck` — `deck` is one of `common`, `shared:R*+R*`, `region-specific:R*`, mirroring the 06/07, 11/12, and 08/09 PDF splits. |

### Plot 06 / 07 — lane layout

Each TF occupies a horizontal lane. Within the lane:

```
│ ┌─────────────────────────────────────────┐ │
│ │           Gaussian ribbon               │ │ ← enrichment
│ │          (centroid + spread)            │ │
│ ├─────────────────────────────────────────┤ │
│ │ │││       │     │   │  (motif ticks)    │ │ ← JASPAR hits
│ └─────────────────────────────────────────┘ │
```

Motif ticks live in a **sub-lane below** the ribbon — not on top, not spanning the whole height. The vertical centroid tick spans only the ribbon, not the motif sub-lane.

### Plot 10 — binding decomposition layout

```
│       ┌───────┐                       │
│      / ∩  ∩   \   signal (logFC-      │ ← upper panel
│     /  •  •    \  weighted Gaussian)  │
│- - -- - -- - -- - - C(x) dashed grey -│ ← coverage overlay
│────|─ R2 R4 ────|──── baseline ───────│ ← y = 0
│    |                                  │
│    ●  R3   ●  R6   negative-lfc        │ ← below 0 but above track
│                    lollipops           │
│    :         :                         │ ← shaded strip (floats below
│    ||  ||    ||   motif ticks          │    the lowest lollipop)
│    ●           ●  predicted events     │ ← called-peak track
│    +120       −800    (labels)         │
```

The called-peak strip is anchored to `min(0, lfc) − 0.06·sig_max`, so depletion lollipops (negative lfc) are never buried under the motif/event lane. When all lfc ≥ 0 this collapses back to "strip sits directly below zero". Dotted connectors run from each predicted event up to the signal value at that x, so a reviewer sees which signal peak each call explains. The dashed grey `C(x)` overlay (scaled to the signal max) shows at a glance where labeling opportunity was strong versus where a call depended on the coverage correction.

**What does a motif-anchored bubble actually show?** Each filled red/orange bubble in the called-peak track is a **surviving motif position** after the per-motif scoring step — not a PWM score readout:

- **Position (x)** = TSS-relative coordinate of a JASPAR (or HOCOMOCO) PWM hit that passed the log-odds cutoff (`threshold_frac = 0.80`) *and* received a non-trivial β. If two or more hits within `merge_dist = 100 bp` both survived, they collapse into one bubble at their amplitude-weighted centroid (see `n_motifs_merged` in the events table).
- **Size** = the coverage-normalized amplitude `β_k = s(m_k) / max(C(m_k), cov_floor · max(C))` (or `Σ β_k` for merged clusters). A strong PWM hit with no CasPEX enrichment around it gets β = 0 and disappears; a weak PWM hit under a strong labelled peak gets a large β and shows up as a big bubble. **It is NOT the PWM match score.**
- **Fill colour**: red/orange (`motif_based = TRUE`) = anchored to a literal PWM hit on the promoter; teal (`motif_based = FALSE`) only appears when the TF has no PWM, in which case the bubble is a local maximum of `s(x)` inside the in-support mask and carries no sequence anchoring.

Only motifs with `β_k > 0.15 · max(β)` become bubbles. All the other PWM hits still appear as tick marks in the sub-lane but don't generate an event call.

---

## 6. Parameter cheat-sheet

| parameter | default | meaning |
|---|---|---|
| `upstream`, `downstream` | 2500, 500 | bp around TSS to fetch/plot |
| `pval_thresh` | 0.05 | significance cutoff |
| `min_regions` | 2 | drop TFs seen in fewer regions |
| `min_lfc` | 0 | spatial-model floor on `logFC`: only `logFC > min_lfc` contributes to centroid/spread/composite. Leave at `0` for the proximity-labeling use case. Plumbed through `run_caspex()` → `run_spatial_model()` → `compute_spatial()`. |
| `top_n` | 25 | lanes in plot 02 |
| `n_common` | 10 | top-composite TFs for motif scan (deck 06/07) |
| `n_shared` | 10 | shared-focal TFs significant in ≥ 2 regions (deck 11/12) |
| `n_specific` | 10 | region-specific TFs per region (deck 08/09) |
| `motif_thresh` | 0.80 | PWM score as fraction of max (JASPAR/HOCOMOCO) |
| `kernel_sigma` | 300 | APEX labeling radius (bp) |
| `min_weight_frac` | 0.15 | minimum normalised weight to keep a peak |
| `min_peak_dist` | 150 | bp between local maxima in the motif-less fallback |
| `merge_dist` | 100 | merge adjacent kept peaks within this distance |
| `weight_mode` | `"mod_t"` | region-weight mode (`mod_t` / `z` / `lfc_pos` / `lfc_signed` / `lfc_x_negp`) — applies to both spatial model and signal |
| `signal_weight` | `NULL` | override `weight_mode` for the signal-building step only (backward-compat alias) |
| `coverage_correct` | `TRUE` | enable coverage-aware per-motif scoring. Leave at `TRUE` for the documented behaviour; the code retains an internal `FALSE` branch for comparison experiments but it is not packaged here. |
| `cov_floor` | 0.05 | relative floor on the denominator; caps inverse-coverage amplification at ~1/cov_floor |
| `edge_guard_frac` | 0.25 | fraction-of-max-coverage floor for the in-support mask. β is trusted only where `c_grid(x) > edge_guard_frac · max(c_grid)`, so zones and peaks cannot form in the low-denominator transition band just inside the `cov_floor` clamp. Must be ≥ `cov_floor`; raise if you still see west-edge calls on sparsely-tiled genes. |
| `tfs_only` | `TRUE` | restrict spatial model to `TFDatabase == "exist"` |

---


## 7. Dependencies

Base R (≥ 4.0) plus:

- `httr`, `jsonlite` — Ensembl + JASPAR REST clients
- `ggplot2`, `dplyr`, `tidyr`, `scales`, `patchwork` — plotting
- `nnls` (optional) — kept in the code path as an internal fallback for a legacy smoothed-signal decomposition; not required for the coverage-aware path.

Install everything in one go:

```r
install.packages(c("httr","jsonlite","ggplot2","dplyr","tidyr",
                   "scales","patchwork","nnls"))
```

---

## 8. Return value

`run_caspex()` invisibly returns a list with:

```
$spatial_df       # per-TF centroid / spread / composite
$long_data        # normalised per-region logFC table
$pos_map          # named numeric vector: region → TSS-relative bp
$gene_info        # Ensembl metadata
$promoter_info    # sequence + coordinate bookkeeping
$motif_results    # list per TF: hits (bp), score threshold, matrix ID
$binding_events   # tidy data.frame of decomposition calls
$weight_mode      # the weight mode requested ("mod_t" by default)
$weight_source    # what was *actually* used — "moderated t from input" or
                  #   "signed z-score (derived from p-value; no `t` column in input)"
$coverage_correct # TRUE for the packaged pipeline
$cov_floor        # 0.05 by default
$edge_guard_frac  # 0.25 by default
$plots            # named list of ggplot objects
  $grna, $track, $heat, $gc, $motif,
  $common, $common_motif,
  $region_specific        # named list: region → ggplot
  $region_specific_motif  # named list: region → ggplot
  $shared                 # ggplot or NULL if no shared-focal TFs
  $shared_motif           # ggplot or NULL
```

At the end of every run, a structured summary block is also printed to the console so you can audit the analytical path without digging into the object:

```
--- Run summary --------------------------------------------------
 Weight mode      : mod_t
 Weight source    : signed z-score (derived from p-value; no `t` column in input)
 Kernel sigma (bp): 300
 Binding path     : coverage-normalized per-motif β=s/C (cov_floor=0.05, edge_guard_frac=0.25)
 Spatial TFs      : 139
 Binding events   : 88 (40 TFs)
------------------------------------------------------------------
```

Interactive usage:

```r
result$plots$region_specific_motif[["R3"]]   # single-region plot
inspect_tf(result, "GATA6")                   # per-TF inspector
```

---

## 9. Files in this folder

- `caspex_analysis.R` — the pipeline (single source of truth, JASPAR-backed).
- `caspex_analysis_hocomoco.R` — **alternative motif backend**: same pipeline wired to HOCOMOCO v11/v12 PWMs. Sources `caspex_analysis.R` automatically. See §11.
- `caspex_chipatlas.R` — optional public-ChIP-seq validation overlay (ChIP-Atlas). Auto-sourced by `caspex_analysis.R` when present. See §11.5.
- `caspex_extras.R` — supplementary result plots and diagnostics (per-TF one-pagers, TF co-occurrence, ranked events, permutation null, σ sensitivity, jackknife, volcanoes, QC, plus coverage-aware diagnostics D.1–D.3). Source *after* `caspex_analysis.R`; see §12.
- `1-try.R` — thin runner for the shipped gene specific 7-gRNA dataset (JASPAR backend, coverage-aware).
- `1-try-hocomoco.R` — same runner using the HOCOMOCO backend.
- `2-run_caspex_extras_in_all.R` — regenerate the supplementary plots for an existing `result` object without re-running the primary pipeline.
- `example_inputs/` — per-region manifest (`grnas.tsv`) and the seven `Region*.txt` tables for the shipped example.
- `caspex_output/` — everything the pipeline writes (PDFs + CSVs above), including an `extras/` subfolder when the extras runner has been invoked.
- `caspex_output_hocomoco/` — default output folder of the HOCOMOCO runner.
- `README.md` — this document (user/developer guide).

---

## 10. HOCOMOCO backend (`caspex_analysis_hocomoco.R`)

An alternative motif source for the *exact same* pipeline. Same spatial model, same binding-event decomposition, same plots — the only thing that changes is where TF position-weight matrices come from.

### Why a HOCOMOCO backend?

- **Offline after one download.** The JASPAR backend issues one REST call per TF (~1–2 s each). The HOCOMOCO backend downloads a single MEME bundle (a few MB) into a persistent cache on first use, then works entirely offline.
- **Different coverage.** HOCOMOCO (<https://hocomoco11.autosome.org>, v12 at <https://hocomoco12.autosome.org>) curates human and mouse TFs with explicit quality grades (A > B > C > D); JASPAR CORE is broader taxonomically but collapses multiple PWMs per TF with its own ranking. Running both is a useful sanity check — TFs whose called events are consistent across databases are the most trustworthy.
- **One best PWM per gene.** HOCOMOCO ships many PWMs per TF (different experimental sources, trimming variants). The backend keeps exactly the highest-rated PWM per gene symbol, so the pipeline sees a clean `TF → PWM` map without needing to deduplicate.

### Quick start

A ready-to-run script (`1-try-hocomoco.R`) mirrors `1-try.R` but calls the HOCOMOCO backend. Source it to reproduce the analysis end-to-end:

```r
source("1-try-hocomoco.R")
```

Or call the wrapper directly:

```r
source("caspex_analysis_hocomoco.R")   # also sources caspex_analysis.R

inputs <- load_caspex_inputs("example_inputs")

result <- run_caspex_hocomoco(
  gene             = "gene_name",
  grnas            = inputs$grnas,
  data_files       = inputs$data_files,
  out_dir          = "caspex_output_hocomoco",
  hocomoco_version = "v12",          # "v12" (default) or "v11"
  hocomoco_species = "human",        # "human" (default) or "mouse"
  coverage_correct = TRUE,
  cov_floor        = 0.05,
  edge_guard_frac  = 0.25
)
```

### What actually happens under the hood

1. `caspex_analysis.R` is sourced first — every non-motif function (loader, Ensembl lookup, spatial model, plots) is identical to the JASPAR build.
2. `download_hocomoco_bundle(version, species)` ensures the MEME file for the requested version/species is present in `tools::R_user_dir("caspex", "cache")`. First call hits the network; subsequent calls are free.
3. `parse_meme_file()` reads the bundle into a list of PPMs. For each matrix, `.extract_gene_symbol()` handles both v11 (`SOX2_HUMAN.H11MO.0.A`) and v12 (`SOX2.H12CORE.0.P.B`) naming, and `.extract_quality()` extracts the A/B/C/D grade. Only the best-rated PWM per gene symbol survives.
4. `run_caspex_hocomoco()` swaps the `run_motif_scan` function with a HOCOMOCO-backed closure for the duration of the call, invokes the unchanged `run_caspex()`, then restores the original on exit.
5. `.rewrite_jaspar_in_plots()` rewrites "JASPAR" → "HOCOMOCO v12" in plot titles/subtitles so the output deck labels its source correctly.
6. The returned list gains three extra fields: `motif_source = "HOCOMOCO"`, `hocomoco_version`, `hocomoco_species`.

### Helpers

```r
list_hocomoco_tfs(version = "v12", species = "human")
check_hocomoco_coverage(c("GATA6","NFKB1","XYZ1"))
clear_hocomoco_memory_cache()
download_hocomoco_bundle("v12", "human", force = TRUE)
```

### If the automatic download is blocked

HOCOMOCO's canonical host is not always reachable from firewalled networks. The downloader tries each known URL in order; if all fail, download the MEME bundle manually from <https://hocomoco12.autosome.org/downloads_v12> (e.g. `H12CORE_meme_format.meme` for human v12) and pass it via `hocomoco_from_file = "/path/to/H12CORE_meme_format.meme"`. Equivalently, prime the cache once and then drop the argument on subsequent calls:

```r
download_hocomoco_bundle("v12", "human",
                         from_file = "/path/to/H12CORE_meme_format.meme",
                         force = TRUE)
```

Cache location is `tools::R_user_dir("caspex", which = "cache")`.

---

## 11.5. ChIP-Atlas validation overlay (optional)

Turn on public ChIP-seq peaks from [ChIP-Atlas](https://chip-atlas.dbcls.jp) as a validation track by passing `chipatlas = TRUE` to `run_caspex()`:

```r
result <- run_caspex(
  gene       = "ATP7B",
  grnas      = inputs$grnas,
  data_files = inputs$data_files,
  chipatlas                 = TRUE,   # fetch peaks for every motif-scanned TF
  chipatlas_threshold       = "05",   # "05" = Q<1e-5 (default), "10", "20"
  chipatlas_max_experiments = 50      # newest-first SRX cap per TF
)
```

Peaks render as a sub-lane under each TF's binding-event bubbles in `10_binding_deconvolution.pdf` (one thin row per SRX experiment) and as union-peak strips in the mini-browser decks (06/07, 08/09, 11/12). You can use this to sanity-check predicted binding events against the published literature without leaving the plot.

**First run downloads** (cached to `tools::R_user_dir("caspex", "cache")/chipatlas/`):

- `experimentList.tab` — ChIP-Atlas's master metadata file (~345 MB). One-time.
- Per-SRX narrowPeak BEDs for every motif-scanned TF — typically a few MB each, capped at `chipatlas_max_experiments` per TF, newest SRX first.

Subsequent runs read from cache and are fully offline. Set `chipatlas = FALSE` (default) to skip the feature entirely and keep the pipeline network-free after the initial JASPAR/HOCOMOCO motif scan.

See `caspex_chipatlas.R` for the backend (cache management via `clear_chipatlas_cache()`, per-TF fetch via `fetch_chipatlas_peaks()`, batch wrapper `run_chipatlas_scan()`).

---

## 12. Supplementary result plots (`caspex_extras.R`)

Source this file *after* `caspex_analysis.R` and a successful `run_caspex()` call. It adds a set of interpretation, robustness, QC, and coverage-diagnostic plots without touching the core pipeline.

```r
source("caspex_analysis.R")
source("caspex_extras.R")

inputs <- load_caspex_inputs("example_inputs")
result <- run_caspex(gene = "gene_name",
                     grnas = inputs$grnas,
                     data_files = inputs$data_files,
                     coverage_correct = TRUE,
                     cov_floor        = 0.05,
                     edge_guard_frac  = 0.25)

extras <- run_caspex_extras(result, out_dir = "caspex_output/extras",
                             n_perm = 500)
```

### A — biological interpretation

| file | what it shows |
|---|---|
| `A1_tf_one_pagers.pdf` | One page per TF: decomposition view on top, per-region logFC bars with p-value stars in the middle, a one-line stats block at the bottom. Default set is top-10 common + all region-specific motif-scanned TFs. |
| `A2_tf_family.pdf` | JASPAR TF-family representation among the top-30 composite hits. |
| `A3_event_density.pdf` | Weighted histogram of all predicted binding events across all TFs on the TSS-relative axis. |
| `A4_composite_vs_specificity.pdf` | Scatter of composite enrichment vs per-TF specificity score. |
| `A5_tf_cooccurrence.pdf` | Symmetric heatmap: for each pair of TFs with called events, the fraction of one TF's events that sit within ±100 bp of an event of the other TF. |
| `A6_ranked_events.pdf` + `.csv` | Top-50 individual binding events ranked by a composite confidence score: per-TF-normalized β plus jackknife survival fraction (from B.3; defaults to 1 if B.3 is skipped). |

### B — statistical robustness

| file | what it shows |
|---|---|
| `B1_permutation_null.pdf` + `.csv` | Observed composite score per TF vs its permutation-null distribution from `n_perm` shuffles of the region-position map. Per-TF empirical p-value and BH-adjusted FDR. |
| `B2_sigma_sensitivity.pdf` | Decomposed event positions at σ ∈ {100, 200, 250, 300, 500} bp for the top-15 motif-scanned TFs. |
| `B3_event_jackknife.pdf` + `.csv` | Leave-one-region-out jackknife: for each called event, how often it survives when a region is dropped (tolerance ±100 bp). |

### C — QC / data sanity

| file | what it shows |
|---|---|
| `C1_volcano_per_region.pdf` | Classic volcano (logFC vs −log10 p) per region; TFs highlighted. |
| `C2_region_correlation.pdf` | Pearson correlation heatmap of logFC vectors across regions. |
| `C3_pval_histograms.pdf` | Per-region p-value histogram — should be roughly uniform with an excess near zero. |

### D — coverage-aware diagnostics

| file | what it shows |
|---|---|
| `D1_coverage_rescue.pdf` | Scatter of event β vs local coverage `C(event_pos)`, coloured by distance-to-nearest-gRNA, triangle-marked where the event sits near the `cov_floor`. The most direct audit of the coverage correction: any "rescued" distal call shows up here as a triangle in the top-left. |
| `D2_covfloor_sensitivity.pdf` + `.csv` | Analog of B.2 but sweeping `cov_floor ∈ {0.02, 0.05, 0.10, 0.20}`. For each motif-scanned TF, tracks the positions of the top-3 events across the sweep. Flat lines = robust; drift to the left = floor-sensitive rescue. Sweep floors are configurable via `cov_floors`. |
| `D3_coverage_stack.pdf` | Three-panel per-TF diagnostic: `s(x)` (signal), `C(x)` (coverage), and `β(x) = s/max(C, floor·max C)` stacked with called events overlaid on the β panel. |

### Skip-list

You can override the steps that are run with an explicit `skip` argument:

```r
run_caspex_extras(result, skip = c("permutation","jackknife","sigma"))
```

Individual functions return ggplot objects so they can be edited interactively (e.g. `p + labs(title = "...")`), and the heavy compute functions are split into `run_*` (returns a data object) and `plot_*` (renders it) so results can be cached across plotting tweaks.

---

## 13. Citing CasPEX

If you use this pipeline in a publication, please cite the canonical CasPEX method paper (Myers et al. 2018) alongside this repository. A formal citation for the pipeline itself will be added here once the accompanying paper is published. The background section above lists the foundational proximity-labeling and dCas9-proximity-labeling references that should typically be cited together with any result produced by this tool.

---

## License

No license file is included at this time — all rights reserved by the author. A formal license will be added alongside the manuscript release. Contact the author for permitted-use inquiries until then.
