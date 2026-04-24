# =============================================================================
# 1-try.R — end-to-end runner for the shipped ATP7B 7-gRNA example.
# Uses coverage-aware per-motif scoring (s/C) with the default edge guard and
# optional ChIP-Atlas public-ChIP-seq validation tracks (see §13 in README.md).
#
# Output layout:
#   caspex_output/            <- primary deck (PDFs + CSVs)
#     extras/                 <- A1-A6, B1-B3, C1-C3, D1-D3 (auto-skips B4, C4)
# =============================================================================

source("caspex_analysis.R")   # auto-sources caspex_chipatlas.R if present

# Per-region gRNA sequences + logFC files live in example_inputs/grnas.tsv
# (edit that file to add/remove regions or change protospacers).
inputs <- load_caspex_inputs("example_inputs")

result <- run_caspex(
  gene             = "ATP7B",                # any HGNC symbol
  grnas            = inputs$grnas,
  data_files       = inputs$data_files,
  out_dir          = "caspex_output",
  coverage_correct = TRUE,                   # coverage-aware per-motif scoring
  cov_floor        = 0.05,                   # relative floor on the denominator
  edge_guard_frac  = 0.25,                   # in-support mask threshold (5x cov_floor)
  # Edge-bleed cap: drop events where either boundary gRNA (westernmost or
  # easternmost guide) contributes >25% of the local Gaussian weight sum.
  # Suppresses tail-inflated events past the outermost guides. NULL disables.
  edge_grna_weight_cap      = 0.25,
  # --- Optional ChIP-Atlas validation overlay (public ChIP-seq peaks) -----
  # Fetches peak BEDs for every motif-scanned TF and renders them as a
  # sub-lane under the Plot-10 binding-event bubbles plus union-peak strips
  # in the mini-browser decks. FIRST RUN downloads experimentList.tab
  # (~345 MB) + per-SRX BEDs (capped at chipatlas_max_experiments per TF);
  # everything is cached under tools::R_user_dir("caspex","cache")/chipatlas/.
  # Set chipatlas=FALSE to keep the run fully offline.
  chipatlas                 = TRUE,
  chipatlas_threshold       = "05",          # "05" (Q<1e-5, default), "10", "20"
  chipatlas_max_experiments = 50
)

source("caspex_extras.R")
extras <- run_caspex_extras(result, out_dir = "caspex_output/extras")
