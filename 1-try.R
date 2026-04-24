# =============================================================================
# 1-try.R — end-to-end runner for the shipped ATP7B 7-gRNA example.
# Uses coverage-aware per-motif scoring (s/C) with the default edge guard.
#
# Output layout:
#   caspex_output/            <- primary deck (PDFs + CSVs)
#     extras/                 <- A1-A6, B1-B3, C1-C3, D1-D3 (auto-skips B4, C4)
# =============================================================================

source("caspex_analysis.R")

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
  edge_guard_frac  = 0.25                    # in-support mask threshold (5x cov_floor)
)

source("caspex_extras.R")
extras <- run_caspex_extras(result, out_dir = "caspex_output/extras")
