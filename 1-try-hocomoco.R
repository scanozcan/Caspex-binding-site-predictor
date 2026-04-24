# =============================================================================
# 1-try-hocomoco.R — end-to-end runner using the HOCOMOCO motif backend.
# Uses coverage-aware per-motif scoring (s/C) with the default edge guard.
#
# Output layout:
#   caspex_output_hocomoco/
#     extras/
# =============================================================================

source("caspex_analysis_hocomoco.R")   # also sources caspex_analysis.R

# Per-region gRNAs + logFC files live in example_inputs/grnas.tsv.
inputs <- load_caspex_inputs("example_inputs")

result <- run_caspex_hocomoco(
  gene             = "ATP7B",
  grnas            = inputs$grnas,
  data_files       = inputs$data_files,
  out_dir          = "caspex_output_hocomoco",
  hocomoco_version = "v12",                  # "v12" (default) or "v11"
  hocomoco_species = "human",                # "human" (default) or "mouse"
  coverage_correct = TRUE,                   # coverage-aware per-motif scoring
  cov_floor        = 0.05,
  edge_guard_frac  = 0.25
)

source("caspex_extras.R")
run_caspex_extras(result, out_dir = "caspex_output_hocomoco/extras")
