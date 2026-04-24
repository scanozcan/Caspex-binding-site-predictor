# =============================================================================
# preprocess_regions.R
# Split Busra's single merged proteomics table (cassiopeia_output.txt) into
# four region-specific text files matching the CasPEX input schema.
#
# Region -> condition mapping (inferred from the dep-<cond>-r<N>.xlsx
# filenames shipped alongside the merged table):
#   R1 = DS   (dep-ds-r1)
#   R2 = g5   (dep-g5-r2)
#   R3 = g3   (dep-g3-r3)
#   R4 = g1   (dep-g1-r4)
#
# Output columns per region file (tab-separated, header row):
#   name, logFC, ave_exp, P.Value_reg<N>__vs__NT, minus_log_pval, TFDatabase
#
# TFDatabase takes values {exist, absent} based on membership in
# TFLibrary.txt (Feng/Zheng library) sitting next to this script.
#
# Usage (from TF_analysis/):
#   Rscript input/preprocess_regions.R
# or from inside the input/ folder:
#   setwd(".../TF_analysis/input"); source("preprocess_regions.R")
# =============================================================================

# Resolve paths relative to *this* script so the driver can be invoked from
# anywhere.
this_dir <- tryCatch(
  dirname(sys.frame(1)$ofile),
  error = function(e) getwd()
)
if (is.null(this_dir) || !nzchar(this_dir)) this_dir <- getwd()

src_file <- file.path(this_dir, "can_busra_caspex_analysis_cassiopeia_output.txt")
tf_file  <- file.path(this_dir, "TFLibrary.txt")

if (!file.exists(src_file))
  stop("merged table not found: ", src_file)
if (!file.exists(tf_file))
  stop("TFLibrary.txt not found next to this script: ", tf_file)

# Region -> (condition token, 1-based region index)
mapping <- list(
  R1 = list(cond = "DS", idx = 1L),
  R2 = list(cond = "g5", idx = 2L),
  R3 = list(cond = "g3", idx = 3L),
  R4 = list(cond = "g1", idx = 4L)
)

# --- load TF library (one symbol per line, strip CR/whitespace) -----------
tf_lines <- readLines(tf_file, warn = FALSE)
tf_set   <- unique(trimws(sub("\r$", "", tf_lines)))
tf_set   <- tf_set[nzchar(tf_set)]
message("TF reference: ", length(tf_set), " unique symbols")

# --- read the merged table ------------------------------------------------
merged <- read.delim(src_file, check.names = FALSE, stringsAsFactors = FALSE)
message("Merged table: ", nrow(merged), " rows x ", ncol(merged), " cols")

# primary gene symbol = first token before ';' in "Gene names"; fall back
# to first Protein ID token if blank (matches MaxQuant output conventions)
primary <- function(s) sub(";.*$", "", trimws(s))
name_raw <- primary(merged[["Gene names"]])
fallback <- primary(merged[["Protein IDs"]])
name_raw[!nzchar(name_raw)] <- fallback[!nzchar(name_raw)]

keep_name <- nzchar(name_raw)
is_tf     <- ifelse(name_raw %in% tf_set, "exist", "absent")

# --- write one Region<N>.txt per mapping entry ---------------------------
for (region in names(mapping)) {
  cond <- mapping[[region]]$cond
  n    <- mapping[[region]]$idx

  lfc  <- suppressWarnings(as.numeric(merged[[paste0("logFC_",   cond, "__vs__NT")]]))
  ave  <- suppressWarnings(as.numeric(merged[[paste0("AveExpr_", cond, "__vs__NT")]]))
  pval <- suppressWarnings(as.numeric(merged[[paste0("P.Value_", cond, "__vs__NT")]]))

  ok <- keep_name & is.finite(lfc) & is.finite(pval)

  # Clamp pval = 0 to 1e-300 so -log10 stays finite; matches the upstream
  # limma convention where tiny p-values occasionally get under/overflowed.
  p_clamped <- pmax(pval[ok], 1e-300)

  out <- data.frame(
    name          = name_raw[ok],
    logFC         = lfc[ok],
    ave_exp       = ave[ok],
    P.Value       = pval[ok],
    minus_log_pval = -log10(p_clamped),
    TFDatabase    = is_tf[ok],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  names(out)[names(out) == "P.Value"] <- paste0("P.Value_reg", n, "__vs__NT")

  out_path <- file.path(this_dir, paste0("Region", n, ".txt"))
  write.table(out, out_path, sep = "\t", quote = FALSE,
              row.names = FALSE, na = "")
  message(sprintf("%s  (cond %-2s)  ->  %s  | rows=%5d  TFs=%4d  dropped=%d",
                  region, cond, basename(out_path),
                  sum(ok), sum(out$TFDatabase == "exist"), sum(!ok)))
}

message("Done. These Region*.txt files are ready to be picked up by ",
        "load_caspex_inputs(\"input\") together with busra_grnas.tsv.")
