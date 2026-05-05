#!/usr/bin/env Rscript

library(data.table)

setwd("/home/mlensink/slimsimulations/ABCslim_trees")

read_rejection <- function(chrom_label, filtered = FALSE) {
  suffix <- if (filtered) "_uncapped_mean_Q50" else ""
  dir <- file.path("output", "ABC_results", paste0("rejection_", chrom_label, suffix))
  path <- file.path(dir, paste0("posterior_summary_", chrom_label, suffix, "_rejection_tol001.csv"))
  dt <- fread(path)
  dt[, `:=`(
    chromosome = chrom_label,
    method = "rejection",
    analysis = if (filtered) "uncapped_mean_Q50" else "original"
  )]
  dt
}

read_abcrf <- function(chrom_label, stat_set = "all_variable", filtered = FALSE) {
  suffix <- if (filtered) "_uncapped_mean_Q50" else ""
  run_label <- paste0(chrom_label, suffix, "_", stat_set)
  path <- file.path(
    "output", "ABC_results", paste0("abcrf_", run_label),
    paste0("posterior_summary_", run_label, ".csv")
  )
  dt <- fread(path)
  setnames(dt, "Parameter", "parameter")
  dt[, `:=`(
    chromosome = chrom_label,
    method = paste0("abcrf_", stat_set),
    analysis = if (filtered) "uncapped_mean_Q50" else "original"
  )]
  dt
}

rows <- list(
  read_rejection("chr1", filtered = FALSE),
  read_rejection("chr1", filtered = TRUE),
  read_rejection("chr5", filtered = FALSE),
  read_rejection("chr5", filtered = TRUE)
)

for (chrom_label in c("chr1", "chr5")) {
  for (stat_set in c("all_variable")) {
    original_path <- file.path(
      "output", "ABC_results", paste0("abcrf_", chrom_label, "_", stat_set),
      paste0("posterior_summary_", chrom_label, "_", stat_set, ".csv")
    )
    filtered_path <- file.path(
      "output", "ABC_results", paste0("abcrf_", chrom_label, "_uncapped_mean_Q50_", stat_set),
      paste0("posterior_summary_", chrom_label, "_uncapped_mean_Q50_", stat_set, ".csv")
    )
    if (file.exists(original_path) && file.exists(filtered_path)) {
      rows <- c(
        rows,
        list(
          read_abcrf(chrom_label, stat_set = stat_set, filtered = FALSE),
          read_abcrf(chrom_label, stat_set = stat_set, filtered = TRUE)
        )
      )
    }
  }
}

comparison <- rbindlist(rows, fill = TRUE)
setcolorder(comparison, c("method", "analysis", "chromosome", "parameter"))

out <- "output/ABC_results/uncapped_mean_Q50_comparison_summary.csv"
fwrite(comparison, out)
cat("Wrote", out, "\n")
