#!/usr/bin/env Rscript

# ABCRF Analysis for ABCslim_trees - Chromosome 5
# This script performs Approximate Bayesian Computation with Random Forests
# to infer mutation rate heterogeneity parameters from genomic data.

library(abcrf)
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(patchwork)
library(scales)

# Set working directory
setwd("/home/mlensink/slimsimulations/ABCslim_trees")

# File paths
priors_file <- "priors/priors_2.6.25.csv"
chrom_dir <- Sys.getenv("ABCRF_CHROM", "chrom5")
if (!grepl("^chrom[0-9A-Za-z_]+$", chrom_dir)) {
  stop("ABCRF_CHROM must look like 'chrom1'")
}
chrom_label <- sub("^chrom", "chr", chrom_dir)
sim_stats_file <- file.path("output", "simstats", paste0("sim_stats_", chrom_dir, ".csv"))
target_stats_file <- file.path("output", "ABC_results", "observed", paste0("target_stats_", chrom_dir, ".csv"))
out_dir <- "output/ABC_results"
uncapped_mean_filter <- Sys.getenv("ABC_UNCAPPED_MEAN_FILTER", "0") == "1"
filter_q <- as.numeric(Sys.getenv("ABC_FILTER_Q", "50"))
filter_cap <- as.numeric(Sys.getenv("ABC_FILTER_CAP", "1"))
filter_label <- if (uncapped_mean_filter) paste0("_uncapped_mean_Q", filter_q) else ""
stat_set <- Sys.getenv("ABCRF_STAT_SET", "all_variable")
valid_stat_sets <- c("base", "core", "core_contrast", "core_sfs", "all_variable")
if (!stat_set %in% valid_stat_sets) {
  stop("ABCRF_STAT_SET must be one of: ", paste(valid_stat_sets, collapse = ", "))
}
run_label <- paste0(chrom_label, filter_label, "_", stat_set)
run_dir <- file.path(out_dir, paste0("abcrf_", run_label))
checkpoint_dir <- file.path(run_dir, "checkpoints")
figure_dir <- "figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
ncores <- min(as.integer(Sys.getenv("ABCRF_NCORES", Sys.getenv("SLURM_CPUS_PER_TASK", "4"))), 4)
checkpoint_version <- 2L
cat("ABCRF cores:", ncores, "\n")
cat("ABCRF chromosome:", chrom_dir, "\n")
cat("ABCRF stat set:", stat_set, "\n")
cat("Uncapped mean filter:", uncapped_mean_filter, "\n")
if (uncapped_mean_filter) {
  cat("Keeping rows with abs(gdfe) *", filter_q, "<=", filter_cap,
      "and abs(idfe) *", filter_q, "<=", filter_cap, "\n")
}

# Load data
cat("Loading data...\n")
priors <- read.csv(priors_file)
sim_stats <- read.csv(sim_stats_file)
target_stats <- read.csv(target_stats_file)

safe_ratio <- function(num, den) {
  out <- num / den
  out[!is.finite(out)] <- NA_real_
  out
}

add_contrasts <- function(dt) {
  dt$snp_ratio <- safe_ratio(dt$g_snps, dt$i_snps)
  dt$pi_ratio <- safe_ratio(dt$g_pi, dt$i_pi)
  dt$theta_ratio <- safe_ratio(dt$g_theta_w, dt$i_theta_w)
  dt$td_diff <- dt$g_td - dt$i_td
  dt
}

sim_stats <- add_contrasts(sim_stats)
target_stats <- add_contrasts(target_stats)

# Merge priors with sim stats
ref_table <- merge(priors, sim_stats, by = "ID")
ref_table <- as.data.table(ref_table)

if (uncapped_mean_filter) {
  before_n <- nrow(ref_table)
  ref_table <- ref_table[
    abs(gdfe) * filter_q <= filter_cap &
      abs(idfe) * filter_q <= filter_cap
  ]
  cat("Rows before uncapped-mean filter:", before_n, "\n")
  cat("Rows after uncapped-mean filter:", nrow(ref_table), "\n")
}

# Define parameters to infer
params <- c("gmu", "imu", "gd", "id", "gdfe", "idfe")

# Define summary statistic sets for sensitivity runs. The site-total columns
# are fixed by chromosome interval design, so they carry no inferential signal.
base_stats <- c("g_snps", "i_snps", "g_td", "i_td")
core_stats <- c(
  "g_snps", "i_snps",
  "g_pi", "i_pi",
  "g_theta_w", "i_theta_w",
  "g_td", "i_td",
  "g_singletons", "i_singletons",
  "g_doubletons", "i_doubletons"
)
sfs_stats <- grep("_sfs_", colnames(sim_stats), value = TRUE)
contrast_stats <- c("snp_ratio", "pi_ratio", "theta_ratio", "td_diff")

stat_cols <- switch(
  stat_set,
  base = base_stats,
  core = core_stats,
  core_contrast = c(core_stats, contrast_stats),
  core_sfs = c(core_stats, sfs_stats),
  all_variable = setdiff(colnames(sim_stats), c("ID", "g_sites_total", "i_sites_total"))
)
stat_cols <- intersect(stat_cols, colnames(sim_stats))
stat_cols <- stat_cols[vapply(ref_table[, ..stat_cols], is.numeric, logical(1))]

target_ok <- vapply(target_stats[, stat_cols], function(x) all(is.finite(x)), logical(1))
if (any(!target_ok)) {
  cat("Dropping non-finite observed stats:", paste(stat_cols[!target_ok], collapse = ", "), "\n")
}
stat_cols <- stat_cols[target_ok]

sim_ok <- vapply(ref_table[, ..stat_cols], function(x) all(is.finite(x)), logical(1))
if (any(!sim_ok)) {
  cat("Dropping non-finite simulated stats:", paste(stat_cols[!sim_ok], collapse = ", "), "\n")
}
stat_cols <- stat_cols[sim_ok]

var_ok <- vapply(ref_table[, ..stat_cols], function(x) sd(x, na.rm = TRUE) > 0, logical(1))
if (any(!var_ok)) {
  cat("Dropping zero-variance stats:", paste(stat_cols[!var_ok], collapse = ", "), "\n")
}
stat_cols <- stat_cols[var_ok]
sumstat_names <- stat_cols
sumstats <- ref_table[, ..sumstat_names]
cat("Summary statistics used:", paste(sumstat_names, collapse = ", "), "\n")

stat_set_labels <- c(
  base = "Base: SNP counts and Tajima's D",
  core = "Core: SNP counts, pi, Watterson's theta, Tajima's D, singleton and doubleton counts",
  core_contrast = "Core + contrasts: core stats plus SNP, pi, theta ratios and Tajima's D difference",
  core_sfs = "Core + SFS: core stats plus folded SFS bins",
  all_variable = "All variable: all non-constant simulated summary statistics"
)
stat_legend <- paste0(
  "ABCRF statistics: ", stat_set_labels[[stat_set]],
  " (", length(sumstat_names), " columns). ",
  "Columns used: ", paste(sumstat_names, collapse = ", ")
)

# Observed summary statistics
obs_sumstats <- target_stats[, stat_cols]
obs_sumstats <- as.data.frame(obs_sumstats)

# Train and predict one ABCRF model at a time. Keeping all six forests in
# memory at once can exceed the Slurm memory limit on the full reference table.
cat("Training ABCRF models and predicting posteriors...\n")
posterior_summary <- data.table(Parameter = character(), NMAE = numeric())
weights_table <- data.table(Parameter = character(), Value = list(), Weight = list())
for (param in params) {
  checkpoint_file <- file.path(checkpoint_dir, paste0(param, "_posterior.rds"))
  model_file <- file.path(checkpoint_dir, paste0(param, "_model.rds"))
  model_cols <- c(param, sumstat_names)
  model_data <- as.data.frame(ref_table[, ..model_cols])

  if (file.exists(checkpoint_file)) {
    cat("Loading checkpoint for", param, "\n")
    checkpoint <- readRDS(checkpoint_file)
  } else {
    checkpoint <- NULL
  }

  needs_refresh <- is.null(checkpoint) ||
    is.null(checkpoint$version) ||
    checkpoint$version < checkpoint_version

  if (needs_refresh) {
    if (!is.null(checkpoint)) {
      cat("Refreshing legacy checkpoint for", param, "\n")
    }

    if (file.exists(model_file)) {
      cat("Loading trained model for", param, "\n")
      model <- readRDS(model_file)
    } else {
      cat("Training model for", param, "\n")
      model <- regAbcrf(
        as.formula(paste(param, "~ .")),
        data = model_data,
        ntree = 500,
        paral = TRUE,
        ncores = ncores
      )
      saveRDS(model, model_file)
    }

    cat("Predicting posterior for", param, "\n")
    posterior <- predict(
      model,
      obs = obs_sumstats,
      training = model_data,
      paral = TRUE,
      ncores = ncores,
      rf.weights = TRUE
    )
    checkpoint <- list(
      version = checkpoint_version,
      parameter = param,
      expectation = posterior$expectation,
      median = posterior$med,
      quantiles = posterior$quantiles,
      post.NMAE.mean = posterior$post.NMAE.mean,
      values = model_data[[param]],
      weights = as.numeric(drop(posterior$weights))
    )
    saveRDS(checkpoint, checkpoint_file)
    rm(model, posterior)
    gc()
  }

  posterior_summary <- rbind(
    posterior_summary,
    data.table(Parameter = param, NMAE = checkpoint$post.NMAE.mean)
  )
  weights_table <- rbind(
    weights_table,
    data.table(Parameter = param, Value = list(checkpoint$values), Weight = list(checkpoint$weights))
  )

  rm(model_data, checkpoint)
  gc()
}

# Set factor levels
my_palette <- c("gmu" = "#1F78B4", "imu" = "#A6CEE3",
                "gd" = "#6A3D9A", "id" = "#CAB2D6",
                "gdfe" = "#E31A1C", "idfe" = "#FB9A99")
my_levels <- c("gmu","imu","gd","id","gdfe","idfe")
posterior_summary[, Parameter := factor(Parameter, levels = my_levels)]
weights_table[, Parameter := factor(Parameter, levels = my_levels)]

# Making the prior/posterior figure
plot_info <- data.table(
  param = c("gmu", "imu", "gd", "id", "gdfe", "idfe"),
  label = c("Mutation Rate\n(Genic)(log10)",
            "Mutation Rate\n(Intergenic)(log10)",
            "Proportion of Deleterious Mutations\n(Genic)",
            "Proportion of Deleterious Mutations\n(Intergenic)",
            "Mean Selection Coefficient\n(Genic)(log10)",
            "Mean Selection Coefficient\n(Intergenic)(log10)"),
  color = c("#1F78B4", "#A6CEE3", "#6A3D9A", "#CAB2D6", "#E31A1C", "#FB9A99"),
  log10_x = c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE),
  neglog_x = c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE),
  alpha = c(0.6, 0.8, 0.5, 0.6, 0.4, 0.4)
)
plotting_reftable <- ref_table[, .(gmu, imu, gd, id, gdfe, idfe)]

make_density_plot <- function(param, label, color, log10_x = FALSE, neglog_x = FALSE, alpha = 0.5) {
  prior_dt <- data.table(value = plotting_reftable[[param]], source = "Prior", weight = NA_real_)
  posterior_value <- weights_table[Parameter == param, Value][[1]]
  posterior_weight <- weights_table[Parameter == param, Weight][[1]]
  posterior_dt <- data.table(
    value = posterior_value,
    source = "Posterior",
    weight = posterior_weight
  )
  
  p <- ggplot() +
    geom_density(data = prior_dt, aes(x = value, fill = "Prior"), alpha = 0.3) +
    geom_density(data = posterior_dt, aes(x = value, weight = weight, fill = "Posterior"), alpha = alpha) +
    labs(x = label, y = "Density") +
    scale_fill_manual(values = c("Prior" = "gray", "Posterior" = color)) +
    theme_minimal() +
    theme(
      legend.title = element_blank(),
      legend.position = "none",
      legend.box.background = element_blank(),
      axis.text = element_text(size = 6),
      axis.title = element_text(size = 7)
    )
  
  if (log10_x) {
    p <- p + scale_x_log10()
  }
  if (neglog_x) {
    p <- p + scale_x_continuous(
      trans = trans_new(
        name = "neg_log10_abs",
        transform = function(x) -log10(abs(x)),
        inverse = function(x) -10^(-x)
      ),
      breaks = c(-0.1, -0.01, -0.001)
    )
  }
  
  return(p)
}

plots <- mapply(
  FUN = make_density_plot,
  param = plot_info$param,
  label = plot_info$label,
  color = plot_info$color,
  log10_x = plot_info$log10_x,
  neglog_x = plot_info$neglog_x,
  alpha = plot_info$alpha,
  SIMPLIFY = FALSE
)

combined_plot <- plots$gmu / plots$gd / plots$gdfe / plots$imu / plots$id / plots$idfe +
  plot_layout(ncol = 3) +
  plot_annotation(
    caption = stat_legend,
    theme = theme(
      plot.caption = element_text(size = 8, hjust = 0, lineheight = 1.05, margin = margin(t = 8))
    )
  )

# Save plot
ggsave(file.path(run_dir, paste0("posterior_prior_combined_", run_label, ".png")), combined_plot, width = 12, height = 8)

# Save data
write.csv(posterior_summary, file.path(run_dir, paste0("posterior_summary_", run_label, ".csv")), row.names = FALSE)
write.csv(
  data.frame(
    stat_set = stat_set,
    stat_set_description = unname(stat_set_labels[[stat_set]]),
    statistic = sumstat_names
  ),
  file.path(run_dir, paste0("summary_statistics_", run_label, ".csv")),
  row.names = FALSE
)

cat("Analysis complete. Results saved in", run_dir, "\n")
