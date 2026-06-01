#!/usr/bin/env Rscript

# ABCRF Analysis for ABCslim_trees - Chromosome 5
# This script performs Approximate Bayesian Computation with Random Forests
# to infer mutation rate heterogeneity parameters from genomic data.

library(abcrf)
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(scales)

# Set working directory
project_root <- Sys.getenv("ABC_PROJECT_ROOT", "/home/mlensink/ABCslim_trees")
if (dir.exists(project_root)) {
  setwd(project_root)
}

abc_data_root <- Sys.getenv("ABC_DATA_ROOT", "/group/gmonroegrp2/mlensink/ABC_data")

# File paths
priors_file <- Sys.getenv("ABC_PRIORS_FILE", "priors/priors_2.6.25.csv")
chrom_dir <- Sys.getenv("ABCRF_CHROM", "chrom5")
if (!grepl("^chrom[0-9A-Za-z_]+$", chrom_dir)) {
  stop("ABCRF_CHROM must look like 'chrom1'")
}
chrom_label <- sub("^chrom", "chr", chrom_dir)
stat_version <- Sys.getenv("ABC_STAT_VERSION", "v2")
valid_stat_versions <- c("v1", "v2")
if (!stat_version %in% valid_stat_versions) {
  stop("ABC_STAT_VERSION must be one of: ", paste(valid_stat_versions, collapse = ", "))
}
default_sim_stats_file <- if (stat_version == "v2") {
  file.path(abc_data_root, chrom_dir, paste0("sim_stats_v2_", chrom_dir, "_merged.csv"))
} else {
  file.path(abc_data_root, chrom_dir, "simstats", paste0("sim_stats_", chrom_dir, ".csv"))
}
target_segment_files <- c(
  chrom1 = "chr1_10000000_10284120.sim_stats_v2.csv",
  chrom2 = "chr2_18000000_18280622.sim_stats_v2.csv",
  chrom3 = "chr3_22000000_22349526.sim_stats_v2.csv",
  chrom4 = "chr4_2000000_2660337.sim_stats_v2.csv",
  chrom5 = "chr5_1_377208.sim_stats_v2.csv"
)
default_target_stats_file <- if (stat_version == "v2" && chrom_dir %in% names(target_segment_files)) {
  file.path(abc_data_root, "real_data_segments", "simstats_v2", target_segment_files[[chrom_dir]])
} else {
  file.path(abc_data_root, chrom_dir, "observed", paste0("target_stats_", chrom_dir, ".csv"))
}
sim_stats_file <- Sys.getenv(
  "ABCRF_SIM_STATS_FILE",
  default_sim_stats_file
)
target_stats_file <- Sys.getenv(
  "ABCRF_TARGET_STATS_FILE",
  default_target_stats_file
)
out_dir <- Sys.getenv(
  "ABC_OUTPUT_DIR",
  file.path(abc_data_root, "randomforestABC", if (stat_version == "v2") "ABC_results_statv2" else "ABC_results")
)
uncapped_mean_filter <- Sys.getenv("ABC_UNCAPPED_MEAN_FILTER", "0") == "1"
filter_q <- as.numeric(Sys.getenv("ABC_FILTER_Q", "50"))
filter_cap <- as.numeric(Sys.getenv("ABC_FILTER_CAP", "1"))
filter_label <- if (uncapped_mean_filter) paste0("_uncapped_mean_Q", filter_q) else ""
analysis_label <- Sys.getenv("ABC_ANALYSIS_LABEL", "")
if (nzchar(analysis_label)) {
  filter_label <- paste0(filter_label, "_", analysis_label)
}
stat_set <- Sys.getenv("ABCRF_STAT_SET", "all_variable")
valid_stat_sets <- c("base", "core", "core_contrast", "core_sfs", "all_variable")
if (!stat_set %in% valid_stat_sets) {
  stop("ABCRF_STAT_SET must be one of: ", paste(valid_stat_sets, collapse = ", "))
}
run_label <- paste0(chrom_label, filter_label, "_", stat_set)
run_dir <- file.path(out_dir, paste0("abcrf_", run_label))
checkpoint_dir <- file.path(run_dir, "checkpoints")
figure_dir <- Sys.getenv(
  "ABC_FIGURE_DIR",
  file.path(abc_data_root, "randomforestABC", if (stat_version == "v2") "figures_statv2" else "figures")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
ncores <- min(as.integer(Sys.getenv("ABCRF_NCORES", Sys.getenv("SLURM_CPUS_PER_TASK", "4"))), 4)
checkpoint_version <- 2L
cat("ABCRF cores:", ncores, "\n")
cat("ABCRF chromosome:", chrom_dir, "\n")
cat("ABCRF stat set:", stat_set, "\n")
cat("Target stats file:", target_stats_file, "\n")
cat("Simulation stats file:", sim_stats_file, "\n")
cat("Priors file:", priors_file, "\n")
cat("Analysis label:", ifelse(nzchar(analysis_label), analysis_label, "(none)"), "\n")
cat("ABC stat version:", stat_version, "\n")
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
  if (stat_version == "v2" &&
      all(c("g_snps_perbp", "i_snps_perbp") %in% names(dt))) {
    dt$snp_ratio <- safe_ratio(dt$g_snps_perbp, dt$i_snps_perbp)
  } else {
    dt$snp_ratio <- safe_ratio(dt$g_snps, dt$i_snps)
  }
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
if (stat_version == "v2") {
  snp_stats <- c("g_snps_perbp", "i_snps_perbp")
  singleton_stats <- c("g_singletons_perbp", "i_singletons_perbp")
  doubleton_stats <- c("g_doubletons_perbp", "i_doubletons_perbp")
  excluded_count_stats <- c(
    "g_snps", "i_snps",
    "g_singletons", "i_singletons",
    "g_doubletons", "i_doubletons"
  )
} else {
  snp_stats <- c("g_snps", "i_snps")
  singleton_stats <- c("g_singletons", "i_singletons")
  doubleton_stats <- c("g_doubletons", "i_doubletons")
  excluded_count_stats <- character()
}

base_stats <- c(snp_stats, "g_td", "i_td")
core_stats <- c(
  snp_stats,
  "g_pi", "i_pi",
  "g_theta_w", "i_theta_w",
  "g_td", "i_td",
  singleton_stats,
  doubleton_stats
)
sfs_stats <- grep("_sfs_", colnames(sim_stats), value = TRUE)
contrast_stats <- c("snp_ratio", "pi_ratio", "theta_ratio", "td_diff")

stat_cols <- switch(
  stat_set,
  base = base_stats,
  core = core_stats,
  core_contrast = c(core_stats, contrast_stats),
  core_sfs = c(core_stats, sfs_stats),
  all_variable = setdiff(
    colnames(sim_stats),
    c("ID", "g_sites_total", "i_sites_total", excluded_count_stats)
  )
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
  base = if (stat_version == "v2") "Base: SNPs per bp and Tajima's D" else "Base: SNP counts and Tajima's D",
  core = if (stat_version == "v2") {
    "Core: SNPs, singletons and doubletons per bp, pi, Watterson's theta and Tajima's D"
  } else {
    "Core: SNP counts, pi, Watterson's theta, Tajima's D, singleton and doubleton counts"
  },
  core_contrast = if (stat_version == "v2") {
    "Core + contrasts: statv2 core stats plus SNP-per-bp, pi, theta ratios and Tajima's D difference"
  } else {
    "Core + contrasts: core stats plus SNP, pi, theta ratios and Tajima's D difference"
  },
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

transform_density_value <- function(value, log10_x, neglog_x) {
  if (log10_x) {
    return(log10(value))
  }
  if (neglog_x) {
    return(-log10(abs(value)))
  }
  value
}

combined_density <- rbindlist(lapply(seq_len(nrow(plot_info)), function(i) {
  info <- plot_info[i]
  prior_dt <- data.table(
    value = plotting_reftable[[info$param]],
    source = "Prior",
    weight = 1
  )
  posterior_dt <- data.table(
    value = weights_table[Parameter == info$param, Value][[1]],
    source = "Posterior",
    weight = weights_table[Parameter == info$param, Weight][[1]]
  )
  dt <- rbind(prior_dt, posterior_dt)
  dt[, panel := info$label]
  dt[, plot_value := transform_density_value(value, info$log10_x, info$neglog_x)]
  dt[is.finite(plot_value)]
}))
combined_density[, panel := factor(panel, levels = plot_info$label)]

combined_plot <- ggplot(combined_density, aes(x = plot_value, weight = weight, fill = source)) +
  geom_density(alpha = 0.35, na.rm = TRUE) +
  facet_wrap(~ panel, scales = "free", ncol = 3) +
  scale_fill_manual(values = c("Prior" = "gray", "Posterior" = "#2166ac")) +
  labs(
    x = "Transformed parameter value",
    y = "Density",
    fill = NULL,
    caption = stat_legend
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    plot.caption = element_text(size = 8, hjust = 0, lineheight = 1.05, margin = margin(t = 8)),
    strip.text = element_text(face = "bold", size = 8),
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 8)
  )

# Save prior/posterior plots in the run directory and the shared figure directory.
ggsave(file.path(run_dir, paste0("posterior_prior_combined_", run_label, ".png")), combined_plot, width = 12, height = 8)
ggsave(
  file.path(figure_dir, paste0("abcrf_", run_label, "_prior_posterior_6facet.png")),
  combined_plot,
  width = 12,
  height = 8,
  dpi = 300
)

# Save mutation-rate posteriors alongside the other ABCRF figures.
mutation_weights <- rbindlist(lapply(c("gmu", "imu"), function(param) {
  data.table(
    Parameter = param,
    Value = weights_table[Parameter == param, Value][[1]],
    Weight = weights_table[Parameter == param, Weight][[1]]
  )
}))
mutation_weights[, Parameter := factor(Parameter, levels = c("gmu", "imu"))]
mutation_means <- mutation_weights[, .(
  mean_value = if (sum(Weight, na.rm = TRUE) > 0) {
    weighted.mean(Value, Weight, na.rm = TRUE)
  } else {
    mean(Value, na.rm = TRUE)
  }
), by = Parameter]
mutation_means[, mean_label := label_scientific(digits = 3)(mean_value)]

p_mrate <- ggplot(mutation_weights, aes(x = Value, weight = Weight)) +
  geom_density(fill = "#4393c3", color = "#2166ac", alpha = 0.45, na.rm = TRUE) +
  geom_vline(
    data = mutation_means,
    aes(xintercept = mean_value),
    inherit.aes = FALSE,
    color = "#b2182b",
    linewidth = 0.8,
    linetype = "dashed"
  ) +
  geom_label(
    data = mutation_means,
    aes(x = mean_value, y = Inf, label = paste0("mean = ", mean_label)),
    inherit.aes = FALSE,
    vjust = 1.2,
    hjust = -0.05,
    size = 3,
    label.size = 0.15
  ) +
  facet_grid(Parameter ~ ., scales = "free_y") +
  scale_x_log10(labels = label_scientific(digits = 2)) +
  labs(
    title = "Mutation-rate random-forest posteriors",
    x = "Posterior mutation rate, log10 scale",
    y = "Density"
  ) +
  theme_bw(base_size = 12) +
  theme(strip.text.y = element_text(angle = 0))

ggsave(
  file.path(figure_dir, paste0("abcrf_", run_label, "_posterior_mutation_rates.png")),
  p_mrate,
  width = 8,
  height = 6.5,
  dpi = 300
)

# Save variable-importance figure for the trained ABCRF models.
importance_panels <- c(
  gmu = "Genic\nMutation rate",
  gd = "Genic\nProportion deleterious",
  gdfe = "Genic\nDFE",
  imu = "Intergenic\nMutation rate",
  id = "Intergenic\nProportion deleterious",
  idfe = "Intergenic\nDFE"
)
importance <- list()
for (param in params) {
  model_file <- file.path(checkpoint_dir, paste0(param, "_model.rds"))
  if (!file.exists(model_file)) {
    next
  }
  model_obj <- readRDS(model_file)
  imp <- model_obj$model.rf$variable.importance
  if (is.null(imp) || length(imp) == 0) {
    next
  }
  dt <- data.table(
    Parameter = param,
    statistic = names(imp),
    importance = as.numeric(imp)
  )
  dt <- dt[order(-importance)][seq_len(min(.N, 12))]
  importance[[param]] <- dt
}

if (length(importance) > 0) {
  imp_dt <- rbindlist(importance)
  imp_dt[, importance_panel := factor(
    importance_panels[Parameter],
    levels = importance_panels[c("gmu", "gd", "gdfe", "imu", "id", "idfe")]
  )]
  imp_dt[, statistic_ranked := reorder(statistic, importance)]

  p_importance <- ggplot(imp_dt, aes(x = importance, y = statistic_ranked)) +
    geom_col(fill = "#2166ac", alpha = 0.85, width = 0.72) +
    facet_wrap(~ importance_panel, ncol = 3, scales = "free") +
    labs(
      title = "Random-forest variable importance",
      subtitle = "Top 12 statistics per parameter; each panel has independent x and y scales",
      x = "Ranger impurity importance",
      y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold"),
      axis.text.y = element_text(size = 7),
      panel.spacing.y = unit(1.1, "lines"),
      panel.spacing.x = unit(1.4, "lines")
    )

  ggsave(
    file.path(figure_dir, paste0("abcrf_", run_label, "_variable_importance.png")),
    p_importance,
    width = 14,
    height = 11,
    dpi = 300
  )
} else {
  warning("No variable importance was available for ", run_label)
}

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
