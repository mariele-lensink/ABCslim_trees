#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)

out_root <- Sys.getenv(
  "ABC_CONV_OUTPUT_DIR",
  "/group/gmonroegrp2/mlensink/ABC_data/convergence_test_chrom5_100_abcrf_all_variable"
)
figure_dir <- Sys.getenv(
  "ABC_CONV_PERF_FIGURE_DIR",
  file.path(
    Sys.getenv("HOME"),
    "ABCslim_trees",
    "figures",
    "abc_convergence_chrom5_100_abcrf_all_variable",
    "performance"
  )
)
performance_dir <- Sys.getenv(
  "ABC_CONV_PERF_OUTPUT_DIR",
  file.path(out_root, "performance")
)
priors_file <- Sys.getenv("ABC_CONV_PRIORS_FILE", "priors/priors_2.6.25.csv")

dir.create(performance_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

param_cols <- c("gmu", "imu", "gd", "id", "gdfe", "idfe")
param_labels <- c(
  gmu = "Genic mutation rate",
  imu = "Intergenic mutation rate",
  gd = "Genic selected-site fraction",
  id = "Intergenic selected-site fraction",
  gdfe = "Genic mean DFE effect",
  idfe = "Intergenic mean DFE effect"
)
param_units <- c(
  gmu = "log10",
  imu = "log10",
  gd = "linear",
  id = "linear",
  gdfe = "log10(abs)",
  idfe = "log10(abs)"
)

run_labels <- c(
  rejection_core4 = "Rejection ABC, core 4",
  abcrf_core4 = "ABC-RF, core 4",
  abcrf_all_variable = "ABC-RF, all variable"
)

transform_value <- function(x, parameter) {
  if (parameter %in% c("gmu", "imu")) {
    out <- log10(x)
  } else if (parameter %in% c("gdfe", "idfe")) {
    out <- log10(abs(x))
  } else {
    out <- x
  }
  out[!is.finite(out)] <- NA_real_
  out
}

read_replicate_summaries <- function() {
  combined_path <- file.path(out_root, "posterior_summary_by_replicate.csv")
  if (file.exists(combined_path)) {
    return(fread(combined_path))
  }

  partial_dir <- file.path(out_root, "_replicates")
  files <- list.files(
    partial_dir,
    pattern = "^posterior_summary_replicate_[0-9]+\\.csv$",
    full.names = TRUE
  )
  if (length(files) == 0) {
    stop("No posterior summary files found in ", partial_dir)
  }
  rbindlist(lapply(files, fread), fill = TRUE)
}

read_stat_manifests <- function() {
  combined_path <- file.path(out_root, "summary_statistics_used_by_replicate.csv")
  if (file.exists(combined_path)) {
    return(fread(combined_path))
  }

  partial_dir <- file.path(out_root, "_replicates")
  files <- list.files(
    partial_dir,
    pattern = "^summary_statistics_replicate_[0-9]+\\.csv$",
    full.names = TRUE
  )
  if (length(files) == 0) {
    return(NULL)
  }
  rbindlist(lapply(files, fread), fill = TRUE)
}

make_run <- function(method, stat_set) {
  paste(method, stat_set, sep = "_")
}

add_transformed_columns <- function(dt) {
  dt <- copy(dt)
  dt[, run := make_run(method, stat_set)]
  dt[, run_label := fifelse(run %in% names(run_labels), run_labels[run], run)]
  dt[, parameter_label := param_labels[parameter]]
  dt[, scale := param_units[parameter]]

  for (col in c("true_value", "mean", "median", "q025", "q975")) {
    new_col <- paste0(col, "_t")
    dt[, (new_col) := transform_value(get(col), parameter), by = parameter]
  }
  dt[, mean_error_t := mean_t - true_value_t]
  dt[, median_error_t := median_t - true_value_t]
  dt[, interval_width_t := abs(q975_t - q025_t)]
  dt[, covered_95_t := true_value_t >= pmin(q025_t, q975_t) &
       true_value_t <= pmax(q025_t, q975_t)]
  dt
}

prior_widths <- function() {
  priors <- fread(priors_file)
  missing_params <- setdiff(param_cols, names(priors))
  if (length(missing_params) > 0) {
    stop("Prior file is missing columns: ", paste(missing_params, collapse = ", "))
  }
  rbindlist(lapply(param_cols, function(param) {
    x_t <- transform_value(priors[[param]], param)
    qs <- quantile(x_t, c(0.025, 0.975), na.rm = TRUE)
    data.table(
      parameter = param,
      prior_q025_t = unname(qs[[1]]),
      prior_q975_t = unname(qs[[2]]),
      prior_width_95_t = unname(qs[[2]] - qs[[1]])
    )
  }))
}

scorecard_metrics <- function(dt) {
  widths <- prior_widths()
  metrics <- dt[, {
    ok <- is.finite(true_value_t) & is.finite(median_t)
    fit <- if (sum(ok) >= 3 && sd(true_value_t[ok]) > 0) {
      lm(median_t[ok] ~ true_value_t[ok])
    } else {
      NULL
    }
    data.table(
      n_replicates = uniqueN(replicate),
      median_abs_error = median(abs(median_error_t), na.rm = TRUE),
      mean_abs_error = mean(abs(median_error_t), na.rm = TRUE),
      rmse = sqrt(mean(median_error_t^2, na.rm = TRUE)),
      signed_bias = mean(median_error_t, na.rm = TRUE),
      coverage_95 = mean(covered_95_t, na.rm = TRUE),
      median_interval_width = median(interval_width_t, na.rm = TRUE),
      mean_interval_width = mean(interval_width_t, na.rm = TRUE),
      true_estimate_cor = if (sum(ok) >= 3) cor(true_value_t[ok], median_t[ok]) else NA_real_,
      recovery_r2 = if (!is.null(fit)) summary(fit)$r.squared else NA_real_,
      recovery_slope = if (!is.null(fit)) coef(fit)[[2]] else NA_real_,
      recovery_intercept = if (!is.null(fit)) coef(fit)[[1]] else NA_real_
    )
  }, by = .(run, run_label, method, stat_set, parameter, parameter_label, scale)]
  metrics <- merge(metrics, widths, by = "parameter", all.x = TRUE, sort = FALSE)
  metrics[, posterior_contraction_95 := median_interval_width / prior_width_95_t]
  setcolorder(metrics, c(
    "run", "run_label", "method", "stat_set", "parameter", "parameter_label", "scale",
    setdiff(names(metrics), c("run", "run_label", "method", "stat_set", "parameter", "parameter_label", "scale"))
  ))
  metrics[]
}

threshold_metrics <- function(dt) {
  mutation_like <- dt[parameter %in% c("gmu", "imu", "gdfe", "idfe")]
  fraction_like <- dt[parameter %in% c("gd", "id")]

  fold_dt <- mutation_like[, .(
    within_2x = mean(abs(median_error_t) <= log10(2), na.rm = TRUE),
    within_5x = mean(abs(median_error_t) <= log10(5), na.rm = TRUE),
    within_10x = mean(abs(median_error_t) <= 1, na.rm = TRUE)
  ), by = .(run, run_label, method, stat_set, parameter, parameter_label, scale)]
  fold_long <- melt(
    fold_dt,
    id.vars = c("run", "run_label", "method", "stat_set", "parameter", "parameter_label", "scale"),
    variable.name = "threshold",
    value.name = "fraction_passing"
  )

  abs_dt <- fraction_like[, .(
    within_0.05 = mean(abs(median_error_t) <= 0.05, na.rm = TRUE),
    within_0.10 = mean(abs(median_error_t) <= 0.10, na.rm = TRUE),
    within_0.20 = mean(abs(median_error_t) <= 0.20, na.rm = TRUE)
  ), by = .(run, run_label, method, stat_set, parameter, parameter_label, scale)]
  abs_long <- melt(
    abs_dt,
    id.vars = c("run", "run_label", "method", "stat_set", "parameter", "parameter_label", "scale"),
    variable.name = "threshold",
    value.name = "fraction_passing"
  )

  rbindlist(list(fold_long, abs_long), fill = TRUE)
}

summarize_variable_importance <- function() {
  run_dirs <- list.dirs(out_root, recursive = FALSE, full.names = TRUE)
  files <- unlist(lapply(run_dirs, function(d) {
    list.files(
      d,
      pattern = "^variable_importance_replicate_[0-9]+_observed_.*\\.csv$",
      full.names = TRUE
    )
  }))
  if (length(files) == 0) {
    return(NULL)
  }
  imp <- rbindlist(lapply(files, fread), fill = TRUE)
  if (!"method" %in% names(imp)) {
    imp[, method := "abcrf"]
  }
  imp[, run := make_run(method, stat_set)]
  imp[, run_label := fifelse(run %in% names(run_labels), run_labels[run], run)]
  imp[, parameter_label := param_labels[parameter]]
  imp[, rank_within_replicate := frank(-importance, ties.method = "average"),
      by = .(run, replicate, parameter)]

  imp_summary <- imp[, .(
    n_replicates = uniqueN(replicate),
    mean_importance = mean(importance, na.rm = TRUE),
    median_importance = median(importance, na.rm = TRUE),
    sd_importance = sd(importance, na.rm = TRUE),
    mean_rank = mean(rank_within_replicate, na.rm = TRUE),
    top5_fraction = mean(rank_within_replicate <= 5, na.rm = TRUE)
  ), by = .(run, run_label, stat_set, parameter, parameter_label, statistic)]
  setorder(imp_summary, run, parameter, -mean_importance)
  imp_summary
}

plot_true_vs_estimated <- function(dt, estimate_col = "median") {
  stopifnot(estimate_col %in% c("median", "mean"))
  estimate_t_col <- paste0(estimate_col, "_t")
  estimate_label <- if (estimate_col == "median") "Posterior median" else "Posterior mean"
  file_suffix <- if (estimate_col == "median") "median" else "mean"

  plot_dt <- copy(dt)
  plot_dt[, estimate_t := get(estimate_t_col)]
  r2_dt <- plot_dt[
    is.finite(true_value_t) & is.finite(estimate_t),
    {
      fit <- if (.N >= 3 && sd(true_value_t) > 0) lm(estimate_t ~ true_value_t) else NULL
      data.table(
        x = min(true_value_t, na.rm = TRUE),
        y = max(estimate_t, na.rm = TRUE),
        label = if (!is.null(fit)) sprintf("r = %.3f\nslope = %.3f", cor(true_value_t, estimate_t), coef(fit)[[2]]) else "r = NA\nslope = NA"
      )
    },
    by = .(run_label, parameter_label, scale)
  ]

  p <- ggplot(
    plot_dt,
    aes(x = true_value_t, y = estimate_t, color = run_label)
  ) +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.35, color = "grey35") +
    geom_point(size = 1.25, alpha = 0.75) +
    geom_text(
      data = r2_dt,
      aes(x = x, y = y, label = label, color = run_label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 2.4,
      show.legend = FALSE
    ) +
    facet_wrap(
      ~ parameter_label + scale,
      scales = "free",
      ncol = 3,
      labeller = label_wrap_gen(width = 24)
    ) +
    labs(
      x = "True parameter value, transformed",
      y = paste0(estimate_label, ", transformed"),
      color = "Run"
    ) +
    theme_bw(base_size = 9) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
  ggsave(
    file.path(figure_dir, paste0("true_vs_posterior_", file_suffix, ".png")),
    p,
    width = 8.5,
    height = 6.5,
    dpi = 300
  )
}

plot_error_distributions <- function(dt) {
  p <- ggplot(dt, aes(x = median_error_t, fill = run_label)) +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey25") +
    geom_density(alpha = 0.35, linewidth = 0.35) +
    facet_wrap(
      ~ parameter_label + scale,
      scales = "free",
      ncol = 3,
      labeller = label_wrap_gen(width = 24)
    ) +
    labs(
      x = "Posterior median error, transformed scale",
      y = "Density",
      fill = "Run"
    ) +
    theme_bw(base_size = 9) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
  ggsave(file.path(figure_dir, "posterior_median_error_distributions.png"), p, width = 8.5, height = 6.5, dpi = 300)
}

plot_scorecard <- function(metrics) {
  plot_dt <- melt(
    metrics,
    id.vars = c("run", "run_label", "parameter", "parameter_label", "scale"),
    measure.vars = c("coverage_95", "median_abs_error", "posterior_contraction_95", "recovery_r2"),
    variable.name = "metric",
    value.name = "value"
  )
  plot_dt[, metric_label := factor(
    metric,
    levels = c("coverage_95", "median_abs_error", "posterior_contraction_95", "recovery_r2"),
    labels = c("95% coverage", "Median abs. error", "Posterior / prior width", "Recovery R2")
  )]

  p <- ggplot(plot_dt, aes(x = parameter_label, y = value, fill = run_label)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.65) +
    facet_wrap(~ metric_label, scales = "free_y", ncol = 2) +
    labs(x = NULL, y = "Metric value", fill = "Run") +
    theme_bw(base_size = 9) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 35, hjust = 1),
      panel.grid.minor = element_blank()
    )
  ggsave(file.path(figure_dir, "performance_scorecard.png"), p, width = 9, height = 6.5, dpi = 300)
}

plot_variable_importance <- function(imp_summary) {
  if (is.null(imp_summary) || nrow(imp_summary) == 0) {
    return(invisible(NULL))
  }
  top_stats <- imp_summary[, head(.SD, 10), by = .(run, parameter)]
  top_stats[, normalized_importance := mean_importance / max(mean_importance, na.rm = TRUE),
            by = .(run, parameter)]
  top_stats[, statistic_rank := frank(normalized_importance, ties.method = "first"),
            by = .(run, parameter)]
  top_stats[, statistic_plot := reorder(paste(parameter, statistic, sep = "__"), statistic_rank)]
  stat_labels <- setNames(
    sub("^.*__", "", levels(top_stats$statistic_plot)),
    levels(top_stats$statistic_plot)
  )
  p <- ggplot(top_stats, aes(x = statistic_plot, y = normalized_importance, fill = top5_fraction)) +
    geom_col(width = 0.7) +
    coord_flip() +
    facet_grid(parameter_label ~ run_label, scales = "free_y", space = "free_y") +
    scale_x_discrete(labels = stat_labels) +
    scale_fill_viridis_c(option = "C", limits = c(0, 1)) +
    labs(
      x = NULL,
      y = "Mean variable importance, normalized within parameter",
      fill = "Top-5 fraction"
    ) +
    theme_bw(base_size = 8) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      strip.text.y = element_text(angle = 0)
    )
  ggsave(file.path(figure_dir, "abcrf_variable_importance_top10.png"), p, width = 9.5, height = 9, dpi = 300)
}

summary_dt <- add_transformed_columns(read_replicate_summaries())
summary_dt <- summary_dt[parameter %in% param_cols]
setorder(summary_dt, run, parameter, replicate)

metrics <- scorecard_metrics(summary_dt)
thresholds <- threshold_metrics(summary_dt)
stat_manifest <- read_stat_manifests()
imp_summary <- summarize_variable_importance()

fwrite(summary_dt, file.path(performance_dir, "posterior_summary_by_replicate_transformed.csv"))
fwrite(metrics, file.path(performance_dir, "performance_scorecard.csv"))
fwrite(thresholds, file.path(performance_dir, "threshold_accuracy.csv"))
if (!is.null(stat_manifest)) {
  fwrite(stat_manifest, file.path(performance_dir, "summary_statistics_used_by_replicate.csv"))
}
if (!is.null(imp_summary)) {
  fwrite(imp_summary, file.path(performance_dir, "abcrf_variable_importance_summary.csv"))
}

plot_true_vs_estimated(summary_dt, estimate_col = "median")
plot_true_vs_estimated(summary_dt, estimate_col = "mean")
plot_error_distributions(summary_dt)
plot_scorecard(metrics)
plot_variable_importance(imp_summary)

cat("Wrote convergence performance outputs to", performance_dir, "\n")
cat("Wrote convergence performance figures to", figure_dir, "\n")
