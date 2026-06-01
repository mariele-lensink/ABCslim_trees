#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)

out_root <- Sys.getenv(
  "ABC_CONV_OUTPUT_DIR",
  "../../../group/gmonroegrp2/mlensink/ABC_data/convergence_test"
)
figure_dir <- Sys.getenv(
  "ABC_CONV_FIGURE_DIR",
  file.path(Sys.getenv("HOME"), "figures", "abc_convergence")
)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

param_labels <- c(
  gmu = "genic mutation rate\nlog10",
  imu = "intergenic mutation rate\nlog10",
  gd = "genic selected-site fraction",
  id = "intergenic selected-site fraction",
  gdfe = "genic mean DFE effect\nlog10(abs)",
  idfe = "intergenic mean DFE effect\nlog10(abs)"
)
param_levels <- names(param_labels)
param_row_groups <- c(
  gmu = "Genic",
  imu = "Intergenic",
  gd = "Genic",
  id = "Intergenic",
  gdfe = "Genic",
  idfe = "Intergenic"
)
param_column_groups <- c(
  gmu = "Mutation rate",
  imu = "Mutation rate",
  gd = "Selected-site fraction",
  id = "Selected-site fraction",
  gdfe = "DFE",
  idfe = "DFE"
)
run_names <- c("rejection_core4", "abcrf_core4", "abcrf_all_variable")
run_labels <- c(
  rejection_core4 = "Rejection ABC, core 4 stats",
  abcrf_core4 = "ABC-RF, core 4 stats",
  abcrf_all_variable = "ABC-RF, all variable stats"
)

combine_replicate_outputs <- function() {
  partial_dir <- file.path(out_root, "_replicates")
  summary_files <- list.files(
    partial_dir,
    pattern = "^posterior_summary_replicate_[0-9]+\\.csv$",
    full.names = TRUE
  )
  stat_files <- list.files(
    partial_dir,
    pattern = "^summary_statistics_replicate_[0-9]+\\.csv$",
    full.names = TRUE
  )

  if (length(summary_files) == 0) {
    return(invisible(FALSE))
  }

  message("Combining ", length(summary_files), " replicate posterior summaries from ", partial_dir)
  write_combined <- Sys.getenv("ABC_CONV_WRITE_COMBINED", "1") != "0"
  summary_dt <- rbindlist(lapply(summary_files, fread), fill = TRUE)
  setorder(summary_dt, replicate, method, stat_set, parameter)

  if (length(stat_files) > 0) {
    stat_dt <- rbindlist(lapply(stat_files, fread), fill = TRUE)
    setorder(stat_dt, replicate, method, stat_set, statistic)
    if (write_combined) {
      fwrite(stat_dt, file.path(out_root, "summary_statistics_used_by_replicate.csv"))
    }
  } else {
    warning("No per-replicate summary-statistic manifest files found in ", partial_dir)
  }

  if (write_combined) {
    fwrite(summary_dt, file.path(out_root, "posterior_summary_by_replicate.csv"))
  }

  accuracy <- summary_dt[, .(
    n_replicates = uniqueN(replicate),
    coverage_95 = mean(true_value >= q025 & true_value <= q975, na.rm = TRUE),
    mean_abs_mean_error = mean(abs(mean_error), na.rm = TRUE)
  ), by = .(method, stat_set, parameter)]
  if (write_combined) {
    fwrite(accuracy, file.path(out_root, "convergence_accuracy_summary.csv"))
  }

  summary_dt[, run := paste(method, stat_set, sep = "_")]
  if (write_combined) for (run_name in intersect(unique(summary_dt$run), run_names)) {
    run_dir <- file.path(out_root, run_name)
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
    run_dt <- summary_dt[run == run_name]
    fwrite(run_dt[, !"run"], file.path(run_dir, "posterior_summary_by_replicate.csv"))
    run_accuracy <- accuracy[paste(method, stat_set, sep = "_") == run_name]
    fwrite(run_accuracy, file.path(run_dir, "convergence_accuracy_summary.csv"))
    if (exists("stat_dt")) {
      run_parts <- strsplit(run_name, "_", fixed = TRUE)[[1]]
      run_method <- run_parts[[1]]
      run_stat_set <- paste(run_parts[-1], collapse = "_")
      fwrite(
        stat_dt[method == run_method & stat_set == run_stat_set],
        file.path(run_dir, "summary_statistics_used_by_replicate.csv")
      )
    }
  }

  summary_dt
}

read_run_summary <- function(run_name) {
  path <- file.path(out_root, run_name, "posterior_summary_by_replicate.csv")
  if (!file.exists(path)) {
    warning("Skipping missing summary: ", path)
    return(NULL)
  }
  dt <- fread(path)
  dt[, run := run_name]
  dt[, run_label := run_labels[[run_name]]]
  dt
}

plot_recovery <- function(dt, out_path, title, estimate_col = "mean") {
  stopifnot(estimate_col %in% c("mean", "median"))
  estimate_label <- paste("Inferred", estimate_col)
  true_value_label <- "True value"
  symbol_labels <- c(true_value_label, estimate_label)
  interval_color_values <- c(
    "Posterior interval covers true value" = "grey62",
    "Posterior interval misses true value" = "#b2182b"
  )
  dt <- copy(dt)
  dt <- dt[parameter %in% param_levels]
  dt[, scale_type := fifelse(
    parameter %in% c("gmu", "imu"),
    "log10",
    fifelse(parameter %in% c("gdfe", "idfe"), "neglog10_abs", "linear")
  )]
  for (col in c("true_value", estimate_col, "q025", "q975")) {
    plot_col <- paste0(col, "_plot")
    dt[, (plot_col) := NA_real_]
    dt[scale_type == "linear", (plot_col) := get(col)]
    dt[scale_type == "log10" & get(col) > 0, (plot_col) := log10(get(col))]
    dt[scale_type == "neglog10_abs" & get(col) != 0, (plot_col) := log10(abs(get(col)))]
  }
  setnames(dt, paste0(estimate_col, "_plot"), "estimate_plot")
  dt <- dt[
    is.finite(true_value_plot) &
      is.finite(estimate_plot) &
      is.finite(q025_plot) &
      is.finite(q975_plot)
  ]
  dt[, x_min_plot := pmin(q025_plot, q975_plot)]
  dt[, x_max_plot := pmax(q025_plot, q975_plot)]
  dt[, coverage_status := fifelse(
    true_value_plot >= x_min_plot & true_value_plot <= x_max_plot,
    "Covered",
    "Failed"
  )]
  dt[, coverage_status := factor(
    coverage_status,
    levels = c("Covered", "Failed"),
    labels = c("Posterior interval covers true value", "Posterior interval misses true value")
  )]
  gmu_order <- unique(dt[parameter == "gmu", .(replicate, gmu_true_plot = true_value_plot)])
  setorder(gmu_order, gmu_true_plot, replicate)
  gmu_order[, observed_y := .I]
  dt <- merge(dt, gmu_order[, .(replicate, observed_y)], by = "replicate", all.x = TRUE, sort = FALSE)
  dt <- dt[!is.na(observed_y)]
  dt[, row_group := factor(param_row_groups[parameter], levels = c("Genic", "Intergenic"))]
  dt[, column_group := factor(
    param_column_groups[parameter],
    levels = c("Mutation rate", "Selected-site fraction", "DFE")
  )]
  dt[, parameter := factor(parameter, levels = param_levels, labels = param_labels)]
  setorder(dt, row_group, column_group, observed_y, replicate)
  dt[, observed_label := paste0("rep ", replicate, " / ID ", observed_ID)]

  facet_layer <- if (uniqueN(dt$run_label) > 1) {
    facet_grid(row_group ~ run_label + column_group, scales = "free_x")
  } else {
    facet_grid(row_group ~ column_group, scales = "free_x")
  }

  y_breaks <- sort(unique(dt$observed_y))
  y_labels <- if (length(y_breaks) <= 20) {
    dt[match(y_breaks, observed_y), observed_label]
  } else {
    fifelse(y_breaks %% 10 == 0 | y_breaks == 1, as.character(y_breaks), "")
  }

  p <- ggplot(dt, aes(y = observed_y)) +
    geom_segment(
      aes(x = x_min_plot, xend = x_max_plot, yend = observed_y, color = coverage_status),
      linewidth = 0.65,
      lineend = "round"
    ) +
    geom_point(
      aes(x = true_value_plot, shape = true_value_label),
      size = 0.55,
      color = "#1b9e77"
    ) +
    geom_point(
      aes(x = estimate_plot, shape = estimate_label),
      size = 0.55,
      color = "#111111"
    ) +
    facet_layer +
    scale_y_continuous(
      breaks = y_breaks,
      labels = y_labels,
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_color_manual(
      values = interval_color_values,
      name = "Posterior interval"
    ) +
    guides(
      color = guide_legend(order = 1, override.aes = list(linewidth = 1.1)),
      shape = guide_legend(
        order = 2,
        override.aes = list(color = c("#1b9e77", "#111111"), size = 1.8)
      )
    ) +
    scale_shape_manual(
      values = setNames(c(16, 16), symbol_labels),
      breaks = symbol_labels,
      name = "Figure symbols"
    ) +
    labs(
      title = title,
      x = "Parameter value; mutation-rate facets use log10, DFE facets use log10(abs), selected-site fractions are linear",
      y = "Convergence replicate, ordered by true genic mutation rate"
    ) +
    theme_bw(base_size = 9) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7.5),
      axis.text.y = element_text(size = 5.5),
      axis.text.x = element_text(size = 6.5),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.spacing.x = unit(0.5, "lines"),
      panel.spacing.y = unit(0.35, "lines"),
      strip.background = element_rect(fill = "grey95", color = "grey70")
    )

  ggsave(out_path, p, width = 10.5, height = 5.8, dpi = 300)
}

combined_summaries <- combine_replicate_outputs()
if (is.data.table(combined_summaries)) {
  summaries <- combined_summaries[run %in% run_names]
  summaries[, run_label := run_labels[run]]
} else {
  summaries <- rbindlist(Filter(Negate(is.null), lapply(run_names, read_run_summary)), fill = TRUE)
}
if (nrow(summaries) == 0) {
  stop("No convergence posterior summaries found in ", out_root)
}

for (run_name in unique(summaries$run)) {
  run_dt <- summaries[run == run_name]
  run_fig_dir <- file.path(figure_dir, run_name)
  dir.create(run_fig_dir, recursive = TRUE, showWarnings = FALSE)
  plot_recovery(
    run_dt,
    file.path(run_fig_dir, "posterior_recovery_intervals_mean.png"),
    paste0(run_labels[[run_name]], ": posterior recovery, mean estimate"),
    estimate_col = "mean"
  )
}

plot_recovery(
  summaries,
  file.path(figure_dir, "posterior_recovery_intervals_all_runs_mean.png"),
  "Convergence test posterior recovery, mean estimate",
  estimate_col = "mean"
)

coverage <- summaries[, .(
  n_replicates = uniqueN(replicate),
  coverage_95 = mean(true_value >= q025 & true_value <= q975, na.rm = TRUE),
  mean_abs_mean_error = mean(abs(mean_error), na.rm = TRUE)
), by = .(run, run_label, parameter)]
fwrite(coverage, file.path(figure_dir, "posterior_recovery_plot_summary.csv"))

cat("Saved posterior recovery figures in", figure_dir, "\n")
