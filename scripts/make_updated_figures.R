#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(scales)
library(grid)

project_root <- Sys.getenv("ABC_PROJECT_ROOT", "/home/mlensink/ABCslim_trees")
if (dir.exists(project_root)) {
  setwd(project_root)
}

abc_data_root <- Sys.getenv("ABC_DATA_ROOT", "/group/gmonroegrp2/mlensink/ABC_data")
abc_results_root <- Sys.getenv(
  "ABC_RESULTS_ROOT",
  file.path(abc_data_root, "randomforestABC", "ABC_results_statv2")
)
out_root <- Sys.getenv(
  "ABC_UPDATED_FIGURE_DIR",
  file.path(abc_data_root, "randomforestABC", "figures_updated_statv2")
)
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "rejection"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "abcrf"), recursive = TRUE, showWarnings = FALSE)

param_cols <- c("gmu", "imu", "gd", "id", "gdfe", "idfe")
param_labels <- c(
  gmu = "Genic mutation rate",
  imu = "Intergenic mutation rate",
  gd = "Genic deleterious proportion",
  id = "Intergenic deleterious proportion",
  gdfe = "Genic DFE mean",
  idfe = "Intergenic DFE mean"
)
mutation_params <- c("gmu", "imu")

param_plot_info <- data.table(
  parameter_raw = param_cols,
  row_group = factor(
    c("Genic", "Intergenic", "Genic", "Intergenic", "Genic", "Intergenic"),
    levels = c("Genic", "Intergenic")
  ),
  column_group = factor(
    c("Mutation rate", "Mutation rate", "Proportion deleterious",
      "Proportion deleterious", "DFE", "DFE"),
    levels = c("Mutation rate", "Proportion deleterious", "DFE")
  ),
  panel_label = factor(
    param_labels[param_cols],
    levels = param_labels[c("gmu", "imu", "gd", "id", "gdfe", "idfe")]
  ),
  transform_type = c("log10", "log10", "identity", "identity", "neglog10_abs", "neglog10_abs")
)
param_plot_info[, importance_panel := factor(
  paste(row_group, column_group, sep = "\n"),
  levels = c(
    "Genic\nMutation rate",
    "Genic\nProportion deleterious",
    "Genic\nDFE",
    "Intergenic\nMutation rate",
    "Intergenic\nProportion deleterious",
    "Intergenic\nDFE"
  )
)]

theme_set(theme_bw(base_size = 12))

read_if_exists <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  fread(path)
}

case_chrom <- function(case_name) {
  hit <- regmatches(case_name, regexpr("chr[0-9]+", case_name))
  if (length(hit) == 0 || hit == "") {
    return(NA_character_)
  }
  hit
}

case_uses_uncapped_filter <- function(case_name) {
  grepl("uncapped_mean_Q50", case_name)
}

simstats_for_chrom <- function(chrom_label) {
  chrom_dir <- sub("^chr", "chrom", chrom_label)
  file.path(abc_data_root, chrom_dir, paste0("sim_stats_v2_", chrom_dir, "_merged.csv"))
}

target_segment_files <- c(
  chr1 = "chr1_10000000_10284120.sim_stats_v2.csv",
  chr2 = "chr2_18000000_18280622.sim_stats_v2.csv",
  chr3 = "chr3_22000000_22349526.sim_stats_v2.csv",
  chr4 = "chr4_2000000_2660337.sim_stats_v2.csv",
  chr5 = "chr5_1_377208.sim_stats_v2.csv"
)

target_for_chrom <- function(chrom_label) {
  if (!chrom_label %in% names(target_segment_files)) {
    stop("No all-sites 1001 target segment is configured for ", chrom_label)
  }
  file.path(abc_data_root, "real_data_segments", "simstats_v2", target_segment_files[[chrom_label]])
}

prior_for_case <- function(case_name) {
  chrom_label <- case_chrom(case_name)
  sim_path <- simstats_for_chrom(chrom_label)
  if (!file.exists(sim_path)) {
    warning("No simstats found for ", case_name, ": ", sim_path)
    return(NULL)
  }
  priors <- fread("priors/priors_2.6.25.csv")
  sim_stats <- fread(sim_path, select = "ID")
  par <- priors[match(sim_stats$ID, ID)]
  par <- par[!is.na(ID)]
  if (case_uses_uncapped_filter(case_name)) {
    par <- par[abs(gdfe) * 50 <= 1 & abs(idfe) * 50 <= 1]
  }
  par[, ..param_cols]
}

long_params <- function(dt, type, weight_col = NULL) {
  keep_cols <- c(param_cols, weight_col)
  keep_cols <- keep_cols[!is.na(keep_cols)]
  long <- melt(
    dt[, ..keep_cols],
    measure.vars = param_cols,
    variable.name = "parameter_raw",
    value.name = "value"
  )
  long[, type := type]
  add_param_plot_fields(long)
}

transform_param_value <- function(value, transform_type) {
  ifelse(
    transform_type == "log10",
    log10(value),
    ifelse(transform_type == "neglog10_abs", -log10(abs(value)), value)
  )
}

add_param_plot_fields <- function(dt) {
  dt <- merge(dt, param_plot_info, by = "parameter_raw", all.x = TRUE, sort = FALSE)
  dt[, plot_value := transform_param_value(value, transform_type)]
  dt <- dt[is.finite(plot_value)]
  dt
}

weighted_median <- function(value, weight) {
  ok <- is.finite(value) & is.finite(weight) & weight > 0
  value <- value[ok]
  weight <- weight[ok]
  if (length(value) == 0 || sum(weight) <= 0) {
    return(NA_real_)
  }
  ord <- order(value)
  value <- value[ord]
  weight <- weight[ord]
  value[which(cumsum(weight) / sum(weight) >= 0.5)[1]]
}

summary_value_table <- function(long_dt, weight_col = NULL, summary_label = "mean") {
  if (is.null(weight_col)) {
    out <- long_dt[, .(summary_value = mean(value, na.rm = TRUE)), by = .(
      parameter_raw, row_group, column_group, panel_label, transform_type
    )]
  } else {
    out <- long_dt[, .(
      summary_value = if (sum(get(weight_col), na.rm = TRUE) > 0) {
        weighted.mean(value, get(weight_col), na.rm = TRUE)
      } else {
        mean(value, na.rm = TRUE)
      }
    ), by = .(parameter_raw, row_group, column_group, panel_label, transform_type)]
  }
  out[, summary_plot_value := transform_param_value(summary_value, transform_type)]
  out <- out[is.finite(summary_plot_value)]
  out[, summary_label := summary_label]
  out[, summary_value_label := label_scientific(digits = 3)(summary_value)]
  out
}

save_plot <- function(path, plot, width = 12, height = 8) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggsave(path, plot, width = width, height = height, dpi = 300)
  message("wrote ", path)
}

plot_prior_posterior <- function(prior_long, posterior_long, posterior_summary, title, out_path) {
  plot_dt <- rbindlist(list(prior_long, posterior_long), fill = TRUE)
  if (!("weight" %in% names(plot_dt))) {
    plot_dt[, weight := 1]
  }
  plot_dt[is.na(weight), weight := 1]
  p <- ggplot(plot_dt, aes(x = plot_value, weight = weight, fill = type, color = type)) +
    geom_density(alpha = 0.35, adjust = 1, na.rm = TRUE) +
    geom_vline(
      data = posterior_summary,
      aes(xintercept = summary_plot_value),
      color = "#b2182b",
      linewidth = 0.8,
      linetype = "dashed"
    ) +
    geom_label(
      data = posterior_summary,
      aes(x = summary_plot_value, y = Inf, label = paste0("mean = ", summary_value_label)),
      inherit.aes = FALSE,
      vjust = 1.2,
      hjust = -0.05,
      size = 3,
      label.size = 0.15
    ) +
    facet_grid(row_group ~ column_group, scales = "free_x") +
    scale_fill_manual(values = c(Prior = "grey70", Posterior = "#2166ac")) +
    scale_color_manual(values = c(Prior = "grey45", Posterior = "#2166ac")) +
    labs(
      title = title,
      x = "Transformed parameter value: log10(mu), proportion deleterious, -log10(abs(DFE mean))",
      y = "Density",
      fill = NULL,
      color = NULL
    ) +
    theme(
      legend.position = "top",
      strip.text = element_text(face = "bold")
    )
  save_plot(out_path, p, width = 13, height = 8.5)
}

save_column_plots <- function(path, plots, width = 15, height = 8.5) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  png(path, width = width, height = height, units = "in", res = 300)
  on.exit(dev.off(), add = TRUE)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(1, length(plots))))
  for (i in seq_along(plots)) {
    print(plots[[i]], vp = viewport(layout.pos.row = 1, layout.pos.col = i))
  }
  message("wrote ", path)
}

weighted_welch_test <- function(x, y, wx = NULL, wy = NULL) {
  weighted_stats <- function(value, weight = NULL) {
    if (is.null(weight)) {
      ok <- is.finite(value)
      value <- value[ok]
      return(list(mean = mean(value), variance = var(value), n_eff = length(value)))
    }
    ok <- is.finite(value) & is.finite(weight) & weight > 0
    value <- value[ok]
    weight <- weight[ok]
    norm_weight <- weight / sum(weight)
    mean_value <- sum(norm_weight * value)
    list(
      mean = mean_value,
      variance = sum(norm_weight * (value - mean_value)^2) / (1 - sum(norm_weight^2)),
      n_eff = sum(weight)^2 / sum(weight^2)
    )
  }

  x_stats <- weighted_stats(x, wx)
  y_stats <- weighted_stats(y, wy)
  se <- sqrt(x_stats$variance / x_stats$n_eff + y_stats$variance / y_stats$n_eff)
  t_stat <- (x_stats$mean - y_stats$mean) / se
  df <- (x_stats$variance / x_stats$n_eff + y_stats$variance / y_stats$n_eff)^2 /
    ((x_stats$variance / x_stats$n_eff)^2 / (x_stats$n_eff - 1) +
       (y_stats$variance / y_stats$n_eff)^2 / (y_stats$n_eff - 1))
  p_value <- 2 * pt(-abs(t_stat), df)
  list(t = t_stat, df = df, p = p_value)
}

mutation_rate_test_label <- function(mut, weight_col = NULL) {
  g <- mut[parameter_raw == "gmu"]
  i <- mut[parameter_raw == "imu"]
  if (nrow(g) == 0 || nrow(i) == 0) {
    return(NULL)
  }
  wx <- if (is.null(weight_col)) NULL else g[[weight_col]]
  wy <- if (is.null(weight_col)) NULL else i[[weight_col]]
  test <- weighted_welch_test(g$value, i$value, wx, wy)
  p_label <- if (is.finite(test$p) && test$p < 0.001) {
    "p < 0.001"
  } else {
    paste0("p = ", label_pvalue(accuracy = 0.001)(test$p))
  }
  paste0("Weighted Welch test for genic vs intergenic posterior mean mutation rates: ", p_label)
}

plot_prior_posterior_column_free_y <- function(prior_long, posterior_long, posterior_summary, title, out_path) {
  plot_dt <- rbindlist(list(prior_long, posterior_long), fill = TRUE)
  if (!("weight" %in% names(plot_dt))) {
    plot_dt[, weight := 1]
  }
  plot_dt[is.na(weight), weight := 1]
  plot_dt[, density_group := fifelse(type == "Prior", "Prior", as.character(column_group))]
  plot_dt[, density_group := factor(density_group, levels = c("Prior", "Mutation rate", "Proportion deleterious", "DFE"))]
  density_colors <- c(
    Prior = "grey70",
    `Mutation rate` = "#1F78B4",
    `Proportion deleterious` = "#6A3D9A",
    DFE = "#E31A1C"
  )

  make_type_plot <- function(column_name, plot_title, x_label) {
    dt <- plot_dt[column_group == column_name]
    med <- copy(posterior_summary[column_group == column_name])
    panel_ranges <- dt[, .(xmin = min(plot_value, na.rm = TRUE), xmax = max(plot_value, na.rm = TRUE)), by = row_group]
    med <- merge(med, panel_ranges, by = "row_group", all.x = TRUE, sort = FALSE)
    med[, label_hjust := fifelse(
      column_group == "Mutation rate",
      -0.05,
      fifelse(summary_plot_value < (xmin + xmax)/2, -0.05, 1.05)
    )]
    x_expand <- if (column_name == "Mutation rate") c(0.06, 0.30) else c(0.20, 0.20)
    ggplot(dt, aes(x = plot_value, weight = weight, fill = density_group, color = density_group)) +
      geom_density(alpha = 0.40, adjust = 1, na.rm = TRUE) +
      geom_vline(
        data = med,
        aes(xintercept = summary_plot_value),
        color = "#b2182b",
        linewidth = 0.9,
        linetype = "dashed"
      ) +
      geom_label(
        data = med,
        aes(x = summary_plot_value, y = Inf, label = paste0("mean = ", summary_value_label), hjust = label_hjust),
        inherit.aes = FALSE,
        vjust = 1.35,
        size = 4.2,
        label.size = 0.15
      ) +
      facet_grid(rows = vars(row_group), scales = "fixed") +
      scale_fill_manual(values = density_colors, drop = FALSE) +
      scale_color_manual(values = density_colors, drop = FALSE) +
      scale_x_continuous(expand = expansion(mult = x_expand)) +
      labs(title = plot_title, x = x_label, y = "Density", fill = NULL, color = NULL) +
      coord_cartesian(clip = "off") +
      theme_minimal(base_size = 15) +
      theme(
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 16),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 13),
        legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5, size = 17),
        plot.margin = margin(5.5, 11, 5.5, 11)
      )
  }

  plots <- list(
    make_type_plot("Mutation rate", "Mutation rates", "log10(mu)"),
    make_type_plot("Proportion deleterious", "Fraction deleterious", "Proportion deleterious"),
    make_type_plot("DFE", "Strength of deleterious effects", "-log10(abs(DFE mean))")
  )
  save_column_plots(out_path, plots, width = 13.2, height = 6.8)
}

plot_mutation_posteriors <- function(posterior_long, posterior_summary, title, out_path, weight_col = NULL) {
  mut <- posterior_long[parameter_raw %in% mutation_params]
  mut_summary <- posterior_summary[parameter_raw %in% mutation_params]
  test_label <- mutation_rate_test_label(mut, weight_col)
  mapping <- aes(x = value)
  if (!is.null(weight_col)) {
    mapping <- aes(x = value, weight = get(weight_col))
  }
  p <- ggplot(mut, mapping) +
    geom_density(fill = "#4393c3", color = "#2166ac", alpha = 0.45, na.rm = TRUE) +
    geom_vline(
      data = mut_summary,
      aes(xintercept = summary_value),
      color = "#b2182b",
      linewidth = 0.8,
      linetype = "dashed"
    ) +
    geom_label(
      data = mut_summary,
      aes(x = summary_value, y = Inf, label = paste0("mean = ", summary_value_label)),
      inherit.aes = FALSE,
      vjust = 1.2,
      hjust = -0.05,
      size = 3,
      label.size = 0.15
    ) +
    facet_grid(row_group ~ ., scales = "free_y") +
    scale_x_log10(labels = label_scientific(digits = 2)) +
    labs(
      title = title,
      subtitle = test_label,
      x = "Posterior mutation rate, log10 scale",
      y = "Density"
    ) +
    theme(
      strip.text.y = element_text(angle = 0, face = "bold", size = 13),
      plot.subtitle = element_text(size = 11),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 11)
    )
  save_plot(out_path, p, width = 8, height = 6.5)
}

plot_rejection_stat_histograms <- function(case_name, post, out_path) {
  chrom_label <- case_chrom(case_name)
  sim_path <- simstats_for_chrom(chrom_label)
  target_path <- target_for_chrom(chrom_label)
  sim_stats <- read_if_exists(sim_path)
  target <- read_if_exists(target_path)
  if (is.null(sim_stats) || is.null(target)) {
    warning("Skipping posterior-stat histograms for ", case_name, ": missing simstats or target stats")
    return(invisible(FALSE))
  }

  post_stats <- sim_stats[post[, .(ID)], on = "ID", nomatch = 0]
  if (nrow(post_stats) == 0) {
    warning("Skipping posterior-stat histograms for ", case_name, ": no posterior IDs matched simstats")
    return(invisible(FALSE))
  }

  stat_cols <- setdiff(names(post_stats), c("ID", "g_sites_total", "i_sites_total"))
  stat_cols <- stat_cols[vapply(post_stats[, ..stat_cols], is.numeric, logical(1))]
  long_stats <- melt(
    post_stats[, c("ID", stat_cols), with = FALSE],
    id.vars = "ID",
    variable.name = "statistic",
    value.name = "value"
  )
  observed <- melt(
    target[, ..stat_cols],
    measure.vars = stat_cols,
    variable.name = "statistic",
    value.name = "observed"
  )

  p <- ggplot(long_stats, aes(x = value)) +
    geom_histogram(bins = 40, fill = "grey70", color = "white", na.rm = TRUE) +
    geom_vline(
      data = observed,
      aes(xintercept = observed),
      color = "#b2182b",
      linewidth = 0.7
    ) +
    facet_wrap(~ statistic, scales = "free", ncol = 4) +
    scale_x_continuous(labels = label_number()) +
    labs(
      title = paste0(case_name, ": accepted rejection posterior summary statistics"),
      subtitle = "Red vertical line marks observed 1001 Genomes target value",
      x = "Posterior accepted simulation statistic value",
      y = "Accepted simulations"
    )
  save_plot(out_path, p, width = 16, height = 18)
}

process_rejection_case <- function(case_dir) {
  case_name <- basename(case_dir)
  sample_file <- list.files(case_dir, pattern = "^posterior_samples_.*rejection_tol001\\.csv$", full.names = TRUE)
  if (length(sample_file) == 0) {
    warning("Skipping rejection case with no posterior samples: ", case_name)
    return(invisible(FALSE))
  }
  post <- fread(sample_file[[1]])
  missing_params <- setdiff(param_cols, names(post))
  if (length(missing_params) > 0) {
    warning("Skipping rejection case with missing parameter columns: ", case_name)
    return(invisible(FALSE))
  }

  prior <- prior_for_case(case_name)
  if (is.null(prior)) {
    return(invisible(FALSE))
  }

  prior_long <- long_params(prior, "Prior")
  post_long <- long_params(post, "Posterior")
  post_summary <- summary_value_table(post_long)

  plot_prior_posterior(
    prior_long,
    post_long,
    post_summary,
    paste0(case_name, ": prior and rejection posterior"),
    file.path(out_root, "rejection", paste0(case_name, "_prior_posterior_6facet.png"))
  )
  plot_prior_posterior_column_free_y(
    prior_long,
    post_long,
    post_summary,
    paste0(case_name, ": prior and rejection posterior"),
    file.path(out_root, "rejection", paste0(case_name, "_prior_posterior_6facet_column_free_y.png"))
  )
  plot_mutation_posteriors(
    post_long,
    post_summary,
    paste0(case_name, ": mutation-rate rejection posteriors"),
    file.path(out_root, "rejection", paste0(case_name, "_posterior_mutation_rates.png"))
  )
  plot_rejection_stat_histograms(
    case_name,
    post,
    file.path(out_root, "rejection", paste0(case_name, "_posterior_summary_stats_histograms.png"))
  )
  invisible(TRUE)
}

rf_posterior_long <- function(checkpoint_dir) {
  pieces <- list()
  means <- list()
  for (param in param_cols) {
    path <- file.path(checkpoint_dir, paste0(param, "_posterior.rds"))
    if (!file.exists(path)) {
      return(NULL)
    }
    obj <- readRDS(path)
    dt <- data.table(
      parameter_raw = param,
      value = as.numeric(obj$values),
      weight = as.numeric(drop(obj$weights))
    )
    dt <- add_param_plot_fields(dt)
    pieces[[param]] <- dt
    summary_value <- if (!is.null(obj$mean)) {
      as.numeric(obj$mean)
    } else if (sum(dt$weight, na.rm = TRUE) > 0) {
      weighted.mean(dt$value, dt$weight, na.rm = TRUE)
    } else {
      mean(dt$value, na.rm = TRUE)
    }
    means[[param]] <- data.table(
      parameter_raw = param,
      row_group = param_plot_info[parameter_raw == param, row_group],
      column_group = param_plot_info[parameter_raw == param, column_group],
      panel_label = param_plot_info[parameter_raw == param, panel_label],
      transform_type = param_plot_info[parameter_raw == param, transform_type],
      summary_value = summary_value,
      summary_plot_value = transform_param_value(summary_value, param_plot_info[parameter_raw == param, transform_type]),
      summary_label = "mean",
      summary_value_label = label_scientific(digits = 3)(summary_value)
    )
  }
  list(post = rbindlist(pieces), summary = rbindlist(means))
}

plot_rf_variable_importance <- function(case_name, checkpoint_dir, out_path) {
  importance <- list()
  for (param in param_cols) {
    path <- file.path(checkpoint_dir, paste0(param, "_model.rds"))
    if (!file.exists(path)) {
      next
    }
    obj <- readRDS(path)
    imp <- obj$model.rf$variable.importance
    if (is.null(imp)) {
      next
    }
    dt <- data.table(
      parameter = param,
      statistic = names(imp),
      importance = as.numeric(imp)
    )
    dt <- dt[order(-importance)][seq_len(min(.N, 12))]
    importance[[param]] <- dt
  }
  if (length(importance) == 0) {
    warning("No variable importance available for ", case_name)
    return(invisible(FALSE))
  }
  imp_dt <- rbindlist(importance)
  setnames(imp_dt, "parameter", "parameter_raw")
  imp_dt <- merge(imp_dt, param_plot_info, by = "parameter_raw", all.x = TRUE, sort = FALSE)
  imp_dt[, statistic_ranked := reorder(statistic, importance)]

  p <- ggplot(imp_dt, aes(x = importance, y = statistic_ranked)) +
    geom_col(fill = "#2166ac", alpha = 0.85, width = 0.72) +
    facet_wrap(~ importance_panel, ncol = 3, scales = "free") +
    labs(
      title = paste0(case_name, ": random-forest variable importance"),
      subtitle = "Top 12 statistics per parameter; each panel has independent x and y scales",
      x = "Ranger impurity importance",
      y = NULL
    ) +
    theme(
      strip.text = element_text(face = "bold"),
      axis.text.y = element_text(size = 7),
      panel.spacing.y = unit(1.1, "lines"),
      panel.spacing.x = unit(1.4, "lines")
    )
  save_plot(out_path, p, width = 14, height = 11)
}

process_abcrf_case <- function(case_dir) {
  case_name <- basename(case_dir)
  checkpoint_dir <- file.path(case_dir, "checkpoints")
  if (!dir.exists(checkpoint_dir)) {
    warning("Skipping ABCRF case with no checkpoints directory: ", case_name)
    return(invisible(FALSE))
  }
  rf <- rf_posterior_long(checkpoint_dir)
  if (is.null(rf)) {
    warning("Skipping ABCRF posterior plots with incomplete posterior RDS files: ", case_name)
    return(invisible(FALSE))
  }

  prior_long <- copy(rf$post)
  prior_long[, type := "Prior"]
  prior_long[, weight := NULL]
  post_long <- copy(rf$post)
  post_long[, type := "Posterior"]

  plot_prior_posterior(
    prior_long,
    post_long,
    rf$summary,
    paste0(case_name, ": prior and random-forest posterior"),
    file.path(out_root, "abcrf", paste0(case_name, "_prior_posterior_6facet.png"))
  )
  plot_prior_posterior_column_free_y(
    prior_long,
    post_long,
    rf$summary,
    paste0(case_name, ": prior and random-forest posterior"),
    file.path(out_root, "abcrf", paste0(case_name, "_prior_posterior_6facet_column_free_y.png"))
  )
  plot_mutation_posteriors(
    post_long,
    rf$summary,
    paste0(case_name, ": mutation-rate random-forest posteriors"),
    file.path(out_root, "abcrf", paste0(case_name, "_posterior_mutation_rates.png")),
    weight_col = "weight"
  )
  plot_rf_variable_importance(
    case_name,
    checkpoint_dir,
    file.path(out_root, "abcrf", paste0(case_name, "_variable_importance.png"))
  )
  invisible(TRUE)
}

abc_dirs <- list.dirs(abc_results_root, recursive = FALSE, full.names = TRUE)
case_grep <- Sys.getenv("ABC_CASE_GREP", "")
if (nzchar(case_grep)) {
  abc_dirs <- abc_dirs[grepl(case_grep, basename(abc_dirs))]
}
rejection_dirs <- abc_dirs[grepl("/rejection_", abc_dirs)]
abcrf_dirs <- abc_dirs[grepl("/abcrf_", abc_dirs) & dir.exists(file.path(abc_dirs, "checkpoints"))]
only_variable_importance <- Sys.getenv("ONLY_VARIABLE_IMPORTANCE", "0") == "1"

if (!only_variable_importance) {
  message("Rejection cases: ", length(rejection_dirs))
  for (case_dir in rejection_dirs) {
    message("Processing ", basename(case_dir))
    tryCatch(process_rejection_case(case_dir), error = function(e) warning(basename(case_dir), ": ", conditionMessage(e)))
  }
}

message("ABCRF cases: ", length(abcrf_dirs))
for (case_dir in abcrf_dirs) {
  message("Processing ", basename(case_dir))
  if (only_variable_importance) {
    checkpoint_dir <- file.path(case_dir, "checkpoints")
    tryCatch(
      plot_rf_variable_importance(
        basename(case_dir),
        checkpoint_dir,
        file.path(out_root, "abcrf", paste0(basename(case_dir), "_variable_importance.png"))
      ),
      error = function(e) warning(basename(case_dir), ": ", conditionMessage(e))
    )
  } else {
    tryCatch(process_abcrf_case(case_dir), error = function(e) warning(basename(case_dir), ": ", conditionMessage(e)))
  }
}

message("Done. Figures are in ", out_root)
