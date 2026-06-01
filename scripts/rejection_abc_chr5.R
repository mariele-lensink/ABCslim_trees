#!/usr/bin/env Rscript

library(abc)
library(data.table)
library(ggplot2)
library(scales)
library(patchwork)

setwd("/home/mlensink/slimsimulations/ABCslim_trees")

chrom_dir <- Sys.getenv("ABCREJ_CHROM", "chrom5")
if (!grepl("^chrom[0-9A-Za-z_]+$", chrom_dir)) {
  stop("ABCREJ_CHROM must look like 'chrom1'")
}
chrom_label <- sub("^chrom", "chr", chrom_dir)
uncapped_mean_filter <- Sys.getenv("ABC_UNCAPPED_MEAN_FILTER", "0") == "1"
filter_q <- as.numeric(Sys.getenv("ABC_FILTER_Q", "50"))
filter_cap <- as.numeric(Sys.getenv("ABC_FILTER_CAP", "1"))
filter_label <- if (uncapped_mean_filter) paste0("_uncapped_mean_Q", filter_q) else ""
analysis_label <- Sys.getenv("ABC_ANALYSIS_LABEL", "")
if (nzchar(analysis_label)) {
  filter_label <- paste0(filter_label, "_", analysis_label)
}
stat_set <- Sys.getenv("ABCREJ_STAT_SET", "all_variable")
valid_stat_sets <- c("base", "perbp_base", "all_variable")
if (!stat_set %in% valid_stat_sets) {
  stop("ABCREJ_STAT_SET must be one of: ", paste(valid_stat_sets, collapse = ", "))
}

target_stats_file <- Sys.getenv(
  "ABC_TARGET_STATS_FILE",
  file.path("output", "ABC_results", "observed", paste0("target_stats_", chrom_dir, ".csv"))
)
sim_stats_file <- Sys.getenv(
  "ABC_SIM_STATS_FILE",
  file.path("output", "simstats", paste0("sim_stats_", chrom_dir, ".csv"))
)
sim_params_file <- Sys.getenv("ABC_PRIORS_FILE", "priors/priors_2.6.25.csv")
out_root <- Sys.getenv("ABC_OUTPUT_DIR", file.path("output", "ABC_results"))
out_dir <- file.path(out_root, paste0("rejection_", chrom_label, filter_label))
figure_dir <- Sys.getenv("ABC_FIGURE_DIR", "figures")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

cat("Rejection ABC chromosome:", chrom_dir, "\n")
cat("Target stats file:", target_stats_file, "\n")
cat("Simulation stats file:", sim_stats_file, "\n")
cat("Priors file:", sim_params_file, "\n")
cat("Analysis label:", ifelse(nzchar(analysis_label), analysis_label, "(none)"), "\n")
cat("Rejection stat set:", stat_set, "\n")
cat("Uncapped mean filter:", uncapped_mean_filter, "\n")
if (uncapped_mean_filter) {
  cat("Keeping rows with abs(gdfe) *", filter_q, "<=", filter_cap,
      "and abs(idfe) *", filter_q, "<=", filter_cap, "\n")
}
target_stats <- fread(target_stats_file)
sim_stats <- fread(sim_stats_file)
sim_params <- fread(sim_params_file)

param_cols <- c("gmu", "imu", "gd", "id", "gdfe", "idfe")
tol <- 0.01
tols_to_compare <- c(0.005, 0.01, 0.02)

safe_ratio <- function(num, den) {
  out <- num / den
  out[!is.finite(out)] <- NA_real_
  out
}

add_contrasts <- function(dt) {
  dt[, snp_ratio := safe_ratio(g_snps, i_snps)]
  dt[, pi_ratio := safe_ratio(g_pi, i_pi)]
  dt[, theta_ratio := safe_ratio(g_theta_w, i_theta_w)]
  dt[, td_diff := g_td - i_td]
  dt
}

sim_stats <- add_contrasts(copy(sim_stats))
target_stats <- add_contrasts(copy(target_stats))

# Match prior rows to the simulations that have summary statistics.
par <- sim_params[match(sim_stats$ID, ID)]

if (anyNA(par$ID)) {
  stop("Some simulation IDs in sim_stats were not found in priors.")
}

if (uncapped_mean_filter) {
  filter_keep <- (abs(par$gdfe) * filter_q <= filter_cap) &
    (abs(par$idfe) * filter_q <= filter_cap)
  cat("Rows before uncapped-mean filter:", nrow(par), "\n")
  cat("Rows after uncapped-mean filter:", sum(filter_keep), "\n")
  sim_stats <- sim_stats[filter_keep]
  par <- par[filter_keep]
}

# Use either a minimal SNP-count/Tajima's D set, the statv2 per-bp SNP/Tajima's
# D set, or the expanded folded-SFS/downsampled stats. For the expanded set,
# exclude ratio-derived summaries before distance calculation.
base_stats <- c("g_snps", "i_snps", "g_td", "i_td")
perbp_base_stats <- c("g_snps_perbp", "i_snps_perbp", "g_td", "i_td")
ratio_stats <- c("snp_ratio", "pi_ratio", "theta_ratio")
sumstat_cols <- switch(
  stat_set,
  base = base_stats,
  perbp_base = perbp_base_stats,
  all_variable = setdiff(
    names(sim_stats),
    c("ID", "g_sites_total", "i_sites_total", ratio_stats)
  )
)
sumstat_cols <- intersect(sumstat_cols, names(sim_stats))
sumstat_cols <- sumstat_cols[
  vapply(sim_stats[, ..sumstat_cols], function(x) sd(x, na.rm = TRUE) > 0, logical(1))
]

keep <- complete.cases(sim_stats[, ..sumstat_cols]) &
  complete.cases(par[, ..param_cols])

sim_stats <- sim_stats[keep]
par <- par[keep]

sim_mat <- as.matrix(sim_stats[, ..sumstat_cols])
sim_scaled <- scale(sim_mat)

target_mat <- as.matrix(target_stats[, ..sumstat_cols])
target_scaled <- scale(
  target_mat,
  center = attr(sim_scaled, "scaled:center"),
  scale = attr(sim_scaled, "scaled:scale")
)

res <- abc(
  target = as.data.frame(target_scaled),
  param = as.data.frame(par[, ..param_cols]),
  sumstat = as.data.frame(sim_scaled),
  tol = tol,
  method = "rejection"
)

# Compute the same z-scaled Euclidean distance explicitly so the accepted IDs
# can be tracked and plotted alongside the full prior.
target_vec <- as.numeric(target_scaled[1, ])
par[, dist := sqrt(rowSums((sim_scaled - matrix(
  target_vec,
  nrow = nrow(sim_scaled),
  ncol = ncol(sim_scaled),
  byrow = TRUE
))^2))]

accepted_n <- ceiling(tol * nrow(par))
eps <- sort(par$dist)[accepted_n]
par[, type := "Prior"]
par[order(dist)[seq_len(accepted_n)], type := "Posterior"]
par[, type := factor(type, levels = c("Prior", "Posterior"))]

post <- copy(par[type == "Posterior", c("ID", param_cols, "dist"), with = FALSE])
post[, type := "Posterior"]

longpost <- melt(
  post,
  id.vars = c("ID", "dist", "type"),
  measure.vars = param_cols,
  variable.name = "parameter",
  value.name = "value"
)

posterior_summary <- longpost[, .(
  mean = mean(value, na.rm = TRUE),
  median = median(value, na.rm = TRUE),
  q025 = quantile(value, 0.025, na.rm = TRUE),
  q975 = quantile(value, 0.975, na.rm = TRUE)
), by = parameter]

tol_summary <- rbindlist(lapply(tols_to_compare, function(tol_i) {
  n_i <- ceiling(tol_i * nrow(par))
  ids_i <- par[order(dist)[seq_len(n_i)], ID]
  post_i <- melt(
    par[ID %in% ids_i, c("ID", param_cols), with = FALSE],
    id.vars = "ID",
    measure.vars = param_cols,
    variable.name = "parameter",
    value.name = "value"
  )
  post_i[, .(
    tol = tol_i,
    accepted_n = n_i,
    mean = mean(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    q025 = quantile(value, 0.025, na.rm = TRUE),
    q975 = quantile(value, 0.975, na.rm = TRUE)
  ), by = parameter]
}))

fwrite(
  posterior_summary,
  file.path(out_dir, paste0("posterior_summary_", chrom_label, filter_label, "_rejection_tol001.csv"))
)
fwrite(
  tol_summary,
  file.path(out_dir, paste0("posterior_summary_", chrom_label, filter_label, "_rejection_tol_sensitivity.csv"))
)
fwrite(
  post,
  file.path(out_dir, paste0("posterior_samples_", chrom_label, filter_label, "_rejection_tol001.csv"))
)

# Posterior distributions of mutation rates.
longsub <- longpost[parameter %in% c("gmu", "imu")]
vline_df <- longsub[, .(
  x = mean(value, na.rm = TRUE),
  label = label_scientific(digits = 2)(mean(value, na.rm = TRUE))
), by = parameter]

p_mrate <- ggplot(longsub, aes(x = value)) +
  geom_density(aes(fill = parameter), alpha = 0.4) +
  geom_vline(
    data = vline_df,
    aes(xintercept = x, color = parameter),
    linetype = "dashed",
    linewidth = 1
  ) +
  geom_text(
    data = vline_df,
    aes(x = x, y = 0, label = label, color = parameter),
    angle = 90,
    vjust = 1.5,
    hjust = -0.1,
    size = 4
  ) +
  scale_color_manual(values = c(gmu = "red", imu = "blue")) +
  scale_fill_manual(values = c(gmu = "pink", imu = "lightblue")) +
  facet_grid(parameter ~ ., scales = "fixed", space = "free") +
  scale_x_log10(labels = label_scientific()) +
  labs(
    title = "Posterior Distributions of Genic and Intergenic Mutation Rates",
    x = "Mutation rate per site per generation, log10",
    y = "Posterior density"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    strip.text.y = element_text(size = 12, face = "bold", angle = 0),
    legend.position = "none"
  )

ggsave(
  file.path(figure_dir, paste0("posterior_mutation_rates_", chrom_label, filter_label, "_rejection_tol001.png")),
  p_mrate,
  width = 7,
  height = 5
)

# Prior/posterior density comparison for all inferred parameters.
longdt <- melt(
  par[, c("ID", "type", param_cols), with = FALSE],
  id.vars = c("ID", "type"),
  measure.vars = param_cols,
  variable.name = "parameter",
  value.name = "value"
)

d_mu <- copy(longdt[parameter %in% c("gmu", "imu")])
d_mu[, facet := factor(parameter, levels = c("gmu", "imu"),
                       labels = c("Genic mu", "Intergenic mu"))]

p_mu <- ggplot(d_mu, aes(x = value, fill = type, color = type)) +
  geom_density(alpha = 0.40, adjust = 1) +
  scale_x_log10(labels = label_scientific(), name = "mu per site per generation") +
  scale_fill_manual(values = c(Prior = "blue", Posterior = "red")) +
  scale_color_manual(values = c(Prior = "blue4", Posterior = "firebrick")) +
  facet_grid(rows = vars(facet), scales = "fixed") +
  labs(title = "Mutation rates") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), strip.text = element_text(face = "bold"))

d_p <- copy(longdt[parameter %in% c("gd", "id")])
d_p[, facet := factor(parameter, levels = c("gd", "id"),
                      labels = c("Genic P(del)", "Intergenic P(del)"))]

p_prop <- ggplot(d_p, aes(x = value, fill = type, color = type)) +
  geom_density(alpha = 0.40, adjust = 1) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = label_percent(accuracy = 1),
    name = "Proportion deleterious"
  ) +
  scale_fill_manual(values = c(Prior = "blue", Posterior = "red")) +
  scale_color_manual(values = c(Prior = "blue4", Posterior = "firebrick")) +
  facet_grid(rows = vars(facet)) +
  labs(title = "Fraction of mutations that are deleterious") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )

d_s <- copy(longdt[parameter %in% c("gdfe", "idfe")])
d_s[, abs_s := pmax(abs(value), .Machine$double.eps)]
d_s[, facet := factor(parameter, levels = c("gdfe", "idfe"),
                      labels = c("Genic |s|", "Intergenic |s|"))]

p_s <- ggplot(d_s, aes(x = abs_s, fill = type, color = type)) +
  geom_density(alpha = 0.40, adjust = 1) +
  scale_x_log10(labels = label_scientific(), name = "|s| per generation") +
  scale_fill_manual(values = c(Prior = "blue", Posterior = "red")) +
  scale_color_manual(values = c(Prior = "blue4", Posterior = "firebrick")) +
  facet_grid(rows = vars(facet)) +
  labs(title = "Strength of deleterious effects") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )

final_density_plot <- ((p_mu | p_prop | p_s) +
  plot_layout(guides = "collect")) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(color = "grey50", fill = NA, linewidth = 0.5)
  )

ggsave(
  file.path(figure_dir, paste0("prior_posterior_density_", chrom_label, filter_label, "_rejection_tol001.png")),
  final_density_plot,
  width = 12,
  height = 8
)

p_dist <- ggplot(par, aes(x = dist)) +
  geom_histogram(bins = 50, fill = "grey70", color = "white") +
  geom_vline(xintercept = eps, color = "red", linetype = "dashed") +
  labs(
    title = "Distribution of Distances",
    subtitle = paste0("Red dashed = epsilon, top ", accepted_n, " simulations"),
    x = "Euclidean distance to observed summary stats, z-scaled",
    y = "Count"
  ) +
  theme_minimal()

ggsave(
  file.path(figure_dir, paste0("distance_histogram_", chrom_label, filter_label, "_rejection_tol001.png")),
  p_dist,
  width = 7,
  height = 5
)

ratio_dt <- rbindlist(list(
  par[, .(dist, type, metric = "gmu / imu", value = safe_ratio(gmu, imu))],
  par[, .(dist, type, metric = "gd / id", value = safe_ratio(gd, id))],
  par[, .(dist, type, metric = "|gdfe| / |idfe|", value = safe_ratio(abs(gdfe), abs(idfe)))]
))
ratio_dt <- ratio_dt[is.finite(dist) & is.finite(value) & dist > 0 & value > 0]
ratio_dt[, metric := factor(metric, levels = c("gmu / imu", "gd / id", "|gdfe| / |idfe|"))]

ratio_means <- ratio_dt[type == "Posterior", .(
  mean_value = mean(value, na.rm = TRUE)
), by = metric]

p_ratios <- ggplot(ratio_dt, aes(x = dist, y = value, color = type)) +
  geom_point(alpha = 0.25, stroke = 0, size = 1) +
  geom_vline(xintercept = eps, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 1, linetype = "dotted", color = "grey40") +
  geom_hline(
    data = ratio_means,
    aes(yintercept = mean_value),
    color = "firebrick",
    linetype = "solid",
    inherit.aes = FALSE
  ) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(values = c(Prior = "grey75", Posterior = "firebrick")) +
  facet_grid(metric ~ ., scales = "free_y") +
  labs(
    title = "Parameter ratios as a function of distance to observed stats",
    subtitle = "Vertical dashed line marks epsilon; dotted line marks ratio = 1; solid red line is posterior mean",
    x = "Distance to observed summary stats, z-scaled log10",
    y = "Parameter ratio, log10"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(
  file.path(figure_dir, paste0("ratio_distance_faceted_", chrom_label, filter_label, "_rejection_tol001.png")),
  p_ratios,
  width = 7,
  height = 10
)

# Posterior summary-statistic histograms for accepted simulations.
post_stats <- sim_stats[post[, .(ID)], on = "ID", nomatch = 0]
if (nrow(post_stats) > 0) {
  hist_stat_cols <- sumstat_cols
  hist_stat_cols <- hist_stat_cols[
    vapply(post_stats[, ..hist_stat_cols], is.numeric, logical(1)) &
      hist_stat_cols %in% names(target_stats)
  ]

  long_stats <- melt(
    post_stats[, c("ID", hist_stat_cols), with = FALSE],
    id.vars = "ID",
    variable.name = "statistic",
    value.name = "value"
  )
  observed_stats <- melt(
    target_stats[, ..hist_stat_cols],
    measure.vars = hist_stat_cols,
    variable.name = "statistic",
    value.name = "observed"
  )

  p_summary_stats <- ggplot(long_stats, aes(x = value)) +
    geom_histogram(bins = 40, fill = "grey70", color = "white", na.rm = TRUE) +
    geom_vline(
      data = observed_stats,
      aes(xintercept = observed),
      inherit.aes = FALSE,
      color = "#b2182b",
      linewidth = 0.7
    ) +
    facet_wrap(~ statistic, scales = "free", ncol = 4) +
    scale_x_continuous(labels = label_number()) +
    labs(
      title = "Accepted rejection posterior summary statistics",
      subtitle = "Red vertical line marks observed 1001 Genomes target value",
      x = "Accepted simulation statistic value",
      y = "Accepted simulations"
    ) +
    theme_bw(base_size = 12)

  ggsave(
    file.path(figure_dir, paste0("posterior_summary_stats_histograms_", chrom_label, filter_label, "_rejection_tol001.png")),
    p_summary_stats,
    width = 16,
    height = 18,
    dpi = 300
  )
} else {
  warning("No accepted posterior IDs matched sim_stats; skipping summary-statistic histograms.")
}

cat("Rejection ABC complete.\n")
cat("Rows in reference table:", nrow(par), "\n")
cat("Summary statistics used:", length(sumstat_cols), "\n")
cat("Accepted simulations at tol=", tol, ": ", accepted_n, "\n", sep = "")
cat("Results written to ", out_dir, "\n", sep = "")
