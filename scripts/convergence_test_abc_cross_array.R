#!/usr/bin/env Rscript

library(data.table)
library(abcrf)

repo_dir <- Sys.getenv("ABC_REPO_DIR", getwd())
setwd(repo_dir)

sim_stats_file <- Sys.getenv(
  "ABC_CONV_SIM_STATS_FILE",
  "../../../group/gmonroegrp2/mlensink/ABC_data/chrom5/sim_stats_v2_chrom5_merged.csv"
)
target_sim_stats_file <- Sys.getenv("ABC_CONV_TARGET_SIM_STATS_FILE", sim_stats_file)
priors_file <- Sys.getenv("ABC_CONV_PRIORS_FILE", "priors/priors_2.6.25.csv")
out_root <- Sys.getenv(
  "ABC_CONV_OUTPUT_DIR",
  "../../../group/gmonroegrp2/mlensink/ABC_data/convergence_test"
)
n_reps <- as.integer(Sys.getenv("ABC_CONV_N_REPS", "10"))
seed <- as.integer(Sys.getenv("ABC_CONV_SEED", "20260526"))
tol <- as.numeric(Sys.getenv("ABC_CONV_REJECTION_TOL", "0.01"))
ntree <- as.integer(Sys.getenv("ABC_CONV_RF_NTREE", "500"))
ncores <- min(as.integer(Sys.getenv("ABCRF_NCORES", Sys.getenv("SLURM_CPUS_PER_TASK", "4"))), 4)
rep_index <- as.integer(Sys.getenv("ABC_CONV_REP_INDEX", Sys.getenv("SLURM_ARRAY_TASK_ID", "1")))

param_cols <- c("gmu", "imu", "gd", "id", "gdfe", "idfe")
core4_stats <- c("g_snps_perbase", "i_snps_perbase", "g_td", "i_td")

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "observed"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "rejection_core4"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "abcrf_core4"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "abcrf_all_variable"), recursive = TRUE, showWarnings = FALSE)

cat("Convergence-test reference sim stats:", sim_stats_file, "\n")
cat("Convergence-test target sim stats:", target_sim_stats_file, "\n")
cat("Convergence-test priors:", priors_file, "\n")
cat("Convergence-test output:", out_root, "\n")
cat("Replicates:", n_reps, "\n")
cat("Seed:", seed, "\n")
cat("Rejection tolerance:", tol, "\n")
cat("ABCRF ntree:", ntree, "\n")
cat("ABCRF cores:", ncores, "\n")
cat("Replicate index:", rep_index, "\n")

safe_ratio <- function(num, den) {
  out <- num / den
  out[!is.finite(out)] <- NA_real_
  out
}

weighted_quantile <- function(x, w, probs) {
  ok <- is.finite(x) & is.finite(w) & w >= 0
  x <- x[ok]
  w <- w[ok]
  if (length(x) == 0 || sum(w) <= 0) {
    return(rep(NA_real_, length(probs)))
  }
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cw <- cumsum(w) / sum(w)
  vapply(probs, function(p) x[which(cw >= p)[1]], numeric(1))
}

summarize_unweighted <- function(dt, value_col, true_dt, method, stat_set, rep_i, obs_id) {
  rbindlist(lapply(param_cols, function(param) {
    value <- dt[[param]]
    true_value <- true_dt[[param]]
    data.table(
      replicate = rep_i,
      observed_ID = obs_id,
      method = method,
      stat_set = stat_set,
      parameter = param,
      true_value = true_value,
      mean = mean(value, na.rm = TRUE),
      median = median(value, na.rm = TRUE),
      q025 = quantile(value, 0.025, na.rm = TRUE),
      q975 = quantile(value, 0.975, na.rm = TRUE),
      mean_error = mean(value, na.rm = TRUE) - true_value,
      median_error = median(value, na.rm = TRUE) - true_value
    )
  }))
}

summarize_weighted <- function(value, weight, true_value, method, stat_set, rep_i, obs_id, param, nmae) {
  qs <- weighted_quantile(value, weight, c(0.025, 0.5, 0.975))
  wmean <- if (sum(weight, na.rm = TRUE) > 0) {
    weighted.mean(value, weight, na.rm = TRUE)
  } else {
    NA_real_
  }
  data.table(
    replicate = rep_i,
    observed_ID = obs_id,
    method = method,
    stat_set = stat_set,
    parameter = param,
    true_value = true_value,
    mean = wmean,
    median = qs[2],
    q025 = qs[1],
    q975 = qs[3],
    mean_error = wmean - true_value,
    median_error = qs[2] - true_value,
    post_NMAE_mean = nmae
  )
}

add_v2_derived_stats <- function(dt) {
  if (all(c("g_snps_copy", "g_callable_sites") %in% names(dt))) {
    dt[, g_snps_perbase := safe_ratio(g_snps_copy, g_callable_sites)]
  }
  if (all(c("i_snps_copy", "i_callable_sites") %in% names(dt))) {
    dt[, i_snps_perbase := safe_ratio(i_snps_copy, i_callable_sites)]
  }
  if (all(c("g_snps_perbase", "i_snps_perbase") %in% names(dt))) {
    dt[, snp_ratio := safe_ratio(g_snps_perbase, i_snps_perbase)]
  }
  if (all(c("g_pi", "i_pi") %in% names(dt))) {
    dt[, pi_ratio := safe_ratio(g_pi, i_pi)]
  }
  if (all(c("g_theta_w", "i_theta_w") %in% names(dt))) {
    dt[, theta_ratio := safe_ratio(g_theta_w, i_theta_w)]
  }
  if (all(c("g_td", "i_td") %in% names(dt))) {
    dt[, td_diff := g_td - i_td]
  }
  dt
}

select_stats <- function(ref, target, stat_set) {
  stat_cols <- switch(
    stat_set,
    core4 = core4_stats,
    all_variable = setdiff(names(ref), c("ID", param_cols))
  )
  stat_cols <- intersect(stat_cols, names(ref))
  stat_cols <- stat_cols[vapply(ref[, ..stat_cols], is.numeric, logical(1))]
  stat_cols <- stat_cols[vapply(ref[, ..stat_cols], function(x) all(is.finite(x)), logical(1))]
  stat_cols <- stat_cols[vapply(target[, ..stat_cols], function(x) all(is.finite(x)), logical(1))]
  stat_cols <- stat_cols[vapply(ref[, ..stat_cols], function(x) sd(x, na.rm = TRUE) > 0, logical(1))]
  if (length(stat_cols) == 0) {
    stop("No usable statistics for stat set: ", stat_set)
  }
  stat_cols
}

run_rejection <- function(ref, target, true_dt, stat_cols, rep_i, obs_id) {
  ref_mat <- as.matrix(ref[, ..stat_cols])
  scaled_ref <- scale(ref_mat)
  scaled_target <- scale(
    as.matrix(target[, ..stat_cols]),
    center = attr(scaled_ref, "scaled:center"),
    scale = attr(scaled_ref, "scaled:scale")
  )
  target_vec <- as.numeric(scaled_target[1, ])
  dist <- sqrt(rowSums((scaled_ref - matrix(
    target_vec,
    nrow = nrow(scaled_ref),
    ncol = ncol(scaled_ref),
    byrow = TRUE
  ))^2))
  accepted_n <- max(1, ceiling(tol * nrow(ref)))
  accepted <- copy(ref[order(dist)[seq_len(accepted_n)], c("ID", param_cols), with = FALSE])
  accepted[, dist := sort(dist)[seq_len(accepted_n)]]

  sample_path <- file.path(
    out_root,
    "rejection_core4",
    sprintf("posterior_samples_replicate_%02d_observed_%s.csv", rep_i, obs_id)
  )
  fwrite(accepted, sample_path)

  summary <- summarize_unweighted(
    accepted,
    value_col = "value",
    true_dt = true_dt,
    method = "rejection",
    stat_set = "core4",
    rep_i = rep_i,
    obs_id = obs_id
  )
  summary[, accepted_n := accepted_n]
  summary[, reference_n := nrow(ref)]
  summary
}

run_abcrf <- function(ref, target, true_dt, stat_cols, stat_set, rep_i, obs_id) {
  obs_sumstats <- as.data.frame(target[, ..stat_cols])
  summaries <- list()
  importances <- list()

  for (param in param_cols) {
    cat("Replicate", rep_i, "ABCRF", stat_set, "parameter", param, "\n")
    model_data <- as.data.frame(ref[, c(param, stat_cols), with = FALSE])
    model <- regAbcrf(
      as.formula(paste(param, "~ .")),
      data = model_data,
      ntree = ntree,
      paral = TRUE,
      ncores = ncores
    )
    posterior <- predict(
      model,
      obs = obs_sumstats,
      training = model_data,
      paral = TRUE,
      ncores = ncores,
      rf.weights = TRUE
    )

    value <- model_data[[param]]
    weight <- as.numeric(drop(posterior$weights))
    summaries[[param]] <- summarize_weighted(
      value = value,
      weight = weight,
      true_value = true_dt[[param]],
      method = "abcrf",
      stat_set = stat_set,
      rep_i = rep_i,
      obs_id = obs_id,
      param = param,
      nmae = posterior$post.NMAE.mean
    )

    imp <- model$model.rf$variable.importance
    if (!is.null(imp) && length(imp) > 0) {
      importances[[param]] <- data.table(
        replicate = rep_i,
        observed_ID = obs_id,
        stat_set = stat_set,
        parameter = param,
        statistic = names(imp),
        importance = as.numeric(imp)
      )
    }
    rm(model, posterior, model_data)
    gc()
  }

  if (length(importances) > 0) {
    fwrite(
      rbindlist(importances),
      file.path(
        out_root,
        paste0("abcrf_", stat_set),
        sprintf("variable_importance_replicate_%02d_observed_%s.csv", rep_i, obs_id)
      )
    )
  }

  rbindlist(summaries)
}

accuracy_summary <- function(summary_dt) {
  summary_dt[, .(
    n_replicates = uniqueN(replicate),
    mean_abs_median_error = mean(abs(median_error), na.rm = TRUE),
    mean_abs_mean_error = mean(abs(mean_error), na.rm = TRUE),
    coverage_95 = mean(true_value >= q025 & true_value <= q975, na.rm = TRUE)
  ), by = .(method, stat_set, parameter)]
}

write_independent_outputs <- function(summary_by_run, stat_manifest_by_run) {
  for (run_name in names(summary_by_run)) {
    if (length(summary_by_run[[run_name]]) == 0) {
      next
    }

    run_summary <- rbindlist(summary_by_run[[run_name]], fill = TRUE)
    fwrite(
      run_summary,
      file.path(out_root, run_name, "posterior_summary_by_replicate.csv")
    )
    fwrite(
      accuracy_summary(run_summary),
      file.path(out_root, run_name, "convergence_accuracy_summary.csv")
    )

    if (length(stat_manifest_by_run[[run_name]]) > 0) {
      fwrite(
        rbindlist(stat_manifest_by_run[[run_name]]),
        file.path(out_root, run_name, "summary_statistics_used_by_replicate.csv")
      )
    }
  }
}

ref_sim_stats <- add_v2_derived_stats(fread(sim_stats_file))
target_sim_stats <- add_v2_derived_stats(fread(target_sim_stats_file))
priors <- fread(priors_file)

missing_params <- setdiff(c("ID", param_cols), names(priors))
if (length(missing_params) > 0) {
  stop("Missing prior columns: ", paste(missing_params, collapse = ", "))
}

common_ids <- Reduce(intersect, list(ref_sim_stats$ID, target_sim_stats$ID, priors$ID))
if (length(common_ids) < n_reps) {
  stop("Need at least ", n_reps, " shared IDs between sim stats and priors.")
}

set.seed(seed)
observed_ids <- sample(common_ids, n_reps)
if (is.na(rep_index) || rep_index < 1 || rep_index > length(observed_ids)) {
  stop("ABC_CONV_REP_INDEX must be between 1 and ", length(observed_ids))
}
observed_table <- merge(
  data.table(replicate = seq_along(observed_ids), observed_ID = observed_ids),
  priors[, c("ID", param_cols), with = FALSE],
  by.x = "observed_ID",
  by.y = "ID",
  sort = FALSE
)
if (rep_index == 1) {
  fwrite(observed_table, file.path(out_root, "observed", "sampled_observed_ids_and_true_parameters.csv"))
}

ref_table <- merge(priors, ref_sim_stats, by = "ID")
ref_table <- add_v2_derived_stats(as.data.table(ref_table))

all_summaries <- list()
stat_manifest <- list()
summary_by_run <- list(
  rejection_core4 = list(),
  abcrf_core4 = list(),
  abcrf_all_variable = list()
)
stat_manifest_by_run <- list(
  rejection_core4 = list(),
  abcrf_core4 = list(),
  abcrf_all_variable = list()
)

for (rep_i in rep_index) {
  obs_id <- observed_ids[[rep_i]]
  cat("Starting replicate", rep_i, "observed ID", obs_id, "\n")
  target <- copy(target_sim_stats[ID == obs_id, setdiff(names(target_sim_stats), "ID"), with = FALSE])
  target[, ID := obs_id]
  true_dt <- priors[ID == obs_id]
  ref <- copy(ref_table[ID != obs_id])

  stat_cols <- select_stats(ref, target, "core4")
  run_name <- "rejection_core4"
  manifest_dt <- data.table(
    replicate = rep_i,
    observed_ID = obs_id,
    method = "rejection",
    stat_set = "core4",
    statistic = stat_cols
  )
  stat_manifest[[length(stat_manifest) + 1]] <- manifest_dt
  stat_manifest_by_run[[run_name]][[length(stat_manifest_by_run[[run_name]]) + 1]] <- manifest_dt
  summary_dt <- run_rejection(
    ref,
    target,
    true_dt,
    stat_cols,
    rep_i,
    obs_id
  )
  all_summaries[[length(all_summaries) + 1]] <- summary_dt
  summary_by_run[[run_name]][[length(summary_by_run[[run_name]]) + 1]] <- summary_dt

  for (stat_set in c("core4", "all_variable")) {
    stat_cols <- select_stats(ref, target, stat_set)
    run_name <- paste0("abcrf_", stat_set)
    manifest_dt <- data.table(
      replicate = rep_i,
      observed_ID = obs_id,
      method = "abcrf",
      stat_set = stat_set,
      statistic = stat_cols
    )
    stat_manifest[[length(stat_manifest) + 1]] <- manifest_dt
    stat_manifest_by_run[[run_name]][[length(stat_manifest_by_run[[run_name]]) + 1]] <- manifest_dt
    summary_dt <- run_abcrf(
      ref,
      target,
      true_dt,
      stat_cols,
      stat_set,
      rep_i,
      obs_id
    )
    all_summaries[[length(all_summaries) + 1]] <- summary_dt
    summary_by_run[[run_name]][[length(summary_by_run[[run_name]]) + 1]] <- summary_dt
  }

  partial_dir <- file.path(out_root, "_replicates")
  dir.create(partial_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(
    rbindlist(all_summaries, fill = TRUE),
    file.path(partial_dir, sprintf("posterior_summary_replicate_%02d.csv", rep_i))
  )
  fwrite(
    rbindlist(stat_manifest),
    file.path(partial_dir, sprintf("summary_statistics_replicate_%02d.csv", rep_i))
  )
}

cat("Convergence replicate", rep_index, "complete. Results saved in", file.path(out_root, "_replicates"), "\n")
