##This script populates a table of parameter values for however many simulations you plan to run
library(data.table)
library(parallel)

# Set a base seed for reproducibility
set.seed(as.integer(Sys.time()))

# Generate 100,000 unique IDs and random parameters
num_simulations <- 100000
selection_min_abs <- 1e-3
# Old full prior range: 1e-3 to 1e-1. For Q25, 1e-1 scales to -2.5
# and many mutation effects hit the SLiM -1 selection-coefficient cap.
selection_max_abs <- 2e-2
out_file <- "/home/mlensink/slimsimulations/ABCslim_trees/priors/priors_Q25_smax2e-2_2026-05-05.csv"

log_uniform <- function(n, min, max) 10^runif(n, log10(min), log10(max))
log_uniform_neg <- function(n, min_abs, max_abs) -10^runif(n, log10(min_abs), log10(max_abs))

params <- data.table(
  ID = 1:num_simulations,
  #unscaled
  gmu = log_uniform(num_simulations, 1e-10, 1e-7),
  imu = log_uniform(num_simulations, 1e-10, 1e-7),
  #unscaled
  gd = runif(num_simulations,0.1,0.75),
  id = runif(num_simulations,0.1,0.75),
  #unscaled
  gdfe = log_uniform_neg(num_simulations, selection_min_abs, selection_max_abs),

  idfe = log_uniform_neg(num_simulations, selection_min_abs, selection_max_abs)
)

# Write to a CSV file
fwrite(params, out_file, sep = ",", col.names = TRUE)
