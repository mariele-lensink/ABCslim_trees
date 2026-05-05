args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) == 0) {
  stop("Unable to determine script path")
}

script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
python_script <- file.path(script_dir, "calc_sim_stats.py")
python_bin <- Sys.getenv(
  "PYTHON_BIN",
  unset = "/home/mlensink/miniconda3/envs/abc_trees_vcf/bin/python"
)

if (!nzchar(python_bin) || !file.exists(python_bin)) {
  stop("Python executable was not found; set PYTHON_BIN to a valid interpreter")
}

if (!file.exists(python_script)) {
  stop(sprintf("Python stats script missing: %s", python_script))
}

status <- system2(python_bin, python_script)
if (!identical(status, 0L)) {
  stop(sprintf("calc_sim_stats.py exited with status %d", status))
}
