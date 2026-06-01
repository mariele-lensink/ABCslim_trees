#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default) {
  hit <- match(flag, args)
  if (is.na(hit) || hit == length(args)) {
    return(default)
  }
  args[[hit + 1]]
}

has_flag <- function(flag) {
  flag %in% args
}

data_dir <- get_arg("--data-dir", "/group/gmonroegrp2/mlensink/1001data")
metadata_csv <- get_arg("--metadata", file.path(data_dir, "1001_admixture.csv"))
accessions_csv <- get_arg("--accessions", file.path(data_dir, "1001genomes_accessions.csv"))
vcf_in <- get_arg("--vcf", file.path(data_dir, "1001genomes_snp-short-indel_only_ACGTN.vcf.gz"))
bcftools <- get_arg("--bcftools", "/home/mlensink/miniconda3/envs/abc_trees_vcf/bin/bcftools")
run_bcftools <- has_flag("--run-bcftools")

meta <- read.csv(metadata_csv, stringsAsFactors = FALSE, check.names = FALSE)
required_cols <- c("id", "country")
missing_cols <- setdiff(required_cols, names(meta))
if (length(missing_cols) > 0) {
  stop("Metadata is missing required columns: ", paste(missing_cols, collapse = ", "))
}

if (file.exists(accessions_csv)) {
  accessions <- read.csv(accessions_csv, stringsAsFactors = FALSE, header = FALSE)
  if (ncol(accessions) >= 11) {
    accession_groups <- data.frame(
      id = as.character(accessions[[1]]),
      admixture_group = accessions[[11]],
      stringsAsFactors = FALSE
    )
    meta <- merge(meta, accession_groups, by = "id", all.x = TRUE, sort = FALSE)
  }
}

if (!("admixture_group" %in% names(meta))) {
  meta$admixture_group <- NA_character_
}

vcf_header <- NULL
con <- gzfile(vcf_in, open = "rt")
on.exit(close(con), add = TRUE)
repeat {
  line <- readLines(con, n = 1)
  if (length(line) == 0) {
    break
  }
  if (startsWith(line, "#CHROM")) {
    vcf_header <- strsplit(line, "\t", fixed = TRUE)[[1]]
    break
  }
}

if (is.null(vcf_header) || length(vcf_header) < 10) {
  stop("Could not find VCF #CHROM header with sample IDs in ", vcf_in)
}

vcf_samples <- vcf_header[10:length(vcf_header)]
sample_set <- unique(vcf_samples)
meta$id <- as.character(meta$id)
meta_vcf <- meta[meta$id %in% sample_set, , drop = FALSE]
meta_vcf <- meta_vcf[match(sample_set, meta_vcf$id), , drop = FALSE]
meta_vcf <- meta_vcf[!is.na(meta_vcf$id), , drop = FALSE]

central_europe_clean_countries <- c("GER", "AUT", "SUI", "CZE", "BEL", "NED")
central_europe_clean_groups <- c("central_europe", "germany")

global_outlier_countries <- c(
  # Iberian relicts
  "ESP", "POR",
  # Other highly diverged relicts described in the 1,135 genomes paper
  "CPV", "LBN",
  # Scandinavian populations
  "SWE", "NOR", "FIN", "DEN", "Finland",
  # Asian populations
  "CHN", "JPN", "Japan", "KAZ", "KGZ", "TJK", "UZB", "AFG", "IND",
  # Caucasus-associated populations
  "GEO", "ARM", "AZE"
)

global_outlier_groups <- c(
  "relict",
  "spain",
  "north_sweden",
  "south_sweden",
  "asia",
  "italy_balkan_caucasus"
)

subsets <- list(
  central_europe_clean = meta_vcf[
    meta_vcf$country %in% central_europe_clean_countries &
      (is.na(meta_vcf$admixture_group) | meta_vcf$admixture_group %in% central_europe_clean_groups),
    ,
    drop = FALSE
  ],
  global_no_outliers = meta_vcf[
    !(meta_vcf$country %in% global_outlier_countries) &
      (is.na(meta_vcf$admixture_group) | !(meta_vcf$admixture_group %in% global_outlier_groups)),
    ,
    drop = FALSE
  ]
)

subset_summary <- data.frame(
  subset = character(),
  n_samples = integer(),
  countries = character(),
  admixture_groups = character(),
  keep_file = character(),
  vcf_file = character(),
  stringsAsFactors = FALSE
)

for (subset_name in names(subsets)) {
  subset_meta <- subsets[[subset_name]]
  subset_dir <- file.path(data_dir, subset_name)
  dir.create(subset_dir, recursive = TRUE, showWarnings = FALSE)

  keep_file <- file.path(subset_dir, paste0(subset_name, ".keep"))
  metadata_out <- file.path(subset_dir, paste0(subset_name, "_metadata.csv"))
  vcf_out <- file.path(subset_dir, paste0(subset_name, ".vcf.gz"))

  writeLines(subset_meta$id, keep_file)
  write.csv(subset_meta, metadata_out, row.names = FALSE, quote = TRUE)

  country_counts <- sort(table(subset_meta$country), decreasing = TRUE)
  country_summary <- paste(paste0(names(country_counts), ":", as.integer(country_counts)), collapse = ";")
  group_counts <- sort(table(subset_meta$admixture_group), decreasing = TRUE)
  group_summary <- paste(paste0(names(group_counts), ":", as.integer(group_counts)), collapse = ";")

  subset_summary <- rbind(
    subset_summary,
    data.frame(
      subset = subset_name,
      n_samples = nrow(subset_meta),
      countries = country_summary,
      admixture_groups = group_summary,
      keep_file = keep_file,
      vcf_file = vcf_out,
      stringsAsFactors = FALSE
    )
  )

  if (run_bcftools) {
    message("Subsetting ", subset_name, " (", nrow(subset_meta), " samples)")
    status <- system2(
      bcftools,
      args = c("view", "-S", keep_file, "-Oz", "-o", vcf_out, vcf_in)
    )
    if (status != 0) {
      stop("bcftools view failed for ", subset_name)
    }
    status <- system2(bcftools, args = c("index", "-t", "--force", vcf_out))
    if (status != 0) {
      stop("bcftools index failed for ", subset_name)
    }
  }
}

summary_out <- file.path(data_dir, "1001_subset_summary.csv")
write.csv(subset_summary, summary_out, row.names = FALSE, quote = TRUE)

print(subset_summary[, c("subset", "n_samples", "keep_file", "vcf_file")], row.names = FALSE)
if (!run_bcftools) {
  message("Keep lists written. Re-run with --run-bcftools to write subset VCFs.")
}
