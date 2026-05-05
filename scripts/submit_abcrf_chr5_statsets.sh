#!/bin/bash
set -euo pipefail

cd /home/mlensink/slimsimulations/ABCslim_trees

chrom="${ABCRF_CHROM:-chrom5}"
chrom_label="${chrom/chrom/chr}"
uncapped_filter="${ABC_UNCAPPED_MEAN_FILTER:-0}"
filter_q="${ABC_FILTER_Q:-50}"
filter_cap="${ABC_FILTER_CAP:-1}"
include_all_variable="${ABCRF_INCLUDE_ALL_VARIABLE:-1}"

for stat_set in base core core_contrast core_sfs all_variable; do
  if [[ "$stat_set" == "all_variable" && "$include_all_variable" != "1" ]]; then
    continue
  fi

  filter_tag=""
  if [[ "$uncapped_filter" == "1" ]]; then
    filter_tag="_uncapped_Q${filter_q}"
  fi

  sbatch \
    --job-name="abcrf_${chrom_label}${filter_tag}_${stat_set}" \
    --export=ALL,ABCRF_CHROM="${chrom}",ABCRF_STAT_SET="${stat_set}",ABCRF_NCORES="${ABCRF_NCORES:-4}",ABC_UNCAPPED_MEAN_FILTER="${uncapped_filter}",ABC_FILTER_Q="${filter_q}",ABC_FILTER_CAP="${filter_cap}" \
    scripts/abcrf_chr5.sbatch
done
