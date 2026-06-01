#!/usr/bin/env bash
#SBATCH -D /home/mlensink/ABCslim_trees/
#SBATCH -o /home/mlensink/ABCslim_trees/logs/observed_1001_segments_v2/%x_%A_%a.out
#SBATCH -e /home/mlensink/ABCslim_trees/logs/observed_1001_segments_v2/%x_%A_%a.err
#SBATCH -J obs1001segv2
#SBATCH -t 1-00:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=1
#SBATCH --account=gmonroegrp
#SBATCH --partition=bml
#SBATCH --array=0-4
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=melensink@ucdavis.edu

set -euo pipefail

source /home/mlensink/miniconda3/etc/profile.d/conda.sh
conda activate "${CONDA_ENV:-/group/gmonroegrp2/mlensink/conda_envs/abc_trees}"

if [[ -z "${VCFTOOLS_BIN:-}" ]]; then
  if command -v vcftools >/dev/null 2>&1; then
    export VCFTOOLS_BIN="$(command -v vcftools)"
  elif command -v module >/dev/null 2>&1; then
    module load vcftools/0.1.16
    export VCFTOOLS_BIN="$(command -v vcftools)"
  else
    echo "Could not find vcftools. Set VCFTOOLS_BIN=/path/to/vcftools before submitting." >&2
    exit 1
  fi
fi

PROJECT_DIR="${PROJECT_DIR:-/home/mlensink/ABCslim_trees}"
SEGMENT_DIR="${SEGMENT_DIR:-/group/gmonroegrp2/mlensink/ABC_data/real_data_segments/vcf_subsets_all_sites}"
OUT_DIR="${OUT_DIR:-/group/gmonroegrp2/mlensink/ABC_data/real_data_segments/simstats_v2}"
TMP_ROOT="${TMP_ROOT:-${SLURM_TMPDIR:-/tmp}}"

cd "${PROJECT_DIR}"
mkdir -p "${OUT_DIR}" /home/mlensink/ABCslim_trees/logs/observed_1001_segments_v2

segments=(
  "chr1_10000000_10284120 chrom1 1 10000000 genomeinfo/gene100_chr1_10000000.csv genomeinfo/intergene100_chr1_10000000.csv"
  "chr2_18000000_18280622 chrom2 2 18000000 genomeinfo/gene100_chr2_18000000.csv genomeinfo/intergene100_chr2_18000000.csv"
  "chr3_22000000_22349526 chrom3 3 22000000 genomeinfo/gene100_chr3_22000000.csv genomeinfo/intergene100_chr3_22000000.csv"
  "chr4_2000000_2660337 chrom4 4 2000000 genomeinfo/gene100_chr4_2000000.csv genomeinfo/intergene100_chr4_2000000.csv"
  "chr5_1_377208 chrom5 5 1 genomeinfo/gene100_chr5_1.csv genomeinfo/intergene100_chr5_1.csv"
)

task_id="${SLURM_ARRAY_TASK_ID:-0}"
if (( task_id < 0 || task_id >= ${#segments[@]} )); then
  echo "Invalid SLURM_ARRAY_TASK_ID=${task_id}; expected 0-4" >&2
  exit 1
fi

read -r SEGMENT CHROM_DIR VCF_CHROM SEGMENT_START GENE_CSV INTERGENE_CSV <<< "${segments[$task_id]}"

export CHROM_DIR
export VCF_CHROM
export VCF_POS_OFFSET="$((SEGMENT_START - 1))"
export GENE_CSV
export INTERGENE_CSV
export VCF_FILE="${VCF_FILE:-${SEGMENT_DIR}/${SEGMENT}.vcf.gz}"
export VCF_FILE_ID="$((task_id + 1))"
export SIM_STATS_OUT="${SIM_STATS_OUT:-${OUT_DIR}/${SEGMENT}.sim_stats_v2.csv}"
export SIM_STATS_CORES="${SLURM_CPUS_PER_TASK:-1}"
export SAMPLE_INDIVIDUALS="${SAMPLE_INDIVIDUALS:-1135}"
export SAMPLE_SEED="${SAMPLE_SEED:-1729}"
export FOLD_SFS="${FOLD_SFS:-1}"
export SNP_ONLY="${SNP_ONLY:-1}"
export REQUIRE_COMPLETE_VCF_SAMPLES="${REQUIRE_COMPLETE_VCF_SAMPLES:-0}"
export MIN_VCF_CALL_FRACTION="${MIN_VCF_CALL_FRACTION:-0.8}"
export ALL_SITES_VCF="${ALL_SITES_VCF:-1}"
export SIM_STATS_TMP_DIR="${SIM_STATS_TMP_DIR:-${TMP_ROOT}/calc_sim_stats_v2_${SLURM_JOB_ID:-manual}_${SEGMENT}}"

echo "Running v2 observed segment stats"
echo "SEGMENT=${SEGMENT}"
echo "VCF_FILE=${VCF_FILE}"
echo "CHROM_DIR=${CHROM_DIR}"
echo "VCF_CHROM=${VCF_CHROM}"
echo "SEGMENT_START=${SEGMENT_START}"
echo "VCF_POS_OFFSET=${VCF_POS_OFFSET}"
echo "GENE_CSV=${GENE_CSV}"
echo "INTERGENE_CSV=${INTERGENE_CSV}"
echo "SIM_STATS_OUT=${SIM_STATS_OUT}"
echo "ALL_SITES_VCF=${ALL_SITES_VCF}"
echo "MIN_VCF_CALL_FRACTION=${MIN_VCF_CALL_FRACTION}"
echo "VCFTOOLS_BIN=${VCFTOOLS_BIN}"

python ./scripts/calc_sim_stats_v2.py
