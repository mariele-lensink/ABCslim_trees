#!/bin/bash
#SBATCH -D /home/mlensink/slimsimulations/ABCslim_trees/
#SBATCH -o /home/mlensink/slimsimulations/ABCslim_trees/logs/subset_simstats/%x_%A_%a.out
#SBATCH -e /home/mlensink/slimsimulations/ABCslim_trees/logs/subset_simstats/%x_%A_%a.err
#SBATCH -J subset_simstats
#SBATCH -t 1-00:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --account=gmonroegrp
#SBATCH --partition=bmh
#SBATCH --array=0-99%25
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=melensink@ucdavis.edu

set -euo pipefail

source /home/mlensink/miniconda3/etc/profile.d/conda.sh
conda activate /home/mlensink/miniconda3/envs/abc_trees_vcf
echo "Python: $(which python)"

: "${TARGET_LABEL:?Set TARGET_LABEL to central_europe_clean or global_no_outliers}"
export CHROM_DIR="${CHROM_DIR:-chrom1}"
export MIN_ID="${MIN_ID:-1}"
export MAX_ID="${MAX_ID:-100000}"
export CHUNK_SIZE="${CHUNK_SIZE:-1000}"
export BASE_OUT="${BASE_OUT:-/home/mlensink/slimsimulations/ABCslim_trees/output}"
export DERIVED_DIR="${DERIVED_DIR:-${BASE_OUT}/constant/${CHROM_DIR}/trees_derived}"
export SIM_STATS_CORES="${SLURM_CPUS_PER_TASK}"
export SAMPLE_SEED="${SAMPLE_SEED:-1729}"
export FOLD_SFS="${FOLD_SFS:-1}"
export SNP_ONLY="${SNP_ONLY:-0}"

if [[ -z "${SAMPLE_INDIVIDUALS:-}" ]]; then
  case "${TARGET_LABEL}" in
    central_europe_clean) export SAMPLE_INDIVIDUALS=187 ;;
    global_no_outliers) export SAMPLE_INDIVIDUALS=635 ;;
    *) echo "Unknown TARGET_LABEL=${TARGET_LABEL}; set SAMPLE_INDIVIDUALS explicitly." >&2; exit 1 ;;
  esac
else
  export SAMPLE_INDIVIDUALS
fi

chunk_start=$((MIN_ID + SLURM_ARRAY_TASK_ID * CHUNK_SIZE))
chunk_end=$((chunk_start + CHUNK_SIZE - 1))
if (( chunk_start > MAX_ID )); then
  echo "Chunk start ${chunk_start} is above MAX_ID ${MAX_ID}; nothing to do."
  exit 0
fi
if (( chunk_end > MAX_ID )); then
  chunk_end="${MAX_ID}"
fi

export START_ID="${chunk_start}"
export END_ID="${chunk_end}"

chunk_dir="${BASE_OUT}/ABC_results/${TARGET_LABEL}/simstats/chunks/${CHROM_DIR}"
export SIM_STATS_OUT="${SIM_STATS_OUT:-${chunk_dir}/sim_stats_${CHROM_DIR}_${START_ID}_${END_ID}.csv}"

mkdir -p "$(dirname "${SIM_STATS_OUT}")" /home/mlensink/slimsimulations/ABCslim_trees/logs/subset_simstats

echo "Running subset-matched calc_sim_stats.py for IDs ${START_ID}-${END_ID}"
echo "TARGET_LABEL=${TARGET_LABEL}"
echo "CHROM_DIR=${CHROM_DIR}"
echo "DERIVED_DIR=${DERIVED_DIR}"
echo "SIM_STATS_CORES=${SIM_STATS_CORES}"
echo "SAMPLE_INDIVIDUALS=${SAMPLE_INDIVIDUALS}"
echo "SAMPLE_SEED=${SAMPLE_SEED}"
echo "FOLD_SFS=${FOLD_SFS}"
echo "SNP_ONLY=${SNP_ONLY}"
echo "Writing output to ${SIM_STATS_OUT}"

python ./scripts/calc_sim_stats.py
