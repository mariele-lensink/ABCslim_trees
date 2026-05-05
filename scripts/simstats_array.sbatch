#!/bin/bash
#SBATCH -D /home/mlensink/slimsimulations/ABCslim_trees/
#SBATCH -o /home/mlensink/slimsimulations/ABCslim_trees/logs/simstats_chunks/%x_%A_%a.out
#SBATCH -e /home/mlensink/slimsimulations/ABCslim_trees/logs/simstats_chunks/%x_%A_%a.err
#SBATCH -J simstats
#SBATCH -t 1-00:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --account=gmonroegrp
#SBATCH --partition=bml
#SBATCH --array=0-99
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=melensink@ucdavis.edu

set -euo pipefail

source /home/mlensink/miniconda3/etc/profile.d/conda.sh
conda activate /home/mlensink/miniconda3/envs/abc_trees_vcf
echo "Python: $(which python)"

export CHROM_DIR="${CHROM_DIR:-chrom1}"
export MIN_ID="${MIN_ID:-1}"
export MAX_ID="${MAX_ID:-100000}"
export CHUNK_SIZE="${CHUNK_SIZE:-1000}"
export BASE_OUT="${BASE_OUT:-/home/mlensink/slimsimulations/ABCslim_trees/output}"
export DERIVED_DIR="${DERIVED_DIR:-${BASE_OUT}/constant/${CHROM_DIR}/trees_derived}"
export VCF_DIR="${VCF_DIR:-${BASE_OUT}/constant/${CHROM_DIR}/vcf}"
export SIM_STATS_CORES="${SLURM_CPUS_PER_TASK}"
export SAMPLE_INDIVIDUALS="${SAMPLE_INDIVIDUALS:-1135}"
export SAMPLE_SEED="${SAMPLE_SEED:-1729}"
export FOLD_SFS="${FOLD_SFS:-1}"
export SNP_ONLY="${SNP_ONLY:-0}"

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
export SIM_STATS_OUT="${SIM_STATS_OUT:-${BASE_OUT}/simstats/chunks/sim_stats_${CHROM_DIR}_${START_ID}_${END_ID}.csv}"

mkdir -p "$(dirname "${SIM_STATS_OUT}")" /home/mlensink/slimsimulations/ABCslim_trees/logs/simstats_chunks

echo "Running calc_sim_stats.py for IDs ${START_ID}-${END_ID}"
echo "CHROM_DIR=${CHROM_DIR}"
echo "DERIVED_DIR=${DERIVED_DIR}"
echo "VCF_DIR=${VCF_DIR}"
echo "SIM_STATS_CORES=${SIM_STATS_CORES}"
echo "SAMPLE_INDIVIDUALS=${SAMPLE_INDIVIDUALS}"
echo "SAMPLE_SEED=${SAMPLE_SEED}"
echo "FOLD_SFS=${FOLD_SFS}"
echo "SNP_ONLY=${SNP_ONLY}"
echo "Writing output to ${SIM_STATS_OUT}"

python ./scripts/calc_sim_stats.py
