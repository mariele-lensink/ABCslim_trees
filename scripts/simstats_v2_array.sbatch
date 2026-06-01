#!/bin/bash
#SBATCH -D /home/mlensink/ABCslim_trees/
#SBATCH -o /home/mlensink/ABCslim_trees/logs/simstats_v2/%x_%A_%a.out
#SBATCH -e /home/mlensink/ABCslim_trees/logs/simstats_v2/%x_%A_%a.err
#SBATCH -J simstats_v2
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
conda activate /group/gmonroegrp2/mlensink/conda_envs/abc_trees

export CHROM_DIR="${CHROM_DIR:-chrom5}"
export MIN_ID="${MIN_ID:-1}"
export MAX_ID="${MAX_ID:-100000}"
export CHUNK_SIZE="${CHUNK_SIZE:-1000}"
export BASE_OUT="${BASE_OUT:-/group/gmonroegrp2/mlensink/ABC_data}"
export DERIVED_DIR="${DERIVED_DIR:-${BASE_OUT}/${CHROM_DIR}/trees_derived}"
export SIM_STATS_CORES="${SLURM_CPUS_PER_TASK}"
export SAMPLE_INDIVIDUALS="${SAMPLE_INDIVIDUALS:-1135}"
export SAMPLE_SEED="${SAMPLE_SEED:-1729}"
export FOLD_SFS="${FOLD_SFS:-1}"
export SNP_ONLY="${SNP_ONLY:-0}"

chunk_start=$((MIN_ID + SLURM_ARRAY_TASK_ID * CHUNK_SIZE))
chunk_end=$((chunk_start + CHUNK_SIZE - 1))
if (( chunk_end > MAX_ID )); then chunk_end="${MAX_ID}"; fi
if (( chunk_start > MAX_ID )); then exit 0; fi

export START_ID="${chunk_start}"
export END_ID="${chunk_end}"

chunk_dir="${BASE_OUT}/${CHROM_DIR}/simstats_2/chunks"
export SIM_STATS_OUT="${chunk_dir}/sim_stats_v2_${CHROM_DIR}_${START_ID}_${END_ID}.csv"

mkdir -p "$chunk_dir" /home/mlensink/ABCslim_trees/logs/simstats_v2

python ./scripts/calc_sim_stats_v2.py
