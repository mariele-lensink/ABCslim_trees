#!/bin/bash
#SBATCH -D /home/mlensink/slimsimulations/ABCslim_trees/
#SBATCH -o /home/mlensink/slimsimulations/ABCslim_trees/logs/simstats/simstats_stdout-%j.txt
#SBATCH -e /home/mlensink/slimsimulations/ABCslim_trees/logs/simstats/simstats_stderr-%j.txt
#SBATCH -J abc_slim
#SBATCH -t 7:00:00
#SBATCH --mem 64G
#SBATCH --cpus-per-task=16
#SBATCH --account=gmonroegrp
#SBATCH --partition=bmh
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=melensink@ucdavis.edu

set -euo pipefail

source /home/mlensink/miniconda3/etc/profile.d/conda.sh
conda activate /home/mlensink/miniconda3/envs/abc_trees_vcf
echo "Python: $(which python)"

export CHROM_DIR="${CHROM_DIR:-chrom1}"
export SIM_STATS_CORES="${SLURM_CPUS_PER_TASK}"
export START_ID="${START_ID:-1}"
export END_ID="${END_ID:-50000}"
export BASE_OUT="${BASE_OUT:-/home/mlensink/slimsimulations/ABCslim_trees/output}"
export DERIVED_DIR="${DERIVED_DIR:-${BASE_OUT}/constant/${CHROM_DIR}/trees_derived}"
export VCF_DIR="${VCF_DIR:-${BASE_OUT}/constant/${CHROM_DIR}/vcf}"
export SAMPLE_INDIVIDUALS="${SAMPLE_INDIVIDUALS:-1135}"
export SAMPLE_SEED="${SAMPLE_SEED:-1729}"
export FOLD_SFS="${FOLD_SFS:-1}"
export SNP_ONLY="${SNP_ONLY:-0}"
export SIM_STATS_OUT="${SIM_STATS_OUT:-${BASE_OUT}/simstats/sim_stats_${CHROM_DIR}_${START_ID}_${END_ID}.csv}"

mkdir -p "$(dirname "${SIM_STATS_OUT}")" /home/mlensink/slimsimulations/ABCslim_trees/logs/simstats

echo "Running calc_sim_stats.py for IDs ${START_ID}-${END_ID}"
echo "CHROM_DIR=${CHROM_DIR}"
echo "DERIVED_DIR=${DERIVED_DIR}"
echo "VCF_DIR=${VCF_DIR}"
echo "SAMPLE_INDIVIDUALS=${SAMPLE_INDIVIDUALS}"
echo "SAMPLE_SEED=${SAMPLE_SEED}"
echo "FOLD_SFS=${FOLD_SFS}"
echo "SNP_ONLY=${SNP_ONLY}"
echo "Writing output to ${SIM_STATS_OUT}"

python ./scripts/calc_sim_stats.py
