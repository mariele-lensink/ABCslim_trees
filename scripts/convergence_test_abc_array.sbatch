#!/bin/bash
#SBATCH --job-name=abc_conv_abcrf
#SBATCH --output=/home/mlensink/ABCslim_trees/logs/ABC_results/%x_%A_%a.out
#SBATCH --error=/home/mlensink/ABCslim_trees/logs/ABC_results/%x_%A_%a.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=192G
#SBATCH --array=1-100%14
#SBATCH --partition=bml
#SBATCH --account=gmonroegrp

set -euo pipefail
source /group/gmonroegrp2/mlensink/miniconda3/etc/profile.d/conda.sh
conda activate /group/gmonroegrp2/mlensink/conda_envs/abc_trees

cd /home/mlensink/ABCslim_trees
mkdir -p logs/ABC_results



export ABC_REPO_DIR="${ABC_REPO_DIR:-/home/mlensink/ABCslim_trees}"
export ABC_CONV_CHROM="${ABC_CONV_CHROM:-chrom1}"
export ABC_CONV_SIM_STATS_FILE="${ABC_CONV_SIM_STATS_FILE:-/group/gmonroegrp2/mlensink/ABC_data/${ABC_CONV_CHROM}/sim_stats_v2_${ABC_CONV_CHROM}_merged.csv}"
export ABC_CONV_PRIORS_FILE="${ABC_CONV_PRIORS_FILE:-priors/priors_2.6.25.csv}"
export ABC_CONV_OUTPUT_DIR="${ABC_CONV_OUTPUT_DIR:-/group/gmonroegrp2/mlensink/ABC_data/convergence_test_${ABC_CONV_CHROM}_100_abcrf_all_variable}"
export ABC_CONV_FIGURE_DIR="${ABC_CONV_FIGURE_DIR:-${HOME}/figures/abc_convergence_${ABC_CONV_CHROM}_100_abcrf_all_variable}"
export ABC_CONV_N_REPS="${ABC_CONV_N_REPS:-100}"
export ABC_CONV_SEED="${ABC_CONV_SEED:-20260526}"
export ABC_CONV_REJECTION_TOL="${ABC_CONV_REJECTION_TOL:-0.01}"
export ABC_CONV_RF_NTREE="${ABC_CONV_RF_NTREE:-500}"
export ABC_CONV_RUNS="${ABC_CONV_RUNS:-abcrf_all_variable}"
export ABCRF_NCORES="${ABCRF_NCORES:-${SLURM_CPUS_PER_TASK:-4}}"

echo "ABC_REPO_DIR=${ABC_REPO_DIR}"
echo "ABC_CONV_CHROM=${ABC_CONV_CHROM}"
echo "ABC_CONV_SIM_STATS_FILE=${ABC_CONV_SIM_STATS_FILE}"
echo "ABC_CONV_PRIORS_FILE=${ABC_CONV_PRIORS_FILE}"
echo "ABC_CONV_OUTPUT_DIR=${ABC_CONV_OUTPUT_DIR}"
echo "ABC_CONV_FIGURE_DIR=${ABC_CONV_FIGURE_DIR}"
echo "ABC_CONV_N_REPS=${ABC_CONV_N_REPS}"
echo "ABC_CONV_SEED=${ABC_CONV_SEED}"
echo "ABC_CONV_REJECTION_TOL=${ABC_CONV_REJECTION_TOL}"
echo "ABC_CONV_RF_NTREE=${ABC_CONV_RF_NTREE}"
echo "ABC_CONV_RUNS=${ABC_CONV_RUNS}"
echo "ABCRF_NCORES=${ABCRF_NCORES}"

export ABC_CONV_REP_INDEX="${SLURM_ARRAY_TASK_ID}"
Rscript scripts/convergence_test_abc_array.R
