#!/bin/bash
#SBATCH --job-name=abc_conv
#SBATCH --output=/home/mlensink/ABCslim_trees/logs/ABC_results/abc_convergence_%j.out
#SBATCH --error=/home/mlensink/ABCslim_trees/logs/ABC_results/abc_convergence_%j.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --partition=bmh
#SBATCH --account=gmonroegrp

set -euo pipefail
source /group/gmonroegrp2/mlensink/miniconda3/etc/profile.d/conda.sh
conda activate /group/gmonroegrp2/mlensink/conda_envs/abc_trees

cd /home/mlensink/ABCslim_trees
mkdir -p logs/ABC_results



export ABC_REPO_DIR="${ABC_REPO_DIR:-/home/mlensink/ABCslim_trees}"
export ABC_CONV_SIM_STATS_FILE="${ABC_CONV_SIM_STATS_FILE:-../../../group/gmonroegrp2/mlensink/ABC_data/chrom5/sim_stats_v2_chrom5_merged.csv}"
export ABC_CONV_PRIORS_FILE="${ABC_CONV_PRIORS_FILE:-priors/priors_2.6.25.csv}"
export ABC_CONV_OUTPUT_DIR="${ABC_CONV_OUTPUT_DIR:-../../../group/gmonroegrp2/mlensink/ABC_data/convergence_test}"
export ABC_CONV_FIGURE_DIR="${ABC_CONV_FIGURE_DIR:-${HOME}/figures/abc_convergence}"
export ABC_CONV_N_REPS="${ABC_CONV_N_REPS:-10}"
export ABC_CONV_SEED="${ABC_CONV_SEED:-20260526}"
export ABC_CONV_REJECTION_TOL="${ABC_CONV_REJECTION_TOL:-0.01}"
export ABC_CONV_RF_NTREE="${ABC_CONV_RF_NTREE:-500}"
export ABCRF_NCORES="${ABCRF_NCORES:-${SLURM_CPUS_PER_TASK:-4}}"

echo "ABC_REPO_DIR=${ABC_REPO_DIR}"
echo "ABC_CONV_SIM_STATS_FILE=${ABC_CONV_SIM_STATS_FILE}"
echo "ABC_CONV_PRIORS_FILE=${ABC_CONV_PRIORS_FILE}"
echo "ABC_CONV_OUTPUT_DIR=${ABC_CONV_OUTPUT_DIR}"
echo "ABC_CONV_FIGURE_DIR=${ABC_CONV_FIGURE_DIR}"
echo "ABC_CONV_N_REPS=${ABC_CONV_N_REPS}"
echo "ABC_CONV_SEED=${ABC_CONV_SEED}"
echo "ABC_CONV_REJECTION_TOL=${ABC_CONV_REJECTION_TOL}"
echo "ABC_CONV_RF_NTREE=${ABC_CONV_RF_NTREE}"
echo "ABCRF_NCORES=${ABCRF_NCORES}"

Rscript scripts/convergence_test_abc.R
Rscript scripts/plot_convergence_recovery.R
