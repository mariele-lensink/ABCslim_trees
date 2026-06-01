#!/bin/bash
#SBATCH --job-name=conv_perf
#SBATCH --output=/home/mlensink/ABCslim_trees/logs/convergence_performance/%x_%j.out
#SBATCH --error=/home/mlensink/ABCslim_trees/logs/convergence_performance/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --partition=bmh
#SBATCH --account=gmonroegrp
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=melensink@ucdavis.edu

set -euo pipefail

: "${PROJECT_ROOT:=/home/mlensink/ABCslim_trees}"
: "${R_SCRIPT:=/group/gmonroegrp2/mlensink/conda_envs/abc_trees/bin/Rscript}"
: "${ABC_CONV_OUTPUT_DIR:=/group/gmonroegrp2/mlensink/ABC_data/convergence_test_chrom5_100_abcrf_all_variable}"
: "${ABC_CONV_PERF_OUTPUT_DIR:=${ABC_CONV_OUTPUT_DIR}/performance}"
: "${ABC_CONV_PERF_FIGURE_DIR:=${PROJECT_ROOT}/figures/abc_convergence_chrom5_100_abcrf_all_variable/performance}"
: "${ABC_CONV_PRIORS_FILE:=${PROJECT_ROOT}/priors/priors_2.6.25.csv}"

mkdir -p "${PROJECT_ROOT}/logs/convergence_performance"

cd "$PROJECT_ROOT"

echo "PROJECT_ROOT=$PROJECT_ROOT"
echo "ABC_CONV_OUTPUT_DIR=$ABC_CONV_OUTPUT_DIR"
echo "ABC_CONV_PERF_OUTPUT_DIR=$ABC_CONV_PERF_OUTPUT_DIR"
echo "ABC_CONV_PERF_FIGURE_DIR=$ABC_CONV_PERF_FIGURE_DIR"
echo "ABC_CONV_PRIORS_FILE=$ABC_CONV_PRIORS_FILE"
echo "R_SCRIPT=$R_SCRIPT"

export ABC_CONV_OUTPUT_DIR
export ABC_CONV_PERF_OUTPUT_DIR
export ABC_CONV_PERF_FIGURE_DIR
export ABC_CONV_PRIORS_FILE

"$R_SCRIPT" scripts/summarize_convergence_performance.R
