#!/bin/bash
#SBATCH --job-name=abcrf
#SBATCH --output=/home/mlensink/ABCslim_trees/logs/ABC_results/abcrf_%j.out
#SBATCH --error=/home/mlensink/ABCslim_trees/logs/ABC_results/abcrf_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --partition=bmh
#SBATCH --account=gmonroegrp

set -euo pipefail
source /group/gmonroegrp2/mlensink/miniconda3/etc/profile.d/conda.sh
conda activate /group/gmonroegrp2/mlensink/conda_envs/abc_trees

# Run the ABCRF script
cd /home/mlensink/ABCslim_trees
echo "ABCRF_STAT_SET=${ABCRF_STAT_SET:-all_variable}"
echo "ABCRF_CHROM=${ABCRF_CHROM:-chrom5}"
echo "ABCRF_NCORES=${ABCRF_NCORES:-${SLURM_CPUS_PER_TASK:-4}}"
echo "ABC_UNCAPPED_MEAN_FILTER=${ABC_UNCAPPED_MEAN_FILTER:-0}"
echo "ABC_STAT_VERSION=${ABC_STAT_VERSION:-v2}"
echo "ABC_DATA_ROOT=${ABC_DATA_ROOT:-/group/gmonroegrp2/mlensink/ABC_data}"
echo "ABC_OUTPUT_DIR=${ABC_OUTPUT_DIR:-/group/gmonroegrp2/mlensink/ABC_data/randomforestABC/ABC_results_statv2}"
echo "ABC_FIGURE_DIR=${ABC_FIGURE_DIR:-/group/gmonroegrp2/mlensink/ABC_data/randomforestABC/figures_statv2}"
echo "ABC_FILTER_Q=${ABC_FILTER_Q:-50}"
echo "ABC_FILTER_CAP=${ABC_FILTER_CAP:-1}"
Rscript scripts/abcrf_chr5.R
