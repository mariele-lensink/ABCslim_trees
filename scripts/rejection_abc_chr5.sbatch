#!/bin/bash
#SBATCH --job-name=rejection_abc
#SBATCH --output=/home/mlensink/slimsimulations/ABCslim_trees/logs/ABC_results/rejection_abc_%j.out
#SBATCH --error=/home/mlensink/slimsimulations/ABCslim_trees/logs/ABC_results/rejection_abc_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --partition=bmh
#SBATCH --account=gmonroegrp

set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh
conda activate abc_trees

cd /home/mlensink/slimsimulations/ABCslim_trees
echo "ABCREJ_CHROM=${ABCREJ_CHROM:-chrom5}"
Rscript scripts/rejection_abc_chr5.R
