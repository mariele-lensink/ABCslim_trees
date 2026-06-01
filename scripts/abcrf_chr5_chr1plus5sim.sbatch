#!/bin/bash
#SBATCH --job-name=abcrf
#SBATCH --output=/home/mlensink/slimsimulations/ABCslim_trees/logs/ABC_results/abcrf_%j.out
#SBATCH --error=/home/mlensink/slimsimulations/ABCslim_trees/logs/ABC_results/abcrf_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --partition=bmh
#SBATCH --account=gmonroegrp

set -euo pipefail
source /group/gmonroegrp2/mlensink/miniconda3/etc/profile.d/conda.sh
conda activate abc_trees

# Run the ABCRF script
cd /home/mlensink/slimsimulations/ABCslim_trees
echo "ABCRF_STAT_SET=${ABCRF_STAT_SET:-all_variable}"
echo "ABCRF_CHROM=${ABCRF_CHROM:-chrom5}"
echo "ABCRF_NCORES=${ABCRF_NCORES:-${SLURM_CPUS_PER_TASK:-4}}"
Rscript scripts/abcrf_chr5_chr1plus5sim.R
