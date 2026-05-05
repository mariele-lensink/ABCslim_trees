#!/bin/bash
#SBATCH -D /home/mlensink/slimsimulations/ABCslim_trees/
#SBATCH -o /home/mlensink/slimsimulations/ABCslim_trees/logs/td_logs/quicktd_stdout-%j.txt
#SBATCH -e /home/mlensink/slimsimulations/ABCslim_trees/logs/td_logs/quicktd_stderr-%j.txt
#SBATCH -J abc_slim
#SBATCH -t 10-00:00:00
#SBATCH --mem 15G # Reduced memory per task
#SBATCH --cpus-per-task=64 # Reduced CPU request
#SBATCH --account=gmonroegrp
#SBATCH --partition=bmh
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=melensink@ucdavis.edu

bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate slimsims

# Run the updated script
./scripts/vcf_calc_td_1.26.2025.sh
