#!/bin/bash
#SBATCH -D /home/mlensink/slimsimulations/ABCslim_trees/
#SBATCH -o /home/mlensink/slimsimulations/ABCslim_trees/logs/ABC_results/observed_1001_%j.out
#SBATCH -e /home/mlensink/slimsimulations/ABCslim_trees/logs/ABC_results/observed_1001_%j.err
#SBATCH -J obs1001stats
#SBATCH -t 1-00:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --account=gmonroegrp
#SBATCH --partition=bml
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=melensink@ucdavis.edu

set -euo pipefail

source /home/mlensink/miniconda3/etc/profile.d/conda.sh
conda activate /home/mlensink/miniconda3/envs/abc_trees_vcf
echo "Python: $(which python)"

export CHROM_DIR="${1:-${CHROM_DIR:-chrom1}}"
if [[ -z "${VCF_CHROM:-}" ]]; then
  case "${CHROM_DIR}" in
    chrom1) export VCF_CHROM="1" ;;
    chrom5) export VCF_CHROM="5" ;;
    *) echo "Set VCF_CHROM for CHROM_DIR=${CHROM_DIR}" >&2; exit 1 ;;
  esac
else
  export VCF_CHROM
fi
export VCF_FILE="${VCF_FILE:-/home/mlensink/slimsimulations/ABCslim_trees/../../rawdata/arabidopsis/1001genomes_snp-short-indel_only_ACGTN.vcf.gz}"
export VCF_FILE_ID="${VCF_FILE_ID:-1001}"
export SIM_STATS_CORES="${SLURM_CPUS_PER_TASK}"
export SAMPLE_INDIVIDUALS="${SAMPLE_INDIVIDUALS:-1135}"
export SAMPLE_SEED="${SAMPLE_SEED:-1729}"
export FOLD_SFS="${FOLD_SFS:-1}"
export SNP_ONLY="${SNP_ONLY:-1}"
export REQUIRE_COMPLETE_VCF_SAMPLES="${REQUIRE_COMPLETE_VCF_SAMPLES:-0}"
export MIN_VCF_CALL_FRACTION="${MIN_VCF_CALL_FRACTION:-0.8}"
export BASE_OUT="${BASE_OUT:-/home/mlensink/slimsimulations/ABCslim_trees/output}"
export SIM_STATS_OUT="${SIM_STATS_OUT:-${BASE_OUT}/ABC_results/observed/target_stats_${CHROM_DIR}.csv}"
export OBS_TMP_DIR="${OBS_TMP_DIR:-${SLURM_TMPDIR:-/tmp}/abc_obs_1001_${SLURM_JOB_ID:-manual}_${CHROM_DIR}}"

mkdir -p "$(dirname "${SIM_STATS_OUT}")" /home/mlensink/slimsimulations/ABCslim_trees/logs/ABC_results
mkdir -p "${OBS_TMP_DIR}"
trap 'rm -rf "${OBS_TMP_DIR}"' EXIT

subset_vcf="${OBS_TMP_DIR}/1001_${CHROM_DIR}.vcf.gz"

echo "Running observed 1001 VCF simstats"
echo "CHROM_DIR=${CHROM_DIR}"
echo "VCF_CHROM=${VCF_CHROM}"
echo "VCF_FILE=${VCF_FILE}"
echo "VCF_FILE_ID=${VCF_FILE_ID}"
echo "SAMPLE_INDIVIDUALS=${SAMPLE_INDIVIDUALS}"
echo "SAMPLE_SEED=${SAMPLE_SEED}"
echo "FOLD_SFS=${FOLD_SFS}"
echo "SNP_ONLY=${SNP_ONLY}"
echo "REQUIRE_COMPLETE_VCF_SAMPLES=${REQUIRE_COMPLETE_VCF_SAMPLES}"
echo "MIN_VCF_CALL_FRACTION=${MIN_VCF_CALL_FRACTION}"
echo "OBS_TMP_DIR=${OBS_TMP_DIR}"
echo "Writing output to ${SIM_STATS_OUT}"

bcftools view -r "${VCF_CHROM}" -Oz -o "${subset_vcf}" "${VCF_FILE}"

export VCF_FILE="${subset_vcf}"
export VCF_CHROM=""
python ./scripts/calc_sim_stats.py
