#!/bin/bash
#SBATCH --job-name=trees_vcf_tajima
#SBATCH --output=/home/mlensink/slimsimulations/ABCslim_trees/logs/trees_vcf_tajima/%x_%A_%a.out
#SBATCH --error=/home/mlensink/slimsimulations/ABCslim_trees/logs/trees_vcf_tajima/%x_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=30:00
#SBATCH --partition=bmh
#SBATCH --account=gmonroegrp
#SBATCH --array=0-49999
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=melensink@ucdavis.edu

set -euo pipefail

# Load conda in a way that works with either miniconda3 or anaconda3 installs.
if [[ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
  source "$HOME/miniconda3/etc/profile.d/conda.sh"
elif [[ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]]; then
  source "$HOME/anaconda3/etc/profile.d/conda.sh"
else
  echo "Could not find conda.sh in ~/miniconda3 or ~/anaconda3"
  exit 1
fi

: "${PIPELINE_ENV:=slim_pipeline}"
conda activate "$PIPELINE_ENV"

: "${PARAM_FILE:=/home/mlensink/slimsimulations/ABCslim_trees/priors/priors_2.6.25.csv}"
: "${OFFSET:=0}"

: "${BASE_OUT:=/home/mlensink/slimsimulations/ABCslim_trees/output}"
: "${CHROM_DIR:=chrom1}"
: "${TREES_DIR:=${BASE_OUT}/constant/${CHROM_DIR}/trees}"
: "${VCF_DIR:=${BASE_OUT}/constant/${CHROM_DIR}/vcf}"
: "${DERIVED_DIR:=${BASE_OUT}/constant/${CHROM_DIR}/trees_derived}"
: "${TAJIMA_DIR:=${BASE_OUT}/constant/${CHROM_DIR}/tajima}"

: "${GENE_CSV:=/home/mlensink/slimsimulations/ABCslim_trees/genomeinfo/gene100_chr5_1.csv}"
: "${INTERGENE_CSV:=/home/mlensink/slimsimulations/ABCslim_trees/genomeinfo/intergene100_chr5_1.csv}"

: "${RECAP:=1}"
: "${NE_RECAP:=6000}"
: "${RECOMB_RECAP:=2e-6}"
: "${Q:=50}"
: "${TAJIMA_WINDOW:=100}"

LOGDIR=/home/mlensink/slimsimulations/ABCslim_trees/logs/trees_vcf_tajima
mkdir -p "$LOGDIR" "$DERIVED_DIR" "$VCF_DIR" "$TAJIMA_DIR"

LINE=$((SLURM_ARRAY_TASK_ID + OFFSET + 2))
data=$(sed -n "${LINE}p" "$PARAM_FILE" || true)

if [[ -z "${data}" ]]; then
  echo "No data at line ${LINE} in ${PARAM_FILE}. Exiting."
  exit 0
fi

# Expected columns: ID,gmu,imu,gd,id,gdfe,idfe
IFS=',' read -r ID gmu imu gd id gdfe idfe <<< "$data"

TREES_IN="${TREES_DIR}/${ID}.trees"
VCF_OUT="${VCF_DIR}/${ID}.vcf"
DERIVED_OUT="${DERIVED_DIR}/${ID}.trees"
TAJIMA_PREFIX="${TAJIMA_DIR}/${ID}"
TAJIMA_OUT="${TAJIMA_PREFIX}.Tajima.D"

echo "JOB=$SLURM_JOB_ID TASK=$SLURM_ARRAY_TASK_ID LINE=$LINE ID=$ID"

# 1) Tree must exist or job ends immediately.
if [[ ! -f "$TREES_IN" ]]; then
  echo "no tree: ${TREES_IN}"
  exit 0
fi

# 2) If Tajima output already exists, no work needed.
if [[ -f "$TAJIMA_OUT" ]]; then
  echo "already done: ${TAJIMA_OUT}"
  exit 0
fi

# 3) Generate derived trees + VCF only if either output is missing.
if [[ ! -f "$DERIVED_OUT" || ! -f "$VCF_OUT" ]]; then
  echo "Running trees_to_vcf for ID=${ID}"

  REC_ARGS=()
  if [[ "$RECAP" == "1" ]]; then
    REC_ARGS=(--recapitate --Ne "$NE_RECAP" --recomb "$RECOMB_RECAP")
  fi

  python /home/mlensink/slimsimulations/ABCslim_trees/scripts/trees_to_vcf.py \
    --trees "$TREES_IN" \
    --out_vcf "$VCF_OUT" \
    --out_trees "$DERIVED_OUT" \
    --gene_csv "$GENE_CSV" \
    --intergene_csv "$INTERGENE_CSV" \
    --Q "$Q" \
    --gmu "$gmu" --imu "$imu" --gd "$gd" --id "$id" \
    --seed "$ID" \
    "${REC_ARGS[@]}"
else
  echo "trees_derived and vcf already exist for ID=${ID}; skipping trees_to_vcf"
fi

# 4) Run TajimaD if not present yet.
if [[ -f "$TAJIMA_OUT" ]]; then
  echo "already done after checks: ${TAJIMA_OUT}"
  exit 0
fi

if [[ ! -f "$VCF_OUT" ]]; then
  echo "VCF missing after trees_to_vcf step: ${VCF_OUT}"
  exit 1
fi

vcftools --vcf "$VCF_OUT" --out "$TAJIMA_PREFIX" --TajimaD "$TAJIMA_WINDOW"

echo "Finished ID=${ID}"
