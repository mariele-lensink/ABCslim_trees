#!/bin/bash
#SBATCH --job-name=trees_vcf
#SBATCH --output=/home/mlensink/slimsimulations/ABCslim_trees/logs/const_trees_to_vcf/%x_%A_%a.out
#SBATCH --error=/home/mlensink/slimsimulations/ABCslim_trees/logs/const_trees_to_vcf/%x_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=8G
#SBATCH --time=7:00:00
#SBATCH --partition=bmh
#SBATCH --account=gmonroegrp
#SBATCH --array=0-49999
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=melensink@ucdavis.edu

set -euo pipefail

source ~/miniconda3/etc/profile.d/conda.sh
conda activate abc_trees_vcf
: "${PARAM_FILE:=/home/mlensink/slimsimulations/ABCslim_trees/priors/priors_2.6.25.csv}"
: "${OFFSET:=0}"   # set to 0 for first half, 50000 for second half, etc.

: "${BASE_OUT:=/home/mlensink/slimsimulations/ABCslim_trees/output}"
: "${CHROM_DIR:=chrom1}"
: "${TREES_DIR:=${BASE_OUT}/constant/${CHROM_DIR}/trees}"
: "${VCF_DIR:=${BASE_OUT}/constant/${CHROM_DIR}/vcf}"
: "${DERIVED_DIR:=${BASE_OUT}/constant/${CHROM_DIR}/trees_derived}"

mkdir -p "$DERIVED_DIR"

# These must point to the interval CSVs your script expects (2-col start,stop; no header)
: "${GENE_CSV:=/home/mlensink/slimsimulations/ABCslim_trees/genomeinfo/gene100_chr5_1.csv}"
: "${INTERGENE_CSV:=/home/mlensink/slimsimulations/ABCslim_trees/genomeinfo/intergene100_chr5_1.csv}"
: "${L:=284122}"

# Recapitation parameters (used only if RECAP=1)
: "${RECAP:=1}"
: "${NE_RECAP:=6000}"
: "${RECOMB_RECAP:=2e-6}"
: "${Q:=50}"
: "${WRITE_VCF:=1}"

LOGDIR=/home/mlensink/slimsimulations/ABCslim_trees/logs/const_trees_to_vcf
mkdir -p "$LOGDIR" "$VCF_DIR"

LINE=$((SLURM_ARRAY_TASK_ID + OFFSET + 2))   # +2 skips header
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

echo "JOB=$SLURM_JOB_ID TASK=$SLURM_ARRAY_TASK_ID LINE=$LINE ID=$ID"
echo "IN=$TREES_IN"
echo "OUT=$VCF_OUT"
echo "DERIVED_OUT=$DERIVED_OUT"
echo "Q=$Q RECAP=$RECAP Ne_recap=$NE_RECAP recomb_recap=$RECOMB_RECAP"
echo "WRITE_VCF=$WRITE_VCF"


REC_ARGS=()
if [[ "$RECAP" == "1" ]]; then
  REC_ARGS=(--recapitate --Ne "$NE_RECAP" --recomb "$RECOMB_RECAP")
fi

VCF_ARGS=()
if [[ "$WRITE_VCF" == "1" ]]; then
  VCF_ARGS=(--out_vcf "$VCF_OUT")
fi

python /home/mlensink/slimsimulations/ABCslim_trees/scripts/trees_to_vcf.py \
  --trees "$TREES_IN" \
  --out_trees "$DERIVED_OUT" \
  --gene_csv "$GENE_CSV" \
  --intergene_csv "$INTERGENE_CSV" \
  --Q "$Q" \
  --L "$L" \
  --gmu "$gmu" --imu "$imu" --gd "$gd" --id "$id" \
  --seed "$ID" \
  "${VCF_ARGS[@]}" \
  "${REC_ARGS[@]}"
