#!/bin/bash
#SBATCH --job-name=derive_trees
#SBATCH --output=/home/mlensink/ABCslim_trees/logs/derive_trees/%x_%A_%a.out
#SBATCH --error=/home/mlensink/ABCslim_trees/logs/derive_trees/%x_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=500M
#SBATCH --time=7:00:00
#SBATCH --partition=low
#SBATCH --account=publicgrp
#SBATCH --array=0-49999
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=melensink@ucdavis.edu

set -euo pipefail

: "${PROJECT_ROOT:=/home/mlensink/ABCslim_trees}"
: "${ABC_DATA_ROOT:=/group/gmonroegrp2/mlensink/ABC_data}"
: "${PYTHON:=/group/gmonroegrp2/mlensink/conda_envs/abc_trees/bin/python}"
: "${PARAM_FILE:=${PROJECT_ROOT}/priors/priors_2.6.25.csv}"
: "${CHROM_DIR:?Set CHROM_DIR to chrom2, chrom3, or chrom4}"
: "${OFFSET:=0}"

: "${RECAP:=1}"
: "${NE_RECAP:=6000}"
: "${RECOMB_RECAP:=2e-6}"
: "${Q:=50}"

case "$CHROM_DIR" in
  chrom2)
    : "${GENE_CSV:=${PROJECT_ROOT}/genomeinfo/gene100_chr2_18000000.csv}"
    : "${INTERGENE_CSV:=${PROJECT_ROOT}/genomeinfo/intergene100_chr2_18000000.csv}"
    ;;
  chrom3)
    : "${GENE_CSV:=${PROJECT_ROOT}/genomeinfo/gene100_chr3_22000000.csv}"
    : "${INTERGENE_CSV:=${PROJECT_ROOT}/genomeinfo/intergene100_chr3_22000000.csv}"
    ;;
  chrom4)
    : "${GENE_CSV:=${PROJECT_ROOT}/genomeinfo/gene100_chr4_2000000.csv}"
    : "${INTERGENE_CSV:=${PROJECT_ROOT}/genomeinfo/intergene100_chr4_2000000.csv}"
    ;;
  *)
    echo "Unsupported CHROM_DIR=$CHROM_DIR; expected chrom2, chrom3, or chrom4" >&2
    exit 2
    ;;
esac

: "${TREES_DIR:=${ABC_DATA_ROOT}/${CHROM_DIR}/trees}"
: "${DERIVED_DIR:=${ABC_DATA_ROOT}/${CHROM_DIR}/trees_derived}"

mkdir -p "$DERIVED_DIR"

LINE=$((SLURM_ARRAY_TASK_ID + OFFSET + 2))
data=$(sed -n "${LINE}p" "$PARAM_FILE" || true)

if [[ -z "${data}" ]]; then
  echo "No data at line ${LINE} in ${PARAM_FILE}. Exiting."
  exit 0
fi

IFS=',' read -r ID gmu imu gd id gdfe idfe <<< "$data"

TREES_IN="${TREES_DIR}/${ID}.trees"
DERIVED_OUT="${DERIVED_DIR}/${ID}.trees"

echo "JOB=$SLURM_JOB_ID TASK=$SLURM_ARRAY_TASK_ID LINE=$LINE ID=$ID"
echo "CHROM_DIR=$CHROM_DIR OFFSET=$OFFSET"
echo "TREES_IN=$TREES_IN"
echo "DERIVED_OUT=$DERIVED_OUT"
echo "GENE_CSV=$GENE_CSV"
echo "INTERGENE_CSV=$INTERGENE_CSV"
echo "Q=$Q RECAP=$RECAP NE_RECAP=$NE_RECAP RECOMB_RECAP=$RECOMB_RECAP"

if [[ ! -s "$TREES_IN" ]]; then
  echo "No non-empty raw tree: $TREES_IN"
  exit 0
fi

if [[ -s "$DERIVED_OUT" ]]; then
  echo "Derived tree already exists: $DERIVED_OUT"
  exit 0
fi

REC_ARGS=()
if [[ "$RECAP" == "1" ]]; then
  REC_ARGS=(--recapitate --Ne "$NE_RECAP" --recomb "$RECOMB_RECAP")
fi

"$PYTHON" "${PROJECT_ROOT}/scripts/trees_to_vcf.py" \
  --trees "$TREES_IN" \
  --out_trees "$DERIVED_OUT" \
  --gene_csv "$GENE_CSV" \
  --intergene_csv "$INTERGENE_CSV" \
  --Q "$Q" \
  --gmu "$gmu" --imu "$imu" --gd "$gd" --id "$id" \
  --seed "$ID" \
  "${REC_ARGS[@]}"
