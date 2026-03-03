#!/bin/bash
#SBATCH --job-name=const_trees
#SBATCH --output=/home/mlensink/slimsimulations/ABCslim_trees/logs/const_trees/%x_%A_%a.out
#SBATCH --error=/home/mlensink/slimsimulations/ABCslim_trees/logs/const_trees/%x_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=8G
#SBATCH --time=2-00:00:00
#SBATCH --partition=bml
#SBATCH --account=gmonroegrp
#SBATCH --array=0-49999
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=melensink@ucdavis.edu

set -euo pipefail

source ~/miniconda3/etc/profile.d/conda.sh
conda activate abc_trees

# ----------------------------
# User-configurable knobs (override via sbatch --export=...)
# ----------------------------
: "${PARAM_FILE:=/home/mlensink/slimsimulations/ABCslim_trees/priors/priors_2.6.25.csv}"
: "${OFFSET:=0}"                  # 0 => first row after header is array task 0
: "${END:=60000}"                 # constant burn-in length (10*Nanc for Nanc=6000)
: "${SIMPLIFY_INTERVAL:=200}"
: "${SELF:=0.98}"

: "${SCRIPT_CONST:=/home/mlensink/slimsimulations/ABCslim_trees/models/ABCtrees.slim}"

BASE_OUT=/home/mlensink/slimsimulations/ABCslim_trees/output/
OUT_CONST="${BASE_OUT}constant/trees"
LOGDIR=/home/mlensink/slimsimulations/ABCslim_trees/logs/const_trees

mkdir -p "$OUT_CONST" "$LOGDIR"

# ----------------------------
# Pick the row for this task
# ----------------------------
# CSV has header on line 1. First data row is line 2.
LINE=$((SLURM_ARRAY_TASK_ID + OFFSET + 2))
data=$(sed -n "${LINE}p" "$PARAM_FILE" || true)

if [[ -z "${data}" ]]; then
  echo "No data at line ${LINE} in ${PARAM_FILE}. Exiting."
  exit 0
fi

# Expected columns: ID,gmu,imu,gd,id,gdfe,idfe
IFS=',' read -r ID gmu imu gd id gdfe idfe <<< "$data"

# ----------------------------
# Run constant model
# ----------------------------
echo "JOB=$SLURM_JOB_ID TASK=$SLURM_ARRAY_TASK_ID LINE=$LINE ID=$ID END=$END"
echo "Params: gmu=$gmu imu=$imu gd=$gd id=$id gdfe=$gdfe idfe=$idfe"
echo "OUTTREES=$OUT_CONST"

# Write SLiM stdout/stderr to a per-ID log file (keeps array .out clean)
slim \
  -d "ID=${ID}" \
  -d "OUTTREES='${OUT_CONST}'" \
  -d "END=${END}" \
  -d "SIMPLIFY_INTERVAL=${SIMPLIFY_INTERVAL}" \
  -d "SELF=${SELF}" \
  -d "gmu=${gmu}" -d "imu=${imu}" -d "gd=${gd}" -d "id=${id}" -d "gdfe=${gdfe}" -d "idfe=${idfe}" \
  "${SCRIPT_CONST}" \
  > "${LOGDIR}/${ID}_const_END${END}.log" 2>&1

# ----------------------------
# Fail loudly if output missing
# ----------------------------
expected="${OUT_CONST}/${ID}.trees"
if [[ ! -s "$expected" ]]; then
  echo "ERROR: expected trees file missing or empty: $expected"
  echo "Last 50 lines of SLiM log:"
  tail -n 50 "${LOGDIR}/${ID}_const_END${END}.log" || true
  exit 1
fi

echo "Done: ID=$ID END=$END trees=${expected}"
