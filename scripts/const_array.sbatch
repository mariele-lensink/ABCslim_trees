#!/bin/bash
#SBATCH --job-name=const_Q25_chr5
#SBATCH --output=/home/mlensink/slimsimulations/ABCslim_trees/logs/const_trees/chrom5_Q25/%x_%A_%a.out
#SBATCH --error=/home/mlensink/slimsimulations/ABCslim_trees/logs/const_trees/chrom5_Q25/%x_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=8G
#SBATCH --time=24:00:00
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
: "${PARAM_FILE:=/home/mlensink/slimsimulations/ABCslim_trees/priors/priors_Q25_smax2e-2_2026-05-05.csv}"
: "${OFFSET:=0}"                  # 0 => first row after header is array task 0
: "${END:=120000}"                # constant burn-in length (10*Nanc for Q=25)
: "${Q:=25}"                      # scaling factor; NPOP defaults to NE_REAL/Q inside SLiM
: "${SIMPLIFY_INTERVAL:=200}"
: "${SELF:=0.98}"
: "${CHUNK_SIZE:=50000}"          # Maximum jobs per array submission
: "${CHUNK_INDEX:=0}"             # Index of the current chunk (0-based)

: "${SCRIPT_CONST:=/home/mlensink/slimsimulations/ABCslim_trees/models/ABCtrees.slim}"

: "${BASE_OUT:=/group/gmonroegrp2/mlensink/ABC_data}"
: "${CHROM_DIR:=chrom5_Q25_constant}"
: "${OUT_CONST:=${BASE_OUT}/${CHROM_DIR}/trees}"
: "${GENEFILE:=/home/mlensink/slimsimulations/ABCslim_trees/genomeinfo/gene100_chr5_1.csv}"
: "${INTERFILE:=/home/mlensink/slimsimulations/ABCslim_trees/genomeinfo/intergene100_chr5_1.csv}"
LOGDIR=/home/mlensink/slimsimulations/ABCslim_trees/logs/const_trees/chrom5_Q25

mkdir -p "$OUT_CONST" "$LOGDIR"

# ----------------------------
# Adjust task ID for chunk. Keep this zero-based so SLURM array task 0 reads
# the first data row after the CSV header.
# ----------------------------
TASK_ID=$((SLURM_ARRAY_TASK_ID + CHUNK_INDEX * CHUNK_SIZE))

# ----------------------------
# Pick the row for this task
# ----------------------------
# CSV has header on line 1. First data row is line 2.
LINE=$((TASK_ID + OFFSET + 2))
data=$(sed -n "${LINE}p" "$PARAM_FILE" || true)

if [[ -z "${data}" ]]; then
  echo "No data at line ${LINE} in ${PARAM_FILE}. Exiting."
  exit 0
fi

# Expected columns: ID,gmu,imu,gd,id,gdfe,idfe
IFS=',' read -r ID gmu imu gd id gdfe idfe <<< "$data"

# ----------------------------
# Check if output file exists
# ----------------------------
expected="${OUT_CONST}/${ID}.trees"
if [[ -s "$expected" ]]; then
  echo "Output file already exists: $expected. Skipping."
  exit 0
fi

# ----------------------------
# Run constant model
# ----------------------------
echo "JOB=$SLURM_JOB_ID ARRAY_TASK=$SLURM_ARRAY_TASK_ID TASK=$TASK_ID LINE=$LINE ID=$ID END=$END Q=$Q"
echo "Params: gmu=$gmu imu=$imu gd=$gd id=$id gdfe=$gdfe idfe=$idfe"
echo "OUTTREES=$OUT_CONST"

SLIM_ARGS=(
  -d "ID=${ID}"
  -d "OUTTREES='${OUT_CONST}'"
  -d "END=${END}"
  -d "Q=${Q}"
  -d "SIMPLIFY_INTERVAL=${SIMPLIFY_INTERVAL}"
  -d "SELF=${SELF}"
  -d "gmu=${gmu}" -d "imu=${imu}" -d "gd=${gd}" -d "id=${id}" -d "gdfe=${gdfe}" -d "idfe=${idfe}"
)

if [[ -n "$GENEFILE" ]]; then
  SLIM_ARGS+=(-d "GENEFILE='${GENEFILE}'")
fi

if [[ -n "$INTERFILE" ]]; then
  SLIM_ARGS+=(-d "INTERFILE='${INTERFILE}'")
fi

# Write SLiM stdout/stderr to a per-ID log file (keeps array .out clean)
slim \
  "${SLIM_ARGS[@]}" \
  "${SCRIPT_CONST}" \
  > "${LOGDIR}/${ID}_const_Q${Q}_END${END}.log" 2>&1

# ----------------------------
# Fail loudly if output missing
# ----------------------------
if [[ ! -s "$expected" ]]; then
  echo "ERROR: expected trees file missing or empty: $expected"
  echo "Last 50 lines of SLiM log:"
  tail -n 50 "${LOGDIR}/${ID}_const_Q${Q}_END${END}.log" || true
  exit 1
fi

echo "Done: ID=$ID Q=$Q END=$END trees=${expected}"
