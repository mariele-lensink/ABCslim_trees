#!/bin/bash
#SBATCH --job-name=bgs_test
#SBATCH --output=/home/mlensink/slimsimulations/ABCslim_trees/logs/bgs_test/%x_%A_%a.out
#SBATCH --error=/home/mlensink/slimsimulations/ABCslim_trees/logs/bgs_test/%x_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=8G
#SBATCH --time=2-00:00:00
#SBATCH --partition=bml
#SBATCH --account=gmonroegrp
#SBATCH --qos=gmonroegrp-bml-qos
#SBATCH --array=0-999%100
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=melensink@ucdavis.edu

set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh
conda activate abc_trees

# ----------------------------
# User-configurable knobs (override via sbatch --export=...)
# ----------------------------
: "${PARAM_FILE:=/home/mlensink/slimsimulations/ABCslim_trees/priors/priors_2.6.25.csv}"
: "${END:=20000}"                 # BGS test plan: 5000, 10000, 20000, 40000
: "${OFFSET:=0}"                  # if you want to start from later in the file
: "${SIMPLIFY_INTERVAL:=200}"
: "${SELF:=0.98}"

# Demography defaults as fractions of END (can override)
: "${BOT_FRAC:=0.2}"              # 80% reduction => keep 20%
: "${BOT_FRAC_TBOT:=0.50}"        # bottleneck time = 0.50*END
: "${BOT_FRAC_TREC:=0.75}"        # recovery time  = 0.75*END
: "${EXP_MULT:=3.0}"              # expansion multiplier at T_EXP
: "${EXP_FRAC_TEXP:=0.75}"        # expansion time = 0.75*END

# Script paths
: "${SCRIPT_CONST:=/home/mlensink/slimsimulations/ABCslim_trees/models/ABCtrees.slim}"
: "${SCRIPT_BOT:=/home/mlensink/slimsimulations/ABCslim_trees/models/ABCtrees_bottleneck.slim}"
: "${SCRIPT_EXP:=/home/mlensink/slimsimulations/ABCslim_trees/models/ABCtrees_expansion.slim}"

BASE_OUT=/home/mlensink/slimsimulations/ABCslim_trees/output/
LOGDIR=/home/mlensink/slimsimulations/ABCslim_trees/logs/bgs_test/
mkdir -p "$LOGDIR"

# ----------------------------
# Pick the row for this task
# ----------------------------
LINE=$((SLURM_ARRAY_TASK_ID + OFFSET + 2))   # +2 skips header; array is 0-based
data=$(sed -n "${LINE}p" "$PARAM_FILE" || true)

if [[ -z "${data}" ]]; then
  echo "No data at line ${LINE} in ${PARAM_FILE}. Exiting."
  exit 0
fi

IFS=',' read -r ID gmu imu gd id gdfe idfe <<< "$data"

# Demography times (integers)
T_BOT=$(python3 - <<PY
import math
END=int(${END})
print(int(round(END*float("${BOT_FRAC_TBOT}"))))
PY
)

T_REC=$(python3 - <<PY
import math
END=int(${END})
print(int(round(END*float("${BOT_FRAC_TREC}"))))
PY
)

T_EXP=$(python3 - <<PY
import math
END=int(${END})
print(int(round(END*float("${EXP_FRAC_TEXP}"))))
PY
)

# Output directories per model + END
OUT_CONST="${BASE_OUT}constant/END${END}"
OUT_BOT="${BASE_OUT}bottleneck/END${END}"
OUT_EXP="${BASE_OUT}expansion/END${END}"
mkdir -p "$OUT_CONST" "$OUT_BOT" "$OUT_EXP"

echo "JOB=$SLURM_JOB_ID TASK=$SLURM_ARRAY_TASK_ID LINE=$LINE ID=$ID END=$END"
echo "Params: gmu=$gmu imu=$imu gd=$gd id=$id gdfe=$gdfe idfe=$idfe"
echo "Demography: T_BOT=$T_BOT T_REC=$T_REC BOT_FRAC=$BOT_FRAC | T_EXP=$T_EXP EXP_MULT=$EXP_MULT"

# ----------------------------
# Run 3 models in parallel (1 CPU each)
# ----------------------------
run_const() {
  slim \
    -d "ID=${ID}" \
    -d "OUTTREES='${OUT_CONST}'" \
    -d END=${END} \
    -d SIMPLIFY_INTERVAL=${SIMPLIFY_INTERVAL} \
    -d SELF=${SELF} \
    -d gmu=${gmu} -d imu=${imu} -d gd=${gd} -d id=${id} -d gdfe=${gdfe} -d idfe=${idfe} \
    "${SCRIPT_CONST}" \
    > "${LOGDIR}/${ID}_const_END${END}.log" 2>&1
}

run_bot() {
  slim \
    -d "ID=${ID}" \
    -d "OUTTREES='${OUT_BOT}'" \
    -d END=${END} \
    -d SIMPLIFY_INTERVAL=${SIMPLIFY_INTERVAL} \
    -d SELF=${SELF} \
    -d T_BOT=${T_BOT} -d T_REC=${T_REC} -d BOT_FRAC=${BOT_FRAC} \
    -d gmu=${gmu} -d imu=${imu} -d gd=${gd} -d id=${id} -d gdfe=${gdfe} -d idfe=${idfe} \
    "${SCRIPT_BOT}" \
    > "${LOGDIR}/${ID}_bottleneck_END${END}.log" 2>&1
}

run_exp() {
  slim \
    -d "ID=${ID}" \
    -d "OUTTREES='${OUT_EXP}'" \
    -d END=${END} \
    -d SIMPLIFY_INTERVAL=${SIMPLIFY_INTERVAL} \
    -d SELF=${SELF} \
    -d T_EXP=${T_EXP} -d EXP_MULT=${EXP_MULT} \
    -d gmu=${gmu} -d imu=${imu} -d gd=${gd} -d id=${id} -d gdfe=${gdfe} -d idfe=${idfe} \
    "${SCRIPT_EXP}" \
    > "${LOGDIR}/${ID}_expansion_END${END}.log" 2>&1
}

run_const &
pid1=$!
run_bot &
pid2=$!
run_exp &
pid3=$!

# Fail the job if any model fails
wait $pid1
wait $pid2
wait $pid3

echo "Done: ID=$ID END=$END"
