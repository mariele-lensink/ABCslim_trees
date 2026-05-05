#!/bin/bash
#SBATCH --job-name=demes_exp
#SBATCH --output=/home/mlensink/slimsimulations/ABCslim_trees/logs/demes_expansion/%x_%A_%a.out
#SBATCH --error=/home/mlensink/slimsimulations/ABCslim_trees/logs/demes_expansion/%x_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=16G
#SBATCH --time=3-00:00:00
#SBATCH --partition=bml
#SBATCH --account=gmonroegrp
#SBATCH --array=0-49999
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=melensink@ucdavis.edu

set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh
conda activate abc_trees

: "${PARAM_FILE:=/home/mlensink/slimsimulations/ABCslim_trees/priors/priors_2.6.25.csv}"
: "${OFFSET:=0}"
: "${END:=6000}"
: "${SIMPLIFY_INTERVAL:=200}"
: "${SELF:=0.98}"
: "${DEMES_EXP_YAML:=/home/mlensink/slimsimulations/ABCslim_trees/demes/expansion.yaml}"
: "${BASE_OUT:=/home/mlensink/slimsimulations/ABCslim_trees/output}"
: "${CHROM_DIR:=chrom1}"
: "${CONST_TREES_DIR:=${BASE_OUT}/constant/${CHROM_DIR}/trees}"
: "${SCRIPT_EXP:=/home/mlensink/slimsimulations/ABCslim_trees/models/ABCtrees_expansion_deme.slim}"

OUT_EXP=/home/mlensink/slimsimulations/ABCslim_trees/output/expansion/trees
LOGDIR=/home/mlensink/slimsimulations/ABCslim_trees/logs/demes_expansion
mkdir -p "$OUT_EXP" "$LOGDIR"

LINE=$((SLURM_ARRAY_TASK_ID + OFFSET + 2))
data=$(sed -n "${LINE}p" "$PARAM_FILE" || true)

if [[ -z "${data}" ]]; then
  echo "No data at line ${LINE} in ${PARAM_FILE}. Exiting."
  exit 0
fi

IFS=',' read -r ID gmu imu gd id gdfe idfe <<< "$data"

CONST_TREE="${CONST_TREES_DIR}/${ID}.trees"
if [[ ! -s "$CONST_TREE" ]]; then
    echo "Missing constant tree file: $CONST_TREE" >&2
    exit 1
fi

if [[ ! -f "${DEMES_EXP_YAML}" ]]; then
    echo "Missing expansion demes YAML: ${DEMES_EXP_YAML}" >&2
    exit 1
fi

mapfile -t DEMES_INFO < <(python3 - "${DEMES_EXP_YAML}" "${END}" <<'PY'
import sys
import yaml

yaml_path = sys.argv[1]
end = int(sys.argv[2])

with open(yaml_path, "r", encoding="utf-8") as fh:
    graph = yaml.safe_load(fh)

demes_list = graph.get("demes", [])
if len(demes_list) != 1:
    raise SystemExit("Expected exactly one deme in YAML.")

epochs = demes_list[0].get("epochs", [])
if not epochs:
    raise SystemExit("YAML must define at least one epoch.")

initial_size = int(round(float(epochs[0]["start_size"])))
event_times = []
event_sizes = []
summary_parts = [f"initial_size={initial_size}"]

for idx, epoch in enumerate(epochs[1:], start=1):
    if "start_time" not in epoch:
        raise SystemExit(f"Epoch {idx} is missing start_time.")
    if "start_size" not in epoch:
        raise SystemExit(f"Epoch {idx} is missing start_size.")

    start_time = float(epoch["start_time"])
    size = int(round(float(epoch["start_size"])))
    forward_time = int(round(end - start_time))

    if forward_time < 1 or forward_time > end:
        raise SystemExit(
            f"Epoch {idx} start_time={start_time} maps to invalid forward generation {forward_time} for END={end}."
        )

    event_times.append(forward_time)
    event_sizes.append(size)
    summary_parts.append(f"gen={forward_time}->N={size}")

if any(event_times[i] <= event_times[i - 1] for i in range(1, len(event_times))):
    raise SystemExit("Epoch start_times must map to strictly increasing forward event times.")

times_expr = "c(" + ",".join(str(x) for x in event_times) + ")" if event_times else "c()"
sizes_expr = "c(" + ",".join(str(x) for x in event_sizes) + ")" if event_sizes else "c()"

print(initial_size)
print(times_expr)
print(sizes_expr)
print("; ".join(summary_parts))
PY
)

INITIAL_SIZE="${DEMES_INFO[0]}"
EVENT_TIMES_EXPR="${DEMES_INFO[1]}"
EVENT_SIZES_EXPR="${DEMES_INFO[2]}"
DEMES_SUMMARY="${DEMES_INFO[3]}"

cat <<EOF
JOB=$SLURM_JOB_ID TASK=$SLURM_ARRAY_TASK_ID LINE=$LINE ID=$ID END=$END MODEL=expansion
Params: gmu=$gmu imu=$imu gd=$gd id=$id gdfe=$gdfe idfe=$idfe
Using constant tree: $CONST_TREE
Demes YAML: $DEMES_EXP_YAML
Demography schedule: $DEMES_SUMMARY
EOF

slim \
  -d "ID='${ID}'" \
  -d "OUTTREES='${OUT_EXP}'" \
  -d "INITTREES='${CONST_TREE}'" \
  -d "INITIAL_SIZE=${INITIAL_SIZE}" \
  -d "EVENT_TIMES=${EVENT_TIMES_EXPR}" \
  -d "EVENT_SIZES=${EVENT_SIZES_EXPR}" \
  -d END=${END} -d SIMPLIFY_INTERVAL=${SIMPLIFY_INTERVAL} -d SELF=${SELF} \
  -d gmu=${gmu} -d imu=${imu} -d gd=${gd} -d id=${id} -d gdfe=${gdfe} -d idfe=${idfe} \
  "${SCRIPT_EXP}" \
  > "${LOGDIR}/${ID}_expansion.log" 2>&1

echo "Done: ID=$ID END=$END MODEL=expansion"
