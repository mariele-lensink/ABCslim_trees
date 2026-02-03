#!/bin/bash
#SBATCH --job-name=slim_array
#SBATCH --array=40000-49999%350
#SBATCH --output=/home/mlensink/slimsimulations/ABCslim/ABC_slim/logs/logs_100k_2/slim_array_%A_%a.out
#SBATCH --error=/home/mlensink/slimsimulations/ABCslim/ABC_slim/logs/logs_100k_2/slim_array_%A_%a.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=5-00:00:00
#SBATCH --partition=bmm
#SBATCH --account=gmonroegrp
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=melensink@ucdavis.edu

source ~/anaconda3/etc/profile.d/conda.sh
conda activate slimsims

echo "SLURM_JOB_ID=$SLURM_JOB_ID"
echo "SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID"

# Default OFFSET if not passed in
: ${OFFSET:=0}

# Load the line corresponding to this task from the parameter file
PARAM_FILE="/home/mlensink/slimsimulations/ABCslim/ABC_slim/data/prior_parameters_100k_half2_9.2.2025.csv"
LINE=$(($SLURM_ARRAY_TASK_ID + $OFFSET + 2))  # Skip header and align 0-based index with data lines  
# Now that we start at line 2, task ID maps directly to the data line
data=$(sed -n "${LINE}p" $PARAM_FILE)

IFS=',' read -r ID gmu imu gd id gdfe idfe <<< "$data"


TAJIMA_OUT="/home/mlensink/slimsimulations/ABCslim/ABC_slim/data/tajima100k_2/${ID}.Tajima.D"
FINAL_VCF="/home/mlensink/slimsimulations/ABCslim/ABC_slim/data/vcf100k_2/${ID}.vcf"
FINAL_VCF_GZ="/home/mlensink/slimsimulations/ABCslim/ABC_slim/data/vcf100k_2/${ID}.vcf.gz"


# CASE 1: Tajima exists but VCF doesn't → both must be redone
if [[ ! -f "$FINAL_VCF" && ! -f "$FINAL_VCF_GZ" && -f "$TAJIMA_OUT" ]]; then
    echo "Tajima D exists but VCF is missing for ID $ID. Deleting Tajima file and rerunning both."
    rm -f "$TAJIMA_OUT"
fi

# Run SLiM if VCF missing
if [[ ! -f "$FINAL_VCF" && ! -f "$FINAL_VCF_GZ" ]]; then
    echo "VCF missing for ID $ID. Running SLiM."
    slim -d ID=$ID \
         -d JOBID=$SLURM_JOB_ID \
         -d gmu=$gmu \
         -d imu=$imu \
         -d gd=$gd \
         -d id=$id \
         -d gdfe=$gdfe \
         -d idfe=$idfe \
         /home/mlensink/slimsimulations/ABCslim/ABC_slim/scripts/ABC.slim 

    if [[ ! -f "$FINAL_VCF" ]]; then
        echo "SLiM failed to produce output for ID $ID. Logging failure."
        mkdir -p logs/failed/tmp
        echo "$ID" > logs/failed/tmp/${ID}.fail
        exit 1
    fi
else
    echo "VCF already exists for ID $ID . Skipping SLiM."
fi

if [[ -f "$FINAL_VCF" ]]; then
    echo "Using uncompressed VCF from data/vcf for ID $ID."
    VCFCMD="vcftools --vcf $FINAL_VCF"
elif [[ -f "$FINAL_VCF_GZ" ]]; then
    echo "Using compressed VCF from data/vcf (streaming) for ID $ID."
    VCFCMD="gunzip -c $FINAL_VCF_GZ | vcftools --vcf -"
else
    echo "ERROR: No VCF file found for ID $ID. Cannot run TajimaD."
    mkdir -p logs/failed/tmp
    echo "$ID" > logs/failed/tmp/${ID}.fail
    exit 1
fi

# Run TajimaD if output is missing
if [[ ! -f "$TAJIMA_OUT" ]]; then
    echo "Running Tajima D for ID $ID."
    eval "$VCFCMD --out /home/mlensink/slimsimulations/ABCslim/ABC_slim/data/tajima100k_2/${ID} --TajimaD 100"
else
    echo "Tajima D already exists for ID $ID. Skipping."
fi


# Compress VCF and copy to permanent storage
if [[ -f "$FINAL_VCF" ]]; then
    echo "Compressing and storing VCF for ID $ID."
    bgzip -c "$FINAL_VCF">"$FINAL_VCF_GZ"
fi

# If last task in array, aggregate failed IDs and clean up
if [[ $SLURM_ARRAY_TASK_ID -eq 49999 ]]; then
    echo "Aggregating failure logs..."
    mkdir -p logs/failed
    cat logs/failed/tmp/*.fail > /home/mlensink/slimsimulations/ABCslim/ABC_slim/logs/failed/slim_failed_ids.txt
    rm -f logs/failed/tmp/*.fail
    rmdir logs/failed/tmp
fi
