#!/bin/bash

PROJECT_DIR="/camp/home/palk/home/users/palk/Projects/012-Ben-Clarke-ERAP1-KD/analysis/04-locating-mutations-in-the-cell-line"
SLURM_DIR=${PROJECT_DIR}/slurm/
FASTQ_DIR="/camp/lab/patanir/data/STPs/babs/inputs/ben.clarke/asf/RN22251/220916_A01366_0279_AHN3GFDSX3/fastq/"
FILE_LIST=$(awk 'BEGIN{FS="/";} {if(NR==FNR){split($NF, arr, "."); gsub("_R1", "", arr[1]); A[arr[1]]=$0; next;}}{split($NF,arr, "."); gsub("_R2", "", arr[1]); if(arr[1] in A){OFS = " "; print arr[1],$0,A[arr[1]]}}' <(ls ${FASTQ_DIR}/*_R1_001.fastq.gz) <(ls ${FASTQ_DIR}/*_R2_001.fastq.gz))
FILE_LIST_LEN=$(awk 'BEGIN{FS="/";} {if(NR==FNR){split($NF, arr, "."); gsub("_R1", "", arr[1]); A[arr[1]]=$0; next;}}{split($NF,arr, "."); gsub("_R2", "", arr[1]); if(arr[1] in A){OFS = " "; print arr[1],$0,A[arr[1]]}}' <(ls ${FASTQ_DIR}/*_R1_001.fastq.gz) <(ls ${FASTQ_DIR}/*_R2_001.fastq.gz) | wc -l)
output_path="/camp/home/palk/home/users/palk/Projects/012-Ben-Clarke-ERAP1-KD/analysis/04-locating-mutations-in-the-cell-line/output_files"
mutscan_csv_mutations="/camp/home/palk/home/users/palk/Genomes/als_mutation_list_mutscan.csv"
ANALYSIS_DIR=${PROJECT_DIR}

echo "#!/bin/bash
#SBATCH --job-name=run_mutscan
#SBATCH --chdir=${PROJECT_DIR}
#SBATCH --get-user-env
#SBATCH --mail-type=ALL
#SBATCH --mail-user=palk@crick.ac.uk
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu
#SBATCH --time=72:00:00
#SBATCH --mem-per-cpu=8GB
#SBATCH --array=1-${FILE_LIST_LEN}%8
#SBATCH --output=${SLURM_DIR}/run_mutscan.%a.stdout
#SBATCH --error=${SLURM_DIR}/run_mutscan.%a.stderr

source /camp/home/palk/home/users/palk/anaconda3/etc/profile.d/conda.sh
conda activate MutScan

declare -A READ_1_LIST
declare -A READ_2_LIST
declare -A OUT_LIST" > ${SLURM_DIR}/run_mutscan.sh

k=1
while read line; do


FILE_ARRAY=(${line})

FILENAME=${FILE_ARRAY[0]}
READ_1=${FILE_ARRAY[1]}
READ_2=${FILE_ARRAY[2]}

echo "
READ_1_LIST[$k]=$READ_1
READ_2_LIST[$k]=$READ_2
OUT_LIST[$k]=$FILENAME" >> ${SLURM_DIR}/run_mutscan.sh

k=$((++k))

done <<< "$FILE_LIST"

echo "

CURRENT_JOB=\${SLURM_ARRAY_TASK_ID}

CURRENT_READ_1=\${READ_1_LIST[\${CURRENT_JOB}]}
CURRENT_READ_2=\${READ_2_LIST[\${CURRENT_JOB}]}
CURRENT_OUT=\${OUT_LIST[\${CURRENT_JOB}]}

mutscan -1 \$CURRENT_READ_1 -2 \$CURRENT_READ_2 -m $mutscan_csv_mutations -h ${output_path}/\${CURRENT_OUT}.html
" >> ${SLURM_DIR}/run_mutscan.sh

sbatch ${SLURM_DIR}/run_mutscan.sh
