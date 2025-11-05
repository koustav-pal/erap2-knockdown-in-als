#!/bin/bash

PROJECT_DIR="/camp/home/palk/home/users/palk/Projects/012-Ben-Clarke-ERAP1-KD/analysis/04-locating-mutations-in-the-cell-line"
SLURM_DIR=${PROJECT_DIR}/slurm/
cd $PROJECT_DIR
FILE_LIST=$(ls ../01-alignment/nfcore/star_rsem/*.markdup.sorted.bam)
FILE_LIST_LEN=$(ls ../01-alignment/nfcore/star_rsem/*.markdup.sorted.bam | wc -l)
GENOME="/camp/home/palk/home/users/palk/Genomes/ENSEMBL_hg38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GTF="/camp/home/palk/home/users/palk/Genomes/ENSEMBL_hg38/Homo_sapiens.GRCh38.104.chr_patch_hapl_scaff.gtf"
ANALYSIS_DIR=${PROJECT_DIR}
IRFINDER_RESULTS_DIR=$ANALYSIS_DIR/irfinder_results
IRFINDER_REF_DIR="../../../0002-Astrocyte-LifeArc-intron-retention/analysis/03-intron-retention-analysis/irfinder_results/reference/"


echo "#!/bin/bash
#SBATCH --job-name=run_rnaseqmut
#SBATCH --chdir=${PROJECT_DIR}
#SBATCH --get-user-env
#SBATCH --mail-type=ALL
#SBATCH --mail-user=palk@crick.ac.uk
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu
#SBATCH --time=72:00:00
#SBATCH --mem-per-cpu=8GB
#SBATCH --array=1-${FILE_LIST_LEN}%8
#SBATCH --output=${SLURM_DIR}/run_ir_finder.%a.stdout
#SBATCH --error=${SLURM_DIR}/run_ir_finder.%a.stderr


export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:~/home/users/palk/anaconda3/lib/

declare -A BAM_LIST
declare -A OUT_LIST" > ${SLURM_DIR}/run_rnaseqmut.sh

k=1
while read line; do

FILE_ARRAY=(${line})
FILENAME=${FILE_ARRAY[0]}

FILEBASE=$(dirname "$FILENAME")
DIR_BASE=$(basename $(dirname "$FILEBASE"))
filename=$(basename -- "$FILENAME")


OUTDIR=${IRFINDER_RESULTS_DIR}/$DIR_BASE/$filename

echo "
BAM_LIST[$k]=$FILENAME
OUT_LIST[$k]=$filename" >> ${SLURM_DIR}/run_rnaseqmut.sh

k=$((++k))

done <<< "$FILE_LIST"

echo "

CURRENT_JOB=\${SLURM_ARRAY_TASK_ID}

CURRENT_BAM=\${BAM_LIST[\${CURRENT_JOB}]}
CURRENT_OUT=\${OUT_LIST[\${CURRENT_JOB}]}

../../../../GitHub/rnaseqmut/bin/rnaseqmut -r ${GENOME} -l input_files/mutation_list -d -k -s 8 -m 1 \${CURRENT_BAM} > output_files/\${CURRENT_OUT}
" >> ${SLURM_DIR}/run_rnaseqmut.sh

sbatch ${SLURM_DIR}/run_rnaseqmut.sh
