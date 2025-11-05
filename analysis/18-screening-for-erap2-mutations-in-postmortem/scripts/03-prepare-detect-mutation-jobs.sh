#!/bin/bash

PROJECT_DIR="/export2/home/koustav/Projects/012-Ben-Clarke-ERAP1-KD/analysis/18-screening-for-erap2-mutations-in-postmortem"
JOB_DIR=${PROJECT_DIR}"/tmp"
REF_FA="/export2/home/koustav/Genomes/ENSEMBL_hg38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
MUTATION_LIST=${PROJECT_DIR}"/input_files/01-list-of-mutations.txt"
PATIENT_LIST=${PROJECT_DIR}"/input_files/02-list-of-patient-samples.txt"
OUT_DIR=${PROJECT_DIR}"/output_files"
IDENTITY_FILE=${PROJECT_DIR}"/id_rsa"
SCRIPT_PATH="/export2/home/koustav/Repositories/rnaseqmut/bin/rnaseqmut"
NUM_FILES=$(cat ${PATIENT_LIST} | awk '{if(NR > 1) print $9}' | wc -l)

echo "#!/usr/bin/env bash
#$ -S /bin/bash
#$ -N detect_mutations
#$ -wd ${PROJECT_DIR}
#$ -o ${JOB_DIR}/detect_mutations_\$TASK_ID.stdout
#$ -e ${JOB_DIR}/detect_mutations_\$TASK_ID.stderr
#$ -V
#$ -l h_rt=168:00:00
#$ -l mem_free=64G
#$ -q all.q
#$ -t 1-${NUM_FILES}

export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/export2/home/koustav/anaconda3/lib/
source /export2/home/koustav/anaconda3/etc/profile.d/conda.sh
declare -a BAM_FILE
declare -a CRAM_FILE
declare -a GSM
declare -a OUT_LIST" > ${JOB_DIR}/detect_mutations.sh

k=1
tail -n +2 "$PATIENT_LIST" | while IFS=$'\t' read -r -a FILE_ARRAY; do

BAM_FILE=${FILE_ARRAY[2]}
GSM=$(basename "${FILE_ARRAY[3]}")
COLUMN_ID=${FILE_ARRAY[8]}

FILENAME=${BAM_FILE%.*}

echo "CRAM_FILE[$k]=${JOB_DIR}/${GSM}/${FILENAME}.cram
BAM_FILE[$k]=${JOB_DIR}/${GSM}/${FILENAME}.bam
GSM[$k]=${GSM}
OUT_LIST[$k]=${OUT_DIR}/${COLUMN_ID}" >> ${JOB_DIR}/detect_mutations.sh
k=$((++k))

done

echo "

TASK_ID=\${SGE_TASK_ID}

CURRENT_GSM=\${GSM[\$TASK_ID]}
CURRENT_CRAM=\${CRAM_FILE[\$TASK_ID]}
CURRENT_BAM=\${BAM_FILE[\$TASK_ID]}
CURRENT_OUT=\${OUT_LIST[\$TASK_ID]}

if [[ ! -f \"${JOB_DIR}/\$CURRENT_GSM\" ]]; then 
	rsync -arv patanilab@10.246.164.99:/var/services/homes/patanilab/2025-crick-projects/proj-patanir-seq/palk/0103-NYGC-consortium/output-files/\$CURRENT_GSM ${JOB_DIR}/
fi

if [[ -f \"\${CURRENT_CRAM}\" ]]; then 
	conda run -n Samtools samtools view -b -1 --reference ${REF_FA} -o \$CURRENT_BAM \$CURRENT_CRAM && rm \$CURRENT_CRAM
fi

if [[ -f \"\${CURRENT_BAM}\" ]]; then 
	conda run -n Samtools samtools index \${CURRENT_BAM}
	${SCRIPT_PATH} -r ${REF_FA} \\
	-l ${MUTATION_LIST} \\
	-d -k -s 8 -m 1 \\
	\${CURRENT_BAM} > \${CURRENT_OUT}
fi

#if [[ -f \"\${CURRENT_OUT}\" ]]; then 
#	rm -rf ${JOB_DIR}/\$CURRENT_GSM
#fi


" >> ${JOB_DIR}/detect_mutations.sh
