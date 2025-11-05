#!/usr/bin/env bash
#$ -S /bin/bash
#$ -N detect_gwas_mutations
#$ -wd /export2/home/koustav/Projects/012-Ben-Clarke-ERAP1-KD/analysis/18-screening-for-erap2-mutations-in-postmortem
#$ -o /export2/home/koustav/Projects/012-Ben-Clarke-ERAP1-KD/analysis/18-screening-for-erap2-mutations-in-postmortem/tmp/detect_gwas_mutations.stdout
#$ -e /export2/home/koustav/Projects/012-Ben-Clarke-ERAP1-KD/analysis/18-screening-for-erap2-mutations-in-postmortem/tmp/detect_gwas_mutations.stderr
#$ -V
#$ -l h_rt=168:00:00
#$ -l mem_free=250G
#$ -q all.q

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/export2/home/koustav/anaconda3/lib/
source /export2/home/koustav/anaconda3/etc/profile.d/conda.sh

# cat input_files/02-list-of-patient-samples.txt | awk -F"\t" '{if(NR == 1){next;} print $9}' | wc -l

ALL_FILES=$(cat input_files/02-list-of-patient-samples.txt | awk -F"\t" '{if(NR == 1){next;} if(NR < 1348){printf "%s ", "gwas_scan_output_files/"$9;} else{printf "%s", "gwas_scan_output_files/"$9} }')

# cat input_files/02-list-of-patient-samples.txt | awk -F"\t" 'BEGIN{a=1; b=1; c=1; d=1}{ if(NR==1){next;} if($23 == "ctrl"){j=a; a++;} else if($23 == "als"){ j = b; b++; } else if($23 == "ftd"){ j = c; c++; } else if($23 == "als+ftd"){ j = d; d++; } print $23"_"j;}' | wc -l

LABELS=$(cat input_files/02-list-of-patient-samples.txt | awk -F"\t" 'BEGIN{a=1;}{ if(NR==1){next;} if(NR < 1348){printf "%s,", $23"_"a;} else{printf "%s", $23"_"a;} a++;}')
CONTROL_GROUPS=$(cat input_files/02-list-of-patient-samples.txt | awk -F"\t" 'BEGIN{a=0;}{ if(NR==1){next;} if($23 == "ctrl"){j=a; if(NR < 1338){printf "%s,",j;} else{printf "%s",j;}} a++;}')

python3 /export2/home/koustav/Repositories/rnaseqmut/script/merge2ndvcf.py -l $LABELS $ALL_FILES > tmp/GWAS_MUTATION_SCREENING_RESULTS.txt

python3 /export2/home/koustav/Repositories/rnaseqmut/script/filtermut.py -c ${CONTROL_GROUPS} -f 0.0 -d 10 -b 0 -l ${LABELS} < tmp/GWAS_MUTATION_SCREENING_RESULTS.txt > gwas_scan_output_files/GWAS_FILTERED_MUTATION_SCREENING_RESULTS.vcf