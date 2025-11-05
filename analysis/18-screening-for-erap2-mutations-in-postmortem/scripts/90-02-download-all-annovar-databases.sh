#!/usr/bin/env bash
#$ -S /bin/bash
#$ -N run_annovar
#$ -wd /export2/home/koustav/Projects/012-Ben-Clarke-ERAP1-KD/analysis/18-screening-for-erap2-mutations-in-postmortem
#$ -o /export2/home/koustav/Projects/012-Ben-Clarke-ERAP1-KD/analysis/18-screening-for-erap2-mutations-in-postmortem/tmp/run_annovar.stdout
#$ -e /export2/home/koustav/Projects/012-Ben-Clarke-ERAP1-KD/analysis/18-screening-for-erap2-mutations-in-postmortem/tmp/run_annovar.stderr
#$ -V
#$ -l h_rt=168:00:00
#$ -l mem_free=250G
#$ -q all.q

source /export2/home/koustav/anaconda3/etc/profile.d/conda.sh

conda activate perl

perl /export2/home/koustav/packages/annovar/annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene databases/humandb/
perl /export2/home/koustav/packages/annovar/annotate_variation.pl -buildver hg38 -downdb -webfrom annovar gnomad41_exome databases/humandb/
perl /export2/home/koustav/packages/annovar/annotate_variation.pl -buildver hg38 -downdb -webfrom annovar clinvar_20250721 databases/humandb/
perl /export2/home/koustav/packages/annovar/annotate_variation.pl -buildver hg38 -downdb -webfrom annovar dbnsfp47a databases/humandb/
perl /export2/home/koustav/packages/annovar/annotate_variation.pl -buildver hg38 -downdb -webfrom annovar logofunc databases/humandb/
perl /export2/home/koustav/packages/annovar/annotate_variation.pl -buildver hg38 -downdb -webfrom annovar polyphen2 databases/humandb/