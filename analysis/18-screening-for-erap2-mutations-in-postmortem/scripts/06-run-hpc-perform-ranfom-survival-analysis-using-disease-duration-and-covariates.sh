#!/usr/bin/env bash
#$ -S /bin/bash
#$ -N detect_mutations
#$ -wd /export2/home/koustav/Projects/012-Ben-Clarke-ERAP1-KD/analysis/18-screening-for-erap2-mutations-in-postmortem
#$ -o /export2/home/koustav/Projects/012-Ben-Clarke-ERAP1-KD/analysis/18-screening-for-erap2-mutations-in-postmortem/tmp/survival_analysis.stdout
#$ -e /export2/home/koustav/Projects/012-Ben-Clarke-ERAP1-KD/analysis/18-screening-for-erap2-mutations-in-postmortem/tmp/survival_analysis.stderr
#$ -V
#$ -l h_rt=168:00:00
#$ -l mem_free=250G
#$ -q all.q

source /export2/home/koustav/anaconda3/etc/profile.d/conda.sh

conda activate r_4.4

Rscript scripts/06-perform-ranfom-survival-analysis-using-disease-duration-and-covariates.R
