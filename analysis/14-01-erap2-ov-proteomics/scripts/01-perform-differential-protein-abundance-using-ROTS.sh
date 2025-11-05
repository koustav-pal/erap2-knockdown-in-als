#!/usr/bin/env bash
#$ -S /bin/bash
#$ -N ov_ROTS
#$ -wd /export2/home/koustav/Projects/012-Ben-Clarke-ERAP1-KD/analysis/14-01-erap2-ov-proteomics/
#$ -o /export2/home/koustav/Projects/012-Ben-Clarke-ERAP1-KD/analysis/14-01-erap2-ov-proteomics/tmp/erap2_ov_ROTS.stdout
#$ -e /export2/home/koustav/Projects/012-Ben-Clarke-ERAP1-KD/analysis/14-01-erap2-ov-proteomics/tmp/erap2_ov_ROTS.stderr
#$ -V
#$ -l h_rt=168:00:00
#$ -l mem_free=64G
#$ -q all.q

source /export2/home/koustav/anaconda3/etc/profile.d/conda.sh

conda activate r_4.4

Rscript scripts/01-perform-differential-protein-abundance-using-ROTS.R