#!/usr/bin/env bash
#$ -S /bin/bash
#$ -N compute_gene_counts
#$ -wd /export2/home/koustav/Projects/012-Ben-Clarke-ERAP1-KD/analysis/11-postmortem-data-analysis/
#$ -o /export2/home/koustav/Projects/012-Ben-Clarke-ERAP1-KD/analysis/11-postmortem-data-analysis/tmp/compute_gene_counts.stdout
#$ -e /export2/home/koustav/Projects/012-Ben-Clarke-ERAP1-KD/analysis/11-postmortem-data-analysis/tmp/compute_gene_counts.stderr
#$ -V
#$ -l h_rt=168:00:00
#$ -l mem_free=250G
#$ -q all.q

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/export2/home/koustav/anaconda3/lib/
source /export2/home/koustav/anaconda3/etc/profile.d/conda.sh

conda activate r_4.4

Rscript scripts/01-save-gene-level-counts.R