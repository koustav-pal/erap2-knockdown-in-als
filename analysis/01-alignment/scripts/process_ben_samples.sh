#!/usr/bin/bash
#SBATCH --job-name=process_ben_erap2
#SBATCH --chdir=/nemo/lab/patanir/home/users/palk/Projects/012-Ben-Clarke-ERAP1-KD/analysis/01-alignment
#SBATCH --output=/nemo/lab/patanir/home/users/palk/Projects/012-Ben-Clarke-ERAP1-KD/analysis/01-alignment/slurm/nextflow/process_ben_erap2.stdout
#SBATCH --error=/nemo/lab/patanir/home/users/palk/Projects/012-Ben-Clarke-ERAP1-KD/analysis/01-alignment/slurm/nextflow/process_ben_erap2.stderr
#SBATCH --get-user-env
#SBATCH --mail-type=ALL
#SBATCH --mail-user=palk@crick.ac.uk
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=ncpu
#SBATCH --time=72:00:00
#SBATCH --mem=32GB


pipeline_version="3.7"

export NXF_SINGULARITY_CACHEDIR=/flask/apps/misc/stp/babs/nf-core/singularity/rnaseq/$pipeline_version/

# Specify variables for alignment

TMP_DIR="/nemo/lab/patanir/home/users/palk/Projects/012-Ben-Clarke-ERAP1-KD/analysis/01-alignment/work/"
OUT_DIR="/nemo/lab/patanir/home/users/palk/Projects/012-Ben-Clarke-ERAP1-KD/analysis/01-alignment/nfcore"

ml purge  
ml Nextflow
ml Singularity
ml CAMP_proxy

nextflow run nf-core/rnaseq \
--input /nemo/lab/patanir/home/users/palk/Projects/012-Ben-Clarke-ERAP1-KD/analysis/01-alignment/input_files/nfcore_sampleinfo.csv \
--outdir $OUT_DIR \
--email koustav.pal@crick.ac.uk \
--fasta /nemo/lab/patanir/home/users/palk/Genomes/ENSEMBL_hg38/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--gtf /nemo/lab/patanir/home/users/palk/Genomes/ENSEMBL_hg38/Homo_sapiens.GRCh38.104.chr_patch_hapl_scaff.gtf \
--aligner star_rsem \
-work-dir $TMP_DIR \
--save_trimmed \
-profile crick \
--extra_star_align_args --outSAMattributes All \
--save_reference \
-r $pipeline_version \
--max_time 72.h \
-resume