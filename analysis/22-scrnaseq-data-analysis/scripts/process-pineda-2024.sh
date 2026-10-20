cd /home/koustav/Projects/012-Ben-Clarke-ERAP1-KD/analysis/22-scrnaseq-data-analysis/

pipeline_version="4.2.0"
TMP_DIR="/data/shared/public_datasets/scrnaseq-als-ftld-motor-pfc-cortex-pineda-et-al-2024/work"
OUT_DIR="/data/shared/public_datasets/scrnaseq-als-ftld-motor-pfc-cortex-pineda-et-al-2024/nfcore"

# export NXF_OFFLINE='true'
# Human
nextflow run nf-core/scrnaseq \
--input ./input_files/pineda-2024-nfcore-sampleinfo.csv \
--outdir $OUT_DIR \
--star_index /data/shared/public_datasets/Genomes/ENSEMBL_hg38/star_index/ \
--aligner star \
--protocol 10XV3 \
-work-dir $TMP_DIR \
-profile sge,apptainer \
--save_reference \
-c ./input_files/star.config \
-r $pipeline_version \
-resume


# # Mouse
# nextflow run nf-core/scrnaseq \
# --input ./input_files/pineda-2024-mouse-nfcore-sampleinfo.csv \
# --outdir $OUT_DIR \
# --fasta /data/shared/public_datasets/Genomes/ENSEMBL_mm10/Mus_musculus.GRCm39.dna.primary_assembly.fa \
# --gtf /data/shared/public_datasets/Genomes/ENSEMBL_mm10/Mus_musculus.GRCm39.104.chr.gtf \
# --aligner star \
# --protocol 10XV3 \
# -work-dir $TMP_DIR \
# -profile sge,apptainer \
# --save_reference \
# -c ./input_files/star.config \
# -r $pipeline_version \
# -resume