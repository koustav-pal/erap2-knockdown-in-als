cd /home/koustav/Projects/012-Ben-Clarke-ERAP1-KD/analysis/22-scrnaseq-data-analysis/

pipeline_version="4.2.0"
TMP_DIR="/data/shared/public_datasets/scrnaseq-als-spinal-cord-zelic-2025/work"
OUT_DIR="/data/shared/public_datasets/scrnaseq-als-spinal-cord-zelic-2025/nfcore"

# export NXF_OFFLINE='true'

nextflow run nf-core/scrnaseq \
--input ./input_files/zelic-2025-nfcore-sampleinfo.csv \
--outdir $OUT_DIR \
--fasta /data/shared/public_datasets/Genomes/ENSEMBL_hg38/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--gtf /data/shared/public_datasets/Genomes/ENSEMBL_hg38/Homo_sapiens.GRCh38.104.chr_patch_hapl_scaff.gtf \
--aligner star \
--protocol 10XV3 \
-work-dir $TMP_DIR \
-profile sge,apptainer \
--save_reference \
-c ./input_files/star.config \
-r $pipeline_version \
-resume
