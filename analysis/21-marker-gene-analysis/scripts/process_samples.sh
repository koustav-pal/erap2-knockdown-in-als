cd /data/shared/lab_internal/012-Ben-Clarke-ERAP1-KD/analysis/21-marker-gene-analysis
export NXF_HOME="/home/koustav/.nextflow"
# source /home/svu/koustav.pal/miniconda3/etc/profile.d/conda.sh
# source activate base
# conda activate nextflow
 
pipeline_version="1.0.0"
TMP_DIR="./work"
OUT_DIR="./nfcore"

# all_files=$(ls /home/koustav/Projects/100-05-ben-clarke-SFPQKD-in-astrocytes-collaboration/analysis/public_data/purified-cortex-human-zhang-2016/nfcore/star_salmon/*.bam)
# source /home/koustav/anaconda3/etc/profile.d/conda.sh
# conda run -n samtools samtools index --threads 8 --bai -M $all_files


nextflow run koustav-pal/nf-featurequant \
-hub github \
-r main \
-latest \
--input ./input_files/nfcore_sampleinfo.csv \
--annotation_rds /data/shared/public_datasets/Genomes/ENSEMBL_hg38/featurequant_ref/feature_ranges.rds \
--organism human \
--aligner star \
--outdir $OUT_DIR \
-work-dir $TMP_DIR \
-profile sge,apptainer \
--pipeline_version $pipeline_version \
--max_time 72.h \
-resume