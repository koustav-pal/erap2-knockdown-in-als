source /home/koustav/anaconda3/etc/profile.d/conda.sh
conda activate kingfisher

mkdir -p /data/shared/public_datasets/scrnaseq-motor-cortex-c9-als-li-et-al-2023/fastq
mkdir -p /data/shared/public_datasets/scrnaseq-motor-cortex-c9-als-li-et-al-2023/input_files

cp /home/koustav/Projects/012-Ben-Clarke-ERAP1-KD/analysis/22-scrnaseq-data-analysis/input_files/sample_ids.csv \
/data/shared/public_datasets/scrnaseq-motor-cortex-c9-als-li-et-al-2023/input_files/

cd /data/shared/public_datasets/scrnaseq-motor-cortex-c9-als-li-et-al-2023/fastq

kingfisher get --run-identifiers-list ../input_files/sampleinfo.csv -m aws-http -f fastq.gz --include-technical --read-layout 10x-gex --read-name-style illumina

kingfisher annotate -p PRJNA907856 -o ../input_files/sampleinfo.csv -f csv -a