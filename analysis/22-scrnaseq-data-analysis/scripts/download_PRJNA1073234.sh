source /home/koustav/anaconda3/etc/profile.d/conda.sh
conda activate kingfisher

mkdir -p /data/shared/public_datasets/scrnaseq-als-ftld-motor-pfc-cortex-pineda-et-al-2024/fastq
cd /data/shared/public_datasets/scrnaseq-als-ftld-motor-pfc-cortex-pineda-et-al-2024/fastq
kingfisher get -p PRJNA1073234 -m aws-http -f fastq.gz --include-technical --read-layout 10x-gex --read-name-style illumina
mkdir ../input_files
kingfisher annotate -p PRJNA1073234 -o ../input_files/sampleinfo.csv -f csv -a
