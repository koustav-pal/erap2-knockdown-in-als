source /home/koustav/anaconda3/etc/profile.d/conda.sh
conda activate kingfisher

mkdir -p /data/shared/public_datasets/scrnaseq-als-spinal-cord-zelic-2025/fastq
cd /data/shared/public_datasets/scrnaseq-als-spinal-cord-zelic-2025/fastq
kingfisher get -p PRJNA1211983 -m aws-http -f fastq.gz --include-technical --read-layout 10x-gex --read-name-style illumina
mkdir ../input_files
kingfisher annotate -p PRJNA1211983 -o ../input_files/sampleinfo.csv -f csv -a