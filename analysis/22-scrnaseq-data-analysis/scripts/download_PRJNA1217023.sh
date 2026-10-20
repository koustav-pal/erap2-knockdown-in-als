source /home/koustav/anaconda3/etc/profile.d/conda.sh
conda activate kingfisher

mkdir -p /data/shared/public_datasets/spatial-als-spinal-cord-zhang-kiskinis-gate-2026/fastq
cd /data/shared/public_datasets/spatial-als-spinal-cord-zhang-kiskinis-gate-2026/fastq
kingfisher get -p PRJNA1217023 -m aws-http -f fastq.gz --include-technical --read-layout 10x-gex --read-name-style illumina
mkdir ../input_files
kingfisher annotate -p PRJNA1217023 -o ../input_files/sampleinfo.csv -f csv -a


mkdir -p /data/shared/public_datasets/scrnaseq-als-spinal-cord-zhang-kiskinis-gate-2026/fastq
cd /data/shared/public_datasets/scrnaseq-als-spinal-cord-zhang-kiskinis-gate-2026/fastq
kingfisher get -p PRJNA1119910 -m aws-http -f fastq.gz --include-technical --read-layout 10x-gex --read-name-style illumina
mkdir ../input_files
kingfisher annotate -p PRJNA1119910 -o ../input_files/sampleinfo.csv -f csv -a