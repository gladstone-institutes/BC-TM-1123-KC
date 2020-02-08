#!/bin/sh
#PBS -l walltime=700:00:00
#PBS -l mem=200gb
set -e           # <-- abort if ANY command returns non-zero value
set -u           # <-- abort on undefined variables
set -o pipefail  # <-- show failed exit codes properly

crangerDir=/data/home/kchoudha/softwares/cellranger-3.1.0
refs=${crangerDir}/References/refdata-cellranger-GRCh38-3.0.0
fqDir=/data/projects/TM-1123/rawfiles/outs/fastq_path/HH7J2DRXX
outDir=/data/projects/TM-1123/analysis_120619/20.cellranger_count

#for i in {1}
#do
i=1
${crangerDir}/cellranger count --id=S$i \
	--fastqs=${fqDir}/$i \
	--sample=$i \
	--transcriptome=${refs} \
	--localcores=16 \
	--localmem=200
#done

