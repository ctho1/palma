#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --partition requeue
#SBATCH --time=0:30:00
#SBATCH --mem=20G
#SBATCH --job-name=sturgeon_pipeline
#SBATCH --mail-type=ALL
#SBATCH --error ./log/%x_%j.err.txt
#SBATCH --output ./log/%x_%j.out.txt
#SBATCH --mail-user=christian.thomas@ukmuenster.de

###################################################
# User Configuration
###################################################
threads=24
bowtie2db=/scratch/tmp/thomachr/references/metaphlandb


###################################################
# Run MetaPhlAn 4
###################################################
mkdir -p bowtie2out
mkdir -p output
ml palma/2022a GCC/11.3.0 OpenMPI/4.1.4 MetaPhlAn/4.0.6 Bowtie2/2.4.5

input_files=$(find ./fastq -type f -name "*R1_001.fastq.gz")

for R1 in $input_files; do

R2=`echo $R1 | sed 's/R1_001/R2_001/'`
R1_base=`basename $R1 .fastq.gz`

metaphlan $R1,$R2 \
	--bowtie2out bowtie2out/${R1_base}metagenome.bowtie2.bz2 \
	--bowtie2db $bowtie2db \
	--nproc $threads \
	--input_type fastq \
	-o output/${R1_base}_profiled_metagenome.txt
	
done