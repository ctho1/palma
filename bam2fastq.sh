#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --partition normal
#SBATCH --time=4:00:00
#SBATCH --mem=20G
#SBATCH --job-name=bam2fastq
#SBATCH --mail-type=ALL
#SBATCH --output log.txt
#SBATCH --mail-user=christian.thomas@ukmuenster.de

ml palma/2019a  GCC/8.2.0-2.31.1
ml SAMtools/1.9

name=`basename $1 .bam`

samtools bam2fq -1 ./fastq/${name}_R1_001.fastq.gz -2 ./fastq/${name}_R2_001.fastq.gz --threads 24 $1