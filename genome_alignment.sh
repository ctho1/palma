#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --partition normal
#SBATCH --time=00:30:00
#SBATCH --mem=50G
#SBATCH --job-name=alignment
#SBATCH --mail-type=ALL
#SBATCH --output log.txt
#SBATCH --mail-user=christian.thomas@ukmuenster.de

## Parameters ############################################################################
RefGenome=/scratch/tmp/thomachr/references/genomes/tbev/sequence.fasta
FASTQ_DIR=./fastq
BAM_DIR=./bam
threads=12

if [ -f "$RefGenome".amb ]; then
    echo "Reference genome index exists."
else 
    echo "Reference genome index does not exist. Creating index with bwa index"
    bwa index $RefGenome
fi


## Alignment #############################################################################
in_files=$(find $FASTQ_DIR -type f -name "*R1*.fastq.gz" -print|sort)

for R1 in $in_files; do
        R2=`echo $R1 | sed 's/R1_001/R2_001/'`
        R1_base=`basename $R1 .fastq.gz`

module unload
ml palma/2019a  GCC/8.2.0-2.31.1 BWA/0.7.17 HTSlib/1.9 SAMtools/1.9
bwa mem -t $threads $RefGenome $R1 $R2 > $BAM_DIR/"$R1_base".sam

samtools view -S -b $BAM_DIR/"$R1_base".sam > $BAM_DIR/"$R1_base".bam
samtools sort -@ 8 -m 2G $BAM_DIR/"$R1_base".bam -o $BAM_DIR/"$R1_base".sorted.bam
samtools index $BAM_DIR/"$R1_base".sorted.bam

rm -f $BAM_DIR/"$R1_base".bam
rm -f $BAM_DIR/"$R1_base".sam

## Coverage Plot ##########################################################################

module unload
ml palma/2020a  GCC/9.3.0  OpenMPI/4.0.3  palma/2020a  iccifort/2020.1.217  impi/2019.7.217  matplotlib/3.2.1-Python-3.8.2
tinycov covplot --skip 100 --res 500 --out $BAM_DIR/"$R1_base".coverage.jpg --text $BAM_DIR/"$R1_base".coverage.txt $BAM_DIR/"$R1_base".sorted.bam

done
