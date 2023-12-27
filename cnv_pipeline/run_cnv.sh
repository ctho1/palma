#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --partition requeue
#SBATCH --time=06:00:00
#SBATCH --mem=30G
#SBATCH --job-name=cnv_pipeline
#SBATCH --mail-type=ALL
#SBATCH --output log.txt
#SBATCH --mail-user=christian.thomas@ukmuenster.de

## Parameters ############################################################################
RefGenome=/scratch/tmp/thomachr/references/hg38/hg38.fa
FASTQ_DIR=./fastq
BAM_DIR=./bam
threads=24
mkdir -p bam
mkdir -p results

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

mkdir -p $BAM_DIR/${R1_base}
mkdir -p $BAM_DIR/${R1_base}/filtered
mkdir -p results/${R1_base}
rm -r results/${R1_base}
module unload
ml palma/2019a  GCC/8.2.0-2.31.1 BWA/0.7.17 HTSlib/1.9 SAMtools/1.9
bwa mem -t $threads $RefGenome $R1 $R2 > $BAM_DIR/${R1_base}/"$R1_base".sam

samtools view -S -b $BAM_DIR/${R1_base}/"$R1_base".sam > $BAM_DIR/${R1_base}/"$R1_base".bam
samtools sort -@ 8 -m 2G $BAM_DIR/${R1_base}/"$R1_base".bam -o $BAM_DIR/${R1_base}/"$R1_base".sorted.bam
samtools index $BAM_DIR/${R1_base}/"$R1_base".sorted.bam

rm -f $BAM_DIR/${R1_base}/"$R1_base".bam
rm -f $BAM_DIR/${R1_base}/"$R1_base".sam

# Remove all regions with coverage > 50X
module unload
ml palma/2021a GCC/10.3.0 BEDTools/2.30.0
bedtools genomecov -ibam $BAM_DIR/${R1_base}/"$R1_base".sorted.bam -bg | awk '$4 > 50' > $BAM_DIR/${R1_base}/"$R1_base".coverage.bed
#bedtools intersect -abam $BAM_DIR/${R1_base}/"$R1_base".sorted.bam  -b $BAM_DIR/${R1_base}/"$R1_base".coverage.bed -v > $BAM_DIR/${R1_base}/filtered/"$R1_base".sorted.filtered.bam

module unload
ml palma/2022a GCC/11.3.0 SAMtools/1.16.1

samtools view -@ 8 -L $BAM_DIR/${R1_base}/"$R1_base".coverage.bed -b $BAM_DIR/${R1_base}/"$R1_base".sorted.bam -U $BAM_DIR/${R1_base}/filtered/"$R1_base".sorted.filtered.bam > $BAM_DIR/${R1_base}/in_target.bam
samtools index $BAM_DIR/${R1_base}/filtered/"$R1_base".sorted.filtered.bam  
samtools index $BAM_DIR/${R1_base}/in_target.bam  

# Remove Coverage Peaks using MACS2
#module unload
#ml palma/2019a  GCC/8.2.0-2.31.1  OpenMPI/3.1.3 MACS2/2.2.6-Python-3.7.2
#macs2 callpeak -f BAMPE -t $BAM_DIR/${R1_base}/"$R1_base".sorted.bam -n ${R1_base} --outdir $BAM_DIR/${R1_base} -B

#module unload
#ml palma/2022a GCC/11.3.0 SAMtools/1.16.1
#samtools view -@ 6 -L $BAM_DIR/${R1_base}/"$R1_base"_peaks.narrowPeak -b $BAM_DIR/${R1_base}/"$R1_base".sorted.bam -U $BAM_DIR/${R1_base}/filtered/"$R1_base".sorted.filtered.bam > $BAM_DIR/${R1_base}/in_target.bam
#samtools index $BAM_DIR/${R1_base}/filtered/"$R1_base".sorted.filtered.bam  

## CNV Calling ###########################################################################

module unload
ml palma/2020b GCC/10.2.0 OpenMPI/4.0.5 R/4.0.3
Rscript cnv_plot.R $BAM_DIR/${R1_base}/filtered results/${R1_base}

done