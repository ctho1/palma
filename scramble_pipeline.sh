#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --partition normal
#SBATCH --time=02:00:00
#SBATCH --mem=30G
#SBATCH --job-name=scramble_pipeline
#SBATCH --mail-type=ALL
#SBATCH --output ./log_scramble_pipeline.txt
#SBATCH --mail-user=ccthomas10@gmail.com

## Parameters ############################################################################
RefGenome=/scratch/tmp/thomachr/references/GRCh37/GRCh37.fa
SCRAMBLE_DIR=/scratch/tmp/thomachr/software/scramble
OUT_DIR=/scratch/tmp/thomachr/scramble/out
FASTQ_DIR=/scratch/tmp/thomachr/scramble/fastq
BAM_DIR=/scratch/tmp/thomachr/scramble/bam
threads=24

## Alignment #############################################################################
in_files=$(find $FASTQ_DIR -type f -name "*R1_001.fastq.gz" -print|sort)

for R1 in $in_files; do
        R2=`echo $R1 | sed 's/R1_001/R2_001/'`
        R1_base=`basename $R1`

module unload
ml palma/2019a  GCC/8.2.0-2.31.1 BWA/0.7.17 HTSlib/1.9 SAMtools/1.9
bwa mem -t $threads $RefGenome $R1 $R2 > $BAM_DIR/"$R1_base".sam

samtools view -S -b $BAM_DIR/"$R1_base".sam > $BAM_DIR/"$R1_base".bam
samtools sort -@ 8 -m 2G $BAM_DIR/"$R1_base".bam -o $BAM_DIR/"$R1_base".sorted.bam
samtools index $BAM_DIR/"$R1_base".sorted.bam

rm -f $BAM_DIR/"$R1_base".bam
rm -f $BAM_DIR/"$R1_base".sam

#done

## Scramble ##############################################################################
module unload
ml palma/2020b  GCC/10.2.0  OpenMPI/4.0.5 R-bundle-Bioconductor/3.12-R-4.0.3 BLAST+/2.11.0

#in_files=$(find $BAM_DIR -type f -name "*.sorted.bam" -print|sort)
#for BAM in $in_files; do

#NAME=`basename $BAM`

$SCRAMBLE_DIR/cluster_identifier/src/build/cluster_identifier $BAM_DIR/"$R1_base".sorted.bam > $OUT_DIR/"$R1_base"_clusters.txt

Rscript --vanilla $SCRAMBLE_DIR/cluster_analysis/bin/SCRAMble.R \
	--out-name $OUT_DIR/"$R1_base" \
	--cluster-file $OUT_DIR/"$R1_base"_clusters.txt \
	--install-dir $SCRAMBLE_DIR/cluster_analysis/bin \
	--mei-refs $SCRAMBLE_DIR/cluster_analysis/resources/MEI_consensus_seqs.fa \
	--pct-align 20 \
	--ref $RefGenome \
	--eval-meis \
	--eval-dels 

#done

## VCF Anno ##############################################################################
module unload
ml palma/2019a  GCC/8.2.0-2.31.1  OpenMPI/3.1.3 VEP/102.0-Perl-5.28.1

#in_files=$(find $OUT_DIR -type f -name "*.vcf" -print|sort)
#for VCF in $in_files; do
#NAME=`basename $VCF`
vep -i $OUT_DIR/"$R1_base".vcf -o ./vcf/"$R1_base".vep.vcf --database --species "human" --assembly GRCh37 --symbol --force --tab

done


##########################################################################################
exit 0
