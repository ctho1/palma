#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=36
#SBATCH --partition requeue
#SBATCH --time=4:00:00
#SBATCH --mem=140G
#SBATCH --job-name=sturgeon_pipeline
#SBATCH --mail-type=ALL
#SBATCH --error ./logs/%x_%j.err.txt
#SBATCH --output ./logs/%x_%j.out.txt
#SBATCH --mail-user=christian.thomas@ukmuenster.de

mkdir -p dorado_output
mkdir -p modkit_output
mkdir -p results

sturgeon_model=/scratch/tmp/thomachr/software/sturgeon/include/models/general.zip
dorado_model=/scratch/tmp/thomachr/software/dorado-0.4.2-linux-x64/models/dna_r10.4.1_e8.2_400bps_hac@v4.2.0

# search for dirs
dirs=$(find ./pod5 -mindepth 1 -maxdepth 1 -type d -print|sort)

for dir in $dirs; do

base=`basename $dir`
mkdir -p tmp/${base}
mkdir -p tmp/${base}/hg38
mkdir -p tmp/${base}/chm13v2
mkdir -p results/${base}
mkdir -p results/${base}/cnv

# Dorado Basecalling
dorado basecaller --modified-bases 5mCG_5hmCG -x cpu \
	$dorado_model $dir > tmp/${base}/${base}_calls.bam

# Dorado Alignment chm13v2
dorado aligner /scratch/tmp/thomachr/references/T2T/chm13v2.0.fa -t 36 \
	tmp/${base}/${base}_calls.bam > tmp/${base}/chm13v2/${base}_chm13v2_alignment.bam
	
# Dorado Alignment hg38
dorado aligner /scratch/tmp/thomachr/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
	-t 36 tmp/${base}/${base}_calls.bam > tmp/${base}/hg38/${base}_hg38_alignment.bam

module unload
ml palma/2020b GCC/10.2.0 OpenMPI/4.0.5 R/4.0.3

Rscript cnv_plot.R tmp/${base}/hg38/ results/${base}/cnv

# Modkit
modkit adjust-mods --convert h m tmp/${base}/chm13v2/${base}_chm13v2_alignment.bam tmp/${base}/${base}_calls_modkit.bam
modkit extract tmp/${base}/${base}_calls_modkit.bam tmp/${base}/${base}_calls_modkit.txt

# Sturgeon
module unload
ml palma/2021a  GCC/10.3.0  OpenMPI/4.1.1 ONNX-Runtime/1.10.0
sturgeon inputtobed -i tmp/${base}/ -o tmp/${base}/ -s modkit

sturgeon predict \
-i tmp/${base} \
-o results/${base} \
--model-files $sturgeon_model \
--plot-results

done
