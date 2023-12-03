#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=36
#SBATCH --partition requeue
#SBATCH --time=4:00:00
#SBATCH --mem=120G
#SBATCH --job-name=sturgeon_pipeline
#SBATCH --mail-type=ALL
#SBATCH --error ./log/%x_%j.err.txt
#SBATCH --output ./log/%x_%j.out.txt
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
mkdir -p modkit_output/${base}
mkdir -p results/${base}

# Dorado
dorado basecaller --modified-bases 5mCG_5hmCG -x cpu \
	$dorado_model $dir > dorado_output/${base}_calls.bam

dorado aligner /scratch/tmp/thomachr/references/T2T/chm13v2.0.fa \
	dorado_output/${base}_calls.bam > dorado_output/${base}_aligned_calls.bam

dorado summary dorado_output/${base}_aligned_calls.bam > dorado_output/${base}_aligned_calls.summary.txt

# Modkit
modkit adjust-mods --convert h m ./dorado_output/${base}_aligned_calls.bam ./modkit_output/${base}/${base}_calls_modkit.bam
modkit extract modkit_output/${base}/${base}_calls_modkit.bam modkit_output/${base}/${base}_calls_modkit.txt

# Sturgeon
ml palma/2021a  GCC/10.3.0  OpenMPI/4.1.1 ONNX-Runtime/1.10.0
sturgeon inputtobed -i modkit_output/${base} -o modkit_output/${base} -s modkit

sturgeon predict \
-i modkit_output/${base} \
-o results/${base} \
--model-files $sturgeon_model \
--plot-results

done
