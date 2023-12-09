#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --partition normal
#SBATCH --time=1:00:00
#SBATCH --mem=40G
#SBATCH --job-name=salmon
#SBATCH --mail-type=ALL
#SBATCH --output log.txt
#SBATCH --mail-user=christian.thomas@ukmuenster.de

R1_base=`basename $1 .fastq.gz`

salmon quant -i /scratch/tmp/thomachr/references/salmon/hg38_partial_sa \
	         -l A \
             -r $1 \
             -p 24 \
	     	 --validateMappings \
	         --seqBias \
	     	 --gcBias \
             -o ./output/${R1_base}