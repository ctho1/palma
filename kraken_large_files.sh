#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --partition normal
#SBATCH --time=4:00:00
#SBATCH --mem=80G
#SBATCH --job-name=kraken2_paired_end
#SBATCH --mail-type=ALL
#SBATCH --output log_kraken2_paired_end.txt
#SBATCH --mail-user=christian.thomas@ukmuenster.de

ml palma/2021b  GCC/11.2.0  OpenMPI/4.1.1
ml Kraken2/2.1.2

# tell bash to be verbose and to abort on error
set -o pipefail
set -x -e -u

## Parameters ############################################################################
KRAKEN_PATH=/scratch/tmp/thomachr/metagenomics/kraken2/
DATABABASE=${KRAKEN_PATH}/kraken_db/k2_standard_20220607
THREADS=8
NAME=`basename $1 .fastq.gz`

## Run Kraken2 ###########################################################################
kraken2 \
	--db $DATABABASE \
	--threads $THREADS \
	--minimum-hit-groups 3 \
	--use-names \
	--report-minimizer-data \
	--output - \
	--report ./report/${NAME}.kraken.report.txt \
	--paired $1 $2

exit 0
