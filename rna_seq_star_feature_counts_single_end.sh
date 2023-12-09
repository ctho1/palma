#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --partition normal
#SBATCH --time=2:00:00
#SBATCH --mem=80G
#SBATCH --job-name=RNAseq_alignment_single_end
#SBATCH --mail-type=ALL
#SBATCH --output log.txt
#SBATCH --mail-user=christian.thomas@ukmuenster.de

cd /scratch/tmp/thomachr/alignment

#!/bin/bash
ml palma/2019a GCC/8.2.0-2.31.1 SAMtools/1.9 Subread/2.0.0

# tell bash to be verbose and to abort on error
set -o pipefail
set -x -e -u

# set arguments
STAR_INDEX_DIR=/scratch/tmp/thomachr/references/STAR_index_hg38_ENSEMBL93
THREADS=24

in_files=$(find ./fastq -type f -name "*.fastq.gz" -print|sort)

for R1 in $in_files; do

R1_base=`basename $R1`

# align FastQ files
/scratch/tmp/thomachr/software/STAR-2.7.9a/source/STAR \
	--genomeDir "$STAR_INDEX_DIR" \
	--genomeLoad NoSharedMemory \
	--readFilesIn "$R1" \
	--readFilesCommand zcat \
	--outFileNamePrefix ./bam/${R1_base}_ \
	--outSAMtype BAM SortedByCoordinate \
	--runThreadN "$THREADS"

samtools index ./bam/${R1_base}_Aligned.sortedByCoord.out.bam

done

# prepare featureCount matrix

featureCounts -T 4 -s 2 -t exon -g gene_id -a /scratch/tmp/thomachr/references/ENSEMBL93.gtf -o ./cts/counts.txt ./bam/*.out.bam


#done
exit 0
