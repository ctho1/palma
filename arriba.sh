#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --partition normal
#SBATCH --time=3:00:00
#SBATCH --mem=80G
#SBATCH --job-name=arriba
#SBATCH --error ./log/%x_%j.err.txt
#SBATCH --output ./log/%x_%j.out.txt

ml palma/2019a
ml GCC/8.2.0-2.31.1
ml OpenMPI/3.1.3
ml SAMtools/1.9
ml R-bundle-Bioconductor/3.9-R-3.6.0
mkdir -p log

# set arguments
BASE_DIR=/scratch/tmp/thomachr/arriba/arriba_v2.4.0
STAR_INDEX_DIR=$BASE_DIR/STAR_index_GRCh37viral_GENCODE19
ANNOTATION_GTF=$BASE_DIR/GENCODE19.gtf
ASSEMBLY_FA=$BASE_DIR/GRCh37viral.fa
BLACKLIST_TSV=$BASE_DIR/database/blacklist_hg19_hs37d5_GRCh37_v2.4.0.tsv.gz
KNOWN_FUSIONS_TSV=$BASE_DIR/database/known_fusions_hg19_hs37d5_GRCh37_v2.4.0.tsv.gz
TAGS_TSV="$KNOWN_FUSIONS_TSV" # different files can be used for filtering and tagging, but the provided one can be used for both
PROTEIN_DOMAINS_GFF3=$BASE_DIR/database/protein_domains_hg19_hs37d5_GRCh37_v2.4.0.gff3
DATABABASE=/scratch/tmp/thomachr/metagenomics/kraken2/kraken_db/k2_standard_20220607 # Kraken2 Standard Database
EUPATHDB=/scratch/tmp/thomachr/metagenomics/kraken2/kraken_db/EuPathDB48 # Kraken2 EuPathDB48 Database
CENTRIFUGE_INDEX=/scratch/tmp/thomachr/software/centrifuge-1.0.4/indices/hpvc
THREADS=24
R1_base=`basename $1 .fastq.gz`
TMP_DIR=./tmp/${R1_base}
OUT_DIR=./out/${R1_base}
mkdir "$TMP_DIR"
mkdir "$OUT_DIR"

# align FastQ files (STAR >=2.7.6a recommended)
/scratch/tmp/thomachr/software/STAR-2.7.10b/source/STAR \
	--runThreadN "$THREADS" \
	--outFileNamePrefix "$TMP_DIR"/ \
	--genomeDir "$STAR_INDEX_DIR" --genomeLoad NoSharedMemory \
	--readFilesIn "$1" "$2" --readFilesCommand zcat \
	--outStd BAM_Unsorted --outSAMtype BAM Unsorted --outSAMunmapped Within --outBAMcompression 0 \
	--outFilterMultimapNmax 50 --peOverlapNbasesMin 10 --alignSplicedMateMapLminOverLmate 0.5 --alignSJstitchMismatchNmax 5 -1 5 5 \
	--chimSegmentMin 10 --chimOutType WithinBAM HardClip --chimJunctionOverhangMin 10 --chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --chimSegmentReadGapMax 3 --chimMultimapNmax 50 |

tee "$TMP_DIR"/${R1_base}_Aligned.out.bam |

# call arriba
# -O "$OUT_DIR"/${R1_base}fusions_discarded.tsv

"$BASE_DIR/arriba" \
	-x /dev/stdin -I \
	-o "$OUT_DIR"/${R1_base}fusions.tsv -f intronic,in_vitro,internal_tandem_duplication \
	-a "$ASSEMBLY_FA" -g "$ANNOTATION_GTF" -b "$BLACKLIST_TSV" -k "$KNOWN_FUSIONS_TSV" -t "$TAGS_TSV" -p "$PROTEIN_DOMAINS_GFF3"

# sorting and indexing is only required for visualization
samtools sort -@ "$THREADS" -m $((40000/THREADS))M -T tmp -O bam "$TMP_DIR"/${R1_base}_Aligned.out.bam > "$TMP_DIR"/${R1_base}_Aligned.sortedByCoord.out.bam
rm -f "$TMP_DIR"_Aligned.out.bam
samtools index "$TMP_DIR"/${R1_base}_Aligned.sortedByCoord.out.bam

# Plot Coverage of NGS Panel Genes
samtools depth -d 0 -b $BASE_DIR/coverage_regions.tsv "$TMP_DIR"/${R1_base}_Aligned.sortedByCoord.out.bam > "$TMP_DIR"/${R1_base}_panel_coverage.txt
Rscript --vanilla $BASE_DIR/coverage_plot.R $BASE_DIR/coverage_regions.tsv "$TMP_DIR"/${R1_base}_panel_coverage.txt "$OUT_DIR"/${R1_base}_panel_coverage.pdf

$BASE_DIR/draw_fusions.R \
    --fusions="$OUT_DIR"/${R1_base}fusions.tsv \
    --alignments="$TMP_DIR"/${R1_base}_Aligned.sortedByCoord.out.bam \
    --output="$OUT_DIR"/${R1_base}_fusions.pdf \
    --annotation="$ANNOTATION_GTF" \
    --cytobands=$BASE_DIR/database/cytobands_hg19_hs37d5_GRCh37_v2.4.0.tsv \
    --proteinDomains="$PROTEIN_DOMAINS_GFF3"

#$BASE_DIR/scripts/extract_fusion-supporting_alignments.sh "$OUT_DIR"fusions.tsv "$TMP_DIR"/${R1_base}_Aligned.sortedByCoord.out.bam "$OUT_DIR"
$BASE_DIR/scripts/quantify_virus_expression.sh "$TMP_DIR"/${R1_base}_Aligned.sortedByCoord.out.bam "$OUT_DIR"/${R1_base}virus_expression.tsv

module unload

ml palma/2021b  GCC/11.2.0  OpenMPI/4.1.1 
ml Kraken2/2.1.2 Bracken/2.7

kraken2 \
	--db $DATABABASE \
	--threads 8 \
	--minimum-hit-groups 3 \
	--use-names \
	--report-minimizer-data \
	--output - \
	--report "$OUT_DIR"/${R1_base}.kraken.report.txt \
	--paired "$1" "$2"

bracken \
	-d $DATABABASE \
	-i "$OUT_DIR"/${R1_base}.kraken.report.txt \
	-r 75 \
	-l G \
	-t 10 \
	-o "$OUT_DIR"/${R1_base}.bracken_genus.txt \
	-w "$OUT_DIR"/${R1_base}.bracken_genus.report.txt      

bracken \
	-d $DATABABASE \
	-i "$OUT_DIR"/${R1_base}.kraken.report.txt \
	-r 75 \
	-l S \
	-t 10 \
	-o "$OUT_DIR"/${R1_base}.bracken_species.txt \
	-w "$OUT_DIR"/${R1_base}.bracken_species.report.txt

kraken2 \
	--db $EUPATHDB \
	--threads 8 \
	--minimum-hit-groups 3 \
	--use-names \
	--report-minimizer-data \
	--output - \
	--report "$OUT_DIR"/${R1_base}.kraken.eupathdb.report.txt \
	--paired "$1" "$2"

bracken \
	-d $EUPATHDB \
	-i "$OUT_DIR"/${R1_base}.kraken.eupathdb.report.txt \
	-r 75 \
	-l G \
	-t 10 \
	-o "$OUT_DIR"/${R1_base}.eupathdb.bracken_genus.txt \
	-w "$OUT_DIR"/${R1_base}.eupathdb.bracken_genus.report.txt        

bracken \
	-d $EUPATHDB \
	-i "$OUT_DIR"/${R1_base}.kraken.eupathdb.report.txt \
	-r 75 \
	-l S \
	-t 10 \
	-o "$OUT_DIR"/${R1_base}.eupathdb.bracken_species.txt \
	-w "$OUT_DIR"/${R1_base}.eupathdb.bracken_species.report.txt 

module unload
ml palma/2019b  GCC/8.3.0  OpenMPI/3.1.4 HISAT2/2.2.1
export PATH=/scratch/tmp/thomachr/software/centrifuge-1.0.4:$PATH

centrifuge -q \
-1 "$1" \
-2 "$2" \
--report-file "$OUT_DIR"/${R1_base}.centrifuge.txt \
--threads $THREADS \
-x $CENTRIFUGE_INDEX

