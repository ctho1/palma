#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=36
#SBATCH --partition normal
#SBATCH --time=1:00:00
#SBATCH --mem=80G
#SBATCH --job-name=ngs_pipeline
#SBATCH --mail-type=ALL
#SBATCH --error ./log/%x_%j.err.txt
#SBATCH --output ./log/%x_%j.out.txt

## Parameters ############################################################################
BASE_DIR=/scratch/tmp/thomachr/arriba/arriba_v2.4.0
RefGenome=/scratch/tmp/thomachr/references/hg19/hg19.fa
STAR_INDEX_DIR=$BASE_DIR/STAR_index_GRCh37viral_GENCODE19
ANNOTATION_GTF=$BASE_DIR/GENCODE19.gtf
ASSEMBLY_FA=$BASE_DIR/GRCh37viral.fa
BLACKLIST_TSV=$BASE_DIR/database/blacklist_hg19_hs37d5_GRCh37_v2.4.0.tsv.gz
KNOWN_FUSIONS_TSV=$BASE_DIR/database/known_fusions_hg19_hs37d5_GRCh37_v2.4.0.tsv.gz
TAGS_TSV="$KNOWN_FUSIONS_TSV" 
PROTEIN_DOMAINS_GFF3=$BASE_DIR/database/protein_domains_hg19_hs37d5_GRCh37_v2.4.0.gff3
DATABABASE=/scratch/tmp/thomachr/references/kraken/k2_standard_20230314
EUPATHDB=/scratch/tmp/thomachr/references/kraken/k2_eupathdb48_20230407
CNV_REF=/scratch/tmp/thomachr/references/panel_v4_reference.cnn
export PATH="/home/t/thomachr/bin:$PATH"
THREADS=36
R1_base=`basename $1 .fastq.gz`
mkdir -p ./tmp
mkdir -p ./output
mkdir -p ./log
TMP_DIR=./tmp/${R1_base}
OUT_DIR=./output/${R1_base}
mkdir -p "$TMP_DIR"
mkdir -p "$OUT_DIR"
mkdir -p "$OUT_DIR"/cnv
mkdir -p "$OUT_DIR"/kraken

### CNV calling with CNVkit ################################################################
module purge
ml palma/2019a GCC/8.2.0-2.31.1 BWA/0.7.17
bwa mem -t $THREADS $RefGenome $1 $2 > "$TMP_DIR"/${R1_base}_aligned.sam

module purge
ml palma/2022a GCC/11.3.0 SAMtools/1.16.1
samtools view -@ "$THREADS" -S -b "$TMP_DIR"/${R1_base}_aligned.sam > "$TMP_DIR"/${R1_base}_aligned.bam
samtools sort -@ "$THREADS" -m $((40000/THREADS))M -T tmp -O bam "$TMP_DIR"/${R1_base}_aligned.bam -o "$TMP_DIR"/${R1_base}_aligned_sorted.bam
samtools index "$TMP_DIR"/${R1_base}_aligned_sorted.bam

module purge
ml palma/2022b foss/2022b CNVkit/0.9.10-R-4.2.2
cnvkit.py batch "$TMP_DIR"/${R1_base}_aligned_sorted.bam --reference $CNV_REF --processes 0 \
	--drop-low-coverage --output-dir $OUT_DIR/cnv --diagram

# all chromosomes
cnvkit.py scatter -s $OUT_DIR/cnv/*aligned_sorted.cn{s,r} --segment-color 'blue' --fig-size 10 4 --y-min "-2" --y-max "4" \
	-o $OUT_DIR/cnv_plot.pdf

# chr1
cnvkit.py scatter -s $OUT_DIR/cnv/*aligned_sorted.cn{s,r} -c chr1 --segment-color 'blue' \
-o "$OUT_DIR"/cnv/${R1_base}_chr1.pdf
# chr2
cnvkit.py scatter -s $OUT_DIR/cnv/*aligned_sorted.cn{s,r} -c chr2 --segment-color 'blue' \
-o "$OUT_DIR"/cnv/${R1_base}_chr2.pdf
# chr3
cnvkit.py scatter -s $OUT_DIR/cnv/*aligned_sorted.cn{s,r} -c chr3 --segment-color 'blue' \
-o "$OUT_DIR"/cnv/${R1_base}_chr3.pdf
# chr4
cnvkit.py scatter -s $OUT_DIR/cnv/*aligned_sorted.cn{s,r} -c chr4 --segment-color 'blue' \
-g PDGFRA,FGFR3 -o "$OUT_DIR"/cnv/${R1_base}_chr4.pdf
# chr5
cnvkit.py scatter -s $OUT_DIR/cnv/*aligned_sorted.cn{s,r} -c chr5 --segment-color 'blue' \
-o "$OUT_DIR"/cnv/${R1_base}_chr5.pdf
# chr6
cnvkit.py scatter -s $OUT_DIR/cnv/*aligned_sorted.cn{s,r} -c chr6 --segment-color 'blue' \
-o "$OUT_DIR"/cnv/${R1_base}_chr6.pdf
# chr7
cnvkit.py scatter -s $OUT_DIR/cnv/*aligned_sorted.cn{s,r} -c chr7 --segment-color 'blue' \
-g EGFR,MET -o "$OUT_DIR"/cnv/${R1_base}_chr7.pdf
cp "$OUT_DIR"/cnv/${R1_base}_chr7.pdf "$OUT_DIR"/EGFR_MET.pdf

# chr8
cnvkit.py scatter -s $OUT_DIR/cnv/*aligned_sorted.cn{s,r} -c chr8 --segment-color 'blue' \
-g FGFR1 -o "$OUT_DIR"/cnv/${R1_base}_chr8.pdf
# chr9
cnvkit.py scatter -s $OUT_DIR/cnv/*aligned_sorted.cn{s,r} -c chr9 --segment-color 'blue' \
-g CDKN2A -o "$OUT_DIR"/cnv/${R1_base}_chr9.pdf
cp "$OUT_DIR"/cnv/${R1_base}_chr9.pdf "$OUT_DIR"/CDKN2A.pdf

# chr10
cnvkit.py scatter -s $OUT_DIR/cnv/*aligned_sorted.cn{s,r} -c chr10 --segment-color 'blue' \
-g PTEN -o "$OUT_DIR"/cnv/${R1_base}_chr10.pdf
# chr11
cnvkit.py scatter -s $OUT_DIR/cnv/*aligned_sorted.cn{s,r} -c chr11 --segment-color 'blue' \
-o "$OUT_DIR"/cnv/${R1_base}_chr11.pdf
# chr12
cnvkit.py scatter -s $OUT_DIR/cnv/*aligned_sorted.cn{s,r} -c chr12 --segment-color 'blue' \
-o "$OUT_DIR"/cnv/${R1_base}_chr12.pdf
# chr13
cnvkit.py scatter -s $OUT_DIR/cnv/*aligned_sorted.cn{s,r} -c chr13 --segment-color 'blue' \
-o "$OUT_DIR"/cnv/${R1_base}_chr13.pdf
# chr14
cnvkit.py scatter -s $OUT_DIR/cnv/*aligned_sorted.cn{s,r} -c chr14 --segment-color 'blue' \
-o "$OUT_DIR"/cnv/${R1_base}_chr14.pdf
# chr15
cnvkit.py scatter -s $OUT_DIR/cnv/*aligned_sorted.cn{s,r} -c chr15 --segment-color 'blue' \
-o "$OUT_DIR"/cnv/${R1_base}_chr15.pdf
# chr16
cnvkit.py scatter -s $OUT_DIR/cnv/*aligned_sorted.cn{s,r} -c chr16 --segment-color 'blue' \
-o "$OUT_DIR"/cnv/${R1_base}_chr16.pdf
# chr17
cnvkit.py scatter -s $OUT_DIR/cnv/*aligned_sorted.cn{s,r} -c chr17 --segment-color 'blue' \
-g TP53 -o "$OUT_DIR"/cnv/${R1_base}_chr17.pdf
# chr18
cnvkit.py scatter -s $OUT_DIR/cnv/*aligned_sorted.cn{s,r} -c chr18 --segment-color 'blue' \
-o "$OUT_DIR"/cnv/${R1_base}_chr18.pdf
# chr19
cnvkit.py scatter -s $OUT_DIR/cnv/*aligned_sorted.cn{s,r} -c chr19 --segment-color 'blue' \
-o "$OUT_DIR"/cnv/${R1_base}_chr19.pdf
# chr20
cnvkit.py scatter -s $OUT_DIR/cnv/*aligned_sorted.cn{s,r} -c chr20 --segment-color 'blue' \
-o "$OUT_DIR"/cnv/${R1_base}_chr20.pdf
# chr21
cnvkit.py scatter -s $OUT_DIR/cnv/*aligned_sorted.cn{s,r} -c chr21 --segment-color 'blue' \
-o "$OUT_DIR"/cnv/${R1_base}_chr21.pdf
# chr22
cnvkit.py scatter -s $OUT_DIR/cnv/*aligned_sorted.cn{s,r} -c chr22 --segment-color 'blue' \
-g SMARCB1 -o "$OUT_DIR"/cnv/${R1_base}_chr22.pdf


### Arriba #################################################################################
module unload
ml palma/2022a GCC/11.3.0 SAMtools/1.16.1 OpenMPI/4.1.4 R-bundle-Bioconductor/3.15-R-4.2.1
/scratch/tmp/thomachr/software/STAR-2.7.10b/source/STAR \
	--runThreadN "$THREADS" \
	--outFileNamePrefix "$TMP_DIR"/ \
	--genomeDir "$STAR_INDEX_DIR" --genomeLoad NoSharedMemory \
	--readFilesIn "$1" "$2" --readFilesCommand zcat \
	--outStd BAM_Unsorted --outSAMtype BAM Unsorted --outSAMunmapped Within --outBAMcompression 0 \
	--outFilterMultimapNmax 50 --peOverlapNbasesMin 10 --alignSplicedMateMapLminOverLmate 0.5 --alignSJstitchMismatchNmax 5 -1 5 5 \
	--chimSegmentMin 10 --chimOutType WithinBAM HardClip --chimJunctionOverhangMin 10 --chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --chimSegmentReadGapMax 3 --chimMultimapNmax 50 |

tee "$TMP_DIR"/${R1_base}_Aligned.out.bam |

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

$BASE_DIR/scripts/quantify_virus_expression.sh "$TMP_DIR"/${R1_base}_Aligned.sortedByCoord.out.bam "$OUT_DIR"/${R1_base}virus_expression.tsv

### Kraken2 #################################################################################
module unload
ml palma/2021b  GCC/11.2.0  OpenMPI/4.1.1 Kraken2/2.1.2 Bracken/2.7

kraken2 \
	--db $DATABABASE \
	--threads 8 \
	--minimum-hit-groups 3 \
	--use-names \
	--report-minimizer-data \
	--output - \
	--report "$OUT_DIR"/kraken/${R1_base}.kraken.report.txt \
	--paired "$1" "$2"

bracken \
	-d $DATABABASE \
	-i "$OUT_DIR"/kraken/${R1_base}.kraken.report.txt \
	-r 75 \
	-l G \
	-t 3 \
	-o "$OUT_DIR"/kraken/${R1_base}.bracken_genus.txt \
	-w "$OUT_DIR"/kraken/${R1_base}.bracken_genus.report.txt      

bracken \
	-d $DATABABASE \
	-i "$OUT_DIR"/kraken/${R1_base}.kraken.report.txt \
	-r 75 \
	-l S \
	-t 3 \
	-o "$OUT_DIR"/kraken/${R1_base}.bracken_species.txt \
	-w "$OUT_DIR"/kraken/${R1_base}.bracken_species.report.txt

kraken2 \
	--db $EUPATHDB \
	--threads 8 \
	--minimum-hit-groups 3 \
	--use-names \
	--report-minimizer-data \
	--output - \
	--report "$OUT_DIR"/kraken/${R1_base}.kraken.eupathdb.report.txt \
	--paired "$1" "$2"

bracken \
	-d $EUPATHDB \
	-i "$OUT_DIR"/kraken/${R1_base}.kraken.eupathdb.report.txt \
	-r 75 \
	-l G \
	-t 3 \
	-o "$OUT_DIR"/kraken/${R1_base}.eupathdb.bracken_genus.txt \
	-w "$OUT_DIR"/kraken/${R1_base}.eupathdb.bracken_genus.report.txt        

bracken \
	-d $EUPATHDB \
	-i "$OUT_DIR"/kraken/${R1_base}.kraken.eupathdb.report.txt \
	-r 75 \
	-l S \
	-t 3 \
	-o "$OUT_DIR"/kraken/${R1_base}.eupathdb.bracken_species.txt \
	-w "$OUT_DIR"/kraken/${R1_base}.eupathdb.bracken_species.report.txt

### Create summary report ###################################################################
module purge
#ml palma/2020b GCC/10.2.0 OpenMPI/4.0.5 R/4.0.3 Pandoc/2.10
ml palma/2023b GCC/13.2.0 R/4.4.1 Pandoc/3.1.2
parent_dir=$(dirname "$1")
sub_dir=$(dirname "$parent_dir")
Rscript /home/t/thomachr/scripts/ngs_report.R $sub_dir/output/${R1_base} ${R1_base}
