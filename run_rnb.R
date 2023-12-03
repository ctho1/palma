library(openxlsx)
library(RnBeads)
parallel.setup(24)
options(bitmapType="cairo")

anno <- read.xlsx("sample_sheet.xlsx")
anno$Biobank.ID <- paste0("#",anno$Biobank.ID)
anno$Age <- as.numeric(anno$Age)
anno <- anno[,c("Biobank.ID","Sentrix_ID","Sentrix_Position","Diagnosis","Group","Age")]
write.csv(anno, file="sample_annotation.csv", quote=FALSE, row.names = FALSE)

data.dir <- "."
idat.dir <- file.path(data.dir, "idat")
sample.annotation <- file.path(data.dir, "sample_annotation.csv")
analysis.dir <- "./analysis"
report.dir <- file.path(analysis.dir, "reports")

rnb.options(analysis.name = "DNA Methylation Profiling of CD8+ T cells of RRMS, PPMS and HD across Age",
            filtering.sex.chromosomes.removal=TRUE,
            identifiers.column="Biobank.ID",
	        normalization.background.method = "none",
            #qc.cnv = TRUE, # causes error 
            inference.age.column = "Age",
	        filtering.snp = "any",
            filtering.cross.reactive = TRUE,
	        filtering.missing.value.quantile = 0.5,
	        differential.adjustment.sva = FALSE,
	        differential.adjustment.celltype = FALSE)

rnb.run.analysis(dir.reports=report.dir, sample.sheet=sample.annotation,
                 data.dir=idat.dir, data.type="infinium.idat.dir")
