.libPaths('/scratch/tmp/thomachr/software/R/')
library(ACE)
library(QDNAseq)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

runACE(inputdir = args[1], 
       outputdir = args[2],
       filetype = 'bam', 
       genome = 'hg38', 
       c(500,1000), 
       ploidies = 2, 
       imagetype = 'pdf', 
       method = 'RMSE', 
       penalty = 0, 
       cap = 4, 
       bottom = 0, 
       trncname = FALSE, 
       printsummaries = TRUE,
       autopick = TRUE)

tmp <- readRDS(paste0(args[2],"/1000kbp.rds"))

pdf(paste0(args[2],"/cnv_plot_1000kbp.pdf"),width=10,height=5)
singleplot(tmp, ploidy=2, QDNAseqobjectsample = 1, onlyautosomes = TRUE) + 
  ylim(0,4)
dev.off()

pdf(paste0(args[2],"/cnv_plot_1000kbp_col.pdf"),width=11,height=6)
ACEcall(tmp, QDNAseqobjectsample = TRUE, ploidy = 2, cap = 5)$calledplot
dev.off()

tmp <- readRDS(paste0(args[2],"/500kbp.rds"))

pdf(paste0(args[2],"/cnv_plot_500kbp.pdf"),width=10,height=5)
singleplot(tmp, ploidy=2, QDNAseqobjectsample = 1, onlyautosomes = TRUE) + 
  ylim(0,4)
dev.off()