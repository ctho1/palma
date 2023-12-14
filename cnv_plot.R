.libPaths('/scratch/tmp/thomachr/software/R/')
library(ACE)

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
       cap = 12, 
       bottom = 0, 
       trncname = FALSE, 
       printsummaries = TRUE,
       autopick = TRUE)
