suppressMessages(library(DecompPipeline))
suppressMessages(library(MeDeCom))
library(openxlsx)
parallel.setup(36)

rnb.options(import.table.separator = ";")
idat.dir <- "./idat"
tmp <- read.xlsx("./sample_sheet.xlsx")

write.table(tmp[tmp$array=="450K",],"samplesheet_450k.csv",sep=";",quote=FALSE)
write.table(tmp[tmp$array=="EPIC",],"samplesheet_EPIC.csv",sep=";",quote=FALSE)

samples_450k <- "./samplesheet_450k.csv"
samples_epic <- "./samplesheet_EPIC.csv"

rnb.450k <- rnb.execute.import(data.source = c(idat.dir,samples_450k),data.type = "infinium.idat.dir")
rnb.epic <- rnb.execute.import(data.source = c(idat.dir,samples_epic),data.type = "infinium.idat.dir")

rnb.set <- rnb.combine.arrays(rnb.450k, rnb.epic, type="common")

md.res <- start.decomp.pipeline(rnb.set=rnb.set,
                                analysis.name="decomp_10k_bmiq_2:10",
                                cores=64,
                                Ks=2:8,
                                lambda.grid=c(0,0.01,0.001,0.0001,0.00001),
                                execute.lump=TRUE,
                                factorviz.outputs=TRUE,
                                marker.selection="var",
                                n.markers=10000,
                                min.n.beads=3,
                                min.int.quant=0.05,
                                max.int.quant=0.95,
                                filter.beads=TRUE,
                                filter.context=TRUE,
                                filter.cross.reactive=TRUE,
                                filter.na=TRUE,
                                filter.snp=TRUE,
                                filter.sex.chromosomes=TRUE,
                                normalization="bmiq")

