#.libPaths( c( .libPaths(), "/gpfs/fs0/home/ew152/Rlibs") )
#library("S4Vectors", lib.loc=.libPaths()[3])
#library("RSQLite", lib.loc=.libPaths()[3])
library(DiffBind)
setwd("/Users/dohlman/Documents/Lihua_ATACSeq/lihua_analysis/")

# filelist <- c("atac.sampletable.3pairs.no1.csv",
#               "atac.sampletable.3pairs.no2.csv",
#               "atac.sampletable.3pairs.no3.csv",
#               "atac.sampletable.3pairs.no4.csv")

filelist <- c("atac.sampletable.liver-colon.allsamples.csv")

# for (filename in filelist) {
  rm(list = ls(pattern = "^atac"))
  tag <- strsplit(filename, '.', fixed=TRUE)[[1]][4]
  atac.list <- dba(sampleSheet=filename)
  atac.count <- dba.count(atac.list)
  atac.contrast <- dba.contrast(atac.count, categories=DBA_CONDITION, block=DBA_REPLICATE, minMembers = 2)
  atac.analyze <- dba.analyze(atac.contrast, method=DBA_ALL_METHODS)
  atac.DB.DESEQ2BLOCK <- dba.report(atac.analyze, file=paste("atac",tag,"DESeq2Block","total",sep='.'), th = 1, method = DBA_DESEQ2_BLOCK)
  atac.DB.EDGERBLOCK <- dba.report(atac.analyze, file = paste("atac",tag,"EdgeRBlock","total",sep='.'), th = 1, method = DBA_EDGER_BLOCK)
  atac.DB.DESEQ2 <- dba.report(atac.analyze, file = paste("atac",tag,"DESeq2","total",sep='.'), th = 1, method = DBA_DESEQ2)
  atac.DB.EDGER <- dba.report(atac.analyze, file = paste("atac",tag,"EdgeR","total",sep='.'), th = 1, method = DBA_EDGER)
  save(atac.list, atac.count, atac.contrast, atac.analyze,
       atac.DB.DESEQ2BLOCK,
       atac.DB.DESEQ2,
       atac.DB.EDGERBLOCK,
       atac.DB.EDGER,
       file = paste("DiffBind","total",tag,"allmethods","RData",sep='.'))
# }

