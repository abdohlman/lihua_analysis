library(DiffBind)
setwd("/Volumes/ShenLabData/Data/HighThroughput/LiHua/ATAC_2017Sep/diffbind/")
atac.list <- dba(sampleSheet="atac.sampletable.3pairs.csv")
atac.count <- dba.count(atac.list)
atac.contrast <- dba.contrast(atac.count,categories=DBA_CONDITION, block=DBA_REPLICATE,minMembers = 2)# block=DBA_REPLICATE for paired samples, minMembers for number of replicats <3
# atac.contrast <- dba.contrast(atac.count,categories=DBA_CONDITION, block=DBA_REPLICATE)# block=DBA_REPLICATE for paired samples, minMembers for number of replicats <3
atac.analyze <- dba.analyze(atac.contrast,method=DBA_ALL_METHODS) 
atac.DB.DESEQ2BLOCK <- dba.report(atac.analyze,file = "atac.P1P3.DESeq2Block.total",th = 1,method = DBA_DESEQ2_BLOCK)
atac.DB.EDGERBLOCK <- dba.report(atac.analyze,file = "atac.P1P3.EdgeRBlock.total",th = 1,method = DBA_EDGER_BLOCK)
atac.DB.DESEQ2 <- dba.report(atac.analyze,file = "atac.P1P3.DESeq2.total",th = 1,method = DBA_DESEQ2)
atac.DB.EDGER <- dba.report(atac.analyze,file = "atac.P1P3.EdgeR.total",th = 1,method = DBA_EDGER)
save(atac.list,atac.count,atac.contrast,atac.analyze,atac.DB.DESEQ2BLOCK,atac.DB.DESEQ2,atac.DB.EDGERBLOCK,atac.DB.EDGER, file = "DiffBind.total.P1P3.allmethods.RData")
