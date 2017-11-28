library(DiffBind)
# setwd("/Volumes/ShenLabData/Data/HighThroughput/LiHua/ATAC_2017Sep/diffbind/")
setwd("/Users/dohlman/Documents/Lihua_ATACSeq/lihua_analysis/")
load("DiffBind.total.no1.allmethods.RData")
# dba.plotMA(atac.analyze)
pdf('PCA.kyc.1.pdf')
dba.plotPCA(atac.analyze,DBA_CONDITION,label = DBA_FACTOR,attributes = DBA_CONDITION,method = DBA_DESEQ2,score = DBA_SCORE_RPKM)
dev.off()
pdf('PCA.kyc.101.pdf')
dba.plotPCA(atac.analyze,DBA_CONDITION,label = DBA_FACTOR,attributes = DBA_CONDITION,method = DBA_DESEQ2,score = DBA_SCORE_SUMMIT)
dev.off()


# dba.plotBox(atac.analyze)
dba.plotHeatmap(atac.analyze,attributes = DBA_CONDITION,method = DBA_DESEQ2,score = DBA_SCORE_RPKM)
dba.plotHeatmap(atac.analyze,attributes = DBA_ID,method = DBA_DESEQ2_BLOCK,score = DBA_SCORE_RPKM)

library(rtracklayer)
# export(atac.DB,"diffbind.bed",format = "bed")
atac.DB<-atac.DB.DESEQ2BLOCK
anadata<-data.frame(mcols(atac.DB))
df <- data.frame(chr=seqnames(atac.DB),
                 starts=start(atac.DB)-1,
                 ends=end(atac.DB),
                 names=c(rep(".", length(atac.DB))),
                 scores=c(rep(".", length(atac.DB))),
                 strands=strand(atac.DB),
                 pval=anadata$p.value,
                 FDR =anadata$FDR,
                 Fold=anadata$Fold,
                 ave=anadata$Conc
)

df$sign<-"0"
df$sign[df$pval<=0.05 & df$Fold> 1]<-  "1"
df$sign[df$pval<=0.05 & df$Fold< -1]<- "-1"

library(ggplot2)
ggplot(data = df)+
  geom_point(shape=21,aes(x=ave,y=Fold,fill=sign,alpha=sign),color="grey90",size=2)+
  scale_fill_manual(values = c("1"="orangered3","0"="grey60","-1"="navyblue"))+
  scale_alpha_manual(values = c("1"=0.8,"0"=0.8,"-1"=0.8))+
  xlab("mean peaks signal")+ylab("fold change log2(colon/liver)")+#ggtitle(sample)+
  ylim(-2.5,3)+
  theme_bw()+
  theme(aspect.ratio = 0.5)
### P.Value<0.05 fold change > 2 or fold change < 0.5
total.df<-df
# total.df[4]<-rownames(total.df)

pos_df<-df[df$sign=="1",]
neg_df<-df[df$sign=="-1",]
write.table(pos_df[,1:9],paste("./homer/diffbind.",samplename,".UP.bed",sep = ""),sep = "\t",quote = FALSE,col.names = FALSE,row.names = FALSE)
write.table(pos_df[,1:4],paste("./homer/diffbind.",samplename,".UP.homer.bed",sep = ""),sep = "\t",quote = FALSE,col.names = FALSE,row.names = FALSE)
write.table(neg_df[,1:9],paste("./homer/diffbind.",samplename,".DOWN.bed",sep = ""),sep = "\t",quote = FALSE,col.names = FALSE,row.names = FALSE)
write.table(neg_df[,1:4],paste("./homer/diffbind.",samplename,".DOWN.homer.bed",sep = ""),sep = "\t",quote = FALSE,col.names = FALSE,row.names = FALSE)
write.table(df,paste("./diffbind.","total.bed",sep = ""),sep = "\t",quote = FALSE,col.names = FALSE,row.names = FALSE)
