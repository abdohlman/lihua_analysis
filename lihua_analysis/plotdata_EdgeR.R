library(DiffBind)
#setwd("/Volumes/ShenLabData/Data/HighThroughput/LiHua/ATAC_2017Sep/diffbind/")
setwd("/Users/dohlman/Documents/Lihua_ATACSeq/lihua_analysis/")
#load("DiffBind.total.3pairs.allmethods.RData")

filelist <- c("DiffBind.total.no1.allmethods.RData",
              "DiffBind.total.no2.allmethods.RData",
              "DiffBind.total.no3.allmethods.RData",
              "DiffBind.total.no4.allmethods.RData")

for (filename in filelist) {

  load(filename)
  samplename <- strsplit(filename, ".", fixed=TRUE)[[1]][3]
  
  # dba.plotMA(atac.analyze)
  pdf(paste0('figures/',samplename,'.PCA.EDGER.RPKM.pdf'))
  dba.plotPCA(atac.analyze,DBA_CONDITION,label = DBA_FACTOR,attributes = DBA_CONDITION,method = DBA_EDGER,score = DBA_SCORE_RPKM)
  dev.off()
  
  pdf(paste0('figures/',samplename,'.PCA.EDGER.SUMMIT.pdf'))
  dba.plotPCA(atac.analyze,DBA_CONDITION,label = DBA_FACTOR,attributes = DBA_CONDITION,method = DBA_EDGER,score = DBA_SCORE_SUMMIT)
  dev.off()
  
  # dba.plotBox(atac.analyze)
  pdf(paste0('figures/',samplename,'.HEATMAP.EDGER.pdf'))
  dba.plotHeatmap(atac.analyze,attributes = DBA_CONDITION,method = DBA_EDGER,score = DBA_SCORE_RPKM)
  dev.off()
  
  pdf(paste0('figures/',samplename,'.HEATMAP.EDGER_BLOCK.pdf'))
  dba.plotHeatmap(atac.analyze,attributes = DBA_ID,method = DBA_EDGER_BLOCK,score = DBA_SCORE_RPKM)
  dev.off()
  
  library(rtracklayer)
  library(ggplot2)
  
  ## EDGER BLOCK MA PLOT
  
  # export(atac.DB,"diffbind.bed",format = "bed")
  atac.DB<-atac.DB.EDGERBLOCK
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
  
  p <- ggplot(data = df) +
    geom_point(shape=21,aes(x=ave,y=Fold,fill=sign,alpha=sign),color="grey90",size=2)+
    scale_fill_manual(values = c("1"="orangered3","0"="grey60","-1"="navyblue"))+
    scale_alpha_manual(values = c("1"=0.8,"0"=0.8,"-1"=0.8))+
    xlab("mean peaks signal")+ylab("fold change log2(colon/liver)")+#ggtitle(sample)+
    ylim(-2.5,3)+
    theme_bw()+
    theme(aspect.ratio = 0.5)
  
  method <- 'EDGER_BLOCK'
  
  ggsave(paste0('figures/',samplename,'.MAplot.',method,'.pdf'), plot = p, width = 8, height = 5)
  
  total.df<-df
  # total.df[4]<-rownames(total.df)
  
  pos_df<-df[df$sign=="1",]
  neg_df<-df[df$sign=="-1",]
  write.table(pos_df[,1:9],paste("./bed/diffbind.",samplename,".",method,".UP.bed",sep = ""),sep = "\t",quote = FALSE,col.names = FALSE,row.names = FALSE)
  write.table(pos_df[,1:4],paste("./bed/diffbind.",samplename,".",method,".UP.homer.bed",sep = ""),sep = "\t",quote = FALSE,col.names = FALSE,row.names = FALSE)
  write.table(neg_df[,1:9],paste("./bed/diffbind.",samplename,".",method,".DOWN.bed",sep = ""),sep = "\t",quote = FALSE,col.names = FALSE,row.names = FALSE)
  write.table(neg_df[,1:4],paste("./bed/diffbind.",samplename,".",method,".DOWN.homer.bed",sep = ""),sep = "\t",quote = FALSE,col.names = FALSE,row.names = FALSE)
  write.table(df,paste("./diffbind.","total.bed",sep = ""),sep = "\t",quote = FALSE,col.names = FALSE,row.names = FALSE)
  
  
  ## EDGER BLOCK MA PLOT
  
  atac.DB<-atac.DB.EDGER
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
  
  p <- ggplot(data = df) +
    geom_point(shape=21,aes(x=ave,y=Fold,fill=sign,alpha=sign),color="grey90",size=2)+
    scale_fill_manual(values = c("1"="orangered3","0"="grey60","-1"="navyblue"))+
    scale_alpha_manual(values = c("1"=0.8,"0"=0.8,"-1"=0.8))+
    xlab("mean peaks signal")+ylab("fold change log2(colon/liver)")+#ggtitle(sample)+
    ylim(-2.5,3)+
    theme_bw()+
    theme(aspect.ratio = 0.5)
  
  method <- 'EDGER'
  
  ggsave(paste0('figures/',samplename,'.MAplot.',method,'.pdf'), plot = p, width = 8, height = 5)
  
  #print(p)
  ### P.Value<0.05 fold change > 2 or fold change < 0.5
  total.df<-df
  # total.df[4]<-rownames(total.df)
  
  pos_df<-df[df$sign=="1",]
  neg_df<-df[df$sign=="-1",]
  write.table(pos_df[,1:9],paste("./bed/diffbind.",samplename,".",method,".UP.bed",sep = ""),sep = "\t",quote = FALSE,col.names = FALSE,row.names = FALSE)
  write.table(pos_df[,1:4],paste("./bed/diffbind.",samplename,".",method,".UP.homer.bed",sep = ""),sep = "\t",quote = FALSE,col.names = FALSE,row.names = FALSE)
  write.table(neg_df[,1:9],paste("./bed/diffbind.",samplename,".",method,".DOWN.bed",sep = ""),sep = "\t",quote = FALSE,col.names = FALSE,row.names = FALSE)
  write.table(neg_df[,1:4],paste("./bed/diffbind.",samplename,".",method,".DOWN.homer.bed",sep = ""),sep = "\t",quote = FALSE,col.names = FALSE,row.names = FALSE)
  write.table(df,paste("./diffbind.","total.bed",sep = ""),sep = "\t",quote = FALSE,col.names = FALSE,row.names = FALSE)
  
}