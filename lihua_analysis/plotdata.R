library(DiffBind)
library(rtracklayer)
library(ggplot2)
#setwd("/Volumes/ShenLabData/Data/HighThroughput/LiHua/ATAC_2017Sep/diffbind/")
setwd("/Users/ergangwang/Downloads/Lihua_ATACSeq/seq2_analysis/")
#load("DiffBind.total.3pairs.allmethods.RData")

# filelist <- c("DiffBind.total.no1.allmethods.RData",
#               "DiffBind.total.no2.allmethods.RData",
#               "DiffBind.total.no3.allmethods.RData",
#               "DiffBind.total.no4.allmethods.RData",
#               "DiffBind.total.allsamples.allmethods.RData")
###################################################
#if do liver.colon  need to change the ylab below
###################################################
filelist <- c("DiffBind.total.liver.allmethods.RData")

methods <- c(DBA_DESEQ2,
             DBA_DESEQ2_BLOCK,
             DBA_EDGER,
             DBA_EDGER_BLOCK)

# methods <- c(atac.DB.DESEQ2BLOCK,
#              atac.DB.DESEQ2BLOCK,
#              atac.DB.EDGER,
#              atac.DB.EDGERBLOCK)

# PCA
# pdf('figures/PCA.DESEQ.RPKM.pdf')
# par(mfrow=c(4,4))

for (filename in filelist) {

  rm(list = ls(pattern = "^atac"))
  load(filename)
  samplename <- strsplit(filename, ".", fixed=TRUE)[[1]][3]
  
  pdf(paste('figures/PCA',samplename,'rpkm','pdf',sep='.'))
  dba.plotPCA(atac.analyze, DBA_CONDITION, label = DBA_REPLICATE, attributes = DBA_CONDITION,score = DBA_SCORE_RPKM)
  dev.off()
  
  pdf(paste('figures/PCA',samplename,'summit', 'pdf',sep='.'))
  dba.plotPCA(atac.analyze, DBA_CONDITION, label = DBA_REPLICATE, attributes = DBA_CONDITION,score = DBA_SCORE_SUMMIT)
  dev.off()
  
  pdf(paste('figures/HEATMAP',samplename,'rpkm','pdf',sep='.'))
  dba.plotHeatmap(atac.analyze,score = DBA_SCORE_RPKM)
  dev.off()
  
  pdf(paste('figures/HEATMAP',samplename,'summit','pdf',sep='.'))
  dba.plotHeatmap(atac.analyze,score = DBA_SCORE_SUMMIT)
  dev.off()
  
  
  for (method in methods) {
    
    print(paste(samplename, method))
    
    atac.DB <- dba.report(atac.analyze,th = 1,method = method)
    anadata<-data.frame(mcols(atac.DB))
    df <- data.frame(chr=seqnames(atac.DB),
                     starts=start(atac.DB)-1,
                     ends=end(atac.DB),
                     names=c(rep(".", length(atac.DB))),
                     scores=c(rep(".", length(atac.DB))),
                     strands=strand(atac.DB),
                     pval=anadata$p.value,
                     FDR=anadata$FDR,
                     Fold=anadata$Fold,
                     ave=anadata$Conc
    )

    df$sign<-"0"
    df$sign[df$pval<=0.05 & df$Fold >= 1.0]<-  "1"
    df$sign[df$pval<=0.05 & df$Fold <= -1.0]<- "-1"

    p <- ggplot(data = df) +
      geom_point(shape=21,aes(x=ave,y=Fold,fill=sign,alpha=sign),color="grey90",size=2)+
      scale_fill_manual(values = c("1"="orangered3","0"="grey60","-1"="navyblue"))+
      scale_alpha_manual(values = c("1"=0.8,"0"=0.8,"-1"=0.8))+
      xlab("mean peaks signal")+ylab("fold change log2(liver/colon)")+#ggtitle(sample)+
      ylim(-2.5,3)+
      theme_bw()+
      theme(aspect.ratio = 0.5)
    
    ggsave(paste('figures/MAPLOT',samplename,method,'pdf',sep='.'), plot = p, width = 8, height = 5)
    
    # total.df<-df
    # total.df[4]<-rownames(total.df)

    pos_df<-df[df$sign=="1",]
    neg_df<-df[df$sign=="-1",]
    
    print(dim(pos_df))
    print(dim(neg_df))
    
    write.table(pos_df[,1:9],paste("bed/diffbind.",samplename,".",method,".UP.liver.bed",sep = ""),sep = "\t",quote = FALSE,col.names = FALSE,row.names = FALSE)
    write.table(pos_df[,1:4],paste("bed/diffbind.",samplename,".",method,".UP.homer.liver.bed",sep = ""),sep = "\t",quote = FALSE,col.names = FALSE,row.names = FALSE)
    write.table(neg_df[,1:9],paste("bed/diffbind.",samplename,".",method,".DOWN.colon.bed",sep = ""),sep = "\t",quote = FALSE,col.names = FALSE,row.names = FALSE)
    write.table(neg_df[,1:4],paste("bed/diffbind.",samplename,".",method,".DOWN.homer.colon.bed",sep = ""),sep = "\t",quote = FALSE,col.names = FALSE,row.names = FALSE)
    write.table(df,paste("./diffbind.","total.bed",sep = ""),sep = "\t",quote = FALSE,col.names = FALSE,row.names = FALSE)
    # total.df[4]<-rownames(total.df)
  }
}

# 
#   # PCA RPKM
#   par(mfrow=c(4,4))
#   
#   # PCA
#   pdf(paste0('figures/',samplename,'.PCA.DESEQ.RPKM.pdf'))
#   dba.plotPCA(atac.analyze,DBA_CONDITION,label = DBA_FACTOR,attributes = DBA_CONDITION,method = DBA_DESEQ2,score = DBA_SCORE_RPKM)
#   dev.off()
#   
#   pdf(paste0('figures/',samplename,'.PCA.DESEQ.SUMMIT.pdf'))
#   dba.plotPCA(atac.analyze,DBA_CONDITION,label = DBA_FACTOR,attributes = DBA_CONDITION,method = DBA_DESEQ2,score = DBA_SCORE_SUMMIT)
#   dev.off()
#   
#   # dba.plotBox(atac.analyze)
#   
#   # HEATMAP
#   pdf(paste0('figures/',samplename,'.HEATMAP.DESEQ.pdf'))
#   dba.plotHeatmap(atac.analyze,attributes = DBA_CONDITION,method = DBA_DESEQ2,score = DBA_SCORE_RPKM)
#   dev.off()
#   
#   pdf(paste0('figures/',samplename,'.HEATMAP.DESEQ_BLOCK.pdf'))
#   dba.plotHeatmap(atac.analyze,attributes = DBA_ID,method = DBA_DESEQ2_BLOCK,score = DBA_SCORE_RPKM)
#   dev.off()
#   
#   library(rtracklayer)
#   library(ggplot2)
#   
#   ## DESEQ BLOCK MA PLOT
#   
#   # export(atac.DB,"diffbind.bed",format = "bed")
#   atac.DB<-atac.DB.DESEQ2BLOCK
#   anadata<-data.frame(mcols(atac.DB))
#   df <- data.frame(chr=seqnames(atac.DB),
#                    starts=start(atac.DB)-1,
#                    ends=end(atac.DB),
#                    names=c(rep(".", length(atac.DB))),
#                    scores=c(rep(".", length(atac.DB))),
#                    strands=strand(atac.DB),
#                    pval=anadata$p.value,
#                    FDR =anadata$FDR,
#                    Fold=anadata$Fold,
#                    ave=anadata$Conc
#   )
#   
#   df$sign<-"0"
#   df$sign[df$pval<=0.05 & df$Fold> 1]<-  "1"
#   df$sign[df$pval<=0.05 & df$Fold< -1]<- "-1"
#   
#   p <- ggplot(data = df) +
#     geom_point(shape=21,aes(x=ave,y=Fold,fill=sign,alpha=sign),color="grey90",size=2)+
#     scale_fill_manual(values = c("1"="orangered3","0"="grey60","-1"="navyblue"))+
#     scale_alpha_manual(values = c("1"=0.8,"0"=0.8,"-1"=0.8))+
#     xlab("mean peaks signal")+ylab("fold change log2(colon/liver)")+#ggtitle(sample)+
#     ylim(-2.5,3)+
#     theme_bw()+
#     theme(aspect.ratio = 0.5)
#   
#   method <- 'DESEQ2_BLOCK'
#   
#   ggsave(paste0('figures/',samplename,'.MAplot.',method,'.pdf'), plot = p, width = 8, height = 5)
#   
#   total.df<-df
#   # total.df[4]<-rownames(total.df)
#   
#   pos_df<-df[df$sign=="1",]
#   neg_df<-df[df$sign=="-1",]
#   write.table(pos_df[,1:9],paste("./bed/diffbind.",samplename,".",method,".UP.bed",sep = ""),sep = "\t",quote = FALSE,col.names = FALSE,row.names = FALSE)
#   write.table(pos_df[,1:4],paste("./bed/diffbind.",samplename,".",method,".UP.homer.bed",sep = ""),sep = "\t",quote = FALSE,col.names = FALSE,row.names = FALSE)
#   write.table(neg_df[,1:9],paste("./bed/diffbind.",samplename,".",method,".DOWN.bed",sep = ""),sep = "\t",quote = FALSE,col.names = FALSE,row.names = FALSE)
#   write.table(neg_df[,1:4],paste("./bed/diffbind.",samplename,".",method,".DOWN.homer.bed",sep = ""),sep = "\t",quote = FALSE,col.names = FALSE,row.names = FALSE)
#   write.table(df,paste("./diffbind.","total.bed",sep = ""),sep = "\t",quote = FALSE,col.names = FALSE,row.names = FALSE)
#   # total.df[4]<-rownames(total.df)
#   
#   ## DESEQ BLOCK MA PLOT
#   
#   atac.DB<-atac.DB.DESEQ2
#   anadata<-data.frame(mcols(atac.DB))
#   df <- data.frame(chr=seqnames(atac.DB),
#                    starts=start(atac.DB)-1,
#                    ends=end(atac.DB),
#                    names=c(rep(".", length(atac.DB))),
#                    scores=c(rep(".", length(atac.DB))),
#                    strands=strand(atac.DB),
#                    pval=anadata$p.value,
#                    FDR =anadata$FDR,
#                    Fold=anadata$Fold,
#                    ave=anadata$Conc
#   )
#   
#   df$sign<-"0"
#   df$sign[df$pval<=0.05 & df$Fold> 1]<-  "1"
#   df$sign[df$pval<=0.05 & df$Fold< -1]<- "-1"
#   
#   p <- ggplot(data = df) +
#     geom_point(shape=21,aes(x=ave,y=Fold,fill=sign,alpha=sign),color="grey90",size=2)+
#     scale_fill_manual(values = c("1"="orangered3","0"="grey60","-1"="navyblue"))+
#     scale_alpha_manual(values = c("1"=0.8,"0"=0.8,"-1"=0.8))+
#     xlab("mean peaks signal")+ylab("fold change log2(colon/liver)")+#ggtitle(sample)+
#     ylim(-2.5,3)+
#     theme_bw()+
#     theme(aspect.ratio = 0.5)
#   
#   method <- 'DESEQ2'
#   
#   ggsave(paste0('figures/',samplename,'.MAplot.',method,'.pdf'), plot = p, width = 8, height = 5)
#   
#   total.df<-df
#   # total.df[4]<-rownames(total.df)
#   
#   pos_df<-df[df$sign=="1",]
#   neg_df<-df[df$sign=="-1",]
#   write.table(pos_df[,1:9],paste("./bed/diffbind.",samplename,".",method,".UP.bed",sep = ""),sep = "\t",quote = FALSE,col.names = FALSE,row.names = FALSE)
#   write.table(pos_df[,1:4],paste("./bed/diffbind.",samplename,".",method,".UP.homer.bed",sep = ""),sep = "\t",quote = FALSE,col.names = FALSE,row.names = FALSE)
#   write.table(neg_df[,1:9],paste("./bed/diffbind.",samplename,".",method,".DOWN.bed",sep = ""),sep = "\t",quote = FALSE,col.names = FALSE,row.names = FALSE)
#   write.table(neg_df[,1:4],paste("./bed/diffbind.",samplename,".",method,".DOWN.homer.bed",sep = ""),sep = "\t",quote = FALSE,col.names = FALSE,row.names = FALSE)
#   write.table(df,paste("./diffbind.","total.bed",sep = ""),sep = "\t",quote = FALSE,col.names = FALSE,row.names = FALSE)
# }