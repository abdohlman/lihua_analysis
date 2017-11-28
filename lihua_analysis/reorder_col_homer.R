bed <- as.data.frame(read.table("/Users/ergangwang/Downloads/Lihua_ATACSeq/lihua_analysis/bed_liver.colon/diffbind.allsamples.edgeRlm.DOWN.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
bed <- bed[, c(1, 2, 3, 4, 9, 6, 7, 8, 4)]
bed$V4 = paste(bed$V2, bed$V3, sep = "-")
bed$V4 = paste(bed$V1, bed$V4, sep = ":")
write.table(bed, "/Users/ergangwang/Downloads/Lihua_ATACSeq/lihua_analysis/bed_liver.colon/diffbind.allsamples.edgeRlm.DOWN.reordername.bed", sep="\t", row.name = F, col.name = F, quote=FALSE)

bed2 <- as.data.frame(read.table("/Users/ergangwang/Downloads/Lihua_ATACSeq/lihua_analysis/bed_liver.colon/diffbind.allsamples.edgeRlm.UP.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
bed2 <- bed2[, c(1, 2, 3, 4, 9, 6, 7, 8, 4)]
bed2$V4 = paste(bed2$V2, bed2$V3, sep = "-")
bed2$V4 = paste(bed2$V1, bed2$V4, sep = ":")
write.table(bed2, "/Users/ergangwang/Downloads/Lihua_ATACSeq/lihua_analysis/bed_liver.colon/diffbind.allsamples.edgeRlm.UP.reordername.bed", sep="\t", row.name = F, col.name = F, quote=FALSE)