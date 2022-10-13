##########################################################################
# File Name:        02.different_expression_select.r
# Author:           wuzengding
# mail:             wuzengding@163.com
# Created Time:     Mon 11 Apr 2022 02:15:19 PM CST
##########################################################################


args <-	commandArgs()
BioType <- args[6]
MatrixCount <- args[7]
MatrixTpm <- args[8]
EnsGeneBio <- args[9]
OutDir <- args[10]

print(paste0("The BioType you input is :",args[6]))
print("\n")


#########Select genes by BioType ########
print("##########Select genes by BioType ##########")
BioTypeItem <- c("lncRNA","protein_coding","miRNA","processed_pseudogene","nonsense_mediated_decay","unprocessed_pseudogene","processed_transcript","TEC","snRNA","misc_RNA","retained_intron","transcribed_unprocessed_pseudogene","IG_V_gene","rRNA_pseudogene","IG_J_gene","IG_C_pseudogene","polymorphic_pseudogene","IG_V_pseudogene","snoRNA","IG_D_gene","IG_C_gene","TR_V_gene","rRNA","IG_J_pseudogene","TR_J_gene","transcribed_processed_pseudogene","unitary_pseudogene","vault_RNA","scaRNA","transcribed_unitary_pseudogene","TR_V_pseudogene","ribozyme","TR_J_pseudogene","translated_unprocessed_pseudogene","Mt_tRNA","non_stop_decay","TR_D_gene","Mt_rRNA","TR_C_gene","scRNA")
if (BioType %in% BioTypeItem){
	print(paste("You select all genes whose transcript biotype is",BioType))
}else{
	print("The transcript biotype you input is illegal,you input biotype in the follow list:")
	print(paste("lncRNA","protein_coding","miRNA","processed_pseudogene","nonsense_mediated_decay","unprocessed_pseudogene","processed_transcript","TEC","snRNA","misc_RNA","retained_intron","transcribed_unprocessed_pseudogene","IG_V_gene","rRNA_pseudogene","IG_J_gene","IG_C_pseudogene","polymorphic_pseudogene","IG_V_pseudogene","snoRNA","IG_D_gene","IG_C_gene","TR_V_gene","rRNA","IG_J_pseudogene","TR_J_gene","transcribed_processed_pseudogene","unitary_pseudogene","vault_RNA","scaRNA","transcribed_unitary_pseudogene","TR_V_pseudogene","ribozyme","TR_J_pseudogene","translated_unprocessed_pseudogene","Mt_tRNA","non_stop_decay","TR_D_gene","Mt_rRNA","TR_C_gene","scRNA",sep=","))
	quit (0)
}

EnsGeneBioDF <- read.csv(EnsGeneBio,header=1,sep="\t")
EnsGeneBioDF <- EnsGeneBioDF[EnsGeneBioDF$Transcript_Biotype==BioType,]
#write.csv(EnsGeneBioDF,file=paste0(OutDir, "/nsemblbiodf.test.txt"))

#########Different Expression Gene Select ################

library(ggplot2)
library(dplyr)
if (TRUE){
library(DESeq2)
library(tidyr)


print("##########DESeq2 analysis##########")
MatrixCountDF <- read.csv(MatrixCount,header=1,row.names=1,sep="\t")
#print(EnsGeneBioDF$ENSEMBL)
#print("gap")
MatrixCountDF <- MatrixCountDF[which(rownames(MatrixCountDF) %in% EnsGeneBioDF$GENE),]
#write.csv(MatrixCountDF,file=paste0(OutDir, "/test.txt"))
MatrixCountDF <- 2^MatrixCountDF -1
MatrixCountCts <- as.matrix(MatrixCountDF)
#write.csv(MatrixCountCts,file=paste0(OutDir, "/MatrixCountDF.test.txt"))
SampAnnoDF <- data.frame(SampleID=colnames(MatrixCountDF))
SampAnnoDF <- SampAnnoDF %>% separate(SampleID,c("Cohort"))
SampAnnoDF <- data.frame(SampleID=colnames(MatrixCountDF),Condition=SampAnnoDF$Cohort)
Nsample <- nrow(SampAnnoDF)
#write.csv(SampAnnoDF,file=paste0(OutDir, "/SampAnnoDF.test.txt"))

DDS <- DESeqDataSetFromMatrix(countData=round(MatrixCountCts), colData=SampAnnoDF, design=~ Condition)
keep <- rowSums(counts(DDS)) >= Nsample*5
DDS <- DDS[keep,]
DDS$Condition <- relevel(DDS$Condition, ref="GTEX")
DDS <- DESeq(DDS)
Result <- results(DDS)
ResultAnno <- Result[order(Result$log2FoldChange, decreasing = c(TRUE)), ]
ResultAnno[which(ResultAnno$log2FoldChange >= 1 & ResultAnno$padj < 0.01),'sig'] <- 'up'
ResultAnno[which(ResultAnno$log2FoldChange <= -1 & ResultAnno$padj < 0.01),'sig'] <- 'down'
ResultAnno[which(abs(ResultAnno$log2FoldChange) <= 1 | ResultAnno$padj >= 0.01),'sig'] <- 'none'
write.csv(ResultAnno,file=paste0(OutDir, "/DESeq2.DEG.",BioType,".csv"))
}

######### DGE select result exhibition by ploting ##############
print("##########Volcano  and Heatmap plot ############")
ResultAnno <- read.csv(file=paste0(OutDir, "/DESeq2.DEG.",BioType,".csv"),row.names=1,header=1)
p_volcano <- ggplot(data = ResultAnno, aes(x = log2FoldChange, y = -log10(padj), color = sig))+
	geom_point(size = 1) + 
	scale_color_manual(values = c('red', 'gray', 'green'), limits = c('up', 'none', 'down'))+
	labs(x = 'log2 Fold Change', y = '-log10 adjust p-value', title = 'control vs treat', color = '')+
	theme(plot.title = element_text(hjust = 0.5, size = 14), panel.grid = element_blank(),panel.background = element_rect(color = 'black', fill = 'transparent'), legend.key = element_rect(fill = 'transparent'))+
	geom_vline(xintercept = c(-1, 1), lty = 3, color = 'black')+
	geom_hline(yintercept = 2, lty = 3, color = 'black')
ggsave(paste0(OutDir, "/DESeq2.DEG.",BioType,".Volcano.pdf"), plot = p_volcano,width=10,height=8)

library(pheatmap)
DegSelect <- ResultAnno[-log10(ResultAnno$padj)> 30,]
DegSelect <- slice_max(DegSelect,order_by=log2FoldChange,n=200)
EnsGeneList <- EnsGeneBioDF[which(EnsGeneBioDF$GENE %in% rownames(DegSelect)),]
#write.csv(EnsGeneList,file=paste0(OutDir, "/TEST2.csv"))
MatrixTpmDF <- read.csv(MatrixTpm,header=1,row.names=1,sep="\t")
MatrixTpmSelectDF <- MatrixTpmDF[which(rownames(MatrixTpmDF) %in% EnsGeneList$ENSEMBL),]
MatrixTpmSelectmergeDF <-merge(MatrixTpmSelectDF,EnsGeneList,by.x=0,by.y="ENSEMBL")
#write.csv(MatrixTpmSelectmergeDF,file=paste0(OutDir, "/TEST.csv"))
rownames(MatrixTpmSelectmergeDF) <- MatrixTpmSelectmergeDF$GENE
MatrixTpmSelectmergeDF$Row.names <- NULL
MatrixTpmSelectmergeDF$GENE <- NULL
MatrixTpmSelectmergeDF$Transcript_Biotype <- NULL
#write.csv(MatrixTpmSelectmergeDF,file=paste0(OutDir, "/TEST3.csv"))
p_heatmap <- pheatmap(MatrixTpmSelectmergeDF,cluster_rows= FALSE,cluster_cols=FALSE,fontsize_col=2,
		xlab='', ylab="", fontsize_row=2,cellheight=2,main= "TPM heatmap gene of top200  Cancer vs Control")
ggsave(paste0(OutDir,"/DESeq2.DEG.",BioType,".TPM.Heatmap.gene.of.top200.pdf"),plot = p_heatmap,width=10,height=8)
dev.off()

DegSelect <- ResultAnno[-log10(ResultAnno$padj)> 30,]
DegSelect <- slice_max(DegSelect,order_by=log2FoldChange,n=100)
EnsGeneList <- EnsGeneBioDF[which(EnsGeneBioDF$GENE %in% rownames(DegSelect)),]
MatrixTpmDF <- read.csv(MatrixTpm,header=1,row.names=1,sep="\t")
MatrixTpmSelectDF <- MatrixTpmDF[which(rownames(MatrixTpmDF) %in% EnsGeneList$ENSEMBL),]
MatrixTpmSelectmergeDF <-merge(MatrixTpmSelectDF,EnsGeneList,by.x=0,by.y="ENSEMBL")
rownames(MatrixTpmSelectmergeDF) <- MatrixTpmSelectmergeDF$GENE
MatrixTpmSelectmergeDF$Row.names <- NULL
MatrixTpmSelectmergeDF$GENE <- NULL
MatrixTpmSelectmergeDF$Transcript_Biotype <- NULL
p_heatmap <- pheatmap(MatrixTpmSelectmergeDF,cluster_rows= FALSE,cluster_cols=FALSE,fontsize_col=2,
		xlab='', ylab="", fontsize_row=3,cellheight=2,main= "TPM heatmap gene of top100  Cancer vs Control")
ggsave(paste0(OutDir,"/DESeq2.DEG.",BioType,".TPM.Heatmap.gene.of.top100.pdf"),plot = p_heatmap,width=10,height=8)
dev.off()
