##########################################################################
# File Name:        03.negative_select_according_tissue.r
# Author:           wuzengding
# mail:             wuzengding@163.com
# Created Time:     Wed 13 Apr 2022 01:56:29 PM CST
##########################################################################



args <-	commandArgs()
DEG <- args[6]
MatrixCanTisTpm <- args[7]
EnsGeneBio <- args[8]
Phenotype <- args[9]
OutDir <- args[10]
Cancer <- args[11]

print(paste0("The DEG you input is :",args[6]))
print(paste0("The MatrixCanTisTpm you input is :",args[7]))
print(paste0("The EnsGeneBio you input is :",args[8]))
print(paste0("The Phenotype you input is :",args[9]))
print(paste0("The OutDir you input is :",args[10]))

library(stringr)
library(dplyr)
#########  Plot TPM boxplot of selected genes  ########
print("########## Plot TPM boxplot of selected genes ##########")
DegSelectDF <- read.csv(DEG,header=1,row.names=1,sep=",")
DegSelectDF <- DegSelectDF[-log10(DegSelectDF$padj)> 30,]
DegSelectDF <- slice_max(DegSelectDF,order_by=log2FoldChange,n=200)

EnsGeneBioDF <- read.csv(EnsGeneBio,header=1,sep="\t")
EnsGeneBioDF <- EnsGeneBioDF[which(EnsGeneBioDF$gene_name %in% rownames(DegSelectDF)),]

MatrixCanTisTpmDF <- read.csv(MatrixCanTisTpm,header=1, row.names=1,sep="\t")
#MatrixCanTisTpmDF <- MatrixCanTisTpmDF[which(rownames(MatrixCanTisTpmDF) %in% EnsGeneBioDF$gene_id),]
MatrixCanTisTpmDF <-merge(MatrixCanTisTpmDF,EnsGeneBioDF,by.x=0,by.y="gene_id")
rownames(MatrixCanTisTpmDF) <- MatrixCanTisTpmDF$gene_name
MatrixCanTisTpmDF$Row.names <- NULL
MatrixCanTisTpmDF$gene_name <- NULL
MatrixCanTisTpmDF$gene_biotype <- NULL
MatrixCanTisTpmDF <- t(MatrixCanTisTpmDF)
MatrixCanTisTpmDF <- 2^MatrixCanTisTpmDF -0.001
#write.csv(MatrixCanTisTpmDF,file=paste0(OutDir, "/MatrixCanTisTpmDF.csv"))

print("yes here1")
PhenotypeDF <- read.csv(Phenotype,header=1, sep="\t")
#write.csv(PhenotypeDF,file=paste0(OutDir, "/PhenotypeDF2.csv"))
#PhenotypeDF <- data.frame(lapply(PhenotypeDF,function(x) {gsub("-",".",x)}))
PhenotypeDF <- PhenotypeDF %>% mutate(sample = str_replace_all(sample,"-","."))
#write.csv(PhenotypeDF,file=paste0(OutDir, "/test2.csv"))
PhenotypeDF <- PhenotypeDF[which((PhenotypeDF$X_study =="GTEX") & (PhenotypeDF$X_sample_type !="Cell Line")),]
PhenotypeDF <- PhenotypeDF[,c("sample","X_primary_site")]
#write.csv(PhenotypeDF,file=paste0(OutDir, "/PhenotypeDF.csv"))

MatrixCanTisTpmDF <- merge(MatrixCanTisTpmDF,PhenotypeDF,by.x=0,by.y="sample",all.x=TRUE)
#write.csv(MatrixCanTisTpmDF,file=paste0(OutDir, "/test.csv"))
MatrixCanTisTpmDF[is.na(MatrixCanTisTpmDF$X_primary_site),c("X_primary_site")] <- Cancer
MatrixCanTisTpmDF <- MatrixCanTisTpmDF[MatrixCanTisTpmDF$X_primary_site!="",]
rownames(MatrixCanTisTpmDF) <- MatrixCanTisTpmDF$Row.names 
MatrixCanTisTpmDF$Row.names <- NULL
colnameid <- colnames(MatrixCanTisTpmDF)
colnames(MatrixCanTisTpmDF) <- str_replace_all(colnameid,"-","_")

library("ggplot2")
select_gene <-c()
for (gene in colnames(MatrixCanTisTpmDF)) {
	print(paste("gene:",gene))
	if ((gene != "X_primary_site") &&(gene != "3830403N18Rik") &&(gene != "4930503L19Rik")){
		print(paste("gene2:",gene))
		ggplot(MatrixCanTisTpmDF, aes_string(x="X_primary_site",y=gene))+theme(axis.text.x = element_text(angle=60,hjust=1))+
		xlab("")+ylab("Expression TPM")+ggtitle(paste("Expression level of",gene,"in different tissues"))+
		theme(plot.title = element_text(hjust = 0.5,size=18))+geom_boxplot()
		ggsave(paste0(OutDir,"/Expression_level_of_",gene,"_in_different_tissues.png"),width=10,height=8)
		a <- MatrixCanTisTpmDF[which(MatrixCanTisTpmDF$X_primary_site %in% c("Heart","Brain","Spleen","Liver")),gene]
		b <- MatrixCanTisTpmDF[which(MatrixCanTisTpmDF$X_primary_site %in% c("Heart","Brain","Spleen","Liver")),gene]
		if ((mean(a) < 0.5)&& mean(b) <2 ){
			select_gene <- append(select_gene,gene)
		}
	}
}
DegSelectDF <- DegSelectDF[which(rownames(DegSelectDF) %in% select_gene),]
write.csv(DegSelectDF,file=paste0(OutDir,"/",Cancer,"_select.gene.csv"))
######### prevalance  accumlate ########