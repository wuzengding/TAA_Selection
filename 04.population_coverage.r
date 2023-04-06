library(dplyr)
library(readr)

# 获取外界传递的参数
args <- commandArgs(trailingOnly = TRUE)
file_a <- args[1]
file_b <- args[2]
outpath <- args[3]
TAA_gene <- args[4]
	#c("COX7B2","BRDT","CT83","CTAG1A","CLDN6","GAGE2A","IGF2BP1","MAGEA1","MAGEA2","MAGEA3","MAGEA4","MAGEA6","MAGEA9B","MAGEB2","MAGEC2","PAGE2")
cancer_type <- args[5]
	#"LUAD"

TAA_gene <- strsplit(TAA_gene,",")[[1]]
print(TAA_gene)

TPM_cutoff <- 35

if (FALSE)
	{
############注释内容
	}
for (arg in args) {
	cat (arg, "\n")
	}
# 读取文件a和文件b
ens2gene_df <- read.table(file_a, header = TRUE, na.strings = "NA", fill = TRUE, stringsAsFactors = FALSE, sep="\t")
tpm_df <- read.table(file_b, header = TRUE, stringsAsFactors = FALSE, sep="\t")

# 合并ens2gene_df和tpm_df
gene_tpm_df <- inner_join(ens2gene_df, tpm_df, by = c("gene_id" = "sample"))

# 过滤gene_name中出现在TAA_gene的行
gene_tpm_TAA_df <- gene_tpm_df %>% filter(gene_name %in% TAA_gene)

# 设置gene_name 为rowname
rownames(gene_tpm_TAA_df) <- gene_tpm_TAA_df$gene_name

# 过滤column name包含"TCGA-"的列
gene_tpm_TAA_TCGA_df <- gene_tpm_TAA_df %>% select(matches("TCGA"))

gene_tpm_TAA_TCGA_transformed_df <- 2^gene_tpm_TAA_TCGA_df -0.001

# 输出结果到指定路径
#write.csv(gene_tpm_TAA_TCGA_transformed_df, file.path(outpath, "gene_tpm_TAA_filtered.csv"),row.names = TRUE)

# 计算每个基因在肿瘤个样本中TPM≥35的样本数
N_sample <- ncol(select(gene_tpm_TAA_TCGA_transformed_df,contains("TCGA")))
N_gene <- nrow(gene_tpm_TAA_TCGA_transformed_df)

print(N_sample)
print(N_gene)

gene_tpm_TAA_TCGA_transformed_df <- gene_tpm_TAA_TCGA_transformed_df %>% 
	mutate(population_cov = round(rowSums(select(., contains("TCGA")) %>% as.matrix >= TPM_cutoff)/N_sample, 4),
		population_num = rowSums(select(., contains("TCGA")) %>% as.matrix >= TPM_cutoff))

# 根据population_cov值排序
gene_tpm_TAA_TCGA_transformed_df <- gene_tpm_TAA_TCGA_transformed_df %>%
  arrange(desc(population_cov))

#write.csv(gene_tpm_TAA_TCGA_transformed_df, file.path(outpath,"gene_tpm_TAA_TPM35_coverage2.csv"),row.names = TRUE)

# 按照基因顺序统计每个基因在前n个基因中TPM≥35的样本数，保存为新列 popu_accum_cov
i_accum_row <- rep(0, N_gene)
for(i in 1:N_gene)
	{
	j_accum_col <- rep(0,N_sample)
	for (j in 1:N_sample) 
		{
		if (any(gene_tpm_TAA_TCGA_transformed_df[1:i,j] >= TPM_cutoff))
			{
			j_accum_col[j] <- 1
			}
		}
	i_accum_row[i] <- sum(j_accum_col)
	}
gene_tpm_TAA_TCGA_transformed_df$popu_accum_cov <- i_accum_row / N_sample
gene_tpm_TAA_TCGA_transformed_df$popu_accum_num <- i_accum_row

write.csv(gene_tpm_TAA_TCGA_transformed_df, file.path(outpath, paste(cancer_type,"gene_tpm_TAA_TPM35_popu_cov.csv", sep="_")),row.names = TRUE)


### plot curve  ###

library(ggplot2)

#gene_tpm_TAA_TCGA_transformed_df <- read.table("/mnt/data2/wuzengding/02.ResDev/11.TAA_selection_human/LAML_gene_tpm_TAA_TPM35_popu_cov.csv", header = TRUE, sep = ",", row.names = 1)

gene_tpm_TAA_TCGA_transformed_df$gene_name <- row.names(gene_tpm_TAA_TCGA_transformed_df)
gene_tpm_TAA_TCGA_transformed_df$gene_name <- factor(gene_tpm_TAA_TCGA_transformed_df$gene_name, levels=rownames(gene_tpm_TAA_TCGA_transformed_df))

#绘制人群覆盖度的点图
plot_cov <- ggplot(gene_tpm_TAA_TCGA_transformed_df, aes(x = gene_name, y = population_cov)) +
	geom_point(aes(fill = "Population Coverage"), shape = 21, size = 4.5,) +
	ylab("Population Coverage") +
	scale_x_discrete(name=NULL) +
	theme(axis.text.x = element_text(face = "bold", angle=45, hjust=1))

# 绘制人群累积覆盖度的点图
plot_accum_cov <- plot_cov+
	geom_point(aes(y = popu_accum_cov, fill = "Accumulative Population Coverage"), shape = 24, size = 3, ) +
	labs(title="Population Coverage and Accumulative Coverage",face = "bold") +
	scale_fill_manual(values = c("Population Coverage" = "red", "Accumulative Population Coverage" = "blue"),
					  name = "",
					  guide = guide_legend(override.aes = list(shape = c(24, 21),size = c(3,4.5)))) +
	theme(plot.title = element_text(size=15, face="bold", hjust=0.5),
		  legend.position = "top",
		  legend.justification = "right",
		  legend.box.just = "right"
			)

ggsave(filename = paste0(outpath, paste0(cancer_type,"_gene_tpm_TAA_TPM35_popu_cov.pdf")), plot = plot_accum_cov, width = 10, height = 8)
