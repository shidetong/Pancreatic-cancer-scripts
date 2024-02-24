rm(list=ls())
setwd("/public/home/shidetong/projects/yf/RNA-seq/20200924/ballgown/")
gene_count <- read.csv("gene_count_matrix.csv",stringsAsFactors = F)
library(stringr)
gene_count_1<-rep(NA,nrow(gene_count))
for (i in 1:nrow(gene_count)){
  gene_count_1[i] <- unlist(str_split(gene_count[i,1], pattern = "\\|"))[1]
}
gene_count$gene_id <- gene_count_1
rownames(gene_count) <- gene_count[,1]
gene_count <-gene_count[,-1]
gene_count_group_1 <- gene_count[, 1:6]
gene_count_group_2 <- gene_count[, c(4,5,6, 7, 8,9)]
gene_count_group_1
gene_count_group_2
write.csv(gene_count_group_1, file = "/public/home/shidetong/projects/yf/RNA-seq/20200924/ballgown/gene_count_group_1.csv", row.names = TRUE)
write.csv(gene_count_group_2, file = "/public/home/shidetong/projects/yf/RNA-seq/20200924/ballgown/gene_count_group_2.csv", row.names = TRUE)
library(DESeq2)
condition_group_1 <- factor(c(rep("Bxpc", 3), rep("H6C7", 3)), levels = c("Bxpc","H6C7"))
condition_group_2 <- factor(c(rep("H6C7", 3), rep("Panc", 3)), levels = c("H6C7","Panc"))
condition_group_1
condition_group_2
colData_group_1 <- data.frame(row.names = colnames(gene_count_group_1), condition_group_1)
colData_group_2 <- data.frame(row.names = colnames(gene_count_group_2), condition_group_2)
colData_group_1
colData_group_2
dds_group_1 <- DESeqDataSetFromMatrix(gene_count_group_1, colData_group_1, design = ~condition_group_1)
dds_group_2 <- DESeqDataSetFromMatrix(gene_count_group_2, colData_group_2, design = ~condition_group_2)
dds_group_1 <- DESeq(dds_group_1)
dds_group_2 <- DESeq(dds_group_2)
dds_group_1
dds_group_2
res_group_1 = results(dds_group_1, contrast=c("condition_group_1","Bxpc","H6C7"))
res_group_2 = results(dds_group_2, contrast=c("condition_group_2","Panc","H6C7"))
res_group_1 = res_group_1[order(res_group_1$pvalue),]
res_group_2 = res_group_2[order(res_group_2$pvalue),]
head(res_group_1)
head(res_group_2)
summary(res_group_1)
summary(res_group_2)
write.csv(res_group_1, file="/public/home/shidetong/projects/yf/RNA-seq/20200924/ballgown/all_different_genes_group_1_genecount.csv")
write.csv(res_group_2, file="/public/home/shidetong/projects/yf/RNA-seq/20200924/ballgown/all_different_genes_group_2_genecount.csv")
table(res_group_1$pvalue<0.05)
table(res_group_2$pvalue<0.05)
significant_pvalue_different_genes_group_1 <- subset(res_group_1, pvalue < 0.05 & abs(log2FoldChange) > 0.584963)
significant_pvalue_different_genes_group_2 <- subset(res_group_2, pvalue < 0.05 & abs(log2FoldChange) > 0.584963)
dim(significant_pvalue_different_genes_group_1)
dim(significant_pvalue_different_genes_group_2)
head(significant_pvalue_different_genes_group_1)
head(significant_pvalue_different_genes_group_2)
write.csv(significant_pvalue_different_genes_group_1, file = "/public/home/shidetong/projects/yf/RNA-seq/20200924/ballgown/significant_pvalue_different_genes_group_1_genecount.csv")
write.csv(significant_pvalue_different_genes_group_2, file = "/public/home/shidetong/projects/yf/RNA-seq/20200924/ballgown/significant_pvalue_different_genes_group_2_genecount.csv")
