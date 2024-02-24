library(AnnotationDbi)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DESeq2)
library(DOSE)
library(stringr)
library(AnnotationDbi)
library(ggplot2)
library(pheatmap)

gene1 = read.table('/public/home/shidetong/projects/yf/hic/HiCHap_workspace/Matrix/pc_gene/gene_name/HvsP_A2B.csv')
gene2 = read.table('/public/home/shidetong/projects/yf/hic/HiCHap_workspace/Matrix/pc_gene/gene_name/HvsB_A2B.csv')
gene <- merge(gene1,gene2,by='V1')
gene3 = read.table('/public/home/shidetong/projects/yf/hic/HiCHap_workspace/Matrix/pc_gene/gene_name/HvsP_B2A.csv')
gene4 = read.table('/public/home/shidetong/projects/yf/hic/HiCHap_workspace/Matrix/pc_gene/gene_name/HvsB_B2A.csv')
gene <- merge(gene3,gene4,by='V1')

A_B = read.table('/public/home/shidetong/projects/yf/hic/HiCHap_workspace/Matrix/pc_gene/gene_name/overlapRNAdif/hig.csv')
B_A = read.table('/public/home/shidetong/projects/yf/hic/HiCHap_workspace/Matrix/pc_gene/gene_name/overlapRNAdif/low.csv')
ego_all<-enrichGO(gene=A_B,
                 OrgDb      = org.Hs.eg.db,
                 keyType    = 'SYMBOL',
                 ont        = "ALL",
                 pAdjustMethod = "fdr",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.2)
barplot(ego_all,showCategory = 10,title="The GO_ALL enrichment analysis of all DEGs ")+
scale_size(range = c(5,20))+
scale_y_discrete(labels=function(ego_all)str_wrap(ego_all,width=5))
ggsave('/public/home/shidetong/projects/yf/hic/HiCHap_workspace/Matrix/pc_gene/gene_name/AB.png')


ggsave('/public/home/shidetong/projects/yf/hic/HiCHap_workspace/Matrix/pc_gene/gene_name/HvsB_A2B.png')


write.table(gene,file='/public/home/shidetong/projects/yf/hic/HiCHap_workspace/Matrix/pc_gene/gene_name/B2A.csv',row.names = F, quote = F, sep="\t")

