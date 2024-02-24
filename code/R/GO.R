setwd("/public/home/shidetong/projects/yf/RNA-seq/20200924/count")
options(stringsAsFactors = FALSE)
Bxpc1<-read.table("Bxpc-1.count",sep = "\t",col.names = c("gene_id","Bxpc_1"))
Bxpc2<-read.table("Bxpc-1.count",sep = "\t",col.names = c("gene_id","Bxpc_2"))
Bxpc3<-read.table("Bxpc-2.count",sep = "\t",col.names = c("gene_id","Bxpc_3"))
H6C71<-read.table("H6C7-1.count",sep = "\t",col.names = c("gene_id","H6C7_1"))
H6C72<-read.table("H6C7-2.count",sep = "\t",col.names = c("gene_id","H6C7_2"))
H6C73<-read.table("H6C7-3.count",sep = "\t",col.names = c("gene_id","H6C7_3"))
Panc1<-read.table("Panc-1.count",sep = "\t",col.names = c("gene_id","Panc_1"))
Panc2<-read.table("Panc-2.count",sep = "\t",col.names = c("gene_id","Panc_2"))
Panc3<-read.table("Panc-3.count",sep = "\t",col.names = c("gene_id","Panc_3"))



raw_count1 <- merge(merge(Bxpc1, Bxpc2, by="gene_id"), merge(Bxpc3, H6C71, by="gene_id"))
raw_count2 <- merge(merge(H6C72, H6C73, by="gene_id"), merge(Panc1, Panc2, by="gene_id"))
raw_count3 <- merge(raw_count1,raw_count2,by='gene_id')
raw_count4 <- merge(raw_count3,Panc3,by='gene_id')


head(raw_count4)
tail(raw_count4)
raw_count_filt <- raw_count4[-1:-5,]
head(raw_count_filt)
ENSEMBL <- gsub("\\.\\d*\\_\\d*", "", raw_count_filt$gene_id) 


row.names(raw_count_filt) <- ENSEMBL
head(raw_count_filt)


write.csv(raw_count_filt,"/public/home/shidetong/projects/yf/RNA-seq/20200924/count/all.count")


##GO分析
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)

sig.gene<-raw_count_filt
head(sig.gene)
gene<-sig.gene[,1]
head(gene)
gene.df<-bitr(gene, fromType = "ENSEMBL", 
              toType = c("SYMBOL","ENTREZID"),
              OrgDb = org.Hs.eg.db)

head(gene.df)

ego_cc<-enrichGO(gene       = gene.df$ENSEMBL,
                 OrgDb      = org.Hs.eg.db,
                 keyType    = 'ENSEMBL',
                 ont        = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.01,
                 qvalueCutoff = 0.05)
ego_bp<-enrichGO(gene       = gene.df$ENSEMBL,
                 OrgDb      = org.Hs.eg.db,
                 keyType    = 'ENSEMBL',
                 ont        = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.01,
                 qvalueCutoff = 0.05)
barplot(ego_bp,showCategory = 18,title="The GO_BP enrichment analysis of all DEGs ")+ 
  scale_size(range=c(2, 12))+
  scale_x_discrete(labels=function(ego_bp) str_wrap(ego_bp,width = 25))
dotplot(ego_bp)
enrichMap(ego_bp)
plotGOgraph(ego_bp)


library(stringr)
kk<-enrichKEGG(gene=gene.df$ENTREZID,
               organism = 'hsa',
               pvalueCutoff = 0.05)
kk[1:30]
barplot(kk,showCategory = 25, title="The KEGG enrichment analysis of all DEGs")+
  scale_size(range=c(2, 12))+
  scale_x_discrete(labels=function(kk) str_wrap(kk,width = 25))
dotplot(kk,showCategory = 25, title="The KEGG enrichment analysis of all DEGs")+
  scale_size(range=c(2, 12))+
  scale_x_discrete(labels=function(kk) str_wrap(kk,width = 25))

browseKEGG(kk, 'hsa04974')







