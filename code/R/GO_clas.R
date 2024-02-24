library(edgeR)
library(AnnotationDbi)
library(clusterProfiler)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(DESeq2)
library(DOSE)
library(stringr)
library(ggplot2)
library(ggrepel)
library(dplyr)




# c1 = read.csv('/public/home/shidetong/projects/yf/ATAC/diffbind/mac/classfilypeak/c1_gene.csv',sep = '\t',header = F)
# c2 = read.csv('/public/home/shidetong/projects/yf/ATAC/diffbind/mac/classfilypeak/c2_gene.csv',sep = '\t',header = F)
# c3 = read.csv('/public/home/shidetong/projects/yf/ATAC/diffbind/mac/classfilypeak/c3_gene.csv',sep = '\t',header = F)


c1 = read.csv('/public/home/shidetong/projects/yf/ATAC/diffbind/mac/classfilypeak/overlapwithRNA/C1gene.csv',sep = '\t',header = F)
c2 = read.csv('/public/home/shidetong/projects/yf/ATAC/diffbind/mac/classfilypeak/overlapwithRNA/C2gene.csv',sep = '\t',header = F)
c3 = read.csv('/public/home/shidetong/projects/yf/ATAC/diffbind/mac/classfilypeak/overlapwithRNA/C3gene.csv',sep = '\t',header = F)

gene = c3$V1

ego_all<-enrichGO(gene,
                 OrgDb      = org.Hs.eg.db,
                 keyType    = 'SYMBOL',
                 ont        = "ALL",
                 pAdjustMethod = "fdr",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.2)

barplot(ego_all,showCategory = 10,title="The GO_ALL enrichment analysis of all DEGs ")+
scale_size(range = c(5,10))+
scale_y_discrete(labels=function(ego_all)str_wrap(ego_all,width=30))+
theme(axis.text.y = element_text(size=15))+
theme(axis.title.x.bottom  = element_text(size = 20))+
scale_color_gradient(low = 'orange',high='red')+
scale_fill_gradient(low = 'orange',high='red')

write.table(ego_all,file = '/public/home/shidetong/projects/yf/ATAC/diffbind/mac/classfilypeak/plot/c3GO.csv',sep = '\t')
ggsave('/public/home/shidetong/projects/yf/ATAC/diffbind/mac/classfilypeak/plot/c2GO_top10.png',dpi=300, width=10, height=5)

ggsave('/public/home/shidetong/projects/yf/ATAC/diffbind/mac/classfilypeak/overlapwithRNA/c3GO_top10.png',dpi=300, width=10, height=5)

#

###特殊类型
ego_all<-enrichGO(gene,
                 OrgDb      = org.Hs.eg.db,
                 keyType    = 'SYMBOL',
                 ont        = "CC",
                 pAdjustMethod = "fdr",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.2)

barplot(ego_all,showCategory = 10,title="The GO_CC enrichment analysis of all DEGs ")+
scale_size(range = c(5,10))+
scale_y_discrete(labels=function(ego_all)str_wrap(ego_all,width=30))+
theme(axis.text.y = element_text(size=15))+
theme(axis.title.x.bottom  = element_text(size = 20))+
scale_color_gradient(low = 'orange',high='red')+
scale_fill_gradient(low = 'orange',high='red')
# scale_color_gradientn(colors = pal(20))+
# scale_fill_gradientn(colors = pal(20))

ggsave('/public/home/shidetong/projects/yf/ATAC/diffbind/mac/classfilypeak/plot/GO_color/c3GO_top10_CC.png',dpi=300, width=10, height=5)


ggplot(ego_all,aes(x=Count, y=Description))+
geom_point()
#自定义颜色

cols<-c('#E64E00','#65B48E','#F29089')
pal<-colorRampPalette(cols)
image(x=1:20,y=1,z=as.matrix(1:20),col=pal(20))

 

ego<- enrichKEGG(gene = gene ,   #基因列表(同GO) 
                # OrgDb = txdb,
                organism = "hsa",  #物种
                keyType = "KEGG",  #指定的基因ID类型，默认为kegg
                # minGSSize = 10, 
                # maxGSSize = 500,
                pvalueCutoff = 0.05,  
                pAdjustMethod = "BH",
                )
p <- barplot(eKEGG,showCategory=50)

###

##        . /public/home/zhangfy/anaconda3/envs/R/bin/R




genes<-bitr(gene, fromType="SYMBOL",  toType="ENTREZID",  OrgDb="org.Hs.eg.db")
eKEGG <- enrichKEGG(gene = genes$ENTREZID,
keyType = "kegg",
organism = "hsa",
pvalueCutoff = 0.05,
qvalueCutoff = 0.2,
pAdjustMethod = "BH")

setwd('/public/home/shidetong/projects/yf/ATAC/diffbind/mac/classfilypeak/plot/KEGG')
pdf("c3.pdf",width=6,height=3)
barplot(eKEGG,showCategory = 5,title="The KEGG enrichment analysis of all DEGs ")+
scale_size(range = c(5,10))+
scale_y_discrete(labels=function(ego_all)str_wrap(ego_all,width=30))+
theme(axis.text.y = element_text(size=15))+
theme(axis.title.x.bottom  = element_text(size = 20))+
scale_color_gradient(low = 'orange',high='red')+
scale_fill_gradient(low = 'orange',high='red')+
theme(panel.grid=element_blank())
dev.off()


png(filename ="/public/home/shidetong/projects/yf/ATAC/diffbind/mac/classfilypeak/plot/KEGG/c1_top5.png",width=500, height=250)
barplot(eKEGG,showCategory = 5,title="The KEGG enrichment analysis of all DEGs ")+
scale_size(range = c(5,10))+
scale_y_discrete(labels=function(ego_all)str_wrap(ego_all,width=30))+
theme(axis.text.y = element_text(size=15))+
theme(axis.title.x.bottom  = element_text(size = 20))+
scale_color_gradient(low = 'orange',high='red')+
scale_fill_gradient(low = 'orange',high='red')
# print(p)
dev.off()




