.libPaths("/public/home/tangy/miniconda3/envs/r413/lib/R/library")
library(ChIPseeker)
library(ggplot2)
library(ggimage)
library("clusterProfiler")
library('ggupset')
library("org.Hs.eg.db")
library("TxDb.Mmusculus.UCSC.mm10.knownGene")
library('TxDb.Hsapiens.UCSC.hg19.knownGene')
# txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene 
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

##绘制barplot
c1 <- readPeakFile('/public/home/shidetong/projects/yf/ATAC/diffbind/mac/classfilypeak/c1.bed')
c2 <- readPeakFile('/public/home/shidetong/projects/yf/ATAC/diffbind/mac/classfilypeak/c2.bed')
c3 <- readPeakFile('/public/home/shidetong/projects/yf/ATAC/diffbind/mac/classfilypeak/c3.bed')

peak=GenomicRanges::GRangesList(C1= c1,C2 = c2,C3 = c3)

peakAnnoList <- lapply(peak, annotatePeak, 
                       TxDb=txdb,tssRegion=c(-3000, 3000))
p <- plotAnnoBar(peakAnnoList)+
theme(axis.text.y = element_text(size = 25),axis.text.x = element_text(size = 20))+
theme(axis.title = element_text(size = 25))+
theme(legend.text = element_text(size = 20))

ggsave('/public/home/shidetong/projects/yf/ATAC/diffbind/mac/classfilypeak/plot/clas_bar_test.png',dpi=300, width=10, height=5)

png( 
    filename = "/public/home/shidetong/projects/yf/ATAC/diffbind/mac/classfilypeak/plot/clas_bar.png", # 文件名称
    width = 640,           # 宽
    height = 320,          # 高
    units = "px",          # 单位
    bg = "white",          # 背景颜色
    res = 72)  
p <- plotAnnoBar(peakAnnoList)+
theme(axis.text.y = element_text(size = 25),axis.text.x = element_text(size = 20))+
theme(axis.title = element_text(size = 25))+
theme(legend.text = element_text(size = 25))
dev.off()

##----------------------------------------------------------
peak=GenomicRanges::GRangesList(CBX6=M,CBX7=P)
covplot(peak)
covplot(peak, weightCol="V5") + facet_grid(chr ~ .id)
covplot(peak, weightCol="V5", chrs=c("chr17", "chr18"), 
        xlim=c(4e7, 5e7)) + facet_grid(chr ~ .id)

x = annotatePeak(M, tssRegion=c(-1000, 1000), TxDb=txdb)
x
p =as.GRanges(x)
as.GRanges(x) #%>% head(10)
tmp=as.data.frame(x)
x2 = annotatePeak(M, tssRegion=c(-1000, 1000), TxDb=txdb, addFlankGeneInfo=TRUE, flankDistance=5000)

peakAnno <- annotatePeak(M, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Mm.eg.db",level="gene")         #针对比对出的ChIP数据peak，注释在peak附近的基因信息
tt = as.GRanges(peakAnno)
write.table(tt,'tt.csv',sep ='\t',row.names = FALSE)
x3 = annotatePeak(M, tssRegion=c(-1000, 1000), TxDb=txdb, 
                  addFlankGeneInfo=TRUE, flankDistance=5000,
                  annoDb = "org.Mm.eg.db")
tmp3=as.data.frame(x3)

peakHeatmap(P, weightCol="V5", TxDb=txdb, 
            upstream=3000, downstream=3000, 
            color=rainbow(length(M)))
###            
promoter <- getPromoters(TxDb=txdb, 
                  upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(f, 
                          windows=promoter)
tagHeatmap(tagMatrix, xlim=c(-3000, 3000), 
           color="red")
###

plotAvgProf2(M, TxDb=txdb, 
             upstream=3000, downstream=3000,
             xlab="Genomic Region (5'->3')", 
             ylab = "Read Count Frequency",
             conf = 0.95, resample = 1000)
#饼图
peakAnno <- annotatePeak(peak, 
                         tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Mm.eg.db")
#饼图

plotAnnoPie(peakAnno)
plotAnnoPie(peakAnno,ndigit = 2,cex = 0.9,legend.position = 'rightside')
#柱状图
plotAnnoBar(peakAnno,ndigit = 2,cex = 1.2)


promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tag <- getTagMatrix(peak,windows = promoter)
x3 = annotatePeak(peak, tssRegion=c(-1000, 1000), TxDb=txdb, 
                  addFlankGeneInfo=TRUE, flankDistance=5000,
                  annoDb = "org.Mm.eg.db",verbose = FALSE)

p =as.GRanges(peakAnno)
# as.GRanges(p)
tmp=as.data.frame(p)
gene = tmp$SYMBOL
ge = unique(gene)
ge = as.data.frame(ge)

ego <- enrichGO(gene = ge, 
                    keyType = "SYMBOL", 
                    OrgDb = org.Mm.eg.db, 
                    ont = "BP", 
                    pAdjustMethod = "BH", 
                    qvalueCutoff = 0.05, 
                    readable = TRUE)

gene <- seq2gene(peak, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)




####
BCF1_B <- readPeakFile('/public/home/shidetong/projects/sdt/chip-seq/haplotype/B6_M/mac/Maternal_B6_peaks.narrowPeak')
BCF1_C <- readPeakFile('/public/home/shidetong/projects/sdt/chip-seq/haplotype/B6_M/mac/Paternal_B6_peaks.narrowPeak')
CBF1_C <- readPeakFile('/public/home/shidetong/projects/sdt/chip-seq/haplotype/Cast_M/mac/Maternal_Cast_peaks.narrowPeak')
CBF1_B <- readPeakFile('/public/home/shidetong/projects/sdt/chip-seq/haplotype/Cast_M/mac/Paternal_Cast_peaks.narrowPeak')

BCF1 <- readPeakFile('/public/home/shidetong/projects/sdt/chip-seq/Parental_transmission/B6_batch2_Rep1/mac/B6_batch2_Chip_CTCF_1_peaks.narrowPeak')
CBF1 <- readPeakFile('/public/home/shidetong/projects/sdt/chip-seq/Parental_transmission/Cast_batch2_Rep1/mac/Cast_batch2_Chip_CTCF_1_peaks.narrowPeak')

peak <- CBF1

peakAnno <- annotatePeak(peak, 
                         tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Mm.eg.db")

p =as.GRanges(peakAnno)
tmp=as.data.frame(p)
gene = tmp$SYMBOL
ge = unique(gene)
g = as.data.frame(ge)

write.table(g,'/public/home/shidetong/projects/sdt/chip-seq/haplotype/chipseeker/CBF1.csv',row.names=  F,sep = '\t',col.names = F,quote=F)

 


setwd('/public/home/shidetong/projects/sdt/chip-seq/haplotype/chipseeker/')
png( 
    filename = "CBF1_pie.png", # 文件名称
    width = 480,           # 宽
    height = 480,          # 高
    units = "px",          # 单位
    bg = "white",          # 背景颜色
    res = 72)  
plotAnnoPie(peakAnno)
dev.off()


png( 
    filename = "CBF1_bar.png", # 文件名称
    width = 480,           # 宽
    height = 480,          # 高
    units = "px",          # 单位
    bg = "white",          # 背景颜色
    res = 72)  
plotAnnoBar(peakAnno)
dev.off()


png( 
    filename = "CBF1_veen.png", # 文件名称
    width = 480,           # 宽
    height = 480,          # 高
    units = "px",          # 单位
    bg = "white",          # 背景颜色
    res = 72)  
vennpie(peakAnno)
dev.off()


png( 
    filename = "BCF1_B_up.png", # 文件名称
    width = 480,           # 宽
    height = 480,          # 高
    units = "px",          # 单位
    bg = "white",          # 背景颜色
    res = 72)  
upsetplot(peakAnno, vennpie=TRUE)
dev.off()

