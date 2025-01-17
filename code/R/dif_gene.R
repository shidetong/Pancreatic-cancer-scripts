gene <- read.table('/public/home/shidetong/projects/yf/Plot/Signal_plot/Gene_list.csv')

BvsH_B <- read.csv('/public/home/shidetong/projects/yf/RNA-seq/20200924/dif/BvsH_B_dif.csv',sep = '\t')
BvsH_H <- read.csv('/public/home/shidetong/projects/yf/RNA-seq/20200924/dif/BvsH_H_dif.csv',sep = '\t')
PvsH_P <- read.csv('/public/home/shidetong/projects/yf/RNA-seq/20200924/dif/PvsH_P_dif.csv',sep = '\t')
PvsH_H <- read.csv('/public/home/shidetong/projects/yf/RNA-seq/20200924/dif/PvsH_H_dif.csv',sep = '\t')

BvsH_B <- rownames(BvsH_B)
BvsH_B <- as.data.frame(BvsH_B)
BvsH_H <- rownames(BvsH_H)
BvsH_H <- as.data.frame(BvsH_H)
PvsH_P <- rownames(PvsH_P)
PvsH_P <- as.data.frame(PvsH_P)
PvsH_H <- rownames(PvsH_H)
PvsH_H <- as.data.frame(PvsH_H)

BH_B <- intersect(gene$V1,BvsH_B$BvsH_B)
BH_H <- intersect(gene$V1,BvsH_H$BvsH_H)
PH_P <- intersect(gene$V1,PvsH_P$PvsH_P)
PH_H <- intersect(gene$V1,PvsH_H$PvsH_H)

setwd('/public/home/shidetong/projects/yf/RNA-seq/20200924/dif/overlap')
write.table(BH_B,'BH_B.csv',sep = '\t',row.names = F,col.names = F)
write.table(BH_H,'BH_H.csv',sep = '\t',row.names = F,col.names = F)
write.table(PH_P,'PH_P.csv',sep = '\t',row.names = F,col.names = F)
write.table(PH_H,'PH_H.csv',sep = '\t',row.names = F,col.names = F)