import pandas as pd 

file = pd.read_csv('/public/home/shidetong/projects/yf/ATAC/diffbind/mac/dif_deseq2.txt',sep = '\t',index_col = 0)


c1 = file[(file['FDR']>0.05)&(abs(file['Fold'])<0.5)]
c2 = file[(file['FDR']<0.05)&(file['Fold']<0.5)]
c3 = file[(file['FDR']<0.05)&(file['Fold']>0.5)]

c1.to_csv('/public/home/shidetong/projects/yf/ATAC/diffbind/mac/classfilypeak/c1.bed',index = None,header = None,sep = '\t')
c2.to_csv('/public/home/shidetong/projects/yf/ATAC/diffbind/mac/classfilypeak/c2.bed',index = None,header = None,sep = '\t')
c3.to_csv('/public/home/shidetong/projects/yf/ATAC/diffbind/mac/classfilypeak/c3.bed',index = None,header = None,sep = '\t')


