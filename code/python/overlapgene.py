import pandas as pd 

c1 = pd.read_csv('/public/home/shidetong/projects/yf/ATAC/diffbind/mac/classfilypeak/c1_gene.csv',sep = '\t',header = None)
c2 = pd.read_csv('/public/home/shidetong/projects/yf/ATAC/diffbind/mac/classfilypeak/c2_gene.csv',sep = '\t',header = None)
c3 = pd.read_csv('/public/home/shidetong/projects/yf/ATAC/diffbind/mac/classfilypeak/c3_gene.csv',sep = '\t',header = None)


nodif_gene = pd.read_csv('/public/home/shidetong/projects/yf/RNA-seq/20200924/dif/nodiff_gene.csv',sep ='\t')
nodif_gene[0] =nodif_gene.index
symbol = nodif_gene[[0]]
symbol.index=range(symbol.shape[0])

Hhig =pd.read_csv('/public/home/shidetong/projects/yf/RNA-seq/20200924/dif/overlap/H_hig.csv',sep = '\t',header = None)
Hlow = pd.read_csv('/public/home/shidetong/projects/yf/RNA-seq/20200924/dif/overlap/H_low.csv',sep = '\t',header = None)



C1 = pd.merge(c1,symbol,on = [0],how = 'inner')
C2 = pd.merge(Hhig,c2,on =[0],how = 'inner')
C3 = pd.merge(Hlow,c3,on = [0],how = 'inner')


C1.to_csv('/public/home/shidetong/projects/yf/ATAC/diffbind/mac/classfilypeak/overlapwithRNA/C1gene.csv',sep = '\t',header = None,index  =None)
C2.to_csv('/public/home/shidetong/projects/yf/ATAC/diffbind/mac/classfilypeak/overlapwithRNA/C2gene.csv',sep = '\t',header = None,index  =None)
C3.to_csv('/public/home/shidetong/projects/yf/ATAC/diffbind/mac/classfilypeak/overlapwithRNA/C3gene.csv',sep = '\t',header = None,index  =None)


