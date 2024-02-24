import pandas as pd 

h6c7_1 = pd.read_csv('/public/home/shidetong/projects/yf/RNA-seq/20200924/FPKM/H6C7_1.txt',sep = '\t',header = None)
h6c7_2 = pd.read_csv('/public/home/shidetong/projects/yf/RNA-seq/20200924/FPKM/H6C7_2.txt',sep = '\t',header = None)
h6c7_3 = pd.read_csv('/public/home/shidetong/projects/yf/RNA-seq/20200924/FPKM/H6C7_3.txt',sep = '\t',header = None)
bxpc_1 = pd.read_csv('/public/home/shidetong/projects/yf/RNA-seq/20200924/FPKM/Bxpc_1.txt',sep = '\t',header = None)
bxpc_2 = pd.read_csv('/public/home/shidetong/projects/yf/RNA-seq/20200924/FPKM/Bxpc_2.txt',sep = '\t',header = None)
bxpc_3 = pd.read_csv('/public/home/shidetong/projects/yf/RNA-seq/20200924/FPKM/Bxpc_3.txt',sep = '\t',header = None)
panc_1 = pd.read_csv('/public/home/shidetong/projects/yf/RNA-seq/20200924/FPKM/Panc_1.txt',sep = '\t',header = None)
panc_2 = pd.read_csv('/public/home/shidetong/projects/yf/RNA-seq/20200924/FPKM/Panc_2.txt',sep = '\t',header = None)
panc_3 = pd.read_csv('/public/home/shidetong/projects/yf/RNA-seq/20200924/FPKM/Panc_3.txt',sep = '\t',header = None)

h6c7_1 = h6c7_1.drop_duplicates()
h6c7_2 = h6c7_2.drop_duplicates()
h6c7_3 = h6c7_3.drop_duplicates()
bxpc_1 = bxpc_1.drop_duplicates()
bxpc_2 = bxpc_2.drop_duplicates()
bxpc_3 = bxpc_3.drop_duplicates()
panc_1 = panc_1.drop_duplicates()
panc_2 = panc_2.drop_duplicates()
panc_3 = panc_3.drop_duplicates()

f1 = pd.merge(h6c7_1,h6c7_2,on=[0],how = 'inner')
f1.drop_duplicates(subset=0,keep='first',inplace=True)
f2 = pd.merge(h6c7_3,bxpc_1,on=[0],how = 'inner')
f2.drop_duplicates(subset=0,keep='first',inplace=True)
f3 = pd.merge(bxpc_2,bxpc_3,on=[0],how = 'inner')
f3.drop_duplicates(subset=0,keep='first',inplace=True)
f4 = pd.merge(panc_1,panc_2,on = [0],how = 'inner')
f4.drop_duplicates(subset=0,keep='first',inplace=True)

f5 = pd.merge(f1,f2,on = [0],how = 'inner')
f5.drop_duplicates(subset=0,keep='first',inplace=True)
f6 = pd.merge(f3,f4,on = [0],how = 'inner')
f6.drop_duplicates(subset=0,keep='first',inplace=True)
f7 = pd.merge(f5,f6,on = [0],how = 'inner')
f7.drop_duplicates(subset=0,keep='first',inplace=True)
f8 = pd.merge(f7,panc_3,on = [0],how = 'inner')
f8.drop_duplicates(subset=0,keep='first',inplace=True)

f8.columns = ['symbol','h6c71','h6c72','h6c73','bxpc1','bxpc2','bxpc3','panc1','panc2','panc3']

df= f8[(f8['h6c71'] != 0) | (f8['h6c72'] != 0)|(f8['h6c73'] != 0)|
(f8['bxpc1'] != 0)|(f8['bxpc2'] != 0)|(f8['bxpc3'] != 0)|
(f8['panc1'] != 0)|(f8['panc2'] != 0)|(f8['panc3'] != 0)]


df['H6C7']=(df['h6c71']+df['h6c72']+df['h6c73'])/3
df['BXPC']=(df['bxpc1']+df['bxpc2']+df['bxpc3'])/3
df['PANC']=(df['panc1']+df['panc2']+df['panc3'])/3
df = df.loc[:,['symbol','H6C7','BXPC','PANC']]
df.to_csv('/public/home/shidetong/projects/yf/RNA-seq/20200924/FPKM/fpkm_all_sig.csv',sep = '\t',index = None)