import pandas as pd 

df = pd.read_csv('/public/home/shidetong/projects/yf/RNA-seq/20200924/count/count_new.csv',sep ='\t')

df['H6C7']=(df['H6C7_1']+df['H6C7_2']+df['H6C7_3'])/3
df['BXPC']=(df['Bxpc_1']+df['Bxpc_2']+df['Bxpc_3'])/3
df['PANC']=(df['Panc_1']+df['Panc_2']+df['Panc_3'])/3
df = df.loc[:,['gene_id','H6C7','BXPC','PANC']]

df.to_csv('/public/home/shidetong/projects/yf/RNA-seq/20200924/count/gecount_nes.csv',sep ='\t',index = None)