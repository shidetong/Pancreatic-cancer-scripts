import pandas as pd 
import matplotlib.pyplot as plt 
import seaborn as sns

fl = pd.read_csv('/public/home/shidetong/projects/yf/RNA-seq/20200924/count/gecount_nes.csv',sep = '\t')
fl = fl.loc[:,['H6C7','BXPC','PANC']]

# fl = pd.read_csv('/public/home/shidetong/projects/yf/RNA-seq/20200924/count/count_new.csv',sep = '\t')
# fl = fl.iloc[:,1:10]
 
cor = fl.corr()

# ax = sns.heatmap(cor,vmax=0.79,annot=True)
# ax.set_title("RNA pearson correlation",size = 20)
# plt.savefig('/public/home/shidetong/projects/yf/RNA-seq/20200924/count/RNA_cor.png')


ax = sns.heatmap(cor,vmax=0.80,vmin = 0.75,annot=False,cmap='GnBu',square=True)
ax.set_title("RNA pearson correlation",size = 20)
plt.savefig('/public/home/shidetong/projects/yf/RNA-seq/20200924/count/RNA_cor_new.pdf')




