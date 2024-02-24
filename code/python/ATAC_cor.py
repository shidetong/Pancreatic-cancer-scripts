import pandas as pd 
import seaborn as sns 
import matplotlib.pyplot as plt 

cor = pd.read_csv('/public/home/shidetong/projects/yf/ATAC/bw/merge_plot/pearsonCorr_readCounts.tab',sep = '\t',skiprows = 2,header = None)
cor = cor.iloc[:,[1,2,3]]


ax = sns.heatmap(cor,vmax=0.88,vmin = 0.81,annot=False,cmap='GnBu',square=True)
ax.set_title("ATAC pearson correlation",size = 20)
plt.savefig('/public/home/shidetong/projects/yf/ATAC/bw/merge_plot/ATAC_cor_new.pdf')
