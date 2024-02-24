from numpy import histogram, histogram_bin_edges
import pandas as pd 
import matplotlib.pyplot as plt 
import seaborn as sns
import numpy as np



H6C7 = pd.read_csv('/public/home/shidetong/projects/yf/hic/HiCHap_workspace/Matrix/H6C7/Cooler/Traditional_TADs/Traditional_TADs_DI_40K.txt',sep = '\t',header = None) 
BXPC = pd.read_csv('/public/home/shidetong/projects/yf/hic/HiCHap_workspace/Matrix/Bxpc/Cooler/Traditional_TADs/Traditional_TADs_DI_40K.txt',sep = '\t',header = None)
PANC = pd.read_csv('/public/home/shidetong/projects/yf/hic/HiCHap_workspace/Matrix/PANC/Cooler/Traditional_TADs/Traditional_TADs_DI_40K.txt',sep = '\t',header = None)

pc = pd.concat([H6C7,BXPC,PANC],axis=1)

pc.columns = ['a','b','c','d','e','f']

IF = pc.loc[:,['b','d','f']]
IF.columns = ['H6C7','BXPC','PANC']

cor = IF.corr()
cor = pd.read_csv('/public/home/shidetong/projects/yf/ATAC/bw/merge_plot/pearsonCorr_readCounts.tab')
# ax = sns.heatmap(cor,vmax=0.4,annot=True)
ax = sns.heatmap(cor,vmax=0.87,vmin = 0.86,annot=False,cmap='GnBu',square=True,linewidths=0.5,linecolor = 'black')
# sns.clustermap(cor,vmax=0.7,annot=True)
# sns.set_style('whitegrid', {'font.sans-serif': ['simhei','FangSong']})
ax.set_title("Hi-C pearson correlation",size = 20)

plt.savefig('/public/home/shidetong/projects/yf/hic/plot/hic_cor.png')



H6C7 = np.load('/public/home/shidetong/projects/yf/hic/HiCHap_workspace/Matrix/H6C7/npz/20200518H6C7_40K.npz')
BXPC = np.load('/public/home/shidetong/projects/yf/hic/HiCHap_workspace/Matrix/Bxpc/npz/20200318Bxpc_40K.npz')

cor=np.vstack(H6C7,BXPC)