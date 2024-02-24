import math
from pickle import TRUE
from string import whitespace
import numpy as np
import pandas as pd 
import csv , copy
import xlrd
from itertools import islice
from sklearn.cluster import KMeans
from scipy import stats
from matplotlib.backends.backend_pdf import PdfPages
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster import hierarchy
from scipy import cluster 
import seaborn as sns
import copy 
import scipy
import scipy.cluster.hierarchy as sch
from itertools import islice  
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from statannotations.Annotator import Annotator


my_cmap = plt.get_cmap('bwr')
my_cmap.set_bad('#2672a1')


# A2B = pd.read_csv('/public/home/shidetong/projects/yf/hic/HiCHap_workspace/Matrix/pc_gene/gene_name/A2B.csv',sep = '\t',header = None)
# B2A = pd.read_csv('/public/home/shidetong/projects/yf/hic/HiCHap_workspace/Matrix/pc_gene/gene_name/B2A.csv',sep = '\t',header = None)
 
A2B = pd.read_csv('/public/home/shidetong/projects/yf/hic/HiCHap_workspace/Matrix/pc_gene/gene_name/A2B.csv',sep = '\t',header = None)
B2A = pd.read_csv('/public/home/shidetong/projects/yf/hic/HiCHap_workspace/Matrix/pc_gene/gene_name/B2A.csv',sep = '\t',header = None)
 

Hhig = pd.read_csv('/public/home/shidetong/projects/yf/RNA-seq/20200924/dif/overlap/H_hig.csv',sep = '\t',header = None)
Hlow = pd.read_csv('/public/home/shidetong/projects/yf/RNA-seq/20200924/dif/overlap/H_low.csv',sep = '\t',header = None)




# B_AH = pd.merge(B2A,Hhig,on=[0],how = 'inner')
# A_BL = pd.merge(A2B,Hlow,on=[0],how = 'inner')


B_AH = pd.merge(B2A,Hlow,on=[0],how = 'inner')
A_BL = pd.merge(A2B,Hhig,on=[0],how = 'inner')


gene = pd.read_csv('/public/home/shidetong/projects/yf/RNA-seq/20200924/count/gene_count.csv',sep = '\t')

fpkm = pd.read_csv('/public/home/shidetong/projects/yf/RNA-seq/20200924/FPKM/fpkm_all_sig.csv',sep = '\t')
gene = fpkm

B_AH.columns = ['symbol']
gene1 = pd.merge(gene,B_AH,on=['symbol'],how='inner')
gene1 = gene1.loc[:,['H6C7']]


B_AH.columns = ['symbol']
gene2 = pd.merge(gene,B_AH,on=['symbol'],how='inner')
gene2 = gene2.loc[:,['BXPC']]

A_BL.columns = ['symbol']
gene3 = pd.merge(gene,A_BL,on=['symbol'],how='inner')
gene3 = gene3.loc[:,['H6C7']]

A_BL.columns = ['symbol']
gene4 = pd.merge(gene,A_BL,on=['symbol'],how='inner')
gene4 = gene4.loc[:,['BXPC']]


B_AH.columns = ['symbol']
gene5 = pd.merge(gene,B_AH,on=['symbol'],how='inner')
gene5 = gene5.loc[:,['PANC']]

A_BL.columns = ['symbol']
gene6 = pd.merge(gene,A_BL,on=['symbol'],how='inner')
gene6 = gene6.loc[:,['PANC']]

#构建绘制箱线图数据类型
c = [gene1['H6C7'],gene2['BXPC'],gene5['PANC'],gene3['H6C7'],gene4['BXPC'],gene6['PANC']]




fl = pd.DataFrame(c)
file = fl.T
file.columns = ['c1' , 'c2' , 'c3','c4','c5','c6']
f = file.melt()
f = f.dropna()

x= 'variable'
y= 'value'
order = ['c1' , 'c2' , 'c3','c4','c5','c6']



#绘制箱线图
fig,ax = plt.subplots(figsize=(7,5),dpi=100,facecolor="w")
my_pal = {"c1": "#FB8402", "c2": "#FB8402", "c3":"#FB8402","c4": "salmon","c5": "salmon","c6": "salmon"}
sns.boxplot(data=f, x=x, y=y, order=order,ax=ax,showfliers=False, palette=my_pal)

pairs=[("c1", "c2"), ("c5", "c4"),('c1','c3'),('c6','c4')]
annotator = Annotator(ax, pairs, data=f, x=x, y=y, order=order)
annotator.configure(test='Mann-Whitney', text_format='star',line_height=0.03,line_width=1)
annotator.apply_and_annotate()
###test:t-test_ind, t-test_welch, t-test_paired, Mann-Whitney, Mann-Whitney-gt, Mann-Whitney-ls, Levene, Wilcoxon, Kruskal, Brunner-Munzel
ax.tick_params(which='major',direction='in',length=3,width=1.,labelsize=14,bottom=False)
for spine in ["top","left","right"]:
    ax.spines[spine].set_visible(False)
# ax.spines['bottom'].set_linewidth(2)
# ax.grid(axis='y',ls='--',c='gray')
sns.despine(top=False, right=False, left=False, bottom=False)
ax.set_axisbelow(True)

plt.savefig('/public/home/shidetong/projects/yf/fig/fig2d_fpkm.png')
