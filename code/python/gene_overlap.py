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

# Our Own Color Map
my_cmap = plt.get_cmap('bwr')
my_cmap.set_bad('#2672a1')

"""
ATAC gene 与 组蛋白gene   分类后的 交集
"""

c1 = pd.read_csv('/public/home/shidetong/projects/yf/ATAC/diffbind/mac/classfilypeak/c1_gene.csv',header=None)
c2 =pd.read_csv('/public/home/shidetong/projects/yf/ATAC/diffbind/mac/classfilypeak/c2_gene.csv',header = None)
c3 = pd.read_csv('/public/home/shidetong/projects/yf/ATAC/diffbind/mac/classfilypeak/c3_gene.csv',header = None)


h6c7h3k4me3 = pd.read_csv('/public/home/shidetong/projects/yf/chip-seq/histon_ATAC_dif/histongene/H6C7H3K4me3_gene.csv',header = None)
h6c7h3k27me3 = pd.read_csv('/public/home/shidetong/projects/yf/chip-seq/histon_ATAC_dif/histongene/H6C7H3K27me3_gene.csv',header = None)
bxpch3k4me3 = pd.read_csv('/public/home/shidetong/projects/yf/chip-seq/histon_ATAC_dif/histongene/BXPCH3K4me3_gene.csv',header = None)
bxpch3k27me3 = pd.read_csv('/public/home/shidetong/projects/yf/chip-seq/histon_ATAC_dif/histongene/BXPCH3K27me3_gene.csv',header = None)
panch3k4me3 = pd.read_csv('/public/home/shidetong/projects/yf/chip-seq/histon_ATAC_dif/histongene/PANCH3K4me3_gene.csv',header = None)
panch3k27me3 = pd.read_csv('/public/home/shidetong/projects/yf/chip-seq/histon_ATAC_dif/histongene/PANCH3K27me3_gene.csv',header = None)


gene = pd.read_csv('/public/home/shidetong/projects/yf/RNA-seq/20200924/count/gene_count.csv',sep = '\t')

fpkm = pd.read_csv('/public/home/shidetong/projects/yf/RNA-seq/20200924/FPKM/fpkm_all_sig.csv',sep = '\t')
gene = fpkm
###H3K4
c1_h6c7 = pd.merge(c1,h6c7h3k4me3,on=[0],how = 'inner')
c2_h6c7 = pd.merge(c2,h6c7h3k4me3,on = [0],how = 'inner')
c3_h6c7 = pd.merge(c3,h6c7h3k4me3,on = [0],how = 'inner')

c1_bxpc = pd.merge(c1,bxpch3k4me3,on = [0],how='inner')
c2_bxpc = pd.merge(c2,bxpch3k4me3,on = [0],how='inner')
c3_bxpc = pd.merge(c3,bxpch3k4me3,on = [0],how='inner')

c1_panc = pd.merge(c1,panch3k4me3,on = [0],how='inner')
c2_panc = pd.merge(c2,panch3k4me3,on = [0],how='inner')
c3_panc = pd.merge(c3,panch3k4me3,on = [0],how='inner')



def merge_gene(g,ge):
    g.columns = ['symbol']
    gene = pd.merge(ge,g,on=['symbol'],how = 'inner')
    gene = gene.loc[:['H6C7']]
    return gene

c1_h6c7.columns = ['symbol']
gene1 = pd.merge(gene,c1_h6c7,on=['symbol'],how='inner')
gene1 = gene1.loc[:,['H6C7']]
# gene1.to_csv('/public/home/shidetong/projects/yf/chip-seq/histon_ATAC_dif/gene_overlap/c1_h6c7_h3k4me3.csv',index = None,sep ='\t')


c2_h6c7.columns = ['symbol']
gene2 = pd.merge(gene,c2_h6c7,on=['symbol'],how='inner')
gene2 = gene2.loc[:,['H6C7']]
# gene2.to_csv('/public/home/shidetong/projects/yf/chip-seq/histon_ATAC_dif/gene_overlap/c2_h6c7_h3k4me3.csv',index = None,sep ='\t')


c3_h6c7.columns = ['symbol']
gene3 = pd.merge(gene,c3_h6c7,on=['symbol'],how='inner')
gene3 = gene3.loc[:,['H6C7']]
# gene3.to_csv('/public/home/shidetong/projects/yf/chip-seq/histon_ATAC_dif/gene_overlap/c3_h6c7_h3k4me3.csv',index = None,sep ='\t')


c1_bxpc.columns = ['symbol']
gene1 = pd.merge(gene,c1_bxpc,on=['symbol'],how='inner')
gene4 = gene1.loc[:,['BXPC']]
# gene1.to_csv('/public/home/shidetong/projects/yf/chip-seq/histon_ATAC_dif/gene_overlap/c1_bxpc_h3k4me3.csv',index = None,sep ='\t')


c2_bxpc.columns = ['symbol']
gene2 = pd.merge(gene,c2_bxpc,on=['symbol'],how='inner')
gene5 = gene2.loc[:,['BXPC']]
# gene2.to_csv('/public/home/shidetong/projects/yf/chip-seq/histon_ATAC_dif/gene_overlap/c2_bxpc_h3k4me3.csv',index = None,sep ='\t')


c3_bxpc.columns = ['symbol']
gene3 = pd.merge(gene,c3_bxpc,on=['symbol'],how='inner')
gene6 = gene3.loc[:,['BXPC']]
# gene3.to_csv('/public/home/shidetong/projects/yf/chip-seq/histon_ATAC_dif/gene_overlap/c3_bxpc_h3k4me3.csv',index = None,sep ='\t')


c1_panc.columns = ['symbol']
gene1 = pd.merge(gene,c1_panc,on=['symbol'],how='inner')
gene7 = gene1.loc[:,['PANC']]
# gene1.to_csv('/public/home/shidetong/projects/yf/chip-seq/histon_ATAC_dif/gene_overlap/c1_panc_h3k4me3.csv',index = None,sep ='\t')

c2_panc.columns = ['symbol']
gene2 = pd.merge(gene,c2_panc,on=['symbol'],how='inner')
gene8 = gene2.loc[:,['PANC']]
# gene2.to_csv('/public/home/shidetong/projects/yf/chip-seq/histon_ATAC_dif/gene_overlap/c2_panc_h3k4me3.csv',index = None,sep ='\t')


c3_panc.columns = ['symbol']
gene3 = pd.merge(gene,c3_panc,on=['symbol'],how='inner')
gene9 = gene3.loc[:,['PANC']]
# gene3.to_csv('/public/home/shidetong/projects/yf/chip-seq/histon_ATAC_dif/gene_overlap/c3_panc_h3k4me3.csv',index = None,sep ='\t')



###plot

def Box_plot(data):                
    left, bottom, width, height = 0.2 , 0.2 , 0.6 , 0.7
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (12, 12),facecolor='white',)
    ax = fig.add_axes(size_axes)
    ax.boxplot(data[0] , positions=[1] , showfliers=False, widths = 0.4 , 
            boxprops={'color': 'darkred','linewidth':2},
            medianprops={'color':'darkred','linewidth':2},
            capprops={'color':'darkred','linewidth':2},
            whiskerprops={'color':'darkred','linewidth':2})
    ax.boxplot(data[1] , positions=[2] , showfliers=False, widths = 0.4 ,
            boxprops={'color': 'darkred','linewidth':2},
            medianprops={'color':'darkred','linewidth':2},
            capprops={'color':'darkred','linewidth':2},
            whiskerprops={'color':'darkred','linewidth':2})
    ax.boxplot(data[2] , positions=[3] , showfliers=False, widths = 0.4 ,
            boxprops={'color': 'darkred','linewidth':2},
            medianprops={'color':'darkred','linewidth':2},
            capprops={'color':'darkred','linewidth':2},
            whiskerprops={'color':'darkred','linewidth':2})
    ax.boxplot(data[3] , positions=[4] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'dodgerblue','linewidth':2},
            medianprops={'color':'dodgerblue','linewidth':2},
            capprops={'color':'dodgerblue','linewidth':2},
            whiskerprops={'color':'dodgerblue','linewidth':2})
    ax.boxplot(data[4] , positions=[5] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'dodgerblue','linewidth':2},
            medianprops={'color':'dodgerblue','linewidth':2},
            capprops={'color':'dodgerblue','linewidth':2},
            whiskerprops={'color':'dodgerblue','linewidth':2})
    ax.boxplot(data[5] , positions=[6] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'dodgerblue','linewidth':2},
            medianprops={'color':'dodgerblue','linewidth':2},
            capprops={'color':'dodgerblue','linewidth':2},
            whiskerprops={'color':'dodgerblue','linewidth':2})
    ax.boxplot(data[6] , positions=[7] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'fuchsia','linewidth':2},
            medianprops={'color':'fuchsia','linewidth':2},
            capprops={'color':'fuchsia','linewidth':2},
            whiskerprops={'color':'fuchsia','linewidth':2})
    ax.boxplot(data[7] , positions=[8] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'fuchsia','linewidth':2},
            medianprops={'color':'fuchsia','linewidth':2},
            capprops={'color':'fuchsia','linewidth':2},
            whiskerprops={'color':'fuchsia','linewidth':2})
    ax.boxplot(data[8] , positions=[9] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'fuchsia','linewidth':2},
            medianprops={'color':'fuchsia','linewidth':2},
            capprops={'color':'fuchsia','linewidth':2},
            whiskerprops={'color':'fuchsia','linewidth':2})
    

    
    # ax.plot([0.5,3.5],[0,0], lw = 1.5, ls = '--', color = 'darkblue')
    ax.set_xticks([1 , 2 , 3,4,5,6,7,8,9])
    ax.set_xticklabels(['c1' , 'c2' , 'c3','c4','c5','c6' ,'c7','c8','c9'] , fontsize = 20)
    ax.set_xlim((0.5 ,9.5))
    ax.text((1+2)*.5,2,'T')


    return fig


data = [gene1['H6C7'],gene2['H6C7'],gene3['H6C7']]


c1 = [gene1['H6C7'],gene4['BXPC'],gene7['PANC']]
c2 = [gene2['H6C7'],gene5['BXPC'],gene8['PANC']]
c3 = [gene3['H6C7'],gene6['BXPC'],gene9['PANC']]

c = [gene1['H6C7'],gene4['BXPC'],gene7['PANC'],
    gene2['H6C7'],gene5['BXPC'],gene8['PANC'],
    gene3['H6C7'],gene6['BXPC'],gene9['PANC']]



###snsplot
import seaborn as sns 

fl = pd.DataFrame(c)
file = fl.T
file.columns = ['c1' , 'c2' , 'c3','c4','c5','c6' ,'c7','c8','c9']
f = file.melt()
f = f.dropna()

x= 'variable'
y= 'value'
order = ['c1' , 'c2' , 'c3','c4','c5','c6' ,'c7','c8','c9']

fig,ax = plt.subplots(figsize=(7,5),dpi=100,facecolor="w")
my_pal = {"c1": "DarkOrange", "c2": "DarkOrange", "c3":"DarkOrange","c4": "SteelBlue", "c5": "SteelBlue", "c6":"SteelBlue","c7": "LightSkyBlue", "c8": "LightSkyBlue", "c9":"LightSkyBlue"}
sns.boxplot(data=f, x=x, y=y, order=order,ax=ax,showfliers=False, palette=my_pal)

pairs=[("c1", "c2"), ("c1", "c3"),('c2','c3'), ("c4", "c5"),('c4','c6'),('c5','c6'),('c7','c8'),('c7','c9'),('c8','c9')]
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

plt.savefig('/public/home/shidetong/projects/yf/chip-seq/histon_ATAC_dif/plot/H3K4me3_dif_bar.png')




#####
sns.boxplot(x = f['variable'],y=f['value'],showfliers=False,palette="Blues")

stat,p_value = scipy.stats.ttest_ind(f[f["variable"]=="c1"]["value"],
                                     f[f["variable"]=="c2"]["value"],
                                     equal_var=False)
p_value



###
p = 0
for i in f.loc[:,:].values:
    if i[0]=='c1' or i[0]=='c2' or i[0]== 'c3':
        f.loc[p,'cla'] = 'H6C7'
    elif i[0]=='c4' or i[0]=='c5' or i[0]=='c6':
        f.loc[p,'cla'] = 'bxpc'
    elif i[0]=='c7' or i[0]=='c8' or i[0]=='c9':
        f.loc[p,'cla'] = 'panc'
    else:
        continue
    p+=1

sns.boxplot(x = f['cla'],y=f['value'],showfliers=False,hue=f['variable'])











#h3k27me3

c1_h6c7 = pd.merge(c1,h6c7h3k27me3,on=[0],how = 'inner')
c2_h6c7 = pd.merge(c2,h6c7h3k27me3,on = [0],how = 'inner')
c3_h6c7 = pd.merge(c3,h6c7h3k27me3,on = [0],how = 'inner')

c1_bxpc = pd.merge(c1,bxpch3k27me3,on = [0],how='inner')
c2_bxpc = pd.merge(c2,bxpch3k27me3,on = [0],how='inner')
c3_bxpc = pd.merge(c3,bxpch3k27me3,on = [0],how='inner')

c1_panc = pd.merge(c1,panch3k27me3,on = [0],how='inner')
c2_panc = pd.merge(c2,panch3k27me3,on = [0],how='inner')
c3_panc = pd.merge(c3,panch3k27me3,on = [0],how='inner')



def merge_gene(g,ge):
    g.columns = ['symbol']
    gene = pd.merge(ge,g,on=['symbol'],how = 'inner')
    gene = gene.loc[:['H6C7']]
    return gene

c1_h6c7.columns = ['symbol']
gene1 = pd.merge(gene,c1_h6c7,on=['symbol'],how='inner')
gene1 = gene1.loc[:,['H6C7']]
# gene1.to_csv('/public/home/shidetong/projects/yf/chip-seq/histon_ATAC_dif/gene_overlap/c1_h6c7_h3k4me3.csv',index = None,sep ='\t')


c2_h6c7.columns = ['symbol']
gene2 = pd.merge(gene,c2_h6c7,on=['symbol'],how='inner')
gene2 = gene2.loc[:,['H6C7']]
# gene2.to_csv('/public/home/shidetong/projects/yf/chip-seq/histon_ATAC_dif/gene_overlap/c2_h6c7_h3k4me3.csv',index = None,sep ='\t')


c3_h6c7.columns = ['symbol']
gene3 = pd.merge(gene,c3_h6c7,on=['symbol'],how='inner')
gene3 = gene3.loc[:,['H6C7']]
# gene3.to_csv('/public/home/shidetong/projects/yf/chip-seq/histon_ATAC_dif/gene_overlap/c3_h6c7_h3k4me3.csv',index = None,sep ='\t')


c1_bxpc.columns = ['symbol']
gene1 = pd.merge(gene,c1_bxpc,on=['symbol'],how='inner')
gene4 = gene1.loc[:,['BXPC']]
# gene1.to_csv('/public/home/shidetong/projects/yf/chip-seq/histon_ATAC_dif/gene_overlap/c1_bxpc_h3k4me3.csv',index = None,sep ='\t')


c2_bxpc.columns = ['symbol']
gene2 = pd.merge(gene,c2_bxpc,on=['symbol'],how='inner')
gene5 = gene2.loc[:,['BXPC']]
# gene2.to_csv('/public/home/shidetong/projects/yf/chip-seq/histon_ATAC_dif/gene_overlap/c2_bxpc_h3k4me3.csv',index = None,sep ='\t')


c3_bxpc.columns = ['symbol']
gene3 = pd.merge(gene,c3_bxpc,on=['symbol'],how='inner')
gene6 = gene3.loc[:,['BXPC']]
# gene3.to_csv('/public/home/shidetong/projects/yf/chip-seq/histon_ATAC_dif/gene_overlap/c3_bxpc_h3k4me3.csv',index = None,sep ='\t')


c1_panc.columns = ['symbol']
gene1 = pd.merge(gene,c1_panc,on=['symbol'],how='inner')
gene7 = gene1.loc[:,['PANC']]
# gene1.to_csv('/public/home/shidetong/projects/yf/chip-seq/histon_ATAC_dif/gene_overlap/c1_panc_h3k4me3.csv',index = None,sep ='\t')

c2_panc.columns = ['symbol']
gene2 = pd.merge(gene,c2_panc,on=['symbol'],how='inner')
gene8 = gene2.loc[:,['PANC']]
# gene2.to_csv('/public/home/shidetong/projects/yf/chip-seq/histon_ATAC_dif/gene_overlap/c2_panc_h3k4me3.csv',index = None,sep ='\t')


c3_panc.columns = ['symbol']
gene3 = pd.merge(gene,c3_panc,on=['symbol'],how='inner')
gene9 = gene3.loc[:,['PANC']]
# gene3.to_csv('/public/home/shidetong/projects/yf/chip-seq/histon_ATAC_dif/gene_overlap/c3_panc_h3k4me3.csv',index = None,sep ='\t')


c = [gene1['H6C7'],gene4['BXPC'],gene7['PANC'],
    gene2['H6C7'],gene5['BXPC'],gene8['PANC'],
    gene3['H6C7'],gene6['BXPC'],gene9['PANC']]



###snsplot
import seaborn as sns 

fl = pd.DataFrame(c)
file = fl.T
file.columns = ['c1' , 'c2' , 'c3','c4','c5','c6' ,'c7','c8','c9']
f = file.melt()
f = f.dropna()

x= 'variable'
y= 'value'
order = ['c1' , 'c2' , 'c3','c4','c5','c6' ,'c7','c8','c9']

fig,ax = plt.subplots(figsize=(7,5),dpi=100,facecolor="w")
my_pal = {"c1": "g", "c2": "g", "c3":"g","c4": "r", "c5": "r", "c6":"r","c7": "y", "c8": "y", "c9":"y"}
sns.boxplot(data=f, x=x, y=y, order=order,ax=ax,showfliers=False, palette=my_pal)

pairs=[("c1", "c2"), ("c1", "c3"),('c2','c3'), ("c4", "c5"),('c4','c6'),('c5','c6'),('c7','c8'),('c7','c9'),('c8','c9')]
annotator = Annotator(ax, pairs, data=f, x=x, y=y, order=order)
annotator.configure(test='Mann-Whitney', text_format='star',line_height=0.03,line_width=1)
annotator.apply_and_annotate()
###test:t-test_ind, t-test_welch, t-test_paired, Mann-Whitney, Mann-Whitney-gt, Mann-Whitney-ls, Levene, Wilcoxon, Kruskal, Brunner-Munzel
ax.tick_params(which='major',direction='in',length=3,width=1.,labelsize=14,bottom=False)
for spine in ["top","left","right"]:
    ax.spines[spine].set_visible(False)
ax.spines['bottom'].set_linewidth(2)
# ax.grid(axis='y',ls='--',c='gray')
ax.set_axisbelow(True)

plt.savefig('/public/home/shidetong/projects/yf/chip-seq/histon_ATAC_dif/plot/H3K27me3_dif_bar.png')
