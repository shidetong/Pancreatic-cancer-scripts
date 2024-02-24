import pandas as pd 
import numpy as np
import seaborn as sns 
import matplotlib.pyplot as plt 
import scipy
from statannotations.Annotator import Annotator



def TAD_size(Domain):
    pc_type = np.dtype({'names':['chr' , 'start','end'] , 
                    'formats':['U8' , np.int , np.int]})
    TAD = np.loadtxt(Domain,pc_type)
    TAD = pd.DataFrame(TAD)
    p = 0
    for i in TAD.loc[:,:].values:
        TAD.loc[p,'diff'] = i[2]-i[1]
        p += 1
    return TAD


H6C7 = TAD_size('/public/home/shidetong/projects/yf/hic/HiCHap_workspace/Matrix/H6C7/Cooler/Traditional_TADs/Traditional_TADs_Domain_40K.txt')
BXPC = TAD_size('/public/home/shidetong/projects/yf/hic/HiCHap_workspace/Matrix/Bxpc/Cooler/Traditional_TADs/Traditional_TADs_Domain_40K.txt')
PANC = TAD_size('/public/home/shidetong/projects/yf/hic/HiCHap_workspace/Matrix/PANC/Cooler/Traditional_TADs/Traditional_TADs_Domain_40K.txt')

data = [H6C7['diff'],BXPC['diff'],PANC['diff']]

fl = pd.DataFrame(data)
file = fl.T
file.columns = ['H6C7','BxPC3','PANC-1']
f = file.melt()
f = f.dropna()
f.columns = ['cells','TAD_size']

x= 'cells'
y= 'TAD_size'
order = ['H6C7','BxPC3','PANC-1']

fig,ax = plt.subplots(figsize=(12,7),dpi=100,facecolor="w")
my_pal = {"H6C7": "#FB8402", "BxPC3": "#219EBC", "PANC-1":"#90C9E6"}
sns.boxplot(data=f, x=x, y=y, order=order,ax=ax,showfliers=False, palette=my_pal)

pairs=[("H6C7", "BxPC3"), ("H6C7", "PANC-1"),('BxPC3','PANC-1')]
annotator = Annotator(ax, pairs, data=f, x=x, y=y, order=order)
annotator.configure(test='Mann-Whitney', text_format='star',line_height=0.03,line_width=1)
annotator.apply_and_annotate()
###test:t-test_ind, t-test_welch, t-test_paired, Mann-Whitney, Mann-Whitney-gt, Mann-Whitney-ls, Levene, Wilcoxon, Kruskal, Brunner-Munzel
ax.tick_params(which='major',direction='in',length=3,width=1.,labelsize=30,bottom=False)
for spine in ["top","left","right"]:
    ax.spines[spine].set_visible(False)
# ax.spines['bottom'].set_linewidth(2)
# ax.grid(axis='y',ls='--',c='gray')
ax.xaxis.label.set_size(30)
ax.yaxis.label.set_size(30)
sns.despine(top=False, right=False, left=False, bottom=False)
ax.set_axisbelow(True)

plt.savefig('/public/home/shidetong/projects/yf/hic/plot/TAD_size_newcolor.png')

