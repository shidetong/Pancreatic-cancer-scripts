from __future__ import division
import math
from string import whitespace
import numpy as np
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

# Our Own Color Map
my_cmap = plt.get_cmap('bwr')
my_cmap.set_bad('#2672a1')

#my_cmap = LinearSegmentedColormap.from_list('interaction',
#                                            ['skyblue' , 'k' , 'yellow'])
#my_cmap.set_bad('#2672a1')

pc_type = np.dtype({'names':['chr' , 'pc'] , 
                    'formats':['U8' , np.float]})
data_type1 = np.dtype({'names':['CCs' , 'NT5' , 'NT6' , 'F35' , 'F40'] , 
                    'formats':[np.float , np.float , np.float , np.float , np.float]})
data_type2 = np.dtype({'names':['CCs' ,  'F35' ] , 
                    'formats':[np.float ,  np.float ]})
data_type3 = np.dtype({'names':['chr','CCs' ,  'F35' ] , 
                    'formats':['U8',np.float ,  np.float ]})

chroms = ['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19']
chroms = ['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19','20','21','22','X']



def covert_pvalue(pvalue):
    if pvalue <= 0.0001:
        return '****'
    elif pvalue <= 0.001:
        return '***'
    elif pvalue <= 0.01:
        return '**'
    elif pvalue <= 0.05:
        return '*'
    return 'ns'

def Load_PC(PC_fil):
    '''
    '''
    data = np.loadtxt(PC_fil , dtype = pc_type)
    data_new = []
    for g in chroms:
        tmp = data[data['chr'] == g]
        for i in tmp:
            data_new.append(i)
    data_new = np.array(data_new , dtype = data.dtype)
    
    return data_new
    
    
    
# def PC_changed(CCS , NT5 , NT6 , F35 , F40):
#     A_B = []
#     B_A = []
#     for i in range(len(CCS)):
#         if CCS[i]['chr'] in chroms:
#             if (CCS[i]['pc'] > 0) and (F35[i]['pc'] < 0) and (F40[i]['pc'] < 0):
#                 A_B.append((CCS[i]['pc'] , NT5[i]['pc'] , NT6[i]['pc'] , F35[i]['pc'] , F40[i]['pc']))
#             elif (CCS[i]['pc'] < 0) and (F35[i]['pc'] > 0) and (F40[i]['pc'] > 0):
#                 B_A.append((CCS[i]['pc'] , NT5[i]['pc'] , NT6[i]['pc'] , F35[i]['pc'] , F40[i]['pc']))
#     A_B = np.array(A_B , dtype = data_type)
#     B_A = np.array(B_A , dtype = data_type)
#     return A_B , B_A

def PC_changed(CCS,F35):
    A_B = []
    B_A = []
    A_A = []
    B_B = []
    for i in range(len(CCS)):
        if CCS[i]['chr'] in chroms:
            if (CCS[i]['pc'] > 0) and (F35[i]['pc'] < 0) :
                A_B.append((CCS[i]['pc'] ,  F35[i]['pc'] ))
            elif (CCS[i]['pc'] < 0) and (F35[i]['pc'] > 0) :
                B_A.append((CCS[i]['pc'] ,  F35[i]['pc'] ))
            elif (CCS[i]['pc'] > 0) and (F35[i]['pc'] > 0):
                A_A.append((CCS[i]['pc'] ,  F35[i]['pc'] ))
            elif (CCS[i]['pc'] < 0) and (F35[i]['pc'] < 0) :
                A_A.append((CCS[i]['pc'] ,  F35[i]['pc']))
    A_B = np.array(A_B , dtype = data_type2)
    B_A = np.array(B_A , dtype = data_type2)
    A_A = np.array(A_A , dtype = data_type2)
    B_B = np.array(B_B , dtype = data_type2)
    print(len(A_B),len(B_A),len(A_A),len(B_B))
    return A_B , B_A , A_A , B_B

def PC_changed(CCS,F35):
    ##查看ABpc所在染色体号
    A_B = []
    B_A = []
    A_A = []
    B_B = []
    for i in range(len(CCS)):
        if CCS[i]['chr'] in chroms:
            if (CCS[i]['pc'] > 0) and (F35[i]['pc'] < 0) :
                A_B.append((CCS[i]['chr'],CCS[i]['pc'] ,  F35[i]['pc']))
            elif (CCS[i]['pc'] < 0) and (F35[i]['pc'] > 0) :
                B_A.append((CCS[i]['pc'] ,  F35[i]['pc'] ))
            elif (CCS[i]['pc'] > 0) and (F35[i]['pc'] > 0):
                A_A.append((CCS[i]['pc'] ,  F35[i]['pc'] ))
            elif (CCS[i]['pc'] < 0) and (F35[i]['pc'] < 0) :
                B_B.append((CCS[i]['pc'] ,  F35[i]['pc']))
    A_B = np.array(A_B , dtype = data_type3)
    B_A = np.array(B_A , dtype = data_type2)
    A_A = np.array(A_A , dtype = data_type2)
    B_B = np.array(B_B , dtype = data_type2)
    print(len(A_B),len(B_A),len(A_A),len(B_B))
    return A_B , B_A , A_A , B_B


def Box_plot(data):                
    left, bottom, width, height = 0.2 , 0.2 , 0.6 , 0.7
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (12, 12),facecolor='white',)
    ax = fig.add_axes(size_axes)
    ax.boxplot(data[0] , positions=[1] , showfliers=False, widths = 0.7 , 
            boxprops={'color': 'darkred','linewidth':2},
            medianprops={'color':'darkred','linewidth':2},
            capprops={'color':'darkred','linewidth':2},
            whiskerprops={'color':'darkred','linewidth':2})
    ax.boxplot(data[1] , positions=[2] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'dodgerblue','linewidth':2},
            medianprops={'color':'dodgerblue','linewidth':2},
            capprops={'color':'dodgerblue','linewidth':2},
            whiskerprops={'color':'dodgerblue','linewidth':2})
    ax.boxplot(data[2] , positions=[3] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'orange','linewidth':2},
            medianprops={'color':'orange','linewidth':2},
            capprops={'color':'orange','linewidth':2},
            whiskerprops={'color':'orange','linewidth':2})
    # ax.boxplot(data[3] , positions=[4] , showfliers=False, widths = 0.7 ,
    #         boxprops={'color': 'green','linewidth':2},
    #         medianprops={'color':'green','linewidth':2},
    #         capprops={'color':'green','linewidth':2},
    #         whiskerprops={'color':'green','linewidth':2})
    # ax.boxplot(data[4] , positions=[5] , showfliers=False, widths = 0.7 ,
    #         boxprops={'color': 'fuchsia','linewidth':2},
    #         medianprops={'color':'fuchsia','linewidth':2},
    #         capprops={'color':'fuchsia','linewidth':2},
    #         whiskerprops={'color':'fuchsia','linewidth':2})
    
    ax.set_facecolor('white')
    
    ax.plot([0.5,3.5],[0,0], lw = 1.5, ls = '--', color = 'darkblue')
    ax.set_xticks([1 , 2 , 3 ])
    ax.set_xticklabels(['Stable' , 'A to B' , 'B to A' ] , fontsize = 20)
    ax.set_yticklabels( fontsize = 20)
    ax.set_xlim((0.5 , 3.5))
    ax.set_ylim((-0.25 , 0.25))
    # ax.set_xlabel('x label',color = 'k')
    # ax.set_ylabel('x label',color = 'k')
    ax = plt.gca() # 获取当前的axes
    ax.spines['right'].set_color('k')
    ax.spines['left'].set_color('k')
    ax.spines['top'].set_color('k')
    ax.spines['bottom'].set_color('k')
    # 设置 ax1 区域背景颜色  
    # plt.rcParams['axes.facecolor']='white'
    # plt.rcParams['savefig.facecolor']='white'
    # ax.tick_params(axis='y',
    #              labelsize=15, # y轴字体大小设置
    #              color='k',    # y轴标签颜色设置  
    #              labelcolor='k', # y轴字体颜色设置
    #              direction='in' # y轴标签方向设置
    #               ) 
    # ax.tick_params(axis='x',
    #              labelsize=25, # y轴字体大小设置
    #              color='k',    # y轴标签颜色设置  
    #              labelcolor='k', # y轴字体颜色设置
    #              direction='in' # y轴标签方向设置
    #               ) 
    return fig
             
def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
    
                   
                   

CCS = Load_PC('/public/home/shidetong/projects/yf/hic/HiCHap_workspace/Matrix/H6C7/Cooler/Traditional_PC/Traditional_PC_Compartment_500K.txt')
F35 = Load_PC('/public/home/shidetong/projects/yf/hic/HiCHap_workspace/Matrix/Bxpc/Cooler/Traditional_PC/Traditional_PC_Compartment_500K.txt')
F35 = Load_PC('/public/home/shidetong/projects/yf/hic/HiCHap_workspace/Matrix/PANC/Cooler/Traditional_PC/Traditional_PC_Compartment_500K.txt')


A_B , B_A , A_A , B_B = PC_changed(CCS,F35)

# A_B , B_A = PC_changed(CCS , NT5 , NT6 , F35 , F40)

data1 = [A_B['CCs'] , A_B['F35'] , A_B['F40']]
data2 = [B_A['CCs'] , B_A['F35'] , B_A['F40']]
data3 = [A_A['CCs'] , A_A['F35'] , A_A['F40']]
data4 = [B_B['CCs'] , B_B['F35'] , B_B['F40']]

data5 = [A_A['CCs'],A_B['CCs'],B_A['CCs']]
data6 = [A_A['F35'],A_B['F35'],B_A['F35']]
data7 = [A_A['F40'],A_B['F40'],B_A['F40']]

fig1 = Box_plot(data1)
fig2 = Box_plot(data2)
fig3 = Box_plot(data3)
fig4 = Box_plot(data4)
fig5 = Box_plot(data5)
fig6 = Box_plot(data6)
fig7 = Box_plot(data7)
run_Plot(fig5 , '/public/home/shidetong/projects/yf/hic/plot/Panc_AB.pdf')
# run_Plot(fig2 , 'F:\\work\\ntESC_3Dreprogramming\\Figures\\Plot_new\\Fig2\\compartment_B_A_200K.pdf')



 