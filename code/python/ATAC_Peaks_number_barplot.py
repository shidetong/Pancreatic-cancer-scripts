

from __future__ import division
import numpy as np
import csv
from itertools import islice
from sklearn.cluster import KMeans
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster import hierarchy
from scipy import cluster 
import seaborn as sns
import copy
import scipy
import scipy.cluster.hierarchy as sch
from itertools import islice
import os

cells = ['CCS' , 'NT5' , 'NT6' ]

p_type = ({'names':['chr' , 'start' , 'end'],
             'formats':['U8' , np.int , np.int]})




def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
    
    
    
datafolder = '/public/home/shidetong/projects/yf/ATAC/mac'
CCS = np.loadtxt(os.path.join(datafolder , 'H6C7_ATAC_2_bed_peaks.narrowPeak') , dtype = p_type , usecols=(0 , 1 , 2))
NT5 = np.loadtxt(os.path.join(datafolder , 'Bxpc_ATAC_1_bed_peaks.narrowPeak') , dtype = p_type , usecols=(0 , 1 , 2))
NT6 = np.loadtxt(os.path.join(datafolder , 'Panc_ATAC_1_bed_peaks.narrowPeak') , dtype = p_type , usecols=(0 , 1 , 2))


p_numbers = {'CCS':len(CCS) , 'NT5':len(NT5) , 'NT6':len(NT6) }



Left = 0.2 ; HB = 0.2 ; width = 0.6 ; HH = 0.6


fig = plt.figure(figsize = (12 , 12))
ax = fig.add_axes([Left  , HB , width , HH])
ax.bar([0 , 1, 2 ] , p_numbers.values() , color = 'c')

ax.set_xlim((-0.5 , 4.5))
# ax.set_ylim((0 , 40000))
ax.set_title('ATAC Peak Numbers')
ax.set_xticks([0 , 1 , 2 ])
ax.set_xticklabels(['CCS' , 'NT5' , 'NT6' ])
ax.text(-0.1 , list(p_numbers.values())[0] + 20 , str(list(p_numbers.values())[0]))
ax.text(0.9 , list(p_numbers.values())[1] + 20 , str(list(p_numbers.values())[1]))
ax.text(1.9 , list(p_numbers.values())[2] + 20 , str(list(p_numbers.values())[2]))


run_Plot(fig , 'F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\ATAC\\ATAC_new\\peaks\\idr_new\\ATAC_Peak_Numbers_barplot_new.pdf')






















