

from __future__ import division
import numpy as np
#from tadlib.calfea.analyze import getmatrix
import matplotlib
# Use a non-interactive backend
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import os
import pyBigWig
import seaborn as sns
from scipy.interpolate import  interp1d
from scipy import stats
import pandas as pd
from itertools import islice

#--------------------------------------------------------------------------
## Matplotlib Settings
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
# Our Own Color Map
my_cmap = LinearSegmentedColormap.from_list('interaction',
                                            ['#FFFFFF','#CD0000'])
my_cmap.set_bad('#2672a1')




def Load_mm10(mm10_fil):
    '''
    mm10_fil: gene length file  , formats: ['chr' , 'length']
    
    '''
    mm10_type = ({'names':['chr' , 'length'],
                  'formats':['U8' , np.int]})
    mm10 = {}
    data = np.loadtxt(mm10_fil , dtype = mm10_type)
    for i in data:
        mm10[i['chr']] = i['length']
        
    return mm10



def Sig_To_100bp(signal_fil , chro):
    """
    signal_fil: bigwig file
    chro: 'chr' or ''
    
    """
    
    signal = pyBigWig.open(signal_fil)
    
    New_Data = {}
    for g in chroms:
        New_Data[g] = {}
        tmp_data = np.array(list(signal.intervals(chro + g)) , dtype = signal_type)
        max_ = tmp_data['end'].max()
        bin_size = max_ // 100 + 1
        New_Data[g] = np.zeros((bin_size,))
        for line in tmp_data:
            start = line['start'] // 100
            # end = line['end'] // 100
            # for i in range(start,end):
            New_Data[g][start] += line['value']
    
    return New_Data




def caxis_S(ax):
    """
    Axis Control for signal plots.
    """
    for spine in ['right', 'top']:
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis = 'x', bottom = False, top = False, left = False,
                   right = False, labelbottom = False, labeltop = False,
                   labelleft = False, labelright = False , labelsize = 12)
    ax.tick_params(axis = 'y', bottom = False, top = False, left = True,
                   right = False, labelbottom = False, labeltop = False,
                   labelleft = True, labelright = False , labelsize = 12)
    ax.spines['left'].set_lw(1.5)
    ax.spines['left'].set_color('black')
    ax.spines['left'].set_alpha(0.9)
    ax.spines['left'].set_linestyle('dotted')
    



def properU(pos):
    """
    Express a genomic position in a proper unit (KB, MB, or both).
    
    """
    i_part = int(pos) // 1000000 # Integer Part
    d_part = (int(pos) % 1000000) // 1000 # Decimal Part
    
    if (i_part > 0) and (d_part > 0):
        return ''.join([str(i_part), 'M', str(d_part), 'K'])
    elif (i_part == 0):
        return ''.join([str(d_part), 'K'])
    else:
        return ''.join([str(i_part), 'M'])
    
    

def standard_axes_lim(N,*figs):
    """
    N: ax numbers
    """    
    lim_min = []
    lim_max = []
    for s_f in figs:
        lim_min.append(s_f.axes[N].get_ylim()[0])
        lim_max.append(s_f.axes[N].get_ylim()[1])
    
    for s_f in figs:
        s_f.axes[N].set_ylim(min(lim_min),max(lim_max))
        s_f.axes[N].set_yticks([round(min(lim_min)/2 , 1) , 0 , round(max(lim_max)/2 , 1)])
        


def Sig_Plot(data,start,end,chro,fig,location,color,label):
    """
    data: signal data 
    
    """
    tmp = data[chro]
    # sig_start = start // R - 50
    # sig_end = end // R + 50
    sig_data = tmp[start:end]
    max_ = sig_data.max()
    

    ax = fig.add_axes(location)
    ax.fill_between(np.arange(len(sig_data)),sig_data, facecolor = color, edgecolor = 'none')

    ax.set_ylabel(label,fontsize = 15,rotation = 'horizontal',labelpad = 50)
    ax.set_xlim((0,len(sig_data)))
    # ax.set_xlabel(chro,fontsize = 5,rotation = 'horizontal',labelpad = 20)
    
    caxis_S(ax)
    return max_
    
    

def Load_gtf(gtfil):
    gtf_type = np.dtype({'names':['gene_id' , 'gene_name' , 'chr' , 'strand' , 'start' , 'end'],
                     'formats':['S64' , 'S64' , 'S8' , 'S4' , np.int , np.int]})
    gtf = open(gtfil , 'r')
    gtf_1 = []
    for i in islice(gtf , 5 , None):
        a = i.strip().split()
        if a[2] == 'gene':
            gene_id = i.strip().split('\"')[1]
            gene_name = i.strip().split('\"')[5]
            chro = a[0]
            strand = a[6]
            start = a[3]
            end = a[4]
            gtf_1.append((gene_id , gene_name , chro , strand , start , end))
    gtf = np.array(gtf_1 , dtype = gtf_type)
    return gtf


    

##-----------------------data_type--------------------------------

signal_type = np.dtype({'names':['start' , 'end' , 'value'] , 
                    'formats':[np.int , np.int , np.float]})


datatype = ({'names':['chr' , 'start' , 'end' , 'score'],
             'formats':['S8' , np.int , np.int , np.float]})

p_type = ({'names':['chr' , 'start' , 'end'],
             'formats':['S8' , np.int , np.int]})

g_type = ({'names':['gene_name' , 'chr' , 'start' , 'end'],
           'formats':['S64' , 'S8' , np.int , np.int]})


R= 1000
chroms = ['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19' , '20' , '21' , '22' , 'X']
cells = ['Bxpc' , 'H6C7' , 'Panc']

##----------------------files----------------------------------------

RNA_signalfolder = '/public/home/shidetong/projects/yf/RNA-seq/20200924/bw'



#Chip_Data


Chip_Bxpc_CTCF = Sig_To_100bp('/public/home/shidetong/projects/yf/chip-seq/chip_2/bw/Bxpc_CTCF-1.bw' , 'chr')
Chip_H6C7_CTCF = Sig_To_100bp('/public/home/shidetong/projects/yf/chip-seq/chip_2/bw/H6C7_CTCF-1.bw' , 'chr')
Chip_Panc_CTCF = Sig_To_100bp('/public/home/shidetong/projects/yf/chip-seq/chip_2/bw/Panc_CTCF-1.bw' , 'chr')


CTCF = {'Bxpc':Chip_Bxpc_CTCF , 
        'H6C7':Chip_H6C7_CTCF , 
        'Panc':Chip_Panc_CTCF}



Chip_Bxpc_H3K4 = Sig_To_100bp('/public/home/shidetong/projects/yf/chip-seq/20201223/Bxpc/bw/1223BxpcH3K4.bw' , 'chr')
Chip_H6C7_H3K4 = Sig_To_100bp('/public/home/shidetong/projects/yf/chip-seq/20210107/H6C7/bw/0107H6C7H3K4.bw' , 'chr')
Chip_Panc_H3K4 = Sig_To_100bp('/public/home/shidetong/projects/yf/chip-seq/20201223/Panc/bw/1223PANCH3K4.bw' , 'chr')


H3K4 = {'Bxpc':Chip_Bxpc_H3K4 , 
        'H6C7':Chip_H6C7_H3K4 , 
        'Panc':Chip_Panc_H3K4}



Chip_Bxpc_Rad21 = Sig_To_100bp('/public/home/shidetong/projects/yf/chip-seq/chip_2/bw/Bxpc_Rad21.bw' , 'chr')
Chip_H6C7_Rad21 = Sig_To_100bp('/public/home/shidetong/projects/yf/chip-seq/chip_2/bw/H6C7_Rad21.bw' , 'chr')
Chip_Panc_Rad21 = Sig_To_100bp('/public/home/shidetong/projects/yf/chip-seq/chip_2/bw/Panc_Rad21.bw' , 'chr')


Rad21 = {'Bxpc':Chip_Bxpc_Rad21 , 
         'H6C7':Chip_H6C7_Rad21 , 
         'Panc':Chip_Panc_Rad21}



Bxpc_ATAC = Sig_To_100bp('/public/home/shidetong/projects/yf/ATAC/bw/Bxpc_ATAC_1.bw' , 'chr')
H6C7_ATAC = Sig_To_100bp('/public/home/shidetong/projects/yf/ATAC/bw/H6C7_ATAC_3.bw' , 'chr')
Panc_ATAC = Sig_To_100bp('/public/home/shidetong/projects/yf/ATAC/bw/Panc_ATAC_1.bw' , 'chr')


ATAC = {'Bxpc':Bxpc_ATAC , 
         'H6C7':H6C7_ATAC , 
         'Panc':Panc_ATAC}




#RNA_Data


RNA_Bxpc = Sig_To_100bp(os.path.join(RNA_signalfolder , 'Bxpc.bw') , 'chr')
RNA_H6C7 = Sig_To_100bp(os.path.join(RNA_signalfolder , 'H6C7.bw') , 'chr')
RNA_Panc = Sig_To_100bp(os.path.join(RNA_signalfolder , 'Panc.bw') , 'chr')


RNA_Data = {'Bxpc':RNA_Bxpc , 
            'H6C7':RNA_H6C7 , 
            'Panc':RNA_Panc}








#----------------------------Gtf_files-----------------------------

hg19 = '/public/home/shidetong/wedata/test/genenum/gencode.v37lift37.annotation.gtf'
gtf = Load_gtf(hg19)



##-----------------------------Gene_list-------------------------


# gene_list = pd.read_csv('/public/home/shidetong/projects/yf/Plot/Signal_plot/gene_list.csv' , header=None)

# genes = [] 
# for i in range(len(gene_list)):
#     gene = gene_list.iloc[i,0].rstrip('?')
#     if gene in gtf['gene_name']:
#         genes.append(gene)

# Genes = {'1_1000':genes[:1000] , '1001_2000':genes[1000:2000] , '2001_3078':genes[2000:3000]}



# gene_list = pd.read_csv('/public/home/shidetong/projects/yf/RNA-seq/20200924/dif/overlap/list/H_low.csv' , header=None)

# genes = [] 
# for i in range(len(gene_list)):
#     gene = gene_list.iloc[i,0].rstrip('?')
#     if gene in gtf['gene_name']:
#         genes.append(gene)

# Genes = {'1_93':genes[:93]}



gene_list = pd.read_csv('/public/home/shidetong/projects/yf/Plot/new_sig/gene_15.csv' , header=None)

genes = [] 
for i in range(len(gene_list)):
    gene = gene_list.iloc[i,0].rstrip('?')
    if gene in gtf['gene_name']:
        genes.append(gene)

Genes = {'1_15':genes[:15]}
# ge = ['CCSER1','Mtap','mir31',
#       'Pde4d','RIPK4',
#       'Smad4','Tmprss2','tnc']

# ge = ['Smad4']

# Genes = {'1_1':ge[:1]}
##   'FAM190A','IFNA1','IFNA2','IFNA5','IFNA13','IFNA17'
# gene_list = pd.read_csv('/public/home/shidetong/projects/yf/RNA-seq/20200924/dif/overlap/list/H_hig.csv' , header=None)

# genes = [] 
# for i in range(len(gene_list)):
#     gene = gene_list.iloc[i,0].rstrip('?')
#     if gene in gtf['gene_name']:
#         genes.append(gene)

# Genes = {'1_65':genes[:65]}




##----------------------------Plot_sigs------------------------------------



size = (12, 20)
Left = 0.15 ; HB = 0.15 ; width = 0.7 ; HH = 0.05




# for cl in ['1_93']:
#     pp = PdfPages('/public/home/shidetong/projects/yf/Plot/new_sig/low_100_' + cl + '.pdf')
#     for i in Genes[cl]:
#         gene = gtf[gtf['gene_name'] == i][0]
#         g = gene['chr'].lstrip('chr')
#         start = gene['start'] // 100 - 10  
#         end = gene['end'] // 100 + 10
#         if start < 0:
#             start = 0
#             ticks = [0 , gene['start'] // 100 , gene['end'] // 100  , end - start]
#         else:
#             ticks = [0 , 10 , gene['end'] // 100 - start , end - start]
        
#         fig = plt.figure(figsize = size)


for cl in ['1_15']:
    pp = PdfPages('/public/home/shidetong/projects/yf/Plot/new_sig/newgene_' + cl + '.pdf')
    for i in Genes[cl]:
        gene = gtf[gtf['gene_name'] == i][0]
        g = gene['chr'].lstrip('chr')
        start = gene['start'] // 100 - 10000  
        end = gene['end'] // 100 + 10000
        if start < 0:
            start = 0
            ticks = [0 , gene['start'] // 100 , gene['end'] // 100  , end - start]
        else:
            ticks = [0 , 10000 , gene['end'] // 100 - start , end - start]
        
        fig = plt.figure(figsize = size)
        #RNA
        n = 0
        max_rna = []
        for c in cells:      
            location1 = [Left , HB + n , width , HH]
            color = 'fuchsia'
            max_ = Sig_Plot(data = RNA_Data[c],
                            start = start,
                            end = end,
                            chro = g,
                            fig = fig,
                            location = location1,
                            color = color,
                            label = 'RNA_' + c)
            max_rna.append(max_)
            n += 0.05
        max_r = max(max_rna) * 1.1
        for N in range(3):
            fig.axes[N].set_ylim(0 , max_r)
            
        
        #ATAC
        max_ATAC = []
        for c in cells:
            location2 = [Left , HB + n , width , HH]
            color = 'blue'
            max_ = Sig_Plot(data = ATAC[c],
                            start = start,
                            end = end,
                            chro = g,
                            fig = fig,
                            location = location2,
                            color = color,
                            label = 'ATAC_' + c)
            max_ATAC.append(max_)
            n += 0.05
        max_p = max(max_ATAC)  * 1.1
        for N in range(3,6):
            fig.axes[N].set_ylim(0 , max_p)
            
        #Rad21
        max_Rad21 = []
        for c in cells:
            location2 = [Left , HB + n , width , HH]
            color = 'darkorange'
            max_ = Sig_Plot(data = Rad21[c],
                            start = start,
                            end = end,
                            chro = g,
                            fig = fig,
                            location = location2,
                            color = color,
                            label = 'Rad21_' + c)
            max_Rad21.append(max_)
            n += 0.05
        max_Rad21 = max(max_Rad21)  * 1.1
        for N in range(6,9):
            fig.axes[N].set_ylim(0 , max_Rad21)
            
        #H3K4
        max_H3K4 = []
        for c in cells:
            location2 = [Left , HB + n , width , HH]
            color = 'darkgreen'
            max_ = Sig_Plot(data = H3K4[c],
                            start = start,
                            end = end,
                            chro = g,
                            fig = fig,
                            location = location2,
                            color = color,
                            label = 'H3K4_' + c)
            max_H3K4.append(max_)
            n += 0.05
        max_h3k4 = max(max_H3K4)  * 1.1
        for N in range(9,12):
            fig.axes[N].set_ylim(0 , max_h3k4)
            

        #CTCF
        max_CTCF = []
        for c in cells:
            location2 = [Left , HB + n , width , HH]
            color = 'maroon'
            max_ = Sig_Plot(data = CTCF[c],
                            start = start,
                            end = end,
                            chro = g,
                            fig = fig,
                            location = location2,
                            color = color,
                            label = 'CTCF_' + c)
            max_CTCF.append(max_)
            n += 0.05
        max_ctcf = max(max_CTCF)  * 1.1
        for N in range(12,15):
            fig.axes[N].set_ylim(0 , max_ctcf)
            
            
        ax1 = fig.add_axes([Left  , HB , width , 0])
        length = end  - start 
        
        pos = [((start + t) * 100) for t in ticks]
        labels = ['' ,  i + '_' + properU(pos[1]) , properU(pos[2]) , '']
        ax1.set_xticks(ticks)
        ax1.set_xticklabels(labels)
        ax1.set_xlabel('chr' + g ,fontsize = 20 ,rotation = 'horizontal')
        pp.savefig(fig)
        plt.close()
    pp.close()
    
    

            

    
    













