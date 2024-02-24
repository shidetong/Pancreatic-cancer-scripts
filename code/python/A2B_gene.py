import pandas as pd 
import numpy as np
from pygments import highlight

H6C7 = pd.read_csv('/public/home/shidetong/projects/yf/hic/HiCHap_workspace/Matrix/H6C7/Cooler/Traditional_PC/Traditional_PC_Compartment_500K.txt',sep = '\t',header = None)
BXPC = pd.read_csv('/public/home/shidetong/projects/yf/hic/HiCHap_workspace/Matrix/Bxpc/Cooler/Traditional_PC/Traditional_PC_Compartment_500K.txt',sep = '\t',header = None)
PANC = pd.read_csv('/public/home/shidetong/projects/yf/hic/HiCHap_workspace/Matrix/PANC/Cooler/Traditional_PC/Traditional_PC_Compartment_500K.txt',sep = '\t',header = None)

H6C7[2] = H6C7.index
BXPC[2] = BXPC.index
PANC[2] = PANC.index


h = H6C7.sort_values(by=[0,2],ascending=[True,True])
b = BXPC.sort_values(by=[0,2],ascending=[True,True])
p = PANC.sort_values(by=[0,2],ascending=[True,True])

h.reset_index(drop=True, inplace=True)
b.reset_index(drop=True, inplace=True)
p.reset_index(drop=True, inplace=True)


h = h.iloc[:,[0,1]]
b = b.iloc[:,[0,1]]
p = p.iloc[:,[0,1]]



###chrom_size
hg19_genome_size = {
    "chr1": 249250621,
    "chr10": 135534747,
    "chr11": 135006516,
    "chr12": 133851895,
    "chr13": 115169878,
    "chr14": 107349540,
    "chr15": 102531392,
    "chr16": 90354753,
    "chr17": 81195210,
    "chr18": 78077248,
    "chr19": 59128983,
    "chr2": 243199373,
    "chr20": 63025520,
    "chr21": 48129895,
    "chr22": 51304566,
    "chr3": 198022430,
    "chr4": 191154276,
    "chr5": 180915260,
    "chr6": 171115067,
    "chr7": 159138663,
    "chr8": 146364022,
    "chr9": 141213431,
    "chrX": 155270560
}

#------------------------------------------------------------
# 定义开始值和结束值
start = 0
end = 500000

# 计算间距并生成等间距点
interval = 5
x = np.arange(start + interval, end + interval, interval)

# 创建一个空矩阵
result_matrix = np.empty((0,2), float)

# 将区间起点和终点添加到矩阵中
for i in range(len(x)):
    segment_start = x[i-1] + 1 if i > 0 else start + 1
    segment_end = x[i]
    result_matrix = np.vstack([result_matrix, [segment_start, segment_end]])

# 将最后一行的结束值改为终点
result_matrix[-1][1] = end

print(result_matrix)
#-------------------------------------------------------------------

merged_dataframe = pd.DataFrame()
for chr,value in hg19_genome_size.items():
    # 定义开始值和结束值
    start = 0
    end = value

    # 计算间距并生成等间距点
    interval = 500000
    x = np.arange(start + interval, end + interval, interval)

    # 创建一个空矩阵
    result_matrix = np.empty((0,2), float)

    # 将区间起点和终点添加到矩阵中
    for i in range(len(x)):
        segment_start = x[i-1] + 1 if i > 0 else start + 1
        segment_end = x[i]
        result_matrix = np.vstack([result_matrix, [segment_start, segment_end]])

    # 将最后一行的结束值改为终点
    result_matrix[-1][1] = end
    a = pd.DataFrame(result_matrix)
    merged_dataframe = pd.concat([merged_dataframe,a])

merged_dataframe.reset_index(drop=True, inplace=True)

H = pd.concat([h,merged_dataframe],axis=1) 
B = pd.concat([b,merged_dataframe],axis=1)
P = pd.concat([p,merged_dataframe],axis=1)
H.columns = ['chr','pc','start','end']
H[['start']] = H[['start']].astype(int)
H[['end']] = H[['end']].astype(int)
B.columns = ['chr','pc','start','end']
B[['start']] = B[['start']].astype(int)
B[['end']] = B[['end']].astype(int)
P.columns = ['chr','pc','start','end']
P[['start']] = P[['start']].astype(int)
P[['end']] = P[['end']].astype(int)
H.to_csv('/public/home/shidetong/projects/yf/hic/HiCHap_workspace/Matrix/pcaddsit/pc_sit_new/H6C7.csv',index = None,sep = '\t',header = None)
B.to_csv('/public/home/shidetong/projects/yf/hic/HiCHap_workspace/Matrix/pcaddsit/pc_sit_new/BXPC.csv',index = None,sep = '\t',header = None)
P.to_csv('/public/home/shidetong/projects/yf/hic/HiCHap_workspace/Matrix/pcaddsit/pc_sit_new/PANC.csv',index = None,sep = '\t',header = None)


####
pc_type = np.dtype({'names':['chr' , 'pc','start','end'] , 
                    'formats':['U8' , np.float,np.int,np.int]})

data_type3 = np.dtype({'names':['chr','CCs' ,  'F35','start','end' ] , 
                    'formats':['U8',np.float ,  np.float,np.int,np.int ]})

chroms = ['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19','20','21','22','X']

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

def PC_changed(CCS,F35):
    ##查看ABpc所在染色体号
    A_B = []
    B_A = []
    A_A = []
    B_B = []
    for i in range(len(CCS)):
        if CCS[i]['chr'] in chroms:
            if (CCS[i]['pc'] > 0) and (F35[i]['pc'] < 0) :
                A_B.append((CCS[i]['chr'],CCS[i]['pc'] ,  F35[i]['pc'],CCS[i]['start'],CCS[i]['end']))
            elif (CCS[i]['pc'] < 0) and (F35[i]['pc'] > 0) :
                B_A.append((CCS[i]['chr'],CCS[i]['pc'] ,  F35[i]['pc'] ,CCS[i]['start'],CCS[i]['end']))
            elif (CCS[i]['pc'] > 0) and (F35[i]['pc'] > 0):
                A_A.append((CCS[i]['chr'],CCS[i]['pc'] ,  F35[i]['pc'] ,CCS[i]['start'],CCS[i]['end']))
            elif (CCS[i]['pc'] < 0) and (F35[i]['pc'] < 0) :
                B_B.append((CCS[i]['chr'],CCS[i]['pc'] ,  F35[i]['pc'],CCS[i]['start'],CCS[i]['end']))
    A_B = np.array(A_B , dtype = data_type3)
    B_A = np.array(B_A , dtype = data_type3)
    A_A = np.array(A_A , dtype = data_type3)
    B_B = np.array(B_B , dtype = data_type3)
    print(len(A_B),len(B_A),len(A_A),len(B_B))
    return A_B , B_A , A_A , B_B

 
CCS = Load_PC('/public/home/shidetong/projects/yf/hic/HiCHap_workspace/Matrix/pcaddsit/pc_sit_new/H6C7.csv')
F35 = Load_PC('/public/home/shidetong/projects/yf/hic/HiCHap_workspace/Matrix/pcaddsit/pc_sit_new/BXPC.csv')
F35 = Load_PC('/public/home/shidetong/projects/yf/hic/HiCHap_workspace/Matrix/pcaddsit/pc_sit_new/PANC.csv')


A_B , B_A , A_A , B_B = PC_changed(CCS,F35)

HvsP_A2B = pd.DataFrame(A_B)
HvsP_B2A = pd.DataFrame(B_A)

HvsP_A2B.to_csv('/public/home/shidetong/projects/yf/hic/HiCHap_workspace/Matrix/pcaddsit/pc_sit_new/PC_AB/HvsP_A2B.csv',sep = '\t',index = None)
HvsP_B2A.to_csv('/public/home/shidetong/projects/yf/hic/HiCHap_workspace/Matrix/pcaddsit/pc_sit_new/PC_AB/HvsP_B2A.csv',sep = '\t',index = None)



    
    





##----------------------
#RNA差异基因与A2B,B2A交集
A2B = pd.read_csv('/public/home/shidetong/projects/yf/hic/HiCHap_workspace/Matrix/pc_gene/gene_name/A2B.csv',sep = '\t',header = None)
B2A = pd.read_csv('/public/home/shidetong/projects/yf/hic/HiCHap_workspace/Matrix/pc_gene/gene_name/B2A.csv',sep = '\t',header = None)
 
Hhig = pd.read_csv('/public/home/shidetong/projects/yf/RNA-seq/20200924/dif/overlap/H_hig.csv',sep = '\t',header = None)
Hlow = pd.read_csv('/public/home/shidetong/projects/yf/RNA-seq/20200924/dif/overlap/H_low.csv',sep = '\t',header = None)

ge = pd.concat([Hhig,Hlow],axis=0)

A_B = pd.merge(A2B,Hlow,on=[0],how = 'inner')
B_A = pd.merge(B2A,Hhig,on = [0],how= 'inner')


A_B = pd.merge(A2B,ge,on=[0],how = 'inner')
B_A = pd.merge(B2A,ge,on=[0],how = 'inner')

A_B.to_csv('/public/home/shidetong/projects/yf/hic/HiCHap_workspace/Matrix/pc_gene/gene_name/overlapRNAdif/low.csv',sep = '\t',header = None,index = None)
B_A.to_csv('/public/home/shidetong/projects/yf/hic/HiCHap_workspace/Matrix/pc_gene/gene_name/overlapRNAdif/hig.csv',sep = '\t',header = None,index = None)