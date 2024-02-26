from copy import PyStringMap
import pandas as pd 
import pysam
import numpy as np
import os 
import matplotlib.pyplot as plt



def  motif_seq(motif_file,outfile):
    '''
    提取motif序列，输出表格(PPM矩阵)，绘制motif序列图，后用R绘图
    '''
    motif = pd.read_csv(motif_file,sep = '\t')
    motif = motif.drop(motif.tail(3).index)
    seq = motif[['matched_sequence']]

    counta = []
    countc = []
    countt = []
    countg = []
    for p in range(20):
        c =0
        a = 0
        t = 0
        g = 0
        for i in seq.iloc[:,:].values:  
            q = i[0][p]
            if q == 'T':
                t = t+1
            elif q == 'A':
                a = a+1
            elif q == 'C':
                c= c+1
            elif q == 'G':
                g = g+1
        counta.append(a)
        countc.append(c)
        countt.append(t)
        countg.append(g)
    A = pd.DataFrame(counta)
    C = pd.DataFrame(countc)
    T = pd.DataFrame(countt)
    G = pd.DataFrame(countg)
    frames = [A,C,T,G]
    aa = pd.concat(frames,axis = 1)
    # aa = end.applymap(lambda x: x / len(seq))
    aa.columns = ['A','C','T','G']
    aa = aa.iloc[3:15,:]
    aa.to_csv(outfile,sep = '\t',index= None)
    return aa




###########
def bam(bamfile,loop_loc,outfile):
    """
    通过bam文件和，位置文件，得到所有的位置中的序列信息
    """
    bf2 = pysam.AlignmentFile(bamfile, 'rb')
    mo = open(outfile,'w')
    
    for i in list:
        for line in bf2.fetch('chr1', 4759853,  4799853):
            #print(line.reference_name,line.pos,line.seq)
            mo.write('>'+line.reference_name+'_')
            mo.write(str(line.pos)+'\n')
            mo.write(line.seq+'\n')
        mo.close()


def get_loc(loopspp):
    """
    得到loop位置信息
    """
    loop_type = np.dtype({'names':['chr' , 'start' , 'end'],
                      'formats':['S8' , np.int , np.int]})
    loops = []
    Loop = np.loadtxt(loopspp, usecols = (0,2,7) , dtype = loop_type, skiprows = 1)
    for i in Loop:
        if i['end'] - i['start'] >= 40000:
            loops.append(i) 
        else:
            continue
    loops = np.array(loops , dtype = loop_type)
    loops = loops.tolist()
    return loops 

###########



def get_loc(loc):
    loops = []
    loop_type = np.dtype({'names':['chr' , 'start' , 'end'],
                      'formats':['S8' , np.int , np.int]})
    Loop = np.loadtxt(loc, usecols = (0,1,2) , dtype = loop_type, skiprows = 1)
    for i in Loop:
        if i['end'] - i['start'] >= 40:
            loops.append(i) 
        else:
            continue
    loops = np.array(loops , dtype = loop_type)
    loops = loops.tolist()
    return loops 

#ssr = get_loc('/public/home/shidetong/wedata/meme/bins_350b/peak/loc_2.csv')

def bam_reads(bam,loc,out):
    """
    通过bam文件和，位置文件，得到所有的位置中的序列信息
    """
    a = get_loc(loc)
    a = a.sort()
    bf = pysam.AlignmentFile(bam, 'rb')
    for i in a :
        mo = open(out + i[0]+'_'+str(i[1])+'_'+str(i[2])+'.txt','w')
        for line in bf.fetch(i[0],i[1],i[2]): 
        #         print(line.reference_name,line.pos,line.seq)
            mo.write('>'+line.reference_name+'_')
            mo.write(str(line.pos)+'\n')
            mo.write(line.seq+'\n')
        mo.close()

##的所有序列信息
def bam_reads(bam,loc,out):
    """
    通过bam文件和，位置文件，得到所有的位置中的所有序列信息
    """
    a = get_loc(loc)
    # a = a.sort()
    bf = pysam.AlignmentFile(bam, 'rb')
    mo = open(out+'B6M.txt','w')
    for i in a :
        # mo = open(out + i[0]+'_'+str(i[1])+'_'+str(i[2])+'.txt','w')
        for line in bf.fetch(i[0],i[1],i[2]): 
        #         print(line.reference_name,line.pos,line.seq)
            mo.write('>'+line.reference_name+'_')
            mo.write(str(line.pos)+'\n')
            mo.write(line.seq+'\n')
    mo.close() 

def motif(file,outfile):
    """
    得到每个motif文件的具体位置file:fimo.tsv
    """
    motif = pd.read_csv(file,sep = '\t') 
    motif = motif.drop(motif.tail(3).index)   #从尾部去除三行注释信息
    seq = motif['sequence_name'].str.split('_')

    sit = []
    chro = []
    for i in seq:
        #print (i[1])
        sit.append(int(i[1]))
        chro.append(i[0])
    pdsit = pd.DataFrame(sit)
    pdsit.columns = ['sit']
    pdchr = pd.DataFrame(chro)
    pdchr.columns = ['chr']

    motif_sit = pd.concat([motif,pdsit],axis=1)
    motif_sit[['start']] = motif_sit[['start']].astype('Int64')
    motif_sit[['stop']] = motif_sit[['stop']].astype('Int64')

    #DataFrame列表中增加列
    motif_sit['chr'] = pdchr
    motif_sit['begin'] = motif_sit['start'] + motif_sit['sit']
    motif_sit['end'] = motif_sit['stop'] + motif_sit['sit']
    motif_IF = motif_sit.loc[:,['chr','begin','end','matched_sequence','p-value']]
    motif_IF.to_csv(outfile,index = None,sep = '\t')
    return motif_IF



###motif 矩阵 PWM
def  motif_seq(motif_file,outfile):
    '''
    提取motif序列，输出表格(PPM矩阵)，绘制motif序列图，后用R绘图
    '''
    motif = pd.read_csv(motif_file,sep = '\t')
    motif = motif.drop(motif.tail(3).index)
    seq = motif[['matched_sequence']]

    counta = []
    countc = []
    countt = []
    countg = []
    for p in range(20):
        c =0
        a = 0
        t = 0
        g = 0
        for i in seq.iloc[:,:].values:  
            q = i[0][p]
            if q == 'T':
                t = t+1
            elif q == 'A':
                a = a+1
            elif q == 'C':
                c= c+1
            elif q == 'G':
                g = g+1
        counta.append(a)
        countc.append(c)
        countt.append(t)
        countg.append(g)
    A = pd.DataFrame(counta)
    C = pd.DataFrame(countc)
    T = pd.DataFrame(countt)
    G = pd.DataFrame(countg)
    frames = [A,C,T,G]
    end = pd.concat(frames,axis = 1)
    aa = end.applymap(lambda x: x / len(seq))
    aa.columns = ['A','C','T','G']
    aa = aa.iloc[3:15,:]
    aa.to_csv(outfile,sep = '\t',index= None)
    return aa

fl = []
def crossover(dir,fl):
    """
    遍历文件夹下所有的.tsv文件
    """
    for i in os.listdir(dir):
        path = os.path.join(dir,i)
        if os.path.isfile(path):
            if os.path.splitext(path)[1] == '.tsv':
                fl.append(i) 
        elif os.path.isdir(path):
            newdir = path
            crossover(newdir,fl)
    return fl


#遍历文件夹，并输出绝对路径。
import pandas as pd 
dfs = pd.DataFrame()
dir  = '/public/home/shidetong/wedata/meme/out/B6M/loc_1'
fl = open('/public/home/shidetong/B6M_loc_1.txt','w')
def fordir(dir,fl):
    for root_dir,sub_dir,files in os.walk(dir):
        for file in files:
            if file.endswith('.tsv'):
                file_name = os.path.join(root_dir,file)
                fl.write(file_name + '\n')
    fl.close()


##获取碱基序列

def get_motif_seq(dir,out):
    """
    遍历全部以.tsv为结尾的文件，得到每个位置的motif序列
    """
    #dir  = '/public/home/shidetong/wedata/meme/out/B6M/loc_2'
    #out = '/public/home/shidetong/wedata/meme/seq/loc_1'
    #fl = open('/public/home/shidetong/B6M_loc_1.txt','w')
    l = []
    fl = []
    for root_dir,sub_dir,files in os.walk(dir):
        for file in files:
            if file.endswith('.tsv'):
                file_name = os.path.join(root_dir,file)
                motif = pd.read_csv(file_name,sep = '\t')
                motif = motif.drop(motif.tail(3).index)
                lenth = len(motif)
                if lenth == 0:
                    g = 0           #判断motif文件是否为空
                    fl.append(g)
                elif lenth != 0:
                    seq = motif[['matched_sequence']]
                    counta = []
                    countc = []
                    countt = []
                    countg = []
                    for p in range(20):
                        c =0
                        a = 0
                        t = 0
                        g = 0
                        for i in seq.iloc[:,:].values:
                            q = i[0][p]
                            if q == 'T':
                                t = t+1
                            elif q == 'A':
                                a = a+1
                            elif q == 'C':
                                c= c+1
                            elif q == 'G':
                                g = g+1
                        counta.append(a)
                        countc.append(c)
                        countt.append(t)
                        countg.append(g)
                    A = pd.DataFrame(counta)
                    C = pd.DataFrame(countc)
                    T = pd.DataFrame(countt)
                    G = pd.DataFrame(countg)
                    frames = [A,C,T,G]
                    end = pd.concat(frames,axis = 1)
                    aa = end.applymap(lambda x: x / len(seq))
                    aa.columns = ['A','C','T','G']
                    seq = aa.idxmax(1)  #取某一行最大值的列名
                    seq = seq.tolist()
                    seq = ''.join(seq)  #列表的多个元素合并为一个字符
                    fl.append(seq)
    df = pd.DataFrame(fl)
    df.to_csv(out,header = None,index = None,sep = '\t')

#concat_loc_all_motif
def read_motif(B6M_loc,B6P_loc,CastM_loc,CastP_loc,out):
    """
    整合四组数据，将motif一一对应
    """
    B6M = pd.read_csv(B6M_loc, sep = '\t')
    B6P = pd.read_csv(B6P_loc, sep = '\t')
    CastM = pd.read_csv(CastM_loc, sep = '\t')
    CastP = pd.read_csv(CastP_loc, sep = '\t')
    B6 = pd.merge(B6M,B6P, on = ['sit'])
    Cast = pd.merge(CastM,CastP,on = ['sit'])
    loc = pd.merge(B6,Cast,on = ['sit'])
    loc.columns = ['sit','B6M','B6P','CastM','CastP']
    #loc = pd.concat([B6M,B6P,CastM,CastP],axis =1)
    loc_T = pd.DataFrame(loc.values.T, index = loc.columns , columns = loc.index)  #矩阵转置
    loc_T.to_csv(out,index = None , sep= '\t')
    return loc_T


def  motif_seq(motif_file,outfile):
    '''
    提取motif序列，逆序输出表格，绘制motif序列图，后用R绘图
    '''
    motif = pd.read_csv(motif_file,sep = '\t')
    motif = motif.drop(motif.tail(3).index)
    seq = motif[['matched_sequence']]

    counta = []
    countc = []
    countt = []
    countg = []
    for p in range(20):
        c =0
        a = 0
        t = 0
        g = 0
        for i in seq.iloc[:,:].values:  
            q = i[0]  ##去掉中括号
            q = q[::-1]  #逆序读取
            q = q[p]
            if q == 'T':
                t = t+1
            elif q == 'A':
                a = a+1
            elif q == 'C':
                c= c+1
            elif q == 'G':
                g = g+1
        counta.append(a)
        countc.append(c)
        countt.append(t)
        countg.append(g)
    A = pd.DataFrame(counta)
    C = pd.DataFrame(countc)
    T = pd.DataFrame(countt)
    G = pd.DataFrame(countg)
    frames = [A,C,T,G]
    end = pd.concat(frames,axis = 1)
    aa = end.applymap(lambda x: x / len(seq))
    aa.applymap(lambda x: '%.2f'%x)##保留两位小数
    aa.columns = ['A','C','T','G']
    aa = aa.iloc[3:15,:]
    aa.to_csv(outfile,sep = '\t',index= None)
    return aa



def get_sig_motif_seq(dir):
    """
    单个motif序列
    """
#dir  = '/public/home/shidetong/wedata/meme/out/B6M/loc_2'
#out = '/public/home/shidetong/wedata/meme/seq/loc_1'
#fl = open('/public/home/shidetong/B6M_loc_1.txt','w')
    l = []
    fl = []
    motif = pd.read_csv(dir,sep = '\t')
    motif = motif.drop(motif.tail(3).index)
    lenth = len(motif)
    if lenth == 0:
        g = 0           #判断motif文件是否为空
        fl.append(g)
    elif lenth != 0:
        seq = motif[['matched_sequence']]
        counta = []
        countc = []
        countt = []
        countg = []
        for p in range(20):
            c =0
            a = 0
            t = 0
            g = 0
            for i in seq.iloc[:,:].values:
                q = i[0][p]
                if q == 'T':
                    t = t+1
                elif q == 'A':
                    a = a+1
                elif q == 'C':
                    c= c+1
                elif q == 'G':
                    g = g+1
            counta.append(a)
            countc.append(c)
            countt.append(t)
            countg.append(g)
        A = pd.DataFrame(counta)
        C = pd.DataFrame(countc)
        T = pd.DataFrame(countt)
        G = pd.DataFrame(countg)
        frames = [A,C,T,G]
        end = pd.concat(frames,axis = 1)
        aa = end.applymap(lambda x: x / len(seq))
        aa.columns = ['A','C','T','G']
        seq = aa.idxmax(1)  #取某一行最大值的列名
        seq = seq.tolist()
        seq = ''.join(seq)  #列表的多个元素合并为一个字符
    return seq


def get_motif_seq(dir,out):
    """
    获取两列，第一列为文件信息，第二列为序列信息
    """
#dir  = '/public/home/shidetong/wedata/meme/out/B6M/loc_2'
#out = '/public/home/shidetong/wedata/meme/seq/loc_1'
#fl = open('/public/home/shidetong/B6M_loc_1.txt','w')
    dic = {}
    l = []
    fl = []
    for root_dir,sub_dir,files in os.walk(dir):
        for file in files:
            if file.endswith('.tsv'):
                ro_fl = root_dir
                file_name = os.path.join(root_dir,file)
                motif = pd.read_csv(file_name,sep = '\t')
                motif = motif.drop(motif.tail(3).index)
                lenth = len(motif)
                if lenth == 0:
                    g = 'NO'          #判断motif文件是否为空
                    fl.append(g)
                    l.append(ro_fl)
                elif lenth != 0:
                    seq = motif[['matched_sequence']]
                    counta = []
                    countc = []
                    countt = []
                    countg = []
                    for p in range(20):
                        c =0
                        a = 0
                        t = 0
                        g = 0
                        for i in seq.iloc[:,:].values:
                            q = i[0][p]
                            if q == 'T':
                                t = t+1
                            elif q == 'A':
                                a = a+1
                            elif q == 'C':
                                c= c+1
                            elif q == 'G':
                                g = g+1
                        counta.append(a)
                        countc.append(c)
                        countt.append(t)
                        countg.append(g)
                    A = pd.DataFrame(counta)
                    C = pd.DataFrame(countc)
                    T = pd.DataFrame(countt)
                    G = pd.DataFrame(countg)
                    frames = [A,C,T,G]
                    end = pd.concat(frames,axis = 1)
                    aa = end.applymap(lambda x: x / len(seq))
                    aa.columns = ['A','C','T','G']
                    seq = aa.idxmax(1)  #取某一行最大值的列名
                    seq = seq.tolist()
                    seq = ''.join(seq)  #列表的多个元素合并为一个字符
                    l.append(ro_fl)
                    fl.append(seq)
    dic = {'sit':l,'seq':fl}
    dic = pd.DataFrame(dic)
    dic['new'] = dic['sit'].str.split('/chr',1).str[1]
    dic = dic.drop('sit',axis =1)
    dic = dic[['new','seq']]
    dic['sit'] = dic['new'].str.split('.',1).str[0]
    dic = dic.drop('new',axis = 1)
    dic = dic[['sit','seq']]
    dic.to_csv(out,sep = '\t',index = None)
    return dic 


def base_ratio(loc_file):
    loc1 = pd.read_csv(loc_file,sep = ',')
    loc = loc1.iloc[6:10,:]#选位置
    locT = loc.T
    locT.columns = [0,1,2,3]
    loc = locT[[1]]   #选位置
    loc = loc.dropna(axis=0,how='all')
    counta = []
    countc = []
    countt = []
    countg = []
    for p in range(12):
        c =0
        a = 0
        t = 0
        g = 0
        for i in loc.iloc[:,:].values:  
            q = i[0][p]
            if q == 'T':
                t = t+1
            elif q == 'A':
                a = a+1
            elif q == 'C':
                c= c+1
            elif q == 'G':
                g = g+1
        counta.append(a)
        countc.append(c)
        countt.append(t)
        countg.append(g)
    A = pd.DataFrame(counta)
    C = pd.DataFrame(countc)
    T = pd.DataFrame(countt)
    G = pd.DataFrame(countg)
    frames = [A,C,T,G]
    end = pd.concat(frames,axis = 1)
    aa = end.applymap(lambda x: x / len(loc))
    aa.applymap(lambda x: '%.2f'%x)    #小数转为百分数
    aa.columns = ['A','C','T','G']
    aa.index = ['C','C','A','C','C','A','G','G','G','G','G','C']
    return aa 


#


import pandas as pd 
##通过pandas得到peak位置
def pd_peak(peakfile):
    peak = pd.read_csv(peakfile,skiprows=1,sep = '\t',header = None,usecols=[0,1,2])
    peak.columns = ['chr','start','end']
    a = []
    for i in peak.iloc[:,:].values:
        T= i.tolist()
        a.append(T)
    return a
    


###读取参考基因组某些位置碱基序列
def base_fa(fa_file,peakfile,outfile):
    peak = pd_peak(peakfile)
    bf = pysam.FastaFile(fa_file)
    mo = open(outfile+'Castbiasseq.txt','w')
    for i in peak:
        chro = i[0]
        start = i[1]
        end = i[2]
        seq = bf.fetch(i[0])[i[1]-1:i[2]]
        # mo.write('>'+chro+'_'+start+'_'+end+'\n'+seq)
        mo.write('>'+chro+'_')
    mo.colse()
      



        
    



