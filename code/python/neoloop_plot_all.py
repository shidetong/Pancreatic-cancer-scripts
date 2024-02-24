#%%
from neoloop.visualize.core import *
import cooler
import numpy as np
from numpy.lib.function_base import median 
import pandas as pd 
#%%
# clr = cooler.Cooler('/public/home/shidetong/projects/yf/hic/test/hic_breakfinder/PANC/inter_10000.cool')
# #%%
# file = open ('/public/home/shidetong/projects/yf/hic/test/hic_breakfinder/PANC/panc.assemblies.txt','r')



# lines = file.readlines()
# len(lines)
# #for i in range(len(lines)):
# #    line = lines[i]
# #    print(line)


# for i in range(len(lines)) :
#     line = lines[i]
#     assembly = line   
# #assembly = 'C8      translocation,10,16000000,+,7,43260000,-        10,15050000     7,44120000'

#     vis = Triangle(clr, assembly, n_rows=9, figsize=(7, 7), track_partition=[5, 1,1,1,1,1,1,1, 0.5], correct='weight')

#     vis.matrix_plot(vmin=0)
#     vis.plot_chromosome_bounds(linewidth=2.5)

#     #vis.plot_loops('/public/home/shidetong/projects/yf/hic/test/hic_breakfinder/20200318Bxpc/Bxpc.neo-loops.txt', face_color='none', marker_size=40, cluster=True)


#     vis.plot_signal('RNA_Panc', '/public/home/shidetong/projects/yf/RNA-seq/20200924/bw/Panc.bw', label_size=10, data_range_size=5, max_value=50, color='#E31A1C')
#     vis.plot_signal('CTCF', '/public/home/shidetong/projects/yf/chip-seq/chip_2/bw/Panc_CTCF-1.bw', label_size=10, data_range_size=5, max_value=50, color='#6A3D9A')
#     vis.plot_signal('H3K27ac', '/public/home/shidetong/projects/yf/chip-seq/20201223/Panc/bw/1223PANCH3K27.bw', label_size=10, data_range_size=5, max_value=50, color='#4B0082')
#     vis.plot_signal('H3K4', '/public/home/shidetong/projects/yf/chip-seq/20201223/Panc/bw/1223PANCH3K4.bw', label_size=10, data_range_size=5, max_value=50, color='#7CFC00')
#     vis.plot_signal('PolII', '/public/home/shidetong/projects/yf/chip-seq/20201223/Panc/bw/1223PANCPolII.bw', label_size=10, data_range_size=5, max_value=50, color='#FF4500')
#     vis.plot_signal('ATAC', '/public/home/shidetong/projects/yf/ATAC/bw/Panc_ATAC_1.bw', label_size=10, data_range_size=5, max_value=50, color='#008000')
#     vis.plot_signal('Rad21', '/public/home/shidetong/projects/yf/chip-seq/chip_2/bw/Panc_Rad21.bw', label_size=10, data_range_size=5, max_value=50, color='#6A3D9A')
#     #vis.plot_genes(filter_=['PRAME','BCRP4', 'RAB36', 'BCR', 'ABL1', 'NUP214'],label_aligns={'PRAME':'right','RAB36':'right'}, fontsize=9)
#     vis.plot_chromosome_bar(name_size=11, coord_size=4.8)

#         #vis.show()
#     # p=i-1
#     vis.outfig('/public/home/shidetong/projects/yf/Plot/new_transplot/Panc/' + 'panc' + '_'+str(i) + '.jpg')
# %%





# ####bxpc
# clr = cooler.Cooler('/public/home/shidetong/projects/yf/hic/test/hic_breakfinder/20200318Bxpc/inter_10000.cool')
# #%%
# file = open ('/public/home/shidetong/projects/yf/hic/test/hic_breakfinder/20200318Bxpc/Bxpc.assemblies.txt','r')



# lines = file.readlines()
# len(lines)
# #for i in range(len(lines)):
# #    line = lines[i]
# #    print(line)


# for i in range(len(lines)) :
#     line = lines[i]
#     assembly = line   
# #assembly = 'C8      translocation,10,16000000,+,7,43260000,-        10,15050000     7,44120000'

#     vis = Triangle(clr, assembly, n_rows=9, figsize=(7, 7), track_partition=[5, 1,1,1,1,1,1,1, 0.5], correct='weight')

#     vis.matrix_plot(vmin=0)
#     vis.plot_chromosome_bounds(linewidth=2.5)

#     #vis.plot_loops('/public/home/shidetong/projects/yf/hic/test/hic_breakfinder/20200318Bxpc/Bxpc.neo-loops.txt', face_color='none', marker_size=40, cluster=True)


#     vis.plot_signal('RNA', '/public/home/shidetong/projects/yf/RNA-seq/20200924/bw/Bxpc.bw', label_size=10, data_range_size=5, max_value=50, color='#E31A1C')
#     vis.plot_signal('CTCF', '/public/home/shidetong/projects/yf/chip-seq/chip_2/bw/Bxpc_CTCF-1.bw', label_size=10, data_range_size=5, max_value=50, color='#6A3D9A')
#     vis.plot_signal('H3K27ac', '/public/home/shidetong/projects/yf/chip-seq/20201223/Bxpc/bw/1223BxpcH3K27.bw', label_size=10, data_range_size=5, max_value=50, color='#4B0082')
#     vis.plot_signal('H3K4', '/public/home/shidetong/projects/yf/chip-seq/20201223/Bxpc/bw/1223BxpcH3K4.bw', label_size=10, data_range_size=5, max_value=50, color='#7CFC00')
#     vis.plot_signal('PolII', '/public/home/shidetong/projects/yf/chip-seq/20201223/Bxpc/bw/1223BxpcPolII.bw', label_size=10, data_range_size=5, max_value=50, color='#FF4500')
#     vis.plot_signal('ATAC', '/public/home/shidetong/projects/yf/ATAC/bw/Bxpc_ATAC_1.bw', label_size=10, data_range_size=5, max_value=50, color='#008000')
#     vis.plot_signal('Rad21', '/public/home/shidetong/projects/yf/chip-seq/chip_2/bw/Bxpc_Rad21.bw', label_size=10, data_range_size=5, max_value=50, color='#6A3D9A')
#     #vis.plot_genes(filter_=['PRAME','BCRP4', 'RAB36', 'BCR', 'ABL1', 'NUP214'],label_aligns={'PRAME':'right','RAB36':'right'}, fontsize=9)
#     vis.plot_chromosome_bar(name_size=11, coord_size=4.8)

#         #vis.show()
#     # p=i-1
#     vis.outfig('/public/home/shidetong/projects/yf/Plot/new_transplot/bxpc/' + 'bxpc' + '_'+str(i) + '.jpg')




# ####panc_h6c7
# clr = cooler.Cooler('/public/home/shidetong/projects/yf/hic/test/hic_breakfinder/20200330H6C7-2/inter_10000.cool')
# #%%
# file = open ('/public/home/shidetong/projects/yf/hic/test/hic_breakfinder/PANC/panc.assemblies.txt','r')



# lines = file.readlines()
# len(lines)
# #for i in range(len(lines)):
# #    line = lines[i]
# #    print(line)


# for i in range(len(lines)) :
#     line = lines[i]
#     assembly = line   
# #assembly = 'C8      translocation,10,16000000,+,7,43260000,-        10,15050000     7,44120000'

#     vis = Triangle(clr, assembly, n_rows=9, figsize=(7, 7), track_partition=[5, 1,1,1,1,1,1,1, 0.5], correct='weight')

#     vis.matrix_plot(vmin=0)
#     vis.plot_chromosome_bounds(linewidth=2.5)

#     #vis.plot_loops('/public/home/shidetong/projects/yf/hic/test/hic_breakfinder/20200318Bxpc/Bxpc.neo-loops.txt', face_color='none', marker_size=40, cluster=True)


#     vis.plot_signal('RNA', '/public/home/shidetong/projects/yf/RNA-seq/20200924/bw/H6C7.bw', label_size=10, data_range_size=5, max_value=50, color='#E31A1C')
#     vis.plot_signal('CTCF', '/public/home/shidetong/projects/yf/chip-seq/chip_2/bw/H6C7_CTCF-1.bw', label_size=10, data_range_size=5, max_value=50, color='#6A3D9A')
#     vis.plot_signal('H3K27ac', '/public/home/shidetong/projects/yf/chip-seq/20210107/H6C7/bw/1217H6C7H3K27.bw', label_size=10, data_range_size=5, max_value=50, color='#4B0082')
#     vis.plot_signal('H3K4', '/public/home/shidetong/projects/yf/chip-seq/20210107/H6C7/bw/0107H6C7H3K4.bw', label_size=10, data_range_size=5, max_value=50, color='#7CFC00')
#     vis.plot_signal('PolII', '/public/home/shidetong/projects/yf/chip-seq/20210107/H6C7/bw/0107H6C7PolII.bw', label_size=10, data_range_size=5, max_value=50, color='#FF4500')
#     vis.plot_signal('ATAC', '/public/home/shidetong/projects/yf/ATAC/bw/H6C7_ATAC_3.bw', label_size=10, data_range_size=5, max_value=50, color='#008000')
#     vis.plot_signal('Rad21', '/public/home/shidetong/projects/yf/chip-seq/chip_2/bw/H6C7_Rad21.bw', label_size=10, data_range_size=5, max_value=50, color='#6A3D9A')
#     #vis.plot_genes(filter_=['PRAME','BCRP4', 'RAB36', 'BCR', 'ABL1', 'NUP214'],label_aligns={'PRAME':'right','RAB36':'right'}, fontsize=9)
#     vis.plot_chromosome_bar(name_size=11, coord_size=4.8)

#         #vis.show()
#     # p=i-1
#     vis.outfig('/public/home/shidetong/projects/yf/Plot/new_transplot/Panc/' + 'h6c7' + '_'+str(i) + '.jpg')



    
# ####bxpc_h6c7
# clr = cooler.Cooler('/public/home/shidetong/projects/yf/hic/test/hic_breakfinder/20200330H6C7-2/inter_10000.cool')
# #%%
# file = open ('/public/home/shidetong/projects/yf/hic/test/hic_breakfinder/20200318Bxpc/Bxpc.assemblies.txt','r')



# lines = file.readlines()
# len(lines)
# #for i in range(len(lines)):
# #    line = lines[i]
# #    print(line)


# for i in range(len(lines)) :
#     line = lines[i]
#     assembly = line   
# #assembly = 'C8      translocation,10,16000000,+,7,43260000,-        10,15050000     7,44120000'

#     vis = Triangle(clr, assembly, n_rows=9, figsize=(7, 7), track_partition=[5, 1,1,1,1,1,1,1, 0.5], correct='weight')

#     vis.matrix_plot(vmin=0)
#     vis.plot_chromosome_bounds(linewidth=2.5)

#     #vis.plot_loops('/public/home/shidetong/projects/yf/hic/test/hic_breakfinder/20200318Bxpc/Bxpc.neo-loops.txt', face_color='none', marker_size=40, cluster=True)


#     vis.plot_signal('RNA', '/public/home/shidetong/projects/yf/RNA-seq/20200924/bw/H6C7.bw', label_size=10, data_range_size=5, max_value=50, color='#E31A1C')
#     vis.plot_signal('CTCF', '/public/home/shidetong/projects/yf/chip-seq/chip_2/bw/H6C7_CTCF-1.bw', label_size=10, data_range_size=5, max_value=50, color='#6A3D9A')
#     vis.plot_signal('H3K27ac', '/public/home/shidetong/projects/yf/chip-seq/20210107/H6C7/bw/1217H6C7H3K27.bw', label_size=10, data_range_size=5, max_value=50, color='#4B0082')
#     vis.plot_signal('H3K4', '/public/home/shidetong/projects/yf/chip-seq/20210107/H6C7/bw/0107H6C7H3K4.bw', label_size=10, data_range_size=5, max_value=50, color='#7CFC00')
#     vis.plot_signal('PolII', '/public/home/shidetong/projects/yf/chip-seq/20210107/H6C7/bw/0107H6C7PolII.bw', label_size=10, data_range_size=5, max_value=50, color='#FF4500')
#     vis.plot_signal('ATAC', '/public/home/shidetong/projects/yf/ATAC/bw/H6C7_ATAC_3.bw', label_size=10, data_range_size=5, max_value=50, color='#008000')
#     vis.plot_signal('Rad21', '/public/home/shidetong/projects/yf/chip-seq/chip_2/bw/H6C7_Rad21.bw', label_size=10, data_range_size=5, max_value=50, color='#6A3D9A')
#     #vis.plot_genes(filter_=['PRAME','BCRP4', 'RAB36', 'BCR', 'ABL1', 'NUP214'],label_aligns={'PRAME':'right','RAB36':'right'}, fontsize=9)
#     vis.plot_chromosome_bar(name_size=11, coord_size=4.8)

#         #vis.show()
#     # p=i-1
#     vis.outfig('/public/home/shidetong/projects/yf/Plot/new_transplot/bxpc/' + 'h6c7' + '_'+str(i) + '.jpg')







#test------------------------------------------------
####panc_h6c7
# clr1 = cooler.Cooler('/public/home/shidetong/projects/yf/hic/test/hic_breakfinder/20200330H6C7-2/inter_10000.cool')
clr1 = cooler.Cooler('/public/home/shidetong/projects/yf/hic/test/hic_breakfinder/20200318Bxpc/inter_10000.cool')
# clr2 = cooler.Cooler('/public/home/shidetong/projects/yf/hic/test/hic_breakfinder/PANC/inter_10000.cool')
#%%
# file = open ('/public/home/shidetong/projects/yf/hic/test/hic_breakfinder/PANC/panc.assemblies.txt','r')
file = open ('/public/home/shidetong/projects/yf/hic/test/hic_breakfinder/20200318Bxpc/Bxpc.assemblies.txt','r')


# lines = file.readlines()
# len(lines)
# #for i in range(len(lines)):
# #    line = lines[i]
# #    print(line)


# for i in range(len(lines)) :
#     line = lines[i]
#     assembly = line   
# #assembly = 'C8      translocation,10,16000000,+,7,43260000,-        10,15050000     7,44120000'

#     vis = Triangle(clr1, assembly, n_rows=10, figsize=(7, 14), track_partition=[5,5, 1,1,1,1,1,1,1, 0.5], correct='weight')
#     vis = Triangle(clr2, assembly, n_rows=10, figsize=(7, 14), track_partition=[5,5, 1,1,1,1,1,1,1, 0.5], correct='weight')
#     vis.matrix_plot(vmin=0)
#     vis.plot_chromosome_bounds(linewidth=2.5)
#     # vis = Triangle(clr2, assembly, n_rows=10, figsize=(7, 14), track_partition=[5,5, 1,1,1,1,1,1,1, 0.5], correct='weight')
#     # vis.matrix_plot(vmin=0)
#     # vis.plot_chromosome_bounds(linewidth=2.5)



#     #vis.plot_loops('/public/home/shidetong/projects/yf/hic/test/hic_breakfinder/20200318Bxpc/Bxpc.neo-loops.txt', face_color='none', marker_size=40, cluster=True)


#     vis.plot_signal('RNA', '/public/home/shidetong/projects/yf/RNA-seq/20200924/bw/H6C7.bw', label_size=10, data_range_size=5, max_value=50, color='#E31A1C')
#     vis.plot_signal('CTCF', '/public/home/shidetong/projects/yf/chip-seq/chip_2/bw/H6C7_CTCF-1.bw', label_size=10, data_range_size=5, max_value=50, color='#6A3D9A')
#     vis.plot_signal('H3K27ac', '/public/home/shidetong/projects/yf/chip-seq/20210107/H6C7/bw/1217H6C7H3K27.bw', label_size=10, data_range_size=5, max_value=50, color='#4B0082')
#     vis.plot_signal('H3K4', '/public/home/shidetong/projects/yf/chip-seq/20210107/H6C7/bw/0107H6C7H3K4.bw', label_size=10, data_range_size=5, max_value=50, color='#7CFC00')
#     vis.plot_signal('PolII', '/public/home/shidetong/projects/yf/chip-seq/20210107/H6C7/bw/0107H6C7PolII.bw', label_size=10, data_range_size=5, max_value=50, color='#FF4500')
#     vis.plot_signal('ATAC', '/public/home/shidetong/projects/yf/ATAC/bw/H6C7_ATAC_3.bw', label_size=10, data_range_size=5, max_value=50, color='#008000')
#     vis.plot_signal('Rad21', '/public/home/shidetong/projects/yf/chip-seq/chip_2/bw/H6C7_Rad21.bw', label_size=10, data_range_size=5, max_value=50, color='#6A3D9A')
#     #vis.plot_genes(filter_=['PRAME','BCRP4', 'RAB36', 'BCR', 'ABL1', 'NUP214'],label_aligns={'PRAME':'right','RAB36':'right'}, fontsize=9)
#     vis.plot_chromosome_bar(name_size=11, coord_size=4.8)

#         #vis.show()
#     # p=i-1
#     vis.outfig('/public/home/shidetong/projects/yf/Plot/test/' + 'h6c7' + '_'+str(i) + '.jpg')





lines = file.readlines()
len(lines)
#for i in range(len(lines)):
#    line = lines[i]
#    print(line)



line = lines[54]
assembly = line   
#assembly = 'C8      translocation,10,16000000,+,7,43260000,-        10,15050000     7,44120000'

# vis = Triangle(clr1, assembly, n_rows=10, figsize=(7, 14), track_partition=[5,5, 1,1,1,1,1,1,1, 0.5], correct='weight')
# vis = Triangle(clr2, assembly, n_rows=14, figsize=(7, 7), track_partition=[5,0.4,0.4,0.4,0.4,0.4,1,1,1,1,1,1,1,0.5], correct='weight')
# vis = Triangle(clr2, assembly, n_rows=10, figsize=(7, 7), track_partition=[5,0.4,1,1,1,1,1,1,1,0.5], correct='weight')
vis = Triangle(clr1, assembly, n_rows=4, figsize=(7, 3), track_partition=[5,0.4,1,0.5], correct='weight')

vis.matrix_plot(vmin=0)
vis.plot_chromosome_bounds(linewidth=2.5)
# vis = Triangle(clr2, assembly, n_rows=10, figsize=(7, 14), track_partition=[5,5, 1,1,1,1,1,1,1, 0.5], correct='weight')
# vis.matrix_plot(vmin=0)
# vis.plot_chromosome_bounds(linewidth=2.5)



#vis.plot_loops('/public/home/shidetong/projects/yf/hic/test/hic_breakfinder/20200318Bxpc/Bxpc.neo-loops.txt', face_color='none', marker_size=40, cluster=True)
# #PANC_10
# vis.plot_genes(filter_=['SSPN', 'RASSF8'],label_aligns={'SSPN':'right','RASSF8':'left'}, fontsize=9)
# vis.plot_genes(filter_=['ITPR2', 'BHLHE41'],label_aligns={'ITPR2':'right','BHLHE41':'left'}, fontsize=9)
# vis.plot_genes(filter_=[ 'AC022509.3','AC022509.1'],label_aligns={'AC022509.1':'right','AC022509.3':'left'}, fontsize=9)
# vis.plot_genes(filter_=[ 'AC022509.5','AC022509.2'],label_aligns={'AC022509.2':'right','AC022509.5':'left'}, fontsize=9)
# vis.plot_genes(filter_=[ 'AC024145.1','AC055720.2'],label_aligns={'AC055720.2':'right','AC024145.1':'left'}, fontsize=9)
#PANC_19

vis.plot_genes(filter_=['ABCG1','TFF2','MX2'],label_aligns={'ABCG1':'right','TFF2':'left','MX2':'left'}, fontsize=9)
# vis.plot_genes(filter_=['HOXB3', 'HOXB4'],label_aligns={'HOXB3':'right','HOXB4':'left'}, fontsize=9)
# vis.plot_genes(filter_=[ 'HOXB5','HOXB9'],label_aligns={'HOXB5':'right','HOXB9':'left'}, fontsize=9)
# vis.plot_genes(filter_=[ 'AC103702.2','AC091179.2'],label_aligns={'AC103702.2':'left','AC091179.2':'right'}, fontsize=9)
# vis.plot_genes(filter_=[ 'RPL9P28','AC055720.2'],label_aligns={'AC055720.2':'right','AC024145.1':'left'}, fontsize=9)
# #PANC_23
# vis.plot_genes(filter_=['RASSF8', 'SSPN'],label_aligns={'RASSF8':'left','SSPN':'right'}, fontsize=9)
# # vis.plot_genes(filter_=['FBXO15', 'TIMM21'],label_aligns={'TIMM21':'left','FBXO15':'right'}, fontsize=9)
# # vis.plot_genes(filter_=['CYB5A', 'FAUP1'],label_aligns={'CYB5A':'right','FAUP1':'left'}, fontsize=9)
# vis.plot_genes(filter_=[ 'DIPK1C','CNDP2'],label_aligns={'DIPK1C':'right','CNDP2':'left'}, fontsize=9)
# # vis.plot_genes(filter_=[ 'CNDP1','ZNF407'],label_aligns={'CNDP1':'right','ZNF407':'left'}, fontsize=9)
# # vis.plot_genes(filter_=[ 'AC024145.1','AC055720.2'],label_aligns={'AC055720.2':'right','AC024145.1':'left'}, fontsize=9)


vis.plot_signal('RNA', '/public/home/shidetong/projects/yf/RNA-seq/20200924/bw/Bxpc.bw', label_size=10, data_range_size=5, max_value=50, color='#E31A1C')


# vis.plot_signal('RNA_Panc', '/public/home/shidetong/projects/yf/RNA-seq/20200924/bw/Panc.bw', label_size=10, data_range_size=5, max_value=50, color='#E31A1C')
# vis.plot_signal('CTCF', '/public/home/shidetong/projects/yf/chip-seq/chip_2/bw/Panc_CTCF-1.bw', label_size=10, data_range_size=5, max_value=50, color='#6A3D9A')
# vis.plot_signal('H3K27ac', '/public/home/shidetong/projects/yf/chip-seq/20201223/Panc/bw/1223PANCH3K27.bw', label_size=10, data_range_size=5, max_value=50, color='#4B0082')
# vis.plot_signal('H3K4', '/public/home/shidetong/projects/yf/chip-seq/20201223/Panc/bw/1223PANCH3K4.bw', label_size=10, data_range_size=5, max_value=50, color='#7CFC00')
# vis.plot_signal('PolII', '/public/home/shidetong/projects/yf/chip-seq/20201223/Panc/bw/1223PANCPolII.bw', label_size=10, data_range_size=5, max_value=50, color='#FF4500')
# vis.plot_signal('ATAC', '/public/home/shidetong/projects/yf/ATAC/bw/Panc_ATAC_1.bw', label_size=10, data_range_size=5, max_value=50, color='#008000')
# vis.plot_signal('Rad21', '/public/home/shidetong/projects/yf/chip-seq/chip_2/bw/Panc_Rad21.bw', label_size=10, data_range_size=5, max_value=50, color='#6A3D9A')
#vis.plot_genes(filter_=['PRAME','BCRP4', 'RAB36', 'BCR', 'ABL1', 'NUP214'],label_aligns={'PRAME':'right','RAB36':'right'}, fontsize=9)
#vis.plot_genes(filter_=['PRAME','BCRP4', 'RAB36', 'BCR', 'ABL1', 'NUP214'],label_aligns={'PRAME':'right','RAB36':'right'}, fontsize=9)
vis.plot_chromosome_bar(name_size=11, coord_size=4.8)

    #vis.show()
# p=i-1
# vis.outfig('/public/home/shidetong/projects/yf/Plot/test/' + 'h6c7' + '_'+str(i) + '.jpg')

vis.outfig('/public/home/shidetong/projects/yf/Plot/test/test_H6C7_54_fig6B.jpg',dpi=300)


# %%
