import matplotlib.pyplot as plt 
import seaborn as sns

# sns.set_palette(sns.color_palette("BuGn_d"))

TAD = [4390,4464,4298]
Loop = [7874,14349,11416]
name = ['H6C7','BxPC3','PANC-1']



#TAD_bar
plt.figure(figsize=(7,5))
plt.bar(range(len(TAD)),TAD,color=['#FB8402','#219EBC','#90C9E6'])
plt.tick_params(labelsize=18)
x = [0,1,2]
lab = ['H6C7','BxPC3','PANC-1']
plt.xticks(x,lab)
plt.savefig('/public/home/shidetong/projects/yf/hic/plot/TAD_bar_newcolor.png',dpi=300)

#Loop_bar
plt.figure(figsize=(7,5))
plt.bar(range(len(Loop)),Loop,color=['#FB8402','#219EBC','#90C9E6'])
plt.tick_params(labelsize=18)
x = [0,1,2]
lab = ['H6C7','BxPC3','PANC-1']
plt.xticks(x,lab)
plt.savefig('/public/home/shidetong/projects/yf/hic/plot/Loop_bar_newcolor.png',dpi=300)