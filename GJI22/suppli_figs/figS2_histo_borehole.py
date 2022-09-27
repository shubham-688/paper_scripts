###
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
import sys
# from brokenaxes import brokenaxes
from mpl_toolkits.axes_grid1 import make_axes_locatable
#plt.ion()

file=open('/Users/Shubham/Research/geological_data_straya/GSSA_basement/cover_thickness_2020.txt','r')
lines=file.readlines()
n_line=len(lines)
x=[]
y=[]
z=[]
for i in range(1,n_line):
    split_line=lines[i].split()
    xyz_t=[]
    x.append(float(split_line[0].rstrip()))
    y.append(float(split_line[1].rstrip()))
    z.append(float(split_line[2].rstrip()))

file.close()

z=np.array(z)

# fig=plt.subplots(figsize=(10,7))
# plt.style.use('seaborn-paper')
# # plt.style.use('seaborn')
# n, bins, patches = plt.hist(z, 41, rwidth=.75, facecolor='DimGray', alpha=0.75)
# plt.ylabel('Number of borehole measurements',fontfamily='serif',fontsize=15)
# plt.grid(True)
# plt.show()

# bin=[0,150,300,450,600,750,1000,2000,3000,4000]
bin_edge=np.linspace(0,4000,21)
bin_centre=np.linspace(100,3900,20)
my_hist = np.histogram(z, bins = bin_edge)[0]

f, (ax, ax2) = plt.subplots(2,1,sharex = True, facecolor = 'w',figsize=(7,5)) # make the axes
###
ax.grid(color='cadetblue', linewidth=.4,alpha=.25)
ax2.grid(color='cadetblue', linewidth=.4,alpha=.25)
ax.bar(bin_centre, my_hist,width=145,color='DarkGray',alpha=.85) # plot on top axes
ax2.bar(bin_centre, my_hist,width=145,color='DarkGray',alpha=.85) # plot on bottom axes
###
ax.set_ylim([25501,25750]) # numbers here are specific to this example
ax2.set_ylim([0, 251]) # numbers here are specific to this example
###
ax.spines['bottom'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax.yaxis.tick_left()
# ax.tick_params(labeltop='off')
ax2.set_xticks(bin_centre)
ax2.set_yticks(np.linspace(0,200,5))
ax.set_yticks(np.linspace(25550,25700,4))
# ax.xaxis.set_visible(False)
plt.subplots_adjust(hspace=0.08)
# ax.grid(color='cadetblue', linewidth=.45,alpha=.35)
# ax2.grid(color='cadetblue', linewidth=.45,alpha=.35)
####
f.text(0.027, 0.5, 'No. of boreholes',fontfamily='serif',fontsize=13, ha='center', va='center', rotation='vertical')
# plt.ylabel('# borehole measurements',fontfamily='serif',fontsize=12)
f.text(0.5, 0.001, 'Basement depth (m)',fontfamily='serif',fontsize=13, ha='center', va='center')

# plt.xlabel('Basement depth (m)',fontfamily='serif',fontsize=12)
labels = ax2.get_xticklabels()
plt.setp(labels, rotation=40, ha="right",rotation_mode="anchor")
ax.tick_params(bottom=False, labelbottom=False)
####
d = .013 # how big to make the diagonal lines in axes coordinates
# arguments to pass plot, just so we don't keep repeating them
kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
ax.plot((1-d,1+d),(-d,+d), **kwargs) # top-left diagonal
ax.plot((-d,+d),(-d,d), **kwargs) # bottom-left diagonal
#
kwargs.update(transform=ax2.transAxes) # switch to the bottom axes
ax2.plot((1-d,1+d),(1-d,1+d), **kwargs) # top-right diagonal
ax2.plot((-d,d),(1-d,1+d), **kwargs) # bottom-right diagonal
plt.savefig('BH_histo.pdf',bbox_inches='tight', pad_inches=0.2)
plt.show()
