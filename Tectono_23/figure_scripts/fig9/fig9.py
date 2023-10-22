import numpy as np
import obspy
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import matplotlib.patches as patches
#########
def great_circ_calc_dist(lat1, lon1, lat2, lon2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2
    c = 2 * asin(sqrt(a))
    km = 6371 * c
    return km
def dist_cal(lat1, lon1, lat2, lon2):
    """
    Determine station distance from cross-section pt
    """
    distance=111*np.sqrt((lat1-lat2)**2+(long1-long2)**2)
    # st.max_P=max_
    return distance
def read_txt(file,station):
    for line in open (file,'r'):
        dt=0
        line=line.split()
        if station==line[4]:
            dt= (float(line[2])-float(line[3]))
            st_lat=float(line[0])
            return dt,st_lat

def read_moho(file,station):
    for line in open (file,'r'):
        line=line.split()
        if station==line[3]:
            moho= float(line[2])
            st_lat=float(line[0])
            return moho,st_lat
## stations_all_single_Pms.txt
## PLOTTING C3
cross_lat=-26
cross_long=133.5


st_C3=['AES01','AES04','AES06','AES14','AES16','AEB04']
st_C2=['AEB15','AEB16','AEB17','AEB18','AEB19','AEB20','AES13','AES15','AES20']
st_C1=['AES07','AES08','AES09','AES11','AES10','AES12','AEB14'] # not taking this.

st_5g=['AES01','AES04','AES16','AEB04','AES20','AEB20']
st_6k=['INGOM','TWINS','MCDOU','MTEBA','BONBO','WHYML','WILGE','BULGU']#'YERDA','COOND',
###
moho_5g_1val=[]
for line in open('moho_oneVal.txt','r'):
    line=line.split()
    if line[3] in st_5g:
        moho_5g_1val.append([float(line[2]),float(line[0]),line[3]])

moho_5g_2val=[]
for line in open('stations_all_moho1_2.txt','r'):
    line=line.split()
    if line[4] in st_5g:
        moho_5g_2val.append([float(line[2]),float(line[3]),float(line[0]),line[4]])

moho_6k_1val=[]
for line in open('6k_stations_Ps.txt','r'):
    line=line.split()
    if line[2] in st_6k:
        moho_6k_1val.append([float(line[3])*8.9,float(line[1]),line[2]])

reflection_moho=[]
for line in open('H-reflect20.xyz','r'):
    line=line.split()
    if -30.5< float(line[1]) <-25 and 133< float(line[2]) <135:
        reflection_moho.append([float(line[3]),float(line[1]),float(line[4])])
bilby_rf=[]
for line in open('H-RecfB.xyz','r'):
    line=line.split()
    if -30.5< float(line[1]) <-24.5 and 132.9< float(line[2]) <135:
        bilby_rf.append([float(line[3]),float(line[1]),float(line[4]),line[0]])
######

fig, ax=plt.subplots(figsize=(12,5))
# plt.style.use('seaborn-whitegrid')
plt.style.use('seaborn-paper')
ax.set_facecolor('whitesmoke')
ax.grid(which='major', axis='both',color='black', linestyle='--',linewidth=.5,
alpha=.25,zorder=1)

#xy, width, height
rect = patches.Rectangle((-25.6, 30), -.55, 30, linewidth=.1, edgecolor='peru',
facecolor='peru',alpha=.2,zorder=2)

rect1 = patches.Rectangle((-27.2, 30), .55,30, linewidth=.1, edgecolor='peru',
facecolor='peru',alpha=.2,zorder=2)

rect2 = patches.Rectangle((-29, 30), .6,30, linewidth=.1, edgecolor='peru',
facecolor='peru',alpha=.2,zorder=2)

# Add the patch to the Axes
ax.add_patch(rect)
ax.add_patch(rect1)
ax.add_patch(rect2)


###
ax.axvline(x=-25.6,ls='--',lw='1.5',c='peru')
ax.axvline(x=-27.2,ls='--',lw='1.5',c='peru')
ax.axvline(x=-29,ls='--',lw='1.5',c='peru')



##
plt.scatter([i[1] for i in reflection_moho],[i[0] for i in reflection_moho],s=95, marker='^',\
c=[i[2] for i in reflection_moho],cmap='Blues',vmin = .4, vmax =.95,
edgecolors='black',linewidths=.25, label='GOMA',zorder=5)

plt.scatter([i[1] for i in bilby_rf],[i[0] for i in bilby_rf],s=95, marker='o',\
color='darkseagreen',vmin = .4, vmax =.85,edgecolors='black',linewidths=.25,label='Bilby',zorder=5)

plt.scatter([i[1] for i in moho_5g_1val],[i[0] for i in moho_5g_1val], marker='+',s=125,\
color='firebrick',alpha=.95, label='5G-6K',zorder=5)

plt.scatter([i[1] for i in moho_6k_1val],[i[0] for i in moho_6k_1val], marker='+',s=125,\
color='firebrick', alpha=.95,zorder=5)

plt.scatter([i[2] for i in moho_5g_2val],[i[0] for i in moho_5g_2val], marker='+',s=125,\
color='firebrick',alpha=.95,zorder=5)
#
plt.scatter([i[2] for i in moho_5g_2val],[i[1] for i in moho_5g_2val], marker='+',s=125,\
color='firebrick',alpha=.95,zorder=5)

plt.ylim(30, 60)
plt.xticks(np.arange(-30.5, -24, .5),fontsize=12)
plt.yticks(np.arange(30, 65, 5),fontsize=12) # was 0,5,.5 for eq
plt.gca().invert_xaxis()
plt.gca().invert_yaxis()

#ax.set_yticks(np.linspace(0,30,7))
plt.xlabel('Latitude',fontsize=14,labelpad=5)
plt.ylabel('Moho depth (km)',fontsize=14)
plt.legend(frameon='True',facecolor='gainsboro')
# fig.text(.368, .8, 'Musgrave Province',fontsize=12, ha='center', va='center',color='indigo')
plt.savefig('Moho_cross_new_.png',bbox_inches='tight', dpi=1000,pad_inches=0.2)
plt.show()
