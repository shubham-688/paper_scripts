# gets data from iris and saves waveforms around earthquake
import matplotlib.pyplot as plt
import obspy
from obspy import read, Stream, UTCDateTime,read_events
from obspy.core.event import Origin
from obspy.taup import TauPyModel
from obspy.geodetics import gps2dist_azimuth
from obspy.clients.fdsn import Client
import numpy as np
import datetime
import os
import sys
from rf import read_rf, rfstats,RFStream
from rf import iter_event_data, IterMultipleComponents
from matplotlib import ticker
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from tqdm import tqdm
from math import cos, sin, asin, sqrt, radians
from matplotlib import cm
# import matplotlib as mpl
from importlib import reload
# reload(hk)
#####################

def calc_distance(lat1, lon1, lat2, lon2):
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

def dist_cal(self,lat,long):
    """
    Determine station distance from cross-section pt
    """
    st = self.stats
    # st.distance=calc_distance(lat,long,st.station_latitude,st.station_longitude)
    st.distance=111*np.sqrt((st.station_latitude-lat)**2+(st.station_longitude-long)**2)
    # st.max_P=max_
    return st.distance
#########################################
# This script uses the refined RF for h-K stacking using functions_hk.
# https://hiperseis.readthedocs.io/en/develop/seismic.receiver_fn.html#module-seismic.receiver_fn.rf_stacking
#########################################
########
# station_5G='OOD'
# for station_5G in ['AES01','AES13','AES15']:
# client = Client("IRIS")
# client_5g=Client('http://auspass.edu.au:80',user='5g',password='grape71')
# inventory = client_5g.get_stations(network="5G",station='*',level='response')

# station_name=[]
# for i in range(len(inventory.get_contents()['stations'])):
    # station_name.append(inventory.get_contents()['stations'][i].split()[0].split('.')[1])

data_stack=[]
file_list=['rf_stack_1Q.h5','rf_stack_5J.h5','rf_stack_7B.h5','rf_stack_AU.h5','rf_stack_1F.h5',\
'rf_stack_3G.h5','rf_stack_6F.h5','rf_stack_7I.h5','rf_stack_S1.h5','rf_stack_1G.h5',\
'rf_stack_5G.h5','rf_stack_6K.h5','rf_stack_7K.h5','rf_stack_YJ.h5']
#

stack_all=RFStream()
for i in range(len(file_list)):
    stack=read(file_list[i],'H5')
    for j in range(len(stack)):
        stack_all.append(stack[j])

print('Len of stack=',len(stack_all))
####
stack_select=RFStream()

station_file=open("P_delay_SA_all_noNan.txt","r")
for line1 in station_file:
    # lat_st=round(float(line1.split()[0]),3)
    # long_st=round(float(line1.split()[1]),3)
    RF_tps=float(line1.split()[2])
    name= line1.split()[3]
    for i in range(len(stack_all)):
        if name==stack_all[i].stats.station:
            stack_all[i].stats.tpdelay=RF_tps
            stack_select.append(stack_all[i])

        # else:
        #     print('Station',name,'not found')

station_file.close()
print('Len of stack=',len(stack_select))
######

#sys.exit()
stack_select.sort(keys=['tpdelay']) # sorts stream wrt to tpdelay
###
dt = UTCDateTime("1970-01-01T00:00:00")
stack_select_try=stack_select.copy()
stack_select_try.trim(dt+22.5,dt+35) # gets 5 sec before and 15 sec after onset
####color map
cmap=cm.plasma# plasma plasma plasma
# cmap=mpl.colormaps['plasma']
#######
########
fig = plt.figure(figsize=(6, 8))
ax1 = fig.add_axes([0.1, 0.07, 0.15, 0.85]) #[left, bottom, width, height]
ax1.set_ylim(59, 0)
ax1.set_xlim(-2.5,10)
ax1.set_yticks(np.linspace(0,59,60))
ax1.tick_params(axis='y',left=False,pad=-2)
# ax1.yaxis.labelpad = -12
ax1.grid(which='both', axis='x',color='LightGrey', linestyle='--',linewidth=.5,alpha=.85)
ax1.xaxis.set_minor_locator(MultipleLocator(2.5))
ax1.xaxis.set_major_locator(MultipleLocator(5))
stack_ax1=stack_select_try[:58]
for i in range(len(stack_ax1)):
    auto=stack_ax1[i]
    time = np.arange(auto.stats.npts) * auto.stats.delta
    time=time-2.5 # shifts onset to 0 sec
    auto.data /= np.max(np.abs(auto.data))
    if np.max(auto.data) == 1:
        auto.data=auto.data*-1
    # dist = dist_cal(auto,cross_lat,cross_long)
    ax1.plot(time,auto.data + i+1,  lw=0.25, color='black')
    ax1.fill_between(time, i+1, auto.data + i+1, lw=0.55,
                      color=cm.plasma(auto.stats.tpdelay/1.32), where=(auto.data < 0),alpha=.8)
    # plt.setp(auto.stats.station,fontsize=3)
#
st_label=[' ']
for i in range(len(stack_ax1)):
    auto=stack_ax1[i]
    st_label.append(auto.stats.station)
    # plt.setp(auto.stats.station,fontsize=3)
plt.setp(ax1.get_yticklabels(),fontsize=4.5)
ax1.yaxis.set_major_formatter(ticker.FixedFormatter(st_label))

########################################## Second column
ax2 = fig.add_axes([0.3, 0.07, 0.15, 0.85]) #[left, bottom, width, height]
ax2.set_ylim(59, 0)
ax2.set_xlim(-2.5,10)
ax2.set_yticks(np.linspace(0,59,60))
ax2.tick_params(axis='y',left=False,labelright=False,pad=-2)
# ax1.yaxis.labelpad = -12
ax2.grid(which='both', axis='x',color='LightGrey', linestyle='--',linewidth=.5,alpha=.85)
ax2.xaxis.set_minor_locator(MultipleLocator(2.5))
ax2.xaxis.set_major_locator(MultipleLocator(5))
#
stack_ax2=stack_select_try[58:116]
for i in range(len(stack_ax2)):
    auto=stack_ax2[i]
    time = np.arange(auto.stats.npts) * auto.stats.delta
    time=time-2.5 # shifts onset to 0 sec
    auto.data /= np.max(np.abs(auto.data))
    if np.max(auto.data) == 1:
        auto.data=auto.data*-1
    # dist = dist_cal(auto,cross_lat,cross_long)
    ax2.plot(time,auto.data + i+1,  lw=0.25, color='black')
    ax2.fill_between(time, i+1, auto.data + i+1, lw=0.55,
                      color=cm.plasma(auto.stats.tpdelay/1.32), where=(auto.data < 0),alpha=.95)
    # plt.setp(auto.stats.station,fontsize=3)
#
st_label=[' ']
for i in range(len(stack_ax2)):
    auto=stack_ax2[i]
    st_label.append(auto.stats.station)
    # plt.setp(auto.stats.station,fontsize=3)
plt.setp(ax2.get_yticklabels(),fontsize=4.5)
ax2.yaxis.set_major_formatter(ticker.FixedFormatter(st_label))
########################################## third column
ax3 = fig.add_axes([0.46, 0.07, 0.15, 0.85]) #[left, bottom, width, height]
ax3.set_ylim(59, 0)
ax3.set_xlim(-2.5,10)
ax3.set_yticks(np.linspace(0,59,60))
ax3.tick_params(axis='y',left=False,labelleft=False,labelright=True,pad=-2)
# ax1.yaxis.labelpad = -12
ax3.grid(which='both', axis='x',color='LightGrey', linestyle='--',linewidth=.5,alpha=.85)
ax3.xaxis.set_minor_locator(MultipleLocator(2.5))
ax3.xaxis.set_major_locator(MultipleLocator(5))
#
stack_ax3=stack_select_try[116:174]
for i in range(len(stack_ax3)):
    auto=stack_ax3[i]
    time = np.arange(auto.stats.npts) * auto.stats.delta
    time=time-2.5 # shifts onset to 0 sec
    auto.data /= np.max(np.abs(auto.data))
    if np.max(auto.data) == 1:
        auto.data=auto.data*-1
    if stack_ax3[i].stats.station == 'CU16' or stack_ax3[i].stats.station == 'E1E1' or stack_ax3[i].stats.station == 'GW08':
        auto.data=auto.data*-1
    if stack_ax3[i].stats.station == 'AEB04' or stack_ax3[i].stats.station == 'CU08' :
        auto.data=auto.data*-1
    # dist = dist_cal(auto,cross_lat,cross_long)
    ax3.plot(time,auto.data + i+1,  lw=0.25, color='black')
    ax3.fill_between(time, i+1, auto.data + i+1, lw=0.55,
                      color=cm.plasma(auto.stats.tpdelay/1.32), where=(auto.data < 0),alpha=.95)
    # plt.setp(auto.stats.station,fontsize=3)
#
st_label=[' ']
for i in range(len(stack_ax3)):
    auto=stack_ax3[i]
    st_label.append(auto.stats.station)
    # plt.setp(auto.stats.station,fontsize=3)
plt.setp(ax3.get_yticklabels(),fontsize=4.5)
ax3.yaxis.set_major_formatter(ticker.FixedFormatter(st_label))
########################################## fourth column
ax4 = fig.add_axes([0.66, 0.07, 0.15, 0.85]) #[left, bottom, width, height]
ax4.set_ylim(59, 0)
ax4.set_xlim(-2.5,10)
ax4.set_yticks(np.linspace(0,59,60))
ax4.tick_params(axis='y',left=False,labelleft=False,labelright=True,pad=-2)
# ax1.yaxis.labelpad = -12
ax4.grid(which='both', axis='x',color='LightGrey', linestyle='--',linewidth=.5,alpha=.85)
ax4.xaxis.set_minor_locator(MultipleLocator(2.5))
ax4.xaxis.set_major_locator(MultipleLocator(5))
#
stack_ax4=stack_select_try[174:231]
for i in range(len(stack_ax4)):
    auto=stack_ax4[i]
    time = np.arange(auto.stats.npts) * auto.stats.delta
    time=time-2.5 # shifts onset to 0 sec
    auto.data /= np.max(np.abs(auto.data))
    if np.max(auto.data) == 1:
        auto.data=auto.data*-1
    if stack_ax4[i].stats.station == 'CU35' or stack_ax4[i].stats.station == 'CU10':
        auto.data=auto.data*-1
    if stack_ax4[i].stats.station == 'AES03' or stack_ax4[i].stats.station == 'CU03' :
        auto.data=auto.data*-1
    # dist = dist_cal(auto,cross_lat,cross_long)
    ax4.plot(time,auto.data + i+1,  lw=0.25, color='black')
    ax4.fill_between(time, i+1, auto.data + i+1, lw=0.55,
                      color=cm.plasma(auto.stats.tpdelay/1.32), where=(auto.data < 0),alpha=.95)
    # plt.setp(auto.stats.station,fontsize=3)
#
st_label=[' ']
for i in range(len(stack_ax4)):
    auto=stack_ax4[i]
    st_label.append(auto.stats.station)
    # plt.setp(auto.stats.station,fontsize=3)
plt.setp(ax4.get_yticklabels(),fontsize=4.5)
ax4.yaxis.set_major_formatter(ticker.FixedFormatter(st_label))
##########################
fig.text(0.45, 0.01, 'Time (s)',fontsize=13, ha='center', va='center')
############
ax6=fig.add_axes([0.25, 0.955, .41, 0.014])
sm = plt.cm.ScalarMappable(norm=plt.Normalize(vmin=0, vmax=1.32),cmap=cm.plasma )
sm.set_array(np.arange(0,1.33))
cbar = plt.colorbar(sm,cax=ax6,orientation='horizontal')
ax6.xaxis.set_minor_locator(MultipleLocator(.1))
ax6.xaxis.set_major_locator(MultipleLocator(.2))
fig.text(0.72, 0.96, 'TPs$_{b}$ (s)',fontsize=10, ha='center', va='center')
fig.text(0.12, 0.025, '(a)',fontsize=9, ha='center', va='center')
fig.text(0.32, 0.025, '(b)',fontsize=9, ha='center', va='center')
fig.text(0.605, 0.025, '(c)',fontsize=9, ha='center', va='center')
fig.text(0.805, 0.025, '(d)',fontsize=9, ha='center', va='center')
# plt.show()
plt.savefig('all_stack_plasma_h_n.pdf',bbox_inches='tight', pad_inches=0.15)
