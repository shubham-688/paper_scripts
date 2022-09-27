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
from importlib import reload
# reload(hk)
#####################
# this func adds tpdelay in stack.stats and returns sorted and trimed stack
def get_select_stack(Pdelay_txt,stack_rf):
    stack_select=RFStream()
    station_file=open(Pdelay_txt,"r")
    for line1 in station_file:
        RF_tps=float(line1.split()[2])
        name= line1.split()[3]
        for i in range(len(stack_rf)):
            if name==stack_rf[i].stats.station:
                stack_rf[i].stats.tpdelay=RF_tps
                stack_select.append(stack_rf[i])
    station_file.close()
    print('Len of stack=',len(stack_select))

    stack_select.sort(keys=['tpdelay']) # sorts stream wrt to tpdelay
    ###
    dt = UTCDateTime("1970-01-01T00:00:00")
    # stack_select_try=stack_select.copy()
    stack_select.trim(dt+22.5,dt+35)
    return stack_select

def add_tps_1hz(stack_rf):
    stack_select=RFStream()
    station_file=open("P_delay_5G_1Hz_all.txt","r")
    for line1 in station_file:
        RF_tps=float(line1.split()[2])
        name= line1.split()[3]
        for i in range(len(stack_rf)):
            if name==stack_rf[i].stats.station:
                stack_rf[i].stats.tpdelay_1Hz=RF_tps
                stack_rf[i].stats.tpdelay_1Hzdiff=RF_tps-stack_rf[i].stats.tpdelay
                stack_select.append(stack_rf[i])
    station_file.close()
    stack_select.sort(keys=['tpdelay_1Hz']) # sorts stream wrt to tpdelay from 1 Hz
    return stack_select
#########################################
########

data_stack=[]
# file_list=['rf_stack_5G_1Hz.h5','rf_stack_5G_25Hz.h5','rf_stack_5G_4Hz.h5']

stack_1hz = read_rf('rf_stack_5G_1Hz.h5', 'H5')
stack_25hz = read_rf('rf_stack_5G_25Hz.h5', 'H5')
stack_4hz = read_rf('rf_stack_5G_4Hz.h5', 'H5')

stack_select_1Hz=get_select_stack("P_delay_5G_1Hz_all.txt",stack_1hz)
stack_select_25Hz=get_select_stack("P_delay_5G_25Hz_all.txt",stack_25hz)
stack_select_4Hz=get_select_stack("P_delay_5G_4Hz_all.txt",stack_4hz)

stack_select_25Hz=add_tps_1hz(stack_select_25Hz) # add tps1Hz and sorts based on it
stack_select_4Hz=add_tps_1hz(stack_select_4Hz)

#######
########
fig = plt.figure(figsize=(6, 10))
ax1 = fig.add_axes([0.1, 0.07, 0.15, 0.85]) #[left, bottom, width, height]
ax1.set_ylim(39, 0)
ax1.set_xlim(-2.5,10)
ax1.set_yticks(np.linspace(0,39,40))
ax1.tick_params(axis='y',left=False,length=.1,pad=2)
# ax1.yaxis.labelpad = -12
ax1.grid(which='both', axis='x',color='LightGrey', linestyle='--',linewidth=.5,alpha=.85)
ax1.xaxis.set_minor_locator(MultipleLocator(2.5))
ax1.xaxis.set_major_locator(MultipleLocator(5))
# ax1.set_xlabel('Time (s)')
# ax1.set_xlabel('Distance (km)')
# ax1.annotate('Vertical ACF AEB09', xy=(.6, 1.05), xycoords='axes fraction')
# ax1.set_title('Profile A-A\'',fontsize=12)

stack_ax1=stack_select_1Hz
for i in range(len(stack_ax1)):
    auto=stack_ax1[i]
    time = np.arange(auto.stats.npts) * auto.stats.delta
    time=time-2.5 # shifts onset to 0 sec
    auto.data /= np.max(np.abs(auto.data))
    if np.max(auto.data) == 1:
        auto.data=auto.data*-1
    if stack_ax1[i].stats.station == 'AES03' or stack_ax1[i].stats.station == 'AEB04':
        auto.data=auto.data*-1
    # dist = dist_cal(auto,cross_lat,cross_long)
    ax1.plot(time,auto.data*.5 + i+1,  lw=0.25, color='black')
    # ax1.fill_between(time, i+1, auto.data + i+1, lw=0.55,
                      # color=cm.plasma(auto.stats.tpdelay/1.32), where=(auto.data < 0),alpha=.8)
    ax1.fill_between(time, i+1, auto.data*.5 + i+1, lw=0.55,
                      color='darkorange', where=(auto.data < 0),alpha=.8)
    # plt.setp(auto.stats.station,fontsize=3)
#
st_label=[' ']
tps_label=[' ']
for i in range(len(stack_ax1)):
    auto=stack_ax1[i]
    # label='{}({})'.format(auto.stats.station,int(auto.stats.RF_Pcorr))
    st_label.append(auto.stats.station)
    tps_label.append(auto.stats.tpdelay)
    # plt.setp(auto.stats.station,fontsize=3)
plt.setp(ax1.get_yticklabels(),fontsize=5.5)
ax1.yaxis.set_major_formatter(ticker.FixedFormatter(st_label))
ax_1T = ax1.twinx()
ax_1T.set_ylim(39,0)
ax_1T.set_yticks(np.linspace(0,39,40))
#ax_1T.tick_params(axis='y',left=False,pad=-2)
ax_1T.tick_params(axis='y',left=False,labelleft=False,labelright=True,pad=2,length=.1)
plt.setp(ax_1T.get_yticklabels(),fontsize=5.5)
ax_1T.yaxis.set_major_formatter(ticker.FixedFormatter(tps_label))
ax_1T.set_ylabel('TPs$_{b}$ (s)')
########################################## Second column
ax2 = fig.add_axes([0.325, 0.07, 0.15, 0.85]) #[left, bottom, width, height]
ax2.set_ylim(39, 0)
ax2.set_xlim(-2.5,10)
ax2.set_yticks(np.linspace(0,39,40))
ax2.tick_params(axis='y',left=False,labelright=False,pad=-2)
# ax1.yaxis.labelpad = -12
ax2.grid(which='both', axis='x',color='LightGrey', linestyle='--',linewidth=.5,alpha=.85)
ax2.xaxis.set_minor_locator(MultipleLocator(2.5))
ax2.xaxis.set_major_locator(MultipleLocator(5))
#
# stack_ax2=stack_select_try[58:116]
stack_ax2=stack_select_25Hz

for i in range(len(stack_ax2)):
    auto=stack_ax2[i]
    time = np.arange(auto.stats.npts) * auto.stats.delta
    time=time-2.5 # shifts onset to 0 sec
    auto.data /= np.max(np.abs(auto.data))
    if np.max(auto.data) == 1:
        auto.data=auto.data*-1
    # dist = dist_cal(auto,cross_lat,cross_long)
    ax2.plot(time,auto.data*.5 + i+1,  lw=0.25, color='black')
    ax2.fill_between(time, i+1, auto.data*.5 + i+1, lw=0.55,
                      color='olivedrab', where=(auto.data < 0),alpha=.8)
    # plt.setp(auto.stats.station,fontsize=3)
#
st_label=[' ']
for i in range(len(stack_ax2)):
    auto=stack_ax2[i]
    # label='{}({})'.format(auto.stats.station,int(auto.stats.RF_Pcorr))
    # st_label.append(label)
    # plt.setp(auto.stats.station,fontsize=3)
plt.setp(ax2.get_yticklabels(),fontsize=4.5)
ax2.yaxis.set_major_formatter(ticker.FixedFormatter(st_label))
########################################## third column
ax3 = fig.add_axes([0.525, 0.07, 0.15, 0.85]) #[left, bottom, width, height]
ax3.set_ylim(39, 0)
ax3.set_xlim(-2.5,10)
ax3.set_yticks(np.linspace(0,39,40))
ax3.tick_params(axis='y',left=False,labelleft=False,labelright=True,pad=-2)
# ax1.yaxis.labelpad = -12
ax3.grid(which='both', axis='x',color='LightGrey', linestyle='--',linewidth=.5,alpha=.85)
ax3.xaxis.set_minor_locator(MultipleLocator(2.5))
ax3.xaxis.set_major_locator(MultipleLocator(5))
#
# stack_ax3=stack_select_try[116:174]
stack_ax3=stack_select_4Hz

for i in range(len(stack_ax3)):
    auto=stack_ax3[i]
    time = np.arange(auto.stats.npts) * auto.stats.delta
    time=time-2.5 # shifts onset to 0 sec
    auto.data /= np.max(np.abs(auto.data))
    if np.max(auto.data) == 1:
        auto.data=auto.data*-1
    # dist = dist_cal(auto,cross_lat,cross_long)
    ax3.plot(time,auto.data*.5 + i+1,  lw=0.25, color='black')
    ax3.fill_between(time, i+1, auto.data*.5 + i+1, lw=0.55,
                      color='palevioletred', where=(auto.data < 0),alpha=.75)
    # plt.setp(auto.stats.station,fontsize=3)
#
tps_label=[' ']
for i in range(len(stack_ax3)):
    auto=stack_ax3[i]
    # label='{}({})'.format(auto.stats.station,int(auto.stats.RF_Pcorr))
    tps_label.append(auto.stats.tpdelay)
    # plt.setp(auto.stats.station,fontsize=3)
plt.setp(ax3.get_yticklabels(),fontsize=5.5)
ax3.yaxis.set_major_formatter(ticker.FixedFormatter([' ']))
##########################
######################################### third column
ax4 = fig.add_axes([0.725, 0.07, 0.15, 0.85]) #[left, bottom, width, height]
ax4.set_ylim(39, 0)
ax4.set_xlim(-.3,.3)
ax4.set_yticks(np.linspace(0,39,40))
ax4.tick_params(axis='y',left=False,labelleft=False,labelright=True,pad=-1)
# ax1.yaxis.labelpad = -12
ax4.grid(which='minor', axis='x',color='Grey', linestyle='--',linewidth=.75,alpha=.85)
ax4.xaxis.set_minor_locator(MultipleLocator(.125))
ax4.xaxis.set_major_locator(MultipleLocator(.25))
#
for i in range(len(stack_ax3)):

    ax4.plot(stack_ax3[i].stats.tpdelay_1Hzdiff,i+1,'o-', ms=5.5,color='mediumvioletred',alpha=.5)
    ax4.plot(stack_ax2[i].stats.tpdelay_1Hzdiff,i+1,'o-', ms=5.5,color='darkolivegreen',alpha=.4)

plt.xticks(fontsize=8)
st_label=[' ']
for i in range(len(stack_ax2)):
    auto=stack_ax2[i]
    # label='{}({})'.format(auto.stats.station,int(auto.stats.RF_Pcorr))
    st_label.append(auto.stats.station)
plt.setp(ax4.get_yticklabels(),fontsize=5.5)
ax4.yaxis.set_major_formatter(ticker.FixedFormatter(st_label))
# stack_ax3=stack_select_try[116:174]
##########################
fig.text(0.475, 0.025, 'Time (s)',fontsize=13, ha='center', va='center')
###########
# ax5= fig.add_axes([0.87, 0.35, 0.025, 0.3]) #[left, bottom, width, height]
# sm = plt.cm.ScalarMappable(norm=plt.Normalize(vmin=0, vmax=1.32),cmap=cm.plasma )
# sm.set_array(np.arange(0,1.33))
# cbar = plt.colorbar(sm,cax=ax5)
# # cbar.set_label('TPs (s)')
# fig.text(0.9, 0.665, 'TPsb (s)',fontsize=10, ha='center', va='center')
# plt.colorbar()
############
# ax6=fig.add_axes([0.25, 0.955, .41, 0.014])
# sm = plt.cm.ScalarMappable(norm=plt.Normalize(vmin=0, vmax=1.32),cmap=cm.plasma )
# sm.set_array(np.arange(0,1.33))
# cbar = plt.colorbar(sm,cax=ax6,orientation='horizontal')
# ax6.xaxis.set_minor_locator(MultipleLocator(.1))
# ax6.xaxis.set_major_locator(MultipleLocator(.2))
fig.text(0.175, 0.935, '0.1-1 Hz',fontsize=12, color='darkorange',ha='center', va='center')
fig.text(0.4, 0.935, '0.1-2.5 Hz',fontsize=12, color='olivedrab',ha='center', va='center')
fig.text(0.595, 0.935, '0.1-4 Hz',fontsize=12, color='palevioletred',ha='center', va='center')
fig.text(0.785, 0.935, '\u03B4 TPs$_{b}$',fontsize=11, color='grey',ha='center', va='center')

fig.text(0.175, 0.96, '(a)',fontsize=13, ha='center', va='center')
fig.text(0.4, 0.96, '(b)',fontsize=13, ha='center', va='center')
fig.text(0.6, 0.96, '(c)',fontsize=13, ha='center', va='center')
fig.text(0.785, 0.96, '(d)',fontsize=13, ha='center', va='center')

plt.savefig('stack_freq_all.png',bbox_inches='tight', pad_inches=0.15,dpi=300)
#plt.show()
