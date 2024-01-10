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
def get_autoC(r0,Dt,len,samp_rate):
    resonance_filter = np.zeros(len)
    resonance_filter[0] = 1
    resonance_filter[int(Dt * samp_rate)] = r0
    return resonance_filter

def get_decay_cos(r0,Dt,time):
    m_t=-r0**(time/Dt)*np.cos(np.pi*time/Dt)
    return m_t
###################
data_stack=[]

file_list='rf_stack_5G.h5'
file_6k='rf_stack_6K_8st.h5'
# file_list='rf_stack_5G_rvr.h5'
file_list_2='rf_stack_5G_rvr_Dtrmv.h5'
file_6k_2='rf_stack_6K_rvr_Dtrmv_8st.h5'

stack_all=RFStream()
stack_all_2=RFStream()

stack=read(file_list,'H5')+read(file_6k,'H5')
stack_2=read(file_list_2,'H5')+read(file_6k_2,'H5')

for j in range(len(stack)):
    stack_all.append(stack[j])

for j in range(len(stack_2)):
    stack_all_2.append(stack_2[j])

# print('Len of stack=',len(stack_all))
####
stack_select=RFStream()
stack_select_2=RFStream()

station_file=open("P_delay_5G_6K_thin_sedi.txt","r")
for line1 in station_file:
    # lat_st=round(float(line1.split()[0]),3)
    # long_st=round(float(line1.split()[1]),3)
    RF_tps=float(line1.split()[2])
    name= line1.split()[3]
    for i in range(len(stack_all)):
        if name==stack_all[i].stats.station:
            stack_all[i].stats.tpdelay=RF_tps
            stack_select.append(stack_all[i])
    for i in range(len(stack_all_2)):
        if name==stack_all_2[i].stats.station:
            stack_all_2[i].stats.tpdelay=RF_tps
            stack_select_2.append(stack_all_2[i])

station_file.close()
print('Len of stack=',len(stack_select))
print('Len of stack2=',len(stack_select_2))

# sys.exit()
######

stack_select.sort(keys=['tpdelay']) # sorts stream wrt to tpdelay
stack_select_2.sort(keys=['tpdelay']) # sorts stream wrt to tpdelay
# sys.exit()

###
dt = UTCDateTime("1970-01-01T00:00:00")
stack_select_try=stack_select.copy()# -7.5 to 30 sec
stack_select_try_2=stack_select_2.copy() # -7.5 to 26 sec

stack_select_try.trim(dt+5,dt+22.5) # gets 2.5 sec before and 15 sec after onset
stack_select_try_2.trim(dt+5,dt+22.5) # gets 2.5 sec before and 15 sec after onset


####color map
# cmap=cm.viridis # plasma viridis cividis
# cmap=mpl.colormaps['plasma']
#######
########
fig = plt.figure(figsize=(5.5, 10))
ax1 = fig.add_axes([0.0, 0.07, 0.45, 0.85]) #[left, bottom, width, height]
ax1.set_ylim(31, 0)
ax1.set_xlim(-2.5,15)
ax1.set_yticks(np.linspace(0,31,32))
ax1.tick_params(axis='y',left=False,pad=1)
# ax1.yaxis.labelpad = -12
ax1.grid(which='major', axis='x',color='LightGrey', linestyle='--',linewidth=.95,alpha=.95)
ax1.xaxis.set_minor_locator(MultipleLocator(2.5))
ax1.xaxis.set_major_locator(MultipleLocator(5))
# ax1.set_xlabel('Time (s)')
# ax1.set_xlabel('Distance (km)')
# ax1.annotate('Vertical ACF AEB09', xy=(.6, 1.05), xycoords='axes fraction')
# ax1.set_title('Profile A-A\'',fontsize=12)
stack_ax1=stack_select_try#[:10]
stack2_ax1=stack_select_try_2#[:10]

for i in range(len(stack_ax1)):
    auto=stack_ax1[i]
    auto2=stack2_ax1[i]

    time = np.arange(auto.stats.npts) * auto.stats.delta
    time=time-2.5 # shifts onset to 0 sec
    auto.data /= np.max(np.abs(auto.data))
    auto2.data /= np.max(np.abs(auto2.data))
    if np.max(auto.data) == 1:
        auto.data=auto.data*-1
    if np.max(auto2.data) == 1:
        auto2.data=auto2.data*-1
    # if stack_ax1[i].stats.station == 'AEB04':
    #     auto.data=auto.data*-1
    # dist = dist_cal(auto,cross_lat,cross_long)
    l1,=ax1.plot(time,auto.data + i+1,  lw=1.15, color='darkkhaki',label='Uncorrected RF')
    l2,=ax1.plot(time,auto2.data + i+1,  lw=1.15, color='navy',alpha=.6,label='Sediment corrected RF')

    # ax1.fill_between(time, i+1, auto.data + i+1, lw=0.55,
    #                   color=cm.viridis(auto.stats.tpdelay/1.32), where=(auto.data < 0),alpha=.8)
    # plt.setp(auto.stats.station,fontsize=3)
#
st_label=[' ']
reso_filter=[' ']
for i in range(len(stack2_ax1)):
    auto=stack2_ax1[i]
    reso_st=[]
    # reso_filter.append('{:.2f}, {:.2f}, {:.2f}'.format(auto.stats.tpdelay,auto.stats.Dt_mean,auto.stats.Dt_std))
    st_label.append(auto.stats.station)
plt.setp(ax1.get_yticklabels(),fontsize=7)
ax1.yaxis.set_major_formatter(ticker.FixedFormatter(st_label))

########################################## Second column sinuodal decay
data_stack=[]
file_list='rf_stack_5G.h5'
file_6k='rf_stack_6K_8st.h5'
# file_list='rf_stack_5G_rvr.h5'
file_list_2='rf_stack_5G_rvr_Dtrmv.h5'
file_6k_2='rf_stack_6K_rvr_Dtrmv_8st.h5'

stack_all=RFStream()
stack_all_2=RFStream()

stack=read(file_list,'H5')+read(file_6k,'H5')
stack_2=read(file_list_2,'H5')+read(file_6k_2,'H5')

for j in range(len(stack)):
    stack_all.append(stack[j])

for j in range(len(stack_2)):
    stack_all_2.append(stack_2[j])

# print('Len of stack=',len(stack_all))
####
stack_select=RFStream()
stack_select_2=RFStream()

station_file=open("P_delay_5G_6K_thin_sedi.txt","r")
for line1 in station_file:
    # lat_st=round(float(line1.split()[0]),3)
    # long_st=round(float(line1.split()[1]),3)
    RF_tps=float(line1.split()[2])
    name= line1.split()[3]
    for i in range(len(stack_all)):
        if name==stack_all[i].stats.station:
            stack_all[i].stats.tpdelay=RF_tps
            stack_select.append(stack_all[i])
    for i in range(len(stack_all_2)):
        if name==stack_all_2[i].stats.station:
            stack_all_2[i].stats.tpdelay=RF_tps
            stack_select_2.append(stack_all_2[i])

station_file.close()
print('Len of stack=',len(stack_select))
print('Len of stack2=',len(stack_select_2))

######
stack_select.sort(keys=['tpdelay']) # sorts stream wrt to tpdelay
stack_select_2.sort(keys=['tpdelay']) # sorts stream wrt to tpdelay
# sys.exit()

###
dt = UTCDateTime("1970-01-01T00:00:00")
stack_select_try=stack_select.copy()# -7.5 to 30 sec
stack_select_try_2=stack_select_2.copy() # -7.5 to 26 sec

stack_select_try.trim(dt+5,dt+22.5) # gets 2.5 sec before and 15 sec after onset
stack_select_try_2.trim(dt+5,dt+22.5)
ax2 = fig.add_axes([0.5, 0.07, 0.45, 0.85]) #[left, bottom, width, height]
ax2.set_ylim(31, 0)
ax2.set_xlim(-2.5,15)
ax2.set_yticks(np.linspace(0,31,32))
ax2.tick_params(axis='y',left=False,labelleft=False,pad=1)
# ax1.yaxis.labelpad = -12
ax2.grid(which='major', axis='x',color='LightGrey', linestyle='--',linewidth=.95,alpha=.95)
ax2.xaxis.set_minor_locator(MultipleLocator(2.5))
ax2.xaxis.set_major_locator(MultipleLocator(5))
# ax2.set_xlabel('Time (s)')
#
stack_ax2=stack_select_try#[10:20]
stack2_ax2=stack_select_try_2#[10:20]

for i in range(len(stack_ax1)):
    auto=stack_ax1[i]
    auto2=stack2_ax1[i]

    time = np.arange(auto.stats.npts) * auto.stats.delta
    time=time-2.5 # shifts onset to 0 sec
    auto.data /= np.max(np.abs(auto.data))
    auto2.data /= np.max(np.abs(auto2.data))
    if np.max(auto.data) == 1:
        auto.data=auto.data*-1
    if np.max(auto2.data) == 1:
        auto2.data=auto2.data*-1
    # if stack_ax1[i].stats.station == 'AEB04':
    #     auto.data=auto.data*-1
    # dist = dist_cal(auto,cross_lat,cross_long)
    ax2.plot(time,auto.data + i+1,  lw=1.15, color='darkkhaki')
    ### r0 Dt business
    if not np.isnan(auto2.stats.r0_mean):
        dt=auto2.stats.tpdelay
        reso_filter=get_autoC(auto2.stats.r0_mean,auto2.stats.Dt_mean,len(auto2.data),auto2.stats.sampling_rate)
        cos_dec=get_decay_cos(auto2.stats.r0_mean,auto2.stats.Dt_mean,time+2.5)
        l3,=ax2.plot(time+2.5+auto.stats.tpdelay,cos_dec + i+1,  lw=1.15, color='maroon',alpha=.6,label='Sinuodal decay')

    # ax1.fill_between(time, i+1, auto.data + i+1, lw=0.55,
    #                   color=cm.viridis(auto.stats.tpdelay/1.32), where=(auto.data < 0),alpha=.8)
    # plt.setp(auto.stats.station,fontsize=3)
#
st_label=[' ']
reso_filter=[' ']
for i in range(len(stack2_ax1)):
    auto=stack2_ax1[i]
    reso_st=[]
    # reso_filter.append('{:.2f}, {:.2f}, {:.2f}'.format(auto.stats.tpdelay,auto.stats.Dt_mean,auto.stats.Dt_std))
    reso_filter.append('{:.2f}, {:.2f}'.format(auto.stats.tpdelay,auto.stats.Dt_mean))
    st_label.append(auto.stats.station)
plt.setp(ax2.get_yticklabels(),fontsize=7)
# ax1.yaxis.set_major_formatter(ticker.FixedFormatter())
ax_1T = ax2.twinx()
ax_1T.set_ylim(31,0)
ax_1T.set_yticks(np.linspace(0,31,32))
ax_1T.tick_params(axis='y',left=False,labelleft=False,labelright=True,pad=2,length=.1)
plt.setp(ax_1T.get_yticklabels(),fontsize=6.5)
ax_1T.yaxis.set_major_formatter(ticker.FixedFormatter(reso_filter))
################
plt.legend([l1,l2,l3],['Uncorrected RF','Sediment corrected RF','Sinusoidal decay'], loc=[-1,1.02],ncol=3,fontsize=9,handletextpad=.5,borderaxespad=1.5,columnspacing=1)
fig.text(0.5, 0.02, 'Time (s)',fontsize=12, ha='center', va='center')
fig.text(-0.04, 0.925, '(a)',fontsize=12, ha='center', va='center')
fig.text(0.99, 0.925, '(b)',fontsize=12, ha='center', va='center')

plt.savefig('5g_stack_rvr_tgthr_2col.pdf',bbox_inches='tight', pad_inches=0.15)
plt.show()

