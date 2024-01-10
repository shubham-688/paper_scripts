# takes the saved rf's (filtered/rever removed/correlated) and creates a stacked RF .h5 file.
# also adds to the staked rf info about r0, Dt (resonance filter parameteres)
import matplotlib.pyplot as plt
import obspy
from obspy import read, Stream, UTCDateTime,read_events
from obspy.core.event import Origin
from obspy.taup import TauPyModel
from obspy.geodetics import gps2dist_azimuth
from obspy.clients.fdsn import Client
from obspy.core.trace import Trace
import numpy as np
import datetime
import os
from scipy.signal import hilbert, correlate
import rf
import sys
from rf import read_rf, rfstats,RFStream
from rf import iter_event_data, IterMultipleComponents
from tqdm import tqdm
from scipy import signal
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from matplotlib.font_manager import FontProperties

#########################################
# This script mkaes RF, auto, sinouidal, corr Rf plot for a station_name

#########################################
def stack(self):
    """
    Return stack of traces with same seed ids.

    Traces with same id need to have the same number of datapoints.
    Each trace in the returned stream will correspond to one unique seed
    id.
    """
    ids = set(tr.id for tr in self)
    tr = self[0]
    traces = []
    for id in ids:
        net, sta, loc, cha = id.split('.')
        data = np.mean([tr.data for tr in self if tr.id == id], axis=0)
        header = {'network': net, 'station': sta, 'location': loc,
                  'channel': cha, 'sampling_rate': tr.stats.sampling_rate}
        for entry in ('phase', 'moveout', 'station_latitude',
                      'station_longitude', 'station_elevation',
                      'processing'):
            if entry in tr.stats:
                header[entry] = tr.stats[entry]
        tr2 = RFTrace(data=data, header=header)
        if 'onset' in tr.stats:
            onset = tr.stats.onset - tr.stats.starttime
            tr2.stats.onset = tr2.stats.starttime + onset
        traces.append(tr2)
    return self.__class__(traces)

####### this gets mean of an array after removing values outside 2.5 std
def get_stdrem_mean(arr):
    # mags=[i[0] for i in mags_st]
    # if len(mags)==0:
    #     return 'nan'
    # elif len(mags)==1:
    #     mean=np.mean(np.array(mags))
    #     return round(float(mean),3)
    # else:
    # mags_ar=np.array(mags)
    mean = np.mean(arr, axis=0)
    sd = np.std(arr, axis=0)

    arr_stdrem = [x for x in arr if (x > mean - 2.5 * sd)]
    arr_stdrem = [x for x in arr_stdrem if (x < mean + 2.5 * sd)]
    mm=np.mean(np.array(arr_stdrem))
    std=np.std(np.array(arr_stdrem))

    # return round(float(mm),3)
    return mm,std

def get_decay_cos(r0,Dt,time):
    m_t=-r0**(time/Dt)*np.cos(np.pi*time/Dt)
    return m_t
###
def apply_reverberation_filter(cha_data):
    result_stream = []
    station_file=open("../P_delay_5G_noNan.txt","r")
    for line1 in station_file:
        RF_tps=float(line1.split()[2])
        name= line1.split()[3]
        if name==cha_data[0].stats.station:
            cha_data[0].stats.update({'pdelay':RF_tps})
    station_file.close()
    auto_c=RFStream()
    for i, tr in enumerate(cha_data):
        lead_time = tr.stats.onset - tr.stats.starttime
        auto_c_temp=Trace()
        relative_time = tr.times() - lead_time
        mask = np.array((relative_time < 0) | (relative_time > 5.0))
        loc = np.argmax(np.ma.masked_array(tr.data, mask=mask))
        dt = relative_time[loc]

        data = correlate(tr.data, tr.data, mode='full')
        data /= np.max(data)
        data = data[len(data) // 2:]
        auto_c_temp.data=data

        r0 = -(np.min(data))
        Dt = np.argmin(data) * 1. / tr.stats.sampling_rate

        tr_copy = tr.copy()
        resonance_filter = np.zeros(len(tr.data))
        resonance_filter[0] = 1
        resonance_filter[int(Dt * tr.stats.sampling_rate)] = r0
        tr_copy.data = np.convolve(tr_copy.data, resonance_filter, mode='full')
        tr_copy.data = tr_copy.data[:len(tr_copy.data) // 2 + 1]

        assert tr.data.shape == tr_copy.data.shape, 'Input/output length mismatch detected in ' \
                                                    'reverberation removal routine'

        tr_copy.stats.update({'t1_offset':dt,
                              't2_offset':Dt - dt,
                              't3_offset':Dt,
                              'r0':r0})
        result_stream.append(tr_copy)
        tr_ccc=tr_copy.copy()
        auto_c_temp.stats=tr_ccc.stats
        auto_c_temp.stats.starttime=tr_ccc.stats.onset # setting startime as onset for autoC
        t=auto_c_temp.stats.starttime
        auto_c.append(auto_c_temp.trim(t,t+30))
    # end for

    return rf.RFStream(result_stream),auto_c
# end func
#############################################################
#############################################################
data = os.path.join('../rf_data/AEB14', '')
rffile_pdelay = data + 'rf_profile_rfs_pdel.h5'
rf_pd=read_rf(rffile_pdelay, 'H5')

stream_rever_rmv=RFStream()
auto_c=RFStream()
stream_rever_rmv,auto_c=apply_reverberation_filter(rf_pd.select(component='R'))
auto_c_stack=auto_c.stack()

rf_stack_5G='rf_stack_5G.h5'
# rf_stack_5G_rvr= 'rf_stack_5G_rvr.h5'
rf_stack_5G_rvr_Dtrmv= 'rf_stack_5G_rvr_Dtrmv.h5'

data_rf = read_rf(rf_stack_5G, 'H5')
# data_rf = read_rf(rf_stack_5G_rvr_Dtrmv, 'H5')

data_rf_rvr_dtrm = read_rf(rf_stack_5G_rvr_Dtrmv, 'H5')

data_rf.select(station='AEB14')[0].stats.back_azimuth=100
data_rf.select(station='AEB14')[0].stats.distance=10

# kw = {'trim': (-2.5, 10),'fig_width':7, 'fillcolors': ('steelblue', 'gainsboro'), \
# 'trace_height': .5, 'show_vlines': 'True','scale':1.25}#,'info':None}
# data_rf.select(station='AEB14').plot_rf(**kw)
# plt.savefig('AEB14_rvr.pdf',bbox_inches='tight', pad_inches=0.2)
# plt.show()
#SeaGreen
# sys.exit()
###################################

fig = plt.figure(figsize=(7, 5))
ax1 = fig.add_axes([0.1, 0.7, 0.8, 0.2]) #[left, bottom, width, height]
ax1.set_xlim(-2.5,10)
ax1.tick_params(axis='y',left=False,labelleft=False,pad=1)
ax1.grid(which='major', axis='x',color='LightGrey', linestyle='--',linewidth=.75,alpha=.85)
# ax1.grid(which='major', axis='y',color='LightGrey', linestyle='--',linewidth=.75,alpha=.85)

ax1.xaxis.set_minor_locator(MultipleLocator(1.25))
ax1.xaxis.set_major_locator(MultipleLocator(2.5))
ax1.tick_params(axis='x',bottom=False,labelbottom=False,pad=1)

auto=data_rf.select(station='AEB14')[0]
time = np.arange(auto.stats.npts) * auto.stats.delta
time=time-7.5 # shifts onset to 0 sec
auto.data /= np.max(np.abs(auto.data))
ax1.plot(time,auto.data*1,  lw=1.25, color='darkSeaGreen')
ax1.fill_between(time, 0, auto.data*1, lw=0.05,color='SeaGreen', where=(auto.data > 0),alpha=.75)

#######################
ax2 = fig.add_axes([0.1, 0.5, 0.8, 0.2]) #[left, bottom, width, height]
ax2.set_xlim(-2.5,10)
ax2.tick_params(axis='y',left=False,labelleft=False,pad=1)
ax2.tick_params(axis='x',bottom=False,labelbottom=False,pad=1)
ax2.xaxis.set_minor_locator(MultipleLocator(1.25))
ax2.xaxis.set_major_locator(MultipleLocator(2.5))
ax2.grid(which='major', axis='x',color='LightGrey', linestyle='--',linewidth=.75,alpha=.85)
# ax2.grid(which='major', axis='y',color='LightGrey', linestyle='--',linewidth=.75,alpha=.85)

auto=auto_c_stack[0]
time_ = np.arange(auto.stats.npts) * auto.stats.delta
ax2.plot(time_,auto,  lw=1.25, color='slategrey')
ax2.fill_between(time_, 0, auto, lw=0.05,color='slategrey', where=(auto.data > 0),alpha=.5)

#######################
ax3 = fig.add_axes([0.1, 0.3, 0.8, 0.2]) #[left, bottom, width, height]
ax3.set_xlim(-2.5,10)
ax3.tick_params(axis='y',left=False,labelleft=False,pad=1)
ax3.tick_params(axis='x',bottom=False,labelbottom=False,pad=1)
ax3.xaxis.set_minor_locator(MultipleLocator(1.25))
ax3.xaxis.set_major_locator(MultipleLocator(2.5))
ax3.grid(which='major', axis='x',color='LightGrey', linestyle='--',linewidth=.75,alpha=.85)
# ax3.grid(which='major', axis='y',color='LightGrey', linestyle='--',linewidth=.75,alpha=.85)

auto=data_rf_rvr_dtrm.select(station='AEB14')[0]
r0=auto.stats.r0_mean
Dt=auto.stats.Dt_mean
cos_dec=get_decay_cos(r0,Dt,time+7.5)

ax3.plot(time+7.5,cos_dec*-1,  lw=1.25, color='darkkhaki')
ax3.fill_between(time+7.5, 0, cos_dec*-1, lw=0.05,color='khaki', where=(cos_dec < 0),alpha=.8)
#plotting r0 Dt
ax2.plot(Dt,-r0,marker='o', color='maroon',markersize=4,alpha=.8)
ax2.plot((Dt,Dt),(0,-r0),color='maroon',ls='--',lw=1.5,alpha=.65)
ax2.plot((0,Dt),(-r0,-r0),color='maroon',ls='--',lw=1.5,alpha=.65)
# fig.text(0.31, 0.53, 'r$_{0}$',fontsize=10, color='maroon',ha='center', va='center')
fig.text(0.305, 0.48, '($\Delta$t, -r$_{0}$)',fontsize=10, color='maroon',ha='center', va='center')
####################### this is last row
ax4 = fig.add_axes([0.1, 0.1, 0.8, 0.2]) #[left, bottom, width, height]
ax4.set_xlim(-2.5,10)
ax4.xaxis.set_minor_locator(MultipleLocator(1.25))
ax4.xaxis.set_major_locator(MultipleLocator(2.5))
ax4.tick_params(axis='y',left=False,labelleft=False,pad=1)
ax4.tick_params(axis='x',bottom=True,labelbottom=True,pad=1)
ax4.grid(which='major', axis='x',color='LightGrey', linestyle='--',linewidth=.75,alpha=.85)
# ax4.grid(axis='y',color='LightGrey', linestyle='--',linewidth=.75,alpha=.85)

auto=data_rf_rvr_dtrm.select(station='AEB14')[0]
time = np.arange(auto.stats.npts) * auto.stats.delta
time=time-7.5 # shifts onset to 0 sec
auto.data /= np.max(np.abs(auto.data))
ax4.plot(time,auto.data*1,  lw=1.25, color='steelblue')
ax4.fill_between(time, 0, auto.data*1, lw=0.05,color='steelblue', where=(auto.data > 0),alpha=.65)
ax4.set_xlabel('Time (s)',fontsize=12,fontname='serif')
fig.text(0.1, .92, '5G.AEB14',fontsize=12,fontname='serif',fontstyle='italic', color='navy',ha='left', va='center')
fig.text(0.89, 0.85, 'Uncorrected RF',fontsize=10,fontname='serif',color='darkGreen',ha='right', va='center')
fig.text(0.89, 0.65, 'RF Autocorrelation',fontsize=10,fontname='serif',color='dimGrey',ha='right', va='center')
fig.text(0.89, 0.45, 'Sinusoidal decay',fontsize=10,fontname='serif',color='olive',ha='right', va='center')
fig.text(0.89, 0.25, 'Sediment corrected RF',fontsize=10,fontname='serif',color='steelblue',ha='right', va='center')

fig.text(0.295, 0.28, 'PS$_{b}$',fontsize=11,fontname='serif',color='maroon',ha='left', va='center')
fig.text(0.59, 0.235, 'PS$_{m}$',fontsize=11,fontname='serif',color='maroon',ha='left', va='center')

plt.savefig('aeb14.pdf',bbox_inches='tight', pad_inches=0.15)

plt.show()

#
