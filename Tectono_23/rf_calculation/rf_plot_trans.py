# takes pre-computed rf's and refine them. Also, applies resonance filter from Yu et al., 2015.
# finally saves the stacks as well.
#rf_data is copied from Rese./south_aus/freq_vari_rf/rf_profile_data
# use rf_data_5g_FDSN.py from that folder to get data again if need be
import matplotlib.pyplot as plt
import obspy
from obspy import read, Stream, UTCDateTime,read_events
from obspy.core.event import Origin
from obspy.taup import TauPyModel
from obspy.geodetics import gps2dist_azimuth
from obspy.clients.fdsn import Client
from obspy.core.trace import Trace
import numpy as np
import scipy.fft as fft
from scipy import signal
from scipy.signal import hilbert, correlate
import datetime
import os
import sys
import rf
from rf import read_rf, rfstats,RFStream
from rf import iter_event_data, IterMultipleComponents
from tqdm import tqdm
from scipy import signal
from matplotlib import ticker
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from math import cos, sin, asin, sqrt, radians
from matplotlib import cm
import rf_plotrf as mm
from importlib import reload
reload(mm)
#########################################
# This script refines the saved RF.
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
##
def max_p_onset(self):
    """
    Determine max amplitude of r-rf around P-onset.
    """

    onset=self.stats.onset
    self_copy = self.copy()
    self_copy.trim(onset-5,onset+ 30)
    self_copy.taper(max_percentage=0.05)
    st = self_copy.stats
    rel_time = getattr(st, 'onset')

    winP_st = rel_time - st.starttime - 0.1    # P onset -1
    winP_end = rel_time - st.starttime + 1.5  # P onset +1
    t = np.arange(st.npts) * 1. / st.sampling_rate
    max_abs=max(abs(self_copy.data[(t >= winP_st) * (t <= winP_end)])) # gets max amp around P (-1,1 sec)
    max_=max(self_copy.data[(t >= winP_st) * (t <= winP_end)])
    self.stats.max_abs_P=max_abs
    self.stats.max_P=max_
    return max_abs,max_
#
# function to stack rfs acoording to bazi
def bazi_stacks(rfs,bin):
    rf_bazi=RFStream()
    rf_=rfs.copy()
    baz_bin=np.arange(0,360,bin)
    for i,b in enumerate(baz_bin):
        bin_rfs=RFStream()
        for j,tr in enumerate(rf_):
             if b < tr.stats.back_azimuth <= b+bin:
                 bin_rfs.append(tr)
        bin_rfs_stack=RFStream()
        if len(bin_rfs)>0:
            bin_rfs_stack=bin_rfs.stack()
            bin_rfs_stack[0].stats.update({'back_azimuth':b+bin/2,'distance':0,'rfs_in_bin':len(bin_rfs)})
            # bin_rfs_stack.stats.network=tr.stats.network
            # bin_rfs_stack.stats.station=tr.stats.station
            rf_bazi.append(bin_rfs_stack[0])

    return rf_bazi
#
########
client = Client("IRIS")
client_5g=Client('http://auspass.edu.au:80',user='5g',password='apple12')
inventory = client_5g.get_stations(network="5G",station='*',level='response')

station_name=[]
for i in range(len(inventory.get_contents()['stations'])):
    station_name.append(inventory.get_contents()['stations'][i].split()[0].split('.')[1])

time_delay=[]
station=[]
# station_name= ['AEB07','AEB13','AEB12','AES08','AES12','AEB15','AES16','AES03','AES06']
no_sedi_st=['AES01','AES13','AES15'] #station with no sedi cover

thick_sedi_st=['AEB12','AEB11','AEB08','AEB10','AEB07','AEB02','AEB03','AEB09'\
,'AEB05','AES07','AES17','AES18','AES19','AES01','AES13','AES15'] # station with thick sedi
# last three are no sedi stations

# thick_sedi_st=['AES07','AES17','AES18','AES19'] # station with thick sedi
station_name_S=[]
for name in station_name:
    if name[2]=='B':
         station_name_S.append(name)

station_name_S=['AES01','AEB15','AES04','AEB19','AEB11']
station_name_S=['AES13','AES15']

for station_5G in station_name_S:
    print('Doing', station_5G)

    # station_5G='ASR12'
    data = os.path.join('rf_data/{}'.format(station_5G), '')
    catfile = data+'rf_profile_events.xml'
    datafile = data + 'rf_profile_data.h5'
    rffile = data + 'rf_profile_rfs.h5' ### change here for diff freq
    rffile_pdelay = data + 'rf_profile_rfs_pdel.h5' ### change here for diff freq
    rffile_ref = data + 'rf_profile_rfs_ref_rvr.h5' ### chnage here for diff freq
#############################################################

    data_rf = read_rf(rffile_pdelay, 'H5')
    # data= read_rf(datafile, 'H5')
    print('No. of RF after direct P corrc. =',len(data_rf)/3)

    kw = {'trim': (-2.5, 25),'fig_width':6, 'fillcolors': ('palevioletred', 'darkgrey'), 'trace_height': 0.1, 'show_vlines': 'True','scale':2.5}
    # data_rf.select(component='T', station='{}'.format(station_5G)).sort(['distance']).plot_rf(**kw)
    # plt.savefig('figs_rf_t/{}_{}_directP.pdf'.format(station_5G,'T'),bbox_inches='tight', pad_inches=0.2)
    plt.close()
    #######
    rfs=data_rf.select(component='T').copy()
    rfs.sort(['back_azimuth'])
    rfs.trim2(-10,35,'onset')
    rfs_stack=bazi_stacks(rfs,10)
    rf_stack_all=rfs.stack()
    rf_stack_all[0].stats.totalrf=len(rfs)

    kw = {'trim': (-2.5, 15.5),'fig_width':5, 'fillcolors': ('palevioletred', 'gainsboro'),\
    'trace_height': 0.1, 'show_vlines': 'True','scale':4.5,'info':None}
    mm.plot_rf_t(rfs_stack,0,0,rf_stack_all,**kw)
    # plt.savefig('figs_baz_st/{}_.1_1_rvrADJ_10.pdf'.format(station_5G),bbox_inches='tight', pad_inches=0.1)
    plt.savefig('figs_rf_t/{}_{}_10.pdf'.format(station_5G,'T'),bbox_inches='tight', pad_inches=0.2)
