# takes pre-computed rf's and refine them. Also, applies resonance filter from Yu et al., 2015.
# finally saves the stacks as well.
#rf_data is copied from Rese./south_aus/freq_vari_rf/rf_profile_data
# use rf_data_5g_FDSN.py from that folder to get data again if need be
# this script was to look at station one by one, curating the RFs with resonable values of Dt.
# the auto_c plot was used to do that.
# Later, used this script to save RFs for dt>.7s (from AEB12 onawards to save rf int0
#cruated rf files for later use.)
import matplotlib.pyplot as plt
import obspy
from obspy import read, Stream, UTCDateTime,read_events
from obspy.core.event import Origin
from obspy.taup import TauPyModel
from obspy.geodetics import gps2dist_azimuth
from obspy.clients.fdsn import Client
from obspy.core.trace import Trace
import numpy as np
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

def _get_residual(x, y):
    """
    Get the residual between x and y and the power of the residual.
    """
    r = x - y
    sumsq = np.dot(r, r)

    return r, sumsq
# end func

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
def apply_reverberation_filter(cha_data):
    result_stream = []
    station_file=open("stacked_rf_tpsb/P_delay_5G_noNan.txt","r")
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
########
client = Client("IRIS")
client_5g=Client('http://auspass.edu.au:80',user='5g',password='___')
inventory = client_5g.get_stations(network="5G",station='*',level='response')

station_name=[]
for i in range(len(inventory.get_contents()['stations'])):
    station_name.append(inventory.get_contents()['stations'][i].split()[0].split('.')[1])

no_sedi_st=['AES01','AES13','AES15'] #station with no sedi cover
thick_sedi_st=['AEB12','AEB11','AEB08','AEB10','AEB07','AEB02','AEB03','AEB09'\
,'AEB05','AES07','AES17','AES18','AES19','AES01','AES13','AES15'] # station with thick sedi
# last three are no sedi stations

station_name_S=[]
for name in station_name:
    if name[2]=='S':
         station_name_S.append(name)

# for station_5G in station_name_S:
for station_5G in ['AES20']:

    print('Doing', station_5G)

    try:
        data = os.path.join('rf_data/{}'.format(station_5G), '')
        # catfile = data+'rf_profile_events.xml'
        # datafile = data + 'rf_profile_data.h5'
        # rffile = data + 'rf_profile_rfs_1Hz.h5' ### change here for diff freq
        rffile_pdelay = data + 'rf_profile_rfs_pdel.h5' ###
        rffile_ref = data + 'rf_profile_rfs_ref_rvr.h5'
        rffile_ref_cura = data + 'rf_profile_rfs_ref_rvr_rmv.h5' ###Rf's with Dt constrains

    #############################################################

        rvr_rf = read_rf(rffile_ref, 'H5')
        print('No. of RF after P corr/rvr filter =',len(rvr_rf))

        if station_5G in thick_sedi_st:
            print('RVR not applied')
            rvr_rf.write(rffile_ref_cura, 'H5')
        else:
            stream_rever_rmv_cur=RFStream()
            for tr in rvr_rf:
                # tr=stream3c.select(component='R')
                # if tr.stats.max_P < 1:#used as another quality check for thick sediment stations
                if .1 < tr.stats.t3_offset < 1 and tr.stats.max_P < 1:#and tr.stats.distance :# change here for different stations based on Dt in auto_c plots
                    stream_rever_rmv_cur.append(tr)
            print('length after removing large DT=', len(stream_rever_rmv_cur))

            stream_rever_rmv_cur.write(rffile_ref_cura, 'H5')

            kw = {'trim': (-2.5, 25),'fig_width':6, 'fillcolors': ('cadetblue', 'gainsboro'), 'trace_height': 0.1, 'show_vlines': 'True','scale':2.5}
            stream_rever_rmv_cur.select(component='R', station='{}'.format(station_5G)).sort(['distance']).plot_rf(**kw)
            plt.savefig('figs_rf_dist/{}_{}_.1_1_rvrD_rmvDt.pdf'.format(station_5G,'R'),bbox_inches='tight', pad_inches=0.2)
            # plt.show()
            plt.close()

    except:
        print('Data/rf not found for station',station_5G)
    print('#############################################################')
