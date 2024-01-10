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
from obspy import read, Stream, Trace
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

def _get_residual(x, y):
    """
    Get the residual between x and y and the power of the residual.
    """
    r = x - y
    sumsq = np.dot(r, r)

    return r, sumsq
# end func
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
#sys.exit()
########
client = Client("IRIS")
client_5g=Client('http://auspass.edu.au:80',user='5g',password='apple12')
inventory = client_5g.get_stations(network="5G",station='*',level='response')

station_name=[]
for i in range(len(inventory.get_contents()['stations'])):
    station_name.append(inventory.get_contents()['stations'][i].split()[0].split('.')[1])

no_sedi_st=['AES01','AES13','AES15'] #station with no sedi cover
thick_sedi_st=['AEB12','AEB11','AEB08','AEB10','AEB07','AEB02','AEB03','AEB09'\
,'AEB05','AES07','AES17','AES18','AES19'] # station with thick sedi

station_name_S=[]
for name in station_name:
    if name[2]=='S':
         station_name_S.append(name)

# station_name=['AEB01']
rf_pcorr_stack_all=read_rf('stacked_rfs/rf_stack_5G.h5','H5')
rf_rvr_rmvDt_stack_all=read_rf('stacked_rfs/rf_stack_5G_rvr_Dtrmv.h5')

station_name_S=['AES04']
for station_5G in station_name_S:
    print('Doing', station_5G)


    try:
        data = os.path.join('rf_data/{}'.format(station_5G), '')
        # catfile = data+'rf_profile_events.xml'
        # datafile = data + 'rf_profile_data.h5'
        # rffile = data + 'rf_profile_rfs_1Hz.h5' ### change here for diff freq
        rffile_pdelay = data + 'rf_profile_rfs_pdel.h5' ###
        rffile_ref = data + 'rf_profile_rfs_ref_rvr.h5' ###
        rffile_ref_cura = data + 'rf_profile_rfs_ref_rvr_rmv.h5' ###Rf's with Dt constrains

    #############################################################
        rf_pcorr_stack=rf_pcorr_stack_all.select(station='{}'.format(station_5G))
        rf_rvr_rmvDt_stack=rf_rvr_rmvDt_stack_all.select(station='{}'.format(station_5G))
        r0=rf_rvr_rmvDt_stack[0].stats.r0_mean
        Dt_stack=rf_rvr_rmvDt_stack[0].stats.Dt_mean
        ###################

        rfs = read_rf(rffile_ref_cura, 'H5')
        print('No. of RFs =',len(rfs))
        rfs.sort(['back_azimuth'])
        rfs.trim2(-10,35,'onset')
        rfs_stack=bazi_stacks(rfs,10)
        rf_stack_all=rfs.stack()
        rf_stack_all[0].stats.totalrf=len(rfs)
        rf_stack_all.trim2(-2.5,15.5,'onset')
        rf_stack_all[0].data /= np.max(np.abs(rf_stack_all[0].data))
        rf_pcorr_stack.trim2(-2.5,15.5,'onset')
        rf_pcorr_stack[0].data /= np.max(np.abs(rf_pcorr_stack[0].data))



        print('Length of rf after bazi stack =',len(rfs_stack))
        dt=0
        Dt=0
        for line in open('st_dt_DT_rvr.txt'):
            line=line.split()
            if line[0]==station_5G:
                dt=float(line[1])
                Dt=float(line[2])
        for line in open('st_dt_rvr_thick.txt'):# loop for thick sedi stations.
            line=line.split()
            if line[0]==station_5G:
                dt=float(line[1])
                Dt=2*dt

        print('dt,Dt,Dt from stack, r0=',dt,Dt,Dt_stack,r0)
        # stream_rever_rmv=RFStream()
        # auto_c=RFStream()

        # kw = {'trim': (-2.5, 27.5),'fig_width':7, 'fillcolors': ('darksalmon', 'gainsboro'),\
        #  'trace_height': 0.1, 'show_vlines': 'True','scale':3,'info':(('back_azimuth', u'baz (°)', 'C0'),('distance', '', 'w'))}
        kw = {'trim': (-2.5, 15.5),'fig_width':5, 'fillcolors': ('darkseagreen', 'gainsboro'),\
        'trace_height': 0.1, 'show_vlines': 'True','scale':4.5,'info':None}
        mm.plot_rf_edi(rfs_stack,dt,Dt_stack,r0,rf_stack_all,rf_pcorr_stack,**kw)
        # plt.savefig('figs_baz_st/{}_.1_1_rvrADJ_10.pdf'.format(station_5G),bbox_inches='tight', pad_inches=0.1)
        plt.savefig('figs_baz/{}_{}_rvrADJ_10.pdf'.format(station_5G,'R'),bbox_inches='tight', pad_inches=0.2)
        # plt.close()
        plt.show()


    except:
        print('Data/rf not found for station',station_5G)
    print('###################################################')
