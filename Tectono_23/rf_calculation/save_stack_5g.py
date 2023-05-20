# takes the saved rf's (filtered/rever removed/correlated) and creates a stacked RF .h5 file.
# also adds to the staked rf info about r0, Dt (resonance filter parameteres)
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
from tqdm import tqdm
from scipy import signal
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

###
#############################################################
#############################################################
client = Client("IRIS")
client_5g=Client('http://auspass.edu.au:80',user='5g',password='apple12')
inventory = client_5g.get_stations(network="5G",station='*',level='response')

station_name=[]
for i in range(len(inventory.get_contents()['stations'])):
    station_name.append(inventory.get_contents()['stations'][i].split()[0].split('.')[1])

stack_5G=RFStream()
rf_stack_5G='rf_stack_5G.h5'
# rf_stack_5G= 'rf_stack_5G_rvr.h5'
# rf_stack_5G= 'rf_stack_5G_rvr_Dtrmv.h5'

###
data_rf_rvr=RFStream()
# station_name=['AEB04']
for station_5G in station_name:
    print('Doing', station_5G)
    try:
        data = os.path.join('rf_data/{}'.format(station_5G), '')
        rffile = data + 'rf_profile_rfs.h5'
        rffile_pdelay = data + 'rf_profile_rfs_pdel.h5' ###
        rffile_ref = data + 'rf_profile_rfs_ref_rvr.h5' # final rfs used for stack
        rffile_ref_cura = data + 'rf_profile_rfs_ref_rvr_rmv.h5' ###Rf's with Dt constrains

        #######
        data_rf_rvr = read_rf(rffile_pdelay, 'H5')
        # data_rf_rvr = read_rf(rffile_ref, 'H5')
        # data_rf_rvr = read_rf(rffile_ref_cura, 'H5')

        data_rf_rvr.sort(['back_azimuth'])
        data_rf_rvr.trim2(-7.5,30,'onset')
        print('No. of RF after directP =',len(data_rf_rvr)/3)
        ##########
	    # take mean of r0, dt
        # Dt=[]
        # r0=[]
        # try:
        #     Dt = np.array([tr.stats.t3_offset for tr in data_rf_rvr])
        #     r0 = np.array([tr.stats.r0 for tr in data_rf_rvr])
        # except:
        #     pass
        #######
        # mean,std=get_stdrem_mean(Dt)# used this to remove large std vals.

        # mean=np.mean(Dt)
        # std=np.std(Dt)
        # print('r0_mean=',np.mean(r0),'Dt_mean=',mean,'Dt_std=',std)
        #
        stack=data_rf_rvr.stack()
        # try:
        #     stack.select(component='R')[0].stats.update({'r0_mean':np.mean(r0),'Dt_mean':mean,'Dt_std':std})
        # except:
        #     pass
        stack_5G.append(stack.select(component='R')[0])

    except:
        print('Data/rf not found for station',station_5G)

print('Len of all stack=', len(stack_5G))
stack_5G.write(rf_stack_5G,'h5')
##

# used the bottom bit for plotting supplimenatry fig for BSSA

# kw = {'trim': (-5, 15),'fig_width':5.5, 'fillcolors': ('SeaGreen', 'DarkGray'), \
# 'trace_height': 0.1, 'show_vlines': 'True','scale':2.5}
# data_rf.select(component='R', station='{}'.format(station_5G)).sort(['back_azimuth']).plot_rf(**kw)
# plt.savefig('AEB07_rf.pdf',bbox_inches='tight', pad_inches=0.2)
# plt.show()
