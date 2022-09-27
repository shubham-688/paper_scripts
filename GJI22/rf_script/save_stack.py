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

#############################################################
#############################################################
client = Client("IRIS")
client_6F=Client('http://auspass.edu.au:80')#,user='6F',password='apple12')
inventory = client_6F.get_stations(network='6F', station='*',level='response', \
minlatitude=-39,maxlatitude=-24.9,minlongitude=128,maxlongitude=141.5)

station_name=[]
for i in range(len(inventory.get_contents()['stations'])):
    station_name.append(inventory.get_contents()['stations'][i].split()[0].split('.')[1])

stack_6F=RFStream()
rf_stack_6F= 'rf_stack_6F.h5'
for station_6F in station_name:
    print('Doing', station_6F)
    try:
        # station_6F='ASR12'
        data = os.path.join('rf_data/{}'.format(station_6F), '')
        # catfile = data+'rf_profile_events.xml'
        # datafile = data + 'rf_profile_data.h5'
        # rffile = data + 'rf_profile_rfs.h5'
        rffile_ref = data + 'rf_profile_rfs_ref.h5' # final rfs used for stack

        # station= data+'perm_st.txt'


    #############################################################

        data_rf = read_rf(rffile_ref, 'H5')
        # data= read_rf(datafile, 'H5')
        # sys.exit()
        print('No. of RF after SNR =',len(data_rf)/3)
        stack=data_rf.stack()
        stack_6F.append(stack.select(component='R')[0])

    except:
        print('Data/rf not found for station',station_6F)

print('Len of all stack=', len(stack_6F))
stack_6F.write(rf_stack_6F,'h5')
