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
from rf.imaging import plot_profile_map
from rf.profile import profile
from rf import iter_event_data, IterMultipleComponents,get_profile_boxes
from tqdm import tqdm
from scipy import signal
from matplotlib import ticker
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from math import cos, sin, asin, sqrt, radians
from matplotlib import cm
import fx_ccp as fx
from importlib import reload
reload(fx)
#########################################

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
###

########
client = Client("IRIS")
client_5g=Client('http://auspass.edu.au:80',user='5g',password='apple12')
inventory = client_5g.get_stations(network="5G",station='*',level='response')

station_name=[]
for i in range(len(inventory.get_contents()['stations'])):
    station_name.append(inventory.get_contents()['stations'][i].split()[0].split('.')[1])

prof_a=['AEB10','AEB06','AEB05','AEB03','AES17','AES18','AES19','AES07','AEB07','AES20']#,'AEB02'
prof_b=['AEB07','AEB08','AEB09','AEB11','AEB13','AEB12']
prof_c1=['AEB02','AES07','AES08','AES09','AES11','AES10','AES12','AEB14']#'AEB01'
prof_c2=['AES13','AEB15','AES15','AEB16','AEB17','AEB18','AEB19','AEB20','AES20']
prof_c3=['AES01','AES04','AES16','AES14','AEB04','AES06']#'AES02','AES03','AES05'

prof_c_big=prof_c1+prof_c2+prof_c3

# inventory = client_5g.get_stations(network="5G",station='AEB07,AEB08,AEB09,AEB11,AEB13,AEB12',level='response')
inventory = client_5g.get_stations(network="5G",station=','.join(prof_c_big),level='response')

rf_all=RFStream()
# station_name=['AES01']
for station_5G in prof_c_big:
    # print('Doing', station_5G)
    try:
        data = os.path.join('../rf_data/{}'.format(station_5G), '')
        # catfile = data+'rf_profile_events.xml'
        # datafile = data + 'rf_profile_data.h5'
        # rffile = data + 'rf_profile_rfs_1Hz.h5' ### change here for diff freq
        rffile_pdelay = data + 'rf_profile_rfs_pdel.h5' ###
        rffile_ref = data + 'rf_profile_rfs_ref_rvr.h5' ###
        rffile_ref_cura = data + 'rf_profile_rfs_ref_rvr_rmv.h5' ###Rf's with Dt constrains
    #############################################################

        rfs = read_rf(rffile_ref_cura, 'H5')
        rf_all.extend(rfs)

    except:
        print('Data/rf not found for station',station_5G)
##
print('{} 5G stations done'.format(len(prof_c_big)))
### aus array stations
# 6K stations with sedi: NILPI, WHYML, WILGE, WIRRE
# 6K stations with no results: COMWH, YALYM, MERNA
AusArray_st=['TWINS','MILCK','BILLA','STUCK','PARAK','MULGA','WITCH'\
'ARCOO','BOSWO','OOAKD','SOGAP','KOOTA']# ADD 'NILPI','MERNA','WIRRE','WILPO'

client_sa = Client('http://auspass.edu.au:8080',user='6k',password='melon93')
inventory_sa = client_sa.get_stations(network='6K', station=','.join(AusArray_st),level='response')

for station_6K in AusArray_st:
# for station_6K in ['WIRRE']:
    try:
        data = os.path.join('../../sa_array_rf/rf_data/{}'.format(station_6K), '')
        rffile_pdelay = data + 'rf_profile_rfs_pdel.h5' ###
        # #rffile_ref = data + 'rf_profile_rfs_ref_rvr.h5' ###
        rffile_ref_cura = data + 'rf_profile_rfs_ref_rvr_rmv.h5' ###Rf's with Dt constrains
    #############################################################
        rfs = read_rf(rffile_ref_cura, 'H5')
        if rfs[0].stats.sampling_rate ==200:
            rfs.decimate(2) # for station STUCK. It has sampling rate of 200.
        rf_all.extend(rfs.select(component='R'))
    except:
        print('Data/rf not found for station',station_5G)
##
print('Len of all rfs',len(AusArray_st+prof_c_big))
####
# rf_all.write('rfs_CC_big.H5', 'H5')
print('Len of rf_all=',len(rf_all)) # 6064!

#rf_all=read_rf('rfs_CC_big.H5', 'H5')
############ CCP bit
#kw = {'tt_model': 'LA_moho35'}
## PLOTTING C, 133.8 -25.5,end 138.5 -31.75 # 833km #147.67 azimuth
# rf_all=read_rf('rfs_C.H5', 'H5')
print('############')
print('Doing Profile C')
stream=rf_all.copy()
# stream=rf_all.select(component='R').copy()
# stream=rf_all.select(network='5G').copy()

ppoints = stream.ppoints(50)
#ppoints(pp_depth, pp_phase=None, model='iasp91')Return coordinates of piercing point calculated by 1D ray tracing.
### 850,171
boxes = get_profile_boxes((-25.5, 133.8), 147.67, np.linspace(0, 850,171), width=220)
#get_profile_boxes(latlon0, azimuth, bins, width=2000)

fig, ax=plt.subplots(figsize=(12,12))
fx.plot_profile_map(boxes,inventory=inventory_sa+inventory, ppoints=ppoints)
plt.savefig('pp_50km_CC_big_n.pdf',bbox_inches='tight', pad_inches=0.1)
plt.box(False)
plt.show()

sys.exit()
######
pstream = profile(tqdm(stream), boxes)
# pstream.write('profile_stream_CC_50km_big.H5', 'H5')
# pstream=read_rf('profile_stream_CC_50km_big.H5', 'H5')
print(pstream)
#####
pstream_temp=pstream.copy()
pstream_4=RFStream()
for tr in pstream:
    if tr.stats.num>4:
        pstream_4.append(tr)
pstream_temp=pstream_4
pstream_temp.trim2(-2.5, 15, 'onset')
# plt.figure(figsize=(30, 10))
fx.plot_profile(profile=pstream_temp.normalize(),scale=7,figsize=(19, 5.5), top='hist',fillcolors=('indianred', 'lightsteelblue'))#,moveout_model='ak135')
plt.gcf()#.set_size_inches(50, 10)
plt.savefig('ccp_50km_CC_>4rfs_n.png',bbox_inches='tight', pad_inches=0.1,dpi=400)
plt.show()
