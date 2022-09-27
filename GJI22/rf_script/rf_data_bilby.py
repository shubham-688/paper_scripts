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
# import functions_decov as dc
from importlib import reload
# reload(dc)

###################################
# This script uses one station data for events and calculates RF.
# Only data with SNR > 2 is saved.
##################################
def get_start_end_time(station,year):
    if year==2018:
        for line in open("/Users/Shubham/Research/Lake_eyre_data/Lake_Eyre_s1_time.txt", "r"):
            if len(line.split()) > 3:
                if line.split()[1]==station:
                    start_t=UTCDateTime(line.split()[6])
                    end_t=UTCDateTime(line.split()[7])
                    break
    if year==2019:
        for line in open("/Users/Shubham/Research/Lake_eyre_data/Lake_Eyre_s2_time.txt", "r"):
            if len(line.split()) > 3:
                if line.split()[0]==station:
                    start_t=UTCDateTime('{} {}'.format(line.split()[3],line.split()[4]))
                    end_t=UTCDateTime('{} {}'.format(line.split()[6],line.split()[7]))
                    break
    return start_t,end_t
################
def signoise(self):#, winsig, winnoise, relative='onset'):
        """
        Determine signal noise ratio by dividing the maximum in the two windows.
        """
        st = self.stats
        self_copy = self.copy()
        self_copy.detrend().taper(max_percentage=0.05)
        self_copy.filter('bandpass', freqmin=0.1, freqmax=1)#,corners=2, zerophase=True)
        winsig=[-5,25] #signal window
        winnoise=[-45,-15] #noise window
        # if relative in ('onset', 'middle'):
            # if relative == 'onset':
        rel_time = getattr(st, 'onset')
            # else:
            #     rel_time = st.starttime + (st.endtime - st.starttime) / 2
        winsig0 = rel_time - st.starttime + winsig[0]
        winsig1 = rel_time - st.starttime + winsig[1]
        winnoise0 = rel_time - st.starttime + winnoise[0]
        winnoise1 = rel_time - st.starttime + winnoise[1]
        #
        t = np.arange(self.stats.npts) * 1. / st.sampling_rate
        datasig = self_copy.data[(t >= winsig0) * (t <= winsig1)]
        datanoise = self_copy.data[(t >= winnoise0) * (t <= winnoise1)]
        # ipshell()
        try:
            st.signoise = max(abs(datasig)) / max(abs(datanoise))
            return st.signoise
        except:
            st.signoise = 0
            return st.signoise
########
# ddir= 'data'
client = Client("IRIS")
client_aus=Client('http://auspass.edu.au:80')#,user='6F',password='apple12')

inventory = client_aus.get_stations(network='6F', station='*',level='response', \
minlatitude=-39,maxlatitude=-24.9,minlongitude=128,maxlongitude=141.5)
# inventory.plot(projection='ortho')

# inventory = client_aus.get_stations(network="6F",station='*',level='response')

station_name=[]
for i in range(len(inventory.get_contents()['stations'])):
    station_name.append(inventory.get_contents()['stations'][i].split()[0].split('.')[1])

print(station_name)
# sys.exit()
for station_6F in station_name:
    print('Doing', station_6F)

    # station_6F='ASR13'
    data = os.path.join('rf_data/{}'.format(station_6F), '')
    catfile = data+'rf_profile_events.xml'
    datafile = data + 'rf_profile_data.h5'
    rffile = data + 'rf_profile_rfs.h5'
    # station= data+'perm_st.txt'

    if not os.path.exists(data):
        os.makedirs(data)

    # client_aus = Client('http://auspass.edu.au:8080',user_agent='shubham.agrawal@anu.edu.au') #,user='username',password='password')
    #st = client_aus.get_waveforms('5G','AES09','','*',starttime=obspy.UTCDateTime(2019,9,10,5,46,4),endtime=obspy.UTCDateTime(2019,9,10,5,46,24))
    st_lat=inventory.get_coordinates('6F.{}..BHZ'.format(station_6F))['latitude']
    st_long=inventory.get_coordinates('6F.{}..BHZ'.format(station_6F))['longitude']
    station_coord=inventory.get_coordinates('6F.{}..BHZ'.format(station_6F))


    # (starttime,endtime)=get_start_end_time('{}'.format(station_6F),2018)

    # starttime=obspy.UTCDateTime(2017,4,25,0,0,0)#2017-04-25
    # endtime=obspy.UTCDateTime(2019,2,28,0,0,0)#2019-02-28

    starttime=inventory.select(station=station_6F)[0][0].start_date
    endtime=inventory.select(station=station_6F)[0][0].end_date
    if int((endtime-starttime)/(24*3600)) > 365:
        endtime=starttime+365*24*3600
    print('Station',station_6F,'active for', int((endtime-starttime)/(24*3600)),'days')
    #
    if not os.path.exists(catfile):
        catalog = client.get_events(starttime=starttime, endtime=endtime,minmagnitude=5.5, \
        latitude=st_lat,longitude=st_long, minradius=30, maxradius=95)
        catalog.write(catfile, 'QUAKEML')

    catalog = read_events(catfile)
    # fig = catalog.plot(label=None)
    print('# of events:',len(catalog))
    # inv_station = client_aus.get_stations(network="6F",station=station_6F,level='response')

    if not os.path.exists(datafile):
        stream = RFStream()
        # kw = {'tt_model': 'LA_moho35'}
        with tqdm() as pbar:
            # try:
            for s in iter_event_data(catalog, inventory.select(station=station_6F), client_aus.get_waveforms, pbar=pbar):#,**kw):
                # s.filter('bandpass', freqmin=0.03, freqmax=5)
                SNR=signoise(s.select(component='Z')[0])
                if SNR >= 1.5:
                    s.detrend('linear')
                    s.detrend('demean').taper(max_percentage=0.05)
                    stream.extend(s)
            # except:
                # print('')
        stream.write(datafile, 'H5')
        print('Len of data per compo after SNR:',len(stream)/3)

    # sys.exit()
    try:
        data = read_rf(datafile, 'H5')
        stream = RFStream()
        # sys.exit()
        for stream3c in tqdm(IterMultipleComponents(data, 'onset', 3)):
            stream3c.detrend('linear')
            stream3c.detrend('demean').taper(max_percentage=0.05)
            stream3c.filter('bandpass', freqmin=0.1, freqmax=1)
            stream3c.trim2(-25, 75, 'onset')
            if len(stream3c) != 3:
                continue
            # kw = {'waterlevel': 0.05, 'gauss':2.5}
            #stream3c.rf(rotate='NE->RT',deconvolve='freq',**kw) # freq domain deconvolve
            stream3c.rf(rotate='NE->RT')
            # stream3c.rotate('NE->RT')

        # try:
        #     a.rf(rotate='NE->RT',deconvolve='func', func=dc.rf_iter_deconv)
        # except:
        #     print('Rf deconvolution didn`t work')
            stream3c.moveout()
            stream.extend(stream3c)


        stream.write(rffile, 'H5')

        print('No. of RF=',len(stream)/3,'for station',station_6F)
        print('--------------------------------------------------\n')
    except:
        print('NOt enough data/rf for station',station_6F)
        print('--------------------------------------------------\n')
sys.exit()
# stream = read_rf(rffile, 'H5')
kw = {'trim': (-5, 20),'fig_width':6, 'fillcolors': ('black', 'gray'), 'trace_height': 0.1, 'show_vlines': 'True','scale':2.5}
# kw = {'trim': (-5, 20), 'fillcolors': ('black', 'gray'), 'trace_height': 0.05, 'show_vlines': 'True','scale':2}
stream.select(component='R', station='{}'.format(station_6F)).sort(['back_azimuth']).plot_rf(**kw)
plt.savefig('{}_{}_.1_1.pdf'.format(station_6F,'R'),bbox_inches='tight', pad_inches=0.2)
#
