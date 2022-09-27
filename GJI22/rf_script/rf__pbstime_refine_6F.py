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

def _get_residual(x, y):
    """
    Get the residual between x and y and the power of the residual.
    """
    r = x - y
    sumsq = np.dot(r, r)

    return r, sumsq
# end func

def calc_cc_conv(data_rf,data):
# rf_all=data_rf.copy()
    rf_all=data_rf
    sig_all=data.copy()
    rf_all.sort(keys=['onset'])
    sig_all.sort(keys=['onset'])
    sig_all.rotate('NE->RT')
    sig_all.trim2(-25, 75, 'onset')
    ###
    cc=[]
    fit=[]
    if len(rf_all) == len(sig_all):
    # print('True')
        for i in range(0,int(len(rf_all)/3)):
        # for i in range(0,2)):
            try:
                rf=rf_all[3*i:3*i+3]
                sig=sig_all[3*i:3*i+3]

                # rf,sig are stream containing same three traces for RF and ZNE
                obs_Z = sig.select(component='Z')[0].copy()
                obs_R = sig.select(component='R')[0].copy()
                obs_rfR = rf.select(component='R')[0].copy()
                sr = obs_Z.stats.sampling_rate

                # Filter using SNR bandpass
                obs_Z.detrend().taper(max_percentage=0.05, max_length=2.)
                obs_R.detrend().taper(max_percentage=0.05, max_length=2.)
                obs_Z.filter('bandpass', freqmin=0.1, freqmax=1., zerophase=True)
                obs_R.filter('bandpass', freqmin=0.1, freqmax=1., zerophase=True)
                obs_rfR.filter('bandpass', freqmin=0.1, freqmax=1., zerophase=True)

                                # Convolve Z with rfR to obtain predicted R
                pred_R = obs_R.copy()
                pred_R.stats.channel = 'PRR'
                st=pred_R.stats.starttime
                time_shift=pred_R.stats.onset-pred_R.stats.starttime # rel P onset
                ind1 = int(np.ceil(time_shift*sr))
                # ind1 = int((len(obs_rfR.data)/2.))
                ind2 = ind1+len(obs_Z.data) # [leadin:leadin + n]
                pred_R.data = np.convolve(obs_Z.data, obs_rfR.data, mode='full')[ind1:ind2]
                # pred_R_1.data = signal.convolve(obs_Z.data, obs_rfR.data, mode='full')[ind1:ind2]
                # trim all traces from 0 to 20. sec following P-wave (fftshift first)
                # obs_Z.data = np.fft.fftshift(obs_Z.data)#[ind1:ind1+int(sr*1.5*10.)]
                obs_Z.trim(st+25,st+35)
                # obs_R.data = np.fft.fftshift(obs_R.data)#[ind1:ind1+int(sr*1.5*10.)]
                obs_R.trim(st+25,st+35)
                # pred_R.data = np.fft.fftshift(pred_R.data)#[ind1:ind1+int(sr*1.5*10.)]
                pred_R.trim(st+25,st+35)

                # resid, sumsq_ip1 = _get_residual(obs_R.normalize().data, pred_R.normalize().data)
                # power = np.dot(obs_R.data, obs_R.data)
                # fit.append(100 * (1 - sumsq_ip1/power))
                # rf.select(component='R')[0].stats.fit= 100 * (1 - sumsq_ip1/power)
                # Get cross correlation coefficient between observed and predicted Q

                rf.select(component='R')[0].stats.cc_conv = np.corrcoef(obs_R.data*obs_R.data, pred_R.data*pred_R.data)[0][1]
                # rf.select(component='R')[0].stats.cc_conv = np.corrcoef(abs(obs_R.data), abs(pred_R.data))[0][1]
                cc.append(np.corrcoef(obs_R.data, pred_R.data)[0][1])

                # print(rf.stats.cc_conv)
                # test = Stream(traces=[obs_Z, obs_R, pred_R])
                # test = Stream(traces=[obs_Z.normalize(), obs_R.normalize(), pred_R.normalize()])
                # test.plot()
                ########
            except:
                print('  ')#do nothing
    return cc,rf_all

def filter_crosscorr_coeff(rf_stream, time_window=(-2, 25), threshold_cc=0.60, min_fraction=0.20, apply_moveout=False):
    """:param rf_stream: Stream of RF traces to filter, should be **for a single component of a single station**
    :param threshold_cc: Threshold cross-correlation coefficient, defaults to 0.70. (0.8 in Sippl)
    :param min_fraction: Minimum fraction of coefficients above threshold_cc, defaults to 0.15.
    :return: Filtered stream of RF traces
    """
    # Early exit if we don't have enough traces for similarity filtering to be meaningful.
    if len(rf_stream) < 3:
        return rf_stream
    # end if

    # Trim RFs to time range so that subsequent cross-correlation computations relate to the
    # relevant period around and after onset.
    data_cc = rf_stream.copy().trim2(*time_window, 'onset')
    if not data_cc:
        return data_cc
    # end if

    # Apply optional moveout
    if apply_moveout:
        data_cc.moveout()
    # end if
    # Gather all RFs into a single array for efficient computation of correlation coefficients
    # between all traces
    data_array = np.array([tr.data for tr in data_cc])
    # Compute cross-correlation coefficients. cc matrix will be symmetric.
    #### Each row of cc indicates the degree of correlation between each other trace. ###### cc is a square matrix..
    cc = np.corrcoef(data_array)
    # Determine mask of which traces meet the similarity filtering criteria
    fraction_above_threshold = np.sum(cc >= threshold_cc, axis=1)/len(data_cc)
    keep_trace_mask = (fraction_above_threshold >= min_fraction)
    # kept_data = rf.RFStream([tr for i, tr in enumerate(rf_stream) if keep_trace_mask[i]])
    kept_data=RFStream()
    kept_data = RFStream([tr for i, tr in enumerate(rf_stream) if keep_trace_mask[i]])
    # kept_data = [tr for i, tr in enumerate(rf_stream) if keep_trace_mask[i]]
    return kept_data

def max_p_onset(self):
    """
    Determine max amplitude of r-rf around P-onset.
    """
    st = self.stats
    self_copy = self.copy()
    self_copy.taper(max_percentage=0.05)
    rel_time = getattr(st, 'onset')

    winP_st = rel_time - st.starttime - 0.5    # P onset -1
    winP_end = rel_time - st.starttime + 2  # P onset +1
    t = np.arange(self.stats.npts) * 1. / st.sampling_rate
    max_abs=max(abs(self_copy.data[(t >= winP_st) * (t <= winP_end)])) # gets max amp around P (-1,1 sec)
    max_=max(self_copy.data[(t >= winP_st) * (t <= winP_end)])
    st.max_abs_P=max_abs
    st.max_P=max_
    return max_abs,max_

########
client = Client("IRIS")
client_6F=Client('http://auspass.edu.au:80')#,user='6F',password='apple12')
inventory = client_6F.get_stations(network='6F', station='*',level='response', \
minlatitude=-39,maxlatitude=-24.9,minlongitude=128,maxlongitude=141.5)

station_name=[]
for i in range(len(inventory.get_contents()['stations'])):
    station_name.append(inventory.get_contents()['stations'][i].split()[0].split('.')[1])

# sys.exit()
time_delay=[]
station=[]

for station_6F in station_name:
    print('Doing', station_6F)
    try:
        # station_6F='ASR12'
        data = os.path.join('rf_data/{}'.format(station_6F), '')
        catfile = data+'rf_profile_events.xml'
        datafile = data + 'rf_profile_data.h5'
        rffile = data + 'rf_profile_rfs.h5'
        rffile_ref = data + 'rf_profile_rfs_ref.h5'
        # station= data+'perm_st.txt'


    #############################################################

        data_rf = read_rf(rffile, 'H5')
        data= read_rf(datafile, 'H5')
        # sys.exit()
        print('No. of RF after SNR =',len(data_rf)/3)

        data_rf_cc = RFStream()
        # sys.exit()
        (cc_conv_all,data_rf_cc)=calc_cc_conv(data_rf,data)
        # sys.exit()
        stream_rf = RFStream()
        for stream3c in tqdm(IterMultipleComponents(data_rf_cc, 'onset', 3)):
            stream3c.taper(max_percentage=0.05)
            max_rf_r_abs=max(abs(stream3c.select(component='R')[0].data))
            max_rf_r=max(stream3c.select(component='R')[0].data)
            cc_conv=stream3c.select(component='R')[0].stats.cc_conv
            stream3c.select(component='R')[0].stats.max_rf_r_abs=max_rf_r_abs
            (max_abs_P,max_P)=max_p_onset(stream3c.select(component='R')[0])
            if max_abs_P==max_rf_r_abs and max_P >0:# and cc_conv>0.4: # if max of Direct P in R-rf is greater than rest of signal
                stream_rf.extend(stream3c)


        stream_rf.write(rffile_ref, 'H5')
        print('No. of RF after direct P corrc. =',len(stream_rf)/3)

        kw = {'trim': (-5, 20),'fig_width':6, 'fillcolors': ('navy', 'grey'), 'trace_height': 0.1, 'show_vlines': 'True','scale':2.5}
        stream_rf.select(component='R', station='{}'.format(station_6F)).sort(['back_azimuth']).plot_rf(**kw)
        # plt.show()
        plt.savefig('figs_rf/{}_{}_.1_1_directP.pdf'.format(station_6F,'R'),bbox_inches='tight', pad_inches=0.2)
        # plt.savefig('{}_{}_.1_1_directP_conv.5.pdf'.format(station_6F,'R'),bbox_inches='tight', pad_inches=0.2)
        plt.close()
        ##################

        stack=stream_rf.stack()
        stack_R=stack.select(component='R')[0].data
        stack_Z=stack.select(component='Z')[0].data

        time_P= (np.argmax(stack_R)-np.argmax(stack_Z))/stack[0].stats.sampling_rate

        print("direct-P Time delay for station", station_6F, '=', time_P,'sec')
        time_delay.append(time_P)
        station.append(station_6F)
    except:
        print('Data/rf not found for station',station_6F)

ff=open('P_delay_lat_long_Bilby.txt','w')

for i in range(len(station)):
    station_6F=station[i]
    lat=inventory.select(station=station_6F)[0][0].latitude
    long=inventory.select(station=station_6F)[0][0].longitude

    ff.write('{} {} {} {}\n'.format(lat,long,float(str(np.around(time_delay[i],2))),station_6F))
ff.close()

sys.exit()
stream_rf_R = RFStream()
for stream3c in tqdm(IterMultipleComponents(stream_rf, 'onset', 3)):
    stream_rf_R.extend(stream3c.select(component='R'))

stream_rf_R_cc = RFStream()
stream_rf_R_cc=filter_crosscorr_coeff(stream_rf_R,threshold_cc=0.6, min_fraction=0.1)

print('No. of final RF =',len(stream_rf_R_cc))

stream_rf_R_cc.write(rffile_R_cc, 'H5')

# stream = read_rf(rffile, 'H5')

kw = {'trim': (-5, 20),'fig_width':6,'fillcolors': ('black', 'gray'),'trace_height': 0.1,\
'stack_height':0.45 ,'show_vlines': 'True','scale':2.5}
stream_rf_R_cc.select(component='R', station='{}'.format(station_6F)).sort(['back_azimuth']).plot_rf(**kw)
# plt.show()
plt.savefig('{}_{}_.1_1_time_cc.6.1_convo.pdf'.format(station_6F,'R'),bbox_inches='tight', pad_inches=0.2)
# sys.exit()

########
