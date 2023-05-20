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

                pred_R = obs_R.copy()
                pred_R.stats.channel = 'PRR'
                st=pred_R.stats.starttime
                time_shift=pred_R.stats.onset-pred_R.stats.starttime # rel P onset
                ind1 = int(np.ceil(time_shift*sr))

                ind2 = ind1+len(obs_Z.data) # [leadin:leadin + n]
                pred_R.data = np.convolve(obs_Z.data, obs_rfR.data, mode='full')[ind1:ind2]

                obs_Z.trim(st+20,st+45)
                obs_R.trim(st+20,st+45)
                pred_R.trim(st+20,st+45)
                rf.select(component='R')[0].stats.cc_conv = np.corrcoef(obs_R.data*obs_R.data, pred_R.data*pred_R.data)[0][1]
                # rf.select(component='R')[0].stats.cc_conv = np.corrcoef(abs(obs_R.data), abs(pred_R.data))[0][1]
                cc.append(np.corrcoef(obs_R.data, pred_R.data)[0][1])
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
def has_reverberations(cha_data, dt_max=0.2):
    """
    Checks if the data shows signs of the presence of reverberations
    :param cha_data: List or iterable of RF traces to use for H-k stacking.
    :type cha_data: Iterable(rf.RFTrace)
    :param dt_max: if the median temporal offset between an RF peak and the onset time > dt_max,
                   this function returns true
    :type dt_max: float
    return: Bool
    """

    dt_array = []
    for i, tr in enumerate(cha_data):
        lead_time = tr.stats.onset - tr.stats.starttime

        relative_time = tr.times() - lead_time
        mask = np.array((relative_time < 0) | (relative_time > 5.0))
        loc = np.argmax(np.ma.masked_array(tr.data, mask=mask))
        dt_array.append(relative_time[loc])
    # end for
    dt_array = np.array(dt_array)

    if(np.median(dt_array) > dt_max): return True

    return False
# end func
#
def apply_reverberation_filter(cha_data):
    result_stream = []
    station_file=open("P_delay_5G_noNan.txt","r")
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
client_5g=Client('http://auspass.edu.au:80',user='5g',password='____')
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

station_name_S=['AEB02']
for station_5G in station_name_S:
    print('Doing', station_5G)
    try:
        # station_5G='ASR12'
        data = os.path.join('rf_data/{}'.format(station_5G), '')
        catfile = data+'rf_profile_events.xml'
        datafile = data + 'rf_profile_data.h5'
        rffile = data + 'rf_profile_rfs.h5' ### change here for diff freq
        rffile_pdelay = data + 'rf_profile_rfs_pdel.h5' ### change here for diff freq
        rffile_ref = data + 'rf_profile_rfs_ref_rvr.h5' ### chnage here for diff freq
    #############################################################

        data_rf = read_rf(rffile, 'H5')
        # data= read_rf(datafile, 'H5')
        print('No. of RF after SNR =',len(data_rf)/3)

        # data_rf_cc = RFStream()
        # (cc_conv_all,data_rf_cc)=calc_cc_conv(data_rf,data)
        stream_rf = RFStream()
        for stream3c in tqdm(IterMultipleComponents(data_rf, 'onset', 3)):
            stream3c.taper(max_percentage=0.05)
            # stream3c.filter('bandpass', freqmin=0.1, freqmax=1)
            max_rf_r_abs=max(abs(stream3c.select(component='R')[0].data))
            max_rf_r=max(stream3c.select(component='R')[0].data)
            # cc_conv=stream3c.select(component='R')[0].stats.cc_conv
            stream3c.select(component='R')[0].stats.max_rf_r_abs=max_rf_r_abs
            (max_abs_P,max_P)=max_p_onset(stream3c.select(component='R')[0])
            # if max_abs_P==max_rf_r_abs:# and max_P >0:# and cc_conv>0.4: # if max of Direct P in R-rf is greater than rest of signal
            if max_P==max_rf_r and max_rf_r_abs <1:# and cc_conv>0.4 :# and max_P >0:# and cc_conv>0.4:
                stream_rf.extend(stream3c)

        stream_rf.write(rffile_pdelay, 'H5')
        print('No. of RF after direct P corrc. =',len(stream_rf)/3)

        kw = {'trim': (-2.5, 25),'fig_width':6, 'fillcolors': ('teal', 'darkgrey'), 'trace_height': 0.1, 'show_vlines': 'True','scale':2.5}
        stream_rf.select(component='R', station='{}'.format(station_5G)).sort(['distance']).plot_rf(**kw)
        plt.savefig('figs_rf_dist/{}_{}_.1_1_directPD.pdf'.format(station_5G,'R'),bbox_inches='tight', pad_inches=0.2)
        plt.close()

        ##
        ## following bit applies resonance filter if the station have sediments.
        #
        stream_rever_rmv=RFStream()
        auto_c=RFStream()
        if station_5G in thick_sedi_st:
            print('No rvr applied')
            stream_rever_rmv=stream_rf.select(component='R')
        else:
            stream_rever_rmv,auto_c=apply_reverberation_filter(stream_rf.select(component='R'))
            ### for plotting auto correlations saved in 'auto_c'
            ##
            fig = plt.figure(figsize=(4,3))
            # plt.style.use('grayscale')
            ax1 = fig.add_axes([0.1, 0.07, 0.85, 0.85]) #[left, bottom, width, height]
            ax1.set_ylim(26, 95)
            ax1.set_xlim(0,10)
            ax1.set_yticks(np.linspace(30,90,3))
            for i,tr in enumerate(auto_c):
                dist=tr.stats.distance
                auto=tr
                time = np.arange(auto.stats.npts) * auto.stats.delta
                auto.data /= np.max(np.abs(auto.data))
                r0=auto.stats.r0
                Dt=auto.stats.t3_offset#*auto.stats.sampling_rate

                ax1.plot(time,5*auto.data + dist,  lw=0.25, color='black')
                ax1.plot(Dt,5*-r0 + dist,marker='o', color='maroon',markersize=1,alpha=.8)

                # ax1.fill_between(time, i+1, auto.data + i+1, lw=0.55,
                                  # color=cm.viridis(auto.stats.tpdelay/1.32), where=(auto.data < 0),alpha=.8)
                # plt.setp(auto.stats.station,fontsize=3)
            ax1.xaxis.set_minor_locator(MultipleLocator(.25))
            ax1.yaxis.set_minor_locator(MultipleLocator(10))
            ax1.xaxis.set_major_locator(MultipleLocator(1))
            plt.setp(ax1.get_xticklabels(),fontsize=8)

            ax1.grid(which='major', axis='x',color='LightGrey', linestyle='--',linewidth=.45,alpha=.75)
            ax1.set_ylabel('Distance')
            ax1.set_xlabel('Time')
            ax1.set_title('AutoC for {}'.format(station_5G),fontsize=11)
            plt.savefig('figs_rf_dist/{}_{}_.1_1_autoC.pdf'.format(station_5G,'R'),bbox_inches='tight', pad_inches=0.2)
            # fig.show()
            plt.close()
            #
            ###
        stream_rever_rmv.write(rffile_ref, 'H5')

        ##################

        kw = {'trim': (-2.5, 25),'fig_width':6, 'fillcolors': ('cadetblue', 'gainsboro'), 'trace_height': 0.1, 'show_vlines': 'True','scale':2.5}
        stream_rever_rmv.select(component='R', station='{}'.format(station_5G)).sort(['distance']).plot_rf(**kw)
        plt.savefig('figs_rf_dist/{}_{}_.1_1_directP_rvrD.pdf'.format(station_5G,'R'),bbox_inches='tight', pad_inches=0.2)
        plt.close()


        # stack=stream_rf.stack()
        # stack.select(component='R')[0].stats.RF_Pcorr=len(stream_rf)/3
        # stack.select(component='R')[0].stats.RF_SNR=len(data_rf)/3
        # stack_5G.append(stack.select(component='R')[0])
        # stack_R=stack.select(component='R')[0].data
        # stack_Z=stack.select(component='Z')[0].data
        #
        # time_P= (np.argmax(stack_R)-np.argmax(stack_Z))/stack[0].stats.sampling_rate
        #
        # print("direct-P Time delay for station", station_5G, '=', time_P,'sec')
        # time_delay.append(time_P)
        # station.append(station_5G)
    except:
        print('Data/rf not found for station',station_5G)
    print('#############################################################')
