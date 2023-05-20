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
