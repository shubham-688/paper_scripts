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
