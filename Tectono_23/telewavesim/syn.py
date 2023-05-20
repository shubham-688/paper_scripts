from telewavesim.utils import Model
from telewavesim import utils as ut
from telewavesim import wiggle as wg

import matplotlib.pyplot as plt
import obspy
from obspy import read, Stream, UTCDateTime,read_events
from obspy.core.event import Origin
from obspy.taup import TauPyModel
import numpy as np
import sys

import funcs as fn
from importlib import reload
reload(fn)

#########
thick = [.02, 42.7, 0]       # Second layer thickness is irrelevant
rho = [2400,3000., 3300.]   # Second rho value is irrelevant as we use a pre-defined elastic tensor
vp = [3.8,6.45, 7.9]         # Likewise for vp
vs = [1.6,3.76, 4.6]        # Likewise for vs
flag = ['iso','iso', 'iso']  # Both layers are isotropic
model = Model(thick, rho, vp, vs, flag)


# sys.exit()
# Define three-layer model with isotropic crust and antigorite upper mantle layer over isotropic half-space
# model = utils.Model([2, 40, 0], [2800., None, 3300.], [4.6, 0, 6.0], [2.6, 0, 3.6], ['iso', 'iso', 'iso'], [0, 0, 0], [0, 0, 0], [0, 0, 0])
# for th in [2,2.2,2.4]:
# thick = [.02, 42.7, 0]

model = Model(thick, rho, vp, vs, flag)
slow = 0.06     # s/km
baz = np.arange(0., 360., 10.)
npts = 3000
dt = 0.025  # s
wvtype = 'P'
# st = utils.run_plane(model, slow, npts, dt,wvtype=wvtype)
######
trR = Stream(); trT = Stream()

for bb in baz:
    # Calculate the plane waves seismograms
    trxyz = ut.run_plane(model, slow, npts, dt, bb, wvtype=wvtype, obs=False)

    # Then the transfer functions in Z-R-T coordinate system
    tfs = ut.tf_from_xyz(trxyz, pvh=False)

    # Append to streams
    trR.append(tfs[0]); trT.append(tfs[1])

# Set frequency corners in Hz
f1 = 0.01
f2 = 1.0

# Filter to get wave-like traces
trR.filter('bandpass',freqmin=f1, freqmax=f2, corners=2, zerophase=True)
trT.filter('bandpass',freqmin=f1, freqmax=f2, corners=2, zerophase=True)

trR_stack, trT_stack = ut.stack_all(trR, trT, pws=True)

# Plot as wiggles
fn.rf_wiggles_baz_r(trR, trT, trR_stack, trT_stack, 'test', btyp='baz',
                  scale=1.e3, tmin=-5., tmax=18., save=True, ftitle='s07_{}km_3.8_1.6_2200'.format(0),
                  wvtype='P')
# plt.show()
