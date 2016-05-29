from numpy import genfromtxt, rad2deg

import matplotlib.pyplot as plt
import mplstereonet

from obspy.imaging.beachball import aux_plane

from focal_mech.lib.classify_mechanism import classify, translate_to_sphharm
from focal_mech.io.read_hash import read_demo, read_hash_solutions

from focal_mech.util.hash_routines import hash_to_classifier
from focal_mech.lib.sph_harm import get_sph_harm
from focal_mech.lib.correlate import corr_dc, _corr_dc

event = 3145744

hash_solns = read_hash_solutions("example1.out")
hash_focal = rad2deg(hash_solns[event])

mechs = genfromtxt("3145744_realizations.out")[:,:3]
polarity_data = read_demo("north1.phase", "scsn.reverse", reverse=True)

fig = plt.figure(facecolor="white")
ax = fig.add_subplot(111, projection='stereonet')


for mech in mechs[:10]:
    strike, dip, rake = mech
    ax.plane(strike-90, dip, '-r', linewidth=2)
    
    strike, dip, rake = aux_plane(*mech)
    ax.plane(strike-90, dip, '-r', linewidth=2)
     
strike, dip, rake = hash_focal
ax.plane(strike-90, dip, '-b', linewidth=2)
    
strike, dip, rake = aux_plane(*hash_focal)
ax.plane(strike-90, dip, '-b', linewidth=2)
                 
azi = rad2deg(polarity_data[event][:,0])
toa = rad2deg(polarity_data[event][:,1])
                    
polarity = polarity_data_nr[event][:,2]
for a, t, p in zip(azi, toa, polarity):
    if p > 0:
        ax.pole(a, t,'o', markeredgecolor='red', markerfacecolor='red')
    else:
        ax.pole(a, t,'o', markeredgecolor='blue', markerfacecolor='white')

fig.suptitle("Event %d" % event)
ax.grid()

plt.show()
        