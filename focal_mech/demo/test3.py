from numpy import genfromtxt, rad2deg

import matplotlib.pyplot as plt
import mplstereonet

from obspy.imaging.beachball import aux_plane

from focal_mech.io.read_hash import read_demo, read_hash_solutions

event = 3146815

hash_solns = read_hash_solutions("example1.out")
hash_focal = rad2deg(hash_solns[event])

mechs = genfromtxt("3146815_realizations.out")[:,:3]
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
                    
polarity = polarity_data[event][:,2]
for a, t, p in zip(azi, toa, polarity):
    if p > 0:
        ax.pole(a, t,'o', markeredgecolor='red', markerfacecolor='red')
    else:
        ax.pole(a, t,'o', markeredgecolor='blue', markerfacecolor='white')

fig.suptitle("Event %d" % event)
ax.grid()

plt.show()
        