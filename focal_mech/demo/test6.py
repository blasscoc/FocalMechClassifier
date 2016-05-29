from numpy import array, rad2deg, pi, mgrid, argmin

from matplotlib.pylab import contour
import matplotlib.pyplot as plt
import mplstereonet

from obspy.imaging.beachball import aux_plane

from focal_mech.lib.classify_mechanism import classify, translate_to_sphharm
from focal_mech.io.read_hash import read_demo, read_hash_solutions

from focal_mech.util.hash_routines import hash_to_classifier
from focal_mech.lib.sph_harm import get_sph_harm
from focal_mech.lib.correlate import corr_shear

hash_solns = read_hash_solutions("example1.out")

# we want solutions that are symetric
polarity_data = read_demo("north1.phase", "scsn.reverse", reverse=True)
inputs = hash_to_classifier(polarity_data, parity=1)

event = 3146815

result = classify(*inputs[event], kernel_degree=2)
Alm = translate_to_sphharm(*result, kernel_degree=2)

coeffs = array([Alm[0,0], 
                Alm[1,-1], Alm[1,0], Alm[1,1], 
                Alm[2,-2], Alm[2,-1], Alm[2,0], Alm[2,1], Alm[2,2]])


svm_soln, f = corr_shear(Alm)

resolution = (200,400)
longi, lati, Z = get_sph_harm(resolution=resolution)
mech = coeffs.dot(Z).real

longi.shape = resolution
lati.shape = resolution
mech.shape = resolution

c = contour(longi, lati, mech, [0])
pth1 = c.collections[0].get_paths()[0].vertices
pth1 = rad2deg(pth1)

pth2 = c.collections[0].get_paths()[1].vertices
pth2 = rad2deg(pth2)

hash_focal = rad2deg(hash_solns[event])


event2 = 3158361


result = classify(*inputs[event2], kernel_degree=2)
Alm = translate_to_sphharm(*result, kernel_degree=2)

coeffs = array([Alm[0,0], 
                Alm[1,-1], Alm[1,0], Alm[1,1], 
                Alm[2,-2], Alm[2,-1], Alm[2,0], Alm[2,1], Alm[2,2]])

svm_soln2, f = corr_shear(Alm)

resolution = (200,400)
longi, lati, Z = get_sph_harm(resolution=resolution)
mech = coeffs.dot(Z).real

longi.shape = resolution
lati.shape = resolution
mech.shape = resolution

c = contour(longi, lati, mech, [0])
pth3 = c.collections[0].get_paths()[0].vertices
pth3 = rad2deg(pth3)

pth4 = c.collections[0].get_paths()[1].vertices
pth4 = rad2deg(pth4)

hash_focal2 = rad2deg(hash_solns[event2])


event3 = 3153955


result = classify(*inputs[event3], kernel_degree=2)
Alm = translate_to_sphharm(*result, kernel_degree=2)

coeffs = array([Alm[0,0], 
                Alm[1,-1], Alm[1,0], Alm[1,1], 
                Alm[2,-2], Alm[2,-1], Alm[2,0], Alm[2,1], Alm[2,2]])

svm_soln3, f = corr_shear(Alm)

resolution = (200,400)
longi, lati, Z = get_sph_harm(resolution=resolution)
mech = coeffs.dot(Z).real

longi.shape = resolution
lati.shape = resolution
mech.shape = resolution

c = contour(longi, lati, mech, [0])
pth5 = c.collections[0].get_paths()[0].vertices
pth5 = rad2deg(pth5)

pth6 = c.collections[0].get_paths()[1].vertices
pth6 = rad2deg(pth6)

hash_focal3 = rad2deg(hash_solns[event3])


fig = plt.figure(facecolor="white", figsize=(10,20))
ax = fig.add_subplot(221, projection='stereonet')

ax.rake(pth1[:,0], pth1[:,1] +90.0, 90.0, ':', color='red', linewidth=3)
ax.rake(pth2[:,0], pth2[:,1] +90.0, 90.0, ':', color='red', linewidth=3)

strike, dip, rake = svm_soln
ax.plane(strike, dip, '-r', linewidth=2)

strike, dip, rake = aux_plane(*svm_soln)
ax.plane(strike, dip, '-r', linewidth=2)

strike, dip, rake = hash_focal
ax.plane(strike-90, dip, 'g-', linewidth=2)

strike, dip, rake = aux_plane(*hash_focal)
ax.plane(strike-90, dip,'g-', linewidth=2)

azi = rad2deg(polarity_data[event][:,0])
toa = rad2deg(polarity_data[event][:,1])

polarity = polarity_data[event][:,2]
for a, t, p in zip(azi, toa, polarity):
    if p > 0:
        ax.pole(a, t,'o', markeredgecolor='red', markerfacecolor='red')
    else:
        ax.pole(a, t,'o', markeredgecolor='blue', markerfacecolor='white')

ax.grid()


ax = fig.add_subplot(222, projection='stereonet')

ax.rake(pth3[:,0], pth3[:,1] +90.0, 90.0, ':', color='red', linewidth=3)
ax.rake(pth4[:,0], pth4[:,1] +90.0, 90.0, ':', color='red', linewidth=3)

strike, dip, rake = svm_soln2
ax.plane(strike, dip, '-r', linewidth=2)

strike, dip, rake = aux_plane(*svm_soln2)
ax.plane(strike, dip, '-r', linewidth=2)

strike, dip, rake = hash_focal2
ax.plane(strike-90, dip, 'g-', linewidth=2)

strike, dip, rake = aux_plane(*hash_focal2)
ax.plane(strike-90, dip,'g-', linewidth=2)

azi = rad2deg(polarity_data[event2][:,0])
toa = rad2deg(polarity_data[event2][:,1])

polarity = polarity_data[event2][:,2]
for a, t, p in zip(azi, toa, polarity):
    if p > 0:
        ax.pole(a, t,'o', markeredgecolor='red', markerfacecolor='red')
    else:
        ax.pole(a, t,'o', markeredgecolor='blue', markerfacecolor='white')

ax.grid()


ax = fig.add_subplot(224, projection='stereonet')

ax.rake(pth5[:,0], pth5[:,1] +90.0, 90.0, ':', color='red', linewidth=3)
ax.rake(pth6[:,0], pth6[:,1] +90.0, 90.0, ':', color='red', linewidth=3)

strike, dip, rake = svm_soln3
ax.plane(strike, dip, '-r', linewidth=2)

strike, dip, rake = aux_plane(*svm_soln3)
ax.plane(strike, dip, '-r', linewidth=2)

strike, dip, rake = hash_focal3
ax.plane(strike-90, dip, 'g-', linewidth=2)

strike, dip, rake = aux_plane(*hash_focal3)
ax.plane(strike-90, dip,'g-', linewidth=2)

azi = rad2deg(polarity_data[event3][:,0])
toa = rad2deg(polarity_data[event3][:,1])

polarity = polarity_data[event3][:,2]
for a, t, p in zip(azi, toa, polarity):
    if p > 0:
        ax.pole(a, t,'o', markeredgecolor='red', markerfacecolor='red')
    else:
        ax.pole(a, t,'o', markeredgecolor='blue', markerfacecolor='white')

ax.grid()


plt.tight_layout(pad=4.0, h_pad=20.0)

plt.show()
