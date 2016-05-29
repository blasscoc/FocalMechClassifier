from numpy import array, rad2deg, logical_and, rad2deg, pi, mgrid, argmin

from matplotlib.pylab import contour
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
hash_solns_nr = read_hash_solutions("example1-noreverse.out")

# we want solutions that are symetric
polarity_data = read_demo("north1.phase", "scsn.reverse", reverse=True)
polarity_data_nr = read_demo("north1.phase", "scsn.reverse", reverse=False)

inputs = hash_to_classifier(polarity_data, parity=1)
inputs_nr = hash_to_classifier(polarity_data_nr, parity=1)


result = classify(*inputs[event], kernel_degree=2)
Alm = translate_to_sphharm(*result, kernel_degree=2)

coeffs = array([Alm[0,0], 
                Alm[1,-1], Alm[1,0], Alm[1,1], 
                Alm[2,-2], Alm[2,-1], Alm[2,0], Alm[2,1], Alm[2,2]])

result_nr = classify(*inputs_nr[event], kernel_degree=2)
Alm_nr = translate_to_sphharm(*result_nr, kernel_degree=2)

coeffs_nr = array([Alm_nr[0,0], 
                Alm_nr[1,-1], Alm_nr[1,0], Alm_nr[1,1], 
                Alm_nr[2,-2], Alm_nr[2,-1], Alm_nr[2,0], Alm_nr[2,1], Alm_nr[2,2]])


X, Y, Z = mgrid[0:2*pi:10j, 0:pi:10j, 0:2*pi:10j]
x0s = zip(X.ravel(),Y.ravel(),Z.ravel())
res = [_corr_dc(x0,coeffs[4:]) for x0 in x0s]
svm_soln, f = corr_dc(Alm, x0=x0s[argmin(res)])

res = [_corr_dc(x0,coeffs_nr[4:]) for x0 in x0s]
svm_soln_nr, f = corr_dc(Alm_nr, x0=x0s[argmin(res)])


resolution = (200,400)
longi, lati, Z = get_sph_harm(resolution=resolution)
mech = coeffs.dot(Z).real

longi.shape = resolution
lati.shape = resolution
mech.shape = resolution

#c = contour(longi, lati, mech, [0])
#pth1 = c.collections[0].get_paths()[0].vertices
#pth1 = rad2deg(pth1)

#pth2 = c.collections[0].get_paths()[1].vertices
#pth2 = rad2deg(pth2)

hash_focal = rad2deg(hash_solns[event])
hash_focal_nr = rad2deg(hash_solns_nr[event])

fig = plt.figure(facecolor="white")
ax = fig.add_subplot(111, projection='stereonet')

indx = mech > 0
# reflect the points in the lower half space
#ax.pole(rad2deg(longi[indx]), rad2deg(lati[indx]))

#ax.rake(pth1[:,0], pth1[:,1] +90.0, 90.0, '-', color='red', linewidth=3)
#ax.rake(pth2[:,0], pth2[:,1] +90.0, 90.0, '-', color='red', linewidth=3)


strike, dip, rake = svm_soln
ax.plane(round(strike), dip, '-r', linewidth=2)

strike, dip, rake = aux_plane(*svm_soln)
ax.plane(strike, dip, '-r', linewidth=2)

strike, dip, rake = svm_soln_nr
ax.plane(round(strike), dip, '--r', linewidth=2)

strike, dip, rake = aux_plane(*svm_soln_nr)
ax.plane(strike, dip, '--r', linewidth=2)



strike, dip, rake = hash_focal
ax.plane(strike-90, dip, 'g-', linewidth=2)

strike, dip, rake = aux_plane(*hash_focal)
ax.plane(strike-90, dip,'g-', linewidth=2)


strike, dip, rake = hash_focal_nr
ax.plane(strike-90, dip, 'g--', linewidth=2)

strike, dip, rake = aux_plane(*hash_focal_nr)
ax.plane(strike-90, dip,'g--', linewidth=2)


azi = rad2deg(polarity_data_nr[event][:,0])
toa = rad2deg(polarity_data_nr[event][:,1])

# https://github.com/markcwill/hashpy/blob/master/hashpy/plotting/focalmechplotter.py #
# Calculate strike azi from direct (dip-pointing) azi 
#azi = azi - 90.
#--- HASH takeoffs are 0-180 from vertical UP!!
#--- Stereonet angles 0-90 inward (dip)
#--- Classic FM's are toa from center???
#indx = logical_and(toa >= 0, toa < 90)
#toa[indx] = 90 - toa[indx]
#indx = logical_and(toa >= 90, toa < 180)
#toa[indx] = 270 - toa[indx]


polarity = polarity_data[event][:,2]
for a, t, p in zip(azi, toa, polarity):
    if p > 0:
        ax.pole(a, t,'o', markeredgecolor='red', markerfacecolor='red')
    else:
        ax.pole(a, t,'o', markeredgecolor='blue', markerfacecolor='white')

fig.suptitle("Event %d" % event)
ax.grid()

plt.show()
