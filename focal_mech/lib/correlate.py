
from numpy import (pi, rad2deg, inner, array, zeros, mgrid,
                       argmin)
from scipy.linalg import norm
from scipy.optimize import minimize

from focal_mech.lib.sph_harm import WignerD2


def _corr_shear(x, alm):
    """
    Rotates the (13),(31) double couple source you can 
    find in Ben-Menahem and Singh. Monopole, dipole 
    terms for this source are zero.
    """
    
    strike, dip, rake = x
    # Wigner is ZYZ Euler rotation, \gamma = -rake
    D = WignerD2(strike, dip, -rake).conjugate()
    # Template Spectrum : glm = (0, -1j, 0, -1j, 0)
    prop = (inner(D[:,3], alm) + inner(D[:,1], alm))
    # Maximize, not minimize.
    return -norm(prop)

def _scan_shear(alm):
    """
    :param alm: Quadrupole components. 

    Scans for the best starting point, the optimization
    can get stuck otherwise.
    """
   
    # scan for a good starting point
    X, Y, Z = mgrid[0:2*pi:10j, 0:pi:10j, 0:2*pi:10j]
    x0s = zip(X.ravel(),Y.ravel(),Z.ravel())
    res = [_corr_shear(x0, alm) for x0 in x0s]
    x0 = x0s[argmin(res)]

    print x0
    
    return x0
    
def corr_shear(Alm):
    """
    :param Alm: the spectrum of the classifier function. 

    Scans for the best starting point, the optimization
    can get stuck otherwise.
    """

    # correlation with iso, dipole is zero.
    alm = array([Alm[2,-2], Alm[2,-1],
                 Alm[2,0],
                 Alm[2,1],Alm[2,2]])
    
    # pick a good starting point.
    x0 = _scan_shear(alm)

    f = lambda x : _corr_shear(x,alm)
    results = minimize(f, x0=x0,
                    bounds=((0,2*pi), (0,pi), (0,2*pi)))

    return rad2deg(results.x), results.fun
    
def corr_tensile(Alm):
    raise Exception("Method not implemented")
