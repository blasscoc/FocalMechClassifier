from numpy import (pi, complex_, empty, mgrid, exp, sin, cos, zeros,
                   sqrt, mod, amin, amax, where)
from scipy.special import sph_harm as scipy_sph_harm


def sph_harm(m, l, longi, colat):
    """
    Adds the Cordon Shortly phase to the sph_harm.

    :parma m: order
    :param l: degree
    :param longi: longitude [0,2*pi]
    :param colat: colatitude [0,pi]
    """
    cordon_shortley = (-1)**m    
    return cordon_shortley * scipy_sph_harm(m, l, longi, colat)

# removes amplitude from a radiation pattern
def beachball(s):
    maxval= amax(amax(s))
    minval= amin(amin(s))
    return where(s>0,1,0)

def get_sph_harm(resolution=(25,25)):
    """ Initialize a grid of template harmonic functions
    """
    Nlong, Nlati = resolution
    
    longi, lati = mgrid[0:2*pi:(Nlati*1j), 0:pi:(Nlong*1j)]

    longi = longi.ravel()
    lati = lati.ravel()
    
    s_harm = empty(shape=(9,Nlong*Nlati), dtype=complex_)

    s_harm[0,:] = sph_harm(0, 0, longi, lati)
    
    s_harm[1,:] = sph_harm(-1, 1, longi, lati)
    s_harm[2,:] = sph_harm(0, 1, longi, lati)
    s_harm[3,:] = sph_harm(1, 1, longi, lati)
    
    s_harm[4,:] = sph_harm(-2, 2, longi, lati)
    s_harm[5,:] = sph_harm(-1, 2, longi, lati)
    s_harm[6,:] = sph_harm(0, 2, longi, lati)
    s_harm[7,:] = sph_harm(1, 2, longi, lati)
    s_harm[8,:] = sph_harm(2, 2, longi, lati)
    
    return longi, lati, s_harm

def WignerD1(alpha, beta, gamma):
    """
    Compute the Wigner-D matrix for l=1.
    
    Parameters - Euler angles    
    :param alpha: strike
    :param beta: dip
    :param gamma: -rake

    Returns:
    --------
    :rtype D: D^1_M'M(alpha,beta,gamma) : for -1 <= M <= 1    
    """

    D = zeros([3,3]) * complex(1,0)

    # exp(-i\alpha), exp(-i\gamma)
    neialp = exp(-complex(0,1)*alpha)
    neigam = exp(-complex(0,1)*gamma)
    # exp(i\alpha), exp(i\gamma)
    peialp = exp(complex(0,1)*alpha)
    peigam = exp(complex(0,1)*gamma)

    D[0,0] = neialp * 0.5 * (1.0 + cos(beta)) * neigam
    D[0,1] = -neialp * sin(beta) / sqrt(2.0)
    D[0,2] = neialp * 0.5 * (1.0 - cos(beta)) * peigam

    D[1,0] = sin(beta) * neigam / sqrt(2.0)
    D[1,1] = cos(beta)
    D[1,2] = -sin(beta) * peigam / sqrt(2.0)

    D[2,0] = peialp * 0.5 * (1.0 - cos(beta)) * neigam
    D[2,1] = peialp * sin(beta) / sqrt(2.0)
    D[2,2] = peialp * 0.5 * (1.0 + cos(beta)) * peigam

    return D

def WignerD2(alpha, beta, gamma):
    """
    Compute the wigner d matrix for l = 2.
    
    Parameters - Euler angles    
    :param alpha: strike
    :param beta: dip
    :param gamma: -rake

    Returns:
    :rtype D: D^2_M'M(alpha,beta,gamma) : for -2<=M<= 2
    """

    D = zeros([5,5]) * complex(1,0)

    # exp(-i\alpha), exp(-i\gamma)
    neialp2 = exp(-complex(0,2)*alpha)
    neigam2 = exp(-complex(0,2)*gamma)

    # exp(-i\alpha), exp(-i\gamma)
    neialp1 = exp(-complex(0,1)*alpha)
    neigam1 = exp(-complex(0,1)*gamma)

    # exp(-i\alpha), exp(-i\gamma)
    peialp1 = exp(complex(0,1)*alpha)
    peigam1 = exp(complex(0,1)*gamma)

    # exp(i\alpha), exp(i\gamma)
    peialp2 = exp(complex(0,2)*alpha)
    peigam2 = exp(complex(0,2)*gamma)

    D[0,0] = neialp2 * 0.25 * (1 + cos(beta)) * (1 + cos(beta)) * neigam2
    D[0,1] = -neialp2 * 0.5 * sin(beta) * (1 + cos(beta)) * neigam1
    D[0,2] = neialp2 * sqrt(3.0/8.0) * sin(beta) * sin(beta)
    D[0,3] = -neialp2 * 0.5 * sin(beta) * (1 - cos(beta)) * peigam1
    D[0,4] = neialp2 * 0.25 * (1 - cos(beta)) * (1 - cos(beta)) * peigam2

    D[1,0] = neialp1 * 0.5 * sin(beta) * (1 + cos(beta)) * neigam2
    D[1,1] = neialp1 * 0.5 * (2.0*cos(beta)**2 + cos(beta) - 1) * neigam1
    D[1,2] = -neialp1 * sqrt(3.0/2.0) * sin(beta) * cos(beta) 
    D[1,3] = -neialp1 * 0.5 * (2.0*cos(beta)**2 - cos(beta) - 1) * peigam1
    D[1,4] = -neialp1 * 0.5 * sin(beta) * (1 - cos(beta)) * peigam2

    D[2,0] = sqrt(3.0/8.0) * sin(beta) * sin(beta) * neigam2
    D[2,1] = sqrt(3.0/2.0) * sin(beta) * cos(beta) * neigam1
    D[2,2] = 1.5 * cos(beta) * cos(beta) - 0.5
    D[2,3] = -sqrt(3.0/2.0) * sin(beta) * cos(beta) * peigam1
    D[2,4] = sqrt(3.0/8.0) * sin(beta) * sin(beta) * peigam2

    D[3,0] = peialp1 * 0.5 * sin(beta) * (1 - cos(beta)) * neigam2
    D[3,1] = -peialp1 * 0.5 * (2.0*cos(beta)**2 - cos(beta) - 1) * neigam1
    D[3,2] = peialp1 * sqrt(3.0/2.0) * sin(beta) * cos(beta)
    D[3,3] = peialp1 * 0.5 * (2.0*cos(beta)**2 + cos(beta) - 1) * peigam1
    D[3,4] = -peialp1 * 0.5 * sin(beta) * (1 + cos(beta)) * peigam2

    D[4,0] = peialp2 * 0.25 * (1 - cos(beta)) * (1 - cos(beta)) * neigam2
    D[4,1] = peialp2 * 0.5 * sin(beta) * (1 - cos(beta)) * neigam1
    D[4,2] = peialp2 * sqrt(3.0/8.0) * sin(beta) * sin(beta)
    D[4,3] = peialp2 * 0.5 * sin(beta) * (1 + cos(beta)) * peigam1
    D[4,4] = peialp2 * 0.25 * (1 + cos(beta)) * (1 + cos(beta)) * peigam2

    return D.conjugate()
