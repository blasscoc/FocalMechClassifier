from numpy import (dot, array, sign, arccos, zeros, pi, c_, arctan2,
                         sum, pi, sqrt, conjugate, logical_and, arcsin)
from math import gamma

from sklearn import svm

from focal_mech.lib.sph_harm import sph_harm


def kernel(x, y, degree=2, coeff=1):
    """
    Evaluates the polynomial kernel used by sklearn, (x.dot(y) + coeff)^degree

    Parameters:
    -----------
    :param x,y: cartesian coords.
    :param degree, coeff: the degree of the kernel, 2 and 1 respectively 
                          is adequate for iso, dc and clvd classification.                             
    """
    return (dot(x,y) + coeff)**degree
    
def classify(*args, **kwargs):
    """Parameters:
    -----------
    :param args: (x, y, z, data) of each observation 

    :param kwargs: kernel_degree and kernel_coeff parameterize the
    kernels, 2 and 1 respectively is adequate for iso, dc and clvd
    classification.

    Returns:
    --------
    :rtype dual_coeff:  dual coefficient 
    :rtype theta: azimuth of the data point
    :rtype phi: polar angle of the data point

    Notes: 
    ------ 
    We're using the sklearn support vector machine to
    find a optimal seperating hyperplanes classifying the polarity of
    points on a sphere. Our particular choice of kernel can be shown
    to generate a expansion in terms of a basis spherical harmonics up
    to order equal to the kernel degree.

    Isotropic, double couple and clvd sources can be represented in
    this basis, which makes this a convenient expansion.

    """
    
    if 'kernel_degree' in kwargs:
        kernel_degree = kwargs['kernel_degree']
    else:
        kernel_degree = 2

    if 'kernel_coeff' in kwargs:
        kernel_coeff = kwargs['kernel_coeff']
    else:
        kernel_coeff = 1

    x, y, z, data = args
    inputs = array([x.ravel(),y.ravel(),z.ravel()]).T
    classes = sign(data.ravel())
    classes[classes<=0] = -1

    poly_svc = svm.SVC(kernel='poly', degree=kernel_degree, 
                       coef0=kernel_coeff, C=1.0).fit(inputs, classes)

    n_support = poly_svc.n_support_
    if len(n_support) > 2:
        raise Exception("Not supported yet.")

    # \beta_0
    intercept = poly_svc.intercept_

    # we only have classes : 1 or 0 (y_i ~ \pm 1)    
    dual_coeff = poly_svc.dual_coef_[0,:]

    num_support = len(dual_coeff)

    #we need the angles of the support vectors:
    # measured from up - the colatitude
    theta = arccos(poly_svc.support_vectors_[:,2])
    
    # the forward half, phi > 0
    phi = zeros(num_support)

    # phi the angle between [-pi,pi] measured as azimuth
    phi = arctan2(poly_svc.support_vectors_[:,1],poly_svc.support_vectors_[:,0])
    
    in_sample = poly_svc.predict(c_[inputs])

    return dual_coeff, phi, theta, intercept, in_sample

def scholkopf_norm(elle,kernel_degree):
    # Scholkopf and Smola, page 113, eq 4.72 coefficient a_l,l > kernel_degree=0
    norm = ((2.0**(kernel_degree+1) *
                 gamma(kernel_degree+1) * gamma(kernel_degree+1)) / 
                   (gamma(kernel_degree+2+elle) * gamma(kernel_degree+1-elle)))
    return norm

def calc_00(intercept, kernel_degree=2):
    """
    The hyper-plane is:
    f(x) = \sum( ... ) + \beta
    beta is the bias, which is the overlap with the 00 mode.
    """
    norm = scholkopf_norm(0,kernel_degree)
    # from the addition theorem
    norm *= 4 * pi * 0.5 * sqrt(1.0/pi)

    return intercept  / norm

def calc_alm(dual_coeff, phi, theta, elle, emm, kernel_degree=2):
    """Computes the spectral content of the classifier function estimated by 
    the support vector machine. 

    Alm = 4\pi/(2l+1) * \sum_{i} \dual_coeff[i] a[elle] Y*_{elle,emm}(phi,theta)

    Parameters:
    -----------

    :param dual_coeff : the lagrange multipliers * class (-1,1) output
                 of the support vector machine.  e.g. see sklearn.svm
                 dual_coef_

    :param phi: the azimuth of each observation
    :param theta: the polar angle of each observation

    :param elle: "l" as in (Y_lm), the degree degree of the spherical harmonic
    :param emm: "m" as in (Y_lm) the order of the spherical harmonic, -l<=m<= l  

    :param kernel_degree: must be an integer, for the kernel used
    (<x,x'> + 1)^d, use the default for focal mechanism
    classification.
                       
    Returns:
    --------
    :rtype Alm: A dict Alm[m,l] containing the spectrum in the sph_harm basis.

    Notes:
    ------    
    A kernel degree of 2 truncates the expansion l=2 and so is
    recommended for clvd, dc, iso classification.
    """

    # l = elle, and the kernel degree is 2
    norm = scholkopf_norm(elle,kernel_degree)
    # from the addition theorem
    norm *= 4 * pi / (2.0 * elle + 1)

    # scipy does y_ml (with -l <= m <= l)
    # rather than y_lm (convention used in Jackson)
    return sum(conjugate(sph_harm(emm, elle, phi, theta) ) * dual_coeff) * norm

def translate_to_sphharm(*args, **kwargs):
    """
    
    Take the results of the classification (using kernels) and
    translate the result to a basis of spherical harmonics.

    :param args: pass the output of classify
    :rtype Alm: the spectrum of the classifier function.
    """

    kernel_degree = kwargs.pop("kernel_degree",2)
    
    dual_coeff, phi, theta, intercept, _ = args
    
    # Scholkopf and Smola, page 113, eq 4.72 coefficient a_l, l>kernel_degree=0
    Alm = {}
    Alm[0,0] = calc_00(intercept)
    Alm[0,0] = Alm[0,0][0]
    for l in range(1,kernel_degree+1):
        for m in range(-l,l+1):            
            Alm[(l,m)] = calc_alm(dual_coeff, phi, theta, l, m)
    
    return Alm

