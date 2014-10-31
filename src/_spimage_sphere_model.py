import numpy
import logging
import scipy.signal
from scipy.optimize import leastsq
import scipy.stats
import spimage
__all__ = ['fit_sphere_diameter', 'fit_sphere_model', 'get_sphere_model_size', 'get_sphere_model_scale']


def fit_sphere_diameter(img, msk, method=None, full_output=False, **kwargs):
    """
    Fit the diameter of a sphere to diffraction data.

    usage:
    ======
    diameter = fit_sphere_diameter(img, msk)
    diameter = fit_sphere_diameter(img, msk, method='pearson', ...)

    """
    # Default method for fitting of diameter
    if method is None: method = 'pearson'

    # Fit the diameter using Pearson correlation
    if method == 'pearson':
        diameter, info = fit_sphere_diameter_pearson(img, msk, full_output=True, **kwargs)
    else:
        diameter, info = [1., "There is no fitting diameter method %s" %method]
        print info

    # Return with/without full information
    if full_output:
        return diameter, info
    else:
        return diameter


def fit_sphere_diameter_pearson(img, msk,x0=0, y0=0, rmax=None, downsampling=1, do_brute=True, full_output=False):
    """
    Fit the diameter of a sphere using pearson correlation.

    """
    # Handle defaults
    if rmax is None: rmax = max(img.shape)/2

    # Downsample image
    X,Y = numpy.meshgrid(numpy.arange(0.,img.shape[1],1.),numpy.arange(0.,img.shape[0],1.))
    msk = msk[::downsampling,::downsampling]
    s = img.shape
    cx = (s[1]-1)/2.+x0
    cy = (s[0]-1)/2.+y0
    print X-cx
    print Y-cy
    print spimage.grid(s, (y0,x0))[1]
    print spimage.grid(s, (y0,x0))[0]
    Rsq = (X-cx)**2+(Y-cy)**2
    Mr = (rmax**2)>=Rsq
    Mr = Mr[::downsampling,::downsampling]
    Xm = X[msk*Mr]
    Ym = Y[msk*Mr]

    # First fit the radius
    # Start with brute force with a sensible range
    # We'll assume at least 20x oversampling
    if(do_brute):
        pixel_size = (D*wavelength/(p * image.shape[0]))
        range = [(pixel_size, pixel_size*image.shape[0]/(40))]
        Ns = 10 * range[0][1]/range[0][0]
        r = scipy.optimize.brute(err, range, Ns=Ns)[0]

    # End with least square
    # r, success = leastsq(err, r, maxfev=maxfev,xtol=1e-3)
    r, cov_r, infodict, mesg, ier = leastsq(err, r, maxfev=maxfev,xtol=1e-3, full_output=1)
    # print infodict
    r = r[0]



#from pylab import *
def fit_sphere_model(image,mask,params,r_max,downsampling=1,do_brute=True):
    # Downsample image
    X,Y = numpy.meshgrid(numpy.arange(0.,image.shape[1],1.),numpy.arange(0.,image.shape[0],1.))
    mask = mask[::downsampling,::downsampling]
    s = image.shape
    cx = (s[1]-1)/2.+params["offCenterX"]
    cy = (s[0]-1)/2.+params["offCenterY"]
    Rsq = (X-cx)**2+(Y-cy)**2
    Mr = (r_max**2)>=Rsq
    Mr = Mr[::downsampling,::downsampling]
    Xm = X[mask*Mr]
    Ym = Y[mask*Mr]
    
    #imsave("img.png",log10(image))
    #imsave("mask.png",mask)
    #imsave("Mr.png",Mr)

    p = params["detectorPixelSizeUM"]*1.E-6
    D = params["detectorDistanceMM"]*1.E-3
    wavelength = params["photonWavelengthNM"]*1.E-9
    h = spimage.DICT_physical_constants['h']
    c = spimage.DICT_physical_constants['c']
    ey_J = h*c/wavelength
    d = params["diameterNM"]*1.E-9
    I0 = params["intensityMJUM2"]*1.E-3*10E12/ey_J
    Mat = spimage.Material(material_type=params["materialType"])
    rho_e = Mat.get_electron_density()
    r = d/2.
    V = 4/3.*numpy.pi*r**3

    fitimg = image[mask*Mr]/params["detectorADUPhoton"]/params["detectorQuantumEfficiency"]

    #q = generate_absqmap(X-cx,Y-cy,p,D,wavelength)
    #I_fit = lambda K,r: I_sphere_diffraction(K,q,r)    
    qm = generate_absqmap(Xm-cx,Ym-cy,p,D,wavelength)
    I_fit_m = lambda K,r: I_sphere_diffraction(K,qm,r)

    # v[0]: K, v[1]: r
    #i_fit = lambda v: I_fit(v[0],v[1])
    i_fit_m = lambda v: I_fit_m(v[0],v[1])

    K = I0 * ( rho_e*p/D*spimage.DICT_physical_constants["re"]*V )**2

    err = lambda v: 1-scipy.stats.pearsonr(i_fit_m([K,v]),fitimg)[0]
    maxfev=1000 

    # First fit the radius
    # Start with brute force with a sensible range
    # We'll assume at least 20x oversampling
    if(do_brute):
        pixel_size = (D*wavelength/(p * image.shape[0]))
        range = [(pixel_size, pixel_size*image.shape[0]/(40))]
        Ns = 10 * range[0][1]/range[0][0]
        r = scipy.optimize.brute(err, range, Ns=Ns)[0]
    # End with least square
    # r, success = leastsq(err, r, maxfev=maxfev,xtol=1e-3)
    r, cov_r, infodict, mesg, ier = leastsq(err, r, maxfev=maxfev,xtol=1e-3, full_output=1)
    # print infodict
    r = r[0]

    # Now fit the intensity
    err2 = lambda v: i_fit_m([v,r])-fitimg
    K, success = leastsq(err2, K, maxfev=maxfev, xtol=1e-3)

    v1 = [K,r]
    Vnew = 4/3.*numpy.pi*v1[1]**3
    I0 = v1[0] / (rho_e*p/D*spimage.DICT_physical_constants["re"]*Vnew)**2
    params["intensityMJUM2"] = I0/1.E-3/1.E12*ey_J
    params["diameterNM"] = v1[1]*2/1.E-9
    #print params
    return params

# scattering amplitude from homogeneous sphere:
# -----------------------------------------------
# Source:
# Feigin 1987
#
# r: sphere radius
#
# F = sqrt(I_0) rho_e p/D r_0 4/3 pi r^3 [ 3 { sin(qr) - qr cos(qr) } / (qr)^3 ]
#   = sqrt(I_0) rho_e p/D r_0 V f(r,qx,qy)
# f = 3 { sin(qr) - qr cos(qr) } / (qr)^3
# K = I_0 (rho_e p/D r_0 V)^2
# S = I_0 rho_e^2 = K / (p/D r_0 V)^2
# ============================================================================================
# I = F^2 = K [ f(r,qx,qy) ]^2
# ============================================================================================
_F_sphere_diffraction = lambda K,q,r: numpy.sqrt(abs(K))*3*(numpy.sin(q*r)-q*r*numpy.cos(q*r))/((q*r)**3+numpy.finfo("float64").eps)
F_sphere_diffraction = lambda K,q,r: ((q*r)**6 < numpy.finfo("float64").resolution)*numpy.sqrt(abs(K)) + ((q*r)**6 >= numpy.finfo("float64").resolution)*_F_sphere_diffraction(K,q,r)
_I_sphere_diffraction = lambda K,q,r: abs(K)*(3*(numpy.sin(q*r)-q*r*numpy.cos(q*r))/((q*r)**3+numpy.finfo("float64").eps))**2
I_sphere_diffraction = lambda K,q,r: ((q*r)**6 < numpy.finfo("float64").resolution)*abs(K) + ((q*r)**6 >= numpy.finfo("float64").resolution)*_I_sphere_diffraction(K,q,r)

def generate_absqmap(X,Y,p,D,wavelength):
    R_Ewald = 2*numpy.pi/(1.*wavelength)
    qx = R_Ewald*p*X/D
    qy = R_Ewald*p*Y/D
    q_map = numpy.sqrt(qx**2+qy**2)
    return q_map      


def get_sphere_model_size(diameter, wavelength, pixelsize, detector_distance):
    """
    returns the size of the sphere model:

    modelsize = (diameter/2.) * (2pi / wavelength) * pixelsize * (1/detector_distance)

    Parameters:
    ===========
    diameter:          The diameter of the sphere in [nm]
    wavelength:        The wavelength of the xray beam in [nm]
    pixelsize:         The size of one detector pixel in [mu m]
    detector_distance: The distance from the interaction to the detector in [mm]
    """    
    # convert to SI units
    r = diameter*1.e-9/2. 
    wl = wavelength*1.e-9
    p = pixelsize*1.e-6
    D = detector_distance*1.e-3
    # compute modelsize
    k = 2*numpy.pi/wl
    modelSize = r*k*p/D
    return modelSize

def get_sphere_model_scale(intensity, diameter, wavelength, pixelsize, detector_distance, detector_efficiency, detector_adu_photon, material):
    """
    returns the scale of the sphere model:

    V = (4/3 * pi * diameter/2 ) **3
    I_0 = (intensity / (hc / wavelength) ) * 1e12
    K = I_0 * ( rho_e * pixelsize / detector_distance * re * V )**2
    scale = K * detector_adu_photon * detector_efficiency

    Parameters:
    ===========
    intensity:           The photon intensity on the sample in [mJ/mu m]
    diameter:            The diameter of the sphere in [nm]
    wavelength:          The wavelength of the xray beam in [nm]
    pixelsize:           The size of one detector pixel in [mu m]
    detector_distance:   The distance from the interaction to the detector in [mm]
    detector_efficiency: The quantum efficieny of the detector
    detector_adu_photon: The nr. ADU units per photon 
    material:            The material of the sample -> defines electron_density rho_e

    Physical constants used:
    ========================
    h:  Planck constant (6.62606957E-34)
    c:  Speed of light (299792458)
    qe: Electron charge (1.60217657E-19)
    re: Electron radius (2.8179403267E-15)
    """
    # convert to SI units
    wl = wavelength*1.e-9
    r = diameter*1.e-9/2.
    i = intensity*1e-3
    p = pixelsize*1.e-6
    D = detector_distance*1.e-3
    # Get some physical constants
    h =  spimage.DICT_physical_constants['h']
    c =  spimage.DICT_physical_constants['c']
    qe = spimage.DICT_physical_constants['e']
    re = spimage.DICT_physical_constants["re"]            
    # compute scale
    ey_J = h*c/wl
    V  = 4/3.*numpy.pi*r**3
    I_0 = i/ey_J*1.E12            
    rho_e = spimage.Material(material_type=material).get_electron_density()
    QE = detector_efficiency
    ADUP = detector_adu_photon
    K = I_0*(rho_e*p/D*re*V)**2
    scale = K * QE * ADUP
    return scale
