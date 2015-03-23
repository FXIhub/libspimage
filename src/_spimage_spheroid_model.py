import numpy
import logging
import scipy.signal
from scipy.optimize import leastsq
import scipy.stats
import spimage
from _spimage_sphere_model import _prepare_for_fitting,sphere_model_convert_intensity_to_scaling,sphere_model_convert_scaling_to_intensity,sphere_model_convert_diameter_to_size,sphere_model_convert_size_to_diameter
__all__ = ['fit_spheroid_shape','fit_spheroid_shape_pearson','fit_spheroid_intensity','fit_spheroid_intensity_pixelwise','I_spheroid_diffraction','fit_full_spheroid_model']

def fit_spheroid_shape(img, msk, diameter_a, diameter_c, phi, intensity, wavelength, pixel_size, detector_distance, method=None, full_output=False, **kwargs):
    """
    Fit the shape of a spheroid to diffraction data.

    usage:
    ======
    diameter_a, diamter_c, phi = fit_spheroid_shape(img, msk)
    diameter_a, diamter_c, phi = fit_spheroid_shape(img, msk, method='pearson', ...)

    """
    if method is None: method = 'pearson'
    if method == 'pearson': diameter_a, diameter_c, phi, info = fit_spheroid_shape_pearson(img, msk, diameter_a, diameter_c, phi, intensity, wavelength, pixel_size, detector_distance, full_output=True, **kwargs)
    else: diameter_a, diameter_c, phi, info = [diameter_a, diameter_c, phi, "There is no fitting dimensions method %s" %method]
    
    diameter_a = abs(diameter_a)
    diameter_c = abs(diameter_c)
    
    if full_output: return diameter_a, diameter_c, phi, info
    else: return diameter_a, diameter_c, phi, info

def fit_spheroid_shape_pearson(img, msk, diameter_a, diameter_c, phi, intensity, wavelength, pixel_size, detector_distance, full_output=False,  x0=0, y0=0, detector_adu_photon=1, detector_quantum_efficiency=1, material='water', rmax=None, downsampling=1, do_brute_evals=0, maxfev=1000, do_photon_counting=False, deltab=0.5, force_prolate=True):
    """
    Fit the shape of a spheroid using pearson correlation.

    """
    Xmc, Ymc, img, msk = _prepare_for_fitting(img, msk, x0, y0, rmax, downsampling, detector_adu_photon, do_photon_counting, pixel_size, detector_distance)
    diameter = (diameter_a**2*diameter_c)**(1/3.)
    S   = sphere_model_convert_intensity_to_scaling(intensity, diameter, wavelength, pixel_size, detector_distance, detector_quantum_efficiency, 1, material)
    size    = lambda d: sphere_model_convert_diameter_to_size(d, wavelength, pixel_size, detector_distance)
    I_fit_m = lambda da,dc,p: I_spheroid_diffraction(S,Xmc,Ymc,size(da),size(dc),0.,p)
    def E_fit_m(p):
        if not (img[msk].std() and I_fit_m(p[0],p[1],p[2]).std()): return numpy.ones(len(p))
        else: return 1-scipy.stats.pearsonr(I_fit_m(p[0],p[1],p[2]),img[msk])[0]*numpy.ones(len(p))
    
    # Start with brute force with a sensible range
    # We'll assume at least 20x oversampling
    if do_brute_evals:
        dmin = sphere_model_convert_size_to_diameter(1./(downsampling*img.shape[0]), wavelength, pixel_size, detector_distance)
        dmax = dmin*downsampling*img.shape[0]/20
        pmin = 0.
        pmax = numpy.pi
        Ns = do_brute_evals
        diameter = scipy.optimize.brute(E_fit_m, [(dmin, dmax),(dmin,dmax),(pmin,pmax)], Ns=Ns)[0]

    # End with least square
    p0 = numpy.array([diameter_a,diameter_c,phi])
    p, cov, infodict, mesg, ier = leastsq(E_fit_m, p0, maxfev=maxfev, xtol=1e-5, full_output=True)
    [diameter_a, diameter_c, phi] = p
    diameter_a = abs(diameter_a)
    diameter_c = abs(diameter_c)
    phi = phi % numpy.pi
    
    if force_prolate and diameter_a > diameter_c:
        return fit_spheroid_shape_pearson(img, msk, diameter_c, diameter_a, phi+numpy.pi/2., intensity, wavelength, pixel_size, detector_distance, full_output,
                                          x0, y0, detector_adu_photon, detector_quantum_efficiency, material, rmax, downsampling, do_brute_evals, maxfev, do_photon_counting, deltab)
    
    # Reduced Chi-squared and standard error
    chisquared = ((I_fit_m(diameter_a, diameter_c, phi) - img[msk])**2).sum()/(img.shape[0]*img.shape[1] - 1)
    if cov is not None:
        pcov = cov[0,0]*chisquared
    else:
        pcov = None
    
    if full_output:
        infodict['error'] = chisquared
        infodict['pcov'] = pcov
        return diameter_a, diameter_c, phi, infodict
    else:
        return diameter_a, diameter_c, phi

def fit_spheroid_intensity(img, msk, diameter_a, diameter_c, phi, intensity, wavelength, pixel_size, detector_distance, method=None, full_output=False, **kwargs):
    """
    Fit the the intensity of a spheroid to diffraction data.

    usage:
    ======
    intensity = fit_spheroid_intensity(img, msk, diameter_a, diameter_c, phi, intensity, wavelength, pixel_size, detector_size)
    intensity = fit_spheroid_intensity(img, msk, diameter_a, diameter_c, phi, intensity, wavelength, pixel_size, detector_size, method='pixelwise', ...)

    """
    if method is None: method = 'pixelwise'
    if method == 'pixelwise': intensity, info = fit_spheroid_intensity_pixelwise(img, msk, diameter_a, diameter_c, phi, intensity, wavelength, pixel_size, detector_distance, full_output=True, **kwargs)
    else: intensity, info = [intensity, "There is no fitting intensity method %s" %method]

    if full_output: return intensity, info
    else: return intensity

def fit_spheroid_intensity_pixelwise(img, msk, diameter_a, diameter_c, phi, intensity, wavelength, pixel_size, detector_distance, full_output=False, x0=0, y0=0, detector_adu_photon=1, detector_quantum_efficiency=1, material='water', rmax=None, downsampling=1, maxfev=1000, do_photon_counting=False):
    """
    Fit the intensity [J / m^2] on a spheroid to diffraction data using least square optimization.
    The cost function is defined as a pixelwise comparison between model and data.

    usage:
    ======
    intensity       = fit_spheroid_intensity(img, msk, diameter_a, diameter_c, phi, intensity, pixel_size, detector_distance)
    intensity, info = fit_spheroid_intensity(img, msk, diameter_a, diameter_c, phi, intensity, pixel_size, detector_distance, full_output=True)
    intensity, info = fit_spheroid_intensity(img, msk, diameter_a, diameter_c, phi, intensity, pixel_size, detector_distance, full_output=True, x0=0, y0=0, adup=1, queff=1, material='water', rmax=None, downsampling=1, maxfev=1000)

    Parameters:
    ===========
    ...

    """
    Xmc, Ymc, img, msk = _prepare_for_fitting(img, msk, x0, y0, rmax, downsampling, detector_adu_photon, do_photon_counting, pixel_size, detector_distance)
    size_a = sphere_model_convert_diameter_to_size(diameter_a, wavelength, pixel_size, detector_distance)
    size_c = sphere_model_convert_diameter_to_size(diameter_c, wavelength, pixel_size, detector_distance)
    diameter = (diameter_a**2*diameter_c)**(1/3.)
    scaling = lambda i: sphere_model_convert_intensity_to_scaling(i, diameter, wavelength, pixel_size, detector_distance, detector_quantum_efficiency, 1, material)
    I_fit_m = lambda i: I_spheroid_diffraction(scaling(i),Xmc,Ymc,size_a,size_c,0.,phi)
    E_fit_m = lambda i: I_fit_m(i) - img[msk]
    #print E_fit_m(intensity)

    if len(img[msk]):
        [intensity], cov, infodict, mesg, ier = leastsq(E_fit_m, intensity, maxfev=maxfev, xtol=1e-3, full_output=True)
        # Reduced Chi-squared and standard error
        chisquared = (E_fit_m(intensity)**2).sum()/(img.shape[0]*img.shape[1] - 1)
        if cov is not None:
            pcov = cov[0,0]*chisquared
        else:
            pcov = None
    else:
        infodict={}
        pcov = None
        chisquared = 0
        
    if full_output:
        infodict['error'] = chisquared
        infodict['pcov'] = pcov
        return intensity, infodict
    else:
        return intensity

def fit_full_spheroid_model(img, msk, diameter_a, diameter_c, phi, intensity, wavelength, pixel_size, detector_distance, full_output=False, x0=0, y0=0, detector_adu_photon=1., detector_quantum_efficiency=1., material='water', rmax=None, downsampling=1, maxfev=1000, deltab=0.2, do_photon_counting=False,n=1):
    diameter_a = max(diameter_a, 1.E-9)
    diameter_c = max(diameter_c, 1.E-9)
    #intensity = max(intensity,1.)
    #x0 = min(x0, img.shape[1])
    #y0 = min(y0, img.shape[0])
    for i in range(n):
        Xm, Ym, img, msk = _prepare_for_fitting(img, msk, x0, y0, rmax, downsampling, detector_adu_photon, do_photon_counting, pixel_size, detector_distance)
        size = lambda d: sphere_model_convert_diameter_to_size(d, wavelength, pixel_size, detector_distance)
        d = lambda da,dc: (da**2*dc)**(1/3.)
        scaling = lambda i,da,dc: sphere_model_convert_intensity_to_scaling(i, d(da,dc), wavelength, pixel_size, detector_distance, detector_quantum_efficiency, 1, material)
        I_fit_m = lambda dx,dy,da,dc,p,i: I_spheroid_diffraction(scaling(i,da,dc), Xm-dx, Ym-dy, size(da), size(dc), 0., p)
        E_fit_m = lambda p: I_fit_m(p[0],p[1],p[2],p[3],p[4],p[5]) - img[msk]
        p0 = numpy.array([0,0,diameter_a,diameter_c,phi,intensity])
        x0_bound = (None, None)
        y0_bound = (None, None)
        d_a_bound  = ((1-deltab)*diameter_a, (1+deltab)*diameter_a)
        d_c_bound  = ((1-deltab)*diameter_c, (1+deltab)*diameter_c)
        p_bound = (None, None)
        i_bound  = (None, None)
        bounds   = numpy.array([x0_bound, y0_bound , d_a_bound, d_c_bound, p_bound, i_bound])
        p, cov, info, mesg, ier = spimage.leastsqbound(E_fit_m, p0, maxfev=maxfev, xtol=1e-5, full_output=True, bounds=bounds)
        [dx0, dy0, diameter_a, diameter_c, phi, intensity] = p
        diameter_a = abs(diameter_a)
        diameter_c = abs(diameter_c)
        x0 += dx0
        y0 += dy0
        phi = phi % numpy.pi
        if (numpy.isfinite(p) == False).sum() > 0:
            break

    # Reduced Chi-squared and standard errors
    chisquared = (E_fit_m(p)**2).sum()/(img.shape[0]*img.shape[1] - len(p))
    #print cov
    if cov is not None:
        pcov = numpy.diag(cov)*chisquared
    else:
        pcov = numpy.array(len(p)*[None])
    if full_output:
        info["error"] = chisquared
        info["pcov"]  = pcov
        return x0, y0, diameter_a, diameter_c, phi, intensity, info
    else:
        return x0, y0, diameter_a, diameter_c, phi, intensity

def I_spheroid_diffraction(A,qX,qY,a,c,theta,phi):
    """
    scattering amplitude from homogeneous spheroid:
    -----------------------------------------------
    Sources:
    Feigin 1987
    Hamzeh,Bragg 1974
    
    f(q,a,c,theta,phi)   = 3 { sin(2pi*|q|*H) - 2pi*q*H * cos(2pi*|q|*H) } / (2pi*q*H)^3
    H(q,a,c,theta,phi)   = sqrt(a^2 sin^2(g)+c^2 cos^2(g))
    g(q,a,c,theta,phi)   = arccos( 2pi * ( q[0] cos(theta) cos(phi) - q[1] cos(theta) sin(phi) ) / (2pi*|q|) )
    I(A,q,a,c,theta,phi) = A * [f]^2 

    Parameters:
    ===========
    q: reciprocal coordinates in [px]
    a: radius perpendicular to the rotation axis of the ellipsoid [1/px]
    c: radius along the rotation axis of the ellipsoid [1/px]
    theta: rotation around x-axis (1st)
    phi: rotation around z-axis (2nd)
    A: Scaling factor in [ADU]
    """
    cos = numpy.cos
    sin = numpy.sin
    arccos = numpy.arccos
    sqrt = numpy.sqrt
    pi = numpy.pi
    eps = numpy.finfo("float64").eps
    res = numpy.finfo("float64").resolution
    #q = lambda qX,qY: sqrt(qX**2+qY**2)
    #g = lambda qX,qY,theta,phi: arccos(2*pi*(qY*cos(theta)*cos(phi)-qX*cos(theta)*sin(phi))/(2*pi*q(qX,qY)+eps))
    #H = lambda qX,qY,a,c,theta,phi: sqrt(a**2*sin(g(qX,qY,theta,phi))**2+c**2*cos(g(qX,qY,theta,phi))**2)
    #qH = lambda qX,qY,a,c,theta,phi: q(qX,qY)*H(qX,qY,a,c,theta,phi)
    #_I_spheroid_diffraction = lambda A,qX,qY,a,c,theta,phi: abs(A) * ( 3.*(sin(2*pi*qH(qX,qY,a,c,theta,phi))-2*pi*qH(qX,qY,a,c,theta,phi)*cos(2*pi*qH(qX,qY,a,c,theta,phi)))/(2*pi*qH(qX,qY,a,c,theta,phi)**3+eps) )**2
    #return ((2*pi*qH(qX,qY,a,c,theta,phi)**6 < res) * abs(A) + (2*pi*qH(qX,qY,a,c,theta,phi)**6 >= res) * _I_spheroid_diffraction(A,qX,qY,a,c,theta,phi))
    q = lambda qX,qY: sqrt(qX**2+qY**2)
    g = lambda qX,qY,theta,phi: arccos((-qX*cos(theta)*sin(phi)+qY*cos(theta)*cos(phi))/(q(qX,qY)+eps))
    H = lambda qX,qY,a,c,theta,phi: sqrt(a**2*sin(g(qX,qY,theta,phi))**2+c**2*cos(g(qX,qY,theta,phi))**2)
    qH = lambda qX,qY,a,c,theta,phi: q(qX,qY)*H(qX,qY,a,c,theta,phi)
    _I = lambda A,qX,qY,a,c,theta,phi: abs(A)*(3.*(sin(qH(qX,qY,a,c,theta,phi))-qH(qX,qY,a,c,theta,phi)*cos(qH(qX,qY,a,c,theta,phi)))/(qH(qX,qY,a,c,theta,phi)**3+eps))**2
    I = lambda A,qX,qY,a,c,theta,phi: (qH(qX,qY,a,c,theta,phi)**6 < res)*abs(A) + (qH(qX,qY,a,c,theta,phi)**6 >= res)*_I(A,qX,qY,a,c,theta,phi)
    return I(A,2*pi*qX,2*pi*qY,a,c,theta,phi)

#_convert_to_non_flat = lambda diameter_a,diameter_c,phi: [diameter_a,diameter_c,phi] if (diameter_c >= diameter_a) else [diameter_c,diameter_a,phi+numpy.pi/2.]
