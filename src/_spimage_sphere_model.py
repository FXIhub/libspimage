import numpy
import logging
import scipy.signal
from scipy.optimize import leastsq
import scipy.stats
import spimage
#from pylab import *
__all__ = ['fit_sphere_diameter', 'fit_sphere_intensity', 'fit_full_sphere_model', 'I_sphere_diffraction', 'sphere_model_convert_diameter_to_size', 'sphere_model_convert_intensity_to_scaling', 'sphere_model_convert_scaling_to_intensity']

def fit_sphere_diameter(img, msk, diameter_nm, intensity_mJ_per_um2, wavelength_nm, pixelsize_um, detector_distance_mm, method=None, full_output=False, **kwargs):
    """
    Fit the diameter of a sphere to diffraction data.

    usage:
    ======
    diameter_nm = fit_sphere_diameter(img, msk)
    diameter_nm = fit_sphere_diameter(img, msk, method='pearson', ...)

    """
    if method is None: method = 'pearson'
    if method == 'pearson': diameter_nm, info = fit_sphere_diameter_pearson(img, msk, diameter_nm, intensity_mJ_per_um2, wavelength_nm, pixelsize_um, detector_distance_mm, full_output=True, **kwargs)
    elif method == 'pixelwise': diameter_nm, info = fit_sphere_diameter_pixelwise(img, msk, diameter_nm, intensity_mJ_per_um2, wavelength_nm, pixelsize_um, detector_distance_mm, full_output=True, **kwargs)
    else: diameter_nm, info = [diameter_nm, "There is no fitting diameter method %s" %method]

    if full_output: return diameter_nm, info
    else: return diameter_nm

def fit_sphere_diameter_pearson(img, msk, diameter_nm, intensity_mJ_per_um2, wavelength_nm, pixelsize_um, detector_distance_mm, full_output=False,  x0=0, y0=0, detector_adu_photon=1, detector_quantum_efficiency=1, material='water', rmax=None, downsampling=1, do_brute_evals=0, maxfev=1000, do_photon_counting=False,**kwargs):
    """
    Fit the diameter of a sphere using pearson correlation.

    """
    Xmc, Ymc, img, msk = _prepare_for_fitting(img, msk, x0, y0, rmax, downsampling, detector_adu_photon, do_photon_counting)
    Rmc = numpy.sqrt(Xmc**2 + Ymc**2)
    S   = sphere_model_convert_intensity_to_scaling(intensity_mJ_per_um2, diameter_nm, wavelength_nm, pixelsize_um, detector_distance_mm, detector_quantum_efficiency, 1, material)
    I_fit_m = lambda d: I_sphere_diffraction(S,Rmc,sphere_model_convert_diameter_to_size(d, wavelength_nm, pixelsize_um, detector_distance_mm))
    def E_fit_m(d):
        if not (img[msk].std() and I_fit_m(d).std()): return 1.
        else: return 1-scipy.stats.pearsonr(I_fit_m(d),img[msk])[0]
        
    # Start with brute force with a sensible range
    # We'll assume at least 20x oversampling
    if do_brute_evals:
        dmin = sphere_model_convert_size_to_diameter(1./(downsampling*img.shape[0]), wavelength_nm, pixelsize_um, detector_distance_mm)
        dmax = dmin*downsampling*img.shape[0]/20
        Ns = do_brute_evals
        diameter_nm = scipy.optimize.brute(E_fit_m, [(dmin, dmax)], Ns=Ns)[0]

    # End with least square
    [diameter_nm], cov, infodict, mesg, ier = leastsq(E_fit_m, diameter_nm, maxfev=maxfev, xtol=1e-5, full_output=True)
    diameter_nm = abs(diameter_nm)

    # Reduced Chi-squared and standard error
    chisquared = ((I_fit_m(diameter_nm) - img[msk])**2).sum()/(img.shape[0]*img.shape[1] - 1)
    if cov is not None:
        pcov = cov[0,0]*chisquared
    else:
        pcov = None
    
    if full_output:
        infodict['error'] = chisquared
        infodict['pcov'] = pcov
        return diameter_nm, infodict
    else:
        return diameter_nm

def fit_sphere_diameter_pixelwise(img, msk, diameter_nm, intensity_mJ_per_um2, wavelength_nm, pixelsize_um, detector_distance_mm, full_output=False, x0=0, y0=0, detector_adu_photon=1, detector_quantum_efficiency=1, material='water', rmax=None, downsampling=1, maxfev=1000, deltab=0.5, do_photon_counting=False):
    """
    Fit the diameter of a sphere minimizing sum of pixelwise difference.
    """
    Xmc, Ymc, img, msk, = _prepare_for_fitting(img, msk, x0, y0, rmax, downsampling, detector_adu_photon, do_photon_counting)
    Rmc = numpy.sqrt(Xmc**2 + Ymc**2)
    S   = sphere_model_convert_intensity_to_scaling(intensity_mJ_per_um2, diameter_nm, wavelength_nm, pixelsize_um, detector_distance_mm, detector_quantum_efficiency, 1, material)
    I_fit_m = lambda d: I_sphere_diffraction(S,Rmc,sphere_model_convert_diameter_to_size(d, wavelength_nm, pixelsize_um, detector_distance_mm))
    E_fit_m = lambda d: I_fit_m(d) - img[msk]
    
<<<<<<< HEAD
    bounds  = numpy.array([(diameter-deltab*diameter, diameter+deltab*diameter)])
    p, cov, infodict, mesg, ier = leastsq(E_fit_m, numpy.array([diameter]), maxfev=maxfev, xtol=1e-5, full_output=True)
    #p, cov, infodict, mesg, ier = spimage.leastsqbound(E_fit_m, numpy.array([diameter]), maxfev=maxfev, xtol=1e-5, full_output=True, bounds=bounds)
    [diameter] = p
=======
    bounds  = numpy.array([(diameter_nm-deltab*diameter_nm, diameter_nm+deltab*diameter_nm)])
    p, cov, infodict, mesg, ier = leastsq(E_fit_m, numpy.array([diameter_nm]), maxfev=maxfev, xtol=1e-5, full_output=True)
    #p, cov, infodict, mesg, ier = spimage.leastsqbound(E_fit_m, numpy.array([diameter_nm]), maxfev=maxfev, xtol=1e-5, full_output=True, bounds=bounds)
    [diameter_nm] = p
>>>>>>> e0e5ecbd33a0e68211178e35bfbde973d79a5437
    
    # Reduced Chi-squared and standard error
    chisquared = ((I_fit_m(diameter_nm) - img[msk])**2).sum()/(img.shape[0]*img.shape[1] - 1)
    if cov is not None:
        pcov = cov[0,0]*chisquared
    else:
        pcov = None
    #print pcov
    if full_output:
        infodict['error'] = chisquared
        infodict['pcov']  = pcov
        return diameter_nm, infodict
    else:
        return diameter_nm
    
def fit_sphere_intensity(img, msk, diameter_nm, intensity_mJ_per_um2, wavelength_nm, pixelsize_um, detector_distance_mm, method=None, full_output=False, **kwargs):
    """
    Fit the the intensity of a sphere to diffraction data.

    usage:
    ======
    intensity_mJ_per_um2 = fit_sphere_intensity(img, msk, diameter_nm, intensity_mJ_per_um2, wavelength_nm, pixelsize_um, detector_size)
    intensity_mJ_per_um2 = fit_sphere_intensity(img, msk, diameter_nm, intensity_mJ_per_um2, wavelength_nm, pixelsize_um, detector_size, method='pixelwise', ...)
    intensity_mJ_per_um2 = fit_sphere_intensity(img, msk, diameter_nm, intensity_mJ_per_um2, wavelength_nm, pixelsize_um, detector_size, method='nrphotons', ...)

    """
    if method is None: method = 'pixelwise'
    if method == 'pixelwise': intensity_mJ_per_um2, info = fit_sphere_intensity_pixelwise(img, msk, diameter_nm, intensity_mJ_per_um2, wavelength_nm, pixelsize_um, detector_distance_mm, full_output=True, **kwargs)
    elif method == 'nrphotons': intensity_mJ_per_um2, info = fit_sphere_intensity_nrphotons(img, msk, diameter_nm, intensity_mJ_per_um2, wavelength_nm, pixelsize_um, detector_distance_mm, full_output=True, **kwargs)
    else: intensity_mJ_per_um2, info = [intensity_mJ_per_um2, "There is no fitting intensity method %s" %method]

    if full_output: return intensity_mJ_per_um2, info
    else: return intensity_mJ_per_um2
        
def fit_sphere_intensity_pixelwise(img, msk, diameter_nm, intensity_mJ_per_um2, wavelength_nm, pixelsize_um, detector_distance_mm, full_output=False, x0=0, y0=0, detector_adu_photon=1, detector_quantum_efficiency=1, material='water', rmax=None, downsampling=1, maxfev=1000, do_photon_counting=False,**kwargs):
    """
    Fit the intensity [mJ / um^2] on a sphere to diffraction data using least square optimization.
    The cost function is defined as a pixelwise comparison between model and data.

    usage:
    ======
    intensity_mJ_per_um2       = fit_sphere_intensity(img, msk, diameter_nm, intensity_mJ_per_um2, pixelsize_um, detector_distance_mm)
    intensity_mJ_per_um2, info = fit_sphere_intensity(img, msk, diameter_nm, intensity_mJ_per_um2, pixelsize_um, detector_distance_mm, full_output=True)
    intensity_mJ_per_um2, info = fit_sphere_intensity(img, msk, diameter_nm, intensity_mJ_per_um2, pixelsize_um, detector_distance_mm, full_output=True, x0=0, y0=0, adup=1, queff=1, material='water', rmax=None, downsampling=1, maxfev=1000)

    Parameters:
    ===========
    ...

    """
    Xmc, Ymc, img, msk = _prepare_for_fitting(img, msk, x0, y0, rmax, downsampling, detector_adu_photon, do_photon_counting)
    Rmc = numpy.sqrt(Xmc**2 + Ymc**2)
    size = sphere_model_convert_diameter_to_size(diameter_nm, wavelength_nm, pixelsize_um, detector_distance_mm)
    I_fit_m = lambda i: I_sphere_diffraction(sphere_model_convert_intensity_to_scaling(i, diameter_nm, wavelength_nm, pixelsize_um, detector_distance_mm, detector_quantum_efficiency, 1, material),Rmc,size)
    E_fit_m = lambda i: I_fit_m(i) - img[msk]
    #print E_fit_m(intensity_mJ_per_um2)

    if len(img[msk]):
        [intensity_mJ_per_um2], cov, infodict, mesg, ier = leastsq(E_fit_m, intensity_mJ_per_um2, maxfev=maxfev, xtol=1e-3, full_output=True)
        # Reduced Chi-squared and standard error
        chisquared = (E_fit_m(intensity_mJ_per_um2)**2).sum()/(img.shape[0]*img.shape[1] - 1)
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
        return intensity_mJ_per_um2, infodict
    else:
        return intensity_mJ_per_um2

def fit_sphere_intensity_nrphotons(img, msk, diameter_nm, intensity_mJ_per_um2, wavelength_nm, pixelsize_um, detector_distance_mm, full_output=False, x0=0, y0=0, detector_adu_photon=1, detector_quantum_efficiency=1, material='water', rmax=None, downsampling=1, maxfev=1000, do_photon_counting=False,**kwargs):
    """
    Fit the intensity [mJ / um^2] on a sphere to diffraction data using least square optimization.
    The cost function is defined as the difference between the nr. of estimated photons scattered in both model and data.
    """
    Xmc, Ymc, img, msk = _prepare_for_fitting(img, msk, x0, y0, rmax, downsampling, detector_adu_photon, do_photon_counting)
    Rmc = numpy.sqrt(Xmc**2 + Ymc**2)
    nr_photons = img[msk].sum()
    size = sphere_model_convert_diameter_to_size(diameter_nm, wavelength_nm, pixelsize_um, detector_distance_mm)
    I_fit_m = lambda i: I_sphere_diffraction(sphere_model_convert_intensity_to_scaling(i, diameter_nm, wavelength_nm, pixelsize_um, detector_distance_mm, detector_quantum_efficiency, 1, material),Rmc,size)
    E_fit_m = lambda i: I_fit_m(i).sum() - nr_photons
    #print E_fit_m(intensity_mJ_per_um2)
    intensity_mJ_per_um2 = scipy.optimize.newton(E_fit_m, intensity_mJ_per_um2, maxiter=maxfev, tol=1e-3)
    [intensity_mJ_per_um2], cov, infodict, mesg, ier = leastsq(E_fit_m, intensity_mJ_per_um2, maxfev=maxfev, xtol=1e-3, full_output=True)
    
    # Reduced Chi-squared and standard error
    chisquared = (E_fit_m(intensity_mJ_per_um2)**2).sum()/(img.shape[0]*img.shape[1] - 1)
    if cov is not None:
        pcov = cov[0,0]*chisquared
    else:
        pcov = None
    
    if full_output:
        infodict['error'] = chisquared
        infodict['pcov'] = pcov
        return intensity_mJ_per_um2, infodict
    else:
        return intensity_mJ_per_um2
    
def fit_full_sphere_model(img, msk, diameter_nm, intensity_mJ_per_um2, wavelength_nm, pixelsize_um, detector_distance_mm, full_output=False, x0=0, y0=0, detector_adu_photon=1., detector_quantum_efficiency=1., material='water', rmax=None, downsampling=1, maxfev=1000, deltab=0.2, do_photon_counting=False):
    diameter_nm = max(diameter_nm, 1.)
    #intensity_mJ_per_um2 = max(intensity_mJ_per_um2,1.)
    #x0 = min(x0, img.shape[1])
    #y0 = min(y0, img.shape[0])
    Xm, Ym, img, msk = _prepare_for_fitting(img, msk, 0, 0, rmax, downsampling, detector_adu_photon, do_photon_counting)
    Rmc     = lambda x,y:     numpy.sqrt((Xm - x)**2 + (Ym - y)**2)
    size    = lambda d:       sphere_model_convert_diameter_to_size(d, wavelength_nm, pixelsize_um, detector_distance_mm)
    scaling = lambda i,d:     sphere_model_convert_intensity_to_scaling(i, d, wavelength_nm, pixelsize_um, detector_distance_mm, detector_quantum_efficiency, 1, material)
    I_fit_m = lambda x,y,d,i: I_sphere_diffraction(scaling(i,d), Rmc(x,y), size(d))
    E_fit_m = lambda p:       I_fit_m(p[0],p[1],p[2],p[3]) - img[msk]
    p0 = numpy.array([x0,y0,diameter_nm,intensity_mJ_per_um2])
    #print E_fit_m(p0)
    x0_bound = (None, None)
    y0_bound = (None, None)
    d_bound  = (diameter_nm-deltab*diameter_nm, diameter_nm+deltab*diameter_nm)
    i_bound = (None, None)
    bounds   = numpy.array([x0_bound, y0_bound , d_bound, i_bound])
<<<<<<< HEAD
    p, cov, infodict, mesg, ier = spimage.leastsqbound(E_fit_m, numpy.array([x0,y0,diameter,intensity]), maxfev=maxfev, xtol=1e-5, full_output=True, bounds=bounds)
    [x0, y0, diameter, intensity] = p
=======
    p, cov, infodict, mesg, ier = spimage.leastsqbound(E_fit_m, numpy.array([x0,y0,diameter_nm,intensity_mJ_per_um2]), maxfev=maxfev, xtol=1e-5, full_output=True, bounds=bounds)
    [x0, y0, diameter_nm, intensity_mJ_per_um2] = p
>>>>>>> e0e5ecbd33a0e68211178e35bfbde973d79a5437

    # Reduced Chi-squared and standard errors
    chisquared = (E_fit_m(p)**2).sum()/(img.shape[0]*img.shape[1] - len(p))
    if cov is not None:
        pcov = numpy.diag(cov)*chisquared
    else:
        pcov = numpy.array(4*[None])
    if full_output:
        infodict["error"] = chisquared
        infodict["pcov"]  = pcov
        return x0, y0, diameter_nm, intensity_mJ_per_um2, infodict
    else:
        return x0, y0, diameter_nm, intensity_mJ_per_um2
    
def _prepare_for_fitting(img, msk, x0, y0, rmax, downsampling, adup, do_photon_counting):
    if rmax is None: rmax = max(img.shape)/2  # Handle default of rmax
    s = img.shape # Shape of image
    Y,X = spimage.grid(s, (0,0)) # Non-centered grid vectors in [px]
    Mr = spimage.rmask(s, rmax, (x0,y0))      # Radial mask
    msk = msk * Mr   # Merge mask and radial mask
    img = img[::downsampling, ::downsampling]
    msk = msk[::downsampling, ::downsampling]
    Xm = X[msk] # Non-centered (masked) grid x vectors in [px]
    Ym = Y[msk] # Non-centered (masked) grid y vectors in [px]
    Xmc = Xm - x0   # Centered (masked) grid x vectors in [px]
    Ymc = Ym - y0   # Centered (masked) grid y vectors in [px]
    img = img / adup
    if do_photon_counting: img = numpy.round(img) * (img > 0)
    return Xmc, Ymc, img, msk

def I_sphere_diffraction(A,q,s):
    """
    scattered intensity from homogeneous sphere.
    --------------------------------------------
    Source:
    Feigin 1987

    f(q,s) = 3 { sin(2pi*q*s) - 2pi*q*s cos(2pi*q*s) } / (2pi*q*s)^3
    I(A,q,s) = A * [f(q,s)]^2

    Parameters:
    ===========
    q: reciprocal coordinates in [px]
    s: size of the sphere in [1/px]
    A: Scaling factor in [ADU]
    """
    _I_sphere_diffraction = lambda A,q,s: abs(A) * ( 3*(numpy.sin(2*numpy.pi*q*s)-2*numpy.pi*q*s*numpy.cos(2*numpy.pi*q*s))/((2*numpy.pi*q*s)**3+numpy.finfo("float64").eps) )**2
    return ((q*s)**6 < numpy.finfo("float64").resolution) * abs(A) + ((q*s)**6 >= numpy.finfo("float64").resolution) * _I_sphere_diffraction(A,q,s)

def sphere_model_convert_diameter_to_size(diameter_nm, wavelength_nm, pixelsize_um, detector_distance_mm):
    """
    converts the sphere diameter [nm] to size of the model [1/px]

    size = ( diameter * pixelsize ) / ( 2 * wavelength * detector_distance )

    Parameters:
    ===========
    diameter_nm:          The diameter of the sphere in [nm]
    wavelength_nm:        The wavelength of the xray beam in [nm]
    pixelsize_um:         The size of one detector pixel in [mu m]
    detector_distance_mm: The distance from the interaction to the detector in [mm]
    """    
    r  = 1.e-9 * diameter_nm/2.
    p  = 1.e-6 * pixelsize_um
    D  = 1.e-3 * detector_distance_mm
    wl = 1.e-9 * wavelength_nm
    return (r*p) / (D*wl)    

def sphere_model_convert_size_to_diameter(size, wavelength_nm, pixelsize_um, detector_distance_mm):
    """
    converts the size of the model [1/px] to sphere diameter [nm]

    diameter = 2 * size * detector_distance * wavelength / pixelsize

    Parameters:
    ===========
    size:              The size of the sphere model in [1/px]
    wavelength_nm:        The wavelength of the xray beam in [nm]
    pixelsize_um:         The size of one detector pixel in [mu m]
    detector_distance_mm: The distance from the interaction to the detector in [mm]
    """    
    p  = 1.e-6 * pixelsize_um
    D  = 1.e-3 * detector_distance_mm
    wl = 1.e-9 * wavelength_nm
    r  = (size*D*wl) / p
    diameter_nm = 1.e9 * r * 2
    return diameter_nm

def sphere_model_convert_intensity_to_scaling(intensity_mJ_per_um2, diameter_nm, wavelength_nm, pixelsize_um, detector_distance_mm, detector_qe=1.0, detector_adu_photon=1.0, material="water"):
    """
    converts the photon intensity [mJ/um^2] on a sphere to the scaling factor [ADU] of the sphere model.

    V = (4/3 * pi * diameter/2 ) **3
    I0 = (intensity / (hc / wavelength) ) * 1e12
    K =  I0 * ( REAL{dn} * pixelsize / detector_distance 2 pi / wavelength**2 * V )**2)
    ( less accurate: K = I0 * ( rho_e * pixelsize / detector_distance * re * V )**2 )
    scaling = K * detector_adu_photon * detector_efficiency

    Parameters:
    ===========
    intensity_mJ_per_um2:    The photon intensity on the sample in [mJ/um^2]
    diameter_nm:             The diameter of the sphere in [nm]
    wavelength_nm:           The wavelength of the xray beam in [nm]
    pixelsize_um:            The size of one detector pixel in [um/px]
    detector_distance_mm:    The distance from the interaction to the detector in [mm]
    detector_efficiency:     The quantum efficieny of the detector (default: 1.0)
    detector_adu_per_photon: The nr. ADU units per photon [ADU / Np] (default: 1.0)
    material:                The material of the sample (default: water) -> defines electron_density rho_e in [px/m^3]

    Physical constants used:
    ========================
    h:  Planck constant (6.62606957E-34)   in [Js]
    c:  Speed of light  (299792458)        in [m/s]
    qe: Electron charge (1.60217657E-19)   in [C]
    re: Electron radius (2.8179403267E-15) in [m]
    """
    # Convert to SI units
    wl = 1.e-9*wavelength_nm        # [m]
    r  = 1.e-9*diameter_nm/2.       # [m]
    i  = 1.e-3*intensity_mJ_per_um2 # [J/um^2]
    p  = 1.e-6*pixelsize_um         # [m/px]
    D  = 1.e-3*detector_distance_mm # [m]
    # Get some physical constants
    h  = spimage.DICT_physical_constants['h']  # [Js]
    c  = spimage.DICT_physical_constants['c']  # [m/s]
    qe = spimage.DICT_physical_constants['e']  # [C]
    re = spimage.DICT_physical_constants["re"] # [m]
    # Convert i to A
    ey_J  = h*c/wl                 # [J] 
    V     = 4/3.*numpy.pi*r**3     # [m^3]
    I0    = i/ey_J*1e12            # [Np/m^2] 
    #rho_e = spimage.Material(material_type=material).get_electron_density() # [px/m^3]
    ph_eV = h*c/wl/qe              # [eV]
    dn    = spimage.Material(material_type=material).get_dn(photon_energy_eV=ph_eV) # []
    QE    = detector_qe
    ADUP  = detector_adu_photon
    #K     = I0*(rho_e*p/D*re*V)**2 # [Np/m^2 (px 1/m^3 m/px 1/m m m^3 )^2]
    K     = I0*(dn.real*p/D*2*numpy.pi/wl**2*V)**2 # [Np/m^2 (px 1/m^3 m/px 1/m m m^3 )^2]
    A     = K*QE*ADUP              # [ADU]  
    return A

def sphere_model_convert_scaling_to_intensity(scaling, diameter_nm, wavelength_nm, pixelsize_um, detector_distance_mm, detector_qe=1.0, detector_adu_photon=1.0, material="water"):
    """
    converts the scaling factor [ADU] of the sphere model to photon intensity [mJ/um^2] on sphere.

    V = (4/3 * pi * diameter/2 )**3 
    K = scaling / detector_adu_photon / detector_efficiency
    I0 = K / ( rho_e * pixelsize / detector_distance * re * V )**2
    intensity = (I0 * 1e-12) * (hc / wavelength)
    
    Parameters:
    ===========
    scaling:             The scaling factor of the sphere model [ADU]
    diameter_nm:            The diameter of the sphere in [nm]
    wavelength_nm:          The wavelength of the xray beam in [nm]
    pixelsize_um:           The size of one detector pixel in [um/px]
    detector_distance_mm:   The distance from the interaction to the detector in [mm]
    detector_efficiency: The quantum efficieny of the detector
    detector_adu_photon: The nr. ADU units per photon [ADU / Np]
    material:            The material of the sample -> defines electron_density rho_e [px/m^3]

    Physical constants used:
    ========================
    h:  Planck constant (6.62606957E-34)     in [Js]
    c:  Speed of light  (299792458)          in [m/s]
    qe: Electron charge (1.60217657E-19)     in [C]
    re: Electron radius (2.8179403267E-15)   in [m]
    """
    # convert to SI units
    wl = 1.e-9*wavelength_nm         # [m]
    r  = 1.e-9*diameter_nm/2.        # [m]
    p  = 1.e-6*pixelsize_um          # [m/px]
    D  = 1.e-3*detector_distance_mm  # [m]
    # Get some physical constants
    h  = spimage.DICT_physical_constants['h']  # [Js]
    c  = spimage.DICT_physical_constants['c']  # [m/s]
    qe = spimage.DICT_physical_constants['e']  # [C]
    re = spimage.DICT_physical_constants["re"] # [m]
    # Convert A to i
    ey_J  = h*c/wl                # [J]
    V     = 4/3.*numpy.pi*r**3    # [m^3]
    #rho_e = spimage.Material(material_type=material).get_electron_density() # [px/m^3]
    ph_eV = h*c/wl/qe             # [eV]
    dn    = spimage.Material(material_type=material).get_dn(photon_energy_eV=ph_eV) # []
    QE    = detector_qe   # []
    ADUP  = detector_adu_photon   # [ADU / Np]
    K     = scaling/QE/ADUP       # [Np / m^2
    #I0    = K/(rho_e*p/D*re*V)**2 # [Np / m^2]
    I0    = K/(dn.real*p/D*2*numpy.pi/wl**2*V)**2 # [Np/m^2]
    i     = I0 * 1e-12 * ey_J     # [J/m^2]
    # Convert i to [mJ/um^2]
    intensity_mJ_per_um2 = 1e3*i
    return intensity_mJ_per_um2
