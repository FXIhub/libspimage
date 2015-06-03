import numpy
import logging
import scipy.signal
from scipy.optimize import leastsq
import scipy.stats
import spimage
#from pylab import *
#__all__ = ['fit_sphere_diameter', 'fit_sphere_diameter_radial', 'fit_sphere_intensity', 'fit_full_sphere_model', 'I_sphere_diffraction', 'sphere_model_convert_diameter_to_size', 'sphere_model_convert_intensity_to_scaling', 'sphere_model_convert_scaling_to_intensity']

def fit_sphere_diameter(img, msk, diameter, intensity, wavelength, pixel_size, detector_distance, method=None, full_output=False, **kwargs):
    """
    Fit the diameter of a sphere to diffraction data.

    usage:
    ======
    diameter = fit_sphere_diameter(img, msk)
    diameter = fit_sphere_diameter(img, msk, method='pearson', ...)

    """
    if method is None: method = 'pearson'
    if method == 'pearson': diameter, info = fit_sphere_diameter_pearson(img, msk, diameter, intensity, wavelength, pixel_size, detector_distance, full_output=True, **kwargs)
    elif method == 'pixelwise': diameter, info = fit_sphere_diameter_pixelwise(img, msk, diameter, intensity, wavelength, pixel_size, detector_distance, full_output=True, **kwargs)
    else: diameter, info = [diameter, "There is no fitting diameter method %s" %method]

    diameter = abs(diameter)

    if full_output: return diameter, info
    else: return diameter

                                             
def fit_sphere_diameter_pearson(img, msk, diameter, intensity, wavelength, pixel_size, detector_distance, full_output=False,  x0=0, y0=0, detector_adu_photon=1, detector_quantum_efficiency=1, material='water', rmax=None, downsampling=1, do_brute_evals=0, maxfev=1000, do_photon_counting=False):
    """
    Fit the diameter of a sphere using pearson correlation.
    """
    Xmc, Ymc, img, msk = _prepare_for_fitting(img, msk, x0, y0, rmax, downsampling, detector_adu_photon, do_photon_counting, pixel_size, detector_distance)
    Rmc = numpy.sqrt(Xmc**2 + Ymc**2)
    S   = sphere_model_convert_intensity_to_scaling(intensity, diameter, wavelength, pixel_size, detector_distance, detector_quantum_efficiency, 1, material)
    I_fit_m = lambda d: I_sphere_diffraction(S,Rmc,sphere_model_convert_diameter_to_size(d, wavelength, pixel_size, detector_distance))
    def E_fit_m(d):
        if not (img[msk].std() and I_fit_m(d).std()): return 1.
        else: return 1-scipy.stats.pearsonr(I_fit_m(d),img[msk])[0]
        
    # Start with brute force with a sensible range
    # We'll assume at least 20x oversampling
    if do_brute_evals:
        dmin = sphere_model_convert_size_to_diameter(1./(downsampling*img.shape[0]), wavelength, pixel_size, detector_distance)
        dmax = dmin*downsampling*img.shape[0]/20
        Ns = do_brute_evals
        diameter = scipy.optimize.brute(E_fit_m, [(dmin, dmax)], Ns=Ns)[0]

    # End with least square
    [diameter], cov, infodict, mesg, ier = leastsq(E_fit_m, diameter, maxfev=maxfev, xtol=1e-5, full_output=True)
    diameter = abs(diameter)

    # Reduced Chi-squared and standard error
    chisquared = ((I_fit_m(diameter) - img[msk])**2).sum()/(img.shape[0]*img.shape[1] - 1)
    if cov is not None:
        pcov = cov[0,0]*chisquared
    else:
        pcov = numpy.nan
    
    if full_output:
        infodict['error'] = chisquared
        infodict['pcov'] = pcov
        return diameter, infodict
    else:
        return diameter

def fit_sphere_diameter_pixelwise(img, msk, diameter, intensity, wavelength, pixel_size, detector_distance, full_output=False, x0=0, y0=0, detector_adu_photon=1, detector_quantum_efficiency=1, material='water', rmax=None, downsampling=1, maxfev=1000, deltab=0.5, do_photon_counting=False):
    """
    Fit the diameter of a sphere minimizing sum of pixelwise difference.
    """
    Xmc, Ymc, img, msk, = _prepare_for_fitting(img, msk, x0, y0, rmax, downsampling, detector_adu_photon, do_photon_counting, pixel_size, detector_distance)
    Rmc = numpy.sqrt(Xmc**2 + Ymc**2)
    S   = sphere_model_convert_intensity_to_scaling(intensity, diameter, wavelength, pixel_size, detector_distance, detector_quantum_efficiency, 1, material)
    I_fit_m = lambda d: I_sphere_diffraction(S,Rmc,sphere_model_convert_diameter_to_size(d, wavelength, pixel_size, detector_distance))
    E_fit_m = lambda d: I_fit_m(d) - img[msk]
    
    bounds  = numpy.array([(diameter-deltab*diameter, diameter+deltab*diameter)])
    p, cov, infodict, mesg, ier = leastsq(E_fit_m, numpy.array([diameter]), maxfev=maxfev, xtol=1e-5, full_output=True)
    [diameter] = p
    diameter = abs(diameter)
    
    # Reduced Chi-squared and standard error
    chisquared = ((I_fit_m(diameter) - img[msk])**2).sum()/(img.shape[0]*img.shape[1] - 1)
    if cov is not None:
        pcov = cov[0,0]*chisquared
    else:
        pcov = numpy.nan

    if full_output:
        infodict['error'] = chisquared
        infodict['pcov']  = pcov
        return diameter, infodict
    else:
        return diameter
    
def fit_sphere_intensity(img, msk, diameter, intensity, wavelength, pixel_size, detector_distance, method=None, full_output=False, **kwargs):
    """
    Fit the the intensity of a sphere to diffraction data.

    usage:
    ======
    intensity = fit_sphere_intensity(img, msk, diameter, intensity, wavelength, pixel_size, detector_size)
    intensity = fit_sphere_intensity(img, msk, diameter, intensity, wavelength, pixel_size, detector_size, method='pixelwise', ...)
    intensity = fit_sphere_intensity(img, msk, diameter, intensity, wavelength, pixel_size, detector_size, method='nrphotons', ...)
    """
    if method is None: method = 'pixelwise'
    if method == 'pixelwise': intensity, info = fit_sphere_intensity_pixelwise(img, msk, diameter, intensity, wavelength, pixel_size, detector_distance, full_output=True, **kwargs)
    elif method == 'nrphotons': intensity, info = fit_sphere_intensity_nrphotons(img, msk, diameter, intensity, wavelength, pixel_size, detector_distance, full_output=True, **kwargs)
    else: intensity, info = [intensity, "There is no fitting intensity method %s" %method]

    if full_output: return intensity, info
    else: return intensity
        
def fit_sphere_intensity_pixelwise(img, msk, diameter, intensity, wavelength, pixel_size, detector_distance, full_output=False, x0=0, y0=0, detector_adu_photon=1, detector_quantum_efficiency=1, material='water', rmax=None, downsampling=1, maxfev=1000, do_photon_counting=False):
    """
    Fit the intensity [J / m^2] on a sphere to diffraction data using least square optimization.
    The cost function is defined as a pixelwise comparison between model and data.

    usage:
    ======
    intensity       = fit_sphere_intensity(img, msk, diameter, intensity, pixel_size, detector_distance)
    intensity, info = fit_sphere_intensity(img, msk, diameter, intensity, pixel_size, detector_distance, full_output=True)
    intensity, info = fit_sphere_intensity(img, msk, diameter, intensity, pixel_size, detector_distance, full_output=True, x0=0, y0=0, adup=1, queff=1, material='water', rmax=None, downsampling=1, maxfev=1000)

    Parameters:
    ===========
    ...

    """
    Xmc, Ymc, img, msk = _prepare_for_fitting(img, msk, x0, y0, rmax, downsampling, detector_adu_photon, do_photon_counting, pixel_size, detector_distance)
    Rmc = numpy.sqrt(Xmc**2 + Ymc**2)
    size = sphere_model_convert_diameter_to_size(diameter, wavelength, pixel_size, detector_distance)
    I_fit_m = lambda i: I_sphere_diffraction(sphere_model_convert_intensity_to_scaling(i, diameter, wavelength, pixel_size, detector_distance, detector_quantum_efficiency, 1, material),Rmc,size)
    E_fit_m = lambda i: I_fit_m(i) - img[msk]

    if len(img[msk]):
        [intensity], cov, infodict, mesg, ier = leastsq(E_fit_m, intensity, maxfev=maxfev, xtol=1e-3, full_output=True)
        # Reduced Chi-squared and standard error
        chisquared = (E_fit_m(intensity)**2).sum()/(img.shape[0]*img.shape[1] - 1)
        if cov is not None:
            pcov = cov[0,0]*chisquared
        else:
            pcov = np.nan
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

def fit_sphere_intensity_nrphotons(img, msk, diameter, intensity, wavelength, pixel_size, detector_distance, full_output=False, x0=0, y0=0, detector_adu_photon=1, detector_quantum_efficiency=1, material='water', rmax=None, downsampling=1, maxfev=1000, do_photon_counting=False):
    """
    Fit the intensity [J / m^2] on a sphere to diffraction data using least square optimization.
    The cost function is defined as the difference between the nr. of estimated photons scattered in both model and data.
    """
    Xmc, Ymc, img, msk = _prepare_for_fitting(img, msk, x0, y0, rmax, downsampling, detector_adu_photon, do_photon_counting, pixel_size, detector_distance)
    Rmc = numpy.sqrt(Xmc**2 + Ymc**2)
    nr_photons = img[msk].sum()
    size = sphere_model_convert_diameter_to_size(diameter, wavelength, pixel_size, detector_distance)
    I_fit_m = lambda i: I_sphere_diffraction(sphere_model_convert_intensity_to_scaling(i, diameter, wavelength, pixel_size, detector_distance, detector_quantum_efficiency, 1, material),Rmc,size)
    E_fit_m = lambda i: I_fit_m(i).sum() - nr_photons
    intensity = scipy.optimize.newton(E_fit_m, intensity, maxiter=maxfev, tol=1e-3)
    [intensity], cov, infodict, mesg, ier = leastsq(E_fit_m, intensity, maxfev=maxfev, xtol=1e-3, full_output=True)
    
    # Reduced Chi-squared and standard error
    chisquared = (E_fit_m(intensity)**2).sum()/(img.shape[0]*img.shape[1] - 1)
    if cov is not None:
        pcov = cov[0,0]*chisquared
    else:
        pcov = np.nan
    
    if full_output:
        infodict['error'] = chisquared
        infodict['pcov'] = pcov
        return intensity, infodict
    else:
        return intensity
    
def fit_full_sphere_model(img, msk, diameter, intensity, wavelength, pixel_size, detector_distance, full_output=False, x0=0, y0=0, detector_adu_photon=1., detector_quantum_efficiency=1., material='water', rmax=None, downsampling=1, maxfev=1000, deltab=0.2, do_photon_counting=False, n=1):
    diameter = max(diameter, 1.E-9)
    for i in range(n):
        Xm, Ym, _img, _msk = _prepare_for_fitting(img, msk, x0, y0, rmax, downsampling, detector_adu_photon, do_photon_counting, pixel_size, detector_distance)
        Rmc     = lambda dx,dy:     numpy.sqrt((Xm - dx)**2 + (Ym - dy)**2)
        size    = lambda d:       sphere_model_convert_diameter_to_size(d, wavelength, pixel_size, detector_distance)
        scaling = lambda i,d:     sphere_model_convert_intensity_to_scaling(i, d, wavelength, pixel_size, detector_distance, detector_quantum_efficiency, 1, material)
        I_fit_m = lambda dx,dy,d,i: I_sphere_diffraction(scaling(i,d), Rmc(dx,dy), size(d))
        E_fit_m = lambda p:       I_fit_m(p[0],p[1],p[2],p[3]) - _img[_msk]
        p0 = numpy.array([0.,0.,diameter,intensity])
        #print E_fit_m(p0)
        x0_bound = (None, None)
        y0_bound = (None, None)
        d_bound  = (diameter-deltab*diameter, diameter+deltab*diameter)
        i_bound = (None, None)
        bounds   = numpy.array([x0_bound, y0_bound , d_bound, i_bound])
        p, cov, infodict, mesg, ier = spimage.leastsqbound(E_fit_m, numpy.array([x0,y0,diameter,intensity]), maxfev=maxfev, xtol=1e-5, full_output=True, bounds=bounds)
        [dx0, dy0, diameter, intensity] = p
        diameter = abs(diameter)
        x0 += dx0
        y0 += dy0
        if (numpy.isfinite(p) == False).sum() > 0:
            break


    # Reduced Chi-squared and standard errors
    chisquared = (E_fit_m(p)**2).sum()/(_img.shape[0]*_img.shape[1] - len(p))
    if cov is not None:
        pcov = numpy.diag(cov)*chisquared
    else:
        pcov = numpy.array(4*[numpy.nan])
    if full_output:
        infodict["error"] = chisquared
        infodict["pcov"]  = pcov
        return x0, y0, diameter, intensity, infodict
    else:
        return x0, y0, diameter, intensity
    
def _prepare_for_fitting(img, msk, x0, y0, rmax, downsampling, adup, do_photon_counting, pixel_size, detector_distance):
    if rmax is None: rmax = max(img.shape)/2  # Handle default of rmax
    s = img.shape # Shape of image
    Mr = spimage.rmask(s, rmax, (x0,y0))      # Radial mask
    msk = msk * Mr   # Merge mask and radial mask
    img = img[::downsampling, ::downsampling]
    msk = msk[::downsampling, ::downsampling]
    Y,X = spimage.grid(s, (0,0)) # Non-centered grid vectors in [px]
    # Center grid
    Xc = X - x0   # Centered grid x vectors in [px]
    Yc = Y - y0   # Centered grid y vectors in [px]
    # Grid to q-space (but we keep units in px)
    Xc = spimage.x_to_qx(Xc,pixel_size,detector_distance)
    Yc = spimage.y_to_qy(Yc,pixel_size,detector_distance)
    Xmc = Xc[msk] # Centered and masked grid x vectors in [px]
    Ymc = Yc[msk] # Centered and masked grid y vectors in [px]
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

def sphere_model_convert_diameter_to_size(diameter, wavelength, pixel_size, detector_distance):
    """
    converts the sphere diameter [nm] to size of the model [1/px]

    size = ( diameter * pixel_size ) / ( 2 * wavelength * detector_distance )

    Parameters:
    ===========
    diameter:          The diameter of the sphere in [m]
    wavelength:        The wavelength of the xray beam in [m]
    pixel_size:        The size of one detector pixel in [m]
    detector_distance: The distance from the interaction to the detector in [m]
    """    
    r  = diameter/2.
    p  = pixel_size
    D  = detector_distance
    wl = wavelength
    return (r*p) / (D*wl)    

def sphere_model_convert_size_to_diameter(size, wavelength, pixel_size, detector_distance):
    """
    converts the size of the model [1/px] to sphere diameter [m]

    diameter = 2 * size * detector_distance * wavelength / pixel_size

    Parameters:
    ===========
    size:              The size of the sphere model in [1/px]
    wavelength:        The wavelength of the xray beam in [m]
    pixel_size:        The size of one detector pixel in [m]
    detector_distance: The distance from the interaction to the detector in [m]
    """    
    p  = pixel_size
    D  = detector_distance
    wl = wavelength
    r  = (size*D*wl) / p
    diameter = r * 2
    return diameter

def sphere_model_convert_intensity_to_scaling(intensity, diameter, wavelength, pixel_size, detector_distance, detector_qe=1.0, detector_adu_photon=1.0, material="water"):
    """
    converts the photon intensity [1/m2] on a sphere to the scaling factor [ADU] of the sphere model.

    V = (4/3 * pi * diameter/2 ) **3
    I0 = (intensity / (hc / wavelength) ) * 1e12
    K =  I0 * ( REAL{dn} * pixel_size / detector_distance 2 pi / wavelength**2 * V )**2)
    ( less accurate: K = I0 * ( rho_e * pixel_size / detector_distance * re * V )**2 )
    scaling = K * detector_adu_photon * detector_efficiency

    Parameters:
    ===========
    intensity:               The photon intensity on the sample in [J/m2]
    diameter:                The diameter of the sphere in [m]
    wavelength:              The wavelength of the xray beam in [m]
    pixel_size:              The size of one detector pixel in [m/px]
    detector_distance:       The distance from the interaction to the detector in [m]
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
    wl = wavelength        # [m]
    r  = diameter/2.       # [m]
    i  = intensity         # [J/m^2]
    p  = pixel_size        # [m/px]
    D  = detector_distance # [m]
    # Get some physical constants
    h  = spimage.DICT_physical_constants['h']  # [Js]
    c  = spimage.DICT_physical_constants['c']  # [m/s]
    qe = spimage.DICT_physical_constants['e']  # [C]
    re = spimage.DICT_physical_constants["re"] # [m]
    # Convert i to A
    ey_J  = h*c/wl                 # [J] 
    V     = 4/3.*numpy.pi*r**3     # [m^3]
    I0    = i/ey_J            # [J/m^2] 
    #rho_e = spimage.Material(material_type=material).get_electron_density() # [px/m^3]
    ph_eV = h*c/wl/qe              # [eV]
    dn    = spimage.Material(material_type=material).get_dn(photon_energy_eV=ph_eV) # []
    QE    = detector_qe
    ADUP  = detector_adu_photon
    #K     = I0*(rho_e*p/D*re*V)**2 # [1/m^2 (px 1/m^3 m/px 1/m m m^3 )^2]
    K     = I0*(dn.real*p/D*2*numpy.pi/wl**2*V)**2 # [Np/m^2 (px 1/m^3 m/px 1/m m m^3 )^2]
    A     = K*QE*ADUP              # [ADU]  
    return A

def sphere_model_convert_scaling_to_intensity(scaling, diameter, wavelength, pixel_size, detector_distance, detector_qe=1.0, detector_adu_photon=1.0, material="water"):
    """
    converts the scaling factor [ADU] of the sphere model to photon intensity [J/m^2] on sphere.

    V = (4/3 * pi * diameter/2 )**3 
    K = scaling / detector_adu_photon / detector_efficiency
    I0 = K / ( rho_e * pixel_size / detector_distance * re * V )**2
    intensity = (I0 * 1e-12) * (hc / wavelength)
    
    Parameters:
    ===========
    scaling:             The scaling factor of the sphere model [ADU]
    diameter:            The diameter of the sphere in [m]
    wavelength:          The wavelength of the xray beam in [m]
    pixel_size:          The size of one detector pixel in [m/px]
    detector_distance:   The distance from the interaction to the detector in [m]
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
    wl = wavelength         # [m]
    r  = diameter/2.        # [m]
    p  = pixel_size         # [m/px]
    D  = detector_distance  # [m]
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
    QE    = detector_qe           # []
    ADUP  = detector_adu_photon   # [ADU / Np]
    K     = scaling/QE/ADUP       # [Np / m^2]
    #I0    = K/(rho_e*p/D*re*V)**2 # [Np / m^2]
    I0    = K/(dn.real*p/D*2*numpy.pi/wl**2*V)**2 # [Np/m^2]
    i     = I0 * ey_J             # [J/m^2]
    return i


def fit_sphere_diameter_radial(r, img_r, msk_r, diameter, intensity, wavelength, pixel_size, detector_distance, full_output=False, detector_adu_photon=1, detector_quantum_efficiency=1, material='water', maxfev=1000, do_brute_evals=0):
    if len(img_r.shape) > 1 or len(msk_r.shape) > 1:
        print "Error: Inputs have to be one-dimensional."
        return
    S = sphere_model_convert_intensity_to_scaling(intensity, diameter, wavelength, pixel_size, detector_distance, detector_quantum_efficiency, 1, material)
    Rm = numpy.array(r, dtype="float")
    Rm = spimage.x_to_qx(Rm, pixel_size, detector_distance)
    Rm = Rm[msk_r]
    I_fit_m = lambda d: I_sphere_diffraction(S,Rm,sphere_model_convert_diameter_to_size(d, wavelength, pixel_size, detector_distance))
    def E_fit_m(d):
        if not (img_r[msk_r].std() and I_fit_m(d).std()): return 1.
        else: return 1-scipy.stats.pearsonr(I_fit_m(d), img_r[msk_r])[0]

    # Start with brute force with a sensible range
    # We'll assume at least 20x oversampling
    if do_brute_evals:
        dmin = sphere_model_convert_size_to_diameter(1./(img_r.size), wavelength, pixel_size, detector_distance)
        dmax = dmin*img.size/20
        Ns = do_brute_evals
        diameter = scipy.optimize.brute(E_fit_m, [(dmin, dmax)], Ns=Ns)[0]

    # End with least square
    [diameter], cov, infodict, mesg, ier = leastsq(E_fit_m, diameter, maxfev=maxfev, xtol=1e-5, full_output=True)
    diameter = abs(diameter)

    # Reduced Chi-squared and standard error
    chisquared = ((I_fit_m(diameter) - img_r[msk_r])**2).sum()/(img_r.size - 1)
    if cov is not None:
        pcov = cov[0,0]*chisquared
    else:
        pcov = numpy.nan
    
    if full_output:
        infodict['error'] = chisquared
        infodict['pcov'] = pcov
        return diameter, infodict
    else:
        return diameter
