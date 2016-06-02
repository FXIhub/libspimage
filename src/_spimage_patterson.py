import numpy,time
import _spimage_utils
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage.morphology import distance_transform_edt
# NOTE: The scipy FFT for real input is much faster than the implementation from numpy
from scipy.fftpack import fftn,fftshift
from scipy.signal import convolve2d

def patterson(image, mask, floor_cut=None, mask_smooth=1., darkfield_x=None, darkfield_y=None, darkfield_sigma=None, normalize_median=False, radial_boost=False, log_boost=False, gauss_damp=False, gauss_damp_sigma=None, gauss_damp_threshold=None, subtract_fourier_kernel=False, log_min=1., full_output=False):
    info = {}
    # Making sure that we leave memory of the input arrays untouched
    M = numpy.array(mask.copy(),dtype=numpy.float64)
    I = numpy.array(image.copy(),dtype=numpy.float64)
    # Clipping
    if floor_cut is not None:
        tmp = I <= floor_cut
        if tmp.sum() > 0:
            I[tmp] = 0.
    I0 = I.copy()
    if radial_boost or gauss_damp:
        Ny, Nx = I.shape
        cx = (Nx-1)/2
        cy = (Ny-1)/2
        X,Y = numpy.meshgrid(numpy.arange(Nx),numpy.arange(Ny))
        Rsq = (X - cx)**2 + (Y - cx)**2
    if radial_boost:
        # Radial boost: Multiply intensities by radius**4
        I = I * Rsq**2
    if gauss_damp:
        if gauss_damp_sigma is None:
            s = None
            if gauss_damp_threshold is None:
                s = (Nx + Ny)/4.
            else:
                Ir = I0.copy()
                Ir[M==0] = numpy.inf
                r,Ir = _spimage_utils.radial_mean(Ir,rout=True)
                inf = numpy.isfinite(Ir)==False
                if inf.sum() > 0:
                    tmp = Ir[inf==False][0]
                    for i in numpy.arange(len(r)):
                        if inf[i]:
                            Ir[i] = tmp
                        else:
                            tmp = Ir[i]
                #import matplotlib.pyplot as pypl
                #fig = pypl.figure()
                #ax = fig.add_subplot(111)
                #ax.plot(r,Ir)
                Ir = gaussian_filter(Ir,len(Ir)/100.)
                #ax.plot(r,Ir)
                tmp = Ir > gauss_damp_threshold
                if tmp.sum() != 0:
                    s = r[tmp][-1]
                #if s is not None:
                #    ax.axvline(s)
                #fig.savefig("./Ir_%i.png" % (numpy.random.randint(1000)))
                #pypl.close(fig)
        else:
            s = gauss_damp_sigma
        if s is not None:
            G = numpy.exp(-Rsq/(2*s**2))    
            I = I * G
        else:
            print "WARNING: Gaussian kernel could not be applied. Sigma undefined."
    if log_boost:
        tmp = I <= log_min
        if tmp.sum() > 0:
            I[tmp] = log_min
        I = numpy.log10(I)
    # Make kernel
    K = kernel(M,smooth=mask_smooth,x=darkfield_x,y=darkfield_y,sigma=darkfield_sigma)
    # Fourier transform
    P = fftshift(fftn(K*I))
    if subtract_fourier_kernel:
        fK = fftshift(fftn(K))
        if full_output:
            info["patterson0"] = P.copy()
            if normalize_median:
                info["patterson0"] /= numpy.median(abs(P))
        from scipy.optimize import leastsq
        err = lambda c: abs(P - c*fK).flatten()
        c = leastsq(err,1.)[0]
        P = P - c * fK
    if normalize_median:
        # Normalization
        P /= numpy.median(abs(P))
    if full_output:
        info["kernel"] = K
        info["intensities"] = I0
        info["intensities_times_kernel"] = K*I
        return P, info
    else:
        return P

def kernel(mask,smooth=5.,x=None,y=None,sigma=None):
    K = distance_transform_edt(mask, sampling=None, return_distances=True, return_indices=False, distances=None, indices=None)
    K = K > (smooth/5.)
    K = gaussian_filter(numpy.array(K,dtype=numpy.float64),smooth)
    if x is not None and y is not None and sigma is not None: 
        # multiply by gaussian darkfield kernel
        Ny, Nx = M.shape
        cx = (Nx-1)/2
        cy = (Ny-1)/2
        X,Y = numpy.ogrid[0:Ny, 0:Nx]
        Rsq = (X - (cx + x))**2 + (Y - (cx + y))**2
        G = numpy.exp(-Rsq/(2*sigma**2))
        K = K*G
    return K

