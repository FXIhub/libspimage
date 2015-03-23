import numpy
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage.morphology import distance_transform_edt

def patterson(image, mask, floor_cut=0., mask_smooth=1., darkfield_x=None, darkfield_y=None, darkfield_sigma=None, normalize_median=False, radial_boost=False, subtract_fourier_kernel=False, full_output=False):
    info = {}
    # Making sure that we leave memory of the input arrays untouched
    M = numpy.array(mask.copy(),dtype=numpy.float64)
    I = numpy.array(image.copy(),dtype=numpy.float64)
    # Clipping
    #I = numpy.clip(I*M, floor_cut, numpy.inf)
    tmp = I <= floor_cut
    if tmp.sum() > 0:
        I[tmp] = 0.
    if radial_boost:
        # Radial boost: Multiply intensities by radius**4
        Ny, Nx = I.shape
        cx = (Nx-1)/2
        cy = (Ny-1)/2
        X,Y = numpy.meshgrid(numpy.arange(Nx),numpy.arange(Ny))
        R = numpy.hypot(X - cx, Y - cx)
        s = 200.
        G = numpy.exp(-R**2/(2*s**2))
        I = I * R**4 * G
    # Make kernel
    K = kernel(M,smooth=mask_smooth,x=darkfield_x,y=darkfield_y,sigma=darkfield_sigma)
    # Fourier transform
    P = numpy.fft.fftshift(numpy.fft.fftn(K*I))
    if subtract_fourier_kernel:
        fK = numpy.fft.fftshift(numpy.fft.fftn(K))
        if full_output:
            info["patterson0"] = P.copy()
            if normalize_median:
                info["patterson0"] /= numpy.median(abs(P))
        s_fK = abs(fK[fK.shape[0]/2,:]).sum() + abs(fK[:,fK.shape[1]/2]).sum()
        s_P = abs(P[P.shape[0]/2,:]).sum() + abs(P[:,P.shape[1]/2]).sum()
        P = abs(P) - abs(fK / s_fK * s_P)
    if normalize_median:
        # Normalization
        P /= numpy.median(abs(P))
    if full_output:
        info["kernel"] = K
        info["intensities"] = K*I
        P_nk = numpy.fft.fftshift(numpy.fft.fft2(numpy.fft.fftshift(numpy.array(image,dtype=numpy.float64))))
        if normalize_median:
            P_nk /= numpy.median(abs(P_nk))
        info["patterson_no_kernel"] = P_nk
        return P, info
    else:
        return P

def kernel(mask,smooth=1.,x=None,y=None,sigma=None):
    K = distance_transform_edt(mask, sampling=None, return_distances=True, return_indices=False, distances=None, indices=None)
    K = K > (smooth * 5)
    K = gaussian_filter(numpy.array(K,dtype=numpy.float64),smooth)
    if x is not None and y is not None and sigma is not None: 
        # multiply by gaussian darkfield kernel
        Ny, Nx = M.shape
        cx = (Nx-1)/2
        cy = (Ny-1)/2
        X,Y = numpy.ogrid[0:Ny, 0:Nx]
        R = numpy.hypot(X - (cx + x), Y - (cx + y))
        G = numpy.exp(-R**2/(2*sigma**2))
        K = K*G
    return K
