import numpy
import scipy.signal
import scipy.ndimage

def patterson(I,M,image_threshold=0.,mask_smooth=1.,mask_threshold=1.,darkfield_x=None,darkfield_y=None,darkfield_sigma=None,normalize_median=False):
    I = numpy.clip(I, image_threshold, numpy.inf)
    K = kernel(M,smooth=mask_smooth,threshold=mask_threshold,x=darkfield_x,y=darkfield_y,sigma=darkfield_sigma)
    P = numpy.fft.fftshift(numpy.fft.fft2(numpy.fft.fftshift(K*I)))
    if normalize_median:
        P /= numpy.median(abs(P))
    return P

def kernel(mask,smooth=1.,threshold=1.,x=None,y=None,sigma=None):
    K = scipy.ndimage.gaussian_filter(numpy.array(mask,dtype="float"),smooth)
    t = K.max()*threshold
    K[K<t] = 0
    K = scipy.ndimage.gaussian_filter(numpy.array(K,dtype="float"),smooth)
    K /= K.max()
    if x!=None and y!=None and sigma!=None: 
        # multiply by gaussian darkfield kernel
        Ny, Nx = mask.shape
        cx = (Nx-1)/2
        cy = (Ny-1)/2
        NN = Nx*Ny
        X, Y = numpy.ogrid[0:Ny, 0:Nx]
        r = numpy.hypot(X - (cx + x), Y - (cx + y))
        G = numpy.exp(-r**2/(2*sigma**2))
        K = K*G
    return K
