import numpy
import logging

def downsample(array,factor,mode="pick",**kwargs):
    d = len(list(array.shape))
    if d == 2:
        return _downsample2d(array,factor,mode,**kwargs)
    elif d == 3:
        return _downsample3d(array,factor,mode,**kwargs)
    else:
        print("ERROR: Downsampling in %i dimensions is not supported.")
        
def _downsample2d(array,factor,mode="pick",**kwargs):
    available_modes = ["pick","integrate"]#,"interpolate"]
    if not mode in available_modes:
        print("ERROR: %s is not a valid mode." % mode)
        return
    mask         = kwargs.get("mask",None)
    bad_bits     = kwargs.get("bad_bits",None)
    min_N_pixels = kwargs.get("min_N_pixels",1)
    factor = int(round(factor))
    if factor == 1:
        if mask is None:
            return array.copy()
        else:
            return [array.copy(),mask.copy()]
    if mask is None:
        mask = None
    else:
        mask = numpy.array(mask,dtype="int16")
    Nx = array.shape[1]
    Ny = array.shape[0]
    if mode == "pick": 
        Y,X = numpy.indices(array.shape)
        pick = ((Y%factor == 0)*(X%factor == 0))
        Ny_new = pick[:,0].sum()
        Nx_new = pick[0,:].sum()
        pick = pick.flatten()
        A = array.flatten().copy()
        array_new = (A[pick]).reshape((Ny_new,Nx_new))
        if mask is not None:
            M = mask.flatten().copy()
            mask_new = (M[pick]).reshape((Ny_new,Nx_new))
            return [array_new,mask_new]
        else:
            return array_new
    elif mode == "integrate": # non-conservative if mask is given
        Nx_new = int(numpy.ceil(Nx/float(factor)))
        Ny_new = int(numpy.ceil(Ny/float(factor)))
        Nx = Nx_new * factor
        Ny = Ny_new * factor
        A = numpy.zeros(shape=(Ny,Nx),dtype=array.dtype)
        A[:array.shape[0],:array.shape[1]] = array[:,:]
        A = A.flat
        Y,X = numpy.indices((Ny,Nx))
        Y = Y.flatten()
        X = X.flatten()
        Y //= factor
        X //= factor
        superp = Y*Nx+X
        superp_order = superp.argsort()
        A = A[superp_order]
        A = A.reshape((Nx_new*Ny_new,factor*factor))
        if mask is None:
            B = A.sum(1)
            return B.reshape((Ny_new,Nx_new))
        if mask is not None:
            AM = numpy.zeros(shape=(Ny,Nx),dtype="int16")
            AM[:mask.shape[0],:mask.shape[1]] = mask[:,:]
            AM = AM.flat
            AM = AM[superp_order]
            AM = AM.reshape((Nx_new*Ny_new,factor*factor))
            if bad_bits is None:
                B = (A*AM).sum(1)
                BN = AM.sum(1)
                BM = BN != 0
            else:
                B = (A*((AM & bad_bits) == 0)).sum(1)
                BN = ((AM & bad_bits) == 0).sum(1)
                BM = AM[:,0]
                for i in range(1,factor*factor):
                    BM |= AM[:,i]
                BM[BN >= min_N_pixels] = BM[BN >= min_N_pixels] & ~bad_bits
            B[BN >= min_N_pixels] = B[BN >= min_N_pixels] * factor/numpy.float64(BN[BN >= min_N_pixels])
            return [B.reshape((Ny_new,Nx_new)),BM.reshape((Ny_new,Nx_new))]


def _downsample3d(array,factor,mode="pick",**kwargs):
    available_modes = ["pick","integrate"]#,"interpolate"]                                                                                                                      
    if not mode in available_modes:
        print("ERROR: %s is not a valid mode." % mode)
        return
    mask         = kwargs.get("mask",None)
    bad_bits     = kwargs.get("bad_bits",None)
    min_N_pixels = kwargs.get("min_N_pixels",1)
    factor = int(round(factor))
    if factor == 1:
        if mask is None:
            return array.copy()
        else:
            return [array.copy(),mask.copy()]
    if mask is None:
        mask = None
    else:
        mask = numpy.array(mask,dtype="int16")
    Nx = array.shape[2]
    Ny = array.shape[1]
    Nz = array.shape[0]
    if mode == "pick":
        Z,Y,X = numpy.indices(array.shape)
        pick = ((Z%factor == 0)*(Y%factor == 0)*(X%factor == 0))
        Nz_new = pick[:,:,0].sum()
        Ny_new = pick[:,0,:].sum()
        Nx_new = pick[0,:,:].sum()
        pick = pick.flatten()
        A = array.flatten().copy()
        array_new = (A[pick]).reshape((Nz_new,Ny_new,Nx_new))
        if mask is not None:
            M = mask.flatten().copy()
            mask_new = (M[pick]).reshape((Nz_new,Ny_new,Nx_new))
            return [array_new,mask_new]
        else:
            return array_new
    elif mode == "integrate": # non-conservative if mask is given                                                                                                               
        Nx_new = int(numpy.ceil(Nx/float(factor)))
        Ny_new = int(numpy.ceil(Ny/float(factor)))
        Nz_new = int(numpy.ceil(Nz/float(factor)))
        Nx = Nx_new * factor
        Ny = Ny_new * factor
        Nz = Nz_new * factor
        A = numpy.zeros(shape=(Nz,Ny,Nx),dtype=array.dtype)
        A[:array.shape[0],:array.shape[1],:array.shape[2]] = array[:,:,:]
        A = A.flat
        Z,Y,X = numpy.indices((Nz,Ny,Nx))
        Z = Z.flatten()
        Y = Y.flatten()
        X = X.flatten()
        Z /= factor
        Y /= factor
        X /= factor
        superp = Z*Ny*Nx+Y*Nx+X
        superp_order = superp.argsort()
        A = A[superp_order]
        A = A.reshape((Nx_new*Ny_new*Nz_new,factor*factor*factor))
        if mask is None:
            B = A.sum(1)
            return B.reshape((Nz_new,Ny_new,Nx_new))
        if mask is not None:
            AM = numpy.zeros(shape=(Nz,Ny,Nx),dtype="int16")
            AM[:mask.shape[0],:mask.shape[1],:mask.shape[2]] = mask[:,:,:]
            AM = AM.flat
            AM = AM[superp_order]
            AM = AM.reshape((Nx_new*Ny_new*Nz_new,factor*factor*factor))
            if bad_bits is None:
                B = (A*AM).sum(1)
                BN = AM.sum(1)
                BM = BN != 0
            else:
                B = (A*((AM & bad_bits) == 0)).sum(1)
                BN = ((AM & bad_bits) == 0).sum(1)
                BM = AM[:,:,0]
                for i in range(1,factor*factor*factor):
                    BM |= AM[:,:,i]
                BM[BN >= min_N_pixels] = BM[BN >= min_N_pixels] & ~bad_bits
            B[BN >= min_N_pixels] = B[BN >= min_N_pixels] * factor/numpy.float64(BN[BN >= min_N_pixels])
            return [B.reshape((Nz_new,Ny_new,Nx_new)),BM.reshape((Nz_new,Ny_new,Nx_new))]

def crop(pattern,cropLength,center='middle',bg=0):
    if cropLength > min(list(pattern.shape)):
        print("WARNING: Could not crop image to larger size.")
        return pattern.copy()
    if (numpy.asarray(list(pattern.shape)) == cropLength).all():
        return numpy.copy(pattern)
    
    temp = pattern.copy()
    if center == 'middle':
        x_center = (pattern.shape[1] - 1)/2.
        y_center = (pattern.shape[0] - 1)/2.
    else:
        x_center = center[1]
        y_center = center[0]

    #x_start = (pattern.shape[1]-cropLength)/2
    #y_start = (pattern.shape[0]-cropLength)/2
    x_start = max([int(round(x_center - cropLength/2)), 0])
    y_start = max([int(round(y_center - cropLength/2)), 0])
    x_stop =  x_start+cropLength
    y_stop =  y_start+cropLength

    patternCropped = numpy.ones(shape=(cropLength,cropLength),dtype=pattern.dtype)*bg
    patternCropped = temp[y_start:y_stop,x_start:x_stop]
    return patternCropped

def binImage(img, binning, msk=None, output_binned_mask=False):
    """ This function bins a 2D image. The image will be cropped before binning if the dimensions are not a multiple of the binning factor. 
    If a mask is provided the mask is applied to the image before binning. 
    Binned pixels from partly masked out regions are scaled up such that the mean value within the bin matches the binned value divided by the binning factor squared.

    Args:
        :img:       Array of the native image
        :binning(ing): Binning factor
    Kwargs:
        :msk(int, bool):  Mask array, True (or 1) means it is a valid pixel
        :output_binned_mask: (bool): Toggle if you want that the binned mask is return besides the binned image
    
    :Authors:
        Max F. Hantke (hantke@xray.bmc.uu.se)
    """
    nx = img.shape[1]
    ny = img.shape[0]
    valid_input = True
    img_new = img
    msk_new = msk
    if binning > nx or binning > ny:
        valid_input = False
        logging.warning("Image with dimensions %i x %i too small to be binned by %i x %i.", ny, nx, binning, binning)
    if msk is not None:
        # Check for matching dimensions
        if msk.shape[0] != img.shape[0] or msk.shape[1] != img.shape[1]:
            logging.error("Dimensions of image (%i x %i) and mask (%i x %i) do not match.", img.shape[0], img.shape[1], msk.shape[0], msk.shape[1])
            valid_input = False
    if valid_input:
        # Crop image such that dimensions are multiples of binning
        nx_new = nx - nx % binning
        ny_new = ny - ny % binning
        img_new = img[:ny_new,:nx_new]
        if msk is not None:
            # Crop mask
            msk_new = msk[:ny_new,:nx_new]
            # Set masked out values in image to zero
            img_new = img_new * msk_new
            # Bin mask
            msk_new = numpy.array(msk_new, dtype="int")
            msk_new = msk_new.reshape(ny_new // binning, binning, nx_new // binning, binning)
            msk_new = msk_new.sum(axis=3).sum(axis=1)
        # New dimensions for binned pixels
        img_new = img_new.reshape(ny_new // binning, binning, nx_new // binning, binning)
        img_new = img_new.sum(axis=3).sum(axis=1)
        if msk is not None:
            img_new *= float(binning**2) / (msk_new + numpy.finfo("float").eps)
    if output_binned_mask:
        return img_new, msk_new
    else:
        return img_new

def _testBinImage(binning,masking=True):
    from scipy import misc
    l1 = misc.lena()
    l1 = l1[:511,:480]
    l1 = l1[:l1.shape[0] - l1.shape[0] % binning,:l1.shape[1] - l1.shape[1] % binning]
    S1 = l1.sum()
    m1 = m2 = None
    if masking:
        m1 = numpy.random.random(l1.shape)
        m1 = numpy.array(numpy.round(m1),dtype="bool")
        l1 *= m1
    l2,m2 = binImage(l1,binning,m1,output_binned_mask=True)
    S2 = l2.sum()
    print("Sum original (cropped) image: %f" % S1)
    print("Sum binned image: %f" % S2)
    print("Sum difference: %.3f %%" % (100.*abs(S1-S2).sum()/2./(S1+S2).sum()))
    return l1,l2,m1,m2
    

def _radialImage(img,mode="mean",cx=None,cy=None,msk=None,output_r=False):
    """ This function calculates a radial representation of a 2D image.
    If a mask is provided the mask is applied to the image before projection. The output of radii that include only masked out pixels is set to nan.
    The center is put into the middle of the respective axis if a center coordinate is not specified (set to None).

    Args:
        :img:       Array of the native image
        :mode(str): Projection mode can be mean, sum, std or median
    Kwargs:
        :cx(float): Center x coordinate
        :cy(float): Center y coordinate
        :msk(int, bool):  Mask array, True (or 1) means it is a valid pixel
        :output_r(book): Set to true if also the array of the radii shall be in the output

    :Authors:
        Max F. Hantke (hantke@xray.bmc.uu.se)
    """
    if mode == "mean": f = numpy.mean
    elif mode == "sum": f = numpy.sum
    elif mode == "std": f = numpy.std
    elif mode == "median": f = numpy.median
    else:
        logging.error("ERROR: No valid mode given for radial projection.")
        return None
    if cx is None: cx = (img.shape[1]-1)/2.0
    if cy is None: cy = (img.shape[0]-1)/2.0
    X,Y = numpy.meshgrid(numpy.arange(img.shape[1]),numpy.arange(img.shape[0]))
    R = numpy.sqrt((X - float(cx))**2 + (Y - float(cy))**2)
    R = R.round()
    if msk is not None:
        if (msk == 0).sum() > 0:
            R[msk == 0] = -1
    radii = numpy.arange(R.min(),R.max()+1,1)
    if radii[0] == -1:
        radii = radii[1:]
    values = numpy.zeros_like(radii)
    for i in range(0,len(radii)):
        tmp = R==radii[i]
        if tmp.sum() > 0:
            values[i] = f(img[tmp])
        else:
            values[i] = numpy.nan
    if (numpy.isfinite(values) == False).sum() > 0:
        tmp = numpy.isfinite(values)
        values = values[tmp]
        radii  = radii[tmp]
    if output_r:
        return radii,values
    else:
        return values

def _radialImage_3D(img,mode="mean",cx=None,cy=None,cz=None,msk=None,output_r=False):
    if mode == "mean": f = numpy.mean
    elif mode == "sum": f = numpy.sum
    elif mode == "std": f = numpy.std
    elif mode == "median": f = numpy.median
    else:
        logging.error("ERROR: No valid mode given for radial projection.")
        return None
    if cx is None: cx = (img.shape[2]-1)/2.0
    if cy is None: cy = (img.shape[1]-1)/2.0
    if cz is None: cz = (img.shape[0]-1)/2.0
    X,Y,Z = numpy.meshgrid(numpy.arange(img.shape[2]),numpy.arange(img.shape[1]),numpy.arange(img.shape[0]))
    R = numpy.sqrt((X - float(cx))**2 + (Y - float(cy))**2 + (Z - float(cz))**2)
    R = R.round()
    if msk is not None:
        if (msk == 0).sum() > 0:
            R[msk == 0] = -1
    radii = numpy.arange(R.min(),R.max()+1,1)
    if radii[0] == -1:
        radii = radii[1:]
    values = numpy.zeros_like(radii)
    for i in range(0,len(radii)):
        tmp = R==radii[i]
        if tmp.sum() > 0:
            values[i] = f(img[tmp])
        else:
            values[i] = numpy.nan
    if (numpy.isfinite(values) == False).sum() > 0:
        tmp = numpy.isfinite(values)
        values = values[tmp]
        radii  = radii[tmp]
    if output_r:
        return radii,values
    else:
        return values


def radialSumImage(img,**kwargs):
    if len(img.shape) == 2:
        return _radialImage(img,mode="sum",**kwargs)

    if len(img.shape) == 3:
        return _radialImage_3D(img,mode="sum",**kwargs)

def radialStd(img,**kwargs):
    if len(img.shape) == 2:
        return _radialImage(img,mode="std",**kwargs)

    if len(img.shape) == 3:
        return _radialImage_3D(img,mode="std",**kwargs)

def radialMeanImage(img,**kwargs):
    if len(img.shape) == 3:
        return _radialImage_3D(img,mode="mean",**kwargs)

    if len(img.shape) == 2:
        return _radialImage(img,mode="mean",**kwargs)

def radialMedianImage(img,**kwargs):
    if len(img.shape) == 2:
        return _radialImage(img,mode="median",**kwargs)

    if len(img.shape) == 3:
        return _radialImage_3D(img,mode="median",**kwargs)

def _testRadialImage(Nx=100,Ny=103,cx=45.3,cy=43.2):
    X,Y = numpy.meshgrid(numpy.arange(Nx),numpy.arange(Ny))
    R = numpy.sqrt((X-cx)**2+(Y-cy)**2)
    img = R.round()
    print("Std-zero-test: Succeeded = %s" % (radialStdImage(img,cx=cx,cy=cy).sum()==0.))
    return img,radialMeanImage(img,cx=cx,cy=cy,output_r=True)
