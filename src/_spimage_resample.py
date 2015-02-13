import numpy

def downsample(array,factor,mode="pick",**kwargs):
    d = len(list(array.shape))
    if d == 2:
        return _downsample2d(array,factor,mode,**kwargs)
    elif d == 3:
        return _downsample2d(array,factor,mode,**kwargs)
    else:
        print "ERROR: Downsampling in %i dimensions is not supported."
        
def _downsample2d(array,factor,mode="pick",**kwargs):
    available_modes = ["pick","integrate"]#,"interpolate"]
    if not mode in available_modes:
        print "ERROR: %s is not a valid mode." % mode
        return
    mask         = kwargs.get("mask",None)
    bad_bits     = kwargs.get("bad_bits",None)
    min_N_pixels = kwargs.get("min_N_pixels",1)
    factor = int(round(factor0))
    if factor == 1:
        if mask == None:
            return array.copy()
        else:
            return [array.copy(),mask.copy()]
    if mask == None:
        mask = None
    else:
        mask = numpy.array(mask,dtype="int16")
    Nx = array.shape[1]
    Ny = array.shape[0]
    if mode == "pick": 
        Y,X = numpy.indices(array2.shape)
        pick = ((Y%factor == 0)*(X%factor == 0))
        print pick.shape
        Ny_new = pick[:,0].sum()
        Nx_new = pick[0,:].sum()
        pick = pick.flatten()
        A = array.flatten().copy()
        array_new = (A[pick]).reshape((Ny_new,Nx_new))
        if mask != None:
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
        Y /= factor
        X /= factor
        superp = Y*Nx+X
        superp_order = superp.argsort()
        A = A[superp_order]
        A = A.reshape((Nx_new*Ny_new,factor*factor))
        if mask == None:
            B = A.sum(1)
            return B.reshape((Ny_new,Nx_new))
        if mask != None:
            AM = numpy.zeros(shape=(Ny,Nx),dtype="int16")
            AM[:mask.shape[0],:mask.shape[1]] = mask[:,:]
            AM = AM.flat
            AM = AM[superp_order]
            AM = AM.reshape((Nx_new*Ny_new,factor*factor))
            if bad_bits == None:
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


def _downsample3d(array,factor,mask):
    available_modes = ["pick"]
    if not mode in available_modes:
        print "ERROR: %s is not a valid mode." % mode
        return
    factor = int(round(factor))
    if factor == 1:
        return array.copy()
    Nz = array.shape[0]
    Ny = array.shape[1]
    Nx = array.shape[2]
    Nz_new = int(numpy.ceil(1.0*Nz/factor))
    Ny_new = int(numpy.ceil(1.0*Ny/factor))
    Nx_new = int(numpy.ceil(1.0*Nx/factor))  
    array_new = numpy.zeros(shape=(Nz_new,Ny_new,Nx_new),dtype=array.dtype)
    for z_new in numpy.arange(0,Nz_new,1):
        for y_new in numpy.arange(0,Ny_new,1):
            for x_new in numpy.arange(0,Nx_new,1):
                z_min = z_new*factor
                z_max = min([(z_new+1)*factor,Nz])
                y_min = y_new*factor
                y_max = min([(y_new+1)*factor,Ny])
                x_min = x_new*factor
                x_max = min([(x_new+1)*factor,Nx])
                array_new[z_new,y_new,x_new] = array[z_min:z_max,y_min:y_max,x_min:x_max].sum()/(1.*mask[z_min:z_max,y_min:y_max,x_min:x_max].sum())
    return array_new
