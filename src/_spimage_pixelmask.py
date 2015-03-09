import numpy
from _spimage_utils import _turn180

class PixelMask:
    # CXI pixelmask bits
    PIXEL_IS_PERFECT = 0
    PIXEL_IS_INVALID = 1
    PIXEL_IS_SATURATED = 2
    PIXEL_IS_HOT = 4
    PIXEL_IS_DEAD = 8
    PIXEL_IS_SHADOWED = 16
    PIXEL_IS_IN_PEAKMASK = 32
    PIXEL_IS_TO_BE_IGNORED = 64
    PIXEL_IS_BAD = 128
    PIXEL_IS_OUT_OF_RESOLUTION_LIMITS = 256
    PIXEL_IS_MISSING = 512
    PIXEL_IS_NOISY = 1024
    PIXEL_IS_ARTIFACT_CORRECTED = 2048
    PIXEL_FAILED_ARTIFACT_CORRECTION = 4096
    PIXEL_IS_PEAK_FOR_HITFINDER = 8192
    PIXEL_IS_PHOTON_BACKGROUND_CORRECTED = 16384
    # Cummulative bit options
    PIXEL_IS_IN_MASK = PIXEL_IS_INVALID |  PIXEL_IS_SATURATED | PIXEL_IS_HOT | PIXEL_IS_DEAD | PIXEL_IS_SHADOWED | PIXEL_IS_IN_PEAKMASK | PIXEL_IS_TO_BE_IGNORED | PIXEL_IS_BAD | PIXEL_IS_MISSING

def select_central_speckle(image, mask, threshold=0.1, radius_extension_factor=1.4, symmetrize_mask=True, cx=None, cy=None, **kwargs): 
    _mask = mask.copy()
    _mask = numpy.array(_mask,dtype="int")
    if cx is None:
        _cx = (mask.shape[1]-1)/2.
    else:
        _cx = cx
    if cy is None:
        _cy = (mask.shape[1]-1)/2.
    else:
        _cy = cy
    if symmetrize_mask:
        _mask *= symmetrize(_mask,_cx,_cy)
    image_threshold = threshold*image[_mask==1].max()
    temp = (image*_mask)>image_threshold
    if temp.sum() > 0:
        _mask *= disk_mask(temp,_cx,_cy,radius_extension_factor)
    return _mask

def symmetrize(mask,cx,cy):
    _mask = mask.copy()
    _mask *= _turn180(_mask,cx,cy)
    return _mask

def disk_mask(mask,cx,cy,radius_extension_factor):
    X,Y = numpy.meshgrid(numpy.arange(mask.shape[1]),numpy.arange(mask.shape[0]))
    X -= cx
    Y -= cy
    R = numpy.sqrt(X**2+Y**2)
    r_max = R[mask==1].max()
    _mask = R <= r_max*radius_extension_factor
    return _mask

    
