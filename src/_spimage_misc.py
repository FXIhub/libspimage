import numpy as np
import scipy as sp
from scipy import constants
__all__ = ['grid', 'rgrid', 'rmask', 'DICT_physical_constants']

def grid(shape, center=None):
    if center is None: center = (0,0)
    yy, xx = np.indices(shape)
    xx = xx.astype(np.float) - ( center[0] + shape[1]/2. - 0.5 )
    yy = yy.astype(np.float) - ( center[1] + shape[0]/2. - 0.5 )
    return yy, xx

def rgrid(shape, center=None):
    if center is None: center = (0,0)
    yy, xx = grid(shape, center)
    rr = np.sqrt( (xx**2) + (yy**2) )
    rr[np.where(rr == 0)] = 1e-20
    return rr

def rmask(shape, r, center=None):
    if center is None: center = (0,0)
    yy, xx = grid(shape, center)
    rmask = r**2 >= (yy**2 + xx**2)
    return rmask

DICT_physical_constants = {'e':constants.e,
                           'c':constants.c,
                           'h':constants.h,
                           're':constants.value("classical electron radius"),
                           'barn':1E-28,
                           'u':constants.value("atomic mass constant")}


