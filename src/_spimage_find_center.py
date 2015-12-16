import numpy as np
import scipy as sp
import spimage
from _spimage_conventions import pos_to_center,center_to_pos
from _spimage_utils import _symmetrize,_turn180,_gaussian_smooth_2d1d

def find_center(img, msk, method=None, errout=False, **kwargs):
    """
    Find the center of a diffraction pattern.

    usage:
    ======
    x,y = find_center(img, msk)
    x,y = find_center(img, msk, method='quadrant',  x0=0, y0=0, dmax=None, threshold=None, solver='L-BFGS-B')

    x,y = find_center(img, msk, method='blurred', x0=0, y0=0, threshold=None, blur_radius=4, dmax=5)
    x,y = find_center(img, msk, method='pixelwise_fast', x0=0, y0=0, dmax=5, rmax=None)
    x,y = find_center(img, msk, method='pixelwise_slow', x0=0, y0=0, dmax=5, rmax=None)
    """

    # Default method for center finding
    if method is None: method = 'quadrant'

    # Find center using "quadrant" method
    if method == 'quadrant':
        x,y,e = find_center_quadrant(img, msk, **kwargs)
    # Find center using fast C implementation of "pixelwise" method
    elif method == 'pixelwise_fast':
        x,y,e = find_center_pixelwise_fast(img, msk, **kwargs)
    # Find center using slow implementation of "pixelwise" method
    elif method == 'pixelwise_slow':
        x,y,e = find_center_pixelwise_slow(img, msk, **kwargs)
    # Find center using blurred version of 'pixelwise' method
    elif method == 'blurred':
        x,y,e = find_center_blurred(img, msk, **kwargs)
    # Return 0,0 if method is not defined
    else:
        x,y,e = (0,0,0)
        print "There is no center finding method %s" %method

    # Check for reasonable numbers
    if abs(x) > img.shape[1]/2: x = 0
    if abs(y) > img.shape[0]/2: y = 0

    if errout:
        return (x,y,e)
    else:
        return (x,y)

def find_center_quadrant(img, msk, dmax=5, x0=0, y0=0, threshold=None, solver='L-BFGS-B'):
    """
    Find the center of a diffraction pattern using the quadrant method.
    For every possible center shift (within - dmax ... + dmax) a centrosymmetric mask is calculated.
    The image is divided into four quadrants A,B,C,D around any possible center position:
    
    +----+----+
    | A  | B  |
    +----+----|
    | C  | D  |
    +----+----+

    Depending of different center positions, the cost function

    E = |sum(w_A*A)-sum(w_D*D)|^2 + |sum(w_B*B) - sum(w_C*C)|^2 ) / sum(wA + wB + wC + wD)

    with w_i being the respective centrosymmetric masks for i=A,B,C,D
    is minimized using a given solver (default = 'L-BFGS-B').

    usage: 
    ======
    x,y = find_center_quadrant(img, msk,  x0=0, y0=0, dmax=None, threshold=None, solver='L-BFGS-B')
    """    
    class CentroSymmetricMask:
        def __init__(self, mask, dx, dy):
            mask = mask.astype(np.bool)
            self.omask = mask
            self.ny, self.nx = mask.shape
            self.dx = dx
            self.dy = dy
            self.define_mask()

        def define_mask(self):
            csmask0 = sp.flipud(sp.fliplr(self.omask))
            csmask = np.zeros((self.ny+4*self.dy, self.nx+4*self.dx)).astype(np.bool)
            csmask[self.dy:self.dy+self.ny,self.dx:self.dx+self.nx] = csmask0
            self.csmask = csmask
                                                                                                                                                                                                                                                                
        def get(self, x0=0, y0=0):
            return self.omask & self.csmask[self.dy-2*y0:self.dy+self.ny-2*y0,self.dx-2*x0:self.dx+self.nx-2*x0]

    class Minimizer:
        def __init__(self,img, msk, x0, y0, maxshift, solver):
            self.csmask = CentroSymmetricMask(msk, 2*maxshift, 2*maxshift)
            self.image = img
            self.Ny, self.Nx = img.shape
            self.x0_initial = x0
            self.y0_initial = y0
            self.maxshift = maxshift
            self.solver = solver

        def manipulate(self, threshold=None):
            if threshold is not None:
                self.image[self.image < threshold] = 0.
            #self.image = np.log(self.image)
            
        def update_center(self, x0,y0):
            self.yt = max(2*y0, 0)
            self.yb = min(self.Ny + 2*y0, self.Ny)
            self.ym = (self.yb + self.yt)/2
            self.xl = max(2*x0, 0)
            self.xr = min(self.Nx + 2*x0, self.Nx)
            self.xm = (self.xl + self.xr)/2

        def update_mask(self, x0,y0):
            self.mask = self.csmask.get(x0,y0)

        def error(self, p):
            [x0, y0] = p
            self.update_center(x0,y0)
            self.update_mask(x0,y0)

            wA = self.mask[self.yt:self.ym, self.xl:self.xm]
            wB = self.mask[self.yt:self.ym, self.xm:self.xr]
            wC = self.mask[self.ym:self.yb, self.xl:self.xm]
            wD = self.mask[self.ym:self.yb, self.xm:self.xr]

            A = self.image[self.yt:self.ym, self.xl:self.xm][wA]
            B = self.image[self.yt:self.ym, self.xm:self.xr][wB]
            C = self.image[self.ym:self.yb, self.xl:self.xm][wC]
            D = self.image[self.ym:self.yb, self.xm:self.xr][wD]
            
            norm = 2*(wA.sum() + wB.sum())
            error = np.sqrt( abs(A.sum() - D.sum())**2 + abs(B.sum() - C.sum())**2 ) / norm
            return error
            
        def error_smooth(self, p):
            [x0, y0] = p
            x0f = np.floor(x0)
            x0c = x0f + 1
            y0f = np.floor(y0)
            y0c = y0f + 1
            err_ff = self.error([x0f, y0f])
            err_cc = self.error([x0c, y0c])
            err_cf = self.error([x0c, y0f])
            err_fc = self.error([x0f, y0c])
            wff = (x0c - x0)  * (y0c - y0)
            wcc = (x0  - x0f) * (y0  - y0f)
            wcf = (x0  - x0f) * (y0c - y0)
            wfc = (x0c - x0)  * (y0  - y0f)
            error = wff*err_ff + wcc*err_cc + wcf*err_cf + wfc*err_fc
            return error

        def error_and_gradient(self,p):
            err = self.error(p)
            dx = self.egrgror([p[0]+1, p[1]]) - err
            dy = self.error([p[0]  , p[1]+1]) - err
            return err, np.array([dx,dy])

        def start(self):
            self.res = sp.optimize.minimize(self.error_smooth, np.array([self.x0_initial,self.y0_initial]), method=self.solver, jac=False, options={'disp':False}, bounds=[(-self.maxshift+1, self.maxshift-1), (-self.maxshift+1, self.maxshift-1)])

        def error_landscape(self):
            XY = np.mgrid[-self.maxshift:self.maxshift:(2*self.maxshift+1)*1j]
            E = np.zeros((2*self.maxshift+1,2*self.maxshift+1))
            for j in range(E.shape[0]):
                for i in range(E.shape[1]):
                    E[j,i] = self.error([XY[i],XY[j]])
            return E

    m = Minimizer(img, msk, x0, y0, dmax, solver)
    m.manipulate(threshold)
    m.start()
    x = m.res["x"][0]
    y = m.res["x"][1]
    e = m.error_smooth(m.res["x"])
    return (x,y,e)

def find_center_pixelwise_fast(img, msk, x0=0, y0=0, dmax=5, rmax=None):
    """
    Find center of diffraction pattern using a pixelwise comparison of centry-symmetric pixels.

    This is a faster C implementation.

    usage:
    ======
    x,y = find_center_pixelwise_fast(img, msk, x0, y0, dmax=5, rmax=None)
    """
    if rmax is not None: msk &= (spimage.rgrid(msk.shape, (x0,y0)) < rmax)
    I = spimage.sp_image_alloc(int(np.ceil(img.shape[1])), int(np.ceil(img.shape[0])), 1)
    I.image[:] = img
    I.mask[:]  = msk
    I.detector.image_center[:] = np.array([center_to_pos(x0,img.shape[1]), center_to_pos(y0,img.shape[0]), 0 ])
    success = spimage.sp_find_center_refine(I, dmax, 0, None)
    x = pos_to_center(I.detector.image_center[0],img.shape[1])
    y = pos_to_center(I.detector.image_center[1],img.shape[0])
    spimage.sp_image_free(I)
    return (x,y,success)

def find_center_pixelwise_slow(img, msk, x0=0, y0=0, dmax=5, rmax=None):
    """
    Find center of diffraction pattern using pixelwise comparison of centro-symmetric pixels.
    
    This is the original python implementation taken from owl/fit.py
    
    usage:
    ======
    x,y = find_center_pixelwise_slow(img, msk, x0, y0, dmax=5, rmax=None)
    """
    s = img.shape
    if rmax is None: rmax = np.sqrt(2)*max(s)
    cx_g = (s[1]-1)/2.+x0
    cy_g = (s[0]-1)/2.+y0
    cx_g = np.round(cx_g*2)/2.
    cy_g = np.round(cy_g*2)/2.
    ddc = 0.5
    N_sam1= int(np.round(2*dmax/ddc))+1
    cx_sam1 = np.linspace(cx_g-dmax,cx_g+dmax,N_sam1)
    cy_sam1 = np.linspace(cy_g-dmax,cy_g+dmax,N_sam1)
    N_sam2= int(np.round(4*dmax/ddc))+1
    cx_sam2 = np.linspace(cx_g-dmax*2,cx_g+dmax*2,N_sam2)
    cy_sam2 = np.linspace(cy_g-dmax*2,cy_g+dmax*2,N_sam2)

    # extend mask so that every pixel has a partner at every possible center
    msk_ext = msk.copy()
    for cy in cy_sam2:
        for cx in cx_sam2:
            msk_ext *= _symmetrize(msk,cx,cy)
    Nme = msk_ext.sum()
    errs = np.zeros(shape=(N_sam1,N_sam1))
    r_max_sq = rmax**2

    X,Y = np.meshgrid(np.arange(img.shape[1]),np.arange(img.shape[0]))    
    for cx,icx in zip(cx_sam1,np.arange(N_sam1)):
        for cy,icy in zip(cy_sam1,np.arange(N_sam1)):
            
            # SLOW CODE
            #for x,y,v1 in zip((Xme-cx),(Yme-cy),imgme):
            #    M = (Xm-cx==-x)*(Ym-cy==-y)
            #    if M.sum() == 1:
            #        v2 = imgm[M==True]
            #        errs[icy,icx] += abs(v1-v2)
            #    else:
            #        print x,y,M.sum()
            
            # FAST CODE (does the same)
            r_sq = ((X-cx)**2+(Y-cy)**2)
            rmsk = r_sq < r_max_sq
            img_turned = _turn180(img,cx,cy)
            diff = abs((img-img_turned)*msk_ext*rmsk)
            errs[icy,icx] = diff.sum()
            #print cx,cy,errs[icy,icx]
  
    errs_sm = _gaussian_smooth_2d1d(errs,dmax)
    i_min = errs.flatten().argmin()
    cxi_min = i_min % N_sam1
    cyi_min = i_min/N_sam1
    #if full_output:
    #    D = {}
    #    D["msk_ext"] = msk_ext
    #    D["img_turned_msk_ext"] = img_turned*msk_ext
    #    D["img_msk_ext"] = img*msk_ext
    #    D["errmap"] = errs/errs.max()
    #    D["errmap_sm"] = gaussian_smooth_2d1d(errs,3.)
    #    D["errmap_sm"] /= D["errmap_sm"].max()
    #    D["errmap_min"] = errs.min()
    #    return cx_sam1[cxi_min],cy_sam1[cyi_min],D
    #else:
    cx_r = cx_sam1[cxi_min]
    cy_r = cy_sam1[cyi_min]
    x = cx_r-(s[1]-1)/2.
    y = cy_r-(s[0]-1)/2.
    return (x,y, errs.flatten()[i_min])

def find_center_blurred(img, msk, x0=0, y0=0, threshold=None, blur_radius=4., dmax=5):
    """
    Find the center using blurred version of 'pixelwise' method.

    usage: 
    ======
    x,y = find_center_blurred(img, msk, x0=0, y0=0, threshold=None, blur_radius=4, dmax=5)
    """
    if threshold is None: threshold = img.min()
    I = spimage.sp_image_alloc(img.shape[1],img.shape[0],1)
    I.image[:] = img * (img >= threshold)
    I.mask[:] = msk[:]
    I.detector.image_center[:] = np.array([center_to_pos(x0,img.shape[1]),
                                           center_to_pos(y0,img.shape[0]), 0 ])
    kernel = spimage.sp_gaussian_kernel(float(blur_radius),int(blur_radius*8+1),int(blur_radius*8+1),1)
    c = spimage.sp_image_convolute_with_mask(I,kernel,np.array([1,1,1]).astype(np.int32))
    ds = spimage.sp_image_alloc(img.shape[1]/4,img.shape[0]/4,1)
    ds.image[:] = c.image[:-3:4,:-3:4]
    ds.mask[:] = c.mask[:-3:4,:-3:4]
    ds.detector.image_center[:] = c.detector.image_center[:] / 4.0
    spimage.sp_find_center_refine_minimal_mask(ds, int(1+dmax/4), 0)
    c.detector.image_center[:] = ds.detector.image_center[:] * 4.0
    spimage.sp_image_free(ds)
    score = spimage.sp_find_center_refine_minimal_mask(c, 4, 0)
    x = pos_to_center(c.detector.image_center[0],img.shape[1])
    y = pos_to_center(c.detector.image_center[1],img.shape[0])
    spimage.sp_image_free(I)
    spimage.sp_image_free(kernel)
    spimage.sp_image_free(c)
    return (x,y,score)
    
