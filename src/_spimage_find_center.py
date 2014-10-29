import numpy as np
import scipy as sp
import spimage


def find_center(img, msk, method=None, **kwargs):
    """
    Find the center of a diffraction pattern.

    usage:

    x,y = find_center(img, msk)
    x,y = find_center(img, msk, method='quadrant',  x0=0, y0=0, dmax=None, solver='L-BFGS-B')
    x,y = find_center(img, msk, method='pixelwise', x0=0, y0=0, dmax=5, rmax=None, downsampling=1)
    """

    # Default method for center finding
    if method is None: method = 'quadrant'

    # Find center using "quadrant" method
    if method == 'quadrant':
        x,y = find_center_quadrant(img, msk, **kwargs)
    # Find center using "pixelwise" method
    elif method == 'pixelwise':
        x,y = find_center_pixelwise(img, msk, **kwargs)
    # Return 0,0 if method is not defined
    else:
        x,y = (0,0)
        print "There is no center finding method %s" %method
    return x,y


def find_center_quadrant(img, msk, x0=0, y0=0, dmax=None, solver='L-BFGS-B'):
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

        def manipulate(self, threshold=13.):
            self.image[self.image < threshold] = 1.
            self.image = np.log(self.image)
            
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
    m.manipulate()
    m.start()
    x = m.res["x"][0]
    y = m.res["x"][1]
    return x,y

def find_center_pixelwise(img, msk, x0, y0, dmax=5, rmax=None, downsampling=1.):
    """
    Find center of diffraction pattern using a pixelwise comparison of centry-symmetric pixels.
    """
    print dmax, rmax, downsampling
    #rmax = None
    #if rmax is not None: msk &= (spimage.rgrid(msk.shape, (x0,y0)) < rmax)
    I = spimage.sp_image_alloc(np.ceil(img.shape[1]/float(downsampling)).astype(int), np.ceil(img.shape[0]/float(downsampling)).astype(int), 1)
    I.image[:] = img[::downsampling,::downsampling]
    I.mask[:]  = msk[::downsampling,::downsampling]
    I.detector.image_center[:] = np.array([x0/float(downsampling) + img.shape[1]/2, y0/float(downsampling) + img.shape[0]/2, 0 ])
    success = spimage.sp_find_center_refine(I, dmax, 0, None)
    print success, I.detector.image_center[0], I.detector.image_center[1]
    x = (I.detector.image_center[0]  - img.shape[1]/2 ) * downsampling
    y = (I.detector.image_center[1]  - img.shape[0]/2 ) * downsampling
    return x,y
