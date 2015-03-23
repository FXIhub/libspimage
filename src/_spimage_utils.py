import numpy,scipy

# ================ #
# Helper functions #
# =================#
def _symmetrize(M,cx,cy):
    M_new = M.copy()
    M_new *= _turn180(M,cx,cy)
    return M_new
        
def _turn180(img,cx=None,cy=None):
    if cx == None:
        cx1 = (img.shape[0]-1)/2
    if cy == None:
        cy1 = (img.shape[0]-1)/2
    cx1 = round(cx*2)/2.
    cy1 = round(cy*2)/2.
    Nx1 = int(2*min([cx1,img.shape[1]-1-cx1]))+1
    Ny1 = int(2*min([cy1,img.shape[0]-1-cy1]))+1
    y_start = int(round(cy1-(Ny1-1)/2.))
    y_stop = int(round(cy1+(Ny1-1)/2.))+1
    x_start = int(round(cx1-(Nx1-1)/2.))
    x_stop = int(round(cx1+(Nx1-1)/2.))+1
    img_new = numpy.zeros(shape=(img.shape[0],img.shape[1]),dtype=img.dtype)
    img_new[y_start:y_stop,x_start:x_stop] = numpy.rot90(numpy.rot90(img[y_start:y_stop,x_start:x_stop]))
    return img_new

def _gaussian_smooth_2d1d(I,sm,precision=1.):
    N = 2*int(numpy.round(precision*sm))+1
    if len(I.shape) == 2:
        kernel = numpy.zeros(shape=(N,N))
        X,Y = numpy.meshgrid(numpy.arange(0,N,1),numpy.arange(0,N,1))
        X = X-N/2
        kernel = numpy.exp(X**2/(2.0*sm**2))
        kernel /= kernel.sum()
        Ism = scipy.signal.convolve2d(I,kernel,mode='same',boundary='wrap')
        return Ism
    elif len(I.shape) == 1:
        print "Error input"
        return []

def get_R_and_Theta_map(Nx,Ny,cx=None,cy=None):
    if not cx:
        cx = (Nx-1)/2.0
    if not cy:
        cy = (Ny-1)/2.0
    x = numpy.arange(0,Nx,1.0)-cx
    y = numpy.arange(0,Ny,1.0)-cy
    X,Y = numpy.meshgrid(x,y)
    R = numpy.sqrt(X**2+Y**2)
    R = R.round()
    Theta = numpy.arctan(-Y/(X+numpy.finfo('float64').eps))
    Theta[X<0] += numpy.pi
    Theta += numpy.pi/2.0
    #numpy.imsave("Theta.png" , Theta)
    #numpy.imsave("X.png" , X)
    #numpy.imsave("Y.png" , Y)
    return [R,Theta]

def _radial(image,mode="mean",**kwargs):
    if mode == "mean": f = numpy.mean
    elif mode == "sum": f = numpy.sum
    elif mode == "std": f = numpy.std
    elif mode == "median": f = numpy.median
    else:
        print "ERROR: No valid mode given for radial projection."
        return
    if 'cx' in kwargs: cx = kwargs['cx']
    else: cx = (image.shape[1]-1)/2.0
    if 'cy' in kwargs: cy = kwargs['cy'] 
    else: cy = (image.shape[0]-1)/2.0
    R = get_R_and_Theta_map(image.shape[1],image.shape[0],cx,cy)[0]
    R = R.round()
    R[numpy.isfinite(image)==False] = -1
    radii = numpy.arange(R.min(),R.max()+1,1)
    if radii[0] == -1:
        radii = radii[1:]
    values = numpy.zeros_like(radii)
    for i in range(0,len(radii)):
        values[i] = f(image[R==radii[i]])
    if 'rout' in kwargs: return numpy.array([radii,values])
    else: return values
def radial_sum(image,**kwargs):
    return _radial(image,"sum",**kwargs)
def radial_std(image,**kwargs):
    return _radial(image,"std",**kwargs)
def radial_mean(image,**kwargs):
    return _radial(image,"mean",**kwargs)
def radial_median(image,**kwargs):
    return _radial(image,"median",**kwargs)

def cone_pixel_average(image,N_theta,cx=None,cy=None):
    [R,Theta] = get_R_and_Theta_map(image.shape[1],image.shape[0],cx,cy)
    R[numpy.isfinite(image) == False] = -1
    radii = numpy.arange(R.min(),R.max()+1,1)
    if radii[0] == -1:
        radii = radii[1:]
    values = numpy.zeros(shape=(N_theta,len(radii)))
    for j in range(0,N_theta):
        theta_min = j/(1.0*N_theta)*2.0*numpy.pi
        theta_max = (j+1)/(1.0*N_theta)*2.0*numpy.pi
        theta_center = (theta_max+theta_min)/2.0
        theta_interval = theta_max-theta_min
        theta_image = image[abs(Theta-theta_center)<=theta_interval/2.0]
        theta_R = R[abs(Theta-theta_center)<=theta_interval/2.0]
        for i in range(0,len(radii)):
            temp = theta_image[theta_R==radii[i]].copy()
            temp = temp[numpy.isfinite(temp)]
            values[j,i] = temp.mean()
    return [radii,values]
