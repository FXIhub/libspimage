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
        print("Error input")
        return []

def _radial(image,f=numpy.mean,shell_thickness=1.0,**kwargs):
    """
    Radial integration in N-dimensions. Assumes the input array has the same size in all dimensions. 
    Default integration method is the mean. 
    """
    n_dim = len(image.shape) 
    im_dim = image.shape[0] 
    im_center = image.shape[0]//2
    num_shells = int((im_dim-im_center)/shell_thickness) 

    c = numpy.array([im_center, im_center, im_center]) 
    
    r_out = numpy.arange(shell_thickness, (num_shells+1) * shell_thickness, shell_thickness) 
    r_in = r_out - 1 
    
    if n_dim == 2: 
        image = image[:,:,None] 
        c = numpy.array([im_center, im_center, 0]) 
    elif n_dim == 1:
        image = image[:,None,None] 
        c = numpy.array([im_center, 0, 0]) 

    x, y, z = numpy.meshgrid(numpy.arange(image.shape[0]), numpy.arange(image.shape[1]), numpy.arange(image.shape[2]), indexing='ij')
    radial_mask = (((x - c[0]) ** 2 + (y - c[1]) ** 2 + (z - c[2]) ** 2 >= r_in[:, None, None, None] ** 2) 
                   & ((x - c[0]) ** 2 + (y - c[1]) ** 2 + (z - c[2]) ** 2 < r_out[:, None, None, None] ** 2))

    prtf_r = numpy.array([f(numpy.squeeze(image[shell])) for shell in radial_mask]) 
    return prtf_r[numpy.isfinite(prtf_r)] 
def radial(image, **kwargs): 
    return _radial(image,**kwargs)

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
