import os,re,sys,h5py,numpy,time,cmath
import spimage

# Set up logger
import logging
logger = logging.getLogger('PRTF')
logger.setLevel("WARNING")

# Toggle here debug mode for testing the functions
debug = False
this_folder = "."#os.path.dirname(os.path.realpath(__file__))
if debug:
    import matplotlib

def prtf(imgs0,msks0,**kwargs):
    """
    This function calculates the phase retrieval transfer function from a stack of real-space images and its supports.
    """

    out = {}

    logger0 = kwargs.get("logger",logger)
    #K = numpy.random.randint(1000)
    K = 0
    enantio = kwargs.get("enantio",True)
    shifted = kwargs.get("shifted",True)
    center_result = kwargs.get("center_result",None)
    pixels_to_exclude = kwargs.get("pixels_to_exclude",None)
    do_maximize_overlap = kwargs.get("maximize_overlap",True)
    do_minimize_phase_ramp = kwargs.get("minimize_phase_ramp",False)
    do_phase_match = kwargs.get("real_space_phase_match",True)
    do_align_com_support = kwargs.get("align_com_support",False)

    if logger0 != None:
        s = "  "
        s += "- "
        if not enantio: "no "
        s += "enantio\n"
        s += "- "
        if not shifted: "not "
        s += "shifted\n"
        s += "- "
        if pixels_to_exclude == None:
            s += "no pixels specified to exclude\n"
        else:
            s += "%i pixels specified to exclude\n" % pixel_to_exclude.sum()
        s += "- "
        if not do_maximize_overlap: s += "not "
        s += "maximizing overlap\n"
        s += "- "
        if not do_minimize_phase_ramp: s += "not "
        s += "minimizing phase ramp\n"
        s += "- "
        if not do_phase_match: s += "no "
        s += "phase matching\n"
        s += "- "
        if not do_align_com_support: s += "not "
        s += "aligning center of mass of support\n"
        s += "- "
        if center_result == None: 
            s += "no centering of the resulting image\n"
        else:
            s += "centering result in mode %s\n" % center_result
        logger0.info("PRTF runs with the folloing configuration:\n %s" % s)

    if debug:
        os.system("rm %s/testdata/prtf*" % this_folder)

    Nx = imgs0.shape[2]
    Ny = imgs0.shape[1]
    cx = kwargs.get("cx",(Nx-1)/2.)
    cy = kwargs.get("cy",(Ny-1)/2.)
    selection = kwargs.get("selection",numpy.ones(imgs0.shape[0],dtype="bool"))
    N = selection.sum()
    imgs = numpy.zeros(shape=(N,imgs0.shape[1],imgs0.shape[2]),dtype=imgs0.dtype)
    msks = numpy.zeros(shape=(N,msks0.shape[1],msks0.shape[2]),dtype="float")
    k = 0
    for i in range(imgs0.shape[0]):
        if selection[i]:
            if shifted:
                imgs[k,:,:] = numpy.fft.fftshift(imgs0[i,:,:])
                msks[k,:,:] = numpy.fft.fftshift(msks0[i,:,:])
            else:
                imgs[k,:,:] = imgs0[i,:,:]
                msks[k,:,:] = msks0[i,:,:]
            k += 1
    if debug:
        out["dbg_imgs"] = imgs.copy()
        out["dbg_msks"] = msks.copy()

    # Average reconstructions
    # superimpose for obtaining the averaged result of the reconstruction
    imgs1 = numpy.zeros(shape=(N,Ny,Nx),dtype=imgs0.dtype)
    msks1 = numpy.zeros(shape=(N,Ny,Nx),dtype="float")
    for i in range(0,N):
        img1 = imgs[i,:,:].copy()
        msk1 = msks[i,:,:].copy()
        img0 = imgs1[0,:,:].copy()
        msk0 = msks1[0,:,:].copy()

        if debug:
            j=0
            matplotlib.pyplot.imsave(this_folder+"/testdata/prtf_%i_%i_I_%i.png" % (K,i,j),abs(img1),vmin=0.,vmax=3.)
            matplotlib.pyplot.imsave(this_folder+"/testdata/prtf_%i_%i_IP_%i.png" % (K,i,j),numpy.angle(img1) % (2.*numpy.pi),vmin=0.,vmax=2.*numpy.pi)
            matplotlib.pyplot.imsave(this_folder+"/testdata/prtf_%i_%i_M_%i.png" % (K,i,j),abs(msk1))
        if do_minimize_phase_ramp:
            [img1,translation] = minimize_phase_ramp(img1,shifted=False,periodic_boundary=True)
            msk1 = numpy.int16(abs(fourier_translation(msk1,translation)).round())
        
        
        if debug:
            j+=1
            matplotlib.pyplot.imsave(this_folder+"/testdata/prtf_%i_%i_I_%i.png" % (K,i,j),abs(img1),vmin=0.,vmax=3.)
            matplotlib.pyplot.imsave(this_folder+"/testdata/prtf_%i_%i_IP_%i.png" % (K,i,j),numpy.angle(img1) % (2.*numpy.pi),vmin=0.,vmax=2.*numpy.pi)
            matplotlib.pyplot.imsave(this_folder+"/testdata/prtf_%i_%i_A_M_%i.png" % (K,i,j),abs(msk1))
        if do_maximize_overlap and i!=0:
            [img1,translation,turned] = maximize_overlap(img0,img1,enantio)
            if turned: msk1 = fft_turn180(msk1)
            msk1 = abs(fourier_translation(msk1,translation))

        if debug:
            j+=1
            matplotlib.pyplot.imsave(this_folder+"/testdata/prtf_%i_%i_I_%i.png" % (K,i,j),abs(img1),vmin=0.,vmax=3.)
            matplotlib.pyplot.imsave(this_folder+"/testdata/prtf_%i_%i_IP_%i.png" % (K,i,j),numpy.angle(img1) % (2.*numpy.pi),vmin=0.,vmax=2.*numpy.pi)
            matplotlib.pyplot.imsave(this_folder+"/testdata/prtf_%i_%i_A_M_%i.png" % (K,i,j),abs(msk1))
        if do_phase_match and i!=0:
            weights = abs(img0)*abs(img1)
            img1 = abs(img1)*numpy.exp(1.j*(numpy.angle(img1)+phase_match(img0,img1,weights)))

        if debug:
            j+=1
            matplotlib.pyplot.imsave(this_folder+"/testdata/prtf_%i_%i_I_%i.png" % (K,i,j),abs(img1),vmin=0.,vmax=3.)
            matplotlib.pyplot.imsave(this_folder+"/testdata/prtf_%i_%i_IP_%i.png" % (K,i,j),numpy.angle(img1) % (2.*numpy.pi),vmin=0.,vmax=2.*numpy.pi)
            matplotlib.pyplot.imsave(this_folder+"/testdata/prtf_%i_%i_A_M_%i.png" % (K,i,j),abs(msk1))
            print "Power: %f" % (abs(img1).sum()/(abs(img0).sum()+numpy.finfo("float32").eps))
            print "Avg. phase: %f,%f" % ((msk0*numpy.angle(img1)).mean(),(msk0*numpy.angle(img0)).mean())
            print "Diff: %f" % (abs(img1-img0).mean()/abs(img0).mean())
        imgs1[i,:,:] = img1[:,:]
        msks1[i,:,:] = numpy.int16(abs(msk1).round())[:,:]
    imgs1_super = imgs1.mean(0)
    msks1_super = msks1.mean(0)
    # Make PRTF
    # go to fourier space
    fimgs = numpy.zeros_like(imgs)
    fimgs1 = numpy.zeros_like(imgs)
    for i in range(N):
        fimgs[i,:,:] = numpy.fft.fftn(imgs[i,:,:])
        fimgs1[i,:,:] = numpy.fft.fftn(imgs1[i,:,:])
    # mask zeros
    PRTF = numpy.zeros_like(imgs)
    tmp = abs(fimgs1) != 0.
    if pixels_to_exclude != None:
        tmp *=  pixels_to_exclude
    PRTF[tmp] = fimgs1[tmp]/abs(fimgs1[tmp])
    PRTF = abs(PRTF.mean(0))
    PRTF[(fimgs == 0).sum(0) != 0] = 0.
    PRTF = numpy.array(PRTF,dtype="float32")

    if debug:
        matplotlib.pyplot.imsave(this_folder+"/testdata/superI%i.png" % numpy.random.randint(1000),abs(imgs1_super),vmin=0,vmax=2.)
        matplotlib.pyplot.imsave(this_folder+"/testdata/superM%i.png" % numpy.random.randint(1000),abs(msks1_super),vmin=0,vmax=1.)

    if center_result != None:
        if center_result == "image":
            CM = center_of_mass(abs(imgs1_super))
        elif center_result == "support_times_image":
            CM = center_of_mass(msks1_super*abs(imgs1_super))
        elif center_result == "support":
            CM = center_of_mass(msks1_super)
        imgs1_super = pixel_translation(imgs1_super,CM,3)
        msks1_super = pixel_translation(msks1_super,CM,3)

    if do_align_com_support:
        com = center_of_mass(msks1_super)
        if com[0] > 0 and com[1] > 0:
            imgs1_super = fft_turn180(imgs1_super)
            msks1_super = abs(fft_turn180(msks1_super))

    if shifted:
        imgs1_super = numpy.fft.fftshift(imgs1_super)
        msks1_super = numpy.fft.fftshift(msks1_super)
        for i in range(N):
            fimgs1[i,:,:] = numpy.fft.fftshift(fimgs1[i,:,:])
            imgs1[i,:,:] = numpy.fft.fftshift(imgs1[i,:,:])
            msks1[i,:,:] = numpy.fft.fftshift(msks1[i,:,:])
        PRTF = numpy.fft.fftshift(PRTF)

    msks1_super = numpy.int16(msks1_super)
    msks1 = numpy.int16(msks1)
    out["prtf"] = PRTF
    out["super_image"] = imgs1_super
    out["super_mask"] = msks1_super
    out["images"] = imgs1
    out["masks"] = msks1
    out["fourier_images"] = fimgs1
    return out


def recover_translation(imgA,imgB,enantio=False):
    imgfA = numpy.fft.fftn(imgA)
    imgfB = numpy.fft.fftn(imgB)
    d = len(list(imgfA.shape))
    f = numpy.indices(imgfA.shape)
    for i in range(d):
        f[i] = f[i]-numpy.ceil(imgfA.shape[i]/2.)
        f[i,:] = numpy.fft.fftshift(f[i,:])[:]
    if enantio == False:
        # Check superposition with image
        imgB_new = imgB
        c = abs(numpy.fft.ifftn(imgfA*imgfB.conj()))
        turned = False
    else:
        imgB_turned = fft_turn180(imgB)
        imgfB_turned = numpy.fft.fftn(imgB_turned)
        # Check superposition with normal and rotated image
        cc = [abs(numpy.fft.ifftn(imgfA*imgfB.conj())),abs(numpy.fft.ifftn(imgfA*imgfB_turned.conj()))]
        if debug:
            os.system("rm %s/testdata/RT*" % this_folder)
            matplotlib.pyplot.imsave(this_folder+"/testdata/RT_imgB_turned.png",abs(numpy.fft.fftshift(imgB_turned)),vmin=0.,vmax=abs(imgB).max())
            matplotlib.pyplot.imsave(this_folder+"/testdata/RT_imgB.png",abs(numpy.fft.fftshift(imgB)),vmin=0.,vmax=abs(imgB).max())
            matplotlib.pyplot.imsave(this_folder+"/testdata/RT_imgA.png",abs(numpy.fft.fftshift(imgA)),vmin=0.,vmax=abs(imgB).max())
            matplotlib.pyplot.imsave(this_folder+"/testdata/RT_imgfB_turned.png",abs(numpy.fft.fftshift(imgfB_turned)),vmin=0.,vmax=abs(imgfB).max())
            matplotlib.pyplot.imsave(this_folder+"/testdata/RT_imgfB.png",abs(numpy.fft.fftshift(imgfB)),vmin=0.,vmax=abs(imgfB).max())
            matplotlib.pyplot.imsave(this_folder+"/testdata/RT_imgfA.png",abs(numpy.fft.fftshift(imgfA)),vmin=0.,vmax=abs(imgfB).max())
            matplotlib.pyplot.imsave(this_folder+"/testdata/RT_CC0.png",cc[0])
            matplotlib.pyplot.imsave(this_folder+"/testdata/RT_CC1.png",cc[1])
        Mcc = numpy.array([cc[0].max(),cc[1].max()])
        i_max = Mcc.argmax()
        if i_max == 0:
            turned = False
            c = cc[0]
        else:
            turned = True
            c = cc[1]
    index_max = c.argmax()
    translation = []
    for i in range(d):
        translation.append(f[i,:].flatten()[index_max])
    translation = numpy.array(translation)
    return [translation,turned]


def fourier_translation(A,t,rotation=False):
    fA = numpy.fft.fftn(A)
    d = len(list(fA.shape))
    f = numpy.indices(fA.shape)
    for i in range(d):
        f[i] = f[i]-numpy.ceil(fA.shape[i]/2.)
        f[i,:] = numpy.fft.fftshift(f[i,:])[:]
    tmp = 0
    for i,ti,fi in zip(range(d),t,f):
        tmp = tmp + 2*numpy.pi*fi[:,:]*ti/fA.shape[i]
    A_translated = numpy.fft.ifftn(fA*numpy.exp(-1.j*tmp))
    #print "%e" % (abs(A).sum()/abs(A_translated).sum())
    return A_translated

def pixel_translation(A,t,order=1):
    if A.dtype == "complex64" or A.dtype == "complex128":
        return (1.*pixel_translation(A.real,t,order)+1.j*pixel_translation(A.imag,t,order))
    from scipy import ndimage
    d = len(list(A.shape))
    g = numpy.indices(A.shape)
    gt = numpy.indices(A.shape)
    for i in range(d):
        gt[i] = (gt[i]+t[i]) % A.shape[i]
    return ndimage.map_coordinates(A, gt, order=order)


# Minimize the difference between the phases of a and b by adding a constant phase to b.
# (transferred into python from libspimage)
def phase_match(imgA,imgB,weights=None): # typically weights = (abs(imgA)*abs(imgB))
    diff = numpy.angle(imgA)-numpy.angle(imgB)
    if weights == None:
        w = 1/(1.*len(diff.flatten()) + numpy.finfo('float64').eps)
    else:
        w = weights / (weights.sum() + numpy.finfo('float64').eps)
    return (diff*w).sum()


# This functions translates image b so that it's phases 
# are as close as possible to a.
# The translation is done in fourier space and both images
# should be in real space
# (transferred into python from libspimage)
def maximize_overlap(imgA0,imgB0,enantio=False):
    imgA = imgA0.copy()
    imgB = imgB0.copy()
    [translation,turned] = recover_translation(imgA,imgB,enantio)
    if turned: imgB = fft_turn180(imgB)
    imgB = fourier_translation(imgB,translation)
    return [imgB,translation,turned]

# in fourier space subtract phase ramp obtained by leastsq, corresponding to translation in real space
# input: real space
# output: real space and fourier space data
def minimize_phase_ramp(img,shifted=False,periodic_boundary=False):
    from scipy.optimize import leastsq,fmin_cobyla
    imgf = numpy.fft.fftn(numpy.fft.fftshift(img))
    pimgf = numpy.angle(imgf)+numpy.pi

    d = len(list(imgf.shape))
    f = numpy.indices(imgf.shape)
    for i in range(d):
        f[i] = f[i]-numpy.ceil(imgf.shape[i]/2.)
        f[i,:] = numpy.fft.fftshift(f[i,:])[:]
    if d == 1:
        p = lambda v: (v[0] + v[1]*f) % (2*numpy.pi)
        v00 = pimgf[f[0]==0][0]
        v01 = pimgf[f[0]==1][0]-v00
        v0 = [v00,v01]
    elif d ==2:
        p = lambda v: (v[0] + v[1]*f[0,:,:] + v[2]*f[1,:,:])  % (2*numpy.pi)
        v00 = pimgf[(f[0]==0)*(f[1]==0)][0]
        v01 = pimgf[(f[0]==1)*(f[1]==0)][0]-v00
        v02 = pimgf[(f[0]==0)*(f[1]==1)][0]-v00
        v0 = [v00,v01,v02]
    elif d ==3:
        p = lambda v: v[0] + (v[1]*f[0,:,:] + v[2]*f[1,:,:] + v[3]*f[2,:,:])  % (2*numpy.pi)
        v00 = pimgf[(f[0]==0)*(f[1]==0)*(f[2]==0)][0]
        v01 = pimgf[(f[0]==1)*(f[1]==0)*(f[2]==0)][0]-v00
        v02 = pimgf[(f[0]==0)*(f[1]==1)*(f[2]==0)][0]-v00
        v03 = pimgf[(f[0]==0)*(f[1]==0)*(f[2]==1)][0]-v00
        v0 = [v00,v01,v02,v03]
    err = lambda v: ((pimgf-p(v))**2).sum()
    v1, success = leastsq(lambda v: numpy.ones(len(v))*err(v),v0)

    v2 = v1.copy()
    if periodic_boundary:
        for i in range(d):
            m1 = v1[i+1]*imgf.shape[i]/(2.*numpy.pi)
            m2 = round(v1[i+1]*imgf.shape[i]/(2.*numpy.pi))
            v2[i+1] = m2/m1*v1[i+1]

    resf = abs(imgf)*numpy.exp(1.j*(pimgf-p(v2)))
    res = numpy.fft.ifftn(resf)

    if shifted:
        resf = numpy.fft.fftshift(resf)
        res = numpy.fft.fftshift(res)

    translation = (v2[1:]-numpy.pi)/(2*numpy.pi)*numpy.array(pimgf.shape)
    if debug:
        matplotlib.pyplot.imsave(this_folder+"/testdata/subtract_phase_ramp_img0.png",abs(img))
        matplotlib.pyplot.imsave(this_folder+"/testdata/subtract_phase_ramp_img1.png",abs(res))
        matplotlib.pyplot.imsave(this_folder+"/testdata/subtract_phase_ramp_img0t.png",abs(fourier_translation(img,translation)))
    return [res,translation]

def center_of_mass(img0,shifted=False):
    img = abs(img0)
    img = img/(1.*img.sum()+numpy.finfo("float32").eps)
    d = len(list(img.shape))
    cm = numpy.zeros(d)
    f = numpy.indices(img.shape)
    for i in range(d):
        f[i] = f[i]-numpy.ceil(img.shape[i]/2.)
        if not shifted:
            f[i,:] = numpy.fft.fftshift(f[i,:])[:]
        cm[i] = (f[i,:]*img[:]).sum()
        if debug:
            matplotlib.pyplot.imsave(this_folder+"/testdata/f%i.png" % i,f[i])
    return cm

def fft_turn180(I):
    return numpy.fft.fftn(numpy.fft.fftn(I))/(1.*I.size)
