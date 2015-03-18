import numpy

def array_to_image(img,msk=None):
    s = img.shape
    d = len(list(s))
    if d == 3:
        sp_img = spimage.sp_image_alloc(s[2],s[1],s[0])
    else:
        sp_img = spimage.sp_image_alloc(s[2],s[1],1)
    sp_img.image[:] = img[:]
    if msk not None:
        sp_img.mask[:] = msk[:]
    return sp_img

def prtf(images_rs,supports,translate=True,enantio=True):
    S = images_rs.shape
    s = list(S)
    N = s.pop(0)
    s = tuple(s)

    image0_rs = images_rs[0]
    image0_fs = numpy.fftn(image0_rs)

    sp_image0_rs = array_to_image(image0_rs)
    sp_image0_fs = array_to_image(image0_fs)

    sp_amp_fs = spimage.sp_image_duplicate(sp_image0_fs,spimage.SP_COPY_ALL)
    spimage.sp_image_dephase(sp_amp_fs)

    spimage.sp_image_free(sp_image0_rs)
    spimage.sp_image_free(sp_image0_fs)

    sum_fs = image0_fs.copy()
    sum_fs[abs(sum_fs) > 0.] /= abs(sum_fs[abs(sum_fs) > 0.])

    sp_sum_fs = array_to_image(sum_fs)

    zeros = numpy.zeros(shape=s,dtype="int")
    zeros[abs(sum_fs) <= 0.] = 1

    sp_avg_img = array_to_image(image0_rs)
    avg_msk = numpy.zeros(shape=s,dtype="float")
    
    images_rs_super = numpy.zeros(shape=S,dtype="complex128")
    masks_rs_super = numpy.zeros(shape=S,dtype="bool")

    for i,img,sup in zip(range(1,N),images_rs[1:],supports[1:]):
        sp_img = array_to_image(img,sup)

        if translate:
            spimage.sp_image_superimpose(sp_avg_img,sp_img,spimage.SpEnantiomorph)
            spimage.sp_image_phase_match(sp_avg_img,sp_img,2)
        spimage.sp_image_add(sp_avg_img,sp_img)
        if sp_img.mask.sum() > 0:
            avg_msk[sp_img.mask] += 1
        
        images_rs_super[i,:] = sp_img.image[:]
        masks_rs_super[i,:] = sp_img.mask[:]
        
        sp_tmp = spimage.sp_image_fft(sp_img)
        sp_tmpamp = spimage.sp_image_duplicate(sp_tmp,spimage.SP_COPY_ALL)
        spimage.sp_image_dephase(sp_tmpamp);
        spimage.sp_image_add(sp_amp_fs,sp_tmpamp)
        
        positive = abs(sp_tmp.image) > 0.
        sp_tmp.image[positive] /= abs(sp_tmp.image)[positive]
        zeros += (positive == False)
        
        spimage.sp_image_add(sp_sum_fs,sp_tmp)
        
        spimage.sp_image_free(sp_img)
        spimage.sp_image_free(sp_tmp)
        spimage.sp_image_free(sp_tmpamp)
  
  sp_prtf = spimage.sp_image_duplicate(sp_sum_fs,spimage.SP_COPY_DATA|spimage.SP_COPY_MASK)
  sp_prtf.image[:] = N
  sp_prtf.image[zeros > 0] = 0.
  sp_image_dephase(prtf)

  avg_img = sp_avg_img.image[:]
  avg_sup = avg_msk > 0
  prtf = sp_prtf.image[:]

  for sp_i in [sp_prtf,sp_avg_img,sp_amp_fs,sp_sum_fs]:
      spimage.sp_image_free(sp_i)  

  out = {}
  out["prtf"] = prtf
  out["super_image"] = avg_img
  if full_out:
      out["prtf_r"] = radial_mean(prtf,cx=s[1]/2,cy=s[0]/2)
      out["super_mask"] = avg_sup
      out["images"] = images_rs_super
      out["masks"] = masks_rs_super

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
