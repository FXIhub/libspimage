#include <stdlib.h>
#include <math.h>
#include <hdf5.h>
#include <tiffio.h>
#include <png.h>
#include <float.h>
#include <ctype.h>
#ifdef _WIN32
#include <string.h>
#else
#include <strings.h>
#endif
#ifdef _USE_DMALLOC
#include <dmalloc.h>
#endif

#include "spimage.h"


Image * sp_image_noise_estimate(Image * intensities,Image * autocorrelation_support){
  Image * a = sp_image_duplicate(intensities,SP_COPY_ALL);
  sp_image_to_intensities(a);
  sp_image_rephase(a,SP_ZERO_PHASE);
  Image * ac = sp_image_ifft(a);
  sp_image_mul_elements(ac,autocorrelation_support);
  Image * filtered_intensities = sp_image_fft(ac);
  sp_image_scale(filtered_intensities,1.0/sp_image_size(ac));
  Image * std_dev = sp_image_alloc(sp_image_x(a),sp_image_y(a),sp_image_z(a));
  real arbitrary_constant = 1;
  for(int i = 0;i<sp_image_size(a);i++){
    /* Here lies a problem. The standard deviation is probably related to this difference, but I have no idea by what constant */
    sp_real(std_dev->image->data[i]) = arbitrary_constant*sp_real(filtered_intensities->image->data[i]) - sp_real(a->image->data[i]);
  }
  return std_dev;
}
