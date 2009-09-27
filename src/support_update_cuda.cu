#include <spimage.h>

static real bezier_map_interpolation(sp_smap * map, real x);
static void support_from_absolute_threshold_cuda(SpPhaser * ph, Image * blur, real abs_threshold);

int sp_support_threshold_update_support_cuda(SpPhaser * ph){
  SpSupportThresholdParameters * params = (SpSupportThresholdParameters *)ph->sup_algorithm->params;
  real radius =  bezier_map_interpolation(params->blur_radius_map,ph->iteration);
  Image * blur = gaussian_blur(ph->g1, radius);
  real rel_threshold = bezier_map_interpolation(params->threshold,ph->iteration);
  real abs_threshold = sp_image_max(blur,NULL,NULL,NULL,NULL)*rel_threshold;  
  support_from_absolute_threshold_cuda(ph,blur,abs_threshold);
  sp_image_free(blur);
  return 0;
}

static void support_from_absolute_threshold_cuda(SpPhaser * ph, Image * blur, real abs_threshold){
  for(int i =0 ;i<ph->image_size;i++){
    if(sp_cabs(sp_image_get_by_index(blur,i)) > abs_threshold){
      ph->pixel_flags->data[i] |= SpPixelInsideSupport;
    }else{
      ph->pixel_flags->data[i] &= ~SpPixelInsideSupport;
    }
  }
}


static real bezier_map_interpolation(sp_smap * map, real x){
  sp_list * keys = sp_smap_get_keys(map);
  sp_list * values = sp_smap_get_values(map);
  unsigned int idx = 0;
  for(idx = 0;idx<sp_list_size(keys);idx++){
    if(x < sp_list_get(keys,idx)){
      break;
    }
  }
  if(idx == 0){
    return sp_list_get(values,0);
  }
  if(idx == sp_list_size(keys)){
    return sp_list_get(values,sp_list_size(keys)-1);
  }
  /* Cubic Bezier curve taken from http://en.wikipedia.org/wiki/Bézier_curve */
  real p0y = sp_list_get(values,idx-1);
  real p0x = sp_list_get(keys,idx-1);
  real p3y = sp_list_get(values,idx);
  real p3x = sp_list_get(keys,idx);
  real t = (2 - 2*(p0x - p3x)*
	    pow(8*p0x*p3x*x - 2*p3x*pow(p0x,2) - 4*x*pow(p0x,2) + 
		2*pow(p0x,3) - 2*p0x*pow(p3x,2) - 4*x*pow(p3x,2) + 
		2*pow(p3x,3) + pow(pow(p0x - p3x,4)*
				   (6*p0x*p3x - 16*p0x*x - 16*p3x*x + 5*pow(p0x,2) +
				    5*pow(p3x,2) + 16*pow(x,2)),0.5f),
		-0.3333333333333333f) +
	    2*pow(p0x - p3x,-1)*
	    pow(8*p0x*p3x*x - 2*p3x*pow(p0x,2) - 4*x*pow(p0x,2) +
		2*pow(p0x,3) - 2*p0x*pow(p3x,2) - 4*x*pow(p3x,2) +
		2*pow(p3x,3) + pow(pow(p0x - p3x,4)*
				   (6*p0x*p3x - 16*p0x*x -
				    16*p3x*x + 5*pow(p0x,2) + 
				     5*pow(p3x,2) + 16*pow(x,2)),0.5f)
		,0.3333333333333333f)
	    )/4.;  
  return 3*p0y*t*pow(1 - t,2) + p0y*pow(1 - t,3) +
    3*p3y*(1 - t)*pow(t,2) + p3y*pow(t,3);
}
