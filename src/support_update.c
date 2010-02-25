#include <spimage.h>

static real bezier_map_interpolation(sp_smap * map, real x);
static void support_from_absolute_threshold(SpPhaser * ph, Image * blur, real abs_threshold);
static int descend_complex_compare(const void * pa,const void * pb);
int sp_support_static_update_cuda(void *ph);
int sp_support_area_update_cuda(void *ph);
int sp_support_threshold_update_cuda(void *ph);
int sp_support_template_update_cuda(void *ph);

SpSupportAlgorithm * sp_support_threshold_alloc(int update_period, sp_smap * blur_radius,sp_smap * threshold){
  SpSupportAlgorithm * ret = sp_malloc(sizeof(SpSupportAlgorithm));
  ret->type = SpSupportThreshold;
  ret->update_period = update_period;
  SpSupportThresholdParameters * params = sp_malloc(sizeof(SpSupportThresholdParameters));
  params->blur_radius_map = blur_radius;
  params->threshold = threshold;
  ret->params = params;
#ifdef _USE_CUDA
  ret->function = sp_support_threshold_update_cuda;
#else
  ret->function = sp_support_threshold_update;
#endif
  return ret;
}

SpSupportAlgorithm * sp_support_area_alloc(int update_period, sp_smap * blur_radius,sp_smap * area){
  SpSupportAlgorithm * ret = sp_malloc(sizeof(SpSupportAlgorithm));
  ret->type = SpSupportArea;
  ret->update_period = update_period;
  SpSupportAreaParameters * params = sp_malloc(sizeof(SpSupportAreaParameters));
  params->blur_radius_map = blur_radius;
  params->area = area;
  ret->params = params;
#ifdef _USE_CUDA
  ret->function = sp_support_area_update_cuda;
#else
  ret->function = sp_support_area_update;
#endif
  return ret;
}

SpSupportAlgorithm * sp_support_template_alloc(int update_period, Image *initial_support, real blur_radius, sp_smap *area){
  SpSupportAlgorithm *ret = sp_malloc(sizeof(SpSupportAlgorithm));
  ret->type = SpSupportTemplate;
  ret->update_period = update_period;
  SpSupportTemplateParameters * params = sp_malloc(sizeof(SpSupportTemplateParameters));
  params->blured = sp_gaussian_blur(initial_support,blur_radius);
  params->area = area;
  params->original_area = 0.0;
  //params->blur_radius = blur_radius;
  for (int i = 0; i < sp_image_size(initial_support); i++) {
    if (sp_real(initial_support->image->data[i]) != 0.0 ||
	sp_imag(initial_support->image->data[i]) != 0.0) {
      params->original_area += 1.0;
    }
  }
  params->original_area /= (real)sp_image_size(initial_support);
  //params->blured = gaussian_blur_sensitive(initial_support,blur_radius);
  params->sorted = sp_image_duplicate(params->blured,SP_COPY_DATA);
  qsort(params->sorted->image->data,sp_c3matrix_size(params->sorted->image),sizeof(Complex),descend_complex_compare);

  real max_threshold = sp_cabs(params->sorted->image->data[(int)(((real)sp_image_size(params->sorted))*sp_smap_max(params->area)*params->original_area)]);
  real min_threshold = sp_cabs(params->sorted->image->data[(int)(((real)sp_image_size(params->sorted))*sp_smap_min(params->area)*params->original_area)]);

  if (max_threshold < 1e-4) {
    fprintf(stderr,"Dangerously large template area with this threshold!\n");
  }
  if (1.0-min_threshold < 1e-4) {
    fprintf(stderr,"Dangerously small template area with this threshold!\n");
  }

#ifdef _USE_CUDA
  ret->function = sp_support_template_update_cuda;
#else
  ret->function = sp_support_template_update;
#endif

  ret->params = params;
  return ret;
}

SpSupportAlgorithm * sp_support_static_alloc(int update_period){
  SpSupportAlgorithm * ret = sp_malloc(sizeof(SpSupportAlgorithm));
  ret->type = SpSupportStatic;
  ret->update_period = update_period;
  SpSupportStaticParameters * params = sp_malloc(sizeof(SpSupportStaticParameters));
  ret->params = params;
#ifdef _USE_CUDA
  ret->function = sp_support_static_update_cuda;
#else
  ret->function = sp_support_static_update;
#endif
  return ret;
}

int sp_support_area_update_support(SpPhaser * ph){
  SpSupportAreaParameters * params = ph->sup_algorithm->params;
  real radius =  bezier_map_interpolation(params->blur_radius_map,ph->iteration);
  Image * tmp = sp_image_duplicate(ph->g1,SP_COPY_DATA);
  sp_image_dephase(tmp);
  Image * blur = sp_gaussian_blur(tmp, radius);
  real area = bezier_map_interpolation(params->area,ph->iteration);
  qsort(blur->image->data,sp_c3matrix_size(blur->image),sizeof(Complex),descend_complex_compare);
  real abs_threshold = sp_cabs(blur->image->data[(int)(sp_image_size(blur)*area)]);
  sp_image_free(blur);
  blur = sp_gaussian_blur(tmp, radius);
  sp_image_free(tmp);
  support_from_absolute_threshold(ph,blur,abs_threshold);
  sp_image_free(blur);
  return 0;
}

int sp_support_threshold_update_support(SpPhaser * ph){
  SpSupportThresholdParameters * params = ph->sup_algorithm->params;
  real radius =  bezier_map_interpolation(params->blur_radius_map,ph->iteration);
  Image * tmp = sp_image_duplicate(ph->g1,SP_COPY_DATA);
  sp_image_dephase(tmp);
  Image * blur = sp_gaussian_blur(tmp, radius);  
  sp_image_free(tmp);
  real rel_threshold = bezier_map_interpolation(params->threshold,ph->iteration);
  real abs_threshold = sp_image_max(blur,NULL,NULL,NULL,NULL)*rel_threshold;  
  support_from_absolute_threshold(ph,blur,abs_threshold);
  sp_image_free(blur);
  return 0;
}

int sp_support_template_update_support(SpPhaser * ph){
  SpSupportTemplateParameters * params = ph->sup_algorithm->params;
  //real radius = params->blur_radius;
  //Image *blur = gaussian_blur(params->template,params->blur_radius);
  //Image *support = sp_image_duplicate(params->blured,SP_COPY_DATA);
  real area = bezier_map_interpolation(params->area,ph->iteration);
  //qsort(blur->image->data,sp_c3matrix_size(blur->image),sizeof(Complex),descend_real_compare);
  //real abs_threshold = sp_cabs(blur->image->data[(int)(sp_image_size(blur)*area)]);area
  real abs_threshold = sp_cabs(params->sorted->image->data[(int)(sp_image_size(params->sorted)*area*params->original_area)]);
  //sp_image_free(blur);
  //support_from_absolute_threshold(ph,support,abs_threshold);
  support_from_absolute_threshold(ph,params->blured,abs_threshold);
  //sp_image_free(support);
  return 0;
}

int sp_support_static_update_support(SpPhaser * ph){
  return 0;
}


static void support_from_absolute_threshold(SpPhaser * ph, Image * blur, real abs_threshold){
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
				    5*pow(p3x,2) + 16*pow(x,2)),0.5),
		-0.3333333333333333) +
	    2*pow(p0x - p3x,-1)*
	    pow(8*p0x*p3x*x - 2*p3x*pow(p0x,2) - 4*x*pow(p0x,2) +
		2*pow(p0x,3) - 2*p0x*pow(p3x,2) - 4*x*pow(p3x,2) +
		2*pow(p3x,3) + pow(pow(p0x - p3x,4)*
				   (6*p0x*p3x - 16*p0x*x -
				    16*p3x*x + 5*pow(p0x,2) + 
				     5*pow(p3x,2) + 16*pow(x,2)),0.5)
		,0.3333333333333333)
	    )/4.;  
  return 3*p0y*t*pow(1 - t,2) + p0y*pow(1 - t,3) +
    3*p3y*(1 - t)*pow(t,2) + p3y*pow(t,3);
}


static int descend_complex_compare(const void * pa,const void * pb){
  Complex a,b;
  a = *((Complex *)pa);
  b = *((Complex *)pb);
  if(sp_cabs(a) < sp_cabs(b)){
    return 1;
  }else if(sp_cabs(a) == sp_cabs(b)){
    return 0;
  }else{
    return -1;
  }
}
