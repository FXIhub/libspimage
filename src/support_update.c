#include <spimage.h>

static real bezier_map_interpolation(sp_smap * map, real x);
static void support_from_absolute_threshold(SpPhaser * ph, Image * blur, real abs_threshold);
static int descend_complex_compare(const void * pa,const void * pb);

SpSupportAlgorithm * sp_support_threshold_alloc(sp_smap * blur_radius,sp_smap * threshold){
  SpSupportAlgorithm * ret = sp_malloc(sizeof(SpSupportAlgorithm));
  ret->type = SpSupportThreshold;
  SpSupportThresholdParameters * params = sp_malloc(sizeof(SpSupportThresholdParameters));
  params->blur_radius_map = blur_radius;
  params->threshold = threshold;
  ret->params = params;
  /*
#ifdef _USE_CUDA
  ret->function = sp_support_threshold_update_support_cuda;
#else
  ret->function = sp_support_threshold_update_support;
#endif
  */
  ret->function = sp_support_threshold_update_support;
  return ret;
}

SpSupportAlgorithm * sp_support_area_alloc(sp_smap * blur_radius,sp_smap * area){
  SpSupportAlgorithm * ret = sp_malloc(sizeof(SpSupportAlgorithm));
  ret->type = SpSupportArea;
  SpSupportAreaParameters * params = sp_malloc(sizeof(SpSupportAreaParameters));
  params->blur_radius_map = blur_radius;
  params->area = area;
  ret->params = params;
  /*
#ifdef _USE_CUDA
  ret->function = sp_support_area_update_support_cuda;
#else
  ret->function = sp_support_area_update_support;
#endif
  */
  ret->function = sp_support_area_update_support;
  return ret;
}

SpSupportAlgorithm * sp_support_template_alloc(Image *initial_support, real blur_radius, sp_smap *area){
  SpSupportAlgorithm *ret = sp_malloc(sizeof(SpSupportAlgorithm));
  ret->type = SpSupportTemplate;
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
  /*
#ifdef _USE_CUDA
  ret->function = sp_support_template_update_support_cuda;
#else
  ret->function = sp_support_template_update_support;
#endif
  */
  ret->function = sp_support_template_update_support;
  ret->params = params;
  return ret;
}

SpSupportAlgorithm * sp_support_static_alloc(){
  SpSupportAlgorithm * ret = sp_malloc(sizeof(SpSupportAlgorithm));
  ret->type = SpSupportStatic;
  SpSupportStaticParameters * params = sp_malloc(sizeof(SpSupportStaticParameters));
  ret->params = params;
  /*
#ifdef _USE_CUDA
  ret->function = sp_support_static_update_support_cuda;
#else
  ret->function = sp_support_static_update_support;
#endif
  */
  ret->function = sp_support_static_update_support;
  return ret;
}

SpSupportAlgorithm * sp_support_close_alloc(int size) {
  SpSupportAlgorithm *ret = sp_malloc(sizeof(SpSupportAlgorithm));
  ret->type = SpSupportClose;
  SpSupportCloseParameters * params = sp_malloc(sizeof(SpSupportCloseParameters));
  ret->params = params;
  params->size = size;
  /*
#ifdef _USE_CUDA
  ret->function = sp_support_close_update_support_cuda;
#else
  ret->function = sp_support_close_update_support;
#endif
  */
  ret->function = sp_support_close_update_support;
  return ret;
}

int sp_support_area_update_support(SpSupportAlgorithm *alg, SpPhaser * ph){
#ifdef _USE_CUDA
  if (ph->engine == SpEngineCUDA) {
    return sp_support_area_update_support_cuda(alg,ph);
  } else {
    return sp_support_area_update_support_cpu(alg,ph);
  }
#else
  return sp_support_area_update_support_cpu(alg,ph);
#endif
}

int sp_support_area_update_support_cpu(SpSupportAlgorithm *alg, SpPhaser * ph){
  SpSupportAreaParameters * params = alg->params;
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

int sp_support_threshold_update_support(SpSupportAlgorithm *alg, SpPhaser *ph){
#ifdef _USE_CUDA
  if (ph->engine == SpEngineCUDA) {
    return sp_support_threshold_update_support_cuda(alg,ph);
  } else {
    return sp_support_threshold_update_support_cpu(alg,ph);
  }
#else
  return sp_support_threshold_update_support_cpu(alg,ph);
#endif
}

int sp_support_threshold_update_support_cpu(SpSupportAlgorithm *alg, SpPhaser * ph){
  SpSupportThresholdParameters * params = alg->params;
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

int sp_support_template_update_support(SpSupportAlgorithm *alg, SpPhaser * ph){
#ifdef _USE_CUDA
  if (ph->engine == SpEngineCUDA) {
    return sp_support_template_update_support_cuda(alg,ph);
  } else {
    return sp_support_template_update_support_cpu(alg,ph);
  }
#else
  return sp_support_template_update_support_cpu(alg,ph);
#endif
}

int sp_support_template_update_support_cpu(SpSupportAlgorithm *alg, SpPhaser * ph){
  SpSupportTemplateParameters * params = alg->params;
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

int sp_support_static_update_support(SpSupportAlgorithm *alg, SpPhaser * ph){
  return 0;
}

int sp_support_close_update_support(SpSupportAlgorithm *alg, SpPhaser * ph){
#ifdef _USE_CUDA
  if (ph->engine == SpEngineCUDA) {
    return sp_support_close_update_support_cuda(alg,ph);
  } else {
    return sp_support_close_update_support_cpu(alg,ph);
  }
#else
  return sp_support_close_update_support_cpu(alg,ph);
#endif
}

int sp_support_close_update_support_cpu(SpSupportAlgorithm *alg, SpPhaser * ph){
  SpSupportCloseParameters *params = (SpSupportCloseParameters *)alg->params;
  int pixels = params->size;
  sp_i3matrix *tmp1 = sp_i3matrix_duplicate(ph->pixel_flags);
  sp_i3matrix *tmp2 = sp_i3matrix_duplicate(ph->pixel_flags);
  sp_i3matrix *foo;

  for (int i = 0; i < pixels; i++) {
    for (int x = 0; x < sp_i3matrix_x(ph->pixel_flags); x++) {
      for (int y = 0; y < sp_i3matrix_y(ph->pixel_flags); y++) {
	for (int z = 0; z < sp_i3matrix_z(ph->pixel_flags); z++) {
	  
	  if ((sp_i3matrix_get(tmp1,x,y,z) == SpPixelInsideSupport) ||
	      ((x != 0 && sp_i3matrix_get(tmp1,x-1,y,z) != SpPixelInsideSupport) ||
	       (x != sp_i3matrix_x(tmp1)-1 && sp_i3matrix_get(tmp1,x+1,y,z) != SpPixelInsideSupport) ||
	       (y != 0 && sp_i3matrix_get(tmp1,x,y-1,z) != SpPixelInsideSupport) ||
	       (y != sp_i3matrix_y(tmp1)-1 && sp_i3matrix_get(tmp1,x,y+1,z) != SpPixelInsideSupport) ||
	       (z != 0 && sp_i3matrix_get(tmp1,x,y,z-1) != SpPixelInsideSupport) ||
	       (z != sp_i3matrix_z(tmp1)-1 && sp_i3matrix_get(tmp1,x,y,z+1) != SpPixelInsideSupport))) {
	    tmp2->data[z*sp_i3matrix_y(tmp2)*sp_i3matrix_x(tmp2) +
		       y*sp_i3matrix_x(tmp2) + x] |= SpPixelInsideSupport;
	  } else {
	    tmp2->data[z*sp_i3matrix_y(tmp2)*sp_i3matrix_x(tmp2) +
		       y*sp_i3matrix_x(tmp2) + x] &= ~SpPixelInsideSupport;
	  }
	}
      }
    }
    foo = tmp1;
    tmp1 = tmp2;
    tmp2 = foo;
  }

  for (int i = 0; i < pixels; i++) {
    for (int x = 0; x < sp_i3matrix_x(ph->pixel_flags); x++) {
      for (int y = 0; y < sp_i3matrix_y(ph->pixel_flags); y++) {
	for (int z = 0; z < sp_i3matrix_z(ph->pixel_flags); z++) {
	  
	  if ((sp_i3matrix_get(tmp1,x,y,z) != SpPixelInsideSupport) ||
	      ((x != 0 && sp_i3matrix_get(tmp1,x-1,y,z) == SpPixelInsideSupport) ||
	       (x != sp_i3matrix_x(tmp1)-1 && sp_i3matrix_get(tmp1,x+1,y,z) == SpPixelInsideSupport) ||
	       (y != 0 && sp_i3matrix_get(tmp1,x,y-1,z) == SpPixelInsideSupport) ||
	       (y != sp_i3matrix_y(tmp1)-1 && sp_i3matrix_get(tmp1,x,y+1,z) == SpPixelInsideSupport) ||
	       (z != 0 && sp_i3matrix_get(tmp1,x,y,z-1) == SpPixelInsideSupport) ||
	       (z != sp_i3matrix_z(tmp1)-1 && sp_i3matrix_get(tmp1,x,y,z+1) == SpPixelInsideSupport))) {
	    tmp2->data[z*sp_i3matrix_y(tmp2)*sp_i3matrix_x(tmp2) +
		       y*sp_i3matrix_x(tmp2) + x] &= ~SpPixelInsideSupport;
	  } else {
	    tmp2->data[z*sp_i3matrix_y(tmp2)*sp_i3matrix_x(tmp2) +
		       y*sp_i3matrix_x(tmp2) + x] |= SpPixelInsideSupport;
	  }
	}
      }
    }
    foo = tmp1;
    tmp1 = tmp2;
    tmp2 = foo;
  }

  for (int i = 0; i < sp_i3matrix_size(ph->pixel_flags); i++) {
    ph->pixel_flags->data[i] = tmp1->data[i];
  }

  sp_i3matrix_free(tmp1);
  sp_i3matrix_free(tmp2);

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

/* Alloc an array with a certain size. After this the elemnets
   of the array needs to be set */
SpSupportArray * sp_support_array_alloc(int size, int update_period){
  SpSupportArray *ret = sp_malloc(sizeof(SpSupportArray));
  ret->size = size;
  ret->algorithms = sp_malloc(ret->size*sizeof(SpSupportAlgorithm *));
  ret->update_period = update_period;
  return ret;
}

/* Insert an elemnet in the array (doesn't change the size) */
void sp_support_array_set(SpSupportArray *array, int i, SpSupportAlgorithm *algorithm){
  array->algorithms[i] = algorithm;
}

/* Create an array with one element */
SpSupportArray * sp_support_array_init(SpSupportAlgorithm *algorithm, int update_period){
  SpSupportArray *ret = sp_support_array_alloc(1,update_period);
  sp_support_array_set(ret,0,algorithm);
  ret->update_period = update_period;
  return ret;
}

/* Append an element to the array */
void sp_support_array_append(SpSupportArray *array, SpSupportAlgorithm *algorithm){
  SpSupportAlgorithm **alg = sp_malloc((array->size+1)*sizeof(SpSupportAlgorithm *));
  memcpy(alg,array->algorithms,array->size*sizeof(SpSupportAlgorithm *));
  alg[array->size] = algorithm;
  sp_free(array->algorithms);
  array->algorithms = alg;
  array->size += 1;
}

int sp_support_array_update(SpSupportArray *array, SpPhaser *ph){
  for (int i = 0; i < array->size; i++) {
    ((int(*)(SpSupportAlgorithm *,SpPhaser *))array->algorithms[i]->function)(array->algorithms[i],ph);
  }
  return 0;
}
