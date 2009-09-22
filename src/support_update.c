#include <spimage.h>


SpSupportAlgorithm * sp_support_fixed_alloc(int update_period, sp_smap * blur_radius,real threshold){
  SpSupportAlgorithm * ret = sp_malloc(sizeof(SpSupportAlgorithm));
  ret->type = SpSupportFixed;
  SpSupportFixedParameters * params = sp_malloc(sizeof(SpSupportFixedParameters));
  params->update_period = update_period;
  params->blur_radius_map = blur_radius;
  params->threshold = threshold;
  ret->params = params;
  return ret;
}

