#include <spimage.h>


SpSupportAlgorithm * sp_support_threshold_alloc(int update_period, sp_smap * blur_radius,sp_smap * threshold){
  SpSupportAlgorithm * ret = sp_malloc(sizeof(SpSupportAlgorithm));
  ret->type = SpSupportThreshold;
  ret->update_period = update_period;
  SpSupportThresholdParameters * params = sp_malloc(sizeof(SpSupportThresholdParameters));
  params->blur_radius_map = blur_radius;
  params->threshold = threshold;
  ret->params = params;
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
  return ret;
}

int sp_support_area_update_support(SpPhaser * ph){
  return 0;
}

int sp_support_threshold_update_support(SpPhaser * ph){
  return 0;
}
