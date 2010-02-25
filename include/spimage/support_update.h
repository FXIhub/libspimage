#ifndef _SP_SUPPORT_UPDATE_H_
#define _SP_SUPPORT_UPDATE_H_ 1

#include "phasing.h"

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

  typedef enum{SpSupportThreshold=1,SpSupportArea,SpSupportTemplate,SpSupportStatic}SpSupportAlgorithmType;

/*! This structure is private */
typedef struct{
  sp_smap * blur_radius_map;
  sp_smap * threshold;
  int close;
}SpSupportThresholdParameters;

/*! This structure is private */
typedef struct{
  sp_smap * blur_radius_map;
  sp_smap * area;
  int close;
}SpSupportAreaParameters;

/*! This structure is private */
typedef struct{
  //real blur_radius;
  Image *blured;
  Image *sorted;
  sp_smap *area;
  real original_area;
}SpSupportTemplateParameters;

typedef struct{
}SpSupportStaticParameters;

/*! This structure is private */
typedef struct{
  SpSupportAlgorithmType type;
  int update_period;
  void * params;
  int (*function)(void *); 
}SpSupportAlgorithm;

  SpSupportAlgorithm * sp_support_threshold_alloc(int update_period, sp_smap * blur_radius,sp_smap * threshold);
  SpSupportAlgorithm * sp_support_area_alloc(int update_period, sp_smap * blur_radius,sp_smap * area);
  SpSupportAlgorithm * sp_support_template_alloc(int update_period, Image *initial_support, real blur_radius, sp_smap *area);
  SpSupportAlgorithm * sp_support_static_alloc(int update_period);

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif
