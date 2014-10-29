#ifndef _SP_SUPPORT_UPDATE_H_
#define _SP_SUPPORT_UPDATE_H_ 1

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

  typedef enum{SpSupportThreshold=1,SpSupportArea,SpSupportCentreImage,SpSupportTemplate,SpSupportStatic,SpSupportClose}SpSupportAlgorithmType;

/*! This structure is private */
typedef struct{
  sp_smap * blur_radius_map;
  sp_smap * threshold;
}SpSupportThresholdParameters;

/*! This structure is private */
typedef struct{
  sp_smap * blur_radius_map;
  sp_smap * area;
}SpSupportAreaParameters;

typedef struct{
  void * dummy;
}SpSupportCentreImageParameters;

/*! This structure is private */
typedef struct{
  //real blur_radius;
  Image *blured;
  Image *sorted;
  sp_smap *area;
  real original_area;
}SpSupportTemplateParameters;

typedef struct{
  void * dummy;
}SpSupportStaticParameters;

typedef struct{
  int size;
}SpSupportCloseParameters;

/*! This structure is private */
typedef struct{
  SpSupportAlgorithmType type;
  void * params;
  void *function;
  //int (*function)(void *); 
  //int (*function)(SpPhaser *);
}SpSupportAlgorithm;

typedef struct{
  int size;
  SpSupportAlgorithm **algorithms;
  int update_period;
}SpSupportArray;

  SpSupportAlgorithm * sp_support_threshold_alloc(sp_smap * blur_radius,sp_smap * threshold);
  SpSupportAlgorithm * sp_support_area_alloc(sp_smap * blur_radius,sp_smap * area);
  SpSupportAlgorithm * sp_support_template_alloc(Image *initial_support, real blur_radius, sp_smap *area);
  SpSupportAlgorithm * sp_support_static_alloc();
  SpSupportAlgorithm * sp_support_close_alloc(int size);
  SpSupportAlgorithm * sp_support_centre_image_alloc();
  
  SpSupportArray * sp_support_array_alloc(int size, int update_period);
  SpSupportArray * sp_support_array_init(SpSupportAlgorithm *algorithm, int update_period);
  void sp_support_array_set(SpSupportArray *array, int i, SpSupportAlgorithm *algorithm);
  void sp_support_array_append(SpSupportArray *array, SpSupportAlgorithm *);

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif
