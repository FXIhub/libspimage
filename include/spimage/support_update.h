#ifndef _SP_SUPPORT_UPDATE_H_
#define _SP_SUPPORT_UPDATE_H_ 1

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

typedef enum{SpSupportFixed=1}SpSupportAlgorithmType;

/*! This structure is private */
typedef struct{
  int update_period;
  sp_smap * blur_radius_map;
  real threshold;
}SpSupportFixedParameters;

/*! This structure is private */
typedef struct{
  SpSupportAlgorithmType type;
  void * params;
}SpSupportAlgorithm;

SpSupportAlgorithm * sp_support_fixed_alloc(int update_period, sp_smap * blur_radius,real threshold);

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif
