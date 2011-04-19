#ifndef _SP_PHASING_H_
#define _SP_PHASING_H_ 1

#include "image.h"
#include "cuda_util.h"
#include "support_update.h"
#include "map.h"

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

typedef enum{SpModelRandomPhases=1,SpModelZeroPhases=2,SpModelRandomValues=4,SpModelMaskedOutZeroed=256}SpModelInitialization;
typedef enum{SpSupportFromPatterson=1}SpSupportInitialization;
typedef enum{SpHIO=1,SpRAAR,SpDiffMap,SpER}SpPhasingAlgorithmType;
  typedef enum{SpNoConstraints=0,SpRealObject=1,SpPositiveRealObject=2,SpPositiveComplexObject=4,SpPositivityFlipping=8,SpCentrosymmetricObject=16}SpPhasingConstraints;
typedef enum{SpEngineAutomatic=0,SpEngineCPU=1,SpEngineCUDA=2}SpPhasingEngine;
typedef enum{SpPixelInsideSupport=1,SpPixelMeasuredAmplitude=2}SpPhasingPixelFlags;
typedef enum{SpRecoverPhases=0,SpRecoverAmplitudes=1}SpPhasingObjective;
/*! This structure is private */
typedef struct{
  sp_smap * beta;
  SpPhasingConstraints constraints;
}SpPhasingHIOParameters;

typedef struct{
  SpPhasingConstraints constraints;
}SpPhasingERParameters;


typedef SpPhasingHIOParameters SpPhasingRAARParameters;

/*! This structure is private */
typedef struct{
  sp_smap * beta;
  real gamma1;
  real gamma2;
  SpPhasingConstraints constraints;
}SpPhasingDiffMapParameters;

/*! This structure is private */
typedef struct{
  SpPhasingAlgorithmType type;
  void * params;
}SpPhasingAlgorithm;

/*! This structure is private */
typedef struct{
  /* amplitudes are used for phase recovery */
  sp_3matrix * amplitudes;
  /* phased_amplitudes are used for amplitude recovery */
  sp_c3matrix * phased_amplitudes;

  sp_i3matrix * pixel_flags;
  SpPhasingObjective phasing_objective;
  SpPhasingAlgorithm * algorithm;
  //SpSupportAlgorithm * sup_algorithm;
  SpSupportArray * sup_algorithm;
  int iteration;
  int image_size;
  int nx;
  int ny;
  int nz;

  /* These are images that are exposed to the user
   when sp_phaser_model() sp_phaser_model_change()
   and sp_phaser_support() are called */
  Image * model;
  int model_iteration;
  Image * old_model;
  int old_model_iteration;
  Image * model_change;
  int model_change_iteration;
  Image * support;
  int support_iteration;
  Image * amplitudes_image;
  Image * fmodel;
  int fmodel_iteration;

  SpPhasingEngine engine;

  Image * g0;
  Image * g1;

#ifdef _USE_CUDA
  cufftHandle cufft_plan;
  float * d_amplitudes;
  cufftComplex * d_phased_amplitudes;
  int * d_pixel_flags;
  cufftComplex * d_g0;
  cufftComplex * d_g1;
  int threads_per_block;
  int number_of_blocks;
#endif
}SpPhaser;


spimage_EXPORT SpPhasingAlgorithm * sp_phasing_hio_alloc(sp_smap * beta, SpPhasingConstraints constraints);
spimage_EXPORT SpPhasingAlgorithm * sp_phasing_raar_alloc(sp_smap * beta, SpPhasingConstraints constraints);
spimage_EXPORT SpPhasingAlgorithm * sp_phasing_diff_map_alloc(sp_smap * beta,real gamma1,real gamma2, SpPhasingConstraints constraints);
spimage_EXPORT SpPhasingAlgorithm * sp_phasing_er_alloc(SpPhasingConstraints constraints);

spimage_EXPORT SpPhaser * sp_phaser_alloc();
spimage_EXPORT void sp_phaser_free(SpPhaser * ph);

spimage_EXPORT const Image * sp_phaser_model(SpPhaser * ph);
spimage_EXPORT const Image * sp_phaser_model_with_support(SpPhaser * ph);
spimage_EXPORT const Image * sp_phaser_old_model(SpPhaser * ph);
  /*! Returns the fourier transform of the current best model */
spimage_EXPORT const Image * sp_phaser_fmodel(SpPhaser * ph);
spimage_EXPORT const Image * sp_phaser_fmodel_with_mask(SpPhaser * ph);
  spimage_EXPORT void sp_phaser_set_model(SpPhaser * ph,const Image * model);
  spimage_EXPORT void sp_phaser_set_support(SpPhaser * ph,const Image * support);
  spimage_EXPORT void sp_phaser_set_phased_amplitudes(SpPhaser * ph,const Image * phased_amplitudes);
  spimage_EXPORT void sp_phaser_set_amplitudes(SpPhaser * ph,const Image * amplitudes);
spimage_EXPORT Image * sp_phaser_model_change(SpPhaser * ph);
spimage_EXPORT const Image * sp_phaser_support(SpPhaser * ph);
spimage_EXPORT const Image * sp_phaser_amplitudes(SpPhaser * ph);
spimage_EXPORT int sp_phaser_init(SpPhaser * ph, SpPhasingAlgorithm * alg, SpSupportArray * sup_alg, SpPhasingEngine engine);
  spimage_EXPORT int sp_phaser_init_model(SpPhaser * ph,const Image * model, int flags);
  spimage_EXPORT int sp_phaser_init_support(SpPhaser * ph,const Image * support, int flags, real value);
spimage_EXPORT int sp_phaser_iterate(SpPhaser * ph, int iterations);
  spimage_EXPORT void sp_phaser_set_objective(SpPhaser * ph, SpPhasingObjective obj);
#ifdef _USE_CUDA
  int phaser_iterate_hio_cuda(SpPhaser * ph,int iterations);  
  int phaser_iterate_raar_cuda(SpPhaser * ph,int iterations);  
  int phaser_iterate_diff_map_cuda(SpPhaser * ph,int iterations);  
  int phaser_iterate_er_cuda(SpPhaser * ph,int iterations);  
  int sp_support_threshold_update_support_cuda(SpSupportAlgorithm *alg, SpPhaser * ph);
  int sp_support_area_update_support_cuda(SpSupportAlgorithm *alg, SpPhaser * ph);
  int sp_support_template_update_support_cuda(SpSupportAlgorithm *alg, SpPhaser * ph);
  int sp_support_static_update_support_cuda(SpSupportAlgorithm *alg, SpPhaser * ph);
  int sp_support_close_update_support_cuda(SpSupportAlgorithm *alg, SpPhaser * ph);
  int sp_support_centre_image_cuda(SpSupportAlgorithm *alg, SpPhaser * ph);
  int sp_proj_module_cuda(Image * a, Image * amp);
#endif

  int sp_support_area_update_support(SpSupportAlgorithm *alg, SpPhaser * ph);
  int sp_support_area_update_support_cpu(SpSupportAlgorithm *alg, SpPhaser * ph);
  int sp_support_threshold_update_support(SpSupportAlgorithm *alg, SpPhaser * ph);
  int sp_support_threshold_update_support_cpu(SpSupportAlgorithm *alg, SpPhaser * ph);
  int sp_support_template_update_support(SpSupportAlgorithm *alg, SpPhaser * ph);
  int sp_support_template_update_support_cpu(SpSupportAlgorithm *alg, SpPhaser * ph);
  int sp_support_static_update_support(SpSupportAlgorithm *alg, SpPhaser * ph);
  int sp_support_close_update_support(SpSupportAlgorithm *alg, SpPhaser * ph);
  int sp_support_close_update_support_cpu(SpSupportAlgorithm *alg, SpPhaser * ph);
  int sp_support_centre_image(SpSupportAlgorithm *alg, SpPhaser * ph);
  int sp_support_centre_image_cpu(SpSupportAlgorithm *alg, SpPhaser * ph);

  int sp_support_array_update(SpSupportArray *array, SpPhaser *ph);

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

  
#endif
