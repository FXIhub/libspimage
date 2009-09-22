#ifndef _SP_PHASING_H_
#define _SP_PHASING_H_ 1

#include "image.h"
#include "cuda_util.h"
#include "support_update.h"

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

typedef enum{SpModelRandomPhases=1,SpModelZeroPhases=2,SpModelRandomValues=4,SpModelMaskedOutZeroed=256}SpModelInitialization;
typedef enum{SpHIO=1,SpRAAR,SpDiffMap}SpPhasingAlgorithmType;
typedef enum{SpNoConstraints=0,SpRealObject=1,SpPositiveRealObject=2,SpPositiveComplexObject=4,SpPositivityFlipping=8}SpPhasingConstraints;
typedef enum{SpOutputNothing=0,SpOutputModel=1,SpOutputModelChange=2}SpPhasingOutput;
typedef enum{SpEngineAutomatic=0,SpEngineCPU=1,SpEngineCUDA=2}SpPhasingEngine;
typedef enum{SpPixelInsideSupport=1,SpPixelMeasuredAmplitude=2}SpPhasingPixelFlags;
/*! This structure is private */
typedef struct{
  real beta;
  SpPhasingConstraints constraints;
}SpPhasingHIOParameters;


typedef SpPhasingHIOParameters SpPhasingRAARParameters;

/*! This structure is private */
typedef struct{
  real beta;
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
  sp_3matrix * amplitudes;
  sp_i3matrix * pixel_flags;
  SpPhasingAlgorithm * algorithm;
  SpSupportAlgorithm * sup_algorithm;
  int iteration;
  int image_size;

  Image * model;
  int model_iteration;
  Image * model_change;
  int model_change_iteration;

  SpPhasingEngine engine;

  Image * g0;
  Image * g1;

#ifdef _USE_CUDA
  cufftHandle cufft_plan;
  cudaStream_t calc_stream;
  cudaStream_t transfer_stream;
  float * d_amplitudes;
  int * d_pixel_flags;
  cufftComplex * d_g0;
  cufftComplex * d_g1;
  cufftComplex * d_g0_transfer;
  cufftComplex * d_g1_transfer;
  int threads_per_block;
  int number_of_blocks;
#endif
}SpPhaser;


spimage_EXPORT SpPhasingAlgorithm * sp_phasing_hio_alloc(real beta, SpPhasingConstraints constraints);
spimage_EXPORT SpPhasingAlgorithm * sp_phasing_raar_alloc(real beta, SpPhasingConstraints constraints);

spimage_EXPORT SpPhaser * sp_phaser_alloc();
spimage_EXPORT void sp_phaser_free(SpPhaser * ph);

spimage_EXPORT Image * sp_phaser_model(const SpPhaser * ph, int * iteration);
spimage_EXPORT Image * sp_phaser_model_change(const SpPhaser * ph, int * iteration);
spimage_EXPORT int sp_phaser_init(SpPhaser * ph, SpPhasingAlgorithm * alg,Image * amplitudes,Image * support, SpPhasingEngine engine);
spimage_EXPORT int sp_phaser_init_model(SpPhaser * ph,const Image * model, int flags);
spimage_EXPORT int sp_phaser_iterate(SpPhaser * ph, int iterations, SpPhasingOutput output);

#ifdef _USE_CUDA
int phaser_iterate_hio_cuda(SpPhaser * ph,int iterations, SpPhasingOutput output);  
int phaser_iterate_raar_cuda(SpPhaser * ph,int iterations, SpPhasingOutput output);  
#endif

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

  
#endif
