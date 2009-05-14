#ifndef _SP_PHASING_H_
#define _SP_PHASING_H_ 1

#include "image.h"

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

typedef enum{SpModelRandomPhases=1,SpModelZeroPhases=2,SpModelRandomValues=4,SpModelMaskedOutZeroed=256}SpModelInitialization;
typedef enum{SpHIO=1,SpRAAR,SpDiffMap}SpPhasingAlgorithmType;
  typedef enum{SpNoConstraints=0,SpRealObject=1,SpPositiveRealObject=2,SpPositiveComplexObject=4,SpPositivityFlipping=8}SpPhasingConstraints;

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
  Image * amplitudes;
  Image * support;
  Image * model;
  Image * model_change;
  SpPhasingAlgorithm * algorithm;
  int iteration;
}SpPhaser;


SpPhasingAlgorithm * sp_phasing_hio_alloc(real beta, SpPhasingConstraints constraints);
SpPhasingAlgorithm * sp_phasing_raar_alloc(real beta, SpPhasingConstraints constraints);

SpPhaser * sp_phaser_alloc();
void sp_phaser_free(SpPhaser * ph);

spimage_EXPORT Image * sp_phaser_model(const SpPhaser * ph);
spimage_EXPORT Image * sp_phaser_model_change(const SpPhaser * ph);
spimage_EXPORT int sp_phaser_init(SpPhaser * ph, SpPhasingAlgorithm * alg,Image * amplitudes,Image * support);
spimage_EXPORT int sp_phaser_init_model(SpPhaser * ph,const Image * model, int flags);
spimage_EXPORT int sp_phaser_iterate(SpPhaser * ph);
#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

  
#endif
