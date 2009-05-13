#ifndef _SP_PHASING_H_
#define _SP_PHASING_H_ 1

#include "image.h"

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

typedef struct{
  Image * amplitudes;
  Image * support;
  Image * prev_model;
  Image * model;
  Image * fmodel;
  Image * prev_fmodel;
  int algorithm;
  real beta;
  int support_size;
  int constraints;
}SpPhaser;

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

  
#endif
