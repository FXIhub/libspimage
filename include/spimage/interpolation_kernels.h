#ifndef _INTERPOLATION_KERNELS_H_
#define _INTERPOLATION_KERNELS_H_ 1

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "mem_util.h"


#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

/** @defgroup Interpolation Interpolation
 *  This modules provides image interpolation
 *  @{
 */


typedef enum{
  Kernel_Real_R2 = 1,
  Kernel_Real_R = 2,
}SP_Kernel_Type;

typedef struct _sp_kernel{
  real * table;
  /* The size of the table */
  int size;
  /* The boundary of the kernel given in r^2 or r depending on the kernel type */
  real boundary;
  real boundary_in_3D;
  SP_Kernel_Type type;  
}sp_kernel;

  /*! 
    \param max_r2 the maximum radius squared of the kernel
    \param tolerance the maximum error allowed
    \return the created kernel
   */
spimage_EXPORT sp_kernel * sp_create_spline2_kernel_table(real max_r2,real tolerance);
spimage_EXPORT sp_kernel * sp_create_spline2_r_kernel_table(real max_r2,int size);
  /*!
    \param k the kernel to sample from
    \param pos the position where the kernel will be sampled
    \return the value of the kernel at the given position
   */
spimage_EXPORT real sp_kernel_table_sample(sp_kernel * k, real pos);
  /*!
    \param k the kernel to free
   */
spimage_EXPORT void sp_kernel_free(sp_kernel * k);


/*@}*/

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

  
#endif
