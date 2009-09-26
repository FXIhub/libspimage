#ifndef _SP_IMAGE_FILTER_CUDA_H_
#define _SP_IMAGE_FILTER_CUDA_H_ 1

#ifdef _USE_CUDA

#include <cuda_runtime_api.h>
#include <cufft.h>

#ifdef __cplusplus
extern "C"
{
#endif
spimage_EXPORT  void sp_create_gaussian_kernel_cuda(cufftComplex * a,int x, int y, int z, float radius);
  spimage_EXPORT  void sp_gaussian_blur_cuda(cufftComplex * in, cufftComplex * out,int x, int y, int z, float radius, cufftHandle plan);

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif /* CUDA */
#endif /* _SP_IMAGE_FILTER_CUDA_H_ */
