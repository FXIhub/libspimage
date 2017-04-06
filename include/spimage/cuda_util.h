#ifndef _SP_CUDA_UTIL_H_
#define _SP_CUDA_UTIL_H_ 1

/** @defgroup CudaUtils CUDA Utils
 *  This module helps to interface with CUDA
 *  It's mostly for internal use.
 *  @{
 */

#ifdef _USE_CUDA

#include <cuda_runtime_api.h>
#include <cufft.h>
#include <spimage/sperror.h>



#ifndef NDEBUG
#define sp_cuda_check_errors() __sp_cuda_check_errors(__FILE__,__LINE__)
#else
#define sp_cuda_check_errors() 
#endif

#define cufftSafeCall(err)           __cufftSafeCall     (err, __FILE__, __LINE__)
static inline void __cufftSafeCall( cufftResult err, const char *file, const int line )
{
    if( CUFFT_SUCCESS != err) {
        fprintf(stderr, "cufftSafeCall() CUFFT error in file <%s>, line %i.\n",
                file, line);
        exit(-1);
    }
}
#define cutilSafeCall(err)           __cudaSafeCall      (err, __FILE__, __LINE__)
static inline void __cudaSafeCall( cudaError_t err, const char *file, const int line )
{
    if( cudaSuccess != err) {
        fprintf(stderr, "cudaSafeCall() Runtime API error in file <%s>, line %i : %s.\n",
                file, line, cudaGetErrorString( err) );
        exit(-1);
    }
}



#define cutilCheckMsg(msg)           __cutilCheckMsg     (msg, __FILE__, __LINE__)




static inline void __cutilCheckMsg( const char *errorMessage, const char *file, const int line )
{   
    cudaError_t err = cudaGetLastError();
    if( cudaSuccess != err) {
        fprintf(stderr, "cutilCheckMsg() CUTIL CUDA error: %s in file <%s>, line %i : %s.\n",
                errorMessage, file, line, cudaGetErrorString( err) );
        exit(-1);
    }
#ifdef _DEBUG
    err = cudaThreadSynchronize();
    if( cudaSuccess != err) {
        fprintf(stderr, "cutilCheckMsg cudaThreadSynchronize error: %s in file <%s>, line %i : %s.\n",
                errorMessage, file, line, cudaGetErrorString( err) );
        exit(-1);
    }
#endif
}


static inline void __sp_cuda_check_errors(const char * file, const int line){
  cudaError_t error = cudaGetLastError();
  if(error != cudaSuccess ){
    fprintf(stderr,"CUDA error \"%s\" at %s:%d\n",cudaGetErrorString(error),file,line);
    abort();
  }
}

/*! Creates an image from a cufftComplex.

   The output image will be 1 dimensional.
 */
spimage_EXPORT Image * sp_get_image_from_cuda(cufftComplex * a, int size);


/*! Calculates the gridSize and blockSize to handle an image of a given size
 */
spimage_EXPORT void sp_cuda_launch_parameters(int size, int * gridSize, int * blockSize);
#endif /* _USE_CUDA */

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

  /*! Types of CUDA devices.
   */
  typedef enum{SpCUDANoDevice=0,SpCUDAEmulatedDevice=1,SpCUDAHardwareDevice=2}SpCUDADeviceType;

  /*! Returns the CUDA device currently being used if any.
   */
spimage_EXPORT  SpCUDADeviceType sp_cuda_get_device_type();

/*@}*/
#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif /* _SP_CUDA_UTIL_H_ */
