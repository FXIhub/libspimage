#ifndef _SP_CUDA_UTIL_H_
#define _SP_CUDA_UTIL_H_ 1

#ifdef _USE_CUDA

#include <cuda_runtime_api.h>
#include <cufft.h>


#define cutilSafeCall(err)           __cudaSafeCall      (err, __FILE__, __LINE__)
#define cufftSafeCall(err)           __cufftSafeCall     (err, __FILE__, __LINE__)
#define cutilCheckMsg(msg)           __cutilCheckMsg     (msg, __FILE__, __LINE__)

static inline void __cudaSafeCall( cudaError_t err, const char *file, const int line )
{
    if( cudaSuccess != err) {
        fprintf(stderr, "cudaSafeCall() Runtime API error in file <%s>, line %i : %s.\n",
                file, line, cudaGetErrorString( err) );
        exit(-1);
    }
}

static inline void __cufftSafeCall( cufftResult err, const char *file, const int line )
{
    if( CUFFT_SUCCESS != err) {
        fprintf(stderr, "cufftSafeCall() CUFFT error in file <%s>, line %i.\n",
                file, line);
        exit(-1);
    }
}


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

#endif /* _USE_CUDA */

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */
  typedef enum{SpCUDANoDevice=0,SpCUDAEmulatedDevice=1,SpCUDAHardwareDevice=2}SpCUDADeviceType;

spimage_EXPORT  SpCUDADeviceType sp_cuda_get_device_type();

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif /* _SP_CUDA_UTIL_H_ */
