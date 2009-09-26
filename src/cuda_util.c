#include <spimage.h>

#ifndef NDEBUG
#define sp_cuda_check_errors() __sp_cuda_check_errors(__FILE__,__LINE__)
#else
#define sp_cuda_check_errors() 
#endif

SpCUDADeviceType sp_cuda_get_device_type(){
#ifdef _USE_CUDA
  int dev_count;
  cudaGetDeviceCount(&dev_count);
  if(dev_count == 0){
    return SpCUDANoDevice;
  }else{
    struct cudaDeviceProp dev_prop;
    cudaGetDeviceProperties(&dev_prop, 0);
    if(dev_prop.major == 9999 && dev_prop.minor == 9999){
      /* emulated device */
      return SpCUDAEmulatedDevice;
    }else{
      return SpCUDAHardwareDevice;
    }
  }
#else
  return SpCUDANoDevice;
#endif
}