#include <spimage.h>
#include <stdint.h>

#ifndef NDEBUG
#define sp_cuda_check_errors() __sp_cuda_check_errors(__FILE__,__LINE__)
#else
#define sp_cuda_check_errors() 
#endif

SpCUDADeviceType sp_cuda_get_device_type(){
#ifdef _USE_CUDA
  int dev_count = 0;
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

#ifdef _USE_CUDA
Image * sp_get_image_from_cuda(cufftComplex * a, int size){
  Image * ret = sp_image_alloc(size,1,1);
  cutilSafeCall(cudaMemcpy(ret->image->data,a,sizeof(cufftComplex)*size,cudaMemcpyDeviceToHost));
  return ret;
}
#endif

void sp_cuda_launch_parameters(int64_t size, int * gridSize, int * blockSize){
#ifdef _USE_CUDA
  int64_t max_blockSize, max_gridSize;
  struct cudaDeviceProp dev_prop;
  cudaGetDeviceProperties(&dev_prop, 0);
  if(dev_prop.major < 3){
    max_blockSize = 512;
    max_gridSize = 65535;
  }else{
    max_blockSize = 1024;
    max_gridSize = 2147483647;
  }
  *blockSize = 64;
  while(size > max_gridSize* (*blockSize)){
    *blockSize *= 2;
  }
  if(*blockSize > max_blockSize){
    sp_error_fatal("Image too large to handle with current CUDA compute capability %d.%d\n. Please use a graphics card with compute capability 3.0 or greater.", dev_prop.major, dev_prop.minor);
  }
  *gridSize = (size+*blockSize-1)/(*blockSize);
#endif
}
