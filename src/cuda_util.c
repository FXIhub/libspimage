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

#ifdef _USE_CUDA
Image * sp_get_image_from_cuda(cufftComplex * a, int size){
  Image * ret = sp_image_alloc(size,1,1);
  cutilSafeCall(cudaMemcpy(ret->image->data,a,sizeof(cufftComplex)*size,cudaMemcpyDeviceToHost));
  return ret;
}
#endif

sp_vector * sp_cuda_get_max_grid_size(){
#ifdef _USE_CUDA
  int dev_count;
  cudaGetDeviceCount(&dev_count);
  if(dev_count == 0){
    return NULL;
  }else{
    struct cudaDeviceProp dev_prop;
    cudaGetDeviceProperties(&dev_prop, 0);
    sp_vector * ret = sp_vector_alloc(3);
    sp_vector_set(ret,0,dev_prop.maxGridSize[0]);
    sp_vector_set(ret,1,dev_prop.maxGridSize[1]);
    sp_vector_set(ret,2,dev_prop.maxGridSize[2]);
    return ret;
  }
#endif
  return NULL;
}

