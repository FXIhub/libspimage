#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include "spimage.h"


#define cutilSafeCall(err)           __cudaSafeCall      (err, __FILE__, __LINE__)
#define cufftSafeCall(err)           __cufftSafeCall     (err, __FILE__, __LINE__)

inline void __cudaSafeCall( cudaError err, const char *file, const int line )
{
    if( cudaSuccess != err) {
        fprintf(stderr, "cudaSafeCall() Runtime API error in file <%s>, line %i : %s.\n",
                file, line, cudaGetErrorString( err) );
        exit(-1);
    }
}

inline void __cufftSafeCall( cufftResult err, const char *file, const int line )
{
    if( CUFFT_SUCCESS != err) {
        fprintf(stderr, "cufftSafeCall() CUFFT error in file <%s>, line %i.\n",
                file, line);
        exit(-1);
    }
}

Image * sp_image_cuda_ifft(Image * img){
  cufftComplex *d_img;
  cufftHandle plan;
  int size = sp_image_size(img);
  Image * out = sp_image_duplicate(img,SP_COPY_DETECTOR);
  cutilSafeCall(cudaMalloc((void**)&d_img, sizeof(cufftComplex)*size));
  cutilSafeCall(cudaMemcpy(d_img, img->image->data, sizeof(cufftComplex)*size, cudaMemcpyHostToDevice));
  if(sp_image_z(img) == 1){
    cufftSafeCall(cufftPlan2d(&plan, sp_image_x(img),sp_image_y(img), CUFFT_C2C));
  }else{
    cufftSafeCall(cufftPlan3d(&plan, sp_image_x(img),sp_image_y(img),sp_image_z(img), CUFFT_C2C));
  }
  cufftSafeCall(cufftExecC2C(plan, d_img, d_img, CUFFT_INVERSE));
  cutilSafeCall(cudaMemcpy(out->image->data, d_img, sizeof(cufftComplex)*size, cudaMemcpyDeviceToHost));
  return out;
}

Image * sp_image_cuda_fft(Image * img){
  cufftComplex *d_img;
  cufftHandle plan;
  int size = sp_image_size(img);
  Image * out = sp_image_duplicate(img,SP_COPY_DETECTOR);
  cutilSafeCall(cudaMalloc((void**)&d_img, sizeof(cufftComplex)*size));
  cutilSafeCall(cudaMemcpy(d_img, img->image->data, sizeof(cufftComplex)*size, cudaMemcpyHostToDevice));
  if(sp_image_z(img) == 1){
    cufftSafeCall(cufftPlan2d(&plan, sp_image_x(img),sp_image_y(img), CUFFT_C2C));
  }else{
    cufftSafeCall(cufftPlan3d(&plan, sp_image_x(img),sp_image_y(img),sp_image_z(img), CUFFT_C2C));
  }
  cufftSafeCall(cufftExecC2C(plan, d_img, d_img, CUFFT_FORWARD));
  cutilSafeCall(cudaMemcpy(out->image->data, d_img, sizeof(cufftComplex)*size, cudaMemcpyDeviceToHost));
  return out;
}
