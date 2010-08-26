#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include "spimage.h"


Image * sp_image_cuda_ifft(const Image * img){
  cufftComplex *d_img;
  cufftHandle plan;
  int size = sp_image_size(img);
  Image * out = sp_image_duplicate(img,SP_COPY_DETECTOR);
  cutilSafeCall(cudaMalloc((void**)&d_img, sizeof(cufftComplex)*size));
  cutilSafeCall(cudaMemcpy(d_img, img->image->data, sizeof(cufftComplex)*size, cudaMemcpyHostToDevice));
  if(sp_image_z(img) == 1){
    cufftSafeCall(cufftPlan2d(&plan, sp_image_y(img),sp_image_x(img), CUFFT_C2C));
  }else{
    cufftSafeCall(cufftPlan3d(&plan, sp_image_z(img),sp_image_y(img),sp_image_x(img), CUFFT_C2C));
  }
  cufftSafeCall(cufftExecC2C(plan, d_img, d_img, CUFFT_INVERSE));
  cutilSafeCall(cudaMemcpy(out->image->data, d_img, sizeof(cufftComplex)*size, cudaMemcpyDeviceToHost));
  cufftSafeCall(cufftDestroy(plan));
  cutilSafeCall(cudaFree(d_img));
  return out;
}

Image * sp_image_cuda_fft(const Image * img){
  cufftComplex *d_img;
  cufftHandle plan;
  int size = sp_image_size(img);
  Image * out = sp_image_duplicate(img,SP_COPY_DETECTOR);
  cutilSafeCall(cudaMalloc((void**)&d_img, sizeof(cufftComplex)*size));
  cutilSafeCall(cudaMemcpy(d_img, img->image->data, sizeof(cufftComplex)*size, cudaMemcpyHostToDevice));
  if(sp_image_z(img) == 1){
    cufftSafeCall(cufftPlan2d(&plan, sp_image_y(img),sp_image_x(img), CUFFT_C2C));
  }else{
    cufftSafeCall(cufftPlan3d(&plan, sp_image_z(img),sp_image_y(img),sp_image_x(img), CUFFT_C2C));
  }
  cufftSafeCall(cufftExecC2C(plan, d_img, d_img, CUFFT_FORWARD));
  cutilSafeCall(cudaMemcpy(out->image->data, d_img, sizeof(cufftComplex)*size, cudaMemcpyDeviceToHost));
  cutilSafeCall(cudaFree(d_img));
  cufftSafeCall(cufftDestroy(plan));
  return out;
}
