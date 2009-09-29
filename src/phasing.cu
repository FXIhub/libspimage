#include <spimage.h>

__global__ void CUDA_module_projection(cufftComplex* g, const float* amp,const int * pixel_flags, const  int size);
__global__ void CUDA_support_projection_hio(cufftComplex* g1, const cufftComplex* g0,const int * pixel_flags,const  int size,const float beta);
__global__ void CUDA_complex_scale(cufftComplex * a, int size ,float scale);
__global__ void CUDA_support_projection_raar(cufftComplex* g1, const cufftComplex* g0,const int * pixel_flags,const  int size,const float beta);
__global__ void CUDA_apply_constraints(cufftComplex* g, const int * pixel_flags,const  int size,const SpPhasingConstraints constraints);

int sp_proj_module_cuda(Image * a, Image * amp){
  cufftComplex * d_a;
  int * d_pixel_flags;
  float * d_amp;
  cutilSafeCall(cudaMalloc((void **)&d_a,sizeof(cufftComplex)*sp_image_size(a)));
  cutilSafeCall(cudaMalloc((void **)&d_pixel_flags,sizeof(int)*sp_image_size(a)));
  cutilSafeCall(cudaMalloc((void **)&d_amp,sizeof(float)*sp_image_size(a)));
  cutilSafeCall(cudaMemcpy(d_a,a->image->data,sizeof(cufftComplex)*sp_image_size(a),cudaMemcpyHostToDevice));
  sp_i3matrix * pixel_flags = sp_i3matrix_alloc(sp_image_x(a),sp_image_y(a),sp_image_z(a));
  sp_3matrix * h_amp = sp_3matrix_alloc(sp_image_x(a),sp_image_y(a),sp_image_z(a));
  for(int i =0 ;i<sp_image_size(a);i++){
    h_amp->data[i] = sp_real(amp->image->data[i]);
    pixel_flags->data[i] = 0;
    if(amp->mask->data[i]){
      pixel_flags->data[i] |= SpPixelMeasuredAmplitude;
    }
  }
  cutilSafeCall(cudaMemcpy(d_amp,h_amp->data,sizeof(float)*sp_image_size(a),cudaMemcpyHostToDevice));
  cutilSafeCall(cudaMemcpy(d_pixel_flags,pixel_flags->data,sizeof(int)*sp_image_size(a),cudaMemcpyHostToDevice));
  int threads_per_block = 64;
  int number_of_blocks = (sp_image_size(a)+threads_per_block-1)/threads_per_block;
  CUDA_module_projection<<<number_of_blocks, threads_per_block>>>(d_a,d_amp,d_pixel_flags,sp_image_size(a));
  cutilSafeCall(cudaMemcpy(a->image->data,d_a,sizeof(cufftComplex)*sp_image_size(a),cudaMemcpyDeviceToHost));
  cutilSafeCall(cudaFree(d_amp));
  cutilSafeCall(cudaFree(d_pixel_flags));
  cutilSafeCall(cudaFree(d_a));
  return 0;
}

int phaser_iterate_hio_cuda(SpPhaser * ph,int iterations){
  SpPhasingHIOParameters * params = (SpPhasingHIOParameters *)ph->algorithm->params;
  const real beta = params->beta;
  for(int i = 0;i<iterations;i++){
    cufftComplex * swap = ph->d_g0;
    ph->d_g0 = ph->d_g1;
    ph->d_g1 = swap;
    /* executes FFT processes */
    cufftSafeCall(cufftExecC2C(ph->cufft_plan, ph->d_g0, ph->d_g1, CUFFT_FORWARD));
    CUDA_module_projection<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_g1,ph->d_amplitudes,ph->d_pixel_flags,ph->image_size);
    sp_cuda_check_errors();
    cufftSafeCall(cufftExecC2C(ph->cufft_plan, ph->d_g1, ph->d_g1, CUFFT_INVERSE));
    /* normalize */
    CUDA_complex_scale<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_g1,ph->image_size, 1.0f / (ph->image_size));
    sp_cuda_check_errors();
    CUDA_support_projection_hio<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_g1,ph->d_g0,ph->d_pixel_flags,ph->image_size,beta);
    sp_cuda_check_errors();
  }
  CUDA_apply_constraints<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_g1,ph->d_pixel_flags,ph->image_size,params->constraints);
  sp_cuda_check_errors();
  ph->iteration += iterations;
  return 0;
}

int phaser_iterate_raar_cuda(SpPhaser * ph,int iterations){
  SpPhasingRAARParameters * params = (SpPhasingRAARParameters *)ph->algorithm->params;
  const real beta = params->beta;
  for(int i = 0;i<iterations;i++){
    cufftComplex * swap = ph->d_g0;
    ph->d_g0 = ph->d_g1;
    ph->d_g1 = swap;
    /* executes FFT processes */
    cufftSafeCall(cufftExecC2C(ph->cufft_plan, ph->d_g0, ph->d_g1, CUFFT_FORWARD));
    CUDA_module_projection<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_g1,ph->d_amplitudes,ph->d_pixel_flags,ph->image_size);
    sp_cuda_check_errors();
    cufftSafeCall(cufftExecC2C(ph->cufft_plan, ph->d_g1, ph->d_g1, CUFFT_INVERSE));
    /* normalize */
    CUDA_complex_scale<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_g1,ph->image_size, 1.0f / (ph->image_size));
    sp_cuda_check_errors();
    CUDA_support_projection_raar<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_g1,ph->d_g0,ph->d_pixel_flags,ph->image_size,beta);
    sp_cuda_check_errors();
  }
  CUDA_apply_constraints<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_g1,ph->d_pixel_flags,ph->image_size,params->constraints);
  sp_cuda_check_errors();
  ph->iteration += iterations;
  return 0;
}
