#include <spimage.h>

__global__ void CUDA_module_projection(cufftComplex* g, const float* amp,const int * pixel_flags, const  int size);
__global__ void CUDA_support_projection_hio(cufftComplex* g1, const cufftComplex* g0,const int * pixel_flags,const  int size,const float beta);
__global__ void CUDA_complex_scale(cufftComplex * a, int size ,float scale);
__global__ void CUDA_support_projection_raar(cufftComplex* g1, const cufftComplex* g0,const int * pixel_flags,const  int size,const float beta);
__global__ void CUDA_apply_constraints(cufftComplex* g, const int * pixel_flags,const  int size,const SpPhasingConstraints constraints);

int phaser_iterate_hio_cuda(SpPhaser * ph,int iterations){
  SpPhasingRAARParameters * params = (SpPhasingRAARParameters *)ph->algorithm->params;
  const real beta = params->beta;
  for(int i = 0;i<iterations;i++){
    cufftComplex * swap = ph->d_g0;
    ph->d_g0 = ph->d_g1;
    ph->d_g1 = swap;
    /* executes FFT processes */
    cufftSafeCall(cufftExecC2C(ph->cufft_plan, ph->d_g0, ph->d_g1, CUFFT_FORWARD));
    CUDA_module_projection<<<ph->number_of_blocks, ph->threads_per_block,0,ph->calc_stream>>>(ph->d_g1,ph->d_amplitudes,ph->d_pixel_flags,ph->image_size);
    sp_cuda_check_errors();
    cufftSafeCall(cufftExecC2C(ph->cufft_plan, ph->d_g1, ph->d_g1, CUFFT_INVERSE));
    /* normalize */
    CUDA_complex_scale<<<ph->number_of_blocks, ph->threads_per_block,0,ph->calc_stream>>>(ph->d_g1,ph->image_size, 1.0f / (ph->image_size));
    sp_cuda_check_errors();
    CUDA_support_projection_hio<<<ph->number_of_blocks, ph->threads_per_block,0,ph->calc_stream>>>(ph->d_g1,ph->d_g0,ph->d_pixel_flags,ph->image_size,beta);
    sp_cuda_check_errors();
  }
  CUDA_apply_constraints<<<ph->number_of_blocks, ph->threads_per_block,0,ph->calc_stream>>>(ph->d_g1,ph->d_pixel_flags,ph->image_size,params->constraints);
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
    CUDA_module_projection<<<ph->number_of_blocks, ph->threads_per_block,0,ph->calc_stream>>>(ph->d_g1,ph->d_amplitudes,ph->d_pixel_flags,ph->image_size);
    sp_cuda_check_errors();
    cufftSafeCall(cufftExecC2C(ph->cufft_plan, ph->d_g1, ph->d_g1, CUFFT_INVERSE));
    /* normalize */
    CUDA_complex_scale<<<ph->number_of_blocks, ph->threads_per_block,0,ph->calc_stream>>>(ph->d_g1,ph->image_size, 1.0f / (ph->image_size));
    sp_cuda_check_errors();
    CUDA_support_projection_raar<<<ph->number_of_blocks, ph->threads_per_block,0,ph->calc_stream>>>(ph->d_g1,ph->d_g0,ph->d_pixel_flags,ph->image_size,beta);
    sp_cuda_check_errors();
  }
  CUDA_apply_constraints<<<ph->number_of_blocks, ph->threads_per_block,0,ph->calc_stream>>>(ph->d_g1,ph->d_pixel_flags,ph->image_size,params->constraints);
  sp_cuda_check_errors();
  ph->iteration += iterations;
  return 0;
}
