#include <spimage.h>

static void phaser_cuda_handle_output_begin(SpPhaser * ph,SpPhasingOutput output);
static void phaser_cuda_handle_output_end(SpPhaser * ph,SpPhasingOutput output);

__global__ void CUDA_module_projection(cufftComplex* g, const float* amp,const int * pixel_flags, const  int size);
__global__ void CUDA_support_projection_hio(cufftComplex* g1, const cufftComplex* g0,const int * pixel_flags,const  int size,const float beta);
__global__ void CUDA_complex_scale(cufftComplex * a, int size ,float scale);
__global__ void CUDA_support_projection_raar(cufftComplex* g1, const cufftComplex* g0,const int * pixel_flags,const  int size,const float beta);
__global__ void CUDA_apply_constraints(cufftComplex* g, const int * pixel_flags,const  int size,const SpPhasingConstraints constraints);

int phaser_iterate_hio_cuda(SpPhaser * ph,int iterations, SpPhasingOutput output){
  SpPhasingRAARParameters * params = (SpPhasingRAARParameters *)ph->algorithm->params;
  const real beta = params->beta;
  cudaStreamSynchronize(ph->calc_stream);
  phaser_cuda_handle_output_begin(ph,output);
  for(int i = 0;i<iterations;i++){
    cufftComplex * swap = ph->d_g0;
    ph->d_g0 = ph->d_g1;
    ph->d_g1 = swap;
    /* executes FFT processes */
    cufftSafeCall(cufftExecC2C(ph->cufft_plan, ph->d_g0, ph->d_g1, CUFFT_FORWARD));
    CUDA_module_projection<<<ph->number_of_blocks, ph->threads_per_block,0,ph->calc_stream>>>(ph->d_g1,ph->d_amplitudes,ph->d_pixel_flags,ph->image_size);
    cutilCheckMsg("module projection execution failed\n");
    cufftSafeCall(cufftExecC2C(ph->cufft_plan, ph->d_g1, ph->d_g1, CUFFT_INVERSE));
    /* normalize */
    CUDA_complex_scale<<<ph->number_of_blocks, ph->threads_per_block,0,ph->calc_stream>>>(ph->d_g1,ph->image_size, 1.0f / (ph->image_size));
    cutilCheckMsg("scale execution failed\n");
        CUDA_support_projection_hio<<<ph->number_of_blocks, ph->threads_per_block,0,ph->calc_stream>>>(ph->d_g1,ph->d_g0,ph->d_pixel_flags,ph->image_size,beta);
    cutilCheckMsg("support projection execution failed\n");
  }
  CUDA_apply_constraints<<<ph->number_of_blocks, ph->threads_per_block,0,ph->calc_stream>>>(ph->d_g1,ph->d_pixel_flags,ph->image_size,params->constraints);
  cudaStreamSynchronize(ph->transfer_stream);
  cudaStreamSynchronize(ph->calc_stream);
  phaser_cuda_handle_output_end(ph,output);
  return 0;
}

int phaser_iterate_raar_cuda(SpPhaser * ph,int iterations, SpPhasingOutput output){
  SpPhasingRAARParameters * params = (SpPhasingRAARParameters *)ph->algorithm->params;
  const real beta = params->beta;
  cudaStreamSynchronize(ph->calc_stream);
  phaser_cuda_handle_output_begin(ph,output);
  for(int i = 0;i<iterations;i++){
    cufftComplex * swap = ph->d_g0;
    ph->d_g0 = ph->d_g1;
    ph->d_g1 = swap;
    /* executes FFT processes */
    cufftSafeCall(cufftExecC2C(ph->cufft_plan, ph->d_g0, ph->d_g1, CUFFT_FORWARD));
    CUDA_module_projection<<<ph->number_of_blocks, ph->threads_per_block,0,ph->calc_stream>>>(ph->d_g1,ph->d_amplitudes,ph->d_pixel_flags,ph->image_size);
    cutilCheckMsg("module projection execution failed\n");
    cufftSafeCall(cufftExecC2C(ph->cufft_plan, ph->d_g1, ph->d_g1, CUFFT_INVERSE));
    /* normalize */
    CUDA_complex_scale<<<ph->number_of_blocks, ph->threads_per_block,0,ph->calc_stream>>>(ph->d_g1,ph->image_size, 1.0f / (ph->image_size));
    cutilCheckMsg("scale execution failed\n");
        CUDA_support_projection_raar<<<ph->number_of_blocks, ph->threads_per_block,0,ph->calc_stream>>>(ph->d_g1,ph->d_g0,ph->d_pixel_flags,ph->image_size,beta);
    cutilCheckMsg("support projection execution failed\n");
  }
  CUDA_apply_constraints<<<ph->number_of_blocks, ph->threads_per_block,0,ph->calc_stream>>>(ph->d_g1,ph->d_pixel_flags,ph->image_size,params->constraints);
  cudaStreamSynchronize(ph->transfer_stream);
  phaser_cuda_handle_output_end(ph,output);
  return 0;
}


static void phaser_cuda_handle_output_begin(SpPhaser * ph,SpPhasingOutput output){
  if(output & SpOutputModel){
    cutilSafeCall(cudaMemcpy(ph->d_g1_transfer,ph->d_g1,sizeof(cufftComplex)*ph->image_size,cudaMemcpyDeviceToDevice));
  }
  if(output & SpOutputModelChange){
    cutilSafeCall(cudaMemcpy(ph->d_g0_transfer,ph->d_g0,sizeof(cufftComplex)*ph->image_size,cudaMemcpyDeviceToDevice));
  }
  if(output & SpOutputModel){
    cutilSafeCall(cudaMemcpyAsync(ph->g1->image->data,ph->d_g1_transfer,sizeof(cufftComplex)*ph->image_size,cudaMemcpyDeviceToHost,ph->transfer_stream));
  }
  if(output & SpOutputModelChange){
    cutilSafeCall(cudaMemcpyAsync(ph->g0->image->data,ph->d_g0_transfer,sizeof(cufftComplex)*ph->image_size,cudaMemcpyDeviceToHost,ph->transfer_stream));
  }
}

static void phaser_cuda_handle_output_end(SpPhaser * ph,SpPhasingOutput output){
  if(output & SpOutputModel){
    sp_image_memcpy(ph->model,ph->g1);
  }
  if(output & SpOutputModelChange){
    sp_image_memcpy(ph->model_change,ph->g1);
    sp_image_sub(ph->model_change,ph->g0);
  }
}


