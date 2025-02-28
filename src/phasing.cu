#include <spimage.h>
#include <thrust/reduce.h>
#include <thrust/device_vector.h>

__global__ void CUDA_module_projection(cufftComplex* g, const float* amp, const float* amp_min, const float* amp_max, const int * pixel_flags,const  int size, const SpPhasingConstraints constraints);
__global__ void CUDA_support_projection_hio(cufftComplex* g1, const cufftComplex* g0, const cufftComplex* gp, const int * pixel_flags,const  int size, const float beta);
__global__ void CUDA_support_projection_er(cufftComplex* g1, cufftComplex *gp, const int * pixel_flags, const  int size);
__global__ void CUDA_complex_scale(cufftComplex * a, int size ,float scale);
__global__ void CUDA_complex_add(cufftComplex * a, int size ,cufftComplex add);
__global__ void CUDA_support_projection_raar(cufftComplex* g1, const cufftComplex* g0, const cufftComplex* gp, const int * pixel_flags,const  int size, const float beta);
__global__ void CUDA_apply_constraints(cufftComplex* g, const int * pixel_flags,const  int size,const SpPhasingConstraints constraints);
__global__ void CUDA_apply_fourier_constraints(cufftComplex* g, const  int size,const SpPhasingConstraints constraints);
__global__ void CUDA_phased_amplitudes_projection(cufftComplex* g, const cufftComplex* phased_amp,const int * pixel_flags, const  int size);
__global__ void CUDA_diff_map_f1(cufftComplex* f1, const cufftComplex* g0,const int * pixel_flags,const float gamma1,const  int size);
__global__ void CUDA_diff_map(cufftComplex* Pi2f1,cufftComplex* Pi2rho, const cufftComplex* g0,cufftComplex* g1,const int * pixel_flags,const float gamma2,const float beta,const  int size);
__global__ void CUDA_random_rephase(cufftComplex * a, float * uniform_random, int size);
__global__ void CUDA_real_to_complex(cufftComplex * out, float * in, int size);
__global__ void CUDA_complex_abs2(cufftComplex * a, int size);
__global__ void CUDA_ereal(cufftComplex * out, const cufftComplex * in, const int * pixel_flags, int size);
__global__ void CUDA_efourier(cufftComplex * out, const cufftComplex * fmodel, const float* amp, const int * pixel_flags, int size);
__global__ void CUDA_FcFo(cufftComplex * out, const cufftComplex * fmodel, const float* amp, const int * pixel_flags, int size);
__global__ void CUDA_pixel_flags_to_complex(cufftComplex * out, const int * pixel_flags, int size);

struct addCufftComplex{   
  __device__ cufftComplex operator()(const cufftComplex lhs, const cufftComplex rhs) { 
    cufftComplex temp = lhs;
#ifdef _STRICT_IEEE_754
    temp.x = __fadd_rn(temp.x,rhs.x);
    temp.y = __fadd_rn(temp.y,rhs.y);
#else
    temp.x += rhs.x;
    temp.y += rhs.y;
#endif
    return temp;
  } 
};

static void random_rephase_cuda(SpPhaser * ph, cufftComplex *  img);

int phaser_iterate_er_cuda(SpPhaser * ph,int iterations){
  SpPhasingERParameters * params = (SpPhasingERParameters *)ph->algorithm->params;
  for(int i = 0;i<iterations;i++){
    cufftComplex * swap = ph->d_g0;
    ph->d_g0 = ph->d_g1;
    ph->d_g1 = swap;
    /* executes FFT processes */
    cufftExecC2C(ph->cufft_plan, ph->d_g0, ph->d_g1, CUFFT_FORWARD);
    CUDA_apply_fourier_constraints<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_g1,ph->image_size,params->constraints);    
    if(ph->phasing_objective == SpRecoverPhases){
      CUDA_module_projection<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_g1,ph->d_amplitudes,ph->d_amplitudes_min,ph->d_amplitudes_max,ph->d_pixel_flags,ph->image_size,params->constraints);
    }else if(ph->phasing_objective == SpRecoverAmplitudes){
      CUDA_phased_amplitudes_projection<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_g1,ph->d_phased_amplitudes,ph->d_pixel_flags,ph->image_size);
    }else{
      abort();
    }
    sp_cuda_check_errors();
    sp_cuda_check_errors();
    cufftExecC2C(ph->cufft_plan, ph->d_g1, ph->d_gp, CUFFT_INVERSE);
    /* normalize */
    CUDA_complex_scale<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_gp,ph->image_size, 1.0f / (ph->image_size));
    sp_cuda_check_errors();
    CUDA_support_projection_er<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_g1, ph->d_gp, ph->d_pixel_flags,ph->image_size);
    sp_cuda_check_errors();
    if(params->constraints != SpNoConstraints){
      CUDA_apply_constraints<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_g1,ph->d_pixel_flags,ph->image_size,params->constraints);
    }
    sp_cuda_check_errors();
  }
  ph->iteration += iterations;
  return 0;
}

int phaser_iterate_hio_cuda(SpPhaser * ph,int iterations){
  SpPhasingHIOParameters * params = (SpPhasingHIOParameters *)ph->algorithm->params;
  for(int i = 0;i<iterations;i++){
    real beta = sp_smap_interpolate(params->beta,ph->iteration);
    cufftComplex * swap = ph->d_g0;
    ph->d_g0 = ph->d_g1;
    ph->d_g1 = swap;
    /* executes FFT processes */
    cufftExecC2C(ph->cufft_plan, ph->d_g0, ph->d_g1, CUFFT_FORWARD);
    /* The fourier constraints have to be applied before the amplitude projection otherwise the algorithm never converges,
     probably because there is a deficit of power after the constraints */
    CUDA_apply_fourier_constraints<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_g1,ph->image_size,params->constraints);
    sp_cuda_check_errors();
    if(ph->phasing_objective == SpRecoverPhases){
      CUDA_module_projection<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_g1,ph->d_amplitudes,ph->d_amplitudes_min,ph->d_amplitudes_max,ph->d_pixel_flags,ph->image_size,params->constraints);
    }else if(ph->phasing_objective == SpRecoverAmplitudes){
      CUDA_phased_amplitudes_projection<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_g1,ph->d_phased_amplitudes,ph->d_pixel_flags,ph->image_size);
    }else{
      abort();
    }
    sp_cuda_check_errors();

    /* The fourier constraints cannot be applied in this location! See comment above */
    cufftExecC2C(ph->cufft_plan, ph->d_g1, ph->d_gp, CUFFT_INVERSE);
    /* normalize */
    CUDA_complex_scale<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_gp,ph->image_size, 1.0f / (ph->image_size));
    sp_cuda_check_errors();
    CUDA_support_projection_hio<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_g1,ph->d_g0,ph->d_gp,ph->d_pixel_flags,ph->image_size,beta);
    sp_cuda_check_errors();
    if(params->constraints != SpNoConstraints){
      CUDA_apply_constraints<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_g1,ph->d_pixel_flags,ph->image_size,params->constraints);
    }
    sp_cuda_check_errors();
    ph->iteration++; 
  }
  return 0;
}

int phaser_iterate_diff_map_cuda(SpPhaser * ph,int iterations){
  SpPhasingDiffMapParameters * params = (SpPhasingDiffMapParameters *)ph->algorithm->params;
  const real gamma1 = params->gamma1;
  const real gamma2 = params->gamma2;
  cufftComplex * f1;
  cudaMalloc((void **)&f1,sizeof(cufftComplex)*ph->image_size);
  for(int i = 0;i<iterations;i++){
    real beta = sp_smap_interpolate(params->beta,ph->iteration);
    cufftComplex * swap = ph->d_g0;
    ph->d_g0 = ph->d_g1;
    ph->d_g1 = swap;
    /* executes FFT processes */
    cufftExecC2C(ph->cufft_plan, ph->d_g0, ph->d_g1, CUFFT_FORWARD);
    CUDA_apply_fourier_constraints<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_g1,ph->image_size,params->constraints);
    CUDA_diff_map_f1<<<ph->number_of_blocks, ph->threads_per_block>>>(f1,ph->d_g0,ph->d_pixel_flags,gamma1,ph->image_size);
    cufftExecC2C(ph->cufft_plan, f1, f1, CUFFT_FORWARD);
    CUDA_module_projection<<<ph->number_of_blocks, ph->threads_per_block>>>(f1,ph->d_amplitudes,ph->d_amplitudes_min,ph->d_amplitudes_max,ph->d_pixel_flags,ph->image_size,params->constraints); 
    cufftExecC2C(ph->cufft_plan, f1, f1, CUFFT_INVERSE);
    CUDA_module_projection<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_g1,ph->d_amplitudes,ph->d_amplitudes_min,ph->d_amplitudes_max,ph->d_pixel_flags,ph->image_size,params->constraints);
    cufftExecC2C(ph->cufft_plan, ph->d_g1, ph->d_gp, CUFFT_INVERSE);
    sp_cuda_check_errors();
    /* normalize */
    CUDA_complex_scale<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_gp,ph->image_size, 1.0f / (ph->image_size));
    CUDA_complex_scale<<<ph->number_of_blocks, ph->threads_per_block>>>(f1,ph->image_size, 1.0f / (ph->image_size));
    sp_cuda_check_errors();
    CUDA_diff_map<<<ph->number_of_blocks, ph->threads_per_block>>>(f1,ph->d_gp,ph->d_g0,ph->d_g1,ph->d_pixel_flags,gamma2,beta,ph->image_size);
    sp_cuda_check_errors();
    if(params->constraints != SpNoConstraints){
      CUDA_apply_constraints<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_g1,ph->d_pixel_flags,ph->image_size,params->constraints);
    }
    sp_cuda_check_errors();
    ph->iteration++; 
  }
  cudaFree(f1);
  return 0;
}

int phaser_iterate_raar_cuda(SpPhaser * ph,int iterations){
  SpPhasingRAARParameters * params = (SpPhasingRAARParameters *)ph->algorithm->params;
  for(int i = 0;i<iterations;i++){
    real beta = sp_smap_interpolate(params->beta,ph->iteration);
    cufftComplex * swap = ph->d_g0;
    ph->d_g0 = ph->d_g1;
    ph->d_g1 = swap;
    /* executes FFT processes */
    cufftExecC2C(ph->cufft_plan, ph->d_g0, ph->d_g1, CUFFT_FORWARD);
    sp_cuda_check_errors();
    CUDA_apply_fourier_constraints<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_g1,ph->image_size,params->constraints);
    sp_cuda_check_errors();
    if(ph->phasing_objective == SpRecoverPhases){
      CUDA_module_projection<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_g1,ph->d_amplitudes,ph->d_amplitudes_min,ph->d_amplitudes_max,ph->d_pixel_flags,ph->image_size,params->constraints);
    }else if(ph->phasing_objective == SpRecoverAmplitudes){
      CUDA_phased_amplitudes_projection<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_g1,ph->d_phased_amplitudes,ph->d_pixel_flags,ph->image_size);
    }else{
      abort();
    }
    sp_cuda_check_errors();
    cufftExecC2C(ph->cufft_plan, ph->d_g1, ph->d_gp, CUFFT_INVERSE);
    /* normalize */
    CUDA_complex_scale<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_gp,ph->image_size, 1.0f / (ph->image_size));
    sp_cuda_check_errors();
    CUDA_support_projection_raar<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_g1,ph->d_g0,ph->d_gp,ph->d_pixel_flags,ph->image_size,beta);
    sp_cuda_check_errors();
    if(params->constraints != SpNoConstraints){
      CUDA_apply_constraints<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_g1,ph->d_pixel_flags,ph->image_size,params->constraints);
    }
    sp_cuda_check_errors();
    ph->iteration++;
  }
  return 0;

}

int sp_phaser_init_model_cuda(SpPhaser * ph, const Image * user_model, int flags){
  if(!ph){
    return -1;
  }
  if(ph->model){
    sp_image_free(ph->model);
  }
  if(ph->model_change){
    sp_image_free(ph->model_change);
  }
  /* allocate GPU memory */
  cutilSafeCall(cudaMalloc((void**)&ph->d_g0, sizeof(cufftComplex)*ph->image_size));
  cutilSafeCall(cudaMalloc((void**)&ph->d_g1, sizeof(cufftComplex)*ph->image_size));
  cutilSafeCall(cudaMalloc((void**)&ph->d_gp, sizeof(cufftComplex)*ph->image_size));
  cutilSafeCall(cudaMalloc((void**)&ph->d_fmodel, sizeof(cufftComplex)*ph->image_size));

  cutilSafeCall(cudaMemset(ph->d_g0, 0, sizeof(cufftComplex)*ph->image_size));
  cutilSafeCall(cudaMemset(ph->d_gp, 0, sizeof(cufftComplex)*ph->image_size));
  if(ph->nz == 1){
    cufftPlan2d(&ph->cufft_plan, ph->ny, ph->nx, CUFFT_C2C);
  }else{
    cufftPlan3d(&ph->cufft_plan, ph->nz, ph->ny, ph->nx, CUFFT_C2C);
  }

  ph->model = sp_image_alloc(ph->nx,ph->ny,ph->nz);
  ph->model->phased = 1;
  if(user_model){
    cutilSafeCall(cudaMemcpy(ph->d_g1, user_model->image->data, sizeof(cufftComplex)*ph->image_size, cudaMemcpyHostToDevice));
  }else if(flags & SpModelRandomPhases){
    CUDA_real_to_complex<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_g1,ph->d_amplitudes, ph->image_size);
    /* randomize phases */
    random_rephase_cuda(ph, ph->d_g1);
    cufftSafeCall(cufftExecC2C(ph->cufft_plan, ph->d_g1, ph->d_g1, CUFFT_INVERSE));
    CUDA_complex_scale<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_g1,ph->image_size, 1.0/sp_image_size(ph->model));
  }else if(flags & SpModelZeroPhases){
    cutilSafeCall(cudaMemcpy(ph->d_g1,  ph->amplitudes->data, sizeof(cufftComplex)*ph->image_size, cudaMemcpyHostToDevice));
    /* randomize phases */
    cufftSafeCall(cufftExecC2C(ph->cufft_plan, ph->d_g1, ph->d_g1, CUFFT_INVERSE));
    CUDA_complex_scale<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_g1,ph->image_size, 1.0/sp_image_size(ph->model));
  }else if(flags & SpModelRandomValues){
    curandSafeCall(curandGenerateUniform(ph->gen, (float *)ph->d_g1, ph->image_size*2));
    cufftComplex add;
    add.x = -0.5;
    add.y = -0.5;
    CUDA_complex_add<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_g1,ph->image_size, add);
    /* Note the result does not follow Parseval's theorem as in the CPU code */
  }else{
    return -3;
  }
  if(flags & SpModelMaskedOutZeroed){
    sp_error_fatal("Not implemented in CUDA");
  }
  ph->model->phased = 1;
  ph->model_change = sp_image_alloc(sp_image_x(ph->model),sp_image_y(ph->model),sp_image_z(ph->model));
  cutilSafeCall(cudaMemcpy(ph->model->image->data, ph->d_g1, sizeof(cufftComplex)*ph->image_size, cudaMemcpyDeviceToHost));  
  return 0;
}

real sp_phaser_ereal_cuda(SpPhaser * ph){
  /* CUDA_ereal takes the model before projection,  d_gp, and calculate the error pixelwise and stores in d_g0 */
  CUDA_ereal<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_g0, ph->d_gp, ph->d_pixel_flags,ph->image_size);
  thrust::device_ptr<cufftComplex> beginc =  thrust::device_pointer_cast(ph->d_g0);
  thrust::device_ptr<cufftComplex> endc =  thrust::device_pointer_cast((cufftComplex *)(ph->d_g0+ph->image_size));
  cufftComplex sum = {0,0};
  sum = thrust::reduce(beginc,endc,sum,addCufftComplex());
  real ereal = sqrt(sum.x / sum.y);
  return ereal;
}

real sp_phaser_support_fraction_cuda(SpPhaser * ph){
  /* CUDA_ereal takes the model before projection,  d_gp, and calculate the error pixelwise and stores in d_g0 */
  CUDA_pixel_flags_to_complex<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_g0, ph->d_pixel_flags,ph->image_size);
  thrust::device_ptr<cufftComplex> beginc =  thrust::device_pointer_cast(ph->d_g0);
  thrust::device_ptr<cufftComplex> endc =  thrust::device_pointer_cast((cufftComplex *)(ph->d_g0+ph->image_size));
  cufftComplex sum = {0,0};
  sum = thrust::reduce(beginc,endc,sum,addCufftComplex());
  real sup_frac = sum.x / real(ph->image_size);
  return sup_frac;  
}

real sp_phaser_efourier_cuda(SpPhaser * ph, real * FcFo){
  /* Calculate fmodel */
  cufftSafeCall(cufftExecC2C(ph->cufft_plan, ph->d_g1, ph->d_fmodel, CUFFT_FORWARD));
  /* CUDA_ereal takes the model before projection,  d_gp, and calculate the error pixelwise and stores in d_g0 */
  CUDA_efourier<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_g0, ph->d_fmodel, ph->d_amplitudes, ph->d_pixel_flags,ph->image_size);
  thrust::device_ptr<cufftComplex> beginc =  thrust::device_pointer_cast(ph->d_g0);
  thrust::device_ptr<cufftComplex> endc =  thrust::device_pointer_cast((cufftComplex *)(ph->d_g0+ph->image_size));
  cufftComplex sum = {0,0};
  sum = thrust::reduce(beginc,endc,sum,addCufftComplex());
  real efourier = sqrt(sum.x / sum.y);
  if(FcFo != NULL){
    /* store the FcFo ratio in this pointer */
    CUDA_FcFo<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_g0, ph->d_fmodel, ph->d_amplitudes, ph->d_pixel_flags,ph->image_size);
    sum.x = 0;
    sum.y = 0;
    sum = thrust::reduce(beginc,endc,sum,addCufftComplex());
    *FcFo = sum.x/sum.y;
  }
  return efourier;  
}

static void random_rephase_cuda(SpPhaser * ph, cufftComplex *  img){
  float * d_uni;
  /* Allocate n floats on device */
  cutilSafeCall(cudaMalloc((void **)&d_uni, ph->image_size*sizeof(float)));
  curandSafeCall(curandGenerateUniform(ph->gen, d_uni, ph->image_size));
  CUDA_random_rephase<<<ph->number_of_blocks, ph->threads_per_block>>>(img,d_uni,ph->image_size);
}
