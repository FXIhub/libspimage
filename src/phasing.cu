#include <spimage.h>

__global__ void CUDA_module_projection(cufftComplex* g, const float* amp,const int * pixel_flags, const  int size);
__global__ void CUDA_support_projection_hio(cufftComplex* g1, const cufftComplex* g0,const int * pixel_flags,const  int size,const float beta);
__global__ void CUDA_support_projection_er(cufftComplex* g1,const int * pixel_flags,const  int size);
__global__ void CUDA_complex_scale(cufftComplex * a, int size ,float scale);
__global__ void CUDA_support_projection_raar(cufftComplex* g1, const cufftComplex* g0,const int * pixel_flags,const  int size,const float beta);
__global__ void CUDA_apply_constraints(cufftComplex* g, const int * pixel_flags,const  int size,const SpPhasingConstraints constraints);
__global__ void CUDA_apply_fourier_constraints(cufftComplex* g, const  int size,const SpPhasingConstraints constraints);
__global__ void CUDA_phased_amplitudes_projection(cufftComplex* g, const cufftComplex* phased_amp,const int * pixel_flags, const  int size);
__global__ void CUDA_diff_map_f1(cufftComplex* f1, const cufftComplex* g0,const int * pixel_flags,const float gamma1,const  int size);
__global__ void CUDA_diff_map(cufftComplex* Pi2f1,cufftComplex* Pi2rho, const cufftComplex* g0,cufftComplex* g1,const int * pixel_flags,const float gamma2,const float beta,const  int size);
__global__ void CUDA_ramp_final(cufftComplex *g, const float ax, const float ay, const int size, const int nx, const int *pixel_flags);

int sp_proj_module_cuda(Image * a, Image * amp){
  cufftComplex * d_a;
  int * d_pixel_flags;
  float * d_amp;
  cudaMalloc((void **)&d_a,sizeof(cufftComplex)*sp_image_size(a));
  cudaMalloc((void **)&d_pixel_flags,sizeof(int)*sp_image_size(a));
  cudaMalloc((void **)&d_amp,sizeof(float)*sp_image_size(a));
  cudaMemcpy(d_a,a->image->data,sizeof(cufftComplex)*sp_image_size(a),cudaMemcpyHostToDevice);
  sp_i3matrix * pixel_flags = sp_i3matrix_alloc(sp_image_x(a),sp_image_y(a),sp_image_z(a));
  sp_3matrix * h_amp = sp_3matrix_alloc(sp_image_x(a),sp_image_y(a),sp_image_z(a));
  for(int i =0 ;i<sp_image_size(a);i++){
    h_amp->data[i] = sp_real(amp->image->data[i]);
    pixel_flags->data[i] = 0;
    if(amp->mask->data[i]){
      pixel_flags->data[i] |= SpPixelMeasuredAmplitude;
    }
  }
  cudaMemcpy(d_amp,h_amp->data,sizeof(float)*sp_image_size(a),cudaMemcpyHostToDevice);
  cudaMemcpy(d_pixel_flags,pixel_flags->data,sizeof(int)*sp_image_size(a),cudaMemcpyHostToDevice);
  int threads_per_block = 64;
  int number_of_blocks = (sp_image_size(a)+threads_per_block-1)/threads_per_block;
  CUDA_module_projection<<<number_of_blocks, threads_per_block>>>(d_a,d_amp,d_pixel_flags,sp_image_size(a));
  cudaMemcpy(a->image->data,d_a,sizeof(cufftComplex)*sp_image_size(a),cudaMemcpyDeviceToHost);
  cudaFree(d_amp);
  cudaFree(d_pixel_flags);
  cudaFree(d_a);
  sp_cuda_check_errors();
  return 0;
}

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
      CUDA_module_projection<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_g1,ph->d_amplitudes,ph->d_pixel_flags,ph->image_size);
    }else if(ph->phasing_objective == SpRecoverAmplitudes){
      CUDA_phased_amplitudes_projection<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_g1,ph->d_phased_amplitudes,ph->d_pixel_flags,ph->image_size);
    }else{
      abort();
    }
    sp_cuda_check_errors();
    sp_cuda_check_errors();
    cufftExecC2C(ph->cufft_plan, ph->d_g1, ph->d_g1, CUFFT_INVERSE);
    /* normalize */
    CUDA_complex_scale<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_g1,ph->image_size, 1.0f / (ph->image_size));
    sp_cuda_check_errors();
    CUDA_support_projection_er<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_g1,ph->d_pixel_flags,ph->image_size);
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
      CUDA_module_projection<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_g1,ph->d_amplitudes,ph->d_pixel_flags,ph->image_size);
    }else if(ph->phasing_objective == SpRecoverAmplitudes){
      CUDA_phased_amplitudes_projection<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_g1,ph->d_phased_amplitudes,ph->d_pixel_flags,ph->image_size);
    }else{
      abort();
    }
    sp_cuda_check_errors();

    /* The fourier constraints cannot be applied in this location! See comment above */
    cufftExecC2C(ph->cufft_plan, ph->d_g1, ph->d_g1, CUFFT_INVERSE);
    /* normalize */
    CUDA_complex_scale<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_g1,ph->image_size, 1.0f / (ph->image_size));
    sp_cuda_check_errors();
    CUDA_support_projection_hio<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_g1,ph->d_g0,ph->d_pixel_flags,ph->image_size,beta);
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
    CUDA_module_projection<<<ph->number_of_blocks, ph->threads_per_block>>>(f1,ph->d_amplitudes,ph->d_pixel_flags,ph->image_size);
    cufftExecC2C(ph->cufft_plan, f1, f1, CUFFT_INVERSE);
    CUDA_module_projection<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_g1,ph->d_amplitudes,ph->d_pixel_flags,ph->image_size);
    cufftExecC2C(ph->cufft_plan, ph->d_g1, ph->d_g1, CUFFT_INVERSE);
    sp_cuda_check_errors();
    /* normalize */
    CUDA_complex_scale<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_g1,ph->image_size, 1.0f / (ph->image_size));
    CUDA_complex_scale<<<ph->number_of_blocks, ph->threads_per_block>>>(f1,ph->image_size, 1.0f / (ph->image_size));
    sp_cuda_check_errors();
    CUDA_diff_map<<<ph->number_of_blocks, ph->threads_per_block>>>(f1,ph->d_g1,ph->d_g0,ph->d_g1,ph->d_pixel_flags,gamma2,beta,ph->image_size);
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
      CUDA_module_projection<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_g1,ph->d_amplitudes,ph->d_pixel_flags,ph->image_size);
    }else if(ph->phasing_objective == SpRecoverAmplitudes){
      CUDA_phased_amplitudes_projection<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_g1,ph->d_phased_amplitudes,ph->d_pixel_flags,ph->image_size);
    }else{
      abort();
    }
    sp_cuda_check_errors();
    cufftExecC2C(ph->cufft_plan, ph->d_g1, ph->d_g1, CUFFT_INVERSE);
    /* normalize */
    CUDA_complex_scale<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_g1,ph->image_size, 1.0f / (ph->image_size));
    sp_cuda_check_errors();
    CUDA_support_projection_raar<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_g1,ph->d_g0,ph->d_pixel_flags,ph->image_size,beta);
    sp_cuda_check_errors();
    if(params->constraints != SpNoConstraints){
      CUDA_apply_constraints<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_g1,ph->d_pixel_flags,ph->image_size,params->constraints);
      if (params->constraints & SpRampObject) {
	apply_ramp_constraint_cuda(ph);
      }
    }
    sp_cuda_check_errors();
    ph->iteration++;
  }
  return 0;

}

int apply_ramp_constraint_cuda(SpPhaser *ph)
{
  float x_tmp, y_tmp;
  float x2,y2,xy,kx,ky;
  Image *tmp = sp_image_alloc(ph->nx,ph->ny,1);
  //int foo = tmp->image->x + tmp->image->y;
  //printf("%d\n",sp_image_y(tmp));
  //printf("%d x %d\n",sp_image_x(tmp),sp_image_y(tmp));
  printf("%d x %d\n",ph->nx,ph->ny);
  cudaMemcpy(tmp->image->data,ph->d_g1,sizeof(cufftComplex)*ph->image_size,cudaMemcpyDeviceToHost);
  cudaMemcpy(ph->pixel_flags->data,ph->d_pixel_flags,sizeof(int)*ph->image_size,cudaMemcpyDeviceToHost);
  int x,y;
  for (x = 0; x < ph->nx/2; x++) {
    x_tmp = (float) x;
    for (y = 0; y < ph->ny/2; y++) {
      if (sp_i3matrix_get(ph->pixel_flags,x,y,0) & SpPixelInsideSupport) {
	y_tmp = (float) y;
	x2 += x_tmp*x_tmp;
	y2 += y_tmp*y_tmp;
	xy += x_tmp*y_tmp;
	kx += x_tmp*sp_carg(sp_image_get(tmp,x,y,0));
	ky += y_tmp*sp_carg(sp_image_get(tmp,x,y,0));
      }
    }
    for (y = ph->ny/2; y < ph->ny; y++) {
      if (sp_i3matrix_get(ph->pixel_flags,x,y,0) & SpPixelInsideSupport) {
	y_tmp = (float) (y-ph->ny);
	x2 += x_tmp*x_tmp;
	y2 += y_tmp*y_tmp;
	xy += x_tmp*y_tmp;
	kx += x_tmp*sp_carg(sp_image_get(tmp,x,y,0));
	ky += y_tmp*sp_carg(sp_image_get(tmp,x,y,0));
      }
    }
  }
  for (x = ph->nx/2; x < ph->nx; x++) {
    x_tmp = (float) (x-ph->nx);
    for (y = 0; y < ph->ny/2; y++) {
      if (sp_i3matrix_get(ph->pixel_flags,x,y,0) & SpPixelInsideSupport) {
	y_tmp = (float) y;
	x2 += x_tmp*x_tmp;
	y2 += y_tmp*y_tmp;
	xy += x_tmp*y_tmp;
	kx += x_tmp*sp_carg(sp_image_get(tmp,x,y,0));
	ky += y_tmp*sp_carg(sp_image_get(tmp,x,y,0));
      }
    }
    for (y = ph->ny/2; y < ph->ny; y++) {
      if (sp_i3matrix_get(ph->pixel_flags,x,y,0) & SpPixelInsideSupport) {
	y_tmp = (float) (y-ph->ny);
	x2 += x_tmp*x_tmp;
	y2 += y_tmp*y_tmp;
	xy += x_tmp*y_tmp;
	kx += x_tmp*sp_carg(sp_image_get(tmp,x,y,0));
	ky += y_tmp*sp_carg(sp_image_get(tmp,x,y,0));
      }
    }
  }
  sp_image_free(tmp);
  /*
  for (int x = 0; x < sp_image_x(new_model); x++) {
    for (int y = 0; y < sp_image_y(new_model); y++) {
      if (sp_i3matrix_get(ph->pixel_flags,x,y,0) & SpPixelInsideSupport) {
	if (x < ph->nx/2) {
	  x_tmp = (real) x;
	} else {
	  x_tmp = (real)( x - ph->nx );
	}
	if (y < ph->ny/2) {
	  y_tmp = (real) y;
	} else {
	  y_tmp = (real)( y - ph->ny );
	}
	x2 += x_tmp*x_tmp;
	y2 += y_tmp*y_tmp;
	xy += x_tmp*y_tmp;
	kx += x_tmp*sp_carg(sp_image_get(new_model,x,y,0));
	ky += y_tmp*sp_carg(sp_image_get(new_model,x,y,0));
      }
    }
  }
  */
  float ax = (kx*y2-ky*xy) / (x2*y2 - xy*xy);
  float ay = (ky*x2-kx*xy) / (x2*y2 - xy*xy);
  printf("\nax = %g\nay = %g\nsize = %d\nnx = %d\n",ax,ay,ph->image_size,ph->nx);

  Image *tmp2 = sp_image_alloc(ph->nx,ph->ny,1);
  cudaMemcpy(tmp2->image->data,ph->d_g1,sizeof(cufftComplex)*ph->image_size,cudaMemcpyDeviceToHost);
  sp_image_write(tmp2,"debug_before_kernel.h5",0);

  CUDA_ramp_final<<<ph->number_of_blocks, ph->threads_per_block>>>(ph->d_g1, ax, ay, ph->image_size, ph->nx, ph->d_pixel_flags);

  cudaMemcpy(tmp2->image->data,ph->d_g1,sizeof(cufftComplex)*ph->image_size,cudaMemcpyDeviceToHost);
  sp_image_write(tmp2,"debug_after_kernel.h5",0);
  sp_image_free(tmp2);

  //sp_image_free(tmp);
  return 0;
}
