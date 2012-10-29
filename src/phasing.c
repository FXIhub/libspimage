#include <spimage.h>


static int phaser_iterate_er(SpPhaser * ph,int iterations);
static int phaser_iterate_hio(SpPhaser * ph, int iterations);
static int phaser_iterate_raar(SpPhaser * ph, int iterations);
static int phaser_iterate_diff_map(SpPhaser * ph, int iterations);
static Image * phaser_iterate_diff_map_f1(Image * real_in,sp_i3matrix * pixel_flags,real gamma1);
static void phaser_apply_constraints(SpPhaser * ph,Image * new_model, SpPhasingConstraints constraints);

static void phaser_module_projection(Image * a, sp_3matrix * amp, sp_i3matrix * pixel_flags);
static void phaser_phased_amplitudes_projection(Image * a, sp_c3matrix * phased_amp, sp_i3matrix * pixel_flags);


static void phaser_check_dimensions(SpPhaser * ph,const Image * a);

SpPhasingAlgorithm * sp_phasing_diff_map_alloc(sp_smap * beta, real gamma1, real gamma2, SpPhasingConstraints constraints){
  SpPhasingAlgorithm * ret = sp_malloc(sizeof(SpPhasingAlgorithm));
  ret->type = SpDiffMap;
  SpPhasingDiffMapParameters * params = sp_malloc(sizeof(SpPhasingDiffMapParameters));
  params->beta = beta;
  if(isinf(gamma1)){    
    /* BUG: At the moment the gammas don't change with beta*/
    real beta = sp_smap_interpolate(params->beta,0);

    /* Automatically calculate the gamma1 value*/
    /* Optimal value according to 
       Veit Elser 2003 "Random projections and the optimization of an algorithm for phase retrieval 
    */
    real sigma = 0.1;
    params->gamma1 = -(4+(2+beta)*sigma + beta*sigma*sigma)/(beta*(4-sigma+sigma*sigma));
  }else{
    params->gamma1 = gamma1;
  }
  if(isinf(gamma1)){
    /* BUG: At the moment the gammas don't change with beta*/
    real beta = sp_smap_interpolate(params->beta,0);
    /* Optimal value according to 
       Veit Elser 2003 "Random projections and the optimization of an algorithm for phase retrieval 
    */
    params->gamma2 = (3-beta)/(2*beta);
  }else{
    params->gamma2 = gamma2;
  }
  params->constraints = constraints;
  ret->params = params;
  return ret;
}

SpPhasingAlgorithm * sp_phasing_raar_alloc(sp_smap * beta, real sigma_noise, SpPhasingConstraints constraints){
  SpPhasingAlgorithm * ret = sp_malloc(sizeof(SpPhasingAlgorithm));
  ret->type = SpRAAR;
  SpPhasingRAARParameters * params = sp_malloc(sizeof(SpPhasingRAARParameters));
  params->beta = beta;
  params->sigma_noise = sigma_noise;
  params->constraints = constraints;
  ret->params = params;
  return ret;
}

SpPhasingAlgorithm * sp_phasing_hio_alloc(sp_smap * beta, real sigma_noise, SpPhasingConstraints constraints){
  SpPhasingAlgorithm * ret = sp_malloc(sizeof(SpPhasingAlgorithm));
  ret->type = SpHIO;
  SpPhasingHIOParameters * params = sp_malloc(sizeof(SpPhasingHIOParameters));
  params->beta = beta;
  params->sigma_noise = sigma_noise;
  params->constraints = constraints;
  ret->params = params;
  return ret;
}

SpPhasingAlgorithm * sp_phasing_er_alloc(SpPhasingConstraints constraints){
  SpPhasingAlgorithm * ret = sp_malloc(sizeof(SpPhasingAlgorithm));
  ret->type = SpER;
  SpPhasingERParameters * params = sp_malloc(sizeof(SpPhasingERParameters));
  params->constraints = constraints;
  ret->params = params;
  return ret;
}

SpPhaser * sp_phaser_alloc(){
  SpPhaser * ret = sp_malloc(sizeof(SpPhaser));
  memset(ret,0,sizeof(SpPhaser));
  ret->model_iteration = -1;
  ret->model_change_iteration = -1;
  ret->support_iteration = -1;
  ret->phasing_objective = SpRecoverPhases;  
  return ret;
}

void sp_phaser_set_objective(SpPhaser * ph, SpPhasingObjective obj){
  ph->phasing_objective = obj;  
}

void sp_phaser_free(SpPhaser * ph){
  if(ph->model){
    sp_image_free(ph->model);
    ph->model = 0;
  }
  if(ph->model_change){
    sp_image_free(ph->model_change);
    ph->model_change = 0;
  }
  if(ph->amplitudes){
    sp_3matrix_free(ph->amplitudes);
    ph->amplitudes = 0;
  }
  if(ph->pixel_flags){
    sp_i3matrix_free(ph->pixel_flags);
    ph->pixel_flags = 0;
  }
  if (ph->phased_amplitudes) {
    sp_c3matrix_free(ph->phased_amplitudes);
    ph->phased_amplitudes = 0;
  }
  if (ph->amplitudes_image) {
    sp_image_free(ph->amplitudes_image);
    ph->amplitudes_image = 0;
  }
  if (ph->old_model) {
    sp_image_free(ph->old_model);
    ph->old_model = 0;
  }
  if (ph->support) {
    sp_image_free(ph->support);
    ph->support = 0;
  }
  if (ph->fmodel) {
    sp_image_free(ph->fmodel);
    ph->fmodel = 0;
  }
#ifdef _USE_CUDA
  if(ph->engine == SpEngineCUDA){

    cudaFree(ph->d_amplitudes);
    ph->d_amplitudes = 0;
    cudaFree(ph->d_pixel_flags);
    ph->d_pixel_flags = 0;
    cudaFree(ph->d_g0);
    ph->d_g0 = 0;
    cudaFree(ph->d_g1);
    ph->d_g1 = 0;
    if(ph->d_phased_amplitudes){
      cudaFree(ph->d_phased_amplitudes);
      ph->d_phased_amplitudes = 0;
    }
    if(ph->cufft_plan){
      cufftDestroy(ph->cufft_plan);
      ph->cufft_plan = 0;
    }
  }else{
    sp_image_free(ph->g0);
    sp_image_free(ph->g1);
  }
#else
  sp_image_free(ph->g0);
  sp_image_free(ph->g1);
#endif
  free(ph);
}


void sp_phaser_set_model(SpPhaser * ph,const Image * model){
  phaser_check_dimensions(ph,model);
  ph->model_iteration = -1;
  if(ph->engine == SpEngineCPU){
    sp_image_memcpy(ph->g1,model);
  }else if(ph->engine == SpEngineCUDA){
#ifdef _USE_CUDA
    /* transfer the model from the graphics card to the main memory */
    cutilSafeCall(cudaMemcpy(ph->d_g1,model->image->data,sizeof(cufftComplex)*ph->image_size,cudaMemcpyHostToDevice));
#else
    abort();
#endif    
  }
}


void sp_phaser_set_support(SpPhaser * ph,const Image * support){
  phaser_check_dimensions(ph,support);
  ph->support_iteration = -1;
  if(!ph->pixel_flags){
    ph->pixel_flags = sp_i3matrix_alloc(ph->nx,ph->ny,ph->nz);
    for(int i= 0;i<ph->image_size;i++){
      ph->pixel_flags->data[i] = 0;
    }
  }
  for(int i = 0;i<ph->image_size;i++){
    if(sp_real(support->image->data[i])){
      ph->pixel_flags->data[i] |= SpPixelInsideSupport;
    }else{
      ph->pixel_flags->data[i] &= ~SpPixelInsideSupport;      
    }
  }
  if(ph->engine == SpEngineCUDA){
#ifdef _USE_CUDA
    if(!ph->d_pixel_flags){
      cudaMalloc((void **)&ph->d_pixel_flags,sizeof(int)*ph->image_size);
    }
    /* transfer the model from the graphics card to the main memory */
    cudaMemcpy(ph->d_pixel_flags,ph->pixel_flags->data,sizeof(int)*ph->image_size,cudaMemcpyHostToDevice);
#else
    abort();
#endif    
  }
}


void sp_phaser_set_phased_amplitudes(SpPhaser * ph,const Image * phased_amplitudes){
  phaser_check_dimensions(ph,phased_amplitudes);
  if(!ph->phased_amplitudes){
    ph->phased_amplitudes = sp_c3matrix_alloc(ph->nx,ph->ny,ph->nz);
  }
  if(!ph->pixel_flags){
    ph->pixel_flags = sp_i3matrix_alloc(ph->nx,ph->ny,ph->nz);
    for(int i= 0;i<ph->image_size;i++){
      ph->pixel_flags->data[i] = 0;
    }
  }
  for(int i = 0;i<ph->image_size;i++){
    ph->phased_amplitudes->data[i] = phased_amplitudes->image->data[i];
    if(phased_amplitudes->mask->data[i]){
      ph->pixel_flags->data[i] |= SpPixelMeasuredAmplitude;
    }else{
      ph->pixel_flags->data[i] &= ~SpPixelMeasuredAmplitude;
    }
  }
  if(ph->engine == SpEngineCUDA){
#ifdef _USE_CUDA
    if(!ph->d_phased_amplitudes){
      cudaMalloc((void **)&ph->d_phased_amplitudes,sizeof(cufftComplex)*ph->image_size);
    }
    if(!ph->d_pixel_flags){
      cudaMalloc((void **)&ph->d_pixel_flags,sizeof(int)*ph->image_size);
    }
    /* transfer the model from the graphics card to the main memory */
    cudaMemcpy(ph->d_pixel_flags,ph->pixel_flags->data,sizeof(int)*ph->image_size,cudaMemcpyHostToDevice);
    cudaMemcpy(ph->d_phased_amplitudes,ph->phased_amplitudes->data,sizeof(cufftComplex)*ph->image_size,cudaMemcpyHostToDevice);
#else
    abort();
#endif    
  }
}

static void phaser_check_dimensions(SpPhaser * ph, const Image * a){
  if(!ph->nx){
    ph->nx = sp_image_x(a);
    ph->ny = sp_image_y(a);
    ph->nz = sp_image_z(a);
    ph->image_size = ph->nx*ph->ny*ph->nz;
#ifdef _USE_CUDA

    /* We will have to be very careful in the future about block size and number of blocks. 
       We should check what's available in the hardware and make our choice accordingly.
       At the moment i'm assuming we have something like an NVIDIA 280 GTX with 512 threads per 
       block maximum and 65535 x 65535 x 1 of maximum block size.
       
     */
    ph->threads_per_block = 64;
    while(sp_image_size(a) > 65535*ph->threads_per_block){
      ph->threads_per_block *= 2;
    }
    if(ph->threads_per_block > 512){
      sp_error_fatal("Image to large to handle with current CUDA code!");
    }
    ph->number_of_blocks = (ph->image_size+ph->threads_per_block-1)/ph->threads_per_block;
#endif
  }
  if(ph->nx != sp_image_x(a)){
    abort();
  }
  if(ph->ny != sp_image_y(a)){
    abort();
  }
  if(ph->nz != sp_image_z(a)){
    abort();
  }
}

void sp_phaser_set_amplitudes(SpPhaser * ph,const Image * amplitudes){
  phaser_check_dimensions(ph,amplitudes);
  if(!ph->amplitudes){
    ph->amplitudes = sp_3matrix_alloc(ph->nx,ph->ny,ph->nz);
  }
  if(!ph->pixel_flags){
    ph->pixel_flags = sp_i3matrix_alloc(ph->nx,ph->ny,ph->nz);
    for(int i= 0;i<ph->image_size;i++){
      ph->pixel_flags->data[i] = 0;
    }
  }
  int masked = 0;
  for(int i = 0;i<ph->image_size;i++){
    ph->amplitudes->data[i] = sp_real(amplitudes->image->data[i]);
    if(amplitudes->mask->data[i]){
      ph->pixel_flags->data[i] |= SpPixelMeasuredAmplitude;
      masked++;
    }else{
      ph->pixel_flags->data[i] &= ~SpPixelMeasuredAmplitude;
    }
  }
  if(masked == 0){
    fprintf(stderr,"Amplitudes mask is all zeros!\n");
  }
  if(ph->engine == SpEngineCUDA){
#ifdef _USE_CUDA
    if(!ph->d_amplitudes){
      cudaMalloc((void **)&ph->d_amplitudes,sizeof(float)*ph->image_size);
    }
    if(!ph->d_pixel_flags){
      cudaMalloc((void **)&ph->d_pixel_flags,sizeof(int)*ph->image_size);
    }
    /* transfer the model from the graphics card to the main memory */
    cutilSafeCall(cudaMemcpy(ph->d_pixel_flags,ph->pixel_flags->data,sizeof(int)*ph->image_size,cudaMemcpyHostToDevice));
    cutilSafeCall(cudaMemcpy(ph->d_amplitudes,ph->amplitudes->data,sizeof(float)*ph->image_size,cudaMemcpyHostToDevice));
#else
    abort();
#endif    
  }
}

const Image * sp_phaser_amplitudes(SpPhaser * ph){
  if(!ph->amplitudes_image){
    ph->amplitudes_image = sp_image_alloc(ph->nx,
					  ph->ny,
					  ph->nz);
    for(int i = 0;i<sp_image_size(ph->amplitudes_image);i++){
      sp_real(ph->amplitudes_image->image->data[i]) = ph->amplitudes->data[i];
      sp_imag(ph->amplitudes_image->image->data[i]) = 0;
      if(ph->pixel_flags->data[i] & SpPixelMeasuredAmplitude){
	ph->amplitudes_image->mask->data[i] = 1;
      }else{
	ph->amplitudes_image->mask->data[i] = 0;
      }
    }
  }
  return ph->amplitudes_image;
}

const Image * sp_phaser_model(SpPhaser * ph){
  if(ph->model_iteration != ph->iteration){
    if(!ph->model){
      ph->model = sp_image_alloc(ph->nx,ph->ny,ph->nz);
    }
    ph->model_iteration = ph->iteration;
    if(ph->engine == SpEngineCPU){
      sp_image_memcpy(ph->model,ph->g1);
    }else if(ph->engine == SpEngineCUDA){
#ifdef _USE_CUDA
      /* transfer the model from the graphics card to the main memory */
      cutilSafeCall(cudaMemcpy(ph->model->image->data,ph->d_g1,sizeof(cufftComplex)*ph->image_size,cudaMemcpyDeviceToHost));
#else
      return NULL;
#endif    
    }
  }
  return ph->model;
}

static Image * sp_phaser_model_non_const(SpPhaser * ph){
  if(ph->model_iteration != ph->iteration){
    if(!ph->model){
      ph->model = sp_image_alloc(ph->nx,ph->ny,ph->nz);
    }
    ph->model_iteration = ph->iteration;
    if(ph->engine == SpEngineCPU){
      sp_image_memcpy(ph->model,ph->g1);
    }else if(ph->engine == SpEngineCUDA){
#ifdef _USE_CUDA
      /* transfer the model from the graphics card to the main memory */
      cutilSafeCall(cudaMemcpy(ph->model->image->data,ph->d_g1,sizeof(cufftComplex)*ph->image_size,cudaMemcpyDeviceToHost));
#else
      return NULL;
#endif    
    }
  }
  return ph->model;
}

const Image * sp_phaser_model_with_support(SpPhaser * ph){
  Image * model = sp_phaser_model_non_const(ph);
  const Image * support = sp_phaser_support(ph);
  sp_image_image_to_mask(support,model);
  return model;
}

const Image * sp_phaser_fmodel(SpPhaser * ph){
  if(ph->fmodel_iteration != ph->iteration){
    if(!ph->fmodel){
      ph->fmodel = sp_image_alloc(ph->nx,ph->ny,ph->nz);
    }

    ph->fmodel_iteration = ph->iteration;
    if(ph->engine == SpEngineCPU){
      sp_image_memcpy(ph->fmodel,ph->g1);
      sp_image_fft_fast(ph->fmodel,ph->fmodel);
    }else if(ph->engine == SpEngineCUDA){
#ifdef _USE_CUDA
      /* transfer the model from the graphics card to the main memory */
      cutilSafeCall(cudaMemcpy(ph->fmodel->image->data,ph->d_g1,sizeof(cufftComplex)*ph->image_size,cudaMemcpyDeviceToHost));
      /* not really efficient here */
      sp_image_fft_fast(ph->fmodel,ph->fmodel);
#else
      return NULL;
#endif    
    }
  }
  ph->fmodel->phased = 1;
  return ph->fmodel;
}

static Image * sp_phaser_fmodel_non_const(SpPhaser * ph){
  if(ph->fmodel_iteration != ph->iteration){
    if(!ph->fmodel){
      ph->fmodel = sp_image_alloc(ph->nx,ph->ny,ph->nz);
    }

    ph->fmodel_iteration = ph->iteration;
    if(ph->engine == SpEngineCPU){
      sp_image_memcpy(ph->fmodel,ph->g1);
      sp_image_fft_fast(ph->fmodel,ph->fmodel);
    }else if(ph->engine == SpEngineCUDA){
#ifdef _USE_CUDA
      /* transfer the model from the graphics card to the main memory */
      cutilSafeCall(cudaMemcpy(ph->fmodel->image->data,ph->d_g1,sizeof(cufftComplex)*ph->image_size,cudaMemcpyDeviceToHost));
      /* not really efficient here */
      sp_image_fft_fast(ph->fmodel,ph->fmodel);
#else
      return NULL;
#endif    
    }
  }
  ph->fmodel->phased = 1;
  return ph->fmodel;
}

const Image * sp_phaser_fmodel_with_mask(SpPhaser * ph){
  Image *fmodel = sp_phaser_fmodel_non_const(ph);
  const Image *amplitudes = sp_phaser_amplitudes(ph);
  sp_image_mask_to_mask(amplitudes,fmodel);
  /*
  for (int i = 0; i < sp_image_size(fmodel); i++) {
    fmodel->mask->data[i] = amplitudes->mask->data[i];
  }
  */
  return fmodel;
}

const Image * sp_phaser_old_model(SpPhaser * ph){
  if(ph->old_model_iteration != ph->iteration){
    if(!ph->old_model){
      ph->old_model = sp_image_alloc(ph->nx,
				  ph->ny,
				  ph->nz);
    }
    ph->old_model_iteration = ph->iteration;
    if(ph->engine == SpEngineCPU){
      sp_image_memcpy(ph->old_model,ph->g0);
    }else if(ph->engine == SpEngineCUDA){
#ifdef _USE_CUDA
      /* transfer the model from the graphics card to the main memory */
      cutilSafeCall(cudaMemcpy(ph->old_model->image->data,ph->d_g0,sizeof(cufftComplex)*ph->image_size,cudaMemcpyDeviceToHost));
#else
      return NULL;
#endif    
    }
  }
  return ph->old_model;
}


Image * sp_phaser_model_change(SpPhaser * ph){
  if(ph->model_iteration != ph->iteration){
    sp_phaser_model(ph);
  }
  if(ph->model_change_iteration != ph->iteration){
    ph->model_change_iteration = ph->iteration;
    if(ph->engine == SpEngineCPU){
      sp_image_memcpy(ph->model_change,ph->model);
      sp_image_sub(ph->model_change,ph->g0);
    }else if(ph->engine == SpEngineCUDA){
#ifdef _USE_CUDA
      /* transfer the model from the graphics card to the main memory */
      cutilSafeCall(cudaMemcpy(ph->model_change->image->data,ph->d_g0,sizeof(cufftComplex)*ph->image_size,cudaMemcpyDeviceToHost));
      sp_image_scale(ph->model_change,-1.0f);
      sp_image_add(ph->model_change,ph->model);
#else
      return NULL;
#endif    
    }
  }
  return ph->model_change;
}
 


int sp_phaser_init(SpPhaser * ph, SpPhasingAlgorithm * alg,SpSupportArray * sup_alg,SpPhasingEngine engine){
  if(!ph){
    fprintf(stderr,"Phaser is NULL!\n");
    return -1;
  }
  if(!alg){
    fprintf(stderr,"Algorithm is NULL!\n");
    return -2;
  }
  ph->algorithm = alg;
  ph->sup_algorithm = sup_alg;
  ph->iteration = 0;
  /* default engine is CPU */
  ph->engine = SpEngineCPU;
#ifdef _USE_CUDA
  if(engine == SpEngineAutomatic || engine == SpEngineCUDA){
    if(sp_cuda_get_device_type() == SpCUDAHardwareDevice){
      ph->engine = SpEngineCUDA;
    }
  }
#endif
  if(engine != SpEngineAutomatic && engine != ph->engine){
    fprintf(stderr,"Requested engine unavailable. Defaulting to CPU engine.\n");
  }
  return 0;
} 

int sp_phaser_init_model(SpPhaser * ph, const Image * user_model, int flags){
  if(!ph){
    return -1;    
  }
  if(ph->model){
    sp_image_free(ph->model);
  }
  if(ph->model_change){
    sp_image_free(ph->model_change);
  }
  if(user_model){
    ph->model = sp_image_duplicate(user_model,SP_COPY_ALL);
  }else if(flags & SpModelRandomPhases){
    Image * tmp = sp_image_alloc(ph->nx,ph->ny,ph->nz);
    for(int i = 0;i<sp_image_size(tmp);i++){
      sp_real(tmp->image->data[i]) = ph->amplitudes->data[i];
    }
    /* randomize phases */
    sp_image_rephase(tmp,SP_RANDOM_PHASE);
    tmp->shifted = 1;
    ph->model = sp_image_ifft(tmp);
    sp_image_free(tmp);
    sp_image_scale(ph->model,1.0/sp_image_size(ph->model));
  }else if(flags & SpModelZeroPhases){
    Image * tmp = sp_image_alloc(ph->nx,ph->ny,ph->nz);
    for(int i = 0;i<sp_image_size(tmp);i++){
      sp_real(tmp->image->data[i]) = ph->amplitudes->data[i];
    }
    sp_image_rephase(tmp,SP_ZERO_PHASE);
    tmp->shifted = 1;
    ph->model = sp_image_ifft(tmp);
    sp_image_free(tmp);
    sp_image_scale(ph->model,1.0/sp_image_size(ph->model));
  }else if(flags & SpModelRandomValues){
    ph->model = sp_image_alloc(ph->nx,ph->ny,ph->nz);
    /* try to start with reasonable random values */
    double sum = 0;
    for(int i= 0;i<sp_image_size(ph->model);i++){
      sum += ph->amplitudes->data[i]*ph->amplitudes->data[i];
    }
    sum /= sp_image_size(ph->model);
    for(int i = 0;i<sp_image_size(ph->model);i++){
      ph->model->image->data[i] = sp_cinit(p_drand48()-0.5,p_drand48()-0.5);
    }
    double model_sum = sp_image_integrate2(ph->model);
    /* Make the square of the norm match */
    sp_image_scale(ph->model,sqrt(sum/model_sum));    
  }else{
    return -3;
  }
  if(flags & SpModelMaskedOutZeroed){
    Image * tmp = sp_image_fft(ph->model);
    for(int i = 0;i<sp_image_size(tmp);i++){
      if((ph->pixel_flags->data[i] & SpPixelMeasuredAmplitude) == 0){
	tmp->image->data[i] = sp_cinit(0,0);
      }
    }
    sp_image_free(ph->model);
    ph->model = sp_image_ifft(tmp);
    sp_image_scale(ph->model,1.0/sp_image_size(ph->model));
    sp_image_free(tmp);
  }
  ph->model->phased = 1;
  ph->model_change = sp_image_alloc(sp_image_x(ph->model),sp_image_y(ph->model),sp_image_z(ph->model));
  sp_image_fill(ph->model_change,sp_cinit(0,0));

  if(ph->engine == SpEngineCPU){
    ph->g0 = sp_image_duplicate(ph->model,SP_COPY_ALL);
    sp_image_fill(ph->g0,sp_cinit(0,0));
    ph->g1 = sp_image_duplicate(ph->model,SP_COPY_ALL);
  }
#ifdef _USE_CUDA
  if(ph->engine == SpEngineCUDA){
    /* allocate GPU memory */
    cudaMalloc((void**)&ph->d_g0, sizeof(cufftComplex)*ph->image_size);
    cudaMalloc((void**)&ph->d_g1, sizeof(cufftComplex)*ph->image_size);
    
    cutilSafeCall(cudaMemcpy(ph->d_g1, ph->model->image->data, sizeof(cufftComplex)*ph->image_size, cudaMemcpyHostToDevice));
    cutilSafeCall(cudaMemset(ph->d_g0, 0, sizeof(cufftComplex)*ph->image_size));
    if(sp_image_z(ph->model) == 1){
      cufftPlan2d(&ph->cufft_plan, sp_image_y(ph->model),sp_image_x(ph->model), CUFFT_C2C);
    }else{
      cufftPlan3d(&ph->cufft_plan, sp_image_z(ph->model),sp_image_y(ph->model),sp_image_x(ph->model), CUFFT_C2C);
    }
  }
#endif 
  return 0;
}


int sp_phaser_init_support(SpPhaser * ph, const Image * support, int flags, real value){
  if(support){
    phaser_check_dimensions(ph, support);
    if(!ph->pixel_flags){
      ph->pixel_flags = sp_i3matrix_alloc(ph->nx,ph->ny,ph->nz);
      for(int i= 0;i<ph->image_size;i++){
	ph->pixel_flags->data[i] = 0;
      }
    }
    /* use given support*/
    for(int i = 0;i<ph->image_size;i++){
      if(sp_real(sp_image_get_by_index(support,i))){
	ph->pixel_flags->data[i] |= SpPixelInsideSupport;
      }else{
	ph->pixel_flags->data[i] &= ~SpPixelInsideSupport;
      }
    }
  }else if(flags & SpSupportFromPatterson){
    Image * tmp = sp_image_alloc(ph->nx,ph->ny,ph->nz);
    for(int i = 0;i<ph->image_size;i++){
      sp_real(tmp->image->data[i]) = ph->amplitudes->data[i]*ph->amplitudes->data[i];
      sp_imag(tmp->image->data[i]) = 0;
    }
    tmp->phased = 1;
    tmp->shifted = 1;
    
    Image * patterson = sp_image_ifft(tmp);
    real abs_threshold = sp_image_max(patterson,0,0,0,0)*value;
    for(int i = 0;i<ph->image_size;i++){
      if(sp_cabs(patterson->image->data[i]) > abs_threshold){
	ph->pixel_flags->data[i] |= SpPixelInsideSupport;
      }else{
	ph->pixel_flags->data[i] &= ~SpPixelInsideSupport;
      }
    }
    sp_image_free(patterson);
  }else{
    return -1;
  }
  if(ph->engine == SpEngineCUDA){
#ifdef _USE_CUDA
    if(!ph->d_pixel_flags){
      cudaMalloc((void **)&ph->d_pixel_flags,sizeof(int)*ph->image_size);
    }
    cudaMemcpy(ph->d_pixel_flags, ph->pixel_flags->data, sizeof(int)*ph->image_size, cudaMemcpyHostToDevice);
#else
    return -1;
#endif
  }
  return 0;
}

const Image * sp_phaser_support(SpPhaser * ph){
  if(!ph->support){
    ph->support = sp_image_alloc(ph->nx,ph->ny,ph->nz);    
  }
  if(ph->iteration != ph->support_iteration){
    ph->support_iteration = ph->iteration;
    sp_image_fill(ph->support,sp_cinit(0,0));
    if(ph->engine == SpEngineCUDA){
#ifdef _USE_CUDA
      cutilSafeCall(cudaMemcpy(ph->pixel_flags->data,ph->d_pixel_flags,sizeof(int)*ph->image_size,cudaMemcpyDeviceToHost));
#else
      return NULL;
#endif
    }
    for(int i = 0;i<ph->image_size;i++){
      if(ph->pixel_flags->data[i] & SpPixelInsideSupport){
	sp_real(ph->support->image->data[i]) = 1;
      }
    }
  }
  return ph->support;
}

int sp_phaser_iterate(SpPhaser * ph, int iterations){
  int (*phaser_iterate_pointer)(SpPhaser *, int) = NULL; 
  // int (*phaser_update_support_pointer)(SpPhaser *) = NULL; 
  if(!ph){
    return -1;
  }
  if(!ph->algorithm){
    return -2;
  }
  if(!ph->model){
    return -3;
  }
  if((!ph->amplitudes && ph->phasing_objective == SpRecoverPhases)
     ||
     (!ph->phased_amplitudes && ph->phasing_objective == SpRecoverAmplitudes)){
    return -4;
  }
  if(!ph->pixel_flags){
    return -5;
  }
  if(!ph->model_change){
    return -6;
  }
  /*
  if(ph->sup_algorithm){
    if(ph->engine == SpEngineCUDA){
#ifdef _USE_CUDA
      if(ph->sup_algorithm->type == SpSupportThreshold){
	phaser_update_support_pointer = sp_support_threshold_update_support_cuda;
      }
      if(ph->sup_algorithm->type == SpSupportArea){
	phaser_update_support_pointer = sp_support_area_update_support_cuda;
      }
      if (ph->sup_algorithm->type == SpSupportTemplate){
	phaser_update_support_pointer = sp_support_template_update_support_cuda;
      }
      if (ph->sup_algorithm->type == SpSupportStatic){
	phaser_update_support_pointer = sp_support_static_update_support_cuda;
      }
#else
      return -8;
#endif
    }else{
      if(ph->sup_algorithm->type == SpSupportThreshold){
	phaser_update_support_pointer = sp_support_threshold_update_support;
      }
      if(ph->sup_algorithm->type == SpSupportArea){
	phaser_update_support_pointer = sp_support_area_update_support;
      }
      if(ph->sup_algorithm->type == SpSupportTemplate){
	phaser_update_support_pointer = sp_support_template_update_support;
      }
      if(ph->sup_algorithm->type == SpSupportStatic){
	phaser_update_support_pointer = sp_support_static_update_support;
      }
    }
  }
  */
  if(ph->algorithm->type == SpHIO){
    if(ph->engine == SpEngineCUDA){
#ifdef _USE_CUDA
      phaser_iterate_pointer = phaser_iterate_hio_cuda;
      /*      return phaser_iterate_hio_cuda(ph,iterations,output); */
#else
      return -8;
#endif
    }else {
      phaser_iterate_pointer = phaser_iterate_hio;
    }
  }else if(ph->algorithm->type == SpRAAR){
    if(ph->engine == SpEngineCUDA){
#ifdef _USE_CUDA
      phaser_iterate_pointer = phaser_iterate_raar_cuda;
#else
      return -8;
#endif
    }else{
      phaser_iterate_pointer = phaser_iterate_raar;
    }
  }else if(ph->algorithm->type == SpDiffMap){
    if(ph->engine == SpEngineCUDA){
#ifdef _USE_CUDA
      //      phaser_iterate_pointer = phaser_iterate_diff_map_cuda;
      phaser_iterate_pointer = phaser_iterate_diff_map_cuda;
#else
      return -8;
#endif
    }else{
      phaser_iterate_pointer = phaser_iterate_diff_map;
    }
  }else if(ph->algorithm->type == SpER){
    if(ph->engine == SpEngineCUDA){
#ifdef _USE_CUDA
      phaser_iterate_pointer = phaser_iterate_er_cuda;
#else
      return -8;
#endif
    }else{
      phaser_iterate_pointer = phaser_iterate_er;
    }
  }
  if(phaser_iterate_pointer == NULL){
    return -7;
  }
  int ret = 0;
  if(ph->sup_algorithm && ph->sup_algorithm->algorithms[0]->function){
    /* iterate up to the point of the support update */
    while(iterations){ 
      int to_support_update = ph->sup_algorithm->update_period-1-(ph->iteration)%ph->sup_algorithm->update_period;
      int to_iterate = sp_min(iterations,to_support_update);
      ret = phaser_iterate_pointer(ph,to_iterate);
      iterations -= to_iterate;
      to_support_update -= to_iterate;
      if(to_support_update == 0 && iterations > 0){
	//phaser_update_support_pointer(ph);
	//ph->sup_algorithm->function(ph);
	//((int(*)(SpPhaser *))ph->sup_algorithm->function)(ph);
	sp_support_array_update(ph->sup_algorithm,ph);
	
	ph->iteration++;
	iterations -= 1;
      }
    }
  }else{
    if(iterations){
      ret = phaser_iterate_pointer(ph,iterations);
    }
  }
  return ret;
}

static void phaser_apply_constraints(SpPhaser * ph,Image * new_model, SpPhasingConstraints constraints){
  /* Apply constraints */
  for(int i =0;i<sp_image_size(new_model);i++){
    if(ph->pixel_flags->data[i] & SpPixelInsideSupport){
    if(constraints & SpRealObject){
      sp_imag(new_model->image->data[i]) = 0;
    }else if(constraints & SpPositiveRealObject){
      if(sp_real(new_model->image->data[i]) < 0){
	if(constraints & SpPositivityFlipping){
	  sp_real(new_model->image->data[i]) = fabs(sp_real(new_model->image->data[i]));
	}else{
	  sp_real(new_model->image->data[i]) = 0;
	}
      }
      sp_imag(new_model->image->data[i]) = 0;	
    }else if(constraints & SpPositiveComplexObject){
      if(sp_real(new_model->image->data[i]) < 0){
	if(constraints & SpPositivityFlipping){
	  sp_real(new_model->image->data[i]) = fabs(sp_real(new_model->image->data[i]));
	}else{
	  sp_real(new_model->image->data[i]) = 0;
	}
      }
      if(sp_imag(new_model->image->data[i]) < 0){
	if(constraints & SpPositivityFlipping){
	  sp_imag(new_model->image->data[i]) = fabs(sp_imag(new_model->image->data[i]));
	}else{
	  sp_imag(new_model->image->data[i]) = 0;
	}
      }

    }
    }
  }
}

static void phaser_apply_fourier_constraints(SpPhaser * ph,Image * new_amplitudes, SpPhasingConstraints constraints){
  /* Apply constraints */
  for(int i =0;i<sp_image_size(new_amplitudes);i++){
    if(constraints & SpCentrosymmetricObject){
      sp_imag(new_amplitudes->image->data[i]) = 0;
    }
  }
}

static void phaser_module_projection(Image * a, sp_3matrix * amp, sp_i3matrix * pixel_flags){
  for(int i = 0;i<sp_image_size(a);i++){
    if(pixel_flags->data[i] & SpPixelMeasuredAmplitude){
      const float m = amp->data[i]/sp_cabs(a->image->data[i]);
      if(isfinite(m)){
	a->image->data[i] = sp_cscale(a->image->data[i],m);
      }else{
	sp_real(a->image->data[i]) = amp->data[i];
	sp_imag(a->image->data[i]) = 0;
      }
    }
  }  
}


static void phaser_phased_amplitudes_projection(Image * a, sp_c3matrix * phased_amp, sp_i3matrix * pixel_flags){
  for(int i = 0;i<sp_image_size(a);i++){
    if(pixel_flags->data[i] & SpPixelMeasuredAmplitude){
      a->image->data[i] = phased_amp->data[i];
    }
  }  
}

static int phaser_iterate_er(SpPhaser * ph,int iterations){
  SpPhasingERParameters * params = ph->algorithm->params;
  for(int i = 0;i<iterations;i++){
    Image * swap = ph->g0;
    ph->g0 = ph->g1;
    ph->g1 = swap;
    sp_image_fft_fast(ph->g0,ph->g1);
    phaser_apply_fourier_constraints(ph,ph->g1,params->constraints);
    if(ph->phasing_objective == SpRecoverPhases){
      phaser_module_projection(ph->g1,ph->amplitudes,ph->pixel_flags);
    }else if(ph->phasing_objective == SpRecoverAmplitudes){
      phaser_phased_amplitudes_projection(ph->g1,ph->phased_amplitudes,ph->pixel_flags);
    }else{
      abort();
    }
    sp_image_ifft_fast(ph->g1,ph->g1);
    sp_image_scale(ph->g1,1.0/sp_image_size(ph->g1));
    for(int i =0;i<sp_image_size(ph->g1);i++){
      if(ph->pixel_flags->data[i] & SpPixelInsideSupport){
	// Nothing to do here 
      }else{
	ph->g1->image->data[i] = sp_cinit(0,0);
      }
    }
    phaser_apply_constraints(ph,ph->g1,params->constraints);
    ph->iteration++;
  }
  return 0;
}

static int phaser_iterate_hio(SpPhaser * ph,int iterations){
  SpPhasingHIOParameters * params = ph->algorithm->params;
  float sigma_noise = params->sigma_noise;
  for(int i = 0;i<iterations;i++){
    real beta = sp_smap_interpolate(params->beta,ph->iteration);
    Image * swap = ph->g0;
    ph->g0 = ph->g1;
    ph->g1 = swap;
    sp_image_fft_fast(ph->g0,ph->g1);
    phaser_apply_fourier_constraints(ph,ph->g1,params->constraints);
    if(ph->phasing_objective == SpRecoverPhases){
      phaser_module_projection(ph->g1,ph->amplitudes,ph->pixel_flags);
    }else if(ph->phasing_objective == SpRecoverAmplitudes){
      phaser_phased_amplitudes_projection(ph->g1,ph->phased_amplitudes,ph->pixel_flags);
    }else{
      abort();
    }
    sp_image_ifft_fast(ph->g1,ph->g1);
    sp_image_scale(ph->g1,1.0/sp_image_size(ph->g1));
    for(int i =0;i<sp_image_size(ph->g1);i++){
      if(ph->pixel_flags->data[i] & SpPixelInsideSupport){
	// Nothing to do here 
      }

      else if(SpNoiseTolerance){
	if (sp_cabs(ph->g1->image->data[i]) > 3*sigma_noise){
	  ph->g1->image->data[i] = sp_csub(ph->g0->image->data[i],sp_cscale(ph->g1->image->data[i],beta));
	}
	else{
	  ph->g1->image->data[i] = sp_cinit(0, 0);
	}
      }
      else{
	ph->g1->image->data[i] = sp_csub(ph->g0->image->data[i],sp_cscale(ph->g1->image->data[i],beta));
      }
    }
    phaser_apply_constraints(ph,ph->g1,params->constraints);
    ph->iteration++;
  }
  return 0;
}


static int phaser_iterate_raar(SpPhaser * ph,int iterations){
  SpPhasingRAARParameters * params = ph->algorithm->params;
  float sigma_noise = params->sigma_noise;
  for(int i = 0;i<iterations;i++){
    real beta = sp_smap_interpolate(params->beta,ph->iteration);
    Image * swap = ph->g0;
    ph->g0 = ph->g1;
    ph->g1 = swap;
    sp_image_fft_fast(ph->g0,ph->g1);
    phaser_apply_fourier_constraints(ph,ph->g1,params->constraints);
    SpPhasingHIOParameters * params = ph->algorithm->params;
    if(ph->phasing_objective == SpRecoverPhases){
      phaser_module_projection(ph->g1,ph->amplitudes,ph->pixel_flags);
    }else if(ph->phasing_objective == SpRecoverAmplitudes){
      phaser_phased_amplitudes_projection(ph->g1,ph->phased_amplitudes,ph->pixel_flags);
    }else{
      abort();
    }
    sp_image_ifft_fast(ph->g1,ph->g1);
    sp_image_scale(ph->g1,1.0/sp_image_size(ph->g1));
    for(int i =0;i<sp_image_size(ph->g1);i++){
      /* A bit of documentation about the equation:
	 
	 Rs = 2*Ps-I; Rm = 2*Pm-I
	 
	 RAAR = 1/2 * beta * (RsRm + I) + (1 - beta) * Pm;    
	 RAAR = 2*beta*Ps*Pm+(1-2*beta)*Pm - beta * (Ps-I)
	 
	 Which reduces to:
	 
	 Inside the support: Pm
	 Outside the support: (1 - 2*beta)*Pm + beta*I
	 
      */    
      if(ph->pixel_flags->data[i] & SpPixelInsideSupport){
	// Nothing to do here 
      }

      else if(ph->pixel_flags->data[i] & SpNoiseTolerance){
	if (sp_cabs(ph->g1->image->data[i]) > 3*sigma_noise){
	  ph->g1->image->data[i] = sp_cadd(sp_cscale(ph->g1->image->data[i],1-2*beta),sp_cscale(ph->g0->image->data[i],beta));
	}
	else{
	  ph->g1->image->data[i] = sp_cinit(0, 0);
	}
      }

      else{
	ph->g1->image->data[i] = sp_cadd(sp_cscale(ph->g1->image->data[i],1-2*beta),sp_cscale(ph->g0->image->data[i],beta));      
      }
    }
    phaser_apply_constraints(ph,ph->g1,params->constraints);
    ph->iteration++;
  }
  return 0;
}


static int phaser_iterate_diff_map(SpPhaser * ph,int iterations){
  SpPhasingDiffMapParameters * params = ph->algorithm->params;
  const real gamma1 = params->gamma1;
  const real gamma2 = params->gamma2;
  for(int i = 0;i<iterations;i++){
    real beta = sp_smap_interpolate(params->beta,ph->iteration);
    Image * swap = ph->g0;
    ph->g0 = ph->g1;
    ph->g1 = swap;
    sp_image_fft_fast(ph->g0,ph->g1);
    phaser_apply_fourier_constraints(ph,ph->g1,params->constraints);
    Image * f1 = phaser_iterate_diff_map_f1(ph->g0,ph->pixel_flags,gamma1);
    sp_image_fft_fast(f1,f1);
    phaser_module_projection(f1,ph->amplitudes,ph->pixel_flags);
    sp_image_ifft_fast(f1,f1);
    Image * Pi2f1 = f1;
    phaser_module_projection(ph->g1,ph->amplitudes,ph->pixel_flags);
    sp_image_ifft_fast(ph->g1,ph->g1);
    Image * Pi2rho = ph->g1;
    
    int size = ph->image_size;
    for(int i = 0;i<size;i++){
      if(ph->pixel_flags->data[i] & SpPixelInsideSupport){
	sp_real(ph->g1->image->data[i]) = sp_real(ph->g0->image->data[i])+(beta)*((1+gamma2)*sp_real(Pi2rho->image->data[i])/size-gamma2*sp_real(ph->g0->image->data[i]));
	sp_imag(ph->g1->image->data[i]) = sp_imag(ph->g0->image->data[i])+(beta)*((1+gamma2)*sp_imag(Pi2rho->image->data[i])/size-gamma2*sp_imag(ph->g0->image->data[i]));
      }else{
	ph->g1->image->data[i] = ph->g0->image->data[i];
      }
      sp_real(ph->g1->image->data[i]) -= beta*sp_real(Pi2f1->image->data[i])/size;
      sp_imag(ph->g1->image->data[i]) -= beta*sp_imag(Pi2f1->image->data[i])/size;
    }
    sp_image_free(Pi2f1);
    phaser_apply_constraints(ph,ph->g1,params->constraints);
    ph->iteration++;
  }
  return 0;
}

Image * phaser_iterate_diff_map_f1(Image * real_in,sp_i3matrix * pixel_flags,real gamma1){
  Image * out = sp_image_duplicate(real_in,SP_COPY_DATA|SP_COPY_MASK);
  for(int i = 0;i<sp_image_size(out);i++){
    if(pixel_flags->data[i] & SpPixelInsideSupport){
      /* inside support */
      /* 1+gamma1-gamma1 is 1 so  do nothing */
    }else{
      /* outside support */
      sp_real(out->image->data[i]) = -gamma1*sp_real(out->image->data[i]);
      sp_imag(out->image->data[i]) = -gamma1*sp_imag(out->image->data[i]);

    }
  }
  return out;
}
