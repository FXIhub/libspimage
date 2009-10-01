#include <spimage.h>


static int phaser_iterate_hio(SpPhaser * ph, int iterations);
static int phaser_iterate_raar(SpPhaser * ph, int iterations);
static void phaser_apply_constraints(SpPhaser * ph,Image * new_model, SpPhasingConstraints constraints);

static void phaser_module_projection(Image * a, sp_3matrix * amp, sp_i3matrix * pixel_flags);

SpPhasingAlgorithm * sp_phasing_raar_alloc(real beta, SpPhasingConstraints constraints){
  SpPhasingAlgorithm * ret = sp_malloc(sizeof(SpPhasingAlgorithm));
  ret->type = SpRAAR;
  SpPhasingRAARParameters * params = sp_malloc(sizeof(SpPhasingRAARParameters));
  params->beta = beta;
  params->constraints = constraints;
  ret->params = params;
  return ret;
}

SpPhasingAlgorithm * sp_phasing_hio_alloc(real beta, SpPhasingConstraints constraints){
  SpPhasingAlgorithm * ret = sp_malloc(sizeof(SpPhasingAlgorithm));
  ret->type = SpHIO;
  SpPhasingHIOParameters * params = sp_malloc(sizeof(SpPhasingHIOParameters));
  params->beta = beta;
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
  return ret;
}

void sp_phaser_free(SpPhaser * ph){
  if(ph->model){
    sp_image_free(ph->model);
  }
  if(ph->model_change){
    sp_image_free(ph->model_change);
  }
#ifdef _USE_CUDA
  if(ph->engine == SpEngineCUDA){

    cudaFree(ph->d_amplitudes);
    cudaFree(ph->d_pixel_flags);
    cudaFree(ph->d_g0);
    cudaFree(ph->d_g1);
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
  ph->model_iteration = -1;
  if(ph->engine == SpEngineCPU){
    sp_image_memcpy(ph->g1,model);
  }else if(ph->engine == SpEngineCUDA){
#ifdef _USE_CUDA
    /* transfer the model from the graphics card to the main memory */
    cutilSafeCall(cudaMemcpy(ph->d_g1,model->image->data,sizeof(cufftComplex)*ph->image_size,cudaMemcpyHostToDevice));
#else
    return NULL;
#endif    
  }
}


void sp_phaser_set_support(SpPhaser * ph,const Image * support){
  ph->support_iteration = -1;
  for(int i = 0;i<ph->image_size;i++){
    if(sp_real(support->image->data[i])){
      ph->pixel_flags->data[i] |= SpPixelInsideSupport;
    }else{
      ph->pixel_flags->data[i] &= ~SpPixelInsideSupport;      
    }
  }
  if(ph->engine == SpEngineCUDA){
#ifdef _USE_CUDA
    /* transfer the model from the graphics card to the main memory */
    cutilSafeCall(cudaMemcpy(ph->d_pixel_flags,ph->pixel_flags->data,sizeof(int)*ph->image_size,cudaMemcpyHostToDevice));
#else
    return NULL;
#endif    
  }
}

const Image * sp_phaser_amplitudes(SpPhaser * ph){
  if(!ph->amplitudes_image){
    ph->amplitudes_image = sp_image_alloc(sp_3matrix_x(ph->amplitudes),
					  sp_3matrix_y(ph->amplitudes),
					  sp_3matrix_z(ph->amplitudes));
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


const Image * sp_phaser_fmodel(SpPhaser * ph){
  if(ph->fmodel_iteration != ph->iteration){
    if(!ph->fmodel){
      ph->fmodel = sp_image_alloc(sp_3matrix_x(ph->amplitudes),
				  sp_3matrix_y(ph->amplitudes),
				  sp_3matrix_z(ph->amplitudes));
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
  return ph->fmodel;
}

const Image * sp_phaser_old_model(SpPhaser * ph){
  if(ph->old_model_iteration != ph->iteration){
    if(!ph->old_model){
      ph->old_model = sp_image_alloc(sp_3matrix_x(ph->amplitudes),
				  sp_3matrix_y(ph->amplitudes),
				  sp_3matrix_z(ph->amplitudes));
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
 

int sp_phaser_init(SpPhaser * ph, SpPhasingAlgorithm * alg,SpSupportAlgorithm * sup_alg, Image * amplitudes,SpPhasingEngine engine){
  if(!ph){
    fprintf(stderr,"Phaser is NULL!\n");
    return -1;
  }
  if(!alg){
    fprintf(stderr,"Algorithm is NULL!\n");
    return -2;
  }
  if(!amplitudes){
    fprintf(stderr,"Amplitudes is NULL!\n");
    return -3;
  }
  ph->image_size = sp_image_size(amplitudes);
  ph->algorithm = alg;
  ph->sup_algorithm = sup_alg;
  ph->amplitudes = sp_3matrix_alloc(sp_image_x(amplitudes),sp_image_y(amplitudes),sp_image_z(amplitudes));
  ph->pixel_flags = sp_i3matrix_alloc(sp_image_x(amplitudes),sp_image_y(amplitudes),sp_image_z(amplitudes));
  for(int i =0 ;i<sp_image_size(amplitudes);i++){
    ph->amplitudes->data[i] = sp_real(sp_image_get_by_index(amplitudes,i));
    ph->pixel_flags->data[i] = 0;
    if(amplitudes->mask->data[i]){
      ph->pixel_flags->data[i] |= SpPixelMeasuredAmplitude;
    }
  }
  ph->iteration = 0;
  int maskedIn = 0;
  /* Do some sanity checks */
  for(int i = 0;i<sp_image_size(amplitudes);i++){
    if(amplitudes->mask->data[i]){
      maskedIn++;
    }
  }
  if(!maskedIn){
    fprintf(stderr,"Amplitudes mask is all zeros!\n");
    return -6;
  }
  /* default engine is CPU */
  ph->engine = SpEngineCPU;
#ifdef _USE_CUDA
  ph->threads_per_block = 64;
  ph->number_of_blocks = (ph->image_size+ph->threads_per_block-1)/ph->threads_per_block;
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
  if(!ph->amplitudes){
    return -2;
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
    Image * tmp = sp_image_alloc(sp_3matrix_x(ph->amplitudes),sp_3matrix_y(ph->amplitudes),sp_3matrix_z(ph->amplitudes));
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
    Image * tmp = sp_image_alloc(sp_3matrix_x(ph->amplitudes),sp_3matrix_y(ph->amplitudes),sp_3matrix_z(ph->amplitudes));
    for(int i = 0;i<sp_image_size(tmp);i++){
      sp_real(tmp->image->data[i]) = ph->amplitudes->data[i];
    }
    sp_image_rephase(tmp,SP_ZERO_PHASE);
    tmp->shifted = 1;
    ph->model = sp_image_ifft(tmp);
    sp_image_free(tmp);
    sp_image_scale(ph->model,1.0/sp_image_size(ph->model));
  }else if(flags & SpModelRandomValues){
    ph->model = sp_image_alloc(sp_3matrix_x(ph->amplitudes),sp_3matrix_y(ph->amplitudes),sp_3matrix_y(ph->amplitudes));
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
    cutilSafeCall(cudaMalloc((void**)&ph->d_amplitudes, sizeof(float)*ph->image_size));
    cutilSafeCall(cudaMalloc((void**)&ph->d_pixel_flags, sizeof(int)*ph->image_size));
    cutilSafeCall(cudaMalloc((void**)&ph->d_g0, sizeof(cufftComplex)*ph->image_size));
    cutilSafeCall(cudaMalloc((void**)&ph->d_g1, sizeof(cufftComplex)*ph->image_size));
    /*    cutilSafeCall(cudaMalloc((void**)&ph->d_g0_transfer, sizeof(cufftComplex)*ph->image_size));
	  cutilSafeCall(cudaMalloc((void**)&ph->d_g1_transfer, sizeof(cufftComplex)*ph->image_size));*/
    /* Used for image transfer */
    /*    ph->g0 = sp_image_duplicate(ph->model,SP_COPY_ALL);
    sp_image_fill(ph->g0,sp_cinit(0,0));
    ph->g1 = sp_image_duplicate(ph->model,SP_COPY_ALL);
    sp_free(ph->g0->image->data);
    sp_free(ph->g1->image->data);*/
    /* allocate page locked memory for image transfer */
    /*    cutilSafeCall(cudaMallocHost((void**)&ph->g0->image->data,ph->image_size*sizeof(cufftComplex)));
	  cutilSafeCall(cudaMallocHost((void**)&ph->g1->image->data,ph->image_size*sizeof(cufftComplex)));*/


    
    /* transfer to GPU memory */
    if(sizeof(real) == sizeof(float)){
      cutilSafeCall(cudaMemcpy(ph->d_amplitudes, ph->amplitudes->data, sizeof(float)*ph->image_size, cudaMemcpyHostToDevice));
    }else{
      /* we have to convert to real to float */
      float * dummy = sp_malloc(sizeof(float)*ph->image_size);
      for(int i = 0;i<ph->image_size;i++){
	dummy[i] = ph->amplitudes->data[i];
      }
      cutilSafeCall(cudaMemcpy(ph->d_amplitudes, dummy, sizeof(float)*ph->image_size, cudaMemcpyHostToDevice));
      sp_free(dummy);
    }
    cutilSafeCall(cudaMemcpy(ph->d_pixel_flags, ph->pixel_flags->data, sizeof(int)*ph->image_size, cudaMemcpyHostToDevice));
    cutilSafeCall(cudaMemcpy(ph->d_g1, ph->model->image->data, sizeof(cufftComplex)*ph->image_size, cudaMemcpyHostToDevice));
    cutilSafeCall(cudaMemset(ph->d_g0, 0, sizeof(cufftComplex)*ph->image_size));
    if(sp_image_z(ph->model) == 1){
      cufftSafeCall(cufftPlan2d(&ph->cufft_plan, sp_image_y(ph->model),sp_image_x(ph->model), CUFFT_C2C));
    }else{
      cufftSafeCall(cufftPlan3d(&ph->cufft_plan, sp_image_z(ph->model),sp_image_y(ph->model),sp_image_x(ph->model), CUFFT_C2C));
    }
  }
#endif 
  return 0;
}


int sp_phaser_init_support(SpPhaser * ph, const Image * support, int flags, real value){
  if(support){
    if(sp_3matrix_x(ph->amplitudes) != sp_image_x(support) || 
       sp_3matrix_y(ph->amplitudes) != sp_image_y(support) ||
       sp_3matrix_z(ph->amplitudes) != sp_image_z(support) ){
      fprintf(stderr,"Support and amplitudes dimensions don't match!\n");
      return -5;
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
    Image * tmp = sp_image_alloc(sp_3matrix_x(ph->amplitudes),sp_3matrix_y(ph->amplitudes),sp_3matrix_z(ph->amplitudes));
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
    cutilSafeCall(cudaMemcpy(ph->d_pixel_flags, ph->pixel_flags->data, sizeof(int)*ph->image_size, cudaMemcpyHostToDevice));
#else
    return -1;
#endif
  }
  return 0;
}

const Image * sp_phaser_support(SpPhaser * ph){
  if(!ph->support){
    ph->support = sp_image_alloc(sp_3matrix_x(ph->amplitudes),sp_3matrix_y(ph->amplitudes),sp_3matrix_z(ph->amplitudes));    
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
  int (*phaser_update_support_pointer)(SpPhaser *) = NULL; 
  if(!ph){
    return -1;
  }
  if(!ph->algorithm){
    return -2;
  }
  if(!ph->model){
    return -3;
  }
  if(!ph->amplitudes){
    return -4;
  }
  if(!ph->pixel_flags){
    return -5;
  }
  if(!ph->model_change){
    return -6;
  }
  if(ph->sup_algorithm){
    if(ph->engine == SpEngineCUDA){
#ifdef _USE_CUDA
      if(ph->sup_algorithm->type == SpSupportThreshold){
	phaser_update_support_pointer = sp_support_threshold_update_support_cuda;
      }
      if(ph->sup_algorithm->type == SpSupportArea){
	phaser_update_support_pointer = sp_support_area_update_support_cuda;
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
    }
  }
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
  }
  if(phaser_iterate_pointer == NULL){
    return -7;
  }
  int ret = 0;
  if(phaser_update_support_pointer){
    /* iterate up to the point of the support update */
    while(iterations){
      int to_support_update = ph->sup_algorithm->update_period-1-(ph->iteration)%ph->sup_algorithm->update_period;
      int to_iterate = sp_min(iterations,to_support_update);
      ret = phaser_iterate_pointer(ph,to_iterate);
      iterations -= to_iterate;
      to_support_update -= to_iterate;
      if(to_support_update == 0 && iterations > 0){
	phaser_update_support_pointer(ph);
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

static int phaser_iterate_hio(SpPhaser * ph,int iterations){
  SpPhasingHIOParameters * params = ph->algorithm->params;
  const real beta = params->beta;
  for(int i = 0;i<iterations;i++){
    Image * swap = ph->g0;
    ph->g0 = ph->g1;
    ph->g1 = swap;
    sp_image_fft_fast(ph->g0,ph->g1);
    phaser_module_projection(ph->g1,ph->amplitudes,ph->pixel_flags);
    sp_image_ifft_fast(ph->g1,ph->g1);
    sp_image_scale(ph->g1,1.0/sp_image_size(ph->g1));
    for(int i =0;i<sp_image_size(ph->g1);i++){
      if(ph->pixel_flags->data[i] & SpPixelInsideSupport){
	// Nothing to do here 
      }else{
	ph->g1->image->data[i] = sp_csub(ph->g0->image->data[i],sp_cscale(ph->g1->image->data[i],beta));
      }
    }
    phaser_apply_constraints(ph,ph->g1,params->constraints);
  }
  ph->iteration += iterations;
  return 0;
}


static int phaser_iterate_raar(SpPhaser * ph,int iterations){
  SpPhasingRAARParameters * params = ph->algorithm->params;
  const real beta = params->beta;
  for(int i = 0;i<iterations;i++){
    Image * swap = ph->g0;
    ph->g0 = ph->g1;
    ph->g1 = swap;
    sp_image_fft_fast(ph->g0,ph->g1);
    SpPhasingHIOParameters * params = ph->algorithm->params;
    phaser_module_projection(ph->g1,ph->amplitudes,ph->pixel_flags);
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
      }else{
	ph->g1->image->data[i] = sp_cadd(sp_cscale(ph->g1->image->data[i],1-2*beta),sp_cscale(ph->g0->image->data[i],beta));      
      }
    }
    phaser_apply_constraints(ph,ph->g1,params->constraints);
  }
  ph->iteration += iterations;
  return 0;
}
