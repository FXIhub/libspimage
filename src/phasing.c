#include <spimage.h>


static int phaser_iterate_hio(SpPhaser * ph);
static int phaser_iterate_raar(SpPhaser * ph);
static void phase_apply_constraints(SpPhaser * ph,Image * new_model, SpPhasingConstraints constraints);

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
  return ret;
}

void sp_phaser_free(SpPhaser * ph){
  if(ph->model){
    sp_image_free(ph->model);
  }
  if(ph->model_change){
    sp_image_free(ph->model_change);
  }
  free(ph);
}

Image * sp_phaser_model(const SpPhaser * ph){
  return ph->model;
}

Image * sp_phaser_model_change(const SpPhaser * ph){
  return ph->model_change;
}
 

int sp_phaser_init(SpPhaser * ph, SpPhasingAlgorithm * alg,Image * amplitudes, Image * support){
  if(!ph){
    return -1;
  }
  if(!alg){
    return -2;
  }
  if(!amplitudes){
    return -3;
  }
  if(!support){
    return -4;
  }
  ph->algorithm = alg;
  ph->amplitudes = amplitudes;
  ph->support = support;
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
    return -5;
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
    Image * tmp = sp_image_duplicate(ph->amplitudes,SP_COPY_ALL);
    /* randomize phases */
    sp_image_rephase(tmp,SP_RANDOM_PHASE);
    ph->model = sp_image_ifft(tmp);
    sp_image_free(tmp);
    sp_image_scale(ph->model,1.0/sp_image_size(ph->model));
  }else if(flags & SpModelZeroPhases){
    Image * tmp = sp_image_duplicate(ph->amplitudes,SP_COPY_ALL);
    sp_image_rephase(tmp,SP_ZERO_PHASE);
    ph->model = sp_image_ifft(tmp);
    sp_image_free(tmp);
    sp_image_scale(ph->model,1.0/sp_image_size(ph->model));
  }else if(flags & SpModelRandomValues){
    ph->model = sp_image_alloc(sp_image_x(ph->amplitudes),sp_image_y(ph->amplitudes),sp_image_y(ph->amplitudes));
    /* try to start with reasonable random values */
    double sum = sp_image_integrate2(ph->amplitudes);
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
      if(ph->amplitudes->mask->data[i] == 0){
	tmp->image->data[i] = sp_cinit(0,0);
      }
    }
    sp_image_free(ph->model);
    ph->model = sp_image_ifft(tmp);
    sp_image_scale(ph->model,1.0/sp_image_size(ph->model));
    sp_image_free(tmp);
  }
  ph->model_change = sp_image_alloc(sp_image_x(ph->amplitudes),sp_image_y(ph->amplitudes),sp_image_y(ph->amplitudes));
  sp_image_fill(ph->model_change,sp_cinit(0,0));
  return 0;
}


int sp_phaser_iterate(SpPhaser * ph){
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
  if(!ph->support){
    return -5;
  }
  if(!ph->model_change){
    return -6;
  }
  if(ph->algorithm->type == SpHIO){
    return phaser_iterate_hio(ph);
  }else if(ph->algorithm->type == SpRAAR){
    return phaser_iterate_raar(ph);
  }
  return -7;
}

static void phase_apply_constraints(SpPhaser * ph,Image * new_model, SpPhasingConstraints constraints){
  /* Apply constraints */
  for(int i =0;i<sp_image_size(new_model);i++){
    if(sp_real(ph->support->image->data[i])){
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

static int phaser_iterate_hio(SpPhaser * ph){
  Image * fmodel = sp_image_fft(ph->model);
  SpPhasingHIOParameters * params = ph->algorithm->params;
  sp_proj_module(fmodel,ph->amplitudes,SpInPlace);
  Image * new_model = sp_image_ifft(fmodel);
  sp_image_free(fmodel);
  sp_image_scale(new_model,1.0/sp_image_size(new_model));
  for(int i =0;i<sp_image_size(new_model);i++){
    if(sp_real(ph->support->image->data[i])){
      // Nothing to do here 
    }else{
      new_model->image->data[i] = sp_csub(ph->model->image->data[i],sp_cscale(new_model->image->data[i],params->beta));
    }
  }
  phase_apply_constraints(ph,new_model,params->constraints);
  sp_image_memcpy(ph->model_change,new_model);
  sp_image_sub(ph->model_change,ph->model);
  sp_image_free(ph->model);
  ph->model = new_model;
  return 0;
}


static int phaser_iterate_raar(SpPhaser * ph){
  Image * fmodel = sp_image_fft(ph->model);
  SpPhasingRAARParameters * params = ph->algorithm->params;
  real beta = params->beta;
  sp_proj_module(fmodel,ph->amplitudes,SpInPlace);
  Image * new_model = sp_image_ifft(fmodel);
  sp_image_free(fmodel);
  sp_image_scale(new_model,1.0/sp_image_size(new_model));
  for(int i =0;i<sp_image_size(new_model);i++){
    /* A bit of documentation about the equation:

     Rs = 2*Ps-I; Rm = 2*Pm-I

     RAAR = 1/2 * beta * (RsRm + I) + (1 - beta) * Pm;    
     RAAR = 2*beta*Ps*Pm+(1-2*beta)*Pm - beta * (Ps-I)

     Which reduces to:

     Inside the support: Pm
     Outside the support: (1 - 2*beta)*Pm + beta*I
     
    */    
    if(sp_real(ph->support->image->data[i])){
      // Nothing to do here 
    }else{
      new_model->image->data[i] = sp_cadd(sp_cscale(new_model->image->data[i],1-2*beta),sp_cscale(ph->model->image->data[i],beta));      
    }
  }
  /* Apply constraints */
  phase_apply_constraints(ph,new_model,params->constraints);
  sp_image_memcpy(ph->model_change,new_model);
  sp_image_sub(ph->model_change,ph->model);
  sp_image_free(ph->model);
  ph->model = new_model;
  return 0;
}
