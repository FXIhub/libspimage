#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include "AllTests.h"



static Image * create_test_image(int size, real oversampling,SpPhasingConstraints c);

int test_sp_phasing_common(CuTest * tc,SpPhasingAlgorithm * alg,int size, real oversampling, SpPhasingConstraints solution_constraints,
			   real beamstop,real tol){
  int attempts = 0;
  double change = 0;
  while(1){
    Image * solution = create_test_image(size,oversampling,solution_constraints);
    Image * f = sp_image_fft(solution);
    sp_image_rephase(f,SP_ZERO_PHASE);
    for(int i = 0;i<sp_image_size(f);i++){
      if(sp_image_dist(f,i,SP_TO_CORNER) < beamstop){
	f->image->data[i] = sp_cinit(0,0);    
	f->mask->data[i] = 0;    
      }else{
	f->mask->data[i] = 1;    
      }
    }
    
    Image * support = sp_image_duplicate(solution,SP_COPY_ALL);
    for(int i =0;i<sp_image_size(support);i++){
      if(sp_cabs(support->image->data[i])){
	support->image->data[i] = sp_cinit(1,0);
      }
    }
    SpPhaser * ph = sp_phaser_alloc();
    CuAssertTrue(tc,sp_phaser_init(ph,alg,NULL,SpEngineAutomatic) == 0);
    sp_phaser_set_amplitudes(ph,f);
    int i =0;
    change = 0;
    int max_iter = 300;
    CuAssertTrue(tc,sp_phaser_init_model(ph,NULL,SpModelRandomPhases) == 0); 
    CuAssertTrue(tc,sp_phaser_init_support(ph,support,0,0) == 0); 
    do{
      CuAssertTrue(tc,sp_phaser_iterate(ph,1) == 0);
      change = sp_image_integrate2(sp_phaser_model_change(ph));
      //      printf("Iter = %d Delta = %g\n",i,change);
      i++;
    }while(change > tol && i < max_iter);
    sp_phaser_free(ph);
    sp_image_free(f);
    sp_image_free(support);    
    if(change <= tol || attempts > 10){
      break;
    }
    attempts++;
  }
  if(attempts > 10){
    return -1;
  }
  return 1;
} 


int test_sp_phasing_cuda_common(CuTest * tc,SpPhasingAlgorithm * alg,int size, real oversampling, SpPhasingConstraints solution_constraints,
				real tol){
  Image * solution = create_test_image(size,oversampling,solution_constraints);
  Image * f = sp_image_fft(solution);
  sp_image_rephase(f,SP_ZERO_PHASE);
  for(int i = 0;i<sp_image_size(f);i++){
    f->mask->data[i] = 1;    
  }
    
  Image * support = sp_image_duplicate(solution,SP_COPY_ALL);
  for(int i =0;i<sp_image_size(support);i++){
    if(sp_cabs(support->image->data[i])){
      support->image->data[i] = sp_cinit(1,0);
    }
  }
  SpPhaser * ph_cpu = sp_phaser_alloc();
  SpPhaser * ph_cuda = sp_phaser_alloc();
  CuAssertTrue(tc,sp_phaser_init(ph_cpu,alg,NULL,SpEngineCPU) == 0);
  sp_phaser_set_amplitudes(ph_cpu,f);
  CuAssertTrue(tc,sp_phaser_init(ph_cuda,alg,NULL,SpEngineCUDA) == 0);
  sp_phaser_set_amplitudes(ph_cuda,f);
  int i =0;
  int max_iter = 10;
  CuAssertTrue(tc,sp_phaser_init_model(ph_cpu,NULL,SpModelRandomPhases) == 0); 
  CuAssertTrue(tc,sp_phaser_init_model(ph_cuda,sp_phaser_model(ph_cpu),0) == 0); 
  CuAssertTrue(tc,sp_phaser_init_support(ph_cpu,support,0,0) == 0); 
  CuAssertTrue(tc,sp_phaser_init_support(ph_cuda,support,0,0) == 0); 
  do{
    CuAssertTrue(tc,sp_phaser_iterate(ph_cpu,1) == 0);
    CuAssertTrue(tc,sp_phaser_iterate(ph_cuda,1) == 0);
    const Image * cpu_model = sp_phaser_model(ph_cpu);
    const Image * cuda_model = sp_phaser_model(ph_cuda);
    for(int i =0 ;i<ph_cpu->image_size;i++){
      CuAssertComplexEquals(tc,cpu_model->image->data[i],
			    cuda_model->image->data[i],sqrt(REAL_EPSILON)+
			    sqrt(REAL_EPSILON)*sp_cabs(cuda_model->image->data[i]));  
    }
    i++;
  }while(i < max_iter);
  sp_phaser_free(ph_cpu);
  sp_phaser_free(ph_cuda);
  sp_image_free(f);
  sp_image_free(solution);
  sp_image_free(support);    
  return 1;
} 

int test_sp_support_common(CuTest * tc,SpPhasingAlgorithm * alg,SpSupportArray * sup_alg, int size, real oversampling, SpPhasingConstraints solution_constraints,
			   real beamstop,real tol){
  int attempts = 0;
  double change = 0;
  const int max_attempts = 20;
  while(1){
    Image * solution = create_test_image(size,oversampling,solution_constraints);
    sp_image_write(solution,"debug_solution.h5",0);
    Image * f = sp_image_fft(solution);
    sp_image_rephase(f,SP_ZERO_PHASE);
    for(int i = 0;i<sp_image_size(f);i++){
      if(sp_image_dist(f,i,SP_TO_CORNER) < beamstop){
	f->image->data[i] = sp_cinit(0,0);    
	f->mask->data[i] = 0;    
      }else{
	f->mask->data[i] = 1;    
      }
    }
    
    Image * support = sp_image_duplicate(solution,SP_COPY_ALL);
    /* start with full support */
    sp_image_fill(support,sp_cinit(1,0));
    /*
    for(int i =0;i<sp_image_size(support);i++){
      if(sp_cabs(support->image->data[i])){
	support->image->data[i] = sp_cinit(1,0);
      }
    }
    */
    SpPhaser * ph = sp_phaser_alloc();
    CuAssertTrue(tc,sp_phaser_init(ph,alg,sup_alg,SpEngineCPU) == 0);
    sp_phaser_set_amplitudes(ph,f);
    int i =0;
    change = 0;
    int max_iter = 1000;
    CuAssertTrue(tc,sp_phaser_init_model(ph,NULL,SpModelRandomPhases) == 0); 
    CuAssertTrue(tc,sp_phaser_init_support(ph,NULL,SpSupportFromPatterson,0.004) == 0); 
    sp_image_write(sp_phaser_support(ph),"debug_init_support.h5",0);
    do{
      CuAssertTrue(tc,sp_phaser_iterate(ph,1) == 0);
      change = sp_image_integrate2(sp_phaser_model_change(ph));
      sp_image_write(sp_phaser_support(ph),"debug_support.h5",0);
      sp_image_write(sp_phaser_model(ph),"debug_model.h5",0);
      //      printf("Iter = %d Delta = %g\n",i,change);
      i++;
    }while(change > tol && i < max_iter);
    sp_phaser_free(ph);
    sp_image_free(f);
    sp_image_free(support);    
    if(change <= tol || attempts > max_attempts){
      break;
    }
    attempts++;
  }
  if(attempts > max_attempts){
    return -1;
  }
  return 1;
} 

void test_sp_support_cuda_common(CuTest * tc,SpSupportArray * sup_alg){
  int size = 16;
  real oversampling = 2;
  Image * solution = create_test_image(size,oversampling,SpNoConstraints);
  Image * f = sp_image_fft(solution);
  sp_image_rephase(f,SP_ZERO_PHASE);
  for(int i = 0;i<sp_image_size(f);i++){
    f->mask->data[i] = 1;    
  }
  sp_smap * beta = sp_smap_create_from_pair(0,0.9);
  SpPhaser * ph_cpu = sp_phaser_alloc();
  SpPhaser * ph_cuda = sp_phaser_alloc();
  SpPhasingAlgorithm * alg = sp_phasing_hio_alloc(beta,0);
  CuAssertTrue(tc,sp_phaser_init(ph_cpu,alg,sup_alg,SpEngineCPU) == 0);
  sp_phaser_set_amplitudes(ph_cpu,f);
  CuAssertTrue(tc,sp_phaser_init_model(ph_cpu,NULL,SpModelRandomPhases) == 0); 
  CuAssertTrue(tc,sp_phaser_init_support(ph_cpu,NULL,SpSupportFromPatterson,0.004) == 0); 

  CuAssertTrue(tc,sp_phaser_init(ph_cuda,alg,sup_alg,SpEngineCUDA) == 0);
  sp_phaser_set_amplitudes(ph_cuda,f);
  /* use previous model */
  Image * cpu_model = sp_image_duplicate(sp_phaser_model(ph_cpu),SP_COPY_ALL);
  CuAssertTrue(tc,sp_phaser_init_model(ph_cuda,cpu_model,0) == 0); 
  CuAssertTrue(tc,sp_phaser_init_support(ph_cuda,NULL,SpSupportFromPatterson,0.004) == 0); 
  sp_image_write(sp_phaser_support(ph_cpu),"debug_init_support_cpu.h5",0);
  sp_image_write(sp_phaser_support(ph_cuda),"debug_init_support_cuda.h5",0);
  const Image * sup_cpu = sp_phaser_support(ph_cpu);
  const Image * sup_cuda = sp_phaser_support(ph_cuda);
  int max_iter = 20;
  for(int i = 0;i<sp_image_size(f);i++){
    CuAssertComplexEquals(tc,sup_cuda->image->data[i],sup_cpu->image->data[i],sp_cabs(sup_cpu->image->data[i])*REAL_EPSILON);
  }
  real tol = REAL_EPSILON;
  int step = 1;
  for(int i = 0;i<max_iter;i+=step){

    CuAssertTrue(tc,sp_phaser_iterate(ph_cpu,step) == 0);
    CuAssertTrue(tc,sp_phaser_iterate(ph_cuda,step) == 0);
    sup_cpu = sp_phaser_support(ph_cpu);
    sup_cuda = sp_phaser_support(ph_cuda);
    sp_image_write(sup_cpu,"debug_support_cpu.h5",0);
    sp_image_write(sup_cuda,"debug_support_cuda.h5",0);
    Image * model_cpu = sp_image_duplicate(sp_phaser_model(ph_cpu),SP_COPY_ALL);
    const Image * model_cuda = sp_phaser_model(ph_cuda);
    sp_image_write(model_cpu,"debug_model_cpu.h5",0);
    sp_image_write(model_cuda,"debug_model_cuda.h5",0);
    
      /* the support will be different because of the the blurring is not binary equivalent.
	 But the area must be the same */
      //      CuAssertComplexEquals(tc,sup_cuda->image->data[j],sup_cpu->image->data[j],sp_cabs(sup_cpu->image->data[j])*REAL_EPSILON);
    CuAssertDblEquals(tc,sp_image_integrate2(sup_cuda),sp_image_integrate2(sup_cpu),REAL_EPSILON);

    sp_image_sub(model_cpu,model_cuda);
    CuAssertDblEquals(tc,sp_image_integrate2(model_cpu),0,sp_image_integrate2(model_cuda)*tol);
    /* make sure there's no cumulative divergence*/
    sp_phaser_set_model(ph_cuda,sp_phaser_model(ph_cpu));
    sp_phaser_set_support(ph_cuda,sp_phaser_support(ph_cpu));
    sp_image_free(model_cpu);

  }
  sp_phaser_free(ph_cpu);
  sp_phaser_free(ph_cuda);
  sp_image_free(f);
} 


int test_sp_phasing_speed_common(CuTest * tc,SpPhasingAlgorithm * alg,SpSupportArray * sup_alg, int size, int oversampling, SpPhasingConstraints solution_constraints,
				 int iterations,SpPhasingEngine engine){
  Image * solution = create_test_image(size,oversampling,solution_constraints);
  Image * f = sp_image_fft(solution);
  sp_image_rephase(f,SP_ZERO_PHASE);
  for(int i = 0;i<sp_image_size(f);i++){
    f->mask->data[i] = 1;    
  }
  
  Image * support = sp_image_duplicate(solution,SP_COPY_ALL);
  for(int i =0;i<sp_image_size(support);i++){
    if(sp_cabs(support->image->data[i])){
      support->image->data[i] = sp_cinit(1,0);
    }
  }
  SpPhaser * ph = sp_phaser_alloc();
  CuAssertTrue(tc,sp_phaser_init(ph,alg,sup_alg,engine) == 0);
  sp_phaser_set_amplitudes(ph,f);
  CuAssertTrue(tc,sp_phaser_init_model(ph,NULL,SpModelRandomPhases) == 0); 
  CuAssertTrue(tc,sp_phaser_init_support(ph,support,0,0) == 0); 
  int timer = sp_timer_start();
  CuAssertTrue(tc,sp_phaser_iterate(ph,iterations) == 0);
  /* retrieve result and make sure the calculations are finished */
  sp_phaser_model(ph);
  int delta_t = sp_timer_stop(timer);
  sp_phaser_free(ph);
  sp_image_free(f);
  sp_image_free(support);    

  return delta_t;
} 



int test_sp_phasing_success_common(CuTest * tc,SpPhasingAlgorithm * alg,Image * solution, real beamstop,real stop_tol,real match_tol){
  Image * f = sp_image_fft(solution);
  sp_image_rephase(f,SP_ZERO_PHASE);
  for(int i = 0;i<sp_image_size(f);i++){
    if(sp_image_dist(f,i,SP_TO_CORNER) < beamstop){
      f->image->data[i] = sp_cinit(0,0);    
      f->mask->data[i] = 0;    
    }else{
      f->mask->data[i] = 1;    
    }
  }

  Image * support = sp_image_duplicate(solution,SP_COPY_ALL);
  for(int i =0;i<sp_image_size(support);i++){
    if(sp_cabs(support->image->data[i])){
      support->image->data[i] = sp_cinit(1,0);
    }
  }
  SpPhaser * ph = sp_phaser_alloc();
  CuAssertTrue(tc,sp_phaser_init(ph,alg,NULL,SpEngineAutomatic) == 0);
  //  CuAssertTrue(tc,sp_phaser_init(ph,alg,NULL,SpEngineCPU) == 0);
  sp_phaser_set_amplitudes(ph,f);

  int i =0;
  double change = 0;
  int max_iter = 300;
  int step = 50;
  CuAssertTrue(tc,sp_phaser_init_model(ph,NULL,SpModelRandomPhases) == 0); 
  CuAssertTrue(tc,sp_phaser_init_support(ph,support,0,0) == 0); 
  do{
    CuAssertTrue(tc,sp_phaser_iterate(ph,step) == 0);
    change = sp_image_integrate2(sp_phaser_model_change(ph));
    //    printf("Iter = %d Delta = %g\n",i,change);
    i++;
  }while(change > stop_tol && i < max_iter);
  Image * model = sp_image_duplicate(sp_phaser_model(ph),SP_COPY_ALL);
  sp_image_superimpose_fractional(solution,model,SpEnantiomorph|SpCorrectPhaseShift,1);
  int match = 1;
  for(int i =0;i<sp_image_size(solution);i++){
    if(sp_real(support->image->data[i])){
      if(sp_cabs(sp_csub(model->image->data[i],solution->image->data[i])) >match_tol){
	match = 0;
	break;
      }
    }
  }
  sp_image_free(model);
  sp_phaser_free(ph);
  sp_image_free(f);
  sp_image_free(support);
  return match;
} 


int test_sp_phasing_noisy_success_common(CuTest * tc,SpPhasingAlgorithm * alg,Image * solution, real beamstop,real stop_tol,real match_tol, real photons_per_pixel,const gsl_rng * r){
  Image * f = sp_image_fft(solution);
  sp_image_rephase(f,SP_ZERO_PHASE);
  for(int i = 0;i<sp_image_size(f);i++){
    if(sp_image_dist(f,i,SP_TO_CORNER) < beamstop){
      f->image->data[i] = sp_cinit(0,0);    
      f->mask->data[i] = 0;    
    }else{
      f->mask->data[i] = 1;    
    }
  }
  /* square the image */
  for(int i = 0;i<sp_image_size(f);i++){
    sp_real(f->image->data[i]) = sp_real(f->image->data[i])*sp_real(f->image->data[i]);
  }
  sp_image_scale(f,sp_image_size(f)*photons_per_pixel/sp_cabs(sp_image_integrate(f)));
  /* Take a poisson value and the square root */
  for(int i = 0;i<sp_image_size(f);i++){
    sp_real(f->image->data[i]) = sqrt(gsl_ran_poisson(r,sp_real(f->image->data[i])));
  }


  Image * support = sp_image_duplicate(solution,SP_COPY_ALL);
  for(int i =0;i<sp_image_size(support);i++){
    if(sp_cabs(support->image->data[i])){
      support->image->data[i] = sp_cinit(1,0);
    }
  }
  sp_image_write(support,"support_noisy.h5",0);
  SpPhaser * ph = sp_phaser_alloc();
  CuAssertTrue(tc,sp_phaser_init(ph,alg,NULL,SpEngineAutomatic) == 0);
  sp_phaser_set_amplitudes(ph,f);

  int i =0;
  double change = 0;
  int max_iter = 300;
  CuAssertTrue(tc,sp_phaser_init_model(ph,NULL,SpModelRandomPhases) == 0); 
  CuAssertTrue(tc,sp_phaser_init_support(ph,support,0,0) == 0); 
  do{
    CuAssertTrue(tc,sp_phaser_iterate(ph,1) == 0);
    change = sp_image_integrate2(sp_phaser_model_change(ph));
    //      printf("Iter = %d Delta = %g\n",i,change);
    i++;
  }while(change > stop_tol && i < max_iter);
  Image * model = sp_image_duplicate(sp_phaser_model(ph),SP_COPY_ALL);
  /* Scale the solution to match the model */
  double model_support_sum = 0;
  double solution_support_sum = 0;
  for(int i = 0;i<sp_image_size(model);i++){
    if(sp_real(support->image->data[i])){
      model_support_sum += sp_cabs(model->image->data[i]);
      solution_support_sum += sp_cabs(solution->image->data[i]);
    }    
  }
  sp_image_scale(model,solution_support_sum/model_support_sum);
  sp_image_write(model,"model_noisy.h5",0);
  sp_image_write(solution,"solution_noisy.h5",0);
  sp_image_superimpose_fractional(solution,model,SpEnantiomorph|SpCorrectPhaseShift,1);
  int match = 1;
  double total_error = 0;
  for(int i =0;i<sp_image_size(solution);i++){
    if(sp_real(support->image->data[i])){
      total_error += sp_cabs2(sp_csub(model->image->data[i],solution->image->data[i]));
    }
  }
  if(total_error/sp_image_integrate2(solution) > match_tol){
    match = 0;
  }
  sp_image_free(model);
  sp_phaser_free(ph);
  sp_image_free(f);
  sp_image_free(support);
  return match;
} 


static Image * create_test_image(int size, real oversampling,SpPhasingConstraints c){
  Image * a = sp_image_alloc(size*oversampling,size*oversampling,1);
  sp_image_fill(a,sp_cinit(0,0));
  for(int x = 0;x<size;x++){
    for(int y = 0;y<size;y++){
      if(c & SpRealObject){
	sp_image_set(a,x,y,0,sp_cinit(p_drand48()-0.5,0));
      }else if(c & SpPositiveRealObject){
	sp_image_set(a,x,y,0,sp_cinit(p_drand48(),0));
      }else if(c & SpPositiveComplexObject){
	sp_image_set(a,x,y,0,sp_cinit(p_drand48(),p_drand48()));
      }else{
	sp_image_set(a,x,y,0,sp_cinit(p_drand48()-0.5,p_drand48()-0.5));
      }
    }
  }
  if(c & SpCentrosymmetricObject){
    size = size*oversampling;
    for(int x = 0;x<size;x++){
      for(int y = 0;y<size;y++){
	
	sp_image_set(a,(size-x)%size,(size-y)%size,0,sp_cconj(sp_image_get(a,x,y,0)));
      }
    }    
    sp_imag(a->image->data[0]) = 0;
  }
  a->phased = 1;
  return a;  
}

void test_sp_phasing_hio_success_rate(CuTest * tc){
  /* Simple phasing example */
  int size = 4;
  int oversampling = 2;
  sp_smap * beta = sp_smap_create_from_pair(0,0.8);
  real stop_tol = 1e-10;
  real match_tol = 1e-4;
  int nruns = 30;
  Image * pada;
  SpPhasingAlgorithm * alg;
  int n_success = 0;
  alg = sp_phasing_hio_alloc(beta,0);
  n_success = 0;
  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,0);
    n_success += test_sp_phasing_success_common(tc,alg,pada,0,stop_tol,match_tol);  
    sp_image_free(pada);
  }
#ifndef NDEBUG
  printf("HIO Complex success rate = %5.2f\n",(real)n_success/nruns);
#endif
  alg = sp_phasing_hio_alloc(beta,SpPositiveComplexObject);
  CuAssertTrue(tc,n_success>0);
  n_success = 0;
  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,SpPositiveComplexObject);  
    n_success += test_sp_phasing_success_common(tc,alg,pada,0,stop_tol,match_tol);  
    sp_image_free(pada);
  }
#ifndef NDEBUG
  printf("HIO Positive Complex success rate = %5.2f\n",(real)n_success/nruns);
#endif
  alg = sp_phasing_hio_alloc(beta,SpPositiveComplexObject| SpPositivityFlipping);
  CuAssertTrue(tc,n_success>0);
  n_success = 0;
  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,SpPositiveComplexObject);  
    n_success += test_sp_phasing_success_common(tc,alg,pada,0,stop_tol,match_tol);  
    sp_image_free(pada);
  }
#ifndef NDEBUG
  printf("HIO Flipping Positive Complex success rate = %5.2f\n",(real)n_success/nruns);
#endif
  CuAssertTrue(tc,n_success>0);
  alg = sp_phasing_hio_alloc(beta,SpRealObject);
  n_success = 0;
  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,SpRealObject);  
    n_success += test_sp_phasing_success_common(tc,alg,pada,0,stop_tol,match_tol);  
    sp_image_free(pada);
  }
#ifndef NDEBUG
  printf("HIO Real success rate = %5.2f\n",(real)n_success/nruns);
#endif
  CuAssertTrue(tc,n_success>0);
  n_success = 0;
  alg = sp_phasing_hio_alloc(beta,SpPositiveRealObject);
  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,SpPositiveRealObject);
    n_success += test_sp_phasing_success_common(tc,alg,pada,0,stop_tol,match_tol);  
    sp_image_free(pada); 
  }
#ifndef NDEBUG  
  printf("HIO Positive Real success rate = %5.2f\n",(real)n_success/nruns);
#endif
  CuAssertTrue(tc,n_success>0);
  n_success = 0;
  alg = sp_phasing_hio_alloc(beta,SpPositiveRealObject| SpPositivityFlipping);
  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,SpPositiveRealObject);
    n_success += test_sp_phasing_success_common(tc,alg,pada,0,stop_tol,match_tol);  
    sp_image_free(pada); 
  }
#ifndef NDEBUG
  printf("HIO Flipping Positive Real success rate = %5.2f\n",(real)n_success/nruns);
#endif
  CuAssertTrue(tc,n_success>0);
  PRINT_DONE;
}



void test_sp_phasing_hio_noisy_success_rate(CuTest * tc){
  /* Simple phasing example */
  int size = 4;
  int oversampling = 2;
  sp_smap * beta = sp_smap_create_from_pair(0,0.8);
  real stop_tol = 1e-10;
  real match_tol = 1e-2;
  //  real photons_per_pixel = 1/(match_tol*match_tol*match_tol);
  real photons_per_pixel = 1/(match_tol*match_tol);
  int nruns = 100;
  Image * pada;
  SpPhasingAlgorithm * alg;
  gsl_rng * r = gsl_rng_alloc(gsl_rng_taus);
  int n_success = 0;
  real beamstop = 0;
  alg = sp_phasing_hio_alloc(beta,0);
  n_success = 0;
  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,0);
    n_success += test_sp_phasing_noisy_success_common(tc,alg,pada,beamstop,stop_tol,match_tol,photons_per_pixel,r);
    sp_image_free(pada);
  }
#ifndef NDEBUG
  printf("HIO Complex success rate with %g photons per pixel = %5.2f\n",photons_per_pixel,(real)n_success/nruns);
#endif

  alg = sp_phasing_hio_alloc(beta,0);
  CuAssertTrue(tc,n_success>0);
  n_success = 0;
  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,SpPositiveRealObject);
    n_success += test_sp_phasing_noisy_success_common(tc,alg,pada,beamstop,stop_tol,match_tol,photons_per_pixel,r);  
    sp_image_free(pada);
  }
#ifndef NDEBUG
  printf("HIO Positive Complex object no constraints success rate with %g photons per pixel = %5.2f\n",photons_per_pixel,(real)n_success/nruns);
#endif

  alg = sp_phasing_hio_alloc(beta,SpPositiveComplexObject);
  CuAssertTrue(tc,n_success>0);
  n_success = 0;
  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,SpPositiveComplexObject);  
    n_success += test_sp_phasing_noisy_success_common(tc,alg,pada,beamstop,stop_tol,match_tol,photons_per_pixel,r);  
    sp_image_free(pada);
  }
#ifndef NDEBUG
  printf("HIO Positive Complex success rate with %g photons per pixel = %5.2f\n",photons_per_pixel,(real)n_success/nruns);
#endif


  alg = sp_phasing_hio_alloc(beta,SpPositiveComplexObject| SpPositivityFlipping);
  CuAssertTrue(tc,n_success>0);
  n_success = 0;
  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,SpPositiveComplexObject);  
    n_success += test_sp_phasing_noisy_success_common(tc,alg,pada,beamstop,stop_tol,match_tol,photons_per_pixel,r);  
    sp_image_free(pada);
  }
#ifndef NDEBUG
  printf("HIO Flipping Positive Complex success rate with %g photons per pixel = %5.2f\n",photons_per_pixel,(real)n_success/nruns);
#endif
  CuAssertTrue(tc,n_success>0);
  alg = sp_phasing_hio_alloc(beta,SpRealObject);
  n_success = 0;
  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,SpRealObject);  
    n_success += test_sp_phasing_noisy_success_common(tc,alg,pada,beamstop,stop_tol,match_tol,photons_per_pixel,r);  
    sp_image_free(pada);
  }
#ifndef NDEBUG
  printf("HIO Real success rate with %g photons per pixel = %5.2f\n",photons_per_pixel,(real)n_success/nruns);
#endif
  CuAssertTrue(tc,n_success>0);
  n_success = 0;
  alg = sp_phasing_hio_alloc(beta,SpPositiveRealObject);
  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,SpPositiveRealObject);
    n_success += test_sp_phasing_noisy_success_common(tc,alg,pada,beamstop,stop_tol,match_tol,photons_per_pixel,r);  
    sp_image_free(pada); 
  }
#ifndef NDEBUG  
  printf("HIO Positive Real success rate with %g photons per pixel = %5.2f\n",photons_per_pixel,(real)n_success/nruns);
#endif
  CuAssertTrue(tc,n_success>0);
  n_success = 0;
  alg = sp_phasing_hio_alloc(beta,SpPositiveRealObject| SpPositivityFlipping);
  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,SpPositiveRealObject);
    n_success += test_sp_phasing_noisy_success_common(tc,alg,pada,beamstop,stop_tol,match_tol,photons_per_pixel,r);  
    sp_image_free(pada); 
  }
#ifndef NDEBUG
  printf("HIO Flipping Positive Real success rate with %g photons per pixel = %5.2f\n",photons_per_pixel,(real)n_success/nruns);
#endif
  CuAssertTrue(tc,n_success>0);

  beamstop = 1;

  n_success = 0;
  alg = sp_phasing_hio_alloc(beta,SpRealObject);

  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,SpRealObject);
    n_success += test_sp_phasing_noisy_success_common(tc,alg,pada,beamstop,stop_tol,match_tol,photons_per_pixel,r);  
    sp_image_free(pada); 
  }
#ifndef NDEBUG
  printf("HIO Real success rate with %g photons per pixel\nand beamstop radius %g = %5.2f\n",photons_per_pixel,beamstop,(real)n_success/nruns);
#endif
  CuAssertTrue(tc,n_success>0);

  n_success = 0;
  alg = sp_phasing_hio_alloc(beta,0);

  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,0);
    n_success += test_sp_phasing_noisy_success_common(tc,alg,pada,beamstop,stop_tol,match_tol,photons_per_pixel,r);  
    sp_image_free(pada); 
  }
#ifndef NDEBUG
  printf("HIO Complex success rate with %g photons per pixel\nand beamstop radius %g = %5.2f\n",photons_per_pixel,beamstop,(real)n_success/nruns);
#endif
  CuAssertTrue(tc,n_success>0);
  PRINT_DONE;
}


void test_sp_phasing_raar_success_rate(CuTest * tc){
  /* Simple phasing example */
  int size = 4;
  int oversampling = 2;
  sp_smap * beta = sp_smap_create_from_pair(0,0.8);
  real stop_tol = 1e-10;
  real match_tol = 1e-4;
  int nruns = 30;
  Image * pada;
  SpPhasingAlgorithm * alg;
  int n_success = 0;
  alg = sp_phasing_raar_alloc(beta,0);
  n_success = 0;
  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,0);
    n_success += test_sp_phasing_success_common(tc,alg,pada,0,stop_tol,match_tol);  
    sp_image_free(pada);
  }
#ifndef NDEBUG
  printf("RAAR Complex success rate = %5.2f\n",(real)n_success/nruns);
#endif
  alg = sp_phasing_raar_alloc(beta,SpPositiveComplexObject);
  CuAssertTrue(tc,n_success>0);
  n_success = 0;
  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,SpPositiveComplexObject);  
    n_success += test_sp_phasing_success_common(tc,alg,pada,0,stop_tol,match_tol);  
    sp_image_free(pada);
  }
#ifndef NDEBUG
  printf("RAAR Positive Complex success rate = %5.2f\n",(real)n_success/nruns);
#endif
  alg = sp_phasing_raar_alloc(beta,SpPositiveComplexObject| SpPositivityFlipping);
  CuAssertTrue(tc,n_success>0);
  n_success = 0;
  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,SpPositiveComplexObject);  
    n_success += test_sp_phasing_success_common(tc,alg,pada,0,stop_tol,match_tol);  
    sp_image_free(pada);
  }
#ifndef NDEBUG
  printf("RAAR Flipping Positive Complex success rate = %5.2f\n",(real)n_success/nruns);
#endif
  alg = sp_phasing_raar_alloc(beta,SpRealObject);
  CuAssertTrue(tc,n_success>0);
  n_success = 0;
  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,SpRealObject);  
    n_success += test_sp_phasing_success_common(tc,alg,pada,0,stop_tol,match_tol);  
    sp_image_free(pada);
  }
#ifndef NDEBUG
  printf("RAAR Real success rate = %5.2f\n",(real)n_success/nruns);
#endif
  CuAssertTrue(tc,n_success>0);

  n_success = 0;
  alg = sp_phasing_raar_alloc(beta,SpPositiveRealObject);
  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,SpPositiveRealObject);
    n_success += test_sp_phasing_success_common(tc,alg,pada,0,stop_tol,match_tol);  
    sp_image_free(pada); 
  }
#ifndef NDEBUG
  printf("RAAR Positive Real success rate = %5.2f\n",(real)n_success/nruns);
#endif
  CuAssertTrue(tc,n_success>0);

  n_success = 0;
  alg = sp_phasing_raar_alloc(beta,SpPositiveRealObject| SpPositivityFlipping);
  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,SpPositiveRealObject);
    n_success += test_sp_phasing_success_common(tc,alg,pada,0,stop_tol,match_tol);  
    sp_image_free(pada); 
  }
#ifndef NDEBUG
  printf("RAAR Flipping Positive Real success rate = %5.2f\n",(real)n_success/nruns);
#endif
  CuAssertTrue(tc,n_success>0);
  PRINT_DONE;
}


void test_sp_phasing_raar_noisy_success_rate(CuTest * tc){
  /* Simple phasing example */
  int size = 4;
  int oversampling = 2;
  sp_smap * beta = sp_smap_create_from_pair(0,0.8);
  real stop_tol = 1e-10;
  real match_tol = 1e-2;
  real photons_per_pixel = 1/(match_tol*match_tol*match_tol);
  int nruns = 30;
  Image * pada;
  SpPhasingAlgorithm * alg;
  gsl_rng * r = gsl_rng_alloc(gsl_rng_taus);
  int n_success = 0;
  real beamstop =0;
  alg = sp_phasing_raar_alloc(beta,0);
  n_success = 0;
  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,0);
    n_success += test_sp_phasing_noisy_success_common(tc,alg,pada,0,stop_tol,match_tol,photons_per_pixel,r);
    sp_image_free(pada);
  }
#ifndef NDEBUG
  printf("RAAR Complex success rate with %g photons per pixel = %5.2f\n",photons_per_pixel,(real)n_success/nruns);
#endif
  alg = sp_phasing_raar_alloc(beta,SpPositiveComplexObject);
  CuAssertTrue(tc,n_success>0);
  n_success = 0;
  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,SpPositiveComplexObject);  
    n_success += test_sp_phasing_noisy_success_common(tc,alg,pada,0,stop_tol,match_tol,photons_per_pixel,r);  
    sp_image_free(pada);
  }
#ifndef NDEBUG
  printf("RAAR Positive Complex success rate with %g photons per pixel = %5.2f\n",photons_per_pixel,(real)n_success/nruns);
#endif
  alg = sp_phasing_raar_alloc(beta,SpPositiveComplexObject| SpPositivityFlipping);
  CuAssertTrue(tc,n_success>0);
  n_success = 0;
  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,SpPositiveComplexObject);  
    n_success += test_sp_phasing_noisy_success_common(tc,alg,pada,0,stop_tol,match_tol,photons_per_pixel,r);  
    sp_image_free(pada);
  }
#ifndef NDEBUG
  printf("RAAR Flipping Positive Complex success rate with %g photons per pixel = %5.2f\n",photons_per_pixel,(real)n_success/nruns);
#endif
  CuAssertTrue(tc,n_success>0);
  alg = sp_phasing_raar_alloc(beta,SpRealObject);
  n_success = 0;
  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,SpRealObject);  
    n_success += test_sp_phasing_noisy_success_common(tc,alg,pada,0,stop_tol,match_tol,photons_per_pixel,r);  
    sp_image_free(pada);
  }
#ifndef NDEBUG
  printf("RAAR Real success rate with %g photons per pixel = %5.2f\n",photons_per_pixel,(real)n_success/nruns);
#endif
  CuAssertTrue(tc,n_success>0);
  n_success = 0;
  alg = sp_phasing_raar_alloc(beta,SpPositiveRealObject);
  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,SpPositiveRealObject);
    n_success += test_sp_phasing_noisy_success_common(tc,alg,pada,0,stop_tol,match_tol,photons_per_pixel,r);  
    sp_image_free(pada); 
  }
#ifndef NDEBUG  
  printf("RAAR Positive Real success rate with %g photons per pixel = %5.2f\n",photons_per_pixel,(real)n_success/nruns);
#endif
  CuAssertTrue(tc,n_success>0);
  n_success = 0;
  alg = sp_phasing_raar_alloc(beta,SpPositiveRealObject| SpPositivityFlipping);
  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,SpPositiveRealObject);
    n_success += test_sp_phasing_noisy_success_common(tc,alg,pada,0,stop_tol,match_tol,photons_per_pixel,r);  
    sp_image_free(pada); 
  }
#ifndef NDEBUG
  printf("RAAR Flipping Positive Real success rate with %g photons per pixel = %5.2f\n",photons_per_pixel,(real)n_success/nruns);
#endif
  CuAssertTrue(tc,n_success>0);

  beamstop = 1;

  n_success = 0;
  alg = sp_phasing_raar_alloc(beta,SpRealObject);

  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,SpRealObject);
    n_success += test_sp_phasing_noisy_success_common(tc,alg,pada,beamstop,stop_tol,match_tol,photons_per_pixel,r);  
    sp_image_free(pada); 
  }
#ifndef NDEBUG
  printf("RAAR Real success rate with %g photons per pixel\nand beamstop radius %g = %5.2f\n",photons_per_pixel,beamstop,(real)n_success/nruns);
#endif
  CuAssertTrue(tc,n_success>0);

  n_success = 0;
  alg = sp_phasing_raar_alloc(beta,0);

  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,0);
    n_success += test_sp_phasing_noisy_success_common(tc,alg,pada,beamstop,stop_tol,match_tol,photons_per_pixel,r);  
    sp_image_free(pada); 
  }
#ifndef NDEBUG
  printf("RAAR Complex success rate with %g photons per pixel\nand beamstop radius %g = %5.2f\n",photons_per_pixel,beamstop,(real)n_success/nruns);
#endif
  CuAssertTrue(tc,n_success>0);
  PRINT_DONE;
}


void test_sp_phasing_diff_map_success_rate(CuTest * tc){
  /* Simple phasing example */
  int size = 4;
  int oversampling = 2;
  sp_smap * beta = sp_smap_create_from_pair(0,0.8);
  real stop_tol = 1e-10;
  real match_tol = 1e-4;
  int nruns = 30;
  Image * pada;
  SpPhasingAlgorithm * alg;
  int n_success = 0;
  alg = sp_phasing_diff_map_alloc(beta,INFINITY,INFINITY,0);
  n_success = 0;
  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,0);
    n_success += test_sp_phasing_success_common(tc,alg,pada,0,stop_tol,match_tol);  
    sp_image_free(pada);
  }
#ifndef NDEBUG
  printf("DIFF_MAP Complex success rate = %5.2f\n",(real)n_success/nruns);
#endif
  alg = sp_phasing_diff_map_alloc(beta,INFINITY,INFINITY,SpPositiveComplexObject);
  CuAssertTrue(tc,n_success>0);
  n_success = 0;
  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,SpPositiveComplexObject);  
    n_success += test_sp_phasing_success_common(tc,alg,pada,0,stop_tol,match_tol);  
    sp_image_free(pada);
  }
#ifndef NDEBUG
  printf("DIFF_MAP Positive Complex success rate = %5.2f\n",(real)n_success/nruns);
#endif
  alg = sp_phasing_diff_map_alloc(beta,INFINITY,INFINITY,SpPositiveComplexObject| SpPositivityFlipping);
  CuAssertTrue(tc,n_success>0);
  n_success = 0;
  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,SpPositiveComplexObject);  
    n_success += test_sp_phasing_success_common(tc,alg,pada,0,stop_tol,match_tol);  
    sp_image_free(pada);
  }
#ifndef NDEBUG
  printf("DIFF_MAP Flipping Positive Complex success rate = %5.2f\n",(real)n_success/nruns);
#endif
  alg = sp_phasing_diff_map_alloc(beta,INFINITY,INFINITY,SpRealObject);
  CuAssertTrue(tc,n_success>0);
  n_success = 0;
  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,SpRealObject);  
    n_success += test_sp_phasing_success_common(tc,alg,pada,0,stop_tol,match_tol);  
    sp_image_free(pada);
  }
#ifndef NDEBUG
  printf("DIFF_MAP Real success rate = %5.2f\n",(real)n_success/nruns);
#endif
  CuAssertTrue(tc,n_success>0);

  n_success = 0;
  alg = sp_phasing_diff_map_alloc(beta,INFINITY,INFINITY,SpPositiveRealObject);
  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,SpPositiveRealObject);
    n_success += test_sp_phasing_success_common(tc,alg,pada,0,stop_tol,match_tol);  
    sp_image_free(pada); 
  }
#ifndef NDEBUG
  printf("DIFF_MAP Positive Real success rate = %5.2f\n",(real)n_success/nruns);
#endif
  CuAssertTrue(tc,n_success>0);

  n_success = 0;
  alg = sp_phasing_diff_map_alloc(beta,INFINITY,INFINITY,SpPositiveRealObject| SpPositivityFlipping);
  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,SpPositiveRealObject);
    n_success += test_sp_phasing_success_common(tc,alg,pada,0,stop_tol,match_tol);  
    sp_image_free(pada); 
  }
#ifndef NDEBUG
  printf("DIFF_MAP Flipping Positive Real success rate = %5.2f\n",(real)n_success/nruns);
#endif
  CuAssertTrue(tc,n_success>0);
  PRINT_DONE;
}

void test_sp_phasing_diff_map_noisy_success_rate(CuTest * tc){
  /* Simple phasing example */
  int size = 4;
  int oversampling = 2;
  sp_smap * beta = sp_smap_create_from_pair(0,0.8);
  real stop_tol = 1e-10;
  real match_tol = 1e-2;
  real photons_per_pixel = 1/(match_tol*match_tol*match_tol);
  int nruns = 30;
  Image * pada;
  SpPhasingAlgorithm * alg;
  gsl_rng * r = gsl_rng_alloc(gsl_rng_taus);
  int n_success = 0;
  real beamstop =0;
  alg = sp_phasing_diff_map_alloc(beta,INFINITY,INFINITY,0);
  n_success = 0;
  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,0);
    n_success += test_sp_phasing_noisy_success_common(tc,alg,pada,0,stop_tol,match_tol,photons_per_pixel,r);
    sp_image_free(pada);
  }
#ifndef NDEBUG
  printf("DIFF MAP Complex success rate with %g photons per pixel = %5.2f\n",photons_per_pixel,(real)n_success/nruns);
#endif
  alg = sp_phasing_diff_map_alloc(beta,INFINITY,INFINITY,SpPositiveComplexObject);
  CuAssertTrue(tc,n_success>0);
  n_success = 0;
  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,SpPositiveComplexObject);  
    n_success += test_sp_phasing_noisy_success_common(tc,alg,pada,0,stop_tol,match_tol,photons_per_pixel,r);  
    sp_image_free(pada);
  }
#ifndef NDEBUG
  printf("DIFF MAP Positive Complex success rate with %g photons per pixel = %5.2f\n",photons_per_pixel,(real)n_success/nruns);
#endif
  alg = sp_phasing_diff_map_alloc(beta,INFINITY,INFINITY,SpPositiveComplexObject| SpPositivityFlipping);
  CuAssertTrue(tc,n_success>0);
  n_success = 0;
  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,SpPositiveComplexObject);  
    n_success += test_sp_phasing_noisy_success_common(tc,alg,pada,0,stop_tol,match_tol,photons_per_pixel,r);  
    sp_image_free(pada);
  }
#ifndef NDEBUG
  printf("DIFF MAP Flipping Positive Complex success rate with %g photons per pixel = %5.2f\n",photons_per_pixel,(real)n_success/nruns);
#endif
  CuAssertTrue(tc,n_success>0);
  alg = sp_phasing_diff_map_alloc(beta,INFINITY,INFINITY,SpRealObject);
  n_success = 0;
  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,SpRealObject);  
    n_success += test_sp_phasing_noisy_success_common(tc,alg,pada,0,stop_tol,match_tol,photons_per_pixel,r);  
    sp_image_free(pada);
  }
#ifndef NDEBUG
  printf("DIFF MAP Real success rate with %g photons per pixel = %5.2f\n",photons_per_pixel,(real)n_success/nruns);
#endif
  CuAssertTrue(tc,n_success>0);
  n_success = 0;
  alg = sp_phasing_diff_map_alloc(beta,INFINITY,INFINITY,SpPositiveRealObject);
  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,SpPositiveRealObject);
    n_success += test_sp_phasing_noisy_success_common(tc,alg,pada,0,stop_tol,match_tol,photons_per_pixel,r);  
    sp_image_free(pada); 
  }
#ifndef NDEBUG  
  printf("DIFF MAP Positive Real success rate with %g photons per pixel = %5.2f\n",photons_per_pixel,(real)n_success/nruns);
#endif
  CuAssertTrue(tc,n_success>0);
  n_success = 0;
  alg = sp_phasing_diff_map_alloc(beta,INFINITY,INFINITY,SpPositiveRealObject| SpPositivityFlipping);
  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,SpPositiveRealObject);
    n_success += test_sp_phasing_noisy_success_common(tc,alg,pada,0,stop_tol,match_tol,photons_per_pixel,r);  
    sp_image_free(pada); 
  }
#ifndef NDEBUG
  printf("DIFF MAP Flipping Positive Real success rate with %g photons per pixel = %5.2f\n",photons_per_pixel,(real)n_success/nruns);
#endif
  CuAssertTrue(tc,n_success>0);

  beamstop = 1;

  n_success = 0;
  alg = sp_phasing_diff_map_alloc(beta,INFINITY,INFINITY,SpRealObject);

  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,SpRealObject);
    n_success += test_sp_phasing_noisy_success_common(tc,alg,pada,beamstop,stop_tol,match_tol,photons_per_pixel,r);  
    sp_image_free(pada); 
  }
#ifndef NDEBUG
  printf("DIFF MAP Real success rate with %g photons per pixel\nand beamstop radius %g = %5.2f\n",photons_per_pixel,beamstop,(real)n_success/nruns);
#endif
  CuAssertTrue(tc,n_success>0);

  n_success = 0;
  alg = sp_phasing_diff_map_alloc(beta,INFINITY,INFINITY,0);

  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,0);
    n_success += test_sp_phasing_noisy_success_common(tc,alg,pada,beamstop,stop_tol,match_tol,photons_per_pixel,r);  
    sp_image_free(pada); 
  }
#ifndef NDEBUG
  printf("DIFF MAP Complex success rate with %g photons per pixel\nand beamstop radius %g = %5.2f\n",photons_per_pixel,beamstop,(real)n_success/nruns);
#endif
  CuAssertTrue(tc,n_success>0);
  PRINT_DONE;
}

void test_sp_phasing_hio(CuTest * tc){
  /* Simple phasing example */
  int size = 4;
  int oversampling = 2;
  sp_smap * beta = sp_smap_create_from_pair(0,0.8);
  real tol = 1e-10;
  SpPhasingAlgorithm * alg = sp_phasing_hio_alloc(beta,0);
  //  CuAssertIntEquals(tc,test_sp_phasing_common(tc,alg,size,oversampling,SpNoConstraints,0,tol),1);

  /*  if(sp_cuda_get_device_type() == SpCUDAHardwareDevice){
    CuAssertIntEquals(tc,test_sp_phasing_cuda_common(tc,alg,size,oversampling,SpNoConstraints,tol),1);
    }
  alg = sp_phasing_hio_alloc(beta,SpPositiveRealObject);
  CuAssertIntEquals(tc,test_sp_phasing_common(tc,alg,size,oversampling,SpPositiveRealObject,0,tol),1);
  alg = sp_phasing_hio_alloc(beta,SpRealObject);
  CuAssertIntEquals(tc,test_sp_phasing_common(tc,alg,size,oversampling,SpRealObject,0,tol),1);
  alg = sp_phasing_hio_alloc(beta,SpPositiveComplexObject);
  CuAssertIntEquals(tc,test_sp_phasing_common(tc,alg,size,oversampling,SpPositiveComplexObject,0,tol),1);
  alg = sp_phasing_hio_alloc(beta,SpPositiveComplexObject);
  CuAssertIntEquals(tc,test_sp_phasing_common(tc,alg,size,oversampling,SpPositiveComplexObject,1,tol),1);
  alg = sp_phasing_hio_alloc(beta,0);
  CuAssertIntEquals(tc,test_sp_phasing_common(tc,alg,size,oversampling,SpNoConstraints,1,tol),1);
  alg = sp_phasing_hio_alloc(beta,SpPositiveRealObject);
  CuAssertIntEquals(tc,test_sp_phasing_common(tc,alg,size,oversampling,SpPositiveRealObject,1,tol),1);
  */

  /* Also try some large sizes */
  if(sp_cuda_get_device_type() == SpCUDAHardwareDevice){
    int size = 1024;
  int oversampling = 2;
    CuAssertIntEquals(tc,test_sp_phasing_cuda_common(tc,alg,size,oversampling,SpNoConstraints,tol),1);
  }
  PRINT_DONE;
}


void test_sp_support_speed(CuTest * tc){
#ifndef NDEBUG
  /* Simple phasing example */
  int size = 512;
  int oversampling = 2;
  sp_smap * beta = sp_smap_create_from_pair(0,0.8);
  SpPhasingAlgorithm * alg = sp_phasing_hio_alloc(beta,0);
  int iterations = 500;
  sp_smap * blur_radius = sp_smap_alloc(2);
  sp_smap_insert(blur_radius,0,3);
  sp_smap_insert(blur_radius,2000,0.7);
  sp_smap * threshold = sp_smap_alloc(1);
  sp_smap_insert(threshold,0,0.15);  
  sp_smap * area = sp_smap_alloc(1);
  sp_smap_insert(area,0,1.3*1.0/(oversampling*oversampling));  
  SpSupportArray * sup_alg;
  if(sp_cuda_get_device_type() == SpCUDAHardwareDevice){
    sup_alg = sp_support_array_init(sp_support_threshold_alloc(blur_radius,threshold),1);
    int delta_t = test_sp_phasing_speed_common(tc,alg,sup_alg,size,oversampling,SpNoConstraints,iterations,SpEngineCUDA);
    printf("CUDA %dx%d threshold support = %g iterations per second\n",size*oversampling,size*oversampling,(1.0e6*iterations)/delta_t);
    sup_alg = sp_support_array_init(sp_support_area_alloc(blur_radius,area),1);
    delta_t = test_sp_phasing_speed_common(tc,alg,sup_alg,size,oversampling,SpNoConstraints,iterations,SpEngineCUDA);
    printf("CUDA %dx%d area support = %g iterations per second\n",size*oversampling,size*oversampling,(1.0e6*iterations)/delta_t);
  }
  iterations = 4;
  sup_alg = sp_support_array_init(sp_support_threshold_alloc(blur_radius,threshold),1);
  int delta_t = test_sp_phasing_speed_common(tc,alg,sup_alg,size,oversampling,SpNoConstraints,iterations,SpEngineCPU);
  printf("CPU %dx%d with threshold support = %g iterations per second\n",size*oversampling,size*oversampling,(1.0e6*iterations)/delta_t);
  sup_alg = sp_support_array_init(sp_support_area_alloc(blur_radius,area),1);
  delta_t = test_sp_phasing_speed_common(tc,alg,sup_alg,size,oversampling,SpNoConstraints,iterations,SpEngineCPU);
  printf("CPU %dx%d with area support = %g iterations per second\n",size*oversampling,size*oversampling,(1.0e6*iterations)/delta_t);
#endif
  PRINT_DONE;
}

void test_sp_phasing_hio_speed(CuTest * tc){
#ifndef NDEBUG
  /* Simple phasing example */
  int size = 512;
  int oversampling = 2;
  sp_smap * beta = sp_smap_create_from_pair(0,0.8);
  SpPhasingAlgorithm * alg = sp_phasing_hio_alloc(beta,0);
  int iterations = 2000;
  sp_smap * blur_radius = sp_smap_alloc(2);
  sp_smap_insert(blur_radius,0,3);
  sp_smap_insert(blur_radius,2000,0.7);
  sp_smap * threshold = sp_smap_alloc(1);
  sp_smap_insert(threshold,0,0.15);  
  sp_smap * area = sp_smap_alloc(1);
  sp_smap_insert(area,0,1.3*1.0/(oversampling*oversampling));  
  SpSupportArray * sup_alg;

  if(sp_cuda_get_device_type() == SpCUDAHardwareDevice){
    int delta_t = test_sp_phasing_speed_common(tc,alg,NULL,size,oversampling,SpNoConstraints,iterations,SpEngineCUDA);
    printf("CUDA HIO %dx%d = %g iterations per second\n",size*oversampling,size*oversampling,(1.0e6*iterations)/delta_t);
    sup_alg = sp_support_array_init(sp_support_threshold_alloc(blur_radius,threshold),20);
    delta_t = test_sp_phasing_speed_common(tc,alg,sup_alg,size,oversampling,SpNoConstraints,iterations,SpEngineCUDA);
    printf("CUDA HIO %dx%d with threshold support every 20 = %g iterations per second\n",size*oversampling,size*oversampling,(1.0e6*iterations)/delta_t);
    sup_alg = sp_support_array_init(sp_support_area_alloc(blur_radius,area),20);
    delta_t = test_sp_phasing_speed_common(tc,alg,sup_alg,size,oversampling,SpNoConstraints,iterations,SpEngineCUDA);
    printf("CUDA HIO %dx%d with area support every 20 = %g iterations per second\n",size*oversampling,size*oversampling,(1.0e6*iterations)/delta_t);




  }
  iterations = 10;
  int delta_t = test_sp_phasing_speed_common(tc,alg,NULL,size,oversampling,SpNoConstraints,iterations,SpEngineCPU);
  printf("CPU HIO %dx%d = %g iterations per second\n",size*oversampling,size*oversampling,(1.0e6*iterations)/delta_t);
  sup_alg = sp_support_array_init(sp_support_threshold_alloc(blur_radius,threshold),20);
  delta_t = test_sp_phasing_speed_common(tc,alg,sup_alg,size,oversampling,SpNoConstraints,iterations,SpEngineCPU);
  printf("CPU HIO %dx%d with threshold support every 20 = %g iterations per second\n",size*oversampling,size*oversampling,(1.0e6*iterations)/delta_t);
  sup_alg = sp_support_array_init(sp_support_area_alloc(blur_radius,area),20);
  delta_t = test_sp_phasing_speed_common(tc,alg,sup_alg,size,oversampling,SpNoConstraints,iterations,SpEngineCPU);
  printf("CPU HIO %dx%d with area support every 20 = %g iterations per second\n",size*oversampling,size*oversampling,(1.0e6*iterations)/delta_t);
#endif
  PRINT_DONE;
}


void test_sp_phasing_hio_speed_by_size(CuTest * tc){
#ifndef NDEBUG
  for(int size = 64;size<2000;size*=2){
  /* Simple phasing example */
    int oversampling = 2;
    sp_smap * beta = sp_smap_create_from_pair(0,0.8);
    SpPhasingAlgorithm * alg = sp_phasing_raar_alloc(beta,0);
    int iterations = 10*1024*1024/(size*size);
    sp_smap * blur_radius = sp_smap_alloc(2);
    sp_smap_insert(blur_radius,0,3);
    sp_smap_insert(blur_radius,2000,0.7);
    sp_smap * threshold = sp_smap_alloc(1);
    sp_smap_insert(threshold,0,0.15);  
    sp_smap * area = sp_smap_alloc(1);
    sp_smap_insert(area,0,1.3*1.0/(oversampling*oversampling));  
    
    if(sp_cuda_get_device_type() == SpCUDAHardwareDevice){
      //      int delta_t = test_sp_phasing_speed_common(tc,alg,NULL,size,oversampling,SpNoConstraints,iterations,SpEngineCUDA);
      int delta_t = test_sp_phasing_speed_common(tc,alg,NULL,size,oversampling,SpNoConstraints,iterations,SpEngineCPU);
      printf("CUDA HIO %dx%d = %g iterations per second\n",size*oversampling,size*oversampling,(1.0e6*iterations)/delta_t);
    }
  }
#endif
  PRINT_DONE;
}



void test_sp_phasing_raar(CuTest * tc){
  /* Simple phasing example */
  int size = 4;
  int oversampling = 2;
  real tol = 1e-10;
  sp_smap * beta = sp_smap_create_from_pair(0,0.8);
  SpPhasingAlgorithm * alg = sp_phasing_raar_alloc(beta,0);
  CuAssertIntEquals(tc,test_sp_phasing_common(tc,alg,size,oversampling,SpNoConstraints,0,tol),1);
  alg = sp_phasing_raar_alloc(beta,SpPositiveRealObject);
  CuAssertIntEquals(tc,test_sp_phasing_common(tc,alg,size,oversampling,SpPositiveRealObject,0,tol),1);
  alg = sp_phasing_raar_alloc(beta,SpRealObject);
  CuAssertIntEquals(tc,test_sp_phasing_common(tc,alg,size,oversampling,SpRealObject,0,tol),1);
  alg = sp_phasing_raar_alloc(beta,SpPositiveComplexObject);
  CuAssertIntEquals(tc,test_sp_phasing_common(tc,alg,size,oversampling,SpPositiveComplexObject,0,tol),1);
  alg = sp_phasing_raar_alloc(beta,SpPositiveComplexObject);
  CuAssertIntEquals(tc,test_sp_phasing_common(tc,alg,size,oversampling,SpPositiveComplexObject,1,tol),1);
  alg = sp_phasing_raar_alloc(beta,SpPositiveRealObject);
  CuAssertIntEquals(tc,test_sp_phasing_common(tc,alg,size,oversampling,SpPositiveRealObject,1,tol),1);
  PRINT_DONE;
}

void test_sp_phasing_raar_speed(CuTest * tc){
#ifndef NDEBUG
  /* Simple phasing example */
  int size = 512;
  int oversampling = 2;
  sp_smap * beta = sp_smap_create_from_pair(0,0.8);
  SpPhasingAlgorithm * alg = sp_phasing_raar_alloc(beta,0);
  int iterations = 2000;
  if(sp_cuda_get_device_type() == SpCUDAHardwareDevice){
    int delta_t = test_sp_phasing_speed_common(tc,alg,NULL,size,oversampling,SpNoConstraints,iterations,SpEngineCUDA);
    printf("CUDA RAAR %dx%d = %g iterations per second\n",size*oversampling,size*oversampling,(1.0e6*iterations)/delta_t);
  }
  iterations = 10;
  int delta_t = test_sp_phasing_speed_common(tc,alg,NULL,size,oversampling,SpNoConstraints,iterations,SpEngineCPU);
  printf("CPU RAAR %dx%d = %g iterations per second\n",size*oversampling,size*oversampling,(1.0e6*iterations)/delta_t);
#endif
  PRINT_DONE;
}


void test_sp_phasing_diff_map_speed(CuTest * tc){
#ifndef NDEBUG
  /* Simple phasing example */
  int size = 512;
  int oversampling = 2;
  sp_smap * beta = sp_smap_create_from_pair(0,0.8);
  SpPhasingAlgorithm * alg = sp_phasing_diff_map_alloc(beta,INFINITY,INFINITY,0);
  int iterations = 2000;
  if(sp_cuda_get_device_type() == SpCUDAHardwareDevice){
    int delta_t = test_sp_phasing_speed_common(tc,alg,NULL,size,oversampling,SpNoConstraints,iterations,SpEngineCUDA);
    printf("CUDA DIFF MAP %dx%d = %g iterations per second\n",size*oversampling,size*oversampling,(1.0e6*iterations)/delta_t);
  }
  iterations = 10;
  int delta_t = test_sp_phasing_speed_common(tc,alg,NULL,size,oversampling,SpNoConstraints,iterations,SpEngineCPU);
  printf("CPU DIFF MAP %dx%d = %g iterations per second\n",size*oversampling,size*oversampling,(1.0e6*iterations)/delta_t);
#endif
  PRINT_DONE;
}


void test_sp_support_hio(CuTest * tc){
  /* Simple phasing example */
  int size = 4;
  int oversampling = 2;
  sp_smap * beta = sp_smap_create_from_pair(0,0.8);
  real tol = 1e-5;
  SpPhasingAlgorithm * alg = sp_phasing_hio_alloc(beta,0);
  sp_smap * blur_radius = sp_smap_alloc(2);
  sp_smap_insert(blur_radius,0,3);
  sp_smap_insert(blur_radius,2000,0.7);
  sp_smap * threshold = sp_smap_alloc(1);
  sp_smap_insert(threshold,0,0.15);  
  sp_smap * area = sp_smap_alloc(1);
  sp_smap_insert(area,0,1.3*1.0/(oversampling*oversampling));  
  SpSupportArray * sup_alg = sp_support_array_init(sp_support_threshold_alloc(blur_radius,threshold),20);
  CuAssertIntEquals(tc,test_sp_support_common(tc,alg,sup_alg,size,oversampling,SpNoConstraints,0,tol),1);
  sup_alg = sp_support_array_init(sp_support_area_alloc(blur_radius,area),20);
  CuAssertIntEquals(tc,test_sp_support_common(tc,alg,sup_alg,size,oversampling,SpNoConstraints,0,tol),1);
  PRINT_DONE;
}

void test_sp_support_raar(CuTest * tc){
  /* Simple phasing example */
  int size = 4;
  int oversampling = 2;
  sp_smap * beta = sp_smap_create_from_pair(0,0.8);
  real tol = 1e-5;
  SpPhasingAlgorithm * alg = sp_phasing_raar_alloc(beta,0);
  sp_smap * blur_radius = sp_smap_alloc(2);
  sp_smap_insert(blur_radius,0,3);
  sp_smap_insert(blur_radius,2000,0.7);
  sp_smap * threshold = sp_smap_alloc(1);
  sp_smap_insert(threshold,0,0.15);  
  sp_smap * area = sp_smap_alloc(1);
  sp_smap_insert(area,0,1.3*1.0/(oversampling*oversampling));  
  SpSupportArray * sup_alg = sp_support_array_init(sp_support_threshold_alloc(blur_radius,threshold),20);
  CuAssertIntEquals(tc,test_sp_support_common(tc,alg,sup_alg,size,oversampling,SpNoConstraints,0,tol),1);
  sup_alg = sp_support_array_init(sp_support_area_alloc(blur_radius,area),20);
  CuAssertIntEquals(tc,test_sp_support_common(tc,alg,sup_alg,size,oversampling,SpNoConstraints,0,tol),1);
  PRINT_DONE;
}

void test_sp_support_diff_map(CuTest * tc){
  /* Simple phasing example */
  int size = 4;
  int oversampling = 2;
  sp_smap * beta = sp_smap_create_from_pair(0,0.8);
  real tol = 1e-5;
  SpPhasingAlgorithm * alg = sp_phasing_diff_map_alloc(beta,INFINITY,INFINITY,0);
  sp_smap * blur_radius = sp_smap_alloc(2);
  sp_smap_insert(blur_radius,0,3);
  sp_smap_insert(blur_radius,2000,0.7);
  sp_smap * threshold = sp_smap_alloc(1);
  sp_smap_insert(threshold,0,0.15);  
  sp_smap * area = sp_smap_alloc(1);
  sp_smap_insert(area,0,1.3*1.0/(oversampling*oversampling));  
  SpSupportArray * sup_alg = sp_support_array_init(sp_support_threshold_alloc(blur_radius,threshold),20);
  CuAssertIntEquals(tc,test_sp_support_common(tc,alg,sup_alg,size,oversampling,SpNoConstraints,0,tol),1);
  sup_alg = sp_support_array_init(sp_support_area_alloc(blur_radius,area),20);
  CuAssertIntEquals(tc,test_sp_support_common(tc,alg,sup_alg,size,oversampling,SpNoConstraints,0,tol),1);
  PRINT_DONE;
}


void test_sp_support_cuda(CuTest * tc){
  if(sp_cuda_get_device_type() == SpCUDAHardwareDevice){
    real oversampling = 2;
    sp_smap * blur_radius = sp_smap_alloc(2);
    sp_smap_insert(blur_radius,0,3);
    sp_smap_insert(blur_radius,2000,0.7);
    sp_smap * threshold = sp_smap_alloc(1);
    sp_smap_insert(threshold,0,0.15);  
    sp_smap * area = sp_smap_alloc(1);
    sp_smap_insert(area,0,1.3*1.0/(oversampling*oversampling));  
    SpSupportArray * sup_alg = sp_support_array_init(sp_support_threshold_alloc(blur_radius,threshold),4);
    test_sp_support_cuda_common(tc,sup_alg);
    sup_alg = sp_support_array_init(sp_support_area_alloc(blur_radius,area),4);
    test_sp_support_cuda_common(tc,sup_alg);
    PRINT_DONE;
  }
}


void test_sp_phasing_fourier_constraints(CuTest * tc){
  /* Simple phasing example */
  int size = 4;
  int oversampling = 3;
  sp_smap * beta = sp_smap_create_from_pair(0,0.8);
  real stop_tol = 1e-10;
  real match_tol = 1e-4;
  int nruns = 30;
  Image * pada;
  SpPhasingAlgorithm * alg;
  int n_success = 0;
  alg = sp_phasing_hio_alloc(beta,SpCentrosymmetricObject);
  n_success = 0;
  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,SpCentrosymmetricObject);
    sp_image_write(pada,"out.png",SpColormapGrayScale);
    n_success += test_sp_phasing_success_common(tc,alg,pada,0,stop_tol,match_tol);  
    sp_image_free(pada);
  }
#ifndef NDEBUG
  printf("Centrosymmetric HIO object success rate = %5.4f\n",(real)n_success/nruns);
#endif
  /*
  alg = sp_phasing_hio_alloc(beta,0);
  n_success = 0;
  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,SpCentrosymmetricObject);
    sp_image_write(pada,"out.png",SpColormapGrayScale);
    n_success += test_sp_phasing_success_common(tc,alg,pada,0,stop_tol,match_tol);  
    sp_image_free(pada);
  }
  //#ifndef NDEBUG
  printf("Complex HIO object success rate = %5.4f\n",(real)n_success/nruns);
  //#endif
  */
}

CuSuite* phasing_get_suite(void)
{
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, test_sp_phasing_hio);
  SUITE_ADD_TEST(suite,test_sp_phasing_hio_speed_by_size);
  SUITE_ADD_TEST(suite, test_sp_support_cuda);
  SUITE_ADD_TEST(suite, test_sp_support_speed);
  SUITE_ADD_TEST(suite, test_sp_phasing_hio_speed);
  SUITE_ADD_TEST(suite, test_sp_support_hio);
  SUITE_ADD_TEST(suite, test_sp_phasing_hio_success_rate);
  SUITE_ADD_TEST(suite, test_sp_phasing_hio_noisy_success_rate);
  SUITE_ADD_TEST(suite, test_sp_phasing_raar);
  SUITE_ADD_TEST(suite, test_sp_phasing_raar_speed);
  SUITE_ADD_TEST(suite, test_sp_support_raar);
  SUITE_ADD_TEST(suite, test_sp_phasing_raar_success_rate);
  SUITE_ADD_TEST(suite, test_sp_phasing_raar_noisy_success_rate);
  SUITE_ADD_TEST(suite, test_sp_phasing_diff_map_speed); 
  SUITE_ADD_TEST(suite, test_sp_support_diff_map);
  SUITE_ADD_TEST(suite, test_sp_phasing_diff_map_success_rate);
  SUITE_ADD_TEST(suite, test_sp_phasing_diff_map_noisy_success_rate);  
  SUITE_ADD_TEST(suite, test_sp_phasing_fourier_constraints);

  return suite;
}
