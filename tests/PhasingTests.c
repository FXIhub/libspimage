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
    CuAssertTrue(tc,sp_phaser_init(ph,alg,f,support,SpEngineAutomatic) == 0);
    int i =0;
    change = 0;
    int max_iter = 300;
    CuAssertTrue(tc,sp_phaser_init_model(ph,NULL,SpModelRandomPhases) == 0); 
    do{
      CuAssertTrue(tc,sp_phaser_iterate(ph,1,SpOutputModel|SpOutputModelChange) == 0);
      change = sp_image_integrate2(sp_phaser_model_change(ph,0));
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


int test_sp_phasing_speed_common(CuTest * tc,SpPhasingAlgorithm * alg,int size, int oversampling, SpPhasingConstraints solution_constraints,
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
  CuAssertTrue(tc,sp_phaser_init(ph,alg,f,support,engine) == 0);
  CuAssertTrue(tc,sp_phaser_init_model(ph,NULL,SpModelRandomPhases) == 0); 
  int timer = sp_timer_start();
  CuAssertTrue(tc,sp_phaser_iterate(ph,iterations,0) == 0);
  /* retrieve result and make sure the calculations are finished */
  CuAssertTrue(tc,sp_phaser_iterate(ph,0,SpOutputModel) == 0);
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
    CuAssertTrue(tc,sp_phaser_init(ph,alg,f,support,SpEngineAutomatic) == 0);
  //    CuAssertTrue(tc,sp_phaser_init(ph,alg,f,support,SpEngineCPU) == 0);

  int i =0;
  double change = 0;
  int max_iter = 300;
  CuAssertTrue(tc,sp_phaser_init_model(ph,NULL,SpModelRandomPhases) == 0); 
  do{
    CuAssertTrue(tc,sp_phaser_iterate(ph,1,SpOutputModel|SpOutputModelChange) == 0);
    change = sp_image_integrate2(sp_phaser_model_change(ph,0));
    //      printf("Iter = %d Delta = %g\n",i,change);
    i++;
  }while(change > stop_tol && i < max_iter);
  Image * model = sp_phaser_model(ph,0);
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
  CuAssertTrue(tc,sp_phaser_init(ph,alg,f,support,SpEngineAutomatic) == 0);

  int i =0;
  double change = 0;
  int max_iter = 300;
  CuAssertTrue(tc,sp_phaser_init_model(ph,NULL,SpModelRandomPhases) == 0); 
  do{
    CuAssertTrue(tc,sp_phaser_iterate(ph,1,SpOutputModel|SpOutputModelChange) == 0);
    change = sp_image_integrate2(sp_phaser_model_change(ph,0));
    //      printf("Iter = %d Delta = %g\n",i,change);
    i++;
  }while(change > stop_tol && i < max_iter);
  Image * model = sp_phaser_model(ph,0);
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
  a->phased = 1;
  return a;  
}

void test_sp_phasing_hio_success_rate(CuTest * tc){
  /* Simple phasing example */
  int size = 4;
  int oversampling = 2;
  real beta = 0.8;
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

}



void test_sp_phasing_hio_noisy_success_rate(CuTest * tc){
  /* Simple phasing example */
  int size = 4;
  int oversampling = 2;
  real beta = 0.8;
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

}


void test_sp_phasing_raar_success_rate(CuTest * tc){
  /* Simple phasing example */
  int size = 4;
  int oversampling = 2;
  real beta = 0.8;
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
}


void test_sp_phasing_raar_noisy_success_rate(CuTest * tc){
  /* Simple phasing example */
  int size = 4;
  int oversampling = 2;
  real beta = 0.8;
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
}


void test_sp_phasing_hio(CuTest * tc){
  /* Simple phasing example */
  int size = 4;
  int oversampling = 2;
  real beta = 0.8;
  real tol = 1e-10;
  SpPhasingAlgorithm * alg = sp_phasing_hio_alloc(beta,0);
  CuAssertIntEquals(tc,test_sp_phasing_common(tc,alg,size,oversampling,SpNoConstraints,0,tol),1);
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
}


void test_sp_phasing_hio_speed(CuTest * tc){
#ifndef NDEBUG
  /* Simple phasing example */
  int size = 512;
  int oversampling = 2;
  real beta = 0.8;
  SpPhasingAlgorithm * alg = sp_phasing_hio_alloc(beta,0);
  int iterations = 2000;
  if(sp_cuda_get_device_type() == SpCUDAHardwareDevice){
    int delta_t = test_sp_phasing_speed_common(tc,alg,size,oversampling,SpNoConstraints,iterations,SpEngineCUDA);
    printf("CUDA HIO %dx%d = %g iterations per second\n",size*oversampling,size*oversampling,(1.0e6*iterations)/delta_t);
  }
  iterations = 10;
  int delta_t = test_sp_phasing_speed_common(tc,alg,size,oversampling,SpNoConstraints,iterations,SpEngineCPU);
  printf("CPU HIO %dx%d = %g iterations per second\n",size*oversampling,size*oversampling,(1.0e6*iterations)/delta_t);
#endif
}



void test_sp_phasing_raar(CuTest * tc){
  /* Simple phasing example */
  int size = 4;
  int oversampling = 2;
  real tol = 1e-10;
  real beta = 0.8;
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
}

void test_sp_phasing_raar_speed(CuTest * tc){
#ifndef NDEBUG
  /* Simple phasing example */
  int size = 512;
  int oversampling = 2;
  real beta = 0.8;
  SpPhasingAlgorithm * alg = sp_phasing_raar_alloc(beta,0);
  int iterations = 2000;
  if(sp_cuda_get_device_type() == SpCUDAHardwareDevice){
    int delta_t = test_sp_phasing_speed_common(tc,alg,size,oversampling,SpNoConstraints,iterations,SpEngineCUDA);
    printf("CUDA RAAR %dx%d = %g iterations per second\n",size*oversampling,size*oversampling,(1.0e6*iterations)/delta_t);
  }
  iterations = 10;
  int delta_t = test_sp_phasing_speed_common(tc,alg,size,oversampling,SpNoConstraints,iterations,SpEngineCPU);
  printf("CPU RAAR %dx%d = %g iterations per second\n",size*oversampling,size*oversampling,(1.0e6*iterations)/delta_t);
#endif
}


CuSuite* phasing_get_suite(void)
{
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, test_sp_phasing_hio_speed);
  SUITE_ADD_TEST(suite, test_sp_phasing_hio);
  SUITE_ADD_TEST(suite, test_sp_phasing_raar);
  SUITE_ADD_TEST(suite, test_sp_phasing_raar_speed);
  SUITE_ADD_TEST(suite,test_sp_phasing_hio_success_rate);
  SUITE_ADD_TEST(suite,test_sp_phasing_hio_noisy_success_rate);
  SUITE_ADD_TEST(suite,test_sp_phasing_raar_success_rate);
  SUITE_ADD_TEST(suite,test_sp_phasing_raar_noisy_success_rate);
  return suite;
}
