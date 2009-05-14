
#include "AllTests.h"


static Image * create_test_image(int size, real oversampling,SpPhasingConstraints c);

int test_sp_phasing_common(CuTest * tc,SpPhasingAlgorithm * alg,int size, real oversampling, SpPhasingConstraints solution_constraints,
			   real beamstop,real tol){
  int attempts = 0;
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
    CuAssertTrue(tc,sp_phaser_init(ph,alg,f,support) == 0);
    int i =0;
    double change = 0;
    int max_iter = 300;
    CuAssertTrue(tc,sp_phaser_init_model(ph,NULL,SpModelRandomPhases) == 0); 
    do{
      CuAssertTrue(tc,sp_phaser_iterate(ph) == 0);
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
  CuAssertTrue(tc,sp_phaser_init(ph,alg,f,support) == 0);

  int i =0;
  double change = 0;
  int max_iter = 300;
  CuAssertTrue(tc,sp_phaser_init_model(ph,NULL,SpModelRandomPhases) == 0); 
  do{
    CuAssertTrue(tc,sp_phaser_iterate(ph) == 0);
    change = sp_image_integrate2(sp_phaser_model_change(ph));
    //      printf("Iter = %d Delta = %g\n",i,change);
    i++;
  }while(change > stop_tol && i < max_iter);
  Image * model = sp_phaser_model(ph);
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

static Image * create_test_image(int size, real oversampling,SpPhasingConstraints c){
  Image * a = sp_image_alloc(size,size,1);
  for(int i = 0;i<sp_image_size(a);i++){
    if(c & SpRealObject){
      a->image->data[i] = sp_cinit(p_drand48()-0.5,0);
    }else if(c & SpPositiveRealObject){
      a->image->data[i] = sp_cinit(p_drand48(),0);
    }else if(c & SpPositiveComplexObject){
      a->image->data[i] = sp_cinit(p_drand48(),p_drand48());
    }else{
      a->image->data[i] = sp_cinit(p_drand48()-0.5,p_drand48()-0.5);
    }
  }
  a->phased = 1;
  Image * pada = sp_image_edge_extend(a,size*(oversampling-1),SP_ZERO_PAD_EDGE,SP_2D);
  pada->phased = 1;
  sp_image_free(a);
  return pada;  
}

void test_sp_phasing_hio_success_rate(CuTest * tc){
  /* Simple phasing example */
  int size = 4;
  int oversampling = 2;
  real beta = 0.8;
  real stop_tol = 1e-10;
  real match_tol = 1e-4;
  int nruns = 10;
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
  printf("HIO Complex success rate = %5.2f\n",(real)n_success/nruns);
  alg = sp_phasing_hio_alloc(beta,SpPositiveComplexObject);
  CuAssertTrue(tc,n_success>0);
  n_success = 0;
  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,SpPositiveComplexObject);  
    n_success += test_sp_phasing_success_common(tc,alg,pada,0,stop_tol,match_tol);  
    sp_image_free(pada);
  }
  printf("HIO Positive Complex success rate = %5.2f\n",(real)n_success/nruns);
  CuAssertTrue(tc,n_success>0);
  alg = sp_phasing_hio_alloc(beta,SpRealObject);
  n_success = 0;
  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,SpRealObject);  
    n_success += test_sp_phasing_success_common(tc,alg,pada,0,stop_tol,match_tol);  
    sp_image_free(pada);
  }
  printf("HIO Real success rate = %5.2f\n",(real)n_success/nruns);
  CuAssertTrue(tc,n_success>0);
  n_success = 0;
  alg = sp_phasing_hio_alloc(beta,SpPositiveRealObject);
  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,SpPositiveRealObject);
    n_success += test_sp_phasing_success_common(tc,alg,pada,0,stop_tol,match_tol);  
    sp_image_free(pada); 
  }
  
  printf("HIO Positive Real success rate = %5.2f\n",(real)n_success/nruns);
  CuAssertTrue(tc,n_success>0);
  n_success = 0;
  alg = sp_phasing_hio_alloc(beta,SpPositiveRealObject| SpPositivityFlipping);
  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,SpPositiveRealObject);
    n_success += test_sp_phasing_success_common(tc,alg,pada,0,stop_tol,match_tol);  
    sp_image_free(pada); 
  }
  printf("HIO Flipping Positive Real success rate = %5.2f\n",(real)n_success/nruns);
  CuAssertTrue(tc,n_success>0);

}


void test_sp_phasing_raar_success_rate(CuTest * tc){
  /* Simple phasing example */
  int size = 4;
  int oversampling = 2;
  real beta = 0.8;
  real stop_tol = 1e-10;
  real match_tol = 1e-4;
  int nruns = 10;
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
  printf("RAAR Complex success rate = %5.2f\n",(real)n_success/nruns);
  alg = sp_phasing_raar_alloc(beta,SpPositiveComplexObject);
  CuAssertTrue(tc,n_success>0);
  n_success = 0;
  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,SpPositiveComplexObject);  
    n_success += test_sp_phasing_success_common(tc,alg,pada,0,stop_tol,match_tol);  
    sp_image_free(pada);
  }
  printf("RAAR Positive Complex success rate = %5.2f\n",(real)n_success/nruns);
  alg = sp_phasing_raar_alloc(beta,SpRealObject);
  CuAssertTrue(tc,n_success>0);
  n_success = 0;
  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,SpRealObject);  
    n_success += test_sp_phasing_success_common(tc,alg,pada,0,stop_tol,match_tol);  
    sp_image_free(pada);
  }
  printf("RAAR Real success rate = %5.2f\n",(real)n_success/nruns);
  CuAssertTrue(tc,n_success>0);

  n_success = 0;
  alg = sp_phasing_raar_alloc(beta,SpPositiveRealObject);
  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,SpPositiveRealObject);
    n_success += test_sp_phasing_success_common(tc,alg,pada,0,stop_tol,match_tol);  
    sp_image_free(pada); 
  }
  printf("RAAR Positive Real success rate = %5.2f\n",(real)n_success/nruns);
  CuAssertTrue(tc,n_success>0);

  n_success = 0;
  alg = sp_phasing_raar_alloc(beta,SpPositiveRealObject| SpPositivityFlipping);
  for(int i = 0;i<nruns;i++){
    pada = create_test_image(size,oversampling,SpPositiveRealObject);
    n_success += test_sp_phasing_success_common(tc,alg,pada,0,stop_tol,match_tol);  
    sp_image_free(pada); 
  }
  printf("RAAR Flipping Positive Real success rate = %5.2f\n",(real)n_success/nruns);
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
  CuAssertIntEquals(tc,test_sp_phasing_common(tc,alg,size,oversampling,SpPositiveComplexObject,2,tol),1);
  alg = sp_phasing_hio_alloc(beta,0);
  CuAssertIntEquals(tc,test_sp_phasing_common(tc,alg,size,oversampling,SpNoConstraints,2,tol),1);
  alg = sp_phasing_hio_alloc(beta,SpPositiveRealObject);
  CuAssertIntEquals(tc,test_sp_phasing_common(tc,alg,size,oversampling,SpPositiveRealObject,2,tol),1);
}



void test_sp_phasing_raar(CuTest * tc){
  /* Simple phasing example */
  int size = 4;
  int oversampling = 2;
  real tol = 1e-12;
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
  CuAssertIntEquals(tc,test_sp_phasing_common(tc,alg,size,oversampling,SpPositiveComplexObject,2,tol),1);
  alg = sp_phasing_raar_alloc(beta,SpPositiveRealObject);
  CuAssertIntEquals(tc,test_sp_phasing_common(tc,alg,size,oversampling,SpPositiveRealObject,2,tol),1);
}

CuSuite* phasing_get_suite(void)
{
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, test_sp_phasing_hio);
  SUITE_ADD_TEST(suite, test_sp_phasing_raar);
  SUITE_ADD_TEST(suite,test_sp_phasing_hio_success_rate);
  SUITE_ADD_TEST(suite,test_sp_phasing_raar_success_rate);
  return suite;
}
