
#include "AllTests.h"


void test_sp_phasing_common(CuTest * tc,SpPhasingAlgorithm * alg,Image * solution, real beamstop,real tol){
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

  while(1){
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
    if(i < max_iter){
      break;
    }
  }
  Image * model = sp_phaser_model(ph);
  sp_image_superimpose_fractional(solution,model,SpEnantiomorph|SpCorrectPhaseShift,1);
  for(int i =0;i<sp_image_size(solution);i++){
    if(sp_real(support->image->data[i])){
      CuAssertComplexEquals(tc,model->image->data[i],solution->image->data[i],1e-2);
    }
  }
  sp_phaser_free(ph);
  sp_image_free(f);
  sp_image_free(support);
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

void test_sp_phasing_hio(CuTest * tc){
  /* Simple phasing example */
  int size = 4;
  int oversampling = 2;
  real beta = 0.8;
  real tol = 1e-10;
  Image * pada = create_test_image(size,oversampling,0);
  SpPhasingAlgorithm * alg = sp_phasing_hio_alloc(beta,0);
  test_sp_phasing_common(tc,alg,pada,0,tol);
  sp_image_free(pada);
  pada = create_test_image(size,oversampling,SpPositiveRealObject);
  alg = sp_phasing_hio_alloc(beta,SpPositiveRealObject);
  test_sp_phasing_common(tc,alg,pada,0,tol);
  sp_image_free(pada);
  pada = create_test_image(size,oversampling,SpRealObject);
  alg = sp_phasing_hio_alloc(beta,SpRealObject);
  test_sp_phasing_common(tc,alg,pada,0,tol);
  sp_image_free(pada);
  pada = create_test_image(size,oversampling,SpPositiveComplexObject);
  alg = sp_phasing_hio_alloc(beta,SpPositiveComplexObject);
  test_sp_phasing_common(tc,alg,pada,0,tol);
  sp_image_free(pada);
  pada = create_test_image(size,oversampling,SpPositiveComplexObject);
  alg = sp_phasing_hio_alloc(beta,SpPositiveComplexObject);
  test_sp_phasing_common(tc,alg,pada,2,tol);
  sp_image_free(pada);
  pada = create_test_image(size,oversampling,0);
  alg = sp_phasing_hio_alloc(beta,0);
  test_sp_phasing_common(tc,alg,pada,2,tol);
  sp_image_free(pada);
  pada = create_test_image(size,oversampling,SpPositiveRealObject);
  alg = sp_phasing_hio_alloc(beta,SpPositiveRealObject);
  test_sp_phasing_common(tc,alg,pada,2,tol);
  sp_image_free(pada);
}



void test_sp_phasing_raar(CuTest * tc){
  /* Simple phasing example */
  int size = 4;
  int oversampling = 2;
  real tol = 1e-12;
  real beta = 0.8;
  Image * a = sp_image_alloc(size,size,1);
  for(int i = 0;i<sp_image_size(a);i++){
    a->image->data[i] = sp_cinit(p_drand48()-0.5,p_drand48()-0.5);
  }
  a->phased = 1;
  Image * pada = sp_image_edge_extend(a,size*(oversampling-1),SP_ZERO_PAD_EDGE,SP_2D);
  pada->phased = 1;
  SpPhasingAlgorithm * alg = sp_phasing_raar_alloc(beta,0);
  test_sp_phasing_common(tc,alg,pada,0,tol);
  sp_image_free(pada);
  sp_image_free(a);
  a = sp_image_alloc(size,size,1);
  for(int i = 0;i<sp_image_size(a);i++){
    a->image->data[i] = sp_cinit(p_drand48(),0);
  }
  a->phased = 1;
  pada = sp_image_edge_extend(a,size*(oversampling-1),SP_ZERO_PAD_EDGE,SP_2D);
  pada->phased = 1;
  alg = sp_phasing_raar_alloc(beta,SpPositiveRealObject);
  test_sp_phasing_common(tc,alg,pada,0,tol);
  sp_image_free(pada);
  sp_image_free(a);
  a = sp_image_alloc(size,size,1);
  for(int i = 0;i<sp_image_size(a);i++){
    a->image->data[i] = sp_cinit(p_drand48()-0.5,0);
  }
  a->phased = 1;
  pada = sp_image_edge_extend(a,size*(oversampling-1),SP_ZERO_PAD_EDGE,SP_2D);
  pada->phased = 1;
  alg = sp_phasing_raar_alloc(beta,SpRealObject);
  test_sp_phasing_common(tc,alg,pada,0,tol);
  sp_image_free(pada);
  sp_image_free(a);
  a = sp_image_alloc(size,size,1);
  for(int i = 0;i<sp_image_size(a);i++){
    a->image->data[i] = sp_cinit(p_drand48(),p_drand48());
  }
  a->phased = 1;
  pada = sp_image_edge_extend(a,size*(oversampling-1),SP_ZERO_PAD_EDGE,SP_2D);
  pada->phased = 1;
  alg = sp_phasing_raar_alloc(beta,SpPositiveComplexObject);
  test_sp_phasing_common(tc,alg,pada,0,tol);
  sp_image_free(pada);
  sp_image_free(a);
  a = sp_image_alloc(size,size,1);
  for(int i = 0;i<sp_image_size(a);i++){
    a->image->data[i] = sp_cinit(p_drand48(),p_drand48());
  }
  a->phased = 1;
  pada = sp_image_edge_extend(a,size*(oversampling-1),SP_ZERO_PAD_EDGE,SP_2D);
  pada->phased = 1;
  alg = sp_phasing_raar_alloc(beta,SpPositiveComplexObject);
  test_sp_phasing_common(tc,alg,pada,2,tol);
  sp_image_free(pada);
  sp_image_free(a);
  a = sp_image_alloc(size,size,1);
  for(int i = 0;i<sp_image_size(a);i++){
    a->image->data[i] = sp_cinit(p_drand48(),0);
  }
  a->phased = 1;
  pada = sp_image_edge_extend(a,size*(oversampling-1),SP_ZERO_PAD_EDGE,SP_2D);
  pada->phased = 1;
  alg = sp_phasing_raar_alloc(beta,SpPositiveRealObject);
  test_sp_phasing_common(tc,alg,pada,2,tol);
  sp_image_free(pada);
  sp_image_free(a);
}

CuSuite* phasing_get_suite(void)
{
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, test_sp_phasing_hio);
  SUITE_ADD_TEST(suite, test_sp_phasing_raar);
  return suite;
}
