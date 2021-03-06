
#include "AllTests.h"

void test_sp_prtf_basic(CuTest* tc){
  int size = 2;
  Image * list[2];
  Image * a = sp_image_alloc(size,size,size);
  for(int i = 0;i<sp_image_size(a);i++){
    a->image->data[i] = sp_cinit(p_drand48(),p_drand48());
  }
  Image * b = sp_image_duplicate(a,SP_COPY_ALL);
  list[0] = sp_image_fft(a);
  list[1] = sp_image_fft(b);
  Image * prtf = sp_prtf_basic(list,2);
  for(int i = 0;i<sp_image_size(prtf);i++){
    /* Check that the magnitude is 1 */
    CuAssertDblEquals(tc,1,sp_cabs(prtf->image->data[i]),fabs(REAL_EPSILON));
  }
  sp_image_free(list[1]);
  for(int i = 0;i<sp_image_size(b);i++){
    b->image->data[i] = sp_cinit(p_drand48(),p_drand48());
  }
  list[1] = sp_image_fft(b);
  sp_image_free(prtf);
  prtf = sp_prtf_basic(list,2);
  for(int i = 0;i<sp_image_size(prtf);i++){
    /* Check that the magnitude is 1 */
    CuAssertTrue(tc,sp_cabs(prtf->image->data[i]) >= 0);
    CuAssertTrue(tc,sp_cabs(prtf->image->data[i]) <= 1);
  }
  sp_image_free(prtf);
  sp_image_free(a);
  sp_image_free(b);
  sp_image_free(list[0]);
  sp_image_free(list[1]);
}


void test_sp_prtf_advanced(CuTest* tc){
  int size = 10;
  Image * list[2];
  Image * a = sp_image_alloc(size,size,1);
  for(int i = 0;i<sp_image_size(a);i++){
    a->image->data[i] = sp_cinit(p_drand48(),p_drand48());
  }
  Image * tmp = sp_image_shift(a);
  sp_image_free(a);
  a = tmp;
  Image * b = sp_image_duplicate(a,SP_COPY_ALL);
  list[0] = sp_image_fft(a);
  list[1] = sp_image_fft(b);
  Image * prtf = sp_prtf_advanced(list,2,SpFourierSpace|SpOutOfPlace);
  for(int i = 0;i<sp_image_size(prtf);i++){
    /* Check that the magnitude is 1 */
    CuAssertDblEquals(tc,1,sp_cabs(prtf->image->data[i]),1000*fabs(REAL_EPSILON));
  }
  sp_image_free(list[0]);
  sp_image_free(list[1]);
  list[0] = sp_image_duplicate(a,SP_COPY_ALL);
  list[1] = sp_image_duplicate(b,SP_COPY_ALL);
  prtf = sp_prtf_advanced(list,2,SpRealSpace|SpOutOfPlace);
  for(int i = 0;i<sp_image_size(prtf);i++){
    /* Check that the magnitude is 1 */
    CuAssertDblEquals(tc,1,sp_cabs(prtf->image->data[i]),1000*fabs(REAL_EPSILON));
  }

  /* test the phase shift correction */
  real phi = p_drand48()*2*M_PI;
  for(int i = 0;i<sp_image_size(b);i++){
    real phase = sp_carg(a->image->data[i])+phi;
    real mag = sp_cabs(a->image->data[i]);
    b->image->data[i] = sp_cinit(mag*cos(phase),mag*sin(phase));
  }
  sp_image_free(list[1]);
  list[1] = sp_image_duplicate(b,SP_COPY_ALL);
  prtf = sp_prtf_advanced(list,2,SpRealSpace|SpOutOfPlace);
  for(int i = 0;i<sp_image_size(prtf);i++){
    /* Check that the magnitude is 1 */
    CuAssertDblEquals(tc,1,sp_cabs(prtf->image->data[i]),fabs(1000*REAL_EPSILON));
  }


  sp_image_free(prtf);
  /* test centrosymmetry correction */
  sp_image_reflect(list[1],1,SP_ORIGO);
  sp_image_conj(list[1]);

  prtf = sp_prtf_advanced(list,2,SpRealSpace|SpOutOfPlace);
  for(int i = 0;i<sp_image_size(prtf);i++){
    /* Check that the magnitude is 1 */
    CuAssertDblEquals(tc,1,sp_cabs(prtf->image->data[i]),fabs(1000*REAL_EPSILON));
  }

  /* test translation correction */
  sp_image_translate(list[1], rand()%sp_image_x(list[1]),rand()%sp_image_y(list[1]),rand()%sp_image_z(list[1]),SP_TRANSLATE_WRAP_AROUND);
  sp_image_free(prtf);

  prtf = sp_prtf_advanced(list,2,SpRealSpace|SpOutOfPlace);
  for(int i = 0;i<sp_image_size(prtf);i++){
    /* Check that the magnitude is 1 */
    CuAssertDblEquals(tc,1,sp_cabs(prtf->image->data[i]),fabs(1000*REAL_EPSILON));
  }  

  /* test fine translation correction with centrosymmetry and phase shift */
  list[1]->phased = 1;
  sp_image_fourier_translate(list[1], 0.5,0,0);
  sp_image_free(prtf);
  /* test centrosymmetry correction */
  sp_image_reflect(list[1],1,SP_ORIGO);
  sp_image_conj(list[1]);
  sp_image_phase_shift(list[1],0.25,1);


  prtf = sp_prtf_advanced(list,2,SpRealSpace|SpOutOfPlace);
  for(int i = 0;i<sp_image_size(prtf);i++){
    /* Check that the magnitude is 1 */
    CuAssertDblEquals(tc,1,sp_cabs(prtf->image->data[i]),fabs(1000*REAL_EPSILON));
  }  


  /* test fine translation correction */
  list[1]->phased = 1;
  sp_image_fourier_translate(list[1], 0.5,0,0);
  sp_image_free(prtf);

  prtf = sp_prtf_advanced(list,2,SpRealSpace|SpOutOfPlace);
  for(int i = 0;i<sp_image_size(prtf);i++){
    /* Check that the magnitude is 1 */
    CuAssertDblEquals(tc,1,sp_cabs(prtf->image->data[i]),fabs(1000*REAL_EPSILON));
  }  

  /* test fine translation correction with in place prtf*/
  list[1]->phased = 1;
  sp_image_fourier_translate(list[1], 0.5,0,0);
  sp_image_free(prtf);

  prtf = sp_prtf_advanced(list,2,SpRealSpace|SpInPlace);
  for(int i = 0;i<sp_image_size(prtf);i++){
    /* Check that the magnitude is 1 */
    CuAssertDblEquals(tc,1,sp_cabs(prtf->image->data[i]),fabs(1000*REAL_EPSILON));
    CuAssertComplexEquals(tc,list[1]->image->data[i],list[0]->image->data[i],fabs(1000*REAL_EPSILON*(sp_cabs(list[1]->image->data[i]))+REAL_EPSILON*1000));  
  }  

  sp_image_free(list[1]);
  for(int i = 0;i<sp_image_size(b);i++){
    b->image->data[i] = sp_cinit(p_drand48(),p_drand48());
  }
  list[1] = sp_image_duplicate(b,SP_COPY_ALL);
  sp_image_free(prtf);
  prtf = sp_prtf_advanced(list,2,SpRealSpace|SpOutOfPlace);
  for(int i = 0;i<sp_image_size(prtf);i++){
    /* Check that the magnitude is 1 */
    CuAssertTrue(tc,sp_cabs(prtf->image->data[i]) >= 0);
    CuAssertTrue(tc,sp_cabs(prtf->image->data[i]) <= 1);
  }
  sp_image_free(prtf);
  sp_image_free(a);
  sp_image_free(b);
  sp_image_free(list[0]);
  sp_image_free(list[1]);
}


CuSuite* prtf_get_suite(void)
{
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, test_sp_prtf_basic);
  SUITE_ADD_TEST(suite, test_sp_prtf_advanced);
  return suite;
}

