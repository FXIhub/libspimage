#include <assert.h>
#include <setjmp.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "spimage.h"

#include "CuTest.h"
#include <complex.h>

/*-------------------------------------------------------------------------*
 * Helper functions
 *-------------------------------------------------------------------------*/


/*-------------------------------------------------------------------------*
 * CuString Test
 *-------------------------------------------------------------------------*/

static inline void CuAssertComplexEquals(CuTest * tc,complex a, complex b, real delta){
  CuAssertTrue(tc,fabs(creal(a - b)) < delta && fabs(cimag(a - b)) < delta);
}

void test_sp_min(CuTest* tc)
{
  int i;
  CuAssertTrue(tc,sp_min(0,1) == 0);
  CuAssertTrue(tc,sp_min(1,0) == 0);
}

void test_sp_max(CuTest* tc)
{
  int i;
  CuAssertTrue(tc,sp_max(0,1) == 1);
  CuAssertTrue(tc,sp_max(1,0) == 1);
}

void test_sp_vector_alloc(CuTest* tc)
{
  int i;
  sp_vector * v = sp_vector_alloc(4);
  CuAssertTrue(tc,sp_vector_size(v) == 4);
  for(i = 0;i<sp_vector_size(v);i++){
    CuAssertDblEquals(tc,sp_vector_get(v,i),0,REAL_EPSILON);
  }
  sp_vector_free(v);
}

void test_sp_cvector_alloc(CuTest* tc)
{
  int i;
  sp_cvector * v = sp_cvector_alloc(4);
  CuAssertTrue(tc,sp_cvector_size(v) == 4);
  for(i = 0;i<sp_cvector_size(v);i++){
    CuAssertComplexEquals(tc,sp_cvector_get(v,i),0,REAL_EPSILON);
  }
  sp_cvector_free(v);
}


void test_sp_matrix_alloc(CuTest* tc)
{
  int i,j;
  sp_matrix * m = sp_matrix_alloc(4,3);
  CuAssertTrue(tc,sp_matrix_size(m) == 4*3);
  CuAssertTrue(tc,sp_matrix_rows(m) == 4);
  CuAssertTrue(tc,sp_matrix_cols(m) == 3);
  for(i = 0;i<sp_matrix_rows(m);i++){
    for(j = 0;j<sp_matrix_cols(m);j++){
      CuAssertDblEquals(tc,sp_matrix_get(m,i,j),0,REAL_EPSILON);
    }
  }
  sp_matrix_free(m);
}

void test_sp_vector_set_get(CuTest* tc){
  sp_vector * v = sp_vector_alloc(4);
  sp_vector_set(v,1,5);
  CuAssertDblEquals(tc,sp_vector_get(v,1),5,REAL_EPSILON);
  sp_vector_set(v,2,-1);
  CuAssertDblEquals(tc,sp_vector_get(v,2),-1,REAL_EPSILON);
  sp_vector_free(v);  
}

void test_sp_cvector_set_get(CuTest* tc){
  sp_cvector * v = sp_cvector_alloc(4);
  sp_cvector_set(v,1,5-4i);
  CuAssertComplexEquals(tc,sp_cvector_get(v,1),5-4i,REAL_EPSILON);
  sp_cvector_set(v,2,-1i);
  CuAssertDblEquals(tc,sp_cvector_get(v,2),-1i,REAL_EPSILON);
  CuAssertDblEquals(tc,creal(sp_cvector_get(v,2)),0,REAL_EPSILON);
  sp_cvector_free(v);  
}

void test_sp_matrix_set_get(CuTest* tc){
  sp_matrix * m = sp_matrix_alloc(4,3);
  sp_matrix_set(m,1,2,5);
  CuAssertDblEquals(tc,sp_matrix_get(m,1,2),5,REAL_EPSILON);
  sp_matrix_set(m,3,1,-1);
  CuAssertDblEquals(tc,sp_matrix_get(m,3,1),-1,REAL_EPSILON);
  sp_matrix_free(m);  
}

void test_sp_vector_memcpy(CuTest * tc){
  sp_vector * v = sp_vector_alloc(4);
  sp_vector * u = sp_vector_alloc(4);
  int i;
  for(i = 0;i<4;i++){
    sp_vector_set(v,i,rand());
  }
  sp_vector_memcpy(u,v);
  for(i = 0;i<4;i++){
    CuAssertDblEquals(tc,sp_vector_get(v,i),sp_vector_get(u,i),fabs(REAL_EPSILON*sp_vector_get(u,i)));
  }
  sp_vector_free(u);
  sp_vector_free(v);  
}

void test_sp_cvector_memcpy(CuTest * tc){
  sp_cvector * v = sp_cvector_alloc(4);
  sp_cvector * u = sp_cvector_alloc(4);
  int i;
  for(i = 0;i<4;i++){
    sp_cvector_set(v,i,rand()+rand()*1i);
  }
  sp_cvector_memcpy(u,v);
  for(i = 0;i<4;i++){
    CuAssertComplexEquals(tc,sp_cvector_get(v,i),sp_cvector_get(u,i),fabs(REAL_EPSILON*sp_cvector_get(u,i)));
  }
  sp_cvector_free(u);
  sp_cvector_free(v);  
}

void test_sp_matrix_memcpy(CuTest* tc)
{
  int i,j;
  sp_matrix * m = sp_matrix_alloc(4,3);
  sp_matrix * n = sp_matrix_alloc(4,3);
  for(i = 0;i<sp_matrix_rows(m);i++){
    for(j = 0;j<sp_matrix_cols(m);j++){
      sp_matrix_set(m,i,j,rand());
    }
  }
  sp_matrix_memcpy(n,m);
  for(i = 0;i<sp_matrix_rows(m);i++){
    for(j = 0;j<sp_matrix_cols(m);j++){
      CuAssertDblEquals(tc,sp_matrix_get(m,i,j),sp_matrix_get(n,i,j),fabs(REAL_EPSILON*sp_matrix_get(n,i,j)));
    }
  }
  sp_matrix_free(m);
  sp_matrix_free(n);
}


void test_sp_vector_add(CuTest * tc){
  sp_vector * v = sp_vector_alloc(4);
  sp_vector * u = sp_vector_alloc(4);
  sp_vector * t = sp_vector_alloc(4);
  int i;
  for(i = 0;i<4;i++){
    sp_vector_set(v,i,rand());
    sp_vector_set(u,i,rand());
  }
  sp_vector_memcpy(t,v);
  sp_vector_add(t,u);
  for(i = 0;i<4;i++){
    CuAssertDblEquals(tc,sp_vector_get(t,i),sp_vector_get(u,i)+sp_vector_get(v,i),fabs(REAL_EPSILON*sp_vector_get(t,i)));
  }
  sp_vector_free(u);
  sp_vector_free(v);  
  sp_vector_free(t);  
}

void test_sp_cvector_add(CuTest * tc){
  sp_cvector * v = sp_cvector_alloc(4);
  sp_cvector * u = sp_cvector_alloc(4);
  sp_cvector * t = sp_cvector_alloc(4);
  int i;
  for(i = 0;i<4;i++){
    sp_cvector_set(v,i,rand()+rand()*1i);
    sp_cvector_set(u,i,rand()+rand()*1i);
  }
  sp_cvector_memcpy(t,v);
  sp_cvector_add(t,u);
  for(i = 0;i<4;i++){
    CuAssertComplexEquals(tc,sp_cvector_get(t,i),sp_cvector_get(v,i)+sp_cvector_get(u,i),fabs(REAL_EPSILON*sp_cvector_get(t,i)));
  }
  sp_cvector_free(u);
  sp_cvector_free(v);  
  sp_cvector_free(t);  
}

void test_sp_matrix_add(CuTest* tc)
{
  int i,j;
  sp_matrix * m = sp_matrix_alloc(4,3);
  sp_matrix * n = sp_matrix_alloc(4,3);
  sp_matrix * p = sp_matrix_alloc(4,3);
  for(i = 0;i<sp_matrix_rows(m);i++){
    for(j = 0;j<sp_matrix_cols(m);j++){
      sp_matrix_set(m,i,j,rand());
      sp_matrix_set(n,i,j,rand());
    }
  }
  sp_matrix_memcpy(p,m);
  sp_matrix_add(p,n);
  for(i = 0;i<sp_matrix_rows(m);i++){
    for(j = 0;j<sp_matrix_cols(m);j++){
      CuAssertDblEquals(tc,sp_matrix_get(p,i,j),sp_matrix_get(m,i,j) + sp_matrix_get(n,i,j),fabs(REAL_EPSILON*sp_matrix_get(p,i,j)));
    }
  }
  sp_matrix_free(m);
  sp_matrix_free(n);
  sp_matrix_free(p);
}


void test_sp_vector_sub(CuTest * tc){
  sp_vector * v = sp_vector_alloc(4);
  sp_vector * u = sp_vector_alloc(4);
  sp_vector * t = sp_vector_alloc(4);
  int i;
  for(i = 0;i<4;i++){
    sp_vector_set(v,i,rand());
    sp_vector_set(u,i,rand());
  }
  sp_vector_memcpy(t,v);
  sp_vector_sub(t,u);
  for(i = 0;i<4;i++){
    CuAssertDblEquals(tc,sp_vector_get(t,i),sp_vector_get(v,i)-sp_vector_get(u,i),fabs(REAL_EPSILON*sp_vector_get(t,i)));
  }
  sp_vector_free(u);
  sp_vector_free(v);  
  sp_vector_free(t);  
}

void test_sp_cvector_sub(CuTest * tc){
  sp_cvector * v = sp_cvector_alloc(4);
  sp_cvector * u = sp_cvector_alloc(4);
  sp_cvector * t = sp_cvector_alloc(4);
  int i;
  for(i = 0;i<4;i++){
    sp_cvector_set(v,i,rand()+rand()*1i);
    sp_cvector_set(u,i,rand()+rand()*1i);
  }
  sp_cvector_memcpy(t,v);
  sp_cvector_sub(t,u);
  for(i = 0;i<4;i++){
    CuAssertComplexEquals(tc,sp_cvector_get(t,i),sp_cvector_get(v,i)-sp_cvector_get(u,i),fabs(REAL_EPSILON*sp_cvector_get(t,i)));
  }
  sp_cvector_free(u);
  sp_cvector_free(v);  
  sp_cvector_free(t);  
}

void test_sp_matrix_sub(CuTest* tc)
{
  int i,j;
  sp_matrix * m = sp_matrix_alloc(4,3);
  sp_matrix * n = sp_matrix_alloc(4,3);
  sp_matrix * p = sp_matrix_alloc(4,3);
  for(i = 0;i<sp_matrix_rows(m);i++){
    for(j = 0;j<sp_matrix_cols(m);j++){
      sp_matrix_set(m,i,j,rand());
      sp_matrix_set(n,i,j,rand());
    }
  }
  sp_matrix_memcpy(p,m);
  sp_matrix_sub(p,n);
  for(i = 0;i<sp_matrix_rows(m);i++){
    for(j = 0;j<sp_matrix_cols(m);j++){
      CuAssertDblEquals(tc,sp_matrix_get(p,i,j),sp_matrix_get(m,i,j) - sp_matrix_get(n,i,j),fabs(REAL_EPSILON*sp_matrix_get(p,i,j)));
    }
  }
  sp_matrix_free(m);
  sp_matrix_free(n);
  sp_matrix_free(p);
}


void test_sp_vector_mul(CuTest * tc){
  sp_vector * v = sp_vector_alloc(4);
  sp_vector * u = sp_vector_alloc(4);
  sp_vector * t = sp_vector_alloc(4);
  int i;
  for(i = 0;i<4;i++){
    sp_vector_set(v,i,rand());
    sp_vector_set(u,i,rand());
  }
  sp_vector_memcpy(t,v);
  sp_vector_mul(t,u);
  for(i = 0;i<4;i++){
    CuAssertDblEquals(tc,sp_vector_get(t,i),sp_vector_get(v,i)*sp_vector_get(u,i),fabs(REAL_EPSILON*sp_vector_get(t,i)));
  }
  sp_vector_free(u);
  sp_vector_free(v);  
  sp_vector_free(t);  
}

void test_sp_cvector_mul(CuTest * tc){
  sp_cvector * v = sp_cvector_alloc(4);
  sp_cvector * u = sp_cvector_alloc(4);
  sp_cvector * t = sp_cvector_alloc(4);
  int i;
  for(i = 0;i<4;i++){
    sp_cvector_set(v,i,rand()+rand()*1i);
    sp_cvector_set(u,i,rand()+rand()*1i);
  }
  sp_cvector_memcpy(t,v);
  sp_cvector_mul(t,u);
  for(i = 0;i<4;i++){
    CuAssertComplexEquals(tc,sp_cvector_get(t,i),sp_cvector_get(v,i)*sp_cvector_get(u,i),fabs(REAL_EPSILON*sp_cvector_get(t,i)));
  }
  sp_cvector_free(u);
  sp_cvector_free(v);  
  sp_cvector_free(t);  
}

void test_sp_matrix_mul_elements(CuTest* tc)
{
  int i,j;
  sp_matrix * m = sp_matrix_alloc(4,3);
  sp_matrix * n = sp_matrix_alloc(4,3);
  sp_matrix * p = sp_matrix_alloc(4,3);
  for(i = 0;i<sp_matrix_rows(m);i++){
    for(j = 0;j<sp_matrix_cols(m);j++){
      sp_matrix_set(m,i,j,rand());
      sp_matrix_set(n,i,j,rand());
    }
  }
  sp_matrix_memcpy(p,m);
  sp_matrix_mul_elements(p,n);
  for(i = 0;i<sp_matrix_rows(m);i++){
    for(j = 0;j<sp_matrix_cols(m);j++){
      CuAssertDblEquals(tc,sp_matrix_get(p,i,j),sp_matrix_get(m,i,j) * sp_matrix_get(n,i,j),fabs(REAL_EPSILON*sp_matrix_get(p,i,j)));
    }
  }
  sp_matrix_free(m);
  sp_matrix_free(n);
  sp_matrix_free(p);
}

void test_sp_vector_div(CuTest * tc){
  sp_vector * v = sp_vector_alloc(4);
  sp_vector * u = sp_vector_alloc(4);
  sp_vector * t = sp_vector_alloc(4);
  int i;
  for(i = 0;i<4;i++){
    sp_vector_set(v,i,rand());
    sp_vector_set(u,i,rand());
  }
  sp_vector_memcpy(t,v);
  sp_vector_div(t,u);
  for(i = 0;i<4;i++){
    CuAssertDblEquals(tc,sp_vector_get(t,i),sp_vector_get(v,i)/sp_vector_get(u,i),fabs(REAL_EPSILON*sp_vector_get(t,i)));
  }
  sp_vector_free(u);
  sp_vector_free(v);  
  sp_vector_free(t);  
}

void test_sp_cvector_div(CuTest * tc){
  sp_cvector * v = sp_cvector_alloc(4);
  sp_cvector * u = sp_cvector_alloc(4);
  sp_cvector * t = sp_cvector_alloc(4);
  int i;
  for(i = 0;i<4;i++){
    sp_cvector_set(v,i,rand()+rand()*1i);
    sp_cvector_set(u,i,rand()+rand()*1i);
  }
  sp_cvector_memcpy(t,v);
  sp_cvector_div(t,u);
  for(i = 0;i<4;i++){
    CuAssertComplexEquals(tc,sp_cvector_get(t,i),sp_cvector_get(v,i)/sp_cvector_get(u,i),fabs(REAL_EPSILON*sp_cvector_get(t,i)));
  }
  sp_cvector_free(u);
  sp_cvector_free(v);  
  sp_cvector_free(t);  
}

void test_sp_matrix_div_elements(CuTest* tc)
{
  int i,j;
  sp_matrix * m = sp_matrix_alloc(4,3);
  sp_matrix * n = sp_matrix_alloc(4,3);
  sp_matrix * p = sp_matrix_alloc(4,3);
  for(i = 0;i<sp_matrix_rows(m);i++){
    for(j = 0;j<sp_matrix_cols(m);j++){
      sp_matrix_set(m,i,j,rand());
      sp_matrix_set(n,i,j,rand());
    }
  }
  sp_matrix_memcpy(p,m);
  sp_matrix_div_elements(p,n);
  for(i = 0;i<sp_matrix_rows(m);i++){
    for(j = 0;j<sp_matrix_cols(m);j++){
      CuAssertDblEquals(tc,sp_matrix_get(p,i,j),sp_matrix_get(m,i,j) / sp_matrix_get(n,i,j),fabs(REAL_EPSILON*sp_matrix_get(p,i,j)));
    }
  }
  sp_matrix_free(m);
  sp_matrix_free(n);
  sp_matrix_free(p);
}


void test_sp_vector_scale(CuTest * tc){
  sp_vector * v = sp_vector_alloc(4);
  sp_vector * u = sp_vector_alloc(4);
  int i;
  real x = rand()%100;
  for(i = 0;i<4;i++){
    sp_vector_set(v,i,rand());
  }
  sp_vector_memcpy(u,v);
  sp_vector_scale(u,x);
  for(i = 0;i<4;i++){
    CuAssertDblEquals(tc,sp_vector_get(u,i),sp_vector_get(v,i)*x,fabs(REAL_EPSILON*sp_vector_get(u,i)));
  }
  sp_vector_free(u);
  sp_vector_free(v);  
}

void test_sp_cvector_scale(CuTest * tc){
  sp_cvector * v = sp_cvector_alloc(4);
  sp_cvector * u = sp_cvector_alloc(4);
  int i;
  real x = rand()%100;
  for(i = 0;i<4;i++){
    sp_cvector_set(v,i,rand()+rand()*1i);
  }
  sp_cvector_memcpy(u,v);
  sp_cvector_scale(u,x);
  for(i = 0;i<4;i++){
    CuAssertComplexEquals(tc,sp_cvector_get(u,i),sp_cvector_get(v,i)*x,fabs(REAL_EPSILON*sp_cvector_get(u,i)));
  }
  sp_cvector_free(u);
  sp_cvector_free(v);  
}

void test_sp_matrix_scale(CuTest * tc){
  sp_matrix * v = sp_matrix_alloc(4,4);
  sp_matrix * u = sp_matrix_alloc(4,4);
  int i,j;
  real x = rand()%100;
  for(i = 0;i<sp_matrix_rows(v);i++){
    for(j = 0;j<sp_matrix_cols(v);j++){
      sp_matrix_set(v,i,j,rand());
    }
  }
  sp_matrix_memcpy(u,v);
  sp_matrix_scale(u,x);
  for(i = 0;i<sp_matrix_rows(v);i++){
    for(j = 0;j<sp_matrix_cols(v);j++){
      CuAssertDblEquals(tc,sp_matrix_get(u,i,j),sp_matrix_get(v,i,j)*x,fabs(REAL_EPSILON*sp_matrix_get(u,i,j)));
    }
  }
  sp_matrix_free(u);
  sp_matrix_free(v);  
}


void test_sp_vector_add_constant(CuTest * tc){
  sp_vector * v = sp_vector_alloc(4);
  sp_vector * u = sp_vector_alloc(4);
  int i;
  real x = rand()%100;
  for(i = 0;i<4;i++){
    sp_vector_set(v,i,rand());
  }
  sp_vector_memcpy(u,v);
  sp_vector_add_constant(u,x);
  for(i = 0;i<4;i++){
    CuAssertDblEquals(tc,sp_vector_get(u,i),sp_vector_get(v,i)+x,fabs(REAL_EPSILON*sp_vector_get(u,i)));
  }
  sp_vector_free(u);
  sp_vector_free(v);  
}

void test_sp_cvector_add_constant(CuTest * tc){
  sp_cvector * v = sp_cvector_alloc(4);
  sp_cvector * u = sp_cvector_alloc(4);
  int i;
  real x = rand()%100;
  for(i = 0;i<4;i++){
    sp_cvector_set(v,i,rand()+rand()*1i);
  }
  sp_cvector_memcpy(u,v);
  sp_cvector_add_constant(u,x);
  for(i = 0;i<4;i++){
    CuAssertComplexEquals(tc,sp_cvector_get(u,i),sp_cvector_get(v,i)+x,fabs(REAL_EPSILON*sp_cvector_get(u,i)));
  }
  sp_cvector_free(u);
  sp_cvector_free(v);  
}

void test_sp_matrix_add_constant(CuTest * tc){
  sp_matrix * v = sp_matrix_alloc(4,4);
  sp_matrix * u = sp_matrix_alloc(4,4);
  int i,j;
  real x = rand()%100;
  for(i = 0;i<sp_matrix_rows(v);i++){
    for(j = 0;j<sp_matrix_cols(v);j++){
      sp_matrix_set(v,i,j,rand());
    }
  }
  sp_matrix_memcpy(u,v);
  sp_matrix_add_constant(u,x);
  for(i = 0;i<sp_matrix_rows(v);i++){
    for(j = 0;j<sp_matrix_cols(v);j++){
      CuAssertDblEquals(tc,sp_matrix_get(u,i,j),sp_matrix_get(v,i,j)+x,fabs(REAL_EPSILON*sp_matrix_get(u,i,j)));
    }
  }
  sp_matrix_free(u);
  sp_matrix_free(v);  
}

void test_sp_matrix_transpose(CuTest * tc){
  sp_matrix * v = sp_matrix_alloc(4,3);
  sp_matrix * u = sp_matrix_alloc(3,4);
  int i,j;
  real x = rand()%100;
  for(i = 0;i<sp_matrix_rows(v);i++){
    for(j = 0;j<sp_matrix_cols(v);j++){
      sp_matrix_set(v,i,j,rand());
      sp_matrix_set(u,j,i,sp_matrix_get(v,i,j));
    }
  }
  sp_matrix_transpose(v);
  for(i = 0;i<sp_matrix_rows(v);i++){
    for(j = 0;j<sp_matrix_cols(v);j++){
      CuAssertDblEquals(tc,sp_matrix_get(u,i,j),sp_matrix_get(v,i,j),fabs(REAL_EPSILON*sp_matrix_get(u,i,j)));
    }
  }
  sp_matrix_free(u);
  sp_matrix_free(v);  
}

CuSuite* linear_alg_get_suite(void)
{
  CuSuite* suite = CuSuiteNew();
  
  SUITE_ADD_TEST(suite, test_sp_min);
  SUITE_ADD_TEST(suite, test_sp_max);
  
  SUITE_ADD_TEST(suite, test_sp_vector_alloc);
  SUITE_ADD_TEST(suite, test_sp_cvector_alloc);
  SUITE_ADD_TEST(suite, test_sp_matrix_alloc);
  
  SUITE_ADD_TEST(suite, test_sp_vector_set_get);
  SUITE_ADD_TEST(suite, test_sp_cvector_set_get);
  SUITE_ADD_TEST(suite, test_sp_matrix_set_get);
  
  SUITE_ADD_TEST(suite, test_sp_vector_memcpy);
  SUITE_ADD_TEST(suite, test_sp_cvector_memcpy);
  SUITE_ADD_TEST(suite, test_sp_matrix_memcpy);
  
  SUITE_ADD_TEST(suite, test_sp_vector_add);
  SUITE_ADD_TEST(suite, test_sp_cvector_add);
  SUITE_ADD_TEST(suite, test_sp_matrix_add);
  
  SUITE_ADD_TEST(suite, test_sp_vector_sub);
  SUITE_ADD_TEST(suite, test_sp_cvector_sub);
  SUITE_ADD_TEST(suite, test_sp_matrix_sub);
  
  SUITE_ADD_TEST(suite, test_sp_vector_mul);
  SUITE_ADD_TEST(suite, test_sp_cvector_mul);
  SUITE_ADD_TEST(suite, test_sp_matrix_mul_elements);
  
  
  SUITE_ADD_TEST(suite, test_sp_vector_div);
  SUITE_ADD_TEST(suite, test_sp_cvector_div);
  SUITE_ADD_TEST(suite, test_sp_matrix_div_elements);
  
  SUITE_ADD_TEST(suite, test_sp_vector_scale);
  SUITE_ADD_TEST(suite, test_sp_cvector_scale);
  SUITE_ADD_TEST(suite, test_sp_matrix_scale);
  
  SUITE_ADD_TEST(suite, test_sp_vector_add_constant);
  SUITE_ADD_TEST(suite, test_sp_cvector_add_constant);
  SUITE_ADD_TEST(suite, test_sp_matrix_add_constant);

  SUITE_ADD_TEST(suite, test_sp_matrix_transpose);
  return suite;
}



void RunAllTests(void)
{
	CuString *output = CuStringNew();
	CuSuite* suite = CuSuiteNew();

	CuSuiteAddSuite(suite, linear_alg_get_suite());

	CuSuiteRun(suite);
	CuSuiteSummary(suite, output);
	CuSuiteDetails(suite, output);
	printf("%s\n", output->buffer);
}

int main(void)
{
	RunAllTests();
}
