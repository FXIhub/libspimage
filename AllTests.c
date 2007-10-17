#include <assert.h>
#include <setjmp.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "spimage.h"

#include "CuTest.h"
#ifdef _USE_DMALLOC
#include <dmalloc.h>
#endif

/*-------------------------------------------------------------------------*
 * Helper functions
 *-------------------------------------------------------------------------*/


/*-------------------------------------------------------------------------*
 * CuString Test
 *-------------------------------------------------------------------------*/

/* FM: I changed from a function to a macro so that the line number would be reported accurately */

static Complex czero = {0,0};

#define CuAssertComplexEquals(__tc,___a,___b,__delta) do{\
    Complex __a = ___a;\
    Complex __b = ___b;\
    CuAssertTrue(__tc,fabs(sp_real(sp_csub(__a, __b))) < __delta && fabs(sp_imag(sp_csub(__a, __b))) < __delta);\
  }while(0)

void test_sp_min(CuTest* tc)
{
  CuAssertTrue(tc,sp_min(0,1) == 0);
  CuAssertTrue(tc,sp_min(1,0) == 0);
}

void test_sp_max(CuTest* tc)
{
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
    CuAssertComplexEquals(tc,sp_cvector_get(v,i),czero,REAL_EPSILON);
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

void test_sp_3matrix_alloc(CuTest* tc)
{
  int i,j,k;
  sp_3matrix * m = sp_3matrix_alloc(4,3,2);
  CuAssertTrue(tc,sp_3matrix_size(m) == 4*3*2);
  CuAssertTrue(tc,sp_3matrix_x(m) == 4);
  CuAssertTrue(tc,sp_3matrix_y(m) == 3);
  CuAssertTrue(tc,sp_3matrix_z(m) == 2);
  for(i = 0;i<sp_3matrix_x(m);i++){
    for(j = 0;j<sp_3matrix_y(m);j++){
      for(k = 0;k<sp_3matrix_z(m);k++){
	CuAssertDblEquals(tc,sp_3matrix_get(m,i,j,k),0,REAL_EPSILON);
      }
    }
  }
  sp_3matrix_free(m);
}

void test_sp_cmatrix_alloc(CuTest* tc)
{
  int i,j;
  sp_cmatrix * m = sp_cmatrix_alloc(4,3);
  CuAssertTrue(tc,sp_cmatrix_size(m) == 4*3);
  CuAssertTrue(tc,sp_cmatrix_rows(m) == 4);
  CuAssertTrue(tc,sp_cmatrix_cols(m) == 3);
  for(i = 0;i<sp_cmatrix_rows(m);i++){
    for(j = 0;j<sp_cmatrix_cols(m);j++){
      CuAssertComplexEquals(tc,sp_cmatrix_get(m,i,j),czero,REAL_EPSILON);
    }
  }
  sp_cmatrix_free(m);
}

void test_sp_c3matrix_alloc(CuTest* tc)
{
  int i,j,k;
  sp_c3matrix * m = sp_c3matrix_alloc(4,3,2);
  CuAssertTrue(tc,sp_c3matrix_size(m) == 4*3*2);
  CuAssertTrue(tc,sp_c3matrix_x(m) == 4);
  CuAssertTrue(tc,sp_c3matrix_y(m) == 3);
  CuAssertTrue(tc,sp_c3matrix_z(m) == 2);
  for(i = 0;i<sp_c3matrix_x(m);i++){
    for(j = 0;j<sp_c3matrix_y(m);j++){
      for(k = 0;k<sp_c3matrix_z(m);k++){
	CuAssertComplexEquals(tc,sp_c3matrix_get(m,i,j,k),czero,REAL_EPSILON);
      }
    }
  }
  sp_c3matrix_free(m);
}

void test_sp_i3matrix_alloc(CuTest* tc)
{
  int i,j,k;
  sp_i3matrix * m = sp_i3matrix_alloc(4,3,2);
  CuAssertTrue(tc,sp_i3matrix_size(m) == 4*3*2);
  CuAssertTrue(tc,sp_i3matrix_x(m) == 4);
  CuAssertTrue(tc,sp_i3matrix_y(m) == 3);
  CuAssertTrue(tc,sp_i3matrix_z(m) == 2);
  for(i = 0;i<sp_i3matrix_x(m);i++){
    for(j = 0;j<sp_i3matrix_y(m);j++){
      for(k = 0;k<sp_i3matrix_z(m);k++){
	CuAssertIntEquals(tc,sp_i3matrix_get(m,i,j,k),0);
      }
    }
  }
  sp_i3matrix_free(m);
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
  sp_cvector_set(v,1,sp_cinit(5,-4));
  CuAssertComplexEquals(tc,sp_cvector_get(v,1),sp_cinit(5,-4),REAL_EPSILON);
  sp_cvector_set(v,2,sp_cinit(0,-1));
  CuAssertComplexEquals(tc,sp_cvector_get(v,2),sp_cinit(0,-1),REAL_EPSILON);
  CuAssertDblEquals(tc,sp_real(sp_cvector_get(v,2)),0,REAL_EPSILON);
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

void test_sp_3matrix_set_get(CuTest* tc){
  sp_3matrix * m = sp_3matrix_alloc(4,3,2);
  sp_3matrix_set(m,1,2,1,5);
  CuAssertDblEquals(tc,sp_3matrix_get(m,1,2,1),5,REAL_EPSILON);
  sp_3matrix_set(m,3,1,0,-1);
  CuAssertDblEquals(tc,sp_3matrix_get(m,3,1,0),-1,REAL_EPSILON);
  sp_3matrix_free(m);  
}

void test_sp_cmatrix_set_get(CuTest* tc){
  sp_cmatrix * m = sp_cmatrix_alloc(4,3);
  sp_cmatrix_set(m,1,2,sp_cinit(5,0));
  CuAssertComplexEquals(tc,sp_cmatrix_get(m,1,2),sp_cinit(5,0),REAL_EPSILON);
  sp_cmatrix_set(m,3,1,sp_cinit(-1,0));
  CuAssertComplexEquals(tc,sp_cmatrix_get(m,3,1),sp_cinit(-1,0),REAL_EPSILON);
  sp_cmatrix_free(m);  
}

void test_sp_c3matrix_set_get(CuTest* tc){
  sp_c3matrix * m = sp_c3matrix_alloc(4,3,2);
  sp_c3matrix_set(m,1,2,1,sp_cinit(5,0));
  CuAssertComplexEquals(tc,sp_c3matrix_get(m,1,2,1),sp_cinit(5,0),REAL_EPSILON);
  sp_c3matrix_set(m,3,1,0,sp_cinit(-1,0));
  CuAssertComplexEquals(tc,sp_c3matrix_get(m,3,1,0),sp_cinit(-1,0),REAL_EPSILON);
  sp_c3matrix_free(m);  
}

void test_sp_imatrix_set_get(CuTest* tc){
  sp_imatrix * m = sp_imatrix_alloc(4,3);
  sp_imatrix_set(m,1,2,5);
  CuAssertIntEquals(tc,sp_imatrix_get(m,1,2),5);
  sp_imatrix_set(m,3,1,-1);
  CuAssertIntEquals(tc,sp_imatrix_get(m,3,1),-1);
  sp_imatrix_free(m);
}

void test_sp_i3matrix_set_get(CuTest* tc){
  sp_i3matrix * m = sp_i3matrix_alloc(4,3,2);
  sp_i3matrix_set(m,1,2,1,5);
  CuAssertIntEquals(tc,sp_i3matrix_get(m,1,2,1),5);
  sp_i3matrix_set(m,3,1,0,-1);
  CuAssertIntEquals(tc,sp_i3matrix_get(m,3,1,0),-1);
  sp_i3matrix_free(m);  
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
    sp_cvector_set(v,i,sp_cinit(rand(),rand()));
  }
  sp_cvector_memcpy(u,v);
  for(i = 0;i<4;i++){
    CuAssertComplexEquals(tc,sp_cvector_get(v,i),sp_cvector_get(u,i),fabs(REAL_EPSILON*sp_cabs(sp_cvector_get(u,i))));
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

void test_sp_3matrix_memcpy(CuTest* tc)
{
  int i,j,k;
  sp_3matrix * m = sp_3matrix_alloc(4,3,2);
  sp_3matrix * n = sp_3matrix_alloc(4,3,2);
  for(i = 0;i<sp_3matrix_x(m);i++){
    for(j = 0;j<sp_3matrix_y(m);j++){
      for(k = 0;k<sp_3matrix_z(m);k++){
	sp_3matrix_set(m,i,j,k,rand());
      }
    }
  }
  sp_3matrix_memcpy(n,m);
  for(i = 0;i<sp_3matrix_x(m);i++){
    for(j = 0;j<sp_3matrix_y(m);j++){
      for(k = 0;k<sp_3matrix_z(m);k++){
      CuAssertDblEquals(tc,sp_3matrix_get(m,i,j,k),sp_3matrix_get(n,i,j,k),fabs(REAL_EPSILON*sp_3matrix_get(n,i,j,k)));
      }
    }
  }
  sp_3matrix_free(m);
  sp_3matrix_free(n);
}

void test_sp_cmatrix_memcpy(CuTest* tc)
{
  int i,j;
  sp_cmatrix * m = sp_cmatrix_alloc(4,3);
  sp_cmatrix * n = sp_cmatrix_alloc(4,3);
  for(i = 0;i<sp_cmatrix_rows(m);i++){
    for(j = 0;j<sp_cmatrix_cols(m);j++){
      sp_cmatrix_set(m,i,j,sp_cinit(rand(),rand()));
    }
  }
  sp_cmatrix_memcpy(n,m);
  for(i = 0;i<sp_cmatrix_rows(m);i++){
    for(j = 0;j<sp_cmatrix_cols(m);j++){
      CuAssertComplexEquals(tc,sp_cmatrix_get(m,i,j),sp_cmatrix_get(n,i,j),REAL_EPSILON*sp_cabs(sp_cmatrix_get(n,i,j)));
    }
  }
  sp_cmatrix_free(m);
  sp_cmatrix_free(n);
}

void test_sp_c3matrix_memcpy(CuTest* tc)
{
  int i,j,k;
  sp_c3matrix * m = sp_c3matrix_alloc(4,3,2);
  sp_c3matrix * n = sp_c3matrix_alloc(4,3,2);
  for(i = 0;i<sp_c3matrix_x(m);i++){
    for(j = 0;j<sp_c3matrix_y(m);j++){
      for(k = 0;k<sp_c3matrix_z(m);k++){
	sp_c3matrix_set(m,i,j,k,sp_cinit(rand(),rand()));
      }
    }
  }
  sp_c3matrix_memcpy(n,m);
  for(i = 0;i<sp_c3matrix_x(m);i++){
    for(j = 0;j<sp_c3matrix_y(m);j++){
      for(k = 0;k<sp_c3matrix_z(m);k++){
      CuAssertComplexEquals(tc,sp_c3matrix_get(m,i,j,k),sp_c3matrix_get(n,i,j,k),fabs(REAL_EPSILON*sp_cabs(sp_c3matrix_get(n,i,j,k))));
      }
    }
  }
  sp_c3matrix_free(m);
  sp_c3matrix_free(n);
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
    sp_cvector_set(v,i,sp_cinit(rand(),rand()));
    sp_cvector_set(u,i,sp_cinit(rand(),rand()));
  }
  sp_cvector_memcpy(t,v);
  sp_cvector_add(t,u);
  for(i = 0;i<4;i++){
    CuAssertComplexEquals(tc,sp_cvector_get(t,i),sp_cadd(sp_cvector_get(v,i),sp_cvector_get(u,i)),fabs(REAL_EPSILON*sp_cabs(sp_cvector_get(t,i))));
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

void test_sp_3matrix_add(CuTest* tc)
{
  int i,j,k;
  sp_3matrix * m = sp_3matrix_alloc(4,3,2);
  sp_3matrix * n = sp_3matrix_alloc(4,3,2);
  sp_3matrix * p = sp_3matrix_alloc(4,3,2);
  for(i = 0;i<sp_3matrix_x(m);i++){
    for(j = 0;j<sp_3matrix_y(m);j++){
      for(k = 0;k<sp_3matrix_z(m);k++){
	sp_3matrix_set(m,i,j,k,rand());
	sp_3matrix_set(n,i,j,k,rand());
      }
    }
  }

  sp_3matrix_memcpy(p,m);
  sp_3matrix_add(p,n);
  for(i = 0;i<sp_3matrix_x(m);i++){
    for(j = 0;j<sp_3matrix_y(m);j++){
      for(k = 0;k<sp_3matrix_z(m);k++){
	CuAssertDblEquals(tc,sp_3matrix_get(p,i,j,k),sp_3matrix_get(m,i,j,k) + sp_3matrix_get(n,i,j,k),fabs(REAL_EPSILON*sp_3matrix_get(p,i,j,k)));
      }
    }
  }
  sp_3matrix_free(m);
  sp_3matrix_free(n);
  sp_3matrix_free(p);
}

void test_sp_cmatrix_add(CuTest* tc)
{
  int i,j;
  sp_cmatrix * m = sp_cmatrix_alloc(4,3);
  sp_cmatrix * n = sp_cmatrix_alloc(4,3);
  sp_cmatrix * p = sp_cmatrix_alloc(4,3);
  for(i = 0;i<sp_cmatrix_rows(m);i++){
    for(j = 0;j<sp_cmatrix_cols(m);j++){
      sp_cmatrix_set(m,i,j,sp_cinit(rand(),rand()));
      sp_cmatrix_set(n,i,j,sp_cinit(rand(),rand()));
    }
  }
  sp_cmatrix_memcpy(p,m);
  sp_cmatrix_add(p,n,NULL);
  for(i = 0;i<sp_cmatrix_rows(m);i++){
    for(j = 0;j<sp_cmatrix_cols(m);j++){
      CuAssertComplexEquals(tc,sp_cmatrix_get(p,i,j),sp_cadd(sp_cmatrix_get(m,i,j),sp_cmatrix_get(n,i,j)),fabs(REAL_EPSILON*sp_cabs(sp_cmatrix_get(p,i,j))));
    }
  }
  sp_cmatrix_free(m);
  sp_cmatrix_free(n);
  sp_cmatrix_free(p);
}

void test_sp_c3matrix_add(CuTest* tc)
{
  int i,j,k;
  sp_c3matrix * m = sp_c3matrix_alloc(4,3,2);
  sp_c3matrix * n = sp_c3matrix_alloc(4,3,2);
  sp_c3matrix * p = sp_c3matrix_alloc(4,3,2);
  for(i = 0;i<sp_c3matrix_x(m);i++){
    for(j = 0;j<sp_c3matrix_y(m);j++){
      for(k = 0;k<sp_c3matrix_z(m);k++){
	sp_c3matrix_set(m,i,j,k,sp_cinit(rand(),rand()));
	sp_c3matrix_set(n,i,j,k,sp_cinit(rand(),rand()));
      }
    }
  }
  sp_c3matrix_memcpy(p,m);
  sp_c3matrix_add(p,n,NULL);
  for(i = 0;i<sp_c3matrix_x(m);i++){
    for(j = 0;j<sp_c3matrix_y(m);j++){
      for(k = 0;k<sp_c3matrix_z(m);k++){
	CuAssertComplexEquals(tc,sp_c3matrix_get(p,i,j,k),sp_cadd(sp_c3matrix_get(m,i,j,k),sp_c3matrix_get(n,i,j,k)),
			      fabs(REAL_EPSILON*sp_cabs(sp_c3matrix_get(p,i,j,k))));
      }
    }
  }
  sp_c3matrix_free(m);
  sp_c3matrix_free(n);
  sp_c3matrix_free(p);
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
    sp_cvector_set(v,i,sp_cinit(rand(),rand()));
    sp_cvector_set(u,i,sp_cinit(rand(),rand()));
  }
  sp_cvector_memcpy(t,v);
  sp_cvector_sub(t,u);
  for(i = 0;i<4;i++){
    CuAssertComplexEquals(tc,sp_cvector_get(t,i),sp_csub(sp_cvector_get(v,i),sp_cvector_get(u,i)),fabs(REAL_EPSILON*sp_cabs(sp_cvector_get(t,i))));
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

void test_sp_3matrix_sub(CuTest* tc)
{
  int i,j,k;
  sp_3matrix * m = sp_3matrix_alloc(4,3,2);
  sp_3matrix * n = sp_3matrix_alloc(4,3,2);
  sp_3matrix * p = sp_3matrix_alloc(4,3,2);
  for(i = 0;i<sp_3matrix_x(m);i++){
    for(j = 0;j<sp_3matrix_y(m);j++){
      for(k = 0;k<sp_3matrix_z(m);k++){
	sp_3matrix_set(m,i,j,k,rand());
	sp_3matrix_set(n,i,j,k,rand());
      }
    }
  }
  sp_3matrix_memcpy(p,m);
  sp_3matrix_sub(p,n);
  for(i = 0;i<sp_3matrix_x(m);i++){
    for(j = 0;j<sp_3matrix_y(m);j++){
      for(k = 0;k<sp_3matrix_z(m);k++){
	CuAssertDblEquals(tc,sp_3matrix_get(p,i,j,k),sp_3matrix_get(m,i,j,k) - sp_3matrix_get(n,i,j,k),fabs(REAL_EPSILON*sp_3matrix_get(p,i,j,k)));
      }
    }
  }
  sp_3matrix_free(m);
  sp_3matrix_free(n);
  sp_3matrix_free(p);
}

void test_sp_cmatrix_sub(CuTest* tc)
{
  int i,j;
  sp_cmatrix * m = sp_cmatrix_alloc(4,3);
  sp_cmatrix * n = sp_cmatrix_alloc(4,3);
  sp_cmatrix * p = sp_cmatrix_alloc(4,3);
  for(i = 0;i<sp_cmatrix_rows(m);i++){
    for(j = 0;j<sp_cmatrix_cols(m);j++){
      sp_cmatrix_set(m,i,j,sp_cinit(rand(),rand()));
      sp_cmatrix_set(n,i,j,sp_cinit(rand(),rand()));
    }
  }
  sp_cmatrix_memcpy(p,m);
  sp_cmatrix_sub(p,n);
  for(i = 0;i<sp_cmatrix_rows(m);i++){
    for(j = 0;j<sp_cmatrix_cols(m);j++){
      CuAssertComplexEquals(tc,sp_cmatrix_get(p,i,j),sp_csub(sp_cmatrix_get(m,i,j),sp_cmatrix_get(n,i,j)),fabs(REAL_EPSILON*sp_cabs(sp_cmatrix_get(p,i,j))));
    }
  }
  sp_cmatrix_free(m);
  sp_cmatrix_free(n);
  sp_cmatrix_free(p);
}

void test_sp_c3matrix_sub(CuTest* tc)
{
  int i,j,k;
  sp_c3matrix * m = sp_c3matrix_alloc(4,3,2);
  sp_c3matrix * n = sp_c3matrix_alloc(4,3,2);
  sp_c3matrix * p = sp_c3matrix_alloc(4,3,2);
  for(i = 0;i<sp_c3matrix_x(m);i++){
    for(j = 0;j<sp_c3matrix_y(m);j++){
      for(k = 0;k<sp_c3matrix_z(m);k++){
	sp_c3matrix_set(m,i,j,k,sp_cinit(rand(),rand()));
	sp_c3matrix_set(n,i,j,k,sp_cinit(rand(),rand()));
      }
    }
  }
  sp_c3matrix_memcpy(p,m);
  sp_c3matrix_sub(p,n);
  for(i = 0;i<sp_c3matrix_x(m);i++){
    for(j = 0;j<sp_c3matrix_y(m);j++){
      for(k = 0;k<sp_c3matrix_z(m);k++){
	CuAssertComplexEquals(tc,sp_c3matrix_get(p,i,j,k),sp_csub(sp_c3matrix_get(m,i,j,k),sp_c3matrix_get(n,i,j,k)),fabs(REAL_EPSILON*sp_cabs(sp_c3matrix_get(p,i,j,k))));
      }
    }
  }
  sp_c3matrix_free(m);
  sp_c3matrix_free(n);
  sp_c3matrix_free(p);
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
    sp_cvector_set(v,i,sp_cinit(rand(),rand()));
    sp_cvector_set(u,i,sp_cinit(rand(),rand()));
  }
  sp_cvector_memcpy(t,v);
  sp_cvector_mul(t,u);
  for(i = 0;i<4;i++){
    CuAssertComplexEquals(tc,sp_cvector_get(t,i),sp_cmul(sp_cvector_get(v,i),sp_cvector_get(u,i)),fabs(REAL_EPSILON*sp_cabs(sp_cvector_get(t,i))));
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

void test_sp_3matrix_mul_elements(CuTest* tc)
{
  int i,j,k;
  sp_3matrix * m = sp_3matrix_alloc(4,3,2);
  sp_3matrix * n = sp_3matrix_alloc(4,3,2);
  sp_3matrix * p = sp_3matrix_alloc(4,3,2);
  for(i = 0;i<sp_3matrix_x(m);i++){
    for(j = 0;j<sp_3matrix_y(m);j++){
      for(k = 0;k<sp_3matrix_z(m);k++){
	sp_3matrix_set(m,i,j,k,rand());
	sp_3matrix_set(n,i,j,k,rand());
      }
    }
  }
  sp_3matrix_memcpy(p,m);
  sp_3matrix_mul_elements(p,n);
  for(i = 0;i<sp_3matrix_x(m);i++){
    for(j = 0;j<sp_3matrix_y(m);j++){
      for(k = 0;k<sp_3matrix_z(m);k++){
	CuAssertDblEquals(tc,sp_3matrix_get(p,i,j,k),sp_3matrix_get(m,i,j,k) * sp_3matrix_get(n,i,j,k),fabs(REAL_EPSILON*sp_3matrix_get(p,i,j,k)));
      }
    }
  }
  sp_3matrix_free(m);
  sp_3matrix_free(n);
  sp_3matrix_free(p);
}

void test_sp_cmatrix_mul_elements(CuTest* tc)
{
  int i,j;
  sp_cmatrix * m = sp_cmatrix_alloc(4,3);
  sp_cmatrix * n = sp_cmatrix_alloc(4,3);
  sp_cmatrix * p = sp_cmatrix_alloc(4,3);
  for(i = 0;i<sp_cmatrix_rows(m);i++){
    for(j = 0;j<sp_cmatrix_cols(m);j++){
      sp_cmatrix_set(m,i,j,sp_cinit(rand(),rand()));
      sp_cmatrix_set(n,i,j,sp_cinit(rand(),rand()));
    }
  }
  sp_cmatrix_memcpy(p,m);
  sp_cmatrix_mul_elements(p,n);
  for(i = 0;i<sp_cmatrix_rows(m);i++){
    for(j = 0;j<sp_cmatrix_cols(m);j++){
      CuAssertComplexEquals(tc,sp_cmatrix_get(p,i,j),sp_cmul(sp_cmatrix_get(m,i,j),sp_cmatrix_get(n,i,j)),fabs(REAL_EPSILON*sp_cabs(sp_cmatrix_get(p,i,j))));
    }
  }
  sp_cmatrix_free(m);
  sp_cmatrix_free(n);
  sp_cmatrix_free(p);
}

void test_sp_c3matrix_mul_elements(CuTest* tc)
{
  int i,j,k;
  sp_c3matrix * m = sp_c3matrix_alloc(4,3,2);
  sp_c3matrix * n = sp_c3matrix_alloc(4,3,2);
  sp_c3matrix * p = sp_c3matrix_alloc(4,3,2);
  for(i = 0;i<sp_c3matrix_x(m);i++){
    for(j = 0;j<sp_c3matrix_y(m);j++){
      for(k = 0;k<sp_c3matrix_z(m);k++){
	sp_c3matrix_set(m,i,j,k,sp_cinit(rand(),rand()));
	sp_c3matrix_set(n,i,j,k,sp_cinit(rand(),rand()));
      }
    }
  }
  sp_c3matrix_memcpy(p,m);
  sp_c3matrix_mul_elements(p,n);
  for(i = 0;i<sp_c3matrix_x(m);i++){
    for(j = 0;j<sp_c3matrix_y(m);j++){
      for(k = 0;k<sp_c3matrix_z(m);k++){
	CuAssertComplexEquals(tc,sp_c3matrix_get(p,i,j,k),sp_cmul(sp_c3matrix_get(m,i,j,k),sp_c3matrix_get(n,i,j,k)),fabs(REAL_EPSILON*sp_cabs(sp_c3matrix_get(p,i,j,k))));
      }
    }
  }
  sp_c3matrix_free(m);
  sp_c3matrix_free(n);
  sp_c3matrix_free(p);
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
    sp_cvector_set(v,i,sp_cinit(rand(),rand()));
    sp_cvector_set(u,i,sp_cinit(rand(),rand()));
  }
  sp_cvector_memcpy(t,v);
  sp_cvector_div(t,u);
  for(i = 0;i<4;i++){
    CuAssertComplexEquals(tc,sp_cvector_get(t,i),sp_cdiv(sp_cvector_get(v,i),sp_cvector_get(u,i)),fabs(REAL_EPSILON*sp_cabs(sp_cvector_get(t,i))));
  }
  sp_cvector_free(u);
  sp_cvector_free(v);  
  sp_cvector_free(t);  
}

void test_sp_matrix_div_elements(CuTest* tc){
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

void test_sp_3matrix_div_elements(CuTest* tc){
  int i,j,k;
  sp_3matrix * m = sp_3matrix_alloc(4,3,2);
  sp_3matrix * n = sp_3matrix_alloc(4,3,2);
  sp_3matrix * p = sp_3matrix_alloc(4,3,2);
  for(i = 0;i<sp_3matrix_x(m);i++){
    for(j = 0;j<sp_3matrix_y(m);j++){
      for(k = 0;k<sp_3matrix_z(m);k++){
      sp_3matrix_set(m,i,j,k,rand());
      sp_3matrix_set(n,i,j,k,rand());
      }
    }
  }
  sp_3matrix_memcpy(p,m);
  sp_3matrix_div_elements(p,n);
  for(i = 0;i<sp_3matrix_x(m);i++){
    for(j = 0;j<sp_3matrix_y(m);j++){
      for(k = 0;k<sp_3matrix_z(m);k++){
      CuAssertDblEquals(tc,sp_3matrix_get(p,i,j,k),sp_3matrix_get(m,i,j,k) / sp_3matrix_get(n,i,j,k),fabs(REAL_EPSILON*sp_3matrix_get(p,i,j,k)));
      }
    }
  }
  sp_3matrix_free(m);
  sp_3matrix_free(n);
  sp_3matrix_free(p);
}

void test_sp_cmatrix_div_elements(CuTest* tc){
  int i,j;
  sp_cmatrix * m = sp_cmatrix_alloc(4,3);
  sp_cmatrix * n = sp_cmatrix_alloc(4,3);
  sp_cmatrix * p = sp_cmatrix_alloc(4,3);
  for(i = 0;i<sp_cmatrix_rows(m);i++){
    for(j = 0;j<sp_cmatrix_cols(m);j++){
      sp_cmatrix_set(m,i,j,sp_cinit(rand(),rand()));
      sp_cmatrix_set(n,i,j,sp_cinit(rand(),rand()));
    }
  }
  sp_cmatrix_memcpy(p,m);
  sp_cmatrix_div_elements(p,n);
  for(i = 0;i<sp_cmatrix_rows(m);i++){
    for(j = 0;j<sp_cmatrix_cols(m);j++){
      CuAssertComplexEquals(tc,sp_cmatrix_get(p,i,j),sp_cdiv(sp_cmatrix_get(m,i,j),sp_cmatrix_get(n,i,j)),fabs(REAL_EPSILON*sp_cabs(sp_cmatrix_get(p,i,j))));
    }
  }
  sp_cmatrix_free(m);
  sp_cmatrix_free(n);
  sp_cmatrix_free(p);
}

void test_sp_c3matrix_div_elements(CuTest* tc){
  int i,j,k;
  sp_c3matrix * m = sp_c3matrix_alloc(4,3,2);
  sp_c3matrix * n = sp_c3matrix_alloc(4,3,2);
  sp_c3matrix * p = sp_c3matrix_alloc(4,3,2);
  for(i = 0;i<sp_c3matrix_x(m);i++){
    for(j = 0;j<sp_c3matrix_y(m);j++){
      for(k = 0;k<sp_c3matrix_z(m);k++){
	sp_c3matrix_set(m,i,j,k,sp_cinit(rand(),rand()));
	sp_c3matrix_set(n,i,j,k,sp_cinit(rand(),rand()));
      }
    }
  }
  sp_c3matrix_memcpy(p,m);
  sp_c3matrix_div_elements(p,n);
  for(i = 0;i<sp_c3matrix_x(m);i++){
    for(j = 0;j<sp_c3matrix_y(m);j++){
      for(k = 0;k<sp_c3matrix_z(m);k++){
      CuAssertComplexEquals(tc,sp_c3matrix_get(p,i,j,k),sp_cdiv(sp_c3matrix_get(m,i,j,k),sp_c3matrix_get(n,i,j,k)),fabs(REAL_EPSILON*sp_cabs(sp_c3matrix_get(p,i,j,k))));
      }
    }
  }
  sp_c3matrix_free(m);
  sp_c3matrix_free(n);
  sp_c3matrix_free(p);
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
    sp_cvector_set(v,i,sp_cinit(rand()%10000,rand()%10000));
  }
  sp_cvector_memcpy(u,v);
  sp_cvector_scale(u,sp_cinit(x,0));
  for(i = 0;i<4;i++){
    CuAssertComplexEquals(tc,sp_cvector_get(u,i),sp_cscale(sp_cvector_get(v,i),x),fabs(REAL_EPSILON*sp_cabs(sp_cvector_get(u,i))));
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

void test_sp_3matrix_scale(CuTest * tc){
  sp_3matrix * v = sp_3matrix_alloc(4,4,4);
  sp_3matrix * u = sp_3matrix_alloc(4,4,4);
  int i,j,k;
  real x = rand()%100;
  for(i = 0;i<sp_3matrix_x(v);i++){
    for(j = 0;j<sp_3matrix_y(v);j++){
      for(k = 0;k<sp_3matrix_z(v);k++){
	sp_3matrix_set(v,i,j,k,rand());
      }
    }
  }
  sp_3matrix_memcpy(u,v);
  sp_3matrix_scale(u,x);
  for(i = 0;i<sp_3matrix_x(v);i++){
    for(j = 0;j<sp_3matrix_y(v);j++){
      for(k = 0;k<sp_3matrix_z(v);k++){
      CuAssertDblEquals(tc,sp_3matrix_get(u,i,j,k),sp_3matrix_get(v,i,j,k)*x,fabs(REAL_EPSILON*sp_3matrix_get(u,i,j,k)));
      }
    }
  }
  sp_3matrix_free(u);
  sp_3matrix_free(v);  
}

void test_sp_cmatrix_scale(CuTest * tc){
  sp_cmatrix * v = sp_cmatrix_alloc(4,4);
  sp_cmatrix * u = sp_cmatrix_alloc(4,4);
  int i,j;
  real x = rand()%100;
  for(i = 0;i<sp_cmatrix_rows(v);i++){
    for(j = 0;j<sp_cmatrix_cols(v);j++){
      sp_cmatrix_set(v,i,j,sp_cinit(rand(),rand()));
    }
  }
  sp_cmatrix_memcpy(u,v);
  sp_cmatrix_scale(u,sp_cinit(x,0));
  for(i = 0;i<sp_cmatrix_rows(v);i++){
    for(j = 0;j<sp_cmatrix_cols(v);j++){
      CuAssertComplexEquals(tc,sp_cmatrix_get(u,i,j),sp_cscale(sp_cmatrix_get(v,i,j),x),fabs(REAL_EPSILON*sp_cabs(sp_cmatrix_get(u,i,j))));
    }
  }
  sp_cmatrix_free(u);
  sp_cmatrix_free(v);  
}

void test_sp_c3matrix_scale(CuTest * tc){
  sp_c3matrix * v = sp_c3matrix_alloc(4,4,4);
  sp_c3matrix * u = sp_c3matrix_alloc(4,4,4);
  int i,j,k;
  real x = rand()%100;
  for(i = 0;i<sp_c3matrix_x(v);i++){
    for(j = 0;j<sp_c3matrix_y(v);j++){
      for(k = 0;k<sp_c3matrix_z(v);k++){
	sp_c3matrix_set(v,i,j,k,sp_cinit(rand(),rand()));
      }
    }
  }
  sp_c3matrix_memcpy(u,v);
  sp_c3matrix_scale(u,sp_cinit(x,0));
  for(i = 0;i<sp_c3matrix_x(v);i++){
    for(j = 0;j<sp_c3matrix_y(v);j++){
      for(k = 0;k<sp_c3matrix_z(v);k++){
      CuAssertComplexEquals(tc,sp_c3matrix_get(u,i,j,k),sp_cscale(sp_c3matrix_get(v,i,j,k),x),fabs(REAL_EPSILON*sp_cabs(sp_c3matrix_get(u,i,j,k))));
      }
    }
  }
  sp_c3matrix_free(u);
  sp_c3matrix_free(v);  
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
    sp_cvector_set(v,i,sp_cinit(rand(),rand()));
  }
  sp_cvector_memcpy(u,v);
  sp_cvector_add_constant(u,sp_cinit(x,0));
  for(i = 0;i<4;i++){
    CuAssertComplexEquals(tc,sp_cvector_get(u,i),sp_cinit(sp_real(sp_cvector_get(v,i))+x,sp_imag(sp_cvector_get(v,i))),fabs(REAL_EPSILON*sp_cabs(sp_cvector_get(u,i))));
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

void test_sp_3matrix_add_constant(CuTest * tc){
  sp_3matrix * v = sp_3matrix_alloc(4,4,4);
  sp_3matrix * u = sp_3matrix_alloc(4,4,4);
  int i,j,k;
  real x = rand()%100;
  for(i = 0;i<sp_3matrix_x(v);i++){
    for(j = 0;j<sp_3matrix_y(v);j++){
      for(k = 0;k<sp_3matrix_z(v);k++){
	sp_3matrix_set(v,i,j,k,rand());
      }
    }
  }
  sp_3matrix_memcpy(u,v);
  sp_3matrix_add_constant(u,x);
  for(i = 0;i<sp_3matrix_x(v);i++){
    for(j = 0;j<sp_3matrix_y(v);j++){
      for(k = 0;k<sp_3matrix_z(v);k++){
	CuAssertDblEquals(tc,sp_3matrix_get(u,i,j,k),sp_3matrix_get(v,i,j,k)+x,fabs(REAL_EPSILON*sp_3matrix_get(u,i,j,k)));
      }
    }
  }
  sp_3matrix_free(u);
  sp_3matrix_free(v);  
}

void test_sp_cmatrix_add_constant(CuTest * tc){
  sp_cmatrix * v = sp_cmatrix_alloc(4,4);
  sp_cmatrix * u = sp_cmatrix_alloc(4,4);
  int i,j;
  Complex x = {rand()%100, (rand()%100)};
  for(i = 0;i<sp_cmatrix_rows(v);i++){
    for(j = 0;j<sp_cmatrix_cols(v);j++){
      sp_cmatrix_set(v,i,j,sp_cinit(rand(),rand()));
    }
  }
  sp_cmatrix_memcpy(u,v);
  sp_cmatrix_add_constant(u,x);
  for(i = 0;i<sp_cmatrix_rows(v);i++){
    for(j = 0;j<sp_cmatrix_cols(v);j++){
      CuAssertComplexEquals(tc,sp_cmatrix_get(u,i,j),sp_cadd(sp_cmatrix_get(v,i,j),x),fabs(REAL_EPSILON*sp_cabs(sp_cmatrix_get(u,i,j))));
    }
  }
  sp_cmatrix_free(u);
  sp_cmatrix_free(v);  
}

void test_sp_c3matrix_add_constant(CuTest * tc){
  sp_c3matrix * v = sp_c3matrix_alloc(4,4,4);
  sp_c3matrix * u = sp_c3matrix_alloc(4,4,4);
  int i,j,k;
  real x = rand()%100;
  for(i = 0;i<sp_c3matrix_x(v);i++){
    for(j = 0;j<sp_c3matrix_y(v);j++){
      for(k = 0;k<sp_c3matrix_z(v);k++){
	sp_c3matrix_set(v,i,j,k,sp_cinit(rand(),rand()));
      }
    }
  }
  sp_c3matrix_memcpy(u,v);
  sp_c3matrix_add_constant(u,sp_cinit(x,0));
  for(i = 0;i<sp_c3matrix_x(v);i++){
    for(j = 0;j<sp_c3matrix_y(v);j++){
      for(k = 0;k<sp_c3matrix_z(v);k++){
	CuAssertComplexEquals(tc,sp_c3matrix_get(u,i,j,k),sp_cadd(sp_c3matrix_get(v,i,j,k),sp_cinit(x,0)),fabs(REAL_EPSILON*sp_cabs(sp_c3matrix_get(u,i,j,k))));
      }
    }
  }
  sp_c3matrix_free(u);
  sp_c3matrix_free(v);  
}

void test_sp_matrix_transpose(CuTest * tc){
  sp_matrix * v = sp_matrix_alloc(4,3);
  sp_matrix * u = sp_matrix_alloc(3,4);
  int i,j;
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

void test_sp_cmatrix_transpose(CuTest * tc){
  sp_cmatrix * v = sp_cmatrix_alloc(4,3);
  sp_cmatrix * u = sp_cmatrix_alloc(3,4);
  int i,j;
  for(i = 0;i<sp_cmatrix_rows(v);i++){
    for(j = 0;j<sp_cmatrix_cols(v);j++){
      sp_cmatrix_set(v,i,j,sp_cinit(rand(),rand()));
      sp_cmatrix_set(u,j,i,sp_cmatrix_get(v,i,j));
    }
  }
  sp_cmatrix_transpose(v);
  for(i = 0;i<sp_cmatrix_rows(v);i++){
    for(j = 0;j<sp_cmatrix_cols(v);j++){
      CuAssertComplexEquals(tc,sp_cmatrix_get(u,i,j),sp_cmatrix_get(v,i,j),fabs(REAL_EPSILON*sp_cabs(sp_cmatrix_get(u,i,j))));
    }
  }
  sp_cmatrix_free(u);
  sp_cmatrix_free(v);  
}

void test_sp_vector_dot_prod(CuTest * tc){
  sp_vector * v = sp_vector_alloc(2);
  sp_vector * u = sp_vector_alloc(2);
  sp_vector_set(v,0,1);
  sp_vector_set(v,1,3);
  sp_vector_set(u,0,2);
  sp_vector_set(u,1,4);
/*  v = [ 1 3 ]; u = [2 4];
    v.u = 14;*/
  CuAssertDblEquals(tc,sp_vector_dot_prod(u,v),14,fabs(REAL_EPSILON*14));
  CuAssertDblEquals(tc,sp_vector_dot_prod(v,u),14,fabs(REAL_EPSILON*14));
  sp_vector_free(u);
  sp_vector_free(v);  
}


void test_sp_cvector_dot_prod(CuTest * tc){
  sp_cvector * v = sp_cvector_alloc(2);
  sp_cvector * u = sp_cvector_alloc(2);
  sp_cvector_set(v,0,sp_cinit(1,0));
  sp_cvector_set(v,1,sp_cinit(3,0));
  sp_cvector_set(u,0,sp_cinit(2,0));
  sp_cvector_set(u,1,sp_cinit(4,0));
/*  v = [ 1 3 ]; u = [2 4];
    v.u = 14;*/
  CuAssertDblEquals(tc,sp_cabs(sp_cvector_dot_prod(u,v)),14,fabs(REAL_EPSILON*14));
  CuAssertDblEquals(tc,sp_cabs(sp_cvector_dot_prod(v,u)),14,fabs(REAL_EPSILON*14));
  sp_cvector_free(u);
  sp_cvector_free(v);  
}


void test_sp_vector_outer_prod(CuTest * tc){
  sp_vector * v = sp_vector_alloc(2);
  sp_vector * u = sp_vector_alloc(3);
  sp_3matrix * m;
  sp_vector_set(v,0,1);
  sp_vector_set(v,1,3);
  sp_vector_set(u,0,2);
  sp_vector_set(u,1,4);
  sp_vector_set(u,2,-1);
/*  v = [1 3];
    u = [2 4 -1];
    v x u = [ 2     4    -1 ]
            [ 6    12    -3 ]
*/
  m = sp_vector_outer_prod(v,u);
  
  CuAssertDblEquals(tc,sp_3matrix_get(m,0,0,0),2,fabs(REAL_EPSILON*2));
  CuAssertDblEquals(tc,sp_3matrix_get(m,0,1,0),4,fabs(REAL_EPSILON*4));
  CuAssertDblEquals(tc,sp_3matrix_get(m,0,2,0),-1,fabs(REAL_EPSILON*-1));
  CuAssertDblEquals(tc,sp_3matrix_get(m,1,0,0),6,fabs(REAL_EPSILON*6));
  CuAssertDblEquals(tc,sp_3matrix_get(m,1,1,0),12,fabs(REAL_EPSILON*12));
  CuAssertDblEquals(tc,sp_3matrix_get(m,1,2,0),-3,fabs(REAL_EPSILON*-3));
  sp_vector_free(u);
  sp_vector_free(v);  
  sp_3matrix_free(m);  
}


void test_sp_vector_norm(CuTest * tc){
  sp_vector * v = sp_vector_alloc(2);
  real n;
  sp_vector_set(v,0,1);
  sp_vector_set(v,1,3);
/*  v = [1 3];
    |v| = sqrt(10);

*/
  n = sp_vector_norm(v);  
  CuAssertDblEquals(tc,n,sqrt(10),fabs(REAL_EPSILON*n));
  sp_vector_free(v);  
}

void test_sp_cvector_norm(CuTest * tc){
  sp_cvector * v = sp_cvector_alloc(2);
  real n;
  sp_cvector_set(v,0,sp_cinit(1,0));
  sp_cvector_set(v,1,sp_cinit(3,0));
/*  v = [1 3];
    |v| = sqrt(10);

*/
  n = sp_cvector_norm(v);  
  CuAssertDblEquals(tc,n,sqrt(10),fabs(REAL_EPSILON*n));
  sp_cvector_free(v);  
}


void test_sp_matrix_vector_prod(CuTest * tc){
  sp_vector * v = sp_vector_alloc(2);
  sp_vector * u;
  sp_matrix * m = sp_matrix_alloc(3,2);

  sp_vector_set(v,0,1);
  sp_vector_set(v,1,3);
  sp_matrix_set(m,0,0,2);
  sp_matrix_set(m,0,1,6);
  sp_matrix_set(m,1,0,4);
  sp_matrix_set(m,1,1,12);
  sp_matrix_set(m,2,0,-1);
  sp_matrix_set(m,2,1,-3);
/*  v = [ 1  3 ];

    m = [  2  6 ]
        [  4 12 ]
        [ -1 -3	]
	
    mv = [ 20 40 -10 ]

*/
  u = sp_matrix_vector_prod(m,v);
  
  CuAssertTrue(tc,sp_vector_size(u)==3);
  CuAssertDblEquals(tc,sp_vector_get(u,0),20,fabs(REAL_EPSILON*20));
  CuAssertDblEquals(tc,sp_vector_get(u,1),40,fabs(REAL_EPSILON*40));
  CuAssertDblEquals(tc,sp_vector_get(u,2),-10,fabs(REAL_EPSILON*-10));

  sp_vector_free(u);
  sp_vector_free(v);  
  sp_matrix_free(m);  
}

void test_sp_cmatrix_cvector_prod(CuTest * tc){
  sp_cvector * v = sp_cvector_alloc(2);
  sp_cvector * u;
  sp_cmatrix * m = sp_cmatrix_alloc(3,2);
  sp_cvector_set(v,0,sp_cinit(1,0));
  sp_cvector_set(v,1,sp_cinit(3,0));
  sp_cmatrix_set(m,0,0,sp_cinit(2,0));
  sp_cmatrix_set(m,0,1,sp_cinit(6,0));
  sp_cmatrix_set(m,1,0,sp_cinit(4,0));
  sp_cmatrix_set(m,1,1,sp_cinit(12,0));
  sp_cmatrix_set(m,2,0,sp_cinit(-1,0));
  sp_cmatrix_set(m,2,1,sp_cinit(-3,0));
/*  v = [ 1  3 ];

    m = [  2  6 ]
        [  4 12 ]
        [ -1 -3	]
	
    mv = [ 20 40 -10 ]

*/
  u = sp_cmatrix_cvector_prod(m,v);
  
  CuAssertTrue(tc,sp_cvector_size(u)==3);
  CuAssertComplexEquals(tc,sp_cvector_get(u,0),sp_cinit(20,0),fabs(REAL_EPSILON*20));
  CuAssertComplexEquals(tc,sp_cvector_get(u,1),sp_cinit(40,0),fabs(REAL_EPSILON*40));
  CuAssertComplexEquals(tc,sp_cvector_get(u,2),sp_cinit(-10,0),fabs(REAL_EPSILON*-10));

  sp_cvector_free(u);
  sp_cvector_free(v);  
  sp_cmatrix_free(m);  
}


void test_sp_matrix_mul(CuTest * tc){

  sp_matrix * m = sp_matrix_alloc(3,2);
  sp_matrix * n = sp_matrix_alloc(2,3);
  sp_matrix * mn;
  sp_matrix * nm;
  sp_matrix_set(m,0,0,2);
  sp_matrix_set(m,0,1,6);
  sp_matrix_set(m,1,0,4);
  sp_matrix_set(m,1,1,12);
  sp_matrix_set(m,2,0,-1);
  sp_matrix_set(m,2,1,-3);

  sp_matrix_set(n,0,0,-2);
  sp_matrix_set(n,0,1,0);
  sp_matrix_set(n,0,2,1);
  sp_matrix_set(n,1,0,3);
  sp_matrix_set(n,1,1,4);
  sp_matrix_set(n,1,2,1);

/*  
    m = [  2  6 ]
        [  4 12 ]
        [ -1 -3	]
	
    n = [ -2 0 1 ]
        [  3 4 1 ]
	
    m*n = 14    24     8
          28    48    16
          -7   -12    -4

    n*m = -5   -15
          21    63
        

*/
  mn = sp_matrix_mul(m,n);
  nm = sp_matrix_mul(n,m);
  

  CuAssertDblEquals(tc,sp_matrix_get(mn,0,0),14,fabs(REAL_EPSILON*14));
  CuAssertDblEquals(tc,sp_matrix_get(mn,0,1),24,fabs(REAL_EPSILON*24));
  CuAssertDblEquals(tc,sp_matrix_get(mn,0,2),8,fabs(REAL_EPSILON*8));
  CuAssertDblEquals(tc,sp_matrix_get(mn,1,0),28,fabs(REAL_EPSILON*28));
  CuAssertDblEquals(tc,sp_matrix_get(mn,1,1),48,fabs(REAL_EPSILON*48));
  CuAssertDblEquals(tc,sp_matrix_get(mn,1,2),16,fabs(REAL_EPSILON*16));
  CuAssertDblEquals(tc,sp_matrix_get(mn,2,0),-7,fabs(REAL_EPSILON*-7));
  CuAssertDblEquals(tc,sp_matrix_get(mn,2,1),-12,fabs(REAL_EPSILON*-12));
  CuAssertDblEquals(tc,sp_matrix_get(mn,2,2),-4,fabs(REAL_EPSILON*-4));


  CuAssertDblEquals(tc,sp_matrix_get(nm,0,0),-5,fabs(REAL_EPSILON*-5));
  CuAssertDblEquals(tc,sp_matrix_get(nm,0,1),-15,fabs(REAL_EPSILON*-15));
  CuAssertDblEquals(tc,sp_matrix_get(nm,1,0),21,fabs(REAL_EPSILON*21));
  CuAssertDblEquals(tc,sp_matrix_get(nm,1,1),63,fabs(REAL_EPSILON*63));

  sp_matrix_free(m);  
  sp_matrix_free(n);  
  sp_matrix_free(mn);
  sp_matrix_free(nm);
}

void test_sp_cmatrix_mul(CuTest * tc){

  sp_cmatrix * m = sp_cmatrix_alloc(3,2);
  sp_cmatrix * n = sp_cmatrix_alloc(2,3);
  sp_cmatrix * mn;
  sp_cmatrix * nm;

  sp_cmatrix_set(m,0,0,sp_cinit(2,0));
  sp_cmatrix_set(m,0,1,sp_cinit(6,0));
  sp_cmatrix_set(m,1,0,sp_cinit(4,0));
  sp_cmatrix_set(m,1,1,sp_cinit(12,0));
  sp_cmatrix_set(m,2,0,sp_cinit(-1,0));
  sp_cmatrix_set(m,2,1,sp_cinit(-3,0));

  sp_cmatrix_set(n,0,0,sp_cinit(-2,0));
  sp_cmatrix_set(n,0,1,sp_cinit(0,0));
  sp_cmatrix_set(n,0,2,sp_cinit(1,0));
  sp_cmatrix_set(n,1,0,sp_cinit(3,0));
  sp_cmatrix_set(n,1,1,sp_cinit(4,0));
  sp_cmatrix_set(n,1,2,sp_cinit(1,0));

/*  
    m = [  2  6 ]
        [  4 12 ]
        [ -1 -3	]
	
    n = [ -2 0 1 ]
        [  3 4 1 ]
	
    m*n = 14    24     8
          28    48    16
          -7   -12    -4

    n*m = -5   -15
          21    63
        

*/
  mn = sp_cmatrix_mul(m,n);
  nm = sp_cmatrix_mul(n,m);
  

  CuAssertComplexEquals(tc,sp_cmatrix_get(mn,0,0),sp_cinit(14,0),fabs(REAL_EPSILON*14));
  CuAssertComplexEquals(tc,sp_cmatrix_get(mn,0,1),sp_cinit(24,0),fabs(REAL_EPSILON*24));
  CuAssertComplexEquals(tc,sp_cmatrix_get(mn,0,2),sp_cinit(8,0),fabs(REAL_EPSILON*8));
  CuAssertComplexEquals(tc,sp_cmatrix_get(mn,1,0),sp_cinit(28,0),fabs(REAL_EPSILON*28));
  CuAssertComplexEquals(tc,sp_cmatrix_get(mn,1,1),sp_cinit(48,0),fabs(REAL_EPSILON*48));
  CuAssertComplexEquals(tc,sp_cmatrix_get(mn,1,2),sp_cinit(16,0),fabs(REAL_EPSILON*16));
  CuAssertComplexEquals(tc,sp_cmatrix_get(mn,2,0),sp_cinit(-7,0),fabs(REAL_EPSILON*-7));
  CuAssertComplexEquals(tc,sp_cmatrix_get(mn,2,1),sp_cinit(-12,0),fabs(REAL_EPSILON*-12));
  CuAssertComplexEquals(tc,sp_cmatrix_get(mn,2,2),sp_cinit(-4,0),fabs(REAL_EPSILON*-4));


  CuAssertComplexEquals(tc,sp_cmatrix_get(nm,0,0),sp_cinit(-5,0),fabs(REAL_EPSILON*-5));
  CuAssertComplexEquals(tc,sp_cmatrix_get(nm,0,1),sp_cinit(-15,0),fabs(REAL_EPSILON*-15));
  CuAssertComplexEquals(tc,sp_cmatrix_get(nm,1,0),sp_cinit(21,0),fabs(REAL_EPSILON*21));
  CuAssertComplexEquals(tc,sp_cmatrix_get(nm,1,1),sp_cinit(63,0),fabs(REAL_EPSILON*63));

  sp_cmatrix_free(m);  
  sp_cmatrix_free(n);  
  sp_cmatrix_free(mn);
  sp_cmatrix_free(nm);
}


void test_sp_matrix_invert(CuTest * tc){
  int i,j;
  sp_matrix * m = sp_matrix_alloc(2,2);
  sp_matrix * n = sp_matrix_alloc(2,2);
  sp_matrix * id;
  for(i = 0;i<sp_matrix_rows(m);i++){
    for(j = 0;j<sp_matrix_cols(m);j++){
      sp_matrix_set(m,i,j,rand()%20);
    }
  }
  sp_matrix_memcpy(n,m);
  sp_matrix_invert(n);
  id = sp_matrix_mul(m,n);
  sp_matrix_set_identity(m);
  
  for(i = 0;i<sp_matrix_rows(m);i++){
    for(j = 0;j<sp_matrix_cols(m);j++){
      CuAssertDblEquals(tc,sp_matrix_get(m,i,j),sp_matrix_get(id,i,j),fabs(REAL_EPSILON*(sp_matrix_get(m,i,j))+sqrt(REAL_EPSILON)));
    }
  }
  sp_matrix_free(m);
  sp_matrix_free(n);
  sp_matrix_free(id);
}

void test_sp_cmatrix_invert(CuTest * tc){
  int i,j;
  sp_cmatrix * m = sp_cmatrix_alloc(2,2);
  sp_cmatrix * n = sp_cmatrix_alloc(2,2);
  sp_cmatrix * id;
  for(i = 0;i<sp_cmatrix_rows(m);i++){
    for(j = 0;j<sp_cmatrix_cols(m);j++){
      sp_cmatrix_set(m,i,j,sp_cinit(rand()%20,(rand()%20)));
    }
  }
  sp_cmatrix_memcpy(n,m);
  sp_cmatrix_invert(n);
  id = sp_cmatrix_mul(m,n);
  sp_cmatrix_set_identity(m);
  
  for(i = 0;i<sp_cmatrix_rows(m);i++){
    for(j = 0;j<sp_cmatrix_cols(m);j++){
      CuAssertComplexEquals(tc,sp_cmatrix_get(m,i,j),sp_cmatrix_get(id,i,j),fabs(REAL_EPSILON*sp_cabs((sp_cmatrix_get(m,i,j)))+REAL_EPSILON*2*2));
    }
  }
  sp_cmatrix_free(m);
  sp_cmatrix_free(n);
  sp_cmatrix_free(id);
}


CuSuite* linear_alg_get_suite(void)
{
  CuSuite* suite = CuSuiteNew();
  

  SUITE_ADD_TEST(suite, test_sp_min);
  SUITE_ADD_TEST(suite, test_sp_max);

  SUITE_ADD_TEST(suite, test_sp_vector_alloc);
  SUITE_ADD_TEST(suite, test_sp_cvector_alloc);
  SUITE_ADD_TEST(suite, test_sp_matrix_alloc);
  SUITE_ADD_TEST(suite, test_sp_3matrix_alloc);
  SUITE_ADD_TEST(suite, test_sp_cmatrix_alloc);
  SUITE_ADD_TEST(suite, test_sp_c3matrix_alloc);
  SUITE_ADD_TEST(suite, test_sp_i3matrix_alloc);

  SUITE_ADD_TEST(suite, test_sp_vector_set_get);
  SUITE_ADD_TEST(suite, test_sp_cvector_set_get);
  SUITE_ADD_TEST(suite, test_sp_matrix_set_get);
  SUITE_ADD_TEST(suite, test_sp_3matrix_set_get);
  SUITE_ADD_TEST(suite, test_sp_cmatrix_set_get);
  SUITE_ADD_TEST(suite, test_sp_c3matrix_set_get);
  SUITE_ADD_TEST(suite, test_sp_imatrix_set_get);
  SUITE_ADD_TEST(suite, test_sp_i3matrix_set_get);

  SUITE_ADD_TEST(suite, test_sp_vector_memcpy);
  SUITE_ADD_TEST(suite, test_sp_cvector_memcpy);
  SUITE_ADD_TEST(suite, test_sp_matrix_memcpy);
  SUITE_ADD_TEST(suite, test_sp_3matrix_memcpy);
  SUITE_ADD_TEST(suite, test_sp_cmatrix_memcpy);
  SUITE_ADD_TEST(suite, test_sp_c3matrix_memcpy);

  SUITE_ADD_TEST(suite, test_sp_vector_add);
  SUITE_ADD_TEST(suite, test_sp_cvector_add);
  SUITE_ADD_TEST(suite, test_sp_matrix_add);
  SUITE_ADD_TEST(suite, test_sp_3matrix_add);
  SUITE_ADD_TEST(suite, test_sp_cmatrix_add);
  SUITE_ADD_TEST(suite, test_sp_c3matrix_add);

  SUITE_ADD_TEST(suite, test_sp_vector_sub);
  SUITE_ADD_TEST(suite, test_sp_cvector_sub);
  SUITE_ADD_TEST(suite, test_sp_matrix_sub);
  SUITE_ADD_TEST(suite, test_sp_3matrix_sub);
  SUITE_ADD_TEST(suite, test_sp_cmatrix_sub);
  SUITE_ADD_TEST(suite, test_sp_c3matrix_sub);

  SUITE_ADD_TEST(suite, test_sp_vector_mul);
  SUITE_ADD_TEST(suite, test_sp_cvector_mul);
  SUITE_ADD_TEST(suite, test_sp_matrix_mul_elements);
  SUITE_ADD_TEST(suite, test_sp_3matrix_mul_elements);
  SUITE_ADD_TEST(suite, test_sp_cmatrix_mul_elements);
  SUITE_ADD_TEST(suite, test_sp_c3matrix_mul_elements);

  
  SUITE_ADD_TEST(suite, test_sp_vector_div);
  SUITE_ADD_TEST(suite, test_sp_cvector_div);
  SUITE_ADD_TEST(suite, test_sp_matrix_div_elements);
  SUITE_ADD_TEST(suite, test_sp_3matrix_div_elements);
  SUITE_ADD_TEST(suite, test_sp_cmatrix_div_elements);
  SUITE_ADD_TEST(suite, test_sp_c3matrix_div_elements);

  SUITE_ADD_TEST(suite, test_sp_vector_scale);
  SUITE_ADD_TEST(suite, test_sp_cvector_scale);
  SUITE_ADD_TEST(suite, test_sp_matrix_scale);
  SUITE_ADD_TEST(suite, test_sp_3matrix_scale);
  SUITE_ADD_TEST(suite, test_sp_cmatrix_scale);
  SUITE_ADD_TEST(suite, test_sp_c3matrix_scale);

  SUITE_ADD_TEST(suite, test_sp_vector_add_constant);
  SUITE_ADD_TEST(suite, test_sp_cvector_add_constant);
  SUITE_ADD_TEST(suite, test_sp_matrix_add_constant);
  SUITE_ADD_TEST(suite, test_sp_3matrix_add_constant);
  SUITE_ADD_TEST(suite, test_sp_cmatrix_add_constant);
  SUITE_ADD_TEST(suite, test_sp_c3matrix_add_constant);

  SUITE_ADD_TEST(suite, test_sp_vector_dot_prod);
  SUITE_ADD_TEST(suite, test_sp_cvector_dot_prod);

  SUITE_ADD_TEST(suite, test_sp_vector_outer_prod);

  SUITE_ADD_TEST(suite, test_sp_vector_norm);
  SUITE_ADD_TEST(suite, test_sp_cvector_norm);

  SUITE_ADD_TEST(suite, test_sp_matrix_vector_prod);
  SUITE_ADD_TEST(suite, test_sp_cmatrix_cvector_prod);

  SUITE_ADD_TEST(suite, test_sp_matrix_transpose);
  SUITE_ADD_TEST(suite, test_sp_cmatrix_transpose);

  SUITE_ADD_TEST(suite, test_sp_matrix_mul);
  SUITE_ADD_TEST(suite, test_sp_cmatrix_mul);

  SUITE_ADD_TEST(suite, test_sp_matrix_invert);
  SUITE_ADD_TEST(suite, test_sp_cmatrix_invert);
  
  return suite;
}


void test_sp_image_edge_extend(CuTest * tc){
  Image * a = sp_image_alloc(2,2,1);
  sp_image_set(a,0,0,0,sp_cinit(1,0));
  sp_image_set(a,1,0,0,sp_cinit(2,0));
  sp_image_set(a,0,1,0,sp_cinit(3,0));
  sp_image_set(a,1,1,0,sp_cinit(4,0));
  Image * b = sp_image_edge_extend(a,1,SP_ZERO_PAD_EDGE,SP_2D);
  for(int x = 0;x<sp_image_x(b);x++){
    for(int y = 0;y<sp_image_y(b);y++){
      if(!x || x == 3 || !y || y == 3){
	CuAssertComplexEquals(tc,sp_image_get(b,x,y,0),sp_cinit(0,0),fabs(REAL_EPSILON*(sp_cabs(sp_image_get(b,x,y,0)))+REAL_EPSILON));
      }
    }
  }
  sp_image_free(b);
  b = sp_image_edge_extend(a,1,SP_SYMMETRIC_EDGE,SP_2D);
  CuAssertComplexEquals(tc,sp_image_get(b,0,0,0),sp_image_get(a,0,0,0),fabs(REAL_EPSILON*(sp_cabs(sp_image_get(b,0,0,0)))+REAL_EPSILON));
  CuAssertComplexEquals(tc,sp_image_get(b,2,0,0),sp_image_get(a,1,0,0),fabs(REAL_EPSILON*sp_cabs(sp_image_get(b,2,0,0))+REAL_EPSILON));
  sp_image_free(b);
  b = sp_image_edge_extend(a,1,SP_REPLICATE_EDGE,SP_2D);
  CuAssertComplexEquals(tc,sp_image_get(b,0,0,0),sp_image_get(a,0,0,0),fabs(REAL_EPSILON*sp_cabs(sp_image_get(b,0,0,0))+REAL_EPSILON));
  CuAssertComplexEquals(tc,sp_image_get(b,2,0,0),sp_image_get(a,1,0,0),fabs(REAL_EPSILON*sp_cabs(sp_image_get(b,2,0,0))+REAL_EPSILON));
  sp_image_free(b);
  b = sp_image_edge_extend(a,1,SP_CIRCULAR_EDGE,SP_2D);
  CuAssertComplexEquals(tc,sp_image_get(b,0,0,0),sp_image_get(a,1,1,0),fabs(REAL_EPSILON*sp_cabs(sp_image_get(b,0,0,0))+REAL_EPSILON));
  CuAssertComplexEquals(tc,sp_image_get(b,2,0,0),sp_image_get(a,1,1,0),fabs(REAL_EPSILON*sp_cabs(sp_image_get(b,2,0,0))+REAL_EPSILON));
  sp_image_free(b);
  sp_image_free(a);
}

void test_sp_bubble_sort(CuTest * tc){
  real array[4] = {4,-5.3, 1000,2};
  sp_bubble_sort(array,4);
  CuAssertDblEquals(tc,array[0],-5.3,fabs(REAL_EPSILON*(-5.3)+REAL_EPSILON));
  CuAssertDblEquals(tc,array[1],2,fabs(REAL_EPSILON*(2)+REAL_EPSILON));
  CuAssertDblEquals(tc,array[2],4,fabs(REAL_EPSILON*(4)+REAL_EPSILON));
  CuAssertDblEquals(tc,array[3],1000,fabs(REAL_EPSILON*(1000)+REAL_EPSILON));  
}


void test_sp_image_median_filter(CuTest * tc){
  Image * a = sp_image_alloc(2,2,1);
  sp_image_set(a,0,0,0,sp_cinit(1,0));
  sp_image_set(a,1,0,0,sp_cinit(2,0));
  sp_image_set(a,0,1,0,sp_cinit(3,0));
  sp_image_set(a,1,1,0,sp_cinit(4,0));
  sp_i3matrix * kernel = sp_i3matrix_alloc(2,2,1);
  sp_i3matrix_add_constant(kernel,2);
  sp_image_median_filter(a,kernel,SP_ZERO_PAD_EDGE,SP_2D);
  CuAssertDblEquals(tc,sp_cabs(sp_image_get(a,0,0,0)),0.0,fabs(REAL_EPSILON*(0)+REAL_EPSILON));
  CuAssertDblEquals(tc,sp_cabs(sp_image_get(a,1,1,0)),2.5,fabs(REAL_EPSILON*(2.5)+REAL_EPSILON));
  sp_image_free(a);
  sp_i3matrix_free(kernel);
}

void test_sp_image_max(CuTest * tc){
  int x = 5;
  int y = 5;
  int z = 5;
  Image * img = sp_image_alloc(x,y,z);
  for(long long i = 0;i<sp_image_size(img);i++){
    img->image->data[i] = sp_cinit(0,0);
  }
  real max = sp_image_max(img,NULL,NULL,NULL,NULL);
  CuAssertDblEquals(tc,max,0,(0+1)*fabs(REAL_EPSILON));
  long long ind = rand()%sp_image_size(img);
  img->image->data[ind] = sp_cinit(1,0);
  long long ind2;
  max = sp_image_max(img,&ind2,NULL,NULL,NULL);
  CuAssertDblEquals(tc,max,1,(1+1)*fabs(REAL_EPSILON));
  CuAssertIntEquals(tc,ind,ind2);
  sp_image_free(img);
}

void test_sp_image_gaussian_blur(CuTest * tc){
  Image * a = sp_image_alloc(11,11,1);
  sp_image_fill(a,sp_cinit(0,0));
  sp_image_set(a,5,5,0,sp_cinit(1,0));
  Image * b = gaussian_blur(a,1);
  /* Check if integral over image is mantained */
  Complex a_sum = {0,0};
  Complex b_sum = {0,0};
  for(int i = 0;i<sp_image_size(a);i++){
    sp_cincr(a_sum,a->image->data[i]);
    sp_cincr(b_sum,b->image->data[i]);
  }
  /* the 0.00015 follows from the fact that 99.7% of the gaussian is within 3
     standard deviation and the gaussian blur only calculates 3 standard deviations. */
  CuAssertComplexEquals(tc,a_sum,b_sum,fabs(0.00015*(sp_cabs(a_sum))));
  long long ind1;
  long long ind2;
  /* Image maximum must remain in the same place */
  sp_image_max(a,&ind1,NULL,NULL,NULL);
  sp_image_max(b,&ind2,NULL,NULL,NULL);
  CuAssertIntEquals(tc,ind1,ind2);
  sp_image_free(a);
  sp_image_free(b);
  a = sp_image_alloc(11,11,11);
  sp_image_fill(a,sp_cinit(0,0));
  sp_image_set(a,5,5,5,sp_cinit(1,0));
  b = gaussian_blur(a,1);
  /* Check if integral over image is mantained */
  a_sum = sp_cinit(0,0);
  b_sum = sp_cinit(0,0);
  for(int i = 0;i<sp_image_size(a);i++){
    sp_cincr(a_sum,a->image->data[i]);
    sp_cincr(b_sum,b->image->data[i]);
  }
  /* the 0.00015 follows from the fact that 99.7% of the gaussian is within 3
     standard deviation and the gaussian blur only calculates 3 standard deviations. */
  CuAssertComplexEquals(tc,a_sum,b_sum,fabs(0.00015*(sp_cabs(a_sum))));
  /* Image maximum must remain in the same place */
  sp_image_max(a,&ind1,NULL,NULL,NULL);
  sp_image_max(b,&ind2,NULL,NULL,NULL);
  CuAssertIntEquals(tc,ind1,ind2);
  sp_image_free(a);
  sp_image_free(b);
  
}

void test_cube_crop(CuTest * tc){
  Image * a = sp_image_alloc(3,3,3);
  for(long long i = 0;i<sp_image_size(a);i++){
    a->image->data[i] = sp_cinit(2,0);
  }
  sp_image_set(a,2,1,1,sp_cinit(3,0));
  Image *b = cube_crop(a,2,1,1,2,1,1);
  CuAssertComplexEquals(tc,sp_image_get(b,0,0,0),sp_cinit(3,0),3*fabs(REAL_EPSILON));  
}
  


void test_sp_image_h5_read_write(CuTest * tc){
  /* First lets check 2D */
  int size = 100;
  Image * a;
  Image * b;
  a = sp_image_alloc(size,size,1);
  a->phased = 1;
  for(int z = 0; z<sp_image_z(a);z++){
    for(int y = 0; y<sp_image_y(a);y++){
      for(int x = 0; x<sp_image_x(a);x++){
	sp_image_set(a,0,0,0,sp_cinit(rand()%50,rand()%50));
      }
    }
  }
  sp_image_write(a,"test.h5",sizeof(float));
  b = sp_image_read("test.h5",0);
  for(int z = 0; z<sp_image_z(a);z++){
    for(int y = 0; y<sp_image_y(a);y++){
      for(int x = 0; x<sp_image_x(a);x++){
	CuAssertComplexEquals(tc,sp_image_get(b,x,y,z),sp_image_get(a,x,y,z),fabs(sp_cabs(sp_image_get(a,x,y,z))*(REAL_EPSILON))+REAL_EPSILON);
      }
    }
  }
  sp_image_free(a);
  sp_image_free(b);
  remove("test.h5");

  a = sp_image_alloc(size,size,1);
  a->phased = 1;
  for(int z = 0; z<sp_image_z(a);z++){
    for(int y = 0; y<sp_image_y(a);y++){
      for(int x = 0; x<sp_image_x(a);x++){
	sp_image_set(a,0,0,0,sp_cinit(rand()%50,rand()%50));
      }
    }
  }

  /* Now with double precision */
  sp_image_write(a,"test.h5",sizeof(double));
  b = sp_image_read("test.h5",0);
  for(int z = 0; z<sp_image_z(a);z++){
    for(int y = 0; y<sp_image_y(a);y++){
      for(int x = 0; x<sp_image_x(a);x++){
	CuAssertComplexEquals(tc,sp_image_get(b,x,y,z),sp_image_get(a,x,y,z),fabs(sp_cabs(sp_image_get(a,x,y,z))*(REAL_EPSILON))+REAL_EPSILON);
      }
    }
  }
  sp_image_free(a);
  sp_image_free(b);
  remove("test.h5");

  /* And now 3D */
  size = 10;
  a = sp_image_alloc(size,size,size);
  a->phased = 1;
  for(int z = 0; z<sp_image_z(a);z++){
    for(int y = 0; y<sp_image_y(a);y++){
      for(int x = 0; x<sp_image_x(a);x++){
	sp_image_set(a,0,0,0,sp_cinit(rand()%50,rand()%50));
      }
    }
  }
  sp_image_write(a,"test.h5",sizeof(float));
  b = sp_image_read("test.h5",0);
  for(int z = 0; z<sp_image_z(a);z++){
    for(int y = 0; y<sp_image_y(a);y++){
      for(int x = 0; x<sp_image_x(a);x++){
	CuAssertComplexEquals(tc,sp_image_get(b,x,y,z),sp_image_get(a,x,y,z),fabs(sp_cabs(sp_image_get(a,x,y,z))*(REAL_EPSILON))+REAL_EPSILON);
      }
    }
  }
  sp_image_free(a);
  sp_image_free(b);
  remove("test.h5");


  /* And in double precision */
  size = 10;
  a = sp_image_alloc(size,size,size);
  a->phased = 1;
  for(int z = 0; z<sp_image_z(a);z++){
    for(int y = 0; y<sp_image_y(a);y++){
      for(int x = 0; x<sp_image_x(a);x++){
	sp_image_set(a,0,0,0,sp_cinit(rand()%50,rand()%50));
      }
    }
  }
  sp_image_write(a,"test.h5",sizeof(double));
  b = sp_image_read("test.h5",0);
  for(int z = 0; z<sp_image_z(a);z++){
    for(int y = 0; y<sp_image_y(a);y++){
      for(int x = 0; x<sp_image_x(a);x++){
	CuAssertComplexEquals(tc,sp_image_get(b,x,y,z),sp_image_get(a,x,y,z),fabs(sp_cabs(sp_image_get(a,x,y,z))*(REAL_EPSILON))+REAL_EPSILON);
      }
    }
  }
  sp_image_free(a);
  sp_image_free(b);
  remove("test.h5");
}


CuSuite* image_get_suite(void)
{
  CuSuite* suite = CuSuiteNew();
  
  SUITE_ADD_TEST(suite, test_sp_image_edge_extend);
  SUITE_ADD_TEST(suite, test_sp_bubble_sort);
  SUITE_ADD_TEST(suite, test_sp_image_median_filter);
  SUITE_ADD_TEST(suite, test_sp_image_max);
  SUITE_ADD_TEST(suite,test_sp_image_gaussian_blur);
  SUITE_ADD_TEST(suite,test_cube_crop);
  SUITE_ADD_TEST(suite,test_sp_image_h5_read_write);
  return suite;
}

int RunAllTests(void)
{
	CuString *output = CuStringNew();
	CuSuite* suite = CuSuiteNew();

	CuSuiteAddSuite(suite, linear_alg_get_suite());
	CuSuiteAddSuite(suite, image_get_suite());

	CuSuiteRun(suite);
	CuSuiteSummary(suite, output);
	CuSuiteDetails(suite, output);
	printf("%s\n", output->buffer);
	return suite->failCount;
}

int main(void)
{
	return RunAllTests();
}
