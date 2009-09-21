
#include "AllTests.h"
#include <gsl/gsl_vector_float.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>

Image * sp_image_cuda_ifft(Image * img);

static Complex czero = {0,0};

void test_sp_clog(CuTest * tc){
  Complex a = sp_cinit(1,1);
  Complex log_a = sp_clog(a);
  CuAssertComplexEquals(tc,log_a,sp_cinit(log(2)/2.0,M_PI/4.0),REAL_EPSILON);  
}

void test_sp_cexp(CuTest * tc){
  Complex a = sp_cinit(1,1);
  Complex exp_a = sp_cexp(a);
  CuAssertComplexEquals(tc,exp_a,sp_cinit(exp(1.0)*cos(1.0),exp(1.0)*sin(1.0)),5*REAL_EPSILON);  
}

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

  /* test performance */
  int total_time = 0;
  v = sp_vector_alloc(1024*1024);
  u = sp_vector_alloc(1024*1024);
  int ntimes = 200;
  for(i = 0;i<sp_vector_size(v);i++){
    sp_vector_set(v,i,rand());
    sp_vector_set(u,i,rand());
  }
  for(int j = 0;j<ntimes;j++){  
    int timer = sp_timer_start();
    //    sp_vector_add(v,u);
    sp_vector_dot_prod(v,u);
    total_time += sp_timer_stop(timer);
  }
  printf("Average time for vector(1024*1024) dot product is %e ms\n",total_time/1000.0/ntimes);

  gsl_vector_float * gsl_v = gsl_vector_float_alloc(1024*1024);
  gsl_vector_float * gsl_u = gsl_vector_float_alloc(1024*1024);
  for(i = 0;i<gsl_u->size;i++){
    gsl_vector_float_set(gsl_v,i,rand());
    gsl_vector_float_set(gsl_u,i,rand());
  }
  total_time = 0;
  for(int j = 0;j<ntimes;j++){     
    int timer = sp_timer_start();
    float res;
    gsl_blas_sdot(gsl_v,gsl_u,&res);
    //    cblas_saxpy(gsl_u->size,2.0,gsl_v->data,1,gsl_u->data,1);
    total_time += sp_timer_stop(timer);
  }
  printf("Average time for gsl vector(1024*1024) dot product is %e ms\n",total_time/1000.0/ntimes);

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
	CuAssertComplexEquals(tc,sp_c3matrix_get(u,i,j,k),sp_cscale(sp_c3matrix_get(v,i,j,k),x),fabs(REAL_EPSILON*sp_cabs(sp_c3matrix_get(u,i,j,k)))+REAL_EPSILON);
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


void test_sp_create_spline2_kernel_table(CuTest * tc){
  real tol = 1e-5;
  sp_kernel * k = sp_create_spline2_kernel_table(4.1,tol);
  real step = 1e-6;
  double max_error = 0;
  double max_error_r = 0;
  for(double r = 0;r<1.5;r+=step){
    double r2 = r*r;
    double sample = sp_kernel_table_sample(k,r2);
    double f = (r+1.5)*(r+1.5) - 3.0*(r+0.5)*(r+0.5);
     if (r > 0.5) f += 3.0*(r-0.5)*(r-0.5);
    if (r > 1.5) f -= (r-1.5)*(r-1.5);	
    double error = fabs(f-sample);
    if(error  > max_error){
      max_error = error;
      max_error_r = r;
    }    
  }
#ifndef NDEBUG
  printf("Spline2 R2 Table Max Error = %g at r = %f for tolerance = %g\n",max_error,max_error_r,tol);
#endif
  CuAssertTrue(tc,max_error < tol*1.1);

  /* Do a control run */
  int t = sp_timer_start();
  step = 1e-5;
  int i = 0;
  /* I need double precision or the small step is gonna be washed in the numerical error */
  for(double r2 = 0;r2<4;r2+=step){
    real sample = 1;
    sample = i;
    i++;
  }
  long long int control_dt = sp_timer_stop(t);


  t = sp_timer_start();
  i = 0;
  /* I need double precision or the small step is gonna be washed in the numerical error */
  for(double r2 = 0;r2<4;r2+=step){
    real sample = sp_kernel_table_sample(k,r2);
    sample = i;
    i++;
  }
  long long int dt = sp_timer_stop(t)-control_dt;
#ifndef NDEBUG
  printf("Spline R2 %d steps in %lld micro seconds %g steps/us\n",i,dt,(real)i/dt);
#endif
}


real spline_interpolation2(Image * a, real x1, real y1,real z1)
{
  int x,y,z;
  real r;
  real w;
  real w_tot = 0.0;
  real res = 0.0;
  
  int x_min = MAX((int)(x1+0.5)-2,0);
  int x_max = MIN((int)(x1+0.5)+2,sp_image_x(a)-1);
  int y_min = MAX((int)(y1+0.5)-2,0);
  int y_max = MIN((int)(y1+0.5)+2,sp_image_y(a)-1);
  int z_min = MAX((int)(z1+0.5)-2,0);
  int z_max = MIN((int)(z1+0.5)+2,sp_image_z(a)-1);
  for (z = z_min; z <= z_max; z++) {
    for (y = y_min; y <= y_max; y++) {
      for (x = x_min; x <= x_max; x++) {
	r = sqrt(((real)x - x1)*((real)x - x1) + ((real)y - y1)*((real)y - y1)+((real)z - z1)*((real)z - z1));
	if (r < 2.0) {
	  w = (r+1.5)*(r+1.5) - 3.0*(r+0.5)*(r+0.5);
	  if (r > 0.5) w += 3.0*(r-0.5)*(r-0.5);
	  if (r > 1.5) w -= (r-1.5)*(r-1.5);
	  res += sp_cabs(sp_image_get(a,x,y,z))*w;
	  w_tot += w;
	}
      }
    }
  }
  if (w_tot != 0.0 && res) return res /= w_tot;
  else return 0.0;
}

void test_sp_c3matrix_kernel_interpolation(CuTest * tc){
  real tolerance = 1e-5;
  Image * a = sp_image_alloc(4,4,4);
  sp_c3matrix * m = a->image;
  int ntests = 1000;
  double max_error = 0;
  sp_kernel * k = sp_create_spline2_kernel_table(3*2.5*2.5,tolerance);
  for(int i = 0;i<sp_c3matrix_size(m);i++){
    m->data[i] = sp_cinit(p_drand48()-0.5,p_drand48()-0.5);
  }
  for(int i = 0;i<ntests;i++){
    real x1 = p_drand48()*3;
    real y1 = p_drand48()*3;
    real z1 = p_drand48()*3;
    real f = sp_c3matrix_kernel_interpolation(m,x1,y1,z1,k);
    real f2 = spline_interpolation2(a,x1,y1,z1);
    if(fabs(f-f2) > max_error){
      max_error = fabs(f-f2);
    }
    CuAssertDblEquals(tc,f,f2,tolerance*(f+f2)/2.0);
  }
#ifndef NDEBUG
  printf("Max error between spline interpolation and c3matrix_kernel_interpolation with a spline kernel = %g using a tolerance = %g\n",max_error,tolerance);
#endif
}


void test_sp_c3matrix_rotate(CuTest * tc){
  sp_c3matrix * a = sp_c3matrix_alloc(2,2,1);
  sp_c3matrix_set(a,0,0,0,sp_cinit(1,0));
  sp_c3matrix_set(a,1,0,0,sp_cinit(2,0));
  sp_c3matrix_set(a,1,1,0,sp_cinit(3,0));
  sp_c3matrix_set(a,0,1,0,sp_cinit(4,0));
  sp_c3matrix * b = sp_c3matrix_rotate(a,sp_ZAxis,sp_90Degrees,0);
  CuAssertComplexEquals(tc,sp_c3matrix_get(b,0,0,0),sp_cinit(2,0),REAL_EPSILON);  
  sp_c3matrix_free(b);
  b = sp_c3matrix_rotate(a,sp_ZAxis,sp_180Degrees,0);
  CuAssertComplexEquals(tc,sp_c3matrix_get(b,0,0,0),sp_cinit(3,0),REAL_EPSILON);  
  sp_c3matrix_free(b);
  b = sp_c3matrix_rotate(a,sp_ZAxis,sp_270Degrees,0);
  CuAssertComplexEquals(tc,sp_c3matrix_get(b,0,0,0),sp_cinit(4,0),REAL_EPSILON);  
  sp_c3matrix_free(b);
  b = sp_c3matrix_rotate(a,sp_ZAxis,sp_90Degrees,0);
  sp_c3matrix * c = sp_c3matrix_rotate(b,sp_ZAxis,sp_270Degrees,0);
  for(int i =0;i<sp_c3matrix_size(a);i++){
    CuAssertComplexEquals(tc,a->data[i],c->data[i],REAL_EPSILON);  
  }
  sp_c3matrix_free(a);
  sp_c3matrix_free(b);
  sp_c3matrix_free(c);
}

void test_sp_matrix_rotate(CuTest * tc){
  sp_matrix * a = sp_matrix_alloc(2,2);
  sp_matrix_set(a,0,0,1);
  sp_matrix_set(a,1,0,2);
  sp_matrix_set(a,1,1,3);
  sp_matrix_set(a,0,1,4);
  sp_matrix * b = sp_matrix_rotate(a,sp_90Degrees,0);
  CuAssertDblEquals(tc,sp_matrix_get(b,0,0),2,REAL_EPSILON);  
  sp_matrix_free(b);
  b = sp_matrix_rotate(a,sp_180Degrees,0);
  CuAssertDblEquals(tc,sp_matrix_get(b,0,0),3,REAL_EPSILON);  
  sp_matrix_free(b);
  b = sp_matrix_rotate(a,sp_270Degrees,0);
  CuAssertDblEquals(tc,sp_matrix_get(b,0,0),4,REAL_EPSILON);  
  sp_matrix_free(b);
  b = sp_matrix_rotate(a,sp_90Degrees,0);
  sp_matrix * c = sp_matrix_rotate(b,sp_270Degrees,0);
  for(int i =0;i<sp_matrix_size(a);i++){
    CuAssertDblEquals(tc,a->data[i],c->data[i],REAL_EPSILON);  
  }
  sp_matrix_free(a);
  sp_matrix_free(b);
  sp_matrix_free(c);
}

void test_sp_image_cuda_ifft(CuTest * tc){
  Image * a = sp_image_alloc(512,512,1);
  for(int i = 0;i<sp_image_size(a);i++){
    a->image->data[0] = sp_cinit((float)rand()/RAND_MAX,(float)rand()/RAND_MAX);
  }
  a->phased = 1;
  a->shifted = 1;
  Image * cuda_out = sp_image_cuda_ifft(a);
  Image * fftw_out = sp_image_ifft(a);
  double error = 0;
  double sum = 0;
  for(int i= 0;i<sp_image_size(a);i++){
    error += sp_cabs(sp_csub(cuda_out->image->data[i],fftw_out->image->data[i]));
    sum += sp_cabs(fftw_out->image->data[i]);
  }
  printf("CUDA ifft tot. error = %e\n",error);
  printf("CUDA ifft rel. error = %e\n",error/sum);
}

void test_sp_image_cuda_fft(CuTest * tc){
  Image * a = sp_image_alloc(2048,2048,1);
  for(int i = 0;i<sp_image_size(a);i++){
    a->image->data[0] = sp_cinit((float)rand()/RAND_MAX,(float)rand()/RAND_MAX);
  }
  a->phased = 1;
  Image * cuda_out = sp_image_cuda_fft(a);
  Image * fftw_out = sp_image_fft(a);
  double error = 0;
  double sum = 0;
  for(int i= 0;i<sp_image_size(a);i++){
    error += sp_cabs(sp_csub(cuda_out->image->data[i],fftw_out->image->data[i]));
    sum += sp_cabs(fftw_out->image->data[i]);
  }
  printf("CUDA fft tot. error = %e\n",error);
  printf("CUDA fft rel. error = %e\n",error/sum);
}

CuSuite* linear_alg_get_suite(void)
{
  CuSuite* suite = CuSuiteNew();

  SUITE_ADD_TEST(suite, test_sp_clog);  
  SUITE_ADD_TEST(suite, test_sp_cexp);  
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

  //  SUITE_ADD_TEST(suite, test_sp_matrix_invert);
  //  SUITE_ADD_TEST(suite, test_sp_cmatrix_invert);

  SUITE_ADD_TEST(suite,test_sp_create_spline2_kernel_table);
  SUITE_ADD_TEST(suite,test_sp_c3matrix_kernel_interpolation);

  SUITE_ADD_TEST(suite,test_sp_c3matrix_rotate);
  SUITE_ADD_TEST(suite,test_sp_matrix_rotate);

  SUITE_ADD_TEST(suite,test_sp_image_cuda_ifft);
  SUITE_ADD_TEST(suite,test_sp_image_cuda_fft);

  return suite;
}
