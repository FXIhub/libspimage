#include <stdlib.h>
#include <string.h>
#include "spimage.h"

sp_vector * sp_vector_alloc(const int size){
  sp_vector * ret = malloc(sizeof(sp_vector));
  ret->size = size;
  ret->data = calloc(size,sizeof(real));
}


sp_cvector * sp_cvector_alloc(const int size){
  sp_cvector * ret = malloc(sizeof(sp_cvector));
  ret->size = size;
  ret->data = calloc(size,sizeof(complex));
}

void sp_vector_free(sp_vector * v){
  free(v->data);
  free(v);
}

void sp_cvector_free(sp_cvector * v){
  free(v->data);
  free(v);
}




sp_matrix * sp_matrix_alloc(unsigned int nrows, unsigned int ncols){
  sp_matrix * res = malloc(sizeof(sp_matrix));
  res->rows = nrows;
  res->cols = ncols;
  res->data = calloc(nrows*ncols,sizeof(real));
  return res;
}

void sp_matrix_free(sp_matrix * a){
  free(a->data);
  free(a);
}




void sp_matrix_invert(sp_matrix * a){
  /* Try to invert a matrix by gaussian elimination */
  int n = a->cols;
  double t,ten;
  int i,j,k,kb,kp1,l,nm1;
  real * ipvt = calloc(n,sizeof(real));
  real * work = calloc(n,sizeof(real));
  
  
  for (k = 0; k < n; k++) {
    sp_matrix_set(a,k,k,1.0 /sp_matrix_get(a,k,k) );
    t = -sp_matrix_get(a,k,k);
    /* dscal(k-1,t,a(1,k),1) */
    for (i = 0; i < k; i++){
      sp_matrix_set(a,i,k,sp_matrix_get(a,i,k)*t);
    }
    kp1 = k + 1;
    if ( n > kp1) {
      for (j = kp1; j < n; j++) {
	t = sp_matrix_get(a,k,j);
	sp_matrix_set(a,k,j,0);
	/* daxpy(k,t,a(1,k),1,a(1,j),1) */
	for (i = 0; i < k+1; i++){
	  sp_matrix_set(a,i,j,sp_matrix_get(a,i,j) + t * sp_matrix_get(a,i,k));
	}
      }
    }
  }
}
