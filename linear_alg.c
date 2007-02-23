#include <stdlib.h>
#include <string.h>
#include "spimage.h"

sp_vector * sp_vector_alloc(const int size){
  sp_vector * ret = malloc(sizeof(sp_vector));
  ret->size = size;
  ret->data = calloc(size,sizeof(real));
  return ret;
}


sp_cvector * sp_cvector_alloc(const int size){
  sp_cvector * ret = malloc(sizeof(sp_cvector));
  ret->size = size;
  ret->data = calloc(size,sizeof(complex));
  return ret;
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






void sp_matrix_invert(sp_matrix * m){
  int n,i,j;
  real x;
  sp_matrix * inv = sp_matrix_alloc(m->rows,m->cols);
  sp_matrix_set_identity(inv);
  /* triangularize the matrix */
  /* For every row */
  for(i = 0;i<m->rows;i++){
    /* set leading element to 1 */
    x = 1.0/sp_matrix_get(m,i,i);
    sp_matrix_scale_row(m,i,x);
    sp_matrix_scale_row(inv,i,x);
    /* For every row below us*/
    for(j = i+1;j<m->rows;j++){
      /* set leading element to 0 */
      x = -sp_matrix_get(m,j,i);
      sp_matrix_row_add_row(m,i,j,x);
      sp_matrix_row_add_row(inv,i,j,x);
    }    
  }

  
  /* Now from the bottom up */
  /* triangularize the matrix */
  /* For every row */
  for(i = m->rows-1;i>=0;i--){
    /* set leading element to 1 */
    x = 1.0/sp_matrix_get(m,i,i);
    sp_matrix_scale_row(m,i,x);
    sp_matrix_scale_row(inv,i,x);
    /* For every row above us*/
    for(j = i-1;j>=0;j--){
      /* set leading element to 0 */
      x = -sp_matrix_get(m,j,i);
      sp_matrix_row_add_row(m,i,j,x);
      sp_matrix_row_add_row(inv,i,j,x);
    }    
  }
  sp_matrix_memcpy(m,inv);
  sp_matrix_free(inv);
}

void sp_matrix_print(sp_matrix * a,FILE * fp){
  int i,j;
  if(!fp){
    fp = stdout;
  }
  for(i = 0;i<a->cols;i++){
    fprintf(fp,"|");
    for(j = 0;j<a->rows;j++){
      fprintf(fp,"\t%e",sp_matrix_get(a,i,j));
    }
    fprintf(fp,",\t|\n");
  }
}
