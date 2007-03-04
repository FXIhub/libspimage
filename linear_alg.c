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
  ret->data = calloc(size,sizeof(Complex));
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

sp_imatrix * sp_imatrix_alloc(unsigned int nrows, unsigned int ncols){
  sp_imatrix * res = malloc(sizeof(sp_imatrix));
  res->rows = nrows;
  res->cols = ncols;
  res->data = calloc(nrows*ncols,sizeof(int));
  return res;
}


sp_cmatrix * sp_cmatrix_alloc(unsigned int nrows, unsigned int ncols){
  sp_cmatrix * res = malloc(sizeof(sp_cmatrix));
  res->rows = nrows;
  res->cols = ncols;
  res->data = calloc(nrows*ncols,sizeof(Complex));
  return res;
}


sp_cmatrix * sp_cmatrix_duplicate(sp_cmatrix * m){
  sp_cmatrix * res = sp_cmatrix_alloc(sp_cmatrix_rows(m),sp_cmatrix_cols(m));
  sp_cmatrix_memcpy(res,m);
  return res;
}

void sp_matrix_free(sp_matrix * a){
  free(a->data);
  free(a);
}


void sp_imatrix_free(sp_imatrix * a){
  free(a->data);
  free(a);
}


void sp_cmatrix_free(sp_cmatrix * a){
  free(a->data);
  free(a);
}






void sp_matrix_invert(sp_matrix * m){
  int i,j;
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


void sp_cmatrix_invert(sp_cmatrix * m){
  int i,j;
  Complex x;
  sp_cmatrix * inv = sp_cmatrix_alloc(m->rows,m->cols);
  sp_cmatrix_set_identity(inv);
  /* triangularize the matrix */
  /* For every row */
  for(i = 0;i<m->rows;i++){
    /* set leading element to 1 */
    x = 1.0/sp_cmatrix_get(m,i,i);
    sp_cmatrix_scale_row(m,i,x);
    sp_cmatrix_scale_row(inv,i,x);
    /* For every row below us*/
    for(j = i+1;j<m->rows;j++){
      /* set leading element to 0 */
      x = -sp_cmatrix_get(m,j,i);
      sp_cmatrix_row_add_row(m,i,j,x);
      sp_cmatrix_row_add_row(inv,i,j,x);
    }    
  }

  
  /* Now from the bottom up */
  /* triangularize the matrix */
  /* For every row */
  for(i = m->rows-1;i>=0;i--){
    /* set leading element to 1 */
    x = 1.0/sp_cmatrix_get(m,i,i);
    sp_cmatrix_scale_row(m,i,x);
    sp_cmatrix_scale_row(inv,i,x);
    /* For every row above us*/
    for(j = i-1;j>=0;j--){
      /* set leading element to 0 */
      x = -sp_cmatrix_get(m,j,i);
      sp_cmatrix_row_add_row(m,i,j,x);
      sp_cmatrix_row_add_row(inv,i,j,x);
    }    
  }
  sp_cmatrix_memcpy(m,inv);
  sp_cmatrix_free(inv);
}

void sp_matrix_print(sp_matrix * a,FILE * fp){
  int i,j;
  if(!fp){
    fp = stdout;
  }
  for(i = 0;i<a->rows;i++){
    fprintf(fp,"|");
    for(j = 0;j<a->cols;j++){
      fprintf(fp,"\t%e",sp_matrix_get(a,i,j));
    }
    fprintf(fp,",\t|\n");
  }
}


void sp_imatrix_print(sp_imatrix * a,FILE * fp){
  int i,j;
  if(!fp){
    fp = stdout;
  }
  for(i = 0;i<a->rows;i++){
    fprintf(fp,"|");
    for(j = 0;j<a->cols;j++){
      fprintf(fp,"\t%d",sp_imatrix_get(a,i,j));
    }
    fprintf(fp,",\t|\n");
  }
}

void sp_cmatrix_print(sp_cmatrix * a,FILE * fp){
  int i,j;
  if(!fp){
    fp = stdout;
  }
  for(i = 0;i<a->rows;i++){
    fprintf(fp,"|");
    for(j = 0;j<a->cols;j++){
      fprintf(fp,"\t%e %ei",creal(sp_cmatrix_get(a,i,j)),cimag(sp_cmatrix_get(a,i,j)));
    }
    fprintf(fp,",\t|\n");
  }
}
