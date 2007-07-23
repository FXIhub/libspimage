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



sp_3matrix * sp_3matrix_alloc(unsigned int nx, unsigned int ny, unsigned int nz){
  sp_3matrix * res = malloc(sizeof(sp_3matrix));
  res->x = nx;
  res->y = ny;
  res->z = nz;
  res->data = calloc(nx*ny*nz,sizeof(real));
  return res;
}

sp_i3matrix * sp_i3matrix_alloc(unsigned int nx, unsigned int ny, unsigned int nz){
  sp_i3matrix * res = malloc(sizeof(sp_i3matrix));
  res->x = nx;
  res->y = ny;
  res->z = nz;
  res->data = calloc(nx*ny*nz,sizeof(int));
  return res;
}

sp_c3matrix * sp_c3matrix_alloc(unsigned int nx, unsigned int ny, unsigned int nz){
  sp_c3matrix * res = malloc(sizeof(sp_c3matrix));
  res->x = nx;
  res->y = ny;
  res->z = nz;
  res->data = calloc(nx*ny*nz,sizeof(Complex));
  return res;
}


sp_c3matrix * sp_c3matrix_duplicate(sp_c3matrix * m){
  sp_c3matrix * res = sp_c3matrix_alloc(sp_c3matrix_x(m),sp_c3matrix_y(m),sp_c3matrix_z(m));
  sp_c3matrix_memcpy(res,m);
  return res;
}

void sp_3matrix_free(sp_3matrix * a){
  free(a->data);
  free(a);
}

void sp_i3matrix_free(sp_i3matrix * a){
  free(a->data);
  free(a);
}

void sp_c3matrix_free(sp_c3matrix * a){
  free(a->data);
  free(a);
}





void sp_matrix_invert(sp_matrix * m){
  int i,j;
  real x;
  sp_matrix * inv = sp_matrix_alloc(m->rows,m->cols);
  if(m->rows == 2 && m->cols == 2){
    double scale = sp_matrix_get(m,0,0)*sp_matrix_get(m,1,1)-sp_matrix_get(m,0,1)*sp_matrix_get(m,1,0);
    assert(scale != 0);
    sp_matrix_set(inv,0,0,sp_matrix_get(m,1,1)/scale);
    sp_matrix_set(inv,0,1,-sp_matrix_get(m,0,1)/scale);
    sp_matrix_set(inv,1,0,-sp_matrix_get(m,1,0)/scale);
    sp_matrix_set(inv,1,1,sp_matrix_get(m,0,0)/scale);
    sp_matrix_memcpy(m,inv);
    sp_matrix_free(inv);
    return;
  }
  sp_matrix_set_identity(inv);
  /* triangularize the matrix */
  /* For every row */
  for(i = 0;i<m->rows;i++){
    /* set leading element to 1 */
    assert(sp_matrix_get(m,i,i) != 0);
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
    assert(1.0/sp_matrix_get(m,i,i) != 0);
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

void sp_3matrix_invert(sp_3matrix * m){
  int i,j;
  real x;
  sp_3matrix * inv = sp_3matrix_alloc(m->x,m->y,0);
  if(m->y == 2 && m->x == 2){
    double scale = sp_3matrix_get(m,0,0,0)*sp_3matrix_get(m,1,1,0)-sp_3matrix_get(m,0,1,0)*sp_3matrix_get(m,1,0,0);
    assert(scale != 0);
    sp_3matrix_set(inv,0,0,0,sp_3matrix_get(m,1,1,0)/scale);
    sp_3matrix_set(inv,0,1,0,-sp_3matrix_get(m,0,1,0)/scale);
    sp_3matrix_set(inv,1,0,0,-sp_3matrix_get(m,1,0,0)/scale);
    sp_3matrix_set(inv,1,1,0,sp_3matrix_get(m,0,0,0)/scale);
    sp_3matrix_memcpy(m,inv);
    sp_3matrix_free(inv);
    return;
  }
  sp_3matrix_set_identity(inv);
  /* triangularize the matrix */
  /* For every row */
  for(i = 0;i<m->y;i++){
    /* set leading element to 1 */
    assert(sp_3matrix_get(m,i,i,0) != 0);
    x = 1.0/sp_3matrix_get(m,i,i,0);
    sp_3matrix_scale_row(m,i,x);
    sp_3matrix_scale_row(inv,i,x);
    /* For every row below us*/
    for(j = i+1;j<m->y;j++){
      /* set leading element to 0 */
      x = -sp_3matrix_get(m,j,i,0);
      sp_3matrix_row_add_row(m,i,j,x);
      sp_3matrix_row_add_row(inv,i,j,x);
    }    
  }

  
  /* Now from the bottom up */
  /* triangularize the matrix */
  /* For every row */
  for(i = m->y-1;i>=0;i--){
    /* set leading element to 1 */
    assert(1.0/sp_3matrix_get(m,i,i,0) != 0);
    x = 1.0/sp_3matrix_get(m,i,i,0);
    sp_3matrix_scale_row(m,i,x);
    sp_3matrix_scale_row(inv,i,x);
    /* For every row above us*/
    for(j = i-1;j>=0;j--){
      /* set leading element to 0 */
      x = -sp_3matrix_get(m,j,i,0);
      sp_3matrix_row_add_row(m,i,j,x);
      sp_3matrix_row_add_row(inv,i,j,x);
    }    
  }
  sp_3matrix_memcpy(m,inv);
  sp_3matrix_free(inv);
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
    assert(1.0/sp_cmatrix_get(m,i,i) != 0);
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
    assert(1.0/sp_cmatrix_get(m,i,i) != 0);
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
