#include <stdlib.h>
#include <string.h>
#ifdef _USE_DMALLOC
#include <dmalloc.h>
#endif
#include "spimage.h"

sp_vector * _sp_vector_alloc(const int size,char * file, int line){
  sp_vector * ret = _sp_malloc(sizeof(sp_vector),file,line);
  ret->size = size;
  ret->data = _sp_calloc(size,sizeof(real),file,line);
  return ret;
}


sp_cvector * _sp_cvector_alloc(const int size,char * file, int line){
  sp_cvector * ret = _sp_malloc(sizeof(sp_cvector),file,line);
  ret->size = size;
  ret->data = _sp_calloc(size,sizeof(Complex),file,line);
  return ret;
}

void _sp_vector_free(sp_vector * v,char * file,int line){
  _sp_free(v->data,file,line);
  _sp_free(v,file,line);
}

void _sp_cvector_free(sp_cvector * v,char * file, int line){
  _sp_free(v->data,file,line);
  _sp_free(v,file,line);
}




sp_matrix * _sp_matrix_alloc(unsigned int nrows, unsigned int ncols,char * file, int line){
  sp_matrix * res = _sp_malloc(sizeof(sp_matrix),file,line);
  res->rows = nrows;
  res->cols = ncols;
  res->data = _sp_calloc(nrows*ncols,sizeof(real),file,line);
  return res;
}

sp_imatrix * _sp_imatrix_alloc(unsigned int nrows, unsigned int ncols,char * file,int line){
  sp_imatrix * res = _sp_malloc(sizeof(sp_imatrix),file,line);
  res->rows = nrows;
  res->cols = ncols;
  res->data = _sp_calloc(nrows*ncols,sizeof(int),file,line);
  return res;
}


sp_cmatrix * _sp_cmatrix_alloc(unsigned int nrows, unsigned int ncols,char * file,int line){
  sp_cmatrix * res = _sp_malloc(sizeof(sp_cmatrix),file,line);
  res->rows = nrows;
  res->cols = ncols;
  res->data = _sp_calloc(nrows*ncols,sizeof(Complex),file,line);
  return res;
}


sp_cmatrix * _sp_cmatrix_duplicate(sp_cmatrix * m,char * file, int line){
  sp_cmatrix * res = _sp_cmatrix_alloc(sp_cmatrix_rows(m),sp_cmatrix_cols(m),file,line);
  sp_cmatrix_memcpy(res,m);
  return res;
}

void _sp_matrix_free(sp_matrix * a,char * file, int line){
  _sp_free(a->data,file,line);
  _sp_free(a,file,line);
}


void _sp_imatrix_free(sp_imatrix * a,char * file, int line){
  _sp_free(a->data,file,line);
  _sp_free(a,file,line);
}


void _sp_cmatrix_free(sp_cmatrix * a,char * file, int line){
  _sp_free(a->data,file,line);
  _sp_free(a,file,line);
}



sp_3matrix * _sp_3matrix_alloc(unsigned int nx, unsigned int ny, unsigned int nz,char * file,int line){
  sp_3matrix * res = _sp_malloc(sizeof(sp_3matrix),file,line);
  res->x = nx;
  res->y = ny;
  res->z = nz;
  res->data = _sp_calloc(nx*ny*nz,sizeof(real),file,line);
  return res;
}

sp_i3matrix * _sp_i3matrix_alloc(unsigned int nx, unsigned int ny, unsigned int nz,char * file, int line){
  sp_i3matrix * res = _sp_malloc(sizeof(sp_i3matrix),file,line);
  res->x = nx;
  res->y = ny;
  res->z = nz;
  res->data = _sp_calloc(nx*ny*nz,sizeof(int),file,line);
  return res;
}

sp_c3matrix * _sp_c3matrix_alloc(unsigned int nx, unsigned int ny, unsigned int nz,char * file, int line){
  sp_c3matrix * res = _sp_malloc(sizeof(sp_c3matrix),file,line);
  res->x = nx;
  res->y = ny;
  res->z = nz;
  res->data = _sp_calloc(nx*ny*nz,sizeof(Complex),file,line);
  return res;
}


sp_c3matrix * _sp_c3matrix_duplicate(sp_c3matrix * m, char * file, int line){
  sp_c3matrix * res = _sp_c3matrix_alloc(sp_c3matrix_x(m),sp_c3matrix_y(m),sp_c3matrix_z(m),file,line);
  sp_c3matrix_memcpy(res,m);
  return res;
}

void _sp_3matrix_free(sp_3matrix * a,char * file, int line){
  _sp_free(a->data,file,line);
  _sp_free(a,file,line);
}

void _sp_i3matrix_free(sp_i3matrix * a,char * file, int line){
  _sp_free(a->data,file,line);
  _sp_free(a,file,line);
}

void _sp_c3matrix_free(sp_c3matrix * a,char * file, int line){
  _sp_free(a->data,file,line);
  _sp_free(a,file,line);
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
    assert(sp_real(sp_cmatrix_get(m,i,i)) != 0 && sp_imag(sp_cmatrix_get(m,i,i)) != 0);
    x = sp_cdiv(sp_cinit(1.0,0.0),sp_cmatrix_get(m,i,i));
    sp_cmatrix_scale_row(m,i,x);
    sp_cmatrix_scale_row(inv,i,x);
    /* For every row below us*/
    for(j = i+1;j<m->rows;j++){
      /* set leading element to 0 */
      x = sp_cmatrix_get(m,j,i);
      sp_real(x) = -sp_real(x);
      sp_imag(x) = -sp_imag(x);
      sp_cmatrix_row_add_row(m,i,j,x);
      sp_cmatrix_row_add_row(inv,i,j,x);
    }    
  }

  
  /* Now from the bottom up */
  /* triangularize the matrix */
  /* For every row */
  for(i = m->rows-1;i>=0;i--){
    /* set leading element to 1 */
    assert(sp_cabs(sp_cmatrix_get(m,i,i)));
    x = sp_cdiv(sp_cinit(1.0,0.0),sp_cmatrix_get(m,i,i));
    sp_cmatrix_scale_row(m,i,x);
    sp_cmatrix_scale_row(inv,i,x);
    /* For every row above us*/
    for(j = i-1;j>=0;j--){
      /* set leading element to 0 */
      x = sp_cmatrix_get(m,j,i);
      sp_real(x) = -sp_real(x);
      sp_imag(x) = -sp_imag(x);
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
      fprintf(fp,"\t%e %ei",sp_real(sp_cmatrix_get(a,i,j)),sp_imag(sp_cmatrix_get(a,i,j)));
    }
    fprintf(fp,",\t|\n");
  }
}


sp_vector * sp_c3matrix_center_of_mass(sp_c3matrix * a){
  sp_vector * res = sp_vector_alloc(3);
  int i = 0;
  real sum= 0;
  for(int z = 0;z<sp_c3matrix_z(a);z++){
    for(int y = 0;y<sp_c3matrix_y(a);y++){
      for(int x = 0;x<sp_c3matrix_x(a);x++){	
	real m = sp_cabs(a->data[i]);
	res->data[0] += m*x;
	res->data[1] += m*y;
	res->data[2] += m*z;
	sum += m;
	i++;
      }
    }
  }
  if(sum){
    res->data[0] /= sum;
    res->data[1] /= sum;
    res->data[2] /= sum;
  }
  return res;
}

int sp_complex_descend_compare(const void * pa,const void * pb){
  Complex a,b;
  a = *((Complex *)pa);
  b = *((Complex *)pb);
  if(sp_cabs(a) < sp_cabs(b)){
    return 1;
  }else if(sp_cabs(a) == sp_cabs(b)){
    return 0;
  }else{
    return -1;
  }
}

int sp_complex_ascend_compare(const void * pa,const void * pb){
  Complex a,b;
  a = *((Complex *)pa);
  b = *((Complex *)pb);
  if(sp_cabs(a) < sp_cabs(b)){
    return -1;
  }else if(sp_cabs(a) == sp_cabs(b)){
    return 0;
  }else{
    return 1;
  }
}

