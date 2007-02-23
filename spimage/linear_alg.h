#ifndef _LINEAR_ALG_H_
#define _LINEAR_ALG_H_ 1

#include <complex.h>
#include <string.h>
#include <math.h>

typedef struct{
  unsigned int size;
  real * data;
} sp_vector;

typedef struct{
  unsigned int size;
  complex * data;
} sp_cvector;


typedef struct{
  unsigned int rows;
  unsigned int cols;
  real * data;
} sp_matrix;


static inline real sp_min(real a,real b){
  if(a< b){
    return a;
  }
  return b;
}

static inline real sp_max(real a,real b){
  if(a > b){
    return a;
  }
  return b;
}


/*! This function allocates memory for a matrix of size nrows rows by ncols columns and initializes all the elements of the matrix to zero.
 *
 */
spimage_EXPORT sp_matrix * sp_matrix_alloc(unsigned int nrows, unsigned int ncols);

/*! This function frees a previously allocated matrix m
 *
 */
spimage_EXPORT void sp_matrix_free(sp_matrix * m);


/*! Creates an empty zero initialized vector of the desired size.
 *
 * This function creates a vector of length n, returning a pointer
 * to a newly initialized vector struct. All vector elements are set
 * to 0.
 */
spimage_EXPORT sp_vector * sp_vector_alloc(const int size);

/*! Creates an empty zero initialized complex vector of the desired size.
 *
 * This function creates a vector of length n, returning a pointer
 * to a newly initialized vector struct. All vector elements are set
 * to 0.
 */
spimage_EXPORT sp_cvector * sp_cvector_alloc(const int size);

/*! Frees a previously allocated vector.
 *
 */
spimage_EXPORT void sp_vector_free(sp_vector * v);

/*! Frees a previously allocated complex vector.
 *
 */
spimage_EXPORT void sp_cvector_free(sp_cvector * v);


/*! This function returns the size of the vector v
 *
 */
static inline unsigned int sp_vector_size(const sp_vector * v){
  return v->size;   
}

/*! This function returns the size of the complex vector v
 *
 */
static inline unsigned int sp_cvector_size(const sp_cvector * v){
  return v->size;   
}

/*! This function returns the (row,col)-th element of a matrix m.
 *
 * row and col must lie in the range of 0 to nrows-1 and 0 to ncols-1.
 */
static inline real sp_matrix_get (const sp_matrix * m, unsigned int row, unsigned int col){
  return m->data[row*m->cols+col];
}

/*! This function sets the (row,col)-th element of a matrix m to x.
 *
 * row and col must lie in the range of 0 to nrows-1 and 0 to ncols-1.
 */
static inline void sp_matrix_set (sp_matrix * m, unsigned int row, unsigned int col, real x){
  m->data[row*m->cols+col] = x;
}


/*! This function returns the nth element of a vector v.
 *
 * n must be in the range of 0 to size-1.
 */
static inline real sp_vector_get (const sp_vector * v, unsigned int n){
  return v->data[n];
}

/*! This function returns the the nth element of a complex vector v.
 *
 * n must be in the range of 0 to size-1.
 */
static inline complex sp_cvector_get(const sp_cvector * v, unsigned int n){
  return v->data[n];
}

/*! This function sets the nth element of a vector v to x.
 *
 * n must be in the range of 0 to size-1.
 */
static inline void sp_vector_set (const sp_vector * v, unsigned int n, real x){
  v->data[n] = x;
}

/*! This function sets the nth element of a complex vector v to x.
 *
 * n must be in the range of 0 to size-1.
 */
static inline void sp_cvector_set (const sp_cvector * v, unsigned int n, complex x){
  v->data[n] = x;
}

/*! This function copies the elements of the matrix b into the matrix a. 
  *
  * The two matrices must have the same dimensions.
  */
static inline void sp_matrix_memcpy(sp_matrix * dest, const sp_matrix * src){
  int i;
  if(src->rows*src->cols < 1024){
    /* avoid function call and make inline possibly useful */
    for(i = 0;i<src->rows*src->cols;i++){
      dest->data[i] = src->data[i];
    }
  }else{
    memcpy(dest->data,src->data,sizeof(real)*src->rows*src->cols);
  }
}


/*! This function adds the elements of vector b to the elements of vector a, a'_i = a_i + b_i. 
 *
 * The two vectors must have the same length.
 */
static inline void sp_vector_add(sp_vector * a, const sp_vector * b){
  int i;
  for(i = 0;i<a->size;i++){
    a->data[i] += b->data[i];
  }
}


/*! This function adds the elements of a complex vector b to the elements of a complex vector a, a'_i = a_i + b_i. 
 *
 * The two vectors must have the same length.
 */
static inline void sp_cvector_add(sp_cvector * a, const sp_cvector * b){
  int i;
  for(i = 0;i<a->size;i++){
    a->data[i] += b->data[i];
  }
}

/*! This function subtracts the elements of vector b to the elements of vector a, a'_i = a_i - b_i. 
 *
 * The two vectors must have the same length.
 */
static inline void sp_vector_sub(sp_vector * a, const sp_vector * b){
  int i;
  for(i = 0;i<a->size;i++){
    a->data[i] -= b->data[i];
  }
}


/*! This function subtracts the elements of complex vector b to the elements of complex vector a, a'_i = a_i - b_i. 
 *
 * The two vectors must have the same length.
 */
static inline void sp_cvector_sub(sp_cvector * a, const sp_cvector * b){
  int i;
  for(i = 0;i<a->size;i++){
    a->data[i] -= b->data[i];
  }
}

/*! This function mutiplies the elements of vector a with the elements of vector b, a'_i = a_i * b_i. 
 *
 * The two vectors must have the same length.
 */
static inline void sp_vector_mul(sp_vector * a, const sp_vector * b){
  int i;
  for(i = 0;i<a->size;i++){
    a->data[i] *= b->data[i];
  }
}

/*! This function mutiplies the elements of complex vector a with the elements of complex vector b, a'_i = a_i * b_i. 
 *
 * The two vectors must have the same length.
 */
static inline void sp_cvector_mul(sp_cvector * a, const sp_cvector * b){
  int i;
  for(i = 0;i<a->size;i++){
    a->data[i] *= b->data[i];
  }
}

/*! This function divides the elements of vector a with the elements of vector b, a'_i = a_i / b_i. 
 *
 * The two vectors must have the same length.
 */
static inline void sp_vector_div(sp_vector * a, const sp_vector * b){
  int i;
  for(i = 0;i<a->size;i++){
    a->data[i] /= b->data[i];
  }
}

/*! This function divides the elements of the complex vector a with the elements of the complex vector b, a'_i = a_i / b_i. 
 *
 * The two vectors must have the same length.
 */
static inline void sp_cvector_div(sp_cvector * a, const sp_cvector * b){
  int i;
  for(i = 0;i<a->size;i++){
    a->data[i] /= b->data[i];
  }
}


/*! This function multiplies the elements of vector a by the constant factor x, a'_i = x a_i.
 *
 */
static inline void sp_vector_scale(sp_vector * a, const real x){
  int i;
  for(i = 0;i<a->size;i++){
    a->data[i] *= x;
  }
}

/*! This function multiplies the elements of a complex vector a by the constant factor x, a'_i = x a_i.
 *
 */
static inline void sp_cvector_scale(sp_cvector * a, const complex x){
  int i;
  for(i = 0;i<a->size;i++){
    a->data[i] *= x;
  }
}


/*! This function adds the constant value x to the elements of the vector a, a'_i = a_i + x.
 *
 */
static inline void sp_vector_add_constant(sp_vector * a, const real x){
  int i;
  for(i = 0;i<a->size;i++){
    a->data[i] += x;
  }
}

/*! This function adds the constant value x to the elements of the complex vector a, a'_i = a_i + x.
 *
 */
static inline void sp_cvector_add_constant(sp_cvector * a, const complex x){
  int i;
  for(i = 0;i<a->size;i++){
    a->data[i] += x;
  }
}


/*! This function calculates the dot product of a and b.
  *
  * The two vectors must have the same length.
  */
static inline real sp_vector_dot_prod(sp_vector * a, const sp_vector * b){
  int i;
  real ret = 0;
  for(i = 0;i<a->size;i++){
    ret += a->data[i]*b->data[i];
  }
  return ret;
}

/*! This function calculates the dot product of a and b.
  *
  * The two vectors must have the same length.
  */
static inline real sp_cvector_dot_prod(sp_cvector * a, const sp_cvector * b){
  int i;
  real ret = 0;
  for(i = 0;i<a->size;i++){
    ret += a->data[i]*b->data[i];
  }
  return ret;
}


/*! This function calculates the outer product of vector a and b.
  *
  */
static inline sp_matrix * sp_vector_outer_prod(sp_vector * a, const sp_vector * b){
  int i,j;
  sp_matrix * ret = sp_matrix_alloc(a->size,b->size);
  for(i = 0;i<a->size;i++){
    for(j = 0;j<b->size;j++){
      sp_matrix_set(ret,i,j,a->data[i]*b->data[j]);
    }
  }
  return ret;
}


/*! This function calculates the norm of the vector a.
  *
  */
static inline real sp_vector_norm(sp_vector * a){
  int i;
  real ret = 0;
  for(i = 0;i<a->size;i++){
    ret += a->data[i]*a->data[i];
  }
  return sqrt(ret);
}

/*! This function calculates the norm of the vector a.
  *
  */
static inline real sp_cvector_norm(sp_cvector * a){
  int i;
  real ret = 0;
  for(i = 0;i<a->size;i++){
    ret += a->data[i]*a->data[i];
  }
  return sqrt(ret);
}


/*! This function copies the elements of the vector src into the vector dest. 
  *
  * The two vectors must have the same length.
  */
static inline void sp_vector_memcpy(sp_vector * dest, const sp_vector * src){
  int i;
  if(src->size < 1024){
    for(i = 0;i<src->size;i++){
      dest->data[i] = src->data[i];
    }
  }else{
    memcpy(dest->data,src->data,sizeof(real)*src->size);
  }

}

/*! This function copies the elements of the complex vector src into the complex vector dest. 
  *
  * The two vectors must have the same length.
  */
static inline void sp_cvector_memcpy(sp_cvector * dest, const sp_cvector * src){
  int i;
  if(src->size < 1024){
    for(i = 0;i<src->size;i++){
      dest->data[i] = src->data[i];
    }
  }else{
    memcpy(dest->data,src->data,sizeof(real)*src->size);
  }
}






/*! This function returns the number of cells in m, that is, rows*cols 
 *
 */
static inline unsigned int sp_matrix_size (const sp_matrix * m){
  return m->rows*m->cols;
}


/*! This function returns the number of rows in m
 *
 */
static inline unsigned int sp_matrix_rows (const sp_matrix * m){
  return m->rows;
}


/*! This function returns the number of colums in m
 *
 */
static inline unsigned int sp_matrix_cols (const sp_matrix * m){
  return m->cols;
}


/*! This function sets the diagonal elements of m to 1 and the rest to 0
 *
 */
static inline void sp_matrix_set_identity(sp_matrix * m){
  int i;
  memset(m->data,0,sizeof(real)*m->rows*m->cols);
  for(i = 0;i<m->rows && i<m->cols;i++){
    sp_matrix_set(m,i,i,1);
  }
}




/*! This function adds the elements of matrix b to the elements of matrix a, a'_ij = a_ij + b_ij. 
 *
 * The two matrices must have the same dimensions.
 */
static inline void sp_matrix_add(sp_matrix * a, const sp_matrix * b){
  int i;
  for(i = 0;i<sp_matrix_size(b);i++){
    a->data[i] += b->data[i];
  }
}

/*! This function subtracts the elements of matrix b to the elements of matrix a, a'_ij = a_ij - b_ij. 
 *
 * The two matrices must have the same dimensions.
 */
static inline void sp_matrix_sub(sp_matrix * a, const sp_matrix * b){
  int i;
  for(i = 0;i<sp_matrix_size(b);i++){
    a->data[i] -= b->data[i];
  }
}

/*! This function mutiplies the elements of matrix a with the elements of matrix b, a'_ij = a_ij * b_ij. 
 *
 * The two matrices must have the same dimensions.
 */
static inline void sp_matrix_mul_elements(sp_matrix * a, const sp_matrix * b){
  int i;
  for(i = 0;i<sp_matrix_size(b);i++){
    a->data[i] *= b->data[i];
  }
}


/*! This function divides the elements of matrix a with the elements of matrix b, a'_ij = a_ij / b_ij. 
 *
 * The two matrices must have the same dimensions.
 */
static inline void sp_matrix_div_elements(sp_matrix * a, const sp_matrix * b){
  int i;
  for(i = 0;i<sp_matrix_size(b);i++){
    a->data[i] /= b->data[i];
  }
}


/*! This function multiplies the elements of the matrix a by the constant factor x, a'_ij = x a_ij.
 *
 */
static inline void sp_matrix_scale(sp_matrix * a, const real x){
  int i;
  for(i = 0;i<sp_matrix_size(a);i++){
    a->data[i] *= x;
  }
}


/*! This function adds the constant value x to the elements of the matrix a, a'_ij = a_ij + x.
 *
 */
static inline void sp_matrix_add_constant(sp_matrix * a, const real x){
  int i;
  for(i = 0;i<sp_matrix_size(a);i++){
    a->data[i] += x;
  }
}

/*! This function transposes matrix a.
 *
 */
static inline void sp_matrix_transpose(sp_matrix * a){
  int i,j;
  /* exchange dimensions */
  sp_matrix * tmp = sp_matrix_alloc(a->cols,a->rows);
  for(i = 0;i<a->rows;i++){
    for(j = 0;j<a->cols;j++){
      sp_matrix_set(tmp,j,i,sp_matrix_get(a,i,j));
    }
  }
  /* copy from tmp the useful things and discard original array */
  a->cols = tmp->cols;
  a->rows = tmp->rows;
  free(a->data);
  a->data = tmp->data;
  free(tmp);
}


/*! This function returns the product of matrix m with vector v.
 *
 *  The size of v must be the same as the number of cols in m.
 */
static inline sp_vector * sp_matrix_vector_prod(const sp_matrix * m, const sp_vector * v){
  int i,j;
  sp_vector * ret = sp_vector_alloc(m->rows);
  for(i = 0;i<m->rows;i++){
    for(j = 0;j<m->cols;j++){
      ret->data[i] += sp_matrix_get(m,i,j)*v->data[j];
    }
  }
  return ret;
}


/*! This function returns the product of matrix a with matrix b.
 *
 *  The size of a must be the same as the size of the transpose of b.
 */
static inline sp_matrix * sp_matrix_mul(const sp_matrix * a, const sp_matrix * b){
  int i,j,k;
  real tmp;
  sp_matrix * ret = sp_matrix_alloc(a->rows,b->cols);
  for(i = 0;i<a->rows;i++){
    for(j = 0;j<b->cols;j++){
      tmp = 0;
      for(k = 0;k<a->cols;k++){
	tmp += sp_matrix_get(a,i,k)*sp_matrix_get(b,k,j);
      }
      sp_matrix_set(ret,i,j,tmp);
    }
  }
  return ret;
}

/*! This function scales all elements of row n of matrix m by x.
 *
 */
static inline void sp_matrix_scale_row(sp_matrix * m, int n, real x){
  int i;
  for(i = 0;i<sp_matrix_cols(m);i++){
    sp_matrix_set(m,n,i,sp_matrix_get(m,n,i)*x);
  }
}

/*! This function add the elements of row from, multiplied by factor, to the elements of row to, of matrix m.
 *
 */
static inline void sp_matrix_row_add_row(sp_matrix * m, int from, int to, real factor){
  int i;
  for(i = 0;i<sp_matrix_cols(m);i++){
    sp_matrix_set(m,to,i,sp_matrix_get(m,from,i)*factor+sp_matrix_get(m,to,i));
  }
}

/*! This functions returns the inverse of matrix a.
 *
 */
spimage_EXPORT void sp_matrix_invert(sp_matrix * a);

/*! This functions tries to print the matrix in a human readable way
 */
spimage_EXPORT void sp_matrix_print(sp_matrix * a, FILE * fp);



#endif
