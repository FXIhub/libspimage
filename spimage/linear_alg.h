#ifndef _LINEAR_ALG_H_
#define _LINEAR_ALG_H_ 1

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "mem_util.h"


#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */


#ifdef _SP_DOUBLE_PRECISION
typedef double real;
#define sp_real(a) (a).re
#define sp_imag(a) (a).im
typedef struct {
     double re, im;
}Complex;
#define REAL_MAX DBL_MAX
#define REAL_EPSILON DBL_EPSILON
#else
typedef float real;
#define sp_real(a) (a).re
#define sp_imag(a) (a).im
typedef struct {
     float re, im;
}Complex;
#define REAL_MAX FLT_MAX
#define REAL_EPSILON FLT_EPSILON
#endif


typedef struct{
  unsigned int size;
  real * data;
} sp_vector;

typedef struct{
  unsigned int size;
  Complex * data;
} sp_cvector;


typedef struct{
  unsigned int rows;
  unsigned int cols;
  real * data;
} sp_matrix;

typedef struct{
  unsigned int rows;
  unsigned int cols;
  int * data;
} sp_imatrix;


typedef struct{
  unsigned int rows;
  unsigned int cols;
  Complex * data;
} sp_cmatrix;

typedef struct{
  unsigned int x;
  unsigned int y;
  unsigned int z;
  real * data;
} sp_3matrix;

typedef struct{
  unsigned int x;
  unsigned int y;
  unsigned int z;
  int * data;
} sp_i3matrix;

typedef struct{
  unsigned int x;
  unsigned int y;
  unsigned int z;
  Complex * data;
} sp_c3matrix;


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

#define sp_swap(a,b,t){ t _temp = a; a = b;b = _temp}


/* You have to be careful when using this macro because it can cause unintended side effects! */
#define sp_cincr(a,b){sp_real(a) += sp_real(b); sp_imag(a) += sp_imag(b);}

static inline Complex sp_cinit(real a, real b){
  Complex ret = {a,b};
  return ret;
}

static inline Complex sp_cconj(Complex a){  
  sp_imag(a) = -sp_imag(a);
  return a;
}

static inline real sp_cabs(Complex a){
  return sqrt(sp_imag(a)*sp_imag(a)+sp_real(a)*sp_real(a));
}


static inline real sp_carg(Complex a){
    return atan2(sp_imag(a),sp_real(a));
}

static inline Complex sp_cadd(Complex a, Complex b){
  sp_real(a) += sp_real(b);
  sp_imag(a) += sp_imag(b);
  return a;
}


static inline Complex sp_csub(Complex a, Complex b){
  sp_real(a) -= sp_real(b);
  sp_imag(a) -= sp_imag(b);
  return a;
}


static inline Complex sp_cscale(Complex a, real b){
  sp_real(a) *= b;
  sp_imag(a) *= b;
  return a;
}

static inline Complex sp_cmul(Complex a, Complex b){
  Complex ret;
  sp_real(ret) = sp_real(a)*sp_real(b)-sp_imag(a)*sp_imag(b);
  sp_imag(ret) = sp_real(a)*sp_imag(b)+sp_real(b)*sp_imag(a);  
  return ret;
}

static inline Complex sp_cdiv(Complex a, Complex b){
  Complex ret;
  sp_real(ret) = (sp_real(a)*sp_real(b)+sp_imag(a)*sp_imag(b))/(sp_real(b)*sp_real(b)+sp_imag(b)*sp_imag(b));
  sp_imag(ret) = (sp_imag(a)*sp_real(b)-sp_real(a)*sp_imag(b))/(sp_real(b)*sp_real(b)+sp_imag(b)*sp_imag(b));
  return ret;
}

/*! This function allocates memory for a 3matrix of size nrows rows by ncols columns and initializes all the elements of the matrix to zero.
 *
 */
spimage_EXPORT sp_3matrix * _sp_3matrix_alloc(unsigned int nx, unsigned int ny, unsigned int nz,char * file, int line);
#define sp_3matrix_alloc(nx,ny,nz) _sp_3matrix_alloc(nx,ny,nz,__FILE__,__LINE__)

/*! This function allocates memory for a i3matrix of size nrows rows by ncols columns and initializes all the elements of the matrix to zero.
 *
 */
spimage_EXPORT sp_i3matrix * _sp_i3matrix_alloc(unsigned int nx, unsigned int ny, unsigned int nz,char * file,int line);
#define sp_i3matrix_alloc(nx,ny,nz) _sp_i3matrix_alloc(nx,ny,nz,__FILE__,__LINE__)

/*! This function allocates memory for a Complex 3matrix of size nrows rows by ncols columns and initializes all the elements of the matrix to zero.
 *
 */
spimage_EXPORT sp_c3matrix * _sp_c3matrix_alloc(unsigned int nx, unsigned int ny, unsigned int nz, char * file, int line);
#define sp_c3matrix_alloc(nx,ny,nz) _sp_c3matrix_alloc(nx,ny,nz,__FILE__,__LINE__)

/*! This function creates a duplicate of it's argument and returns a pointer to it
 *
 */
spimage_EXPORT sp_c3matrix * _sp_c3matrix_duplicate(sp_c3matrix * m, char * file, int line);
#define sp_c3matrix_duplicate(m) _sp_c3matrix_duplicate(m,__FILE__,__LINE__)

/*! This function frees a previously allocated matrix m
 *
 */
spimage_EXPORT void _sp_3matrix_free(sp_3matrix * m, char * file, int line);
#define sp_3matrix_free(m) _sp_3matrix_free(m,__FILE__,__LINE__)

/*! This function frees a previously allocated matrix m
 *
 */
spimage_EXPORT void _sp_i3matrix_free(sp_i3matrix * m, char * file, int line);
#define sp_i3matrix_free(m) _sp_i3matrix_free(m,__FILE__,__LINE__)

/*! This function frees a previously allocated Complex matrix m
 *
 */
spimage_EXPORT void _sp_c3matrix_free(sp_c3matrix * m, char * file, int line);
#define sp_c3matrix_free(m) _sp_c3matrix_free(m,__FILE__,__LINE__)



/*! This function allocates memory for a matrix of size nrows rows by ncols columns and initializes all the elements of the matrix to zero.
 *
 */
spimage_EXPORT sp_matrix * _sp_matrix_alloc(unsigned int nrows, unsigned int ncols, char * file, int line);
#define sp_matrix_alloc(nrows,ncols) _sp_matrix_alloc(nrows,ncols,__FILE__,__LINE__)

/*! This function allocates memory for a Complex matrix of size nrows rows by ncols columns and initializes all the elements of the matrix to zero.
 *
 */
spimage_EXPORT sp_cmatrix * _sp_cmatrix_alloc(unsigned int nrows, unsigned int ncols, char * file, int line);
#define sp_cmatrix_alloc(nrows,ncols) _sp_cmatrix_alloc(nrows,ncols,__FILE__,__LINE__)

/*! This function creates a duplicate of it's argument and returns a pointer to it
 *
 */
spimage_EXPORT sp_cmatrix * _sp_cmatrix_duplicate(sp_cmatrix * m, char * file, int line);
#define sp_cmatrix_duplicate(m) _sp_cmatrix_duplicate(m,__FILE__,__LINE__)

/*! This function allocates memory for an Integer matrix of size nrows rows by ncols columns and initializes all the elements of the matrix to zero.
 *
 */
spimage_EXPORT sp_imatrix * _sp_imatrix_alloc(unsigned int nrows, unsigned int ncols, char * file, int line);
#define sp_imatrix_alloc(nrows,ncols) _sp_imatrix_alloc(nrows,ncols,__FILE__,__LINE__)

/*! This function frees a previously allocated matrix m
 *
 */
spimage_EXPORT void _sp_matrix_free(sp_matrix * m, char * file, int line);
#define sp_matrix_free(m) _sp_matrix_free(m,__FILE__,__LINE__)

/*! This function frees a previously allocated Complex matrix m
 *
 */
spimage_EXPORT void _sp_cmatrix_free(sp_cmatrix * m, char * file, int line);
#define sp_cmatrix_free(m) _sp_cmatrix_free(m,__FILE__,__LINE__)

/*! This function frees a previously allocated Complex matrix m
 *
 */
spimage_EXPORT void _sp_imatrix_free(sp_imatrix * m, char * file, int line);
#define sp_imatrix_free(m) _sp_imatrix_free(m,__FILE__,__LINE__)


/*! Creates an empty zero initialized vector of the desired size.
 *
 * This function creates a vector of length n, returning a pointer
 * to a newly initialized vector struct. All vector elements are set
 * to 0.
 */
spimage_EXPORT sp_vector * _sp_vector_alloc(const int size, char * file, int line);
#define sp_vector_alloc(size) _sp_vector_alloc(size,__FILE__,__LINE__)

/*! Creates an empty zero initialized Complex vector of the desired size.
 *
 * This function creates a vector of length n, returning a pointer
 * to a newly initialized vector struct. All vector elements are set
 * to 0.
 */
spimage_EXPORT sp_cvector * _sp_cvector_alloc(const int size, char * file, int line);
#define sp_cvector_alloc(size) _sp_cvector_alloc(size,__FILE__,__LINE__)
 
/*! Frees a previously allocated vector.
 *
 */
spimage_EXPORT void _sp_vector_free(sp_vector * v, char * file, int line);
#define sp_vector_free(v) _sp_vector_free(v,__FILE__,__LINE__)

/*! Frees a previously allocated Complex vector.
 *
 */
spimage_EXPORT void _sp_cvector_free(sp_cvector * v, char * file, int line);
#define sp_cvector_free(v) _sp_cvector_free(v,__FILE__,__LINE__)

/*! This function returns the size of the vector v
 *
 */
static inline unsigned int sp_vector_size(const sp_vector * v){
  return v->size;   
}

/*! This function returns the size of the Complex vector v
 *
 */
static inline unsigned int sp_cvector_size(const sp_cvector * v){
  return v->size;   
}



/*! This function returns the (x,y,z)-th element of a matrix m.
 *
 * x, y and z must lie in the range of 0 to x-1, 0 to y-1 and 0 to z-1.
 */
static inline real sp_3matrix_get (const sp_3matrix * m, unsigned int x, unsigned int y, unsigned int z){
  return m->data[z*m->y*m->x+y*m->x+x];
}

/*! This function returns the (x,y,z)-th element of a matrix m.
 *
 * x, y and z must lie in the range of 0 to x-1, 0 to y-1 and 0 to z-1.
 */
static inline int sp_i3matrix_get (const sp_i3matrix * m, unsigned int x, unsigned int y, unsigned int z){
  return m->data[z*m->y*m->x+y*m->x+x];
}

/*! This function returns the (x,y,z)-th element of a Complex matrix m.
 *
 * x, y and z must lie in the range of 0 to x-1, 0 to y-1 and 0 to z-1.
 */
static inline Complex sp_c3matrix_get (const sp_c3matrix * m, unsigned int x, unsigned int y, unsigned int z){
  return m->data[z*m->y*m->x+y*m->x+x];
}

/*! This function sets the (x,y,z)-th element of a matrix m to entry.
 *
 * x, y and z must lie in the range of 0 to x-1, 0 to y-1 and 0 to z-1.
 */
static inline void sp_3matrix_set (sp_3matrix * m, unsigned int x, unsigned int y, unsigned int z, real entry){
  m->data[z*m->y*m->x+y*m->x+x] = entry;
}

/*! This function sets the (x,y,z)-th element of a matrix m to entry.
 *
 * x, y and z must lie in the range of 0 to x-1, 0 to y-1 and 0 to z-1.
 */
static inline void sp_i3matrix_set (sp_i3matrix * m, unsigned int x, unsigned int y, unsigned int z, int entry){
  m->data[z*m->y*m->x+y*m->x+x] = entry;
}

/*! This function sets the (x,y,z)-th element of an Complex matrix m to entry.
 *
 * x, y and z must lie in the range of 0 to x-1, 0 to y-1 and 0 to z-1.
 */
static inline void sp_c3matrix_set (sp_c3matrix * m, unsigned int x, unsigned int y, unsigned int z, Complex entry){
  m->data[z*m->y*m->x+y*m->x+x] = entry;
}

/*! This function returns the (row,col)-th element of a matrix m.
 *
 * Row and col must lie in the range of 0 to nrows-1 and 0 to ncols-1.
 */
static inline real sp_matrix_get (const sp_matrix * m, unsigned int row, unsigned int col){
  return m->data[col*m->rows+row];
}

/*! This function returns the (row,col)-th element of an Integer matrix m.
 *
 * row and col must lie in the range of 0 to nrows-1 and 0 to ncols-1.
 */
static inline int sp_imatrix_get (const sp_imatrix * m, unsigned int row, unsigned int col){
  return m->data[col*m->rows+row];
}


/*! This function returns the (row,col)-th element of a Complex matrix m.
 *
 * row and col must lie in the range of 0 to nrows-1 and 0 to ncols-1.
 */
static inline Complex sp_cmatrix_get (const sp_cmatrix * m, unsigned int row, unsigned int col){
  return m->data[col*m->rows+row];
}

/*! This function sets the (row,col)-th element of a matrix m to x.
 *
 * row and col must lie in the range of 0 to nrows-1 and 0 to ncols-1.
 */
static inline void sp_matrix_set (sp_matrix * m, unsigned int row, unsigned int col, real x){
  m->data[col*m->rows+row] = x;
}

/*! This function sets the (row,col)-th element of an Integer matrix m to x.
 *
 * row and col must lie in the range of 0 to nrows-1 and 0 to ncols-1.
 */
static inline void sp_imatrix_set (sp_imatrix * m, unsigned int row, unsigned int col, int n){
  m->data[col*m->rows+row] = n;
}

/*! This function sets the (row,col)-th element of a Complex matrix m to x.
 *
 * row and col must lie in the range of 0 to nrows-1 and 0 to ncols-1.
 */
static inline void sp_cmatrix_set (sp_cmatrix * m, unsigned int row, unsigned int col, Complex x){
  m->data[col*m->rows+row] = x;
}

/*! This function returns the nth element of a vector v.
 *
 * n must be in the range of 0 to size-1.
 */
static inline real sp_vector_get (const sp_vector * v, unsigned int n){
  return v->data[n];
}

/*! This function returns the the nth element of a Complex vector v.
 *
 * n must be in the range of 0 to size-1.
 */
static inline Complex sp_cvector_get(const sp_cvector * v, unsigned int n){
  return v->data[n];
}

/*! This function sets the nth element of a vector v to x.
 *
 * n must be in the range of 0 to size-1.
 */
static inline void sp_vector_set (const sp_vector * v, unsigned int n, real x){
  v->data[n] = x;
}

/*! This function sets the nth element of a Complex vector v to x.
 *
 * n must be in the range of 0 to size-1.
 */
static inline void sp_cvector_set (const sp_cvector * v, unsigned int n, Complex x){
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


/*! This function copies the elements of the integer matrix b into the integer matrix a. 
  *
  * The two matrices must have the same dimensions.
  */
static inline void sp_imatrix_memcpy(sp_imatrix * dest, const sp_imatrix * src){
  int i;
  if(src->rows*src->cols < 1024){
    /* avoid function call and make inline possibly useful */
    for(i = 0;i<src->rows*src->cols;i++){
      dest->data[i] = src->data[i];
    }
  }else{
    memcpy(dest->data,src->data,sizeof(int)*src->rows*src->cols);
  }
}


/*! This function copies the elements of the Complex matrix b into the Complex matrix a. 
  *
  * The two matrices must have the same dimensions.
  */
static inline void sp_cmatrix_memcpy(sp_cmatrix * dest, const sp_cmatrix * src){
  int i;
  if(src->rows*src->cols < 1024){
    /* avoid function call and make inline possibly useful */
    for(i = 0;i<src->rows*src->cols;i++){
      dest->data[i] = src->data[i];
    }
  }else{
    memcpy(dest->data,src->data,sizeof(Complex)*src->rows*src->cols);
  }
}

/*! This function copies the elements of the matrix b into the matrix a. 
  *
  * The two matrices must have the same dimensions.
  */
static inline void sp_3matrix_memcpy(sp_3matrix * dest, const sp_3matrix * src){
  int i;
  if(src->x*src->y*src->z < 1024){
    /* avoid function call and make inline possibly useful */
    for(i = 0;i<src->x*src->y*src->z;i++){
      dest->data[i] = src->data[i];
    }
  }else{
    memcpy(dest->data,src->data,sizeof(real)*src->x*src->y*src->z);
  }
}

/*! This function copies the elements of the matrix b into the matrix a. 
  *
  * The two matrices must have the same dimensions.
  */
static inline void sp_i3matrix_memcpy(sp_i3matrix * dest, const sp_i3matrix * src){
  int i;
  if(src->x*src->y*src->z < 1024){
    /* avoid function call and make inline possibly useful */
    for(i = 0;i<src->x*src->y*src->z;i++){
      dest->data[i] = src->data[i];
    }
  }else{
    memcpy(dest->data,src->data,sizeof(int)*src->x*src->y*src->z);
  }
}

/*! This function copies the elements of the Complex matrix b into the Complex matrix a. 
  *
  * The two matrices must have the same dimensions.
  */
static inline void sp_c3matrix_memcpy(sp_c3matrix * dest, const sp_c3matrix * src){
  int i;
  if(src->x*src->y*src->z < 1024){
    /* avoid function call and make inline possibly useful */
    for(i = 0;i<src->x*src->y*src->z;i++){
      dest->data[i] = src->data[i];
    }
  }else{
    memcpy(dest->data,src->data,sizeof(Complex)*src->x*src->y*src->z);
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


/*! This function adds the elements of a Complex vector b to the elements of a Complex vector a, a'_i = a_i + b_i. 
 *
 * The two vectors must have the same length.
 */
static inline void sp_cvector_add(sp_cvector * a, const sp_cvector * b){
  int i;
  for(i = 0;i<a->size;i++){
    a->data[i] = sp_cadd(a->data[i],b->data[i]);
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


/*! This function subtracts the elements of Complex vector b to the elements of Complex vector a, a'_i = a_i - b_i. 
 *
 * The two vectors must have the same length.
 */
static inline void sp_cvector_sub(sp_cvector * a, const sp_cvector * b){
  int i;
  for(i = 0;i<a->size;i++){
    a->data[i] = sp_csub(a->data[i],b->data[i]);
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

/*! This function mutiplies the elements of Complex vector a with the elements of Complex vector b, a'_i = a_i * b_i. 
 *
 * The two vectors must have the same length.
 */
static inline void sp_cvector_mul(sp_cvector * a, const sp_cvector * b){
  int i;
  for(i = 0;i<a->size;i++){
    a->data[i] = sp_cmul(a->data[i],b->data[i]);
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

/*! This function divides the elements of the Complex vector a with the elements of the Complex vector b, a'_i = a_i / b_i. 
 *
 * The two vectors must have the same length.
 */
static inline void sp_cvector_div(sp_cvector * a, const sp_cvector * b){
  int i;
  for(i = 0;i<a->size;i++){
    a->data[i] = sp_cdiv(a->data[i],b->data[i]);
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

/*! This function multiplies the elements of a Complex vector a by the constant factor x, a'_i = x a_i.
 *
 */
static inline void sp_cvector_scale(sp_cvector * a, const Complex x){
  int i;
  for(i = 0;i<a->size;i++){
    a->data[i] = sp_cmul(a->data[i],x);
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

/*! This function adds the constant value x to the elements of the Complex vector a, a'_i = a_i + x.
 *
 */
static inline void sp_cvector_add_constant(sp_cvector * a, const Complex x){
  int i;
  for(i = 0;i<a->size;i++){
    a->data[i] = sp_cadd(a->data[i],x);
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
static inline Complex sp_cvector_dot_prod(sp_cvector * a, const sp_cvector * b){
  int i;
  Complex ret = {0,0};
  for(i = 0;i<a->size;i++){
    ret = sp_cadd(ret,sp_cmul(a->data[i],sp_cconj(b->data[i])));
  }
  return ret;
}


/*! This function calculates the outer product of vector a and b.
  *
  */
static inline sp_3matrix * sp_vector_outer_prod(sp_vector * a, const sp_vector * b){
  int i,j;
  sp_3matrix * ret = sp_3matrix_alloc(a->size,b->size,1);
  for(i = 0;i<a->size;i++){
    for(j = 0;j<b->size;j++){
      sp_3matrix_set(ret,i,j,0,a->data[i]*b->data[j]);
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
    ret += sp_cabs(a->data[i])*sp_cabs(a->data[i]);
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

/*! This function copies the elements of the Complex vector src into the Complex vector dest. 
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
static inline int sp_matrix_size (const sp_matrix * m){
  return m->rows*m->cols;
}


/*! This function returns the number of cells in m, that is, rows*cols 
 *
 */
static inline int sp_imatrix_size (const sp_imatrix * m){
  return m->rows*m->cols;
}

/*! This function returns the number of cells in m, that is, rows*cols 
 *
 */
static inline int sp_cmatrix_size (const sp_cmatrix * m){
  return m->rows*m->cols;
}


/*! This function returns the number of rows in m
 *
 */
static inline int sp_matrix_rows (const sp_matrix * m){
  return m->rows;
}

/*! This function returns the number of rows in m
 *
 */
static inline int sp_imatrix_rows (const sp_imatrix * m){
  return m->rows;
}

/*! This function returns the number of rows in m
 *
 */
static inline int sp_cmatrix_rows (const sp_cmatrix * m){
  return m->rows;
}


/*! This function returns the number of colums in m
 *
 */
static inline int sp_matrix_cols (const sp_matrix * m){
  return m->cols;
}

/*! This function returns the number of colums in m
 *
 */
static inline int sp_imatrix_cols (const sp_imatrix * m){
  return m->cols;
}

/*! This function returns the number of colums in m
 *
 */
static inline int sp_cmatrix_cols (const sp_cmatrix * m){
  return m->cols;
}

/*! This function returns the size in x of m
 *
 */
static inline int sp_3matrix_x (const sp_3matrix * m){
  return m->x;
}

/*! This function returns the size in y of m
 *
 */
static inline int sp_3matrix_y (const sp_3matrix * m){
  return m->y;
}

/*! This function returns the size in z of m
 *
 */
static inline int sp_3matrix_z (const sp_3matrix * m){
  return m->z;
}

/*! This function returns the size in x of m
 *
 */
static inline int sp_i3matrix_x (const sp_i3matrix * m){
  return m->x;
}

/*! This function returns the size in y of m
 *
 */
static inline int sp_i3matrix_y (const sp_i3matrix * m){
  return m->y;
}

/*! This function returns the size in z of m
 *
 */
static inline int sp_i3matrix_z (const sp_i3matrix * m){
  return m->z;
}

/*! This function returns the size in x of m
 *
 */
static inline int sp_c3matrix_x (const sp_c3matrix * m){
  return m->x;
}

/*! This function returns the size in y of m
 *
 */
static inline int sp_c3matrix_y (const sp_c3matrix * m){
  return m->y;
}

/*! This function returns the size in z of m
 *
 */
static inline int sp_c3matrix_z (const sp_c3matrix * m){
  return m->z;
}

/*! This function returns the number of cells in m, that is, x*y*z 
 *
 */
static inline long long sp_3matrix_size (const sp_3matrix * m){
  long long s;
  s = m->x*m->y*m->z;
  return s;
}

/*! This function returns the number of cells in m, that is, x*y*z 
 *
 */
static inline long long sp_i3matrix_size (const sp_i3matrix * m){
  return m->x*m->y*m->z;
}

/*! This function returns the number of cells in m, that is, x*y*z 
 *
 */
static inline long long sp_c3matrix_size (const sp_c3matrix * m){
  return m->x*m->y*m->z;
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

static inline void sp_3matrix_set_identity(sp_3matrix * m){
  int i;
  memset(m->data,0,sizeof(real)*m->x*m->y);
  for(i = 0;i<m->x && i<m->y;i++){
    sp_3matrix_set(m,i,i,0,1);
  }
}

/*! This function sets the diagonal elements of m to 1 and the rest to 0
 *
 */
static inline void sp_imatrix_set_identity(sp_imatrix * m){
  int i;
  memset(m->data,0,sizeof(int)*m->rows*m->cols);
  for(i = 0;i<m->rows && i<m->cols;i++){
    sp_imatrix_set(m,i,i,1);
  }
}

/*! This function sets the diagonal elements of m to 1 and the rest to 0
 *
 */
static inline void sp_cmatrix_set_identity(sp_cmatrix * m){
  int i;
  Complex one = {1,0};
  memset(m->data,0,sizeof(Complex)*m->rows*m->cols);
  for(i = 0;i<m->rows && i<m->cols;i++){
    sp_cmatrix_set(m,i,i,one);
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

/*! This function adds the elements of integer matrix b to the elements of integer matrix a, a'_ij = a_ij + b_ij. 
 *
 * The two matrices must have the same dimensions.
 */
static inline void sp_imatrix_add(sp_imatrix * a, const sp_imatrix * b){
  int i;
  for(i = 0;i<sp_imatrix_size(b);i++){
    a->data[i] += b->data[i];
  }
}

/*! This function adds the elements of Complex matrix b scaled by x to the elements of Complex matrix a
 *  a'_ij = a_ij + b_ij * x. 
 *
 * The two matrices must have the same dimensions.
 */
static inline void sp_cmatrix_add(sp_cmatrix * a, const sp_cmatrix * b, Complex * x){
  int i;
  if(x && ((sp_real(*x) != 1) || (sp_imag(*x) != 0))){
    for(i = 0;i<sp_cmatrix_size(b);i++){
      a->data[i] = sp_cadd(a->data[i],sp_cmul(b->data[i],(*x)));
    }
  }else{
    for(i = 0;i<sp_cmatrix_size(b);i++){
      a->data[i] = sp_cadd(a->data[i],b->data[i]);
    }
  }
}

/*! This function adds the elements of matrix b to the elements of matrix a, a'_ijk = a_ijk + b_ijk. 
 *
 * The two matrices must have the same dimensions.
 */
static inline void sp_3matrix_add(sp_3matrix * a, const sp_3matrix * b){
  int i;
  for(i = 0;i<sp_3matrix_size(b);i++){
    a->data[i] += b->data[i];
  }
}

/*! This function adds the elements of integer matrix b to the elements of integer matrix a, a'_ijk = a_ijk + b_ijk. 
 *
 * The two matrices must have the same dimensions.
 */
static inline void sp_i3matrix_add(sp_i3matrix * a, const sp_i3matrix * b){
  int i;
  for(i = 0;i<sp_i3matrix_size(b);i++){
    a->data[i] += b->data[i];
  }
}

/*! This function adds the elements of Complex matrix b scaled by x to the elements of Complex matrix a
 *  a'_ijk = a_ijk + b_ijk * x. 
 *
 * The two matrices must have the same dimensions.
 */
static inline void sp_c3matrix_add(sp_c3matrix * a, const sp_c3matrix * b, Complex * x){
  int i;
  if(x && ((sp_real(*x) != 1) || (sp_imag(*x) != 0))){
    for(i = 0;i<sp_c3matrix_size(b);i++){
      a->data[i] = sp_cadd(a->data[i],sp_cmul(b->data[i],(*x)));
    }
  }else{
    for(i = 0;i<sp_c3matrix_size(b);i++){
      a->data[i] = sp_cadd(a->data[i],b->data[i]);
    }
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


/*! This function subtracts the elements of integer matrix b to the elements of integer matrix a, a'_ij = a_ij - b_ij. 
 *
 * The two matrices must have the same dimensions.
 */
static inline void sp_imatrix_sub(sp_imatrix * a, const sp_imatrix * b){
  int i;
  for(i = 0;i<sp_imatrix_size(b);i++){
    a->data[i] -= b->data[i];
  }
}

/*! This function subtracts the elements of Complex matrix b to the elements of Complex matrix a, a'_ij = a_ij - b_ij. 
 *
 * The two matrices must have the same dimensions.
 */
static inline void sp_cmatrix_sub(sp_cmatrix * a, const sp_cmatrix * b){
  int i;
  for(i = 0;i<sp_cmatrix_size(b);i++){
    a->data[i] = sp_csub(a->data[i],b->data[i]);
  }
}


/*! This function subtracts the elements of matrix b to the elements of matrix a, a'_ij = a_ij - b_ij. 
 *
 * The two matrices must have the same dimensions.
 */
static inline void sp_3matrix_sub(sp_3matrix * a, const sp_3matrix * b){
  int i;
  for(i = 0;i<sp_3matrix_size(b);i++){
    a->data[i] -= b->data[i];
  }
}


/*! This function subtracts the elements of integer matrix b to the elements of integer matrix a, a'_ij = a_ij - b_ij. 
 *
 * The two matrices must have the same dimensions.
 */
static inline void sp_i3matrix_sub(sp_i3matrix * a, const sp_i3matrix * b){
  int i;
  for(i = 0;i<sp_i3matrix_size(b);i++){
    a->data[i] -= b->data[i];
  }
}

/*! This function subtracts the elements of Complex matrix b to the elements of Complex matrix a, a'_ij = a_ij - b_ij. 
 *
 * The two matrices must have the same dimensions.
 */
static inline void sp_c3matrix_sub(sp_c3matrix * a, const sp_c3matrix * b){
  int i;
  for(i = 0;i<sp_c3matrix_size(b);i++){
    a->data[i] = sp_csub(a->data[i],b->data[i]);
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

/*! This function mutiplies the elements of integer matrix a with the elements of integer matrix b, a'_ij = a_ij * b_ij. 
 *
 * The two matrices must have the same dimensions.
 */
static inline void sp_imatrix_mul_elements(sp_imatrix * a, const sp_imatrix * b){
  int i;
  for(i = 0;i<sp_imatrix_size(b);i++){
    a->data[i] *= b->data[i];
  }
}

/*! This function mutiplies the elements of Complex matrix a with the elements of Complex matrix b, a'_ij = a_ij * b_ij. 
 *
 * The two matrices must have the same dimensions.
 */
static inline void sp_cmatrix_mul_elements(sp_cmatrix * a, const sp_cmatrix * b){
  int i;
  for(i = 0;i<sp_cmatrix_size(b);i++){
    a->data[i] = sp_cmul(a->data[i],b->data[i]);
  }
}

/*! This function mutiplies the elements of matrix a with the elements of matrix b, a'_ijk = a_ijk * b_ijk. 
 *
 * The two matrices must have the same dimensions.
 */
static inline void sp_3matrix_mul_elements(sp_3matrix * a, const sp_3matrix * b){
  int i;
  for(i = 0;i<sp_3matrix_size(b);i++){
    a->data[i] *= b->data[i];
  }
}

/*! This function mutiplies the elements of integer matrix a with the elements of integer matrix b, a'_ijk = a_ijk * b_ijk. 
 *
 * The two matrices must have the same dimensions.
 */
static inline void sp_i3matrix_mul_elements(sp_i3matrix * a, const sp_i3matrix * b){
  int i;
  for(i = 0;i<sp_i3matrix_size(b);i++){
    a->data[i] *= b->data[i];
  }
}

/*! This function mutiplies the elements of Complex matrix a with the elements of Complex matrix b, a'_ijk = a_ijk * b_ijk. 
 *
 * The two matrices must have the same dimensions.
 */
static inline void sp_c3matrix_mul_elements(sp_c3matrix * a, const sp_c3matrix * b){
  int i;
  for(i = 0;i<sp_c3matrix_size(b);i++){
    a->data[i] = sp_cmul(a->data[i],b->data[i]);
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

static inline void sp_3matrix_div_elements(sp_3matrix * a, const sp_3matrix * b){
  int i;
  for(i = 0;i<sp_3matrix_size(b);i++){
    a->data[i] /= b->data[i];
  }
}

/*! This function divides the elements of integer matrix a with the elements of integer matrix b, a'_ij = a_ij / b_ij. 
 *
 * The two matrices must have the same dimensions.
 */
static inline void sp_imatrix_div_elements(sp_imatrix * a, const sp_imatrix * b){
  int i;
  for(i = 0;i<sp_imatrix_size(b);i++){
    a->data[i] /= b->data[i];
  }
}

static inline void sp_i3matrix_div_elements(sp_i3matrix * a, const sp_i3matrix * b){
  int i;
  for(i = 0;i<sp_i3matrix_size(b);i++){
    a->data[i] /= b->data[i];
  }
}

/*! This function divides the elements of Complex matrix a with the elements of Complex matrix b, a'_ij = a_ij / b_ij. 
 *
 * The two matrices must have the same dimensions.
 */
static inline void sp_cmatrix_div_elements(sp_cmatrix * a, const sp_cmatrix * b){
  int i;
  for(i = 0;i<sp_cmatrix_size(b);i++){
    assert(sp_cabs(b->data[i]) != 0);
    a->data[i] = sp_cdiv(a->data[i],b->data[i]);
  }
}

static inline void sp_c3matrix_div_elements(sp_c3matrix * a, const sp_c3matrix * b){
  int i;
  for(i = 0;i<sp_c3matrix_size(b);i++){
    assert(sp_cabs(b->data[i]) != 0);
    a->data[i] = sp_cdiv(a->data[i],b->data[i]);
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

/*! This function multiplies the elements of the integer matrix a by the constant factor x, a'_ij = x a_ij.
 *
 */
static inline void sp_imatrix_scale(sp_imatrix * a, const real x){
  int i;
  for(i = 0;i<sp_imatrix_size(a);i++){
    a->data[i] *= x;
  }
}


/*! This function multiplies the elements of the Complex matrix a by the constant factor x, a'_ij = x a_ij.
 *
 */
static inline void sp_cmatrix_scale(sp_cmatrix * a, const Complex x){
  int i;
  for(i = 0;i<sp_cmatrix_size(a);i++){
    a->data[i] = sp_cmul(a->data[i],x);
  }
}

/*! This function multiplies the elements of the matrix a by the constant factor x, a'_ij = x a_ij.
 *
 */
static inline void sp_3matrix_scale(sp_3matrix * a, const real x){
  int i;
  for(i = 0;i<sp_3matrix_size(a);i++){
    a->data[i] *= x;
  }
}

/*! This function multiplies the elements of the integer matrix a by the constant factor x, a'_ij = x a_ij.
 *
 */
static inline void sp_i3matrix_scale(sp_i3matrix * a, const real x){
  int i;
  for(i = 0;i<sp_i3matrix_size(a);i++){
    a->data[i] *= x;
  }
}


/*! This function multiplies the elements of the Complex matrix a by the constant factor x, a'_ij = x a_ij.
 *
 */
static inline void sp_c3matrix_scale(sp_c3matrix * a, const Complex x){
  int i;
  for(i = 0;i<sp_c3matrix_size(a);i++){
    a->data[i] = sp_cmul(a->data[i],x);
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
}static inline void sp_3matrix_add_constant(sp_3matrix * a, const real x){
  int i;
  for(i = 0;i<sp_3matrix_size(a);i++){
    a->data[i] += x;
  }
}

/*! This function adds the constant value x to the elements of the integer matrix a, a'_ij = a_ij + x.
 *
 */
static inline void sp_imatrix_add_constant(sp_imatrix * a, const int x){
  int i;
  for(i = 0;i<sp_imatrix_size(a);i++){
    a->data[i] += x;
  }
}
static inline void sp_i3matrix_add_constant(sp_i3matrix * a, const int x){
  int i;
  for(i = 0;i<sp_i3matrix_size(a);i++){
    a->data[i] += x;
  }
}


/*! This function adds the constant value x to the elements of the Complex matrix a, a'_ij = a_ij + x.
 *
 */
static inline void sp_cmatrix_add_constant(sp_cmatrix * a, const Complex x){
  int i;
  for(i = 0;i<sp_cmatrix_size(a);i++){
    a->data[i] = sp_cadd(a->data[i],x);
  }
}
static inline void sp_c3matrix_add_constant(sp_c3matrix * a, const Complex x){
  int i;
  for(i = 0;i<sp_c3matrix_size(a);i++){
    a->data[i] = sp_cadd(a->data[i],x);
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
  sp_free(a->data);
  a->data = tmp->data;
  sp_free(tmp);
}

static inline void sp_3matrix_transpose(sp_3matrix * a){
  int i,j,k;
  /* exchange dimensions */
  sp_3matrix * tmp = sp_3matrix_alloc(a->x,a->y,a->z);
  for(i = 0;i<a->x;i++){
    for(j = 0;j<a->y;j++){
      for(k = 0;k<a->z;k++){
      sp_3matrix_set(tmp,j,i,k,sp_3matrix_get(a,i,j,k));
      }
    }
  }
  /* copy from tmp the useful things and discard original array */
  a->x = tmp->x;
  a->y = tmp->y;
  a->z = tmp->z;
  sp_free(a->data);
  a->data = tmp->data;
  sp_free(tmp);
}

/*! This function transposes matrix a.
 *
 */
static inline void sp_imatrix_transpose(sp_imatrix * a){
  int i,j;
  /* exchange dimensions */
  sp_imatrix * tmp = sp_imatrix_alloc(a->cols,a->rows);
  for(i = 0;i<a->rows;i++){
    for(j = 0;j<a->cols;j++){
      sp_imatrix_set(tmp,j,i,sp_imatrix_get(a,i,j));
    }
  }
  /* copy from tmp the useful things and discard original array */
  a->cols = tmp->cols;
  a->rows = tmp->rows;
  sp_free(a->data);
  a->data = tmp->data;
  sp_free(tmp);
}

static inline void sp_i3matrix_transpose(sp_i3matrix * a){
  int i,j,k;
  /* exchange dimensions */
  sp_i3matrix * tmp = sp_i3matrix_alloc(a->x,a->y,a->z);
  for(i = 0;i<a->x;i++){
    for(j = 0;j<a->y;j++){
      for(k = 0;k<a->z;k++){
	sp_i3matrix_set(tmp,i,j,k,sp_i3matrix_get(a,j,i,k));
      }
    }
  }
  /* copy from tmp the useful things and discard original array */
  a->x = tmp->x;
  a->y = tmp->y;
  a->z = tmp->z;
  sp_free(a->data);
  a->data = tmp->data;
  sp_free(tmp);
}


/*! This function transposes Complex matrix a.
 *
 * The Complex conjugate of the cells is not calculated
 */
static inline void sp_cmatrix_transpose(sp_cmatrix * a){
  int i,j;
  /* exchange dimensions */
  sp_cmatrix * tmp = sp_cmatrix_alloc(a->cols,a->rows);
  for(i = 0;i<a->rows;i++){
    for(j = 0;j<a->cols;j++){
      sp_cmatrix_set(tmp,j,i,sp_cmatrix_get(a,i,j));
    }
  }
  /* copy from tmp the useful things and discard original array */
  a->cols = tmp->cols;
  a->rows = tmp->rows;
  sp_free(a->data);
  a->data = tmp->data;
  sp_free(tmp);
}

static inline void sp_c3matrix_transpose_xy(sp_c3matrix * a){
  int i,j,k;
  /* exchange dimensions */
  sp_c3matrix * tmp = sp_c3matrix_alloc(a->x,a->y,a->z);
  for(i = 0;i<a->x;i++){
    for(j = 0;j<a->y;j++){
      for(k = 0;k<a->z;k++){
	sp_c3matrix_set(tmp,i,j,k,sp_c3matrix_get(a,j,i,k));
      }
    }
  }
  /* copy from tmp the useful things and discard original array */
  a->x = tmp->x;
  a->y = tmp->y;
  a->z = tmp->z;
  sp_free(a->data);
  a->data = tmp->data;
  sp_free(tmp);
}

static inline void sp_i3matrix_transpose_xy(sp_i3matrix * a){
  int i,j,k;
  /* exchange dimensions */
  sp_i3matrix * tmp = sp_i3matrix_alloc(a->x,a->y,a->z);
  for(i = 0;i<a->x;i++){
    for(j = 0;j<a->y;j++){
      for(k = 0;k<a->z;k++){
	sp_i3matrix_set(tmp,i,j,k,sp_i3matrix_get(a,j,i,k));
      }
    }
  }
  /* copy from tmp the useful things and discard original array */
  a->x = tmp->x;
  a->y = tmp->y;
  a->z = tmp->z;
  sp_free(a->data);
  a->data = tmp->data;
  sp_free(tmp);
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

static inline sp_vector * sp_3matrix_vector_prod(const sp_3matrix * m, const sp_vector * v){
  int i,j;
  sp_vector * ret = sp_vector_alloc(m->y);
  for(i = 0;i<m->y;i++){
    for(j = 0;j<m->x;j++){
      ret->data[i] += sp_3matrix_get(m,j,i,0)*v->data[j];
    }
  }
  return ret;
}

/*! This function returns the product of Complex matrix m with Complex vector v.
 *
 *  The size of v must be the same as the number of cols in m.
 */
static inline sp_cvector * sp_cmatrix_cvector_prod(const sp_cmatrix * m, const sp_cvector * v){
  int i,j;
  sp_cvector * ret = sp_cvector_alloc(m->rows);
  for(i = 0;i<m->rows;i++){
    for(j = 0;j<m->cols;j++){
      ret->data[i] = sp_cadd(ret->data[i],sp_cmul(sp_cmatrix_get(m,i,j),v->data[j]));
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

static inline sp_3matrix * sp_3matrix_mul(const sp_3matrix * a, const sp_3matrix * b){
  int i,j,k;
  real tmp;
  sp_3matrix * ret = sp_3matrix_alloc(a->x,b->y,0);
  for(i = 0;i<a->x;i++){
    for(j = 0;j<b->y;j++){
      tmp = 0;
      for(k = 0;k<a->y;k++){
	tmp += sp_3matrix_get(a,k,i,0)*sp_3matrix_get(b,j,k,0);
      }
      sp_3matrix_set(ret,i,j,0,tmp);
    }
  }
  return ret;
}

/*! This function returns the product of Complex matrix a with Complex matrix b.
 *
 *  The size of a must be the same as the size of the transpose of b.
 */
static inline sp_cmatrix * sp_cmatrix_mul(const sp_cmatrix * a, const sp_cmatrix * b){
  int i,j,k;
  Complex tmp;
  sp_cmatrix * ret = sp_cmatrix_alloc(a->rows,b->cols);
  for(i = 0;i<a->rows;i++){
    for(j = 0;j<b->cols;j++){
      sp_real(tmp) = 0;
      sp_imag(tmp) = 0;
      for(k = 0;k<a->cols;k++){
	tmp = sp_cadd(tmp,sp_cmul(sp_cmatrix_get(a,i,k),sp_cmatrix_get(b,k,j)));
      }
      sp_cmatrix_set(ret,i,j,tmp);
    }
  }
  return ret;
}

static inline sp_c3matrix * sp_c3matrix_mul(const sp_c3matrix * a, const sp_c3matrix * b){
  int i,j,k;
  Complex tmp;
  sp_c3matrix * ret = sp_c3matrix_alloc(a->x,b->y,0);
  for(i = 0;i<a->x;i++){
    for(j = 0;j<b->y;j++){
      sp_real(tmp) = 0;
      sp_imag(tmp) = 0;
      for(k = 0;k<a->y;k++){
	tmp = sp_cadd(tmp,sp_cmul(sp_c3matrix_get(a,k,i,0),sp_c3matrix_get(b,j,k,0)));
      }
      sp_c3matrix_set(ret,i,j,0,tmp);
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

static inline void sp_3matrix_scale_row(sp_3matrix * m, int n, real x){
  int i;
  for(i = 0;i<sp_3matrix_x(m);i++){
    sp_3matrix_set(m,n,i,0,sp_3matrix_get(m,n,i,0)*x);
  }
}

/*! This function scales all elements of row n of Complex matrix m by x.
 *
 */
static inline void sp_cmatrix_scale_row(sp_cmatrix * m, int n, Complex x){
  int i;
  for(i = 0;i<sp_cmatrix_cols(m);i++){
    sp_cmatrix_set(m,n,i,sp_cmul(sp_cmatrix_get(m,n,i),x));
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

static inline void sp_3matrix_row_add_row(sp_3matrix * m, int from, int to, real factor){
  int i;
  for(i = 0;i<sp_3matrix_x(m);i++){
    sp_3matrix_set(m,to,i,0,sp_3matrix_get(m,from,i,0)*factor+sp_3matrix_get(m,to,i,0));
  }
}

/*! This function add the elements of row from, multiplied by factor, to the elements of row to, of Complex matrix m.
 *
 */
static inline void sp_cmatrix_row_add_row(sp_cmatrix * m, int from, int to, Complex factor){
  int i;
  for(i = 0;i<sp_cmatrix_cols(m);i++){
    sp_cmatrix_set(m,to,i,sp_cadd(sp_cmul(sp_cmatrix_get(m,from,i),factor),sp_cmatrix_get(m,to,i)));
  }
}

/*! This functions returns the inverse of matrix a.
 *
 */
spimage_EXPORT void sp_matrix_invert(sp_matrix * a);

spimage_EXPORT void sp_3matrix_invert(sp_3matrix * a);

/*! This functions tries to print the matrix in a human readable way
 */
spimage_EXPORT void sp_matrix_print(sp_matrix * a, FILE * fp);


/*! This functions returns the inverse of Complex matrix a.
 *
 */
spimage_EXPORT void sp_cmatrix_invert(sp_cmatrix * a);

/*! This functions tries to print the Complex matrix in a human readable way
 */
spimage_EXPORT void sp_cmatrix_print(sp_cmatrix * a, FILE * fp);


/*! This function returns a cvector numerically similar to v
 *
 * The elements of v are copied to the real part of the result
 */
static inline sp_cvector * sp_vector_to_cvector(const sp_vector * v){
  int i;
  int size = sp_vector_size(v);
  sp_cvector * ret = sp_cvector_alloc(size);
  for(i = 0;i<size;i++){
    sp_real(ret->data[i]) = v->data[i];
    sp_imag(ret->data[i]) = 0;
  }
  return ret;
}

/*! This function returns a cmatrix numerically similar to m
 *
 * The elements of m are copied to the real part of the result
 */
static inline sp_cmatrix * sp_matrix_to_cmatrix(const sp_matrix * m){
  int i,j;
  int rows = sp_matrix_rows(m);
  int cols = sp_matrix_cols(m);
  sp_cmatrix * ret = sp_cmatrix_alloc(rows,cols);
  for(i = 0;i<rows;i++){
    for(j = 0;j<cols;j++){
      Complex tmp = {sp_matrix_get(m,i,j),0};
      sp_cmatrix_set(ret,i,j,tmp);
    }
  }
  return ret;
}



/*! This function returns the index of a given row and column combination
 *
 */
static inline int sp_matrix_get_index(const sp_matrix * m, const int row, const int col){
  return col*m->rows+row;
}


/*! This function returns the index of a given row and column combination
 *
 */
static inline int sp_imatrix_get_index(const sp_imatrix * m, const int row, const int col){
  return col*m->rows+row;
}

/*! This function returns the index of a given row and column combination
 *
 */
static inline int sp_cmatrix_get_index(const sp_cmatrix * m, const int row, const int col){
  return col*m->rows+row;
}

/*! This function returns the index of a given (x,y,z) combination
 *
 */
static inline long long sp_3matrix_get_index(const sp_3matrix * m, const int x, const int y, const int z){
  return z*m->y*m->x+y*m->x+x;
}


/*! This function returns the index of a given (x,y,z) combination
 *
 */
static inline long long sp_i3matrix_get_index(const sp_i3matrix * m, const int x, const int y, const int z){
  return z*m->y*m->x+y*m->x+x;
}

/*! This function returns the index of a given row and column combination
 *
 */
static inline long long sp_c3matrix_get_index(const sp_c3matrix * m, const int x, const int y, const int z){
  return z*m->y*m->x+y*m->x+x;
}

/*! This function changes m to it's complex conjugate
 *
 */
static inline void sp_cmatrix_conj(sp_cmatrix * m){
  int i;
  for(i = 0;i<sp_cmatrix_size(m);i++){
    m->data[i] = sp_cconj(m->data[i]);
  }
}

static inline void sp_c3matrix_conj(sp_c3matrix * m){
  int i;
  for(i = 0;i<sp_c3matrix_size(m);i++){
    m->data[i] = sp_cconj(m->data[i]);
  }
}

/*! Returns the cabs value of the element with the smallest cabs
 *
 */
static inline real sp_cmatrix_min(const sp_cmatrix * m, int * index){
  real min = sp_cabs(m->data[0]);
  int i,ii;
  for(i = 1;i<sp_cmatrix_size(m);i++){
    if(sp_cabs(m->data[i]) < min){
      min = sp_cabs(m->data[i]);
      ii = i;
    }
  }
  if(index){
    *index = ii;
  }
  return min;
}

static inline real sp_c3matrix_min(const sp_c3matrix * m, long long * index){
  real min = sp_cabs(m->data[0]);
  int i,ii;
  for(i = 1;i<sp_c3matrix_size(m);i++){
    if(sp_cabs(m->data[i]) < min){
      min = sp_cabs(m->data[i]);
      ii = i;
    }
  }
  if(index){
    *index = ii;
  }
  return min;
}

/*! Returns the cabs value of the element with the biggest cabs
 *
 */
static inline real sp_cmatrix_max(const sp_cmatrix * m, int * index){
  real max = sp_cabs(m->data[0]);
  int i;
  int i_max = 0;
  for(i = 1;i<sp_cmatrix_size(m);i++){
    if(sp_cabs(m->data[i]) > max){
      max = sp_cabs(m->data[i]);
      i_max = i;
    }
  }
  if(index){
    *index = i_max;
  }
  return max;
}
static inline real sp_c3matrix_max(const sp_c3matrix * m, long long * index){
  real max = sp_cabs(m->data[0]);
  long long i;
  long long i_max = 0;
  for(i = 1;i<sp_c3matrix_size(m);i++){
    if(sp_cabs(m->data[i]) > max){
      max = sp_cabs(m->data[i]);
      i_max = i;
    }
  }
  if(index){
    *index = i_max;
  }
  return max;
}

/*! Returns the interpolated value of m at the floating point row frow and column fcol
 *
 */
static inline Complex sp_cmatrix_interp(const sp_cmatrix * m, real frow, real fcol){
  int x = fcol;
  int y = frow;
  real u = fcol-x;
  real v = frow-y;
  Complex res = {0,0};
  if(x >= sp_cmatrix_cols(m)-1){
    x = sp_cmatrix_cols(m)-2;
    u = 1;
  }
  if(y >= sp_cmatrix_rows(m)-1){
    y = sp_cmatrix_rows(m)-2;
    v = 1;
  }
  res = sp_cadd(sp_cadd(sp_cadd(sp_cscale(sp_cscale(sp_cmatrix_get(m,y,x),(1-u)),(1-v)),
			sp_cscale(sp_cscale(sp_cmatrix_get(m,y,x+1),(u)),(1-v))),
			sp_cscale(sp_cscale(sp_cmatrix_get(m,y+1,x),(1-u)),(v))),
		sp_cscale(sp_cscale(sp_cmatrix_get(m,y+1,x+1),(u)),(v)));
  return res;
}

static inline Complex sp_c3matrix_interp(const sp_c3matrix * m, real fx, real fy, real fz){
  int x = fx;
  int y = fy;
  int z = fz;
  real u = fx-x;
  real v = fy-y;
  real w = fz-z;
  Complex res = {0,0};

  if(x >= sp_c3matrix_x(m)-1){
    x = sp_c3matrix_x(m)-2;
    u = 1;
  }
  if(y >= sp_c3matrix_y(m)-1){
    y = sp_c3matrix_y(m)-2;
    v = 1;
  }
  if(sp_c3matrix_z(m) > 1){
    if(z >= sp_c3matrix_z(m)-1){
      z = sp_c3matrix_z(m)-2;
      w = 1;
    }
  }else{
    z = 0;
    w = 0;
  }
  res = sp_cscale(sp_cscale(sp_cscale(sp_c3matrix_get(m,x,y,z),(1-u)),(1-v)),(1-w));
  if(u){res = sp_cadd(res,sp_cscale(sp_cscale(sp_cscale(sp_c3matrix_get(m,x+1,y,z),(u)),(1-v)),(1-w)));}
  if(v){res = sp_cadd(res,sp_cscale(sp_cscale(sp_cscale(sp_c3matrix_get(m,x,y+1,z),(1-u)),(v)),(1-w)));}
  if(w){res = sp_cadd(res,sp_cscale(sp_cscale(sp_cscale(sp_c3matrix_get(m,x,y,z+1),(1-u)),(1-v)),(w)));}
  if(u && v){res = sp_cadd(res,sp_cscale(sp_cscale(sp_cscale(sp_c3matrix_get(m,x+1,y+1,z),(u)),(v)),(1-w)));}
  if(u && w){res = sp_cadd(res,sp_cscale(sp_cscale(sp_cscale(sp_c3matrix_get(m,x+1,y,z+1),(u)),(1-v)),(w)));}
  if(v && w){res = sp_cadd(res,sp_cscale(sp_cscale(sp_cscale(sp_c3matrix_get(m,x,y+1,z+1),(1-u)),(v)),(w)));}
  if(u && v && w){res = sp_cadd(res,sp_cscale(sp_cscale(sp_cscale(sp_c3matrix_get(m,x+1,y+1,z+1),(u)),(v)),(w)));}
  return res;
}

/*! Returns the interpolated value of m at the floating point row frow and column fcol
 *
 */
static inline real sp_matrix_interp(const sp_matrix * m, real frow, real fcol){
  int x = fcol;
  int y = frow;
  real u = fcol-x;
  real v = frow-y;
  real res = 0;
  if(x >= sp_matrix_cols(m)-1){
    x = sp_matrix_cols(m)-2;
    u = 1;
  }
  if(y >= sp_matrix_rows(m)-1){
    y = sp_matrix_rows(m)-2;
    v = 1;
  }
  res = sp_matrix_get(m,y,x)*(1-u)*(1-v)+
    sp_matrix_get(m,y,x+1)*(u)*(1-v)+
    sp_matrix_get(m,y+1,x)*(1-u)*(v)+
    sp_matrix_get(m,y+1,x+1)*(u)*(v);
  return res;
}


/*! Returns the interpolated value of m at the floating point row frow and column fcol
 *
 */
static inline int sp_imatrix_interp(const sp_imatrix * m, real frow, real fcol){
  int x = fcol;
  int y = frow;
  real u = fcol-x;
  real v = frow-y;
  int res = 0;
  if(x >= sp_imatrix_cols(m)-1){
    x = sp_imatrix_cols(m)-2;
    u = 1;
  }
  if(y >= sp_imatrix_rows(m)-1){
    y = sp_imatrix_rows(m)-2;
    v = 1;
  }
  res = sp_imatrix_get(m,y,x)*(1-u)*(1-v)+
    sp_imatrix_get(m,y,x+1)*(u)*(1-v)+
    sp_imatrix_get(m,y+1,x)*(1-u)*(v)+
    sp_imatrix_get(m,y+1,x+1)*(u)*(v) + 0.5;
  return res;
}
static inline real sp_3matrix_interp(const sp_3matrix * m, real fx, real fy, real fz){
  int x = fx;
  int y = fy;
  real u = fx-x;
  real v = fy-y;
  real res = 0;
  if(fz == 0){
    if(x >= sp_3matrix_x(m)-1){
      x = sp_3matrix_y(m)-2;
      u = 1;
    }
    if(y >= sp_3matrix_y(m)-1){
      y = sp_3matrix_y(m)-2;
      v = 1;
    }
    res = sp_3matrix_get(m,x,y,0)*(1-u)*(1-v)+
      sp_3matrix_get(m,x+1,y,0)*(u)*(1-v)+
      sp_3matrix_get(m,x,y+1,0)*(1-u)*(v)+
      sp_3matrix_get(m,x+1,y+1,0)*(u)*(v);
  }else{
    int z = fz;
    real w = fz-z;
    if(x >= sp_3matrix_x(m)-1){
      x = sp_3matrix_x(m)-2;
      u = 1;
    }
    if(y >= sp_3matrix_y(m)-1){
      y = sp_3matrix_y(m)-2;
      v = 1;
    }
    if(z >= sp_3matrix_z(m)-1){
      z = sp_3matrix_z(m)-2;
      w = 1;
    }
    res = sp_3matrix_get(m,x,y,z)*(1-u)*(1-v)*(1-w)+
      sp_3matrix_get(m,x+1,y,z)*(u)*(1-v)*(1-w)+
      sp_3matrix_get(m,x,y+1,z)*(1-u)*(v)*(1-w)+
      sp_3matrix_get(m,x,y,z+1)*(1-u)*(1-v)*(w)+
      sp_3matrix_get(m,x+1,y+1,z)*(u)*(v)*(1-w)+
      sp_3matrix_get(m,x+1,y,z+1)*(u)*(1-v)*(w)+
      sp_3matrix_get(m,x,y+1,z+1)*(1-u)*(v)*(w)+
      sp_3matrix_get(m,x+1,y+1,z+1)*(u)*(v)*(w);
  }
  return res;
}


static inline int sp_i3matrix_interp(const sp_i3matrix * m, real fx, real fy, real fz){
  int x = fx;
  int y = fy;
  real u = fx-x;
  real v = fy-y;
  int res = 0;
  if(fz == 0){
    if(x >= sp_i3matrix_x(m)-1){
      x = sp_i3matrix_y(m)-2;
      u = 1;
    }
    if(y >= sp_i3matrix_y(m)-1){
      y = sp_i3matrix_y(m)-2;
      v = 1;
    }
    res = sp_i3matrix_get(m,x,y,0)*(1-u)*(1-v)+
      sp_i3matrix_get(m,x+1,y,0)*(u)*(1-v)+
      sp_i3matrix_get(m,x,y+1,0)*(1-u)*(v)+
      sp_i3matrix_get(m,x+1,y+1,0)*(u)*(v);
  }else{
    int z = fz;
    real w = fz-z;
    if(x >= sp_i3matrix_x(m)-1){
      x = sp_i3matrix_x(m)-2;
      u = 1;
    }
    if(y >= sp_i3matrix_y(m)-1){
      y = sp_i3matrix_y(m)-2;
      v = 1;
    }
    if(z >= sp_i3matrix_z(m)-1){
      z = sp_i3matrix_z(m)-2;
      w = 1;
    }
    res = sp_i3matrix_get(m,x,y,z)*(1-u)*(1-v)*(1-w)+
      sp_i3matrix_get(m,x+1,y,z)*(u)*(1-v)*(1-w)+
      sp_i3matrix_get(m,x,y+1,z)*(1-u)*(v)*(1-w)+
      sp_i3matrix_get(m,x,y,z+1)*(1-u)*(1-v)*(w)+
      sp_i3matrix_get(m,x+1,y+1,z)*(u)*(v)*(1-w)+
      sp_i3matrix_get(m,x+1,y,z+1)*(u)*(1-v)*(w)+
      sp_i3matrix_get(m,x,y+1,z+1)*(1-u)*(v)*(w)+
      sp_i3matrix_get(m,x+1,y+1,z+1)*(u)*(v)*(w);
  }
  return res;
}

/*! Resizes complex matrix m to the desired size. 
 *
 *  The content of the matrix will be destroyed.
 */
static inline void _sp_cmatrix_realloc(sp_cmatrix * m, int row, int col, char * file, int line){
  m->rows = row;
  m->cols = col;
  m->data = (Complex *)_sp_realloc(m->data,sizeof(Complex)*sp_cmatrix_size(m), file, line);
}
static inline void _sp_c3matrix_realloc(sp_c3matrix * m, int x, int y, int z, char * file, int line){
  m->x = x;
  m->y = y;
  m->z = z;
  m->data = (Complex *)_sp_realloc(m->data,sizeof(Complex)*sp_c3matrix_size(m),file,line);
}

/*! Resizes integer matrix m to the desired size. 
 *
 *  The content of the matrix will be destroyed.
 */
static inline void _sp_imatrix_realloc(sp_imatrix * m, int row, int col, char * file, int line){
  m->rows = row;
  m->cols = col;
  m->data = (int *)_sp_realloc(m->data,sizeof(int)*sp_imatrix_size(m),file,line);
}
static inline void _sp_i3matrix_realloc(sp_i3matrix * m, int x, int y, int z,char * file, int line){
  m->x = x;
  m->y = y;
  m->z = z;
  m->data = (int *)_sp_realloc(m->data,sizeof(Complex)*sp_i3matrix_size(m),file,line);
}

/*! Resizes matrix m to the desired size. 
 *
 *  The content of the matrix will be destroyed.
 */
static inline void _sp_matrix_realloc(sp_matrix * m, int row, int col,char * file, int line){
  m->rows = row;
  m->cols = col;
  m->data = (real *)_sp_realloc(m->data,sizeof(real)*sp_matrix_size(m),file,line);
}
static inline void _sp_3matrix_realloc(sp_3matrix * m, int x, int y, int z, char * file, int line){
  m->x = x;
  m->y = y;
  m->z = z;
  m->data = (real *)_sp_realloc(m->data,sizeof(real)*sp_3matrix_size(m),file,line);
}


/*! This function returns the index of a given row and column combination
 *
 */
static inline void sp_cmatrix_get_row_col(const sp_cmatrix * m, int index, int * row, int * col){
  *row = index%m->rows;
  *col = index/m->rows;
}

static inline void sp_c3matrix_get_xyz(const sp_c3matrix * m, long long index, int * x, int * y, int * z){
  *x = index%m->z%m->y;
  *y = index/m->x%m->z;
  *z = index/m->y/m->x;
}


/*! This function sets the real part to the cabs of each element and the complex part to 0
 *
 */
static inline void sp_cmatrix_to_real(const sp_cmatrix * m){
  int i;
  for(i = 0;i<sp_cmatrix_size(m);i++){
    sp_real(m->data[i]) = sp_cabs(m->data[i]);
    sp_imag(m->data[i]) = 0;
  }
}


/*! This function returns the Froenius inner product of complex matrix a and b
 *
 * The Froenius inner product is defined as the sum of a element by element multiplication
 * of the matrix elements. It's the basis for the Froenius norm.
 * Both matrices must obviously have the same dimensions.
 */
static inline Complex sp_cmatrix_froenius_prod(const sp_cmatrix * a, const sp_cmatrix * b){
  Complex ret = {0,0}; 
  int i;
  for(i = 0;i<sp_cmatrix_size(a);i++){
    ret = sp_cadd(ret,sp_cmul(a->data[i],sp_cconj(b->data[i])));
  }
  return ret;
}

static inline Complex sp_c3matrix_froenius_prod(const sp_c3matrix * a, const sp_c3matrix * b){
  Complex ret = {0,0};
  int i;
  for(i = 0;i<sp_c3matrix_size(a);i++){
    ret = sp_cadd(ret,sp_cmul(a->data[i],sp_cconj(b->data[i])));
  }
  return ret;
}



/*! Calculates the center of mass of the c3matrix a.
 *
 *  The formula used is R = 1/M * Sum(m_i * r_i) Where m_i is equal to sp_cabs(a->data[i]) and r_i is the position of i
 */
spimage_EXPORT sp_vector * sp_c3matrix_center_of_mass(sp_c3matrix * a);




#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

  
#endif
