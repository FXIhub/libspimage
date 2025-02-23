%module spimage_pybackend

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
%{
  /* Includes the header in the wrapper code */
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "../include/spimage/colormap.h"
#include "../include/spimage/cuda_util.h"
#include "../include/spimage/fft.h"
#include "../include/spimage/hashtable.h"
#include "../include/spimage/image_filter_cuda.h"
#include "../include/spimage/image_filter.h"
#include "../include/spimage/image.h"
#include "../include/spimage/image_io.h"
#include "../include/spimage/image_noise.h"
#include "../include/spimage/image_sphere.h"
#include "../include/spimage/image_util.h"
#include "../include/spimage/interpolation_kernels.h"
#include "../include/spimage/linear_alg.h"
#include "../include/spimage/list.h"
#include "../include/spimage/map.h"
#include "../include/spimage/mem_util.h"
#include "../include/spimage/phasing.h"
#include "../include/spimage/prtf.h"
#include "../include/spimage/sperror.h"
#include "../include/spimage/statistics.h"
#include "../include/spimage/support_update.h"	
#include "../include/spimage/time_util.h"
  //#include <numpy/npy_common.h>
#include <numpy/arrayobject.h>
  
  %}

%init %{
        import_array();
%}

/* Parse the header file to generate wrappers */
%typemap(out) Complex {
    $result = PyComplex_FromDoubles($1.re,$1.im);
}
/* Parse the header file to generate wrappers */
%typemap(in) Complex {
    $1.re = PyComplex_RealAsDouble($input);
    $1.im = PyComplex_ImagAsDouble($input);
}


%typemap(out) sp_i3matrix * {
  if($1->z == 1){
    /* Swap the order of the dimensions so we can plot things easily in imshow.
       This is bound to create confusion at some point in the future. */
    npy_intp dims[2] = {$1->y,$1->x};
    $result = PyArray_SimpleNewFromData(2, dims, NPY_INT, $1->data);
  }else{
    npy_intp dims[3] = {$1->z,$1->y,$1->x};
    $result = PyArray_SimpleNewFromData(3, dims, NPY_INT, $1->data);
  }
}

%typemap(in) sp_i3matrix * {
  /* Make sure input is a NumPy array and an Int */
  PyArrayObject *arr;
  if(PyArray_Check($input) == 0 || PyArray_ISINTEGER((PyArrayObject *)$input) == 0){
    PyErr_SetString( PyExc_TypeError, "not an int array" );
    return NULL;
  } 
  arr = (PyArrayObject *)($input);
  //  if ((arr = (PyArrayObject *) PyArray_ContiguousFromObject($input,NPY_INT, 0, 0)) == NULL) return NULL;
  npy_intp * dim = PyArray_DIMS(arr);
  $1->x = dim[0];
  $1->y = dim[1];
  $1->z = dim[2];
  //  $1->data = (int *)arr->data;
  $1->data = (int *) PyArray_DATA (arr);
  $input = (PyObject *)arr;
}

%typemap(out) sp_c3matrix * {
  if($1->z == 1){
    /* Swap the order of the dimensions so we can plot things easily in imshow.
       This is bound to create confusion at some point in the future. */
    npy_intp dims[2] = {$1->y,$1->x};
    $result = PyArray_SimpleNewFromData(2, dims, NPY_CFLOAT, $1->data);
  }else{
    npy_intp dims[3] = {$1->z,$1->y,$1->x};
    $result = PyArray_SimpleNewFromData(3, dims, NPY_CFLOAT, $1->data);
  }
}

%typemap(out) sp_3matrix * {
  if($1->z == 1){
    /* Swap the order of the dimensions so we can plot things easily in imshow.
       This is bound to create confusion at some point in the future. */
    npy_intp dims[2] = {$1->y,$1->x};
    $result = PyArray_SimpleNewFromData(2, dims, NPY_FLOAT, $1->data);
  }else{
    npy_intp dims[3] = {$1->z,$1->y,$1->x};
    $result = PyArray_SimpleNewFromData(3, dims, NPY_FLOAT, $1->data);
  }
}

%typemap(out) float[3] {
  /*$result = PyTuple_New(3);
  PyTuple_SetItem($result, 0, PyFloat_FromDouble($1[0]));
  PyTuple_SetItem($result, 1, PyFloat_FromDouble($1[1]));
  PyTuple_SetItem($result, 2, PyFloat_FromDouble($1[2]));*/
  npy_intp dims[1] = {3};
  $result = PyArray_SimpleNewFromData(1, dims, NPY_FLOAT, $1);
}


%typemap(in) int * {
  PyArrayObject *arr;
  if(PyArray_Check($input) == 0 || PyArray_ISINTEGER((PyArrayObject *)$input) == 0 
     || PyArray_TYPE((PyArrayObject *)$input) != NPY_INT){
    PyErr_SetString( PyExc_TypeError, "int * argument not an int32 numpy array" );
    return NULL;
  } 
  arr = (PyArrayObject *)($input);
  //  $1 = (int *)arr->data;
  $1 = (int *) PyArray_DATA (arr);
}

%typemap(in) float * {
  PyArrayObject *arr;
  if($input == Py_None){
    $1 = (float *)NULL;
  }else{ 
    if(PyArray_Check($input) == 0 || PyArray_ISFLOAT((PyArrayObject *)$input) == 0 
     || PyArray_TYPE((PyArrayObject *)$input) != NPY_FLOAT32){
    PyErr_SetString( PyExc_TypeError, "float * argument not an float32 numpy array" );
      return NULL;
    } 
    arr = (PyArrayObject *)($input);
    //  $1 = (int *)arr->data;
    $1 = (float *) PyArray_DATA (arr);
  }
}

%typemap(in) Image ** {
  $1 = NULL;
  if (PyList_Check($input)) {
    const size_t size = PyList_Size($input);
    $1 = (Image**)malloc((size+1) * sizeof(Image*));
    for (int i = 0; i < size; ++i) {
      void *argp = 0 ;
      const int res = SWIG_ConvertPtr(PyList_GetItem($input, i), &argp, $*1_descriptor, 0);
      if (!SWIG_IsOK(res)) {
        SWIG_exception_fail(SWIG_ArgError(res), "in method '" "$symname" "', argument " "$argnum"" of type '" "$1_type""'");
      }
      $1[i] = (Image*)(argp);
    }
    $1[size] = NULL;
  }
  else {
    // Raise exception
    SWIG_exception_fail(SWIG_TypeError, "Expected list in $symname");
  }
}


%include "../include/spimage/image.h"
%include "../include/spimage/colormap.h"
%include "../include/spimage/cuda_util.h"
%include "../include/spimage/fft.h"
%include "../include/spimage/hashtable.h"
%include "../include/spimage/image_filter_cuda.h"
%include "../include/spimage/image_filter.h"
%include "../include/spimage/image_io.h"
%include "../include/spimage/image_noise.h"
%include "../include/spimage/image_sphere.h"
%include "../include/spimage/image_util.h"
%include "../include/spimage/interpolation_kernels.h"
%include "../include/spimage/linear_alg.h"
%include "../include/spimage/list.h"
%include "../include/spimage/map.h"
%include "../include/spimage/mem_util.h"
%include "../include/spimage/phasing.h"
%include "../include/spimage/prtf.h"
%include "../include/spimage/sperror.h"
%include "../include/spimage/statistics.h"
%include "../include/spimage/support_update.h"	
%include "../include/spimage/time_util.h"
%include "../include/spimage/find_center.h"

/* %include "../include/spimage/image.h" */
/* %include "../include/spimage/fft.h" */
/* %include "../include/spimage/image_util.h" */
/* %include "../include/spimage/linear_alg.h" */
/* %include "../include/spimage/mem_util.h" */
/* %include "../include/spimage/image_sphere.h" */
/* %include "../include/spimage/sperror.h" */
/* %include "../include/spimage/statistics.h" */
/* %include "../include/spimage/image_noise.h" */
/* %include "../include/spimage/hashtable.h" */
/* %include "../include/spimage/interpolation_kernels.h" */
/* %include "../include/spimage/time_util.h" */
/* %include "../include/spimage/list.h" */
/* %include "../include/spimage/map.h" */
/* %include "../include/spimage/prtf.h" */
/* %include "../include/spimage/phasing.h" */
