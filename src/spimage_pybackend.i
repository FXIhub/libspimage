%module spimage_pybackend
%{
  /* Includes the header in the wrapper code */
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
    $result = PyArray_SimpleNewFromData(2, dims, PyArray_INT, $1->data);
  }else{
    npy_intp dims[3] = {$1->z,$1->y,$1->x};
    $result = PyArray_SimpleNewFromData(3, dims, PyArray_INT, $1->data);
  }
}

%typemap(in) sp_i3matrix * {
  /* Make sure input is a NumPy array and an Int */
  PyArrayObject *arr;
  if(PyArray_Check($input) == 0 || PyArray_ISINTEGER($input) == 0){
    PyErr_SetString( PyExc_TypeError, "not an int array" );
    return NULL;
  } 
  arr = (PyArrayObject *)($input);
  //  if ((arr = (PyArrayObject *) PyArray_ContiguousFromObject($input,PyArray_INT, 0, 0)) == NULL) return NULL;
  npy_intp * dim = PyArray_DIMS(arr);
  $1->x = dim[0];
  $1->y = dim[1];
  $1->z = dim[2];
  $1->data = (int *)arr->data;
  $input = (PyObject *)arr;
}

%typemap(out) sp_c3matrix * {
  if($1->z == 1){
    /* Swap the order of the dimensions so we can plot things easily in imshow.
       This is bound to create confusion at some point in the future. */
    npy_intp dims[2] = {$1->y,$1->x};
    $result = PyArray_SimpleNewFromData(2, dims, PyArray_CFLOAT, $1->data);
  }else{
    npy_intp dims[3] = {$1->z,$1->y,$1->x};
    $result = PyArray_SimpleNewFromData(3, dims, PyArray_CFLOAT, $1->data);
  }
}

%typemap(out) sp_3matrix * {
  if($1->z == 1){
    /* Swap the order of the dimensions so we can plot things easily in imshow.
       This is bound to create confusion at some point in the future. */
    npy_intp dims[2] = {$1->y,$1->x};
    $result = PyArray_SimpleNewFromData(2, dims, PyArray_FLOAT, $1->data);
  }else{
    npy_intp dims[3] = {$1->z,$1->y,$1->x};
    $result = PyArray_SimpleNewFromData(3, dims, PyArray_FLOAT, $1->data);
  }
}

%typemap(out) float[3] {
  /*$result = PyTuple_New(3);
  PyTuple_SetItem($result, 0, PyFloat_FromDouble($1[0]));
  PyTuple_SetItem($result, 1, PyFloat_FromDouble($1[1]));
  PyTuple_SetItem($result, 2, PyFloat_FromDouble($1[2]));*/
  npy_intp dims[1] = {3};
  $result = PyArray_SimpleNewFromData(1, dims, PyArray_FLOAT, $1);
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
