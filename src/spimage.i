%module spimage
 %{
 /* Includes the header in the wrapper code */
#include "../include/spimage/image.h"
#include "../include/spimage/fft.h"
#include "../include/spimage/image_util.h"
#include "../include/spimage/linear_alg.h"
#include "../include/spimage/mem_util.h"
#include "../include/spimage/image_sphere.h"
#include "../include/spimage/sperror.h"
#include "../include/spimage/statistics.h"
#include "../include/spimage/image_noise.h"
#include "../include/spimage/hashtable.h"
#include "../include/spimage/interpolation_kernels.h"
#include "../include/spimage/time_util.h"
#include "../include/spimage/list.h"
#include "../include/spimage/map.h"
#include "../include/spimage/prtf.h"
#include "../include/spimage/phasing.h"
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


%include "../include/spimage/image.h"
%include "../include/spimage/fft.h"
%include "../include/spimage/image_util.h"
%include "../include/spimage/linear_alg.h"
%include "../include/spimage/mem_util.h"
%include "../include/spimage/image_sphere.h"
%include "../include/spimage/sperror.h"
%include "../include/spimage/statistics.h"
%include "../include/spimage/image_noise.h"
%include "../include/spimage/hashtable.h"
%include "../include/spimage/interpolation_kernels.h"
%include "../include/spimage/time_util.h"
%include "../include/spimage/list.h"
%include "../include/spimage/map.h"
%include "../include/spimage/prtf.h"
%include "../include/spimage/phasing.h"


