/*

A few remarks about coordinate systems.

The top left pixel of a detector is called the pixel (0,0) and the lower right the pixel (max_x,max_y)
The top left corner of a detector has the following coordinates: (-0.5, -0.5)
The bottom right corner of a detector has the following coordinates: (max_x+0.5, max_y+0.5)
This results from the finite size of a pixel.

The coordinate units are pixels.

So it follows that the center of the detector is (max_x/2,max_y/2) and its dimensions are (max_x+1,max_y+1)

*/

#ifndef _CIMAGE_H_
#define _CIMAGE_H_

#include <float.h>
#include <complex.h>

#ifdef _SP_DOUBLE_PRECISION
typedef double real;
typedef _Complex double Complex;
#define REAL_MAX DBL_MAX
#define REAL_EPSILON DBL_EPSILON
#define conjr(a) conj(a)
#define cabsr(a) cabs(a)
#define cargr(a) carg(a)
#else
typedef float real;
typedef _Complex float Complex;
#define REAL_MAX FLT_MAX
#define REAL_EPSILON FLT_EPSILON
#define cabsr(a) cabsf(a)
#define conjr(a) conjf(a)
#define cargr(a) cargf(a)
#endif


#if defined(_WIN32) || defined(WIN32) /* Win32 version */
#ifdef spimage_EXPORTS
#  define spimage_EXPORT __declspec(dllexport)
#else
#  define spimage_EXPORT __declspec(dllimport)
#endif
#else
/* unix needs nothing */
#define spimage_EXPORT
#endif

typedef enum{SP_1D=1,SP_2D=2,SP_3D=3} Dimensions;

#include "linear_alg.h"

/*! Structure that keeps all the information about the Detector type.

*/

typedef struct{
  real image_center[3];
  real pixel_size[3];
  real detector_distance;
  real lambda;
}Detector;


/*! Main structure that keeps all the information about an image.

*/


typedef struct{
  /* this flag tell us whether we should try to access
   the real and complex part of the image seperately or
   just the amplitudes */
  int phased;  
  /* this flag tell us whether we should try to access
   the amplitudes (when scaled) or the intensities
  (when unscalled) */
  int scaled;
  /* The actual image */
  sp_c3matrix * image;
  /* The integer mask */
  sp_i3matrix * mask;
  Detector * detector;
  /* this flag tells wether the image is shifted 
     (in FFTW format, center on the corners) or not*/
  int shifted;
  sp_3matrix * rec_coords;
  Dimensions num_dimensions;

}Image;

#endif
