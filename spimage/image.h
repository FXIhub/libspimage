#ifndef _CIMAGE_H_
#define _CIMAGE_H_

#include <float.h>

#ifdef DOUBLE
typedef double real;
typedef _Complex double complex;
#define REAL_MAX DBL_MAX
#define REAL_EPSILON DBL_EPSILON
#else
typedef float real;
typedef _Complex float complex;
#define REAL_MAX FLT_MAX
#define REAL_EPSILON FLT_EPSILON
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

/*! Structure that keeps all the information about the Detector type.

*/

typedef struct{
  int size[2];
  real image_center[2];
  real pixel_size;
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
  /* Intensity of the image */
  real * intensities;
  /* amplitude of the image */
  real * amplitudes;
  /* This always point to the correct
   array depending on scaled (meaning it points
  to intensities when scaled is on and vice versa)*/
  real * image;
  /* real part of the image */
  real * r;
  /* complex part */
  real * c;
  real * mask;
  Detector * detector;
  /* this flag tells wether the image is shifted 
     (in FFTW format, center on the corners) or not*/
  int shifted;
}Image;

#endif
