#ifndef _CIMAGE_H_
#define _CIMAGE_H_


#ifdef DOUBLE
typedef double real;
#else
typedef float real;
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
