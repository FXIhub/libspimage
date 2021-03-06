#ifndef _IMAGE_H_
#define _IMAGE_H_ 1

#include <float.h>

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */


#if defined(_WIN32) || defined(WIN32) /* Win32 version */
#ifdef SPIMAGE_NO_DLL
  /* Nothing happens */
  //#  define spimage_EXPORT __declspec(dllimport)
#define spimage_EXPORT
#else
  /* export as DLL */
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
  /*! Center of the image in pixels */
  real image_center[3]; 
  /*! Size of the pixels in meters*/
  real pixel_size[3]; 
  /*! Distance between the scatterer and the detector in meters   */
  real detector_distance;
  /*! Photon wavelength used in meters */
  real wavelength;
  /*! The orientation of the detector */
  SpRotation * orientation;
}Detector;


/*! Main structure that keeps all the information about an image.

A few remarks about coordinate systems.

The top left pixel of a detector is called the pixel (0,0) and the lower right the pixel (max_x,max_y)
The top left corner of a detector has the following coordinates: (-0.5, -0.5)
The bottom right corner of a detector has the following coordinates: (max_x+0.5, max_y+0.5)
This results from the finite size of a pixel.

The coordinate units are pixels.

So it follows that the center of the detector is (max_x/2,max_y/2) and its dimensions are (max_x+1,max_y+1)

Please note that max_x == sp_image_x()-1.

*/
typedef struct{
  /*! tells us whether we should try to access
   the real and complex part of the image seperately or
   just the amplitudes */
  int phased;  
  /*! tells us whether we should try to access
   the amplitudes (when scaled) or the intensities
  (when unscalled) */
  int scaled;
  /*! The actual image */
  sp_c3matrix * image;
  /*! The integer mask */
  sp_i3matrix * mask;
  Detector * detector;
  /*! this flag tells wether the image is shifted 
     (in FFTW format, center on the corners) or not*/
  int shifted;
  /*! if set, contain the coordinates in reciprocal space, in \f$ m^{-1} \f$ */
  sp_3matrix * rec_coords;
  /*! the dimensionality of the image */
  Dimensions num_dimensions;
}Image;

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */


#endif
