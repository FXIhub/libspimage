#ifndef _COLORMAP_H_
#define _COLORMAP_H_

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

#include "spimage.h"

typedef struct{
  real r;
  real g;
  real b;
}sp_rgb;

/*! 
  Calculates the color corresponding to the given value using a certain colormap.
  value must fall in the interval [0,1[.
*/
spimage_EXPORT sp_rgb sp_colormap_rgb_from_value(real value, int colormap);


spimage_EXPORT void sp_colormap_create_table(sp_rgb color_table[256],int colormap);

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */


#endif 
