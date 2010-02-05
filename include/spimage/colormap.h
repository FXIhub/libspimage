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

  typedef enum{SpColormapFirstColorScheme=1,SpColormapGrayScale=1,SpColormapTraditional=2,
	     SpColormapHot=4,SpColormapRainbow=8,SpColormapJet=16,
	       SpColormapWheel=32,SpColormapLastColorScheme=64,SpColormapLogScale=128,SpColormapPhase=256,
	       SpColormapWeightedPhase=512,SpColormapMask=1024,SpColormapShadedMask=2048}SpColormap;

/*! 
  Calculates the color corresponding to the given value using a certain colormap.
  value must fall in the interval [0,1[.
*/
spimage_EXPORT sp_rgb sp_colormap_rgb_from_value(real value, int colormap);


spimage_EXPORT void sp_colormap_create_table(sp_rgb color_table[256],int colormap);

  spimage_EXPORT real sp_colormap_scale_value(real input,int colormap,real max_v, real min_v);

  spimage_EXPORT void sp_colormap_write_rgb(unsigned char * out,Image * img, int colormap,sp_rgb * color_table,real max_v, real min_v, int x, int y, int z, int red_blue_swap);
#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */


#endif 
