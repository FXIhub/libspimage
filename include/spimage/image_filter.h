#ifndef _IMAGE_FILTER_H_
#define _IMAGE_FILTER_H_ 1

#include "image.h"
#include "linear_alg.h"

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

/** @defgroup ImageFilter Image Filtering
 *  Convolutions, Correlations and Filtering functions
 *  @{
 */

/*! Convolutes the images with a gaussian function of a given radius in pixels .
 * 
 * The filter function is given by:
 * \f$ f(x,y) = (2 \pi \times radius)^{-0.5} \times \exp{\frac{-(x^2+y^2)}{(2 \times radius^2)}} \f$
 *  The mask will not be blurred.
 */
spimage_EXPORT Image * gaussian_blur(Image * in, real radius);
/*! Multiplies the image with a gaussian centered in the image center of a given radius 
 *
 *   The mask will not be multiplied.
 */
spimage_EXPORT Image * gaussian_filter(Image * in, real radius, int in_place);

/*! Convolutes Image a with Image b 
 *
 *  Note that b must fit inside a or the other way around.
 *  The convolution only acts on the amplitudes, 
 *  the real and complex parts of the Image will be set to 0.
 *  The resulting mask will be the same as the one from image a.
 */
spimage_EXPORT Image * sp_image_convolute(Image * a, Image * b, int * size);

/*! Convolutes Image a with Image b 
 *
 *  Note that b must fit inside a or the other way around.
 *  The convolution only acts on the amplitudes, 
 *  the real and complex parts of the Image will be set to 0.
 *  The resulting mask will be the same as the one from image a.
 *
 * The convolution is calculated with better than pixel precision. 
 * precision is the inverse of the precision of the output. precision == 2 
 * corresponds to output with 1/2 pixel precision.
 * Flag tells whether to calculate convolution (flag == 0) or correlation (flag != 0)
 * The other parameters are similar to the normal convolution
 */
 spimage_EXPORT Image * sp_image_convolute_fractional(Image * a, Image * b, int * size,int precision, int flag);
  
/*! Low pass filter using a centered square window of side edge_size 
 *
 *  The mask is also low pass filtered.
 */

spimage_EXPORT Image * low_pass_square_filter(Image * in, int edge_size);
/*! Low pass filter using a centered gaussian window of side edge_size 
 *
 * The mask is also low pass filtered.
 */
spimage_EXPORT Image * low_pass_gaussian_filter(Image * in, int edge_size);

/*! Correlates Image a with Image b 
 *
 *  Note that b must fit inside a or the other way around.
 *  The correlation only acts on the amplitudes, 
 *  the real and complex parts of the Image will be set to 0.
 *  The resulting mask will be the same as the one from image a.
 *  If size is NULL no padding is done
 */
spimage_EXPORT Image * sp_image_cross_correlate(Image * a, Image * b, int * size);

/*! Returns a rectangular window centered on the middle of the image
 * with the given width and height. If shifted, the center is in the upper left corner
 * and the window wraps around the image
 */
spimage_EXPORT Image * rectangular_window(int image_x, int image_y, int width, int height, int shifted);

spimage_EXPORT Image * cube_window(int image_x, int image_y, int image_z, int dx, int dy, int dz, int shifted); 

/*! Returns a circular window centered on the image center
 * with the given radius
 */
spimage_EXPORT Image * circular_window(int x, int y, int radius, int shifted);
 
spimage_EXPORT Image * spherical_window(int x, int y, int z, int radius, int shifted); 

/*! Convolutes Image a with Image b at point i
 *
 *  The convolution only acts on the amplitudes. 
 *  This convolution, is not FFT based, and as such there is
 *  no wrap around effect. That is points of b which are outside
 *  of a are simply discarded.
 */
spimage_EXPORT real sp_point_convolute(const Image * a,const Image * b, int index);

/*! Returns the variance of the image in respect to the vicinity defined by window
 *
 * This function convolutes img with the window (making sure there's no aliasing)
 * and then returns the difference between this averaged image and the input img
 */
spimage_EXPORT Image * sp_image_local_variance(Image * img, Image * window);
/*@}*/

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif
