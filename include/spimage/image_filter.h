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


/*! Convolute the image with a square window.
 *  The filter function is given by:
 *  
 *  f(x,y) = 1/((2*radius+1)^2)) 
 *
 *  If type is SP_3D the blur is done in 3D.
 *
 * \param in input image
 * \param radius distance from the center of the square/cube to the center of the face
 * \param type if given SP_3D the blur is done in 3D. Otherwise it's done in 2D.
 * \return the blurred image
 */
  Image * sp_square_blur(Image * in, real radius, int type);

/*! Convolutes the images with a gaussian function of a given radius in pixels .
 * 
 * \param in input image
 * \param radius standard deviation of the gaussian kernel in pixels
 * \return the blurred image
 * 
 * The filter function is given by:
 * \f$ f(x,y) = (2 \pi \times radius)^{-0.5} \times \exp{\frac{-(x^2+y^2)}{(2 \times radius^2)}} \f$
 *  The mask will not be blurred.
 */
spimage_EXPORT Image * sp_gaussian_blur(Image * in, real radius);

/*! Multiplies the image with a gaussian centered in the image center of a given radius 
 *
 * \param in input image
 * \param radius radius of the gaussian function used. This is defined as the point where the gaussian drops below 0.0009
 * \param in_place if true the multiplication is done in place
 * \return the input image multiplied with a gaussian image
 *
 *   The mask will not be multiplied.
 */
spimage_EXPORT Image * sp_gaussian_filter(Image * in, real radius, int in_place);

/*! Convolutes Image a with Image b 
 * 
 * \param a input image a
 * \param b input image b
 * \param size a 3 element array describing the size to which the images will be zero padded before the convolution. If NULL then the size of the largest input image is used.
 * \return the convolution of a with b after properly padding the images.
 *
 *  Note that b must fit inside of a or the other way around.
 *  The convolution only acts on the amplitudes, 
 *  the real and complex parts of the Image will be set to 0.
 *  The resulting mask will be the same as the one from image a.
 */
spimage_EXPORT Image * sp_image_convolute(Image * a, Image * b, int * size);

/*! Convolutes Image a with Image b 
 *
 * \param a input image a
 * \param b input image b
 * \param size a 3 element array describing the size to which the images will be zero padded before the convolution. If NULL then the size of the largest input image is used.
 * \param precision multiply the size of the output image by precision before back fourier transforming, which gives better than pixel precision by interpolation
 * \param flag if true do a correlation, otherwise do a convolution
 * \return the convolution of a with b after properly padding the images.
 *
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
 *
 * \deprecated Please use the sp_image_convolute and sp_image_cross_correlate with the size parameter to achieve the same result
 * \sa sp_image_convolute sp_image_cross_correlate
 */
 spimage_EXPORT Image * sp_image_convolute_fractional(Image * a, Image * b, int * size,int precision, int flag);
  
/*! Low pass filter using a centered square window of side edge_size 
 *
 * \param in the input image
 * \param edge_size the side of the square window in pixels
 * \return the input image multiplied by a square window
 *
 *  The mask is also low pass filtered.
 *
 */
spimage_EXPORT Image * sp_low_pass_square_filter(Image * in, int edge_size);

/*! Low pass filter using a centered gaussian window of side edge_size 
 *
 *
 * \param in the input image
 * \param edge_size the side of the gaussian window in pixels
 * \return the input image multiplied by a gaussian window
 *
 * The edge of the gaussian is defined as the point where it's value drops below 0.0009
 * The mask is also low pass filtered.
 */
spimage_EXPORT Image * sp_low_pass_gaussian_filter(Image * in, int edge_size);

/*! Correlates Image a with Image b 
 *
 * \param a input image
 * \param b input image
 * \param size a 3 element array describing the size to which the images will be zero padded before the convolution. If NULL then the size of the largest input image is used.
 * \return image a correlated with b using the requested zero padding
 *
 *  Note that b must fit inside a or the other way around.
 *  The correlation only acts on the amplitudes, 
 *  the real and complex parts of the Image will be set to 0.
 *  The resulting mask will be the same as the one from image a.
 *  If size is NULL no padding is done
 */
spimage_EXPORT Image * sp_image_cross_correlate(Image * a, Image * b, int * size);

/*! Returns a rectangular window centered on the middle of the image
 * with the given width and height. 
 *
 * \param image_x the width of the output image
 * \param image_y the height of the output image
 * \param width the width of the region where the filter is 1
 * \param height the height of the region where the filter is 1
 * \param shifted if true the center is in the upper left corner and the window wraps around the image
 * \return a filter function that is 1 inside a given rectangle and 0 outside, possibly shifted
 *
 * \sa sp_cube_window
 */
spimage_EXPORT Image * sp_rectangular_window(int image_x, int image_y, int width, int height, int shifted);

/*! Returns a cubic window centered on the middle of the image
 * with the given width and height and depth. 
 *
 * \param image_x the width of the output image
 * \param image_y the height of the output image
 * \param image_z the depth of the output image
 * \param width the width of the region where the filter is 1
 * \param height the height of the region where the filter is 1
 * \param depth the depth of the region where the filter is 1
 * \param shifted if true the center is in the upper left corner and the window wraps around the image
 * \return a filter function that is 1 inside a given cube and 0 outside, possibly shifted
 *
 * \sa sp_rectangular_window
 */
spimage_EXPORT Image * sp_cube_window(int image_x, int image_y, int image_z, int width, int height, int depth, int shifted); 

/*! Returns a circular window centered on the image center
 * with the given radius.
 * \param x the width of the output image
 * \param y the height of the output image
 * \param radius the radius of the circle where the filter is 1
 * \param shifted if true the center is in the upper left corner and the window wraps around the image
 * \return a filter function that is 1 inside a given circle and 0 outside, possibly shifted
 *
 * \sa sp_spherical_window
 */
spimage_EXPORT Image * sp_circular_window(int x, int y, int radius, int shifted);

/*! Returns a circular window centered on the image center
 * with the given radius.
 * \param x the width of the output image
 * \param y the height of the output image
 * \param z the depth of the output image
 * \param radius the radius of the circle where the filter is 1
 * \param shifted if true the center is in the upper left corner and the window wraps around the image
 * \return a filter function that is 1 inside a given sphere and 0 outside, possibly shifted
 *
 * \sa sp_circular_window
 */
spimage_EXPORT Image * sp_spherical_window(int x, int y, int z, int radius, int shifted); 

/*! Convolutes Image a with Image b at point i
 *
 * \param a input image
 * \param b input image
 * \param index index of the pixel where to calculate the convolution
 * \return the value of the convoluted pixels
 *
 *  The convolution only acts on the amplitudes. 
 *  This convolution, is not FFT based, and as such there is
 *  no wrap around effect. That is points of b which are outside
 *  of a are simply discarded.
 */
spimage_EXPORT real sp_point_convolute(const Image * a,const Image * b, int index);

/*! Returns the variance of the image in respect to the vicinity defined by window
 *
 * \param img the input image
 * \param window the filter that the input is convoluted with
 * \return The absolute difference between the input image and a smoothed version of itself, that is \f$out =  |img - img \star window| \f$
 *
 * This function convolutes img with the window (making sure there's no aliasing)
 * and then returns the absolute difference between this averaged image and the input img
 */
spimage_EXPORT Image * sp_image_local_variance(Image * img, Image * window);
/*@}*/

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif
