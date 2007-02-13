#ifndef _IMAGE_UTIL_H_
#define _IMAGE_UTIL_H_

#include "image.h"

#define TSIZE(a) (a->detector->size[0]*a->detector->size[1])

#define COLOR_GRAYSCALE 1
#define COLOR_TRADITIONAL 2
#define COLOR_HOT 4
#define COLOR_RAINBOW 8
#define COLOR_JET 16
#define LOG_SCALE 32

#define OUT_OF_PLACE 0
#define IN_PLACE 1

#ifndef MAX
#define MAX(a,b) ((a)>(b)?(a):(b))
#endif

#ifndef MIN
#define MIN(a,b) ((a)<(b)?(a):(b))
#endif

/** @defgroup Distance
 *  Calculates several kinds of distances in an image
 *  @{
 */


/*! Calculates the distance from the pixel corresponding to index i to the closest corner.  
 */
spimage_EXPORT real dist_to_corner(int i, Image * in);
/*! Calculates the distance from the pixel corresponding to index i to the image center,
  as defined in the Detector.  
 */
spimage_EXPORT real dist_to_center(int i, Image * in);
/*! Calculates the distance from the pixel corresponding to index i to closest axis going
through the image center as defined in the Detector.  
 */
spimage_EXPORT real dist_to_axis(int i,Image * in);
/*! Calculates the "Manhattan distance" from the pixel corresponding to index i to the image center,
  as defined in the Detector.  
 */
spimage_EXPORT real square_dist_to_center(int i,Image * in);
/*@}*/



/** @defgroup IO Input-Output
 *  Image input output routines
 *  @{
 */

/*! Read image from a file in spimage file format.
 *
 * It alocates and returns an Image structre will all the information
 * read from the file filename. It returns a NULL pointer
 * if it cannot open the file, or the file is not an spimage file.
 */

spimage_EXPORT Image * read_imagefile(const char * filename);

/*! Write an image from to a file in spimage format.
 *
 * It creates an spimage file with all the informations contained
 * in the Image structure with the specified precision.
 * output_precision must be 4 or 8 and specifies the number of bytes
 * used for storing each floating point number.
 */
spimage_EXPORT void write_img(Image * img,char * filename,int output_precision);

/*! Create an image using the values from the TIFF file
*/
spimage_EXPORT Image * read_tiff(const char * filename);
/*! Write a TIFF using the values from the image file
*/
spimage_EXPORT void write_tiff(Image * img, char * filename);

/*! Write the norm of the image to a file in png format, using 
  the specified color map
*/
spimage_EXPORT int write_png(Image * img, char * filename, int color);

/* Read png file and convert to image */
spimage_EXPORT Image*  read_png(char * filename);

/*! Write the norm of the image to a file in VTK ASCII format as a 2D Structured grid
 */
spimage_EXPORT int write_vtk(Image * img, char * filename);

/*! Write the mask of the image to a file in png format, using 
  the specified color map
*/
spimage_EXPORT int write_mask_to_png(Image * img, char * filename, int color);
/*@}*/



/** @defgroup MM Memory Management
 *  Create, copy and free images
 *  @{
 */

/*! Create a copy of Image in
 */
spimage_EXPORT Image * imgcpy(Image * in);
/*! Delete image Image in
 */
spimage_EXPORT void freeimg(Image * in);

/*! Creates an Image with the same dimensions as a.
 *  The value of the pixels is all set to 0
 */
spimage_EXPORT Image * create_empty_img(Image * a);

/*! Creates an empty image with given dimensions. Mask is set to 0 everywhere.
 *  Image in unscaled and unphased.
 */
spimage_EXPORT Image * create_new_img(int x, int y);

/*! Copies the contents of src to dst.
 *
 * Both images must have the same dimensions!
 */
spimage_EXPORT void sp_image_memcpy(Image * dst,Image *src);

/*@}*/



/** @defgroup ImageCrop Image Cropping
 *  Crops parts of the image
 *  @{
 */

/*! Masks out all pixels whose distance to center
  is smaller than radius. It also puts their values to 0.
 */
spimage_EXPORT void remove_lowres(Image * in, real radius);
/*! Masks out all pixels whose "manhattan distance" to center
  is bigger than radius. It also puts their values to 0.
 */
spimage_EXPORT Image * limit_resolution(Image * img, int resolution);
/*! Creates a new image from the rectange define by the upper left corner (x1,y1)
  and lower right corner (x2,y2).
 */
spimage_EXPORT Image * rectangle_crop(Image * in, int x1, int y1, int x2, int y2); 

/*! Returns a ray from the image along a given direction

   Returns the values of the image along a radial
   sector in a certain direction(in radians)
   samples defines how many values we return
   
   The origin of the vector is given by point,
   or by the image center if point is NULL.

   intersection old the position of the image where the 
   sector vector intersect the image border.

   The values are obtained by linear interpolation
 */
spimage_EXPORT Image * get_image_radial_sector(Image * img, real * point, real direction, int samples, real * intersection);

/*! Gives the distance to the border on a given direction

  Returns the value of the distance between
  "point" (Or the img center in case "point" is NULL)
  and the border of the image in a certain direction
  specified by an angle (in radians)
  The point of intersection is returned
  in "intersection".
*/
spimage_EXPORT real get_image_radial_distance_to_border(Image * img, real * point, real direction, real * intersection);

/*! Creates a circular symmetric window by rotating a sector around center

  Takes a sector and rotates it around the center to create an image.
  It assumes the sector is pixel scaled (each bin 1 pixel ).
  That is it does not try to stretch the sector to cover the image.
*/
spimage_EXPORT Image * circular_image_from_sector(Image * sector, int * size, real * center);


/*@}*/


/** @defgroup ImageAlgebra Image Algebra
 *  Basic algebra with images
 *  @{
 */
/*! Add Image b to Image a
 */
spimage_EXPORT void add_image(Image * a, Image * b);
/*! Subtract Image b from Image a
 */
spimage_EXPORT void sp_image_sub(Image * a, Image * b);

/*! Mutiplies Image a with Image b element by element
 */
spimage_EXPORT Image * sp_image_mul_elements(Image * a, Image * b);

/*! Transforms an image in its complex conjugate
 */
spimage_EXPORT void sp_image_complex_conj(Image * a);

/*! Sets a to 1/a. 
 */
spimage_EXPORT int sp_image_invert(Image * a);


/*! Calculates the dot product of a and b
 *
 * The real part of the result is return in *r and the imaginary part in *j
 * A return value of 0 means success, otherwise it indicates an error.
 */
spimage_EXPORT int sp_image_dot_prod(Image * a, Image * b, real * r, real * j);

/*! Returns the phase of a complexed value pixel.
 */
spimage_EXPORT real get_phase_angle(Image * img, int i);

/*! Calculates the \f$ \cos() \f$ of the complexed value pixel i.
 *  Can also be understood as \f$ \frac{Re(f(i))}{|f(i)|} \f$
 */
spimage_EXPORT real real_factor(Image * img, int i);
/*! Calculates the \f$ \sin() \f$ of the complexed value pixel i.
 *  Can also be understood as \f$ \frac{Im(f(i))}{|f(i)|} \f$
 */
spimage_EXPORT real complex_factor(Image * img, int i);
/*! Calculates the norm of the complexed value pixel i.
 */
spimage_EXPORT real norm(Image * img, int i);
/*! Calculates the integral of all the image, independently of the mask.
 */
spimage_EXPORT real integrated_intensity(Image * a);

/*! Turns a real valued image into a complexed valued one,
  mantaining the real part and setting the imaginary part to 0.
 */
spimage_EXPORT void rephase(Image *  img);

/*! Turns a complexed valued image into a real valued one, by
  taking the norm of each pixel.
 */
spimage_EXPORT void dephase(Image *  img);

/*! Turns a real valued image into a complexed valued one,
  mantaining the norm and using a random phase.
 */
spimage_EXPORT void random_rephase(Image * img);

/*! Multiply Image img by a scalar value.
 */
spimage_EXPORT void multiply_scalar_by_image(Image * img, real value);
/*! Returns the index of the image maximum, and set *max to the
  value of that pixel
 */
spimage_EXPORT int get_image_maximum(Image * img, int * x, int * y, real * max);

/*! Returns the phase of the image on each pixel, in radians
 *
 * Returns the phase of a in radians, from [0, 2*M_PI[ or NULL if there's an error
 */
spimage_EXPORT Image * sp_image_get_phases(Image * a);

/*@}*/


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
spimage_EXPORT Image * convolute_img(Image * a, Image * b, int * size);

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
 */
spimage_EXPORT Image * cross_correlate_img(Image * a, Image * b, int * size);

/*! Returns a rectangular window centered on the middle of the image
 * with the given width and height. If shifted, the center is in the upper left corner
 * and the window wraps around the image
 */
spimage_EXPORT Image * rectangular_window(int image_x, int image_y, int width, int height, int shifted);


/*! Returns a circular window centered on the image center
 * with the given radius
 */
spimage_EXPORT Image * circular_window(int x, int y, int radius, int shifted);

/*! Convolutes Image a with Image b at point i
 *
 *  The convolution only acts on the amplitudes. 
 *  This convolution, is not FFT based, and as such there is
 *  no wrap around effect. That is points of b which are outside
 *  of a are simply discarded.
 */
spimage_EXPORT real point_convolute_img(Image * a, Image * b, int i);


/*! Returns the variance of the image in respect to the vicinity defined by window
 *
 * This function convolutes img with the window (making sure there's no aliasing)
 * and then returns the difference between this averaged image and the input img
 */

spimage_EXPORT Image * image_local_variance(Image * img, Image * window);
/*@}*/


/** @defgroup Projectors Projector Operators
 *  Projector Operators for fixed point iterative algorithms
 *  @{
 */

/*! Module projector
 *
 *  Returns an image that combines the phases of a with the amplitudes of b
 *  Both images are assumed to be in reciprocal space
 */
spimage_EXPORT Image * sp_proj_module(Image * a, Image * b);

/*! Support projector
 *
 *  Sets to zero all regions of the image a for which image b is 0
 *  and keeps the image for regions for which image b is 1
 *  Both images are assumed to be in real space
 *
 *  Warning: Image b must only take values 0 and 1
 */
spimage_EXPORT Image * sp_proj_support(Image * a, Image * b);

/*@}*/
spimage_EXPORT void add_gaussian_noise(Image * in, real level);
spimage_EXPORT float box_muller(float m, float s);
spimage_EXPORT real lin_image_interpol(Image * img, real x, real y);
spimage_EXPORT real centro_sym_value(int index,Image * img);
spimage_EXPORT int centro_sym_index(int index,Image * img);
spimage_EXPORT Image * centro_sym_correlation(Image  * img);
spimage_EXPORT Image * shift_quadrants(Image * img);
spimage_EXPORT Image * make_shifted_image_square(Image * in);
spimage_EXPORT Image * make_unshifted_image_square(Image * in);
spimage_EXPORT Image * average_centrosymetry(Image * in);
spimage_EXPORT real r_factor(Image * fobs, Image *fcalc,real low_intensity_cutoff);
spimage_EXPORT Image * reflect_xy(Image * in, int in_place);
spimage_EXPORT Image * reflect_x(Image * in, int in_place);
spimage_EXPORT Image * reflect_y(Image * in, int in_place);
spimage_EXPORT Image * make_real(Image * in, int in_place);
spimage_EXPORT Image * centrosym_convolve(Image * in);
spimage_EXPORT Image * get_phase_angle_image(Image * img);
spimage_EXPORT Image * zero_pad_image(Image * a, int newx, int newy, int pad_mask);
spimage_EXPORT Image * shift_center_to_top_left(Image * a);
spimage_EXPORT real I_divergenge(Image * a, Image * b);
spimage_EXPORT Image * square_blur(Image * in, real radius);
spimage_EXPORT Image * patterson_function(Image * a);
spimage_EXPORT void find_center(Image * img, real * center_x, real * center_y);
spimage_EXPORT int pixel_to_index(Image * img, real * point);
spimage_EXPORT void smooth_edges(Image * img, int resolution);
spimage_EXPORT Image * fourier_rescale(Image * img, int x, int y);
spimage_EXPORT real bilinear_interpol_img(Image * img, real * data, real v_x, real v_y);
spimage_EXPORT Image * bilinear_rescale(Image * img, int new_x, int new_y);
spimage_EXPORT void resize_empty_image(Image * a,int newx, int newy);
spimage_EXPORT Image * get_mask_from_image(Image * a);
spimage_EXPORT int quadrant_shift_index(Image * a, int index);

spimage_EXPORT Image * normalize_image(Image * in);

spimage_EXPORT real p_drand48();

#endif
