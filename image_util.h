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
real dist_to_corner(int i, Image * in);
/*! Calculates the distance from the pixel corresponding to index i to the image center,
  as defined in the Detector.  
 */
real dist_to_center(int i, Image * in);
/*! Calculates the distance from the pixel corresponding to index i to closest axis going
through the image center as defined in the Detector.  
 */
real dist_to_axis(int i,Image * in);
/*! Calculates the "Manhattan distance" from the pixel corresponding to index i to the image center,
  as defined in the Detector.  
 */
real square_dist_to_center(int i,Image * in);
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

Image * read_imagefile(const char * filename);

/*! Write an image from to a file in spimage format.
 *
 * It creates an spimage file with all the informations contained
 * in the Image structure with the specified precision.
 * output_precision must be 4 or 8 and specifies the number of bytes
 * used for storing each floating point number.
 */
void write_img(Image * img,char * filename,int output_precision);

/*! Create an image using the values from the TIFF file
*/
Image * read_tiff(char * filename);
/*! Write a TIFF using the values from the image file
*/
void write_tiff(Image * img, char * filename);

/*! Write the norm of the image to a file in png format, using 
  the specified color map
*/
int write_png(Image * img, char * filename, int color);

/*! Write the norm of the image to a file in VTK ASCII format as a 2D Structured grid
 */
int write_vtk(Image * img, char * filename);

/*! Write the mask of the image to a file in png format, using 
  the specified color map
*/
int write_mask_to_png(Image * img, char * filename, int color);
/*@}*/



/** @defgroup MM Memory Management
 *  Create, copy and free images
 *  @{
 */

/*! Create a copy of Image in
 */
Image * imgcpy(Image * in);
/*! Delete image Image in
 */
void freeimg(Image * in);

/*! Creates an Image with the same dimensions as a.
 *  The value of the pixels is all set to 0
 */
Image * create_empty_img(Image * a);

/*! Creates an empty image with given dimensions. Mask is set to 0 everywhere.
 *  Image in unscaled and unphased.
 */
Image * create_new_img(int x, int y);
/*@}*/

/** @defgroup ImageCrop Image Cropping
 *  Crops parts of the image
 *  @{
 */

/*! Masks out all pixels whose distance to center
  is smaller than radius. It also puts their values to 0.
 */
void remove_lowres(Image * in, real radius);
/*! Masks out all pixels whose "manhattan distance" to center
  is bigger than radius. It also puts their values to 0.
 */
Image * limit_resolution(Image * img, int resolution);
/*! Creates a new image from the rectange define by the upper left corner (x1,y1)
  and lower right corner (x2,y2).
 */
Image * rectangle_crop(Image * in, int x1, int y1, int x2, int y2); 
/*@}*/


/** @defgroup ImageAlgebra Image Algebra
 *  Basic algebra with images
 *  @{
 */
/*! Add Image b to Image a
 */
void add_image(Image * a, Image * b);
/*! Subtract Image b from Image a
 */
void sub_image(Image * a, Image * b);
/*! Returns the phase of a complexed value pixel.
 */
real get_phase_angle(Image * img, int i);

/*! Calculates the \f$ \cos() \f$ of the complexed value pixel i.
 *  Can also be understood as \f$ \frac{Re(f(i))}{|f(i)|} \f$
 */
real real_factor(Image * img, int i);
/*! Calculates the \f$ \sin() \f$ of the complexed value pixel i.
 *  Can also be understood as \f$ \frac{Im(f(i))}{|f(i)|} \f$
 */
real complex_factor(Image * img, int i);
/*! Calculates the norm of the complexed value pixel i.
 */
real norm(Image * img, int i);
/*! Calculates the integral of all the image, independently of the mask.
 */
real integrated_intensity(Image * a);

/*! Turns a real valued image into a complexed valued one,
  mantaining the real part and setting the imaginary part to 0.
 */
void rephase(Image *  img);

/*! Turns a complexed valued image into a real valued one, by
  taking the norm of each pixel.
 */
void dephase(Image *  img);

/*! Turns a real valued image into a complexed valued one,
  mantaining the norm and using a random phase.
 */
void random_rephase(Image * img);

/*! Multiply Image img by a scalar value.
 */
void multiply_scalar_by_image(Image * img, real value);
/*! Returns the index of the image maximum, and set *max to the
  value of that pixel
 */
int get_image_maximum(Image * img, int * x, int * y, real * max);

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
Image * gaussian_blur(Image * in, real radius);
/*! Multiplies the image with a gaussian centered in the image center of a given radius 
 *
 *   The mask will not be multiplied.
 */
Image * gaussian_filter(Image * in, real radius, int in_place);

/*! Convolutes Image a with Image b 
 *
 *  Note that b must fit inside a or the other way around.
 *  The convolution only acts on the amplitudes, 
 *  the real and complex parts of the Image will be set to 0.
 *  The resulting mask will be the same as the one from image a.
 */
Image * convolute_img(Image * a, Image * b);

/*! Low pass filter using a centered square window of side edge_size 
 *
 *  The mask is also low pass filtered.
 */

Image * low_pass_square_filter(Image * in, int edge_size);
/*! Low pass filter using a centered gaussian window of side edge_size 
 *
 * The mask is also low pass filtered.
 */
Image * low_pass_gaussian_filter(Image * in, int edge_size);

/*! Correlates Image a with Image b 
 *
 *  Note that b must fit inside a or the other way around.
 *  The correlation only acts on the amplitudes, 
 *  the real and complex parts of the Image will be set to 0.
 *  The resulting mask will be the same as the one from image a.
 */
Image * cross_correlate_img(Image * a, Image * b);

/*! Returns a rectangular window centered on the image center
 * with the given width and height 
 */
Image * rectangular_window(Image * a, int width, int height);


/*! Returns a circular window centered on the image center
 * with the given radius
 */
Image * circular_window(Image * a, int radius);

/*! Convolutes Image a with Image b at point i
 *
 *  The convolution only acts on the amplitudes. 
 *  This convolution, is not FFT based, and as such there is
 *  no wrap around effect. That is points of b which are outside
 *  of a are simply discarded.
 */
real point_convolute_img(Image * a, Image * b, int i);

/*@}*/

void add_gaussian_noise(Image * in, real level);
float box_muller(float m, float s);
real lin_image_interpol(Image * img, real x, real y);
real centro_sym_value(int index,Image * img);
int centro_sym_index(int index,Image * img);
Image * centro_sym_correlation(Image  * img);
Image * shift_quadrants(Image * img);
Image * make_shifted_image_square(Image * in);
Image * make_unshifted_image_square(Image * in);
Image * average_centrosymetry(Image * in);
real r_factor(Image * fobs, Image *fcalc,real low_intensity_cutoff);
Image * reflect_xy(Image * in, int in_place);
Image * reflect_x(Image * in, int in_place);
Image * reflect_y(Image * in, int in_place);
Image * make_real(Image * in, int in_place);
Image * centrosym_convolve(Image * in);
Image * get_phase_image(Image * img);
Image * get_phase_angle_image(Image * img);
Image * zero_pad_image(Image * a, int newx, int newy, int pad_mask);
Image * shift_center_to_top_left(Image * a);
real I_divergenge(Image * a, Image * b);
Image * square_blur(Image * in, real radius);
Image * patterson_function(Image * a);
void find_center(Image * img, real * center_x, real * center_y);
int pixel_to_index(Image * img, real * point);
void smooth_edges(Image * img, int resolution);
Image * fourier_rescale(Image * img, int x, int y);
real bilinear_interpol_img(Image * img, real * data, real v_x, real v_y);
Image * bilinear_rescale(Image * img, int new_x, int new_y);
void resize_empty_image(Image * a,int newx, int newy);
Image * get_mask_from_image(Image * a);
int quadrant_shift_index(Image * a, int index);
#endif
