#ifndef _IMAGE_UTIL_H_
#define _IMAGE_UTIL_H_

#include "image.h"
#include "linear_alg.h"

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */


#define TSIZE(a) (a->detector->size[0]*a->detector->size[1])

#define COLOR_GRAYSCALE 1
#define COLOR_TRADITIONAL 2
#define COLOR_HOT 4
#define COLOR_RAINBOW 8
#define COLOR_JET 16
#define LOG_SCALE 32
#define COLOR_PHASE 64

#define OUT_OF_PLACE 0
#define IN_PLACE 1

#ifndef MAX
#define MAX(a,b) ((a)>(b)?(a):(b))
#endif

#ifndef MIN
#define MIN(a,b) ((a)<(b)?(a):(b))
#endif

#define SP_ZERO_PHASE 1
#define SP_RANDOM_PHASE 2
#define SP_GAUSSIAN 3

#define SP_TO_AXIS 4
#define SP_TO_CENTER 5
#define SP_TO_CENTER2 6
#define SP_TO_CORNER 7

#define SP_COPY_DATA 1
#define SP_COPY_DETECTOR 0
#define SP_COPY_MASK 2

#define SP_AXIS_XY 0
#define SP_AXIS_X 1
#define SP_AXIS_Y 2
#define SP_ORIGO 3

#define SP_ZERO_PAD_EDGE 1
#define SP_SYMMETRIC_EDGE 2
#define SP_REPLICATE_EDGE 3
#define SP_CIRCULAR_EDGE 4

#define SP_TRANSLATE_WRAP_AROUND 1
#define SP_TRANSLATE_DISCARD_OUTSIDE 2

#define SP_ENANTIOMORPH 1

/** @defgroup Distance
 *  Calculates several kinds of distances in an image
 *  @{
 */


/*! Calculates the distance from the pixel corresponding to index i
 *  
 * The type flag determines to which point and how the distance is calculated.
 * The possible values for the flag is:
 *
 * SP_TO_AXIS - Distance to the closest axis going through the middle of the image.
 *
 * SP_TO_CENTER - Distance to center of the image.
 *
 * SP_TO_CENTER2 - Same as SP_TO_CENTER squared.
 *
 * SP_TO_CORNER - Distance to the closest corner.
 *
 */
spimage_EXPORT real sp_image_dist(Image * in, int i, int type);
/*@}*/



/** @defgroup IO Input-Output
 *  Image input output routines
 *  @{
 */

/*! Write the given file to the specified filename
 * 
 * The file type is infered from the filename extension.
 * Supported extensions are:
 *
 * .h5 - Create an hdf5 file
 * .tiff or .tif - Create a TIFF file
 * .png - Create a png file
 * .vtk - Create a VTK file
 * .csv - Create a CSV file
 *
 * It creates an spimage file with all the informations contained
 * in the Image structure with the specified precision.
 *
 * The meaning of the flags depends on the type of file created:
 * There are only flags for .h5 and .png files, for the others this
 * value is ignored 
 * 
 *   File type    Flag Value          Meaning
 *
 *      .h5
 *              sizeof(float)      Data should be written in single precision
 *              sizeof(double)     Data should be written in double precision
 *
 *      .png
 *              COLOR_GRAYSCALE    Image will use a grayscale color palete
 *              COLOR_TRADITIONAL  Image will use the traditional color map
 *              COLOR_HOT          Image will use the hot color map
 *              COLOR_RAINBOW      Image will use the rainbow color map
 *              COLOR_JET          Image will use the jet color map
 *              LOG_SCALE          Image will be written in log scale
 *
 */
spimage_EXPORT void sp_image_write(Image * img,const char * filename, int flags);


/*! Reads an image from the specified filename
 * 
 * The file type is infered from the filename extension.
 * Supported extensions are:
 *
 * .h5 - Reads a hdf5 file
 * .tiff or .tif - Reads a TIFF file
 * .png - Reads a png file
 *
 *
 * The flags value is currently ignored.
 *
 */
spimage_EXPORT Image * _sp_image_read(const char * filename,int flags, char * file, int line );
#define sp_image_read(filename,flags) _sp_image_read(filename,flags,__FILE__,__LINE__)

/*! Write the mask of the image to a file in png format, using 
  the specified color map
*/
spimage_EXPORT int write_mask_to_png(Image * img, char * filename, int color);
/*@}*/



/** @defgroup MM Memory Management
 *  Create, copy and free images
 *  @{
 */

/*! Allocates an image of width x and height y
 */
spimage_EXPORT Image * _sp_image_alloc(int x, int y, int z, char * file, int line);
#define sp_image_alloc(x,y,z) _sp_image_alloc(x,y,z,__FILE__,__LINE__)

/*! Create a copy of Image in
 *
 * Possible values for flags are a bitwise combination of the following options
 *
 * SP_COPY_DATA - memcpy the image data from in to the newly created image
 * SP_COPY_MASK - memcpy the mask from in to the newly created image
 */
spimage_EXPORT Image * _sp_image_duplicate(Image * in,int flags, char * file, int line);
#define sp_image_duplicate(in, flags) _sp_image_duplicate(in,flags,__FILE__,__LINE__)

/*! Delete image Image in
 */
spimage_EXPORT void _sp_image_free(Image * in, char * file, int line);
#define sp_image_free(in) _sp_image_free(in,__FILE__,__LINE__)

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
spimage_EXPORT void sp_image_high_pass(Image * in, real radius, int type);
/*! Masks out all pixels whose "manhattan distance" to center
  is bigger than radius. It also puts their values to 0.
 */
spimage_EXPORT Image * sp_image_low_pass(Image * img, int resolution, int type);
/*! Creates a new image from the rectange define by the upper left corner (x1,y1)
  and lower right corner (x2,y2).
 */
spimage_EXPORT Image * rectangle_crop(Image * in, int x1, int y1, int x2, int y2); 
spimage_EXPORT Image * cube_crop(Image * in, int x1, int y1, int z1, int x2, int y2, int z2);

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
spimage_EXPORT real sp_image_distance_to_edge(Image * img, real * point, real direction, real * intersection);

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

/*! Sets the point x,y,z (or x,y,0 for 2D) to the value v
 */
static inline void sp_image_set(Image * a,int x,int y,int z,Complex v){
  sp_c3matrix_set(a->image,x,y,z,v);
}

/*! Returns the image value at point x,y,z (or x,y,0 for 2D)
 */
static inline Complex sp_image_get(Image * a,int x,int y,int z){
  return sp_c3matrix_get(a->image,x,y,z);
}

/*! Returns x*y*z size of image a
 */
static inline long long sp_image_size(Image * a){
  return sp_c3matrix_size(a->image);
}

/*! Returns the x length of image a
 */
static inline int sp_image_x(const Image * a){
  return sp_c3matrix_x(a->image);
}

/*! Returns the y length of image a
 */
static inline int sp_image_y(const Image * a){
  return sp_c3matrix_y(a->image);
}

/*! Returns the z length of image a
 */
static inline int sp_image_z(const Image * a){
  return sp_c3matrix_z(a->image);
}

/*! Returns the index of point x,y,z (or x,y,0)  of image a
 */
static inline  long long sp_image_get_index(const Image * a, int x, int y, int z){
  return sp_c3matrix_get_index(a->image,x,y,z);
}


/*! Fill the entire image with value
 */
spimage_EXPORT void sp_image_fill(Image * a, Complex value);

/*! Add Image b to Image a
 */
spimage_EXPORT void sp_image_add(Image * a, Image * b);
/*! Subtract Image b from Image a
 */
spimage_EXPORT void sp_image_sub(Image * a, Image * b);

/*! Mutiplies Image a with Image b element by element
 */
spimage_EXPORT void sp_image_mul_elements(Image * a, Image * b);

/*! Transforms an image in its complex conjugate
 */
spimage_EXPORT void sp_image_conj(Image * a);

/*! Sets a to 1/a. 
 */
spimage_EXPORT int sp_image_invert(Image * a);


/*! Calculates the dot product of a and b as if they were vectors 
 *
 */
spimage_EXPORT Complex sp_image_dot_prod(Image * a, Image * b);



spimage_EXPORT real integrated_intensity(Image * a);

/*! Turns a real valued image into a complexed valued one.
 *
 * The possible values for type are:
 *
 * SP_ZERO_PHASE - mantain the real part and setting the imaginary part to 0.
 * 
 * SP_RANDOM_PHASE - The phase of the result is randomly distributed and 
 * the absolute value is equal to the real value of the input.
 */
spimage_EXPORT void sp_image_rephase(Image *  img, int type);

/*! Turns a complexed valued image into a real valued one, by
  taking the norm of each pixel.
 */
spimage_EXPORT void sp_image_dephase(Image *  img);


/*! Transform image into intensities
 *
 * Multiplies the image by the complex conjugate of itself if the image is scaled.
 * If the image is not scaled it does nothing.
 */
static inline void sp_image_to_intensities(Image *  img){
  if(img->scaled){
    for(long long i = 0;i<sp_image_size(img);i++){
      img->image->data[i] = sp_cmul(img->image->data[i],sp_cconj(img->image->data[i]));
    }
    img->scaled = 0;
  }
}

/*! Transform image into amplitudes
 *
 * Returns the square root of each pixel of the image if the image is not scaled.
 * If the image is already scaled it does nothing.
 */
static inline  void sp_image_to_amplitudes(Image *  img){
  if(!img->scaled){
    for(long long i = 0;i<sp_image_size(img);i++){
      sp_real(img->image->data[i]) = sqrt(sp_real(img->image->data[i]));
      sp_imag(img->image->data[i]) = 0;
    }
    img->scaled = 1;
  }
}


/*! Multiply Image img by a scalar value.
 */
spimage_EXPORT void sp_image_scale(Image * img, real value);
/*! Returns the maximum value of the image and sets *index to the
  index of that pixel along iwth *x and *y
 */
spimage_EXPORT real sp_image_max(Image * img, long long * index, int * x, int * y, int * z);

/*! Returns the phase of the image on each pixel, in radians
 *
 * Returns the phase of a in radians, from [0, 2*M_PI[ or NULL if there's an error
 */
spimage_EXPORT Image * sp_image_get_phases(Image * a);

/*! Transposes the images
*/
spimage_EXPORT void sp_image_tranpose(Image * in);

/*! Returns the linearly interpolated value of the image at v_x, v_y
*/
spimage_EXPORT real sp_image_interp(Image * img, real v_x, real v_y, real v_z);

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
spimage_EXPORT Image * sp_image_convolute(Image * a, Image * b, int * size);

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
spimage_EXPORT void sp_add_noise(Image * in, real level,int type);
spimage_EXPORT float box_muller(float m, float s);
spimage_EXPORT real lin_image_interpol(Image * img, real x, real y);
spimage_EXPORT real centro_sym_value(long long index,Image * img);
spimage_EXPORT int centro_sym_index(long long index,Image * img);
spimage_EXPORT Image * centro_sym_correlation(Image  * img);
spimage_EXPORT Image * sp_image_shift(Image * img);
spimage_EXPORT Image * make_shifted_image_square(Image * in);
spimage_EXPORT Image * make_unshifted_image_square(Image * in);
spimage_EXPORT Image * average_centrosymetry(Image * in);
spimage_EXPORT real r_factor(Image * fobs, Image *fcalc,real low_intensity_cutoff);
spimage_EXPORT Image * sp_image_reflect(Image * in, int in_place,int axis);
spimage_EXPORT Image * centrosym_convolve(Image * in);
spimage_EXPORT Image * get_phase_angle_image(Image * img);
spimage_EXPORT Image * zero_pad_image(Image * a, int newx, int newy, int newz, int pad_mask);
spimage_EXPORT Image * shift_center_to_top_left(Image * a);
spimage_EXPORT real I_divergenge(Image * a, Image * b);
spimage_EXPORT Image * square_blur(Image * in, real radius, int type);
spimage_EXPORT Image * sp_image_patterson(Image * a);
spimage_EXPORT void find_center(Image * img, real * center_x, real * center_y, real * center_z);
spimage_EXPORT long long pixel_to_index(Image * img, real * point);
spimage_EXPORT void sp_image_smooth_edges(Image * img, sp_i3matrix * mask, int flags, real * value);
spimage_EXPORT Image * fourier_rescale(Image * img, int x, int y, int z);
spimage_EXPORT real bilinear_interpol_img(Image * img, real * data, real v_x, real v_y);
spimage_EXPORT Image * bilinear_rescale(Image * img, int new_x, int new_y, int new_z);
spimage_EXPORT void resize_empty_image(Image * a,int newx, int newy);
spimage_EXPORT Image * get_mask_from_image(Image * a);
spimage_EXPORT int sp_image_shift_index(Image * a, long long index);
spimage_EXPORT void sp_image_normalize(Image * in);

spimage_EXPORT real p_drand48();

/*! Resizes the input image to the desired size. It's content is discarded
 *
 */
spimage_EXPORT void _sp_image_realloc(Image * img, int new_x, int new_y, int new_z, char * file, int line);
#define sp_image_realloc(img,x,y,z) _sp_image_realloc(img,x,y,z,__FILE__,__LINE__)

/*! Filters the input image with a median filter
 *
 * The kernels tells how big the window is and how to weight the pixels.
 * The center of the center is equal to its dimensions/2.
 * The edge_flags correspond to the sp_image_edge_extend flags().
 */
spimage_EXPORT void sp_image_median_filter(Image * a,sp_i3matrix * kernel, int edge_flags, int type);


/*! Extend an image outside its borders by radius pixels on each
 *  direction.
 *
 *
 * Shifted images are padded in the middle.
 * The edge_flags determines the kind of extension done.
 *
 * Possible values are:
 *
 * SP_ZERO_PAD_EDGE - The values outside the border are set to 0
 * SP_SYMMETRIC_EDGE - The values outside the bounds of the image
 * are computed by mirror-reflecting the image across the image border.
 * SP_REPLICATE_EDGE - The values outside the bounds of the image are 
 * assumed to equal the nearest image border value.
 * SP_CIRCULAR_EDGE - The values outside the bounds of the image are
 * computed by implicitly assuming the input image is periodic.
 */
spimage_EXPORT Image * sp_image_edge_extend(Image * a, int radius, int edge_flags, int type);

/*! Inserts image from into image to at the position at_x, at_y
 *
 *  If from doesn't fit in to, the image is silently clipped.
 */
spimage_EXPORT void sp_image_insert(Image * to, Image * from, int at_x, int at_y, int at_z);

/*! Implements the beloved bubble sort for a list of reals
 *
 * The array will be sorted in ascending order
 */
spimage_EXPORT void sp_bubble_sort(real * a, int n);

/*! Contrast stretches an image
 *
 * The image is partitioned in x_div*y_div divisions
 * and each division is contrast stretched 
 * The scaling is linear interpolated from division to division
 */
spimage_EXPORT void sp_image_adaptative_constrast_stretch(Image * a,int x_div, int y_div);


/*! Calculates the fourier space coordinates of each pixel of the image 
 *
 * The matrices k_x k_y and k_z are filled with the coordinates.
 * If any of these matrices is NULL then it's not filed in .
 * These matrices should be prealocated.
 *
 * To access the x fourier coordinate at pixel_x = 1 pixel_y = 2
 * the following should be used: sp_matrix_get(k_x,1,2) and not the other way around!
 *
 * The wavelength distance to detector and pixel size of the image must be correctly set for
 * this function to work!
 */
spimage_EXPORT void sp_image_fourier_coords(Image * in, sp_3matrix * k_x, sp_3matrix * k_y, sp_3matrix * k_z);

/*! Initialize the random number generator with a given seed
 *
 */
spimage_EXPORT void sp_srand(int i);


/*! Superimposes image b on top of image a 
 *
 *  flags is a bitwise combination of the following:
 *
 *  SP_ENANTIOMORPH - allow to try to superimpose not only b but also
 *  the "mirror image" of b [b(-x)].
 *
 */
spimage_EXPORT void sp_image_superimpose(Image * a,Image * b, int flags);


/*! Translates an image by x,y,z.
 * 
 *  What happens to the data that falls outside of the boundaries of the image depends on the flag.
 *
 *  SP_TRANSLATE_WRAP_AROUND - The data is brought back inside the image assuming a periodic image.
 *  SP_TRANSLATE_DISCARD_OUTSIDE - The data that falls outside is simply discarded.
 */
spimage_EXPORT void sp_image_translate(Image * a,int x,int y,int z, int flags);


/*! Calculates the Real Space R factor between a and b
 *
 *  This routine does not scale or translate the images! This must be done before hand
 *  
 *  The formula used is Sum(|Fa|-|Fb|)/Sum(|Fa|+|Fb|)
 */
spimage_EXPORT real sp_image_rs_r_factor(Image *a, Image *b);


/*! Calculates the correlation coefficient of image a and b
 * 
 */
spimage_EXPORT real sp_image_correlation_coefficient(Image * a,Image * b);


spimage_EXPORT sp_vector * sp_image_center_of_mass(Image * a);

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif
