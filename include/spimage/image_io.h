#ifndef _IMAGE_IO_H_
#define _IMAGE_IO_H_ 1

#include "image.h"
#include "linear_alg.h"

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

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
 *              SpColormapGrayScale    Image will use a grayscale color palete
 *              SpColormapTraditional  Image will use the traditional color map
 *              SpColormapHot          Image will use the hot color map
 *              SpColormapRainbow      Image will use the rainbow color map
 *              SpColormapJet          Image will use the jet color map
 *              SpColormapWheel        Image will use the color wheel color map
 *              SpColormapLogScale     Image will be written in log scale
 *
 */
spimage_EXPORT void sp_image_write(const Image * img,const char * filename, int flags);


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
spimage_EXPORT Image * _sp_image_read(const char * filename,int flags, const char * file, int line );
#ifdef SWIG
  Image * sp_image_read(const char * filename,int flags){
    return _sp_image_read(filename,flags,__FILE__,__LINE__);
  }
#else
#define sp_image_read(filename,flags) _sp_image_read(filename,flags,__FILE__,__LINE__)
#endif

/*! Write the mask of the image to a file in png format, using 
  the specified color map
*/
spimage_EXPORT int write_mask_to_png(const Image * img, char * filename, int color);
/*@}*/


#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif
