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
 * \param img The image to be output
 * \param filename The name of the file to be created
 * \param flags describe the colormap, data precision, iteration number, and type of functionality
 *
 * The file type is infered from the filename extension.
 * Supported extensions are:
 *
 *   - .h5 - Create an hdf5 file
 *   - .cxi - Create a CXIDB file
 *   - .tiff or .tif - Create a TIFF file
 *   - .png - Create a png file
 *   - .vtk - Create a VTK file
 *   - .csv - Create a CSV file
 *
 * It creates an spimage file with all the informations contained
 * in the Image structure with the specified precision.
 *
 * The meaning of the flags depends on the type of file created:
 * There are only flags for .h5, .cxi and .png files, for the others this
 * value is ignored. The last 32 bits were added to flags for .cxi support
 * and contain the iteration number of the reconstructed image.
 * 
 <table>
 <tr> <td>File type</td>    <td>Flag Value</td><td>          Meaning</td></tr>
 *
 *<tr><td>       .h5</td></tr>
 *<tr><td></td><td>              sizeof(float)</td><td>      Data should be written in single precision</td></tr>
 *<tr><td></td><td>              sizeof(double)</td><td>     Data should be written in double precision</td></tr>
 *<tr><td>       .cxi</td></tr>
 *<tr><td></td><td>              create</td><td>             Data should be written to new file</td></tr>
 *<tr><td></td><td>              append entry</td><td>       Data should be appended to existing file as new entry</td></tr>
 *<tr><td></td><td>              append image</td><td>       Data should be appended to existing entry as new image</td></tr>
 *<tr><td></td><td>              append data</td><td>        Data should be appended to existing image</td></tr>
 *<tr><td></td><td>              iteration number</td><td>   The last 32 bits are interpreted as the iteration number of the image</td></tr> 
 *<tr><td>       .png</td></tr>
 *<tr><td></td><td>              SpColormapGrayScale    </td><td>Image will use a grayscale color palete</td></tr>
 *<tr><td></td><td>              SpColormapTraditional  </td><td>Image will use the traditional color map</td></tr>
 *<tr><td></td><td>              SpColormapHot          </td><td>Image will use the hot color map</td></tr>
 *<tr><td></td><td>              SpColormapRainbow      </td><td>Image will use the rainbow color map</td></tr>
 *<tr><td></td><td>              SpColormapJet          </td><td>Image will use the jet color map</td></tr>
 *<tr><td></td><td>              SpColormapWheel        </td><td>Image will use the color wheel color map</td></tr>
 *<tr><td></td><td>              SpColormapLogScale     </td><td>Image will be written in log scale</td></tr>
 *  </table>
 */
spimage_EXPORT void sp_image_write(const Image *img, const char *filename, long long flags);


/*! Reads an image from the specified filename
 * 
 *
 * \param filename The name of the file to be read
 * \param flags currently unused
 * \return The image that was read, or NULL if there was an error
 * The file type is infered from the filename extension.
 * Supported extensions are:
 *
 * .h5 - Reads a hdf5 file
 * .cxi - Reads a CXIDB file
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
  
 * \param img The image to be output
 * \param filename The name of the file to be created
 * \param color describe the colormap
*/
spimage_EXPORT int write_mask_to_png(const Image * img, char * filename, int color);
/*@}*/


#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif
