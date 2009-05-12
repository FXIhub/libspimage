#ifndef _PRTF_H_
#define _PRTF_H_ 1

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */


/*! Calculates a prtf image from the list(of size n) of input images
 *  
 * The prtf is defined as the complex sum of the images divided
 * by the sum of the absolute value of the images
 *
 * The images are assumed to be in fourier space
 *
 * Keep in mind that this function does not try to correct for phase
 * shift or centrosymmetry!
 */
spimage_EXPORT Image * sp_prtf_basic(Image ** list, int n);


/*! Calculates a prtf image from the list(of size n) of input images
 * after trying to correct for phase shifts and centrosymmetry
 *
 * flags is the logical OR of:
 * 
 *  SpFourierSpace - Indicates that the input is in diffraction space
 *  SpRealSpace - Indicates that the input is in real space
 *  One of these two flags must be specified!
 *
 *  SpInPlace - Allows the function to actually modify the input list 
 *              so that they will end superimposed
 *  SpOutOfPlace - Do not modify input images
 *  One of these two flags must be specified!
 *
 * The prtf is defined as the complex sum of the images divided
 * by the sum of the absolute value of the images
 *
 */
spimage_EXPORT Image * sp_prtf_advanced(Image ** list,int n,int flags);


  /*! Calculates the radial average of the prtf image. It returns a list with the average prtf 
     with respect to the distance to the center in pixels.
   */
  sp_list * sp_prtf_by_resolution(Image * prtf);
#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif


