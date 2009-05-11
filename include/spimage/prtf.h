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
 * space indicates in which space the input images are
 *
 * The prtf is defined as the complex sum of the images divided
 * by the sum of the absolute value of the images
 *
 */
spimage_EXPORT Image * sp_prtf_advanced(Image ** list,int n,SpSpace space);
#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif


