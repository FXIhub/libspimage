#ifndef _IMAGE_NOISE_H_
#define _IMAGE_NOISE_H_ 1

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

/*! This function tries to determine the standard deviation of each pixel using
 *  the constraint that the autocorrelation of the intensities is space limited.
 *
 *  This function calculates the discrepancy of the intensities and the intensities
 *  convoluted with the fourier transform of the autocorrelation support
 *  If the image was noiseless this 2 images should be the same, but in the presence of
 *  noise they will defer. I should write some paper about it to describe this in more
 *  detail.
 */
spimage_EXPORT Image * sp_image_noise_estimate(Image * intensities, Image * autocorrelation_support);


#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif
