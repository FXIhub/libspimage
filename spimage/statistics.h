#ifndef _STATISTICS_H_
#define _STATISTICS_H_ 1

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */


/*! Returns the inverse of the cumulative distribution function of the normal distribution
 *
 * This function is taken from the GNU Scientific Library.
 */
spimage_EXPORT double embedded_gsl_cdf_gaussian_Pinv (const double P, const double sigma);


/*! Returns the inverse of the complement of the cumulative distribution function of the normal distribution
 *
 * This function is taken from the GNU Scientific Library.
 */
spimage_EXPORT double embedded_gsl_cdf_gaussian_Qinv (const double Q, const double sigma);



#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif
