#ifndef _FFT_H_
#define _FFT_H_

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

/** @defgroup FFT Fast Fourier Transform
 *  Forward and Backward fourier transform

 *  The forward transform is define as \f$ F(h) = \int_{-\infty}^{\infty}f(r) \exp(- 2 \pi i h \cdot r) {\textrm d} r \f$
 *  @{
 */

#ifndef M_PI
#define M_PI 3.1415926535897932846
#endif

#define FFTW3
/* #define FFTW2 */

#if defined FFTW3 && !defined __CUDACC__
#include <fftw3.h>
#define sp_image_ifft(a) sp_image_ifftw3(a)
#define sp_image_fft(a) sp_image_fftw3(a) 
#define sp_image_ifft_fast(a,b) sp_image_ifftw3_fast(a,b)
#define sp_image_fft_fast(a,b) sp_image_fftw3_fast(a,b) 
  //#define sp_cmatrix_ifft(a) sp_cmatrix_ifftw3(a)
#define sp_image_1d_fft(a,b) sp_image_1d_fftw3(a,b)
#define sp_image_1d_ifft(a,b) sp_image_1d_fftw3(a,b)
//#define sp_cmatrix_fft(a) sp_cmatrix_fftw3(a) 
#define sp_c3matrix_ifft(a) sp_c3matrix_ifftw3(a)
#define sp_c3matrix_fft(a) sp_c3matrix_fftw3(a)


#ifdef _SP_DOUBLE_PRECISION
  typedef fftw_complex fftwr_complex;
  typedef fftw_plan fftwr_plan;
  #define fftwr_execute(a) fftw_execute(a)
  #define fftwr_destroy_plan(a) fftw_destroy_plan(a)
  #define fftwr_malloc(a) fftw_malloc(a)
  #define fftwr_free(a) fftw_free(a)
  #define fftwr_plan_dft_2d(a,b,c,d,e,f) fftw_plan_dft_2d(a,b,c,d,e,f)
  #define fftwr_plan_dft_r2c_2d(a,b,c,d,e) fftw_plan_dft_r2c_2d(a,b,c,d,e)
  #define fftwr_plan_dft_3d(a,b,c,d,e,f,g) fftw_plan_dft_3d(a,b,c,d,e,f,g)
  #define fftwr_plan_dft_r2c_3d(a,b,c,d,e,f) fftw_plan_dft_r2c_3d(a,b,c,d,e,f)
  #define fftwr_plan_guru_split_dft(a,b,c,d,e,f,g,h,i) fftw_plan_guru_split_dft(a,b,c,d,e,f,g,h,i)
  #define fftwr_init_threads() fftw_init_threads()
  #define fftwr_plan_with_nthreads(a) fftw_plan_with_nthreads(a)
  #define fftwr_plan_many_dft(a,b,c,d,e,f,g,h,i,j,k,l,m) fftw_plan_many_dft(a,b,c,d,e,f,g,h,i,j,k,l,m)
#else
  typedef fftwf_complex fftwr_complex;
  typedef fftwf_plan fftwr_plan;
  #define fftwr_execute(a) fftwf_execute(a)
  #define fftwr_destroy_plan(a) fftwf_destroy_plan(a)
  #define fftwr_malloc(a) fftwf_malloc(a)
  #define fftwr_free(a) fftwf_free(a)
  #define fftwr_plan_dft_2d(a,b,c,d,e,f) fftwf_plan_dft_2d(a,b,c,d,e,f) 
  #define fftwr_plan_dft_r2c_2d(a,b,c,d,e) fftwf_plan_dft_r2c_2d(a,b,c,d,e) 
  #define fftwr_plan_dft_3d(a,b,c,d,e,f,g) fftwf_plan_dft_3d(a,b,c,d,e,f,g)
  #define fftwf_plan_dft_r2c_3d(a,b,c,d,e,f) fftwf_plan_dft_r2c_3d(a,b,c,d,e,f)
  #define fftwr_plan_guru_split_dft(a,b,c,d,e,f,g,h,i) fftwf_plan_guru_split_dft(a,b,c,d,e,f,g,h,i)
  #define fftwr_init_threads() fftwf_init_threads()
  #define fftwr_plan_with_nthreads(a) fftwf_plan_with_nthreads(a)
  #define fftwr_plan_many_dft(a,b,c,d,e,f,g,h,i,j,k,l,m) fftwf_plan_many_dft(a,b,c,d,e,f,g,h,i,j,k,l,m)
#endif

#endif
#ifdef FFTW2
#ifdef _SP_DOUBLE_PRECISION
#include <dfftw.h>
#else 
#include <sfftw.h>
#endif
#define image_guru_rev_fft(a) image_rev_fftw2(a)
#define image_rev_fft(a) image_rev_fftw2(a)
#define image_fft(a) image_fftw2(a) 
#define image_guru_fft(a) image_fftw2(a) 
#endif
  
/*! Returns the forward FFT of img.
 *
 * The mask is copied unchanged
 * The center of the image is set to the middle of the image
 */
spimage_EXPORT Image * sp_image_fft(const Image * img);

/*! Simply does the fft of the input data and puts it into out */
spimage_EXPORT void sp_image_fft_fast(const Image * in, Image * out);

  /*!Returns the 1d fourier transform of img along the specified axis. 
   * axis = 0,1,2 corresponds to x,y,z respectively. (axis = 1 only
   * works for 2d images.
   */
spimage_EXPORT Image * sp_image_1d_fft(const Image *img,int axis);

  /*!Returns the 1d back fourier transform of img along the specified axis. 
   * axis = 0,1,2 corresponds to x,y,z respectively. (axis = 1 only
   * works for 2d images.
   */
spimage_EXPORT Image * sp_image_1d_ifft(const Image *img,int axis);
/*! Returns the forward FFT of m.
 */
//spimage_EXPORT sp_cmatrix * sp_cmatrix_fft(sp_cmatrix * m);
/*! Returns the forward FFT of m.
 */
spimage_EXPORT sp_c3matrix * sp_c3matrix_fft(const sp_c3matrix * m);
/*! Returns the backward FFT of img.
 *
 * The mask is copied unchanged
 */
spimage_EXPORT Image * sp_image_ifft(const Image * img);
/*! Simply does the ifft of the input data and puts it into out */
spimage_EXPORT void sp_image_ifft_fast(const Image * in, Image * out);


/*! Returns the backward FFT of m.
 */
//spimage_EXPORT sp_cmatrix * sp_cmatrix_ifft(sp_cmatrix * img);
/*! Returns the backward FFT of m.
 */
spimage_EXPORT sp_c3matrix * sp_c3matrix_ifft(const sp_c3matrix * img);
/*! Initializes the fft routine and tells it to use 
 *   nthreads threads
 */
spimage_EXPORT int sp_init_fft(int nthreads);

#ifdef _USE_CUDA
  spimage_EXPORT Image * sp_image_cuda_ifft(const Image * img);
  spimage_EXPORT Image * sp_image_cuda_fft(const Image * img);
#endif

/*@}*/
#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif
