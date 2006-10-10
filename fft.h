#ifndef _FFT_H_
#define _FFT_H_


#define FFTW3
/* #define FFTW2 */
/* #define DOUBLE */

#ifdef FFTW3
#include <fftw3.h>
#define image_guru_rev_fft(a) image_rev_fftw3(a)
#define image_rev_fft(a) image_rev_fftw3(a)
#define image_fft(a) image_fftw3(a) 
#define image_guru_fft(a) image_fftw3(a)

#ifdef DOUBLE
  typedef fftw_complex fftwr_complex;
  typedef fftw_plan fftwr_plan;
  #define fftwr_execute(a) fftw_execute(a)
  #define fftwr_destroy_plan(a) fftw_destroy_plan(a)
  #define fftwr_malloc(a) fftw_malloc(a)
  #define fftwr_free(a) fftw_free(a)
  #define fftwr_plan_dft_2d(a,b,c,d,e,f) fftw_plan_dft_2d(a,b,c,d,e,f)
  #define fftwr_plan_dft_r2c_2d(a,b,c,d,e) fftw_plan_dft_r2c_2d(a,b,c,d,e)
  #define fftwr_plan_guru_split_dft(a,b,c,d,e,f,g,h,i) fftw_plan_guru_split_dft(a,b,c,d,e,f,g,h,i)
#else
  typedef fftwf_complex fftwr_complex;
  typedef fftwf_plan fftwr_plan;
  #define fftwr_execute(a) fftwf_execute(a)
  #define fftwr_destroy_plan(a) fftwf_destroy_plan(a)
  #define fftwr_malloc(a) fftwf_malloc(a)
  #define fftwr_free(a) fftwf_free(a)
  #define fftwr_plan_dft_2d(a,b,c,d,e,f) fftwf_plan_dft_2d(a,b,c,d,e,f) 
  #define fftwr_plan_dft_r2c_2d(a,b,c,d,e) fftwf_plan_dft_r2c_2d(a,b,c,d,e) 
  #define fftwr_plan_guru_split_dft(a,b,c,d,e,f,g,h,i) fftwf_plan_guru_split_dft(a,b,c,d,e,f,g,h,i)
#endif

#endif
#ifdef FFTW2
#ifdef DOUBLE
#include <dfftw.h>
#else 
#include <sfftw.h>
#endif
#define image_guru_rev_fft(a) image_rev_fftw2(a)
#define image_rev_fft(a) image_rev_fftw2(a)
#define image_fft(a) image_fftw2(a) 
#define image_guru_fft(a) image_fftw2(a) 
#endif

Image * image_fft(Image * img);
Image * image_rev_fft(Image * img);
Image * real_image_fft(Image * img);
Image * image_guru_fft(Image * img);
Image * image_guru_rev_fft(Image * img);

void create_fft_bench_table(int iter, FILE * f);

#endif
