#include <stdlib.h>
/* #include <sys/time.h>*/
#include <time.h>
#include <limits.h>
#include "spimage.h"



real * sp_image_fft_shift(real * fftout, Image * a){
  int x,y;
  int nx = sp_cmatrix_cols(a->image);
  int ny = sp_cmatrix_rows(a->image);
  
  real * img_d = malloc(sizeof(real)*(ny/2+1)*nx);
  for(x = 0;x<nx;x++){
    for(y = 0;y<ny/2+1;y++){
      img_d[y+x*(ny/2+1)] = fftout[y+((nx/2+x)%nx)*(ny/2+1)];
    }
  }
  return img_d;
}


int sp_init_fft(int nthreads){
  int ret = 0; 
  if(nthreads == 1){
    /* No need to do anything*/
    return 0;
  }
#ifdef FFTW3
  ret = fftwr_init_threads();  
  if(!ret){
    perror("Error initializing parallel fftw!\n");
  }else{
    fftwr_plan_with_nthreads(nthreads);
  }  
#endif
  return ret;
}

#ifdef FFTW3

Image * sp_image_ifftw3(Image * img){
  fftwr_complex *out; 
  fftwr_complex *in; 
  fftwr_plan plan;
  Image * res;
  if(!img->phased){
    fprintf(stderr,"Error: Trying reverse fft an unphased image!\n");
    abort();
  }
  if(!img->shifted){
    fprintf(stderr,"Error: Trying to rev_fft an unshifted image!\n");
    abort();
  }

  res = sp_image_duplicate(img,SP_COPY_DETECTOR);
  /* We're gonna rely on the binary compatibility between fftwr_complex and Complex type */
  /*
  in = (fftwr_complex*) fftwr_malloc(sizeof(fftwr_complex) * (sp_cmatrix_rows(img->image))*(sp_cmatrix_cols(img->image)));  
  out = (fftwr_complex*) fftwr_malloc(sizeof(fftwr_complex) * (sp_cmatrix_rows(img->image))*(sp_cmatrix_cols(img->image)));
  */
  in = img->image->data;
  out = res->image->data;
  plan = fftwr_plan_dft_2d(sp_cmatrix_cols(img->image),sp_cmatrix_rows(img->image),in,out, FFTW_BACKWARD,FFTW_ESTIMATE);

  fftwr_execute(plan);
  fftwr_destroy_plan(plan);
  res->shifted = 0;
  return res;
}


sp_cmatrix * sp_cmatrix_ifftw3(sp_cmatrix * m){
  fftwr_complex *out; 
  fftwr_complex *in; 
  fftwr_plan plan;
  sp_cmatrix * res;
  res = sp_cmatrix_alloc(sp_cmatrix_rows(m),sp_cmatrix_cols(m));
  in = m->data;
  out = res->data;
  plan = fftwr_plan_dft_2d(sp_cmatrix_cols(m),sp_cmatrix_rows(m),in,out, FFTW_BACKWARD,FFTW_ESTIMATE);

  fftwr_execute(plan);
  fftwr_destroy_plan(plan);
  return res;
}

#endif

#ifdef FFTW2

Image * sp_image_ifftw2(Image * img){
  fftw_complex *out; 
  fftw_complex *in; 
  fftwnd_plan plan;
  int i,j;
  if(!img->phased){
    fprintf(stderr,"Error: Trying reverse fft an unphased image!\n");
    abort();
  }
  if(!img->shifted){
    fprintf(stderr,"Error: Trying to rev_fft an unshifted image!\n");
    abort();
  }

  Image * res = create_empty_img(img);
  /* We're gonna rely on the binary compatibility between fftwr_complex and Complex type */
  in = img->image->data;
  out = res->image->data;
  plan = fftw2d_create_plan(sp_cmatrix_cols(img->image),sp_cmatrix_rows(img->image), FFTW_BACKWARD,FFTW_ESTIMATE);

  j = 0;
  fftwnd_one(plan,in,out);
  fftwnd_destroy_plan(plan);
  res->shifted = 0;
  return res;
}
#endif


#ifdef FFTW3

Image * sp_image_fftw3(Image * img){
  fftwr_complex *out; 
  fftwr_complex *in; 
  fftwr_plan plan;
  Image * res = sp_image_duplicate(img,SP_COPY_DETECTOR);

  sp_image_rephase(res,SP_ZERO_PHASE);
  /* Rely on binary compatibility of the Complex type */
  in = img->image->data;
  out = res->image->data;
  plan = fftwr_plan_dft_2d(sp_cmatrix_cols(img->image),sp_cmatrix_rows(img->image),in,out,FFTW_FORWARD,FFTW_ESTIMATE);

  fftwr_execute(plan);
  fftwr_destroy_plan(plan);
  res->shifted = 1;
  res->detector->image_center[0] = (sp_cmatrix_cols(res->image)-1)/2.0;
  res->detector->image_center[1] = (sp_cmatrix_rows(res->image)-1)/2.0;
  return res;
}

sp_cmatrix * sp_cmatrix_fftw3(sp_cmatrix * m){
  fftwr_complex *out; 
  fftwr_complex *in; 
  fftwr_plan plan;
  sp_cmatrix * res = sp_cmatrix_alloc(sp_cmatrix_rows(m),sp_cmatrix_cols(m));

  /* Rely on binary compatibility of the Complex type */
  in = m->data;
  out = res->data;
  plan = fftwr_plan_dft_2d(sp_cmatrix_cols(m),sp_cmatrix_rows(m),in,out,FFTW_FORWARD,FFTW_ESTIMATE);

  fftwr_execute(plan);
  fftwr_destroy_plan(plan);
  return res;
}

#endif

#ifdef FFTW2

sp_cmatrix * sp_cmatrix_fftw2(sp_cmatrix * m){
  fftw_complex *out; 
  fftw_complex *in; 
  fftwnd_plan plan;
  int i;
  sp_cmatrix * res = sp_cmatrix_alloc(sp_cmatrix_rows(m),sp_cmatrix_cols(m));

  rephase(res);
  /* Rely on binary compatibility of the Complex type */
  in = m->data;
  out = res->data;
  plan = fftw2d_create_plan(sp_cmatrix_cols(m),sp_cmatrix_rows(m),FFTW_FORWARD,FFTW_ESTIMATE);

  fftwnd_one(plan,in,out);
  fftwnd_destroy_plan(plan);
  return res;
}

#endif


