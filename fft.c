#include <stdlib.h>
/* #include <sys/time.h>*/
#include <time.h>
#include <limits.h>
#ifdef _USE_DMALLOC
#include <dmalloc.h>
#endif

#include "spimage.h"



real * sp_image_fft_shift(real * fftout, Image * a){
  int x,y,z;
  int nx = sp_c3matrix_x(a->image);
  int ny = sp_c3matrix_y(a->image);
  int nz = sp_c3matrix_z(a->image);
  
  real * img_d = sp_malloc(sizeof(real)*(nx/2+1)*ny*nz);
  for(x = 0;x<nx;x++){
    for(y = 0;y<ny;y++){
      for(z = 0;z<nz;z++){
	img_d[z*(nx/2+1)*(ny/2+1)+y*(nx/2+1)+x] = 
	  fftout[((nz/2+z)%nz)*(nx/2+1)*(ny/2+1)+((ny/2+y)%y)*(nx/2+1)+x];
	/* Not sure of this, didn't investigete to really understand
	   The old 2D version looked like:
	   img_d[y+x*(ny/2+1)] = fftout[y+((nx/2+x)%nx)*(ny/2+1)]; */
      }
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
  
  /* Kinda buggy FM: I'm gonna assume the compiler doesn't make strange things to our structure and everything is binary compatible */
  in = (fftwr_complex *)img->image->data;
  out = (fftwr_complex *)res->image->data;
  plan = fftwr_plan_dft_3d(sp_c3matrix_x(img->image),sp_c3matrix_y(img->image),sp_c3matrix_z(img->image),in,out, FFTW_BACKWARD,FFTW_ESTIMATE);

  fftwr_execute(plan);
  fftwr_destroy_plan(plan);
  res->shifted = 0;
  return res;
}

sp_c3matrix * sp_c3matrix_ifftw3(sp_c3matrix * m){
  fftwr_complex *out; 
  fftwr_complex *in; 
  fftwr_plan plan;
  sp_c3matrix * res;
  res = sp_c3matrix_alloc(sp_c3matrix_x(m),sp_c3matrix_y(m),sp_c3matrix_z(m));
  in = (fftwr_complex *)m->data;
  out = (fftwr_complex *)res->data;
  plan = fftwr_plan_dft_3d(sp_c3matrix_x(m),sp_c3matrix_y(m),sp_c3matrix_z(m),in,out, FFTW_BACKWARD,FFTW_ESTIMATE);

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
  plan = fftw3d_create_plan(sp_c3matrix_x(img->image),sp_c3matrix_y(img->image),sp_c3matrix_z(img->image), FFTW_BACKWARD,FFTW_ESTIMATE);

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
  in = (fftwr_complex *)img->image->data;
  out = (fftwr_complex *)res->image->data;
  plan = fftwr_plan_dft_3d(sp_c3matrix_x(img->image),sp_c3matrix_y(img->image),sp_c3matrix_z(img->image),in,out,FFTW_FORWARD,FFTW_ESTIMATE);

  fftwr_execute(plan);
  fftwr_destroy_plan(plan);
  res->shifted = 1;
  /*changed from
    res->detector->image_center[0] = (sp_c3matrix_x(res->image)-1)/2.0;
   */
  res->detector->image_center[0] = (sp_c3matrix_x(res->image))/2.0;
  res->detector->image_center[1] = (sp_c3matrix_y(res->image))/2.0;
  res->detector->image_center[2] = (sp_c3matrix_z(res->image))/2.0;
  return res;
}

Image * sp_image_1d_fftw3(Image * img, int axis) {
  fftwr_complex *out;
  fftwr_complex *in;
  fftwr_plan plan;
  Image *res = sp_image_duplicate(img,SP_COPY_DETECTOR);

  sp_image_rephase(res, SP_ZERO_PHASE);
  in = (fftwr_complex *) img->image->data;
  out = (fftwr_complex *) res->image->data;
  int inembed[1], onembed[1];
  inembed[0] = sp_image_size(img); onembed[0] = sp_image_size(img);
  int idist, odist, istride, ostride;
  int *n;
  int howmany;
  if (axis == 0) {
    idist = sp_image_x(img);
    istride = 1;
    n = malloc(sp_image_y(img)*sp_image_z(img)*sizeof(int));
    for (int i = 0; i < sp_image_y(img)*sp_image_z(img); i++) n[i] = sp_image_x(img);
    howmany = sp_image_y(img)*sp_image_z(img);
  } else if (axis == 1) {
    if (sp_image_z(img) != 1) printf("trying to 1d fft along the second y-axis of a 3d image is not yet implemented");
    idist = 1;
    istride = sp_image_x(img);
    n = malloc(sp_image_x(img)*sizeof(int));
    for (int i = 0; i < sp_image_x(img); i++) n[i] = sp_image_y(img);
    howmany = sp_image_x(img);
  } else if (axis == 2) {
    idist = 1;
    istride = sp_image_x(img)*sp_image_y(img);
    n = malloc(sp_image_x(img)*sp_image_y(img)*sizeof(int));
    for (int i = 0; i < sp_image_x(img)*sp_image_y(img); i++) n[i] = sp_image_z(img);
    howmany = sp_image_x(img)*sp_image_y(img);
  }
  odist = idist;
  ostride = istride;
  plan = fftwr_plan_many_dft(1, n, sp_image_y(img)*sp_image_z(img),
			     in, inembed, istride, idist, out, onembed, ostride, odist,
			     FFTW_FORWARD, FFTW_ESTIMATE);
  fftwr_execute(plan);
  fftwr_destroy_plan(plan);
  return res;
}

Image * sp_image_1d_ifftw3(Image * img, int axis) {
  fftwr_complex *out;
  fftwr_complex *in;
  fftwr_plan plan;
  Image *res = sp_image_duplicate(img,SP_COPY_DETECTOR);

  sp_image_rephase(res, SP_ZERO_PHASE);
  in = (fftwr_complex *) img->image->data;
  out = (fftwr_complex *) res->image->data;
  int inembed[1], onembed[1];
  inembed[0] = sp_image_size(img); onembed[0] = sp_image_size(img);
  int idist, odist, istride, ostride;
  int *n;
  int howmany;
  if (axis == 0) {
    idist = sp_image_x(img);
    istride = 1;
    n = malloc(sp_image_y(img)*sp_image_z(img)*sizeof(int));
    for (int i = 0; i < sp_image_y(img)*sp_image_z(img); i++) n[i] = sp_image_x(img);
    howmany = sp_image_y(img)*sp_image_z(img);
  } else if (axis == 1) {
    if (sp_image_z(img) != 1) printf("trying to 1d fft along the second y-axis of a 3d image is not yet implemented");
    idist = 1;
    istride = sp_image_x(img);
    n = malloc(sp_image_x(img)*sizeof(int));
    for (int i = 0; i < sp_image_x(img); i++) n[i] = sp_image_y(img);
    howmany = sp_image_x(img);
  } else if (axis == 2) {
    idist = 1;
    istride = sp_image_x(img)*sp_image_y(img);
    n = malloc(sp_image_x(img)*sp_image_y(img)*sizeof(int));
    for (int i = 0; i < sp_image_x(img)*sp_image_y(img); i++) n[i] = sp_image_z(img);
    howmany = sp_image_x(img)*sp_image_y(img);
  }
  odist = idist;
  ostride = istride;
  plan = fftwr_plan_many_dft(1, n, sp_image_y(img)*sp_image_z(img),
			     in, inembed, istride, idist, out, onembed, ostride, odist,
			     FFTW_BACKWARD, FFTW_ESTIMATE);
  fftwr_execute(plan);
  fftwr_destroy_plan(plan);
  return res;
}


sp_c3matrix * sp_c3matrix_fftw3(sp_c3matrix * m){
  fftwr_complex *out; 
  fftwr_complex *in; 
  fftwr_plan plan;
  sp_c3matrix * res = sp_c3matrix_alloc(sp_c3matrix_x(m),sp_c3matrix_y(m),
					sp_c3matrix_z(m));

  /* Rely on binary compatibility of the Complex type */
  in = (fftwr_complex *)m->data;
  out = (fftwr_complex *)res->data;
  plan = fftwr_plan_dft_3d(sp_c3matrix_x(m),sp_c3matrix_y(m),sp_c3matrix_z(m),in,out,FFTW_FORWARD,FFTW_ESTIMATE);

  fftwr_execute(plan);
  fftwr_destroy_plan(plan);
  return res;
}

#endif

#ifdef FFTW2

sp_c3matrix * sp_c3matrix_fftw2(sp_c3matrix * m){
  fftw_complex *out; 
  fftw_complex *in; 
  fftwnd_plan plan;
  int i;
  sp_c3matrix * res = sp_c3matrix_alloc(sp_c3matrix_x(m),sp_c3matrix_y(m),
					sp_c3matrix_z(m));

  rephase(res);
  /* Rely on binary compatibility of the Complex type */
  in = m->data;
  out = res->data;
  plan = fftw3d_create_plan(sp_c3matrix_x(m),sp_c3matrix_y(m),sp_c3matrix_z(m),
			    FFTW_FORWARD,FFTW_ESTIMATE);
  fftwnd_one(plan,in,out);
  fftwnd_destroy_plan(plan);
  return res;
}

#endif


