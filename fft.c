#include <stdlib.h>
/* #include <sys/time.h>*/
#include <time.h>
#include <limits.h>
#include "spimage.h"



real * fft_shift(real * fftout, Detector *det){
  int x,y;
  real * img_d = malloc(sizeof(real)*(det->size[1]/2+1)*det->size[0]);
  for(x = 0;x<det->size[0];x++){
    for(y = 0;y<det->size[1]/2+1;y++){
      img_d[y+x*(det->size[1]/2+1)] = fftout[y+((det->size[0]/2+x)%det->size[0])*(det->size[1]/2+1)];
    }
  }
  return img_d;
}


int init_fft(int nthreads){
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

Image * image_rev_fftw3(Image * img){
  fftwr_complex *out; 
  fftwr_complex *in; 
  fftwr_plan plan;
  Image * res;
  int i,j;
  if(!img->phased){
    fprintf(stderr,"Error: Trying reverse fft an unphased image!\n");
    abort();
  }
  if(!img->shifted){
    fprintf(stderr,"Error: Trying to rev_fft an unshifted image!\n");
    abort();
  }

  res = create_empty_img(img);
  in = (fftwr_complex*) fftwr_malloc(sizeof(fftwr_complex) * (img->detector->size[1])*(img->detector->size[0]));  
  out = (fftwr_complex*) fftwr_malloc(sizeof(fftwr_complex) * (img->detector->size[1])*(img->detector->size[0]));
  plan = fftwr_plan_dft_2d(img->detector->size[0],img->detector->size[1],in,out, FFTW_BACKWARD,FFTW_ESTIMATE);

  j = 0;
  for(i = 0;i<TSIZE(img);i++){
      in[i][0] = img->r[i];
      in[i][1] = img->c[i];
  }
  fftwr_execute(plan);
  fftwr_destroy_plan(plan);
  fftwr_free(in);
  for(i = 0;i<TSIZE(img);i++){
    res->r[i] = out[i][0];
    res->c[i] = out[i][1];
    res->image[i] = norm(res,i);
  }
  fftwr_free(out);
  res->shifted = 0;
  return res;
}

#endif

#ifdef FFTW2

Image * image_rev_fftw2(Image * img){
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
  in = (fftw_complex*) malloc(sizeof(fftw_complex) * (img->detector->size[1])*(img->detector->size[0]));  
  out = (fftw_complex*) malloc(sizeof(fftw_complex) * (img->detector->size[1])*(img->detector->size[0]));
  plan = fftw2d_create_plan(img->detector->size[0],img->detector->size[1], FFTW_BACKWARD,FFTW_ESTIMATE);

  j = 0;
  for(i = 0;i<TSIZE(img);i++){
      in[i].re = img->r[i];
      in[i].im = img->c[i];
  }
  fftwnd_one(plan,in,out);
  fftwnd_destroy_plan(plan);
  free(in);
  for(i = 0;i<TSIZE(img);i++){
    res->r[i] = out[i].re;
    res->c[i] = out[i].im;
    res->image[i] = norm(res,i);
  }
  free(out);
  res->shifted = 0;
  return res;
}
#endif


#ifdef FFTW3

Image * image_guru_rev_fftw3(Image * img){
  fftwr_plan plan;
  int i;
  fftw_iodim dims[2];
  Image * res = create_empty_img(img);

  dims[0].n = img->detector->size[0];
  dims[0].is = img->detector->size[1];
  dims[0].os = img->detector->size[1];
  dims[1].n = img->detector->size[1];
  dims[1].is = 1;
  dims[1].os = 1;

  if(!img->phased){
    fprintf(stderr,"Error: Trying reverse fft an unphased image!\n");
    abort();
  }
  if(!img->shifted){
    fprintf(stderr,"Error: Trying to rev_fft an unshifted image!\n");
    abort();
  }

  if(!res->c || !res->r){
    fprintf(stderr,"out of memory?\n");
    abort();    
  }

  /* The exchange in complex and real part are on purpose. That's how you get the back fourier transform */
  plan = fftwr_plan_guru_split_dft(2,dims,0,dims,img->c,img->r,res->c,res->r,FFTW_ESTIMATE);
/*  if(!((problem_dft *)plan->prb)->io){
    fprintf(stderr,"FFTW bug? Plan problem with NULL io!\n");
    abort();
  }*/

  fftwr_execute(plan);
  fftwr_destroy_plan(plan);
  for(i = 0;i<TSIZE(img);i++){
    res->image[i] = norm(res,i);
  }
  res->shifted = 0;
  return res;
}

#endif

#ifdef FFTW3

Image * image_fftw3(Image * img){
  fftwr_complex *out; 
  fftwr_complex *in; 
  fftwr_plan plan;
  int i;
  Image * res = create_empty_img(img);

/*  if(!img->scaled){
    fprintf(stderr,"Error: Trying fft an unscaled image!\n");
    abort();
  }*/

  rephase(res);
  in = (fftwr_complex*) fftwr_malloc(sizeof(fftwr_complex) * (img->detector->size[1])*(img->detector->size[0]));
  out = (fftwr_complex*) fftwr_malloc(sizeof(fftwr_complex) * (img->detector->size[1])*(img->detector->size[0]));
  plan = fftwr_plan_dft_2d(img->detector->size[0],img->detector->size[1],in,out,FFTW_FORWARD,FFTW_ESTIMATE);

  if(img->phased){
    for(i = 0;i<TSIZE(img);i++){
      in[i][0] = img->r[i];
      in[i][1] = img->c[i];
    }
  }else{
    for(i = 0;i<TSIZE(img);i++){
      in[i][0] = img->image[i];
      in[i][1] = 0;
    }
  }
  fftwr_execute(plan);
  for(i = 0;i<(img->detector->size[1])*(img->detector->size[0]);i++){
   res->r[i] = out[i][0];
   res->c[i] =  out[i][1];
   res->image[i] =  norm(res,i);
  }
  fftwr_destroy_plan(plan);
  fftwr_free(in);
  fftwr_free(out);
  res->shifted = 1;
  res->detector->image_center[0] = (res->detector->size[0]-1)/2.0;
  res->detector->image_center[1] = (res->detector->size[1]-1)/2.0;
  return res;
}
#endif

#ifdef FFTW2

Image * image_fftw2(Image * img){
  fftw_complex *out; 
  fftw_complex *in; 
  fftwnd_plan plan;
  int i;
  Image * res = create_empty_img(img);

/*  if(!img->scaled){
    fprintf(stderr,"Error: Trying fft an unscaled image!\n");
    abort();
  }*/

  rephase(res);
  in = (fftw_complex*) malloc(sizeof(fftw_complex) * (img->detector->size[1])*(img->detector->size[0]));
  out = (fftw_complex*) malloc(sizeof(fftw_complex) * (img->detector->size[1])*(img->detector->size[0]));
  plan = fftw2d_create_plan(img->detector->size[0],img->detector->size[1],FFTW_FORWARD,FFTW_ESTIMATE);

  if(img->phased){
    for(i = 0;i<TSIZE(img);i++){
      in[i].re = img->r[i];
      in[i].im = img->c[i];
    }
  }else{
    for(i = 0;i<TSIZE(img);i++){
      in[i].re = img->image[i];
      in[i].im = 0;
    }
  }
  fftwnd_one(plan,in,out);
  for(i = 0;i<(img->detector->size[1])*(img->detector->size[0]);i++){
   res->r[i] = out[i].re;
   res->c[i] =  out[i].im;
   res->image[i] =  norm(res,i);
  }
  fftwnd_destroy_plan(plan);
  free(in);
  free(out);
  res->shifted = 1;
  res->detector->image_center[0] = (res->detector->size[0]-1)/2.0;
  res->detector->image_center[1] = (res->detector->size[1]-1)/2.0;
  return res;
}

#endif

#ifdef FFTW3

Image * image_guru_fftw3(Image * img){
  fftwr_plan plan;
  fftw_iodim dims[2];
  Image * res;
  real * zero;
  int i;

  dims[0].n = img->detector->size[0];
  dims[0].is = img->detector->size[1];
  dims[0].os = img->detector->size[1];
  dims[1].n = img->detector->size[1];
  dims[1].is = 1;
  dims[1].os = 1;
  res = create_empty_img(img);
  zero = NULL;

  if(!res->phased){
    res->c = malloc(sizeof(real)*TSIZE(res));
    res->r = malloc(sizeof(real)*TSIZE(res));
    res->phased = 1;
  }
  if(!res->c || !res->r){
    perror("out of memory?!?\n");
    abort();
  }

  if(img->phased){
    plan = fftwr_plan_guru_split_dft(2,dims,0,dims,img->r,img->c,res->r,res->c,FFTW_ESTIMATE);
  }else{
    zero = calloc(TSIZE(res),sizeof(real));
    plan = fftwr_plan_guru_split_dft(2,dims,0,dims,img->image,zero,res->r,res->c,FFTW_ESTIMATE);
  }
  fftwr_execute(plan);
  if(!img->phased){
    free(zero);
  }
  for(i = 0;i<(img->detector->size[1])*(img->detector->size[0]);i++){
    res->image[i] =  norm(res,i);
  }
  fftwr_destroy_plan(plan);
  res->shifted = 1;
  res->detector->image_center[0] = (res->detector->size[0]-1)/2.0;
  res->detector->image_center[1] = (res->detector->size[1]-1)/2.0;
  return res;
}



#endif

#ifdef FFTW3

Image * real_image_fft(Image * img){
  fftwr_complex *out; 
  fftwr_plan plan;
  int x,y;
  Image * res = create_empty_img(img);

  if(!img->scaled){
    fprintf(stderr,"Error: Trying fft an unscaled image!\n");
    abort();
  }

  rephase(res);
  out = (fftwr_complex*) fftwr_malloc(sizeof(fftwr_complex) * (img->detector->size[1]/2+1)*(img->detector->size[0]));
  plan = fftwr_plan_dft_r2c_2d(img->detector->size[0],img->detector->size[1],img->amplitudes,out,FFTW_ESTIMATE);
  fftwr_execute(plan);
  for(x = 0;x<img->detector->size[0];x++){
    for(y = 0;y<img->detector->size[1]/2+1;y++){
      res->r[x*img->detector->size[1]+y] = out[x*(img->detector->size[1]/2+1)+y][0];
      res->c[x*img->detector->size[1]+y] = out[x*(img->detector->size[1]/2+1)+y][1];
      res->amplitudes[x*img->detector->size[1]+y] = norm(res,x*img->detector->size[1]+y);

      res->r[(img->detector->size[0]-x-1)*img->detector->size[1]+(img->detector->size[1]-y-1)] = 
	out[x*(img->detector->size[1]/2+1)+y][0];
      res->c[(img->detector->size[0]-x-1)*img->detector->size[1]+(img->detector->size[1]-y-1)] = 
	out[x*(img->detector->size[1]/2+1)+y][1];
      res->amplitudes[(img->detector->size[0]-x-1)*img->detector->size[1]+(img->detector->size[1]-y-1)] =
	res->amplitudes[x*img->detector->size[1]+y];
    }
  }

  fftwr_destroy_plan(plan);
  fftwr_free(out);
  res->shifted = 1;
  res->detector->image_center[0] = (res->detector->size[0]-1)/2.0;
  res->detector->image_center[1] = (res->detector->size[1]-1)/2.0;
  return res;
}

#endif


