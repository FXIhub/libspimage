#ifdef __MINGW32__
/* resolve conflict with ssize_t*/
#define _NO_OLDNAMES
typedef long off_t;
#endif

#ifndef PNG_DEBUG
#  define PNG_DEBUG 3
#endif

#include <stdlib.h>
#include <math.h>
#include <hdf5.h>
#include <tiffio.h>
#include <png.h>
#include <float.h>
#include "spimage.h"


static Image * zero_pad_shifted_image(Image * a, int newx, int newy, int pad_mask);
static Image * zero_pad_unshifted_image(Image * a, int newx, int newy, int pad_mask);

real p_drand48(){
	real ret = ((double)rand())/(double)RAND_MAX;
  return (ret);
}

/* Returns the patterson from diffraction image a */
Image * patterson_function(Image * a){
  Image * b = imgcpy(a);
  Image * c;
  Image * d;
  int i;

  if(b->scaled){
    for(i = 0;i<TSIZE(a);i++){
      b->image[i] *= b->image[i];
    }
  }
  c = image_fft(b);

  freeimg(b);
  d = shift_quadrants(c);
  freeimg(c);
  return d; 
}

/* 
   Returns in the image correlated with it's
   centrosymetric image

   Clevearly implement using a single convolution
   (after a bug hunt)
*/
Image * centrosym_convolve(Image * in){
/*  Image *csym = reflect_xy(in,OUT_OF_PLACE);*/
  Image * res;
  res = convolute_img(in,in,NULL);/*csym);*/
/*  freeimg(csym);*/
  return res;
}

/* Make all phases 0, also known
   as making the image real */
Image * make_real(Image * in, int in_place){
  Image * out;
  int i;
  if(!in->phased){
    fprintf(stderr,"Trying to run make_real on unphased image!");
    abort();
  }
  if(in_place){
    out = in;
  }else{
    out = imgcpy(in);
  }
  for(i = 0;i<TSIZE(in);i++){
    out->r[i] = out->image[i];
    out->c[i] = 0;
  }
  return out;
  
}

/* reflect an image through the center, also 
 known as rotating by 180 degrees */
Image * reflect_xy(Image * in, int in_place){
  return reflect_x(reflect_y(in,in_place),1);
}

/* reflect an image on the x axis, also 
 known as vertical flip */
Image * reflect_x(Image * in, int in_place){
  real tmp;
  Image * out;
  int x,y,y2;
  if(in_place){
    out = in;
  }else{
    out = imgcpy(in);
  }
  for(x = 0;x<in->detector->size[0];x++){
    for(y = 0;y<in->detector->size[1]/2.0;y++){
      y2 = in->detector->size[1]-y-1;
      tmp = in->image[x*in->detector->size[1]+y];
      out->image[x*in->detector->size[1]+y] =in->image[x*in->detector->size[1]+y2];
      out->image[x*in->detector->size[1]+y2] = tmp;
      if(in->phased){
	tmp = in->r[x*in->detector->size[1]+y];
	out->r[x*in->detector->size[1]+y] =in->r[x*in->detector->size[1]+y2];
	out->r[x*in->detector->size[1]+y2] = tmp;
	tmp = in->c[x*in->detector->size[1]+y];
	out->c[x*in->detector->size[1]+y] =in->c[x*in->detector->size[1]+y2];
	out->c[x*in->detector->size[1]+y2] = tmp;	  
      }
    }
  }
  return out;
}

/* reflect an image on the y axis, also 
 known as horizontal flip */
Image * reflect_y(Image * in, int in_place){
  real tmp;
  Image * out;
  int x,y,x2;
  if(in_place){
    out = in;
  }else{
    out = imgcpy(in);
  }
  for(x = 0;x<in->detector->size[0]/2.0;x++){
    x2 = in->detector->size[0]-x-1;
    for(y = 0;y<in->detector->size[1];y++){
      tmp = in->image[x*in->detector->size[1]+y];
      out->image[x*in->detector->size[1]+y] =in->image[x2*in->detector->size[1]+y];
      out->image[x2*in->detector->size[1]+y] = tmp;
      if(in->phased){
	tmp = in->r[x*in->detector->size[1]+y];
	out->r[x*in->detector->size[1]+y] =in->r[x2*in->detector->size[1]+y];
	out->r[x2*in->detector->size[1]+y] = tmp;
	tmp = in->c[x*in->detector->size[1]+y];
	out->c[x*in->detector->size[1]+y] =in->c[x2*in->detector->size[1]+y];
	out->c[x2*in->detector->size[1]+y] = tmp;	  
      }
    }
  }
  return out;
}

/* Image should be shifted and squared */
Image * average_centrosymetry(Image * in){
  int x,y;
  real noise = 0;
  int ind1,ind2;
  int k = 0;
  Image * out = imgcpy(in);
  
  if(!in->shifted){
    fprintf(stderr,"Error: Using average_centrosymetry on unshifted image!\n");
    exit(1);
  }
  if(in->detector->size[0] != in->detector->size[1]){
    fprintf(stderr,"Error: Using average_centrosymetry with a non square image!\n");
    exit(1);
  }

  for(x = 0;x<in->detector->size[0];x++){
    for(y = 0;y<in->detector->size[1];y++){
      ind1 = (in->detector->size[0]-1-x)*in->detector->size[1]+(in->detector->size[1]-1-y);
      ind2 = x*in->detector->size[1]+y;
      if(dist_to_corner(x*in->detector->size[1]+y,in) < 1500){
	if((in->image[ind1] + in->image[ind2]) && in->mask[ind1] && in->mask[ind2]){
	  noise += fabs(in->image[ind1] - in->image[ind2])/
 	  ((in->image[ind1] + in->image[ind2])/2);
	  k++;
	  out->image[ind2] = (in->image[ind1] + in->image[ind2])/2;
	}else{
	  out->image[ind2] = 0;
	}
      }
    }
  }
  fprintf(stderr,"Noise = %f\n",noise/k);
  fprintf(stderr,"Noise2 = %f\n",noise);
  return out;
}


Image * make_shifted_image_square(Image * in){
  Image * out;
  int x,y;
  int xout,yout;
  out = imgcpy(in);

  if(!in->shifted){
    fprintf(stderr,"Trying to call make_shifted_image_square with an unshifted image|\n");
    exit(1);
  }

  /* make it a square by removing a line or a column from the center */
  if(in->detector->size[0] > in->detector->size[1]){
    /* remove a column */
    out->detector->size[0] = in->detector->size[1];
    free(out->image);
    out->image = malloc(sizeof(real)*TSIZE(out));
    if(out->scaled){
      out->amplitudes = out->image;
    }else{
      out->intensities = out->image;
    }
    free(out->mask);
    out->mask = malloc(sizeof(real)*TSIZE(out));
    yout = 0;
    xout = 0;
    for(x = 0;x<in->detector->size[0];x++){
      for(y = 0;y<in->detector->size[1];y++){
	if(fabs(x - (in->detector->size[0]-1)/2.0) < (in->detector->size[0]-out->detector->size[1])/2.0){
	  continue;
	}
	out->image[xout*out->detector->size[1]+yout] = in->image[x*in->detector->size[1]+y];
	out->mask[xout*out->detector->size[1]+yout] = in->mask[x*in->detector->size[1]+y];
	yout = (yout+1)%out->detector->size[1];
	if(yout == 0){
	  xout++;
	}
      }
    }
  }else if(in->detector->size[0] < in->detector->size[1]){
    /* remove a line */
    out->detector->size[1] = in->detector->size[0];
    free(out->image);
    out->image = malloc(sizeof(real)*TSIZE(out));      
    if(out->scaled){
      out->amplitudes = out->image;
    }else{
      out->intensities = out->image;
    }
    free(out->mask);
    out->mask = malloc(sizeof(real)*TSIZE(out));
    yout = 0;
    xout = 0;

    for(x = 0;x<in->detector->size[0];x++){
      for(y = 0;y<in->detector->size[1];y++){
	if(fabs(x - (in->detector->size[1]-1)/2.0) < (in->detector->size[1]-out->detector->size[0])/2.0){
	  continue;
	}
	out->image[xout*out->detector->size[1]+yout] = in->image[x*in->detector->size[1]+y];
	out->mask[xout*out->detector->size[1]+yout] = in->mask[x*in->detector->size[1]+y];
	yout = (yout+1)%out->detector->size[1];
	if(yout == 0){
	  xout++;
	}

      }
    }
  }
  return out;
}



/* Make an unshifted image square by padding with zeroes on the smaller dimension */
Image * make_unshifted_image_square(Image * in){
  int size = MAX(in->detector->size[0],in->detector->size[1]);
  return zero_pad_image(in,size,size,0);
}

Image * old_make_unshifted_image_square(Image * in){
  Image * out;
  int x,y;
  int xout,yout;
  int llimit,ulimit;
  xout = 0;
  yout = 0;
  if(in->shifted){
    fprintf(stderr,"Trying to call make_unshifted_image_square with a shifted image|\n");
    exit(1);
  }

  out = imgcpy(in);
  /* make it a square by removing a line or a column from the sides */
  if(in->detector->size[0] > in->detector->size[1]){
    /* reposition the center */
    out->detector->image_center[0] *= in->detector->size[1];
    out->detector->image_center[0] /= out->detector->size[0];
    /* remove columns */
    out->detector->size[0] = in->detector->size[1];
    free(out->image);
    out->image = malloc(sizeof(real)*TSIZE(out));
    if(out->scaled){
      out->amplitudes = out->image;
    }else{
      out->intensities = out->image;
    }
    
    free(out->mask);
    out->mask = malloc(sizeof(real)*TSIZE(out));
    yout = 0;
    xout = 0;
    llimit = (in->detector->size[0]-out->detector->size[1])/2.0;
    ulimit = in->detector->size[0] - ((in->detector->size[0]-out->detector->size[1]) - llimit);
    for(x = 0;x<in->detector->size[0];x++){
      for(y = 0;y<in->detector->size[1];y++){
	if(x < llimit || x >= ulimit){
	  continue;
	}
	out->image[xout*out->detector->size[1]+yout] = in->image[x*in->detector->size[1]+y];
	out->mask[xout*out->detector->size[1]+yout] = in->mask[x*in->detector->size[1]+y];
	yout = (yout+1)%out->detector->size[1];
	if(yout == 0){
	  xout++;
	}
      }
    }
  }else if(in->detector->size[0] < in->detector->size[1]){
    /* reposition the center */
    out->detector->image_center[1] *= in->detector->size[0];
    out->detector->image_center[1] /= out->detector->size[1];
    /* remove lines */
    out->detector->size[1] = in->detector->size[0];
    free(out->image);
      out->image = malloc(sizeof(real)*TSIZE(out));
    if(out->scaled){
      out->amplitudes = out->image;
    }else{
      out->intensities = out->image;
    }
    free(out->mask);
    out->mask = malloc(sizeof(real)*TSIZE(out));
    yout = 0;
    xout = 0;

    llimit = (in->detector->size[1]-out->detector->size[0])/2.0;
    ulimit = in->detector->size[1] - ((in->detector->size[1]-out->detector->size[0]) - llimit);

    for(x = 0;x<in->detector->size[0];x++){
      for(y = 0;y<in->detector->size[1];y++){
	if(y < llimit || y >= ulimit){
	  continue;
	}
	out->image[xout*out->detector->size[1]+yout] = in->image[x*in->detector->size[1]+y];
	out->mask[xout*out->detector->size[1]+yout] = in->mask[x*in->detector->size[1]+y];
	yout = (yout+1)%out->detector->size[1];
	if(yout == 0){
	  xout++;
	}

      }
    }
  }
  return out;
}


/*
  For unshifted images it shifted the quadrants around image_center and 
  extra zero pad is used on the smaller quadrants to make them the same size
  as the biggest quadrant.
  For shifted images it shifted the quadrants around (size[]-1)/2 

*/
Image * shift_quadrants(Image * img){
  Image * out;
  int i;
  int index1,index2;
  int x,y;
  int newx,newy;
  int max_x,max_y;

  /* fft shift the image */
  out = imgcpy(img);
  if(!img->shifted){
    max_x = MAX(img->detector->image_center[0],img->detector->size[0]-img->detector->image_center[0]);
    max_y = MAX(img->detector->image_center[1],img->detector->size[1]-img->detector->image_center[1]);
    resize_empty_image(out,2*max_x,2*max_y);
  }

		   
  for(i = 0;i<TSIZE(out);i++){
    out->image[i] = 0;
    out->mask[i] = 0;
  }
  if(img->shifted){
    out->detector->image_center[0] = (img->detector->size[0])/2.0;
    out->detector->image_center[1] = (img->detector->size[1])/2.0;
    img->detector->image_center[0] = (img->detector->size[0])/2.0;
    img->detector->image_center[1] = (img->detector->size[0])/2.0;
  }else{
    out->detector->image_center[0] = 0;
    out->detector->image_center[1] = 0;
  }
  out->shifted = !out->shifted;
  /* shift quadrants */
  for(x = 0;x<img->detector->size[0];x++){
    for(y = 0;y<img->detector->size[1];y++){
      index1 = x*img->detector->size[1]+y;
      index2 = 0;
      if(y < img->detector->image_center[1]){
	if(x < img->detector->image_center[0]){
	  newx = out->detector->size[0]-(img->detector->image_center[0]-x);
	  newy = out->detector->size[1]-(img->detector->image_center[1]-y);
	  if(newx < img->detector->size[0]/2.0 ||
	     newy < img->detector->size[1]/2.0){
	    index2 = -1;
	  }
	}else{
	  newx = x-img->detector->image_center[0];
	  newy = out->detector->size[1]-(img->detector->image_center[1]-y);
	  if(newx >= img->detector->size[0]/2.0 ||
	     newy < img->detector->size[1]/2.0){
	    index2 = -1;
	  }
	}
      }else{	
	if(x < img->detector->image_center[0]){
	  newx = out->detector->size[0]-(img->detector->image_center[0]-x);
	  newy = y-img->detector->image_center[1];
	  if(newx < img->detector->size[0]/2.0 ||
	     newy >= img->detector->size[1]/2.0){
	    index2 = -1;
	  }
	}else{
	  newx = x-img->detector->image_center[0];
	  newy = y-img->detector->image_center[1];
	  if(newx >= img->detector->size[0]/2.0 ||
	     newy >= img->detector->size[1]/2.0){
	    index2 = -1;
	  }

	}
      }
      if(index2 != -1){
	index2 = (newx)*out->detector->size[1]+newy;
      }
/*      if(newx >= img->detector->size[0]){
	index2 = -1;
      }
      if(newy >= img->detector->size[1]){
	index2 = -1;
      }
*/
      if(index2 != -1){
	if(img->phased){
	  out->r[index2] = img->r[index1];
	  out->c[index2] = img->c[index1];
	}
	out->image[index2] = img->image[index1];
	out->mask[index2] = img->mask[index1];
      }
    }
  }

  
  return out;
}


/* resolution given in pixels from center. 
   Image should be shifted. Square window. */
Image * limit_resolution(Image * img, int resolution){
  Image * res = imgcpy(img);
  int x,y,nx,ny;
  int dx,dy;
  if(img->shifted == 0){
    fprintf(stderr,"Error: Trying to limit resolution on an unshifted image\n");
  }
  if(resolution*2 > res->detector->size[0] || resolution*2 > res->detector->size[1]){
    return imgcpy(img);
  }  
  res->detector->size[0] = resolution*2;
  res->detector->size[1] = resolution*2;
  free(res->image);
  free(res->mask);
  res->image = malloc(sizeof(real)*TSIZE(res));
  if(res->scaled){
    res->amplitudes = res->image;
  }else{
    res->intensities = res->image;
  }
  if(res->phased){
    free(res->c);
    free(res->r);
    res->c = malloc(sizeof(real)*TSIZE(res));
    res->r = malloc(sizeof(real)*TSIZE(res));
  }
  res->mask = malloc(sizeof(real)*TSIZE(res));
  nx = 0;
  ny = 0;
  for(x = 0;x<img->detector->size[0];x++){
    dx = x;
    if(img->detector->size[0]-x-1 < dx){
      dx = img->detector->size[0]-x-1;
    }
    for(y = 0;y<img->detector->size[1];y++){
      dy = y;
      if(img->detector->size[1]-y-1 < dy){
	dy = img->detector->size[1]-y-1;
      }
      if(dx < resolution && dy < resolution){
	res->image[nx*res->detector->size[1]+ny] = img->image[x*img->detector->size[1]+y];
	if(res->phased){
	  res->c[nx*res->detector->size[1]+ny] = img->c[x*img->detector->size[1]+y];
	  res->r[nx*res->detector->size[1]+ny] = img->r[x*img->detector->size[1]+y];
	}
	res->mask[nx*res->detector->size[1]+ny] = img->mask[x*img->detector->size[1]+y];
	ny++;
	if(ny == res->detector->size[1]){
	  nx++;
	  ny = 0;
	}
      }      
    }    
  }
  return res;
}



/* resolution given in pixels from center. 
   Image should be shifted. Square window. */
void smooth_edges(Image * img, int resolution){
  real cutoff = 0.80;
  real tmp;
  int i;
  if(img->shifted == 0){
    fprintf(stderr,"Error: Trying to smooth edges on an unshifted image\n");
  }
  for(i =0 ;i<TSIZE(img);i++){
    if(dist_to_center(i,img) > resolution*cutoff){
      if(dist_to_center(i,img) <= resolution){
	/* tmp goes from 0 to 1 with higher resolution */
	tmp = (1/(1.0-cutoff))*(dist_to_center(i,img)/resolution-cutoff);
	tmp *= 1;
	tmp *= tmp;
	img->image[i] *= exp(-(tmp));
	if(img->phased){
	  img->r[i] *= exp(-(tmp));
	  img->c[i] *= exp(-(tmp));
	}
      }else{
	img->image[i] = 0;
	if(img->phased){
	  img->r[i] = 0;
	  img->c[i] = 0;
	}

      }
    }
  }
}


real dist_to_axis(int i, Image * in){
  int x = i/in->detector->size[1];
  int y = i%in->detector->size[1];
  real dx,dy;
  if(in->shifted){
    dx = MIN(x,in->detector->size[0]-x);
    dy = MIN(y,in->detector->size[1]-y);
  }else{
    dx = fabs(x-in->detector->image_center[0]);
    dy = fabs(y-in->detector->image_center[1]);
  }
  return MIN(dx,dy);
}

real centro_sym_value(int index,Image * img){
  int x,y;
  real nx,ny;
  x = index/img->detector->size[1];
  y = index%img->detector->size[1];
  nx = 2*img->detector->image_center[0]-x;
  ny = 2*img->detector->image_center[1]-y;
  if(nx < 1 || nx >= img->detector->size[0]-2 ||
     ny < 1 || ny >= img->detector->size[1]-2){
    return -1;
  }
  return lin_image_interpol(img, nx, ny);
}

int centro_sym_index(int index,Image * img){
  int x,y;
  int nx,ny;
  x = index/img->detector->size[1];
  y = index%img->detector->size[1];
  nx = 2*img->detector->image_center[0]-x;
  ny = 2*img->detector->image_center[1]-y;
  if(nx < 0 || nx >= img->detector->size[0] ||
     ny < 0 || ny >= img->detector->size[1]){
    return index;
  }
  return nx*img->detector->size[1]+ny;
}

Image * centro_sym_correlation(Image  * img){
  int x,y;
  Image * res = imgcpy(img);
  int index = 0;
  real csvalue;
  for(x = 0;x<img->detector->size[0];x++){
    for(y = 0;y<img->detector->size[1];y++){
      csvalue = centro_sym_value(index,img);
      if(!img->mask[index] || csvalue == -1 || fabs(img->image[index])+fabs(csvalue) < 1){
	res->image[index] = 1.0;
      }else{
	res->image[index] = 1.0 - fabs(img->image[index]-csvalue)/(fabs(img->image[index])+fabs(csvalue));      
      }
      if(res->image[index] < 0 || res->image[index] > 1){
	/* Houston we have a problem */
	exit(1);
      }
      index++;
    }
  }
  return res;
}


real lin_image_interpol(Image * img, real x, real y){
  int ix = x;
  int iy = y;
  real t = (x-ix);
  real u = (y-iy);
  return (1.0-t) * (1.0-u) * img->image[ix*img->detector->size[1]+iy] +
    t * (1.0-u) * img->image[(ix+1)*img->detector->size[1]+iy] +
    t * u * img->image[(ix+1)*img->detector->size[1]+iy+1] +
    (1.0-t) * u * img->image[ix*img->detector->size[1]+iy+1];
}


real dist_to_center(int i, Image * in){
  int x = i/in->detector->size[1];
  int y = i%in->detector->size[1];
  real dx,dy;
  if(in->shifted){
    dx = MIN(x,in->detector->size[0]-x);
    dy = MIN(y,in->detector->size[1]-y);
  }else{
    dx = x-in->detector->image_center[0];
    dy = y-in->detector->image_center[1];
  }
  return sqrt(dx*dx+dy*dy);
}

real square_dist_to_center(int i, Image * in){
  int x = i/in->detector->size[1];
  int y = i%in->detector->size[1];
  real dx,dy;
  if(in->shifted){
    dx = MIN(x,in->detector->size[0]-x);
    dy = MIN(y,in->detector->size[1]-y);
  }else{
    dx = fabs(x-in->detector->image_center[0]);
    dy = fabs(y-in->detector->image_center[1]);
  }
  return MAX(dx,dy);

}

real dist_to_corner(int i, Image * in){
  int x = i/in->detector->size[1];
  int y = i%in->detector->size[1];
  real dx,dy;
  if(in->detector->size[0]-1-x < x){
    dx = in->detector->size[0]-1-x;
  }else{
    dx = x;
  }
  if(in->detector->size[1]-1-y < y){
    dy = in->detector->size[0]-1-y;
  }else{
    dy = y;
  }
  return sqrt(dx*dx+dy*dy);
}


void dephase(Image *  img){
  if(img->phased){
    free(img->r);
    free(img->c);
    img->phased = 0;    
  }
}

void rephase(Image *  img){ 
/*  if(!img->scaled){
    fprintf(stderr,"Error: Trying to phase unscaled image!\n");
    exit(1);
  }*/
  if(!img->phased){
    img->r = calloc(TSIZE(img),sizeof(real));
    img->c = calloc(TSIZE(img),sizeof(real));
    memcpy(img->r,img->image,TSIZE(img)*sizeof(real));
    img->phased = 1;    
  }
}

void random_rephase(Image *  img){  
  int i;
  real phase;
  if(!img->scaled){
    fprintf(stderr,"Error: Trying to phase unscaled image!\n");
    exit(1);
  }
  if(!img->phased){
    img->r = calloc(TSIZE(img),sizeof(real));
    img->c = calloc(TSIZE(img),sizeof(real));
  }
  for(i = 0;i<TSIZE(img);i++){
    phase = p_drand48()*M_PI*2;
    img->r[i] = cos(phase)*img->image[i];
    img->c[i] = sin(phase)*img->image[i];
  }
  img->phased = 1;    
}


/* returns the phase of the image from [0, 2*M_PI[ or -1 if there's an error */
Image * sp_image_get_phases(Image * img){
  Image * res;
  int i;
  if(!img->phased){
    return NULL;
  }
  res = create_empty_img(img);
  dephase(res);
  for(i = 0;i<TSIZE(res);i++){
    res->image[i] = get_phase_angle(img,i);
  }
  return res;
}


/* returns the phase of index i from [0, 2*M_PI[ or -1 if there's an error */
Image * get_phase_image(Image * img){
  Image * res;
  int i;
  if(!img->phased){
    return NULL;
  }
  res = imgcpy(img);
  for(i = 0;i<TSIZE(res);i++){
    res->r[i] /= res->image[i];
    res->c[i] /= res->image[i];
    res->image[i] = 1;
  }
  return res;
}


/* returns the phase of index i from ]-M_PI, M_PI] or -1 if there's an error */
real get_phase_angle(Image * img, int i){
  if(!img->phased){
    return -1;
  }
  if(img->c[i]){
    return atan(img->r[i]/img->c[i]);
  }else{
    if(img->r[i] < 0){
      return M_PI;
    }else{
      return 0;
    }
  }

}


real real_factor(Image * img, int i){
  if(img->image[i]){
    return img->r[i]/img->image[i];
  }else{
    return 1;
  }
}

real complex_factor(Image * img, int i){
  if(img->image[i]){
    return img->c[i]/img->image[i];
  }else{
    return 0;
  }
}



inline real norm(Image * img, int i){
  return sqrt(img->r[i]*img->r[i]+img->c[i]*img->c[i]);
}


void add_image(Image * a, Image * b){
  int i;
  if(a->phased && b->phased){
    for(i = 0;i<TSIZE(a);i++){
      a->r[i] = a->r[i]+b->r[i];
      a->c[i] = a->c[i]+b->c[i];
      a->image[i] = norm(a,i);
    }
  }else{
    for(i = 0;i<TSIZE(a);i++){
      a->image[i] = a->image[i]+b->image[i];
    }
  }
}

void sp_image_sub(Image * a, Image * b){
  int i;
  if(a->phased && b->phased){
    for(i = 0;i<TSIZE(a);i++){
      a->r[i] = a->r[i]-b->r[i];
      a->c[i] = a->c[i]-b->c[i];
      a->image[i] = norm(a,i);
    }
  }else{
    for(i = 0;i<TSIZE(a);i++){
      a->image[i] = a->image[i]-b->image[i];
    }
  }
}

/* mean m, standard deviation s */
/* gaussian distribution generator */
float box_muller(float m, float s){        
  float x1, x2, w, y1;
  static float y2;
  static int use_last = 0;
  
  if (use_last){        /* use value from previous call */
    y1 = y2;
    use_last = 0;
  }  else {
    do {
      x1 = 2.0 * p_drand48() - 1.0;
      x2 = 2.0 * p_drand48() - 1.0;
      w = x1 * x1 + x2 * x2;
    } while ( w >= 1.0 );
    
    w = sqrt( (-2.0 * log( w ) ) / w );
    y1 = x1 * w;
    y2 = x2 * w;
    use_last = 1;
  }  
  return( m + y1 * s );
}


void add_gaussian_noise(Image * in, real level){
  int i;
  if(level <= 0){
    return;
  }
  if(!in->scaled){    
    fprintf(stderr,"Error: Trying add gaussian noise to unscaled image!\n");
    exit(1);
  }

  for(i = 0;i<TSIZE(in);i++){
    if(in->phased){
      in->r[i] *= (1.0+box_muller(0,level));
      in->c[i] *= (1.0+box_muller(0,level));
      in->image[i] = norm(in,i);      
    }else{
      in->image[i] *= (1.0+box_muller(0,level));
    }
  }
}

void remove_lowres(Image * in, real radius){
  int x,y;
  real dist,dx,dy;
  if(radius <= 0){
    return;
  }
  for(x = 0;x<in->detector->size[0];x++){
    for(y = 0;y<in->detector->size[1];y++){
      if(x > in->detector->size[0]/2.0){
	dx = in->detector->size[0]-x;
      }else{
	dx = x;
      }
      if(y > in->detector->size[1]/2.0){
	dy = in->detector->size[1]-y;
      }else{
	dy = y;
      }
      dist = sqrt(dx*dx+dy*dy);      
      if(dist <= radius){
	in->mask[x*in->detector->size[1]+y] = 0;
	if(in->phased){
	  in->r[x*in->detector->size[1]+y] = 0;
	  in->c[x*in->detector->size[1]+y] = 0;
	}
	in->image[x*in->detector->size[1]+y] = 0;
      }
    }
  }
}

void freeimg(Image * in){
  free(in->detector);
  free(in->image);
  if(in->phased){
    free(in->c);
    free(in->r);
  }
  free(in->mask);
  free(in);
}

Image * imgcpy(Image * in){
  Image  *res = malloc(sizeof(Image));
  if(!res){
    perror("Out of memory!\n");
    abort();
  }
  memcpy(res,in,sizeof(Image));
  res->detector = malloc(sizeof(Detector));
  if(!res->detector){
    perror("Out of memory!\n");
    abort();
  }

  memcpy(res->detector,in->detector,sizeof(Detector));
  res->image = malloc(sizeof(real)*TSIZE(res));
  if(!res->image){
    perror("Out of memory!\n");
    abort();
  }
  memcpy(res->image,in->image,sizeof(real)*TSIZE(res)); 
  if(in->scaled){
    res->amplitudes = res->image;
  }else{
    res->intensities = res->image;
  }
  res->mask = malloc(sizeof(real)*TSIZE(res));
  memcpy(res->mask,in->mask,sizeof(real)*TSIZE(res));
  if(in->phased){
    res->c = malloc(sizeof(real)*TSIZE(res));
    if(!res->c){
      perror("Out of memory!\n");
      abort();
    }
    memcpy(res->c,in->c,sizeof(real)*TSIZE(res));
    res->r = malloc(sizeof(real)*TSIZE(res));
    if(!res->r){
      perror("Out of memory!\n");
      abort();
    }
    memcpy(res->r,in->r,sizeof(real)*TSIZE(res));
  }
  return res;
}

/* creates an image similar to a but empty */
Image * create_empty_img(Image * in){
  Image  *res = malloc(sizeof(Image));
  if(!res){
    perror("Out of memory!\n");
    abort();
  }
  memcpy(res,in,sizeof(Image));
  res->detector = malloc(sizeof(Detector));
  if(!res->detector){
    perror("Out of memory!\n");
    abort();
  }
  memcpy(res->detector,in->detector,sizeof(Detector));
  res->mask = calloc(TSIZE(res),sizeof(real));
  res->image = calloc(TSIZE(res),sizeof(real));
  if(!res->image){
    perror("Out of memory!\n");
    abort();
  }  
  if(in->scaled){
    res->amplitudes = res->image;
  }else{
    res->intensities = res->image;
  }
  if(in->phased){
    res->c = calloc(TSIZE(res),sizeof(real));
    if(!res->c){
      perror("Out of memory!\n");
      abort();
    }
    res->r = calloc(TSIZE(res),sizeof(real));
    if(!res->r){
      perror("Out of memory!\n");
      abort();
    }
  }
  return res;
}

Image * create_new_img(int x, int y){
  Image  *res = malloc(sizeof(Image));
  if(!res){
    perror("Out of memory!\n");
    abort();
  }
  res->detector = malloc(sizeof(Detector));
  res->detector->size[0] = x;
  res->detector->size[1] = y;
  res->detector->image_center[0] = x/2;
  res->detector->image_center[1] = y/2;
  if(!res->detector){
    perror("Out of memory!\n");
    abort();
  }
  res->mask = calloc(TSIZE(res),sizeof(real));
  res->scaled = 0;
  res->intensities = calloc(TSIZE(res),sizeof(real));
  if(!res->intensities){
    perror("Out of memory!\n");
    abort();
  }
  res->image = res->intensities;
  res->phased = 0;
  return res;
}



void write_img(Image * img,char * filename, int output_precision){
  hid_t dataspace_id;
  hid_t dataset_id;
  hid_t file_id;
  int status;
  hsize_t  dims[2];
  real values[2];
  real * phases;
  int i;
  hid_t out_type_id = 0;
  hid_t mem_type_id = 0;
  hid_t plist;
  hsize_t chunk_size[2] = {img->detector->size[0],img->detector->size[1]};
  if(output_precision == sizeof(double)){
    out_type_id = H5T_NATIVE_DOUBLE;
  }else if(output_precision == sizeof(float)){
    out_type_id = H5T_NATIVE_FLOAT;
  }else{
    abort();
  }
  if(sizeof(real) == sizeof(float)){
    mem_type_id = H5T_NATIVE_FLOAT;
  }else if(sizeof(real) == sizeof(double)){
    mem_type_id = H5T_NATIVE_DOUBLE;
  }else{
    abort();
  }

  dims[0] = img->detector->size[0];
  dims[1] = img->detector->size[1];
  file_id = H5Fcreate(filename,  H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  dataspace_id = H5Screate_simple( 2, dims, NULL );

  plist = H5Pcreate (H5P_DATASET_CREATE);
  H5Pset_chunk(plist,2,chunk_size);
  H5Pset_deflate(plist,6);

  if(img->scaled){
    dataset_id = H5Dcreate(file_id, "/amplitudes", out_type_id,
			   dataspace_id, plist);
    status = H5Dwrite(dataset_id,mem_type_id , H5S_ALL, H5S_ALL,
		      H5P_DEFAULT, img->amplitudes);
    status = H5Dclose(dataset_id);
  }else{
    dataset_id = H5Dcreate(file_id, "/intensities", out_type_id,
			   dataspace_id, plist);
    status = H5Dwrite(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		      H5P_DEFAULT, img->intensities);
    status = H5Dclose(dataset_id);
  }


  dataset_id = H5Dcreate(file_id, "/mask", out_type_id,
			 dataspace_id, plist);
  status = H5Dwrite(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, img->mask);
  status = H5Dclose(dataset_id);
  if(img->phased){
    phases = malloc(sizeof(real)*TSIZE(img));
    for(i = 0;i<TSIZE(img);i++){
      if(img->r[i] < 0){
	phases[i] = atan(img->c[i]/img->r[i])+M_PI;
      }else if(img->r[i] > 0){
	phases[i] = atan(img->c[i]/img->r[i]);	
      }else{
	if(img->c[i] > 0){
	  phases[i] = M_PI/2.0;
	}else if(img->c[i] < 0){
	  phases[i] = -M_PI/2.0;
	}else{
	  /* 0 amplitude, arbitrary 0 phase */
	  phases[i] = 0;
	}
      }
    }

    dataset_id = H5Dcreate(file_id, "/phases", out_type_id,
			   dataspace_id, plist);
    status = H5Dwrite(dataset_id,mem_type_id, H5S_ALL, H5S_ALL,
		      H5P_DEFAULT, phases);
    free(phases);
    status = H5Dclose(dataset_id);	
    dataset_id = H5Dcreate(file_id, "/real", out_type_id,
			   dataspace_id, plist);
    status = H5Dwrite(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		      H5P_DEFAULT, img->r);
    status = H5Dclose(dataset_id);
    dataset_id = H5Dcreate(file_id, "/complex",out_type_id,
			   dataspace_id, plist);
    status = H5Dwrite(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, img->c);
    status = H5Dclose(dataset_id);
  }
  status = H5Sclose(dataspace_id);
  dims[0] = 2;
  dataspace_id = H5Screate_simple( 1, dims, NULL );
  dataset_id = H5Dcreate(file_id, "/image_center",out_type_id ,
			 dataspace_id, H5P_DEFAULT);
  values[0] = img->detector->image_center[0];
  values[1] = img->detector->image_center[1];
  status = H5Dwrite(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, values);
  status = H5Dclose(dataset_id);
  status = H5Sclose(dataspace_id);

  dims[0] = 1;
  dataspace_id = H5Screate_simple( 1, dims, NULL );
  values[0] = img->phased;
  dataset_id = H5Dcreate(file_id, "/phased", out_type_id,
			 dataspace_id, H5P_DEFAULT);
  status = H5Dwrite(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, values);
  status = H5Dclose(dataset_id);

  values[0] = img->shifted;
  dataset_id = H5Dcreate(file_id, "/shifted", out_type_id,
			 dataspace_id, H5P_DEFAULT);
  status = H5Dwrite(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, values);
  status = H5Dclose(dataset_id);

  values[0] = img->detector->lambda;
  dataset_id = H5Dcreate(file_id, "/lambda", out_type_id,
			 dataspace_id, H5P_DEFAULT);
  status = H5Dwrite(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, values);
  status = H5Dclose(dataset_id);

  values[0] = img->detector->pixel_size;
  dataset_id = H5Dcreate(file_id, "/pixel_size", out_type_id,
			 dataspace_id, H5P_DEFAULT);
  status = H5Dwrite(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, values);
  status = H5Dclose(dataset_id);

  values[0] = img->detector->detector_distance;
  dataset_id = H5Dcreate(file_id, "/detector_distance", out_type_id,
			 dataspace_id, H5P_DEFAULT);
  status = H5Dwrite(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, values);
  status = H5Dclose(dataset_id);

  values[0] = img->scaled;
  dataset_id = H5Dcreate(file_id, "/scaled", out_type_id,
			 dataspace_id, H5P_DEFAULT);
  status = H5Dwrite(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, values);
  status = H5Dclose(dataset_id);
  status = H5Sclose(dataspace_id);


  status = H5Fclose(file_id);
}

Image * read_imagefile(const char * filename){
  Image * res = malloc(sizeof(Image));
  int file_id,dataset_id,space;
  int status;
  hsize_t dims[2];
  hid_t mem_type_id = 0;
  real values[2];
  if(sizeof(real) == sizeof(float)){
    mem_type_id = H5T_NATIVE_FLOAT;
  }else if(sizeof(real) == sizeof(double)){
    mem_type_id = H5T_NATIVE_DOUBLE;
  }else{
    abort();
  }

  res->detector = malloc(sizeof(Detector));

  file_id = H5Fopen(filename,H5F_ACC_RDONLY,H5P_DEFAULT);
  dataset_id = H5Dopen(file_id, "/mask");
  space = H5Dget_space(dataset_id);
  H5Sget_simple_extent_dims(space,dims,NULL);
  res->detector->size[0] = dims[0];
  res->detector->size[1] = dims[1];  

  res->mask = malloc(sizeof(real)*dims[0]*dims[1]);
  status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		   H5P_DEFAULT, res->mask);
  status = H5Dclose(dataset_id);

  dataset_id = H5Dopen(file_id, "/image_center");
  status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		   H5P_DEFAULT, values);
  status = H5Dclose(dataset_id);
  res->detector->image_center[0] = values[0];
  res->detector->image_center[1] = values[1];

  dataset_id = H5Dopen(file_id, "/phased");
  status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		   H5P_DEFAULT, values);
  status = H5Dclose(dataset_id);
  res->phased = values[0];

  dataset_id = H5Dopen(file_id, "/shifted");
  status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		   H5P_DEFAULT, values);
  status = H5Dclose(dataset_id);
  res->shifted = values[0];

  dataset_id = H5Dopen(file_id, "/scaled");
  status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		   H5P_DEFAULT, values);
  status = H5Dclose(dataset_id);
  res->scaled = values[0];

  dataset_id = H5Dopen(file_id, "/detector_distance");
  status = H5Dread(dataset_id,  mem_type_id, H5S_ALL, H5S_ALL,
		   H5P_DEFAULT, values);
  status = H5Dclose(dataset_id);
  res->detector->detector_distance = values[0];

  dataset_id = H5Dopen(file_id, "/lambda");
  status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		   H5P_DEFAULT, values);
  status = H5Dclose(dataset_id);
  res->detector->lambda = values[0];

  dataset_id = H5Dopen(file_id, "/pixel_size");
  status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		   H5P_DEFAULT, values);
  status = H5Dclose(dataset_id);
  res->detector->pixel_size = values[0];


  if(res->phased){
    res->c = malloc(sizeof(real)*dims[0]*dims[1]);
    res->r = malloc(sizeof(real)*dims[0]*dims[1]);
    dataset_id = H5Dopen(file_id, "/complex");
    status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		     H5P_DEFAULT, res->c);
    status = H5Dclose(dataset_id);

    dataset_id = H5Dopen(file_id, "/real");
    status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		     H5P_DEFAULT, res->r);
    status = H5Dclose(dataset_id);
  }

  if(res->scaled){
    dataset_id = H5Dopen(file_id, "/amplitudes");
    space = H5Dget_space(dataset_id);
    H5Sget_simple_extent_dims(space,dims,NULL);
    res->detector->size[0] = dims[0];
    res->detector->size[1] = dims[1];  
    res->amplitudes = malloc(sizeof(real)*dims[0]*dims[1]);  
    status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		     H5P_DEFAULT, res->amplitudes);
    status = H5Dclose(dataset_id);    
    res->image = res->amplitudes;
  }else{
    dataset_id = H5Dopen(file_id, "/intensities");
    space = H5Dget_space(dataset_id);
    H5Sget_simple_extent_dims(space,dims,NULL);
    res->detector->size[0] = dims[0];
    res->detector->size[1] = dims[1];  
    res->intensities = malloc(sizeof(real)*dims[0]*dims[1]);  
    status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		     H5P_DEFAULT, res->intensities);
    status = H5Dclose(dataset_id);
    res->image = res->intensities;
  }
  status = H5Fclose(file_id);
  return res;
  
}


Image * read_tiff(const char * filename){
  Image * res = malloc(sizeof(Image));
  float * img;
  int nstrips;
  int stripsize;
  TIFF * tif;
  int i,x,y;;


  res->detector = malloc(sizeof(Detector));

  tif = TIFFOpen(filename, "r");  
  if(!tif){
    return NULL;
  }
  
  nstrips = TIFFNumberOfStrips(tif);
  stripsize = TIFFStripSize(tif);
  img = malloc(nstrips*stripsize);
  for(i = 0;i<nstrips;i++){
    TIFFReadEncodedStrip(tif,i,&(img[i*stripsize/4]),stripsize);
  }
  TIFFClose(tif);
  res->detector->size[0] = stripsize/sizeof(float);
  res->detector->size[1] = nstrips;  
  res->scaled = 0;
  res->intensities = malloc(sizeof(real)*nstrips*stripsize/sizeof(float));
  /* Transpose image */
  for(x = 0;(unsigned int)x<(stripsize/sizeof(float));x++){
    for(y = 0;y<nstrips;y++){
      res->intensities[x*nstrips+y] = img[y*(stripsize/sizeof(float))+x];
    }
    
  }
  free(img);
  
  res->mask = malloc(sizeof(real)*TSIZE(res));
  for(i = 0;i<TSIZE(res);i++){
    res->mask[i] = 1;
  }
  res->detector->image_center[0] =  (res->detector->size[0]-1)/2.0;
  res->detector->image_center[1] = (res->detector->size[1]-1)/2.0;
  res->phased = 0;
  res->image = res->intensities;
  return res;  
}


void write_tiff(Image * img, char * filename){
  float * data;
  int nstrips;
  int stripsize;
  TIFF * tif;
  int x,y;;



  tif = TIFFOpen(filename, "w");  
  nstrips = img->detector->size[1];
  stripsize = img->detector->size[0]*sizeof(float);
  TIFFSetField(tif,TIFFTAG_IMAGEWIDTH,stripsize/sizeof(float));
  TIFFSetField(tif,TIFFTAG_ROWSPERSTRIP,nstrips);
  TIFFSetField(tif,TIFFTAG_IMAGELENGTH,nstrips);
  TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 32);
  TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);
  TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
  TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
  TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
  data = malloc(nstrips*stripsize);

  /* Transpose image */
  for(x = 0;(unsigned int)x<(stripsize/sizeof(float));x++){
    for(y = 0;y<nstrips;y++){
      data[x+(stripsize/sizeof(float))*y] = img->image[x*img->detector->size[1]+y];
    }
    
  }

  TIFFWriteEncodedStrip(tif, 0, data, nstrips*stripsize);
  /*  for(i = 0;i<nstrips;i++){
    TIFFWriteEncodedStrip(tif,i,&(data[i*stripsize/4]),stripsize);
    }*/
  TIFFClose(tif);
  free(data);
}




/* b must fit inside a or the other way around */

/* If wrap arond is on, a is taken to be periodic and the outout size will
   be as big as the input 
   
   Otherwise the output size is equal to a+b-1

   If size is not a NULL pointer it will be the size of transform used.
   A size smaller than the size of a+b will results in wrap around effects.
   If size is NULL the size of a will be used (resulting in wrap around effects).

*/
Image * cross_correlate_img(Image * a, Image * b, int * size){
  int x;
  int y;  
  int i;
  Image * res;
  Image * a_ft;
  Image * b_ft;
  Image * tmp;

  if(a->detector->size[0] < b->detector->size[0]){
    /* swap */
    res = a;
    a = b;
    b = res;
  }

  if(size){
    x = size[0];
    y = size[1];
  }else{
    x = a->detector->size[0];
    y = a->detector->size[1];
  }

  tmp = zero_pad_image(a,x,y,1);
  a_ft = image_fft(tmp);
  freeimg(tmp);

  tmp = zero_pad_image(b,x,y,1);
  b_ft = image_fft(tmp);
  freeimg(tmp);

  tmp = create_empty_img(a_ft);
  tmp->shifted = 1;
  rephase(tmp);
  /* Now do the multiplication in fourier space */
  /* Using the Convolution Theorem */
  for(i = 0;i<x*y;i++){
    tmp->r[i] = a_ft->r[i]*b_ft->r[i]+a_ft->c[i]*b_ft->c[i];
    tmp->c[i] = -a_ft->r[i]*b_ft->c[i]+a_ft->c[i]*b_ft->r[i];
  }

  freeimg(a_ft);
  freeimg(b_ft);
  /* Backtransform */
  res = image_rev_fft(tmp);
  freeimg(tmp);

  /* should be all real */
  for(i = 0;i<TSIZE(res);i++){
    res->r[i] /= TSIZE(res);
    res->c[i] /= TSIZE(res);
    res->image[i] /= TSIZE(res);
  }

  return res;  
}

/* if wrap around is on b must fit inside a or the other way around */

/* Convolute only acts on the amplitudes 
   real and complex part will be set to 0 
*/
/* If size is not a NULL pointer it will be the size of transform used.
   A size smaller than the size of a+b will results in wrap around effects.
   If size is NULL the size of a will be used (resulting in wrap around effects).
*/
Image * convolute_img(Image * a, Image * b, int * size){
  int x;
  int y;  
  int i;
  Image * res;
  Image * a_ft;
  Image * b_ft;
  Image * tmp;

  if(a->detector->size[0] < b->detector->size[0]){
    /* swap */
    res = a;
    a = b;
    b = res;
  }
  
  if(!size){
    x = a->detector->size[0]/*+(b->detector->size[0]-1)/2*/;
    y = a->detector->size[1]/*+(b->detector->size[1]-1)/2*/;
  }else{
    x = size[0];
    y = size[1];
  }

  tmp = zero_pad_image(a,x,y,1);
  dephase(tmp);
  a_ft = image_fft(tmp);
  freeimg(tmp);

  tmp = zero_pad_image(b,x,y,1);
  dephase(tmp);
  b_ft = image_fft(tmp);
  freeimg(tmp);

  tmp = create_empty_img(a_ft);
  tmp->shifted = 1;
  rephase(tmp);
  /* Now do the multiplication in fourier space */
  /* Using the Convolution Theorem */
  for(i = 0;i<x*y;i++){
    tmp->r[i] = a_ft->r[i]*b_ft->r[i]-a_ft->c[i]*b_ft->c[i];
    tmp->c[i] = a_ft->r[i]*b_ft->c[i]+a_ft->c[i]*b_ft->r[i];
  }

  freeimg(a_ft);
  freeimg(b_ft);
  /* Backtransform */
  res = image_rev_fft(tmp);
  freeimg(tmp);

  dephase(res);
  /* should be all real */
  for(i = 0;i<TSIZE(res);i++){
    res->image[i] /= TSIZE(res);
  }

  return res;  
}


/* Convolute the image with a gaussian filter.
 The filter function is given by:

f(x,y) = 1/sqrt(2*M_PI*radius) * exp(-(x^2+y^2)/(2*radius^2)) */
Image * gaussian_blur(Image * in, real radius){
  /* Lets make this convolution using a fourier transform shallw we... good....*/
  int x,y;
  int i,j;
  int filter_side = ceil(radius)*3*2+1;
  real * filter = malloc(sizeof(real)*filter_side*filter_side);
  real total_filter = 0;
  Image * filter_img = create_empty_img(in);
  Image * centered_filter;
  Image * res;
  filter_img->detector->size[0] = filter_side;
  filter_img->detector->size[1] = filter_side;
  filter_img->detector->image_center[0] = (filter_side-1)/2.0;
  filter_img->detector->image_center[1] = (filter_side-1)/2.0;
  free(filter_img->image);
  filter_img->image = filter;
  
  if(filter_img->scaled){
    filter_img->amplitudes = filter_img->image;
  }else{
    filter_img->intensities = filter_img->image;
  }
  dephase(filter_img);
  for(x = -ceil(radius)*3;x<=ceil(radius)*3;x++){
    i = x+ceil(radius)*3;
    for(y = -ceil(radius)*3;y<=ceil(radius)*3;y++){
      j = y+ceil(radius)*3;
      filter[i*filter_side+j] = 1/sqrt(2*M_PI*radius) * exp(-(x*x+y*y)/(2*radius*radius));
      total_filter += filter[i*filter_side+j];
    }
  }
  for(i = 0;i<TSIZE(filter_img);i++){
    filter_img->image[i] /= total_filter;
  }
  centered_filter = shift_center_to_top_left(filter_img);
  centered_filter->shifted = 1;
  freeimg(filter_img);
  res = convolute_img(in, centered_filter,NULL);
  freeimg(centered_filter);
  return res;
}


/* Convolute the image with a square window.
 The filter function is given by:

f(x,y) = 1/((2*radius+1)^2)) */
Image * square_blur(Image * in, real radius){
  /* Lets make this convolution using a fourier transform shallw we... good....*/
  int x,y;
  int i,j;
  int filter_side = ceil(radius)*3*2+1;
  real * filter = malloc(sizeof(real)*filter_side*filter_side);
  real total_filter = 0;
  Image * filter_img = create_empty_img(in);
  Image * centered_filter;
  Image * res;
  filter_img->detector->size[0] = filter_side;
  filter_img->detector->size[1] = filter_side;
  filter_img->detector->image_center[0] = (filter_side-1)/2.0;
  filter_img->detector->image_center[1] = (filter_side-1)/2.0;
  free(filter_img->image);
  filter_img->image = filter;
  if(filter_img->scaled){
    filter_img->amplitudes = filter_img->image;
  }else{
    filter_img->intensities = filter_img->image;
  }
  dephase(filter_img);
  for(x = -ceil(radius)*3;x<=ceil(radius)*3;x++){
    i = x+ceil(radius)*3;
    for(y = -ceil(radius)*3;y<=ceil(radius)*3;y++){
      j = y+ceil(radius)*3;
      filter[i*filter_side+j] = 1.0/((2*radius+1)*(2*radius+1));
      total_filter += filter[i*filter_side+j];
    }
  }
  for(i = 0;i<TSIZE(filter_img);i++){
    filter_img->image[i] /= total_filter;
  }
  centered_filter = shift_center_to_top_left(filter_img);
  freeimg(filter_img);
  res = convolute_img(in, centered_filter,NULL);
  freeimg(centered_filter);
  return res;
}


int write_mask_to_png(Image * img, char * filename, int color){
  Image  * res = imgcpy(img);
  int ret;
  memcpy(res->image,res->mask,TSIZE(img)*sizeof(real));
  ret = write_png(res,filename,color);
  freeimg(res);
  return ret;
}


#ifdef _WIN32
  /* png_init_io seems to crash in windows using GnuWin32 libpng-1.2.8*/

#  define READFILE(file, data, length, check) \
     check=(png_size_t)fread(data,(png_size_t)1,length,file)
#  define WRITEFILE(file, data, length, check) \
     check=(png_size_t)fwrite(data,(png_size_t)1, length, file)
#define NEAR_BUF_SIZE 1024

static void
pngtest_write_data(png_structp png_ptr, png_bytep data, png_size_t length)
{
   png_uint_32 check;

   WRITEFILE((FILE *)png_ptr->io_ptr,  data, length, check);
   if (check != length)
   {
      png_error(png_ptr, "Write Error");
   }
}

static void
pngtest_read_data(png_structp png_ptr, png_bytep data, png_size_t length)
{
   png_size_t check;

   /* fread() returns 0 on error, so it is OK to store this in a png_size_t
    * instead of an int, which is what fread() actually returns.
    */
   READFILE((png_FILE_p)png_ptr->io_ptr, data, length, check);

   if (check != length)
   {
      png_error(png_ptr, "Read Error!");
   }
}
#endif 

Image * read_png(char * filename){
 FILE *fp = fopen(filename, "rb");
 int i,j,index;
 png_uint_32 width,height;
 int bit_depth,color_type,interlace_type,compression_type,filter_method;
 png_structp png_ptr = png_create_read_struct
   (PNG_LIBPNG_VER_STRING, (png_voidp)NULL,
    NULL,NULL);
 png_infop info_ptr = png_create_info_struct(png_ptr);
 png_infop end_info = png_create_info_struct(png_ptr);
 png_byte ** row_pointers;
 Image * res;
#ifdef _WIN32
 png_set_read_fn(png_ptr, (png_voidp)fp, pngtest_read_data);
#else
 png_init_io(png_ptr, fp);
#endif
 png_read_info(png_ptr, info_ptr);
 png_get_IHDR(png_ptr, info_ptr, &width, &height,
	      &bit_depth, &color_type, &interlace_type,
	      &compression_type, &filter_method);

 if(color_type == 6){
   bit_depth *= 4;
 }else if(color_type == 4){
   bit_depth *= 2;
 }else if(color_type == 2){
   bit_depth *= 3;
 }
 row_pointers = malloc(sizeof(png_byte *)*height);
 for(i = 0;i<height;i++){
   row_pointers[i] = malloc(sizeof(png_byte)*width*bit_depth/8);
 }
 png_read_image(png_ptr, row_pointers);
 res = create_new_img(width,height);
 for(i = 0;i<height;i++){
   for(j = 0;j<width;j++){
     res->image[j*height+i] = row_pointers[i][(int)(j*bit_depth/8)];
   }
 }
 return res;
}




int write_png(Image * img, char * filename, int color){

  FILE *fp = fopen(filename, "wb");
  
  png_structp png_ptr; 
  png_infop info_ptr;
  int bit_depth = 8;
  int color_type;
  int interlace_type = PNG_INTERLACE_NONE;
  int compression_type = PNG_COMPRESSION_TYPE_DEFAULT;
  int filter_method = PNG_FILTER_TYPE_DEFAULT;
  int png_transforms = PNG_TRANSFORM_IDENTITY/*|PNG_TRANSFORM_INVERT_MONO*/;
  int pixel_size = 0;
  int i,x,y;
  real log_of_2;
  real color_table[3][256];
  real scale,offset,max_v,min_v,value;
  png_byte ** row_pointers;

  /*fclose(fp);
  return 0;*/
  max_v = 0;
  min_v = REAL_MAX;

/* Fill color tables */

  for(i = 0;i<256;i++){
    value = i/255.0;
    if(color & COLOR_GRAYSCALE){       
      color_table[0][i] = value;
      color_table[1][i] = value;
      color_table[2][i] = value;
    }else if(color & COLOR_TRADITIONAL){       
      color_table[0][i] = sqrt(value);
      color_table[1][i] = value*value*value;
      color_table[2][i] = sin(value*2*M_PI);
    }else if(color & COLOR_HOT){
      color_table[0][i] = 3*value;
      color_table[1][i] = 3*value-1;
      color_table[2][i] = 3*value-2;
    }else if(color & COLOR_RAINBOW){
      color_table[0][i] = fabs(2*value-0.5);
      color_table[1][i] = sin(value*M_PI);
     color_table[2][i] = cos(value*M_PI/2);
    }else if(color & COLOR_JET){
      if(value < 1/8.0){
	color_table[0][i] = 0;
	color_table[1][i] = 0;
	color_table[2][i] = (value+1.0/8.0)*4;	   
      }else if(value < 3/8.0){
	color_table[0][i] = 0;
	color_table[1][i] = (value-1.0/8.0)*4;
	color_table[2][i] = 1;
      }else if(value < 5/8.0){
	color_table[0][i] = (value-3.0/8.0)*4;
	color_table[1][i] = 1;
	color_table[2][i] = 1-(value-3.0/8.0)*4;
      }else if(value < 7/8.0){
	color_table[0][i] = 1;
	color_table[1][i] = 1-(value-5.0/8.0)*4;
	color_table[2][i] = 0;
      }else if(value <= 1.01){
	color_table[0][i] = 1-(value-7.0/8.0)*4;;
	color_table[1][i] = 0;
	color_table[2][i] = 0;
      }
    }
    color_table[0][i] = MIN(1,color_table[0][i]);
    color_table[1][i] = MIN(1,color_table[1][i]);
    color_table[2][i] = MIN(1,color_table[2][i]);
    color_table[0][i] *= 255;
    color_table[1][i] *= 255;
    color_table[2][i] *= 255;
  }
  if (!fp){
    perror("Couldn't open file!\n");
    abort();
    return -1;
  }
  png_ptr = png_create_write_struct
    (PNG_LIBPNG_VER_STRING, (png_voidp)NULL/*user_error_ptr*/,
     NULL/*user_error_fn*/, NULL/*user_warning_fn*/);
  if (!png_ptr){
    perror("Couldn't allocate write structure!\n");
    abort();
    return -1;
  }
  info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr){
    png_destroy_write_struct(&png_ptr,
			     (png_infopp)NULL);
    perror("Couldn't allocate info structure!\n");
    abort();
    return (-1);
  }
  if (setjmp(png_jmpbuf(png_ptr))){
    png_destroy_write_struct(&png_ptr, &info_ptr);
    fclose(fp);
    perror("Couldn't setjmp!\n");
    abort();
    return (-1);
  }
  #ifdef _WIN32
    png_set_write_fn(png_ptr, (png_voidp)fp,  pngtest_write_data,NULL);
  #else
   png_init_io(png_ptr, fp);
  #endif

  
  color_type = PNG_COLOR_TYPE_RGB;
  /* 8 bits 3 channels */
  pixel_size = 3*1;
   /* png_set_compression_level(png_ptr,Z_BEST_COMPRESSION); */
  png_set_IHDR(png_ptr, info_ptr, img->detector->size[0], img->detector->size[1],
	       bit_depth, color_type, interlace_type,
	       compression_type, filter_method);
  
  

  row_pointers = png_malloc(png_ptr,img->detector->size[1]*sizeof(png_byte *));
  for (i=0; i<img->detector->size[1]; i++){
    row_pointers[i] = png_malloc(png_ptr,img->detector->size[0]*pixel_size*sizeof(png_byte));
  }
  
  /* We're gonna scale the image so that it fits on the 8 bits */
  for(i = 0;i<TSIZE(img);i++){
    if(img->image[i] < min_v){
      min_v = img->image[i];
    }
    if(img->image[i] > max_v){
      max_v = img->image[i];
    }
  }
  if(max_v-min_v){
    scale = 1/(max_v-min_v);
  }else{
    scale = 1;
  }
  offset = min_v;
  i = 0;
  log_of_2 = log(2.0);
  for(x = 0;x<img->detector->size[0];x++){
    for(y = 0;y<img->detector->size[1];y++){



      /* traditional color scale taken from gnuplot manual */
      if(color & LOG_SCALE){
	value = log((img->image[i]-offset)*scale+1)/log_of_2;
      }else{
	value = ((img->image[i]-offset)*scale);
      }
      value *= 255;
      row_pointers[y][x*3] =  color_table[0][(int)value];
      row_pointers[y][x*3+1] = color_table[1][(int)value];
      row_pointers[y][x*3+2] = color_table[2][(int)value];
      i++;
    }
  }
  png_set_rows(png_ptr, info_ptr, row_pointers);
  
  png_write_png(png_ptr, info_ptr, png_transforms, NULL);
  png_write_flush(png_ptr);
  /* png_write_end(png_ptr, info_ptr);*/
  for(i=0; i<img->detector->size[1]; i++){
    png_free(png_ptr,row_pointers[i]);
  }
  png_free(png_ptr,row_pointers);
  png_destroy_write_struct(&png_ptr, &info_ptr);
  fflush(fp);
  fclose(fp);
  return 0;
}

real r_factor(Image * fobs, Image *fcalc, real low_intensity_cutoff){
  real den = 0;
  real num = 0;
  int i;
  for(i = 0;i<TSIZE(fobs);i++){
    if(!fobs->mask[i] || fobs->image[i] < low_intensity_cutoff){
      continue;
    }
    num +=  fabs(fobs->image[i]-fcalc->image[i]);
    den += fobs->image[i];
  }
  return num/den;
}



/* Low pass filter using a centered square window of side edge_size */
Image * low_pass_square_filter(Image * in, int edge_size){
  Image * fft_img = image_fft(in);
  Image * res;
  Image * tmp;
  int i = 0;
  for(i = 0;i<TSIZE(in);i++){
    if(square_dist_to_center(i,in) > edge_size/2.0){
      fft_img->image[i] = 0;
      fft_img->r[i] = 0;
      fft_img->c[i] = 0;
    }
  }
  tmp = imgcpy(fft_img);
  for(i = 0;i<TSIZE(tmp);i++){
    tmp->image[i] = log(tmp->image[i]+1);
  }
  write_png(tmp,"low_pass.png",COLOR_JET);
  freeimg(tmp);

  res = image_rev_fft(fft_img);
  /* scale appropriately */
  for(i = 0;i<TSIZE(res);i++){
    res->image[i] /= TSIZE(res);
  }
  freeimg(fft_img);
  fft_img = image_fft(res);
  tmp = imgcpy(fft_img);
  for(i = 0;i<TSIZE(tmp);i++){
    tmp->image[i] = log(tmp->image[i]+1);
  }
  write_png(tmp,"after_low_pass.png",COLOR_JET);
  freeimg(tmp);

  
  return res;
}

/* Low pass filter using a centered gaussian window of side edge_size */
Image * low_pass_gaussian_filter(Image * in, int edge_size){
  Image * fft_img = image_fft(in);
  Image * res;
  Image * mask;
  int i = 0;
  gaussian_filter(fft_img,edge_size/2.0,1);
  res = image_rev_fft(fft_img);
  freeimg(fft_img);
  /* scale appropriately */
  for(i = 0;i<TSIZE(res);i++){
    res->image[i] /= TSIZE(res);
  }
  if(!in->phased){
    dephase(res);
  }
  /* Also low pass filter the mask */
  mask = imgcpy(in);
  memcpy(mask->image,in->mask,TSIZE(in)*sizeof(real)); 
  fft_img = image_fft(mask);
  gaussian_filter(fft_img,edge_size/2.0,1);
  freeimg(mask);
  mask = image_rev_fft(fft_img);
  freeimg(fft_img);
  /* scale appropriately */
  for(i = 0;i<TSIZE(mask);i++){
    /* if the mask is not really want then we have unkown information and we'll make it 0 */
    if(mask->image[i] < TSIZE(res)-1){
      res->mask[i] = 0;
    }else{
      res->mask[i] = 1;
    }
  }
  freeimg(mask);
  
  return res;
}

/* Filter using a centered gaussian window of side edge_size */
Image * gaussian_filter(Image * in, real radius,int in_place){
  Image * res;
  const real scaling_factor = 7;
  int i;
  real scaling;
  if(in_place){
    res = in;
  }else{
    res = imgcpy(in);
  }
  /* the scaling of 7 makes sure the pattern is scaled by 0.0009 at the edge */

  for(i = 0;i<TSIZE(in);i++){
    scaling = (dist_to_center(i,in)/(radius))*(dist_to_center(i,in)/(radius));
    res->image[i] *= exp(-scaling*scaling_factor);
    res->c[i] *= exp(-scaling*scaling_factor);
    res->r[i] *= exp(-scaling*scaling_factor);
  }
  return res;
}


Image * zero_pad_image(Image * a, int newx, int newy, int pad_mask){
  if(a->shifted){
    return zero_pad_shifted_image(a,newx,newy,pad_mask);
  }else{
    return zero_pad_unshifted_image(a,newx,newy,pad_mask);
  }
  return NULL;
}


/* Zero pads an image in the following way 

     1 1 2 2
in = 1 1 2 2 newx = 6 newy = 6
     3 3 4 4
     3 3 4 4

out = 1 1 0 0 2 2
      1 1 0 0 2 2
      0 0 0 0 0 0
      0 0 0 0 0 0
      3 3 0 0 4 4
      3 3 0 0 4 4

The in is split through the middle. In case of odd in dimensions
the upper left corner gets the extra row/column.

*/
static Image * zero_pad_shifted_image(Image * a, int newx, int newy,int pad_mask){
  Image * out;
  int x,y;
  int sourcex;
  int sourcey;
  if(newx < a->detector->size[0] || 
     newy < a->detector->size[1]){
    fprintf(stderr,"Negative padding!\n");
    abort();
  }else if(newx == a->detector->size[0] &&
	   newy == a->detector->size[1]){
    return imgcpy(a);
  }
  out = create_empty_img(a);
  out->detector->size[0] = newx;
  out->detector->size[1] = newy;
  free(out->image);
  out->image = malloc(sizeof(real)*TSIZE(out));
  free(out->mask);
  out->mask = malloc(sizeof(real)*TSIZE(out));
  if(out->scaled){
    out->amplitudes = out->image;
  }else{
    out->intensities = out->image;
  }
  if(out->phased){
    free(out->c);
    free(out->r);
    out->c = malloc(sizeof(real)*TSIZE(out));
    out->r = malloc(sizeof(real)*TSIZE(out));
  }

  for(x = 0;x<out->detector->size[0];x++){
    if(x < a->detector->size[0]/2.0){
      sourcex = x;
    }else if(out->detector->size[0]-x-1 < (a->detector->size[0]-1)/2.0){
      sourcex = (a->detector->size[0]-1)-(out->detector->size[0]-x-1);
    }else{
      sourcex = -1;
    }

    for(y = 0;y<out->detector->size[1];y++){
      if(y < a->detector->size[1]/2.0){
	sourcey = y;
      }else if(out->detector->size[1]-y-1 < (a->detector->size[1]-1)/2.0){
	sourcey = (a->detector->size[0]-1)-(out->detector->size[1]-y-1);
      }else{
	sourcey = -1;
      }
      if(sourcey == -1 || sourcex == -1){
	out->image[x*out->detector->size[1]+y] = 0;
	out->mask[x*out->detector->size[1]+y] = pad_mask;
	if(out->phased){
	  out->r[x*out->detector->size[1]+y] = 0;
	  out->c[x*out->detector->size[1]+y] = 0;
	}
      }else{
	out->image[x*out->detector->size[1]+y] = a->image[sourcex*a->detector->size[1]+sourcey];
	out->mask[x*out->detector->size[1]+y] = a->mask[sourcex*a->detector->size[1]+sourcey];
	if(out->phased){
	  out->r[x*out->detector->size[1]+y] = a->r[sourcex*a->detector->size[1]+sourcey];
	  out->c[x*out->detector->size[1]+y] = a->c[sourcex*a->detector->size[1]+sourcey];	  
	}
      }
    }
  }
  return out;
}


/* Zero pads an image in the following way 

     1 1 2 2
in = 1 1 2 2 newx = 6 newy = 6
     3 3 4 4
     3 3 4 4

out = 1 1 2 2 0 0
      1 1 2 2 0 0
      3 3 4 4 0 0
      3 3 4 4 0 0
      0 0 0 0 0 0
      0 0 0 0 0 0

*/
static Image * zero_pad_unshifted_image(Image * a, int newx, int newy, int pad_mask){
  Image * out;
  int x,y;
  int sourcex;
  int sourcey;
  if(newx < a->detector->size[0] || 
     newy < a->detector->size[1]){
    fprintf(stderr,"Negative padding!\n");
    abort();
  }else if(newx == a->detector->size[0] &&
	   newy == a->detector->size[1]){
    return imgcpy(a);
  }
  out = create_empty_img(a);
  out->detector->size[0] = newx;
  out->detector->size[1] = newy;
  free(out->image);
  out->image = malloc(sizeof(real)*TSIZE(out));
  free(out->mask);
  out->mask = malloc(sizeof(real)*TSIZE(out));
  if(out->scaled){
    out->amplitudes = out->image;
  }else{
    out->intensities = out->image;
  }
  if(out->phased){
    free(out->c);
    free(out->r);
    out->c = malloc(sizeof(real)*TSIZE(out));
    out->r = malloc(sizeof(real)*TSIZE(out));
  }

  for(x = 0;x<out->detector->size[0];x++){
    if(x < a->detector->size[0]){
      sourcex = x;
    }else{
      sourcex = -1;
    }
    
    for(y = 0;y<out->detector->size[1];y++){
      if(y < a->detector->size[1]){
	sourcey = y;
      }else{
	sourcey = -1;
      }
      if(sourcey == -1 || sourcex == -1){
	out->image[x*out->detector->size[1]+y] = 0;
	out->mask[x*out->detector->size[1]+y] = pad_mask;
	if(out->phased){
	  out->r[x*out->detector->size[1]+y] = 0;
	  out->c[x*out->detector->size[1]+y] = 0;
	}
      }else{
	out->image[x*out->detector->size[1]+y] = a->image[sourcex*a->detector->size[1]+sourcey];
	out->mask[x*out->detector->size[1]+y] = a->mask[sourcex*a->detector->size[1]+sourcey];
	if(out->phased){
	  out->r[x*out->detector->size[1]+y] = a->r[sourcex*a->detector->size[1]+sourcey];
	  out->c[x*out->detector->size[1]+y] = a->c[sourcex*a->detector->size[1]+sourcey];	  
	}
      }
    }
  }
  return out;
}



/* Shifts center to top left corner in the following way 

     1 1 2 2
in = 1 1 2 2 
     3 3 4 4
     3 3 4 4

out = 4 4 3 3
      4 4 3 3
      2 2 1 1
      2 2 1 1

Center is taken from image_center[]
*/
Image * shift_center_to_top_left(Image * a){
  Image * out = imgcpy(a);
  int x,y;
  int destx;
  int desty;
  for(x = 0;x<out->detector->size[0];x++){
    destx = (int)(x+out->detector->size[0]-out->detector->image_center[0])%out->detector->size[0];
    for(y = 0;y<out->detector->size[1];y++){
      desty = (int)(y+out->detector->size[1]-out->detector->image_center[1])%out->detector->size[1];
      out->image[destx*out->detector->size[1]+desty] = a->image[x*a->detector->size[1]+y];
      if(out->phased){
	out->r[destx*out->detector->size[1]+desty] = a->r[x*a->detector->size[1]+y];
	out->c[destx*out->detector->size[1]+desty] = a->c[x*a->detector->size[1]+y];	  
      }
    }
  }
  return out;
}


/* Csiszar I-divergence */
real I_divergenge(Image * a, Image * b){
  double suma = 0;
  double sumb = 0;
  int i;
  for(i = 0;i<TSIZE(a);i++){
    suma += a->image[i] - b->image[i];
    if(b->image[i]){
      sumb += b->image[i] * log(b->image[i]/(a->image[i]+FLT_EPSILON));
    }
  }
  return suma+sumb;
}


real integrated_intensity(Image * a){
  double sum = 0;
  int i;
  for(i = 0;i<TSIZE(a);i++){
    sum += a->image[i];
  }
  return sum;
}

int write_vtk(Image * img, char * filename){
  FILE * f = fopen(filename,"w");
  int i,x,y;
  if(!f){
    perror("Bad file in write_vtk!");
    abort();
  }
  fprintf(f,"# vtk DataFile Version 2.0\n");
  fprintf(f,"Generated by image_util write_vtk()\n");
  fprintf(f,"ASCII\n");
  fprintf(f,"DATASET STRUCTURED_POINTS\n");
  fprintf(f,"DIMENSIONS %d %d 1\n",img->detector->size[0],img->detector->size[1]);
  fprintf(f,"ORIGIN 0 %d 0\n",img->detector->size[1]);
  fprintf(f,"SPACING 1 -1 1\n");
  fprintf(f,"POINT_DATA %d\n",TSIZE(img));
  fprintf(f,"SCALARS amplitudes float 1\n");
  fprintf(f,"LOOKUP_TABLE default\n");
  fprintf(f,"%6g",img->image[0]);
  for(y = 0;y<img->detector->size[1];y++){
    for(x = 0;x<img->detector->size[0];x++){
      fprintf(f," %g",img->image[x*img->detector->size[1]+y]);    
    }
  }
/*  for(i = 1;i<TSIZE(img);i++){
    fprintf(f," %g",img->image[i]);    
  }*/
  fprintf(f,"\n");
  fflush(f);
  fclose(f);
  return 0;
}


Image * rectangle_crop(Image * in, int x1, int y1, int x2, int y2){
  Image * cropped;
  int i;
  /* x1,y1 should be upper left, x2,y2 lower right */
  if(x1 > x2 || y1 > y2){
    return NULL;
  }
  cropped = create_empty_img(in);
  cropped->detector->image_center[0] -= x1;
  cropped->detector->image_center[1] -= y1;
  cropped->detector->size[0] = x2-x1+1;
  cropped->detector->size[1] = y2-y1+1;
  free(cropped->image);
  cropped->image = malloc(sizeof(real)*TSIZE(cropped));
  free(cropped->mask);
  cropped->mask = malloc(sizeof(real)*TSIZE(cropped));
  if(cropped->scaled){
    cropped->amplitudes = cropped->image;
  }else{
    cropped->intensities = cropped->image;
  }
  for(i = x1;i<= x2;i++){
    memcpy(&cropped->image[(i-x1)*cropped->detector->size[1]],&in->image[(i)*in->detector->size[1]+y1],cropped->detector->size[1]*sizeof(real));
    memcpy(&cropped->mask[(i-x1)*cropped->detector->size[1]],&in->mask[(i)*in->detector->size[1]+y1],cropped->detector->size[1]*sizeof(real));
    if(in->phased){
      memcpy(&cropped->r[(i-x1)*cropped->detector->size[1]],&in->r[(i)*in->detector->size[1]+y1],cropped->detector->size[1]*sizeof(real));
      memcpy(&cropped->c[(i-x1)*cropped->detector->size[1]],&in->c[(i)*in->detector->size[1]+y1],cropped->detector->size[1]*sizeof(real));
    }
  }
  return cropped;
}


void find_center(Image * img, real * center_x, real * center_y){
  int x,y;
  float bx = -1;
  float by = -1;
  Image * a = centrosym_convolve(img);
  real max = 0;
  for(x = 0;x<img->detector->size[0];x++){
    for(y = 0;y<img->detector->size[1];y++){
      if(a->image[x*img->detector->size[1]+y] > max){
	max = a->image[x*img->detector->size[1]+y];
	bx = x;
	by = y;
      }      
    }
  }
  if(bx < img->detector->size[0]/2.0){
    bx = (img->detector->size[0])/2.0+bx/2.0;
  }else{
    bx = (img->detector->size[0])/2.0-(img->detector->size[0]-bx)/2.0;
  }
  if(by < img->detector->size[1]/2.0){
    by = (img->detector->size[1])/2.0+by/2.0;
  }else{
    by = (img->detector->size[1])/2.0-(img->detector->size[1]-by)/2.0;
  }
  printf("Center x - %f y - %f\n",bx,by);
  write_png(a,"centrosym_convolve.png",COLOR_JET|LOG_SCALE);
  *center_x  = bx;
  *center_y = by;
  freeimg(a);
}

int pixel_to_index(Image * img, real * point){
  return ((int)point[0])*img->detector->size[1]+point[1];
}


/* This doesn't really rescale the mask which is a problem and doesn't really downscale */
Image * fourier_rescale(Image * img, int new_x, int new_y){
  Image * res = image_fft(img);
  Image * tmp;
  int i;
  real inv_size;
  if(new_x < img->detector->size[0] ||
     new_y < img->detector->size[1]){
    perror("fourier_scale doesn't downscale");
    abort();
  }
  tmp = zero_pad_image(res,new_x,new_y,1);
  freeimg(res);
  res = image_rev_fft(tmp);
  freeimg(tmp);
  res->detector->image_center[0] = img->detector->image_center[0]*(new_x/img->detector->size[0]);
  res->detector->image_center[1] = img->detector->image_center[1]*(new_x/img->detector->size[1]);
  inv_size = 1.0/TSIZE(img);
  for(i = 0;i<TSIZE(res);i++){
    res->image[i] *= inv_size;
  }
  if(res->phased){
    for(i = 0;i<TSIZE(res);i++){
      res->r[i] *= inv_size;
      res->c[i] *= inv_size;
    }    
  }
  return res;
}


Image * bilinear_rescale(Image * img, int new_x, int new_y){
  Image * res = imgcpy(img);

  real virtual_x;
  real virtual_y;
  int x,y;
  res->detector->size[0] = new_x;
  res->detector->size[1] = new_y;
  res->detector->image_center[0] *= new_x/img->detector->size[0];
  res->detector->image_center[1] *= new_y/img->detector->size[1];
  res->image = realloc(res->image,TSIZE(res)*sizeof(real));
  if(res->scaled){
    res->amplitudes = res->image;
  }else{
    res->intensities = res->image;
  }
  res->mask = realloc(res->mask,TSIZE(res)*sizeof(real));
  if(res->phased){
    res->r = realloc(res->r,TSIZE(res)*sizeof(real));
    res->c = realloc(res->c,TSIZE(res)*sizeof(real));
  }
  

  for(x = 0; x<res->detector->size[0]; x++){
    virtual_x = (real)x*img->detector->size[0]/res->detector->size[0];
    for(y = 0; y<res->detector->size[1]; y++){
      virtual_y = (real)y*img->detector->size[1]/res->detector->size[1];
      if(res->phased){
	res->r[x*res->detector->size[1]+y] = bilinear_interpol_img(img,img->r,virtual_x,virtual_y);
	res->c[x*res->detector->size[1]+y] = bilinear_interpol_img(img,img->c,virtual_x,virtual_y);
	res->image[x*res->detector->size[1]+y] = bilinear_interpol_img(img,img->image,virtual_x,virtual_y);
      }else{
	res->image[x*res->detector->size[1]+y] = bilinear_interpol_img(img,img->image,virtual_x,virtual_y);
      }
      res->mask[x*res->detector->size[1]+y] = bilinear_interpol_img(img,img->mask,virtual_x,virtual_y);
      /* If we don't have the full information we don't have any information */
      if(res->mask[x*res->detector->size[1]+y] < 0.99){
	res->mask[x*res->detector->size[1]+y] = 0;
      }else{
	res->mask[x*res->detector->size[1]+y] = 1;
      }
    }
  } 
  return res;
}
  

real bilinear_interpol_img(Image * img, real * data, real v_x, real v_y){
  int x = v_x;
  int y = v_y;
  real u = v_x-x;
  real v = v_y-y;
  real res = 0;
  if(x >= img->detector->size[0]-1){
    x = img->detector->size[0]-2;
    u = 1;
  }
  if(y >= img->detector->size[1]-1){
    y = img->detector->size[1]-2;
    v = 1;
  }
  res = data[x*img->detector->size[1]+y]*(1-u)*(1-v)+
    data[(x+1)*img->detector->size[1]+y]*(u)*(1-v)+
    data[x*img->detector->size[1]+(y+1)]*(1-u)*(v)+
    data[(x+1)*img->detector->size[1]+(y+1)]*(u)*(v);
  return res;
}

void multiply_scalar_by_image(Image * img, real value){
  int i;
  for(i = 0;i<TSIZE(img);i++){
    img->image[i] *= value;
    if(img->phased){
      img->r[i] *= value;
      img->c[i] *= value;
    }
  }
}


int get_image_maximum(Image * img, int * x, int * y, real * max){
  real tmp = -1;
  int i;
  int index = -1;
  for(i = 0;i<TSIZE(img);i++){
    if(img->image[i] > tmp){
      tmp = img->image[i];
      index = i;
    }
  }
  *max = tmp;
  *y = index%img->detector->size[1];
  *x = index/img->detector->size[1];
  return index;
}


void resize_empty_image(Image * img, int new_x, int new_y){
  Image * res = img;

  res->detector->size[0] = new_x;
  res->detector->size[1] = new_y;
  res->detector->image_center[0] *= new_x/img->detector->size[0];
  res->detector->image_center[1] *= new_y/img->detector->size[1];
  res->image = realloc(res->image,TSIZE(res)*sizeof(real));
  if(res->scaled){
    res->amplitudes = res->image;
  }else{
    res->intensities = res->image;
  }
  res->mask = realloc(res->mask,TSIZE(res)*sizeof(real));
  if(res->phased){
    res->r = realloc(res->r,TSIZE(res)*sizeof(real));
    res->c = realloc(res->c,TSIZE(res)*sizeof(real));
  }
}

Image * rectangular_window(int image_x, int image_y, int width, int height, int shifted){
  Image * res = create_new_img(image_x,image_y);
  int x,y,i;
  int center[2];
  i = 0;
  
  if(shifted){
    for(x = 0;x<image_x;x++){
      for(y = 0;y<image_y;y++){
	if((fabs(x) < width/2 || fabs(image_x-x) < width/2 )&&
	   (fabs(y) < height/2 || fabs(image_y-y) < height/2 )){
	  res->image[i] = 1;
	}else{
	  res->image[i] = 0;
	}
	i++;
      }
    }
  }else{
    center[0] = image_x/2;
    center[1] = image_y/2;
    for(x = 0;x<image_x;x++){
      for(y = 0;y<image_y;y++){
	if(fabs(x-center[0]) < width/2 &&
	   fabs(y-center[1]) < height/2){
	  res->image[i] = 1;
	}else{
	  res->image[i] = 0;
	}
	i++;
      }
    }
  }
  rephase(res);
  return res;
}


Image * circular_window(int x, int y, int radius, int shifted){
  Image * res = create_new_img(x,y);
  int i;
  if(shifted){
    res->detector->image_center[0] = 0;
    res->detector->image_center[1] = 0;
  }else{
    res->detector->image_center[0] = x/2;
    res->detector->image_center[1] = y/2;
  }
  for(i = 0;i<x*y;i++){
    if(dist_to_center(i,res) < radius){
      res->image[i] = 1;
    }else{
      res->image[i] = 0;
    }
  }
  rephase(res);
  return res;
}


Image * normalize_image(Image * in){
  Image * res = imgcpy(in);
  double integral = 0;
  int i;
  for(i = 0;i<TSIZE(in);i++){
    integral += in->image[i];
  }
  for(i = 0;i<TSIZE(in);i++){
    res->image[i] /= integral;
  }
  if(res->phased){
    res->r[i] /= integral;
    res->c[i] /= integral;
  }
  return res;
}


Image * get_mask_from_image(Image * a){
  Image * res = imgcpy(a);
  int i;
  for(i = 0;i<TSIZE(a);i++){
    res->image[i] = res->mask[i];
  }
  return res;
}

real point_convolute_img(Image * a, Image * b, int index){
  int index_x, index_y,x,y;
  double out = 0;
  int ai,bi;
  Image * tmp;
  if(a->shifted){
    tmp = shift_quadrants(a);
    index = quadrant_shift_index(a,index);
    a = tmp;
  }
  if(b->shifted){
    fprintf(stderr,"Point convoluting with a shifted function is not currently defined!\n");
    return 0;
  }
  index_x = index/a->detector->size[1]-a->detector->image_center[0];
  index_y = index%a->detector->size[1]-a->detector->image_center[1];
  for(x = -b->detector->image_center[0];x<b->detector->size[0]-b->detector->image_center[0];x++){
    for(y = -b->detector->image_center[1];y<b->detector->size[1]-b->detector->image_center[1];y++){
      if(x+index_x < -a->detector->image_center[0] || 
	 x+index_x >=  a->detector->size[0]-a->detector->image_center[0] ||
	 y+index_y < -a->detector->image_center[1] || 
	 y+index_y >=  a->detector->size[1]-a->detector->image_center[1]){
	/* we're outside of image a */
	continue;
      }
      ai = (index_x+x+a->detector->image_center[0])*a->detector->size[1]+index_y+y+a->detector->image_center[1];
      bi = (x+b->detector->image_center[0])*b->detector->size[1]+y+b->detector->image_center[1];
      out += a->image[ai]*b->image[bi];
    }
  }
  if(a->shifted){
    freeimg(tmp);
  }
  return out;
}

int quadrant_shift_index(Image * a, int index){
  int x = index/a->detector->size[1];
  int y = index%a->detector->size[1];
  if(a->shifted){
    x = (x+a->detector->size[0]/2)%a->detector->size[0];
    y = (y+a->detector->size[1]/2)%a->detector->size[1];
  }else{
    x -= a->detector->image_center[0];
    y -= a->detector->image_center[1];
    x += a->detector->size[0];
    x %= a->detector->size[0];
    y += a->detector->size[1];
    y %= a->detector->size[1];
  }
  return x*a->detector->size[1]+y;
}

/*
  Returns the value of the distance between
  "point" (Or the img center in case "point" is NULL)
  and the border of the image in a certain direction
  specified by an angle (in radians)
  The point of intersection is returned
  in "intersection".
*/

real get_image_radial_distance_to_border(Image * img, real * point, real direction, real * intersection){
  real d1,d2;
  real x1,y1,x2,y2;
  real * center;
  int * dim;
  if(point){
    center = point;
  }else{
    center = img->detector->image_center;
  }
  dim = img->detector->size;
  if(cos(direction)>0){
    /* d1 assumes that the intersection is with the right border */
    y1 = center[1]-tan(direction)*(dim[0]-center[0]);
    x1 = dim[0]-1; 
  }else{
    /* d1 assumes that the intersection is with the left border */
    y1 = center[1]+(tan(direction)*(center[0]));
    x1 = 0; 
  }
  d1 = sqrt((x1-center[0])*(x1-center[0])+(y1-center[1])*(y1-center[1]));
  if(sin(direction) > 0){
    /* d2 assumes that the intersection is with the upper border */
    x2 = center[0]+(1.0/tan(direction))*(center[1]);
    y2 = 0;
  }else{
    /* d2 assumes that the intersection is with the lower border */
    x2 = center[0]-(1.0/tan(direction))*(dim[1]-center[1]);
    y2 = dim[1]-1;
  }
  d2 = sqrt((x2-center[0])*(x2-center[0])+(y2-center[1])*(y2-center[1]));
  if(d1 < d2){
    if(intersection){
      intersection[0] = x1;
      intersection[1] = y1;
    }
    return d1;
  }else{
    if(intersection){
      intersection[0] = x2;
      intersection[1] = y2;
    }
    return d2;
  }              
}

/* 
   returns the values of the image along a radial
   sector in a certain direction(in radians)
   samples defines how many values we return
   
   The origin of the vector is given by point,
   or by the image center if point is NULL.

   intersection old the position of the image where the 
   sector vector intersect the image border.

   The values are obtained by linear interpolation
 */
Image * get_image_radial_sector(Image * img, real * point, real direction, int samples, real * intersection){
  int i;
  Image  *ret = create_new_img(samples,1);
  real fpixel[2];
  real d_to_border;
  int px0,py0;
  real t,u;
  real * center;
  int * dim;
  if(point){
    center = point;
  }else{
    center = img->detector->image_center;
  }
  dim = img->detector->size;
  d_to_border = get_image_radial_distance_to_border(img,center,direction, intersection);
  for(i = 0;i<samples;i++){
    fpixel[0] = center[0]+cos(direction)*d_to_border*((real)i/samples);
    /* The - happens because the pixel 0,0 on an image is the upper left not the lower left */
    fpixel[1] = center[1]-sin(direction)*d_to_border*((real)i/samples);
    /* bilinear interpolation around fpixel */
    px0 = (int)fpixel[0];
    py0 = (int)fpixel[1];
    t = (fpixel[0]-px0);
    u = (fpixel[1]-py0);
    /* border problems */
    if(px0 >= img->detector->size[0]-1){
      px0 = img->detector->size[0]-2;
      t = 1;
    }
    if(py0 >= img->detector->size[1]-1){
      py0 = img->detector->size[1]-2;
      u = 1;
    }
    ret->image[i] =  (1.0-t) * (1.0-u) * img->image[px0*dim[1]+py0]+
	  t * (1.0-u) * img->image[(px0+1)*dim[1]+py0]+
	  t * u * img->image[(px0+1)*dim[1]+py0+1]+
	  (1.0-t)* u * img->image[(px0)*dim[1]+py0+1];
    ret->mask[i] =  (1.0-t) * (1.0-u) * img->mask[px0*dim[1]+py0]+
	  t * (1.0-u) * img->mask[(px0+1)*dim[1]+py0]+
	  t * u * img->mask[(px0+1)*dim[1]+py0+1]+
	  (1.0-t)* u * img->mask[(px0)*dim[1]+py0+1];
  }
  return ret;  
}


/*
  Takes a sector and rotates it around the center to create an image.
  It assumes the sector is pixel scaled (each bin 1 pixel ).
  That is it does not try to stretch the sector to cover the image.
*/

Image * circular_image_from_sector(Image * sector, int * img_size, real * center){
  int x,y;
  Image * ret = create_new_img(img_size[0],img_size[1]);
  int bin;
  real d;
  real u;
  for(x = 0;x<img_size[0];x++){
    for(y = 0;y<img_size[1];y++){
      /* calculate distance to the center and sample from sector accordingly */
      d = sqrt((x-center[0])*(x-center[0])+(y-center[1])*(y-center[1]));
      bin = d;
      u = d-bin;
      if(d+1 < sector->detector->size[0]-1){
	ret->image[x*img_size[1]+y] = (1.0-u)*sector->image[bin];
	ret->image[x*img_size[1]+y] += (u)*sector->image[bin+1];
	ret->mask[x*img_size[1]+y] = 1;
      }else{
	ret->image[x*img_size[1]+y] = 0;
	ret->mask[x*img_size[1]+y] = 0;
      }
    }
  }
  return ret;
}



Image * image_local_variance(Image * img, Image * window){
  int i;
  int size[2] = {img->detector->size[0]+window->detector->size[0]-1,img->detector->size[1]+window->detector->size[1]-1};
  Image * norm_window = normalize_image(window);
  Image * ra = convolute_img(img,norm_window,size);
  Image * res = rectangle_crop(ra,0,0,img->detector->size[0]-1,img->detector->size[1]-1);
  write_png(img,"non_averaged.png",COLOR_JET);
  write_png(ra,"total_averaged.png",COLOR_JET);
  write_png(res,"crop_total_averaged.png",COLOR_JET);
  for(i = 0;i<TSIZE(res);i++){
    res->image[i] = fabs(res->image[i]-img->image[i]);
  }
  freeimg(ra);
  ra = convolute_img(res,norm_window,size);
  freeimg(res);
  res = rectangle_crop(ra,0,0,img->detector->size[0]-1,img->detector->size[1]-1);
  freeimg(ra);
  freeimg(norm_window);

  return res;
}



int sp_image_dot_prod(Image * a, Image * b, real * r, real * j){
  int i;
  *r = 0;
  *j = 0;
  if(TSIZE(a) != TSIZE(b)){
    return -1;
  }
  if(a->phased && b->phased){
    for(i = 0;i<TSIZE(a);i++){
      *r += a->r[i]*b->r[i] - a->c[i]*b->c[i];
      *j += a->c[i]*b->r[i] + a->r[i]*b->c[i];
    }
  }else if(a->phased && !b->phased){
    for(i = 0;i<TSIZE(a);i++){
      *r += a->r[i]*b->image[i];
      *j += a->c[i]*b->image[i];
    }
  }else if(!a->phased && b->phased){
    for(i = 0;i<TSIZE(a);i++){
      *r += b->r[i]*a->image[i];
      *j += b->c[i]*a->image[i];
    }
  }else if(!a->phased && !b->phased){
    for(i = 0;i<TSIZE(a);i++){
      *r += b->image[i]*a->image[i];
    }
  }else{
    fprintf(stderr,"Cannot get here!\n");
    return -1;
  }
  return 0;
}

Image * sp_proj_module(Image * a, Image * b){
  Image * ret = imgcpy(a);
  int i;
  for(i = 0;i<TSIZE(a);i++){
    ret->r[i] = b->image[i]*ret->r[i]/ret->image[i];
    ret->c[i] = b->image[i]*ret->c[i]/ret->image[i];
    ret->image[i] = b->image[i];
  }
  return ret;
}


Image * sp_proj_support(Image * a, Image * b){
  Image * ret = imgcpy(a);
  int i;
  for(i = 0;i<TSIZE(a);i++){
    ret->r[i] *= b->image[i];
    ret->c[i] *= b->image[i];
    ret->image[i] *= b->image[i];
  }
  return ret;
}


int sp_image_invert(Image * a){
  int i;
  if(a->phased){
    for(i = 0;i<TSIZE(a);i++){
      a->r[i] = 1/a->r[i];
      a->c[i] = 1/a->c[i];
      a->image[i] = norm(a,i);
    }
  }else{
    for(i = 0;i<TSIZE(a);i++){
      a->image[i] = 1/a->image[i];
    }
  }
  return 0;
}


Image * sp_image_mul_elements(Image * a, Image * b){
  Image * ret = NULL;
  int i;
  if(a->phased && b->phased){
    ret = imgcpy(a);
    for(i = 0;i<TSIZE(a);i++){
      ret->r[i] = a->r[i]*b->r[i] - a->c[i]*b->c[i];
      ret->c[i] = a->c[i]*b->r[i] + a->r[i]*b->c[i];
      ret->image[i] = norm(ret,i);
    }
  }else if(a->phased && !b->phased){
    ret = imgcpy(a);
    for(i = 0;i<TSIZE(a);i++){
      ret->r[i] = a->r[i]*b->image[i];
      ret->c[i] = a->c[i]*b->image[i];
      ret->image[i] *= b->image[i];
    }
  }else if(!a->phased && b->phased){
    ret = imgcpy(b);
    for(i = 0;i<TSIZE(a);i++){
      ret->r[i] = b->r[i]*a->image[i];
      ret->c[i] = b->c[i]*a->image[i];
      ret->image[i] *= a->image[i];
    }
  }else if(!a->phased && !b->phased){
    ret = imgcpy(a);
    for(i = 0;i<TSIZE(a);i++){
      ret->image[i] = b->image[i]*a->image[i];
    }
  }else{
    fprintf(stderr,"Cannot get here!\n");
    return NULL;
  }
  return ret;
}

void sp_image_complex_conj(Image * a){  
  int i;
  if(a->phased){
    for(i = 0;i<TSIZE(a);i++){
      a->c[i] -= a->c[i];
    }
  }
}


void sp_image_memcpy(Image * dst,Image * src){  
  if(src->phased && dst->phased){
    memcpy(dst->c,src->c,sizeof(real)*TSIZE(src));
    memcpy(dst->r,src->r,sizeof(real)*TSIZE(src));
  }
  memcpy(dst->mask,src->mask,sizeof(real)*TSIZE(src));
  memcpy(dst->image,src->image,sizeof(real)*TSIZE(src));
  memcpy(dst->detector,src->detector,sizeof(Detector));
}
