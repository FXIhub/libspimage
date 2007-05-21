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
#include <ctype.h>
#include <strings.h>
#include "spimage.h"


static Image * zero_pad_shifted_image(Image * a, int newx, int newy, int pad_mask);
static Image * zero_pad_unshifted_image(Image * a, int newx, int newy, int pad_mask);
static real dist_to_axis(int i, Image * in);
static real dist_to_center(int i, Image * in);
static real square_dist_to_center(int i, Image * in);
static real dist_to_corner(int i, Image * in);
static void random_rephase(Image *  img);
static Image * reflect_xy(Image * in, int in_place);
static Image * reflect_x(Image * in, int in_place);
static Image * reflect_y(Image * in, int in_place);
static void write_h5_img(Image * img,const char * filename, int output_precision);
static Image * read_imagefile(const char * filename);
static Image * read_tiff(const char * filename);
static  void write_tiff(Image * img,const char * filename);
static  void write_csv(Image * img,const char * filename);
static Image * read_png(const char * filename);
static int write_png(Image * img,const char * filename, int color);
static int write_vtk(Image * img,const char * filename);

real p_drand48(){
	real ret = ((double)rand())/(double)RAND_MAX;
  return (ret);
}

/* Returns the patterson from diffraction image a */
Image * sp_image_patterson(Image * a){
  Image * b = sp_image_duplicate(a,SP_COPY_DATA|SP_COPY_MASK);
  Image * c;
  Image * d;

  if(b->scaled){
    sp_cmatrix_mul_elements(b->image,b->image);
  }
  c = sp_image_fft(b);

  sp_image_free(b);
  d = sp_image_shift(c);
  sp_image_free(c);
  return d; 
}


Image * sp_image_reflect(Image * in, int in_place, int axis){
  if(axis == SP_AXIS_XY){
    return reflect_xy(in,in_place);
  }else if(axis == SP_AXIS_X){
    return reflect_x(in,in_place);
  }else if(axis == SP_AXIS_Y){
    return reflect_y(in,in_place);
  }
  return NULL;
}

/* reflect an image through the center, also 
 known as rotating by 180 degrees */
static Image * reflect_xy(Image * in, int in_place){
  return reflect_x(reflect_y(in,in_place),1);
}

/* reflect an image on the x axis, also 
 known as vertical flip */
static Image * reflect_x(Image * in, int in_place){
  Complex tmp;
  Image * out;
  int x,y,y2;
  if(in_place){
    out = in;
  }else{
    out = sp_image_duplicate(in,SP_COPY_DATA|SP_COPY_MASK);
  }
  for(x = 0;x<sp_cmatrix_cols(in->image);x++){
    for(y = 0;y<sp_cmatrix_rows(in->image)/2.0;y++){
      y2 = sp_cmatrix_rows(in->image)-y-1;
      tmp = sp_cmatrix_get(in->image,y,x);
      sp_cmatrix_set(out->image,y,x,sp_cmatrix_get(in->image,y2,x));
      sp_cmatrix_set(out->image,y2,x, tmp);
    }
  }
  return out;
}

/* reflect an image on the y axis, also 
 known as horizontal flip */
static Image * reflect_y(Image * in, int in_place){
  Complex tmp;
  Image * out;
  int x,y,x2;
  if(in_place){
    out = in;
  }else{
    out = sp_image_duplicate(in,SP_COPY_DATA|SP_COPY_DATA);
  }
  for(x = 0;x<sp_cmatrix_cols(in->image)/2.0;x++){
    x2 = sp_cmatrix_cols(in->image)-x-1;
    for(y = 0;y<sp_cmatrix_rows(in->image);y++){
      tmp = sp_cmatrix_get(in->image,y,x);
      sp_cmatrix_set(out->image,y,x,sp_cmatrix_get(in->image,y,x2));
      sp_cmatrix_set(out->image,y,x2,tmp);
    }
  }
  return out;
}

/* Image should be shifted and squared */
Image * sp_average_centrosymetry(Image * in){
  int x,y;
  real noise = 0;
  int ind1,ind2;
  int k = 0;
  Image * out = sp_image_duplicate(in,SP_COPY_DATA|SP_COPY_MASK);
  
  if(!in->shifted){
    fprintf(stderr,"Error: Using average_centrosymetry on unshifted image!\n");
    exit(1);
  }
  if(sp_cmatrix_cols(in->image) != sp_cmatrix_rows(in->image)){
    fprintf(stderr,"Error: Using average_centrosymetry with a non square image!\n");
    exit(1);
  }

  for(x = 0;x<sp_cmatrix_cols(in->image);x++){
    for(y = 0;y<sp_cmatrix_rows(in->image);y++){
      ind1 = sp_cmatrix_get_index(in->image,sp_cmatrix_rows(in->image)-1-y,sp_cmatrix_cols(in->image)-1-x);
      ind2= sp_cmatrix_get_index(in->image,y,x);
      if(dist_to_corner(x*sp_cmatrix_rows(in->image)+y,in) < 1500){
	if((in->image->data[ind1] + in->image->data[ind2]) && in->mask->data[ind1] && in->mask->data[ind2]){
	  noise += cabs(in->image->data[ind1] - in->image->data[ind2])/
 	  ((in->image->data[ind1] + in->image->data[ind2])/2);
	  k++;
	  out->image->data[ind2] = (in->image->data[ind1] + in->image->data[ind2])/2;
	}else{
	  out->image->data[ind2] = 0;
	}
      }
    }
  }
  fprintf(stderr,"Noise = %f\n",noise/k);
  fprintf(stderr,"Noise2 = %f\n",noise);
  return out;
}


Image * sp_make_shifted_image_square(Image * in){
  Image * out;
  int x,y;
  int xout = 0;
  int yout = 0;
  out = sp_image_duplicate(in,SP_COPY_DATA|SP_COPY_MASK);

  if(!in->shifted){
    fprintf(stderr,"Trying to call make_shifted_image_square with an unshifted image|\n");
    exit(1);
  }

  /* make it a square by removing a line or a column from the center */
  if(sp_cmatrix_cols(in->image) > sp_cmatrix_rows(in->image)){
    /* remove a column */
    sp_cmatrix_free(out->image);
    out->image = sp_cmatrix_alloc(sp_cmatrix_rows(in->image),sp_cmatrix_rows(in->image));
    sp_imatrix_free(out->mask);
    out->mask = sp_imatrix_alloc(sp_cmatrix_rows(in->image),sp_cmatrix_rows(in->image));
    yout = 0;
    xout = 0;
    for(x = 0;x<sp_cmatrix_cols(in->image);x++){
      for(y = 0;y<sp_cmatrix_rows(in->image);y++){
	if(fabs(x - (sp_cmatrix_cols(in->image)-1)/2.0) < (sp_cmatrix_cols(in->image)-sp_cmatrix_rows(out->image))/2.0){
	  continue;
	}
	sp_cmatrix_set(out->image,yout,xout,sp_cmatrix_get(in->image,y,x));
	sp_imatrix_set(out->mask,yout,xout,sp_imatrix_get(in->mask,y,x));
	yout = (yout+1)%sp_cmatrix_rows(out->image);
	if(yout == 0){
	  xout++;
	}
      }
    }
  }else if(sp_cmatrix_cols(in->image) < sp_cmatrix_rows(in->image)){
    /* remove a line */
    sp_cmatrix_free(out->image);
    out->image = sp_cmatrix_alloc(sp_cmatrix_cols(in->image),sp_cmatrix_cols(in->image));
    sp_imatrix_free(out->mask);
    out->mask = sp_imatrix_alloc(sp_cmatrix_cols(in->image),sp_cmatrix_cols(in->image));

    for(x = 0;x<sp_cmatrix_cols(in->image);x++){
      for(y = 0;y<sp_cmatrix_rows(in->image);y++){
	if(fabs(x - (sp_cmatrix_rows(in->image)-1)/2.0) < (sp_cmatrix_rows(in->image)-sp_cmatrix_cols(out->image))/2.0){
	  continue;
	}
	sp_cmatrix_set(out->image,yout,xout,sp_cmatrix_get(in->image,y,x));
	sp_imatrix_set(out->mask,yout,xout,sp_imatrix_get(in->mask,y,x));
	yout = (yout+1)%sp_cmatrix_rows(out->image);
	if(yout == 0){
	  xout++;
	}

      }
    }
  }
  return out;
}



/* Make an unshifted image square by padding with zeroes on the smaller dimension */
Image * sp_make_unshifted_image_square(Image * in){
  int size = sp_max(sp_cmatrix_cols(in->image),sp_cmatrix_rows(in->image));
  return zero_pad_image(in,size,size,0);
}



/*
  For unshifted images it shifted the quadrants around image_center and 
  extra zero pad is used on the smaller quadrants to make them the same size
  as the biggest quadrant.
  For shifted images it shifted the quadrants around (size[]-1)/2 

*/
Image * sp_image_shift(Image * img){
  Image * out;
  int i;
  int index1,index2;
  int x,y;
  int newx,newy;
  int max_x,max_y;

  /* fft shift the image */
  out = sp_image_duplicate(img,SP_COPY_DATA|SP_COPY_MASK);
  if(!img->shifted){
    max_x = sp_max(img->detector->image_center[0],sp_cmatrix_cols(img->image)-img->detector->image_center[0]);
    max_y = sp_max(img->detector->image_center[1],sp_cmatrix_rows(img->image)-img->detector->image_center[1]);
    sp_image_realloc(out,2*max_x,2*max_y);
  }

		   
  for(i = 0;i<sp_image_size(out);i++){
    out->image->data[i] = 0;
    out->mask->data[i] = 0;
  }
  if(img->shifted){
    out->detector->image_center[0] = (sp_cmatrix_cols(img->image))/2.0;
    out->detector->image_center[1] = (sp_cmatrix_rows(img->image))/2.0;
    img->detector->image_center[0] = (sp_cmatrix_cols(img->image))/2.0;
    img->detector->image_center[1] = (sp_cmatrix_cols(img->image))/2.0;
  }else{
    out->detector->image_center[0] = 0;
    out->detector->image_center[1] = 0;
  }
  out->shifted = !out->shifted;
  /* shift quadrants */
  for(x = 0;x<sp_cmatrix_cols(img->image);x++){
    for(y = 0;y<sp_cmatrix_rows(img->image);y++){
      index1 = x*sp_cmatrix_rows(img->image)+y;
      index2 = 0;
      if(y < img->detector->image_center[1]){
	if(x < img->detector->image_center[0]){
	  newx = sp_cmatrix_cols(out->image)-(img->detector->image_center[0]-x);
	  newy = sp_cmatrix_rows(out->image)-(img->detector->image_center[1]-y);
	  if(newx < sp_cmatrix_cols(img->image)/2.0 ||
	     newy < sp_cmatrix_rows(img->image)/2.0){
	    index2 = -1;
	  }
	}else{
	  newx = x-img->detector->image_center[0];
	  newy = sp_cmatrix_rows(out->image)-(img->detector->image_center[1]-y);
	  if(newx >= sp_cmatrix_cols(img->image)/2.0 ||
	     newy < sp_cmatrix_rows(img->image)/2.0){
	    index2 = -1;
	  }
	}
      }else{	
	if(x < img->detector->image_center[0]){
	  newx = sp_cmatrix_cols(out->image)-(img->detector->image_center[0]-x);
	  newy = y-img->detector->image_center[1];
	  if(newx < sp_cmatrix_cols(img->image)/2.0 ||
	     newy >= sp_cmatrix_rows(img->image)/2.0){
	    index2 = -1;
	  }
	}else{
	  newx = x-img->detector->image_center[0];
	  newy = y-img->detector->image_center[1];
	  if(newx >= sp_cmatrix_cols(img->image)/2.0 ||
	     newy >= sp_cmatrix_rows(img->image)/2.0){
	    index2 = -1;
	  }

	}
      }
      if(index2 != -1){
	index2 = sp_cmatrix_get_index(out->image,newy,newx);
      }

      if(index2 != -1){
	out->image->data[index2] = img->image->data[index1];
	out->mask->data[index2] = img->mask->data[index1];
      }
    }
  }

  
  return out;
}


/* resolution given in pixels from center. 
   Image should be shifted. Square window. */
Image * sp_image_low_pass(Image * img, int resolution){
  Image * res = sp_image_duplicate(img,SP_COPY_DATA|SP_COPY_MASK);
  int x,y,nx,ny;
  int dx,dy;
  if(img->shifted == 0){
    fprintf(stderr,"Error: Trying to limit resolution on an unshifted image\n");
  }
  if(resolution*2 > sp_cmatrix_cols(res->image) || resolution*2 > sp_cmatrix_rows(res->image)){
    return sp_image_duplicate(img,SP_COPY_DATA|SP_COPY_MASK);
  }  
  sp_cmatrix_free(res->image);
  sp_imatrix_free(res->mask);
  res->image = sp_cmatrix_alloc(resolution*2,resolution*2);
  res->mask = sp_imatrix_alloc(resolution*2,resolution*2);
  nx = 0;
  ny = 0;
  for(x = 0;x<sp_cmatrix_cols(img->image);x++){
    dx = x;
    if(sp_cmatrix_cols(img->image)-x-1 < dx){
      dx = sp_cmatrix_cols(img->image)-x-1;
    }
    for(y = 0;y<sp_cmatrix_rows(img->image);y++){
      dy = y;
      if(sp_cmatrix_rows(img->image)-y-1 < dy){
	dy = sp_cmatrix_rows(img->image)-y-1;
      }
      if(dx < resolution && dy < resolution){
	sp_cmatrix_set(res->image,ny,nx,sp_cmatrix_get(img->image,y,x));
	sp_imatrix_set(res->mask,ny,nx,sp_imatrix_get(img->mask,y,x));
	ny++;
	if(ny == sp_cmatrix_rows(res->image)){
	  nx++;
	  ny = 0;
	}
      }      
    }    
  }
  return res;
}



/*! This functoin multiplies it's input image with the mask after convoluting
 *   it with a selected filter
 *
 *  The mask should be set to 1 in the region of the image with meaningful image
 *  The possible values for the flag are:
 *
 *  SP_GAUSSIAN - *value corresponds to the standard deviation of the gaussian in pixels
 *
 */
void sp_image_smooth_edges(Image * img, sp_imatrix * mask, int flags, real * value){
  int i;
  Image * tmp = sp_image_duplicate(img,SP_COPY_DATA|SP_COPY_MASK);
  if(flags & SP_GAUSSIAN){
    for(i = 0;i<sp_image_size(tmp);i++){
      tmp->image->data[i] = mask->data[i];
    }
    Image * blur_mask = gaussian_blur(tmp,*value);

    /* eat out the edge of the mask*/
    for(i = 0;i<sp_image_size(blur_mask);i++){
      if(creal(blur_mask->image->data[i]) < 0.99){
	blur_mask->image->data[i] = 0;
      }
    }
      
    sp_image_free(tmp);
    tmp = blur_mask;
    blur_mask = gaussian_blur(tmp,*value);

    sp_cmatrix_mul_elements(img->image,blur_mask->image);
  }  
}


real sp_centro_sym_value(int index,Image * img){
  int x,y;
  real nx,ny;
  x = index/sp_cmatrix_rows(img->image);
  y = index%sp_cmatrix_rows(img->image);
  nx = 2*img->detector->image_center[0]-x;
  ny = 2*img->detector->image_center[1]-y;
  if(nx < 1 || nx >= sp_cmatrix_cols(img->image)-2 ||
     ny < 1 || ny >= sp_cmatrix_rows(img->image)-2){
    return -1;
  }
  return sp_image_interp(img, nx, ny);
}

int sp_centro_sym_index(int index,Image * img){
  int x,y;
  int nx,ny;
  x = index/sp_cmatrix_rows(img->image);
  y = index%sp_cmatrix_rows(img->image);
  nx = 2*img->detector->image_center[0]-x;
  ny = 2*img->detector->image_center[1]-y;
  if(nx < 0 || nx >= sp_cmatrix_cols(img->image) ||
     ny < 0 || ny >= sp_cmatrix_rows(img->image)){
    return index;
  }
  return nx*sp_cmatrix_rows(img->image)+ny;
}

Image * sp_centro_sym_correlation(Image  * img){
  int x,y;
  Image * res = sp_image_duplicate(img,SP_COPY_DATA|SP_COPY_MASK);
  int index = 0;
  real csvalue;
  for(x = 0;x<sp_cmatrix_cols(img->image);x++){
    for(y = 0;y<sp_cmatrix_rows(img->image);y++){
      csvalue = sp_centro_sym_value(index,img);
      if(!img->mask->data[index] || csvalue == -1 || cabs(img->image->data[index])+fabs(csvalue) < 1){
	res->image->data[index] = 1.0;
      }else{
	res->image->data[index] = 1.0 - cabs(img->image->data[index]-csvalue)/(cabs(img->image->data[index])+fabs(csvalue));      
      }
      if(creal(res->image->data[index]) < 0 || creal(res->image->data[index]) > 1){
	/* Houston we have a problem */
	exit(1);
      }
      index++;
    }
  }
  return res;
}


real sp_image_dist(Image * in, int i, int type){
  if(type == SP_TO_AXIS){
    return dist_to_axis(i,in);
  }if(type == SP_TO_CENTER){
    return dist_to_center(i,in);
  }if(type == SP_TO_CENTER2){
    return square_dist_to_center(i,in);
  }if(type == SP_TO_CORNER){
    return  dist_to_corner(i, in);
  }
  return -1;
}

static real dist_to_axis(int i, Image * in){
  int x = i/sp_cmatrix_rows(in->image);
  int y = i%sp_cmatrix_rows(in->image);
  real dx,dy;
  if(in->shifted){
    dx = MIN(x,sp_cmatrix_cols(in->image)-x);
    dy = MIN(y,sp_cmatrix_rows(in->image)-y);
  }else{
    dx = fabs(x-in->detector->image_center[0]);
    dy = fabs(y-in->detector->image_center[1]);
  }
  return MIN(dx,dy);
}



static real dist_to_center(int i, Image * in){
  int x = i/sp_cmatrix_rows(in->image);
  int y = i%sp_cmatrix_rows(in->image);
  real dx,dy;
  if(in->shifted){
    dx = MIN(x,sp_cmatrix_cols(in->image)-x);
    dy = MIN(y,sp_cmatrix_rows(in->image)-y);
  }else{
    dx = x-in->detector->image_center[0];
    dy = y-in->detector->image_center[1];
  }
  return sqrt(dx*dx+dy*dy);
}

static real square_dist_to_center(int i, Image * in){
  int x = i/sp_cmatrix_rows(in->image);
  int y = i%sp_cmatrix_rows(in->image);
  real dx,dy;
  if(in->shifted){
    dx = MIN(x,sp_cmatrix_cols(in->image)-x);
    dy = MIN(y,sp_cmatrix_rows(in->image)-y);
  }else{
    dx = fabs(x-in->detector->image_center[0]);
    dy = fabs(y-in->detector->image_center[1]);
  }
  return MAX(dx,dy);

}

static real dist_to_corner(int i, Image * in){
  int x = i/sp_cmatrix_rows(in->image);
  int y = i%sp_cmatrix_rows(in->image);
  real dx,dy;
  if(sp_cmatrix_cols(in->image)-1-x < x){
    dx = sp_cmatrix_cols(in->image)-1-x;
  }else{
    dx = x;
  }
  if(sp_cmatrix_rows(in->image)-1-y < y){
    dy = sp_cmatrix_cols(in->image)-1-y;
  }else{
    dy = y;
  }
  return sqrt(dx*dx+dy*dy);
}


void sp_image_dephase(Image *  img){
  int i = 0;
  img->phased = 0;    
  for(i = 0;i<sp_image_size(img);i++){
    img->image->data[i] = cabs(img->image->data[i]);
  }
}

void sp_image_rephase(Image *  img, int type){ 
  img->phased = 1;
  if(type == SP_ZERO_PHASE){
    return;
  }else if(type == SP_RANDOM_PHASE){
    random_rephase(img);
  }
}

static void random_rephase(Image *  img){  
  int i;
  real phase;
  for(i = 0;i<sp_image_size(img);i++){
    phase = p_drand48()*M_PI*2;
    img->image->data[i] = cos(phase)*img->image->data[i]+ sin(phase)*img->image->data[i]*I;
  }
  img->phased = 1;    
}


/* returns the phase of the image from [-PI, PI[ or -1 if there's an error */
Image * sp_image_get_phases(Image * img){
  Image * res;
  int i;
  if(!img->phased){
    return NULL;
  }
  res = sp_image_duplicate(img,SP_COPY_DETECTOR);
  sp_image_dephase(res);
  for(i = 0;i<sp_image_size(res);i++){
    res->image->data[i] = carg(img->image->data[i]);
  }
  return res;
}





void sp_image_add(Image * a, Image * b){
  sp_cmatrix_add(a->image,b->image,NULL);
}

void sp_image_sub(Image * a, Image * b){
  sp_cmatrix_sub(a->image,b->image);
}


/* mean m, standard deviation s */
/* gaussian distribution generator */
real sp_box_muller(real m, real s){        
  float x1, x2, w, y1;
  static real y2;
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


void sp_add_noise(Image * in, real level, int type){
  int i;
  if(type == SP_GAUSSIAN){
    if(level <= 0){
      return;
    }
    for(i = 0;i<sp_image_size(in);i++){
      in->image->data[i] = (1.0+sp_box_muller(0,level))*creal(in->image->data[i])+(1.0+sp_box_muller(0,level))*cimag(in->image->data[i])*I;
    }
  }
}

void sp_image_high_pass(Image * in, real radius){
  int x,y;
  real dist,dx,dy;
  if(radius <= 0){
    return;
  }
  for(x = 0;x<sp_cmatrix_cols(in->image);x++){
    for(y = 0;y<sp_cmatrix_rows(in->image);y++){
      if(x > sp_cmatrix_cols(in->image)/2.0){
	dx = sp_cmatrix_cols(in->image)-x;
      }else{
	dx = x;
      }
      if(y > sp_cmatrix_rows(in->image)/2.0){
	dy = sp_cmatrix_rows(in->image)-y;
      }else{
	dy = y;
      }
      dist = sqrt(dx*dx+dy*dy);      
      if(dist <= radius){
	in->mask->data[x*sp_cmatrix_rows(in->image)+y] = 0;
	in->image->data[x*sp_cmatrix_rows(in->image)+y] = 0;
      }
    }
  }
}

void sp_image_free(Image * in){
  free(in->detector);
  sp_cmatrix_free(in->image);
  sp_imatrix_free(in->mask);
  free(in);
}

Image * sp_image_duplicate(Image * in, int flags){
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
  res->image = sp_cmatrix_alloc(sp_cmatrix_rows(in->image),sp_cmatrix_cols(in->image));
  if(!res->image){
    perror("Out of memory!\n");
    abort();
  }
  if(flags & SP_COPY_DATA){
    sp_cmatrix_memcpy(res->image,in->image);
  }

  res->mask = sp_imatrix_alloc(sp_cmatrix_rows(in->image),sp_cmatrix_cols(in->image));
  if(!res->mask){
    perror("Out of memory!\n");
    abort();
  }
  if(flags & SP_COPY_MASK){
    sp_imatrix_memcpy(res->mask,in->mask);
  }
  return res;
}


Image * sp_image_alloc(int x, int y){
  Image  *res = malloc(sizeof(Image));
  if(!res){
    perror("Out of memory!\n");
    abort();
  }
  res->detector = malloc(sizeof(Detector));
  res->detector->image_center[0] = x/2;
  res->detector->image_center[1] = y/2;
  if(!res->detector){
    perror("Out of memory!\n");
    abort();
  }
  res->mask = sp_imatrix_alloc(y,x);
  if(!res->mask){
    perror("Out of memory!\n");
    abort();
  }  
  res->scaled = 0;
  res->image = sp_cmatrix_alloc(y,x);
  if(!res->image){
    perror("Out of memory!\n");
    abort();
  }  
  res->phased = 0;
  return res;
}



void sp_image_write(Image * img, const char * filename, int flags){
  char buffer[1024];
  strcpy(buffer,filename);
  for(int i = 0;i<strlen(buffer);i++){
    buffer[i] = tolower(buffer[i]);
  }
  /* select the correct function depending on the buffer extension */
  if(rindex(buffer,'.') && strcmp(rindex(buffer,'.'),".h5") == 0){
    /* we have an h5 file */
    write_h5_img(img,filename,sizeof(real));
  }else if(rindex(buffer,'.') && strcmp(rindex(buffer,'.'),".png") == 0){
    write_png(img,filename,flags);
  }else if(rindex(buffer,'.') && strcmp(rindex(buffer,'.'),".vtk") == 0){
    write_vtk(img,filename);
  }else if(rindex(buffer,'.') && (strcmp(rindex(buffer,'.'),".tif") == 0 ||strcmp(rindex(buffer,'.'),".tiff") == 0 )){
    write_tiff(img,filename);
  }else if(rindex(buffer,'.') && (strcmp(rindex(buffer,'.'),".csv") == 0)){
    write_csv(img,filename);
  }else{
    fprintf(stderr,"Unsupported file type: %s\n",filename);
  }
}

Image * sp_image_read(const char * filename, int flags){
  char buffer[1024];
  strcpy(buffer,filename);
  for(int i = 0;i<strlen(buffer);i++){
    buffer[i] = tolower(buffer[i]);
  }
  /* select the correct function depending on the filename extension */
  if(rindex(buffer,'.') && strcmp(rindex(buffer,'.'),".h5") == 0){
    /* we have an h5 file */
    return read_imagefile(filename);
  }else if(rindex(buffer,'.') && strcmp(rindex(buffer,'.'),".png") == 0){
    return read_png(filename);
  }else if(rindex(buffer,'.') && strcmp(rindex(buffer,'.'),".vtk") == 0){
    fprintf(stderr,"Cannot read VTK files!\n");
    return NULL;
  }else if(rindex(buffer,'.') && (strcmp(rindex(buffer,'.'),".tif") == 0 ||strcmp(rindex(buffer,'.'),".tiff") == 0 )){
    return read_tiff(filename);
  }else{
    fprintf(stderr,"Unsupported file type: %s\n",filename);
  }
  return NULL;
}


static void write_h5_img(Image * img,const char * filename, int output_precision){
  hid_t dataspace_id;
  hid_t dataset_id;
  hid_t file_id;
  int status;
  int version;
  hsize_t  dims[2];
  real values[2];
  sp_matrix * tmp;
  int i;
  hid_t out_type_id = 0;
  hid_t mem_type_id = 0;
  hid_t plist;
  hsize_t chunk_size[2] = {sp_cmatrix_cols(img->image),sp_cmatrix_rows(img->image)};
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

  dims[0] = sp_cmatrix_cols(img->image);
  dims[1] = sp_cmatrix_rows(img->image);
  file_id = H5Fcreate(filename,  H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  dataspace_id = H5Screate_simple( 2, dims, NULL );

  plist = H5Pcreate (H5P_DATASET_CREATE);
  H5Pset_chunk(plist,2,chunk_size);
  H5Pset_deflate(plist,6);

  dataset_id = H5Dcreate(file_id, "/mask", H5T_NATIVE_INT,
			 dataspace_id, plist);
  status = H5Dwrite(dataset_id,H5T_NATIVE_INT , H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, img->mask->data);
  status = H5Dclose(dataset_id);

  tmp = sp_matrix_alloc(sp_cmatrix_rows(img->image),sp_cmatrix_cols(img->image));
  for(i = 0;i<sp_image_size(img);i++){
    tmp->data[i] = creal(img->image->data[i]);
  }

  dataset_id = H5Dcreate(file_id, "/real", out_type_id,
			 dataspace_id, plist);
  status = H5Dwrite(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, tmp->data);
  status = H5Dclose(dataset_id);
  sp_matrix_free(tmp);

  if(img->phased){
    tmp = sp_matrix_alloc(sp_cmatrix_rows(img->image),sp_cmatrix_cols(img->image));
    for(i = 0;i<sp_image_size(img);i++){
      tmp->data[i] = cimag(img->image->data[i]);
    }

    dataset_id = H5Dcreate(file_id, "/imag",out_type_id,
			   dataspace_id, plist);
    status = H5Dwrite(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, tmp->data);
    status = H5Dclose(dataset_id);
    sp_matrix_free(tmp);

  }
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


  version = 2;
  dataset_id = H5Dcreate(file_id, "/version", H5T_NATIVE_INT,
			 dataspace_id, H5P_DEFAULT);
  status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, &version);
  status = H5Dclose(dataset_id);


  status = H5Sclose(dataspace_id);


  status = H5Fclose(file_id);
}

Image * read_imagefile(const char * filename){
  Image * res = malloc(sizeof(Image));
  int file_id,dataset_id,space;
  int status,i;
  int version;
  hsize_t dims[2];
  hid_t mem_type_id = 0;
  H5E_auto_t func;
  void * client_data;
  real values[2];
  sp_matrix * tmp;
  if(sizeof(real) == sizeof(float)){
    mem_type_id = H5T_NATIVE_FLOAT;
  }else if(sizeof(real) == sizeof(double)){
    mem_type_id = H5T_NATIVE_DOUBLE;
  }else{
    abort();
  }



  res->detector = malloc(sizeof(Detector));

  
  file_id = H5Fopen(filename,H5F_ACC_RDONLY,H5P_DEFAULT);

  H5Eget_auto(&func,&client_data);
  /* turn off warning to check version because it might not exist */
  H5Eset_auto(NULL,NULL);
  dataset_id = H5Dopen(file_id, "/version");
  H5Eset_auto(func,client_data);
  if(dataset_id>=0){
    status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
		     H5P_DEFAULT, &version);
    if(version == 2){
      dataset_id = H5Dopen(file_id, "/mask");
      space = H5Dget_space(dataset_id);
      H5Sget_simple_extent_dims(space,dims,NULL);
      res->image = sp_cmatrix_alloc(dims[1],dims[0]);
      
      res->mask = sp_imatrix_alloc(dims[1],dims[0]);
      status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
		       H5P_DEFAULT, res->mask->data);
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
	tmp = sp_matrix_alloc(sp_imatrix_rows(res->mask),sp_imatrix_cols(res->mask));
	dataset_id = H5Dopen(file_id, "/imag");
	status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
			 H5P_DEFAULT, tmp->data);
	status = H5Dclose(dataset_id);
	for(i = 0;i<sp_matrix_size(tmp);i++){
	  res->image->data[i] = tmp->data[i]*I;
	}
	sp_matrix_free(tmp);
      }
      
      tmp = sp_matrix_alloc(sp_imatrix_rows(res->mask),sp_imatrix_cols(res->mask));
      dataset_id = H5Dopen(file_id, "/real");
      status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		       H5P_DEFAULT, tmp->data);
      status = H5Dclose(dataset_id);
      for(i = 0;i<sp_matrix_size(tmp);i++){
	res->image->data[i] += tmp->data[i];
      }
      sp_matrix_free(tmp);
      
      status = H5Fclose(file_id);
    }
  }else{
    dataset_id = H5Dopen(file_id, "/mask");
    space = H5Dget_space(dataset_id);
    H5Sget_simple_extent_dims(space,dims,NULL);
    res->image = sp_cmatrix_alloc(dims[1],dims[0]);
    
    res->mask = sp_imatrix_alloc(dims[1],dims[0]);
    tmp = sp_matrix_alloc(sp_imatrix_rows(res->mask),sp_imatrix_cols(res->mask));

    status = H5Dread(dataset_id,mem_type_id , H5S_ALL, H5S_ALL,
		     H5P_DEFAULT, tmp->data);
    for(i = 0;i<sp_matrix_size(tmp);i++){
      res->mask->data[i] = tmp->data[i];
    }
    sp_matrix_free(tmp);
    
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
      tmp = sp_matrix_alloc(sp_imatrix_rows(res->mask),sp_imatrix_cols(res->mask));
      dataset_id = H5Dopen(file_id, "/complex");
      status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		       H5P_DEFAULT, tmp->data);
      status = H5Dclose(dataset_id);
      for(i = 0;i<sp_matrix_size(tmp);i++){
	res->image->data[i] = tmp->data[i]*I;
      }
      sp_matrix_free(tmp);
      tmp = sp_matrix_alloc(sp_imatrix_rows(res->mask),sp_imatrix_cols(res->mask));
      dataset_id = H5Dopen(file_id, "/real");
      status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		       H5P_DEFAULT, tmp->data);
      status = H5Dclose(dataset_id);
      for(i = 0;i<sp_matrix_size(tmp);i++){
	res->image->data[i] += tmp->data[i];
      }
      sp_matrix_free(tmp);
      
    }else{
      if(!res->scaled){
	tmp = sp_matrix_alloc(sp_imatrix_rows(res->mask),sp_imatrix_cols(res->mask));
	dataset_id = H5Dopen(file_id, "/intensities");
	status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
			 H5P_DEFAULT, tmp->data);
	status = H5Dclose(dataset_id);
	for(i = 0;i<sp_matrix_size(tmp);i++){
	  res->image->data[i] += tmp->data[i];
	}
	sp_matrix_free(tmp);
      }else{
	tmp = sp_matrix_alloc(sp_imatrix_rows(res->mask),sp_imatrix_cols(res->mask));
	dataset_id = H5Dopen(file_id, "/amplitudes");
	status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
			 H5P_DEFAULT, tmp->data);
	status = H5Dclose(dataset_id);
	for(i = 0;i<sp_matrix_size(tmp);i++){
	  res->image->data[i] += tmp->data[i];
	}
	sp_matrix_free(tmp);	 
      }
    }        
    status = H5Fclose(file_id);    
  }
  return res;
  
}



Image * read_tiff(const char * filename){
  Image * out = malloc(sizeof(Image));
  out->detector = malloc(sizeof(Detector));
  int bpp = 4;  
  int datatype = 0;
  int width,height;
  unsigned char * tmpuc;
  int nstrips;
  int stripsize;
  int i,x,y;
  unsigned char * img;
  float * tmpf;
  unsigned short * tmpi;
  unsigned short * tmpui;
  TIFF * tif; 

  tif = TIFFOpen(filename, "r");
  if(TIFFGetField(tif,TIFFTAG_BITSPERSAMPLE,&bpp)){
    bpp /= 8;
  }
  if(!TIFFGetField(tif,TIFFTAG_SAMPLEFORMAT,&datatype)){
    if(bpp == 1){
      datatype = SAMPLEFORMAT_VOID;
    }else if(bpp == 2){
      datatype = SAMPLEFORMAT_UINT;
    }
  }
  
  if(!TIFFGetField(tif,TIFFTAG_IMAGELENGTH,&height)){
    perror("Could not get image height!\n");
    return NULL;
  }
  if(!TIFFGetField(tif,TIFFTAG_IMAGEWIDTH,&width)){
    perror("Could not get image width!\n");
    return NULL;
  }
  
  nstrips = TIFFNumberOfStrips(tif);
  stripsize = TIFFStripSize(tif);
  img = malloc(nstrips*stripsize);
  for(i = 0;i<nstrips;i++){
    TIFFReadEncodedStrip(tif,i,img+i*stripsize,stripsize);
  }
  TIFFClose(tif);
  
  /* Transpose image, because TIFF is saved row by row (which they call strips)
     unlike Hawk that saves column by column */
  out->image = sp_cmatrix_alloc(height,width);
  out->mask = sp_imatrix_alloc(height,width);
  if(datatype == SAMPLEFORMAT_UINT){
    tmpui = (unsigned short *)img;
    for(x = 0;x<width;x++){
      for(y = 0;y<height;y++){
	tmpui = (unsigned short *)(img+y*width*bpp+x*bpp);
	out->image->data[x*height+y] = *tmpui;
      }    
    }
  }else if(datatype == SAMPLEFORMAT_IEEEFP){
    for(x = 0;x<width;x++){
      for(y = 0;y<height;y++){
	tmpf = (float *)(img+y*width*bpp+x*bpp);
	out->image->data[x*height+y] = *tmpf;
      }    
    }
  }else if(datatype == SAMPLEFORMAT_VOID){
    for(x = 0;x<width;x++){
      for(y = 0;y<height;y++){
	tmpuc = img+y*width*bpp+x*bpp;
	out->image->data[x*height+y] = *tmpuc;
      }    
    }
  }else if(datatype == SAMPLEFORMAT_INT){
    for(x = 0;x<width;x++){
      for(y = 0;y<height;y++){
	tmpi = (unsigned short *)(img+y*width*bpp+x*bpp);
	out->image->data[x*height+y] = (*tmpi);
      }    
    }
  }
  for(i = 0;i<sp_cmatrix_size(out->image);i++){
    out->mask->data[i] = 1;
  }
  free(img);
  out->scaled = 0;
  out->phased = 0;
  out->shifted = 0;
  out->detector->image_center[0] = width/2;
  out->detector->image_center[1] = height/2;
  return out;
}


void write_tiff(Image * img,const char * filename){
  float * data;
  int nstrips;
  int stripsize;
  TIFF * tif;
  int x,y;
  int width = sp_image_width(img);
  int height = sp_image_height(img);


  tif = TIFFOpen(filename, "w");  
  nstrips = height;
  stripsize = width*sizeof(float);

  TIFFSetField(tif,TIFFTAG_IMAGEWIDTH,width);
  TIFFSetField(tif,TIFFTAG_ROWSPERSTRIP,1);
  TIFFSetField(tif,TIFFTAG_IMAGELENGTH,height);
  TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 32);
  TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);
  TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
  TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
  TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
  data = malloc(nstrips*stripsize);
  for(y = 0;y<sp_image_height(img);y++){
    for(x = 0;x<sp_image_width(img);x++){
      data[x] =cabsr(sp_image_get(img,x,y));      
    }
    TIFFWriteEncodedStrip(tif,y,data,stripsize);
  }

  /*  for(i = 0;i<nstrips;i++){
    TIFFWriteEncodedStrip(tif,i,&(data[i*stripsize/4]),stripsize);
    }*/
  TIFFClose(tif);
  free(data);
}


/*! Write an image to CSV format 
 */
void write_csv(Image * img,const char * filename){
  FILE * f = fopen(filename,"w");
  int x,y;
  if(!f){
    perror("Could not write CSV");
    exit(0);
  }
  fprintf(f,"x,y,amplitude,phase,real,imaginary\n");
  for(y = 0;y<sp_image_height(img);y++){
    for(x = 0;x<sp_image_width(img);x++){
      fprintf(f,"%d,%d,%f,%f,%f,%f\n",x,y,cabsr(sp_image_get(img,x,y)),carg(sp_image_get(img,x,y)),creal(sp_image_get(img,x,y)),cimag(sp_image_get(img,x,y)));
    }
  }
  fclose(f);
}

/* b must fit inside a or the other way around */

/* If wrap arond is on, a is taken to be periodic and the outout size will
   be as big as the input 
   
   Otherwise the output size is equal to a+b-1

   If size is not a NULL pointer it will be the size of transform used.
   A size smaller than the size of a+b will results in wrap around effects.
   If size is NULL the size of a will be used (resulting in wrap around effects).

*/
Image * sp_image_cross_correlate(Image * a, Image * b, int * size){
  int x;
  int y;  
  int i;
  Image * res;
  Image * a_ft;
  Image * b_ft;
  Image * tmp;

  if(sp_cmatrix_cols(a->image) < sp_cmatrix_cols(b->image)){
    /* swap */
    res = a;
    a = b;
    b = res;
  }

  if(size){
    x = size[0];
    y = size[1];
  }else{
    x = sp_cmatrix_cols(a->image);
    y = sp_cmatrix_rows(a->image);
  }

  tmp = zero_pad_image(a,x,y,1);
  a_ft = sp_image_fft(tmp);
  sp_image_free(tmp);

  tmp = zero_pad_image(b,x,y,1);
  b_ft = sp_image_fft(tmp);
  sp_image_free(tmp);

  tmp = sp_image_duplicate(a_ft,SP_COPY_DETECTOR);
  tmp->shifted = 1;
  sp_image_rephase(tmp,SP_ZERO_PHASE);
  /* Now do the multiplication in fourier space */
  /* Using the Convolution Theorem */
  for(i = 0;i<x*y;i++){
    tmp->image->data[i] = a_ft->image->data[i]*conj(b_ft->image->data[i]);
  }

  sp_image_free(a_ft);
  sp_image_free(b_ft);
  /* Backtransform */
  res = sp_image_ifft(tmp);
  sp_image_free(tmp);

  /* should be all real */
  for(i = 0;i<sp_image_size(res);i++){
    res->image->data[i] /= sp_image_size(res);
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
Image * sp_image_convolute(Image * a, Image * b, int * size){
  int x;
  int y;  
  int i;
  Image * res;
  Image * a_ft;
  Image * b_ft;
  Image * tmp;

  if(sp_cmatrix_cols(a->image) < sp_cmatrix_cols(b->image)){
    /* swap */
    res = a;
    a = b;
    b = res;
  }
  
  if(!size){
    x = sp_cmatrix_cols(a->image)/*+(sp_cmatrix_cols(b->image)-1)/2*/;
    y = sp_cmatrix_rows(a->image)/*+(sp_cmatrix_rows(b->image)-1)/2*/;
  }else{
    x = size[0];
    y = size[1];
  }

  tmp = zero_pad_image(a,x,y,1);
/*  sp_image_dephase(tmp);*/
  a_ft = sp_image_fft(tmp);
  sp_image_free(tmp);

  tmp = zero_pad_image(b,x,y,1);
/*  sp_image_dephase(tmp);*/
  b_ft = sp_image_fft(tmp);
  sp_image_free(tmp);

  tmp = sp_image_duplicate(a_ft,SP_COPY_DETECTOR);
  tmp->shifted = 1;
/*  sp_image_rephase(tmp,SP_ZERO_PHASE);*/
  /* Now do the multiplication in fourier space */
  /* Using the Convolution Theorem */
  for(i = 0;i<x*y;i++){
    tmp->image->data[i] = a_ft->image->data[i]*b_ft->image->data[i];
  }

  sp_image_free(a_ft);
  sp_image_free(b_ft);
  /* Backtransform */
  res = sp_image_ifft(tmp);
  sp_image_free(tmp);

/*  sp_image_dephase(res);*/
  /* should be all real */
  for(i = 0;i<sp_image_size(res);i++){
    res->image->data[i] /= sp_image_size(res);
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
  real total_filter = 0;
  Image * filter_img = sp_image_alloc(filter_side,filter_side);
  Image * centered_filter;
  Image * res;
  Image * tmp;
  filter_img->detector->image_center[0] = (filter_side-1)/2.0;
  filter_img->detector->image_center[1] = (filter_side-1)/2.0;
  
  sp_image_dephase(filter_img);
  for(x = -ceil(radius)*3;x<=ceil(radius)*3;x++){
    i = x+ceil(radius)*3;
    for(y = -ceil(radius)*3;y<=ceil(radius)*3;y++){
      j = y+ceil(radius)*3;
      filter_img->image->data[i*filter_side+j] = 1/sqrt(2*M_PI*radius) * exp(-(x*x+y*y)/(2*radius*radius));
      /* Make the filter symmetric in the imaginary part */
/*      filter_img->image->data[i*filter_side+j] = filter_img->image->data[i*filter_side+j] + filter_img->image->data[i*filter_side+j]*I;*/
      total_filter += filter_img->image->data[i*filter_side+j];
    }
  }
  for(i = 0;i<sp_image_size(filter_img);i++){
    filter_img->image->data[i] /= total_filter;
  }
  centered_filter = shift_center_to_top_left(filter_img);
  centered_filter->shifted = 1;
  sp_image_free(filter_img);
  res = sp_image_convolute(in, centered_filter,NULL);
  sp_image_free(centered_filter);
  /* we should crop the result if it's bigger than the input */
  if(sp_image_size(res) > sp_image_size(in)){
    tmp = rectangle_crop(res, (sp_cmatrix_cols(res->image)-sp_cmatrix_cols(in->image))/2,
			 (sp_cmatrix_rows(res->image)-sp_cmatrix_rows(in->image))/2, 
			 sp_cmatrix_cols(in->image)/2-1+(sp_cmatrix_cols(res->image)-sp_cmatrix_cols(in->image))/2,
			 sp_cmatrix_cols(in->image)/2-1+(sp_cmatrix_rows(res->image)-sp_cmatrix_rows(in->image))/2);
    sp_image_free(res);
    res = tmp;
  }
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
  sp_cmatrix * filter = sp_cmatrix_alloc(filter_side,filter_side);
  real total_filter = 0;
  Image * filter_img = sp_image_duplicate(in,SP_COPY_DETECTOR);
  Image * centered_filter;
  Image * res;
  filter_img->detector->image_center[0] = (filter_side-1)/2.0;
  filter_img->detector->image_center[1] = (filter_side-1)/2.0;
  sp_cmatrix_free(filter_img->image);
  filter_img->image = filter;
  sp_image_dephase(filter_img);
  for(x = -ceil(radius)*3;x<=ceil(radius)*3;x++){
    i = x+ceil(radius)*3;
    for(y = -ceil(radius)*3;y<=ceil(radius)*3;y++){
      j = y+ceil(radius)*3;
      filter->data[i*filter_side+j] = 1.0/((2*radius+1)*(2*radius+1));
      total_filter += filter->data[i*filter_side+j];
    }
  }
  for(i = 0;i<sp_image_size(filter_img);i++){
    filter_img->image->data[i] /= total_filter;
  }
  centered_filter = shift_center_to_top_left(filter_img);
  sp_image_free(filter_img);
  res = sp_image_convolute(in, centered_filter,NULL);
  sp_image_free(centered_filter);
  return res;
}


int write_mask_to_png(Image * img, char * filename, int color){
  Image  * res = sp_image_duplicate(img,SP_COPY_DATA|SP_COPY_MASK);
  int ret;
  int i;
  for(i = 0;i<sp_image_size(img);i++){
    res->image->data[i] = res->mask->data[i];
  }
  ret = write_png(res,filename,color);
  sp_image_free(res);
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

Image * read_png(const char * filename){
 FILE *fp = fopen(filename, "rb");
 int i,j;
 png_uint_32 width,height;
 int bit_depth,color_type,interlace_type,compression_type,filter_method;
 png_structp png_ptr = png_create_read_struct
   (PNG_LIBPNG_VER_STRING, (png_voidp)NULL,
    NULL,NULL);
 png_infop info_ptr = png_create_info_struct(png_ptr);
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
 res = sp_image_alloc(width,height);
 for(i = 0;i<height;i++){
   for(j = 0;j<width;j++){
     res->image->data[j*height+i] = row_pointers[i][(int)(j*bit_depth/8)];
   }
 }
 return res;
}




int write_png(Image * img,const char * filename, int color){

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
  real phase;
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
  png_set_IHDR(png_ptr, info_ptr, sp_cmatrix_cols(img->image), sp_cmatrix_rows(img->image),
	       bit_depth, color_type, interlace_type,
	       compression_type, filter_method);
  
  

  row_pointers = png_malloc(png_ptr,sp_cmatrix_rows(img->image)*sizeof(png_byte *));
  for (i=0; i<sp_cmatrix_rows(img->image); i++){
    row_pointers[i] = png_malloc(png_ptr,sp_cmatrix_cols(img->image)*pixel_size*sizeof(png_byte));
  }
  
  /* We're gonna scale the image so that it fits on the 8 bits */
  min_v = sp_cmatrix_min(img->image,NULL);
  max_v = sp_cmatrix_max(img->image,NULL);
  if(max_v-min_v){
    scale = 1/(max_v-min_v);
  }else{
    scale = 1;
  }
  offset = min_v;
  i = 0;
  log_of_2 = log(2.0);
  /* this is a special kind of color */
  for(x = 0;x<sp_cmatrix_cols(img->image);x++){
    for(y = 0;y<sp_cmatrix_rows(img->image);y++){
      /* traditional color scale taken from gnuplot manual */
      if(color & LOG_SCALE){
	value = log((cabs(img->image->data[i])-offset)*scale+1)/log_of_2;
      }else{
	value = ((cabs(img->image->data[i])-offset)*scale);
      }
      if(color & COLOR_PHASE){
	phase = (256*(cargr(img->image->data[i])+3.1416)/(2*3.1416));
	row_pointers[y][x*3] =  sqrt(value)*color_table[0][(int)phase];
	row_pointers[y][x*3+1] = sqrt(value)*color_table[1][(int)phase];
	row_pointers[y][x*3+2] = sqrt(value)*color_table[2][(int)phase];
      }else{
	value *= 255;
	row_pointers[y][x*3] =  color_table[0][(int)value];
	row_pointers[y][x*3+1] = color_table[1][(int)value];
	row_pointers[y][x*3+2] = color_table[2][(int)value];
      }
      i++;
    }
  }
  png_set_rows(png_ptr, info_ptr, row_pointers);
  
  png_write_png(png_ptr, info_ptr, png_transforms, NULL);
  png_write_flush(png_ptr);
  /* png_write_end(png_ptr, info_ptr);*/
  for(i=0; i<sp_cmatrix_rows(img->image); i++){
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
  for(i = 0;i<sp_image_size(fobs);i++){
    if(!fobs->mask->data[i] || creal(fobs->image->data[i]) < low_intensity_cutoff){
      continue;
    }
    num +=  cabs(fobs->image->data[i]-fcalc->image->data[i]);
    den += fobs->image->data[i];
  }
  return num/den;
}



/* Low pass filter using a centered square window of side edge_size */
Image * low_pass_square_filter(Image * in, int edge_size){
  Image * fft_img = sp_image_fft(in);
  Image * res;
  Image * tmp;
  int i = 0;
  for(i = 0;i<sp_image_size(in);i++){
    if(square_dist_to_center(i,in) > edge_size/2.0){
      fft_img->image->data[i] = 0;
    }
  }
  tmp = sp_image_duplicate(fft_img,SP_COPY_DATA|SP_COPY_MASK);
  for(i = 0;i<sp_image_size(tmp);i++){
    tmp->image->data[i] = log(tmp->image->data[i]+1);
  }
  write_png(tmp,"low_pass.png",COLOR_JET);
  sp_image_free(tmp);

  res = sp_image_ifft(fft_img);
  /* scale appropriately */
  for(i = 0;i<sp_image_size(res);i++){
    res->image->data[i] /= sp_image_size(res);
  }
  sp_image_free(fft_img);
  fft_img = sp_image_fft(res);
  tmp = sp_image_duplicate(fft_img,SP_COPY_DATA|SP_COPY_MASK);
  for(i = 0;i<sp_image_size(tmp);i++){
    tmp->image->data[i] = log(tmp->image->data[i]+1);
  }
  write_png(tmp,"after_low_pass.png",COLOR_JET);
  sp_image_free(tmp);

  
  return res;
}

/* Low pass filter using a centered gaussian window of side edge_size */
Image * low_pass_gaussian_filter(Image * in, int edge_size){
  Image * fft_img = sp_image_fft(in);
  Image * res;
  Image * mask;
  int i = 0;
  gaussian_filter(fft_img,edge_size/2.0,1);
  res = sp_image_ifft(fft_img);
  sp_image_free(fft_img);
  /* scale appropriately */
  for(i = 0;i<sp_image_size(res);i++){
    res->image->data[i] /= sp_image_size(res);
  }
  if(!in->phased){
    sp_image_dephase(res);
  }
  /* Also low pass filter the mask */
  mask = sp_image_duplicate(in,SP_COPY_DATA|SP_COPY_MASK);
  memcpy(mask->image,in->mask,sp_image_size(in)*sizeof(real)); 
  fft_img = sp_image_fft(mask);
  gaussian_filter(fft_img,edge_size/2.0,1);
  sp_image_free(mask);
  mask = sp_image_ifft(fft_img);
  sp_image_free(fft_img);
  /* scale appropriately */
  for(i = 0;i<sp_image_size(mask);i++){
    /* if the mask is not really want then we have unkown information and we'll make it 0 */
    if(cabs(mask->image->data[i]) < sp_image_size(res)-1){
      res->mask->data[i] = 0;
    }else{
      res->mask->data[i] = 1;
    }
  }
  sp_image_free(mask);
  
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
    res = sp_image_duplicate(in,SP_COPY_DATA|SP_COPY_MASK);
  }
  /* the scaling of 7 makes sure the pattern is scaled by 0.0009 at the edge */

  for(i = 0;i<sp_image_size(in);i++){
    scaling = (dist_to_center(i,in)/(radius))*(dist_to_center(i,in)/(radius));
    res->image->data[i] *= exp(-scaling*scaling_factor);
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
Image * zero_pad_shifted_image(Image * a, int newx, int newy,int pad_mask){
  Image * out;
  int x,y;
  int sourcex;
  int sourcey;
  if(newx < sp_cmatrix_cols(a->image) || 
     newy < sp_cmatrix_rows(a->image)){
    fprintf(stderr,"Negative padding!\n");
    abort();
  }else if(newx == sp_cmatrix_cols(a->image) &&
	   newy == sp_cmatrix_rows(a->image)){
    return sp_image_duplicate(a,SP_COPY_DATA|SP_COPY_MASK);
  }
  out = sp_image_duplicate(a,SP_COPY_DETECTOR);
  sp_cmatrix_free(out->image);
  out->image = sp_cmatrix_alloc(newy,newx);
  sp_imatrix_free(out->mask);
  out->mask = sp_imatrix_alloc(newy,newx);

  for(x = 0;x<sp_cmatrix_cols(out->image);x++){
    if(x < sp_cmatrix_cols(a->image)/2.0){
      sourcex = x;
    }else if(sp_cmatrix_cols(out->image)-x-1 < (sp_cmatrix_cols(a->image)-1)/2.0){
      sourcex = (sp_cmatrix_cols(a->image)-1)-(sp_cmatrix_cols(out->image)-x-1);
    }else{
      sourcex = -1;
    }

    for(y = 0;y<sp_cmatrix_rows(out->image);y++){
      if(y < sp_cmatrix_rows(a->image)/2.0){
	sourcey = y;
      }else if(sp_cmatrix_rows(out->image)-y-1 < (sp_cmatrix_rows(a->image)-1)/2.0){
	sourcey = (sp_cmatrix_cols(a->image)-1)-(sp_cmatrix_rows(out->image)-y-1);
      }else{
	sourcey = -1;
      }
      if(sourcey == -1 || sourcex == -1){
	out->image->data[x*sp_cmatrix_rows(out->image)+y] = 0;
	out->mask->data[x*sp_cmatrix_rows(out->image)+y] = pad_mask;
      }else{
	out->image->data[x*sp_cmatrix_rows(out->image)+y] = a->image->data[sourcex*sp_cmatrix_rows(a->image)+sourcey];
	out->mask->data[x*sp_cmatrix_rows(out->image)+y] = a->mask->data[sourcex*sp_cmatrix_rows(a->image)+sourcey];
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
Image * zero_pad_unshifted_image(Image * a, int newx, int newy, int pad_mask){
  Image * out;
  int x,y;
  int sourcex;
  int sourcey;
  if(newx < sp_cmatrix_cols(a->image) || 
     newy < sp_cmatrix_rows(a->image)){
    fprintf(stderr,"Negative padding!\n");
    abort();
  }else if(newx == sp_cmatrix_cols(a->image) &&
	   newy == sp_cmatrix_rows(a->image)){
    return sp_image_duplicate(a,SP_COPY_DATA|SP_COPY_MASK);
  }
  out = sp_image_duplicate(a,SP_COPY_DETECTOR);
  sp_cmatrix_free(out->image);
  out->image = sp_cmatrix_alloc(newy,newx);
  sp_imatrix_free(out->mask);
  out->mask = sp_imatrix_alloc(newy,newx);

  for(x = 0;x<sp_cmatrix_cols(out->image);x++){
    if(x < sp_cmatrix_cols(a->image)){
      sourcex = x;
    }else{
      sourcex = -1;
    }
    
    for(y = 0;y<sp_cmatrix_rows(out->image);y++){
      if(y < sp_cmatrix_rows(a->image)){
	sourcey = y;
      }else{
	sourcey = -1;
      }
      if(sourcey == -1 || sourcex == -1){
	out->image->data[x*sp_cmatrix_rows(out->image)+y] = 0;
	out->mask->data[x*sp_cmatrix_rows(out->image)+y] = pad_mask;
      }else{
	out->image->data[x*sp_cmatrix_rows(out->image)+y] = a->image->data[sourcex*sp_cmatrix_rows(a->image)+sourcey];
	out->mask->data[x*sp_cmatrix_rows(out->image)+y] = a->mask->data[sourcex*sp_cmatrix_rows(a->image)+sourcey];
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
  Image * out = sp_image_duplicate(a,SP_COPY_DATA|SP_COPY_MASK);
  int x,y;
  int destx;
  int desty;
  for(x = 0;x<sp_cmatrix_cols(out->image);x++){
    destx = (int)(x+sp_cmatrix_cols(out->image)-out->detector->image_center[0])%sp_cmatrix_cols(out->image);
    for(y = 0;y<sp_cmatrix_rows(out->image);y++){
      desty = (int)(y+sp_cmatrix_rows(out->image)-out->detector->image_center[1])%sp_cmatrix_rows(out->image);
      out->image->data[destx*sp_cmatrix_rows(out->image)+desty] = a->image->data[x*sp_cmatrix_rows(a->image)+y];
    }
  }
  return out;
}


/* Csiszar I-divergence */
real I_divergenge(Image * a, Image * b){
  double suma = 0;
  double sumb = 0;
  int i;
  for(i = 0;i<sp_image_size(a);i++){
    suma += a->image->data[i] - b->image->data[i];
    if(b->image->data[i]){
      sumb += b->image->data[i] * log(b->image->data[i]/(a->image->data[i]+FLT_EPSILON));
    }
  }
  return suma+sumb;
}


real integrated_intensity(Image * a){
  double sum = 0;
  int i;
  for(i = 0;i<sp_image_size(a);i++){
    sum += a->image->data[i];
  }
  return sum;
}

int write_vtk(Image * img, const char * filename){
  FILE * f = fopen(filename,"w");
  int x,y;
  if(!f){
    perror("Bad file in write_vtk!");
    abort();
  }
  fprintf(f,"# vtk DataFile Version 2.0\n");
  fprintf(f,"Generated by image_util write_vtk()\n");
  fprintf(f,"ASCII\n");
  fprintf(f,"DATASET STRUCTURED_POINTS\n");
  fprintf(f,"DIMENSIONS %d %d 1\n",sp_cmatrix_cols(img->image),sp_cmatrix_rows(img->image));
  fprintf(f,"ORIGIN 0 %d 0\n",sp_cmatrix_rows(img->image));
  fprintf(f,"SPACING 1 -1 1\n");
  fprintf(f,"POINT_DATA %d\n",sp_image_size(img));
  fprintf(f,"SCALARS amplitudes float 1\n");
  fprintf(f,"LOOKUP_TABLE default\n");
  fprintf(f,"%6g",cabs(img->image->data[0]));
  for(y = 0;y<sp_cmatrix_rows(img->image);y++){
    for(x = 0;x<sp_cmatrix_cols(img->image);x++){
      fprintf(f," %g",cabs(img->image->data[x*sp_cmatrix_rows(img->image)+y]));    
    }
  }
/*  for(i = 1;i<sp_image_size(img);i++){
    fprintf(f," %g",img->image->data[i]);    
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
  cropped = sp_image_duplicate(in,SP_COPY_DETECTOR);
  cropped->detector->image_center[0] -= x1;
  cropped->detector->image_center[1] -= y1;
  sp_cmatrix_free(cropped->image);
  cropped->image = sp_cmatrix_alloc(y2-y1+1,x2-x1+1);
  sp_imatrix_free(cropped->mask);
  cropped->mask = sp_imatrix_alloc(y2-y1+1,x2-x1+1);


  for(i = x1;i<= x2;i++){
    memcpy(&cropped->image->data[(i-x1)*sp_cmatrix_rows(cropped->image)],&in->image->data[(i)*sp_cmatrix_rows(in->image)+y1],sp_cmatrix_rows(cropped->image)*sizeof(real));
    memcpy(&cropped->mask->data[(i-x1)*sp_cmatrix_rows(cropped->image)],&in->mask->data[(i)*sp_cmatrix_rows(in->image)+y1],sp_cmatrix_rows(cropped->image)*sizeof(real));
  }
  return cropped;
}


void find_center(Image * img, real * center_x, real * center_y){
  int x,y;
  float bx = -1;
  float by = -1;
  Image * a = sp_image_convolute(img,img,NULL);
  real max = 0;
  int index;
  max = sp_cmatrix_max(a->image,&index);
  sp_cmatrix_get_row_col(a->image,index,&y,&x);
  bx = x;
  by = y;
  if(bx < sp_cmatrix_cols(img->image)/2.0){
    bx = (sp_cmatrix_cols(img->image))/2.0+bx/2.0;
  }else{
    bx = (sp_cmatrix_cols(img->image))/2.0-(sp_cmatrix_cols(img->image)-bx)/2.0;
  }
  if(by < sp_cmatrix_rows(img->image)/2.0){
    by = (sp_cmatrix_rows(img->image))/2.0+by/2.0;
  }else{
    by = (sp_cmatrix_rows(img->image))/2.0-(sp_cmatrix_rows(img->image)-by)/2.0;
  }
  printf("Center x - %f y - %f\n",bx,by);
  write_png(a,"centrosym_convolve.png",COLOR_JET|LOG_SCALE);
  *center_x  = bx;
  *center_y = by;
  sp_image_free(a);
}

int pixel_to_index(Image * img, real * point){
  return ((int)point[0])*sp_cmatrix_rows(img->image)+point[1];
}


/* This doesn't really rescale the mask which is a problem and doesn't really downscale */
Image * fourier_rescale(Image * img, int new_x, int new_y){
  Image * res = sp_image_fft(img);
  Image * tmp;
  int i;
  real inv_size;
  if(new_x < sp_cmatrix_cols(img->image) ||
     new_y < sp_cmatrix_rows(img->image)){
    perror("fourier_scale doesn't downscale");
    abort();
  }
  tmp = zero_pad_image(res,new_x,new_y,1);
  sp_image_free(res);
  res = sp_image_ifft(tmp);
  sp_image_free(tmp);
  res->detector->image_center[0] = img->detector->image_center[0]*(new_x/sp_cmatrix_cols(img->image));
  res->detector->image_center[1] = img->detector->image_center[1]*(new_x/sp_cmatrix_rows(img->image));
  inv_size = 1.0/sp_image_size(img);
  for(i = 0;i<sp_image_size(res);i++){
    res->image->data[i] *= inv_size;
  }
  return res;
}


Image * bilinear_rescale(Image * img, int new_x, int new_y){
  Image * res = sp_image_duplicate(img,SP_COPY_DATA|SP_COPY_MASK);

  real virtual_x;
  real virtual_y;
  int x,y;
  res->detector->image_center[0] *= new_x/sp_cmatrix_cols(img->image);
  res->detector->image_center[1] *= new_y/sp_cmatrix_rows(img->image);
  sp_cmatrix_free(res->image);
  res->image = sp_cmatrix_alloc(new_y,new_x);
  sp_imatrix_free(res->mask);
  res->mask = sp_imatrix_alloc(new_y,new_x);
  

  for(x = 0; x<sp_cmatrix_cols(res->image); x++){
    virtual_x = (real)x*sp_cmatrix_cols(img->image)/sp_cmatrix_cols(res->image);
    for(y = 0; y<sp_cmatrix_rows(res->image); y++){
      virtual_y = (real)y*sp_cmatrix_rows(img->image)/sp_cmatrix_rows(res->image);
      res->image->data[x*sp_cmatrix_rows(res->image)+y] = sp_cmatrix_interp(img->image,virtual_y,virtual_x);	
      res->mask->data[x*sp_cmatrix_rows(res->image)+y] = sp_imatrix_interp(img->mask,virtual_y,virtual_x);
    }
  } 
  return res;
}
  

real sp_image_interp(Image * img, real v_x, real v_y){
  return creal(sp_cmatrix_interp(img->image,v_y,v_x));
}

void sp_image_scale(Image * img, real value){
  sp_cmatrix_scale(img->image,value);
}


real sp_image_max(Image * img, int * index,int * x, int * y){
  real ret;
  ret = sp_cmatrix_max(img->image,index);
  if(index){
    if(x){
      *x = *index/sp_image_height(img);
    }
    if(y){
      *y = *index%sp_image_height(img);
    }
  }
  return ret;
}


void sp_image_realloc(Image * img, int new_x, int new_y){
  sp_cmatrix_realloc(img->image,new_y,new_x);
  sp_imatrix_realloc(img->mask,new_y,new_x);
}

Image * rectangular_window(int image_x, int image_y, int width, int height, int shifted){
  Image * res = sp_image_alloc(image_x,image_y);
  int x,y,i;
  int center[2];
  i = 0;
  
  if(shifted){
    for(x = 0;x<image_x;x++){
      for(y = 0;y<image_y;y++){
	if((fabs(x) < width/2 || fabs(image_x-x) < width/2 )&&
	   (fabs(y) < height/2 || fabs(image_y-y) < height/2 )){
	  res->image->data[i] = 1;
	}else{
	  res->image->data[i] = 0;
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
	  res->image->data[i] = 1;
	}else{
	  res->image->data[i] = 0;
	}
	i++;
      }
    }
  }
  sp_image_rephase(res,SP_ZERO_PHASE);
  return res;
}


Image * circular_window(int x, int y, int radius, int shifted){
  Image * res = sp_image_alloc(x,y);
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
      res->image->data[i] = 1;
    }else{
      res->image->data[i] = 0;
    }
  }
  sp_image_rephase(res,SP_ZERO_PHASE);
  return res;
}


void sp_image_normalize(Image * in){
  double integral = 0;
  int i;
  for(i = 0;i<sp_image_size(in);i++){
    integral += cabs(in->image->data[i]);
  }
  for(i = 0;i<sp_image_size(in);i++){
    in->image->data[i] /= integral;
  }
}


Image * sp_image_get_mask(Image * a){
  Image * res = sp_image_duplicate(a,SP_COPY_DATA|SP_COPY_MASK);
  int i;
  for(i = 0;i<sp_image_size(a);i++){
    res->image->data[i] = res->mask->data[i];
  }
  return res;
}

real sp_point_convolute(Image * a, Image * b, int index){
  int index_x, index_y,x,y;
  double out = 0;
  int ai,bi;
  Image * tmp = NULL;
  if(a->shifted){
    tmp = sp_image_shift(a);
    index = sp_image_shift_index(a,index);
    a = tmp;
  }
  if(b->shifted){
    fprintf(stderr,"Point convoluting with a shifted function is not currently defined!\n");
    return 0;
  }
  index_x = index/sp_cmatrix_rows(a->image)-a->detector->image_center[0];
  index_y = index%sp_cmatrix_rows(a->image)-a->detector->image_center[1];
  for(x = -b->detector->image_center[0];x<sp_cmatrix_cols(b->image)-b->detector->image_center[0];x++){
    for(y = -b->detector->image_center[1];y<sp_cmatrix_rows(b->image)-b->detector->image_center[1];y++){
      if(x+index_x < -a->detector->image_center[0] || 
	 x+index_x >=  sp_cmatrix_cols(a->image)-a->detector->image_center[0] ||
	 y+index_y < -a->detector->image_center[1] || 
	 y+index_y >=  sp_cmatrix_rows(a->image)-a->detector->image_center[1]){
	/* we're outside of image a */
	continue;
      }
      ai = (index_x+x+a->detector->image_center[0])*sp_cmatrix_rows(a->image)+index_y+y+a->detector->image_center[1];
      bi = (x+b->detector->image_center[0])*sp_cmatrix_rows(b->image)+y+b->detector->image_center[1];
      out += a->image->data[ai]*b->image->data[bi];
    }
  }
  if(a->shifted){
    sp_image_free(tmp);
  }
  return out;
}

int sp_image_shift_index(Image * a, int index){
  int x = index/sp_cmatrix_rows(a->image);
  int y = index%sp_cmatrix_rows(a->image);
  if(a->shifted){
    x = (x+sp_cmatrix_cols(a->image)/2)%sp_cmatrix_cols(a->image);
    y = (y+sp_cmatrix_rows(a->image)/2)%sp_cmatrix_rows(a->image);
  }else{
    x -= a->detector->image_center[0];
    y -= a->detector->image_center[1];
    x += sp_cmatrix_cols(a->image);
    x %= sp_cmatrix_cols(a->image);
    y += sp_cmatrix_rows(a->image);
    y %= sp_cmatrix_rows(a->image);
  }
  return x*sp_cmatrix_rows(a->image)+y;
}

/*
  Returns the value of the distance between
  "point" (Or the img center in case "point" is NULL)
  and the border of the image in a certain direction
  specified by an angle (in radians)
  The point of intersection is returned
  in "intersection".
*/

real sp_image_distance_to_edge(Image * img, real * point, real direction, real * intersection){
  real d1,d2;
  real x1,y1,x2,y2;
  real * center;
  int dim[2];
  if(point){
    center = point;
  }else{
    center = img->detector->image_center;
  }
  dim[0] = sp_cmatrix_cols(img->image);
  dim[1] = sp_cmatrix_rows(img->image);
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
Image * sp_image_radial_sector(Image * img, real * point, real direction, int samples, real * intersection){
  int i;
  Image  *ret = sp_image_alloc(samples,1);
  real fpixel[2];
  real d_to_border;
  real * center;
  int dim[2];
  if(point){
    center = point;
  }else{
    center = img->detector->image_center;
  }
  dim[0] = sp_cmatrix_cols(img->image);
  dim[1] = sp_cmatrix_rows(img->image);

  d_to_border = sp_image_distance_to_edge(img,center,direction, intersection);
  for(i = 0;i<samples;i++){
    fpixel[0] = center[0]+cos(direction)*d_to_border*((real)i/samples);
    /* The - happens because the pixel 0,0 on an image is the upper left not the lower left */
    fpixel[1] = center[1]-sin(direction)*d_to_border*((real)i/samples);
    /* bilinear interpolation around fpixel */
    ret->image->data[i] = sp_cmatrix_interp(img->image,fpixel[1],fpixel[0]);
    /* bilinear interpolation around fpixel */
    ret->mask->data[i] = sp_imatrix_interp(img->mask,fpixel[1],fpixel[0]);
  }
  return ret;  
}


/*
  Takes a sector and rotates it around the center to create an image.
  It assumes the sector is pixel scaled (each bin 1 pixel ).
  That is it does not try to stretch the sector to cover the image.
*/

Image * sp_image_create_from_sector(Image * sector, int * img_size, real * center){
  int x,y;
  Image * ret = sp_image_alloc(img_size[0],img_size[1]);
  int bin;
  real d;
  real u;
  for(x = 0;x<img_size[0];x++){
    for(y = 0;y<img_size[1];y++){
      /* calculate distance to the center and sample from sector accordingly */
      d = sqrt((x-center[0])*(x-center[0])+(y-center[1])*(y-center[1]));
      bin = d;
      u = d-bin;
      if(d+1 < sp_cmatrix_cols(sector->image)-1){
	ret->image->data[x*img_size[1]+y] = (1.0-u)*sector->image->data[bin];
	ret->image->data[x*img_size[1]+y] += (u)*sector->image->data[bin+1];
	ret->mask->data[x*img_size[1]+y] = 1;
      }else{
	ret->image->data[x*img_size[1]+y] = 0;
	ret->mask->data[x*img_size[1]+y] = 0;
      }
    }
  }
  return ret;
}



Image * sp_image_local_variance(Image * img, Image * window){
  int i;
  int size[2] = {sp_cmatrix_cols(img->image)+sp_cmatrix_cols(window->image)-1,sp_cmatrix_rows(img->image)+sp_cmatrix_rows(window->image)-1};
  Image * norm_window = sp_image_duplicate(window,SP_COPY_DATA|SP_COPY_MASK);
  
  sp_image_normalize(norm_window);
  Image * ra = sp_image_convolute(img,norm_window,size);
  Image * res = rectangle_crop(ra,0,0,sp_cmatrix_cols(img->image)-1,sp_cmatrix_rows(img->image)-1);
  write_png(img,"non_averaged.png",COLOR_JET);
  write_png(ra,"total_averaged.png",COLOR_JET);
  write_png(res,"crop_total_averaged.png",COLOR_JET);
  for(i = 0;i<sp_image_size(res);i++){
    res->image->data[i] = fabs(res->image->data[i]-img->image->data[i]);
  }
  sp_image_free(ra);
  ra = sp_image_convolute(res,norm_window,size);
  sp_image_free(res);
  res = rectangle_crop(ra,0,0,sp_cmatrix_cols(img->image)-1,sp_cmatrix_rows(img->image)-1);
  sp_image_free(ra);
  sp_image_free(norm_window);

  return res;
}



Complex sp_image_dot_prod(Image * a, Image * b){
  return sp_cmatrix_froenius_prod(a->image,b->image);
}

Image * sp_proj_module(Image * a, Image * b){
  Image * ret = sp_image_duplicate(a,SP_COPY_DATA|SP_COPY_MASK);
  int i;
  for(i = 0;i<sp_image_size(a);i++){
    if(b->mask->data[i]){
      ret->image->data[i] *= cabs(b->image->data[i])/cabs(a->image->data[i]);
    }else{
      ret->image->data[i] = a->image->data[i];
    }
  }
  return ret;
}


Image * sp_proj_support(Image * a, Image * b){
  Image * ret = sp_image_duplicate(a,SP_COPY_DATA|SP_COPY_MASK);
  int i;
  for(i = 0;i<sp_image_size(a);i++){
    if(!b->image->data[i]){
      ret->image->data[i] = 0;
    }
  }
  return ret;
}


int sp_image_invert(Image * a){
  int i;
  for(i = 0;i<sp_image_size(a);i++){
    a->image->data[i] = 1.0/a->image->data[i];
  }
  return 0;
}


void sp_image_mul_elements(Image * a, Image * b){
  sp_cmatrix_mul_elements(a->image,b->image);  
}

void sp_image_conj(Image * a){  
  sp_cmatrix_conj(a->image);
}


void sp_image_memcpy(Image * dst,Image * src){  
  sp_imatrix_memcpy(dst->mask,src->mask);
  sp_cmatrix_memcpy(dst->image,src->image);
  memcpy(dst->detector,src->detector,sizeof(Detector));
}

void sp_image_transpose(Image * a){
  sp_cmatrix_transpose(a->image);
  sp_imatrix_transpose(a->mask);
}

/*! Inserts image from into image to at the position at_x, at_y
 *
 *  If from doesn't fit in to, the image is silently clipped.
 */
void sp_image_insert(Image * to, Image * from, int at_x, int at_y){
 int x;
  for(x = 0;x<sp_min(sp_image_width(from),sp_image_width(to)-at_x);x++){
    memcpy(&(to->image->data[sp_image_get_index(to,x+at_x,at_y)]),
	   &from->image->data[sp_image_get_index(from,x,0)],
	   sp_min(sp_image_height(from),sp_image_height(to)-at_x)*sizeof(Complex));
  }
}

/*! Extend an image outside its borders by radius pixels on each
 *  direction.
 *
 *
 * Shifted images are padded in the middle.
 * The edge_flags determines the kind of extension done.
 *
 * Possible values are:
 *
 * SP_ZERO_PAD_EDGE - The values outside the border are set to 0
 * SP_SYMMETRIC_EDGE - The values outside the bounds of the image
 * are computed by mirror-reflecting the image across the image border.
 * SP_REPLICATE_EDGE - The values outside the bounds of the image are 
 * assumed to equal the nearest image border value.
 * SP_CIRCULAR_EDGE - The values outside the bounds of the image are
 * computed by implicitly assuming the input image is periodic.
 */
Image * sp_image_edge_extend(Image * a, int radius, int edge_flags){
  Image * res = sp_image_alloc(sp_cmatrix_cols(a->image)+radius*2,sp_cmatrix_rows(a->image)+radius*2);
  int x,y,x0,y0;
  /* first put there the initial image */
  sp_image_insert(res,a,radius,radius);
  /* now fill the edges */
  if(edge_flags == SP_ZERO_PAD_EDGE){
    return res;
  }else if(edge_flags == SP_SYMMETRIC_EDGE){
    /* First do the four corners */
    for(x = 0,x0=x;x<radius;x++){
      for(y = 0,y0=y,y0=y;y<radius;y++){
	sp_image_set(res,x,y,sp_image_get(a,radius-x-1,radius-y-1));
      }
    }
    for(x = radius+sp_image_width(a),x0=x;x<sp_image_width(res);x++){
      for(y = 0,y0=y;y<radius;y++){
	sp_image_set(res,x,y,sp_image_get(a,sp_image_width(a)-(x-x0)-1,radius-y-1));
      }
    }

    for(x = radius+sp_image_width(a),x0=x;x<sp_image_width(res);x++){
      for(y = radius+sp_image_height(a),y0=y;y<sp_image_height(res);y++){
	sp_image_set(res,x,y,sp_image_get(a,sp_image_width(a)-(x-x0)-1,sp_image_height(a)-(y-y0)-1));
      }
    }

    for(x = 0,x0=x;x<radius;x++){
      for(y = radius+sp_image_height(a),y0=y;y<sp_image_height(res);y++){
	sp_image_set(res,x,y,sp_image_get(a,radius-(x-x0)-1,sp_image_height(a)-(y-y0)-1));
      }
    }
    /* And now the four sides */
    for(x = radius,x0=x;x<radius+sp_image_width(a);x++){
      for(y = 0,y0=y;y<radius;y++){
	sp_image_set(res,x,y,sp_image_get(a,(x-x0),radius-(y-y0)-1));
      }
    }
    for(x = radius,x0=x;x<radius+sp_image_width(a);x++){
      for(y = radius+sp_image_height(a),y0=y;y<sp_image_height(res);y++){
	sp_image_set(res,x,y,sp_image_get(a,(x-x0),sp_image_height(a)-(y-y0)-1));
      }
    }
    for(x = 0,x0=x;x<radius;x++){
      for(y = radius,y0=y;y<radius+sp_image_height(a);y++){
	sp_image_set(res,x,y,sp_image_get(a,radius-(x-x0)-1,(y-y0)));
      }
    }
    for(x = radius+sp_image_width(a),x0=x;x<sp_image_width(a);x++){
      for(y = radius,y0=y;y<radius+sp_image_height(a);y++){
	sp_image_set(res,x,y,sp_image_get(a,sp_image_width(a)-(x-x0)-1,(y-y0)));
      }
    }
  }else if(edge_flags == SP_REPLICATE_EDGE){
    /* First do the four corners */
    for(x = 0,x0=x;x<radius;x++){
      for(y = 0,y0=y,y0=y;y<radius;y++){
	sp_image_set(res,x,y,sp_image_get(a,0,0));
      }
    }
    for(x = radius+sp_image_width(a),x0=x;x<sp_image_width(res);x++){
      for(y = 0,y0=y;y<radius;y++){
	sp_image_set(res,x,y,sp_image_get(a,sp_image_width(a)-1,0));
      }
    }

    for(x = radius+sp_image_width(a),x0=x;x<sp_image_width(res);x++){
      for(y = radius+sp_image_height(a),y0=y;y<sp_image_height(res);y++){
	sp_image_set(res,x,y,sp_image_get(a,sp_image_width(a)-1,sp_image_height(a)-1));
      }
    }

    for(x = 0,x0=x;x<radius;x++){
      for(y = radius+sp_image_height(a),y0=y;y<sp_image_height(res);y++){
	sp_image_set(res,x,y,sp_image_get(a,0,sp_image_height(a)-1));
      }
    }
    /* And now the four sides */
    for(x = radius,x0=x;x<radius+sp_image_width(a);x++){
      for(y = 0,y0=y;y<radius;y++){
	sp_image_set(res,x,y,sp_image_get(a,(x-x0),0));
      }
    }
    for(x = radius,x0=x;x<radius+sp_image_width(a);x++){
      for(y = radius+sp_image_height(a),y0=y;y<sp_image_height(res);y++){
	sp_image_set(res,x,y,sp_image_get(a,(x-x0),sp_image_height(a)-1));
      }
    }
    for(x = 0,x0=x;x<radius;x++){
      for(y = radius,y0=y;y<radius+sp_image_height(a);y++){
	sp_image_set(res,x,y,sp_image_get(a,0,(y-y0)));
      }
    }
    for(x = radius+sp_image_width(a),x0=x;x<sp_image_width(a);x++){
      for(y = radius,y0=y;y<radius+sp_image_height(a);y++){
	sp_image_set(res,x,y,sp_image_get(a,sp_image_width(a)-1,(y-y0)));
      }
    }
    
  }else if(edge_flags == SP_CIRCULAR_EDGE){
    /* First do the four corners */
    for(x = 0,x0=x;x<radius;x++){
      for(y = 0,y0=y,y0=y;y<radius;y++){
	sp_image_set(res,x,y,sp_image_get(a,sp_image_width(a)-radius+(x-x0),sp_image_height(a)-radius+(y-y0)));
      }
    }
    for(x = radius+sp_image_width(a),x0=x;x<sp_image_width(res);x++){
      for(y = 0,y0=y;y<radius;y++){
	sp_image_set(res,x,y,sp_image_get(a,(x-x0),sp_image_height(a)-radius+(y-y0)));
      }
    }

    for(x = radius+sp_image_width(a),x0=x;x<sp_image_width(res);x++){
      for(y = radius+sp_image_height(a),y0=y;y<sp_image_height(res);y++){
	sp_image_set(res,x,y,sp_image_get(a,(x-x0),(y-y0)));
      }
    }

    for(x = 0,x0=x;x<radius;x++){
      for(y = radius+sp_image_height(a),y0=y;y<sp_image_height(res);y++){
	sp_image_set(res,x,y,sp_image_get(a,sp_image_width(a)-radius+(x-x0),(y-y0)));
      }
    }
    /* And now the four sides */
    for(x = radius,x0=x;x<radius+sp_image_width(a);x++){
      for(y = 0,y0=y;y<radius;y++){
	sp_image_set(res,x,y,sp_image_get(a,(x-x0),sp_image_height(a)-radius+(y-y0)));
      }
    }
    for(x = radius,x0=x;x<radius+sp_image_width(a);x++){
      for(y = radius+sp_image_height(a),y0=y;y<sp_image_height(res);y++){
	sp_image_set(res,x,y,sp_image_get(a,(x-x0),(y-y0)));
      }
    }
    for(x = 0,x0=x;x<radius;x++){
      for(y = radius,y0=y;y<radius+sp_image_height(a);y++){
	sp_image_set(res,x,y,sp_image_get(a,sp_image_height(a)-radius+(x-x0),(y-y0)));
      }
    }
    for(x = radius+sp_image_width(a),x0=x;x<sp_image_width(a);x++){
      for(y = radius,y0=y;y<radius+sp_image_height(a);y++){
	sp_image_set(res,x,y,sp_image_get(a,(x-x0),(y-y0)));
      }
    }
  }
  return  res;  
}

/* Filters the input image with a median filter
 *
 * The kernels tells how big the window is and how to weight the pixels.
 * The center of the center is equal to its dimensions/2.
 * The edge_flags correspond to the sp_image_edge_extend flags().
 */
void sp_image_median_filter(Image * a,sp_imatrix * kernel, int edge_flags){
  int integral = 0;
  int n,i,x,y;
  int kx,ky;
  int kcx = sp_imatrix_cols(kernel)/2;
  int kcy = sp_imatrix_rows(kernel)/2;

  for(i =0 ;i<sp_imatrix_size(kernel);i++){
    integral += kernel->data[i];
  }
  real * buffer = malloc(sizeof(real)*integral);
  int radius = sp_max((sp_imatrix_cols(kernel)+1)/2,(sp_imatrix_rows(kernel)+1)/2);
  Image * work = sp_image_edge_extend(a,radius,edge_flags);
  
  for(x = radius;x<radius+sp_image_width(a);x++){
    for(y = radius;y<radius+sp_image_height(a);y++){
      n = 0;
      for(kx = -kcx;kx<sp_imatrix_cols(kernel)-kcx;kx++){
	for(ky = -kcy;ky<sp_imatrix_rows(kernel)-kcy;ky++){
	  for(i = 0;i<sp_imatrix_get(kernel,ky+kcy,kx+kcx);i++){
	    buffer[n++] = sp_image_get(work,x+kx,y+ky);
	  }
	}
      }
      sp_bubble_sort(buffer,n);
      if(n%2){
	/* odd n take the middle one */
	sp_image_set(a,x-radius,y-radius,buffer[n/2]);
      }else{
	/* even n take the average */
	sp_image_set(a,x-radius,y-radius,(buffer[n/2]+buffer[n/2-1])/2);
      }
    }
  }    
}
 

/*! Contrast stretches an image
 *
 * The image is partitioned in x_div*y_div points. x_div and y_div must be at least 2
 * Each point is contrast stretched in relation to its neighbourhood
 * The scaling is linear interpolated from division to division
 */
void sp_image_adaptative_constrast_stretch(Image * a,int x_div, int y_div){
  int dx,dy,x,y;
  int x_div_size = sp_image_width(a)/(x_div-1);
  int y_div_size = sp_image_height(a)/(y_div-1);
  sp_matrix * div_max = sp_matrix_alloc(x_div,y_div);
  sp_matrix * div_min = sp_matrix_alloc(x_div,y_div);
  sp_matrix * div_mean = sp_matrix_alloc(x_div,y_div);
  sp_matrix * div_std_dev = sp_matrix_alloc(x_div,y_div);
  sp_matrix * offsets = sp_matrix_alloc(sp_image_width(a),sp_image_height(a));
  sp_matrix * factors = sp_matrix_alloc(sp_image_width(a),sp_image_height(a));
  int radius = sp_max(x_div_size/2,y_div_size/2);
  Image * work = sp_image_edge_extend(a, radius, SP_SYMMETRIC_EDGE);
  for(dx = 0;dx <x_div;dx++){
    for(dy = 0;dy <y_div;dy++){
      sp_matrix_set(div_max,dx,dy,sp_image_get(work,radius+dx*x_div_size,radius+dy*y_div_size));
      sp_matrix_set(div_min,dx,dy,sp_image_get(work,radius+dx*x_div_size,radius+dy*y_div_size));
      for(x = radius+dx*x_div_size;x<radius+(dx+1)*x_div_size;x++){	
	if(x == sp_image_width(work)){
	  break;
	}

	for(y = radius+dy*y_div_size;y<radius+(dy+1)*y_div_size;y++){
	  if(y == sp_image_height(work)){
	    break;
	  }
	  if(cabs(sp_image_get(work,x,y)) < sp_matrix_get(div_min,dx,dy)){
	    sp_matrix_set(div_min,dx,dy,cabs(sp_image_get(work,x,y)));
	  }
	  if(cabs(sp_image_get(work,x,y)) > sp_matrix_get(div_max,dx,dy)){
	    sp_matrix_set(div_max,dx,dy,cabs(sp_image_get(work,x,y)));
	  }
	  sp_matrix_set(div_mean,dx,dy,sp_image_get(work,x,y)+sp_matrix_get(div_mean,dx,dy));
	}
      }
      sp_matrix_set(div_mean,dx,dy,sp_matrix_get(div_mean,dx,dy)/(x_div_size*y_div_size));
      /* Second pass to calculate standard deviation */
      for(x = radius+dx*x_div_size;x<radius+(dx+1)*x_div_size;x++){	
	if(x == sp_image_width(work)){
	  break;
	}

	for(y = radius+dy*y_div_size;y<radius+(dy+1)*y_div_size;y++){
	  if(y == sp_image_height(work)){
	    break;
	  }
	  sp_matrix_set(div_std_dev,dx,dy,sp_matrix_get(div_std_dev,dx,dy)+
			cabs(sp_image_get(work,x,y)*sp_image_get(work,x,y)-
			     sp_matrix_get(div_mean,dx,dy)*sp_matrix_get(div_mean,dx,dy)));
	}
      }
      sp_matrix_set(div_std_dev,dx,dy,sqrt(sp_matrix_get(div_std_dev,dx,dy)/((x_div_size*y_div_size)-1)));
    }
  }
  sp_image_free(work);
  for(x = 0;x<sp_image_width(a);x++){
    for(y = 0;y<sp_image_height(a);y++){
      sp_matrix_set(offsets,x,y,-sp_matrix_interp(div_mean,sp_min((float)x/x_div_size,x_div-0.001),
						  sp_min((float)y/y_div_size,y_div-0.001))); 
      /* allow 3 standard deviation on each side of the mean */
      sp_matrix_set(factors,x,y,1.0/(6*sp_matrix_interp(div_std_dev,sp_min((float)x/x_div_size,x_div-0.001),
						      sp_min((float)y/y_div_size,y_div-0.001)))); 
    }
  }
  for(x = 0;x<sp_image_width(a);x++){
    for(y = 0;y<sp_image_height(a);y++){
      /* cap extreme values */
      if(creal(sp_image_get(a,x,y)) < -sp_matrix_get(offsets,x,y)-1/sp_matrix_get(factors,x,y)){
	sp_image_set(a,x,y,-sp_matrix_get(offsets,x,y)-1/sp_matrix_get(factors,x,y));
      }else if(creal(sp_image_get(a,x,y)) > -sp_matrix_get(offsets,x,y)+1/sp_matrix_get(factors,x,y)){
	sp_image_set(a,x,y,-sp_matrix_get(offsets,x,y)+1/sp_matrix_get(factors,x,y));
      }
      sp_image_set(a,x,y,(sp_image_get(a,x,y)+sp_matrix_get(offsets,x,y))*sp_matrix_get(factors,x,y));
    }
  }
  sp_matrix_free(offsets);
  sp_matrix_free(factors);
  sp_matrix_free(div_max);
  sp_matrix_free(div_min);
}


void sp_bubble_sort(real * a, int n){
  int i,j;
  int flag = 0;
  real tmp;

  for(i = 0;i<n-1;i++){
    flag = 0;
    for(j = 0;j<n-1-i;j++){
      if(a[j] > a[j+1]){
	tmp = a[j+1];
	a[j+1] = a[j];
	a[j] = tmp;
	flag = 1;
      }
    }
    if(!flag){
      break;
    }
  }
}


void sp_image_fourier_coords(Image * in, sp_matrix * k_x, sp_matrix * k_y, sp_matrix * k_z){
  /* We need to get the wavelength and detector size */
  /* Calculate the fourier coordinates of each pixel in the detector */
  /* First we project the pixels on the ewald sphere and then we calculate it's coordinates */

  /* number of pixels */
  int nx, ny;
  /* pixel index */
  int x,y;
  /* physical location of pixel*/
  real px,py;
  /* reciprocal coordinates */
  real rx,ry;
  real real_to_reciprocal = 1.0/(in->detector->detector_distance*in->detector->lambda);
  real ewald_radius = 1.0/in->detector->lambda;
  real distance_to_ewald_sphere_center;

  real det_width = in->detector->pixel_size * sp_image_width(in);
  real det_height = in->detector->pixel_size * sp_image_height(in);
  
  nx = sp_image_width(in);
  ny = sp_image_height(in);

  for(x = 0;x<nx;x++){
    for(y = 0;y<ny;y++){
      /* 
	 Calculate the pixel coordinates in reciprocal space 	 
	 by dividing the physical position by detector_distance*wavelength.
	 
	 CCD center at image_center(nx-1)/2,(ny-1)/2

	 Upper left corner of the detector with negative x and positive y
      */
      px = ((x-in->detector->image_center[0])/nx)*det_width;
      py = ((in->detector->image_center[1]-y)/ny)*det_height;

      rx = px*real_to_reciprocal;
      ry = py*real_to_reciprocal;
      /* Project pixel into Ewald sphere. */
      distance_to_ewald_sphere_center = sqrt(rx*rx+ry*ry+ewald_radius*ewald_radius);
      if(k_x){      
	sp_matrix_set(k_x,x,y,rx * ewald_radius/distance_to_ewald_sphere_center);
      }
      if(k_y){
	sp_matrix_set(k_y,x,y,ry * ewald_radius/distance_to_ewald_sphere_center);
      }
      if(k_z){
	sp_matrix_set(k_z,x,y,ewald_radius-(ewald_radius * ewald_radius/distance_to_ewald_sphere_center));
      }
    }
  }
}
