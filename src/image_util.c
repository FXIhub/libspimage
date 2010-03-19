/* #ifdef __MINGW32__ */
/* /\* resolve conflict with ssize_t*\/ */
/* #define _NO_OLDNAMES */
/* typedef long off_t; */

/* #ifndef _SSIZE_T_ */
/* #define _SSIZE_T_ */
/* #warning "defining ssize_t" */
/* typedef unsigned int ssize_t; */
/* #endif */

//#endif

#ifndef PNG_DEBUG
#  define PNG_DEBUG 3
#endif

#define _XOPEN_SOURCE 500

#include <stdlib.h>
#include <math.h>
#include <hdf5.h>
#include <tiffio.h>
#include <png.h>
#include <float.h>
#include <ctype.h>
#include <strings.h>
#ifdef _USE_DMALLOC
#include <dmalloc.h>
#endif

#include "spimage.h"




static Image * zero_pad_shifted_image(Image * a, int newx, int newy, int newz, int pad_mask);
static Image * zero_pad_unshifted_image(Image * a, int newx, int newy, int newz, int pad_mask);
static real dist_to_axis(int i, Image * in);
static real dist_to_center(int i, Image * in);
static real square_dist_to_center(int i, Image * in);
static real dist_to_corner(int i, Image * in);
static void random_rephase(Image *  img);
static Image * reflect_xy(Image * in, int in_place);
static Image * reflect_x(Image * in, int in_place);
static Image * reflect_y(Image * in, int in_place);
static Image * reflect_origo(Image * in, int in_place);
/*static void hsv_to_rgb(float H,float S,float V,float * R,float *G,float *B);

static void hsv_to_rgb(float H,float S,float V,float * R,float *G,float *B){
  if( V == 0 ){ *R 
= 0; *G = 0; *B = 0; 
  }else if( S == 0 ) {                                                                   
    *R = V;                                                            
    *G = V;                                                            
    *B = V;                                                            
  } else {                                                                   
    const double hf = H / 60.0;                                       
    const int    i  = (int) floor( hf );                              
    const double f  = hf - i;                                         
    const double pv  = V * ( 1 - S );                                 
    const double qv  = V * ( 1 - S * f );                             
    const double tv  = V * ( 1 - S * ( 1 - f ) );                     
    switch( i ){                                                               
    case 0: 
      *R = V; 
      *G = tv;
      *B = pv;
      break; 
    case 1:
      *R = qv;
      *G = V;
      *B = pv;
      break;
    case 2:
      *R = pv; 
      *G = V;
      *B = tv;
      break; 
    case 3: 
      *R = pv;
      *G = qv;
      *B = V;
      break;
    case 4:  
      *R = tv; 
      *G = pv;
      *B = V;
      break;  
    case 5:
      *R = V;
      *G = pv;
      *B = qv; 
      break;
    case 6: 
      *R = V;
      *G = tv;    
      *B = pv; 
      break; 
    case -1:  
      *R = V;
      *G = pv; 
      *B = qv;
      break;
    default:
      sp_error_fatal("i Value error in HSV to *R*G*B conversion, Value is %d",i);
      break;
    }									
  }									
  *R *= 255.0F;                                                        
  *G *= 255.0F;                                                        
  *B *= 255.0F;  
}
*/





void sp_srand(int i){
  srand(i);
}

real p_drand48(){
  real ret = ((double)rand())/(double)RAND_MAX;
  return (ret);
}


int sp_image_is_valid(const Image * a){
  int valid = 1;
  for(int i = 0;i<sp_image_size(a);i++){
    if(!isfinite(sp_real(a->image->data[i])) ||
       !isfinite(sp_imag(a->image->data[i])) ||
       !isfinite(sp_cabs(a->image->data[i]))){
      valid = 0;
      break;
    }
  }
  return valid;
}

/* Returns the patterson from diffraction image a */
Image * sp_image_patterson(Image * a){
  Image * b = sp_image_duplicate(a,SP_COPY_DATA|SP_COPY_MASK);
  Image * c;
  Image * d;

  if(b->scaled){
    sp_c3matrix_mul_elements(b->image,b->image);
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
  }else if(axis == SP_ORIGO){
    if(sp_image_z(in) == 1){
      return reflect_xy(in,in_place);
    }else{
      return reflect_origo(in,in_place);
    }
  }
  return NULL;
}

Image * sp_image_rotate(Image * in, SpAxis axis, SpAngle angleDef, int in_place){
  double angle = 0;
  if(angleDef == sp_0Degrees){
    return sp_image_duplicate(in,SP_COPY_DATA|SP_COPY_MASK); 
  }
  sp_matrix * rot = sp_matrix_alloc(2,2);
  if(angleDef == sp_90Degrees){
    angle = M_PI/2;
  }
  if(angleDef == sp_180Degrees){
    angle = M_PI;
  }  
  if(angleDef == sp_270Degrees){
    angle = 3*M_PI/2;
  }
  if(sp_image_x(in) != sp_image_y(in)){
    sp_error_fatal("Cannot rotate non square images, sorry.");
  }
  if(axis == sp_XAxis){
    sp_error_fatal("X axis rotation not implement yet, sorry.");
  }
  if(axis == sp_YAxis){
    sp_error_fatal("Y axis rotation not implement yet, sorry.");
  }
  sp_matrix_set(rot,0,0,cos(angle));
  sp_matrix_set(rot,0,1,sin(angle));
  sp_matrix_set(rot,1,0,-sin(angle));
  sp_matrix_set(rot,1,1,cos(angle));
  Image * out = sp_image_duplicate(in,SP_COPY_DATA|SP_COPY_MASK);
  sp_matrix * newx = sp_matrix_alloc(sp_image_x(in),sp_image_y(in));
  sp_matrix * newy = sp_matrix_alloc(sp_image_x(in),sp_image_y(in));
  int min_x = 1e9;
  int min_y = 1e9;
  for(int x = 0;x < sp_image_x(in);x++){
    for(int y = 0;y < sp_image_y(in);y++){
      int new_x = x*sp_matrix_get(rot,0,0)+y*sp_matrix_get(rot,0,1);
      int new_y = x*sp_matrix_get(rot,1,0)+y*sp_matrix_get(rot,1,1);
      sp_matrix_set(newx,x,y,new_x);
      sp_matrix_set(newy,x,y,new_y);
      if(min_x > new_x){
	min_x = new_x;
      }
      if(min_y > new_y){
	min_y = new_y;
      }
    }
  }
  for(int i = 0;i<sp_matrix_size(newx);i++){
    newx->data[i] -= min_x;
    newy->data[i] -= min_y;
  }
  for(int x = 0;x < sp_image_x(in);x++){
    for(int y = 0;y < sp_image_y(in);y++){
      sp_image_set(out,sp_matrix_get(newx,x,y),sp_matrix_get(newy,x,y),0,sp_image_get(in,x,y,0));
    }
  }
  sp_matrix_free(newx);
  sp_matrix_free(newy);
  sp_matrix_free(rot);
  return out;
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
  for(x = 0;x<sp_c3matrix_x(in->image);x++){
    for(y = 0;y<sp_c3matrix_y(in->image)/2.0;y++){
      y2 = sp_c3matrix_y(in->image)-y-1;
      tmp = sp_c3matrix_get(in->image,x,y,0);
      sp_c3matrix_set(out->image,x,y,0,sp_c3matrix_get(in->image,x,y2,0));
      sp_c3matrix_set(out->image,x,y2,0, tmp);
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
    out = sp_image_duplicate(in,SP_COPY_DATA|SP_COPY_MASK);
  }
  for(x = 0;x<sp_c3matrix_x(in->image)/2.0;x++){
    x2 = sp_c3matrix_x(in->image)-x-1;
    for(y = 0;y<sp_c3matrix_y(in->image);y++){
      tmp = sp_c3matrix_get(in->image,x,y,0);
      sp_c3matrix_set(out->image,x,y,0,sp_c3matrix_get(in->image,x2,y,0));
      sp_c3matrix_set(out->image,x2,y,0,tmp);
    }
  }
  return out;
}

static Image * reflect_origo(Image * in, int in_place){
  Complex tmp;
  Image * out;
  int x,y,z;
  if(in_place){
    out = in;
  }else{
    out = sp_image_duplicate(in,SP_COPY_DATA|SP_COPY_MASK);
  }
  for(z = 0;z<sp_c3matrix_z(in->image)/2.0;z++){
    for(y = 0; y<sp_c3matrix_y(in->image);y++){
      for(x = 0; x<sp_c3matrix_x(in->image);x++){
	tmp = sp_c3matrix_get(in->image,x,y,z);
	sp_c3matrix_set(out->image,x,y,z,sp_c3matrix_get(in->image,sp_c3matrix_x(in->image)-x-1,sp_c3matrix_y(in->image)-y-1,sp_c3matrix_z(in->image)-z-1));
	sp_c3matrix_set(out->image,sp_c3matrix_x(in->image)-x-1,sp_c3matrix_y(in->image)-y-1,sp_c3matrix_z(in->image)-z-1,tmp);
      }
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
  if(sp_c3matrix_x(in->image) != sp_c3matrix_y(in->image)){
    fprintf(stderr,"Error: Using average_centrosymetry with a non square image!\n");
    exit(1);
  }
  if(sp_c3matrix_z(in->image) != 1){
    fprintf(stderr,"Error: Using average_centrosymetry on 3D image");
  }

  for(x = 0;x<sp_c3matrix_x(in->image);x++){
    for(y = 0;y<sp_c3matrix_y(in->image);y++){
      ind1 = sp_c3matrix_get_index(in->image,sp_c3matrix_y(in->image)-1-y,sp_c3matrix_x(in->image)-1-x,0);
      ind2= sp_c3matrix_get_index(in->image,x,y,0);
      if(dist_to_corner(x*sp_c3matrix_y(in->image)+y,in) < 1500){
	if(sp_cabs(sp_cadd(in->image->data[ind1],in->image->data[ind2])) && in->mask->data[ind1] && in->mask->data[ind2]){
	  noise += sp_cabs(sp_csub(in->image->data[ind1],in->image->data[ind2]))/(sp_cabs(sp_cscale(sp_cadd(in->image->data[ind1],in->image->data[ind2]),0.5)));
	  k++;
	  out->image->data[ind2] = sp_cscale(sp_cadd(in->image->data[ind1],in->image->data[ind2]),0.5);
	}else{
	  sp_real(out->image->data[ind2]) = 0;
	  sp_imag(out->image->data[ind2]) = 0;
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
  if(sp_c3matrix_x(in->image) > sp_c3matrix_y(in->image)){
    /* remove a column */
    sp_c3matrix_free(out->image);
    out->image = sp_c3matrix_alloc(sp_c3matrix_y(in->image),sp_c3matrix_y(in->image),1);
    sp_i3matrix_free(out->mask);
    out->mask = sp_i3matrix_alloc(sp_c3matrix_y(in->image),sp_c3matrix_y(in->image),1);
    yout = 0;
    xout = 0;
    for(x = 0;x<sp_c3matrix_y(in->image);x++){
      for(y = 0;y<sp_c3matrix_y(in->image);y++){
	if(fabs(x - (sp_c3matrix_x(in->image)-1)/2.0) < (sp_c3matrix_x(in->image)-sp_c3matrix_y(out->image))/2.0){
	  continue;
	}
	sp_c3matrix_set(out->image,xout,yout,0,sp_c3matrix_get(in->image,x,y,0));
	sp_i3matrix_set(out->mask,xout,yout,0,sp_i3matrix_get(in->mask,x,y,0));
	yout = (yout+1)%sp_c3matrix_y(out->image);
	if(yout == 0){
	  xout++;
	}
      }
    }
  }else if(sp_c3matrix_x(in->image) < sp_c3matrix_y(in->image)){
    /* remove a line */
    sp_c3matrix_free(out->image);
    out->image = sp_c3matrix_alloc(sp_c3matrix_x(in->image),sp_c3matrix_x(in->image),1);
    sp_i3matrix_free(out->mask);
    out->mask = sp_i3matrix_alloc(sp_c3matrix_x(in->image),sp_c3matrix_x(in->image),1);

    for(x = 0;x<sp_c3matrix_x(in->image);x++){
      for(y = 0;y<sp_c3matrix_y(in->image);y++){
	if(fabs(x - (sp_c3matrix_y(in->image)-1)/2.0) < (sp_c3matrix_y(in->image)-sp_c3matrix_x(out->image))/2.0){
	  continue;
	}
	sp_c3matrix_set(out->image,xout,yout,0,sp_c3matrix_get(in->image,x,y,0));
	sp_i3matrix_set(out->mask,xout,yout,0,sp_i3matrix_get(in->mask,x,y,0));
	yout = (yout+1)%sp_c3matrix_y(out->image);
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
  if(sp_c3matrix_z(in->image) != 1){
    fprintf(stderr,"error: Trying to square a 3D image");
  }
  int size = sp_max(sp_c3matrix_x(in->image),sp_c3matrix_y(in->image));
  return zero_pad_image(in,size,size,1,0);
}


/*! This function calculates the size of the input vector
  after being shifted and possibly padded (depending on the value
  of the flag pad) */
static int shift_size(int size, int shift_origin, int pad){
  if(pad){
    return sp_max((size-shift_origin)*2,shift_origin*2);
  }else{
    return size;
  }
}
/*! This function simply takes as argument the index you want to shift (i),
  the size of the vector you want to shift (size), and the new origin (shift_origin).
  If you wish the new vector to have as many components to the right of the shift_origin
  as to the left of the shift_origin set pad to 1. 
*/
static int shift_coordinate(int i, int size, int shift_origin, int pad){
  if(!pad){
    return (i-shift_origin+size)%size;
  }else{
    int new_size = shift_size(size,shift_origin,pad);
    if(i<shift_origin){
      return new_size+(i-shift_origin);
    }else{
      return i-shift_origin;
    }
  }
  return -1;
}

Image * sp_image_shift(Image * img){
  Image * out;
  const int small_constant = 0.001;    
  int pad = 1;
  int pad_z = 1;
  if(img->num_dimensions == 2){
    pad_z = 0;
  }
  int new_origin[3];

  if(img->shifted){
    new_origin[0] = ceil((sp_image_x(img)-1.0)/2.0-small_constant);
    new_origin[1] = ceil((sp_image_y(img)-1.0)/2.0-small_constant);
    new_origin[2] = ceil((sp_image_z(img)-1.0)/2.0-small_constant);
  }else{
    new_origin[0] = ceil(img->detector->image_center[0]-small_constant);
    new_origin[1] = ceil(img->detector->image_center[1]-small_constant);
    new_origin[2] = ceil(img->detector->image_center[2]-small_constant);
  }
  int new_size[3] = {shift_size(sp_image_x(img),new_origin[0],pad),
		     shift_size(sp_image_y(img),new_origin[1],pad),
		     shift_size(sp_image_z(img),new_origin[2],pad_z)};
  /* Try to duplicate the original image as faithfully as possible */
  out = sp_image_duplicate(img,SP_COPY_ALL);
  sp_image_realloc(out,new_size[0],new_size[1],new_size[2]);
  out->shifted = !img->shifted;

  for(int i = 0;i<sp_image_size(out);i++){
    out->image->data[i] = sp_cinit(0,0);
    out->mask->data[i] = 1;
  }
  /* We're going to shift the image in all 3 dimensions by shifting each dimension individually */
  for(int z = 0;z<sp_image_z(img);z++){
    int new_z = shift_coordinate(z,sp_image_z(img),new_origin[2],pad_z);
    for(int y = 0;y<sp_image_y(img);y++){
      int new_y = shift_coordinate(y,sp_image_y(img),new_origin[1],pad);
      for(int x = 0;x<sp_image_x(img);x++){
	int new_x = shift_coordinate(x,sp_image_x(img),new_origin[0],pad);
	sp_image_set(out,new_x,new_y,new_z,sp_image_get(img,x,y,z));
	sp_image_mask_set(out,new_x,new_y,new_z,sp_image_mask_get(img,x,y,z));
      }
    }
  }		       
  if(img->shifted){
    out->detector->image_center[0] = (sp_image_x(img)-1.0)/2.0;
    out->detector->image_center[1] = (sp_image_y(img)-1.0)/2.0;
    out->detector->image_center[2] = (sp_image_z(img)-1.0)/2.0;
  }else{
    out->detector->image_center[0] = 0;
    out->detector->image_center[1] = 0;
    out->detector->image_center[2] = 0;
  }

  return out;
}

/*
  For unshifted images it shifted the quadrants around image_center and 
  extra zero pad is used on the smaller quadrants to make them the same size
  as the biggest quadrant.
  For shifted images it shifted the quadrants around (size[]-1)/2 
  This function is generalized for 3D patterns as well.
*/
Image * sp_image_shift2(Image * img){
  Image * out;
  int i;
  int index1,index2;
  int x,y,z;
  int newx,newy,newz;
  real max_x,max_y,max_z;//changed from int
  const int small_constant = 0.001;    
  int new_origin[3] = {ceil(img->detector->image_center[0]-small_constant),
		       ceil(img->detector->image_center[1]-small_constant),
		       ceil(img->detector->image_center[2]-small_constant)};
  


  /* for purposes of shifting the image the pixels which goes to the upper left of the image is
     the ceil(img->detector->image_center) pixel 
     
     A small constant is subtracted from img->detector->image_center to deal with numerical errors
     when the image_center is an integer (say 50).
  */

  /* fft shift the image */
  out = sp_image_duplicate(img,SP_COPY_DATA|SP_COPY_MASK);
  if(!img->shifted){
    max_x = sp_max(new_origin[0],sp_c3matrix_x(img->image)-1-new_origin[0]);
    max_y = sp_max(new_origin[1],sp_c3matrix_y(img->image)-1-new_origin[1]);
    max_z = sp_max(new_origin[2],sp_c3matrix_z(img->image)-1-new_origin[2]);
    //was 2*max_x+1 before (other way of defining center)
    /* FM BUG: There was a bug here for 2D images with max_z so I added this hack */
    /* With the new way of defining the top left pixel this is no longer required */
    /*    if(max_z == 0){
      max_z = 0.5;
      img->detector->image_center[2] = 0.5;
    }
    */
    sp_image_realloc(out,2*max_x,2*max_y,2*max_z);
  }

		   
  for(i = 0;i<sp_image_size(out);i++){
    sp_real(out->image->data[i]) = 0;
    sp_imag(out->image->data[i]) = 0;
    out->mask->data[i] = 0;
  }
  //doesn't work for shifted images with z=1. Could be helped by adding
  //a -1 when seting the center.
  //also a problem with different output size for shifted (i-1) and
  //unshifted(i) images (i is input side length). (even for 3D)
  if(img->shifted){
    out->detector->image_center[0] = (sp_c3matrix_x(img->image))/2.0;
    out->detector->image_center[1] = (sp_c3matrix_y(img->image))/2.0;
    out->detector->image_center[2] = (sp_c3matrix_z(img->image))/2.0;
    img->detector->image_center[0] = (sp_c3matrix_x(img->image))/2.0;
    img->detector->image_center[1] = (sp_c3matrix_y(img->image))/2.0;
    img->detector->image_center[2] = (sp_c3matrix_z(img->image))/2.0;
  }else{
    out->detector->image_center[0] = 0;
    out->detector->image_center[1] = 0;
    out->detector->image_center[2] = 0;
  }

  out->shifted = !out->shifted;
  /* shift quadrants */

  if(img->num_dimensions == 2){
    for(x = 0;x<sp_c3matrix_x(img->image);x++){
      for(y = 0;y<sp_c3matrix_y(img->image);y++){
	index1 = y*sp_c3matrix_x(img->image) + x;
	index2 = 0;
	if(y < img->detector->image_center[1]){
	  if(x < img->detector->image_center[0]){
	    newx = sp_c3matrix_x(out->image)-(img->detector->image_center[0]-x);
	    newy = sp_c3matrix_y(out->image)-(img->detector->image_center[1]-y);
	    if(newx < sp_c3matrix_x(img->image)/2 ||
	       newy < sp_c3matrix_y(img->image)/2){
	      index2 = -1;
	    }
	  }else{
	    newx = x-img->detector->image_center[0];
	    newy = sp_c3matrix_y(out->image)-(img->detector->image_center[1]-y);
	    if(newx >= sp_c3matrix_x(img->image)/2 ||
	       newy < sp_c3matrix_y(img->image)/2){
	      index2 = -1;
	    }
	  }
	}else{	
	  if(x < img->detector->image_center[0]){
	    newx = sp_c3matrix_x(out->image)-(img->detector->image_center[0]-x);
	    newy = y-img->detector->image_center[1];
	    if(newx < sp_c3matrix_x(img->image)/2 ||
	       newy >= sp_c3matrix_y(img->image)/2){
	      index2 = -1;
	    }
	  }else{
	    newx = x-img->detector->image_center[0];
	    newy = y-img->detector->image_center[1];
	    if(newx >= sp_c3matrix_x(img->image)/2 ||
	       newy >= sp_c3matrix_y(img->image)/2){
	      index2 = -1;
	    }	      
	  }
	}
	if(index2 != -1){
	  index2 = sp_c3matrix_get_index(out->image,newx,newy,0);
	}
	
	if(index2 != -1){
	  out->image->data[index2] = img->image->data[index1];
	  out->mask->data[index2] = img->mask->data[index1];
	}
      }
    }
  }else if(img->num_dimensions == 3){
    for(x = 0;x<sp_c3matrix_x(img->image);x++){
      for(y = 0;y<sp_c3matrix_y(img->image);y++){
	for(z = 0;z<sp_c3matrix_z(img->image);z++){
	  index1 = z*sp_c3matrix_y(img->image)*sp_c3matrix_x(img->image) +
	    y*sp_c3matrix_x(img->image) + x;
	  index2 = 0;
	  if(z < img->detector->image_center[2]){
	    if(y < img->detector->image_center[1]){
	      if(x < img->detector->image_center[0]){
		newx = sp_c3matrix_x(out->image)-(img->detector->image_center[0]-x);
		newy = sp_c3matrix_y(out->image)-(img->detector->image_center[1]-y);
		newz = sp_c3matrix_z(out->image)-(img->detector->image_center[2]-z);
		if(newx < sp_c3matrix_x(img->image)/2 ||
		   newy < sp_c3matrix_y(img->image)/2 ||
		   newz < sp_c3matrix_z(img->image)/2){
		  index2 = -1;
		}
	      }else{
		newx = x-img->detector->image_center[0];
		newy = sp_c3matrix_y(out->image)-(img->detector->image_center[1]-y);
		newz = sp_c3matrix_z(out->image)-(img->detector->image_center[2]-z);
		if(newx >= sp_c3matrix_x(img->image)/2 ||
		   newy < sp_c3matrix_y(img->image)/2 ||
		   newz < sp_c3matrix_z(img->image)/2){
		  index2 = -1;
		}
	      }
	    }else{	
	      if(x < img->detector->image_center[0]){
		newx = sp_c3matrix_x(out->image)-(img->detector->image_center[0]-x);
		newy = y-img->detector->image_center[1];
		newz = sp_c3matrix_z(out->image)-(img->detector->image_center[2]-z);
		if(newx < sp_c3matrix_x(img->image)/2 ||
		   newy >= sp_c3matrix_y(img->image)/2 ||
		   newz < sp_c3matrix_z(img->image)/2){
		  index2 = -1;
		}
	      }else{
		newx = x-img->detector->image_center[0];
		newy = y-img->detector->image_center[1];
		newz = sp_c3matrix_z(out->image)-(img->detector->image_center[2]-z);
		if(newx >= sp_c3matrix_x(img->image)/2 ||
		   newy >= sp_c3matrix_y(img->image)/2 ||
		   newz < sp_c3matrix_z(img->image)/2){
		  index2 = -1;
		}	      
	      }
	    }
	  }else{
	    if(y < img->detector->image_center[1]){
	      if(x < img->detector->image_center[0]){
		newx = sp_c3matrix_x(out->image)-(img->detector->image_center[0]-x);
		newy = sp_c3matrix_y(out->image)-(img->detector->image_center[1]-y);
		newz = z-img->detector->image_center[2];
		if(newx < sp_c3matrix_x(img->image)/2 ||
		   newy < sp_c3matrix_y(img->image)/2 ||
		   newz >= sp_c3matrix_z(img->image)/2){
		  index2 = -1;
		}
	      }else{
		newx = x-img->detector->image_center[0];
		newy = sp_c3matrix_y(out->image)-(img->detector->image_center[1]-y);
		newz = z-img->detector->image_center[2];
		if(newx >= sp_c3matrix_x(img->image)/2 ||
		   newy < sp_c3matrix_y(img->image)/2 ||
		   newz >= sp_c3matrix_z(img->image)/2){
		  index2 = -1;
		}
	      }
	    }else{	
	      if(x < img->detector->image_center[0]){
		newx = sp_c3matrix_x(out->image)-(img->detector->image_center[0]-x);
		newy = y-img->detector->image_center[1];
		newz = z-img->detector->image_center[2];
		if(newx < sp_c3matrix_x(img->image)/2 ||
		   newy >= sp_c3matrix_y(img->image)/2 ||
		   newz >= sp_c3matrix_z(img->image)/2){
		  index2 = -1;
		}
	      }else{
		newx = x-img->detector->image_center[0];
		newy = y-img->detector->image_center[1];
		newz = z-img->detector->image_center[2];
		if(newx >= sp_c3matrix_x(img->image)/2 ||
		   newy >= sp_c3matrix_y(img->image)/2 ||
		   newz >= sp_c3matrix_z(img->image)/2){
		  index2 = -1;
		}	      
	      }
	    }
	  }
	  if(index2 != -1){
	    index2 = sp_c3matrix_get_index(out->image,newx,newy,newz);
	  }
	  
	  if(index2 != -1){
	    out->image->data[index2] = img->image->data[index1];
	    out->mask->data[index2] = img->mask->data[index1];
	  }
	}
      }
    }
  }  
  return out;
}


/* resolution given in pixels from center. 
   Image should be shifted. Square window. */
Image * sp_image_low_pass(Image * img, int resolution, int type){
  Image * res = sp_image_duplicate(img,SP_COPY_DATA|SP_COPY_MASK);
  int x,y,z,nx,ny,nz;
  int dx,dy,dz;
  if(img->shifted == 0){
    fprintf(stderr,"Error: Trying to limit resolution on an unshifted image\n");
  }
  if(img->num_dimensions == SP_2D){
    if(resolution*2 > sp_c3matrix_x(res->image) || resolution*2 > sp_c3matrix_y(res->image)){
      return sp_image_duplicate(img,SP_COPY_DATA|SP_COPY_MASK);
    }  
    sp_c3matrix_free(res->image);
    sp_i3matrix_free(res->mask);
    res->image = sp_c3matrix_alloc(resolution*2,resolution*2,1);
    res->mask = sp_i3matrix_alloc(resolution*2,resolution*2,1);
    nx = 0;
    ny = 0;
    for(x = 0;x<sp_c3matrix_x(img->image);x++){
      dx = x;
      if(sp_c3matrix_x(img->image)-x-1 < dx){
      dx = sp_c3matrix_x(img->image)-x-1;
      }
      for(y = 0;y<sp_c3matrix_y(img->image);y++){
	dy = y;
	if(sp_c3matrix_y(img->image)-y-1 < dy){
	  dy = sp_c3matrix_y(img->image)-y-1;
	}
	if(dx < resolution && dy < resolution){
	  sp_c3matrix_set(res->image,nx,ny,0,sp_c3matrix_get(img->image,x,y,0));
	  sp_i3matrix_set(res->mask,nx,ny,0,sp_i3matrix_get(img->mask,x,y,0));
	  ny++;
	  if(ny == sp_c3matrix_y(res->image)){
	    nx++;
	    ny = 0;
	  }
	}  
      }
    }
  }else{
    if(resolution*2 > sp_c3matrix_x(res->image) ||
       resolution*2 > sp_c3matrix_y(res->image) ||
       resolution*2 > sp_c3matrix_z(res->image)){
      return sp_image_duplicate(img,SP_COPY_DATA|SP_COPY_MASK);
    }
    sp_c3matrix_free(res->image);
    sp_i3matrix_free(res->mask);
    res->image = sp_c3matrix_alloc(resolution*2,resolution*2,resolution*2);
    res->mask = sp_i3matrix_alloc(resolution*2,resolution*2,resolution*2);
    nx = 0;
    ny = 0;
    nz = 0;
    for(x = 0;x<sp_c3matrix_x(img->image);x++){
      dx = x;
      if(sp_c3matrix_x(img->image)-x-1 < dx){
	dx = sp_c3matrix_x(img->image)-x-1;
      }
      for(y = 0;y<sp_c3matrix_y(img->image);y++){
	dy = y;
	if(sp_c3matrix_y(img->image)-y-1 < dy){
	  dy = sp_c3matrix_y(img->image)-y-1;
	}
	for(z = 0;z<sp_c3matrix_z(img->image);z++){
	  dz = z;
	  if(sp_c3matrix_z(img->image)-z-1 < dz){
	    dz = sp_c3matrix_z(img->image)-z-1;
	  }
	  if(dx < resolution && dy < resolution && dz < resolution){
	    sp_c3matrix_set(res->image,nx,ny,nz,sp_c3matrix_get(img->image,x,y,z));
	    sp_i3matrix_set(res->mask,nx,ny,nz,sp_i3matrix_get(img->mask,x,y,z));
	    nz++;
	    if(nz == sp_c3matrix_z(res->image)){
	      ny++;
	      nz = 0;
	      if(ny == sp_c3matrix_y(res->image)){
		nx++;
		ny = 0;
	      }
	    }
	  }
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
void sp_image_smooth_edges(Image * img, sp_i3matrix * mask, int flags, real * value){
  int i;
  Image * tmp = sp_image_duplicate(img,SP_COPY_DATA|SP_COPY_MASK);
  if(flags & SP_GAUSSIAN){
    for(i = 0;i<sp_image_size(tmp);i++){
      sp_real(tmp->image->data[i]) = mask->data[i];
      sp_imag(tmp->image->data[i]) = 0;
    }
    printf("img (%i,%i,%i)\n",sp_image_x(img),sp_image_y(img),sp_image_z(img));
    Image * blur_mask = sp_gaussian_blur(tmp,*value);

    printf("blur_mask (%i,%i,%i)\n",sp_image_x(blur_mask),sp_image_y(blur_mask),sp_image_z(blur_mask));
    /* eat out the edge of the mask*/
    for(i = 0;i<sp_image_size(blur_mask);i++){
      if(sp_real(blur_mask->image->data[i]) < 0.99){
	sp_real(blur_mask->image->data[i]) = 0;
	sp_imag(blur_mask->image->data[i]) = 0;
      }
    }
      
    sp_image_free(tmp);
    tmp = blur_mask;
    blur_mask = sp_gaussian_blur(tmp,*value);
    printf("blur_mask (%i,%i,%i)\n",sp_image_x(blur_mask),sp_image_y(blur_mask),sp_image_z(blur_mask));

    sp_c3matrix_mul_elements(img->image,blur_mask->image);
  }  
}

real sp_centro_sym_value(long long index,Image * img){
  int x,y,z;
  real nx,ny,nz;
  x = index/sp_c3matrix_y(img->image)/sp_c3matrix_z(img->image);
  y = index/sp_c3matrix_z(img->image)%sp_c3matrix_y(img->image);
  z = index%sp_c3matrix_z(img->image)%sp_c3matrix_y(img->image);
  nx = 2*img->detector->image_center[0]-x;
  ny = 2*img->detector->image_center[1]-y;
  nz = 2*img->detector->image_center[2]-z;
  if(nx < 1 || nx >= sp_c3matrix_x(img->image)-2 ||
     ny < 1 || ny >= sp_c3matrix_y(img->image)-2 ||
     nz < 1 || nz >= sp_c3matrix_z(img->image)-2){
    return -1;
  }
  return sp_image_interp(img, nx, ny, nz);
}

int sp_centro_sym_index(long long index,Image * img){
  int x,y,z;
  int nx,ny,nz;
  x = index/sp_c3matrix_y(img->image)/sp_c3matrix_z(img->image);
  y = index/sp_c3matrix_z(img->image)%sp_c3matrix_y(img->image);
  z = index%sp_c3matrix_z(img->image)%sp_c3matrix_y(img->image);
  nx = 2*img->detector->image_center[0]-x;
  ny = 2*img->detector->image_center[1]-y;
  nz = 2*img->detector->image_center[2]-z;
  if(nx < 0 || nx >= sp_c3matrix_x(img->image) ||
     ny < 0 || ny >= sp_c3matrix_y(img->image) ||
     nz < 0 || nz >= sp_c3matrix_z(img->image)){
    return index;
  }
  return nz*sp_c3matrix_x(img->image)*sp_c3matrix_y(img->image)+ny*sp_c3matrix_x(img->image)+nx;
}

Image * sp_centro_sym_correlation(Image  * img){
  int x,y,z;
  Image * res = sp_image_duplicate(img,SP_COPY_DATA|SP_COPY_MASK);
  int index = 0;
  real csvalue;
  for(x = 0;x<sp_c3matrix_x(img->image);x++){
    for(y = 0;y<sp_c3matrix_y(img->image);y++){
      for(z = 0;z<sp_c3matrix_z(img->image);z++){
	csvalue = sp_centro_sym_value(index,img);
	if(!img->mask->data[index] || csvalue == -1 || sp_cabs(img->image->data[index])+fabs(csvalue) < 1){
	  sp_real(res->image->data[index]) = 1.0;
	  sp_imag(res->image->data[index]) = 0;
	}else{
	  Complex tmp = img->image->data[index];
	  sp_real(tmp) -= csvalue;
	  sp_real(res->image->data[index]) = 1.0 - sp_cabs(tmp)/(sp_cabs(img->image->data[index])+fabs(csvalue));      
	  sp_imag(res->image->data[index]) = 0;
	}
	if(sp_real(res->image->data[index]) < 0 || sp_real(res->image->data[index]) > 1){
	  /* Houston we have a problem */
	  exit(1);
	}
	index++;
      }
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
  float x,y,z;
  sp_image_get_coords_from_index(in,i,&x,&y,&z,SpTopLeftCorner);
  real dx,dy,dz;
  if(in->shifted){
    dx = MIN(x,sp_c3matrix_x(in->image)-x);
    dy = MIN(y,sp_c3matrix_y(in->image)-y);
    dz = MIN(z,sp_c3matrix_z(in->image)-z);
  }else{
    dx = fabs(x-in->detector->image_center[0]);
    dy = fabs(y-in->detector->image_center[1]);
    dz = fabs(z-in->detector->image_center[2]);
  }
  return MIN(dx,MIN(dy,dz));
}



static real dist_to_center(int i, Image * in){
  float x,y,z;
  sp_image_get_coords_from_index(in,i,&x,&y,&z,SpTopLeftCorner);
  real dx,dy,dz;
  if(in->shifted){
    dx = MIN(x,sp_c3matrix_x(in->image)-x);
    dy = MIN(y,sp_c3matrix_y(in->image)-y);
    dz = MIN(z,sp_c3matrix_z(in->image)-z);
  }else{
    dx = x-in->detector->image_center[0];
    dy = y-in->detector->image_center[1];
    dz = z-in->detector->image_center[2];
  }
  return sqrt(dx*dx+dy*dy+dz*dz);
}

static real square_dist_to_center(int i, Image * in){
  float x,y,z;
  sp_image_get_coords_from_index(in,i,&x,&y,&z,SpTopLeftCorner);
  real dx,dy,dz;
  if(in->shifted){
    dx = MIN(x,sp_c3matrix_x(in->image)-x);
    dy = MIN(y,sp_c3matrix_y(in->image)-y);
    dz = MIN(z,sp_c3matrix_z(in->image)-z);
  }else{
    dx = fabs(x-in->detector->image_center[0]);
    dy = fabs(y-in->detector->image_center[1]);
    dz = fabs(z-in->detector->image_center[2]);
  }
  return MAX(dx,MAX(dy,dz));

}

static real dist_to_corner(int i, Image * in){
  float x,y,z;
  sp_image_get_coords_from_index(in,i,&x,&y,&z,SpTopLeftCorner);
  real dx,dy,dz;
  if(sp_c3matrix_x(in->image)-1-x < x){
    dx = sp_c3matrix_x(in->image)-1-x;
  }else{
    dx = x;
  }
  if(sp_c3matrix_y(in->image)-1-y < y){
    dy = sp_c3matrix_y(in->image)-1-y;
  }else{
    dy = y;
  }
  if(sp_c3matrix_z(in->image)-1-z < z){
    dz = sp_c3matrix_z(in->image)-1-z;
  }else{
    dz = z;
  }
  return sqrt(dx*dx+dy*dy+dz*dz);
}


void sp_image_dephase(Image *  img){
  int i = 0;
  img->phased = 0;    
  for(i = 0;i<sp_image_size(img);i++){
    sp_real(img->image->data[i]) = sp_cabs(img->image->data[i]);
    sp_imag(img->image->data[i]) = 0;
  }
}

void sp_image_rephase(Image *  img, int type){ 
  img->phased = 1;
  if(type == SP_ZERO_PHASE){
    for(int i = 0;i<sp_image_size(img);i++){
      img->image->data[i] = sp_cinit(sp_cabs(img->image->data[i]),0);
    }
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
    real amp = sp_cabs(img->image->data[i]);
    sp_real(img->image->data[i]) = cos(phase)*amp;
    sp_imag(img->image->data[i]) = sin(phase)*amp;
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
    sp_real(res->image->data[i]) = sp_carg(img->image->data[i]);
    sp_imag(res->image->data[i]) = 0;
  }
  return res;
}





void sp_image_add(Image * a, const Image * b){
  sp_c3matrix_add(a->image,b->image,NULL);
}

void sp_image_sub(Image * a,const Image * b){
  sp_c3matrix_sub(a->image,b->image);
}


void sp_image_fill(Image * a, Complex value){
  for(int i = 0;i<sp_image_size(a);i++){
    a->image->data[i] = value;
  }
}

Complex sp_image_integrate(const Image * a){
  size_t size = sp_image_size(a);
  Complex ret = {0,0};
  for(size_t i = 0;i<size;i++){
    sp_cincr(ret,a->image->data[i]);
  }
  return ret;
}

real sp_image_integrate2(const Image * a){
  size_t size = sp_image_size(a);
  double ret = 0;
  for(size_t i = 0;i<size;i++){
    ret += sp_cabs2(a->image->data[i]);
  }
  return ret;
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
      sp_real(in->image->data[i]) = (1.0+sp_box_muller(0,level))*sp_real(in->image->data[i]);
      sp_imag(in->image->data[i]) = (1.0+sp_box_muller(0,level))*sp_imag(in->image->data[i]);
    }
  }
}

void sp_image_high_pass(Image * in, real radius, int type){
  int x,y,z;
  real dist,dx,dy,dz;
  if(radius <= 0){
    return;
  }
  if(type == SP_3D){
    for(x = 0;x<sp_c3matrix_x(in->image);x++){
      for(y = 0;y<sp_c3matrix_y(in->image);y++){
	for(z = 0;z<sp_c3matrix_z(in->image);z++){
	  if(x > sp_c3matrix_x(in->image)/2.0){
	    dx = sp_c3matrix_x(in->image)-x;
	  }else{
	    dx = x;
	  }
	  if(y > sp_c3matrix_y(in->image)/2.0){
	    dy = sp_c3matrix_y(in->image)-y;
	  }else{
	    dy = y;
	  }
	  if(z > sp_c3matrix_z(in->image)/2.0){
	    dz = sp_c3matrix_z(in->image)-z;
	  }else{
	    dz = z;
	  }
	  dist = sqrt(dx*dx+dy*dy+dz*dz);
	  if(dist <= radius){
	    in->mask->data[z*sp_c3matrix_x(in->image)*sp_c3matrix_y(in->image)+
			   y*sp_c3matrix_x(in->image)+x] = 0;
	    sp_real(in->image->data[z*sp_c3matrix_x(in->image)*sp_c3matrix_y(in->image)+
				    y*sp_c3matrix_x(in->image)+y]) = 0;
	    sp_imag(in->image->data[z*sp_c3matrix_x(in->image)*sp_c3matrix_y(in->image)+
				    y*sp_c3matrix_x(in->image)+y]) = 0;
	  }
	}
      }
    }
  }else if(in->num_dimensions == SP_2D){
    for(x = 0;x<sp_c3matrix_x(in->image);x++){
      for(y = 0;y<sp_c3matrix_y(in->image);y++){
	if(x > sp_c3matrix_x(in->image)/2.0){
	  dx = sp_c3matrix_x(in->image)-x;
	}else{
	  dx = x;
	}
	if(y > sp_c3matrix_y(in->image)/2.0){
	  dy = sp_c3matrix_y(in->image)-y;
	}else{
	  dy = y;
	}
	dist = sqrt(dx*dx+dy*dy);
	if(dist <= radius){
	  in->mask->data[y*sp_c3matrix_x(in->image)+x] = 0;
	  sp_real(in->image->data[y*sp_c3matrix_x(in->image)+x]) = 0;
	  sp_imag(in->image->data[y*sp_c3matrix_x(in->image)+x]) = 0;
	}
      }
    }
  }
}

void _sp_image_free(Image * in, const char * file, int line){
  _sp_c3matrix_free(in->image,file,line);
  _sp_i3matrix_free(in->mask,file,line);
#ifdef _DEBUG_MEM
  _sp_free(in->detector,file,line);
  _sp_free(in,file,line);
#else
  sp_free(in->detector);
  sp_free(in);
#endif

}

Image * _sp_image_duplicate(const Image * in, int flags,const char * file, int line){
  Image  *res = sp_malloc(sizeof(Image));
  if(!res){
    sp_error_fatal("Out of memory!");
  }
  memcpy(res,in,sizeof(Image));
  res->detector = sp_malloc(sizeof(Detector));
  if(!res->detector){
    sp_error_fatal("Out of memory!");
  }

  memcpy(res->detector,in->detector,sizeof(Detector));
  res->image = _sp_c3matrix_alloc(sp_c3matrix_x(in->image),sp_c3matrix_y(in->image),
				 sp_c3matrix_z(in->image),file,line);
  if(!res->image){
    sp_error_fatal("Out of memory!");
  }
  if(flags & SP_COPY_DATA){
    sp_c3matrix_memcpy(res->image,in->image);
  }

  res->mask = _sp_i3matrix_alloc(sp_c3matrix_x(in->image),sp_c3matrix_y(in->image),
				 sp_c3matrix_z(in->image),file,line);
  if(!res->mask){
    sp_error_fatal("Out of memory!");
  }
  if(flags & SP_COPY_MASK){
    sp_i3matrix_memcpy(res->mask,in->mask);
  }
  return res;
}


Image * _sp_image_alloc(int x, int y, int z,const char * file, int line){
  Image  *res = sp_malloc(sizeof(Image));
  if(!res){
    perror("Out of memory!\n");
    abort();
  }
  res->detector = sp_malloc(sizeof(Detector));
  res->detector->image_center[0] = (x-1)/2.0;
  res->detector->image_center[1] = (y-1)/2.0;
  res->detector->image_center[2] = (z-1)/2.0;
  res->detector->orientation = NULL;
  res->mask = _sp_i3matrix_alloc(x,y,z,file,line);
  res->scaled = 0;
  res->shifted = 0;
  res->image = _sp_c3matrix_alloc(x,y,z,file,line);
  res->phased = 0;
  res->rec_coords = 0;

  if(z == 1){
    /*assume we have a 2D image */
    res->num_dimensions = SP_2D;
  }else if(z > 1){
    res->num_dimensions = SP_3D;
  }else{
    perror("Dimensions must be positive integers");
    abort();
  }
  return res;
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
  int z;
  int i;
  Image * res;
  Image * a_ft;
  Image * b_ft;
  Image * tmp;

  if(sp_c3matrix_x(a->image) < sp_c3matrix_x(b->image)){
    /* swap */
    res = a;
    a = b;
    b = res;
  }

  if(size){
    x = size[0];
    y = size[1];
    z = size[2];
  }else{
    x = sp_c3matrix_x(a->image);
    y = sp_c3matrix_y(a->image);
    z = sp_c3matrix_z(a->image);
  }

  tmp = zero_pad_image(a,x,y,z,1);
  a_ft = sp_image_fft(tmp);
  sp_image_free(tmp);

  tmp = zero_pad_image(b,x,y,z,1);
  b_ft = sp_image_fft(tmp);
  sp_image_free(tmp);

  tmp = sp_image_duplicate(a_ft,SP_COPY_DETECTOR);
  tmp->shifted = 1;
  sp_image_rephase(tmp,SP_ZERO_PHASE);
  /* Now do the multiplication in fourier space */
  /* Using the Convolution Theorem */
  for(i = 0;i<x*y*z;i++){
    tmp->image->data[i] = sp_cmul(a_ft->image->data[i],sp_cconj(b_ft->image->data[i]));
  }

  sp_image_free(a_ft);
  sp_image_free(b_ft);
  /* Backtransform */
  res = sp_image_ifft(tmp);
  sp_image_free(tmp);

  /* should be all real */
  for(i = 0;i<sp_image_size(res);i++){
    res->image->data[i] = sp_cscale(res->image->data[i],1.0/sp_image_size(res));
    
  }
  return res;  
}


/* Calculates the convolution (or correlation if flag != 0) of a with b.
 *
 * The convolution is calculated with better than pixel precision. 
 * precision is the inverse of the precision of the output. precision == 2 
 * corresponds to output with 1/2 pixel precision.
 * Flag tells whether to calculate convolution (flag == 0) or correlation (flag != 0)
 * The other parameters are similar to the normal convolution
 */
Image * sp_image_convolute_fractional(Image * a, Image * b, int * size,int precision, int flag){
  int x;
  int y;
  int z;
  int i;
  Image * res;
  Image * a_ft;
  Image * b_ft;
  Image * tmp;

  if(sp_c3matrix_x(a->image) < sp_c3matrix_x(b->image)){
    /* swap */
    res = a;
    a = b;
    b = res;
  }

  /* special treatment of 2D images */
  int precision_z = precision;
  if(a->num_dimensions == SP_2D){
    precision_z = 1;
  }

  if(size){
    x = size[0];
    y = size[1];
    z = size[2];
  }else{
    x = sp_c3matrix_x(a->image);
    y = sp_c3matrix_y(a->image);
    z = sp_c3matrix_z(a->image);
  }

  tmp = zero_pad_image(a,x,y,z,1);
  a_ft = sp_image_fft(tmp);
  sp_image_free(tmp);
  tmp = zero_pad_image(a_ft,x*precision,y*precision,z*precision_z,1);
  sp_image_free(a_ft);
  a_ft = tmp;
  

  tmp = zero_pad_image(b,x,y,z,1);
  b_ft = sp_image_fft(tmp);
  sp_image_free(tmp);
  tmp = zero_pad_image(b_ft,x*precision,y*precision,z*precision_z,1);
  sp_image_free(b_ft);
  b_ft = tmp;


  tmp = sp_image_duplicate(a_ft,SP_COPY_DETECTOR);
  tmp->shifted = 1;
  sp_image_rephase(tmp,SP_ZERO_PHASE);
  /* Now do the multiplication in fourier space */
  /* Using the Convolution Theorem */
  for(i = 0;i<sp_image_size(tmp);i++){
    if(flag){
      tmp->image->data[i] = sp_cmul(a_ft->image->data[i],sp_cconj(b_ft->image->data[i]));
    }else{
      tmp->image->data[i] = sp_cmul(a_ft->image->data[i],b_ft->image->data[i]);
    }
  }

  sp_image_free(a_ft);
  sp_image_free(b_ft);
  /* Backtransform */
  res = sp_image_ifft(tmp);
  sp_image_free(tmp);

  /* should be all real */
  /* Note well that we're using the size of the INPUT to scale! */
  for(i = 0;i<sp_image_size(res);i++){
    res->image->data[i] = sp_cscale(res->image->data[i],1.0/sp_image_size(a));    
  }

  res->detector->image_center[0] = 0;
  res->detector->image_center[1] = 1;
  res->detector->image_center[2] = 2;
  res->shifted = 1;

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

/* FM: BUG! THIS FUNCTION IS BADLY CODED AND NEEDS TO BE REDONE! */
Image * sp_image_convolute(Image * a, Image * b, int * size){
  int x;
  int y;  
  int z;
  int i;
  Image * res;
  Image * a_ft;
  Image * b_ft;
  Image * tmp;

  
  /*if(sp_c3matrix_x(a->image) < sp_c3matrix_x(b->image)){*/
    /* swap */
  /*    res = a;
    a = b;
    b = res;
  }
  
  if(!size){
    x = sp_c3matrix_x(a->image);
    y = sp_c3matrix_y(a->image);
    z = sp_c3matrix_z(a->image);
  }else{
    x = size[0];
    y = size[1];
    z = size[2];
  }
  
  tmp = zero_pad_image(a,x,y,z,1);
  */

  x = sp_max(sp_image_x(a),sp_image_x(b));
  y = sp_max(sp_image_y(a),sp_image_y(b));
  z = sp_max(sp_image_z(a),sp_image_z(b));

  tmp = zero_pad_image(a,x,y,z,1);

/*  sp_image_dephase(tmp);*/
  a_ft = sp_image_fft(tmp);
  sp_image_free(tmp);

  tmp = zero_pad_image(b,x,y,z,1);
/*  sp_image_dephase(tmp);*/
  b_ft = sp_image_fft(tmp);
  sp_image_free(tmp);

  tmp = sp_image_duplicate(a_ft,SP_COPY_DETECTOR);
  tmp->shifted = 1;
/*  sp_image_rephase(tmp,SP_ZERO_PHASE);*/
  /* Now do the multiplication in fourier space */
  /* Using the Convolution Theorem */
  for(i = 0;i<x*y*z;i++){
    tmp->image->data[i] = sp_cmul(a_ft->image->data[i],b_ft->image->data[i]);
  }

  sp_image_free(a_ft);
  sp_image_free(b_ft);
  /* Backtransform */
  res = sp_image_ifft(tmp);
  sp_image_free(tmp);

  /*crop to make the image the original size*/
/*  tmp = cube_crop(res,
		     (sp_image_x(res)-sp_image_x(a))/2,
		     (sp_image_y(res)-sp_image_y(a))/2,
		     (sp_image_z(res)-sp_image_z(a))/2,
		     (sp_image_x(res)+sp_image_x(a))/2-1,
		     (sp_image_y(res)+sp_image_y(a))/2-1,
		     (sp_image_z(res)+sp_image_z(a))/2-1);

  sp_image_free(res);
  res = tmp;*/

/*  sp_image_dephase(res);*/
  /* should be all real */
  for(i = 0;i<sp_image_size(res);i++){
    res->image->data[i] = sp_cscale(res->image->data[i],1.0/sp_image_size(res));
  }

  return res;  
}



real r_factor(Image * fobs, Image *fcalc, real low_intensity_cutoff){
  real den = 0;
  real num = 0;
  int i;
  for(i = 0;i<sp_image_size(fobs);i++){
    if(!fobs->mask->data[i] || sp_real(fobs->image->data[i]) < low_intensity_cutoff){
      continue;
    }
    num +=  sp_cabs(sp_csub(fobs->image->data[i],fcalc->image->data[i]));
    den += sp_cabs(fobs->image->data[i]);
  }
  return num/den;
}





Image * zero_pad_image(Image * a, int newx, int newy, int newz, int pad_mask){
  if(a->shifted){
    return zero_pad_shifted_image(a,newx,newy,newz,pad_mask);
  }else{
    return zero_pad_unshifted_image(a,newx,newy,newz,pad_mask);
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
Image * zero_pad_shifted_image(Image * a, int newx, int newy, int newz,int pad_mask){
  Image * out;
  int x,y,z;
  int sourcex;
  int sourcey;
  int sourcez;
  if(newx < sp_c3matrix_x(a->image) || 
     newy < sp_c3matrix_y(a->image) ||
     newz < sp_c3matrix_z(a->image)){
    fprintf(stderr,"Negative padding!\n");
    printf("shifted (%i,%i,%i) - (%i,%i,%i)\n",newx,newz,newy,sp_c3matrix_x(a->image),sp_c3matrix_y(a->image),sp_c3matrix_z(a->image));
    abort();
  }else if(newx == sp_c3matrix_x(a->image) &&
	   newy == sp_c3matrix_y(a->image) &&
	   newz == sp_c3matrix_z(a->image)){
    return sp_image_duplicate(a,SP_COPY_DATA|SP_COPY_MASK);
  }
  out = sp_image_duplicate(a,SP_COPY_DETECTOR);
  sp_c3matrix_free(out->image);
  out->image = sp_c3matrix_alloc(newx,newy,newz);
  sp_i3matrix_free(out->mask);
  out->mask = sp_i3matrix_alloc(newx,newy,newz);

  for(x = 0;x<sp_c3matrix_x(out->image);x++){
    if(x < sp_c3matrix_x(a->image)/2.0){
      sourcex = x;
    }else if(sp_c3matrix_x(out->image)-x-1 < (sp_c3matrix_x(a->image)-1)/2.0){
      sourcex = (sp_c3matrix_x(a->image)-1)-(sp_c3matrix_x(out->image)-x-1);
    }else{
      sourcex = -1;
    }

    for(y = 0;y<sp_c3matrix_y(out->image);y++){
      if(y < sp_c3matrix_y(a->image)/2.0){
	sourcey = y;
      }else if(sp_c3matrix_y(out->image)-y-1 < (sp_c3matrix_y(a->image)-1)/2.0){
	sourcey = (sp_c3matrix_y(a->image)-1)-(sp_c3matrix_y(out->image)-y-1);
      }else{
	sourcey = -1;
      }
      for(z = 0; z<sp_c3matrix_z(out->image);z++){
	if(z < sp_c3matrix_z(a->image)/2.0){
	  sourcez = z;
	}else if(sp_c3matrix_z(out->image)-z-1 <(sp_c3matrix_z(a->image)-1)/2.0){
	  sourcez = (sp_c3matrix_z(a->image)-1)-(sp_c3matrix_z(out->image)-z-1);
	}else{
	  sourcez = -1;
	}
	if(sourcex == -1 || sourcey == -1 || sourcez == -1){
	  out->image->data[z*sp_c3matrix_y(out->image)*sp_c3matrix_x(out->image)+
			   y*sp_c3matrix_x(out->image)+x] = sp_cinit(0,0);
	  out->mask->data[z*sp_c3matrix_y(out->image)*sp_c3matrix_x(out->image)+
			  y*sp_c3matrix_x(out->image)+x] = pad_mask;
	}else{
	  out->image->data[z*sp_c3matrix_y(out->image)*sp_c3matrix_x(out->image)+y*sp_c3matrix_x(out->image)+x] = a->image->data[sourcez*sp_c3matrix_y(a->image)*sp_c3matrix_x(a->image)+sourcey*sp_c3matrix_x(a->image)+sourcex];
	  out->mask->data[z*sp_c3matrix_y(out->image)*sp_c3matrix_x(out->image)+y*sp_c3matrix_x(out->image)+x] = a->mask->data[sourcez*sp_c3matrix_y(a->image)*sp_c3matrix_x(a->image)+sourcey*sp_c3matrix_x(a->image)+sourcex];
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
Image * zero_pad_unshifted_image(Image * a, int newx, int newy, int newz, int pad_mask){
  Image * out;
  int x,y,z;
  int sourcex;
  int sourcey;
  int sourcez;
  if(newx < sp_c3matrix_x(a->image) || 
     newy < sp_c3matrix_y(a->image) ||
     newz < sp_c3matrix_z(a->image)){
    fprintf(stderr,"Negative padding!\n");
    printf("unshifted (%i,%i,%i) - (%i,%i,%i)\n",newx,newz,newy,sp_c3matrix_x(a->image),sp_c3matrix_y(a->image),sp_c3matrix_z(a->image));
    abort();
  }else if(newx == sp_c3matrix_x(a->image) &&
	   newy == sp_c3matrix_y(a->image) &&
	   newz == sp_c3matrix_z(a->image)){
    return sp_image_duplicate(a,SP_COPY_DATA|SP_COPY_MASK);
  }
  out = sp_image_duplicate(a,SP_COPY_DETECTOR);
  sp_c3matrix_free(out->image);
  out->image = sp_c3matrix_alloc(newx,newy,newz);
  sp_i3matrix_free(out->mask);
  out->mask = sp_i3matrix_alloc(newx,newy,newz);

  for(x = 0;x<sp_c3matrix_x(out->image);x++){
    if(x < sp_c3matrix_x(a->image)){
      sourcex = x;
    }else{
      sourcex = -1;
    }
    
    for(y = 0;y<sp_c3matrix_y(out->image);y++){
      if(y < sp_c3matrix_y(a->image)){
	sourcey = y;
      }else{
	sourcey = -1;
      }

      for(z = 0;z<sp_c3matrix_z(out->image);z++){
	if(z < sp_c3matrix_z(a->image)){
	  sourcez = z;
	}else{
	  sourcez = -1;
	}
	if(sourcey == -1 || sourcex == -1 || sourcez == -1){
	out->image->data[z*sp_c3matrix_y(out->image)*sp_c3matrix_x(out->image)+
			 y*sp_c3matrix_x(out->image)+x] = sp_cinit(0,0);
	out->mask->data[z*sp_c3matrix_y(out->image)*sp_c3matrix_x(out->image)+
			y*sp_c3matrix_x(out->image)+x] = pad_mask;
	}else{
	  out->image->data[z*sp_c3matrix_y(out->image)*sp_c3matrix_x(out->image)+y*sp_c3matrix_x(out->image)+x] = a->image->data[sourcez*sp_c3matrix_y(a->image)*sp_c3matrix_x(a->image)+sourcey*sp_c3matrix_x(a->image)+sourcex];
	  out->mask->data[z*sp_c3matrix_y(out->image)*sp_c3matrix_x(out->image)+y*sp_c3matrix_x(out->image)+x] = a->mask->data[sourcez*sp_c3matrix_y(a->image)*sp_c3matrix_x(a->image)+sourcey*sp_c3matrix_x(a->image)+sourcex];
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
  Image * out = sp_image_duplicate(a,SP_COPY_DATA|SP_COPY_MASK);
  int x,y,z;
  int destx;
  int desty;
  int destz;
  for(x = 0;x<sp_c3matrix_x(out->image);x++){
    destx = (int)(x+sp_c3matrix_x(out->image)-out->detector->image_center[0])%sp_c3matrix_x(out->image);
    for(y = 0;y<sp_c3matrix_y(out->image);y++){
      desty = (int)(y+sp_c3matrix_y(out->image)-out->detector->image_center[1])%sp_c3matrix_y(out->image);
      for(z = 0;z<sp_c3matrix_z(out->image);z++){
	destz = (int)(z+sp_c3matrix_z(out->image)-out->detector->image_center[2])%sp_c3matrix_z(out->image);
	out->image->data[destz*sp_c3matrix_y(out->image)*sp_c3matrix_x(out->image)+desty*sp_c3matrix_x(out->image)+destx] = a->image->data[z*sp_c3matrix_y(a->image)*sp_c3matrix_x(a->image)+y*sp_c3matrix_x(a->image)+x];
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
  for(i = 0;i<sp_image_size(a);i++){
    suma += sp_cabs(sp_csub(a->image->data[i], b->image->data[i]));
    if(sp_cabs(b->image->data[i])){
      sumb += sp_cabs(b->image->data[i]) * log(sp_cabs(b->image->data[i])/(sp_cabs(a->image->data[i])+FLT_EPSILON));
    }
  }
  return suma+sumb;
}


real integrated_intensity(Image * a){
  double sum = 0;
  int i;
  for(i = 0;i<sp_image_size(a);i++){
    sum += sp_cabs(a->image->data[i]);
  }
  return sum;
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
  sp_c3matrix_free(cropped->image);
  cropped->image = sp_c3matrix_alloc(x2-x1+1,y2-y1+1,sp_c3matrix_z(in->image));
  sp_i3matrix_free(cropped->mask);
  cropped->mask = sp_i3matrix_alloc(x2-x1+1,y2-y1+1,sp_c3matrix_z(in->image));

  for(i = y1;i<= y2;i++){
    memcpy(&cropped->image->data[(i-y1)*sp_c3matrix_x(cropped->image)],&in->image->data[(i)*sp_c3matrix_x(in->image)+x1],sp_c3matrix_x(cropped->image)*sizeof(Complex));
    memcpy(&cropped->mask->data[(i-y1)*sp_c3matrix_x(cropped->image)],&in->mask->data[(i)*sp_c3matrix_x(in->image)+x1],sp_c3matrix_x(cropped->image)*sizeof(int));
  }
  return cropped;
}

/* Crops the image so that both the plane x1 and x2 is in the cube.
 * For example cube_crop(a,1,2,3,1,2,3) will return the spot (1,2,3)
 */
Image * cube_crop(Image * in, int x1, int y1, int z1, int x2, int y2, int z2){
  Image * cropped;
  int i,j;
  /* x1,y1 should be upper left, x2,y2 lower right */
  if(x1 > x2 || y1 > y2 || z1 > z2){
    return NULL;
  }
  cropped = sp_image_duplicate(in,SP_COPY_DETECTOR);
  cropped->detector->image_center[0] -= x1;
  cropped->detector->image_center[1] -= y1;
  cropped->detector->image_center[2] -= z1;
  sp_c3matrix_free(cropped->image);
  cropped->image = sp_c3matrix_alloc(x2-x1+1,y2-y1+1,z2-z1+1);
  sp_i3matrix_free(cropped->mask);
  cropped->mask = sp_i3matrix_alloc(x2-x1+1,y2-y1+1,z2-z1+1);


  for(i = z1;i<=z2;i++){
    for(j = y1;j<=y2;j++){
      memcpy(&cropped->image->data[(i-z1)*sp_c3matrix_y(cropped->image)*sp_c3matrix_x(cropped->image)+(j-y1)*sp_c3matrix_x(cropped->image)],
	     &in->image->data[(i)*sp_c3matrix_y(in->image)*sp_c3matrix_x(in->image)+(j)*sp_c3matrix_x(in->image)+x1],sp_c3matrix_x(cropped->image)*2*sizeof(real));
      memcpy(&cropped->mask->data[(i-z1)*sp_c3matrix_y(cropped->image)*sp_c3matrix_x(cropped->image)+(j-y1)*sp_c3matrix_x(cropped->image)],
	     &in->mask->data[(i)*sp_c3matrix_y(in->image)*sp_c3matrix_x(in->image)+(j)*sp_c3matrix_x(in->image)+x1],sp_c3matrix_x(cropped->image)*sizeof(int));
    }
  }
  return cropped;
}


void sp_find_center(Image * img, real * center_x, real * center_y, real * center_z){
  int x,y,z;
  float bx = -1;
  float by = -1;
  float bz = -1;
  Image * a = sp_image_convolute(img,img,NULL);
  real max = 0;
  long long index;
  max = sp_c3matrix_max(a->image,&index);
  sp_c3matrix_get_xyz(a->image,index,&x,&y,&z);
  bx = x;
  by = y;
  bz = z;
  if(bx < sp_c3matrix_x(img->image)/2.0){
    bx = (sp_c3matrix_x(img->image))/2.0+bx/2.0;
  }else{
    bx = (sp_c3matrix_x(img->image))/2.0-(sp_c3matrix_x(img->image)-bx)/2.0;
  }
  if(by < sp_c3matrix_y(img->image)/2.0){
    by = (sp_c3matrix_y(img->image))/2.0+by/2.0;
  }else{
    by = (sp_c3matrix_y(img->image))/2.0-(sp_c3matrix_y(img->image)-by)/2.0;
  }
  if(bz < sp_c3matrix_z(img->image)/2.0){
    bz = (sp_c3matrix_z(img->image))/2.0+bz/2.0;
  }else{
    bz = (sp_c3matrix_z(img->image))/2.0-(sp_c3matrix_z(img->image)-bz)/2.0;
  }
  printf("Center x - %f y - %f z - %f\n",bx,by,bz);
  if(a->num_dimensions == SP_2D){
    //    sp_image_write(a,"centrosym_convolve.png",SpColormapJet|SpColormapLogScale);
  }
  *center_x = bx;
  *center_y = by;
  *center_z = bz;
  sp_image_free(a);
}

long long pixel_to_index(Image * img, real * point){
  return ((int)point[2])*sp_c3matrix_y(img->image)*sp_c3matrix_x(img->image)+
    (int)point[1]*sp_c3matrix_x(img->image)+(int)point[0];
}


/* This doesn't really rescale the mask which is a problem and doesn't really downscale */
Image * fourier_rescale(Image * img, int new_x, int new_y, int new_z){
  Image * res = sp_image_fft(img);
  Image * tmp;
  int i;
  real inv_size;
  if(new_x < sp_c3matrix_x(img->image) ||
     new_y < sp_c3matrix_y(img->image) ||
     new_z < sp_c3matrix_z(img->image)){
    perror("fourier_scale doesn't downscale");
    abort();
  }
  tmp = zero_pad_image(res,new_x,new_y,new_z,1);
  sp_image_free(res);
  res = sp_image_ifft(tmp);
  sp_image_free(tmp);
  res->detector->image_center[0] = img->detector->image_center[0]*(new_x/sp_c3matrix_x(img->image));
  res->detector->image_center[1] = img->detector->image_center[1]*(new_y/sp_c3matrix_y(img->image));
  res->detector->image_center[2] = img->detector->image_center[2]*(new_z/sp_c3matrix_z(img->image));
  inv_size = 1.0/sp_image_size(img);
  for(i = 0;i<sp_image_size(res);i++){
    res->image->data[i] = sp_cscale(res->image->data[i],inv_size);
  }
  return res;
}


Image * bilinear_rescale(Image * img, int new_x, int new_y, int new_z){
  Image * res = sp_image_duplicate(img,SP_COPY_DATA|SP_COPY_MASK);

  real virtual_x;
  real virtual_y;
  real virtual_z;
  int x,y,z;
  res->detector->image_center[0] *= new_x/(real)sp_c3matrix_x(img->image);
  res->detector->image_center[1] *= new_y/(real)sp_c3matrix_y(img->image);
  res->detector->image_center[2] *= new_z/(real)sp_c3matrix_z(img->image);
  res->detector->pixel_size[0] *= sp_c3matrix_x(img->image)/(real)new_x;
  res->detector->pixel_size[1] *= sp_c3matrix_y(img->image)/(real)new_y;
  res->detector->pixel_size[2] *= sp_c3matrix_z(img->image)/(real)new_z;
  sp_c3matrix_free(res->image);
  res->image = sp_c3matrix_alloc(new_x,new_y,new_z);
  sp_i3matrix_free(res->mask);
  res->mask = sp_i3matrix_alloc(new_x,new_y,new_z);
  

  for(x = 0; x<sp_image_x(res); x++){
    virtual_x = (real)x*(sp_image_x(img)-1)/sp_image_x(res);
    for(y = 0; y<sp_image_y(res); y++){
      virtual_y = (real)y*(sp_image_y(img)-1)/sp_image_y(res);
      for(z = 0; z<sp_image_z(res); z++){
	virtual_z = (real)z*(sp_image_z(img)-1)/sp_image_z(res);
	sp_image_set(res,x,y,z,sp_c3matrix_interp(img->image,virtual_x,virtual_y,virtual_z));
	//	sp_i3matrix_set(res->mask,x,y,z,round(sp_i3matrix_interp(img->mask,virtual_x,virtual_y,virtual_z)));
	sp_i3matrix_set(res->mask,x,y,z,round(sp_i3matrix_nearest_neighbour_interp(img->mask,virtual_x,virtual_y,virtual_z)));
      }
    }
  }
  return res;
}
  

real sp_image_interp(const Image * img, real v_x, real v_y, real v_z){
  return sp_real(sp_c3matrix_interp(img->image,v_x,v_y,v_z));
}

void sp_image_scale(Image * img, real value){
  Complex tmp = {value,0};
  sp_c3matrix_scale(img->image,tmp);
}


real sp_image_max(Image * img, long long * index,int * x, int * y, int * z){
  real ret;
  float fx,fy,fz;
  ret = sp_c3matrix_max(img->image,index);
  if(index){
    sp_image_get_coords_from_index(img,*index,&fx,&fy,&fz,SpTopLeftCorner);
    if(x){
      *x = (int)round(fx);
      /*      *x = *index%(sp_image_z(img)*sp_image_y(img));*/
    }
    if(y){
      *y = (int)round(fy);
    }
    if(z){
      *z = (int)round(fz);
      /*      *z = *index%sp_image_x(img)%sp_image_y(img);*/
    }
  }
  return ret;
}


void _sp_image_realloc(Image * img, int new_x, int new_y, int new_z,const char * file, int line){
  _sp_c3matrix_realloc(img->image,new_x,new_y,new_z,file,line);
  _sp_i3matrix_realloc(img->mask,new_x,new_y,new_z,file,line);
}



void sp_image_normalize(Image * in){
  double integral = 0;
  int i;
  for(i = 0;i<sp_image_size(in);i++){
    integral += sp_cabs(in->image->data[i]);
  }
  for(i = 0;i<sp_image_size(in);i++){
    in->image->data[i] = sp_cscale(in->image->data[i],1.0/integral);
  }
}


Image * sp_image_get_mask(const Image * a){
  Image * res = sp_image_duplicate(a,SP_COPY_DATA|SP_COPY_MASK);
  int i;
  for(i = 0;i<sp_image_size(a);i++){
    res->image->data[i] = sp_cinit(res->mask->data[i],0);
  }
  return res;
}


int sp_image_shift_index(Image * a, long long index){
  int x,y,z;
  real fx,fy,fz;
  sp_image_get_coords_from_index(a,index,&fx,&fy,&fz,SpTopLeftCorner);
  x = (int)round(fx);
  y = (int)round(fy);
  z = (int)round(fz);
  /*  int x = index%sp_c3matrix_z(a->image)%sp_c3matrix_y(a->image);
  int y = index/sp_c3matrix_x(a->image)%sp_c3matrix_z(a->image);
  int z = index/sp_c3matrix_x(a->image)/sp_c3matrix_y(a->image);*/
  if(a->shifted){
    x = (x+sp_c3matrix_x(a->image)/2)%sp_c3matrix_x(a->image);
    y = (y+sp_c3matrix_y(a->image)/2)%sp_c3matrix_y(a->image);
    z = (z+sp_c3matrix_z(a->image)/2)%sp_c3matrix_z(a->image);
  }else{
    x -= a->detector->image_center[0];
    y -= a->detector->image_center[1];
    z -= a->detector->image_center[2];
    x += sp_c3matrix_x(a->image);
    x %= sp_c3matrix_x(a->image);
    y += sp_c3matrix_y(a->image);
    y %= sp_c3matrix_y(a->image);
    z += sp_c3matrix_z(a->image);
    z %= sp_c3matrix_z(a->image);
  }
  long long ret = z*sp_c3matrix_y(a->image)*sp_c3matrix_x(a->image)+
    y*sp_c3matrix_x(a->image)+x;
  return ret;

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
  dim[0] = sp_c3matrix_x(img->image);
  dim[1] = sp_c3matrix_y(img->image);
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
  Image  *ret = sp_image_alloc(samples,1,1);
  real fpixel[2];
  real d_to_border;
  real * center;
  int dim[2];
  if(point){
    center = point;
  }else{
    center = img->detector->image_center;
  }
  dim[0] = sp_c3matrix_x(img->image);
  dim[1] = sp_c3matrix_y(img->image);

  d_to_border = sp_image_distance_to_edge(img,center,direction, intersection);
  for(i = 0;i<samples;i++){
    fpixel[0] = center[0]+cos(direction)*d_to_border*((real)i/samples);
    /* The - happens because the pixel 0,0 on an image is the upper left not the lower left */
    fpixel[1] = center[1]-sin(direction)*d_to_border*((real)i/samples);
    /* bilinear interpolation around fpixel */
    ret->image->data[i] = sp_c3matrix_interp(img->image,fpixel[1],fpixel[0],0);
    /* bilinear interpolation around fpixel */
    ret->mask->data[i] = sp_i3matrix_interp(img->mask,fpixel[1],fpixel[0],0);
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
  Image * ret = sp_image_alloc(img_size[0],img_size[1],1);
  int bin;
  real d;
  real u;
  for(x = 0;x<img_size[0];x++){
    for(y = 0;y<img_size[1];y++){
      /* calculate distance to the center and sample from sector accordingly */
      d = sqrt((x-center[0])*(x-center[0])+(y-center[1])*(y-center[1]));
      bin = d;
      u = d-bin;
      if(d+1 < sp_c3matrix_x(sector->image)-1){
	ret->image->data[y*img_size[0]+x] = sp_cscale(sector->image->data[bin],(1.0-u));
	ret->image->data[y*img_size[0]+x] = sp_cadd(ret->image->data[y*img_size[0]+x],sp_cscale(sector->image->data[bin+1],u));
	ret->mask->data[y*img_size[0]+x] = 1;
      }else{
	ret->image->data[y*img_size[0]+x] = sp_cinit(0,0);
	ret->mask->data[y*img_size[0]+x] = 0;
      }
    }
  }
  return ret;
}





Complex sp_image_dot_prod(Image * a, Image * b){
  return sp_c3matrix_froenius_prod(a->image,b->image);
}

Image * sp_proj_module(Image * a, Image * b,SpPlace place){
  Image * ret;
  if(place & SpInPlace){
    ret = a;
  }else{
    ret = sp_image_duplicate(a,SP_COPY_ALL);
  }
  int i;
  for(i = 0;i<sp_image_size(a);i++){
    if(b->mask->data[i]){
      if(sp_cabs(a->image->data[i])){
	  ret->image->data[i] = sp_cscale(ret->image->data[i],sp_cabs(b->image->data[i])/sp_cabs(a->image->data[i]));
      }else{
	ret->image->data[i] = b->image->data[i];
      }
    }else{
      ret->image->data[i] = a->image->data[i];
    }
  }
  return ret;
}


Image * sp_proj_support(Image * a, Image * b, SpPlace place){
  Image * ret;
  if(place & SpInPlace){
    ret = a;
  }else{
    ret = sp_image_duplicate(a,SP_COPY_ALL);
  }
  int i;
  for(i = 0;i<sp_image_size(a);i++){
    if(!sp_cabs(b->image->data[i])){
      ret->image->data[i] = sp_cinit(0,0);
    }
  }
  return ret;
}



typedef struct{
  real value;
  int index;    
}ValueIndexPair;

static int compare_ValueIndexPair(const void * a, const void * b){
  if((*(ValueIndexPair *)a).value < (*(ValueIndexPair *)b).value){
    return -1;
  }else if((*(ValueIndexPair *)a).value == (*(ValueIndexPair *)b).value){
    return 0;
  }else{
    return 1;
  }  
}

Image * sp_proj_module_histogram(Image * a, Image * exp, Image * std_dev){
  ValueIndexPair *  norm_sq_int = sp_malloc(sizeof(ValueIndexPair)*sp_image_size(a));
  Image * ret = sp_image_duplicate(a,SP_COPY_DATA|SP_COPY_MASK);
  int valid_points = 0;
  for(int i = 0;i<sp_image_size(a);i++){
    if(exp->mask->data[i]){
      norm_sq_int[valid_points].value = sp_cabs2(a->image->data[i])-sp_cabs2(exp->image->data[i]);
      norm_sq_int[valid_points].value /= sp_real(std_dev->image->data[i]);
      norm_sq_int[valid_points].index = i;
      valid_points++;
    }
  }
  qsort(norm_sq_int,valid_points,sizeof(ValueIndexPair),compare_ValueIndexPair);
  int vpi = 0;
  for(int i = 0;i<sp_image_size(a);i++){
    if(exp->mask->data[i]){
      int index = norm_sq_int[vpi].index;
      double sigma = sp_real(std_dev->image->data[index]);      
      double new_int = embedded_gsl_cdf_gaussian_Pinv((vpi+1.0)/(valid_points+1.0),sigma)+sp_cabs2(exp->image->data[index]);
      /* If the standard deviation is so big that the projected intensity would be negative we make it 0 */
      if(new_int < 0){
	new_int = 0;
      }
      /* We have to guard against having 0 amplitude*/
      if(sp_cabs(a->image->data[index])){
	ret->image->data[index] = sp_cscale(ret->image->data[index],sqrt(new_int)/sp_cabs(a->image->data[index]));
      }else{
	ret->image->data[index] = sp_cinit(sqrt(new_int),0);
      }
      vpi++;
    }else{
      ret->image->data[i] = a->image->data[i];
    }
  }      
  sp_free(norm_sq_int);
  return ret;
}

int sp_image_invert(Image * a){
  int i;
  for(i = 0;i<sp_image_size(a);i++){
    sp_real(a->image->data[i]) = 1.0/sp_real(a->image->data[i]);
  }
  return 0;
}


void sp_image_mul_elements(Image * a, const Image * b){
  sp_c3matrix_mul_elements(a->image,b->image);  
}

void sp_image_conj(Image * a){  
  sp_c3matrix_conj(a->image);
}


void sp_image_memcpy(Image * dst,const Image * src){  
  sp_i3matrix_memcpy(dst->mask,src->mask);
  sp_c3matrix_memcpy(dst->image,src->image);
  memcpy(dst->detector,src->detector,sizeof(Detector));
}

void sp_image_transpose(Image * a){
  sp_c3matrix_transpose_xy(a->image);
  sp_i3matrix_transpose_xy(a->mask);
}

/*! Inserts image from into image to at the position at_x, at_y, at_z
 *
 *  If from doesn't fit in to, the image is silently clipped.
 */
void sp_image_insert(Image * to, Image * from, int at_x, int at_y, int at_z){
 int y,z;
  for(z = 0;z<sp_min(sp_image_z(from),sp_image_z(to)-at_z);z++){
    for(y = 0;y<sp_min(sp_image_y(from),sp_image_y(to)-at_y);y++){
      memcpy(&(to->image->data[sp_image_get_index(to,at_x,y+at_y,z+at_z)]),
	     &from->image->data[sp_image_get_index(from,0,y,z)],
	     sp_min(sp_image_x(from),sp_image_x(to)-at_x)*sizeof(Complex));
    }
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
 * 
 * type is either SP_2D or SP_3D
 */
Image * sp_image_edge_extend(Image * a, int radius, int edge_flags, int type){
  Image * res;
  int x,y,z,x0,y0,z0;
  /* first put there the initial image */
  if(a->num_dimensions == SP_2D){
    res = sp_image_alloc(sp_c3matrix_x(a->image)+radius*2,sp_c3matrix_y(a->image)+radius*2,1);
    memset(res->image->data,0,sp_image_size(res)*sizeof(Complex));
    sp_image_insert(res,a,radius,radius,0);
  }else if(a->num_dimensions == SP_3D){
    res = sp_image_alloc(sp_c3matrix_x(a->image)+radius*2,sp_c3matrix_y(a->image)+radius*2,sp_c3matrix_z(a->image)+radius*2);
    sp_image_insert(res,a,radius,radius,radius);
  }else{
    res = NULL;
    abort();
  }
  /* now fill the edges */
  if(edge_flags == SP_ZERO_PAD_EDGE){
    return res;
  }else if(edge_flags == SP_SYMMETRIC_EDGE){
    /* First do the four/eight corners */
    for(x = 0,x0=x;x<radius;x++){
      for(y = 0,y0=y;y<radius;y++){
	if(a->num_dimensions == SP_2D){
	  sp_image_set(res,x,y,0,sp_image_get(a,radius-x-1,radius-y-1,0));
	}else if(a->num_dimensions == SP_3D){
	  for(z = 0,z0=z;z<radius;z++){
	    sp_image_set(res,x,y,z,sp_image_get(a,radius-x-1,radius-y-1,
						radius-z-1));
	  }
	}
      }
    }
    for(x = radius+sp_image_x(a),x0=x;x<sp_image_x(res);x++){
      for(y = 0,y0=y;y<radius;y++){
	if(a->num_dimensions == SP_2D){
	  sp_image_set(res,x,y,0,sp_image_get(a,sp_image_x(a)-(x-x0)-1,radius-y-1,
					      0));
	}else if(a->num_dimensions == SP_3D){
	  for(z = 0, z0=z;z<radius;z++){
	    sp_image_set(res,x,y,z,sp_image_get(a,sp_image_x(a)-(x-x0)-1,
						radius-y-1,radius-z-1));
	  }
	}
      }
    }
    
    for(x = radius+sp_image_x(a),x0=x;x<sp_image_x(res);x++){
      for(y = radius+sp_image_y(a),y0=y;y<sp_image_y(res);y++){
	if(a->num_dimensions == SP_2D){
	  sp_image_set(res,x,y,0,sp_image_get(a,sp_image_x(a)-(x-x0)-1,sp_image_y(a)-(y-y0)-1,0));
	}else if(a->num_dimensions == SP_3D){
	  for(z = 0,z0=z;z<radius;z++){
	    sp_image_set(res,x,y,z,sp_image_get(a,sp_image_x(a)-(x-x0)-1,sp_image_y(a)-(y-y0)-1,radius-z-1));
	  }
	}
      }
    }
    
    for(x = 0,x0=x;x<radius;x++){
      for(y = radius+sp_image_y(a),y0=y;y<sp_image_y(res);y++){
	if(a->num_dimensions == SP_2D){
	  sp_image_set(res,x,y,0,sp_image_get(a,radius-(x-x0)-1,sp_image_y(a)-(y-y0)-1,0));
	}else if(a->num_dimensions == SP_3D){
	  for(z = 0,z0=z;z<radius;z++){
	    sp_image_set(res,x,y,z,sp_image_get(a,radius-(x-x0)-1,sp_image_y(a)-(y-y0)-1,radius-z-1));
	  }
	}
      }
    }
    if(a->num_dimensions  == SP_3D){
      for(x = 0,x0=x;x<radius;x++){
	for(y = 0,y0=y;y<radius;y++){
	  for(z = radius+sp_image_z(a),z0=z;z<sp_image_z(res);z++){
	    sp_image_set(res,x,y,z,sp_image_get(a,radius-x-1,radius-y-1,
						sp_image_z(a)-(z-z0)-1));
	  }
	}
      }
      for(x = radius+sp_image_x(a),x0=x;x<sp_image_x(res);x++){
	for(y = 0,y0=y;y<radius;y++){
	  for(z = radius+sp_image_z(a),z0=z;z<sp_image_z(res);z++){
	    sp_image_set(res,x,y,z,sp_image_get(a,sp_image_x(a)-(x-x0)-1,
						radius-y-1,sp_image_z(a)-(z-z0)-1));
	  }
	}
      }
      for(x = radius+sp_image_x(a),x0=x;x<sp_image_x(res);x++){
	for(y = radius+sp_image_y(a),y0=y;y<sp_image_y(res);y++){
	  for(z = radius+sp_image_z(a),z0=z;z<sp_image_z(res);z++){
	    sp_image_set(res,x,y,z,sp_image_get(a,sp_image_x(a)-(x-x0)-1,sp_image_y(a)-(y-y0)-1,sp_image_z(a)-(z-z0)-1));
	  }
	}
      }
      for(x = 0,x0=x;x<radius;x++){
	for(y = radius+sp_image_y(a),y0=y;y<sp_image_y(res);y++){
	  for(z = radius+sp_image_z(a),z0=z;z<sp_image_z(res);z++){
	    sp_image_set(res,x,y,z,sp_image_get(a,radius-(x-x0)-1,sp_image_y(a)-(y-y0)-1,sp_image_z(a)-(z-z0)-1));
	  }
	}
      }
    }
    /* And now the four/six sides */
    for(x = radius,x0=x;x<radius+sp_image_x(a);x++){
      for(y = 0,y0=y;y<radius;y++){
	if(a->num_dimensions == SP_2D){
	  sp_image_set(res,x,y,0,sp_image_get(a,(x-x0),radius-(y-y0)-1,0));
	}else if(a->num_dimensions == SP_3D){
	  for(z = radius,z0=z;z<radius+sp_image_z(a);z++){
	    sp_image_set(res,x,y,z,sp_image_get(a,(x-x0),radius-(y-y0)-1,(z-z0)));
	  }
	}
      }
    }
    for(x = radius,x0=x;x<radius+sp_image_x(a);x++){
      for(y = radius+sp_image_y(a),y0=y;y<sp_image_y(res);y++){
	if(a->num_dimensions == SP_2D){
	  sp_image_set(res,x,y,0,sp_image_get(a,(x-x0),sp_image_y(a)-(y-y0)-1,0));
	}else if(a->num_dimensions == SP_3D){
	  for(z = radius,z0=z;z<radius+sp_image_z(a);z++){
	    sp_image_set(res,x,y,z,sp_image_get(a,(x-x0),sp_image_y(a)-(y-y0)-1,(z-z0)));
	  }
	}
      }
    }
    for(x = 0,x0=x;x<radius;x++){
      for(y = radius,y0=y;y<radius+sp_image_y(a);y++){
	if(a->num_dimensions == SP_2D){
	  sp_image_set(res,x,y,0,sp_image_get(a,radius-(x-x0)-1,(y-y0),0));
	}else if(a->num_dimensions == SP_3D){
	  for(z = radius,z0=z;z<radius+sp_image_z(a);z++){
	    sp_image_set(res,x,y,z,sp_image_get(a,radius-(x-x0)-1,(y-y0),(z-z0)));
	  }
	}
      }
    }
    for(x = radius+sp_image_x(a),x0=x;x<sp_image_x(a);x++){
      for(y = radius,y0=y;y<radius+sp_image_y(a);y++){
	if(a->num_dimensions == SP_2D){
	  sp_image_set(res,x,y,0,sp_image_get(a,sp_image_x(a)-(x-x0)-1,(y-y0),0));
	}else if(a->num_dimensions == SP_3D){
	  for(z = radius,z0=z;z<radius+sp_image_z(a);z++){
	    sp_image_set(res,x,y,z,sp_image_get(a,sp_image_x(a)-(x-x0)-1,(y-y0),(z-z0)));
	  }
	}
      }
    }
    if(a->num_dimensions == SP_2D){
      for(x = radius,x0=x;x<radius+sp_image_x(a);x++){
	for(y = radius,y0=y;y<radius+sp_image_y(a);y++){
	  for(z = 0,z0=z;z<radius;z++){
	    sp_image_set(res,x,y,z,sp_image_get(a,(x-x0),(y-y0),radius-(z-z0)-1));
	  }
	}
      }
      for(x = radius,x0=x;x<radius+sp_image_x(a);x++){
	for(y = radius,y0=y;y<radius+sp_image_y(a);y++){
	  for(z = radius+sp_image_z(a),z0=z;z<sp_image_z(a);z++){
	    sp_image_set(res,x,y,z,sp_image_get(a,(x-x0),(y-y0),sp_image_z(a)-(z-z0)));
	  }
	}
      }
    }
  }else if(edge_flags == SP_CIRCULAR_EDGE){

    /* Lets make this a bit more simple */
    for(x = 0;x<sp_image_x(res);x++){
      /* calculate x on the original image */
      x0 = (x + sp_image_x(a)-(radius%sp_image_x(a)))%sp_image_x(a);
      for(y = 0;y<sp_image_y(res);y++){
	y0 = (y + sp_image_y(a)-(radius%sp_image_y(a)))%sp_image_y(a);
	if(a->num_dimensions == SP_2D){
	  sp_image_set(res,x,y,0,sp_image_get(a,x0,y0,0));
	}else if(a->num_dimensions == SP_3D){
	  for(z = 0;x<sp_image_z(res);x++){
	    z0 = (z + sp_image_z(a)-(radius%sp_image_z(a)))%sp_image_z(a);
	    sp_image_set(res,x,y,z,sp_image_get(a,x0,y0,z0));
	  }
	}
      }
    }
  }else if(edge_flags == SP_REPLICATE_EDGE){
    /* First do the four/eight corners */
    for(x = 0,x0=x;x<radius;x++){
      for(y = 0,y0=y,y0=y;y<radius;y++){
	if(a->num_dimensions == SP_2D){
	  sp_image_set(res,x,y,0,sp_image_get(a,0,0,0));
	}else if(a->num_dimensions == SP_3D){
	  for(z = 0,z0=z;z<radius;z++){
	    sp_image_set(res,x,y,z,sp_image_get(a,0,0,0));
	  }
	}
      }
    }
    for(x = radius+sp_image_x(a),x0=x;x<sp_image_x(res);x++){
      for(y = 0,y0=y;y<radius;y++){
	if(a->num_dimensions == SP_2D){
	sp_image_set(res,x,y,0,sp_image_get(a,sp_image_x(a)-1,0,0));
	}else if(a->num_dimensions == SP_3D){
	  for(z = 0,z0=z;z<radius;z++){
	    sp_image_set(res,x,y,z,sp_image_get(a,sp_image_x(a)-1,0,0));
	  }
	}
      }
    }

    for(x = radius+sp_image_x(a),x0=x;x<sp_image_x(res);x++){
      for(y = radius+sp_image_y(a),y0=y;y<sp_image_y(res);y++){
	if(type == SP_2D){
	sp_image_set(res,x,y,0,sp_image_get(a,sp_image_x(a)-1,sp_image_y(a)-1,0));
	}else{
	  for(z = 0,z0=z;z<radius;z++){
	    sp_image_set(res,x,y,z,sp_image_get(a,sp_image_x(a)-1,sp_image_y(a)-1,0));
	  }
	}
      }
    }

    for(x = 0,x0=x;x<radius;x++){
      for(y = radius+sp_image_y(a),y0=y;y<sp_image_y(res);y++){
	if(type == SP_2D){
	sp_image_set(res,x,y,0,sp_image_get(a,0,sp_image_y(a)-1,0));
	}else{
	  for(z = 0,z0=z;z<radius;z++){
	    sp_image_set(res,x,y,z,sp_image_get(a,0,sp_image_y(a)-1,0));
	  }
	}
      }
    }
    if(type != SP_2D){
      for(x = 0,x0=x;x<radius;x++){
	for(y = 0,y0=y,y0=y;y<radius;y++){
	  for(z = radius+sp_image_z(a),z0=z;z<sp_image_z(res);z++){
	    sp_image_set(res,x,y,z,sp_image_get(a,0,0,sp_image_z(a)-1));
	  }
	}
      }
      for(x = radius+sp_image_x(a),x0=x;x<sp_image_x(res);x++){
	for(y = 0,y0=y;y<radius;y++){
	  for(z = radius+sp_image_z(a),z0=z;z<sp_image_z(res);z++){
	    sp_image_set(res,x,y,z,sp_image_get(a,0,0,sp_image_z(a)-1));
	  }
	}
      }
      
      for(x = radius+sp_image_x(a),x0=x;x<sp_image_x(res);x++){
	for(y = radius+sp_image_y(a),y0=y;y<sp_image_y(res);y++){
	  for(z = radius+sp_image_z(a),z0=z;z<sp_image_z(res);z++){
	    sp_image_set(res,x,y,z,sp_image_get(a,sp_image_x(a)-1,sp_image_y(a)-1,sp_image_z(a)-1));
	  }
	}
      }
      
      for(x = 0,x0=x;x<radius;x++){
	for(y = radius+sp_image_y(a),y0=y;y<sp_image_y(res);y++){
	  for(z = radius+sp_image_z(a),z0=z;z<sp_image_z(res);z++){
	    sp_image_set(res,x,y,z,sp_image_get(a,0,sp_image_y(a)-1,sp_image_z(a)-1));
	  }
	}
      }
    }
    /* And now the four/six sides */
    for(x = radius,x0=x;x<radius+sp_image_x(a);x++){
      for(y = 0,y0=y;y<radius;y++){
	if(a->num_dimensions == SP_2D){
	  sp_image_set(res,x,y,0,sp_image_get(a,(x-x0),0,0));
	}else{
	  for(z = radius,z0=z;z<radius+sp_image_z(a);z++){
	    sp_image_set(res,x,y,z,sp_image_get(a,(x-x0),0,(z-z0)));
	  }
	}
      }
    }
    for(x = radius,x0=x;x<radius+sp_image_x(a);x++){
      for(y = radius+sp_image_y(a),y0=y;y<sp_image_y(res);y++){
	if(a->num_dimensions == SP_2D){
	  sp_image_set(res,x,y,0,sp_image_get(a,(x-x0),sp_image_y(a),0));
	}else{
	  for(z = radius,z0=z;z<radius+sp_image_z(a);z++){
	    sp_image_set(res,x,y,z,sp_image_get(a,(x-x0),sp_image_y(a)-1,(z-z0)));
	  }
	}
      }
    }
    for(x = 0,x0=x;x<radius;x++){
      for(y = radius,y0=y;y<radius+sp_image_y(a);y++){
	if(a->num_dimensions == SP_2D){
	  sp_image_set(res,x,y,0,sp_image_get(a,sp_image_x(a)-radius+(x-x0),(y-y0),0));
	}else{
	  for(z = radius,z0=z;z<radius+sp_image_z(a);z++){
	    sp_image_set(res,x,y,0,sp_image_get(a,0,(y-y0),(z-z0)));
	  }
	}
      }
    }
    for(x = radius+sp_image_x(a),x0=x;x<sp_image_x(res);x++){
      for(y = radius,y0=y;y<radius+sp_image_y(a);y++){
	if(a->num_dimensions == SP_2D){
	  sp_image_set(res,x,y,0,sp_image_get(a,(x-x0),(y-y0),0));
	}else{
	  for(z = radius,z0=z;z<radius+sp_image_z(a);z++){
	    sp_image_set(res,x,y,z,sp_image_get(a,sp_image_x(a)-1,(y-y0),(z-z0)));
	  }
	}
      }
    }
    if(type != SP_2D){
      for(x = radius,x0=x;x<radius+sp_image_x(a);x++){
	for(y = radius,y0=y;y<radius+sp_image_y(a);y++){
	  for(z = 0,z0=z;z<radius;z++){
	    sp_image_set(res,x,y,z,sp_image_get(a,(x-x0),(y-y0),0));
	  }
	}
      }
      for(x = radius,x0=x;x<radius+sp_image_x(a);x++){
	for(y = radius,y0=y;y<radius+sp_image_y(a);y++){
	  for(z = radius+sp_image_z(a),z0=z;z<sp_image_z(res);z++){
	    sp_image_set(res,x,y,z,sp_image_get(a,(x-x0),(y-y0),sp_image_z(a)-1));
	  }
	}
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
void sp_image_median_filter(Image * a,sp_i3matrix * kernel, int edge_flags, int type){
  int integral = 0;
  int n,i,x,y,z;
  int kx,ky,kz;
  int kcx = sp_i3matrix_x(kernel)/2;
  int kcy = sp_i3matrix_y(kernel)/2;
  int kcz = sp_i3matrix_z(kernel)/2;
  int radius;
  Image * work;

  for(i =0 ;i<sp_i3matrix_size(kernel);i++){
    integral += kernel->data[i];
  }
  real * buffer = sp_malloc(sizeof(real)*integral);
  if(a->num_dimensions == SP_2D){
    radius = sp_max((sp_i3matrix_x(kernel)+1)/2,
			(sp_i3matrix_y(kernel)+1)/2);
  }else{
    radius = sp_max((sp_i3matrix_x(kernel)+1)/2,
			sp_max((sp_i3matrix_y(kernel)+1)/2,
			       (sp_i3matrix_z(kernel)+1)/2));
  }
  if(a->num_dimensions == SP_2D){
    work = sp_image_edge_extend(a,radius,edge_flags,SP_2D);
  }else{
    work = sp_image_edge_extend(a,radius,edge_flags,SP_3D);
  }

  if(a->num_dimensions == SP_2D){
    z = 0;
    for(x = radius;x<radius+sp_image_x(a);x++){
      for(y = radius;y<radius+sp_image_y(a);y++){
	n = 0;
	for(kx = -kcx;kx<sp_i3matrix_x(kernel)-kcx;kx++){
	  for(ky = -kcy;ky<sp_i3matrix_y(kernel)-kcy;ky++){
	    for(i = 0;i<sp_i3matrix_get(kernel,ky+kcy,kx+kcx,0);i++){
	      buffer[n++] = sp_cabs(sp_image_get(work,x+kx,y+ky,0));
	    }
	  }
	}
	sp_bubble_sort(buffer,n);
	Complex tmp = sp_cinit(0,0);
	if(n%2){	  
	  /* odd n take the middle one */
	  sp_real(tmp) = buffer[n/2];
	  if(a->num_dimensions == SP_2D){
	    sp_image_set(a,x-radius,y-radius,0,tmp);
	  }else{
	    sp_image_set(a,x-radius,y-radius,z-radius,tmp);
	  }
	}else{
	  /* even n take the average */
	  sp_real(tmp) = (buffer[n/2]+buffer[n/2-1])/2;
	  sp_image_set(a,x-radius,y-radius,0,tmp);
	}
      }
    }
  }else{
    for(x = radius;x<radius+sp_image_x(a);x++){
      for(y = radius;y<radius+sp_image_y(a);y++){
	for(z = radius;z<radius+sp_image_z(a);z++){
	  n = 0;
	  for(kx = -kcx;kx<sp_i3matrix_x(kernel)-kcx;kx++){
	    for(ky = -kcy;ky<sp_i3matrix_y(kernel)-kcy;ky++){
	      for(kz = -kcz;kz<sp_i3matrix_z(kernel)-kcz;kz++){
		for(i = 0;i<sp_i3matrix_get(kernel,ky+kcy,kx+kcx,kz+kcz);i++){
		  buffer[n++] = sp_cabs(sp_image_get(work,x+kx,y+ky,z+kz));
		}
	      }
	    }
	  }
	  sp_bubble_sort(buffer,n);
	  Complex tmp = sp_cinit(0,0);
	  if(n%2){
	    sp_real(tmp) = buffer[n/2];
	    /* odd n take the middle one */
	    if(a->num_dimensions == SP_2D){
	      sp_image_set(a,x-radius,y-radius,0,tmp);
	    }else{
	      sp_image_set(a,x-radius,y-radius,z-radius,tmp);
	    }
	  }else{
	    sp_real(tmp) =(buffer[n/2]+buffer[n/2-1])/2; 
	    /* even n take the average */
	    if(a->num_dimensions == SP_2D){
	      sp_image_set(a,x-radius,y-radius,0,tmp);
	    }else{
	      sp_image_set(a,x-radius,y-radius,z-radius,tmp);
	    }
	  }
	}
      }
    }
  }
  sp_free(buffer);
  sp_image_free(work);
}

/*! Contrast stretches an image
 *
 * The image is partitioned in x_div*y_div points. x_div and y_div must be at least 2
 * Each point is contrast stretched in relation to its neighbourhood
 * The scaling is linear interpolated from division to division
 */
void sp_image_adaptative_constrast_stretch(Image * a,int x_div, int y_div){
  int dx,dy,x,y;
  int x_div_size = sp_image_x(a)/(x_div-1);
  int y_div_size = sp_image_y(a)/(y_div-1);
  sp_3matrix * div_max = sp_3matrix_alloc(x_div,y_div,1);
  sp_3matrix * div_min = sp_3matrix_alloc(x_div,y_div,1);
  sp_3matrix * div_mean = sp_3matrix_alloc(x_div,y_div,1);
  sp_3matrix * div_std_dev = sp_3matrix_alloc(x_div,y_div,1);
  sp_3matrix * offsets = sp_3matrix_alloc(sp_image_x(a),sp_image_y(a),1);
  sp_3matrix * factors = sp_3matrix_alloc(sp_image_x(a),sp_image_y(a),1);
  int radius = sp_max(x_div_size/2,y_div_size/2);
  Image * work = sp_image_edge_extend(a, radius, SP_SYMMETRIC_EDGE,SP_2D);
  for(dx = 0;dx <x_div;dx++){
    for(dy = 0;dy <y_div;dy++){
      sp_3matrix_set(div_max,dx,dy,0,sp_cabs(sp_image_get(work,radius+dx*x_div_size,radius+dy*y_div_size,0)));
      sp_3matrix_set(div_min,dx,dy,0,sp_cabs(sp_image_get(work,radius+dx*x_div_size,radius+dy*y_div_size,0)));
      for(x = radius+dx*x_div_size;x<radius+(dx+1)*x_div_size;x++){	
	if(x == sp_image_x(work)){
	  break;
	}

	for(y = radius+dy*y_div_size;y<radius+(dy+1)*y_div_size;y++){
	  if(y == sp_image_y(work)){
	    break;
	  }
	  if(sp_cabs(sp_image_get(work,x,y,0)) < sp_3matrix_get(div_min,dx,dy,0)){
	    sp_3matrix_set(div_min,dx,dy,0,sp_cabs(sp_image_get(work,x,y,0)));
	  }
	  if(sp_cabs(sp_image_get(work,x,y,0)) > sp_3matrix_get(div_max,dx,dy,0)){
	    sp_3matrix_set(div_max,dx,dy,0,sp_cabs(sp_image_get(work,x,y,0)));
	  }
	  sp_3matrix_set(div_mean,dx,dy,0,sp_cabs(sp_image_get(work,x,y,0))+sp_3matrix_get(div_mean,dx,dy,0));
	}
      }
      sp_3matrix_set(div_mean,dx,dy,0,sp_3matrix_get(div_mean,dx,dy,0)/(x_div_size*y_div_size));
      /* Second pass to calculate standard deviation */
      for(x = radius+dx*x_div_size;x<radius+(dx+1)*x_div_size;x++){	
	if(x == sp_image_x(work)){
	  break;
	}

	for(y = radius+dy*y_div_size;y<radius+(dy+1)*y_div_size;y++){
	  if(y == sp_image_y(work)){
	    break;
	  }
	  sp_3matrix_set(div_std_dev,dx,dy,0,sp_3matrix_get(div_std_dev,dx,dy,0)+
			fabs(sp_cabs(sp_image_get(work,x,y,0))*sp_cabs(sp_image_get(work,x,y,0))-
			     sp_3matrix_get(div_mean,dx,dy,0)*sp_3matrix_get(div_mean,dx,dy,0)));
	}
      }
      sp_3matrix_set(div_std_dev,dx,dy,0,sqrt(sp_3matrix_get(div_std_dev,dx,dy,0)/((x_div_size*y_div_size)-1)));
    }
  }
  sp_image_free(work);
  for(x = 0;x<sp_image_x(a);x++){
    for(y = 0;y<sp_image_y(a);y++){
      sp_3matrix_set(offsets,x,y,0,-sp_3matrix_interp(div_mean,sp_min((float)x/x_div_size,x_div-0.001),sp_min((float)y/y_div_size,y_div-0.001),0)); 
      /* allow 3 standard deviation on each side of the mean */
      sp_3matrix_set(factors,x,y,0,1.0/(6*sp_3matrix_interp(div_std_dev,sp_min((float)x/x_div_size,x_div-0.001),sp_min((float)y/y_div_size,y_div-0.001),0))); 
    }
  }
  for(x = 0;x<sp_image_x(a);x++){
    for(y = 0;y<sp_image_y(a);y++){
      /* cap extreme values */
      Complex tmp = sp_cinit(0,0);
      if(sp_real(sp_image_get(a,x,y,0)) < -sp_3matrix_get(offsets,x,y,0)-1/sp_3matrix_get(factors,x,y,0)){	
	sp_real(tmp) = -sp_3matrix_get(offsets,x,y,0)-1/sp_3matrix_get(factors,x,y,0);
	sp_image_set(a,x,y,0,tmp);
      }else if(sp_real(sp_image_get(a,x,y,0)) > -sp_3matrix_get(offsets,x,y,0)+1/sp_3matrix_get(factors,x,y,0)){
	sp_real(tmp) = -sp_3matrix_get(offsets,x,y,0)+1/sp_3matrix_get(factors,x,y,0);
	sp_image_set(a,x,y,0,tmp);
      }
      tmp = sp_image_get(a,x,y,0);
      sp_real(tmp) += sp_3matrix_get(offsets,x,y,0);
      tmp = sp_cscale(tmp,sp_3matrix_get(factors,x,y,0));
      sp_image_set(a,x,y,0,tmp);
    }
  }
  sp_3matrix_free(offsets);
  sp_3matrix_free(factors);
  sp_3matrix_free(div_max);
  sp_3matrix_free(div_min);
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
 
 
void sp_image_fourier_coords(Image * in, sp_3matrix * k_x, sp_3matrix * k_y, sp_3matrix * k_z){
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
  real real_to_reciprocal = 1.0/(in->detector->detector_distance*in->detector->wavelength);
  real ewald_radius = 1.0/in->detector->wavelength;
  real distance_to_ewald_sphere_center;

  real det_x = in->detector->pixel_size[0] * sp_image_x(in);
  real det_y = in->detector->pixel_size[1] * sp_image_y(in);
  
  nx = sp_image_x(in);
  ny = sp_image_y(in);

  for(x = 0;x<nx;x++){
    for(y = 0;y<ny;y++){
      /* 
	 Calculate the pixel coordinates in reciprocal space 	 
	 by dividing the physical position by detector_distance*wavelength.
	 
	 CCD center at image_center(nx-1)/2,(ny-1)/2

	 Upper left corner of the detector with negative x and positive y
      */
      px = ((x-in->detector->image_center[0])/nx)*det_x;
      py = ((in->detector->image_center[1]-y)/ny)*det_y;
      
      rx = px*real_to_reciprocal;
      ry = py*real_to_reciprocal;
      /* Project pixel into Ewald sphere. */
      distance_to_ewald_sphere_center = sqrt(rx*rx+ry*ry+ewald_radius*ewald_radius);
      if(k_x){      
	sp_3matrix_set(k_x,x,y,0,rx * ewald_radius/distance_to_ewald_sphere_center);
      }
      if(k_y){
	sp_3matrix_set(k_y,x,y,0,ry * ewald_radius/distance_to_ewald_sphere_center);
      }
      if(k_z){
	sp_3matrix_set(k_z,x,y,0,ewald_radius-(ewald_radius * ewald_radius/distance_to_ewald_sphere_center));
      }
    }
  }
}





/*! Superimposes image b on top of image a 
 *
 *  flags is a bitwise combination of the following:
 *
 *  SP_ENANTIOMORPH - allow to try to superimpose not only b but also
 *  the "mirror image" of b [b(-x)].
 *
*/
void sp_image_superimpose(const Image * _a,Image * _b, SpSuperimposeFlags flags){
  int x,y,z;
  long long index;
  int center_invert = 0;
  real max;
  /* check maximum overlap of a and b */
  Image * a = sp_image_duplicate(_a,SP_COPY_DATA);
  Image * b = sp_image_duplicate(_b,SP_COPY_DATA);
  
  sp_image_dephase(a);
  sp_image_dephase(b);
  Image * direct_overlap = sp_image_cross_correlate(a,b,NULL);
  max = sp_image_max(direct_overlap,&index,&x,&y,&z);
  sp_image_free(direct_overlap);
  if(flags & SpEnantiomorph){
    int x2,y2,z2;
    long long index2;
    Image * enantio_overlap = sp_image_convolute(a,b,NULL);
    real max2 = sp_image_max(enantio_overlap,&index2,&x2,&y2,&z2);
    sp_image_free(enantio_overlap);
    if(max2 > max){
      center_invert = 1;
      max = max2;
      x = x2+1;
      y = y2+1;
      z = z2+1;
      sp_image_reflect(_b,IN_PLACE,SP_ORIGO);
    }
  }
  sp_image_free(a);
  sp_image_free(b);
  sp_image_translate(_b,x,y,z,SP_TRANSLATE_WRAP_AROUND);
}



/*! Superimposes image b on top of image a with fractional pixel precision
 *
 *  flags is a bitwise combination of the following:
 *
 *  SpEnantiomorph - allow to try to superimpose not only b but also
 *  the "mirror image" of b [b(-x)].
 *  SpCorrectPhaseShift - allow to try to superimpose not only b but also
 *  the "mirror image" of b [b(-x)].
 *
 *  A precision==2 corresponds to superpositions with 1/2 pixels precision
 *  precision==3 corresponds to superpositions with 1/3 pixels precision and so forth
 *  The image is padded with zeroes so as to becomes precision*original so the 
 *  run time is proportional to the precision to the power of the image dimension.
 *
 *
*/
void sp_image_superimpose_fractional(Image * _a,Image * _b, SpSuperimposeFlags flags, int precision){
  int x,y,z;
  long long index;
  int center_invert = 0;
  real max;
  /* check maximum overlap of a and b */
  Image * a = sp_image_duplicate(_a,SP_COPY_DATA);
  Image * b = sp_image_duplicate(_b,SP_COPY_DATA);
  real phase_correction = 0;

  Image * direct_overlap = sp_image_convolute_fractional(a,b,NULL,precision,1);
  /* It's important that this is the maximum of the absolute value
   Because in the case of phase shifted input that maximum is gonna have the
  same value in abosolute value as |(a+bi)*((a-bi)*exp(phi))| = a^2+b^2 */
  max = sp_image_max(direct_overlap,&index,&x,&y,&z);
  phase_correction = sp_carg(direct_overlap->image->data[index]);
  sp_image_free(direct_overlap);
  if(flags & SpEnantiomorph){
    int x2,y2,z2;
    long long index2;

    Image * enantio_overlap = sp_image_convolute_fractional(a,b,NULL,precision,0);

    real max2 = sp_image_max(enantio_overlap,&index2,&x2,&y2,&z2);
    if(max2 > max){
      center_invert = 1;
      max = max2;
      x = x2+precision;
      y = y2+precision;
      z = z2+precision;
      phase_correction = sp_carg(enantio_overlap->image->data[index2]);
      /* enantiomorph is both reflected and conjugated */
      sp_image_reflect(_b,IN_PLACE,SP_ORIGO);
      sp_image_conj(_b);
    }
    sp_image_free(enantio_overlap);
  }
  sp_image_free(a);
  sp_image_free(b);
  sp_image_fourier_translate(_b,(real)x/precision,(real)y/precision,(real)z/precision);
  if(flags & SpCorrectPhaseShift){
    sp_image_phase_shift(_b,phase_correction,1);
  }
}


void sp_image_superimpose_fractional2(Image * _a,Image * _b, SpSuperimposeFlags flags, int precision){
  int x,y,z;
  long long index;
  int center_invert = 0;
  real max;
  /* check maximum overlap of a and b */
  Image * a = sp_image_duplicate(_a,SP_COPY_DATA);
  Image * b = sp_image_duplicate(_b,SP_COPY_DATA);
  real phase_correction = 0;

  Image * direct_overlap = sp_image_convolute_fractional(a,b,NULL,1,1);
  /* It's important that this is the maximum of the absolute value
   Because in the case of phase shifted input that maximum is gonna have the
  same value in abosolute value as |(a+bi)*((a-bi)*exp(phi))| = a^2+b^2 */
  max = sp_image_max(direct_overlap,&index,&x,&y,&z);
  phase_correction = sp_carg(direct_overlap->image->data[index]);
  sp_image_free(direct_overlap);
  if(flags & SpEnantiomorph){
    int x2,y2,z2;
    long long index2;

    Image * enantio_overlap = sp_image_convolute_fractional(a,b,NULL,precision,0);

    real max2 = sp_image_max(enantio_overlap,&index2,&x2,&y2,&z2);
    if(max2 > max){
      center_invert = 1;
      max = max2;
      x = x2+precision;
      y = y2+precision;
      z = z2+precision;
      phase_correction = sp_carg(enantio_overlap->image->data[index2]);
      /* enantiomorph is both reflected and conjugated */
      sp_image_reflect(_b,IN_PLACE,SP_ORIGO);
      sp_image_conj(_b);
    }
    sp_image_free(enantio_overlap);
  }
  sp_image_free(a);
  sp_image_free(b);
  sp_image_fourier_translate(_b,(real)x/1,(real)y/1,(real)z/1);
  if(flags & SpCorrectPhaseShift){
    sp_image_phase_shift(_b,phase_correction,1);
  }
  max = 0;
  real  max_dx = 0;
  real  max_dy = 0;
  real  max_dz = 0;
  if(b->num_dimensions != SP_2D){
    for(real dz = -0.5;dz<0.5;dz+= 1.0/precision){
      for(real dy = -0.5;dy<0.5;dy+= 1.0/precision){
	for(real dx = -0.5;dx<0.5;dx+= 1.0/precision){
	  Image * b = sp_image_duplicate(_b,SP_COPY_DATA);
	  sp_image_fourier_translate(b,dx,dy,dz);
	  real sum = sp_cabs(sp_c3matrix_froenius_prod(_a->image,b->image));
	  if(sum > max){
	    max = sum;
	    max_dx = dx;
	    max_dy = dy;
	    max_dz = dz;
	  }
	  sp_image_free(b);
	}
      }
    }
  }else{
    for(real dy = -1+1.0/precision;dy<1;dy+= 1.0/precision){
      for(real dx = -1+1.0/precision;dx<1;dx+= 1.0/precision){
	Image * b = sp_image_duplicate(_b,SP_COPY_DATA);
	sp_image_fourier_translate(b,dx,dy,0);
	real sum = sp_cabs(sp_c3matrix_froenius_prod(_a->image,b->image));
	if(sum > max){
	  max = sum;
	  max_dx = dx;
	  max_dy = dy;
	}
	sp_image_free(b);
      }
    }
  }
  sp_image_fourier_translate(_b,max_dx,max_dy,max_dz);
}


/*! Minimize the difference between the phases of a and b by adding a constant phase to b.
 *
 * The returned value is the phase factor in radians.
 * The method used to minimize the phase difference is to take the average 
 * of the vectors representing the phase difference between a and b.
 * The constant phase is then simply the angle of the average vector.
 * If wieghted is 1 the magnitude of each pixel is used as a weighting
 * for the averaging. If it's 2 then the square of the magnitude is used.
 * The weight is taken from image a. 
 * Both images are assumed to have the same dimensions.
 */
real sp_image_phase_match(Image * a, Image * b,int weighted){
  int i;
  Complex average = sp_cinit(0,0);
  real angle;
  real magnitude;
  Complex v;
  real phi_constant;
  for(i = 0;i<sp_image_size(a);i++){
    angle = sp_carg(a->image->data[i]) - sp_carg(b->image->data[i]);
    if(weighted == 1){
      magnitude = sp_cabs(a->image->data[i]);
    }else if(weighted == 2){
      magnitude = sp_cabs2(a->image->data[i]);
    }else{
      magnitude = 1;
    }
    v = sp_cinit(cos(angle)*magnitude,sin(angle)*magnitude);
    average = sp_cadd(average,v);
  }
  phi_constant = sp_carg(average);
  for(i = 0;i<sp_image_size(b);i++){
    angle = sp_carg(b->image->data[i])+phi_constant;
    magnitude = sp_cabs(b->image->data[i]);
    b->image->data[i] = sp_cinit(cos(angle)*magnitude,sin(angle)*magnitude);
  }
  return phi_constant;
}

void sp_image_translate(Image * a, int x,int y,int z,int flags){
  int image_z = sp_image_z(a);
  int image_y = sp_image_y(a);
  int image_x = sp_image_x(a);
  Image * tmp = sp_image_alloc(image_x,image_y,image_z);
  if(flags & SP_TRANSLATE_WRAP_AROUND){
    for(int pz = 0;pz<image_z;pz++){
      int nz = pz+z;
      while(nz < 0){
	nz += image_z;
      }
      nz = nz % image_z;
      for(int py = 0;py<sp_image_y(a);py++){
	int ny = py+y;
	while(ny < 0){
	  ny += image_y;
	}
	ny = ny % image_y;	
	for(int px = 0;px<sp_image_x(a);px++){
	  int nx = px+x;
	  while(nx < 0){
	    nx += image_x;
	  }
	  nx = nx % image_x;	
	  sp_image_set(tmp,nx,ny,nz,sp_image_get(a,px,py,pz));
	  sp_i3matrix_set(tmp->mask,nx,ny,nz,sp_i3matrix_get(a->mask,px,py,pz));	  	  
	}
      }
    }
  }else if(flags & SP_TRANSLATE_DISCARD_OUTSIDE){
    for(int pz = 0;pz<image_z;pz++){
      int nz = pz+z;
      if(nz < 0 || nz >= image_z){
	continue;
      }
      for(int py = 0;py<sp_image_y(a);py++){
	int ny = py+y;
	if(ny < 0 || ny >= image_y){
	  continue;
	}
	for(int px = 0;px<sp_image_x(a);px++){
	  int nx = px+x;
	  if(nx < 0 || nx >= image_x){
	    continue;
	  }
	  sp_image_set(tmp,nx,ny,nz,sp_image_get(a,px,py,pz));
	  sp_i3matrix_set(tmp->mask,nx,ny,nz,sp_i3matrix_get(a->mask,px,py,pz));	  	  
	}
      }
    }
  }
  sp_c3matrix * tmp2 = a->image;
  sp_i3matrix * tmp3 = a->mask;
  a->image = tmp->image;
  a->mask = tmp->mask;
  tmp->image = tmp2;
  tmp->mask = tmp3;
  sp_image_free(tmp);
}
  

void sp_image_fourier_translate(Image * ra, real t_x, real t_y, real t_z){
  Image * a = sp_image_fft(ra);
  /* fourier frequency x*/
  int f_x;
  /* fourier frequency y*/
  int f_y;
  /* fourier frequency z*/
  int f_z;
  long long i = 0;
  real nx_inv = 1.0/sp_image_x(a);
  real ny_inv = 1.0/sp_image_y(a);
  real nz_inv = 1.0/sp_image_z(a);
  int x,y,z;
  real two_pi = 2*M_PI;
  for(z = 0;z<sp_image_z(a);z++){    
    if(z < sp_image_z(a)/2){
      f_z = z;
    }else{
      f_z = -(sp_image_z(a)-z);
    }
    for(y = 0;y<sp_image_y(a);y++){
      if(y < sp_image_y(a)/2){
	f_y = y;
      }else{
	f_y = -(sp_image_y(a)-y);
      }
      for(x = 0;x<sp_image_x(a);x++){
	if(x < sp_image_x(a)/2){
	  f_x = x;
	}else{
	  f_x = -(sp_image_x(a)-x);
	}
	a->image->data[i] = sp_cmul(a->image->data[i],
				    sp_cinit(cos(f_x * t_x * nx_inv * two_pi+
						 f_y * t_y * ny_inv * two_pi+
						 f_z * t_z * nz_inv * two_pi),
					     -sin(f_x * t_x * nx_inv * two_pi+
						  f_y * t_y * ny_inv * two_pi+
						  f_z * t_z * nz_inv * two_pi)));
	i++;
      }
    }
  }
  Image * tmp = sp_image_ifft(a);
  sp_image_scale(tmp,1.0/sp_image_size(tmp));
  sp_image_memcpy(ra,tmp);
  sp_image_free(a);
  sp_image_free(tmp);
}

real sp_image_rs_r_factor(Image * a, Image * b){
  real sum_dif = 0;
  real sum_sum = 0;
  for(int i = 0; i < sp_image_size(a);i++){
    sum_sum += sp_cabs(a->image->data[i])+sp_cabs(b->image->data[i]);
    sum_dif += fabs(sp_cabs(a->image->data[i])-sp_cabs(b->image->data[i]));
  }
  return sum_dif/sum_sum;
}

real sp_image_correlation_coefficient(const Image * a,const Image * b){
  real sum_sq_x = 0;
  real sum_sq_y = 0;
  real sum_coproduct = 0;
  real mean_x = sp_cabs(a->image->data[0]);
  real mean_y = sp_cabs(b->image->data[0]);
  for(int i= 1; i < sp_image_size(a);i++){
    real sweep = i/(i+1.0);
    real delta_x = sp_cabs(a->image->data[i]) - mean_x;
    real delta_y = sp_cabs(b->image->data[i]) - mean_y;
    sum_sq_x += delta_x * delta_x * sweep;
    sum_sq_y += delta_y * delta_y * sweep;
    sum_coproduct += delta_x * delta_y * sweep;
    mean_x += delta_x / i;
    mean_y += delta_y / i;
  }
  real pop_sd_x = sqrt( sum_sq_x / sp_image_size(a) );
  real  pop_sd_y = sqrt( sum_sq_y / sp_image_size(b) );
  real cov_x_y = sum_coproduct / sp_image_size(a);
  real correlation = cov_x_y / (pop_sd_x * pop_sd_y);
  return correlation;
}

sp_vector * sp_image_center_of_mass(Image * a){
  return sp_c3matrix_center_of_mass(a->image);
}


int sp_image_get_coords_from_index(const Image * in,int index,real * x, real * y, real * z, SpOrigin origin){
  int nx,ny,nz;
  nx = sp_image_x(in);
  ny = sp_image_y(in);
  nz = sp_image_z(in);
  if(origin == SpTopLeftCorner){    
    *z = index/(ny*nx);
    *y = (index%(ny*nx))/nx;
    *x = index%(nx);
    return 0;
  }else if(origin == SpImageCenter){
    *z = index/(ny*nx);
    *y = (index%(ny*nx))/nx;
    *x = index%(nx);
    *z -= in->detector->image_center[2];
    *y -= in->detector->image_center[1];
    *x -= in->detector->image_center[0];
    return 0;
  }else{
    return -1;
  }
  return -1;
}


Image * sp_background_adaptative_mesh(Image * a,int cols, int rows, int slices){
  Image * cell_min;
  if(cols > sp_image_x(a)){
    cols = sp_image_x(a);
  }
  if(rows > sp_image_y(a)){
    rows = sp_image_y(a);
  }
  if(slices > sp_image_z(a)){
    slices = sp_image_z(a);
  }
  if(cols < 1){
    cols = 1;
  }
  if(rows < 1){
    rows = 1;
  }

  if(slices < 1){
    slices = 1;
  }
  cell_min = sp_image_alloc(cols,rows,slices);
  sp_image_rephase(cell_min,SP_ZERO_PHASE);


  /* Do not be scared by the deep indentation!
     The first three for loops just loop around the different cells
     The second group of three for loops loops around the pixels inside the cell to find the minimum of the cell
   */
  for(int c = 0;c<cols;c++){
    for(int r = 0;r<rows;r++){
      for(int s = 0;s<slices;s++){
	sp_image_set(cell_min,c,r,s,sp_cinit(1e9,0));	
	for(int x = c*sp_image_x(a)/cols;x<(c+1)*sp_image_x(a)/cols;x++){
	  for(int y = r*sp_image_y(a)/rows;y<(r+1)*sp_image_y(a)/rows;y++){
	    for(int z = s*sp_image_z(a)/slices;z<(s+1)*sp_image_z(a)/slices;z++){

	      if(sp_cabs(sp_image_get(a,x,y,z))	< sp_cabs(sp_image_get(cell_min,c,r,s))){
	 	sp_image_set(cell_min,c,r,s,sp_image_get(a,x,y,z));
	      }
	    }
	  }
	}
      }
    }
  }
  sp_image_write(cell_min,"cell_min.tif",0);
  Image * ret = bilinear_rescale(cell_min,sp_image_x(a),sp_image_y(a),sp_image_z(a));
  sp_image_free(cell_min);
  return ret;
}


Image * sp_image_phase_shift(Image * a, real phi, int in_place){
  Image * out;
  if(in_place){
    out = a;
  }else{
    out = sp_image_duplicate(a,SP_COPY_ALL);
  }
  for(int i = 0;i<sp_image_size(a);i++){
    real mag = sp_cabs(a->image->data[i]);
    real arg = sp_carg(a->image->data[i]);
    sp_real(out->image->data[i]) = mag*cos(arg+phi);
    sp_imag(out->image->data[i]) = mag*sin(arg+phi);
  }
  return out;
}

int sp_image_contains_coordinates(const Image *a, real x, real y, real z){
  if(x < 0 || x>=sp_image_x(a) || y< 0 || y>=sp_image_y(a) || z<0 || z>=sp_image_z(a)){
    return 0;
  }
  return 1;
}

void sp_image_image_to_mask(Image *in, Image *out)
{
  if (sp_image_x(in) != sp_image_x(out) ||
      sp_image_y(in) != sp_image_y(out) ||
      sp_image_z(in) != sp_image_z(out)) {
    printf("incompatible sizes in sp_image_image_to_mask\n");
    return;
  }
  const int i_max = sp_image_size(in);
  for (int i = 0; i < i_max; i++) {
    if (sp_real(in->image->data[i]) == 0.0 &&
	sp_imag(in->image->data[i]) == 0.0) {
      out->mask->data[i] = 0;
    } else {
      out->mask->data[i] = 1;
    }
  }
}
/* out will get a new image that is zero when in->mask is zero
     and (1.0,0.0) where it is nonzero */
void sp_image_mask_to_image(Image *in, Image *out)
{
  if (sp_image_x(in) != sp_image_x(out) ||
      sp_image_y(in) != sp_image_y(out) ||
      sp_image_z(in) != sp_image_z(out)) {
    printf("incompatible sizes in sp_image_mask_to_image\n");
    return;
  }
  const int i_max = sp_image_size(in);
  for (int i = 0; i < i_max; i++) {
    if (in->mask->data[i] == 0) {
      out->image->data[i] = sp_cinit(0.0,0.0);
    } else {
      out->image->data[i] = sp_cinit(1.0,0.0);
    }
  }
}
/* copies the image from in to out */
void sp_image_image_to_image(Image *in, Image *out)
{
  if (sp_image_x(in) != sp_image_x(out) ||
      sp_image_y(in) != sp_image_y(out) ||
      sp_image_z(in) != sp_image_z(out)) {
    printf("incompatible sizes in sp_image_image_to_image\n");
    return;
  }
  const int i_max = sp_image_size(in);
  for (int i = 0; i < i_max; i++) {
    out->image->data[i] = in->image->data[i];
  }
}
/* copies the mask from in to out */
void sp_image_mask_to_mask(Image *in, Image *out)
{
  if (sp_image_x(in) != sp_image_x(out) ||
      sp_image_y(in) != sp_image_y(out) ||
      sp_image_z(in) != sp_image_z(out)) {
    printf("incompatible sizes in sp_image_mask_to_mask\n");
    return;
  }
  const int i_max = sp_image_size(in);
  for (int i = 0; i < i_max; i++) {
    out->mask->data[i] = in->mask->data[i];
  }
}

/* out gets a new mask that is zero where in->image is nonzero
   and one wher it is zero */
void sp_image_invert_mask(Image *in, Image *out)
{
  if (sp_image_x(in) != sp_image_x(out) ||
      sp_image_y(in) != sp_image_y(out) ||
      sp_image_z(in) != sp_image_z(out)) {
    printf("incompatible sizes in sp_image_mask_to_mask\n");
    return;
  }
  const int i_max = sp_image_size(in);
  for (int i = 0; i < i_max; i++) {
    if (in->mask->data[i] == 0) {
      out->mask->data[i] = 1;
    } else {
      out->mask->data[i] = 0;
    }
  }
}

void sp_image_grow_mask(Image *in, int pixels){
  Image *tmp1 = sp_image_duplicate(in,SP_COPY_DATA|SP_COPY_MASK|SP_COPY_DETECTOR);
  Image *tmp2 = sp_image_duplicate(in,SP_COPY_DATA|SP_COPY_MASK|SP_COPY_DETECTOR);
  Image *foo;

  for (int i = 0; i < pixels; i++) {
    for (int x = 0; x < sp_image_x(in); x++) {
      for (int y = 0; y < sp_image_y(in); y++) {
	for (int z = 0; z < sp_image_z(in); z++) {
	  
	  if ((sp_image_mask_get(tmp1,x,y,z) == 1) ||
	      ((x != 0 && sp_image_mask_get(tmp1,x-1,y,z) != 0) ||
	       (x != sp_image_x(tmp1)-1 && sp_image_mask_get(tmp1,x+1,y,z) != 0) ||
	       (y != 0 && sp_image_mask_get(tmp1,x,y-1,z) != 0) ||
	       (y != sp_image_y(tmp1)-1 && sp_image_mask_get(tmp1,x,y+1,z) != 0) ||
	       (z != 0 && sp_image_mask_get(tmp1,x,y,z-1) != 0) ||
	       (z != sp_image_z(tmp1)-1 && sp_image_mask_get(tmp1,x,y,z+1) != 0))) {
	    sp_image_mask_set(tmp2,x,y,z,1);
	  } else {
	    sp_image_mask_set(tmp2,x,y,z,0);
	  }
	}
      }
    }
    foo = tmp1;
    tmp1 = tmp2;
    tmp2 = foo;
  }
  sp_image_mask_to_mask(tmp1,in);
  sp_image_free(tmp1);
  sp_image_free(tmp2);
}

void sp_image_shrink_mask(Image *in, int pixels){
  Image *tmp1 = sp_image_duplicate(in,SP_COPY_DATA|SP_COPY_MASK|SP_COPY_DETECTOR);
  Image *tmp2 = sp_image_duplicate(in,SP_COPY_DATA|SP_COPY_MASK|SP_COPY_DETECTOR);
  Image *foo;

  for (int i = 0; i < pixels; i++) {
    for (int x = 0; x < sp_image_x(in); x++) {
      for (int y = 0; y < sp_image_y(in); y++) {
	for (int z = 0; z < sp_image_z(in); z++) {
	  if ((sp_image_mask_get(tmp1,x,y,z) == 0) ||
	      ((x != 0 && sp_image_mask_get(tmp1,x-1,y,z) == 0) ||
	       (x != sp_image_x(tmp1)-1 && sp_image_mask_get(tmp1,x+1,y,z) == 0) ||
	       (y != 0 && sp_image_mask_get(tmp1,x,y-1,z) == 0) ||
	       (y != sp_image_y(tmp1)-1 && sp_image_mask_get(tmp1,x,y+1,z) == 0) ||
	       (z != 0 && sp_image_mask_get(tmp1,x,y,z-1) == 0) ||
	       (z != sp_image_z(tmp1)-1 && sp_image_mask_get(tmp1,x,y,z+1) == 0))) {
	    sp_image_mask_set(tmp2,x,y,z,0);
	  } else {
	    sp_image_mask_set(tmp2,x,y,z,1);
	  }
	}
      }
    }
    foo = tmp1;
    tmp1 = tmp2;
    tmp2 = foo;
  }
  sp_image_mask_to_mask(tmp1,in);
  sp_image_free(tmp1);
  sp_image_free(tmp2);
}
