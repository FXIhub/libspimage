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
static void write_h5_img(Image * img,const char * filename, int output_precision);
static Image * _read_imagefile(const char * filename, char * file, int line);
static Image * read_tiff(const char * filename);
static  void write_tiff(Image * img,const char * filename);
static  void write_csv(Image * img,const char * filename);
static Image * read_png(const char * filename);
static int write_png(Image * img,const char * filename, int color);
static int write_vtk(Image * img,const char * filename);
static Image * read_smv(const char * filename);


void sp_srand(int i){
  srand(i);
}

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
    return reflect_origo(in,in_place);
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
      sp_c3matrix_set(out->image,x,y,0,sp_c3matrix_get(in->image,x,y,0));
      sp_c3matrix_set(out->image,x,y,0,tmp);
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



/*
  For unshifted images it shifted the quadrants around image_center and 
  extra zero pad is used on the smaller quadrants to make them the same size
  as the biggest quadrant.
  For shifted images it shifted the quadrants around (size[]-1)/2 
  This function is generalized for 3D patterns as well.
*/
Image * sp_image_shift(Image * img){
  Image * out;
  int i;
  int index1,index2;
  int x,y,z;
  int newx,newy,newz;
  real max_x,max_y,max_z;//changed from int

  /* fft shift the image */
  out = sp_image_duplicate(img,SP_COPY_DATA|SP_COPY_MASK);
  if(!img->shifted){
    max_x = sp_max(img->detector->image_center[0],sp_c3matrix_x(img->image)-1-img->detector->image_center[0]);
    max_y = sp_max(img->detector->image_center[1],sp_c3matrix_y(img->image)-1-img->detector->image_center[1]);
    max_z = sp_max(img->detector->image_center[2],sp_c3matrix_z(img->image)-1-img->detector->image_center[2]);
    //was 2*max_x+1 before (other way of defining center)
    /* FM BUG: There was a bug here for 2D images with max_z so I added this hack */
    if(max_z == 0){
      max_z = 0.5;
    }
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
    if(resolution*2-1 > sp_c3matrix_x(res->image) || resolution*2-1 > sp_c3matrix_y(res->image)){
      return sp_image_duplicate(img,SP_COPY_DATA|SP_COPY_MASK);
    }  
    sp_c3matrix_free(res->image);
    sp_i3matrix_free(res->mask);
    res->image = sp_c3matrix_alloc(resolution*2-1,resolution*2-1,1);
    res->mask = sp_i3matrix_alloc(resolution*2-1,resolution*2-1,1);
    nx = 0;
    ny = 0;
    for(x = 0;x<sp_c3matrix_x(img->image);x++){
      dx = x;
      if(sp_c3matrix_x(img->image)-x-1 < dx){
      dx = sp_c3matrix_x(img->image)-x;
      }
      for(y = 0;y<sp_c3matrix_y(img->image);y++){
	dy = y;
	if(sp_c3matrix_y(img->image)-y-1 < dy){
	  dy = sp_c3matrix_y(img->image)-y;
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
    if(resolution*2-1 > sp_c3matrix_x(res->image) ||
       resolution*2-1 > sp_c3matrix_y(res->image) ||
       resolution*2-1 > sp_c3matrix_z(res->image)){
      return sp_image_duplicate(img,SP_COPY_DATA|SP_COPY_MASK);
    }
    sp_c3matrix_free(res->image);
    sp_i3matrix_free(res->mask);
    res->image = sp_c3matrix_alloc(resolution*2-1,resolution*2-1,resolution*2-1);
    res->mask = sp_i3matrix_alloc(resolution*2-1,resolution*2-1,resolution*2-1);
    nx = 0;
    ny = 0;
    nz = 0;
    for(x = 0;x<sp_c3matrix_x(img->image);x++){
      dx = x;
      if(sp_c3matrix_x(img->image)-x-1 < dx){
	dx = sp_c3matrix_x(img->image)-x;
      }
      for(y = 0;y<sp_c3matrix_y(img->image);y++){
	dy = y;
	if(sp_c3matrix_y(img->image)-y-1 < dy){
	  dy = sp_c3matrix_y(img->image)-y;
	}
	for(z = 0;z<sp_c3matrix_z(img->image);z++){
	  dz = z;
	  if(sp_c3matrix_z(img->image)-z-1 < dz){
	    dz = sp_c3matrix_z(img->image)-z;
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
    Image * blur_mask = gaussian_blur(tmp,*value);

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
    blur_mask = gaussian_blur(tmp,*value);
    printf("blur_mask (%i,%i,%i)\n",sp_image_x(blur_mask),sp_image_y(blur_mask),sp_image_z(blur_mask));

    sp_c3matrix_mul_elements(img->image,blur_mask->image);
  }  
}

real sp_centro_sym_value(int index,Image * img){
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

int sp_centro_sym_index(int index,Image * img){
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
  int x = i/sp_c3matrix_y(in->image)/sp_c3matrix_z(in->image);
  int y = i/sp_c3matrix_z(in->image)%sp_c3matrix_y(in->image);
  int z = i%sp_c3matrix_z(in->image)%sp_c3matrix_y(in->image);
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
  int x = i/sp_c3matrix_z(in->image)/sp_c3matrix_y(in->image);
  int y = i/sp_c3matrix_z(in->image)%sp_c3matrix_y(in->image);
  int z = i%sp_c3matrix_z(in->image)%sp_c3matrix_y(in->image);
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
  int x = i/sp_c3matrix_z(in->image)/sp_c3matrix_y(in->image);
  int y = i/sp_c3matrix_z(in->image)%sp_c3matrix_y(in->image);
  int z = i%sp_c3matrix_z(in->image)%sp_c3matrix_y(in->image);
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
  int x = i%sp_c3matrix_z(in->image)%sp_c3matrix_y(in->image);
  int y = i/sp_c3matrix_x(in->image)%sp_c3matrix_z(in->image);
  int z = i/sp_c3matrix_x(in->image)/sp_c3matrix_y(in->image);
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
    sp_real(img->image->data[i]) = cos(phase)*sp_cabs(img->image->data[i]);
    sp_imag(img->image->data[i]) = sin(phase)*sp_cabs(img->image->data[i]);
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





void sp_image_add(Image * a, Image * b){
  sp_c3matrix_add(a->image,b->image,NULL);
}

void sp_image_sub(Image * a, Image * b){
  sp_c3matrix_sub(a->image,b->image);
}


void sp_image_fill(Image * a, Complex value){
  for(int i = 0;i<sp_image_size(a);i++){
    a->image->data[i] = value;
  }
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

void _sp_image_free(Image * in, char * file, int line){
  _sp_free(in->detector,file,line);
  _sp_c3matrix_free(in->image,file,line);
  _sp_i3matrix_free(in->mask,file,line);
  _sp_free(in,file,line);
}

Image * _sp_image_duplicate(Image * in, int flags,char * file, int line){
  Image  *res = _sp_malloc(sizeof(Image),file,line);
  if(!res){
    sp_error_fatal("Out of memory!");
  }
  memcpy(res,in,sizeof(Image));
  res->detector = _sp_malloc(sizeof(Detector),file,line);
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


Image * _sp_image_alloc(int x, int y, int z,char * file, int line){
  Image  *res = _sp_malloc(sizeof(Image),file,line);
  if(!res){
    perror("Out of memory!\n");
    abort();
  }
  res->detector = _sp_malloc(sizeof(Detector),file,line);
  res->detector->image_center[0] = x/2.0;
  res->detector->image_center[1] = y/2.0;
  res->detector->image_center[2] = z/2.0;
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
    if(img->num_dimensions == SP_3D){
      fprintf(stderr,"Cannot export 3D file to csv");
    }
    write_csv(img,filename);
  }else{
    fprintf(stderr,"Unsupported file type: %s\n",filename);
  }
}

Image * _sp_image_read(const char * filename, int flags, char * file, int line){
  char buffer[1024];
  strcpy(buffer,filename);
  for(int i = 0;i<strlen(buffer);i++){
    buffer[i] = tolower(buffer[i]);
  }
  /* select the correct function depending on the filename extension */
  if(rindex(buffer,'.') && strcmp(rindex(buffer,'.'),".h5") == 0){
    /* we have an h5 file */
    return _read_imagefile(filename,file,line);
  }else if(rindex(buffer,'.') && strcmp(rindex(buffer,'.'),".png") == 0){
    /* we  have a png file */
    return read_png(filename);
  }else if(rindex(buffer,'.') && strcmp(rindex(buffer,'.'),".vtk") == 0){
    /* we have a vtk file */
    fprintf(stderr,"Cannot read VTK files!\n");
    return NULL;
  }else if(rindex(buffer,'.') && (strcmp(rindex(buffer,'.'),".tif") == 0 ||strcmp(rindex(buffer,'.'),".tiff") == 0 )){
    /* we have a tiff file */
    return read_tiff(filename);
  }else if(rindex(buffer,'.') && (strcmp(rindex(buffer,'.'),".smv") == 0 ||strcmp(rindex(buffer,'.'),".SMV") == 0 )){
    /* we have an smv file */
    return read_smv(filename);
  }else{
    fprintf(stderr,"Unsupported file type: %s\n",filename);
  }
  return NULL;
}


/* superseeded by the new one. Just leave it here in case strange bugs show up in the new one */
/*
static void write_h5_img(Image * img,const char * filename, int output_precision){
  hid_t dataspace_id;
  hid_t dataset_id;
  hid_t file_id;
  int status;
  int version;
  hsize_t  dims[2];
  real values[2];
  sp_3matrix * tmp;
  int i;
  hid_t out_type_id = 0;
  hid_t mem_type_id = 0;
  hid_t plist;
  hsize_t chunk_size[2] = {sp_c3matrix_x(img->image),sp_c3matrix_y(img->image)};
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

  dims[0] = sp_c3matrix_x(img->image);
  dims[1] = sp_c3matrix_y(img->image);
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

  tmp = sp_3matrix_alloc(sp_c3matrix_x(img->image),sp_c3matrix_y(img->image),1);
  for(i = 0;i<sp_image_size(img);i++){
    tmp->data[i] = sp_real(img->image->data[i]);
  }

  dataset_id = H5Dcreate(file_id, "/real", out_type_id,
			 dataspace_id, plist);
  status = H5Dwrite(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, tmp->data);
  status = H5Dclose(dataset_id);
  sp_3matrix_free(tmp);

  if(img->phased){
    tmp = sp_3matrix_alloc(sp_c3matrix_x(img->image),sp_c3matrix_y(img->image),1);
    for(i = 0;i<sp_image_size(img);i++){
      tmp->data[i] = sp_imag(img->image->data[i]);
    }

    dataset_id = H5Dcreate(file_id, "/imag",out_type_id,
			   dataspace_id, plist);
    status = H5Dwrite(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, tmp->data);
    status = H5Dclose(dataset_id);
    sp_3matrix_free(tmp);

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

  values[0] = img->detector->pixel_size[0];
  values[1] = img->detector->pixel_size[1];
  values[2] = img->detector->pixel_size[2];
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
*/

static void write_h5_img(Image * img,const char * filename, int output_precision){
  hid_t dataspace_id;
  hid_t dataset_id;
  hid_t file_id;
  int status;
  int version;
  hsize_t  dims[3];
  real values[3];
  sp_3matrix * tmp;
  int i;
  hid_t out_type_id = 0;
  hid_t mem_type_id = 0;
  hid_t plist;
  hsize_t chunk_size[3] = {sp_c3matrix_x(img->image),sp_c3matrix_y(img->image),sp_c3matrix_z(img->image)};
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

  dims[0] = sp_c3matrix_x(img->image);
  dims[1] = sp_c3matrix_y(img->image);
  dims[2] = sp_c3matrix_z(img->image);
  file_id = H5Fcreate(filename,  H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  dataspace_id = H5Screate_simple( 3, dims, NULL );

  plist = H5Pcreate (H5P_DATASET_CREATE);
  H5Pset_chunk(plist,3,chunk_size);
  H5Pset_deflate(plist,6);

  dataset_id = H5Dcreate(file_id, "/mask", H5T_NATIVE_INT,
			 dataspace_id, plist);
  status = H5Dwrite(dataset_id,H5T_NATIVE_INT , H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, img->mask->data);
  status = H5Dclose(dataset_id);

  tmp = sp_3matrix_alloc(sp_c3matrix_x(img->image),sp_c3matrix_y(img->image),
			 sp_c3matrix_z(img->image));
  for(i = 0;i<sp_image_size(img);i++){
    tmp->data[i] = sp_real(img->image->data[i]);
  }

  dataset_id = H5Dcreate(file_id, "/real", out_type_id,
			 dataspace_id, plist);
  status = H5Dwrite(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, tmp->data);
  status = H5Dclose(dataset_id);
  sp_3matrix_free(tmp);

  if(img->phased){
    tmp = sp_3matrix_alloc(sp_c3matrix_x(img->image),sp_c3matrix_y(img->image),
			   sp_c3matrix_z(img->image));
    for(i = 0;i<sp_image_size(img);i++){
      tmp->data[i] = sp_imag(img->image->data[i]);
    }

    dataset_id = H5Dcreate(file_id, "/imag",out_type_id,
			   dataspace_id, plist);
    status = H5Dwrite(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, tmp->data);
    status = H5Dclose(dataset_id);
    sp_3matrix_free(tmp);

  }
  dims[0] = 3;
  dataspace_id = H5Screate_simple( 1, dims, NULL );
  dataset_id = H5Dcreate(file_id, "/image_center",out_type_id ,
			 dataspace_id, H5P_DEFAULT);
  values[0] = img->detector->image_center[0];
  values[1] = img->detector->image_center[1];
  values[2] = img->detector->image_center[2];
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

  values[0] = img->detector->pixel_size[0];
  values[1] = img->detector->pixel_size[1];
  values[2] = img->detector->pixel_size[2];
  dataset_id = H5Dcreate(file_id, "/pixel_size", out_type_id,
			 dataspace_id, H5P_DEFAULT);
  status = H5Dwrite(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, values);
  status = H5Dclose(dataset_id);

  values[0] = img->num_dimensions;
  dataset_id = H5Dcreate(file_id, "/num_dimensions", out_type_id,
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


Image * _read_imagefile(const char * filename, char * file, int line){
  Image * res = _sp_malloc(sizeof(Image),file, line);
  int file_id,dataset_id,space;
  int status,i;
  int version;
  hsize_t dims[3];
  hid_t mem_type_id = 0;
  H5E_auto_t func;
  void * client_data;
  real values[3];
  sp_3matrix * tmp;
  if(sizeof(real) == sizeof(float)){
    mem_type_id = H5T_NATIVE_FLOAT;
  }else if(sizeof(real) == sizeof(double)){
    mem_type_id = H5T_NATIVE_DOUBLE;
  }else{
    abort();
  }
  
  
  
  res->detector = _sp_malloc(sizeof(Detector),file,line);
  
  
  H5Eget_auto(&func,&client_data);
  /* turn off warning to check file and version because they might not exist */
  H5Eset_auto(NULL,NULL);  

  file_id = H5Fopen(filename,H5F_ACC_RDONLY,H5P_DEFAULT);
  if(file_id < 0){
    sp_error_fatal("Unable to open %s",filename);
  }
  

  dataset_id = H5Dopen(file_id, "/version");
  /* File includes version information */
  if(dataset_id>=0){
    status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
		     H5P_DEFAULT, &version);
    if(status < 0){
      sp_error_fatal("Unable to read dataset");
    }
    if(version == 2){
      dataset_id = H5Dopen(file_id, "/mask");
      if(dataset_id < 0){
	sp_error_fatal("Unable to open dataset");
      }
      space = H5Dget_space(dataset_id);
      H5Sget_simple_extent_dims(space,dims,NULL);
      if(H5Sget_simple_extent_ndims(space) == 3){
	res->image = _sp_c3matrix_alloc(dims[0],dims[1],dims[2],file,line);
	res->mask = _sp_i3matrix_alloc(dims[0],dims[1],dims[2],file,line);
      }else{
	res->image = _sp_c3matrix_alloc(dims[0],dims[1],1,file,line);
	res->mask = _sp_i3matrix_alloc(dims[0],dims[1],1,file,line);
      }
      
      status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
		       H5P_DEFAULT, res->mask->data);
      if(status < 0){
	sp_error_fatal("Unable to read dataset");
      }


      status = H5Dclose(dataset_id);
      
      dataset_id = H5Dopen(file_id, "/image_center");
      if(dataset_id < 0){
	sp_error_fatal("Unable to open dataset");
      }
      status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		       H5P_DEFAULT, values);
      if(status < 0){
	sp_error_fatal("Unable to read dataset");
      }
      status = H5Dclose(dataset_id);
      res->detector->image_center[0] = values[0];
      res->detector->image_center[1] = values[1];
      if(values[2]){
	res->detector->image_center[2] = values[2];
      }else{
	res->detector->image_center[2] = 0;
      }
      
      dataset_id = H5Dopen(file_id, "/phased");
      if(dataset_id < 0){
	sp_error_fatal("Unable to open dataset");
      }

      status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		       H5P_DEFAULT, values);
      if(status < 0){
	sp_error_fatal("Unable to read dataset");
      }

      status = H5Dclose(dataset_id);
      res->phased = values[0];
      
      dataset_id = H5Dopen(file_id, "/shifted");
      if(dataset_id < 0){
	sp_error_fatal("Unable to open dataset");
      }

      status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		       H5P_DEFAULT, values);
      if(status < 0){
	sp_error_fatal("Unable to read dataset");
      }

      status = H5Dclose(dataset_id);
      res->shifted = values[0];
      
      dataset_id = H5Dopen(file_id, "/scaled");
      if(dataset_id < 0){
	sp_error_fatal("Unable to open dataset");
      }

      status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		       H5P_DEFAULT, values);
      if(status < 0){
	sp_error_fatal("Unable to read dataset");
      }

      status = H5Dclose(dataset_id);
      res->scaled = values[0];
      
      dataset_id = H5Dopen(file_id, "/detector_distance");
      if(dataset_id < 0){
	sp_error_fatal("Unable to open dataset");
      }

      status = H5Dread(dataset_id,  mem_type_id, H5S_ALL, H5S_ALL,
		       H5P_DEFAULT, values);
      if(status < 0){
	sp_error_fatal("Unable to read dataset");
      }

      status = H5Dclose(dataset_id);
      res->detector->detector_distance = values[0];
      
      dataset_id = H5Dopen(file_id, "/lambda");
      if(dataset_id < 0){
	sp_error_fatal("Unable to open dataset");
      }

      status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		       H5P_DEFAULT, values);
      if(status < 0){
	sp_error_fatal("Unable to read dataset");
      }

      status = H5Dclose(dataset_id);
      res->detector->lambda = values[0];
      
      dataset_id = H5Dopen(file_id, "/pixel_size");
      if(dataset_id < 0){
	sp_error_fatal("Unable to open dataset");
      }

      status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		       H5P_DEFAULT, values);
      if(status < 0){
	sp_error_fatal("Unable to read dataset");
      }

      status = H5Dclose(dataset_id);
      res->detector->pixel_size[0] = values[0];
      res->detector->pixel_size[1] = values[1];
      res->detector->pixel_size[2] = values[2];
      

      H5Eget_auto(&func,&client_data);
      /* turn off warning to check num_dimensions because it might not exist */
      H5Eset_auto(NULL,NULL);
      dataset_id = H5Dopen(file_id, "/num_dimensions");
      H5Eset_auto(func,client_data);
      if(dataset_id>=0){
	status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
			 H5P_DEFAULT, values);
	if(status < 0){
	  sp_error_fatal("Unable to read dataset");
	}

	status = H5Dclose(dataset_id);
	res->num_dimensions = values[0];
      }else{
	/* we'll try to guess the dimensions */
	res->num_dimensions = 0;
	if(sp_image_x(res) > 1){
	  res->num_dimensions++;
	}
	if(sp_image_y(res) > 1){
	  res->num_dimensions++;
	}
	if(sp_image_z(res) > 1){
	  res->num_dimensions++;
	}	  
      }
      
      if(res->phased){
	tmp = _sp_3matrix_alloc(sp_i3matrix_x(res->mask),
			       sp_i3matrix_y(res->mask),
			       sp_i3matrix_z(res->mask),file,line);
	dataset_id = H5Dopen(file_id, "/imag");
	if(dataset_id < 0){
	  sp_error_fatal("Unable to open dataset");
	}

	status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
			 H5P_DEFAULT, tmp->data);
      if(status < 0){
	sp_error_fatal("Unable to read dataset");
      }

	status = H5Dclose(dataset_id);
	for(i = 0;i<sp_3matrix_size(tmp);i++){
	  sp_imag(res->image->data[i]) = tmp->data[i];
	  sp_real(res->image->data[i]) = 0;
	}
	sp_3matrix_free(tmp);
      }
      
      tmp = _sp_3matrix_alloc(sp_i3matrix_x(res->mask),sp_i3matrix_y(res->mask),
			     sp_i3matrix_z(res->mask),file,line);
      dataset_id = H5Dopen(file_id, "/real");
      if(dataset_id < 0){
	sp_error_fatal("Unable to open dataset");
      }

      status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		       H5P_DEFAULT, tmp->data);
      if(status < 0){
	sp_error_fatal("Unable to read dataset");
      }

      status = H5Dclose(dataset_id);
      for(i = 0;i<sp_3matrix_size(tmp);i++){
	sp_real(res->image->data[i]) += tmp->data[i];
      }
      sp_3matrix_free(tmp);
      
      status = H5Fclose(file_id);
    }
  }else{
    /* File does *NOT* includes version information */
    dataset_id = H5Dopen(file_id, "/mask");
      if(dataset_id < 0){
	sp_error_fatal("Unable to open dataset");
      }

    space = H5Dget_space(dataset_id);
    res->num_dimensions = H5Sget_simple_extent_ndims(space);
    H5Sget_simple_extent_dims(space,dims,NULL);
    if(H5Sget_simple_extent_ndims(space) == 3){
      res->image = _sp_c3matrix_alloc(dims[0],dims[1],dims[2],file,line);
      res->mask = _sp_i3matrix_alloc(dims[0],dims[1],dims[2],file,line);
    }else if(H5Sget_simple_extent_ndims(space) == 2){
      res->image = _sp_c3matrix_alloc(dims[0],dims[1],1,file,line);
      res->mask = _sp_i3matrix_alloc(dims[0],dims[1],1,file,line);
    }else{
      fprintf(stderr,"File has unsupported number of dimensions!\n");
      abort();
    }
    tmp = _sp_3matrix_alloc(sp_i3matrix_x(res->mask),sp_i3matrix_y(res->mask),
			   sp_i3matrix_z(res->mask),file,line);
    
    status = H5Dread(dataset_id,mem_type_id , H5S_ALL, H5S_ALL,
		     H5P_DEFAULT, tmp->data);
      if(status < 0){
	sp_error_fatal("Unable to read dataset");
      }

    for(i = 0;i<sp_3matrix_size(tmp);i++){
      res->mask->data[i] = tmp->data[i];
    }
    sp_3matrix_free(tmp);
    
    status = H5Dclose(dataset_id);
    
    dataset_id = H5Dopen(file_id, "/image_center");
      if(dataset_id < 0){
	sp_error_fatal("Unable to open dataset");
      }

    status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		     H5P_DEFAULT, values);
      if(status < 0){
	sp_error_fatal("Unable to read dataset");
      }

    status = H5Dclose(dataset_id);
    res->detector->image_center[0] = values[0];
    res->detector->image_center[1] = values[1];
    if(values[2]){
      res->detector->image_center[2] = values[2];
    }else{
      res->detector->image_center[2] = 0;
    }
    
    dataset_id = H5Dopen(file_id, "/phased");
      if(dataset_id < 0){
	sp_error_fatal("Unable to open dataset");
      }

    status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		     H5P_DEFAULT, values);
      if(status < 0){
	sp_error_fatal("Unable to read dataset");
      }

    status = H5Dclose(dataset_id);
    res->phased = values[0];
    
    dataset_id = H5Dopen(file_id, "/shifted");
      if(dataset_id < 0){
	sp_error_fatal("Unable to open dataset");
      }

    status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		     H5P_DEFAULT, values);
      if(status < 0){
	sp_error_fatal("Unable to read dataset");
      }

    status = H5Dclose(dataset_id);
    res->shifted = values[0];
    
    dataset_id = H5Dopen(file_id, "/scaled");
      if(dataset_id < 0){
	sp_error_fatal("Unable to open dataset");
      }

    status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		     H5P_DEFAULT, values);
      if(status < 0){
	sp_error_fatal("Unable to read dataset");
      }

    status = H5Dclose(dataset_id);
    res->scaled = values[0];
    
    dataset_id = H5Dopen(file_id, "/detector_distance");
      if(dataset_id < 0){
	sp_error_fatal("Unable to open dataset");
      }

    status = H5Dread(dataset_id,  mem_type_id, H5S_ALL, H5S_ALL,
		     H5P_DEFAULT, values);
      if(status < 0){
	sp_error_fatal("Unable to read dataset");
      }

    status = H5Dclose(dataset_id);
    res->detector->detector_distance = values[0];
    
    dataset_id = H5Dopen(file_id, "/lambda");
      if(dataset_id < 0){
	sp_error_fatal("Unable to open dataset");
      }

    status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		     H5P_DEFAULT, values);
      if(status < 0){
	sp_error_fatal("Unable to read dataset");
      }

    status = H5Dclose(dataset_id);
    res->detector->lambda = values[0];
    
    dataset_id = H5Dopen(file_id, "/pixel_size");
      if(dataset_id < 0){
	sp_error_fatal("Unable to open dataset");
      }

    status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		     H5P_DEFAULT, values);
      if(status < 0){
	sp_error_fatal("Unable to read dataset");
      }

    status = H5Dclose(dataset_id);
    res->detector->pixel_size[0] = values[0];
    res->detector->pixel_size[1] = values[0];
    res->detector->pixel_size[2] = values[0];

    
    if(res->phased){
      tmp = _sp_3matrix_alloc(sp_i3matrix_x(res->mask),sp_i3matrix_y(res->mask),
			     sp_i3matrix_z(res->mask),file,line);
      dataset_id = H5Dopen(file_id, "/complex");
      if(dataset_id < 0){
	sp_error_fatal("Unable to open dataset");
      }

      status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		       H5P_DEFAULT, tmp->data);
      if(status < 0){
	sp_error_fatal("Unable to read dataset");
      }

      status = H5Dclose(dataset_id);
      for(i = 0;i<sp_3matrix_size(tmp);i++){
	sp_imag(res->image->data[i]) = tmp->data[i];
	sp_real(res->image->data[i]) = 0;
      }
      sp_3matrix_free(tmp);
      tmp = _sp_3matrix_alloc(sp_i3matrix_x(res->mask),sp_i3matrix_y(res->mask),
			     sp_i3matrix_z(res->mask),file,line);
      dataset_id = H5Dopen(file_id, "/real");
      if(dataset_id < 0){
	sp_error_fatal("Unable to open dataset");
      }

      status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
		       H5P_DEFAULT, tmp->data);
      if(status < 0){
	sp_error_fatal("Unable to read dataset");
      }

      status = H5Dclose(dataset_id);
      for(i = 0;i<sp_3matrix_size(tmp);i++){
	sp_real(res->image->data[i]) += tmp->data[i];
      }
      sp_3matrix_free(tmp);
      
    }else{
      if(!res->scaled){
	tmp = _sp_3matrix_alloc(sp_i3matrix_x(res->mask),sp_i3matrix_y(res->mask),
			       sp_i3matrix_z(res->mask),file,line);
	dataset_id = H5Dopen(file_id, "/intensities");
      if(dataset_id < 0){
	sp_error_fatal("Unable to open dataset");
      }

	status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
			 H5P_DEFAULT, tmp->data);
      if(status < 0){
	sp_error_fatal("Unable to read dataset");
      }

	status = H5Dclose(dataset_id);
	for(i = 0;i<sp_3matrix_size(tmp);i++){
	  sp_real(res->image->data[i]) += tmp->data[i];
	}
	sp_3matrix_free(tmp);
      }else{
	tmp = _sp_3matrix_alloc(sp_i3matrix_x(res->mask),sp_i3matrix_y(res->mask),
			       sp_i3matrix_z(res->mask),file,line);
	dataset_id = H5Dopen(file_id, "/amplitudes");
      if(dataset_id < 0){
	sp_error_fatal("Unable to open dataset");
      }

	status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL,
			 H5P_DEFAULT, tmp->data);
      if(status < 0){
	sp_error_fatal("Unable to read dataset");
      }

	status = H5Dclose(dataset_id);
	for(i = 0;i<sp_3matrix_size(tmp);i++){
	  sp_real(res->image->data[i]) += tmp->data[i];
	}
	sp_3matrix_free(tmp);	 
      }
    }        
    status = H5Fclose(file_id);    
  }
  H5Eset_auto(func,client_data);
  return res;
  
}



Image * read_tiff(const char * filename){
  Image * out = sp_malloc(sizeof(Image));
  out->detector = sp_malloc(sizeof(Detector));
  int bpp = 4;  
  int datatype = 0;
  int width,height;
  int nstrips;
  int stripsize;
  int i;
  unsigned char * img;
  float * tmpf;
  short * tmpi;
  unsigned short * tmpui;
  unsigned char * tmpuc;

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
  img = sp_malloc(nstrips*stripsize);
  for(i = 0;i<nstrips;i++){
    TIFFReadEncodedStrip(tif,i,img+i*stripsize,stripsize);
  }
  TIFFClose(tif);
  
  out->image = sp_c3matrix_alloc(width,height,1);
  out->mask = sp_i3matrix_alloc(width,height,1);

  if(datatype == SAMPLEFORMAT_UINT){
    tmpui = (unsigned short *)img;
    for(i = 0;i<sp_c3matrix_size(out->image);i++){
      sp_real(out->image->data[i])= tmpui[i];
      sp_imag(out->image->data[i]) = 0;
    }
  }else if(datatype == SAMPLEFORMAT_IEEEFP){
    tmpf = (float *)img;
    for(i = 0;i<sp_c3matrix_size(out->image);i++){
      sp_real(out->image->data[i]) = tmpf[i];
      sp_imag(out->image->data[i]) = 0;
    }
  }else if(datatype == SAMPLEFORMAT_VOID){
    tmpuc = (unsigned char *)img;
    for(i = 0;i<sp_c3matrix_size(out->image);i++){
      sp_real(out->image->data[i]) = tmpuc[i];
      sp_imag(out->image->data[i]) = 0;
    }
  }else if(datatype == SAMPLEFORMAT_INT){
    tmpi = (short *)(img);
    for(i = 0;i<sp_c3matrix_size(out->image);i++){
      sp_real(out->image->data[i]) = tmpi[i];
      sp_imag(out->image->data[i]) = 0;
    }
  }



  for(i = 0;i<sp_c3matrix_size(out->image);i++){
    out->mask->data[i] = 1;
  }
  sp_free(img);
  out->scaled = 0;
  out->phased = 0;
  out->shifted = 0;
  out->detector->image_center[0] = width/2;
  out->detector->image_center[1] = height/2;
  out->num_dimensions = SP_2D;
  return out;
}


void write_tiff(Image * img,const char * filename){
  float * data;
  int nstrips;
  int stripsize;
  TIFF * tif;
  int x,y;
  int width = sp_image_x(img);
  int height = sp_image_y(img);


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
  data = sp_malloc(nstrips*stripsize);
  for(y = 0;y<sp_image_x(img);y++){
    for(x = 0;x<sp_image_y(img);x++){
      data[x] =sp_cabs(sp_image_get(img,x,y,0));      
    }
    TIFFWriteEncodedStrip(tif,y,data,stripsize);
  }

  /*  for(i = 0;i<nstrips;i++){
    TIFFWriteEncodedStrip(tif,i,&(data[i*stripsize/4]),stripsize);
    }*/
  TIFFClose(tif);
  sp_free(data);
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
  for(y = 0;y<sp_image_y(img);y++){
    for(x = 0;x<sp_image_x(img);x++){
      fprintf(f,"%d,%d,%f,%f,%f,%f\n",x,y,sp_cabs(sp_image_get(img,x,y,0)),sp_carg(sp_image_get(img,x,y,0)),sp_real(sp_image_get(img,x,y,0)),sp_imag(sp_image_get(img,x,y,0)));
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


/* Convolute the image with a gaussian filter.
 The filter function is given by:

f(x,y) = 1/sqrt(2*M_PI*radius) * exp(-(x^2+y^2)/(2*radius^2)) */
Image * gaussian_blur(Image * in, real radius){
  /* Lets make this convolution using a fourier transform shallw we... good....*/
  int x,y,z;
  int i,j,k;
  int filter_side = ceil(radius)*3*2+1;
  real radius_z;
  if(in->num_dimensions == SP_2D){
    radius_z = 0;
  }else if(in->num_dimensions == SP_3D){
    radius_z = radius;
  }else{
    radius_z = 0;
    abort();
  }
  real total_filter = 0;
  Image * filter_img = sp_image_alloc(filter_side,filter_side, ceil(radius_z)*3*2+1);

  Image * centered_filter;
  Image * res;
  Image * tmp;
  filter_img->detector->image_center[0] = (sp_image_x(filter_img)-1)/2.0;
  filter_img->detector->image_center[1] = (sp_image_y(filter_img)-1)/2.0;
  filter_img->detector->image_center[2] = (sp_image_z(filter_img)-1)/2.0;
  
  sp_image_dephase(filter_img);
  for(x = -ceil(radius)*3;x<=ceil(radius)*3;x++){
    i = x+ceil(radius)*3;
    for(y = -ceil(radius)*3;y<=ceil(radius)*3;y++){
      j = y+ceil(radius)*3;
      for(z = -ceil(radius_z)*3;z<=ceil(radius_z)*3;z++){
	k = z+ceil(radius_z)*3;
	sp_real(filter_img->image->data[k*filter_side*filter_side+j*filter_side+i]) = 1/sqrt(2*M_PI*radius) * exp(-(x*x+y*y+z*z)/(2*radius*radius));
	sp_imag(filter_img->image->data[k*filter_side*filter_side+j*filter_side+i]) = 0;
	/* Make the filter symmetric in the imaginary part */
	/*      filter_img->image->data[i*filter_side+j] = filter_img->image->data[i*filter_side+j] + filter_img->image->data[i*filter_side+j]*I;*/
	total_filter += sp_cabs(filter_img->image->data[k*filter_side*filter_side+j*filter_side+i]);
      }
    }
  }
  for(i = 0;i<sp_image_size(filter_img);i++){
    filter_img->image->data[i] = sp_cscale(filter_img->image->data[i],1.0/total_filter);
  }
  centered_filter = shift_center_to_top_left(filter_img);
  centered_filter->shifted = 1;
  sp_image_free(filter_img);
  res = sp_image_convolute(in, centered_filter,NULL);
  sp_image_free(centered_filter);
  /* we should crop the result if it's bigger than the input */
  if(sp_image_size(res) > sp_image_size(in)){
    /*tmp = cube_crop(res, (sp_c3matrix_x(res->image)-sp_c3matrix_x(in->image))/2,
		    (sp_c3matrix_y(res->image)-sp_c3matrix_y(in->image))/2,
		    (sp_c3matrix_z(res->image)-sp_c3matrix_z(in->image))/2,
		    sp_c3matrix_x(in->image)/2-1+(sp_c3matrix_x(res->image)-sp_c3matrix_x(in->image))/2,
		    sp_c3matrix_x(in->image)/2-1+(sp_c3matrix_y(res->image)-sp_c3matrix_y(in->image))/2,
		    sp_c3matrix_x(in->image)/2-1+(sp_c3matrix_z(res->image)-sp_c3matrix_z(in->image))/2);*/
    tmp = cube_crop(res, (sp_c3matrix_x(res->image)-sp_c3matrix_x(in->image))/2,
		    (sp_c3matrix_y(res->image)-sp_c3matrix_y(in->image))/2,
		    (sp_c3matrix_z(res->image)-sp_c3matrix_z(in->image))/2,
		    (sp_c3matrix_x(res->image)+sp_c3matrix_x(in->image))/2-1,
		    (sp_c3matrix_y(res->image)+sp_c3matrix_y(in->image))/2-1,
		    (sp_c3matrix_z(res->image)+sp_c3matrix_z(in->image))/2-1);
    sp_image_free(res);
    res = tmp;
  }
  return res;
}


/* Convolute the image with a square window.
 The filter function is given by:

f(x,y) = 1/((2*radius+1)^2)) */
Image * square_blur(Image * in, real radius, int type){
  /* Lets make this convolution using a fourier transform shallw we... good....*/
  int x,y,z;
  int i,j,k;
  int filter_side = ceil(radius)*3*2+1;
  Image * filter_img = sp_image_duplicate(in,SP_COPY_DETECTOR);
  sp_c3matrix * filter = sp_c3matrix_alloc(filter_side,filter_side,(in->num_dimensions == SP_2D)? 1:filter_side);
  filter_img->detector->image_center[2] = (filter_side-1)/2.0;
  if(in->num_dimensions == SP_2D){
    filter_img->detector->image_center[2] = 0;
  }
  real total_filter = 0;
  Image * centered_filter;
  Image * res;
  filter_img->detector->image_center[0] = (filter_side-1)/2.0;
  filter_img->detector->image_center[1] = (filter_side-1)/2.0;
  sp_c3matrix_free(filter_img->image);
  filter_img->image = filter;
  sp_image_dephase(filter_img);
  for(x = -ceil(radius)*3;x<=ceil(radius)*3;x++){
    i = x+ceil(radius)*3;
    for(y = -ceil(radius)*3;y<=ceil(radius)*3;y++){
      j = y+ceil(radius)*3;
      for(z = -ceil(radius)*3;z<=ceil(radius)*3;z++){
	if(type == SP_3D){k = z+ceil(radius)*3;}else{k = 0;}
	filter->data[k*filter_side*filter_side+j*filter_side+i] = sp_cinit(1.0/((2*radius+1)*(2*radius+1)),0);
	total_filter += sp_real(filter->data[k*filter_side*filter_side+j*filter_side+i]);
      }
    }
  }
  for(i = 0;i<sp_image_size(filter_img);i++){
    filter_img->image->data[i] = sp_cscale(filter_img->image->data[i],1.0/total_filter);
  }
  centered_filter = shift_center_to_top_left(filter_img);
  sp_image_free(filter_img);
  res = sp_image_convolute(in, centered_filter,NULL);
  sp_image_free(centered_filter);
  return res;
}


int write_mask_to_png(Image * img, char * filename, int color){
  Image  * res = sp_image_duplicate(img,SP_COPY_DATA|SP_COPY_MASK);
  int ret = 1;
  int i;
  if(sp_i3matrix_z(img->mask) != 1){
    sp_image_free(res);
    fprintf(stderr,"Can't write 3D mask to png");
  }else{
    for(i = 0;i<sp_image_size(img);i++){
      res->image->data[i] = sp_cinit(res->mask->data[i],0);
    }
    ret = write_png(res,filename,color);
    sp_image_free(res);
  }
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
 row_pointers = sp_malloc(sizeof(png_byte *)*height);
 for(i = 0;i<height;i++){
   row_pointers[i] = sp_malloc(sizeof(png_byte)*width*bit_depth/8);
 }
 png_read_image(png_ptr, row_pointers);
 res = sp_image_alloc(width,height,1);
 for(i = 0;i<height;i++){
   for(j = 0;j<width;j++){
     res->image->data[i*width+j] = sp_cinit(row_pointers[i][(int)(j*bit_depth/8)],0);
     res->mask->data[i*width+j] = 1;
   }
 }
 return res;
}




int write_png(Image * img,const char * filename, int color){

  if(img->num_dimensions != SP_2D){
    fprintf(stderr,"Can only write png of 2D images in write_png!\n");
    abort();
  }

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
  png_set_IHDR(png_ptr, info_ptr, sp_c3matrix_x(img->image), sp_c3matrix_y(img->image),
	       bit_depth, color_type, interlace_type,
	       compression_type, filter_method);
  
  

  row_pointers = png_malloc(png_ptr,sp_c3matrix_y(img->image)*sizeof(png_byte *));
  for (i=0; i<sp_c3matrix_y(img->image); i++){
    row_pointers[i] = png_malloc(png_ptr,sp_c3matrix_x(img->image)*pixel_size*sizeof(png_byte));
  }
  
  /* We're gonna scale the image so that it fits on the 8 bits */
  min_v = sp_c3matrix_min(img->image,NULL);
  max_v = sp_c3matrix_max(img->image,NULL);
  if(max_v-min_v){
    scale = 1/(max_v-min_v);
  }else{
    scale = 1;
  }
  offset = min_v;
  i = 0;
  log_of_2 = log(2.0);
  /* this is a special kind of color */
  for(y = 0;y<sp_c3matrix_y(img->image);y++){
    for(x = 0;x<sp_c3matrix_x(img->image);x++){
      /* traditional color scale taken from gnuplot manual */
      if(color & LOG_SCALE){
	value = 1-log((sp_cabs(img->image->data[i])-offset)*scale+FLT_EPSILON)/log(FLT_EPSILON);
      }else{
	value = ((sp_cabs(img->image->data[i])-offset)*scale);
      }
      if(color & COLOR_PHASE){
	phase = (256*(sp_carg(img->image->data[i])+3.1416)/(2*3.1416));
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
  for(i=0; i<sp_c3matrix_y(img->image); i++){
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
    if(!fobs->mask->data[i] || sp_real(fobs->image->data[i]) < low_intensity_cutoff){
      continue;
    }
    num +=  sp_cabs(sp_csub(fobs->image->data[i],fcalc->image->data[i]));
    den += sp_cabs(fobs->image->data[i]);
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
      fft_img->image->data[i] = sp_cinit(0,0);
    }
  }
  tmp = sp_image_duplicate(fft_img,SP_COPY_DATA|SP_COPY_MASK);
  for(i = 0;i<sp_image_size(tmp);i++){
    tmp->image->data[i] = sp_cinit(log(sp_cabs(tmp->image->data[i])+1),0);
  }
  write_png(tmp,"low_pass.png",COLOR_JET);
  sp_image_free(tmp);

  res = sp_image_ifft(fft_img);
  /* scale appropriately */
  for(i = 0;i<sp_image_size(res);i++){
    res->image->data[i] = sp_cscale(res->image->data[i],sp_image_size(res));
  }
  sp_image_free(fft_img);
  fft_img = sp_image_fft(res);
  tmp = sp_image_duplicate(fft_img,SP_COPY_DATA|SP_COPY_MASK);
  for(i = 0;i<sp_image_size(tmp);i++){
    tmp->image->data[i] = sp_cinit(log(sp_cabs(tmp->image->data[i])+1),0);
  }
  //  write_png(tmp,"after_low_pass.png",COLOR_JET); //not compatible with 3D
  sp_image_free(tmp);
  
  return res;
}

//I am here
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
    res->image->data[i] = sp_cscale(res->image->data[i],1.0/sp_image_size(res));
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
    if(sp_cabs(mask->image->data[i]) < sp_image_size(res)-1){
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
    res->image->data[i] = sp_cscale(res->image->data[i],exp(-scaling*scaling_factor));
  }
  return res;
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

/* Superseeded by the new write_vtk */
/*
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
  fprintf(f,"DIMENSIONS %d %d 1\n",sp_c3matrix_x(img->image),sp_c3matrix_y(img->image));
  fprintf(f,"ORIGIN 0 %d 0\n",sp_c3matrix_y(img->image));
  fprintf(f,"SPACING 1 -1 1\n");
  fprintf(f,"POINT_DATA %lld\n",sp_image_size(img));
  fprintf(f,"SCALARS amplitudes float 1\n");
  fprintf(f,"LOOKUP_TABLE default\n");
  fprintf(f,"%6g",sp_cabs(img->image->data[0]));
  for(y = 0;y<sp_c3matrix_y(img->image);y++){
    for(x = 0;x<sp_c3matrix_x(img->image);x++){
      fprintf(f," %g",sp_cabs(img->image->data[y*sp_c3matrix_x(img->image)+x]));    
    }
  }
*/
/*  for(i = 1;i<sp_image_size(img);i++){
    fprintf(f," %g",img->image->data[i]);    
  }*/
/*
  fprintf(f,"\n");
  fflush(f);
  fclose(f);
  return 0;
}
*/

int write_vtk(Image * img, const char * filename){
  FILE * f = fopen(filename,"w");
  int x,y,z;
  if(!f){
    perror("Bad file in write_vtk!");
    abort();
  }
  fprintf(f,"# vtk DataFile Version 2.0\n");
  fprintf(f,"Generated by image_util write_vtk()\n");
  fprintf(f,"ASCII\n");
  fprintf(f,"DATASET STRUCTURED_POINTS\n");
  fprintf(f,"DIMENSIONS %d %d %d\n",sp_c3matrix_x(img->image),
	  sp_c3matrix_y(img->image),sp_c3matrix_z(img->image));
  fprintf(f,"ORIGIN 0 0 0\n");//changed from 0 y 0
  fprintf(f,"SPACING 1 1 1\n");//changed from 1 -1 1 when going to 3d ??
  fprintf(f,"POINT_DATA %lld\n",sp_image_size(img));
  fprintf(f,"SCALARS amplitudes float 1\n");
  fprintf(f,"LOOKUP_TABLE default\n");
  fprintf(f,"%6g",sp_cabs(img->image->data[0]));
  for(z = 0;z<sp_c3matrix_z(img->image);z++){
    for(y = 0;y<sp_c3matrix_y(img->image);y++){
      for(x = 0;x<sp_c3matrix_x(img->image);x++){
	fprintf(f," %g",sp_cabs(img->image->data[z*sp_c3matrix_x(img->image)*sp_c3matrix_y(img->image)+y*sp_c3matrix_x(img->image)+x]));
      }
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
  sp_c3matrix_free(cropped->image);
  cropped->image = sp_c3matrix_alloc(y2-y1+1,x2-x1+1,sp_c3matrix_z(in->image));
  sp_i3matrix_free(cropped->mask);
  cropped->mask = sp_i3matrix_alloc(y2-y1+1,x2-x1+1,sp_c3matrix_z(in->image));


  for(i = y1;i<= y2;i++){
    memcpy(&cropped->image->data[(i-y1)*sp_c3matrix_x(cropped->image)],&in->image->data[(i)*sp_c3matrix_x(in->image)+x1],sp_c3matrix_x(cropped->image)*sizeof(real));
    memcpy(&cropped->mask->data[(i-y1)*sp_c3matrix_x(cropped->image)],&in->mask->data[(i)*sp_c3matrix_x(in->image)+x1],sp_c3matrix_x(cropped->image)*sizeof(real));
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
  cropped->image = sp_c3matrix_alloc(y2-y1+1,x2-x1+1,z2-z1+1);
  sp_i3matrix_free(cropped->mask);
  cropped->mask = sp_i3matrix_alloc(y2-y1+1,x2-x1+1,z2-z1+1);


  for(i = z1;i<=z2;i++){
    for(j = y1;j<=y2;j++){
      memcpy(&cropped->image->data[(i-z1)*sp_c3matrix_y(cropped->image)*sp_c3matrix_x(cropped->image)+(j-y1)*sp_c3matrix_x(cropped->image)],
	     &in->image->data[(i)*sp_c3matrix_y(in->image)*sp_c3matrix_x(in->image)+(j)*sp_c3matrix_x(in->image)+x1],sp_c3matrix_x(cropped->image)*sizeof(real));
      memcpy(&cropped->mask->data[(i-z1)*sp_c3matrix_y(cropped->image)*sp_c3matrix_x(cropped->image)+(j-y1)*sp_c3matrix_x(cropped->image)],
	     &in->mask->data[(i)*sp_c3matrix_y(in->image)*sp_c3matrix_x(in->image)+(j)*sp_c3matrix_x(in->image)+x1],sp_c3matrix_x(cropped->image)*sizeof(real));
    }
  }
  return cropped;
}


void find_center(Image * img, real * center_x, real * center_y, real * center_z){
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
    write_png(a,"centrosym_convolve.png",COLOR_JET|LOG_SCALE);
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
  

  for(x = 0; x<sp_c3matrix_x(res->image); x++){
    virtual_x = (real)x*sp_c3matrix_x(img->image)/sp_c3matrix_x(res->image);
    for(y = 0; y<sp_c3matrix_y(res->image); y++){
      virtual_y = (real)y*sp_c3matrix_y(img->image)/sp_c3matrix_y(res->image);
      for(z = 0; z<sp_c3matrix_z(res->image); z++){
	virtual_z = (real)z*sp_c3matrix_z(img->image)/sp_c3matrix_z(res->image);
	res->image->data[z*sp_c3matrix_y(res->image)*sp_c3matrix_x(res->image)+
			 y*sp_c3matrix_x(res->image)+x] =
	  sp_c3matrix_interp(img->image,virtual_x,virtual_y,virtual_z);	
	res->mask->data[z*sp_c3matrix_y(res->image)*sp_c3matrix_x(res->image)+
			y*sp_c3matrix_x(res->image)+x] = 
	  sp_i3matrix_interp(img->mask,virtual_x,virtual_y,virtual_z);
      }
    }
  }
  return res;
}
  

real sp_image_interp(Image * img, real v_x, real v_y, real v_z){
  return sp_real(sp_c3matrix_interp(img->image,v_x,v_y,v_z));
}

void sp_image_scale(Image * img, real value){
  Complex tmp = {value,0};
  sp_c3matrix_scale(img->image,tmp);
}


real sp_image_max(Image * img, long long * index,int * x, int * y, int * z){
  real ret;
  ret = sp_c3matrix_max(img->image,index);
  if(index){
    if(x){
      *x = *index%sp_image_z(img)%sp_image_y(img);
    }
    if(y){
      *y = *index/sp_image_x(img)%sp_image_z(img);
    }
    if(z){
      *z = *index/sp_image_x(img)/sp_image_y(img);
    }
  }
  return ret;
}


void _sp_image_realloc(Image * img, int new_x, int new_y, int new_z, char * file, int line){
  _sp_c3matrix_realloc(img->image,new_x,new_y,new_z,file,line);
  _sp_i3matrix_realloc(img->mask,new_x,new_y,new_z,file,line);
}

Image * rectangular_window(int image_x, int image_y, int width, int height, int shifted){
  Image * res = sp_image_alloc(image_x,image_y,1);
  int x,y,i;
  int center[2];
  i = 0;
  
  if(shifted){
    for(x = 0;x<image_x;x++){
      for(y = 0;y<image_y;y++){
	if((fabs(x) < width/2 || fabs(image_x-x) < width/2 )&&
	   (fabs(y) < height/2 || fabs(image_y-y) < height/2 )){
	  res->image->data[i] = sp_cinit(1,0);
	}else{
	  res->image->data[i] = sp_cinit(0,0);
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
	  res->image->data[i] = sp_cinit(1,0);
	}else{
	  res->image->data[i] = sp_cinit(0,0);
	}
	i++;
      }
    }
  }
  sp_image_rephase(res,SP_ZERO_PHASE);
  return res;
}

Image * cube_window(int image_x, int image_y, int image_z, int dx, int dy, int dz, int shifted){
  Image * res = sp_image_alloc(image_x,image_y,image_z);
  int x,y,z,i;
  int center[3];
  i = 0;
  
  if(shifted){
    for(x = 0;x<image_x;x++){
      for(y = 0;y<image_y;y++){
	for(z = 0;z<image_z;z++){
	  if((fabs(x) < dx/2 || fabs(image_x-x) < dx/2 )&&
	     (fabs(y) < dy/2 || fabs(image_y-y) < dy/2 )&&
	     (fabs(z) < dz/2 || fabs(image_z-z) < dz/2 )){
	    res->image->data[i] = sp_cinit(1,0);
	  }else{
	    res->image->data[i] = sp_cinit(0,0);
	  }
	  i++;
	}
      }
    }
  }else{
    center[0] = image_x/2;
    center[1] = image_y/2;
    center[2] = image_z/2;
    for(x = 0;x<image_x;x++){
      for(y = 0;y<image_y;y++){
	for(z = 0;z<image_z;z++){
	  if(fabs(x-center[0]) < dx/2 &&
	     fabs(y-center[1]) < dy/2 &&
	     fabs(z-center[2]) < dz/2){
	    res->image->data[i] = sp_cinit(1,0);
	  }else{
	    res->image->data[i] = sp_cinit(0,0);
	  }
	  i++;
	}
      }
    }
  }
  sp_image_rephase(res,SP_ZERO_PHASE);
  return res;
}

Image * circular_window(int x, int y, int radius, int shifted){
  Image * res = sp_image_alloc(x,y,1);
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
      res->image->data[i] = sp_cinit(1,0);
    }else{
      res->image->data[i] = sp_cinit(0,0);
    }
  }
  sp_image_rephase(res,SP_ZERO_PHASE);
  return res;
}

Image * spherical_window(int x, int y, int z, int radius, int shifted){
  Image * res = sp_image_alloc(x,y,z);
  int i;
  if(shifted){
    res->detector->image_center[0] = 0;
    res->detector->image_center[1] = 0;
    res->detector->image_center[2] = 0;
  }else{
    res->detector->image_center[0] = x/2;
    res->detector->image_center[1] = y/2;
    res->detector->image_center[2] = z/2;
  }
  for(i = 0;i<x*y*z;i++){
    if(dist_to_center(i,res) < radius){
      res->image->data[i] = sp_cinit(1,0);
    }else{
      res->image->data[i] = sp_cinit(0,0);
    }
  }
  sp_image_rephase(res,SP_ZERO_PHASE);
  return res;
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


Image * sp_image_get_mask(Image * a){
  Image * res = sp_image_duplicate(a,SP_COPY_DATA|SP_COPY_MASK);
  int i;
  for(i = 0;i<sp_image_size(a);i++){
    res->image->data[i] = sp_cinit(res->mask->data[i],0);
  }
  return res;
}

real sp_point_convolute(Image * a, Image * b, int index){
  int index_x, index_y, index_z,x,y,z;
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
  index_x = index%sp_c3matrix_z(a->image)%sp_c3matrix_y(a->image)-
    a->detector->image_center[0];
  index_y = index/sp_c3matrix_x(a->image)%sp_c3matrix_z(a->image)-
    a->detector->image_center[1];
  index_z = index/sp_c3matrix_x(a->image)/sp_c3matrix_y(a->image)-
    a->detector->image_center[2];
  for(x = -b->detector->image_center[0];x<sp_c3matrix_x(b->image)-b->detector->image_center[0];x++){
    for(y = -b->detector->image_center[1];y<sp_c3matrix_y(b->image)-b->detector->image_center[1];y++){
      for(z = -b->detector->image_center[2];z<sp_c3matrix_z(b->image)-b->detector->image_center[2];z++){
	if(x+index_x < -a->detector->image_center[0] || 
	   x+index_x >=  sp_c3matrix_x(a->image)-a->detector->image_center[0] ||
	   y+index_y < -a->detector->image_center[1] || 
	   y+index_y >=  sp_c3matrix_y(a->image)-a->detector->image_center[1] ||
	   z+index_z < -a->detector->image_center[2] ||
	   z+index_z >=  sp_c3matrix_z(a->image)-a->detector->image_center[2]){
	  /* we're outside of image a */
	  continue;
	}
	ai = (index_z+z+a->detector->image_center[2])*sp_c3matrix_y(a->image)*sp_c3matrix_x(a->image)+(index_y+y+a->detector->image_center[1])*sp_c3matrix_x(a->image)+(index_x+x+a->detector->image_center[0]);
	bi = (z+b->detector->image_center[2])*sp_c3matrix_y(b->image)*sp_c3matrix_x(b->image)+(y+b->detector->image_center[1])*sp_c3matrix_x(b->image)+(x+b->detector->image_center[0]);
	out += sp_cabs(a->image->data[ai])*sp_cabs(b->image->data[bi]);
      }
    }
  }
  if(a->shifted){
    sp_image_free(tmp);
  }
  return out;
  }

int sp_image_shift_index(Image * a, long long index){
  int x = index%sp_c3matrix_z(a->image)%sp_c3matrix_y(a->image);
  int y = index/sp_c3matrix_x(a->image)%sp_c3matrix_z(a->image);
  int z = index/sp_c3matrix_x(a->image)/sp_c3matrix_y(a->image);
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



Image * sp_image_local_variance(Image * img, Image * window){
  int i;
  int size[3] = {sp_c3matrix_x(img->image)+sp_c3matrix_x(window->image)-1,sp_c3matrix_y(img->image)+sp_c3matrix_y(window->image)-1,sp_c3matrix_z(img->image)+sp_c3matrix_z(window->image)-1};
  Image * norm_window = sp_image_duplicate(window,SP_COPY_DATA|SP_COPY_MASK);
  
  sp_image_normalize(norm_window);
  Image * ra = sp_image_convolute(img,norm_window,size);
  Image * res = cube_crop(ra,0,0,0,sp_c3matrix_x(img->image)-1,sp_c3matrix_y(img->image)-1,sp_c3matrix_z(img->image)-1);
  /*  write_png(img,"non_averaged.png",COLOR_JET);
  write_png(ra,"total_averaged.png",COLOR_JET);
  write_png(res,"crop_total_averaged.png",COLOR_JET);*/
  for(i = 0;i<sp_image_size(res);i++){
    res->image->data[i] = sp_cinit(sp_cabs(sp_csub(res->image->data[i],img->image->data[i])),0);
  }
  sp_image_free(ra);
  ra = sp_image_convolute(res,norm_window,size);
  sp_image_free(res);
  res = cube_crop(ra,0,0,0,sp_c3matrix_x(img->image)-1,sp_c3matrix_y(img->image)-1,sp_c3matrix_z(img->image)-1);
  sp_image_free(ra);
  sp_image_free(norm_window);

  return res;
}


Complex sp_image_dot_prod(Image * a, Image * b){
  return sp_c3matrix_froenius_prod(a->image,b->image);
}

Image * sp_proj_module(Image * a, Image * b){
  Image * ret = sp_image_duplicate(a,SP_COPY_DATA|SP_COPY_MASK);
  int i;
  for(i = 0;i<sp_image_size(a);i++){
    if(b->mask->data[i]){
      ret->image->data[i] = sp_cscale(ret->image->data[i],sp_cabs(b->image->data[i])/sp_cabs(a->image->data[i]));
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
    if(!sp_cabs(b->image->data[i])){
      ret->image->data[i] = sp_cinit(0,0);
    }
  }
  return ret;
}


int sp_image_invert(Image * a){
  int i;
  for(i = 0;i<sp_image_size(a);i++){
    sp_real(a->image->data[i]) = 1.0/sp_real(a->image->data[i]);
  }
  return 0;
}


void sp_image_mul_elements(Image * a, Image * b){
  sp_c3matrix_mul_elements(a->image,b->image);  
}

void sp_image_conj(Image * a){  
  sp_c3matrix_conj(a->image);
}


void sp_image_memcpy(Image * dst,Image * src){  
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
  real real_to_reciprocal = 1.0/(in->detector->detector_distance*in->detector->lambda);
  real ewald_radius = 1.0/in->detector->lambda;
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


static Image * read_smv(const char * filename){
  Image * res = NULL;
  FILE * fp = fopen(filename,"rb");
  char buffer[1024];
  int header_size = 0;
  int x_size = 0;
  int y_size = 0;
  while(fgets(buffer,1024,fp)){
    /* stop loop when we find a line with original folder */
    if(strstr(buffer,"OriginalFolder")){
      break;
    }
    char * p;
    if(strstr(buffer,"HEADER_BYTES=")){
      p = strstr(buffer,"HEADER_BYTES=")+strlen("HEADER_BYTES=");
      header_size = atoi(p);
    }
    if(strstr(buffer,"SIZE1=")){
      p = strstr(buffer,"SIZE1=")+strlen("SIZE1=");
	x_size = atoi(p);
    }
    if(strstr(buffer,"SIZE2=")){
      p = strstr(buffer,"SIZE2=")+strlen("SIZE2=");
      y_size = atoi(p);
    }
  }
  if(!x_size || !y_size || !header_size){
    return NULL;
  }
  res = sp_image_alloc(x_size,y_size,1);
  fseek(fp,header_size,SEEK_SET);
  unsigned short * data = sp_malloc(sizeof(unsigned short)*x_size*y_size);
  fread((void *)data,sizeof(unsigned short),x_size*y_size,fp);  
  for(int x = 0; x < x_size;x++){
    for(int y = 0; y < y_size;y++){
      sp_image_set(res,x,y,0,sp_cinit(data[y*x_size+x],0));
      sp_i3matrix_set(res->mask,x,y,0,1);
    }
  }
  return res;
}



/*! Superimposes image b on top of image a 
 *
 *  flags is a bitwise combination of the following:
 *
 *  SP_ENANTIOMORPH - allow to try to superimpose not only b but also
 *  the "mirror image" of b [b(-x)].
 *
*/
void sp_image_superimpose(Image * a,Image * b, int flags){
  int x,y,z;
  long long index;
  int center_invert = 0;
  real max;
  /* check maximum overlap of a and b */
  Image * direct_overlap = sp_image_cross_correlate(a,b,NULL);
  max = sp_image_max(direct_overlap,&index,&x,&y,&z);
  sp_image_free(direct_overlap);
  if(flags & SP_ENANTIOMORPH){
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
      sp_image_reflect(b,IN_PLACE,SP_ORIGO);
    }
  }

  sp_image_translate(b,x,y,z,SP_TRANSLATE_WRAP_AROUND);
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
  
real sp_image_rs_r_factor(Image * a, Image * b){
  real sum_dif = 0;
  real sum_sum = 0;
  for(int i = 0; i < sp_image_size(a);i++){
    sum_sum += sp_cabs(a->image->data[i])+sp_cabs(b->image->data[i]);
    sum_dif += fabs(sp_cabs(a->image->data[i])-sp_cabs(b->image->data[i]));
  }
  return sum_dif/sum_sum;
}

real sp_image_correlation_coefficient(Image * a,Image * b){
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
