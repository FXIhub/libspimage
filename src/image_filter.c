#ifdef _USE_DMALLOC
#include <dmalloc.h>
#endif

#include "spimage.h"



/* Convolute the image with a gaussian filter.
 The filter function is given by:

f(x,y) = 1/sqrt(2*M_PI*radius) * exp(-(x^2+y^2)/(2*radius^2)) */
Image * sp_gaussian_blur_old(Image * in, real radius){
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

Image * sp_gaussian_blur(Image * in, real radius){
  int x, y, z;
  Image *fourier_image = sp_image_fft(in);
  real coordinate_rad;
  for (int i = 0; i < sp_image_size(fourier_image); i++) {
    coordinate_rad = sp_image_dist(fourier_image, i, SP_TO_TOP_LEFT) / ((real)sp_image_x(fourier_image));
    fourier_image->image->data[i] = sp_cmul(fourier_image->image->data[i], sp_cinit(exp(-2.*pow(M_PI,2)*pow(coordinate_rad,2)*pow(radius,2))/((real) sp_image_size(fourier_image)), 0.));
  }

  Image *ret = sp_image_ifft(fourier_image);
  return ret;
}


/* Convolute the image with a square window.
 The filter function is given by:

f(x,y) = 1/((2*radius+1)^2)) */
Image * sp_square_blur(Image * in, real radius, int type){
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


/* Low pass filter using a centered square window of side edge_size */
Image * sp_low_pass_square_filter(Image * in, int edge_size){
  Image * fft_img = sp_image_fft(in);
  Image * res;
  Image * tmp;
  int i = 0;
  for(i = 0;i<sp_image_size(in);i++){
    if(sp_image_dist(in,i,SP_TO_CENTER2) > edge_size/2.0){
      fft_img->image->data[i] = sp_cinit(0,0);
    }
  }
  tmp = sp_image_duplicate(fft_img,SP_COPY_DATA|SP_COPY_MASK);
  for(i = 0;i<sp_image_size(tmp);i++){
    tmp->image->data[i] = sp_cinit(log(sp_cabs(tmp->image->data[i])+1),0);
  }
  //  sp_image_write(tmp,"low_pass.png",SpColormapJet);
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
Image * sp_low_pass_gaussian_filter(Image * in, int edge_size){
  Image * fft_img = sp_image_fft(in);
  Image * res;
  Image * mask;
  int i = 0;
  sp_gaussian_filter(fft_img,edge_size/2.0,1);
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
  sp_gaussian_filter(fft_img,edge_size/2.0,1);
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
Image * sp_gaussian_filter(Image * in, real radius,int in_place){
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
    scaling = (sp_image_dist(in,i,SP_TO_CENTER)/(radius))*(sp_image_dist(in,i,SP_TO_CENTER)/(radius));
    res->image->data[i] = sp_cscale(res->image->data[i],exp(-scaling*scaling_factor));
  }
  return res;
}


Image * sp_rectangular_window(int image_x, int image_y, int width, int height, int shifted){
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

Image * sp_cube_window(int image_x, int image_y, int image_z, int dx, int dy, int dz, int shifted){
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

Image * sp_circular_window(int x, int y, int radius, int shifted){
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
    if(sp_image_dist(res,i,SP_TO_CENTER) < radius){
      res->image->data[i] = sp_cinit(1,0);
    }else{
      res->image->data[i] = sp_cinit(0,0);
    }
  }
  sp_image_rephase(res,SP_ZERO_PHASE);
  return res;
}

Image * sp_spherical_window(int x, int y, int z, int radius, int shifted){
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
    if(sp_image_dist(res,i,SP_TO_CENTER) < radius){
      res->image->data[i] = sp_cinit(1,0);
    }else{
      res->image->data[i] = sp_cinit(0,0);
    }
  }
  sp_image_rephase(res,SP_ZERO_PHASE);
  return res;
} 

real sp_point_convolute(const Image * a,const Image * b, int index){
  real index_x, index_y, index_z;
  int x,y,z;
  real out = 0;
  int ai,bi;
  /* do shifted convolutions as normal convolutions */
  /*
  if(a->shifted){
    tmp = sp_image_shift(a);
    index = sp_image_shift_index(a,index);
    a = tmp;
  }
  
  if(b->shifted){
    fprintf(stderr,"Point convoluting with a shifted function is not currently defined!\n");
    return -1;
  }
  */
  sp_image_get_coords_from_index(a,index,&index_x,&index_y,&index_z,SpTopLeftCorner);
  /*  index_x = index%sp_c3matrix_z(a->image)%sp_c3matrix_y(a->image)-
    a->detector->image_center[0];
  index_y = index/sp_c3matrix_x(a->image)%sp_c3matrix_z(a->image)-
    a->detector->image_center[1];
  index_z = index/sp_c3matrix_x(a->image)/sp_c3matrix_y(a->image)-
  a->detector->image_center[2];*/
  /*  for(x = -b->detector->image_center[0];x<sp_c3matrix_x(b->image)-b->detector->image_center[0];x++){
    for(y = -b->detector->image_center[1];y<sp_c3matrix_y(b->image)-b->detector->image_center[1];y++){
    for(z = -b->detector->image_center[2];z<sp_c3matrix_z(b->image)-b->detector->image_center[2];z++){*/
  for(x = 0;x<sp_c3matrix_x(b->image);x++){
    for(y = 0;y<sp_c3matrix_y(b->image);y++){
      for(z = 0;z<sp_c3matrix_z(b->image);z++){
	if(!sp_image_contains_coordinates(a,x+index_x-b->detector->image_center[0],y+index_y-b->detector->image_center[1],z+index_z-b->detector->image_center[2])){
	  /* we're outside of image a */
	  continue;
	}
	ai = sp_image_get_index(a,x+index_x-b->detector->image_center[0],y+index_y-b->detector->image_center[1],z+index_z-b->detector->image_center[2]);
	bi = sp_image_get_index(b,x,y,z);
	out += sp_cabs(a->image->data[ai])*sp_cabs(b->image->data[bi]);
      }
    }
  }
  return out;
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

