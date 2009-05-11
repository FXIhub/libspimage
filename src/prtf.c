#include "spimage.h"


void sp_image_fourier_translate(Image * a, real t_x, real t_y, real t_z){
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
				       sp_cinit(cos(f_x * t_x * nx_inv * two_pi+f_y * t_y * ny_inv * two_pi+f_z * t_z * nz_inv * two_pi),
						-sin(f_x * t_x * nx_inv * two_pi+f_y * t_y * ny_inv * two_pi+f_z * t_z * nz_inv * two_pi)));
	i++;
      }
    }
  }
}


Image * sp_prtf_basic(Image ** list,int n){
  if(n < 1){
    return NULL;
  }
  for(int i = 0;i<n;i++){
    if(sp_image_x(list[0]) != sp_image_x(list[i]) ||
       sp_image_y(list[0]) != sp_image_y(list[i]) ||
       sp_image_z(list[0]) != sp_image_z(list[i])){
      return NULL;
    }
       
  }
  Image * ret = sp_image_alloc(sp_image_x(list[0]),sp_image_y(list[0]),sp_image_z(list[0]));
  int size = sp_image_size(list[0]);
  sp_image_fill(ret,sp_cinit(0,0));
  for(int i = 0;i<n;i++){
    sp_image_add(ret,list[i]);
  }
  for(int i = 0;i<size;i++){
    real div = 0;
    for(int j = 0;j<n;j++){
      div += sp_cabs(list[j]->image->data[i]);
    }
    if(div){
      ret->image->data[i] = sp_cscale(ret->image->data[i],1.0/div);
    }
  }
  return ret;
}



Image * sp_prtf_advanced(Image ** list,int n,SpSpace space){
  Image ** real_image = sp_malloc(sizeof(Image *)*n);
  Image ** fourier_image = sp_malloc(sizeof(Image *)*n);
  /* We'll do the superposition and phase match in real space
     but the prtf in reciprocal space
  */
  if(space == RealSpace){
    for(int i = 0;i<n;i++){
      real_image[i] = sp_image_duplicate(list[i],SP_COPY_ALL);
    }
  }else{
    for(int i = 0;i<n;i++){
      real_image[i] = sp_image_ifft(list[i]);
    }
  }
  sp_vector * com_master = sp_image_center_of_mass(real_image[0]);
  for(int i = 1;i<n;i++){
    /* rough alignment */
    sp_image_superimpose(real_image[0],real_image[i],SP_ENANTIOMORPH);
    /* precise alignment */
    sp_vector * com = sp_image_center_of_mass(real_image[i]);
    Image * tmp = sp_image_fft(real_image[i]);
    sp_image_fourier_translate(tmp,sp_vector_get(com_master,0)-sp_vector_get(com,0),
			       sp_vector_get(com_master,1)-sp_vector_get(com,1),
			       sp_vector_get(com_master,2)-sp_vector_get(com,2));
    sp_image_free(real_image[i]);
    real_image[i] = sp_image_ifft(tmp);
    sp_image_free(tmp);
    sp_image_scale(real_image[i],1.0/sp_image_size(real_image[i]));

    real phi = sp_image_phase_match(real_image[0],real_image[i],2);
  }
  for(int i = 0;i<n;i++){
    fourier_image[i] = sp_image_fft(real_image[i]);
  }
  Image * ret = sp_prtf_basic(fourier_image,n);
  for(int i = 0;i<n;i++){
    sp_image_free(real_image[i]);
    sp_image_free(fourier_image[i]);
  }
  sp_free(fourier_image);
  sp_free(real_image);
  return ret;
}

sp_list * sp_prtf_by_resolution(Image * prtf){
  int size = sp_image_x(prtf)+sp_image_y(prtf)+sp_image_z(prtf);
  real * sum = sp_malloc(sizeof(real)*size);
  int * div = sp_malloc(sizeof(int)*size);
  for(int i = 0;i<size;i++){
    sum[i] = 0;
    div[i] = 0;
  }
  for(int i = 0;i<sp_image_size(prtf);i++){
    int d = sp_image_dist(prtf,i,SP_TO_CORNER);
    sum[d] += sp_cabs(prtf->image->data[i]);
    div[d]++;
  }
  int max_size = 0;
  for(int i = 0;i<size;i++){
    if(div[i]){
      max_size = i;
    }
  }
  sp_list * ret = sp_list_alloc(max_size);
  for(int i = 0;i<max_size;i++){
    if(div[i]){
      sp_list_set(ret,i,sum[i]/div[i]);
    }else{
      sp_list_set(ret,i,0);
    }
  }
  sp_free(sum);
  sp_free(div);
  return ret;
}
